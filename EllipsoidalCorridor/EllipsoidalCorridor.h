#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>

#include "LagrangePolynomial.h"

namespace CorrGen{
    struct EllipsoidalCorridorOptions{
        EllipsoidalCorridorOptions(){}
        Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor

        LagrangePolynomial<double> reference_path; // The reference path for the corridor
        int N_q; // The number of dimensions for the point cloud.

        int N_ctrl; // The number of control points for the ellipsoid
        int N_samples; // The number of times to sample the ellipsoid. 
        // The constraints are enforced at these points.

        double r_diagonal = 10.0; // The diagonal dominance for the ellipsoid matrix
        double r_f = 100.0; // The diagonal dominance for the ellipsoid matrix
    };


    template<typename T>
    struct EllipsoidalCorridor{
        
        LagrangePolynomial<T> E_; // Ellipsoid matrix
        LagrangePolynomial<T> d_; // Ellipsoid vector
        LagrangePolynomial<T> f_; // Ellipsoid scalar
        LagrangePolynomial<double> reference_path_; // Reference point for the ellipsoid

        EllipsoidalCorridor( LagrangePolynomial<T> E_path, LagrangePolynomial<T>d_path, LagrangePolynomial<T> f_path, LagrangePolynomial<double> ref_path)
            : E_(E_path), d_(d_path), f_(f_path), reference_path_(ref_path) {
                // drake::log()->info("Ellipsoid created");
            }

        T evaluate(double eps, Eigen::VectorXd point) const {
            // Evaluates if 'point' is in the ellipsoidal cross section of the corridor.
            // Corridor(e) = {q | (q-c_e)E_e(q-c_e) - 2*d_e(q-c_e) + f_e <= 1}

            Eigen::VectorXd refp = reference_path_.value(eps);

            auto o_point = point - refp; // Assuming c is the center of the ellipsoid
            // drake::log()->info("o_point: \n{}", o_point);
            // drake::log()->info("point: \n{}", point);
            // drake::log()->info("reference_path_: \n{}", reference_path_.value(eps));
            
            auto E = E_.value(eps);
            auto d = d_.value(eps);
            auto f = f_.value(eps);

            // drake::log()->info("E: \n{}", E);
            // drake::log()->info("d: \n{}", d);
            // drake::log()->info("f: \n{}", f);

            auto quadratic_term = o_point.transpose() * E * o_point;
            // drake::log()->info("quadratic_term: {}", quadratic_term);
            auto linear_term = - 2 * d.transpose() * o_point;
            // drake::log()->info("linear_term: {}", linear_term);
            auto constant_term = f;
            // drake::log()->info("constant_term: {}", constant_term);

            return (quadratic_term + linear_term + constant_term)(0,0) - 1.0;
        }

        Eigen::VectorX<T> center(double eps) const {
            // Calculate the center of the ellipsoid
            Eigen::MatrixX<T> E_inv = E_.value(eps).inverse();
            return reference_path_.value(eps) + (E_inv * d_.value(eps));
        }

    };

    EllipsoidalCorridor<double> GenerateEllipsoidalCorridor(
        const EllipsoidalCorridorOptions& options) {

        int N_q = options.N_q;
        int N_ctrl = options.N_ctrl;
        int N_samples = options.N_samples;
        double r_diagonal = options.r_diagonal;
        double r_f = options.r_f;
        auto point_cloud = options.point_cloud; // Each col is a point
        auto ref_path = options.reference_path; // Reference point for the ellipsoid

        drake::solvers::MathematicalProgram prog;

        std::vector<Eigen::MatrixX<drake::symbolic::Variable>> E_vars(N_ctrl); 
        std::vector<Eigen::MatrixX<drake::symbolic::Variable>> d_vars(N_ctrl); 
        std::vector<Eigen::MatrixX<drake::symbolic::Variable>> f_vars(N_ctrl);

        std::vector<Eigen::MatrixX<drake::symbolic::Expression>> E_expr(N_ctrl); 
        std::vector<Eigen::MatrixX<drake::symbolic::Expression>> d_expr(N_ctrl); 
        std::vector<Eigen::MatrixX<drake::symbolic::Expression>> f_expr(N_ctrl);
        std::vector<double> ctrl_breaks(N_ctrl);

        std::vector<double> sample_breaks(N_samples);
        for (int k = 0; k < N_samples; ++k) {
            sample_breaks[k] = static_cast<double>(k) / (N_samples-1); // Example control breaks
        }

        for (int k = 0; k < N_ctrl; ++k) {
            E_vars[k] = (prog.NewContinuousVariables(N_q, N_q, "E_" + std::to_string(k)));
            d_vars[k] = (prog.NewContinuousVariables(N_q, 1, "d_" + std::to_string(k)));
            f_vars[k] = (prog.NewContinuousVariables(1, 1, "f_" + std::to_string(k)));
            E_expr[k] = E_vars[k];
            d_expr[k] = d_vars[k];
            f_expr[k] = f_vars[k];
            ctrl_breaks[k] = static_cast<double>(k) / (N_ctrl-1); // Example control breaks
        }

        EllipsoidalCorridor<drake::symbolic::Expression> corr_C(
            LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks, E_expr),
            LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks, d_expr),
            LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks, f_expr),
            ref_path);

        for(int k = 0; k < N_samples; k++){
            double eps = sample_breaks[k];
            Eigen::MatrixX<drake::symbolic::Expression> E = corr_C.E_.value(eps);
            Eigen::VectorX<drake::symbolic::Expression> d = corr_C.d_.value(eps);
            Eigen::VectorX<drake::symbolic::Expression> f = corr_C.f_.value(eps);

            Eigen::VectorXd ref_point = corr_C.reference_path_.value(eps);

            for (int i = 0; i < N_q; ++i) {
                for (int j = 0; j < N_q; ++j) {
                    if (i != j) {
                        // Add constraints for diagonal dominance
                        prog.AddLinearConstraint(E(i, i) >= r_diagonal * E(i, j));
                        prog.AddLinearConstraint(E(i, i) >= -r_diagonal * E(i, j));
                    
                        // Add constraints for symmetry
                        prog.AddLinearConstraint(E(i, j) == E(j, i)); // Symmetry constraint
                    }
                }
                prog.AddLinearConstraint(E(i, i) >= 0.001);
                prog.AddLinearConstraint(E(i, i) <= 100.0); // Ensure E is bounded

                // prog.AddLinearConstraint(d(i) <= 2.5); // Ensure d is bounded
                // prog.AddLinearConstraint(d(i) >= -2.5); // Ensure d is bounded
            }

            prog.AddLinearConstraint(f[0] == 0.0); // Ensure f is bounded
            // prog.AddLinearConstraint(f[0] >= ); // Ensure f is bounded
            

            
            // Add non-membership constraints for each point in the point cloud
            for (int p_index = 0; p_index < options.point_cloud.cols(); ++p_index) {
                Eigen::VectorXd point = options.point_cloud.col(p_index);

                // Non-membership constraint
                prog.AddLinearConstraint(corr_C.evaluate(eps, point) >= 0);
            }

            // Add a cost to minimize the trace of the ellipsoid matrix E and the scalar f
            prog.AddLinearCost(E.trace() + N_q * r_f * f[0]);

            // prog.AddLinearConstraint(f[0] == 0); // Ensure f is bounded
        }

        auto result = drake::solvers::Solve(prog);

        if (!result.is_success()) {
            throw std::runtime_error("Ellipsoid generation failed");
        }

        // Extract the optimized values
        std::vector<Eigen::MatrixX<double>> E_opt(N_ctrl);
        std::vector<Eigen::MatrixX<double>> d_opt(N_ctrl);
        std::vector<Eigen::MatrixX<double>> f_opt(N_ctrl);

        for (int k = 0; k < N_ctrl; ++k) {
            E_opt[k] = result.GetSolution(E_vars[k]);
            d_opt[k] = result.GetSolution(d_vars[k]);
            f_opt[k] = result.GetSolution(f_vars[k]);
        }
        
        // Create the ellipsoid object
        EllipsoidalCorridor<double> ellipsoid_corridor(
            LagrangePolynomial<double>(ctrl_breaks, E_opt),
            LagrangePolynomial<double>(ctrl_breaks, d_opt),
            LagrangePolynomial<double>(ctrl_breaks, f_opt),
            ref_path
        );
    
        // drake::log()->info("Ellipsoid generated successfully");

        return ellipsoid_corridor;
    }
}   // namespace CorrGen