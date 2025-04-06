#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>

#include "InflatedEllipsoid.h"

#include "LagrangePolynomial.h"

namespace CorrGen{
    struct EllipsoidalCorridorOptions{
        EllipsoidalCorridorOptions(){}
        Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor

        LagrangePolynomial<double> reference_path; // The reference path for the corridor
        int N_q; // The number of dimensions for the point cloud.

        LagrangePolynomial<double> cost_path; // The path for the cost imposed on the ellipsoid

        int N_ctrl; // The number of control points for the ellipsoid
        int N_samples; // The number of times to sample the ellipsoid. 
        // The constraints are enforced at these points.

        double r_f = 1.0; // Cost on 'f'

        double r_min = 0.01; // Minimum radius of the ellipsoid

        double max_d = 0.25; // Maximum distance from the center of the ellipsoid to the reference in any direction
    };


    template<typename T>
    struct EllipsoidalCorridor{
        
        LagrangePolynomial<T> W_; // Ellipsoid matrix
        LagrangePolynomial<double> reference_path_; // Reference point for the ellipsoid

        EllipsoidalCorridor( LagrangePolynomial<T> W_path, LagrangePolynomial<double> ref_path)
            : W_(W_path), reference_path_(ref_path) 
            {
            }

        Eigen::MatrixX<T> W(double t) const {
            // Evaluate the ellipsoid matrix at time t
            return W_.value(t);
        }

        Eigen::VectorX<T> center(double t) const {
            // Calculate the center of the ellipsoid
            Eigen::MatrixX<T> W = W_(t);
            Eigen::MatrixX<T> E = W.topLeftCorner(W.rows() - 1, W.cols() - 1);
            Eigen::VectorX<T> d = W.topRightCorner(W.rows() - 1, 1);

            Eigen::MatrixX<T> E_inv = E.inverse();
            Eigen::VectorX<T> q_tilde = reference_path_.value(t) + (E_inv * d);
            return q_tilde;
        }

    };

    EllipsoidalCorridor<double> GenerateEllipsoidalCorridor(
        const EllipsoidalCorridorOptions& options) {

        int N_q = options.N_q;
        int N_ctrl = options.N_ctrl;
        int N_samples = options.N_samples;
        auto point_cloud = options.point_cloud; // Each col is a point
        auto ref_path = options.reference_path; // Reference point for the ellipsoid
        auto cost_path = options.cost_path; // 
        auto d_max = options.max_d; // Maximum distance from the center of the ellipsoid to the reference in any direction
        auto r_min = options.r_min; // Minimum radius of the ellipsoid


        drake::solvers::MathematicalProgram prog;

        std::vector<Eigen::MatrixX<drake::symbolic::Variable>> W_vars(N_ctrl);
        
        std::vector<Eigen::MatrixX<drake::symbolic::Expression>> W_expr(N_ctrl); 
        std::vector<double> ctrl_breaks(N_ctrl);
        std::vector<double> sample_breaks(N_samples);
        
        for (int k = 0; k < N_ctrl; ++k) {
            W_vars[k] = (prog.NewContinuousVariables(N_q+1, N_q+1, "W_" + std::to_string(k)));
            W_expr[k] = W_vars[k];
            ctrl_breaks[k] = static_cast<double>(k) / (N_ctrl-1); // Example control breaks
        }
        for (int k = 0; k < N_samples; ++k) {
            sample_breaks[k] = static_cast<double>(k) / (N_samples-1); // Example control breaks
        }


        EllipsoidalCorridor<drake::symbolic::Expression> corr_C(
            LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks, W_expr),
            ref_path);

        std::vector< InflatedEllipsoid<drake::symbolic::Expression> > corr_Cs(N_samples);

        for(int k = 0; k < N_samples; k++){
            double eps = sample_breaks[k];
            Eigen::MatrixX<drake::symbolic::Expression> W = corr_C.W(eps);
            
            Eigen::VectorXd ref_point = corr_C.reference_path_.value(eps);
            
            Eigen::VectorXd d_hat = cost_path.value(eps).normalized();

            corr_Cs[k] = InflatedEllipsoid<drake::symbolic::Expression>(W, ref_point);            

            auto& cal_C = corr_Cs[k];

            cal_C.EnforceIdentity(prog);
            cal_C.ApproximatelyMaximizeRadius(prog, d_hat);
            // cal_C.EnforceOffsetNorm(prog, d_max*d_max);
            // cal_C.ApproximatelyEnforceMinRadius(prog, 0.1, d_hat);
            cal_C.ApproximatelyEnforceOffsetNorm(prog, d_max, d_hat);
            
            for (int i = 0; i < point_cloud.cols(); ++i) {
                Eigen::VectorXd point = point_cloud.col(i);
                cal_C.AddNonMembershipConstraint(prog, point);
            }
        }

        corr_Cs[0].AddMembershipConstraint(prog, (Eigen::VectorXd)ref_path.value(0));
        corr_Cs[N_samples-1].AddMembershipConstraint(prog, (Eigen::VectorXd)ref_path.value(1));

        auto result = drake::solvers::Solve(prog);
        auto result_solver_details = result.get_solution_result();

        drake::log()->info("Solver details: {}", result_solver_details);

        if (!result.is_success()) {
            throw std::runtime_error("Ellipsoid generation failed");
        }

        // Extract the optimized values
        std::vector<Eigen::MatrixX<double>> W_opt(N_ctrl);

        for (int k = 0; k < N_ctrl; ++k) {
            W_opt[k] = result.GetSolution(W_vars[k]);
        }
        
        // Create the ellipsoid object
        EllipsoidalCorridor<double> ellipsoid_corridor(
            LagrangePolynomial<double>(ctrl_breaks, W_opt),
            ref_path
        );
    
        // drake::log()->info("Ellipsoid generated successfully");

        return ellipsoid_corridor;
    }
}   // namespace CorrGen