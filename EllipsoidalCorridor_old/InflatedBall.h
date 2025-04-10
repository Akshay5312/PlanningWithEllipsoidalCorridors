#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>

#include "../PathParameterization/LagrangePolynomial.h"

namespace CorrGen{
    struct InflatedBallOptions{
        Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor
        Eigen::MatrixXd contain_points; // The points that must be contained in the ball

        Eigen::VectorXd reference_point; // The reference point for the ball
        int N_q; // The number of dimensions for the point cloud.

        double r_f = 1.0; // The diagonal dominance for the ball matrix
    };


    template<typename T>
    struct InflatedBall{
        
        Eigen::VectorX<T> d_; // Ball vector
        T f_; // Ball rad
        T r_; // Ball rad
        
        InflatedBall(Eigen::VectorX<T> d, T f, T r)
            : d_(d), f_(f), r_(r) {}

            T evaluate(Eigen::VectorXd point) const {
                // Evaluates if 'point' is in the Ball cross section of the corridor.
                
                T quadratic_term = (point.transpose() * point);
                
                T linear_term = ref_p.transpose() * point;
                
                return (quadratic_term + linear_term + f_) - r;
            }

        Eigen::VectorX<T> center() const {
            // Calculate the center of the Ball
            Eigen::MatrixX<T> E_inv = E_.inverse();
            return reference_point_ + (E_inv * d_);
        }

        std::vector<Eigen::VectorXd> sample(int N_samples) const {
            // Sample points within the Ball
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(E_);
            Eigen::MatrixXd C = eigen_solver.operatorSqrt();
            Eigen::MatrixXd C_inv = C.inverse();
            // Eigen::MatrixX<T> E_inv = E_.inverse();

            Eigen::VectorXd q_tilde = center();

            double o = (q_tilde - reference_point_).transpose() * E_ * (q_tilde - reference_point_);

            std::vector<Eigen::VectorXd> samples(N_samples);
            for (int i = 0; i < N_samples; ++i) {
                Eigen::VectorXd random_point = Eigen::VectorXd::Random(E_.rows());
                random_point.normalize();
                samples[i] = sqrt(1+o) * C_inv * random_point + q_tilde;
            }
            return samples;
        }
    };

    InflatedBall<double> GenerateInflatedBall(
        const InflatedBallOptions& options) {

        // drake::log()->info("Generating Ball");

        int N_q = options.N_q;
        
        double r_f = options.r_f;
        auto point_cloud = options.point_cloud; // Each col is a point
        auto contain_points = options.contain_points; // Each col is a point
        auto ref_point = options.reference_point; // Reference point for the Ball

        drake::solvers::MathematicalProgram prog;

        Eigen::VectorX<drake::symbolic::Variable> d_var = prog.NewContinuousVariables(N_q, 1, "d");
        drake::symbolic::Variable f_var = prog.NewContinuousVariables(1, 1, "f")[0];
        drake::symbolic::Variable r_var = prog.NewContinuousVariables(1, 1, "r")[0];

        InflatedBall<drake::symbolic::Expression> cal_C(d_var, f_var, r_var);

        // Add non-membership constraints for each point in the point cloud
        for (int p_index = 0; p_index < options.point_cloud.cols(); ++p_index) {
            Eigen::VectorXd point = options.point_cloud.col(p_index);

            // Non-membership constraint
            prog.AddLinearConstraint(cal_C.evaluate(point) >= 0);
        }

        // Add a cost to minimize the trace of the Ball matrix E and the scalar f
        prog.AddLinearCost(E.trace());

        // prog.AddLinearConstraint(f[0] == 0); // Ensure f is bounded


        auto result = drake::solvers::Solve(prog);

        if (!result.is_success()) {
            throw std::runtime_error("Ball generation failed");
        }

        // Extract the optimized values
        Eigen::MatrixX<double> E_opt = result.GetSolution(E_var);
        Eigen::VectorX<double> d_opt = result.GetSolution(d_var);

        // Create the Ball object
        InflatedBall<double> ball(E_opt, d_opt, ref_point);
        drake::log()->info("Ball generated successfully");

        return ball;
    }
}   // namespace CorrGen