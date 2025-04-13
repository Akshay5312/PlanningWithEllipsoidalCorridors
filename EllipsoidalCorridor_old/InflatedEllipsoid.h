#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>

#include "../PathParameterization/LagrangePolynomial.h"

namespace CorrGen{
    struct InflatedEllipsoidOptions{
        Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor
        Eigen::MatrixXd contain_points; // The points that must be contained in the ellipsoid

        Eigen::VectorXd reference_point; // The reference point for the ellipsoid
        int N_q; // The number of dimensions for the point cloud.

        Eigen::VectorXd d_hat; // The direction of the ellipsoid

        double d_max; // The maximum offset for the ellipsoid

        double r_min = 0.1; // The minimum radius of the ellipsoid
    };


    template<typename T>
    struct InflatedEllipsoid{
        
        Eigen::MatrixX<T> W_; // Ellipsoid matrix
        Eigen::VectorX<double> reference_point_; // Reference point for the ellipsoid
        int N_q; // The number of dimensions for the point cloud.

        InflatedEllipsoid() : N_q(0) {}
        InflatedEllipsoid(Eigen::MatrixX<T> W, Eigen::VectorX<double> reference_point)
            : W_(W), reference_point_(reference_point), N_q(W.rows() - 1) {
            }

            
            T evaluate(Eigen::VectorXd point) const {
                // Evaluates if 'point' is in the ellipsoidal cross section of the corridor.
                
                // Corridor(e) = {q | [q-c; -1] * W * [q-c; -1] <= 1}
                // where c is the reference point. W is PSD.

                Eigen::VectorXd lifted_point = Eigen::VectorXd::Zero(point.size() + 1);
                lifted_point.head(point.size()) = point - reference_point_;
                lifted_point[point.size()] = -1.0;
                Eigen::VectorX<T> quadratic_term = lifted_point.transpose() * W_ * lifted_point;

                return (quadratic_term)[0] - 1.0;
            }

            drake::symbolic::Expression evaluate(Eigen::VectorX<drake::symbolic::Expression> point) const {
                // Evaluates if 'point' is in the ellipsoidal cross section of the corridor.
                
                // Corridor(e) = {q | [q-c; -1] * W * [q-c; -1] <= 1}
                // where c is the reference point. W is PSD.

                Eigen::VectorX<drake::symbolic::Expression> lifted_point 
                    = Eigen::VectorX<drake::symbolic::Expression>::Zero(point.size() + 1);
                lifted_point.head(point.size()) = point - reference_point_;
                lifted_point[point.size()] = -1.0;
                Eigen::VectorX<drake::symbolic::Expression> quadratic_term = lifted_point.transpose() * W_ * lifted_point;

                return (quadratic_term)[0] - 1.0;
            }

            void EnforceSymmetry(drake::solvers::MathematicalProgram& prog){
                prog.AddPositiveSemidefiniteConstraint(W_);
                for(int i = 0; i < N_q+1; ++i){
                    for(int j = 0; j < N_q+1; ++j){
                        if(i != j){
                            // Add constraints for symmetry
                            prog.AddLinearConstraint(W_(i, j) == W_(j, i)); // Symmetry constraint
                        }
                    }
                }
            }


            void EnforceOffsetBoxBounds(drake::solvers::MathematicalProgram& prog, double d_max){
                for(int i = 0; i < N_q; ++i){
                    prog.AddLinearConstraint(W_(N_q,i) <= d_max);
                    prog.AddLinearConstraint(W_(N_q,i) >= -d_max);
                }
            }

            void EnforceOffsetNorm(drake::solvers::MathematicalProgram& prog, double d_max){
                // Enforce that the norm of d is less than d_max
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                prog.AddLinearConstraint(d.transpose()*d <= d_max*d_max);
            }

            void ApproximatelyEnforceOffsetNorm(drake::solvers::MathematicalProgram& prog, double d_max, Eigen::VectorXd d_hat){
                // Enforce that the norm of d is less than d_max
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                prog.AddLinearConstraint((d_hat.transpose() * d) <= d_max);
            }
                
            // void AddCostOnOffset(drake::solvers::MathematicalProgram& prog, Eigen::VectorXd d_cost){
            //     // Add cost on d
            //     Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
            //     prog.AddCost(-d_cost.transpose() * d);
            // }

            void ApproximatelyMaximizeRadius(drake::solvers::MathematicalProgram& prog, Eigen::VectorXd d_hat){
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                prog.AddLinearCost( -(1 + d_hat.transpose() * d - f) );
            }

            void MaximizeRadius(drake::solvers::MathematicalProgram& prog){
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                prog.AddCost( -(d.transpose() * d - f + 1.0) );
            }

            void EnforceDiagonality(drake::solvers::MathematicalProgram& prog){
                for(int i = 0; i < N_q; ++i){
                    for(int j = 0; j < N_q; ++j){
                        if(i != j){
                            prog.AddLinearConstraint(W_(i, j) == 0);
                        }
                    }
                }
            }

            void EnforceMaxEigRatio(drake::solvers::MathematicalProgram& prog, double rho){
                for(int i = 0; i < N_q; ++i){
                    prog.AddLinearConstraint(W_(i,i) <= 1);
                    prog.AddLinearConstraint(W_(i,i) >= 1 / rho);
                }
            }

            void EnforceIdentity(drake::solvers::MathematicalProgram& prog){
                for(int i = 0; i < N_q; ++i){
                    prog.AddLinearConstraint(W_(i,i) == 1);
                }
                EnforceDiagonality(prog);
            }

            void ApproximatelyEnforceMinRadius(drake::solvers::MathematicalProgram& prog, double r_min, Eigen::VectorX<T> d_hat){
                // Enforce that the ellipsoid is at least r_min away from the reference point
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                prog.AddLinearConstraint( 1 + (d_hat.transpose() * d) - f >= r_min*r_min );
            }

            void EnforceConicMembershipOfOffset(drake::solvers::MathematicalProgram& prog, Eigen::VectorX<T> d_hat, double conic_rad){
                // Enforce that the offset is in the conic hull of the ellipsoid
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                // hmm how would you do this? I think I have an idea but
            }

            void EnforceNonlinearMinRadius(drake::solvers::MathematicalProgram& prog, double r_min){
                // Enforce that the ellipsoid is at least r_min away from the reference point
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                prog.AddConstraint( 1 + (d.transpose()*d) - f >= r_min*r_min );
            }
                        
            void AddNonMembershipConstraint(drake::solvers::MathematicalProgram& prog, Eigen::VectorXd point){
                prog.AddLinearConstraint(evaluate(point) >= 0);
            }

            void AddMembershipConstraint(drake::solvers::MathematicalProgram& prog, Eigen::VectorX<drake::symbolic::Expression> point){
                prog.AddConstraint(evaluate(point) <= 0);
            }

            void AddBoundaryConstraint(drake::solvers::MathematicalProgram& prog, Eigen::VectorX<drake::symbolic::Expression> point){
                prog.AddConstraint(evaluate(point) <= 0);
            }

            void AddMembershipConstraint(drake::solvers::MathematicalProgram& prog, Eigen::VectorXd point){
                prog.AddLinearConstraint(evaluate(point) <= 0);
            }
            
            void AddRadiusCost(drake::solvers::MathematicalProgram& prog, double r_f){
                prog.AddLinearCost(r_f * W_(N_q, N_q));
            }

            Eigen::VectorX<T> center() const {
                // Calculate the center of the ellipsoid
            
                Eigen::MatrixX<T> E = W_.topLeftCorner(N_q, N_q);
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);

                Eigen::MatrixX<T> E_inv = E.inverse();

                Eigen::VectorX<T> q_tilde = reference_point_ + (E_inv * d);

                return q_tilde;
            }

            std::vector<Eigen::VectorX<T>> sample(int N_samples) {
                // Sample points within the ellipsoid
                Eigen::VectorX<T> q_tilde = center();

                drake::log()->info("projecting to W: \n{}" , W_);

                auto E = W_.topLeftCorner(N_q, N_q);
                Eigen::VectorX<T> d = W_.topRightCorner(N_q, 1);
                T f = W_(N_q, N_q);

                T radius = d.transpose() * E.inverse() * d - f + 1.0;
                radius = sqrt(radius);

                drake::log()->info("radius: {}", radius);

                
                std::vector<Eigen::VectorX<T>> samples(N_samples);
                for (int i = 0; i < N_samples; ++i) {
                    Eigen::VectorX<T> random_point = Eigen::VectorX<T>::Random(N_q);
                    random_point.normalize();

                    random_point = random_point * 5.0; // Scale the random point

                    drake::solvers::MathematicalProgram proj_prog;
                    auto proj_sol_vars = proj_prog.NewContinuousVariables(N_q, "proj_sol_vars");

                    Eigen::VectorX<drake::symbolic::Expression> proj_sol_expr = proj_sol_vars;
                    AddBoundaryConstraint(proj_prog, proj_sol_expr);
                    
                    proj_prog.AddQuadraticCost(
                        (proj_sol_vars - random_point).squaredNorm());

                    auto result = drake::solvers::Solve(proj_prog);
                    if (!result.is_success()) {
                        // drake::log()->warn("Projection failed");
                        samples[i] = (q_tilde);
                        continue;
                    }
                    Eigen::VectorX<T> projected_point = result.GetSolution(proj_sol_vars);
                    
                    samples[i] = projected_point;
                }
                return samples;
                
            }

    };

    InflatedEllipsoid<double> GenerateInflatedEllipsoid(
        const InflatedEllipsoidOptions& options) {

        int N_q = options.N_q;
        auto point_cloud = options.point_cloud; // Each col is a point
        auto contain_points = options.contain_points; // Each col is a point
        auto ref_point = options.reference_point; // Reference point for the ellipsoid
        auto d_hat = options.d_hat; // Direction of the ellipsoid

        auto d_max = options.d_max; // Maximum offset for the ellipsoid

        double r_min = options.r_min; // Minimum radius of the ellipsoid

        drake::solvers::MathematicalProgram prog;
        Eigen::MatrixX<drake::symbolic::Variable> W_var = prog.NewContinuousVariables(N_q+1, N_q+1, "W");
        Eigen::MatrixX<drake::symbolic::Expression> W = W_var;
        InflatedEllipsoid<drake::symbolic::Expression> cal_C(W, ref_point);

        // drake::symbolic::Variable phi = prog.NewContinuousVariables(1, 1, "phi")(0,0);

        cal_C.EnforceIdentity(prog);
        cal_C.ApproximatelyMaximizeRadius(prog, d_hat);
        // cal_C.EnforceOffsetNorm(prog, d_max);
        cal_C.ApproximatelyEnforceMinRadius(prog, r_min, d_hat);
        cal_C.ApproximatelyEnforceOffsetNorm(prog, d_max, d_hat);
        cal_C.EnforceSymmetry(prog);
        

        // For now, enforce the exact directionality of the ellipsoid.
        // prog.AddLinearConstraint(W_var.topRightCorner(N_q, 1) == phi * d_hat);
        // prog.AddLinearConstraint(phi >= r_min);

        // Add non-membership constraints for each point in the point cloud
        for (int p_index = 0; p_index < options.point_cloud.cols(); ++p_index) {
            Eigen::VectorXd point = options.point_cloud.col(p_index);
            // Non-membership constraint
            cal_C.AddNonMembershipConstraint(prog, point);
        }

        for(int p_index = 0; p_index < options.contain_points.cols(); ++p_index) {
            Eigen::VectorXd point = options.contain_points.col(p_index);
            // membership constraint
            cal_C.AddMembershipConstraint(prog, point);
        }

        auto result = drake::solvers::Solve(prog);

        if (!result.is_success()) {
            throw std::runtime_error("Ellipsoid generation failed");
        }

        // Extract the optimized values
        Eigen::MatrixX<double> W_opt = result.GetSolution(W_var);

        // Create the ellipsoid object
        InflatedEllipsoid<double> ellipsoid(W_opt, ref_point);
        drake::log()->info("Ellipsoid generated successfully");

        drake::log()->info("W: \n{}", W_opt);

        auto center = ellipsoid.center();
        auto evaluate_center = ellipsoid.evaluate(center);
        drake::log()->info("center: \n{}", center.transpose());
        drake::log()->info("evaluate center: {}", evaluate_center);

        return ellipsoid;
    }
}   // namespace CorrGen