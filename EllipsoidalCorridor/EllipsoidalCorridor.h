#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>

#include "LagrangePolynomial.h"

namespace CorrGen{
    template<typename ReferencePathType>
    // This struct is used to define the options for calculating the Ellipsoidal Corridor
    struct EllipsoidalCorridorOptions{
        Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor
        ReferencePathType reference_path; // The initial guess for the path

        int N_q; // The number of dimensions for the point cloud.

        int num_ctr_pts = 5; // The number of control points for the B-spline trajectory
        int num_samples = 100; // The number of samples for the trajectory

        double r_diagonal = 10.0; // The diagonal dominance for the ellipsoid matrix
        double r_f = 1.0; // The diagonal dominance for the ellipsoid matrix


        // The following are NOT IMPLEMENTED, but may be nice to have.
        // Enforce that the corridor contains the reference path. 
        // This is done by constraining f_eps <= 1 \forall epsilon in [0, 1].
        bool contain_reference_path = false;
        // Enforce that the corridor "tapers" at the boundaries.
        bool taper = false; 
        // Constrain the maximum diffeence between elements defining two consequtive ellipsoids.
        double ellipsoidal_velocity_constraint = 1e6;

        double time_horizon = 1;

        // C^-1*E^-1*d approximates the offset of the center of the ellipsoid from the reference trajectory.
        // Assuming that E is significantly diagonally dominant, and the ellipsoid is approximately a ball (i.e. E \approx I),
        // then "d" dictates the direction in which the center of the ellipsoid is offset from the reference trajectory.
        // This can be exploited to create corridors with differing homotopy by adding some random weight to "d" in the cost. 
        // The central trajectory will be in (approximately) the opposing direction of the weight.
        bool randomize_homotopy = false;
    };


    template<typename T>
    class EllipsoidalCorridor{
        public:
            EllipsoidalCorridor(CorrGen::LagrangePolynomial<T> E, CorrGen::LagrangePolynomial<T> d, CorrGen::LagrangePolynomial<T> f, drake::trajectories::PiecewisePolynomial<double> reference_path)
                : E_(E), d_(d), f_(f), c_(reference_path) {
                    drake::log()->info("Ellipsoidal Corridor created");
                }

            // Accessors for the ellipsoidal corridor properties
            CorrGen::LagrangePolynomial<T> E() const { return E_; }
            CorrGen::LagrangePolynomial<T> d() const { return d_; }
            CorrGen::LagrangePolynomial<T> f() const { return f_; }
            drake::trajectories::PiecewisePolynomial<double> c() const { return c_; }

            // Evaluate the ellipsoidal corridor properties at a given epsilon
            // where epsilon is a parameter that varies along the path from 0 to 1
            Eigen::MatrixX<T> E(double epsilon) const { return E_.value(epsilon); }
        
            Eigen::VectorX<T> d(double epsilon) const { return d_.value(epsilon); }

            Eigen::VectorX<T> f(double epsilon) const { return f_.value(epsilon); }

            Eigen::VectorXd c(double epsilon) const { return c_.value(epsilon); }

            T evaluate(double epsilon, Eigen::VectorXd point){
                // Evaluates if 'point' is in the ellipsoidal cross section of the corridor at 'e'.
                // Corridor(e) = {q | (q-c_e)E_e(q-c_e) - 2*d_e(q-c_e) + f_e <= 1}
                auto o_point = point - c(epsilon);
                return (o_point.transpose()*E(epsilon)*o_point - 2*d(epsilon).transpose()*o_point) + f(epsilon)[0] - 1;
            }

            // This is the path that the center of the ellipsoid traces.
            Eigen::VectorX<T> ellipsoid_center(double epsilon) const { 
                // ONLY IMPLEMENT IF ISSAME<T, DOUBLE>
                auto E_eps = E(epsilon);
                auto d_eps = d(epsilon);
                auto c_eps = c(epsilon);

                // Get the inverse of E and C
                Eigen::MatrixXd E_inv = E_eps.inverse();
                // Calculate the center of the ellipsoid

                Eigen::VectorXd center = c_eps + (E_inv * d_eps);

                return center;
            }
            
        private:
            CorrGen::LagrangePolynomial<T> E_; // Ellipsoid coridor matrix path
            CorrGen::LagrangePolynomial<T> d_; // Ellipsoid corridor vector path
            CorrGen::LagrangePolynomial<T> f_; // Ellipsoid corridor scalar path
            drake::trajectories::PiecewisePolynomial<double> c_; // Reference path for the corridor

            // An ellipsoidal corridor is defined as 
            // Corridor(e) = {q | (q-c_e)E_e(q-c_e) - 2*d_e(q-c_e) + f_e <= 1}
    };



    EllipsoidalCorridor<double> GenerateEllipsoidalCorridor(
        const EllipsoidalCorridorOptions<drake::trajectories::PiecewisePolynomial<double>>& options){

        drake::log()->info("Generating Ellipsoidal Corridor");

        // Generate the ellipsoidal corridor based on the provided options
        // For now, let us solve the linear programming relaxation

        int N_samples = options.num_samples;
        int N_ctr_pts = options.num_ctr_pts;
        int N_q = options.N_q;
        Eigen::MatrixXd point_cloud = options.point_cloud; // Each col is a point
        drake::trajectories::PiecewisePolynomial<double> reference_path = options.reference_path;

        double r_diagonal = options.r_diagonal;
        double r_f = options.r_f;

        // Initialize the Linear Program
        drake::solvers::MathematicalProgram corridor_solver_lp;

        drake::log()->info("Creating the corridor solver");

        // Define the variables for the ellipsoidal corridor
        // E is the ellipsoid matrix, d is the vector, f is the scalar
        // We want to define control points for each path.

        std::vector<drake::solvers::MatrixXDecisionVariable> E_ctr_vars(N_ctr_pts);
        std::vector<drake::solvers::VectorXDecisionVariable> d_ctr_vars(N_ctr_pts);
        // std::vector<drake::solvers::VectorXDecisionVariable> f_ctr_vars(N_ctr_pts);
        
        std::vector<double> ctrl_breaks_dbl(N_ctr_pts);

        // Initialize the control points for the ellipsoidal corridor
        // These are the variables we want to optimize over
        for(int k = 0; k < N_ctr_pts; k++) {
            E_ctr_vars[k] = corridor_solver_lp.NewContinuousVariables(N_q, N_q, "E_" + std::to_string(k));
            d_ctr_vars[k] = corridor_solver_lp.NewContinuousVariables(N_q, 1, "d_" + std::to_string(k));
            // f_ctr_vars[k] = corridor_solver_lp.NewContinuousVariables(1, 1, "f_" + std::to_string(k));
            // ctrl_breaks[k] = static_cast<double>(k) / (N_ctr_pts - 1);
            ctrl_breaks_dbl[k] = static_cast<double>(k) / (N_ctr_pts - 1);
        }

        drake::log()->info("Control points initialized");

        std::vector<Eigen::MatrixX<drake::symbolic::Expression>> E_ctr_expr(N_ctr_pts);
        std::vector<Eigen::MatrixX<drake::symbolic::Expression>> d_ctr_expr(N_ctr_pts);
        // std::vector<Eigen::MatrixX<drake::symbolic::Expression>> f_ctr_expr(N_ctr_pts);

        for (int k = 0; k < N_ctr_pts; ++k) {
            E_ctr_expr[k] = E_ctr_vars[k];
            d_ctr_expr[k] = d_ctr_vars[k];
            // f_ctr_expr[k] = f_ctr_vars[k];
        }


        // Get the times to evaluate the ellipsoidal corridor at
        auto sample_breaks = std::vector<double>(N_samples);
        for (int k = 0; k < N_samples; k++) {
            sample_breaks[k] = static_cast<double>(k) / (N_samples - 1);
        }

        auto E_traj = 
            CorrGen::LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks_dbl, E_ctr_expr);
        auto d_traj =
            CorrGen::LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks_dbl, d_ctr_expr);
        // auto f_traj =
        //     CorrGen::LagrangePolynomial<drake::symbolic::Expression>(ctrl_breaks_dbl, f_ctr_expr);

        drake::log()->info("Trajectories created");


        // SET UP THE LP PROBLEM
        for(int k = 0; k < N_samples; k++){
            drake::log()->info("Sample: {}", k);
            double eps = sample_breaks[k];
            drake::log()->info("eps: {}", eps);
            
            // Get the ellipsoid matrix, vector, and scalar at this sample point
            auto E_eps = E_traj.value(eps);

            drake::log()->info("E_eps: \n{}", E_eps);
            
            Eigen::VectorX<drake::symbolic::Expression> d_eps = d_traj.value(eps);

            drake::log()->info("d_eps: \n{}", d_eps);
            
            // drake::symbolic::Expression f_eps = f_traj.value(eps)(0,0);

            // Add cost to the prog.
            corridor_solver_lp.AddLinearCost( E_eps.trace() );
        
            // Add non-membership constriant for each point
            for(int p_index = 0; p_index < point_cloud.cols(); p_index++){
                
                auto point = point_cloud.col(p_index);

                auto o_point = point - reference_path.value(eps);
                // drake::log()->info("o_point: {}", o_point);

                drake::symbolic::Expression non_membership = 
                    (o_point.transpose()*E_eps*o_point - 2*d_eps.transpose()*o_point) - 1;
                corridor_solver_lp.AddLinearConstraint(non_membership >= 0);
            }


            // Add symmetry constraint
            for(int i = 0; i < N_q; i++){
                for(int j = i+1; j < N_q; j++){
                    corridor_solver_lp.AddLinearEqualityConstraint( E_eps(i,j) == E_eps(j,i) );
                    // corridor_solver_lp.AddLinearConstraint( E_eps(i,j) <= 10000 );
                }
                corridor_solver_lp.AddLinearConstraint( E_eps(i,i) >= 0.1 );
                corridor_solver_lp.AddLinearConstraint( E_eps(i,i) <= 1e6 );
                corridor_solver_lp.AddLinearConstraint( d_eps(i) <= 0.1);
                corridor_solver_lp.AddLinearConstraint( d_eps(i) >= 0.1);
            }

            // corridor_solver_lp.AddLinearConstraint(f_eps == 0);

            // Add diagonal dominance
            for(int i = 0; i < N_q; i++){
                for(int j = i+1; j < N_q; j++){
                    corridor_solver_lp.AddLinearConstraint( E_eps(i,i) >= r_diagonal*E_eps(j,i) );
                    corridor_solver_lp.AddLinearConstraint( E_eps(i,i) >= -r_diagonal*E_eps(j,i) );
                    
                }
            }


        }

        // Restrict the initial and final ellipsoids to contain c_0, c_T. 
        // Let us do this (for now) by constraining f_0 = f_T = 0.
        // We can also state f_0=f_T=1, d_0=d_T=0. This would constrain the corridor to have point valued ellipsoids at eps=0,1.
        // corridor_solver_lp.AddLinearEqualityConstraint(f_traj.value(0)(0,0) == 0 );
        // corridor_solver_lp.AddLinearEqualityConstraint(f_traj.value(1)(0,0)== 0 );

        drake::log()->info("Constraints added. Solving");

        auto result = drake::solvers::Solve(corridor_solver_lp);

        drake::log()->info("Solved");

        if(!result.is_success()) {
            drake::log()->error("Failed to solve the corridor optimization problem. Status: {}", result.get_solution_result());
            throw std::runtime_error("Ellipsoidal Corridor generation failed.");
        }


        std::vector<Eigen::MatrixXd> E_ctr(N_ctr_pts);
        std::vector<Eigen::MatrixXd> d_ctr(N_ctr_pts);
        std::vector<Eigen::MatrixXd> f_ctr(N_ctr_pts);

        for(int k = 0; k < N_ctr_pts; ++k) {
            E_ctr[k] = result.GetSolution(E_ctr_vars[k]);
            d_ctr[k] = result.GetSolution(d_ctr_vars[k]);
            f_ctr[k] = Eigen::VectorXd::Zero(1);

            drake::log()->info("Control point {}: E = \n{}", k, E_ctr[k]);
            drake::log()->info("Control point {}: d = \n{}", k, d_ctr[k]);
            // drake::log()->info("Control point {}: f = \n{}", k, f_ctr[k]);

        }

        EllipsoidalCorridor<double> corridor_result(
            CorrGen::LagrangePolynomial(ctrl_breaks_dbl, E_ctr),
            CorrGen::LagrangePolynomial(ctrl_breaks_dbl, d_ctr),
            CorrGen::LagrangePolynomial(ctrl_breaks_dbl, f_ctr),
            reference_path
        );

        return corridor_result;
    }
};