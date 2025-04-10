#pragma once

#include "../BallalCorridor/BallCorridor.h"
#include "TimeVaryingProjectionToBall.h"

namespace CorrPlanning{
    using namespace CorrGen;

    struct LinearizedRobotDynamics{
        Eigen::MatrixXd A; // Linearized dynamics matrix
        Eigen::MatrixXd B; // Linearized input matrix

        Eigen::MatrixXd R; // Cost on the input
        Eigen::MatrixXd P; // Cost on the offset from the state
    };

    class CorridorPlanningProblem{

        public:
        CorridorPlanningProblem(
            BallCorridor<double>& corridor,
            LinearizedRobotDynamics dynamics,
            int N_samples = 20,
            int N_y = 20
        ) : N_q_(corridor.d(0).rows()), 
            N_y_(N_y), 
            N_samples_(N_samples), 
            corridor_(corridor), 
            dynamics_(dynamics),
            projection_(N_y, N_q_, 4, N_samples)
            {
                int N_x = N_q_*2;
                Eigen::MatrixXd A = dynamics_.A;
                Eigen::MatrixXd R = dynamics_.R;
                if(R.rows() != N_x || R.cols() != N_x){
                    drake::log()->warn("R is not square, setting R to identity");
                    R = Eigen::MatrixXd::Identity(N_x, N_x);
                }

                sqrt_Q_ = Eigen::MatrixXd::Identity(N_x*N_samples, N_y);
                b_ = Eigen::VectorXd::Zero(N_x*N_samples);
                // we want to minimize ||R*(A*(r*\bar{L}*(y) + [q_tilde;d_q_tilde]) - (r*\bar{dL}*y+[d_q_tilde;dd_q_tilde])) ||^2 = ||sqrt(Q)*y + b||. 
                // where r*\bar{L}*y + [q_tilde;0] is the state of the system, r*L*y + q_tilde is the configuration, and \bar{L} = [L ;dL]
                // q_tilde is the center point of the ball corridor, r is the radius of the ball corridor. 
                // L projects a the unit ball in Ny to a unit ball in Nx, centered at the origin. Thus, r*L*y+q_tilde projects to a point in the scaled/offset ballal corridor.
                for(int i = 0; i < N_samples_; i++){
                    double eps = static_cast<double>(i) / (N_samples_-1);
                    Eigen::VectorXd q_tilde = corridor_.center(eps);
                    Eigen::VectorXd d_q_tilde = corridor_.d_center(eps);
                    Eigen::VectorXd dd_q_tilde = corridor_.dd_center(eps);

                    double radius = sqrt(corridor_.getBall(eps).squared_radius());
                    Eigen::MatrixXd L = radius*projection_.L(eps);
                    Eigen::MatrixXd dL = radius*projection_.dL(eps);
                    Eigen::MatrixXd ddL = radius*projection_.ddL(eps);  
                    assert(L.rows() == N_q_);
                    assert(L.cols() == N_y_);
                    assert(q_tilde.rows() == N_q_);

                    Eigen::MatrixXd L_bar = Eigen::MatrixXd(N_x, N_y_);
                    L_bar.block(0, 0, N_q_, N_y_) = L;
                    L_bar.block(N_q_, 0, N_q_, N_y_) = dL;

                    Eigen::MatrixXd dL_bar = Eigen::MatrixXd(N_x, N_y_);
                    dL_bar.block(0, 0, N_q_, N_y_) = dL;
                    dL_bar.block(N_q_, 0, N_q_, N_y_) = ddL;
                    
                    drake::log()->info("L_bar: \n{}", L_bar);
                    drake::log()->info("dL_bar: \n{}", dL_bar);

                    Eigen::MatrixXd sqrt_Q_i = R*(A*L_bar - dL_bar);
                    sqrt_Q_.block(i*N_x, 0, N_x, N_y) = sqrt_Q_i;
                    Eigen::VectorXd bar_q_tilde = Eigen::VectorXd::Zero(N_x);
                    bar_q_tilde.head(N_q_) = q_tilde;
                    bar_q_tilde.tail(N_q_) = d_q_tilde;

                    Eigen::VectorXd bar_d_q_tilde = Eigen::VectorXd::Zero(N_x);
                    bar_d_q_tilde.head(N_q_) = d_q_tilde;
                    bar_d_q_tilde.tail(N_q_) = dd_q_tilde;

                    Eigen::VectorXd b_i = (A*bar_q_tilde - bar_d_q_tilde);

                    b_.segment(i*N_x, N_x) = b_i;
                }

                drake::log()->info("sqrt_Q_: \n{}", sqrt_Q_);
                drake::log()->info("b_: \n{}", b_);

            }

        Eigen::VectorXd x(double eps, Eigen::VectorXd y){
            // x = r*L*y + q_tilde
            Eigen::VectorXd q_tilde = corridor_.center(eps);
            double radius = sqrt(corridor_.getBall(eps).squared_radius());
            Eigen::MatrixXd L = radius*projection_.L(eps);
            Eigen::VectorXd q = L*y + q_tilde;

            Eigen::MatrixXd dL = radius*projection_.dL(eps);
            Eigen::VectorXd dq = radius*dL*y;

            Eigen::VectorXd x = Eigen::VectorXd::Zero(N_q_*2);
            x.head(N_q_) = q;
            x.tail(N_q_) = dq;

            return x;
        }

        // void CreateQcQp(drake::solvers::MathematicalProgram& prog, const Eigen::VectorX<drake::symbolic::Variable>& y){
        //     prog.AddQuadraticCost( (sqrt_Q_*y + b_).squaredNorm() );
            
        //     prog.AddQuadraticConstraint(
        //         2*Eigen::MatrixXd::Identity(N_y_,N_y_), y, 
        //         Eigen::VectorX<drake::symbolic::Expression>::Zero(N_y_),
        //         Eigen::VectorX<drake::symbolic::Expression>::Ones(N_y_)
        //     );
        // }

        Eigen::VectorXd ApproximateSolutionByProjection(){
            // Find the pinverse of sqrt_Q_
            // drake::log()->info("sqrt_Q_: \n{}", sqrt_Q_);
            // drake::log()->info("b_: \n{}", b_);
            Eigen::MatrixXd sqrt_Q_pinv = sqrt_Q_.completeOrthogonalDecomposition().pseudoInverse();
            // drake::log()->info("sqrt_Q_pinv: \n{}", sqrt_Q_pinv);
            // Find the solution to the linear system
            Eigen::VectorXd y_opt = sqrt_Q_pinv * b_;

            double y_opt_norm = y_opt.norm();
               
            // Project this to a valid collision free solution
            Eigen::VectorXd y_proj = y_opt;
            if(y_opt_norm > 1.0){
                y_proj = y_opt / y_opt_norm;             
            }
            
            return y_proj;
        }

        private:

        int N_q_;
        int N_y_;
        int N_samples_;
        
        // Eigen::MatrixXd Q_; // The cost matrix for the dynamics
        Eigen::MatrixXd sqrt_Q_; // The cost matrix for the dynamics
        Eigen::VectorXd b_;

        BallCorridor<double> corridor_;
        LinearizedRobotDynamics dynamics_;
        VaryingProjectionToBall projection_;
    };
}