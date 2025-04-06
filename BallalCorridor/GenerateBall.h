#pragma once
#include "Ball.h"

namespace CorrGen{
struct BallGenerationOptions{
    int N_q;
    
    Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor
    
    Eigen::VectorXd reference_point; // The reference point for the ball

    bool SolveLinearApproximation = false;

    struct LinearApproximationOptions{
        Eigen::VectorXd d_hat; // The direction of offset between the center of the ball and the reference point

        double d_max; // The maximum offset for the ball

        double r_min = 0.05; // The minimum radius of the ball

        // We enforce the center of the ball to be within some 'cone' from the offset, centered at d_hat
        // This factor relates to the 'width' of the cone. 
        double conic_factor = 1e-1; 
    };

    LinearApproximationOptions approx_options;
};

Ball<double> GenerateInflatedBall(const BallGenerationOptions& options){
    drake::solvers::MathematicalProgram prog;
    
    int N_q = options.N_q;

    Eigen::VectorX<drake::symbolic::Variable> d_var = prog.NewContinuousVariables(N_q, 1);
    drake::symbolic::Variable f_var = prog.NewContinuousVariables(1, 1)(0,0);

    Eigen::VectorX<drake::symbolic::Expression> d_expr = d_var;
    drake::symbolic::Expression f_expr = f_var;

    Ball<drake::symbolic::Expression> ball_C(d_var, f_var, options.reference_point);

    if(options.SolveLinearApproximation){
        auto d_hat = options.approx_options.d_hat.normalized();
        auto d_max = options.approx_options.d_max;
        auto r_min = options.approx_options.r_min;
        auto conic_factor = options.approx_options.conic_factor;

        ball_C.ApproxEnforceMinRadius(prog, r_min, d_hat, d_max);
        ball_C.ApproxAddMaxRadiusCost(prog, d_hat, (r_min + d_max)/2);
        ball_C.ApproxEnforceMaxOffset(prog, d_hat, d_max);
        ball_C.EnforceConicDirectionality(prog, d_hat, conic_factor*d_max);
    }else{
        ball_C.EnforceMinRadius(prog, options.approx_options.r_min);
        ball_C.AddMaxRadiusCost(prog);
        ball_C.EnforceMaxOffset(prog, options.approx_options.d_max);
    }

    // Enforce the points to be outside the ball
    ball_C.EnforceNonMembership(prog, options.point_cloud);
    
    // Solve the program
    auto result = drake::solvers::Solve(prog);
    if (!result.is_success()) {
        drake::log()->warn("Ball generation failed. Solver result : {}", result.get_solution_result());
        return Ball<double>(Eigen::VectorXd::Zero(N_q), 0.0, Eigen::VectorXd::Zero(N_q));
    }

    Eigen::VectorX<double> d = result.GetSolution(d_var);
    double f = result.GetSolution(f_var);
    drake::log()->info("Ball center: {}", d.transpose());
    drake::log()->info("radius: {}", sqrt(d.transpose()*d - f + 1.0));

    // Return the ball
    return Ball<double>(d, f, options.reference_point); // Return the ball with the center and radius


}



}