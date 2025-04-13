#pragma once
#include "Ball.h"
#include "BallCorridor.h"

namespace CorrGen{
struct BallCorridorGenerationOptions{
    int N_q;
    
    Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor
    
    LagrangePolynomial<double> reference_path; // The reference path for the corridor

    int N_ctrl; // The number of control points for the corridor
    int N_samples; // The number of times to sample the corridor.

    bool SolveLinearApproximation = false;

    bool enforce_taper = false; // Enforce the corridor to be 'tapered' s.t. it contains only the initial and final points at the extremes.

    struct LinearApproximationOptions{
        LagrangePolynomial<double> d_hat; // The direction of offset between the center of the ball and the reference point

        double d_max; // The maximum offset for the ball

        double r_min = 0.05; // The minimum radius of the ball

        // We enforce the center of the ball to be within some 'cone' from the offset, centered at d_hat
        // This factor relates to the 'width' of the cone. 
        double conic_factor = 1e-1; 
    };

    LinearApproximationOptions approx_options;
};

BallCorridor<double> GenerateBallCorridor(const BallCorridorGenerationOptions& options, bool* success = nullptr){
    drake::solvers::MathematicalProgram prog;
    
    int N_q = options.N_q;
    int N_ctrl = options.N_ctrl;
    int N_samples = options.N_samples;

    std::vector<Eigen::VectorX<drake::symbolic::Variable>> d_vars(N_ctrl);
    std::vector<drake::symbolic::Variable> f_vars(N_ctrl);
    std::vector<double> break_times(N_ctrl);
    for(int i = 0; i < N_ctrl; i++){
        d_vars[i] = prog.NewContinuousVariables(N_q, 1);
        f_vars[i] = prog.NewContinuousVariables(1, 1)(0,0);
        break_times[i] = static_cast<double>(i) / (N_ctrl-1); //s control breaks
    }

    std::vector<Eigen::MatrixX<drake::symbolic::Expression>> d_exprs(N_ctrl);
    std::vector<Eigen::MatrixX<drake::symbolic::Expression>> f_exprs(N_ctrl);
    for(int i = 0; i < N_ctrl; i++){
        d_exprs[i] = d_vars[i];
        f_exprs[i] = Eigen::MatrixX<drake::symbolic::Expression>(1,1);
        f_exprs[i](0,0) = f_vars[i];
    }

    BallCorridor<drake::symbolic::Expression> corridor(
        LagrangePolynomial<drake::symbolic::Expression>(break_times, d_exprs),
        LagrangePolynomial<drake::symbolic::Expression>(break_times, f_exprs),
        options.reference_path
    );

    // Sample the corridor at uniformly distant times
    std::vector<double> sample_breaks(N_samples);
    for(int k = 0; k < N_samples; k++){
        sample_breaks[k] = static_cast<double>(k) / (N_samples-1); // Example control breaks
    }
    std::vector<Ball<drake::symbolic::Expression>> balls;

    for(int k = 0; k < N_samples; k++){
        double eps = sample_breaks[k]; // sample breaks
        
        balls.push_back(corridor.getBall(eps));

        if(options.SolveLinearApproximation){
            auto d_hat = options.approx_options.d_hat.value(eps).normalized();
            auto d_max = options.approx_options.d_max;
            auto r_min = options.approx_options.r_min;
            auto conic_factor = options.approx_options.conic_factor;

            if((!options.enforce_taper) || (( eps >= break_times[1]) && ( eps <= break_times[N_ctrl-2]))) {         
                balls[k].ApproxEnforceMinRadius(prog, r_min, d_hat, (r_min + d_max)/2);   
            }

            balls[k].ApproxAddMaxRadiusCost(prog, d_hat, (r_min + d_max)/2);
            balls[k].ApproxEnforceMaxOffset(prog, d_hat, d_max);
            balls[k].EnforceConicDirectionality(prog, d_hat, conic_factor*d_max);
        }else{
            balls[k].EnforceMinRadius(prog, options.approx_options.r_min);
            balls[k].AddMaxRadiusCost(prog);
            balls[k].EnforceMaxOffset(prog, options.approx_options.d_max);
        }

        // Enforce the points to be outside the ball
        balls[k].EnforceNonMembership(prog, options.point_cloud);
        prog.AddLinearConstraint(
            balls[k].center() <= Eigen::VectorXd::Ones(N_q) * 1.0);
        prog.AddLinearConstraint(
            balls[k].center() >= Eigen::VectorXd::Ones(N_q) * 0.0);
    }
    
    if(options.enforce_taper){
        //   If we 'enforce taper', then we cannot enforce a 'minimum radius'
        prog.AddLinearConstraint(
            balls[0].d() == Eigen::VectorXd::Zero(N_q));
        prog.AddLinearConstraint(
            balls[0].f() == 1.0);

        prog.AddLinearConstraint(
            balls[N_samples-1].d() == Eigen::VectorXd::Zero(N_q));
        prog.AddLinearConstraint(
            balls[N_samples-1].f() == 1.0);
        // prog.AddLinearConstraint(
        //     balls[0].f() <= 1.0);
        // prog.AddLinearConstraint(
        //     balls[N_samples-1].f() <= 1.0);
    }else{
        // auto initial_q = options.reference_path.value(0);
        // auto final_q = options.reference_path.value(1);
        prog.AddLinearConstraint(
            balls[0].f() <= 1.0);
        prog.AddLinearConstraint(
            balls[N_samples-1].f() <= 1.0);
    }

    // Solve the program
    auto result = drake::solvers::Solve(prog);
    
    if (success != nullptr) {
        *success = result.is_success();
    }

    if (!result.is_success()) {
        drake::log()->warn("Ball corridor generation failed. Solver result : {}", result.get_solution_result());
        drake::log()->info("is success? {}", *success);
        return BallCorridor<double>();
    }

    std::vector<Eigen::MatrixX<double>> d(N_ctrl);
    std::vector<Eigen::MatrixX<double>> f(N_ctrl);

    for(int i = 0; i < N_ctrl; i++){
        d[i] = result.GetSolution(d_vars[i]);
        f[i] = Eigen::MatrixXd(1,1);
        f[i](0,0) =
            result.GetSolution(f_vars[i]);
    }

    auto retval = BallCorridor<double>(
        LagrangePolynomial<double>(break_times, d),
        LagrangePolynomial<double>(break_times, f),
        options.reference_path
    );

    return retval; // Return the ballal corridor with the center and radius
}



}