#include "../BallalCorridor/BallCorridor.h"
#include "../BallalCorridor/GenerateBallCorridor.h"
#include "TimeVaryingProjectionToBall.h"
#include "CorridorPlanningProblem.h"

#include <fstream>


namespace CorrPlanning{
    using namespace CorrGen;

    struct ProblemParameters{
        LinearizedRobotDynamics dynamics;
        Eigen::VectorXd q0;
        Eigen::VectorXd qf;
        Eigen::MatrixXd collision_points;
        int N_ctr;
        int N_samples;
        int N_y;
        int NL_ctrl;

        double r_min;
        double d_max;
        double conic_factor;
    };
    
    class IterativePlanningProblem{
        public:
        IterativePlanningProblem(
            ProblemParameters& params,
            int N_corridors = 10
        ) : N_corridors_(N_corridors),
            params_(params),
            L_(params.N_y, params.q0.rows(), params.NL_ctrl, params.NL_ctrl)
        {
            options_.N_q = params.q0.rows();
            options_.N_ctrl = params.N_ctr;
            options_.N_samples = params.N_samples;
            // options_.point_cloud = Eigen::MatrixXd::Zero(options_.N_q, params.collision_points.size());
            // for(int i = 0; i < params.collision_points.size(); i++){
            //     options_.point_cloud.col(i) = params.collision_points[i];
            // }
            options_.point_cloud = params.collision_points;
            options_.SolveLinearApproximation = false;
            options_.enforce_taper = true;

            options_.reference_path = LagrangePolynomial<double>(
                {0, 1}, 
                {params.q0, params.qf}
            );

            options_.approx_options.d_hat = LagrangePolynomial<double>(
                {0, 1}, 
                {Eigen::VectorXd::Zero(options_.N_q), 
                    Eigen::VectorXd::Zero(options_.N_q)}
            );

            while(corridors_.size() < N_corridors_){
                // Generate a random corridor
                // Start with a random d_hat
                auto random_d_hat = Eigen::VectorXd::Random(options_.N_q).normalized();
                options_.approx_options.d_hat = LagrangePolynomial<double>(
                    {0, 1}, 
                    {random_d_hat, random_d_hat}
                );


                options_.approx_options.d_max = params.d_max; // Maximum distance from the center to the reference in any direction
                options_.approx_options.r_min = params.r_min; // Minimum radius
                options_.approx_options.conic_factor = params.conic_factor;

                // Generate the corridor
                bool is_success = false;
                auto corridor =  GenerateBallCorridor(options_, &is_success);
                if(is_success){
                    corridors_.push_back(corridor);
                    drake::log()->info("corridor {} found.", corridors_.size());
                }else{
                    drake::log()->info("failed.");
                    continue;
                }
            }
        }

        // Solve the planning problem
        Eigen::VectorXd ApproximateInEachCorridor(int& solution_corridor_idx){
            std::vector<Eigen::VectorXd> solutions_(N_corridors_);
            std::vector<double> costs_(N_corridors_);
            // Iterate over the corridors
            int min_cost_index = 0;
            double min_cost = std::numeric_limits<double>::max();
            for(int i = 0; i < N_corridors_; i++){
                auto corridor = corridors_[i];
                // Solve the planning problem
                auto planning_problem = CorridorPlanningProblem(
                    corridor, 
                    params_.dynamics, 
                    L_,
                    params_.N_samples
                );
                // Solve the planning problem

                auto approximate_solution = planning_problem.ApproximateSolutionByProjection();
                auto cost = planning_problem.solution_cost(approximate_solution);
                // Store the solution
                solutions_[i] = approximate_solution;
                costs_[i] = cost;
                if(cost < min_cost){
                    min_cost = cost;
                    min_cost_index = i;
                }
            }

            solution_corridor_idx = min_cost_index;

            // Return the best solution
            return solutions_[min_cost_index];
        }

        void SaveAllCorridorsToFile(std::string file_name){
            std::ofstream file(file_name);
            if (file.is_open()) {
                for(int i = 0; i < corridors_.size(); i++){
                    file << corridors_[i].toStringOfSamples();
                }
                file.close();
            } else {
                // std::cerr << "Unable to open file: " << file_name << std::endl;
            }
        }

        std::vector<Eigen::VectorXd> SampleSolution(Eigen::VectorXd y, int corridor_idx){
            auto corridor = corridors_[corridor_idx];
            std::vector<Eigen::VectorXd>  sol_traj(params_.N_samples);
            for(int i = 0; i < params_.N_samples; i++){
                double eps = static_cast<double>(i)/(params_.N_samples);
                auto ball = corridor.getBall(eps);
                double r = sqrt(ball.squared_radius());
                auto L = r*L_.L(eps);
                Eigen::VectorXd q = r*L*y + ball.center();

                drake::log()->info("output q({}): \n{}",eps, q.transpose());
                
                sol_traj[i] = q;
            }

            return sol_traj;
        }

        private:

        int N_corridors_;
        BallCorridorGenerationOptions options_;
        VaryingProjectionToBall L_;
        ProblemParameters params_;
        std::vector<BallCorridor<double>> corridors_;
    };

}