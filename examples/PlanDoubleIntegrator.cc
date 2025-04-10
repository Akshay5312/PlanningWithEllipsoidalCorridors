
#include "../BallalCorridor/BallCorridor.h"
#include "../CorridorPlanning/CorridorPlanningProblem.h"
#include "../CorridorPlanning/TimeVaryingProjectionToBall.h"

#include <drake/common/text_logging.h>
#include <drake/solvers/mathematical_program.h>

#include <iostream>
#include <fstream>

int main(){
    using namespace CorrGen;
    using namespace CorrPlanning;
    
    int Ny = 20;
    int N_samples= 20;

    std::string current_folder = __FILE__;
    current_folder = current_folder.substr(0, current_folder.find_last_of("/"));
    std::string corridor_file = current_folder + "/result_corridors/corridor_1.txt";

    // Read the corridor from a file
    std::ifstream file(corridor_file);
    std::string corridor_string;
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            corridor_string += line + "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << corridor_file << std::endl;
    }
    // Parse the corridor string to create a BallCorridor object
    auto corridor = BallCorridor<double>::FromString(corridor_string);

    // Create the linearized robot dynamics
    LinearizedRobotDynamics dynamics;
    
    // DOUBLE INTEGRATOR DYNAMICS
    dynamics.A = Eigen::MatrixXd::Zero(4, 4);
    dynamics.A(0, 2) = 1.0;
    dynamics.A(1, 3) = 1.0;

    dynamics.B = Eigen::MatrixXd::Zero(4, 2);
    dynamics.B(2, 0) = 1.0;
    dynamics.B(3, 1) = 1.0;

    dynamics.R = Eigen::MatrixXd::Identity(4, 4); // Cost on the input
    dynamics.P = Eigen::MatrixXd::Identity(4, 4); // Cost on the offset from the state
    dynamics.R(0, 0) = 0.0;
    dynamics.R(1, 1) = 0.0;
    dynamics.R(2, 2) = 100.0;
    dynamics.R(3, 3) = 100.0;


    // Create the corridor planning problem
    CorridorPlanningProblem planning_problem(corridor, dynamics, N_samples, Ny);

    // Solve the planning problem

    auto start_time = std::chrono::high_resolution_clock::now();
    auto approximate_solution = planning_problem.ApproximateSolutionByProjection();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    drake::log()->info("Solving took {} milliseconds", elapsed_time.count() * 1000);
    drake::log()->info("Approximate solution: \n{}", approximate_solution.transpose());
    drake::log()->info("Approximate solution norm: {}", approximate_solution.norm());


    std::string result_path_file = current_folder + "/sample_result_data/path_1.txt";
    std::ofstream result_file(result_path_file);
    for(int i = 0; i < N_samples; i++){
        double eps = static_cast<double>(i) / (N_samples-1);
        result_file << planning_problem.x(eps, approximate_solution).head(2).transpose() << "\n";
    }

}