
#include "../BallalCorridor/BallCorridor.h"
#include "../CorridorPlanning/CorridorPlanningProblem.h"
#include "../CorridorPlanning/TimeVaryingProjectionToBall.h"
#include "../CorridorPlanning/IterativelySolve.h"

#include <drake/common/text_logging.h>
#include <drake/solvers/mathematical_program.h>
#include "helpers/read_point_cloud.h"

#include <iostream>
#include <fstream>

int main(){
    using namespace CorrGen;
    using namespace CorrPlanning;
    
    int Ny = 50;
    int N_samples= 50;
    int N_ctrl = 5;
    int NL_ctrl = 4;

    int Nq;
    // Get the current folder location
    std::string current_folder = __FILE__;
    current_folder = current_folder.substr(0, current_folder.find_last_of("/"));

    // Read the point cloud data provided
    std::string point_cloud_file = current_folder + "/sample_point_cloud_data/point_cloud_1.txt";
    Eigen::MatrixXd point_cloud = readPointCloud(point_cloud_file, Nq);


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

    ProblemParameters params;
    params.dynamics = dynamics;
    params.N_ctr = N_ctrl;
    params.N_y = Ny;
    params.NL_ctrl = NL_ctrl;
    params.collision_points = point_cloud;

    params.N_samples = N_samples;
    params.q0 = Eigen::VectorXd::Zero(2);
    params.q0[0] = 0.1;
    params.q0[1] = 0.1;
    
    params.qf = Eigen::VectorXd::Zero(2);
    params.qf[0] = 0.9;
    params.qf[1] = 0.9;

    params.r_min = 0.08;
    params.d_max = 0.75;
    params.conic_factor = 0.1;

    IterativePlanningProblem planning_problem(params, 10);

    int solution_corridor_index = 0;
    Eigen::VectorXd y_sol = planning_problem.ApproximateInEachCorridor( solution_corridor_index );    

    std::string corridor_file = current_folder + "/result_corridors/corridor_1.txt";
    planning_problem.SaveAllCorridorsToFile(corridor_file);

    std::string trajectory_file = current_folder + "/result_corridors/trajectory_1.txt";
    std::vector<Eigen::VectorXd> trajectory = planning_problem.SampleSolution(y_sol, solution_corridor_index);
    savePointCloud(trajectory_file, trajectory);
    
}