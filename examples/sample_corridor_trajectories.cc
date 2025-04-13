


#include "../BallalCorridor/BallCorridor.h"
#include "../CorridorPlanning/CorridorPlanningProblem.h"
#include "../CorridorPlanning/TimeVaryingProjectionToBall.h"

#include <drake/common/text_logging.h>
#include <drake/solvers/mathematical_program.h>
#include "helpers/read_point_cloud.h"

#include <iostream>
#include <fstream>


int main(){
    using namespace CorrGen;
    using namespace CorrPlanning;
    
    int Ny = 10;
    int N_samples= 100;
    int N_b_samples = 50;
    
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

    int Nq=corridor.d(0).rows();
    VaryingProjectionToBall projection(Ny, Nq, 10,10);

    std::vector<Eigen::VectorXd> ball_samples(N_b_samples);
    std::vector<std::vector<Eigen::VectorXd>> paths(N_b_samples);
    for(int i = 0; i < N_b_samples; i++){

        if(i < Ny){
            ball_samples[i] = Eigen::VectorXd::Zero(Ny);
            ball_samples[i][i] = -1.0;
        }else{
            if(i < 2*Ny){
                ball_samples[i] = Eigen::VectorXd::Zero(Ny);
                ball_samples[i][i-Ny] = 1.0;
            }else{
                ball_samples[i] = 2*(Eigen::VectorXd::Random(Ny) + Eigen::VectorXd::Random(Ny) + Eigen::VectorXd::Random(Ny) + Eigen::VectorXd::Random(Ny)) - 4*Eigen::VectorXd::Ones(Ny);
                ball_samples[i] = ball_samples[i]/ball_samples[i].norm();
            }
        }
        std::vector<Eigen::VectorXd> path_samples(N_samples);

        for(int j = 0; j < N_samples; j++){
            double eps = static_cast<double>(j) / (N_samples-1);

            Eigen::VectorXd q_tilde = corridor.center(eps);
            double radius = sqrt(corridor.getBall(eps).squared_radius());
            Eigen::MatrixXd L = radius*projection.L(eps);
            Eigen::VectorXd q = L*ball_samples[i] + q_tilde;

            path_samples[j]=(q);
        }
        paths[i]=(path_samples);
    }


    std::string result_paths_file = current_folder + "/sample_result_data/paths_in_1.txt";
    std::ofstream result_file(result_paths_file);
    for(int i = 0; i < N_b_samples; i++){
        for(int j = 0; j < N_samples; j++){
            double eps = static_cast<double>(i) / (N_samples-1);
            result_file << paths[i][j].transpose() << "\n";
        }
    }

    std::vector<double> sample_times(N_samples);
    for (int k = 0; k < N_samples; ++k) {
        sample_times[k] = static_cast<double>(k) / (N_samples-1); // Example control breaks
    }    
    
    std::vector<Eigen::VectorXd> corridor_samples;
    for(int k = 0; k < N_samples; k++){
        double eps = sample_times[k]; // sample breaks

        auto ball_eps = corridor.getBall(eps);

        auto samples_at_eps = ball_eps.sample(1000);
        for(int i = 0; i < samples_at_eps.size(); ++i) {
            corridor_samples.push_back(samples_at_eps[i]);
        }

        drake::log()->info("eps: {}", eps);
        
        drake::log()->info("radius: \n{}", sqrt(ball_eps.squared_radius()));
        drake::log()->info("center: \n{}", ball_eps.center().transpose());
    }

    std::string ellipsoid_cloud_file = current_folder + "/sample_result_data/corridor_1.txt";
    savePointCloud(ellipsoid_cloud_file, corridor_samples);


}