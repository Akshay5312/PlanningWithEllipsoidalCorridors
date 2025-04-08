#include <iostream>
#include "../BallalCorridor/Ball.h"
#include "../BallalCorridor/GenerateBall.h"
#include "../BallalCorridor/GenerateBallCorridor.h"
#include "helpers/read_point_cloud.h"

int main(){
    using namespace CorrGen;
    int Nq;
    int N_ctrl = 5;
    int N_samples = 100;

    // Get the current folder location
    std::string current_folder = __FILE__;
    current_folder = current_folder.substr(0, current_folder.find_last_of("/"));

    // Read the point cloud data provided
    std::string point_cloud_file = current_folder + "/sample_point_cloud_data/point_cloud_1.txt";
    Eigen::MatrixXd point_cloud = readPointCloud(point_cloud_file, Nq);

    // Let us plan from 0,0,... to 1,1,...
    // The point cloud data can be generated in the box 0 < q < 1
    Eigen::VectorXd start_point(2);
    start_point << 0.1, 0.1;
    Eigen::VectorXd end_point(2);
    end_point << 0.9, 0.9;

    LagrangePolynomial<double> reference_path({0, 1}, {start_point, end_point});
    Eigen::VectorXd d_hat_reference(2);
    d_hat_reference << -1.0, 1.0;
    d_hat_reference = d_hat_reference.normalized();
    LagrangePolynomial<double> d_hat_path({0, 1}, {d_hat_reference, d_hat_reference});

    CorrGen::BallCorridorGenerationOptions corridor_options;
    corridor_options.point_cloud = point_cloud;
    corridor_options.reference_path = reference_path;
    corridor_options.N_q = Nq;
    corridor_options.N_samples = N_samples;
    corridor_options.N_ctrl = N_ctrl;
    corridor_options.enforce_taper = false; // Enforce the corridor to be 'tapered' s.t. it contains only the initial and final points at the extremes

    corridor_options.SolveLinearApproximation = true;
    corridor_options.approx_options.d_hat = d_hat_path;
    corridor_options.approx_options.d_max = 0.750; // Maximum distance from the center to the reference in any direction
    corridor_options.approx_options.r_min = 0.07; // Minimum radius
    corridor_options.approx_options.conic_factor = 0.1; // The conic factor

    // Create the corridor.
    BallCorridor<double> corridor = GenerateBallCorridor(corridor_options);

    // The corridor is now created.

    // Save the corridor to a file
    std::string corridor_file = current_folder + "/result_corridors/corridor_1.txt";
    std::ofstream file(corridor_file);
    if (file.is_open()) {
        file << corridor.toString();
        file.close();
    } else {
        std::cerr << "Unable to open file: " << corridor_file << std::endl;
    }


    // //VISUALIZATION
    // //option 1. Meshcat
    // // Create a meshcat visualizer, and add the points in the point cloud to it. 
    // // Sample the surface and interior of all ellipsoids.
    // // Create an alphashape mesh. This is the corridor.
    // // Add the alphashape to meshcat to visualize it.

    // // option 2. Matlab. 
    // // Sample the corridor to find ellipsoids at times between 0 and 1.

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
        
        auto d_hat_eps = corridor_options.approx_options.d_hat.value(eps).normalized();

        drake::log()->info("radius: \n{}", sqrt(ball_eps.squared_radius()));
        drake::log()->info("APPROXIMATE RADIUS: \n{}", ball_eps.approximate_squared_radius(d_hat_eps));
        drake::log()->info("center: \n{}", ball_eps.center().transpose());
    }

    std::string ellipsoid_cloud_file = current_folder + "/sample_result_data/corridor_1.txt";
    savePointCloud(ellipsoid_cloud_file, corridor_samples);

    

}