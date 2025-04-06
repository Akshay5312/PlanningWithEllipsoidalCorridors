#include "../EllipsoidalCorridor/EllipsoidalCorridor.h"
#include "../EllipsoidalCorridor/InflatedEllipsoid.h"
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
    start_point << 0.0, 0.0;
    Eigen::VectorXd end_point(2);
    end_point << 1.0, 1.0;

    LagrangePolynomial<double> reference_path({0, 1}, {start_point, end_point});
    Eigen::VectorXd cost_vector = Eigen::VectorXd::Random(Nq).normalized();
    LagrangePolynomial<double> cost_path({0, 1}, {cost_vector, cost_vector});

    CorrGen::EllipsoidalCorridorOptions corridor_options;
    corridor_options.point_cloud = point_cloud;
    corridor_options.cost_path = cost_path;
    corridor_options.reference_path = reference_path;
    corridor_options.N_q = Nq;
    // corridor_options.r_diagonal = 10.0;
    corridor_options.N_samples = N_samples;
    corridor_options.N_ctrl = N_ctrl;
    corridor_options.r_f = 0.01; // Cost on the radius of the ellipsoid
    corridor_options.max_d = 1.00; // Maximum distance from the center of the ellipsoid to the reference in any direction
    corridor_options.r_min = 0.1; // Minimum radius of the ellipsoid

    // Create the corridor.
    EllipsoidalCorridor<double> corridor = GenerateEllipsoidalCorridor(corridor_options);

    // The corridor is now created.


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
        Eigen::MatrixXd W_eps = corridor.W(sample_times[k]);

        InflatedEllipsoid<double> cal_C(W_eps, corridor.reference_path_.value(sample_times[k]));

        auto samples_at_eps = cal_C.sample(1000);
        for(int i = 0; i < samples_at_eps.size(); ++i) {
            corridor_samples.push_back(samples_at_eps[i]);
        }
    }

    std::string ellipsoid_cloud_file = current_folder + "/sample_result_data/ellipsoid_corridor_1.txt";
    savePointCloud(ellipsoid_cloud_file, corridor_samples);

}