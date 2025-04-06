#include "../EllipsoidalCorridor/InflatedEllipsoid.h"
#include "helpers/read_point_cloud.h"



int main(){
    using namespace CorrGen;
    int Nq;
    // Get the current folder location
    std::string current_folder = __FILE__;
    current_folder = current_folder.substr(0, current_folder.find_last_of("/"));


    // Read the point cloud data provided
    std::string point_cloud_file = current_folder + "/sample_point_cloud_data/point_cloud_1.txt";
    Eigen::MatrixXd point_cloud = readPointCloud(point_cloud_file, Nq);
    Eigen::MatrixXd contain_points = Eigen::MatrixXd::Zero(Nq, 1);
    contain_points.col(0) << 0.75, 0.75;
    // contain_points.col(0) << 0.25, 0.75;

    // Let us plan from 0,0,... to 1,1,...
    // The point cloud data can be generated in the box 0 < q < 1
    Eigen::VectorXd reference_point(Nq);
    reference_point << 0.0, 0.0;

    InflatedEllipsoidOptions ellipsoid_options;
    ellipsoid_options.point_cloud = point_cloud;
    ellipsoid_options.reference_point = reference_point;
    ellipsoid_options.contain_points = contain_points;
    ellipsoid_options.N_q = Nq;
    ellipsoid_options.d_hat = Eigen::VectorXd::Ones(Nq);
    ellipsoid_options.d_hat = ellipsoid_options.d_hat.normalized();
    ellipsoid_options.d_max = 1.0; // The maximum offset for the ellipsoid

    
    // Create the corridor.
    InflatedEllipsoid<double> cal_C = GenerateInflatedEllipsoid(ellipsoid_options);

    // // The corridor is now created.

    // auto ellipsoid_samples = corridor.sample(1500);


    // //VISUALIZATION
    // //option 1. Meshcat
    // // Create a meshcat visualizer, and add the points in the point cloud to it. 
    // // Sample the surface and interior of all ellipsoids.
    // // Create an alphashape mesh. This is the corridor.
    // // Add the alphashape to meshcat to visualize it.

    // // option 2. Matlab. 
    // // Sample the corridor to find ellipsoids at times between 0 and 1.
    // // First, find samples in a ball.
    int N_ball_samples = 1000;

    auto ellipsoid_samples = cal_C.sample(N_ball_samples);

    auto E_eps = cal_C.W_.topLeftCorner(Nq, Nq);
    Eigen::VectorXd d_eps = cal_C.W_.topRightCorner(Nq, 1);
    double f_eps = cal_C.W_(Nq, Nq);

    drake::log()->info("E_eps: \n{}", E_eps);
    drake::log()->info("d_eps: \n{}", d_eps.transpose());
    drake::log()->info("d_hat: \n{}", ellipsoid_options.d_hat.transpose());
    drake::log()->info("f_eps: {}", f_eps);
    drake::log()->info("W: \n{}", cal_C.W_);

    auto approx_radius = sqrt(ellipsoid_options.d_hat.transpose() * d_eps - f_eps + 1.0);
    drake::log()->info("approx radius: {}", approx_radius);

    auto actual_radius = sqrt(d_eps.transpose() * d_eps - f_eps + 1.0);
    drake::log()->info("actual radius: {}", actual_radius);

    // Save the points to an "ellipsoid cloud" file.
    std::string ellipsoid_cloud_file = current_folder + "/sample_result_data/ellipsoid_1.txt";
    savePointCloud(ellipsoid_cloud_file, ellipsoid_samples);

    // for( int k = 0; k < point_cloud.cols(); k++){
    //     drake::log()->info("point cloud: {}", point_cloud.col(k).transpose());
        
    //     drake::log()->info("point eval: {}", corridor.evaluate( (Eigen::VectorXd)point_cloud.col(k) ));
    // }

}