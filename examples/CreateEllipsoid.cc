#include "../EllipsoidalCorridor/SingleEllipsoid.h"
#include "helpers/read_point_cloud.h"



int main(){
    int Nq;
    // Get the current folder location
    std::string current_folder = __FILE__;
    current_folder = current_folder.substr(0, current_folder.find_last_of("/"));


    // Read the point cloud data provided
    std::string point_cloud_file = current_folder + "/sample_point_cloud_data/point_cloud_1.txt";
    Eigen::MatrixXd point_cloud = readPointCloud(point_cloud_file, Nq);

    // Let us plan from 0,0,... to 1,1,...
    // The point cloud data can be generated in the box 0 < q < 1
    Eigen::VectorXd reference_point = Eigen::VectorXd::Ones(Nq) * 0.25; 

    InflatedEllipsoidOptions ellipsoid_options;
    ellipsoid_options.point_cloud = point_cloud;
    ellipsoid_options.reference_point = reference_point;
    ellipsoid_options.N_q = Nq;
    
    // Create the corridor.
    InflatedEllipsoid<double> corridor = GenerateInflatedEllipsoid(ellipsoid_options);

    // The corridor is now created.


    //VISUALIZATION
    //option 1. Meshcat
    // Create a meshcat visualizer, and add the points in the point cloud to it. 
    // Sample the surface and interior of all ellipsoids.
    // Create an alphashape mesh. This is the corridor.
    // Add the alphashape to meshcat to visualize it.

    // option 2. Matlab. 
    // Sample the corridor to find ellipsoids at times between 0 and 1.
    // First, find samples in a ball.
    int N_ball_samples = 1000;
    Eigen::MatrixXd ball_samples(Nq, N_ball_samples);
    for (int k = 0; k < N_ball_samples; ++k) {
        ball_samples.col(k) = Eigen::VectorXd::Random(Nq);
        ball_samples.col(k) = ball_samples.col(k) / ball_samples.col(k).norm();
        // ball_samples.col(k) = ball_samples.col(k) * Eigen::VectorXd::Random(1);
    }


    std::vector<Eigen::VectorXd> ellipsoid_samples(N_ball_samples);

    // Sample points within the ellipsoids
    auto E_eps = corridor.E_;
    auto d_eps = corridor.d_;
    auto f_eps = corridor.f_;
    
    drake::log()->info("E_eps: {}", E_eps);
    drake::log()->info("d_eps: {}", d_eps);
    drake::log()->info("f_eps: {}", f_eps);
    // Add cost to the prog.
    

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(E_eps);
    Eigen::MatrixXd C_eps = eigen_solver.operatorSqrt();
    drake::log()->info("C_eps: \n{}", C_eps);
    Eigen::MatrixXd C_eps_inv = C_eps.inverse();

    Eigen::VectorXd q_tilde = corridor.center();

    
    // drake::log()->info("q_tilde: \n{}", C_eps);

    for(int i = 0; i < N_ball_samples; ++i){
        double o = (
            f_eps + 
            (2*d_eps.transpose() * reference_point
            + reference_point.transpose() * E_eps * reference_point 
            - q_tilde.transpose() * E_eps * q_tilde)[0] );

        drake::log()->info("o: {}", o);

        ellipsoid_samples[i] = sqrt(1 - o) *(C_eps_inv * (ball_samples.col(i))) + q_tilde; 
    }

    // Save the points to an "ellipsoid cloud" file.
    std::string ellipsoid_cloud_file = current_folder + "/sample_result_data/ellipsoid_1.txt";
    savePointCloud(ellipsoid_cloud_file, ellipsoid_samples);

}