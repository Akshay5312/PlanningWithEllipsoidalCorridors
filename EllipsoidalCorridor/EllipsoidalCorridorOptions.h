#pragma once

#include <Eigen/Dense>

template<typename ReferencePathType>
// This struct is used to define the options for calculating the Ellipsoidal Corridor
struct EllipsoidalCorridorOptions{
    Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor
    ReferencePathType reference_path; // The initial guess for the path

    int N_q; // The number of dimensions for the point cloud.

    int num_ctr_pts = 5; // The number of control points for the B-spline trajectory
    int num_samples = 100; // The number of samples for the trajectory

    double r_diagonal = 10.0; // The diagonal dominance for the ellipsoid matrix
    double r_f = 1.0; // The diagonal dominance for the ellipsoid matrix


    // The following are NOT IMPLEMENTED, but may be nice to have.
    // Enforce that the corridor contains the reference path. 
    // This is done by constraining f_eps <= 1 \forall epsilon in [0, 1].
    bool contain_reference_path = false;
    // Enforce that the corridor "tapers" at the boundaries.
    bool taper = false; 
    // Constrain the maximum diffeence between elements defining two consequtive ellipsoids.
    double ellipsoidal_velocity_constraint = 1e6;

    double time_horizon = 1;

    // C^-1*E^-1*d approximates the offset of the center of the ellipsoid from the reference trajectory.
    // Assuming that E is significantly diagonally dominant, and the ellipsoid is approximately a ball (i.e. E \approx I),
    // then "d" dictates the direction in which the center of the ellipsoid is offset from the reference trajectory.
    // This can be exploited to create corridors with differing homotopy by adding some random weight to "d" in the cost. 
    // The central trajectory will be in (approximately) the opposing direction of the weight.
    bool randomize_homotopy = false;
};