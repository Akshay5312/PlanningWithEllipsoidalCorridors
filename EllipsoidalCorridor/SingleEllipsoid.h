#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>

#include "EllipsoidalCorridorOptions.h"
#include "LagrangePolynomial.h"


struct InflatedEllipsoidOptions{
    Eigen::MatrixXd point_cloud; // The point cloud used to define the corridor

    Eigen::VectorXd reference_point; // The reference point for the ellipsoid
    int N_q; // The number of dimensions for the point cloud.

    double r_diagonal = 10.0; // The diagonal dominance for the ellipsoid matrix
    double r_f = 1.0; // The diagonal dominance for the ellipsoid matrix
};


template<typename T>
struct InflatedEllipsoid{
    Eigen::MatrixX<T> E_; // Ellipsoid matrix
    Eigen::VectorX<T> d_; // Ellipsoid vector
    T f_; // Ellipsoid scalar
    Eigen::VectorX<double> reference_point_; // Reference point for the ellipsoid

    InflatedEllipsoid(Eigen::MatrixX<T> E, Eigen::VectorX<T> d, T f, Eigen::VectorX<double> reference_point)
        : E_(E), d_(d), f_(f), reference_point_(reference_point) {
            drake::log()->info("Ellipsoid created");
        }

    template<typename Tpoint>
    T evaluate(Eigen::VectorX<Tpoint> point) const {
        // Evaluates if 'point' is in the ellipsoidal cross section of the corridor.
        // Corridor(e) = {q | (q-c_e)E_e(q-c_e) - 2*d_e(q-c_e) + f_e <= 1}

        auto o_point = point - reference_point_; // Assuming c is the center of the ellipsoid
        return (o_point.transpose() * E_ * o_point - 2 * d_.transpose() * o_point) + f_ - 1;
    }

    Eigen::VectorX<T> center() const {
        // Calculate the center of the ellipsoid
        Eigen::MatrixX<T> E_inv = E_.inverse();
        return reference_point_ + (E_inv * d_);
    }

};

InflatedEllipsoid<double> GenerateInflatedEllipsoid(
    const InflatedEllipsoidOptions& options) {

    drake::log()->info("Generating Ellipsoid");

    int N_q = options.N_q;
    double r_diagonal = options.r_diagonal;
    double r_f = options.r_f;
    auto point_cloud = options.point_cloud; // Each col is a point
    auto ref_point = options.reference_point; // Reference point for the ellipsoid

    drake::log()->info("Ellipsoid generated");

    drake::solvers::MathematicalProgram prog;

    Eigen::MatrixX<drake::symbolic::Variable> E_var = prog.NewContinuousVariables(N_q, N_q, "E");
    Eigen::VectorX<drake::symbolic::Variable> d_var = prog.NewContinuousVariables(N_q, 1, "d");
    Eigen::VectorX<drake::symbolic::Variable> f_var = prog.NewContinuousVariables(1, 1, "f");

    Eigen::MatrixX<drake::symbolic::Expression> E = E_var;
    Eigen::VectorX<drake::symbolic::Expression> d = d_var;
    Eigen::VectorX<drake::symbolic::Expression> f = f_var;

    InflatedEllipsoid<drake::symbolic::Expression> cal_C(E, d, f[0], ref_point);

    for (int i = 0; i < N_q; ++i) {
        for (int j = 0; j < N_q; ++j) {
            if (i != j) {

                // Add constraints for diagonal dominance
                prog.AddLinearConstraint(E(i, i) <= r_diagonal * E(i, j));
                prog.AddLinearConstraint(E(i, i) >= -r_diagonal * E(i, j));
            
                // Add constraints for symmetry
                prog.AddLinearConstraint(E(i, j) == E(j, i)); // Symmetry constraint
            }
        }
        prog.AddLinearConstraint(E(i, i) >= 0.1);
    }
    
    // The reference must be in the interior of the ellipsoid, so we add a constraint for that
    // Eigen::VectorXd zero_point = Eigen::VectorXd::Zero(N_q);
    prog.AddLinearConstraint(cal_C.evaluate(ref_point) <= 0);
    // Add non-membership constraints for each point in the point cloud
    for (int p_index = 0; p_index < options.point_cloud.cols(); ++p_index) {
        Eigen::VectorXd point = options.point_cloud.col(p_index);

        // Non-membership constraint
        prog.AddLinearConstraint(cal_C.evaluate(point) >= 0);
    }

    // Add a cost to minimize the trace of the ellipsoid matrix E and the scalar f
    prog.AddLinearCost(E.trace() + N_q * r_f * f[0]);

    prog.AddLinearConstraint(f[0] == 0); // Ensure f is bounded


    auto result = drake::solvers::Solve(prog);

    if (!result.is_success()) {
        throw std::runtime_error("Ellipsoid generation failed");
    }

    // Extract the optimized values
    Eigen::MatrixX<double> E_opt = result.GetSolution(E_var);
    Eigen::VectorX<double> d_opt = result.GetSolution(d_var);
    Eigen::VectorX<double> f_opt = result.GetSolution(f_var);

    // Create the ellipsoid object
    InflatedEllipsoid<double> ellipsoid(E_opt, d_opt, f_opt[0], ref_point);
    drake::log()->info("Ellipsoid generated successfully");

    return ellipsoid;
}