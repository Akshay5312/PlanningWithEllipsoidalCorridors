#pragma once
#include "../PathParameterization/LagrangePolynomial.h"


namespace CorrPlanning{
class VaryingProjectionToBall{

    
    public:
    VaryingProjectionToBall(
        int N_y,
        int N_q,
        int control_grain,
        int enforce_orthonormality_grain
    ){
        N_y_ = N_y;
        N_q_ = N_q;

        // Create a RANDOM projection to a ball in N_q.
        
        std::vector<Eigen::MatrixXd> onorm_ctrl_points(control_grain);
        
        
        std::vector<double> onorm_breaks(control_grain);
        for(int i = 0; i < control_grain; i++){
            Eigen::MatrixXd random_matrix = Eigen::MatrixXd::Random(N_y, N_y);
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(random_matrix);
            Eigen::MatrixXd orthonormal_matrix = (qr.householderQ());
            orthonormal_matrix = orthonormal_matrix.leftCols(N_y);
            onorm_ctrl_points[i] = orthonormal_matrix;
            onorm_breaks[i] = static_cast<double>(i) / (control_grain-1);
        }

        // Create the Lagrange polynomial for the orthonormality constraint
        drake::trajectories::PiecewisePolynomial<double> onorm_path = 
        drake::trajectories::PiecewisePolynomial<double>::CubicWithContinuousSecondDerivatives(
            onorm_breaks,
            onorm_ctrl_points
        );

        std::vector<double> onorm_new_breaks(enforce_orthonormality_grain);

        std::vector<Eigen::MatrixXd> onorm_new_ctrl_points(enforce_orthonormality_grain);
        for(int i = 0; i < enforce_orthonormality_grain; i++){
            double eps = static_cast<double>(i) / (enforce_orthonormality_grain-1);
            Eigen::MatrixXd onorm_matrix = onorm_path.value(eps);
            // project to orthonormality using QR 
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(onorm_matrix);
            Eigen::MatrixXd orthonormal_matrix = qr.householderQ();
            orthonormal_matrix = orthonormal_matrix.leftCols(N_y);
            onorm_new_ctrl_points[i] = orthonormal_matrix;
            onorm_new_breaks[i] = eps;
        }

        // Create the Lagrange polynomial
        drake::trajectories::PiecewisePolynomial<double> onorm_new_path = 
        drake::trajectories::PiecewisePolynomial<double>::CubicWithContinuousSecondDerivatives(
            onorm_new_breaks,
            onorm_new_ctrl_points
        );

        onorm_path_ = onorm_new_path;
    }

    int N_y() const {return N_y_;}
    int N_q() const {return N_q_;}
    

    // drake::trajectories::PiecewisePolynomial<double><double>& path(){
    //     return onorm_path_;
    // }

    Eigen::MatrixXd L(double eps) const {
        // Evaluate the path at a given epsilon
        return onorm_path_.value(eps);
    }
    Eigen::MatrixXd dL(double eps) const {
        // Evaluate the derivative of the path at a given epsilon
        return onorm_path_.EvalDerivative(eps, 1);
    }
    Eigen::MatrixXd ddL(double eps) const {
        // Evaluate the second derivative of the path at a given epsilon
        return onorm_path_.EvalDerivative(eps, 2);
    }

    private:
    int N_y_;
    int N_q_;
    drake::trajectories::PiecewisePolynomial<double> onorm_path_;
};
} // namespace CorrPlanning