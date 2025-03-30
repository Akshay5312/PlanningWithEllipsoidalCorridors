#pragma once

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <eigen3/Eigen/Dense>

namespace CorrGen{
    template<typename T>
    class LagrangePolynomial{
        public:
        LagrangePolynomial(std::vector<double> break_times, std::vector<Eigen::MatrixX<T>> control_points)
            : break_times_(break_times), ctrl_points_(Eigen::MatrixX<Eigen::VectorX<T>>(control_points[0].rows(), control_points[0].cols()))
        {
            
            assert(break_times.size() == control_points.size());
            N_breaks = break_times.size();

            assert(break_times[0] == 0);
            // assert(N_breaks == 1);

            for(int i = 1; i < N_breaks-1; i++){
                assert(break_times[i-1] < break_times[i] );
                assert(break_times[i] > 0);
                assert(break_times[i] < 1);

                drake::log()->info("break_times[{}]: {}", i, break_times[i]);
            }

            assert(break_times[N_breaks-1] == 1);

            out_rows = control_points[0].rows();
            out_cols = control_points[0].cols();

            for(int i = 0; i < out_rows; i++){
                for(int j = 0; j < out_cols; j++){
                    ctrl_points_(i,j) = Eigen::VectorX<T>( N_breaks );

                    for(int k = 0; k < N_breaks; k++){
                        ctrl_points_(i,j)[k] = control_points[k](i,j);
                    }
                }
            }
        }

        Eigen::VectorXd basis(double eps) const {
            // Lagrange interpolating polynomial basis
            // L(t) = \prod_{i=0}^{N-1} (t - t_i)/(t_j - t_i)
            // where t_i are the break times and t_j are the break times not equal to t_i
            // This is a vector of size N_breaks
            Eigen::VectorXd L(N_breaks);

            for(int i = 0; i < N_breaks; i++){
                L[i] = 1;
                for(int j = 0; j < N_breaks; j++){
                    if(i != j){
                        L[i] *= (eps - break_times_[j])/(break_times_[i] - break_times_[j]);
                    }
                }
            }

            return L;
        }

        Eigen::MatrixX<T> value(double eps) const {
            Eigen::VectorXd L = basis(eps);
            Eigen::MatrixX<T> result(out_rows, out_cols);

            for(int i = 0; i < out_rows; i++){
                for(int j = 0; j < out_cols; j++){
                    result(i,j) = L.transpose() * ctrl_points_(i,j);
                }
            }

            return result;
        } 


        private:
        int N_breaks;
        std::vector<double> break_times_;
        Eigen::MatrixX<Eigen::VectorX<T>> ctrl_points_;

        int out_rows = 1;
        int out_cols = 1;
    };
}