#pragma once
#include <Eigen/Dense>
#include "Ball.h"
#include "../PathParameterization/LagrangePolynomial.h"

namespace CorrGen{
    template<typename T>
    class BallCorridor{
        public:
        BallCorridor(const LagrangePolynomial<T>& d, const LagrangePolynomial<T>& f, const LagrangePolynomial<double>& ref)
            : d_(d), f_(f), ref_(ref) {
            Nq_ = d_.value(0).rows();
        }

        BallCorridor(){}

        static BallCorridor FromString(const std::string& str){
            // Parse the string to initialize the object
            std::istringstream iss(str);
            std::string line;
            std::vector<double> eps_values;
            std::vector<Eigen::MatrixX<T>> d_values;
            std::vector<Eigen::MatrixX<T>> f_values;
            std::vector<Eigen::MatrixX<double>> ref_values;
            int Nq = 0;
                    
            while (std::getline(iss, line)) {
                drake::log()->info("line: {}", line);
                if (line.find("eps[") != std::string::npos) {
                    size_t pos = line.find("=");
                    eps_values.push_back(std::stod(line.substr(pos + 1)));
                } else if (line.find("d[") != std::string::npos) {
                    size_t pos = line.find("=");
                    size_t end_pos = line.find("]", pos);
                    std::istringstream vec_stream(line.substr(pos + 1));
                    
                    // Find N_q
                    if(Nq == 0){
                        std::istringstream temp_stream(line.substr(pos + 1));
                        T temp_value;
                        while (temp_stream >> temp_value) {
                            ++Nq;
                        }
                        temp_stream.clear();
                        temp_stream.seekg(0, std::ios::beg);
                    }
                    Eigen::VectorX<T> vec(Nq);
                    
                    for (int i = 0; i < Nq; ++i) {
                        vec_stream >> vec(i);
                    }
                    d_values.push_back(vec);
                    drake::log()->info("d is {}", vec.transpose());
                } else if (line.find("ref[") != std::string::npos) {
                    size_t pos = line.find("=");
                    std::istringstream vec_stream(line.substr(pos + 1));
                    Eigen::VectorX<double> vec(Nq);
                    for (int i = 0; i < Nq; ++i) {
                        vec_stream >> vec(i);
                    }
                    drake::log()->info("ref[]: {}", vec.transpose());
                    ref_values.push_back(vec);
                } else if (line.find("f[") != std::string::npos) {
                    size_t pos = line.find("=");
                    Eigen::MatrixX<T> f_val = Eigen::MatrixX<T>::Zero(1, 1);
                    f_val(0, 0) = std::stod(line.substr(pos + 1));
                    f_values.push_back(f_val);
                    drake::log()->info("f is {}", f_val);
                } 
            }

            return BallCorridor(
                LagrangePolynomial<T>(eps_values, d_values),
                LagrangePolynomial<T>(eps_values, f_values),
                LagrangePolynomial<double>(eps_values, ref_values)
            );
        }

        Eigen::VectorX<T> center(double eps) const {
            return d_.value(eps) + ref_.value(eps);
        }
        
        Eigen::VectorX<T> d_center(double eps) const {
            return d_.d_value(eps);
        }

        Eigen::VectorX<T> dd_center(double eps) const {
            return d_.dd_value(eps);
        }

        Eigen::VectorX<T> d(double eps) const {
            return d_.value(eps);
        }

        T f(double eps) const {
            return f_.value(eps)(0,0);
        }

        Ball<T> getBall(double eps) const {
            return Ball<T>(d_.value(eps), f_.value(eps)(0,0), ref_.value(eps));
        }

        std::string toString() const {
            std::string str = "BallCorridor:\n";
            auto break_times = d_.break_times();
            for(int i = 0; i < d_.break_times().size(); i++){
                str += "eps[" + std::to_string(i) + "] = " + std::to_string(d_.break_times()[i]) + "\n";
            }
            for(int i = 0; i < d_.break_times().size(); i++){
                auto d_val = d_.value(d_.break_times()[i]);
                std::stringstream ss;
                ss << d_val.transpose();
                std::string d_as_str = ss.str();
                str += "d[" + std::to_string(i) + "] = " +  d_as_str + "\n";
            }
            for(int i = 0; i < d_.break_times().size(); i++){
                str += "f[" + std::to_string(i) + "] = " + std::to_string(f_.value(d_.break_times()[i])(0,0)) + "\n";
            }
            for(int i = 0; i < d_.break_times().size(); i++){
                auto ref_val = ref_.value(d_.break_times()[i]);
                std::stringstream ss;
                ss << ref_val.transpose();
                std::string ref_as_str = ss.str();
                str += "ref[" + std::to_string(i) + "] = " +  ref_as_str + "\n";
            }
            return str;
        }

        std::string toStringOfSamples(int N_time_samples = 50){
            std::string str = "";
            std::vector<double> sample_times(N_time_samples);
            for (int k = 0; k < N_time_samples; ++k) {
                sample_times[k] = static_cast<double>(k) / (N_time_samples-1);
                auto ball_eps = getBall(sample_times[k]);
                auto ball_samples = ball_eps.sample(1000);
                for(int i = 0; i < ball_samples.size(); ++i) {
                    for(int j = 0; j < ball_samples[i].rows(); ++j) {
                        str += std::to_string(ball_samples[i][j]) + " ";
                    }
                    str += "\n";
                }
            }

            return str;
        }



        LagrangePolynomial<T> d_;
        LagrangePolynomial<T> f_;
        LagrangePolynomial<double> ref_;
        int Nq_;
    };
}


