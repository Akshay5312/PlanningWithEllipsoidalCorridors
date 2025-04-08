#pragma once
#include <Eigen/Dense>
#include "Ball.h"
#include "LagrangePolynomial.h"

namespace CorrGen{
    template<typename T>
    class BallCorridor{
        public:
        BallCorridor(const LagrangePolynomial<T>& d, const LagrangePolynomial<T>& f, const LagrangePolynomial<double>& ref)
            : d_(d), f_(f), ref_(ref) {
            Nq_ = d_.value(0).rows();
        }

        static BallCorridor FromString(const std::string& str){
            // Parse the string to initialize the object
            std::istringstream iss(str);
            std::string line;
            std::vector<double> eps_values;
            std::vector<Eigen::VectorX<T>> d_values;
            std::vector<T> f_values;
            std::vector<Eigen::VectorX<double>> ref_values;

            while (std::getline(iss, line)) {
                if (line.find("eps[") != std::string::npos) {
                    size_t pos = line.find("=");
                    eps_values.push_back(std::stod(line.substr(pos + 1)));
                } else if (line.find("d[") != std::string::npos) {
                    size_t pos = line.find("=");
                    std::istringstream vec_stream(line.substr(pos + 1));
                    Eigen::VectorX<T> vec(Nq_);
                    for (int i = 0; i < Nq_; ++i) {
                        vec_stream >> vec(i);
                    }
                    d_values.push_back(vec);
                } else if (line.find("f[") != std::string::npos) {
                    size_t pos = line.find("=");
                    f_values.push_back(static_cast<T>(std::stod(line.substr(pos + 1))));
                } else if (line.find("ref[") != std::string::npos) {
                    size_t pos = line.find("=");
                    std::istringstream vec_stream(line.substr(pos + 1));
                    Eigen::VectorX<double> vec(Nq_);
                    for (int i = 0; i < Nq_; ++i) {
                        vec_stream >> vec(i);
                    }
                    ref_values.push_back(vec);
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
            for(int i = 0; i < d_.break_times().size(); i++){
                str += "eps[" + std::to_string(i) + "] = " + std::to_string(d_.break_times()[i]) + "\n";
            }
            for(int i = 0; i < d_.break_times().size(); i++){
                str += "d[" + std::to_string(i) + "] = " + d_.value(d_.break_times()[i]).transpose().format(Eigen::IOFormat(Eigen::FullPrecision)) + "\n";
            }
            for(int i = 0; i < f_.break_times().size(); i++){
                str += "f[" + std::to_string(i) + "] = " + std::to_string(f_.value(f_.break_times()[i])(0,0)) + "\n";
            }
            for(int i = 0; i < ref_.break_times().size(); i++){
                str += "ref[" + std::to_string(i) + "] = " + ref_.value(ref_.break_times()[i]).transpose().format(Eigen::IOFormat(Eigen::FullPrecision)) + "\n";
            }
            return str;
        }



        LagrangePolynomial<T> d_;
        LagrangePolynomial<T> f_;
        LagrangePolynomial<double> ref_;
        int Nq_;
    };
}


