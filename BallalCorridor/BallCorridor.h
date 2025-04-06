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

        LagrangePolynomial<T> d_;
        LagrangePolynomial<T> f_;
        LagrangePolynomial<double> ref_;
        int Nq_;
    };
}


