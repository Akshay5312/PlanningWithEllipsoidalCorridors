#pragma once
#include <Eigen/Dense>

#include <drake/common/trajectories/piecewise_polynomial.h>
#include <drake/solvers/solve.h>
#include <drake/solvers/mathematical_program.h>

#include <drake/common/text_logging.h>


namespace CorrGen{


template<typename T>
class Ball{
    public:
    Ball(const Eigen::VectorX<T>& center_offset, const T& f, const Eigen::VectorX<T>& ref)
        : d_(center_offset), f_(f), ref_(ref), Nq_(center_offset.rows()) {
            assert(r.rows() == d_.rows());
    }

    Eigen::VectorX<T> center(){
        return d_ + ref_;
    }

    T evaluate(const Eigen::VectorXd& point){
        Eigen::VectorX<T> point_offset = point - ref_;
        return (( point_offset.transpose() * point_offset) - 2*(point_offset.transpose() * d_))[0] + f_;
    }

    Eigen::VectorX<T>& d(){
        return d_;
    }

    T& f(){
        return f_;
    }

    Eigen::VectorX<T>& ref(){
        return ref_;
    }

    void EnforceNonMembership(drake::solvers::MathematicalProgram& prog, const Eigen::VectorXd& point){
        prog.AddLinearConstraint(evaluate(point) >= 1);
    }
    void EnforceNonMembership(drake::solvers::MathematicalProgram& prog, const std::vector<Eigen::VectorXd>& points){
        for(const auto& point : points ){
            EnforceNonMembership(prog, point);
        }
    }
    void EnforceNonMembership(drake::solvers::MathematicalProgram& prog, const Eigen::MatrixXd& points){
        for(int i = 0; i < points.cols(); i++){
            EnforceNonMembership(prog, Eigen::VectorXd(points.col(i)));
        }
    }

    // NLP methods

    T squared_radius(){
        return 1.0 + d_.transpose() * d_ - f_;
    }

    void EnforceMinRadius(drake::solvers::MathematicalProgram& prog, double r_min){
        prog.AddConstraint(squared_radius() >= r_min*r_min);
    }

    void AddMaxRadiusCost(drake::solvers::MathematicalProgram& prog){
        prog.AddCost(-squared_radius());
    }

    void AddMinRadiusCost(drake::solvers::MathematicalProgram& prog){
        prog.AddCost(squared_radius());
    }

    void EnforceMaxOffset(drake::solvers::MathematicalProgram& prog, double offset_max){
        prog.AddQuadraticConstraint(d_.squaredNorm(), 0, offset_max*offset_max);
    }


    // APPROXIMATE THE NLP WITH AN LP
    // We approximate ||d||^2 with some <d_hat,d> * r_guess
    // r_guess is the guess for ||d||. For now let us set it to some offset_max.
    // 

    T approximate_squared_radius(const Eigen::VectorXd& d_hat, double r_guess = 1.0){
        return 1.0 + r_guess*d_hat.transpose() * d_ - f_;
    }

    void EnforceConicDirectionality(drake::solvers::MathematicalProgram& prog, const Eigen::VectorXd& d_hat, double r=0){
        // (d_hat.transpose()*d) * d_hat = d & ||d_hat|| = 1? d is along d_hat
        // We enforce that 'd_' is along some _linear_ cone with axis d_hat. 
        auto d_hat_ = d_hat.normalized();
        r = abs(r);
        prog.AddLinearConstraint( (d_hat_*d_hat_.transpose()) * d_ - d_ <= Eigen::VectorXd::Ones(Nq_)*r );
        prog.AddLinearConstraint( (d_hat_*d_hat_.transpose()) * d_ - d_ >= -Eigen::VectorXd::Ones(Nq_)*r );
    }


    void ApproxEnforceMinRadius(drake::solvers::MathematicalProgram& prog, double r_min, const Eigen::VectorXd& d_hat, double r_guess = 1.0){
        prog.AddLinearConstraint(approximate_squared_radius(d_hat, r_guess) >= r_min*r_min);
    }

    void ApproxAddMaxRadiusCost(drake::solvers::MathematicalProgram& prog, const Eigen::VectorXd& d_hat, double r_guess = 1.0){
        prog.AddLinearCost(-approximate_squared_radius(d_hat, r_guess));
    }

    void ApproxAddMinRadiusCost(drake::solvers::MathematicalProgram& prog, const Eigen::VectorXd& d_hat, double r_guess = 1.0){
        prog.AddLinearCost(approximate_squared_radius(d_hat, r_guess));
    }

    void ApproxEnforceMaxOffset(drake::solvers::MathematicalProgram& prog, const Eigen::VectorXd& d_hat, double offset_max){
        prog.AddLinearConstraint(d_hat.normalized().transpose()*d_ <= offset_max);
    }


    std::vector<Eigen::VectorX<T>> sample(int N_samples) {
        // Sample points within the ellipsoid
        Eigen::VectorX<T> q_tilde = center();

        T radius_sq = squared_radius();
        T radius = sqrt(radius_sq);

        std::vector<Eigen::VectorX<T>> samples(N_samples);
        for (int i = 0; i < N_samples; ++i) {
            Eigen::VectorXd random_point = Eigen::VectorXd::Random(Nq_);
            random_point.normalize();
            samples[i] = random_point * radius + q_tilde;
        }
        return samples;
    }


    private:
    Eigen::VectorX<T> d_;
    Eigen::VectorX<T> ref_;
    T f_;
    int Nq_;
};


}