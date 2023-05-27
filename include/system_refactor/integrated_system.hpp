#include <iostream>
#include <vector>
#include <type_traits>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

#include "base_system.hpp"

class IntegratedSystem
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using state_t = Eigen::VectorXd;
    using input_t = Eigen::VectorXd;
    using output_t = Eigen::VectorXd;

    IntegratedSystem(const std::shared_ptr<BaseSystem> &sys, const double &dt)
        : sys_(sys), dt_(dt), t_(0.0){};
    IntegratedSystem(const IntegratedSystem &other) : dt_(other.dt_), t_(other.t_){};
    ~IntegratedSystem() = default;

    virtual IntegratedSystem *clone() const { new IntegratedSystem(*this); };

    void update(void)
    {
        auto func = std::bind(&BaseSystem::dynamics, &(*sys_), std::placeholders::_2, std::placeholders::_2, std::placeholders::_3);
        stepper.do_step(func, *sys_, t_, dt_);
    }
    static void interface(IntegratedSystem &controller, const input_t &u, output_t &y)
    {
        controller.sys_->setInput(u);
        controller.update();
        controller.sys_->getOutput(y);
    }

protected:
    boost::numeric::odeint::runge_kutta4<state_t> stepper;
    std::shared_ptr<BaseSystem> sys_;
    double dt_, t_;
};
