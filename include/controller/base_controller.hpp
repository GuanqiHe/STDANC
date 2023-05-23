#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

class BaseController
{
public:
    typedef std::vector<double> state_t;

    BaseController(const state_t &initial_state, const double &dt) : x_(initial_state), dt_(dt), t_(0.0){};
    ~BaseController(){};

    virtual void dynamics(const state_t &x, state_t &dx, double t) = 0;
    virtual void setInput(const double &y) = 0;
    virtual double getOutput() = 0;
    void update(void)
    {
        auto func = std::bind(&BaseController::dynamics, &(*this), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        stepper.do_step(func, x_, t_, dt_);
    }
    static void interface(BaseController &controller, const double &y, double &u)
    {
        controller.setInput(y);
        controller.update();
        u = controller.getOutput();
    }

private:
    boost::numeric::odeint::runge_kutta_dopri5<state_t> stepper;
    state_t x_;
    double dt_, t_;
};
