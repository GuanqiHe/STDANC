#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

class BaseSystem
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef Eigen::VectorXd state_t;
    typedef Eigen::VectorXd input_t;
    typedef Eigen::VectorXd output_t;
    typedef double time_t;

    BaseSystem(unsigned int x_dim, unsigned int u_dim, unsigned int y_dim, const state_t &x0)
        : x_(x_dim), u_(u_dim), y_(y_dim)
    {
        x_ = x0;
        u_.setZero();
        y_.setZero();
    };
    BaseSystem(const BaseSystem &other) : x_(other.x_), u_(other.u_), y_(other.y_){};
    ~BaseSystem() = default;

    virtual BaseSystem *clone() const = 0;

    void dynamics_int(const state_t &x, state_t &dx, time_t t){
        dynamics(x_, u_, t, dx)};

    void setInput(const input_t &u) { u_ = u; }

    virtual void dynamics(const state_t &x, const input_t &u, const double &t, state_t &dx) const = 0;
    virtual void getOutput(const state_t &x, output_t &y) const = 0;

protected:
    state_t x_;
    input_t u_;
    output_t y_;
};
