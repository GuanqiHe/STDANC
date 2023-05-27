#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

#include "base_system.hpp"

class MTController : public BaseSystem
{
public:
    using BaseSystem::state_t;

    MTController(const state_t &x_init, double w_star, double g1, double g2, double dt, double bound) : BaseSystem(x_init, dt), w_star_(w_star), bound_(bound)
    {
        G << g1, g2;
        S << 0, w_star_,
            -w_star_, 0;
    }

    void dynamics(const state_t &x, state_t &dx, double t) override
    {
        Eigen::Vector2d w_hat;
        w_hat << x[0], x[1];
        // if (w_hat.norm() >= bound_)
        // {
        //     w_hat = w_hat / w_hat.norm() * bound_;
        // }

        Eigen::Vector2d dot_w_hat = S * w_hat + G * (input);

        dx[0] = dot_w_hat(0);
        dx[1] = dot_w_hat(1);
    }

    void setInput(const double &y) override { input = y; }

    double getOutput() override
    {
        return x_[0];
    }

protected:
    double input, w_star_, bound_;
    Eigen::Matrix<double, 2, 1> G;
    Eigen::Matrix2d S;
};
