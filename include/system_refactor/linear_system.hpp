#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

#include "base_system.hpp"

class LinearSystem : virtual public BaseSystem
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Base = BaseSystem;
    using Base::input_t;
    using Base::output_t;
    using Base::state_t;
    using Base::time_t;
    using matrix_t = Eigen::MatrixXd;

    LinearSystem(unsigned int x_dim, unsigned int u_dim, unsigned int y_dim,
                 const matrix_t &A, const matrix_t &B, const matrix_t &C, const state_t &x0)
        : Base(x_dim, u_dim, y_dim, x0), A_(x_dim, x_dim), B_(x_dim, u_dim), C_(y_dim, x_dim)
    {
        A_ = A;
        B_ = B;
        C_ = C;
    }

    virtual LinearSystem *clone() const override { return new LinearSystem(*this); };

    void dynamics(const state_t &x, state_t &dx, double t) override
    {
        dx = A * x + B * u_;
    }

protected:
    matrix_t A_, B_, C_;
};
