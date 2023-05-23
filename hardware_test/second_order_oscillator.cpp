#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include <iomanip>
#include <msgpack.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

struct datapack_t
{
    std::vector<std::string> state_name;
    std::vector<state_type> state_data;
    MSGPACK_DEFINE_MAP(state_name, state_data);
};

void equations(const state_type &y, state_type &dy, double _x)
{
    const double pi(boost::math::constants::pi<double>()), w_star(50 * 2 * pi);
    Eigen::MatrixXd S(2, 2);
    S << 0, w_star,
        -w_star, 0;

    Eigen::Vector2d w;
    w << y[0], y[1];

    Eigen::Vector2d x;
    x << y[2], y[3];

    Eigen::Vector2d w_hat;
    w_hat << y[4], y[5];

    Eigen::Vector2d dot_w = S * w;

    Eigen::MatrixXd A(2, 2), B(2, 1), C(1, 2);
    A << 0, 1,
        -20, -9;
    B << 0, 1;
    C << -10, 10;

    double u = w(1) + w_hat(1);
    Eigen::Vector2d dot_x = A * x + B * u;
    double out = (C * x).value();

    Eigen::MatrixXd G(2,1);
    G<< 1, 0;

    Eigen::Vector2d dot_w_hat = S * w_hat + G*(-10*1*out);

    dy[0] = dot_w(0);
    dy[1] = dot_w(1);
    dy[2] = dot_x(0);
    dy[3] = dot_x(1);
    dy[4] = dot_w_hat(0);
    dy[5] = dot_w_hat(1);
}

int main(int argc, char **argv)
{
    std::ofstream stream("output.bin", std::ios::binary);
    msgpack::packer<std::ofstream> packer(stream);

    datapack_t data;
    data.state_name = {"t", "w1", "w2", "x1", "x2"};

    const double dt = 0.0001, T = 100.0;
    runge_kutta_dopri5<state_type> stepper;

    data.state_data.reserve(1.0 / dt * T * sizeof(state_type(7)));

    state_type y(6);
    // initial values
    y[0] = 0.0; //  Theta 1
    y[1] = 2.0; // dTheta 1
    y[2] = 0.0;
    y[3] = 0.0;
    y[4] = 0.0;
    y[5] = 0.0;

    for (double t(0.0); t <= T; t += dt)
    {
        state_type state = {t, y[0], y[1], y[2], y[3], y[4], y[5]};

        data.state_data.push_back(state);

        stepper.do_step(equations, y, t, dt);
    }

    packer.pack(data);

    stream.close();
}