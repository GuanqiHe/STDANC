#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <thread>
#include <cmath>

#include "yaml.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include "msgpack.hpp"

#include "utils.h"

class Wang2016Single
{
public:
    typedef std::vector<double> state_type;
    // initial state
    state_type x_state;
    // control gain
    double vep, rho, beta, h, r1, r2;
    Eigen::Vector2d G;
    Eigen::RowVector2d Gamma;
    // omega_star
    double omega_star;
    Eigen::Matrix2d S;
    int a;
    double t, y, dt;

    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;

    Wang2016Single(double delta_t, double vep_param, double rho_param, int a_param, double beta_param, double h_param);
    double Outputs();
    void Update();
    void Derivatives(const state_type &x_state, state_type &x_dot_state, double t);
    void Integrals();
    void setControllerInput(double controller_input);
};

void *controllerInit(int argc, char *argv[]);

double controllerCompute(void *ctrl_ptr, double u);

void controllerFinish(void *ctrl_ptr);
