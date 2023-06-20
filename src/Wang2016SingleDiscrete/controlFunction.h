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

class Wang2016SingleDiscrete
{
public:
    // initial state
    Eigen::VectorXd x;
    // control gain
    double vep, rho, beta, h, r1, r2;
    Eigen::Vector2d G;
    Eigen::RowVector2d Gamma;
    // omega_star
    double omega_star;
    Eigen::Matrix2d S;
    int a, k;
    double y;

    Wang2016SingleDiscrete(double vep_param, double rho_param, int a_param, double beta_param, double h_param);
    double Outputs();
    void Update();
    void setControllerInput(double controller_input);
};

void *controllerInit(int argc, char *argv[]);

double controllerCompute(void *ctrl_ptr, double u);

void controllerFinish(void *ctrl_ptr);
