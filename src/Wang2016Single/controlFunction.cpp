#include "controlFunction.h"
#include "logger/data_logging.hpp"

data_logger logger;

#define M_PI 3.14159265358979323846

template <typename XprType, typename RowFactorType, typename ColFactorType>
auto repelem(const XprType &xpr, RowFactorType row_factor, ColFactorType col_factor)
{
    using namespace Eigen;

    const int RowFactor = internal::get_fixed_value<RowFactorType>::value;
    const int ColFactor = internal::get_fixed_value<ColFactorType>::value;
    const int NRows = XprType::RowsAtCompileTime == Dynamic || RowFactor == Dynamic ? Dynamic : XprType::RowsAtCompileTime * RowFactor;
    const int NCols = XprType::ColsAtCompileTime == Dynamic || ColFactor == Dynamic ? Dynamic : XprType::ColsAtCompileTime * ColFactor;
    const int nrows = internal::get_runtime_value(row_factor) * xpr.rows();
    const int ncols = internal::get_runtime_value(col_factor) * xpr.cols();

    return xpr(
        Array<int, NRows, 1>::LinSpaced(nrows, 0, xpr.rows() - 1),
        Array<int, NCols, 1>::LinSpaced(ncols, 0, xpr.cols() - 1));
}

Wang2016Single::Wang2016Single(double delta_t, double vep_param, double rho_param, int a_param, double beta_param, double h_param)
{
    vep = vep_param;
    rho = rho_param;
    beta = beta_param;
    h = h_param;
    r1 = 0.1;
    r2 = 3;
    omega_star = 50 * 2 * M_PI;
    G << 1, 0;
    Gamma << 1, 0;
    S << 0, omega_star,
        -omega_star, 0;
    t = 0;
    y = 0;
    a = a_param;
    dt = delta_t;
    Eigen::VectorXd chi_hat_init(6);
    chi_hat_init << 1, 0, -1, 0, 0, -1;

    Eigen::VectorXd x(19);
    x.setZero();
    x.segment(6, 6) << chi_hat_init;
    x.segment(17, 2) << 2 * Eigen::Vector2d::Random();
    state_type x_temp(x.data(), x.data() + x.size());
    x_state = x_temp;
}

double Wang2016Single::Outputs()
{
    double u;
    Eigen::Vector2d v_hat;
    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x_state.data(), x_state.size());
    v_hat << x.segment(17, 2);
    u = Gamma * v_hat;
    return u;
}

void Wang2016Single::Update()
{
    Eigen::Vector3d J;
    Eigen::Index min_ind;
    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x_state.data(), x_state.size());
    J << x.segment(14, 3);
    J(a) = J(a) - h;
    J.minCoeff(&min_ind);
    a = min_ind;
    // a = 1;
}

void Wang2016Single::Derivatives(const state_type &x_state, state_type &x_dot_state, double t)
{
    Eigen::Vector2d xi1_hat, xi1_hat_dot, v_hat, v_hat_dot, chi1, chi2, chi3, vph1, vph2, vph3, chi_hat_dot1, chi_hat_dot2, chi_hat_dot3;
    Eigen::Vector3d J, e, y1_hat, y_til, m_sr, J_dot;
    Eigen::VectorXd zeta_hat(6), chi_hat(6), zeta_hat_dot(6), chi_hat_dot(6), xi1_hat_replicate(6), y_til_repelem(6), m_sr_repelem(6), vph(6), x_dot(19);
    Eigen::Matrix2d F, Ical1, Ical2, Ical3;
    Eigen::MatrixXd chi_hat_b(2, 3), zeta_hat_b(2, 3), blkdiag_Gamma(3, 6), blkdiag_S(6, 6), blkdiag_G(6, 3), vph_b(2, 3);
    double e_a;
    state_type x_temp = x_state;
    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x_temp.data(), x_temp.size());

    zeta_hat << x.segment(0, 6);
    chi_hat << x.segment(6, 6);
    xi1_hat << x.segment(12, 2);
    J << x.segment(14, 3);
    v_hat << x.segment(17, 2);

    F << S - vep * G * Gamma;

    // Multiple - model estimator and switching controller
    zeta_hat_b << zeta_hat.reshaped(2, 3);
    chi_hat_b << chi_hat.reshaped(2, 3);
    e << -vep * chi_hat_b.cwiseProduct(zeta_hat_b).colwise().sum().transpose();
    e_a = e(a);

    blkdiag_Gamma.setZero();
    blkdiag_Gamma(0, 0) = 1;
    blkdiag_Gamma(1, 2) = 1;
    blkdiag_Gamma(2, 4) = 1;
    y1_hat << blkdiag_Gamma * zeta_hat;

    blkdiag_S.setZero();
    blkdiag_S.topLeftCorner(2, 2) = S;
    blkdiag_S.block(2, 2, 2, 2) = S;
    blkdiag_S.bottomRightCorner(2, 2) = S;

    blkdiag_G << blkdiag_Gamma.transpose();
    y_til << y1_hat.array() - y;

    zeta_hat_dot << blkdiag_S * zeta_hat + chi_hat * e_a - vep * blkdiag_G * y_til;

    // Multiple - model update law
    m_sr << 1 + xi1_hat.squaredNorm() + y_til.array().square();

    xi1_hat_replicate << xi1_hat.replicate(3, 1);
    y_til_repelem << repelem(y_til, 2, 1);
    m_sr_repelem << repelem(m_sr, 2, 1);

    vph << -rho * pow(vep, 2) * xi1_hat_replicate.cwiseProduct(y_til_repelem).cwiseQuotient(m_sr_repelem);
    vph_b << vph.reshaped(2, 3);

    chi1 << chi_hat_b.col(0);
    chi2 << chi_hat_b.col(1);
    chi3 << chi_hat_b.col(2);

    vph1 << vph_b.col(0);
    vph2 << vph_b.col(1);
    vph3 << vph_b.col(2);

    if (sqrt(3) * chi1(0) + chi1(1) < r1 || chi1(0) / sqrt(3) + chi1(1) < 0 || chi1(0) < 0)
    {
        Ical1.setIdentity();
    }
    else if (chi1.norm() > r2)
    {
        Ical1 << chi1 * chi1.transpose() / (chi1.transpose() * chi1);
    }
    else
    {
        Ical1.setZero();
    }

    if (-sqrt(3) * chi2(0) + chi2(1) < r1 || -chi2(0) / sqrt(3) + chi2(1) < 0 || chi2(0) > 0)
    {
        Ical2.setIdentity();
    }
    else if (chi2.norm() > r2)
    {
        Ical2 << chi2 * chi2.transpose() / (chi2.transpose() * chi2);
    }
    else
    {
        Ical2.setZero();
    }

    if (chi3(0) / sqrt(3) + chi3(1) > 0 || -chi3(0) / sqrt(3) + chi3(1) > 0 || chi3(1) > -r1 / 2)
    {
        Ical3.setIdentity();
    }
    else if (chi3.norm() > r2)
    {
        Ical3 << chi3 * chi3.transpose() / (chi3.transpose() * chi3);
    }
    else
    {
        Ical3.setZero();
    }

    chi_hat_dot1 << vph1 - Ical1 * vph1;
    chi_hat_dot2 << vph2 - Ical2 * vph2;
    chi_hat_dot3 << vph3 - Ical3 * vph3;

    chi_hat_dot << chi_hat_dot1, chi_hat_dot2, chi_hat_dot3;
    xi1_hat_dot << F.transpose() * xi1_hat + G * e_a;

    // Hysteresis switching logic
    J_dot << beta * y_til.array().square();

    // Certainty - equivalence estimator
    v_hat_dot << S * v_hat + G * e_a;

    x_dot << zeta_hat_dot, chi_hat_dot, xi1_hat_dot, J_dot, v_hat_dot;
    state_type dx_temp(x_dot.data(), x_dot.data() + x_dot.size());
    x_dot_state = dx_temp;
}

void Wang2016Single::Integrals()
{
    auto func = std::bind(&Wang2016Single::Derivatives, &(*this), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    stepper.do_step(func, x_state, t, dt);
    t += dt;
}

void Wang2016Single::setControllerInput(double controller_input)
{
    y = controller_input;
}

void *controllerInit(int argc, char *argv[])
{
    std::string config_path = "config.yaml";
    if (argc >= 3)
    {
        config_path = std::string(argv[2]);
    }
    else if (argc == 2)
    {
        config_path = std::string(argv[1]);
    }

    YAML::Node config = YAML::LoadFile(config_path);

    const std::string log_folder = config["log_folder"].as<std::string>();
    const std::string controller_log_path = log_folder + '/' + config["controller_log_path"].as<std::string>();
    // const double w_star = config["dist_freq"].as<double>() * M_PI * 2;
    const double dt = 1 / config["sample_fs"].as<double>();
    const double vep_param = config["vep_param"].as<double>();
    const double rho_param = config["rho_param"].as<double>();
    const int a_param = config["a_param"].as<int>();
    const double beta_param = config["beta_param"].as<double>();
    const double h_param = config["h_param"].as<double>();

    Wang2016Single *ctrl_ptr = new Wang2016Single(dt, vep_param, rho_param, a_param, beta_param, h_param);

    int sample_len = (config["run_time"].as<double>()) * config["sample_fs"].as<double>();

    logger.init({"t", "y", "a", "chi_hat_1", "chi_hat_2", "J_1", "J_2", "J_3"}, sample_len, controller_log_path);

    std::cout << "Controller log path: " << controller_log_path << std::endl;

    return (void *)(ctrl_ptr);
}

double controllerCompute(void *ctrl_ptr, double controller_input)
{
    Wang2016Single *ptr = (Wang2016Single *)(ctrl_ptr);
    double u;

    u = ptr->Outputs();
    ptr->setControllerInput(controller_input);
    ptr->Integrals();
    ptr->Update();

    logger.log({ptr->t, ptr->y, double(ptr->a), ptr->x_state[6 + 2 * ptr->a], ptr->x_state[7 + 2 * ptr->a], ptr->x_state[14], ptr->x_state[15], ptr->x_state[16]});

    return u;
}

void controllerFinish(void *ctrl_ptr)
{
    logger.write();
    Wang2016Single *ptr = (Wang2016Single *)(ctrl_ptr);
    delete ptr;
}