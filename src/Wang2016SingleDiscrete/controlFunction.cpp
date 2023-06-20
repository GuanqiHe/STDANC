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

Wang2016SingleDiscrete::Wang2016SingleDiscrete(double vep_param, double rho_param, int a_param, double beta_param, double h_param)
{
    vep = vep_param;
    rho = rho_param;
    beta = beta_param;
    h = h_param;
    r1 = 0.1;
    r2 = 3;
    omega_star = 50 * 2 * M_PI / 5000;
    G << 1, 0;
    Gamma << 1, 0;
    S << std::cos(omega_star), std::sin(omega_star),
        -std::sin(omega_star), std::cos(omega_star);
    k = 0;
    y = 0;
    a = a_param;
    Eigen::VectorXd chi_hat_init(6);
    chi_hat_init << 1, 0, -1, 0, 0, -1;

    x.setZero(19);
    x.segment(6, 6) << chi_hat_init;
    x.segment(17, 2) << 2 * Eigen::Vector2d::Random();
}

double Wang2016SingleDiscrete::Outputs()
{
    double u;
    Eigen::Vector2d v_hat;
    v_hat << x.segment(17, 2);
    u = Gamma * v_hat;
    return u;
}

void Wang2016SingleDiscrete::Update()
{
    Eigen::Vector2d xi1_hat, v_hat, chi1, chi2, chi3, vph1, vph2, vph3;
    Eigen::Vector3d J, e, y1_hat, y_til, m_sr, y_til_sq;
    Eigen::VectorXd zeta_hat(6), chi_hat(6), xi1_hat_replicate(6), y_til_repelem(6), m_sr_repelem(6), vph(6);
    Eigen::Matrix2d F, Ical1, Ical2, Ical3;
    Eigen::MatrixXd chi_hat_b(2, 3), zeta_hat_b(2, 3), blkdiag_Gamma(3, 6), blkdiag_S(6, 6), blkdiag_G(6, 3), vph_b(2, 3);
    double e_a;
    Eigen::Index min_ind;

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

    zeta_hat << blkdiag_S * zeta_hat + chi_hat * e_a - vep * blkdiag_G * y_til;

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

    chi1 << chi1 + vph1 - Ical1 * vph1;
    chi2 << chi2 + vph2 - Ical2 * vph2;
    chi3 << chi3 + vph3 - Ical3 * vph3;

    chi_hat << chi1, chi2, chi3;
    xi1_hat << F.transpose() * xi1_hat + G * e_a;

    // Hysteresis switching logic
    y_til_sq << y_til.array().square();
    J << J + beta * y_til_sq;

    // Certainty - equivalence estimator
    v_hat << S * v_hat + G * e_a;

    x << zeta_hat, chi_hat, xi1_hat, J, v_hat;

    // Update switching controller
    J(a) = J(a) - h;
    J.minCoeff(&min_ind);
    a = min_ind;

    k += 1;
}

void Wang2016SingleDiscrete::setControllerInput(double controller_input)
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
    const double vep_param = config["vep_param"].as<double>();
    const double rho_param = config["rho_param"].as<double>();
    const int a_param = config["a_param"].as<int>();
    const double beta_param = config["beta_param"].as<double>();
    const double h_param = config["h_param"].as<double>();

    Wang2016SingleDiscrete *ctrl_ptr = new Wang2016SingleDiscrete(vep_param, rho_param, a_param, beta_param, h_param);

    int sample_len = (config["run_time"].as<double>()) * config["sample_fs"].as<double>();

    logger.init({"k", "y", "a", "chi_hat_1", "chi_hat_2", "J_1", "J_2", "J_3"}, sample_len, controller_log_path);

    std::cout << "Controller log path: " << controller_log_path << std::endl;

    return (void *)(ctrl_ptr);
}

double controllerCompute(void *ctrl_ptr, double controller_input)
{
    Wang2016SingleDiscrete *ptr = (Wang2016SingleDiscrete *)(ctrl_ptr);
    double u;

    u = ptr->Outputs();
    ptr->setControllerInput(controller_input);
    ptr->Update();

    logger.log({double(ptr->k), ptr->y, double(ptr->a), ptr->x(6 + 2 * ptr->a), ptr->x(7 + 2 * ptr->a), ptr->x(14), ptr->x(15), ptr->x(16)});

    return u;
}

void controllerFinish(void *ctrl_ptr)
{
    logger.write();
    Wang2016SingleDiscrete *ptr = (Wang2016SingleDiscrete *)(ctrl_ptr);
    delete ptr;
}