#include "controlFunction.h"
#include "logger/data_logging.hpp"

data_logger logger;

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

    double w_star = config["dist_freq"].as<double>() * M_PI * 2;
    double theta1 = config["theta1"].as<double>();
    double theta2 = config["theta2"].as<double>();
    double alpha = config["alpha"].as<double>();
    double epsilon = config["epsilon"].as<double>();
    double gamma = config["gamma"].as<double>();
    double dt = 1 / config["sample_fs"].as<double>();

    SwitchingBasedAFC *ctrl_ptr = new SwitchingBasedAFC(w_star, dt, theta1, theta2, alpha, epsilon, gamma);

    int sample_len = (config["run_time"].as<double>()) * config["sample_fs"].as<double>();

    logger.init({"t", "y", "w0", "w1", "zeta0", "zeta1", "xi0", "xi1", "eta0", "eta1"}, sample_len, config["controller_log_path"].as<std::string>());

    return (void *)(ctrl_ptr);
}

double controllerCompute(void *ctrl_ptr, double y)
{
    SwitchingBasedAFC *ptr = (SwitchingBasedAFC *)(ctrl_ptr);

    logger.log({ptr->tw, ptr->out, ptr->w[0], ptr->w[1], ptr->w[2], ptr->w[3], ptr->w[4], ptr->w[5], ptr->w[6], ptr->w[7]});

    ptr->setInput(y);
    return ptr->computeOutput();
}

void controllerFinish(void *ctrl_ptr)
{
    logger.write();
    SwitchingBasedAFC *ptr = (SwitchingBasedAFC *)(ctrl_ptr);
    delete ptr;
}
