#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>

#include "NIDAQmx.h"
#include "yaml.h"
#include "msgpack.hpp"
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
//#include <boost/numeric/odeint.hpp>

#include "controlFunction.h"

typedef std::vector<double> array_t;

struct datapack_t
{
	array_t t;
	array_t y;
	array_t d;
	array_t u;
	MSGPACK_DEFINE_MAP(t, y, d, u);
};

typedef std::vector<double> state_type;
void equations(state_type &y, double _x);
double model_input = 0.0;
double model_output = 0.0;
//double omega1 = 0.0628;
//double w_star = 50 * 2 * M_PI;
//double v[2] = {0,2};
//double Gamma[2] = {1,0};
//double S[2][2] = {
  //  {std::cos(omega1), std::sin(omega1)}, {-std::sin(omega1), std::cos(omega1)}
//};
namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
	std::string config_path = "config.yaml";
	if (argc != 1)
	{
		config_path = std::string(argv[1]);
	}

	// read parameter from config file
	YAML::Node config = YAML::LoadFile(config_path);
	const float64 sampleFs = config["sample_fs"].as<float64>();
	const float64 dist_freq = config["dist_freq"].as<float64>();
	const float64 dist_gain = config["dist_gain"].as<float64>();
	const float64 runTime = config["run_time"].as<float64>(); // sec
	const float64 warmUp = config["warm_up"].as<float64>();	  // sec
	// const std::string log_path =  config["log_path"].as<std::string>();
	const std::string log_folder = config["log_folder"].as<std::string>();
	const std::string device_log_path = log_folder + '/' + config["device_log_path"].as<std::string>();


	try
	{
		if (!fs::exists(log_folder))
		{
			if (fs::create_directory(log_folder))
			{
				printf("created log dir %s \n", log_folder.c_str());
			}
			else
			{
				printf("Failed to create log dir!\n");
				return 1;
			}
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << std::endl;
	}



	printf("config file path: %s \n", config_path.c_str());
	printf("log file path: %s \n", device_log_path.c_str());

	const int32 samplesPerChan = 1;

	const int32 writeNumChan = 2;
	float64 write_origin[samplesPerChan * writeNumChan + 100] = {0};
	float64 *writeChan0 = &write_origin[0], *writeChan1 = &write_origin[samplesPerChan];

	const int32 totalNumSamples = (runTime + warmUp) * sampleFs + sampleFs;

	std::cout << "run time: " << runTime << " "
			  << "total data points: " << totalNumSamples << " "
			  << "Hz: " << sampleFs / samplesPerChan << std::endl;

	int64 index = 0;
	float64 globalTime = 0;
	void *ctrl_ptr = controllerInit(argc, argv);

	std::ofstream stream(device_log_path, std::ios::binary);
	msgpack::packer<std::ofstream> packer(stream);
	datapack_t data;
	data.t.reserve(totalNumSamples);
	data.d.reserve(totalNumSamples);
	data.u.reserve(totalNumSamples);
	data.y.reserve(totalNumSamples);

	//boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
	state_type x(2);
	x = {0.0, 0.0};

	printf("start warm up\n");
	for (int32 i = 0; i < warmUp * sampleFs; i++)
	{
		//float64 d = 0;
        //double result[2] = {0};
        //for (int m = 0; m < 2; ++m) {
          //  for (int n = 0; n < 2; ++n) {
            //    result[m] += S[m][n] * v[n];
          //  }
        //}
        //for (int m = 0; m < 2; ++m) {
          //  d += Gamma[m] * result[m];
        //}
        //for (int m = 0; m < 2; ++m) {
          //  v[m] = result[m];
        //}
		float64 d = dist_gain * sin(globalTime * dist_freq * 2 * M_PI);
		float64 u = 0;
		float64 y = model_output;

		writeChan0[0] = d;
		writeChan1[0] = u;

		index += 1;
		globalTime += 1 / sampleFs;

		data.t.push_back(globalTime);
		data.d.push_back(d);
		data.y.push_back(y);
		data.u.push_back(u);

		model_input = u + d;
		//stepper.do_step(equations, x, globalTime, 1 / sampleFs);

	    Eigen::Vector2d x;

	    //Eigen::MatrixXd A(2, 2), B(2, 1), C(1, 2);
	    //A << 1.4, -0.728, 0.728, 0;
	    //B << 1, 0;
		//C << 2, -3.0219;
		Eigen::MatrixXd A(2, 2), B(2, 1), C(1, 2);
	    A << 0.999996, 0.00019982,
		-0.004, 0.9982;
	    B << 0.0000019988, 0.02;
	    C << -10, 10;


	    double ud = model_input;
	
	    x = A * x + ud * B;	
	    model_output = (C * x).value();

		if (i >= warmUp * sampleFs)
		{
			break;
		}
	}

	printf("start run main loop\n");
	for (int i = 0; i < runTime * sampleFs; i++)
	{
	//	float64 d = 0;
      //  double result[2] = {0};
      //  for (int m = 0; m < 2; ++m) {
     //       for (int n = 0; n < 2; ++n) {
      //          result[m] += S[m][n] * v[n];
      //      }
     //   }
     //   for (int m = 0; m < 2; ++m) {
       //     d += Gamma[m] * result[m];
      //  }
      //  for (int m = 0; m < 2; ++m) {
      //      v[m] = result[m];
     //   }
		float64 d = dist_gain * sin(globalTime * dist_freq * 2 * M_PI);
		float64 y = model_output;
		float64 u = controllerCompute(ctrl_ptr, y);

		 u = std::min(std::max(u, -4.0), 4.0);

		writeChan0[0] = d;
		writeChan1[0] = u;

		index += 1;
		globalTime += 1 / sampleFs;

		data.t.push_back(globalTime);
		data.d.push_back(d);
		data.y.push_back(y);
		data.u.push_back(u);

		model_input = u + d;
		//stepper.do_step(equations, x, globalTime, 1 / sampleFs);


	    Eigen::Vector2d x;

	    //Eigen::MatrixXd A(2, 2), B(2, 1), C(1, 2);
	    //A << 1.4, -0.728,
		//0.728, 0;
	    //B << 1, 0;
	    //C << 2, -3.0219;

		Eigen::MatrixXd A(2, 2), B(2, 1), C(1, 2);
	    A << 0.999996, 0.00019982,
		-0.004, 0.9982;
	    B << 0.0000019988, 0.02;
	    C << -10, 10;




	    double ud = model_input;
	
	    x = A * x + B * ud;	
	    model_output = (C * x).value();
	    Eigen::MatrixXd G(2,1);
	    G << 1, 0;

		if (i >= runTime * sampleFs)
		{
			break;
		}
	}

	controllerFinish(ctrl_ptr);

	printf("End of program\n");
	packer.pack(data);
	stream.close();

	return 0;
}