#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>

#include "NIDAQmx.h"
#include "yaml.h"
#include "msgpack.hpp"
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

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
void equations(const state_type &y, state_type &dy, double _x);
double model_input = 0.0;
double model_output = 0.0;
double w_star = 50 * 2 * M_PI;

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
	const float64 runTime = config["run_time"].as<float64>();	  // sec
	const float64 warmUp = config["warm_up"].as<float64>();		  // sec
	const float64 startWait = config["start_wait"].as<float64>(); // sec
	const float64 analog_read_bias = config["analog_read_bias"].as<float64>();
	const float64 analog_read_gain = config["analog_read_gain"].as<float64>();
	// const std::string log_path =  config["log_path"].as<std::string>();

	printf("config file path: %s \n", config_path.c_str());
	// printf("log file path: %s \n", log_path.c_str());

	int32 error = 0;
	char errBuff[2048] = {'\0'};
	TaskHandle taskHandle_w = 0, taskHandle_r = 0;
	char trigName[256] = {0};

	const int32 samplesPerChan = 1;

	const int32 writeNumChan = 2;
	float64 write_origin[samplesPerChan * writeNumChan + 100] = {0};
	float64 *writeChan0 = &write_origin[0], *writeChan1 = &write_origin[samplesPerChan];

	const int32 readNumChan = 1;
	float64 read_origin[samplesPerChan * readNumChan] = {0};
	float64 *readChan0 = &read_origin[0];

	const int32 totalNumSamples = (runTime + warmUp) * sampleFs + sampleFs;

	std::cout << "run time: " << runTime << " "
			  << "total data points: " << totalNumSamples << " "
			  << "Hz: " << sampleFs / samplesPerChan << std::endl;

	int64 index = 0;
	float64 globalTime = 0;
	float64 timeout = 10.0 / sampleFs;
	void *ctrl_ptr = controllerInit(argc, argv);

	std::ofstream stream(config["data_log_path"].as<std::string>(), std::ios::binary);
	msgpack::packer<std::ofstream> packer(stream);
	datapack_t data;
	data.t.reserve(totalNumSamples);
	data.d.reserve(totalNumSamples);
	data.u.reserve(totalNumSamples);
	data.y.reserve(totalNumSamples);


    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
	state_type x(2);
	x = {0.0, 0.0};

	printf("start warm up\n");
	for (int32 i = 0; i < warmUp * sampleFs; i++)
	{

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
		stepper.do_step(equations, x, globalTime, 1 / sampleFs);

		if (i >= warmUp * sampleFs)
		{
			break;
		}
	}

	printf("start run main loop\n");
	for (int i = 0; i < runTime * sampleFs; i++)
	{

		float64 d = dist_gain * sin(globalTime * dist_freq * 2 * M_PI);
		float64 y = model_output;
		float64 u = controllerCompute(ctrl_ptr, y);

		// u = std::min(std::max(u, -4.0), 4.0);

		writeChan0[0] = d;
		writeChan1[0] = u;

		index += 1;
		globalTime += 1 / sampleFs;

		data.t.push_back(globalTime);
		data.d.push_back(d);
		data.y.push_back(y);
		data.u.push_back(u);

		model_input = u + d;
		stepper.do_step(equations, x, globalTime, 1 / sampleFs);

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

void equations(const state_type &y, state_type &dy, double _x)
{
	Eigen::Vector2d x;
	x << y[0], y[1];

	Eigen::MatrixXd A(2, 2), B(2, 1), C(1, 2);
	A << 0, 1,
		-20, -9;
	B << 0, 100;
	C << -10, 10;

	double u = model_input;
	Eigen::Vector2d dot_x = A * x + B * u;
	model_output = (C * x).value();

	Eigen::MatrixXd G(2, 1);
	G << 1, 0;

	dy[0] = dot_x(0);
	dy[1] = dot_x(1);
}
