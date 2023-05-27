#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>

#include "NIDAQmx.h"
#include "msgpack.hpp"

#define DAQmxErrChk(functionCall)            \
	if (DAQmxFailed(error = (functionCall))) \
		goto Error;                          \
	else

#define PI 3.1415926535
inline float64 MIN_MAX(float64 val, float64 min, float64 max)
{
	if (val > max)
	{
		return max;
	}
	else if (val < min)
	{
		return min;
	}
	return val;
}

typedef std::vector<double> array_t;

struct datapack_t
{
	array_t t;
	array_t y;
	array_t d;
	array_t u;
	// array_t u_a;
	// array_t y_d;
	// array_t eta1;
	// array_t eta2;
	// array_t xi1;
	// array_t xi2;
	// array_t zeta1;
	// array_t zeta2;
	// MSGPACK_DEFINE_MAP(t, y, d, u, u_a, y_d, eta1, eta2, xi1, xi2, zeta1, zeta2);
	MSGPACK_DEFINE_MAP(t, y, d, u);

};

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);

class AFCController
{
public:
	typedef std::vector<double> state_type;
	float64 out, tw, w_star_, dt_, bound_;

	Eigen::Matrix<double, 2, 1> G;
	Eigen::Matrix2d S;

	Eigen::Vector2d theta_hat;
	float64 alpha, epsilon, gamma, d1, d2;

	float64 u_a, y_d;
	Eigen::Vector2d eta_, xi_, zeta_;

	state_type w;
	boost::numeric::odeint::runge_kutta4<state_type> stepper;

	AFCController(float64 w_star, float64 dt, float64 bound = 3.95) : w(8), w_star_(w_star),
																	  dt_(dt),
																	  bound_(bound)
	{
		tw = 0.0;
		out = 0.0;

		G << 1, 0;
		S << 0, w_star_,
			-w_star_, 0;

		theta_hat << -1.0, 1.0;

		alpha = 3.0;
		epsilon = 3.0;
		gamma = 1.0;
		d1 = 0.1;
		d2 = 3;

		eta_ = theta_hat;

		w[6] = theta_hat(0);
		w[7] = theta_hat(1);
	}

	Eigen::Vector2d exoCpyDynamics(const Eigen::Vector2d &w_hat, const float64 u_a) const
	{
		return S * w_hat + G * (u_a);
	}

	float64 ydiffCalc(const Eigen::Vector2d &zeta_hat) const
	{
		return G.transpose() * zeta_hat - out;
	}

	float64 uaCalc(const Eigen::Vector2d &zeta_hat) const
	{
		return -epsilon * theta_hat.transpose() * zeta_hat;
	}

	Eigen::Vector2d zetaDynamics(const Eigen::Vector2d &zeta_hat) const
	{
		return S * zeta_hat + theta_hat * uaCalc(zeta_hat) - alpha * G * ydiffCalc(zeta_hat);
	}

	Eigen::Vector2d xiDynamics(const Eigen::Vector2d &xi) const
	{
		Eigen::Matrix2d Malpha;
		Malpha << alpha, 0, 0, 0;
		return (S - Malpha).transpose() * xi + G * u_a;
	}

	Eigen::Vector2d etaDynamics(const Eigen::Vector2d &eta, const Eigen::Vector2d &xi) const
	{
		return -gamma * xi * (y_d - (theta_hat - eta).transpose() * xi);
	}

	void dot_projection(const Eigen::Vector2d &x, float64 bound, Eigen::Vector2d &dot_x) const
	{
		dot_x = (x.norm() > d2 || (x.norm() == d2 && dot_x.dot(x) > 0)) ? (dot_x - x * x.transpose() / (x.norm() * x.norm()) * dot_x) : dot_x;
	}

	void projection(Eigen::Vector2d &x, float64 bound) const
	{
		x = x.norm() > bound ? x / x.norm() * bound : x;
	}

	void equations(const state_type &y, state_type &dy, double _x)
	{
		Eigen::Vector2d w_hat, dot_w_hat, zeta_hat, dot_zeta_hat, eta, dot_eta, xi, dot_xi;
		w_hat << y[0], y[1];
		zeta_hat << y[2], y[3];
		xi << y[4], y[5];
		eta << y[6], y[7];

		eta_ = eta;
		xi_ = xi;
		zeta_ = zeta_hat;

		u_a = uaCalc(zeta_hat);
		y_d = ydiffCalc(zeta_hat);
		dot_zeta_hat = zetaDynamics(zeta_hat);
		dot_w_hat = exoCpyDynamics(w_hat, u_a);
		dot_xi = xiDynamics(xi);
		dot_eta = etaDynamics(eta, xi);

		// dot_projection(w_hat, bound_, dot_w_hat);
		// dot_projection(eta, d2, dot_eta);

		dy[0] = dot_w_hat(0);
		dy[1] = dot_w_hat(1);
		dy[2] = dot_zeta_hat(0);
		dy[3] = dot_zeta_hat(1);
		dy[4] = dot_xi(0);
		dy[5] = dot_xi(1);
		dy[6] = dot_eta(0);
		dy[7] = dot_eta(1);
	}

	void setInput(float64 y) { out = y; }

	float64 computeOutput()
	{
		auto func = std::bind(&AFCController::equations, &(*this), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		stepper.do_step(func, w, tw, dt_);

		Eigen::Vector2d w_hat, eta;
		w_hat << w[0], w[1];
		eta << w[6], w[7];
		projection(w_hat, bound_);
		projection(eta, d2);
		w[0] = w_hat(0);
		w[1] = w_hat(1);
		w[6] = eta(0);
		w[7] = eta(1);

		return w[0];
		// return MIN_MAX(w[0], -4.0, 4.0);
	}
};

void afcController(void *ptr_controller, const float64 &y, float64 &u)
{
	AFCController *ptr = (AFCController *)(ptr_controller);
	ptr->setInput(y);
	u = ptr->computeOutput();
}

void (*controlTaskRun)(const float64 &y, float64 &u) = NULL;

float64 input2dB(float64 data) { return 20 * log10f64(data * 1000 / 50 / (2 * 1e-5)); }

void controlTaskInit()
{
}

int main(void)
{
	int32 error = 0;
	char errBuff[2048] = {'\0'};
	TaskHandle taskHandle_w = 0, taskHandle_r = 0;

	const int32 samplesPerChan = 1000;
	const float64 sampleFs = 25000.0, dist_freq = 35.0;

	const int32 writeNumChan = 2;
	float64 write_origin[samplesPerChan * writeNumChan] = {0};
	float64 *writeChan0 = &write_origin[0], *writeChan1 = &write_origin[samplesPerChan];

	const int32 readNumChan = 1;
	float64 read_origin[samplesPerChan * readNumChan] = {0};
	float64 *readChan0 = &read_origin[0];

	const float runTime = 20, warmUp = 1.0; // sec
	const int32 totalNumSamples = (runTime + warmUp) * sampleFs + sampleFs;
	std::cout << "run time: " << runTime << " "
			  << "total data points: " << totalNumSamples << std::endl;

	std::ofstream stream("distRejTemp.bin", std::ios::binary);
	msgpack::packer<std::ofstream> packer(stream);
	datapack_t data;
	data.d.reserve(totalNumSamples);
	data.t.reserve(totalNumSamples);
	data.u.reserve(totalNumSamples);
	data.y.reserve(totalNumSamples);
	// data.u_a.reserve(totalNumSamples);
	// data.y_d.reserve(totalNumSamples);
	// data.eta1.reserve(totalNumSamples);
	// data.eta2.reserve(totalNumSamples);
	// data.xi1.reserve(totalNumSamples);
	// data.xi2.reserve(totalNumSamples);
	// data.zeta1.reserve(totalNumSamples);
	// data.zeta2.reserve(totalNumSamples);


	int64 index = 0;
	float64 globalTime = 0;

	AFCController controller(dist_freq * 2 * PI, 1.0 / sampleFs);

	/*********************************************/
	// DAQmx Configure Code
	/*********************************************/
	DAQmxErrChk(DAQmxCreateTask("", &taskHandle_w));
	DAQmxErrChk(DAQmxCreateAOVoltageChan(taskHandle_w, "Mod1/ao0:1", "", -4.0, 4.0, DAQmx_Val_Volts, NULL));
	DAQmxErrChk(DAQmxCfgSampClkTiming(taskHandle_w, "OnboardClock", sampleFs, DAQmx_Val_Rising, DAQmx_Val_ContSamps, samplesPerChan));
	DAQmxErrChk(DAQmxRegisterDoneEvent(taskHandle_w, 0, DoneCallback, NULL));

	DAQmxErrChk(DAQmxCreateTask("", &taskHandle_r));
	DAQmxErrChk(DAQmxCreateAIVoltageChan(taskHandle_r, "Mod2/ai0", "", DAQmx_Val_PseudoDiff, -4.0, 4.0, DAQmx_Val_Volts, NULL));
	DAQmxErrChk(DAQmxCfgSampClkTiming(taskHandle_r, "OnboardClock", sampleFs, DAQmx_Val_Rising, DAQmx_Val_ContSamps, samplesPerChan));

	/*********************************************/
	// DAQmx Write Code
	/*********************************************/
	DAQmxErrChk(DAQmxWriteAnalogF64(taskHandle_w, samplesPerChan, 0, 10.0, DAQmx_Val_GroupByChannel, write_origin, NULL, NULL));

	/*********************************************/
	// DAQmx Start Code
	/*********************************************/
	DAQmxErrChk(DAQmxStartTask(taskHandle_w));
	DAQmxErrChk(DAQmxStartTask(taskHandle_r));

	for (int i = 0; i < warmUp * sampleFs / samplesPerChan; i++)
	{

		for (int i = 0; i < samplesPerChan; i++)
		{

			float64 d = 2.0 * sin(globalTime * dist_freq * 2 * PI);
			float64 u = 0;
			float64 y = (read_origin[i] - 0.00025) * 8000;

			data.t.push_back(globalTime);
			data.d.push_back(d);
			data.u.push_back(u);
			data.y.push_back(y);
			// data.u_a.push_back(controller.u_a);
			// data.y_d.push_back(controller.y_d);
			// data.eta1.push_back(controller.eta_(0));
			// data.eta2.push_back(controller.eta_(1));
			// data.xi1.push_back(controller.xi_(0));
			// data.xi2.push_back(controller.xi_(1));
			// data.zeta1.push_back(controller.zeta_(0));
			// data.zeta2.push_back(controller.zeta_(1));


			writeChan0[i] = d;
			writeChan1[i] = u;

			index += 1;
			globalTime += 1 / sampleFs;
		}

		DAQmxErrChk(DAQmxWriteAnalogF64(taskHandle_w, samplesPerChan, 0, 10.0, DAQmx_Val_GroupByChannel, write_origin, NULL, NULL));
		DAQmxErrChk(DAQmxReadAnalogF64(taskHandle_r, samplesPerChan, 10.0, DAQmx_Val_GroupByChannel, read_origin, 1000, NULL, NULL));
	}

	for (int i = 0; i < runTime * sampleFs / samplesPerChan; i++)
	{

		for (int i = 0; i < samplesPerChan; i++)
		{

			float64 d = 2.0 * sin(globalTime * dist_freq * 2 * PI);
			float64 u = 0;
			float64 y = (read_origin[i] - 0.00025) * 8000;

			afcController(&controller, y, u);

			data.t.push_back(globalTime);
			data.d.push_back(d);
			data.u.push_back(u);
			data.y.push_back(y);
			// data.u_a.push_back(controller.u_a);
			// data.y_d.push_back(controller.y_d);
			// data.eta1.push_back(controller.eta_(0));
			// data.eta2.push_back(controller.eta_(1));
			// data.xi1.push_back(controller.xi_(0));
			// data.xi2.push_back(controller.xi_(1));
			// data.zeta1.push_back(controller.zeta_(0));
			// data.zeta2.push_back(controller.zeta_(1));

			writeChan0[i] = d;
			writeChan1[i] = u;

			index += 1;
			globalTime += 1 / sampleFs;
		}

		DAQmxErrChk(DAQmxWriteAnalogF64(taskHandle_w, samplesPerChan, 0, 10.0, DAQmx_Val_GroupByChannel, write_origin, NULL, NULL));
		DAQmxErrChk(DAQmxReadAnalogF64(taskHandle_r, samplesPerChan, 10.0, DAQmx_Val_GroupByChannel, read_origin, 1000, NULL, NULL));
	}

	// printf("Generating voltage continuously. Press Enter to interrupt\n");
	// getchar();

Error:
	if (DAQmxFailed(error))
		DAQmxGetExtendedErrorInfo(errBuff, 2048);
	if (taskHandle_w != 0)
	{
		/*********************************************/
		// DAQmx Stop Code
		/*********************************************/
		DAQmxStopTask(taskHandle_w);
		DAQmxClearTask(taskHandle_w);
	}
	if (DAQmxFailed(error))
		printf("DAQmx Error: %s\n", errBuff);
	printf("End of program, press Enter key to quit\n");
	getchar();
	packer.pack(data);
	stream.close();
	return 0;
}

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData)
{
	int32 error = 0;
	char errBuff[2048] = {'\0'};

	// Check to see if an error stopped the task.
	DAQmxErrChk(status);

Error:
	if (DAQmxFailed(error))
	{
		DAQmxGetExtendedErrorInfo(errBuff, 2048);
		DAQmxClearTask(taskHandle);
		printf("DAQmx Error: %s\n", errBuff);
	}
	return 0;
}
