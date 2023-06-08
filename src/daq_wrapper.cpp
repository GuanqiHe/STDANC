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

#include "controlFunction.h"
#include "logger/data_logging.hpp"

using namespace std::chrono;

struct DAQwrapper_t
{
	const int32 sampleFs;
	const int32 samplesPerChan; // 1
	const float64 timeout;

	TaskHandle writeTaskHandle;
	const int32 writeNumChan;
	float64 writeData[2];
	float64 *writeChannel0;
	float64 *writeChannel1;

	TaskHandle readTaskHandle;
	const int32 readNumChan;
	float64 readData[2];
	float64 *readChannel0;
	float64 *readChannel1;

	int32 error;

	DAQwrapper_t(int32 sample_frequency) : sampleFs(sample_frequency),
										   samplesPerChan(1),
										   timeout(10.0 / sampleFs),

										   writeTaskHandle(0),
										   writeNumChan(2),
										   writeData{0, 0},
										   writeChannel0(&writeData[0]),
										   writeChannel1(&writeData[1]),

										   readTaskHandle(0),
										   readNumChan(1),
										   readData{0, 0},
										   readChannel0(&readData[0]),
										   readChannel1(&readData[1]),

										   error(0)
	{
	}
};

class DAQmxFailedException : public std::exception
{
public:
	DAQwrapper_t *wrapper;
	DAQmxFailedException(DAQwrapper_t *wrapper_) : std::exception(), wrapper(wrapper_) {}
	virtual ~DAQmxFailedException() throw() {}
};

inline void DAQmxErrChk(DAQwrapper_t *wrapper, int functionCall)
{
	if (DAQmxFailed(wrapper->error = functionCall))
	{
		throw DAQmxFailedException(wrapper);
	}
}

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);

void DAQmxFailedErrorHandler(TaskHandle taskHandle, int32 error);

void *hardwareInit(int sampleFs, double startWait);

void hardwareRun(void *wrapper_);

void hardwareFinish(void *wrapper_);

inline int64 timeConversion(int64 t) { return t + 2082844800; }

int main(int argc, char *argv[])
{
	std::string config_path = "config.yaml";
	if (argc != 1)
	{
		config_path = std::string(argv[1]);
	}

	// read parameter from config file
	YAML::Node config = YAML::LoadFile(config_path);
	const int32 sampleFs = config["sample_fs"].as<float64>();
	const float64 dist_freq = config["dist_freq"].as<float64>();
	const float64 dist_gain = config["dist_gain"].as<float64>();
	const float64 runTime = config["run_time"].as<float64>();	  // sec
	const float64 warmUp = config["warm_up"].as<float64>();		  // sec
	const float64 startWait = config["start_wait"].as<float64>(); // sec
	const float64 analog_read_bias = config["analog_read_bias"].as<float64>();
	const float64 analog_read_gain = config["analog_read_gain"].as<float64>();
	// const std::string log_path =  config["log_path"].as<std::string>();
	const int sample_len = (config["run_time"].as<double>() + config["warm_up"].as<double>()) * config["sample_fs"].as<double>() + 1;
	const int32 totalNumSamples = (runTime + warmUp) * sampleFs + sampleFs;

	std::cout << "run time: " << runTime << " "
			  << "total data points: " << totalNumSamples << " "
			  << "Hz: " << sampleFs << std::endl;

	data_logger logger;
	logger.init({"t", "y", "d", "u"}, sample_len, config["data_log_path"].as<std::string>());

	printf("config file path: %s \n", config_path.c_str());
	// printf("log file path: %s \n", log_path.c_str());

	int64 index = 0;
	float64 globalTime = 0;
	void *ctrl_ptr = controllerInit(argc, argv);

	try
	{
		DAQwrapper_t *hardware_wrapper = (DAQwrapper_t *)hardwareInit(sampleFs, startWait);

		_V2::system_clock::time_point start, end;

		start = system_clock::now();

		printf("start warm up\n");
		for (int32 i = 0; i < warmUp * sampleFs; i++)
		{

			float64 d = dist_gain * sin(globalTime * dist_freq * 2 * M_PI);
			float64 u = 0;
			float64 y = *hardware_wrapper->readChannel0 * analog_read_gain + analog_read_bias;

			*hardware_wrapper->writeChannel0 = d;
			*hardware_wrapper->writeChannel1 = u;

			index += 1;
			globalTime += 1 / sampleFs;

			logger.log({globalTime, y, d, u});

			hardwareRun(hardware_wrapper);

			if (i >= warmUp * sampleFs)
			{
				break;
			}
		}
		end = system_clock::now();
		std::cout << "Warm up time: " << double(duration_cast<microseconds>(end - start).count()) * microseconds::period::num / microseconds::period::den << std::endl;

		start = system_clock::now();

		printf("start run main loop\n");
		for (int i = 0; i < runTime * sampleFs; i++)
		{

			float64 d = dist_gain * sin(globalTime * dist_freq * 2 * M_PI);
			float64 y = *hardware_wrapper->readChannel0 * analog_read_gain + analog_read_bias;
			float64 u = controllerCompute(ctrl_ptr, y);

			u = std::min(std::max(u, -4.0), 4.0);

			*hardware_wrapper->writeChannel0 = d;
			*hardware_wrapper->writeChannel1 = u;

			index += 1;
			globalTime += 1 / sampleFs;

			logger.log({globalTime, y, d, u});

			hardwareRun(hardware_wrapper);

			if (i >= runTime * sampleFs)
			{
				break;
			}
		}
		end = system_clock::now();
		std::cout << "Run time: " << double(duration_cast<microseconds>(end - start).count()) * microseconds::period::num / microseconds::period::den << std::endl;

		controllerFinish(ctrl_ptr);

		hardwareFinish(hardware_wrapper);
	}
	catch (DAQmxFailedException &e)
	{
		DAQmxFailedErrorHandler(e.wrapper->readTaskHandle, e.wrapper->error);
		DAQmxFailedErrorHandler(e.wrapper->readTaskHandle, e.wrapper->error);
	}

	logger.write();

	return 0;
}

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData)
{
	int32 error = 0;

	if (DAQmxFailed(error = status))
	{
		DAQmxFailedErrorHandler(taskHandle, error);
	}

	return 0;
}

void DAQmxFailedErrorHandler(TaskHandle taskHandle, int32 error)
{
	char errBuff[2048] = {'\0'};

	if (DAQmxFailed(error))
	{
		DAQmxGetExtendedErrorInfo(errBuff, 2048);
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
		printf("DAQmx Error: %s\n", errBuff);
	}
}

void *hardwareInit(int sampleFs, double startWait)
{

	DAQwrapper_t *wrapper = new DAQwrapper_t(sampleFs);

	CVIAbsoluteTime t;

	/*********************************************/
	// DAQmx Configure Code
	/*********************************************/
	DAQmxErrChk(wrapper, DAQmxCreateTask("write", &wrapper->writeTaskHandle));
	DAQmxErrChk(wrapper, DAQmxCreateAOVoltageChan(wrapper->writeTaskHandle, "Mod1/ao0:1", "", -4.0, 4.0, DAQmx_Val_Volts, NULL));
	DAQmxErrChk(wrapper, DAQmxCfgSampClkTiming(wrapper->writeTaskHandle, "OnboardClock", sampleFs, DAQmx_Val_Rising, DAQmx_Val_ContSamps, wrapper->samplesPerChan));
	DAQmxErrChk(wrapper, DAQmxRegisterDoneEvent(wrapper->writeTaskHandle, 0, DoneCallback, NULL));

	DAQmxErrChk(wrapper, DAQmxCreateTask("read", &wrapper->readTaskHandle));
	DAQmxErrChk(wrapper, DAQmxCreateAIVoltageChan(wrapper->readTaskHandle, "Mod2/ai0", "", DAQmx_Val_PseudoDiff, -4.0, 4.0, DAQmx_Val_Volts, NULL));
	DAQmxErrChk(wrapper, DAQmxCfgSampClkTiming(wrapper->readTaskHandle, "OnboardClock", sampleFs, DAQmx_Val_Rising, DAQmx_Val_ContSamps, wrapper->samplesPerChan));
	DAQmxErrChk(wrapper, DAQmxRegisterDoneEvent(wrapper->readTaskHandle, 0, DoneCallback, NULL));

	t.cviTime.lsb = 0;
	t.cviTime.msb = timeConversion(time(NULL) + int(startWait));
	DAQmxCfgTimeStartTrig(wrapper->writeTaskHandle, t, DAQmx_Val_HostTime);
	DAQmxCfgTimeStartTrig(wrapper->readTaskHandle, t, DAQmx_Val_HostTime);

	/*********************************************/
	// DAQmx Write Code
	/*********************************************/
	float64 firstWrite[100] = {0};
	DAQmxErrChk(wrapper, DAQmxWriteAnalogF64(wrapper->writeTaskHandle, 50, 0, 10.0, DAQmx_Val_GroupByChannel, firstWrite, NULL, NULL));
	// DAQmxErrChk(wrapper, DAQmxWriteAnalogScalarF64(wrapper->writeTaskHandle, 0, 10.0, 0.0, NULL));

	/*********************************************/
	// DAQmx Start Code
	/*********************************************/
	DAQmxErrChk(wrapper, DAQmxStartTask(wrapper->writeTaskHandle));
	DAQmxErrChk(wrapper, DAQmxStartTask(wrapper->readTaskHandle));

	DAQmxErrChk(wrapper, DAQmxReadAnalogScalarF64(wrapper->readTaskHandle, 10.0, wrapper->readChannel0, NULL));

	return (void *)(wrapper);
}

void hardwareRun(void *wrapper_)
{
	DAQwrapper_t *wrapper = (DAQwrapper_t *)(wrapper_);
	DAQmxErrChk(wrapper, DAQmxReadAnalogScalarF64(wrapper->readTaskHandle, wrapper->timeout, wrapper->readChannel0, NULL));
	DAQmxErrChk(wrapper, DAQmxWriteAnalogF64(wrapper->writeTaskHandle, wrapper->samplesPerChan, 0, wrapper->timeout, DAQmx_Val_GroupByChannel, wrapper->writeData, NULL, NULL));
}

void hardwareFinish(void *wrapper_)
{
	DAQwrapper_t *wrapper = (DAQwrapper_t *)(wrapper_);

	DAQmxStopTask(wrapper->writeTaskHandle);
	DAQmxClearTask(wrapper->writeTaskHandle);
	DAQmxStopTask(wrapper->readTaskHandle);
	DAQmxClearTask(wrapper->readTaskHandle);

	printf("End of program\n");

	delete wrapper;
}
