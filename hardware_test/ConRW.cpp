/*********************************************************************
 *
 * ANSI C Example program:
 *    ContGen-ExtClk.c
 *
 * Example Category:
 *    AO
 *
 * Description:
 *    This example demonstrates how to output a continuous periodic
 *    waveform using an external clock.
 *
 * Instructions for Running:
 *    1. Select the Physical Channel to correspond to where your
 *       signal is output on the DAQ device.
 *    2. Enter the Minimum and Maximum Voltage Ranges.
 *    3. Specify the external sample clock source (typically a PFI or
 *       RTSI pin) in the timing function.
 *
 * Steps:
 *    1. Create a task.
 *    2. Create an Analog Output Voltage Channel.
 *    3. Define the parameters for an External Clock Source.
 *       Additionally, define the sample mode to be Continuous.
 *    4. Write the waveform to the output buffer.
 *    5. Call the Start function.
 *    6. Wait until the user presses the Stop button.
 *    7. Call the Clear Task function to clear the Task.
 *    8. Display an error if any.
 *
 * I/O Connections Overview:
 *    Make sure your signal output terminal matches the Physical
 *    Channel I/O Control. Also, make sure your external clock
 *    terminal matches the Clock Source Control. For further
 *    connection information, refer to your hardware reference manual.
 *
 *********************************************************************/
#include "NIDAQmx.h"
#include <cstdio>
#include <cmath>
#include "msgpack.hpp"
#include <iostream>
#include <fstream>

#define DAQmxErrChk(functionCall)            \
	if (DAQmxFailed(error = (functionCall))) \
		goto Error;                          \
	else

#define PI 3.1415926535

typedef float64 array_t[10000];

struct datapack_t
{
	// array_t u;
	// array_t d;
	array_t t;
	array_t y;
	MSGPACK_DEFINE_MAP(t, y);
};

datapack_t save;

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);

int main(void)
{
	std::ofstream stream("save.bin", std::ios::binary);
    msgpack::packer<std::ofstream> packer(stream);
	
	
	int32 error = 0;
	TaskHandle taskHandle_w = 0, taskHandle_r = 0;
	const int32 samplesPerChan = 1000, numChan = 2;
	float64 data[samplesPerChan * numChan] = {0};
	float64 data_init[samplesPerChan * numChan] = {0};
	char errBuff[2048] = {'\0'};
	float64 sampleFs = 25000.0, frequency = 110.0;

	int i_c = 0, i_n = 0, t_i = 0;

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
	DAQmxErrChk(DAQmxWriteAnalogF64(taskHandle_w, samplesPerChan, 0, 10.0, DAQmx_Val_GroupByChannel, data_init, NULL, NULL));

	/*********************************************/
	// DAQmx Start Code
	/*********************************************/
	DAQmxErrChk(DAQmxStartTask(taskHandle_w));
	DAQmxErrChk(DAQmxStartTask(taskHandle_r));

	for (int i = 0; i < 10; i++)
	{
		for (int i = 0; i < samplesPerChan; i++)
		{
			data[i] = 2.0 * sin((double)(i_c + i) * frequency * PI / sampleFs);
			data[samplesPerChan + i] = 2.0 * sin((double)(i_c + i) * frequency * PI / sampleFs);
			save.t[t_i++] = t_i / sampleFs;
		}
		i_c += samplesPerChan;

		DAQmxErrChk(DAQmxWriteAnalogF64(taskHandle_w, samplesPerChan, 0, 10.0, DAQmx_Val_GroupByChannel, data, NULL, NULL));
		DAQmxErrChk(DAQmxReadAnalogF64(taskHandle_r, samplesPerChan, 10.0, DAQmx_Val_GroupByChannel, &save.y[i*1000], 1000, NULL, NULL));
	}

    packer.pack(save);
    stream.close();
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
