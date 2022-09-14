// FluidEngineDemo.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <array>
#include <chrono>
#include <thread>
#include "ParticleSystem.h"

const size_t kBufferSize = 80;
const double M_PI = 3.141592653f;

const char* kGrayScaleTable = " .:-=+*#%@";
const size_t kGrayScaleTableSize = sizeof(kGrayScaleTable) / sizeof(char);

using namespace std;
using namespace chrono;

void updateWaves(const double timeInterval, double& x, double& speed)
{
    x += timeInterval * speed;

    //Boundary reflection
    if (x > 1.f)
    {
        speed *= -1.f;
        x += 1.f + timeInterval * speed;
    }
    else if (x < 0.f)
    {
        speed *= -1.f;
        x = timeInterval * speed;
    }
}

void accumulateWaveToHeightField(const double x, const double waveLength, const double maxHeight, array<double, kBufferSize>& heightField)
{
    const double quarterWaveLength = 0.25f * waveLength;
    const int start = static_cast<int>((x - quarterWaveLength) * kBufferSize);
    const int end = static_cast<int>((x + quarterWaveLength) * kBufferSize);

    for (int i = start; i < end; i++)
    {
        int iNew = i;
        if (i < 0)
        {
            iNew = -i - 1;
        }
        else if (i >= static_cast<int>(kBufferSize))
        {
            iNew = 2 * kBufferSize - i - 1;
        }

        double distance = fabs((i + 0.5f) / kBufferSize - x);
        double height = maxHeight * 0.5f * (cos(min(distance * M_PI / quarterWaveLength, M_PI)) + 1.0);
        heightField[iNew] += height;
    }
}

void draw(const array<double, kBufferSize>& heightField)
{
	string buffer(kBufferSize, ' ');

	// Convert height field to grayscale
	for (size_t i = 0; i < kBufferSize; ++i) {
		double height = heightField[i];
		size_t tableIndex = min(
			static_cast<size_t>(floor(kGrayScaleTableSize * height)),
			kGrayScaleTableSize - 1);
		buffer[i] = kGrayScaleTable[tableIndex];

	}

	// Clear old prints
	for (size_t i = 0; i < kBufferSize; ++i) {
		printf("\b");

	}

	// Draw new buffer
	printf("%s", buffer.c_str());
	fflush(stdout);
}

void DemoDraw()
{
    const double waveLengthX = 0.8f;
    const double waveLengthY = 1.2f;

    const double maxHeightX = 0.5f;
    const double maxHeightY = 0.4f;

    double x = 0.f;
    double y = 1.f;
    double speedX = 1.f;
    double speedY = -0.5f;

    const int fps = 100;
    const double timeIntreval = 1.f / fps;

    array<double, kBufferSize> heightField;

    for (int i = 0; i < 10000; i++)
    {
        // Update waves
        updateWaves(timeIntreval, x, speedX);
        updateWaves(timeIntreval, y, speedY);

        // Clear height field
        for (double& height : heightField) {
            height = 0.0;
        }

        // Accumulate waves for each center point
        accumulateWaveToHeightField(
            x, waveLengthX, maxHeightX, heightField);
        accumulateWaveToHeightField(
            y, waveLengthY, maxHeightY, heightField);

        // Draw height field
        draw(heightField);

        // Wait
        this_thread::sleep_for(milliseconds(1000 / fps));
    }

    printf("\n");
    fflush(stdout);
}

void DemoDraw_PBD()
{
    PciSphSystemSolver3 pciSphSysSolver;

    const int fps = 100;
    const double timeIntreval = 1.f / fps;

}

int main()
{
    //DemoDraw();

    DemoDraw_PBD();

    return 0;
}
