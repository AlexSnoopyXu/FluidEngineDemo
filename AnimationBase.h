#pragma once
#include <iostream>
#include "MathFoundationLib.h"
#include "MathGeoLib-1.5/src/Math/float3.h"

using namespace std;

struct Frame final
{
	unsigned int index = 0;
	double timeIntervalInSeconds = 1 / 60.0;

	double timeInSeconds()
	{
		return index * timeIntervalInSeconds;
	}

	void advance()
	{
		++index;
	}

	void advance(unsigned int delta)
	{
		index += delta;
	}
};

class AnimationBase
{
public:
	void update(const Frame& frame);

protected:
	virtual void onUpdate(const Frame& frame) = 0;
};

class PhysicsAnimation : public AnimationBase
{
protected:
	virtual void onAdvanceTimStep(double timeIntervalInSeconds) = 0;

private:
	virtual void onUpdate(const Frame& frame) final override;
	void advanceTimeStep(double timeIntervalInSeconds);

private:
	Frame _currentFrame;

};