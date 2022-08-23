#pragma once
#include <iostream>
#include "MathFoundationLib.h"
#include <vector>
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

class SimpleMassSpringAnimation : public PhysicsAnimation
{
public:
	struct Edge
	{
		size_t first;
		size_t second;
	};

	vector<Vector3f> positions;
	vector<Vector3f> velocities;
	vector<Vector3f> forces;
	vector<Edge> edges;

	SimpleMassSpringAnimation(size_t numOfPoints = 10)
	{
		size_t numOfEdges = numOfPoints - 1;
		positions.resize(numOfPoints);
		velocities.resize(numOfPoints);
		forces.resize(numOfPoints);
		edges.resize(numOfEdges);

		for (size_t i = 0; i < numOfPoints; i++)
		{
			positions[i].x = static_cast<float>(i);
		}

		for (size_t i = 0; i < numOfEdges; i++)
		{
			edges[i] = { i, i + 1 };
		}
	}

protected:
	virtual void onAdvanceTimStep(double timeIntervalInSeconds) override;
};

