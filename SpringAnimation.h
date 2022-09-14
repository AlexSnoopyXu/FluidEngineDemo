#pragma once

#include "AnimationBase.h"

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

private:
	const Vector3f gravity = Vector3f(0.0f, 9.8f, 0.0f);
	const float _mass = 0.5f;
	const float _restLength = 1.0f;
	const float _stiffness = 500.f;
	const float _dampingCofficient = 1.0f;
	const float _dragCofficient = 0.1;
};

