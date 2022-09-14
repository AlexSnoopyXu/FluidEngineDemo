#include "SpringAnimation.h"

void SimpleMassSpringAnimation::onAdvanceTimStep(double timeIntervalInSeconds)
{
	size_t numberOfPoints = positions.size();
	size_t numberOfEdges = edges.size();


	for (size_t i = 0; i < numberOfPoints; i++)
	{
		// compute gravity forces
		forces[i] = _mass * gravity;
		// compute drag forces
		forces[i] += -_dragCofficient * velocities[i];
	}

	for (size_t i = 0; i < numberOfEdges; i++)
	{
		// compute spring forces
		size_t pointIndex0 = edges[i].first;
		size_t pointIndex1 = edges[i].second;
		Vector3f pos0 = positions[pointIndex0];
		Vector3f pos1 = positions[pointIndex1];
		Vector3f r = pos0 - pos1;
		float distance = r.Length();
		if (distance)
		{
			Vector3f force = -_stiffness * (distance - _restLength) * r.Normalized();
			forces[pointIndex0] += force;
			forces[pointIndex1] -= force;
		}

		// compute damping forces
		Vector3f vel0 = velocities[pointIndex0];
		Vector3f vel1 = velocities[pointIndex1];
		Vector3f damping = -_dampingCofficient * (vel0 - vel1);
		forces[pointIndex0] += damping;
		forces[pointIndex1] -= damping;
	}

	// upate states
	for (size_t i = 0; i < numberOfPoints; i++)
	{
		Vector3f newAcceleration = forces[i] / _mass;
		Vector3f newVelocity = velocities[i] + newAcceleration * timeIntervalInSeconds;
		Vector3f newPosition = positions[i] + newVelocity * timeIntervalInSeconds;

		velocities[i] = newVelocity;
		positions[i] = newPosition;
	}

}