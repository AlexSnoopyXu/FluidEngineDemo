#include "ParticleSystem.h"

ParticleSystemData3::ParticleSystemData3()
{

}

void ParticleSystemData3::resize(size_t newNumberOfParticles)
{

}

size_t ParticleSystemData3::numOfParticles() const
{

}

const Vector3f* const ParticleSystemData3::positions() const
{
	return _positions.cbegin()._Ptr;
}

const Vector3f* const ParticleSystemData3::velocities() const
{
	return _velocities.cbegin()._Ptr;
}

const Vector3f* const ParticleSystemData3::forces() const
{
	return _forces.cbegin()._Ptr;
}

const float ParticleSystemData3::mass() const 
{
	return _mass;
}

const float ParticleSystemData3::radius() const
{
	return _radius;
}

void ParticleSystemData3::addParticle(const Vector3f& newPosition, const Vector3f& newVelocity = Vector3f(), const Vector3f& newForce = Vector3f())
{

}

void ParticleSystemData3::addParticles(const VectorArray& newPositions, const VectorArray& newVelocities = VectorArray(), const VectorArray& newForce = VectorArray())
{

}

void ParticleSystemData3::setMass(float mass)
{
	_mass = mass;
}

/* Solver begin */

ParticleSystemSolver3::ParticleSystemSolver3()
{
	ParticleSystemData3::ParticleSystemData3();
	_particleSystemData = std::make_shared<ParticleSystemData3>();
	_wind = std::make_shared<ConstantVectorField3>();
}

void ParticleSystemSolver3::onAdvanceTimStep(double timeIntervalInSeconds)
{
	beginAdvanceTimeStep();

	accumulateForces(timeIntervalInSeconds);
	timeIntegration(timeIntervalInSeconds);
	resolveCollision();

	endAdvanceTimeStep();
}

void ParticleSystemSolver3::accumulateForces(double timeStepInSeconds)
{
	accumulateExternalForces();
}

void ParticleSystemSolver3::resolveCollision()
{
	auto positions = _particleSystemData->positions();
	auto velocities = _particleSystemData->velocities();

	if (_collider != nullptr)
	{
		size_t n = _particleSystemData->numOfParticles();
		const float radius = _particleSystemData->radius();
		for (size_t i = 0; i < n; i++)
		{
			_collider->resolveCollision(positions[i], velocities[i], radius, _restitutionCoeffcient, _newPositions[i], _newVelocities[i]);
		}
	}

}

void ParticleSystemSolver3::beginAdvanceTimeStep()
{
	size_t n = _particleSystemData->numOfParticles();
	_newPositions.resize(n);
	_newVelocities.resize(n);

	// Clear forces
	auto forces = const_cast<Vector3f*>(_particleSystemData->forces());
	for (size_t i = 0; i < n; i++)
	{
		forces[i] = Vector3f();
	}

	// For subclass
	onBeginAdvanceTimeStep();
}

void ParticleSystemSolver3::endAdvanceTimeStep()
{
	size_t n = _particleSystemData->numOfParticles();

	auto positions = const_cast<Vector3f*>(_particleSystemData->positions());
	auto velocities = const_cast<Vector3f*>(_particleSystemData->velocities());
	for (size_t i = 0; i < n; i++)
	{
		positions[i] = _newPositions[i];
		velocities[i] = _newVelocities[i];
	}

	// For subclass
	onEndAdvanceTimeStep();
}

void ParticleSystemSolver3::timeIntegration(double timeIntervalInSeconds)
{
	size_t n = _particleSystemData->numOfParticles();
	auto positions = _particleSystemData->positions();
	auto forces = _particleSystemData->forces();
	auto velocities = _particleSystemData->velocities();
	const float mass = _particleSystemData->mass();

	for (size_t i = 0; i < n; i++)
	{
		// Integrate velocity first
		Vector3f& newVelocity = _newVelocities[i];
		newVelocity = velocities[i] + timeIntervalInSeconds * forces[i] / mass;

		// Integrate position
		Vector3f& newPosition = _newPositions[i];
		newPosition = positions[i] + newVelocity * timeIntervalInSeconds;
	}
}

void ParticleSystemSolver3::accumulateExternalForces()
{
	size_t n = _particleSystemData->numOfParticles();
	auto positions = _particleSystemData->positions();
	auto forces = const_cast<Vector3f*>(_particleSystemData->forces());
	auto velocities = _particleSystemData->velocities();
	const float mass = _particleSystemData->mass();

	for (size_t i = 0; i < n; i++)
	{
		// Gravity
		Vector3f force = mass * _gravity;

		// Wind force
		Vector3f relativeVel = velocities[i] - _wind->sample(positions[i]);
		force += -_dragCoefficient * relativeVel;

		forces[i] += force;
	}

}