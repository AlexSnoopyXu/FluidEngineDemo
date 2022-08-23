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

}

const Vector3f* const ParticleSystemData3::velocities() const
{

}

const Vector3f* const ParticleSystemData3::forces() const
{

}

void ParticleSystemData3::addParticle(const Vector3f& newPosition, const Vector3f& newVelocity = Vector3f(), const Vector3f& newForce = Vector3f())
{

}

void ParticleSystemData3::addParticles(const VectorArray& newPositions, const VectorArray& newVelocities = VectorArray(), const VectorArray& newForce = VectorArray())
{

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

}

void ParticleSystemSolver3::resolveCollision()
{

}

void ParticleSystemSolver3::beginAdvanceTimeStep()
{

}

void ParticleSystemSolver3::endAdvanceTimeStep()
{

}

void ParticleSystemSolver3::timeIntegration(double timeIntervalInSeconds)
{

}