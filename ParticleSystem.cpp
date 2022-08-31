#include "ParticleSystem.h"

ParticleSystemData3::ParticleSystemData3()
{

}

void ParticleSystemData3::resize(size_t newNumberOfParticles)
{

}

size_t ParticleSystemData3::numOfParticles() const
{
	return 0;
}

const VectorArray& const ParticleSystemData3::positions() const
{
	return _positions;
}

const VectorArray& const ParticleSystemData3::velocities() const
{
	return _velocities;
}

const VectorArray& const ParticleSystemData3::forces() const
{
	return _forces;
}

const float ParticleSystemData3::mass() const 
{
	return _mass;
}

const float ParticleSystemData3::radius() const
{
	return _radius;
}

void ParticleSystemData3::addParticle(const Vector3f& newPosition, const Vector3f& newVelocity, const Vector3f& newForce)
{

}

void ParticleSystemData3::addParticles(const VectorArray& newPositions, const VectorArray& newVelocities, const VectorArray& newForce)
{

}

void ParticleSystemData3::setMass(float mass)
{
	_mass = mass;
}

void ParticleSystemData3::buildNearbySearcher(float maxSearchRadius)
{
	_pointNeighborSearcher = std::make_shared<PointHashGridSearcher3>(kDefaultHashGridResolution, kDefaultHashGridResolution, kDefaultHashGridResolution, 2 * maxSearchRadius);
	_pointNeighborSearcher->build(positions());
}

void ParticleSystemData3::buildNearbyList(float maxSearchRadius)
{
	_neighborLists.resize(numOfParticles());

	auto points = positions();
	for (size_t i = 0; i < numOfParticles(); i++)
	{
		Vector3f origin = points[i];
		_neighborLists[i].clear();
		_pointNeighborSearcher->forEachNearbyPoint(
			origin,
			maxSearchRadius,
			[&](size_t j, const Vector3f&) {
			if (j != i)
			{
				_neighborLists[i].push_back(j);
			}
		}
		);
	}
}

/* Sph system*/

void SphSystemData3::updateDensities()
{
	auto& p = positions();
	auto& d = densities();

	for (size_t i = 0; i < numOfParticles(); i++)
	{
		d[i] = mass() * sumOfKernelNearby(p[i]);
	}
}

float SphSystemData3::sumOfKernelNearby(const Vector3f& position) const
{
	float sum = 0.f;
	SphStdkernel3 kernel(radius());
	
	neighborSearcher()->forEachNearbyPoint(position, radius(),
		[&](size_t, const Vector3f& neighborPosition) {
		float dist = position.Distance(neighborPosition);
		sum += kernel(dist);
	});

	return sum;
}

Vector3f SphSystemData3::interpolate(const Vector3f& origin, const VectorArray& values) const
{
	Vector3f sum;
	auto& d = densities();
	SphStdkernel3 kernel(radius());
	neighborSearcher()->forEachNearbyPoint(origin, radius(),
		[&](size_t i, const Vector3f& neighborPosition) {
		float dist = origin.Distance(neighborPosition);
		float weight = mass() / d[i] * kernel(dist);
		sum += weight * values[i];
	});

	return sum;
}

/* Solver begin */

ParticleSystemSolver3::ParticleSystemSolver3()
{
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
	auto forces = const_cast<VectorArray&>(_particleSystemData->forces());
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

	auto positions = const_cast<VectorArray&>(_particleSystemData->positions());
	auto velocities = const_cast<VectorArray&>(_particleSystemData->velocities());
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
	auto forces = const_cast<VectorArray&>(_particleSystemData->forces());
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