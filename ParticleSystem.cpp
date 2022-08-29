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

/* Foundations function of particles */

PointHashGridSearcher3::PointHashGridSearcher3(const Vector3f& resolution, float gridSpacing)
{

}

PointHashGridSearcher3::PointHashGridSearcher3(size_t resolutionX, size_t resolutionY, size_t resolutionZ, float gridSpacing)
{

}

void PointHashGridSearcher3::build(const ParticleSystemData3::VectorArray& points)
{
	_buckets.clear();
	_points.clear();

	if (points.size() == 0)
	{
		return;
	}

	if (_resolution.x <= 0 || _resolution.y <= 0 || _resolution.z <= 0)
	{
		return;
	}

	_buckets.resize(_resolution.x * _resolution.y * _resolution.z);
	_points.resize(points.size());

	for (size_t i = 0; i < points.size(); i++)
	{
		_points[i] = points[i];
		size_t key = getHashKeyFromPosition(points[i]);
		_buckets[key].push_back(i);
	}
}

void PointHashGridSearcher3::forEachNearbyPoint(const Vector3f& origin, float radius, const ForEachNeighborPointFunc& callback) const
{
	_resolution.x;
}

size_t PointHashGridSearcher3::getHashKeyFromPosition(const Vector3f& position)
{
	Vector3f bucketIndex = getBucketIndex(position);
	return getHashKeyFromBucketIndex(bucketIndex);
}

Vector3f PointHashGridSearcher3::getBucketIndex(const Vector3f& position)
{
	Vector3f index;
	index.x = static_cast<size_t>(std::floor(position.x / _gridSpacing));
	index.y = static_cast<size_t>(std::floor(position.y / _gridSpacing));
	index.z = static_cast<size_t>(std::floor(position.z / _gridSpacing));
	return index;
}

size_t PointHashGridSearcher3::getHashKeyFromBucketIndex(const Vector3f& bucketIndex)
{
	Vector3f wrappedIndex = bucketIndex;
	wrappedIndex.x = static_cast<size_t>(bucketIndex.x) % static_cast<size_t>(_resolution.x);
	wrappedIndex.y = static_cast<size_t>(bucketIndex.y) % static_cast<size_t>(_resolution.y);
	wrappedIndex.z = static_cast<size_t>(bucketIndex.z) % static_cast<size_t>(_resolution.z);

	if (wrappedIndex.x < 0)
	{
		wrappedIndex.x += static_cast<size_t>(_resolution.x);
	}
	if (wrappedIndex.y < 0)
	{
		wrappedIndex.y += static_cast<size_t>(_resolution.y);
	}
	if (wrappedIndex.z < 0)
	{
		wrappedIndex.z += static_cast<size_t>(_resolution.z);
	}

	return static_cast<size_t>((wrappedIndex.z * _resolution.y + wrappedIndex.y) * _resolution.x + wrappedIndex.x);
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