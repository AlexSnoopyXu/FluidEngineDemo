#include "ParticleSystem.h"

ParticleSystemData3::ParticleSystemData3()
{

}

void ParticleSystemData3::resize(size_t newNumberOfParticles)
{
	_positions.resize(newNumberOfParticles);
	_velocities.resize(newNumberOfParticles);
	_forces.resize(newNumberOfParticles);
}

size_t ParticleSystemData3::numOfParticles() const
{
	return _positions.size();
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

const std::vector<vector<size_t>>& ParticleSystemData3::neighborLists() const
{
	return _neighborLists;
}

void ParticleSystemData3::addParticle(const Vector3f& newPosition, const Vector3f& newVelocity, const Vector3f& newForce)
{
	_positions.push_back(newPosition);
	_velocities.push_back(newVelocity);
	_forces.push_back(newForce);
}

void ParticleSystemData3::addParticles(const VectorArray& newPositions, const VectorArray& newVelocities, const VectorArray& newForce)
{
	if (newPositions.size() != newVelocities.size() || newVelocities.size() != newForce.size())
	{
		return;
	}
	for (size_t i = 0; i < newPositions.size(); i++)
	{
		_positions.push_back(newPositions[i]);
		_velocities.push_back(newVelocities[i]);
		_forces.push_back(newForce[i]);
	}
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

float SphSystemData3::sumOfKernelNearby(const Vector3f& position, SPHKernelType kernelType) const
{
	float sum = 0.f;
	SphKernelBase3 kernel;
	if (kernelType == SPHKernelType::Std)
	{
		kernel = SphStdKernel3(radius());
	}
	else
	{
		kernel = SphSpikyKernel3(radius());
	}
	
	neighborSearcher()->forEachNearbyPoint(position, radius(),
		[&](size_t, const Vector3f& neighborPosition) {
		float dist = position.Distance(neighborPosition);
		sum += kernel(dist);
	});

	return sum;
}

Vector3f SphSystemData3::interpolate(const Vector3f& origin, const VectorArray& values, SPHKernelType kernelType) const
{
	Vector3f sum;
	auto& d = densities();
	SphKernelBase3 kernel;
	if (kernelType == SPHKernelType::Std)
	{
		kernel = SphStdKernel3(radius());
	}
	else
	{
		kernel = SphSpikyKernel3(radius());
	}
	neighborSearcher()->forEachNearbyPoint(origin, radius(),
		[&](size_t i, const Vector3f& neighborPosition) {
		float dist = origin.Distance(neighborPosition);
		float weight = mass() / d[i] * kernel(dist);
		sum += weight * values[i];
	});

	return sum;
}

Vector3f SphSystemData3::gradientAt(size_t i, const FloatArray& values, SPHKernelType kernelType) const
{
	Vector3f sum;

	auto p = positions();
	auto d = densities();
	const auto& neighbors = neighborLists()[i];
	const Vector3f& origin = p[i];
	SphKernelBase3 kernel;
	if (kernelType == SPHKernelType::Std)
	{
		kernel = SphStdKernel3(radius());
	}
	else
	{
		kernel = SphSpikyKernel3(radius());
	}

	for (size_t j : neighbors)
	{
		const Vector3f& neighborPosition = p[j];
		float dist = origin.Distance(neighborPosition);
		if (dist > 0)
		{
			const Vector3f& dir = (neighborPosition - origin).Normalized();
			sum += d[i] * mass() * (values[i] / powf(d[i], 2) + values[j] / powf(d[j], 2)) * kernel.gradient(dist, dir);
		}
	}

	return sum;
}

float SphSystemData3::laplaciantAt(size_t i, const FloatArray& values, SPHKernelType kernelType) const
{
	float sum = 0.f;

	auto p = positions();
	auto d = densities();
	const auto& neighbors = neighborLists()[i];
	const Vector3f& origin = p[i];
	SphKernelBase3 kernel;
	if (kernelType == SPHKernelType::Std)
	{
		kernel = SphStdKernel3(radius());
	}
	else
	{
		kernel = SphSpikyKernel3(radius());
	}	

	for (size_t j : neighbors)
	{
		const Vector3f& neighborPosition = p[j];
		float dist = origin.Distance(neighborPosition);
		sum += mass() * (values[j] - values[j]) / d[j] * kernel.secondDerivative(dist);
	}

	return sum;
}

/* Solver begin */

ParticleSystemSolver3::ParticleSystemSolver3()
{
	_particleSystemData = std::make_shared<ParticleSystemData3>();
	_wind = std::make_shared<ConstantVectorField3>(Vector3f(1.0f, 1.0f, 1.0f));
	_collider = std::make_shared<Collider3>();
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
			_collider->resolveCollision(positions[i], velocities[i], radius, _restitutionCoefficient, _newPositions[i], _newVelocities[i]);
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


/*
*	SphSystemSolver
*/


SphSystemSolver3::SphSystemSolver3() {
	_sphSystemData3 = std::make_shared<SphSystemData3>();
}

void SphSystemSolver3::accumulateForces(double timeStepInSeconds)
{
	accumulateNonPressureForces(timeStepInSeconds);
	accumulatePressureForces(timeStepInSeconds);
}

void SphSystemSolver3::onBeginAdvanceTimeStep()
{
	auto particles = sphSystemData3();
	particles->buildNearbySearcher(_kMaxSearchRadius);
	particles->buildNearbyList(_kMaxSearchRadius);
	particles->updateDensities();
}

void SphSystemSolver3::onEndAdvanceTimeStep()
{
	computePseudeViscosity();
}

void SphSystemSolver3::accumulateNonPressureForces(float timeStepInSeconds)
{
	ParticleSystemSolver3::accumulateForces(timeStepInSeconds);
	accumulateViscosityForece();
}

void SphSystemSolver3::accumulatePressureForces(float timeStepInSeconds)
{
	auto particles = sphSystemData3();
	auto x = particles->positions();
	auto d = particles->densities();
	auto f = particles->forces();
	auto p = particles->pressures();

	computePressure();
	
	accumulatePressureForces(x, d, p, _pressureForces);
}

void SphSystemSolver3::accumulatePressureForces(VectorArray& positions, FloatArray& densities, FloatArray& pressures, VectorArray& pressureForces)
{
	auto particles = sphSystemData3();
	computePressure();

	size_t numberOfParticls = particles->numOfParticles();
	const float massSquared = Square(particles->mass());
	const SphSpikyKernel3 kernel(particles->radius());

	for (size_t i = 0; i < numberOfParticls; i++)
	{
		const auto& neighbors = particles->neighborLists()[i];
		for (size_t j : neighbors)
		{
			double dist = positions[i].Distance(positions[j]);
			if (dist)
			{
				Vector3f dir = (positions[j] - positions[i]).Normalized();
				pressureForces[i] -= massSquared * ((pressures[i] / Square(densities[i])) + (pressures[j] / Square(densities[j]))) * kernel.gradient(dist, dir);
			}
		}
	}
}

void SphSystemSolver3::computePressure()
{
	auto particles = sphSystemData3();
	size_t numberOfParticls = particles->numOfParticles();
	auto p = particles->pressures();
	auto d = particles->densities();

	const float targetDensity = particles->targetDensity();
	const float eosScale = targetDensity * powf(_kSpeedOfSound, 2) / _eosExponent;

	for (size_t i = 0; i < numberOfParticls; i++)
	{
		p[i] = computePressureFromEOS(d[i], targetDensity, eosScale, _eosExponent);
	}
}

void SphSystemSolver3::accumulateViscosityForece()
{
	auto particles = sphSystemData3();
	size_t numberOfParticls = particles->numOfParticles();
	auto x = particles->positions();
	auto d = particles->densities();
	auto f = particles->forces();
	auto v = particles->velocities();

	const float massSquared = Square(particles->mass());
	const SphSpikyKernel3 kernel(particles->radius());

	for (size_t i = 0; i < numberOfParticls; i++)
	{
		const auto& neighbors = particles->neighborLists()[i];
		for (size_t j : neighbors)
		{
			double dist = x[i].Distance(x[j]);
			f[i] += viscosityCoefficient() * massSquared * (v[j] - v[i]) / d[j] * kernel.secondDerivative(dist);
		}
	}
}

void SphSystemSolver3::computePseudeViscosity()
{

}

float SphSystemSolver3::computePressureFromEOS(float density, float targetDensity, float eosScale, float eosExponent, float negativePressureScale)
{
	float p = eosScale / eosExponent * std::powf((density / targetDensity) - 1, eosExponent);
	if (p < 0)
	{
		p *= negativePressureScale;
	}
	return p;
}

void PciSphSystemSolver3::resolveCollision(VectorArray& oriPositions, VectorArray& oriVelocities, VectorArray& newPositions, VectorArray& newVelocity)
{
	auto particles = sphSystemData3();
	auto positions = oriPositions;
	auto velocities = oriVelocities;

	if (collider() != nullptr)
	{
		size_t n = particles->numOfParticles();
		const float radius = particles->radius();
		for (size_t i = 0; i < n; i++)
		{
			collider()->resolveCollision(positions[i], velocities[i], radius, restitutionCoefficient(), newPositions[i], newVelocity[i]);
		}
	}
}

void PciSphSystemSolver3::accumulatePressureForces(float timeStepInSeconds)
{
	auto particles = sphSystemData3();
	size_t numberOfParticls = particles->numOfParticles();
	const float targetDensity = particles->targetDensity();
	const float mass = particles->mass();
	const float delta = computeDelta(timeStepInSeconds);

	auto p = particles->pressures();
	auto x = particles->positions();
	auto v = particles->velocities();
	auto f = particles->forces();

	FloatArray ds(numberOfParticls, 0.0f);
	SphStdKernel3 kernel(particles->radius());

	_tempPositions.resize(numberOfParticls);
	_tempVelocities.resize(numberOfParticls);
	_pressureForces.resize(numberOfParticls);
	_densityErrors.resize(numberOfParticls);

	for (size_t i = 0; i < numberOfParticls; i++)
	{
		p[i] = 0.f;
		_pressureForces[i] = Vector3f();
		_densityErrors[i] = 0.f;
	}

	for (size_t i = 0; i < _maxNumberOfIterations; i++)
	{
		// Predict position and velocity
		for (size_t i = 0; i < numberOfParticls; i++)
		{
			_tempVelocities[i] = v[i] + timeStepInSeconds / mass * (f[i] + _pressureForces[i]);
			_tempPositions[i] = x[i] + timeStepInSeconds * _tempVelocities[i];
		}

		// Resolve collision
		resolveCollision(_tempPositions, _tempVelocities, _tempPositions, _tempVelocities);;

		// Compute pressure from density error
		for (size_t i = 0; i < numberOfParticls; i++)
		{
			float weightSum = 0.f;
			const auto& neighbors = particles->neighborLists()[i];
			for (size_t j : neighbors)
			{
				float dist = _tempPositions[j].Distance(_tempPositions[i]);
				weightSum += kernel(dist);
			}
			weightSum += kernel(0);

			float density = weightSum * mass;
			float densityError = density - targetDensity;
			float pressure = delta * densityError;

			if (pressure < 0)
			{
				pressure *= negativePressureScale();
				densityError *= negativePressureScale();
			}

			p[i] += pressure;
			ds[i] = density;
			_densityErrors[i] = densityError;
		}

		// Compute pressure gradient force
		SphSystemSolver3::accumulatePressureForces(x, ds, p, _pressureForces);

		// Compute max density error
		float maxDensityError = 0.f;
		for (size_t i = 0; i < numberOfParticls; i++)
		{
			maxDensityError = abs(max(maxDensityError, _densityErrors[i]));
		}

		float densityErrorRatio = maxDensityError / targetDensity;
		if (abs(densityErrorRatio) < _maxDensityErrorRatio)
		{
			break;
		}
	}

	// Accumulate pressure force
	for (size_t i = 0; i < numberOfParticls; i++)
	{
		f[i] += _pressureForces[i];
	}
}