#pragma once
#include "AnimationBase.h"
#include <functional>
#include "FoundationalFunctionLib.h"

#define Square(x) MathFoundationLib::square(x)

struct ParticleSystemData3 {
public:

	ParticleSystemData3();
	ParticleSystemData3(float mass, float radius) : _mass(mass), _radius(radius) {}
	virtual ~ParticleSystemData3() {}

	void resize(size_t newNumberOfParticles);
	size_t numOfParticles() const;

	const VectorArray& const positions() const;
	const VectorArray& const velocities() const;
	const VectorArray& const forces() const;
	const float mass() const;
	const float radius() const;
	const std::vector<vector<size_t>>& neighborLists() const;

	void addParticle(const Vector3f& newPosition, const Vector3f& newVelocity = Vector3f(), const Vector3f& newForce = Vector3f());
	void addParticles(const VectorArray& newPositions, const VectorArray& newVelocities = VectorArray(), const VectorArray& newForce = VectorArray());
	void setMass(float mass);

	void buildNearbySearcher(float maxSearchRadius);
	void buildNearbyList(float maxSearchRadius);

protected:
	std::shared_ptr<PointHashGridSearcher3> neighborSearcher() const { return _pointNeighborSearcher; }

private:
	VectorArray _positions;
	VectorArray _velocities;
	VectorArray _forces;
	float _mass;
	float _radius;

	const size_t kDefaultHashGridResolution = 100;

	std::shared_ptr<PointHashGridSearcher3> _pointNeighborSearcher;
	std::vector<vector<size_t>> _neighborLists;
};

struct SphSystemData3 : public ParticleSystemData3
{
public:
	SphSystemData3() {}
	virtual ~SphSystemData3() {}

	const FloatArray& densities() const { return _densities; }
	FloatArray& densities() { return _densities; }

	const float targetDensity() const { return 10.f; }

	const FloatArray& pressures() const { return _pressures; }
	FloatArray& pressures() { return _pressures; }

	void updateDensities();
	float sumOfKernelNearby(const Vector3f& position, SPHKernelType kernelType = SPHKernelType::Std) const;

	Vector3f interpolate(const Vector3f& origin, const VectorArray& values, SPHKernelType kernelType = SPHKernelType::Std) const;

	Vector3f gradientAt(size_t i, const FloatArray& values, SPHKernelType kernelType = SPHKernelType::Std) const;
	float laplaciantAt(size_t i, const FloatArray& values, SPHKernelType kernelType = SPHKernelType::Std) const;
private:
	FloatArray _densities;
	FloatArray _pressures;
};

class ParticleSystemSolver3 : public PhysicsAnimation 
{
public:
	ParticleSystemSolver3();
	virtual ~ParticleSystemSolver3() {}

protected:
	virtual void onAdvanceTimStep(double timeIntervalInSeconds) override;
	virtual void accumulateForces(double timeStepInSeconds);
	void resolveCollision();
	float viscosityCoefficient() const { return _viscosityCoefficient; }

	virtual void onBeginAdvanceTimeStep() = 0;
	virtual void onEndAdvanceTimeStep() = 0;
	float negativePressureScale() { return -1.f; }

	shared_ptr<Collider3> collider() { return _collider; }
	float restitutionCoefficient() { return _restitutionCoefficient; }

private:
	shared_ptr<ParticleSystemData3> _particleSystemData;
	shared_ptr<ConstantVectorField3> _wind;

	float _dragCoefficient = 1e-4;
	Vector3f _gravity = Vector3f(0.f, 9.8f, 0.f);
	float _restitutionCoefficient = 5.f;
	float _viscosityCoefficient = 5.f;

	VectorArray _newPositions;
	VectorArray _newVelocities;
	shared_ptr<Collider3> _collider;

private:
	void beginAdvanceTimeStep();
	void endAdvanceTimeStep();
	void timeIntegration(double timeIntervalInSeconds);
	void accumulateExternalForces();
};

class SphSystemSolver3 : public ParticleSystemSolver3
{
public:
	SphSystemSolver3();
	virtual ~SphSystemSolver3(){}

	shared_ptr<SphSystemData3> sphSystemData3() const { return _sphSystemData3; }

protected:
	virtual void accumulateForces(double timeStepInSeconds) override;
	virtual void onBeginAdvanceTimeStep() override;
	virtual void onEndAdvanceTimeStep() override;

	virtual void accumulateNonPressureForces(float timeStepInSeconds);
	virtual void accumulatePressureForces(float timeStepInSeconds);
	void accumulatePressureForces(VectorArray& positions, FloatArray& densities, FloatArray& pressures, VectorArray& pressureForces);

	void computePressure();
	void accumulateViscosityForece();
	void computePseudeViscosity();

	float computePressureFromEOS(float density, float targetDensity, float eosScale, float eosExponent, float negativePressureScale = 0);

protected:
	VectorArray _pressureForces;

private:
	shared_ptr<SphSystemData3> _sphSystemData3;

	const float _kMaxSearchRadius = 10.f;
	const float _kSpeedOfSound = 340.f;
	const float _eosExponent = 2.f;
};

class PciSphSystemSolver3 : public SphSystemSolver3 
{
	PciSphSystemSolver3() {}
	virtual ~PciSphSystemSolver3(){}

protected:
	virtual void accumulatePressureForces(float timeStepInSeconds) override;

private:
	float computeDelta(float timeStepInSeconds) { return 23.3f; }
	void resolveCollision(VectorArray& oriPositions, VectorArray& oriVelocities, VectorArray& newPositions, VectorArray& newVelocity);
private:
	float _maxDensityErrorRatio = 0.01f;
	size_t _maxNumberOfIterations = 5;

	VectorArray _tempPositions;
	VectorArray _tempVelocities;	
	FloatArray _densityErrors;
};