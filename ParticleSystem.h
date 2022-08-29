#pragma once
#include "AnimationBase.h"
#include <functional>

struct ParticleSystemData3 {
public:
	typedef std::vector<Vector3f> VectorArray;

	ParticleSystemData3();
	ParticleSystemData3(float mass, float radius) : _mass(mass), _radius(radius) {}
	virtual ~ParticleSystemData3() {}

	void resize(size_t newNumberOfParticles);
	size_t numOfParticles() const;

	const Vector3f* const positions() const;
	const Vector3f* const velocities() const;
	const Vector3f* const forces() const;
	const float mass() const;
	const float radius() const;

	void addParticle(const Vector3f& newPosition, const Vector3f& newVelocity = Vector3f(), const Vector3f& newForce = Vector3f());
	void addParticles(const VectorArray& newPositions, const VectorArray& newVelocities = VectorArray(), const VectorArray& newForce = VectorArray());
	void setMass(float mass);
private:
	VectorArray _positions;
	VectorArray _velocities;
	VectorArray _forces;
	float _mass;
	float _radius;
};

class PointNeighborSearcher3
{
public:
	typedef std::function<void(size_t, const Vector3f&)> ForEachNeighborPointFunc;

	PointNeighborSearcher3();
	virtual ~PointNeighborSearcher3(){}

	virtual void build(const ParticleSystemData3::VectorArray& points) = 0;

	virtual void forEachNearbyPoint(const Vector3f& origin, float radius, const ForEachNeighborPointFunc& callback) const = 0;
};

class PointHashGridSearcher3 final : public PointNeighborSearcher3
{
public:
	PointHashGridSearcher3(const Vector3f& resolution, float gridSpacing);
	PointHashGridSearcher3(size_t resolutionX, size_t resolutionY, size_t resolutionZ, float gridSpacing);

	virtual void build(const ParticleSystemData3::VectorArray& points) override;
	virtual void forEachNearbyPoint(const Vector3f& origin, float radius, const ForEachNeighborPointFunc& callback) const override;

private:
	float _gridSpacing;
	Vector3f _resolution;
	ParticleSystemData3::VectorArray _points;
	std::vector<std::vector<size_t>> _buckets;

private:
	size_t getHashKeyFromPosition(const Vector3f& position);
	Vector3f getBucketIndex(const Vector3f& position);
	size_t getHashKeyFromBucketIndex(const Vector3f& bucketIndex);
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

	virtual void onBeginAdvanceTimeStep() = 0;
	virtual void onEndAdvanceTimeStep() = 0;

private:
	shared_ptr<ParticleSystemData3> _particleSystemData;
	shared_ptr<ConstantVectorField3> _wind;

	float _dragCoefficient = 1e-4;
	Vector3f _gravity = Vector3f(0.f, 9.8f, 0.f);
	float _restitutionCoeffcient = 5.f;

	ParticleSystemData3::VectorArray _newPositions;
	ParticleSystemData3::VectorArray _newVelocities;
	shared_ptr<Collider3> _collider;

private:
	void beginAdvanceTimeStep();
	void endAdvanceTimeStep();
	void timeIntegration(double timeIntervalInSeconds);
	void accumulateExternalForces();
};

