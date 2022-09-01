#pragma once
#include "MathFoundationLib.h"
#include <functional>
#include <vector>

enum class SPHKernelType {
	None = 0,

	Std = 1,
	Spiky = 2,

	Max
};

class Collider3 {
public:
	Collider3() {}
	virtual ~Collider3() {}

	void resolveCollision(const Vector3f& currentPosition, const Vector3f& currentVelocity, float radius, float restitutionCoefficent, Vector3f& newPosition, Vector3f& newVelocity);
};

class PointNeighborSearcher3
{
public:
	typedef std::function<void(size_t, const Vector3f&)> ForEachNeighborPointFunc;

	PointNeighborSearcher3() {}
	virtual ~PointNeighborSearcher3() {}

	virtual void build(const VectorArray& points) = 0;

	virtual void forEachNearbyPoint(const Vector3f& origin, float radius, const ForEachNeighborPointFunc& callback) const = 0;
};

class PointHashGridSearcher3 final : public PointNeighborSearcher3
{
public:
	PointHashGridSearcher3(const Vector3f& resolution, float gridSpacing);
	PointHashGridSearcher3(size_t resolutionX, size_t resolutionY, size_t resolutionZ, float gridSpacing);

	virtual void build(const VectorArray& points) override;
	virtual void forEachNearbyPoint(const Vector3f& origin, float radius, const ForEachNeighborPointFunc& callback) const override;

private:
	float _gridSpacing = 0;
	Vector3f _resolution;
	VectorArray _points;
	std::vector<std::vector<size_t>> _buckets;

private:
	size_t getHashKeyFromPosition(const Vector3f& position);
	Vector3f getBucketIndex(const Vector3f& position) const;
	size_t getHashKeyFromBucketIndex(const Vector3f& bucketIndex) const;
	void getNearbyKeys(const Vector3f& origin, size_t* nearbyKeys) const;
};

struct SphKernelBase3
{
	SphKernelBase3() {}
	virtual ~SphKernelBase3(){}
	virtual float operator()(float distance) const { return 0.f; }

	virtual float firstDerivative(float distance) const { return 0.f; }
	virtual Vector3f gradient(float distance, const Vector3f& direction) const { return Vector3f(); }

	virtual float secondDerivative(float distance) const { return 0.f; }
};

struct SphStdKernel3 : public SphKernelBase3
{
	float h, h2, h3, h5;

	SphStdKernel3() : h(0), h2(0), h3(0), h5(0) {}
	explicit SphStdKernel3(float kernelRadius) : h(kernelRadius), h2(h * h), h3(h2 * h), h5(h3 * h2) {}
	SphStdKernel3(const SphStdKernel3& other) : h(other.h), h2(other.h2), h3(other.h3), h5(other.h5) {}
	virtual float operator()(float distance) const override;

	/*
	* Differential operators
	*/
	virtual float firstDerivative(float distance) const override;
	virtual Vector3f gradient(float distance, const Vector3f& direction) const override;

	virtual float secondDerivative(float distance) const override;
};

struct SphSpikyKernel3 : public SphKernelBase3
{
	float h, h2, h3, h4, h5;

	SphSpikyKernel3() : h(0), h2(0), h3(0), h4(0), h5(0) {}
	explicit SphSpikyKernel3(float kernelRadius) : h(kernelRadius), h2(h* h), h3(h2* h), h4(h3 * h), h5(h3* h2) {}
	SphSpikyKernel3(const SphSpikyKernel3& other) : h(other.h), h2(other.h2), h3(other.h3), h4(other.h4), h5(other.h5) {}
	virtual float operator()(float distance) const override;

	/*
	* Differential operators
	*/
	virtual float firstDerivative(float distance) const override;
	virtual Vector3f gradient(float distance, const Vector3f& direction) const override;

	virtual float secondDerivative(float distance) const override;
};