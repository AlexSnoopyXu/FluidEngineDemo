#pragma once
#include "MathFoundationLib.h"
#include <functional>
#include <vector>

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

