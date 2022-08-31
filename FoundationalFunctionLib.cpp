/* Foundations function of particles */
#include "FoundationalFunctionLib.h"

void Collider3::resolveCollision(const Vector3f& currentPosition, const Vector3f& currentVelocity, float radius, float restitutionCoefficent, Vector3f& newPosition, Vector3f& newVelocity)
{

}

PointHashGridSearcher3::PointHashGridSearcher3(const Vector3f& resolution, float gridSpacing) : _gridSpacing(gridSpacing), _resolution(resolution)
{

}

PointHashGridSearcher3::PointHashGridSearcher3(size_t resolutionX, size_t resolutionY, size_t resolutionZ, float gridSpacing) : _gridSpacing(gridSpacing), _resolution(resolutionX, resolutionY, resolutionZ)
{

}

void PointHashGridSearcher3::build(const VectorArray& points)
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
	if (_buckets.empty())
	{
		return;
	}

	size_t nearbyKeys[8];
	getNearbyKeys(origin, nearbyKeys);

	const float queryRadiusSquared = radius * radius;
	for (size_t i = 0; i < 8; i++)
	{
		const auto& bucket = _buckets[nearbyKeys[i]];
		size_t numberOfPointsInBucket = bucket.size();

		for (size_t j = 0; j < numberOfPointsInBucket; j++)
		{
			size_t pointIndex = bucket[j];
			float rSquared = _points[pointIndex].DistanceSq(origin);
			if (rSquared <= queryRadiusSquared)
			{
				callback(pointIndex, _points[pointIndex]);
			}
		}
	}
}

size_t PointHashGridSearcher3::getHashKeyFromPosition(const Vector3f& position)
{
	Vector3f bucketIndex = getBucketIndex(position);
	return getHashKeyFromBucketIndex(bucketIndex);
}

Vector3f PointHashGridSearcher3::getBucketIndex(const Vector3f& position) const
{
	Vector3f index;
	index.x = static_cast<size_t>(floor(position.x / _gridSpacing));
	index.y = static_cast<size_t>(floor(position.y / _gridSpacing));
	index.z = static_cast<size_t>(floor(position.z / _gridSpacing));
	return index;
}

size_t PointHashGridSearcher3::getHashKeyFromBucketIndex(const Vector3f& bucketIndex) const
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


void PointHashGridSearcher3::getNearbyKeys(const Vector3f& origin, size_t* nearbyKeys) const
{
	Vector3f originBucketIndex = getBucketIndex(origin);
	Vector3f nearbyBucketIndeics[8];

	for (size_t i = 0; i < 8; i++)
	{
		nearbyBucketIndeics[i] = originBucketIndex;
	}

	if ((originBucketIndex.x + 0.5f) * _gridSpacing <= origin.x)
	{
		nearbyBucketIndeics[4].x += 1;
		nearbyBucketIndeics[5].x += 1;
		nearbyBucketIndeics[6].x += 1;
		nearbyBucketIndeics[7].x += 1;
	}
	else
	{
		nearbyBucketIndeics[4].x -= 1;
		nearbyBucketIndeics[5].x -= 1;
		nearbyBucketIndeics[6].x -= 1;
		nearbyBucketIndeics[7].x -= 1;
	}

	if ((originBucketIndex.y + 0.5f) * _gridSpacing <= origin.y)
	{
		nearbyBucketIndeics[2].y += 1;
		nearbyBucketIndeics[3].y += 1;
		nearbyBucketIndeics[6].y += 1;
		nearbyBucketIndeics[7].y += 1;
	}
	else
	{
		nearbyBucketIndeics[2].y -= 1;
		nearbyBucketIndeics[3].y -= 1;
		nearbyBucketIndeics[6].y -= 1;
		nearbyBucketIndeics[7].y -= 1;
	}

	if ((originBucketIndex.z + 0.5f) * _gridSpacing <= origin.z)
	{
		nearbyBucketIndeics[1].z += 1;
		nearbyBucketIndeics[3].z += 1;
		nearbyBucketIndeics[5].z += 1;
		nearbyBucketIndeics[7].z += 1;
	}
	else
	{
		nearbyBucketIndeics[1].z -= 1;
		nearbyBucketIndeics[3].z -= 1;
		nearbyBucketIndeics[5].z -= 1;
		nearbyBucketIndeics[7].z -= 1;
	}

	for (size_t i = 0; i < 8; i++)
	{
		nearbyKeys[i] = getHashKeyFromBucketIndex(nearbyBucketIndeics[i]);
	}
}

float SphStdkernel3::operator()(float distance) const
{
	float distanceSquared = powf(distance, 2);
	if (distanceSquared > h2)
	{
		return 0.f;
	}
	else
	{
		float x = 1 - distanceSquared / h2;
		return 315.f / (64 * kPiD * h3) * powf(x, 3);
	}
	return 0.f;
}