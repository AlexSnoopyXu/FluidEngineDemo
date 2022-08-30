#pragma once
#include "math.h"
#include <vector>
#include "MathGeoLib-1.5/src/Math/float3.h"

typedef float3 Vector3f;
typedef std::vector<Vector3f> VectorArray;

class MathFoundationLib
{
};

struct Vector3D
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	Vector3D() = default;
	Vector3D(double ix, double iy, double iz) : x(ix), y(iy), z(iz) {}

	double size()
	{
		return sqrt(sizeSquare());
	}

	double sizeSquare()
	{
		return x * x + y * y + z * z;
	}

	Vector3D normalized()
	{
		const double s = size();
		return Vector3D(x / s, y / x, z / s);
	}

	// dot
	double operator|(const Vector3D& v)
	{
		return x * v.x + y * v.y;
	}

	// cross
	Vector3D operator*(const Vector3D& b)
	{
		return Vector3D(y * b.z - b.y * z, z * b.x - x * b.z, x * b.y - b.x * y);
	}

	Vector3D operator+(const Vector3D& v)
	{
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}

	Vector3D operator-(const Vector3D& v)
	{
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}
};

class Field3 {
public:
	Field3() {}
	virtual ~Field3() {}
};

class ScalarField3 : public Field3 {
public:
	ScalarField3() {}
	virtual ~ScalarField3() {}

	virtual float sample(const Vector3f& x) const = 0;
};

class VectorField3 : public Field3 {
public:
	VectorField3() {}
	virtual ~VectorField3() {}

	virtual Vector3f sample(const Vector3f& x) const = 0;
};

class ConstantVectorField3 : public VectorField3 {
public:
	ConstantVectorField3() {}
	ConstantVectorField3(const Vector3f& value) : _value(value) {}
	virtual ~ConstantVectorField3() {}

	virtual Vector3f sample(const Vector3f& x) const override
	{
		return _value;
	}

private:
	const Vector3f _value = Vector3f(1.f, 1.f, 1.f);
};

class Collider3 {
public:
	Collider3() {}
	virtual ~Collider3() {}

	void resolveCollision(const Vector3f& currentPosition, const Vector3f& currentVelocity, float radius, float restitutionCoefficent, Vector3f& newPosition, Vector3f& newVelocity);
};