#pragma once
#include "math.h"
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