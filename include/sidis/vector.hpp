#ifndef __SIDIS_VECTOR_HPP__
#define __SIDIS_VECTOR_HPP__

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

/**
 * A spatial 3-vector.
 */
struct Vec3 {
	static Vec3 const ZERO;
	static Vec3 const X;
	static Vec3 const Y;
	static Vec3 const Z;

	Real x;
	Real y;
	Real z;

	Vec3() : x(0.), y(0.), z(0.) { }
	Vec3(Real x, Real y, Real z) : x(x), y(y), z(z) { }

	// Arithmentic operations.
	Vec3& operator+=(Vec3 const& rhs) {
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	Vec3& operator-=(Vec3 const& rhs) {
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
	Vec3& operator*=(Real const& rhs) {
		x *= rhs;
		y *= rhs;
		z *= rhs;
		return *this;
	}
	Vec3& operator/=(Real const& rhs) {
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}
	Vec3 operator-() const {
		return Vec3(-x, -y, -z);
	}

	// Magnitude operations.
	Real norm_sq() const;
	Real norm() const;
	Vec3 unit() const;

	// Projection operations.
	Vec3 par(Vec3 const& v) const;
	Vec3 perp(Vec3 const& v) const;
};

inline Vec3 operator+(Vec3 lhs, Vec3 const& rhs) {
	lhs += rhs;
	return lhs;
}
inline Vec3 operator-(Vec3 lhs, Vec3 const& rhs) {
	lhs -= rhs;
	return lhs;
}
inline Vec3 operator*(Vec3 lhs, Real const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline Vec3 operator*(Real const& lhs, Vec3 rhs) {
	rhs *= lhs;
	return rhs;
}
inline Vec3 operator/(Vec3 lhs, Real const& rhs) {
	lhs /= rhs;
	return lhs;
}

inline Real dot(Vec3 const& lhs, Vec3 const& rhs) {
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
inline Vec3 cross(Vec3 const& lhs, Vec3 const& rhs) {
	return Vec3(
		lhs.y * rhs.z - lhs.z * rhs.y,
		lhs.z * rhs.x - lhs.x * rhs.z,
		lhs.x * rhs.y - lhs.y * rhs.x);
}
Real angle_between(Vec3 const& vec_a, Vec3 const& vec_b);

/**
 * A 4-vector for representing 4-momenta of particles.
 */
struct Vec4 {
	static Vec4 const ZERO;
	static Vec4 const T;
	static Vec4 const X;
	static Vec4 const Y;
	static Vec4 const Z;

	Real t;
	Vec3 r;

	Vec4() : t(0.), r(0., 0., 0.) { }
	Vec4(Real t, Real x, Real y, Real z) : t(t), r(x, y, z) { }
	Vec4(Real t, Vec3 r) : t(t), r(r) { }

	// Arithmentic operations.
	Vec4& operator+=(Vec4 const& rhs) {
		t += rhs.t;
		r += rhs.r;
		return *this;
	}
	Vec4& operator-=(Vec4 const& rhs) {
		t -= rhs.t;
		r -= rhs.r;
		return *this;
	}
	Vec4& operator*=(Real const& rhs) {
		t *= rhs;
		r *= rhs;
		return *this;
	}
	Vec4& operator/=(Real const& rhs) {
		t /= rhs;
		r /= rhs;
		return *this;
	}
	Vec4 operator-() const {
		return Vec4(-t, -r);
	}

	// Magnitude operations.
	Real norm_sq() const;
	Real norm() const;
	Vec4 unit() const;

	// Projection operations.
	Vec4 par(Vec4 const& v) const;
	Vec4 perp(Vec4 const& v) const;
};

inline Vec4 operator+(Vec4 lhs, Vec4 const& rhs) {
	lhs += rhs;
	return lhs;
}
inline Vec4 operator-(Vec4 lhs, Vec4 const& rhs) {
	lhs -= rhs;
	return lhs;
}
inline Vec4 operator*(Vec4 lhs, Real const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline Vec4 operator*(Real const& lhs, Vec4 rhs) {
	rhs *= lhs;
	return rhs;
}
inline Vec4 operator/(Vec4 lhs, Real const& rhs) {
	lhs /= rhs;
	return lhs;
}

inline Real dot(Vec4 const& lhs, Vec4 const& rhs) {
	return lhs.t * rhs.t - dot(lhs.r, rhs.r);
}
Vec4 cross(Vec4 const& a, Vec4 const& b, Vec4 const& c);

}
}

#endif

