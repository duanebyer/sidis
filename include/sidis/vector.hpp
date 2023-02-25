#ifndef SIDIS_VECTOR_HPP
#define SIDIS_VECTOR_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

/**
 * \defgroup LinAlgGroup Linear algebra
 *
 * Types and functions for linear algebra in 3 and 4 dimensions.
 */
/// \{

/**
 * A spatial 3-vector.
 */
struct Vec3 {
	Real x;
	Real y;
	Real z;

	Vec3() : x(0.), y(0.), z(0.) { }
	Vec3(Real x, Real y, Real z) : x(x), y(y), z(z) { }

	// Equality operations.
	bool operator==(Vec3 const& rhs) const {
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}
	bool operator!=(Vec3 const& rhs) const {
		return !(*this == rhs);
	}

	// Arithmetic operations.
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

Vec3 const VEC3_ZERO(0., 0., 0.);
Vec3 const VEC3_X(1., 0., 0.);
Vec3 const VEC3_Y(0., 1., 0.);
Vec3 const VEC3_Z(0., 0., 1.);

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
Real angle_between(Vec3 const& v_a, Vec3 const& v_b);

/**
 * A 4-vector for representing 4-momenta of particles.
 */
struct Vec4 {
	Real t;
	Real x;
	Real y;
	Real z;

	Vec4() : t(0.), x(0.), y(0.), z(0.) { }
	Vec4(Real t, Real x, Real y, Real z) : t(t), x(x), y(y), z(z) { }
	Vec4(Real t, Vec3 r) : t(t), x(r.x), y(r.y), z(r.z) { }

	static Vec4 from_length_and_r(Real m, Vec3 const& p);
	static Vec4 from_length_and_t(Real m, Real t, Vec3 const& dir);

	Vec3 r() const {
		return Vec3(x, y, z);
	}

	// Equality operations.
	bool operator==(Vec4 const& rhs) const {
		return t == rhs.t && x == rhs.x && y == rhs.y && z == rhs.z;
	}
	bool operator!=(Vec4 const& rhs) const {
		return !(*this == rhs);
	}

	// Arithmetic operations.
	Vec4& operator+=(Vec4 const& rhs) {
		t += rhs.t;
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	Vec4& operator-=(Vec4 const& rhs) {
		t -= rhs.t;
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
	Vec4& operator*=(Real const& rhs) {
		t *= rhs;
		x *= rhs;
		y *= rhs;
		z *= rhs;
		return *this;
	}
	Vec4& operator/=(Real const& rhs) {
		t /= rhs;
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}
	Vec4 operator-() const {
		return Vec4(-t, -x, -y, -z);
	}

	// Magnitude operations.
	Real norm_sq() const;
	Real norm() const;
	Vec4 unit() const;
	int sign() const;

	// Projection operations.
	Vec4 par(Vec4 const& v) const;
	Vec4 perp(Vec4 const& v) const;
};

Vec4 const VEC4_ZERO(0., 0., 0., 0.);
Vec4 const VEC4_T(1., 0., 0., 0.);
Vec4 const VEC4_X(0., 1., 0., 0.);
Vec4 const VEC4_Y(0., 0., 1., 0.);
Vec4 const VEC4_Z(0., 0., 0., 1.);

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
	return lhs.t * rhs.t - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z;
}
Vec4 cross(Vec4 const& a, Vec4 const& b, Vec4 const& c);

/// \}

}
}

#endif

