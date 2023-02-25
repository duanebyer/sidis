#ifndef SIDIS_TRANSFORM_HPP
#define SIDIS_TRANSFORM_HPP

#include "sidis/numeric.hpp"
#include "sidis/vector.hpp"

namespace sidis {
namespace math {

struct Transform4;

/**
 * \addtogroup LinAlgGroup
 */
/// \{

/**
 * A linear transformation acting on 3-vectors.
 */
struct Transform3 {
	static Transform3 rotate(Vec3 const& dir, Real angle);
	static Transform3 rotate_to(Vec3 const& dir_old, Vec3 const& dir_new);
	static Transform3 rotate_to(Vec3 const& z_new);
	static Transform3 rotate_basis(Vec3 const& z_axis, Vec3 const& y_up);
	static Transform3 scale(Vec3 const& dir, Real scale);
	static Transform3 project(Vec3 const& dir);
	static Transform3 bivec(Vec3 const& dir);

	Vec3 x;
	Vec3 y;
	Vec3 z;

	Transform3() : x(VEC3_ZERO), y(VEC3_ZERO), z(VEC3_ZERO) { }
	Transform3(Vec3 const& x, Vec3 const& y, Vec3 const& z)
		: x(x), y(y), z(z) { }
	Transform3(
		Real xx, Real xy, Real xz,
		Real yx, Real yy, Real yz,
		Real zx, Real zy, Real zz) :
		x(xx, xy, xz),
		y(yx, yy, yz),
		z(zx, zy, zz) { }
	explicit Transform3(Real s) : Transform3(
		s, 0., 0.,
		0., s, 0.,
		0., 0., s) { };

	Transform3 transpose() const;
	Real trace() const;
	Real det() const;
	Transform3 inv() const;

	// Equality operations.
	bool operator==(Transform3 const& rhs) const {
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}
	bool operator!=(Transform3 const& rhs) const {
		return !(*this == rhs);
	}

	// Arithmetic operations.
	Transform3& operator+=(Transform3 const& rhs) {
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	Transform3& operator-=(Transform3 const& rhs) {
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
	Transform3& operator*=(Real const& rhs) {
		x *= rhs;
		y *= rhs;
		z *= rhs;
		return *this;
	}
	Transform3& operator/=(Real const& rhs) {
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}
	Transform3& operator*=(Transform3 const& rhs) {
		Transform3 rhs_tr = rhs.transpose();
		x = Vec3(dot(x, rhs_tr.x), dot(x, rhs_tr.y), dot(x, rhs_tr.z));
		y = Vec3(dot(y, rhs_tr.x), dot(y, rhs_tr.y), dot(y, rhs_tr.z));
		z = Vec3(dot(z, rhs_tr.x), dot(z, rhs_tr.y), dot(z, rhs_tr.z));
		return *this;
	}
	Transform3 operator-() const {
		return Transform3(-x, -y, -z);
	}

	// Apply coordinate transformation to objects.
	Vec3 transform(Vec3 const& other) const;
	Vec4 transform(Vec4 const& other) const;
	Transform3 transform(Transform3 const& other) const;
	Transform4 transform(Transform4 const& other) const;
};

Transform3 const TRANSFORM3_ZERO(
	0., 0., 0.,
	0., 0., 0.,
	0., 0., 0.);
Transform3 const TRANSFORM3_ID(
	1., 0., 0.,
	0., 1., 0.,
	0., 0., 1.);

inline Transform3 operator+(Transform3 lhs, Transform3 const& rhs) {
	lhs += rhs;
	return lhs;
}
inline Transform3 operator-(Transform3 lhs, Transform3 const& rhs) {
	lhs -= rhs;
	return lhs;
}
inline Transform3 operator*(Transform3 lhs, Real const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline Transform3 operator*(Real const& lhs, Transform3 rhs) {
	rhs *= lhs;
	return rhs;
}
inline Transform3 operator/(Transform3 lhs, Real const& rhs) {
	lhs /= rhs;
	return lhs;
}
inline Transform3 operator*(Transform3 lhs, Transform3 const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline Vec3 operator*(Transform3 lhs, Vec3 rhs) {
	return Vec3(
		dot(lhs.x, rhs),
		dot(lhs.y, rhs),
		dot(lhs.z, rhs));
}

inline Vec3 dot(Transform3 const& lhs, Vec3 const& rhs) {
	return lhs * rhs;
}
inline Vec3 dot(Vec3 const& lhs, Transform3 const& rhs) {
	return rhs.transpose() * lhs;
}
inline Transform3 dot(Transform3 const& lhs, Transform3 const& rhs) {
	return lhs * rhs;
}
inline Transform3 outer(Vec3 const& lhs, Vec3 const& rhs) {
	return Transform3(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}
inline Transform3 cross(Vec3 const& v) {
	return Transform3(
		0., -v.z, v.y,
		v.z, 0., -v.x,
		-v.y, v.x, 0.);
}
inline Transform3 asym_prod(Vec3 const& v1, Vec3 const& v2) {
	return outer(v1, v2) - outer(v2, v1);
}
inline Transform3 sym_prod(Vec3 const& v1, Vec3 const& v2) {
	return outer(v1, v2) + outer(v1, v2);
}

/**
 * A linear transformation acting on 4-vectors.
 */
struct Transform4 {
	static Transform4 rotate(Vec3 const& dir, Real angle);
	static Transform4 rotate_to(Vec3 const& dir_old, Vec3 const& dir_new);
	static Transform4 rotate_to(Vec3 const& z_axis);
	static Transform4 rotate_basis(Vec3 const& z_axis, Vec3 const& y_up);
	static Transform4 scale(Vec3 const& dir, Real scale);
	static Transform4 project(Vec3 const& dir);
	static Transform4 project(Vec4 const& dir);
	static Transform4 boost(Vec3 const& dir, Real rapidity);
	static Transform4 boost_to(Vec4 const& t_new);
	static Transform4 transform_to(Vec4 const& dir_old, Vec4 const& dir_new);
	static Transform4 transform_to(Vec4 const& t_new);
	static Transform4 bivec(Vec3 const& v, Vec3 const& bv);

	Vec4 t;
	Vec4 x;
	Vec4 y;
	Vec4 z;

	Transform4() : t(VEC4_ZERO), x(VEC4_ZERO), y(VEC4_ZERO), z(VEC4_ZERO) { }
	Transform4(Vec4 const& t, Vec4 const& x, Vec4 const& y, Vec4 const& z)
		: t(t), x(x), y(y), z(z) { }
	Transform4(
		Real tt, Real tx, Real ty, Real tz,
		Real xt, Real xx, Real xy, Real xz,
		Real yt, Real yx, Real yy, Real yz,
		Real zt, Real zx, Real zy, Real zz) :
		t(tt, -tx, -ty, -tz),
		x(xt, -xx, -xy, -xz),
		y(yt, -yx, -yy, -yz),
		z(zt, -zx, -zy, -zz) { }
	explicit Transform4(Real s) : Transform4(
		s, 0., 0., 0.,
		0., s, 0., 0.,
		0., 0., s, 0.,
		0., 0., 0., s) { };
	Transform4(Transform3 const& tr) :
		Transform4(
			VEC4_T,
			Vec4(0., -tr.x),
			Vec4(0., -tr.y),
			Vec4(0., -tr.z)) { }

	Transform4 transpose() const;
	Real trace() const;
	Real det() const;
	Transform4 inv() const;

	// Equality operations.
	bool operator==(Transform4 const& rhs) const {
		return t == rhs.t && x == rhs.x && y == rhs.y && z == rhs.z;
	}
	bool operator!=(Transform4 const& rhs) const {
		return !(*this == rhs);
	}

	// Arithmetic operations.
	Transform4& operator+=(Transform4 const& rhs) {
		t += rhs.t;
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	Transform4& operator-=(Transform4 const& rhs) {
		t -= rhs.t;
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
	Transform4& operator*=(Real const& rhs) {
		t *= rhs;
		x *= rhs;
		y *= rhs;
		z *= rhs;
		return *this;
	}
	Transform4& operator/=(Real const& rhs) {
		t /= rhs;
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}
	Transform4& operator*=(Transform4 const& rhs) {
		Transform4 rhs_tr = rhs.transpose();
		t = Vec4(dot(t, rhs_tr.t), dot(t, rhs_tr.x), dot(t, rhs_tr.y), dot(t, rhs_tr.z));
		x = Vec4(dot(x, rhs_tr.t), dot(x, rhs_tr.x), dot(x, rhs_tr.y), dot(x, rhs_tr.z));
		y = Vec4(dot(y, rhs_tr.t), dot(y, rhs_tr.x), dot(y, rhs_tr.y), dot(y, rhs_tr.z));
		z = Vec4(dot(z, rhs_tr.t), dot(z, rhs_tr.x), dot(z, rhs_tr.y), dot(z, rhs_tr.z));
		return *this;
	}
	Transform4 operator-() const {
		return Transform4(-t, -x, -y, -z);
	}

	// Apply coordinate transformation to objects.
	Vec4 transform(Vec4 const& other) const;
	Transform4 transform(Transform3 const& other) const;
	Transform4 transform(Transform4 const& other) const;
};

Transform4 const TRANSFORM4_ZERO(
	0., 0., 0., 0.,
	0., 0., 0., 0.,
	0., 0., 0., 0.,
	0., 0., 0., 0.);
Transform4 const TRANSFORM4_ID(
	1., 0., 0., 0.,
	0., 1., 0., 0.,
	0., 0., 1., 0.,
	0., 0., 0., 1.);

inline Transform4 operator+(Transform4 lhs, Transform4 const& rhs) {
	lhs += rhs;
	return lhs;
}
inline Transform4 operator-(Transform4 lhs, Transform4 const& rhs) {
	lhs -= rhs;
	return lhs;
}
inline Transform4 operator*(Transform4 lhs, Real const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline Transform4 operator*(Real const& lhs, Transform4 rhs) {
	rhs *= lhs;
	return rhs;
}
inline Transform4 operator/(Transform4 lhs, Real const& rhs) {
	lhs /= rhs;
	return lhs;
}
inline Transform4 operator*(Transform4 lhs, Transform4 const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline Vec4 operator*(Transform4 lhs, Vec4 rhs) {
	return Vec4(
		dot(lhs.t, rhs),
		dot(lhs.x, rhs),
		dot(lhs.y, rhs),
		dot(lhs.z, rhs));
}

inline Vec4 dot(Transform4 const& lhs, Vec4 const& rhs) {
	return lhs * rhs;
}
inline Vec4 dot(Vec4 const& lhs, Transform4 const& rhs) {
	return rhs.transpose() * lhs;
}
inline Transform4 dot(Transform4 const& lhs, Transform4 const& rhs) {
	return lhs * rhs;
}
inline Transform4 outer(Vec4 const& lhs, Vec4 const& rhs) {
	return Transform4(lhs.t * rhs, lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}
inline Transform4 cross(Vec4 const& a, Vec4 const& b) {
	return Transform4(
		0., a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x,
		a.z * b.y - a.y * b.z, 0., a.t * b.z - a.z * b.t, a.y * b.t - a.t * b.y,
		a.x * b.z - a.z * b.x, a.z * b.t - a.t * b.z, 0., a.t * b.x - a.x * b.t,
		a.y * b.x - a.x * b.y, a.t * b.y - a.y * b.t, a.x * b.t - a.t * b.x, 0.);
}
inline Transform4 asym_prod(Vec4 const& v1, Vec4 const& v2) {
	return outer(v1, v2) - outer(v2, v1);
}
inline Transform4 sym_prod(Vec4 const& v1, Vec4 const& v2) {
	return outer(v1, v2) + outer(v1, v2);
}

/// \}

}
}

#endif

