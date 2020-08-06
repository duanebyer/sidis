#include "sidis/vector.hpp"

#include <algorithm>
#include <cmath>

using namespace sidis;
using namespace sidis::math;

Vec3 const Vec3::ZERO = Vec3();
Vec3 const Vec3::X = Vec3(1., 0., 0.);
Vec3 const Vec3::Y = Vec3(0., 1., 0.);
Vec3 const Vec3::Z = Vec3(0., 0., 1.);

Real Vec3::norm_sq() const {
	return dot(*this, *this);
}

Real Vec3::norm() const {
	// This approach ensures that we don't overflow/underflow for large vectors.
	Real max = std::max({
		std::abs(this->x),
		std::abs(this->y),
		std::abs(this->z) });
	Real x_r = this->x / max;
	Real y_r = this->y / max;
	Real z_r = this->z / max;
	return max * std::sqrt(x_r * x_r + y_r * y_r + z_r * z_r);
}

Vec3 Vec3::unit() const {
	Real norm = this->norm();
	if (norm == 0.) {
		return Vec3::ZERO;
	} else {
		return *this / norm;
	}
}

Vec3 Vec3::par(Vec3 const& v) const {
	Real norm_sq = v.norm_sq();
	if (norm_sq == 0.) {
		return Vec3::ZERO;
	} else {
		return dot(*this, v) * v / norm_sq;
	}
}

Vec3 Vec3::perp(Vec3 const& v) const {
	return *this - this->par(v);
}

Real math::angle_between(Vec3 const& vec_a, Vec3 const& vec_b) {
	return std::acos(dot(vec_a.unit(), vec_b.unit()));
}

Vec4 const Vec4::ZERO = Vec4();
Vec4 const Vec4::T = Vec4(1., 0., 0., 0.);
Vec4 const Vec4::X = Vec4(0., 1., 0., 0.);
Vec4 const Vec4::Y = Vec4(0., 0., 1., 0.);
Vec4 const Vec4::Z = Vec4(0., 0., 0., 1.);

Vec4 math::cross(Vec4 const& a, Vec4 const& b, Vec4 const& c) {
	return Vec4(
		+a.r.x * b.r.y * c.r.z + a.r.y * b.r.z * c.r.x + a.r.z * b.r.x * c.r.y
		-a.r.x * b.r.z * c.r.y - a.r.y * b.r.x * c.r.z - a.r.z * b.r.y * c.r.x,
		+  a.t * b.r.z * c.r.x + a.r.y *   b.t * c.r.z + a.r.z * b.r.y *   c.t
		-  a.t * b.r.x * c.r.z - a.r.y * b.r.z *   c.t - a.r.z *   b.t * c.r.y,
		+  a.t * b.r.x * c.r.z + a.r.x * b.r.z *   c.t + a.r.z *   b.t * c.r.x
		-  a.t * b.r.z * c.r.x - a.r.x *   b.t * c.r.z - a.r.z * b.r.x *   c.t,
		+  a.t * b.r.y * c.r.x + a.r.x *   b.t * c.r.y + a.r.y * b.r.x *   c.t
		-  a.t * b.r.x * c.r.y - a.r.x * b.r.y *   c.t - a.r.y *   b.t * c.r.x);
}

Real Vec4::norm_sq() const {
	return dot(*this, *this);
}

Real Vec4::norm() const {
	Real max = std::max({
		std::abs(this->t),
		std::abs(this->r.x),
		std::abs(this->r.y),
		std::abs(this->r.z) });
	Real t_r = this->t / max;
	Real x_r = this->r.x / max;
	Real y_r = this->r.y / max;
	Real z_r = this->r.z / max;
	return max * std::sqrt(t_r * t_r - (x_r * x_r + y_r * y_r + z_r * z_r));
}

Vec4 Vec4::unit() const {
	Real norm = this->norm();
	if (norm == 0.) {
		return Vec4::ZERO;
	} else {
		return *this / norm;
	}
}

Vec4 Vec4::par(Vec4 const& v) const {
	Real norm_sq = v.norm_sq();
	if (norm_sq == 0.) {
		return Vec4::ZERO;
	} else {
		return dot(*this, v) * v / norm_sq;
	}
}

Vec4 Vec4::perp(Vec4 const& v) const {
	return *this - this->par(v);
}

