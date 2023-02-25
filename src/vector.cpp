#include "sidis/vector.hpp"

#include <algorithm>
#include <cmath>

using namespace sidis;
using namespace sidis::math;

Real Vec3::norm_sq() const {
	return dot(*this, *this);
}

Real Vec3::norm() const {
	// This approach ensures that we don't overflow/underflow for large vectors.
	Real max = std::max({
		std::abs(x),
		std::abs(y),
		std::abs(z) });
	if (max == 0.) {
		return 0.;
	} else {
		Real x_r = x / max;
		Real y_r = y / max;
		Real z_r = z / max;
		return max * std::sqrt(x_r * x_r + y_r * y_r + z_r * z_r);
	}
}

Vec3 Vec3::unit() const {
	Real n = norm();
	if (n == 0.) {
		return VEC3_ZERO;
	} else {
		return *this / n;
	}
}

Vec3 Vec3::par(Vec3 const& v) const {
	Real n_sq = v.norm_sq();
	if (n_sq == 0.) {
		return *this;
	} else {
		return dot(*this, v) * v / n_sq;
	}
}

Vec3 Vec3::perp(Vec3 const& v) const {
	Real n_sq = v.norm_sq();
	if (n_sq == 0.) {
		return VEC3_ZERO;
	} else {
		return cross(cross(v, *this), v) / n_sq;
	}
}

Real math::angle_between(Vec3 const& v_a, Vec3 const& v_b) {
	return std::acos(dot(v_a.unit(), v_b.unit()));
}

Vec4 Vec4::from_length_and_r(Real m, Vec3 const& p) {
	Real max = std::max({
		std::abs(m),
		std::abs(p.x),
		std::abs(p.y),
		std::abs(p.z) });
	Real m_r = m / max;
	Real x_r = p.x / max;
	Real y_r = p.y / max;
	Real z_r = p.z / max;
	Real energy = max * std::sqrt(m_r * m_r + x_r * x_r + y_r * y_r + z_r * z_r);
	return Vec4(energy, p);
}

Vec4 Vec4::from_length_and_t(Real m, Real t, Vec3 const& dir) {
	Real max = std::max({ std::abs(m), std::abs(t) });
	Real m_r = m / max;
	Real t_r = t / max;
	Real p = max * std::sqrt(t_r * t_r - m_r * m_r);
	return Vec4(t, p * dir.unit());
}

Vec4 math::cross(Vec4 const& a, Vec4 const& b, Vec4 const& c) {
	return Vec4(
		+a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y
		-a.x * b.z * c.y - a.y * b.x * c.z - a.z * b.y * c.x,
		+a.t * b.z * c.y + a.y * b.t * c.z + a.z * b.y * c.t
		-a.t * b.y * c.z - a.y * b.z * c.t - a.z * b.t * c.y,
		+a.t * b.x * c.z + a.x * b.z * c.t + a.z * b.t * c.x
		-a.t * b.z * c.x - a.x * b.t * c.z - a.z * b.x * c.t,
		+a.t * b.y * c.x + a.x * b.t * c.y + a.y * b.x * c.t
		-a.t * b.x * c.y - a.x * b.y * c.t - a.y * b.t * c.x);
}

Real Vec4::norm_sq() const {
	return dot(*this, *this);
}

Real Vec4::norm() const {
	Real max = std::max({ std::abs(t), std::abs(x), std::abs(y), std::abs(z) });
	if (max == 0.) {
		return 0.;
	} else {
		Real t_r = t / max;
		Real x_r = x / max;
		Real y_r = y / max;
		Real z_r = z / max;
		return max * std::sqrt(
			std::abs(t_r * t_r - (x_r * x_r + y_r * y_r + z_r * z_r)));
	}
}

Vec4 Vec4::unit() const {
	Real n = norm();
	if (n == 0.) {
		return VEC4_ZERO;
	} else {
		return *this / n;
	}
}

int Vec4::sign() const {
	Real n_sq = norm_sq();
	return (0. < n_sq) - (n_sq < 0.);
}

Vec4 Vec4::par(Vec4 const& v) const {
	Real n_sq = v.norm_sq();
	if (n_sq == 0.) {
		return *this;
	} else {
		return dot(*this, v) * v / n_sq;
	}
}

Vec4 Vec4::perp(Vec4 const& v) const {
	return *this - par(v);
}

