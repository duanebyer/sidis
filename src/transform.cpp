#include "sidis/transform.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

Transform3 Transform3::rotate(Vec3 const& dir, Real angle) {
	Vec3 dir_unit = dir.unit();
	Real cos = std::cos(angle);
	Real sin = std::sin(angle);
	return cos * TRANSFORM3_ID
		+ sin * cross(dir_unit)
		+ (1. - cos) * outer(dir_unit, dir_unit);
}

Transform3 Transform3::rotate_to(Vec3 const& dir_old, Vec3 const& dir_new) {
	Vec3 dir_old_unit = dir_old.unit();
	Vec3 dir_new_unit = dir_new.unit();
	Vec3 dir = cross(dir_old_unit, dir_new_unit).unit();
	Real cos = dot(dir_old_unit, dir_new_unit);
	Real sin = cross(dir_old_unit, dir_new_unit).norm();
	return cos * TRANSFORM3_ID
		+ sin * cross(dir)
		+ (1. - cos) * outer(dir, dir);
}

Transform3 Transform3::rotate_to(Vec3 const& z_new) {
	return Transform3::rotate_to(VEC3_Z, z_new);
}

Transform3 Transform3::rotate_basis(Vec3 const& z_axis, Vec3 const& y_up) {
	Vec3 z = z_axis.unit();
	Vec3 x = cross(y_up, z_axis).unit();
	Vec3 y = cross(z, x).unit();
	return Transform3(x, y, z);
}

Transform3 Transform3::scale(Vec3 const& dir, Real scale) {
	Vec3 dir_unit = dir.unit();
	return TRANSFORM3_ID + (scale - 1.) * outer(dir_unit, dir_unit);
}

Transform3 Transform3::project(Vec3 const& dir) {
	Real n_sq = dir.norm_sq();
	return outer(dir, dir) / n_sq;
}

Transform3 Transform3::bivec(Vec3 const& dir) {
	return cross(dir);
}

Transform3 Transform3::transpose() const {
	return Transform3(
		Vec3(x.x, y.x, z.x),
		Vec3(x.y, y.y, z.y),
		Vec3(x.z, y.z, z.z));
}

Real Transform3::trace() const {
	return x.x + y.y + z.z;
}

Real Transform3::det() const {
	return x.x * y.y * z.z + x.y * y.z * z.x + x.z * y.x * z.y
		- x.x * y.z * z.y - x.z * y.y * z.x - x.y * y.x * z.z;
}

Transform3 Transform3::inv() const {
	return 1. / det() * Transform3(
		y.y * z.z - y.z * z.y, x.z * z.y - x.y * z.z, x.y * y.z - x.z * y.y,
		y.z * z.x - y.x * z.z, x.x * z.z - x.z * z.x, x.z * y.x - x.x * y.z,
		y.x * z.y - y.y * z.x, x.y * z.x - x.x * z.y, x.x * y.y - x.y * y.x);
}

Vec3 Transform3::transform(Vec3 const& other) const {
	return dot(*this, other);
}

Vec4 Transform3::transform(Vec4 const& other) const {
	return Transform4(*this).transform(other);
}

Transform3 Transform3::transform(Transform3 const& other) const {
	return dot(dot(*this, other), this->transpose());
}

Transform4 Transform3::transform(Transform4 const& other) const {
	return Transform4(*this).transform(other);
}

Transform4 Transform4::rotate(Vec3 const& dir, Real angle) {
	return Transform4(Transform3::rotate(dir, angle));
}
Transform4 Transform4::rotate_to(Vec3 const& dir_old, Vec3 const& dir_new) {
	return Transform4(Transform3::rotate_to(dir_old, dir_new));
}
Transform4 Transform4::rotate_to(Vec3 const& z_axis) {
	return Transform4(Transform3::rotate_to(z_axis));
}
Transform4 Transform4::rotate_basis(Vec3 const& z_axis, Vec3 const& y_up) {
	return Transform4(Transform3::rotate_basis(z_axis, y_up));
}
Transform4 Transform4::scale(Vec3 const& dir, Real scale) {
	return Transform4(Transform3::scale(dir, scale));
}
Transform4 Transform4::project(Vec3 const& dir) {
	return Transform4(Transform3::project(dir));
}

Transform4 Transform4::boost(Vec3 const& dir, Real rapidity) {
	Vec4 dir_unit = Vec4(0., dir.unit());
	Real cosh = std::cosh(rapidity);
	Real sinh = std::sinh(rapidity);
	return TRANSFORM4_ID
		+ sinh * (outer(dir_unit, VEC4_T) - outer(VEC4_T, dir_unit))
		+ (1. - cosh) * (outer(dir_unit, dir_unit) - outer(VEC4_T, VEC4_T));
}

Transform4 Transform4::boost_to(Vec4 const& t_new) {
	return Transform4::transform_to(t_new);
}

Transform4 Transform4::transform_to(Vec4 const& dir_old, Vec4 const& dir_new) {
	Vec4 dir_old_unit = dir_old.unit();
	Vec4 dir_new_unit = dir_new.unit();
	Real cos = dot(dir_old_unit, dir_new_unit);
	int s_old = dir_old_unit.sign();
	int s = dir_new_unit.sign();
	if (s_old != s || s == 0) {
		return std::numeric_limits<Real>::quiet_NaN() * TRANSFORM4_ID;
	}
	Transform4 sym = outer(dir_old_unit, dir_old_unit)
		+ outer(dir_new_unit, dir_new_unit);
	Transform4 asym = outer(dir_old_unit, dir_new_unit)
		- outer(dir_new_unit, dir_old_unit);
	Transform4 transport = 2. * s * cos * outer(dir_new_unit, dir_old_unit);
	return TRANSFORM4_ID - (sym + asym - transport) / (s + cos);
}

Transform4 Transform4::transform_to(Vec4 const& t_new) {
	return Transform4::transform_to(VEC4_T, t_new);
}

Transform4 Transform4::project(Vec4 const& dir) {
	Real n_sq = dir.norm_sq();
	return outer(dir, dir) / n_sq;
}

Transform4 Transform4::bivec(Vec3 const& v, Vec3 const& bv) {
	return Transform4(
		0., -v.x, -v.y, -v.z,
		-v.x, 0., -bv.z, bv.y,
		-v.y, bv.z, 0., -bv.x,
		-v.z, -bv.y, bv.x, 0.);

}

Transform4 Transform4::transpose() const {
	return Transform4(
		Vec4(t.t, x.t, y.t, z.t),
		Vec4(t.x, x.x, y.x, z.x),
		Vec4(t.y, x.y, y.y, z.y),
		Vec4(t.z, x.z, y.z, z.z));
}

Real Transform4::trace() const {
	return t.t - x.x - y.y - z.z;
}

Real Transform4::det() const {
	return -t.t * x.x * y.y * z.z - t.t * x.y * y.z * z.x - t.t * x.z * y.x * z.y
		+ t.t * x.x * y.z * z.y + t.t * x.z * y.y * z.x + t.t * x.y * y.x * z.z
		+ t.x * x.t * y.y * z.z + t.x * x.y * y.z * z.t + t.x * x.z * y.t * z.y
		- t.x * x.t * y.z * z.y - t.x * x.z * y.y * z.t - t.x * x.y * y.t * z.z
		+ t.y * x.x * y.t * z.z + t.y * x.t * y.z * z.x + t.y * x.z * y.x * z.t
		- t.y * x.x * y.z * z.t - t.y * x.z * y.t * z.x - t.y * x.t * y.x * z.z
		+ t.z * x.x * y.y * z.t + t.z * x.y * y.t * z.x + t.z * x.t * y.x * z.y
		- t.z * x.x * y.t * z.y - t.z * x.t * y.y * z.x - t.z * x.y * y.x * z.t;
}

Transform4 Transform4::inv() const {
	return 1. / det() * Transform4(
		-(x.x * y.y * z.z + x.y * y.z * z.x + x.z * y.x * z.y - x.x * y.z * z.y - x.z * y.y * z.x - x.y * y.x * z.z),
		-(t.x * y.z * z.y + t.z * y.y * z.x + t.y * y.x * z.z - t.x * y.y * z.z - t.y * y.z * z.x - t.z * y.x * z.y),
		-(x.x * t.z * z.y + x.z * t.y * z.x + x.y * t.x * z.z - x.x * t.y * z.z - x.y * t.z * z.x - x.z * t.x * z.y),
		-(x.x * y.z * t.y + x.z * y.y * t.x + x.y * y.x * t.z - x.x * y.y * t.z - x.y * y.z * t.x - x.z * y.x * t.y),
		x.t * y.z * z.y + x.z * y.y * z.t + x.y * y.t * z.z - x.t * y.y * z.z - x.y * y.z * z.t - x.z * y.t * z.y,
		t.t * y.y * z.z + t.y * y.z * z.t + t.z * y.t * z.y - t.t * y.z * z.y - t.z * y.y * z.t - t.y * y.t * z.z,
		t.t * x.z * z.y + t.z * x.y * z.t + t.y * x.t * z.z - t.t * x.y * z.z - t.y * x.z * z.t - t.z * x.t * z.y,
		t.t * y.z * x.y + t.z * y.y * x.t + t.y * y.t * x.z - t.t * y.y * x.z - t.y * y.z * x.t - t.z * y.t * x.y,
		x.x * y.z * z.t + x.z * y.t * z.x + x.t * y.x * z.z - x.x * y.t * z.z - x.t * y.z * z.x - x.z * y.x * z.t,
		y.x * t.z * z.t + y.z * t.t * z.x + y.t * t.x * z.z - y.x * t.t * z.z - y.t * t.z * z.x - y.z * t.x * z.t,
		x.x * t.t * z.z + x.t * t.z * z.x + x.z * t.x * z.t - x.x * t.z * z.t - x.z * t.t * z.x - x.t * t.x * z.z,
		x.x * t.z * y.t + x.z * t.t * y.x + x.t * t.x * y.z - x.x * t.t * y.z - x.t * t.z * y.x - x.z * t.x * y.t,
		x.x * y.t * z.y + x.t * y.y * z.x + x.y * y.x * z.t - x.x * y.y * z.t - x.y * y.t * z.x - x.t * y.x * z.y,
		z.x * y.t * t.y + z.t * y.y * t.x + z.y * y.x * t.t - z.x * y.y * t.t - z.y * y.t * t.x - z.t * y.x * t.y,
		x.x * z.t * t.y + x.t * z.y * t.x + x.y * z.x * t.t - x.x * z.y * t.t - x.y * z.t * t.x - x.t * z.x * t.y,
		x.x * y.y * t.t + x.y * y.t * t.x + x.t * y.x * t.y - x.x * y.t * t.y - x.t * y.y * t.x - x.y * y.x * t.t);
}

Vec4 Transform4::transform(Vec4 const& other) const {
	return dot(*this, other);
}

Transform4 Transform4::transform(Transform3 const& other) const {
	return transform(Transform4(other));
}

Transform4 Transform4::transform(Transform4 const& other) const {
	return dot(dot(*this, other), this->transpose());
}

