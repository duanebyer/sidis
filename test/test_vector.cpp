#include <catch2/catch.hpp>

#include <cmath>
#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/constant.hpp>
#include <sidis/transform.hpp>
#include <sidis/vector.hpp>

// Declare stream operators before including templates, since ADL fails us here.
namespace {
std::istream& operator>>(std::istream& in, sidis::math::Vec3& vec) {
	in >> vec.x;
	in >> vec.y;
	in >> vec.z;
	return in;
}
std::istream& operator>>(std::istream& in, sidis::math::Vec4& vec) {
	in >> vec.t;
	in >> vec.x;
	in >> vec.y;
	in >> vec.z;
	return in;
}
}

#include "abs_matcher.hpp"
#include "rel_matcher.hpp"
#include "stream_generator.hpp"

using namespace sidis;

TEST_CASE(
		"3-vector basic operations checks",
		"[vec]") {
	math::Vec3 vec = GENERATE(
		from_stream<math::Vec3>(
			std::move(std::ifstream("data/vec3_vals.dat")),
			true));

	std::stringstream ss;
	ss << "v = (" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
	INFO(ss.str());

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	Real norm = vec.norm();
	if (norm == 0.) {
		CHECK((vec.x == 0. && vec.y == 0. && vec.z == 0.));
	} else {
		CHECK_THAT(
			vec.unit().norm(),
			RelMatcher<Real>(1., prec));
	}
	CHECK_THAT(
		(2. * vec).norm(),
		RelMatcher<Real>(2. * norm, prec));
	CHECK_THAT(
		(0.5 * vec).norm(),
		RelMatcher<Real>(0.5 * norm, prec));
	CHECK_THAT(
		(vec + 0.5 * vec).norm(),
		RelMatcher<Real>(1.5 * norm, prec));
	CHECK_THAT(
		dot(vec, vec),
		RelMatcher<Real>(norm * norm, prec));
	CHECK_THAT(
		cross(vec, vec).norm(),
		AbsMatcher<Real>(0., prec * norm));
	CHECK_THAT(
		vec.perp(vec).norm(),
		AbsMatcher<Real>(0., prec * norm));
	CHECK_THAT(
		vec.par(vec).x,
		RelMatcher<Real>(vec.x, prec));
	CHECK_THAT(
		vec.par(vec).y,
		RelMatcher<Real>(vec.y, prec));
	CHECK_THAT(
		vec.par(vec).z,
		RelMatcher<Real>(vec.z, prec));
	CHECK_THAT(
		vec.par(10. * math::VEC3_X).norm(),
		RelMatcher<Real>(std::abs(vec.x), prec));
	CHECK_THAT(
		vec.par(10. * math::VEC3_Y).norm(),
		RelMatcher<Real>(std::abs(vec.y), prec));
	CHECK_THAT(
		vec.par(10. * math::VEC3_Z).norm(),
		RelMatcher<Real>(std::abs(vec.z), prec));
}

TEST_CASE(
		"3-vector rotation checks",
		"[vec]") {
	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();
	math::Vec3 test;

	math::Transform3 rotate_x = math::Transform3::rotate(
		math::Vec3(1., 0., 0.),
		PI / 2.);
	math::Vec3 dir(1., 1., 0.);
	math::Transform3 rotate_dir = math::Transform3::rotate(
		dir,
		-PI / 3.);

	// Apply simple `x` rotation to some vectors.
	test = rotate_x * math::VEC3_X;
	CHECK_THAT(test.x, RelMatcher<Real>(1., prec));
	CHECK_THAT(test.y, RelMatcher<Real>(0., prec));
	CHECK_THAT(test.z, RelMatcher<Real>(0., prec));

	test = rotate_x * math::Vec3(1., 2., 3.);
	CHECK_THAT(test.x, RelMatcher<Real>(1., prec));
	CHECK_THAT(test.y, RelMatcher<Real>(-3., prec));
	CHECK_THAT(test.z, RelMatcher<Real>(2., prec));
	CHECK_THAT(test.norm_sq(), RelMatcher<Real>(14., prec));

	// Check that the direction of rotation is an eigenvector.
	test = rotate_dir * dir;
	CHECK_THAT(test.x, RelMatcher<Real>(dir.x, prec));
	CHECK_THAT(test.y, RelMatcher<Real>(dir.y, prec));
	CHECK_THAT(test.z, RelMatcher<Real>(dir.z, prec));

	// Check determinant of 1.
	CHECK_THAT(rotate_dir.det(), RelMatcher<Real>(1., prec));

	// Check that the rotation matrix is orthogonal.
	math::Transform3 rotate_dir_inv = rotate_dir.inv();
	math::Transform3 rotate_dir_tr = rotate_dir.transpose();
	CHECK_THAT(rotate_dir_inv.x.x, RelMatcher<Real>(rotate_dir_tr.x.x, prec));
	CHECK_THAT(rotate_dir_inv.x.y, RelMatcher<Real>(rotate_dir_tr.x.y, prec));
	CHECK_THAT(rotate_dir_inv.x.z, RelMatcher<Real>(rotate_dir_tr.x.z, prec));
	CHECK_THAT(rotate_dir_inv.y.x, RelMatcher<Real>(rotate_dir_tr.y.x, prec));
	CHECK_THAT(rotate_dir_inv.y.y, RelMatcher<Real>(rotate_dir_tr.y.y, prec));
	CHECK_THAT(rotate_dir_inv.y.z, RelMatcher<Real>(rotate_dir_tr.y.z, prec));
	CHECK_THAT(rotate_dir_inv.z.x, RelMatcher<Real>(rotate_dir_tr.z.x, prec));
	CHECK_THAT(rotate_dir_inv.z.y, RelMatcher<Real>(rotate_dir_tr.z.y, prec));
	CHECK_THAT(rotate_dir_inv.z.z, RelMatcher<Real>(rotate_dir_tr.z.z, prec));
}

TEST_CASE(
		"3-vector rotate-to checks",
		"[vec]") {
	math::Vec3 vec = GENERATE(
		from_stream<math::Vec3>(
			std::move(std::ifstream("data/vec3_vals.dat")),
			true));

	std::stringstream ss;
	ss << "v = (" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
	INFO(ss.str());

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	if (vec.norm_sq() != 0.) {
		math::Vec3 vec_from(-0.34, 0.96, 1.2);
		math::Transform3 rotate_1 = math::Transform3::rotate_to(vec);
		math::Transform3 rotate_2 = math::Transform3::rotate_to(vec_from, vec);

		// Apply to the `from` vector.
		math::Vec3 result_1 = rotate_1 * math::VEC3_Z;
		math::Vec3 result_2 = rotate_2 * vec_from.unit();
		CHECK_THAT(result_1.x, AbsMatcher<Real>(vec.unit().x, prec * vec.norm()));
		CHECK_THAT(result_1.y, AbsMatcher<Real>(vec.unit().y, prec * vec.norm()));
		CHECK_THAT(result_1.z, AbsMatcher<Real>(vec.unit().z, prec * vec.norm()));
		CHECK_THAT(result_2.x, AbsMatcher<Real>(vec.unit().x, prec * vec.norm()));
		CHECK_THAT(result_2.y, AbsMatcher<Real>(vec.unit().y, prec * vec.norm()));
		CHECK_THAT(result_2.z, AbsMatcher<Real>(vec.unit().z, prec * vec.norm()));

		// Check composition of transformations.
		math::Vec3 test(1., 2., 3.);
		math::Vec3 apply_1 = (rotate_1 * rotate_2) * test;
		math::Vec3 apply_2 = rotate_1 * (rotate_2 * test);
		CHECK_THAT(apply_1.x, RelMatcher<Real>(apply_2.x, prec));
		CHECK_THAT(apply_1.y, RelMatcher<Real>(apply_2.y, prec));
		CHECK_THAT(apply_1.z, RelMatcher<Real>(apply_2.z, prec));

		math::Vec3 apply_3 = dot(test, dot(rotate_1, rotate_2));
		math::Vec3 apply_4 = dot(dot(test, rotate_1), rotate_2);
		CHECK_THAT(apply_3.x, RelMatcher<Real>(apply_4.x, prec));
		CHECK_THAT(apply_3.y, RelMatcher<Real>(apply_4.y, prec));
		CHECK_THAT(apply_3.z, RelMatcher<Real>(apply_4.z, prec));

		// Check determinant of 1.
		CHECK_THAT(rotate_1.det(), RelMatcher<Real>(1., prec));
		CHECK_THAT(rotate_2.det(), RelMatcher<Real>(1., prec));

		// Check that the rotation matrix is orthogonal.
		math::Transform3 rotate_1_inv = rotate_1.inv();
		math::Transform3 rotate_1_tr = rotate_1.transpose();
		CHECK_THAT(rotate_1_inv.x.x, RelMatcher<Real>(rotate_1_tr.x.x, prec));
		CHECK_THAT(rotate_1_inv.x.y, RelMatcher<Real>(rotate_1_tr.x.y, prec));
		CHECK_THAT(rotate_1_inv.x.z, RelMatcher<Real>(rotate_1_tr.x.z, prec));
		CHECK_THAT(rotate_1_inv.y.x, RelMatcher<Real>(rotate_1_tr.y.x, prec));
		CHECK_THAT(rotate_1_inv.y.y, RelMatcher<Real>(rotate_1_tr.y.y, prec));
		CHECK_THAT(rotate_1_inv.y.z, RelMatcher<Real>(rotate_1_tr.y.z, prec));
		CHECK_THAT(rotate_1_inv.z.x, RelMatcher<Real>(rotate_1_tr.z.x, prec));
		CHECK_THAT(rotate_1_inv.z.y, RelMatcher<Real>(rotate_1_tr.z.y, prec));
		CHECK_THAT(rotate_1_inv.z.z, RelMatcher<Real>(rotate_1_tr.z.z, prec));
		math::Transform3 rotate_2_inv = rotate_2.inv();
		math::Transform3 rotate_2_tr = rotate_2.transpose();
		CHECK_THAT(rotate_2_inv.x.x, RelMatcher<Real>(rotate_2_tr.x.x, prec));
		CHECK_THAT(rotate_2_inv.x.y, RelMatcher<Real>(rotate_2_tr.x.y, prec));
		CHECK_THAT(rotate_2_inv.x.z, RelMatcher<Real>(rotate_2_tr.x.z, prec));
		CHECK_THAT(rotate_2_inv.y.x, RelMatcher<Real>(rotate_2_tr.y.x, prec));
		CHECK_THAT(rotate_2_inv.y.y, RelMatcher<Real>(rotate_2_tr.y.y, prec));
		CHECK_THAT(rotate_2_inv.y.z, RelMatcher<Real>(rotate_2_tr.y.z, prec));
		CHECK_THAT(rotate_2_inv.z.x, RelMatcher<Real>(rotate_2_tr.z.x, prec));
		CHECK_THAT(rotate_2_inv.z.y, RelMatcher<Real>(rotate_2_tr.z.y, prec));
		CHECK_THAT(rotate_2_inv.z.z, RelMatcher<Real>(rotate_2_tr.z.z, prec));
	}
}

TEST_CASE(
		"4-vector basic operations checks",
		"[vec]") {
	math::Vec4 vec = GENERATE(
		from_stream<math::Vec4>(
			std::move(std::ifstream("data/vec4_vals.dat")),
			true));

	std::stringstream ss;
	ss << "v = ("
		<< vec.t << ", "
		<< vec.x << ", "
		<< vec.y << ", "
		<< vec.z << ")" << std::endl;
	INFO(ss.str());

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	Real norm = vec.norm();
	if (norm == 0.) {
		CHECK_THAT(
			vec.r().norm(),
			RelMatcher<Real>(vec.t, prec));
	} else {
		CHECK_THAT(
			vec.unit().norm(),
			RelMatcher<Real>(1., prec));
	}
	CHECK_THAT(
		(2. * vec).norm(),
		RelMatcher<Real>(2. * norm, prec));
	CHECK_THAT(
		(0.5 * vec).norm(),
		RelMatcher<Real>(0.5 * norm, prec));
	CHECK_THAT(
		(vec + 0.5 * vec).norm(),
		RelMatcher<Real>(1.5 * norm, prec));
	CHECK_THAT(
		dot(vec, vec),
		RelMatcher<Real>(norm * norm, prec));
}

TEST_CASE(
		"4-vector boost checks",
		"[vec]") {
	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();
	math::Vec4 test;

	math::Transform4 boost_x = math::Transform4::boost(
		math::Vec3(1., 0., 0.),
		std::log(2.));
	math::Vec3 dir(1., 1., 0.);
	math::Transform4 boost_dir = math::Transform4::boost(
		dir,
		std::log(1.5));

	// Apply simple `x` boost to some vectors.
	test = boost_x * math::VEC4_T;
	CHECK_THAT(test.t, RelMatcher<Real>(1.25, prec));
	CHECK_THAT(test.x, RelMatcher<Real>(0.75, prec));
	CHECK_THAT(test.y, RelMatcher<Real>(0., prec));
	CHECK_THAT(test.z, RelMatcher<Real>(0., prec));

	test = boost_x * math::Vec4(10., 3., 4., 5.);
	CHECK_THAT(test.t, RelMatcher<Real>(10. * 1.25 + 3. * 0.75, prec));
	CHECK_THAT(test.x, RelMatcher<Real>(3. * 1.25 + 10. * 0.75, prec));
	CHECK_THAT(test.y, RelMatcher<Real>(4., prec));
	CHECK_THAT(test.z, RelMatcher<Real>(5., prec));
	CHECK_THAT(test.norm_sq(), RelMatcher<Real>(50., prec));

	// Check that the light path in the direction of the boost is an
	// eigenvector.
	test = boost_dir * math::Vec4(dir.norm(), dir);
	CHECK_THAT(test.x / test.t, RelMatcher<Real>(dir.x / dir.norm(), prec));
	CHECK_THAT(test.y / test.t, RelMatcher<Real>(dir.y / dir.norm(), prec));
	CHECK_THAT(test.z / test.t, RelMatcher<Real>(dir.z / dir.norm(), prec));

	// Check determinant of 1.
	CHECK_THAT(boost_dir.det(), RelMatcher<Real>(1., prec));

	// Check that the boost matrix is orthogonal.
	math::Transform4 boost_dir_inv = boost_dir.inv();
	math::Transform4 boost_dir_tr = boost_dir.transpose();
	CHECK_THAT(boost_dir_inv.t.t, RelMatcher<Real>(boost_dir_tr.t.t, prec));
	CHECK_THAT(boost_dir_inv.t.x, RelMatcher<Real>(boost_dir_tr.t.x, prec));
	CHECK_THAT(boost_dir_inv.t.y, RelMatcher<Real>(boost_dir_tr.t.y, prec));
	CHECK_THAT(boost_dir_inv.t.z, RelMatcher<Real>(boost_dir_tr.t.z, prec));
	CHECK_THAT(boost_dir_inv.x.t, RelMatcher<Real>(boost_dir_tr.x.t, prec));
	CHECK_THAT(boost_dir_inv.x.x, RelMatcher<Real>(boost_dir_tr.x.x, prec));
	CHECK_THAT(boost_dir_inv.x.y, RelMatcher<Real>(boost_dir_tr.x.y, prec));
	CHECK_THAT(boost_dir_inv.x.z, RelMatcher<Real>(boost_dir_tr.x.z, prec));
	CHECK_THAT(boost_dir_inv.y.t, RelMatcher<Real>(boost_dir_tr.y.t, prec));
	CHECK_THAT(boost_dir_inv.y.x, RelMatcher<Real>(boost_dir_tr.y.x, prec));
	CHECK_THAT(boost_dir_inv.y.y, RelMatcher<Real>(boost_dir_tr.y.y, prec));
	CHECK_THAT(boost_dir_inv.y.z, RelMatcher<Real>(boost_dir_tr.y.z, prec));
	CHECK_THAT(boost_dir_inv.z.t, RelMatcher<Real>(boost_dir_tr.z.t, prec));
	CHECK_THAT(boost_dir_inv.z.x, RelMatcher<Real>(boost_dir_tr.z.x, prec));
	CHECK_THAT(boost_dir_inv.z.y, RelMatcher<Real>(boost_dir_tr.z.y, prec));
	CHECK_THAT(boost_dir_inv.z.z, RelMatcher<Real>(boost_dir_tr.z.z, prec));
}

TEST_CASE(
		"4-vector boost-to checks",
		"[vec]") {
	math::Vec4 vec = GENERATE(
		from_stream<math::Vec4>(
			std::move(std::ifstream("data/vec4_vals.dat")),
			true));

	std::stringstream ss;
	ss << "v = ("
		<< vec.t << ", "
		<< vec.x << ", "
		<< vec.y << ", "
		<< vec.z << ")" << std::endl;
	INFO(ss.str());

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	math::Vec4 vec_from(2.1, -0.34, 0.96, 1.2);
	math::Transform4 boost_1 = math::Transform4::boost_to(vec);
	math::Transform4 boost_2 = math::Transform4::transform_to(vec_from, vec);

	// Apply to the `from` vector.
	math::Vec4 result_1 = boost_1 * math::VEC4_T;
	math::Vec4 result_2 = boost_2 * vec_from.unit();
	CHECK_THAT(result_1.t, AbsMatcher<Real>(vec.unit().t, prec * vec.t));
	CHECK_THAT(result_1.x, AbsMatcher<Real>(vec.unit().x, prec * vec.t));
	CHECK_THAT(result_1.y, AbsMatcher<Real>(vec.unit().y, prec * vec.t));
	CHECK_THAT(result_1.z, AbsMatcher<Real>(vec.unit().z, prec * vec.t));
	CHECK_THAT(result_2.t, AbsMatcher<Real>(vec.unit().t, prec * vec.t));
	CHECK_THAT(result_2.x, AbsMatcher<Real>(vec.unit().x, prec * vec.t));
	CHECK_THAT(result_2.y, AbsMatcher<Real>(vec.unit().y, prec * vec.t));
	CHECK_THAT(result_2.z, AbsMatcher<Real>(vec.unit().z, prec * vec.t));

	// Check composition of transformations.
	math::Vec4 test(1., 0.1, 0.2, 0.3);
	math::Vec4 apply_1 = (boost_1 * boost_2) * test;
	math::Vec4 apply_2 = boost_1 * (boost_2 * test);
	CHECK_THAT(apply_1.t, RelMatcher<Real>(apply_2.t, prec));
	CHECK_THAT(apply_1.x, RelMatcher<Real>(apply_2.x, prec));
	CHECK_THAT(apply_1.y, RelMatcher<Real>(apply_2.y, prec));
	CHECK_THAT(apply_1.z, RelMatcher<Real>(apply_2.z, prec));

	math::Vec4 apply_3 = dot(test, dot(boost_1, boost_2));
	math::Vec4 apply_4 = dot(dot(test, boost_1), boost_2);
	CHECK_THAT(apply_3.t, RelMatcher<Real>(apply_4.t, prec));
	CHECK_THAT(apply_3.x, RelMatcher<Real>(apply_4.x, prec));
	CHECK_THAT(apply_3.y, RelMatcher<Real>(apply_4.y, prec));
	CHECK_THAT(apply_3.z, RelMatcher<Real>(apply_4.z, prec));

	// Check that the light path in the direction of the boost is an
	// eigenvector.
	if (vec.r().norm_sq() != 0.) {
		math::Vec4 light = math::Vec4(vec.r().norm(), vec.r());
		math::Vec4 light_1 = boost_1 * light;
		CHECK_THAT(light_1.x / light_1.t, RelMatcher<Real>(light.x / light.t, prec));
		CHECK_THAT(light_1.y / light_1.t, RelMatcher<Real>(light.y / light.t, prec));
		CHECK_THAT(light_1.z / light_1.t, RelMatcher<Real>(light.z / light.t, prec));
	}

	// Check determinant of 1.
	CHECK_THAT(boost_1.det(), RelMatcher<Real>(1., prec));
	CHECK_THAT(boost_2.det(), RelMatcher<Real>(1., prec));

	// Check that the boost matrix is orthogonal.
	math::Transform4 boost_1_inv = boost_1.inv();
	math::Transform4 boost_1_tr = boost_1.transpose();
	CHECK_THAT(boost_1_inv.t.t, RelMatcher<Real>(boost_1_tr.t.t, prec));
	CHECK_THAT(boost_1_inv.t.x, RelMatcher<Real>(boost_1_tr.t.x, prec));
	CHECK_THAT(boost_1_inv.t.y, RelMatcher<Real>(boost_1_tr.t.y, prec));
	CHECK_THAT(boost_1_inv.t.z, RelMatcher<Real>(boost_1_tr.t.z, prec));
	CHECK_THAT(boost_1_inv.x.t, RelMatcher<Real>(boost_1_tr.x.t, prec));
	CHECK_THAT(boost_1_inv.x.x, RelMatcher<Real>(boost_1_tr.x.x, prec));
	CHECK_THAT(boost_1_inv.x.y, RelMatcher<Real>(boost_1_tr.x.y, prec));
	CHECK_THAT(boost_1_inv.x.z, RelMatcher<Real>(boost_1_tr.x.z, prec));
	CHECK_THAT(boost_1_inv.y.t, RelMatcher<Real>(boost_1_tr.y.t, prec));
	CHECK_THAT(boost_1_inv.y.x, RelMatcher<Real>(boost_1_tr.y.x, prec));
	CHECK_THAT(boost_1_inv.y.y, RelMatcher<Real>(boost_1_tr.y.y, prec));
	CHECK_THAT(boost_1_inv.y.z, RelMatcher<Real>(boost_1_tr.y.z, prec));
	CHECK_THAT(boost_1_inv.z.t, RelMatcher<Real>(boost_1_tr.z.t, prec));
	CHECK_THAT(boost_1_inv.z.x, RelMatcher<Real>(boost_1_tr.z.x, prec));
	CHECK_THAT(boost_1_inv.z.y, RelMatcher<Real>(boost_1_tr.z.y, prec));
	CHECK_THAT(boost_1_inv.z.z, RelMatcher<Real>(boost_1_tr.z.z, prec));
	math::Transform4 boost_2_inv = boost_2.inv();
	math::Transform4 boost_2_tr = boost_2.transpose();
	CHECK_THAT(boost_2_inv.t.t, RelMatcher<Real>(boost_2_tr.t.t, prec));
	CHECK_THAT(boost_2_inv.t.x, RelMatcher<Real>(boost_2_tr.t.x, prec));
	CHECK_THAT(boost_2_inv.t.y, RelMatcher<Real>(boost_2_tr.t.y, prec));
	CHECK_THAT(boost_2_inv.t.z, RelMatcher<Real>(boost_2_tr.t.z, prec));
	CHECK_THAT(boost_2_inv.x.t, RelMatcher<Real>(boost_2_tr.x.t, prec));
	CHECK_THAT(boost_2_inv.x.x, RelMatcher<Real>(boost_2_tr.x.x, prec));
	CHECK_THAT(boost_2_inv.x.y, RelMatcher<Real>(boost_2_tr.x.y, prec));
	CHECK_THAT(boost_2_inv.x.z, RelMatcher<Real>(boost_2_tr.x.z, prec));
	CHECK_THAT(boost_2_inv.y.t, RelMatcher<Real>(boost_2_tr.y.t, prec));
	CHECK_THAT(boost_2_inv.y.x, RelMatcher<Real>(boost_2_tr.y.x, prec));
	CHECK_THAT(boost_2_inv.y.y, RelMatcher<Real>(boost_2_tr.y.y, prec));
	CHECK_THAT(boost_2_inv.y.z, RelMatcher<Real>(boost_2_tr.y.z, prec));
	CHECK_THAT(boost_2_inv.z.t, RelMatcher<Real>(boost_2_tr.z.t, prec));
	CHECK_THAT(boost_2_inv.z.x, RelMatcher<Real>(boost_2_tr.z.x, prec));
	CHECK_THAT(boost_2_inv.z.y, RelMatcher<Real>(boost_2_tr.z.y, prec));
	CHECK_THAT(boost_2_inv.z.z, RelMatcher<Real>(boost_2_tr.z.z, prec));
}

TEST_CASE(
		"4-vector spacelike transform-to checks",
		"[vec]") {
	math::Vec4 vec = GENERATE(
		from_stream<math::Vec4>(
			std::move(std::ifstream("data/vec4_spacelike_vals.dat")),
			true));

	std::stringstream ss;
	ss << "v = ("
		<< vec.t << ", "
		<< vec.x << ", "
		<< vec.y << ", "
		<< vec.z << ")" << std::endl;
	INFO(ss.str());

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	math::Vec4 vec_from(0.63, -0.34, 0.96, 1.2);
	math::Transform4 boost = math::Transform4::transform_to(vec_from, vec);

	// Apply to the `from` vector.
	math::Vec4 result = boost * vec_from.unit();
	CHECK_THAT(result.t, AbsMatcher<Real>(vec.unit().t, prec * vec.r().norm()));
	CHECK_THAT(result.x, AbsMatcher<Real>(vec.unit().x, prec * vec.r().norm()));
	CHECK_THAT(result.y, AbsMatcher<Real>(vec.unit().y, prec * vec.r().norm()));
	CHECK_THAT(result.z, AbsMatcher<Real>(vec.unit().z, prec * vec.r().norm()));

	// Check composition of transformations.
	math::Vec4 test(1., 0.1, 0.2, 0.3);
	math::Transform4 transform = math::Transform4::boost(
		math::Vec3(1., 1., 0.), 1.2);
	math::Vec4 apply_1 = (boost * transform) * test;
	math::Vec4 apply_2 = boost * (transform * test);
	CHECK_THAT(apply_1.t, RelMatcher<Real>(apply_2.t, prec));
	CHECK_THAT(apply_1.x, RelMatcher<Real>(apply_2.x, prec));
	CHECK_THAT(apply_1.y, RelMatcher<Real>(apply_2.y, prec));
	CHECK_THAT(apply_1.z, RelMatcher<Real>(apply_2.z, prec));

	math::Vec4 apply_3 = dot(test, dot(boost, transform));
	math::Vec4 apply_4 = dot(dot(test, boost), transform);
	CHECK_THAT(apply_3.t, RelMatcher<Real>(apply_4.t, prec));
	CHECK_THAT(apply_3.x, RelMatcher<Real>(apply_4.x, prec));
	CHECK_THAT(apply_3.y, RelMatcher<Real>(apply_4.y, prec));
	CHECK_THAT(apply_3.z, RelMatcher<Real>(apply_4.z, prec));

	// Check determinant of 1.
	CHECK_THAT(boost.det(), RelMatcher<Real>(1., prec));

	// Check that the boost matrix is orthogonal.
	math::Transform4 boost_inv = boost.inv();
	math::Transform4 boost_tr = boost.transpose();
	CHECK_THAT(boost_inv.t.t, RelMatcher<Real>(boost_tr.t.t, prec));
	CHECK_THAT(boost_inv.t.x, RelMatcher<Real>(boost_tr.t.x, prec));
	CHECK_THAT(boost_inv.t.y, RelMatcher<Real>(boost_tr.t.y, prec));
	CHECK_THAT(boost_inv.t.z, RelMatcher<Real>(boost_tr.t.z, prec));
	CHECK_THAT(boost_inv.x.t, RelMatcher<Real>(boost_tr.x.t, prec));
	CHECK_THAT(boost_inv.x.x, RelMatcher<Real>(boost_tr.x.x, prec));
	CHECK_THAT(boost_inv.x.y, RelMatcher<Real>(boost_tr.x.y, prec));
	CHECK_THAT(boost_inv.x.z, RelMatcher<Real>(boost_tr.x.z, prec));
	CHECK_THAT(boost_inv.y.t, RelMatcher<Real>(boost_tr.y.t, prec));
	CHECK_THAT(boost_inv.y.x, RelMatcher<Real>(boost_tr.y.x, prec));
	CHECK_THAT(boost_inv.y.y, RelMatcher<Real>(boost_tr.y.y, prec));
	CHECK_THAT(boost_inv.y.z, RelMatcher<Real>(boost_tr.y.z, prec));
	CHECK_THAT(boost_inv.z.t, RelMatcher<Real>(boost_tr.z.t, prec));
	CHECK_THAT(boost_inv.z.x, RelMatcher<Real>(boost_tr.z.x, prec));
	CHECK_THAT(boost_inv.z.y, RelMatcher<Real>(boost_tr.z.y, prec));
	CHECK_THAT(boost_inv.z.z, RelMatcher<Real>(boost_tr.z.z, prec));
}

