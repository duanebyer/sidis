#include "phase_space_generator.hpp"

#include <random>

#include <sidis/bound.hpp>
#include <sidis/extra/math.hpp>

using namespace sidis;

namespace {

std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<Real> dist(0., 1.);

}

bool PhaseSpaceGenerator::next() {
	Real point[6] = {
		dist(rng), dist(rng), dist(rng),
		dist(rng), dist(rng), dist(rng),
	};
	cut::take(_ps, _S, point, &_value, nullptr);
	return true;
}

bool PhaseSpaceRadGenerator::next() {
	Real point[9] = {
		dist(rng), dist(rng), dist(rng),
		dist(rng), dist(rng), dist(rng),
		dist(rng), dist(rng), dist(rng),
	};
	cut::take(_ps, _S, point, &_value, nullptr);
	return true;
}

bool PhaseSpaceSurfaceGenerator::next() {
	// Choose a random variable from `x`, `y`, `z`, `ph_t_sq` to put near the
	// surface.
	Real select_float = dist(rng);
	int select;
	if (select_float < 1./4.) {
		select = 0;
	} else if (select_float < 2./4.) {
		select = 1;
	} else if (select_float < 3./4.) {
		select = 2;
	} else {
		select = 3;
	}
	Real value = dist(rng) < 0.5 ? 1. - _skin_depth : _skin_depth;
	Real x = cut::x_bound(_ps, _S).lerp(select == 0 ? value : dist(rng));
	Real y = cut::y_bound(_ps, _S, x).lerp(select == 1 ? value : dist(rng));
	Real z = cut::z_bound(_ps, _S, x, y).lerp(select == 2 ? value : dist(rng));
	Real ph_t_sq = cut::ph_t_sq_bound(_ps, _S, x, y, z).lerp(select == 3 ? value : dist(rng));
	Real phi_h = math::Bound(-PI, PI).lerp(dist(rng));
	Real phi = math::Bound(-PI, PI).lerp(dist(rng));
	kin::PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, phi };
	_value = ph_space;
	return true;
}

bool PhaseSpaceRadSurfaceGenerator::next() {
	// Choose a random variable from `x`, `y`, `z`, `ph_t_sq`, `tau`, `R` to put
	// near the surface.
	Real select_float = dist(rng);
	int select;
	if (select_float < 1./6.) {
		select = 0;
	} else if (select_float < 2./6.) {
		select = 1;
	} else if (select_float < 3./6.) {
		select = 2;
	} else if (select_float < 4./6.) {
		select = 3;
	} else if (select_float < 5./6.) {
		select = 4;
	} else {
		select = 5;
	}
	Real value = dist(rng) < 0.5 ? 1. - _skin_depth : _skin_depth;
	Real x = cut::x_bound(_ps, _S).lerp(select == 0 ? value : dist(rng));
	Real y = cut::y_bound(_ps, _S, x).lerp(select == 1 ? value : dist(rng));
	Real z = cut::z_bound(_ps, _S, x, y).lerp(select == 2 ? value : dist(rng));
	Real ph_t_sq = cut::ph_t_sq_bound(_ps, _S, x, y, z).lerp(select == 3 ? value : dist(rng));
	Real phi_h = math::Bound(-PI, PI).lerp(dist(rng));
	Real phi = math::Bound(-PI, PI).lerp(dist(rng));
	kin::PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, phi };
	kin::Kinematics kin(_ps, _S, ph_space);
	Real tau = cut::tau_bound(kin).lerp(select == 4 ? value : dist(rng));
	Real phi_k = math::Bound(-PI, PI).lerp(dist(rng));
	Real R = cut::R_bound(kin, tau, phi_k).lerp(select == 5 ? value : dist(rng));
	kin::PhaseSpaceRad ph_space_rad {
		x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R,
	};
	_value = ph_space_rad;
	return true;
}

Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpace> gen_phase_space(
		part::Particles ps, Real S) {
	return Catch::Generators::GeneratorWrapper<kin::PhaseSpace>(
		std::unique_ptr<Catch::Generators::IGenerator<kin::PhaseSpace> >(
			new PhaseSpaceGenerator(ps, S)));
}

Catch::Generators::GeneratorWrapper<kin::PhaseSpaceRad> gen_phase_space_rad(
		part::Particles ps, Real S) {
	return Catch::Generators::GeneratorWrapper<kin::PhaseSpaceRad>(
		std::unique_ptr<Catch::Generators::IGenerator<kin::PhaseSpaceRad> >(
			new PhaseSpaceRadGenerator(ps, S)));
}

Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpace> gen_phase_space_surface(
		part::Particles ps, Real S, Real skin_depth) {
	return Catch::Generators::GeneratorWrapper<kin::PhaseSpace>(
		std::unique_ptr<Catch::Generators::IGenerator<kin::PhaseSpace> >(
			new PhaseSpaceSurfaceGenerator(ps, S, skin_depth)));
}

Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpaceRad> gen_phase_space_rad_surface(
		part::Particles ps, Real S, Real skin_depth) {
	return Catch::Generators::GeneratorWrapper<kin::PhaseSpaceRad>(
		std::unique_ptr<Catch::Generators::IGenerator<kin::PhaseSpaceRad> >(
			new PhaseSpaceRadSurfaceGenerator(ps, S, skin_depth)));
}

