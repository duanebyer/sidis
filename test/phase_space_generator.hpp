#ifndef SIDIS_TEST_PHASE_SPACE_GENERATOR_HPP
#define SIDIS_TEST_PHASE_SPACE_GENERATOR_HPP

#include <catch2/catch.hpp>

#include <memory>

#include <sidis/sidis.hpp>

/*
 * Generates random points in phase space.
 */
class PhaseSpaceGenerator final :
		public Catch::Generators::IGenerator<sidis::kin::PhaseSpace> {
	sidis::part::Particles _ps;
	sidis::Real _S;
	sidis::kin::PhaseSpace _value;

public:
	PhaseSpaceGenerator(sidis::part::Particles ps, sidis::Real S) :
			_ps(ps),
			_S(S) {
		next();
	}

	bool next() override;
	sidis::kin::PhaseSpace const& get() const override {
		return _value;
	}
};

/*
 * Generates random points in radiative phase space.
 */
class PhaseSpaceRadGenerator final :
		public Catch::Generators::IGenerator<sidis::kin::PhaseSpaceRad> {
	sidis::part::Particles _ps;
	sidis::Real _S;
	sidis::kin::PhaseSpaceRad _value;

public:
	PhaseSpaceRadGenerator(sidis::part::Particles ps, sidis::Real S) :
			_ps(ps),
			_S(S) {
		next();
	}

	bool next() override;
	sidis::kin::PhaseSpaceRad const& get() const override {
		return _value;
	}
};

/*
 * Generates random points near the surface of phase space.
 */
class PhaseSpaceSurfaceGenerator final :
		public Catch::Generators::IGenerator<sidis::kin::PhaseSpace> {
	sidis::part::Particles _ps;
	sidis::Real _S;
	sidis::Real _skin_depth;
	sidis::kin::PhaseSpace _value;

public:
	PhaseSpaceSurfaceGenerator(
			sidis::part::Particles ps,
			sidis::Real S,
			sidis::Real skin_depth) :
			_ps(ps),
			_S(S),
			_skin_depth(skin_depth) {
		next();
	}

	bool next() override;
	sidis::kin::PhaseSpace const& get() const override {
		return _value;
	}
};

/*
 * Generates random points near the surface of radiative phase space.
 */
class PhaseSpaceRadSurfaceGenerator final :
		public Catch::Generators::IGenerator<sidis::kin::PhaseSpaceRad> {
	sidis::part::Particles _ps;
	sidis::Real _S;
	sidis::Real _skin_depth;
	sidis::kin::PhaseSpaceRad _value;

public:
	PhaseSpaceRadSurfaceGenerator(
			sidis::part::Particles ps,
			sidis::Real S,
			sidis::Real skin_depth) :
			_ps(ps),
			_S(S),
			_skin_depth(skin_depth) {
		next();
	}

	bool next() override;
	sidis::kin::PhaseSpaceRad const& get() const override {
		return _value;
	}
};

Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpace> gen_phase_space(
		sidis::part::Particles ps, sidis::Real S);
Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpaceRad> gen_phase_space_rad(
		sidis::part::Particles ps, sidis::Real S);
Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpace> gen_phase_space_surface(
		sidis::part::Particles ps, sidis::Real S, sidis::Real skin_depth);
Catch::Generators::GeneratorWrapper<sidis::kin::PhaseSpaceRad> gen_phase_space_rad_surface(
		sidis::part::Particles ps, sidis::Real S, sidis::Real skin_depth);

#endif

