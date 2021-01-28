#ifndef SIDIS_CONSTANT_HPP
#define SIDIS_CONSTANT_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace constant {

static Real const PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446L;
static Real const E = 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174136L;
extern Real const INF;

static Real const ALPHA = 7.2973525664e-3L;

static Real const MASS_E = 0.0005109989461L;
static Real const MASS_MU = 0.1056583745L;
static Real const MASS_TAU = 1.77686L;

static Real const MASS_P = 0.9382720813L;
static Real const MASS_N = 0.939565413L;
static Real const MASS_D = 1.875612928L;

static Real const MASS_PI_0 = 0.1349770L;
static Real const MASS_PI = 0.13957061L;
static Real const MASS_K_0 = 0.497611L;
static Real const MASS_K = 0.493677L;
static Real const MASS_PHI = 1.019461L;
static Real const MASS_J_PSI = 3.096900L;

enum class Lepton {
	E, MU, TAU,
};

enum class Nucleus {
	P,
	N,
	D,
};

enum class Hadron {
	PI_0,
	PI_P,
	PI_M,
	K_0,
	K_P,
	K_M,
	PHI,
	J_PSI,
};

enum class Quark {
	U, D, S, C, B, T,
	U_B, D_B, S_B, C_B, B_B, T_B,
};

inline Real mass(Lepton lepton) {
	switch (lepton) {
	case Lepton::E:
		return MASS_E;
	case Lepton::MU:
		return MASS_MU;
	case Lepton::TAU:
		return MASS_TAU;
	}
	return 0.;
}

inline Real mass(Nucleus nucleus) {
	switch (nucleus) {
	case Nucleus::P:
		return MASS_P;
	case Nucleus::N:
		return MASS_N;
	case Nucleus::D:
		return MASS_D;
	}
	return 0.;
}

inline Real mass(Hadron hadron) {
	switch (hadron) {
	case Hadron::PI_0:
		return MASS_PI_0;
	case Hadron::PI_P:
	case Hadron::PI_M:
		return MASS_PI;
	case Hadron::K_0:
		return MASS_K_0;
	case Hadron::K_P:
	case Hadron::K_M:
		return MASS_K;
	case Hadron::PHI:
		return MASS_PHI;
	case Hadron::J_PSI:
		return MASS_J_PSI;
	}
	return 0.;
}

inline Real charge(Lepton lepton) {
	switch (lepton) {
	case Lepton::E:
	case Lepton::MU:
	case Lepton::TAU:
		return -1.;
	}
	return 0.;
}

inline Real charge(Nucleus nucleus) {
	switch (nucleus) {
	case Nucleus::P:
	case Nucleus::D:
		return +1.;
	case Nucleus::N:
		return 0.;
	}
	return 0.;
}

inline Real charge(Hadron hadron) {
	switch (hadron) {
	case Hadron::PI_0:
	case Hadron::K_0:
	case Hadron::PHI:
	case Hadron::J_PSI:
		return 0.;
	case Hadron::PI_P:
	case Hadron::K_P:
		return +1.;
	case Hadron::PI_M:
	case Hadron::K_M:
		return -1.;
	}
	return 0.;
}

inline Real charge(Quark quark) {
	switch (quark) {
	case Quark::U:
	case Quark::C:
	case Quark::T:
		return +2./3.;
	case Quark::D:
	case Quark::S:
	case Quark::B:
		return -1./3.;
	case Quark::U_B:
	case Quark::C_B:
	case Quark::T_B:
		return -2./3.;
	case Quark::D_B:
	case Quark::S_B:
	case Quark::B_B:
		return +1./3.;
	}
	return 0.;
}

inline char const* name(Lepton lepton) {
	switch (lepton) {
	case Lepton::E:
		return "e";
	case Lepton::MU:
		return "μ";
	case Lepton::TAU:
		return "τ";
	}
	return "<error>";
}

inline char const* name(Nucleus nucleus) {
	switch (nucleus) {
	case Nucleus::P:
		return "p";
	case Nucleus::N:
		return "n";
	case Nucleus::D:
		return "d";
	}
	return "<error>";
}

inline char const* name(Hadron hadron) {
	switch (hadron) {
	case Hadron::PI_0:
		return "π⁰";
	case Hadron::PI_P:
		return "π⁺";
	case Hadron::PI_M:
		return "π⁻";
	case Hadron::K_0:
		return "K⁰";
	case Hadron::K_P:
		return "K⁺";
	case Hadron::K_M:
		return "K⁻";
	case Hadron::PHI:
		return "φ";
	case Hadron::J_PSI:
		return "J/ψ";
	}
	return "<error>";
}

inline char const* name(Quark quark) {
	switch (quark) {
	case Quark::U:
		return "u";
	case Quark::D:
		return "d";
	case Quark::S:
		return "s";
	case Quark::C:
		return "c";
	case Quark::B:
		return "b";
	case Quark::T:
		return "t";
	case Quark::U_B:
		return "u bar";
	case Quark::D_B:
		return "d bar";
	case Quark::S_B:
		return "s bar";
	case Quark::C_B:
		return "c bar";
	case Quark::B_B:
		return "b bar";
	case Quark::T_B:
		return "t bar";
	}
	return "<error>";
}

}
}

#endif

