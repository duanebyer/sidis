#ifndef SIDIS_CONSTANT_HPP
#define SIDIS_CONSTANT_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace constant {

static Real const PI = 3.141592653589793238462643383279502797479068098137295573004504331874296718662975536062731407582759857177734375;
static Real const E = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466392;

static Real const ALPHA = 7.2973525664e-3;

static Real const MASS_E = 0.0005109989461;
static Real const MASS_MU = 0.1056583745;
static Real const MASS_TAU = 1.77686;

static Real const MASS_P = 0.9382720813;
static Real const MASS_N = 0.939565413;
static Real const MASS_D = 1.875612928;

static Real const MASS_PI_0 = 0.1349770;
static Real const MASS_PI = 0.13957061;
static Real const MASS_K_0 = 0.497611;
static Real const MASS_K = 0.493677;

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

inline Real mass(Nucleus nucleon) {
	switch (nucleon) {
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

inline Real charge(Nucleus nucleon) {
	switch (nucleon) {
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

}
}

#endif

