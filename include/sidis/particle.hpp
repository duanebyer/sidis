#ifndef SIDIS_PARTICLE_HPP
#define SIDIS_PARTICLE_HPP

#include "sidis/constant.hpp"

namespace sidis {
namespace part {

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

}
}

#endif

