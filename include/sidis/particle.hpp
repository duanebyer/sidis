#ifndef SIDIS_PARTICLE_HPP
#define SIDIS_PARTICLE_HPP

#include "sidis/constant.hpp"

namespace sidis {
namespace part {

/**
 * \defgroup ParticleInfo Particle info
 * Masses and charges of various types of particles.
 */
/// \{

/// Beam leptons types.
enum class Lepton {
	E,          ///< Electron \f$e^{-}\f$
	MU,         ///< Muon \f$\mu^{-}\f$
	TAU,        ///< Tau lepton \f$\tau^{-}\f$
	MASSLESS_E, ///< "Massless" electron \f$e^{-}\f$
};

/// Target nuclei types.
enum class Nucleus {
	P, ///< Proton \f$p\f$
	N, ///< Neutron \f$n\f$
	D, ///< Deuteron \f$d\f$
};

/// Detected hadron types.
enum class Hadron {
	PI_0,  ///< Pion \f$\pi^{0}\f$
	PI_P,  ///< Pion \f$\pi^{+}\f$
	PI_M,  ///< Pion \f$\pi^{-}\f$
	K_0,   ///< Kaon \f$K^{0}\f$
	K_P,   ///< Kaon \f$K^{+}\f$
	K_M,   ///< Kaon \f$K^{-}\f$
	PHI,   ///< Phi meson \f$\phi\f$
	J_PSI, ///< J/psi meson \f$J/\psi\f$
};

/// A full description of the particles involved in a SIDIS process. Contains a
/// target Nucleus, a Lepton beam, and a leading Hadron.
struct Particles {
	/// Target nucleus type.
	part::Nucleus const target;
	/// Lepton beam type.
	part::Lepton const beam;
	/// Detected hadron type.
	part::Hadron const hadron;
	/// Mass of target.
	Real const M;
	/// Mass of beam lepton.
	Real const m;
	/// Mass of detected hadron.
	Real const mh;
	/// Threshold mass for SIDIS process.
	Real const Mth;

	/// Bundle together a set of particles describing a setup for a SIDIS
	/// experiment.
	Particles(
		part::Nucleus target,
		part::Lepton beam,
		part::Hadron hadron,
		Real Mth);
};

/// \name Particle masses
/// \{

/// Return the mass of a particle.
inline Real mass(Lepton lepton) {
	switch (lepton) {
	case Lepton::E:
		return MASS_E;
	case Lepton::MU:
		return MASS_MU;
	case Lepton::TAU:
		return MASS_TAU;
	case Lepton::MASSLESS_E:
		return 0.;
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
/// \}

/// \name Particle charges
/// \{

/// Returns the charge of a particle.
inline Real charge(Lepton lepton) {
	switch (lepton) {
	case Lepton::E:
	case Lepton::MU:
	case Lepton::TAU:
	case Lepton::MASSLESS_E:
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
/// \}

/// \name Particle names
/// \{

/// Returns a string naming the particle.
inline char const* name(Lepton lepton) {
	switch (lepton) {
	case Lepton::E:
		return "e";
	case Lepton::MU:
		return "μ";
	case Lepton::TAU:
		return "τ";
	case Lepton::MASSLESS_E:
		return "e";
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
/// \}

/// \}

}
}

#endif

