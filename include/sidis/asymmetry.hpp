#ifndef SIDIS_ASYMMETRY_HPP
#define SIDIS_ASYMMETRY_HPP

#include "sidis/integ_params.hpp"
#include "sidis/numeric.hpp"

namespace sidis {

namespace part {
	struct Particles;
}
namespace sf {
	class SfSet;
}

namespace asym {

math::IntegParams const DEFAULT_INTEG_PARAMS {
	math::IntegMethod::CUBATURE,
	100000,
	10000,
};
math::IntegParams const DEFAULT_INTEG_PARAMS_ASYM {
	math::IntegMethod::CUBATURE,
	100000,
	10000,
};

/**
 * \defgroup AsymGroup Transverse single-spin asymmetries
 * Functions used to calculate transverse single-spin asymmetries.
 *
 * These methods make use of the approximation \f$\phi_S\approx\phi\f$, which
 * may not be valid at low energies.
 */

/// Integral of unpolarized cross-section over \f$\phi_h}\f$ and \f$\phi\f$.
math::EstErr uu_integ(
	sf::SfSet const& sf_set,
	part::Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
	bool include_rc=false,
	math::IntegParams params=DEFAULT_INTEG_PARAMS);

/// Integral of transverse part of cross-section \f$\sigma^{UT}-\sigma^{UU}\f$.
///
/// The integral is weighted by a \f$\sin\f$ function with integer coefficients
/// for \f$\phi_S\f$ and \f$\phi_h\f$, and an offset in quarter-circles.
math::EstErr ut_integ(
	sf::SfSet const& sf_set,
	part::Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
	int phi_s_coeff, int phi_h_coeff, int offset,
	bool include_rc=false,
	math::IntegParams params=DEFAULT_INTEG_PARAMS);

/// Transverse single-spin asymmetry.
///
/// The asymmetry is chosen by a \f$\sin\f$ function with integer coefficients
/// for \f$\phi_S\f$ and \f$\phi_h\f$, and an offset in quarter-circles.
math::EstErr ut_asym(
	sf::SfSet const& sf_set,
	part::Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
	int phi_s_coeff, int phi_h_coeff, int offset,
	bool include_rc=false,
	math::IntegParams params=DEFAULT_INTEG_PARAMS_ASYM);

/// Sivers asymmetry.
math::EstErr sivers(
	sf::SfSet const& sf_set,
	part::Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
	bool include_rc=false,
	math::IntegParams params=DEFAULT_INTEG_PARAMS_ASYM);

/// Collins asymmetry.
math::EstErr collins(
	sf::SfSet const& sf_set,
	part::Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
	bool include_rc=false,
	math::IntegParams params=DEFAULT_INTEG_PARAMS_ASYM);

}
}

#endif

