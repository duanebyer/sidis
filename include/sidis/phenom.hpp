#ifndef SIDIS_PHENOM_HPP
#define SIDIS_PHENOM_HPP

#include "sidis/numeric.hpp"

namespace sidis {

namespace kin {
	struct Kinematics;
}

namespace ph {

/**
 * \defgroup PhenomGroup Phenomenlogical inputs
 * Types and functions used to calculate various phenomenological inputs needed
 * for radiative corrections to the cross-section.
 */
/// \{

/// QED coupling constant at electron mass.
Real const ALPHA_QED_0 = 1.L / 137.03599L;
/// Calculates the QED coupling constant at a certain \f$Q^2\f$ value.
Real alpha_qed(Real Q_sq=0.);
/// Factor \f$\frac{\alpha}{\pi}\delta_{\text{vac}}^{\text{had}}\sigma_{B}\f$
/// gives the vacuum polarization cross-section due to hadron loops.
Real delta_vac_had(kin::Kinematics const& kin);

/**
 * Bundle of phenomenological inputs for the cross section, at a specific
 * kinematics. This structure caches the phenomenological inputs so they can be
 * re-used throughout various cross-section calculations. It also allows for
 * custom phenomenological inputs to be provided, e.g., an improved hadron
 * vacuum polarization calculation.
 */
struct Phenom {
	/// Electromagnetic coupling constant, evolved to the scale of the process.
	Real alpha_qed;
	/// Contribution to vacuum polarization from hadrons.
	Real delta_vac_had;

	/// Computes the phenomenological inputs at the provided kinematics.
	explicit Phenom(kin::Kinematics const& kin);
	/// Choose the electromagnetic coupling constant, but compute the other
	/// phenomenological inputs with the default method at the provided
	/// kinematics.
	Phenom(Real alpha_qed, kin::Kinematics const& kin);
	/// Provide custom values for the phenomenological inputs.
	Phenom(Real alpha_qed, Real delta_vac_had);
};
/// \}

}
}

#endif

