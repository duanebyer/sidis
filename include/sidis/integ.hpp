#ifndef SIDIS_INTEG_HPP
#define SIDIS_INTEG_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

/**
 * An integral estimate combined with a symmetric error bar.
 *
 * \sa xs::nrad_integ()
 * \sa xs::rad_f_integ()
 * \sa xs::rad_integ()
 */
struct EstErr {
	Real val;
	Real err;
};

/**
 * Parameters to be used for numerical integration.
 *
 * \sa xs::nrad_integ()
 * \sa xs::rad_f_integ()
 * \sa xs::rad_integ()
 */
struct IntegParams {
	/// Maximum number of evaluations of the integrand.
	unsigned max_evals;
	/// Relative precision goal for the integral.
	Real err_rel;
	/// Absolute precision goal for the integral.
	Real err_abs;
};

}
}

#endif

