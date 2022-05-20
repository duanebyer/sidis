#ifndef SIDIS_INTEG_PARAMS_HPP
#define SIDIS_INTEG_PARAMS_HPP

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
 * Available methods to use for numerical integration.
 */
enum class IntegMethod {
	CUBATURE,
	MC_PLAIN,
	MISER,
	VEGAS,
};

/**
 * Parameters to use for numerical integration.
 *
 * \sa xs::nrad_integ()
 * \sa xs::rad_f_integ()
 * \sa xs::rad_integ()
 */
struct IntegParams {
	/// Method of integration to use.
	IntegMethod method;
	/// Number of evaluations of the integrand.
	unsigned num_evals;
	/// Relative precision goal for the integral.
	Real err_rel;
	/// Absolute precision goal for the integral.
	Real err_abs;
};

}
}

#endif

