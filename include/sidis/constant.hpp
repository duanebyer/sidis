#ifndef SIDIS_CONSTANT_HPP
#define SIDIS_CONSTANT_HPP

#include "sidis/numeric.hpp"

namespace sidis {

/**
 * \defgroup Constant Constants
 * \addtogroup Constant
 */
/// \{

/// Mathematical constant \f$\pi\f$.
static Real const PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446L;
/// Mathematical constant \f$e\f$.
static Real const E = 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174136L;
/// Positive infinity.
extern Real const INF;

/// Fine structure constant.
static Real const ALPHA = 7.2973525664e-3L;

/// \name Particle masses
/// \{
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
/// \}

/// \}

}

#endif

