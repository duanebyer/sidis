# SIDIS event generator

This is an event generator for polarized semi-inclusive deep inelastic
scattering (SIDIS) with QED radiative corrections, following closely work done
in [1]. Currently under construction.

## Quick start

Many demonstrations of the different features of the `sidis` library can be
found in the `examples` folder. To get started quickly:

```cpp
#include <iostream>
#include <sidis/sidis.hpp>
#include <sidis/extra/vector.hpp>

#define PI (3.14159)

int main() {
	sidis::kin::Initial initial_state(
		sidis::constant::Nucleus::P, // Target nucleus.
		sidis::constant::Lepton::E,  // Beam lepton.
		10.6                         // Beam energy.
	);
	sidis::kin::PhaseSpace phase_space {
		0.2,      // Bjorken x.
		0.9,      // Bjorken y.
		0.3,      // Bjorken z.
		2.,       // Transverse momentum of hadron, squared.
		0.5 * PI, // Azimuthal angle of hadron.
		0.,       // Azimuthal angle of transverse target polarization.
	};
	kin::Kinematics kin(
		initial_state,
		phase_space,
		sidis::constant::Hadron::PI_P,                     // Leading hadron.
		sidis::constant::MASS_P + sidis::constant::MASS_PI // Threshold mass.
	);
	sidis::kin::Final final_state(initial_state, kin);
	sidis::Real beam_pol = 0.;
	sidis::math::Vec3 target_pol(0., 0., 0.);
	// Compute structure functions with WW-type approximation.
	sidis::sf::model::WW ww;
	sidis::Real born_xs = sidis::xs::born(beam_pol, target_pol, kin, ww);
	std::cout << "Born unpolarized cross-section is " << born_xs << std::endl;
	return 0;
}
```

## Build

The `sidis` library uses CMake for building. To get started quickly, run the
following commands (on Linux):

```bash
# Clone the project including submodules.
git clone https://github.com/duanebyer/sidis.git
cd sidis
git submodule init
git submodule update
# Build the project.
mkdir build
cd build
cmake ..
make
# Optionally install the library files to your system.
make install
```

The following CMake build options may be of use:
* `Sidis_REAL_TYPE`: The floating point type to use for all cross-section
  calculations. Only the standard C++ floating point types `float`, `double`,
  and `long double` are supported.
* `BUILD_TESTING`: Whether to build the tests.
* `CMAKE_INSTALL_PREFIX`: The directory to which to install all library files.

## Acknowledgements

This software makes use of the following libraries:
* [WW-SIDIS](https://github.com/prokudin/WW-SIDIS): Structure functions for
  proton using WW-type approximation [2].
* [mstwpdf](https://mstwpdf.hepforge.org/): Parton distribution functions for
  proton [3].
* [cubature](https://github.com/stevengj/cubature): Multi-dimensional numerical
  integration.
* [cog](https://nedbatchelder.com/code/cog/): Code generation with Python.

## References

[1] I. Akushevich and A. Ilyichev, 2019. Lowest order QED radiative effects in
    polarized SIDIS. Phys. Rev. D 100(3), 033005.

[2] S. Bastami, H. Avakian, A. V. Efremov, A. Kotzinian, B. U. Musch, B.
    Parsamyan, A. Prokudin, M. Schlegel, P. Schweitzer. Semi-inclusive deep-
	inelastic scattering in Wandzura-Wilczek-type approximation. J. High Energy
	Phys. 1906(2019), 007.

[3] A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt. Parton distributions
    for the LHC. Eur. Phys. J. C. 63(2009), 189-285.

