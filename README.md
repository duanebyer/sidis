# SIDIS event generator

This is an event generator for polarized semi-inclusive deep inelastic
scattering (SIDIS) with QED radiative corrections, following closely work done
in [1]. Currently under construction.

## Build

The `sidis` library uses CMake for building. To get started quickly, run the
following commands (using bash):

```
# First, clone the project including submodules.
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
* [WW-SIDIS](https://github.com/prokudin/WW-SIDIS) structure functions for
  proton using WW-type approximation [2].
* [mstwpdf](https://mstwpdf.hepforge.org/) parton distribution functions for
  proton [3].

## References

[1] I. Akushevich and A. Ilyichev, 2019. Lowest order QED radiative effects in
    polarized SIDIS. Phys. Rev. D 100(3), 033005.

[2] S. Bastami, H. Avakian, A. V. Efremov, A. Kotzinian, B. U. Musch, B.
    Parsamyan, A. Prokudin, M. Schlegel, P. Schweitzer. Semi-inclusive deep-
	inelastic scattering in Wandzura-Wilczek-type approximation. J. High Energy
	Phys. 1906(2019), 007.

[3] A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt. Parton distributions
    for the LHC. Eur. Phys. J. C. 63(2009), 189-285.

