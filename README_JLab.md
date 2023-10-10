# SIDIS-RC EvGen on JLab ifarm

This is an event generator for polarized semi-inclusive deep inelastic
scattering (SIDIS) with QED radiative corrections, following closely work done
in [1]. The project is divided into a library component `sidis` for calculating
SIDIS cross-sections, and a binary component `sidisgen` for generating events.

A detailed description of the generator can be found on
[the arxiv](https://arxiv.org/abs/2210.03785).

## Build

This is a specific introduction for building the generator on the JLab ifarm machine. For building the generator, ROOT must be installed on your system. Optional: for building
the documentation, Doxygen must be installed; for building the tests, Catch2
must be installed.

At the beginning, check the environment.
```bash
source /path/to/thisroot.sh
module load gcc/latest
module load cmake/latest
```
Use the correct version of the software before the building process. 
```bash
# Clone the project including submodules.
git clone https://github.com/duanebyer/sidis.git
cd sidis
git submodule init
git submodule update
```
The git submodule process setup the external directory. Check the content in external directory, if it is the same as on [the github](https://github.com/duanebyer/sidis/tree/master/external). If not, download the corresponding repositories. 
```bash
# At the external directory
git clone git@github.com:duanebyer/bubble.git
git clone git@github.com:duanebyer/cubature-cpp.git
git clone git@github.com:duanebyer/mstwpdf.git
```
Then go back to the sidis directory

```bash
# Make a build directory.
mkdir build
cd build
# Configure the build.
# Optionally install the library and binary files to your system. Use
# `CMAKE_INSTALL_PREFIX` during the configure to choose the install location.
cmake ..
# Build.
make
make install
```
Now the sidisgen is installed in build/app/sidisgen. To make it more convenient to access the generator, you could store its path in a file (eg. setup_sidis). Then, when you want to use the generator, you can simply source the setup_sidis file. 
```bash
source path/to/thisroot.sh
module load gcc/latest
#Add the path to the generator
export SIDISGEN_PATH=/path/to/sidis/build/app
export PATH=${PATH}:$SIDISGEN_PATH/sidisgen
```

The following CMake configuration options may be of use:
* `Sidis_REAL_TYPE`: The floating point type to use for all cross-section
  calculations. Only the standard C++ floating point types `float`, `double`,
  and `long double` are supported.
* `Sidis_BUILD_TESTS`: Whether to build the tests.
* `Sidis_BUILD_EXAMPLES`: Whether to build the examples.
* `Sidis_BUILD_DOXYGEN`: Whether to build the Doxygen documentation.
* `Sidis_BUILD_APPS`: Whether to build the `sidisgen` binary.
* `Sidis_ENABLE_IPO`: Whether to build with interprocedural optimization.

## Quick start

### Generator

The `sidisgen` generator uses a parameter file. The allowed options are listed
by `sidisgen --help`. As an example:

```csv
mc.num_events        10000
file.gen             gen.root
file.event           event.root
phys.sf_set          prokudin
phys.rc_method       approx
phys.soft_threshold  0.01
phys.mass_threshold  1.073249081
setup.beam           e
setup.target         p
setup.hadron         pi+
setup.beam_energy    10.6
setup.beam_pol       0
setup.target_pol     0 0 0
```

Then, initialize a FOAM to approximate the differential cross-section in the
specified kinematic region.

```bash
sidisgen --initialize <param-file>
```

Finally, use the FOAM for Monte-Carlo event generation.

```bash
sidisgen --generate <param-file>
```

To process the resulting ROOT file, see `examples/process_events.cpp`.

### Job submission
The job submission system using here on JLab ifram is swif2. Here is an introduction for [JLab Scientific Computing](https://scicomp.jlab.org/scicomp/home). Check if you have a  [Slurm User Account](https://scicomp.jlab.org/scicomp/slurmJob/slurmAccount).
First, some executable files are needed 
```bash
#!/bin/bash -l
cd where/you/want/to/start
source environment/setup/file 
sidisgen --initialize path/to/input
sidisgen --generate path/to/input
```
Then, you could check if you executable file can run under interactive mode
```bash
srun --pty bash
./executable.sh
```
To submit the jobs, a json file for this job is needed. 
```javascript
{ "jobs":[{
  "account":"Your account",
  "command":["/path/to/your/executable/file_1.sh"],
  "constraint":"centos79",
  "cpu_cores":1,
  "diskBytes":20000000000,
  "name":"job_name",
  "partition":"production",
  "ramBytes":5000000000,
  "stderr":"/where/do/you/want/to/save/out_1.out",
  "stdout":"/where/do/you/want/to/save/err_1.err",
  "timeSecs":28800
  },
  {
  "account":"Your account",
  "command":["/path/to/your/executable/file_2.sh"],
  "constraint":"centos79",
  "cpu_cores":1,
  "diskBytes":20000000000,
  "name":"job_name",
  "partition":"production",
  "ramBytes":5000000000,
  "stderr":"/where/do/you/want/to/save/out_2.out",
  "stdout":"/where/do/you/want/to/save/err_2.err",
  "timeSecs":28800
    }
],
"name":"job_name"
  }
```
To submit the jobs, use swif2 to submit the json file and run the jobs
```bash
swif2 import -file jobs.json
swif2 run job_name
```
There are some commends to check the job status
```bash
#Check all your jobs
swif2 list
#Cancel one job
swif2 cancel job_name -delete
```



### Library

Many demonstrations of the different features of the `sidis` library can be
found in the `examples` folder. To get started quickly:

```cpp
#include <iostream>
#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>

sidis::Real const PI = sidis::PI;
sidis::Real const M_TH = sidis::MASS_P + sidis::MASS_PI_0;

int main() {
	sidis::part::Particles particles(
		sidis::part::Nucleus::P,   // Target nucleus.
		sidis::part::Lepton::E,    // Beam lepton.
		sidis::part::Hadron::PI_P, // Leading hadron.
		M_TH                       // Threshold mass of undetected part.
	);
	sidis::Real S = 2. * 10.6 * particles.M; // Kinematic variable `S = 2 p k1`.
	sidis::kin::PhaseSpace phase_space {
		0.2,      // Bjorken x.
		0.9,      // Bjorken y.
		0.3,      // Bjorken z.
		2.,       // Transverse momentum of hadron, squared.
		0.5 * PI, // Azimuthal angle of hadron.
		0.,       // Azimuthal angle of transverse target polarization.
	};
	sidis::kin::Kinematics kin(particles, S, phase_space);
	sidis::Real beam_pol = 0.;
	sidis::math::Vec3 target_pol(0., 0., 0.);
	// Compute structure functions with WW-type approximation.
	sidis::sf::set::ProkudinSfSet sf;
	sidis::Real born_xs = sidis::xs::born(kin, sf, beam_pol, target_pol);
	std::cout << "Born unpolarized cross-section is " << born_xs << std::endl;
	return 0;
}
```


## Acknowledgements

This software makes use of the following libraries:
* [WW-SIDIS](https://github.com/prokudin/WW-SIDIS): Structure functions for
  proton using WW-type approximation [2].
* [mstwpdf](https://mstwpdf.hepforge.org/): Parton distribution functions for
  proton [3].
* [ROOT](https://root.cern/): Plotting and data analysis.
* [FOAM](http://jadach.web.cern.ch/jadach/Foam/Index.html): Monte-Carlo event
  generation using spatial partitioning.
* [cubature](https://github.com/stevengj/cubature): Multi-dimensional numerical
  integration.
* [GSL](https://www.gnu.org/software/gsl/): Multi-dimensional Monte-Carlo
  integration.

## References

[1] I. Akushevich and A. Ilyichev, 2019. Lowest order QED radiative effects in
    polarized SIDIS. Phys. Rev. D 100(3), 033005.

[2] S. Bastami, H. Avakian, A. V. Efremov, A. Kotzinian, B. U. Musch, B.
    Parsamyan, A. Prokudin, M. Schlegel, P. Schweitzer. Semi-inclusive deep-
	inelastic scattering in Wandzura-Wilczek-type approximation. J. High Energy
	Phys. 1906(2019), 007.

[3] A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt. Parton distributions
    for the LHC. Eur. Phys. J. C. 63(2009), 189-285.

