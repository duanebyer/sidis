#ifndef SIDISGEN_PARAMS_HPP
#define SIDISGEN_PARAMS_HPP

#include <istream>
#include <ostream>
#include <string>

#include <TFile.h>

#include <sidis/numeric.hpp>
#include <sidis/constant.hpp>
#include <sidis/extra/vector.hpp>

// Keeps track of the various parameters that can be used for a run of the
// generator. This structure can be read/written to both a ROOT file or a plain
// text file.
struct Params {
	std::string event_file_name;
	std::string foam_file_name;
	sidis::Real beam_energy;
	sidis::constant::Nucleus target;
	sidis::constant::Lepton beam;
	sidis::math::Vec3 target_pol;
	sidis::Real beam_pol;
	sidis::Real k0_cut;

	Params() :
		event_file_name("gen.root"),
		foam_file_name("foam.root"),
		beam_energy(10.),
		target(sidis::constant::Nucleus::P),
		beam(sidis::constant::Lepton::E),
		target_pol(0., 0., 0.),
		beam_pol(0.),
		k0_cut(0.01) { }

	void write_root(TFile& file) const;
	void read_root(TFile& file);

	void write(std::ostream& file) const;
	void read(std::istream& file);
};

#endif

