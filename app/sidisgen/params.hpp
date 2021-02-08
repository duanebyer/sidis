#ifndef SIDISGEN_PARAMS_HPP
#define SIDISGEN_PARAMS_HPP

#include <istream>
#include <ostream>
#include <string>

#include <TFile.h>

#include <sidis/numeric.hpp>
#include <sidis/constant.hpp>
#include <sidis/extra/bounds.hpp>
#include <sidis/extra/vector.hpp>

// Within a major version, there is forward compatibility (e.x. 1.4 is
// forward-compatible with 1.5). Between major versions, there is no
// compatibility.
#define SIDIS_PARAMS_VERSION_MAJOR 1
#define SIDIS_PARAMS_VERSION_MINOR 2

// Keeps track of the various parameters that can be used for a run of the
// generator. This structure can be read/written to both a ROOT file or a plain
// text file.
struct Params {
	Int_t version_major;
	Int_t version_minor;
	std::string event_file;
	std::string foam_nrad_file;
	std::string foam_rad_file;
	std::string sf_set;
	Long_t num_events;
	Long_t num_init;
	Int_t seed;
	Int_t seed_init;
	sidis::Real beam_energy;
	sidis::constant::Lepton beam;
	sidis::constant::Nucleus target;
	sidis::constant::Hadron hadron;
	sidis::Real mass_threshold;
	sidis::math::Vec3 target_pol;
	sidis::Real beam_pol;
	sidis::Real k_0_bar;
	sidis::math::Bounds x_cut;
	sidis::math::Bounds y_cut;
	sidis::math::Bounds z_cut;
	sidis::math::Bounds ph_t_sq_cut;
	sidis::math::Bounds phi_h_cut;
	sidis::math::Bounds phi_cut;

	Params() :
		version_major(SIDIS_PARAMS_VERSION_MAJOR),
		version_minor(SIDIS_PARAMS_VERSION_MINOR),
		event_file("gen.root"),
		foam_nrad_file("foam-nrad.root"),
		foam_rad_file("foam-rad.root"),
		sf_set("ww"),
		num_events(10000),
		num_init(1000),
		seed(0),
		seed_init(0),
		beam_energy(10.),
		beam(sidis::constant::Lepton::E),
		target(sidis::constant::Nucleus::P),
		hadron(sidis::constant::Hadron::PI_P),
		mass_threshold(sidis::constant::MASS_P + sidis::constant::MASS_PI_0),
		target_pol(0., 0., 0.),
		beam_pol(0.),
		k_0_bar(0.01),
		x_cut(sidis::math::Bounds::UNIT),
		y_cut(sidis::math::Bounds::UNIT),
		z_cut(sidis::math::Bounds::UNIT),
		ph_t_sq_cut(sidis::math::Bounds::POSITIVE),
		phi_h_cut(-sidis::constant::PI, sidis::constant::PI),
		phi_cut(-sidis::constant::PI, sidis::constant::PI) { }

	void write_root(TFile& file) const;
	void read_root(TFile& file);

	void write(std::ostream& file) const;
	void read(std::istream& file);

	// Check whether a `TFoam` generated for one set of parameters can be used
	// for another set.
	bool compatible_foam(Params const& foam_params) const;
};

#endif

