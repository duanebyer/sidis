#ifndef SIDISGEN_PARAMS_HPP
#define SIDISGEN_PARAMS_HPP

#include <istream>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>

#include <TFile.h>

#include <sidis/sidis.hpp>

#include "event_type.hpp"

// Within a major version, there is forward compatibility (e.x. 1.4 is
// forward-compatible with 1.5). Between major versions, there is no
// compatibility.
#define SIDIS_PARAMS_VERSION_MAJOR 4
#define SIDIS_PARAMS_VERSION_MINOR 1

struct Toggle {
	bool on;
	Toggle() : on(false) { }
	Toggle(bool on) : on(on) { }
	operator bool() const {
		return on;
	}
};

enum class RcMethod {
	NONE,
	APPROX,
	EXACT,
};

struct Version {
	int v_major;
	int v_minor;
	Version(
		int v_major=SIDIS_PARAMS_VERSION_MAJOR,
		int v_minor=SIDIS_PARAMS_VERSION_MINOR) :
		v_major(v_major),
		v_minor(v_minor) { }
	bool operator==(Version const& rhs) const {
		return v_major == rhs.v_major && v_minor == rhs.v_minor;
	}
	bool operator!=(Version const& rhs) const {
		return !(*this == rhs);
	}
};

// Represents a single parameter for the generator. Stores a name, and an
// optional value.
template<typename T>
class Param {
	std::string _name;
	T _value;
	bool _occupied;

public:
	Param(std::string name, T const& value) :
		_name(name),
		_value(value),
		_occupied(true) { }
	Param(std::string name) :
		_name(name),
		_value(),
		_occupied(false) { }

	T const& operator*() const {
		if (!_occupied) {
			throw std::runtime_error(
				"No value to get from parameter '" + _name + "'.");
		} else {
			return _value;
		}
	}
	T& operator*() {
		if (!_occupied) {
			throw std::runtime_error(
				"No value to get from parameter '" + _name + "'.");
		} else {
			return _value;
		}
	}
	T const* operator->() const {
		return &(operator*());
	}
	T* operator->() {
		return &(operator*());
	}
	T const& get() const {
		return operator*();
	}
	T& get() {
		return operator*();
	}
	T const& get_or(T const& default_value) const {
		if (_occupied) {
			return _value;
		} else {
			return default_value;
		}
	}
	T& get_or_insert(T const& default_value) {
		if (!_occupied) {
			_value = default_value;
			_occupied = true;
		}
		return _value;
	}
	void take_first(Param<T> const& other) {
		if (!_occupied) {
			if (other._occupied) {
				_value = other._value;
				_occupied = true;
			}
		}
	}
	void reset() {
		_value = T();
		_occupied = false;
	}
	void reset(T const& value) {
		_value = value;
		_occupied = true;
	}

	char const* name() const {
		return _name.c_str();
	}
	bool occupied() const {
		return _occupied;
	}
	explicit operator bool() const {
		return _occupied;
	}

	bool operator==(Param<T> const& rhs) const {
		if (_occupied && rhs._occupied) {
			return _value == rhs._value;
		} else {
			return _occupied == rhs._occupied;
		}
	}
	bool operator!=(Param<T> const& rhs) const {
		return !(*this == rhs);
	}
};

// Keeps track of the various parameters that can be used for a run of the
// generator. This structure can be read/written to both a ROOT file or a plain
// text file.
struct Params {
	Param<Version> version;
	Param<Toggle> strict;
	Param<std::string> event_file;
	Param<RcMethod> rc_method;
	Param<Toggle> nrad_gen;
	Param<Double_t> nrad_rej_scale;
	Param<Int_t> nrad_max_cells;
	Param<Double_t> nrad_target_eff;
	Param<Double_t> nrad_scale_exp;
	Param<Toggle> rad_gen;
	Param<Double_t> rad_rej_scale;
	Param<Int_t> rad_max_cells;
	Param<Double_t> rad_target_eff;
	Param<Double_t> rad_scale_exp;
	Param<Toggle> write_momenta;
	Param<Toggle> write_photon;
	Param<Toggle> write_sf_set;
	Param<Toggle> write_mc_coords;
	Param<std::string> foam_file;
	Param<std::string> sf_set;
	Param<Long_t> num_events;
	Param<Double_t> rej_weight;
	Param<std::multiset<Int_t> > seed;
	Param<Int_t> nrad_seed_init;
	Param<Int_t> rad_seed_init;
	Param<sidis::Real> beam_energy;
	Param<sidis::part::Lepton> beam;
	Param<sidis::part::Nucleus> target;
	Param<sidis::part::Hadron> hadron;
	Param<sidis::Real> Mth;
	Param<sidis::math::Vec3> target_pol;
	Param<sidis::Real> beam_pol;
	Param<sidis::Real> k_0_bar;
	Param<sidis::math::Bound> cut_x;
	Param<sidis::math::Bound> cut_y;
	Param<sidis::math::Bound> cut_z;
	Param<sidis::math::Bound> cut_ph_t_sq;
	Param<sidis::math::Bound> cut_phi_h;
	Param<sidis::math::Bound> cut_phi;
	Param<sidis::math::Bound> cut_Q_sq;
	Param<sidis::math::Bound> cut_t;
	Param<sidis::math::Bound> cut_W_sq;
	Param<sidis::math::Bound> cut_r;
	Param<sidis::math::Bound> cut_mx_sq;
	Param<sidis::math::Bound> cut_qt_to_Q;
	Param<sidis::math::Bound> cut_lab_mom_q;
	Param<sidis::math::Bound> cut_lab_mom_k2;
	Param<sidis::math::Bound> cut_lab_mom_h;
	Param<sidis::math::Bound> cut_lab_theta_q;
	Param<sidis::math::Bound> cut_lab_theta_k2;
	Param<sidis::math::Bound> cut_lab_theta_h;
	Param<sidis::math::Bound> cut_tau;
	Param<sidis::math::Bound> cut_phi_k;
	Param<sidis::math::Bound> cut_R;
	Param<sidis::math::Bound> cut_k_0_bar;
	Param<sidis::math::Bound> cut_lab_mom_k;
	Param<sidis::math::Bound> cut_lab_theta_k;

	// TODO: Re-order these to be compatible with the help command.
	Params() :
		version("version"),
		strict("strict"),
		event_file("file.event_out"),
		rc_method("phys.rc_method"),
		nrad_gen("mc.nrad.gen"),
		nrad_rej_scale("mc.nrad.gen.rej_scale"),
		nrad_max_cells("mc.nrad.init.max_cells"),
		nrad_target_eff("mc.nrad.init.target_eff"),
		nrad_scale_exp("mc.nrad.init.scale_exp"),
		rad_gen("mc.rad.gen"),
		rad_rej_scale("mc.rad.gen.rej_scale"),
		rad_max_cells("mc.rad.init.max_cells"),
		rad_target_eff("mc.rad.init.target_eff"),
		rad_scale_exp("mc.rad.init.scale_exp"),
		write_momenta("file.write_momenta"),
		write_photon("file.write_photon"),
		write_sf_set("file.write_sf_set"),
		write_mc_coords("file.write_mc_coords"),
		foam_file("file.foam_out"),
		sf_set("phys.sf_set"),
		num_events("mc.num_events"),
		rej_weight("mc.rej_weight"),
		seed("mc.seed"),
		nrad_seed_init("mc.nrad.init.seed"),
		rad_seed_init("mc.rad.init.seed"),
		beam_energy("setup.beam_energy"),
		beam("setup.beam"),
		target("setup.target"),
		hadron("setup.hadron"),
		Mth("phys.mass_threshold"),
		target_pol("setup.target_pol"),
		beam_pol("setup.beam_pol"),
		k_0_bar("phys.soft_threshold"),
		cut_x("cut.x"),
		cut_y("cut.y"),
		cut_z("cut.z"),
		cut_ph_t_sq("cut.ph_t_sq"),
		cut_phi_h("cut.phi_h"),
		cut_phi("cut.phi"),
		cut_Q_sq("cut.Q_sq"),
		cut_t("cut.t"),
		cut_W_sq("cut.W_sq"),
		cut_r("cut.r"),
		cut_mx_sq("cut.mx_sq"),
		cut_qt_to_Q("cut.qt_to_Q"),
		cut_lab_mom_q("cut.lab_mom_q"),
		cut_lab_mom_k2("cut.lab_mom_k2"),
		cut_lab_mom_h("cut.lab_mom_h"),
		cut_lab_theta_q("cut.lab_theta_q"),
		cut_lab_theta_k2("cut.lab_theta_k2"),
		cut_lab_theta_h("cut.lab_theta_h"),
		cut_tau("cut.tau"),
		cut_phi_k("cut.phi_k"),
		cut_R("cut.R"),
		cut_k_0_bar("cut.k_0_bar"),
		cut_lab_mom_k("cut.lab_mom_k"),
		cut_lab_theta_k("cut.lab_theta_k") { }

	void write_root(TFile& file) const;
	void read_root(TFile& file);

	void write_stream(std::ostream& file) const;
	void read_stream(std::istream& file);

	// Takes the supplied parameters and fills in missing ones to make a valid
	// run card. If unable to do so, throw an exception. If `strict` is enabled,
	// then `fill_defaults` will never change a parameter that has been set by
	// the user (for example, disabling `write_photon` when no radiative
	// corrections are being applied).
	void fill_defaults();
	bool valid() const {
		Params other = *this;
		other.fill_defaults();
		return other == *this;
	}

	// Check whether a FOAM produced using a set of parameters can be used to
	// generate events with this set of parameters. If not compatible, throws an
	// exception with a message describing the incompatibility.
	void compatible_with_foam(EventType type, Params const& foam_params) const;
	// Checks whether two parameter files used for event generation can be
	// merged together.
	void compatible_with_merge(Params const& other_params) const;
	// Combine two parameter files together. If the parameter files cannot be
	// combined, then throws an exception.
	void merge(Params const& params);

	bool operator==(Params const& rhs) const;
	bool operator!=(Params const& rhs) const {
		return !(*this == rhs);
	}
};

#endif

