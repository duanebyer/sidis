#ifndef SIDISGEN_PARAMS_HPP
#define SIDISGEN_PARAMS_HPP

#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

#include <TFile.h>

#include <sidis/sidis.hpp>

// Within a major version, there is forward compatibility (e.x. 1.4 is
// forward-compatible with 1.5). Between major versions, there is no
// compatibility.
#define SIDIS_PARAMS_VERSION_MAJOR 2
#define SIDIS_PARAMS_VERSION_MINOR 0

enum class RcMethod {
	NONE,
	APPROX,
	EXACT,
};

struct Version {
	int major;
	int minor;
	Version(
		int major=SIDIS_PARAMS_VERSION_MAJOR,
		int minor=SIDIS_PARAMS_VERSION_MINOR) :
		major(major),
		minor(minor) { }
	bool operator==(Version const& rhs) const {
		return major == rhs.major && minor == rhs.minor;
	}
	bool operator!=(Version const& rhs) const {
		return !(*this == rhs);
	}
};

// Represents a single parameter for the generator. Stores a name, and an
// optional value.
template<typename T>
class Param {
	std::string const _name;
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
	operator bool() const {
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
	Param<std::string> event_file;
	Param<RcMethod> rc_method;
	Param<bool> gen_nrad;
	Param<bool> gen_rad;
	Param<bool> write_photon;
	Param<std::string> foam_nrad_file;
	Param<std::string> foam_rad_file;
	Param<std::string> sf_set;
	Param<Long_t> num_events;
	Param<Long_t> num_init;
	Param<Int_t> seed;
	Param<Int_t> seed_init;
	Param<sidis::Real> beam_energy;
	Param<sidis::part::Lepton> beam;
	Param<sidis::part::Nucleus> target;
	Param<sidis::part::Hadron> hadron;
	Param<sidis::Real> mass_threshold;
	Param<sidis::math::Vec3> target_pol;
	Param<sidis::Real> beam_pol;
	Param<sidis::Real> k_0_bar;
	Param<sidis::math::Bound> x_cut;
	Param<sidis::math::Bound> y_cut;
	Param<sidis::math::Bound> z_cut;
	Param<sidis::math::Bound> ph_t_sq_cut;
	Param<sidis::math::Bound> phi_h_cut;
	Param<sidis::math::Bound> phi_cut;
	Param<sidis::math::Bound> Q_sq_cut;
	Param<sidis::math::Bound> t_cut;
	Param<sidis::math::Bound> w_cut;
	Param<sidis::math::Bound> r_cut;
	Param<sidis::math::Bound> mx_sq_cut;
	Param<sidis::math::Bound> q_0_cut;
	Param<sidis::math::Bound> k2_0_cut;
	Param<sidis::math::Bound> ph_0_cut;
	Param<sidis::math::Bound> theta_q_cut;
	Param<sidis::math::Bound> theta_k2_cut;
	Param<sidis::math::Bound> theta_h_cut;
	Param<sidis::math::Bound> tau_cut;
	Param<sidis::math::Bound> phi_k_cut;
	Param<sidis::math::Bound> k_0_bar_cut;
	Param<sidis::math::Bound> k_0_cut;
	Param<sidis::math::Bound> theta_k_cut;

	Params() :
		version("version"),
		event_file("event-file"),
		rc_method("rc-method"),
		gen_nrad("gen-nrad"),
		gen_rad("gen-rad"),
		write_photon("write-photon"),
		foam_nrad_file("foam-nrad-file"),
		foam_rad_file("foam-rad-file"),
		sf_set("sf-set"),
		num_events("num-events"),
		num_init("num-init"),
		seed("seed"),
		seed_init("seed-init"),
		beam_energy("beam-energy"),
		beam("beam"),
		target("target"),
		hadron("hadron"),
		mass_threshold("mass-threshold"),
		target_pol("target-pol"),
		beam_pol("beam-pol"),
		k_0_bar("soft-threshold"),
		x_cut("x-cut"),
		y_cut("y-cut"),
		z_cut("z-cut"),
		ph_t_sq_cut("ph-t-sq-cut"),
		phi_h_cut("phi-h-cut"),
		phi_cut("phi-cut"),
		Q_sq_cut("Q-sq-cut"),
		t_cut("t-cut"),
		w_cut("w-cut"),
		r_cut("r-cut"),
		mx_sq_cut("mx-sq-cut"),
		q_0_cut("q-0-cut"),
		k2_0_cut("k2-0-cut"),
		ph_0_cut("ph-0-cut"),
		theta_q_cut("theta-q-cut"),
		theta_k2_cut("theta-k2-cut"),
		theta_h_cut("theta-h-cut"),
		tau_cut("tau-cut"),
		phi_k_cut("phi-k-cut"),
		k_0_bar_cut("k-0-bar-cut"),
		k_0_cut("k-0-cut"),
		theta_k_cut("theta-k-cut") { }

	void write_root(TFile& file) const;
	void read_root(TFile& file);

	void write_stream(std::ostream& file) const;
	void read_stream(std::istream& file);

	// Takes the supplied parameters and fills in missing ones to make a valid
	// run card. If unable to do so, throw an exception. If `strict` is enabled,
	// then `make_valid` will never change a parameter that has been set by the
	// user (for example, disabling `write-photon` when no radiative corrections
	// are being applied).
	void make_valid(bool strict=true);
	bool valid() const {
		Params other = *this;
		other.make_valid();
		return other == *this;
	}

	// Check whether a `TFoam` generated for one set of parameters can be used
	// for another set. If not compatible, throws an exception with a message
	// describing the incompatibility.
	void compatible_with_foam(Params const& foam_params) const;

	bool operator==(Params const& rhs) const;
	bool operator!=(Params const& rhs) const {
		return !(*this == rhs);
	}
};

#endif

