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
	Param<Toggle> gen_nrad;
	Param<Toggle> gen_rad;
	Param<Toggle> write_photon;
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
	Param<sidis::Real> Mth;
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
	Param<sidis::math::Bound> W_sq_cut;
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
	Param<sidis::math::Bound> R_cut;
	Param<sidis::math::Bound> k_0_bar_cut;
	Param<sidis::math::Bound> k_0_cut;
	Param<sidis::math::Bound> theta_k_cut;

	Params() :
		version("version"),
		strict("strict"),
		event_file("event_file"),
		rc_method("rc_method"),
		gen_nrad("gen_nrad"),
		gen_rad("gen_rad"),
		write_photon("write_photon"),
		foam_nrad_file("foam_nrad_file"),
		foam_rad_file("foam_rad_file"),
		sf_set("sf_set"),
		num_events("num_events"),
		num_init("num_init"),
		seed("seed"),
		seed_init("seed_init"),
		beam_energy("beam_energy"),
		beam("beam"),
		target("target"),
		hadron("hadron"),
		Mth("mass_threshold"),
		target_pol("target_pol"),
		beam_pol("beam_pol"),
		k_0_bar("soft_threshold"),
		x_cut("x_cut"),
		y_cut("y_cut"),
		z_cut("z_cut"),
		ph_t_sq_cut("ph_t_sq_cut"),
		phi_h_cut("phi_h_cut"),
		phi_cut("phi_cut"),
		Q_sq_cut("Q_sq_cut"),
		t_cut("t_cut"),
		W_sq_cut("W_sq_cut"),
		r_cut("r_cut"),
		mx_sq_cut("mx_sq_cut"),
		q_0_cut("q_0_cut"),
		k2_0_cut("k2_0_cut"),
		ph_0_cut("ph_0_cut"),
		theta_q_cut("theta_q_cut"),
		theta_k2_cut("theta_k2_cut"),
		theta_h_cut("theta_h_cut"),
		tau_cut("tau_cut"),
		phi_k_cut("phi_k_cut"),
		R_cut("R_cut"),
		k_0_bar_cut("k_0_bar_cut"),
		k_0_cut("k_0_cut"),
		theta_k_cut("theta_k_cut") { }

	void write_root(TFile& file) const;
	void read_root(TFile& file);

	void write_stream(std::ostream& file) const;
	void read_stream(std::istream& file);

	// Takes the supplied parameters and fills in missing ones to make a valid
	// run card. If unable to do so, throw an exception. If `strict` is enabled,
	// then `make_valid` will never change a parameter that has been set by the
	// user (for example, disabling `write_photon` when no radiative corrections
	// are being applied).
	void make_valid();
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

