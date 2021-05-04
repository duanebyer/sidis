#include "params.hpp"

#include <algorithm>
#include <cctype>
#include <ios>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <type_traits>
#include <unordered_map>

#include <TArrayD.h>
#include <TArrayI.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TString.h>
#include <TVector3.h>

using namespace sidis;
using namespace sidis::math;

namespace {

template<typename T, typename=void>
struct RootParser {
	static void write_root(TFile&, Param<T> const&) {
	}
	static void read_root(TFile&, Param<T>&) {
	}
};

// Remove whitespace from beginning and end of string.
std::string trim(std::string str) {
	auto not_whitespace = [](char c) {
		return !std::isspace(c);
	};
	str.erase(str.begin(), std::find_if(str.begin(), str.end(), not_whitespace));
	str.erase(std::find_if(str.rbegin(), str.rend(), not_whitespace).base(), str.end());
	return str;
}

// Search for comment character '#' and remove anything following it.
std::string trim_comment(std::string str) {
	str.erase(std::find(str.begin(), str.end(), '#'), str.end());
	return str;
}

template<typename T>
void write_param_root(TFile& file, Param<T> const& param) {
	if (param.occupied()) {
		RootParser<T>::write_root(file, param);
	}
}
template<typename T>
void read_param_root(TFile& file, Param<T>& param) {
	RootParser<T>::read_root(file, param);
}

template<typename T>
void write_param_stream(std::ostream& os, Param<T> const& param, bool force=false) {
	if (param.occupied() || force) {
		os << param.name() << " " << *param << std::endl;
	}
	if (!os) {
		throw std::runtime_error(
			std::string("Could not write parameter '")
			+ param.name() + "' to stream.");
	}
}

template<typename T>
void read_param_stream(std::istream& is, Param<T>& param) {
	T value{};
	is >> value;
	param.reset(value);
	if (!is) {
		throw std::runtime_error(
			std::string("Could not read parameter '")
			+ param.name() + "' from stream.");
	}
}

template<>
void read_param_stream<std::string>(std::istream& is, Param<std::string>& param) {
	std::string value;
	std::getline(is, value);
	param.reset(trim(value));
	if (!is) {
		throw std::runtime_error(
			std::string("Could not read parameter '")
			+ param.name() + "' from stream.");
	}
}

template<typename T>
void consume_param_from_map(
		std::unordered_map<std::string, std::string>& map,
		Param<T>& param) {
	param.reset();
	auto p = map.find(param.name());
	if (p != map.end()) {
		// Strip all trailing whitespace or comments from param.
		std::istringstream ss(p->second);
		map.erase(p);
		try {
			read_param_stream(ss, param);
			std::string rem;
			std::getline(ss, rem);
			if (!rem.empty()) {
				throw std::runtime_error("");
			}
		} catch (std::exception const& e) {
			throw std::runtime_error(
				std::string("Failed to parse parameter '")
				+ param.name() + "' from stream.");
		}
	}
}

// Parsers for input and output of various types to ROOT files.
// Version.
template<>
struct RootParser<Version> {
	static void write_root(TFile& file, Param<Version> const& param) {
		Int_t vals[2] = { param->v_major, param->v_minor };
		TArrayI object(2, vals);
		Int_t result = file.WriteObject<TArrayI>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root(TFile& file, Param<Version>& param) {
		auto* object = file.Get<TArrayI>(param.name());
		if (object != nullptr && object->GetSize() == 2) {
			param.reset(Version(object->At(0), object->At(1)));
		} else {
			param.reset();
		}
	}
};
// Number types.
template<typename T>
struct RootParser<T, typename std::enable_if<std::is_arithmetic<T>::value>::type> {
	static void write_root(TFile& file, Param<T> const& param) {
		TParameter<T> object("", *param);
		Int_t result = file.WriteObject<TParameter<T> >(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root(TFile& file, Param<T>& param) {
		auto* object = file.Get<TParameter<T> >(param.name());
		if (object != nullptr) {
			param.reset(object->GetVal());
		} else {
			param.reset();
		}
	}
};
// Enums.
template<typename T>
using is_scoped_enum = std::integral_constant<
	bool,
	std::is_enum<T>::value && !std::is_convertible<T, int>::value>;
template<typename T>
struct RootParser<T, typename std::enable_if<is_scoped_enum<T>::value>::type> {
	using Integral = typename std::underlying_type<T>::type;
	static void write_root(TFile& file, Param<T> const& param) {
		TParameter<Integral> object("", static_cast<Integral>(*param));
		Int_t result = file.WriteObject<decltype(object)>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root(TFile& file, Param<T>& param) {
		auto object = file.Get<TParameter<Integral> >(param.name());
		if (object != nullptr) {
			param.reset(static_cast<T>(object->GetVal()));
		} else {
			param.reset();
		}
	}
};
// Strings.
template<>
struct RootParser<std::string> {
	static void write_root(TFile& file, Param<std::string> const& param) {
		TObjString object(param->c_str());
		Int_t result = file.WriteObject<TObjString>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root(TFile& file, Param<std::string>& param) {
		auto object = file.Get<TObjString>(param.name());
		if (object != nullptr) {
			param.reset(object->GetString().Data());
		} else {
			param.reset();
		}
	}
};
// Math types.
template<>
struct RootParser<Vec3> {
	static void write_root(TFile& file, Param<Vec3> const& param) {
		TVector3 object(param->x, param->y, param->z);
		Int_t result = file.WriteObject<TVector3>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root(TFile& file, Param<Vec3>& param) {
		auto object = file.Get<TVector3>(param.name());
		if (object != nullptr) {
			param.reset(Vec3(object->X(), object->Y(), object->Z()));
		} else {
			param.reset();
		}
	}
};
template<>
struct RootParser<Bound> {
	static void write_root(TFile& file, Param<Bound> const& param) {
		Double_t vals[2] = { param->min(), param->max() };
		TArrayD object(2, vals);
		Int_t result = file.WriteObject<TArrayD>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root(TFile& file, Param<Bound>& param) {
		auto object = file.Get<TArrayD>(param.name());
		if (object != nullptr && object->GetSize() == 2 && object->At(0) <= object->At(1)) {
			param.reset(Bound(object->At(0), object->At(1)));
		} else {
			param.reset();
		}
	}
};

// Parsers for streams.
// Version.
std::ostream& operator<<(std::ostream& os, Version const& version) {
	return os << version.v_major << "." << version.v_minor;
}
std::istream& operator>>(std::istream& is, Version& version) {
	is >> version.v_major;
	char period;
	is >> period;
	if (period != '.') {
		is.setstate(std::ios_base::failbit);
	}
	is >> version.v_minor;
	return is;
}
// Enums.
std::ostream& operator<<(std::ostream& os, RcMethod const& rc_method) {
	switch (rc_method) {
	case RcMethod::NONE:
		return os << "none";
	case RcMethod::APPROX:
		return os << "approx";
	case RcMethod::EXACT:
		return os << "exact";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, RcMethod& rc_method) {
	std::string name;
	is >> name;
	if (name == "none") {
		rc_method = RcMethod::NONE;
	} else if (name == "approx") {
		rc_method = RcMethod::APPROX;
	} else if (name == "exact") {
		rc_method = RcMethod::EXACT;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}
std::ostream& operator<<(std::ostream& os, part::Nucleus const& nucleus) {
	switch (nucleus) {
	case part::Nucleus::P:
		return os << "p";
	case part::Nucleus::N:
		return os << "n";
	case part::Nucleus::D:
		return os << "d";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, part::Nucleus& nucleus) {
	std::string name;
	is >> name;
	if (name == "p" || name == "proton") {
		nucleus = part::Nucleus::P;
	} else if (name == "n" || name == "neutron") {
		nucleus = part::Nucleus::N;
	} else if (name == "d" || name == "deuteron") {
		nucleus = part::Nucleus::D;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}
std::ostream& operator<<(std::ostream& os, part::Lepton const& lepton) {
	switch (lepton) {
	case part::Lepton::E:
		return os << "e";
	case part::Lepton::MU:
		return os << "mu";
	case part::Lepton::TAU:
		return os << "tau";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, part::Lepton& lepton) {
	std::string name;
	is >> name;
	if (name == "e" || name == "electron") {
		lepton = part::Lepton::E;
	} else if (name == "mu" || name == "muon") {
		lepton = part::Lepton::MU;
	} else if (name == "tau" || name == "tauon") {
		lepton = part::Lepton::TAU;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}
std::ostream& operator<<(std::ostream& os, part::Hadron const& hadron) {
	switch (hadron) {
	case part::Hadron::PI_0:
		return os << "pi0";
	case part::Hadron::PI_P:
		return os << "pi+";
	case part::Hadron::PI_M:
		return os << "pi-";
	case part::Hadron::K_0:
		return os << "K0";
	case part::Hadron::K_P:
		return os << "K+";
	case part::Hadron::K_M:
		return os << "K-";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, part::Hadron& hadron) {
	std::string name;
	is >> name;
	if (name == "pi0" || name == "pion0") {
		hadron = part::Hadron::PI_0;
	} else if (name == "pi+" || name == "pion+") {
		hadron = part::Hadron::PI_P;
	} else if (name == "pi-" || name == "pion-") {
		hadron = part::Hadron::PI_M;
	} else if (name == "K0" || name == "kaon0") {
		hadron = part::Hadron::K_0;
	} else if (name == "K+" || name == "kaon+") {
		hadron = part::Hadron::K_P;
	} else if (name == "K-" || name == "kaon-") {
		hadron = part::Hadron::K_M;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}
// Math types.
std::ostream& operator<<(std::ostream& os, Vec3 const& vec) {
	return os << vec.x << " " << vec.y << " " << vec.z;
}
std::istream& operator>>(std::istream& is, Vec3& vec) {
	is >> vec.x >> vec.y >> vec.z;
	return is;
}
std::ostream& operator<<(std::ostream& os, Bound const& bounds) {
	os << bounds.min() << " " << bounds.max();
	return os;
}
std::istream& operator>>(std::istream& is, Bound& bounds) {
	Real min;
	Real max;
	is >> min >> max;
	if (min <= max) {
		bounds = Bound(min, max);
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}

}

void Params::write_root(TFile& file) const {
	file.cd();
	write_param_root(file, version);
	write_param_root(file, event_file);
	write_param_root(file, rc_method);
	write_param_root(file, gen_nrad);
	write_param_root(file, gen_rad);
	write_param_root(file, write_photon);
	write_param_root(file, foam_nrad_file);
	write_param_root(file, foam_rad_file);
	write_param_root(file, sf_set);
	write_param_root(file, num_events);
	write_param_root(file, num_init);
	write_param_root(file, seed);
	write_param_root(file, seed_init);
	write_param_root(file, beam_energy);
	write_param_root(file, beam);
	write_param_root(file, target);
	write_param_root(file, hadron);
	write_param_root(file, mass_threshold);
	write_param_root(file, target_pol);
	write_param_root(file, beam_pol);
	write_param_root(file, k_0_bar);
	write_param_root(file, x_cut);
	write_param_root(file, y_cut);
	write_param_root(file, z_cut);
	write_param_root(file, ph_t_sq_cut);
	write_param_root(file, phi_h_cut);
	write_param_root(file, phi_cut);
	write_param_root(file, Q_sq_cut);
	write_param_root(file, t_cut);
	write_param_root(file, W_sq_cut);
	write_param_root(file, r_cut);
	write_param_root(file, mx_sq_cut);
	write_param_root(file, q_0_cut);
	write_param_root(file, k2_0_cut);
	write_param_root(file, ph_0_cut);
	write_param_root(file, theta_q_cut);
	write_param_root(file, theta_k2_cut);
	write_param_root(file, theta_h_cut);
	write_param_root(file, tau_cut);
	write_param_root(file, phi_k_cut);
	write_param_root(file, R_cut);
	write_param_root(file, k_0_bar_cut);
	write_param_root(file, k_0_cut);
	write_param_root(file, theta_k_cut);
}

void Params::read_root(TFile& file) {
	read_param_root(file, version);
	if (version->v_major != SIDIS_PARAMS_VERSION_MAJOR
			|| version->v_minor > SIDIS_PARAMS_VERSION_MINOR) {
		throw std::runtime_error(
			std::string("Cannot read parameters from version ")
			+ std::to_string(version->v_major) + "."
			+ std::to_string(version->v_minor) + ".");
	}
	read_param_root(file, event_file);
	read_param_root(file, rc_method);
	read_param_root(file, gen_nrad);
	read_param_root(file, gen_rad);
	read_param_root(file, write_photon);
	read_param_root(file, foam_nrad_file);
	read_param_root(file, foam_rad_file);
	read_param_root(file, sf_set);
	read_param_root(file, num_events);
	read_param_root(file, num_init);
	read_param_root(file, seed);
	read_param_root(file, seed_init);
	read_param_root(file, beam_energy);
	read_param_root(file, beam);
	read_param_root(file, target);
	read_param_root(file, hadron);
	read_param_root(file, mass_threshold);
	read_param_root(file, target_pol);
	read_param_root(file, beam_pol);
	read_param_root(file, k_0_bar);
	read_param_root(file, x_cut);
	read_param_root(file, y_cut);
	read_param_root(file, z_cut);
	read_param_root(file, ph_t_sq_cut);
	read_param_root(file, phi_h_cut);
	read_param_root(file, phi_cut);
	read_param_root(file, Q_sq_cut);
	read_param_root(file, t_cut);
	read_param_root(file, W_sq_cut);
	read_param_root(file, r_cut);
	read_param_root(file, mx_sq_cut);
	read_param_root(file, q_0_cut);
	read_param_root(file, k2_0_cut);
	read_param_root(file, ph_0_cut);
	read_param_root(file, theta_q_cut);
	read_param_root(file, theta_k2_cut);
	read_param_root(file, theta_h_cut);
	read_param_root(file, tau_cut);
	read_param_root(file, phi_k_cut);
	read_param_root(file, R_cut);
	read_param_root(file, k_0_bar_cut);
	read_param_root(file, k_0_cut);
	read_param_root(file, theta_k_cut);
}

void Params::write_stream(std::ostream& file) const {
	write_param_stream(file, version, true);
	write_param_stream(file, event_file);
	write_param_stream(file, rc_method);
	write_param_stream(file, gen_nrad);
	write_param_stream(file, gen_rad);
	write_param_stream(file, write_photon);
	write_param_stream(file, foam_nrad_file);
	write_param_stream(file, foam_rad_file);
	write_param_stream(file, sf_set);
	write_param_stream(file, num_events);
	write_param_stream(file, num_init);
	write_param_stream(file, seed);
	write_param_stream(file, seed_init);
	write_param_stream(file, beam_energy);
	write_param_stream(file, beam);
	write_param_stream(file, target);
	write_param_stream(file, hadron);
	write_param_stream(file, mass_threshold);
	write_param_stream(file, target_pol);
	write_param_stream(file, beam_pol);
	write_param_stream(file, k_0_bar);
	write_param_stream(file, x_cut);
	write_param_stream(file, y_cut);
	write_param_stream(file, z_cut);
	write_param_stream(file, ph_t_sq_cut);
	write_param_stream(file, phi_h_cut);
	write_param_stream(file, phi_cut);
	write_param_stream(file, Q_sq_cut);
	write_param_stream(file, t_cut);
	write_param_stream(file, W_sq_cut);
	write_param_stream(file, r_cut);
	write_param_stream(file, mx_sq_cut);
	write_param_stream(file, q_0_cut);
	write_param_stream(file, k2_0_cut);
	write_param_stream(file, ph_0_cut);
	write_param_stream(file, theta_q_cut);
	write_param_stream(file, theta_k2_cut);
	write_param_stream(file, theta_h_cut);
	write_param_stream(file, tau_cut);
	write_param_stream(file, phi_k_cut);
	write_param_stream(file, R_cut);
	write_param_stream(file, k_0_bar_cut);
	write_param_stream(file, k_0_cut);
	write_param_stream(file, theta_k_cut);
}
void Params::read_stream(std::istream& file) {
	std::unordered_map<std::string, std::string> map;
	while (file) {
		std::string line;
		std::getline(file, line);
		std::stringstream ss(trim(trim_comment(line)));
		std::string key;
		std::string value;
		ss >> key;
		std::getline(ss, value);
		value = trim(value);
		if (!key.empty()) {
			if (map.find(key) != map.end()) {
				throw std::runtime_error("Duplicate parameter '" + key + "'.");
			}
			map[key] = value;
		}
	}

	consume_param_from_map(map, version);
	version.get_or_insert(Version());
	if (version->v_major != SIDIS_PARAMS_VERSION_MAJOR
			|| version->v_minor > SIDIS_PARAMS_VERSION_MINOR) {
		throw std::runtime_error(
			std::string("Cannot read parameters from version ")
			+ std::to_string(version->v_major) + "."
			+ std::to_string(version->v_minor) + ".");
	}
	consume_param_from_map(map, event_file);
	consume_param_from_map(map, rc_method);
	consume_param_from_map(map, gen_nrad);
	consume_param_from_map(map, gen_rad);
	consume_param_from_map(map, write_photon);
	consume_param_from_map(map, foam_nrad_file);
	consume_param_from_map(map, foam_rad_file);
	consume_param_from_map(map, sf_set);
	consume_param_from_map(map, num_events);
	consume_param_from_map(map, num_init);
	consume_param_from_map(map, seed);
	consume_param_from_map(map, seed_init);
	consume_param_from_map(map, beam_energy);
	consume_param_from_map(map, beam);
	consume_param_from_map(map, target);
	consume_param_from_map(map, hadron);
	consume_param_from_map(map, mass_threshold);
	consume_param_from_map(map, target_pol);
	consume_param_from_map(map, beam_pol);
	consume_param_from_map(map, k_0_bar);
	consume_param_from_map(map, x_cut);
	consume_param_from_map(map, y_cut);
	consume_param_from_map(map, z_cut);
	consume_param_from_map(map, ph_t_sq_cut);
	consume_param_from_map(map, phi_h_cut);
	consume_param_from_map(map, phi_cut);
	consume_param_from_map(map, Q_sq_cut);
	consume_param_from_map(map, t_cut);
	consume_param_from_map(map, W_sq_cut);
	consume_param_from_map(map, r_cut);
	consume_param_from_map(map, mx_sq_cut);
	consume_param_from_map(map, q_0_cut);
	consume_param_from_map(map, k2_0_cut);
	consume_param_from_map(map, ph_0_cut);
	consume_param_from_map(map, theta_q_cut);
	consume_param_from_map(map, theta_k2_cut);
	consume_param_from_map(map, theta_h_cut);
	consume_param_from_map(map, tau_cut);
	consume_param_from_map(map, phi_k_cut);
	consume_param_from_map(map, R_cut);
	consume_param_from_map(map, k_0_bar_cut);
	consume_param_from_map(map, k_0_cut);
	consume_param_from_map(map, theta_k_cut);

	if (!map.empty()) {
		std::ostringstream ss_err;
		ss_err << "Unrecognized parameters";
		while (!map.empty()) {
			auto p = map.begin();
			ss_err << " '" << p->first << "'";
			map.erase(p);
		}
		throw std::runtime_error(ss_err.str());
	}
}

void Params::make_valid(bool strict) {
	// This function is somewhat complicated because it tries to consider many
	// possible combinations of options and sort through them to produce a
	// "normalized" form of any parameter file.
	version.get_or_insert(Version());
	if (!num_events.occupied()) {
		throw std::runtime_error("Must specify number of events to generate.");
	}
	if (!beam_energy.occupied()) {
		throw std::runtime_error("Must specify beam energy.");
	}
	if (!beam.occupied()) {
		throw std::runtime_error("Must specify beam lepton flavor.");
	}
	if (!target.occupied()) {
		throw std::runtime_error("Must specify target nucleus.");
	}
	if (!hadron.occupied()) {
		throw std::runtime_error("Must specify leading hadron.");
	}
	if (!mass_threshold.occupied()) {
		throw std::runtime_error("Must specify mass threshold.");
	}
	if (!sf_set.occupied()) {
		throw std::runtime_error("Must specify structure function set.");
	}
	event_file.get_or_insert("gen.root");
	rc_method.get_or_insert(RcMethod::APPROX);
	num_init.get_or_insert(1000);
	seed.get_or_insert(0);
	seed_init.get_or_insert(0);
	target_pol.get_or_insert(Vec3::ZERO);
	beam_pol.get_or_insert(0.);
	// RC methods.
	// The process for determining whether radiative and non-radiative events
	// are generated looks something like this:
	//  * If RC are disabled, then no radiative events are generated.
	//  * If `k_0_bar_cut` is above the soft threshold, then only radiative
	//    events are generated.
	//  * `k_0_bar_cut` cannot be partway between zero and the soft threshold,
	//    because it is impossible to generate events according to that.
	if (*rc_method == RcMethod::APPROX || *rc_method == RcMethod::EXACT) {
		gen_rad.get_or_insert(true);
		k_0_bar.get_or_insert(0.01);
		if (gen_nrad.occupied() && !k_0_bar_cut.occupied()) {
			if (*gen_nrad) {
				k_0_bar_cut.reset(Bound::POSITIVE);
			} else {
				k_0_bar_cut.reset(Bound::POSITIVE + *k_0_bar);
			}
		} else {
			k_0_bar_cut.get_or_insert(Bound::POSITIVE);
		}
	} else {
		if (gen_rad.occupied() && *gen_rad) {
			if (strict) {
				throw std::runtime_error(
					"Cannot generate radiative events without "
					"radiative corrections enabled.");
			}
		}
		gen_rad.reset(false);
		if (k_0_bar.occupied()) {
			if (strict) {
				throw std::runtime_error(
					"Cannot provide soft photon threshold when not using any "
					"radiative cross-sections.");
			}
			k_0_bar.reset();
		}
	}
	if (*gen_rad) {
		if (k_0_bar_cut->min() <= 0.) {
			gen_nrad.get_or_insert(true);
			if (!*gen_nrad) {
				if (strict) {
					throw std::runtime_error(
						"Cannot generate radiative events only with "
						"`k_0_bar_cut min. == 0`.");
				}
			}
		} else if (*k_0_bar <= k_0_bar_cut->min()) {
			if (gen_nrad.occupied() && *gen_nrad) {
				if (strict) {
					throw std::runtime_error(
						"Cannot generate non-radiative events with "
						"`soft_threshold < k_0_bar_cut min.`.");
				}
			}
			gen_nrad.reset(false);
		} else {
			throw std::runtime_error(
				"Invalid `k_0_bar_cut`: minimum must be zero "
				"or larger than `soft_threshold`.");
		}
	} else {
		gen_nrad.get_or_insert(true);
	}
	if (!*gen_nrad && !*gen_rad) {
		throw std::runtime_error(
			"Cannot disable all event types `(nrad, rad)`.");
	}
	// Basic options associated with radiative and non-radiative events in
	// particular.
	if (*gen_nrad) {
		foam_nrad_file.get_or_insert("foam-nrad.root");
	} else if (foam_nrad_file.occupied()) {
		if (strict) {
			throw std::runtime_error(
				"Cannot provide `foam_nrad_file` when no "
				"non-radiative events are being generated.");
		} else {
			foam_nrad_file.reset();
		}
	}
	if (*gen_rad) {
		foam_rad_file.get_or_insert("foam-rad.root");
		write_photon.get_or_insert(true);
	} else {
		if (foam_rad_file.occupied()) {
			if (strict) {
				throw std::runtime_error(
					"Cannot provide `foam_rad_file` when no "
					"radiative events are being generated.");
			} else {
				foam_rad_file.reset();
			}
		}
		if (write_photon.occupied()) {
			if (strict) {
				throw std::runtime_error(
					"Cannot enable `write_photon` when no "
					"radiative events are being generated.");
			} else {
				write_photon.reset();
			}
		}
	}
	// Cuts.
	if (*gen_rad && *gen_nrad) {
		if (tau_cut.occupied()
				|| phi_k_cut.occupied()
				|| R_cut.occupied()
				|| k_0_cut.occupied()
				|| theta_k_cut.occupied()) {
			throw std::runtime_error(
				"Cannot apply radiative cuts to non-radiative events.");
		}
	}
	// Verify that cuts make sense. This isn't comprehensive, but is primarily
	// important to avoid cuts on the azimuthal angles larger than 360 degrees.
	if (strict && x_cut.occupied() && !Bound::UNIT.contains(*x_cut)) {
		throw std::runtime_error(
			"Cut on x must lie between 0 and 1.");
	}
	if (strict && y_cut.occupied() && !Bound::UNIT.contains(*y_cut)) {
		throw std::runtime_error(
			"Cut on y must lie between 0 and 1.");
	}
	if (strict && z_cut.occupied() && !Bound::UNIT.contains(*z_cut)) {
		throw std::runtime_error(
			"Cut on z must lie between 0 and 1.");
	}
	if (phi_h_cut.occupied() && phi_h_cut->size() >= 360.) {
		throw std::runtime_error(
			"Cut on φ_h must be smaller than 360 degrees.");
	}
	if (phi_cut.occupied() && phi_cut->size() >= 360.) {
		throw std::runtime_error(
			"Cut on φ must be smaller than 360 degrees.");
	}
	if (phi_k_cut.occupied() && phi_k_cut->size() >= 360.) {
		throw std::runtime_error(
			"Cut on φ_k must be smaller than 360 degrees.");
	}
	if (strict && theta_q_cut.occupied() && !Bound(0., 180.).contains(*theta_q_cut)) {
		throw std::runtime_error(
			"Cut on θ_q must lie between 0 and 180 degrees.");
	}
	if (strict && theta_h_cut.occupied() && !Bound(0., 180.).contains(*theta_h_cut)) {
		throw std::runtime_error(
			"Cut on θ_h must lie between 0 and 180 degrees.");
	}
	if (strict && theta_k_cut.occupied() && !Bound(0., 180.).contains(*theta_k_cut)) {
		throw std::runtime_error(
			"Cut on θ_k must lie between 0 and 180 degrees.");
	}
}

void Params::compatible_with_foam(Params const& foam_params) const {
	if (version->v_major != foam_params.version->v_major
			|| version->v_minor > foam_params.version->v_minor) {
		throw std::runtime_error("Incompatible versions.");
	} else if (*foam_params.rc_method != *rc_method) {
		throw std::runtime_error("Different RC methods.");
	} else if (*gen_nrad && !*foam_params.gen_nrad) {
		throw std::runtime_error("No non-radiative FOAM available.");
	} else if (*gen_rad && !*foam_params.gen_rad) {
		throw std::runtime_error("No radiative FOAM available.");
	} else if (*sf_set != *foam_params.sf_set) {
		throw std::runtime_error("Different SF set.");
	} else if (*num_init < *foam_params.num_init) {
		throw std::runtime_error("Insufficient initialization sampling.");
	} else if (*seed_init != 0 && *foam_params.seed_init != *seed_init) {
		throw std::runtime_error("Different initialization seed.");
	} else if (*beam_energy != *foam_params.beam_energy) {
		throw std::runtime_error("Different beam energies.");
	} else if (*beam != *foam_params.beam) {
		throw std::runtime_error("Different beam lepton type.");
	} else if (*target != *foam_params.target) {
		throw std::runtime_error("Different target nucleus type.");
	} else if (*hadron != *foam_params.hadron) {
		throw std::runtime_error("Different leading hadron type.");
	} else if (*mass_threshold != *foam_params.mass_threshold) {
		throw std::runtime_error("Different mass thresholds.");
	} else if (*target_pol != *foam_params.target_pol) {
		throw std::runtime_error("Different target polarizations.");
	} else if (*beam_pol != *foam_params.beam_pol) {
		throw std::runtime_error("Different beam polarizations.");
	} else if ((*rc_method == RcMethod::APPROX || *rc_method == RcMethod::EXACT)
			&& *k_0_bar != *foam_params.k_0_bar) {
		throw std::runtime_error("Different soft photon cutoffs.");
	} else if (x_cut.get_or(Bound::UNIT) != foam_params.x_cut.get_or(Bound::UNIT)) {
		throw std::runtime_error("Different cuts on x.");
	} else if (y_cut.get_or(Bound::UNIT) != foam_params.y_cut.get_or(Bound::UNIT)) {
		throw std::runtime_error("Different cuts on y.");
	} else if (z_cut.get_or(Bound::UNIT) != foam_params.z_cut.get_or(Bound::UNIT)) {
		throw std::runtime_error("Different cuts on z.");
	} else if (ph_t_sq_cut != foam_params.ph_t_sq_cut) {
		throw std::runtime_error("Different cuts on ph_t².");
	} else if (phi_h_cut != foam_params.phi_h_cut) {
		throw std::runtime_error("Different cuts on φ_h.");
	} else if (phi_cut != foam_params.phi_cut) {
		throw std::runtime_error("Different cuts on φ.");
	} else if (Q_sq_cut != foam_params.Q_sq_cut) {
		throw std::runtime_error("Different cuts on Q².");
	} else if (t_cut != foam_params.t_cut) {
		throw std::runtime_error("Different cuts on t.");
	} else if (W_sq_cut != foam_params.W_sq_cut) {
		throw std::runtime_error("Different cuts on W².");
	} else if (r_cut != foam_params.r_cut) {
		throw std::runtime_error("Different cuts on r.");
	} else if (mx_sq_cut != foam_params.mx_sq_cut) {
		throw std::runtime_error("Different cuts on mx².");
	} else if (q_0_cut != foam_params.q_0_cut) {
		throw std::runtime_error("Different cuts on q_0.");
	} else if (k2_0_cut != foam_params.k2_0_cut) {
		throw std::runtime_error("Different cuts on k2_0.");
	} else if (ph_0_cut != foam_params.ph_0_cut) {
		throw std::runtime_error("Different cuts on ph_0.");
	} else if (theta_q_cut != foam_params.theta_q_cut) {
		throw std::runtime_error("Different cuts on θ_q.");
	} else if (theta_k2_cut != foam_params.theta_k2_cut) {
		throw std::runtime_error("Different cuts on θ_k2.");
	} else if (theta_h_cut != foam_params.theta_h_cut) {
		throw std::runtime_error("Different cuts on θ_h.");
	} else if (*gen_rad && tau_cut != foam_params.tau_cut) {
		throw std::runtime_error("Different cuts on τ.");
	} else if (*gen_rad && phi_k_cut != foam_params.phi_k_cut) {
		throw std::runtime_error("Different cuts on φ_k.");
	} else if (*gen_rad && R_cut != foam_params.R_cut) {
		throw std::runtime_error("Different cuts on R.");
	} else if (*gen_rad && k_0_bar_cut != foam_params.k_0_bar_cut) {
		throw std::runtime_error("Different cuts on radiated photon energy.");
	} else if (*gen_rad && k_0_cut != foam_params.k_0_cut) {
		throw std::runtime_error("Different cuts on radiated photon energy.");
	} else if (*gen_rad && theta_k_cut != foam_params.theta_k_cut) {
		throw std::runtime_error("Different cuts on radiated photon energy.");
	}
}

bool Params::operator==(Params const& rhs) const {
	return version == rhs.version
		&& event_file == rhs.event_file
		&& rc_method == rhs.rc_method
		&& gen_nrad == rhs.gen_nrad
		&& gen_rad == rhs.gen_rad
		&& write_photon == rhs.write_photon
		&& foam_nrad_file == rhs.foam_nrad_file
		&& foam_rad_file == rhs.foam_rad_file
		&& sf_set == rhs.sf_set
		&& num_events == rhs.num_events
		&& num_init == rhs.num_init
		&& seed == rhs.seed
		&& seed_init == rhs.seed_init
		&& beam_energy == rhs.beam_energy
		&& beam == rhs.beam
		&& target == rhs.target
		&& hadron == rhs.hadron
		&& mass_threshold == rhs.mass_threshold
		&& target_pol == rhs.target_pol
		&& beam_pol == rhs.beam_pol
		&& k_0_bar == rhs.k_0_bar
		&& x_cut == rhs.x_cut
		&& y_cut == rhs.y_cut
		&& z_cut == rhs.z_cut
		&& ph_t_sq_cut == rhs.ph_t_sq_cut
		&& phi_h_cut == rhs.phi_h_cut
		&& phi_cut == rhs.phi_cut
		&& Q_sq_cut == rhs.Q_sq_cut
		&& t_cut == rhs.t_cut
		&& W_sq_cut == rhs.W_sq_cut
		&& r_cut == rhs.r_cut
		&& mx_sq_cut == rhs.mx_sq_cut
		&& q_0_cut == rhs.q_0_cut
		&& k2_0_cut == rhs.k2_0_cut
		&& ph_0_cut == rhs.ph_0_cut
		&& theta_q_cut == rhs.theta_q_cut
		&& theta_k2_cut == rhs.theta_k2_cut
		&& theta_h_cut == rhs.theta_h_cut
		&& tau_cut == rhs.tau_cut
		&& phi_k_cut == rhs.phi_k_cut
		&& R_cut == rhs.R_cut
		&& k_0_bar_cut == rhs.k_0_bar_cut
		&& k_0_cut == rhs.k_0_cut
		&& theta_k_cut == rhs.theta_k_cut;
}

