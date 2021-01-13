#include "params.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <TObjString.h>
#include <TParameter.h>
#include <TString.h>
#include <TVector3.h>

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::math;

#define WRITE_PARAM_ROOT(file, param) \
	{ \
		auto p = convert(param); \
		file.WriteObject(&p, #param); \
	}

#define READ_PARAM_ROOT(file, param) \
	{ \
		auto p = file.Get<decltype(convert(param))>(#param); \
		if (p != nullptr) { \
			param = static_cast<decltype(param)>(convert(*p)); \
		} \
	}

#define WRITE_PARAM(file, param) \
	file << #param << " " << param << std::endl;

#define READ_PARAM(params, param) \
	{ \
		auto p = params.find(#param); \
		if (p != params.end()) { \
			std::istringstream ss(p->second); \
			params.erase(p); \
			ss >> param; \
			if (!ss) { \
				throw std::runtime_error("Failed to parse parameter '" #param "'"); \
			} \
		} \
	}

template<typename T>
using is_scoped_enum = std::integral_constant<
	bool,
	std::is_enum<T>::value && !std::is_convertible<T, int>::value>;

template<
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
T convert(TParameter<T> const& in) {
	return in.GetVal();
}
template<
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
TParameter<T> convert(T const& in) {
	return TParameter<T>("", in);
}
template<
	typename T,
	typename = typename std::enable_if<is_scoped_enum<T>::value>::type>
T convert(TParameter<typename std::underlying_type<T>::type> const& in) {
	return static_cast<T>(in.GetVal());
}
template<
	typename T,
	typename = typename std::enable_if<is_scoped_enum<T>::value>::type>
TParameter<typename std::underlying_type<T>::type> convert(T const& in) {
	return TParameter<typename std::underlying_type<T>::type>(
		"", static_cast<typename std::underlying_type<T>::type>(in));
}
TObjString convert(std::string str) {
	return TObjString(str.c_str());
}
std::string convert(TObjString str) {
	return str.GetString().Data();
}
TVector3 convert(Vec3 const& vec) {
	return TVector3(vec.x, vec.y, vec.z);
}
Vec3 convert(TVector3 const& vec) {
	return Vec3(vec.X(), vec.Y(), vec.Z());
}

std::ostream& operator<<(std::ostream& os, Nucleus const& nucleus) {
	switch (nucleus) {
	case Nucleus::P:
		return os << "p";
	case Nucleus::N:
		return os << "n";
	case Nucleus::D:
		return os << "d";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, Nucleus& nucleus) {
	std::string name;
	is >> name;
	if (name == "p" || name == "proton") {
		nucleus = Nucleus::P;
	} else if (name == "n" || name == "neutron") {
		nucleus = Nucleus::N;
	} else if (name == "d" || name == "deuteron") {
		nucleus = Nucleus::D;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}

std::ostream& operator<<(std::ostream& os, Lepton const& lepton) {
	switch (lepton) {
	case Lepton::E:
		return os << "e";
	case Lepton::MU:
		return os << "mu";
	case Lepton::TAU:
		return os << "tau";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, Lepton& lepton) {
	std::string name;
	is >> name;
	if (name == "e" || name == "electron") {
		lepton = Lepton::E;
	} else if (name == "mu" || name == "muon") {
		lepton = Lepton::MU;
	} else if (name == "tau" || name == "tauon") {
		lepton = Lepton::TAU;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}

std::ostream& operator<<(std::ostream& os, Hadron const& hadron) {
	switch (hadron) {
	case Hadron::PI_0:
		return os << "pi0";
	case Hadron::PI_P:
		return os << "pi+";
	case Hadron::PI_M:
		return os << "pi-";
	case Hadron::K_0:
		return os << "K0";
	case Hadron::K_P:
		return os << "K+";
	case Hadron::K_M:
		return os << "K-";
	default:
		os.setstate(std::ios_base::failbit);
		return os;
	}
}
std::istream& operator>>(std::istream& is, Hadron& hadron) {
	std::string name;
	is >> name;
	if (name == "pi0" || name == "pion0") {
		hadron = Hadron::PI_0;
	} else if (name == "pi+" || name == "pion+") {
		hadron = Hadron::PI_P;
	} else if (name == "pi-" || name == "pion-") {
		hadron = Hadron::PI_M;
	} else if (name == "K0" || name == "kaon0") {
		hadron = Hadron::K_0;
	} else if (name == "K+" || name == "kaon+") {
		hadron = Hadron::K_P;
	} else if (name == "K-" || name == "kaon-") {
		hadron = Hadron::K_M;
	} else {
		is.setstate(std::ios_base::failbit);
	}
	return is;
}

std::ostream& operator<<(std::ostream& os, Vec3 const& vec) {
	os << vec.x << " " << vec.y << " " << vec.z;
	return os;
}
std::istream& operator>>(std::istream& is, Vec3& vec) {
	is >> vec.x >> vec.y >> vec.z;
	return is;
}

void Params::write_root(TFile& file) const {
	file.cd();
	WRITE_PARAM_ROOT(file, event_file);
	WRITE_PARAM_ROOT(file, foam_nrad_file);
	WRITE_PARAM_ROOT(file, foam_rad_file);
	WRITE_PARAM_ROOT(file, num_events);
	WRITE_PARAM_ROOT(file, num_init);
	WRITE_PARAM_ROOT(file, seed);
	WRITE_PARAM_ROOT(file, seed_init);
	WRITE_PARAM_ROOT(file, beam_energy);
	WRITE_PARAM_ROOT(file, beam);
	WRITE_PARAM_ROOT(file, target);
	WRITE_PARAM_ROOT(file, hadron);
	WRITE_PARAM_ROOT(file, target_pol);
	WRITE_PARAM_ROOT(file, beam_pol);
	WRITE_PARAM_ROOT(file, k0_cut);
}

void Params::read_root(TFile& file) {
	READ_PARAM_ROOT(file, event_file);
	READ_PARAM_ROOT(file, foam_nrad_file);
	READ_PARAM_ROOT(file, foam_rad_file);
	READ_PARAM_ROOT(file, num_events);
	READ_PARAM_ROOT(file, num_init);
	READ_PARAM_ROOT(file, seed);
	READ_PARAM_ROOT(file, seed_init);
	READ_PARAM_ROOT(file, beam_energy);
	READ_PARAM_ROOT(file, beam);
	READ_PARAM_ROOT(file, target);
	READ_PARAM_ROOT(file, hadron);
	READ_PARAM_ROOT(file, target_pol);
	READ_PARAM_ROOT(file, beam_pol);
	READ_PARAM_ROOT(file, k0_cut);
}

void Params::write(std::ostream& file) const {
	// TODO: Set floating point precision so that all digits are kept, then
	// unset after.
	WRITE_PARAM(file, event_file);
	WRITE_PARAM(file, foam_nrad_file);
	WRITE_PARAM(file, foam_rad_file);
	WRITE_PARAM(file, num_events);
	WRITE_PARAM(file, num_init);
	WRITE_PARAM(file, seed);
	WRITE_PARAM(file, seed_init);
	WRITE_PARAM(file, beam_energy);
	WRITE_PARAM(file, beam);
	WRITE_PARAM(file, target);
	WRITE_PARAM(file, hadron);
	WRITE_PARAM(file, target_pol);
	WRITE_PARAM(file, beam_pol);
	WRITE_PARAM(file, k0_cut);
}
void Params::read(std::istream& file) {
	std::unordered_map<std::string, std::string> params;
	while (file) {
		std::string key;
		std::string value;
		file >> key;
		auto whitespace = [](char c) {
			return std::isspace(c);
		};
		if (std::all_of(key.begin(), key.end(), whitespace)) {
			continue;
		}
		std::getline(file, value);
		params[key] = value;
	}

	READ_PARAM(params, event_file);
	READ_PARAM(params, foam_nrad_file);
	READ_PARAM(params, foam_rad_file);
	READ_PARAM(params, num_events);
	READ_PARAM(params, num_init);
	READ_PARAM(params, seed);
	READ_PARAM(params, seed_init);
	READ_PARAM(params, beam_energy);
	READ_PARAM(params, beam);
	READ_PARAM(params, target);
	READ_PARAM(params, hadron);
	READ_PARAM(params, target_pol);
	READ_PARAM(params, beam_pol);
	READ_PARAM(params, k0_cut);

	if (!params.empty()) {
		std::ostringstream ss_err;
		ss_err << "Unrecognized parameters";
		while (!params.empty()) {
			auto p = params.begin();
			ss_err << " '" << p->first << "'";
			params.erase(p);
		}
		throw std::runtime_error(ss_err.str());
	}
}

bool Params::compatible_foam(Params const& foam_params) const {
	return foam_params.num_init >= num_init
		&& (foam_params.seed_init == seed_init || 0 == seed_init)
		&& foam_params.beam_energy == beam_energy
		&& foam_params.beam == beam
		&& foam_params.target == target
		&& foam_params.hadron == hadron
		&& foam_params.target_pol.x == target_pol.x
		&& foam_params.target_pol.y == target_pol.y
		&& foam_params.target_pol.z == target_pol.z
		&& foam_params.beam_pol == beam_pol
		&& foam_params.k0_cut == k0_cut;
}

