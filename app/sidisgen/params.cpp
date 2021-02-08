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
	do { \
		auto p = root_from_param(param); \
		file.WriteObject(&p, #param); \
	} while (false)

#define READ_PARAM_ROOT(file, param) \
	do { \
		auto p = file.Get<decltype(root_from_param(param))>(#param); \
		if (p != nullptr) { \
			param = param_from_root<decltype(param)>(*p); \
		} \
	} while (false)

#define WRITE_PARAM(file, param) \
	do { \
		file << #param << " " << param << std::endl; \
	} while (false)

#define READ_PARAM(params, param) \
	do { \
		auto p = params.find(#param); \
		if (p != params.end()) { \
			std::istringstream ss(p->second); \
			params.erase(p); \
			ss >> param; \
			if (!ss) { \
				throw std::runtime_error("Failed to parse parameter '" #param "'"); \
			} \
		} \
	} while (false)

template<typename T>
using is_scoped_enum = std::integral_constant<
	bool,
	std::is_enum<T>::value && !std::is_convertible<T, int>::value>;

template<typename T, typename R>
T param_from_root(R const& in) {
	throw std::runtime_error("Unsupported parameter conversion");
}
template<>
Int_t param_from_root<Int_t>(TParameter<Int_t> const& in) {
	return in.GetVal();
}
template<>
Long_t param_from_root<Long_t>(TParameter<Long_t> const& in) {
	return in.GetVal();
}
template<>
Real param_from_root<Real>(TParameter<Real> const& in) {
	return in.GetVal();
}
template<>
Nucleus param_from_root<Nucleus>(TParameter<int> const& in) {
	return static_cast<Nucleus>(in.GetVal());
}
template<>
Lepton param_from_root<Lepton>(TParameter<int> const& in) {
	return static_cast<Lepton>(in.GetVal());
}
template<>
Hadron param_from_root<Hadron>(TParameter<int> const& in) {
	return static_cast<Hadron>(in.GetVal());
}
template<>
std::string param_from_root<std::string>(TObjString const& in) {
	return in.GetString().Data();
}
template<>
Vec3 param_from_root<Vec3>(TVector3 const& in) {
	return Vec3(in.X(), in.Y(), in.Z());
}
template<>
Bounds param_from_root<Bounds>(TVector2 const& in) {
	return Bounds(in.X(), in.Y());
}

template<
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
TParameter<T> root_from_param(T const& in) {
	return TParameter<T>("", in);
}
template<
	typename T,
	typename = typename std::enable_if<is_scoped_enum<T>::value>::type>
TParameter<typename std::underlying_type<T>::type> root_from_param(T const& in) {
	return TParameter<typename std::underlying_type<T>::type>(
		"", static_cast<typename std::underlying_type<T>::type>(in));
}
TObjString root_from_param(std::string const& str) {
	return TObjString(str.c_str());
}
TVector3 root_from_param(Vec3 const& vec) {
	return TVector3(vec.x, vec.y, vec.z);
}
TVector2 root_from_param(Bounds const& bounds) {
	return TVector2(bounds.min(), bounds.max());
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

std::ostream& operator<<(std::ostream& os, Bounds const& bounds) {
	os << bounds.min() << " " << bounds.max();
	return os;
}
std::istream& operator>>(std::istream& is, Bounds& bounds) {
	Real min;
	Real max;
	is >> min >> max;
	bounds = Bounds(min, max);
	return is;
}

void Params::write_root(TFile& file) const {
	file.cd();
	WRITE_PARAM_ROOT(file, version_major);
	WRITE_PARAM_ROOT(file, version_minor);
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
	WRITE_PARAM_ROOT(file, mass_threshold);
	WRITE_PARAM_ROOT(file, target_pol);
	WRITE_PARAM_ROOT(file, beam_pol);
	WRITE_PARAM_ROOT(file, k_0_bar);
	WRITE_PARAM_ROOT(file, x_cut);
	WRITE_PARAM_ROOT(file, y_cut);
	WRITE_PARAM_ROOT(file, z_cut);
	WRITE_PARAM_ROOT(file, ph_t_sq_cut);
	WRITE_PARAM_ROOT(file, phi_h_cut);
	WRITE_PARAM_ROOT(file, phi_cut);
}

void Params::read_root(TFile& file) {
	READ_PARAM_ROOT(file, version_major);
	READ_PARAM_ROOT(file, version_minor);
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
	READ_PARAM_ROOT(file, mass_threshold);
	READ_PARAM_ROOT(file, target_pol);
	READ_PARAM_ROOT(file, beam_pol);
	READ_PARAM_ROOT(file, k_0_bar);
	READ_PARAM_ROOT(file, x_cut);
	READ_PARAM_ROOT(file, y_cut);
	READ_PARAM_ROOT(file, z_cut);
	READ_PARAM_ROOT(file, ph_t_sq_cut);
	READ_PARAM_ROOT(file, phi_h_cut);
	READ_PARAM_ROOT(file, phi_cut);
}

void Params::write(std::ostream& file) const {
	// TODO: Set floating point precision so that all digits are kept, then
	// unset after.
	WRITE_PARAM(file, version_major);
	WRITE_PARAM(file, version_minor);
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
	WRITE_PARAM(file, mass_threshold);
	WRITE_PARAM(file, target_pol);
	WRITE_PARAM(file, beam_pol);
	WRITE_PARAM(file, k_0_bar);
	WRITE_PARAM(file, x_cut);
	WRITE_PARAM(file, y_cut);
	WRITE_PARAM(file, z_cut);
	WRITE_PARAM(file, ph_t_sq_cut);
	WRITE_PARAM(file, phi_h_cut);
	WRITE_PARAM(file, phi_cut);
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

	READ_PARAM(params, version_major);
	READ_PARAM(params, version_minor);
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
	READ_PARAM(params, mass_threshold);
	READ_PARAM(params, target_pol);
	READ_PARAM(params, beam_pol);
	READ_PARAM(params, k_0_bar);
	READ_PARAM(params, x_cut);
	READ_PARAM(params, y_cut);
	READ_PARAM(params, z_cut);
	READ_PARAM(params, ph_t_sq_cut);
	READ_PARAM(params, phi_h_cut);
	READ_PARAM(params, phi_cut);

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
	return foam_params.version_major == version_major
		&& foam_params.version_minor <= version_minor
		&& foam_params.num_init >= num_init
		&& (foam_params.seed_init == seed_init || 0 == seed_init)
		&& foam_params.beam_energy == beam_energy
		&& foam_params.beam == beam
		&& foam_params.target == target
		&& foam_params.hadron == hadron
		&& foam_params.mass_threshold == mass_threshold
		&& foam_params.target_pol.x == target_pol.x
		&& foam_params.target_pol.y == target_pol.y
		&& foam_params.target_pol.z == target_pol.z
		&& foam_params.beam_pol == beam_pol
		&& foam_params.k_0_bar == k_0_bar
		&& foam_params.x_cut.min() == x_cut.min()
		&& foam_params.x_cut.max() == x_cut.max()
		&& foam_params.y_cut.min() == y_cut.min()
		&& foam_params.y_cut.max() == y_cut.max()
		&& foam_params.z_cut.min() == z_cut.min()
		&& foam_params.z_cut.max() == z_cut.max()
		&& foam_params.ph_t_sq_cut.min() == ph_t_sq_cut.min()
		&& foam_params.ph_t_sq_cut.max() == ph_t_sq_cut.max()
		&& foam_params.phi_h_cut.min() == phi_h_cut.min()
		&& foam_params.phi_h_cut.max() == phi_h_cut.max()
		&& foam_params.phi_cut.min() == phi_cut.min()
		&& foam_params.phi_cut.max() == phi_cut.max();
}

