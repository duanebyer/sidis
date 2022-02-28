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
#include <TDirectory.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TString.h>
#include <TVector3.h>

#include "exception.hpp"

using namespace sidis;
using namespace sidis::math;

namespace {

template<typename T, typename=void>
struct RootParser {
	static void write_root_dir(TDirectory&, Param<T> const&) {
	}
	static void read_root_dir(TDirectory&, Param<T>&) {
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
void write_param_root_dir(TDirectory& file, Param<T> const& param) {
	if (param.occupied()) {
		RootParser<T>::write_root_dir(file, param);
	}
}
template<typename T>
void read_param_root_dir(TDirectory& file, Param<T>& param) {
	RootParser<T>::read_root_dir(file, param);
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
	static void write_root_dir(TDirectory& file, Param<Version> const& param) {
		Int_t vals[2] = { param->v_major, param->v_minor };
		TArrayI object(2, vals);
		Int_t result = file.WriteObject<TArrayI>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<Version>& param) {
		auto* object = file.Get<TArrayI>(param.name());
		if (object != nullptr && object->GetSize() == 2) {
			param.reset(Version(object->At(0), object->At(1)));
		} else {
			param.reset();
		}
	}
};
// Toggle.
template<>
struct RootParser<Toggle> {
	static void write_root_dir(TDirectory& file, Param<Toggle> const& param) {
		TParameter<Bool_t> object("", *param);
		Int_t result = file.WriteObject<TParameter<Bool_t> >(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<Toggle>& param) {
		auto* object = file.Get<TParameter<Bool_t> >(param.name());
		if (object != nullptr) {
			param.reset(object->GetVal());
		} else {
			param.reset();
		}
	}
};
// Number types.
template<typename T>
struct RootParser<T, typename std::enable_if<std::is_arithmetic<T>::value>::type> {
	static void write_root_dir(TDirectory& file, Param<T> const& param) {
		TParameter<T> object("", *param);
		Int_t result = file.WriteObject<TParameter<T> >(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<T>& param) {
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
	static void write_root_dir(TDirectory& file, Param<T> const& param) {
		TParameter<Integral> object("", static_cast<Integral>(*param));
		Int_t result = file.WriteObject<decltype(object)>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<T>& param) {
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
	static void write_root_dir(TDirectory& file, Param<std::string> const& param) {
		TObjString object(param->c_str());
		Int_t result = file.WriteObject<TObjString>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<std::string>& param) {
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
	static void write_root_dir(TDirectory& file, Param<Vec3> const& param) {
		TVector3 object(param->x, param->y, param->z);
		Int_t result = file.WriteObject<TVector3>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<Vec3>& param) {
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
	static void write_root_dir(TDirectory& file, Param<Bound> const& param) {
		Double_t vals[2] = { param->min(), param->max() };
		TArrayD object(2, vals);
		Int_t result = file.WriteObject<TArrayD>(&object, param.name());
		if (result == 0) {
			throw std::runtime_error(
				std::string("Could not write parameter '")
				+ param.name() + "' to ROOT file.");
		}
	}
	static void read_root_dir(TDirectory& file, Param<Bound>& param) {
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
// Toggle.
std::ostream& operator<<(std::ostream& os, Toggle const& toggle) {
	return os << (toggle ? "on" : "off");
}
std::istream& operator>>(std::istream& is, Toggle& toggle) {
	std::string str;
	is >> str;
	if (str == "1" || str == "on" || str == "true" || str == "yes") {
		toggle.on = true;
	} else if (str == "0" || str == "off" || str == "false" || str == "no") {
		toggle.on = false;
	} else {
		is.setstate(std::ios_base::failbit);
	}
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

// Common part to `compatible_with_foam` and `compatible_with_event`.
void compatible_common(Params const& p1, Params const& p2) {
	if (p1.version->v_major != p2.version->v_major) {
		throw std::runtime_error("Incompatible versions.");
	} else if (*p1.rc_method != *p2.rc_method) {
		throw std::runtime_error("Different rc methods.");
	} else if (*p1.sf_set != *p2.sf_set) {
		throw std::runtime_error("Different sf set.");
	} else if (*p1.beam_energy != *p2.beam_energy) {
		throw std::runtime_error("Different beam energies.");
	} else if (*p1.beam != *p2.beam) {
		throw std::runtime_error("Different beam lepton type.");
	} else if (*p1.target != *p2.target) {
		throw std::runtime_error("Different target nucleus type.");
	} else if (*p1.hadron != *p2.hadron) {
		throw std::runtime_error("Different leading hadron type.");
	} else if (*p1.Mth != *p2.Mth) {
		throw std::runtime_error("Different mass thresholds.");
	} else if (*p1.target_pol != *p2.target_pol) {
		throw std::runtime_error("Different target polarizations.");
	} else if (*p1.beam_pol != *p2.beam_pol) {
		throw std::runtime_error("Different beam polarizations.");
	} else if ((*p1.rc_method == RcMethod::APPROX || *p1.rc_method == RcMethod::EXACT)
			&& *p1.k_0_bar != *p2.k_0_bar) {
		throw std::runtime_error("Different soft photon cutoffs.");
	} else if (p1.cut_x.get_or(Bound::UNIT) != p2.cut_x.get_or(Bound::UNIT)) {
		throw std::runtime_error("Different cuts on x.");
	} else if (p1.cut_y.get_or(Bound::UNIT) != p2.cut_y.get_or(Bound::UNIT)) {
		throw std::runtime_error("Different cuts on y.");
	} else if (p1.cut_z.get_or(Bound::UNIT) != p2.cut_z.get_or(Bound::UNIT)) {
		throw std::runtime_error("Different cuts on z.");
	} else if (p1.cut_ph_t_sq != p2.cut_ph_t_sq) {
		throw std::runtime_error("Different cuts on ph_t².");
	} else if (p1.cut_phi_h != p2.cut_phi_h) {
		throw std::runtime_error("Different cuts on φ_h.");
	} else if (p1.cut_phi != p2.cut_phi) {
		throw std::runtime_error("Different cuts on φ.");
	} else if (p1.cut_Q_sq != p2.cut_Q_sq) {
		throw std::runtime_error("Different cuts on Q².");
	} else if (p1.cut_t != p2.cut_t) {
		throw std::runtime_error("Different cuts on t.");
	} else if (p1.cut_W_sq != p2.cut_W_sq) {
		throw std::runtime_error("Different cuts on W².");
	} else if (p1.cut_r != p2.cut_r) {
		throw std::runtime_error("Different cuts on r.");
	} else if (p1.cut_mx_sq != p2.cut_mx_sq) {
		throw std::runtime_error("Different cuts on mx².");
	} else if (p1.cut_q_0 != p2.cut_q_0) {
		throw std::runtime_error("Different cuts on q_0.");
	} else if (p1.cut_k2_0 != p2.cut_k2_0) {
		throw std::runtime_error("Different cuts on k2_0.");
	} else if (p1.cut_ph_0 != p2.cut_ph_0) {
		throw std::runtime_error("Different cuts on ph_0.");
	} else if (p1.cut_theta_q != p2.cut_theta_q) {
		throw std::runtime_error("Different cuts on θ_q.");
	} else if (p1.cut_theta_k2 != p2.cut_theta_k2) {
		throw std::runtime_error("Different cuts on θ_k2.");
	} else if (p1.cut_theta_h != p2.cut_theta_h) {
		throw std::runtime_error("Different cuts on θ_h.");
	}
}

void compatible_common_rad(Params const& p1, Params const& p2) {
	if (*p1.rad_gen && *p2.rad_gen) {
		if (p1.cut_tau != p2.cut_tau) {
			throw std::runtime_error("Different cuts on τ.");
		} else if (p1.cut_phi_k != p2.cut_phi_k) {
			throw std::runtime_error("Different cuts on φ_k.");
		} else if (p1.cut_R != p2.cut_R) {
			throw std::runtime_error("Different cuts on R.");
		} else if (p1.cut_k_0_bar != p2.cut_k_0_bar) {
			throw std::runtime_error(
				std::string("Different '") + p1.cut_k_0_bar.name() + "'.");
		} else if (p1.cut_k_0 != p2.cut_k_0) {
			throw std::runtime_error(
				std::string("Different '") + p1.cut_k_0.name() + "'.");
		} else if (p1.cut_theta_k != p2.cut_theta_k) {
			throw std::runtime_error("Different cuts on θ_k.");
		}
	}
}

}

void Params::write_root(TFile& file) const {
	TDirectory* dir = file.mkdir("params", "params");
	if (dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create directory 'params' in ROOT file."));
	}
	dir->cd();
	write_param_root_dir(*dir, version);
	write_param_root_dir(*dir, strict);
	write_param_root_dir(*dir, event_file);
	write_param_root_dir(*dir, rc_method);
	write_param_root_dir(*dir, nrad_gen);
	write_param_root_dir(*dir, nrad_max_cells);
	write_param_root_dir(*dir, nrad_target_eff);
	write_param_root_dir(*dir, nrad_scale_exp);
	write_param_root_dir(*dir, rad_gen);
	write_param_root_dir(*dir, rad_max_cells);
	write_param_root_dir(*dir, rad_target_eff);
	write_param_root_dir(*dir, rad_scale_exp);
	write_param_root_dir(*dir, write_momenta);
	write_param_root_dir(*dir, write_photon);
	write_param_root_dir(*dir, foam_file);
	write_param_root_dir(*dir, sf_set);
	write_param_root_dir(*dir, num_events);
	write_param_root_dir(*dir, rej_weight);
	write_param_root_dir(*dir, seed);
	write_param_root_dir(*dir, nrad_seed_init);
	write_param_root_dir(*dir, rad_seed_init);
	write_param_root_dir(*dir, beam_energy);
	write_param_root_dir(*dir, beam);
	write_param_root_dir(*dir, target);
	write_param_root_dir(*dir, hadron);
	write_param_root_dir(*dir, Mth);
	write_param_root_dir(*dir, target_pol);
	write_param_root_dir(*dir, beam_pol);
	write_param_root_dir(*dir, k_0_bar);
	write_param_root_dir(*dir, cut_x);
	write_param_root_dir(*dir, cut_y);
	write_param_root_dir(*dir, cut_z);
	write_param_root_dir(*dir, cut_ph_t_sq);
	write_param_root_dir(*dir, cut_phi_h);
	write_param_root_dir(*dir, cut_phi);
	write_param_root_dir(*dir, cut_Q_sq);
	write_param_root_dir(*dir, cut_t);
	write_param_root_dir(*dir, cut_W_sq);
	write_param_root_dir(*dir, cut_r);
	write_param_root_dir(*dir, cut_mx_sq);
	write_param_root_dir(*dir, cut_q_0);
	write_param_root_dir(*dir, cut_k2_0);
	write_param_root_dir(*dir, cut_ph_0);
	write_param_root_dir(*dir, cut_theta_q);
	write_param_root_dir(*dir, cut_theta_k2);
	write_param_root_dir(*dir, cut_theta_h);
	write_param_root_dir(*dir, cut_tau);
	write_param_root_dir(*dir, cut_phi_k);
	write_param_root_dir(*dir, cut_R);
	write_param_root_dir(*dir, cut_k_0_bar);
	write_param_root_dir(*dir, cut_k_0);
	write_param_root_dir(*dir, cut_theta_k);
	file.cd();
}

void Params::read_root(TFile& file) {
	TDirectory* dir = file.Get<TDirectory>("params");
	if (dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Couldn't find directory 'params' in ROOT file."));
	}
	read_param_root_dir(*dir, version);
	if (version->v_major != SIDIS_PARAMS_VERSION_MAJOR
			|| version->v_minor > SIDIS_PARAMS_VERSION_MINOR) {
		throw std::runtime_error(
			std::string("Cannot read parameters from version ")
			+ std::to_string(version->v_major) + "."
			+ std::to_string(version->v_minor) + ".");
	}
	read_param_root_dir(*dir, strict);
	read_param_root_dir(*dir, event_file);
	read_param_root_dir(*dir, rc_method);
	read_param_root_dir(*dir, nrad_gen);
	read_param_root_dir(*dir, nrad_max_cells);
	read_param_root_dir(*dir, nrad_target_eff);
	read_param_root_dir(*dir, nrad_scale_exp);
	read_param_root_dir(*dir, rad_gen);
	read_param_root_dir(*dir, rad_max_cells);
	read_param_root_dir(*dir, rad_target_eff);
	read_param_root_dir(*dir, rad_scale_exp);
	read_param_root_dir(*dir, write_momenta);
	read_param_root_dir(*dir, write_photon);
	read_param_root_dir(*dir, foam_file);
	read_param_root_dir(*dir, sf_set);
	read_param_root_dir(*dir, num_events);
	read_param_root_dir(*dir, rej_weight);
	read_param_root_dir(*dir, seed);
	read_param_root_dir(*dir, nrad_seed_init);
	read_param_root_dir(*dir, rad_seed_init);
	read_param_root_dir(*dir, beam_energy);
	read_param_root_dir(*dir, beam);
	read_param_root_dir(*dir, target);
	read_param_root_dir(*dir, hadron);
	read_param_root_dir(*dir, Mth);
	read_param_root_dir(*dir, target_pol);
	read_param_root_dir(*dir, beam_pol);
	read_param_root_dir(*dir, k_0_bar);
	read_param_root_dir(*dir, cut_x);
	read_param_root_dir(*dir, cut_y);
	read_param_root_dir(*dir, cut_z);
	read_param_root_dir(*dir, cut_ph_t_sq);
	read_param_root_dir(*dir, cut_phi_h);
	read_param_root_dir(*dir, cut_phi);
	read_param_root_dir(*dir, cut_Q_sq);
	read_param_root_dir(*dir, cut_t);
	read_param_root_dir(*dir, cut_W_sq);
	read_param_root_dir(*dir, cut_r);
	read_param_root_dir(*dir, cut_mx_sq);
	read_param_root_dir(*dir, cut_q_0);
	read_param_root_dir(*dir, cut_k2_0);
	read_param_root_dir(*dir, cut_ph_0);
	read_param_root_dir(*dir, cut_theta_q);
	read_param_root_dir(*dir, cut_theta_k2);
	read_param_root_dir(*dir, cut_theta_h);
	read_param_root_dir(*dir, cut_tau);
	read_param_root_dir(*dir, cut_phi_k);
	read_param_root_dir(*dir, cut_R);
	read_param_root_dir(*dir, cut_k_0_bar);
	read_param_root_dir(*dir, cut_k_0);
	read_param_root_dir(*dir, cut_theta_k);
}

void Params::write_stream(std::ostream& file) const {
	write_param_stream(file, version, true);
	write_param_stream(file, strict);
	write_param_stream(file, event_file);
	write_param_stream(file, rc_method);
	write_param_stream(file, nrad_gen);
	write_param_stream(file, nrad_max_cells);
	write_param_stream(file, nrad_target_eff);
	write_param_stream(file, nrad_scale_exp);
	write_param_stream(file, rad_gen);
	write_param_stream(file, rad_max_cells);
	write_param_stream(file, rad_target_eff);
	write_param_stream(file, rad_scale_exp);
	write_param_stream(file, write_momenta);
	write_param_stream(file, write_photon);
	write_param_stream(file, foam_file);
	write_param_stream(file, sf_set);
	write_param_stream(file, num_events);
	write_param_stream(file, rej_weight);
	write_param_stream(file, seed);
	write_param_stream(file, nrad_seed_init);
	write_param_stream(file, rad_seed_init);
	write_param_stream(file, beam_energy);
	write_param_stream(file, beam);
	write_param_stream(file, target);
	write_param_stream(file, hadron);
	write_param_stream(file, Mth);
	write_param_stream(file, target_pol);
	write_param_stream(file, beam_pol);
	write_param_stream(file, k_0_bar);
	write_param_stream(file, cut_x);
	write_param_stream(file, cut_y);
	write_param_stream(file, cut_z);
	write_param_stream(file, cut_ph_t_sq);
	write_param_stream(file, cut_phi_h);
	write_param_stream(file, cut_phi);
	write_param_stream(file, cut_Q_sq);
	write_param_stream(file, cut_t);
	write_param_stream(file, cut_W_sq);
	write_param_stream(file, cut_r);
	write_param_stream(file, cut_mx_sq);
	write_param_stream(file, cut_q_0);
	write_param_stream(file, cut_k2_0);
	write_param_stream(file, cut_ph_0);
	write_param_stream(file, cut_theta_q);
	write_param_stream(file, cut_theta_k2);
	write_param_stream(file, cut_theta_h);
	write_param_stream(file, cut_tau);
	write_param_stream(file, cut_phi_k);
	write_param_stream(file, cut_R);
	write_param_stream(file, cut_k_0_bar);
	write_param_stream(file, cut_k_0);
	write_param_stream(file, cut_theta_k);
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
	consume_param_from_map(map, strict);
	consume_param_from_map(map, event_file);
	consume_param_from_map(map, rc_method);
	consume_param_from_map(map, nrad_gen);
	consume_param_from_map(map, nrad_max_cells);
	consume_param_from_map(map, nrad_target_eff);
	consume_param_from_map(map, nrad_scale_exp);
	consume_param_from_map(map, rad_gen);
	consume_param_from_map(map, rad_max_cells);
	consume_param_from_map(map, rad_target_eff);
	consume_param_from_map(map, rad_scale_exp);
	consume_param_from_map(map, write_momenta);
	consume_param_from_map(map, write_photon);
	consume_param_from_map(map, foam_file);
	consume_param_from_map(map, sf_set);
	consume_param_from_map(map, num_events);
	consume_param_from_map(map, rej_weight);
	consume_param_from_map(map, seed);
	consume_param_from_map(map, nrad_seed_init);
	consume_param_from_map(map, rad_seed_init);
	consume_param_from_map(map, beam_energy);
	consume_param_from_map(map, beam);
	consume_param_from_map(map, target);
	consume_param_from_map(map, hadron);
	consume_param_from_map(map, Mth);
	consume_param_from_map(map, target_pol);
	consume_param_from_map(map, beam_pol);
	consume_param_from_map(map, k_0_bar);
	consume_param_from_map(map, cut_x);
	consume_param_from_map(map, cut_y);
	consume_param_from_map(map, cut_z);
	consume_param_from_map(map, cut_ph_t_sq);
	consume_param_from_map(map, cut_phi_h);
	consume_param_from_map(map, cut_phi);
	consume_param_from_map(map, cut_Q_sq);
	consume_param_from_map(map, cut_t);
	consume_param_from_map(map, cut_W_sq);
	consume_param_from_map(map, cut_r);
	consume_param_from_map(map, cut_mx_sq);
	consume_param_from_map(map, cut_q_0);
	consume_param_from_map(map, cut_k2_0);
	consume_param_from_map(map, cut_ph_0);
	consume_param_from_map(map, cut_theta_q);
	consume_param_from_map(map, cut_theta_k2);
	consume_param_from_map(map, cut_theta_h);
	consume_param_from_map(map, cut_tau);
	consume_param_from_map(map, cut_phi_k);
	consume_param_from_map(map, cut_R);
	consume_param_from_map(map, cut_k_0_bar);
	consume_param_from_map(map, cut_k_0);
	consume_param_from_map(map, cut_theta_k);

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

void Params::fill_defaults() {
	// This function is somewhat complicated because it tries to consider many
	// possible combinations of options and sort through them to produce a
	// "normalized" form of any parameter file.
	version.get_or_insert(Version());
	strict.get_or_insert(true);
	if (!num_events.occupied()) {
		throw std::runtime_error(
			std::string("Must specify number of events to generate (")
			+ num_events.name() + ").");
	}
	if (!beam_energy.occupied()) {
		throw std::runtime_error(
			std::string("Must specify beam energy (")
			+ beam_energy.name() + ").");
	}
	if (!beam.occupied()) {
		throw std::runtime_error(
			std::string("Must specify beam lepton flavor (")
			+ beam.name() + ").");
	}
	if (!target.occupied()) {
		throw std::runtime_error(
			std::string("Must specify target nucleus (")
			+ target.name() + ").");
	}
	if (!hadron.occupied()) {
		throw std::runtime_error(
			std::string("Must specify leading hadron (")
			+ hadron.name() + ").");
	}
	if (!Mth.occupied()) {
		throw std::runtime_error(
			std::string("Must specify mass threshold (")
			+ Mth.name() + ").");
	}
	if (!sf_set.occupied()) {
		throw std::runtime_error(
			std::string("Must specify structure function set (")
			+ sf_set.name() + ").");
	}
	event_file.get_or_insert("gen.root");
	rc_method.get_or_insert(RcMethod::APPROX);
	rej_weight.get_or_insert(1.1);
	seed.get_or_insert(0);
	target_pol.get_or_insert(Vec3::ZERO);
	beam_pol.get_or_insert(0.);
	write_momenta.get_or_insert(true);
	// RC methods.
	// The process for determining whether radiative and non-radiative events
	// are generated looks something like this:
	//  * If RC are disabled, then no radiative events are generated.
	//  * If `cut_k_0_bar` is above the soft threshold, then only radiative
	//    events are generated.
	//  * `cut_k_0_bar` cannot be partway between zero and the soft threshold,
	//    because it is impossible to generate events according to that.
	if (*rc_method == RcMethod::APPROX || *rc_method == RcMethod::EXACT) {
		rad_gen.get_or_insert(true);
		k_0_bar.get_or_insert(0.01);
		if (nrad_gen.occupied() && !cut_k_0_bar.occupied()) {
			if (*nrad_gen) {
				cut_k_0_bar.reset(Bound::POSITIVE);
			} else {
				cut_k_0_bar.reset(Bound::POSITIVE + *k_0_bar);
			}
		} else {
			cut_k_0_bar.get_or_insert(Bound::POSITIVE);
		}
	} else {
		if (rad_gen.occupied() && *rad_gen) {
			if (*strict) {
				throw std::runtime_error(
					std::string("Cannot generate radiative events without "
						"radiative corrections enabled (")
					+ rad_gen.name() + ", " + rc_method.name() + ").");
			}
		}
		rad_gen.reset(false);
		if (k_0_bar.occupied()) {
			if (*strict) {
				throw std::runtime_error(
					std::string("Cannot provide soft photon threshold without "
						"radiative corrections enabled (")
					+ k_0_bar.name() + ", " + rc_method.name() + ").");
			}
			k_0_bar.reset();
		}
	}
	if (*rad_gen) {
		if (cut_k_0_bar->min() <= 0.) {
			nrad_gen.get_or_insert(true);
			if (!*nrad_gen) {
				if (*strict) {
					throw std::runtime_error(
						std::string("Cannot generate radiative events only "
							"with photon energy minimum <= 0 (")
						+ rad_gen.name() + ", " + nrad_gen.name()
						+ ", " + cut_k_0_bar.name() + ").");
				}
			}
		} else if (*k_0_bar <= cut_k_0_bar->min()) {
			if (nrad_gen.occupied() && *nrad_gen) {
				if (*strict) {
					throw std::runtime_error(
						std::string("Cannot generate non-radiative events with "
							"photon energy minimum > soft threshold (")
						+ nrad_gen.name() + ", " + k_0_bar.name() + ", "
						+ cut_k_0_bar.name() + ").");
				}
			}
			nrad_gen.reset(false);
		} else {
			throw std::runtime_error(
				std::string("Cannot have photon energy cut minimum < soft "
					"threshold without generating non-radiative events (")
				+ nrad_gen.name() + ", " + k_0_bar.name() + ", "
				+ cut_k_0_bar.name() + ").");
		}
	} else {
		nrad_gen.get_or_insert(true);
	}
	if (!*nrad_gen && !*rad_gen) {
		throw std::runtime_error(
			std::string("Cannot disable all event types (")
			+ nrad_gen.name() + ", " + rad_gen.name() + ").");
	}
	// Basic options associated with radiative and non-radiative events in
	// particular.
	if (*rad_gen) {
		if (write_momenta.occupied() && *write_momenta) {
			write_photon.get_or_insert(true);
		} else {
			if (write_photon.occupied() && *write_photon) {
				if (*strict) {
					throw std::runtime_error(
						std::string("Cannot enable '") + write_photon.name()
						+ "' when '" + write_momenta.name() + "' is disabled.");
				} else {
					write_photon.reset(false);
				}
			}
		}
		rad_max_cells.get_or_insert(1048576);
		rad_target_eff.get_or_insert(0.65);
		rad_scale_exp.get_or_insert(0.20);
		rad_seed_init.get_or_insert(0);
	} else {
		if (write_photon.occupied() && *write_photon) {
			if (*strict) {
				throw std::runtime_error(
					std::string("Cannot enable '") + write_photon.name()
					+ "' when no radiative events are being generated.");
			} else {
				write_photon.reset();
			}
		}
		if (*strict && rad_max_cells.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + rad_max_cells.name()
				+ "' when no radiative events are being generated.");
		}
		if (*strict && rad_target_eff.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + rad_target_eff.name()
				+ "' when no radiative events are being generated.");
		}
		if (*strict && rad_scale_exp.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + rad_scale_exp.name()
				+ "' when no radiative events are being generated.");
		}
		if (*strict && rad_seed_init.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + rad_seed_init.name()
				+ "' when no radiative events are being generated.");
		}
	}
	if (*nrad_gen) {
		nrad_max_cells.get_or_insert(1048576);
		nrad_target_eff.get_or_insert(0.95);
		nrad_scale_exp.get_or_insert(0.60);
		nrad_seed_init.get_or_insert(0);
	} else {
		if (*strict && nrad_max_cells.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + nrad_max_cells.name()
				+ "' when no non-radiative events are being generated.");
		}
		if (*strict && nrad_target_eff.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + nrad_target_eff.name()
				+ "' when no non-radiative events are being generated.");
		}
		if (*strict && nrad_scale_exp.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + nrad_scale_exp.name()
				+ "' when no non-radiative events are being generated.");
		}
		if (*strict && nrad_seed_init.occupied()) {
			throw std::runtime_error(
				std::string("Cannot specify '") + nrad_seed_init.name()
				+ "' when no non-radiative events are being generated.");
		}
	}
	// Cuts.
	if (*rad_gen && *nrad_gen) {
		if (cut_tau.occupied()
				|| cut_phi_k.occupied()
				|| cut_R.occupied()
				|| cut_k_0.occupied()
				|| cut_theta_k.occupied()) {
			throw std::runtime_error(
				"Cannot apply radiative cuts to non-radiative events.");
		}
	}
	// Verify that cuts make sense. This isn't comprehensive, but is primarily
	// important to avoid cuts on the azimuthal angles larger than 360 degrees.
	if (*strict && cut_x.occupied() && !Bound::UNIT.contains(*cut_x)) {
		throw std::runtime_error(
			std::string("Cut on x must lie between 0 and 1 (")
			+ cut_x.name() + ").");
	}
	if (*strict && cut_y.occupied() && !Bound::UNIT.contains(*cut_y)) {
		throw std::runtime_error(
			std::string("Cut on y must lie between 0 and 1 (")
			+ cut_y.name() + ").");
	}
	if (*strict && cut_z.occupied() && !Bound::UNIT.contains(*cut_z)) {
		throw std::runtime_error(
			std::string("Cut on z must lie between 0 and 1 (")
			+ cut_z.name() + ").");
	}
	if (cut_phi_h.occupied() && cut_phi_h->size() >= 360.) {
		throw std::runtime_error(
			std::string("Cut on φ_h must be smaller than 360 degrees (")
			+ cut_phi_h.name() + ").");
	}
	if (cut_phi.occupied() && cut_phi->size() >= 360.) {
		throw std::runtime_error(
			std::string("Cut on φ must be smaller than 360 degrees (")
			+ cut_phi.name() + ").");
	}
	if (cut_phi_k.occupied() && cut_phi_k->size() >= 360.) {
		throw std::runtime_error(
			std::string("Cut on φ_k must be smaller than 360 degrees (")
			+ cut_phi_k.name() + ").");
	}
	if (*strict && cut_theta_q.occupied() && !Bound(0., 180.).contains(*cut_theta_q)) {
		throw std::runtime_error(
			std::string("Cut on θ_q must lie between 0 and 180 degrees (")
			+ cut_theta_q.name() + ").");
	}
	if (*strict && cut_theta_h.occupied() && !Bound(0., 180.).contains(*cut_theta_h)) {
		throw std::runtime_error(
			std::string("Cut on θ_h must lie between 0 and 180 degrees (")
			+ cut_theta_h.name() + ").");
	}
	if (*strict && cut_theta_k.occupied() && !Bound(0., 180.).contains(*cut_theta_k)) {
		throw std::runtime_error(
			std::string("Cut on θ_k must lie between 0 and 180 degrees (")
			+ cut_theta_k.name() + ").");
	}
}

void Params::compatible_with_foam(EventType type, Params const& foam_params) const {
	if (strict && !foam_params.strict) {
		throw std::runtime_error("Incompatible strictness levels.");
	} else if (type == EventType::NRAD && *nrad_gen && !*foam_params.nrad_gen) {
		throw std::runtime_error("No non-radiative FOAM available.");
	} else if (type == EventType::RAD && *rad_gen && !*foam_params.rad_gen) {
		throw std::runtime_error("No radiative FOAM available.");
	} else if (type == EventType::NRAD && *nrad_seed_init != 0 && (*foam_params.nrad_seed_init != *nrad_seed_init || *foam_params.nrad_max_cells != *nrad_max_cells || *foam_params.nrad_target_eff != *nrad_target_eff || *foam_params.nrad_scale_exp != *nrad_scale_exp)) {
		throw std::runtime_error("Different non-radiative initialization seed/conditions.");
	} else if (type == EventType::RAD && *rad_seed_init != 0 && (*foam_params.rad_seed_init != *rad_seed_init || *foam_params.rad_max_cells != *rad_max_cells || *foam_params.rad_target_eff != *rad_target_eff || *foam_params.rad_scale_exp != *rad_scale_exp)) {
		throw std::runtime_error("Different radiative initializalation seed/conditions.");
	} else if (type == EventType::NRAD && *nrad_gen && !*foam_params.nrad_gen) {
		throw std::runtime_error("FOAM doesn't support non-radiative events.");
	} else if (type == EventType::RAD && *rad_gen && !*foam_params.rad_gen) {
		throw std::runtime_error("FOAM doesn't support radiative events.");
	} else if (type == EventType::NRAD && *nrad_max_cells > *foam_params.nrad_max_cells) {
		throw std::runtime_error("Non-radiative FOAM doesn't have sufficient max cells.");
	} else if (type == EventType::RAD && *rad_max_cells > *foam_params.rad_max_cells) {
		throw std::runtime_error("Radiative FOAM doesn't have sufficient max cells.");
	} else if (type == EventType::NRAD && *nrad_target_eff > *foam_params.nrad_target_eff) {
		throw std::runtime_error("Non-radiative FOAM doesn't have sufficient target efficiency.");
	} else if (type == EventType::RAD && *rad_target_eff > *foam_params.rad_target_eff) {
		throw std::runtime_error("Radiative FOAM doesn't have sufficient target efficiency.");
	} else if (type == EventType::NRAD && *nrad_scale_exp != *foam_params.nrad_scale_exp) {
		throw std::runtime_error("Non-radiative FOAM has different scaling exponent.");
	} else if (type == EventType::RAD && *rad_scale_exp != *foam_params.rad_scale_exp) {
		throw std::runtime_error("Radiative FOAM has different scaling exponent.");
	}
	compatible_common(*this, foam_params);
	if (type == EventType::RAD) {
		compatible_common_rad(*this, foam_params);
	}
}

void Params::compatible_with_event(Params const& params) const {
	if (*seed == *params.seed && *seed != 0) {
		throw std::runtime_error("Can't merge parameters with same seed.");
	}
	if (*rej_weight != *params.rej_weight) {
		throw std::runtime_error("Can't merge parameters with different rejection weights.");
	}
	compatible_common(*this, params);
	compatible_common_rad(*this, params);
}

void Params::merge(Params const& params) {
	compatible_with_event(params);
	version->v_minor = std::max(version->v_minor, params.version->v_minor);
	*strict = *strict && *params.strict;
	event_file.reset("<undefined>");
	*nrad_gen = *nrad_gen || *params.nrad_gen;
	*rad_gen = *rad_gen || *params.rad_gen;
	*write_momenta = *write_momenta && *params.write_momenta;
	if (*write_momenta) {
		*write_photon = *write_photon && *params.write_photon;
	}
	if (foam_file != params.foam_file) {
		foam_file.reset("<undefined>");
	}
	*num_events += *params.num_events;
	if (rej_weight != params.rej_weight) {
		rej_weight.reset(0.);
	}
	*seed = 0;
	if (*nrad_gen && nrad_seed_init != params.nrad_seed_init) {
		nrad_seed_init.reset(0.);
	}
	if (*rad_gen && rad_seed_init != params.rad_seed_init) {
		rad_seed_init.reset(0.);
	}
	if (*rad_gen ^ *params.rad_gen) {
		cut_tau.take_first(params.cut_tau);
		cut_phi_k.take_first(params.cut_phi_k);
		cut_R.take_first(params.cut_R);
		cut_k_0_bar.take_first(params.cut_k_0_bar);
		cut_k_0.take_first(params.cut_k_0);
		cut_theta_k.take_first(params.cut_theta_k);
	}
}

bool Params::operator==(Params const& rhs) const {
	return version == rhs.version
		&& strict == rhs.strict
		&& event_file == rhs.event_file
		&& rc_method == rhs.rc_method
		&& nrad_gen == rhs.nrad_gen
		&& nrad_max_cells == rhs.nrad_max_cells
		&& nrad_target_eff == rhs.nrad_target_eff
		&& nrad_scale_exp == rhs.nrad_scale_exp
		&& rad_gen == rhs.rad_gen
		&& rad_max_cells == rhs.rad_max_cells
		&& rad_target_eff == rhs.rad_target_eff
		&& rad_scale_exp == rhs.rad_scale_exp
		&& write_momenta == rhs.write_momenta
		&& write_photon == rhs.write_photon
		&& foam_file == rhs.foam_file
		&& sf_set == rhs.sf_set
		&& num_events == rhs.num_events
		&& rej_weight == rhs.rej_weight
		&& seed == rhs.seed
		&& rad_seed_init == rhs.rad_seed_init
		&& nrad_seed_init == rhs.nrad_seed_init
		&& beam_energy == rhs.beam_energy
		&& beam == rhs.beam
		&& target == rhs.target
		&& hadron == rhs.hadron
		&& Mth == rhs.Mth
		&& target_pol == rhs.target_pol
		&& beam_pol == rhs.beam_pol
		&& k_0_bar == rhs.k_0_bar
		&& cut_x == rhs.cut_x
		&& cut_y == rhs.cut_y
		&& cut_z == rhs.cut_z
		&& cut_ph_t_sq == rhs.cut_ph_t_sq
		&& cut_phi_h == rhs.cut_phi_h
		&& cut_phi == rhs.cut_phi
		&& cut_Q_sq == rhs.cut_Q_sq
		&& cut_t == rhs.cut_t
		&& cut_W_sq == rhs.cut_W_sq
		&& cut_r == rhs.cut_r
		&& cut_mx_sq == rhs.cut_mx_sq
		&& cut_q_0 == rhs.cut_q_0
		&& cut_k2_0 == rhs.cut_k2_0
		&& cut_ph_0 == rhs.cut_ph_0
		&& cut_theta_q == rhs.cut_theta_q
		&& cut_theta_k2 == rhs.cut_theta_k2
		&& cut_theta_h == rhs.cut_theta_h
		&& cut_tau == rhs.cut_tau
		&& cut_phi_k == rhs.cut_phi_k
		&& cut_R == rhs.cut_R
		&& cut_k_0_bar == rhs.cut_k_0_bar
		&& cut_k_0 == rhs.cut_k_0
		&& cut_theta_k == rhs.cut_theta_k;
}

