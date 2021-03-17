#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <regex>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <TArrayD.h>
#include <TBranch.h>
#include <TClass.h>
#include <TError.h>
#include <TFile.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>

#include <sidis/sidis.hpp>
#include <sidis/extra/math.hpp>
#include <sidis/sf_set/mask.hpp>
#include <sidis/sf_set/prokudin.hpp>
#include <sidis/sf_set/test.hpp>

#include "params.hpp"
#include "utility.hpp"

using namespace sidis;

namespace {

int const SUCCESS = 0;
int const ERROR_ARG_PARSE = -1;
int const ERROR_FILE_NOT_FOUND = -2;
int const ERROR_FILE_NOT_CREATED = -3;
int const ERROR_PARAMS_PARSE = -4;
int const ERROR_PARAMS_INVALID = -5;
int const ERROR_FOAM_INCOMPATIBLE = -6;
int const ERROR_FOAM_NOT_FOUND = -7;
int const ERROR_STRUCTURE_FUNCTIONS_NOT_FOUND = -8;
int const ERROR_STRUCTURE_FUNCTIONS_PARSE = -9;

// In the future, more types of events (such as exclusive) may be included.
enum class EventType {
	NRAD,
	RAD,
};

// All of the relevant information about one kind of event.
struct EventStats {
	EventType type;
	std::unique_ptr<TFile> foam_file;
	std::unique_ptr<TFoamIntegrand> rho;
	TFoam* foam;
	Double_t xs;
	Double_t xs_err;
	ULong_t num_events;
};

// Converts between the `sidis` 4-vector type and the `ROOT` 4-vector type.
TLorentzVector convert_vec4(math::Vec4 v) {
	return TLorentzVector(v.x, v.y, v.z, v.t);
}

// Fills out cuts from parameters.
void cuts(Params params, cut::Cut* cut_out, cut::CutRad* cut_rad_out) {
	if (cut_out != nullptr) {
		*cut_out = cut::Cut();
		cut_out->x = params.x_cut.get_or(math::Bound::INVALID);
		cut_out->y = params.y_cut.get_or(math::Bound::INVALID);
		cut_out->z = params.z_cut.get_or(math::Bound::INVALID);
		cut_out->ph_t_sq = params.ph_t_sq_cut.get_or(math::Bound::INVALID);
		cut_out->phi_h = params.phi_h_cut.get_or(math::Bound::INVALID);
		cut_out->phi = params.phi_cut.get_or(math::Bound::INVALID);
		cut_out->Q_sq = params.Q_sq_cut.get_or(math::Bound::INVALID);
		cut_out->t = params.t_cut.get_or(math::Bound::INVALID);
		cut_out->w = params.w_cut.get_or(math::Bound::INVALID);
		cut_out->r = params.r_cut.get_or(math::Bound::INVALID);
		cut_out->mx_sq = params.mx_sq_cut.get_or(math::Bound::INVALID);
		cut_out->q_0 = params.q_0_cut.get_or(math::Bound::INVALID);
		cut_out->k2_0 = params.k2_0_cut.get_or(math::Bound::INVALID);
		cut_out->ph_0 = params.ph_0_cut.get_or(math::Bound::INVALID);
		cut_out->theta_q = params.theta_q_cut.get_or(math::Bound::INVALID);
		cut_out->theta_k2 = params.theta_k2_cut.get_or(math::Bound::INVALID);
		cut_out->theta_h = params.theta_h_cut.get_or(math::Bound::INVALID);
	}
	if (cut_rad_out != nullptr) {
		*cut_rad_out = cut::CutRad();
		if (*params.gen_rad) {
			cut_rad_out->tau = params.tau_cut.get_or(math::Bound::INVALID);
			cut_rad_out->phi_k = params.phi_k_cut.get_or(math::Bound::INVALID);
			// The `k_0_bar` cut is mandatory.
			cut_rad_out->k_0_bar = *params.k_0_bar_cut;
			cut_rad_out->k_0 = params.k_0_cut.get_or(math::Bound::INVALID);
			cut_rad_out->theta_k = params.theta_k_cut.get_or(math::Bound::INVALID);
		}
	}
}

// Logical union of two boolean arrays.
void zip_and(bool* begin_1, bool* end_1, bool const* begin_2) {
	bool const* it_2 = begin_2;
	for (bool* it_1 = begin_1; it_1 != end_1; ++it_1) {
		*it_1 &= *it_2;
		++it_2;
	}
}

// Allocates the memory for the structure functions.
int alloc_sf(
		Params params,
		std::unique_ptr<sf::SfSet>* sf_out,
		std::unique_ptr<sf::TmdSet>* tmd_out) {
	// Structure function description comes as a series of words. For example:
	// `leading uu prokudin`. The total structure function is built starting
	// from the end backwards. So in this case, a ProkudinSfSet wrapped in two
	// MaskSfSet.
	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	// Load words into an array.
	std::vector<std::string> parts;
	std::stringstream ss(*params.sf_set);
	std::string next_part;
	while (ss >> next_part) {
		parts.push_back(next_part);
	}
	if (parts.empty()) {
		std::cerr << "No structure function provided." << std::endl;
		return ERROR_STRUCTURE_FUNCTIONS_PARSE;
	}

	// The last part is the base structure function.
	std::string base = parts.back();
	parts.pop_back();
	if (base == "prokudin") {
		std::cout << "Using Prokudin structure functions." << std::endl;
		tmd_out->reset();
		sf_out->reset(new sf::set::ProkudinSfSet());
	} else if (base == "test") {
		std::cout << "Using test structure functions." << std::endl;
		tmd_out->reset();
		sf_out->reset(new sf::set::TestSfSet(*params.target));
	} else {
		// TODO: Make this work on Windows as well, and make providing the
		// extension optional.
		std::string file_name = base + ".so";
		if (gSystem->Load(file_name.c_str()) != 0) {
			std::cerr
				<< "Failed to load structure function from shared library file "
				<< "'" << file_name << "'." << std::endl;
			return ERROR_FILE_NOT_FOUND;
		}
		TClass* sf_class = TClass::GetClass(base.c_str());
		if (sf_class->InheritsFrom("sidis::sf::SfSet")) {
			std::cout
				<< "Using structure functions from '"
				<< file_name << "'." << std::endl;
			tmd_out->reset();
			sf_out->reset(static_cast<sf::SfSet*>(sf_class->New()));
		} else if (sf_class->InheritsFrom("sidis::sf::TmdSet")) {
			std::cout
				<< "Using TMDs and FFs from '"
				<< file_name << "'." << std::endl;
			sf::TmdSet* tmd_ptr = static_cast<sf::TmdSet*>(sf_class->New());
			tmd_out->reset(tmd_ptr);
			sf_out->reset(new sf::TmdSfSet(*tmd_ptr));
		} else if (sf_class->InheritsFrom("sidis::sf::GaussianTmdSet")) {
			std::cout
				<< "Using Gaussian TMDs and FFs from '"
				<< file_name << "'." << std::endl;
			sf::GaussianTmdSet* tmd_ptr = static_cast<sf::GaussianTmdSet*>(sf_class->New());
			tmd_out->reset(tmd_ptr);
			sf_out->reset(new sf::GaussianTmdSfSet(*tmd_ptr));
		} else if (sf_class->InheritsFrom("sidis::sf::WwTmdSet")) {
			std::cout
				<< "Using WW-type TMDs and FFs from '"
				<< file_name << "'." << std::endl;
			sf::WwTmdSet* tmd_ptr = static_cast<sf::WwTmdSet*>(sf_class->New());
			tmd_out->reset(tmd_ptr);
			sf_out->reset(new sf::WwTmdSfSet(*tmd_ptr));
		} else if (sf_class->InheritsFrom("sidis::sf::GaussianWwTmdSet")) {
			std::cout
				<< "Using Gaussian WW-type TMDs and FFs from '"
				<< *params.sf_set << "'." << std::endl;
			sf::GaussianWwTmdSet* tmd_ptr = static_cast<sf::GaussianWwTmdSet*>(sf_class->New());
			tmd_out->reset(tmd_ptr);
			sf_out->reset(new sf::GaussianWwTmdSfSet(*tmd_ptr));
		} else {
			std::cerr
				<< "Couldn't find structure functions in file '"
				<< *params.sf_set << ".so'" << std::endl;
			return ERROR_STRUCTURE_FUNCTIONS_NOT_FOUND;
		}
	}

	// Apply filters to the base structure function based on every other part,
	// in reverse order.
	std::regex mask_regex("select([0-9]+)");
	std::smatch match;
	bool mask[sf::set::NUM_SF];
	std::fill_n(mask, sf::set::NUM_SF, true);
	for (auto it = parts.rbegin(); it != parts.rend(); ++it) {
		std::string part = *it;
		if (part == "leading") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_LEADING);
		} else if (part == "subleading") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_SUBLEADING);
		} else if (part == "uu") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_UU);
		} else if (part == "ul") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_UL);
		} else if (part == "ut") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_UT);
		} else if (part == "up") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_UP);
		} else if (part == "ux") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_UX);
		} else if (part == "lu") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_LU);
		} else if (part == "ll") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_LL);
		} else if (part == "lt") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_LT);
		} else if (part == "lp") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_LP);
		} else if (part == "lx") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_LX);
		} else if (part == "xu") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_XU);
		} else if (part == "xl") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_XL);
		} else if (part == "xt") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_XT);
		} else if (part == "xp") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_XP);
		} else if (part == "xx") {
			zip_and(mask, mask + sf::set::NUM_SF, sf::set::MASK_XX);
		} else if (std::regex_match(part, match, mask_regex)) {
			std::size_t idx = std::stoul(match[1]);
			if (idx > sf::set::NUM_SF) {
				std::cerr
					<< "Cannot filter on structure function index "
					<< idx << ": out of bounds." << std::endl;
				return ERROR_STRUCTURE_FUNCTIONS_PARSE;
			}
			bool select_mask[sf::set::NUM_SF] = { false };
			select_mask[idx] = true;
			zip_and(mask, mask + sf::set::NUM_SF, select_mask);
		} else {
			std::cerr
				<< "Unrecognized structure function filter '"
				<< part << "'." << std::endl;
			return ERROR_STRUCTURE_FUNCTIONS_PARSE;
		}
	}
	bool mask_full = true;
	for (std::size_t idx = 0; idx < sf::set::NUM_SF; ++idx) {
		if (!mask[idx]) {
			mask_full = false;
			break;
		}
	}
	if (!mask_full) {
		sf_out->reset(new sf::set::MaskSfSet(mask, std::move(*sf_out)));
	}
	return SUCCESS;
}

struct XsNRad : public TFoamIntegrand {
	Params params;
	cut::Cut cut;
	part::Particles ps;
	Real S;
	sf::SfSet const& sf;

	explicit XsNRad(Params params, cut::Cut cut, sf::SfSet const& sf) :
		params(params),
		cut(cut),
		ps(*params.target, *params.beam, *params.hadron, *params.mass_threshold),
		S(2. * mass(*params.target) * *params.beam_energy),
		sf(sf) { }

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 6) {
			return 0.;
		}
		kin::Kinematics kin;
		Real jacobian;
		if (!cut::take(cut, ps, S, vec, &kin, &jacobian)) {
			return 0.;
		}
		math::Vec3 eta = frame::hadron_from_target(kin) * *params.target_pol;
		// TODO: Evaluate when it is a good approximation to say that
		// `nrad ~ nrad_ir`. This happens because for small `k_0_bar`, the
		// contribution of `rad_f` integrated up to `k_0_bar` becomes vanishingly
		// small, so it can be neglected. However, this must be balanced with
		// choosing `k_0_bar` to be non-zero to avoid the infrared divergence in
		// the radiative part of the cross-section.
		Real xs;
		switch (*params.rc_method) {
		case RcMethod::NONE:
			xs = xs::born(kin, sf, *params.beam_pol, eta);
			break;
		case RcMethod::APPROX:
			xs = xs::nrad_ir(kin, sf, *params.beam_pol, eta, *params.k_0_bar);
			break;
		case RcMethod::EXACT:
			xs = xs::nrad(kin, sf, *params.beam_pol, eta, *params.k_0_bar);
			break;
		default:
			xs = 0.;
		}
		// Some kinematic regions will be out of range for the structure
		// functions, so return 0 in those cases.
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

struct XsRad : public TFoamIntegrand {
	Params params;
	cut::Cut cut;
	cut::CutRad cut_rad;
	part::Particles ps;
	Real S;
	sf::SfSet const& sf;

	explicit XsRad(
		Params params,
		cut::Cut cut,
		cut::CutRad cut_rad,
		sf::SfSet const& sf) :
		params(params),
		cut(cut),
		cut_rad(cut_rad),
		ps(*params.target, *params.beam, *params.hadron, *params.mass_threshold),
		S(2. * mass(*params.target) * *params.beam_energy),
		sf(sf) {
	}

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 9) {
			return 0.;
		}
		kin::KinematicsRad kin_rad;
		Real jacobian;
		if (!cut::take(cut, cut_rad, ps, S, vec, &kin_rad, &jacobian)) {
			return 0.;
		}
		kin::Kinematics kin = kin_rad.project();
		math::Vec3 eta = frame::hadron_from_target(kin) * *params.target_pol;
		Real xs = xs::rad(kin_rad, sf, *params.beam_pol, eta);
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

int command_help() {
	std::cout
		<< "Usage:"                                          << std::endl
		<< "  Prepare FOAM for Monte-Carlo generation"       << std::endl
		<< "    sidisgen --initialize <parameter file>"      << std::endl
		<< "  Generate events"                               << std::endl
		<< "    sidisgen --generate <parameter file>"        << std::endl
		<< "  List parameters used to produce file"          << std::endl
		<< "    sidisgen --inspect <output file>"            << std::endl
		<< "  Show parameter file format information"        << std::endl
		<< "    sidisgen --help-params"                      << std::endl;
	return SUCCESS;
}

int command_help_params() {
	std::cout
		<< "Parameter file format summary."                  << std::endl
		<< "For more detailed information, see docs."        << std::endl
		<< std::endl
		<< "event-file     <ROOT file>"                      << std::endl
		<< "rc-method      <none, approx, exact>"            << std::endl
		<< "gen-nrad       <true, false>"                    << std::endl
		<< "gen-rad        <true, false>"                    << std::endl
		<< "write-photon   <true, false>"                    << std::endl
		<< "foam-nrad-file <ROOT file>"                      << std::endl
		<< "foam-rad-file  <ROOT file>"                      << std::endl
		<< "sf-set         <prokudin, test, ROOT dict.>"     << std::endl
		<< "num-events     <integer>"                        << std::endl
		<< "num-init       <integer>"                        << std::endl
		<< "seed           <integer>"                        << std::endl
		<< "seed-init      <integer>"                        << std::endl
		<< "beam-energy    <energy (GeV)>"                   << std::endl
		<< "beam           <pid>"                            << std::endl
		<< "target         <pid>"                            << std::endl
		<< "mass-threshold <mass (GeV)>"                     << std::endl
		<< "hadron         <pid>"                            << std::endl
		<< "beam-pol       <real in [0, 1]>"                 << std::endl
		<< "target-pol     <vector in unit sphere>"          << std::endl
		<< "soft-threshold <energy (GeV)>"                   << std::endl
		<< "k-0-bar-cut    <min> <max>"                      << std::endl
		<< "x-cut          <min> <max>"                      << std::endl
		<< "y-cut          <min> <max>"                      << std::endl
		<< "z-cut          <min> <max>"                      << std::endl
		<< "ph-t-sq-cut    <min> <max>"                      << std::endl
		<< "phi-h-cut      <min> <max>"                      << std::endl
		<< "phi-cut        <min> <max>"                      << std::endl
		<< "tau-cut        <min> <max>"                      << std::endl
		<< "phi-k-cut      <min> <max>"                      << std::endl
		<< "Q-sq-cut       <min> <max>"                      << std::endl
		<< "t-cut          <min> <max>"                      << std::endl
		<< "w-cut          <min> <max>"                      << std::endl
		<< "r-cut          <min> <max>"                      << std::endl
		<< "mx-sq-cut      <min> <max>"                      << std::endl
		<< "q-0-cut        <min> <max>"                      << std::endl
		<< "k2-0-cut       <min> <max>"                      << std::endl
		<< "ph-0-cut       <min> <max>"                      << std::endl
		<< "k-0-cut        <min> <max>"                      << std::endl
		<< "theta-q-cut    <min> <max>"                      << std::endl
		<< "theta-k2-cut   <min> <max>"                      << std::endl
		<< "theta-ph-cut   <min> <max>"                      << std::endl
		<< "theta-k-cut    <min> <max>"                      << std::endl;
	return SUCCESS;
}

int command_version() {
	std::cout << "sidisgen "
		<< SIDIS_VERSION_MAJOR << "."
		<< SIDIS_VERSION_MINOR << "."
		<< SIDIS_VERSION_PATCH << "."
		<< SIDIS_VERSION_TWEAK << std::endl;
	return SUCCESS;
}

int command_inspect(std::string output_file_name) {
	Params params;
	TFile file(output_file_name.c_str());
	if (file.IsZombie()) {
		std::cerr
			<< "Output file '" << output_file_name
			<< "' not found." << std::endl;
		return ERROR_FILE_NOT_FOUND;
	}
	params.read_root(file);
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout
		<< std::scientific
		<< std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	params.write_stream(std::cout);
	std::cout.flags(flags);
	return SUCCESS;
}

int command_initialize(std::string params_file_name) {
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		std::cerr
			<< "Parameter file '" << params_file_name
			<< "' not found." << std::endl;
		return ERROR_FILE_NOT_FOUND;
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params;
	try {
		params.read_stream(params_file);
	} catch (std::exception const& e) {
		std::cerr
			<< "Failed to parse parameter file '" << params_file_name << "': "
			<< e.what() << std::endl;
		return ERROR_PARAMS_PARSE;
	}
	std::cout << std::endl;
	params.write_stream(std::cout);
	std::cout << std::endl;
	try {
		params.make_valid();
	} catch (std::exception const& e) {
		std::cerr
			<< "Invalid options in parameter file '" << params_file_name << "': "
			<< e.what() << std::endl;
		return ERROR_PARAMS_INVALID;
	}

	ULong_t N_init_nrad = *params.num_init >= 1 ? *params.num_init : 1;
	ULong_t N_init_rad  = *params.num_init >= 1 ? *params.num_init : 1;
	// Note that `TRandom3` uses the time as the seed if zero is provided.
	UInt_t seed = *params.seed_init >= 0 ? *params.seed_init : 0;
	TRandom3 random(seed);

	cut::Cut cut;
	cut::CutRad cut_rad;
	cuts(params, &cut, &cut_rad);

	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	int err = alloc_sf(params, &sf, &tmd);
	if (err != SUCCESS) {
		return err;
	}

	if (*params.gen_nrad) {
		std::cout << "Creating non-radiative FOAM file." << std::endl;
		TFile foam_nrad_file(params.foam_nrad_file->c_str(), "RECREATE");
		if (foam_nrad_file.IsZombie()) {
			std::cerr
				<< "Couldn't open or create file '" << *params.foam_nrad_file
				<< "'." << std::endl;
			return ERROR_FILE_NOT_CREATED;
		}
		foam_nrad_file.cd();
		params.write_root(foam_nrad_file);

		std::cout << "Non-radiative FOAM initialization." << std::endl;
		TFoam foam_nrad("FoamNRad");
		XsNRad xs_nrad(params, cut, *sf);
		foam_nrad.SetChat(0);
		foam_nrad.SetkDim(6);
		foam_nrad.SetRho(&xs_nrad);
		foam_nrad.SetPseRan(&random);
		foam_nrad.SetnSampl(N_init_nrad);
		foam_nrad.Initialize();
		foam_nrad.Write("FoamNRad");
	}

	if (*params.gen_rad && *params.rc_method != RcMethod::NONE) {
		std::cout << "Creating radiative FOAM file." << std::endl;
		TFile foam_rad_file(params.foam_rad_file->c_str(), "RECREATE");
		if (foam_rad_file.IsZombie()) {
			std::cerr
				<< "Couldn't open or create file '" << *params.foam_rad_file
				<< "'." << std::endl;
			return ERROR_FILE_NOT_CREATED;
		}
		foam_rad_file.cd();
		params.write_root(foam_rad_file);

		std::cout << "Radiative FOAM initialization." << std::endl;
		TFoam foam_rad("FoamRad");
		XsRad xs_rad(params, cut, cut_rad, *sf);
		foam_rad.SetChat(0);
		foam_rad.SetkDim(9);
		foam_rad.SetRho(&xs_rad);
		foam_rad.SetPseRan(&random);
		foam_rad.SetnSampl(N_init_rad);
		foam_rad.Initialize();
		foam_rad.Write("FoamRad");
	}

	std::cout << "Finished!" << std::endl;
	return SUCCESS;
}

int command_generate(std::string params_file_name) {
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		std::cerr
			<< "Parameter file '" << params_file_name
			<< "' not found." << std::endl;
		return ERROR_FILE_NOT_FOUND;
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params;
	try {
		params.read_stream(params_file);
	} catch (std::exception const& e) {
		std::cerr
			<< "Failed to parse parameter file '" << params_file_name << "': "
			<< e.what() << std::endl;
		return ERROR_PARAMS_PARSE;
	}
	std::cout << std::endl;
	params.write_stream(std::cout);
	std::cout << std::endl;
	try {
		params.make_valid();
	} catch (std::exception const& e) {
		std::cerr
			<< "Invalid options in parameter file '" << params_file_name << "': "
			<< e.what() << std::endl;
		return ERROR_PARAMS_INVALID;
	}

	// Fill out cut information.
	cut::Cut cut;
	cut::CutRad cut_rad;
	cuts(params, &cut, &cut_rad);

	// Load the structure functions.
	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	int err = alloc_sf(params, &sf, &tmd);
	if (err != SUCCESS) {
		return err;
	}

	std::cout
		<< "Opening event output file '" << *params.event_file
		<< "'." << std::endl;
	TFile event_file(params.event_file->c_str(), "RECREATE");
	if (event_file.IsZombie()) {
		std::cerr
			<< "Couldn't create file '" << *params.event_file
			<< "'." << std::endl;
		return ERROR_FILE_NOT_CREATED;
	}

	UInt_t seed = *params.seed >= 0 ? *params.seed : 0;
	TRandom3 random(seed);

	// Fill out the information for each type of event.
	std::vector<EventStats> event_stats;
	if (*params.gen_nrad) {
		std::cout
			<< "Reading non-radiative FOAM from file '" << *params.foam_nrad_file
			<< "'." << std::endl;
		std::unique_ptr<TFile> foam_nrad_file(
			new TFile(params.foam_nrad_file->c_str(), "OPEN"));
		if (foam_nrad_file->IsZombie()) {
			std::cerr
				<< "Couldn't open file '" << *params.foam_nrad_file
				<< "'." << std::endl;
			return ERROR_FILE_NOT_CREATED;
		}
		Params foam_nrad_params;
		foam_nrad_params.read_root(*foam_nrad_file);
		try {
			if (!foam_nrad_params.valid()) {
				throw std::runtime_error("Invalid FOAM parameters.");
			}
			params.compatible_with_foam(foam_nrad_params);
		} catch (std::exception const& e) {
			std::cerr
				<< "Couldn't use non-radiative FOAM from '" << *params.foam_nrad_file
				<< "' because it uses parameters incompatible with the "
				<< "provided parameter file '" << params_file_name
				<< "': " << e.what() << std::endl;
			return ERROR_FOAM_INCOMPATIBLE;
		}
		TFoam* foam_nrad = foam_nrad_file->Get<TFoam>("FoamNRad");
		if (foam_nrad == nullptr) {
			std::cerr
				<< "Failed to load non-radiative FOAM from file '"
				<< *params.foam_nrad_file << "'.";
			return ERROR_FOAM_NOT_FOUND;
		}
		std::unique_ptr<TFoamIntegrand> rho(new XsNRad(params, cut, *sf));
		foam_nrad->SetPseRan(&random);
		foam_nrad->ResetRho(rho.get());
		event_stats.push_back({
			EventType::NRAD,
			std::move(foam_nrad_file),
			std::move(rho),
			foam_nrad,
			0.,
			0.,
			0,
		});
	}

	if (*params.gen_rad) {
		std::cout
			<< "Reading radiative FOAM from file '" << *params.foam_rad_file
			<< "'." << std::endl;
		std::unique_ptr<TFile> foam_rad_file(
			new TFile(params.foam_rad_file->c_str(), "OPEN"));
		if (foam_rad_file->IsZombie()) {
			std::cerr
				<< "Couldn't open file '" << *params.foam_rad_file
				<< "'." << std::endl;
			return ERROR_FILE_NOT_CREATED;
		}
		Params foam_rad_params;
		foam_rad_params.read_root(*foam_rad_file);
		try {
			if (!foam_rad_params.valid()) {
				throw std::runtime_error("Invalid FOAM parameters.");
			}
			params.compatible_with_foam(foam_rad_params);
		} catch (std::exception const& e) {
			std::cerr
				<< "Couldn't use radiative FOAM from '" << *params.foam_rad_file
				<< "' because it uses parameters incompatible with the "
				<< "provided parameter file '" << params_file_name
				<< "': " << e.what() << std::endl;
			return ERROR_FOAM_INCOMPATIBLE;
		}
		TFoam* foam_rad = foam_rad_file->Get<TFoam>("FoamRad");
		if (foam_rad == nullptr) {
			std::cerr
				<< "Failed to load radiative FOAM from file '"
				<< *params.foam_rad_file << "'.";
			return ERROR_FOAM_NOT_FOUND;
		}
		std::unique_ptr<TFoamIntegrand> rho(new XsRad(params, cut, cut_rad, *sf));
		foam_rad->SetPseRan(&random);
		foam_rad->ResetRho(rho.get());
		event_stats.push_back({
			EventType::RAD,
			std::move(foam_rad_file),
			std::move(rho),
			foam_rad,
			0.,
			0.,
			0,
		});
	}

	part::Nucleus target = *params.target;
	part::Lepton beam = *params.beam;
	part::Hadron hadron = *params.hadron;
	math::Vec3 target_pol = *params.target_pol;
	part::Particles ps(target, beam, hadron, *params.mass_threshold);
	Real S = 2.*(*params.beam_energy)*ps.M;

	kin::Initial init(ps, *params.beam_energy);
	ULong_t N_gen = *params.num_events >= 0 ? *params.num_events : 0;

	event_file.cd();
	TTree events("Events", "Events");
	Int_t type;
	Double_t weight;
	Double_t jacobian;
	Double_t x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R;
	TLorentzVector p, k1, q, k2, ph, k;
	events.Branch("type", &type);
	events.Branch("weight", &weight);
	events.Branch("jacobian", &jacobian);
	events.Branch("x", &x);
	events.Branch("y", &y);
	events.Branch("z", &z);
	events.Branch("ph_t_sq", &ph_t_sq);
	events.Branch("phi_h", &phi_h);
	events.Branch("phi", &phi);
	events.Branch("tau", &tau);
	events.Branch("phi_k", &phi_k);
	events.Branch("R", &R);
	events.Branch("p", "TLorentzVector", &p);
	events.Branch("k1", "TLorentzVector", &k1);
	events.Branch("q", "TLorentzVector", &q);
	events.Branch("k2", "TLorentzVector", &k2);
	events.Branch("ph", "TLorentzVector", &ph);
	if (*params.gen_rad && *params.write_photon) {
		events.Branch("k", "TLorentzVector", &k);
	}
	params.write_root(event_file);

	std::cout << "Generating events." << std::endl;
	bool update_progress = true;
	ULong_t percent = 0;
	ULong_t next_percent_rem = N_gen % 100;
	ULong_t next_percent = N_gen / 100;
	for (ULong_t event_idx = 0; event_idx < N_gen; ++event_idx) {
		while (event_idx >= next_percent + (next_percent_rem != 0)) {
			percent += 1;
			next_percent = math::prod_div(N_gen, percent + 1, 100, next_percent_rem);
			update_progress = true;
		}
		if (update_progress) {
			write_progress_bar(std::cout, percent);
			std::cout << '\r';
			std::cout.flush();
			update_progress = false;
		}
		// Estimate the total radiative and non-radiative cross-sections and
		// generate a radiative/non-radiative event accordingly. The total
		// cross-section estimates are improved as more events are generated.
		Double_t total_xs = 0.;
		Double_t total_xs_err = 0.;
		for (EventStats& stats : event_stats) {
			stats.foam->GetIntegMC(stats.xs, stats.xs_err);
			if (!std::isfinite(stats.xs) || stats.xs == 0.) {
				stats.xs = 0.;
				stats.xs_err = std::numeric_limits<Double_t>::max();
			}
			total_xs += stats.xs;
			total_xs_err += stats.xs_err;
		}
		std::size_t choose_event_type;
		if (event_idx == 0) {
			// On the first event, we don't know anything about the total cross-
			// sections, so choose the event type arbitrarily.
			choose_event_type = 0;
		} else {
			// We choose the type of event to generate as that for which the
			// # of events generated / total # of events is furthest from the
			// cross-section ratio of the two event types.
			Double_t ratio_min = std::numeric_limits<Double_t>::infinity();
			for (std::size_t idx = 0; idx < event_stats.size(); ++idx) {
				EventStats const& stats = event_stats[idx];
				Double_t target = (stats.xs + stats.xs_err) / (total_xs + total_xs_err);
				Double_t ratio = static_cast<Double_t>(stats.num_events) / event_idx
					/ target;
				if (ratio <= ratio_min) {
					ratio_min = ratio;
					choose_event_type = idx;
				}
			}
		}

		// The event vector can store up to the number of dimensions of any of
		// the FOAMs.
		Double_t event_vec[9];
		weight = event_stats[choose_event_type].foam->MCgenerate(event_vec);
		type = static_cast<Int_t>(event_stats[choose_event_type].type);
		event_stats[choose_event_type].num_events += 1;

		// Fill in the branches.
		kin::Kinematics kin;
		kin::KinematicsRad kin_rad;
		switch (event_stats[choose_event_type].type) {
		case EventType::NRAD:
			// Non-radiative event.
			if (cut::take(cut, ps, S, event_vec, &kin, &jacobian)) {
				kin::Final fin(init, target_pol, kin);
				// Fill in branches.
				x = kin.x;
				y = kin.y;
				z = kin.z;
				ph_t_sq = kin.ph_t_sq;
				phi_h = kin.phi_h;
				phi = kin.phi;
				p = convert_vec4(init.p);
				k1 = convert_vec4(init.k1);
				q = convert_vec4(fin.q);
				k2 = convert_vec4(fin.k2);
				ph = convert_vec4(fin.ph);
			} else {
				// Make sure invalid data isn't written to the events.
				weight = 0.;
			}
			break;
		case EventType::RAD:
			// Radiative event.
			if (cut::take(cut, cut_rad, ps, S, event_vec, &kin_rad, &jacobian)) {
				kin::FinalRad fin(init, target_pol, kin_rad);
				// Fill in branches.
				x = kin_rad.x;
				y = kin_rad.y;
				z = kin_rad.z;
				ph_t_sq = kin_rad.ph_t_sq;
				phi_h = kin_rad.phi_h;
				phi = kin_rad.phi;
				tau = kin_rad.tau;
				phi_k = kin_rad.phi_k;
				R = kin_rad.R;
				p = convert_vec4(init.p);
				k1 = convert_vec4(init.k1);
				q = convert_vec4(fin.q);
				k2 = convert_vec4(fin.k2);
				ph = convert_vec4(fin.ph);
				k = convert_vec4(fin.k);
			} else {
				weight = 0.;
			}
			break;
		default:
			throw std::runtime_error("Unrecognized event type.");
		}
		events.Fill();
	}
	write_progress_bar(std::cout, 100);
	std::cout << std::endl;
	std::cout << "Writing events to file." << std::endl;
	event_file.cd();
	events.Write();

	std::cout << "Statistics:" << std::endl;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::scientific << std::setprecision(6);
	TArrayD total_xs(2);
	TArrayD total_xs_err(2);
	for (EventStats& stats : event_stats) {
		Int_t idx = static_cast<Int_t>(stats.type);
		stats.foam->GetIntegMC(stats.xs, stats.xs_err);
		total_xs.SetAt(idx, stats.xs);
		total_xs_err.SetAt(idx, stats.xs_err);
		switch (stats.type) {
		case EventType::NRAD:
			std::cout << "\tNon-radiative events:" << std::endl;
			break;
		case EventType::RAD:
			std::cout << "\tRadiative events:" << std::endl;
			break;
		default:
			throw std::runtime_error("Unrecognized event type.");
		}
		std::cout << "\t\tCount:         " << stats.num_events << std::endl;
		std::cout << "\t\tCross-section: " << stats.xs << " ± " << stats.xs_err << std::endl;
	}
	std::cout.flags(flags);

	// Write total cross-sections to file.
	event_file.cd();
	event_file.WriteObject(&total_xs, "xs_total");
	event_file.WriteObject(&total_xs_err, "xs_total_err");

	return SUCCESS;
}

}

int main(int argc, char** argv) {
#ifdef NDEBUG
	gErrorIgnoreLevel = kFatal;
	gErrorAbortLevel = kFatal;
#else
	gErrorIgnoreLevel = kWarning;
	gErrorAbortLevel = kError;
#endif
	// Parse arguments.
	if (argc <= 1) {
		std::cout << "Try `sidisgen --help`." << std::endl;
		return SUCCESS;
	}
	std::string command = argv[1];
	if (command == "--help" || command == "-?") {
		return command_help();
	} else if (command == "--help-params") {
		return command_help_params();
	} else if (command == "--version" || command == "-v") {
		return command_version();
	} else if (command == "--inspect" || command == "-i") {
		if (argc > 3) {
			std::cerr
				<< "Unexpected argument " << argv[3] << "." << std::endl;
			return ERROR_ARG_PARSE;
		} else if (argc < 3) {
			std::cerr
				<< "Expected ROOT file argument for inspection." << std::endl;
			return ERROR_ARG_PARSE;
		}
		return command_inspect(argv[2]);
	} else if (command == "--initialize" || command == "-i") {
		if (argc > 3) {
			std::cerr
				<< "Unexpected argument " << argv[3] << "." << std::endl;
			return ERROR_ARG_PARSE;
		} else if (argc < 3) {
			std::cerr
				<< "Expected parameter file argument." << std::endl;
			return ERROR_ARG_PARSE;
		}
		return command_initialize(argv[2]);
	} else if (command == "--generate" || command == "-g") {
		if (argc > 3) {
			std::cerr
				<< "Unexpected argument " << argv[3] << "." << std::endl;
			return ERROR_ARG_PARSE;
		} else if (argc < 3) {
			std::cerr
				<< "Expected parameter file argument." << std::endl;
			return ERROR_ARG_PARSE;
		}
		return command_generate(argv[2]);
	} else {
		std::cerr << "Unrecognized command " << command << "." << std::endl;
		return ERROR_ARG_PARSE;
	}
}

