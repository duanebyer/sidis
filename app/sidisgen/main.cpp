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
#include <TArrayL.h>
#include <TBranch.h>
#include <TChain.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>

#include <sidis/sidis.hpp>
#include <sidis/extra/math.hpp>
#include <sidis/sf_set/mask.hpp>
#include <sidis/sf_set/prokudin.hpp>
#include <sidis/sf_set/test.hpp>

#include "event_generator.hpp"
#include "event_stats.hpp"
#include "exception.hpp"
#include "params.hpp"
#include "utility.hpp"

using namespace sidis;

namespace {

// Converts between the `sidis` 4-vector type and the ROOT 4-vector type.
TLorentzVector convert_vec4(math::Vec4 vec) {
	return TLorentzVector(vec.x, vec.y, vec.z, vec.t);
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
void alloc_sf(
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
		throw Exception(
			ERROR_STRUCTURE_FUNCTIONS_PARSE,
			"No structure function provided.");
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
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				std::string("Failed to load structure function from shared ")
				+ "library file '" + file_name + "'.");
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
			throw Exception(
				ERROR_STRUCTURE_FUNCTIONS_NOT_FOUND,
				std::string("Couldn't find structure functions in file '")
				+ *params.sf_set + ".so'.");
		}
	}

	// Apply filters to the base structure function based on every other part,
	// in reverse order.
	std::regex mask_regex("select([0-9]+)", std::regex_constants::basic);
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
				throw Exception(
					ERROR_STRUCTURE_FUNCTIONS_PARSE,
					std::string("Cannot filter on structure function index ")
					+ std::to_string(idx) + " because out of bounds.");
			}
			bool select_mask[sf::set::NUM_SF] = { false };
			select_mask[idx] = true;
			zip_and(mask, mask + sf::set::NUM_SF, select_mask);
		} else {
			throw Exception(
				ERROR_STRUCTURE_FUNCTIONS_PARSE,
				std::string("Unrecognized structure function filter '") + part
				+ "'.");
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
}

int command_help() {
	std::cout
		<< "Usage:"                                                << std::endl
		<< "  Prepare FOAM for Monte-Carlo generation"             << std::endl
		<< "    sidisgen --initialize <parameter file>"            << std::endl
		<< "  Generate events"                                     << std::endl
		<< "    sidisgen --generate <parameter file>"              << std::endl
		<< "  List parameters used to produce file"                << std::endl
		<< "    sidisgen --inspect <output file>"                  << std::endl
		<< "  Merge multiple event files into one"                 << std::endl
		<< "    sidisgen --merge-soft <output file> <input files>" << std::endl
		<< "    sidisgen --merge-hard <output file> <input files>" << std::endl
		<< "  Show parameter file format information"              << std::endl
		<< "    sidisgen --help-params"                            << std::endl;
	return SUCCESS;
}

int command_help_params() {
	std::cout
		<< "Parameter file format summary."                  << std::endl
		<< "For more detailed information, see docs."        << std::endl
		<< std::endl
		<< "event_file     <ROOT file>"                      << std::endl
		<< "rc_method      <none, approx, exact>"            << std::endl
		<< "gen_nrad       <on, off>"                        << std::endl
		<< "gen_rad        <on, off>"                        << std::endl
		<< "write_momenta  <on, off>"                        << std::endl
		<< "write_photon   <on, off>"                        << std::endl
		<< "foam_nrad_file <ROOT file>"                      << std::endl
		<< "foam_rad_file  <ROOT file>"                      << std::endl
		<< "sf_set         <prokudin, test, ROOT dict.>"     << std::endl
		<< "num_events     <integer>"                        << std::endl
		<< "num_init       <integer>"                        << std::endl
		<< "rej_weight     <real in [1, ∞)>"                 << std::endl
		<< "seed           <integer>"                        << std::endl
		<< "seed_init      <integer>"                        << std::endl
		<< "beam_energy    <energy (GeV)>"                   << std::endl
		<< "beam           <pid>"                            << std::endl
		<< "target         <pid>"                            << std::endl
		<< "mass_threshold <mass (GeV)>"                     << std::endl
		<< "hadron         <pid>"                            << std::endl
		<< "beam_pol       <real in [0, 1]>"                 << std::endl
		<< "target_pol     <vector in unit sphere>"          << std::endl
		<< "soft_threshold <energy (GeV)>"                   << std::endl
		<< "k_0_bar_cut    <min> <max>"                      << std::endl
		<< "x_cut          <min> <max>"                      << std::endl
		<< "y_cut          <min> <max>"                      << std::endl
		<< "z_cut          <min> <max>"                      << std::endl
		<< "ph_t_sq_cut    <min> <max>"                      << std::endl
		<< "phi_h_cut      <min> <max>"                      << std::endl
		<< "phi_cut        <min> <max>"                      << std::endl
		<< "tau_cut        <min> <max>"                      << std::endl
		<< "phi_k_cut      <min> <max>"                      << std::endl
		<< "R_cut          <min> <max>"                      << std::endl
		<< "Q_sq_cut       <min> <max>"                      << std::endl
		<< "t_cut          <min> <max>"                      << std::endl
		<< "W_sq_cut       <min> <max>"                      << std::endl
		<< "r_cut          <min> <max>"                      << std::endl
		<< "mx_sq_cut      <min> <max>"                      << std::endl
		<< "q_0_cut        <min> <max>"                      << std::endl
		<< "k2_0_cut       <min> <max>"                      << std::endl
		<< "ph_0_cut       <min> <max>"                      << std::endl
		<< "k_0_cut        <min> <max>"                      << std::endl
		<< "theta_q_cut    <min> <max>"                      << std::endl
		<< "theta_k2_cut   <min> <max>"                      << std::endl
		<< "theta_ph_cut   <min> <max>"                      << std::endl
		<< "theta_k_cut    <min> <max>"                      << std::endl;
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

int command_inspect(char const* file_name) {
	Params params;
	TFile file(file_name, "OPEN");
	if (file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("File '") + file_name + "' not found.");
	}
	params.read_root(file);
	std::cout << "Parameters:" << std::endl;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout
		<< std::scientific
		<< std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	params.write_stream(std::cout);
	std::cout.flags(flags);
	std::cout << std::endl;

	TArrayD* xs = file.Get<TArrayD>("stats/xs");
	TArrayD* xs_err = file.Get<TArrayD>("stats/xs_err");
	TArrayD* norm = file.Get<TArrayD>("stats/norm");
	TArrayD* efficiency = file.Get<TArrayD>("stats/efficiency");
	TArrayL* num_events = file.Get<TArrayL>("stats/num_events");
	if (
			xs == nullptr
			|| xs_err == nullptr
			|| norm == nullptr
			|| efficiency == nullptr
			|| num_events == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Couldn't find statistics in file '") + file_name
			+ "'.");
	}
	std::cout << "Statistics:" << std::endl;
	flags = std::cout.flags();
	std::cout << std::scientific << std::setprecision(6);
	for (std::size_t type_idx = 0; type_idx < NUM_TYPES + 1; ++type_idx) {
		if (num_events->At(type_idx) == 0) {
			continue;
		}
		if (type_idx == NUM_TYPES) {
			std::cout << "\ttotal:" << std::endl;
		} else {
			EventType type = static_cast<EventType>(type_idx);
			std::cout << "\t" << event_type_name(type) << ":" << std::endl;
		}
		std::cout << "\t\tcount:         " << num_events->At(type_idx) << std::endl;
		std::cout << "\t\tcross-section: " << xs->At(type_idx) << " ± " << xs_err->At(type_idx) << std::endl;
		std::cout << "\t\tnorm:          " << norm->At(type_idx) << std::endl;
		std::cout << "\t\tefficiency:    " << efficiency->At(type_idx) << std::endl;
	}
	std::cout.flags(flags);

	return SUCCESS;
}

int command_initialize(char const* params_file_name) {
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Parameter file '") + params_file_name + "' not "
			+ "found.");
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params;
	try {
		params.read_stream(params_file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_PARSE,
			std::string("Failed to parse parameter file '")
			+ params_file_name + "': " + e.what());
	}
	std::cout << std::endl;
	params.write_stream(std::cout);
	std::cout << std::endl;
	try {
		params.make_valid();
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			std::string("Invalid parameter file '") + params_file_name + "': "
			+ e.what());
	}

	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	alloc_sf(params, &sf, &tmd);
	UInt_t seed = *params.seed_init >= 0 ? *params.seed_init : 0;

	if (*params.gen_nrad) {
		TRandom3 random(seed);
		std::cout
			<< "Initializing non-radiative FOAM in file '"
			<< *params.foam_nrad_file << "'." << std::endl;
		EventGenerator::write(EventType::NRAD, params, *sf, &random);
	}
	if (*params.gen_rad) {
		TRandom3 random(seed);
		std::cout
			<< "Initializing radiative FOAM in file '"
			<< *params.foam_rad_file << "'." << std::endl;
		EventGenerator::write(EventType::RAD, params, *sf, &random);
	}
	std::cout << "Finished!" << std::endl;
	return SUCCESS;
}

int command_generate(char const* params_file_name) {
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Parameter file '") + params_file_name + "' not "
			+ "found.");
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params;
	try {
		params.read_stream(params_file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_PARSE,
			std::string("Failed to parse parameter file '")
			+ params_file_name + "': " + e.what());
	}
	std::cout << std::endl;
	params.write_stream(std::cout);
	std::cout << std::endl;
	try {
		params.make_valid();
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			std::string("Invalid parameter file '") + params_file_name + "': "
			+ e.what());
	}

	// Load the structure functions.
	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	alloc_sf(params, &sf, &tmd);

	std::cout
		<< "Opening event output file '" << *params.event_file << "'."
		<< std::endl;
	TFile event_file(params.event_file->c_str(), "RECREATE");
	if (event_file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create file '") + *params.event_file + "'.");
	}

	UInt_t seed = *params.seed >= 0 ? *params.seed : 0;
	TRandom3 random(seed);

	// Load the event generators from file.
	std::vector<EventGenerator> event_generators;
	if (*params.gen_nrad) {
		std::cout
			<< "Loading non-radiative FOAM from file '"
			<< *params.foam_nrad_file << "'." << std::endl;
		event_generators.push_back(EventGenerator::read(
			EventType::NRAD,
			params,
			*sf,
			&random));
	}
	if (*params.gen_rad) {
		std::cout
			<< "Loading radiative FOAM from file '"
			<< *params.foam_rad_file << "'." << std::endl;
		event_generators.push_back(EventGenerator::read(
			EventType::RAD,
			params,
			*sf,
			&random));
	}

	part::Nucleus target = *params.target;
	part::Lepton beam = *params.beam;
	part::Hadron hadron = *params.hadron;
	math::Vec3 target_pol = *params.target_pol;
	part::Particles ps(target, beam, hadron, *params.Mth);
	Real S = 2.*(*params.beam_energy)*ps.M;

	kin::Initial init(ps, *params.beam_energy);
	ULong_t N_gen = *params.num_events >= 0 ? *params.num_events : 0;

	event_file.cd();
	TTree events("events", "events");
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
	if (*params.gen_rad && *params.write_momenta) {
		events.Branch("p", "TLorentzVector", &p);
		events.Branch("k1", "TLorentzVector", &k1);
		events.Branch("q", "TLorentzVector", &q);
		events.Branch("k2", "TLorentzVector", &k2);
		events.Branch("ph", "TLorentzVector", &ph);
		if (*params.gen_rad && *params.write_photon) {
			events.Branch("k", "TLorentzVector", &k);
		}
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
		// Choose a type of event (ex. radiative or non-radiative) to generate.
		std::size_t choose_event_type;
		if (event_idx == 0) {
			// On the first event, we don't know anything about the total cross-
			// sections, so choose the event type arbitrarily.
			choose_event_type = 0;
		} else {
			// We choose the type of event to generate as that for which the
			// # of events generated / total # of events is furthest from the
			// cross-section ratio of the two event types.
			Double_t ratio_max = 0.;
			for (std::size_t idx = 0; idx < event_generators.size(); ++idx) {
				EventStats stats = event_generators[idx].stats();
				Double_t ratio = stats.xs / stats.weight_total;
				if (!std::isfinite(ratio)) {
					ratio = std::numeric_limits<Double_t>::infinity();
				}
				if (ratio > ratio_max) {
					ratio_max = ratio;
					choose_event_type = idx;
				}
			}
		}

		// The event vector can store up to the number of dimensions of any of
		// the FOAMs.
		Double_t event_vec[9];
		weight = event_generators[choose_event_type].generate(event_vec);
		type = static_cast<Int_t>(event_generators[choose_event_type].type());

		// Fill in the branches.
		kin::Kinematics kin;
		kin::KinematicsRad kin_rad;
		switch (event_generators[choose_event_type].type()) {
		case EventType::NRAD:
			// Non-radiative event.
			{
				// Fill in branches.
				x = event_vec[0];
				y = event_vec[1];
				z = event_vec[2];
				ph_t_sq = event_vec[3];
				phi_h = event_vec[4];
				phi = event_vec[5];
				tau = 0.;
				phi_k = 0.;
				R = 0.;
				if (*params.write_momenta) {
					kin::PhaseSpace ph_space {
						x, y, z,
						ph_t_sq, phi_h, phi,
					};
					kin::Kinematics kin(ps, S, ph_space);
					kin::Final fin(init, target_pol, kin);
					p = convert_vec4(init.p);
					k1 = convert_vec4(init.k1);
					q = convert_vec4(fin.q);
					k2 = convert_vec4(fin.k2);
					ph = convert_vec4(fin.ph);
				}
			}
			break;
		case EventType::RAD:
			// Radiative event.
			{
				// Fill in branches.
				x = event_vec[0];
				y = event_vec[1];
				z = event_vec[2];
				ph_t_sq = event_vec[3];
				phi_h = event_vec[4];
				phi = event_vec[5];
				tau = event_vec[6];
				phi_k = event_vec[7];
				R = event_vec[8];
				if (*params.write_momenta) {
					kin::PhaseSpaceRad ph_space {
						x, y, z,
						ph_t_sq, phi_h, phi,
						tau, phi_k, R,
					};
					kin::KinematicsRad kin(ps, S, ph_space);
					kin::FinalRad fin(init, target_pol, kin);
					p = convert_vec4(init.p);
					k1 = convert_vec4(init.k1);
					q = convert_vec4(fin.q);
					k2 = convert_vec4(fin.k2);
					ph = convert_vec4(fin.ph);
					k = convert_vec4(fin.k);
				}
			}
			break;
		case EventType::EXCL:
			throw std::runtime_error("Exclusive events not supported.");
		default:
			throw std::runtime_error("Unrecognized event type.");
		}
		events.Fill();
	}
	write_progress_bar(std::cout, 100);
	std::cout << std::endl;
	std::cout << "Writing events to file." << std::endl;
	event_file.WriteObject(&events, events.GetName());

	// Produce statistics and write them to file.
	std::cout << "Statistics:" << std::endl;
	TDirectory* stats_dir = event_file.mkdir("stats", "stats");
	if (stats_dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create directory 'stats' in ROOT file."));
	}
	stats_dir->cd();
	TArrayD xs(NUM_TYPES + 1);
	TArrayD xs_err(NUM_TYPES + 1);
	TArrayD weight_total(NUM_TYPES + 1);
	TArrayD weight_sq_total(NUM_TYPES + 1);
	TArrayD norm(NUM_TYPES + 1);
	TArrayD efficiency(NUM_TYPES + 1);
	TArrayL num_events(NUM_TYPES + 1);
	xs_err.Reset(std::numeric_limits<Double_t>::infinity());
	std::vector<EventStats> stats;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::scientific << std::setprecision(6);
	for (EventGenerator& gen : event_generators) {
		Int_t idx = static_cast<Int_t>(gen.type());
		stats.push_back(gen.stats());
		xs.SetAt(gen.stats().xs, idx);
		xs_err.SetAt(gen.stats().xs_err, idx);
		weight_total.SetAt(gen.stats().weight_total, idx);
		weight_sq_total.SetAt(gen.stats().weight_sq_total, idx);
		norm.SetAt(gen.stats().norm(), idx);
		efficiency.SetAt(gen.stats().efficiency(), idx);
		num_events.SetAt(gen.stats().num_events, idx);
		std::cout << "\t" << event_type_name(gen.type()) << " events:" << std::endl;
		std::cout << "\t\tcount:         " << gen.stats().num_events << std::endl;
		std::cout << "\t\tcross-section: " << gen.stats().xs << " ± " << gen.stats().xs_err << std::endl;
		std::cout << "\t\tnorm:          " << gen.stats().norm() << std::endl;
		std::cout << "\t\tefficiency:    " << gen.stats().efficiency() << " ± " << gen.stats().efficiency_err() << std::endl;
	}
	EventStats total_stats = EventStats::total(stats.begin(), stats.end());
	xs.SetAt(total_stats.xs, NUM_TYPES);
	xs_err.SetAt(total_stats.xs_err, NUM_TYPES);
	weight_total.SetAt(total_stats.weight_total, NUM_TYPES);
	weight_sq_total.SetAt(total_stats.weight_sq_total, NUM_TYPES);
	norm.SetAt(total_stats.norm(), NUM_TYPES);
	efficiency.SetAt(total_stats.efficiency(), NUM_TYPES);
	num_events.SetAt(total_stats.num_events, NUM_TYPES);
	std::cout << "\ttotal:" << std::endl;
		std::cout << "\t\tcount:         " << total_stats.num_events << std::endl;
		std::cout << "\t\tcross-section: " << total_stats.xs << " ± " << total_stats.xs_err << std::endl;
		std::cout << "\t\tnorm:          " << total_stats.norm() << std::endl;
		std::cout << "\t\tefficiency:    " << total_stats.efficiency() << std::endl;
	std::cout.flags(flags);
	stats_dir->WriteObject(&xs, "xs");
	stats_dir->WriteObject(&xs_err, "xs_err");
	stats_dir->WriteObject(&weight_total, "weight_total");
	stats_dir->WriteObject(&weight_sq_total, "weight_sq_total");
	stats_dir->WriteObject(&norm, "norm");
	stats_dir->WriteObject(&efficiency, "efficiency");
	stats_dir->WriteObject(&num_events, "num_events");

	return SUCCESS;
}

int command_merge_soft(
		char const* file_out_name,
		std::vector<char const*> file_names) {
	TFile file_out(file_out_name, "CREATE");
	if (file_out.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create file '") + file_out_name + "'.");
	}

	std::cout << "Merging parameters from files." << std::endl;
	Params params_out;
	bool first = true;
	for (char const* file_name : file_names) {
		TFile file(file_name, "OPEN");
		Params params;
		params.read_root(file);
		if (!params.valid()) {
			throw Exception(
				ERROR_PARAMS_INVALID,
				std::string("Invalid parameters in '") + file_name + "'.");
		}
		if (first) {
			params_out = params;
			first = false;
		} else {
			params_out.merge(params);
		}
	}
	params_out.write_root(file_out);

	std::cout << "Merging statistics from files." << std::endl;
	std::vector<EventStats> stats[NUM_TYPES + 1];
	for (char const* file_name : file_names) {
		TFile file(file_name, "OPEN");
		TArrayD* xs = file.Get<TArrayD>("stats/xs");
		TArrayD* xs_err = file.Get<TArrayD>("stats/xs_err");
		TArrayD* weight_total = file.Get<TArrayD>("stats/weight_total");
		TArrayD* weight_sq_total = file.Get<TArrayD>("stats/weight_sq_total");
		TArrayL* num_events = file.Get<TArrayL>("stats/num_events");
		for (std::size_t type_idx = 0; type_idx < NUM_TYPES + 1; ++type_idx) {
			EventStats next_stats;
			next_stats.xs = xs->At(type_idx);
			next_stats.xs_err = xs_err->At(type_idx);
			next_stats.weight_total = weight_total->At(type_idx);
			next_stats.weight_sq_total = weight_sq_total->At(type_idx);
			next_stats.num_events = num_events->At(type_idx);
			stats[type_idx].push_back(next_stats);
		}
	}
	EventStats stats_out[NUM_TYPES + 1];
	TDirectory* stats_dir = file_out.mkdir("stats");
	if (stats_dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create directory 'stats' in ROOT file."));
	}
	stats_dir->cd();
	TArrayD xs_out(NUM_TYPES + 1);
	TArrayD xs_err_out(NUM_TYPES + 1);
	TArrayD weight_total_out(NUM_TYPES + 1);
	TArrayD weight_sq_total_out(NUM_TYPES + 1);
	TArrayD norm_out(NUM_TYPES + 1);
	TArrayD efficiency_out(NUM_TYPES + 1);
	TArrayL num_events_out(NUM_TYPES + 1);
	for (std::size_t type_idx = 0; type_idx < NUM_TYPES + 1; ++type_idx) {
		stats_out[type_idx] = EventStats::average(
			stats[type_idx].begin(),
			stats[type_idx].end());
		xs_out.SetAt(stats_out[type_idx].xs, type_idx);
		xs_err_out.SetAt(stats_out[type_idx].xs_err, type_idx);
		weight_total_out.SetAt(stats_out[type_idx].weight_total, type_idx);
		weight_sq_total_out.SetAt(stats_out[type_idx].weight_sq_total, type_idx);
		norm_out.SetAt(stats_out[type_idx].norm(), type_idx);
		efficiency_out.SetAt(stats_out[type_idx].efficiency(), type_idx);
		num_events_out.SetAt(stats_out[type_idx].num_events, type_idx);
	}
	stats_dir->WriteObject(&xs_out, "xs");
	stats_dir->WriteObject(&xs_err_out, "xs_err");
	stats_dir->WriteObject(&weight_total_out, "weight_total");
	stats_dir->WriteObject(&weight_sq_total_out, "weight_sq_total");
	stats_dir->WriteObject(&norm_out, "norm");
	stats_dir->WriteObject(&efficiency_out, "efficiency");
	stats_dir->WriteObject(&num_events_out, "num_events");

	std::cout << "Merging events from files." << std::endl;
	file_out.cd();
	TChain chain("events", "events");
	for (char const* file_name : file_names) {
		std::cout << "\t" << file_name << std::endl;
		chain.Add(file_name);
	}
	file_out.WriteObject(&chain, chain.GetName());

	return SUCCESS;
}

int command_merge_hard(
		char const* file_out_name,
		std::vector<char const*> file_names) {
	static_cast<void>(file_out_name);
	static_cast<void>(file_names);
	throw Exception(
		ERROR_UNIMPLEMENTED,
		"Hard merges are not yet supported.");
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
	try {
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
				throw Exception(
					ERROR_ARG_PARSE,
					std::string("Unexpected argument '") + argv[3] + "'.");
			} else if (argc < 3) {
				throw Exception(
					ERROR_ARG_PARSE,
					"Expected ROOT file argument for inspection.");
			}
			return command_inspect(argv[2]);
		} else if (command == "--initialize" || command == "-i") {
			if (argc > 3) {
				throw Exception(
					ERROR_ARG_PARSE,
					std::string("Unexpected argument '") + argv[3] + "'.");
			} else if (argc < 3) {
				throw Exception(
					ERROR_ARG_PARSE,
					"Expected parameter file argument.");
			}
			return command_initialize(argv[2]);
		} else if (command == "--generate" || command == "-g") {
			if (argc > 3) {
				throw Exception(
					ERROR_ARG_PARSE,
					std::string("Unexpected argument '") + argv[3] + "'.");
			} else if (argc < 3) {
				throw Exception(
					ERROR_ARG_PARSE,
					"Expected parameter file argument.");
			}
			return command_generate(argv[2]);
		} else if (command == "--merge-soft") {
			if (argc < 4) {
				throw Exception(
					ERROR_ARG_PARSE,
					"Expected target file and source file arguments.");
			}
			std::vector<char const*> files;
			for (int idx = 3; idx < argc; ++idx) {
				files.push_back(argv[idx]);
			}
			return command_merge_soft(argv[2], files);
		} else if (command == "--merge-hard") {
			if (argc < 4) {
				throw Exception(
					ERROR_ARG_PARSE,
					"Expected target file and source file arguments.");
			}
			std::vector<char const*> files;
			for (int idx = 3; idx < argc; ++idx) {
				files.push_back(argv[idx]);
			}
			return command_merge_hard(argv[2], files);
		} else {
			throw Exception(
				ERROR_ARG_PARSE,
				std::string("Unrecognized command '") + command + "'.");
		}
	} catch (Exception const& e) {
		std::cerr << "Fatal error: " << e.what() << std::endl;
		return e.error_code;
	} catch (std::exception const& e) {
		std::cerr << "Fatal error: " << e.what() << std::endl;
		return ERROR;
	} catch (...) {
		std::cerr << "Fatal error: An unknown error occurred." << std::endl;
		return ERROR;
	}
}

