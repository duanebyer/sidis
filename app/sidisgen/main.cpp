#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <regex>
#include <string>
#include <sstream>
#include <stdexcept>
#include <type_traits>
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
#include <TSystem.h>
#include <TTree.h>

#include <bubble.hpp>

#include <sidis/sidis.hpp>
#include <sidis/extra/math.hpp>
#include <sidis/sf_set/mask.hpp>
#include <sidis/sf_set/prokudin.hpp>
#include <sidis/sf_set/test.hpp>

#include "generator.hpp"
#include "event_type.hpp"
#include "exception.hpp"
#include "params.hpp"
#include "utility.hpp"

using namespace sidis;

static_assert(
	std::is_same<Real, Double_t>::value,
	"ROOT `Double_t` type must be same as C++ `double` type.");

namespace {

struct DrawProgressBar final :
		public bubble::ExploreProgressReporter<Real>,
		public bubble::TuneProgressReporter<Real> {
	void operator()(bubble::ExploreProgress<Real> progress) override {
		Real percent = 100. * progress.progress;
		write_progress_bar(std::cout, static_cast<unsigned>(percent));
		std::cout << '\r';
		std::cout << std::flush;
	}
	void operator()(bubble::TuneProgress<Real> progress) override {
		Real percent = 100. * progress.progress;
		write_progress_bar(std::cout, static_cast<unsigned>(percent));
		std::cout << '\r';
		std::cout << std::flush;
	}
};

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
		<< "Usage:"                                              << std::endl
		<< "  Prepare FOAM for Monte-Carlo generation"           << std::endl
		<< "    sidisgen initialize <parameter file>"            << std::endl
		<< "  Generate events"                                   << std::endl
		<< "    sidisgen generate <parameter file>"              << std::endl
		<< "  List parameters used to produce file"              << std::endl
		<< "    sidisgen inspect <output file>"                  << std::endl
		<< "  Merge multiple event files into one"               << std::endl
		<< "    sidisgen merge-soft <output file> <input files>" << std::endl
		<< "    sidisgen merge-hard <output file> <input files>" << std::endl
		<< "  Show parameter file format information"            << std::endl
		<< "    sidisgen help-params"                            << std::endl;
	return SUCCESS;
}

int command_help_params() {
	std::cout
		<< "Parameter file format summary."                    << std::endl
		<< "For more detailed information, see docs."          << std::endl
		<< std::endl
		<< "file.event_out      <ROOT file>"                   << std::endl
		<< "file.write_momenta  <on, off>"                     << std::endl
		<< "file.write_photon   <on, off>"                     << std::endl
		<< "file.foam_out       <ROOT file>"                   << std::endl
		<< std::endl
		<< "mc.gen_nrad    <on, off>"                          << std::endl
		<< "mc.gen_rad     <on, off>"                          << std::endl
		<< "mc.num_events  <integer>"                          << std::endl
		<< "mc.rej_weight  <real in [1, ∞)>"                   << std::endl
		<< "mc.seed        <integer>"                          << std::endl
		<< "mc.seed_init   <integer>"                          << std::endl
		<< std::endl
		<< "setup.beam_energy  <energy (GeV)>"                 << std::endl
		<< "setup.beam         <pid>"                          << std::endl
		<< "setup.target       <pid>"                          << std::endl
		<< "setup.hadron       <pid>"                          << std::endl
		<< "setup.beam_pol     <real in [-1, 1]>"              << std::endl
		<< "setup.target_pol   <vector in unit sphere>"        << std::endl
		<< std::endl
		<< "phys.sf_set          <prokudin, test, ROOT dict.>" << std::endl
		<< "phys.rc_method       <none, approx, exact>"        << std::endl
		<< "phys.mass_threshold  <mass (GeV)>"                 << std::endl
		<< "phys.soft_threshold  <energy (GeV)>"               << std::endl
		<< std::endl
		<< "cut.k_0_bar   <min> <max>"                         << std::endl
		<< "cut.x         <min> <max>"                         << std::endl
		<< "cut.y         <min> <max>"                         << std::endl
		<< "cut.z         <min> <max>"                         << std::endl
		<< "cut.ph_t_sq   <min> <max>"                         << std::endl
		<< "cut.phi_h     <min> <max>"                         << std::endl
		<< "cut.phi       <min> <max>"                         << std::endl
		<< "cut.tau       <min> <max>"                         << std::endl
		<< "cut.phi_k     <min> <max>"                         << std::endl
		<< "cut.R         <min> <max>"                         << std::endl
		<< "cut.Q_sq      <min> <max>"                         << std::endl
		<< "cut.t         <min> <max>"                         << std::endl
		<< "cut.W_sq      <min> <max>"                         << std::endl
		<< "cut.r         <min> <max>"                         << std::endl
		<< "cut.mx_sq     <min> <max>"                         << std::endl
		<< "cut.q_0       <min> <max>"                         << std::endl
		<< "cut.k2_0      <min> <max>"                         << std::endl
		<< "cut.ph_0      <min> <max>"                         << std::endl
		<< "cut.k_0       <min> <max>"                         << std::endl
		<< "cut.theta_q   <min> <max>"                         << std::endl
		<< "cut.theta_k2  <min> <max>"                         << std::endl
		<< "cut.theta_ph  <min> <max>"                         << std::endl
		<< "cut.theta_k   <min> <max>"                         << std::endl;
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

	TArrayD* prime_arr = file.Get<TArrayD>("stats/prime");
	TArrayD* weight_moms_arr = file.Get<TArrayD>("stats/weight_mom");
	TArrayD* weight_max_arr = file.Get<TArrayD>("stats/weight_max");
	TArrayL* num_events_arr = file.Get<TArrayL>("stats/num_events");
	TArrayD* norm_arr = file.Get<TArrayD>("stats/norm");
	if (
			prime_arr == nullptr
			|| weight_moms_arr == nullptr
			|| weight_max_arr == nullptr
			|| num_events_arr == nullptr
			|| norm_arr == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Couldn't find statistics in file '") + file_name
			+ "'.");
	}
	std::cout << "Statistics:" << std::endl;
	flags = std::cout.flags();
	std::cout << std::scientific << std::setprecision(6);
	for (std::size_t type_idx = 0; type_idx < NUM_EVENT_TYPES + 1; ++type_idx) {
		auto count = num_events_arr->At(type_idx);
		if (count == 0) {
			continue;
		}
		if (type_idx == NUM_EVENT_TYPES) {
			std::cout << "\ttotal:" << std::endl;
		} else {
			EventType type = static_cast<EventType>(type_idx);
			std::cout << "\t" << event_type_name(type) << ":" << std::endl;
		}

		Real prime = prime_arr->At(type_idx);
		Real norm = norm_arr->At(type_idx);
		std::array<Real, 4> moms = {
			weight_moms_arr->At(4 * type_idx + 0),
			weight_moms_arr->At(4 * type_idx + 1),
			weight_moms_arr->At(4 * type_idx + 2),
			weight_moms_arr->At(4 * type_idx + 3),
		};
		Real max = weight_max_arr->At(type_idx);
		Stats stats(moms, max, count);
		Real mean = stats.mean();
		Real mean_err = std::sqrt(stats.ratio_var1_to_max2()) * stats.max1()
			/ std::sqrt(stats.count());
		Real xs = prime * mean;
		Real xs_err = prime * mean_err;

		std::cout << "\t\tcount:         " << count << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\tefficiency:    " << mean << " ± " << mean_err << std::endl;
	}
	std::cout.flags(flags);

	return SUCCESS;
}

int command_initialize(char const* params_file_name) {
	// Load parameters.
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

	BuilderParams builder_params;
	std::minstd_rand seed_rnd(*params.seed_init >= 0 ? *params.seed_init : 0);
	std::uniform_int_distribution<Seed> seed_dist;

	// Create FOAM file.
	std::string file_name = *params.foam_file;
	std::cout << "Creating ROOT file '" << file_name << "'." << std::endl;
	TFile file(file_name.c_str(), "CREATE");
	if (file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create file '") + file_name + "'.");
	}
	params.write_root(file);

	// Build FOAMs and write to file.
	std::vector<EventType> gen_types;
	if (*params.gen_nrad) {
		gen_types.push_back(EventType::NRAD);
	}
	if (*params.gen_rad) {
		gen_types.push_back(EventType::RAD);
	}
	for (EventType gen_type : gen_types) {
		char const* gen_name = event_type_name(gen_type);
		char const* gen_key = event_type_short_name(gen_type);
		std::cout << "Building " << gen_name << " FOAM." << std::endl;
		DrawProgressBar progress_reporter;
		builder_params.explore_progress_reporter = &progress_reporter;
		builder_params.tune_progress_reporter = &progress_reporter;
		builder_params.seed = seed_dist(seed_rnd);
		Builder builder(gen_type, builder_params, params, *sf);
		try {
			std::cout << "Exploration phase." << std::endl;
			write_progress_bar(std::cout, 0);
			std::cout << '\r';
			std::cout << std::flush;
			builder.explore();
			write_progress_bar(std::cout, 100);
			std::cout << std::endl;

			std::cout << "Tuning phase." << std::endl;
			write_progress_bar(std::cout, 0);
			std::cout << '\r';
			std::cout << std::flush;
			builder.tune();
			write_progress_bar(std::cout, 100);
			std::cout << std::endl;
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_BUILDING_FOAM,
				std::string("Error while building ") + gen_name + " FOAM: "
				+ e.what());
		}
		Real rel_var_err;
		Real rel_var = builder.rel_var(&rel_var_err);
		std::cout << "Constructed " << gen_name << " FOAM with relative "
			<< "variance " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "Writing " << gen_name << " FOAM to file." << std::endl;
		try {
			std::ostringstream os;
			builder.write(os);
			std::string data = os.str();
			file.WriteObject(&data, gen_key);
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_WRITING_FOAM,
				std::string("Failed to write ") + gen_name + " non-radiative "
				+ "FOAM to file '" + file_name + "': " + e.what());
		}
	}

	std::cout << "Finished!" << std::endl;
	return SUCCESS;
}

int command_generate(char const* params_file_name) {
	// Load parameters.
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
	TFile event_file(params.event_file->c_str(), "CREATE");
	if (event_file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create file '") + *params.event_file + "'.");
	}

	GeneratorParams generator_params;
	std::minstd_rand seed_rnd(*params.seed >= 0 ? *params.seed : 0);
	std::uniform_int_distribution<Seed> seed_dist;

	// Load the event generators from file.
	std::vector<EventType> gen_types;
	if (*params.gen_nrad) {
		gen_types.push_back(EventType::NRAD);
	}
	if (*params.gen_rad) {
		gen_types.push_back(EventType::RAD);
	}
	std::vector<Generator> gens;
	std::string foam_file_name = *params.foam_file;
	std::cout << "Opening FOAM file '" << foam_file_name << "'." << std::endl;
	TFile foam_file(foam_file_name.c_str(), "OPEN");
	for (EventType gen_type : gen_types) {
		char const* gen_name = event_type_name(gen_type);
		char const* gen_key = event_type_short_name(gen_type);
		std::cout << "Loading " << gen_name << " FOAM from file." << std::endl;
		try {
			std::string* data = foam_file.Get<std::string>(gen_key);
			if (data == nullptr) {
				throw std::runtime_error(
					std::string("Couldn't find key '") + gen_key + "' in "
					+ "file '" + foam_file_name + "'.");
			}
			std::istringstream is(*data);
			generator_params.seed = seed_dist(seed_rnd);
			gens.emplace_back(gen_type, generator_params, params, *sf, is);
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_READING_FOAM,
				std::string("Failed to read ") + gen_name + " FOAM from file '"
				+ foam_file_name + "': " + e.what());
		}
	}
	foam_file.Close();

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
	Real weight;
	Real jacobian;
	Real x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R;
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
			std::cout << std::flush;
			update_progress = false;
		}
		// Choose a type of event (ex. radiative or non-radiative) to generate.
		std::size_t choose_event_type;
		if (event_idx == 0) {
			// On the first event, we don't know anything about the total cross-
			// sections, so choose the event type arbitrarily.
			choose_event_type = 0;
		} else {
			// We want to generate events so that the total weights contributed
			// by events of each "type" have the same ratio as the cross-
			// sections of each "type".
			Real ratio_max = 0.;
			for (std::size_t idx = 0; idx < gens.size(); ++idx) {
				Generator const& gen = gens[idx];
				// Equivalent to cross-section divided by total weight.
				Real ratio = gen.prime() / gen.weights().count();
				if (!std::isfinite(ratio)) {
					ratio = std::numeric_limits<Real>::infinity();
				}
				if (ratio > ratio_max) {
					ratio_max = ratio;
					choose_event_type = idx;
				}
			}
		}

		// The event vector can store up to the number of dimensions of any of
		// the FOAMs.
		Real event_vec[9];
		weight = gens[choose_event_type].generate(event_vec);
		type = static_cast<Int_t>(gens[choose_event_type].type());

		// Fill in the branches.
		kin::Kinematics kin;
		kin::KinematicsRad kin_rad;
		switch (gens[choose_event_type].type()) {
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
	TArrayD prime_arr(NUM_EVENT_TYPES + 1);
	TArrayD weight_moms_arr(4 * (NUM_EVENT_TYPES + 1));
	TArrayD weight_max_arr(NUM_EVENT_TYPES + 1);
	TArrayL num_events_arr(NUM_EVENT_TYPES + 1);
	TArrayD norm_arr(NUM_EVENT_TYPES + 1);
	Stats stats_total;
	Real prime_total = 0.;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::scientific << std::setprecision(6);
	for (Generator& gen : gens) {
		Stats stats = gen.weights();
		Real prime = gen.prime();
		stats_total += stats;
		prime_total += prime;
		Real norm = prime / stats.count();
		Real mean = stats.mean();
		Real mean_err = std::sqrt(stats.ratio_var1_to_max2()) * stats.max1()
			/ std::sqrt(stats.count());
		Real xs = prime * mean;
		Real xs_err = prime * mean_err;

		Int_t type_idx = static_cast<Int_t>(gen.type());
		prime_arr.SetAt(prime, type_idx);
		weight_moms_arr.SetAt(stats.ratio_m1_to_max1(), 4 * type_idx + 0);
		weight_moms_arr.SetAt(stats.ratio_m2_to_max2(), 4 * type_idx + 1);
		weight_moms_arr.SetAt(stats.ratio_m3_to_max3(), 4 * type_idx + 2);
		weight_moms_arr.SetAt(stats.ratio_m4_to_max4(), 4 * type_idx + 3);
		weight_max_arr.SetAt(stats.max1(), type_idx);
		num_events_arr.SetAt(stats.count(), type_idx);
		norm_arr.SetAt(norm, type_idx);

		std::cout << "\t" << event_type_name(gen.type()) << " events:" << std::endl;
		std::cout << "\t\tcount:         " << stats.count() << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\tefficiency:    " << mean << " ± " << mean_err << std::endl;
	}

	Real norm = prime_total / stats_total.count();
	Real mean = stats_total.mean();
	Real mean_err = std::sqrt(stats_total.ratio_var1_to_max2()) * stats_total.max1()
		/ std::sqrt(stats_total.count());
	Real xs = prime_total * mean;
	Real xs_err = prime_total * mean_err;

	prime_arr.SetAt(prime_total, NUM_EVENT_TYPES);
	weight_moms_arr.SetAt(stats_total.ratio_m1_to_max1(), 4 * NUM_EVENT_TYPES + 0);
	weight_moms_arr.SetAt(stats_total.ratio_m2_to_max2(), 4 * NUM_EVENT_TYPES + 1);
	weight_moms_arr.SetAt(stats_total.ratio_m3_to_max3(), 4 * NUM_EVENT_TYPES + 2);
	weight_moms_arr.SetAt(stats_total.ratio_m4_to_max4(), 4 * NUM_EVENT_TYPES + 3);
	weight_max_arr.SetAt(stats_total.max1(), NUM_EVENT_TYPES);
	num_events_arr.SetAt(stats_total.count(), NUM_EVENT_TYPES);
	norm_arr.SetAt(norm, NUM_EVENT_TYPES);

	std::cout << "\ttotal:" << std::endl;
		std::cout << "\t\tcount:         " << stats_total.count() << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\tefficiency:    " << mean << " ± " << mean_err << std::endl;

	std::cout.flags(flags);
	stats_dir->WriteObject(&prime_arr, "prime");
	stats_dir->WriteObject(&weight_moms_arr, "weight_mom");
	stats_dir->WriteObject(&weight_max_arr, "weight_max");
	stats_dir->WriteObject(&num_events_arr, "num_events");
	stats_dir->WriteObject(&norm_arr, "norm");

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
		if (file.IsZombie()) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				std::string("File '") + file_name + "' not found.");
		}
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
	Real primes[NUM_EVENT_TYPES + 1];
	Stats stats_total[NUM_EVENT_TYPES + 1];
	first = true;
	for (char const* file_name : file_names) {
		// TODO: Ensure that the merged statistics come from the same underlying
		// FOAM.
		TFile file(file_name, "OPEN");
		if (file.IsZombie()) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				std::string("File '") + file_name + "' not found.");
		}
		TArrayD* prime_arr = file.Get<TArrayD>("stats/prime");
		TArrayD* weight_moms_arr = file.Get<TArrayD>("stats/weight_mom");
		TArrayD* weight_max_arr = file.Get<TArrayD>("stats/weight_max");
		TArrayL* num_events_arr = file.Get<TArrayL>("stats/num_events");
		TArrayD* norm_arr = file.Get<TArrayD>("stats/norm");
		if (
				prime_arr == nullptr
				|| weight_moms_arr == nullptr
				|| weight_max_arr == nullptr
				|| num_events_arr == nullptr
				|| norm_arr == nullptr) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				std::string("Couldn't find statistics in file '") + file_name
				+ "'.");
		}
		for (std::size_t type_idx = 0; type_idx < NUM_EVENT_TYPES + 1; ++type_idx) {
			Real prime = prime_arr->At(type_idx);
			auto count = num_events_arr->At(type_idx);
			std::array<Real, 4> moms = {
				weight_moms_arr->At(4 * type_idx + 0),
				weight_moms_arr->At(4 * type_idx + 1),
				weight_moms_arr->At(4 * type_idx + 2),
				weight_moms_arr->At(4 * type_idx + 3),
			};
			Real max = weight_max_arr->At(type_idx);
			Stats stats(moms, max, count);
			stats_total[type_idx] += stats;
			if (first) {
				primes[type_idx] = prime;
			} else {
				if (primes[type_idx] != prime) {
					throw Exception(
						ERROR_FOAM_INCOMPATIBLE,
						std::string("FOAM from file '") + file_name
						+ "' has incompatible prime.");
				}
			}
			first = false;
		}
	}
	TDirectory* stats_dir = file_out.mkdir("stats");
	if (stats_dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create directory 'stats' in ROOT file."));
	}
	stats_dir->cd();
	TArrayD prime_arr_out(NUM_EVENT_TYPES + 1);
	TArrayD weight_moms_arr_out(NUM_EVENT_TYPES + 1);
	TArrayD weight_max_arr_out(NUM_EVENT_TYPES + 1);
	TArrayL num_events_arr_out(NUM_EVENT_TYPES + 1);
	TArrayD norm_arr_out(NUM_EVENT_TYPES + 1);
	for (std::size_t type_idx = 0; type_idx < NUM_EVENT_TYPES + 1; ++type_idx) {
		Real norm = primes[type_idx] / stats_total[type_idx].count();
		prime_arr_out.SetAt(primes[type_idx], type_idx);
		weight_moms_arr_out.SetAt(stats_total[type_idx].ratio_m1_to_max1(), 4 * type_idx + 0);
		weight_moms_arr_out.SetAt(stats_total[type_idx].ratio_m2_to_max2(), 4 * type_idx + 1);
		weight_moms_arr_out.SetAt(stats_total[type_idx].ratio_m3_to_max3(), 4 * type_idx + 2);
		weight_moms_arr_out.SetAt(stats_total[type_idx].ratio_m4_to_max4(), 4 * type_idx + 3);
		weight_max_arr_out.SetAt(stats_total[type_idx].max1(), type_idx);
		num_events_arr_out.SetAt(stats_total[type_idx].count(), type_idx);
		norm_arr_out.SetAt(norm, type_idx);
	}
	stats_dir->WriteObject(&prime_arr_out, "prime");
	stats_dir->WriteObject(&weight_moms_arr_out, "weight_mom");
	stats_dir->WriteObject(&weight_max_arr_out, "weight_max");
	stats_dir->WriteObject(&num_events_arr_out, "num_events");
	stats_dir->WriteObject(&norm_arr_out, "norm");

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
			std::cout << "Try `sidisgen help`." << std::endl;
			return SUCCESS;
		}
		std::string command = argv[1];
		if (command[0] == '-' && command[1] == '-' && command[2] != '-') {
			command = command.substr(2);
		}
		if (command == "help" || command == "-?") {
			return command_help();
		} else if (command == "help-params") {
			return command_help_params();
		} else if (command == "version" || command == "-v") {
			return command_version();
		} else if (command == "inspect" || command == "-i") {
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
		} else if (command == "initialize" || command == "-i") {
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
		} else if (command == "generate" || command == "-g") {
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
		} else if (command == "merge-soft") {
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
		} else if (command == "merge-hard") {
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

