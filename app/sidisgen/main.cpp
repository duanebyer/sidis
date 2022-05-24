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

int const OUTPUT_STATS_PRECISION = 3;

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

// Estimates the quantity `sqrt(E[X^2])`, including first-order correction.
Real est_sqrt_m2(Stats const& stats, Real* err_out=nullptr) {
	Real ratio_m2_to_max2 = stats.ratio_m2_to_max2();
	if (ratio_m2_to_max2 < 0.) {
		ratio_m2_to_max2 = 0.;
	}
	Real est_0 = stats.max() * std::sqrt(ratio_m2_to_max2);
	Real correction = (1. / 8) * stats.ratio_var2_to_m2p2() / stats.count();
	if (!(-0.5 < correction && correction < 0.5)) {
		correction = 0.;
	}
	if (err_out != nullptr) {
		Real ratio_var2_to_max2_m2 = stats.ratio_var2_to_max2_m2();
		if (ratio_var2_to_max2_m2 < 0.) {
			ratio_var2_to_max2_m2 = 0.;
		}
		*err_out = stats.max() * std::sqrt(
			0.25 * ratio_var2_to_max2_m2 / stats.count());
	}
	Real est = est_0 * (1. + correction);
	return est;
}
Real est_mean(Stats const& stats, Real* err_out=nullptr) {
	Real est = stats.mean();
	Real ratio_var1_to_max2 = stats.ratio_var1_to_max2();
	if (ratio_var1_to_max2 < 0.) {
		ratio_var1_to_max2 = 0.;
	}
	if (err_out != nullptr) {
		*err_out = stats.max() * std::sqrt(ratio_var1_to_max2 / stats.count());
	}
	return est;
}
Real est_rel_var(Stats const& stats, Real* err_out=nullptr) {
	Real est_0 = stats.ratio_var1_to_m1p2();
	Real correction = (2. * stats.ratio_skew1_to_var1_m1() - 3. * est_0)
		/ stats.count();
	if (!(-0.5 < correction && correction < 0.5)) {
		correction = 0.;
	}
	Real est = est_0 * (1. + correction);
	if (err_out != nullptr) {
		Real kurt_term = stats.ratio_kurt1_to_var1p2();
		Real skew_term = stats.ratio_skew1_to_var1_m1();
		Real var_term = stats.ratio_var1_to_m1p2();
		Real term = kurt_term - 4. * skew_term + 4. * var_term - 1.;
		if (term < 0.) {
			term = 0.;
		}
		*err_out = est * std::sqrt(term / stats.count());
	}
	return est;
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
		<< "file.event_out        <ROOT file>"                 << std::endl
		<< "file.write_momenta    <on, off>"                   << std::endl
		<< "file.write_photon     <on, off>"                   << std::endl
		<< "file.write_sf_set     <on, off>"                   << std::endl
		<< "file.write_mc_coords  <on, off>"                   << std::endl
		<< "file.foam_out         <ROOT file>"                 << std::endl
		<< std::endl
		<< "mc.nrad.gen              <on, off>"                << std::endl
		<< "mc.nrad.init.seed        <integer>"                << std::endl
		<< "mc.nrad.init.max_cells   <integer>"                << std::endl
		<< "mc.nrad.init.target_eff  <real in [0, 1]>"         << std::endl
		<< "mc.nrad.init.scale_exp   <real>"                   << std::endl
		<< "mc.rad.gen               <on, off>"                << std::endl
		<< "mc.rad.init.seed         <integer>"                << std::endl
		<< "mc.rad.init.max_cells    <integer>"                << std::endl
		<< "mc.rad.init.target_eff   <real in [0, 1]>"         << std::endl
		<< "mc.rad.init.scale_exp    <real>"                   << std::endl
		<< "mc.num_events            <integer>"                << std::endl
		<< "mc.rej_weight            <real in [1, ∞)>"         << std::endl
		<< "mc.seed                  <integer>"                << std::endl
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
		<< "cut.k_0_bar       <min> <max>"                     << std::endl
		<< "cut.x             <min> <max>"                     << std::endl
		<< "cut.y             <min> <max>"                     << std::endl
		<< "cut.z             <min> <max>"                     << std::endl
		<< "cut.ph_t_sq       <min> <max>"                     << std::endl
		<< "cut.phi_h         <min> <max>"                     << std::endl
		<< "cut.phi           <min> <max>"                     << std::endl
		<< "cut.tau           <min> <max>"                     << std::endl
		<< "cut.phi_k         <min> <max>"                     << std::endl
		<< "cut.R             <min> <max>"                     << std::endl
		<< "cut.Q_sq          <min> <max>"                     << std::endl
		<< "cut.t             <min> <max>"                     << std::endl
		<< "cut.W_sq          <min> <max>"                     << std::endl
		<< "cut.r             <min> <max>"                     << std::endl
		<< "cut.mx_sq         <min> <max>"                     << std::endl
		<< "cut.qt_to_Q       <min> <max>"                     << std::endl
		<< "cut.lab_mom_q     <min> <max>"                     << std::endl
		<< "cut.lab_mom_k2    <min> <max>"                     << std::endl
		<< "cut.lab_mom_h     <min> <max>"                     << std::endl
		<< "cut.lab_mom_k     <min> <max>"                     << std::endl
		<< "cut.lab_theta_q   <min> <max>"                     << std::endl
		<< "cut.lab_theta_k2  <min> <max>"                     << std::endl
		<< "cut.lab_theta_h   <min> <max>"                     << std::endl
		<< "cut.lab_theta_k   <min> <max>"                     << std::endl;
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
			|| prime_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| weight_moms_arr == nullptr
			|| weight_moms_arr->GetSize() != 4 * (NUM_EVENT_TYPES + 1)
			|| weight_max_arr == nullptr
			|| weight_max_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| num_events_arr == nullptr
			|| num_events_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| norm_arr == nullptr
			|| norm_arr->GetSize() != NUM_EVENT_TYPES + 1) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Couldn't find statistics in file '") + file_name
			+ "'.");
	}
	std::cout << "Statistics:" << std::endl;
	flags = std::cout.flags();
	std::cout << std::scientific << std::setprecision(OUTPUT_STATS_PRECISION);
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
		Real mean_err;
		Real mean = est_mean(stats, &mean_err);
		Real xs = prime * mean;
		Real xs_err = prime * mean_err;
		Real rel_var_err;
		Real rel_var = est_rel_var(stats, &rel_var_err);

		std::cout << "\t\tcount:         " << count << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tprime:         " << prime << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "\t\tefficiency:    " << 1. / std::sqrt(1. + rel_var) << std::endl;
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
		params.fill_defaults();
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

	BuilderReporters builder_reporters;

	// Create FOAM file.
	std::string file_name = *params.foam_file;
	std::cout << "Creating ROOT file '" << file_name << "'." << std::endl;
	TFile file(file_name.c_str(), "CREATE");
	if (file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create file '") + file_name + "'.");
	}

	// Build FOAMs and write to file.
	std::vector<EventType> gen_types;
	if (*params.nrad_gen) {
		gen_types.push_back(EventType::NRAD);
	}
	if (*params.rad_gen) {
		gen_types.push_back(EventType::RAD);
	}
	for (EventType gen_type : gen_types) {
		char const* gen_name = event_type_name(gen_type);
		char const* gen_key = event_type_short_name(gen_type);
		std::cout << "Building " << gen_name << " FOAM." << std::endl;
		DrawProgressBar progress_reporter;
		builder_reporters.explore_progress = &progress_reporter;
		builder_reporters.tune_progress = &progress_reporter;
		Builder builder(gen_type, builder_reporters, params, *sf);
		// Overwrite the seed parameter with the chosen seed.
		if (gen_type == EventType::NRAD) {
			params.nrad_seed_init.reset(builder.seed());
		} else if (gen_type == EventType::RAD) {
			params.rad_seed_init.reset(builder.seed());
		} else if (gen_type == EventType::EXCL) {
			// Not yet implemented.
		}
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
		std::cout << "Constructed " << gen_name << " FOAM with size "
			<< builder.size() << "." << std::endl;
		std::cout << "\tRelative variance: " << rel_var
			<< " ± " << rel_var_err << std::endl;
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
	params.write_root(file);

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
		params.fill_defaults();
		if (params.seed->size() != 1) {
			throw std::runtime_error("Exactly one seed must be provided.");
		}
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

	// Load the event generators from file.
	std::vector<EventType> gen_types;
	if (*params.nrad_gen) {
		gen_types.push_back(EventType::NRAD);
	}
	if (*params.rad_gen) {
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
			gens.emplace_back(gen_type, params, *sf, is);
			// Overwrite the seed parameter with the actual seed that was chosen
			// in the case that the seed was asked to be chosen randomly.
			params.seed.reset({ gens.back().seed() });
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_READING_FOAM,
				std::string("Failed to read ") + gen_name + " FOAM from file '"
				+ foam_file_name + "': " + e.what());
		}
	}
	foam_file.Close();

	// Setup initial conditions.
	part::Nucleus target = *params.target;
	part::Lepton beam = *params.beam;
	part::Hadron hadron = *params.hadron;
	math::Vec3 target_pol = *params.target_pol;
	part::Particles ps(target, beam, hadron, *params.Mth);
	Real S = 2.*(*params.beam_energy)*ps.M;

	kin::Initial init(ps, *params.beam_energy);
	ULong_t N_gen = *params.num_events >= 0 ? *params.num_events : 0;

	// Prepare branches in output ROOT file.
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
	if (*params.write_momenta) {
		events.Branch("p", "TLorentzVector", &p);
		events.Branch("k1", "TLorentzVector", &k1);
		events.Branch("q", "TLorentzVector", &q);
		events.Branch("k2", "TLorentzVector", &k2);
		events.Branch("ph", "TLorentzVector", &ph);
		if (*params.rad_gen && *params.write_photon) {
			events.Branch("k", "TLorentzVector", &k);
		}
	}
	sf::SfXX sf_out;
	if (*params.write_sf_set) {
		// TODO: Right now, this is depending on the `SfXX` structure having a
		// very specific format. This isn't guaranteed to be true in the future,
		// if things get reorganized. Not sure what a better approach is, as
		// ROOT doesn't have a better way of writing plain C-structs into trees.
		events.Branch("sf", &sf_out,
			"F_UUL/D:F_UUT/D:F_UU_cos_phih/D:F_UU_cos_2phih/D:F_UL_sin_phih/D:F_UL_sin_2phih/D:F_UTL_sin_phih_m_phis/D:F_UTT_sin_phih_m_phis/D:F_UT_sin_2phih_m_phis/D:F_UT_sin_3phih_m_phis/D:F_UT_sin_phis/D:F_UT_sin_phih_p_phis/D:F_LU_sin_phih/D:F_LL/D:F_LL_cos_phih/D:F_LT_cos_phih_m_phis/D:F_LT_cos_2phih_m_phis/D:F_LT_cos_phis/D");
	}
	Real mc_coords[9];
	Real ph_coords[9];
	if (*params.write_mc_coords) {
		events.Branch("mc_coords", &mc_coords, "mc_coords[9]/D");
	}
	// Write parameter file.
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
				Stats stats = gen.weights();
				// Equivalent to cross-section divided by total weight.
				Real ratio = gen.prime() * est_sqrt_m2(stats) / stats.count();
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
		weight = gens[choose_event_type].generate(ph_coords, mc_coords);
		type = static_cast<Int_t>(gens[choose_event_type].type());

		// Fill in the branches.
		kin::Kinematics kin;
		kin::KinematicsRad kin_rad;
		switch (gens[choose_event_type].type()) {
		case EventType::NRAD:
			// Non-radiative event.
			{
				// Fill in branches.
				x = ph_coords[0];
				y = ph_coords[1];
				z = ph_coords[2];
				ph_t_sq = ph_coords[3];
				phi_h = ph_coords[4];
				phi = ph_coords[5];
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
					k = TLorentzVector();
				}
				if (*params.write_sf_set) {
					sf_out = sf->sf(hadron, x, z, S * x * y, ph_t_sq);
				}
			}
			break;
		case EventType::RAD:
			// Radiative event.
			{
				// Fill in branches.
				x = ph_coords[0];
				y = ph_coords[1];
				z = ph_coords[2];
				ph_t_sq = ph_coords[3];
				phi_h = ph_coords[4];
				phi = ph_coords[5];
				tau = ph_coords[6];
				phi_k = ph_coords[7];
				R = ph_coords[8];
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
				if (*params.write_sf_set) {
					sf_out = sf->sf(hadron, x, z, S * x * y, ph_t_sq);
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
	std::size_t count_total = 0;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::scientific << std::setprecision(OUTPUT_STATS_PRECISION);
	for (Generator& gen : gens) {
		// Need these ahead of time, for computing the weight.
		prime_total += gen.prime();
		count_total += gen.weights().count();
	}
	for (Generator& gen : gens) {
		Stats stats = gen.weights();
		Real prime = gen.prime();
		std::size_t count = stats.count();
		Real mean_err;
		Real mean = est_mean(stats, &mean_err);
		Real xs = prime * mean;
		Real xs_err = prime * mean_err;
		Real rel_var_err;
		Real rel_var = est_rel_var(stats, &rel_var_err);

		// The weight is chosen so that the total cross-section averages out
		// correctly. Normally it is near one for all types of events.
		Real weight = (prime / count) / (prime_total / count_total);
		Real norm = prime / count;
		stats_total += weight * stats;

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
		std::cout << "\t\tweight:        " << weight << std::endl;
		std::cout << "\t\tcount:         " << stats.count() << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tprime:         " << prime << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "\t\tefficiency:    " << 1. / std::sqrt(1. + rel_var) << std::endl;
	}

	Real norm = prime_total / count_total;
	Real mean_err;
	Real mean = est_mean(stats_total, &mean_err);
	Real xs = prime_total * mean;
	Real xs_err = prime_total * mean_err;
	Real rel_var_err;
	Real rel_var = est_rel_var(stats_total, &rel_var_err);

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
		std::cout << "\t\tprime:         " << prime_total << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "\t\tefficiency:    " << 1. / std::sqrt(1. + rel_var) << std::endl;

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
	Real primes[NUM_EVENT_TYPES + 1] = {};
	Stats stats_total[NUM_EVENT_TYPES + 1] = {};
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
				|| prime_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| weight_moms_arr == nullptr
				|| weight_moms_arr->GetSize() != 4 * (NUM_EVENT_TYPES + 1)
				|| weight_max_arr == nullptr
				|| weight_max_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| num_events_arr == nullptr
				|| num_events_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| norm_arr == nullptr
				|| norm_arr->GetSize() != NUM_EVENT_TYPES + 1) {
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
		}
		first = false;
	}
	TDirectory* stats_dir = file_out.mkdir("stats");
	if (stats_dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't create directory 'stats' in ROOT file."));
	}
	stats_dir->cd();
	TArrayD prime_arr_out(NUM_EVENT_TYPES + 1);
	TArrayD weight_moms_arr_out(4 * (NUM_EVENT_TYPES + 1));
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

