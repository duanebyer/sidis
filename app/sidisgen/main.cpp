#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <random>
#include <regex>
#include <string>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

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
#include "exception.hpp"
#include "params.hpp"
#include "params_format.hpp"
#include "terminal.hpp"
#include "utility.hpp"

using namespace sidis;

static_assert(
	std::is_same<Real, Double>::value,
	"ROOT `Double_t` type must be same as libsidis `Real` type.");

namespace {

int const OUTPUT_STATS_PRECISION = 3;

struct DrawProgressBar final :
		public bubble::ExploreProgressReporter<Double>,
		public bubble::TuneProgressReporter<Double> {
	void operator()(bubble::ExploreProgress<Double> progress) override {
		Double percent = 100. * progress.progress;
		write_progress_bar(std::cout, static_cast<unsigned>(percent));
		std::cout << '\r';
		std::cout << std::flush;
	}
	void operator()(bubble::TuneProgress<Double> progress) override {
		Double percent = 100. * progress.progress;
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
		Params& params,
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
	std::string sf_set_name = params["phys.sf_set"].any();
	std::stringstream ss(sf_set_name);
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
		sf_out->reset(
			new sf::set::TestSfSet(params["setup.target"].any()));
	} else {
		// TODO: Make this work on Windows as well, and make providing the
		// extension optional.
		std::string file_name = base + ".so";
		if (gSystem->Load(file_name.c_str()) != 0) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				"Failed to load structure function from shared library file '"
				+ file_name + "'.");
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
				<< sf_set_name << "'." << std::endl;
			sf::GaussianWwTmdSet* tmd_ptr = static_cast<sf::GaussianWwTmdSet*>(sf_class->New());
			tmd_out->reset(tmd_ptr);
			sf_out->reset(new sf::GaussianWwTmdSfSet(*tmd_ptr));
		} else {
			throw Exception(
				ERROR_STRUCTURE_FUNCTIONS_NOT_FOUND,
				"Couldn't find structure functions in file '" + sf_set_name
				+ ".so'.");
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
					"Cannot filter on structure function index "
					+ std::to_string(idx) + " because out of bounds.");
			}
			bool select_mask[sf::set::NUM_SF] = { false };
			select_mask[idx] = true;
			zip_and(mask, mask + sf::set::NUM_SF, select_mask);
		} else {
			throw Exception(
				ERROR_STRUCTURE_FUNCTIONS_PARSE,
				"Unrecognized structure function filter '" + part + "'.");
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
		<< "    sidisgen help params"                            << std::endl;
	return SUCCESS;
}

int command_help_params() {
	Params params = PARAMS_STD_FORMAT;
	std::cout
		<< "Parameter file format summary."                      << std::endl
		<< "For more details, try `sidisgen help <param name>`." << std::endl
		<< std::endl;
	std::vector<std::string> table_entries;
	for (std::string const& name : params.names()) {
		table_entries.push_back(name);
		table_entries.push_back(params.usage(name.c_str()));
	}
	std::vector<unsigned> table_widths = { 24, 42 };
	write_table(std::cout, table_entries, table_widths, 4);
	return SUCCESS;
}

int command_help_param(std::string param_name) {
	Params params = PARAMS_STD_FORMAT;
	if (params.names().count(param_name) != 0) {
		std::vector<std::string> table_entries {
			"Summary", params.brief(param_name),
			"Key", param_name,
			"Value", params.usage(param_name),
			"Info", params.doc(param_name),
		};
		std::vector<unsigned> table_widths = { 7, 59 };
		write_table(std::cout, table_entries, table_widths, 4);
		return SUCCESS;
	} else {
		return command_help();
	}
}

int command_version() {
	std::cout << "sidisgen "
		<< SIDIS_VERSION_MAJOR << "."
		<< SIDIS_VERSION_MINOR << "."
		<< SIDIS_VERSION_PATCH << "."
		<< SIDIS_VERSION_TWEAK << std::endl;
	return SUCCESS;
}

int command_inspect(std::string file_name) {
	Params params = PARAMS_STD_FORMAT;
	TFile file(file_name.c_str(), "OPEN");
	if (file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"File '" + file_name + "' not found.");
	}
	params.read_root(file);
	std::cout << "Parameters:" << std::endl;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::setprecision(std::numeric_limits<Double>::digits10 + 1);
	params.write_stream(std::cout);
	std::cout.flags(flags);
	std::cout << std::endl;

	RootArrayD* prime_arr = file.Get<RootArrayD>("stats/prime");
	RootArrayD* weight_moms_arr = file.Get<RootArrayD>("stats/weight_mom");
	RootArrayD* weight_max_arr = file.Get<RootArrayD>("stats/weight_max");
	RootArrayD* num_events_arr = file.Get<RootArrayD>("stats/num_events");
	RootArrayD* num_events_acc_arr = file.Get<RootArrayD>("stats/num_events_acc");
	RootArrayD* norm_arr = file.Get<RootArrayD>("stats/norm");
	if (
			prime_arr == nullptr
			|| prime_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| weight_moms_arr == nullptr
			|| weight_moms_arr->GetSize() != 4 * (NUM_EVENT_TYPES + 1)
			|| weight_max_arr == nullptr
			|| weight_max_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| num_events_arr == nullptr
			|| num_events_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| num_events_acc_arr == nullptr
			|| num_events_acc_arr->GetSize() != NUM_EVENT_TYPES + 1
			|| norm_arr == nullptr
			|| norm_arr->GetSize() != NUM_EVENT_TYPES + 1) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Couldn't find statistics in file '" + file_name + "'.");
	}
	std::cout << "Statistics:" << std::endl;
	flags = std::cout.flags();
	std::cout << std::scientific << std::setprecision(OUTPUT_STATS_PRECISION);
	for (std::size_t ev_idx = 0; ev_idx < NUM_EVENT_TYPES + 1; ++ev_idx) {
		auto count = num_events_arr->At(ev_idx);
		auto count_acc = num_events_acc_arr->At(ev_idx);
		if (count == 0) {
			continue;
		}
		if (ev_idx == NUM_EVENT_TYPES) {
			std::cout << "\ttotal:" << std::endl;
		} else {
			EventType ev_type = static_cast<EventType>(ev_idx);
			std::cout << "\t" << event_type_name(ev_type) << ":" << std::endl;
		}

		Double prime = prime_arr->At(ev_idx);
		Double norm = norm_arr->At(ev_idx);
		Double acceptance = static_cast<Double>(count_acc) / count;
		std::array<Double, 4> moms = {
			weight_moms_arr->At(4 * ev_idx + 0),
			weight_moms_arr->At(4 * ev_idx + 1),
			weight_moms_arr->At(4 * ev_idx + 2),
			weight_moms_arr->At(4 * ev_idx + 3),
		};
		Double max = weight_max_arr->At(ev_idx);

		Stats stats(moms, max, count);
		Stats stats_acc = stats;
		stats_acc.rescale_count(count_acc);
		Double mean_err;
		Double mean = stats.est_mean(&mean_err);
		Double xs = prime * mean;
		Double xs_err = prime * mean_err;
		Double rel_var_err;
		Double rel_var = stats_acc.est_rel_var(&rel_var_err);

		std::cout << "\t\tcount:         " << count_acc << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tprime:         " << prime << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\tacceptance:    " << acceptance << std::endl;
		std::cout << "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "\t\tefficiency:    " << 1. / std::sqrt(1. + rel_var) << std::endl;
	}
	std::cout.flags(flags);

	return SUCCESS;
}

int command_initialize(std::string params_file_name) {
	// Load parameters.
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Parameter file '" + params_file_name + "' not found.");
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params = PARAMS_STD_FORMAT;
	try {
		params.read_stream(params_file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_PARSE,
			"Failed to parse parameter file '" + params_file_name + "': "
			+ e.what());
	}
	std::cout << std::endl;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::setprecision(std::numeric_limits<Double>::digits10 + 1);
	params.write_stream(std::cout);
	std::cout.flags(flags);
	std::cout << std::endl;

	// Load the structure functions.
	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	alloc_sf(params, &sf, &tmd);

	BuilderReporters builder_reporters;

	// Create FOAM file.
	std::string file_name = params["file.foam"].any();
	std::cout << "Creating ROOT file '" << file_name << "'." << std::endl;
	TFile file(file_name.c_str(), "CREATE");
	if (file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Couldn't create file '" + file_name + "'.");
	}

	// Build FOAMs and write to file.
	std::vector<EventType> ev_types = p_enabled_event_types(params);
	std::random_device rnd_dev;
	for (EventType ev_type : ev_types) {
		// Choose specific seeds, if they haven't already been supplied.
		if (!params.is_set(p_name_init_seed(ev_type))) {
			params.set(p_name_init_seed(ev_type), new ValueInt(rnd_dev()));
		}
	}
	// Use a queue here so that the builders can be easily destructed in a FIFO
	// order after they have been initialized and written to file.
	std::queue<Builder> builders;
	for (EventType ev_type : ev_types) {
		char const* ev_name = event_type_name(ev_type);
		std::cout << "Building " << ev_name << " FOAM." << std::endl;
		DrawProgressBar progress_reporter;
		builder_reporters.explore_progress = &progress_reporter;
		builder_reporters.tune_progress = &progress_reporter;
		builders.emplace(ev_type, builder_reporters, params, *sf);
	}
	// Check that all provided parameters were used.
	try {
		params.filter("init"_F).check_complete();
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			"Invalid parameter file '" + params_file_name + "': " + e.what());
	}

	while (!builders.empty()) {
		EventType ev_type = builders.front().ev_type();
		std::string ev_name = event_type_name(ev_type);
		std::string ev_key = event_type_short_name(ev_type);
		try {
			std::cout << "Exploration phase." << std::endl;
			write_progress_bar(std::cout, 0);
			std::cout << '\r';
			std::cout << std::flush;
			builders.front().explore();
			write_progress_bar(std::cout, 100);
			std::cout << std::endl;

			std::cout << "Tuning phase." << std::endl;
			write_progress_bar(std::cout, 0);
			std::cout << '\r';
			std::cout << std::flush;
			builders.front().tune();
			write_progress_bar(std::cout, 100);
			std::cout << std::endl;
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_BUILDING_FOAM,
				"Error while building " + ev_name + " FOAM: " + e.what());
		}
		Double rel_var_err;
		Double rel_var = builders.front().rel_var(&rel_var_err);
		std::cout << "Constructed " << ev_name << " FOAM with size "
			<< builders.front().size() << "." << std::endl;
		std::cout << "\tRelative variance: " << rel_var
			<< " ± " << rel_var_err << std::endl;
		std::cout << "Writing " << ev_name << " FOAM to file." << std::endl;
		try {
			std::ostringstream os;
			std::size_t hash = builders.front().write(os);
			std::string data = os.str();
			params.set(p_name_init_hash(ev_type), new ValueSize(hash));
			file.WriteObject(&data, ev_key.c_str());
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_WRITING_FOAM,
				"Failed to write " + ev_name + " non-radiative FOAM to file '"
				+ file_name + "': " + e.what());
		}
		builders.pop();
	}
	params.write_root(file);

	std::cout << "Finished!" << std::endl;
	return SUCCESS;
}

int command_generate(std::string params_file_name) {
	// Load parameters.
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Parameter file '" + params_file_name + "' not found.");
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params = PARAMS_STD_FORMAT;
	try {
		params.read_stream(params_file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_PARSE,
			"Failed to parse parameter file '" + params_file_name + "': "
			+ e.what());
	}
	std::cout << std::endl;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::setprecision(std::numeric_limits<Double>::digits10 + 1);
	params.write_stream(std::cout);
	std::cout.flags(flags);
	std::cout << std::endl;

	// Load the structure functions.
	std::unique_ptr<sf::SfSet> sf;
	std::unique_ptr<sf::TmdSet> tmd;
	alloc_sf(params, &sf, &tmd);

	// Load the event generators from file.
	std::vector<EventType> ev_types = p_enabled_event_types(params);
	std::vector<Generator> gens;
	std::string foam_file_name = params["file.foam"].any();
	std::cout << "Opening FOAM file '" << foam_file_name << "'." << std::endl;
	TFile foam_file(foam_file_name.c_str(), "OPEN");
	if (foam_file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Couldn't find file '" + foam_file_name + "'.");
	}
	Params params_foam = PARAMS_STD_FORMAT;
	params_foam.read_root(foam_file);
	std::random_device rnd_dev;
	if (!params.is_set("mc.seed")) {
		params.set("mc.seed", new ValueSeedGen(rnd_dev()));
	} else if (params.get<ValueSeedGen>("mc.seed").val.seeds.size() != 1) {
		throw std::runtime_error("Exactly one seed must be provided.");
	}
	check_can_provide_foam(params_foam, params);
	for (EventType ev_type : ev_types) {
		std::string ev_name = event_type_name(ev_type);
		std::string ev_key = event_type_short_name(ev_type);
		std::cout << "Loading " << ev_name << " FOAM from file." << std::endl;
		try {
			std::string* data = foam_file.Get<std::string>(ev_key.c_str());
			if (data == nullptr) {
				throw std::runtime_error(
					"Couldn't find key '" + ev_key + "' in file '"
					+ foam_file_name + "'.");
			}
			std::istringstream is(*data);
			gens.emplace_back(ev_type, params, *sf, is);
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_READING_FOAM,
				"Failed to read " + ev_name + " FOAM from file '"
				+ foam_file_name + "': " + e.what());
		}
	}
	foam_file.Close();

	std::string event_file_name = params["file.event"].any();
	std::cout
		<< "Opening event output file '" << event_file_name << "'."
		<< std::endl;
	TFile event_file(event_file_name.c_str(), "CREATE");
	if (event_file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Couldn't create file '" + event_file_name + "'.");
	}

	// Setup initial conditions.
	part::Nucleus target = params["setup.target"].any();
	part::Lepton beam = params["setup.beam"].any();
	part::Hadron hadron = params["setup.hadron"].any();
	math::Vec3 target_pol = params["setup.target_pol"].any();
	Double beam_energy = params["setup.beam_energy"].any();
	Double M_th = params["phys.mass_threshold"].any();
	part::Particles ps(target, beam, hadron, M_th);
	Double S = 2.*beam_energy*ps.M;

	kin::Initial init(ps, beam_energy);
	std::size_t num_events = std::max<std::size_t>(0, params["mc.num_events"].any());

	// Prepare branches in output ROOT file.
	event_file.cd();
	TTree events("events", "events");
	Int_t ev_idx;
	Double weight;
	Double jacobian;
	Double x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R;
	TLorentzVector p, k1, q, k2, ph, k;
	events.Branch("type", &ev_idx);
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
	bool write_momenta = params["file.write_momenta"].any();
	bool write_photon = false;
	bool write_sf_set = params["file.write_sf_set"].any();
	bool write_mc_coords = params["file.write_mc_coords"].any();
	if (write_momenta) {
		events.Branch("p", "TLorentzVector", &p);
		events.Branch("k1", "TLorentzVector", &k1);
		events.Branch("q", "TLorentzVector", &q);
		events.Branch("k2", "TLorentzVector", &k2);
		events.Branch("ph", "TLorentzVector", &ph);
		if (params["mc.rad.enable"].any()) {
			write_photon = params["file.write_photon"].any();
			if (write_photon) {
				events.Branch("k", "TLorentzVector", &k);
			}
		}
	}
	sf::SfXX sf_out;
	if (write_sf_set) {
		// TODO: Right now, this is depending on the `SfXX` structure having a
		// very specific format. This isn't guaranteed to be true in the future,
		// if things get reorganized. Not sure what a better approach is, as
		// ROOT doesn't have a better way of writing plain C-structs into trees.
		events.Branch("sf", &sf_out,
			"F_UUL/D:F_UUT/D:F_UU_cos_phih/D:F_UU_cos_2phih/D:F_UL_sin_phih/D:F_UL_sin_2phih/D:F_UTL_sin_phih_m_phis/D:F_UTT_sin_phih_m_phis/D:F_UT_sin_2phih_m_phis/D:F_UT_sin_3phih_m_phis/D:F_UT_sin_phis/D:F_UT_sin_phih_p_phis/D:F_LU_sin_phih/D:F_LL/D:F_LL_cos_phih/D:F_LT_cos_phih_m_phis/D:F_LT_cos_2phih_m_phis/D:F_LT_cos_phis/D");
	}
	Double mc_coords[9];
	Double ph_coords[9];
	if (write_mc_coords) {
		events.Branch("mc_coords", &mc_coords, "mc_coords[9]/D");
	}
	// Check that all provided parameters were used.
	try {
		params.filter("gen"_F).check_complete();
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			"Invalid parameter file '" + params_file_name + "': " + e.what());
	}
	// Write parameter file.
	params.write_root(event_file);

	std::cout << "Generating events." << std::endl;
	bool update_progress = true;
	std::size_t percent = 0;
	std::size_t next_percent_rem = num_events % 100;
	std::size_t next_percent = num_events / 100;
	for (std::size_t event_idx = 0; event_idx < num_events; ++event_idx) {
		while (event_idx >= next_percent + (next_percent_rem != 0)) {
			percent += 1;
			next_percent = math::prod_div(num_events, percent + 1, 100, next_percent_rem);
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
			Double ratio_max = 0.;
			for (std::size_t idx = 0; idx < gens.size(); ++idx) {
				Generator const& gen = gens[idx];
				Stats stats = gen.weights();
				Stats stats_acc = gen.weights_acc();
				// Equivalent to cross-section divided by total weight.
				Double ratio = gen.prime() * stats_acc.est_sqrt_m2() / stats.count();
				if (!std::isfinite(ratio)) {
					ratio = std::numeric_limits<Double>::infinity();
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
		EventType ev_type = gens[choose_event_type].ev_type();
		ev_idx = static_cast<Int_t>(ev_type);

		// Fill in the branches.
		kin::Kinematics kin;
		kin::KinematicsRad kin_rad;
		switch (gens[choose_event_type].ev_type()) {
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
				if (write_momenta) {
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
				if (write_sf_set) {
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
				if (write_momenta) {
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
				if (write_sf_set) {
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
			"Couldn't create directory 'stats' in ROOT file.");
	}
	stats_dir->cd();
	RootArrayD prime_arr(NUM_EVENT_TYPES + 1);
	RootArrayD weight_moms_arr(4 * (NUM_EVENT_TYPES + 1));
	RootArrayD weight_max_arr(NUM_EVENT_TYPES + 1);
	RootArrayD num_events_arr(NUM_EVENT_TYPES + 1);
	RootArrayD num_events_acc_arr(NUM_EVENT_TYPES + 1);
	RootArrayD norm_arr(NUM_EVENT_TYPES + 1);
	Stats stats_total;
	Stats stats_acc_total;
	Double prime_total = 0.;
	// TODO: We treat counts as floating point numbers for the statistics, check
	// if this is valid later.
	Double count_total = 0;
	Double count_acc_total = 0;
	flags = std::cout.flags();
	std::cout << std::scientific << std::setprecision(OUTPUT_STATS_PRECISION);
	for (Generator& gen : gens) {
		// Need these ahead of time, for computing the weights and norms.
		prime_total += gen.prime();
		count_total += gen.count();
		count_acc_total += gen.count_acc();
	}
	for (Generator& gen : gens) {
		Stats stats = gen.weights();
		Stats stats_acc = gen.weights_acc();
		Double prime = gen.prime();
		Double mean_err;
		Double mean = stats.est_mean(&mean_err);
		Double xs = prime * mean;
		Double xs_err = prime * mean_err;
		Double rel_var_err;
		Double rel_var = stats_acc.est_rel_var(&rel_var_err);
		Double acceptance = gen.acceptance();

		// The weight is chosen so that the total cross-section averages out
		// correctly. Normally it is near one for all types of events.
		Double weight = (prime / gen.count()) / (prime_total / count_total);
		// This comes from `norm = weight * (prime_total / count_total)`.
		Double norm = prime / gen.count();
		stats_total += weight * stats;
		stats_acc_total += weight * stats_acc;

		Int_t ev_idx = static_cast<Int_t>(gen.ev_type());
		prime_arr.SetAt(prime, ev_idx);
		weight_moms_arr.SetAt(stats.ratio_m1_to_max1(), 4 * ev_idx + 0);
		weight_moms_arr.SetAt(stats.ratio_m2_to_max2(), 4 * ev_idx + 1);
		weight_moms_arr.SetAt(stats.ratio_m3_to_max3(), 4 * ev_idx + 2);
		weight_moms_arr.SetAt(stats.ratio_m4_to_max4(), 4 * ev_idx + 3);
		weight_max_arr.SetAt(stats.max1(), ev_idx);
		num_events_arr.SetAt(gen.count(), ev_idx);
		num_events_acc_arr.SetAt(gen.count_acc(), ev_idx);
		norm_arr.SetAt(norm, ev_idx);

		std::cout << "\t" << event_type_name(gen.ev_type()) << " events:" << std::endl;
		std::cout << "\t\tweight:        " << weight << std::endl;
		std::cout << "\t\tcount:         " << gen.count_acc() << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tprime:         " << prime << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\tacceptance:    " << acceptance << std::endl;
		std::cout << "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "\t\tefficiency:    " << 1. / std::sqrt(1. + rel_var) << std::endl;
	}

	Double norm = prime_total / count_total;
	Double mean_err;
	Double mean = stats_total.est_mean(&mean_err);
	Double xs = prime_total * mean;
	Double xs_err = prime_total * mean_err;
	Double rel_var_err;
	Double rel_var = stats_acc_total.est_rel_var(&rel_var_err);
	Double acceptance = static_cast<Double>(count_acc_total) / count_total;

	prime_arr.SetAt(prime_total, NUM_EVENT_TYPES);
	weight_moms_arr.SetAt(stats_total.ratio_m1_to_max1(), 4 * NUM_EVENT_TYPES + 0);
	weight_moms_arr.SetAt(stats_total.ratio_m2_to_max2(), 4 * NUM_EVENT_TYPES + 1);
	weight_moms_arr.SetAt(stats_total.ratio_m3_to_max3(), 4 * NUM_EVENT_TYPES + 2);
	weight_moms_arr.SetAt(stats_total.ratio_m4_to_max4(), 4 * NUM_EVENT_TYPES + 3);
	weight_max_arr.SetAt(stats_total.max1(), NUM_EVENT_TYPES);
	num_events_arr.SetAt(count_total, NUM_EVENT_TYPES);
	num_events_acc_arr.SetAt(count_acc_total, NUM_EVENT_TYPES);
	norm_arr.SetAt(norm, NUM_EVENT_TYPES);

	std::cout << "\ttotal:" << std::endl;
		std::cout << "\t\tcount:         " << count_acc_total << std::endl;
		std::cout << "\t\tcross-section: " << xs << " ± " << xs_err << std::endl;
		std::cout << "\t\tprime:         " << prime_total << std::endl;
		std::cout << "\t\tnorm:          " << norm << std::endl;
		std::cout << "\t\tacceptance:    " << acceptance << std::endl;
		std::cout << "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl;
		std::cout << "\t\tefficiency:    " << 1. / std::sqrt(1. + rel_var) << std::endl;

	std::cout.flags(flags);
	stats_dir->WriteObject(&prime_arr, "prime");
	stats_dir->WriteObject(&weight_moms_arr, "weight_mom");
	stats_dir->WriteObject(&weight_max_arr, "weight_max");
	stats_dir->WriteObject(&num_events_arr, "num_events");
	stats_dir->WriteObject(&num_events_acc_arr, "num_events_acc");
	stats_dir->WriteObject(&norm_arr, "norm");

	return SUCCESS;
}

int command_merge_soft(
		std::string file_out_name,
		std::vector<std::string> file_names) {
	std::cout << "Merging parameters from files." << std::endl;
	Params params_out;
	bool first = true;
	for (std::string const& file_name : file_names) {
		TFile file(file_name.c_str(), "OPEN");
		if (file.IsZombie()) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				"File '" + file_name + "' not found.");
		}
		Params params = PARAMS_STD_FORMAT;
		params.read_root(file);
		if (first) {
			params_out = params;
			first = false;
		} else {
			params_out = merge_params(params, params_out);
		}
	}

	std::cout << "Merging statistics from files." << std::endl;
	Double primes[NUM_EVENT_TYPES + 1] = {};
	Stats stats_total[NUM_EVENT_TYPES + 1] = {};
	Double count_total[NUM_EVENT_TYPES + 1] = {};
	Double count_acc_total[NUM_EVENT_TYPES + 1] = {};
	first = true;
	for (std::string const& file_name : file_names) {
		// TODO: Ensure that the merged statistics come from the same underlying
		// FOAM.
		TFile file(file_name.c_str(), "OPEN");
		if (file.IsZombie()) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				"File '" + file_name + "' not found.");
		}
		RootArrayD* prime_arr = file.Get<RootArrayD>("stats/prime");
		RootArrayD* weight_moms_arr = file.Get<RootArrayD>("stats/weight_mom");
		RootArrayD* weight_max_arr = file.Get<RootArrayD>("stats/weight_max");
		RootArrayD* num_events_arr = file.Get<RootArrayD>("stats/num_events");
		RootArrayD* num_events_acc_arr = file.Get<RootArrayD>("stats/num_events_acc");
		RootArrayD* norm_arr = file.Get<RootArrayD>("stats/norm");
		if (
				prime_arr == nullptr
				|| prime_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| weight_moms_arr == nullptr
				|| weight_moms_arr->GetSize() != 4 * (NUM_EVENT_TYPES + 1)
				|| weight_max_arr == nullptr
				|| weight_max_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| num_events_arr == nullptr
				|| num_events_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| num_events_acc_arr == nullptr
				|| num_events_acc_arr->GetSize() != NUM_EVENT_TYPES + 1
				|| norm_arr == nullptr
				|| norm_arr->GetSize() != NUM_EVENT_TYPES + 1) {
			throw Exception(
				ERROR_FILE_NOT_FOUND,
				"Couldn't find statistics in file '" + file_name + "'.");
		}
		for (std::size_t ev_idx = 0; ev_idx < NUM_EVENT_TYPES + 1; ++ev_idx) {
			Double prime = prime_arr->At(ev_idx);
			auto count = num_events_arr->At(ev_idx);
			auto count_acc = num_events_acc_arr->At(ev_idx);
			std::array<Double, 4> moms = {
				weight_moms_arr->At(4 * ev_idx + 0),
				weight_moms_arr->At(4 * ev_idx + 1),
				weight_moms_arr->At(4 * ev_idx + 2),
				weight_moms_arr->At(4 * ev_idx + 3),
			};
			Double max = weight_max_arr->At(ev_idx);
			Stats stats(moms, max, count);
			stats_total[ev_idx] += stats;
			count_total[ev_idx] += count;
			count_acc_total[ev_idx] += count_acc;
			if (first) {
				primes[ev_idx] = prime;
			} else {
				if (primes[ev_idx] != prime) {
					throw Exception(
						ERROR_FOAM_INCOMPATIBLE,
						"FOAM from file '" + file_name + "' has incompatible "
						"prime.");
				}
			}
		}
		first = false;
	}

	TFile file_out(file_out_name.c_str(), "CREATE");
	if (file_out.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Couldn't create file '" + file_out_name + "'.");
	}
	params_out.write_root(file_out);
	TDirectory* stats_dir = file_out.mkdir("stats");
	if (stats_dir == nullptr) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Couldn't create directory 'stats' in ROOT file.");
	}
	stats_dir->cd();
	RootArrayD prime_arr_out(NUM_EVENT_TYPES + 1);
	RootArrayD weight_moms_arr_out(4 * (NUM_EVENT_TYPES + 1));
	RootArrayD weight_max_arr_out(NUM_EVENT_TYPES + 1);
	RootArrayD num_events_arr_out(NUM_EVENT_TYPES + 1);
	RootArrayD num_events_acc_arr_out(NUM_EVENT_TYPES + 1);
	RootArrayD norm_arr_out(NUM_EVENT_TYPES + 1);
	for (std::size_t ev_idx = 0; ev_idx < NUM_EVENT_TYPES + 1; ++ev_idx) {
		Double norm = primes[ev_idx] / stats_total[ev_idx].count();
		prime_arr_out.SetAt(primes[ev_idx], ev_idx);
		weight_moms_arr_out.SetAt(stats_total[ev_idx].ratio_m1_to_max1(), 4 * ev_idx + 0);
		weight_moms_arr_out.SetAt(stats_total[ev_idx].ratio_m2_to_max2(), 4 * ev_idx + 1);
		weight_moms_arr_out.SetAt(stats_total[ev_idx].ratio_m3_to_max3(), 4 * ev_idx + 2);
		weight_moms_arr_out.SetAt(stats_total[ev_idx].ratio_m4_to_max4(), 4 * ev_idx + 3);
		weight_max_arr_out.SetAt(stats_total[ev_idx].max1(), ev_idx);
		num_events_arr_out.SetAt(count_total[ev_idx], ev_idx);
		num_events_acc_arr_out.SetAt(count_acc_total[ev_idx], ev_idx);
		norm_arr_out.SetAt(norm, ev_idx);
	}
	stats_dir->WriteObject(&prime_arr_out, "prime");
	stats_dir->WriteObject(&weight_moms_arr_out, "weight_mom");
	stats_dir->WriteObject(&weight_max_arr_out, "weight_max");
	stats_dir->WriteObject(&num_events_arr_out, "num_events");
	stats_dir->WriteObject(&num_events_acc_arr_out, "num_events_acc");
	stats_dir->WriteObject(&norm_arr_out, "norm");

	std::cout << "Merging events from files." << std::endl;
	file_out.cd();
	TChain chain("events", "events");
	for (std::string const& file_name : file_names) {
		std::cout << "\t" << file_name << std::endl;
		chain.Add(file_name.c_str());
	}
	file_out.WriteObject(&chain, chain.GetName());

	return SUCCESS;
}

int command_merge_hard(
		std::string file_out_name,
		std::vector<std::string> file_names) {
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
			if (argc >= 3) {
				std::string topic = argv[2];
				if (topic == "params") {
					return command_help_params();
				} else {
					return command_help_param(topic);
				}
			} else {
				return command_help();
			}
		} else if (command == "version" || command == "-v") {
			return command_version();
		} else if (command == "inspect" || command == "-i") {
			if (argc > 3) {
				throw Exception(
					ERROR_ARG_PARSE,
					"Unexpected argument '" + std::string(argv[3]) + "'.");
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
					"Unexpected argument '" + std::string(argv[3]) + "'.");
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
					"Unexpected argument '" + std::string(argv[3]) + "'.");
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
			std::vector<std::string> files;
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
			std::vector<std::string> files;
			for (int idx = 3; idx < argc; ++idx) {
				files.push_back(argv[idx]);
			}
			return command_merge_hard(argv[2], files);
		} else {
			throw Exception(
				ERROR_ARG_PARSE,
				"Unrecognized command '" + command + "'.");
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

