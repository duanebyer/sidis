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

#include <TArrayC.h>
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

// Converts between the `sidis` 4-vector type and the ROOT 4-vector type.
TLorentzVector convert_vec4(math::Vec4 vec) {
	return TLorentzVector(vec.x, vec.y, vec.z, vec.t);
}

// Convenience structure for array of integrators for each event type.
struct IntegratorArray final {
	Integrator map[NUM_EVENT_TYPES];
	IntegratorArray() : map{} { }
	Integrator& operator[](EventType ev_type) {
		return map[static_cast<std::size_t>(ev_type)];
	}
	Integrator const& operator[](EventType ev_type) const {
		return map[static_cast<std::size_t>(ev_type)];
	}

	// Sums the integrators for each different event type together.
	Integrator total() const {
		Double prime = 0.;
		std::size_t count_acc = 0.;
		Stats weights;
		for (Integrator const& integ : map) {
			prime += integ.prime();
			count_acc += integ.count_acc();
			weights += integ.weights();
		}
		return Integrator(1. / prime, count_acc, weights);
	}

	IntegratorArray& operator+=(IntegratorArray const& rhs) {
		for (std::size_t ev_idx = 0; ev_idx < NUM_EVENT_TYPES; ++ev_idx) {
			map[ev_idx] += rhs.map[ev_idx];
		}
		return *this;
	}
};

void root_write_integs(TDirectory& dir, IntegratorArray const& integs) {
	std::string file_name = dir.GetName();
	// Create directory to hold statistics info.
	TDirectory* stats_dir = dir.mkdir("stats", "stats");
	if (stats_dir == nullptr) {
		throw Exception(
			ERROR_WRITING_STATS,
			"Could not create directory 'stats' in file '" + file_name + "'.");
	}
	stats_dir->cd();

	// Define arrays.
	RootArrayD prime(NUM_EVENT_TYPES + 1);
	RootArrayD weight_mom(4 * (NUM_EVENT_TYPES + 1));
	RootArrayD weight_max(NUM_EVENT_TYPES + 1);
	RootArrayD count(NUM_EVENT_TYPES + 1);
	RootArrayD count_acc(NUM_EVENT_TYPES + 1);
	RootArrayD norm(NUM_EVENT_TYPES + 1);

	// Fill arrays.
	for (std::size_t arr_idx = 0; arr_idx < NUM_EVENT_TYPES + 1; ++arr_idx) {
		Integrator const& integ = arr_idx == 0 ? integs.total() : integs.map[arr_idx - 1];
		prime.SetAt(integ.prime(), arr_idx);
		weight_mom.SetAt(integ.weights().ratio_m1_to_max1(), 4 * arr_idx + 0);
		weight_mom.SetAt(integ.weights().ratio_m2_to_max2(), 4 * arr_idx + 1);
		weight_mom.SetAt(integ.weights().ratio_m3_to_max3(), 4 * arr_idx + 2);
		weight_mom.SetAt(integ.weights().ratio_m4_to_max4(), 4 * arr_idx + 3);
		weight_max.SetAt(integ.weights().max1(), arr_idx);
		// TODO: We treat counts as floating point numbers for the statistics,
		// check if this is valid later.
		count.SetAt(integ.count(), arr_idx);
		count_acc.SetAt(integ.count_acc(), arr_idx);
		norm.SetAt(integ.norm(), arr_idx);
	}

	// Write arrays.
	if (stats_dir->WriteObject(&prime, "prime") == 0
			|| stats_dir->WriteObject(&weight_mom, "weight_mom") == 0
			|| stats_dir->WriteObject(&weight_max, "weight_max") == 0
			|| stats_dir->WriteObject(&count, "num_events") == 0
			|| stats_dir->WriteObject(&count_acc, "num_events_acc") == 0
			|| stats_dir->WriteObject(&norm, "norm") == 0) {
		throw Exception(
			ERROR_WRITING_STATS,
			"Could not write statistics to file '" + file_name + "'.");
	}
}

IntegratorArray root_read_integs(TDirectory& dir) {
	std::string file_name = dir.GetName();
	IntegratorArray integs;
	// Read arrays.
	RootArrayD* prime = dir.Get<RootArrayD>("stats/prime");
	RootArrayD* weight_mom = dir.Get<RootArrayD>("stats/weight_mom");
	RootArrayD* weight_max = dir.Get<RootArrayD>("stats/weight_max");
	RootArrayD* count = dir.Get<RootArrayD>("stats/num_events");
	RootArrayD* count_acc = dir.Get<RootArrayD>("stats/num_events_acc");
	RootArrayD* norm = dir.Get<RootArrayD>("stats/norm");
	if (
			prime == nullptr
			|| prime->GetSize() != NUM_EVENT_TYPES + 1
			|| weight_mom == nullptr
			|| weight_mom->GetSize() != 4 * (NUM_EVENT_TYPES + 1)
			|| weight_max == nullptr
			|| weight_max->GetSize() != NUM_EVENT_TYPES + 1
			|| count == nullptr
			|| count->GetSize() != NUM_EVENT_TYPES + 1
			|| count_acc == nullptr
			|| count_acc->GetSize() != NUM_EVENT_TYPES + 1
			|| norm == nullptr
			|| norm->GetSize() != NUM_EVENT_TYPES + 1) {
		throw Exception(
			ERROR_READING_STATS,
			"Could not read statistics from file '" + file_name + "'.");
	}

	// Fill integrators.
	for (std::size_t arr_idx = 1; arr_idx < NUM_EVENT_TYPES + 1; ++arr_idx) {
		Stats weights = Stats(
			{
				weight_mom->At(4 * arr_idx + 0),
				weight_mom->At(4 * arr_idx + 1),
				weight_mom->At(4 * arr_idx + 2),
				weight_mom->At(4 * arr_idx + 3),
			},
			weight_max->At(arr_idx),
			count->At(arr_idx));
		integs.map[arr_idx - 1] = Integrator(
			1. / prime->At(arr_idx),
			count_acc->At(arr_idx),
			weights);
	}
	return integs;
}

void stream_write_integ(std::ostream& os, std::string header, Integrator const& integ) {
	Double xs_err;
	Double xs = integ.integ(&xs_err);
	Double rel_var_err;
	Double rel_var = integ.rel_var(&rel_var_err);
	Double eff = std::exp(-0.5 * std::log1p(rel_var));

	std::ios_base::fmtflags flags(os.flags());
	os << std::scientific << std::setprecision(OUTPUT_STATS_PRECISION);
	os << "\t" << header << std::endl
		<< "\t\tcount:         " << integ.count_acc() << std::endl
		<< "\t\tcross-section: " << xs << " ± " << xs_err << std::endl
		<< "\t\tprime:         " << integ.prime() << std::endl
		<< "\t\tnorm:          " << integ.norm() << std::endl
		<< std::fixed << std::setprecision(1)
		<< "\t\tacceptance:    " << 100. * integ.acc() << '%' << std::endl
		<< std::scientific << std::setprecision(OUTPUT_STATS_PRECISION)
		<< "\t\trel. variance: " << rel_var << " ± " << rel_var_err << std::endl
		<< std::fixed << std::setprecision(1)
		<< "\t\tefficiency:    " << 100. * eff << '%' << std::endl;
	os.flags(flags);
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
		part::Nucleus target = params["setup.target"].any();
		tmd_out->reset();
		sf_out->reset(new sf::set::TestSfSet(target));
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
				"Could not find structure functions in file '" + sf_set_name
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

	try {
		params.read_root(file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_READING_PARAMS,
			"Failed to read parameters from file '" + file_name + "': "
			+ e.what());
	}

	std::cout << "Parameters:" << std::endl;
	std::ios_base::fmtflags flags(std::cout.flags());
	std::cout << std::setprecision(std::numeric_limits<Double>::digits10 + 1);
	params.write_stream(std::cout);
	std::cout.flags(flags);
	std::cout << std::endl;

	std::cout << "Statistics:" << std::endl;
	try {
		IntegratorArray integs = root_read_integs(file);
		std::size_t num_event_types = 0;
		for (int idx = 0; idx < NUM_EVENT_TYPES; ++idx) {
			EventType ev_type = static_cast<EventType>(idx);
			std::string header = event_type_name(ev_type) + std::string(" events");
			Integrator const& integ = integs[ev_type];
			if (integ.count() != 0) {
				num_event_types += 1;
				stream_write_integ(std::cout, header, integ);
			}
		}
		if (num_event_types > 1) {
			stream_write_integ(std::cout, "total", integs.total());
		}
	} catch (Exception const& e) {
		if (e.error_code == ERROR_READING_STATS) {
			// FOAM files don't store statistics.
			std::cout << "Could not find any statistics." << std::endl;
		} else {
			throw e;
		}
	}

	return SUCCESS;
}

int command_initialize(std::string params_file_name) {
	// Load parameters.
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Could not open parameter file '" + params_file_name + "'.");
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

	// Create generator file.
	std::string file_name = params["file.gen"].any();
	std::cout << "Creating generator file '" << file_name << "'." << std::endl;
	TFile file(file_name.c_str(), "CREATE");
	if (file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Could not create generator file '" + file_name + "'.");
	}

	// Extract parameters needed for building the Monte-Carlo generators.
	std::vector<EventType> ev_types = p_enabled_event_types(params);
	if (ev_types.empty()) {
		throw Exception(
			ERROR_NO_EVENT_TYPES_ENABLED,
			"No event types enabled in parameter file '" + params_file_name
			+ "'.");
	}
	struct BuilderTuple {
		Density density;
		DistParams dist_params;
	};
	std::vector<BuilderTuple> builders;
	std::cout << "Checking parameters." << std::endl;
	for (EventType ev_type : ev_types) {
		try {
			builders.emplace_back(BuilderTuple {
				Density(ev_type, params, *sf),
				DistParams(ev_type, params)
			});
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_PARAMS_INVALID,
				"Invalid parameter file '" + params_file_name + "': "
				+ e.what());
		}
	}

	// Check that all provided parameters were used.
	try {
		params.filter("init"_F).check_complete();
	} catch (std::exception const& e) {
		if (params["strict"].any()) {
			throw Exception(
				ERROR_PARAMS_INVALID,
				"Invalid parameter file '" + params_file_name + "': "
				+ e.what());
		} else {
			std::cout << "Warning: " << e.what() << std::endl;
		}
	}

	// Generate UID for each generator type.
	std::random_device rnd_dev;
	std::mt19937_64 rnd(rnd_dev());
	// Long is guaranteed to be able to store 63 bits at least.
	std::uniform_int_distribution<Long> uid_dist(
		//-0x7FFFFFFFFFFFFFFF,
		0x00000000000000000,
		0x7FFFFFFFFFFFFFFF);

	// Build the generators and write to file.
	for (BuilderTuple const& builder : builders) {
		EventType ev_type = builder.density.event_type;
		std::string ev_name = event_type_name(ev_type);
		std::string ev_key = event_type_short_name(ev_type);
		params.set(p_name_init_uid(ev_type), new ValueLong(uid_dist(rnd)));
		std::cout << "Building " << ev_name << " generator." << std::endl;
		Generator gen(builder.density);
		try {
			gen.build_dist(builder.dist_params);
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_BUILDING_FOAM,
				"Error while building " + ev_name + " generator: " + e.what());
		}
		// TODO: Write more detailed generator statistics to output. Right now,
		// only the prime is displayed. But, we could also show the efficiency,
		// number of cells, relative variance, etc.
		std::ios_base::fmtflags flags(std::cout.flags());
		std::cout << std::scientific << std::setprecision(OUTPUT_STATS_PRECISION);
		std::cout << "Finished " << ev_name << " generator with prime " << gen.prime() << "." << std::endl;
		std::cout.flags(flags);
		std::cout << "Writing " << ev_name << " generator to file." << std::endl;
		// Serialize the generator. This is a little convoluted:
		// * stringstream -> unique_ptr<char[]> -> TArrayC -> ROOT file
		try {
			std::stringstream ss;
			if (!Generator::write_dist(ss, gen)) {
				throw std::runtime_error("Could not write to stream.");
			}
			ss.seekg(0, std::ios_base::end);
			std::streamsize data_len = ss.tellg();
			ss.seekg(0, std::ios_base::beg);
			if (!ss || data_len > std::numeric_limits<Int_t>::max()) {
				throw std::runtime_error("Could not get size of stream.");
			}
			std::unique_ptr<char[]> data_ptr = std::unique_ptr<char[]>(new char[data_len]);
			ss.read(&data_ptr[0], data_len);
			if (!ss || ss.gcount() != data_len) {
				throw std::runtime_error("Could not copy stream to buffer.");
			}
			TArrayC data;
			data.Adopt(data_len, data_ptr.release());
			if (file.WriteObject(&data, ev_key.c_str()) == 0) {
				throw std::runtime_error("Could not write to file.");
			}
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_WRITING_FOAM,
				"Failed to write " + ev_name + " generator to file '"
				+ file_name + "': " + e.what());
		}
	}

	// Write parameters.
	try {
		params.write_root(file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_WRITING_PARAMS,
			"Failed to write parameters to generator file '" + file_name
			+ "'." + e.what());
	}

	std::cout << "Finished!" << std::endl;
	return SUCCESS;
}

int command_generate(std::string params_file_name) {
	// Load parameters.
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Could not open parameter file '" + params_file_name + "'.");
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

	// Open the generator file.
	std::vector<EventType> ev_types = p_enabled_event_types(params);
	if (ev_types.empty()) {
		throw Exception(
			ERROR_NO_EVENT_TYPES_ENABLED,
			"No event types enabled in parameter file '" + params_file_name
			+ "'.");
	}
	std::string foam_file_name = params["file.gen"].any();
	std::cout << "Opening generator file '" << foam_file_name << "'." << std::endl;
	TFile foam_file(foam_file_name.c_str(), "OPEN");
	if (foam_file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			"Could not find generator file '" + foam_file_name + "'.");
	}

	// Load the parameters from the generator file.
	Params params_foam = PARAMS_STD_FORMAT;
	try {
		params_foam.read_root(foam_file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_READING_PARAMS,
			"Failed to read parameters from generator file '" + foam_file_name
			+ "': " + e.what());
	}

	// Create random number engine.
	std::random_device rnd_dev;
	if (!params.is_set("mc.seed")) {
		params.set("mc.seed", new ValueSeedGen(rnd_dev()));
	}
	SeedGen seed_gen = params["mc.seed"].any();
	if (seed_gen.seeds.size() != 1) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			"Invalid parameter file '" + params_file_name + "': Parameter "
			+ "'mc.seed' must provide exactly one seed.");
	}
	RndEngine rnd(*seed_gen.seeds.begin());

	// Check whether the generator is able to provide the events according to
	// what the user requested.
	try {
		check_can_provide_foam(params_foam, params);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_FOAM_INCOMPATIBLE,
			"Generator from '" + foam_file_name + "' is unable to provide "
			+ "events compatible with parameters from '" + params_file_name
			+ "': " + e.what());
	}

	// Keeps all of the relevant data for each generator together.
	struct GenTuple {
		Generator gen;
		IntegratorAccum integ;
		Double rej_scale;
	};
	std::vector<GenTuple> gens;

	// Deserialize the generators.
	for (EventType ev_type : ev_types) {
		std::string ev_name = event_type_name(ev_type);
		std::string ev_key = event_type_short_name(ev_type);
		// Record the UID from the generator.
		params.set_from(params_foam.filter(Filter(ev_key) & "init"_F & "uid"_F));
		std::cout << "Loading " << ev_name << " generator from file." << std::endl;
		try {
			TArrayC* data = foam_file.Get<TArrayC>(ev_key.c_str());
			if (data == nullptr) {
				throw std::runtime_error("Could not find key '" + ev_key + "'.");
			}
			std::stringstream ss;
			ss.write(data->GetArray(), data->GetSize());
			if (!ss) {
				throw std::runtime_error("Could not copy to buffer.");
			}
			Generator gen(Density(ev_type, params, *sf));
			if (!Generator::read_dist(ss, gen)) {
				throw std::runtime_error("Could not read from buffer.");
			}
			Double rej_scale = params[p_name_gen_rej_scale(ev_type)].any();
			IntegratorAccum integ(1. / gen.prime());
			gens.emplace_back(GenTuple {
				Generator(std::move(gen)),
				integ,
				rej_scale,
			});
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_READING_FOAM,
				"Failed to read " + ev_name + " generator from file '"
				+ foam_file_name + "': " + e.what());
		}
	}

	foam_file.Close();

	// Open the event file.
	std::string event_file_name = params["file.event"].any();
	std::cout << "Opening event file '" << event_file_name << "'." << std::endl;
	TFile event_file(event_file_name.c_str(), "CREATE");
	if (event_file.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Could not create event file '" + event_file_name + "'.");
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
	Long num_events = params["mc.num_events"].any();
	if (num_events <= 0) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			"Parameter 'mc.num_events' must be greater than zero.");
	}

	// Prepare branches in output ROOT file.
	event_file.cd();
	TTree events("events", "events");
	Int ev_idx;
	Double weight;
	Double x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R;
	TLorentzVector p, k1, q, k2, ph, k;
	events.Branch("type", &ev_idx);
	events.Branch("weight", &weight);
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
	sf::SfLP sf_out;
	if (write_sf_set) {
		// TODO: Right now, this is depending on the `SfXX` structure having a
		// very specific format. This isn't guaranteed to be true in the future,
		// if things get reorganized. Not sure what a better approach is, as
		// ROOT doesn't have another way of writing plain C-structs into trees.
		events.Branch("sf", &sf_out,
			"F_UUL/D:F_UUT/D:F_UU_cos_phih/D:F_UU_cos_2phih/D:F_UL_sin_phih/D:F_UL_sin_2phih/D:F_UTL_sin_phih_m_phis/D:F_UTT_sin_phih_m_phis/D:F_UT_sin_2phih_m_phis/D:F_UT_sin_3phih_m_phis/D:F_UT_sin_phis/D:F_UT_sin_phih_p_phis/D:F_LU_sin_phih/D:F_LL/D:F_LL_cos_phih/D:F_LT_cos_phih_m_phis/D:F_LT_cos_2phih_m_phis/D:F_LT_cos_phis/D");
	}
	Double mc_coords[9];
	Double jacobian;
	if (write_mc_coords) {
		events.Branch("mc_coords", &mc_coords, "mc_coords[9]/D");
		events.Branch("jac", &jacobian);
		throw Exception(
			ERROR_UNIMPLEMENTED,
			"Parameter 'file.write_mc_coords' not yet implemented.");
	}

	// Check that all provided parameters were used.
	try {
		params.filter("gen"_F).check_complete();
	} catch (std::exception const& e) {
		if (params["strict"].any()) {
			throw Exception(
				ERROR_PARAMS_INVALID,
				"Invalid parameter file '" + params_file_name + "': "
				+ e.what());
		} else {
			std::cout << "Warning: " << e.what() << std::endl;
		}
	}

	// Write parameters.
	try {
		params.write_root(event_file);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_WRITING_PARAMS,
			"Failed to write parameters to event file '" + event_file_name
			+ "': " + e.what());
	}

	// Generate events.
	std::cout << "Generating events." << std::endl;
	bool update_progress = true;
	std::size_t percent = 0;
	Long next_percent_rem = num_events % 100;
	Long next_percent = num_events / 100;
	for (Long event_idx = 0; event_idx < num_events; ++event_idx) {
		// Update the progress bar.
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
		std::size_t chosen_gen_idx;
		if (event_idx == 0) {
			// On the first event, we don't know anything about the total cross-
			// sections, so choose the event type arbitrarily.
			chosen_gen_idx = 0;
		} else {
			// We want to generate events so that the total weights contributed
			// by events of each type have the same ratio as the cross-sections
			// of each type.
			Double ratio_max = 0.;
			for (std::size_t idx = 0; idx < gens.size(); ++idx) {
				Integrator integ = gens[idx].integ.total_fast();
				// Equivalent to cross-section divided by total weight.
				Double ratio = integ.weights_acc().est_sqrt_m2() / integ.norm();
				if (!std::isfinite(ratio)) {
					ratio = std::numeric_limits<Double>::infinity();
				}
				if (ratio > ratio_max) {
					ratio_max = ratio;
					chosen_gen_idx = idx;
				}
			}
		}
		Generator const& gen = gens[chosen_gen_idx].gen;
		EventType ev_type = gen.event_type();

		// Generate an event.
		Event event;
		do {
			// Draw from the current generator.
			event = gen.draw(rnd);
			IntegratorAccum& integ = gens[chosen_gen_idx].integ;
			Double rej_scale = gens[chosen_gen_idx].rej_scale;

			// Apply rejection sampling through reweighting events, where
			// "rejected" events are simply reweighted to zero. The scaling used
			// here assures that the average weight will remain unchanged after
			// rescaling.
			if (rej_scale != 0.) {
				std::uniform_real_distribution<Double> dist;
				Double rej = dist(rnd);
				if (event.weight < rej * rej_scale) {
					// Event is rejected.
					event.weight = 0.;
				} else if (event.weight < rej_scale) {
					// Event is accepted.
					event.weight = rej_scale;
				} else {
					// Event overflows rejection scale. Leave it as is.
				}
			}

			// Update the integrator.
			integ += event.weight;
		} while (event.weight == 0.);

		// Fill in the branches.
		ev_idx = static_cast<Int>(ev_type) + 1;
		weight = event.weight;
		switch (event.event_type) {
		case EventType::NRAD:
			// Non-radiative event.
			{
				x = event.kin.nrad.x;
				y = event.kin.nrad.y;
				z = event.kin.nrad.z;
				ph_t_sq = event.kin.nrad.ph_t_sq;
				phi_h = event.kin.nrad.phi_h;
				phi = event.kin.nrad.phi;
				tau = 0.;
				phi_k = 0.;
				R = 0.;
				if (write_momenta) {
					kin::Final fin(init, target_pol, event.kin.nrad);
					p = convert_vec4(init.p);
					k1 = convert_vec4(init.k1);
					q = convert_vec4(fin.q);
					k2 = convert_vec4(fin.k2);
					ph = convert_vec4(fin.ph);
					k = TLorentzVector();
				}
				if (write_sf_set) {
					sf_out = sf->sf_lp(hadron, x, z, S * x * y, ph_t_sq);
				}
			}
			break;
		case EventType::RAD:
			// Radiative event.
			{
				x = event.kin.rad.x;
				y = event.kin.rad.y;
				z = event.kin.rad.z;
				ph_t_sq = event.kin.rad.ph_t_sq;
				phi_h = event.kin.rad.phi_h;
				phi = event.kin.rad.phi;
				tau = event.kin.rad.tau;
				phi_k = event.kin.rad.phi_k;
				R = event.kin.rad.R;
				if (write_momenta) {
					kin::FinalRad fin(init, target_pol, event.kin.rad);
					p = convert_vec4(init.p);
					k1 = convert_vec4(init.k1);
					q = convert_vec4(fin.q);
					k2 = convert_vec4(fin.k2);
					ph = convert_vec4(fin.ph);
					k = convert_vec4(fin.k);
				}
				if (write_sf_set) {
					sf_out = sf->sf_lp(hadron, x, z, S * x * y, ph_t_sq);
				}
			}
			break;
		default:
			UNREACHABLE();
		}
		events.Fill();
	}
	write_progress_bar(std::cout, 100);
	std::cout << std::endl;

	// Write events to file.
	std::cout << "Writing events to file." << std::endl;
	if (event_file.WriteObject(&events, events.GetName()) == 0) {
		throw Exception(
			ERROR_WRITING_EVENTS,
			"Could not write events to event file '" + event_file_name + "'.");
	}

	// Handle statistics.
	std::cout << "Statistics:" << std::endl;
	IntegratorArray integs;
	for (std::size_t idx = 0; idx < gens.size(); ++idx) {
		Generator const& gen = gens[idx].gen;
		Integrator integ = gens[idx].integ.total();
		std::string header = event_type_name(gen.event_type()) + std::string(" events");
		// Show statistics to user.
		stream_write_integ(std::cout, header, integ);
		// Add integrator to array to keep track of totals.
		integs[gen.event_type()] = integ;
	}
	if (gens.size() > 1) {
		stream_write_integ(std::cout, "total", integs.total());
	}
	// Write stats to file.
	root_write_integs(event_file, integs);

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
		try {
			params.read_root(file);
		} catch (std::exception const& e) {
			throw Exception(
				ERROR_READING_PARAMS,
				"Failed to read parameters from event file '" + file_name
				+ "': " + e.what());
		}
		if (first) {
			params_out = params;
			first = false;
		} else {
			try {
				params_out = merge_params(params, params_out);
			} catch (std::exception const& e) {
				throw Exception(
					ERROR_MERGING_PARAMS,
					"Failed to merge parameters from event file '" + file_name
					+ "': " + e.what());
			}
		}
	}

	std::cout << "Merging statistics from files." << std::endl;
	IntegratorArray integs;
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
		IntegratorArray integs_next = root_read_integs(file);
		if (first) {
			integs = integs_next;
		} else {
			integs += integs_next;
		}
		first = false;
	}

	TFile file_out(file_out_name.c_str(), "CREATE");
	if (file_out.IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			"Could not create file '" + file_out_name + "'.");
	}
	try {
		params_out.write_root(file_out);
	} catch (std::exception const& e) {
		throw Exception(
			ERROR_WRITING_PARAMS,
			"Failed to write parameters to file '" + file_out_name + "': "
			+ e.what());
	}
	root_write_integs(file_out, integs);

	std::cout << "Merging events from files." << std::endl;
	file_out.cd();
	TChain chain("events", "events");
	for (std::string const& file_name : file_names) {
		std::cout << "\t" << file_name << std::endl;
		chain.Add(file_name.c_str());
	}
	if (file_out.WriteObject(&chain, chain.GetName()) == 0) {
		throw Exception(
			ERROR_WRITING_EVENTS,
			"Could not write event chain to file '" + file_out_name + "'.");
	}

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

