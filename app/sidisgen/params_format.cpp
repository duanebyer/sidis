#include "params_format.hpp"

#include <algorithm>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace sidis;

// Within a major version, there is forward compatibility (e.x. 1.4 is
// forward-compatible with 1.5). Between major versions, there is no
// compatibility.
#define PARAMS_VERSION_MAJOR 5
#define PARAMS_VERSION_MINOR 0

/*
Tag overview.
* "meta": Meta-level information (like version) that is needed to ensure
  compatibility, but doesn't affect results.
* "file": File names.
* "write": Describes the event file format.
* "init": Used during the initialization phase, when building the generator.
* "gen": Used in the generation phase.
* "xs": Modifies the differential cross-section, excluding cuts.
* "cut": Cuts on kinematic variables.
* "dist": Modifies the distribution of weights in generated events. Note that if
  a parameter has the "xs" or "cut" tags, it should probably also have this tag.
* "nrad": Applies to non-radiative events.
* "rad": Applies to radiative events.
* "excl": Applies to exclusive events.
* "seed": For random-number-generation seeds.
* "uid": For unique identifiers.
* "num": Special tag exclusively for "mc.num_events" parameter, as this
  parameter must be treated specially.
* "cut-no-nrad": Special tag for cuts which are incompatible with non-radiative
* events. Most cuts on kinematics that depend on the radiated photon should have
* this tag.
*/
extern Params const PARAMS_STD_FORMAT = []() {
	Params params;
	// Parameters are added one-by-one to the `params` object.
	params.add_param(
		"version", new ValueVersion(PARAMS_VERSION_MAJOR, PARAMS_VERSION_MINOR),
		{ "meta" },
		"<major>.<minor>", "version of parameter file",
		"Forward compatibility within a major version (e.x. 1.4 -> 1.5), but "
		"not between major versions.");
	params.add_param(
		"strict", new ValueBool(true),
		{ "meta" },
		"<on/off>", "strictness level of parameter file",
		"If on, unused parameters are an error; otherwise, they are a warning. "
		"Default 'on'.");
	params.add_param(
		"file.event", TypeString::INSTANCE,
		{ "gen", "file", "nrad", "rad", "excl" },
		"<file>", "ROOT file for generated events",
		"Path to ROOT file to be created to hold generated events. Will give "
		"error instead of overwriting an existing file.");
	params.add_param(
		"file.gen", TypeString::INSTANCE,
		{ "init", "gen", "file", "nrad", "rad", "excl" },
		"<file>", "ROOT file for generator",
		"Path to ROOT file to save the generator to after initialization. Will "
		"give error instead of overwriting an existing file.");
	params.add_param(
		"file.write_momenta", new ValueBool(false),
		{ "gen", "write", "nrad", "rad", "excl" },
		"<on/off>", "write particle momenta to file",
		"Should particle momenta be written to output ROOT file? Default "
		"'off'.");
	params.add_param(
		"file.write_photon", new ValueBool(true),
		{ "gen", "write", "rad", "excl" },
		"<on/off>", "write radiated photon momentum to file",
		"Should radiative momentum be written to output ROOT file? Default "
		"'on'.");
	params.add_param(
		"file.write_sf_set", new ValueBool(false),
		{ "gen", "write", "nrad", "rad", "excl" },
		"<on/off>", "write model structure functions to file",
		"Should the model structure functions evaluated for each event be "
		"written to output ROOT file? Default 'off'.");
	params.add_param(
		"file.write_mc_coords", new ValueBool(false),
		{ "gen", "write", "nrad", "rad", "excl" },
		"<on/off>", "write internal MC coordinates to file",
		"Should the internal Monte-Carlo coordinates be written to output ROOT "
		"file? Default 'off'.");
	params.add_param(
		"mc.nrad.enable", new ValueBool(true),
		{ "init", "gen", "nrad" },
		"<on/off>", "generate non-radiative events",
		"Should non-radiative events be generated? Non-radiative events are "
		"those for which no real photon is radiated (including events with "
		"photon energies below the soft threshold). Default 'on'.");
	params.add_param(
		"mc.nrad.gen.rej_scale", new ValueDouble(0.),
		{ "gen", "dist", "nrad" },
		"<real>", "rejection sampling for non-radiative events",
		"The rescaling factor used for non-radiative event weights during "
		"rejection sampling. Larger values make generation slower, but improve "
		"efficiency. Suggested between 0 and 2. Default '0'.");
	params.add_param(
		"mc.nrad.init.uid", TypeLong::INSTANCE,
		{ "init", "uid", "nrad" },
		"<uid>", "UID of the non-radiative generator",
		"Internally used UID of the generator produced during non-radiative "
		"initialization. Should not be set manually.");
	params.add_param(
		"mc.nrad.init.max_cells", new ValueInt(262144),
		{ "init", "dist", "nrad" },
		"<int>", "max number of cells in non-radiative FOAM",
		"Maximum number of cells that will be created during construction of "
		"the non-radiative FOAM. Default '262144'.");
	params.add_param(
		"mc.nrad.init.target_eff", new ValueDouble(0.95),
		{ "init", "dist", "nrad" },
		"<real in [0,1]>", "efficiency for non-radiative FOAM initialization",
		"Efficiency which the non-radiative FOAM will be constructed to "
		"achieve. Larger values may cause the initialization process to take "
		"much longer. Default value '0.95'.");
	params.add_param(
		"mc.nrad.init.scale_exp", new ValueDouble(0.50),
		{ "init", "dist", "nrad" },
		"<real>", "scaling exponent for non-radiative FOAM initialization",
		"Estimate of the scaling exponent relating non-radiative FOAM "
		"efficiency to number of cells. Accurate value allows for faster "
		"construction of FOAM. Suggested between 0 and 2. Default '0.50'.");
	params.add_param(
		"mc.rad.enable", new ValueBool(true),
		{ "init", "gen", "rad" },
		"<on/off>", "generate radiative events",
		"Should radiative events be generated? Radiative events are those for "
		"which a real photon is radiated with energy above the soft threshold. "
		"Default 'on'.");
	params.add_param(
		"mc.rad.gen.rej_scale", new ValueDouble(0.),
		{ "gen", "dist", "rad" },
		"<real>", "rejection sampling for radiative events",
		"The rescaling factor used for radiative event weights during "
		"rejection sampling. Larger values make generation slower, but improve "
		"efficiency. Suggested between 0 and 2. Default '0'.");
	params.add_param(
		"mc.rad.init.uid", TypeLong::INSTANCE,
		{ "init", "uid", "rad" },
		"<uid>", "UID of the radiative generator",
		"Internally used UID of the generator produced during radiative "
		"initialization. Should not be set manually.");
	params.add_param(
		"mc.rad.init.max_cells", new ValueInt(262144),
		{ "init", "dist", "rad" },
		"<int>", "max number of cells in radiative FOAM",
		"Maximum number of cells that will be created during construction of "
		"the radiative FOAM. Default '262144'.");
	params.add_param(
		"mc.rad.init.target_eff", new ValueDouble(0.50),
		{ "init", "dist", "rad" },
		"<real in [0,1]>", "efficiency for radiative FOAM initialization",
		"Efficiency which the radiative FOAM will be constructed to achieve. "
		"Larger values may cause the initialization process to take much "
		"longer. Default value '0.50'.");
	params.add_param(
		"mc.rad.init.scale_exp", new ValueDouble(0.18),
		{ "init", "dist", "rad" },
		"<real>", "scaling exponent for radiative FOAM initialization",
		"Estimate of the scaling exponent relating radiative FOAM efficiency "
		"to number of cells. Accurate value allows for faster construction of "
		"FOAM. Suggested between 0 and 2. Default '0.18'.");
	params.add_param(
		"mc.num_events", TypeLong::INSTANCE,
		{ "gen", "num" },
		"<int>", "number of events to generate",
		"Total number of events that should be generated. Events are randomly "
		"chosen to be a mixture of 'non-radiative' and 'radiative' events.");
	params.add_param(
		"mc.seed", new ValueSeedGen(),
		{ "gen", "seed", "nrad", "rad", "excl" },
		"<int>, any", "seed for generated events",
		"The seed that will be used for generating events. If 'any', seed is "
		"chosen randomly. Default 'any'.");
	params.add_param(
		"setup.beam_energy", TypeDouble::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"<real (GeV)>", "lepton beam energy in the target rest frame",
		"Energy (in GeV) of the lepton beam in the target rest frame.");
	params.add_param(
		"setup.beam", TypeLepton::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"e, mu, tau", "lepton beam particle id",
		"Identifier for the type of particle in the lepton beam.");
	params.add_param(
		"setup.target", TypeNucleus::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"p, n, d", "target nucleus particle id",
		"Identifier for the type of the target nucleus.");
	params.add_param(
		"setup.hadron", TypeHadron::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"pi0, pi+, pi-, K0, K+, K-", "detected hadron particle id",
		"Identifier for the type of the detected hadron.");
	params.add_param(
		"setup.beam_pol", TypeDouble::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"<real in [-1,1]>", "longitudinal polarization of beam",
		"Longitudinal polarization of the lepton beam.");
	params.add_param(
		"setup.target_pol", TypeVec3::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"<vector in unit sphere>", "polarization of target nucleus",
		"Polarization of the target nucleus in the target rest frame.");
	params.add_param(
		"phys.sf_set", TypeString::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"prokudin, <ROOT dict.>", "structure function parameterization",
		"Parameterization to use for structure functions. For energies around "
		"10 GeV, try 'prokudin'. Provide custom parameterization by specifying "
		"a ROOT dictionary in shared library (*.so) form.");
	params.add_param(
		"phys.rc_method", TypeRcMethod::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"none, approx, exact", "radiative corrections method",
		"Method to use for radiative corrections. Methods are: 'none' for Born "
		"cross-section only; 'approx' for an efficient approach that is very "
		"accurate provided the soft threshold is small enough; 'exact' for the "
		"complete NLO radiative corrections. Suggested 'approx'.");
	params.add_param(
		"phys.mass_threshold", TypeDouble::INSTANCE,
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"<real (GeV)>", "mass threshold for SIDIS process",
		"Kinematic threshold mass (in GeV) of undetected particles for the "
		"SIDIS process. For example, for pion production on proton, the "
		"threshold is the mass of proton + mass of pion.");
	params.add_param(
		"phys.soft_threshold", new ValueDouble(0.01),
		{ "init", "gen", "xs", "dist", "nrad", "rad", "excl" },
		"<real (GeV)>", "radiated photon soft threshold",
		"Radiated photon energy (in GeV) dividing non-radiative events from "
		"radiative events. Specified in reference frame `p + q - p_h = 0`. "
		"Should generally be chosen as small as possible. If too small, "
		"radiative cross-section becomes invalid. If too large, RC method "
		"'approx' is no longer accurate. Default '0.01'.");
	params.add_param(
		"cut.k_0_bar", new ValueBound(math::BOUND_POSITIVE),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV)> <max (GeV)>", "bound on photon energy",
		"Defined in frame `p + q - p_h = 0`, same as 'phys.soft_threshold'. "
		"This cut is a pre-requisite for all other cuts on the radiated "
		"photon.");
	params.add_param(
		"cut.x", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min> <max>", "bound on x", "");
	params.add_param(
		"cut.y", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min> <max>", "bound on y", "");
	params.add_param(
		"cut.z", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min> <max>", "bound on z", "");
	params.add_param(
		"cut.ph_t_sq", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV^2)> <max (GeV^2)>", "bound on p_ht^2", "");
	params.add_param(
		"cut.phi_h", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (deg.)> <max (deg.)>", "bound on φ_h", "");
	params.add_param(
		"cut.phi", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (deg.)> <max (deg.)>", "bound on φ", "");
	params.add_param(
		"cut.Q_sq", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV^2)> <max (GeV^2)>", "bound on Q^2", "");
	params.add_param(
		"cut.W_sq", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV^2)> <max (GeV^2)>", "bound on W^2", "");
	params.add_param(
		"cut.mx_sq", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV^2)> <max (GeV^2)>", "bound on m_x^2 (also known as W'^2)", "");
	params.add_param(
		"cut.t", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV^2)> <max (GeV^2)>", "bound on t", "");
	params.add_param(
		"cut.r", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min> <max>", "bound on r", "");
	params.add_param(
		"cut.qt_to_Q", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min> <max>", "bound on qt/Q", "");
	params.add_param(
		"cut.tau", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "rad", "excl", "cut-no-nrad" },
		"<min> <max>", "bound on τ", "");
	params.add_param(
		"cut.phi_k", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "rad", "excl", "cut-no-nrad" },
		"<min> <max>", "bound on φ_k", "");
	params.add_param(
		"cut.R", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "rad", "cut-no-nrad" },
		"<min (GeV^2)> <max (GeV^2)>", "bound on R", "");
	params.add_param(
		"cut.lab.mom_q", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV)> <max (GeV)>", "bound on q momentum in lab frame", "");
	params.add_param(
		"cut.lab.mom_k2", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV)> <max (GeV)>", "bound on k_2 momentum in lab frame", "");
	params.add_param(
		"cut.lab.mom_h", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (GeV)> <max (GeV)>", "bound on p_h momentum in lab frame", "");
	params.add_param(
		"cut.lab.mom_k", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "rad", "cut-no-nrad" },
		"<min (GeV)> <max (GeV)>", "bound on k momentum in lab frame", "");
	params.add_param(
		"cut.lab.theta_q", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (deg.)> <max (deg.)>", "bound on q polar angle in lab frame", "");
	params.add_param(
		"cut.lab.theta_k2", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (deg.)> <max (deg.)>", "bound on k_2 polar angle in lab frame", "");
	params.add_param(
		"cut.lab.theta_h", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "nrad", "rad", "excl" },
		"<min (deg.)> <max (deg.)>", "bound on p_h polar angle in lab frame", "");
	params.add_param(
		"cut.lab.theta_k", new ValueBound(math::BOUND_INVALID),
		{ "init", "gen", "cut", "dist", "rad", "excl", "cut-no-nrad" },
		"<min (deg.)> <max (deg.)>", "bound on k polar angle in lab frame", "");
	return params;
}();

namespace {

template<typename T, typename F>
void params_merge_value(
		Params params_1,
		Params params_2,
		std::string const& name,
		F merge_op,
		Params* params_out) {
	if (!params_1.is_set(name) && !params_2.is_set(name)) {
		// If both are default, let the result also be default.
		params_out->set(name, nullptr);
	} else {
		// Otherwise merge.
		auto val_1 = params_1.get<T>(name).val;
		auto val_2 = params_2.get<T>(name).val;
		params_out->set(name, new T(merge_op(val_1, val_2)));
	}
}

void params_merge_bool_and(
		Params const& params_1,
		Params const& params_2,
		std::string const& name,
		Params* params_out) {
	return params_merge_value<ValueBool>(
		params_1, params_2, name, std::logical_and<bool>(), params_out); 
}

void params_merge_bool_or(
		Params const& params_1,
		Params const& params_2,
		std::string const& name,
		Params* params_out) {
	return params_merge_value<ValueBool>(
		params_1, params_2, name, std::logical_or<bool>(), params_out); 
}

void params_merge_version(
		Params const& params_1,
		Params const& params_2,
		std::string const& name,
		Params* params_out) {
	return params_merge_value<ValueVersion>(
		params_1, params_2, name,
		[](Version a, Version b) {
			if (a.v_major == b.v_major) {
				return a.v_minor > b.v_minor ? a : b;
			} else {
				return a.v_major > b.v_major ? a : b;
			}
		},
		params_out);
}

void params_merge_file(
		Params const& params_1,
		Params const& params_2,
		std::string const& name,
		Params* params_out) {
	return params_merge_value<ValueString>(
		params_1, params_2, name,
		[](std::string a, std::string b) { return a == b ? a : "<unknown>"; },
		params_out);
}

void params_merge_seed_gen(
		Params const& params_1,
		Params const& params_2,
		std::string const& name,
		Params* params_out) {
	return params_merge_value<ValueSeedGen>(
		params_1, params_2, name,
		[](SeedGen a, SeedGen b) { return SeedGen(a, b); },
		params_out);
}

void params_merge_count(
		Params const& params_1,
		Params const& params_2,
		std::string const& name,
		Params* params_out) {
	return params_merge_value<ValueLong>(
		params_1, params_2, name, std::plus<std::size_t>(), params_out);
}

std::runtime_error make_incompatible_param_error(
		std::string const& name,
		Value const& val_src, Value const& val_dst) {
	return std::runtime_error(
		"Parameter '" + name + "' is incompatible between source (value "
		+ val_src.to_string() + ") and dest (value " + val_dst.to_string()
		+ ").");
}

}

std::vector<EventType> p_enabled_event_types(Params& params) {
	RcMethod rc_method = params["phys.rc_method"].any();
	if (rc_method == RcMethod::NONE) {
		if (params[p_name_enable(EventType::NRAD)].any()) {
			return { EventType::NRAD };
		} else {
			return { };
		}
	} else {
		std::vector<EventType> result;
		for (int idx = 0; idx < NUM_EVENT_TYPES; ++idx) {
			EventType ev_type = static_cast<EventType>(idx);
			if (params[p_name_enable(ev_type)].any()) {
				result.push_back(ev_type);
			}
		}
		return result;
	}
}

void check_can_provide_foam(
		Params& params_foam,
		Params& params_gen) {
	params_foam.check_format(PARAMS_STD_FORMAT);
	params_gen.check_format(PARAMS_STD_FORMAT);
	// Check that versions are compatible.
	Version version_foam = params_foam.get<ValueVersion>("version");
	Version version_gen = params_gen.get<ValueVersion>("version");
	// TODO: Check versions can be handled by sidisgen.
	if (!(version_foam.v_major == version_gen.v_major
			&& version_foam.v_minor <= version_gen.v_minor)) {
		throw make_incompatible_param_error(
			"version",
			ValueVersion(version_foam),
			ValueVersion(version_gen));
	}
	// Check that all event types needed for generation are included.
	std::vector<EventType> ev_types = p_enabled_event_types(params_gen);
	Filter filter_ev_type = Filter::REJECT;
	for (EventType ev_type : ev_types) {
		filter_ev_type |= Filter(event_type_short_name(ev_type));
		if (!params_foam[p_name_enable(ev_type)].any()) {
			throw make_incompatible_param_error(
				p_name_enable(ev_type),
				ValueBool(params_foam[p_name_enable(ev_type)].any()),
				ValueBool(params_gen[p_name_enable(ev_type)].any()));
		}
	}
	// Check that the foam parameters used in initialization are equal.
	Filter filter_dist = "init"_F & ("xs"_F | "cut"_F | "dist"_F);
	Params params_dist_foam = params_foam.filter(filter_dist & filter_ev_type);
	Params params_dist_gen = params_gen.filter(filter_dist & filter_ev_type);
	params_dist_foam.check_equivalent(params_dist_gen);
	// Check that the UIDs are compatible.
	for (EventType ev_type : ev_types) {
		std::string uid_name = p_name_init_uid(ev_type);
		Long uid_foam = params_foam[uid_name].any();
		if (params_gen.is_set(uid_name)) {
			Long uid_gen = params_gen[uid_name].any();
			if (uid_foam != uid_gen) {
				throw make_incompatible_param_error(
					uid_name,
					ValueLong(uid_foam),
					ValueLong(uid_gen));
			}
		}
	}
}

Params merge_params(Params& params_1, Params& params_2) {
	params_1.check_format(PARAMS_STD_FORMAT);
	params_2.check_format(PARAMS_STD_FORMAT);
	Params result = PARAMS_STD_FORMAT;
	// Check that versions are compatible, meaning same major version.
	Version version_1 = params_1.get<ValueVersion>("version");
	Version version_2 = params_2.get<ValueVersion>("version");
	// TODO: Check version can be handled by sidisgen.
	if (version_1.v_major != version_2.v_major) {
		throw make_incompatible_param_error(
			"version",
			ValueVersion(version_1),
			ValueVersion(version_2));
	}
	// Find subset of event types shared by both sets of parameters.
	std::vector<EventType> ev_types_1 = p_enabled_event_types(params_1);
	std::vector<EventType> ev_types_2 = p_enabled_event_types(params_2);
	std::vector<EventType> ev_types_shared;
	std::vector<EventType> ev_types;
	std::set_intersection(
		ev_types_1.begin(), ev_types_1.end(),
		ev_types_2.begin(), ev_types_2.end(),
		std::back_inserter(ev_types_shared));
	std::set_union(
		ev_types_1.begin(), ev_types_1.end(),
		ev_types_2.begin(), ev_types_2.end(),
		std::back_inserter(ev_types));

	Filter filter_ev_type_1 = Filter::REJECT;
	Filter filter_ev_type_2 = Filter::REJECT;
	Filter filter_ev_type_shared = Filter::REJECT;
	for (EventType ev_type : ev_types_1) {
		filter_ev_type_1 |= Filter(event_type_short_name(ev_type));
	}
	for (EventType ev_type : ev_types_2) {
		filter_ev_type_2 |= Filter(event_type_short_name(ev_type));
	}
	for (EventType ev_type : ev_types_shared) {
		filter_ev_type_shared |= Filter(event_type_short_name(ev_type));
	}

	// Check that all parameters affecting the event distribution are equal.
	Filter filter_dist = ("xs"_F | "cut"_F | "dist"_F);
	Params params_dist_1 = params_1.filter(filter_dist & filter_ev_type_shared);
	Params params_dist_2 = params_2.filter(filter_dist & filter_ev_type_shared);
	params_dist_1.check_equivalent(params_dist_2);
	// Check that the UIDs are compatible.
	for (EventType ev_type : ev_types_shared) {
		std::string uid_name = p_name_init_uid(ev_type);
		Long uid_1 = params_1[uid_name].any();
		Long uid_2 = params_2[uid_name].any();
		if (uid_1 != uid_2) {
			throw make_incompatible_param_error(
				uid_name,
				ValueLong(uid_1),
				ValueLong(uid_2));
		}
	}
	// Check that the generation seeds are compatible, meaning non-overlapping,
	// so that each set of events are independent.
	if (params_1.is_set("mc.seed") && params_2.is_set("mc.seed")) {
		SeedGen seed_gen_1 = params_1["mc.seed"].any();
		SeedGen seed_gen_2 = params_2["mc.seed"].any();
		for (int seed : seed_gen_1.seeds) {
			if (seed_gen_2.seeds.find(seed) != seed_gen_2.seeds.end()) {
				throw make_incompatible_param_error(
					"mc.seed",
					ValueSeedGen(seed_gen_1),
					ValueSeedGen(seed_gen_2));
			}
		}
	}
	// Merge meta information.
	params_merge_version(params_1, params_2, "version", &result);
	params_merge_bool_and(params_1, params_2, "strict", &result);
	// Merge file parameters.
	Filter filter_file = "file"_F;
	for (std::string const& name : params_1.filter(filter_file).names()) {
		params_merge_file(params_1, params_2, name, &result);
	}
	// Merge write parameters.
	params_merge_bool_and(params_1, params_2, "file.write_momenta", &result);
	params_merge_bool_and(params_1, params_2, "file.write_photon", &result);
	params_merge_bool_and(params_1, params_2, "file.write_sf_set", &result);
	params_merge_bool_and(params_1, params_2, "file.write_mc_coords", &result);
	// Merge event type enables.
	for (EventType ev_type : ev_types) {
		params_merge_bool_or(params_1, params_2, p_name_enable(ev_type), &result);
	}
	// Merge seeds.
	params_merge_seed_gen(params_1, params_2, "mc.seed", &result);
	// Merge counts.
	params_merge_count(params_1, params_2, "mc.num_events", &result);
	// All other information must be equal between `params_dist_1` and
	// `params_dist_2` from the earlier call to `equivalent()`, so no need to
	// merge, just copy value direct from `params_dist_1`.
	Filter filter_set = (filter_dist | "uid"_F | "seed"_F);
	result.set_from(params_1.filter(filter_set & filter_ev_type_1));
	result.set_from(params_2.filter(filter_set & filter_ev_type_2));
	// TODO: This check verifies that every parameter set in `params_1` or
	// `params_2` is also set in `result`, and vice versa. It shouldn't ever get
	// tripped if there aren't any bugs in the merge, so enable it for debug
	// only.
	for (std::string const& name_1 : params_1.names()) {
		if (params_1.is_set(name_1) && !result.is_set(name_1)) {
			throw std::runtime_error(
				"Parameter '" + name_1 + "' failed to merge.");
		}
	}
	for (std::string const& name_2 : params_2.names()) {
		if (params_2.is_set(name_2) && !result.is_set(name_2)) {
			throw std::runtime_error(
				"Parameter '" + name_2 + "' failed to merge.");
		}
	}
	for (std::string const& name : result.names()) {
		if (result.is_set(name)
				&& !params_1.is_set(name)
				&& !params_2.is_set(name)) {
			throw std::runtime_error(
				"Parameter '" + name + "' contaminated merge.");
		}
	}
	return result;
}

#define ESC(...) __VA_ARGS__

#define VALUE_TYPE_DEFINE_SINGLETON(RType) \
	RType const RType::INSTANCE;

#define VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(RType, Wrapped) \
	Wrapped RType::read_stream_base(std::istream& is) const { \
		Wrapped val; \
		is >> val; \
		return val; \
	} \
	void RType::write_stream_base(std::ostream& os, Wrapped const& val) const { \
		os << val; \
	}

#define VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(RType, Wrapped, enum_count, enum_vals_list, names_list) \
	Wrapped RType::read_stream_base(std::istream& is) const { \
		Wrapped enum_vals[enum_count] = enum_vals_list; \
		std::vector<char const*> names[enum_count] = names_list; \
		std::string name_in; \
		is >> name_in; \
		for (std::size_t idx = 0; idx < (enum_count); ++idx) { \
			for (char const* name : names[idx]) { \
				if (name_in == name) { \
					return enum_vals[idx]; \
				} \
			} \
		} \
		is.setstate(std::ios_base::failbit); \
		return enum_vals[0]; \
	} \
	void RType::write_stream_base(std::ostream& os, Wrapped const& val) const { \
		Wrapped enum_vals[enum_count] = enum_vals_list; \
		std::vector<char const*> names[enum_count] = names_list; \
		for (std::size_t idx = 0; idx < (enum_count); ++idx) { \
			if (enum_vals[idx] == val) { \
				os << names[idx][0]; \
				return; \
			} \
		} \
		os.setstate(std::ios_base::failbit); \
	}

#define VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(RType, Wrapped, Number) \
	TParameter<Number> RType::convert_to_root_base(Wrapped const& val) const { \
		return TParameter<Number>("", static_cast<Number>(val)); \
	} \
	Wrapped RType::convert_from_root_base(TParameter<Number>& val) const { \
		return static_cast<Wrapped>(val.GetVal()); \
	}

// Singleton definitions.
// Version.
VALUE_TYPE_DEFINE_SINGLETON(TypeVersion)
// Numbers.
VALUE_TYPE_DEFINE_SINGLETON(TypeDouble)
VALUE_TYPE_DEFINE_SINGLETON(TypeInt)
VALUE_TYPE_DEFINE_SINGLETON(TypeLong)
VALUE_TYPE_DEFINE_SINGLETON(TypeBool)
// Strings.
VALUE_TYPE_DEFINE_SINGLETON(TypeString)
// Random number seeds.
VALUE_TYPE_DEFINE_SINGLETON(TypeSeedGen)
// Enums.
VALUE_TYPE_DEFINE_SINGLETON(TypeRcMethod)
VALUE_TYPE_DEFINE_SINGLETON(TypeNucleus)
VALUE_TYPE_DEFINE_SINGLETON(TypeLepton)
VALUE_TYPE_DEFINE_SINGLETON(TypeHadron)
// Math types.
VALUE_TYPE_DEFINE_SINGLETON(TypeVec3)
VALUE_TYPE_DEFINE_SINGLETON(TypeBound)

// Stream IO definitions.
// Numbers.
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(TypeDouble, Double)
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(TypeInt, Int)
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(TypeLong, Long)
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeRcMethod, RcMethod, 3,
	ESC({ RcMethod::NONE, RcMethod::APPROX, RcMethod::EXACT }),
	ESC({ { "none" }, { "approx" }, { "exact" } }))
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeNucleus, part::Nucleus, 3,
	ESC({ part::Nucleus::P, part::Nucleus::N, part::Nucleus::D }),
	ESC({ { "p", "proton" }, { "n", "neutron" }, { "d", "deuteron" } }))
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeLepton, part::Lepton, 3,
	ESC({ part::Lepton::E, part::Lepton::MU, part::Lepton::TAU }),
	ESC({ { "e", "electron" }, { "mu", "muon" }, { "tau", "tauon" } }))
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeHadron, part::Hadron, 6,
	ESC({ part::Hadron::PI_0, part::Hadron::PI_P, part::Hadron::PI_M, part::Hadron::K_0, part::Hadron::K_P, part::Hadron::K_M }),
	ESC({ { "pi0", "pion0" }, { "pi+", "pion+" }, { "pi-", "pion-" }, { "K0", "kaon0" }, { "K+", "kaon+" }, { "K-", "kaon-" } }))

// ROOT conversion definitions.
// Numbers.
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeDouble, Double, Double)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeInt, Int, Int)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeLong, Long, Long)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeBool, bool, bool)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeRcMethod, RcMethod, int)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeNucleus, part::Nucleus, int)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeLepton, part::Lepton, int)
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeHadron, part::Hadron, int)

// Read/write to stream.

// Version.
Version TypeVersion::read_stream_base(std::istream& is) const {
	Version version;
	is >> version.v_major;
	char period;
	is >> period;
	if (period != '.') {
		is.setstate(std::ios_base::failbit);
	}
	is >> version.v_minor;
	return version;
}
void TypeVersion::write_stream_base(std::ostream& os, Version const& version) const {
	os << version.v_major << "." << version.v_minor;
}

// Numbers.
bool TypeBool::read_stream_base(std::istream& is) const {
	std::string str;
	is >> str;
	if (str == "1" || str == "on" || str == "true" || str == "yes") {
		return true;
	} else if (str == "0" || str == "off" || str == "false" || str == "no") {
		return false;
	} else {
		is.setstate(std::ios_base::failbit);
		return false;
	}
}
void TypeBool::write_stream_base(std::ostream& os, bool const& val) const {
	os << (val ? "on" : "off");
}

// Strings.
std::string TypeString::read_stream_base(std::istream& is) const {
	std::string str;
	std::getline(is, str);
	return str;
}
void TypeString::write_stream_base(std::ostream& os, std::string const& str) const {
	os << str;
}

// Random number seeds.
SeedGen TypeSeedGen::read_stream_base(std::istream& is) const {
	SeedGen seed;
	while (is && !is.eof()) {
		int val;
		is >> val;
		seed.seeds.insert(val);
	}
	return seed;
}
void TypeSeedGen::write_stream_base(std::ostream& os, SeedGen const& seed) const {
	bool first = true;
	for (int val : seed.seeds) {
		if (!first) {
			os << ' ';
			first = false;
		}
		os << val;
	}
}

// Math types.
math::Vec3 TypeVec3::read_stream_base(std::istream& is) const {
	math::Vec3 vec;
	is >> vec.x >> vec.y >> vec.z;
	return vec;
}
void TypeVec3::write_stream_base(std::ostream& os, math::Vec3 const& vec) const {
	os << vec.x << ' ' << vec.y << ' ' << vec.z;
}
math::Bound TypeBound::read_stream_base(std::istream& is) const {
	Real min, max;
	is >> min >> max;
	if (!(min <= max)) {
		is.setstate(std::ios_base::failbit);
	}
	return math::Bound(min, max);
}
void TypeBound::write_stream_base(std::ostream& os, math::Bound const& bound) const {
	os << bound.min() << ' ' << bound.max();
}

// Read/write to ROOT.

// Version.
RootArrayI TypeVersion::convert_to_root_base(Version const& version) const {
	Int vals[2] = { version.v_major, version.v_minor };
	return RootArrayI(2, vals);
}
Version TypeVersion::convert_from_root_base(RootArrayI& version) const {
	if (version.GetSize() != 2) {
		throw std::runtime_error("Wrong number of array elements.");
	}
	return Version(version.At(0), version.At(1));
}
// Strings.
TObjString TypeString::convert_to_root_base(std::string const& str) const {
	return TObjString(str.c_str());
}
std::string TypeString::convert_from_root_base(TObjString& str) const {
	return str.GetString().Data();
}
// Random number seeds.
RootArrayI TypeSeedGen::convert_to_root_base(SeedGen const& seed) const {
	std::vector<Int> vals(seed.seeds.begin(), seed.seeds.end());
	return RootArrayI(vals.size(), vals.data());
}
SeedGen TypeSeedGen::convert_from_root_base(RootArrayI& seed) const {
	SeedGen result;
	for (Int_t i = 0; i < seed.GetSize(); ++i) {
		result.seeds.insert(seed.At(i));
	}
	return result;
}
// Math types.
RootArrayD TypeVec3::convert_to_root_base(math::Vec3 const& vec) const {
	Double_t vals[3] = { vec.x, vec.y, vec.z };
	return RootArrayD(3, vals);
}
math::Vec3 TypeVec3::convert_from_root_base(RootArrayD& vec) const {
	if (vec.GetSize() != 3) {
		throw std::runtime_error("Wrong number of array elements.");
	}
	return math::Vec3(vec.At(0), vec.At(1), vec.At(2));
}
RootArrayD TypeBound::convert_to_root_base(math::Bound const& bound) const {
	Double_t vals[3] = { bound.min(), bound.max() };
	return RootArrayD(2, vals);
}
math::Bound TypeBound::convert_from_root_base(RootArrayD& bound) const {
	if (bound.GetSize() != 2 || !(bound.At(0) <= bound.At(1))) {
		throw std::runtime_error("Wrong number of array elements.");
	}
	return math::Bound(bound.At(0), bound.At(1));
}

