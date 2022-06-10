#include "params_format.hpp"

#include <ios>
#include <istream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace sidis;

// Within a major version, there is forward compatibility (e.x. 1.4 is
// forward-compatible with 1.5). Between major versions, there is no
// compatibility.
#define PARAMS_VERSION_MAJOR 4
#define PARAMS_VERSION_MINOR 1

Params params_format() {
	Params params;
	// Parameters are added one-by-one to the `params` object.
	params.add_param(
		"version", new ValueVersion(PARAMS_VERSION_MAJOR, PARAMS_VERSION_MINOR),
		{ "meta" },
		"<major>.<minor>", "version of parameter file",
		"Forward compatibility within a major version (e.x. 1.4 -> 1.5), but "
		"not between major versions.");
	params.add_param(
		"strict", new ValueBool(ProvideBool::LT_EQ, true),
		{ "meta" },
		"<on/off>", "strictness level of parameter file",
		"If on, unused parameters are an error; otherwise, they are a warning. "
		"Default 'on'.");
	params.add_param(
		"file.event", TypeString::instance(ProvideSimple::ANY),
		{ "file" },
		"<file>", "ROOT file for generated events",
		"Path to ROOT file to be created to hold generated events. Will give "
		"error instead of overwriting an existing file.");
	params.add_param(
		"file.foam", TypeString::instance(ProvideSimple::ANY),
		{ "file" },
		"<file>", "ROOT file for FOAM",
		"Path to the ROOT file to hold the FOAM during initialization. Will "
		"give error instead of overwriting an existing file.");
	params.add_param(
		"file.write_momenta", new ValueBool(false),
		{ "write" },
		"<on/off>", "write particle momenta to file",
		"Should particle momenta be written to output ROOT file? Default "
		"'off'.");
	params.add_param(
		"file.write_photon", new ValueBool(true),
		{ "write" },
		"<on/off>", "write radiated photon momentum to file",
		"Should radiative momentum be written to output ROOT file? Default "
		"'on'.");
	params.add_param(
		"file.write_sf_set", new ValueBool(false),
		{ "write" },
		"<on/off>", "write model structure functions to file",
		"Should the model structure functions evaluated for each event be "
		"written to output ROOT file? Default 'off'.");
	params.add_param(
		"file.write_mc_coords", new ValueBool(false),
		{ "write" },
		"<on/off>", "write internal MC coordinates to file",
		"Should the internal Monte-Carlo coordinates be written to output ROOT "
		"file? Default 'off'.");
	params.add_param(
		"mc.nrad.gen", new ValueBool(ProvideBool::LT_EQ, true),
		{ "nrad", "gen" },
		"<on/off>", "generate non-radiative events",
		"Should non-radiative events be generated? Non-radiative events are "
		"those for which no real photon is radiated (including events with "
		"photon energies below the soft threshold). Default 'on'.");
	params.add_param(
		"mc.nrad.gen.rej_scale", new ValueDouble(0.),
		{ "nrad", "gen" },
		"<real>", "rejection sampling for non-radiative events",
		"The rescaling factor used for non-radiative event weights during "
		"rejection sampling. Larger values make generation slower, but improve "
		"efficiency. Suggested between 0 and 2. Default '0'.");
	params.add_param(
		"mc.nrad.init.seed", new ValueSeedInit(),
		{ "nrad", "init" },
		"<int>, any", "seed for non-radiative FOAM initialization",
		"The seed that will be used for constructing the non-radiative FOAM. "
		"If 'any', seed is chosen randomly. Default 'any'.");
	params.add_param(
		"mc.nrad.init.max_cells", new ValueInt(ProvideOrder::LT_EQ, 262144),
		{ "nrad", "init" },
		"<int>", "max number of cells in non-radiative FOAM",
		"Maximum number of cells that will be created during construction of "
		"the non-radiative FOAM. Default '262144'.");
	params.add_param(
		"mc.nrad.init.target_eff", new ValueDouble(ProvideOrder::LT_EQ, 0.95),
		{ "nrad", "init" },
		"<real in [0,1]>", "efficiency for non-radiative FOAM initialization",
		"Efficiency which the non-radiative FOAM will be constructed to "
		"achieve. Larger values may cause the initialization process to take "
		"much longer. Default value '0.95'.");
	params.add_param(
		"mc.nrad.init.scale_exp", new ValueDouble(0.50),
		{ "nrad", "init" },
		"<real>", "scaling exponent for non-radiative FOAM initialization",
		"Estimate of the scaling exponent relating non-radiative FOAM "
		"efficiency to number of cells. Accurate value allows for faster "
		"construction of FOAM. Suggested between 0 and 2. Default '0.50'.");
	params.add_param(
		"mc.rad.gen", new ValueBool(ProvideBool::LT_EQ, true),
		{ "rad", "gen" },
		"<on/off>", "generate radiative events",
		"Should radiative events be generated? Radiative events are those for "
		"which a real photon is radiated with energy above the soft threshold. "
		"Default 'on'.");
	params.add_param(
		"mc.rad.gen.rej_scale", new ValueDouble(0.),
		{ "rad", "gen" },
		"<real>", "rejection sampling for radiative events",
		"The rescaling factor used for radiative event weights during "
		"rejection sampling. Larger values make generation slower, but improve "
		"efficiency. Suggested between 0 and 2. Default '0'.");
	params.add_param(
		"mc.rad.init.seed", new ValueSeedInit(),
		{ "rad", "init" },
		"<int>, any", "seed for radiative FOAM initialization",
		"The seed that will be used for constructing the radiative FOAM. If "
		"'any', seed is chosen randomly. Default 'any'.");
	params.add_param(
		"mc.rad.init.max_cells", new ValueInt(ProvideOrder::LT_EQ, 262144),
		{ "rad", "init" },
		"<int>", "max number of cells in radiative FOAM",
		"Maximum number of cells that will be created during construction of "
		"the radiative FOAM. Default '262144'.");
	params.add_param(
		"mc.rad.init.target_eff", new ValueDouble(ProvideOrder::LT_EQ, 0.50),
		{ "rad", "init" },
		"<real in [0,1]>", "efficiency for radiative FOAM initialization",
		"Efficiency which the radiative FOAM will be constructed to achieve. "
		"Larger values may cause the initialization process to take much "
		"longer. Default value '0.50'.");
	params.add_param(
		"mc.rad.init.scale_exp", new ValueDouble(0.18),
		{ "rad", "init" },
		"<real>", "scaling exponent for radiative FOAM initialization",
		"Estimate of the scaling exponent relating radiative FOAM efficiency "
		"to number of cells. Accurate value allows for faster construction of "
		"FOAM. Suggested between 0 and 2. Default '0.18'.");
	params.add_param(
		"mc.num_events", TypeInt::instance(ProvideOrder::ANY),
		{ "num" },
		"<int>", "number of events to generate",
		"Total number of events that should be generated. Events are randomly "
		"chosen to be a mixture of 'non-radiative' and 'radiative' events.");
	params.add_param(
		"mc.seed", new ValueSeedGen(),
		{ "nrad", "rad", "gen" },
		"<int>, any", "seed for generated events",
		"The seed that will be used for generating events. If 'any', seed is "
		"chosen randomly. Default 'any'.");
	params.add_param(
		"setup.beam_energy", TypeDouble::instance(),
		{ "nrad", "rad", "dist" },
		"<real (GeV)>", "lepton beam energy in the target rest frame",
		"Energy (in GeV) of the lepton beam in the target rest frame.");
	params.add_param(
		"setup.beam", TypeLepton::instance(),
		{ "nrad", "rad", "dist" },
		"e, mu, tau", "lepton beam particle id",
		"Identifier for the type of particle in the lepton beam.");
	params.add_param(
		"setup.target", TypeNucleus::instance(),
		{ "nrad", "rad", "dist" },
		"p, n, d", "target nucleus particle id",
		"Identifier for the type of the target nucleus.");
	params.add_param(
		"setup.hadron", TypeHadron::instance(),
		{ "nrad", "rad", "dist" },
		"pi0, pi+, pi-, K0, K+, K-", "detected hadron particle id",
		"Identifier for the type of the detected hadron.");
	params.add_param(
		"setup.beam_pol", TypeDouble::instance(),
		{ "nrad", "rad", "dist" },
		"<real in [-1,1]>", "longitudinal polarization of beam",
		"Longitudinal polarization of the lepton beam.");
	params.add_param(
		"setup.target_pol", TypeVec3::instance(),
		{ "nrad", "rad", "dist" },
		"<vector in unit sphere>", "polarization of target nucleus",
		"Polarization of the target nucleus in the target rest frame.");
	params.add_param(
		"phys.sf_set", TypeString::instance(),
		{ "nrad", "rad", "dist" },
		"prokudin, <ROOT dict.>", "structure function parameterization",
		"Parameterization to use for structure functions. For energies around "
		"10 GeV, try 'prokudin'. Provide custom parameterization by specifying "
		"a ROOT dictionary in shared library (*.so) form.");
	params.add_param(
		"phys.rc_method", TypeRcMethod::instance(),
		{ "nrad", "rad", "dist" },
		"none, approx, exact", "radiative corrections method",
		"Method to use for radiative corrections. Methods are: 'none' for Born "
		"cross-section only; 'approx' for an efficient approach that is very "
		"accurate provided the soft threshold is small enough; 'exact' for the "
		"complete NLO radiative corrections. Suggested 'approx'.");
	params.add_param(
		"phys.mass_threshold", TypeDouble::instance(),
		{ "nrad", "rad", "dist" },
		"<real (GeV)>", "mass threshold for SIDIS process",
		"Kinematic threshold mass (in GeV) of undetected particles for the "
		"SIDIS process. For example, for pion production on proton, the "
		"threshold is the mass of proton + mass of pion.");
	params.add_param(
		"phys.soft_threshold", new ValueDouble(0.01),
		{ "nrad", "rad", "dist" },
		"<real (GeV)>", "radiated photon soft threshold",
		"Radiated photon energy (in GeV) dividing non-radiative events from "
		"radiative events. Specified in reference frame `p + q - p_h = 0`. "
		"Should generally be chosen as small as possible. If too small, "
		"radiative cross-section becomes invalid. If too large, RC method "
		"'approx' is no longer accurate. Default '0.01'.");
	params.add_param(
		"cut.k_0_bar", new ValueBound(math::Bound::POSITIVE),
		{ "nrad", "rad", "cut" },
		"<min (GeV)> <max (GeV)>", "bound on photon energy",
		"Defined in frame `p + q - p_h = 0`, same as 'phys.soft_threshold'. "
		"This cut is a pre-requisite for all other cuts on the radiated "
		"photon.");
	params.add_param(
		"cut.x", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min> <max>", "bound on x", "");
	params.add_param(
		"cut.y", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min> <max>", "bound on y", "");
	params.add_param(
		"cut.z", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min> <max>", "bound on z", "");
	params.add_param(
		"cut.ph_t_sq", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV²)> <max (GeV²)>", "bound on p_ht²", "");
	params.add_param(
		"cut.phi_h", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (deg.)> <max (deg.)>", "bound on φ_h", "");
	params.add_param(
		"cut.phi", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (deg.)> <max (deg.)>", "bound on φ", "");
	params.add_param(
		"cut.Q_sq", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV²)> <max (GeV²)>", "bound on Q²", "");
	params.add_param(
		"cut.W_sq", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV²)> <max (GeV²)>", "bound on W²", "");
	params.add_param(
		"cut.mx_sq", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV²)> <max (GeV²)>", "bound on m_x² (also known as W'²)", "");
	params.add_param(
		"cut.t", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV²)> <max (GeV²)>", "bound on t", "");
	params.add_param(
		"cut.r", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min> <max>", "bound on r", "");
	params.add_param(
		"cut.qt_to_Q", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min> <max>", "bound on qt/Q", "");
	params.add_param(
		"cut.tau", new ValueBound(math::Bound::INVALID), { "rad", "cut" },
		"<min> <max>", "bound on τ", "");
	params.add_param(
		"cut.phi_k", new ValueBound(math::Bound::INVALID), { "rad", "cut" },
		"<min> <max>", "bound on φ_k", "");
	params.add_param(
		"cut.R", new ValueBound(math::Bound::INVALID), { "rad", "cut" },
		"<min (GeV²)> <max (GeV²)>", "bound on R", "");
	params.add_param(
		"cut.lab.mom_q", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV)> <max (GeV)>", "bound on q momentum in lab frame", "");
	params.add_param(
		"cut.lab.mom_k2", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV)> <max (GeV)>", "bound on k_2 momentum in lab frame", "");
	params.add_param(
		"cut.lab.mom_h", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (GeV)> <max (GeV)>", "bound on p_h momentum in lab frame", "");
	params.add_param(
		"cut.lab.mom_k", new ValueBound(math::Bound::INVALID), { "rad", "cut" },
		"<min (GeV)> <max (GeV)>", "bound on k momentum in lab frame", "");
	params.add_param(
		"cut.lab.theta_q", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (deg.)> <max (deg.)>", "bound on q polar angle in lab frame", "");
	params.add_param(
		"cut.lab.theta_k2", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (deg.)> <max (deg.)>", "bound on k_2 polar angle in lab frame", "");
	params.add_param(
		"cut.lab.theta_h", new ValueBound(math::Bound::INVALID), { "nrad", "rad", "cut" },
		"<min (deg.)> <max (deg.)>", "bound on p_h polar angle in lab frame", "");
	params.add_param(
		"cut.lab.theta_k", new ValueBound(math::Bound::INVALID), { "rad", "cut" },
		"<min (deg.)> <max (deg.)>", "bound on k polar angle in lab frame", "");
	return params;
}

#define ESC(...) __VA_ARGS__

#define VALUE_TYPE_DEFINE_SINGLETON(RType, Provide, provide_count) \
	RType const& RType::instance(Provide provide) { \
		static std::unique_ptr<RType> singleton[provide_count]; \
		int index = static_cast<int>(provide); \
		if (!(index >= 0 && index < provide_count)) { \
			throw std::runtime_error("Invalid " #Provide "."); \
		} \
		if (singleton[index] == nullptr) { \
			singleton[index] = std::unique_ptr<RType>(new RType(provide)); \
		} \
		return *singleton[index]; \
	}

#define VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(RType, Wrapped) \
	bool RType::can_provide_base(Wrapped const& prod, Wrapped const& cons) const { \
		switch (_provide) { \
		case ProvideSimple::EQ: \
			return prod == cons; \
		case ProvideSimple::ANY: \
			return true; \
		default: \
			throw std::runtime_error("Invalid ProvideOrder."); \
		} \
	}

#define VALUE_TYPE_DEFINE_PROVIDE_ORDER(RType, Wrapped) \
	bool RType::can_provide_base(Wrapped const& prod, Wrapped const& cons) const { \
		switch (_provide) { \
		case ProvideOrder::EQ: \
			return prod == cons; \
		case ProvideOrder::ANY: \
			return true; \
		case ProvideOrder::LT_EQ: \
			return prod >= cons; \
		case ProvideOrder::GT_EQ: \
			return prod <= cons; \
		default: \
			throw std::runtime_error("Invalid ProvideOrder."); \
		} \
	}

#define VALUE_TYPE_DEFINE_PROVIDE_BOOL(RType, Wrapped) \
	bool RType::can_provide_base(Wrapped const& prod, Wrapped const& cons) const { \
		switch (_provide) { \
		case ProvideBool::EQ: \
			return prod == cons; \
		case ProvideBool::ANY: \
			return true; \
		case ProvideBool::AND: \
			return prod && cons; \
		case ProvideBool::OR: \
			return prod || cons; \
		case ProvideBool::GT_EQ: \
			return !(!prod && cons); \
		case ProvideBool::LT_EQ: \
			return !(prod && !cons); \
		default: \
			throw std::runtime_error("Invalid ProvideBool."); \
		} \
	}

#define VALUE_TYPE_DEFINE_PROVIDE_BOUND(RType, Wrapped) \
	bool RType::can_provide_base(Wrapped const& prod, Wrapped const& cons) const { \
		switch (_provide) { \
		case ProvideBound::EQ: \
			return prod == cons; \
		case ProvideBound::ANY: \
			return true; \
		case ProvideBound::IN: \
			return prod.contains(cons); \
		case ProvideBound::OUT: \
			return cons.contains(prod); \
		default: \
			throw std::runtime_error("Invalid ProvideBound."); \
		} \
	}

#define VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(RType, Wrapped) \
	Wrapped RType::read_stream_base(std::istream& is) const { \
		Wrapped value; \
		is >> value; \
		return value; \
	} \
	void RType::write_stream_base(std::ostream& os, Wrapped const& value) const { \
		os << value; \
	}

#define VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(RType, Wrapped, enum_count, enum_vals_list, names_list) \
	Wrapped RType::read_stream_base(std::istream& is) const { \
		Wrapped enum_vals[enum_count] = enum_vals_list; \
		std::vector<char const*> names[enum_count] = names_list; \
		std::string name_in; \
		is >> name_in; \
		for (std::size_t idx = 0; idx < enum_count; ++idx) { \
			for (char const* name : names[idx]) { \
				if (name_in == name) { \
					return enum_vals[idx]; \
				} \
			} \
		} \
		is.setstate(std::ios_base::failbit); \
		return enum_vals[0]; \
	} \
	void RType::write_stream_base(std::ostream& os, Wrapped const& value) const { \
		Wrapped enum_vals[enum_count] = enum_vals_list; \
		std::vector<char const*> names[enum_count] = names_list; \
		for (std::size_t idx = 0; idx < enum_count; ++idx) { \
			if (enum_vals[idx] == value) { \
				os << names[idx][0]; \
				return; \
			} \
		} \
		os.setstate(std::ios_base::failbit); \
	}

#define VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(RType, Wrapped, Number) \
	TParameter<Number> RType::convert_to_root_base(Wrapped const& value) const { \
		return TParameter<Number>("", static_cast<Number>(value)); \
	} \
	Wrapped RType::convert_from_root_base(TParameter<Number>& value) const { \
		return static_cast<Wrapped>(value.GetVal()); \
	}

// Singleton definitions.
// Version.
VALUE_TYPE_DEFINE_SINGLETON(TypeVersion, ProvideSpecial, 3);
// Numbers.
VALUE_TYPE_DEFINE_SINGLETON(TypeDouble, ProvideOrder, 4);
VALUE_TYPE_DEFINE_SINGLETON(TypeInt, ProvideOrder, 4);
VALUE_TYPE_DEFINE_SINGLETON(TypeBool, ProvideBool, 6);
// Strings.
VALUE_TYPE_DEFINE_SINGLETON(TypeString, ProvideSimple, 2);
// Random number seeds.
VALUE_TYPE_DEFINE_SINGLETON(TypeSeedInit, ProvideSpecial, 3);
VALUE_TYPE_DEFINE_SINGLETON(TypeSeedGen, ProvideSpecial, 3);
// Enums.
VALUE_TYPE_DEFINE_SINGLETON(TypeRcMethod, ProvideSimple, 2);
VALUE_TYPE_DEFINE_SINGLETON(TypeNucleus, ProvideSimple, 2);
VALUE_TYPE_DEFINE_SINGLETON(TypeLepton, ProvideSimple, 2);
VALUE_TYPE_DEFINE_SINGLETON(TypeHadron, ProvideSimple, 2);
// Math types.
VALUE_TYPE_DEFINE_SINGLETON(TypeVec3, ProvideSimple, 2);
VALUE_TYPE_DEFINE_SINGLETON(TypeBound, ProvideBound, 4);

// Can provide definitions.
// Numbers.
VALUE_TYPE_DEFINE_PROVIDE_ORDER(TypeDouble, double);
VALUE_TYPE_DEFINE_PROVIDE_ORDER(TypeInt, int);
VALUE_TYPE_DEFINE_PROVIDE_BOOL(TypeBool, bool);
// Strings.
VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(TypeString, std::string);
// Enums.
VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(TypeRcMethod, RcMethod);
VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(TypeNucleus, part::Nucleus);
VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(TypeLepton, part::Lepton);
VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(TypeHadron, part::Hadron);
// Math types.
VALUE_TYPE_DEFINE_PROVIDE_SIMPLE(TypeVec3, math::Vec3);
VALUE_TYPE_DEFINE_PROVIDE_BOUND(TypeBound, math::Bound);

// Stream IO definitions.
// Numbers.
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(TypeDouble, double);
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_SIMPLE(TypeInt, int);
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeRcMethod, RcMethod, 3,
	ESC({ RcMethod::NONE, RcMethod::APPROX, RcMethod::EXACT }),
	ESC({ { "none" }, { "approx" }, { "exact" } }));
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeNucleus, part::Nucleus, 3,
	ESC({ part::Nucleus::P, part::Nucleus::N, part::Nucleus::D }),
	ESC({ { "p", "proton" }, { "n", "neutron" }, { "d", "deuteron" } }));
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeLepton, part::Lepton, 3,
	ESC({ part::Lepton::E, part::Lepton::MU, part::Lepton::TAU }),
	ESC({ { "e", "electron" }, { "mu", "muon" }, { "tau", "tauon" } }));
VALUE_TYPE_DEFINE_READ_WRITE_STREAM_ENUM(
	TypeHadron, part::Hadron, 6,
	ESC({ part::Hadron::PI_0, part::Hadron::PI_P, part::Hadron::PI_M, part::Hadron::K_0, part::Hadron::K_P, part::Hadron::K_M }),
	ESC({ { "pi0", "pion0" }, { "pi+", "pion+" }, { "pi-", "pion-" }, { "K0", "kaon0" }, { "K+", "kaon+" }, { "K-", "kaon-" } }));

// ROOT conversion definitions.
// Numbers.
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeDouble, double, double);
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeInt, int, int);
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeBool, bool, bool);
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeRcMethod, RcMethod, int);
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeNucleus, part::Nucleus, int);
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeLepton, part::Lepton, int);
VALUE_TYPE_DEFINE_CONVERT_ROOT_NUMBER(TypeHadron, part::Hadron, int);

// Can provide.

// Version.
bool TypeVersion::can_provide_base(Version const& prod, Version const& cons) const {
	switch (_provide) {
	case ProvideSpecial::EQ:
		return prod == cons;
	case ProvideSpecial::ANY:
		return true;
	case ProvideSpecial::SPECIAL:
		if (prod.v_major != cons.v_major) {
			return false;
		} else if (prod.v_minor < cons.v_minor) {
			return false;
		} else {
			return true;
		}
	default:
		throw std::runtime_error("Invalid ProvideSpecial.");
	}
}

// Random number seeds.
bool TypeSeedInit::can_provide_base(SeedInit const& prod, SeedInit const& cons) const {
	switch (_provide) {
	case ProvideSpecial::EQ:
		return prod == cons;
	case ProvideSpecial::ANY:
		return true;
	case ProvideSpecial::SPECIAL:
		if (cons.any) {
			return true;
		} else if (prod.any) {
			return false;
		} else {
			return prod.seed == cons.seed;
		}
	default:
		throw std::runtime_error("Invalid ProvideSpecial.");
	}
}
bool TypeSeedGen::can_provide_base(SeedGen const& prod, SeedGen const& cons) const {
	switch (_provide) {
	case ProvideSpecial::EQ:
		return prod == cons;
	case ProvideSpecial::ANY:
		return true;
	case ProvideSpecial::SPECIAL:
		if (cons.any) {
			return true;
		} else if (prod.any) {
			return false;
		} else {
			return prod.seeds == cons.seeds;
		}
	default:
		throw std::runtime_error("Invalid ProvideSpecial.");
	}
}

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
void TypeBool::write_stream_base(std::ostream& os, bool const& value) const {
	os << (value ? "on" : "off");
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
SeedInit TypeSeedInit::read_stream_base(std::istream& is) const {
	std::string str;
	is >> str;
	if (str == "any") {
		return SeedInit();
	} else {
		std::stringstream is_str(str);
		int seed;
		is_str >> seed;
		if (!is_str) {
			is.setstate(std::ios_base::failbit);
		}
		return SeedInit(seed);
	}
}
void TypeSeedInit::write_stream_base(std::ostream& os, SeedInit const& seed) const {
	if (seed.any) {
		os << "any";
	} else {
		os << seed.seed;
	}
}
SeedGen TypeSeedGen::read_stream_base(std::istream& is) const {
	std::string str;
	std::getline(is, str);
	if (str == "any") {
		return SeedGen();
	} else {
		std::stringstream is_str(str);
		SeedGen seed;
		seed.any = false;
		while (is_str && !is_str.eof()) {
			int value;
			is >> value;
			seed.seeds.insert(value);
		}
		if (!is_str) {
			is.setstate(std::ios_base::failbit);
		}
		return seed;
	}
}
void TypeSeedGen::write_stream_base(std::ostream& os, SeedGen const& seed) const {
	if (seed.any) {
		os << "any";
	} else {
		bool first = true;
		for (int value : seed.seeds) {
			if (!first) {
				os << ' ';
				first = false;
			}
			os << value;
		}
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
TArrayI TypeVersion::convert_to_root_base(Version const& version) const {
	Int_t vals[2] = { version.v_major, version.v_minor };
	return TArrayI(2, vals);
}
Version TypeVersion::convert_from_root_base(TArrayI& version) const {
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
TParameter<Int_t> TypeSeedInit::convert_to_root_base(SeedInit const& seed) const {
	if (seed.any) {
		throw std::runtime_error("Can't write seed 'any' to ROOT.");
	} else {
		return TParameter<Int_t>("", seed.seed);
	}
}
SeedInit TypeSeedInit::convert_from_root_base(TParameter<Int_t>& seed) const {
	return SeedInit(seed.GetVal());
}
TArrayI TypeSeedGen::convert_to_root_base(SeedGen const& seed) const {
	if (seed.any) {
		throw std::runtime_error("Can't write seed 'any' to ROOT.");
	} else {
		std::vector<int> vals(seed.seeds.begin(), seed.seeds.end());
		return TArrayI(vals.size(), vals.data());
	}
}
SeedGen TypeSeedGen::convert_from_root_base(TArrayI& seed) const {
	SeedGen result;
	result.any = false;
	for (Int_t i = 0; i < seed.GetSize(); ++i) {
		result.seeds.insert(seed.At(i));
	}
	return result;
}
// Math types.
TArrayD TypeVec3::convert_to_root_base(math::Vec3 const& vec) const {
	Double_t vals[3] = { vec.x, vec.y, vec.z };
	return TArrayD(3, vals);
}
math::Vec3 TypeVec3::convert_from_root_base(TArrayD& vec) const {
	if (vec.GetSize() != 3) {
		throw std::runtime_error("Wrong number of array elements.");
	}
	return math::Vec3(vec.At(0), vec.At(1), vec.At(2));
}
TArrayD TypeBound::convert_to_root_base(math::Bound const& bound) const {
	Double_t vals[3] = { bound.min(), bound.max() };
	return TArrayD(2, vals);
}
math::Bound TypeBound::convert_from_root_base(TArrayD& bound) const {
	if (bound.GetSize() != 2 || !(bound.At(0) <= bound.At(1))) {
		throw std::runtime_error("Wrong number of array elements.");
	}
	return math::Bound(bound.At(0), bound.At(1));
}

