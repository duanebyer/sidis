#ifndef SIDISGEN_PARAMS_FORMAT_HPP
#define SIDISGEN_PARAMS_FORMAT_HPP

#include <set>
#include <stdexcept>
#include <string>

#include <TArrayI.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TObjString.h>
#include <TParameter.h>

#include <sidis/bound.hpp>
#include <sidis/particle.hpp>
#include <sidis/vector.hpp>

#include "params.hpp"
#include "utility.hpp"

// The standard format describing the parameter file used by `sidisgen`. The
// format makes use of the types declared below. Any changes to the parameters
// used by `sidisgen` (e.x. adding a new paramter) should be done with this
// variable.
extern Params const PARAMS_STD_FORMAT;

void check_can_provide_foam(Params& params_foam, Params& params_gen);
Params merge_params(Params& params_1, Params& params_2);

// Represents a major.minor version of a parameter file.
struct Version {
	int v_major;
	int v_minor;
	Version(
		int v_major=0,
		int v_minor=0) :
		v_major(v_major),
		v_minor(v_minor) { }
	bool operator==(Version const& rhs) const {
		return v_major == rhs.v_major && v_minor == rhs.v_minor;
	}
	bool operator!=(Version const& rhs) const {
		return !(*this == rhs);
	}
};

// Seed used in the initialization process.
struct SeedInit {
	// If `_any` is true, then the initialization process can choose any val
	// it wants for the seed, ignoring the specified seed.
	bool any;
	int seed;

	SeedInit() : any(true) { }
	SeedInit(int seed) :
		any(false),
		seed(seed) { }
	bool operator==(SeedInit const& rhs) const {
		if (any && rhs.any) {
			return true;
		} else if (!any && !rhs.any) {
			return seed == rhs.seed;
		} else {
			return false;
		}
	}
	bool operator!=(SeedInit const& rhs) const {
		return !(*this == rhs);
	}
};

// Seeds used in the generation process.
struct SeedGen {
	bool any;
	// Set of all seeds used for the generation process.
	std::multiset<int> seeds;

	SeedGen() : any(true) { }
	SeedGen(int seed) :
		any(false),
		seeds{ seed } { }
	// Merge two `SeedGen`s together.
	SeedGen(SeedGen const& seed_1, SeedGen const& seed_2) :
			any(seed_1.any || seed_2.any),
			seeds() {
		// Can't use `set_union` because of unusual behaviour with multisets.
		if (!any) {
			for (int seed : seed_1.seeds) {
				seeds.insert(seed);
			}
			for (int seed : seed_2.seeds) {
				seeds.insert(seed);
			}
		}
	}
	bool operator==(SeedGen const& rhs) const {
		if (any && rhs.any) {
			return true;
		} else if (!any && !rhs.any) {
			return seeds == rhs.seeds;
		} else {
			return false;
		}
	}
	bool operator!=(SeedGen const& rhs) const {
		return !(*this == rhs);
	}
};

// These functions are helpful for getting names of parameters that depend on
// the type of event, e.x. "mc.nrad.enable", "mc.rad.enable", "mc.excl.enable".
inline std::string p_name_enable(EventType ev_type) {
	return std::string("mc.") + event_type_short_name(ev_type) + ".enable";
}
inline std::string p_name_gen_rej_scale(EventType ev_type) {
	return std::string("mc.") + event_type_short_name(ev_type) + ".gen.rej_scale";
}
inline std::string p_name_init_seed(EventType ev_type) {
	return std::string("mc.") + event_type_short_name(ev_type) + ".init.seed";
}
inline std::string p_name_init_target_eff(EventType ev_type) {
	return std::string("mc.") + event_type_short_name(ev_type) + ".init.target_eff";
}
inline std::string p_name_init_scale_exp(EventType ev_type) {
	return std::string("mc.") + event_type_short_name(ev_type) + ".init.scale_exp";
}
inline std::string p_name_init_max_cells(EventType ev_type) {
	return std::string("mc.") + event_type_short_name(ev_type) + ".init.max_cells";
}

std::set<EventType> p_active_event_types(Params const& params);

// Convenience macros for declaring new types.
#define VALUE_TYPE_DECLARE(RType, RValue, Wrapped, WrappedRoot) \
	class RType final : public Type { \
	public: \
		static RType const INSTANCE; \
		bool equivalent_base(Wrapped const& val_1, Wrapped const& val_2) const; \
		Wrapped convert_from_root_base(WrappedRoot& obj) const; \
		WrappedRoot convert_to_root_base(Wrapped const& obj) const; \
		Wrapped read_stream_base(std::istream& is) const; \
		void write_stream_base(std::ostream& os, Wrapped const& val) const; \
		\
		bool equivalent(Value const& val_1, Value const& val_2) const override; \
		std::unique_ptr<Value> read_root(TDirectory& dir, std::string const& name) const override; \
		void write_root(TDirectory& dir, std::string const& name, Value const& val) const override; \
		std::unique_ptr<Value> read_stream(std::istream& is) const override; \
		void write_stream(std::ostream& os, Value const& val) const override; \
	}; \
	class RValue final : public Value { \
	public: \
		Wrapped val; \
		RValue(Wrapped val) : \
			Value(RType::INSTANCE), \
			val(val) { } \
		template<typename... Ts> \
		RValue(Ts... args) : \
			Value(RType::INSTANCE), \
			val(args...) { } \
		Any any() const override { \
			return Any(val); \
		} \
		operator Wrapped&() { \
			return val; \
		} \
		operator Wrapped const&() const { \
			return val; \
		} \
	}; \
	inline bool RType::equivalent_base(Wrapped const& val_1, Wrapped const& val_2) const { \
		return val_1 == val_2; \
	} \
	inline bool RType::equivalent(Value const& val_1, Value const& val_2) const { \
		return equivalent_base(val_1.as<RValue>().val, val_2.as<RValue>().val); \
	} \
	inline std::unique_ptr<Value> RType::read_root(TDirectory& dir, std::string const& name) const { \
		WrappedRoot* obj = dir.Get<WrappedRoot>(name.c_str()); \
		if (obj == nullptr) { \
			throw std::runtime_error("Couldn't find ROOT object."); \
		} else { \
			return std::make_unique<RValue>(convert_from_root_base(*obj)); \
		} \
	} \
	inline void RType::write_root(TDirectory& dir, std::string const& name, Value const& val) const { \
		WrappedRoot obj = convert_to_root_base(val.as<RValue>().val); \
		dir.WriteObject(&obj, name.c_str()); \
	} \
	inline std::unique_ptr<Value> RType::read_stream(std::istream& is) const { \
		return std::make_unique<RValue>(read_stream_base(is)); \
	} \
	inline void RType::write_stream(std::ostream& os, Value const& val) const { \
		RValue const& value_c = dynamic_cast<RValue const&>(val); \
		return write_stream_base(os, value_c.val); \
	}

// Built-in val types.
// Version.
VALUE_TYPE_DECLARE(TypeVersion, ValueVersion, Version, TArrayI);
// Numbers.
VALUE_TYPE_DECLARE(TypeDouble, ValueDouble, Double, TParameter<Double>);
VALUE_TYPE_DECLARE(TypeInt, ValueInt, Int, TParameter<Int>);
VALUE_TYPE_DECLARE(TypeLong, ValueLong, Long, TParameter<Long>);
VALUE_TYPE_DECLARE(TypeBool, ValueBool, bool, TParameter<bool>);
// Strings.
VALUE_TYPE_DECLARE(TypeString, ValueString, std::string, TObjString);
// Random number seeds.
VALUE_TYPE_DECLARE(TypeSeedInit, ValueSeedInit, SeedInit, TParameter<int>);
VALUE_TYPE_DECLARE(TypeSeedGen, ValueSeedGen, SeedGen, TArrayI);
// Enums.
VALUE_TYPE_DECLARE(TypeRcMethod, ValueRcMethod, RcMethod, TParameter<int>);
VALUE_TYPE_DECLARE(TypeNucleus, ValueNucleus, sidis::part::Nucleus, TParameter<int>);
VALUE_TYPE_DECLARE(TypeLepton, ValueLepton, sidis::part::Lepton, TParameter<int>);
VALUE_TYPE_DECLARE(TypeHadron, ValueHadron, sidis::part::Hadron, TParameter<int>);
// Math types.
VALUE_TYPE_DECLARE(TypeVec3, ValueVec3, sidis::math::Vec3, RootArrayD);
VALUE_TYPE_DECLARE(TypeBound, ValueBound, sidis::math::Bound, RootArrayD);

#endif

