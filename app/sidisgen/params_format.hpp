#ifndef SIDISGEN_PARAMS_FORMAT_HPP
#define SIDISGEN_PARAMS_FORMAT_HPP

#include "params.hpp"

#include <set>
#include <stdexcept>
#include <string>

#include <TArrayD.h>
#include <TArrayI.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TObjString.h>
#include <TParameter.h>

#include <sidis/bound.hpp>
#include <sidis/particle.hpp>
#include <sidis/vector.hpp>

// Returns the standard format describing the parameter file used by `sidisgen`.
// The format makes use of the types declared below. Any changes to the
// parameters used by `sidisgen` (e.x. adding a new paramter) should be done in
// this function.
Params params_format();

// These enums are used for defining the `can_provide` relationship between
// values. For example, a producer with `ProvideOrder::LT` can only provide an
// object to a consumer with a smaller value than itself. An example of this
// relationship is FOAM efficiency: a FOAM with high efficiency can be used to
// generate events for a parameter file requesting a lower efficiency.
//
// In the context of `sidisgen`, the 'producer' creates a FOAM for the
// 'consumer' to use to generate events. If `can_provide(producer, consumer)` is
// true, then the FOAM is compatible with the consumer.

// The basic relationships between producers and consumers. Either the producer
// can only provide to consumers with the same value, or the producer can
// provide to any consumer.
enum class ProvideSimple {
	EQ,
	ANY,
};
// For types which can't really be described by any of the other relationships
// here. There is a 'special' relationship that works uniquely for this type.
// Example: Versions have their own custom compatibility rules.
enum class ProvideSpecial {
	EQ,
	ANY,
	SPECIAL,
};
// Producer can provide to consumer depending on the ordering between them.
// Example: Produced FOAM can only be consumed if requested efficiency is
// smaller (or equal) to what it was produced with.
enum class ProvideOrder {
	EQ,
	ANY,
	LT_EQ,
	GT_EQ,
};
// Producer can provide to consumer depending on a boolean relationship.
// Example: Strictness level of consumer must be either less strict or as strict
// as producer (LT_EQ relationship).
enum class ProvideBool {
	EQ,
	ANY,
	AND,
	OR,
	LT_EQ,
	GT_EQ,
};
// Producer can provide to consumer depending on the nesting of their intervals.
// Example: producer with cuts on kinematics can only provide to consumers with
// tighter cuts.
enum class ProvideBound {
	EQ,
	ANY,
	IN,
	OUT,
};

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
	// If `_any` is true, then the initialization process can choose any value
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
		} else {
			return seed == rhs.seed;
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
	bool operator==(SeedGen const& rhs) const {
		if (any && rhs.any) {
			return true;
		} else {
			return seeds == rhs.seeds;
		}
	}
	bool operator!=(SeedGen const& rhs) const {
		return !(*this == rhs);
	}
};

enum class RcMethod {
	NONE,
	APPROX,
	EXACT,
};

// Convenience macros for declaring new types.
#define VALUE_TYPE_DECLARE(RType, RValue, Wrapped, WrappedRoot, Provide, PROVIDE_DEFAULT) \
	class RType final : public Type { \
		Provide _provide; \
		RType(Provide provide) : _provide(provide) { } \
	public: \
		static RType const& instance(Provide provide=Provide::PROVIDE_DEFAULT); \
		bool equivalent_base(Wrapped const& value_1, Wrapped const& value_2) const; \
		bool can_provide_base(Wrapped const& prod, Wrapped const& cons) const; \
		Wrapped convert_from_root_base(WrappedRoot& obj) const; \
		WrappedRoot convert_to_root_base(Wrapped const& obj) const; \
		Wrapped read_stream_base(std::istream& is) const; \
		void write_stream_base(std::ostream& os, Wrapped const& value) const; \
		\
		bool equivalent(Value const& value_1, Value const& value_2) const override; \
		bool can_provide(Value const& prod, Value const& cons) const override; \
		std::unique_ptr<Value> read_root(TDirectory& dir, char const* name) const override; \
		void write_root(TDirectory& dir, char const* name, Value const& value) const override; \
		std::unique_ptr<Value> read_stream(std::istream& is) const override; \
		void write_stream(std::ostream& os, Value const& value) const override; \
	}; \
	class RValue final : public Value { \
	public: \
		Wrapped value; \
		RValue(Provide provide, Wrapped value) : \
			Value(RType::instance(provide)), \
			value(value) { } \
		RValue(Wrapped value) : \
			Value(RType::instance(Provide::PROVIDE_DEFAULT)), \
			value(value) { } \
		template<typename... Ts> \
		RValue(Provide provide, Ts... args) : \
			Value(RType::instance(provide)), \
			value(args...) { } \
		template<typename... Ts> \
		RValue(Ts... args) : \
			Value(RType::instance(Provide::PROVIDE_DEFAULT)), \
			value(args...) { } \
	}; \
	inline bool RType::equivalent_base(Wrapped const& value_1, Wrapped const& value_2) const { \
		return value_1 == value_2; \
	} \
	inline bool RType::equivalent(Value const& value_1, Value const& value_2) const { \
		return equivalent_base(value_1.as<RValue>().value, value_2.as<RValue>().value); \
	} \
	inline bool RType::can_provide(Value const& prod, Value const& cons) const { \
		return can_provide_base(prod.as<RValue>().value, cons.as<RValue>().value); \
	} \
	inline std::unique_ptr<Value> RType::read_root(TDirectory& dir, char const* name) const { \
		WrappedRoot* obj = dir.Get<WrappedRoot>(name); \
		if (obj == nullptr) { \
			throw std::runtime_error("Couldn't find ROOT object."); \
		} else { \
			return std::make_unique<RValue>(convert_from_root_base(*obj)); \
		} \
	} \
	inline void RType::write_root(TDirectory& dir, char const* name, Value const& value) const { \
		WrappedRoot obj = convert_to_root_base(value.as<RValue>().value); \
		dir.WriteObject(&obj, name); \
	} \
	inline std::unique_ptr<Value> RType::read_stream(std::istream& is) const { \
		return std::make_unique<RValue>(read_stream_base(is)); \
	} \
	inline void RType::write_stream(std::ostream& os, Value const& value) const { \
		RValue const& value_c = dynamic_cast<RValue const&>(value); \
		return write_stream_base(os, value_c.value); \
	}

// Built-in value types.
// Version.
VALUE_TYPE_DECLARE(TypeVersion, ValueVersion, Version, TArrayI, ProvideSpecial, SPECIAL);
// Numbers.
VALUE_TYPE_DECLARE(TypeDouble, ValueDouble, double, TParameter<double>, ProvideOrder, EQ);
VALUE_TYPE_DECLARE(TypeInt, ValueInt, int, TParameter<int>, ProvideOrder, EQ);
VALUE_TYPE_DECLARE(TypeBool, ValueBool, bool, TParameter<bool>, ProvideBool, EQ);
// Strings.
VALUE_TYPE_DECLARE(TypeString, ValueString, std::string, TObjString, ProvideSimple, EQ);
// Random number seeds.
VALUE_TYPE_DECLARE(TypeSeedInit, ValueSeedInit, SeedInit, TParameter<int>, ProvideSpecial, SPECIAL);
VALUE_TYPE_DECLARE(TypeSeedGen, ValueSeedGen, SeedGen, TArrayI, ProvideSpecial, SPECIAL);
// Enums.
VALUE_TYPE_DECLARE(TypeRcMethod, ValueRcMethod, RcMethod, TParameter<int>, ProvideSimple, EQ);
VALUE_TYPE_DECLARE(TypeNucleus, ValueNucleus, sidis::part::Nucleus, TParameter<int>, ProvideSimple, EQ);
VALUE_TYPE_DECLARE(TypeLepton, ValueLepton, sidis::part::Lepton, TParameter<int>, ProvideSimple, EQ);
VALUE_TYPE_DECLARE(TypeHadron, ValueHadron, sidis::part::Hadron, TParameter<int>, ProvideSimple, EQ);
// Math types.
VALUE_TYPE_DECLARE(TypeVec3, ValueVec3, sidis::math::Vec3, TArrayD, ProvideSimple, EQ);
VALUE_TYPE_DECLARE(TypeBound, ValueBound, sidis::math::Bound, TArrayD, ProvideBound, EQ);

#endif

