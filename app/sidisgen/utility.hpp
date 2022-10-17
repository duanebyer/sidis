#ifndef SIDISGEN_UTILITY_HPP
#define SIDISGEN_UTILITY_HPP

#include <limits>
#include <type_traits>

#include <TArrayD.h>
#include <TArrayI.h>
#include <TArrayL.h>
#include <TArrayL64.h>

// Numeric data types. We want to ensure that integers have 32 bits and longs
// have 64 bits.
using Double = Double_t;
using Int = typename std::conditional<
	std::numeric_limits<Int_t>::digits >= 32,
	Int_t, Long_t>::type;
using Long = typename std::conditional<
	std::numeric_limits<Long_t>::digits >= 64,
	Long_t, Long64_t>::type;

using UInt = typename std::make_unsigned<Int>::type;
using ULong = typename std::make_unsigned<Long>::type;

using RootArrayD = TArrayD;
using RootArrayI = typename std::conditional<
	std::numeric_limits<Int_t>::digits >= 32,
	TArrayI, TArrayL>::type;
using RootArrayL = typename std::conditional<
	std::numeric_limits<Long_t>::digits >= 64,
	TArrayL, TArrayL64>::type;

// All methods for handling radiative corrections.
enum class RcMethod {
	// No radiative corrections.
	NONE,
	// Approximate the leading order radiative corrections.
	APPROX,
	// Calculate the leading order radiative corrections exactly.
	EXACT,
};

// All allowed event types.
enum class EventType {
	// Non-radiative (without photon emission).
	NRAD,
	// Radiative
	RAD,
	// Exclusive
	EXCL,
};

int const NUM_EVENT_TYPES = 3;

// Short identifying names for each type of event.
inline char const* event_type_short_name(EventType type) {
	switch(type) {
	case EventType::NRAD:
		return "nrad";
	case EventType::RAD:
		return "rad";
	case EventType::EXCL:
		return "excl";
	default:
		return "<error>";
	}
}

// Longer identifying names (ex. for error message) for each type of event.
inline char const* event_type_name(EventType type) {
	switch (type) {
	case EventType::NRAD:
		return "non-radiative";
	case EventType::RAD:
		return "radiative";
	case EventType::EXCL:
		return "exclusive";
	default:
		return "<error>";
	}
}

#endif

