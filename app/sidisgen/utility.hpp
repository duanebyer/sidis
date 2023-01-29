#ifndef SIDISGEN_UTILITY_HPP
#define SIDISGEN_UTILITY_HPP

#include <cstdlib>
#include <limits>
#include <type_traits>

#include <TArrayD.h>
#include <TArrayI.h>
#include <TArrayL.h>
#include <TArrayL64.h>

// TODO: Have to include this header for the Stats type, would like to move the
// the Stats type into an external library instead.
#include <bubble.hpp>

#if defined(__GNUC__)
#define UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#define UNREACHABLE() __assume(0)
#else
#define UNREACHABLE() exit(128)
#endif

// Numeric data types. We want to ensure that integers have 32 bits and longs
// have 64 bits.
using Double = Double_t;
using Int = typename std::conditional<
	std::numeric_limits<Int_t>::digits >= 32 - 1,
	Int_t, Long_t>::type;
using Long = typename std::conditional<
	std::numeric_limits<Long_t>::digits >= 64 - 1,
	Long_t, Long64_t>::type;

using RootArrayD = TArrayD;
using RootArrayI = typename std::conditional<
	std::numeric_limits<Int_t>::digits >= 32 - 1,
	TArrayI, TArrayL>::type;
using RootArrayL = typename std::conditional<
	std::numeric_limits<Long_t>::digits >= 64 - 1,
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
};

int const NUM_EVENT_TYPES = 2;

inline constexpr std::size_t event_type_dimension(EventType ev_type) {
	if (ev_type == EventType::NRAD) {
		return 6;
	} else if (ev_type == EventType::RAD) {
		return 9;
	} else {
		return 0;
	}
}

// Short identifying names for each type of event.
inline char const* event_type_short_name(EventType type) {
	switch(type) {
	case EventType::NRAD:
		return "nrad";
	case EventType::RAD:
		return "rad";
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
	default:
		return "<error>";
	}
}

// Store statistics (count and moments) of the generated events.
using Stats = bubble::Stats<Double>;
using StatsAccum = bubble::StatsAccum<Double>;

// Accumulates events drawn from a generator to calculate the integral of the
// density function. Allows events from different generators with the same
// density to be combined. Provides an event normalization that can be used when
// computing weighted integrals over the density.
class Integrator final {
	// Normalization, such that total integral = mean event weight / prime.
	Double _prime_inv;
	// Total number of generated events.
	std::size_t _count;
	// Number of accepted events. Events with weight zero are considered
	// "rejected", and all others are considered "accepted".
	std::size_t _count_acc;
	// Statistics of all generated event weights, both accepted and rejected.
	Stats _weights;

public:
	Integrator(Double prime_inv, std::size_t count_acc, Stats weights) :
		_prime_inv(prime_inv),
		_count(weights.count()),
		_count_acc(count_acc),
		_weights(weights) { }
	Integrator() :
		_prime_inv(std::numeric_limits<Double>::quiet_NaN()),
		_count(0), _count_acc(0),
		_weights() { }

	// Combines two integrators together. Integrators with different primes can
	// be combined as well, with the resulting integrator having a prime
	// somewhere in the middle.
	Integrator& operator+=(Integrator const& other) {
		if (_count == 0) {
			*this = other;
		} else if (other._count == 0) {
			// No modification needed.
		} else {
			// The rule for prime combination is `N/P = N1/P1 + N2/P2`. This
			// ensures that expectation values remain unchanged.
			Double lambda = static_cast<Double>(other._count) / (_count + other._count);
			_prime_inv += lambda * (other._prime_inv - _prime_inv);
			_count += other._count;
			_count_acc += other._count_acc;
			_weights += other._weights;
		}
		return *this;
	}

	// Prime of the underlying generator.
	Double prime() const {
		return 1. / _prime_inv;
	}
	// Total number of events.
	std::size_t count() const {
		return _count;
	}
	// Number of accepted events (those having non-zero weight).
	std::size_t count_acc() const {
		return _count_acc;
	}
	Double acc() const {
		return static_cast<Double>(_count_acc) / _count;
	}
	// Normalization constant for events.
	Double norm() const {
		return _count * _prime_inv;
	}
	// Statistics describing the event weights.
	Stats weights() const {
		return _weights;
	}
	// Statistics describing the weights of accepted events.
	Stats weights_acc() const {
		Stats weights_scaled = _weights;
		weights_scaled.rescale_count(_count_acc);
		return weights_scaled;
	}
	// Estimate of integral. Note: error estimate assumes independent events.
	Double integ(Double* err=nullptr) const {
		Double mean_err;
		Double mean = _weights.est_mean(&mean_err);
		if (err != nullptr) {
			*err = mean_err / _prime_inv;
		}
		return mean / _prime_inv;
	}
	// Estimate of relative variance of event weights, which is used as a
	// measure of efficiency. Note: error estimate assumes independent events.
	Double rel_var(Double* err=nullptr) const {
		return weights_acc().est_rel_var(err);
	}
};
inline Integrator operator+(Integrator lhs, Integrator const& rhs) {
	lhs += rhs;
	return lhs;
}

// Accumulates events drawn one-by-one from a generator to estimate the integral
// of the density function. After acquiring all of the events, can be totalled
// to produce an `Integrator`.
class IntegratorAccum final {
	// Members similar to `Integrator`.
	Double _prime_inv;
	std::size_t _count;
	std::size_t _count_acc;
	// Use an accumulator for the statistics.
	StatsAccum _weights;

public:
	// Creates a new integrator for an event generator with a given prime. The
	// prime is usually determined during the process of constructing the event
	// generator, see `build_dist_approx`.
	IntegratorAccum(Double prime_inv) :
		_prime_inv(prime_inv),
		_count(0),
		_count_acc(0),
		_weights() { }

	// Adds a new event from the event generator to the integral. The event must
	// be drawn from a generator with the same prime as provided in the
	// constructor.
	IntegratorAccum& operator+=(Double weight) {
		if (weight != 0.) {
			_count_acc += 1;
		}
		_count += 1;
		_weights += weight;
		return *this;
	}

	// Total up the measurements and get the result.
	Integrator total() const {
		return Integrator(_prime_inv, _count_acc, _weights.total());
	}
	// A faster but less accurate version.
	Integrator total_fast() const {
		return Integrator(_prime_inv, _count_acc, _weights.total_fast());
	}

	Double prime() const {
		return 1. / _prime_inv;
	}
	std::size_t count() const {
		return _count;
	}
	std::size_t count_acc() const {
		return _count_acc;
	}
	Double acc() const {
		return static_cast<Double>(_count_acc) / _count;
	}
	Double norm() const {
		return _count * _prime_inv;
	}
};

#endif

