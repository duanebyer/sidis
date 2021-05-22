#ifndef SIDISGEN_EVENT_STATS
#define SIDISGEN_EVENT_STATS

#include <cmath>

#include <Rtypes.h>

enum class EventType {
	NRAD,
	RAD,
	EXCL,
};

static std::size_t const NUM_TYPES = 3;

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

inline unsigned event_type_dim(EventType type) {
	switch (type) {
	case EventType::NRAD:
		return 6;
	case EventType::RAD:
		return 9;
	case EventType::EXCL:
		return 8;
	default:
		return 0;
	}
}

// Structure for storing statistics related to generated events.
struct EventStats {
	Double_t xs;
	Double_t xs_err;
	Double_t weight_total;
	Double_t weight_sq_total;
	Long_t num_events;

	EventStats();

	// Efficiency is the ratio of the error bars in an "ideal" unweighted MC
	// integration to the error bars of the actual weighted MC integration.
	Double_t efficiency() const;
	// The normalization for events from this generator (to be combined with
	// individual event weights) to ensure that the total integral is correct.
	Double_t norm() const;

	template<typename It>
	static EventStats average(It begin, It end) {
		EventStats result;
		Double_t var_inv_total = 0.;
		Double_t xs_var_inv = 0.;
		for (auto it = begin; it != end; ++it) {
			EventStats stats = *it;
			Double_t var_inv = 1. / (stats.xs_err * stats.xs_err);
			var_inv_total += var_inv;
			xs_var_inv += var_inv * stats.xs;
			result.weight_total += stats.weight_total;
			result.weight_sq_total += stats.weight_sq_total;
			result.num_events += stats.num_events;
		}
		result.xs = xs_var_inv / var_inv_total;
		result.xs_err = std::sqrt(1. / var_inv_total);
		return result;
	}

	template<typename It>
	static EventStats total(It begin, It end) {
		EventStats result;
		Double_t xs_err_sq = 0.;
		for (auto it = begin; it != end; ++it) {
			EventStats stats = *it;
			result.xs += stats.xs;
			xs_err_sq += stats.xs_err * stats.xs_err;
			result.weight_total += stats.weight_total;
			result.weight_sq_total += stats.weight_sq_total;
			result.num_events += stats.num_events;
		}
		result.xs_err = std::sqrt(xs_err_sq);
		return result;
	}
};

#endif

