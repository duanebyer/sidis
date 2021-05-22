#include "event_stats.hpp"

#include <cmath>
#include <limits>

EventStats::EventStats() :
	xs(0.), xs_err(std::numeric_limits<Double_t>::infinity()),
	weight_total(0.),
	weight_sq_total(0.),
	num_events(0) { }

Double_t EventStats::efficiency() const {
	return std::sqrt(
		(weight_total * weight_total)
		/ (weight_sq_total * num_events));
}
Double_t EventStats::norm() const {
	return xs / weight_total;
}

