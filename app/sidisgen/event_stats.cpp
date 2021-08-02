#include "event_stats.hpp"

#include <cmath>
#include <limits>

EventStats::EventStats() :
	xs(0.), xs_err(std::numeric_limits<Double_t>::infinity()),
	weight_total(0.),
	weight_sq_total(0.),
	num_events(0) { }

Double_t EventStats::efficiency() const {
	// TODO: The efficiency estimate is bad for certain weight distributions. It
	// will tend to overestimate the efficiency.
	return std::sqrt((weight_total * weight_total)
		/ (num_events * weight_sq_total));
}
Double_t EventStats::efficiency_err() const {
	Long_t n = num_events;
	Double_t m1 = weight_total / n;
	Double_t m2 = weight_sq_total / n;
	Double_t m3 = weight_p3_total / n;
	Double_t m4 = weight_p4_total / n;
	Double_t cm2 = m2 - m1 * m1;
	Double_t cm4 = m4 - 4. * m1 * m3 + 6. * m1 * m1 * m2 - 3. * m1 * m1 * m1 * m1;
	Double_t eff = efficiency();
	Double_t rel_err = std::sqrt((cm4 / (cm2 * cm2) + 4. * cm2 / (m1 * m1) - 1.) / n);
	// TODO: This error estimate is really not very good! It substantially
	// underestimates the error, especially for low n.
	Double_t err = rel_err * cm2 / (m1 * m1);
	return 0.5 * eff * eff * eff * err;
}
Double_t EventStats::norm() const {
	return xs / weight_total;
}

