#include "utility.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

bool write_progress_bar(std::ostream& os, double fraction, unsigned width) {
	if (fraction < 0. || std::isnan(fraction)) {
		fraction = 0.;
	} else if (fraction > 1.) {
		fraction = 1.;
	}
	int percentage = static_cast<int>(fraction * 100.);
	std::string edge_left_str = "[";
	std::string edge_right_str = "]";
	std::string spacing_str = " ";
	std::ostringstream percentage_ss;
	percentage_ss << std::setw(3) << percentage << '%';
	std::string percentage_str = percentage_ss.str();
	unsigned percentage_width = percentage_str.size();
	unsigned edge_width = edge_left_str.size() + edge_right_str.size();
	unsigned space_width = spacing_str.size();
	if (width <= percentage_width + edge_width + space_width) {
		return false;
	}
	unsigned bar_width = width - percentage_width - edge_width - space_width;
	unsigned bar_filled = static_cast<unsigned>(fraction * bar_width);
	os << edge_left_str;
	for (unsigned idx = 0; idx < bar_width; ++idx) {
		if (idx < bar_filled) {
			os << '=';
		} else {
			os << ' ';
		}
	}
	os << edge_right_str;
	os << spacing_str;
	os << percentage_str;
	return os.good();
}

