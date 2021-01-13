#include "utility.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

bool write_progress_bar(std::ostream& os, unsigned percent, unsigned width) {
	if (percent > 100) {
		percent = 100;
	}
	std::string edge_left_str = "[";
	std::string edge_right_str = "]";
	std::string spacing_str = " ";
	std::ostringstream percent_ss;
	percent_ss << std::setw(3) << percent << '%';
	std::string percent_str = percent_ss.str();
	unsigned percent_width = percent_str.size();
	unsigned edge_width = edge_left_str.size() + edge_right_str.size();
	unsigned space_width = spacing_str.size();
	if (width <= percent_width + edge_width + space_width) {
		return false;
	}
	unsigned bar_width = width - percent_width - edge_width - space_width;
	unsigned bar_filled = bar_width * percent / 100;
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
	os << percent_str;
	return os.good();
}

