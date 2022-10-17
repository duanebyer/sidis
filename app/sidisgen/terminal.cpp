#include "terminal.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

std::vector<std::string> split_lines(std::string input, unsigned width) {
	std::vector<std::string> result;
	result.push_back("");
	// Loop over words.
	unsigned line_width = 0;
	auto word_begin = input.begin();
	while (word_begin != input.end() && *word_begin == ' ') {
		++word_begin;
	}
	while (word_begin != input.end()) {
		auto word_end = std::find(word_begin, input.end(), ' ');
		bool first_word = (line_width == 0);
		unsigned word_width = word_end - word_begin + (first_word ? 0 : 1);
		if (line_width + word_width > width) {
			if (first_word) {
				// Truncate the word.
				word_end = word_begin + width;
			} else {
				// Move to next line.
				line_width = 0;
				result.push_back("");
				first_word = true;
				word_width -= 1;
			}
		}
		// Add word to the current line.
		std::string word = std::string(word_begin, word_end);
		if (!first_word) {
			word = ' ' + word;
		}
		result.back() += word;
		line_width += word_width;
		// Move to next word.
		word_begin = word_end;
		while (word_begin != input.end() && *word_begin == ' ') {
			++word_begin;
		}
	}
	return result;
}

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

bool write_table(
		std::ostream& os,
		std::vector<std::string> entries,
		std::vector<unsigned> widths,
		unsigned spacing) {
	unsigned col_count = widths.size();
	if (col_count == 0 || entries.size() % col_count != 0) {
		return false;
	}
	unsigned row_count = entries.size() / col_count;
	for (unsigned row_idx = 0; row_idx < row_count; ++row_idx) {
		std::vector<std::vector<std::string> > cells;
		unsigned subrow_count = 0;
		for (unsigned col_idx = 0; col_idx < col_count; ++col_idx) {
			cells.push_back(split_lines(
				entries[row_idx * col_count + col_idx],
				widths[col_idx]));
			if (cells.back().size() > subrow_count) {
				subrow_count = cells.back().size();
			}
		}
		for (unsigned subrow_idx = 0; subrow_idx < subrow_count; ++subrow_idx) {
			for (unsigned col_idx = 0; col_idx < col_count; ++col_idx) {
				std::string next = "";
				if (subrow_idx < cells[col_idx].size()) {
					next = cells[col_idx][subrow_idx];
				}
				unsigned excess = widths[col_idx] - next.size()
					+ (col_idx != col_count - 1 ? spacing : 0);
				os << next << std::string(excess, ' ');
			}
			os << std::endl;
		}
	}
	return os.good();
}

