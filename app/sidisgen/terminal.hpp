#ifndef SIDISGEN_TERMINAL_HPP
#define SIDISGEN_TERMINAL_HPP

#include <ostream>
#include <string>
#include <vector>

// Splits a string into lines to prevent going over certain width.
std::vector<std::string> split_lines(std::string input, unsigned width=70);
// Draws a progress bar in the terminal.
bool write_progress_bar(std::ostream& os, unsigned percent, unsigned width=70);
// Draws a table to the terminal.
bool write_table(
	std::ostream& os,
	std::vector<std::string> entries,
	std::vector<unsigned> widths,
	unsigned spacing);

#endif

