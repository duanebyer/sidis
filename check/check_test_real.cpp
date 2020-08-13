#include "sidis/numeric.hpp"

#include <limits>

using Limits = std::numeric_limits<sidis::Real>;

// Require that the `Real` type is at least as capable as the IEEE 754 64-bit
// interchange formats.
static_assert(
	Limits::radix == 2 || Limits::radix == 10,
	"sidis::Real has invalid radix");
static_assert(
	(Limits::radix == 2 && Limits::digits >= 53)
	|| (Limits::radix == 10 && Limits::digits >= 16),
	"sidis::Real has insufficient precision");
static_assert(
	(Limits::radix == 2
		&& Limits::min_exponent <= -1021 && Limits::max_exponent >= 1024)
	|| (Limits::radix == 10
		&& Limits::min_exponent <= -382 && Limits::max_exponent >= 385),
	"sidis::Real has insufficient exponent range");

int main(int argc, char** argv) {
	return 0;
}

