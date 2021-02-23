#include <limits>

using Limits = std::numeric_limits<float>;

// Require that the `float` type is at least as capable as the IEEE 754 32-bit
// interchange formats.
static_assert(
	Limits::radix == 2 || Limits::radix == 10,
	"float has invalid radix");
static_assert(
	(Limits::radix == 2 && Limits::digits >= 24)
	|| (Limits::radix == 10 && Limits::digits >= 7),
	"float has insufficient precision");
static_assert(
	(Limits::radix == 2
		&& Limits::min_exponent <= -125 && Limits::max_exponent >= 128)
	|| (Limits::radix == 10
		&& Limits::min_exponent <= -94 && Limits::max_exponent >= 97),
	"float has insufficient exponent range");
static_assert(
	Limits::is_iec559,
	"float does not satisfy IEEE 754");

int main(int argc, char** argv) {
	return 0;
}

