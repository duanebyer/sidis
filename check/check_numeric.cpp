#include "sidis/numeric.hpp"

#include <type_traits>

static_assert(
	std::is_same<sidis::Real, float>::value
	|| std::is_same<sidis::Real, double>::value
	|| std::is_same<sidis::Real, long double>::value,
	"sidis::Real must be an allowed floating point type");

int main(int argc, char** argv) {
	return 0;
}

