#include "sidis/tmd.hpp"

#include <cmath>
#include <initializer_list>

#include "sidis/extra/exception.hpp"

using namespace sidis;
using namespace sidis::sf;

void FlavorVec::raise_out_of_range(unsigned size) const {
	throw FlavorOutOfRange(size);
}

FlavorVec::FlavorVec(std::initializer_list<Real> list) {
	auto size = list.size();
	if (!(size <= MAX_FLAVOR_VEC_SIZE)) {
		raise_out_of_range(size);
	}
	std::copy(list.begin(), list.end(), _arr);
	_size = size;
}

FlavorVec sf::sqrt_vec(FlavorVec vec) {
	for (unsigned fl = 0; fl < vec.size(); ++fl) {
		vec[fl] = std::sqrt(vec[fl]);
	}
	return vec;
}
FlavorVec sf::sq_vec(FlavorVec vec) {
	for (unsigned fl = 0; fl < vec.size(); ++fl) {
		vec[fl] = vec[fl]*vec[fl];
	}
	return vec;
}
FlavorVec sf::exp_vec(FlavorVec vec) {
	for (unsigned fl = 0; fl < vec.size(); ++fl) {
		vec[fl] = std::exp(vec[fl]);
	}
	return vec;
}
FlavorVec sf::pow_vec(Real lhs, FlavorVec rhs) {
	for (unsigned fl = 0; fl < rhs.size(); ++fl) {
		rhs[fl] = std::pow(lhs, rhs[fl]);
	}
	return rhs;
}
FlavorVec sf::pow_vec(FlavorVec lhs, Real rhs) {
	for (unsigned fl = 0; fl < lhs.size(); ++fl) {
		lhs[fl] = std::pow(lhs[fl], rhs);
	}
	return rhs;
}
FlavorVec sf::pow_vec(FlavorVec lhs, FlavorVec const& rhs) {
	for (unsigned fl = 0; fl < lhs.size(); ++fl) {
		lhs[fl] = std::pow(lhs[fl], rhs[fl]);
	}
	return lhs;
}

FlavorVec sf::tmd_gaussian_factor(FlavorVec var, Real k_perp_sq) {
	return exp_vec(-k_perp_sq/var)/(PI*var);
}

