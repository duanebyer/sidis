#ifndef __SIDIS_HADRONIC_COEFF_HPP__
#define __SIDIS_HADRONIC_COEFF_HPP__

#include "sidis/kinematics.hpp"
#include "sidis/numeric.hpp"
#include "sidis/structure_function.hpp"

namespace sidis {
namespace had {

struct HadUU {
	Real H_1;
	Real H_2;
	Real H_3;
	Real H_4;

	HadUU(kin::Kinematics kin, sf::SfUU sf);
};

struct HadUL {
	Real H_6;
	Real H_8;

	HadUL(kin::Kinematics kin, sf::SfUL sf);
};

struct HadUT1 {
	Real H_6;
	Real H_8;

	HadUT1(kin::Kinematics kin, sf::SfUT sf);
};

struct HadUT2 {
	Real H_1;
	Real H_2;
	Real H_3;
	Real H_4;

	HadUT2(kin::Kinematics kin, sf::SfUT sf);
};

struct HadLU {
	Real H_5;

	HadLU(kin::Kinematics kin, sf::SfLU sf);
};

struct HadLL {
	Real H_7;
	Real H_9;

	HadLL(kin::Kinematics kin, sf::SfLL sf);
};

struct HadLT1 {
	Real H_7;
	Real H_9;

	HadLT1(kin::Kinematics kin, sf::SfLT sf);
};

struct HadLT2 {
	Real H_5;

	HadLT2(kin::Kinematics kin, sf::SfLT sf);
};

}
}

#endif

