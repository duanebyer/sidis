#ifndef SIDIS_LEPTONIC_COEFF_HPP
#define SIDIS_LEPTONIC_COEFF_HPP

#include "sidis/numeric.hpp"

namespace sidis {

namespace kin {
	struct Kinematics;
}

namespace lep {

struct LepBornUU {
	Real theta_1;
	Real theta_2;
	Real theta_3;
	Real theta_4;

	explicit LepBornUU(kin::Kinematics kin);
};

struct LepBornUP {
	Real theta_6;
	Real theta_8;

	explicit LepBornUP(kin::Kinematics kin);
};

struct LepBornLU {
	Real theta_5;

	explicit LepBornLU(kin::Kinematics kin);
};

struct LepBornLP {
	Real theta_7;
	Real theta_9;

	explicit LepBornLP(kin::Kinematics kin);
};

}
}

#endif

