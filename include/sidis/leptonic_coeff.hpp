#ifndef SIDIS_LEPTONIC_COEFF_HPP
#define SIDIS_LEPTONIC_COEFF_HPP

#include "sidis/numeric.hpp"

namespace sidis {

namespace kin {
	struct Kinematics;
	struct KinematicsRad;
}

namespace lep {

// Born coefficients.
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

// AMM coefficients.
struct LepAmmUU {
	Real theta_1;
	Real theta_2;
	Real theta_3;
	Real theta_4;

	explicit LepAmmUU(kin::Kinematics kin);
};
struct LepAmmUP {
	Real theta_6;
	Real theta_8;

	explicit LepAmmUP(kin::Kinematics kin);
};
struct LepAmmLU {
	Real theta_5;

	explicit LepAmmLU(kin::Kinematics kin);
};
struct LepAmmLP {
	Real theta_7;
	Real theta_9;

	explicit LepAmmLP(kin::Kinematics kin);
};

// Radiative coefficients.
struct LepRadUU {
	Real theta_011;
	Real theta_012;
	Real theta_013;
	Real theta_021;
	Real theta_022;
	Real theta_023;
	Real theta_031;
	Real theta_032;
	Real theta_033;
	Real theta_041;
	Real theta_042;
	Real theta_043;

	explicit LepRadUU(kin::KinematicsRad kin);
};
struct LepRadUP {
	Real theta_061;
	Real theta_062;
	Real theta_063;
	Real theta_064;
	Real theta_081;
	Real theta_082;
	Real theta_083;
	Real theta_084;

	explicit LepRadUP(kin::KinematicsRad kin);
};
struct LepRadLU {
	Real theta_051;
	Real theta_052;
	Real theta_053;
	Real theta_151;
	Real theta_152;
	Real theta_153;

	explicit LepRadLU(kin::KinematicsRad kin);
};
struct LepRadLP {
	Real theta_071;
	Real theta_072;
	Real theta_073;
	Real theta_074;
	Real theta_171;
	Real theta_172;
	Real theta_173;
	Real theta_174;
	Real theta_091;
	Real theta_092;
	Real theta_093;
	Real theta_094;
	Real theta_191;
	Real theta_192;
	Real theta_193;
	Real theta_194;

	explicit LepRadLP(kin::KinematicsRad kin);
};

}
}

#endif

