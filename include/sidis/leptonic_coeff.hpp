#ifndef SIDIS_LEPTONIC_COEFF_HPP
#define SIDIS_LEPTONIC_COEFF_HPP

#include "sidis/numeric.hpp"


namespace sidis {

namespace kin {
	struct Kinematics;
	struct KinematicsRad;
}

namespace lep {

/**
 * \defgroup LepGroup Lepton coefficients
 * Coefficients for describing the lepton part of the a SIDIS process. The
 * coefficients come from the various \f$\theta_{i}\f$ coefficients in
 * \cite akushevich2019sidis, including for the Born, AMM, and radiative
 * cross-sections.
 *
 * The coefficients are labeled by a beam and target polarization (see
 * \ref PolSection), but there are some inconsistencies in how the target
 * polarization is labeled. The target `P` polarization is used to indicate a
 * polarization "in-plane", either in the `L` direction or the `T1` direction.
 * The target `U` polarization is used to indicate a target polarization
 * "out-of-plane", either unpolarized or in the `T2` direction. See the
 * \ref XsGroup "cross-section" functions for details about what coefficient
 * sets are needed to compute which cross-sections.
 */
/// \{

/**
 * \defgroup LepBornGroup Born lepton coefficients
 * Coefficients \f$\theta_{i}^{B}\f$ for the Born SIDIS diagram.
 */
/// \{
struct LepBornBaseUU {
	Real theta_1;
	Real theta_2;
	Real theta_3;
	Real theta_4;
	LepBornBaseUU() = default;
	explicit LepBornBaseUU(kin::Kinematics const& kin);
};
struct LepBornBaseUP {
	Real theta_6;
	Real theta_8;
	LepBornBaseUP() = default;
	explicit LepBornBaseUP(kin::Kinematics const& kin);
};
struct LepBornBaseLU {
	Real theta_5;
	LepBornBaseLU() = default;
	explicit LepBornBaseLU(kin::Kinematics const& kin);
};
struct LepBornBaseLP {
	Real theta_7;
	Real theta_9;
	LepBornBaseLP() = default;
	explicit LepBornBaseLP(kin::Kinematics const& kin);
};

struct LepBornUU {
	LepBornBaseUU uu;
	LepBornUU() = default;
	explicit LepBornUU(kin::Kinematics const& kin) : uu(kin) { }
};
struct LepBornUP {
	LepBornBaseUU uu;
	LepBornBaseUP up;
	LepBornUP() = default;
	explicit LepBornUP(kin::Kinematics const& kin) : uu(kin), up(kin) { }
};
struct LepBornLU {
	LepBornBaseUU uu;
	LepBornBaseLU lu;
	LepBornLU() = default;
	explicit LepBornLU(kin::Kinematics const& kin) : uu(kin), lu(kin) { }
};
struct LepBornLP {
	LepBornBaseUU uu;
	LepBornBaseUP up;
	LepBornBaseLU lu;
	LepBornBaseLP lp;
	LepBornLP() = default;
	explicit LepBornLP(kin::Kinematics const& kin) : uu(kin), up(kin), lu(kin), lp(kin) { }
};
/// \}

/**
 * \defgroup LepAmmGroup AMM lepton coefficients
 * Coefficients \f$\theta_{i}^{AMM}\f$ for the AMM cross-section.
 */
/// \{
struct LepAmmBaseUU {
	Real theta_1;
	Real theta_2;
	Real theta_3;
	Real theta_4;
	LepAmmBaseUU() = default;
	explicit LepAmmBaseUU(kin::Kinematics const& kin);
};
struct LepAmmBaseUP {
	Real theta_6;
	Real theta_8;
	LepAmmBaseUP() = default;
	explicit LepAmmBaseUP(kin::Kinematics const& kin);
};
struct LepAmmBaseLU {
	Real theta_5;
	LepAmmBaseLU() = default;
	explicit LepAmmBaseLU(kin::Kinematics const& kin);
};
struct LepAmmBaseLP {
	Real theta_7;
	Real theta_9;
	LepAmmBaseLP() = default;
	explicit LepAmmBaseLP(kin::Kinematics const& kin);
};

struct LepAmmUU {
	LepAmmBaseUU uu;
	LepAmmUU() = default;
	explicit LepAmmUU(kin::Kinematics const& kin) : uu(kin) { }
};
struct LepAmmUP {
	LepAmmBaseUU uu;
	LepAmmBaseUP up;
	LepAmmUP() = default;
	explicit LepAmmUP(kin::Kinematics const& kin) : uu(kin), up(kin) { }
};
struct LepAmmLU {
	LepAmmBaseUU uu;
	LepAmmBaseLU lu;
	LepAmmLU() = default;
	explicit LepAmmLU(kin::Kinematics const& kin) : uu(kin), lu(kin) { }
};
struct LepAmmLP {
	LepAmmBaseUU uu;
	LepAmmBaseUP up;
	LepAmmBaseLU lu;
	LepAmmBaseLP lp;
	LepAmmLP() = default;
	explicit LepAmmLP(kin::Kinematics const& kin) : uu(kin), up(kin), lu(kin), lp(kin) { }
};
/// \}

/**
 * \defgroup LepNradGroup Non-radiative lepton coefficients
 * Coefficients for the non-radiative (Born, vertex correction, vacuum
 * polarization, and soft photon emission) SIDIS diagrams.
 */
/// \{
struct LepNradBaseUU {
	LepBornBaseUU born;
	LepAmmBaseUU amm;
	LepNradBaseUU() = default;
	explicit LepNradBaseUU(kin::Kinematics const& kin);
	LepNradBaseUU(LepBornBaseUU const& born, LepAmmBaseUU const& amm);
};
struct LepNradBaseUP {
	LepBornBaseUP born;
	LepAmmBaseUP amm;
	LepNradBaseUP() = default;
	explicit LepNradBaseUP(kin::Kinematics const& kin);
	LepNradBaseUP(LepBornBaseUP const& born, LepAmmBaseUP const& amm);
};
struct LepNradBaseLU {
	LepBornBaseLU born;
	LepAmmBaseLU amm;
	LepNradBaseLU() = default;
	explicit LepNradBaseLU(kin::Kinematics const& kin);
	LepNradBaseLU(LepBornBaseLU const& born, LepAmmBaseLU const& amm);
};
struct LepNradBaseLP {
	LepBornBaseLP born;
	LepAmmBaseLP amm;
	LepNradBaseLP() = default;
	explicit LepNradBaseLP(kin::Kinematics const& kin);
	LepNradBaseLP(LepBornBaseLP const& born, LepAmmBaseLP const& amm);
};

struct LepNradUU {
	LepNradBaseUU uu;
	LepNradUU() = default;
	explicit LepNradUU(kin::Kinematics const& kin) : uu(kin) { }
	LepNradUU(LepBornUU const& born, LepAmmUU const& amm) : uu(born.uu, amm.uu) { }
};
struct LepNradUP {
	LepNradBaseUU uu;
	LepNradBaseUP up;
	LepNradUP() = default;
	explicit LepNradUP(kin::Kinematics const& kin) : uu(kin), up(kin) { }
	LepNradUP(LepBornUP const& born, LepAmmUP const& amm) : uu(born.uu, amm.uu), up(born.up, amm.up) { }
};
struct LepNradLU {
	LepNradBaseUU uu;
	LepNradBaseLU lu;
	LepNradLU() = default;
	explicit LepNradLU(kin::Kinematics const& kin) : uu(kin), lu(kin) { }
	LepNradLU(LepBornLU const& born, LepAmmLU const& amm) : uu(born.uu, amm.uu), lu(born.lu, amm.lu) { }
};
struct LepNradLP {
	LepNradBaseUU uu;
	LepNradBaseUP up;
	LepNradBaseLU lu;
	LepNradBaseLP lp;
	LepNradLP() = default;
	explicit LepNradLP(kin::Kinematics const& kin) : uu(kin), up(kin), lu(kin), lp(kin) { }
	LepNradLP(LepBornLP const& born, LepAmmLP const& amm) : uu(born.uu, amm.uu), up(born.up, amm.up), lu(born.lu, amm.lu), lp(born.lp, amm.lp) { }
};
/// \}

/**
 * \defgroup LepRadGroup Radiative lepton coefficients
 * Coefficients \f$\theta_{i}^{\,j}\f$ for the radiative SIDIS diagrams.
 */
/// \{
struct LepRadBaseUU {
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
	LepRadBaseUU() = default;
	explicit LepRadBaseUU(kin::KinematicsRad const& kin);
};
struct LepRadBaseUP {
	Real theta_061;
	Real theta_062;
	Real theta_063;
	Real theta_064;
	Real theta_081;
	Real theta_082;
	Real theta_083;
	Real theta_084;
	LepRadBaseUP() = default;
	explicit LepRadBaseUP(kin::KinematicsRad const& kin);
};
struct LepRadBaseLU {
	Real theta_051;
	Real theta_052;
	Real theta_053;
	Real theta_151;
	Real theta_152;
	Real theta_153;
	LepRadBaseLU() = default;
	explicit LepRadBaseLU(kin::KinematicsRad const& kin);
};
struct LepRadBaseLP {
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
	LepRadBaseLP() = default;
	explicit LepRadBaseLP(kin::KinematicsRad const& kin);
};

struct LepRadUU {
	LepRadBaseUU uu;
	LepRadUU() = default;
	explicit LepRadUU(kin::KinematicsRad const& kin) : uu(kin) { }
};
struct LepRadUP {
	LepRadBaseUU uu;
	LepRadBaseUP up;
	LepRadUP() = default;
	explicit LepRadUP(kin::KinematicsRad const& kin) : uu(kin), up(kin) { }
};
struct LepRadLU {
	LepRadBaseUU uu;
	LepRadBaseLU lu;
	LepRadLU() = default;
	explicit LepRadLU(kin::KinematicsRad const& kin) : uu(kin), lu(kin) { }
};
struct LepRadLP {
	LepRadBaseUU uu;
	LepRadBaseUP up;
	LepRadBaseLU lu;
	LepRadBaseLP lp;
	LepRadLP() = default;
	explicit LepRadLP(kin::KinematicsRad const& kin) : uu(kin), up(kin), lu(kin), lp(kin) { }
};
/// \}
/// \}

}
}

#endif

