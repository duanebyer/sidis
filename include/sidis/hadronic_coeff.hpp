#ifndef SIDIS_HADRONIC_COEFF_HPP
#define SIDIS_HADRONIC_COEFF_HPP

#include "sidis/kinematics.hpp"
#include "sidis/numeric.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/vector.hpp"

namespace sidis {
namespace had {

/**
 * \defgroup HadGroup Hadron coefficients
 * Coefficients for describing the hadron part of the SIDIS process. The
 * coefficients come from the \f$\mathcal{H}_{i}\f$ coefficients in
 * \cite akushevich2019sidis.
 *
 * See the \ref XsGroup "cross-section" functions for details about what
 * coefficient sets are needed to compute which cross-sections.
 */
/// \{

/**
 * \defgroup HadStandardGroup Standard hadron coefficients
 * Standard coefficients \f$\mathcal{H}_{i}\f$ for a non-radiative SIDIS
 * process.
 *
 * Structure functions from sf::SfSet are used to calculate these coefficients.
 */
/// \{
struct HadBaseUU {
	Real H_10;
	Real H_20;
	Real H_30;
	Real H_40;
	HadBaseUU() = default;
	HadBaseUU(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseUL {
	Real H_63;
	Real H_83;
	HadBaseUL() = default;
	HadBaseUL(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseUT {
	Real H_12;
	Real H_22;
	Real H_32;
	Real H_42;
	Real H_61;
	Real H_81;
	HadBaseUT() = default;
	HadBaseUT(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseUP {
	HadBaseUL ul;
	HadBaseUT ut;
	HadBaseUP(HadBaseUL const& ul, HadBaseUT const& ut) : ul(ul), ut(ut) { }
	HadBaseUP() = default;
	HadBaseUP(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseLU {
	Real H_50;
	HadBaseLU() = default;
	HadBaseLU(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseLL {
	Real H_73;
	Real H_93;
	HadBaseLL() = default;
	HadBaseLL(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseLT {
	Real H_52;
	Real H_71;
	Real H_91;
	HadBaseLT() = default;
	HadBaseLT(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadBaseLP {
	HadBaseLL ll;
	HadBaseLT lt;
	HadBaseLP(HadBaseLL const& ll, HadBaseLT const& lt) : ll(ll), lt(lt) { }
	HadBaseLP() = default;
	HadBaseLP(kin::Kinematics const& kin, sf::SfSet const& sf);
};

struct HadUU {
	HadBaseUU uu;
	HadUU() = default;
	HadUU(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadUL {
	HadBaseUU uu;
	HadBaseUL ul;
	HadUL() = default;
	HadUL(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadUT {
	HadBaseUU uu;
	HadBaseUT ut;
	HadUT() = default;
	HadUT(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadUP {
	HadBaseUU uu;
	HadBaseUL ul;
	HadBaseUT ut;
	HadUP() = default;
	HadUP(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadLU {
	HadBaseUU uu;
	HadBaseLU lu;
	HadLU() = default;
	HadLU(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadLL {
	HadBaseUU uu;
	HadBaseUL ul;
	HadBaseLU lu;
	HadBaseLL ll;
	HadLL() = default;
	HadLL(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadLT {
	HadBaseUU uu;
	HadBaseUT ut;
	HadBaseLU lu;
	HadBaseLT lt;
	HadLT() = default;
	HadLT(kin::Kinematics const& kin, sf::SfSet const& sf);
};
struct HadLP {
	HadBaseUU uu;
	HadBaseUL ul;
	HadBaseUT ut;
	HadBaseLU lu;
	HadBaseLL ll;
	HadBaseLT lt;
	HadLP() = default;
	HadLP(kin::Kinematics const& kin, sf::SfSet const& sf);
};
/// \}

/** \defgroup HadRadGroup Radiative hadron coefficients
 * These are the shifted versions of the coefficients
 * \f$\tilde{\mathcal{H}}_{i}\f$ calculated through the substitution
 * \f$q \rightarrow q - k\f$. They are needed for computing the radiative
 * cross-section. The radiative cross-section has an infrared divergence in the
 * \f$k \rightarrow 0\f$ case. For removing this divergence, see
 * \ref HadRadFGroup "HadRadFBaseXX"-type coefficients.
 *
 * "Shifted" structure functions from sf::SfSet are used to calculate these
 * coefficients.
 */
/// \{
struct HadRadBaseUU {
	Real H_10;
	Real H_20;
	Real H_30;
	Real H_40;
	HadRadBaseUU() = default;
	HadRadBaseUU(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
struct HadRadBaseUP {
	math::Vec3 H_1;
	math::Vec3 H_2;
	math::Vec3 H_3;
	math::Vec3 H_4;
	math::Vec3 H_6;
	math::Vec3 H_8;
	HadRadBaseUP() = default;
	HadRadBaseUP(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
struct HadRadBaseLU {
	Real H_50;
	HadRadBaseLU() = default;
	HadRadBaseLU(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
struct HadRadBaseLP {
	math::Vec3 H_5;
	math::Vec3 H_7;
	math::Vec3 H_9;
	HadRadBaseLP() = default;
	HadRadBaseLP(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};

struct HadRadUU {
	HadRadBaseUU uu;
	HadRadUU() = default;
	HadRadUU(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
struct HadRadUP {
	HadRadBaseUU uu;
	HadRadBaseUP up;
	HadRadUP() = default;
	HadRadUP(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
struct HadRadLU {
	HadRadBaseUU uu;
	HadRadBaseLU lu;
	HadRadLU() = default;
	HadRadLU(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
struct HadRadLP {
	HadRadBaseUU uu;
	HadRadBaseUP up;
	HadRadBaseLU lu;
	HadRadBaseLP lp;
	HadRadLP() = default;
	HadRadLP(kin::KinematicsRad const& kin, sf::SfSet const& sf);
};
/// \}

/**
 * \defgroup HadRadFGroup Infrared-divergent-free radiative hadron coefficients
 * A variation of the \ref HadRadGroup "HadRadBaseXX"-type coefficients, in
 * which the infrared divergence at \f$k \rightarrow 0\f$ is removed through
 * taking the difference with the unshifted \ref HadStandardGroup "HadXX"-type
 * coefficients
 * (\f$\mathcal{H}_{i}^{F} = \frac{1}{R}(\tilde{\mathcal{H}}_{i} - \mathcal{H}_{i})\f$).
 *
 * "Shifted" structure functions from sf::SfSet are used to calculate these
 * coefficients. It is recommended to supply an sf::SfSet instance directly so
 * that the constructor can calculate the required shift itself.
 */
/// \{
struct HadRadFBaseUU {
	Real H_10;
	Real H_20;
	Real H_30;
	Real H_40;
	Real H_10_diff;
	Real H_20_diff;
	Real H_30_diff;
	Real H_40_diff;
	HadRadFBaseUU() = default;
	HadRadFBaseUU(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadBaseUU const& had_0);
	HadRadFBaseUU(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFBaseUU(kin, sf, HadBaseUU(kin.project(), sf)) { }
};
struct HadRadFBaseUP {
	math::Vec3 H_1;
	math::Vec3 H_2;
	math::Vec3 H_3;
	math::Vec3 H_4;
	math::Vec3 H_6;
	math::Vec3 H_8;
	math::Vec3 H_1_diff;
	math::Vec3 H_2_diff;
	math::Vec3 H_3_diff;
	math::Vec3 H_4_diff;
	math::Vec3 H_6_diff;
	math::Vec3 H_8_diff;
	HadRadFBaseUP() = default;
	HadRadFBaseUP(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadBaseUP const& had_0);
	HadRadFBaseUP(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFBaseUP(kin, sf, HadBaseUP(kin.project(), sf)) { }
};
struct HadRadFBaseLU {
	Real H_50;
	Real H_50_diff;
	HadRadFBaseLU() = default;
	HadRadFBaseLU(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadBaseLU const& had_0);
	HadRadFBaseLU(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFBaseLU(kin, sf, HadBaseLU(kin.project(), sf)) { }
};
struct HadRadFBaseLP {
	math::Vec3 H_5;
	math::Vec3 H_7;
	math::Vec3 H_9;
	math::Vec3 H_5_diff;
	math::Vec3 H_7_diff;
	math::Vec3 H_9_diff;
	HadRadFBaseLP() = default;
	HadRadFBaseLP(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadBaseLP const& had_0);
	HadRadFBaseLP(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFBaseLP(kin, sf, HadBaseLP(kin.project(), sf)) { }
};

struct HadRadFUU {
	HadRadFBaseUU uu;
	HadRadFUU() = default;
	HadRadFUU(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadUU const& had_0);
	HadRadFUU(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFUU(kin, sf, HadUU(kin.project(), sf)) { }
};
struct HadRadFUP {
	HadRadFBaseUU uu;
	HadRadFBaseUP up;
	HadRadFUP() = default;
	HadRadFUP(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadUP const& had_0);
	HadRadFUP(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFUP(kin, sf, HadUP(kin.project(), sf)) { }
};
struct HadRadFLU {
	HadRadFBaseUU uu;
	HadRadFBaseLU lu;
	HadRadFLU() = default;
	HadRadFLU(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadLU const& had_0);
	HadRadFLU(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFLU(kin, sf, HadLU(kin.project(), sf)) { }
};
struct HadRadFLP {
	HadRadFBaseUU uu;
	HadRadFBaseUP up;
	HadRadFBaseLU lu;
	HadRadFBaseLP lp;
	HadRadFLP() = default;
	HadRadFLP(kin::KinematicsRad const& kin, sf::SfSet const& sf, HadLP const& had_0);
	HadRadFLP(kin::KinematicsRad const& kin, sf::SfSet const& sf) :
		HadRadFLP(kin, sf, HadLP(kin.project(), sf)) { }
};
/// \}
/// \}

}
}

#endif

