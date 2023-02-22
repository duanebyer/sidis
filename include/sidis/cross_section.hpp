#ifndef SIDIS_CROSS_SECTION_HPP
#define SIDIS_CROSS_SECTION_HPP

#include "sidis/constant.hpp"
#include "sidis/integ_params.hpp"
#include "sidis/numeric.hpp"

namespace sidis {

namespace math {
	struct Vec3;
}
namespace kin {
	struct Kinematics;
	struct KinematicsRad;
}
namespace ph {
	struct Phenom;
}
namespace lep {
	struct LepBornBaseUU;
	struct LepBornBaseUP;
	struct LepBornBaseLU;
	struct LepBornBaseLP;

	struct LepAmmBaseUU;
	struct LepAmmBaseUP;
	struct LepAmmBaseLU;
	struct LepAmmBaseLP;

	struct LepNradBaseUU;
	struct LepNradBaseUP;
	struct LepNradBaseLU;
	struct LepNradBaseLP;

	struct LepRadBaseUU;
	struct LepRadBaseUP;
	struct LepRadBaseLU;
	struct LepRadBaseLP;
}
namespace had {
	struct HadBaseUU;
	struct HadBaseUL;
	struct HadBaseUT;
	struct HadBaseLU;
	struct HadBaseLL;
	struct HadBaseLT;

	struct HadRadBaseUU;
	struct HadRadBaseUP;
	struct HadRadBaseLU;
	struct HadRadBaseLP;

	struct HadRadFBaseUU;
	struct HadRadFBaseUP;
	struct HadRadFBaseLU;
	struct HadRadFBaseLP;
}
namespace sf {
	class SfSet;
}

namespace xs {

math::IntegParams const DEFAULT_INTEG_PARAMS {
	math::IntegMethod::CUBATURE,
	100000,
	10000,
};

/**
 * \defgroup XsGroup Cross-sections
 * Functions used to calculate differential SIDIS cross-sections. The various
 * kinds of cross-sections correspond to the following diagrams:
 *
 * \image html sidis_rc_feynman_diagrams.svg "Leading order diagrams for radiative corrections to the SIDIS cross-section."
 * \image latex sidis_rc_feynman_diagrams.ps "Leading order diagrams for radiative corrections to the SIDIS cross-section."
 *
 * The kinds of cross-sections are then:
 * * \f$\sigma_{B}\f$: The %Born cross-section.
 * * \f$\sigma_{AMM}\f$: The anomalous magnetic moment cross-section,
 *   related to the vertex correction diagram.
 * * \f$\frac{\alpha}{\pi}\delta_{\text{vac}} \sigma_{B}\f$: A correction factor
 *   to %Born cross-section, corresponding to vacuum polarization diagram.
 * * \f$\frac{\alpha}{\pi}\delta_{\text{vert}} \sigma_{B}\f$: A correction factor
 *   to %Born cross-section, related to vertex correction diagram.
 * * \f$\sigma_{R}\f$: The radiative cross-section, corresponding to radiative
 *   diagrams with photon energy above the soft photon cutoff. Has an infrared
 *   divergence, allowing it to be split into an infrared-divergent part
 *   \f$\sigma_{R}^{IR}\f$ and an infrared-divergent-free part
 *   \f$\sigma_{R}^{F}\f$.
 * * \f$\sigma_{\text{nrad}}\f$: The non-radiative cross-section, corresponding
 *   to %Born, vertex correction, vacuum polarization, and radiative diagrams
 *   with photon energy below the soft photon cutoff.
 * * \f$\sigma_{\text{nrad}}^{IR}\f$: The non-radiative cross-section,
 *   neglecting the difficult-to-compute contribution from \f$\sigma_{R}^{F}\f$,
 *   which is a good approximation to \f$\sigma_{\text{nrad}}\f$ so long as the
 *   soft photon cutoff is small enough.
 *
 * The cross-sections are computed as:
 *
 * \f{equation}{
 *     \sigma = \sigma^{UU} + \pmb{\eta}\cdot\pmb{\sigma}^{UP} + \lambda_e\left(\sigma^{LU} + \pmb{\eta}\cdot\pmb{\sigma}^{LP}\right)
 * \f}
 *
 * Where:
 *
 * \f{equation}{
 *     \pmb{\sigma}^{XP} = \pmb{e}_x\sigma^{XT_1} + \pmb{e}_y\sigma^{XT_2} + \pmb{e}_z\sigma^{XL}
 * \f}
 *
 * For computing the complete cross-section \f$\sigma\f$, use the general
 * cross-section functions like xs::born(). For computing the cross-section
 * parts, such as \f$\sigma^{UT_1}\f$, use the
 * \ref XsBornGroup "base Born cross-section functions".
 */

/**
 * \defgroup GeneralXsGroup General cross-section functions
 * Functions for doing various kinds of cross-section calculations. They take a
 * kin::Kinematics, an sf::SfSet for structure functions, and the beam and
 * target polarizations. Additionally, a ph::Phenom object may (optionally) be
 * provided, to use custom phenomenological inputs.
 * \ingroup XsGroup
 */
/// \{

/// %Born cross-section \f$\sigma_{B}\f$.
Real born(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// \copydoc born()
Real born(kin::Kinematics const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// Anomalous magnetic moment cross-section \f$\sigma_{AMM}\f$, related to
/// vertex correction diagram.
Real amm(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// \copydoc amm()
Real amm(kin::Kinematics const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// Non-radiative cross-section neglecting infrared-divergent-free soft photon
/// part \f$\sigma_{\text{nrad}}^{IR}\f$.
Real nrad_ir(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF);
/// \copydoc nrad_ir()
Real nrad_ir(kin::Kinematics const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF);
/// Non-radiative cross-section \f$\sigma_{\text{nrad}}\f$, integrated over the
/// radiated photon with energy below soft cutoff \p k_0_bar (if \p k_0_bar is
/// set to infinity (default), then the entire radiative part is integrated
/// over).
math::EstErr nrad_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// \copydoc nrad_integ()
math::EstErr nrad_integ(kin::Kinematics const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// Radiative cross-section with infrared divergence removed
/// \f$\sigma_{R}^{F}\f$.
Real rad_f(kin::KinematicsRad const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// \copydoc rad_f()
Real rad_f(kin::KinematicsRad const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// Radiative cross-section \f$\sigma_{R}\f$.
Real rad(kin::KinematicsRad const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// \copydoc rad()
Real rad(kin::KinematicsRad const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);

/// Radiative cross-section with infrared divergence removed
/// \f$\sigma_{R}^{F}\f$, integrated over the radiated photon with energy below
/// soft cutoff \p k_0_bar.
math::EstErr rad_f_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// \copydoc rad_f_integ()
math::EstErr rad_f_integ(kin::Kinematics const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// Radiative cross-section \f$\sigma_{R}\f$, integrated over the radiated
/// photon with energy above soft cutoff \p k_0_bar.
math::EstErr rad_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// \copydoc rad_integ()
math::EstErr rad_integ(kin::Kinematics const& kin, ph::Phenom const& phenom, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// \}

/// \name Born correction factors
/// These correction factors to the %Born cross-section give the contribution
/// from vacuum polarization and from soft radiated photon.
/// \{

/// Factor \f$\frac{\alpha}{\pi}\delta_{VS}\sigma_{B}\f$ gives part of the
/// vertex correction together with the integral over the soft photon part of
/// \f$\sigma_{R}^{IR}\f$. The infrared divergences in these two terms cancel
/// for positive soft photon threshold, giving a finite result.
Real delta_vert_rad_ir(kin::Kinematics const& kin, Real k_0_bar=INF);
/// Factor \f$\frac{\alpha}{\pi}\delta_{H}\sigma_{B}\f$ gives the integral over
/// the hard photon part of \f$\sigma_{R}^{IR}\f$. For the most part, this
/// function isn't very useful, but it is provided for completeness with
/// xs::delta_vert_rad_ir.
Real delta_rad_ir_hard(kin::Kinematics const& kin, Real k_0_bar=INF);
/// Factor \f$\frac{\alpha}{\pi}\delta_{\text{vac}}^{\text{lep}}\sigma_{B}\f$
/// gives the vacuum polarization cross-section due to lepton loops.
Real delta_vac_lep(kin::Kinematics const& kin);
/// Factor \f$\frac{\alpha}{\pi}\delta_{\text{vac}}^{\text{had}}\sigma_{B}\f$
/// gives the vacuum polarization cross-section due to hadron loops.
Real delta_vac_had(kin::Kinematics const& kin);
/// \}

/**
 * \defgroup XsBornGroup Base Born cross-sections
 * %Born cross-section base functions.
 * \sa LepBornGroup
 * \sa HadGroup
 * \ingroup XsGroup
 */
/// \{

struct Born {
	Real coeff;
	Born(kin::Kinematics const& kin, ph::Phenom const& phenom);
};

Real born_base_uu(Born const& b, lep::LepBornBaseUU const& lep, had::HadBaseUU const& had);
Real born_base_ul(Born const& b, lep::LepBornBaseUP const& lep, had::HadBaseUL const& had);
Real born_base_ut1(Born const& b, lep::LepBornBaseUP const& lep, had::HadBaseUT const& had);
Real born_base_ut2(Born const& b, lep::LepBornBaseUU const& lep, had::HadBaseUT const& had);
Real born_base_lu(Born const& b, lep::LepBornBaseLU const& lep, had::HadBaseLU const& had);
Real born_base_ll(Born const& b, lep::LepBornBaseLP const& lep, had::HadBaseLL const& had);
Real born_base_lt1(Born const& b, lep::LepBornBaseLP const& lep, had::HadBaseLT const& had);
Real born_base_lt2(Born const& b, lep::LepBornBaseLU const& lep, had::HadBaseLT const& had);
/// \}

/**
 * \defgroup XsAmmGroup Base AMM cross-sections
 * AMM cross-section base functions.
 * \sa LepAmmGroup
 * \sa HadGroup
 * \ingroup XsGroup
 */
/// \{

struct Amm {
	Real coeff;
	Amm(kin::Kinematics const& kin, ph::Phenom const& phenom);
};

Real amm_base_uu(Amm const& b, lep::LepAmmBaseUU const& lep, had::HadBaseUU const& had);
Real amm_base_ul(Amm const& b, lep::LepAmmBaseUP const& lep, had::HadBaseUL const& had);
Real amm_base_ut1(Amm const& b, lep::LepAmmBaseUP const& lep, had::HadBaseUT const& had);
Real amm_base_ut2(Amm const& b, lep::LepAmmBaseUU const& lep, had::HadBaseUT const& had);
Real amm_base_lu(Amm const& b, lep::LepAmmBaseLU const& lep, had::HadBaseLU const& had);
Real amm_base_ll(Amm const& b, lep::LepAmmBaseLP const& lep, had::HadBaseLL const& had);
Real amm_base_lt1(Amm const& b, lep::LepAmmBaseLP const& lep, had::HadBaseLT const& had);
Real amm_base_lt2(Amm const& b, lep::LepAmmBaseLU const& lep, had::HadBaseLT const& had);
/// \}

/**
 * \defgroup XsNradGroup Base non-radiative cross-sections
 * Non-radiative cross-section base functions. They are used to compute
 * \f$\sigma_{\text{nrad}}^{IR}\f$, so they neglect the soft photon contribution
 * from \f$\sigma_{R}^{F}\f$.
 * \sa LepNradGroup
 * \sa HadGroup
 * \ingroup XsGroup
 */
/// \{

struct Nrad {
	Real coeff_born;
	Real coeff_amm;
	Nrad(kin::Kinematics const& kin, ph::Phenom const& phenom, Real k_0_bar);
};

Real nrad_ir_base_uu(Nrad const& b, lep::LepNradBaseUU const& lep, had::HadBaseUU const& had);
Real nrad_ir_base_ul(Nrad const& b, lep::LepNradBaseUP const& lep, had::HadBaseUL const& had);
Real nrad_ir_base_ut1(Nrad const& b, lep::LepNradBaseUP const& lep, had::HadBaseUT const& had);
Real nrad_ir_base_ut2(Nrad const& b, lep::LepNradBaseUU const& lep, had::HadBaseUT const& had);
Real nrad_ir_base_lu(Nrad const& b, lep::LepNradBaseLU const& lep, had::HadBaseLU const& had);
Real nrad_ir_base_ll(Nrad const& b, lep::LepNradBaseLP const& lep, had::HadBaseLL const& had);
Real nrad_ir_base_lt1(Nrad const& b, lep::LepNradBaseLP const& lep, had::HadBaseLT const& had);
Real nrad_ir_base_lt2(Nrad const& b, lep::LepNradBaseLU const& lep, had::HadBaseLT const& had);
/// \}

/**
 * \defgroup XsRadGroup Base radiative cross-sections
 * Radiative cross-section base functions. They do not allow the target
 * polarization part \f$\sigma^{XP}\f$ to be divided into pieces
 * \f$\sigma^{XL}\f$, \f$\sigma^{XT_1}\f$, and \f$\sigma^{XT_2}\f$.
 * \sa LepRadGroup
 * \sa HadRadGroup
 * \sa HadRadFGroup
 * \ingroup XsGroup
 */
/// \{

struct Rad {
	Real coeff;
	Real R;
	Rad(kin::KinematicsRad const& kin, ph::Phenom const& phenom);
};

Real rad_base_uu(Rad const& b, lep::LepRadBaseUU const& lep, had::HadRadBaseUU const& had);
math::Vec3 rad_base_up(Rad const& b, lep::LepRadBaseUU const& lep_uu, lep::LepRadBaseUP const& lep_up, had::HadRadBaseUP const& had);
Real rad_base_lu(Rad const& b, lep::LepRadBaseLU const& lep, had::HadRadBaseLU const& had);
math::Vec3 rad_base_lp(Rad const& b, lep::LepRadBaseLU const& lep_lu, lep::LepRadBaseLP const& lep_lp, had::HadRadBaseLP const& had);

Real rad_f_base_uu(Rad const& b, lep::LepRadBaseUU const& lep, had::HadRadFBaseUU const& had);
math::Vec3 rad_f_base_up(Rad const& b, lep::LepRadBaseUU const& lep_uu, lep::LepRadBaseUP const& lep_up, had::HadRadFBaseUP const& had);
Real rad_f_base_lu(Rad const& b, lep::LepRadBaseLU const& lep, had::HadRadFBaseLU const& had);
math::Vec3 rad_f_base_lp(Rad const& b, lep::LepRadBaseLU const& lep_lu, lep::LepRadBaseLP const& lep_lp, had::HadRadFBaseLP const& had);
/// \}

}
}

#endif

