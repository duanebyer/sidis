#ifndef SIDIS_CROSS_SECTION_HPP
#define SIDIS_CROSS_SECTION_HPP

#include "sidis/constant.hpp"
#include "sidis/integ.hpp"
#include "sidis/numeric.hpp"

namespace sidis {

namespace math {
	struct Vec3;
}
namespace kin {
	struct Kinematics;
	struct KinematicsRad;
}
namespace lep {
	struct LepBornUU;
	struct LepBornUP;
	struct LepBornUX;
	struct LepBornLU;
	struct LepBornLP;
	struct LepBornLX;
	struct LepBornXX;

	struct LepAmmUU;
	struct LepAmmUP;
	struct LepAmmUX;
	struct LepAmmLU;
	struct LepAmmLP;
	struct LepAmmLX;
	struct LepAmmXX;

	struct LepNradUU;
	struct LepNradUP;
	struct LepNradUX;
	struct LepNradLU;
	struct LepNradLP;
	struct LepNradLX;
	struct LepNradXX;

	struct LepRadUU;
	struct LepRadUP;
	struct LepRadUX;
	struct LepRadLU;
	struct LepRadLP;
	struct LepRadLX;
	struct LepRadXX;
}
namespace had {
	struct HadUU;
	struct HadUL;
	struct HadUT;
	struct HadUP;
	struct HadLU;
	struct HadLL;
	struct HadLT;
	struct HadLP;
	struct HadXX;

	struct HadRadUU;
	struct HadRadUP;
	struct HadRadLU;
	struct HadRadLP;
	struct HadRadXX;

	struct HadRadFUU;
	struct HadRadFUP;
	struct HadRadFLU;
	struct HadRadFLP;
	struct HadRadFXX;
}
namespace sf {
	struct SfSet;
}

namespace xs {

math::IntegParams const DEFAULT_INTEG_PARAMS { 1000000, 1e-6, 0. };

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
 * target polarizations.
 * \ingroup XsGroup
 */
/// \{

/// %Born cross-section \f$\sigma_{B}\f$.
Real born(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// Anomalous magnetic moment cross-section \f$\sigma_{AMM}\f$, related to
/// vertex correction diagram.
Real amm(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// Non-radiative cross-section neglecting infrared-divergent-free soft photon
/// part \f$\sigma_{\text{nrad}}^{IR}\f$.
Real nrad_ir(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF);
/// Non-radiative cross-section \f$\sigma_{\text{nrad}}\f$, integrated over the
/// radiated photon with energy below soft cutoff \p k_0_bar (if \p k_0_bar is
/// set to infinity (default), then the entire radiative part is integrated
/// over).
math::EstErr nrad_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// Radiative cross-section with infrared divergence removed
/// \f$\sigma_{R}^{F}\f$.
Real rad_f(kin::KinematicsRad const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
/// Radiative cross-section \f$\sigma_{R}\f$.
Real rad(kin::KinematicsRad const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);

/// Radiative cross-section with infrared divergence removed
/// \f$\sigma_{R}^{F}\f$, integrated over the radiated photon with energy above
/// soft cutoff \p k_0_bar.
math::EstErr rad_f_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
/// Radiative cross-section \f$\sigma_{R}\f$, integrated over the radiated
/// photon with energy above soft cutoff \p k_0_bar.
math::EstErr rad_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, math::IntegParams params=DEFAULT_INTEG_PARAMS);
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
	explicit Born(kin::Kinematics const& kin);
};

Real born_uu_base(Born const& b, lep::LepBornUU const& lep, had::HadUU const& had);
Real born_ul_base(Born const& b, lep::LepBornUP const& lep, had::HadUL const& had);
Real born_ut1_base(Born const& b, lep::LepBornUP const& lep, had::HadUT const& had);
Real born_ut2_base(Born const& b, lep::LepBornUU const& lep, had::HadUT const& had);
math::Vec3 born_up_base(Born const& b, lep::LepBornUX const& lep, had::HadUP const& had);
Real born_lu_base(Born const& b, lep::LepBornLU const& lep, had::HadLU const& had);
Real born_ll_base(Born const& b, lep::LepBornLP const& lep, had::HadLL const& had);
Real born_lt1_base(Born const& b, lep::LepBornLP const& lep, had::HadLT const& had);
Real born_lt2_base(Born const& b, lep::LepBornLU const& lep, had::HadLT const& had);
math::Vec3 born_lp_base(Born const& b, lep::LepBornLX const& lep, had::HadLP const& had);
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
	explicit Amm(kin::Kinematics const& kin);
};

Real amm_uu_base(Amm const& b, lep::LepAmmUU const& lep, had::HadUU const& had);
Real amm_ul_base(Amm const& b, lep::LepAmmUP const& lep, had::HadUL const& had);
Real amm_ut1_base(Amm const& b, lep::LepAmmUP const& lep, had::HadUT const& had);
Real amm_ut2_base(Amm const& b, lep::LepAmmUU const& lep, had::HadUT const& had);
math::Vec3 amm_up_base(Amm const& b, lep::LepAmmUX const& lep, had::HadUP const& had);
Real amm_lu_base(Amm const& b, lep::LepAmmLU const& lep, had::HadLU const& had);
Real amm_ll_base(Amm const& b, lep::LepAmmLP const& lep, had::HadLL const& had);
Real amm_lt1_base(Amm const& b, lep::LepAmmLP const& lep, had::HadLT const& had);
Real amm_lt2_base(Amm const& b, lep::LepAmmLU const& lep, had::HadLT const& had);
math::Vec3 amm_lp_base(Amm const& b, lep::LepAmmLX const& lep, had::HadLP const& had);
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
	explicit Nrad(kin::Kinematics const& kin, Real k_0_bar);
};

Real nrad_ir_uu_base(Nrad const& b, lep::LepNradUU const& lep, had::HadUU const& had);
Real nrad_ir_ul_base(Nrad const& b, lep::LepNradUP const& lep, had::HadUL const& had);
Real nrad_ir_ut1_base(Nrad const& b, lep::LepNradUP const& lep, had::HadUT const& had);
Real nrad_ir_ut2_base(Nrad const& b, lep::LepNradUU const& lep, had::HadUT const& had);
math::Vec3 nrad_ir_up_base(Nrad const& b, lep::LepNradUX const& lep, had::HadUP const& had);
Real nrad_ir_lu_base(Nrad const& b, lep::LepNradLU const& lep, had::HadLU const& had);
Real nrad_ir_ll_base(Nrad const& b, lep::LepNradLP const& lep, had::HadLL const& had);
Real nrad_ir_lt1_base(Nrad const& b, lep::LepNradLP const& lep, had::HadLT const& had);
Real nrad_ir_lt2_base(Nrad const& b, lep::LepNradLU const& lep, had::HadLT const& had);
math::Vec3 nrad_ir_lp_base(Nrad const& b, lep::LepNradLX const& lep, had::HadLP const& had);
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
	explicit Rad(kin::KinematicsRad const& kin);
};

Real rad_uu_base(Rad const& b, lep::LepRadUU const& lep, had::HadRadUU const& had);
math::Vec3 rad_up_base(Rad const& b, lep::LepRadUX const& lep, had::HadRadUP const& had);
Real rad_lu_base(Rad const& b, lep::LepRadLU const& lep, had::HadRadLU const& had);
math::Vec3 rad_lp_base(Rad const& b, lep::LepRadLX const& lep, had::HadRadLP const& had);

Real rad_f_uu_base(Rad const& b, lep::LepRadUU const& lep, had::HadRadFUU const& had);
math::Vec3 rad_f_up_base(Rad const& b, lep::LepRadUX const& lep, had::HadRadFUP const& had);
Real rad_f_lu_base(Rad const& b, lep::LepRadLU const& lep, had::HadRadFLU const& had);
math::Vec3 rad_f_lp_base(Rad const& b, lep::LepRadLX const& lep, had::HadRadFLP const& had);
/// \}

}
}

#endif

