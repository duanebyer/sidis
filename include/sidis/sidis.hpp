#ifndef SIDIS_SIDIS_HPP
#define SIDIS_SIDIS_HPP

#include "sidis/asymmetry.hpp"
#include "sidis/bound.hpp"
#include "sidis/constant.hpp"
#include "sidis/cross_section.hpp"
#include "sidis/cut.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/integ_params.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/numeric.hpp"
#include "sidis/particle.hpp"
#include "sidis/phenom.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/tmd.hpp"
#include "sidis/transform.hpp"
#include "sidis/vector.hpp"
#include "sidis/version.hpp"

namespace sidis {

/**
 * \mainpage
 *
 * The `sidis` library can be used to calculate the polarized SIDIS
 * cross-section with lowest order radiative corrections, following the work in
 * \cite akushevich2019sidis. To get started quickly see the
 * \ref ExampleSection "example" below. Then, read the documentation for the
 * xs::born(), xs::nrad_integ(), and xs::rad() functions. For a thorough
 * introduction to the `sidis` library, read this page.
 *
 * \section IntroductionSection Introduction to the SIDIS process
 *
 * The `sidis` library is used to compute cross-sections for the semi-inclusive
 * deep inelastic scattering (SIDIS) process, in which a lepton scatters from a
 * quark inside of a proton---destroying it---and the hadron formed from the
 * struck quark is detected. The SIDIS process is a useful probe of the quark
 * and gluon distributions inside of the proton. Since the `sidis` library
 * includes leading order radiative corrections to the SIDIS process, it
 * includes the following Feynman diagrams in its calculations:
 *
 * \image html sidis_rc_feynman_diagrams.svg "Leading order diagrams for radiative corrections to the SIDIS cross-section."
 * \image latex sidis_rc_feynman_diagrams.ps "Leading order diagrams for radiative corrections to the SIDIS cross-section."
 *
 * These diagrams give the following cross-sections:
 * * \f$\sigma_{B}\f$: The six-fold Born cross-section.
 * * \f$\sigma_{AMM}\f$: The six-fold anomalous magnetic moment (AMM)
 *   cross-section, related to the vertex correction diagram.
 * * \f$\sigma_{\text{vac}}\f$: The six-fold vacuum polarization cross-section.
 *   It has an infrared divergence.
 * * \f$\sigma_{R}\f$: The nine-fold radiative cross-section, corresponding to
 *   semi-inclusive diagrams in which a real photon is emitted. It also has an
 *   infrared divergence that cancels with the vacuum polarization divergence.
 * * \f$\sigma_{R,\text{ex}}\f$: The eight-fold exclusive part of the radiative
 *   cross-section. Currently, the `sidis` library does not support computation
 *   of this part of the cross-section.
 *
 * \subsection RcSection Radiative corrections
 *
 * To handle the infrared divergence, the radiative part of the cross-section is
 * first divided into an infrared-divergent part \f$\sigma_{R}^{IR}\f$ and an
 * infrared-divergent-free part \f$\sigma_{R}^{F}\f$ so that:
 *
 * \f{equation}{
 *     \sigma_{R} = \sigma_{R}^{IR} + \sigma_{R}^{F}
 *     \tag{1}
 * \f}
 *
 * A soft-photon cutoff \f$\bar{k}_{0,\text{cutoff}}\f$ is introduced. All
 * radiative events with photon energy below this cutoff (in a certain frame)
 * are combined with the Born, AMM, and vacuum polarization cross-sections to
 * give the total non-radiative cross-section.
 *
 * \f{equation}{
 *     \sigma_{\text{nrad}} = \sigma_{B} + \sigma_{AMM} + \sigma_{\text{vac}} + \int_0^{\bar{k}_{0,\text{cutoff}}} d\bar{k}_0 \, \int d\Omega_{\bar{k}} \, \sigma_{R}
 *     \tag{2}
 * \f}
 *
 * The integration over \f$\sigma_{R}\f$ can be done in two parts using
 * \f$\sigma_{R}^{IR}\f$ and \f$\sigma_{R}^{F}\f$. The integration over
 * \f$\sigma_{R}^{IR}\f$ can be evaluated analytically and then combined with
 * \f$\sigma_{\text{vac}}\f$ to give an infrared-divergent-free six-fold
 * cross-section \f$\sigma_{VR}\f$. The remaining integral over \f$\sigma_{R}^{F}\f$
 * must be evaluated numerically, but can be neglected if the soft-photon cutoff
 * is small enough (due to the absence of an infrared divergence).
 *
 * \f{eqnarray*}{
 *     \sigma_{\text{nrad}} &=& \sigma_{B} + \sigma_{AMM} + \sigma_{VR} + \int_0^{\bar{k}_{0,\text{cutoff}}} d\bar{k}_0 \, \int d\Omega_{\bar{k}} \, \sigma_{R}^{F} \newline
 *     &\approx& \sigma_{B} + \sigma_{AMM} + \sigma_{VR}
 *     \tag{3}
 * \f}
 *
 * Then, the SIDIS process has been reduced to two cross-sections: a six-fold
 * non-radiative cross-section \f$\sigma_{\text{nrad}}\f$, and a nine-fold
 * radiative cross-section \f$\sigma_{R}\f$ (which is now restricted to the
 * region of phase space with \f$\bar{k}_0 > \bar{k}_{0,\text{cutoff}}\f$).
 *
 * The expressions for these cross-sections are derived in detail in
 * \cite akushevich2019sidis.
 *
 * \subsection KinSection Kinematics
 *
 * The non-radiative SIDIS process, with target \f$n\f$, detected hadron
 * \f$h\f$, and undetected fragments \f$x\f$ is:
 *
 * \f{equation}{
 *     e(k_1) + n(p) \rightarrow e(k_2) + h(p_h) + x(p_x)
 *     \tag{4}
 * \f}
 *
 * \image html sidis_kinematics_diagram.svg "SIDIS kinematic diagram."
 * \image latex sidis_kinematics_diagram.ps "SIDIS kinematic diagram."
 *
 * The mass of the fragments \f$m_x\f$ is allowed to vary (although constrained
 * to be larger than the threshold mass \f$M_{\text{th}}\f$), leaving six
 * degrees of freedom. The `sidis` library chooses them as follows:
 *
 * \f{equation}{
 *     x = -\frac{q^2}{2 q p},\qquad y = \frac{q p}{k_1 p},\qquad z = \frac{p_h p}{p q},\qquad p_{t}^2 = \frac{|\pmb{p}_h\times\pmb{q}|^2}{|\pmb{q}|^2},\qquad \phi_h,\qquad \phi
 *     \tag{5}
 * \f}
 *
 * Where \f$q = k_1 - k_2\f$, \f$\phi_h\f$ is the azimuthal angle of the hadron
 * about \f$\pmb{q}\f$, and \f$\phi\f$ is the azimuthal angle of the scattered
 * lepton about the target polarization direction. (Note that in the deep
 * inelastic regime, \f$\phi\approx\phi_S\f$). For the radiative case, when
 * a photon is emitted with 4-momentum \f$k\f$, three additional kinematic
 * variables are used.
 *
 * \f{equation}{
 *     R=2 k p,\qquad \tau=\frac{k q}{k p},\qquad \phi_k
 *     \tag{6}
 * \f}
 *
 * Where \f$\phi_k\f$ is the azimuthal angle of the photon about \f$\pmb{q}\f$.
 *
 * In the exclusive radiative case, \f$m_x\f$ is fixed to a specific value,
 * resulting in only eight degrees of freedom. This case is not handled by the
 * `sidis` library.
 *
 * \section ComputingSection Computing the cross-section
 *
 * The details of the cross-section computation are not necessary to use the
 * `sidis` library, although they can help with getting the most possible
 * performance. At the most general level, the cross-sections are calculated
 * using the process in \cite akushevich2019sidis.
 *
 * \f{equation}{
 *     \sigma = C\sum_{i} \theta_{i}\mathcal{H}_i
 *     \tag{7}
 * \f}
 *
 * Where \f$\theta_{i}\f$ are coefficients describing the lepton part of the
 * interaction, and \f$\mathcal{H}_i\f$ describe the hadron part of the
 * interaction (using the standard SIDIS structure functions).
 *
 * The process by which the cross sections are calculated by the `sidis` library
 * is:
 * * Set up the initial state, including particle IDs, the beam energy, and the
 *   threshold mass \f$M_{\text{th}}\f$ (see part::Particles and kin::Initial).
 * * Choose a point in phase space using the kinematic variables from equations
 *   (5) or (6) (see kin::PhaseSpace and kin::PhaseSpaceRad).
 * * Calculate various kinematic quantities and store the results for later use
 *   (see kin::PhaseSpace and kin::Kinematics).
 * * Calculate the "leptonic coefficients". These coefficients are related to
 *   the \f$\theta_i\f$ coefficients(see \ref LepGroup "lepton coefficients").
 * * Calculate the structure functions (see sf::SfSet).
 * * Use the structure functions to calculate the "hadronic coefficients". These
 *   coefficients are related to the \f$\mathcal{H}_i\f$ (see
 *   \ref HadGroup "hadron coefficients").
 * * Combine the leptonic and hadronic coefficients into a cross-section using
 *   equation (7) (see \ref XsGroup "cross-section" namespace).
 *
 * \subsection PolSection Polarization notation
 *
 * For each different target and beam polarization, and for each type of
 * cross-section, there is a different subset of lepton and hadron coefficients
 * that must be computed. For example, the lep::LepAmmBaseUU coefficients are
 * used when computing the AMM cross-section for unpolarized beam and target.
 * Any polarized cross-section for a lepton beam with longitudinal polarization
 * \f$\lambda_e\f$ and a target with polarization \f$\pmb{\eta}\f$ can be
 * decomposed as:
 *
 * \f{equation}{
 *     \sigma = \sigma^{UU} + \pmb{\eta}\cdot\pmb{\sigma}^{UP} + \lambda_e\left(\sigma^{LU} + \pmb{\eta}\cdot\pmb{\sigma}^{LP}\right)
 *     \tag{8}
 * \f}
 *
 * Where:
 *
 * \f{equation}{
 *     \pmb{\sigma}^{XP} = \pmb{e}_x\sigma^{XT_1} + \pmb{e}_y\sigma^{XT_2} + \pmb{e}_z\sigma^{XL}
 *     \tag{9}
 * \f}
 *
 * The following notation is used by the `sidis` library for polarizations:
 * * Beam
 *   * U: Unpolarized
 *   * L: Longitudinally polarized
 * * Target
 *   * U: Unpolarized
 *   * L: Longitudinally polarized
 *   * T1: Transversely polarized in the lepton plane
 *   * T2: Transversely polarized out of the lepton plane
 *   * T: Transversely polarized
 *   * P: Polarized in any direction
 *
 * For example, the had::HadBaseLP structure is used to calculate the
 * \f$\pmb{\sigma}^{LP}\f$ cross-section, which corresponds to a longitudinally
 * polarized beam and a target polarized in any direction.
 *
 * \subsection CalcKinSection Kinematics
 *
 * The kinematics are determined by the six SIDIS variables. The results of the
 * full kinematic calculations are stored in the kin::Kinematics structure. The
 * six SIDIS variables have kinematic bounds on what values they can take. These
 * bounds can be determined using methods from the \ref CutGroup namespace, such
 * as cut::x_bound(). It may also be useful to apply cuts to various kinematic
 * variables (e.g. when generating events). The cut::Cut and cut::CutRad
 * structures can be used for this purpose.
 *
 * \subsection CalcLepCoeffsSection Lepton coefficients
 *
 * The lepton coefficients are specialized for every type of cross-section
 * computation: lep::LepBorn for Born, lep::LepAmm for AMM, lep::LepNrad for
 * non-radiative (including Born, AMM, and vacuum polarization), and lep::LepRad
 * for radiative. See the \ref LepGroup "lepton coefficients" for more
 * information.
 *
 * \subsection CalcHadCoeffsSection Hadron coefficients
 *
 * There are three types of hadron coefficients: had::Had for non-radiative
 * cross-sections, had::HadRad for radiative cross-sections, and had::HadRadF
 * for the infrared-divergent-free part of the radiative cross-sections
 * \f$\sigma_{R}^{F}\f$ (from equation (1)). See the
 * \ref HadGroup "hadron coefficients" for more information.
 *
 * \subsection CalcXsSection Cross-sections
 *
 * As per equation (8), the cross-sections can be calculated in a number of
 * parts. For example, the `xs::born_base_ut1()` function computes the term
 * \f$\sigma_{B}^{UT_1}\f$. For ease of use, the xs::born() function can combine
 * all of these pieces with the provided polarizations to give a cross-section
 * (and similarly with the other types of cross-sections). See the
 * \ref XsGroup "cross-section page" for more information.
 *
 * \section ExampleSection Example
 *
 * The following example shows how to easily get started using the `sidis`
 * library to compute cross-sections.
 *
 * \include quick_start.cpp
 *
 * To extend this example, try replacing xs::born() with another kind of
 * cross-section, such as xs::amm(), xs::nrad_integ(), or even xs::rad(). You
 * can also experiment with the integrated radiative cross-sections
 * xs::rad_integ() and xs::rad_f_integ().
 *
 * \example quick_start.cpp
 * Demonstrates how to compute a cross-section with the `sidis` library.
 *
 * \example cross_section.cpp
 * A simple utility that computes the cross-section with and without radiative
 * corrections at the provided point in phase space.
 */

}

#endif

