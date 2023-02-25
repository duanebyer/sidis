#ifndef SIDIS_KINEMATICS_HPP
#define SIDIS_KINEMATICS_HPP

#include "sidis/numeric.hpp"
#include "sidis/particle.hpp"
#include "sidis/vector.hpp"

namespace sidis {
namespace kin {

/**
 * \defgroup KinGroup Kinematics info
 * Types and functions related to kinematics calculations.
 */
/// \{

/**
 * Point in phase space. These six kinematic variables uniquely pick out a point
 * in the SIDIS phase space.
 */
struct PhaseSpace {
	/// Bjorken variable \f$x = -\frac{q^2}{2 q p}\f$.
	Real x;
	/// Fraction of incident lepton energy transferred to virtual photon
	/// \f$y = \frac{q p}{k_1 p}\f$.
	Real y;
	/// Fraction of virtual photon energy transferred to hadron
	/// \f$z = \frac{p_h p}{p q}\f$.
	Real z;
	/// Detected hadron momentum transverse to virtual photon.
	Real ph_t_sq;
	/// Azimuthal angle of hadron about virtual photon.
	Real phi_h;
	/// Azimuthal angle of scattered lepton about target polarization.
	Real phi;
};

/**
 * Point in radiative phase space. These nine kinematic variables uniquely pick
 * out a point in the radiative SIDIS phase space.
 */
struct PhaseSpaceRad {
	/// \copydoc PhaseSpace::x
	Real x;
	/// \copydoc PhaseSpace::y
	Real y;
	/// \copydoc PhaseSpace::z
	Real z;
	/// \copydoc PhaseSpace::ph_t_sq
	Real ph_t_sq;
	/// \copydoc PhaseSpace::phi_h
	Real phi_h;
	/// \copydoc PhaseSpace::phi
	Real phi;
	/// Ratio of virtual photon and radiated photon energies
	/// \f$\tau=\frac{kq}{kp}\f$.
	Real tau;
	/// Azimuthal angle of radiated photon about virtual photon.
	Real phi_k;
	/// Related to energy of radiated photon \f$R = 2kp\f$.
	Real R;

	/// Discard the radiative kinematic variables \f$\tau\f$, \f$\phi_k\f$, and
	/// \f$R\f$ to get a point in the non-radiative PhaseSpace.
	PhaseSpace project() const {
		return { x, y, z, ph_t_sq, phi_h, phi };
	}
};

/**
 * Describes the kinematics of a SIDIS process. This structure is used to cache
 * various kinematic quantities so they can be re-used throughout the
 * cross-section calculations.
 */
struct Kinematics {
	/// \name %Initial state description
	/// \{

	/// \copydoc part::Particles::target
	part::Nucleus target;
	/// \copydoc part::Particles::beam
	part::Lepton beam;
	/// \copydoc part::Particles::hadron
	part::Hadron hadron;

	/// Related to the beam energy. Given by \f$S = 2 k_1 p\f$.
	Real S;
	/// \copydoc part::Particles::M
	Real M;
	/// \copydoc part::Particles::m
	Real m;
	/// \copydoc part::Particles::mh
	Real mh;
	/// \copydoc part::Particles::Mth
	Real Mth;
	/// \}

	/// \name Base kinematic variables
	/// A minimum set of kinematic variables needed to describe the SIDIS
	/// process.
	/// \sa PhaseSpace
	/// \{

	/// \copydoc PhaseSpace::x
	Real x;
	/// \copydoc PhaseSpace::y
	Real y;
	/// \copydoc PhaseSpace::z
	Real z;
	/// \copydoc PhaseSpace::ph_t_sq
	Real ph_t_sq;
	/// \copydoc PhaseSpace::phi_h
	Real phi_h;
	/// \copydoc PhaseSpace::phi
	Real phi;
	/// \}

	/// \name Additional kinematic variables.
	/// Additional kinematic variables that are needed for cross-section
	/// calculations.
	/// \{

	/// Squared mass of the virtual photon (which is always space-like).
	Real Q_sq;
	/// Mass of the virtual photon.
	Real Q;
	/// Momentum transfer from virtual photon into hadron \f$t = (q - p_h)^2\f$.
	Real t;
	/// Center-of-mass energy of the virtual photon and hadron system
	/// \f$W^2 = (p + q)^2\f$.
	Real W_sq;
	/// Related to energy of scattered lepton \f$X = 2 p k_2\f$.
	Real X;
	/// Related to energy of virtual photon \f$S_x = 2 p q\f$.
	Real S_x;
	/// Defined as \f$S_p = 2 p (k_1 + k_2)\f$.
	Real S_p;
	/// Defined as \f$V_1 = 2 k_1 p_h \f$.
	Real V_1;
	/// Defined as \f$V_2 = 2 k_2 p_h \f$.
	Real V_2;
	/// Defined as \f$V_{-} = 2 q p_h\f$.
	Real V_m;
	/// Defined as \f$V_{+} = 2 (k_1 + k_2) p_h\f$.
	Real V_p;

	/// Defined as \f$\lambda_S = 4 M^2 |\pmb{k}_1|^2\f$.
	Real lambda_S;
	/// Defined as \f$\lambda_X = 4 M^2 |\pmb{k}_2|^2\f$.
	Real lambda_X;
	/// Defined as \f$\lambda_Y = 4 M^2 |\pmb{q}|^2\f$.
	Real lambda_Y;
	/// Defined as \f$\lambda_1 = 4 M^2 q_t^2|\pmb{k_1}|^2\f$.
	Real lambda_1;
	/// Defined as \f$\lambda_2 = V_{-}^2 + m_h^2 Q^2\f$.
	Real lambda_2;
	/// Defined as \f$\lambda_3 = V_{-}^2 + z Q^2\f$.
	Real lambda_3;

	/// Cached \f$\sqrt{\lambda_S}\f$.
	Real lambda_S_sqrt;
	/// Cached \f$\sqrt{\lambda_X}\f$.
	Real lambda_X_sqrt;
	/// Cached \f$\sqrt{\lambda_Y}\f$.
	Real lambda_Y_sqrt;
	/// Cached \f$\sqrt{\lambda_1}\f$.
	Real lambda_1_sqrt;
	/// Cached \f$\sqrt{\lambda_2}\f$.
	Real lambda_2_sqrt;
	/// Cached \f$\sqrt{\lambda_3}\f$.
	Real lambda_3_sqrt;

	/// Squared mass of the undetected fragments \f$m_x^2 = p_x^2\f$.
	Real mx_sq;
	/// Mass of the undetected fragments.
	Real mx;
	/// Volume element associated with the azimthal angle \f$\phi_h\f$. Given by
	/// \f$\varepsilon^{\mu\nu\rho\sigma}p_{h\mu}p_{\nu}k_{1\rho}q_{\sigma}\f$.
	Real vol_phi_h;

	/// Useful coefficient for structure function calculations.
	Real C_1;
	/// \}

	/// \name 4-momentum components
	/// Energy, and momentum components of the particles in various frames. For
	/// the definitions of the frames themselves, see the
	/// \ref FrameGroup "frame" namespace.
	/// \{

	/// Hadron energy in the lepton frame.
	Real ph_0;
	/// Hadron transverse momentum in the lepton frame.
	Real ph_t;
	/// Hadron longitudinal momentum in the lepton frame.
	Real ph_l;
	/// Cached \f$\cos\phi_h\f$.
	Real cos_phi_h;
	/// Cached \f$\sin\phi_h\f$.
	Real sin_phi_h;

	/// Virtual photon energy in the target frame.
	Real q_0;
	/// Virtual photon tranverse momentum in the target frame.
	Real q_t;
	/// Virtual photon longitudinal momentum in the target frame.
	Real q_l;
	/// The azimuthal angle of the virtual photon about the incoming lepton.
	Real phi_q;
	/// Cached \f$\cos\phi_q\f$.
	Real cos_phi_q;
	/// Cached \f$\sin\phi_q\f$.
	Real sin_phi_q;

	// Components of scattered lepton in target frame.
	/// Scattered lepton energy in the target frame.
	Real k2_0;
	/// Scattered lepton transverse momentum in the target frame.
	Real k2_t;
	/// Scattered lepton longitudinal momentum in the target frame.
	Real k2_l;

	/// Incident lepton transverse momentum in the target frame.
	Real k1_t;
	/// Cached \f$\cos\phi\f$.
	Real cos_phi;
	/// Cached \f$\sin\phi\f$.
	Real sin_phi;
	/// \}

	/// Initialize an empty Kinematics in an invalid state.
	Kinematics() = default;
	/// Fill in a Kinematics corresponding to particles \p ps, with beam energy
	/// given by \p S (with \f$S = 2 p k_1\f$), and at a PhaseSpace point
	/// \p ph_space.
	Kinematics(part::Particles const& ps, Real S, PhaseSpace const& ph_space);
};

/**
 * Describes the kinematics of a radiative SIDIS process. This structure is used
 * to cache various kinematic quantities so they can be re-used throughout the
 * cross-section calculations.
 */
struct KinematicsRad {
	/// \name %Initial state description
	/// \{
	/// \copydoc Kinematics::target
	part::Nucleus target;
	/// \copydoc Kinematics::beam
	part::Lepton beam;
	/// \copydoc Kinematics::hadron
	part::Hadron hadron;

	/// \copydoc Kinematics::S
	Real S;
	/// \copydoc Kinematics::M
	Real M;
	/// \copydoc Kinematics::m
	Real m;
	/// \copydoc Kinematics::mh
	Real mh;
	/// \copydoc Kinematics::Mth
	Real Mth;
	/// \}

	/// \name Base non-radiative kinematic variables
	/// A minimum set of kinematic variables needed to describe the SIDIS
	/// process.
	/// \sa PhaseSpace
	/// \{

	/// \copydoc Kinematics::x
	Real x;
	/// \copydoc Kinematics::y
	Real y;
	/// \copydoc Kinematics::z
	Real z;
	/// \copydoc Kinematics::ph_t_sq
	Real ph_t_sq;
	/// \copydoc Kinematics::phi_h
	Real phi_h;
	/// \copydoc Kinematics::phi
	Real phi;
	/// \}

	/// \name Additional non-radiative kinematic variables.
	/// Additional kinematic variables that are needed for cross-section
	/// calculations.
	/// \{

	/// \copydoc Kinematics::Q_sq
	Real Q_sq;
	/// \copydoc Kinematics::Q
	Real Q;
	/// \copydoc Kinematics::t
	Real t;
	/// \copydoc Kinematics::W_sq
	Real W_sq;
	/// \copydoc Kinematics::X
	Real X;
	/// \copydoc Kinematics::S_x
	Real S_x;
	/// \copydoc Kinematics::S_p
	Real S_p;
	/// \copydoc Kinematics::V_1
	Real V_1;
	/// \copydoc Kinematics::V_2
	Real V_2;
	/// \copydoc Kinematics::V_m
	Real V_m;
	/// \copydoc Kinematics::V_p
	Real V_p;

	/// \copydoc Kinematics::lambda_S
	Real lambda_S;
	/// \copydoc Kinematics::lambda_X
	Real lambda_X;
	/// \copydoc Kinematics::lambda_Y
	Real lambda_Y;
	/// \copydoc Kinematics::lambda_1
	Real lambda_1;
	/// \copydoc Kinematics::lambda_2
	Real lambda_2;
	/// \copydoc Kinematics::lambda_3
	Real lambda_3;
	/// \copydoc Kinematics::lambda_S_sqrt
	Real lambda_S_sqrt;
	/// \copydoc Kinematics::lambda_X_sqrt
	Real lambda_X_sqrt;
	/// \copydoc Kinematics::lambda_Y_sqrt
	Real lambda_Y_sqrt;
	/// \copydoc Kinematics::lambda_1_sqrt
	Real lambda_1_sqrt;
	/// \copydoc Kinematics::lambda_2_sqrt
	Real lambda_2_sqrt;
	/// \copydoc Kinematics::lambda_3_sqrt
	Real lambda_3_sqrt;

	/// \copydoc Kinematics::mx_sq
	Real mx_sq;
	/// \copydoc Kinematics::mx
	Real mx;
	/// \copydoc Kinematics::vol_phi_h
	Real vol_phi_h;

	/// \copydoc Kinematics::C_1
	Real C_1;
	/// \}

	/// \name Non-radiative 4-momentum components
	/// Energy, and momentum components of the particles in various frames. For
	/// the definitions of the frames themselves, see the
	/// \ref FrameGroup "frame" namespace.
	/// \{

	/// \copydoc Kinematics::ph_0
	Real ph_0;
	/// \copydoc Kinematics::ph_t
	Real ph_t;
	/// \copydoc Kinematics::ph_l
	Real ph_l;
	/// \copydoc Kinematics::cos_phi_h
	Real cos_phi_h;
	/// \copydoc Kinematics::sin_phi_h
	Real sin_phi_h;
	/// \copydoc Kinematics::q_0
	Real q_0;
	/// \copydoc Kinematics::q_t
	Real q_t;
	/// \copydoc Kinematics::q_l
	Real q_l;
	/// \copydoc Kinematics::phi_q
	Real phi_q;
	/// \copydoc Kinematics::cos_phi_q
	Real cos_phi_q;
	/// \copydoc Kinematics::sin_phi_q
	Real sin_phi_q;
	/// \copydoc Kinematics::k2_0
	Real k2_0;
	/// \copydoc Kinematics::k2_t
	Real k2_t;
	/// \copydoc Kinematics::k2_l
	Real k2_l;
	/// \copydoc Kinematics::k1_t
	Real k1_t;
	/// \copydoc Kinematics::cos_phi
	Real cos_phi;
	/// \copydoc Kinematics::sin_phi
	Real sin_phi;
	/// \}

	/// \name Base radiative kinematic variables
	/// A minimum set of kinematic variables needed to describe the radiative
	/// SIDIS process.
	/// \sa PhaseSpaceRad
	/// \{

	/// \copydoc PhaseSpaceRad::tau
	Real tau;
	/// \copydoc PhaseSpaceRad::phi_k
	Real phi_k;
	/// \copydoc PhaseSpaceRad::R
	Real R;
	/// \}

	/// \name Additional radiative kinematic variables.
	/// Additional kinematic variables that are needed for cross-section
	/// calculations.
	/// \{

	/// The lower kinematic bound on \f$\tau\f$.
	Real tau_min;
	/// The upper kinematic bound on \f$\tau\f$.
	Real tau_max;
	/// The ratio \f$\frac{k p_h}{k p}\f$.
	Real mu;
	/// The ratio \f$\frac{k_1 k}{p k}\f$.
	Real z_1;
	/// The ratio \f$\frac{k_2 k}{p k}\f$.
	Real z_2;

	/// Defined as \f$\lambda_V = 4 M^2 \pmb{p}_h\cdot\pmb{q}\f$.
	Real lambda_V;
	/// Defined as \f$\lambda_{RV} = 4 M^2 \pmb{k}\cdot\pmb{p}_h\f$.
	Real lambda_RV;
	/// Defined as \f$\lambda_{RY} = 4 M^2 \pmb{k}\cdot\pmb{q}\f$.
	Real lambda_RY;
	/// Defined as \f$\lambda_H = 4 M^2 |\pmb{p}_h|^2\f$.
	Real lambda_H;
	/// Defined as \f$\lambda_z = (\tau_{\text{min}} - \tau)(\tau - \tau_{\text{max}})\lambda_1\f$.
	Real lambda_z;
	/// Cached \f$\sqrt{\lambda_z}\f$.
	Real lambda_z_sqrt;

	/// The radiated photon energy in the frame \f$\pmb{p}_h = 0\f$.
	Real k_0_bar;

	/// Volume element associated with the azimthal angle \f$\phi_k\f$. Given by
	/// \f$\varepsilon^{\mu\nu\rho\sigma}k_{h\mu}p_{\nu}k_{1\rho}q_{\sigma}/R\f$.
	Real vol_phi_k_R;
	/// Volume element. Given by
	/// \f$\varepsilon^{\mu\nu\rho\sigma}k_{\mu}p_{\nu}p_{h\rho}q_{\sigma}\f$
	Real vol_phi_hk;

	/// Defined as \f$F_{22} = \frac{1}{z_2^2}\f$.
	Real F_22;
	/// Defined as \f$F_{21} = \frac{1}{z_1^2}\f$.
	Real F_21;
	/// Defined as \f$F_{2+} = \frac{1}{z_2^2} + \frac{1}{z_1^2}\f$.
	Real F_2p;
	/// Defined as \f$F_{2-} = \frac{1}{z_2^2} - \frac{1}{z_1^2}\f$.
	Real F_2m;
	/// Defined as \f$F_d = \frac{1}{z_1 z_2}\f$.
	Real F_d;
	/// Defined as \f$F_{1+} = \frac{1}{z_1} + \frac{1}{z_2}\f$.
	Real F_1p;
	/// Defined as \f$F_{IR} = m^2 F_{2+} - (Q^2 + 2 m^2)F_d\f$.
	Real F_IR;
	/// \}

	/// \name Radiative 4-momentum components
	/// Energy, and momentum components of the radiated photon in the lepton
	/// frame. For the definition of the lepton frame, see the
	/// \ref FrameGroup "frame" namespace.
	/// \{

	// Components of real photon 4-momentum in lepton frame.
	/// Radiated photon energy in the lepton frame.
	Real k_0;
	/// Radiated photon transverse momentum in the lepton frame.
	Real k_t;
	/// Radiated photon longitudinal momentum in the lepton frame.
	Real k_l;
	/// Cached \f$\cos\phi_k\f$.
	Real cos_phi_k;
	/// Cached \f$\sin\phi_k\f$.
	Real sin_phi_k;
	/// \}

	/// \name Shifted kinematic variables
	/// The following kinematic variables are related to their unshifted
	/// counterparts in the following way: A shifted variable is related to its
	/// unshifted counterpart through the substitution \f$q\rightarrow q - k\f$.
	/// For example, for the Bjorken variable \f$x\f$:
	///
	/// \f{eqnarray*}{
	///     x &=& -\frac{q^2}{2 q p} = \frac{Q^2}{S_x} \newline
	///     \tilde{x} &=& -\frac{(q - k)^2}{2 (q - k) p}
	///         = -\frac{q^2 - 2 q k}{2 q p - 2 k p}
	///         = \frac{Q^2 + R\tau}{S_x - R}
	///         = \frac{\tilde{Q}^2}{\tilde{S}_x}
	/// \f}
	/// \{

	/// Shifted kinematic variable.
	Real shift_x;
	Real shift_y;
	Real shift_z;
	Real shift_ph_t_sq;
	Real shift_phi_h;
	Real shift_phi_q;

	Real shift_cos_phi_h;
	Real shift_sin_phi_h;
	Real shift_cos_phi_q;
	Real shift_sin_phi_q;

	Real shift_Q_sq;
	Real shift_Q;
	Real shift_t;
	Real shift_W_sq;
	Real shift_S_x;
	Real shift_V_m;

	Real shift_lambda_Y;
	Real shift_lambda_1;
	Real shift_lambda_2;
	Real shift_lambda_3;
	Real shift_lambda_Y_sqrt;
	Real shift_lambda_1_sqrt;
	Real shift_lambda_2_sqrt;
	Real shift_lambda_3_sqrt;

	Real shift_mx_sq;
	Real shift_mx;
	Real shift_vol_phi_h;

	Real shift_C_1;

	Real shift_ph_t;
	Real shift_ph_l;
	Real shift_q_0;
	Real shift_q_t;
	Real shift_q_l;
	Real shift_k1_t;
	/// \}

	/// Discard the radiative kinematic variables to get a Kinematics describing
	/// a point in the non-radiative PhaseSpace.
	Kinematics project() const;
	/// Use the shifted kinematic variables to get a Kinematics describing a
	/// shifted point in the non-radiative PhaseSpace.
	Kinematics project_shift() const;

	/// Initialize an empty KinematicsRad in an invalid state.
	KinematicsRad() = default;
	/// Fill in a KinematicsRad corresponding to particles \p ps, with beam
	/// energy given by \p S (with \f$S = 2 p k_1\f$), and at a PhaseSpaceRad
	/// point \p ph_space.
	KinematicsRad(part::Particles const& ps, Real S, PhaseSpaceRad const& ph_space) :
		KinematicsRad(
			Kinematics(
				ps,
				S,
				{
					ph_space.x,
					ph_space.y,
					ph_space.z,
					ph_space.ph_t_sq,
					ph_space.phi_h,
					ph_space.phi,
				}),
			ph_space.tau,
			ph_space.phi_k,
			ph_space.R) { }
	/// Take an existing non-radiative Kinematics, and add the additional base
	/// radiative variables \f$\tau\f$, \f$\phi_k\f$, and \f$R\f$, to form a
	/// radiative KinematicsRad.
	KinematicsRad(Kinematics const& kin, Real tau, Real phi_k, Real R);
};

/**
 * %Initial state of the system before the SIDIS process. Contains the target
 * and beam particle types as well as the initial 4-momenta of both.
 */
struct Initial {
	/// \copydoc part::Particles::target
	part::Nucleus const target;
	/// \copydoc part::Particles::beam
	part::Lepton const beam;
	/// The initial 4-momentum of the target.
	math::Vec4 const p;
	/// The initial 4-momentum of the beam.
	math::Vec4 const k1;

	/// Prepare the initial state with particles \p ps, initial target momentum
	/// \p p, and initial beam momentum \p k1.
	Initial(
		part::Particles const& ps, math::Vec3 p, math::Vec3 k1) :
		target(ps.target),
		beam(ps.beam),
		p(math::Vec4::from_length_and_r(ps.M, p)),
		k1(math::Vec4::from_length_and_r(ps.m, k1)) { }

	/// Prepare the initial state with particles \p ps, with the target at rest,
	/// and with the beam pointed along \f$\pmb{e}_z\f$ with energy
	/// \p beam_energy.
	Initial(part::Particles const& ps, Real beam_energy) :
		target(ps.target),
		beam(ps.beam),
		p(math::Vec4(ps.M, 0., 0., 0.)),
		k1(math::Vec4::from_length_and_t(ps.m, beam_energy, math::VEC3_Z)) { }

	/// Prepare the initial state with information from kinematics object
	/// \p kin, with the target at rest, and with the beam pointed along
	/// \f$\pmb{e}_z\f$.
	Initial(kin::Kinematics const& kin) :
		target(kin.target),
		beam(kin.beam),
		p(math::Vec4(kin.M, 0., 0., 0.)),
		k1(math::Vec4::from_length_and_t(kin.m, kin.S / (2. * kin.M), math::VEC3_Z)) { }
};

/**
 * %Final state of the system after a SIDIS process. Contains the scattered
 * lepton 4-momentum, the virtual photon 4-momentum, and the detected hadron
 * 4-momentum.
 */
struct Final {
	/// \copydoc part::Particles::beam
	part::Lepton beam;
	/// \copydoc part::Particles::hadron
	part::Hadron hadron;
	/// The 4-momentum of the virtual photon.
	math::Vec4 q;
	/// The final 4-momentum of the beam.
	math::Vec4 k2;
	/// The 4-momentum of the detected hadron.
	math::Vec4 ph;

	/// Calculate the final state from the Initial state \p init, the target
	/// polarization in its rest frame \p target_pol, and a Kinematics
	/// describing a SIDIS process.
	Final(Initial const& init, math::Vec3 target_pol, Kinematics const& kin);
};

/**
 * %Final state of the system after a radiative SIDIS process. Similar to Final,
 * but with an added radiated photon 4-momentum.
 * \sa Final
 */
struct FinalRad {
	/// \copydoc Final::beam
	part::Lepton beam;
	/// \copydoc Final::hadron
	part::Hadron hadron;
	/// \copydoc Final::q
	math::Vec4 q;
	/// \copydoc Final::k2
	math::Vec4 k2;
	/// \copydoc Final::ph
	math::Vec4 ph;
	/// The 4-momentum of the radiated photon.
	math::Vec4 k;

	/// Calculate the final state from the Initial state \p init, the target
	/// polarization in its rest frame \p target_pol, and a KinematicsRad
	/// describing a radiative SIDIS process.
	FinalRad(Initial const& init, math::Vec3 target_pol, KinematicsRad const& kin);
};
/// \}

}
}

#endif

