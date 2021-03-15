#ifndef SIDIS_FRAME_HPP
#define SIDIS_FRAME_HPP

namespace sidis {

namespace kin {
	struct Initial;
	struct Kinematics;
	struct KinematicsRad;
}
namespace math {
	struct Transform3;
	struct Transform4;
	struct Vec3;
}
namespace frame {

/**
 * \defgroup FrameGroup Frame transformations
 * The functions here provide transformations between various useful frames for
 * the SIDIS process:
 * * `lab`: The lab frame is the frame in which the initial 4-momenta of all
 *   particles are defined.
 * * `target`: The target frame has the target at rest, the incident lepton
 *   along the \f$\pmb{e}_z\f$ axis, and any transverse polarization of the
 *   target along the \f$\pmb{e}_y\f$ axis.
 * * `lepton`: The lepton frame has its \f$\pmb{e}_z\f$ axis along the virtual
 *   photon momentum and \f$\pmb{e}_x\f$ axis along the transverse part of the
 *   scattered lepton momentum.
 * * `hadron`: The hadron frame has its \f$\pmb{e}_z\f$ axis along the virtual
 *   photon momentum and \f$\pmb{e}_x\f$ axis along the transverse part of the
 *   detected hadron momentum.
 * * `virt_photon`: The virtual photon frame shares its \f$\pmb{e}_x\f$ and
 *   \f$\pmb{e}_y\f$ axes with the hadron frame, but has an additional boost
 *   along \f$\pmb{e}_z\f$ to bring the virtual photon to rest.
 *
 * For radiative kinematics only, the following frames are also used:
 * * `real_photon`: The radiated photon frame has its \f$\pmb{e}_z\f$ axis along
 *   the radiated photon momentum and \f$\pmb{e}_x\f$ axis along the transverse
 *   part of the radiated photon momentum.
 * * `shift`: The shifted form of the hadron frame after applying the
 *   transformation \f$q \rightarrow q - k\f$.
 *
 * To use these functions, generally start by transforming from the lab frame to
 * the target frame using frame::target_from_lab. From the target frame, you can
 * then transform to the various other frames using a kin::Kinematics or
 * kin::KinematicsRad.
 */
/// \{

/// Transformation from lab frame to target frame.
math::Transform4 target_from_lab(kin::Initial const& init, math::Vec3 pol);
/// Transformation from target frame to lab frame.
math::Transform4 lab_from_target(kin::Initial const& init, math::Vec3 pol);

/// Transformation from lepton frame to target frame.
math::Transform3 target_from_lepton(kin::Kinematics const& kin);
/// Transformation from target frame to lepton frame.
math::Transform3 lepton_from_target(kin::Kinematics const& kin);

/// Transformation from hadron frame to target frame.
math::Transform3 target_from_hadron(kin::Kinematics const& kin);
/// Transformation from target frame to hadron frame.
math::Transform3 hadron_from_target(kin::Kinematics const& kin);

/// Transformation from virtual photon frame to target frame.
math::Transform4 target_from_virt_photon(kin::Kinematics const& kin);
/// Transformation from target frame to virtual photon frame.
math::Transform4 virt_photon_from_target(kin::Kinematics const& kin);

/// Transformation from radiated photon frame to target frame.
math::Transform3 target_from_real_photon(kin::KinematicsRad const& kin);
/// Transformation from target frame to radiated photon frame.
math::Transform3 real_photon_from_target(kin::KinematicsRad const& kin);

/// Transformation from shift frame to target frame.
math::Transform3 target_from_shift(kin::KinematicsRad const& kin);
/// Transformation from target frame to shift frame.
math::Transform3 shift_from_target(kin::KinematicsRad const& kin);
/// Transformation from shift frame to hadron frame.
math::Transform3 hadron_from_shift(kin::KinematicsRad const& kin);
/// Transformation from hadron frame to shift frame.
math::Transform3 shift_from_hadron(kin::KinematicsRad const& kin);
/// \}

}
}

#endif

