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

/*
 * The following functions provide transformations between various useful frames
 * for the SIDIS process:
 * * `lab`: The lab frame is the frame in which the initial 4-momenta of all
 *   particles are defined.
 * * `target`: The target frame has the target at rest, the lepton beam along
 *   the `z` axis, and any transverse polarization of the target along the `y`
 *   axis.
 * * `lepton`: The lepton frame (see [1.A5]) has its `z` axis along the virtual
 *   photon momentum and `x` axis along the transverse part of the outgoing
 *   lepton momentum.
 * * `hadron`: The hadron frame (see [1.A2]) has its `z` axis along the virtual
 *   photon momentum and `x` axis along the transverse part of the detected
 *   hadron momentum.
 * * `virt_photon`: The virtual photon frame (see [1.A1]) shares its `x` and `y`
 *   axes with the virtual hadron frame, but has an additional boost along `z`
 *   to bring the virtual photon to rest.
 * * `real_photon`: The real photon frame has its `z` axis along the virtual
 *   photon momentum and `x` axis along the transverse part of the real photon
 *   momentum.
 *
 * For radiative kinematics only, the shifted frame is additionally given as:
 * * `shift`: The shifted form of the hadron frame after applying the
 *   transformation `q -> q - k`.
 *
 * To use these functions, generally start by transforming from the lab frame to
 * the target frame. From the target frame, you can then transform to the
 * various other frames using only the kinematics information in `kin::Kin`.
 */

math::Transform4 target_from_lab(kin::Initial const& init, math::Vec3 pol);
math::Transform4 lab_from_target(kin::Initial const& init, math::Vec3 pol);

math::Transform3 target_from_lepton(kin::Kinematics const& kin);
math::Transform3 lepton_from_target(kin::Kinematics const& kin);

math::Transform3 target_from_hadron(kin::Kinematics const& kin);
math::Transform3 hadron_from_target(kin::Kinematics const& kin);

math::Transform4 target_from_virt_photon(kin::Kinematics const& kin);
math::Transform4 virt_photon_from_target(kin::Kinematics const& kin);

math::Transform3 target_from_real_photon(kin::KinematicsRad const& kin);
math::Transform3 real_photon_from_target(kin::KinematicsRad const& kin);

math::Transform3 target_from_shift(kin::KinematicsRad const& kin);
math::Transform3 shift_from_target(kin::KinematicsRad const& kin);
math::Transform3 hadron_from_shift(kin::KinematicsRad const& kin);
math::Transform3 shift_from_hadron(kin::KinematicsRad const& kin);

}
}

#endif

