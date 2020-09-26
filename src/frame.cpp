#include "sidis/frame.hpp"

#include <iostream>
#include <cmath>

#include "sidis/kinematics.hpp"
#include "sidis/extra/math.hpp"
#include "sidis/extra/transform.hpp"

using namespace sidis;
using namespace sidis::frame;
using namespace sidis::kin;
using namespace sidis::math;

Transform4 frame::lab_from_target(Initial init, Vec3 pol) {
	Transform4 boost = Transform4::boost_to(init.p);
	Vec3 y_axis = pol.perp((boost * init.k1).r()).unit();
	// If the target polarization is parallel to the beam, then choose an
	// arbitrary direction to be the y-axis. One of these is guaranteed to be
	// non-zero.
	if (y_axis.norm_sq() == 0.) {
		y_axis = cross(Vec3::X, (boost * init.k1).r()).unit();
	}
	if (y_axis.norm_sq() == 0.) {
		y_axis = cross(Vec3::Y, (boost * init.k1).r()).unit();
	}
	if (y_axis.norm_sq() == 0.) {
		y_axis = cross(Vec3::Z, (boost * init.k1).r()).unit();
	}
	Transform4 rotate = Transform3::rotate_basis((boost * init.k1).r(), y_axis);
	return rotate * boost;
}
Transform4 frame::target_from_lab(Initial init, Vec3 pol) {
	return lab_from_target(init, pol).transpose();
}

Transform4 frame::target_from_lepton(Kinematics kin) {
	return lepton_from_target(kin).transpose();
}
Transform4 frame::lepton_from_target(Kinematics kin) {
	// Equation [1.A5].
	Real q_norm = kin.lambda_Y_sqrt/(2.*kin.M);
	Vec3 ez = 1./q_norm*Vec3(
		kin.q_t*kin.sin_phi,
		-kin.q_t*kin.cos_phi,
		kin.q_l);
	Vec3 ey = Vec3(-kin.cos_phi, -kin.sin_phi, 0.);
	Vec3 ex = 1./q_norm*Vec3(
		-kin.sin_phi*kin.q_l,
		kin.cos_phi*kin.q_l,
		kin.q_t);
	return Transform4(Transform3(ex, ey, ez));
}

Transform4 frame::target_from_hadron(Kinematics kin) {
	// Equation [1.A2].
	Transform3 rotate(
		kin.cos_phi_h, -kin.sin_phi_h, 0.,
		kin.sin_phi_h, kin.cos_phi_h, 0.,
		0., 0., 1.);
	return target_from_lepton(kin) * Transform4(rotate);
}
Transform4 frame::hadron_from_target(Kinematics kin) {
	return target_from_hadron(kin).transpose();
}

Transform4 frame::target_from_virt_photon(Kinematics kin) {
	// Equation [1.A3].
	Transform4 boost = 1./(2.*kin.M*kin.Q)*Transform4(
		kin.lambda_Y_sqrt, 0., 0., kin.S_x,
		0., 1., 0., 0.,
		0., 0., 1., 0.,
		kin.S_x, 0., 0., kin.lambda_Y_sqrt);
	return target_from_hadron(kin) * boost;
}
Transform4 frame::virt_photon_from_target(Kinematics kin) {
	return target_from_virt_photon(kin).transpose();
}

Transform4 frame::target_from_real_photon(KinematicsRad kin) {
	// Equation [1.A2].
	Transform3 rotate(
		kin.cos_phi_k, -kin.sin_phi_k, 0.,
		kin.sin_phi_k, kin.cos_phi_k, 0.,
		0., 0., 1.);
	return target_from_lepton(kin.project()) * Transform4(rotate);
}
Transform4 frame::real_photon_from_target(KinematicsRad kin) {
	return target_from_real_photon(kin).transpose();
}

Transform4 frame::target_from_shift(KinematicsRad kin) {
	// TODO: See if this can be calculated more accurately by doing it directly.
	return target_from_hadron(kin.project()) * hadron_from_shift(kin);
}
Transform4 frame::shift_from_target(KinematicsRad kin) {
	return target_from_shift(kin).transpose();
}

Transform4 frame::hadron_from_shift(KinematicsRad kin) {
	return shift_from_hadron(kin).transpose();
}
Transform4 frame::shift_from_hadron(KinematicsRad kin) {
	// TODO: Fill in equation number from derivations.
	Vec3 ex(
		1./(kin.shift_lambda_Y*kin.shift_ph_t)*(
			kin.lambda_Y*kin.ph_t
			- 1./(4.*sq(kin.M)*kin.lambda_Y*kin.ph_t)*(
				(sq(kin.R) - 2.*kin.lambda_RY)*(
					sq(kin.lambda_V) - kin.lambda_H*kin.lambda_Y)
				+ (kin.lambda_V - kin.lambda_RV)*(
					kin.lambda_RY*kin.lambda_V - kin.lambda_RV*kin.lambda_Y))),
		1./(kin.lambda_Y_sqrt*kin.shift_lambda_Y*kin.ph_t*kin.shift_ph_t)*(
			2.*(kin.lambda_V - kin.lambda_RV)*kin.vol_phi_hk),
		1./(2.*kin.M*kin.lambda_Y_sqrt*kin.shift_lambda_Y*kin.shift_ph_t)*(
			(sq(kin.R) - kin.lambda_RY)*kin.lambda_V
			+ (kin.lambda_Y - kin.lambda_RY)*kin.lambda_RV));
	Vec3 ey(
		-1./(kin.lambda_Y*kin.shift_lambda_Y_sqrt*kin.ph_t*kin.shift_ph_t)*(
			2.*kin.vol_phi_hk*kin.lambda_V),
		1./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t)*(
			kin.lambda_Y_sqrt*kin.ph_t
			+ 1./(4.*sq(kin.M)*kin.lambda_Y_sqrt*kin.ph_t)*(
				kin.lambda_V*kin.lambda_RV - kin.lambda_H*kin.lambda_RY)),
		1./(kin.lambda_Y_sqrt*kin.shift_lambda_Y_sqrt*kin.shift_ph_t)*(
			4.*kin.M*kin.vol_phi_hk));
	Vec3 ez(
		1./(2.*kin.M*kin.lambda_Y*kin.shift_lambda_Y_sqrt*kin.ph_t)*(
			kin.lambda_RY*kin.lambda_V - kin.lambda_RV*kin.lambda_Y),
		-1./(kin.lambda_Y_sqrt*kin.shift_lambda_Y_sqrt*kin.ph_t)*(
			4.*kin.M*kin.vol_phi_hk),
		1./kin.shift_lambda_Y_sqrt*(
			kin.lambda_Y_sqrt - kin.lambda_RY/kin.lambda_Y_sqrt));
	return Transform4(Transform3(ex, ey, ez));
}

