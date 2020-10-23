#include <array>
#include <iostream>

#include <cubature.hpp>

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>

#include <sidis/sidis.hpp>
#include <sidis/sf_model/ww.hpp>
#include <sidis/extra/bounds.hpp>
#include <sidis/extra/transform.hpp>
#include <sidis/extra/vector.hpp>

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::kin;
using namespace sidis::math;

namespace {

// Enum for how many radiative corrections to apply.
enum class RC {
	None,
	NonRadOnly,
	All,
};

sf::model::WW const model;
Hadron const hadron = Hadron::PI_P;
Real const M_th = MASS_P + MASS_PI_0;

// Integrates the unpolarized cross-section over the azimuthal angles `phi` and
// `phi_h`.
Real xs_uu_integ(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		RC include_rc=RC::None) {
	// The `phi` integration only contributes a factor of `2 pi`, so we don't
	// need to evaluate it. This leaves the `phi_h` integration.
	Real S = 2. * dot(initial_state.p, initial_state.k1);
	Real Q_sq = S * x * y;
	// Calculate the structure functions ahead of time, because they are
	// constant over the integrals we'll be doing.
	sf::SfUU sf = model.sf_uu(hadron, x, z, Q_sq, ph_t_sq);
	cubature::EstErr<Real> nrad_xs_integ { 0., 0. };
	cubature::EstErr<Real> rad_xs_integ { 0., 0. };
	if (include_rc == RC::None) {
		// Without radiative corrections, there is only the Born cross-section.
		nrad_xs_integ = cubature::cubature([&](Real phi_h) {
			PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
			Kinematics kin(initial_state, phase_space, hadron, M_th);
			xs::Born born(kin);
			lep::LepBornUU lep(kin);
			had::HadUU had(kin, sf);
			return xs::born_uu_base(born, lep, had);
		},
		-PI, PI,
		10000, 0., 1e-4);
	} else {
		// Evaluate the "non-radiative" and "radiative" parts separately.
		nrad_xs_integ = cubature::cubature([&](Real phi_h) {
			PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
			Kinematics kin(initial_state, phase_space, hadron, M_th);
			xs::NRad b(kin);
			lep::LepBornUU lep_born(kin);
			lep::LepAmmUU lep_amm(kin);
			had::HadUU had(kin, sf);
			return xs::nrad_uu_base(b, lep_born, lep_amm, had);
		},
		-PI, PI,
		10000, 0., 1e-4);
		if (include_rc == RC::All) {
			rad_xs_integ = cubature::cubature([&](cubature::Point<4, Real> ph) {
				Bounds phi_h_b(0., 2. * PI);
				Real phi_h = phi_h_b.lerp(ph[0]);
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);

				// Rescale the integral so it is evaluated on the unit 4-cube.
				Bounds tau_b = tau_bounds(kin);
				Real tau = tau_b.lerp(ph[1]);
				Bounds phi_k_b(0., 2. * PI);
				Real phi_k = phi_k_b.lerp(ph[2]);
				Bounds R_b = R_bounds(kin, tau, phi_k);
				Real R = R_b.lerp(ph[3]);
				Real jacobian = phi_h_b.size()
					* tau_b.size()
					* phi_k_b.size()
					* R_b.size();

				KinematicsRad kin_rad(kin, tau, phi_k, R);
				xs::Rad b(kin_rad);
				sf::SfUU shift_sf = model.sf_uu(
					hadron,
					kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
				lep::LepRadUU lep(kin_rad);
				had::HadUU had(kin, sf);
				had::HadUU shift_had(kin_rad.project_shift(), shift_sf);

				// Compute hard cross-section.
				Real uu_h = xs::rad_hard_uu_base(b, lep, shift_had);
				// Compute soft cross-section.
				Real uu_s;
				if (std::abs(R) < std::abs(R_b.max) * xs::SMALL_R_REL) {
					uu_s = xs::rad_soft_0(0., Vec3::ZERO, kin_rad, model);
				} else {
					uu_s = xs::rad_soft_uu_base(b, lep, had, shift_had);
				}
				Real xs = uu_h + uu_s;
				if (std::isnan(xs)) {
					xs = 0.;
				}
				return jacobian * xs;
			},
			cubature::Point<4, Real>{ 0., 0., 0., 0. },
			cubature::Point<4, Real>{ 1., 1., 1., 1. },
			10000, 0., 1e-4);
		}
	}

	Real xs_integ = nrad_xs_integ.val + rad_xs_integ.val;
	Real err = nrad_xs_integ.err + rad_xs_integ.err;
	return 2. * PI * xs_integ;
}

// Integrates the transverse part of the cross-section `sigma_UT - sigma_UU`
// over `phi_h` with fixed values for `eta_1` and `eta_2`. This is used by
// `xs_ut_integ` to compute the integral over both `phi` and `phi_h`.
Real xs_ut_integ_h(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		Real eta_1=0., int phi_h_coeff_1=0, int offset_1=0,
		Real eta_2=0., int phi_h_coeff_2=0, int offset_2=0,
		RC include_rc=RC::None) {
	Real S = 2. * dot(initial_state.p, initial_state.k1);
	Real Q_sq = S * x * y;
	sf::SfUP sf = model.sf_up(hadron, x, z, Q_sq, ph_t_sq);
	cubature::EstErr<Real> nrad_xs_integ { 0., 0. };
	cubature::EstErr<Real> rad_xs_integ { 0., 0. };
	if (include_rc == RC::None) {
		nrad_xs_integ = cubature::cubature([&](Real phi_h) {
			PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
			Kinematics kin(initial_state, phase_space, hadron, M_th);
			xs::Born b(kin);
			lep::LepBornUX lep(kin);
			had::HadUP had(kin, sf);
			Vec3 eta(
				eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
				eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
				0.);
			return eta.x * xs::born_ut1_base(b, lep, had)
				+ eta.y * xs::born_ut2_base(b, lep, had);
		},
		-PI, PI,
		10000, 0., 1e-4);
	} else {
		nrad_xs_integ = cubature::cubature([&](Real phi_h) {
			PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
			Kinematics kin(initial_state, phase_space, hadron, M_th);
			xs::NRad b(kin);
			lep::LepBornUX lep_born(kin);
			lep::LepAmmUX lep_amm(kin);
			had::HadUP had(kin, sf);
			Vec3 eta(
				eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
				eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
				0.);
			return eta.x * xs::nrad_ut1_base(b, lep_born, lep_amm, had)
				+ eta.y * xs::nrad_ut2_base(b, lep_born, lep_amm, had);
		},
		-PI, PI,
		10000, 0., 1e-4);
		if (include_rc == RC::All) {
			rad_xs_integ = cubature::cubature([&](cubature::Point<4, Real> ph) {
				Bounds phi_h_b(-PI, PI);
				Real phi_h = phi_h_b.lerp(ph[0]);
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);

				Bounds tau_b = tau_bounds(kin);
				Real tau = tau_b.lerp(ph[1]);
				Bounds phi_k_b(-PI, PI);
				Real phi_k = phi_k_b.lerp(ph[2]);
				Bounds R_b = R_bounds(kin, tau, phi_k);
				Real R = R_b.lerp(ph[3]);
				Real jacobian = phi_h_b.size()
					* tau_b.size()
					* phi_k_b.size()
					* R_b.size();

				KinematicsRad kin_rad(kin, tau, phi_k, R);
				xs::Rad b(kin_rad);
				Transform3 shift_rot = frame::hadron_from_shift(kin_rad);
				sf::SfUP shift_sf = model.sf_up(
					hadron,
					kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
				lep::LepRadUX lep(kin_rad);
				had::HadUP had(kin, sf);
				had::HadUP shift_had(kin_rad.project_shift(), shift_sf);

				Vec3 eta(
					eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
					eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
					0.);
				Real up_h = dot(
					eta,
					xs::rad_hard_up_base(b, lep, shift_had, shift_rot));
				Real up_s;
				if (std::abs(R) < std::abs(R_b.max) * xs::SMALL_R_REL) {
					up_s = xs::rad_soft_0(0., eta, kin_rad, model);
				} else {
					up_s = dot(
						eta,
						xs::rad_soft_up_base(b, lep, had, shift_had, shift_rot));
				}
				Real xs = up_h + up_s;
				if (std::isnan(xs)) {
					xs = 0.;
				}
				return jacobian * xs;
			},
			cubature::Point<4, Real>{ 0., 0., 0., 0. },
			cubature::Point<4, Real>{ 1., 1., 1., 1. },
			10000, 0., 1e-4);
		}
	}

	Real xs_integ = nrad_xs_integ.val + rad_xs_integ.val;
	Real err = nrad_xs_integ.err + rad_xs_integ.err;
	return xs_integ;
}

// Integrates the transverse part of the cross-section over the azimuthal angles
// `phi` and `phi_h`.
Real xs_ut_integ(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		int phi_s_coeff=0, int phi_h_coeff=0, int offset=0,
		RC include_rc=RC::None) {
	// Evaluate the `phi - phi_h` integral on the `eta` components. To do this
	// integral, we can approximate `phi_s` as `phi`. Then, analytically we get:
	if (phi_s_coeff == 0 || phi_s_coeff > 1 || phi_s_coeff < -1) {
		return 0.;
	}
	Real eta_1 = 0.5 * phi_s_coeff;
	Real eta_2 = 0.5;
	return 2. * PI * xs_ut_integ_h(
		initial_state, x, y, z, ph_t_sq,
		eta_1, phi_s_coeff + phi_h_coeff, offset + 1,
		eta_2, phi_s_coeff + phi_h_coeff, offset,
		include_rc);
}

// Calculates a transverse asymmetry associated with the angle:
// ```
// sin(phi_s_coeff * phi + phi_h_coeff * phi_h + 0.5 * offset * PI)
// ```
Real asymmetry(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		int phi_s_coeff, int phi_h_coeff, int offset,
		RC include_rc=RC::None) {
	Real xs_ut = xs_ut_integ(
		initial_state, x, y, z, ph_t_sq,
		phi_s_coeff, phi_h_coeff, offset,
		include_rc);
	Real xs_uu = xs_uu_integ(initial_state, x, y, z, ph_t_sq, include_rc);
	return 2. * xs_ut / xs_uu;
}

Real collins(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		RC include_rc=RC::None) {
	return asymmetry(initial_state, x, y, z, ph_t_sq, 1, 1, 0, include_rc);
}

Real sivers(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		RC include_rc=RC::None) {
	return asymmetry(initial_state, x, y, z, ph_t_sq, -1, 1, 0, include_rc);
}

}

// Computes and plots the Sivers and Collins asymmetries, with and without
// radiative corrections. The computation demonstrates some techniques for
// getting better performance out of the `sidis` library for cases where the
// target and beam polarizations are known.
int main(int argc, char** argv) {
	TApplication app("Transverse asymmetry plots", &argc, argv);

	Initial initial_state(Nucleus::P, Lepton::E, 10.6);
	Real x = 0.2;
	Real y = 0.754;

	// The values of `ph_t` to evaluate the transverse asymmetries.
	std::array<Real, 3> ph_t_vals { 0.1, 0.4, 0.7 };
	// The bounds on `z` must be slightly more restrictive than `[0, 1]` for
	// kinematic reasons.
	Real z_min = 0.1;
	Real z_max = 0.8;

	TF1 f_sivers("Sivers Born", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return sivers(initial_state, x, y, z, ph_t_sq, RC::None);
	}, z_min, z_max, 1);
	TF1 f_sivers_rc("Sivers RC", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return sivers(initial_state, x, y, z, ph_t_sq, RC::NonRadOnly);
	}, z_min, z_max, 1);
	TF1 f_collins("Collins Born", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return collins(initial_state, x, y, z, ph_t_sq, RC::None);
	}, z_min, z_max, 1);
	TF1 f_collins_rc("Collins RC", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return collins(initial_state, x, y, z, ph_t_sq, RC::NonRadOnly);
	}, z_min, z_max, 1);

	TCanvas canvas_sivers("Sivers asymmetry");
	for (Real ph_t : ph_t_vals) {
		f_sivers.SetParameter(0, ph_t);
		f_sivers_rc.SetParameter(0, ph_t);
		f_sivers.Draw();
		f_sivers_rc.Draw();
	}

	TCanvas canvas_collins("Collins asymmetry");
	for (Real ph_t : ph_t_vals) {
		f_collins.SetParameter(0, ph_t);
		f_collins_rc.SetParameter(0, ph_t);
		f_collins.Draw();
		f_collins_rc.Draw();
	}

	app.Run();

	return 0;
}

