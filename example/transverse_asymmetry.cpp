#include <array>
#include <iostream>
#include <vector>

#include <cubature.hpp>

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>

#include <sidis/sidis.hpp>
#include <sidis/extra/bounds.hpp>
#include <sidis/extra/vector.hpp>
#include <sidis/sf_set/ww.hpp>

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::kin;
using namespace sidis::math;

namespace {

// Struct with information about how to apply radiative corrections.
struct RC {
	bool apply_rc;
	Real k0_cut;
	std::size_t num_evals;
	Real prec;

	static RC none() {
		return { false, 0., 100000, 1e-6 };
	}
	static RC cut(Real k0_cut) {
		return { true, k0_cut, 100000, 1e-6 };
	}
	static RC all() {
		return { true, INF, 100000, 1e-6 };
	}
};

sf::model::WW const model;
Hadron const hadron = Hadron::PI_P;
Real const M_th = MASS_P + MASS_PI_0;

// Integrates the unpolarized cross-section over the azimuthal angles `phi` and
// `phi_h`.
Real xs_uu_integ(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		RC rc_info=RC::none()) {
	// The `phi` integration only contributes a factor of `2 pi`, so we don't
	// need to evaluate it. This leaves the `phi_h` integration.
	Real S = 2. * dot(initial_state.p, initial_state.k1);
	Real Q_sq = S * x * y;
	// Calculate the structure functions ahead of time, because they are
	// constant over the integrals we'll be doing.
	sf::SfUU sf = model.sf_uu(hadron, x, z, Q_sq, ph_t_sq);
	if (!rc_info.apply_rc) {
		// Without radiative corrections, there is only the Born cross-section.
		cubature::EstErr<Real> xs_born_integ;
		xs_born_integ = cubature::cubature(
			[&](Real phi_h) {
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);
				xs::Born born(kin);
				lep::LepBornUU lep(kin);
				had::HadUU had(kin, sf);
				return xs::born_uu_base(born, lep, had);
			},
			-PI, PI,
			rc_info.num_evals, 0., rc_info.prec);
		return 2. * PI * xs_born_integ.val;
	} else {
		cubature::EstErr<Real> xs_nrad_ir_integ;
		cubature::EstErr<Real> xs_rad_f_integ;
		// Evaluate the "non-radiative" and "radiative" parts separately so we
		// can take combine the `phi_h` integration with the photon degrees of
		// freedom.
		xs_nrad_ir_integ = cubature::cubature(
			[&](Real phi_h) {
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);
				xs::NRadIR b(kin, rc_info.k0_cut);
				lep::LepBornUU lep_born(kin);
				lep::LepAmmUU lep_amm(kin);
				had::HadUU had(kin, sf);
				return xs::nrad_ir_uu_base(b, lep_born, lep_amm, had);
			},
			-PI, PI,
			rc_info.num_evals, 0., rc_info.prec);
		xs_rad_f_integ = cubature::cubature(
			[&](cubature::Point<4, Real> ph) {
				Bounds phi_h_b(0., 2. * PI);
				Real phi_h = phi_h_b.lerp(ph[0]);
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);

				// Rescale the integral so it is evaluated on the unit 4-cube.
				Bounds tau_b = tau_bounds(kin);
				Real tau = tau_b.lerp(ph[1]);
				Bounds phi_k_b(0., 2. * PI);
				Real phi_k = phi_k_b.lerp(ph[2]);
				Bounds R_b = R_bounds_soft(kin, tau, phi_k, rc_info.k0_cut);
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
				had::HadRadFUU had(kin_rad, sf, shift_sf);
				Real xs = xs::rad_f_uu_base(b, lep, had);
				if (std::isnan(xs)) {
					xs = 0.;
				}
				return jacobian * xs;
			},
			cubature::Point<4, Real>{ 0., 0., 0., 0. },
			cubature::Point<4, Real>{ 1., 1., 1., 1. },
			rc_info.num_evals, 0., rc_info.prec);
		return 2. * PI * (xs_nrad_ir_integ.val + xs_rad_f_integ.val);
	}
}

// Integrates the transverse part of the cross-section `sigma_UT - sigma_UU`
// over `phi_h` with fixed values for `eta_1` and `eta_2`. This is used by
// `xs_ut_integ` to compute the integral over both `phi` and `phi_h`.
Real xs_ut_integ_h(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		Real eta_1=0., int phi_h_coeff_1=0, int offset_1=0,
		Real eta_2=0., int phi_h_coeff_2=0, int offset_2=0,
		RC rc_info=RC::none()) {
	Real S = 2. * dot(initial_state.p, initial_state.k1);
	Real Q_sq = S * x * y;
	sf::SfUP sf = model.sf_up(hadron, x, z, Q_sq, ph_t_sq);
	if (!rc_info.apply_rc) {
		cubature::EstErr<Real> xs_born_integ;
		xs_born_integ = cubature::cubature(
			[&](Real phi_h) {
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
			rc_info.num_evals, 0., rc_info.prec);
		return xs_born_integ.val;
	} else {
		cubature::EstErr<Real> xs_nrad_ir_integ;
		cubature::EstErr<Real> xs_rad_f_integ;
		xs_nrad_ir_integ = cubature::cubature(
			[&](Real phi_h) {
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);
				xs::NRadIR b(kin, rc_info.k0_cut);
				lep::LepBornUX lep_born(kin);
				lep::LepAmmUX lep_amm(kin);
				had::HadUP had(kin, sf);
				Vec3 eta(
					eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
					eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
					0.);
				return eta.x * xs::nrad_ir_ut1_base(b, lep_born, lep_amm, had)
					+ eta.y * xs::nrad_ir_ut2_base(b, lep_born, lep_amm, had);
			},
			-PI, PI,
			rc_info.num_evals, 0., rc_info.prec);
		xs_rad_f_integ = cubature::cubature(
			[&](cubature::Point<4, Real> ph) {
				Bounds phi_h_b(-PI, PI);
				Real phi_h = phi_h_b.lerp(ph[0]);
				PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(initial_state, phase_space, hadron, M_th);

				Bounds tau_b = tau_bounds(kin);
				Real tau = tau_b.lerp(ph[1]);
				Bounds phi_k_b(-PI, PI);
				Real phi_k = phi_k_b.lerp(ph[2]);
				Bounds R_b = R_bounds_soft(kin, tau, phi_k, rc_info.k0_cut);
				Real R = R_b.lerp(ph[3]);
				Real jacobian = phi_h_b.size()
					* tau_b.size()
					* phi_k_b.size()
					* R_b.size();

				KinematicsRad kin_rad(kin, tau, phi_k, R);
				xs::Rad b(kin_rad);
				sf::SfUP shift_sf = model.sf_up(
					hadron,
					kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
				lep::LepRadUX lep(kin_rad);
				had::HadRadFUP had(kin_rad, sf, shift_sf);

				Vec3 eta(
					eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
					eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
					0.);
				Real xs = dot(
					eta,
					xs::rad_f_up_base(b, lep, had));
				if (std::isnan(xs)) {
					xs = 0.;
				}
				return jacobian * xs;
			},
			cubature::Point<4, Real>{ 0., 0., 0., 0. },
			cubature::Point<4, Real>{ 1., 1., 1., 1. },
			rc_info.num_evals, 0., rc_info.prec);
		return xs_nrad_ir_integ.val + xs_rad_f_integ.val;
	}
}

// Integrates the transverse part of the cross-section over the azimuthal angles
// `phi` and `phi_h`.
Real xs_ut_integ(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		int phi_s_coeff=0, int phi_h_coeff=0, int offset=0,
		RC rc_info=RC::none()) {
	// Evaluate the `phi - phi_h` integral on the `eta` components. To do this
	// integral, we can approximate `phi_s` as `phi`. Then, analytically we get:
	if (phi_s_coeff == 0 || phi_s_coeff > 1 || phi_s_coeff < -1) {
		return 0.;
	}
	Real eta_1 = 0.5;
	Real eta_2 = 0.5 * phi_s_coeff;
	return 2. * PI * xs_ut_integ_h(
		initial_state, x, y, z, ph_t_sq,
		eta_1, phi_s_coeff + phi_h_coeff, offset,
		eta_2, phi_s_coeff + phi_h_coeff, offset + 1,
		rc_info);
}

// Calculates a transverse asymmetry associated with the angle:
// ```
// sin(phi_s_coeff * phi + phi_h_coeff * phi_h + 0.5 * offset * PI)
// ```
Real asymmetry(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		int phi_s_coeff, int phi_h_coeff, int offset,
		RC rc_info=RC::none()) {
	Real xs_ut = xs_ut_integ(
		initial_state, x, y, z, ph_t_sq,
		phi_s_coeff, phi_h_coeff, offset,
		rc_info);
	Real xs_uu = xs_uu_integ(initial_state, x, y, z, ph_t_sq, rc_info);
	return 2. * xs_ut / xs_uu;
}

Real collins(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		RC rc_info=RC::none()) {
	return asymmetry(initial_state, x, y, z, ph_t_sq, 1, 1, 0, rc_info);
}

Real sivers(
		Initial initial_state, Real x, Real y, Real z, Real ph_t_sq,
		RC rc_info=RC::none()) {
	return asymmetry(initial_state, x, y, z, ph_t_sq, -1, 1, 0, rc_info);
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
	std::array<Real, 3> ph_t_vals { 0.7, 0.4, 0.1 };
	// The bounds on `z` must be slightly more restrictive than `[0, 1]` for
	// kinematic reasons.
	Real z_min = 0.1;
	Real z_max = 0.8;
	Real k0_cut = INF;

	TF1 f_sivers("Sivers Born", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return sivers(initial_state, x, y, z, ph_t_sq, RC::none());
	}, z_min, z_max, 1);
	TF1 f_sivers_rc("Sivers RC", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return sivers(initial_state, x, y, z, ph_t_sq, RC::cut(k0_cut));
	}, z_min, z_max, 1);
	TF1 f_collins("Collins Born", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return collins(initial_state, x, y, z, ph_t_sq, RC::none());
	}, z_min, z_max, 1);
	TF1 f_collins_rc("Collins RC", [&](Double_t* xs, Double_t* ps) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(ps[0] * ps[0]);
		return collins(initial_state, x, y, z, ph_t_sq, RC::cut(k0_cut));
	}, z_min, z_max, 1);

	Color_t color = 1;
	Style_t style = 0;
	Style_t style_rc = 1;
	f_sivers.SetLineColor(color);
	f_sivers.SetLineStyle(style);
	f_sivers_rc.SetLineColor(color);
	f_sivers_rc.SetLineStyle(style_rc);
	f_collins.SetLineColor(color);
	f_collins.SetLineStyle(style);
	f_collins_rc.SetLineColor(color);
	f_collins_rc.SetLineStyle(style_rc);

	TCanvas canvas_sivers("Sivers asymmetry");
	canvas_sivers.DrawFrame(0., 0., 1., 0.15, "Sivers asymmetry");
	for (Real ph_t : ph_t_vals) {
		f_sivers.SetParameter(0, ph_t);
		f_sivers_rc.SetParameter(0, ph_t);
		f_sivers.DrawCopy("LSAME");
		f_sivers_rc.DrawCopy("LSAME");
	}

	TCanvas canvas_collins("Collins asymmetry");
	canvas_collins.DrawFrame(0., 0., 1., 0.04, "Collins asymmetry");
	for (Real ph_t : ph_t_vals) {
		f_collins.SetParameter(0, ph_t);
		f_collins_rc.SetParameter(0, ph_t);
		f_collins.DrawCopy("LSAME");
		f_collins_rc.DrawCopy("LSAME");
	}

	TFile output("out.root", "RECREATE");
	f_sivers.Write();
	f_sivers_rc.Write();
	f_collins.Write();
	f_collins_rc.Write();
	canvas_sivers.Write();
	canvas_collins.Write();
	output.Close();

	app.Run();

	return 0;
}

