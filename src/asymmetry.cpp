#include "sidis/asymmetry.hpp"

#include <cmath>

#include "sidis/cross_section.hpp"
#include "sidis/cut.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/particle.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/vector.hpp"
#include "sidis/extra/integrate.hpp"

using namespace sidis;
using namespace sidis::asym;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::part;

namespace {

// Integrates over `phi_h` only, with fixed values for `eta_1` and `eta_2`. This
// is used by `ut_integ` to compute the full integral over `phi` and `phi_h`.
EstErr ut_integ_h(
		sf::SfSet const& sf_set,
		Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
		Real eta_1, int phi_h_coeff_1, int offset_1,
		Real eta_2, int phi_h_coeff_2, int offset_2,
		bool include_rc,
		IntegParams params) {
	Real Q_sq = S * x * y;
	sf::SfUP sf = sf_set.sf_up(ps.hadron, x, z, Q_sq, ph_t_sq);
	if (!include_rc) {
		EstErr born_integ = integrate(
			[&](Real phi_h) {
				PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(ps, S, ph_space);
				xs::Born b(kin);
				lep::LepBornUX lep(kin);
				had::HadUP had(kin, sf);
				Vec3 eta(
					eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
					eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
					0.);
				Real result = eta.x * xs::born_ut1_base(b, lep, had)
					+ eta.y * xs::born_ut2_base(b, lep, had);
				return result;
			},
			-PI, PI,
			params);
		return born_integ;
	} else {
		EstErr nrad_ir_integ = integrate(
			[&](Real phi_h) {
				PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(ps, S, ph_space);
				xs::Nrad b(kin, INF);
				lep::LepNradUX lep(kin);
				had::HadUP had(kin, sf);
				Vec3 eta(
					eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
					eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
					0.);
				return eta.x * xs::nrad_ir_ut1_base(b, lep, had)
					+ eta.y * xs::nrad_ir_ut2_base(b, lep, had);
			},
			-PI, PI,
			params);
		EstErr rad_f_integ = integrate(
			[&](std::array<Real, 4> ph) {
				Bound phi_h_b(-PI, PI);
				Real phi_h = phi_h_b.lerp(ph[0]);
				PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(ps, S, ph_space);

				KinematicsRad kin_rad;
				Real jacobian;
				Real ph_rad[3] = { ph[1], ph[2], ph[3] };
				if (!cut::take(kin, ph_rad, &kin_rad, &jacobian)) {
					return 0.;
				}
				jacobian *= phi_h_b.size();

				xs::Rad b(kin_rad);
				sf::SfUP shift_sf = sf_set.sf_up(
					ps.hadron,
					kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
				lep::LepRadUX lep(kin_rad);
				had::HadRadFUP had(kin_rad, sf, shift_sf);

				Vec3 eta(
					eta_1 * std::sin(phi_h_coeff_1 * phi_h + 0.5 * offset_1 * PI),
					eta_2 * std::sin(phi_h_coeff_2 * phi_h + 0.5 * offset_2 * PI),
					0.);
				Real xs = dot(eta, xs::rad_f_up_base(b, lep, had));
				if (std::isnan(xs)) {
					xs = 0.;
				}
				return jacobian * xs;
			},
			std::array<Real, 4>{ 0., 0., 0., 0. },
			std::array<Real, 4>{ 1., 1., 1., 1. },
			params);
		Real val = nrad_ir_integ.val + rad_f_integ.val;
		Real err = std::hypot(nrad_ir_integ.err, rad_f_integ.err);
		return { val, err };
	}
}

}

EstErr asym::uu_integ(
		sf::SfSet const& sf_set,
		Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
		bool include_rc,
		IntegParams params) {
	// The `phi` integration only contributes a factor of `2 pi`, so we don't
	// need to evaluate it. This leaves the `phi_h` integration.
	Real Q_sq = S * x * y;
	// Calculate the structure functions ahead of time, because they are
	// constant over the integrals we'll be doing.
	sf::SfUU sf = sf_set.sf_uu(ps.hadron, x, z, Q_sq, ph_t_sq);
	if (!include_rc) {
		// Without radiative corrections, there is only the Born cross-section.
		EstErr born_integ = integrate(
			[&](Real phi_h) {
				PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(ps, S, ph_space);
				xs::Born born(kin);
				lep::LepBornUU lep(kin);
				had::HadUU had(kin, sf);
				return xs::born_uu_base(born, lep, had);
			},
			-PI, PI,
			params);
		Real val = 2. * PI * born_integ.val;
		Real err = 2. * PI * born_integ.err;
		return { val, err };
	} else {
		// Evaluate the "non-radiative" and "radiative" parts separately so we
		// can take combine the `phi_h` integration with the photon degrees of
		// freedom.
		EstErr nrad_ir_integ = integrate(
			[&](Real phi_h) {
				PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(ps, S, ph_space);
				xs::Nrad b(kin, INF);
				lep::LepNradUU lep(kin);
				had::HadUU had(kin, sf);
				return xs::nrad_ir_uu_base(b, lep, had);
			},
			-PI, PI,
			params);
		EstErr rad_f_integ = integrate(
			[&](std::array<Real, 4> ph) {
				Bound phi_h_b(-PI, PI);
				Real phi_h = phi_h_b.lerp(ph[0]);
				PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, 0. };
				Kinematics kin(ps, S, ph_space);

				KinematicsRad kin_rad;
				Real jacobian;
				Real ph_rad[3] = { ph[1], ph[2], ph[3] };
				if (!cut::take(kin, ph_rad, &kin_rad, &jacobian)) {
					return 0.;
				}
				jacobian *= phi_h_b.size();

				xs::Rad b(kin_rad);
				sf::SfUU shift_sf = sf_set.sf_uu(
					ps.hadron,
					kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
				lep::LepRadUU lep(kin_rad);
				had::HadRadFUU had(kin_rad, sf, shift_sf);
				Real xs = xs::rad_f_uu_base(b, lep, had);
				if (std::isnan(xs)) {
					xs = 0.;
				}
				return jacobian * xs;
			},
			std::array<Real, 4>{ 0., 0., 0., 0. },
			std::array<Real, 4>{ 1., 1., 1., 1. },
			params);
		// It's unclear whether these two errors really are independent, but
		// it's likely a good enough approximation.
		Real val = 2. * PI * (nrad_ir_integ.val + rad_f_integ.val);
		Real err = 2. * PI * std::hypot(nrad_ir_integ.err, rad_f_integ.err);
		return { val, err };
	}
}

EstErr asym::ut_integ(
		sf::SfSet const& sf_set,
		Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
		int phi_s_coeff, int phi_h_coeff, int offset,
		bool include_rc,
		IntegParams params) {
	// Evaluate the `phi - phi_h` integral on the `eta` components. To do this
	// integral, we can approximate `phi_s` as `phi`. Then, analytically we get:
	if (phi_s_coeff == 0 || phi_s_coeff > 1 || phi_s_coeff < -1) {
		return { 0., 0. };
	}
	Real eta_1 = 0.5;
	Real eta_2 = 0.5 * phi_s_coeff;
	EstErr integ_h = ut_integ_h(
		sf_set,
		ps, S, x, y, z, ph_t_sq,
		eta_1, phi_s_coeff + phi_h_coeff, offset,
		eta_2, phi_s_coeff + phi_h_coeff, offset + 1,
		include_rc,
		params);
	Real val = 2. * PI * integ_h.val;
	Real err = 2. * PI * integ_h.err;
	return { val, err };
}

EstErr asym::ut_asym(
		sf::SfSet const& sf_set,
		part::Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
		int phi_s_coeff, int phi_h_coeff, int offset,
		bool include_rc,
		IntegParams params) {
	// Factor used to adjust uncertainty. This should be an over-estimate of the
	// true asymmetry value. It's guaranteed to be smaller than 1, but in
	// practice it will be even smaller than that, so leave it at 0.5.
	Real adjust = 0.5;
	IntegParams uu_params {
		params.method,
		params.num_evals,
		std::max(params.err_rel, 2. * adjust * params.err_abs),
		0.,
	};
	EstErr uu = uu_integ(
		sf_set,
		ps, S, x, y, z, ph_t_sq,
		include_rc,
		uu_params);
	IntegParams ut_params {
		params.method,
		params.num_evals,
		params.err_rel,
		2. * uu.val * params.err_abs,
	};
	EstErr ut = ut_integ(
		sf_set,
		ps, S, x, y, z, ph_t_sq,
		phi_s_coeff, phi_h_coeff, offset,
		include_rc,
		ut_params);
	Real val = 2. * ut.val / uu.val;
	Real err = std::hypot(ut.err / ut.val, uu.err / uu.val) * std::abs(val);
	return { val, err };
}

EstErr asym::sivers(
		sf::SfSet const& sf_set,
		Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
		bool include_rc,
		IntegParams params) {
	return ut_asym(sf_set, ps, S, x, y, z, ph_t_sq, -1, 1, 0, include_rc, params);
}
EstErr asym::collins(
		sf::SfSet const& sf_set,
		Particles ps, Real S, Real x, Real y, Real z, Real ph_t_sq,
		bool include_rc,
		IntegParams params) {
	return ut_asym(sf_set, ps, S, x, y, z, ph_t_sq, 1, -1, 0, include_rc, params);
}

