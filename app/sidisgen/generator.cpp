#include "generator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "params_format.hpp"

using namespace sidis;

namespace {

cut::Cut cut_from_params(Params& params) {
	Double const DEG = PI / 180.;
	cut::Cut result;
	result.x = params["cut.x"].any();
	result.y = params["cut.y"].any();
	result.z = params["cut.z"].any();
	result.ph_t_sq = params["cut.ph_t_sq"].any();
	result.phi_h = DEG * params["cut.phi_h"].any().as<math::Bound>();
	result.phi = DEG * params["cut.phi"].any().as<math::Bound>();
	result.Q_sq = params["cut.Q_sq"].any();
	result.t = params["cut.t"].any();
	result.W_sq = params["cut.W_sq"].any();
	result.r = params["cut.r"].any();
	result.mx_sq = params["cut.mx_sq"].any();
	result.qt_to_Q = params["cut.qt_to_Q"].any();
	result.lab_mom_q = params["cut.lab.mom_q"].any();
	result.lab_mom_k2 = params["cut.lab.mom_k2"].any();
	result.lab_mom_h = params["cut.lab.mom_h"].any();
	result.lab_theta_q = DEG * params["cut.lab.theta_q"].any().as<math::Bound>();
	result.lab_theta_k2 = DEG * params["cut.lab.theta_k2"].any().as<math::Bound>();
	result.lab_theta_h = DEG * params["cut.lab.theta_h"].any().as<math::Bound>();
	return result;
}

cut::CutRad cut_rad_from_params(Params& params) {
	Double const DEG = PI / 180.;
	cut::CutRad result;
	if (params["mc.rad.enable"].any()) {
		result.tau = params["cut.tau"].any();
		result.phi_k = DEG * params["cut.phi_k"].any().as<math::Bound>();
		result.R = params["cut.R"].any();
		// The `k_0_bar` cut is mandatory.
		result.k_0_bar = params["cut.k_0_bar"].any().as<math::Bound>()
			& math::Bound(params["phys.soft_threshold"].any(), INF);
		result.lab_mom_k = params["cut.lab.mom_k"].any();
		result.lab_theta_k = DEG * params["cut.lab.theta_k"].any().as<math::Bound>();
	}
	return result;
}

}

NradDensity::NradDensity(Params& params, sf::SfSet const& sf) :
		_cut(cut_from_params(params)),
		_sf(sf),
		_rc_method(params["phys.rc_method"].any()),
		_soft_threshold(params["phys.soft_threshold"].any()),
		_ps(
			params["setup.target"].any(),
			params["setup.beam"].any(),
			params["setup.hadron"].any(),
			params["phys.mass_threshold"].any()),
		_S(2. * mass(_ps.target) * params["setup.beam_energy"].any().as<Double>()),
		_beam_pol(params["setup.beam_pol"].any()),
		_target_pol(params["setup.target_pol"].any()) {
	if (_rc_method == RcMethod::APPROX || _rc_method == RcMethod::EXACT) {
		// Validate that if there is a `k_0_bar` cut, that its range completely
		// encompasses the non-radiative part from 0 to the soft threshold.
		if (params.is_set("cut.k_0_bar")) {
			math::Bound k_0_bar = params.get_soft<ValueBound>("cut.k_0_bar");
			if (!(k_0_bar.min() <= 0.) || !(k_0_bar.max() >= _soft_threshold)) {
				throw std::runtime_error(
					"Parameter 'cut.k_0_bar' does not encompass the full "
					"non-radiative range from 0 to 'phys.soft_threshold'.");
			}
		}
		// Check that no cuts incompatible with non-radiative events are
		// included. These include mostly cuts on the final state radiated
		// photon, which is unknown in the non-radiative case.
		Params params_cut_rad = params.filter("cut-no-nrad"_F);
		for (std::string const& name : params_cut_rad.names()) {
			if (params_cut_rad.is_set(name)) {
				throw std::runtime_error(
					"Parameter '" + name + "' is incompatible with "
					"non-radiative events.");
			}
		}
	}
}

Double NradDensity::transform(
		Point<6> const& unit_vec,
		kin::Kinematics* kin) const noexcept {
	Double jacobian;
	if (!cut::take(_cut, _ps, _S, unit_vec.data(), kin, &jacobian)) {
		jacobian = 0.;
	}
	return jacobian;
}

Double NradDensity::eval(Point<6> const& unit_vec, kin::Kinematics* kin) const noexcept {
	Double jacobian;
	if (!cut::take(_cut, _ps, _S, unit_vec.data(), kin, &jacobian)) {
		return 0.;
	}
	math::Vec3 eta = frame::hadron_from_target(*kin) * _target_pol;
	// TODO: Evaluate when it is a good approximation to say that
	// `nrad ~ nrad_ir`. This happens because for small `k_0_bar`, the
	// contribution of `rad_f` integrated up to `k_0_bar` becomes vanishingly
	// small, so it can be neglected. However, this must be balanced with
	// choosing `k_0_bar` to be non-zero to avoid the infrared divergence in the
	// radiative part of the cross-section.
	Double xs;
	switch (_rc_method) {
	case RcMethod::NONE:
		xs = xs::born(*kin, _sf, _beam_pol, eta);
		break;
	case RcMethod::APPROX:
		xs = xs::nrad_ir(*kin, _sf, _beam_pol, eta, _soft_threshold);
		break;
	case RcMethod::EXACT:
		xs = xs::nrad_integ(*kin, _sf, _beam_pol, eta, _soft_threshold).val;
		break;
	default:
		UNREACHABLE();
	}
	// Some kinematic regions will be out of range for the structure functions,
	// so return 0 in those cases.
	// TODO: Find a way of notifying when situations like this occur.
	if (!std::isfinite(xs) || !(xs >= 0.)) {
		return 0.;
	} else {
		return jacobian * xs;
	}
}

Double NradDensity::eval(Point<6> const& unit_vec) const noexcept {
	kin::Kinematics kin;
	return eval(unit_vec, &kin);
}

RadDensity::RadDensity(Params& params, sf::SfSet const& sf) :
		_cut(cut_from_params(params)),
		_cut_rad(cut_rad_from_params(params)),
		_sf(sf),
		_rc_method(params["phys.rc_method"].any()),
		_soft_threshold(params["phys.soft_threshold"].any()),
		_ps(
			params["setup.target"].any(),
			params["setup.beam"].any(),
			params["setup.hadron"].any(),
			params["phys.mass_threshold"].any()),
		_S(2. * mass(_ps.target) * params["setup.beam_energy"].any().as<Double>()),
		_beam_pol(params["setup.beam_pol"].any()),
		_target_pol(params["setup.target_pol"].any()) {
	if (_rc_method == RcMethod::NONE) {
		throw std::runtime_error(
			"Cannot enable radiative events while parameter 'phys.rc_method' "
			"has a value of 'none'.");
	}
	// Validate that the `k_0_bar` maximum is above the soft threshold.
	math::Bound k_0_bar = params["cut.k_0_bar"].any();
	if (!(k_0_bar.max() > _soft_threshold)) {
		throw std::runtime_error(
			"Parameter 'cut.k_0_bar' does not encompass any of the radiative "
			"range above 'phys.soft_threshold'.");
	}
}

Double RadDensity::transform(
		Point<9> const& unit_vec,
		kin::KinematicsRad* kin) const noexcept {
	Double jacobian;
	if (!cut::take(_cut, _cut_rad, _ps, _S, unit_vec.data(), kin, &jacobian)) {
		jacobian = 0.;
	}
	return jacobian;
}

Double RadDensity::eval(Point<9> const& unit_vec, kin::KinematicsRad* kin_rad) const noexcept {
	Double jacobian;
	if (!cut::take(_cut, _cut_rad, _ps, _S, unit_vec.data(), kin_rad, &jacobian)) {
		return 0.;
	}
	kin::Kinematics kin = kin_rad->project();
	math::Vec3 eta = frame::hadron_from_target(kin) * _target_pol;
	Double xs = xs::rad(*kin_rad, _sf, _beam_pol, eta);
	if (!std::isfinite(xs) || !(xs >= 0.)) {
		return 0.;
	} else {
		return jacobian * xs;
	}
}

Double RadDensity::eval(Point<9> const& unit_vec) const noexcept {
	kin::KinematicsRad kin_rad;
	return eval(unit_vec, &kin_rad);
}

DistParams::DistParams(EventType event_type, Params& params_full) {
	dist_type = DistType::BUBBLE;
	Double target_eff = params_full[p_name_init_target_eff(event_type)].any();
	Double scale_exp = params_full[p_name_init_scale_exp(event_type)].any();
	Int max_cells = params_full[p_name_init_max_cells(event_type)].any();
	if (max_cells <= 0) {
		max_cells = 1;
	}
	new (&params.bubble) BubbleParams();
	params.bubble.check_samples = 16384;
	params.bubble.target_rel_var = std::expm1(-2. * std::log(target_eff));
	params.bubble.scale_exp_est = scale_exp;
	params.bubble.min_cell_explore_samples = 512;
	params.bubble.hist_num_per_bin = 2;
	params.bubble.max_explore_cells = static_cast<std::size_t>(max_cells);
}

Generator::Generator(Density density) :
		_event_type(density.event_type),
		_dist_valid(false) {
	switch (_event_type) {
	case EventType::NRAD:
		new (&_density.nrad) NradDensity(density.density.nrad);
		new (&_dist.nrad) Dist<6>(DistType::UNIFORM);
		_dist_valid = true;
		break;
	case EventType::RAD:
		new (&_density.rad) RadDensity(density.density.rad);
		new (&_dist.rad) Dist<9>(DistType::UNIFORM);
		_dist_valid = true;
		break;
	default:
		UNREACHABLE();
	}
}

Generator::Generator(Generator&& other) noexcept :
		_event_type(other._event_type),
		_dist_valid(false) {
	switch (_event_type) {
	case EventType::NRAD:
		new (&_density.nrad) NradDensity(std::move(other._density.nrad));
		new (&_dist.nrad) Dist<6>(std::move(other._dist.nrad));
		_dist_valid = true;
		break;
	case EventType::RAD:
		new (&_density.rad) RadDensity(std::move(other._density.rad));
		new (&_dist.rad) Dist<9>(std::move(other._dist.rad));
		_dist_valid = true;
		break;
	default:
		UNREACHABLE();
	}
}

Generator::~Generator() {
	if (_dist_valid) {
		switch (_event_type) {
		case EventType::NRAD:
			_dist_valid = false;
			_dist.nrad.~Dist<6>();
			break;
		case EventType::RAD:
			_dist_valid = false;
			_dist.rad.~Dist<9>();
			break;
		default:
			UNREACHABLE();
		}
	}
}

Double Generator::prime() const {
	switch (_event_type) {
	case EventType::NRAD:
		return _dist.nrad.prime();
	case EventType::RAD:
		return _dist.rad.prime();
	default:
		UNREACHABLE();
	}
}

Event Generator::draw(RndEngine& rnd) const {
	Event event;
	event.event_type = _event_type;
	switch (_event_type) {
	case EventType::NRAD:
		{
			UnitEvent<6> unit_event_nrad = _dist.nrad.draw(rnd);
			event.weight = unit_event_nrad.weight * _density.nrad.eval(unit_event_nrad.vec, &event.kin.nrad);
		}
		break;
	case EventType::RAD:
		{
			UnitEvent<9> unit_event_rad = _dist.rad.draw(rnd);
			event.weight = unit_event_rad.weight * _density.rad.eval(unit_event_rad.vec, &event.kin.rad);
		}
		break;
	default:
		UNREACHABLE();
	}
	return event;
}

void Generator::build_dist(DistParams dist_params) {
	switch (_event_type) {
	case EventType::NRAD:
		_dist_valid = false;
		_dist.nrad.~Dist<6>();
		new (&_dist.nrad) Dist<6>(build_dist_approx<6>(dist_params, _density.nrad));
		_dist_valid = true;
		break;
	case EventType::RAD:
		_dist_valid = false;
		_dist.rad.~Dist<9>();
		new (&_dist.rad) Dist<9>(build_dist_approx<9>(dist_params, _density.rad));
		_dist_valid = true;
		break;
	default:
		UNREACHABLE();
	}
}

std::ostream& Generator::write_dist(std::ostream& os, Generator const& gen) {
	switch (gen._event_type) {
	case EventType::NRAD:
		return Dist<6>::write(os, gen._dist.nrad);
	case EventType::RAD:
		return Dist<9>::write(os, gen._dist.rad);
	default:
		UNREACHABLE();
	}
}

std::istream& Generator::read_dist(std::istream& is, Generator& gen) {
	switch (gen._event_type) {
	case EventType::NRAD:
		return Dist<6>::read(is, gen._dist.nrad);
	case EventType::RAD:
		return Dist<9>::read(is, gen._dist.rad);
	default:
		UNREACHABLE();
	}
}

