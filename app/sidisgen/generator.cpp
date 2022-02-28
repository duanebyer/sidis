#include "generator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

using namespace sidis;

namespace {

cut::Cut get_cut_from_params(Params const& params) {
	Real const DEG = PI / 180.;
	cut::Cut result;
	result.x = params.cut_x.get_or(math::Bound::INVALID);
	result.y = params.cut_y.get_or(math::Bound::INVALID);
	result.z = params.cut_z.get_or(math::Bound::INVALID);
	result.ph_t_sq = params.cut_ph_t_sq.get_or(math::Bound::INVALID);
	result.phi_h = DEG * params.cut_phi_h.get_or(math::Bound::INVALID);
	result.phi = DEG * params.cut_phi.get_or(math::Bound::INVALID);
	result.Q_sq = params.cut_Q_sq.get_or(math::Bound::INVALID);
	result.t = params.cut_t.get_or(math::Bound::INVALID);
	result.W_sq = params.cut_W_sq.get_or(math::Bound::INVALID);
	result.r = params.cut_r.get_or(math::Bound::INVALID);
	result.mx_sq = params.cut_mx_sq.get_or(math::Bound::INVALID);
	result.q_0 = params.cut_q_0.get_or(math::Bound::INVALID);
	result.k2_0 = params.cut_k2_0.get_or(math::Bound::INVALID);
	result.ph_0 = params.cut_ph_0.get_or(math::Bound::INVALID);
	result.theta_q = DEG * params.cut_theta_q.get_or(math::Bound::INVALID);
	result.theta_k2 = DEG * params.cut_theta_k2.get_or(math::Bound::INVALID);
	result.theta_h = DEG * params.cut_theta_h.get_or(math::Bound::INVALID);
	return result;
}

cut::CutRad get_cut_rad_from_params(Params const& params) {
	Real const DEG = PI / 180.;
	cut::CutRad result;
	if (*params.rad_gen) {
		result.tau = params.cut_tau.get_or(math::Bound::INVALID);
		result.phi_k = DEG * params.cut_phi_k.get_or(math::Bound::INVALID);
		result.R = params.cut_R.get_or(math::Bound::INVALID);
		// The `k_0_bar` cut is mandatory.
		result.k_0_bar = *params.cut_k_0_bar
			& math::Bound(*params.k_0_bar, INF);
		result.k_0 = params.cut_k_0.get_or(math::Bound::INVALID);
		result.theta_k = DEG * params.cut_theta_k.get_or(math::Bound::INVALID);
	}
	return result;
}

std::random_device rnd_dev;

}

NradDensity::NradDensity(Params const& params, sf::SfSet const& sf) :
	_params(params),
	_cut(get_cut_from_params(params)),
	_ps(*params.target, *params.beam, *params.hadron, *params.Mth),
	_S(2. * mass(*params.target) * *params.beam_energy),
	_sf(sf) { }

Real NradDensity::transform(Point<6> const& vec, PointFull* ph_out) const noexcept {
	Real jacobian;
	kin::Kinematics kin;
	if (!cut::take(_cut, _ps, _S, vec.data(), &kin, &jacobian)) {
		jacobian = 0.;
	}
	(*ph_out)[0] = kin.x;
	(*ph_out)[1] = kin.y;
	(*ph_out)[2] = kin.z;
	(*ph_out)[3] = kin.ph_t_sq;
	(*ph_out)[4] = kin.phi_h;
	(*ph_out)[5] = kin.phi;
	return jacobian;
}

Real NradDensity::operator()(Point<6> const& vec) const noexcept {
	kin::Kinematics kin;
	Real jacobian;
	if (!cut::take(_cut, _ps, _S, vec.data(), &kin, &jacobian)) {
		return 0.;
	}
	math::Vec3 eta = frame::hadron_from_target(kin) * *_params.target_pol;
	// TODO: Evaluate when it is a good approximation to say that
	// `nrad ~ nrad_ir`. This happens because for small `k_0_bar`, the
	// contribution of `rad_f` integrated up to `k_0_bar` becomes vanishingly
	// small, so it can be neglected. However, this must be balanced with
	// choosing `k_0_bar` to be non-zero to avoid the infrared divergence in the
	// radiative part of the cross-section.
	Real xs;
	switch (*_params.rc_method) {
	case RcMethod::NONE:
		xs = xs::born(kin, _sf, *_params.beam_pol, eta);
		break;
	case RcMethod::APPROX:
		xs = xs::nrad_ir(kin, _sf, *_params.beam_pol, eta, *_params.k_0_bar);
		break;
	case RcMethod::EXACT:
		xs = xs::nrad_integ(kin, _sf, *_params.beam_pol, eta, *_params.k_0_bar).val;
		break;
	default:
		xs = 0.;
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

RadDensity::RadDensity(Params const& params, sf::SfSet const& sf) :
	_params(params),
	_cut(get_cut_from_params(params)),
	_cut_rad(get_cut_rad_from_params(params)),
	_ps(*params.target, *params.beam, *params.hadron, *params.Mth),
	_S(2. * mass(*params.target) * *params.beam_energy),
	_sf(sf) { }

Real RadDensity::transform(Point<9> const& vec, PointFull* ph_out) const noexcept {
	Real jacobian;
	kin::KinematicsRad kin;
	if (!cut::take(_cut, _cut_rad, _ps, _S, vec.data(), &kin, &jacobian)) {
		jacobian = 0.;
	}
	(*ph_out)[0] = kin.x;
	(*ph_out)[1] = kin.y;
	(*ph_out)[2] = kin.z;
	(*ph_out)[3] = kin.ph_t_sq;
	(*ph_out)[4] = kin.phi_h;
	(*ph_out)[5] = kin.phi;
	(*ph_out)[6] = kin.tau;
	(*ph_out)[7] = kin.phi_k;
	(*ph_out)[8] = kin.R;
	return jacobian;
}

Real RadDensity::operator()(Point<9> const& vec) const noexcept {
	kin::KinematicsRad kin_rad;
	Real jacobian;
	if (!cut::take(_cut, _cut_rad, _ps, _S, vec.data(), &kin_rad, &jacobian)) {
		return 0.;
	}
	kin::Kinematics kin = kin_rad.project();
	math::Vec3 eta = frame::hadron_from_target(kin) * *_params.target_pol;
	Real xs = xs::rad(kin_rad, _sf, *_params.beam_pol, eta);
	if (!std::isfinite(xs) || !(xs >= 0.)) {
		return 0.;
	} else {
		return jacobian * xs;
	}
}

Builder::Builder(
		EventType type,
		BuilderReporters const& builder_reporters,
		Params const& params,
		sf::SfSet const& sf) :
		_type(type) {
	std::minstd_rand seed_rnd(rnd_dev());
	std::minstd_rand seed_rnd_fixed;
	std::uniform_int_distribution<Seed> seed_dist;
	switch (_type) {
	case EventType::NRAD:
		seed_rnd_fixed.seed(*params.nrad_seed_init);
		new (&_builder.nrad) NradBuilder(
			NradDensity(params, sf),
			*params.nrad_seed_init == 0 ? seed_dist(seed_rnd) : seed_dist(seed_rnd_fixed));
		_builder.nrad.explore_progress_reporter = builder_reporters.explore_progress;
		_builder.nrad.tune_progress_reporter = builder_reporters.tune_progress;
		_builder.nrad.check_samples = 16384;
		_builder.nrad.target_rel_var =
			std::expm1(-2. * std::log(*params.nrad_target_eff));
		_builder.nrad.scale_exp_est = *params.nrad_scale_exp;
		_builder.nrad.min_cell_explore_samples = 32;
		_builder.nrad.hist_num_per_bin = 2;
		_builder.nrad.max_explore_cells = *params.nrad_max_cells;
		break;
	case EventType::RAD:
		seed_rnd_fixed.seed(*params.rad_seed_init);
		new (&_builder.rad) RadBuilder(
			RadDensity(params, sf),
			*params.rad_seed_init == 0 ? seed_dist(seed_rnd) : seed_dist(seed_rnd_fixed));
		_builder.rad.explore_progress_reporter = builder_reporters.explore_progress;
		_builder.rad.tune_progress_reporter = builder_reporters.tune_progress;
		_builder.rad.check_samples = 16384;
		_builder.rad.target_rel_var =
			std::expm1(-2. * std::log(*params.rad_target_eff));
		_builder.rad.scale_exp_est = *params.rad_scale_exp;
		_builder.rad.min_cell_explore_samples = 32;
		_builder.rad.hist_num_per_bin = 2;
		_builder.rad.max_explore_cells = *params.rad_max_cells;
		break;
	case EventType::EXCL:
		new (&_builder.excl) ExclBuilder(
			ExclDensity(params, sf),
			seed_dist(seed_rnd));
		_builder.excl.explore_progress_reporter = builder_reporters.explore_progress;
		_builder.excl.tune_progress_reporter = builder_reporters.tune_progress;
		_builder.excl.check_samples = 16384;
		_builder.excl.target_rel_var = 0.25;
		_builder.excl.scale_exp_est = 0.25;
		_builder.excl.min_cell_explore_samples = 2;
		_builder.excl.hist_num_per_bin = 2;
		_builder.excl.max_explore_cells = 16777216;
		break;
	}
}

Builder::Builder(Builder&& other) :
		_type(other._type) {
	switch (_type) {
	case EventType::NRAD:
		new (&_builder.nrad) NradBuilder(std::move(other._builder.nrad));
		break;
	case EventType::RAD:
		new (&_builder.rad) RadBuilder(std::move(other._builder.rad));
		break;
	case EventType::EXCL:
		new (&_builder.excl) ExclBuilder(std::move(other._builder.excl));
		break;
	}
}

Builder::~Builder() {
	switch (_type) {
	case EventType::NRAD:
		_builder.nrad.~NradBuilder();
		break;
	case EventType::RAD:
		_builder.rad.~RadBuilder();
		break;
	case EventType::EXCL:
		_builder.excl.~ExclBuilder();
		break;
	}
}

void Builder::explore() {
	switch (_type) {
	case EventType::NRAD:
		_builder.nrad.explore();
		break;
	case EventType::RAD:
		_builder.rad.explore();
		break;
	case EventType::EXCL:
		_builder.excl.explore();
		break;
	}
}

void Builder::tune() {
	switch (_type) {
	case EventType::NRAD:
		_builder.nrad.tune();
		break;
	case EventType::RAD:
		_builder.rad.tune();
		break;
	case EventType::EXCL:
		_builder.excl.tune();
		break;
	}
}

void Builder::write(std::ostream& os) {
	switch (_type) {
	case EventType::NRAD:
		_builder.nrad.write(os);
		break;
	case EventType::RAD:
		_builder.rad.write(os);
		break;
	case EventType::EXCL:
		_builder.excl.write(os);
		break;
	}
}

Real Builder::rel_var(Real* err_out) const {
	switch (_type) {
	case EventType::NRAD:
		return _builder.nrad.rel_var(err_out);
	case EventType::RAD:
		return _builder.rad.rel_var(err_out);
	case EventType::EXCL:
		return _builder.excl.rel_var(err_out);
	default:
		return std::numeric_limits<Real>::quiet_NaN();
	}
}

Generator::Generator(
		EventType type,
		Params const& params,
		sf::SfSet const& sf,
		std::istream& is) :
		_type(type) {
	std::minstd_rand seed_rnd(*params.seed >= 0 ? *params.seed : rnd_dev());
	std::uniform_int_distribution<Seed> seed_dist;
	switch (_type) {
	case EventType::NRAD:
		new (&_generator.nrad) NradGenerator(
			NradDensity(params, sf),
			seed_dist(seed_rnd));
		_generator.nrad.read(is);
		break;
	case EventType::RAD:
		new (&_generator.rad) RadGenerator(
			RadDensity(params, sf),
			seed_dist(seed_rnd));
		_generator.rad.read(is);
		break;
	case EventType::EXCL:
		new (&_generator.excl) ExclGenerator(
			ExclDensity(params, sf),
			seed_dist(seed_rnd));
		_generator.excl.read(is);
		break;
	}
}

Generator::Generator(Generator&& other) :
		_type(other._type),
		_weights(other._weights) {
	switch (_type) {
	case EventType::NRAD:
		new (&_generator.nrad) NradGenerator(std::move(other._generator.nrad));
		break;
	case EventType::RAD:
		new (&_generator.rad) RadGenerator(std::move(other._generator.rad));
		break;
	case EventType::EXCL:
		new (&_generator.excl) ExclGenerator(std::move(other._generator.excl));
		break;
	}
}

Generator::~Generator() {
	switch (_type) {
	case EventType::NRAD:
		_generator.nrad.~NradGenerator();
		break;
	case EventType::RAD:
		_generator.rad.~RadGenerator();
		break;
	case EventType::EXCL:
		_generator.excl.~ExclGenerator();
		break;
	}
}

Real Generator::generate(Real* out) {
	Real weight;
	PointFull ph_space;
	switch (_type) {
	case EventType::NRAD:
		{
			Point<6> vec;
			_generator.nrad.generate(&weight, &vec);
			_generator.nrad.func().transform(vec, &ph_space);
		}
		break;
	case EventType::RAD:
		{
			Point<9> vec;
			_generator.rad.generate(&weight, &vec);
			_generator.rad.func().transform(vec, &ph_space);
		}
		break;
	case EventType::EXCL:
		{
			Point<8> vec;
			_generator.excl.generate(&weight, &vec);
			_generator.excl.func().transform(vec, &ph_space);
		}
		break;
	}
	std::copy(ph_space.begin(), ph_space.end(), out);
	_weights += weight;
	return weight;
}

Real Generator::prime() const {
	switch (_type) {
	case EventType::NRAD:
		return _generator.nrad.prime();
	case EventType::RAD:
		return _generator.rad.prime();
	case EventType::EXCL:
		return _generator.excl.prime();
	default:
		return std::numeric_limits<Real>::quiet_NaN();
	}
}

