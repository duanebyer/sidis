#include "generator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
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

unsigned rand_simple(unsigned prev, std::size_t steps) {
	for (std::size_t idx = 0; idx < steps; ++idx) {
		prev = (1664525 * prev + 1013904223) & 0xFFFF;
	}
	return prev;
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

Double NradDensity::transform(Point<6> const& unit_vec, Point<6>* ph_vec) const noexcept {
	Double jacobian;
	kin::Kinematics kin;
	if (!cut::take(_cut, _ps, _S, unit_vec.data(), &kin, &jacobian)) {
		jacobian = 0.;
	}
	(*ph_vec)[0] = kin.x;
	(*ph_vec)[1] = kin.y;
	(*ph_vec)[2] = kin.z;
	(*ph_vec)[3] = kin.ph_t_sq;
	(*ph_vec)[4] = kin.phi_h;
	(*ph_vec)[5] = kin.phi;
	return jacobian;
}

Double NradDensity::operator()(Point<6> const& vec) const noexcept {
	kin::Kinematics kin;
	Double jacobian;
	if (!cut::take(_cut, _ps, _S, vec.data(), &kin, &jacobian)) {
		return 0.;
	}
	math::Vec3 eta = frame::hadron_from_target(kin) * _target_pol;
	// TODO: Evaluate when it is a good approximation to say that
	// `nrad ~ nrad_ir`. This happens because for small `k_0_bar`, the
	// contribution of `rad_f` integrated up to `k_0_bar` becomes vanishingly
	// small, so it can be neglected. However, this must be balanced with
	// choosing `k_0_bar` to be non-zero to avoid the infrared divergence in the
	// radiative part of the cross-section.
	Double xs;
	switch (_rc_method) {
	case RcMethod::NONE:
		xs = xs::born(kin, _sf, _beam_pol, eta);
		break;
	case RcMethod::APPROX:
		xs = xs::nrad_ir(kin, _sf, _beam_pol, eta, _soft_threshold);
		break;
	case RcMethod::EXACT:
		xs = xs::nrad_integ(kin, _sf, _beam_pol, eta, _soft_threshold).val;
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
	// Validate that the `k_0_bar` maximum is above the soft threshold.
	math::Bound k_0_bar = params["cut.k_0_bar"].any();
	if (!(k_0_bar.max() > _soft_threshold)) {
		throw std::runtime_error(
			"Parameter 'cut.k_0_bar' does not encompass any of the radiative "
			"range above 'phys.soft_threshold'.");
	}
}

Double RadDensity::transform(Point<9> const& unit_vec, Point<9>* ph_vec) const noexcept {
	Double jacobian;
	kin::KinematicsRad kin;
	if (!cut::take(_cut, _cut_rad, _ps, _S, unit_vec.data(), &kin, &jacobian)) {
		jacobian = 0.;
	}
	(*ph_vec)[0] = kin.x;
	(*ph_vec)[1] = kin.y;
	(*ph_vec)[2] = kin.z;
	(*ph_vec)[3] = kin.ph_t_sq;
	(*ph_vec)[4] = kin.phi_h;
	(*ph_vec)[5] = kin.phi;
	(*ph_vec)[6] = kin.tau;
	(*ph_vec)[7] = kin.phi_k;
	(*ph_vec)[8] = kin.R;
	return jacobian;
}

Double RadDensity::operator()(Point<9> const& vec) const noexcept {
	kin::KinematicsRad kin_rad;
	Double jacobian;
	if (!cut::take(_cut, _cut_rad, _ps, _S, vec.data(), &kin_rad, &jacobian)) {
		return 0.;
	}
	kin::Kinematics kin = kin_rad.project();
	math::Vec3 eta = frame::hadron_from_target(kin) * _target_pol;
	Double xs = xs::rad(kin_rad, _sf, _beam_pol, eta);
	if (!std::isfinite(xs) || !(xs >= 0.)) {
		return 0.;
	} else {
		return jacobian * xs;
	}
}

Builder::Builder(
		EventType ev_type,
		BuilderReporters const& reporters,
		Params& params,
		sf::SfSet const& sf) :
		_ev_type(ev_type) {
	SeedInit seed_init = params[p_name_init_seed(_ev_type)].any();
	if (seed_init.any) {
		throw std::runtime_error("Must choose specific seed for Builder.");
	}
	_seed = seed_init.seed;
	// Mix up the seed based on event type, to avoid duplicate seeds for
	// different event types.
	std::minstd_rand seed_rnd(rand_simple(_seed, static_cast<int>(_ev_type)));
	std::uniform_int_distribution<Seed> seed_dist;

	// TODO: Double check a lot of these dynamic casts, they could break on
	// other platforms with differently sized integers/doubles. This might get
	// fixed when we refactor the generator build parameter extraction process.
	bubble::CellBuilderParams<Double> builder_params = {};
	builder_params.check_samples = 16384;
	builder_params.target_rel_var = std::expm1(
		-2. * std::log(params[p_name_init_target_eff(_ev_type)].any().as<Double>()));
	builder_params.scale_exp_est = params[p_name_init_scale_exp(_ev_type)].any();
	builder_params.min_cell_explore_samples = 512;
	builder_params.hist_num_per_bin = 2;
	builder_params.max_explore_cells = params[p_name_init_max_cells(_ev_type)].any();

	switch (_ev_type) {
	case EventType::NRAD:
		new (&_builder.nrad) NradBuilder(
			NradDensity(params, sf),
			seed_dist(seed_rnd));
		_builder.nrad.par = builder_params;
		_builder.nrad.explore_progress_reporter = reporters.explore_progress;
		_builder.nrad.tune_progress_reporter = reporters.tune_progress;
		break;
	case EventType::RAD:
		new (&_builder.rad) RadBuilder(
			RadDensity(params, sf),
			seed_dist(seed_rnd));
		_builder.rad.par = builder_params;
		_builder.rad.explore_progress_reporter = reporters.explore_progress;
		_builder.rad.tune_progress_reporter = reporters.tune_progress;
		break;
	case EventType::EXCL:
		new (&_builder.excl) ExclBuilder(
			ExclDensity(params, sf),
			seed_dist(seed_rnd));
		_builder.excl.par = builder_params;
		_builder.excl.explore_progress_reporter = reporters.explore_progress;
		_builder.excl.tune_progress_reporter = reporters.tune_progress;
		break;
	}
}

Builder::Builder(Builder&& other) :
		_ev_type(other._ev_type) {
	switch (_ev_type) {
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
	switch (_ev_type) {
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
	switch (_ev_type) {
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
	switch (_ev_type) {
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
	switch (_ev_type) {
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

Double Builder::rel_var(Double* err_out) const {
	switch (_ev_type) {
	case EventType::NRAD:
		return _builder.nrad.rel_var(err_out);
	case EventType::RAD:
		return _builder.rad.rel_var(err_out);
	case EventType::EXCL:
		return _builder.excl.rel_var(err_out);
	default:
		return std::numeric_limits<Double>::quiet_NaN();
	}
}

std::size_t Builder::size() const {
	switch (_ev_type) {
	case EventType::NRAD:
		return _builder.nrad.tree().size();
	case EventType::RAD:
		return _builder.rad.tree().size();
	case EventType::EXCL:
		return _builder.excl.tree().size();
	default:
		return 0;
	}
}

Generator::Generator(
		EventType ev_type,
		Params& params,
		sf::SfSet const& sf,
		std::istream& is) :
		_ev_type(ev_type),
		_seed(),
		_rej_scale(0.),
		_weights(),
		_count(0),
		_count_acc(0) {
	SeedGen seed_gen = params["mc.seed"].any();
	if (seed_gen.any || seed_gen.seeds.size() != 1) {
		throw std::runtime_error("Must choose specific seed for Generator.");
	}
	_seed = *seed_gen.seeds.begin();
	std::minstd_rand seed_rnd(rand_simple(_seed, static_cast<int>(_ev_type)));
	std::uniform_int_distribution<Seed> seed_dist;

	// TODO: Double check this dynamic cast as well.
	_rej_scale = params[p_name_gen_rej_scale(_ev_type)].any();

	switch (_ev_type) {
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
		_ev_type(other._ev_type),
		_seed(other._seed),
		_rej_scale(other._rej_scale),
		_weights(other._weights),
		_count(other._count),
		_count_acc(other._count_acc) {
	switch (_ev_type) {
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
	switch (_ev_type) {
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

Double Generator::generate(Double* ph_out, Double* unit_out) {
	Double weight = 0.;
	while (weight == 0.) {
		switch (_ev_type) {
		case EventType::NRAD:
			{
				Point<6> unit_vec, ph_vec;
				_generator.nrad.generate_rej(_rej_scale, &weight, &unit_vec);
				_generator.nrad.func().transform(unit_vec, &ph_vec);
				std::copy(ph_vec.begin(), ph_vec.end(), ph_out);
				std::copy(unit_vec.begin(), unit_vec.end(), unit_out);
			}
			break;
		case EventType::RAD:
			{
				Point<9> unit_vec, ph_vec;
				_generator.rad.generate_rej(_rej_scale, &weight, &unit_vec);
				_generator.rad.func().transform(unit_vec, &ph_vec);
				std::copy(ph_vec.begin(), ph_vec.end(), ph_out);
				std::copy(unit_vec.begin(), unit_vec.end(), unit_out);
			}
			break;
		case EventType::EXCL:
			{
				Point<8> unit_vec, ph_vec;
				_generator.excl.generate_rej(_rej_scale, &weight, &unit_vec);
				_generator.excl.func().transform(unit_vec, &ph_vec);
				std::copy(ph_vec.begin(), ph_vec.end(), ph_out);
				std::copy(unit_vec.begin(), unit_vec.end(), unit_out);
			}
			break;
		}
		_weights += weight;
		_count += 1;
	}
	_count_acc += 1;
	return weight;
}

Double Generator::prime() const {
	switch (_ev_type) {
	case EventType::NRAD:
		return _generator.nrad.prime();
	case EventType::RAD:
		return _generator.rad.prime();
	case EventType::EXCL:
		return _generator.excl.prime();
	default:
		return std::numeric_limits<Double>::quiet_NaN();
	}
}

