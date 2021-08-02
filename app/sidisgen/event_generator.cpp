#include "event_generator.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <TDirectory.h>

#include "exception.hpp"

using namespace sidis;

namespace {

char const* get_foam_file_name(EventType type, Params const& params) {
	switch (type) {
	case EventType::NRAD:
		return params.foam_nrad_file->c_str();
	case EventType::RAD:
		return params.foam_rad_file->c_str();
	case EventType::EXCL:
		throw Exception(
			ERROR_UNIMPLEMENTED,
			"Exclusive tail contribution is not available.");
	default:
		return "<error>";
	}
}

cut::Cut get_cut_from_params(Params const& params) {
	Real const DEG = PI / 180.;
	cut::Cut result;
	result.x = params.x_cut.get_or(math::Bound::INVALID);
	result.y = params.y_cut.get_or(math::Bound::INVALID);
	result.z = params.z_cut.get_or(math::Bound::INVALID);
	result.ph_t_sq = params.ph_t_sq_cut.get_or(math::Bound::INVALID);
	result.phi_h = DEG * params.phi_h_cut.get_or(math::Bound::INVALID);
	result.phi = DEG * params.phi_cut.get_or(math::Bound::INVALID);
	result.Q_sq = params.Q_sq_cut.get_or(math::Bound::INVALID);
	result.t = params.t_cut.get_or(math::Bound::INVALID);
	result.W_sq = params.W_sq_cut.get_or(math::Bound::INVALID);
	result.r = params.r_cut.get_or(math::Bound::INVALID);
	result.mx_sq = params.mx_sq_cut.get_or(math::Bound::INVALID);
	result.q_0 = params.q_0_cut.get_or(math::Bound::INVALID);
	result.k2_0 = params.k2_0_cut.get_or(math::Bound::INVALID);
	result.ph_0 = params.ph_0_cut.get_or(math::Bound::INVALID);
	result.theta_q = DEG * params.theta_q_cut.get_or(math::Bound::INVALID);
	result.theta_k2 = DEG * params.theta_k2_cut.get_or(math::Bound::INVALID);
	result.theta_h = DEG * params.theta_h_cut.get_or(math::Bound::INVALID);
	return result;
}

cut::CutRad get_cut_rad_from_params(Params const& params) {
	Real const DEG = PI / 180.;
	cut::CutRad result;
	if (*params.gen_rad) {
		result.tau = params.tau_cut.get_or(math::Bound::INVALID);
		result.phi_k = DEG * params.phi_k_cut.get_or(math::Bound::INVALID);
		result.R = params.R_cut.get_or(math::Bound::INVALID);
		// The `k_0_bar` cut is mandatory.
		result.k_0_bar = *params.k_0_bar_cut
			& math::Bound(*params.k_0_bar, INF);
		result.k_0 = params.k_0_cut.get_or(math::Bound::INVALID);
		result.theta_k = DEG * params.theta_k_cut.get_or(math::Bound::INVALID);
	}
	return result;
}

class NRadDensity : public EventDensity {
private:
	Params params;
	cut::Cut cut;
	part::Particles ps;
	Real S;
	sf::SfSet const& sf;

public:
	NRadDensity(Params const& params, sf::SfSet const& sf) :
		EventDensity(6),
		params(params),
		cut(get_cut_from_params(params)),
		ps(*params.target, *params.beam, *params.hadron, *params.Mth),
		S(2. * mass(*params.target) * *params.beam_energy),
		sf(sf) { }

	Double_t Transform(Double_t const* vec_unit, Double_t* vec_out) const override {
		Real jacobian;
		kin::Kinematics kin;
		if (!cut::take(cut, ps, S, vec_unit, &kin, &jacobian)) {
			jacobian = 0.;
		}
		vec_out[0] = kin.x;
		vec_out[1] = kin.y;
		vec_out[2] = kin.z;
		vec_out[3] = kin.ph_t_sq;
		vec_out[4] = kin.phi_h;
		vec_out[5] = kin.phi;
		return jacobian;
	}
	Double_t Density(int dim, Double_t* vec_unit) override {
		if (dim != this->dim) {
			throw std::runtime_error(
				"Evaluating integrand of dimension " + std::to_string(this->dim)
				+ " with vector of dimension " + std::to_string(dim) + ".");
		}
		kin::Kinematics kin;
		Real jacobian;
		if (!cut::take(cut, ps, S, vec_unit, &kin, &jacobian)) {
			return 0.;
		}
		math::Vec3 eta = frame::hadron_from_target(kin) * *params.target_pol;
		// TODO: Evaluate when it is a good approximation to say that
		// `nrad ~ nrad_ir`. This happens because for small `k_0_bar`, the
		// contribution of `rad_f` integrated up to `k_0_bar` becomes
		// vanishingly small, so it can be neglected. However, this must be
		// balanced with choosing `k_0_bar` to be non-zero to avoid the infrared
		// divergence in the radiative part of the cross-section.
		Real xs;
		switch (*params.rc_method) {
		case RcMethod::NONE:
			xs = xs::born(kin, sf, *params.beam_pol, eta);
			break;
		case RcMethod::APPROX:
			xs = xs::nrad_ir(kin, sf, *params.beam_pol, eta, *params.k_0_bar);
			break;
		case RcMethod::EXACT:
			xs = xs::nrad(kin, sf, *params.beam_pol, eta, *params.k_0_bar);
			break;
		default:
			throw std::runtime_error(
				"Unsupported radiative corrections method with index "
				+ std::to_string(static_cast<int>(*params.rc_method)) + ".");
		}
		// Some kinematic regions will be out of range for the structure
		// functions, so return 0 in those cases.
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

class RadDensity : public EventDensity {
private:
	Params params;
	cut::Cut cut;
	cut::CutRad cut_rad;
	part::Particles ps;
	Real S;
	sf::SfSet const& sf;

public:
	RadDensity(Params params, sf::SfSet const& sf) :
		EventDensity(9),
		params(params),
		cut(get_cut_from_params(params)),
		cut_rad(get_cut_rad_from_params(params)),
		ps(*params.target, *params.beam, *params.hadron, *params.Mth),
		S(2. * mass(*params.target) * *params.beam_energy),
		sf(sf) { }

	Double_t Transform(Double_t const* vec_unit, Double_t* vec_out) const override {
		Real jacobian;
		kin::KinematicsRad kin;
		if (!cut::take(cut, cut_rad, ps, S, vec_unit, &kin, &jacobian)) {
			jacobian = 0.;
		}
		vec_out[0] = kin.x;
		vec_out[1] = kin.y;
		vec_out[2] = kin.z;
		vec_out[3] = kin.ph_t_sq;
		vec_out[4] = kin.phi_h;
		vec_out[5] = kin.phi;
		vec_out[6] = kin.tau;
		vec_out[7] = kin.phi_k;
		vec_out[8] = kin.R;
		return jacobian;
	}
	Double_t Density(int dim, Double_t* vec_unit) override {
		if (dim != 9) {
			return 0.;
		}
		kin::KinematicsRad kin_rad;
		Real jacobian;
		if (!cut::take(cut, cut_rad, ps, S, vec_unit, &kin_rad, &jacobian)) {
			return 0.;
		}
		kin::Kinematics kin = kin_rad.project();
		math::Vec3 eta = frame::hadron_from_target(kin) * *params.target_pol;
		Real xs = xs::rad(kin_rad, sf, *params.beam_pol, eta);
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

}

EventDensity::EventDensity(int dim) : dim(dim) {
	if (dim <= 0) {
		throw std::runtime_error("Cannot have integrand dimension <= 0.");
	}
}

EventGenerator::EventGenerator(
	EventType type,
	std::unique_ptr<TFile> foam_file,
	std::unique_ptr<EventDensity> density,
	TRandom3* random,
	TFoam* foam) :
	_type(type),
	_foam_file(std::move(foam_file)),
	_density(std::move(density)),
	_random(random),
	_foam(foam),
	_stats(),
	_vec_unit(_density->dim, 0.) { }

std::unique_ptr<EventDensity> EventGenerator::alloc_density(
		EventType type,
		Params const& params,
		sf::SfSet const& sf) {
	switch (type) {
	case EventType::NRAD:
		return std::unique_ptr<EventDensity>(new NRadDensity(params, sf));
	case EventType::RAD:
		return std::unique_ptr<EventDensity>(new RadDensity(params, sf));
	case EventType::EXCL:
		throw Exception(
			ERROR_UNIMPLEMENTED,
			"Exclusive tail contribution is not available.");
	default:
		throw std::runtime_error(
			"Unrecognized event type with index "
			+ std::to_string(static_cast<int>(type)) + ".");
	}
}

EventGenerator EventGenerator::write(
		EventType type,
		Params const& params,
		sf::SfSet const& sf,
		TRandom3* random) {
	if (!params.valid()) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			std::string("Invalid parameters."));
	}
	char const* file_name = get_foam_file_name(type, params);
	std::string foam_name = std::string("foam_") + event_type_name(type);
	ULong_t num_init  = *params.num_init >= 1 ? *params.num_init : 1;
	ULong_t num_cells  = *params.num_cells >= 1 ? *params.num_cells : 1;

	std::unique_ptr<TFile> file(new TFile(file_name, "RECREATE"));
	if (file->IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_CREATED,
			std::string("Couldn't open or create file '") + file_name + "'.");
	}
	params.write_root(*file);
	file->cd();

	std::unique_ptr<EventDensity> density =
		EventGenerator::alloc_density(type, params, sf);
	Bool_t opt_rej = *params.rej_weight >= 1.;
	TFoam* foam = new TFoam(foam_name.c_str());
	foam->SetChat(0);
	foam->SetkDim(density->dim);
	foam->SetRho(density.get());
	foam->SetPseRan(random);
	foam->SetnSampl(num_init);
	foam->SetnCells(num_cells);
	foam->SetOptRej(opt_rej ? 1 : 0);
	// Use variance reduction because it typically does a better job improving
	// the efficiency, as well as because of a bug in the maximum weight
	// reduction that causes a crash for a large number of cells.
	foam->SetOptDrive(1);
	if (opt_rej) {
		foam->SetMaxWtRej(*params.rej_weight);
	}
	foam->Initialize();
	foam->Write(foam_name.c_str());

	return EventGenerator(
		type,
		std::move(file),
		std::move(density),
		random,
		foam);
}

EventGenerator EventGenerator::read(
		EventType type,
		Params const& params,
		sidis::sf::SfSet const& sf,
		TRandom3* random) {
	Params foam_params;
	char const* file_name = get_foam_file_name(type, params);
	std::string foam_name = std::string("foam_") + event_type_name(type);
	std::unique_ptr<TFile> file(new TFile(file_name, "OPEN"));
	if (file->IsZombie()) {
		throw Exception(
			ERROR_FILE_NOT_FOUND,
			std::string("Couldn't open file '") + file_name + "'.");
	}
	foam_params.read_root(*file);
	if (!foam_params.valid()) {
		throw Exception(
			ERROR_PARAMS_INVALID,
			std::string("Invalid FOAM parameter in file '") + file_name + "'.");
	}
	params.compatible_with_foam(type, foam_params);
	file->cd();

	TFoam* foam = file->Get<TFoam>(foam_name.c_str());
	if (foam == nullptr) {
		throw Exception(
			ERROR_FOAM_NOT_FOUND,
			std::string("Failed to load FOAM from file '") + file_name + "'.");
	}

	std::unique_ptr<EventDensity> density =
		EventGenerator::alloc_density(type, params, sf);
	Bool_t opt_rej = *params.rej_weight >= 1.;
	foam->SetPseRan(random);
	foam->ResetRho(density.get());
	foam->SetOptRej(opt_rej ? 1 : 0);
	if (opt_rej) {
		foam->SetMaxWtRej(*params.rej_weight);
	}

	return EventGenerator(
		type,
		std::move(file),
		std::move(density),
		random,
		foam);
}

Double_t EventGenerator::generate(Double_t* vec) {
	Double_t weight = _foam->MCgenerate(_vec_unit.data());
	_density->Transform(_vec_unit.data(), vec);
	_stats.weight_total += weight;
	_stats.weight_sq_total += weight * weight;
	_stats.weight_p3_total += weight * weight * weight;
	_stats.weight_p4_total += weight * weight * weight * weight;
	_stats.num_events += 1;
	_foam->GetIntegMC(_stats.xs, _stats.xs_err);
	if (!std::isfinite(_stats.xs) || _stats.xs == 0.) {
		_stats.xs = 0.;
		_stats.xs_err = std::numeric_limits<Double_t>::infinity();
	}
	if (!std::isfinite(_stats.xs_err)) {
		_stats.xs_err = std::numeric_limits<Double_t>::infinity();
	}
	return weight;
}

