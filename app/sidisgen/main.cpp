#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>

#include <TBranch.h>
#include <TFile.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TTree.h>

#include <sidis/sidis.hpp>
#include <sidis/extra/bounds.hpp>
#include <sidis/extra/transform.hpp>
#include <sidis/extra/vector.hpp>
#include <sidis/sf_model/ww.hpp>

using namespace sidis;

namespace {

sf::model::WW const ww;
constant::Hadron const hadron = constant::Hadron::PI_P;
Real const M_th = constant::MASS_P + constant::MASS_PI_0;
Real const PI = constant::PI;

// Converts between the `sidis` 4-vector type and the `ROOT` 4-vector type.
TLorentzVector convert_vec4(math::Vec4 v) {
	return TLorentzVector(v.x, v.y, v.z, v.t);
}

struct XsNRad : public TFoamIntegrand {
	kin::Initial const initial_state;
	Real const beam_pol;
	math::Vec3 const target_pol;
	Real const k0_cut;
	math::Bounds const x_cut;
	math::Bounds const Q_sq_cut;

	XsNRad(
		kin::Initial initial_state,
		Real beam_pol,
		math::Vec3 target_pol,
		Real k0_cut,
		Real x_max=0.1,
		Real Q_sq_min=1.) :
		initial_state(initial_state),
		beam_pol(beam_pol),
		target_pol(target_pol),
		k0_cut(k0_cut),
		x_cut(0., x_max),
		Q_sq_cut(Q_sq_min, 1000.) { }

	// Converts from the unit hyper-cube to a point in phase space, also
	// providing the associated jacobian.
	kin::Kinematics GetKinematics(Double_t const* vec, Real* jacobian) {
		Real S = 2.*math::dot(initial_state.p, initial_state.k1);
		math::Bounds x_bounds = kin::x_bounds(initial_state) & x_cut;
		Real x = x_bounds.lerp(vec[0]);
		math::Bounds y_cut(Q_sq_cut.min / (S * x), Q_sq_cut.max / (S * x));
		math::Bounds y_bounds = kin::y_bounds(initial_state, x) & y_cut;
		Real y = y_bounds.lerp(vec[1]);
		math::Bounds z_bounds = kin::z_bounds(initial_state, hadron, M_th, x, y);
		Real z = z_bounds.lerp(vec[2]);
		math::Bounds ph_t_sq_bounds = kin::ph_t_sq_bounds(initial_state, hadron, M_th, x, y, z);
		Real ph_t_sq = ph_t_sq_bounds.lerp(vec[3]);
		math::Bounds phi_h_bounds = math::Bounds(0., 2. * PI);
		Real phi_h = phi_h_bounds.lerp(vec[4]);
		math::Bounds phi_bounds = math::Bounds(0., 2. * PI);
		Real phi = phi_bounds.lerp(vec[5]);
		*jacobian = x_bounds.size() * y_bounds.size() * z_bounds.size()
			* ph_t_sq_bounds.size() * phi_h_bounds.size() * phi_bounds.size();

		kin::PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, phi };
		kin::Kinematics kin(initial_state, phase_space, hadron, M_th);
		return kin;
	}

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 6) {
			return 0.;
		}
		Real jacobian;
		kin::Kinematics kin = GetKinematics(vec, &jacobian);
		math::Vec3 eta = frame::hadron_from_target(kin) * target_pol;
		// TODO: Evaluate when it is a good approximation to say that
		// `nrad ~ nrad_ir`. This happens because for small `k0_cut`, the
		// contribution of `rad_f` integrated up to `k0_cut` becomes vanishingly
		// small, so it can be neglected. However, this must be balanced with
		// choosing `k0_cut` to be non-zero to avoid the infrared divergence in
		// the radiative part of the cross-section.
		Real xs = xs::nrad_ir(beam_pol, eta, kin, ww, k0_cut);
		// Some kinematic regions will be out of range for the structure
		// functions, so return 0 in those cases.
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

struct XsRad : public TFoamIntegrand {
	kin::Initial const initial_state;
	Real const beam_pol;
	math::Vec3 const target_pol;
	Real const k0_cut;
	math::Bounds const x_cut;
	math::Bounds const Q_sq_cut;

	XsRad(
		kin::Initial initial_state,
		Real beam_pol,
		math::Vec3 target_pol,
		Real k0_cut,
		Real x_max=0.1,
		Real Q_sq_min=1.) :
		initial_state(initial_state),
		beam_pol(beam_pol),
		target_pol(target_pol),
		k0_cut(k0_cut),
		x_cut(0., x_max),
		Q_sq_cut(Q_sq_min, 1000.) { }

	kin::KinematicsRad GetKinematics(Double_t const* vec, Real* jacobian) {
		Real S = 2.*math::dot(initial_state.p, initial_state.k1);
		math::Bounds x_bounds = kin::x_bounds(initial_state) & x_cut;
		Real x = x_bounds.lerp(vec[0]);
		math::Bounds y_cut(Q_sq_cut.min / (S * x), Q_sq_cut.max / (S * x));
		math::Bounds y_bounds = kin::y_bounds(initial_state, x) & y_cut;
		Real y = y_bounds.lerp(vec[1]);
		math::Bounds z_bounds = kin::z_bounds(initial_state, hadron, M_th, x, y);
		Real z = z_bounds.lerp(vec[2]);
		math::Bounds ph_t_sq_bounds = kin::ph_t_sq_bounds(initial_state, hadron, M_th, x, y, z);
		Real ph_t_sq = ph_t_sq_bounds.lerp(vec[3]);
		math::Bounds phi_h_bounds = math::Bounds(0., 2. * PI);
		Real phi_h = phi_h_bounds.lerp(vec[4]);
		math::Bounds phi_bounds = math::Bounds(0., 2. * PI);
		Real phi = phi_bounds.lerp(vec[5]);

		kin::PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, phi };
		kin::Kinematics kin(initial_state, phase_space, hadron, M_th);

		math::Bounds tau_bounds = kin::tau_bounds(kin);
		Real tau = tau_bounds.lerp(vec[6]);
		math::Bounds phi_k_bounds = math::Bounds(0., 2. * PI);
		Real phi_k = phi_k_bounds.lerp(vec[7]);
		math::Bounds R_bounds = kin::R_bounds_hard(kin, tau, phi_k, k0_cut);
		Real R = R_bounds.lerp(vec[8]);
		*jacobian = x_bounds.size() * y_bounds.size() * z_bounds.size()
			* ph_t_sq_bounds.size() * phi_h_bounds.size() * phi_bounds.size()
			* tau_bounds.size() * phi_k_bounds.size() * R_bounds.size();

		kin::KinematicsRad kin_rad(kin, tau, phi_k, R);
		return kin_rad;
	}

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 9) {
			return 0.;
		}
		Real jacobian;
		kin::KinematicsRad kin_rad = GetKinematics(vec, &jacobian);
		math::Vec3 eta = frame::hadron_from_target(kin_rad.project()) * target_pol;
		Real xs = xs::rad(beam_pol, eta, kin_rad, ww);
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

// Convenience function for drawing a progress bar in the terminal.
bool write_progress_bar(std::ostream& os, double fraction, unsigned width=70) {
	if (fraction < 0. || std::isnan(fraction)) {
		fraction = 0.;
	} else if (fraction > 1.) {
		fraction = 1.;
	}
	int percentage = static_cast<int>(fraction * 100.);
	std::string edge_left_str = "[";
	std::string edge_right_str = "]";
	std::string spacing_str = " ";
	std::stringstream percentage_ss;
	percentage_ss << std::setw(3) << percentage << '%';
	std::string percentage_str = percentage_ss.str();
	unsigned percentage_width = percentage_str.size();
	unsigned edge_width = edge_left_str.size() + edge_right_str.size();
	unsigned space_width = spacing_str.size();
	if (width <= percentage_width + edge_width + space_width) {
		return false;
	}
	unsigned bar_width = width - percentage_width - edge_width - space_width;
	unsigned bar_filled = static_cast<unsigned>(fraction * bar_width);
	os << edge_left_str;
	for (unsigned idx = 0; idx < bar_width; ++idx) {
		if (idx < bar_filled) {
			os << '=';
		} else {
			os << ' ';
		}
	}
	os << edge_right_str;
	os << spacing_str;
	os << percentage_str;
	return os.good();
}

}

int main(int argc, char** argv) {
	Real k0_cut = 0.001;
	Real beam_energy = 11.;
	Real beam_pol = 0.;
	constant::Nucleus target = constant::Nucleus::P;
	constant::Lepton beam = constant::Lepton::E;
	math::Vec3 target_pol(0., 0., 0.);
	kin::Initial initial_state(target, beam, beam_energy);
	unsigned N_init_nrad = 1000;
	unsigned N_init_rad  = 1000;
	unsigned N_gen = 100000;
	// These keep track of how many radiative/non-radiative events have been
	// generated so far, so that the ratio of radiative to non-radiative events
	// generated can be kept roughly the same as the ratio of the total
	// radiative/non-radiative cross-sections.
	unsigned N_gen_nrad = 0;
	unsigned N_gen_rad = 0;

	TFile file_out("gen.root", "Recreate");
	TRandom3 random(0);

	std::cout << "Non-radiative FOAM initialization." << std::endl;
	TFoam foam_nrad("FoamNRad");
	XsNRad xs_nrad(initial_state, beam_pol, target_pol, k0_cut);
	foam_nrad.SetChat(0);
	foam_nrad.SetkDim(6);
	foam_nrad.SetRho(&xs_nrad);
	foam_nrad.SetPseRan(&random);
	foam_nrad.SetnSampl(N_init_nrad);
	foam_nrad.Initialize();
	foam_nrad.Write();

	std::cout << "Radiative FOAM initialization." << std::endl;
	TFoam foam_rad("FoamRad");
	XsRad xs_rad(initial_state, beam_pol, target_pol, k0_cut);
	foam_rad.SetChat(0);
	foam_rad.SetkDim(9);
	foam_rad.SetRho(&xs_rad);
	foam_rad.SetPseRan(&random);
	foam_rad.SetnSampl(N_init_rad);
	foam_rad.Initialize();
	foam_rad.Write();

	TTree events("Events", "Events");
	Bool_t is_rad;
	Double_t weight;
	Double_t x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R;
	Double_t Q_sq;
	TLorentzVector p, k1, q, k2, ph, k;
	events.Branch("is_rad", &is_rad);
	events.Branch("weight", &weight);
	events.Branch("x", &x);
	events.Branch("y", &y);
	events.Branch("z", &z);
	events.Branch("ph_t_sq", &ph_t_sq);
	events.Branch("phi_h", &phi_h);
	events.Branch("phi", &phi);
	events.Branch("tau", &tau);
	events.Branch("phi_k", &phi_k);
	events.Branch("R", &R);
	events.Branch("Q_sq", &Q_sq);
	events.Branch("p", "TLorentzVector", &p);
	events.Branch("k1", "TLorentzVector", &k1);
	events.Branch("q", "TLorentzVector", &q);
	events.Branch("k2", "TLorentzVector", &k2);
	events.Branch("ph", "TLorentzVector", &ph);
	events.Branch("k", "TLorentzVector", &k);

	std::cout << "Writing parameters." << std::endl;
	TParameter<int> p_target("target_pid", static_cast<int>(target));
	TParameter<int> p_beam("beam_pid", static_cast<int>(beam));
	TParameter<Double_t> p_beam_energy("beam_energy", beam_energy);
	p_target.Write();
	p_beam.Write();

	std::cout << "Generating events." << std::endl;
	Double_t total_nrad, total_nrad_err;
	Double_t total_rad, total_rad_err;
	for (unsigned event_idx = 0; event_idx < N_gen; ++event_idx) {
		if (event_idx % (N_gen / 100) == 0) {
			write_progress_bar(std::cout, static_cast<double>(event_idx) / N_gen);
			std::cout << '\r';
			std::cout.flush();
		}
		// Estimate the total radiative and non-radiative cross-sections and
		// generate a radiative/non-radiative event accordingly. The total
		// cross-section estimates are improved as more events are generated.
		foam_nrad.GetIntegMC(total_nrad, total_nrad_err);
		foam_rad.GetIntegMC(total_rad, total_rad_err);
		if (!std::isfinite(total_nrad) || total_nrad == 0.) {
			total_nrad = 0.;
			total_nrad_err = std::numeric_limits<Double_t>::max();
		}
		if (!std::isfinite(total_rad) || total_rad == 0.) {
			total_rad = 0.;
			total_rad_err = std::numeric_limits<Double_t>::max();
		}
		bool choose_nrad;
		if (event_idx == 0) {
			// On the first event, we don't know anything about the total cross-
			// sections, so always generate a non-radiative event (since the
			// non-radiative cross-section is guaranteed non-zero).
			choose_nrad = true;
		} else {
			Double_t total_nrad_max = total_nrad + total_nrad_err;
			Double_t total_rad_max = total_rad + total_rad_err;
			Double_t target_ratio = total_nrad_max
				/ (total_nrad_max + total_rad_max);
			Double_t ratio = static_cast<Double_t>(N_gen_nrad)
				/ (N_gen_nrad + N_gen_rad);
			// The reason we choose whether to generate a radiative event in
			// this way is because if the total radiative cross-section has a
			// large uncertainty, it will cause radiative events to be generated
			// more than they should be, which will in turn reduce the
			// uncertainty in the total radiative cross-section to a reasonable
			// size.
			choose_nrad = (ratio < target_ratio);
		}

		if (choose_nrad) {
			// Generate a non-radiative event.
			is_rad = false;
			N_gen_nrad += 1;
			Double_t event_vec[6];
			weight = foam_nrad.MCgenerate(event_vec);
			Real jacobian;
			kin::Kinematics kin = xs_nrad.GetKinematics(event_vec, &jacobian);
			kin::Final final_state(initial_state, target_pol, kin);
			// Fill in branches.
			x = kin.x;
			y = kin.y;
			z = kin.z;
			ph_t_sq = kin.ph_t_sq;
			phi_h = kin.phi_h;
			phi = kin.phi;
			Q_sq = kin.Q_sq;
			p = convert_vec4(initial_state.p);
			k1 = convert_vec4(initial_state.k1);
			q = convert_vec4(final_state.q);
			k2 = convert_vec4(final_state.k2);
			ph = convert_vec4(final_state.ph);
		} else {
			// Generate a radiative event.
			is_rad = true;
			N_gen_rad += 1;
			Double_t event_vec[9];
			weight = foam_rad.MCgenerate(event_vec);
			Real jacobian;
			kin::KinematicsRad kin = xs_rad.GetKinematics(event_vec, &jacobian);
			kin::FinalRad final_state(initial_state, target_pol, kin);
			// Fill in branches.
			x = kin.x;
			y = kin.y;
			z = kin.z;
			ph_t_sq = kin.ph_t_sq;
			phi_h = kin.phi_h;
			phi = kin.phi;
			tau = kin.tau;
			phi_k = kin.phi_k;
			R = kin.R;
			Q_sq = kin.Q_sq;
			p = convert_vec4(initial_state.p);
			k1 = convert_vec4(initial_state.k1);
			q = convert_vec4(final_state.q);
			k2 = convert_vec4(final_state.k2);
			ph = convert_vec4(final_state.ph);
			k = convert_vec4(final_state.k);
		}
		events.Fill();
	}
	write_progress_bar(std::cout, 1.);
	std::cout << std::endl;
	std::cout << "Writing events to file." << std::endl;
	events.Write();

	foam_nrad.GetIntegMC(total_nrad, total_nrad_err);
	foam_rad.GetIntegMC(total_rad, total_rad_err);

	// Write total cross-sections to file.
	TParameter<Double_t> p_xs_total_nrad("xs_total_nrad", total_nrad);
	TParameter<Double_t> p_xs_total_nrad_err("xs_total_nrad_err", total_nrad_err);
	TParameter<Double_t> p_xs_total_rad("xs_total_rad", total_rad);
	TParameter<Double_t> p_xs_total_rad_err("xs_total_rad_err", total_rad_err);
	p_xs_total_nrad.Write();
	p_xs_total_nrad_err.Write();
	p_xs_total_rad.Write();
	p_xs_total_rad_err.Write();

	file_out.Close();

	std::cout << "Statistics:" << std::endl;
	std::cout << "\tNumber of non-radiative events:    " << N_gen_nrad << std::endl;
	std::cout << "\tNumber of radiative events:        " << N_gen_rad << std::endl;
	std::cout << "\tTotal non-radiative cross-section: "
		<< std::scientific << total_nrad << " ± " << total_nrad_err << std::endl;
	std::cout << "\tTotal radiative cross-section:     "
		<< std::scientific << total_rad << " ± " << total_rad_err << std::endl;
	return 0;
}

