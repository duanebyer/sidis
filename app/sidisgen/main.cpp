#include <cmath>
#include <iostream>

#include <TBranch.h>
#include <TFile.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TLorentzVector.h>
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
	math::Bounds const x_cut;
	math::Bounds const Q_sq_cut;

	XsNRad(
		kin::Initial initial_state,
		Real beam_pol,
		math::Vec3 target_pol,
		Real x_max=0.1,
		Real Q_sq_min=1.) :
		initial_state(initial_state),
		beam_pol(beam_pol),
		target_pol(target_pol),
		x_cut(0., x_max),
		Q_sq_cut(Q_sq_min, 1000.) { }

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 6) {
			return 0.;
		}
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
		Real jacobian = x_bounds.size() * y_bounds.size() * z_bounds.size()
			* ph_t_sq_bounds.size() * phi_h_bounds.size() * phi_bounds.size();

		kin::PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, phi };
		kin::Kinematics kin(initial_state, phase_space, hadron, M_th);
		math::Vec3 eta = frame::hadron_from_target(kin) * target_pol;
		Real xs = xs::nrad(beam_pol, eta, kin, ww);
		// Some kinematic regions will be out of range for the structure
		// functions, so return 0 in those cases.
		if (!std::isfinite(xs)) {
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
	math::Bounds const x_cut;
	math::Bounds const Q_sq_cut;

	XsRad(
		kin::Initial initial_state,
		Real beam_pol,
		math::Vec3 target_pol,
		Real x_max=0.1,
		Real Q_sq_min=1.) :
		initial_state(initial_state),
		beam_pol(beam_pol),
		target_pol(target_pol),
		x_cut(0., x_max),
		Q_sq_cut(Q_sq_min, 1000.) { }

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 9) {
			return 0.;
		}
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
		math::Bounds tau_bounds = kin::tau_bounds(initial_state, x, y);
		Real tau = tau_bounds.lerp(vec[6]);
		math::Bounds phi_k_bounds = math::Bounds(0., 2. * PI);
		Real phi_k = phi_k_bounds.lerp(vec[7]);
		math::Bounds R_bounds = kin::R_bounds(initial_state, hadron, M_th, x, y, z, ph_t_sq, phi_h, tau, phi_k);
		Real R = R_bounds.lerp(vec[8]);
		Real jacobian = x_bounds.size() * y_bounds.size() * z_bounds.size()
			* ph_t_sq_bounds.size() * phi_h_bounds.size() * phi_bounds.size()
			* tau_bounds.size() * phi_k_bounds.size() * R_bounds.size();

		kin::PhaseSpaceRad phase_space { x, y, z, ph_t_sq, phi_h, phi, tau, phi_k, R };
		kin::KinematicsRad kin(initial_state, phase_space, hadron, M_th);
		math::Vec3 eta = frame::hadron_from_target(kin.project()) * target_pol;
		Real xs = xs::rad(beam_pol, eta, kin, ww);
		if (!std::isfinite(xs)) {
			return 0.;
		} else {
			// TODO: Fit this with a negative sign until we understand why the
			// radiative cross-section is returning negative values.
			return -jacobian * xs;
		}
	}
};

}

int main(int argc, char** argv) {
	Real E_b = 11.;
	Real beam_pol = 0.;
	math::Vec3 target_pol(0., 0., 0.);
	kin::Initial initial_state(constant::Nucleus::P, constant::Lepton::E, E_b);
	unsigned N_init_nrad = 1000;
	unsigned N_init_rad  = 1000;
	unsigned N_gen = 100000;

	TFile file_out("gen.root", "Recreate");
	TRandom3 random(0);

	std::cout << "Initializing non-radiative foam...";
	TFoam foam_nrad("FoamNRad");
	XsNRad xs_nrad(initial_state, beam_pol, target_pol);
	foam_nrad.SetkDim(6);
	foam_nrad.SetRho(&xs_nrad);
	foam_nrad.SetPseRan(&random);
	foam_nrad.SetnSampl(N_init_nrad);
	foam_nrad.Initialize();
	foam_nrad.Write();
	std::cout << "complete!" << std::endl;
	Double_t total_nrad, total_err_nrad;
	// TODO: The total integral isn't calculated until some events have been
	// generated.
	foam_nrad.GetIntegMC(total_nrad, total_err_nrad);
	std::cout << "σ_nrad = "
		<< std::scientific << total_nrad << "GeV⁻² ± "
		<< std::fixed << 100. * (total_err_nrad / total_nrad) << "%" << std::endl;

	std::cout << "Initializing radiative foam...";
	TFoam foam_rad("FoamRad");
	XsRad xs_rad(initial_state, beam_pol, target_pol);
	foam_rad.SetkDim(9);
	foam_rad.SetRho(&xs_rad);
	foam_rad.SetPseRan(&random);
	foam_rad.SetnSampl(N_init_rad);
	foam_rad.Initialize();
	foam_rad.Write();
	std::cout << "complete!" << std::endl;
	Double_t total_rad, total_err_rad;
	// TODO: The total integral isn't calculated until some events have been
	// generated.
	foam_rad.GetIntegMC(total_rad, total_err_rad);
	std::cout << "σ_rad = "
		<< std::scientific << total_rad << "GeV⁻² ± "
		<< std::fixed << 100. * (total_err_rad / total_rad) << "%" << std::endl;

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

	std::cout << "Generating events." << std::endl;
	for (unsigned event_idx = 0; event_idx < N_gen; ++event_idx) {
		if (random.Uniform(0., total_rad + total_nrad) < total_nrad) {
			// Generate a non-radiative event.
			is_rad = false;
			Double_t event_vec[6];
			weight = foam_nrad.MCgenerate(event_vec);
			kin::PhaseSpace phase_space {
				event_vec[0], event_vec[1], event_vec[2],
				event_vec[3], event_vec[4], event_vec[5],
			};
			kin::Kinematics kin(initial_state, phase_space, hadron, M_th);
			kin::Final final_state(initial_state, target_pol, kin);
			// Fill in branches.
			x = phase_space.x;
			y = phase_space.y;
			z = phase_space.z;
			ph_t_sq = phase_space.ph_t_sq;
			phi_h = phase_space.phi_h;
			phi = phase_space.phi;
			Q_sq = kin.Q_sq;
			p = convert_vec4(initial_state.p);
			k1 = convert_vec4(initial_state.k1);
			q = convert_vec4(final_state.q);
			k2 = convert_vec4(final_state.k2);
			ph = convert_vec4(final_state.ph);
		} else {
			// Generate a radiative event.
			is_rad = true;
			Double_t event_vec[9];
			weight = foam_nrad.MCgenerate(event_vec);
			kin::PhaseSpaceRad phase_space {
				event_vec[0], event_vec[1], event_vec[2],
				event_vec[3], event_vec[4], event_vec[5],
				event_vec[6], event_vec[7], event_vec[8],
			};
			kin::KinematicsRad kin(initial_state, phase_space, hadron, M_th);
			kin::FinalRad final_state(initial_state, target_pol, kin);
			// Fill in branches.
			x = phase_space.x;
			y = phase_space.y;
			z = phase_space.z;
			ph_t_sq = phase_space.ph_t_sq;
			phi_h = phase_space.phi_h;
			phi = phase_space.phi;
			tau = phase_space.tau;
			phi_k = phase_space.phi_k;
			R = phase_space.R;
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
	events.Write();

	file_out.Close();

	return 0;
}

