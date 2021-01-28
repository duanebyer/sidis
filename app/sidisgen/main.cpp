#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <TBranch.h>
#include <TError.h>
#include <TFile.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TTree.h>

#include <sidis/sidis.hpp>
#include <sidis/extra/bounds.hpp>
#include <sidis/extra/math.hpp>
#include <sidis/extra/vector.hpp>
#include <sidis/sf_set/ww.hpp>

#include "params.hpp"
#include "utility.hpp"

using namespace sidis;

namespace {

int const SUCCESS = 0;
int const ERROR_ARG_PARSE = -1;
int const ERROR_FILE_NOT_FOUND = -2;
int const ERROR_FILE_NOT_CREATED = -3;
int const ERROR_PARAMS_PARSE = -4;
int const ERROR_FOAM_INCOMPATIBLE = -5;
int const ERROR_FOAM_NOT_FOUND = -6;

sf::model::WW const ww;
Real const PI = constant::PI;

// Converts between the `sidis` 4-vector type and the `ROOT` 4-vector type.
TLorentzVector convert_vec4(math::Vec4 v) {
	return TLorentzVector(v.x, v.y, v.z, v.t);
}

struct XsNRad : public TFoamIntegrand {
	Params const params;
	kin::Particles const ps;
	Real const S;

	explicit XsNRad(Params params) :
		params(params),
		ps(params.target, params.beam, params.hadron, params.mass_threshold),
		S(2. * mass(params.target) * params.beam_energy) { }

	// Converts from the unit hyper-cube to a point in phase space, also
	// providing the associated jacobian.
	kin::Kinematics GetKinematics(Double_t const* vec, Real* jacobian) {
		math::Bounds x_bounds = kin::x_bounds(ps, S);
		Real x = x_bounds.lerp(vec[0]);
		math::Bounds y_bounds = kin::y_bounds(ps, S, x);
		Real y = y_bounds.lerp(vec[1]);
		math::Bounds z_bounds = kin::z_bounds(ps, S, x, y);
		Real z = z_bounds.lerp(vec[2]);
		math::Bounds ph_t_sq_bounds = kin::ph_t_sq_bounds(ps, S, x, y, z);
		Real ph_t_sq = ph_t_sq_bounds.lerp(vec[3]);
		math::Bounds phi_h_bounds = math::Bounds(0., 2. * PI);
		Real phi_h = phi_h_bounds.lerp(vec[4]);
		math::Bounds phi_bounds = math::Bounds(0., 2. * PI);
		Real phi = phi_bounds.lerp(vec[5]);
		*jacobian = x_bounds.size() * y_bounds.size() * z_bounds.size()
			* ph_t_sq_bounds.size() * phi_h_bounds.size() * phi_bounds.size();

		kin::PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, phi };
		kin::Kinematics kin(ps, S, phase_space);
		return kin;
	}

	Double_t Density(int dim, Double_t* vec) override {
		if (dim != 6) {
			return 0.;
		}
		Real jacobian;
		kin::Kinematics kin = GetKinematics(vec, &jacobian);
		// TODO: Remove this condition, once we feel confident enough in the
		// kinematic bounds functions.
		if (!kin::valid(kin)) {
			return 0.;
		}
		math::Vec3 eta = frame::hadron_from_target(kin) * params.target_pol;
		// TODO: Evaluate when it is a good approximation to say that
		// `nrad ~ nrad_ir`. This happens because for small `k0_cut`, the
		// contribution of `rad_f` integrated up to `k0_cut` becomes vanishingly
		// small, so it can be neglected. However, this must be balanced with
		// choosing `k0_cut` to be non-zero to avoid the infrared divergence in
		// the radiative part of the cross-section.
		Real xs = xs::nrad_ir(params.beam_pol, eta, kin, ww, params.k0_cut);
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
	Params const params;
	kin::Particles const ps;
	Real const S;

	explicit XsRad(Params params) :
		params(params),
		ps(params.target, params.beam, params.hadron, params.mass_threshold),
		S(2. * mass(params.target) * params.beam_energy) { }

	kin::KinematicsRad GetKinematics(Double_t const* vec, Real* jacobian) {
		math::Bounds x_bounds = kin::x_bounds(ps, S);
		Real x = x_bounds.lerp(vec[0]);
		math::Bounds y_bounds = kin::y_bounds(ps, S, x);
		Real y = y_bounds.lerp(vec[1]);
		math::Bounds z_bounds = kin::z_bounds(ps, S, x, y);
		Real z = z_bounds.lerp(vec[2]);
		math::Bounds ph_t_sq_bounds = kin::ph_t_sq_bounds(ps, S, x, y, z);
		Real ph_t_sq = ph_t_sq_bounds.lerp(vec[3]);
		math::Bounds phi_h_bounds = math::Bounds(0., 2. * PI);
		Real phi_h = phi_h_bounds.lerp(vec[4]);
		math::Bounds phi_bounds = math::Bounds(0., 2. * PI);
		Real phi = phi_bounds.lerp(vec[5]);

		kin::PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, phi };
		kin::Kinematics kin(ps, S, phase_space);

		math::Bounds tau_bounds = kin::tau_bounds(kin);
		Real tau = tau_bounds.lerp(vec[6]);
		math::Bounds phi_k_bounds = math::Bounds(0., 2. * PI);
		Real phi_k = phi_k_bounds.lerp(vec[7]);
		math::Bounds R_bounds = kin::R_bounds_hard(kin, tau, phi_k, params.k0_cut);
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
		kin::Kinematics kin = kin_rad.project();
		if (!kin::valid(kin_rad)) {
			return 0.;
		}
		math::Vec3 eta = frame::hadron_from_target(kin) * params.target_pol;
		Real xs = xs::rad(params.beam_pol, eta, kin_rad, ww);
		if (std::isnan(xs)) {
			return 0.;
		} else {
			return jacobian * xs;
		}
	}
};

int command_help() {
	std::cout
		<< "Usage:"                                        << std::endl
		<< "  Prepare FOAM for Monte-Carlo generation"     << std::endl
		<< "    sidisgen --initialize <parameter file>"    << std::endl
		<< "  Generate events"                             << std::endl
		<< "    sidisgen --generate <parameter file>"      << std::endl
		<< "  List parameters used to produce file"        << std::endl
		<< "    sidisgen --inspect <output file>"          << std::endl
		<< std::endl
		<< "Parameter file format summary (see docs):"     << std::endl
		<< std::endl
		<< "foam_nrad_file <ROOT file>"                      << std::endl
		<< "foam_rad_file  <ROOT file>"                      << std::endl
		<< "event_file     <ROOT file>"                      << std::endl
		<< "num_events     <integer>"                        << std::endl
		<< "num_init       <integer>"                        << std::endl
		<< "seed           <integer>"                        << std::endl
		<< "seed_init      <integer>"                        << std::endl
		<< "beam_energy    <energy>"                         << std::endl
		<< "beam           <pid>"                            << std::endl
		<< "target         <pid>"                            << std::endl
		<< "hadron         <pid>"                            << std::endl
		<< "beam_pol       <U, L>"                           << std::endl
		<< "target_pol     <U, L, T>"                        << std::endl
		// TODO: We could add an intermediate method here that approximates
		// the full 3-d integral over `rad_f` with a 2-d integral over
		// `phi_k` and `tau` only.
		<< "UNIMPLEMENTED OPTIONS"                         << std::endl
		<< "rc_method      <none, full>"                     << std::endl
		<< "nrad_method    <none, ir_only, full>"            << std::endl
		<< "soft_k_cut     <energy>"                         << std::endl
		<< "sf_method      <ww>"                             << std::endl
		<< "x_cut          <min> <max>"                      << std::endl
		<< "Q2_cut         <min> <max>"                      << std::endl
		<< "y_cut          <min> <max>"                      << std::endl
		<< "z_cut          <min> <max>"                      << std::endl
		<< "ph_t_cut       <min> <max>"                      << std::endl;
	return SUCCESS;
}

int command_version() {
	// TODO: Output correct version and build information.
	std::cout << "sidisgen "
		<< SIDIS_VERSION_MAJOR << "."
		<< SIDIS_VERSION_MINOR << "."
		<< SIDIS_VERSION_PATCH << "."
		<< SIDIS_VERSION_TWEAK << std::endl;
	return SUCCESS;
}

int command_inspect(std::string output_file_name) {
	Params params;
	TFile file(output_file_name.c_str());
	if (file.IsZombie()) {
		std::cerr
			<< "Output file '" << output_file_name
			<< "' not found." << std::endl;
		return ERROR_FILE_NOT_FOUND;
	}
	params.read_root(file);
	params.write(std::cout);
	return SUCCESS;
}

int command_initialize(std::string params_file_name) {
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		std::cerr
			<< "Parameter file '" << params_file_name
			<< "' not found." << std::endl;
		return ERROR_FILE_NOT_FOUND;
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params;
	try {
		params.read(params_file);
	} catch (std::runtime_error const& e) {
		std::cerr
			<< "Failed to parse parameter file '" << params_file_name << "': "
			<< e.what() << "." << std::endl;
		return ERROR_PARAMS_PARSE;
	}
	std::cout << std::endl;
	params.write(std::cout);
	std::cout << std::endl;

	ULong_t N_init_nrad = params.num_init >= 1 ? params.num_init : 1;
	ULong_t N_init_rad  = params.num_init >= 1 ? params.num_init : 1;
	// Note that `TRandom3` uses the time as the seed if zero is provided.
	UInt_t seed = params.seed_init >= 0 ? params.seed_init : 0;
	TRandom3 random(seed);

	std::cout << "Creating non-radiative FOAM file." << std::endl;
	TFile foam_nrad_file(params.foam_nrad_file.c_str(), "RECREATE");
	if (foam_nrad_file.IsZombie()) {
		std::cerr
			<< "Couldn't open or create file '" << params.foam_nrad_file
			<< "'." << std::endl;
		return ERROR_FILE_NOT_CREATED;
	}
	foam_nrad_file.cd();
	params.write_root(foam_nrad_file);

	std::cout << "Non-radiative FOAM initialization." << std::endl;
	TFoam foam_nrad("FoamNRad");
	XsNRad xs_nrad(params);
	foam_nrad.SetChat(0);
	foam_nrad.SetkDim(6);
	foam_nrad.SetRho(&xs_nrad);
	foam_nrad.SetPseRan(&random);
	foam_nrad.SetnSampl(N_init_nrad);
	foam_nrad.Initialize();
	foam_nrad.Write("FoamNRad");
	foam_nrad_file.Close();

	std::cout << "Creating radiative FOAM file." << std::endl;
	TFile foam_rad_file(params.foam_rad_file.c_str(), "RECREATE");
	if (foam_rad_file.IsZombie()) {
		std::cerr
			<< "Couldn't open or create file '" << params.foam_rad_file
			<< "'." << std::endl;
		return ERROR_FILE_NOT_CREATED;
	}
	foam_rad_file.cd();
	params.write_root(foam_rad_file);

	std::cout << "Radiative FOAM initialization." << std::endl;
	TFoam foam_rad("FoamRad");
	XsRad xs_rad(params);
	foam_rad.SetChat(0);
	foam_rad.SetkDim(9);
	foam_rad.SetRho(&xs_rad);
	foam_rad.SetPseRan(&random);
	foam_rad.SetnSampl(N_init_rad);
	foam_rad.Initialize();
	foam_rad.Write("FoamRad");
	foam_rad_file.Close();

	std::cout << "Finished!" << std::endl;
	return SUCCESS;
}

int command_generate(std::string params_file_name) {
	std::ifstream params_file(params_file_name);
	if (!params_file) {
		std::cerr
			<< "Parameter file '" << params_file_name
			<< "' not found." << std::endl;
		return ERROR_FILE_NOT_FOUND;
	}
	std::cout << "Reading parameter file '" << params_file_name << "'." << std::endl;
	Params params;
	try {
		params.read(params_file);
	} catch (std::runtime_error const& e) {
		std::cerr
			<< "Failed to parse parameter file '" << params_file_name << "': "
			<< e.what() << "." << std::endl;
		return ERROR_PARAMS_PARSE;
	}
	std::cout << std::endl;
	params.write(std::cout);
	std::cout << std::endl;

	TFile event_file(params.event_file.c_str(), "RECREATE");
	if (event_file.IsZombie()) {
		std::cerr
			<< "Couldn't create file '" << params.event_file
			<< "'." << std::endl;
		return ERROR_FILE_NOT_CREATED;
	}
	TFile foam_nrad_file(params.foam_nrad_file.c_str(), "OPEN");
	if (foam_nrad_file.IsZombie()) {
		std::cerr
			<< "Couldn't open file '" << params.foam_nrad_file
			<< "'." << std::endl;
		return ERROR_FILE_NOT_CREATED;
	}
	TFile foam_rad_file(params.foam_rad_file.c_str(), "OPEN");
	if (foam_rad_file.IsZombie()) {
		std::cerr
			<< "Couldn't open file '" << params.foam_rad_file
			<< "'." << std::endl;
		return ERROR_FILE_NOT_CREATED;
	}

	std::cout << "Verifying FOAMs use the same parameters." << std::endl;
	Params foam_nrad_params;
	foam_nrad_params.read_root(foam_nrad_file);
	if (!params.compatible_foam(foam_nrad_params)) {
		std::cerr
			<< "Couldn't use non-radiative FOAM from '" << params.foam_nrad_file
			<< "' because it was generated with parameters incompatible with "
			<< "the provided parameter file '" << params_file_name
			<< "'." << std::endl;
		return ERROR_FOAM_INCOMPATIBLE;
	}
	Params foam_rad_params;
	foam_rad_params.read_root(foam_rad_file);
	if (!params.compatible_foam(foam_rad_params)) {
		std::cerr
			<< "Couldn't use radiative FOAM from '" << params.foam_rad_file
			<< "' because it was generated with different parameters than "
			<< "provided parameter file '" << params_file_name
			<< "'." << std::endl;
		return ERROR_FOAM_INCOMPATIBLE;
	}

	std::cout << "Reading non-radiative FOAM from file." << std::endl;
	TFoam* foam_nrad = foam_nrad_file.Get<TFoam>("FoamNRad");
	if (foam_nrad == nullptr) {
		std::cerr
			<< "Failed to load non-radiative FOAM from file '"
			<< params.foam_nrad_file << "'.";
		return ERROR_FOAM_NOT_FOUND;
	}

	std::cout << "Reading radiative FOAM from file." << std::endl;
	TFoam* foam_rad = foam_rad_file.Get<TFoam>("FoamRad");
	if (foam_rad == nullptr) {
		std::cerr
			<< "Failed to load radiative FOAM from file '"
			<< params.foam_rad_file << "'.";
		return ERROR_FOAM_NOT_FOUND;
	}

	UInt_t seed = params.seed >= 0 ? params.seed : 0;
	TRandom3 random(seed);

	XsNRad xs_nrad(params);
	foam_nrad->SetPseRan(&random);
	foam_nrad->ResetRho(&xs_nrad);

	XsRad xs_rad(params);
	foam_rad->SetPseRan(&random);
	foam_rad->ResetRho(&xs_rad);

	constant::Nucleus target = params.target;
	constant::Lepton beam = params.beam;
	constant::Hadron hadron = params.hadron;
	math::Vec3 target_pol = params.target_pol;
	kin::Particles ps(target, beam, hadron, params.mass_threshold);

	kin::Initial initial_state(ps, params.beam_energy);
	ULong_t N_gen = params.num_events >= 0 ? params.num_events : 0;
	// These keep track of how many radiative/non-radiative events have been
	// generated so far, so that the ratio of radiative to non-radiative events
	// generated can be kept roughly the same as the ratio of the total
	// radiative/non-radiative cross-sections.
	ULong_t N_gen_nrad = 0;
	ULong_t N_gen_rad = 0;

	event_file.cd();
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
	params.write_root(event_file);

	std::cout << "Generating events." << std::endl;
	Double_t total_nrad, total_nrad_err;
	Double_t total_rad, total_rad_err;

	bool update_progress = true;
	ULong_t percent = 0;
	ULong_t next_percent_rem = N_gen % 100;
	ULong_t next_percent = N_gen / 100;
	for (ULong_t event_idx = 0; event_idx < N_gen; ++event_idx) {
		while (event_idx >= next_percent + (next_percent_rem != 0)) {
			percent += 1;
			next_percent = math::prod_div(N_gen, percent, 100, next_percent_rem);
			update_progress = true;
		}
		if (update_progress) {
			write_progress_bar(std::cout, percent);
			std::cout << '\r';
			std::cout.flush();
			update_progress = false;
		}
		// Estimate the total radiative and non-radiative cross-sections and
		// generate a radiative/non-radiative event accordingly. The total
		// cross-section estimates are improved as more events are generated.
		foam_nrad->GetIntegMC(total_nrad, total_nrad_err);
		foam_rad->GetIntegMC(total_rad, total_rad_err);
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
			weight = foam_nrad->MCgenerate(event_vec);
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
			weight = foam_rad->MCgenerate(event_vec);
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
	write_progress_bar(std::cout, 100);
	std::cout << std::endl;
	std::cout << "Writing events to file." << std::endl;
	event_file.cd();
	events.Write();

	foam_nrad->GetIntegMC(total_nrad, total_nrad_err);
	foam_rad->GetIntegMC(total_rad, total_rad_err);

	// Write total cross-sections to file.
	event_file.cd();
	TParameter<Double_t> p_xs_total_nrad("xs_total_nrad", total_nrad);
	TParameter<Double_t> p_xs_total_nrad_err("xs_total_nrad_err", total_nrad_err);
	TParameter<Double_t> p_xs_total_rad("xs_total_rad", total_rad);
	TParameter<Double_t> p_xs_total_rad_err("xs_total_rad_err", total_rad_err);
	p_xs_total_nrad.Write();
	p_xs_total_nrad_err.Write();
	p_xs_total_rad.Write();
	p_xs_total_rad_err.Write();

	event_file.Close();
	foam_nrad_file.Close();
	foam_rad_file.Close();

	std::cout << "Statistics:" << std::endl;
	std::cout << "\tNumber of non-radiative events:    " << N_gen_nrad << std::endl;
	std::cout << "\tNumber of radiative events:        " << N_gen_rad << std::endl;
	std::cout << "\tTotal non-radiative cross-section: "
		<< std::scientific << total_nrad << " ± " << total_nrad_err << std::endl;
	std::cout << "\tTotal radiative cross-section:     "
		<< std::scientific << total_rad << " ± " << total_rad_err << std::endl;
	return SUCCESS;
}

}

int main(int argc, char** argv) {
#ifdef NDEBUG
	gErrorIgnoreLevel = kFatal;
	gErrorAbortLevel = kFatal;
#else
	gErrorIgnoreLevel = kWarning;
	gErrorAbortLevel = kError;
#endif
	// Parse arguments.
	if (argc <= 1) {
		std::cout << "Try `sidisgen --help`." << std::endl;
		return SUCCESS;
	}
	std::string command = argv[1];
	if (command == "--help" || command == "-?") {
		return command_help();
	} else if (command == "--version" || command == "-v") {
		return command_version();
	} else if (command == "--inspect" || command == "-i") {
		if (argc > 3) {
			std::cerr
				<< "Unexpected argument " << argv[3] << "." << std::endl;
			return ERROR_ARG_PARSE;
		} else if (argc < 3) {
			std::cerr
				<< "Expected ROOT file argument for inspection." << std::endl;
			return ERROR_ARG_PARSE;
		}
		return command_inspect(argv[2]);
	} else if (command == "--initialize" || command == "-i") {
		if (argc > 3) {
			std::cerr
				<< "Unexpected argument " << argv[3] << "." << std::endl;
			return ERROR_ARG_PARSE;
		} else if (argc < 3) {
			std::cerr
				<< "Expected parameter file argument." << std::endl;
			return ERROR_ARG_PARSE;
		}
		return command_initialize(argv[2]);
	} else if (command == "--generate" || command == "-g") {
		if (argc > 3) {
			std::cerr
				<< "Unexpected argument " << argv[3] << "." << std::endl;
			return ERROR_ARG_PARSE;
		} else if (argc < 3) {
			std::cerr
				<< "Expected parameter file argument." << std::endl;
			return ERROR_ARG_PARSE;
		}
		return command_generate(argv[2]);
	} else {
		std::cerr << "Unrecognized command " << command << "." << std::endl;
		return ERROR_ARG_PARSE;
	}
}

