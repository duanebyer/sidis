#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <sidis/sidis.hpp>

using namespace sidis;
using namespace sidis::cut;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::part;

// Plots the kinematically allowed phase space region. This program can generate
// plots in 1, 2, or 3 kinematic variables.
int main(int argc, char** argv) {
	int argc_root = 1;

	Real Mth = MASS_P + MASS_PI_0;
	part::Lepton beam = part::Lepton::E;
	part::Nucleus target = part::Nucleus::P;
	part::Hadron hadron = part::Hadron::PI_P;

	// Parse command line arguments.
	Real beam_energy;
	bool radiative;
	std::vector<std::string> axes;
	try {
		if (argc < 4 || argc > 6) {
			throw std::invalid_argument(
				"Unexpected number of command line arguments");
		}
		beam_energy = std::stold(argv[1]);
		std::string radiative_str = argv[2];
		if (radiative_str == "true"
				|| radiative_str == "on"
				|| radiative_str == "rad") {
			radiative = true;
		} else if (radiative_str == "false"
				|| radiative_str == "off"
				|| radiative_str == "nrad") {
			radiative = false;
		} else {
			throw std::out_of_range(
				"Must select radiative (rad) or "
				"non-radiative (nrad) phase space point");
		}

		for (int idx = 3; idx < argc; ++idx) {
			std::string axis_var = argv[idx];
			if (axis_var == "x"
					|| axis_var == "y" || axis_var == "Q_sq"
					|| axis_var == "z"
					|| axis_var == "ph_t" || axis_var == "ph_t_sq"
					|| axis_var == "phi_h"
					|| axis_var == "phi"
					|| (radiative && axis_var == "tau")
					|| (radiative && axis_var == "phi_k")
					|| (radiative && axis_var == "R")) {
				axes.push_back(axis_var);
			} else {
				throw std::out_of_range(
					"Unrecognized kinematic variable '" + axis_var + "'");
			}
		}
	} catch (std::exception const& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << "Usage: "
			<< "plot_phase_space <E_b> <rad,nrad> "
			<< "<kin. var. x> [kin. var. y] [kin. var. z]"
			<< std::endl;
		return 1;
	}

	std::vector<std::vector<Real> > data;
	std::vector<Real> data_mins;
	std::vector<Real> data_maxs;
	for (std::size_t idx = 0; idx < axes.size(); ++idx) {
		data.push_back(std::vector<Real>());
		data_maxs.push_back(std::numeric_limits<Real>::lowest());
		data_mins.push_back(std::numeric_limits<Real>::max());
	}

	// Generate events randomly in phase space.
	Particles ps(target, beam, hadron, Mth);
	Real S = 2. * ps.M * beam_energy;
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<Real> dist(0., 1.);
	std::cout << "Generating points in phase space." << std::endl;
	std::size_t num_points = 0;
	if (axes.size() == 1) {
		num_points = 100000;
	} else if (axes.size() == 2) {
		num_points = 1000000;
	} else if (axes.size() == 3) {
		num_points = 10000000;
	}
	for (std::size_t n = 0; n < num_points; ++n) {
		Real point[9] = {
			dist(rng), dist(rng), dist(rng),
			dist(rng), dist(rng), dist(rng),
			dist(rng), dist(rng), dist(rng),
		};
		Kinematics kin;
		KinematicsRad kin_rad;
		while (!cut::take(ps, S, point, &kin, nullptr)) { }
		while (!cut::take(kin, point + 6, &kin_rad, nullptr)) { }

		for (std::size_t idx = 0; idx < axes.size(); ++idx) {
			std::string axis_var = axes[idx];
			Real value = 0.;
			if (axis_var == "x") {
				value = kin.x;
			} else if (axis_var == "y") {
				value = kin.y;
			} else if (axis_var == "Q_sq") {
				value = kin.Q_sq;
			} else if (axis_var == "z") {
				value = kin.z;
			} else if (axis_var == "ph_t") {
				value = kin.ph_t;
			} else if (axis_var == "ph_t_sq") {
				value = kin.ph_t_sq;
			} else if (axis_var == "phi_h") {
				value = kin.phi_h;
			} else if (axis_var == "phi") {
				value = kin.phi;
			} else if (axis_var == "tau") {
				value = kin_rad.tau;
			} else if (axis_var == "phi_k") {
				value = kin_rad.phi_k;
			} else if (axis_var == "R") {
				value = kin_rad.R;
			}
			if (value < data_mins[idx]) {
				data_mins[idx] = value;
			}
			if (value > data_maxs[idx]) {
				data_maxs[idx] = value;
			}
			data[idx].push_back(value);
		}
	}

	// Plot the events.
	std::cout << "Plotting." << std::endl;
	TApplication app("Phase space plots", &argc_root, argv);
	TCanvas canvas("Phase space canvas");
	std::unique_ptr<TH1> hist;
	if (axes.size() == 1) {
		TH1D* h = new TH1D(
			"Phase space", "Phase space",
			100, data_mins[0], data_maxs[0]);
		h->FillN(data[0].size(), data[0].data(), nullptr);
		h->GetXaxis()->SetTitle(axes[0].c_str());
		h->Draw("hist");
		hist.reset(h);
	} else if (axes.size() == 2) {
		TH2D* h =new TH2D(
			"Phase space", "Phase space",
			100, data_mins[0], data_maxs[0],
			100, data_mins[1], data_maxs[1]);
		h->FillN(data[0].size(), data[0].data(), data[1].data(), nullptr);
		h->GetXaxis()->SetTitle(axes[0].c_str());
		h->GetYaxis()->SetTitle(axes[1].c_str());
		h->Draw("colz");
		hist.reset(h);
	} else if (axes.size() == 3) {
		TH3D* h = new TH3D(
			"Phase space", "Phase space",
			100, data_mins[0], data_maxs[0],
			100, data_mins[1], data_maxs[1],
			100, data_mins[2], data_maxs[2]);
		for (std::size_t idx = 0; idx < data[0].size(); ++idx) {
			h->Fill(data[0][idx], data[1][idx], data[2][idx]);
		}
		h->GetXaxis()->SetTitle(axes[0].c_str());
		h->GetYaxis()->SetTitle(axes[1].c_str());
		h->GetZaxis()->SetTitle(axes[2].c_str());
		h->Draw("box");
		hist.reset(h);
	}

	app.Run();

	return 0;
}

