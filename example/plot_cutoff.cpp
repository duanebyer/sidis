#include <iostream>
#include <memory>
#include <vector>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLine.h>

#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>

using namespace sidis;

int main(int argc, char** argv) {
	TApplication app("Soft photon cutoff plots", &argc, argv);

	Real Mth = MASS_P + MASS_PI_0;
	part::Lepton beam = part::Lepton::E;
	part::Nucleus target = part::Nucleus::P;
	part::Hadron hadron = part::Hadron::PI_P;

	// Read input parameters from command line.
	Real beam_energy;
	std::unique_ptr<sidis::sf::SfSet> sf(new sf::set::ProkudinSfSet());
	Real x, y, z, ph_t_sq, phi_h, phi;
	math::IntegParams params { math::IntegMethod::CUBATURE, 10000, 1000 };
	//math::IntegParams params { math::IntegMethod::CUBATURE, 100000, 10000 };

	if (argc != 8) {
		std::cout << "Wrong number of arguments." << std::endl;
		std::cout << "Usage" << std::endl
			<< "plot_cutoff "
			<< "<E_b> "
			<< "<x> <y> <z> <ph_t²> <φ_h> <φ>" << std::endl;
		return 1;
	}

	beam_energy = std::stold(argv[1]);
	x = std::stold(argv[2]);
	y = std::stold(argv[3]);
	z = std::stold(argv[4]);
	ph_t_sq = std::stold(argv[5]);
	phi_h = std::stold(argv[6]);
	phi = std::stold(argv[7]);

	part::Particles ps(target, beam, hadron, Mth);
	Real S = 2. * ps.M * beam_energy;
	kin::PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, phi };
	kin::Kinematics kin(ps, S, ph_space);

	math::EstErr xs = xs::nrad_integ(kin, *sf, 0., math::VEC3_ZERO, INF, params);
	std::cout << "Total cross-section: "
		<< xs.val << " ± " << xs.err << std::endl;

	std::vector<Real> k0_cut_arr;
	std::vector<Real> xs_rad_f_arr;
	std::vector<Real> xs_rad_f_err_arr;

	Real k0_min_pow = -3.;
	Real k0_max_pow = 0.;
	std::size_t k0_steps = 100;
	for (std::size_t idx = 0; idx < k0_steps; ++idx) {
		std::cout << "Progress: " << (idx + 1) << " / " << k0_steps << std::endl;
		Real lambda = static_cast<Real>(idx) / (k0_steps - 1);
		Real k0_pow = lambda * k0_min_pow + (1. - lambda) * k0_max_pow;
		Real k0_cut = std::pow(10., k0_pow);
		math::EstErr xs_rad_f = xs::rad_f_integ(kin, *sf, 0., math::VEC3_ZERO, k0_cut, params);

		Real ratio = xs_rad_f.val / xs.val;
		Real err_rel = std::abs(xs_rad_f.err / xs_rad_f.val) + std::abs(xs.err / xs.val);
		Real err = err_rel * ratio;
		k0_cut_arr.push_back(k0_cut);
		xs_rad_f_arr.push_back(-100. * ratio);
		xs_rad_f_err_arr.push_back(100. * err);
	}
	std::cout << "Finished!" << std::endl;

	TCanvas canvas("Effect of soft photon cutoff");
	canvas.SetLogx();

	TGraphErrors graph(
		k0_steps,
		k0_cut_arr.data(), xs_rad_f_arr.data(),
		nullptr, xs_rad_f_err_arr.data());
	graph.SetTitle("");
	graph.GetXaxis()->SetTitle("Soft photon cutoff  #bar{k}_{0} (GeV)");
	graph.GetYaxis()->SetTitle("Relative IR-free part  #minus#sigma_{R}^{F} / #sigma_{SIDIS}^{nrad} (%)");
	graph.GetXaxis()->SetLimits(std::pow(10., k0_min_pow), std::pow(10., k0_max_pow));
	graph.SetFillColor(6);
	graph.SetLineColor(6);
	graph.SetLineWidth(2);
	graph.SetFillStyle(1001);
	graph.Draw("AL3");

	app.Run();

	return 0;
}

