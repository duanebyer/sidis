#include <array>

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>

#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>

using namespace sidis;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::part;

// Computes and plots the Sivers and Collins asymmetries, with and without
// radiative corrections.
int main(int argc, char** argv) {
	TApplication app("Transverse asymmetry plots", &argc, argv);

	sf::set::ProkudinSfSet const sf_set;

	Real E_b = 10.6;
	Particles ps(
		Nucleus::P,
		Lepton::E,
		Hadron::PI_P,
		MASS_P + MASS_PI_0);
	Real S = 2. * ps.M * E_b;
	Real x = 0.2;
	Real y = 0.754;

	// The values of `ph_t` to evaluate the transverse asymmetries.
	std::array<Real, 3> ph_t_vals { 0.7, 0.4, 0.1 };
	// The bounds on `z` must be slightly more restrictive than `[0, 1]` for
	// kinematic reasons.
	Real z_min = 0.1;
	Real z_max = 0.8;

	IntegParams params_sivers { IntegMethod::VEGAS, 10000, 1000 };
	IntegParams params_collins { IntegMethod::VEGAS, 10000, 1000 };
	TF1 f_sivers("Sivers Born", [&](Double_t* xs, Double_t* arr) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(arr[0] * arr[0]);
		return asym::sivers(sf_set, ps, S, x, y, z, ph_t_sq, false, params_sivers).val;
	}, z_min, z_max, 1);
	TF1 f_sivers_rc("Sivers RC", [&](Double_t* xs, Double_t* arr) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(arr[0] * arr[0]);
		return asym::sivers(sf_set, ps, S, x, y, z, ph_t_sq, true, params_sivers).val;
	}, z_min, z_max, 1);
	TF1 f_collins("Collins Born", [&](Double_t* xs, Double_t* arr) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(arr[0] * arr[0]);
		return asym::collins(sf_set, ps, S, x, y, z, ph_t_sq, false, params_collins).val;
	}, z_min, z_max, 1);
	TF1 f_collins_rc("Collins RC", [&](Double_t* xs, Double_t* arr) {
		Real z = static_cast<Real>(xs[0]);
		Real ph_t_sq = static_cast<Real>(arr[0] * arr[0]);
		return asym::collins(sf_set, ps, S, x, y, z, ph_t_sq, true, params_collins).val;
	}, z_min, z_max, 1);

	Color_t color = 1;
	Style_t style = 0;
	Style_t style_rc = 1;
	f_sivers.SetLineColor(color);
	f_sivers.SetLineStyle(style);
	f_sivers_rc.SetLineColor(color);
	f_sivers_rc.SetLineStyle(style_rc);
	f_collins.SetLineColor(color);
	f_collins.SetLineStyle(style);
	f_collins_rc.SetLineColor(color);
	f_collins_rc.SetLineStyle(style_rc);

	TCanvas canvas_sivers("Sivers asymmetry");
	canvas_sivers.DrawFrame(0., 0., 1., 0.15, "Sivers asymmetry");
	for (Real ph_t : ph_t_vals) {
		f_sivers.SetParameter(0, ph_t);
		f_sivers_rc.SetParameter(0, ph_t);
		f_sivers.DrawCopy("LSAME");
		f_sivers_rc.DrawCopy("LSAME");
	}

	TCanvas canvas_collins("Collins asymmetry");
	canvas_collins.DrawFrame(0., 0., 1., 0.04, "Collins asymmetry");
	for (Real ph_t : ph_t_vals) {
		f_collins.SetParameter(0, ph_t);
		f_collins_rc.SetParameter(0, ph_t);
		f_collins.DrawCopy("LSAME");
		f_collins_rc.DrawCopy("LSAME");
	}

	TFile output("out.root", "RECREATE");
	f_sivers.Write();
	f_sivers_rc.Write();
	f_collins.Write();
	f_collins_rc.Write();
	canvas_sivers.Write();
	canvas_collins.Write();
	output.Close();

	app.Run();

	return 0;
}

