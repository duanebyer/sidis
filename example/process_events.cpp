#include <ios>
#include <iostream>
#include <string>

#include <TApplication.h>
#include <TArrayD.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TParameter.h>
#include <TTree.h>

Double_t const MASS_P = 0.9382720813;

// Reads events from a ROOT file produced by `sidisgen`, integrates the
// cross-section, and makes an x-Q2 phase space plot.
int main(int argc, char** argv) {
	int argc_root = 1;
	if (argc != 2) {
		std::cerr << "Usage: process_events <ROOT file>" << std::endl;
		return 1;
	}

	std::string file_name = argv[1];
	TApplication app("Processing events", &argc_root, argv);

	// Load events from file.
	TFile file(file_name.c_str());
	TTree* events = file.Get<TTree>("events");
	// Load the normalizations for events.
	TArrayD* norm_arr = file.Get<TArrayD>("stats/norm");
	// Load parameters (beam energy and target type).
	TParameter<Double_t>* beam_energy_param = file.Get<TParameter<Double_t> >("params/setup.beam_energy");
	TParameter<Int_t>* target_param = file.Get<TParameter<Int_t> >("params/setup.target");
	if (events == nullptr
			|| norm_arr == nullptr
			|| beam_energy_param == nullptr
			|| target_param == nullptr) {
		std::cerr << "Error: failed to load some objects from file '"
			<< file_name << "'." << std::endl;
		return 1;
	}
	Double_t beam_energy = beam_energy_param->GetVal();
	Int_t target = target_param->GetVal();
	Double_t target_mass = 0.;

	if (target == 0) {
		// Proton target.
		target_mass = MASS_P;
	} else {
		// Only handle proton targets for now, if we get something else, quit.
		std::cerr << "Error: only proton targets are supported." << std::endl;
		return 1;
	}
	// Related to COM energy.
	Double_t S = 2. * beam_energy * target_mass;

	// Create histogram.
	TH2D hist = TH2D(
		"hist", "Phase space coverage",
		100, 0., 1.,
		100, 0., S);
	hist.GetXaxis()->SetTitle("x");
	hist.GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
	hist.SetStats(0);

	// Set up branches.
	Double_t x;
	Double_t y;
	Double_t weight;
	Int_t type;
	events->SetBranchAddress("x", &x);
	events->SetBranchAddress("y", &y);
	events->SetBranchAddress("weight", &weight);
	events->SetBranchAddress("type", &type);

	// Keep track of some running totals.
	Double_t weight_total = 0.;
	Double_t weight_sq_total = 0.;

	Long64_t num_events = events->GetEntries();
	for (Long64_t idx = 0; idx < num_events; ++idx) {
		events->GetEntry(idx);
		// The norm must be looked up separately for each type of event.
		Double_t norm = norm_arr->At(type);

		// Fill the histogram. Note that the correct event weight is given by
		// the norm and the weight together. The norm provides an overall
		// scaling for all events of a certaint type (ex. non-radiative or
		// radiative), whereas the weight changes event-by-event.
		Double_t Q_sq = S * x * y;
		hist.Fill(x, Q_sq, norm * weight);

		weight_total += weight * norm;
		weight_sq_total += TMath::Sq(weight * norm);
	}

	// Calculate total cross-section with error.
	Double_t xs = weight_total;
	Double_t xs_err = TMath::Sqrt(weight_sq_total);
	std::cout << std::scientific << std::setprecision(4);
	std::cout << "Cross-section: " << xs
		<< " ± " << xs_err
		<< " GeV⁻²" << std::endl;

	// Draw the histogram.
	TCanvas canvas("canvas");
	hist.Draw("col");

	app.Run();

	return 0;
}

