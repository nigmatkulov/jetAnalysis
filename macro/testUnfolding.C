R__LOAD_LIBRARY(/Users/gnigmat/work/RooUnfold/build/libRooUnfold.dylib)

#include "/Users/gnigmat/work/RooUnfold/build/RooUnfoldResponse.h"
#include "/Users/gnigmat/work/RooUnfold/build/RooUnfoldBayes.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TCanvas.h"

void testUnfolding() {
    // Toy Monte Carlo parameters
    const int nbins = 20;
    const double xmin = 0.0;
    const double xmax = 10.0;

    // True histogram (truth distribution)
    TH1D *hTrue = new TH1D("hTrue", "True distribution", nbins, xmin, xmax);

    // Measured histogram (simulated data)
    TH1D *hMeasured = new TH1D("hMeasured", "Measured distribution", nbins, xmin, xmax);

    // Response matrix (migration matrix)
    RooUnfoldResponse response(hMeasured, hTrue);

    // Generate toy MC data
    TRandom3 rnd(0);
    for (int i = 0; i < 1000000; ++i) {
        double trueValue = rnd.Uniform(xmin, xmax);
        double measuredValue = trueValue + rnd.Gaus(0.0, 0.2);  // Smearing
        hTrue->Fill(trueValue);
        hMeasured->Fill(measuredValue);
        response.Fill(measuredValue, trueValue);  // Fill response matrix
    }

    // Unfold using RooUnfoldBayes
    // RooUnfoldBayes unfold(&response, hMeasured, 4, kFALSE);  // Number of iterations = 4
    // TH1D *hUnfolded = (TH1D*)unfold.Hunfold();

    // Create inverted unfolding
    RooUnfoldInvert unfold(&response, hMeasured);

    // Retrieve inverted distribution
    TH1D* hUnfolded = (TH1D*)unfold.Hunfold();

    // Plotting
    TCanvas *canvas = new TCanvas("canvas", "RooUnfold Example", 800, 600);
    canvas->Divide(2, 2);

    canvas->cd(1);
    hTrue->SetLineColor(kBlue);
    hTrue->Draw();
    hTrue->SetTitle("True Distribution");

    canvas->cd(2);
    hMeasured->SetLineColor(kRed);
    hMeasured->Draw();
    hMeasured->SetTitle("Measured Distribution");

    canvas->cd(3);
    TH2D *hResponse = (TH2D*)response.Hresponse();
    hResponse->Draw("colz");
    hResponse->SetTitle("Response Matrix");

    canvas->cd(4);
    hUnfolded->SetLineColor(kGreen);
    hUnfolded->Draw();
    hUnfolded->SetTitle("Unfolded Distribution");

    //canvas->SaveAs("RooUnfoldExample.png");
}

// CMakeLists.txt
// -find_package( ROOT COMPONENTS Tree Unfold Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser REQUIRED)
// +find_package( ROOT COMPONENTS Tree Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser REQUIRED)

// cmake/PlainROOT.cmake
// -find_package( ROOT COMPONENTS Tree Unfold Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser)
// +find_package( ROOT COMPONENTS Tree Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser)
