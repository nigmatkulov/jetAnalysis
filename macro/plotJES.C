#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TLine.h"
#include "TSystem.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TLatex.h"
#include "THnSparse.h"
#include "TLegend.h"

// C++ headers
#include <iostream>
#include <vector>
#include <sys/stat.h>

//________________
bool directoryExists(const char* directoryPath) {
    TSystemFile file(directoryPath, "");
    return file.IsDirectory();
}

//________________
void createDirectory(const char* directoryPath) {
    // Create the directory with read, write, and execute permissions for owner, group, and others
    if (mkdir(directoryPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
        std::cout << "Directory created successfully." << std::endl;
    } else {
        std::cerr << "Failed to create directory." << std::endl;
    }
}

//________________
void plotCMSHeader(int collSystem = 1, double energy = 8.16) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Simulation}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, Form("%s #sqrt{s_{NN}} = %3.2f TeV", collSystemStr.Data(), energy) );
    t.SetTextSize(0.05);
}

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.2;
    Int_t lineWidth = 2;
    Int_t color = 2;
    if (type == 0) {
        color = 2;        // red
        markerStyle = 20; // filled circle
    }
    else if (type == 1) {
        color = 4;        // blue
        markerStyle = 24; // open circle
    }
    else if (type == 2) {
        color = 1;        // black
        markerStyle = 20; // filled triangle
    }
    else if (type == 3) {
        color = 2;        // red
        markerStyle = 24; // open circle
    }
    else if (type == 4) {
        color = 4;        // blue
        markerStyle = 20; // filled circle
    }
    else {
        color = 6;        // magenta
        markerStyle = 30; // open star
    }

    h->SetLineWidth( lineWidth );
    h->SetLineColor( color );
    
    h->SetMarkerStyle( markerStyle );
    h->SetMarkerColor( color );
    h->SetMarkerSize( markerSize );

    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void set2DStyle(TH2* h) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
    h->GetYaxis()->SetTitleOffset(1.0);
}

//________________
void retrieveJESandJER(TH2D *h2D, TH1D*& hJES, TH1D*& hJER, int type = 0, bool performFit = false) {
    // Initialize pointers to nullptr
    // hJES = nullptr;
    // hJER = nullptr;
    
    // Check if input histogram is valid
    if (!h2D) {
        std::cerr << "[ERROR] Input histogram is nullptr" << std::endl;
        return;
    }
    
    // Check if histogram has sufficient entries for fitting
    if (h2D->GetEntries() == 0) {
        std::cerr << "[WARNING] Histogram " << h2D->GetName() << " has insufficient entries (" 
                  << h2D->GetEntries() << ") for fitting" << std::endl;
        return;
    }
    
    // Retrieve JES and JER
    h2D->RebinX(2); // Rebin for better statistics
    h2D->FitSlicesY();

    // Jet energy scale
    hJES = dynamic_cast<TH1D*>(gDirectory->Get( Form("%s_1", h2D->GetName()) ));
    if (!hJES) {
        std::cerr << "[ERROR] Could not retrieve JES histogram " << Form("%s_1", h2D->GetName()) << std::endl;
        return;
    }
    set1DStyle(hJES, 2);
    hJES->GetYaxis()->SetTitle("JES");
    set1DStyle(hJES, type);

    // Jet energy resolution
    hJER = dynamic_cast<TH1D*>(gDirectory->Get( Form("%s_2", h2D->GetName()) ));
    if (!hJER) {
        std::cerr << "[ERROR] Could not retrieve JER histogram " << Form("%s_2", h2D->GetName()) << std::endl;
        hJES = nullptr; // Reset JES pointer since we failed to get JER
        return;
    }
    set1DStyle(hJER, 2);
    hJER->GetYaxis()->SetTitle("JER");
    set1DStyle(hJER, type);

    TF1 *jerFit;
    if ( performFit ) {
        jerFit = new TF1( Form("%s_fit", hJER->GetName()), "sqrt([0] + [1] / x)", 30., 600.);
        jerFit->SetParameters(0.002, 1.0);
        jerFit->SetLineColor( hJER->GetLineColor() );
        jerFit->SetLineWidth(2);
        hJER->Fit(jerFit, "MRE");
    }
}

//________________
void plotJESandJER(TCanvas *c, TH2 *h2D_trkMax, TH2 *h2D_noCut = nullptr, TH2 *h2D_jetId = nullptr,
                   TH1D *h1JES_trkMax = nullptr, TH1D *h1JES_noCut = nullptr, TH1D *h1JES_jetId = nullptr,
                   TH1D *h1JER_trkMax = nullptr, TH1D *h1JER_noCut = nullptr, TH1D *h1JER_jetId = nullptr,
                   double lowVal = -1.6, double hiVal = 1.6,
                   int collSystem = 0, double energy = 5.02,
                   bool isPt = true) {

    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.35;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);

    // Set xAxis range
    double xRange[2] = {10., 800.};
    double yMuRange[2] = {0.95, 1.05};
    double ySigmaRange[2] = {0., 0.25};
    
    if ( !isPt ) {
        xRange[0] = -4.0;
        xRange[1] = 4.0;
    }
    else {
        // h2D_trkMax->RebinX(2);
        // if ( h2D_noCut ) {
        //     h2D_noCut->RebinX(2);
        // }
        // if ( h2D_jetId ) {
        //     h2D_jetId->RebinX(2);
        // }
    }

    c->Clear();

    bool use2rows = {false};
    if ( !h2D_noCut && !h2D_jetId ) {
        c->Divide(3, 1);
    }
    else {
        c->Divide(3, 2);
        use2rows = {true};
    }

    // Plot 2D distribution
    c->cd(1);
    setPadStyle();
    h2D_trkMax->Draw("colz");
    h2D_trkMax->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h2D_trkMax->GetYaxis()->SetRangeUser(0., 2.);
    // gPad->SetLogz();
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, Form("TrkMax"));

    if ( use2rows ) {

        // Plot 2D distributions
        if ( h2D_noCut ) {
            c->cd(2);
            setPadStyle();
            h2D_noCut->Draw("colz");
            h2D_noCut->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            h2D_noCut->GetYaxis()->SetRangeUser(0., 2.);
            // gPad->SetLogz();
            plotCMSHeader(collSystem, energy);
            t.DrawLatexNDC(xTextPosition, yTextPosition, Form("No Cuts"));
        }

        if ( h2D_jetId ) {
            c->cd(3);
            setPadStyle();
            h2D_jetId->Draw("colz");
            h2D_jetId->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            h2D_jetId->GetYaxis()->SetRangeUser(0., 2.);
            // gPad->SetLogz();
            plotCMSHeader(collSystem, energy);
            t.DrawLatexNDC(xTextPosition, yTextPosition, Form("JetId"));
        }
    } // if ( use2rows)

    // Plot JES
    if ( !use2rows ) {
        c->cd(2);
    }
    else {
        c->cd(4);
    }
    setPadStyle();
    gPad->SetGrid(1, 1);
    h1JES_trkMax->Draw();
    if ( h1JES_noCut && h1JES_noCut->GetEntries() > 0 ) {
        h1JES_noCut->Draw("same");
    }
    if ( h1JES_jetId && h1JES_jetId->GetEntries() > 0 ) {
        h1JES_jetId->Draw("same");
    }
    h1JES_trkMax->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h1JES_trkMax->GetYaxis()->SetRangeUser(yMuRange[0], yMuRange[1]);
    if ( isPt ) {
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%2.1f < #eta^{jet} < %2.1f ", lowVal, hiVal) );
    }
    else {
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%4.0f < p_{T}^{jet} (GeV) < %4.0f ", lowVal, hiVal) );
    }
    plotCMSHeader(collSystem, energy);

    // Plot JER
    if ( !use2rows ) {
        c->cd(3);
    }
    else {
        c->cd(5);
    }
    setPadStyle();
    gPad->SetGrid(1, 1);
    h1JER_trkMax->Draw();
    if ( h1JER_noCut && h1JER_noCut->GetEntries() > 0 ) {
        h1JER_noCut->Draw("same");
    }
    if ( h1JER_jetId && h1JER_jetId->GetEntries() > 0 ) {
        h1JER_jetId->Draw("same");
    }
    h1JER_trkMax->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h1JER_trkMax->GetYaxis()->SetRangeUser(ySigmaRange[0], ySigmaRange[1]);
    plotCMSHeader(collSystem, energy);

    if ( use2rows ) {
        c->cd(6);
        setPadStyle();
        TLegend *leg = new TLegend(0.3, 0.35, 0.7, 0.75);
        leg->SetTextFont(42);
        leg->SetTextSize(0.05);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        if ( isPt ) {
            leg->SetHeader(Form("%2.1f < #eta^{jet} < %2.1f ", lowVal, hiVal));
        }
        else {
            leg->SetHeader(Form("%4.0f < p_{T}^{jet} (GeV) < %4.0f ", lowVal, hiVal));
        }
        leg->AddEntry(h1JER_trkMax, "trkMax", "p");
        if ( h1JER_noCut && h1JER_noCut->GetEntries() > 0 ) {
            leg->AddEntry(h1JER_noCut, "noCut", "p");
        }
        if ( h1JER_jetId && h1JER_jetId->GetEntries() > 0 ) {
            leg->AddEntry(h1JER_jetId, "jetId", "p");
        }
        leg->Draw();
    }
}

//________________
void processJESandJER(TFile *fTrkMax, TFile* fNoCut = nullptr, TFile *fJetId = nullptr,
                   int collSystem = 1, double energy = 8.16, TString date = "20250129") {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    bool performFit = true;

    // fNoCut = nullptr; // Disable noCut for now
    // fJetId = nullptr; // Disable jetId for now

    // Retrieve JES hInclusiveJetJESGenPtGenEtaPtHatWeighted hInclusiveJetJESRecoPtRecoEtaPtHatWeighted
    THnSparseD *hJESPars_trkMax = (THnSparseD*)fTrkMax->Get("hInclusiveJetJESGenPtGenEtaPtHatWeighted");
    hJESPars_trkMax->SetName("hJESPars_trkMax");
    if ( hJESPars_trkMax == nullptr ) {
        std::cerr << Form( "[ERROR] Could not retrieve hInclusiveJetJESGenPtGenEtaPtHatWeighted from the file: %s", fTrkMax->GetName() ) << std::endl;
        return;
    }

    THnSparseD *hJESPars_noCut = nullptr;
    if ( fNoCut ) {
        hJESPars_noCut = (THnSparseD*)fNoCut->Get("hInclusiveJetJESGenPtGenEtaPtHatWeighted");
        hJESPars_noCut->SetName("hJESPars_noCut");
        if ( hJESPars_noCut == nullptr ) {
            std::cerr << Form( "[ERROR] Could not retrieve hInclusiveJetJESGenPtGenEtaPtHatWeighted from the file: %s", fNoCut->GetName() ) << std::endl;
            return;
        }
    }

    THnSparseD *hJESPars_jetId = nullptr;
    if ( fJetId ) {
        hJESPars_jetId = (THnSparseD*)fJetId->Get("hInclusiveJetJESGenPtGenEtaPtHatWeighted");
        hJESPars_jetId->SetName("hJESPars_jetId");
        if ( hJESPars_jetId == nullptr ) {
            std::cerr << Form( "[ERROR] Could not retrieve hInclusiveJetJESGenPtGenEtaPtHatWeighted from the file: %s", fJetId->GetName() ) << std::endl;
            return;
        }
    }

    // Event ptHat binning
    double ptHatVals[] = { 15.};
    int sizeOfPtHatBins = sizeof(ptHatVals) / sizeof(ptHatVals[0]);

    // Jet eta binning
    // double etaVals[] = { -1.6, 1.6 }; // Simplified for this example
    double etaVals[] = { -5.2, -3.6, -3.0, -2.8, -2.6, -2.4, -2.0, -1.6, 0., 1.6, 2.0, 2.4, 2.6, 2.8, 3.0, 3.6, 5.2 };
    int sizeOfEtaBins = sizeof(etaVals) / sizeof(etaVals[0]);

    // Jet pT binning
    double ptVals[] = { 25., 55., 95., 125., 200., 800. };
    int sizeOfPtBins = sizeof(ptVals) / sizeof(ptVals[0]);

    std::cout << Form("Number of ptHat bins: %d, Number of eta bins: %d, Number of pT bins: %d", 
                      sizeOfPtHatBins, sizeOfEtaBins, sizeOfPtBins) << std::endl;

    // Create 2D and 1D histograms
    TH2D *h2D_trkMax{nullptr}, *h2D_noCut{nullptr}, *h2D_jetId{nullptr};
    TH1D *h1JES_trkMax{nullptr}, *h1JES_noCut{nullptr}, *h1JES_jetId{nullptr};
    TH1D *h1JER_trkMax{nullptr}, *h1JER_noCut{nullptr}, *h1JER_jetId{nullptr};

    TCanvas *c = nullptr;
    if ( !fNoCut && !fJetId ) {
        c = new TCanvas( "c", "cJESandJER", 1500, 500 );
    }
    else {
        c = new TCanvas( "c", "cJESandJER", 1500, 1000 );
    }

    // Loop over ptHat bins
    for (unsigned int i{0}; i<sizeOfPtHatBins; i++) {

        double ptHatLow = ptHatVals[i];
        int ptHatBinLow = hJESPars_trkMax->GetAxis(3)->FindBin(ptHatLow);

        std::cout << Form("ptHat low: %4f hi: %4f\n", 
                          hJESPars_trkMax->GetAxis(3)->GetBinLowEdge(ptHatBinLow),
                          hJESPars_trkMax->GetAxis(3)->GetBinUpEdge( hJESPars_trkMax->GetAxis(3)->GetNbins() )  );

        // Set ptHat binning
        hJESPars_trkMax->GetAxis(3)->SetRange( ptHatBinLow, hJESPars_trkMax->GetAxis(3)->GetNbins() );

        if ( hJESPars_noCut ) {
            hJESPars_noCut->GetAxis(3)->SetRange( ptHatBinLow, hJESPars_noCut->GetAxis(3)->GetNbins() );
        }

        if ( hJESPars_jetId ) {
            hJESPars_jetId->GetAxis(3)->SetRange( ptHatBinLow, hJESPars_jetId->GetAxis(3)->GetNbins() );
        }

        //
        // JES vs. jet pT for different eta bins
        //
        for (unsigned int j{0}; j<(sizeOfEtaBins-1); j++) {

            h2D_trkMax = nullptr; h2D_noCut = nullptr; h2D_jetId = nullptr;
            h1JES_trkMax = nullptr; h1JES_noCut = nullptr; h1JES_jetId = nullptr;
            h1JER_trkMax = nullptr; h1JER_noCut = nullptr; h1JER_jetId = nullptr;

            double etaLow = etaVals[j];
            double etaHi = etaVals[j+1];
            int etaBinLow = hJESPars_trkMax->GetAxis(2)->FindBin(etaLow);
            int etaBinHigh = hJESPars_trkMax->GetAxis(2)->FindBin(etaHi) - 1;

            // std::cout << Form("eta low: %.2f hi: %.2f for bins (%d, %d)", 
            //               hJESPars_trkMax->GetAxis(2)->GetBinLowEdge(etaBinLow),
            //               hJESPars_trkMax->GetAxis(2)->GetBinUpEdge( etaBinHigh ),
            //               etaBinLow, etaBinHigh ) << std::endl;

            // Set axis limits for THnSparse
            hJESPars_trkMax->GetAxis(2)->SetRange( etaBinLow, etaBinHigh );

            // Create 2D histogram
            // if ( h2D_trkMax ) { h2D_trkMax->Clear(); }
            h2D_trkMax = dynamic_cast<TH2D*>( hJESPars_trkMax->Projection(0, 1) );
            h2D_trkMax->SetName( Form("hJES_trkMax_pt_%d_%d", i, j) );
            set2DStyle(h2D_trkMax);

            // Retrieve JES and JER
            retrieveJESandJER(h2D_trkMax, h1JES_trkMax, h1JER_trkMax, 0, performFit);

            if ( fNoCut ) {
                // Set axis limits for THnSparse
                hJESPars_noCut->GetAxis(2)->SetRange( etaBinLow, etaBinHigh );

                // Create 2D histogram
                h2D_noCut = dynamic_cast<TH2D*>( hJESPars_noCut->Projection(0, 1) );
                h2D_noCut->SetName( Form("hJES_noCut_pt_%d_%d", i, j) );
                set2DStyle(h2D_noCut);

                // Retrieve JES and JER
                retrieveJESandJER(h2D_noCut, h1JES_noCut, h1JER_noCut, 1, performFit);

                // Protect against empty histograms
                if ( h2D_noCut == nullptr || h2D_noCut->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h2D_noCut = nullptr;
                }
                if ( h1JES_noCut == nullptr || h1JES_noCut->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JES_noCut = nullptr;
                }
                if ( h1JER_noCut == nullptr || h1JER_noCut->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JER_noCut = nullptr;
                }
            } // if ( fNoCut ) 

            if ( fJetId ) {
                // Set axis limits for THnSparse
                hJESPars_jetId->GetAxis(2)->SetRange( etaBinLow, etaBinHigh );

                // Create 2D histogram
                h2D_jetId = dynamic_cast<TH2D*>( hJESPars_jetId->Projection(0, 1) );
                h2D_jetId->SetName( Form("hJES_jetId_pt_%d_%d", i, j) );
                set2DStyle(h2D_jetId);

                // Retrieve JES and JER
                retrieveJESandJER(h2D_jetId, h1JES_jetId, h1JER_jetId, 2, performFit);

                // Protect against empty histograms
                if ( h2D_jetId == nullptr || h2D_jetId->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h2D_jetId = nullptr;
                }
                if ( h1JES_jetId == nullptr || h1JES_jetId->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JES_jetId = nullptr;
                }
                if ( h1JER_jetId == nullptr || h1JER_jetId->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JER_jetId = nullptr;
                }
            } // if ( fJetId ) 

            plotJESandJER(c, h2D_trkMax, h2D_noCut, h2D_jetId,
                          h1JES_trkMax, h1JES_noCut, h1JES_jetId,
                          h1JER_trkMax, h1JER_noCut, h1JER_jetId,
                          etaLow, etaHi, collSystem, energy, true);

            c->SaveAs( Form("%s/JES_vs_pt_ptHat_%d_eta_%s_%s.pdf", date.Data(), (int)ptHatLow, 
                        ((etaLow>=0) ? Form("%.1f", etaLow) : Form("m%.1f", -etaLow)), 
                        ((etaHi>=0) ? Form("%.1f", etaHi) : Form("m%.1f", -etaHi)) ) );

        } // for (unsigned int j{0}; j<(sizeOfEtaBins-1); j++)

        // Restore eta binning
        hJESPars_trkMax->GetAxis(2)->SetRange( 1, hJESPars_trkMax->GetAxis(2)->GetNbins() ); 
        if ( hJESPars_noCut ) {
            hJESPars_noCut->GetAxis(2)->SetRange( 1, hJESPars_noCut->GetAxis(2)->GetNbins() );
        }
        if ( hJESPars_jetId ) {
            hJESPars_jetId->GetAxis(2)->SetRange( 1, hJESPars_jetId->GetAxis(2)->GetNbins() );
        }

        //
        // JES vs. jet eta for different pt bins
        //
        for (unsigned int j{0}; j<(sizeOfPtBins-1); j++) {

            h2D_trkMax = nullptr; h2D_noCut = nullptr; h2D_jetId = nullptr;
            h1JES_trkMax = nullptr; h1JES_noCut = nullptr; h1JES_jetId = nullptr;
            h1JER_trkMax = nullptr; h1JER_noCut = nullptr; h1JER_jetId = nullptr;

            double ptLow = ptVals[j];
            double ptHi = ptVals[j+1];
            int ptBinLow = hJESPars_trkMax->GetAxis(1)->FindBin(ptLow);
            int ptBinHigh = hJESPars_trkMax->GetAxis(1)->FindBin(ptHi) - 1;

            // std::cout << Form("pt low: %.2f hi: %.2f for bins (%d, %d)", 
            //               hJESPars_trkMax->GetAxis(1)->GetBinLowEdge(ptBinLow),
            //               hJESPars_trkMax->GetAxis(1)->GetBinUpEdge( ptBinHigh ),
            //               ptBinLow, ptBinHigh ) << std::endl;

            // Set axis limits for THnSparse
            hJESPars_trkMax->GetAxis(1)->SetRange( ptBinLow, ptBinHigh );

            // Create 2D histogram
            h2D_trkMax = dynamic_cast<TH2D*>( hJESPars_trkMax->Projection(0, 2) );
            h2D_trkMax->SetName( Form("hJES_trkMax_eta_%d_%d", i, j) );
            set2DStyle(h2D_trkMax);

            // Retrieve JES and JER
            retrieveJESandJER(h2D_trkMax, h1JES_trkMax, h1JER_trkMax, 0, false);

            if ( fNoCut ) {
                // Set axis limits for THnSparse
                hJESPars_noCut->GetAxis(1)->SetRange( ptBinLow, ptBinHigh );

                // Create 2D histogram
                h2D_noCut = dynamic_cast<TH2D*>( hJESPars_noCut->Projection(0, 2) );
                h2D_noCut->SetName( Form("hJES_noCut_eta_%d_%d", i, j) );
                set2DStyle(h2D_noCut);

                // Retrieve JES and JER
                retrieveJESandJER(h2D_noCut, h1JES_noCut, h1JER_noCut, 1, false);
                
                // Protection from empty histograms
                if ( h2D_noCut == nullptr || h2D_noCut->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h2D_noCut = nullptr;
                }
                if ( h1JES_noCut == nullptr || h1JES_noCut->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JES_noCut = nullptr;
                }
                if ( h1JER_noCut == nullptr || h1JER_noCut->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JER_noCut = nullptr;
                }
            } // if ( fNoCut )

            if ( fJetId ) {
                // Set axis limits for THnSparse
                hJESPars_jetId->GetAxis(1)->SetRange( ptBinLow, ptBinHigh );

                // Create 2D histogram
                h2D_jetId = dynamic_cast<TH2D*>( hJESPars_jetId->Projection(0, 2) );
                h2D_jetId->SetName( Form("hJES_jetId_eta_%d_%d", i, j) );
                set2DStyle(h2D_jetId);

                // Retrieve JES and JER
                retrieveJESandJER(h2D_jetId, h1JES_jetId, h1JER_jetId, 2, false);

                // Protection from empty histograms
                if ( h2D_jetId == nullptr || h2D_jetId->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h2D_jetId = nullptr;
                }
                if ( h1JES_jetId == nullptr || h1JES_jetId->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JES_jetId = nullptr;
                }
                if ( h1JER_jetId == nullptr || h1JER_jetId->GetEntries() == 0 ) {
                    // Handle empty histogram case
                    h1JER_jetId = nullptr;
                }
            } // if ( fJetId )

            plotJESandJER(c, h2D_trkMax, h2D_noCut, h2D_jetId,
                          h1JES_trkMax, h1JES_noCut, h1JES_jetId,
                          h1JER_trkMax, h1JER_noCut, h1JER_jetId,
                          ptLow, ptHi, collSystem, energy, false);
            c->SaveAs( Form("%s/JES_vs_eta_ptHat_%d_pt_%d_%d.pdf", date.Data(), (int)ptHatLow, (int)ptLow, (int)ptHi) );
        } // for (unsigned int j{0}; j<(sizeOfPtBins-1); j++)

    } // for (unsigned int i{0}; i<sizeOfPtHatBins; i++)

    if ( c ) { delete c; c = nullptr; }
}

//________________
void plotSimpleJES(TFile *f, int collSystem = 0, double energy = 5.02) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    // The function reads the JES (reco pt/gen pt vs. gen pt for |eta|<1.4) from the file

    double xTextPosition = 0.7;
    double yTextPosition = 0.8;
    TLatex t;

    // Retrieve JES
    TH2D *h2D = dynamic_cast<TH2D*>( f->Get("hInclusiveJetJESVsPtGen") );
    set2DStyle(h2D);
    h2D->FitSlicesY();
    // Jet energy scale
    TH1D *hJESMean = (TH1D*)gDirectory->Get("hInclusiveJetJESVsPtGen_1");
    hJESMean->SetName("hJESMean");
    set1DStyle(hJESMean, 2);
    hJESMean->GetYaxis()->SetTitle("JES");

    // Jet energy resolution
    TH1D *hJESSigma = (TH1D*)gDirectory->Get("hInclusiveJetJESVsPtGen_2");
    hJESSigma->SetName("hJESSigma");
    set1DStyle(hJESSigma, 2);
    hJESSigma->GetYaxis()->SetTitle("JER");

    //
    // Create canvas
    //
    TCanvas *cSimpleJES = new TCanvas( "cSimpleJES", "cSimpleJES", 1500, 500 );
    cSimpleJES->Divide(3, 1);

    cSimpleJES->cd(1);
    setPadStyle();
    h2D->Draw("colz");
    //gPad->SetLogz();
    h2D->GetXaxis()->SetRangeUser(0, 800);
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");

    cSimpleJES->cd(2);
    setPadStyle();
    hJESMean->Draw();
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");
    hJESMean->GetXaxis()->SetRangeUser(0, 800);
    hJESMean->GetYaxis()->SetRangeUser(0.98, 1.3);

    cSimpleJES->cd(3);
    setPadStyle();
    hJESSigma->Draw();
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");
    hJESSigma->GetXaxis()->SetRangeUser(0, 800);
    hJESSigma->GetYaxis()->SetRangeUser(-0.005, 0.17);
}

//________________
void plotJES() {

    // Base style
    // gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);


    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Date extraction
    TDatime dt; 
    TString date { Form( "%d",dt.GetDate() ) };

    if ( !directoryExists( date.Data() ) ) {
        createDirectory( date.Data() );
    } 

    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    int collisionSystem = 1;         // 0 - pp, 1 - pPb, 2 - pPb5020, 3 - pPb8160
    double collisionEnergy = 8.16;   // 8.16 TeV
    int direction = 1;               // 0-p-going, 1-Pb-going, 2 - combined
    TString directionStr = (direction == 0) ? "pgoing" : ((direction == 1) ? "Pbgoing" : "");
    int dataTrigger = 0;               // 0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100
    TString dataStr = (dataTrigger == 0) ? "MB" : ((dataTrigger == 1) ? "Jet60" : ((dataTrigger == 2) ? "Jet80" : ((dataTrigger == 3) ? "Jet100" : "unknownData")));
    TString dataDirectionStr = (direction == 0) ? "Pbgoing" : ((direction == 1) ? "pgoing" : "");
    int genType = 0; // 0 - embedding, 1 - pythia
    TString mcTypeStr = (genType == 0) ? "Embedding" : ((genType == 1) ? "Pythia" : "unknownGen");
    TString mcTypeStrLower = mcTypeStr;
    mcTypeStrLower.ToLower();

    // // pp5020 PYTHIA
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", username.Data()) );
    // if ( !inputFile ) {
    //     std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", username.Data()) << std::endl;
    //     return;
    // }
    // collSystem = 0;
    // energy = 5.02;

    //
    // pPb8160 embedding
    //

    TFile *trkMaxFile = nullptr;
    if ( direction < 2 ) {
        trkMaxFile = TFile::Open( Form("~/cernbox/ana/pPb8160/%s/%s/o%s_%s_def_ak4_eta20.root", 
                                    mcTypeStrLower.Data(), directionStr.Data(), mcTypeStr.Data(), 
                                    directionStr.Data()) );
        if ( !trkMaxFile || trkMaxFile->IsZombie() ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/%s/%s/o%s_%s_def_ak4_eta20.root", 
                              mcTypeStrLower.Data(), directionStr.Data(), mcTypeStr.Data(), 
                              directionStr.Data()) << std::endl;
            return;
        }
    }
    else {
        trkMaxFile = TFile::Open( Form("~/cernbox/ana/pPb8160/%s/o%s_pPb8160_def_ak4_eta20.root", 
                                     mcTypeStrLower.Data(), mcTypeStr.Data()) );
        if ( !trkMaxFile || trkMaxFile->IsZombie() ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/%s/o%s_pPb8160_def_ak4_eta20.root", 
                              mcTypeStrLower.Data(), mcTypeStr.Data()) << std::endl;
            return;
        }
    }

    TFile *noCutFile = nullptr;
    if ( direction < 2 ) {
        noCutFile = TFile::Open( Form("~/cernbox/ana/pPb8160/%s/%s/o%s_%s_def_ak4_noCut_eta20.root", 
            mcTypeStrLower.Data(), directionStr.Data(), mcTypeStr.Data(), directionStr.Data()) );
        if ( !noCutFile || noCutFile->IsZombie() ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/%s/%s/o%s_%s_def_ak4_noCut_eta20.root", 
                              mcTypeStrLower.Data(), directionStr.Data(), mcTypeStr.Data(), directionStr.Data()) << std::endl;
            return;
        }
    }
    else {
        noCutFile = TFile::Open( Form("~/cernbox/ana/pPb8160/%s/o%s_pPb8160_def_ak4_noCut_eta20.root",
                                         mcTypeStrLower.Data(), mcTypeStr.Data()) );
        if ( !noCutFile || noCutFile->IsZombie() ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/%s/o%s_pPb8160_def_ak4_noCut_eta20.root",
                              mcTypeStrLower.Data(), mcTypeStr.Data()) << std::endl;
            return;
        }
    }

    TFile *jetIdFile = nullptr;
    if ( direction < 2 ) {
        jetIdFile = TFile::Open( Form("~/cernbox/ana/pPb8160/%s/%s/o%s_%s_def_ak4_jetId_eta20.root", 
                                    mcTypeStrLower.Data(), directionStr.Data(), mcTypeStr.Data(), directionStr.Data()) );
        if ( !jetIdFile || jetIdFile->IsZombie() ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/%s/%s/o%s_%s_def_ak4_jetId_eta20.root", 
                              mcTypeStrLower.Data(), directionStr.Data(), mcTypeStr.Data(), directionStr.Data()) << std::endl;
            return;
        }
    }
    else {
        jetIdFile = TFile::Open( Form("~/cernbox/ana/pPb8160/%s/o%s_pPb8160_def_ak4_jetId_eta20.root", 
                                         mcTypeStrLower.Data(), mcTypeStr.Data()) );
        if ( !jetIdFile || jetIdFile->IsZombie() ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/%s/o%s_pPb8160_def_ak4_jetId_eta20.root", 
                              mcTypeStrLower.Data(), mcTypeStr.Data()) << std::endl;
            return;
        }
    }

    if ( trkMaxFile ) {
        std::cout << Form("TrkMax file opened: %s", trkMaxFile->GetName()) << std::endl;
    }
    if ( noCutFile ) {
        std::cout << Form("NoCut file opened: %s", noCutFile->GetName()) << std::endl;
    }
    if ( jetIdFile ) {
        std::cout << Form("JetId file opened: %s", jetIdFile->GetName()) << std::endl;
    }

    //
    // Plot simple JES
    //
    // plotSimpleJES( trkMaxFile, collisionSystem, collisionEnergy );

    //
    // Plot JES and JER for different ptHat, eta, and pT selections
    //
    processJESandJER( trkMaxFile, noCutFile, jetIdFile, collisionSystem, collisionEnergy, date );
    // processJESandJER( trkMaxFile, noCutFile, nullptr, collisionSystem, collisionEnergy, date );
}
