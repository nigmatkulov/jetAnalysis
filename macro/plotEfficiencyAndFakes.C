// ROOT headers
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TLine.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"

// C++ headers
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
void plotCMSHeader(int collSystem = 0, double energy = 5.02) {
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
void setRatioPadStyle(TPad *pad, bool isUpper = true) {
    if ( isUpper ) {
        pad->SetTopMargin(0.1);
        pad->SetBottomMargin(0.2);
        pad->SetRightMargin(0.1);
        pad->SetLeftMargin(0.2);
        pad->SetFrameLineWidth(2);
    }
    else {
        pad->SetTopMargin(0.1);
        pad->SetBottomMargin(0.2);
        pad->SetRightMargin(0.1);
        pad->SetLeftMargin(0.2);
        pad->SetFrameLineWidth(2);
    }
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 0.9;
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

    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);    
    h->GetYaxis()->SetTitleOffset(1.2);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void setGraphStyle(TGraph *gr, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 0.9;
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

    gr->SetLineWidth( lineWidth );
    gr->SetLineColor( color );
    
    gr->SetMarkerStyle( markerStyle );
    gr->SetMarkerColor( color );
    gr->SetMarkerSize( markerSize );

    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetNdivisions(205);
    gr->GetYaxis()->SetNdivisions(205);    
    gr->GetYaxis()->SetTitleOffset(1.1);
}

//________________
void set2DStyle(TH2* h) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(210);
    h->GetYaxis()->SetNdivisions(210);    
    h->GetYaxis()->SetTitleOffset(1.0);
}

//________________
void plotEfficiency(TFile *f, int collSystem = 0, double energy = 5.02, int jetType = 0, TString date = "20250606") {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    // jetType: 0 = inclusive, 1 = lead, 2 = sublead

    TString inputFileName( f->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "combined";
    }

    double xTextPosition = 0.3;
    double yTextPosition = 0.85;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );

    int canvXSize = 800;
    int canvYSize = 800;

    const int nJetTypes = 3; // Inclusive, Lead, SubLead
    const char* jetTypes[3] = {"Inclusive", "Lead", "SubLead"};
    TH3D *hRefPtEtaPtHat[nJetTypes] = {nullptr};
    TH3D *hGenPtEtaPtHat[nJetTypes] = {nullptr};
    TH2D *hEfficiency2D[nJetTypes] = {nullptr};
    TH2D *hRef2D[nJetTypes] = {nullptr};
    TH2D *hGen2D[nJetTypes] = {nullptr};

    TCanvas *cEfficiency2D = new TCanvas( "cEfficiency2D", "cEfficiency2D", canvXSize, canvYSize );
    setPadStyle();

    // Read all histograms from the file and check if they exist
    for (size_t i = 0; i < nJetTypes; ++i) {
        TString tJetType = jetTypes[i];

        hRefPtEtaPtHat[i] = dynamic_cast<TH3D*>( f->Get( Form("hRef%sJetPtEtaPtHat", tJetType.Data()) ) );
        if ( !hRefPtEtaPtHat[i] ) {
            std::cerr << "Histogram hRef" << tJetType.Data() << "JetPtEtaPtHat not found in file." << std::endl;
            return;
        }
        hGenPtEtaPtHat[i] = dynamic_cast<TH3D*>( f->Get( Form("hGen%sJetPtEtaPtHat", tJetType.Data()) ) );
        if ( !hGenPtEtaPtHat[i] ) {
            std::cerr << "Histogram hGen" << tJetType.Data() << "JetPtEtaPtHat not found in file." << std::endl;
            return;
        }

        //
        // Create 2D efficiency histogram
        //
        hEfficiency2D[i] = dynamic_cast<TH2D*>( hRefPtEtaPtHat[i]->Project3D("yx") );
        if ( !hEfficiency2D[i] ) {
            std::cerr << "Failed to project hRef" << tJetType.Data() << "JetPtEtaPtHat to 2D." << std::endl;
            return;
        }
        hEfficiency2D[i]->SetName( Form("hEfficiency2D_%s", tJetType.Data()) );

        hRef2D[i] = dynamic_cast<TH2D*>( hRefPtEtaPtHat[i]->Project3D("yx") );
        if ( !hRef2D[i] ) {
            std::cerr << "Failed to project hRef" << tJetType.Data() << "JetPtEtaPtHat to 2D." << std::endl;
            return;
        }
        hRef2D[i]->SetName( Form("hRef2D_%s", tJetType.Data()) );

        hGen2D[i] = dynamic_cast<TH2D*>( hGenPtEtaPtHat[i]->Project3D("yx") );
        if ( !hGen2D[i] ) {
            std::cerr << "Failed to project hGen" << tJetType.Data() << "JetPtEtaPtHat to 2D." << std::endl;
            return;
        }
        hGen2D[i]->SetName( Form("hGen2D_%s", tJetType.Data()) );

        // std::cout << Form("%s refJet integral: %f, genJet integral: %f", 
        //                   tJetType.Data(), hRef2D[i]->Integral(), hGen2D[i]->Integral() ) << std::endl;

        // Create 2D efficiency histogram
        hEfficiency2D[i]->Divide( hGen2D[i] );

        cEfficiency2D->Clear();
        hEfficiency2D[i]->Draw("colz");
        hEfficiency2D[i]->GetXaxis()->SetRangeUser(-3.6, 3.6);
        hEfficiency2D[i]->GetYaxis()->SetRangeUser(20., 400.);
        hEfficiency2D[i]->SetMaximum(1.0);
        hEfficiency2D[i]->SetMinimum(0.85);
        hEfficiency2D[i]->GetXaxis()->SetTitle("#eta");
        hEfficiency2D[i]->GetYaxis()->SetTitle("p_{T} (GeV)");
        plotCMSHeader(collSystem, energy);
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%s jet efficiency", tJetType.Data()));

        tJetType.ToLower();
        cEfficiency2D->SaveAs( Form("%s/%s%d_%s_efficiency2D_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")), 
                                    int(energy * 1000), direction.Data(), tJetType.Data() ) );
    } // for (size_t i = 0; i < nJetTypes; ++i)

    // Projections in eta and pt to look at fakes
    const int nEtaRanges = 10;  // Number of eta ranges (9 ranges, 10 boundaries)
    double etaRanges[nEtaRanges] = {-5.0, -3.0, -2.4, -2.0, -1.0, 1.0, 2.0, 2.4, 3.0, 5.0};
    const int nPtRanges = 8;   // Number of pT ranges (7 ranges, 8 boundaries)
    double ptRanges[nPtRanges] = {5., 25., 35., 45., 65., 95., 150., 1505.};

    // Efficiency as a function of pT for each jet type and eta range
    TH1D *hEfficiencyPt[nJetTypes][nEtaRanges - 1];
    TH1D *hGenPt[nJetTypes][nEtaRanges - 1];
    // Efficiency as a function of eta for each jet type and pT range
    TH1D *hEfficiencyEta[nJetTypes][nPtRanges - 1];
    TH1D *hGenEta[nJetTypes][nPtRanges - 1];

    TCanvas *cEfficiencyIndividual = new TCanvas( "cEfficiencyIndividual", "cEfficiencyIndividual", canvXSize, canvYSize );
    setPadStyle();
    TCanvas *cEfficiencyCombined = new TCanvas( "cEfficiencyCombined", "cEfficiencyCombined", canvXSize, canvYSize );
    setPadStyle();

    // Plot projections of matched and fakes on eta as a function of pT
    for (size_t i = 0; i < nEtaRanges - 1; ++i) {
        double etaMin = etaRanges[i];
        double etaMax = etaRanges[i + 1];

        for (size_t j = 0; j < nJetTypes; ++j) {
            TString tJetType = jetTypes[j];
            tJetType.ToLower();

            // std::cout << Form("%s refJet integral in eta range %.1f-%.1f: %f, genJet integral: %f", 
            //                   tJetType.Data(), etaMin, etaMax, 
            //                   hRef2D[j]->Integral(hRef2D[j]->GetXaxis()->FindBin(etaMin), hRef2D[j]->GetXaxis()->FindBin(etaMax), 1, 
            //                                       hRef2D[j]->GetYaxis()->GetNbins()), 
            //                   hGen2D[j]->Integral(hGen2D[j]->GetXaxis()->FindBin(etaMin), hGen2D[j]->GetXaxis()->FindBin(etaMax), 1, 
            //                                       hGen2D[j]->GetYaxis()->GetNbins()) ) << std::endl;

            // Project efficiency in eta for the current range
            hEfficiencyPt[j][i] = dynamic_cast<TH1D*>(hRef2D[j]->ProjectionY( Form("hEfficiencyPt_%d_%d", int(j), int(i)), 
                                                      hRef2D[j]->GetXaxis()->FindBin(etaMin), 
                                                      hRef2D[j]->GetXaxis()->FindBin(etaMax) ));
            if ( !hEfficiencyPt[j][i] ) {
                std::cerr << "Failed to project hEfficiency2D" << tJetType.Data() << " to 1D." << std::endl;
                return;
            }
            set1DStyle(hEfficiencyPt[j][i], j);
            hEfficiencyPt[j][i]->GetXaxis()->SetTitle("p_{T} (GeV)");
            hEfficiencyPt[j][i]->GetYaxis()->SetTitle("Efficiency");

            hGenPt[j][i] = dynamic_cast<TH1D*>(hGen2D[j]->ProjectionY( Form("hGenPt_%d_%d", int(j), int(i)), 
                                               hGen2D[j]->GetXaxis()->FindBin(etaMin), 
                                               hGen2D[j]->GetXaxis()->FindBin(etaMax) ));
            if ( !hGenPt[j][i] ) {
                std::cerr << "Failed to project hGen2D" << tJetType.Data() << " to 1D." << std::endl;
                return;
            }

            // std::cout << Form("%s efficiency in eta range %.1f-%.1f ref: %f gen: %f", 
            //                   tJetType.Data(), etaMin, etaMax, 
            //                   hEfficiencyPt[j][i]->Integral(), hGenPt[j][i]->Integral()) << std::endl;

            //
            hEfficiencyPt[j][i]->Divide(hEfficiencyPt[j][i], hGenPt[j][i], 1., 1., "B"); 
            // std::cout << Form("%s efficiency in eta range %.1f-%.1f after division: %f", 
            //                   tJetType.Data(), etaMin, etaMax, 
            //                   hEfficiencyPt[j][i]->Integral() ) << std::endl;

            cEfficiencyIndividual->Clear();
            cEfficiencyIndividual->cd();
            setPadStyle();
            hEfficiencyPt[j][i]->Draw();
            hEfficiencyPt[j][i]->GetXaxis()->SetRangeUser(15., 400.);
            hEfficiencyPt[j][i]->GetYaxis()->SetRangeUser(0.65, 1.01);
            hEfficiencyPt[j][i]->GetYaxis()->SetTitle("Efficiency");
            plotCMSHeader(collSystem, energy);
            t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%.1f < #eta < %.1f", etaMin, etaMax));
            cEfficiencyIndividual->SaveAs( Form("%s/%s%d_%s_%s_efficiency_eta_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                        int(energy * 1000), direction.Data(), tJetType.Data(),
                                        ((etaMin < 0) ? Form("m%.1f", -etaMin) : Form("%.1f", etaMin)), 
                                        ((etaMax < 0) ? Form("m%.1f", -etaMax) : Form("%.1f", etaMax)) ) );

            // Combine all jet types in eta range
            cEfficiencyCombined->cd();
            if (j == 0) {
                cEfficiencyCombined->Clear();
                hEfficiencyPt[j][i]->Draw();
                hEfficiencyPt[j][i]->GetXaxis()->SetRangeUser(15., 400.);
                hEfficiencyPt[j][i]->GetYaxis()->SetRangeUser(0.35, 1.3);
                hEfficiencyPt[j][i]->GetYaxis()->SetTitle("Efficiency");
                plotCMSHeader(collSystem, energy);
                t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%.1f < #eta < %.1f", etaMin, etaMax));
            } 
            else {
                hEfficiencyPt[j][i]->Draw("same");
            }
        } // for (size_t j = 0; j < nJetTypes; ++j)

        // Finalize combined efficiency plot
        cEfficiencyCombined->SaveAs( Form("%s/%s%d_%s_combined_efficiency_eta_%s_%s.pdf", 
                                        date.Data(), 
                                        (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                        int(energy * 1000), direction.Data(),
                                        ((etaMin < 0) ? Form("m%.1f", -etaMin) : Form("%.1f", etaMin)),
                                        ((etaMax < 0) ? Form("m%.1f", -etaMax) : Form("%.1f", etaMax)) ) );
    } // for (size_t i = 0; i < nEtaRanges - 1; ++i)


    // Clean up
    for (size_t i = 0; i < nJetTypes; ++i) {
        if (hRefPtEtaPtHat[i]) delete hRefPtEtaPtHat[i];
        if (hGenPtEtaPtHat[i]) delete hGenPtEtaPtHat[i];
        if (hEfficiency2D[i]) delete hEfficiency2D[i];
        if (hRef2D[i]) delete hRef2D[i];
        if (hGen2D[i]) delete hGen2D[i];
        for (size_t j = 0; j < nEtaRanges - 1; ++j) {
            if (hEfficiencyPt[i][j]) delete hEfficiencyPt[i][j];
            if (hGenPt[i][j]) delete hGenPt[i][j];
        }
    }
    if (cEfficiencyIndividual) delete cEfficiencyIndividual;
    if (cEfficiencyCombined) delete cEfficiencyCombined;
    if (cEfficiency2D) delete cEfficiency2D;
}

//________________
void plotFakes(TFile *f, int collSystem = 0, double energy = 5.02, int jetType = 0, TString date = "20250606") {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    // jetType: 0 = inclusive, 1 = lead, 2 = sublead

    TString inputFileName( f->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "combined";
    }

    double xTextPosition = 0.4;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );

    const int nJetTypes = 3; // Inclusive, Lead, SubLead
    const char* jetTypes[3] = {"Inclusive", "Lead", "SubLead"};
    const int nHistTitles = 3; // AllJetPtEta, MatchedJetPtEta, UnmatchedJetPtEta
    const char* histTitles[3] = {"AllJetPtEta", "MatchedJetPtEta", "UnmatchedJetPtEta"};

    // Histograms to read from file
    TH2D *hRecoAllJetPtEta[nJetTypes] = {nullptr};
    TH2D *hRecoMatchedJetPtEta[nJetTypes] = {nullptr};
    TH2D *hRecoUnmatchedJetPtEta[nJetTypes] = {nullptr};

    // Ratios
    TH2D *hRecoMatchedFraction[nJetTypes] = {nullptr};
    TH2D *hRecoUnmatchedFraction[nJetTypes] = {nullptr};

    TCanvas *cFakes2D = new TCanvas( "cFakes2D", "cFakes2D", 800, 800 );
    setPadStyle();

    // Read all histograms from the file and check if they exist
    for (int i = 0; i < nJetTypes; ++i) {

        TString tJetType = jetTypes[i];
        for (int j = 0; j < nHistTitles; ++j) {
            TString histName = Form("hReco%s%s", tJetType.Data(), histTitles[j]);
            if (j == 0) {
                hRecoAllJetPtEta[i] = dynamic_cast<TH2D*>( f->Get(histName) );
                if ( !hRecoAllJetPtEta[i] ) {
                    std::cerr << "Histogram " << histName.Data() << " not found in file." << std::endl;
                    return;
                }
            }
            else if (j == 1) {
                hRecoMatchedJetPtEta[i] = dynamic_cast<TH2D*>( f->Get(histName) );
                if ( !hRecoMatchedJetPtEta[i] ) {
                    std::cerr << "Histogram " << histName.Data() << " not found in file." << std::endl;
                    return;
                }
            }
            else if (j == 2) {
                hRecoUnmatchedJetPtEta[i] = dynamic_cast<TH2D*>( f->Get(histName) );
                if ( !hRecoUnmatchedJetPtEta[i] ) {
                    std::cerr << "Histogram " << histName.Data() << " not found in file." << std::endl;
                    return;
                }
            }
        } // for (int j = 0; j < nHistTitles; ++j)

        // Plot matched jet fraction
        hRecoMatchedFraction[i] = dynamic_cast<TH2D*>( hRecoMatchedJetPtEta[i]->Clone( Form("hReco%sMatchedFraction", tJetType.Data()) ) );
        hRecoMatchedFraction[i]->Divide( hRecoAllJetPtEta[i] );
        cFakes2D->Clear();
        hRecoMatchedFraction[i]->Draw("colz");
        hRecoMatchedFraction[i]->GetXaxis()->SetRangeUser(-3.6, 3.6);
        hRecoMatchedFraction[i]->GetYaxis()->SetRangeUser(5.0, 255.0);
        hRecoMatchedFraction[i]->GetZaxis()->SetRangeUser(0.0, 1.0);
        t.DrawLatexNDC(0.3, yTextPosition, Form("%s Matched Jets", tJetType.Data()));
        plotCMSHeader(collSystem, energy);
        cFakes2D->SaveAs( Form("%s/%s%d_%s_%s_matched2D.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                int(energy * 1000), direction.Data(), tJetType.Data() ) );

        
        // Plot unmatched jet fraction
        hRecoUnmatchedFraction[i] = dynamic_cast<TH2D*>( hRecoUnmatchedJetPtEta[i]->Clone( Form("hReco%sUnmatchedFraction", tJetType.Data()) ) );
        hRecoUnmatchedFraction[i]->Divide( hRecoAllJetPtEta[i] );
        cFakes2D->Clear();
        hRecoUnmatchedFraction[i]->Draw("colz");
        hRecoUnmatchedFraction[i]->GetXaxis()->SetRangeUser(-3.6, 3.6);
        hRecoUnmatchedFraction[i]->GetYaxis()->SetRangeUser(5.0, 255.0);
        hRecoUnmatchedFraction[i]->GetZaxis()->SetRangeUser(0.0, 1.0);
        t.DrawLatexNDC(0.2, yTextPosition, Form("%s Unmatched Jets", tJetType.Data()));
        plotCMSHeader(collSystem, energy);
        cFakes2D->SaveAs( Form("%s/%s%d_%s_%s_fakes2D.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                int(energy * 1000), direction.Data(), tJetType.Data() ) );
    } // for (int i = 0; i < nJetTypes; ++i)

    // Projections in eta and pt to look at fakes
    const int nEtaRanges = 10;  // Number of eta ranges (9 ranges, 10 boundaries)
    double etaRanges[nEtaRanges] = {-5.0, -3.0, -2.4, -2.0, -1.0, 1.0, 2.0, 2.4, 3.0, 5.0};
    const int nPtRanges = 8;   // Number of pT ranges (7 ranges, 8 boundaries)
    double ptRanges[nPtRanges] = {5., 25., 35., 45., 65., 95., 150., 1505.};

    // Declare 1D histograms for all, matched and unmatched jets for eta ranges
    TH1D *hRecoAllEta[nJetTypes][nEtaRanges - 1];
    TH1D *hRecoMatchedEta[nJetTypes][nEtaRanges - 1];
    TH1D *hRecoUnmatchedEta[nJetTypes][nEtaRanges - 1];

    // Declare 1D histograms for all, matched and unmatched jets for pT ranges
    TH1D *hRecoAllPt[nJetTypes][nPtRanges - 1];
    TH1D *hRecoMatchedPt[nJetTypes][nPtRanges - 1];
    TH1D *hRecoUnmatchedPt[nJetTypes][nPtRanges - 1];

    TLegend *leg = nullptr;
    TLegend *legMatched = nullptr;
    TLegend *legUnmatched = nullptr;

    TCanvas *cMatched = new TCanvas( "cMatched", "cMatched", 800, 800 );
    setPadStyle();
    TCanvas *cUnmatched = new TCanvas( "cUnmatched", "cUnmatched", 800, 800 );
    setPadStyle();
    TCanvas *cCombined = new TCanvas( "cCombined", "cCombined", 800, 800 );
    setPadStyle();

    // Plot projections of matched and fakes on eta as a function of pT
    for (size_t i = 0; i < nEtaRanges - 1; ++i) {
        double etaMin = etaRanges[i];
        double etaMax = etaRanges[i + 1];

        for (size_t j = 0; j < nJetTypes; ++j) {
            TString tJetType = jetTypes[j];
            tJetType.ToLower();

            // Project matched and unmatched jets in eta for the current range
            hRecoMatchedEta[j][i] = dynamic_cast<TH1D*>(hRecoMatchedJetPtEta[j]->ProjectionY( Form("hMatchedEta_%d_%d", int(j), int(i)), 
                                                                           hRecoMatchedJetPtEta[j]->GetXaxis()->FindBin(etaMin), 
                                                                           hRecoMatchedJetPtEta[j]->GetXaxis()->FindBin(etaMax) ));
            hRecoUnmatchedEta[j][i] = dynamic_cast<TH1D*>(hRecoUnmatchedJetPtEta[j]->ProjectionY( Form("hUnmatchedEta_%d_%d", int(j), int(i)), 
                                                                               hRecoUnmatchedJetPtEta[j]->GetXaxis()->FindBin(etaMin), 
                                                                               hRecoUnmatchedJetPtEta[j]->GetXaxis()->FindBin(etaMax) ));
            hRecoAllEta[j][i] = dynamic_cast<TH1D*>(hRecoAllJetPtEta[j]->ProjectionY( Form("hAllEta_%d_%d", int(j), int(i)), 
                                                                     hRecoAllJetPtEta[j]->GetXaxis()->FindBin(etaMin), 
                                                                     hRecoAllJetPtEta[j]->GetXaxis()->FindBin(etaMax) ));

            // Normalize matched and unmatched jets
            hRecoMatchedEta[j][i]->Divide( hRecoMatchedEta[j][i], hRecoAllEta[j][i], 1.0, 1.0, "B" ); // Normalize matched
            hRecoUnmatchedEta[j][i]->Divide( hRecoUnmatchedEta[j][i], hRecoAllEta[j][i], 1.0, 1.0, "B" ); // Normalize unmatched

            set1DStyle(hRecoMatchedEta[j][i], 0);
            set1DStyle(hRecoUnmatchedEta[j][i], 1);

            //
            // Draw matched and fakes on the same canvas
            //
            cCombined->Clear();
            cCombined->cd();
            hRecoMatchedEta[j][i]->Draw();
            hRecoUnmatchedEta[j][i]->Draw("same");
            hRecoMatchedEta[j][i]->GetXaxis()->SetRangeUser(5., 150.);
            hRecoMatchedEta[j][i]->GetYaxis()->SetRangeUser(1e-5, 1.1);
            hRecoMatchedEta[j][i]->GetYaxis()->SetTitle("Fraction");
            cCombined->SetLogy();
            plotCMSHeader(collSystem, energy);
            t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%.1f < #eta < %.1f", etaMin, etaMax));
            leg = new TLegend(0.45, 0.5, 0.8, 0.65);
            leg->AddEntry(hRecoMatchedEta[j][i], Form("%s matched", tJetType.Data()), "p");
            leg->AddEntry(hRecoUnmatchedEta[j][i], Form("%s fakes", tJetType.Data()), "p");
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.04);
            leg->Draw();
            // Save the canvas with the date and collision system information
            cCombined->SaveAs( Form("%s/%s%d_%s_combined_%s_eta_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                    int(energy * 1000), direction.Data(), tJetType.Data(), 
                                    ((etaMin<0) ? Form("m%.1f", -etaMin) : Form("%.1f", etaMin)),
                                    ((etaMax<0) ? Form("m%.1f", -etaMax) : Form("%.1f", etaMax)) ) );

            //
            // Draw matched jets on a separate canvas
            //
            set1DStyle(hRecoMatchedEta[j][i], j);
            cMatched->cd();
            if (j == 0) {
                hRecoMatchedEta[j][i]->Draw();
                hRecoMatchedEta[j][i]->GetXaxis()->SetRangeUser(5., 150.);
                hRecoMatchedEta[j][i]->GetYaxis()->SetRangeUser(0.3, 1.1);
                hRecoMatchedEta[j][i]->GetYaxis()->SetTitle("Fraction");

                plotCMSHeader(collSystem, energy);
                t.DrawLatexNDC(xTextPosition, 0.85, Form("%.1f < #eta < %.1f", etaMin, etaMax));

                legMatched = new TLegend(0.55, 0.4, 0.84, 0.6);
                legMatched->SetFillColor(0);
                legMatched->SetBorderSize(0);
                legMatched->SetTextSize(0.04);
                legMatched->SetHeader("Matched Jets");
            }
            else {
                hRecoMatchedEta[j][i]->Draw("same");
            }

            legMatched->AddEntry(hRecoMatchedEta[j][i], Form("%s", tJetType.Data()), "p");
            if (j == (nJetTypes - 1)) {
                legMatched->Draw();
            }

            //
            // Draw unmatched jets on a separate canvas
            //
            set1DStyle(hRecoUnmatchedEta[j][i], j);
            cUnmatched->cd();
            if (j == 0) {
                hRecoUnmatchedEta[j][i]->Draw();
                hRecoUnmatchedEta[j][i]->GetXaxis()->SetRangeUser(5., 150.);
                hRecoUnmatchedEta[j][i]->GetYaxis()->SetRangeUser(1e-5, 1.1);
                hRecoUnmatchedEta[j][i]->GetYaxis()->SetTitle("Fraction");
                plotCMSHeader(collSystem, energy);
                t.DrawLatexNDC(xTextPosition, 0.85, Form("%.1f < #eta < %.1f", etaMin, etaMax));
                legUnmatched = new TLegend(0.55, 0.55, 0.84, 0.75);
                legUnmatched->SetFillColor(0);
                legUnmatched->SetBorderSize(0);
                legUnmatched->SetTextSize(0.04);
                legUnmatched->SetHeader("Fake Jets");
                cUnmatched->SetLogy();
            }
            else {
                hRecoUnmatchedEta[j][i]->Draw("same");
            }

            legUnmatched->AddEntry(hRecoUnmatchedEta[j][i], Form("%s", tJetType.Data()), "p");
            if (j == (nJetTypes - 1)) {
                legUnmatched->Draw();
            }

        } // for (size_t j = 0; j < nJetTypes; ++j)

        cMatched->SaveAs( Form("%s/%s%d_%s_matched_eta_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                int(energy * 1000), direction.Data(), 
                                ((etaMin<0) ? Form("m%.1f", -etaMin) : Form("%.1f", etaMin)),
                                ((etaMax<0) ? Form("m%.1f", -etaMax) : Form("%.1f", etaMax)) ) );
        cMatched->Clear();

        cUnmatched->SaveAs( Form("%s/%s%d_%s_fakes_eta_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                  int(energy * 1000), direction.Data(), 
                                  ((etaMin<0) ? Form("m%.1f", -etaMin) : Form("%.1f", etaMin)),
                                  ((etaMax<0) ? Form("m%.1f", -etaMax) : Form("%.1f", etaMax)) ) );
        cUnmatched->Clear();
    } // for (size_t i = 0; i < nEtaRanges - 1; ++i)

    // Plot projections of matched and fakes on pT as a function of eta
    for (size_t i = 0; i < nPtRanges - 1; ++i) {
        double ptMin = ptRanges[i];
        double ptMax = ptRanges[i + 1];

        for (size_t j = 0; j < nJetTypes; ++j) {
            TString tJetType = jetTypes[j];
            tJetType.ToLower();

            // Project matched and unmatched jets in pT for the current range
            hRecoMatchedPt[j][i] = dynamic_cast<TH1D*>(hRecoMatchedJetPtEta[j]->ProjectionX( Form("hMatchedPt_%d_%d", int(j), int(i)), 
                                                                           hRecoMatchedJetPtEta[j]->GetYaxis()->FindBin(ptMin), 
                                                                           hRecoMatchedJetPtEta[j]->GetYaxis()->FindBin(ptMax) ));
            hRecoUnmatchedPt[j][i] = dynamic_cast<TH1D*>(hRecoUnmatchedJetPtEta[j]->ProjectionX( Form("hUnmatchedPt_%d_%d", int(j), int(i)), 
                                                                               hRecoUnmatchedJetPtEta[j]->GetYaxis()->FindBin(ptMin), 
                                                                               hRecoUnmatchedJetPtEta[j]->GetYaxis()->FindBin(ptMax) ));
            hRecoAllPt[j][i] = dynamic_cast<TH1D*>(hRecoAllJetPtEta[j]->ProjectionX( Form("hAllPt_%d_%d", int(j), int(i)), 
                                                                     hRecoAllJetPtEta[j]->GetYaxis()->FindBin(ptMin), 
                                                                     hRecoAllJetPtEta[j]->GetYaxis()->FindBin(ptMax) ));

            // Normalize matched and unmatched jets
            hRecoMatchedPt[j][i]->Divide( hRecoMatchedPt[j][i], hRecoAllPt[j][i], 1.0, 1.0, "B" ); // Normalize matched
            hRecoUnmatchedPt[j][i]->Divide( hRecoUnmatchedPt[j][i], hRecoAllPt[j][i], 1.0, 1.0, "B" ); // Normalize unmatched

            set1DStyle(hRecoMatchedPt[j][i], 0);
            set1DStyle(hRecoUnmatchedPt[j][i], 1);

            //
            // Draw matched and fakes on the same canvas
            //
            cCombined->Clear();
            cCombined->cd();
            hRecoMatchedPt[j][i]->Draw();
            hRecoUnmatchedPt[j][i]->Draw("same");
            hRecoMatchedPt[j][i]->GetXaxis()->SetRangeUser(-3.6, 3.6);
            hRecoMatchedPt[j][i]->GetYaxis()->SetRangeUser(1e-5, 1.1);
            hRecoMatchedPt[j][i]->GetYaxis()->SetTitle("Fraction");
            cCombined->SetLogy();
            plotCMSHeader(collSystem, energy);
            t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%.1f < p_{T} < %.1f", ptMin, ptMax));
            leg = new TLegend(0.45, 0.5, 0.8, 0.65);
            leg->AddEntry(hRecoMatchedPt[j][i], Form("%s matched", tJetType.Data()), "p");
            leg->AddEntry(hRecoUnmatchedPt[j][i], Form("%s fakes", tJetType.Data()), "p");
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.04);
            leg->Draw();
            // Save the canvas with the date and collision system information  
            cCombined->SaveAs( Form("%s/%s%d_%s_combined_%s_pt_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                    int(energy * 1000), direction.Data(), tJetType.Data(), 
                                    ((ptMin<0) ? Form("m%.1f", -ptMin) : Form("%.1f", ptMin)),
                                    ((ptMax<0) ? Form("m%.1f", -ptMax) : Form("%.1f", ptMax)) ) );
            //
            // Draw matched jets on a separate canvas
            //
            set1DStyle(hRecoMatchedPt[j][i], j);
            cMatched->cd();
            if (j == 0) {
                hRecoMatchedPt[j][i]->Draw();
                hRecoMatchedPt[j][i]->GetXaxis()->SetRangeUser(-3.6, 3.6);
                hRecoMatchedPt[j][i]->GetYaxis()->SetRangeUser(0.3, 1.1);
                hRecoMatchedPt[j][i]->GetYaxis()->SetTitle("Fraction");

                plotCMSHeader(collSystem, energy);
                t.DrawLatexNDC(xTextPosition, 0.85, Form("%.1f < p_{T} < %.1f", ptMin, ptMax));

                legMatched = new TLegend(0.55, 0.4, 0.84, 0.6);
                legMatched->SetFillColor(0);
                legMatched->SetBorderSize(0);
                legMatched->SetTextSize(0.04);
                legMatched->SetHeader("Matched Jets");
            }
            else {
                hRecoMatchedPt[j][i]->Draw("same");
            }
            legMatched->AddEntry(hRecoMatchedPt[j][i], Form("%s", tJetType.Data()), "p");
            if (j == (nJetTypes - 1)) {
                legMatched->Draw();
            }
            //
            //
            // Draw unmatched jets on a separate canvas
            //
            set1DStyle(hRecoUnmatchedPt[j][i], j);
            cUnmatched->cd();
            if (j == 0) {
                hRecoUnmatchedPt[j][i]->Draw();
                hRecoUnmatchedPt[j][i]->GetXaxis()->SetRangeUser(-3.6, 3.6);
                hRecoUnmatchedPt[j][i]->GetYaxis()->SetRangeUser(1e-5, 1.1);
                hRecoUnmatchedPt[j][i]->GetYaxis()->SetTitle("Fraction");
                plotCMSHeader(collSystem, energy);
                t.DrawLatexNDC(xTextPosition, 0.85, Form("%.1f < p_{T} < %.1f", ptMin, ptMax));
                legUnmatched = new TLegend(0.55, 0.55, 0.84, 0.75);
                legUnmatched->SetFillColor(0);
                legUnmatched->SetBorderSize(0);
                legUnmatched->SetTextSize(0.04);
                legUnmatched->SetHeader("Fake Jets");
                cUnmatched->SetLogy();
            }
            else {
                hRecoUnmatchedPt[j][i]->Draw("same");
            }
            legUnmatched->AddEntry(hRecoUnmatchedPt[j][i], Form("%s", tJetType.Data()), "p");
            if (j == (nJetTypes - 1)) {
                legUnmatched->Draw();
            }
        } // for (size_t j = 0; j < nJetTypes; ++j)
        cMatched->SaveAs( Form("%s/%s%d_%s_matched_pt_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                int(energy * 1000), direction.Data(), 
                                ((ptMin<0) ? Form("m%.1f", -ptMin) : Form("%.1f", ptMin)),
                                ((ptMax<0) ? Form("m%.1f", -ptMax) : Form("%.1f", ptMax)) ) );
        cMatched->Clear();
        cUnmatched->SaveAs( Form("%s/%s%d_%s_fakes_pt_%s_%s.pdf", date.Data(), (collSystem == 0 ? "pp" : (collSystem == 1 ? "pPb" : "PbPb")),
                                  int(energy * 1000), direction.Data(), 
                                  ((ptMin<0) ? Form("m%.1f", -ptMin) : Form("%.1f", ptMin)),
                                  ((ptMax<0) ? Form("m%.1f", -ptMax) : Form("%.1f", ptMax)) ) );
        cUnmatched->Clear();
    } // for (size_t i = 0; i < nPtRanges - 1; ++i)
    // Clean up
    for (int i = 0; i < nJetTypes; ++i) {
        for (int j = 0; j < nEtaRanges - 1; ++j) {
            if (hRecoAllEta[i][j]) delete hRecoAllEta[i][j];
            if (hRecoMatchedEta[i][j]) delete hRecoMatchedEta[i][j];
            if (hRecoUnmatchedEta[i][j]) delete hRecoUnmatchedEta[i][j];
        }
        for (int j = 0; j < nPtRanges - 1; ++j) {
            if (hRecoAllPt[i][j]) delete hRecoAllPt[i][j];
            if (hRecoMatchedPt[i][j]) delete hRecoMatchedPt[i][j];
            if (hRecoUnmatchedPt[i][j]) delete hRecoUnmatchedPt[i][j];
        }
    }
    for (int i = 0; i < nJetTypes; ++i) {
        if (hRecoAllJetPtEta[i]) delete hRecoAllJetPtEta[i];
        if (hRecoMatchedJetPtEta[i]) delete hRecoMatchedJetPtEta[i];
        if (hRecoUnmatchedJetPtEta[i]) delete hRecoUnmatchedJetPtEta[i];
    }
    if (cCombined) delete cCombined;
    if (cMatched) delete cMatched;
    if (cUnmatched) delete cUnmatched;
    if (leg) delete leg;
    if (legMatched) delete legMatched;
    if (legUnmatched) delete legUnmatched;
    if (cFakes2D) delete cFakes2D;
}

//________________
void plotEfficiencyAndFakes() {

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

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    int collisionSystem = 1;         // 0 - pp, 1 - pPb, 2 - pPb5020, 3 - pPb8160
    double collisionEnergy = 8.16;   // 8.16 TeV
    int direction = 1;               // 0-p-going, 1-Pb-going, 2 - combined
    TString directionStr = (direction == 0) ? "pgoing" : ((direction == 1) ? "Pbgoing" : "");
    int dataTrigger = 0;               // 0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100
    TString dataStr = (dataTrigger == 0) ? "MB" : ((dataTrigger == 1) ? "Jet60" : ((dataTrigger == 2) ? "Jet80" : ((dataTrigger == 3) ? "Jet100" : "unknownData")));
    TString dataDirectionStr = (direction == 0) ? "Pbgoing" : ((direction == 1) ? "pgoing" : "");
    int jetType = 2; // 0 - inclusive, 1 - lead, 2 - sublead
    int matchType = 0; // 0 - inclusive, 1 - matched, 2 - unmatched

    // Monte Carlo file
    TFile *f = nullptr;
    if ( direction < 2 ) {
        f = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_eta20.root", uname.Data(), directionStr.Data(), directionStr.Data()) );
        if ( !f ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_eta20.root", uname.Data(), directionStr.Data(), directionStr.Data()) << std::endl;
            return;
        }
    }
    else {
        f = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_eta20.root", uname.Data()) );
        if ( !f ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_eta20.root", uname.Data()) << std::endl;
            return;
        }
    }

    // Plot efficiency
    plotEfficiency(f, collisionSystem, collisionEnergy, jetType, date);

    // Plot fakes
    // plotFakes(f, collisionSystem, collisionEnergy, jetType, date);
}
