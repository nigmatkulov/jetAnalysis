// ROOT headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TSystemDirectory.h"
#include <sys/stat.h>

// C++ headers
#include <iostream>

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
void setStyle() {
    auto myStyle = new TStyle("myStyle","Gregs style");

    myStyle->SetCanvasBorderMode(0);
    myStyle->SetPadBorderMode(0);
    myStyle->SetPadColor(0);
    myStyle->SetCanvasColor(0);
    myStyle->SetTitleColor(0);
    myStyle->SetStatColor(0);
    myStyle->SetOptStat(0);
    myStyle->SetFrameLineWidth(2);
    myStyle->SetTitleSize(0.08, "XY");
    myStyle->SetLabelSize(0.08, "XY");
    myStyle->SetLabelOffset(0.007, "XY");
    myStyle->SetLabelFont(72, "XY");
    myStyle->cd();
}

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.1;
    Int_t lineWidth = 2;
    Int_t color = 2;
    if (type == 0) {
        color = 2;
        markerStyle = 20;
    }
    else if (type == 1) {
        color = 4;
        markerStyle = 21;
    }
    else if (type == 2) {
        color = 1;
        markerStyle = 22;
    }
    else if (type == 3) {
        color = 2;
        markerStyle = 24;
    }
    else if (type == 4) {
        color = 4;
        markerStyle = 25;
    }
    else if (type == 5) {
        color = 1;
        markerStyle = 26;
    }
    else {
        color = 9;
        markerStyle = 30;
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
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void set2DStyle(TH2* h, Bool_t doRenorm = kFALSE) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./ h->Integral() );
    }
    
}

//________________
void rescaleEta(TH1* h, Bool_t doRenorm = kTRUE) {
    for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++) {
        Double_t val = h->GetBinContent( iBin );
        Double_t valErr = h->GetBinError( iBin );
        Double_t binWidth = h->GetBinWidth( iBin );
        h->SetBinContent( iBin, val / binWidth );
        h->SetBinError( iBin, valErr / binWidth );
    }
    if ( doRenorm ) {
        h->Scale( 1. / h->Integral() );
    }
}

//________________
void rescaleEta(TH2* h) {
    for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++) {
        for (Int_t jBin=1; jBin<=h->GetNbinsY(); jBin++) {
            Double_t val = h->GetBinContent( iBin, jBin );
            Double_t valErr = h->GetBinError( iBin, jBin );
            Double_t binWidthX = h->GetXaxis()->GetBinWidth( iBin );
            Double_t binWidthY = h->GetYaxis()->GetBinWidth( jBin );
            h->SetBinContent( iBin, jBin, val / (binWidthX * binWidthY) );
            h->SetBinError( iBin, jBin, valErr / (binWidthX * binWidthY) );
        } // for (Int_t jBin=1; jBin<=h->GetNbinsY(); jBin++)
    } // for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++)
    h->Scale( 1. / h->Integral() );
}

//________________
void plotPtHat(TFile *inFile, TString date, Int_t jetBranch = 0) {

    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    TString branchName;
    if ( jetBranch == 0 ) {
        branchName = "akCs4";
    }
    else {
        branchName = "ak4";
    }

    //
    // Plot ptHat distribution
    //
    TH1D* hPtHat = (TH1D*)inFile->Get("hPtHatWeighted");
    hPtHat->SetName("hPtHat");

    TCanvas *cPtHat = new TCanvas("cPtHat","cPtHat", 1200, 800);
    setPadStyle();
    set1DStyle(hPtHat, 2);
    hPtHat->Draw();
    gPad->SetLogy(1);
    cPtHat->SaveAs( Form("%s/pPb8160_%s_ptHat_%s.pdf", 
                         date.Data(), direction.Data(), branchName.Data()));
}

//________________
void plotEfficiency(TFile *inFile, TString date, Int_t jetBranch = 0) {

    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    TString branchName;
    if ( jetBranch == 0 ) {
        branchName = "akCs4";
    }
    else {
        branchName = "ak4";
    }

    // Rebinning
    Int_t rebinX{1}, rebinY{1};

    // Plotting options
    Bool_t plotKineCut{kTRUE};
    Bool_t plotJetIdCut{kTRUE};

    // Make latex
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    //
    // Efficiency over the acceptance
    //

    // Retrieve gen jet acceptance
    TH2D *hGenPtVsEta = (TH2D*)inFile->Get("hGenInclusiveJetPtEta");
    hGenPtVsEta->SetName("hGenPtVsEta");

    // Retrieve ref jet acceptance (for reco passed kinematic cuts only)
    TH2D *hRefPtVsEtaKineCut = (TH2D*)inFile->Get("hRecoInclusiveJetRefPtVsEtaKineCut");
    hRefPtVsEtaKineCut->SetName("hRefPtVsEtaKineCut");

    // Retrieve ref jet acceptance (for reco passed kinematic+trkMax cuts)
    TH2D *hRefPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveJetRefPtVsEtaTrkMaxCut");
    hRefPtVsEtaTrkMaxCut->SetName("hRefPtVsEtaTrkMaxCut");

    // Retrieve ref jet acceptance (for reco passed kinematic+jetId cuts)
    TH2D *hRefPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveJetRefPtVsEtaJetIdCut");
    hRefPtVsEtaJetIdCut->SetName("hRefPtVsEtaJetIdCut");

    // Perform rebinning
    hGenPtVsEta->RebinX( rebinX );
    hGenPtVsEta->RebinY( rebinY );
    hRefPtVsEtaKineCut->RebinX( rebinX );
    hRefPtVsEtaKineCut->RebinY( rebinY );
    hRefPtVsEtaTrkMaxCut->RebinX( rebinX );
    hRefPtVsEtaTrkMaxCut->RebinY( rebinY );
    hRefPtVsEtaJetIdCut->RebinX( rebinX );
    hRefPtVsEtaJetIdCut->RebinY( rebinY );

    // Apply 2D style
    set2DStyle(hGenPtVsEta);
    set2DStyle(hRefPtVsEtaKineCut);
    set2DStyle(hRefPtVsEtaTrkMaxCut);
    set2DStyle(hRefPtVsEtaJetIdCut);

    // Calculate efficiencies
    hRefPtVsEtaKineCut->Divide(hRefPtVsEtaKineCut, hGenPtVsEta, 1., 1., "b");
    hRefPtVsEtaTrkMaxCut->Divide(hRefPtVsEtaTrkMaxCut, hGenPtVsEta, 1., 1., "b");
    hRefPtVsEtaJetIdCut->Divide(hRefPtVsEtaJetIdCut, hGenPtVsEta, 1., 1., "b");

    // Create canvas to store 1D distributions
    TCanvas *c = new TCanvas("c", "c", 1200, 800);

    // Plot and save 2D efficiencies
    if ( plotKineCut ) {
        c->cd();
        setPadStyle();
        hRefPtVsEtaKineCut->Draw("colz");
        c->SaveAs( Form("%s/pPb8160_%s_2D_efficiency_kineCut_%s.pdf", 
                        date.Data(), direction.Data(), branchName.Data()) );
    }

    c->cd();
    setPadStyle();
    hRefPtVsEtaTrkMaxCut->Draw("colz");
    c->SaveAs( Form("%s/pPb8160_%s_2D_efficiency_trkMaxCut_%s.pdf", 
                    date.Data(), direction.Data(), branchName.Data()) );

    if ( plotJetIdCut ) {
        c->cd();
        setPadStyle();
        hRefPtVsEtaJetIdCut->Draw("colz");
        c->SaveAs( Form("%s/pPb8160_%s_2D_efficiency_jetIdCut_%s.pdf", 
                        date.Data(), direction.Data(), branchName.Data()) );
    }


    // Retrieve 1D binning
    Int_t fPtBins = hRefPtVsEtaKineCut->GetNbinsY();
    Double_t fPtRange[2] {hRefPtVsEtaKineCut->GetYaxis()->GetBinLowEdge(1), 
                          hRefPtVsEtaKineCut->GetYaxis()->GetBinUpEdge(fPtBins) };
    Int_t fEtaBins = hRefPtVsEtaKineCut->GetNbinsX();
    Double_t fEtaRange[2] { hRefPtVsEtaKineCut->GetXaxis()->GetBinLowEdge(1),
                            hRefPtVsEtaKineCut->GetXaxis()->GetBinUpEdge(fEtaBins)};
    Double_t ptStep = (fPtRange[1]-fPtRange[0]) / fPtBins;
    Double_t etaStep = (fEtaRange[1]-fEtaRange[0]) / fEtaBins;

    // Prepare histograms with 1D efficiency projections and its ratios 
    // to trkMaxCut (because trkMaxCut is default)
    TH1D *hEtaEfficiencyKineCut[fPtBins];
    TH1D *hPtEfficiencyKineCut[fEtaBins];
    TH1D *hEtaEfficiencyTrkMaxCut[fPtBins];
    TH1D *hPtEfficiencyTrkMaxCut[fEtaBins];
    TH1D *hEtaEfficiencyJetIdCut[fPtBins];
    TH1D *hPtEfficiencyJetIdCut[fEtaBins];

    TH1D *hEtaEfficiencyRatioKineCut[fPtBins];
    TH1D *hPtEfficiencyRatioKineCut[fEtaBins];
    TH1D *hEtaEfficiencyRatioJetIdCut[fPtBins];
    TH1D *hPtEfficiencyRatioJetIdCut[fEtaBins];

    // Create legends
    TLegend *hEtaLegend[fPtBins];
    TLegend *hPtLegend[fEtaBins];

    TCanvas *cEtaEfficiency = new TCanvas("cEtaEfficiency", "cEtaEfficiency", 1600, 800);
    cEtaEfficiency->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtEfficiency = new TCanvas("cPtEfficiency", "cPtEfficiency", 1600, 800);
    cPtEfficiency->Divide(5, ( (fEtaBins % 5) == 0 ) ? (fEtaBins / 5) : (fEtaBins / 5 + 1) );

    //
    // Make projections on eta
    //
    for (Int_t i{1}; i<=fPtBins; i++) {

        // TrkMax
        hEtaEfficiencyTrkMaxCut[i] = (TH1D*)hRefPtVsEtaTrkMaxCut->ProjectionX(Form("hEtaEfficiencyTrkMaxCut_%d", i), i, i);
        hEtaEfficiencyTrkMaxCut[i]->SetNameTitle(Form("hEtaEfficiencyTrkMaxCut_%d", i), ";#eta;Efficiency");
        set1DStyle(hEtaEfficiencyTrkMaxCut[i], 2);
        hEtaEfficiencyTrkMaxCut[i]->SetMarkerSize(0.7);

        // Kine
        if ( plotKineCut ) {
            hEtaEfficiencyKineCut[i] = (TH1D*)hRefPtVsEtaKineCut->ProjectionX(Form("hEtaEfficiencyKineCut_%d", i), i, i);
            hEtaEfficiencyKineCut[i]->SetNameTitle(Form("hEtaEfficiencyKineCut_%d", i), ";#eta;Efficiency");
            set1DStyle(hEtaEfficiencyKineCut[i], 0);
            hEtaEfficiencyKineCut[i]->SetMarkerSize(0.7);

            hEtaEfficiencyRatioKineCut[i] = (TH1D*)hEtaEfficiencyKineCut[i]->Clone(Form("hEtaEfficiencyRatioKineCut_%d",i));
            hEtaEfficiencyRatioKineCut[i]->Divide( hEtaEfficiencyRatioKineCut[i], hEtaEfficiencyTrkMaxCut[i], 1., 1., "b");
            hEtaEfficiencyRatioKineCut[i]->GetYaxis()->SetTitle("Ratio to trkMax");
        }

        // JetId
        if ( plotJetIdCut ) {
            hEtaEfficiencyJetIdCut[i] = (TH1D*)hRefPtVsEtaJetIdCut->ProjectionX(Form("hEtaEfficiencyJetIdCut_%d", i), i, i);
            hEtaEfficiencyJetIdCut[i]->SetNameTitle(Form("hEtaEfficiencyJetIdCut_%d", i), ";#eta;Efficiency");
            set1DStyle(hEtaEfficiencyJetIdCut[i], 1);
            hEtaEfficiencyJetIdCut[i]->SetMarkerSize(0.7);

            hEtaEfficiencyRatioJetIdCut[i] = (TH1D*)hEtaEfficiencyJetIdCut[i]->Clone(Form("hEtaEfficiencyRatioJetIdCut_%d",i));
            hEtaEfficiencyRatioJetIdCut[i]->Divide( hEtaEfficiencyRatioJetIdCut[i], hEtaEfficiencyTrkMaxCut[i], 1., 1., "b");
            hEtaEfficiencyRatioJetIdCut[i]->GetYaxis()->SetTitle("Ratio to trkMax");
        }

        // Plot individual distributions
        c->cd();
        setPadStyle();
        hEtaEfficiencyTrkMaxCut[i]->Draw();
        hEtaEfficiencyTrkMaxCut[i]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotKineCut ) {
            hEtaEfficiencyKineCut[i]->Draw("same");
        }
        if ( plotJetIdCut ) {
            hEtaEfficiencyJetIdCut[i]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );

        hEtaLegend[i] = new TLegend(0.35, 0.35, 0.65, 0.55);
        hEtaLegend[i]->SetTextSize(0.05);
        hEtaLegend[i]->SetLineWidth(0);
        hEtaLegend[i]->AddEntry(hEtaEfficiencyTrkMaxCut[i], Form("TrkMax"), "p");
        if ( plotKineCut ) {
            hEtaLegend[i]->AddEntry(hEtaEfficiencyKineCut[i], Form("KineOnly"), "p");
        }
        if ( plotJetIdCut ) {
            hEtaLegend[i]->AddEntry(hEtaEfficiencyJetIdCut[i], Form("JetId"), "p");
        }
        hEtaLegend[i]->Draw();
        c->SaveAs( Form("%s/pPb8160_%s_eta_efficiency_pt_%d_%d_%s.pdf", 
                        date.Data(), direction.Data(), 
                        (Int_t)(fPtRange[0] + (i-1) * ptStep),
                        (Int_t)(fPtRange[0] + i * ptStep),
                        branchName.Data()) );

        // Plot ratios
        if ( plotKineCut || plotJetIdCut ) {
            c->cd();
            setPadStyle();
            if ( plotKineCut ) {
                hEtaEfficiencyRatioKineCut[i]->Draw();
                hEtaEfficiencyRatioKineCut[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
            }
            if ( plotJetIdCut ) {
                if ( plotKineCut ) {
                    hEtaEfficiencyRatioJetIdCut[i]->Draw("same");
                }
                else {
                    hEtaEfficiencyRatioJetIdCut[i]->Draw();
                    hEtaEfficiencyRatioJetIdCut[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
                }
            }
            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                        fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );
            hEtaLegend[i]->Draw();
            c->SaveAs( Form("%s/pPb8160_%s_eta_efficiency_rat2trkMax_pt_%d_%d_%s.pdf", 
                            date.Data(), direction.Data(), 
                            (Int_t)(fPtRange[0] + (i-1) * ptStep),
                            (Int_t)(fPtRange[0] + i * ptStep),
                            branchName.Data()) );
        } // if ( plotKineCut || plotJetIdCut )

        // Plot all in one canvas
        cEtaEfficiency->cd(i);
        setPadStyle();
        hEtaEfficiencyTrkMaxCut[i]->Draw();
        hEtaEfficiencyTrkMaxCut[i]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotKineCut ) {
            hEtaEfficiencyKineCut[i]->Draw("same");
        }
        if ( plotJetIdCut ) {
            hEtaEfficiencyJetIdCut[i]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );
        hEtaLegend[i]->Draw();
    } // for (Int_t i{1}; i<=fPtBins; i++)
    cEtaEfficiency->SaveAs( Form("%s/pPb8160_%s_eta_efficiency_projections_%s.pdf", 
                                 date.Data(), direction.Data(), branchName.Data()) );

    //
    // Projections on pT axis
    //
    for (Int_t i{1}; i<=fEtaBins; i++) {

        // TrkMax
        hPtEfficiencyTrkMaxCut[i] = (TH1D*)hRefPtVsEtaTrkMaxCut->ProjectionY(Form("hPtEfficiencyTrkMaxCut_%d", i), i, i);
        hPtEfficiencyTrkMaxCut[i]->SetNameTitle(Form("hPtEfficiencyTrkMaxCut_%d", i), ";p_{T} (GeV/c);Efficiency");
        set1DStyle(hPtEfficiencyTrkMaxCut[i], 2);
        hPtEfficiencyTrkMaxCut[i]->SetMarkerSize(0.7);

        // Kine
        if ( plotKineCut ) {
            hPtEfficiencyKineCut[i] = (TH1D*)hRefPtVsEtaKineCut->ProjectionY(Form("hPtEfficiencyKineCut_%d", i), i, i);
            hPtEfficiencyKineCut[i]->SetNameTitle(Form("hPtEfficiencyKineCut_%d", i), ";p_{T} (GeV/c);Efficiency");
            set1DStyle(hPtEfficiencyKineCut[i], 0);
            hPtEfficiencyKineCut[i]->SetMarkerSize(0.7);

            hPtEfficiencyRatioKineCut[i] = (TH1D*)hPtEfficiencyKineCut[i]->Clone(Form("hPtEfficiencyRatioKineCut_%d",i));
            hPtEfficiencyRatioKineCut[i]->Divide( hPtEfficiencyRatioKineCut[i], hPtEfficiencyTrkMaxCut[i], 1., 1., "b");
            hPtEfficiencyRatioKineCut[i]->GetYaxis()->SetTitle("Ratio to trkMax");
        }

        // JetId
        if ( plotJetIdCut ) {
            hPtEfficiencyJetIdCut[i] = (TH1D*)hRefPtVsEtaJetIdCut->ProjectionY(Form("hPtEfficiencyJetIdCut_%d", i), i, i);
            hPtEfficiencyJetIdCut[i]->SetNameTitle(Form("hPtEfficiencyJetIdCut_%d", i), ";p_{T} (GeV/c);Efficiency");
            set1DStyle(hPtEfficiencyJetIdCut[i], 1);
            hPtEfficiencyJetIdCut[i]->SetMarkerSize(0.7);

            hPtEfficiencyRatioJetIdCut[i] = (TH1D*)hPtEfficiencyJetIdCut[i]->Clone(Form("hPtEfficiencyRatioJetIdCut_%d",i));
            hPtEfficiencyRatioJetIdCut[i]->Divide( hPtEfficiencyRatioJetIdCut[i], hPtEfficiencyTrkMaxCut[i], 1., 1., "b");
            hPtEfficiencyRatioJetIdCut[i]->GetYaxis()->SetTitle("Ratio to trkMax");
        }

        // Plot individual distributions
        c->cd();
        setPadStyle();
        hPtEfficiencyTrkMaxCut[i]->Draw();
        hPtEfficiencyTrkMaxCut[i]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotKineCut ) {
            hPtEfficiencyKineCut[i]->Draw("same");
        }
        if ( plotJetIdCut ) {
            hPtEfficiencyJetIdCut[i]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                       fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );

        hPtLegend[i] = new TLegend(0.35, 0.35, 0.65, 0.55);
        hPtLegend[i]->SetTextSize(0.05);
        hPtLegend[i]->SetLineWidth(0);
        hPtLegend[i]->AddEntry(hPtEfficiencyTrkMaxCut[i], Form("TrkMax"), "p");
        if ( plotKineCut ) {
            hPtLegend[i]->AddEntry(hPtEfficiencyKineCut[i], Form("KineOnly"), "p");
        }
        if ( plotJetIdCut ) {
            hPtLegend[i]->AddEntry(hPtEfficiencyJetIdCut[i], Form("JetId"), "p");
        }
        hPtLegend[i]->Draw();
        c->SaveAs( Form("%s/pPb8160_%s_pt_efficiency_eta_%d_%d_%s.pdf", 
                        date.Data(), direction.Data(), 
                        (Int_t)( (fEtaRange[0] + (i-1) * etaStep) * 10 ),
                        (Int_t)( (fEtaRange[0] + i * etaStep) * 10),
                        branchName.Data()) );

        // Plot ratios
        if ( plotKineCut || plotJetIdCut ) {
            c->cd();
            setPadStyle();
            if ( plotKineCut ) {
                hPtEfficiencyRatioKineCut[i]->Draw();
                hPtEfficiencyRatioKineCut[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
            }
            if ( plotJetIdCut ) {
                if ( plotKineCut ) {
                    hPtEfficiencyRatioJetIdCut[i]->Draw("same");
                }
                else {
                    hPtEfficiencyRatioJetIdCut[i]->Draw();
                    hPtEfficiencyRatioJetIdCut[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
                }
            }
            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                           fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );
            hPtLegend[i]->Draw();
            c->SaveAs( Form("%s/pPb8160_%s_pt_efficiency_rat2trkMax_eta_%d_%d_%s.pdf", 
                            date.Data(), direction.Data(), 
                            (Int_t)( (fEtaRange[0] + (i-1) * etaStep) * 10),
                            (Int_t)( (fEtaRange[0] + i * etaStep) * 10),
                            branchName.Data()) );
        } // if ( plotKineCut || plotJetIdCut )

        // Plot all in one canvas
        cPtEfficiency->cd(i);
        setPadStyle();
        hPtEfficiencyTrkMaxCut[i]->Draw();
        hPtEfficiencyTrkMaxCut[i]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotKineCut ) {
            hPtEfficiencyKineCut[i]->Draw("same");
        }
        if ( plotJetIdCut ) {
            hPtEfficiencyJetIdCut[i]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                       fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );
        hPtLegend[i]->Draw();
    } // for (Int_t i{1}; i<=fEtaBins; i++)
    cPtEfficiency->SaveAs( Form("%s/pPb8160_%s_pt_efficiency_projections_%s.pdf", 
                                date.Data(), direction.Data(), branchName.Data()) );
}

//________________
void plotEtaDijetCorrelation(TFile *inFile, TString date) {

    TH2D *hDijetEtaRefVsReco = (TH2D*)inFile->Get("hRefDijetEtaVsRecoDijetEta");
    hDijetEtaRefVsReco->SetName("hDijetEtaRefVsReco");

    // Plot efficiency
    TCanvas *cDijetEtaRefVsReco = new TCanvas("cDijetEtaRefVsReco","cDijetEtaRefVsReco", 800, 800);
    setPadStyle();
    hDijetEtaRefVsReco->Draw("colz");
    gPad->SetLogz(1);

    cDijetEtaRefVsReco->SaveAs(Form("%s/pPb8160_ref_vs_reco_responce.pdf", date.Data()));
}

//________________
void plotJetIdHistos(TFile *inFile, TString date) {

    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    // Retrieve histograms
    TH1D *hNHF[4];
    TH1D *hNEmF[4];
    TH1D *hNumOfConst[4];
    TH1D *hMUF[4];
    TH1D *hCHF[4];
    TH1D *hChargedMult[4];
    TH1D *hCEmF[4];
    TH1D *hNumOfNeutPart[4];

    TCanvas *c[4];
    TString text;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Loop over 4 eta bins
    for (Int_t i{0}; i<4; i++) {

        if      ( i == 0 ) { text = "-2.4 #leq #eta #leq 2.4"; }
        else if ( i == 1 ) { text = "2.4 < |#eta| #leq 2.7"; }
        else if ( i == 2 ) { text = "2.7 < |#eta| #leq 3.0"; }
        else               { text = "|#eta| > 3.0"; }


        hNHF[i] = (TH1D*)inFile->Get(Form("hNHF_%d",i));
        set1DStyle(hNHF[i], 0, kTRUE);

        hNEmF[i] = (TH1D*)inFile->Get(Form("hNEmF_%d",i));
        set1DStyle(hNEmF[i], 0, kTRUE);

        hNumOfConst[i] = (TH1D*)inFile->Get(Form("hNumOfConst_%d",i));
        set1DStyle(hNumOfConst[i], 0, kTRUE);

        hMUF[i] = (TH1D*)inFile->Get(Form("hMUF_%d",i));
        set1DStyle(hMUF[i], 0, kTRUE);

        hCHF[i] = (TH1D*)inFile->Get(Form("hCHF_%d",i));
        set1DStyle(hCHF[i], 0, kTRUE);

        hChargedMult[i] = (TH1D*)inFile->Get(Form("hChargedMult_%d",i));
        set1DStyle(hChargedMult[i], 0, kTRUE);

        hCEmF[i] = (TH1D*)inFile->Get(Form("hCEmF_%d",i));
        set1DStyle(hCEmF[i], 0, kTRUE);

        hNumOfNeutPart[i] = (TH1D*)inFile->Get(Form("hNumOfNeutPart_%d",i));
        set1DStyle(hNumOfNeutPart[i], 0, kTRUE);

        c[i] = new TCanvas(Form("c%d",i), Form("c%d",i), 1200, 800);
        c[i]->Divide(4, 2);

        c[i]->cd(1);
        setPadStyle();
        hNEmF[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(2);
        setPadStyle();
        hNHF[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(3);
        setPadStyle();
        hNumOfConst[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(4);
        setPadStyle();
        hMUF[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(5);
        setPadStyle();
        hCHF[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(6);
        setPadStyle();
        hChargedMult[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(7);
        setPadStyle();
        hCEmF[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->cd(8);
        setPadStyle();
        hNumOfNeutPart[i]->Draw();
        gPad->SetLogy();
        t.DrawLatexNDC(0.4, 0.93, text.Data() );

        c[i]->SaveAs(Form("%s/pPb8160_%s_jetId_%d.pdf", date.Data(), direction.Data(), i) );
    }
}

//________________
// jetBranch: 0 - akCs4PF, 1 - ak4PF
void plotJESandJER(TFile *inFile, TString date, Int_t jetBranch = 0) {

    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    Int_t rebinX{1}, rebinY{1};

    // Retrieve JES histogram
    THnSparseD *hJESPtEtaPhi = (THnSparseD*)inFile->Get("hJESInclusiveJetPtEtaPhiWeighted");
    hJESPtEtaPhi->SetName("hJESPtEtaPhi");

    TString branchName;
    if ( jetBranch == 0 ) {
        branchName = "akCs4";
    }
    else {
        branchName = "ak4";
    }

    // Get JES vs ref pT
    hJESPtEtaPhi->GetAxis(2)->SetRange(20, 32); // -1.2 < eta < 1.2
    TH2D* hJESvsPt = (TH2D*)hJESPtEtaPhi->Projection(0, 1);
    hJESvsPt->SetName("hJESvsPt");
    hJESvsPt->RebinX( rebinX );
    hJESvsPt->RebinY( rebinY );
    set2DStyle( hJESvsPt );

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Create canvas
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    setPadStyle();
    hJESvsPt->Draw("colz");
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    canv->SaveAs( Form("%s/pPb8160_%s_jes_2d_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );

    // Create histograms for JES and JER
    TH1D *hJESMean;
    TH1D *hJESSigma;

    // Fit 2D with gaussian
    hJESvsPt->FitSlicesY();
    // Obtain JES
    hJESMean = (TH1D*)gDirectory->Get("hJESvsPt_1");
    hJESMean->SetName("hJESMean");
    set1DStyle( hJESMean, 2 );
    // Obtain JER
    hJESSigma = (TH1D*)gDirectory->Get("hJESvsPt_2");
    hJESSigma->SetName("hJESSigma");
    set1DStyle( hJESSigma, 2 );

    hJESMean->Draw();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    hJESMean->GetYaxis()->SetRangeUser(0.95, 1.1);
    hJESMean->GetYaxis()->SetTitle("JES");
    hJESMean->GetXaxis()->SetRangeUser(20, 200);
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    canv->SaveAs( Form("%s/pPb8160_%s_jes_mu_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );

    hJESSigma->Draw();
    hJESSigma->GetXaxis()->SetRangeUser(20, 800);
    hJESSigma->GetYaxis()->SetTitle("JER");
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    canv->SaveAs( Form("%s/pPb8160_%s_jes_sigma_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );

    // Return eta to original area
    hJESPtEtaPhi->GetAxis(2)->SetRange(1, 52);

    // Plot 2D, mu and sigma on a single canvas
    TCanvas *canv2 = new TCanvas("canv2", "canv2", 1200, 400);
    canv2->Divide(3, 1);

    canv2->cd(1);
    setPadStyle();
    hJESvsPt->Draw("colz");
    gPad->SetLogz();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );

    canv2->cd(2);
    setPadStyle();
    hJESMean->Draw();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );

    canv2->cd(3);
    setPadStyle();
    hJESSigma->Draw();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    
    canv2->SaveAs( Form("%s/pPb8160_%s_jes_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );

    // Make JES and JER as a function of eta for few pT bins
    Int_t ptStep{10};
    Int_t ptFirst{20};
    std::vector<Int_t> ptLow {2, 10, 20};
    std::vector<Int_t> ptHi  {6, 14, 27};
    TH2D *hJESvsEta[ ptLow.size() ];
    TH1D *hJESvsEtaMean[ ptLow.size() ];
    TH1D *hJESvsEtaSigma[ ptLow.size() ];
    TCanvas *canv3[3];
    // Loop over pT bins
    for (Int_t i=0; i<ptLow.size(); i++) {

        // Create canvas for the given pT bin
        canv3[i] = new TCanvas(Form("canv3_%d",i), Form("canv3_%d",i), 1200, 400);
        canv3[i]->Divide(3, 1);

        // Set pT range
        hJESPtEtaPhi->GetAxis(1)->SetRange( ptLow.at(i), ptHi.at(i) );
        hJESvsEta[i] = (TH2D*)hJESPtEtaPhi->Projection(0, 2);
        hJESvsEta[i]->SetName( Form("hJESvsEta_%d",i) );
        hJESvsEta[i]->GetXaxis()->SetTitle("#eta");
        hJESvsEta[i]->RebinX(2);
        hJESvsEta[i]->RebinY(2);
        set2DStyle( hJESvsEta[i] );

        // Fit slices
        hJESvsEta[i]->FitSlicesY();
        hJESvsEtaMean[i] = (TH1D*)gDirectory->Get( Form("hJESvsEta_%d_1", i) );
        hJESvsEtaMean[i]->SetName( Form("hJESvsEtaMean_%d", i) );
        hJESvsEtaMean[i]->GetYaxis()->SetTitle("JES");
        set1DStyle(hJESvsEtaMean[i], 2);
        hJESvsEtaSigma[i] = (TH1D*)gDirectory->Get( Form("hJESvsEta_%d_2", i) );
        hJESvsEtaSigma[i]->SetName( Form("hJESvsEtaSigma_%d", i) );
        hJESvsEtaSigma[i]->GetYaxis()->SetTitle("JER");
        set1DStyle(hJESvsEtaSigma[i],2);

        canv->cd();
        setPadStyle();
        hJESvsEta[i]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );
        canv->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_pt_%d_%d_%s.pdf", 
                           date.Data(),  direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                           ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );

        canv->cd();
        setPadStyle();
        hJESvsEtaMean[i]->Draw();
        hJESvsEtaMean[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );
        canv->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_mu_pt_%d_%d_%s.pdf", 
                           date.Data(),  direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                           ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );

        canv->cd();
        setPadStyle();
        hJESvsEtaSigma[i]->Draw();
        hJESvsEtaSigma[i]->GetYaxis()->SetRangeUser(0., 0.2);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );
        canv->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_mu_pt_%d_%d_%s.pdf", 
                           date.Data(), direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                           ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );

        canv3[i]->cd(1);
        setPadStyle();
        hJESvsEta[i]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );

        canv3[i]->cd(2);
        setPadStyle();
        hJESvsEtaMean[i]->Draw();
        hJESvsEtaMean[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );

        canv3[i]->cd(3);
        setPadStyle();
        hJESvsEtaSigma[i]->Draw();
        hJESvsEtaSigma[i]->GetYaxis()->SetRangeUser(0., 0.2);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );

        canv3[i]->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_all_pt_%d_%d_%s.pdf", 
                               date.Data(), direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                               ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );
    }

}

//________________
void plotDijetDistributions(TFile *inFile, TString date) {

    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};
    Int_t refSelType{2};
    Bool_t doRenorm{kFALSE};

    TH1D *hRecoDijetEta = (TH1D*)inFile->Get("hRecoDijetEta");
    if ( !hRecoDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No recoDijetEta found\n";
    }
    TH1D *hRefDijetEta = (TH1D*)inFile->Get("hRefDijetEta");
    if ( !hRefDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No refDijetEta found\n";
    }
    TH1D *hGenDijetEta = (TH1D*)inFile->Get("hGenDijetEta");
    if ( !hGenDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No genDijetEta found\n";
    }
    TH1D *hRefSelDijetEta = (TH1D*)inFile->Get("hRefSelDijetEta");
    if ( !hRefSelDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No refSelDijetEta found\n";
    }

    set1DStyle(hRecoDijetEta,   recoType, doRenorm);
    set1DStyle(hRefDijetEta,    refType, doRenorm);
    set1DStyle(hGenDijetEta,    genType, doRenorm);
    set1DStyle(hRefSelDijetEta, refSelType, doRenorm);

    rescaleEta( hRecoDijetEta );
    rescaleEta( hRefDijetEta );
    rescaleEta( hGenDijetEta );
    rescaleEta( hRefSelDijetEta );

    TH1D *hReco2Gen = new TH1D("hReco2Gen", "hReco2Gen;#eta_{dijet};#frac{reco}{gen}",
                                hRecoDijetEta->GetNbinsX(), 
                                hRecoDijetEta->GetXaxis()->GetBinLowEdge(1),
                                hRecoDijetEta->GetXaxis()->GetBinUpEdge( hRecoDijetEta->GetNbinsX() ) );
    hReco2Gen->Sumw2();
    set1DStyle(hReco2Gen, recoType);
    TH1D *hRef2Gen = new TH1D("hRef2Gen", "hRef2Gen;#eta_{dijet};#frac{ref}{gen}",
                                hRefDijetEta->GetNbinsX(), 
                                hRefDijetEta->GetXaxis()->GetBinLowEdge(1),
                                hRefDijetEta->GetXaxis()->GetBinUpEdge( hRefDijetEta->GetNbinsX() ) );
    hRef2Gen->Sumw2();
    set1DStyle(hRef2Gen, refType);
    TH1D *hRefSel2Gen = new TH1D("hRefSel2Gen", "hRefSel2Gen;#eta_{dijet};#frac{ref sel}{gen}",
                                hRefDijetEta->GetNbinsX(), 
                                hRefDijetEta->GetXaxis()->GetBinLowEdge(1),
                                hRefDijetEta->GetXaxis()->GetBinUpEdge( hRefDijetEta->GetNbinsX() ) );
    hRefSel2Gen->Sumw2();
    set1DStyle(hRefSel2Gen, refSelType);

    hReco2Gen->Divide(hRecoDijetEta, hGenDijetEta, 1., 1., "b");
    hRef2Gen->Divide(hRefDijetEta, hGenDijetEta, 1., 1., "b");
    hRefSel2Gen->Divide(hRefSelDijetEta, hGenDijetEta, 1., 1., "b");

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TCanvas *c1DEta = new TCanvas("c1DEta", "c1DEta", 1000, 1000);
    c1DEta->Divide(1, 2);

    c1DEta->cd(1);
    setPadStyle();
    hGenDijetEta->Draw();
    hRecoDijetEta->Draw("same");
    hRefDijetEta->Draw("same");
    hRefSelDijetEta->Draw("same");
    hGenDijetEta->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    TLegend *leg = new TLegend(0.65, 0.68, 0.8, 0.92);
    leg->SetLineWidth(0);
    leg->AddEntry(hGenDijetEta,Form("Gen"), "p");
    leg->AddEntry(hRefDijetEta,Form("Ref"), "p");
    leg->AddEntry(hRecoDijetEta,Form("Reco"), "p");
    leg->AddEntry(hRefSelDijetEta,Form("Ref Sel"), "p");
    leg->SetTextSize(0.06);
    leg->Draw();

    c1DEta->cd(2);
    setPadStyle();
    hReco2Gen->Draw();
    hRef2Gen->Draw("same");
    hRefSel2Gen->Draw("same");
    TLegend *leg2 = new TLegend(0.65, 0.68, 0.8, 0.92);
    leg2->SetLineWidth(0);
    leg2->AddEntry(hRef2Gen,Form("Ref/Gen"), "p");
    leg2->AddEntry(hReco2Gen,Form("Reco/Gen"), "p");
    leg2->AddEntry(hRefSel2Gen,Form("RefSel/Gen"), "p");
    leg2->SetTextSize(0.06);
    leg2->Draw();
    hReco2Gen->GetYaxis()->SetRangeUser(0.8, 1.2);
    hReco2Gen->GetYaxis()->SetTitle("Ratio to Gen");
    TLine *l = new TLine(hReco2Gen->GetXaxis()->GetBinLowEdge(1), 1., 
                         hReco2Gen->GetXaxis()->GetBinUpEdge(hReco2Gen->GetNbinsX()), 1.);
    l->SetLineColor(kMagenta);
    l->SetLineWidth(3);
    l->SetLineStyle(3);
    l->Draw();

    c1DEta->SaveAs( Form("%s/pPb8160_%s_eta_dijet_comparison.pdf", date.Data(), direction.Data()) );
}

//________________
void compareInclusiveJetPtSpectra(TFile *inFile, TString date) {

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};

    // Reco inclusive jet pT
    TH1D *hReco = (TH1D*)inFile->Get("hRecoInclusiveJetPt");
    // Reco inclusive jet pT
    TH1D *hRef = (TH1D*)inFile->Get("hRefInclusiveJetPt");
    // Reco inclusive jet pT
    TH1D *hGen = (TH1D*)inFile->Get("hGenInclusiveJetPt");

    // Set style
    set1DStyle(hReco, recoType, 0);
    set1DStyle(hRef, refType, 0);
    set1DStyle(hGen, genType, 0);

    // Create histograms to fill with ratios
    TH1D* hReco2GenRatio = new TH1D("hReco2GenRatio","reco / gen ratio;p_{T} (GeV/c);#frac{reco}{gen}",
                                    hGen->GetNbinsX(), 
                                    hGen->GetXaxis()->GetBinLowEdge(1),
                                    hGen->GetXaxis()->GetBinUpEdge( hGen->GetNbinsX() ) );
    hReco2GenRatio->Sumw2();
    set1DStyle(hReco2GenRatio, recoType);
    TH1D* hRef2GenRatio = new TH1D("hRef2GenRatio","ref / gen ratio;p_{T} (GeV/c);#frac{ref}{gen}",
                                    hGen->GetNbinsX(), 
                                    hGen->GetXaxis()->GetBinLowEdge(1),
                                    hGen->GetXaxis()->GetBinUpEdge( hGen->GetNbinsX() ) );
    hRef2GenRatio->Sumw2();
    set1DStyle(hRef2GenRatio, genType);


    // Make ratios
    hReco2GenRatio->Divide(hReco, hGen, 1., 1., "b");
    hRef2GenRatio->Divide(hRef, hGen, 1., 1., "b");

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TCanvas *cPtSpectra = new TCanvas("cPtSpectra", "cPtSpectra", 600, 1200);

    cPtSpectra->Divide(1, 2);

    cPtSpectra->cd(1);
    setPadStyle();
    hReco->Draw();
    hRef->Draw("same");
    hGen->Draw("same");
    hReco->GetXaxis()->SetTitle("Inclusive jet p_{T} (GeV/c)");
    gPad->SetLogy(1);

    TLegend *leg = new TLegend(0.7, 0.68, 0.85, 0.92);
    leg->SetLineWidth(0);
    leg->AddEntry(hReco,Form("Reco"), "p");
    leg->AddEntry(hRef,Form("Ref"), "p");
    leg->AddEntry(hGen,Form("Gen"), "p");
    leg->SetTextSize(0.06);
    leg->Draw();

    cPtSpectra->cd(2);
    setPadStyle();
    hReco2GenRatio->Draw();
    hRef2GenRatio->Draw("same");
    //hGen->Draw("same");
    hReco2GenRatio->GetYaxis()->SetRangeUser(0.87, 1.07);
    hReco2GenRatio->GetYaxis()->SetTitle("Ratio to gen");
    TLegend *leg2 = new TLegend(0.65, 0.75, 0.82, 0.92);
    leg2->SetLineWidth(0);
    leg2->AddEntry(hReco2GenRatio,Form("Reco/Gen"), "p");
    leg2->AddEntry(hRef2GenRatio,Form("Ref/Gen"), "p");
    leg2->SetTextSize(0.06);
    leg2->Draw();
}

//________________
void makeProjectionsFrom2D(TH2D *h2D, TH1D *hProjX[], TH1D *hProjY[], 
                           Int_t nBinsX = 1, Int_t nBinsY = 1, 
                           const Char_t *hNameX = "hProjX", const Char_t *hNameY = "hProjY",
                           Int_t style = 1) {

    std::cout << "nBinsX: " << nBinsX << " nBinsY: " << nBinsY << std::endl;
    // Make projections on X axis
    for (Int_t i{0}; i<nBinsY; i++) {
        hProjX[i] = dynamic_cast<TH1D*>( h2D->ProjectionX( Form("%s_%d", hNameX, i), i+1, i+1) );
        hProjX[i]->SetNameTitle( Form("%s_%d", hNameX, i), ";#eta");
        set1DStyle( hProjX[i], style );
        hProjX[i]->SetMarkerSize(0.7);
    } // for (Int_t i{1}; i<=nBinsX; i++)

    // Make projections on Y axis
    for (Int_t i{0}; i<nBinsX; i++) {
        hProjY[i] = dynamic_cast<TH1D*>( h2D->ProjectionY( Form("%s_%d", hNameY, i), i+1, i+1) );
        hProjY[i]->SetNameTitle( Form("%s_%d", hNameY, i), ";p_{T} (GeV/c)");
        set1DStyle( hProjY[i], style );
        hProjY[i]->SetMarkerSize(0.7);
        std::cout << " hProjY: " << hProjY[i]->GetName() << std::endl;
    } // for (Int_t i{1}; i<=nBinsY; i++)
}

//________________
void make1DRatio(TH1D *hRat, TH1D *hDen, const Char_t *ratioName = "Ratio to TrkMax", Int_t style = 0) {

    if ( !hDen ) {
        std::cout << "Denominator does not exist" << std::endl;
    }

    hRat->Divide( hRat, hDen, 1., 1., "b" );
    hRat->GetYaxis()->SetTitle( ratioName );
    hRat->GetYaxis()->SetRangeUser(0.85, 1.15);
    set1DStyle(hRat, style);
}

//________________
void plotRecoAndFakes(TFile *inFile, TString date, Int_t jetBranch = 0) {

    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    TString branchName;
    if ( jetBranch == 0 ) {
        branchName = "akCs4";
    }
    else {
        branchName = "ak4";
    }

    // Rebinning
    Int_t rebinX{1}, rebinY{1};

    // Plotting options
    Bool_t plotKineCut{kTRUE};
    Bool_t plotJetIdCut{kTRUE};

    // Make latex
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Color style
    Int_t kineStyle{0};
    Int_t trkMaxStyle{2};
    Int_t jetIdStyle{1};

    //
    // Retrieve kine selection histograms
    //
    TH2D *hJetPtVsEtaKineCut = (TH2D*)inFile->Get("hRecoInclusiveJetPtVsEtaKineCut");
    TH2D *hJetPtVsEtaMatchedKineCut = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEtaKineCut");
    TH2D *hJetPtVsEtaUnmatchedKineCut = (TH2D*)inFile->Get("hRecoInclusiveUnmatchedJetPtVsEtaKineCut");
    // Rebin 2D histograms
    hJetPtVsEtaKineCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaMatchedKineCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaUnmatchedKineCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hJetPtVsEtaMatchedKineCut->Divide(hJetPtVsEtaMatchedKineCut, hJetPtVsEtaKineCut, 1., 1., "b");
    hJetPtVsEtaUnmatchedKineCut->Divide(hJetPtVsEtaUnmatchedKineCut, hJetPtVsEtaKineCut, 1., 1., "b");
    // Set style
    set2DStyle(hJetPtVsEtaMatchedKineCut);
    set2DStyle(hJetPtVsEtaUnmatchedKineCut);

    //
    // Retrieve trkMax selection histograms
    //
    TH2D *hJetPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveJetPtVsEtaTrkMaxCut");
    TH2D *hJetPtVsEtaMatchedTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut");
    TH2D *hJetPtVsEtaUnmatchedTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut");
    // Rebin 2D histograms
    hJetPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaMatchedTrkMaxCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaUnmatchedTrkMaxCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hJetPtVsEtaMatchedTrkMaxCut->Divide(hJetPtVsEtaMatchedTrkMaxCut, hJetPtVsEtaTrkMaxCut, 1., 1., "b");
    hJetPtVsEtaUnmatchedTrkMaxCut->Divide(hJetPtVsEtaUnmatchedTrkMaxCut, hJetPtVsEtaTrkMaxCut, 1., 1., "b");
    // Set style
    set2DStyle(hJetPtVsEtaMatchedTrkMaxCut);
    set2DStyle(hJetPtVsEtaUnmatchedTrkMaxCut);

    TH2D *hLeadJetAllPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoLeadJetAllPtVsEta");
    TH2D *hLeadJetMatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoLeadJetMatchedPtVsEta");
    TH2D *hLeadJetUnmatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoLeadJetUnmatchedPtVsEta");
    TH2D *hSubLeadJetAllPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoSubLeadJetAllPtVsEta");
    TH2D *hSubLeadJetMatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoSubLeadJetMatchedPtVsEta");
    TH2D *hSubLeadJetUnmatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoSubLeadJetUnmatchedPtVsEta");
    // Rebin 2D histograms
    hLeadJetAllPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hLeadJetMatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hLeadJetUnmatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hLeadJetMatchedPtVsEtaTrkMaxCut->Divide(hLeadJetMatchedPtVsEtaTrkMaxCut, hLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    hLeadJetUnmatchedPtVsEtaTrkMaxCut->Divide(hLeadJetUnmatchedPtVsEtaTrkMaxCut, hLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    // Set style
    set2DStyle(hLeadJetMatchedPtVsEtaTrkMaxCut);
    set2DStyle(hLeadJetUnmatchedPtVsEtaTrkMaxCut);
    // Rebin 2D histograms
    hSubLeadJetAllPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetMatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetUnmatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hSubLeadJetMatchedPtVsEtaTrkMaxCut->Divide(hSubLeadJetMatchedPtVsEtaTrkMaxCut, hSubLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    hSubLeadJetUnmatchedPtVsEtaTrkMaxCut->Divide(hSubLeadJetUnmatchedPtVsEtaTrkMaxCut, hSubLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    // Set style
    set2DStyle(hSubLeadJetMatchedPtVsEtaTrkMaxCut);
    set2DStyle(hSubLeadJetUnmatchedPtVsEtaTrkMaxCut);


    //
    // Retrieve JetId selection histograms
    //
    TH2D *hJetPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveJetPtVsEtaJetIdCut");
    TH2D *hJetPtVsEtaMatchedJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEtaJetIdCut");
    TH2D *hJetPtVsEtaUnmatchedJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut");
    // Rebin 2D histograms
    hJetPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaMatchedJetIdCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaUnmatchedJetIdCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hJetPtVsEtaMatchedJetIdCut->Divide(hJetPtVsEtaMatchedJetIdCut, hJetPtVsEtaJetIdCut, 1., 1., "b");
    hJetPtVsEtaUnmatchedJetIdCut->Divide(hJetPtVsEtaUnmatchedJetIdCut, hJetPtVsEtaJetIdCut, 1., 1., "b");
    // Set style
    set2DStyle(hJetPtVsEtaMatchedJetIdCut);
    set2DStyle(hJetPtVsEtaUnmatchedJetIdCut);

    TH2D *hLeadJetAllPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoLeadJetAllPtVsEtaJetIdCut");
    TH2D *hLeadJetMatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoLeadJetMatchedPtVsEtaJetIdCut");
    TH2D *hLeadJetUnmatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoLeadJetUnmatchedPtVsEtaJetIdCut");
    TH2D *hSubLeadJetAllPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoSubLeadJetAllPtVsEtaJetIdCut");
    TH2D *hSubLeadJetMatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoSubLeadJetMatchedPtVsEtaJetIdCut");
    TH2D *hSubLeadJetUnmatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut");
    // Rebin 2D histograms
    hLeadJetAllPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hLeadJetMatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hLeadJetUnmatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hLeadJetMatchedPtVsEtaJetIdCut->Divide(hLeadJetMatchedPtVsEtaJetIdCut, hLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    hLeadJetUnmatchedPtVsEtaJetIdCut->Divide(hLeadJetUnmatchedPtVsEtaJetIdCut, hLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    // Set style
    set2DStyle(hLeadJetMatchedPtVsEtaJetIdCut);
    set2DStyle(hLeadJetUnmatchedPtVsEtaJetIdCut);
    // Rebin 2D histograms
    hSubLeadJetAllPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetMatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetUnmatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    // Divide 2D histograms
    hSubLeadJetMatchedPtVsEtaJetIdCut->Divide(hSubLeadJetMatchedPtVsEtaJetIdCut, hSubLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    hSubLeadJetUnmatchedPtVsEtaJetIdCut->Divide(hSubLeadJetUnmatchedPtVsEtaJetIdCut, hSubLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    // Set style
    set2DStyle(hSubLeadJetMatchedPtVsEtaJetIdCut);
    set2DStyle(hSubLeadJetUnmatchedPtVsEtaJetIdCut);


    // Retrieve 1D binning
    Int_t fPtBins = hJetPtVsEtaTrkMaxCut->GetNbinsY();
    Double_t fPtRange[2] {hJetPtVsEtaTrkMaxCut->GetYaxis()->GetBinLowEdge(1), 
                          hJetPtVsEtaTrkMaxCut->GetYaxis()->GetBinUpEdge(fPtBins) };
    Int_t fEtaBins = hJetPtVsEtaTrkMaxCut->GetNbinsX();
    Double_t fEtaRange[2] { hJetPtVsEtaTrkMaxCut->GetXaxis()->GetBinLowEdge(1),
                            hJetPtVsEtaTrkMaxCut->GetXaxis()->GetBinUpEdge(fEtaBins)};
    Double_t ptStep = (fPtRange[1]-fPtRange[0]) / fPtBins;
    Double_t etaStep = (fEtaRange[1]-fEtaRange[0]) / fEtaBins;

    std::cout << "fPtBins:  " << fPtBins << std::endl;
    std::cout << "fEtaBins: " << fEtaBins << std::endl;

    // Reserve histograms for projections

    // Kine selection
    TH1D *hJetMatchedPtKineCut[ fEtaBins ];
    TH1D *hJetMatchedEtaKineCut[ fPtBins ];
    TH1D *hJetUnmatchedPtKineCut[ fEtaBins ];
    TH1D *hJetUnmatchedEtaKineCut[ fPtBins ];

    TH1D *hJetMatchedPtRatioKineCut[ fEtaBins ];
    TH1D *hJetMatchedEtaRatioKineCut[ fPtBins ];
    TH1D *hJetUnmatchedPtRatioKineCut[ fEtaBins ];
    TH1D *hJetUnmatchedEtaRatioKineCut[ fPtBins ];

    // TrkMax selection
    TH1D *hJetMatchedPtTrkMaxCut[ fEtaBins ];
    TH1D *hJetMatchedEtaTrkMaxCut[ fPtBins ];
    TH1D *hJetUnmatchedPtTrkMaxCut[ fEtaBins ];
    TH1D *hJetUnmatchedEtaTrkMaxCut[ fPtBins ];

    TH1D *hLeadJetMatchedPtTrkMaxCut[ fEtaBins ];
    TH1D *hLeadJetMatchedEtaTrkMaxCut[ fPtBins ];
    TH1D *hLeadJetUnmatchedPtTrkMaxCut[ fEtaBins ];
    TH1D *hLeadJetUnmatchedEtaTrkMaxCut[ fPtBins ];

    TH1D *hSubLeadJetMatchedPtTrkMaxCut[ fEtaBins ];
    TH1D *hSubLeadJetMatchedEtaTrkMaxCut[ fPtBins ];
    TH1D *hSubLeadJetUnmatchedPtTrkMaxCut[ fEtaBins ];
    TH1D *hSubLeadJetUnmatchedEtaTrkMaxCut[ fPtBins ];

    // JetId selection
    TH1D *hJetMatchedPtJetIdCut[ fEtaBins ];
    TH1D *hJetMatchedEtaJetIdCut[ fPtBins ];
    TH1D *hJetUnmatchedPtJetIdCut[ fEtaBins ];
    TH1D *hJetUnmatchedEtaJetIdCut[ fPtBins ];

    TH1D *hJetMatchedPtRatioJetIdCut[ fEtaBins ];
    TH1D *hJetMatchedEtaRatioJetIdCut[ fPtBins ];
    TH1D *hJetUnmatchedPtRatioJetIdCut[ fEtaBins ];
    TH1D *hJetUnmatchedEtaRatioJetIdCut[ fPtBins ];

    TH1D *hLeadJetMatchedPtJetIdCut[ fEtaBins ];
    TH1D *hLeadJetMatchedEtaJetIdCut[ fPtBins ];
    TH1D *hLeadJetUnmatchedPtJetIdCut[ fEtaBins ];
    TH1D *hLeadJetUnmatchedEtaJetIdCut[ fPtBins ];

    TH1D *hSubLeadJetMatchedPtJetIdCut[ fEtaBins ];
    TH1D *hSubLeadJetMatchedEtaJetIdCut[ fPtBins ];
    TH1D *hSubLeadJetUnmatchedPtJetIdCut[ fEtaBins ];
    TH1D *hSubLeadJetUnmatchedEtaJetIdCut[ fPtBins ];


    //
    // Make projections
    //

    // Kine selection
    if ( plotKineCut ) {
        makeProjectionsFrom2D(hJetPtVsEtaMatchedKineCut, hJetMatchedEtaKineCut, hJetMatchedPtKineCut, fEtaBins, fPtBins, "hJetMatchedEtaKineCut", "hJetMatchedPtKineCut", kineStyle);
        makeProjectionsFrom2D(hJetPtVsEtaUnmatchedKineCut, hJetUnmatchedEtaKineCut, hJetUnmatchedPtKineCut, fEtaBins, fPtBins, "hJetUnmatchedEtaKineCut", "hJetUnmatchedPtKineCut", kineStyle+3);
    } // if ( plotKineCut )

    // TrkMax selection
    makeProjectionsFrom2D(hJetPtVsEtaMatchedTrkMaxCut, hJetMatchedEtaTrkMaxCut, hJetMatchedPtTrkMaxCut, fEtaBins, fPtBins, "hJetMatchedEtaTrkMaxCut", "hJetMatchedPtTrkMaxCut", trkMaxStyle);
    makeProjectionsFrom2D(hJetPtVsEtaUnmatchedTrkMaxCut, hJetUnmatchedEtaTrkMaxCut, hJetUnmatchedPtTrkMaxCut, fEtaBins, fPtBins, "hJetUnmatchedEtaTrkMaxCut", "hJetUnmatchedPtTrkMaxCut", trkMaxStyle+3);

    makeProjectionsFrom2D(hLeadJetMatchedPtVsEtaTrkMaxCut, hLeadJetMatchedEtaTrkMaxCut, hLeadJetMatchedPtTrkMaxCut, fEtaBins, fPtBins, "hLeadJetMatchedEtaTrkMaxCut", "hLeadJetMatchedPtTrkMaxCut", trkMaxStyle);
    makeProjectionsFrom2D(hLeadJetUnmatchedPtVsEtaTrkMaxCut, hLeadJetUnmatchedEtaTrkMaxCut, hLeadJetUnmatchedPtTrkMaxCut, fEtaBins, fPtBins, "hLeadJetUnmatchedEtaTrkMaxCut", "hLeadJetUnmatchedPtTrkMaxCut", trkMaxStyle+3);

    makeProjectionsFrom2D(hSubLeadJetMatchedPtVsEtaTrkMaxCut, hSubLeadJetMatchedEtaTrkMaxCut, hSubLeadJetMatchedPtTrkMaxCut, fEtaBins, fPtBins, "hSubLeadJetMatchedEtaTrkMaxCut", "hSubLeadJetMatchedPtTrkMaxCut", trkMaxStyle);
    makeProjectionsFrom2D(hSubLeadJetUnmatchedPtVsEtaTrkMaxCut, hSubLeadJetUnmatchedEtaTrkMaxCut, hSubLeadJetUnmatchedPtTrkMaxCut, fEtaBins, fPtBins, "hSubLeadJetUnmatchedEtaTrkMaxCut", "hSubLeadJetUnmatchedPtTrkMaxCut", trkMaxStyle+3);

    // JetId selection
    if ( plotJetIdCut ) {
        makeProjectionsFrom2D(hJetPtVsEtaMatchedJetIdCut, hJetMatchedEtaJetIdCut, hJetMatchedPtJetIdCut, fEtaBins, fPtBins, "hJetMatchedEtaJetIdCut", "hJetMatchedPtJetIdCut", jetIdStyle);
        makeProjectionsFrom2D(hJetPtVsEtaUnmatchedJetIdCut, hJetUnmatchedEtaJetIdCut, hJetUnmatchedPtJetIdCut, fEtaBins, fPtBins, "hJetUnmatchedEtaJetIdCut", "hJetUnmatchedPtJetIdCut", jetIdStyle+3);

        makeProjectionsFrom2D(hLeadJetMatchedPtVsEtaJetIdCut, hLeadJetMatchedEtaJetIdCut, hLeadJetMatchedPtJetIdCut, fEtaBins, fPtBins, "hLeadJetMatchedEtaJetIdCut", "hLeadJetMatchedPtJetIdCut", jetIdStyle);
        makeProjectionsFrom2D(hLeadJetUnmatchedPtVsEtaJetIdCut, hLeadJetUnmatchedEtaJetIdCut, hLeadJetUnmatchedPtJetIdCut, fEtaBins, fPtBins, "hLeadJetUnmatchedEtaJetIdCut", "hLeadJetUnmatchedPtJetIdCut", jetIdStyle+3);

        makeProjectionsFrom2D(hSubLeadJetMatchedPtVsEtaJetIdCut, hSubLeadJetMatchedEtaJetIdCut, hSubLeadJetMatchedPtJetIdCut, fEtaBins, fPtBins, "hSubLeadJetMatchedEtaJetIdCut", "hSubLeadJetMatchedPtJetIdCut", jetIdStyle);
        makeProjectionsFrom2D(hSubLeadJetUnmatchedPtVsEtaJetIdCut, hSubLeadJetUnmatchedEtaJetIdCut, hSubLeadJetUnmatchedPtJetIdCut, fEtaBins, fPtBins, "hSubLeadJetUnmatchedEtaJetIdCut", "hSubLeadJetUnmatchedPtJetIdCut", jetIdStyle+3);
    } // if ( plotJetIdCut )


    //
    // Make canvases
    //

    TCanvas *c = new TCanvas("c", "c", 1200, 800);

    TCanvas *cEtaProj = new TCanvas("cEtaProj", "cEtaProj", 1200, 800);
    cEtaProj->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtProj = new TCanvas("cPtProj", "cPtProj", 1200, 800);
    cPtProj->Divide(5, ( (fEtaBins % 5) == 0 ) ? (fEtaBins / 5) : (fEtaBins / 5 + 1) );

    TCanvas *cEtaProjRat = new TCanvas("cEtaProjRat", "cEtaProjRat", 1200, 800);
    cEtaProjRat->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtProjRat = new TCanvas("cPtProjRat", "cPtProjRat", 1200, 800);
    cPtProjRat->Divide(5, ( (fEtaBins % 5) == 0 ) ? (fEtaBins / 5) : (fEtaBins / 5 + 1) );

    // Create legends
    TLegend *hEtaLegend[fPtBins];
    TLegend *hPtLegend[fEtaBins];

    //
    // Plotting eta distributions
    //
    for (Int_t i{1}; i<=fPtBins; i++) {

        // Make ratios for the kine selection to trkMax one
        if ( plotKineCut ) {
            hJetMatchedEtaRatioKineCut[i-1] = dynamic_cast<TH1D*>( hJetMatchedEtaKineCut[i-1]->Clone( Form("hJetMatchedEtaRatioKineCut_%d", i-1) ) );
            make1DRatio(hJetMatchedEtaRatioKineCut[i-1], hJetMatchedEtaTrkMaxCut[i-1], "Ratio to TrkMax", kineStyle);
            hJetUnmatchedEtaRatioKineCut[i-1] = dynamic_cast<TH1D*>( hJetUnmatchedEtaKineCut[i-1]->Clone( Form("hJetUnmatchedEtaRatioKineCut_%d", i-1) ) );
            make1DRatio(hJetUnmatchedEtaRatioKineCut[i-1], hJetUnmatchedEtaTrkMaxCut[i-1], "Ratio to TrkMax", kineStyle+3);
        }

        // Make ratios for the jetId selection to trkMax one
        if ( plotJetIdCut ) {
            hJetMatchedEtaRatioJetIdCut[i-1] = dynamic_cast<TH1D*>( hJetMatchedEtaJetIdCut[i-1]->Clone( Form("hJetMatchedEtaRatioJetIdCut_%d", i-1) ) );
            make1DRatio(hJetMatchedEtaRatioJetIdCut[i-1], hJetMatchedEtaTrkMaxCut[i-1], "Ratio to TrkMax", jetIdStyle);
            hJetUnmatchedEtaRatioJetIdCut[i-1] = dynamic_cast<TH1D*>( hJetUnmatchedEtaJetIdCut[i-1]->Clone( Form("hJetUnmatchedEtaRatioJetIdCut_%d", i-1) ) );
            make1DRatio(hJetUnmatchedEtaRatioJetIdCut[i-1], hJetUnmatchedEtaTrkMaxCut[i-1], "Ratio to TrkMax", jetIdStyle+3);
        }

        // Plot individual distributions

        // Plot fakes and matched
        c->cd();
        setPadStyle();
        hJetMatchedEtaTrkMaxCut[i-1]->Draw();
        hJetMatchedEtaTrkMaxCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotJetIdCut ) {
            hJetMatchedEtaJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetMatchedEtaKineCut[i-1]->Draw("same");
        }
        hJetUnmatchedEtaTrkMaxCut[i-1]->Draw("same");
        if ( plotJetIdCut ) {
            hJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetUnmatchedEtaKineCut[i-1]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );

        hEtaLegend[i-1] = new TLegend(0.65, 0.35, 0.75, 0.55);
        hEtaLegend[i-1]->SetTextSize(0.04);
        hEtaLegend[i-1]->SetLineWidth(0);
        hEtaLegend[i-1]->AddEntry(hJetMatchedEtaTrkMaxCut[i-1], Form("TrkMax M"), "p");
        if ( plotKineCut ) {
            hEtaLegend[i-1]->AddEntry(hJetMatchedEtaKineCut[i-1], Form("KineOnly M"), "p");
        }
        if ( plotJetIdCut ) {
            hEtaLegend[i-1]->AddEntry(hJetMatchedEtaJetIdCut[i-1], Form("JetId M"), "p");
        }
        hEtaLegend[i-1]->AddEntry(hJetUnmatchedEtaTrkMaxCut[i-1], Form("TrkMax Un"), "p");
        if ( plotKineCut ) {
            hEtaLegend[i-1]->AddEntry(hJetUnmatchedEtaKineCut[i-1], Form("KineOnly Un"), "p");
        }
        if ( plotJetIdCut ) {
            hEtaLegend[i-1]->AddEntry(hJetUnmatchedEtaJetIdCut[i-1], Form("JetId Un"), "p");
        }
        hEtaLegend[i-1]->Draw();
        c->SaveAs( Form("%s/pPb8160_%s_eta_fakes_pt_%d_%d_%s.pdf", 
                        date.Data(), direction.Data(), 
                        (Int_t)(fPtRange[0] + (i-1) * ptStep),
                        (Int_t)(fPtRange[0] + i * ptStep),
                        branchName.Data()) );

        // All projections
        cEtaProj->cd(i);
        setPadStyle();
        hJetMatchedEtaTrkMaxCut[i-1]->Draw();
        hJetMatchedEtaTrkMaxCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotJetIdCut ) {
            hJetMatchedEtaJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetMatchedEtaKineCut[i-1]->Draw("same");
        }
        hJetUnmatchedEtaTrkMaxCut[i-1]->Draw("same");
        if ( plotJetIdCut ) {
            hJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetUnmatchedEtaKineCut[i-1]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );
        hEtaLegend[i-1]->Draw();


        // Plot ratios
        if ( plotJetIdCut ) {
            // Individual
            c->cd();
            setPadStyle();
            hJetMatchedEtaRatioJetIdCut[i-1]->Draw();
            hJetMatchedEtaRatioJetIdCut[i-1]->GetYaxis()->SetRangeUser(0.75, 1.25);
            if ( plotKineCut ) {
                hJetMatchedEtaRatioKineCut[i-1]->Draw("same");
            }
            hJetUnmatchedEtaRatioJetIdCut[i-1]->Draw("same");
            if ( plotKineCut ) {
                hJetUnmatchedEtaRatioKineCut[i-1]->Draw("same");
            }
            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                           fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );
            hEtaLegend[i-1]->Draw();
            c->SaveAs( Form("%s/pPb8160_%s_eta_fakes_ratio2trkMax_pt_%d_%d_%s.pdf", 
                       date.Data(), direction.Data(), 
                       (Int_t)(fPtRange[0] + (i-1) * ptStep),
                       (Int_t)(fPtRange[0] + i * ptStep),
                       branchName.Data()) );

            //  All projections
            cEtaProjRat->cd(i);
            setPadStyle();
            hJetMatchedEtaRatioJetIdCut[i-1]->Draw();
            hJetMatchedEtaRatioJetIdCut[i-1]->GetYaxis()->SetRangeUser(0.75, 1.25);
            if ( plotKineCut ) {
                hJetMatchedEtaRatioKineCut[i-1]->Draw("same");
            }
            hJetUnmatchedEtaRatioJetIdCut[i-1]->Draw("same");
            if ( plotKineCut ) {
                hJetUnmatchedEtaRatioKineCut[i-1]->Draw("same");
            }
            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                           fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );
            hEtaLegend[i-1]->Draw();            
        } // if ( plotJetIdCut )

    } // for (Int_t i{1}; i<fPtBins; i++)


    //
    // Plotting pt distributions
    //
    std::cout << "Eta bins: " << fEtaBins << std::endl;
    for (Int_t i{1}; i<=fEtaBins; i++) {
        std::cout << "Loop starts" << std::endl;

        // Make ratios for the kine selection to trkMax one
        if ( plotKineCut ) {
            std::cout << "Kine i: " << i << " ";
            hJetMatchedPtRatioKineCut[i-1] = dynamic_cast<TH1D*>( hJetMatchedPtKineCut[i-1]->Clone( Form("hJetMatchedPtRatioKineCut_%d", i-1) ) );
            std::cout << hJetMatchedPtRatioKineCut[i-1]->GetName();
            make1DRatio(hJetMatchedPtRatioKineCut[i-1], hJetMatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", kineStyle);
            std::cout << " " << hJetMatchedPtRatioKineCut[i-1]->Integral() << std::endl;
            hJetUnmatchedPtRatioKineCut[i-1] = dynamic_cast<TH1D*>( hJetUnmatchedPtKineCut[i-1]->Clone( Form("hJetUnmatchedPtRatioKineCut_%d", i-1) ) );
            std::cout << hJetUnmatchedPtRatioKineCut[i-1]->GetName();
            make1DRatio(hJetUnmatchedPtRatioKineCut[i-1], hJetUnmatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", kineStyle+3);
            std::cout << " " << hJetUnmatchedPtRatioKineCut[i-1]->Integral() << std::endl;
        }

        // Make ratios for the jetId selection to trkMax one
        if ( plotJetIdCut ) {
            std::cout << "JetId i: " << i << " ";
            hJetMatchedPtRatioJetIdCut[i-1] = dynamic_cast<TH1D*>( hJetMatchedPtJetIdCut[i-1]->Clone( Form("hJetMatchedPtRatioJetIdCut_%d", i-1) ) );
            std::cout << hJetMatchedPtRatioJetIdCut[i-1]->GetName();
            make1DRatio(hJetMatchedPtRatioJetIdCut[i-1], hJetMatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", jetIdStyle);
            std::cout << " " << hJetMatchedPtRatioJetIdCut[i-1]->Integral() << std::endl;
            hJetUnmatchedPtRatioJetIdCut[i-1] = dynamic_cast<TH1D*>( hJetUnmatchedPtJetIdCut[i-1]->Clone( Form("hJetUnmatchedPtRatioJetIdCut_%d", i-1) ) );
            make1DRatio(hJetUnmatchedPtRatioJetIdCut[i-1], hJetUnmatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", jetIdStyle+3);
        }

        // Plot individual distributions

        // Plot fakes and matched
        c->cd();
        setPadStyle();
        hJetMatchedPtTrkMaxCut[i-1]->Draw();
        hJetMatchedPtTrkMaxCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotJetIdCut ) {
            hJetMatchedPtJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetMatchedPtKineCut[i-1]->Draw("same");
        }
        hJetUnmatchedPtTrkMaxCut[i-1]->Draw("same");
        if ( plotJetIdCut ) {
            hJetUnmatchedPtJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetUnmatchedPtKineCut[i-1]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                       fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );

        hPtLegend[i-1] = new TLegend(0.65, 0.35, 0.75, 0.55);
        hPtLegend[i-1]->SetTextSize(0.04);
        hPtLegend[i-1]->SetLineWidth(0);
        hPtLegend[i-1]->AddEntry(hJetMatchedPtTrkMaxCut[i-1], Form("TrkMax M"), "p");
        if ( plotKineCut ) {
            hPtLegend[i-1]->AddEntry(hJetMatchedPtKineCut[i-1], Form("KineOnly M"), "p");
        }
        if ( plotJetIdCut ) {
            hPtLegend[i-1]->AddEntry(hJetMatchedPtJetIdCut[i-1], Form("JetId M"), "p");
        }
        hPtLegend[i-1]->AddEntry(hJetUnmatchedPtTrkMaxCut[i-1], Form("TrkMax Un"), "p");
        if ( plotKineCut ) {
            hPtLegend[i-1]->AddEntry(hJetUnmatchedPtKineCut[i-1], Form("KineOnly Un"), "p");
        }
        if ( plotJetIdCut ) {
            hPtLegend[i-1]->AddEntry(hJetUnmatchedPtJetIdCut[i-1], Form("JetId Un"), "p");
        }
        hPtLegend[i-1]->Draw();
        c->SaveAs( Form("%s/pPb8160_%s_pt_fakes_eta_%d_%d_%s.pdf", 
                        date.Data(), direction.Data(), 
                        (Int_t)( (fEtaRange[0] + (i-1) * etaStep) * 10),
                        (Int_t)( (fEtaRange[0] + i * etaStep) * 10 ),
                        branchName.Data()) );

        // All projections
        cPtProj->cd(i);
        setPadStyle();
        hJetMatchedPtTrkMaxCut[i-1]->Draw();
        hJetMatchedPtTrkMaxCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);
        if ( plotJetIdCut ) {
            hJetMatchedPtJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetMatchedPtKineCut[i-1]->Draw("same");
        }
        hJetUnmatchedPtTrkMaxCut[i-1]->Draw("same");
        if ( plotJetIdCut ) {
            hJetUnmatchedPtJetIdCut[i-1]->Draw("same");
        }
        if ( plotKineCut ) {
            hJetUnmatchedPtKineCut[i-1]->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                       fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );
        hPtLegend[i-1]->Draw();


        // Plot ratios
        if ( plotJetIdCut ) {
            // Individual
            c->cd();
            setPadStyle();
            hJetMatchedPtRatioJetIdCut[i-1]->Draw();
            hJetMatchedPtRatioJetIdCut[i-1]->GetYaxis()->SetRangeUser(0.75, 1.25);
            if ( plotKineCut ) {
                hJetMatchedPtRatioKineCut[i-1]->Draw("same");
            }
            hJetUnmatchedPtRatioJetIdCut[i-1]->Draw("same");
            if ( plotKineCut ) {
                hJetUnmatchedPtRatioKineCut[i-1]->Draw("same");
            }
            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                           fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );
            hPtLegend[i-1]->Draw();
            c->SaveAs( Form("%s/pPb8160_%s_pt_fakes_ratio2trkMax_eta_%d_%d_%s.pdf", 
                       date.Data(), direction.Data(), 
                       (Int_t)( (fEtaRange[0] + (i-1) * etaStep) * 10),
                       (Int_t)( (fEtaRange[0] + i * etaStep) * 10),
                       branchName.Data()) );

            //  All projections
            cPtProjRat->cd(i);
            setPadStyle();
            hJetMatchedPtRatioJetIdCut[i-1]->Draw();
            hJetMatchedPtRatioJetIdCut[i-1]->GetYaxis()->SetRangeUser(0.75, 1.25);
            if ( plotKineCut ) {
                hJetMatchedPtRatioKineCut[i-1]->Draw("same");
            }
            hJetUnmatchedPtRatioJetIdCut[i-1]->Draw("same");
            if ( plotKineCut ) {
                hJetUnmatchedPtRatioKineCut[i-1]->Draw("same");
            }
            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                           fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );
            hPtLegend[i-1]->Draw();

            std::cout << "Loop end" << std::endl;
        } // if ( plotJetIdCut )

    } // for (Int_t i{1}; i<fEtaBins; i++)



    cEtaProj->SaveAs( Form("%s/pPb8160_%s_eta_fake_projections_all_%s.pdf", 
                           date.Data(), direction.Data(), branchName.Data()) );
    cEtaProjRat->SaveAs( Form("%s/pPb8160_%s_eta_fake_ratios_all_%s.pdf", 
                              date.Data(), direction.Data(), branchName.Data()) );
    cPtProj->SaveAs( Form("%s/pPb8160_%s_pt_fake_projections_all_%s.pdf", 
                           date.Data(), direction.Data(), branchName.Data()) );
    cPtProjRat->SaveAs( Form("%s/pPb8160_%s_pt_fake_ratios_all_%s.pdf", 
                              date.Data(), direction.Data(), branchName.Data()) );

}

//________________
void pPb_embedding_qa(const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing_ak4.root") {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TDatime dt;
    TString date { Form( "%d",dt.GetDate() ) };
    TString inputFileName(inFileName);

    if ( directoryExists( date.Data() ) ) {
        //std::cout << "Directory exists." << std::endl;
    } 
    else {
        createDirectory( date.Data() );
    }

    TFile *inFile = TFile::Open(inFileName);
    if ( !inFile ) {
        std::cout << "[ERROR] Input file does not exist\n";
        return;
    }
    if ( inFile->IsZombie() ) {
        std::cout << "[ERROR] Beware of zombies!!!!\n";
        return;
    }

    Int_t branchId{0};
    if ( inputFileName.Contains("akCs4") ) {
        branchId = {0};
    }
    else {
        branchId = {1};
    }

    // Plot ptHat distribution
    //plotPtHat(inFile, date, branchId);

    // Compare inclusive reco, ref and gen transverse momentum spectra
    //compareInclusiveJetPtSpectra(inFile, date);

    // Plot jet reconstruction efficiency as a function of acceptance (pT vs eta)
    // plotEfficiency(inFile, date, branchId);

    // Plot dijet distributions
    //plotDijetDistributions(inFile, date);

    // Plot reco, reco with matching and calculate fakes
    plotRecoAndFakes(inFile, date, branchId);

    // Plot correlation between ref and reco dijet eta
    //plotEtaDijetCorrelation(inFile, date);

    // Plot distributions for jetId
    //plotJetIdHistos(inFile, date);

    // Plot JES and JER
    //plotJESandJER(inFile, date, branchId);
}