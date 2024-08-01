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
void set2DStyle(TH2* h, Bool_t doRenorm = kFALSE) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
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
void copy2DHisto(TH2* h, TH2* hOrig, Int_t xLow = 1, Int_t xHi = -1, Int_t yLow = 1, Int_t yHi = -1) {
    Int_t xBinLow{ xLow }, xBinHi{ hOrig->GetNbinsX() };
    if ( xHi > 0 && xHi < xBinHi ) {
        xBinHi = xHi;
    }

    Int_t yBinLow{ yLow }, yBinHi{ hOrig->GetNbinsY() };
    if ( yHi > 0 && yHi < yBinHi ) {
        yBinHi = yHi;
    }

    for (Int_t i=xBinLow; i<=xBinHi; i++) {
        for (Int_t j=yBinLow; j<=yBinHi; j++) {
            h->SetBinContent( i, j, hOrig->GetBinContent(i, j) );
            h->SetBinError( i, j, hOrig->GetBinError(i, j) );
        }
    }
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
    Bool_t plotKineCut{kFALSE};
    Bool_t plotJetIdCut{kFALSE};

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
void plotDijetResponseMatrices(TFile *embFile, TFile *mbFile, TFile *jet60File,
                               TFile *jet80File, TFile *jet100File, TString date) {

    TString inputFileName( embFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "pPbgoing";
    }

    TString frame = "lab";

    Int_t ptStep {5};
    Int_t ptLow {30};
    // Bin numbers
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194};
    // pTave values
    std::vector<Int_t> ptDijetPtValLow{};
    std::vector<Int_t> ptDijetPtValHi{};

    Int_t ptBins = ptDijetLow.size();

    for (UInt_t i{0}; i<ptDijetLow.size(); i++) {
        ptDijetPtValLow.push_back( ptLow + (ptDijetLow.at(i)-1) * ptStep );
        ptDijetPtValHi.push_back( ptLow + ptDijetHi.at(i) * ptStep );
    }

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};
    Int_t refSelType{6};
    Int_t dataType{2};
    Bool_t doRenorm{kFALSE};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Retrieve distributions for embedding
    THnSparseD *hReco2RefDijet = dynamic_cast< THnSparseD* >( embFile->Get("hRecoDijetPtEtaRefDijetPtEtaWeighted") );
    hReco2RefDijet->SetName("hReco2RefDijet");
    TH3D *hGenDijetPtEtaDphi = dynamic_cast< TH3D* >( embFile->Get("hGenDijetPtEtaDphiWeighted") );
    hGenDijetPtEtaDphi->SetName("hGenDijetPtEtaDphi");

    // Retrieve distributions for experimental data
    TH3D* hRecoDataMB = dynamic_cast< TH3D* > ( mbFile->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataMB->SetName("hRecoDataMB");
    TH3D* hRecoDataJet60 = dynamic_cast< TH3D* > ( jet60File->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataJet60->SetName("hRecoDataJet60");
    TH3D* hRecoDataJet80 = dynamic_cast< TH3D* > ( jet80File->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataJet80->SetName("hRecoDataJet80");
    TH3D* hRecoDataJet100 = dynamic_cast< TH3D* > ( jet100File->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataJet100->SetName("hRecoDataJet100");

    // Prepare histograms
    TH2D *hRecoPtVsRefPt;
    TH2D *hEmbRecoPtVsRecoEta;
    TH2D *hGenPtVsGenEta;
    TH2D *hDataRecoPtVsRecoEta[ ptDijetLow.size() + 1 ]; // pT bins + integrated
    TH2D *hRecoEtaVsRefEta[ ptDijetLow.size() + 1 ];     // pT bins + integrated
    TH2D *hEmbRefPtVsRefEta[ ptDijetLow.size() + 1 ];    // pT bins + integrated

    TH2D *hRecoPtVsRecoEtaMB = dynamic_cast<TH2D*> ( hRecoDataMB->Project3D("xy") );
    hRecoPtVsRecoEtaMB->SetName("hRecoPtVsRecoEtaMB");
    rescaleEta( hRecoPtVsRecoEtaMB );
    set2DStyle( hRecoPtVsRecoEtaMB, kTRUE );

    TH2D *hRecoPtVsRecoEtaJet60 = dynamic_cast<TH2D*> ( hRecoDataJet60->Project3D("xy") );
    hRecoPtVsRecoEtaJet60->SetName("hRecoPtVsRecoEtaJet60");
    rescaleEta( hRecoPtVsRecoEtaJet60 );
    set2DStyle( hRecoPtVsRecoEtaJet60, kTRUE );

    TH2D *hRecoPtVsRecoEtaJet80 = dynamic_cast<TH2D*> ( hRecoDataJet80->Project3D("xy") );
    hRecoPtVsRecoEtaJet80->SetName("hRecoPtVsRecoEtaJet80");
    rescaleEta( hRecoPtVsRecoEtaJet80 );
    set2DStyle( hRecoPtVsRecoEtaJet80, kTRUE );

    TH2D *hRecoPtVsRecoEtaJet100 = dynamic_cast<TH2D*> ( hRecoDataJet100->Project3D("xy") );
    hRecoPtVsRecoEtaJet100->SetName("hRecoPtVsRecoEtaJet100");
    rescaleEta( hRecoPtVsRecoEtaJet100 );
    set2DStyle( hRecoPtVsRecoEtaJet100, kTRUE );

    // Create canvases
    Int_t sizeX{1200};
    Int_t sizeY{800};

    TCanvas *canv = new TCanvas("canv", "canv", sizeX, sizeY);

    TCanvas *cRecoEtaVsRefEta = new TCanvas("cRecoEtaVsRefEta", "cRecoEtaVsRefEta", sizeX, sizeY);
    cRecoEtaVsRefEta->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cEmbRefPtVsRefEta = new TCanvas("cEmbRefPtVsRefEta", "cEmbRefPtVsRefEta", sizeX, sizeY);
    cEmbRefPtVsRefEta->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Fill and save reco pt vs ref pt for embedding
    hRecoPtVsRefPt = dynamic_cast<TH2D*> ( hReco2RefDijet->Projection(0, 2) );
    set2DStyle(hRecoPtVsRefPt);

    canv->cd();
    setPadStyle();
    hRecoPtVsRecoEtaMB->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/data_pPb8160_recoPtVsRecoEtaMB_%s.pdf", date.Data(), frame.Data() ) );

    canv->cd();
    setPadStyle();
    hRecoPtVsRecoEtaJet60->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/data_pPb8160_recoPtVsRecoEtaJet60_%s.pdf", date.Data(), frame.Data() ) );

    canv->cd();
    setPadStyle();
    hRecoPtVsRecoEtaJet80->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/data_pPb8160_recoPtVsRecoEtaJet80_%s.pdf", date.Data(), frame.Data() ) );

    canv->cd();
    setPadStyle();
    hRecoPtVsRecoEtaJet100->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/data_pPb8160_recoPtVsRecoEtaJet100_%s.pdf", date.Data(), frame.Data() ) );

    canv->cd();
    setPadStyle();
    hRecoPtVsRefPt->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/emb_pPb8160_recoPtVsRefPt_%s.pdf", date.Data(), frame.Data() ) );

    // Fill and save reco pt vs ref eta for embedding
    hEmbRecoPtVsRecoEta = dynamic_cast<TH2D*> ( hReco2RefDijet->Projection(0, 1) );
    set2DStyle(hEmbRecoPtVsRecoEta);

    canv->cd();
    setPadStyle();
    hEmbRecoPtVsRecoEta->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/emb_pPb8160_recoPtVsRecoEta_%s.pdf", date.Data(), frame.Data() ) );

    // Fill and save gen pt vs gen eta for embedding
    hGenPtVsGenEta = dynamic_cast<TH2D*> ( hGenDijetPtEtaDphi->Project3D("xy") );
    hGenPtVsGenEta->SetName( "hGenPtVsGenEta" );
    set2DStyle(hGenPtVsGenEta);

    canv->cd();
    setPadStyle();
    hGenPtVsGenEta->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/emb_pPb8160_genPtVsGenEta_%s.pdf", date.Data(), frame.Data() ) );

    // Loop over pTave bins
    for ( Int_t iPt{0}; iPt<ptDijetLow.size(); iPt++ ) {  

        // Select pTave bin for embedding
        hReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(iPt), ptDijetHi.at(iPt) );
        hRecoEtaVsRefEta[iPt] = dynamic_cast<TH2D*> ( hReco2RefDijet->Projection(3, 1) );
        hRecoEtaVsRefEta[iPt]->SetName( Form("hRecoEtaVsRefEta_%d", iPt) );
        hRecoEtaVsRefEta[iPt]->GetXaxis()->SetTitle("Reco #eta");
        hRecoEtaVsRefEta[iPt]->GetYaxis()->SetTitle("Ref #eta");
        set2DStyle(hRecoEtaVsRefEta[iPt]);

        hEmbRefPtVsRefEta[iPt] = dynamic_cast<TH2D*> ( hReco2RefDijet->Projection(2, 3) );
        hEmbRefPtVsRefEta[iPt]->SetName( Form("hEmbRefPtVsRefEta_%d", iPt) );
        hEmbRefPtVsRefEta[iPt]->GetXaxis()->SetTitle("Ref #eta");
        set2DStyle(hEmbRefPtVsRefEta[iPt]);

        if ( iPt == 0 ) {
            hRecoEtaVsRefEta[ ptDijetLow.size() ] = dynamic_cast<TH2D*> ( hRecoEtaVsRefEta[iPt]->Clone( Form("hRecoEtaVsRefEta_%d", ptDijetLow.size()) ) );
            hEmbRefPtVsRefEta[ ptDijetLow.size() ] = dynamic_cast<TH2D*> ( hEmbRefPtVsRefEta[iPt]->Clone( Form("hEmbRefPtVsRefEta_%d", ptDijetLow.size()) ) );
        }
        else {
            hRecoEtaVsRefEta[ ptDijetLow.size() ]->Add( hRecoEtaVsRefEta[iPt] );
            hEmbRefPtVsRefEta[ ptDijetLow.size() ]->Add( hEmbRefPtVsRefEta[iPt] );
        }

        // Rescale 2D distributions
        rescaleEta( hRecoEtaVsRefEta[iPt] );
        rescaleEta( hEmbRefPtVsRefEta[iPt] );

        // Retrieve experimental distributions
        if ( ptDijetPtValLow.at(iPt) < 80 ) {
            hDataRecoPtVsRecoEta[iPt] = dynamic_cast<TH2D*> ( hRecoPtVsRecoEtaMB->Clone( Form("hDataRecoPtVsRecoEta_%d", iPt) ) );
            hDataRecoPtVsRecoEta[iPt]->Reset();
            copy2DHisto(hDataRecoPtVsRecoEta[iPt], hRecoPtVsRecoEtaMB, 1, hRecoPtVsRecoEtaMB->GetNbinsX(), ptDijetLow.at(iPt), ptDijetHi.at(iPt) );
        }
        else if ( ptDijetPtValLow.at(iPt) < 100 ) {
            hDataRecoPtVsRecoEta[iPt] = dynamic_cast<TH2D*> ( hRecoPtVsRecoEtaJet60->Clone( Form("hDataRecoPtVsRecoEta_%d", iPt) ) );
            hDataRecoPtVsRecoEta[iPt]->Reset();
            copy2DHisto(hDataRecoPtVsRecoEta[iPt], hRecoPtVsRecoEtaJet60, 1, hRecoPtVsRecoEtaJet60->GetNbinsX(), ptDijetLow.at(iPt), ptDijetHi.at(iPt) );
        }
        else if ( ptDijetPtValLow.at(iPt) < 120 ) {
            hDataRecoPtVsRecoEta[iPt] = dynamic_cast<TH2D*> ( hRecoPtVsRecoEtaJet80->Clone( Form("hDataRecoPtVsRecoEta_%d", iPt) ) );
            hDataRecoPtVsRecoEta[iPt]->Reset();
            copy2DHisto(hDataRecoPtVsRecoEta[iPt], hRecoPtVsRecoEtaJet80, 1, hRecoPtVsRecoEtaJet80->GetNbinsX(), ptDijetLow.at(iPt), ptDijetHi.at(iPt) );
        }
        else {
            hDataRecoPtVsRecoEta[iPt] = dynamic_cast<TH2D*> ( hRecoPtVsRecoEtaJet100->Clone( Form("hDataRecoPtVsRecoEta_%d", iPt) ) );
            hDataRecoPtVsRecoEta[iPt]->Reset();
            copy2DHisto(hDataRecoPtVsRecoEta[iPt], hRecoPtVsRecoEtaJet100, 1, hRecoPtVsRecoEtaJet100->GetNbinsX(), ptDijetLow.at(iPt), ptDijetHi.at(iPt) );
        }
        hDataRecoPtVsRecoEta[iPt]->SetName( Form("hDataRecoPtVsRecoEta_%d", iPt) );
        set2DStyle(hDataRecoPtVsRecoEta[iPt]);

        if ( iPt == 0 ) {
            hDataRecoPtVsRecoEta[ ptDijetLow.size() ] = dynamic_cast<TH2D*> ( hDataRecoPtVsRecoEta[iPt]->Clone( Form("hDataRecoPtVsRecoEta_%d", ptDijetLow.size()) ) );
        }
        else {
            hDataRecoPtVsRecoEta[ ptDijetLow.size() ]->Add( hDataRecoPtVsRecoEta[iPt] );
        }

        // Rescale 2D distributions
        rescaleEta( hDataRecoPtVsRecoEta[iPt] );

        // Store individual ref eta vs reco eta projections for embedding (for the given reco pT)
        canv->cd();
        setPadStyle();
        hRecoEtaVsRefEta[iPt]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/emb_pPb8160_refEtaVsRecoEta_pt_%d_%d_%s.pdf", date.Data(),
                      ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                      ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

        // Store individual ref pt vs ref eta projections for embedding (for the given reco pT)
        canv->cd();
        setPadStyle();
        hEmbRefPtVsRefEta[iPt]->Draw("colz");
        hEmbRefPtVsRefEta[iPt]->GetYaxis()->SetRangeUser( (Double_t)ptDijetPtValLow.at(iPt) / 2, (Double_t)ptDijetPtValHi.at(iPt) * 2);
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/emb_pPb8160_refPtVsRefEta_pt_%d_%d_%s.pdf", date.Data(),
                      ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                      ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

        // Store individual reco pt vs reco eta projections for data (for the given reco pT)
        canv->cd();
        setPadStyle();
        hDataRecoPtVsRecoEta[iPt]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/data_pPb8160_recoPtVsRecoEta_pt_%d_%d_%s.pdf", date.Data(),
                      ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                      ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

        // Fill individual ref eta vs reco eta projections for embedding on one canvas (for the given reco pT)
        cRecoEtaVsRefEta->cd( iPt+1 );
        setPadStyle();
        hRecoEtaVsRefEta[iPt]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );

        // Fill individual ref pt vs ref eta projections for embedding (for the given reco pT)
        cEmbRefPtVsRefEta->cd( iPt+1 );
        setPadStyle();
        hEmbRefPtVsRefEta[iPt]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    } // for ( Int_t iPt{0}; iPt<ptDijetLow.size(); iPt++ )

    // Write the reco pt vs reco eta for the data
    canv->cd();
    setPadStyle();
    rescaleEta( hDataRecoPtVsRecoEta[ ptDijetLow.size() ] );
    hDataRecoPtVsRecoEta[ ptDijetLow.size() ]->Draw("colz");
    t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    canv->SaveAs( Form("%s/data_pPb8160_recoPtVsRecoEta_%s.pdf", date.Data(), frame.Data() ) );
    
    cRecoEtaVsRefEta->SaveAs( Form("%s/emb_pPb8160_recoEtaVsRefEta_all_%s.pdf", date.Data(), frame.Data() ) );
    cEmbRefPtVsRefEta->SaveAs( Form("%s/emb_pPb8160_refPtVsRefEta_all_%s.pdf", date.Data(), frame.Data() ) );
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
        direction = "Pb-going";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "p-going";
    }
    else {
        direction = "pPbgoing";
    }

    Int_t rebinX{1}, rebinY{2};

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

    Int_t etaLow{20}, etaHi{32};
    Double_t etaStep = 0.2;
    Double_t etaStart{-5.2};

    // Get JES vs ref pT
    hJESPtEtaPhi->GetAxis(2)->SetRange(etaLow, etaHi); // -1.2 < eta < 1.2
    TH2D* hJESvsPt = (TH2D*)hJESPtEtaPhi->Projection(0, 1);
    hJESvsPt->SetName("hJESvsPt");
    hJESvsPt->RebinX( rebinX );
    hJESvsPt->RebinY( rebinY );
    hJESvsPt->GetYaxis()->SetTitle("p_{T}^{reco} / p_{T}^{ref}");
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
    t.DrawLatexNDC(0.65, 0.75, direction.Data() );
    t.DrawLatexNDC(0.6, 0.65, Form("%2.1f<#eta<%2.1f", etaStart + etaStep * etaLow, etaStart + etaStep * etaHi) );
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
    hJESMean->GetYaxis()->SetRangeUser(0.95, 1.05);
    hJESMean->GetYaxis()->SetTitle("JES");
    //hJESMean->GetXaxis()->SetRangeUser(20, 200);
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    t.DrawLatexNDC(0.65, 0.75, direction.Data() );
    t.DrawLatexNDC(0.6, 0.65, Form("%2.1f<#eta<%2.1f", etaStart + etaStep * etaLow, etaStart + etaStep * etaHi) );
    canv->SaveAs( Form("%s/pPb8160_%s_jes_mu_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );

    hJESSigma->Draw();
    //hJESSigma->GetXaxis()->SetRangeUser(20, 800);
    hJESSigma->GetYaxis()->SetTitle("JER");
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    t.DrawLatexNDC(0.65, 0.75, direction.Data() );
    t.DrawLatexNDC(0.6, 0.65, Form("%2.1f<#eta<%2.1f", etaStart + etaStep * etaLow, etaStart + etaStep * etaHi) );
    canv->SaveAs( Form("%s/pPb8160_%s_jes_sigma_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );

    // Return eta to original area
    hJESPtEtaPhi->GetAxis(2)->SetRange(1, 52);

    // Plot 2D, mu and sigma on a single canvas
    TCanvas *canv2 = new TCanvas("canv2", "canv2", 1200, 400);
    canv2->Divide(3, 1);

    canv2->cd(1);
    setPadStyle();
    hJESvsPt->Draw("colz");
    hJESvsPt->GetYaxis()->SetTitle("p_{T}^{reco} / p_{T}^{ref}");
    //gPad->SetLogz();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    t.DrawLatexNDC(0.65, 0.75, direction.Data() );
    t.DrawLatexNDC(0.6, 0.65, Form("%2.1f<#eta<%2.1f", etaStart + etaStep * etaLow, etaStart + etaStep * etaHi) );

    canv2->cd(2);
    setPadStyle();
    hJESMean->Draw();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    t.DrawLatexNDC(0.65, 0.75, direction.Data() );
    t.DrawLatexNDC(0.6, 0.65, Form("%2.1f<#eta<%2.1f", etaStart + etaStep * etaLow, etaStart + etaStep * etaHi) );

    canv2->cd(3);
    setPadStyle();
    hJESSigma->Draw();
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    t.DrawLatexNDC(0.65, 0.75, direction.Data() );
    t.DrawLatexNDC(0.6, 0.65, Form("%2.1f<#eta<%2.1f", etaStart + etaStep * etaLow, etaStart + etaStep * etaHi) );
    
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
        hJESvsEta[i]->GetYaxis()->SetTitle("p_{T}^{reco} / p_{T}^{ref}");
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
        t.DrawLatexNDC(0.65, 0.75, direction.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );
        canv->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_pt_%d_%d_%s.pdf", 
                           date.Data(),  direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                           ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );

        canv->cd();
        setPadStyle();
        hJESvsEtaMean[i]->Draw();
        hJESvsEtaMean[i]->GetYaxis()->SetRangeUser(0.95, 1.05);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.65, 0.75, direction.Data() );
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
        t.DrawLatexNDC(0.65, 0.75, direction.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );
        canv->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_mu_pt_%d_%d_%s.pdf", 
                           date.Data(), direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                           ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );

        canv3[i]->cd(1);
        setPadStyle();
        hJESvsEta[i]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.65, 0.75, direction.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );

        canv3[i]->cd(2);
        setPadStyle();
        hJESvsEtaMean[i]->Draw();
        hJESvsEtaMean[i]->GetYaxis()->SetRangeUser(0.95, 1.05);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.65, 0.75, direction.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );

        canv3[i]->cd(3);
        setPadStyle();
        hJESvsEtaSigma[i]->Draw();
        hJESvsEtaSigma[i]->GetYaxis()->SetRangeUser(0., 0.2);
        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.65, 0.75, direction.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%d<p_{T} (GeV/c)<%d", ptFirst + (ptLow.at(i)-1) * ptStep, ptFirst + ptHi.at(i) * ptStep ) );

        canv3[i]->SaveAs( Form("%s/pPb8160_%s_jes_vs_eta_all_pt_%d_%d_%s.pdf", 
                               date.Data(), direction.Data(), ptFirst + (ptLow.at(i)-1) * ptStep, 
                               ptFirst + ptHi.at(i) * ptStep, branchName.Data()) );
    }

}

//________________
void plotDijetDistributions(TFile *embFile, TFile *mbFile, TFile *jet60File,
                            TFile *jet80File, TFile *jet100File, TString date) {

    bool useUncrt{true};
    TFile *uncrtFile;
    if ( useUncrt ) {
        uncrtFile = TFile::Open("freezeSyst/oSystematics.root");
    }

    TString inputFileName( embFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "pPbgoing";
    }

    TString frame = "lab";

    Int_t ptStep {5};
    Int_t ptLow {30};
    // Bin numbers
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194};
    // pTave values
    std::vector<Int_t> ptDijetPtValLow{};
    std::vector<Int_t> ptDijetPtValHi{};

    for (UInt_t i{0}; i<ptDijetLow.size(); i++) {
        ptDijetPtValLow.push_back( ptLow + (ptDijetLow.at(i)-1) * ptStep );
        ptDijetPtValHi.push_back( ptLow + ptDijetHi.at(i) * ptStep );
    }
    // std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    // std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194};

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};
    Int_t refSelType{6};
    Int_t dataType{2};
    Bool_t doRenorm{kFALSE};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Retrieve distributions for embedding
    THnSparseD *hReco2RefDijet = dynamic_cast< THnSparseD* >( embFile->Get("hRecoDijetPtEtaRefDijetPtEtaWeighted") );
    hReco2RefDijet->SetName("hReco2RefDijet");
    TH3D *hGenDijetPtEtaDphi = dynamic_cast< TH3D* >( embFile->Get("hGenDijetPtEtaDphiWeighted") );
    hGenDijetPtEtaDphi->SetName("hGenDijetPtEtaDphi");

    // Retrieve distributions for experimental data
    TH3D* hRecoDataMB = dynamic_cast< TH3D* > ( mbFile->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataMB->SetName("hRecoDataMB");
    TH3D* hRecoDataJet60 = dynamic_cast< TH3D* > ( jet60File->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataJet60->SetName("hRecoDataJet60");
    TH3D* hRecoDataJet80 = dynamic_cast< TH3D* > ( jet80File->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataJet80->SetName("hRecoDataJet80");
    TH3D* hRecoDataJet100 = dynamic_cast< TH3D* > ( jet100File->Get("hRecoDijetPtEtaDphiWeighted") );
    hRecoDataJet100->SetName("hRecoDataJet100");

    // Create 1D dijet eta histograms
    TH1D *hRecoEta[ ptDijetLow.size() ];
    TH1D *hRecoDataEta[ ptDijetLow.size() ];
    TH1D *hRefEta[ ptDijetLow.size() ];
    TH1D *hGenEta[ ptDijetLow.size() ];
    TH1D *hRecoEtaDivNeighborEmb[ ptDijetLow.size() - 1 ];
    TH1D *hRecoEtaDivNeighborData[ ptDijetLow.size() - 1 ];

    TH1D *hReco2GenRatio[ ptDijetLow.size() ];
    TH1D *hRef2GenRatio[ ptDijetLow.size() ];
    TH1D *hRecoData2MCRatio[ ptDijetLow.size() ];
    TH1D *hRelSystUncrt[ ptDijetLow.size() ];

    // Create line and legend
    TLine *line;
    TLegend *leg;

    // Create canvases and pads
    Int_t sizeX{1200};
    Int_t sizeY{800};

    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRatData2Mc = new TCanvas("cRatData2Mc", "cRatData2Mc", sizeX, sizeY);
    cRatData2Mc->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRatNeighbor = new TCanvas("cRatNeighbor", "cRatNeighbor", sizeX, sizeY);
    cRatNeighbor->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over pTave bins
    for ( Int_t iPt{0}; iPt<ptDijetLow.size(); iPt++ ) {

        // Retrieve gen distribution
        hGenEta[iPt] = dynamic_cast<TH1D*>( hGenDijetPtEtaDphi->ProjectionY( Form("hGenEta_%d", iPt), ptDijetLow.at(iPt), ptDijetHi.at(iPt) ) );
        hGenEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        rescaleEta( hGenEta[iPt] );
        set1DStyle( hGenEta[iPt], genType );


        // Retrieve reco and ref eta distributions for the given reco pTave bin
        hReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(iPt), ptDijetHi.at(iPt) );

        // Project on reco eta
        hRecoEta[iPt] = dynamic_cast<TH1D*> ( hReco2RefDijet->Projection(1) );
        hRecoEta[iPt]->SetName( Form("hRecoEta_%d", iPt) );
        hRecoEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        rescaleEta( hRecoEta[iPt] );
        set1DStyle( hRecoEta[iPt], recoType );

        // Project on ref eta
        hRefEta[iPt] = dynamic_cast<TH1D*> ( hReco2RefDijet->Projection(3) );
        hRefEta[iPt]->SetName( Form("hRefEta_%d", iPt) );
        hRefEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        rescaleEta( hRefEta[iPt] );
        set1DStyle( hRefEta[iPt], refType );

        // Make ratio of reco to gen distribution
        hReco2GenRatio[iPt] = dynamic_cast<TH1D*>( hRecoEta[iPt]->Clone( Form("hReco2GenRatio_%d", iPt) ) );
        hReco2GenRatio[iPt]->Divide( hReco2GenRatio[iPt], hGenEta[iPt], 1., 1. /*, "b" */ );
        hReco2GenRatio[iPt]->GetXaxis()->SetTitle("#eta^{dijet}");
        hReco2GenRatio[iPt]->GetYaxis()->SetTitle("Ratio to gen");
        
        if ( useUncrt ) {
            hRelSystUncrt[iPt] = dynamic_cast<TH1D*> ( uncrtFile->Get( Form("hTotalSystRel_%d", iPt) ) );
            hRelSystUncrt[iPt]->GetXaxis()->SetTitle("#eta^{dijet}");
            hRelSystUncrt[iPt]->GetYaxis()->SetTitle("Ratio to gen");
        }

        // Make ratio of ref to gen distribution
        hRef2GenRatio[iPt] = dynamic_cast<TH1D*>( hRefEta[iPt]->Clone( Form("hRef2GenRatio_%d", iPt) ) );
        hRef2GenRatio[iPt]->Divide( hRef2GenRatio[iPt], hGenEta[iPt], 1., 1. /*, "b" */ );
        hRef2GenRatio[iPt]->GetXaxis()->SetTitle("#eta^{dijet}");
        hRef2GenRatio[iPt]->GetYaxis()->SetTitle("Ratio to gen");

        // Retrieve experimental distributions
        if (ptDijetPtValLow.at(iPt) < 80) {
            hRecoDataEta[iPt] = dynamic_cast<TH1D*> ( hRecoDataMB->ProjectionY( Form("hRecoDataEta_%d", iPt), ptDijetLow.at(iPt), ptDijetHi.at(iPt) ) );
            hRecoDataEta[iPt]->GetXaxis()->SetTitle( "#eta^{dijet}" );
            hRecoDataEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
            rescaleEta( hRecoDataEta[iPt] );
            set1DStyle( hRecoDataEta[iPt], dataType );
            // hRecoDataEta[iPt]->Draw();

        }
        else if (ptDijetPtValLow.at(iPt) < 100) {
            hRecoDataEta[iPt] = dynamic_cast<TH1D*> ( hRecoDataJet60->ProjectionY( Form("hRecoDataEta_%d", iPt), ptDijetLow.at(iPt), ptDijetHi.at(iPt) ) );
            hRecoDataEta[iPt]->GetXaxis()->SetTitle( "#eta^{dijet}" );
            hRecoDataEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
            rescaleEta( hRecoDataEta[iPt] );
            set1DStyle( hRecoDataEta[iPt], dataType );
        }
        else if (ptDijetPtValLow.at(iPt) < 120) {
            hRecoDataEta[iPt] = dynamic_cast<TH1D*> ( hRecoDataJet80->ProjectionY( Form("hRecoDataEta_%d", iPt), ptDijetLow.at(iPt), ptDijetHi.at(iPt) ) );
            hRecoDataEta[iPt]->GetXaxis()->SetTitle( "#eta^{dijet}" );
            hRecoDataEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
            rescaleEta( hRecoDataEta[iPt] );
            set1DStyle( hRecoDataEta[iPt], dataType );
        }
        else {
            hRecoDataEta[iPt] = dynamic_cast<TH1D*> ( hRecoDataJet100->ProjectionY( Form("hRecoDataEta_%d", iPt), ptDijetLow.at(iPt), ptDijetHi.at(iPt) ) );
            hRecoDataEta[iPt]->GetXaxis()->SetTitle( "#eta^{dijet}" );
            hRecoDataEta[iPt]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
            rescaleEta( hRecoDataEta[iPt] );
            set1DStyle( hRecoDataEta[iPt], dataType );
        }

        // Make ratio of reco data to reco MC
        hRecoData2MCRatio[iPt] = dynamic_cast<TH1D*> ( hRecoDataEta[iPt]->Clone( Form("hRecoData2MCRatio_%d", iPt) ) );
        hRecoData2MCRatio[iPt]->Divide(hRecoData2MCRatio[iPt], hRecoEta[iPt], 1., 1.);
        std::cout << "data nbins: " <<  hRecoData2MCRatio[iPt]->GetNbinsX() << " embed nbins: " << hRecoEta[iPt]->GetNbinsX() << std::endl;
        hRecoData2MCRatio[iPt]->GetXaxis()->SetTitle( "Reco #eta^{dijet}" );
        set1DStyle(hRecoData2MCRatio[iPt], dataType);

        // Make ratio of the reco distribution to the neighbor (previous one)
        if ( iPt != 0 ) {
            hRecoEtaDivNeighborEmb[iPt-1] = dynamic_cast<TH1D*> ( hRecoEta[iPt]->Clone( Form("hRecoEtaDivNeighborEmb_%d", iPt-1)  ) );
            hRecoEtaDivNeighborEmb[iPt-1]->Divide( hRecoEtaDivNeighborEmb[iPt-1], hRecoEta[iPt-1], 1., 1.);
            hRecoEtaDivNeighborEmb[iPt-1]->GetYaxis()->SetTitle("Ratio to neighboring p_{T} bin");
            set1DStyle(hRecoEtaDivNeighborEmb[iPt-1], recoType);

            hRecoEtaDivNeighborData[iPt-1] = dynamic_cast<TH1D*> ( hRecoDataEta[iPt]->Clone( Form("hRecoEtaDivNeighborData_%d", iPt-1)  ) );
            hRecoEtaDivNeighborData[iPt-1]->Divide( hRecoEtaDivNeighborData[iPt-1], hRecoDataEta[iPt-1], 1., 1.);
            hRecoEtaDivNeighborData[iPt-1]->GetYaxis()->SetTitle("Ratio to neighboring p_{T} bin");
            set1DStyle(hRecoEtaDivNeighborData[iPt-1], dataType);
        }

        // Plot comparison of gen, reco and ref
        canv->cd();
        setPadStyle();
        hGenEta[iPt]->Draw();
        hRecoEta[iPt]->Draw("same");
        hRefEta[iPt]->Draw("same");
        hRecoDataEta[iPt]->Draw("same");
        hGenEta[iPt]->GetYaxis()->SetRangeUser(0., 0.15);
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hGenEta[iPt], Form("Gen"), "p");
        leg->AddEntry(hRecoEta[iPt], Form("Reco"), "p");
        leg->AddEntry(hRefEta[iPt], Form("Ref"), "p");
        leg->AddEntry(hRecoDataEta[iPt], Form("Data"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/emb_pPb8160_etaDijet_RecoRefGenComp_pt_%d_%d_%s.pdf", date.Data(),
                      ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                      ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

        // Plot ratios of reco and ref to gen
        canv->cd();
        setPadStyle();
        hRelSystUncrt[iPt]->Draw("E2");
        hReco2GenRatio[iPt]->Draw("same");
        hRef2GenRatio[iPt]->Draw("same");
        hRelSystUncrt[iPt]->GetYaxis()->SetRangeUser(0.8, 1.5);
        hRelSystUncrt[iPt]->GetXaxis()->SetRangeUser(-3., 3.);
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hReco2GenRatio[iPt], Form("Reco / Gen"), "p");
        leg->AddEntry(hRef2GenRatio[iPt], Form("Ref / Gen"), "p");
        leg->Draw();
        line = new TLine(-3.01, 1., 3.01, 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/emb_pPb8160_etaDijet_RecoRefGenRat_pt_%d_%d_%s.pdf", date.Data(),
                      ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                      ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

        // Plot ratios of reco data to reco embedding
        canv->cd();
        setPadStyle();
        hRelSystUncrt[iPt]->Draw("E2");
        hRecoData2MCRatio[iPt]->Draw("same");
        hReco2GenRatio[iPt]->Draw("same");
        hRef2GenRatio[iPt]->Draw("same");
        hRelSystUncrt[iPt]->GetYaxis()->SetRangeUser(0.7, 1.5);
        hRelSystUncrt[iPt]->GetXaxis()->SetRangeUser(-3., 3.);
        hRecoData2MCRatio[iPt]->GetYaxis()->SetTitle("Data / Embed");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hRecoData2MCRatio[iPt], Form("Reco(data) / Reco(MC)"), "p");
        leg->AddEntry(hReco2GenRatio[iPt], Form("Reco(MC) / Gen"), "p");
        leg->AddEntry(hRef2GenRatio[iPt], Form("Ref(Ref axis) / Gen"), "p");
        leg->AddEntry(hRelSystUncrt[iPt], Form("Syst. uncrt."), "f");
        leg->Draw();
        line = new TLine(-3.01, 1., 3.01, 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_data2embed_pt_%d_%d_%s.pdf", date.Data(),
                      ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                      ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

        // Plot comparisons on a single canvas
        cComp->cd( iPt+1 );
        setPadStyle();
        hGenEta[iPt]->Draw();
        hRecoEta[iPt]->Draw("same");
        hRefEta[iPt]->Draw("same");
        hRecoDataEta[iPt]->Draw("same");
        hGenEta[iPt]->GetYaxis()->SetRangeUser(0., 0.15);
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hGenEta[iPt], Form("Gen"), "p");
        leg->AddEntry(hRecoEta[iPt], Form("Reco"), "p");
        leg->AddEntry(hRefEta[iPt], Form("Ref"), "p");
        leg->AddEntry(hRecoDataEta[iPt], Form("Data"), "p");
        leg->Draw();

        // Plot ratios or reco emb to gen and ref emb to gen on a single canvas
        cRat->cd( iPt+1 );
        setPadStyle();
        hReco2GenRatio[iPt]->Draw();
        hRef2GenRatio[iPt]->Draw("same");
        hReco2GenRatio[iPt]->GetYaxis()->SetRangeUser(0.8, 1.5);
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hReco2GenRatio[iPt], Form("Reco / Gen"), "p");
        leg->AddEntry(hRef2GenRatio[iPt], Form("Ref / Gen"), "p");
        leg->Draw();
        line = new TLine(hReco2GenRatio[iPt]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hReco2GenRatio[iPt]->GetXaxis()->GetBinUpEdge(hReco2GenRatio[iPt]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();

        // Plot ratios of experimental data over reco MC
        cRatData2Mc->cd( iPt+1 );
        setPadStyle();
        hRelSystUncrt[iPt]->Draw("E2");
        hRecoData2MCRatio[iPt]->Draw("same");
        hReco2GenRatio[iPt]->Draw("same");
        hRef2GenRatio[iPt]->Draw("same");
        hRelSystUncrt[iPt]->GetYaxis()->SetRangeUser(0.7, 1.2);
        hRelSystUncrt[iPt]->GetXaxis()->SetRangeUser(-3., 3.);
        hRecoData2MCRatio[iPt]->GetYaxis()->SetTitle("Data / MC");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.4);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hRecoData2MCRatio[iPt], Form("Reco(data) / Reco(MC)"), "p");
        leg->AddEntry(hReco2GenRatio[iPt], Form("Reco(MC) / Gen"), "p");
        leg->AddEntry(hRef2GenRatio[iPt], Form("Ref(Ref axis) / Gen"), "p");
        leg->AddEntry(hRelSystUncrt[iPt], Form("Syst. uncrt."), "f");
        leg->Draw();
        line = new TLine(-3., 1., 3., 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();

        // Make ratio of the reco distribution to the neighbor (previous one) in one canvas
        if ( iPt != 0 && iPt != (ptDijetLow.size()-1) ) {

            // Individual distributions
            canv->cd();
            setPadStyle();
            hRecoEtaDivNeighborEmb[iPt-1]->Draw();
            hRecoEtaDivNeighborData[iPt-1]->Draw("same");
            hRecoEtaDivNeighborEmb[iPt-1]->GetYaxis()->SetRangeUser(0.6, 1.3);
            t.DrawLatexNDC(0.15, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d / %d < p_{T}^{ave} (GeV) < %d", 
                           ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep,
                           ptLow + (ptDijetLow.at(iPt-1) - 1) * ptStep, ptLow + ptDijetHi.at(iPt-1) * ptStep) );
            leg = new TLegend(0.2, 0.65, 0.4, 0.85);
            leg->SetTextSize(0.04);
            leg->SetLineWidth(0);
            leg->AddEntry(hRecoEtaDivNeighborEmb[iPt-1], Form("Reco embed"), "p");
            leg->AddEntry(hRecoEtaDivNeighborData[iPt-1], Form("Reco data"), "p");
            leg->Draw();
            line = new TLine(hRecoEtaDivNeighborEmb[iPt-1]->GetXaxis()->GetBinLowEdge(1), 1., 
                             hRecoEtaDivNeighborEmb[iPt-1]->GetXaxis()->GetBinUpEdge(hRecoEtaDivNeighborEmb[iPt-1]->GetNbinsX()), 1.);
            line->SetLineColor(kMagenta);
            line->SetLineWidth(3);
            line->SetLineStyle(3);
            line->Draw();
            canv->SaveAs( Form("%s/emb_pPb8160_etaDijet_recoRatNeighbor_pt_%d_%d_%s.pdf", date.Data(),
                          ptLow + (ptDijetLow.at(iPt)-1) * ptStep, 
                          ptLow + ptDijetHi.at(iPt) * ptStep, frame.Data() ) );

            // Plot all on one canvas
            cRatNeighbor->cd( iPt );
            setPadStyle();
            hRecoEtaDivNeighborEmb[iPt-1]->Draw();
            hRecoEtaDivNeighborData[iPt-1]->Draw("same");
            hRecoEtaDivNeighborEmb[iPt-1]->GetYaxis()->SetRangeUser(0.6, 1.3);
            t.DrawLatexNDC(0.15, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d / %d < p_{T}^{ave} (GeV) < %d", 
                           ptLow + (ptDijetLow.at(iPt) - 1) * ptStep, ptLow + ptDijetHi.at(iPt) * ptStep,
                           ptLow + (ptDijetLow.at(iPt-1) - 1) * ptStep, ptLow + ptDijetHi.at(iPt-1) * ptStep) );
            leg = new TLegend(0.2, 0.65, 0.4, 0.85);
            leg->SetTextSize(0.04);
            leg->SetLineWidth(0);
            leg->AddEntry(hRecoEtaDivNeighborEmb[iPt-1], Form("Reco embed"), "p");
            leg->AddEntry(hRecoEtaDivNeighborData[iPt-1], Form("Reco data"), "p");
            leg->Draw();
            line = new TLine(hRecoEtaDivNeighborEmb[iPt-1]->GetXaxis()->GetBinLowEdge(1), 1., 
                             hRecoEtaDivNeighborEmb[iPt-1]->GetXaxis()->GetBinUpEdge(hRecoEtaDivNeighborEmb[iPt-1]->GetNbinsX()), 1.);
            line->SetLineColor(kMagenta);
            line->SetLineWidth(3);
            line->SetLineStyle(3);
            line->Draw();
        }   
    } // for ( Int_t iPt{0}; iPt<ptDijetLow.size(); iPt++ )

    cComp->SaveAs( Form("%s/emb_pPb8160_etaDijet_RecoRefGenComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/emb_pPb8160_etaDijet_RecoRefGenRat_all_%s.pdf", date.Data(), frame.Data() ) );
    cRatNeighbor->SaveAs( Form("%s/emb_pPb8160_etaDijet_recoRatNeighbor_all_%s.pdf", date.Data(), frame.Data() ) );
    cRatData2Mc->SaveAs( Form("%s/pPb8160_etaDijet_data2embed_all_%s.pdf", date.Data(), frame.Data() ) );
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
    Int_t rebinX{2}, rebinY{2};

    // Plotting options
    Bool_t plotKineCut{kFALSE};
    Bool_t plotJetIdCut{kFALSE};

    // Make latex
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Color style
    Int_t kineStyle{1};
    Int_t trkMaxStyle{2};
    Int_t jetIdStyle{0};

    //
    // Retrieve kine selection histograms
    //
    TH2D *hJetPtVsEtaKineCut = (TH2D*)inFile->Get("hRecoInclusiveJetPtVsEtaKineCut");
    TH2D *hJetPtVsEtaMatchedKineCut = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEtaKineCut");
    TH2D *hJetPtVsEtaUnmatchedKineCut = (TH2D*)inFile->Get("hRecoInclusiveUnmatchedJetPtVsEtaKineCut");

    //
    // Retrieve trkMax selection histograms
    //
    TH2D *hJetPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveJetPtVsEtaTrkMaxCut");
    TH2D *hJetPtVsEtaMatchedTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut");
    TH2D *hJetPtVsEtaUnmatchedTrkMaxCut = (TH2D*)inFile->Get("hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut");
    TH2D *hLeadJetAllPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoLeadJetAllPtVsEta");
    TH2D *hLeadJetMatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoLeadJetMatchedPtVsEta");
    TH2D *hLeadJetUnmatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoLeadJetUnmatchedPtVsEta");
    TH2D *hSubLeadJetAllPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoSubLeadJetAllPtVsEta");
    TH2D *hSubLeadJetMatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoSubLeadJetMatchedPtVsEta");
    TH2D *hSubLeadJetUnmatchedPtVsEtaTrkMaxCut = (TH2D*)inFile->Get("hRecoSubLeadJetUnmatchedPtVsEta");

    //
    // Retrieve JetId selection histograms
    //
    TH2D *hJetPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveJetPtVsEtaJetIdCut");
    TH2D *hJetPtVsEtaMatchedJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEtaJetIdCut");
    TH2D *hJetPtVsEtaUnmatchedJetIdCut = (TH2D*)inFile->Get("hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut");
    TH2D *hLeadJetAllPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoLeadJetAllPtVsEtaJetIdCut");
    TH2D *hLeadJetMatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoLeadJetMatchedPtVsEtaJetIdCut");
    TH2D *hLeadJetUnmatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoLeadJetUnmatchedPtVsEtaJetIdCut");
    TH2D *hSubLeadJetAllPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoSubLeadJetAllPtVsEtaJetIdCut");
    TH2D *hSubLeadJetMatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoSubLeadJetMatchedPtVsEtaJetIdCut");
    TH2D *hSubLeadJetUnmatchedPtVsEtaJetIdCut = (TH2D*)inFile->Get("hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut");

    //
    // Rebin 2D histograms
    //
    hJetPtVsEtaKineCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaMatchedKineCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaUnmatchedKineCut->Rebin2D( rebinX, rebinY );

    hJetPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaMatchedTrkMaxCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaUnmatchedTrkMaxCut->Rebin2D( rebinX, rebinY );
    hLeadJetAllPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hLeadJetMatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hLeadJetUnmatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetAllPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetMatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetUnmatchedPtVsEtaTrkMaxCut->Rebin2D( rebinX, rebinY );

    hJetPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaMatchedJetIdCut->Rebin2D( rebinX, rebinY );
    hJetPtVsEtaUnmatchedJetIdCut->Rebin2D( rebinX, rebinY );
    hLeadJetAllPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hLeadJetMatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hLeadJetUnmatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetAllPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetMatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );
    hSubLeadJetUnmatchedPtVsEtaJetIdCut->Rebin2D( rebinX, rebinY );

    //
    // Make 1D eta distributions to compare
    //
    TH1D *hJetEtaKineCut = dynamic_cast<TH1D*>(hJetPtVsEtaKineCut->ProjectionX( Form("hJetEtaKineCut") ) );
    TH1D *hJetEtaTrkMaxCut = dynamic_cast<TH1D*>(hJetPtVsEtaTrkMaxCut->ProjectionX( Form("hJetEtaTrkMaxCut") ) );
    TH1D *hJetEtaJetIdCut = dynamic_cast<TH1D*>(hJetPtVsEtaJetIdCut->ProjectionX( Form("hJetEtaJetIdCut") ) );

    set1DStyle(hJetEtaKineCut, kineStyle, kTRUE);
    set1DStyle(hJetEtaTrkMaxCut, trkMaxStyle, kTRUE);
    set1DStyle(hJetEtaJetIdCut, jetIdStyle, kTRUE);

    //
    // Divide 2D histograms
    //
    hJetPtVsEtaMatchedKineCut->Divide(hJetPtVsEtaMatchedKineCut, hJetPtVsEtaKineCut, 1., 1., "b");
    hJetPtVsEtaUnmatchedKineCut->Divide(hJetPtVsEtaUnmatchedKineCut, hJetPtVsEtaKineCut, 1., 1., "b");

    hJetPtVsEtaMatchedTrkMaxCut->Divide(hJetPtVsEtaMatchedTrkMaxCut, hJetPtVsEtaTrkMaxCut, 1., 1., "b");
    hJetPtVsEtaUnmatchedTrkMaxCut->Divide(hJetPtVsEtaUnmatchedTrkMaxCut, hJetPtVsEtaTrkMaxCut, 1., 1., "b");
    hLeadJetMatchedPtVsEtaTrkMaxCut->Divide(hLeadJetMatchedPtVsEtaTrkMaxCut, hLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    hLeadJetUnmatchedPtVsEtaTrkMaxCut->Divide(hLeadJetUnmatchedPtVsEtaTrkMaxCut, hLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    hSubLeadJetMatchedPtVsEtaTrkMaxCut->Divide(hSubLeadJetMatchedPtVsEtaTrkMaxCut, hSubLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");
    hSubLeadJetUnmatchedPtVsEtaTrkMaxCut->Divide(hSubLeadJetUnmatchedPtVsEtaTrkMaxCut, hSubLeadJetAllPtVsEtaTrkMaxCut, 1., 1., "b");

    hJetPtVsEtaMatchedJetIdCut->Divide(hJetPtVsEtaMatchedJetIdCut, hJetPtVsEtaJetIdCut, 1., 1., "b");
    hJetPtVsEtaUnmatchedJetIdCut->Divide(hJetPtVsEtaUnmatchedJetIdCut, hJetPtVsEtaJetIdCut, 1., 1., "b");
    hLeadJetMatchedPtVsEtaJetIdCut->Divide(hLeadJetMatchedPtVsEtaJetIdCut, hLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    hLeadJetUnmatchedPtVsEtaJetIdCut->Divide(hLeadJetUnmatchedPtVsEtaJetIdCut, hLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    hSubLeadJetMatchedPtVsEtaJetIdCut->Divide(hSubLeadJetMatchedPtVsEtaJetIdCut, hSubLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");
    hSubLeadJetUnmatchedPtVsEtaJetIdCut->Divide(hSubLeadJetUnmatchedPtVsEtaJetIdCut, hSubLeadJetAllPtVsEtaJetIdCut, 1., 1., "b");

    //
    // Set style
    //
    set2DStyle(hJetPtVsEtaMatchedKineCut);
    set2DStyle(hJetPtVsEtaUnmatchedKineCut);

    set2DStyle(hJetPtVsEtaMatchedTrkMaxCut);
    set2DStyle(hJetPtVsEtaUnmatchedTrkMaxCut);
    set2DStyle(hLeadJetMatchedPtVsEtaTrkMaxCut);
    set2DStyle(hLeadJetUnmatchedPtVsEtaTrkMaxCut);
    set2DStyle(hSubLeadJetMatchedPtVsEtaTrkMaxCut);
    set2DStyle(hSubLeadJetUnmatchedPtVsEtaTrkMaxCut);

    set2DStyle(hJetPtVsEtaMatchedJetIdCut);
    set2DStyle(hJetPtVsEtaUnmatchedJetIdCut);
    set2DStyle(hLeadJetMatchedPtVsEtaJetIdCut);
    set2DStyle(hLeadJetUnmatchedPtVsEtaJetIdCut);
    set2DStyle(hSubLeadJetMatchedPtVsEtaJetIdCut);
    set2DStyle(hSubLeadJetUnmatchedPtVsEtaJetIdCut);

    //
    // Retrieve 1D binning
    //
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

    //
    // Reserve histograms for projections
    //

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

    Int_t sizeX{800};
    Int_t sizeY{1200};
    TCanvas *cEtaProj = new TCanvas("cEtaProj", "cEtaProj", sizeX, sizeY);
    cEtaProj->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtProj = new TCanvas("cPtProj", "cPtProj", sizeX, sizeY);
    cPtProj->Divide(5, ( (fEtaBins % 5) == 0 ) ? (fEtaBins / 5) : (fEtaBins / 5 + 1) );

    TCanvas *cEtaProjRat = new TCanvas("cEtaProjRat", "cEtaProjRat", sizeX, sizeY);
    cEtaProjRat->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtProjRat = new TCanvas("cPtProjRat", "cPtProjRat", sizeX, sizeY);
    cPtProjRat->Divide(5, ( (fEtaBins % 5) == 0 ) ? (fEtaBins / 5) : (fEtaBins / 5 + 1) );

    TCanvas *cEtaProjSep  = new TCanvas("cEtaProjSep", "cEtaProjSep", sizeX, sizeY);
    cEtaProjSep->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtProjSep = new TCanvas("cPtProjSep", "cPtProjSep", sizeX, sizeY);
    cPtProjSep->Divide(5, ( (fEtaBins % 5) == 0 ) ? (fEtaBins / 5) : (fEtaBins / 5 + 1) );

    // Create legends
    TLegend *hEtaLegend[fPtBins];
    TLegend *hEtaLegend2[fPtBins];
    TLegend *hPtLegend[fEtaBins];
    TLegend *hPtLegend2[fEtaBins];

    // Plot integrated eta distributions
    c->cd();
    setPadStyle();
    hJetEtaTrkMaxCut->Draw("same");
    if ( plotJetIdCut ) {
        hJetEtaJetIdCut->Draw("same");
    }
    if ( plotKineCut ) {
        hJetEtaKineCut->Draw("same");
    }
    
    t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
    t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
    hEtaLegend[0] = new TLegend(0.45, 0.35, 0.65, 0.55);
    hEtaLegend[0]->SetTextSize(0.04);
    hEtaLegend[0]->SetLineWidth(0);
    hEtaLegend[0]->AddEntry(hJetEtaTrkMaxCut, Form("TrkMax"), "p");
    if ( plotKineCut ) {
        hEtaLegend[0]->AddEntry(hJetEtaKineCut, Form("KineOnly"), "p");
    }
    if ( plotJetIdCut ) {
        hEtaLegend[0]->AddEntry(hJetEtaJetIdCut, Form("JetId"), "p");
    }
    hEtaLegend[0]->Draw();
    c->SaveAs( Form("%s/pPb8160_%s_eta_distributions_%s.pdf", 
                    date.Data(), direction.Data(), branchName.Data()) );

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
        hEtaLegend[i-1]->AddEntry(hJetMatchedEtaTrkMaxCut[i-1], Form("TrkMax matched"), "p");
        if ( plotKineCut ) {
            hEtaLegend[i-1]->AddEntry(hJetMatchedEtaKineCut[i-1], Form("KineOnly matched"), "p");
        }
        if ( plotJetIdCut ) {
            hEtaLegend[i-1]->AddEntry(hJetMatchedEtaJetIdCut[i-1], Form("JetId matched"), "p");
        }
        hEtaLegend[i-1]->AddEntry(hJetUnmatchedEtaTrkMaxCut[i-1], Form("TrkMax fakes"), "p");
        if ( plotKineCut ) {
            hEtaLegend[i-1]->AddEntry(hJetUnmatchedEtaKineCut[i-1], Form("KineOnly fakes"), "p");
        }
        if ( plotJetIdCut ) {
            hEtaLegend[i-1]->AddEntry(hJetUnmatchedEtaJetIdCut[i-1], Form("JetId fakes"), "p");
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
            set1DStyle(hJetMatchedEtaRatioJetIdCut[i-1], jetIdStyle);
            if ( plotKineCut ) {
                hJetMatchedEtaRatioKineCut[i-1]->Draw("same");
            }
            hJetUnmatchedEtaRatioJetIdCut[i-1]->Draw("same");
            set1DStyle(hJetUnmatchedEtaRatioJetIdCut[i-1], jetIdStyle+3);
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

        // With leading and subleading jets for trkMax selection
        c->cd();
        setPadStyle();
        hJetMatchedEtaTrkMaxCut[i-1]->Draw();
        set1DStyle(hJetMatchedEtaTrkMaxCut[i-1], trkMaxStyle);
        hLeadJetMatchedEtaTrkMaxCut[i-1]->Draw("same");
        set1DStyle(hLeadJetMatchedEtaTrkMaxCut[i-1], kineStyle);
        hSubLeadJetMatchedEtaTrkMaxCut[i-1]->Draw("same");
        set1DStyle(hSubLeadJetMatchedEtaTrkMaxCut[i-1], jetIdStyle);

        hJetUnmatchedEtaTrkMaxCut[i-1]->Draw("same");
        set1DStyle(hJetUnmatchedEtaTrkMaxCut[i-1], trkMaxStyle+3);
        hLeadJetUnmatchedEtaTrkMaxCut[i-1]->Draw("same");
        set1DStyle(hLeadJetUnmatchedEtaTrkMaxCut[i-1], kineStyle+3);
        hSubLeadJetUnmatchedEtaTrkMaxCut[i-1]->Draw("same");
        set1DStyle(hSubLeadJetUnmatchedEtaTrkMaxCut[i-1], jetIdStyle+3);

        hJetMatchedEtaTrkMaxCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);

        t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
        t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
        t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );

        hEtaLegend2[i-1] = new TLegend(0.65, 0.35, 0.75, 0.55);
        hEtaLegend2[i-1]->SetTextSize(0.04);
        hEtaLegend2[i-1]->SetLineWidth(0);
        hEtaLegend2[i-1]->AddEntry("type", "TrkMax sel", "");
        hEtaLegend2[i-1]->AddEntry(hJetMatchedEtaTrkMaxCut[i-1],        Form("Incl matched"), "p");
        hEtaLegend2[i-1]->AddEntry(hLeadJetMatchedEtaTrkMaxCut[i-1],    Form("Lead matched"), "p");
        hEtaLegend2[i-1]->AddEntry(hSubLeadJetMatchedEtaTrkMaxCut[i-1], Form("SubLead matched"), "p");
        hEtaLegend2[i-1]->AddEntry(hJetUnmatchedEtaTrkMaxCut[i-1],        Form("Incl fakes"), "p");
        hEtaLegend2[i-1]->AddEntry(hLeadJetUnmatchedEtaTrkMaxCut[i-1],    Form("Lead fakes"), "p");
        hEtaLegend2[i-1]->AddEntry(hSubLeadJetUnmatchedEtaTrkMaxCut[i-1], Form("SubLead fakes"), "p");

        hEtaLegend2[i-1]->Draw();
        c->SaveAs( Form("%s/pPb8160_%s_eta_fakes_lead_sub_trkMax_pt_%d_%d_%s.pdf", 
                        date.Data(), direction.Data(), 
                        (Int_t)(fPtRange[0] + (i-1) * ptStep),
                        (Int_t)(fPtRange[0] + i * ptStep),
                        branchName.Data()) );

        // With leading and subleading jets for jetId selection
        if ( plotJetIdCut ) {
            c->cd();
            setPadStyle();
            hJetMatchedEtaJetIdCut[i-1]->Draw();
            set1DStyle(hJetMatchedEtaJetIdCut[i-1], jetIdStyle);
            hLeadJetMatchedEtaJetIdCut[i-1]->Draw("same");
            set1DStyle(hLeadJetMatchedEtaJetIdCut[i-1], kineStyle);
            hSubLeadJetMatchedEtaJetIdCut[i-1]->Draw("same");
            set1DStyle(hSubLeadJetMatchedEtaJetIdCut[i-1], trkMaxStyle);

            hJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
            set1DStyle(hJetUnmatchedEtaJetIdCut[i-1], jetIdStyle+3);
            hLeadJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
            set1DStyle(hLeadJetUnmatchedEtaJetIdCut[i-1], kineStyle+3);
            hSubLeadJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
            set1DStyle(hSubLeadJetUnmatchedEtaJetIdCut[i-1], trkMaxStyle+3);

            hJetMatchedEtaJetIdCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);

            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                        fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );

            hEtaLegend2[i-1] = new TLegend(0.65, 0.35, 0.75, 0.55);
            hEtaLegend2[i-1]->SetTextSize(0.04);
            hEtaLegend2[i-1]->SetLineWidth(0);
            hEtaLegend2[i-1]->AddEntry("type", "JetId sel", "");
            hEtaLegend2[i-1]->AddEntry(hJetMatchedEtaJetIdCut[i-1],        Form("Incl matched"), "p");
            hEtaLegend2[i-1]->AddEntry(hLeadJetMatchedEtaJetIdCut[i-1],    Form("Lead matched"), "p");
            hEtaLegend2[i-1]->AddEntry(hSubLeadJetMatchedEtaJetIdCut[i-1], Form("SubLead matched"), "p");
            hEtaLegend2[i-1]->AddEntry(hJetUnmatchedEtaJetIdCut[i-1],        Form("Incl fakes"), "p");
            hEtaLegend2[i-1]->AddEntry(hLeadJetUnmatchedEtaJetIdCut[i-1],    Form("Lead fakes"), "p");
            hEtaLegend2[i-1]->AddEntry(hSubLeadJetUnmatchedEtaJetIdCut[i-1], Form("SubLead fakes"), "p");

            hEtaLegend2[i-1]->Draw();
            c->SaveAs( Form("%s/pPb8160_%s_eta_fakes_lead_sub_jetId_pt_%d_%d_%s.pdf", 
                            date.Data(), direction.Data(), 
                            (Int_t)(fPtRange[0] + (i-1) * ptStep),
                            (Int_t)(fPtRange[0] + i * ptStep),
                            branchName.Data()) );

            // Make for all projections
            cEtaProjSep->cd(i);
            setPadStyle();
            hJetMatchedEtaJetIdCut[i-1]->Draw();
            hLeadJetMatchedEtaJetIdCut[i-1]->Draw("same");
            hSubLeadJetMatchedEtaJetIdCut[i-1]->Draw("same");

            hJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
            hLeadJetUnmatchedEtaJetIdCut[i-1]->Draw("same");
            hSubLeadJetUnmatchedEtaJetIdCut[i-1]->Draw("same");

            hJetMatchedEtaJetIdCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);

            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%3.0f<p_{T} (GeV/c)<%3.0f", 
                           fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep ) );
            hEtaLegend2[i-1]->Draw();
        } // if ( plotJetIdCut )

    } // for (Int_t i{1}; i<fPtBins; i++)


    //
    // Plotting pt distributions
    //
    for (Int_t i{1}; i<=fEtaBins; i++) {

        // Make ratios for the kine selection to trkMax one
        if ( plotKineCut ) {
            hJetMatchedPtRatioKineCut[i-1] = dynamic_cast<TH1D*>( hJetMatchedPtKineCut[i-1]->Clone( Form("hJetMatchedPtRatioKineCut_%d", i-1) ) );
            make1DRatio(hJetMatchedPtRatioKineCut[i-1], hJetMatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", kineStyle);
            hJetUnmatchedPtRatioKineCut[i-1] = dynamic_cast<TH1D*>( hJetUnmatchedPtKineCut[i-1]->Clone( Form("hJetUnmatchedPtRatioKineCut_%d", i-1) ) );
            make1DRatio(hJetUnmatchedPtRatioKineCut[i-1], hJetUnmatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", kineStyle+3);
        }

        // Make ratios for the jetId selection to trkMax one
        if ( plotJetIdCut ) {
            hJetMatchedPtRatioJetIdCut[i-1] = dynamic_cast<TH1D*>( hJetMatchedPtJetIdCut[i-1]->Clone( Form("hJetMatchedPtRatioJetIdCut_%d", i-1) ) );
            make1DRatio(hJetMatchedPtRatioJetIdCut[i-1], hJetMatchedPtTrkMaxCut[i-1], "Ratio to TrkMax", jetIdStyle);
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

            // With leading and subleading jets for jetId selection
            c->cd();
            setPadStyle();
            hJetMatchedPtJetIdCut[i-1]->Draw();
            set1DStyle(hJetMatchedPtJetIdCut[i-1], jetIdStyle);
            hLeadJetMatchedPtJetIdCut[i-1]->Draw("same");
            set1DStyle(hLeadJetMatchedPtJetIdCut[i-1], kineStyle);
            hSubLeadJetMatchedPtJetIdCut[i-1]->Draw("same");
            set1DStyle(hSubLeadJetMatchedPtJetIdCut[i-1], trkMaxStyle);

            hJetUnmatchedPtJetIdCut[i-1]->Draw("same");
            set1DStyle(hJetUnmatchedPtJetIdCut[i-1], jetIdStyle+3);
            hLeadJetUnmatchedPtJetIdCut[i-1]->Draw("same");
            set1DStyle(hLeadJetUnmatchedPtJetIdCut[i-1], kineStyle+3);
            hSubLeadJetUnmatchedPtJetIdCut[i-1]->Draw("same");
            set1DStyle(hSubLeadJetUnmatchedPtJetIdCut[i-1], trkMaxStyle+3);

            hJetMatchedPtJetIdCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);

            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                           fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );

            hPtLegend2[i-1] = new TLegend(0.65, 0.35, 0.75, 0.55);
            hPtLegend2[i-1]->SetTextSize(0.04);
            hPtLegend2[i-1]->SetLineWidth(0);
            hPtLegend2[i-1]->AddEntry("type", "JetId sel", "");
            hPtLegend2[i-1]->AddEntry(hJetMatchedPtJetIdCut[i-1],        Form("Incl matched"), "p");
            hPtLegend2[i-1]->AddEntry(hLeadJetMatchedPtJetIdCut[i-1],    Form("Lead matched"), "p");
            hPtLegend2[i-1]->AddEntry(hSubLeadJetMatchedPtJetIdCut[i-1], Form("SubLead matched"), "p");
            hPtLegend2[i-1]->AddEntry(hJetUnmatchedPtJetIdCut[i-1],        Form("Incl fakes"), "p");
            hPtLegend2[i-1]->AddEntry(hLeadJetUnmatchedPtJetIdCut[i-1],    Form("Lead fakes"), "p");
            hPtLegend2[i-1]->AddEntry(hSubLeadJetUnmatchedPtJetIdCut[i-1], Form("SubLead fakes"), "p");

            hPtLegend2[i-1]->Draw();
            c->SaveAs( Form("%s/pPb8160_%s_pt_fakes_lead_sub_jetId_eta_%d_%d_%s.pdf", 
                            date.Data(), direction.Data(), 
                            (Int_t)( (fEtaRange[0] + (i-1) * etaStep) * 10),
                            (Int_t)( (fEtaRange[0] + i * etaStep) * 10),
                            branchName.Data()) );

            // Make for all projections
            cPtProjSep->cd(i);
            setPadStyle();
            hJetMatchedPtJetIdCut[i-1]->Draw();
            hLeadJetMatchedPtJetIdCut[i-1]->Draw("same");
            hSubLeadJetMatchedPtJetIdCut[i-1]->Draw("same");

            hJetUnmatchedPtJetIdCut[i-1]->Draw("same");
            hLeadJetUnmatchedPtJetIdCut[i-1]->Draw("same");
            hSubLeadJetUnmatchedPtJetIdCut[i-1]->Draw("same");

            hJetMatchedPtJetIdCut[i-1]->GetYaxis()->SetRangeUser(0., 1.05);

            t.DrawLatexNDC(0.35, 0.93, "PYTHIA8+EPOS" );
            t.DrawLatexNDC(0.75, 0.85, branchName.Data() );
            t.DrawLatexNDC(0.3, 0.2, Form("%2.1f<#eta<%2.1f", 
                           fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep ) );
            hPtLegend2[i-1]->Draw();
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
    if ( plotJetIdCut ) {
        cEtaProjSep->SaveAs( Form("%s/pPb8160_%s_eta_fake_lead_sub_jetId_projections_all_%s.pdf", 
                                  date.Data(), direction.Data(), branchName.Data()) );
        cPtProjSep->SaveAs( Form("%s/pPb8160_%s_pt_fake_lead_sub_jetId_projections_all_%s.pdf", 
                                  date.Data(), direction.Data(), branchName.Data()) );
    }

}

//________________
bool isGoodFile(TFile *f) {
    bool isGood{true};
    if ( !f ) {
        std::cout << Form("[ERROR] Input file %s does not exist\n", f->GetName() );
        isGood = {false};
    }
    if ( f->IsZombie() ) {
        std::cout << Form("[ERROR] File %s is zombie\n", f->GetName() );
        isGood = {false};
    }
    return isGood;
}

//________________
void print2DarrayFromHisto(TH2* h) {
    Int_t nBinsX = h->GetNbinsX();
    Int_t nBinsY = h->GetNbinsY();
    std::cout << Form("Correction factor for %s \n", h->GetName() );
    std::cout << Form("Double_t nCorr[%d][%d] = {\n", nBinsX, nBinsY);
    for (Int_t i=1; i<=nBinsX; i++) {
        std::cout << "\t{";
        for (Int_t j=1; j<=nBinsY; j++) {

            Double_t val = h->GetBinContent(i, j);
            if ( val == 0 ) val = {1.};
            if ( j != nBinsY ) {
                std::cout << val << ", ";
            }
            else {
                std::cout << val << " ";
            }
        } // for (Int_t j=1; j<=nBinsY; j++)
        if ( i != nBinsX ) {
            std::cout << "},\n";
        }
        else {
            std::cout << "}\n";
        }
    } // for (Int_t i=1; i<=nBinsX; i++)
    std::cout << "}\n";
}

//________________
void makeMcWeightingMaps(TFile *embFile, TFile *mbFile, TFile *jet60File,
                         TFile *jet80File, TFile *jet100File, TString date) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    Int_t rebin = 2;

    Bool_t doReweighting{kFALSE};
    TString embFileName = embFile->GetName();
    if ( embFileName.Contains("weight") ) {
        doReweighting = {kTRUE};
    }

    // Retrieve distributions for embedding
    TH2D *hEmbRecoPtLeadPtSublead = dynamic_cast< TH2D* >( embFile->Get("hRecoPtLeadPtSublead") );
    hEmbRecoPtLeadPtSublead->SetName("hEmbRecoPtLeadPtSublead");
    hEmbRecoPtLeadPtSublead->Rebin2D( rebin, rebin );
    hEmbRecoPtLeadPtSublead->Scale( 1./hEmbRecoPtLeadPtSublead->Integral() );
    set2DStyle(hEmbRecoPtLeadPtSublead);

    
    TH2D* hEmbRecoPtLeadPtSubleadMcReweight{nullptr};
    if ( doReweighting ) {
        hEmbRecoPtLeadPtSubleadMcReweight = dynamic_cast< TH2D* >( embFile->Get("hRecoPtLeadPtSubleadMcReweight") );
        hEmbRecoPtLeadPtSubleadMcReweight->SetName("hEmbRecoPtLeadPtSubleadMcReweight");
        hEmbRecoPtLeadPtSubleadMcReweight->Rebin2D( rebin, rebin );
        hEmbRecoPtLeadPtSubleadMcReweight->Scale( 1./hEmbRecoPtLeadPtSubleadMcReweight->Integral() );
        set2DStyle(hEmbRecoPtLeadPtSubleadMcReweight);
    }
    
    // Retrieve distributions for MB
    TH2D* hMBRecoPtLeadPtSublead = dynamic_cast< TH2D* > ( mbFile->Get("hRecoPtLeadPtSublead") );
    hMBRecoPtLeadPtSublead->SetName("hMBRecoPtLeadPtSublead");
    hMBRecoPtLeadPtSublead->Rebin2D( rebin, rebin );
    hMBRecoPtLeadPtSublead->Scale( 1./hMBRecoPtLeadPtSublead->Integral() );
    set2DStyle(hMBRecoPtLeadPtSublead);

    TH2D* hMBRecoPtLeadPtSubleadMcReweight{nullptr};
    if ( doReweighting ) {
        hMBRecoPtLeadPtSubleadMcReweight = dynamic_cast< TH2D* > ( mbFile->Get("hRecoPtLeadPtSubleadMcReweight") );
        hMBRecoPtLeadPtSubleadMcReweight->SetName("hMBRecoPtLeadPtSubleadMcReweight");
        hMBRecoPtLeadPtSubleadMcReweight->Rebin2D( rebin, rebin );
        hMBRecoPtLeadPtSubleadMcReweight->Scale( 1./hMBRecoPtLeadPtSubleadMcReweight->Integral() );
        set2DStyle(hMBRecoPtLeadPtSubleadMcReweight);
    }

    
    // Retrieve distributions for Jet60
    TH2D* hJet60RecoPtLeadPtSublead = dynamic_cast< TH2D* > ( jet60File->Get("hRecoPtLeadPtSublead") );
    hJet60RecoPtLeadPtSublead->SetName("hJet60RecoPtLeadPtSublead");
    hJet60RecoPtLeadPtSublead->Rebin2D( rebin, rebin );
    hJet60RecoPtLeadPtSublead->Scale( 1./hJet60RecoPtLeadPtSublead->Integral() );
    set2DStyle(hJet60RecoPtLeadPtSublead);

    TH2D* hJet60RecoPtLeadPtSubleadMcReweight{nullptr};
    if ( doReweighting ) {
        hJet60RecoPtLeadPtSubleadMcReweight = dynamic_cast< TH2D* > ( jet60File->Get("hRecoPtLeadPtSubleadMcReweight") );
        hJet60RecoPtLeadPtSubleadMcReweight->SetName("hJet60RecoPtLeadPtSubleadMcReweight");
        hJet60RecoPtLeadPtSubleadMcReweight->Rebin2D( rebin, rebin );
        hJet60RecoPtLeadPtSubleadMcReweight->Scale( 1./hJet60RecoPtLeadPtSubleadMcReweight->Integral() );
        set2DStyle(hJet60RecoPtLeadPtSubleadMcReweight);
    }
    
    // Retrieve distributions for Jet80
    TH2D* hJet80RecoPtLeadPtSublead = dynamic_cast< TH2D* > ( jet80File->Get("hRecoPtLeadPtSublead") );
    hJet80RecoPtLeadPtSublead->SetName("hJet80RecoPtLeadPtSublead");
    hJet80RecoPtLeadPtSublead->Rebin2D( rebin, rebin );
    hJet80RecoPtLeadPtSublead->Scale( 1./hJet80RecoPtLeadPtSublead->Integral() );
    set2DStyle(hJet80RecoPtLeadPtSublead);

    TH2D* hJet80RecoPtLeadPtSubleadMcReweight{nullptr};
    if ( doReweighting ) {
        hJet80RecoPtLeadPtSubleadMcReweight = dynamic_cast< TH2D* > ( jet80File->Get("hRecoPtLeadPtSubleadMcReweight") );
        hJet80RecoPtLeadPtSubleadMcReweight->SetName("hJet80RecoPtLeadPtSublead");
        hJet80RecoPtLeadPtSubleadMcReweight->Rebin2D( rebin, rebin );
        hJet80RecoPtLeadPtSubleadMcReweight->Scale( 1./hJet80RecoPtLeadPtSubleadMcReweight->Integral() );
        set2DStyle(hJet80RecoPtLeadPtSubleadMcReweight);
    }
    
    // Retrieve distributions for Jet100
    TH2D* hJet100RecoPtLeadPtSublead = dynamic_cast< TH2D* > ( jet100File->Get("hRecoPtLeadPtSublead") );
    hJet100RecoPtLeadPtSublead->SetName("hJet100RecoPtLeadPtSublead");
    hJet100RecoPtLeadPtSublead->Rebin2D( rebin, rebin );
    hJet100RecoPtLeadPtSublead->Scale( 1./hJet100RecoPtLeadPtSublead->Integral() );
    set2DStyle(hJet100RecoPtLeadPtSublead);

    TH2D* hJet100RecoPtLeadPtSubleadMcReweight{nullptr};
    if ( doReweighting ) {
        hJet100RecoPtLeadPtSubleadMcReweight = dynamic_cast< TH2D* > ( jet100File->Get("hRecoPtLeadPtSubleadMcReweight") );
        hJet100RecoPtLeadPtSubleadMcReweight->SetName("hJet100RecoPtLeadPtSubleadMcReweight");
        hJet100RecoPtLeadPtSubleadMcReweight->Rebin2D( rebin, rebin );
        hJet100RecoPtLeadPtSubleadMcReweight->Scale( 1./hJet100RecoPtLeadPtSubleadMcReweight->Integral() );
        set2DStyle(hJet100RecoPtLeadPtSubleadMcReweight);
    }


    // Vz distributions
    TH1D* hEmbVzWeighted = dynamic_cast< TH1D* >( embFile->Get("hVzWeighted") );
    hEmbVzWeighted->SetName("hEmbVzWeighted");
    TH1D* hMBVzWeighted = dynamic_cast< TH1D* >( mbFile->Get("hVzWeighted") );
    hMBVzWeighted->SetName("hMBVzWeighted");
    TH1D* hJet60VzWeighted = dynamic_cast< TH1D* >( jet60File->Get("hVzWeighted") );
    hJet60VzWeighted->SetName("hJet60VzWeighted");
    TH1D* hJet80VzWeighted = dynamic_cast< TH1D* >( jet80File->Get("hVzWeighted") );
    hJet80VzWeighted->SetName("hJet80VzWeighted");
    TH1D* hJet100VzWeighted = dynamic_cast< TH1D* >( jet100File->Get("hVzWeighted") );
    hJet100VzWeighted->SetName("hJet100VzWeighted");

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};
    Int_t refSelType{6};
    Int_t dataType{2};

    set1DStyle(hEmbVzWeighted, genType);
    set1DStyle(hMBVzWeighted, recoType);
    set1DStyle(hJet60VzWeighted, refType);
    set1DStyle(hJet80VzWeighted, refSelType);
    set1DStyle(hJet100VzWeighted, dataType);

    Double_t nEmbNevents = hEmbVzWeighted->Integral();       hEmbVzWeighted->Scale(1./nEmbNevents);
    Double_t nMBNevents = hMBVzWeighted->Integral();         hMBVzWeighted->Scale(1./nMBNevents);
    Double_t nJet60Nevents = hJet60VzWeighted->Integral();   hJet60VzWeighted->Scale(1./nJet60Nevents);
    Double_t nJet80Nevents = hJet80VzWeighted->Integral();   hJet80VzWeighted->Scale(1./nJet80Nevents);
    Double_t nJet100Nevents = hJet100VzWeighted->Integral(); hJet100VzWeighted->Scale(1./nJet100Nevents);

    TCanvas *cVz = new TCanvas("cVz","cVz", 1200, 800);
    cVz->cd();
    setPadStyle();
    hJet100VzWeighted->Draw();
    hJet80VzWeighted->Draw("same");
    hJet60VzWeighted->Draw("same");
    hMBVzWeighted->Draw("same");
    hEmbVzWeighted->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.65, 0.4, 0.85);
    leg->SetTextSize(0.04);
    leg->SetLineWidth(0);
    leg->AddEntry(hEmbVzWeighted, "Embedding", "p");
    leg->AddEntry(hMBVzWeighted, "Min. bias.", "p");
    leg->AddEntry(hJet60VzWeighted, "Jet60", "p");
    leg->AddEntry(hJet80VzWeighted, "Jet80", "p");
    leg->AddEntry(hJet100VzWeighted, "Jet100", "p");
    leg->Draw();

    t.DrawLatexNDC(0.65, 0.85, Form("N_{events}(Emb): %.5f", nEmbNevents) );
    t.DrawLatexNDC(0.65, 0.8, Form("N_{events}(MB): %.f", nMBNevents) );
    t.DrawLatexNDC(0.65, 0.75, Form("N_{events}(Jet60): %.f", nJet60Nevents) );
    t.DrawLatexNDC(0.65, 0.7, Form("N_{events}(Jet80): %.f", nJet80Nevents) );
    t.DrawLatexNDC(0.65, 0.65, Form("N_{events}(Jet100): %.f", nJet100Nevents) );

    // Divide data / MC -> correction factor (multiplication) for MC
    hJet100RecoPtLeadPtSublead->Divide(hJet100RecoPtLeadPtSublead, hEmbRecoPtLeadPtSublead, 1., 1.);
    hJet80RecoPtLeadPtSublead->Divide(hJet80RecoPtLeadPtSublead, hEmbRecoPtLeadPtSublead, 1., 1.);
    hJet60RecoPtLeadPtSublead->Divide(hJet60RecoPtLeadPtSublead, hEmbRecoPtLeadPtSublead, 1., 1.);
    hMBRecoPtLeadPtSublead->Divide(hMBRecoPtLeadPtSublead, hEmbRecoPtLeadPtSublead, 1., 1.);

    //print2DarrayFromHisto(hMBRecoPtLeadPtSublead);
    // print2DarrayFromHisto(hJet60RecoPtLeadPtSublead);
    // print2DarrayFromHisto(hJet80RecoPtLeadPtSublead);
    print2DarrayFromHisto(hJet100RecoPtLeadPtSublead);

    // Before reweighting
    TCanvas *cMaps = new TCanvas("cMaps", "cMaps", 1500, 300);
    cMaps->Divide(5, 1);

    cMaps->cd(1);
    setPadStyle();
    hMBRecoPtLeadPtSublead->Draw("colz");
    t.DrawLatexNDC(0.15, 0.93, Form("Correction factor min. bias") );

    cMaps->cd(2);
    setPadStyle();
    hJet60RecoPtLeadPtSublead->Draw("colz");
    t.DrawLatexNDC(0.15, 0.93, Form("Correction factor jet60") );

    cMaps->cd(3);
    setPadStyle();
    hJet80RecoPtLeadPtSublead->Draw("colz");
    t.DrawLatexNDC(0.15, 0.93, Form("Correction factor jet80") );

    cMaps->cd(4);
    setPadStyle();
    hJet100RecoPtLeadPtSublead->Draw("colz");
    t.DrawLatexNDC(0.15, 0.93, Form("Correction factor jet100") );

    cMaps->cd(5);
    setPadStyle();
    hEmbRecoPtLeadPtSublead->Draw("colz");
    t.DrawLatexNDC(0.15, 0.93, Form("Pure embedding") );

    // After reweighting
    TCanvas *cMaps2 = new TCanvas("cMaps2", "cMaps2", 1500, 300);
    if ( doReweighting ) {
        cMaps2->Divide(5, 1);

        cMaps->cd(1);
        setPadStyle();
        hMBRecoPtLeadPtSublead->Draw("colz");
        t.DrawLatexNDC(0.15, 0.93, Form("Correction factor min. bias") );

        cMaps2->cd(2);
        setPadStyle();
        hJet60RecoPtLeadPtSublead->Draw("colz");
        t.DrawLatexNDC(0.15, 0.93, Form("Correction factor jet60") );

        cMaps2->cd(3);
        setPadStyle();
        hJet80RecoPtLeadPtSublead->Draw("colz");
        t.DrawLatexNDC(0.15, 0.93, Form("Correction factor jet80") );

        cMaps2->cd(4);
        setPadStyle();
        hJet100RecoPtLeadPtSublead->Draw("colz");
        t.DrawLatexNDC(0.15, 0.93, Form("Correction factor jet100") );

        cMaps2->cd(5);
        setPadStyle();
        hEmbRecoPtLeadPtSublead->Draw("colz");
        t.DrawLatexNDC(0.15, 0.93, Form("Pure embedding") );
    }

    cVz->SaveAs( Form("%s/pPb8160_vz_EmbMbJet60Jet80Jet100.pdf", date.Data() ) );
    cMaps->SaveAs( Form("%s/pPb8160_ptLeadPtSublead_ReweightingMaps.pdf", date.Data() ) );
    if ( doReweighting ) {
        cMaps2->SaveAs( Form("%s/pPb8160_ptLeadPtSublead_ReweightingMaps_afterCorr.pdf", date.Data() ) );
    }
}


//________________
void pPb_qa() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // File names
    const Char_t *embeddingFileName = "/Users/gnigmat/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_ak4.root";
    const Char_t *mbFileName = "/Users/gnigmat/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4.root";
    const Char_t *jet60FileName = "/Users/gnigmat/cernbox/ana/pPb8160/exp/Jet60_pPb8160_ak4.root";
    const Char_t *jet80FileName = "/Users/gnigmat/cernbox/ana/pPb8160/exp/Jet80_pPb8160_ak4.root";
    const Char_t *jet100FileName = "/Users/gnigmat/cernbox/ana/pPb8160/exp/Jet100_pPb8160_ak4.root";

    // Files
    TFile *embFile = TFile::Open(embeddingFileName);
    TFile *mbFile = TFile::Open(mbFileName);
    TFile *jet60File = TFile::Open(jet60FileName);
    TFile *jet80File = TFile::Open(jet80FileName);
    TFile *jet100File = TFile::Open(jet100FileName);

    // Check files are good
    if ( !isGoodFile( embFile ) ) return;
    if ( !isGoodFile( mbFile ) ) return;
    if ( !isGoodFile( jet60File ) ) return;
    if ( !isGoodFile( jet80File ) ) return;
    if ( !isGoodFile( jet100File ) ) return;

    // Date extraction
    TDatime dt; 
    TString date { Form( "%d",dt.GetDate() ) };

    if ( !directoryExists( date.Data() ) ) {
        createDirectory( date.Data() );
    } 


    // Int_t branchId{0};
    // if ( inputFileName.Contains("akCs4") ) {
    //     branchId = {0};
    // }
    // else {
    //     branchId = {1};
    // }

    // Plot ptHat distribution
    //plotPtHat(embFile, date);

    // Compare inclusive reco, ref and gen transverse momentum spectra
    //compareInclusiveJetPtSpectra(embFile, date);

    // Plot jet reconstruction efficiency as a function of acceptance (pT vs eta)
    // plotEfficiency(embFile, date);

    // Plot dijet distributions
    // plotDijetDistributions(embFile, mbFile, jet60File, jet80File, jet100File, date);

    // Plot reco, reco with matching and calculate fakes
    // plotRecoAndFakes(embFile, date);

    // Plot correlation between ref and reco dijet eta
    // plotEtaDijetCorrelation(embFile, date);

    // Plot distributions for jetId
    // plotJetIdHistos(embFile, date);

    // Plot JES and JER
    // plotJESandJER(embFile, date);

    // Plot various correlation matrices
    // plotDijetResponseMatrices(embFile, mbFile, jet60File, jet80File, jet100File, date);

    // Build reweighting matrices for the MC w.r.t. data
    makeMcWeightingMaps(embFile, mbFile, jet60File, jet80File, jet100File, date);
}