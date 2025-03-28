// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TRatioPlot.h"
#include "TPad.h"

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
void plotCMSHeader(int collSystem = 0, double energy = 5.02) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, Form("%s #sqrt{s_{NN}} = %3.2f TeV", collSystemStr.Data(), energy) );
    t.SetTextSize(0.05);
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
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
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
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
}

//________________
void rescaleEta(TH1* h) {
    for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++) {
        Double_t val = h->GetBinContent( iBin );
        Double_t valErr = h->GetBinError( iBin );
        Double_t binWidth = h->GetBinWidth( iBin );
        h->SetBinContent( iBin, val / binWidth );
        h->SetBinError( iBin, valErr / binWidth );
    }
    h->Scale( 1. / h->Integral() );
}

//________________
void drawJecComparison(TCanvas *c, TH1D* h1, TH1D* h2, 
    double ptLow = 20., double ptHigh = 1500.,
    double ptHatLow=15.,
    const char* h1Name="Reco", const char* h2Name="Gen",
    int collSystem = 0, double energy = 5.02,
    bool isCM = false,
    bool isJet = true,
    bool isPt = false,
    double etaLow = -1.6, double etaHi = 1.6) {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // Energy in TeV

    // Text
    TLatex t;
    TString frameT = (isCM) ? "CM" : "Lab";
    TString jetType = (isJet) ? "jet" : "dijet";
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    int maximumBin = h1->GetMaximumBin();
    double maximumVal = h1->GetBinContent(maximumBin);

    // Number for plotting position
    double xRange[2] = {-3.2, 3.2};
    if (isPt)
    {
        xRange[0] = 0;
        xRange[1] = 500;
    }
    double yRange[2] = {0.0000001, maximumVal * 1.25};
    double legX[2] = {0.5, 0.7};
    double legY[2] = {0.2, 0.35};
    double ratioYRange[2] = {0.8, 1.2};
    if (isPt)
    {
        ratioYRange[0] = 0.9;
        ratioYRange[1] = 1.3;
    }

    if (!isPt)
    {
        h1->GetXaxis()->SetTitle(Form("#eta^{%s}", jetType.Data()));
        h1->GetYaxis()->SetTitle(Form("1/N_{%s} dN/d#eta^{%s}", jetType.Data(), jetType.Data()));
        h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    }
    else
    {
        h1->GetXaxis()->SetTitle(Form("p_{T}^{%s} (GeV)", jetType.Data()));
        h1->GetYaxis()->SetTitle(Form("1/N_{%s} dN/dp_{T}^{%s}", jetType.Data(), jetType.Data()));
    }

    std::vector<double> gridLineValues{0.95, 1.0, 1.05};

    // Create ratio plot
    TRatioPlot *ratioPlot = new TRatioPlot(h1, h2);
    ratioPlot->SetH1DrawOpt("E");
    ratioPlot->SetH2DrawOpt("E");
    ratioPlot->SetGridlines(gridLineValues);
    ratioPlot->Draw();

    // Set pad parameters
    ratioPlot->GetUpperPad()->SetFrameLineWidth(2);
    ratioPlot->GetLowerPad()->SetFrameLineWidth(2);
    ratioPlot->SetLeftMargin(0.15);
    ratioPlot->SetRightMargin(0.05);
    ratioPlot->SetSeparationMargin(0.03);

    // Lower plot style and titles
    ratioPlot->GetLowerRefXaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefXaxis()->SetTitleOffset(0.85);
    ratioPlot->GetLowerRefYaxis()->SetTitle(Form("%s / %s", h1Name, h2Name));
    ratioPlot->GetLowerRefYaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefYaxis()->SetTitleOffset(1.2);

    ratioPlot->GetLowerRefGraph()->SetMarkerStyle(20);
    ratioPlot->GetLowerRefGraph()->SetMarkerSize(1.2);
    ratioPlot->GetLowerRefGraph()->SetMarkerColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineWidth(2);
    ratioPlot->GetLowYaxis()->SetNdivisions(205);

    ratioPlot->GetLowerRefGraph()->SetMinimum(ratioYRange[0]);
    ratioPlot->GetLowerRefGraph()->SetMaximum(ratioYRange[1]);
    ratioPlot->GetLowerRefXaxis()->SetRangeUser(xRange[0], xRange[1]);

    ratioPlot->GetUpperPad()->cd();
    t.DrawLatexNDC(0.17, 0.84, Form("%s frame", frameT.Data()));
    plotCMSHeader(collSystem, energy);
    if (!isPt) {
        t.DrawLatexNDC(0.4, 0.84, Form("%4.0f < p_{T}^{%s} < %4.0f GeV", ptLow, jetType.Data(), ptHigh));
        t.DrawLatexNDC(0.75, 0.84, Form("#hat{p}_{T} > %3.0f GeV", ptHatLow));
    }
    else{
        t.DrawLatexNDC(0.4, 0.84, Form("%2.1f< #eta^{jet} < %2.1f", etaLow, etaHi));
        t.DrawLatexNDC(0.75, 0.84, Form("#hat{p}_{T} > %3.0f GeV", ptHatLow));
    }
    if (isPt)
    {
        ratioPlot->GetUpperPad()->SetLogy();
    }

    TLegend *leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h1, h1Name, "p");
    leg->AddEntry(h2, h2Name, "p");
    leg->Draw();
}

//________________
void drawEtaPtComparison(TCanvas *c, TH1D* h1, TH1D* h2, 
    int ptLow=15, int ptHigh = 45,
    double etaLow = -1.6, double etaHigh = 1.6,
    int ptHatLow=15,
    const char* h1Name="Reco (embed)", const char* h2Name="Reco (pythia)",
    int collSystem = 1, double energy = 8.16,
    bool isCM = false,
    bool isJet = true,
    bool isPt = false) {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // Energy in TeV

    // Text
    TLatex t;
    TString frameT = (isCM) ? "CM" : "Lab";
    TString jetType = (isJet) ? "jet" : "dijet";
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    int maximumBin = h1->GetMaximumBin();
    double maximumVal = h1->GetBinContent(maximumBin);

    // Number for plotting position
    double xRange[2] = {-3.2, 3.2};
    if (isPt) {
        xRange[0] = 0;
        xRange[1] = 500;
    }
    double yRange[2] = {0.0000001, maximumVal * 1.25};
    double legX[2] = {0.5, 0.7};
    double legY[2] = {0.2, 0.35};
    double ratioYRange[2] = {0.8, 1.2};
    if (isPt) {
        ratioYRange[0] = 0.9;
        ratioYRange[1] = 1.3;
    }

    if (!isPt) {
        h1->GetXaxis()->SetTitle(Form("#eta^{%s}", jetType.Data()));
        h1->GetYaxis()->SetTitle(Form("1/N_{%s} dN/d#eta^{%s}", jetType.Data(), jetType.Data()));
        h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    }
    else {
        h1->GetXaxis()->SetTitle(Form("p_{T}^{%s} (GeV)", jetType.Data()));
        h1->GetYaxis()->SetTitle(Form("1/N_{%s} dN/dp_{T}^{%s}", jetType.Data(), jetType.Data()));
    }

    std::vector<double> gridLineValues{0.95, 1.0, 1.05};

    // Create ratio plot
    TRatioPlot *ratioPlot = new TRatioPlot(h1, h2);
    ratioPlot->SetH1DrawOpt("E");
    ratioPlot->SetH2DrawOpt("E");
    ratioPlot->SetGridlines(gridLineValues);
    ratioPlot->Draw();

    // Set pad parameters
    ratioPlot->GetUpperPad()->SetFrameLineWidth(2);
    ratioPlot->GetLowerPad()->SetFrameLineWidth(2);
    ratioPlot->SetLeftMargin(0.15);
    ratioPlot->SetRightMargin(0.05);
    ratioPlot->SetSeparationMargin(0.03);

    // Lower plot style and titles
    ratioPlot->GetLowerRefXaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefXaxis()->SetTitleOffset(0.85);
    ratioPlot->GetLowerRefYaxis()->SetTitle(Form("%s / %s", h1Name, h2Name));
    ratioPlot->GetLowerRefYaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefYaxis()->SetTitleOffset(1.2);

    ratioPlot->GetLowerRefGraph()->SetMarkerStyle(20);
    ratioPlot->GetLowerRefGraph()->SetMarkerSize(1.2);
    ratioPlot->GetLowerRefGraph()->SetMarkerColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineWidth(2);
    ratioPlot->GetLowYaxis()->SetNdivisions(205);

    ratioPlot->GetLowerRefGraph()->SetMinimum(ratioYRange[0]);
    ratioPlot->GetLowerRefGraph()->SetMaximum(ratioYRange[1]);
    ratioPlot->GetLowerRefXaxis()->SetRangeUser(xRange[0], xRange[1]);

    ratioPlot->GetUpperPad()->cd();

    t.DrawLatexNDC(0.17, 0.84, Form("%s frame", frameT.Data()));
    plotCMSHeader(collSystem, energy);
    if (!isPt) {
        if (isJet) {
            t.DrawLatexNDC(0.4,  0.84, Form("%d < p_{T}^{%s} < %d GeV", ptLow, jetType.Data(), ptHigh));
            t.DrawLatexNDC(0.75, 0.84, Form("#hat{p}_{T} > %d GeV", ptHatLow));
        }
        else {
            t.DrawLatexNDC(0.4, 0.84, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHigh));
        }
    }
    else {
        if (isJet) {
            t.DrawLatexNDC(0.4,  0.84, Form("%2.1f< #eta^{%s} < %2.1f", etaLow, jetType.Data(), etaHigh));
            t.DrawLatexNDC(0.75, 0.84, Form("#hat{p}_{T} > %d GeV", ptHatLow));
        }
        else {
            t.DrawLatexNDC(0.4, 0.84, Form("%2.1f< #eta^{dijet} < %2.1f", etaLow, etaHigh));
        }
    }

    if (isPt){
        ratioPlot->GetUpperPad()->SetLogy();
    }

    TLegend *leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h1, h1Name, "p");
    leg->AddEntry(h2, h2Name, "p");
    leg->Draw();
}

//________________
void rescaleForwardBackward(TH1D *hForward, TH1D *hBackward) {
    // Rescale forward and backward histograms
    // to have the same integral
    double intForward = hForward->Integral();
    double intBackward = hBackward->Integral();
    double scaleFactor = intForward + intBackward;
    hForward->Scale( intForward / scaleFactor );
    hBackward->Scale( intBackward / scaleFactor );
}

//________________
// Compare reco, ref, refSel to gen inclusive jet eta distributions
void drawSingleJetToGenComparison(TCanvas *c, TH1D *hReco, TH1D *hRef = nullptr, 
                                  TH1D *hGen = nullptr, TH1D *hRefSel = nullptr,
                                  TH1D *hData = nullptr,
                                  int ptLow = 50, int ptHi = 60,
                                  int ptHatLow = 15,
                                  bool isCM = false, bool isFB = false,
                                  int collisionSystem = 1, double energy = 8.16) {
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.04 );
    TLegend *leg;

    double xRange[2] = { -3.2, 3.2 };
    double yRange[2] = {0.0000001, 0.08 };

    if ( isCM ) {
        if ( !isFB ) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] = 0.0000001; yRange[1] = 0.08;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.5;
            yRange[0] = 0.8; yRange[1] = 1.2;
        }
    }

    c->cd();
    setPadStyle();
    hReco->Draw();
    hGen->Draw("same");
    if ( hRef ) {
        hRef->Draw("same");
    }
    if ( hData ) {
        hData->Draw("same");
    }
    if ( hRefSel ) {
        hRefSel->Draw("same");
    }
    hReco->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hReco->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    hReco->GetYaxis()->SetTitle("1/N_{jet} dN/d#eta_{jet}");
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);        
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T} (GeV) < %d", ptLow, ptHi) );
    if ( isFB ) {
        leg = new TLegend(0.2, 0.25, 0.4, 0.4);
    }
    else {
        leg = new TLegend(0.2, 0.55, 0.4, 0.7);
    }
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    if ( hData ) {
        leg->AddEntry( hData, "Data", "p" );
    }
    leg->AddEntry( hReco, "Reco", "p" );
    leg->AddEntry( hGen, "Gen", "p" );
    if ( hRef ) {
        leg->AddEntry( hRef, "Ref", "p" );        
    }
    if ( hRefSel ) {
        leg->AddEntry( hRefSel, "RefSel", "p" );
    }
    leg->Draw();

}

//________________
// Plot reco, ref, refSel to gen inclusive jet eta ratio distributions
void drawSingleJetToGenRatio(TCanvas *c, TH1D *hReco2Gen, TH1D *hRef2Gen = nullptr, TH1D *hRefSel2Gen = nullptr,
                             TH1D *hData2Gen = nullptr,
                             int ptLow = 50, int ptHi = 60,
                             int ptHatLow = 15,
                             bool isCM = false, bool isFB = false,
                             int collisionSystem = 1, double energy = 8.16) {
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);
    TLegend *leg;

    double xRange[2] = {-3.2, 3.2};
    double yRange[2] = {0.75, 1.25};

    if (isCM) {
        if (!isFB) {
            xRange[0] = -3.2;
            xRange[1] = 3.2;
            yRange[0] = 0.75;
            yRange[1] = 1.25;
        }
        else {
            xRange[0] = 0;
            xRange[1] = 2.5;
            yRange[0] = 0.8;
            yRange[1] = 1.2;
        }
    }

    c->cd();
    setPadStyle();
    hReco2Gen->Draw();
    if (hRef2Gen) {
        hRef2Gen->Draw("same");
    }
    if (hData2Gen) {
        hData2Gen->Draw("same");
    }
    if (hRefSel2Gen) {
        hRefSel2Gen->Draw("same");
    }
    hReco2Gen->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hReco2Gen->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    hReco2Gen->GetYaxis()->SetTitle("Ratio to gen");
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T} (GeV) < %d", ptLow, ptHi));
    if (isFB) {
        leg = new TLegend(0.35, 0.25, 0.4, 0.4);
    }
    else {
        leg = new TLegend(0.35, 0.25, 0.65, 0.4);
    }
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->AddEntry(hReco2Gen, "Reco / Gen", "p");
    if (hData2Gen) {
        leg->AddEntry(hData2Gen, "Data / Gen", "p");
    }
    if (hRef2Gen) {
        leg->AddEntry(hRef2Gen, "Ref / Gen", "p");
    }
    if (hRefSel2Gen) {
        leg->AddEntry(hRefSel2Gen, "RefSel / Gen", "p");
    }
    leg->Draw();
}



//________________
// Draw comparison of reco, ref, gen and refSel dijet eta distributions
void drawDijetToGenComparison(TCanvas *c, TH1D *hReco, TH1D *hRef = nullptr, 
                         TH1D *hGen = nullptr, TH1D *hRefSel = nullptr,
                         TH1D *hData = nullptr,
                         int ptLow = 50, int ptHi = 60, 
                         bool isCM = false, bool isFB = false,
                         int collisionSystem = 1, double energy = 8.16) {

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );
    TLegend *leg;

    double xRange[2] = { -3.2, 3.2 };
    double yRange[2] = {0.0000001, 0.08 };

    if ( isCM ) {
        if ( !isFB ) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] = 0.0000001; yRange[1] = 0.08;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.5;
            yRange[0] = 0.8; yRange[1] = 1.2;
        }
    }

    c->cd();
    setPadStyle();
    hReco->Draw();
    hGen->Draw("same");
    if ( hRef ) {
        hRef->Draw("same");
    }
    if ( hData ) {
        hData->Draw("same");
    }
    if ( hRefSel ) {
        hRefSel->Draw("same");
    }
    hReco->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hReco->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);        
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi) );
    if ( isFB ) {
        leg = new TLegend(0.2, 0.25, 0.4, 0.4);
    }
    else {
        leg = new TLegend(0.2, 0.55, 0.4, 0.7);
    }
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    if ( hData ) {
        leg->AddEntry( hData, "Data", "p" );
    }
    leg->AddEntry( hReco, "Reco", "p" );
    leg->AddEntry( hGen, "Gen", "p" );
    if ( hRef ) {
        leg->AddEntry( hRef, "Ref", "p" );        
    }
    if ( hRefSel ) {
        leg->AddEntry( hRefSel, "RefSel", "p" );
    }
    leg->Draw();
}

//________________
// Draw reco/gen, ref/gen, refSel/gen and data/gen dijet eta ratios
void drawDijetToGenRatio(TCanvas *c, TH1D *hReco2Gen, TH1D *hRef2Gen = nullptr,
                         TH1D *hRefSel2Gen = nullptr, TH1D *hData = nullptr,
                         int ptLow = 50, int ptHi = 60, 
                         bool isCM = false, bool isFB = false,
                         int collisionSystem = 1, double energy = 8.16) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TLegend *leg;

    double xRange[2] = {-3.2, 3.2};
    double yRange[2] = {0.7, 1.3};

    if (isCM) {
        if (!isFB) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] =0.7; yRange[1] = 1.3;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.5;
            yRange[0] = 0.7; yRange[1] = 1.3;
        }
    }

    c->cd();
    setPadStyle();
    hReco2Gen->Draw();
    if (hRef2Gen) {
        hRef2Gen->Draw("same");
    }
    if (hRefSel2Gen) {
        hRefSel2Gen->Draw("same");
    }
    if (hData) {
        hData->Draw("same");
    }
    hReco2Gen->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hReco2Gen->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi));
    leg = new TLegend(0.2, 0.25, 0.4, 0.4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    if (hData) {
        leg->AddEntry(hData, "Data / Gen", "p");
    }
    leg->AddEntry(hReco2Gen, "Reco / Gen", "p");
    if (hRef2Gen) {
        leg->AddEntry(hRef2Gen, "Ref / Gen", "p");
    }
    if (hRefSel2Gen) {
        leg->AddEntry(hRefSel2Gen, "RefSel / Gen", "p");
    }
    leg->Draw();
}

//________________
// Plot dijet eta comparison of reco, ref and refSel to gen
void plotDijetClosures(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    
    // collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Dijet ptAve binning
    int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
                             100, 110,  120, 130, 140,
                             150, 160,  180, 200, 250, 
                             300, 500};
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);

    std::vector<int> dijetPtVals; 
    dijetPtVals.assign(dijetPtNewVals, dijetPtNewVals + sizeOfPtVals);

    std::cout << "Number of pT bins: " << dijetPtVals.size() << std::endl;

    // Bins for projections from 3D
    std::vector<int> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55 };
    std::vector<int> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94 };

    TH1D *h1DProj[ ptDijetBinLow.size() ];
    TH1D *h1DDirect[ ptDijetBinLow.size() ];

    //
    // Lab frame
    //
    TH1D *hRecoDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hRefSelDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hReco2GenDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hRef2GenDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hRefSel2GenDijetEta1DLab[ ptDijetBinLow.size() ];

    //
    // CM frame
    //
    TH1D *hRecoDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hRefSelDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hReco2GenDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hRef2GenDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hRefSel2GenDijetEta1DCM[ ptDijetBinLow.size() ];

    //
    // Forward/backward ratios
    //

    TH1D *hRecoDijetEtaCMForward1D[ ptDijetBinLow.size() ];
    TH1D *hRecoDijetEtaCMBackward1D[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEtaCMForward1D[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEtaCMBackward1D[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEtaCMForward1D[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEtaCMBackward1D[ ptDijetBinLow.size() ];

    TH1D *hRecoDijetFBEtaCM1D[ ptDijetBinLow.size() ];
    TH1D *hGenDijetFBEtaCM1D[ ptDijetBinLow.size() ];
    TH1D *hRefDijetFBEtaCM1D[ ptDijetBinLow.size() ];

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );
    TLegend *leg;
    TCanvas *c = new TCanvas( "c", "c", 1000, 1000 );

    // Loop over dijet ptAve bins
    for (unsigned int i = 0; i < ptDijetBinLow.size(); i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{1000}, canvY{1000};

        //
        // Lab frame
        //
        hRecoDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1DWeighted_%d", i) ) );
        hRecoDijetEta1DLab[i]->SetName( Form("hRecoDijetEta1DLab_%d", i) );
        set1DStyle( hRecoDijetEta1DLab[i], 0 );
        rescaleEta( hRecoDijetEta1DLab[i] );
        hGenDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1DWeighted_%d", i) ) );
        hGenDijetEta1DLab[i]->SetName( Form("hGenDijetEta1DLab_%d", i) );
        set1DStyle( hGenDijetEta1DLab[i], 4 );
        rescaleEta( hGenDijetEta1DLab[i] );
        hRefDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1DWeighted_%d", i) ) );
        hRefDijetEta1DLab[i]->SetName( Form("hRefDijetEta1DLab_%d", i) );
        set1DStyle( hRefDijetEta1DLab[i], 1 );
        rescaleEta( hRefDijetEta1DLab[i] );
        hRefSelDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefSelDijetEta1DWeighted_%d", i) ) );
        hRefSelDijetEta1DLab[i]->SetName( Form("hRefSelDijetEta1DLab_%d", i) );
        set1DStyle( hRefSelDijetEta1DLab[i], 2 );
        rescaleEta( hRefSelDijetEta1DLab[i] );

        hReco2GenDijetEta1DLab[i] = dynamic_cast<TH1D*>( hRecoDijetEta1DLab[i]->Clone( Form("hReco2GenDijetEta1DLab_%d", i) ) );
        hReco2GenDijetEta1DLab[i]->Divide( hReco2GenDijetEta1DLab[i], hGenDijetEta1DLab[i], 1., 1. );
        hRef2GenDijetEta1DLab[i] = dynamic_cast<TH1D*>( hRefDijetEta1DLab[i]->Clone( Form("hRef2GenDijetEta1DLab_%d", i) ) );
        hRef2GenDijetEta1DLab[i]->Divide( hRef2GenDijetEta1DLab[i], hGenDijetEta1DLab[i], 1., 1., "b" );
        hRefSel2GenDijetEta1DLab[i] = dynamic_cast<TH1D*>( hRefSelDijetEta1DLab[i]->Clone( Form("hRefSel2GenDijetEta1DLab_%d", i) ) );
        hRefSel2GenDijetEta1DLab[i]->Divide( hRefSel2GenDijetEta1DLab[i], hGenDijetEta1DLab[i], 1., 1., "b" );

        //
        // CM frame
        //
        hRecoDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1DCMWeighted_%d", i) ) );
        hRecoDijetEta1DCM[i]->SetName( Form("hRecoDijetEta1DCM_%d", i) );
        set1DStyle( hRecoDijetEta1DCM[i], 0 );
        rescaleEta( hRecoDijetEta1DCM[i] );
        hGenDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1DCMWeighted_%d", i) ) );
        hGenDijetEta1DCM[i]->SetName( Form("hGenDijetEta1DCM_%d", i) );
        set1DStyle( hGenDijetEta1DCM[i], 4 );
        rescaleEta( hGenDijetEta1DCM[i] );
        hRefDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1DCMWeighted_%d", i) ) );
        hRefDijetEta1DCM[i]->SetName( Form("hRefDijetEta1DCM_%d", i) );
        set1DStyle( hRefDijetEta1DCM[i], 1 );
        rescaleEta( hRefDijetEta1DCM[i] );
        hRefSelDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefSelDijetEta1DCMWeighted_%d", i) ) );
        hRefSelDijetEta1DCM[i]->SetName( Form("hRefSelDijetEta1DCM_%d", i) );
        set1DStyle( hRefSelDijetEta1DCM[i], 2 );
        rescaleEta( hRefSelDijetEta1DCM[i] );
        
        hReco2GenDijetEta1DCM[i] = dynamic_cast<TH1D*>( hRecoDijetEta1DCM[i]->Clone( Form("hReco2GenDijetEta1DCM_%d", i) ) );
        hReco2GenDijetEta1DCM[i]->Divide( hReco2GenDijetEta1DCM[i], hGenDijetEta1DCM[i], 1., 1. /* , "b" */ );
        hRef2GenDijetEta1DCM[i] = dynamic_cast<TH1D*>( hRefDijetEta1DCM[i]->Clone( Form("hRef2GenDijetEta1DCM_%d", i) ) );
        hRef2GenDijetEta1DCM[i]->Divide( hRef2GenDijetEta1DCM[i], hGenDijetEta1DCM[i], 1., 1., "b" );
        hRefSel2GenDijetEta1DCM[i] = dynamic_cast<TH1D*>( hRefSelDijetEta1DCM[i]->Clone( Form("hRefSel2GenDijetEta1DCM_%d", i) ) );
        hRefSel2GenDijetEta1DCM[i]->Divide( hRefSel2GenDijetEta1DCM[i], hGenDijetEta1DCM[i], 1., 1., "b" );

        //
        // Forward/backward ratios
        //
        hRecoDijetEtaCMForward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMForward1DWeighted_%d", i) ) );
        hRecoDijetEtaCMForward1D[i]->SetName( Form("hRecoDijetEtaCMForward1D_%d", i) );
        hRecoDijetEtaCMBackward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMBackward1DWeighted_%d", i) ) );
        hRecoDijetEtaCMBackward1D[i]->SetName( Form("hRecoDijetEtaCMBackward1D_%d", i) );
        hGenDijetEtaCMForward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMForward1DWeighted_%d", i) ) );
        hGenDijetEtaCMForward1D[i]->SetName( Form("hGenDijetEtaCMForward1D_%d", i) );
        hGenDijetEtaCMBackward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMBackward1DWeighted_%d", i) ) );
        hGenDijetEtaCMBackward1D[i]->SetName( Form("hGenDijetEtaCMBackward1D_%d", i) );
        hRefDijetEtaCMForward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMForward1DWeighted_%d", i) ) );
        hRefDijetEtaCMForward1D[i]->SetName( Form("hRefDijetEtaCMForward1D_%d", i) );
        hRefDijetEtaCMBackward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMBackward1DWeighted_%d", i) ) );
        hRefDijetEtaCMBackward1D[i]->SetName( Form("hRefDijetEtaCMBackward1D_%d", i) );

        rescaleForwardBackward( hRecoDijetEtaCMForward1D[i], hRecoDijetEtaCMBackward1D[i] );
        rescaleForwardBackward( hGenDijetEtaCMForward1D[i], hGenDijetEtaCMBackward1D[i] );
        rescaleForwardBackward( hRefDijetEtaCMForward1D[i], hRefDijetEtaCMBackward1D[i] );
        set1DStyle( hRecoDijetEtaCMForward1D[i], 0 );
        set1DStyle( hGenDijetEtaCMForward1D[i], 2 );
        set1DStyle( hRefDijetEtaCMForward1D[i], 1 );

        hRecoDijetFBEtaCM1D[i] = dynamic_cast<TH1D*>( hRecoDijetEtaCMForward1D[i]->Clone( Form("hRecoDijetFBEtaCM1D_%d", i) ) );
        hRecoDijetFBEtaCM1D[i]->Divide( hRecoDijetEtaCMBackward1D[i] );
        hRecoDijetFBEtaCM1D[i]->GetYaxis()->SetTitle("Forward/Backward");
        hGenDijetFBEtaCM1D[i] = dynamic_cast<TH1D*>( hGenDijetEtaCMForward1D[i]->Clone( Form("hGenDijetFBEtaCM1D_%d", i) ) );
        hGenDijetFBEtaCM1D[i]->Divide( hGenDijetEtaCMBackward1D[i] );
        hGenDijetFBEtaCM1D[i]->GetYaxis()->SetTitle("Forward/Backward");
        hRefDijetFBEtaCM1D[i] = dynamic_cast<TH1D*>( hRefDijetEtaCMForward1D[i]->Clone( Form("hRefDijetFBEtaCM1D_%d", i) ) );
        hRefDijetFBEtaCM1D[i]->Divide( hRefDijetEtaCMBackward1D[i] );
        hRefDijetFBEtaCM1D[i]->GetYaxis()->SetTitle("Forward/Backward");

        //
        // Plot comparisons in the lab frame
        //
        drawDijetToGenComparison(c, hRecoDijetEta1DLab[i], hRefDijetEta1DLab[i], hGenDijetEta1DLab[i], hRefSelDijetEta1DLab[i], nullptr,
                            dijetPtVals[i], dijetPtVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaLab_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );        

        //
        // Plot ratios in the lab frame
        //
        drawDijetToGenRatio(c, hReco2GenDijetEta1DLab[i], hRef2GenDijetEta1DLab[i], hRefSel2GenDijetEta1DLab[i], nullptr,
                       dijetPtVals[i], dijetPtVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Plot comparisons in the CM frame
        //

        drawDijetToGenComparison(c, hRecoDijetEta1DCM[i], hRefDijetEta1DCM[i], hGenDijetEta1DCM[i], hRefSelDijetEta1DCM[i], nullptr,
                            dijetPtVals[i], dijetPtVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaCM_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Plot ratios in the CM frame
        //
        drawDijetToGenRatio(c, hReco2GenDijetEta1DCM[i], hRef2GenDijetEta1DCM[i], hRefSel2GenDijetEta1DCM[i], nullptr,
                       dijetPtVals[i], dijetPtVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaCM_RecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Forward/backward ratios
        //
        drawDijetToGenComparison(c, hRecoDijetFBEtaCM1D[i], hRefDijetFBEtaCM1D[i], hGenDijetFBEtaCM1D[i], nullptr, nullptr,
                            dijetPtVals[i], dijetPtVals[i+1], true, true, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaFB_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );
            
    } // end loop over dijet pt bins
}

//________________
// Plot eta and pT distributions of reco and gen jets to look at the closures
// after the JECs are applied. Corrections can be studies in the bins of ptHat.
//
void plotInclusiveJetJECClosures(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.6;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Create vector of ptHat and jet pT bins for projections
    int ptHatStart = 15;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins{1}; // 30

    // Jet pT binning
    int jetPtStart = 5;
    int jetPtStep = 10; // Starting from 5 GeV: jetPtStart + (jetPtBinsLow(i) - 1) * jetPtStep
    int jetPtBinsMax = 150;
    std::vector<int> jetPtBinsLow{3,  5,  7,   9,  12,   3,   5,   7}; // 45, 35, 45, 55, 65, 75, 85, 125
    std::vector<int> jetPtBinsHigh {5,  7,  9,  12, 150, 150, 150, 150}; // 65, 35, 45, 55, 65, 75, 85, 125

    // Eta binning
    // 52 bins from (-5.2, 5.2)
    int nEtaBins = 52;
    double etaStep = 0.2;
    std::vector<int> jetEtaBinsLow{19, 12, 38};
    std::vector<int> jetEtaBinsHigh{35, 16, 42};

    // Retrieve histograms
    TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRecoInclusiveJetPtEtaPtHat"));
    // TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRecoMatchedJetPtEtaPtHat"));
    TH3D *hGenPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hGenInclusiveJetPtEtaPtHat"));
    TH3D *hRefPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRefInclusiveJetPtEtaPtHat"));
    TH3D *hRefSelPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRefSelInclusiveJetPtEtaPtHat"));

    // Check that each histogram exist
    if ( !hRecoPtEtaPtHat ) {
        std::cerr << "Histogram hRecoPtEtaPtHat not found in file." << std::endl; return;
    }
    if ( !hGenPtEtaPtHat ) {
        std::cerr << "Histogram hGenPtEtaPtHat not found in file." << std::endl; return;
    }
    if ( !hRefPtEtaPtHat ) {
        std::cerr << "Histogram hRefPtEtaPtHat not found in file." << std::endl; return;
    }
    if ( !hRefSelPtEtaPtHat ) {
        std::cerr << "Histogram hRefSelPtEtaPtHat not found in file." << std::endl; return;
    }


    ptHatBinsMax = hRecoPtEtaPtHat->GetZaxis()->GetNbins();

    // //
    // // 2D distributions
    // //
    // TH2D *hRecoPtVsPtHat = dynamic_cast<TH2D *>(hRecoPtEtaPtHat->Project3D("yz"));
    // hRecoPtVsPtHat->SetName("hRecoPtVsPtHat");
    // set2DStyle(hRecoPtVsPtHat);
    // TH2D *hGenPtVsPtHat = dynamic_cast<TH2D *>(hGenPtEtaPtHat->Project3D("yz"));
    // hGenPtVsPtHat->SetName("hGenPtVsPtHat");
    // set2DStyle(hGenPtVsPtHat);

    // double ptRange[2] = {0, 250};
    // TCanvas *cJetPtVsPtHat = new TCanvas("cJetPtVsPtHat", "cJetPtVsPtHat", 1000, 500);
    // cJetPtVsPtHat->Divide(2, 1);
    // cJetPtVsPtHat->cd(1);
    // setPadStyle();
    // hRecoPtVsPtHat->Draw("colz");
    // hRecoPtVsPtHat->GetXaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    // hRecoPtVsPtHat->GetYaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    // // gPad->SetLogz();
    // hRecoPtVsPtHat->GetXaxis()->SetTitle("#hat{p}_{T} (GeV)");
    // hRecoPtVsPtHat->GetYaxis()->SetTitle("Reco p_{T}^{jet} (GeV)");
    // plotCMSHeader(collSystem, energy);
    // cJetPtVsPtHat->cd(2);
    // setPadStyle();
    // hGenPtVsPtHat->Draw("colz");
    // hGenPtVsPtHat->GetXaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    // hGenPtVsPtHat->GetYaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    // // gPad->SetLogz();
    // hGenPtVsPtHat->GetXaxis()->SetTitle("#hat{p}_{T} (GeV)");
    // hGenPtVsPtHat->GetYaxis()->SetTitle("Gen p_{T}^{jet} (GeV)");
    // plotCMSHeader(collSystem, energy);

    //
    // Declare canvases and histograms
    //
    TCanvas *cPtVsEta[ptHatBins.size()];
    TCanvas *cClosureEta[ptHatBins.size()][jetPtBinsLow.size()];
    TCanvas *cClosurePt[ptHatBins.size()][jetEtaBinsLow.size()];

    TH2D *hRecoPtVsEta[ptHatBins.size()];
    TH2D *hGenPtVsEta[ptHatBins.size()];
    TH1D *hRecoPt[ptHatBins.size()][jetEtaBinsLow.size()];
    TH1D *hGenPt[ptHatBins.size()][jetEtaBinsLow.size()];

    TH1D *hRecoEta[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hGenEta[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hRefEta[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hRefSelEta[ptHatBins.size()][jetPtBinsLow.size()];

    TH1D *hReco2GenEta[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hRef2GenEta[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hRefSel2GenEta[ptHatBins.size()][jetPtBinsLow.size()];

    TCanvas *c = new TCanvas("c", "c", 1200, 1000);


    //
    // Perform analysis and build comparisons for different ptHat intervals.
    // First part loops over jet pT, second part over jet eta.
    //

    // Loop over ptHat bins
    for (unsigned int i = 0; i < ptHatBins.size(); i++) {

        double ptHatLow = hRecoPtEtaPtHat->GetZaxis()->GetBinLowEdge(ptHatBins[i]);

        //
        // Set ptHat range
        //
        hRecoPtEtaPtHat->GetZaxis()->SetRange(ptHatBins[i],  ptHatBinsMax);
        hGenPtEtaPtHat->GetZaxis()->SetRange(ptHatBins[i],  ptHatBinsMax);
        hRefPtEtaPtHat->GetZaxis()->SetRange(ptHatBins[i],  ptHatBinsMax);
        hRefSelPtEtaPtHat->GetZaxis()->SetRange(ptHatBins[i],  ptHatBinsMax);
        

        // // Make 2D histograms
        // hRecoPtVsEta[i] = dynamic_cast<TH2D *>(hRecoPtEtaPtHat->Project3D("yx"));
        // hRecoPtVsEta[i]->SetName(Form("hRecoPtVsEta_%d", i));
        // set2DStyle(hRecoPtVsEta[i]);

        // hGenPtVsEta[i] = dynamic_cast<TH2D *>(hGenPtEtaPtHat->Project3D("yx"));
        // hGenPtVsEta[i]->SetName(Form("hGenPtVsEta_%d", i));
        // set2DStyle(hGenPtVsEta[i]);

        // cPtVsEta[i] = new TCanvas(Form("cPtVsEta_%d", i), Form("cPtVsEta_%d", i), 1000, 500);
        // cPtVsEta[i]->Divide(2, 1);
        // cPtVsEta[i]->cd(1);
        // setPadStyle();
        // hRecoPtVsEta[i]->Draw("colz");
        // hRecoPtVsEta[i]->GetXaxis()->SetRangeUser(-5.2, 5.2);
        // hRecoPtVsEta[i]->GetYaxis()->SetRangeUser(0, 120);
        // hRecoPtVsEta[i]->GetXaxis()->SetTitle("Reco #eta^{jet}");
        // hRecoPtVsEta[i]->GetYaxis()->SetTitle("Reco p_{T}^{jet} (GeV)");
        // plotCMSHeader(collisionSystem, collisionEnergy);
        // gPad->SetLogz();
        // t.DrawLatexNDC(xTextPosition, yTextPosition, Form("#hat{p}_{T} > %4.0f GeV", ptHatLow));
        // cPtVsEta[i]->cd(2);
        // setPadStyle();
        // hGenPtVsEta[i]->Draw("colz");
        // hGenPtVsEta[i]->GetXaxis()->SetRangeUser(-5.2, 5.2);
        // hGenPtVsEta[i]->GetYaxis()->SetRangeUser(0, 120);
        // hGenPtVsEta[i]->GetXaxis()->SetTitle("Gen #eta^{jet}");
        // hGenPtVsEta[i]->GetYaxis()->SetTitle("Gen p_{T}^{jet} (GeV)");
        // plotCMSHeader(collisionSystem, collisionEnergy);
        // gPad->SetLogz();
        // t.DrawLatexNDC(xTextPosition, yTextPosition, Form("#hat{p}_{T} > %4.0f GeV", ptHatLow));

        //
        // Loop over jet pT bins
        //
        for (unsigned int j = 0; j < jetPtBinsLow.size(); j++) {

            //
            // Make projections of the 3D histograms
            //

            double ptLow = hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[j]);
            double ptHi = hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[j]);

            // Reco jets
            hRecoEta[i][j] = dynamic_cast<TH1D *>(hRecoPtEtaPtHat->ProjectionX(Form("hRecoEta_%d_%d", i, j),
                                                                               jetPtBinsLow[j], jetPtBinsHigh[j],
                                                                               ptHatBins[i], ptHatBinsMax));
            set1DStyle(hRecoEta[i][j], 0, kTRUE);

            // Gen jets
            hGenEta[i][j] = dynamic_cast<TH1D *>(hGenPtEtaPtHat->ProjectionX(Form("hGenEta_%d_%d", i, j),
                                                 jetPtBinsLow[j], jetPtBinsHigh[j],
                                                 ptHatBins[i], ptHatBinsMax));
            set1DStyle(hGenEta[i][j], 5, kTRUE);

            // Ref jets
            hRefEta[i][j] = dynamic_cast<TH1D *>(hRefPtEtaPtHat->ProjectionX(Form("hRefEta_%d_%d", i, j),
                                                 jetPtBinsLow[j], jetPtBinsHigh[j],
                                                 ptHatBins[i], ptHatBinsMax));
            set1DStyle(hRefEta[i][j], 1, kTRUE);
        
            // RefSel jets
            hRefSelEta[i][j] = dynamic_cast<TH1D *>(hRefSelPtEtaPtHat->ProjectionX(Form("hRefSelEta_%d_%d", i, j),
                                                 jetPtBinsLow[j], jetPtBinsHigh[j],
                                                 ptHatBins[i], ptHatBinsMax));
            set1DStyle(hRefSelEta[i][j], 2, kTRUE);

            //
            // Ratios of reco, ref and refSel to gen
            //
            hReco2GenEta[i][j] = dynamic_cast<TH1D *>(hRecoEta[i][j]->Clone(Form("hReco2GenEta_%d_%d", i, j)));
            hReco2GenEta[i][j]->Divide(hReco2GenEta[i][j], hGenEta[i][j], 1., 1., "b");
            hRef2GenEta[i][j] = dynamic_cast<TH1D *>(hRefEta[i][j]->Clone(Form("hRef2GenEta_%d_%d", i, j)));
            hRef2GenEta[i][j]->Divide(hRef2GenEta[i][j], hGenEta[i][j], 1., 1., "b");
            hRefSel2GenEta[i][j] = dynamic_cast<TH1D *>(hRefSelEta[i][j]->Clone(Form("hRefSel2GenEta_%d_%d", i, j)));
            hRefSel2GenEta[i][j]->Divide(hRefSel2GenEta[i][j], hGenEta[i][j], 1., 1., "b");

            //
            // Plot comparisons
            //
            drawSingleJetToGenComparison(c, hRecoEta[i][j], hRefEta[i][j], hGenEta[i][j], hRefSelEta[i][j], nullptr,
                                         hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[j]), 
                                         hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[j]),
                                         ptHatLow, false, false, collisionSystem, collisionEnergy);
            c->SaveAs(Form("%s/%s_jetEta_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), 
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[j]), 
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[j])));

            //
            // Plot ratios
            //
            drawSingleJetToGenRatio(c, hReco2GenEta[i][j], hRef2GenEta[i][j], hRefSel2GenEta[i][j], nullptr,
                                    hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[j]), 
                                    hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[j]),
                                    ptHatLow, false, false, collisionSystem, collisionEnergy);
            c->SaveAs(Form("%s/%s_jetEta_RecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), 
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[j]), 
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[j])));


            // // Create canvas
            // cClosureEta[i][j] = new TCanvas(Form("cClosureEta_%d_%d", i, j), Form("cClosureEta_%d_%d", i, j), 700, 800);

            // // Plot distributions
            // drawJecComparison(cClosureEta[i][j], hRecoEta[i][j], hGenEta[i][j],
            //                   hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge( jetPtBinsLow[j] ), 
            //                   hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge( jetPtBinsHigh[j] ),
            //                   ptHatLow,
            //                   "Reco", "Gen", collisionSystem, collisionEnergy, false, true, false,
            //                   hRecoPtEtaPtHat->GetXaxis()->GetBinLowEdge( 1 ),
            //                   hRecoPtEtaPtHat->GetXaxis()->GetBinUpEdge( hRecoPtEtaPtHat->GetXaxis()->GetNbins() )
            //                 );

        } // for (unsigned int j = 0; j < jetPtBinsLow.size(); j++)


        // //
        // // Loop over jet eta bins
        // //
        // for (unsigned int j = 0; j < jetEtaBinsLow.size(); j++) {

        //     // Create canvas
        //     cClosurePt[i][j] = new TCanvas(Form("cClosurePt_%d_%d", i, j), Form("cClosurePt_%d_%d", i, j), 700, 800);
        //     // Make projection of the 3D histograms
        //     hRecoPt[i][j] = dynamic_cast<TH1D *>(hRecoPtEtaPtHat->ProjectionY(Form("hRecoPt_%d_%d", i, j),
        //                                                                       jetEtaBinsLow[j], jetEtaBinsHigh[j],
        //                                                                       ptHatBins[i], ptHatBinsMax));
        //     set1DStyle(hRecoPt[i][j], 0, kFALSE);
        //     hGenPt[i][j] = dynamic_cast<TH1D *>(hGenPtEtaPtHat->ProjectionY(Form("hGenPt_%d_%d", i, j),
        //                                                                     jetEtaBinsLow[j], jetEtaBinsHigh[j],
        //                                                                     ptHatBins[i], ptHatBinsMax));
        //     set1DStyle(hGenPt[i][j], 1, kFALSE);

        //     // Plot distributions
        //     drawJecComparison(cClosurePt[i][j], hRecoPt[i][j], hGenPt[i][j],
        //                       hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge( 1 ),
        //                       hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge( hRecoPtEtaPtHat->GetYaxis()->GetNbins() ),
        //                       ptHatLow,
        //                       "Reco", "Gen", collisionSystem, collisionEnergy, false, true, true,
        //                       hRecoPtEtaPtHat->GetXaxis()->GetBinLowEdge( jetEtaBinsLow[j] ), 
        //                       hRecoPtEtaPtHat->GetXaxis()->GetBinUpEdge( jetEtaBinsHigh[j] ));

        // } // for (unsigned int j = 0; j < jetEtaBinsLow.size(); j++)
    } // for (unsigned int i = 0; i < ptHatBins.size(); i++)
}

//________________
//
// Plot comparison of inclusive jet distributions for different pT and eta bins
// Compare both reco data to reco MC and reco data to gen MC
//
void plotData2McInclusiveJetComparison(TFile *fData, TFile *fMc, int collSystem = 1, double energy = 8.16, TString date = "20250224") {
    
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString prefix = "mb_";
    TString dataFName = fData->GetName();
    if ( dataFName.Contains("jet60", TString::kIgnoreCase) ) {
        prefix = "jet60";
    }
    else if ( dataFName.Contains("jet80", TString::kIgnoreCase) ) {
        prefix = "jet80";
    }
    else if ( dataFName.Contains("jet100", TString::kIgnoreCase) ) {
        prefix = "jet100";
    }

    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", (int)(1000 * energy));
    collSystemStr += Form("_%s", prefix.Data());

    TH3D *hRecoDataPtEtaPtHat = dynamic_cast<TH3D *>(fData->Get("hRecoInclusiveJetPtEtaPtHat"));
    hRecoDataPtEtaPtHat->SetName("hRecoDataPtEtaPtHat");
    TH3D *hRecoMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hRecoInclusiveJetPtEtaPtHat"));
    hRecoMcPtEtaPtHat->SetName("hRecoMcPtEtaPtHat");
    // TH3D *hRecoMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hRecoMatchedJetPtEtaPtHat"));
    // hRecoMcPtEtaPtHat->SetName("hRecoMcPtEtaPtHat");
    TH3D *hGenMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hGenInclusiveJetPtEtaPtHat"));
    hGenMcPtEtaPtHat->SetName("hGenMcPtEtaPtHat");

    // Create vector of ptHat and jet pT bins for projections
    int ptHatStart = 10;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins{4}; // 30

    // Jet pT binning
    int jetPtStart = 5;
    int jetPtStep = 10; // Starting from 5 GeV: jetPtStart + (jetPtBinsLow(i) - 1) * jetPtStep
    int jetPtBinsMax = 150;
    std::vector<int> jetPtBinsLow {4,  5,  6,  7,  8,  9,  12,   4,   5,   7,   9}; 
    std::vector<int> jetPtBinsHigh{5,  6,  7,  8,  9, 12,  50, 150, 150, 150, 150};

    // Eta binning
    // 52 bins from (-5.2, 5.2)
    int nEtaBins = 52;
    double etaStep = 0.2;
    std::vector<int> jetEtaBinsLow{19, 12, 38};
    std::vector<int> jetEtaBinsHigh{35, 16, 42};

    TH1D *hRecoDataEta[jetPtBinsLow.size()];
    TH1D *hRecoMcEta[jetPtBinsLow.size()];
    TH1D *hGenMcEta[jetPtBinsLow.size()];

    TCanvas *cClosureEtaReco2Reco[jetPtBinsLow.size()];
    TCanvas *cClosureEtaReco2Gen[jetPtBinsLow.size()];

    // Jet pT binning
    for (unsigned int i = 0; i < jetPtBinsLow.size(); i++) {

        //
        // Reco data 2 Reco MC comparison
        //

        // Create canvas
        cClosureEtaReco2Reco[i] = new TCanvas(Form("cClosureEtaReco2Reco_%d", i), Form("cClosureEtaReco2Reco_%d", i), 700, 800);

        // Make projection of the 3D histograms
        hRecoDataEta[i] = dynamic_cast<TH1D *>(hRecoDataPtEtaPtHat->ProjectionX(Form("hRecoDataEta_%d", i),
                                                                                 jetPtBinsLow[i], jetPtBinsHigh[i],
                                                                                 1, -1));
        set1DStyle(hRecoDataEta[i], 0, kTRUE);

        hRecoMcEta[i] = dynamic_cast<TH1D *>(hRecoMcPtEtaPtHat->ProjectionX(Form("hRecoMcEta_%d", i),
                                                                              jetPtBinsLow[i], jetPtBinsHigh[i],
                                                                              ptHatBins[0], hRecoMcPtEtaPtHat->GetNbinsZ()));
        set1DStyle(hRecoMcEta[i], 1, kTRUE);

        // Plot distributions
        drawJecComparison(cClosureEtaReco2Reco[i], hRecoDataEta[i], hRecoMcEta[i],
            hRecoDataPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[i]),
            hRecoDataPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[i]),
            ptHatStart + (ptHatBins[0] - 1) * ptHatStep,
            "Data", "Reco MC", collSystem, energy, false, true,
            hRecoDataPtEtaPtHat->GetXaxis()->GetBinLowEdge(1),
            hRecoDataPtEtaPtHat->GetXaxis()->GetBinUpEdge(hRecoDataPtEtaPtHat->GetXaxis()->GetNbins()));
        cClosureEtaReco2Reco[i]->SaveAs( Form("%s/%s_closureEta_Reco2Reco_pt_%d__%d.pdf", 
                                              date.Data(), collSystemStr.Data(), 
                                              int(hRecoDataPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[i])),
                                              int(hRecoDataPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsHigh[i])) ) );


        //
        // Reco data 2 Gen MC comparison
        //

        // Create canvas
        cClosureEtaReco2Gen[i] = new TCanvas(Form("cClosureEtaReco2Gen_%d", i), Form("cClosureEtaReco2Gen_%d", i), 700, 800);
        // Make projection of the 3D histograms
        hGenMcEta[i] = dynamic_cast<TH1D *>(hGenMcPtEtaPtHat->ProjectionX(Form("hGenMcEta_%d", i),
                                                                            jetPtBinsLow[i], jetPtBinsHigh[i],
                                                                            ptHatBins[0], hGenMcPtEtaPtHat->GetNbinsZ()));
        set1DStyle(hGenMcEta[i], 1, kTRUE);

        // Plot distributions
        drawJecComparison(cClosureEtaReco2Gen[i], hRecoDataEta[i], hGenMcEta[i],
                          hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[i]),
                          hGenMcPtEtaPtHat->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[i]),
                          ptHatStart + (ptHatBins[0] - 1) * ptHatStep,
                          "Data", "Gen MC", collSystem, energy, false, true,
                          hGenMcPtEtaPtHat->GetXaxis()->GetBinLowEdge(1),
                          hGenMcPtEtaPtHat->GetXaxis()->GetBinUpEdge(hGenMcPtEtaPtHat->GetXaxis()->GetNbins())
                        );

        // Save canvas
        cClosureEtaReco2Gen[i]->SaveAs( Form("%s/%s_closureEta_Reco2Gen_pt_%d_%d.pdf", 
                                             date.Data(), collSystemStr.Data(), 
                                             int(hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsLow[i])) ,
                                             int(hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge(jetPtBinsHigh[i])) )
                                            );

    } // for (unsigned int i = 0; i < jetPtBinsLow.size(); i++)
}

//________________
//
// Plot comparison of dijet eta distributions for different pTave bins
// Compare both reco data to reco MC and reco data to gen MC
//
void plotData2McDijetComparison(TFile *fData, TFile *fMc, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250224") {

    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString prefix = "mb";
    TString dataFName = fData->GetName();
    if ( dataFName.Contains("jet60", TString::kIgnoreCase) ) {
        prefix = "jet60";
    }
    else if ( dataFName.Contains("jet80", TString::kIgnoreCase) ) {
        prefix = "jet80";
    }
    else if ( dataFName.Contains("jet100", TString::kIgnoreCase) ) {
        prefix = "jet100";
    }

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", (int)(1000 * collisionEnergy));
    collSystemStr += Form("_%s", prefix.Data());

    const int dijetPtBins{16};
    double dijetPtVals[dijetPtBins+1] {  50.,  60.,   70.,  80.,  90.,
                                         100., 110.,  120., 130., 140.,
                                         150., 160.,  180., 200., 250., 
                                         300., 500.};

    //
    // Declare histograms
    //

    // Data
    TH1D *hRecoDataDijetEta[dijetPtBins];
    TH1D *hRecoDataDijetEtaCM[dijetPtBins];
    TH1D *hRecoDataDijetEtaForward[dijetPtBins];
    TH1D *hRecoDataDijetEtaBackward[dijetPtBins];
    TH1D *hRecoDataDijetFBRatio[dijetPtBins];

    // MC
    TH1D *hRecoMcDijetEta[dijetPtBins];
    TH1D *hRecoMcDijetEtaCM[dijetPtBins];
    TH1D *hRecoMcDijetEtaForward[dijetPtBins];
    TH1D *hRecoMcDijetEtaBackward[dijetPtBins];
    TH1D *hRecoMcDijetFBRatio[dijetPtBins];

    TH1D *hGenMcDijetEta[dijetPtBins];
    TH1D *hGenMcDijetEtaCM[dijetPtBins];
    TH1D *hGenMcDijetEtaForward[dijetPtBins];
    TH1D *hGenMcDijetEtaBackward[dijetPtBins];
    TH1D *hGenMcDijetFBRatio[dijetPtBins];

    TH1D *hRefMcDijetEta[dijetPtBins];
    TH1D *hRefMcDijetEtaCM[dijetPtBins];
    TH1D *hRefMcDijetEtaForward[dijetPtBins];
    TH1D *hRefMcDijetEtaBackward[dijetPtBins];
    TH1D *hRefMcDijetFBRatio[dijetPtBins];

    // Ratio of full distributions to gen
    TH1D *hRecoData2GenDijetEta[dijetPtBins];
    TH1D *hRecoMc2GenDijetEta[dijetPtBins];
    TH1D *hRef2GenDijetEta[dijetPtBins];

    TH1D *hRecoData2GenDijetEtaCM[dijetPtBins];
    TH1D *hRecoMc2GenDijetEtaCM[dijetPtBins];
    TH1D *hRef2GenDijetEtaCM[dijetPtBins];

    //
    // Loop over the dijet pT bins and read histograms from the files
    //
    for (int i = 0; i < dijetPtBins; ++i) {

        double integralForward{0};
        double integralBackward{0};

        //
        // Data
        //

        // Eta full (lab frame)
        hRecoDataDijetEta[i] = dynamic_cast<TH1D*>( fData->Get( Form("hRecoDijetEta1DWeighted_%d", i) ) );
        if ( !hRecoDataDijetEta[i] ) {
            std::cerr << Form("Data histogram not found: hRecoDijetEta1DWeighted_%d", i) << std::endl; return;
        }
        hRecoDataDijetEta[i]->SetName( Form("hRecoDataDijetEta_%d", i) );
        set1DStyle( hRecoDataDijetEta[i], 2 );
        rescaleEta( hRecoDataDijetEta[i] );

        // Eta full (CM frame)
        hRecoDataDijetEtaCM[i] = dynamic_cast<TH1D*>( fData->Get( Form("hRecoDijetEta1DCMWeighted_%d", i) ) );
        if ( !hRecoDataDijetEtaCM[i] ) {
            std::cerr << Form("Data histogram not found: hRecoDijetEta1DCMWeighted_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaCM[i]->SetName( Form("hRecoDataDijetEtaCM_%d", i) );
        set1DStyle( hRecoDataDijetEtaCM[i], 2 );
        rescaleEta( hRecoDataDijetEtaCM[i] );

        // Eta forward (CM frame)
        hRecoDataDijetEtaForward[i] = dynamic_cast<TH1D*>( fData->Get( Form("hRecoDijetEtaCMForward1DWeighted_%d", i) ) );
        if ( !hRecoDataDijetEtaForward[i] ) {
            std::cerr << Form("Data histogram not found: hRecoDijetEtaCMForward1DWeighted_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaForward[i]->SetName( Form("hRecoDataDijetEtaForward_%d", i) );
        set1DStyle( hRecoDataDijetEtaForward[i], 2 );
        integralForward = hRecoDataDijetEtaForward[i]->Integral();

        // Eta backward (CM frame)
        hRecoDataDijetEtaBackward[i] = dynamic_cast<TH1D*>( fData->Get( Form("hRecoDijetEtaCMBackward1DWeighted_%d", i) ) );
        if ( !hRecoDataDijetEtaBackward[i] ) {
            std::cerr << Form("Data histogram not found: hRecoDijetEtaCMBackward1DWeighted_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaBackward[i]->SetName( Form("hRecoDataDijetEtaBackward_%d", i) );
        set1DStyle( hRecoDataDijetEtaBackward[i], 2 );
        integralBackward = hRecoDataDijetEtaBackward[i]->Integral();

        // Normalize the forward and backward distributions
        hRecoDataDijetEtaForward[i]->Scale( integralForward / (integralForward + integralBackward) );
        hRecoDataDijetEtaBackward[i]->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hRecoDataDijetFBRatio[i] = dynamic_cast<TH1D*>( hRecoDataDijetEtaForward[i]->Clone( Form("hRecoDataDijetFBRatio_%d", i) ) );
        hRecoDataDijetFBRatio[i]->SetName( Form("hRecoDataDijetFBRatio_%d", i) );
        hRecoDataDijetFBRatio[i]->Divide( hRecoDataDijetEtaBackward[i] );

        //
        // MC (reco)
        //

        // Eta full (lab frame)
        hRecoMcDijetEta[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRecoDijetEta1DWeighted_%d", i) ) );
        if ( !hRecoMcDijetEta[i] ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEta1DWeighted_%d", i) << std::endl; return;
        }
        hRecoMcDijetEta[i]->SetName( Form("hRecoMcDijetEta_%d", i) );
        set1DStyle( hRecoMcDijetEta[i], 0 );
        rescaleEta( hRecoMcDijetEta[i] );

        // Eta full (CM frame)
        hRecoMcDijetEtaCM[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRecoDijetEta1DCMWeighted_%d", i) ) );
        if ( !hRecoMcDijetEtaCM[i] ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEta1DCMWeighted_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaCM[i]->SetName( Form("hRecoMcDijetEtaCM_%d", i) );
        set1DStyle( hRecoMcDijetEtaCM[i], 0 );
        rescaleEta( hRecoMcDijetEtaCM[i] );

        // Eta forward (CM frame)
        hRecoMcDijetEtaForward[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRecoDijetEtaCMForward1DWeighted_%d", i) ) );
        if ( !hRecoMcDijetEtaForward[i] ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEtaCMForward1DWeighted_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaForward[i]->SetName( Form("hRecoMcDijetEtaForward_%d", i) );
        set1DStyle( hRecoMcDijetEtaForward[i], 0 );
        integralForward = hRecoMcDijetEtaForward[i]->Integral();

        // Eta backward (CM frame)
        hRecoMcDijetEtaBackward[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRecoDijetEtaCMBackward1DWeighted_%d", i) ) );
        if ( !hRecoMcDijetEtaBackward[i] ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEtaCMBackward1DWeighted_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaBackward[i]->SetName( Form("hRecoMcDijetEtaBackward_%d", i) );
        set1DStyle( hRecoMcDijetEtaBackward[i], 0 );
        integralBackward = hRecoMcDijetEtaBackward[i]->Integral();

        // Normalize the forward and backward distributions
        hRecoMcDijetEtaForward[i]->Scale( integralForward / (integralForward + integralBackward) );
        hRecoMcDijetEtaBackward[i]->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hRecoMcDijetFBRatio[i] = dynamic_cast<TH1D*>( hRecoMcDijetEtaForward[i]->Clone( Form("hRecoMcDijetFBRatio_%d", i) ) );
        hRecoMcDijetFBRatio[i]->SetName( Form("hRecoMcDijetFBRatio_%d", i) );
        hRecoMcDijetFBRatio[i]->Divide( hRecoMcDijetEtaBackward[i] );

        //
        // MC (gen)
        //

        // Eta full (lab frame)
        hGenMcDijetEta[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hGenDijetEta1DWeighted_%d", i) ) );
        if ( !hGenMcDijetEta[i] ) {
            std::cerr << Form("MC histogram not found: hGenDijetEta1DWeighted_%d", i) << std::endl; return;
        }
        hGenMcDijetEta[i]->SetName( Form("hGenMcDijetEta_%d", i) );
        set1DStyle( hGenMcDijetEta[i], 5 );
        rescaleEta( hGenMcDijetEta[i] );

        // Eta full (CM frame)
        hGenMcDijetEtaCM[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hGenDijetEta1DCMWeighted_%d", i) ) );
        if ( !hGenMcDijetEtaCM[i] ) {
            std::cerr << Form("MC histogram not found: hGenDijetEta1DCMWeighted_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaCM[i]->SetName( Form("hGenMcDijetEtaCM_%d", i) );
        set1DStyle( hGenMcDijetEtaCM[i], 5 );
        rescaleEta( hGenMcDijetEtaCM[i] );

        // Eta forward (CM frame)
        hGenMcDijetEtaForward[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hGenDijetEtaCMForward1DWeighted_%d", i) ) );
        if ( !hGenMcDijetEtaForward[i] ) {
            std::cerr << Form("MC histogram not found: hGenDijetEtaCMForward1DWeighted_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaForward[i]->SetName( Form("hGenMcDijetEtaForward_%d", i) );
        set1DStyle( hGenMcDijetEtaForward[i], 5 );
        integralForward = hGenMcDijetEtaForward[i]->Integral();

        // Eta backward (CM frame)
        hGenMcDijetEtaBackward[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hGenDijetEtaCMBackward1DWeighted_%d", i) ) );
        if ( !hGenMcDijetEtaBackward[i] ) {
            std::cerr << Form("MC histogram not found: hGenDijetEtaCMBackward1DWeighted_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaBackward[i]->SetName( Form("hGenMcDijetEtaBackward_%d", i) );
        set1DStyle( hGenMcDijetEtaBackward[i], 5 );
        integralBackward = hGenMcDijetEtaBackward[i]->Integral();

        // Normalize the forward and backward distributions
        hGenMcDijetEtaForward[i]->Scale( integralForward / (integralForward + integralBackward) );
        hGenMcDijetEtaBackward[i]->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hGenMcDijetFBRatio[i] = dynamic_cast<TH1D*>( hGenMcDijetEtaForward[i]->Clone( Form("hGenMcDijetFBRatio_%d", i) ) );
        hGenMcDijetFBRatio[i]->SetName( Form("hGenMcDijetFBRatio_%d", i) );
        hGenMcDijetFBRatio[i]->Divide( hGenMcDijetEtaBackward[i] );

        //
        // MC (ref)
        //

        // Eta full (lab frame)
        hRefMcDijetEta[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRefDijetEta1DWeighted_%d", i) ) );
        if ( !hRefMcDijetEta[i] ) {
            std::cerr << Form("MC histogram not found: hRefDijetEta1DWeighted_%d", i) << std::endl; return;
        }
        hRefMcDijetEta[i]->SetName( Form("hRefMcDijetEta_%d", i) );
        set1DStyle( hRefMcDijetEta[i], 1 );
        rescaleEta( hRefMcDijetEta[i] );

        // Eta full (CM frame)
        hRefMcDijetEtaCM[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRefDijetEta1DCMWeighted_%d", i) ) );
        if ( !hRefMcDijetEtaCM[i] ) {
            std::cerr << Form("MC histogram not found: hRefDijetEta1DCMWeighted_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaCM[i]->SetName( Form("hRefMcDijetEtaCM_%d", i) );
        set1DStyle( hRefMcDijetEtaCM[i], 1 );
        rescaleEta( hRefMcDijetEtaCM[i] );

        // Eta forward (CM frame)
        hRefMcDijetEtaForward[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRefDijetEtaCMForward1DWeighted_%d", i) ) );
        if ( !hRefMcDijetEtaForward[i] ) {
            std::cerr << Form("MC histogram not found: hRefDijetEtaCMForward1DWeighted_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaForward[i]->SetName( Form("hRefMcDijetEtaForward_%d", i) );
        set1DStyle( hRefMcDijetEtaForward[i], 1 );
        integralForward = hRefMcDijetEtaForward[i]->Integral();

        // Eta backward (CM frame)
        hRefMcDijetEtaBackward[i] = dynamic_cast<TH1D*>( fMc->Get( Form("hRefDijetEtaCMBackward1DWeighted_%d", i) ) );
        if ( !hRefMcDijetEtaBackward[i] ) {
            std::cerr << Form("MC histogram not found: hRefDijetEtaCMBackward1DWeighted_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaBackward[i]->SetName( Form("hRefMcDijetEtaBackward_%d", i) );
        set1DStyle( hRefMcDijetEtaBackward[i], 1 );
        integralBackward = hRefMcDijetEtaBackward[i]->Integral();

        // Normalize the forward and backward distributions
        hRefMcDijetEtaForward[i]->Scale( integralForward / (integralForward + integralBackward) );
        hRefMcDijetEtaBackward[i]->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hRefMcDijetFBRatio[i] = dynamic_cast<TH1D*>( hRefMcDijetEtaForward[i]->Clone( Form("hRefMcDijetFBRatio_%d", i) ) );
        hRefMcDijetFBRatio[i]->SetName( Form("hRefMcDijetFBRatio_%d", i) );
        hRefMcDijetFBRatio[i]->Divide( hRefMcDijetEtaBackward[i] );

        //
        // Ratios of full distributions to gen 
        //

        hRecoData2GenDijetEta[i] = dynamic_cast<TH1D*>( hRecoDataDijetEta[i]->Clone( Form("hRecoData2GenDijetEta_%d", i) ) );
        hRecoData2GenDijetEta[i]->SetName( Form("hRecoData2GenDijetEta_%d", i) );
        hRecoData2GenDijetEta[i]->Divide( hRecoData2GenDijetEta[i], hGenMcDijetEta[i], 1., 1. );

        hRecoMc2GenDijetEta[i] = dynamic_cast<TH1D*>( hRecoMcDijetEta[i]->Clone( Form("hRecoMc2GenDijetEta_%d", i) ) );
        hRecoMc2GenDijetEta[i]->SetName( Form("hRecoMc2GenDijetEta_%d", i) );
        hRecoMc2GenDijetEta[i]->Divide( hRecoMc2GenDijetEta[i], hGenMcDijetEta[i], 1., 1., "b" );

        hRef2GenDijetEta[i] = dynamic_cast<TH1D*>( hRefMcDijetEta[i]->Clone( Form("hRef2GenDijetEta_%d", i) ) );
        hRef2GenDijetEta[i]->SetName( Form("hRef2GenDijetEta_%d", i) );
        hRef2GenDijetEta[i]->Divide( hRef2GenDijetEta[i], hGenMcDijetEta[i], 1., 1., "b" );

        hRecoData2GenDijetEtaCM[i] = dynamic_cast<TH1D*>( hRecoDataDijetEtaCM[i]->Clone( Form("hRecoData2GenDijetEtaCM_%d", i) ) );
        hRecoData2GenDijetEtaCM[i]->SetName( Form("hRecoData2GenDijetEtaCM_%d", i) );
        hRecoData2GenDijetEtaCM[i]->Divide( hRecoData2GenDijetEtaCM[i], hGenMcDijetEtaCM[i], 1., 1. );

        hRecoMc2GenDijetEtaCM[i] = dynamic_cast<TH1D*>( hRecoMcDijetEtaCM[i]->Clone( Form("hRecoMc2GenDijetEtaCM_%d", i) ) );
        hRecoMc2GenDijetEtaCM[i]->SetName( Form("hRecoMc2GenDijetEtaCM_%d", i) );
        hRecoMc2GenDijetEtaCM[i]->Divide( hRecoMc2GenDijetEtaCM[i], hGenMcDijetEtaCM[i], 1., 1., "b" );

        hRef2GenDijetEtaCM[i] = dynamic_cast<TH1D*>( hRefMcDijetEtaCM[i]->Clone( Form("hRef2GenDijetEtaCM_%d", i) ) );
        hRef2GenDijetEtaCM[i]->SetName( Form("hRef2GenDijetEtaCM_%d", i) );
        hRef2GenDijetEtaCM[i]->Divide( hRef2GenDijetEtaCM[i], hGenMcDijetEtaCM[i], 1., 1., "b" );
    } // for (int i = 0; i < dijetPtBins; ++i)

    //
    // Plot comparisons and ratios
    //

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );
    TLegend *leg;
    TCanvas *c = new TCanvas("c", "c", 1200, 1200);

    //
    // Loop over the dijet pT average bins and plot the distributions
    //
    for (int i = 0; i < dijetPtBins; ++i) {

        TLatex t;

        //
        // Full eta distribution and f/b comparisons 
        //

        // Lab frame
        drawDijetToGenComparison(c, hRecoMcDijetEta[i], hGenMcDijetEta[i], hRefMcDijetEta[i], nullptr, hRecoDataDijetEta[i], 
                            (int)dijetPtVals[i], (int)dijetPtVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaLab_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        c->cd();
        setPadStyle();
        hRecoDataDijetEta[i]->Draw();
        hRecoDataDijetEta[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        hRecoDataDijetEta[i]->GetYaxis()->SetRangeUser(0.0, 0.12);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );
        c->SaveAs( Form("%s/%s_dijetEtaLab_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        // CM frame
        drawDijetToGenComparison(c, hRecoMcDijetEtaCM[i], hGenMcDijetEtaCM[i], hRefMcDijetEtaCM[i], nullptr, hRecoDataDijetEtaCM[i], 
                            (int)dijetPtVals[i], (int)dijetPtVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaCM_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );
        
        c->cd();
        setPadStyle();
        hRecoDataDijetEtaCM[i]->Draw();
        hRecoDataDijetEtaCM[i]->GetXaxis()->SetRangeUser(-2., 2.);
        hRecoDataDijetEtaCM[i]->GetYaxis()->SetRangeUser(0., 0.12);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );
        c->SaveAs( Form("%s/%s_dijetEtaCM_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        // Forward-backward ratio
        drawDijetToGenComparison(c, hRecoMcDijetFBRatio[i], hGenMcDijetFBRatio[i], hRefMcDijetFBRatio[i], nullptr, hRecoDataDijetFBRatio[i], 
                            (int)dijetPtVals[i], (int)dijetPtVals[i+1], true, true, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaFBRatio_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        c->cd();
        setPadStyle();
        hRecoDataDijetFBRatio[i]->Draw();
        hRecoDataDijetFBRatio[i]->GetXaxis()->SetRangeUser(0., 2.);
        hRecoDataDijetFBRatio[i]->GetYaxis()->SetRangeUser(0.8, 1.2);

        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );
        c->SaveAs( Form("%s/%s_dijetEtaFBRatio_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        //
        // Full eta distribution ratios
        //

        // Lab frame
        drawDijetToGenRatio(c, hRecoMc2GenDijetEta[i], hRef2GenDijetEta[i], nullptr, hRecoData2GenDijetEta[i],
                       (int)dijetPtVals[i], (int)dijetPtVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaLab_DataRecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        // CM frame
        drawDijetToGenRatio(c, hRecoMc2GenDijetEtaCM[i], hRef2GenDijetEtaCM[i], nullptr, hRecoData2GenDijetEtaCM[i],
                       (int)dijetPtVals[i], (int)dijetPtVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaCM_DataRecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );
    }
}

//________________
// Plot eta and pT distributions of inclusive jets from Pythia and embedding.
// Compare dijet distributions.
//
void plotPythia2EmbeddingComparisons(TFile *fEmbedding, TFile *fPythia, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250305") {
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", (int)(1000 * collisionEnergy));
    collSystemStr += Form("_emb2pythia");

    double xTextPosition = 0.6;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Create vector of ptHat and jet pT bins for projections
    int ptHatStart = 15;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins{1, 4, 6, 7, 8, 10, 30, 50};

    // Jet pT binning
    int jetPtStart = 5;
    int jetPtStep = 10; // Starting from 5 GeV: jetPtStart + (jetPtBins(i) - 1) * jetPtStep
    int jetPtBinsMax = 150;
    std::vector<int> jetPtBinsLow {1, 5,  8, 13,  20}; // 
    std::vector<int> jetPtBinsHigh{4, 7, 12, 19, 150}; // 


    // Eta binning
    // 52 bins from (-5.2, 5.2)
    int nEtaBins = 52;
    double etaStep = 0.2;
    std::vector<int> jetEtaBinsLow{19, 12, 38};
    std::vector<int> jetEtaBinsHigh{35, 16, 42};

    const int dijetPtBins{16};
    double dijetPtVals[dijetPtBins+1] {  50.,  60.,   70.,  80.,  90.,
                                         100., 110.,  120., 130., 140.,
                                         150., 160.,  180., 200., 250., 
                                         300., 500.};

    //
    // Retrieve inclusive jet histograms
    //

    // Embedding
    TH3D *hRecoPtEtaPtHatEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRecoInclusiveJetPtEtaPtHat"));
    hRecoPtEtaPtHatEmb->SetName("hRecoPtEtaPtHatEmb");
    TH3D *hGenPtEtaPtHatEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hGenInclusiveJetPtEtaPtHat"));
    hGenPtEtaPtHatEmb->SetName("hGenPtEtaPtHatEmb");
    TH1D *hRecoEtaEmb[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hGenEtaEmb[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hRecoPtEmb[ptHatBins.size()][jetEtaBinsLow.size()];
    TH1D *hGenPtEmb[ptHatBins.size()][jetEtaBinsLow.size()];

    // Pythia
    TH3D *hRecoPtEtaPtHatPythia = dynamic_cast<TH3D *>(fPythia->Get("hRecoInclusiveJetPtEtaPtHat"));
    hRecoPtEtaPtHatPythia->SetName("hRecoPtEtaPtHatPythia");
    TH3D *hGenPtEtaPtHatPythia = dynamic_cast<TH3D *>(fPythia->Get("hGenInclusiveJetPtEtaPtHat"));
    hGenPtEtaPtHatPythia->SetName("hGenPtEtaPtHatPythia");
    TH1D *hRecoEtaPythia[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hGenEtaPythia[ptHatBins.size()][jetPtBinsLow.size()];
    TH1D *hRecoPtPythia[ptHatBins.size()][jetEtaBinsLow.size()];
    TH1D *hGenPtPythia[ptHatBins.size()][jetEtaBinsLow.size()];

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();

    // Loop over ptHat bins
    for (size_t iPtHat = 0; iPtHat < ptHatBins.size(); ++iPtHat) {

        int ptHatLow = (int)hRecoPtEtaPtHatEmb->GetZaxis()->GetBinLowEdge(ptHatBins[iPtHat]);

        // Loop over ptHat bins and make projections of eta
        for (size_t iJetPt = 0; iJetPt < jetPtBinsLow.size(); ++iJetPt) {

            int ptLow = (int)hRecoPtEtaPtHatEmb->GetYaxis()->GetBinLowEdge(jetPtBinsLow[iJetPt]);
            int ptHigh = (int)hRecoPtEtaPtHatEmb->GetYaxis()->GetBinUpEdge(jetPtBinsHigh[iJetPt]);

            // Embedding
            hRecoEtaEmb[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hRecoPtEtaPtHatEmb->ProjectionX(Form("hRecoEtaEmb_%zu_%zu", iPtHat, iJetPt), 
                                                                                               jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleEta(hRecoEtaEmb[iPtHat][iJetPt]);
            set1DStyle(hRecoEtaEmb[iPtHat][iJetPt], 0);
            hGenEtaEmb[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hGenPtEtaPtHatEmb->ProjectionX(Form("hGenEtaEmb_%zu_%zu", iPtHat, iJetPt), 
                                                                                             jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleEta(hGenEtaEmb[iPtHat][iJetPt]);
            set1DStyle(hGenEtaEmb[iPtHat][iJetPt], 0);

            // Pythia
            hRecoEtaPythia[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hRecoPtEtaPtHatPythia->ProjectionX(Form("hRecoEtaPythia_%zu_%zu", iPtHat, iJetPt), 
                                                                                                   jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleEta(hRecoEtaPythia[iPtHat][iJetPt]);
            set1DStyle(hRecoEtaPythia[iPtHat][iJetPt], 1);
            hGenEtaPythia[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hGenPtEtaPtHatPythia->ProjectionX(Form("hGenEtaPythia_%zu_%zu", iPtHat, iJetPt), 
                                                                                                 jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleEta(hGenEtaPythia[iPtHat][iJetPt]);
            set1DStyle(hGenEtaPythia[iPtHat][iJetPt], 1);

            // Plot comparisons
            c->Clear();
            setPadStyle();
            drawEtaPtComparison(c, hRecoEtaEmb[iPtHat][iJetPt], hRecoEtaPythia[iPtHat][iJetPt], 
                                ptLow, ptHigh, -5.2, 5.2, ptHatLow, 
                                "Reco (embed)", "Reco (pythia)", 1, 8.16, 
                                false, true, false);
            c->SaveAs(Form("%s/%s_reco_eta_ptHat_ptHat%d_jetPt_%d_%d.pdf", date.Data(), collSystemStr.Data(), ptHatLow, ptLow, ptHigh));

            c->Clear();
            setPadStyle();
            drawEtaPtComparison(c, hGenEtaEmb[iPtHat][iJetPt], hGenEtaPythia[iPtHat][iJetPt], 
                                ptLow, ptHigh, -5.2, 5.2, ptHatLow, 
                                "Gen (embed)", "Gen (pythia)", 1, 8.16, 
                                false, true, false);
                                c->SaveAs(Form("%s/%s_gen_eta_ptHat_ptHat%d_jetPt_%d_%d.pdf", date.Data(), collSystemStr.Data(), ptHatLow, ptLow, ptHigh));
            
        } // for (size_t iJetPt = 0; iJetPt < jetPtBinsLow.size(); ++iJetPt)

    } // for (size_t iPtHat = 0; iPtHat < ptHatBins.size(); ++iPtHat)

    //
    // Retrieve dijet histograms
    //

    // Embedding
    TH1D *hGenDijet1DCMEmb[dijetPtBins];
    for (int i = 0; i < dijetPtBins; ++i) {
        hGenDijet1DCMEmb[i] = dynamic_cast<TH1D *>(fEmbedding->Get(Form("hGenDijetEta1DWeighted_%d", i)));
        hGenDijet1DCMEmb[i]->SetName(Form("hGenDijet1DCMEmb_%d", i));
    }

    // Pythia
    TH1D *hGenDijet1DCMPythia[dijetPtBins];
    for (int i = 0; i < dijetPtBins; ++i) {
        hGenDijet1DCMPythia[i] = dynamic_cast<TH1D *>(fPythia->Get(Form("hGenDijetEta1DWeighted_%d", i)));
        hGenDijet1DCMPythia[i]->SetName(Form("hGenDijet1DCMPythia_%d", i));
    }


}

//________________
void plotMcClosures() {

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

    int collisionSystem = 0;  // 0 - pp, 1 - pPb, 2 - pPb5020, 3 - pPb8160
    double collisionEnergy = 5.02;  // 5.02 TeV

    //
    // pPb8160
    //

    collisionSystem = 1;
    collisionEnergy = 8.16;
    int direction = 1; // 0-p-going, 1-Pb-going, 2 - combined // TODO: implement combined
    TString directionStr = (direction == 0) ? "pgoing" : "Pbgoing";


    // MC p-going direction new (coincides with the pPb5020)
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_eta20.root", uname.Data(), directionStr.Data(), directionStr.Data()) );
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_eta20.root", uname.Data()) );
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/build/oTest_pPb8160_dijet_ptHat_50_80_noTrkMax.root", uname.Data()) );
    TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/macro/pPb8160_Pbgoing_noTrkMax/oEmbedding_pPb8160_Pbgoing_noTrkMax.root", uname.Data()) );
    if ( !pPb8160EmbedFile ) {
        std::cerr << Form("File not found: /Users/%s/work/cms/soft/jetAnalysis/macro/pPb8160_Pbgoing_noTrkMax/oEmbedding_pPb8160_Pbgoing_noTrkMax.root", uname.Data()) << std::endl;
        return;
    }

    TFile *pPb8160PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_eta20.root", uname.Data(), directionStr.Data(), directionStr.Data()) );
    // TFile *pPb8160PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_eta20.root", uname.Data()) );
    if ( !pPb8160PythiaFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/pythia/oPythia_%s_def_ak4_eta20.root", uname.Data(), directionStr.Data()) << std::endl;
        return;
    }

    // TFile *pPb8160DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_jerDef_ak4_eta25.root", uname.Data(), directionStr.Data(), directionStr.Data()) );
    // TFile *pPb8160DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/exp/%s/PAEGJet60_%s_ak4_eta20.root", uname.Data(), directionStr.Data(), directionStr.Data()) );
    TFile *pPb8160DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4_eta20.root", uname.Data()) );
    if ( !pPb8160DataFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/exp/%s/PAEGJet60_%s_ak4_eta20.root", uname.Data(), directionStr.Data(), directionStr.Data()) << std::endl;
        return;
    }

    //
    // Plot for inclusive jets JEC closures (scan in eta and pT)
    //
    plotInclusiveJetJECClosures(pPb8160EmbedFile, collisionSystem, collisionEnergy, date);

    //
    // Comparison of dijet reco and ref to gen distributions
    //
    plotDijetClosures( pPb8160EmbedFile, collisionSystem, collisionEnergy, date );


    //
    // Plot comparison of inclusive jet eta distributions to check/validate the JEC
    //
    // plotData2McInclusiveJetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of dijet reco and ref to gen distributions
    //
    // plotData2McDijetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of inclusive jets and dijets for embedding and PYTHIA
    //
    // plotPythia2EmbeddingComparisons(pPb8160EmbedFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);
}
