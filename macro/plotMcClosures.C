// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
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
// Rescale a 1D histogram by its area-weighted integral
// This function rescales a 1D histogram by its area-weighted integral.
// It calculates the total area-weighted sum of the histogram and then normalizes each bin content and error accordingly.
// The function assumes that the histogram has non-equidistant axes.
void rescaleHisto1D(TH1* h1) {
    if (!h1) {
        std::cerr << "Null histogram pointer!" << std::endl;
        return;
    }

    const int nBins = h1->GetNbinsX();
    TAxis* xAxis = h1->GetXaxis();

    double total = 0.0;

    // First pass: calculate total sum (content × bin width)
    for (int i = 1; i <= nBins; ++i) {
        double width = xAxis->GetBinWidth(i);
        double content = h1->GetBinContent(i);
        total += content * width;
    }

    // Second pass: normalize content and error
    if (total > 0) {
        for (int i = 1; i <= nBins; ++i) {
            double content = h1->GetBinContent(i);
            double error   = h1->GetBinError(i);

            h1->SetBinContent(i, content / total);
            h1->SetBinError(i, error / total);
        }
    } 
    else {
        std::cerr << "Warning: total area-normalized sum is zero!" << std::endl;
    }
    h1->Scale(1.0 / h1->Integral()); // Normalize to unity
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

    double xRange[2] = { -3., 3. };
    // double xRange[2] = { 0., 3. };
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
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T} < %d GeV", ptLow, ptHi) );
    t.DrawLatexNDC(0.35, 0.78, Form("#hat{p}_{T} > %d GeV", ptHatLow));
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

    double xRange[2] = {-3.0, 3.0};
    // double xRange[2] = {0.0, 3.0};
    double yRange[2] = {0.75, 1.25};

    if (isCM) {
        if (!isFB) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] = 0.75; yRange[1] = 1.25;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.5;
            yRange[0] = 0.8; yRange[1] = 1.2;
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
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T} < %d GeV", ptLow, ptHi));
    t.DrawLatexNDC(0.35, 0.78, Form("#hat{p}_{T} > %d GeV", ptHatLow));
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

    double xRange[2] = { -3., 3. };
    double yRange[2] = {0.0000001, 0.12 };

    if ( isCM ) {
        if ( !isFB ) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] = 0.0000001; yRange[1] = 0.12;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.4;
            yRange[0] = 0.8; yRange[1] = 1.2;
        }
    }

    c->cd();
    setPadStyle();
    hReco->Draw();
    if ( hGen ) {
        hGen->Draw("same");
    }
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
    // gPad->SetLogy(); // Set log scale for y-axis
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);        
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T}^{ave} GeV < %d", ptLow, ptHi) );
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
    if ( hGen ) {
        leg->AddEntry( hGen, "Gen", "p" );
    }
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

    double xRange[2] = {-3.0, 3.0};
    double yRange[2] = {0.85, 1.15};

    if (isCM) {
        if (!isFB) {
            // xRange[0] = -3.0; xRange[1] = 3.0;
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] = 0.85; yRange[1] = 1.15;
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
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
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
// Function calculates forward, backward and ratio of forward/backward distribution from
// the full 1D histograms in the CM frame
void recalculateFBRatioFromFullDistribution(TH1* hRatio, TH1* hForward, TH1* hBackward, TH1* hFull) {
    if (!hRatio || !hForward || !hBackward || !hFull) {
        std::cerr << "Error: Null histogram pointer passed." << std::endl;
        return;
    }

    if (hRatio->GetNbinsX() != hForward->GetNbinsX() || 
        hForward->GetNbinsX() != hBackward->GetNbinsX()) {
        std::cerr << "Error: Histograms have different number of bins." << std::endl;
        return;
    }

    // std::cout << Form("Number of bins Full: %d Forward: %d Backward: %d Ratio: %d\n",
    //                   hFull->GetNbinsX(), hForward->GetNbinsX(), hBackward->GetNbinsX(), hRatio->GetNbinsX());

    int middleBin = hFull->FindBin(0.);

    // Make forward
    for (int i=middleBin; i<=hFull->GetNbinsX(); i++) {
        double binVal = hFull->GetBinContent(i);
        double binError = hFull->GetBinError(i);
        hForward->SetBinContent(i-middleBin+1, binVal);
        hForward->SetBinError(i-middleBin+1, binError);
    } // for (int i=middleBin; i<=hFull->GetNbinsX(); i++)

    // Make backward
    for (int i=middleBin-1; i>0; i--) {
        double binVal = hFull->GetBinContent(i);
        double binError = hFull->GetBinError(i);
        hBackward->SetBinContent(middleBin-i, binVal);
        hBackward->SetBinError(middleBin-i, binError);
    } // for (int i=1; i<middleBin; i++)

    hRatio->Divide(hForward, hBackward, 1., 1.);
    for (int i=1; i<=hRatio->GetNbinsX(); i++) {
        double binVal = hRatio->GetBinContent(i);
        double binError = hRatio->GetBinError(i);
    }
}

//________________
// Plot dijet eta comparison of reco, ref and refSel to gen
void dijetClosures(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    
    // collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // Dijet ptAve binning
    double dijetPtVals[] {  50.,  60.,   70.,  80.,  90.,
                            100., 110.,  120., 130., 140.,
                            150., 160.,  180., 200., 250., 
                            300., 500.};
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    TH1D *h1DProj{nullptr};
    TH1D *h1DDirect{nullptr};

    //
    // Lab frame
    //
    TH1D *hRecoDijetEta1DLab{nullptr};
    TH1D *hGenDijetEta1DLab{nullptr};
    TH1D *hRefDijetEta1DLab{nullptr};
    TH1D *hRefSelDijetEta1DLab{nullptr};
    TH1D *hReco2GenDijetEta1DLab{nullptr};
    TH1D *hRef2GenDijetEta1DLab{nullptr};
    TH1D *hRefSel2GenDijetEta1DLab{nullptr};

    //
    // CM frame
    //
    TH1D *hRecoDijetEta1DCM{nullptr};
    TH1D *hGenDijetEta1DCM{nullptr};
    TH1D *hRefDijetEta1DCM{nullptr};
    TH1D *hRefSelDijetEta1DCM{nullptr};
    TH1D *hReco2GenDijetEta1DCM{nullptr};
    TH1D *hRef2GenDijetEta1DCM{nullptr};
    TH1D *hRefSel2GenDijetEta1DCM{nullptr};

    //
    // Forward/backward ratios
    //

    TH1D *hRecoDijetEtaCMForward1D{nullptr};
    TH1D *hRecoDijetEtaCMBackward1D{nullptr};
    TH1D *hGenDijetEtaCMForward1D{nullptr};
    TH1D *hGenDijetEtaCMBackward1D{nullptr};
    TH1D *hRefDijetEtaCMForward1D{nullptr};
    TH1D *hRefDijetEtaCMBackward1D{nullptr};

    TH1D *hRecoDijetFBEtaCM1D{nullptr};
    TH1D *hGenDijetFBEtaCM1D{nullptr};
    TH1D *hRefDijetFBEtaCM1D{nullptr};

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );
    TLegend *leg;
    TCanvas *c = new TCanvas( "c", "c", 1000, 1000 );

    // Loop over dijet ptAve bins
    for (int i = 0; i < sizeOfPtVals-1; i++) {

        double ptLow = dijetPtVals[i];
        double ptHi = dijetPtVals[i+1];

        //
        // Lab frame
        //
        hRecoDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1D_%d", i) ) );
        if (!hRecoDijetEta1DLab) {
            std::cout << Form("Failed to retrieve hRecoDijetEta1D_%d", i) << std::endl;
        }
        hRecoDijetEta1DLab->SetName( Form("hRecoDijetEta1DLab_%d", i) );
        rescaleHisto1D( hRecoDijetEta1DLab );
        set1DStyle( hRecoDijetEta1DLab, 0 );
        
        hGenDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1D_%d", i) ) );
        if (!hGenDijetEta1DLab) {
            std::cout << Form("Failed to retrieve hGenDijetEta1D_%d", i) << std::endl;
        }
        hGenDijetEta1DLab->SetName( Form("hGenDijetEta1DLab_%d", i) );
        rescaleHisto1D( hGenDijetEta1DLab );
        set1DStyle( hGenDijetEta1DLab, 4 );

        hRefDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1D_%d", i) ) );
        if (!hRefDijetEta1DLab) {
            std::cout << Form("Failed to retrieve hRefDijetEta1D_%d", i) << std::endl;
        }
        hRefDijetEta1DLab->SetName( Form("hRefDijetEta1DLab_%d", i) );
        rescaleHisto1D( hRefDijetEta1DLab );
        set1DStyle( hRefDijetEta1DLab, 1 );

        hRefSelDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hRefSelDijetEta1D_%d", i) ) );
        if (!hRefSelDijetEta1DLab) {
            std::cout << Form("Failed to retrieve hRefSelDijetEta1D_%d", i) << std::endl;
        }
        hRefSelDijetEta1DLab->SetName( Form("hRefSelDijetEta1DLab_%d", i) );
        rescaleHisto1D( hRefSelDijetEta1DLab );
        set1DStyle( hRefSelDijetEta1DLab, 2 );

        int rebinFactor = 2;
        hRecoDijetEta1DLab->Rebin( rebinFactor );
        hGenDijetEta1DLab->Rebin( rebinFactor );
        hRefDijetEta1DLab->Rebin( rebinFactor );
        hRefSelDijetEta1DLab->Rebin( rebinFactor );

        hReco2GenDijetEta1DLab = dynamic_cast<TH1D*>( hRecoDijetEta1DLab->Clone( Form("hReco2GenDijetEta1DLab_%d", i) ) );
        hReco2GenDijetEta1DLab->Divide( hReco2GenDijetEta1DLab, hGenDijetEta1DLab, 1., 1., "b" );
        hRef2GenDijetEta1DLab = dynamic_cast<TH1D*>( hRefDijetEta1DLab->Clone( Form("hRef2GenDijetEta1DLab_%d", i) ) );
        hRef2GenDijetEta1DLab->Divide( hRef2GenDijetEta1DLab, hGenDijetEta1DLab, 1., 1., "b" );
        hRefSel2GenDijetEta1DLab = dynamic_cast<TH1D*>( hRefSelDijetEta1DLab->Clone( Form("hRefSel2GenDijetEta1DLab_%d", i) ) );
        hRefSel2GenDijetEta1DLab->Divide( hRefSel2GenDijetEta1DLab, hGenDijetEta1DLab, 1., 1., "b" );

        //
        // CM frame
        //
        hRecoDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1DCM_%d", i) ) );
        if (!hRecoDijetEta1DCM) {
            std::cout << Form("Failed to retrieve hRecoDijetEta1DCM_%d", i) << std::endl;
        }
        hRecoDijetEta1DCM->SetName( Form("hRecoDijetEta1DCM_%d", i) );
        set1DStyle( hRecoDijetEta1DCM, 0 );
        rescaleHisto1D( hRecoDijetEta1DCM );

        hGenDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1DCM_%d", i) ) );
        if (!hGenDijetEta1DCM) {
            std::cout << Form("Failed to retrieve hGenDijetEta1DCM_%d", i) << std::endl;
        }
        hGenDijetEta1DCM->SetName( Form("hGenDijetEta1DCM_%d", i) );
        set1DStyle( hGenDijetEta1DCM, 4 );
        rescaleHisto1D( hGenDijetEta1DCM );

        hRefDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1DCM_%d", i) ) );
        if (!hRefDijetEta1DCM) {
            std::cout << Form("Failed to retrieve hRefDijetEta1DCM_%d", i) << std::endl;
        }
        hRefDijetEta1DCM->SetName( Form("hRefDijetEta1DCM_%d", i) );
        set1DStyle( hRefDijetEta1DCM, 1 );
        rescaleHisto1D( hRefDijetEta1DCM );

        hRefSelDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hRefSelDijetEta1DCM_%d", i) ) );
        if (!hRefSelDijetEta1DCM) {
            std::cout << Form("Failed to retrieve hRefSelDijetEta1DCM_%d", i) << std::endl;
        }
        hRefSelDijetEta1DCM->SetName( Form("hRefSelDijetEta1DCM_%d", i) );
        set1DStyle( hRefSelDijetEta1DCM, 2 );
        rescaleHisto1D( hRefSelDijetEta1DCM );

        hRecoDijetEta1DCM->Rebin( rebinFactor );
        hGenDijetEta1DCM->Rebin( rebinFactor );
        hRefDijetEta1DCM->Rebin( rebinFactor );
        hRefSelDijetEta1DCM->Rebin( rebinFactor );

        hReco2GenDijetEta1DCM = dynamic_cast<TH1D*>( hRecoDijetEta1DCM->Clone( Form("hReco2GenDijetEta1DCM_%d", i) ) );
        hReco2GenDijetEta1DCM->Divide( hReco2GenDijetEta1DCM, hGenDijetEta1DCM, 1., 1., "b" );
        hRef2GenDijetEta1DCM = dynamic_cast<TH1D*>( hRefDijetEta1DCM->Clone( Form("hRef2GenDijetEta1DCM_%d", i) ) );
        hRef2GenDijetEta1DCM->Divide( hRef2GenDijetEta1DCM, hGenDijetEta1DCM, 1., 1., "b" );
        hRefSel2GenDijetEta1DCM = dynamic_cast<TH1D*>( hRefSelDijetEta1DCM->Clone( Form("hRefSel2GenDijetEta1DCM_%d", i) ) );
        hRefSel2GenDijetEta1DCM->Divide( hRefSel2GenDijetEta1DCM, hGenDijetEta1DCM, 1., 1., "b" );

        //
        // Forward/backward ratios
        //
        hRecoDijetEtaCMForward1D = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMForward1D_%d", i) ) );
        if (!hRecoDijetEtaCMForward1D) {
            std::cout << Form("Failed to retrieve hRecoDijetEtaCMForward1D_%d", i) << std::endl;
        }
        hRecoDijetEtaCMForward1D->SetName( Form("hRecoDijetEtaCMForward1D_%d", i) );
        hRecoDijetEtaCMBackward1D = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMBackward1D_%d", i) ) );
        if (!hRecoDijetEtaCMBackward1D) {
            std::cout << Form("Failed to retrieve hRecoDijetEtaCMBackward1D_%d", i) << std::endl;
        }
        hRecoDijetEtaCMBackward1D->SetName( Form("hRecoDijetEtaCMBackward1D_%d", i) );
        hGenDijetEtaCMForward1D = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMForward1D_%d", i) ) );
        if (!hGenDijetEtaCMForward1D) {
            std::cout << Form("Failed to retrieve hGenDijetEtaCMForward1D_%d", i) << std::endl;
        }
        hGenDijetEtaCMForward1D->SetName( Form("hGenDijetEtaCMForward1D_%d", i) );
        hGenDijetEtaCMBackward1D = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMBackward1D_%d", i) ) );
        if (!hGenDijetEtaCMBackward1D) {
            std::cout << Form("Failed to retrieve hGenDijetEtaCMBackward1D_%d", i) << std::endl;
        }
        hGenDijetEtaCMBackward1D->SetName( Form("hGenDijetEtaCMBackward1D_%d", i) );
        hRefDijetEtaCMForward1D = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMForward1D_%d", i) ) );
        if (!hRefDijetEtaCMForward1D) {
            std::cout << Form("Failed to retrieve hRefDijetEtaCMForward1D_%d", i) << std::endl;
        }
        hRefDijetEtaCMForward1D->SetName( Form("hRefDijetEtaCMForward1D_%d", i) );
        hRefDijetEtaCMBackward1D = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMBackward1D_%d", i) ) );
        if (!hRefDijetEtaCMBackward1D) {
            std::cout << Form("Failed to retrieve hRefDijetEtaCMBackward1D_%d", i) << std::endl;
        }
        hRefDijetEtaCMBackward1D->SetName( Form("hRefDijetEtaCMBackward1D_%d", i) );

        // rescaleForwardBackward( hRecoDijetEtaCMForward1D, hRecoDijetEtaCMBackward1D );
        // rescaleForwardBackward( hGenDijetEtaCMForward1D, hGenDijetEtaCMBackward1D );
        // rescaleForwardBackward( hRefDijetEtaCMForward1D, hRefDijetEtaCMBackward1D );
        set1DStyle( hRecoDijetEtaCMForward1D, 0 );
        set1DStyle( hGenDijetEtaCMForward1D, 2 );
        set1DStyle( hRefDijetEtaCMForward1D, 1 );

        hRecoDijetFBEtaCM1D = dynamic_cast<TH1D*>( hRecoDijetEtaCMForward1D->Clone( Form("hRecoDijetFBEtaCM1D_%d", i) ) );
        hRecoDijetFBEtaCM1D->Divide( hRecoDijetEtaCMBackward1D );
        hRecoDijetFBEtaCM1D->GetYaxis()->SetTitle("Forward/Backward");
        hGenDijetFBEtaCM1D = dynamic_cast<TH1D*>( hGenDijetEtaCMForward1D->Clone( Form("hGenDijetFBEtaCM1D_%d", i) ) );
        hGenDijetFBEtaCM1D->Divide( hGenDijetEtaCMBackward1D );
        hGenDijetFBEtaCM1D->GetYaxis()->SetTitle("Forward/Backward");
        hRefDijetFBEtaCM1D = dynamic_cast<TH1D*>( hRefDijetEtaCMForward1D->Clone( Form("hRefDijetFBEtaCM1D_%d", i) ) );
        hRefDijetFBEtaCM1D->Divide( hRefDijetEtaCMBackward1D );
        hRefDijetFBEtaCM1D->GetYaxis()->SetTitle("Forward/Backward");

        //
        // Plot comparisons in the lab frame
        //
        drawDijetToGenComparison(c, hRecoDijetEta1DLab, hRefDijetEta1DLab, hGenDijetEta1DLab, hRefSelDijetEta1DLab, nullptr,
                                 ptLow, ptHi, false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        (int)ptLow, (int)ptHi) );

        //
        // Plot ratios in the lab frame
        //
        drawDijetToGenRatio(c, hReco2GenDijetEta1DLab, hRef2GenDijetEta1DLab, hRefSel2GenDijetEta1DLab, nullptr,
                            ptLow, ptHi, false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        (int)ptLow, (int)ptHi) );

        //
        // Plot comparisons in the CM frame
        //

        drawDijetToGenComparison(c, hRecoDijetEta1DCM, hRefDijetEta1DCM, hGenDijetEta1DCM, hRefSelDijetEta1DCM, nullptr,
                                 (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        (int)ptLow, (int)ptHi) );

        //
        // Plot ratios in the CM frame
        //
        drawDijetToGenRatio(c, hReco2GenDijetEta1DCM, hRef2GenDijetEta1DCM, hRefSel2GenDijetEta1DCM, nullptr,
                            (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        (int)ptLow, (int)ptHi) );

        //
        // Forward/backward ratios
        //
        drawDijetToGenComparison(c, hRecoDijetFBEtaCM1D, hRefDijetFBEtaCM1D, hGenDijetFBEtaCM1D, nullptr, nullptr,
                                 (int)ptLow, (int)ptHi, true, true, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaFB_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        (int)ptLow, (int)ptHi) );

    } // end loop over dijet pt bins

    if (c) { delete c; c = nullptr; }
}

//________________
// This function retrieves 3D histograms (eta, pT)  for reconstructed, generated, 
// and reference jets, and computes the ratios (to gen).
// f is the input TFile containing the histograms.
// collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
// collisionEnergy: energy in TeV (default is 8.16 TeV for pPb)
// date: date string for saving the plots (default is "20250129")
void dijetClosuresFrom2D(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // Dijet ptAve binning
    // int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
    //                          100, 110,  120, 130, 140,
    //                          150, 160,  180, 200, 250, 
    //                          300, 500 };
    double dijetPtVals[] { 50.,  90.,   120.,  180., 200., 250., 300., 350., 500. };
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    // Laboratory frame
    TH2D *hRecoDijetPtEtaLab = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaWeighted"));
    // TH2D *hRecoDijetPtEtaLab = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaMatched"));
    if (!hRecoDijetPtEtaLab) {
        std::cerr << "Error: hRecoDijetPtEtaWeighted not found in file." << std::endl;
        return;
    }
    hRecoDijetPtEtaLab->SetName("hRecoDijetPtEtaLab");
    TH2D *hGenDijetPtEtaLab = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaWeighted"));
    if (!hGenDijetPtEtaLab) {
        std::cerr << "Error: hGenDijetPtEtaWeighted not found in file." << std::endl;
        return;
    }
    hGenDijetPtEtaLab->SetName("hGenDijetPtEtaLab");
    TH2D *hRefDijetPtEtaLab = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaWeighted"));
    if (!hRefDijetPtEtaLab) {
        std::cerr << "Error: hRefDijetPtEtaWeighted not found in file." << std::endl;
        return;
    }
    hRefDijetPtEtaLab->SetName("hRefDijetPtEtaLab");
    TH2D *hRefSelDijetPtEtaLab = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaWeighted"));
    if (!hRefSelDijetPtEtaLab) {
        std::cerr << "Error: hRefSelDijetPtEtaWeighted not found in file." << std::endl;
        return;
    }
    hRefSelDijetPtEtaLab->SetName("hRefSelDijetPtEtaLab");
    TH2D *hRecoDijetPtEtaCMInLab = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMInLab"));
    if (!hRecoDijetPtEtaCMInLab) {
        std::cerr << "Error: hRecoDijetPtEtaCMInLab not found in file." << std::endl;
        return;
    }
    hRecoDijetPtEtaCMInLab->SetName("hRecoDijetPtEtaCMInLab");
    TH2D *hGenDijetPtEtaCMInLab = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaCMInLab"));
    if (!hGenDijetPtEtaCMInLab) {
        std::cerr << "Error: hGenDijetPtEtaCMInLab not found in file." << std::endl;
        return;
    }
    hGenDijetPtEtaCMInLab->SetName("hGenDijetPtEtaCMInLab");
    TH2D *hRefDijetPtEtaCMInLab = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaCMInLab"));
    if (!hRefDijetPtEtaCMInLab) {
        std::cerr << "Error: hRefDijetPtEtaCMInLab not found in file." << std::endl;
        return;
    }
    hRefDijetPtEtaCMInLab->SetName("hRefDijetPtEtaCMInLab");
    TH2D *hRefSelDijetPtEtaCMInLab = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaCMInLab"));
    if (!hRefSelDijetPtEtaCMInLab) {
        std::cerr << "Error: hRefSelDijetPtEtaCMInLab not found in file." << std::endl;
        return;
    }
    hRefSelDijetPtEtaCMInLab->SetName("hRefSelDijetPtEtaCMInLab");

    // Center of mass frame
    TH2D *hRecoDijetPtEtaCM = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMWeighted"));
    // TH2D *hRecoDijetPtEtaCM = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMMatched"));
    if (!hRecoDijetPtEtaCM) {
        std::cerr << "Error: hRecoDijetPtEtaCMWeighted not found in file." << std::endl;
        return;
    }
    hRecoDijetPtEtaCM->SetName("hRecoDijetPtEtaCM");
    TH2D *hGenDijetPtEtaCM = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaCMWeighted"));
    if (!hGenDijetPtEtaCM) {
        std::cerr << "Error: hGenDijetPtEtaCMWeighted not found in file." << std::endl;
        return;
    }
    hGenDijetPtEtaCM->SetName("hGenDijetPtEtaCM");
    TH2D *hRefDijetPtEtaCM = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaCMWeighted"));
    if (!hRefDijetPtEtaCM) {
        std::cerr << "Error: hRefDijetPtEtaCMWeighted not found in file." << std::endl;
        return;
    }
    hRefDijetPtEtaCM->SetName("hRefDijetPtEtaCM");
    TH2D *hRefSelDijetPtEtaCM = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaCMWeighted"));
    if (!hRefSelDijetPtEtaCM) {
        std::cerr << "Error: hRefSelDijetPtEtaCMWeighted not found in file." << std::endl;
        return;
    }
    hRefSelDijetPtEtaCM->SetName("hRefSelDijetPtEtaCM");
    TH2D *hRecoDijetPtEtaLabInCM = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaLabInCM"));
    if (!hRecoDijetPtEtaLabInCM) {
        std::cerr << "Error: hRecoDijetPtEtaLabInCM not found in file." << std::endl;
        return;
    }
    hRecoDijetPtEtaLabInCM->SetName("hRecoDijetPtEtaLabInCM");
    TH2D *hGenDijetPtEtaLabInCM = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaLabInCM"));
    if (!hGenDijetPtEtaLabInCM) {
        std::cerr << "Error: hGenDijetPtEtaLabInCM not found in file." << std::endl;
        return;
    }
    hGenDijetPtEtaLabInCM->SetName("hGenDijetPtEtaLabInCM");
    TH2D *hRefDijetPtEtaLabInCM = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaLabInCM"));
    if (!hRefDijetPtEtaLabInCM) {
        std::cerr << "Error: hRefDijetPtEtaLabInCM not found in file." << std::endl;
        return;
    }
    hRefDijetPtEtaLabInCM->SetName("hRefDijetPtEtaLabInCM");
    TH2D *hRefSelDijetPtEtaLabInCM = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaLabInCM"));
    if (!hRefSelDijetPtEtaLabInCM) {
        std::cerr << "Error: hRefSelDijetPtEtaLabInCM not found in file." << std::endl;
        return;
    }
    hRefSelDijetPtEtaLabInCM->SetName("hRefSelDijetPtEtaLabInCM");

    TH2D *hRecoDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMForwardWeighted"));
    if (!hRecoDijetPtEtaCMForward) {
        std::cerr << "Error: hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    hRecoDijetPtEtaCMForward->SetName("hRecoDijetPtEtaCMForward");
    TH2D *hGenDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaCMForwardWeighted"));
    if (!hGenDijetPtEtaCMForward) {
        std::cerr << "Error: hGenDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    hGenDijetPtEtaCMForward->SetName("hGenDijetPtEtaCMForward");
    TH2D *hRefDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaCMForwardWeighted"));
    if (!hRefDijetPtEtaCMForward) {
        std::cerr << "Error: hRefDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    hRefDijetPtEtaCMForward->SetName("hRefDijetPtEtaCMForward");
    // TH2D *hRefSelDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaCMForwardWeighted"));
    // if (!hRefSelDijetPtEtaCMForward) {
    //     std::cerr << "Error: hRefSelDijetPtEtaCMForwardWeighted not found in file." << std::endl;
    //     return;
    // }
    
    TH2D *hRecoDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMBackwardWeighted"));
    if (!hRecoDijetPtEtaCMBackward) {
        std::cerr << "Error: hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    hRecoDijetPtEtaCMBackward->SetName("hRecoDijetPtEtaCMBackward");
    TH2D *hGenDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaCMBackwardWeighted"));
    if (!hGenDijetPtEtaCMBackward) {
        std::cerr << "Error: hGenDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    hGenDijetPtEtaCMBackward->SetName("hGenDijetPtEtaCMBackward");
    TH2D *hRefDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaCMBackwardWeighted"));
    if (!hRefDijetPtEtaCMBackward) {
        std::cerr << "Error: hRefDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    hRefDijetPtEtaCMBackward->SetName("hRefDijetPtEtaCMBackward");
    // TH2D *hRefSelDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaCMBackwardWeighted"));
    // if (!hRefSelDijetPtEtaCMBackward) {
    //     std::cerr << "Error: hRefSelDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
    //     return;
    // }

    // Distributions for laboratory frame
    TH1D *hRecoEtaLab{nullptr};
    TH1D *hGenEtaLab{nullptr};
    TH1D *hRefEtaLab{nullptr};
    TH1D *hRefSelEtaLab{nullptr};
    TH1D *hReco2GenEtaLab{nullptr};
    TH1D *hRef2GenEtaLab{nullptr};
    TH1D *hRefSel2GenEtaLab{nullptr};

    TH1D *hRecoEtaCMInLab{nullptr};
    TH1D *hGenEtaCMInLab{nullptr};
    TH1D *hRefEtaCMInLab{nullptr};
    TH1D *hRefSelEtaCMInLab{nullptr};
    TH1D *hReco2GenEtaCMInLab{nullptr};
    TH1D *hRef2GenEtaCMInLab{nullptr};
    TH1D *hRefSel2GenEtaCMInLab{nullptr};

    TH1D *hRecoPhiLab{nullptr};
    TH1D *hGenPhiLab{nullptr};
    TH1D *hRefPhiLab{nullptr};
    TH1D *hRefSelPhiLab{nullptr};
    TH1D *hReco2GenPhiLab{nullptr};
    TH1D *hRef2GenPhiLab{nullptr};
    TH1D *hRefSel2GenPhiLab{nullptr};

    // Distributions for center of mass frame
    TH1D *hRecoEtaCM{nullptr};
    TH1D *hGenEtaCM{nullptr};
    TH1D *hRefEtaCM{nullptr};
    TH1D *hRefSelEtaCM{nullptr};
    TH1D *hReco2GenEtaCM{nullptr};
    TH1D *hRef2GenEtaCM{nullptr};
    TH1D *hRefSel2GenEtaCM{nullptr};

    TH1D *hRecoEtaLabInCM{nullptr};
    TH1D *hGenEtaLabInCM{nullptr};
    TH1D *hRefEtaLabInCM{nullptr};
    TH1D *hRefSelEtaLabInCM{nullptr};
    TH1D *hReco2GenEtaLabInCM{nullptr};
    TH1D *hRef2GenEtaLabInCM{nullptr};
    TH1D *hRefSel2GenEtaLabInCM{nullptr};

    TH1D *hRecoEtaCMForward{nullptr};
    TH1D *hRecoEtaCMBackward{nullptr};
    TH1D *hGenEtaCMForward{nullptr};
    TH1D *hGenEtaCMBackward{nullptr};
    TH1D *hRefEtaCMForward{nullptr};
    TH1D *hRefEtaCMBackward{nullptr};
    // TH1D *hRefSelEtaCMForward{nullptr};
    // TH1D *hRefSelEtaCMBackward{nullptr};

    TH1D *hRecoEtaFBCM{nullptr};
    TH1D *hGenEtaFBCM{nullptr};
    TH1D *hRefEtaFBCM{nullptr};
    // TH1D *hRefSelEtaFBCM{nullptr};

    TH1D *hRecoPhiCM{nullptr};
    TH1D *hGenPhiCM{nullptr};
    TH1D *hRefPhiCM{nullptr};
    TH1D *hRefSelPhiCM{nullptr};
    TH1D *hReco2GenPhiCM{nullptr};
    TH1D *hRef2GenPhiCM{nullptr};
    TH1D *hRefSel2GenPhiCM{nullptr};

    // Forward, backward, and ratio distributions recalculated from full distribution in CM
    TH1D *hRecoEtaCMForwardRecalc {nullptr};
    TH1D *hRecoEtaCMBackwardRecalc{nullptr};
    TH1D *hGenEtaCMForwardRecalc{nullptr};
    TH1D *hGenEtaCMBackwardRecalc{nullptr};
    TH1D *hRefEtaCMForwardRecalc {nullptr};
    TH1D *hRefEtaCMBackwardRecalc{nullptr};

    TH1D *hRecoFBEtaCMRecalc{nullptr};
    TH1D *hGenFBEtaCMRecalc{nullptr};
    TH1D *hRefFBEtaCMRecalc{nullptr};

    TH1D *hRecoFBDoubleRatio{nullptr};
    TH1D *hGenFBDoubleRatio{nullptr};
    TH1D *hRefFBDoubleRatio{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    TLegend *leg{nullptr};
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Loop over dijet ptAve bins
    for (int i = 0; i < sizeOfPtVals - 1; ++i) {
        double ptLow = dijetPtVals[i];
        double ptHigh = dijetPtVals[i + 1];
        int ptBinLow = hRecoDijetPtEtaLab->GetXaxis()->FindBin(ptLow);
        int ptBinHigh = (hRecoDijetPtEtaLab->GetXaxis()->FindBin(ptHigh)) - 1;

        //
        // Laboratory frame
        //

        // Dijet pseudorapidity
        hRecoEtaLab = dynamic_cast<TH1D *>(hRecoDijetPtEtaLab->ProjectionY(Form("hRecoEtaLab_%d", i), ptBinLow, ptBinHigh));
        hGenEtaLab = dynamic_cast<TH1D *>(hGenDijetPtEtaLab->ProjectionY(Form("hGenEtaLab_%d", i), ptBinLow, ptBinHigh));
        hRefEtaLab = dynamic_cast<TH1D *>(hRefDijetPtEtaLab->ProjectionY(Form("hRefEtaLab_%d", i), ptBinLow, ptBinHigh));
        hRefSelEtaLab = dynamic_cast<TH1D *>(hRefSelDijetPtEtaLab->ProjectionY(Form("hRefSelEtaLab_%d", i), ptBinLow, ptBinHigh));
        hRecoEtaCMInLab = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMInLab->ProjectionY(Form("hRecoEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        hGenEtaCMInLab = dynamic_cast<TH1D *>(hGenDijetPtEtaCMInLab->ProjectionY(Form("hGenEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        hRefEtaCMInLab = dynamic_cast<TH1D *>(hRefDijetPtEtaCMInLab->ProjectionY(Form("hRefEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        hRefSelEtaCMInLab = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCMInLab->ProjectionY(Form("hRefSelEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        set1DStyle(hRecoEtaLab, 0);
        set1DStyle(hGenEtaLab, 5);
        set1DStyle(hRefEtaLab, 1);
        set1DStyle(hRefSelEtaLab, 2);
        set1DStyle(hRecoEtaCMInLab, 0);
        set1DStyle(hGenEtaCMInLab, 3);
        set1DStyle(hRefEtaCMInLab, 1);
        set1DStyle(hRefSelEtaCMInLab, 2);

        rescaleHisto1D(hRecoEtaLab);
        rescaleHisto1D(hGenEtaLab);
        rescaleHisto1D(hRefEtaLab);
        rescaleHisto1D(hRefSelEtaLab);
        rescaleHisto1D(hRecoEtaCMInLab);
        rescaleHisto1D(hGenEtaCMInLab);
        rescaleHisto1D(hRefEtaCMInLab);
        rescaleHisto1D(hRefSelEtaCMInLab);

        // Compute ratios for laboratory frame
        hReco2GenEtaLab = dynamic_cast<TH1D *>(hRecoEtaLab->Clone(Form("hReco2GenEtaLab_%d", i)));
        hReco2GenEtaLab->Divide(hReco2GenEtaLab, hGenEtaLab, 1., 1., "b");
        hRef2GenEtaLab = dynamic_cast<TH1D *>(hRefEtaLab->Clone(Form("hRef2GenEtaLab_%d", i)));
        hRef2GenEtaLab->Divide(hRef2GenEtaLab, hGenEtaLab, 1., 1., "b");
        hRefSel2GenEtaLab = dynamic_cast<TH1D *>(hRefSelEtaLab->Clone(Form("hRefSel2GenEtaLab_%d", i)));
        hRefSel2GenEtaLab->Divide(hRefSel2GenEtaLab, hGenEtaLab, 1., 1., "b");
        hReco2GenEtaCMInLab = dynamic_cast<TH1D *>(hRecoEtaCMInLab->Clone(Form("hReco2GenEtaCMInLab_%d", i)));
        hReco2GenEtaCMInLab->Divide(hReco2GenEtaCMInLab, hGenEtaCMInLab, 1., 1., "b");
        hRef2GenEtaCMInLab = dynamic_cast<TH1D *>(hRefEtaCMInLab->Clone(Form("hRef2GenEtaCMInLab_%d", i)));
        hRef2GenEtaCMInLab->Divide(hRef2GenEtaCMInLab, hGenEtaCMInLab, 1., 1., "b");
        hRefSel2GenEtaCMInLab = dynamic_cast<TH1D *>(hRefSelEtaCMInLab->Clone(Form("hRefSel2GenEtaCMInLab_%d", i)));
        hRefSel2GenEtaCMInLab->Divide(hRefSel2GenEtaCMInLab, hGenEtaCMInLab, 1., 1., "b");

        // std::cout << "pAve bin: " << i 
        //           << " ptLowBin: " << hRecoDijetPtEtaLab->GetXaxis()->FindBin(ptLow) 
        //           << " ptLow: " << hRecoDijetPtEtaLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaLab->GetXaxis()->FindBin(ptLow) )   
        //           << " ptHiBin: " << hRecoDijetPtEtaLab->GetXaxis()->FindBin(ptHi)-1
        //           << " ptHi: " << hRecoDijetPtEtaLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaLab->GetXaxis()->FindBin(ptHi)-1) << std::endl;

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoEtaLab, hRefEtaLab, hGenEtaLab, hRefSelEtaLab, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        drawDijetToGenComparison(c, hRecoEtaCMInLab, hRefEtaCMInLab, hGenEtaCMInLab, hRefSelEtaCMInLab, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                 false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCMInLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenEtaLab, hRef2GenEtaLab, hRefSel2GenEtaLab, nullptr,
                            (int)ptLow, (int)ptHigh,
                            false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        drawDijetToGenRatio(c, hReco2GenEtaCMInLab, hRef2GenEtaCMInLab, hRefSel2GenEtaCMInLab, nullptr,
                            (int)ptLow, (int)ptHigh,
                            false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCMInLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

                       /*
        // Dijet azimuthal angle
        hRecoPhiLab = dynamic_cast<TH1D *>(hRecoDijetPtEtaLab->ProjectionZ(Form("hRecoPhiLab_%d", i),
                                           ptBinLow, ptBinHigh));
        hGenPhiLab = dynamic_cast<TH1D *>(hGenDijetPtEtaLab->ProjectionZ(Form("hGenPhiLab_%d", i), 
                                          ptBinLow, ptBinHigh));
        hRefPhiLab = dynamic_cast<TH1D *>(hRefDijetPtEtaLab->ProjectionZ(Form("hRefPhiLab_%d", i), 
                                          ptBinLow, ptBinHigh));
        hRefSelPhiLab = dynamic_cast<TH1D *>(hRefSelDijetPtEtaLab->ProjectionZ(Form("hRefSelPhiLab_%d", i), 
                                             ptBinLow, ptBinHigh));
        set1DStyle(hRecoPhiLab, 0);
        set1DStyle(hGenPhiLab, 5);
        set1DStyle(hRefPhiLab, 1);
        set1DStyle(hRefSelPhiLab, 2);

        rescaleHisto1D(hRecoPhiLab);
        rescaleHisto1D(hGenPhiLab);
        rescaleHisto1D(hRefPhiLab);
        rescaleHisto1D(hRefSelPhiLab);

        // Compute ratios for laboratory frame
        hReco2GenPhiLab = dynamic_cast<TH1D *>(hRecoPhiLab->Clone(Form("hReco2GenPhiLab_%d", i)));
        hReco2GenPhiLab->Divide(hReco2GenPhiLab, hGenPhiLab, 1., 1., "b");
        hRef2GenPhiLab = dynamic_cast<TH1D *>(hRefPhiLab->Clone(Form("hRef2GenPhiLab_%d", i)));
        hRef2GenPhiLab->Divide(hRef2GenPhiLab, hGenPhiLab, 1., 1., "b");
        hRefSel2GenPhiLab = dynamic_cast<TH1D *>(hRefSelPhiLab->Clone(Form("hRefSel2GenPhiLab_%d", i)));
        hRefSel2GenPhiLab->Divide(hRefSel2GenPhiLab, hGenPhiLab, 1., 1., "b");

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoPhiLab, hRefPhiLab, hGenPhiLab, hRefSelPhiLab, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                 false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));
        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenPhiLab, hRef2GenPhiLab, hRefSel2GenPhiLab, nullptr,
                            (int)ptLow, (int)ptHigh,
                            false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

                       */

        
        //
        // Center of mass frame
        //

        // Dijet pseudorapidity
        hRecoEtaCM = dynamic_cast<TH1D *>(hRecoDijetPtEtaCM->ProjectionY(Form("hRecoEtaCM_%d", i), 
                                          ptBinLow, ptBinHigh));
        hGenEtaCM = dynamic_cast<TH1D *>(hGenDijetPtEtaCM->ProjectionY(Form("hGenEtaCM_%d", i), 
                                         ptBinLow, ptBinHigh));
        hRefEtaCM = dynamic_cast<TH1D *>(hRefDijetPtEtaCM->ProjectionY(Form("hRefEtaCM_%d", i), 
                                         ptBinLow, ptBinHigh));
        hRefSelEtaCM = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCM->ProjectionY(Form("hRefSelEtaCM_%d", i), 
                                            ptBinLow, ptBinHigh));
        hRecoEtaLabInCM = dynamic_cast<TH1D *>(hRecoDijetPtEtaLabInCM->ProjectionY(Form("hRecoEtaLabInCM_%d", i), 
                                               ptBinLow, ptBinHigh));
        hGenEtaLabInCM = dynamic_cast<TH1D *>(hGenDijetPtEtaLabInCM->ProjectionY(Form("hGenEtaLabInCM_%d", i), 
                                              ptBinLow, ptBinHigh));
        hRefEtaLabInCM = dynamic_cast<TH1D *>(hRefDijetPtEtaLabInCM->ProjectionY(Form("hRefEtaLabInCM_%d", i), 
                                              ptBinLow, ptBinHigh));
        hRefSelEtaLabInCM = dynamic_cast<TH1D *>(hRefSelDijetPtEtaLabInCM->ProjectionY(Form("hRefSelEtaLabInCM_%d", i), 
                                                 ptBinLow, ptBinHigh));
                                            
        set1DStyle(hRecoEtaCM, 0);
        set1DStyle(hGenEtaCM, 5);
        set1DStyle(hRefEtaCM, 1);
        set1DStyle(hRefSelEtaCM, 2);
        set1DStyle(hRecoEtaLabInCM, 0);
        set1DStyle(hGenEtaLabInCM, 3);
        set1DStyle(hRefEtaLabInCM, 1);
        set1DStyle(hRefSelEtaLabInCM, 2);

        rescaleHisto1D(hRecoEtaCM);
        rescaleHisto1D(hGenEtaCM);
        rescaleHisto1D(hRefEtaCM);
        rescaleHisto1D(hRefSelEtaCM);
        rescaleHisto1D(hRecoEtaLabInCM);
        rescaleHisto1D(hGenEtaLabInCM);
        rescaleHisto1D(hRefEtaLabInCM);
        rescaleHisto1D(hRefSelEtaLabInCM);

        // Compute ratios for center of mass frame
        hReco2GenEtaCM = dynamic_cast<TH1D *>(hRecoEtaCM->Clone(Form("hReco2GenEtaCM_%d", i)));
        hReco2GenEtaCM->Divide(hReco2GenEtaCM, hGenEtaCM, 1., 1., "b");
        hRef2GenEtaCM = dynamic_cast<TH1D *>(hRefEtaCM->Clone(Form("hRef2GenEtaCM_%d", i)));
        hRef2GenEtaCM->Divide(hRef2GenEtaCM, hGenEtaCM, 1., 1., "b");
        hRefSel2GenEtaCM = dynamic_cast<TH1D *>(hRefSelEtaCM->Clone(Form("hRefSel2GenEtaCM_%d", i)));
        hRefSel2GenEtaCM->Divide(hRefSel2GenEtaCM, hGenEtaCM, 1., 1., "b");
        hReco2GenEtaLabInCM = dynamic_cast<TH1D *>(hRecoEtaLabInCM->Clone(Form("hReco2GenEtaLabInCM_%d", i)));
        hReco2GenEtaLabInCM->Divide(hReco2GenEtaLabInCM, hGenEtaLabInCM, 1., 1., "b");
        hRef2GenEtaLabInCM = dynamic_cast<TH1D *>(hRefEtaLabInCM->Clone(Form("hRef2GenEtaLabInCM_%d", i)));
        hRef2GenEtaLabInCM->Divide(hRef2GenEtaLabInCM, hGenEtaLabInCM, 1., 1., "b");
        hRefSel2GenEtaLabInCM = dynamic_cast<TH1D *>(hRefSelEtaLabInCM->Clone(Form("hRefSel2GenEtaLabInCM_%d", i)));
        hRefSel2GenEtaLabInCM->Divide(hRefSel2GenEtaLabInCM, hGenEtaLabInCM, 1., 1., "b");

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoEtaCM, hRefEtaCM, hGenEtaCM, hRefSelEtaCM, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        drawDijetToGenComparison(c, hRecoEtaLabInCM, hRefEtaLabInCM, hGenEtaLabInCM, hRefSelEtaLabInCM, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLabInCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenEtaCM, hRef2GenEtaCM, hRefSel2GenEtaCM, nullptr,
                            (int)ptLow, (int)ptHigh,
                            true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_RecoRefGenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        drawDijetToGenRatio(c, hReco2GenEtaLabInCM, hRef2GenEtaLabInCM, hRefSel2GenEtaLabInCM, nullptr,
                            (int)ptLow, (int)ptHigh,
                            true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLabInCM_RecoRefGenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));


                       /*
        // Dijet azimuthal angle
        hRecoPhiCM = dynamic_cast<TH1D *>(hRecoDijetPtEtaCM->ProjectionZ(Form("hRecoPhiCM_%d", i), 
                                          ptBinLow, ptBinHigh));
        hGenPhiCM = dynamic_cast<TH1D *>(hGenDijetPtEtaCM->ProjectionZ(Form("hGenPhiCM_%d", i), 
                                         ptBinLow, ptBinHigh));
        hRefPhiCM = dynamic_cast<TH1D *>(hRefDijetPtEtaCM->ProjectionZ(Form("hRefPhiCM_%d", i), 
                                         ptBinLow, ptBinHigh));
        hRefSelPhiCM = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCM->ProjectionZ(Form("hRefSelPhiCM_%d", i), 
                                            ptBinLow, ptBinHigh));
        set1DStyle(hRecoPhiCM, 0);
        set1DStyle(hGenPhiCM, 5);
        set1DStyle(hRefPhiCM, 1);
        set1DStyle(hRefSelPhiCM, 2);

        rescaleHisto1D(hRecoPhiCM);
        rescaleHisto1D(hGenPhiCM);
        rescaleHisto1D(hRefPhiCM);
        rescaleHisto1D(hRefSelPhiCM);

        // Compute ratios for center of mass frame
        hReco2GenPhiCM = dynamic_cast<TH1D *>(hRecoPhiCM->Clone(Form("hReco2GenPhiCM_%d", i)));
        hReco2GenPhiCM->Divide(hReco2GenPhiCM, hGenPhiCM, 1., 1., "b");
        hRef2GenPhiCM = dynamic_cast<TH1D *>(hRefPhiCM->Clone(Form("hRef2GenPhiCM_%d", i)));
        hRef2GenPhiCM->Divide(hRef2GenPhiCM, hGenPhiCM, 1., 1., "b");
        hRefSel2GenPhiCM = dynamic_cast<TH1D *>(hRefSelPhiCM->Clone(Form("hRefSel2GenPhiCM_%d", i)));
        hRefSel2GenPhiCM->Divide(hRefSel2GenPhiCM, hGenPhiCM, 1., 1., "b");

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoPhiCM, hRefPhiCM, hGenPhiCM, hRefSelPhiCM, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));
        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenPhiCM, hRef2GenPhiCM, hRefSel2GenPhiCM, nullptr,
                            (int)ptLow, (int)ptHigh,
                            true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiCM_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        */


        //
        // Forward and backward dijet pseudorapidity
        //
        hRecoEtaCMForward = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMForward->ProjectionY(Form("hRecoEtaCMForward_%d", i), 
                                                 ptBinLow, ptBinHigh));
        hGenEtaCMForward = dynamic_cast<TH1D *>(hGenDijetPtEtaCMForward->ProjectionY(Form("hGenEtaCMForward_%d", i), 
                                                ptBinLow, ptBinHigh));
        hRefEtaCMForward = dynamic_cast<TH1D *>(hRefDijetPtEtaCMForward->ProjectionY(Form("hRefEtaCMForward_%d", i), 
                                                ptBinLow, ptBinHigh));
        // hRefSelEtaCMForward = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCMForward->ProjectionY(Form("hRefSelEtaCMForward_%d", i), 
        //                                            (int)ptBinLow, (int)ptBinHigh));

        hRecoEtaCMBackward = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMBackward->ProjectionY(Form("hRecoEtaCMBackward_%d", i),
                                                  ptBinLow, ptBinHigh));
        hGenEtaCMBackward = dynamic_cast<TH1D *>(hGenDijetPtEtaCMBackward->ProjectionY(Form("hGenEtaCMBackward_%d", i), 
                                                ptBinLow, ptBinHigh));
        hRefEtaCMBackward = dynamic_cast<TH1D *>(hRefDijetPtEtaCMBackward->ProjectionY(Form("hRefEtaCMBackward_%d", i), 
                                                ptBinLow, ptBinHigh));
        // hRefSelEtaCMBackward = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCMBackward->ProjectionY(Form("hRefSelEtaCMBackward_%d", i), 
        //                                             ptBinLow, ptBinHigh));

        set1DStyle(hRecoEtaCMForward, 0);
        set1DStyle(hGenEtaCMForward, 2);
        set1DStyle(hRefEtaCMForward, 1);
        // set1DStyle(hRefSelEtaCMForward, 2);

        set1DStyle(hRecoEtaCMBackward, 0);
        set1DStyle(hGenEtaCMBackward, 2);
        set1DStyle(hRefEtaCMBackward, 1);
        // set1DStyle(hRefSelEtaCMBackward, 2);

        // For forward and backward distributions number of bins is doubled (for default binning)
        int rebinFactor{1};
        hRecoEtaCMForward->Rebin(rebinFactor);
        hGenEtaCMForward->Rebin(rebinFactor);
        hRefEtaCMForward->Rebin(rebinFactor);
        // hRefSelEtaCMForward->Rebin(rebinFactor);

        hRecoEtaCMBackward->Rebin(rebinFactor);
        hGenEtaCMBackward->Rebin(rebinFactor);
        hRefEtaCMBackward->Rebin(rebinFactor);
        // hRefSelEtaCMBackward->Rebin(rebinFactor);

        // rescaleForwardBackward(hRecoEtaCMForward, hRecoEtaCMBackward);
        // rescaleForwardBackward(hGenEtaCMForward, hGenEtaCMBackward);
        // rescaleForwardBackward(hRefEtaCMForward, hRefEtaCMBackward);
        // rescaleForwardBackward(hRefSelEtaCMForward, hRefSelEtaCMBackward);

        // Compute ratios for forward and backward center of mass frame
        hRecoEtaFBCM = dynamic_cast<TH1D *>(hRecoEtaCMForward->Clone(Form("hRecoEtaFBCM_%d", i)));
        hRecoEtaFBCM->Divide(hRecoEtaFBCM, hRecoEtaCMBackward, 1., 1.);
        hGenEtaFBCM = dynamic_cast<TH1D *>(hGenEtaCMForward->Clone(Form("hGenEtaFBCM_%d", i)));
        hGenEtaFBCM->Divide(hGenEtaFBCM, hGenEtaCMBackward, 1., 1.);
        hRefEtaFBCM = dynamic_cast<TH1D *>(hRefEtaCMForward->Clone(Form("hRefEtaFBCM_%d", i)));
        hRefEtaFBCM->Divide(hRefEtaFBCM, hRefEtaCMBackward, 1., 1.);
        // hRefSelEtaFBCM = dynamic_cast<TH1D *>(hRefSelEtaCMForward->Clone(Form("hRefSelEtaFBCM_%d", i)));
        // hRefSelEtaFBCM->Divide(hRefSelEtaFBCM, hRefSelEtaCMBackward, 1., 1.);

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoEtaFBCM, hRefEtaFBCM, hGenEtaFBCM, nullptr, nullptr,
                                 (int)ptLow, (int)ptHigh,
                                 true, true, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaFBRatio_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), 
                       (int)ptLow, (int)ptHigh));

        //
        // Forward/backward manually computed from full 1D in CM
        //

        hRecoEtaCMForwardRecalc = dynamic_cast<TH1D*>( hRecoEtaCMForward->Clone( Form("hRecoEtaCMForwardRecalc_%d", i) ) );
        hRecoEtaCMForwardRecalc->GetYaxis()->SetTitle("Forward");
        hRecoEtaCMBackwardRecalc = dynamic_cast<TH1D*>( hRecoEtaCMBackward->Clone( Form("hRecoEtaCMBackwardRecalc_%d", i) ) );
        hRecoEtaCMBackwardRecalc->GetYaxis()->SetTitle("Backward"); 
        hGenEtaCMForwardRecalc = dynamic_cast<TH1D*>( hGenEtaCMForward->Clone( Form("hGenEtaCMForwardRecalc_%d", i) ) );
        hGenEtaCMForwardRecalc->GetYaxis()->SetTitle("Forward");
        hGenEtaCMBackwardRecalc = dynamic_cast<TH1D*>( hGenEtaCMBackward->Clone( Form("hGenEtaCMBackwardRecalc_%d", i) ) );
        hGenEtaCMBackwardRecalc->GetYaxis()->SetTitle("Backward"); 
        hRefEtaCMForwardRecalc = dynamic_cast<TH1D*>( hRefEtaCMForward->Clone( Form("hRefEtaCMForwardRecalc_%d", i) ) );
        hRefEtaCMForwardRecalc->GetYaxis()->SetTitle("Forward");
        hRefEtaCMBackwardRecalc = dynamic_cast<TH1D*>( hRefEtaCMBackward->Clone( Form("hRefEtaCMBackwardRecalc_%d", i) ) );
        hRefEtaCMBackwardRecalc->GetYaxis()->SetTitle("Backward"); 

        hRecoFBEtaCMRecalc = dynamic_cast<TH1D*>( hRecoEtaFBCM->Clone( Form("hRecoFBEtaCMRecalc_%d", i) ) );
        hRecoFBEtaCMRecalc->GetYaxis()->SetTitle("Forward/Backward"); 
        hGenFBEtaCMRecalc = dynamic_cast<TH1D*>( hGenEtaFBCM->Clone( Form("hGenFBEtaCMRecalc_%d", i) ) );
        hGenFBEtaCMRecalc->GetYaxis()->SetTitle("Forward/Backward"); 
        hRefFBEtaCMRecalc = dynamic_cast<TH1D*>( hRefEtaFBCM->Clone( Form("hRefFBEtaCMRecalc_%d", i) ) );
        hRefFBEtaCMRecalc->GetYaxis()->SetTitle("Forward/Backward"); 

        recalculateFBRatioFromFullDistribution(hRecoFBEtaCMRecalc, hRecoEtaCMForwardRecalc, hRecoEtaCMBackwardRecalc, hRecoEtaCM);
        recalculateFBRatioFromFullDistribution(hGenFBEtaCMRecalc, hGenEtaCMForwardRecalc, hGenEtaCMBackwardRecalc, hGenEtaCM);
        recalculateFBRatioFromFullDistribution(hRefFBEtaCMRecalc, hRefEtaCMForwardRecalc, hRefEtaCMBackwardRecalc, hRefEtaCM);

        drawDijetToGenRatio(c, hRecoFBEtaCMRecalc, hRefFBEtaCMRecalc, hGenFBEtaCMRecalc, nullptr, 
                                 (int)ptLow, (int)ptHigh,
                                 true, true, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaFBRatio_recalculated_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHigh));

        // Double ratios of FB (automated/manual)
        hRecoFBDoubleRatio = dynamic_cast<TH1D *>(hRecoEtaFBCM->Clone(Form("hRecoFBDoubleRatio_%d", i)));
        hRecoFBDoubleRatio->Divide( hRecoFBDoubleRatio, hRecoFBEtaCMRecalc, 1., 1., "b");
        hRecoFBDoubleRatio->GetYaxis()->SetTitle("Forward/Backward (Double Ratio)");
        hGenFBDoubleRatio = dynamic_cast<TH1D *>(hGenEtaFBCM->Clone(Form("hGenFBDoubleRatio_%d", i)));
        hGenFBDoubleRatio->Divide( hGenFBDoubleRatio, hGenFBEtaCMRecalc, 1., 1., "b");
        hGenFBDoubleRatio->GetYaxis()->SetTitle("Forward/Backward (Double Ratio)");
        hRefFBDoubleRatio = dynamic_cast<TH1D *>(hRefEtaFBCM->Clone(Form("hRefFBDoubleRatio_%d", i)));
        hRefFBDoubleRatio->Divide( hRefFBDoubleRatio, hRefFBEtaCMRecalc, 1., 1., "b");
        hRefFBDoubleRatio->GetYaxis()->SetTitle("Forward/Backward (Double Ratio)");

        drawDijetToGenRatio(c, hRecoFBDoubleRatio, hRefFBDoubleRatio, hGenFBDoubleRatio, nullptr,
                            (int)ptLow, (int)ptHigh,
                            true, true, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaFB_doubleRatio_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        (int)ptLow, (int)ptHigh));

    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)

    if (c) { delete c; c = nullptr; }
}

//_________________
// Function plots eta and pT distribution comparisons and ratios of reco, ref, refSel to gen jets
// in bins of ptHat.
// It retrieves histograms for reconstructed, generated, and reference jets, and computes the ratios.
// f is the input TFile containing the histograms.
// collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
// collisionEnergy: energy in TeV (default is 8.16 TeV for pPb)
// jetType: 0 = Inclusive, 1 = Lead, 2 = SubLead
// date: date string for saving the plots (default is "20250129")
void inclusiveJetJECClosures(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, int jetType = 0, int matchType = 0, TString date = "20250129") {
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.6;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TString tMatched = "All";
    tMatched = (matchType == 0) ? "All" : ((matchType == 1) ? "Matched" : "Unmatched");

    TString tJetType = "Inclusive"; 
    tJetType = (jetType == 0) ? "Inclusive" : ((jetType == 1) ? "Lead" : "SubLead");

    std::cout << "Jet type: " << tJetType.Data() << std::endl;
    std::cout << "Match type: " << tMatched.Data() << std::endl;

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } 
    else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } 
    else {
        directionStr = "comb";
    }
    

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // ptHat binning
    double ptHatVals[] { 15., 60., 90. };
    int sizeOfPtHatVals = sizeof(ptHatVals)/sizeof(ptHatVals[0]);

    // jet pT and eta binning
    double_t jetPtVals[] = {30., 50., 90., 120., 200., 500., 1500.};
    int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

    double jetEtaVals[] = {-3.6, -2.4, -1.6, 1.6, 2.4, 3.6};
    int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);

    std::cout << Form("Reco histogram to read: ") << Form("hReco%s%sJetPtEtaPtHat", tJetType.Data(), tMatched.Data() ) << std::endl;
    std::cout << Form("Gen histogram to read: ") << Form("hGen%sJetPtEtaPtHat", tJetType.Data() ) << std::endl;
    std::cout << Form("Ref histogram to read: ") << Form("hRef%sJetPtEtaPtHat", tJetType.Data() ) << std::endl;
    std::cout << Form("RefSel histogram to read: ") << Form("hRefSel%sJetPtEtaPtHat", tJetType.Data() ) << std::endl;

    // Retrieve histograms
    TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hReco%s%sJetPtEtaPtHat", tJetType.Data(), tMatched.Data() ) ) );
    if ( !hRecoPtEtaPtHat ) {
        std::cerr << "Histogram hRecoPtEtaPtHat not found in file." << std::endl; return;
    }
    // TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRecoMatchedJetPtEtaPtHat"));
    TH3D *hGenPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hGen%sJetPtEtaPtHat", tJetType.Data() ) ) );
    if ( !hGenPtEtaPtHat ) {
        std::cerr << "Histogram hGenPtEtaPtHat not found in file." << std::endl; return;
    }
    TH3D *hRefPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hRef%sJetPtEtaPtHat", tJetType.Data() ) ) );
    if ( !hRefPtEtaPtHat ) {
        std::cerr << "Histogram hRefPtEtaPtHat not found in file." << std::endl; return;
    }
    TH3D *hRefSelPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hRefSel%sJetPtEtaPtHat", tJetType.Data() ) ) );
    if ( !hRefSelPtEtaPtHat ) {
        std::cerr << "Histogram hRefSelPtEtaPtHat not found in file." << std::endl; return;
    }

    //
    // Declare canvases and histograms
    //

    TH2D *hRecoPtVsEta{nullptr};
    TH2D *hGenPtVsEta{nullptr};
    TH1D *hRecoPt{nullptr};
    TH1D *hGenPt{nullptr};

    TH1D *hRecoEta{nullptr};
    TH1D *hGenEta{nullptr};
    TH1D *hRefEta{nullptr};
    TH1D *hRefSelEta{nullptr};

    TH1D *hReco2GenEta{nullptr};
    TH1D *hRef2GenEta{nullptr};
    TH1D *hRefSel2GenEta{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid(1, 1);

    //
    // Perform analysis and build comparisons for different ptHat intervals.
    // First part loops over jet pT, second part over jet eta.
    //

    // Loop over ptHat bins
    for (int i = 0; i < sizeOfPtHatVals; i++) {

        double ptHatLow = ptHatVals[i];
        double ptHatHigh = hRecoPtEtaPtHat->GetZaxis()->GetBinUpEdge( hRecoPtEtaPtHat->GetNbinsZ() );
        int ptHatBinLow = hRecoPtEtaPtHat->GetZaxis()->FindBin(ptHatLow);
        int ptHatBinHigh = hRecoPtEtaPtHat->GetNbinsZ();

        //
        // Loop over jet pT bins
        //
        for (int j = 0; j < sizeOfJetPtVals-1; j++) {

            //
            // Make projections of the 3D histograms
            //

            double jetPtLow = jetPtVals[j];
            double jetPtHigh = jetPtVals[j+1];
            int jetPtBinLow = hRecoPtEtaPtHat->GetYaxis()->FindBin(jetPtLow);
            int jetPtBinHigh = hRecoPtEtaPtHat->GetYaxis()->FindBin(jetPtHigh) - 1;


            // Reco jets
            hRecoEta = dynamic_cast<TH1D *>(hRecoPtEtaPtHat->ProjectionX(Form("hRecoEta_%d_%d", i, j),
                                                                               jetPtBinLow, jetPtBinHigh,
                                                                               ptHatBinLow, ptHatBinHigh));
            if ( !hRecoEta ) {
                std::cout << "hRecoEta_" << i << "_" << j << " does not exist" << std::endl;
                return;
            }
            rescaleHisto1D(hRecoEta);
            set1DStyle(hRecoEta, 0);

            // Gen jets
            hGenEta = dynamic_cast<TH1D *>(hGenPtEtaPtHat->ProjectionX(Form("hGenEta_%d_%d", i, j),
                                           jetPtBinLow, jetPtBinHigh,
                                           ptHatBinLow, ptHatBinHigh));
            if ( !hGenEta ) {
                std::cout << "hGenEta_" << i << "_" << j << " does not exist" << std::endl;
                return;
            }
            rescaleHisto1D(hGenEta);
            set1DStyle(hGenEta, 5);

            // Ref jets
            hRefEta = dynamic_cast<TH1D *>(hRefPtEtaPtHat->ProjectionX(Form("hRefEta_%d_%d", i, j),
                                           jetPtBinLow, jetPtBinHigh,
                                           ptHatBinLow, ptHatBinHigh));
            if ( !hRefEta ) {
                std::cout << "hRefEta_" << i << "_" << j << " does not exist" << std::endl;
                return;
            }
            rescaleHisto1D(hRefEta);
            set1DStyle(hRefEta, 1);
        
            // RefSel jets
            hRefSelEta = dynamic_cast<TH1D *>(hRefSelPtEtaPtHat->ProjectionX(Form("hRefSelEta_%d_%d", i, j),
                                              jetPtBinLow, jetPtBinHigh,
                                              ptHatBinLow, ptHatBinHigh));
            if ( !hRefSelEta ) {
                std::cout << "hRefSelEta_" << i << "_" << j << " does not exist" << std::endl;
                return;
            }
            set1DStyle(hRefSelEta, 2);
            rescaleHisto1D(hRefSelEta);

            int rebinFactor{2};
            hRecoEta->Rebin(rebinFactor);
            hGenEta->Rebin(rebinFactor);
            hRefEta->Rebin(rebinFactor);
            hRefSelEta->Rebin(rebinFactor);

            //
            // Ratios of reco, ref and refSel to gen
            //
            hReco2GenEta = dynamic_cast<TH1D *>(hRecoEta->Clone(Form("hReco2GenEta_%d_%d", i, j)));
            hReco2GenEta->Divide(hReco2GenEta, hGenEta, 1., 1., "b");
            hRef2GenEta = dynamic_cast<TH1D *>(hRefEta->Clone(Form("hRef2GenEta_%d_%d", i, j)));
            hRef2GenEta->Divide(hRef2GenEta, hGenEta, 1., 1., "b");
            hRefSel2GenEta = dynamic_cast<TH1D *>(hRefSelEta->Clone(Form("hRefSel2GenEta_%d_%d", i, j)));
            hRefSel2GenEta->Divide(hRefSel2GenEta, hGenEta, 1., 1., "b");
            
            //
            // Plot comparisons
            //
            drawSingleJetToGenComparison(c, hRecoEta, hRefEta, hGenEta, hRefSelEta, nullptr,
                                         jetPtLow, jetPtHigh, ptHatLow, 
                                         false, false, collisionSystem, collisionEnergy);
            c->SaveAs(Form("%s/%s_%s_jetEta_RecoRefGenComp_ptAve_%d_%d_ptHatLow_%d.pdf", 
                      date.Data(), collSystemStr.Data(), directionStr.Data(),
                      (int)jetPtLow, (int)jetPtHigh, (int)ptHatLow));

            //
            // Plot ratios
            //
            drawSingleJetToGenRatio(c, hReco2GenEta, hRef2GenEta, hRefSel2GenEta, nullptr,
                                    jetPtLow, jetPtHigh, ptHatLow, false, false, collisionSystem, collisionEnergy);
            c->SaveAs(Form("%s/%s_%s_jetEta_RecoRef2GenRatio_ptAve_%d_%d_ptHatLow_%d.pdf", 
                      date.Data(), collSystemStr.Data(), directionStr.Data(),
                      (int)jetPtLow, (int)jetPtHigh, (int)ptHatLow));


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

    if (c) { delete c; c = nullptr; }
}

//________________
//
// Plot comparison of inclusive jet distributions for different pT and eta bins
// Compare both reco data to reco MC and reco data to gen MC
//
void data2mcInclusiveJetComparison(TFile *fData, TFile *fMc, int collSystem = 1, double energy = 8.16, TString date = "20250224") {
    
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

    // Direction 
    TString directionStr;
    TString filename = fMc->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", (int)(1000 * energy));
    collSystemStr += Form("_%s", prefix.Data());

    TH2D *hRecoDataPtEta = dynamic_cast<TH2D *>(fData->Get("hRecoInclusiveAllJetPtEtaCM"));
    if (!hRecoDataPtEta) {
        std::cerr << "Error: hRecoDataPtEta not found in data file!" << std::endl;
        return;
    }
    hRecoDataPtEta->SetName("hRecoDataPtEta");
    TH2D *hRecoMcPtEta = dynamic_cast<TH2D *>(fMc->Get("hRecoInclusiveAllJetPtEtaCM"));
    if (!hRecoMcPtEta) {
        std::cerr << "Error: hRecoMcPtEta not found in MC file!" << std::endl;
        return;
    }
    hRecoMcPtEta->SetName("hRecoMcPtEta");
    // TH3D *hRecoMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hRecoMatchedJetPtEtaPtHat"));
    // hRecoMcPtEtaPtHat->SetName("hRecoMcPtEtaPtHat");
    TH2D *hGenMcPtEta = dynamic_cast<TH2D *>(fMc->Get("hGenInclusiveJetPtEtaCM"));
    if (!hGenMcPtEta) {
        std::cerr << "Error: hGenMcPtEta not found in MC file!" << std::endl;
        return;
    }
    hGenMcPtEta->SetName("hGenMcPtEta");

    // jet pT and eta binning
    double jetPtVals[] = {40., 80., 120., 180., 250., 500.};
    int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

    double jetEtaVals[] = {-3.6, -2.4, -1.6, 1.6, 2.4, 3.6};
    int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);


    TH1D *hRecoDataEta{nullptr};
    TH1D *hRecoMcEta{nullptr};
    TH1D *hGenMcEta{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);

    // Jet pT binning
    for ( int i = 0; i < sizeOfJetPtVals-1; i++) {

        double jetPtLow = jetPtVals[i];
        double jetPtHigh = jetPtVals[i+1];
        int jetPtBinLow = hRecoDataPtEta->GetYaxis()->FindBin(jetPtLow);
        int jetPtBinHigh = hRecoDataPtEta->GetYaxis()->FindBin(jetPtHigh)-1;

        //
        // Reco data 2 Reco MC comparison
        //
        // Make projection of the 3D histograms
        hRecoDataEta = dynamic_cast<TH1D *>(hRecoDataPtEta->ProjectionX(Form("hRecoDataEta_%d", i),
                                            jetPtBinLow, jetPtBinHigh));
        if (!hRecoDataEta) {
            std::cerr << Form("Error: hRecoDataEta_%d not found!", i) << std::endl;
            return;
        }
        rescaleHisto1D(hRecoDataEta);
        set1DStyle(hRecoDataEta, 0);

        hRecoMcEta = dynamic_cast<TH1D *>(hRecoMcPtEta->ProjectionX(Form("hRecoMcEta_%d", i),
                                          jetPtBinLow, jetPtBinHigh));
        if (!hRecoMcEta) {
            std::cerr << Form("Error: hRecoMcEta_%d not found!", i) << std::endl;
            return;
        }
        rescaleHisto1D(hRecoMcEta);
        set1DStyle(hRecoMcEta, 1);

        // Plot distributions
        drawJecComparison(c, hRecoDataEta, hRecoMcEta,
            jetPtLow, jetPtHigh, 15,
            "Data", "Reco MC", collSystem, energy, false, true,
            hRecoDataPtEta->GetXaxis()->GetBinLowEdge(1),
            hRecoDataPtEta->GetXaxis()->GetBinUpEdge(hRecoDataPtEta->GetXaxis()->GetNbins()));

        c->SaveAs( Form("%s/%s_closureEta_Reco2Reco_pt_%d_%d.pdf", 
            date.Data(), collSystemStr.Data(), (int)jetPtLow, (int)jetPtHigh ) );


        //
        // Reco data 2 Gen MC comparison
        //

        c->cd(); c->Clear(); setPadStyle();
        // Make projection of the 3D histograms
        hGenMcEta = dynamic_cast<TH1D *>(hGenMcPtEta->ProjectionX(Form("hGenMcEta_%d", i),
            jetPtBinLow, jetPtBinHigh));
        if (!hGenMcEta) {
            std::cerr << Form("Error: hGenMcEta_%d not found!", i) << std::endl;
            return;
        }
        rescaleHisto1D(hGenMcEta);
        set1DStyle(hGenMcEta, 5);

        // Plot distributions
        drawJecComparison(c, hRecoDataEta, hGenMcEta,
                          jetPtLow, jetPtHigh, 15,
                          "Data", "Gen MC", collSystem, energy, false, true,
                          hGenMcPtEta->GetXaxis()->GetBinLowEdge(1),
                          hGenMcPtEta->GetXaxis()->GetBinUpEdge(hGenMcPtEta->GetXaxis()->GetNbins())
                        );

        // Save canvas
        c->SaveAs( Form("%s/%s_closureEta_Reco2Gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), (int)jetPtLow, (int)jetPtHigh) );
    } // for (unsigned int i = 0; i < jetPtBinsLow.size(); i++)

    if (c) { delete c; c = nullptr; }
}

//________________
//
// Plot comparison of dijet eta distributions for different pTave bins
// Compare both reco data to reco MC and reco data to gen MC
//
void data2mcDijetComparison(TFile *fData, TFile *fMc, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250224") {

    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

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

    // Direction
    TString directionStr;
    TString filename = fMc->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // double dijetPtVals[] {  50.,  60.,   70.,  80.,  90.,
    //                         100., 110.,  120., 130., 140.,
    //                         150., 160.,  180., 200., 250., 
    //                         300., 500.};
    double dijetPtVals[] {  50., 80., 100., 120., 180., 200., 250., 300., 350., 500.};
    int sizeOfDijetPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    //
    // Declare histograms
    //

    // Data
    TH3D *hRecoDataDijetPtEtaLab = dynamic_cast<TH3D *>(fData->Get("hRecoDijetPtEtaWeighted"));
    if ( !hRecoDataDijetPtEtaLab ) {
        std::cerr << "Data histogram hRecoDijetPtEtaWeighted not found in file." << std::endl; return;
    }
    hRecoDataDijetPtEtaLab->SetName("hRecoDataDijetPtEtaLab");
    TH3D *hRecoDataDijetPtEtaCM = dynamic_cast<TH3D *>(fData->Get("hRecoDijetPtEtaCMWeighted"));
    if ( !hRecoDataDijetPtEtaCM ) {
        std::cerr << "Data histogram hRecoDijetPtEtaCMWeighted not found in file." << std::endl; return;
    }
    hRecoDataDijetPtEtaCM->SetName("hRecoDataDijetPtEtaCM");
    TH2D *hRecoDataDijetPtEtaCMForward = dynamic_cast<TH2D *>(fData->Get("hRecoDijetPtEtaCMForwardWeighted"));
    if ( !hRecoDataDijetPtEtaCMForward ) {
        std::cerr << "Data histogram hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
    }
    hRecoDataDijetPtEtaCMForward->SetName("hRecoDataDijetPtEtaCMForward");
    TH2D *hRecoDataDijetPtEtaCMBackward = dynamic_cast<TH2D *>(fData->Get("hRecoDijetPtEtaCMBackwardWeighted"));
    if ( !hRecoDataDijetPtEtaCMBackward ) {
        std::cerr << "Data histogram hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
    }
    hRecoDataDijetPtEtaCMBackward->SetName("hRecoDataDijetPtEtaCMBackward");

    TH1D *hRecoDataDijetEtaLab{nullptr};
    TH1D *hRecoDataDijetEtaCM{nullptr};
    TH1D *hRecoDataDijetEtaForward{nullptr};
    TH1D *hRecoDataDijetEtaBackward{nullptr};
    TH1D *hRecoDataDijetFBRatio{nullptr};

    // MC

    // Reco
    TH3D *hRecoMcDijetPtEtaLab = dynamic_cast<TH3D *>(fMc->Get("hRecoDijetPtEtaWeighted"));
    if ( !hRecoMcDijetPtEtaLab ) {
        std::cerr << "MC histogram hRecoDijetPtEtaWeighted not found in file." << std::endl; return;
    }
    hRecoMcDijetPtEtaLab->SetName("hRecoMcDijetPtEtaLab");
    TH3D *hRecoMcDijetPtEtaCM = dynamic_cast<TH3D *>(fMc->Get("hRecoDijetPtEtaCMWeighted"));
    if ( !hRecoMcDijetPtEtaCM ) {
        std::cerr << "MC histogram hRecoDijetPtEtaCMWeighted not found in file." << std::endl; return;
    }
    hRecoMcDijetPtEtaCM->SetName("hRecoMcDijetPtEtaCM");
    TH2D *hRecoMcDijetPtEtaCMForward = dynamic_cast<TH2D *>(fMc->Get("hRecoDijetPtEtaCMForwardWeighted"));
    if ( !hRecoMcDijetPtEtaCMForward ) {
        std::cerr << "MC histogram hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
    }
    hRecoMcDijetPtEtaCMForward->SetName("hRecoMcDijetPtEtaCMForward");
    TH2D *hRecoMcDijetPtEtaCMBackward = dynamic_cast<TH2D *>(fMc->Get("hRecoDijetPtEtaCMBackwardWeighted"));
    if ( !hRecoMcDijetPtEtaCMBackward ) {
        std::cerr << "MC histogram hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
    }
    hRecoMcDijetPtEtaCMBackward->SetName("hRecoMcDijetPtEtaCMBackward");

    TH1D *hRecoMcDijetEtaLab{nullptr};
    TH1D *hRecoMcDijetEtaCM{nullptr};
    TH1D *hRecoMcDijetEtaForward{nullptr};
    TH1D *hRecoMcDijetEtaBackward{nullptr};
    TH1D *hRecoMcDijetFBRatio{nullptr};

    // Gen
    TH3D *hGenMcDijetPtEtaLab = dynamic_cast<TH3D *>(fMc->Get("hGenDijetPtEtaWeighted"));
    if ( !hGenMcDijetPtEtaLab ) {
        std::cerr << "MC histogram hGenDijetPtEtaWeighted not found in file." << std::endl; return;
    }
    hGenMcDijetPtEtaLab->SetName("hGenMcDijetPtEtaLab");
    TH3D *hGenMcDijetPtEtaCM = dynamic_cast<TH3D *>(fMc->Get("hGenDijetPtEtaCMWeighted"));
    if ( !hGenMcDijetPtEtaCM ) {
        std::cerr << "MC histogram hGenDijetPtEtaCMWeighted not found in file." << std::endl; return;
    }
    hGenMcDijetPtEtaCM->SetName("hGenMcDijetPtEtaCM");
    TH2D *hGenMcDijetPtEtaCMForward = dynamic_cast<TH2D *>(fMc->Get("hGenDijetPtEtaCMForwardWeighted"));
    if ( !hGenMcDijetPtEtaCMForward ) {
        std::cerr << "MC histogram hGenDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
    }
    hGenMcDijetPtEtaCMForward->SetName("hGenMcDijetPtEtaCMForward");
    TH2D *hGenMcDijetPtEtaCMBackward =  dynamic_cast<TH2D *>(fMc->Get("hGenDijetPtEtaCMBackwardWeighted"));
    if ( !hGenMcDijetPtEtaCMBackward ) {
        std::cerr << "MC histogram hGenDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
    }
    hGenMcDijetPtEtaCMBackward->SetName("hGenMcDijetPtEtaCMBackward");

    TH1D *hGenMcDijetEtaLab{nullptr};
    TH1D *hGenMcDijetEtaCM{nullptr};
    TH1D *hGenMcDijetEtaForward{nullptr};
    TH1D *hGenMcDijetEtaBackward{nullptr};
    TH1D *hGenMcDijetFBRatio{nullptr};

    // Ref
    TH3D *hRefMcDijetPtEtaLab = dynamic_cast<TH3D *>(fMc->Get("hRefDijetPtEtaWeighted"));
    if ( !hRefMcDijetPtEtaLab ) {
        std::cerr << "MC histogram hRefDijetPtEtaWeighted not found in file." << std::endl; return;
    }
    hRefMcDijetPtEtaLab->SetName("hRefMcDijetPtEtaLab");
    TH3D *hRefMcDijetPtEtaCM = dynamic_cast<TH3D *>(fMc->Get("hRefDijetPtEtaCMWeighted"));
    if ( !hRefMcDijetPtEtaCM ) {
        std::cerr << "MC histogram hRefDijetPtEtaCMWeighted not found in file." << std::endl; return;
    }
    hRefMcDijetPtEtaCM->SetName("hRefMcDijetPtEtaCM");
    TH2D *hRefMcDijetPtEtaCMForward = dynamic_cast<TH2D *>(fMc->Get("hRefDijetPtEtaCMForwardWeighted"));
    if ( !hRefMcDijetPtEtaCMForward ) {
        std::cerr << "MC histogram hRefDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
    }
    hRefMcDijetPtEtaCMForward->SetName("hRefMcDijetPtEtaCMForward");
    TH2D *hRefMcDijetPtEtaCMBackward = dynamic_cast<TH2D *>(fMc->Get("hRefDijetPtEtaCMBackwardWeighted"));
    if ( !hRefMcDijetPtEtaCMBackward ) {
        std::cerr << "MC histogram hRefDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
    }
    hRefMcDijetPtEtaCMBackward->SetName("hRefMcDijetPtEtaCMBackward");

    TH1D *hRefMcDijetEtaLab{nullptr};
    TH1D *hRefMcDijetEtaCM{nullptr};
    TH1D *hRefMcDijetEtaForward{nullptr};
    TH1D *hRefMcDijetEtaBackward{nullptr};
    TH1D *hRefMcDijetFBRatio{nullptr};

    // Ratio of full distributions to gen
    TH1D *hRecoData2GenDijetEtaLab{nullptr};
    TH1D *hRecoMc2GenDijetEtaLab{nullptr};
    TH1D *hRef2GenDijetEtaLab{nullptr};

    TH1D *hRecoData2GenDijetEtaCM{nullptr};
    TH1D *hRecoMc2GenDijetEtaCM{nullptr};
    TH1D *hRef2GenDijetEtaCM{nullptr};

    //
    // Canvas 
    //
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TCanvas *c = new TCanvas("c", "c", 1200, 1200);
    setPadStyle();

    //
    // Loop over the dijet pTave bins, make projections and plot distributions
    //
    for (int i = 0; i < sizeOfDijetPtVals - 1; ++i) {

        double ptLow = dijetPtVals[i];
        double ptHi = dijetPtVals[i + 1];
        int ptBinLow = hRecoDataDijetPtEtaLab->GetXaxis()->FindBin(ptLow);
        int ptBinHigh = hRecoDataDijetPtEtaLab->GetXaxis()->FindBin(ptHi) - 1;

        //
        // Data
        //

        hRecoDataDijetEtaLab = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaLab->ProjectionY( Form("hRecoDataDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoDataDijetEtaLab ) {
            std::cerr << Form("Data histogram not found: hRecoDataDijetEtaLab_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaLab->SetName( Form("hRecoDataDijetEtaLab_%d", i) );
        set1DStyle( hRecoDataDijetEtaLab, 2 );
        rescaleHisto1D( hRecoDataDijetEtaLab );

        hRecoDataDijetEtaCM = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaCM->ProjectionY( Form("hRecoDataDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoDataDijetEtaCM ) {
            std::cerr << Form("Data histogram not found: hRecoDataDijetEtaCM_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaCM->SetName( Form("hRecoDataDijetEtaCM_%d", i) );
        set1DStyle( hRecoDataDijetEtaCM, 2 );
        rescaleHisto1D( hRecoDataDijetEtaCM );

        hRecoDataDijetEtaForward = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaCMForward->ProjectionY( Form("hRecoDataDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoDataDijetEtaForward ) {
            std::cerr << Form("Data histogram not found: hRecoDataDijetEtaForward_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaForward->SetName( Form("hRecoDataDijetEtaForward_%d", i) );
        set1DStyle( hRecoDataDijetEtaForward, 2 );
        
        hRecoDataDijetEtaBackward = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaCMBackward->ProjectionY( Form("hRecoDataDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoDataDijetEtaBackward ) {
            std::cerr << Form("Data histogram not found: hRecoDataDijetEtaBackward_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaBackward->SetName( Form("hRecoDataDijetEtaBackward_%d", i) );
        set1DStyle( hRecoDataDijetEtaBackward, 2 );

        // Normalize forward and backward distributions
        // double integralForward = hRecoDataDijetEtaForward->Integral();
        // double integralBackward = hRecoDataDijetEtaBackward->Integral();
        // hRecoDataDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
        // hRecoDataDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hRecoDataDijetFBRatio = dynamic_cast<TH1D*>( hRecoDataDijetEtaForward->Clone( Form("hRecoDataDijetFBRatio_%d", i) ) );
        hRecoDataDijetFBRatio->SetName( Form("hRecoDataDijetFBRatio_%d", i) );
        hRecoDataDijetFBRatio->Divide( hRecoDataDijetFBRatio, hRecoDataDijetEtaBackward, 1., 1. );

        c->cd(); c->Clear(); setPadStyle();
        gPad->SetGrid(1, 1);
        hRecoDataDijetEtaLab->Draw();
        hRecoDataDijetEtaLab->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hRecoDataDijetEtaLab->GetYaxis()->SetRangeUser(0.00001, 0.15);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        c->cd(); c->Clear(); setPadStyle();
        gPad->SetGrid(1, 1);
        hRecoDataDijetEtaCM->Draw();
        hRecoDataDijetEtaCM->GetXaxis()->SetRangeUser(-2.5, 2.5);
        hRecoDataDijetEtaCM->GetYaxis()->SetRangeUser(0.00001, 0.15);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

        c->cd(); c->Clear(); setPadStyle();
        gPad->SetGrid(1, 1);
        hRecoDataDijetFBRatio->Draw();
        hRecoDataDijetFBRatio->GetXaxis()->SetRangeUser(-2.5, 2.5);
        hRecoDataDijetFBRatio->GetYaxis()->SetRangeUser(0.75, 1.25);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
        c->SaveAs( Form("%s/%s_%s_dijetEtaFBRatio_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

        //
        // MC (reco)
        //

        // Eta full (lab frame)
        hRecoMcDijetEtaLab = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaLab->ProjectionY( Form("hRecoMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoMcDijetEtaLab ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEta1DLab_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaLab->SetName( Form("hRecoMcDijetEtaLab_%d", i) );
        set1DStyle( hRecoMcDijetEtaLab, 0 );
        rescaleHisto1D( hRecoMcDijetEtaLab );

        // Eta full (CM frame)
        hRecoMcDijetEtaCM = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaCM->ProjectionY( Form("hRecoMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoMcDijetEtaCM ) {
            std::cerr << Form("MC histogram not found: hRecoMcDijetEtaCM_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaCM->SetName( Form("hRecoMcDijetEtaCM_%d", i) );
        set1DStyle( hRecoMcDijetEtaCM, 0 );
        rescaleHisto1D( hRecoMcDijetEtaCM );

        // Eta forward (CM frame)
        hRecoMcDijetEtaForward = dynamic_cast<TH1D*>(hRecoMcDijetPtEtaCMForward->ProjectionY( Form("hRecoMcDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoMcDijetEtaForward ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEtaCMForward1D_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaForward->SetName( Form("hRecoMcDijetEtaForward_%d", i) );
        set1DStyle( hRecoMcDijetEtaForward, 0 );

        // Eta backward (CM frame)
        hRecoMcDijetEtaBackward = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaCMBackward->ProjectionY( Form("hRecoMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoMcDijetEtaBackward ) {
            std::cerr << Form("MC histogram not found: hRecoMcDijetEtaBackward_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaBackward->SetName( Form("hRecoMcDijetEtaBackward_%d", i) );
        set1DStyle( hRecoMcDijetEtaBackward, 0 );

        // Normalize the forward and backward distributions
        // integralForward = hRecoMcDijetEtaForward->Integral();
        // integralBackward = hRecoMcDijetEtaBackward->Integral();
        // hRecoMcDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
        // hRecoMcDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hRecoMcDijetFBRatio = dynamic_cast<TH1D*>( hRecoMcDijetEtaForward->Clone( Form("hRecoMcDijetFBRatio_%d", i) ) );
        hRecoMcDijetFBRatio->SetName( Form("hRecoMcDijetFBRatio_%d", i) );
        hRecoMcDijetFBRatio->Divide( hRecoMcDijetFBRatio, hRecoMcDijetEtaBackward, 1., 1. );

        //
        // MC (gen)
        //

        // Eta full (lab frame)
        hGenMcDijetEtaLab = dynamic_cast<TH1D*>( hGenMcDijetPtEtaLab->ProjectionY( Form("hGenMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hGenMcDijetEtaLab ) {
            std::cerr << Form("MC histogram not found: hGenMcDijetEtaLab_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaLab->SetName( Form("hGenMcDijetEtaLab_%d", i) );
        set1DStyle( hGenMcDijetEtaLab, 5 );
        rescaleHisto1D( hGenMcDijetEtaLab );    

        // Eta full (CM frame)
        hGenMcDijetEtaCM = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCM->ProjectionY( Form("hGenMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hGenMcDijetEtaCM ) {
            std::cerr << Form("MC histogram not found: hGenMcDijetEtaCM_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaCM->SetName( Form("hGenMcDijetEtaCM_%d", i) );
        set1DStyle( hGenMcDijetEtaCM, 5 );
        rescaleHisto1D( hGenMcDijetEtaCM );
        
        // Eta forward (CM frame)
        hGenMcDijetEtaForward = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCMForward->ProjectionY( Form("hGenMcDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hGenMcDijetEtaForward ) {
            std::cerr << Form("MC histogram not found: hGenMcDijetEtaForward_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaForward->SetName( Form("hGenMcDijetEtaForward_%d", i) );
        set1DStyle( hGenMcDijetEtaForward, 5 );

        // Eta backward (CM frame)
        hGenMcDijetEtaBackward = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCMBackward->ProjectionY( Form("hGenMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hGenMcDijetEtaBackward ) {
            std::cerr << Form("MC histogram not found: hGenMcDijetEtaBackward_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaBackward->SetName( Form("hGenMcDijetEtaBackward_%d", i) );
        set1DStyle( hGenMcDijetEtaBackward, 5 );
        
        // Normalize the forward and backward distributions
        // integralForward = hGenMcDijetEtaForward->Integral();
        // integralBackward = hGenMcDijetEtaBackward->Integral();
        // hGenMcDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
        // hGenMcDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hGenMcDijetFBRatio = dynamic_cast<TH1D*>( hGenMcDijetEtaForward->Clone( Form("hGenMcDijetFBRatio_%d", i) ) );
        hGenMcDijetFBRatio->SetName( Form("hGenMcDijetFBRatio_%d", i) );
        hGenMcDijetFBRatio->Divide( hGenMcDijetFBRatio, hGenMcDijetEtaBackward, 1., 1. );

        //
        // MC (ref)
        //

        // Eta full (lab frame)
        hRefMcDijetEtaLab = dynamic_cast<TH1D*>( hRefMcDijetPtEtaLab->ProjectionY( Form("hRefMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRefMcDijetEtaLab ) {
            std::cerr << Form("MC histogram not found: hRefMcDijetEtaLab_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaLab->SetName( Form("hRefMcDijetEtaLab_%d", i) );
        set1DStyle( hRefMcDijetEtaLab, 1 );
        rescaleHisto1D( hRefMcDijetEtaLab );

        // Eta full (CM frame)
        hRefMcDijetEtaCM = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCM->ProjectionY( Form("hRefMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRefMcDijetEtaCM ) {
            std::cerr << Form("MC histogram not found: hRefMcDijetEtaCM_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaCM->SetName( Form("hRefMcDijetEtaCM_%d", i) );
        set1DStyle( hRefMcDijetEtaCM, 1 );
        rescaleHisto1D( hRefMcDijetEtaCM );

        // Eta forward (CM frame)
        hRefMcDijetEtaForward = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCMForward->ProjectionY( Form("hRefMcDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRefMcDijetEtaForward ) {
            std::cerr << Form("MC histogram not found: hRefMcDijetEtaForward_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaForward->SetName( Form("hRefMcDijetEtaForward_%d", i) );
        set1DStyle( hRefMcDijetEtaForward, 1 );
        

        // Eta backward (CM frame)
        hRefMcDijetEtaBackward = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCMBackward->ProjectionY( Form("hRefMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRefMcDijetEtaBackward ) {
            std::cerr << Form("MC histogram not found: hRefMcDijetEtaBackward_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaBackward->SetName( Form("hRefMcDijetEtaBackward_%d", i) );
        set1DStyle( hRefMcDijetEtaBackward, 1 );

        // Normalize the forward and backward distributions
        // integralForward = hRefMcDijetEtaForward->Integral();
        // integralBackward = hRefMcDijetEtaBackward->Integral();
        // hRefMcDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
        // hRefMcDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

        // Ratios
        hRefMcDijetFBRatio = dynamic_cast<TH1D*>( hRefMcDijetEtaForward->Clone( Form("hRefMcDijetFBRatio_%d", i) ) );
        hRefMcDijetFBRatio->SetName( Form("hRefMcDijetFBRatio_%d", i) );
        hRefMcDijetFBRatio->Divide( hRefMcDijetFBRatio, hRefMcDijetEtaBackward, 1., 1. );

        //
        // Ratios of full distributions to gen in the laboratory 
        //
        hRecoData2GenDijetEtaLab = dynamic_cast<TH1D*>( hRecoDataDijetEtaLab->Clone( Form("hRecoData2GenDijetEtaLab_%d", i) ) );
        hRecoData2GenDijetEtaLab->SetName( Form("hRecoData2GenDijetEtaLab_%d", i) );
        hRecoData2GenDijetEtaLab->GetYaxis()->SetTitle( "Ratio to Gen" );
        hRecoData2GenDijetEtaLab->Divide( hRecoData2GenDijetEtaLab, hGenMcDijetEtaLab, 1., 1. );

        hRecoMc2GenDijetEtaLab = dynamic_cast<TH1D*>( hRecoMcDijetEtaLab->Clone( Form("hRecoMc2GenDijetEtaLab_%d", i) ) );
        hRecoMc2GenDijetEtaLab->SetName( Form("hRecoMc2GenDijetEtaLab_%d", i) );
        hRecoMc2GenDijetEtaLab->GetYaxis()->SetTitle( "Ratio to Gen" );
        hRecoMc2GenDijetEtaLab->Divide( hRecoMc2GenDijetEtaLab, hGenMcDijetEtaLab, 1., 1., "b" );

        hRef2GenDijetEtaLab = dynamic_cast<TH1D*>( hRefMcDijetEtaLab->Clone( Form("hRef2GenDijetEtaLab_%d", i) ) );
        hRef2GenDijetEtaLab->SetName( Form("hRef2GenDijetEtaLab_%d", i) );
        hRef2GenDijetEtaLab->GetYaxis()->SetTitle( "Ratio to Gen" );
        hRef2GenDijetEtaLab->Divide( hRef2GenDijetEtaLab, hGenMcDijetEtaLab, 1., 1., "b" );

        //
        // Ratios of full distributions to gen in the center-of-mass frame
        //
        hRecoData2GenDijetEtaCM = dynamic_cast<TH1D*>( hRecoDataDijetEtaCM->Clone( Form("hRecoData2GenDijetEtaCM_%d", i) ) );
        hRecoData2GenDijetEtaCM->SetName( Form("hRecoData2GenDijetEtaCM_%d", i) );
        hRecoData2GenDijetEtaCM->GetYaxis()->SetTitle( "Ratio to Gen" );
        hRecoData2GenDijetEtaCM->Divide( hRecoData2GenDijetEtaCM, hGenMcDijetEtaCM, 1., 1. );

        hRecoMc2GenDijetEtaCM = dynamic_cast<TH1D*>( hRecoMcDijetEtaCM->Clone( Form("hRecoMc2GenDijetEtaCM_%d", i) ) );
        hRecoMc2GenDijetEtaCM->SetName( Form("hRecoMc2GenDijetEtaCM_%d", i) );
        hRecoMc2GenDijetEtaCM->GetYaxis()->SetTitle( "Ratio to Gen" );
        hRecoMc2GenDijetEtaCM->Divide( hRecoMc2GenDijetEtaCM, hGenMcDijetEtaCM, 1., 1., "b" );

        hRef2GenDijetEtaCM = dynamic_cast<TH1D*>( hRefMcDijetEtaCM->Clone( Form("hRef2GenDijetEtaCM_%d", i) ) );
        hRef2GenDijetEtaCM->SetName( Form("hRef2GenDijetEtaCM_%d", i) );
        hRef2GenDijetEtaCM->GetYaxis()->SetTitle( "Ratio to Gen" );
        hRef2GenDijetEtaCM->Divide( hRef2GenDijetEtaCM, hGenMcDijetEtaCM, 1., 1., "b" );

        //
        // Plot comparisons
        //


        // Lab frame
        c->cd();
        drawDijetToGenComparison(c, hRecoMcDijetEtaLab, hGenMcDijetEtaLab, hRefMcDijetEtaLab, nullptr, hRecoDataDijetEtaLab,
                            (int)ptLow, (int)ptHi, false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );


        // CM frame
        drawDijetToGenComparison(c, hRecoMcDijetEtaCM, hGenMcDijetEtaCM, hRefMcDijetEtaCM, nullptr, hRecoDataDijetEtaCM,
                            (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );



        // Forward-backward ratio
        drawDijetToGenComparison(c, hRecoMc2GenDijetEtaCM, nullptr, hRef2GenDijetEtaCM, nullptr, hRecoData2GenDijetEtaCM, 
                            (int)ptLow, (int)ptHi, true, true, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaFBRatio_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

        //
        // Plot ratios
        //

        // Lab frame
        drawDijetToGenRatio(c, hRecoMc2GenDijetEtaLab, hRef2GenDijetEtaLab, nullptr, hRecoData2GenDijetEtaLab,
                       (int)ptLow, (int)ptHi, false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_DataRecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

        // CM frame
        drawDijetToGenRatio(c, hRecoMc2GenDijetEtaCM, hRef2GenDijetEtaCM, nullptr, hRecoData2GenDijetEtaCM,
                       (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_DataRecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );
      
    } // for (int i = 0; i < sizeOfDijetPtVals - 1; ++i)

    if (c) { delete c; c = nullptr; }
}

//________________
// Plot eta and pT distributions of inclusive jets from Pythia and embedding.
// Compare dijet distributions.
void pythia2embeddingSingleJet(TFile *fEmbedding, TFile *fPythia, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250305") {
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", (int)(1000 * collisionEnergy));
    collSystemStr += Form("_emb2pythia");

    // Event ptHat binning
    double ptHatVals[] { 15, 45, 90 };
    int sizeOfPtHatVals = sizeof(ptHatVals)/sizeof(ptHatVals[0]);

    // Jet pT binning
    double jetPtVals[] = {30., 50., 90., 120., 200., 500., 1500.};
    int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

    // Jet eta binning
    // double jetEtaVals[] = {-3.0, -2.4, -1.6, 1.6, 2.4, 3.0};
    // int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);

    //
    // Retrieve inclusive jet histograms
    //

    // Embedding
    TH3D *hRecoPtEtaPtHatEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRecoInclusiveJetPtEtaPtHat"));
    if (!hRecoPtEtaPtHatEmb) {
        std::cerr << "Error: hRecoPtEtaPtHatEmb not found!" << std::endl;
        return;
    }
    hRecoPtEtaPtHatEmb->SetName("hRecoPtEtaPtHatEmb");
    TH3D *hGenPtEtaPtHatEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hGenInclusiveJetPtEtaPtHat"));
    if (!hGenPtEtaPtHatEmb) {
        std::cerr << "Error: hGenPtEtaPtHatEmb not found!" << std::endl;
        return;
    }
    hGenPtEtaPtHatEmb->SetName("hGenPtEtaPtHatEmb");

    TH1D *hRecoEtaEmb{nullptr};
    TH1D *hGenEtaEmb{nullptr};
    TH1D *hRecoPtEmb{nullptr};
    TH1D *hGenPtEmb{nullptr};

    // Pythia
    TH3D *hRecoPtEtaPtHatPythia = dynamic_cast<TH3D *>(fPythia->Get("hRecoInclusiveJetPtEtaPtHat"));
    if (!hRecoPtEtaPtHatPythia) {
        std::cerr << "Error: hRecoPtEtaPtHatPythia not found!" << std::endl;
        return;
    }
    hRecoPtEtaPtHatPythia->SetName("hRecoPtEtaPtHatPythia");
    TH3D *hGenPtEtaPtHatPythia = dynamic_cast<TH3D *>(fPythia->Get("hGenInclusiveJetPtEtaPtHat"));
    if (!hGenPtEtaPtHatPythia) {
        std::cerr << "Error: hGenPtEtaPtHatPythia not found!" << std::endl;
        return;
    }
    hGenPtEtaPtHatPythia->SetName("hGenPtEtaPtHatPythia");

    TH1D *hRecoEtaPythia{nullptr};
    TH1D *hGenEtaPythia{nullptr};
    TH1D *hRecoPtPythia{nullptr};
    TH1D *hGenPtPythia{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();

    // Loop over ptHat bins
    for (size_t iPtHat = 0; iPtHat < sizeOfPtHatVals-1; ++iPtHat) {

        double ptHatLow = ptHatVals[iPtHat];
        double ptHatHigh = hRecoPtEtaPtHatEmb->GetZaxis()->GetBinUpEdge( hRecoPtEtaPtHatEmb->GetNbinsZ() );
        int ptHatBinLow = hRecoPtEtaPtHatEmb->GetZaxis()->FindBin(ptHatLow);
        int ptHatBinHigh = hRecoPtEtaPtHatEmb->GetZaxis()->FindBin(ptHatHigh)-1;

        // Loop over ptHat bins and make projections of eta
        for (size_t iJetPt = 0; iJetPt < sizeOfJetPtVals-1; ++iJetPt) {

            double jetPtLow = jetPtVals[iJetPt];
            double jetPtHigh = jetPtVals[iJetPt+1];
            int jetPtBinLow = hRecoPtEtaPtHatEmb->GetYaxis()->FindBin(jetPtLow);
            int jetPtBinHigh = hRecoPtEtaPtHatEmb->GetYaxis()->FindBin(jetPtHigh)-1;

            // Embedding
            hRecoEtaEmb = dynamic_cast<TH1D *>(hRecoPtEtaPtHatEmb->ProjectionX( Form("hRecoEtaEmb_%zu_%zu", iPtHat, iJetPt), 
                                                                                jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
            rescaleHisto1D(hRecoEtaEmb);
            set1DStyle(hRecoEtaEmb, 0);

            hGenEtaEmb = dynamic_cast<TH1D *>(hGenPtEtaPtHatEmb->ProjectionX( Form("hGenEtaEmb_%zu_%zu", iPtHat, iJetPt), 
                                                                                    jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
            rescaleHisto1D(hGenEtaEmb);
            set1DStyle(hGenEtaEmb, 0);

            // Pythia
            hRecoEtaPythia = dynamic_cast<TH1D *>(hRecoPtEtaPtHatPythia->ProjectionX(Form("hRecoEtaPythia_%zu_%zu", iPtHat, iJetPt), 
                                                                                          jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
            rescaleHisto1D(hRecoEtaPythia);
            set1DStyle(hRecoEtaPythia, 1);

            hGenEtaPythia = dynamic_cast<TH1D *>(hGenPtEtaPtHatPythia->ProjectionX(Form("hGenEtaPythia_%zu_%zu", iPtHat, iJetPt), 
                                                                                                 jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
            rescaleHisto1D(hGenEtaPythia);
            set1DStyle(hGenEtaPythia, 1);

            // Plot comparisons
            c->cd(); c->Clear(); setPadStyle();
            drawEtaPtComparison(c, hRecoEtaEmb, hRecoEtaPythia, 
                                jetPtLow, jetPtHigh, -5.2, 5.2, ptHatLow, 
                                "Reco (embed)", "Reco (pythia)", collisionSystem, collisionEnergy, 
                                false, true, false);
            c->SaveAs(Form("%s/%s_reco_eta_ptHat_ptHat%d_jetPt_%d_%d.pdf", 
                      date.Data(), collSystemStr.Data(), (int)ptHatLow, (int)jetPtLow, (int)jetPtHigh));

            c->cd(); c->Clear(); setPadStyle();
            drawEtaPtComparison(c, hGenEtaEmb, hGenEtaPythia, 
                                jetPtLow, jetPtHigh, -5.2, 5.2, ptHatLow, 
                                "Gen (embed)", "Gen (pythia)", collisionSystem, collisionEnergy, 
                                false, true, false);
            c->SaveAs(Form("%s/%s_gen_eta_ptHat_ptHat%d_jetPt_%d_%d.pdf", 
                      date.Data(), collSystemStr.Data(), (int)ptHatLow, (int)jetPtLow, (int)jetPtHigh));

        } // for (size_t iJetPt = 0; iJetPt < sizeOfJetPtVals-1; ++iJetPt)
    } // for (size_t iPtHat = 0; iPtHat < sizeOfPtHatVals-1; ++iPtHat)

    if (c) { delete c; c = nullptr; }
}

//________________
// Compare dijet distributions from embedding and Pythia.
// Compare reco, ref, and gen dijet eta distributions.
void pythia2embeddingDijetComparison(TFile *fEmbedding, TFile *fPythia, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250305") {
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", (int)(1000 * collisionEnergy));
    collSystemStr += Form("_emb2pythia");

    double xTextPosition = 0.35;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Determine the direction based on the filename
    TString directionStr;
    TString fileName = fEmbedding->GetName();
    if (fileName.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (fileName.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // Dijet ptAve binning
    // double dijetPtVals[] { 50.,  60.,   70.,  80.,  90.,
    //                           100., 110.,  120., 130., 140.,
    //                           150., 160.,  180., 200., 250., 
    //                           300., 500. };
    int dijetPtVals[] { 50,  90,   120,  180, 200, 250, 300, 350, 500 };
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    // Dijet dijstributions for embedding
    TH3D *hRecoDijetPtEtaCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRecoDijetPtEtaCMWeighted"));
    if (!hRecoDijetPtEtaCMEmb) {
        std::cerr << "Error: hRecoDijetPtEtaCMWeighted for embedding not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hGenDijetPtEtaCMWeighted"));
    if (!hGenDijetPtEtaCMEmb) {
        std::cerr << "Error: hGenDijetPtEtaCMWeighted for embedding not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRefDijetPtEtaCMWeighted"));
    if (!hRefDijetPtEtaCMEmb) {
        std::cerr << "Error: hRefDijetPtEtaCMWeighted for embedding not found in file." << std::endl;
        return;
    }

    // Dijet dijstributions for pythia
    TH3D *hRecoDijetPtEtaCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hRecoDijetPtEtaCMWeighted"));
    if (!hRecoDijetPtEtaCMPythia) {
        std::cerr << "Error: hRecoDijetPtEtaCMWeighted for pythia not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hGenDijetPtEtaCMWeighted"));
    if (!hGenDijetPtEtaCMPythia) {
        std::cerr << "Error: hGenDijetPtEtaCMWeighted for pythia not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hRefDijetPtEtaCMWeighted"));
    if (!hRefDijetPtEtaCMPythia) {
        std::cerr << "Error: hRefDijetPtEtaCMWeighted for pythia not found in file." << std::endl;
        return;
    }

    // Create 1D histograms for dijet eta distributions
    TH1D *hRecoDijetEtaCMEmb{nullptr};
    TH1D *hGenDijetEtaCMEmb{nullptr};
    TH1D *hRefDijetEtaCMEmb{nullptr};

    TH1D *hRecoDijetEtaCMPythia{nullptr};
    TH1D *hGenDijetEtaCMPythia{nullptr};
    TH1D *hRefDijetEtaCMPythia{nullptr};

    TH1D *hRecoDijetEtaCM_emb2pythia{nullptr};
    TH1D *hGenDijetEtaCM_emb2pythia{nullptr};
    TH1D *hRefDijetEtaCM_emb2pythia{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    TLegend *leg{nullptr};

    // Loop over dijet ptAve bins
    for (int i = 0; i < sizeOfPtVals - 1; ++i) {
        double ptLow = dijetPtVals[i];
        double ptHi = dijetPtVals[i + 1];
        int ptBinLow = hRecoDijetPtEtaCMEmb->GetXaxis()->FindBin(ptLow);
        int ptBinHi = hRecoDijetPtEtaCMEmb->GetXaxis()->FindBin(ptHi)-1;

        // Create 1D histograms for embedding
        hRecoDijetEtaCMEmb = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMEmb->ProjectionY(Form("hRecoDijetEtaCMEmb_%d", i), 
                                                  ptBinLow, ptBinHi));
        if (!hRecoDijetEtaCMEmb) {
            std::cerr << "Error: hRecoDijetEtaCMEmb_" << i << " not found." << std::endl;
            return;
        }
        set1DStyle(hRecoDijetEtaCMEmb, 0);
        rescaleHisto1D(hRecoDijetEtaCMEmb);
        hGenDijetEtaCMEmb = dynamic_cast<TH1D *>(hGenDijetPtEtaCMEmb->ProjectionY(Form("hGenDijetEtaCMEmb_%d", i), 
                                                 ptBinLow, ptBinHi));
        if (!hGenDijetEtaCMEmb) {
            std::cerr << "Error: hGenDijetEtaCMEmb_" << i << " not found." << std::endl;
            return;
        }
        set1DStyle(hGenDijetEtaCMEmb, 0);
        rescaleHisto1D(hGenDijetEtaCMEmb);
        hRefDijetEtaCMEmb = dynamic_cast<TH1D *>(hRefDijetPtEtaCMEmb->ProjectionY(Form("hRefDijetEtaCMEmb_%d", i),
                                                  ptBinLow, ptBinHi));
        if (!hRefDijetEtaCMEmb) {
            std::cerr << "Error: hRefDijetEtaCMEmb_" << i << " not found." << std::endl;
            return;
        }
        set1DStyle(hRefDijetEtaCMEmb, 0);
        rescaleHisto1D(hRefDijetEtaCMEmb);

        // Create 1D histograms for pythia
        hRecoDijetEtaCMPythia = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMPythia->ProjectionY(Form("hRecoDijetEtaCMPythia_%d", i), 
                                                     ptBinLow, ptBinHi));
        if (!hRecoDijetEtaCMPythia) {
            std::cerr << "Error: hRecoDijetEtaCMPythia_" << i << " not found." << std::endl;
            return;
        }
        set1DStyle(hRecoDijetEtaCMPythia, 1);
        rescaleHisto1D(hRecoDijetEtaCMPythia);
        hGenDijetEtaCMPythia = dynamic_cast<TH1D *>(hGenDijetPtEtaCMPythia->ProjectionY(Form("hGenDijetEtaCMPythia_%d", i), 
                                                  ptBinLow, ptBinHi));
        if (!hGenDijetEtaCMPythia) {
            std::cerr << "Error: hGenDijetEtaCMPythia_" << i << " not found." << std::endl;
            return;
        }
        set1DStyle(hGenDijetEtaCMPythia, 1);
        rescaleHisto1D(hGenDijetEtaCMPythia);
        hRefDijetEtaCMPythia = dynamic_cast<TH1D *>(hRefDijetPtEtaCMPythia->ProjectionY(Form("hRefDijetEtaCMPythia_%d", i), 
                                                  ptBinLow, ptBinHi));
        if (!hRefDijetEtaCMPythia) {
            std::cerr << "Error: hRefDijetEtaCMPythia_" << i << " not found." << std::endl;
            return;
        }
        set1DStyle(hRefDijetEtaCMPythia, 1);
        rescaleHisto1D(hRefDijetEtaCMPythia);

        // Draw and save the comparisons
        c->cd(); c->Clear(); setPadStyle();
        drawDijetToGenComparison(c, hRecoDijetEtaCMEmb, hRecoDijetEtaCMPythia, nullptr, nullptr, nullptr,
                                 (int)ptLow, (int)ptHi, 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_recoComp_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));

        c->cd(); c->Clear(); setPadStyle();
        drawDijetToGenComparison(c, hGenDijetEtaCMEmb, hGenDijetEtaCMPythia, nullptr, nullptr, nullptr,
                                 (int)ptLow, (int)ptHi,
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_genComp_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));

        c->cd(); c->Clear(); setPadStyle();
        drawDijetToGenComparison(c, hRefDijetEtaCMEmb, hRefDijetEtaCMPythia, nullptr, nullptr, nullptr,
                                 (int)ptLow, (int)ptHi,
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_refComp_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));

        // Create 1D histograms for embedding to pythia ratios
        hRecoDijetEtaCM_emb2pythia = dynamic_cast<TH1D *>(hRecoDijetEtaCMEmb->Clone(Form("hRecoDijetEtaCM_emb2pythia_%d", i)));
        if (!hRecoDijetEtaCM_emb2pythia) {
            std::cerr << "Error: hRecoDijetEtaCM_emb2pythia_" << i << " not found." << std::endl;
            return;
        }
        hRecoDijetEtaCM_emb2pythia->Divide(hRecoDijetEtaCMEmb, hRecoDijetEtaCMPythia, 1., 1., "b");
        set1DStyle(hRecoDijetEtaCM_emb2pythia, 0);
        
        hGenDijetEtaCM_emb2pythia = dynamic_cast<TH1D *>(hGenDijetEtaCMEmb->Clone(Form("hGenDijetEtaCM_emb2pythia_%d", i)));
        if (!hGenDijetEtaCM_emb2pythia) {
            std::cerr << "Error: hGenDijetEtaCM_emb2pythia_" << i << " not found." << std::endl;
            return;
        }
        hGenDijetEtaCM_emb2pythia->Divide(hGenDijetEtaCMEmb, hGenDijetEtaCMPythia, 1., 1., "b");
        set1DStyle(hGenDijetEtaCM_emb2pythia, 4);
        
        hRefDijetEtaCM_emb2pythia = dynamic_cast<TH1D *>(hRefDijetEtaCMEmb->Clone(Form("hRefDijetEtaCM_emb2pythia_%d", i)));
        if (!hRefDijetEtaCM_emb2pythia) {
            std::cerr << "Error: hRefDijetEtaCM_emb2pythia_" << i << " not found." << std::endl;
            return;


    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
        hRefDijetEtaCM_emb2pythia->Divide(hRefDijetEtaCMEmb, hRefDijetEtaCMPythia, 1., 1., "b");
        set1DStyle(hRefDijetEtaCM_emb2pythia, 2);

        // Draw and save the ratios
        c->cd(); c->Clear(); setPadStyle();
        hRecoDijetEtaCM_emb2pythia->Draw();
        hGenDijetEtaCM_emb2pythia->Draw("same");
        hRefDijetEtaCM_emb2pythia->Draw("same");
        hRecoDijetEtaCM_emb2pythia->GetXaxis()->SetRangeUser(-2.4, 2.4);
        hRecoDijetEtaCM_emb2pythia->GetYaxis()->SetRangeUser(0.8, 1.2);
        hRecoDijetEtaCM_emb2pythia->GetYaxis()->SetTitle("Embedding/Pythia");
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%d < p_{T}^{ave} < %d (GeV)", (int)ptLow, (int)ptHi));
        leg = new TLegend(0.35, 0.25, 0.65, 0.4);
        leg->SetTextSize(0.04);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hRecoDijetEtaCM_emb2pythia, "Reco", "p");
        leg->AddEntry(hGenDijetEtaCM_emb2pythia, "Gen", "p");
        leg->AddEntry(hRefDijetEtaCM_emb2pythia, "Ref", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_ratio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));
    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
    if (c) { delete c; c = nullptr; }
}


//________________
// Plot dijet distributions for data by projecting those from 2D (ptAve, eta) histogram
void dijetDistributionsForData(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    TString triggerName;
    if (filename.Contains("Jet100", TString::kIgnoreCase)) {
        triggerName = "Jet100";
    } else if (filename.Contains("Jet80", TString::kIgnoreCase)) {
        triggerName = "Jet80";
    } else if (filename.Contains("Jet60", TString::kIgnoreCase)) {
        triggerName = "Jet60";
    } else {
        triggerName = "MB";
    }

    // Dijet ptAve binning
    // int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
    //                          100, 110,  120, 130, 140,
    //                          150, 160,  180, 200, 250, 
    //                          300, 500 };
    double dijetPtNewVals[] {  50., 80., 120., 180., 200., 250., 300., 350., 500. };
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);



    // Retrieve 2D distributions from data
    TH2D *hRecoDijetPtEtaWeighted = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaWeighted"));
    if (!hRecoDijetPtEtaWeighted) {
        std::cerr << "Error: hRecoDijetPtEtaWeighted not found in file." << std::endl;
        return;
    }
    TH2D *hRecoDijetPtEtaCMWeighted = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMWeighted"));
    if (!hRecoDijetPtEtaCMWeighted) {
        std::cerr << "Error: hRecoDijetPtEtaCMWeighted not found in file." << std::endl;
        return;
    }

}

//________________
void usefullDistributions(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // Dijet ptAve binning
    // int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
    //                          100, 110,  120, 130, 140,
    //                          150, 160,  180, 200, 250, 
    //                          300, 500 };
    int dijetPtNewVals[9] {  50,  90,   120,  180, 200, 250, 300, 350, 500 };
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);


}

//________________
void plotSingleJetData2McComparison(TFile *fExp, TFile *fMc, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    // TString directionStr;
    // TString filename = f->GetName();
    // if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
    //     directionStr = "Pbgoing";
    // } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
    //     directionStr = "pgoing";
    // } else {
    //     directionStr = "combined";
    // }
}

//_________________
void plotForwardBackwardRatios(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    TH2D *hForward2D = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMForwardWeighted"));
    if (!hForward2D) {
        std::cerr << "Error: hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    hForward2D->SetName("hForward2D");
    hForward2D->SetDirectory(0);
    TH2D *hBackward2D = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMBackwardWeighted"));
    if (!hBackward2D) {
        std::cerr << "Error: hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    hBackward2D->SetName("hBackward2D");
    hBackward2D->SetDirectory(0);


    // int dijetPtVals[] { 50, 80, 100, 140 };
    int dijetPtVals[] { 120, 150, 180, 250, 500 };
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);
    for (int i = 0; i < sizeOfPtVals - 1; ++i) {
        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i + 1];
        int ptBinLow = hForward2D->GetXaxis()->FindBin(ptLow);
        int ptBinHi = hForward2D->GetXaxis()->FindBin(ptHi) - 1;


        // Project to 1D
        TH1D *hForward1D = dynamic_cast<TH1D *>( hForward2D->ProjectionY(Form("hForward1D_%d_%d", ptLow, ptHi)) );
        set1DStyle(hForward1D, 0, true);
        TH1D *hBackward1D = dynamic_cast<TH1D *>( hBackward2D->ProjectionY(Form("hBackward1D_%d_%d", ptLow, ptHi)) );
        set1DStyle(hBackward1D, 1, true);

        // Create ratio
        TH1D *hFwdBwdRatio = dynamic_cast<TH1D *>(hForward1D->Clone(Form("hFwdBwdRatio_%d_%d", ptLow, ptHi)));
        hFwdBwdRatio->Divide(hForward1D, hBackward1D, 1., 1., "b");
        set1DStyle(hFwdBwdRatio, 2);

        TCanvas *c = new TCanvas("c", "c", 800, 800);
        setPadStyle();
        gPad->SetGrid();
        hFwdBwdRatio->Draw();
        hFwdBwdRatio->GetXaxis()->SetRangeUser(0., 2.);
        hFwdBwdRatio->GetYaxis()->SetRangeUser(0.75, 1.25);
        hFwdBwdRatio->GetYaxis()->SetTitle("Forward / Backward");
        TLatex t;
        t.SetTextFont(42);
        t.SetTextSize(0.05);
        t.DrawLatexNDC(0.3, 0.8, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));
        delete c;   
    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
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

    //
    // pPb8160
    //

    int collisionSystem = 1;         // 0 - pp, 1 - pPb, 2 - pPb5020, 3 - pPb8160
    double collisionEnergy = 8.16;   // 8.16 TeV
    int direction = 2;               // 0-p-going, 1-Pb-going, 2 - combined
    TString directionStr = (direction == 0) ? "pgoing" : ((direction == 1) ? "Pbgoing" : "");
    int dataTrigger = 3;               // 0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100
    TString dataStr = (dataTrigger == 0) ? "MB" : ((dataTrigger == 1) ? "Jet60" : ((dataTrigger == 2) ? "Jet80" : ((dataTrigger == 3) ? "Jet100" : "unknownData")));
    TString dataDirectionStr = (direction == 0) ? "Pbgoing" : ((direction == 1) ? "pgoing" : "");
    int jetType = 0; // 0 - inclusive, 1 - lead, 2 - sublead
    int matchType = 0; // 0 - inclusive, 1 - matched, 2 - unmatched
    int cutType = 2; // 0 - no cuts, 1 - trkMax, 2 - jetId
    TString cutTypeStr = (cutType == 0) ? "noCut_" : ((cutType == 1) ? "trkMax_" : ((cutType == 2) ? "jetId_" : "unknownCut"));
    int etaRange = 15;

    //
    // Embedding
    //
    TFile *pPb8160EmbedFile = nullptr;
    if ( direction < 2 ) {
        pPb8160EmbedFile = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_%seta%d.root", 
                                        directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160EmbedFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_%seta%d.root", 
                               directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    else {
        pPb8160EmbedFile = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_%seta%d.root", 
                                        cutTypeStr.Data(), etaRange) );
        if ( !pPb8160EmbedFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_%seta%d.root", 
                               cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/macro/pPb8160_Pbgoing_noTrkMax/oEmbedding_pPb8160_Pbgoing_80.root", uname.Data()) );

    //
    // Pythia
    //
    TFile *pPb8160PythiaFile = nullptr;
    if ( direction < 2 ) {
        pPb8160PythiaFile = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_%seta%d.root", 
                                         directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160PythiaFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_%seta%d.root", 
                              directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    else {
        pPb8160PythiaFile = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_%seta%d.root", 
                                         cutTypeStr.Data(), etaRange) );
        // if ( !pPb8160PythiaFile ) {
        //     std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_%seta%d.root", 
        //                       cutTypeStr.Data(), etaRange) << std::endl;
        //     return;
        // }
    }

    //
    // Data
    //
    TFile *pPb8160DataFile = nullptr;
    if ( direction < 2 ) {
        pPb8160DataFile = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/%s/%s_%s_ak4_%seta%d.root", 
            dataDirectionStr.Data(), dataStr.Data(), dataDirectionStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160DataFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/exp/%s/%s_%s_ak4_%seta%d.root", 
                dataDirectionStr.Data(), dataStr.Data(), dataDirectionStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    else {
        pPb8160DataFile = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/%s_pPb8160_ak4_%seta%d.root", 
            dataStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160DataFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/exp/%s_pPb8160_ak4_%seta%d.root", 
                dataStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }


    // Function plots eta and pT distribution comparisons and ratios of reco, ref, refSel to gen jets
    // in bins of ptHat.
    // It retrieves histograms for reconstructed, generated, and reference jets, and computes the ratios.
    // f is the input TFile containing the histograms.
    // collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // collisionEnergy: energy in TeV (default is 8.16 TeV for pPb)
    // jetType: 0 = Inclusive, 1 = Lead, 2 = SubLead
    // date: date string for saving the plots (default is "20250129")
    // inclusiveJetJECClosures(pPb8160EmbedFile, collisionSystem, collisionEnergy, jetType, matchType, date);

    //
    // Comparison of dijet reco and ref to gen distributions
    //
    // dijetClosures( pPb8160EmbedFile, collisionSystem, collisionEnergy, date );

    //
    // Comparison of dijet reco and ref to gen distributions from 3D histograms
    //
    //dijetClosuresFrom2D( pPb8160EmbedFile, collisionSystem, collisionEnergy, date );
    // dijetClosuresFrom2D( pPb8160PythiaFile, collisionSystem, collisionEnergy, date );

    //
    // Plot comparison of inclusive jet eta distributions to check/validate the JEC
    //
    // data2mcInclusiveJetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of dijet reco and ref to gen distributions
    //
    // data2mcDijetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);
    // data2mcDijetComparison(pPb8160DataFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of inclusive jets and dijets for embedding and PYTHIA
    //
    // pythia2embeddingSingleJet(pPb8160EmbedFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of dijet distributions from embedding and Pythia
    //
    // pythia2embeddingDijetComparison(pPb8160EmbedFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);


    //
    // Plot distributions for the reco dijets
    //
    // recoDijetDistributions(pPb8160EmbedFile, collisionSystem, collisionEnergy, date, true);

    plotForwardBackwardRatios(pPb8160DataFile, collisionSystem, collisionEnergy, date);
}
