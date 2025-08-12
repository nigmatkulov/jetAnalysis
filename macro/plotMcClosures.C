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
// Rescale a 2D histogram by its area-weighted integral
// This function rescales a 2D histogram by its area-weighted integral.
// It calculates the total area-weighted sum of the histogram and then normalizes each bin content and error accordingly.
// The function assumes that the histogram has non-equidistant axes.
void rescaleHisto2D(TH2* h2) {
    if (!h2) {
        std::cerr << "Null histogram pointer!" << std::endl;
        return;
    }

    const int nBinsX = h2->GetNbinsX();
    const int nBinsY = h2->GetNbinsY();
    TAxis* xAxis = h2->GetXaxis();
    TAxis* yAxis = h2->GetYaxis();

    double total = 0.0;

    // First pass: calculate total bin-area-weighted sum
    for (int ix = 1; ix <= nBinsX; ++ix) {
        double dx = xAxis->GetBinWidth(ix);
        for (int iy = 1; iy <= nBinsY; ++iy) {
            double dy = yAxis->GetBinWidth(iy);
            double content = h2->GetBinContent(ix, iy);
            total += content * dx * dy;
        }
    }

    // Second pass: normalize content and error
    if (total > 0) {
        for (int ix = 1; ix <= nBinsX; ++ix) {
            for (int iy = 1; iy <= nBinsY; ++iy) {
                double content = h2->GetBinContent(ix, iy);
                double error   = h2->GetBinError(ix, iy);

                double normContent = content / total;
                double normError   = error / total;

                h2->SetBinContent(ix, iy, normContent);
                h2->SetBinError(ix, iy, normError);
            }
        }
    } else {
        std::cerr << "Warning: total area-normalized sum is zero!" << std::endl;
    }
}

//________________
// Rescale a 3D histogram by its volume-weighted integral
// This function rescales a 3D histogram by its volume-weighted integral.
// It calculates the total volume-weighted sum of the histogram and then normalizes each bin content and error accordingly.
// The function assumes that the histogram has non-equidistant axes.
void rescaleHisto3D(TH3D* h3) {
    if (!h3) return;

    TAxis* xAxis = h3->GetXaxis();
    TAxis* yAxis = h3->GetYaxis();  // Non-equidistant
    TAxis* zAxis = h3->GetZaxis();

    int nX = xAxis->GetNbins();
    int nY = yAxis->GetNbins();
    int nZ = zAxis->GetNbins();

    double total = 0.0;

    // First pass: calculate total volume-weighted sum
    for (int i = 1; i <= nX; ++i) {
        double dx = xAxis->GetBinWidth(i);
        for (int j = 1; j <= nY; ++j) {
            double dy = yAxis->GetBinWidth(j);
            for (int k = 1; k <= nZ; ++k) {
                double dz = zAxis->GetBinWidth(k);
                double content = h3->GetBinContent(i, j, k);
                total += content * dx * dy * dz;
            }
        }
    }

    if (total == 0.0) {
        std::cerr << "Warning: Histogram has zero total volume-weighted integral." << std::endl;
        return;
    }

    // Second pass: normalize bin content and errors
    for (int i = 1; i <= nX; ++i) {
        double dx = xAxis->GetBinWidth(i);
        for (int j = 1; j <= nY; ++j) {
            double dy = yAxis->GetBinWidth(j);
            for (int k = 1; k <= nZ; ++k) {
                double dz = zAxis->GetBinWidth(k);
                double volume = dx * dy * dz;

                double content = h3->GetBinContent(i, j, k);
                double error = h3->GetBinError(i, j, k);

                double normFactor = 1.0 / total;
                double newContent = content * normFactor;
                double newError = error * normFactor;

                h3->SetBinContent(i, j, k, newContent);
                h3->SetBinError(i, j, k, newError);
            }
        }
    }
}

//________________
void checkIntegral(TH1* h) {
    if (h->Integral() <= 0) {
        std::cout << "Warning: " << h->GetName() << " has zero or negative integral." << std::endl;
    }
    if (h->GetEntries() <= 0) {
        std::cout << "Warning: " << h->GetName() << " has zero entries." << std::endl;
    }
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

    double xRange[2] = {0.0, 2.8};
    double yRange[2] = {0.75, 1.25};

    if (isCM) {
        if (!isFB) {
            xRange[0] = 0.0;
            xRange[1] = 2.8;
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
            xRange[0] = -3.0; xRange[1] = 3.0;
            // xRange[0] = -2.5; xRange[1] = 2.5;
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
// Ratio of two histograms with uncertainty propagation
// This function computes the ratio of two histograms A and B, where A is not a subset of B.
// It calculates the ratio R = A / B and propagates the uncertainties using the formula:
// sigma_R = R * sqrt((sigma_A^2 / A^2) + (sigma_B^2 / B^2)), where sigma_A and sigma_B are the
// uncertainties of A and B, respectively. If B is zero, the ratio is set to zero and the uncertainty is also set to zero.
// The function takes three histograms as input: h1 (numerator), h2 (denominator), and hRatio (output histogram for the ratio).
// The function checks for null pointers and ensures that all histograms have the same number of bins.
// It also checks that the denominator histogram (h2) is not zero before performing the ratio calculation.
void computeNonBinomialRatio(TH1D* h1, TH1D* h2, TH1D* hRatio) {
    if (!h1 || !h2 || !hRatio) {
        std::cerr << "Error: Null histogram pointer passed." << std::endl;
        return;
    }

    int nbins = h1->GetNbinsX();
    if (nbins != h2->GetNbinsX() || nbins != hRatio->GetNbinsX()) {
        std::cerr << "Error: Histograms have different number of bins." << std::endl;
        return;
    }

    for (int i = 1; i <= nbins; ++i) {
        double A = h1->GetBinContent(i);
        double B = h2->GetBinContent(i);
        double sigmaA = h1->GetBinError(i);
        double sigmaB = h2->GetBinError(i);

        if (B > 0) {
            double R = A / B;
            double sigmaR = R * std::sqrt(
                (A > 0 ? (sigmaA * sigmaA) / (A * A) : 0) +
                (sigmaB * sigmaB) / (B * B)
            );

            hRatio->SetBinContent(i, R);
            hRatio->SetBinError(i, sigmaR);
        } else {
            // Set zero or sentinel values when B is 0
            hRatio->SetBinContent(i, 0);
            hRatio->SetBinError(i, 0);
        }
    }
}

//________________
// Ratio of two histograms with uncertainty propagation
// This function computes the ratio of two histograms A and B, where A is a subset of B.
// It calculates the ratio R = A / B and propagates the uncertainties using the formula:
// sigma_R = sqrt(R * (1 - R) / B), where R is the ratio and B is the denominator histogram.
// If B is zero, the ratio is set to zero and the uncertainty is also set to zero.
// The function takes three histograms as input: hNumerator (numerator), hDenominator (denominator),
// and hRatio (output histogram for the ratio).
// The function checks for null pointers and ensures that all histograms have the same number of bins.
// It also checks that the denominator histogram (hDenominator) is not zero before performing the ratio calculation.
// The function is specifically designed for cases where the numerator is a subset of the denominator.
// It uses the binomial ratio formula for uncertainty propagation.
void computeBinomialRatio(TH1D* hNumerator, TH1D* hDenominator, TH1D* hRatio) {
    if (!hNumerator || !hDenominator || !hRatio) {
        std::cerr << "Error: Null histogram pointer passed." << std::endl;
        return;
    }

    int nbins = hNumerator->GetNbinsX();
    if (nbins != hDenominator->GetNbinsX() || nbins != hRatio->GetNbinsX()) {
        std::cerr << "Error: Histograms have different number of bins." << std::endl;
        return;
    }

    for (int i = 1; i <= nbins; ++i) {
        double A = hNumerator->GetBinContent(i);
        double B = hDenominator->GetBinContent(i);

        if (B > 0 && A <= B) {
            double R = A / B;
            double sigmaR = std::sqrt(R * (1.0 - R) / B);

            hRatio->SetBinContent(i, R);
            hRatio->SetBinError(i, sigmaR);
        } else {
            hRatio->SetBinContent(i, 0);
            hRatio->SetBinError(i, 0);
        }
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
    int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
                             100, 110,  120, 130, 140,
                             150, 160,  180, 200, 250, 
                             300, 500};
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);


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

        int ptLow = dijetPtNewVals[i];
        int ptHi = dijetPtNewVals[i+1];
        int canvX{1000}, canvY{1000};

        //
        // Lab frame
        //
        hRecoDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1D_%d", i) ) );
        hRecoDijetEta1DLab->SetName( Form("hRecoDijetEta1DLab_%d", i) );
        rescaleHisto1D( hRecoDijetEta1DLab );
        set1DStyle( hRecoDijetEta1DLab, 0 );
        hGenDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1D_%d", i) ) );
        hGenDijetEta1DLab->SetName( Form("hGenDijetEta1DLab_%d", i) );
        rescaleHisto1D( hGenDijetEta1DLab );
        set1DStyle( hGenDijetEta1DLab, 4 );
        hRefDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1D_%d", i) ) );
        hRefDijetEta1DLab->SetName( Form("hRefDijetEta1DLab_%d", i) );
        rescaleHisto1D( hRefDijetEta1DLab );
        set1DStyle( hRefDijetEta1DLab, 1 );
        hRefSelDijetEta1DLab = dynamic_cast<TH1D*>( f->Get( Form("hRefSelDijetEta1D_%d", i) ) );
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
        hRecoDijetEta1DCM->SetName( Form("hRecoDijetEta1DCM_%d", i) );
        set1DStyle( hRecoDijetEta1DCM, 0 );
        rescaleHisto1D( hRecoDijetEta1DCM );
        hGenDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1DCM_%d", i) ) );
        hGenDijetEta1DCM->SetName( Form("hGenDijetEta1DCM_%d", i) );
        set1DStyle( hGenDijetEta1DCM, 4 );
        rescaleHisto1D( hGenDijetEta1DCM );
        hRefDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1DCM_%d", i) ) );
        hRefDijetEta1DCM->SetName( Form("hRefDijetEta1DCM_%d", i) );
        set1DStyle( hRefDijetEta1DCM, 1 );
        rescaleHisto1D( hRefDijetEta1DCM );
        hRefSelDijetEta1DCM = dynamic_cast<TH1D*>( f->Get( Form("hRefSelDijetEta1DCM_%d", i) ) );
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
        hRecoDijetEtaCMForward1D->SetName( Form("hRecoDijetEtaCMForward1D_%d", i) );
        hRecoDijetEtaCMBackward1D = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMBackward1D_%d", i) ) );
        hRecoDijetEtaCMBackward1D->SetName( Form("hRecoDijetEtaCMBackward1D_%d", i) );
        hGenDijetEtaCMForward1D = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMForward1D_%d", i) ) );
        hGenDijetEtaCMForward1D->SetName( Form("hGenDijetEtaCMForward1D_%d", i) );
        hGenDijetEtaCMBackward1D = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMBackward1D_%d", i) ) );
        hGenDijetEtaCMBackward1D->SetName( Form("hGenDijetEtaCMBackward1D_%d", i) );
        hRefDijetEtaCMForward1D = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMForward1D_%d", i) ) );
        hRefDijetEtaCMForward1D->SetName( Form("hRefDijetEtaCMForward1D_%d", i) );
        hRefDijetEtaCMBackward1D = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMBackward1D_%d", i) ) );
        hRefDijetEtaCMBackward1D->SetName( Form("hRefDijetEtaCMBackward1D_%d", i) );

        rescaleForwardBackward( hRecoDijetEtaCMForward1D, hRecoDijetEtaCMBackward1D );
        rescaleForwardBackward( hGenDijetEtaCMForward1D, hGenDijetEtaCMBackward1D );
        rescaleForwardBackward( hRefDijetEtaCMForward1D, hRefDijetEtaCMBackward1D );
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
                                 dijetPtNewVals[i], dijetPtNewVals[i+1],
                                 false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        dijetPtNewVals[i], dijetPtNewVals[i+1]) );

        //
        // Plot ratios in the lab frame
        //
        drawDijetToGenRatio(c, hReco2GenDijetEta1DLab, hRef2GenDijetEta1DLab, hRefSel2GenDijetEta1DLab, nullptr,
                            dijetPtNewVals[i], dijetPtNewVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        dijetPtNewVals[i], dijetPtNewVals[i+1]) );

        //
        // Plot comparisons in the CM frame
        //

        drawDijetToGenComparison(c, hRecoDijetEta1DCM, hRefDijetEta1DCM, hGenDijetEta1DCM, hRefSelDijetEta1DCM, nullptr,
                                 dijetPtNewVals[i], dijetPtNewVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        dijetPtNewVals[i], dijetPtNewVals[i+1]) );

        //
        // Plot ratios in the CM frame
        //
        drawDijetToGenRatio(c, hReco2GenDijetEta1DCM, hRef2GenDijetEta1DCM, hRefSel2GenDijetEta1DCM, nullptr,
                            dijetPtNewVals[i], dijetPtNewVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        dijetPtNewVals[i], dijetPtNewVals[i+1]) );

        //
        // Forward/backward ratios
        //
        drawDijetToGenComparison(c, hRecoDijetFBEtaCM1D, hRefDijetFBEtaCM1D, hGenDijetFBEtaCM1D, nullptr, nullptr,
                                 dijetPtNewVals[i], dijetPtNewVals[i+1], true, true, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaFB_RecoRefGenComp_ptAve_%d_%d.pdf", 
                        date.Data(), collSystemStr.Data(), directionStr.Data(), 
                        dijetPtNewVals[i], dijetPtNewVals[i+1]) );

    } // end loop over dijet pt bins
}


//________________
// This function retrieves 3D histograms (eta, pT, phi)  for reconstructed, generated, 
// and reference jets, and computes the ratios (to gen).
// f is the input TFile containing the histograms.
// collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
// collisionEnergy: energy in TeV (default is 8.16 TeV for pPb)
// date: date string for saving the plots (default is "20250129")
void dijetClosuresFrom3D(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    
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
    int dijetPtNewVals[9] { 50,  90,   120,  180, 200, 250, 300, 350, 500 };
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);

    // Laboratory frame
    TH3D *hRecoDijetPtEtaPhiLab = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiWeighted"));
    // TH3D *hRecoDijetPtEtaPhiLab = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiMatched"));
    if (!hRecoDijetPtEtaPhiLab) {
        std::cerr << "Error: hRecoDijetPtEtaPhiWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaPhiLab = dynamic_cast<TH3D *>(f->Get("hGenDijetPtEtaPhiWeighted"));
    if (!hGenDijetPtEtaPhiLab) {
        std::cerr << "Error: hGenDijetPtEtaPhiWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaPhiLab = dynamic_cast<TH3D *>(f->Get("hRefDijetPtEtaPhiWeighted"));
    if (!hRefDijetPtEtaPhiLab) {
        std::cerr << "Error: hRefDijetPtEtaPhiWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hRefSelDijetPtEtaPhiLab = dynamic_cast<TH3D *>(f->Get("hRefSelDijetPtEtaPhiWeighted"));
    if (!hRefSelDijetPtEtaPhiLab) {
        std::cerr << "Error: hRefSelDijetPtEtaPhiWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hRecoDijetPtEtaPhiCMInLab = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiCMInLab"));
    if (!hRecoDijetPtEtaPhiCMInLab) {
        std::cerr << "Error: hRecoDijetPtEtaPhiCMInLab not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaPhiCMInLab = dynamic_cast<TH3D *>(f->Get("hGenDijetPtEtaPhiCMInLab"));
    if (!hGenDijetPtEtaPhiCMInLab) {
        std::cerr << "Error: hGenDijetPtEtaPhiCMInLab not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaPhiCMInLab = dynamic_cast<TH3D *>(f->Get("hRefDijetPtEtaPhiCMInLab"));
    if (!hRefDijetPtEtaPhiCMInLab) {
        std::cerr << "Error: hRefDijetPtEtaPhiCMInLab not found in file." << std::endl;
        return;
    }
    TH3D *hRefSelDijetPtEtaPhiCMInLab = dynamic_cast<TH3D *>(f->Get("hRefSelDijetPtEtaPhiCMInLab"));
    if (!hRefSelDijetPtEtaPhiCMInLab) {
        std::cerr << "Error: hRefSelDijetPtEtaPhiCMInLab not found in file." << std::endl;
        return;
    }

    // Rescale the histograms by their volume-weighted integral
    rescaleHisto3D(hRecoDijetPtEtaPhiLab);
    rescaleHisto3D(hGenDijetPtEtaPhiLab);
    rescaleHisto3D(hRefDijetPtEtaPhiLab);
    rescaleHisto3D(hRefSelDijetPtEtaPhiLab);
    rescaleHisto3D(hRecoDijetPtEtaPhiCMInLab);
    rescaleHisto3D(hGenDijetPtEtaPhiCMInLab);
    rescaleHisto3D(hRefDijetPtEtaPhiCMInLab);
    rescaleHisto3D(hRefSelDijetPtEtaPhiCMInLab);

    // Center of mass frame
    TH3D *hRecoDijetPtEtaPhiCM = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiCMWeighted"));
    // TH3D *hRecoDijetPtEtaPhiCM = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiCMMatched"));
    if (!hRecoDijetPtEtaPhiCM) {
        std::cerr << "Error: hRecoDijetPtEtaPhiCMWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaPhiCM = dynamic_cast<TH3D *>(f->Get("hGenDijetPtEtaPhiCMWeighted"));
    if (!hGenDijetPtEtaPhiCM) {
        std::cerr << "Error: hGenDijetPtEtaPhiCMWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaPhiCM = dynamic_cast<TH3D *>(f->Get("hRefDijetPtEtaPhiCMWeighted"));
    if (!hRefDijetPtEtaPhiCM) {
        std::cerr << "Error: hRefDijetPtEtaPhiCMWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hRefSelDijetPtEtaPhiCM = dynamic_cast<TH3D *>(f->Get("hRefSelDijetPtEtaPhiCMWeighted"));
    if (!hRefSelDijetPtEtaPhiCM) {
        std::cerr << "Error: hRefSelDijetPtEtaPhiCMWeighted not found in file." << std::endl;
        return;
    }
    TH3D *hRecoDijetPtEtaPhiLabInCM = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiLabInCM"));
    if (!hRecoDijetPtEtaPhiLabInCM) {
        std::cerr << "Error: hRecoDijetPtEtaPhiLabInCM not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaPhiLabInCM = dynamic_cast<TH3D *>(f->Get("hGenDijetPtEtaPhiLabInCM"));
    if (!hGenDijetPtEtaPhiLabInCM) {
        std::cerr << "Error: hGenDijetPtEtaPhiLabInCM not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaPhiLabInCM = dynamic_cast<TH3D *>(f->Get("hRefDijetPtEtaPhiLabInCM"));
    if (!hRefDijetPtEtaPhiLabInCM) {
        std::cerr << "Error: hRefDijetPtEtaPhiLabInCM not found in file." << std::endl;
        return;
    }
    TH3D *hRefSelDijetPtEtaPhiLabInCM = dynamic_cast<TH3D *>(f->Get("hRefSelDijetPtEtaPhiLabInCM"));
    if (!hRefSelDijetPtEtaPhiLabInCM) {
        std::cerr << "Error: hRefSelDijetPtEtaPhiLabInCM not found in file." << std::endl;
        return;
    }

    TH2D *hRecoDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMForwardWeighted"));
    if (!hRecoDijetPtEtaCMForward) {
        std::cerr << "Error: hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    TH2D *hGenDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaCMForwardWeighted"));
    if (!hGenDijetPtEtaCMForward) {
        std::cerr << "Error: hGenDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    TH2D *hRefDijetPtEtaCMForward = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaCMForwardWeighted"));
    if (!hRefDijetPtEtaCMForward) {
        std::cerr << "Error: hRefDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
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
    TH2D *hGenDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hGenDijetPtEtaCMBackwardWeighted"));
    if (!hGenDijetPtEtaCMBackward) {
        std::cerr << "Error: hGenDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    TH2D *hRefDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hRefDijetPtEtaCMBackwardWeighted"));
    if (!hRefDijetPtEtaCMBackward) {
        std::cerr << "Error: hRefDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    // TH2D *hRefSelDijetPtEtaCMBackward = dynamic_cast<TH2D *>(f->Get("hRefSelDijetPtEtaCMBackwardWeighted"));
    // if (!hRefSelDijetPtEtaCMBackward) {
    //     std::cerr << "Error: hRefSelDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
    //     return;
    // }

    // Rescale the histograms by their volume-weighted integral
    rescaleHisto3D(hRecoDijetPtEtaPhiCM);
    rescaleHisto3D(hGenDijetPtEtaPhiCM);
    rescaleHisto3D(hRefDijetPtEtaPhiCM);
    rescaleHisto3D(hRefSelDijetPtEtaPhiCM);
    rescaleHisto3D(hRecoDijetPtEtaPhiLabInCM);
    rescaleHisto3D(hGenDijetPtEtaPhiLabInCM);
    rescaleHisto3D(hRefDijetPtEtaPhiLabInCM);
    rescaleHisto3D(hRefSelDijetPtEtaPhiLabInCM);

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

    TH1D *hRecEtaFBCM{nullptr};
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

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    TLegend *leg{nullptr};
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Loop over dijet ptAve bins
    for (int i = 0; i < sizeOfPtVals - 1; ++i) {
        int ptLow = dijetPtNewVals[i];
        int ptHi = dijetPtNewVals[i + 1];
        int ptBinLow = hRecoDijetPtEtaPhiLab->GetYaxis()->FindBin(ptLow);
        int ptBinHigh = hRecoDijetPtEtaPhiLab->GetYaxis()->FindBin(ptHi) - 1;

        //
        // Laboratory frame
        //

        // Dijet pseudorapidity
        hRecoEtaLab = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiLab->ProjectionY(Form("hRecoEtaLab_%d", i), ptBinLow, ptBinHigh));
        hGenEtaLab = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiLab->ProjectionY(Form("hGenEtaLab_%d", i), ptBinLow, ptBinHigh));
        hRefEtaLab = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiLab->ProjectionY(Form("hRefEtaLab_%d", i), ptBinLow, ptBinHigh));
        hRefSelEtaLab = dynamic_cast<TH1D *>(hRefSelDijetPtEtaPhiLab->ProjectionY(Form("hRefSelEtaLab_%d", i), ptBinLow, ptBinHigh));
        hRecoEtaCMInLab = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiCMInLab->ProjectionY(Form("hRecoEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        hGenEtaCMInLab = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiCMInLab->ProjectionY(Form("hGenEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        hRefEtaCMInLab = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiCMInLab->ProjectionY(Form("hRefEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        hRefSelEtaCMInLab = dynamic_cast<TH1D *>(hRefSelDijetPtEtaPhiCMInLab->ProjectionY(Form("hRefSelEtaCMInLab_%d", i), ptBinLow, ptBinHigh));
        set1DStyle(hRecoEtaLab, 0, true);
        set1DStyle(hGenEtaLab, 5, true);
        set1DStyle(hRefEtaLab, 1, true);
        set1DStyle(hRefSelEtaLab, 2, true);
        set1DStyle(hRecoEtaCMInLab, 0, true);
        set1DStyle(hGenEtaCMInLab, 3, true);
        set1DStyle(hRefEtaCMInLab, 1, true);
        set1DStyle(hRefSelEtaCMInLab, 2, true);

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
        //           << " ptLowBin: " << hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) 
        //           << " ptLow: " << hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) )   
        //           << " ptHiBin: " << hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1
        //           << " ptHi: " << hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1) << std::endl;

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoEtaLab, hRefEtaLab, hGenEtaLab, hRefSelEtaLab, nullptr,
                                 (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) ), 
                                 (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1), 
                                false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        drawDijetToGenComparison(c, hRecoEtaCMInLab, hRefEtaCMInLab, hGenEtaCMInLab, hRefSelEtaCMInLab, nullptr,
                                 (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) ), 
                                 (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1), 
                                 false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCMInLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenEtaLab, hRef2GenEtaLab, hRefSel2GenEtaLab, nullptr,
                            (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) ), 
                            (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1), 
                            false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        drawDijetToGenRatio(c, hReco2GenEtaCMInLab, hRef2GenEtaCMInLab, hRefSel2GenEtaCMInLab, nullptr,
                            (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) ), 
                            (int)hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1), 
                            false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCMInLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        // Dijet azimuthal angle
        hRecoPhiLab = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiLab->ProjectionZ(Form("hRecoPhiLab_%d", i),
                                           hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow), 
                                           hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1));
        hGenPhiLab = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiLab->ProjectionZ(Form("hGenPhiLab_%d", i), 
                                          hGenDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow), 
                                          hGenDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1));
        hRefPhiLab = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiLab->ProjectionZ(Form("hRefPhiLab_%d", i), 
                                          hRefDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow), 
                                          hRefDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1));
        hRefSelPhiLab = dynamic_cast<TH1D *>(hRefSelDijetPtEtaPhiLab->ProjectionZ(Form("hRefSelPhiLab_%d", i), 
                                             hRefSelDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow), 
                                             hRefSelDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1));
        set1DStyle(hRecoPhiLab, 0, true);
        set1DStyle(hGenPhiLab, 5, true);
        set1DStyle(hRefPhiLab, 1, true);
        set1DStyle(hRefSelPhiLab, 2, true);

        // Compute ratios for laboratory frame
        hReco2GenPhiLab = dynamic_cast<TH1D *>(hRecoPhiLab->Clone(Form("hReco2GenPhiLab_%d", i)));
        hReco2GenPhiLab->Divide(hReco2GenPhiLab, hGenPhiLab, 1., 1., "b");
        hRef2GenPhiLab = dynamic_cast<TH1D *>(hRefPhiLab->Clone(Form("hRef2GenPhiLab_%d", i)));
        hRef2GenPhiLab->Divide(hRef2GenPhiLab, hGenPhiLab, 1., 1., "b");
        hRefSel2GenPhiLab = dynamic_cast<TH1D *>(hRefSelPhiLab->Clone(Form("hRefSel2GenPhiLab_%d", i)));
        hRefSel2GenPhiLab->Divide(hRefSel2GenPhiLab, hGenPhiLab, 1., 1., "b");

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoPhiLab, hRefPhiLab, hGenPhiLab, hRefSelPhiLab, nullptr,
                                 hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) ), 
                                 hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1), 
                                 false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiLab_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));
        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenPhiLab, hRef2GenPhiLab, hRefSel2GenPhiLab, nullptr,
                            hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptLow) ), 
                            hRecoDijetPtEtaPhiLab->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiLab->GetXaxis()->FindBin(ptHi)-1), 
                            false, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        
        //
        // Center of mass frame
        //

        // Dijet pseudorapidity
        hRecoEtaCM = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiCM->ProjectionY(Form("hRecoEtaCM_%d", i), 
                                          hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                          hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hGenEtaCM = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiCM->ProjectionY(Form("hGenEtaCM_%d", i), 
                                         hGenDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                         hGenDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hRefEtaCM = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiCM->ProjectionY(Form("hRefEtaCM_%d", i), 
                                         hRefDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                         hRefDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hRefSelEtaCM = dynamic_cast<TH1D *>(hRefSelDijetPtEtaPhiCM->ProjectionY(Form("hRefSelEtaCM_%d", i), 
                                            hRefSelDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                            hRefSelDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hRecoEtaLabInCM = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiLabInCM->ProjectionY(Form("hRecoEtaLabInCM_%d", i), 
                                               hRecoDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptLow), 
                                               hRecoDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptHi)-1));
        hGenEtaLabInCM = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiLabInCM->ProjectionY(Form("hGenEtaLabInCM_%d", i), 
                                              hGenDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptLow), 
                                              hGenDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptHi)-1));
        hRefEtaLabInCM = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiLabInCM->ProjectionY(Form("hRefEtaLabInCM_%d", i), 
                                              hRefDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptLow), 
                                              hRefDijetPtEtaPhiCMInLab->GetXaxis()->FindBin(ptHi)-1));
        hRefSelEtaLabInCM = dynamic_cast<TH1D *>(hRefSelDijetPtEtaPhiLabInCM->ProjectionY(Form("hRefSelEtaLabInCM_%d", i), 
                                                 hRefSelDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptLow), 
                                                 hRefSelDijetPtEtaPhiLabInCM->GetXaxis()->FindBin(ptHi)-1));
                                            
        set1DStyle(hRecoEtaCM, 0, true);
        set1DStyle(hGenEtaCM, 5, true);
        set1DStyle(hRefEtaCM, 1, true);
        set1DStyle(hRefSelEtaCM, 2, true);
        set1DStyle(hRecoEtaLabInCM, 0, true);
        set1DStyle(hGenEtaLabInCM, 3, true);
        set1DStyle(hRefEtaLabInCM, 1, true);
        set1DStyle(hRefSelEtaLabInCM, 2, true);

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
                                 (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow) ), 
                                 (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1), 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        drawDijetToGenComparison(c, hRecoEtaLabInCM, hRefEtaLabInCM, hGenEtaLabInCM, hRefSelEtaLabInCM, nullptr,
                                 (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow) ), 
                                 (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1), 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLabInCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenEtaCM, hRef2GenEtaCM, hRefSel2GenEtaCM, nullptr,
                            (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow) ), 
                            (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1), 
                            true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_RecoRefGenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        drawDijetToGenRatio(c, hReco2GenEtaLabInCM, hRef2GenEtaLabInCM, hRefSel2GenEtaLabInCM, nullptr,
                            (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow) ), 
                            (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1), 
                            true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaLabInCM_RecoRefGenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        // Dijet azimuthal angle
        hRecoPhiCM = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiCM->ProjectionZ(Form("hRecoPhiCM_%d", i), 
                                          hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                          hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hGenPhiCM = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiCM->ProjectionZ(Form("hGenPhiCM_%d", i), 
                                         hGenDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                         hGenDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hRefPhiCM = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiCM->ProjectionZ(Form("hRefPhiCM_%d", i), 
                                         hRefDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                         hRefDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        hRefSelPhiCM = dynamic_cast<TH1D *>(hRefSelDijetPtEtaPhiCM->ProjectionZ(Form("hRefSelPhiCM_%d", i), 
                                              hRefSelDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow), 
                                              hRefSelDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1));
        set1DStyle(hRecoPhiCM, 0, true);
        set1DStyle(hGenPhiCM, 5, true);
        set1DStyle(hRefPhiCM, 1, true);
        set1DStyle(hRefSelPhiCM, 2, true);

        // Compute ratios for center of mass frame
        hReco2GenPhiCM = dynamic_cast<TH1D *>(hRecoPhiCM->Clone(Form("hReco2GenPhiCM_%d", i)));
        hReco2GenPhiCM->Divide(hReco2GenPhiCM, hGenPhiCM, 1., 1., "b");
        hRef2GenPhiCM = dynamic_cast<TH1D *>(hRefPhiCM->Clone(Form("hRef2GenPhiCM_%d", i)));
        hRef2GenPhiCM->Divide(hRef2GenPhiCM, hGenPhiCM, 1., 1., "b");
        hRefSel2GenPhiCM = dynamic_cast<TH1D *>(hRefSelPhiCM->Clone(Form("hRefSel2GenPhiCM_%d", i)));
        hRefSel2GenPhiCM->Divide(hRefSel2GenPhiCM, hGenPhiCM, 1., 1., "b");

        // Plot comparisons
        drawDijetToGenComparison(c, hRecoPhiCM, hRefPhiCM, hGenPhiCM, hRefSelPhiCM, nullptr,
                                 (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow) ), 
                                 (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1), 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiCM_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));
        // Plot ratios
        drawDijetToGenRatio(c, hReco2GenPhiCM, hRef2GenPhiCM, hRefSel2GenPhiCM, nullptr,
                            (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptLow) ), 
                            (int)hRecoDijetPtEtaPhiCM->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaPhiCM->GetXaxis()->FindBin(ptHi)-1), 
                            true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetPhiCM_RecoRef2GenRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));


        //
        // Forward and backward dijet pseudorapidity
        //
        hRecoEtaCMForward = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMForward->ProjectionY(Form("hRecoEtaCMForward_%d", i), 
                                                hRecoDijetPtEtaCMForward->GetXaxis()->FindBin(ptLow), 
                                                hRecoDijetPtEtaCMForward->GetXaxis()->FindBin(ptHi)-1));
        hGenEtaCMForward = dynamic_cast<TH1D *>(hGenDijetPtEtaCMForward->ProjectionY(Form("hGenEtaCMForward_%d", i), 
                                             hGenDijetPtEtaCMForward->GetXaxis()->FindBin(ptLow), 
                                             hGenDijetPtEtaCMForward->GetXaxis()->FindBin(ptHi)-1));
        hRefEtaCMForward = dynamic_cast<TH1D *>(hRefDijetPtEtaCMForward->ProjectionY(Form("hRefEtaCMForward_%d", i), 
                                             hRefDijetPtEtaCMForward->GetXaxis()->FindBin(ptLow), 
                                             hRefDijetPtEtaCMForward->GetXaxis()->FindBin(ptHi)-1));
        // hRefSelEtaCMForward = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCMForward->ProjectionY(Form("hRefSelEtaCMForward_%d", i), 
        //                                           hRefSelDijetPtEtaCMForward->GetXaxis()->FindBin(ptLow), 
        //                                           hRefSelDijetPtEtaCMForward->GetXaxis()->FindBin(ptHi)-1));

        hRecoEtaCMBackward = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMBackward->ProjectionY(Form("hRecoEtaCMBackward_%d", i),
                                                  hRecoDijetPtEtaCMBackward->GetXaxis()->FindBin(ptLow), 
                                                  hRecoDijetPtEtaCMBackward->GetXaxis()->FindBin(ptHi)-1));
        hGenEtaCMBackward = dynamic_cast<TH1D *>(hGenDijetPtEtaCMBackward->ProjectionY(Form("hGenEtaCMBackward_%d", i), 
                                                hGenDijetPtEtaCMBackward->GetXaxis()->FindBin(ptLow), 
                                                hGenDijetPtEtaCMBackward->GetXaxis()->FindBin(ptHi)-1));
        hRefEtaCMBackward = dynamic_cast<TH1D *>(hRefDijetPtEtaCMBackward->ProjectionY(Form("hRefEtaCMBackward_%d", i), 
                                                hRefDijetPtEtaCMBackward->GetXaxis()->FindBin(ptLow), 
                                                hRefDijetPtEtaCMBackward->GetXaxis()->FindBin(ptHi)-1));
        // hRefSelEtaCMBackward = dynamic_cast<TH1D *>(hRefSelDijetPtEtaCMBackward->ProjectionY(Form("hRefSelEtaCMBackward_%d", i), 
        //                                               hRefSelDijetPtEtaCMBackward->GetXaxis()->FindBin(ptLow), 
        //                                               hRefSelDijetPtEtaCMBackward->GetXaxis()->FindBin(ptHi)-1));

        set1DStyle(hRecoEtaCMForward, 0);
        set1DStyle(hGenEtaCMForward, 2);
        set1DStyle(hRefEtaCMForward, 1);
        // set1DStyle(hRefSelEtaCMForward, 2);

        set1DStyle(hRecoEtaCMBackward, 0);
        set1DStyle(hGenEtaCMBackward, 2);
        set1DStyle(hRefEtaCMBackward, 1);
        // set1DStyle(hRefSelEtaCMBackward, 2);

        rescaleForwardBackward(hRecoEtaCMForward, hRecoEtaCMBackward);
        rescaleForwardBackward(hGenEtaCMForward, hGenEtaCMBackward);
        rescaleForwardBackward(hRefEtaCMForward, hRefEtaCMBackward);
        // rescaleForwardBackward(hRefSelEtaCMForward, hRefSelEtaCMBackward);

        // Compute ratios for forward and backward center of mass frame
        hRecEtaFBCM = dynamic_cast<TH1D *>(hRecoEtaCMForward->Clone(Form("hRecEtaFBCM_%d", i)));
        hRecEtaFBCM->Divide(hRecEtaFBCM, hRecoEtaCMBackward, 1., 1.);
        hGenEtaFBCM = dynamic_cast<TH1D *>(hGenEtaCMForward->Clone(Form("hGenEtaFBCM_%d", i)));
        hGenEtaFBCM->Divide(hGenEtaFBCM, hGenEtaCMBackward, 1., 1.);
        hRefEtaFBCM = dynamic_cast<TH1D *>(hRefEtaCMForward->Clone(Form("hRefEtaFBCM_%d", i)));
        hRefEtaFBCM->Divide(hRefEtaFBCM, hRefEtaCMBackward, 1., 1.);
        // hRefSelEtaFBCM = dynamic_cast<TH1D *>(hRefSelEtaCMForward->Clone(Form("hRefSelEtaFBCM_%d", i)));
        // hRefSelEtaFBCM->Divide(hRefSelEtaFBCM, hRefSelEtaCMBackward, 1., 1.);

        // Plot comparisons
        drawDijetToGenComparison(c, hRecEtaFBCM, hRefEtaFBCM, hGenEtaFBCM, nullptr, nullptr,
                                 (int)hRecoDijetPtEtaCMForward->GetXaxis()->GetBinLowEdge( hRecoDijetPtEtaCMForward->GetXaxis()->FindBin(ptLow) ), 
                                 (int)hRecoDijetPtEtaCMForward->GetXaxis()->GetBinUpEdge( hRecoDijetPtEtaCMForward->GetXaxis()->FindBin(ptHi)-1), 
                                 true, true, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaFBRatio_RecoRefGenComp_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), 
                       ptLow, ptHi));

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
    int ptHatVals[] { 15, 60, 90 };
    int sizeOfPtHatVals = sizeof(ptHatVals)/sizeof(ptHatVals[0]);

    // jet pT and eta binning
    int jetPtVals[] = {30, 50, 90, 120, 200, 500, 1500};
    int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

    double jetEtaVals[] = {-3.6, -2.4, -1.6, 1.6, 2.4, 3.6};
    int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);

    std::cout << Form("Reco histogram to read: ") << Form("hReco%s%sJetPtEtaPtHat", tJetType.Data(), tMatched.Data() ) << std::endl;
    std::cout << Form("Gen histogram to read: ") << Form("hGen%sJetPtEtaPtHat", tJetType.Data() ) << std::endl;
    std::cout << Form("Ref histogram to read: ") << Form("hRef%sJetPtEtaPtHat", tJetType.Data() ) << std::endl;
    std::cout << Form("RefSel histogram to read: ") << Form("hRefSel%sJetPtEtaPtHat", tJetType.Data() ) << std::endl;

    // Retrieve histograms
    TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hReco%s%sJetPtEtaPtHat", tJetType.Data(), tMatched.Data() ) ) );
    // TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRecoMatchedJetPtEtaPtHat"));
    TH3D *hGenPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hGen%sJetPtEtaPtHat", tJetType.Data() ) ) );
    TH3D *hRefPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hRef%sJetPtEtaPtHat", tJetType.Data() ) ) );
    TH3D *hRefSelPtEtaPtHat = dynamic_cast<TH3D *>(f->Get( Form("hRefSel%sJetPtEtaPtHat", tJetType.Data() ) ) );

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


    //
    // Perform analysis and build comparisons for different ptHat intervals.
    // First part loops over jet pT, second part over jet eta.
    //

    // Loop over ptHat bins
    for (int i = 0; i < sizeOfPtHatVals-1; i++) {

        int ptHatBinLow = hRecoPtEtaPtHat->GetZaxis()->FindBin(ptHatVals[i]);
        int ptHatBinHigh = hRecoPtEtaPtHat->GetNbinsZ();
        double ptHatLow = hRecoPtEtaPtHat->GetZaxis()->GetBinLowEdge( hRecoPtEtaPtHat->GetZaxis()->FindBin( ptHatVals[i]));
        double ptHatHigh = hRecoPtEtaPtHat->GetZaxis()->GetBinUpEdge( hRecoPtEtaPtHat->GetNbinsZ() );

        //
        // Set ptHat range
        //
        hRecoPtEtaPtHat->GetZaxis()->SetRange(ptHatLow,  ptHatHigh);
        hGenPtEtaPtHat->GetZaxis()->SetRange(ptHatLow,  ptHatHigh);
        hRefPtEtaPtHat->GetZaxis()->SetRange(ptHatLow,  ptHatHigh);
        hRefSelPtEtaPtHat->GetZaxis()->SetRange(ptHatLow,  ptHatHigh);


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
        for (int j = 0; j < sizeOfJetPtVals-1; j++) {

            //
            // Make projections of the 3D histograms
            //

            int ptBinLow = hRecoPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[j]);
            int ptBinHigh = hRecoPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[j+1]) - 1;
            double ptLow = hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(ptBinLow);
            double ptHigh = hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(ptBinHigh);

            // Reco jets
            hRecoEta = dynamic_cast<TH1D *>(hRecoPtEtaPtHat->ProjectionX(Form("hRecoEta_%d_%d", i, j),
                                                                               ptBinLow, ptBinHigh,
                                                                               ptHatBinLow, ptHatBinHigh));
            set1DStyle(hRecoEta, 0, true);
            checkIntegral(hRecoEta);
            rescaleHisto1D(hRecoEta);

            // Gen jets
            hGenEta = dynamic_cast<TH1D *>(hGenPtEtaPtHat->ProjectionX(Form("hGenEta_%d_%d", i, j),
                                           ptBinLow, ptBinHigh,
                                           ptHatBinLow, ptHatBinHigh));
            set1DStyle(hGenEta, 5, true);
            checkIntegral(hGenEta);
            rescaleHisto1D(hGenEta);


            // Ref jets
            hRefEta = dynamic_cast<TH1D *>(hRefPtEtaPtHat->ProjectionX(Form("hRefEta_%d_%d", i, j),
                                           ptBinLow, ptBinHigh,
                                           ptHatBinLow, ptHatBinHigh));
            set1DStyle(hRefEta, 1, true);
            checkIntegral(hRefEta);
            rescaleHisto1D(hRefEta);
        
            // RefSel jets
            hRefSelEta = dynamic_cast<TH1D *>(hRefSelPtEtaPtHat->ProjectionX(Form("hRefSelEta_%d_%d", i, j),
                                                 ptBinLow, ptBinHigh,
                                                 ptHatBinLow, ptHatBinHigh));
            set1DStyle(hRefSelEta, 2, true);
            checkIntegral(hRefSelEta);
            rescaleHisto1D(hRefSelEta);

            //
            // Ratios of reco, ref and refSel to gen
            //
            hReco2GenEta = dynamic_cast<TH1D *>(hRecoEta->Clone(Form("hReco2GenEta_%d_%d", i, j)));
            // hReco2GenEta->Reset();
            // computeNonBinomialRatio(hRecoEta, hGenEta, hReco2GenEta);
            hReco2GenEta->Divide(hReco2GenEta, hGenEta, 1., 1. /*, "b" */);
            hRef2GenEta = dynamic_cast<TH1D *>(hRefEta->Clone(Form("hRef2GenEta_%d_%d", i, j)));
            // hRef2GenEta->Reset();
            // computeBinomialRatio(hRefEta, hGenEta, hRef2GenEta);
            hRef2GenEta->Divide(hRef2GenEta, hGenEta, 1., 1., "b");
            hRefSel2GenEta = dynamic_cast<TH1D *>(hRefSelEta->Clone(Form("hRefSel2GenEta_%d_%d", i, j)));
            // hRefSel2GenEta->Reset();
            // computeBinomialRatio(hRefSelEta, hGenEta, hRefSel2GenEta);
            hRefSel2GenEta->Divide(hRefSel2GenEta, hGenEta, 1., 1., "b");

            //
            // Plot comparisons
            //
            drawSingleJetToGenComparison(c, hRecoEta, hRefEta, hGenEta, hRefSelEta, nullptr,
                                         hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(ptBinLow), 
                                         hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(ptBinHigh),
                                         ptHatLow, false, false, collisionSystem, collisionEnergy);
            c->SaveAs(Form("%s/%s_%s_jetEta_RecoRefGenComp_ptAve_%d_%d_ptHatLow_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(),
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(ptBinLow),
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(ptBinHigh),
                           (int)ptHatLow));

            //
            // Plot ratios
            //
            drawSingleJetToGenRatio(c, hReco2GenEta, hRef2GenEta, hRefSel2GenEta, nullptr,
                                    hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(ptBinLow), 
                                    hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(ptBinHigh),
                                    ptHatLow, false, false, collisionSystem, collisionEnergy);
            c->SaveAs(Form("%s/%s_%s_jetEta_RecoRef2GenRatio_ptAve_%d_%d_ptHatLow_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), 
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinLowEdge(ptBinLow), 
                           (int)hRecoPtEtaPtHat->GetYaxis()->GetBinUpEdge(ptBinHigh), 
                           (int)ptHatLow));


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

    TH3D *hRecoDataPtEtaPtHat = dynamic_cast<TH3D *>(fData->Get("hRecoInclusiveJetPtEtaPtHat"));
    hRecoDataPtEtaPtHat->SetName("hRecoDataPtEtaPtHat");
    TH3D *hRecoMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hRecoInclusiveJetPtEtaPtHat"));
    hRecoMcPtEtaPtHat->SetName("hRecoMcPtEtaPtHat");
    // TH3D *hRecoMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hRecoMatchedJetPtEtaPtHat"));
    // hRecoMcPtEtaPtHat->SetName("hRecoMcPtEtaPtHat");
    TH3D *hGenMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hGenInclusiveJetPtEtaPtHat"));
    hGenMcPtEtaPtHat->SetName("hGenMcPtEtaPtHat");

    // ptHat binning
    int ptHatVals[] { 15, 60, 90 };
    int sizeOfPtHatVals = sizeof(ptHatVals)/sizeof(ptHatVals[0]);

    // jet pT and eta binning
    int jetPtVals[] = {30, 50, 90, 120, 200};
    int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

    double jetEtaVals[] = {-3.6, -2.4, -1.6, 1.6, 2.4, 3.6};
    int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);


    TH1D *hRecoDataEta{nullptr};
    TH1D *hRecoMcEta{nullptr};
    TH1D *hGenMcEta{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);

    // Jet pT binning
    for ( int i = 0; i < sizeOfJetPtVals-1; i++) {

        //
        // Reco data 2 Reco MC comparison
        //
        // Make projection of the 3D histograms
        hRecoDataEta = dynamic_cast<TH1D *>(hRecoDataPtEtaPtHat->ProjectionX(Form("hRecoDataEta_%d", i),
                                            hRecoDataPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i]), 
                                            hRecoDataPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i+1]),
                                            hRecoDataPtEtaPtHat->GetZaxis()->FindBin(ptHatVals[0]), 
                                            hRecoDataPtEtaPtHat->GetNbinsZ()));
        set1DStyle(hRecoDataEta, 0, kTRUE);

        hRecoMcEta = dynamic_cast<TH1D *>(hRecoMcPtEtaPtHat->ProjectionX(Form("hRecoMcEta_%d", i),
                                                                              hRecoMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i]),
                                                                              hRecoMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i+1]),
                                                                              hRecoMcPtEtaPtHat->GetZaxis()->FindBin(ptHatVals[0]),
                                                                              hRecoMcPtEtaPtHat->GetNbinsZ()));
        set1DStyle(hRecoMcEta, 1, kTRUE);

        // Plot distributions
        drawJecComparison(c, hRecoDataEta, hRecoMcEta,
            hRecoDataPtEtaPtHat->GetYaxis()->GetBinLowEdge(hRecoDataPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i])),
            hRecoDataPtEtaPtHat->GetYaxis()->GetBinUpEdge(hRecoDataPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i+1])),
            ptHatVals[0],
            "Data", "Reco MC", collSystem, energy, false, true,
            hRecoDataPtEtaPtHat->GetXaxis()->GetBinLowEdge(1),
            hRecoDataPtEtaPtHat->GetXaxis()->GetBinUpEdge(hRecoDataPtEtaPtHat->GetXaxis()->GetNbins()));

        c->SaveAs( Form("%s/%s_closureEta_Reco2Reco_pt_%d_%d.pdf", 
            date.Data(), collSystemStr.Data(), 
            int(hRecoDataPtEtaPtHat->GetYaxis()->GetBinLowEdge(hRecoDataPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i]))),
            int(hRecoDataPtEtaPtHat->GetYaxis()->GetBinLowEdge(hRecoDataPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i+1]))) ) );


        //
        // Reco data 2 Gen MC comparison
        //

        c->Clear();
        setPadStyle();

        // Make projection of the 3D histograms
        hGenMcEta = dynamic_cast<TH1D *>(hGenMcPtEtaPtHat->ProjectionX(Form("hGenMcEta_%d", i),
            hGenMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i]),
            hGenMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i+1]),
            hGenMcPtEtaPtHat->GetZaxis()->FindBin( hGenMcPtEtaPtHat->GetZaxis()->FindBin(ptHatVals[0])),
            hGenMcPtEtaPtHat->GetNbinsZ()));

        set1DStyle(hGenMcEta, 5, kTRUE);

        // Plot distributions
        drawJecComparison(c, hRecoDataEta, hGenMcEta,
                          hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge( hGenMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i])),
                          hGenMcPtEtaPtHat->GetYaxis()->GetBinUpEdge( hGenMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i])),
                          hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge( hGenMcPtEtaPtHat->GetZaxis()->FindBin(ptHatVals[0])),
                          "Data", "Gen MC", collSystem, energy, false, true,
                          hGenMcPtEtaPtHat->GetXaxis()->GetBinLowEdge(1),
                          hGenMcPtEtaPtHat->GetXaxis()->GetBinUpEdge(hGenMcPtEtaPtHat->GetXaxis()->GetNbins())
                        );

        // Save canvas
        c->SaveAs( Form("%s/%s_closureEta_Reco2Gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), 
                    int(hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge(hGenMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i]))),
                    int(hGenMcPtEtaPtHat->GetYaxis()->GetBinLowEdge(hGenMcPtEtaPtHat->GetYaxis()->FindBin(jetPtVals[i+1]))) ) );
    } // for (unsigned int i = 0; i < jetPtBinsLow.size(); i++)

    delete c; c = nullptr;
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
    double dijetPtVals[] {  50., 90., 120., 180., 200., 250., 300., 350., 500.};
    int sizeOfDijetPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    //
    // Declare histograms
    //

    // Data
    TH3D *hRecoDataDijetPtEtaPhiLab = dynamic_cast<TH3D *>(fData->Get("hRecoDijetPtEtaPhiWeighted"));
    if ( !hRecoDataDijetPtEtaPhiLab ) {
        std::cerr << "Data histogram hRecoDijetPtEtaPhiWeighted not found in file." << std::endl; return;
    }
    hRecoDataDijetPtEtaPhiLab->SetName("hRecoDataDijetPtEtaPhiLab");
    TH3D *hRecoDataDijetPtEtaPhiCM = dynamic_cast<TH3D *>(fData->Get("hRecoDijetPtEtaPhiCMWeighted"));
    if ( !hRecoDataDijetPtEtaPhiCM ) {
        std::cerr << "Data histogram hRecoDijetPtEtaPhiCMWeighted not found in file." << std::endl; return;
    }
    hRecoDataDijetPtEtaPhiCM->SetName("hRecoDataDijetPtEtaPhiCM");
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

    rescaleHisto3D(hRecoDataDijetPtEtaPhiLab);
    rescaleHisto3D(hRecoDataDijetPtEtaPhiCM);
    rescaleHisto2D(hRecoDataDijetPtEtaCMForward);
    rescaleHisto2D(hRecoDataDijetPtEtaCMBackward);

    TH1D *hRecoDataDijetEtaLab{nullptr};
    TH1D *hRecoDataDijetEtaCM{nullptr};
    TH1D *hRecoDataDijetEtaForward{nullptr};
    TH1D *hRecoDataDijetEtaBackward{nullptr};
    TH1D *hRecoDataDijetFBRatio{nullptr};

    // MC

    // Reco
    TH3D *hRecoMcDijetPtEtaPhiLab = dynamic_cast<TH3D *>(fMc->Get("hRecoDijetPtEtaPhiWeighted"));
    if ( !hRecoMcDijetPtEtaPhiLab ) {
        std::cerr << "MC histogram hRecoDijetPtEtaPhiWeighted not found in file." << std::endl; return;
    }
    hRecoMcDijetPtEtaPhiLab->SetName("hRecoMcDijetPtEtaPhiLab");
    TH3D *hRecoMcDijetPtEtaPhiCM = dynamic_cast<TH3D *>(fMc->Get("hRecoDijetPtEtaPhiCMWeighted"));
    if ( !hRecoMcDijetPtEtaPhiCM ) {
        std::cerr << "MC histogram hRecoDijetPtEtaPhiCMWeighted not found in file." << std::endl; return;
    }
    hRecoMcDijetPtEtaPhiCM->SetName("hRecoMcDijetPtEtaPhiCM");
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

    rescaleHisto3D(hRecoMcDijetPtEtaPhiLab);
    rescaleHisto3D(hRecoMcDijetPtEtaPhiCM);
    rescaleHisto2D(hRecoMcDijetPtEtaCMForward);
    rescaleHisto2D(hRecoMcDijetPtEtaCMBackward);

    TH1D *hRecoMcDijetEtaLab{nullptr};
    TH1D *hRecoMcDijetEtaCM{nullptr};
    TH1D *hRecoMcDijetEtaForward{nullptr};
    TH1D *hRecoMcDijetEtaBackward{nullptr};
    TH1D *hRecoMcDijetFBRatio{nullptr};

    // Gen
    TH3D *hGenMcDijetPtEtaPhiLab = dynamic_cast<TH3D *>(fMc->Get("hGenDijetPtEtaPhiWeighted"));
    if ( !hGenMcDijetPtEtaPhiLab ) {
        std::cerr << "MC histogram hGenDijetPtEtaPhiWeighted not found in file." << std::endl; return;
    }
    hGenMcDijetPtEtaPhiLab->SetName("hGenMcDijetPtEtaPhiLab");
    TH3D *hGenMcDijetPtEtaPhiCM = dynamic_cast<TH3D *>(fMc->Get("hGenDijetPtEtaPhiCMWeighted"));
    if ( !hGenMcDijetPtEtaPhiCM ) {
        std::cerr << "MC histogram hGenDijetPtEtaPhiCMWeighted not found in file." << std::endl; return;
    }
    hGenMcDijetPtEtaPhiCM->SetName("hGenMcDijetPtEtaPhiCM");
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

    rescaleHisto3D(hGenMcDijetPtEtaPhiLab);
    rescaleHisto3D(hGenMcDijetPtEtaPhiCM);
    rescaleHisto2D(hGenMcDijetPtEtaCMForward);
    rescaleHisto2D(hGenMcDijetPtEtaCMBackward);

    TH1D *hGenMcDijetEtaLab{nullptr};
    TH1D *hGenMcDijetEtaCM{nullptr};
    TH1D *hGenMcDijetEtaForward{nullptr};
    TH1D *hGenMcDijetEtaBackward{nullptr};
    TH1D *hGenMcDijetFBRatio{nullptr};

    // Ref
    TH3D *hRefMcDijetPtEtaPhiLab = dynamic_cast<TH3D *>(fMc->Get("hRefDijetPtEtaPhiWeighted"));
    if ( !hRefMcDijetPtEtaPhiLab ) {
        std::cerr << "MC histogram hRefDijetPtEtaPhiWeighted not found in file." << std::endl; return;
    }
    hRefMcDijetPtEtaPhiLab->SetName("hRefMcDijetPtEtaPhiLab");
    TH3D *hRefMcDijetPtEtaPhiCM = dynamic_cast<TH3D *>(fMc->Get("hRefDijetPtEtaPhiCMWeighted"));
    if ( !hRefMcDijetPtEtaPhiCM ) {
        std::cerr << "MC histogram hRefDijetPtEtaPhiCMWeighted not found in file." << std::endl; return;
    }
    hRefMcDijetPtEtaPhiCM->SetName("hRefMcDijetPtEtaPhiCM");
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

    rescaleHisto3D(hRefMcDijetPtEtaPhiLab);
    rescaleHisto3D(hRefMcDijetPtEtaPhiCM);
    rescaleHisto2D(hRefMcDijetPtEtaCMForward);
    rescaleHisto2D(hRefMcDijetPtEtaCMBackward);

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
        int ptBinLow = hRecoDataDijetPtEtaPhiLab->GetYaxis()->FindBin(ptLow);
        int ptBinHigh = hRecoDataDijetPtEtaPhiLab->GetYaxis()->FindBin(ptHi) - 1;


        //
        // Data
        //

        hRecoDataDijetEtaLab = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaPhiLab->ProjectionY( Form("hRecoDataDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoDataDijetEtaLab ) {
            std::cerr << Form("Data histogram not found: hRecoDataDijetEtaLab_%d", i) << std::endl; return;
        }
        hRecoDataDijetEtaLab->SetName( Form("hRecoDataDijetEtaLab_%d", i) );
        set1DStyle( hRecoDataDijetEtaLab, 2 );
        rescaleHisto1D( hRecoDataDijetEtaLab );

        hRecoDataDijetEtaCM = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaPhiCM->ProjectionY( Form("hRecoDataDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
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

        // Ratios
        hRecoDataDijetFBRatio = dynamic_cast<TH1D*>( hRecoDataDijetEtaForward->Clone( Form("hRecoDataDijetFBRatio_%d", i) ) );
        hRecoDataDijetFBRatio->SetName( Form("hRecoDataDijetFBRatio_%d", i) );
        hRecoDataDijetFBRatio->Divide( hRecoDataDijetFBRatio, hRecoDataDijetEtaBackward, 1., 1. );

        //
        // MC (reco)
        //

        // Eta full (lab frame)
        hRecoMcDijetEtaLab = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaPhiLab->ProjectionY( Form("hRecoMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoMcDijetEtaLab ) {
            std::cerr << Form("MC histogram not found: hRecoDijetEta1DLab_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaLab->SetName( Form("hRecoMcDijetEtaLab_%d", i) );
        set1DStyle( hRecoMcDijetEtaLab, 0 );
        rescaleHisto1D( hRecoMcDijetEtaLab );

        // Eta full (CM frame)
        hRecoMcDijetEtaCM = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaPhiCM->ProjectionY( Form("hRecoMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
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
        // integralForward = hRecoMcDijetEtaForward->Integral();

        // Eta backward (CM frame)
        hRecoMcDijetEtaBackward = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaCMBackward->ProjectionY( Form("hRecoMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRecoMcDijetEtaBackward ) {
            std::cerr << Form("MC histogram not found: hRecoMcDijetEtaBackward_%d", i) << std::endl; return;
        }
        hRecoMcDijetEtaBackward->SetName( Form("hRecoMcDijetEtaBackward_%d", i) );
        set1DStyle( hRecoMcDijetEtaBackward, 0 );
        // double integralBackward = hRecoMcDijetEtaBackward->Integral();

        // Normalize the forward and backward distributions
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
        hGenMcDijetEtaLab = dynamic_cast<TH1D*>( hGenMcDijetPtEtaPhiLab->ProjectionY( Form("hGenMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hGenMcDijetEtaLab ) {
            std::cerr << Form("MC histogram not found: hGenMcDijetEtaLab_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaLab->SetName( Form("hGenMcDijetEtaLab_%d", i) );
        set1DStyle( hGenMcDijetEtaLab, 5 );
        rescaleHisto1D( hGenMcDijetEtaLab );    

        // Eta full (CM frame)
        hGenMcDijetEtaCM = dynamic_cast<TH1D*>( hGenMcDijetPtEtaPhiCM->ProjectionY( Form("hGenMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
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
        // double integralForward = hGenMcDijetEtaForward->Integral();

        // Eta backward (CM frame)
        hGenMcDijetEtaBackward = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCMBackward->ProjectionY( Form("hGenMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hGenMcDijetEtaBackward ) {
            std::cerr << Form("MC histogram not found: hGenMcDijetEtaBackward_%d", i) << std::endl; return;
        }
        hGenMcDijetEtaBackward->SetName( Form("hGenMcDijetEtaBackward_%d", i) );
        set1DStyle( hGenMcDijetEtaBackward, 5 );
        // double integralBackward = hGenMcDijetEtaBackward->Integral();
        
        // Normalize the forward and backward distributions
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
        hRefMcDijetEtaLab = dynamic_cast<TH1D*>( hRefMcDijetPtEtaPhiLab->ProjectionY( Form("hRefMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRefMcDijetEtaLab ) {
            std::cerr << Form("MC histogram not found: hRefMcDijetEtaLab_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaLab->SetName( Form("hRefMcDijetEtaLab_%d", i) );
        set1DStyle( hRefMcDijetEtaLab, 1 );
        rescaleHisto1D( hRefMcDijetEtaLab );

        // Eta full (CM frame)
        hRefMcDijetEtaCM = dynamic_cast<TH1D*>( hRefMcDijetPtEtaPhiLab->ProjectionY( Form("hRefMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
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
        // integralForward = hRefMcDijetEtaForward->Integral();

        // Eta backward (CM frame)
        hRefMcDijetEtaBackward = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCMBackward->ProjectionY( Form("hRefMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
        if ( !hRefMcDijetEtaBackward ) {
            std::cerr << Form("MC histogram not found: hRefMcDijetEtaBackward_%d", i) << std::endl; return;
        }
        hRefMcDijetEtaBackward->SetName( Form("hRefMcDijetEtaBackward_%d", i) );
        set1DStyle( hRefMcDijetEtaBackward, 1 );
        // integralBackward = hRefMcDijetEtaBackward->Integral();

        // Normalize the forward and backward distributions
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

        c->cd(); c->Clear(); setPadStyle();
        setPadStyle();
        hRecoDataDijetEtaLab->Draw();
        hRecoDataDijetEtaLab->GetXaxis()->SetRangeUser(-2.5, 2.5);
        hRecoDataDijetEtaLab->GetYaxis()->SetRangeUser(0.00001, 0.45);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
        c->SaveAs( Form("%s/%s_%s_dijetEtaLab_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

        // CM frame
        drawDijetToGenComparison(c, hRecoMcDijetEtaCM, hGenMcDijetEtaCM, hRefMcDijetEtaCM, nullptr, hRecoDataDijetEtaCM,
                            (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

        c->cd(); c->Clear(); setPadStyle();
        hRecoDataDijetEtaCM->Draw();
        hRecoDataDijetEtaCM->GetXaxis()->SetRangeUser(-2., 2.);
        hRecoDataDijetEtaCM->GetYaxis()->SetRangeUser(0.00001, 0.45);
        t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
        c->SaveAs( Form("%s/%s_%s_dijetEtaCM_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

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
    }

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
            rescaleHisto1D(hRecoEtaEmb[iPtHat][iJetPt]);
            set1DStyle(hRecoEtaEmb[iPtHat][iJetPt], 0);
            hGenEtaEmb[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hGenPtEtaPtHatEmb->ProjectionX(Form("hGenEtaEmb_%zu_%zu", iPtHat, iJetPt), 
                                                                                             jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleHisto1D(hGenEtaEmb[iPtHat][iJetPt]);
            set1DStyle(hGenEtaEmb[iPtHat][iJetPt], 0);

            // Pythia
            hRecoEtaPythia[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hRecoPtEtaPtHatPythia->ProjectionX(Form("hRecoEtaPythia_%zu_%zu", iPtHat, iJetPt), 
                                                                                                   jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleHisto1D(hRecoEtaPythia[iPtHat][iJetPt]);
            set1DStyle(hRecoEtaPythia[iPtHat][iJetPt], 1);
            hGenEtaPythia[iPtHat][iJetPt] = dynamic_cast<TH1D *>(hGenPtEtaPtHatPythia->ProjectionX(Form("hGenEtaPythia_%zu_%zu", iPtHat, iJetPt), 
                                                                                                 jetPtBinsLow[iJetPt], jetPtBinsHigh[iJetPt], ptHatBins[iPtHat], ptHatBinsMax));
            rescaleHisto1D(hGenEtaPythia[iPtHat][iJetPt]);
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
    int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
                             100, 110,  120, 130, 140,
                             150, 160,  180, 200, 250, 
                             300, 500 };
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);

    // Dijet dijstributions for embedding
    TH3D *hRecoDijetPtEtaPhiCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRecoDijetPtEtaPhiCMWeighted"));
    if (!hRecoDijetPtEtaPhiCMEmb) {
        std::cerr << "Error: hRecoDijetPtEtaPhiCMWeighted for embedding not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaPhiCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hGenDijetPtEtaPhiCMWeighted"));
    if (!hGenDijetPtEtaPhiCMEmb) {
        std::cerr << "Error: hGenDijetPtEtaPhiCMWeighted for embedding not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaPhiCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRefDijetPtEtaPhiCMWeighted"));
    if (!hRefDijetPtEtaPhiCMEmb) {
        std::cerr << "Error: hRefDijetPtEtaPhiCMWeighted for embedding not found in file." << std::endl;
        return;
    }

    // Dijet dijstributions for pythia
    TH3D *hRecoDijetPtEtaPhiCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hRecoDijetPtEtaPhiCMWeighted"));
    if (!hRecoDijetPtEtaPhiCMPythia) {
        std::cerr << "Error: hRecoDijetPtEtaPhiCMWeighted for pythia not found in file." << std::endl;
        return;
    }
    TH3D *hGenDijetPtEtaPhiCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hGenDijetPtEtaPhiCMWeighted"));
    if (!hGenDijetPtEtaPhiCMPythia) {
        std::cerr << "Error: hGenDijetPtEtaPhiCMWeighted for pythia not found in file." << std::endl;
        return;
    }
    TH3D *hRefDijetPtEtaPhiCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hRefDijetPtEtaPhiCMWeighted"));
    if (!hRefDijetPtEtaPhiCMPythia) {
        std::cerr << "Error: hRefDijetPtEtaPhiCMWeighted for pythia not found in file." << std::endl;
        return;
    }

    // Rescaling of 3D histograms
    rescaleHisto3D(hRecoDijetPtEtaPhiCMEmb);
    rescaleHisto3D(hGenDijetPtEtaPhiCMEmb);
    rescaleHisto3D(hRefDijetPtEtaPhiCMEmb);
    rescaleHisto3D(hRecoDijetPtEtaPhiCMPythia);
    rescaleHisto3D(hGenDijetPtEtaPhiCMPythia);
    rescaleHisto3D(hRefDijetPtEtaPhiCMPythia);

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
        int ptLow = dijetPtNewVals[i];
        int ptHi = dijetPtNewVals[i + 1];

        // Create 1D histograms for embedding
        hRecoDijetEtaCMEmb = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiCMEmb->ProjectionY(Form("hRecoDijetEtaCMEmb_%d", i), 
                                                  hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow), 
                                                  hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1));
        if (!hRecoDijetEtaCMEmb) {
            std::cerr << "Error: hRecoDijetEtaCMEmb_" << i << " not found." << std::endl;
            return;
        }
        hGenDijetEtaCMEmb = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiCMEmb->ProjectionY(Form("hGenDijetEtaCMEmb_%d", i), 
                                                  hGenDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow), 
                                                  hGenDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1));
        if (!hGenDijetEtaCMEmb) {
            std::cerr << "Error: hGenDijetEtaCMEmb_" << i << " not found." << std::endl;
            return;
        }
        hRefDijetEtaCMEmb = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiCMEmb->ProjectionY(Form("hRefDijetEtaCMEmb_%d", i),
                                                  hRefDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow), 
                                                  hRefDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1));
        if (!hRefDijetEtaCMEmb) {
            std::cerr << "Error: hRefDijetEtaCMEmb_" << i << " not found." << std::endl;
            return;
        }

        set1DStyle(hRecoDijetEtaCMEmb, 0, true);
        set1DStyle(hGenDijetEtaCMEmb, 0, true);
        set1DStyle(hRefDijetEtaCMEmb, 0, true);

        // Create 1D histograms for pythia
        hRecoDijetEtaCMPythia = dynamic_cast<TH1D *>(hRecoDijetPtEtaPhiCMPythia->ProjectionY(Form("hRecoDijetEtaCMPythia_%d", i), 
                                                  hRecoDijetPtEtaPhiCMPythia->GetXaxis()->FindBin(ptLow), 
                                                  hRecoDijetPtEtaPhiCMPythia->GetXaxis()->FindBin(ptHi)-1));
        if (!hRecoDijetEtaCMPythia) {
            std::cerr << "Error: hRecoDijetEtaCMPythia_" << i << " not found." << std::endl;
            return;
        }
        hGenDijetEtaCMPythia = dynamic_cast<TH1D *>(hGenDijetPtEtaPhiCMPythia->ProjectionY(Form("hGenDijetEtaCMPythia_%d", i), 
                                                  hGenDijetPtEtaPhiCMPythia->GetXaxis()->FindBin(ptLow), 
                                                  hGenDijetPtEtaPhiCMPythia->GetXaxis()->FindBin(ptHi)-1));
        if (!hGenDijetEtaCMPythia) {
            std::cerr << "Error: hGenDijetEtaCMPythia_" << i << " not found." << std::endl;
            return;
        }
        hRefDijetEtaCMPythia = dynamic_cast<TH1D *>(hRefDijetPtEtaPhiCMPythia->ProjectionY(Form("hRefDijetEtaCMPythia_%d", i), 
                                                  hRefDijetPtEtaPhiCMPythia->GetXaxis()->FindBin(ptLow), 
                                                  hRefDijetPtEtaPhiCMPythia->GetXaxis()->FindBin(ptHi)-1));
        if (!hRefDijetEtaCMPythia) {
            std::cerr << "Error: hRefDijetEtaCMPythia_" << i << " not found." << std::endl;
            return;
        }

        set1DStyle(hRecoDijetEtaCMPythia, 1, true);
        set1DStyle(hGenDijetEtaCMPythia, 1, true);
        set1DStyle(hRefDijetEtaCMPythia, 1, true);

        // Draw and save the comparisons
        c->Clear();
        setPadStyle();
        drawDijetToGenComparison(c, hRecoDijetEtaCMEmb, hRecoDijetEtaCMPythia, nullptr, nullptr, nullptr,
                                 (int)hRecoDijetPtEtaPhiCMEmb->GetXaxis()->GetBinLowEdge(hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow)), 
                                 (int)hRecoDijetPtEtaPhiCMEmb->GetXaxis()->GetBinUpEdge(hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1), 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_recoComp_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));
    
        c->Clear();
        setPadStyle();
        drawDijetToGenComparison(c, hGenDijetEtaCMEmb, hGenDijetEtaCMPythia, nullptr, nullptr, nullptr,
                                 (int)hGenDijetPtEtaPhiCMEmb->GetXaxis()->GetBinLowEdge(hGenDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow)), 
                                 (int)hGenDijetPtEtaPhiCMEmb->GetXaxis()->GetBinUpEdge(hGenDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1), 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_genComp_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        c->Clear();
        setPadStyle();
        drawDijetToGenComparison(c, hRefDijetEtaCMEmb, hRefDijetEtaCMPythia, nullptr, nullptr, nullptr,
                                 (int)hRefDijetPtEtaPhiCMEmb->GetXaxis()->GetBinLowEdge(hRefDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow)), 
                                 (int)hRefDijetPtEtaPhiCMEmb->GetXaxis()->GetBinUpEdge(hRefDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1), 
                                 true, false, collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_refComp_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

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
        c->Clear();
        setPadStyle();
        hRecoDijetEtaCM_emb2pythia->Draw();
        hGenDijetEtaCM_emb2pythia->Draw("same");
        hRefDijetEtaCM_emb2pythia->Draw("same");
        hRecoDijetEtaCM_emb2pythia->GetXaxis()->SetRangeUser(-2.4, 2.4);
        hRecoDijetEtaCM_emb2pythia->GetYaxis()->SetRangeUser(0.8, 1.2);
        hRecoDijetEtaCM_emb2pythia->GetYaxis()->SetTitle("Embedding/Pythia");
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%d < p_{T}^{ave} < %d GeV", 
                       (int)hRecoDijetPtEtaPhiCMEmb->GetXaxis()->GetBinLowEdge(hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow)), 
                       (int)hRecoDijetPtEtaPhiCMEmb->GetXaxis()->GetBinUpEdge(hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1)) );
        leg = new TLegend(0.35, 0.25, 0.65, 0.4);
        leg->SetTextSize(0.04);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hRecoDijetEtaCM_emb2pythia, "Reco", "p");
        leg->AddEntry(hGenDijetEtaCM_emb2pythia, "Gen", "p");
        leg->AddEntry(hRefDijetEtaCM_emb2pythia, "Ref", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_ratio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(),
                       (int)hRecoDijetPtEtaPhiCMEmb->GetXaxis()->GetBinLowEdge(hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptLow)), 
                       (int)hRecoDijetPtEtaPhiCMEmb->GetXaxis()->GetBinUpEdge(hRecoDijetPtEtaPhiCMEmb->GetXaxis()->FindBin(ptHi)-1)));
    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
}

//________________
void recoDijetDistributions(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129", bool isMc = false) {
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

    // Reco dijet info [9 dimensions]
    // 0 - dijet pt ave, 1 - dijet eta, 2 - dijet phi,
    // 3 - lead pt, 4 - lead eta, 5 - lead phi,
    // 6 - sublead pt, 7 - sublead eta, 8 - sublead phi
    THnSparseD *hRecoDijetInfo = dynamic_cast<THnSparseD *>(f->Get("hRecoDijetInfo"));
    if (!hRecoDijetInfo) {
        std::cerr << "Error: hRecoDijetInfo not found in file." << std::endl;
        return;
    }

    // Plot dijet ptAve distribution
    TH1D *hDijetPtAve = dynamic_cast<TH1D *>(hRecoDijetInfo->Projection(0));
    set1DStyle(hDijetPtAve, 0, true);

    // Next distributions will be plotted for the given dijet ptAve binning
    TH1D *hDijetEta{nullptr};
    TH2D *hLeadPtVsSubLeadPt{nullptr};
    TH2D *hLeadEtaVsSubLeadEta{nullptr};

    TH3D *hRecoDijetPtEtaPhiLab{nullptr};
    TH3D *hRecoDijetPtEtaPhiLabMatched{nullptr};
    TH3D *hRecoDijetPtEtaPhiCM{nullptr};
    TH3D *hRecoDijetPtEtaPhiCMMatched{nullptr};
    TH2D *hDijetPtAveVsEtaLabInclusive{nullptr};
    TH2D *hDijetPtAveVsEtaLabMatched{nullptr};
    TH2D *hDijetPtAveVsEtaCMInclusive{nullptr};
    TH2D *hDijetPtAveVsEtaCMMatched{nullptr};

    if ( isMc ) {
        hRecoDijetPtEtaPhiLab = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiWeighted"));
        if (!hRecoDijetPtEtaPhiLab) {
            std::cerr << "Error: hRecoDijetPtEtaPhiWeighted not found in file." << std::endl;
            return;
        }
        hRecoDijetPtEtaPhiLab->SetName("hRecoDijetPtEtaPhiLab");

        hRecoDijetPtEtaPhiLabMatched = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiMatched"));
        if (!hRecoDijetPtEtaPhiLabMatched) {
            std::cerr << "Error: hRecoDijetPtEtaPhiMatched not found in file." << std::endl;
            return;
        }
        hRecoDijetPtEtaPhiLabMatched->SetName("hRecoDijetPtEtaPhiLabMatched");

        hRecoDijetPtEtaPhiCM = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiCMWeighted"));
        if (!hRecoDijetPtEtaPhiCM) {
            std::cerr << "Error: hRecoDijetPtEtaPhiCMWeighted not found in file." << std::endl;
            return;
        }
        hRecoDijetPtEtaPhiCM->SetName("hRecoDijetPtEtaPhiCM");

        hRecoDijetPtEtaPhiCMMatched = dynamic_cast<TH3D *>(f->Get("hRecoDijetPtEtaPhiCMMatched"));
        if (!hRecoDijetPtEtaPhiCMMatched) {
            std::cerr << "Error: hRecoDijetPtEtaPhiCMMatched not found in file." << std::endl;
            return;
        }
        hRecoDijetPtEtaPhiCMMatched->SetName("hRecoDijetPtEtaPhiCMMatched");

        hDijetPtAveVsEtaLabInclusive = dynamic_cast<TH2D *>( hRecoDijetPtEtaPhiLab->Project3D("xy") );
        hDijetPtAveVsEtaLabInclusive->SetName("hDijetPtAveVsEtaLabInclusive");
        set2DStyle(hDijetPtAveVsEtaLabInclusive);

        hDijetPtAveVsEtaLabMatched = dynamic_cast<TH2D *>( hRecoDijetPtEtaPhiLabMatched->Project3D("xy") );
        hDijetPtAveVsEtaLabMatched->SetName("hDijetPtAveVsEtaLabMatched");
        set2DStyle(hDijetPtAveVsEtaLabMatched);

        hDijetPtAveVsEtaCMInclusive = dynamic_cast<TH2D *>( hRecoDijetPtEtaPhiCM->Project3D("xy") );
        hDijetPtAveVsEtaCMInclusive->SetName("hDijetPtAveVsEtaCMInclusive");
        set2DStyle(hDijetPtAveVsEtaCMInclusive);

        hDijetPtAveVsEtaCMMatched = dynamic_cast<TH2D *>( hRecoDijetPtEtaPhiCMMatched->Project3D("xy") );
        hDijetPtAveVsEtaCMMatched->SetName("hDijetPtAveVsEtaCMMatched");
        set2DStyle(hDijetPtAveVsEtaCMMatched);
    } // if ( isMc )

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    TLegend *leg{nullptr};
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    hDijetPtAve->Draw();
    hDijetPtAve->GetYaxis()->SetTitle("1/N dN/dp_{T}^{ave}");
    hDijetPtAve->GetXaxis()->SetRangeUser(40, hDijetPtAve->GetXaxis()->GetBinUpEdge(hDijetPtAve->GetXaxis()->GetNbins()));
    c->SetLogy(1);
    c->SetLogx(1);
    c->SetGrid(1, 1);
    c->SaveAs(Form("%s/%s_%s_recoDijet_ptAve.pdf", date.Data(), collSystemStr.Data(), directionStr.Data()));
    c->SetLogy(0);
    c->SetLogx(0);
    c->SetGrid(0, 0);

    if ( isMc ) {
        c->Clear();
        setPadStyle();
        hDijetPtAveVsEtaLabInclusive->Draw("colz");
        hDijetPtAveVsEtaLabInclusive->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetPtAveVsEtaLabInclusive->GetYaxis()->SetRangeUser(35, 505);
        t.DrawLatexNDC(0.3, 0.8, "Inclusive dijets (lab)");
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_ptAveVsEtaLabInclusive.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data()));

        c->Clear();
        setPadStyle();
        hDijetPtAveVsEtaLabMatched->Draw("colz");
        hDijetPtAveVsEtaLabMatched->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetPtAveVsEtaLabMatched->GetYaxis()->SetRangeUser(35, 505);
        t.DrawLatexNDC(0.3, 0.8, "Matched dijets (lab)");
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_ptAveVsEtaLabMatched.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data()));

        c->Clear();
        setPadStyle();
        hDijetPtAveVsEtaCMInclusive->Draw("colz");
        hDijetPtAveVsEtaCMInclusive->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetPtAveVsEtaCMInclusive->GetYaxis()->SetRangeUser(35, 505);
        t.DrawLatexNDC  (0.3, 0.8, "Inclusive dijets (CM)");
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_ptAveVsEtaCMInclusive.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data()));

        c->Clear();
        setPadStyle();
        hDijetPtAveVsEtaCMMatched->Draw("colz");
        hDijetPtAveVsEtaCMMatched->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetPtAveVsEtaCMMatched->GetYaxis()->SetRangeUser(35, 505);
        t.DrawLatexNDC(0.3, 0.8, "Matched dijets (CM)");
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_ptAveVsEtaCMMatched.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data()));

        
        c->Clear();
        setPadStyle();
        hDijetPtAveVsEtaLabMatched->Divide(hDijetPtAveVsEtaLabMatched, hDijetPtAveVsEtaLabInclusive, 1., 1., "b");
        hDijetPtAveVsEtaLabMatched->Draw("colz");
        hDijetPtAveVsEtaLabMatched->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetPtAveVsEtaLabMatched->GetYaxis()->SetRangeUser(35, 505);
        hDijetPtAveVsEtaLabMatched->GetZaxis()->SetRangeUser(0.01, 1.05);
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_ptAveVsEtaLabMatched_fraction.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data()));

        c->Clear();
        setPadStyle();
        hDijetPtAveVsEtaCMMatched->Divide(hDijetPtAveVsEtaCMMatched, hDijetPtAveVsEtaCMInclusive, 1., 1., "b");
        hDijetPtAveVsEtaCMMatched->Draw("colz");
        hDijetPtAveVsEtaCMMatched->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetPtAveVsEtaCMMatched->GetYaxis()->SetRangeUser(35, 505);
        hDijetPtAveVsEtaCMMatched->GetZaxis()->SetRangeUser(0.01, 1.05);
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_ptAveVsEtaCMMatched_fraction.pdf", date.Data(),
                       collSystemStr.Data(), directionStr.Data()));
    }

    // Loop over dijet ptAve bins
    for (int i = 0; i < sizeOfPtVals - 1; ++i) {
        int ptLow = dijetPtNewVals[i];
        int ptHi = dijetPtNewVals[i + 1];
        int ptBinLow = hRecoDijetInfo->GetAxis(0)->FindBin(ptLow);
        int ptBinHi = hRecoDijetInfo->GetAxis(0)->FindBin(ptHi) - 1;

        // Set dijet ptAve range
        hRecoDijetInfo->GetAxis(0)->SetRange(ptBinLow, ptBinHi);

        // Dijet eta distribution
        c->Clear();
        setPadStyle();
        hDijetEta = dynamic_cast<TH1D *>( hRecoDijetInfo->Projection(1) );
        hDijetEta->SetName(Form("hDijetEta_%d_%d", ptLow, ptHi));
        set1DStyle(hDijetEta, 0, true);
        hDijetEta->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hDijetEta->GetYaxis()->SetRangeUser(0.000001, 0.12);
        hDijetEta->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        hDijetEta->Draw();
        t.DrawLatexNDC(0.35, 0.8, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_eta_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        // Dijet lead vs sublead pt distribution
        c->Clear();
        setPadStyle();
        hLeadPtVsSubLeadPt = dynamic_cast<TH2D *>( hRecoDijetInfo->Projection(3, 6) );
        hLeadPtVsSubLeadPt->SetName(Form("hLeadPtVsSubLeadPt_%d_%d", ptLow, ptHi));
        set2DStyle(hLeadPtVsSubLeadPt);
        hLeadPtVsSubLeadPt->GetXaxis()->SetRangeUser(30, 530);
        hLeadPtVsSubLeadPt->GetYaxis()->SetRangeUser(30, 530);
        //hLeadPtVsSubLeadPt->SetRange(0., 0.12);
        hLeadPtVsSubLeadPt->Draw("colz");
        t.DrawLatexNDC(0.35, 0.8, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_leadVsSubLeadPt_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));

        // Dijet lead vs sublead eta distribution
        c->Clear();
        setPadStyle();
        hLeadEtaVsSubLeadEta = dynamic_cast<TH2D *>( hRecoDijetInfo->Projection(4, 7) );
        hLeadEtaVsSubLeadEta->SetName(Form("hLeadEtaVsSubLeadEta_%d_%d", ptLow, ptHi));
        set2DStyle(hLeadEtaVsSubLeadEta);
        hLeadEtaVsSubLeadEta->GetXaxis()->SetRangeUser(-3.0, 3.0);
        hLeadEtaVsSubLeadEta->GetYaxis()->SetRangeUser(-3.0, 3.0);
        //hLeadEtaVsSubLeadEta->SetRange(0., 0.12);
        hLeadEtaVsSubLeadEta->Draw("colz");
        t.DrawLatexNDC(0.35, 0.8, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_recoDijet_leadVsSubLeadEta_ptAve_%d_%d.pdf", date.Data(), 
                       collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));
    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)

    c->Clear();
    delete c;
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
    int dataTrigger = 0;               // 0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100
    TString dataStr = (dataTrigger == 0) ? "MB" : ((dataTrigger == 1) ? "Jet60" : ((dataTrigger == 2) ? "Jet80" : ((dataTrigger == 3) ? "Jet100" : "unknownData")));
    TString dataDirectionStr = (direction == 0) ? "Pbgoing" : ((direction == 1) ? "pgoing" : "");
    int jetType = 0; // 0 - inclusive, 1 - lead, 2 - sublead
    int matchType = 0; // 0 - inclusive, 1 - matched, 2 - unmatched
    int cutType = 2; // 0 - no cuts, 1 - trkMax, 2 - jetId
    TString cutTypeStr = (cutType == 0) ? "noCut_" : ((cutType == 1) ? "" : ((cutType == 2) ? "jetId_" : "unknownCut"));

    //
    // Embedding
    //
    TFile *pPb8160EmbedFile = nullptr;
    if ( direction < 2 ) {
        pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_%seta19.root", uname.Data(), directionStr.Data(), directionStr.Data(), cutTypeStr.Data()) );
        if ( !pPb8160EmbedFile ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_%seta19.root", uname.Data(), directionStr.Data(), directionStr.Data(), cutTypeStr.Data()) << std::endl;
            return;
        }
    }
    else {
        pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_%seta19.root", uname.Data(), cutTypeStr.Data()) );
        if ( !pPb8160EmbedFile ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_%seta19.root", uname.Data(), cutTypeStr.Data()) << std::endl;
            return;
        }
    }
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/macro/pPb8160_Pbgoing_noTrkMax/oEmbedding_pPb8160_Pbgoing_noTrkMax.root", uname.Data()) );
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/macro/pPb8160_Pbgoing_noTrkMax/oEmbedding_pPb8160_Pbgoing_80.root", uname.Data()) );

    //
    // Pythia
    //
    TFile *pPb8160PythiaFile = nullptr;
    if ( direction < 2 ) {
        pPb8160PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_%seta19.root", uname.Data(), directionStr.Data(), directionStr.Data(), cutTypeStr.Data()) );
        if ( !pPb8160PythiaFile ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_%seta19.root", uname.Data(), directionStr.Data(), directionStr.Data(), cutTypeStr.Data()) << std::endl;
            return;
        }
    }
    else {
        pPb8160PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_%seta19.root", uname.Data(), cutTypeStr.Data()) );
        if ( !pPb8160PythiaFile ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_%seta19.root", uname.Data(), cutTypeStr.Data()) << std::endl;
            return;
        }
    }

    //
    // Data
    //
    TFile *pPb8160DataFile = nullptr;
    if ( direction < 2 ) {
        pPb8160DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/exp/%s/%s_%s_ak4_%seta19.root", uname.Data(), dataDirectionStr.Data(), dataStr.Data(), dataDirectionStr.Data(), cutTypeStr.Data()) );
        if ( !pPb8160DataFile ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/exp/%s/%s_%s_ak4_%seta19.root", uname.Data(), dataDirectionStr.Data(), dataStr.Data(), dataDirectionStr.Data(), cutTypeStr.Data()) << std::endl;
            return;
        }
    }
    else {
        pPb8160DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/exp/%s_pPb8160_ak4_%seta19.root", uname.Data(), dataStr.Data(), cutTypeStr.Data()) );
        if ( !pPb8160DataFile ) {
            std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/exp/%s_pPb8160_ak4_%seta19.root", uname.Data(), dataStr.Data(), cutTypeStr.Data()) << std::endl;
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
    // dijetClosuresFrom3D( pPb8160EmbedFile, collisionSystem, collisionEnergy, date );
    // dijetClosuresFrom3D( pPb8160PythiaFile, collisionSystem, collisionEnergy, date );

    //
    // Plot comparison of inclusive jet eta distributions to check/validate the JEC
    //
    // data2mcInclusiveJetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of dijet reco and ref to gen distributions
    //
    data2mcDijetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);
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
}
