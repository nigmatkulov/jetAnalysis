#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TH1.h"
#include "TStyle.h"

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
    gPad->GetFrame()->SetLineWidth(2);
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
    else if (type == 5) {
        color = 8;        // green
        markerStyle = 39; // radiation
    }
    else if ( type == 6 ) {
        color = 9;        // purple
        markerStyle = 20; // circle
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
void setSystUncrtStyle(TH1* h, Int_t type = 0) {
    Int_t color = 29;        // Green-like
    Int_t markerColor = 1;
    Int_t lineWidth = 0;
    Int_t markerStyle = 20; // Full circle
    if ( type == 0 ) {
        color = 29;          // Green-like
        lineWidth = 1;
        markerColor = 1;     // Black
        markerStyle = 20;    // Full circle
    }
    else if ( type == 1) {
        color = 2;          // Red
        lineWidth = 1;
        markerColor = 2;    // Red           
        markerStyle = 21;   // Full circle
    }
    else if ( type == 2 ) {
        color = 4;           // Blue
        lineWidth = 3;
        markerColor = 4;     // Blue
        // markerStyle = 24;    // Open circle
        markerStyle = 71;    // Open circle
    }
    else {
        color = 46;          // Dark red-like
        lineWidth = 1;
        markerStyle = 22;
    }
    h->SetFillColorAlpha(color, 0.35);
    h->SetLineColor(color);
    h->SetMarkerColor(markerColor);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(1.8);
    h->SetLineWidth(lineWidth);
    h->SetFillStyle(1001);           // Solid fill

    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
    h->GetYaxis()->SetTitleOffset(1.1);
}

//________________
void plotCuts(Bool_t isCM = kFALSE) {

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
    c1->cd();
    setPadStyle();

    const Char_t *dijet_dphi = "#frac{5#pi}{6}";

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);
    t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
    t.DrawLatexNDC(0.2, 0.6, "p_{T}^{Subleading} > 40 GeV");
    if ( isCM ) {
        t.DrawLatexNDC(0.2, 0.45, "|#eta_{CM}| < 2.5");
    }
    else {
        t.DrawLatexNDC(0.2, 0.45, "|#eta| < 3");
    }
    t.DrawLatexNDC(0.2, 0.3, Form("#Delta#phi^{dijet} > %s", dijet_dphi));
}

//________________
void plotLegend(Bool_t isRatioToMc = kFALSE) {
    

    TH1D *hData = new TH1D("hData", "Data", 1, 0, 1);
    hData->Sumw2();
    setSystUncrtStyle(hData, 0);

    TH1D *hSyst = new TH1D("hSyst", "Syst. Uncrt.", 1, 0, 1);
    hSyst->Sumw2();
    setSystUncrtStyle(hSyst, 0);

    TH1D *hMc = new TH1D("hMc", "MC", 1, 0, 1);
    hMc->Sumw2();
    setSystUncrtStyle(hMc, 2);

    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 1000);
    c2->cd();
    setPadStyle();

    Double_t xLegend[2] {0.2, 0.7};
    Double_t yLegend[2] {0.5, 0.85};
    if ( isRatioToMc ) {
        yLegend[0] = 0.45;
        yLegend[1] = 0.85;
    }

    TLegend *legend = new TLegend(xLegend[0], yLegend[0], xLegend[1], yLegend[1]);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.06);
    if ( isRatioToMc ) {
        legend->AddEntry(hData, "Data / PYTHIA", "p");
        legend->AddEntry(hMc, "MC / PYTHIA", "p");
        legend->AddEntry(hSyst, "Syst. uncrt.", "f");
    }
    else {
        legend->AddEntry(hData, "Data", "p");
        legend->AddEntry(hSyst, "Syst. uncrt.", "f");
    }
    legend->Draw();
}

//________________
void plotCanvWithInfo() {
    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    Bool_t isRatioToMc = kFALSE;
    Bool_t isCM = kFALSE;

    plotCuts(isCM);
    plotLegend(isRatioToMc);
}