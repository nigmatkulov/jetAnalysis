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

// C++ headers
#include <iostream>
#include <vector>

//________________
void plotCMSHeader() {
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.65, 0.93, "pp 5.02 TeV");
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
void plotComparison(TCanvas *c, TH1D* h1, TH1D* h2, 
                    int ptLow=55, int ptHi=75,
                    const char* h1Name="h1", const char* h2Name="h2") {

    // Text 
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.06);

    // Number for plotting position
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0., 0.12};
    Double_t legX[2] = {0.4, 0.65};
    Double_t legY[2] = {0.2, 0.35};

    // Make ratios
    TH1D *hRatio = dynamic_cast<TH1D*>( h2->Clone("hRatio") );
    hRatio->Divide( h1 );

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot comparison
    //

    c->cd(1);

    setPadStyle();
    // Plot distributions
    h1->Draw();
    h2->Draw("same");
    // Set ranges
    h1->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    h1->GetYaxis()->SetTitle("dN/d#eta^{dijet}");
    h1->GetXaxis()->SetTitle("#eta^{dijet}");

    t.DrawLatexNDC(0.3, 0.85, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi));
    plotCMSHeader();
    

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    leg->AddEntry(h1, Form("%s", h1Name), "p");
    leg->AddEntry(h2, Form("%s", h2Name), "p");
    leg->Draw();

    //
    // Plot ratio
    //

    c->cd(2);
    setPadStyle();
    hRatio->Draw();
    hRatio->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetYaxis()->SetTitle( Form( "%s / %s", h2Name, h1Name ) );
    hRatio->GetXaxis()->SetTitle("#eta^{dijet}");

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    //leg->AddEntry(hRatio, Form( "%s / %s", h2Name, h1Name ), "p");
    leg->Draw();

    // Line at unity
    line = new TLine(xRange[0], 1.0, xRange[1], 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kMagenta);
    line->Draw();
}

//________________
void pp5020_compare2published() {
    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Dijet ptAve binning
    int dijetPtOldVals[7] {25, 55, 75, 95, 115, 150, 400};
    std::vector<int> dijetPtVals(dijetPtOldVals, dijetPtOldVals + sizeof(dijetPtOldVals) / sizeof(int));

    // Files
    TFile *pp5020PubFile = TFile::Open("pPb5020/cms_dijet_eta_5TeV_pub.root");
    if ( !pp5020PubFile ) {
        std::cerr << "File not found: pPb5020/cms_dijet_eta_5TeV_pub.root" << std::endl;
        return;
    }
    TFile *pp5020DataFile = TFile::Open("/Users/nigmatkulov/cernbox/ana/pp5020/exp/pp5020_2017.root");
    if ( !pp5020DataFile ) {
        std::cerr << "File not found: /Users/nigmatkulov/cernbox/ana/pp5020/exp/pp5020_2017.root" << std::endl;
        return;
    }
    TFile *pp5020PythiaFile = TFile::Open("/Users/nigmatkulov/cernbox/ana/pp5020/pythia/pp5020_pythia8.root");
    if ( !pp5020PythiaFile ) {
        std::cerr << "File not found: /Users/nigmatkulov/cernbox/ana/pp5020/pythia/pp5020_pythia8.root" << std::endl;
        return;
    }

    // Histograms 
    TH1D *hPubDijetEta[6];
    TH1D *hDataDijetEta[6];
    TH1D *hPythiaRecoDijetEta[6];
    TH1D *hPythiaGenDijetEta[6];

    // Loop over dijet pT bins
    for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++) {
        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];

        // Initialize histograms with nullptr first
        hPubDijetEta[i] = nullptr;
        hDataDijetEta[i] = nullptr;
        hPythiaRecoDijetEta[i] = nullptr;
        hPythiaGenDijetEta[i] = nullptr;

        // h1lished data
        if ( i != 0 ) {
            hPubDijetEta[i] = dynamic_cast<TH1D*>( pp5020PubFile->Get( Form("ppEta_pt_%d_%d", ptLow, ptHi) ) );
            hPubDijetEta[i]->SetName( Form("hPubDijetEta_%d", i) );
            set1DStyle( hPubDijetEta[i], 2 );
        }

        hDataDijetEta[i] = dynamic_cast<TH1D*>( pp5020DataFile->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hDataDijetEta[i]->SetName( Form("hDataDijetEta_%d", i) );
        rescaleEta( hDataDijetEta[i] );
        set1DStyle( hDataDijetEta[i], 0 );

        hPythiaRecoDijetEta[i] = dynamic_cast<TH1D*>( pp5020PythiaFile->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hPythiaRecoDijetEta[i]->SetName( Form("hPythiaRecoDijetEta_%d", i) );
        rescaleEta( hPythiaRecoDijetEta[i] );
        set1DStyle( hPythiaRecoDijetEta[i], 1 );

        hPythiaGenDijetEta[i] = dynamic_cast<TH1D*>( pp5020PythiaFile->Get( Form("hGenDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hPythiaGenDijetEta[i]->SetName( Form("hPythiaGenDijetEta_%d", i) );
        rescaleEta( hPythiaGenDijetEta[i] );
        set1DStyle( hPythiaGenDijetEta[i], 3 );
    } // for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++)

    // Plot comparisons

    TCanvas *cData2DataComparison[5];
    for (int i{1}; i < dijetPtVals.size() - 2; i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{500}, canvY{1000};
        cData2DataComparison[i-1] = new TCanvas( Form("cData2DataComparison_%d", i-1), 
                                                 Form("cData2DataComparison_%d", i-1), 
                                                 canvX, canvY );

        plotComparison(cData2DataComparison[i-1], hPubDijetEta[i], hDataDijetEta[i], 
                       ptLow, ptHi, "pp5020 pub.", "pp5020 my");
    } // for (int i{1}; i < dijetPtVals.size() - 2; i++)
}