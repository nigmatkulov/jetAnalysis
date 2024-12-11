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
    t.DrawLatexNDC(0.65, 0.93, "pPb 5.02 TeV");
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
void plotComparison(TCanvas *c, TH1D* pub, TH1D* runB, TH1D* runD = nullptr,
                    int ptLow=55, int ptHi=75, TString jetAlgo="akCs4PF") {

    // Text 
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.06);

    // Number for plotting position
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0., 0.12};
    Double_t legX[2] = {0.4, 0.65};
    Double_t legY[2] = {0.2, 0.35};

    // Set style for the data points
    Int_t pubStyle{2};
    Int_t runBStyle{0};
    Int_t runDStyle{1};
    set1DStyle(pub, pubStyle);
    set1DStyle(runB, runBStyle);
    if ( runD ) {
        set1DStyle(runD, runDStyle);
    }

    // Make ratios
    TH1D *ratioB = dynamic_cast<TH1D*>( runB->Clone("ratioB") );
    ratioB->Divide( pub );

    TH1D *ratioD = nullptr;
    if ( runD ) {
        ratioD = dynamic_cast<TH1D*>( runD->Clone("ratioD") );
        ratioD->Divide( pub );
    }

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot comparison
    //
    c->cd(1);
    // Set pad style
    setPadStyle();
    // Plot distributions
    pub->Draw();
    runB->Draw("same");
    if ( runD ) {
        runD->Draw("same");
    }
    // Set ranges
    pub->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    pub->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    pub->GetYaxis()->SetTitle("dN/d#eta^{dijet}");
    pub->GetXaxis()->SetTitle("#eta^{dijet}");

    t.DrawLatexNDC(0.3, 0.85, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi));
    t.DrawLatexNDC(0.75, 0.85, Form("%s", jetAlgo.Data()));
    plotCMSHeader();
    

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    leg->AddEntry(pub, "pPb pub.", "p");
    leg->AddEntry(runB, "runB", "p");
    if ( runD ) {
        leg->AddEntry(runD, "runD", "p");
    }
    leg->Draw();

    c->cd(2);
    setPadStyle();
    ratioB->Draw();
    if ( runD ) {
        ratioD->Draw("same");
    }
    ratioB->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    ratioB->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratioB->GetYaxis()->SetTitle("Ratio to pub.");
    ratioB->GetXaxis()->SetTitle("#eta^{dijet}");

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    leg->AddEntry(ratioB, "runB", "p");
    if ( runD ) {
        leg->AddEntry(ratioD, "runD", "p");
    }
    leg->Draw();

    // Line at unity
    line = new TLine(xRange[0], 1.0, xRange[1], 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kMagenta);
    line->Draw();
}

//________________
void pPb5020_compare2published() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    int ptBinLow{6};
    int ptBinHi{10};

    TString jetAlgo = "akCs4PF";
    // TString jetAlgo = "ak4PF";

    TFile *pubFile = TFile::Open("pPb5020/cms_dijet_eta_5TeV_pub.root");
    TFile *runBFile = TFile::Open( Form("/Users/nigmatkulov/cernbox/ana/pPb5020/exp/PAEGJet_pPb5020_%s.root", jetAlgo.Data()) );
    // TFile *runBFile = TFile::Open( Form("/Users/nigmatkulov/cernbox/ana/pPb5020/exp/RunB/PAEGJet_RunB_pPb5020_%s.root", jetAlgo.Data()) );

    TFile *runDFile = TFile::Open( Form("/Users/nigmatkulov/cernbox/ana/pPb5020/exp/RunD/PAEGJet_RunD_pPb5020_%s.root", jetAlgo.Data()) );
    // TFile *runDFile = TFile::Open( Form("/Users/nigmatkulov/cernbox/ana/pPb5020/exp/PAEGJet_pPb5020_%s.root", jetAlgo.Data()) ); ;

    TH1D *hPubDijetEta_55_75 = dynamic_cast<TH1D*>( pubFile->Get("pPbEta_pt_55_75") );
    hPubDijetEta_55_75->SetName("hPubDijetEta_55_75");
    set1DStyle( hPubDijetEta_55_75, 2);

    TH3D *hRunBPtEtaDphi = dynamic_cast<TH3D*>( runBFile->Get("hRecoDijetPtEtaDphiWeighted") );
    hRunBPtEtaDphi->SetName("hRunBPtEtaDphi");
    TH1D *hRunBDijetEta_55_75 = dynamic_cast<TH1D*>( hRunBPtEtaDphi->ProjectionY("hRunBDijetEta_55_75", ptBinLow, ptBinHi) );
    rescaleEta( hRunBDijetEta_55_75 );
    set1DStyle( hRunBDijetEta_55_75, 0);

    TH3D *hRunDPtEtaDphi = dynamic_cast<TH3D*>( runDFile->Get("hRecoDijetPtEtaDphiWeighted") );
    hRunDPtEtaDphi->SetName("hRunDPtEtaDphi");
    TH1D *hRunDDijetEta_55_75 = dynamic_cast<TH1D*>( hRunDPtEtaDphi->ProjectionY("hRunDDijetEta_55_75", ptBinLow, ptBinHi) );
    rescaleEta( hRunDDijetEta_55_75 );
    set1DStyle( hRunDDijetEta_55_75, 1);

    TCanvas *c = new TCanvas("c", "c", 600, 1200);
    c->Divide(1, 2);

    plotComparison(c, hPubDijetEta_55_75, hRunBDijetEta_55_75, nullptr, 55, 75, jetAlgo);
    // plotComparison(c, hPubDijetEta_55_75, hRunBDijetEta_55_75, hRunDDijetEta_55_75, 55, 75, jetAlgo);

}