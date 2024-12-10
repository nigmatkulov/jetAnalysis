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
void pPb5020_compare2published() {

    int ptBinLow{6};
    int ptBinHi{10};

    TFile *pubFile = TFile::Open("pPb5020/cms_dijet_eta_5TeV_pub.root");
    TFile *runBFile = TFile::Open("/Users/nigmatkulov/cernbox/ana/pPb5020/exp/RunB/PAEGJet_RunB_pPb5020_1.root");
    TFile *runDFile = TFile::Open("/Users/nigmatkulov/cernbox/ana/pPb5020/exp/RunD/PAEGJet_RunD_pPb5020_1.root");

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

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    setPadStyle();
    hPubDijetEta_55_75->Draw();
    hRunBDijetEta_55_75->Draw("same");
    hRunDDijetEta_55_75->Draw("same");
}