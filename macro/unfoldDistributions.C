// ROOT headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THnSparseD.h"

// RooUnfold
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// C++ headers
#include <iostream>

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.10);
    gPad->SetLeftMargin(0.15);
}

//________________
void set1DStyle(TH1 *h, Int_t weight = 0) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.2;
    Int_t lineWidth = 2;
    Int_t color = 2;
    if (weight) {
        color = 4;
        markerStyle = 24;
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
    h->GetYaxis()->SetTitleOffset(1.0);
}

//________________
void set2DStyle(TH2* h) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.0);
}

//________________
void performUnfolding(TFile *inFile, TString date) {

    // Ref
    TH1D *hRefDijetEta = (TH1D*)inFile->Get("hRefDijetEta");
    // Reco
    TH1D *hRecoDijetEta = (TH1D*)inFile->Get("hRecoDijetEta");
    // Response matrix
    TH2D *hDijetEtaRefVsReco = (TH2D*)inFile->Get("hRefDijetEtaVsRecoDijetEta");
    hDijetEtaRefVsReco->SetName("hDijetEtaRefVsReco");

    // Gen
    THnSparse *hGen = (THnSparseD*)inFile->Get("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted");
    TH1D* hGenDijetEta = (TH1D*)hGen->Projection(1);
    // Miss = Gen - Ref
    TH1D* hMiss = new TH1D("hMiss","hMiss;Miss #eta^{dijet}",
                           hDijetEtaRefVsReco->GetNbinsX(), 
                           hDijetEtaRefVsReco->GetXaxis()->GetBinLowEdge(1),
                           hDijetEtaRefVsReco->GetXaxis()->GetBinUpEdge( hDijetEtaRefVsReco->GetNbinsX() ) );
    hMiss->Sumw2();
    hMiss->Add(hGenDijetEta, hRefDijetEta, 1., -1.);

    // Create response matrix

}

//________________
void unfoldDistributions(const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing.root") {

    gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TString date {"20240412"};
    TFile *inFile = TFile::Open(inFileName);

}