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
void plotEfficiency(TFile *inFile, TString date) {

    // Rebinning
    Int_t rebinX{1}, rebinY{1};

    // Read gen jet acceptance
    TH2D *hGenPtVsEta = (TH2D*)inFile->Get("hGenInclusiveJetPtEta");
    hGenPtVsEta->SetName("hGenPtVsEta");

    // Read reco that matched gen acceptance
    TH2D *hRecoPtVsEta = (TH2D*)inFile->Get("hRecoMatchedPtEta");
    hRecoPtVsEta->SetName("hRecoPtVsEta");

    // Preform rebinning
    hGenPtVsEta->RebinX( rebinX );
    hGenPtVsEta->RebinY( rebinY );

    hRecoPtVsEta->RebinX( rebinX );
    hRecoPtVsEta->RebinY( rebinY );

    // Create plot for efficiency
    TH2D *hEfficiency = new TH2D("hEfficiency","Inclusive Reco (matched) / Gen;#eta;p_{T} (GeV/c)",
                                 hGenPtVsEta->GetNbinsX(), 
                                 hGenPtVsEta->GetXaxis()->GetBinLowEdge(1),
                                 hGenPtVsEta->GetXaxis()->GetBinUpEdge( hGenPtVsEta->GetNbinsX() ),
                                 hGenPtVsEta->GetNbinsY(), 
                                 hGenPtVsEta->GetYaxis()->GetBinLowEdge(1),
                                 hGenPtVsEta->GetYaxis()->GetBinUpEdge( hGenPtVsEta->GetNbinsX() ) );
    hEfficiency->Sumw2();
    hEfficiency->Divide(hRecoPtVsEta, hGenPtVsEta, 1., 1.);
    set2DStyle( hEfficiency );

    // Plot efficiency
    TCanvas *cEfficiency = new TCanvas("cEfficiency","cEfficiency", 1200, 800);
    setPadStyle();
    hEfficiency->Draw("colz");
    gPad->SetLogz(1);
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
}

//________________
void pPb_embedding_qa(const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing.root") {

    gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TString date {"20240412"};
    TFile *inFile = TFile::Open(inFileName);

    // Plot jet reconstruction efficiency as a function of acceptance (pT vs eta)
    plotEfficiency(inFile, date);

    // Plot correlation between ref and reco dijet eta
    plotEtaDijetCorrelation(inFile, date);
}