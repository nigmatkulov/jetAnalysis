// ROOT headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THnSparse.h"

// RooUnfold
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// C++ headers
#include <iostream>
#include <vector>

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.10);
    gPad->SetLeftMargin(0.15);
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.2;
    Int_t lineWidth = 2;
    Int_t color = 2;
    if (type == 0) {
        color = 2;
        markerStyle = 20;
    }
    else if (type == 1) {
        color = 4;
        markerStyle = 24;
    }
    else if (type == 2) {
        color = 1;
        markerStyle = 22;
    }
    else if (type == 3) {
        color = 6;
        markerStyle = 26;
    }
    else if (type == 4) {
        color = 3;
        markerStyle = 29;
    }
    else {
        color = 9;
        markerStyle = 30;
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
void unfold1D(TH1D* hReco, TH1D *hRef, TH2D* hResponse, TH1D* hGen = nullptr, TH1D* hUnfold = nullptr, 
              TString date = "20240418", TString name = "unfold", Int_t nIter = 4) {
    Int_t recoType{0};
    Int_t refType{1};
    Int_t unfoldType{2};
    Int_t genType{3};

    // Create response
    RooUnfoldResponse response( hReco, hRef, hResponse, 
                                Form("%s", name.Data()), Form("%s", name.Data()) );

    // Create unfolding procedure 
    // (const RooUnfoldResponseT< Hist, Hist2D > *res, const Hist *meas, Int_t niter=4, 
    //  Bool_t smoothit=false, Bool_t handleFakes=false, const char *name=0, const char *title=0)
    RooUnfoldBayes unfold( &response, hReco, nIter, kFALSE, 
                           Form("%s_UnfoldBayes", name.Data()), Form("%s_UnfoldBayes", name.Data()) );

    std::cout << "-------------------------------------------" << std::endl;
    std::cout << Form("Unfolding for \t%s \n", name.Data() );

    // Create unfolded histogram

    hUnfold = (TH1D*)unfold.Hunfold();
    hUnfold->SetNameTitle( Form("hUnfold_%s", name.Data()), 
                           Form("hUnfold_%s", name.Data()) );

    // Create histograms for ratios to ref
    TH1D *hReco2RefRatio = new TH1D( Form("hReco2RefRatio_%s", name.Data()), 
                                     Form("hReco2RefRatio_%s", name.Data()), 
                                     hReco->GetNbinsX(), 
                                     hReco->GetXaxis()->GetBinLowEdge(1),
                                     hReco->GetXaxis()->GetBinUpEdge( hReco->GetNbinsX() ) );
    hReco2RefRatio->Sumw2();
    TH1D *hUnfold2RefRatio = new TH1D( Form("hUnfold2RefRatio_%s", name.Data()), 
                                       Form("hUnfold2RefRatio_%s", name.Data()), 
                                       hUnfold->GetNbinsX(), 
                                       hUnfold->GetXaxis()->GetBinLowEdge(1),
                                       hUnfold->GetXaxis()->GetBinUpEdge( hUnfold->GetNbinsX() ) );
    hUnfold2RefRatio->Sumw2();

    TH1D *hGen2RefRatio = new TH1D( Form("hGen2RefRatio_%s", name.Data()), 
                                    Form("hGen2RefRatio_%s", name.Data()), 
                                    hGen->GetNbinsX(), 
                                    hGen->GetXaxis()->GetBinLowEdge(1),
                                    hGen->GetXaxis()->GetBinUpEdge( hGen->GetNbinsX() ) );
    hGen2RefRatio->Sumw2();

    // Fill ratios
    hReco2RefRatio->Divide(hReco, hRef, 1., 1.);
    hUnfold2RefRatio->Divide(hUnfold, hRef, 1., 1.);
    if ( hGen ) {
        hGen2RefRatio->Divide(hGen, hRef, 1., 1.);
    }

    // Create canvas for plotting
    TCanvas *canv = new TCanvas( Form("pPb8160_%s", name.Data()), 
                                 Form("pPb8160_%s", name.Data()),
                                 800, 800);
    canv->Divide(1, 2);

    // Plot distributions
    canv->cd(1);
    setPadStyle();
    set1DStyle(hReco, recoType);
    set1DStyle(hRef, refType);
    set1DStyle(hUnfold, unfoldType);
    if (hGen) {
        set1DStyle(hGen, genType);
    }
    hReco->Draw();
    hRef->Draw("same");
    hUnfold->Draw("same");
    if ( hGen ) {
        hGen->Draw("same");
    }
    hReco->GetYaxis()->SetRangeUser(0., 0.0005);

    TLegend *leg = new TLegend(0.75, 0.7, 0.85, 0.9);
    leg->SetLineWidth(0);
    leg->AddEntry(hReco,Form("Reco"), "p");
    leg->AddEntry(hRef,Form("Ref"), "p");
    leg->AddEntry(hUnfold,Form("Unfold"), "p");
    if ( hGen ) {
        leg->AddEntry(hGen,Form("Gen"), "p");
    }
    leg->SetTextSize(0.06);
    leg->Draw();

    // Plot ratios to ref
    canv->cd(2);
    setPadStyle();
    set1DStyle(hReco2RefRatio, recoType);
    set1DStyle(hUnfold2RefRatio, unfoldType);
    if ( hGen ) {
        set1DStyle(hGen2RefRatio, genType);
    }
    hReco2RefRatio->Draw();
    hUnfold2RefRatio->Draw("same");
    if ( hGen ) {
        hGen2RefRatio->Draw("same");
    }

    hReco2RefRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    hReco2RefRatio->GetYaxis()->SetTitle("Ratio to ref");

    canv->SaveAs( Form("%s/pPb8160_%s.pdf", date.Data(), name.Data()) );

    auto* R = response.Hresponse();
    auto* cResponse = new TCanvas(Form("cResponse_%s", name.Data()), 
                                  Form("cResponse_%s", name.Data()), 
                                  800, 800);
    setPadStyle();
    R->SetStats(0);
    R->Draw("colz");
    gPad->SetLogz(1);
    cResponse->SaveAs(Form("%s/pPb8160_%s_response.pdf", date.Data(), name.Data()));
}

//________________
void unfold1DPtHat(TH1D* hReco, TH1D *hRef, TH2D* hResponse, TH1D* hGen = nullptr, TH1D* hUnfold = nullptr, 
                   TString date = "20240418", TString name = "unfold", Int_t nIter = 4, Int_t ptHat=50, 
                   TH1D *hRefOrig = nullptr) {
    
    Int_t recoType{0};
    Int_t refType{1};
    Int_t unfoldType{2};
    Int_t genType{3};

    // Create response
    RooUnfoldResponse response( hReco, hRef, hResponse, 
                                Form("%s", name.Data()), Form("%s", name.Data()) );

    // Create unfolding procedure 
    // (const RooUnfoldResponseT< Hist, Hist2D > *res, const Hist *meas, Int_t niter=4, 
    //  Bool_t smoothit=false, Bool_t handleFakes=false, const char *name=0, const char *title=0)
    RooUnfoldBayes unfold( &response, hReco, nIter, kFALSE, 
                           Form("%s_UnfoldBayes", name.Data()), Form("%s_UnfoldBayes", name.Data()) );

    std::cout << "-------------------------------------------" << std::endl;
    std::cout << Form("Unfolding for \t%s \n", name.Data() );

    // Create unfolded histogram

    hUnfold = (TH1D*)unfold.Hunfold();
    hUnfold->SetNameTitle( Form("hUnfold_%s", name.Data()), 
                           Form("hUnfold_%s", name.Data()) );

    // Create histograms for ratios to ref
    TH1D *hReco2RefRatio = new TH1D( Form("hReco2RefRatio_%s", name.Data()), 
                                     Form("hReco2RefRatio_%s", name.Data()), 
                                     hReco->GetNbinsX(), 
                                     hReco->GetXaxis()->GetBinLowEdge(1),
                                     hReco->GetXaxis()->GetBinUpEdge( hReco->GetNbinsX() ) );
    hReco2RefRatio->Sumw2();
    TH1D *hUnfold2RefRatio = new TH1D( Form("hUnfold2RefRatio_%s", name.Data()), 
                                       Form("hUnfold2RefRatio_%s", name.Data()), 
                                       hUnfold->GetNbinsX(), 
                                       hUnfold->GetXaxis()->GetBinLowEdge(1),
                                       hUnfold->GetXaxis()->GetBinUpEdge( hUnfold->GetNbinsX() ) );
    hUnfold2RefRatio->Sumw2();
    TH1D *hGen2RefRatio = new TH1D( Form("hGen2RefRatio_%s", name.Data()), 
                                    Form("hGen2RefRatio_%s", name.Data()), 
                                    hGen->GetNbinsX(), 
                                    hGen->GetXaxis()->GetBinLowEdge(1),
                                    hGen->GetXaxis()->GetBinUpEdge( hGen->GetNbinsX() ) );
    hGen2RefRatio->Sumw2();

    // Fill ratios
    hReco2RefRatio->Divide(hReco, hRef, 1., 1.);
    hUnfold2RefRatio->Divide(hUnfold, hRef, 1., 1.);
    if ( hGen ) {
        hGen2RefRatio->Divide(hGen, hRef, 1., 1.);
    }

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Create canvas for plotting
    TCanvas *canv = new TCanvas( Form("pPb8160_%s", name.Data()), 
                                 Form("pPb8160_%s", name.Data()),
                                 800, 800);
    canv->Divide(1, 2);

    // Plot distributions
    canv->cd(1);
    setPadStyle();
    set1DStyle(hReco, recoType);
    set1DStyle(hRef, refType);
    set1DStyle(hUnfold, unfoldType);
    if ( hGen ) {
        set1DStyle(hGen, genType);
    }
    hReco->Draw();
    hRef->Draw("same");
    hUnfold->Draw("same");
    if ( hGen ) {
        hGen->Draw("same");
    }
    hReco->GetYaxis()->SetRangeUser(0., 0.0005);
    if ( hRefOrig ) {
        set1DStyle(hRefOrig, 4);
        hRefOrig->Draw("same");
    }

    TLegend *leg = new TLegend(0.7, 0.68, 0.85, 0.92);
    leg->SetLineWidth(0);
    leg->AddEntry(hReco,Form("Reco"), "p");
    leg->AddEntry(hRef,Form("Ref"), "p");
    leg->AddEntry(hUnfold,Form("Unfold"), "p");
    if ( hGen ) {
        leg->AddEntry(hGen,Form("Gen"), "p");
    }
    if (hRefOrig) {
        leg->AddEntry(hRefOrig,Form("Ref (all #hat{p_{T}})"), "p");
    }
    leg->SetTextSize(0.06);
    leg->Draw();

    t.DrawLatexNDC(0.2, 0.9, Form("#hat{p_{T}}>%d GeV/c", ptHat));

    // Plot ratios to ref
    canv->cd(2);
    setPadStyle();
    set1DStyle(hReco2RefRatio, recoType);
    set1DStyle(hUnfold2RefRatio, unfoldType);
    if ( hGen ) {
        set1DStyle(hGen2RefRatio, genType);
    }
    hReco2RefRatio->Draw();
    hUnfold2RefRatio->Draw("same");
    if ( hGen ) {
        hGen2RefRatio->Draw("same");
    }

    hReco2RefRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    hReco2RefRatio->GetYaxis()->SetTitle("Ratio to ref");

    canv->SaveAs( Form("%s/pPb8160_%s.pdf", date.Data(), name.Data()) );

    auto* R = response.Hresponse();
    auto* cResponse = new TCanvas(Form("cResponse_%s", name.Data()), 
                                  Form("cResponse_%s", name.Data()), 
                                  800, 800);
    setPadStyle();
    R->SetStats(0);
    R->Draw("colz");
    gPad->SetLogz(1);
    cResponse->SaveAs(Form("%s/pPb8160_%s_response.pdf", date.Data(), name.Data()));
    t.DrawLatexNDC(0.2, 0.9, Form("All #hat{p_{T}}"));
}

//________________
void unfoldDijetEta1D(TFile *inFile, TString date) {


    // Retrieve THnSparse for reco 2 ref
    THnSparseD *hReco2RefDijet = (THnSparseD*)inFile->Get("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    // Retrieve THnSparse for gen
    THnSparseD *hGenDijet = (THnSparseD*)inFile->Get("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted");

    // Define single-jet binning
    Int_t jetPtBins{50}, jetEtaBins{50};      // fPtRange{20., 520.}, fEtaRange{-5.0, 5.0}
    Double_t jetPtStep{10}, jetEtaStep{0.2};
    // Define dijet binning
    Int_t dijetPtBins{120}, dijetEtaBins{50};      // fDijetPtRange{20., 620.},, fEtaRange{-5.0, 5.0}
    Double_t dijetPtStep{10}, dijetEtaStep{0.2};   // fDijetEtaRange{-5., 5.}

    // Define bins to analyze

    // > 50 GeV/c
    std::vector<Int_t> ptLeadLow{4};  
    std::vector<Int_t> ptLeadHi{50};

    // > 30 GeV/c
    std::vector<Int_t> ptSubLeadLow{2};
    std::vector<Int_t> ptSubLeadHi{50};

    // Integrated (>20 GeV/c)
    std::vector<Int_t> ptDijetLow{1};
    std::vector<Int_t> ptDijetHi{50};

    // Refine reco2ref distribution

    // Reco dijet pT
    hReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(0), ptDijetHi.at(0) );
    // Reco leading jet pT
    hReco2RefDijet->GetAxis(2)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Reco subleading jet pT
    hReco2RefDijet->GetAxis(4)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );


    // Refine gen distribution

    // Gen leading jet pT
    hGenDijet->GetAxis(3)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Gen subleading jet pT
    hGenDijet->GetAxis(6)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    //
    // Make projections
    //

    //
    // Eta dijet
    //

    // Reco and ref
    TH1D *hRecoDijetEta = (TH1D*)hReco2RefDijet->Projection(1);
    TH1D *hRefDijetEta = (TH1D*)hReco2RefDijet->Projection(7);
    TH2D *hRef2RecoDijetEta = (TH2D*)hReco2RefDijet->Projection(7, 1);
    // Gen
    TH1D *hGenDijetEta = (TH1D*)hGenDijet->Projection(1);
    TH1D *hUnfoldEta = new TH1D();
    TString name = "TestEta";
    unfold1D(hRecoDijetEta, hRefDijetEta, hRef2RecoDijetEta, hGenDijetEta, hUnfoldEta, date, name, 4);

    //
    // pT dijet
    //

    // Reco and ref
    TH1D *hRecoDijetPt = (TH1D*)hReco2RefDijet->Projection(0);
    TH1D *hRefDijetPt = (TH1D*)hReco2RefDijet->Projection(6);
    TH2D *hRef2RecoDijetPt = (TH2D*)hReco2RefDijet->Projection(6, 0);
    // Gen
    TH1D *hGenDijetPt = (TH1D*)hGenDijet->Projection(0);
    TH1D *hUnfoldPt = new TH1D();
    name = "TestPt";
    unfold1D(hRecoDijetPt, hRefDijetPt, hRef2RecoDijetPt, hGenDijetPt, hUnfoldPt, date, name, 4);
}

//________________
void performUnfolding(TFile *inFile, TString date) {

    Int_t recoType{0};
    Int_t refType{1};
    Int_t unfoldType{2};
    Int_t genType{3};

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

    TH1D* hReco2RefRatio = new TH1D("hReco2RefRatio","Reco / Ref; #eta^{dijet}; Ratio to Ref",
                                    hDijetEtaRefVsReco->GetNbinsX(), 
                                    hDijetEtaRefVsReco->GetXaxis()->GetBinLowEdge(1),
                                    hDijetEtaRefVsReco->GetXaxis()->GetBinUpEdge( hDijetEtaRefVsReco->GetNbinsX() ) );
    hReco2RefRatio->Sumw2();
    TH1D* hUnfold2RefRatio = new TH1D("hUnfold2RefRatio","Unfold / Ref; #eta^{dijet}; Ratio to Ref",
                                    hDijetEtaRefVsReco->GetNbinsX(), 
                                    hDijetEtaRefVsReco->GetXaxis()->GetBinLowEdge(1),
                                    hDijetEtaRefVsReco->GetXaxis()->GetBinUpEdge( hDijetEtaRefVsReco->GetNbinsX() ) );
    hUnfold2RefRatio->Sumw2();

    // Create response matrix
    //RooUnfoldResponse response( hRecoDijetEta, hMiss, hDijetEtaRefVsReco, "response", "Response matrix");
    RooUnfoldResponse response( hRecoDijetEta, hGenDijetEta, hDijetEtaRefVsReco, "response", "Response matrix");

    // Create unfolding procedure 
    // (const RooUnfoldResponseT< Hist, Hist2D > *res, const Hist *meas, Int_t niter=4, 
    //  Bool_t smoothit=false, Bool_t handleFakes=fales, const char *name=0, const char *title=0)
    RooUnfoldBayes unfold( &response, hRecoDijetEta, 4);

    // Create unfolded histogram
    TH1D* hUnfold = (TH1D*)unfold.Hunfold();
    hUnfold->SetNameTitle("hUnfold","Unfolded distribution;Unfolded #eta^{dijet}");


    // Create canvas to present results
    TCanvas *cUndold = new TCanvas("cUnfold", "Unfolding", 800, 800);
    cUndold->Divide(1, 2);

    cUndold->cd(1);
    setPadStyle();
    set1DStyle(hRecoDijetEta, recoType);
    set1DStyle(hRefDijetEta, refType);
    set1DStyle(hUnfold, unfoldType);

    hRecoDijetEta->Draw();
    hRecoDijetEta->GetXaxis()->SetTitle("#eta^{dijet}");
    hRecoDijetEta->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    hRefDijetEta->Draw("same");
    hUnfold->Draw("same");

    TLegend *leg = new TLegend(0.75, 0.7, 0.85, 0.9);
    leg->SetLineWidth(0);
    leg->AddEntry(hRecoDijetEta,Form("Reco"), "p");
    leg->AddEntry(hRefDijetEta,Form("Ref"), "p");
    leg->AddEntry(hUnfold,Form("Unfold"), "p");
    leg->SetTextSize(0.05);
    leg->Draw();

    cUndold->cd(2);
    setPadStyle();
    hReco2RefRatio->Divide(hRecoDijetEta, hRefDijetEta);
    set1DStyle(hReco2RefRatio, recoType);
    hUnfold2RefRatio->Divide(hUnfold, hRefDijetEta);
    set1DStyle(hUnfold2RefRatio, unfoldType);

    hReco2RefRatio->Draw();
    hUnfold2RefRatio->Draw("same");
    hReco2RefRatio->GetYaxis()->SetRangeUser(0.8, 1.2);

    unfold.PrintTable(std::cout, hRefDijetEta);

    cUndold->SaveAs(Form("%s/pPb8160_etaDijet_1D_unfolding.pdf", date.Data()));
}

//________________
void unfoldDifferentPtHat(TFile *inFile, TFile *inFilePtHat, Int_t ptHat, TString date) {

    // Retrieve THnSparse for reco 2 ref
    THnSparseD *hReco2RefDijet = (THnSparseD*)inFile->Get("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    hReco2RefDijet->SetName("hReco2RefDijet");
    // Retrieve THnSparse for gen
    THnSparseD *hReco2RefDijetPtHat = (THnSparseD*)inFilePtHat->Get("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    // Retrieve THnSparse for gen
    THnSparseD *hGenDijetPtHat = (THnSparseD*)inFilePtHat->Get("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted");

    // Define single-jet binning
    Int_t jetPtBins{50}, jetEtaBins{50};      // fPtRange{20., 520.}, fEtaRange{-5.0, 5.0}
    Double_t jetPtStep{10}, jetEtaStep{0.2};
    // Define dijet binning
    Int_t dijetPtBins{120}, dijetEtaBins{50};      // fDijetPtRange{20., 620.},, fEtaRange{-5.0, 5.0}
    Double_t dijetPtStep{10}, dijetEtaStep{0.2};   // fDijetEtaRange{-5., 5.}

    // Define bins to analyze

    // > 50 GeV/c
    std::vector<Int_t> ptLeadLow{4};  
    std::vector<Int_t> ptLeadHi{50};

    // > 30 GeV/c
    std::vector<Int_t> ptSubLeadLow{2};
    std::vector<Int_t> ptSubLeadHi{50};

    // Integrated (>20 GeV/c)
    std::vector<Int_t> ptDijetLow{1};
    std::vector<Int_t> ptDijetHi{50};

    // Refine reco2ref distribution for all ptHat

    // Reco dijet pT
    hReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(0), ptDijetHi.at(0) );
    // Reco leading jet pT
    hReco2RefDijet->GetAxis(2)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Reco subleading jet pT
    hReco2RefDijet->GetAxis(4)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    // Refine reco2ref distribution for a single ptHat

    // Reco dijet pT
    hReco2RefDijetPtHat->GetAxis(0)->SetRange( ptDijetLow.at(0), ptDijetHi.at(0) );
    // Reco leading jet pT
    hReco2RefDijetPtHat->GetAxis(2)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Reco subleading jet pT
    hReco2RefDijetPtHat->GetAxis(4)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    // Refine gen distribution for a single ptHat

    // TODO: Add range for pT of dijet!!!!

    // Gen leading jet pT
    hGenDijetPtHat->GetAxis(3)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Gen subleading jet pT
    hGenDijetPtHat->GetAxis(6)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    //
    // Make projections
    //

    //
    // Eta dijet
    //

    // Reco and ref
    TH1D *hRecoDijetPtHatEta = (TH1D*)hReco2RefDijetPtHat->Projection(1);
    TH1D *hRefDijetPtHatEta = (TH1D*)hReco2RefDijetPtHat->Projection(7);
    TH2D *hRef2RecoDijetEta = (TH2D*)hReco2RefDijet->Projection(7, 1);
    TH1D *hRefDijetEta = (TH1D*)hReco2RefDijet->Projection(7);
    // Gen
    TH1D *hGenDijetPtHatEta = (TH1D*)hGenDijetPtHat->Projection(1);
    TH1D *hUnfoldEta = new TH1D();
    TString name = "TestEtaPtHat50";
    unfold1DPtHat(hRecoDijetPtHatEta, hRefDijetPtHatEta, hRef2RecoDijetEta, 
                  hGenDijetPtHatEta, hUnfoldEta, date, name, 4, ptHat, hRefDijetEta);
}

//________________
void unfoldDistributions(const Char_t *dateToday = "20240418") {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TString date {dateToday};
    
    const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing.root";
    TFile *inFile = TFile::Open(inFileName);

    const Char_t *inFileNameWithPtHat = "../build/oEmbedding_pPb8160_Pbgoing_50.root";
    TFile *inFileWithPtHat = TFile::Open(inFileNameWithPtHat);

    // Run 1D unfolding
    //performUnfolding(inFile, date);

    // Run 1D unfolding with semi-automated regime
    //unfoldDijetEta1D(inFile, date);

    // Run 1D unfolding for ptHat distributions with the response matrix
    // for total ptHat
    unfoldDifferentPtHat(inFile, inFileWithPtHat, 50, date);
}