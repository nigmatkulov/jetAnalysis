// ROOT headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TSystemDirectory.h"
#include <sys/stat.h>
#include "TLine.h"
#include "TLatex.h"

// RooUnfold
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// C++ headers
#include <iostream>
#include <vector>

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
void rescaleEta(TH1* h) {
    for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++) {
        Double_t val = h->GetBinContent( iBin );
        Double_t valErr = h->GetBinError( iBin );
        Double_t binWidth = h->GetBinWidth( iBin );
        h->SetBinContent( iBin, val / binWidth );
        h->SetBinError( iBin, valErr / binWidth );
    } // for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++)
    h->Scale( 1. / h->Integral() );
}

//________________
void rescaleEta(TH2* h) {
    for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++) {
        for (Int_t jBin=1; jBin<=h->GetNbinsY(); jBin++) {
            Double_t val = h->GetBinContent( iBin, jBin );
            Double_t valErr = h->GetBinError( iBin, jBin );
            Double_t binWidthX = h->GetXaxis()->GetBinWidth( iBin );
            Double_t binWidthY = h->GetYaxis()->GetBinWidth( jBin );
            h->SetBinContent( iBin, jBin, val / (binWidthX * binWidthY) );
            h->SetBinError( iBin, jBin, valErr / (binWidthX * binWidthY) );
        } // for (Int_t jBin=1; jBin<=h->GetNbinsY(); jBin++)
    } // for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++)
    h->Scale( 1. / h->Integral() );
}

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.1;
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
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void set2DStyle(TH2* h, Bool_t doRenorm = kFALSE) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./ h->Integral() );
    }
    
}

//________________
void unfold1D(TH1D* hReco, TH1D *hRef, TH2D* hResponse, TH1D* hGen = nullptr, TH1D* hUnfold = nullptr, 
              TString date = "20240418", TString name = "unfold", Int_t nIter = 4) {
    Int_t recoType{0};
    Int_t refType{1};
    Int_t unfoldType{2};
    Int_t genType{3};

    set1DStyle(hReco, recoType);
    set1DStyle(hRef, refType);
    if ( hGen ) {
        set1DStyle(hGen, genType);
    }

    set2DStyle(hResponse);

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
    set1DStyle(hUnfold, unfoldType);

    std::cout << "ref: " << hRef->GetNbinsX() << " reco: " << hReco->GetNbinsX()
    << " gen: " << hGen->GetNbinsX() << " unfold: " << hUnfold->GetNbinsX() << std::endl;

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

    set1DStyle(hReco2RefRatio, recoType);
    set1DStyle(hUnfold2RefRatio, unfoldType);
    if ( hGen ) {
        set1DStyle(hGen2RefRatio, genType);
    }

    // Fill ratios
    hReco2RefRatio->Divide(hReco, hRef, 1., 1., "b");
    hUnfold2RefRatio->Divide(hUnfold, hRef, 1., 1., "b");
    if ( hGen ) {
        hGen2RefRatio->Divide(hGen, hRef, 1., 1., "b");
    }

    // hReco2RefRatio->Divide(hReco, hRef, 1., 1.);
    // hUnfold2RefRatio->Divide(hUnfold, hRef, 1., 1.);
    // if ( hGen ) {
    //     hGen2RefRatio->Divide(hGen, hRef, 1., 1.);
    // }

    // Create canvas for plotting
    TCanvas *canv = new TCanvas( Form("pPb8160_%s", name.Data()), 
                                 Form("pPb8160_%s", name.Data()),
                                 800, 800);
    canv->Divide(1, 2);

    // Plot distributions
    canv->cd(1);
    setPadStyle();
    // set1DStyle(hReco, recoType);
    // set1DStyle(hRef, refType);
    // set1DStyle(hUnfold, unfoldType);
    // if (hGen) {
    //     set1DStyle(hGen, genType);
    // }
    hReco->Draw();
    hRef->Draw("same");
    hUnfold->Draw("same");
    if ( hGen ) {
        hGen->Draw("same");
    }
    //hReco->GetYaxis()->SetRangeUser(0., 0.0005);

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

    set1DStyle(hReco, recoType);
    set1DStyle(hRef, refType);
    if ( hGen ) {
        set1DStyle(hGen, genType);
    }

    set2DStyle(hResponse);

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
    set1DStyle(hUnfold, unfoldType);

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

    set1DStyle(hReco2RefRatio, recoType);
    set1DStyle(hUnfold2RefRatio, unfoldType);
    if ( hGen ) {
        set1DStyle(hGen2RefRatio, genType);
    }

    // Fill ratios
    hReco2RefRatio->Divide(hReco, hRef, 1., 1., "b");
    hUnfold2RefRatio->Divide(hUnfold, hRef, 1., 1., "b");
    if ( hGen ) {
        hGen2RefRatio->Divide(hGen, hRef, 1., 1., "b");
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

    hReco->Draw();
    hRef->Draw("same");
    hUnfold->Draw("same");
    if ( hGen ) {
        hGen->Draw("same");
    }
    //hReco->GetYaxis()->SetRangeUser(0., 0.0005);
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
    t.DrawLatexNDC(0.2, 0.9, Form("All #hat{p_{T}}"));
    cResponse->SaveAs(Form("%s/pPb8160_%s_response.pdf", date.Data(), name.Data()));
    
}

//________________
void unfoldDijetEta1D(TFile *inFile, TString date) {


    // Retrieve THnSparse for reco 2 ref
    THnSparseD *hReco2RefDijet = (THnSparseD*)inFile->Get("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    hReco2RefDijet->SetName("hReco2RefDijet");
    // Retrieve THnSparse for reco 2 ref with selection on ref
    THnSparseD *hRefSelReco2RefDijet = (THnSparseD*)inFile->Get("hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    hRefSelReco2RefDijet->SetName("hRefSelReco2RefDijet");
    // Retrieve THnSparse for gen
    THnSparseD *hGenDijet = (THnSparseD*)inFile->Get("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted");
    hGenDijet->SetName("hGenDijet");

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

    // Refine reco2ref with ref selection

    // Reco dijet pT
    hRefSelReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(0), ptDijetHi.at(0) );
    // Reco leading jet pT
    hRefSelReco2RefDijet->GetAxis(2)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Reco subleading jet pT
    hRefSelReco2RefDijet->GetAxis(4)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );


    // Refine gen distribution

    // Gen dijet pT
    hGenDijet->GetAxis(0)->SetRange( ptDijetLow.at(0), ptDijetHi.at(0) );
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
    TH1D *hRefSelDijetEta = (TH1D*)hRefSelReco2RefDijet->Projection(7); // Test with ref selection

    TH2D *hRef2RecoDijetEta = (TH2D*)hReco2RefDijet->Projection(7, 1);
    // Gen
    TH1D *hGenDijetEta = (TH1D*)hGenDijet->Projection(1);
    TH1D *hUnfoldEta = new TH1D();
    TString name = "TestEta";
    //unfold1D(hRecoDijetEta, hRefDijetEta, hRef2RecoDijetEta, hGenDijetEta, hUnfoldEta, date, name, 4);
    unfold1D(hRecoDijetEta, hRefDijetEta, hRef2RecoDijetEta, hGenDijetEta, hUnfoldEta, date, name, 4);

    //
    // pT dijet
    //

    // // Reco dijet pT (important. It is extremely hardcoded!!!)
    // hReco2RefDijet->GetAxis(0)->SetRange( 1, 121 );

    // // Reco and ref
    // TH1D *hRecoDijetPt = (TH1D*)hReco2RefDijet->Projection(0);
    // TH1D *hRefDijetPt = (TH1D*)hReco2RefDijet->Projection(6);
    // TH2D *hRef2RecoDijetPt = (TH2D*)hReco2RefDijet->Projection(6, 0);
    // // Gen
    // TH1D *hGenDijetPt = (TH1D*)hGenDijet->Projection(0);
    // TH1D *hUnfoldPt = new TH1D();
    // name = "TestPt";
    // unfold1D(hRecoDijetPt, hRefDijetPt, hRef2RecoDijetPt, hGenDijetPt, hUnfoldPt, date, name, 4);
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
    TString name = Form("TestEtaPtHatBinomial%d",ptHat);
    unfold1DPtHat(hRecoDijetPtHatEta, hRefDijetPtHatEta, hRef2RecoDijetEta, 
                  hGenDijetPtHatEta, hUnfoldEta, date, name, 4, ptHat, hRefDijetEta);
}

//________________
void plotDijetDistributions(TFile *inFile, TString date, Int_t jetBranch = 0) {
    
    TString inputFileName( inFile->GetName() );
    TString direction;
    if ( inputFileName.Contains("Pbgoing") ) {
        direction = "Pbgoing";
    }
    else if ( inputFileName.Contains("pgoing") ) {
        direction = "pgoing";
    }
    else {
        direction = "unknownDir";
    }

    TString branchName;
    if ( jetBranch == 0 ) {
        branchName = "akCs4";
    }
    else {
        branchName = "ak4";
    }

    // Retrieve THnSparse for reco 2 ref
    THnSparseD *hReco2RefDijet = (THnSparseD*)inFile->Get("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    hReco2RefDijet->SetName("hReco2RefDijet");
    // Retrieve THnSparse for reco 2 ref with selection on ref
    THnSparseD *hRefSelReco2RefDijet = (THnSparseD*)inFile->Get("hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted");
    hRefSelReco2RefDijet->SetName("hRefSelReco2RefDijet");
    // Retrieve THnSparse for gen
    THnSparseD *hGenDijet = (THnSparseD*)inFile->Get("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted");
    hGenDijet->SetName("hGenDijet");

    TH3D *hGenDijetPtEtaDphi = (TH3D*)inFile->Get("hGenDijetPtEtaDphi");
    TH3D *hRecoDijetPtEtaDphi = (TH3D*)inFile->Get("hRecoDijetPtEtaDphi");
    TH3D *hRefDijetPtEtaDphi = (TH3D*)inFile->Get("hRefDijetPtEtaDphi");
    TH3D *hRefSelDijetPtEtaDphi = (TH3D*)inFile->Get("hRefSelDijetPtEtaDphi");
    TH3D *hRecoDijetEtaRefDijetEtaRecoDijetPt = (TH3D*)inFile->Get("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt");
    TH3D* hRecoDijetPtEtaDphiJetId = (TH3D*)inFile->Get("hRecoDijetPtEtaDphiJetId");

    // Double_t dijetPtVals[dijetPtBins+1] {  40.,  50.,   60.,  70.,  80.,
    //                                        90., 100.,  110., 120., 130.,
    //                                       140., 150.,  160., 180., 200., 
    //                                       240., 300., 1000.};
    // fDijetPtBins{194}, fDijetPtRange{30., 1000.}


    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // fPtBins{150}, fPtRange{20., 1520.}

    // Leading jet pT selection (> 50 GeV/c)
    std::vector<Int_t> ptLeadLow{4};  
    std::vector<Int_t> ptLeadHi{150};

    // SubLeading jet pT selection (> 30 GeV/c)
    std::vector<Int_t> ptSubLeadLow{2};
    std::vector<Int_t> ptSubLeadHi{150};

    //
    // Create canvases to fill
    //

    TCanvas *canv = new TCanvas("canv", "canv", 1200, 900);
    TCanvas *canv2 = new TCanvas("canv2", "canv2", 600, 1200);
    canv2->Divide(1, 2);
    TCanvas *cRefVsReco = new TCanvas("cRefVsReco", "cRefVsReco", 1600, 800);
    cRefVsReco->Divide( 5, ( (ptDijetLow.size() % 5) == 0 ) ? (ptDijetLow.size() / 5) : (ptDijetLow.size() / 5 + 1) );
    TCanvas *cEtaComp = new TCanvas("cEtaComp", "cEtaComp", 1600, 800);
    cEtaComp->Divide( 5, ( (ptDijetLow.size() % 5) == 0 ) ? (ptDijetLow.size() / 5) : (ptDijetLow.size() / 5 + 1) );
    TCanvas *cEtaRatComp = new TCanvas("cEtaRatComp", "cEtaRatComp", 1600, 800);
    cEtaRatComp->Divide( 5, ( (ptDijetLow.size() % 5) == 0 ) ? (ptDijetLow.size() / 5) : (ptDijetLow.size() / 5 + 1) );

    //
    // Apply leading and subleading jet pT selection
    //

    // Reco leading jet pT
    hReco2RefDijet->GetAxis(2)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Reco subleading jet pT
    hReco2RefDijet->GetAxis(4)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    // Reco leading jet pT
    hRefSelReco2RefDijet->GetAxis(2)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Reco subleading jet pT
    hRefSelReco2RefDijet->GetAxis(4)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    // Gen leading jet pT
    hGenDijet->GetAxis(3)->SetRange( ptLeadLow.at(0), ptLeadHi.at(0) );
    // Gen subleading jet pT
    hGenDijet->GetAxis(6)->SetRange( ptSubLeadLow.at(0), ptSubLeadHi.at(0) );

    //
    // Create histograms with projection
    //

    TH1D *hRecoDijetEta[ ptDijetLow.size() ];
    TH1D *hRefDijetEta[ ptDijetLow.size() ];
    TH1D *hRefSelDijetEta[ ptDijetLow.size() ];
    TH2D *hRef2RecoDijetEta[ ptDijetLow.size() ];
    TH1D *hGenDijetEta[ ptDijetLow.size() ];

    TH1D *hRecoDijetEtaPure[ ptDijetLow.size() ];
    TH1D *hRefDijetEtaPure[ ptDijetLow.size() ];
    TH1D *hRefSelDijetEtaPure[ ptDijetLow.size() ];
    TH1D *hGenDijetEtaPure[ ptDijetLow.size() ];

    // jetId
    TH1D *hRecoDijetEtaPureJetId[ ptDijetLow.size() ];

    TH1D *hReco2GenDijetEta[ ptDijetLow.size() ];
    TH1D *hRef2GenDijetEta[ ptDijetLow.size() ];
    TH1D *hRefSel2GenDijetEta[ ptDijetLow.size() ];
    TH1D *hReco2GenDijetEtaJetId[ ptDijetLow.size() ];

    TLegend *leg[ ptDijetLow.size() ];
    TLegend *leg2[ ptDijetLow.size() ];
    TLegend *leg3[ ptDijetLow.size() ];
    TLine   *line[ ptDijetLow.size() ];
    TLine   *line2[ ptDijetLow.size() ];

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};
    Int_t refSelType{2};
    Int_t jetIdType{5};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Loop over dijet pT bins
    for (Int_t i=0; i<ptDijetLow.size(); i++) {

        std::cout << "Working on bin: " << i << std::endl;
        // Set dijet pT selection
        
        // Reco dijet pT
        hReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(i), ptDijetHi.at(i) );
        // Reco dijet pT
        hRefSelReco2RefDijet->GetAxis(0)->SetRange( ptDijetLow.at(i), ptDijetHi.at(i) );
        // Gen dijet pT
        hGenDijet->GetAxis(0)->SetRange( ptDijetLow.at(i), ptDijetHi.at(i) );

        // Make 1D and 2D projections

        // From THnSparseD
        hRecoDijetEta[i] = (TH1D*)hReco2RefDijet->Projection(1);
        hRecoDijetEta[i]->SetNameTitle( Form("hRecoDijetEta_%d", i), ";#eta^{dijet};1/N dN/d#eta^{dijet}");
        hRefDijetEta[i] = (TH1D*)hReco2RefDijet->Projection(7);
        hRefDijetEta[i]->SetNameTitle( Form("hRefDijetEta_%d", i),";#eta^{dijet};1/N dN/d#eta^{dijet}" );
        hRefSelDijetEta[i] = (TH1D*)hRefSelReco2RefDijet->Projection(7);
        hRefSelDijetEta[i]->SetNameTitle( Form("hRefSelDijetEta_%d", i),";#eta^{dijet};1/N dN/d#eta^{dijet}" );
        //hRef2RecoDijetEta[i] = (TH2D*)hReco2RefDijet->Projection(7, 1);
        //hRef2RecoDijetEta[i]->SetNameTitle( Form("hRef2RecoDijetEta_%d", i),";#eta^{dijet}_{reco};#eta^{dijet}_{ref}" );
        hGenDijetEta[i] = (TH1D*)hGenDijet->Projection(1);
        hGenDijetEta[i]->SetNameTitle( Form("hGenDijetEta_%d", i),";#eta^{dijet};1/N dN/d#eta^{dijet}" );

        // From 3D
        hRecoDijetEtaPure[i] = (TH1D*)hRecoDijetPtEtaDphi->ProjectionY( Form("hRecoDijetEtaPure_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        hRecoDijetEtaPure[i]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        hRefDijetEtaPure[i] = (TH1D*)hRefDijetPtEtaDphi->ProjectionY( Form("hRefDijetEtaPure_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        hRefDijetEtaPure[i]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        hRefSelDijetEtaPure[i] = (TH1D*)hRefSelDijetPtEtaDphi->ProjectionY( Form("hRefSelDijetEtaPure_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        hRefSelDijetEtaPure[i]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
        hGenDijetEtaPure[i] = (TH1D*)hGenDijetPtEtaDphi->ProjectionY( Form("hGenDijetEtaPure_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        hGenDijetEtaPure[i]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");

        // JetId
        hRecoDijetEtaPureJetId[i] = (TH1D*)hRecoDijetPtEtaDphiJetId->ProjectionY( Form("hRecoDijetEtaPureJetId_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        hRecoDijetEtaPureJetId[i]->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");

        // 2D from 3D
        hRecoDijetEtaRefDijetEtaRecoDijetPt->GetZaxis()->SetRange(ptDijetLow.at(i), ptDijetHi.at(i));
        hRef2RecoDijetEta[i] = (TH2D*)hRecoDijetEtaRefDijetEtaRecoDijetPt->Project3D("yx");
        hRef2RecoDijetEta[i]->SetNameTitle( Form("hRef2RecoDijetEta_%d", i),";#eta^{dijet}_{reco};#eta^{dijet}_{ref}" );

        // Rescale eta distributions
        rescaleEta( hRecoDijetEta[i] );
        rescaleEta( hRefDijetEta[i] );
        rescaleEta( hRefSelDijetEta[i] );
        rescaleEta( hGenDijetEta[i] );

        rescaleEta( hRecoDijetEtaPure[i] );
        rescaleEta( hRefDijetEtaPure[i] );
        rescaleEta( hRefSelDijetEtaPure[i] );
        rescaleEta( hGenDijetEtaPure[i] );

        // jetId
        rescaleEta( hRecoDijetEtaPureJetId[i] );

        rescaleEta( hRef2RecoDijetEta[i] );

        // Create ratios of reco, ref and refSel to gen

        // hReco2GenDijetEta[i] = (TH1D*)hRecoDijetEta[i]->Clone( Form("hReco2GenDijetEta_%d", i) );
        // hReco2GenDijetEta[i]->Divide(hRecoDijetEta[i], hGenDijetEta[i], 1., 1., "b");
        // hRef2GenDijetEta[i] = (TH1D*)hRefDijetEta[i]->Clone( Form("hRef2GenDijetEta_%d", i) );
        // hRef2GenDijetEta[i]->Divide(hRef2GenDijetEta[i], hGenDijetEta[i], 1., 1., "b");
        // hRefSel2GenDijetEta[i] = (TH1D*)hRefSelDijetEta[i]->Clone( Form("hRefSel2GenDijetEta_%d", i) );
        // hRefSel2GenDijetEta[i]->Divide( hRefSel2GenDijetEta[i], hGenDijetEta[i], 1., 1., "b");

        hReco2GenDijetEta[i] = (TH1D*)hRecoDijetEtaPure[i]->Clone( Form("hReco2GenDijetEta_%d", i) );
        hReco2GenDijetEta[i]->Divide(hReco2GenDijetEta[i], hGenDijetEtaPure[i], 1., 1., "b");
        hReco2GenDijetEta[i]->GetYaxis()->SetTitle("Ratio to Gen");
        hRef2GenDijetEta[i] = (TH1D*)hRefDijetEtaPure[i]->Clone( Form("hRef2GenDijetEta_%d", i) );
        hRef2GenDijetEta[i]->Divide(hRef2GenDijetEta[i], hGenDijetEtaPure[i], 1., 1., "b");
        hRef2GenDijetEta[i]->GetYaxis()->SetTitle("Ratio to Gen");
        hRefSel2GenDijetEta[i] = (TH1D*)hRefSelDijetEtaPure[i]->Clone( Form("hRefSel2GenDijetEta_%d", i) );
        hRefSel2GenDijetEta[i]->Divide( hRefSel2GenDijetEta[i], hGenDijetEtaPure[i], 1., 1., "b");
        hRefSel2GenDijetEta[i]->GetYaxis()->SetTitle("Ratio to Gen");
        
        // JetId
        hReco2GenDijetEtaJetId[i] = (TH1D*)hRecoDijetEtaPureJetId[i]->Clone( Form("hReco2GenDijetEtaJetId_%d", i) );
        hReco2GenDijetEtaJetId[i]->Divide( hReco2GenDijetEtaJetId[i], hGenDijetEtaPure[i], 1., 1., "b" );
        hReco2GenDijetEtaJetId[i]->GetYaxis()->SetTitle("Ratio to Gen");

        // Set style
        set1DStyle( hRecoDijetEta[i], recoType );
        set1DStyle( hRefDijetEta[i], refType );
        set1DStyle( hRefSelDijetEta[i], refSelType );
        set1DStyle( hGenDijetEta[i], genType );

        set1DStyle( hRecoDijetEtaPure[i], recoType );
        set1DStyle( hRefDijetEtaPure[i], refType );
        set1DStyle( hRefSelDijetEtaPure[i], refSelType );
        set1DStyle( hGenDijetEtaPure[i], genType );
        // JetId
        set1DStyle( hRecoDijetEtaPureJetId[i], jetIdType );

        set1DStyle( hReco2GenDijetEta[i], recoType );
        set1DStyle( hRef2GenDijetEta[i], refType );
        set1DStyle( hRefSel2GenDijetEta[i], refSelType );

        // JetId
        set1DStyle( hReco2GenDijetEtaJetId[i], jetIdType );

        set2DStyle( hRef2RecoDijetEta[i] );

        // 2D correlation
        canv->cd();
        setPadStyle();
        hRef2RecoDijetEta[i]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i)-1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        canv->SaveAs( Form("%s/pPb8160_%s_eta_RefVsReco_pT_%d_%d_%s.pdf", date.Data(), direction.Data(),
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep, branchName.Data()) );

        // 1D projections
        canv2->cd(1);
        setPadStyle();
        // hRecoDijetEta[i]->Draw();
        // hRefDijetEta[i]->Draw("same");
        // hRefSelDijetEta[i]->Draw("same");
        // hGenDijetEta[i]->Draw("same");

        hRecoDijetEtaPure[i]->Draw();
        hRefDijetEtaPure[i]->Draw("same");
        hRefSelDijetEtaPure[i]->Draw("same");
        hGenDijetEtaPure[i]->Draw("same");
        hRecoDijetEtaPureJetId[i]->Draw("same");

        leg[i] = new TLegend(0.75, 0.65, 0.83, 0.87);
        leg[i]->SetLineWidth(0);
        leg[i]->AddEntry(hRecoDijetEtaPure[i],Form("Reco trkMax"), "p");
        leg[i]->AddEntry(hRefDijetEtaPure[i],Form("Ref"), "p");
        leg[i]->AddEntry(hRefSelDijetEtaPure[i],Form("RefSel"), "p");
        leg[i]->AddEntry(hGenDijetEtaPure[i],Form("Gen"), "p");
        leg[i]->AddEntry(hRecoDijetEtaPureJetId[i],Form("Reco jetId"), "p");
        leg[i]->SetTextSize(0.05);
        leg[i]->Draw();
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );

        canv2->cd(2);
        setPadStyle();
        hReco2GenDijetEta[i]->Draw();
        hRef2GenDijetEta[i]->Draw("same");
        hRefSel2GenDijetEta[i]->Draw("same");
        hReco2GenDijetEtaJetId[i]->Draw("same");
        hReco2GenDijetEta[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
        leg2[i] = new TLegend(0.4, 0.65, 0.6, 0.87);
        leg2[i]->SetLineWidth(0);
        leg2[i]->AddEntry(hReco2GenDijetEta[i],Form("Reco/Gen"), "p");
        leg2[i]->AddEntry(hRef2GenDijetEta[i],Form("Ref/Gen"), "p");
        leg2[i]->AddEntry(hRefSel2GenDijetEta[i],Form("RefSel/Gen"), "p");
        leg2[i]->AddEntry(hReco2GenDijetEtaJetId[i],Form("Reco/Gen jetId"), "p");
        leg2[i]->SetTextSize(0.05);
        leg2[i]->Draw();
        line[i] = new TLine(hReco2GenDijetEta[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                            hReco2GenDijetEta[i]->GetXaxis()->GetBinUpEdge(hReco2GenDijetEta[i]->GetNbinsX()), 1.);
        line[i]->SetLineColor(kMagenta);
        line[i]->SetLineWidth(3);
        line[i]->SetLineStyle(3);
        line[i]->Draw();
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        canv2->SaveAs( Form("%s/pPb8160_%s_etaDijet_pt_%d_%d_%s.pdf", date.Data(), direction.Data(),
                       ptLow + (ptDijetLow.at(i)-1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep, branchName.Data()) );

        // Fill all 2D distributions
        cRefVsReco->cd(i+1);
        setPadStyle();
        hRef2RecoDijetEta[i]->Draw("colz");
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );

        cEtaComp->cd(i+1);
        setPadStyle();
        hRecoDijetEta[i]->Draw();
        hRefDijetEta[i]->Draw("same");
        hRefSelDijetEta[i]->Draw("same");
        hGenDijetEta[i]->Draw("same");
        hRecoDijetEtaPureJetId[i]->Draw("same");
        leg[i] = new TLegend(0.75, 0.65, 0.83, 0.87);
        leg[i]->SetLineWidth(0);
        leg[i]->AddEntry(hRecoDijetEta[i],Form("Reco"), "p");
        leg[i]->AddEntry(hRefDijetEta[i],Form("Ref"), "p");
        leg[i]->AddEntry(hRefSelDijetEta[i],Form("RefSel"), "p");
        leg[i]->AddEntry(hGenDijetEta[i],Form("Gen"), "p");
        leg[i]->AddEntry(hRecoDijetEtaPureJetId[i],Form("Reco jetId"), "p");
        leg[i]->SetTextSize(0.05);
        leg[i]->Draw();
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );

        cEtaRatComp->cd(i+1);
        setPadStyle();
        hReco2GenDijetEta[i]->Draw();
        hRef2GenDijetEta[i]->Draw("same");
        hRefSel2GenDijetEta[i]->Draw("same");
        hReco2GenDijetEtaJetId[i]->Draw("same");
        hReco2GenDijetEta[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
        leg3[i] = new TLegend(0.4, 0.65, 0.6, 0.87);
        leg3[i]->SetLineWidth(0);
        leg3[i]->AddEntry(hReco2GenDijetEta[i],Form("Reco/Gen"), "p");
        leg3[i]->AddEntry(hRef2GenDijetEta[i],Form("Ref/Gen"), "p");
        leg3[i]->AddEntry(hRefSel2GenDijetEta[i],Form("RefSel/Gen"), "p");
        leg3[i]->AddEntry(hReco2GenDijetEtaJetId[i],Form("Reco/Gen jetId"), "p");
        leg3[i]->SetTextSize(0.05);
        leg3[i]->Draw();
        line2[i] = new TLine(hReco2GenDijetEta[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                            hReco2GenDijetEta[i]->GetXaxis()->GetBinUpEdge(hReco2GenDijetEta[i]->GetNbinsX()), 1.);
        line2[i]->SetLineColor(kMagenta);
        line2[i]->SetLineWidth(3);
        line2[i]->SetLineStyle(3);
        line2[i]->Draw();
        t.DrawLatexNDC(0.35, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
    } // for (Int_t i=0; i<dijetPtLow.size(); i++)

    cRefVsReco->SaveAs( Form("%s/pPb8160_%s_eta_RefVsReco_all_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );
    cEtaComp->SaveAs( Form("%s/pPb8160_%s_eta_all_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );
    cEtaRatComp->SaveAs( Form("%s/pPb8160_%s_etaRatios_all_%s.pdf", date.Data(), direction.Data(), branchName.Data()) );
}

//________________
void unfoldDistributions() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing_ak4_etaSel.root";
    TFile *inFile = TFile::Open(inFileName);

    TDatime dt;
    TString date { Form( "%d",dt.GetDate() ) };
    TString inputFileName(inFileName);

    if ( directoryExists( date.Data() ) ) {
        //std::cout << "Directory exists." << std::endl;
    } 
    else {
        createDirectory( date.Data() );
    }

    Int_t branchId{0};
    if ( inputFileName.Contains("akCs4") ) {
        branchId = {0};
    }
    else {
        branchId = {1};
    }

    Int_t ptHat = 50; // 50, 120, 370
    const Char_t *inFileNameWithPtHat = Form("../build/oEmbedding_pPb8160_Pbgoing_%d.root", ptHat);
    TFile *inFileWithPtHat = TFile::Open(inFileNameWithPtHat);

    // Run 1D unfolding
    //performUnfolding(inFile, date);

    // Run 1D unfolding with semi-automated regime
    //unfoldDijetEta1D(inFile, date);

    // Plot various dijet distributions
    plotDijetDistributions(inFile, date, branchId);

    // Run 1D unfolding for ptHat distributions with the response matrix
    // for total ptHat
    //unfoldDifferentPtHat(inFile, inFileWithPtHat, ptHat, date);
}