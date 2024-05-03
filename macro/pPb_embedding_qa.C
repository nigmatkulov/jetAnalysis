// ROOT headers
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
#include "TEfficiency.h"
#include "TLine.h"

// C++ headers
#include <iostream>

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
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.15);
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.2;
    Int_t lineWidth = 2;
    Int_t color = 2;
    if ( !h ) {
        std::cout << "[ERROR] set1DStyle - No histogram passed\n";
    }
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

    h->SetLineColor( color );
    h->SetLineWidth( lineWidth );
    
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
    h->GetYaxis()->SetTitleOffset(1.0);

    if (doRenorm) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void plotEfficiency(TFile *inFile, TString date) {

    // Rebinning
    Int_t rebinX{1}, rebinY{1};

    //
    // Plot ptHat distribution
    //
    TH1D* hPtHat = (TH1D*)inFile->Get("hPtHatWeighted");
    hPtHat->SetName("hPtHat");

    TCanvas *cPtHat = new TCanvas("cPtHat","cPtHat", 1200, 800);
    setPadStyle();
    set1DStyle(hPtHat, 2);
    hPtHat->Draw();
    gPad->SetLogy(1);
    cPtHat->SaveAs(Form("%s/pPb8160_ptHat.pdf", date.Data()));


    //
    // Efficiency over the acceptance
    //

    // Read gen jet acceptance
    TH2D *hGenPtVsEta = (TH2D*)inFile->Get("hGenInclusiveJetPtEta");
    hGenPtVsEta->SetName("hGenPtVsEta");
    // Read reco that matched gen acceptance
    TH2D *hRefPtVsEta = (TH2D*)inFile->Get("hRefInclusiveJetPtEta");
    //TH2D *hRefPtVsEta = (TH2D*)inFile->Get("hRecoMatchedPtEta");
    TH2D *hRecoPtVsEta = (TH2D*)inFile->Get("hRecoMatchedPtEta");
    hRefPtVsEta->SetName("hRefPtVsEta");

    // Perform rebinning
    hGenPtVsEta->RebinX( rebinX );
    hGenPtVsEta->RebinY( rebinY );

    hRefPtVsEta->RebinX( rebinX );
    hRefPtVsEta->RebinY( rebinY );

    // Create plot for efficiency
    TH2D *hEfficiency = new TH2D("hEfficiency","Inclusive ref(reco matched)/gen;#eta;p_{T} (GeV/c)",
                                 hGenPtVsEta->GetNbinsX(), 
                                 hGenPtVsEta->GetXaxis()->GetBinLowEdge(1),
                                 hGenPtVsEta->GetXaxis()->GetBinUpEdge( hGenPtVsEta->GetNbinsX() ),
                                 hGenPtVsEta->GetNbinsY(), 
                                 hGenPtVsEta->GetYaxis()->GetBinLowEdge(1),
                                 hGenPtVsEta->GetYaxis()->GetBinUpEdge( hGenPtVsEta->GetNbinsX() ) );
    hEfficiency->Sumw2();
    hEfficiency->Divide(hRefPtVsEta, hGenPtVsEta, 1., 1., "b");
    set2DStyle( hEfficiency );

    // Plot efficiency
    TCanvas *cEfficiency = new TCanvas("cEfficiency","cEfficiency", 1200, 800);
    setPadStyle();
    hEfficiency->Draw("colz");
    gPad->SetLogz(1);

    cEfficiency->SaveAs(Form("%s/pPb8160_pt_vs_eta_efficiency.pdf", date.Data()));


    Int_t fPtBins{50};
    Double_t fPtRange[2] {20., 520.};
    Int_t fEtaBins{50}; 
    Double_t fEtaRange[2] {-5.0, 5.0};
    Double_t ptStep = (fPtRange[1]-fPtRange[0]) / fPtBins;
    Double_t etaStep = (fEtaRange[1]-fEtaRange[0]) / fEtaBins;

    TH1D *hEtaEfficiency[fPtBins];
    TH1D *hPtEfficiency[fEtaBins];

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    TCanvas *cEtaEfficiency = new TCanvas("cEtaEfficiency", "cEtaEfficiency", 1600, 800);
    cEtaEfficiency->Divide(10, 5);

    TCanvas *cPtEfficiency = new TCanvas("cPtEfficiency", "cPtEfficiency", 1600, 800);
    cPtEfficiency->Divide(10, 5);

    // Make projections bin-by-bin
    for (Int_t i{1}; i<=fPtBins; i++) {
        // Project on eta
        hEtaEfficiency[i] = (TH1D*)hEfficiency->ProjectionX(Form("hEtaEfficiency_%d", i), i, i);
        hEtaEfficiency[i]->SetNameTitle(Form("hEtaEfficiency_%d", i), ";#eta;Efficiency");
        set1DStyle(hEtaEfficiency[i], 2);
        hEtaEfficiency[i]->SetMarkerSize(0.7);
        cEtaEfficiency->cd(i);
        setPadStyle();
        hEtaEfficiency[i]->Draw();
        if (i==1) {
            hEtaEfficiency[i]->GetYaxis()->SetRangeUser(0.35, 0.85);
        }
        else {
            hEtaEfficiency[i]->GetYaxis()->SetRangeUser(0.9, 1.05);
        }
        t.DrawLatexNDC(0.2, 0.9, Form("%4.1f < p_{T} (GeV/c) < %4.1f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep) );

        // Project on pT
        hPtEfficiency[i] = (TH1D*)hEfficiency->ProjectionY(Form("hPtEfficiency_%d", i), i, i);
        hPtEfficiency[i]->SetNameTitle(Form("hPtEfficiency_%d", i), ";p_{T} (GeV/c);Efficiency");
        set1DStyle(hPtEfficiency[i], 2);
        hPtEfficiency[i]->SetMarkerSize(0.7);
        cPtEfficiency->cd(i);
        setPadStyle();
        hPtEfficiency[i]->Draw();
        hPtEfficiency[i]->GetYaxis()->SetRangeUser(0.9, 1.02);
        t.DrawLatexNDC(0.2, 0.9, Form("%2.1f < #eta < %2.1f", 
                       fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep) );
    } // for (Int_t i{1}; i<=fPtBins; i++)

    cEtaEfficiency->SaveAs( Form("%s/pPb8160_eta_efficiency_projections.pdf", date.Data()) );
    cPtEfficiency->SaveAs( Form("%s/pPb8160_pt_efficiency_projections.pdf", date.Data()) );


    // // Reco inclusive jet pT
    // TH1D *hReco = (TH1D*)inFile->Get("hRecoInclusiveJetPt");
    // // Reco inclusive jet pT
    // TH1D *hRef = (TH1D*)inFile->Get("hRefInclusiveJetPt");
    // // Reco inclusive jet pT
    // TH1D *hGen = (TH1D*)inFile->Get("hGenInclusiveJetPt");

    // TH1D *hGenPtProj = (TH1D*)hGenPtVsEta->ProjectionY();
    // TH1D *hRefPtProj = (TH1D*)hRefPtVsEta->ProjectionY();
    // TH1D *hRecoPtProj = (TH1D*)hRecoPtVsEta->ProjectionY();

    // TCanvas *cInclPt = new TCanvas("cInclPt", "cInclPt", 1200, 800);
    // cInclPt->Divide(3, 1);

    // cInclPt->cd(1);
    // setPadStyle();
    // set1DStyle(hGen, 0);
    // set1DStyle(hGenPtProj, 1);
    // hGen->Draw();
    // set1DStyle(hGen, 0);
    // hGenPtProj->Draw("same");
    // gPad->SetLogy(1);

    // cInclPt->cd(2);
    // setPadStyle();
    // set1DStyle(hRef, 0);
    // set1DStyle(hRefPtProj, 1);
    // hRef->Draw();
    // hRefPtProj->Draw("same");
    // gPad->SetLogy(1);

    // cInclPt->cd(3);
    // setPadStyle();
    // set1DStyle(hReco, 0);
    // set1DStyle(hRecoPtProj, 1);
    // hReco->Draw();
    // hRecoPtProj->Draw("same");
    // gPad->SetLogy(1);

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

    cDijetEtaRefVsReco->SaveAs(Form("%s/pPb8160_ref_vs_reco_responce.pdf", date.Data()));
}


//________________
void plotDijetDistributions(TFile *inFile, TString date) {

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};
    Int_t refSelType{2};
    Bool_t doRenorm{kTRUE};

    TH1D *hRecoDijetEta = (TH1D*)inFile->Get("hRecoDijetEta");
    if ( !hRecoDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No recoDijetEta found\n";
    }
    TH1D *hRefDijetEta = (TH1D*)inFile->Get("hRefDijetEta");
    if ( !hRefDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No refDijetEta found\n";
    }
    TH1D *hGenDijetEta = (TH1D*)inFile->Get("hGenDijetEta");
    if ( !hGenDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No genDijetEta found\n";
    }
    TH1D *hRefSelDijetEta = (TH1D*)inFile->Get("hRefSelDijetEta");
    if ( !hRefSelDijetEta ) {
        std::cout << "[WARNING] plotDijetDistributions - No refSelDijetEta found\n";
    }

    set1DStyle(hRecoDijetEta,   recoType, doRenorm);
    set1DStyle(hRefDijetEta,    refType, doRenorm);
    set1DStyle(hGenDijetEta,    genType, doRenorm);
    set1DStyle(hRefSelDijetEta, refSelType, doRenorm);

    TH1D *hReco2Gen = new TH1D("hReco2Gen", "hReco2Gen;#eta_{dijet};#frac{reco}{gen}",
                                hRecoDijetEta->GetNbinsX(), 
                                hRecoDijetEta->GetXaxis()->GetBinLowEdge(1),
                                hRecoDijetEta->GetXaxis()->GetBinUpEdge( hRecoDijetEta->GetNbinsX() ) );
    hReco2Gen->Sumw2();
    set1DStyle(hReco2Gen, recoType);
    TH1D *hRef2Gen = new TH1D("hRef2Gen", "hRef2Gen;#eta_{dijet};#frac{ref}{gen}",
                                hRefDijetEta->GetNbinsX(), 
                                hRefDijetEta->GetXaxis()->GetBinLowEdge(1),
                                hRefDijetEta->GetXaxis()->GetBinUpEdge( hRefDijetEta->GetNbinsX() ) );
    hRef2Gen->Sumw2();
    set1DStyle(hRef2Gen, refType);
    TH1D *hRefSel2Gen = new TH1D("hRefSel2Gen", "hRefSel2Gen;#eta_{dijet};#frac{ref sel}{gen}",
                                hRefDijetEta->GetNbinsX(), 
                                hRefDijetEta->GetXaxis()->GetBinLowEdge(1),
                                hRefDijetEta->GetXaxis()->GetBinUpEdge( hRefDijetEta->GetNbinsX() ) );
    hRefSel2Gen->Sumw2();
    set1DStyle(hRefSel2Gen, refSelType);

    hReco2Gen->Divide(hRecoDijetEta, hGenDijetEta, 1., 1., "b");
    hRef2Gen->Divide(hRefDijetEta, hGenDijetEta, 1., 1., "b");
    hRefSel2Gen->Divide(hRefSelDijetEta, hGenDijetEta, 1., 1., "b");

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TCanvas *c1DEta = new TCanvas("c1DEta", "c1DEta", 1000, 1000);
    c1DEta->Divide(1, 2);

    c1DEta->cd(1);
    setPadStyle();
    hGenDijetEta->Draw();
    hRecoDijetEta->Draw("same");
    hRefDijetEta->Draw("same");
    hRefSelDijetEta->Draw("same");
    hGenDijetEta->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    TLegend *leg = new TLegend(0.65, 0.68, 0.8, 0.92);
    leg->SetLineWidth(0);
    leg->AddEntry(hGenDijetEta,Form("Gen"), "p");
    leg->AddEntry(hRefDijetEta,Form("Ref"), "p");
    leg->AddEntry(hRecoDijetEta,Form("Reco"), "p");
    leg->AddEntry(hRefSelDijetEta,Form("Ref Sel"), "p");
    leg->SetTextSize(0.06);
    leg->Draw();

    c1DEta->cd(2);
    setPadStyle();
    hReco2Gen->Draw();
    hRef2Gen->Draw("same");
    hRefSel2Gen->Draw("same");
    TLegend *leg2 = new TLegend(0.65, 0.68, 0.8, 0.92);
    leg2->SetLineWidth(0);
    leg2->AddEntry(hRef2Gen,Form("Ref/Gen"), "p");
    leg2->AddEntry(hReco2Gen,Form("Reco/Gen"), "p");
    leg2->AddEntry(hRefSel2Gen,Form("RefSel/Gen"), "p");
    leg2->SetTextSize(0.06);
    leg2->Draw();
    hReco2Gen->GetYaxis()->SetRangeUser(0.8, 1.2);
    hReco2Gen->GetYaxis()->SetTitle("Ratio to Gen");
    TLine *l = new TLine(hReco2Gen->GetXaxis()->GetBinLowEdge(1), 1., 
                         hReco2Gen->GetXaxis()->GetBinUpEdge(hReco2Gen->GetNbinsX()), 1.);
    l->SetLineColor(kMagenta);
    l->SetLineWidth(3);
    l->SetLineStyle(3);
    l->Draw();

    c1DEta->SaveAs( Form("%s/pPb8160_eta_dijet_comparison.pdf", date.Data()) );
}

//________________
void compareInclusiveJetPtSpectra(TFile *inFile, TString date) {

    Int_t recoType{0};
    Int_t refType{1};
    Int_t genType{3};

    // Reco inclusive jet pT
    TH1D *hReco = (TH1D*)inFile->Get("hRecoInclusiveJetPt");
    // Reco inclusive jet pT
    TH1D *hRef = (TH1D*)inFile->Get("hRefInclusiveJetPt");
    // Reco inclusive jet pT
    TH1D *hGen = (TH1D*)inFile->Get("hGenInclusiveJetPt");

    // Set style
    set1DStyle(hReco, recoType, 0);
    set1DStyle(hRef, refType, 0);
    set1DStyle(hGen, genType, 0);

    // Create histograms to fill with ratios
    TH1D* hReco2GenRatio = new TH1D("hReco2GenRatio","reco / gen ratio;p_{T} (GeV/c);#frac{reco}{gen}",
                                    hGen->GetNbinsX(), 
                                    hGen->GetXaxis()->GetBinLowEdge(1),
                                    hGen->GetXaxis()->GetBinUpEdge( hGen->GetNbinsX() ) );
    hReco2GenRatio->Sumw2();
    set1DStyle(hReco2GenRatio, recoType);
    TH1D* hRef2GenRatio = new TH1D("hRef2GenRatio","ref / gen ratio;p_{T} (GeV/c);#frac{ref}{gen}",
                                    hGen->GetNbinsX(), 
                                    hGen->GetXaxis()->GetBinLowEdge(1),
                                    hGen->GetXaxis()->GetBinUpEdge( hGen->GetNbinsX() ) );
    hRef2GenRatio->Sumw2();
    set1DStyle(hRef2GenRatio, genType);


    // Make ratios
    hReco2GenRatio->Divide(hReco, hGen, 1., 1., "b");
    hRef2GenRatio->Divide(hRef, hGen, 1., 1., "b");

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TCanvas *cPtSpectra = new TCanvas("cPtSpectra", "cPtSpectra", 600, 1200);

    cPtSpectra->Divide(1, 2);

    cPtSpectra->cd(1);
    setPadStyle();
    hReco->Draw();
    hRef->Draw("same");
    hGen->Draw("same");
    hReco->GetXaxis()->SetTitle("Inclusive jet p_{T} (GeV/c)");
    gPad->SetLogy(1);

    TLegend *leg = new TLegend(0.7, 0.68, 0.85, 0.92);
    leg->SetLineWidth(0);
    leg->AddEntry(hReco,Form("Reco"), "p");
    leg->AddEntry(hRef,Form("Ref"), "p");
    leg->AddEntry(hGen,Form("Gen"), "p");
    leg->SetTextSize(0.06);
    leg->Draw();

    cPtSpectra->cd(2);
    setPadStyle();
    hReco2GenRatio->Draw();
    hRef2GenRatio->Draw("same");
    //hGen->Draw("same");
    hReco2GenRatio->GetYaxis()->SetRangeUser(0.87, 1.07);
    hReco2GenRatio->GetYaxis()->SetTitle("Ratio to gen");
    TLegend *leg2 = new TLegend(0.65, 0.75, 0.82, 0.92);
    leg2->SetLineWidth(0);
    leg2->AddEntry(hReco2GenRatio,Form("Reco/Gen"), "p");
    leg2->AddEntry(hRef2GenRatio,Form("Ref/Gen"), "p");
    leg2->SetTextSize(0.06);
    leg2->Draw();
}

//________________
void plotRecoAndFakes(TFile *inFile, TString date) {
    TH2D* hRecoInclusiveAllJetPtVsEta = (TH2D*)inFile->Get("hRecoInclusiveAllJetPtVsEta");
    TH2D* hRecoInclusiveMatchedJetPtVsEta = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEta");
    TH2D* hNumberOfFakes = new TH2D("hNumberOfFakes","Number of fakes (reco - recoMatched);#eta;p_{T} (GeV/c)",
                                    hRecoInclusiveAllJetPtVsEta->GetNbinsX(), hRecoInclusiveAllJetPtVsEta->GetXaxis()->GetXbins()->GetArray(),
                                    hRecoInclusiveAllJetPtVsEta->GetNbinsY(), hRecoInclusiveAllJetPtVsEta->GetYaxis()->GetXbins()->GetArray());
    hNumberOfFakes->Sumw2();
    hNumberOfFakes->Add(hRecoInclusiveAllJetPtVsEta, hRecoInclusiveMatchedJetPtVsEta, 1., -1.);

    TCanvas *cFakes2D = new TCanvas("cFakes2D","cFakes2D", 1300, 400);
    cFakes2D->Divide(3, 1);

    cFakes2D->cd(1);
    setPadStyle();
    hRecoInclusiveAllJetPtVsEta->Draw("colz");

    cFakes2D->cd(2);
    setPadStyle();
    hRecoInclusiveMatchedJetPtVsEta->Draw("colz");

    cFakes2D->cd(3);
    setPadStyle();
    cFakes2D->Draw("colz");
}

//________________
void pPb_embedding_qa(const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing.root") {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TString date {"20240503"};
    TFile *inFile = TFile::Open(inFileName);
    if ( !inFile ) {
        std::cout << "[ERROR] Input file does not exist\n";
        return;
    }
    if ( inFile->IsZombie() ) {
        std::cout << "[ERROR] Beware of zombies!!!!\n";
        return;
    }

    // Compare inclusive reco, ref and gen transverse momentum spectra
    //compareInclusiveJetPtSpectra(inFile, date);

    // Plot jet reconstruction efficiency as a function of acceptance (pT vs eta)
    //plotEfficiency(inFile, date);

    // Plot dijet distributions
    //plotDijetDistributions(inFile, date);

    // Plot various dijet distributions
    //plotDijetDistributions(inFile, date);

    // Plot reco, reco with matching and calculate fakes
    plotRecoAndFakes(inFile, date);

    // Plot correlation between ref and reco dijet eta
    //plotEtaDijetCorrelation(inFile, date);
}