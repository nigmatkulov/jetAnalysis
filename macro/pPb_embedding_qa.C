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
#include "TLine.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TSystemDirectory.h"
#include <sys/stat.h>

// C++ headers
#include <iostream>

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
    //gPad->SetLogz(1);

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
    } // for (Int_t i{1}; i<=fPtBins; i++)

    for (Int_t i{1}; i<=fEtaBins; i++) {
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
    } // for (Int_t i{1}; i<=fEtaBins; i++)

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
    // Retrieve reco
    TH2D* hRecoInclusiveAllJetPtVsEta = (TH2D*)inFile->Get("hRecoInclusiveAllJetPtVsEta");
    // Retrieve reco that matched gen
    TH2D* hRecoInclusiveMatchedJetPtVsEta = (TH2D*)inFile->Get("hRecoInclusiveMatchedJetPtVsEta");
    // Create a histogram to calculate fakes
    TH2D* hNumberOfFakes = (TH2D*)hRecoInclusiveAllJetPtVsEta->Clone("hNumberOfFakes");
    hNumberOfFakes->SetTitle("Number of fakes (reco - recoMatched);#eta;p_{T} (GeV/c)");
    hNumberOfFakes->Add(hRecoInclusiveMatchedJetPtVsEta, -1.);

    // z axis scaling
    Double_t zAxisRange[2] {0., hRecoInclusiveAllJetPtVsEta->GetMaximum()};
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Create canvas to plot 2D distributions
    TCanvas *cFakes2D = new TCanvas("cFakes2D","cFakes2D", 1300, 400);
    cFakes2D->Divide(3, 1);

    // Plot 2D reco
    cFakes2D->cd(1);
    setPadStyle();
    hRecoInclusiveAllJetPtVsEta->Draw("colz");
    t.DrawLatexNDC(0.46, 0.96, Form("Reco") );

    // Plot 2D reco that matched gen
    cFakes2D->cd(2);
    setPadStyle();
    hRecoInclusiveMatchedJetPtVsEta->Draw("colz");
    hRecoInclusiveMatchedJetPtVsEta->GetZaxis()->SetRangeUser(zAxisRange[0], zAxisRange[1]);
    t.DrawLatexNDC(0.4, 0.96, Form("Reco matched") );

    // Plot 2D fakes
    cFakes2D->cd(3);
    setPadStyle();
    hNumberOfFakes->Draw("colz");
    hNumberOfFakes->GetZaxis()->SetRangeUser(zAxisRange[0], zAxisRange[1]);
    t.DrawLatexNDC(0.46, 0.96, Form("Fakes") );

    //
    // Plot projections of fakes as a function of eta and pT
    //

    TH2D *hRecoMatchedJetFrac = (TH2D*)hRecoInclusiveMatchedJetPtVsEta->Clone("hRecoMatchedJetFrac");
    TH2D *hRecoFakeJetFrac = (TH2D*)hNumberOfFakes->Clone("hRecoFakeJetFrac");
    hRecoMatchedJetFrac->Divide(hRecoMatchedJetFrac, hRecoInclusiveAllJetPtVsEta, 1., 1., "b");
    hRecoFakeJetFrac->Divide(hRecoFakeJetFrac, hRecoInclusiveAllJetPtVsEta, 1., 1., "b");

    Int_t fPtBins = hNumberOfFakes->GetNbinsY();
    Double_t fPtRange[2] {hNumberOfFakes->GetYaxis()->GetBinLowEdge(1), 
                          hNumberOfFakes->GetYaxis()->GetBinUpEdge(fPtBins)};
    Int_t fEtaBins = hNumberOfFakes->GetNbinsX();
    Double_t fEtaRange[2] {hNumberOfFakes->GetXaxis()->GetBinLowEdge(1), 
                           hNumberOfFakes->GetXaxis()->GetBinUpEdge(fEtaBins)};
    Double_t ptStep = (fPtRange[1]-fPtRange[0]) / fPtBins;
    Double_t etaStep = (fEtaRange[1]-fEtaRange[0]) / fEtaBins;

    TH1D *hEtaFakes[fPtBins];
    TH1D *hPtFakes[fEtaBins];
    TH1D *hEtaMatched[fPtBins];
    TH1D *hPtMatched[fEtaBins];

    TCanvas *cEtaFakes = new TCanvas("cEtaFakes", "cEtaFakes", 1600, 800);
    cEtaFakes->Divide(5, ( (fPtBins % 5) == 0 ) ? (fPtBins / 5) : (fPtBins / 5 + 1) );

    TCanvas *cPtFakes = new TCanvas("cPtFakes", "cPtFakes", 1600, 800);
    cPtFakes->Divide(10, 5);

    // Make projections bin-by-bin on eta
    for (Int_t i{1}; i<=fPtBins; i++) {
        // Project on eta
        hEtaFakes[i] = (TH1D*)hRecoFakeJetFrac->ProjectionX(Form("hEtaFakes_%d", i), i, i);
        hEtaFakes[i]->SetNameTitle(Form("hEtaFakes_%d", i), ";#eta;Efficiency");
        set1DStyle(hEtaFakes[i], 2);
        hEtaFakes[i]->SetMarkerSize(0.7);

        hEtaMatched[i] = (TH1D*)hRecoMatchedJetFrac->ProjectionX(Form("hEtaMatched_%d", i), i, i);
        hEtaMatched[i]->SetNameTitle(Form("hEtaMatched_%d", i), ";#eta;Efficiency");
        set1DStyle(hEtaMatched[i], 1);
        hEtaMatched[i]->SetMarkerSize(0.7);

        cEtaFakes->cd(i);
        setPadStyle();
        hEtaFakes[i]->Draw();
        hEtaMatched[i]->Draw("same");
        hEtaFakes[i]->GetYaxis()->SetRangeUser(0., 1.01);
        t.DrawLatexNDC(0.35, 0.95, Form("%4.1f < p_{T} (GeV/c) < %4.1f", 
                       fPtRange[0] + (i-1) * ptStep, fPtRange[0] + i * ptStep) );
    } // for (Int_t i{1}; i<=fPtBins; i++)

    // Make projections bin-by-bin on pT
    for (Int_t i{1}; i<=fEtaBins; i++) {
        // Project on pT
        hPtFakes[i] = (TH1D*)hRecoFakeJetFrac->ProjectionY(Form("hPtFakes_%d", i), i, i);
        hPtFakes[i]->SetNameTitle(Form("hPtFakes_%d", i), ";p_{T} (GeV/c);Efficiency");
        set1DStyle(hPtFakes[i], 2);
        hPtFakes[i]->SetMarkerSize(0.7);

        hPtMatched[i] = (TH1D*)hRecoMatchedJetFrac->ProjectionY(Form("hPtMatched_%d", i), i, i);
        hPtMatched[i]->SetNameTitle(Form("hPtMatched_%d", i), ";p_{T} (GeV/c);Efficiency");
        set1DStyle(hPtMatched[i], 1);
        hPtMatched[i]->SetMarkerSize(0.7);

        cPtFakes->cd(i);
        setPadStyle();
        hPtFakes[i]->Draw();
        hPtMatched[i]->Draw("same");
        hPtFakes[i]->GetYaxis()->SetRangeUser(0., 1.01);
        t.DrawLatexNDC(0.35, 0.95, Form("%2.1f < #eta < %2.1f", 
                       fEtaRange[0] + (i-1) * etaStep, fEtaRange[0] + i * etaStep) );
    }

    cEtaFakes->SaveAs( Form("%s/pPb8160_eta_fakes_projections.pdf", date.Data()) );
    cPtFakes->SaveAs( Form("%s/pPb8160_pt_fakes_projections.pdf", date.Data()) );
}

//________________
void pPb_embedding_qa(const Char_t *inFileName = "../build/oEmbedding_pPb8160_Pbgoing_akCs4.root") {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TDatime dt;
    TString date { Form( "%d",dt.GetDate() ) };

    if ( directoryExists( date.Data() ) ) {
        //std::cout << "Directory exists." << std::endl;
    } 
    else {
        createDirectory( date.Data() );
    }


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

    // Plot reco, reco with matching and calculate fakes
    plotRecoAndFakes(inFile, date);

    // Plot correlation between ref and reco dijet eta
    //plotEtaDijetCorrelation(inFile, date);
}