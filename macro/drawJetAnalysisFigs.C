// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TGraph.h"
#include "TMultiGraph.h"

// C++ headers
#include <iostream>
#include <vector>

using namespace std;

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
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.15);
}

//________________
void calculateCumulativeIntegral(TGraph* gr, TH1* h, Int_t type=0) {

    Double_t integral = h->Integral();
    Double_t cumulativeInegtral{0};
    for (Int_t i{1}; i<=h->GetNbinsX(); i++) {
        Double_t x = h->GetXaxis()->GetBinCenter( i );
        cumulativeInegtral += h->GetBinContent( i );
        gr->AddPoint(x, cumulativeInegtral / integral);
    }

    // type: 0 - reco, 1 - gen, 2 - ref
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.2;
    Int_t lineWidth{2}; 
    Int_t color{2};
    if (type == 0) {

    }
    else if (type == 1) { // gen
        color = 4;
        markerStyle = 24;
    }
    else if (type == 2) { // ref
        color = 8;
        markerStyle = 33;
        markerSize = 1.5;
    }
    else if (type == 3) {
        color = 6;
        markerStyle = 35;
        markerSize = 1.3;
    }
    else if ( type == 4 ) {
        color = 7;
        markerStyle = 37;
        markerSize = 1.3;  
    }
    else if ( type == 5 ) {
        color = 11;
        markerStyle = 40;
        markerSize = 1.3;  
    }
    else if ( type == 6 ) {
        color = 15;
        markerStyle = 41;
        markerSize = 1.3;  
    }
    else {
        color = 9;
        markerStyle = 38;
        markerSize = 1.3;
    }

    gr->SetLineWidth( lineWidth );
    gr->SetLineColor( color );
    
    gr->SetMarkerStyle( markerStyle );
    gr->SetMarkerColor( color );
    gr->SetMarkerSize( markerSize );
}

//________________
void setJetStyle(TH1 *h, Int_t type = 0) {
    // type: 0 - reco, 1 - gen, 2 - ref
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.2;
    Int_t lineWidth{2}; 
    Int_t color{2};
    if (type == 0) {

    }
    else if (type == 1) { // gen
        color = 4;
        markerStyle = 24;
    }
    else if (type == 2) { // ref
        color = 8;
        markerStyle = 33;
        markerSize = 1.5;
    }
    else if (type == 3) {
        color = 6;
        markerStyle = 35;
        markerSize = 1.3;
    }
    else if ( type == 4 ) {
        color = 7;
        markerStyle = 37;
        markerSize = 1.3;  
    }
    else if ( type == 5 ) {
        color = 11;
        markerStyle = 40;
        markerSize = 1.3;  
    }
    else if ( type == 6 ) {
        color = 15;
        markerStyle = 41;
        markerSize = 1.3;  
    }
    else {
        color = 9;
        markerStyle = 38;
        markerSize = 1.3;
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
void updateVectorLimits(std::vector<Int_t> &vec) {
    Int_t firstHistBin{1};
    for (Int_t i{0}; i<vec.size(); i++) {
        vec.at(i) += firstHistBin;
    }
}

//________________
void drawEventQuantities(TFile *inFile) {
    // Event histograms
    TH1D *hVz = (TH1D*)inFile->Get("hVz");
    TH1D *hVzWeighted = (TH1D*)inFile->Get("hVzWeighted");
    TH1D *hHiBin = (TH1D*)inFile->Get("hHiBin");
    TH1D *hHiBinWeighted = (TH1D*)inFile->Get("hHiBinWeighted");
    TH1D *hPtHat = (TH1D*)inFile->Get("hPtHat");
    TH1D *hPtHatWeighted = (TH1D*)inFile->Get("hPtHatWeighted");
    TH1D *hPtHatWeight = (TH1D*)inFile->Get("hPtHatWeight");

    set1DStyle(hPtHatWeight, 0);
    set1DStyle(hVz, 0);
    set1DStyle(hHiBin, 0);
    set1DStyle(hPtHat, 0);

    set1DStyle(hVzWeighted, 1);
    set1DStyle(hHiBinWeighted, 1);
    set1DStyle(hPtHatWeighted, 1); 

    // Event information
    TCanvas *canvEvent = new TCanvas("canvEvent","canvEvent", 1200, 800);
    canvEvent->Divide(2, 2);

    canvEvent->cd(1);
    setPadStyle();
    hPtHatWeight->Draw();
    gPad->SetLogy(1);

    canvEvent->cd(2);
    setPadStyle();
    hVz->Draw();
    hVzWeighted->Draw("same");
    gPad->SetLogy(1); 

    canvEvent->cd(3);
    setPadStyle();
    hHiBin->Draw();
    hHiBinWeighted->Draw("same");
    gPad->SetLogy(1);

    canvEvent->cd(4);
    setPadStyle();
    hPtHat->Draw();
    hPtHatWeighted->Draw("same");
    gPad->SetLogy(1);

    canvEvent->Draw();
}

//________________
void drawNumberOfJets(TFile *inFile) {

    const Int_t nCuts = 5;  // pT > 0, 20, 50, 80, 120 GeV/c
    Int_t cutVals[5] = {0, 20, 50, 80, 120};
    TH1D* hNRecoJets[nCuts];
    TH1D* hNGenJets[nCuts];
    TH1D* hNRefJets[nCuts];
    for (Int_t i{0}; i<nCuts; i++) {
        hNRecoJets[i] = (TH1D*)inFile->Get( Form("hNRecoJets_%d",i) );
        setJetStyle(hNRecoJets[i], 0);

        hNGenJets[i] = (TH1D*)inFile->Get( Form("hNGenJets_%d",i) );
        setJetStyle(hNGenJets[i], 1);

        hNRefJets[i] = (TH1D*)inFile->Get( Form("hNRefJets_%d",i) );
        setJetStyle(hNRefJets[i], 2);
    }

    TCanvas *nJetsCanv[5];
    for (Int_t i{0}; i<nCuts; i++) {
        nJetsCanv[i] = new TCanvas( Form("nJetsCanv_%d",i), 
                                    Form("nJetsCanv_%d",i), 
                                    800, 800);
        setPadStyle();
        hNRecoJets[i]->GetXaxis()->SetTitle("Number of jets");
        hNRecoJets[i]->GetYaxis()->SetTitle("Entries");
        hNRecoJets[i]->Draw();
        hNGenJets[i]->Draw("same");
        hNRefJets[i]->Draw("same");
        gPad->SetLogy(1);

        TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->SetHeader("Jet types", "C");
        leg->AddEntry(hNRecoJets[i],Form("Reco jets with p_{T}>%d GeV/c",cutVals[i]), "p");
        leg->AddEntry(hNGenJets[i],Form("Gen jets with p_{T}>%d GeV/c",cutVals[i]), "p");
        leg->AddEntry(hNRefJets[i],Form("Ref jets with p_{T}>%d GeV/c",cutVals[i]), "p");

        leg->Draw();
    }
}

//________________
void drawReco2GenSpectraPtHatComparison(TFile *inFile) {

    // TODO:
    // Should be inclusive vs. gen!!
    THnSparseD *hRecoJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentWeighted");
    THnSparseD *hGenJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hGenJetPtFlavPtHatCentWeighted");

    // TCanvas *cc = new TCanvas("cc","d",800,800);
    // TH1D* hCheck = (TH1D*)hRecoJetPtFlavPtHatCentWeighted->Projection(0);
    // hCheck->Draw();

    // Centrality projections. Centrality bins correspond 0: -10-0; 1: 0-10; 2: 10-20, ...
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    // vector<Int_t> centLow{1};
    // vector<Int_t> centHi {2};

    TH1D* hRecoJetPt[ centLow.size() ];
    TH1D* hGenJetPt[ centLow.size() ];
    TH1D* hReco2GenPtRatio[ centLow.size() ];
    //TRatioPlot* hRatioPlot[ centLow.size() ];

    TCanvas *compCanv = new TCanvas("compCanv","compCanv", 800, 600);
    compCanv->Divide(centLow.size(), 2);

    std::cout << "Here" << std::endl;
    // For pp
    // vector<Int_t> centLow{0};
    // vector<Int_t> centHi {0};

    Double_t ptHatStep = 5.;
    vector<Int_t> ptHatLow{0, 1, 2, 3, 4};
    vector<Int_t> ptHatHi {119, 119, 119, 119, 119};
    Int_t ptHatBin = 0;
    Double_t ptHatVal[2] = {15. + ptHatLow.at(ptHatBin), 615.};
    updateVectorLimits( ptHatLow );
    updateVectorLimits( ptHatHi );

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    Int_t rebinFactor = 2;

    for (UInt_t iCent=0; iCent<centLow.size(); iCent++) {
        // Reconstructed jets
        std::cout << centLow.at(iCent) << " " << centHi.at(iCent) << " " << ptHatLow.at(ptHatBin) << " " << ptHatHi.at(ptHatBin) << std::endl;
        hRecoJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(iCent), centHi.at(iCent) );
        hRecoJetPtFlavPtHatCentWeighted->GetAxis(2)->SetRange( ptHatLow.at(ptHatBin), ptHatHi.at(ptHatBin) );
        hRecoJetPt[iCent] = (TH1D*)hRecoJetPtFlavPtHatCentWeighted->Projection(0);
        hRecoJetPt[iCent]->SetName(Form("hRecoJetPt_%d",iCent));
        hRecoJetPt[iCent]->Rebin(rebinFactor);
        setJetStyle(hRecoJetPt[iCent], 0);
        hRecoJetPt[iCent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

        // Generated jets
        hGenJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(iCent), centHi.at(iCent) );
        hGenJetPtFlavPtHatCentWeighted->GetAxis(2)->SetRange( ptHatLow.at(ptHatBin), ptHatHi.at(ptHatBin) );
        hGenJetPt[iCent] = (TH1D*)hGenJetPtFlavPtHatCentWeighted->Projection(0);
        hGenJetPt[iCent]->SetName(Form("hGenJetPt_%d",iCent));
        hGenJetPt[iCent]->Rebin(rebinFactor);
        setJetStyle(hGenJetPt[iCent], 1);
        hGenJetPt[iCent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

        hReco2GenPtRatio[iCent] = new TH1D(Form("hReco2GenPtRatio_%d",iCent),";p_{T} (GeV/c);Reco/Gen", 
                                            hRecoJetPt[iCent]->GetNbinsX(), 
                                            hRecoJetPt[iCent]->GetXaxis()->GetBinLowEdge(1), 
                                            hRecoJetPt[iCent]->GetXaxis()->GetBinUpEdge( hRecoJetPt[iCent]->GetNbinsX() ));
        hReco2GenPtRatio[iCent]->Divide( hRecoJetPt[iCent], hGenJetPt[iCent], 1., 1. );
        hReco2GenPtRatio[iCent]->SetName(Form("hReco2GenPtRatio_%d", iCent));
        setJetStyle(hReco2GenPtRatio[iCent], 2);

        compCanv->cd(iCent+1);
        setPadStyle();
        hRecoJetPt[iCent]->Draw();
        hGenJetPt[iCent]->Draw("same");
        gPad->SetLogy(1);
        t.DrawLatexNDC(0.25, 0.9, Form("%d-%d%% Pb+Pb #sqrt{s_{NN}}=5020 GeV",(centLow.at(iCent)-1)*10, (centHi.at(iCent))*10));

        compCanv->cd(centLow.size() + iCent+1);
        setPadStyle();
        hReco2GenPtRatio[iCent]->Draw();
        hReco2GenPtRatio[iCent]->GetXaxis()->SetRangeUser(20., 300.);
        hReco2GenPtRatio[iCent]->GetYaxis()->SetRangeUser(0.5, 3.);

    } // for (UInt_t iCent=0; iCent<centLow.size(); iCent++)

}

//________________
void drawPtHatVsLeadingJet(TFile *inFile) {
    // Retrieve histogram
    THnSparseD *hRecoLeadJetPtFlavPtHatCent = (THnSparseD*)inFile->Get("hRecoLeadJetPtFlavPtHatCent");
    THnSparseD *hRecoLeadJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoLeadJetPtFlavPtHatCentWeighted");

    // Centrality projections. Centrality bins correspond 0: -10-0; 1: 0-10; 2: 10-20, ...
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    TH2D* hPtHatVsLeadJetPt[ centLow.size() ];
    TH2D* hPtHatVsLeadJetPtWeighted[ centLow.size() ];
    TCanvas *canv = new TCanvas("canv","canv", 1600, 800);
    canv->Divide(centLow.size(), 2);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    Int_t rebinX = 1;
    Int_t rebinY = 2;

    // Loop over centralities
    for (Int_t i{0}; i<centLow.size(); i++) {

        // Without ptHatW reweighting
        hRecoLeadJetPtFlavPtHatCent->GetAxis(3)->SetRange(centLow.at(i), centHi.at(i));
        hPtHatVsLeadJetPt[i] = (TH2D*)hRecoLeadJetPtFlavPtHatCent->Projection(2,0);
        hPtHatVsLeadJetPt[i]->SetName( Form("hPtHatVsLeadJetPt_%d",i) );
        set2DStyle(hPtHatVsLeadJetPt[i]);
        hPtHatVsLeadJetPt[i]->GetXaxis()->SetTitle("p_{T}^{Lead} (GeV/c)");

        // With ptHatW reweighting
        hRecoLeadJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange(centLow.at(i), centHi.at(i));
        hPtHatVsLeadJetPtWeighted[i] = (TH2D*)hRecoLeadJetPtFlavPtHatCentWeighted->Projection(2,0);
        hPtHatVsLeadJetPtWeighted[i]->SetName( Form("hPtHatVsLeadJetPtWeighted_%d",i) );
        set2DStyle(hPtHatVsLeadJetPtWeighted[i]);
        hPtHatVsLeadJetPtWeighted[i]->GetXaxis()->SetTitle("p_{T}^{Lead} (GeV/c)");

        // Without ptHatW reweighting
        canv->cd(i+1);
        setPadStyle();
        hPtHatVsLeadJetPt[i]->Draw("colz");
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.3, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        t.DrawLatexNDC(0.2, 0.8, "Before reweighting");

        // With ptHatW reweighting
        canv->cd(i+1+centLow.size());
        setPadStyle();
        hPtHatVsLeadJetPtWeighted[i]->Draw("colz");
        hPtHatVsLeadJetPtWeighted[i]->GetXaxis()->SetRangeUser(15., 315.);
        hPtHatVsLeadJetPtWeighted[i]->GetYaxis()->SetRangeUser(15., 315.);
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.3, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        t.DrawLatexNDC(0.25, 0.8, "After reweighting");
    } // for (Int_t i{0}; i<centLow.size(); i++) 
}

//________________
void drawFindPtHatCut(TFile *inFile) {
    // Retrieve gen jet histogram
    THnSparseD *hGenJetPtFlavPtHatCent = (THnSparseD*)inFile->Get("hGenJetPtFlavPtHatCent");

    // Centrality projections. Centrality bins correspond 0: -10-0; 1: 0-10; 2: 10-20, ...
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    // PtHat projections starting 15 GeV with 5 GeV step
    Int_t ptHatMin = 15;
    Int_t pHatStep = 5;
    vector<Int_t> ptHatLow{0, 1, 2, 3, 4, 5, 6};
    vector<Int_t> ptHatHi{119, 119, 119, 119, 119, 119, 119};
    updateVectorLimits( ptHatLow );
    updateVectorLimits( ptHatHi );
    Int_t ptHatBins = ptHatLow.size();

    TH1D *hGenPt[ptHatBins];
    TH1D *hGenPtRatio[ptHatBins];

    TCanvas *canv = new TCanvas("canv","canv", 1200, 800);
    canv->Divide(2, 1);

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (Int_t i{0}; i<ptHatBins; i++) {

        // Retrieve pt spectrum of gen jet pt
        hGenJetPtFlavPtHatCent->GetAxis(2)->SetRange( ptHatLow.at(i), ptHatHi.at(i) );

        //
        // Plot gen jet pT with different ptHat on the same plot
        //
        canv->cd(1);
        setPadStyle();
        hGenPt[i] = (TH1D*)hGenJetPtFlavPtHatCent->Projection(0);
        hGenPt[i]->SetName( Form("hGenPt_%d", i) );
        hGenPt[i]->GetYaxis()->SetTitle("Entries");
        setJetStyle(hGenPt[i], i);
        hGenPt[i]->Draw("same");
        gPad->SetLogy(1);

        leg->AddEntry(hGenPt[i],Form("#hat{p_{T}}>%d GeV/c", ptHatMin + i*pHatStep), "p");
        if ( i == (ptHatBins-1) ) {
            leg->SetLineWidth(0);
            leg->Draw();
        }


        //
        // Create and plot ratios of gen jet pt spectra with ptHat>x to that with ptHat>15 GeV/c
        //
        hGenPtRatio[i] = new TH1D( Form("hGenPtRatio_%d",i),"Ratio of gen jet p_{T} distributions with different #hat{p_{T}} to #hat{p_{T}}>15 GeV/c;p_{T}^{gen} (GeV/c);p_{T} (#hat{p_{T}}>x) / p_{T} (#hat{p_{T}}>15)",
                                   hGenPt[0]->GetNbinsX(), 
                                   hGenPt[0]->GetXaxis()->GetBinLowEdge(1),
                                   hGenPt[0]->GetXaxis()->GetBinUpEdge( hGenPt[0]->GetNbinsX() ) );
        setJetStyle(hGenPtRatio[i], i);

        canv->cd(2);
        setPadStyle();
        if (i == 0) {
            hGenPtRatio[i]->Divide(hGenPt[i], hGenPt[i]);
        }
        else {
            hGenPtRatio[i]->Divide(hGenPt[i], hGenPt[0]);
        }
        hGenPtRatio[i]->Draw("same");
        hGenPtRatio[i]->GetYaxis()->SetRangeUser(0., 1.25);
        hGenPtRatio[i]->GetXaxis()->SetRangeUser(15., 215.);
    }
}

//________________
void drawPtRawCorrGen(TFile *inFile) {
    auto hRecoJetRawPtCorrPtGenPtCent = (THnSparseD*)inFile->Get("hRecoJetRawPtCorrPtGenPtCent");

    // Centrality projections. Centrality bins correspond 0: -10-0; 1: 0-10; 2: 10-20, ...
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    TH2D* hPtCorrVsPtRaw[ centLow.size() ];
    TH2D* hPtCorrVsPtGen[ centLow.size() ];
    TCanvas *canv = new TCanvas("canv","canv", 1200, 800);
    canv->Divide(centLow.size(), 2);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    Int_t rebinX = 1;
    Int_t rebinY = 1;

    for (Int_t i{0}; i<centLow.size(); i++) {
        hRecoJetRawPtCorrPtGenPtCent->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );

        // PtCorr vs. PtRaw
        hPtCorrVsPtRaw[i] = (TH2D*)hRecoJetRawPtCorrPtGenPtCent->Projection(1,0);
        hPtCorrVsPtRaw[i]->SetName( Form("hPtCorrVsPtRaw%d",i) );
        set2DStyle(hPtCorrVsPtRaw[i]);

        // PtCorr vs. PtGen
        hPtCorrVsPtGen[i] = (TH2D*)hRecoJetRawPtCorrPtGenPtCent->Projection(1,2);
        hPtCorrVsPtGen[i]->SetName( Form("hPtCorrVsPtGen%d",i) );
        set2DStyle(hPtCorrVsPtGen[i]);

        // PtCorr vs. PtRaw
        canv->cd(i+1);
        setPadStyle();
        hPtCorrVsPtRaw[i]->Draw("colz");
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.3, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));

        // PtCorr vs. PtGen
        canv->cd(i+1+centLow.size());
        setPadStyle();
        hPtCorrVsPtGen[i]->Draw("colz");
        gPad->SetLogz(1);
    } // for (Int_t i{0}; i<centLow.size(); i++)    
}

//________________
void drawCompareGenPP2PbPb(TFile *pbpbFile, TFile *ppFile) {
    // Read data for PbPb
    auto hGenJetPtEtaPhiCentPbPb = (THnSparseD*)pbpbFile->Get("hGenJetPtEtaPhiCentWeighted");
    hGenJetPtEtaPhiCentPbPb->SetName("hGenJetPtEtaPhiCentPbPb");
    auto hPtHatPbPb = (TH1D*)pbpbFile->Get("hPtHat");
    hPtHatPbPb->SetName("hPtHatPbPb");
    // Read data for pp
    auto hGenJetPtEtaPhiCentPP = (THnSparseD*)ppFile->Get("hGenJetPtEtaPhiCentWeighted");
    hGenJetPtEtaPhiCentPP->SetName("hGenJetPtEtaPhiCentPP");
    auto hPtHatPP = (TH1D*)pbpbFile->Get("hPtHat");
    hPtHatPP->SetName("hPtHatPP");

    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    vector<Int_t> ptBinLow{5};
    vector<Int_t> ptBinHi{49};
    vector<Int_t> etaBinLow{9};
    vector<Int_t> etaBinHi{39};

    updateVectorLimits( ptBinLow );
    updateVectorLimits( ptBinHi );
    updateVectorLimits( etaBinLow );
    updateVectorLimits( etaBinHi );

    // Set common pt and eta limits
    hGenJetPtEtaPhiCentPbPb->GetAxis(0)->SetRange( ptBinLow.at(0), ptBinHi.at(0) );
    hGenJetPtEtaPhiCentPbPb->GetAxis(1)->SetRange( etaBinLow.at(0), etaBinHi.at(0) );

    hGenJetPtEtaPhiCentPP->GetAxis(0)->SetRange( ptBinLow.at(0), ptBinHi.at(0) );
    hGenJetPtEtaPhiCentPP->GetAxis(1)->SetRange( etaBinLow.at(0), etaBinHi.at(0) );

    // Define histograms
    Int_t rebinPt = 1;
    Int_t rebinEta = 1;
    Int_t rebinPhi = 1;

    // Pythia + Hydjet
    TH1D* hGenJetPtPbPb[ centLow.size() ];
    TH1D* hGenJetEtaPbPb[ centLow.size() ];
    TH1D* hGenJetPhiPbPb[ centLow.size() ];
    TH2D* hGenJetPhiVsEtaPbPb[ centLow.size() ];
    
    // Pythia
    TH1D* hGenJetPtPP = (TH1D*)hGenJetPtEtaPhiCentPP->Projection(0);
    hGenJetPtPP->Scale(1./hGenJetPtPP->Integral());
    hGenJetPtPP->Rebin( rebinPt );
    TH1D* hGenJetEtaPP = (TH1D*)hGenJetPtEtaPhiCentPP->Projection(1);
    hGenJetEtaPP->Scale(1./hGenJetEtaPP->Integral());
    hGenJetEtaPP->Rebin( rebinEta );
    TH1D* hGenJetPhiPP = (TH1D*)hGenJetPtEtaPhiCentPP->Projection(2);
    hGenJetPhiPP->Scale(1./hGenJetPhiPP->Integral());
    hGenJetPhiPP->Rebin( rebinPhi );
    TH2D* hGenJetPhiVsEtaPP = (TH2D*)hGenJetPtEtaPhiCentPP->Projection(2,1);
    set2DStyle(hGenJetPhiVsEtaPP);
    TH1D* hPbPb2PPRaio[centLow.size()];
    set1DStyle(hGenJetPtPP, 1);
    set1DStyle(hGenJetEtaPP, 1);
    set1DStyle(hGenJetPhiPP, 1);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Set canvases
    TCanvas *canv = new TCanvas("canv","canv", 1600, 800);
    canv->Divide( centLow.size(), 2);

    TCanvas *canv2 = new TCanvas("canv2","canv2", 1200, 800);
    canv2->Divide( centLow.size()+1, 3);
    
    canv2->cd(1);
    setPadStyle();
    hGenJetPhiVsEtaPP->Draw("colz");
    gPad->SetLogz(1);
    t.DrawLatexNDC(0.3, 0.85, "PYTHIA");

    TLegend *leg;

    for (Int_t i{0}; i<centLow.size(); i++) {
        hGenJetPtEtaPhiCentPbPb->GetAxis(3)->SetRange(centLow.at(i), centHi.at(i));
        hGenJetPtPbPb[i] = (TH1D*)hGenJetPtEtaPhiCentPbPb->Projection(0);
        hGenJetPtPbPb[i]->SetName(Form("%s_%d",hGenJetPtEtaPhiCentPbPb->GetName(),i));
        hGenJetPtPbPb[i]->Rebin( rebinPt );
        set1DStyle(hGenJetPtPbPb[i], 0);
        hGenJetEtaPbPb[i] = (TH1D*)hGenJetPtEtaPhiCentPbPb->Projection(1);
        hGenJetEtaPbPb[i]->SetName(Form("%s_%d",hGenJetPtEtaPhiCentPbPb->GetName(),i));
        hGenJetEtaPbPb[i]->Rebin( rebinEta );
        set1DStyle(hGenJetEtaPbPb[i], 0);
        hGenJetPhiPbPb[i] = (TH1D*)hGenJetPtEtaPhiCentPbPb->Projection(2);
        hGenJetPhiPbPb[i]->SetName(Form("%s_%d",hGenJetPtEtaPhiCentPbPb->GetName(),i));
        hGenJetPhiPbPb[i]->Rebin( rebinPhi );
        set1DStyle(hGenJetPhiPbPb[i], 0);


        hGenJetPhiVsEtaPbPb[i] = (TH2D*)hGenJetPtEtaPhiCentPbPb->Projection(2,1);
        hGenJetPhiVsEtaPbPb[i]->SetName(Form("hGenJetPhiVsEtaPbPb_%d",i));
        set2DStyle(hGenJetPhiVsEtaPbPb[i]);

        canv->cd(i+1);
        setPadStyle();
        hGenJetPtPbPb[i]->Scale(1./hGenJetPtPbPb[i]->Integral());
        hGenJetPtPbPb[i]->Draw();
        hGenJetPtPP->Draw("same");
        gPad->SetLogy(1);

        if ( i == 0) {
            leg = new TLegend(0.5, 0.65, 0.9, 0.75);
            leg->AddEntry(hGenJetPtPbPb[i],"PYTHIA+HYDJET","p");
            leg->AddEntry(hGenJetPtPP,"PYTHIA","p");
            leg->SetBorderSize(0);
            leg->SetTextSize(0.05);
            leg->Draw();
            t.DrawLatexNDC(0.55, 0.78, "#hat{p_{T}} > 50 GeV/c");
        }
        t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        

        canv->cd(i+1+centLow.size());
        setPadStyle();
        hPbPb2PPRaio[i] = new TH1D( Form("hPbPb2PPRaio_%d",i), "Ratio of gen jet p_{T} spectra in PbPb to pp;p_{T}^{gen} (GeV/c);dN^{PbPb}/dp_{T} / dN^{pp}/dp_{T}",
                                    hGenJetPtPP->GetNbinsX(),
                                    hGenJetPtPP->GetXaxis()->GetBinLowEdge(1),
                                    hGenJetPtPP->GetXaxis()->GetBinUpEdge( hGenJetPtPP->GetNbinsX() ) );
        hPbPb2PPRaio[i]->Sumw2();
        hPbPb2PPRaio[i]->Divide(hGenJetPtPbPb[i], hGenJetPtPP);
        setJetStyle(hPbPb2PPRaio[i], 2);
        hPbPb2PPRaio[i]->Draw();

        canv2->cd(i+2);
        setPadStyle();
        hGenJetPhiVsEtaPbPb[i]->Draw("colz");
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.3, 0.85, Form("PYTHIA+HYDJET %d-%d%%", (centLow.at(i)-2)*10,  (centHi.at(i)-1)*10) );

        canv2->cd(i+2+centLow.size()+1);
        setPadStyle();
        hGenJetEtaPbPb[i]->Scale(1./hGenJetEtaPbPb[i]->Integral());
        hGenJetEtaPbPb[i]->Draw();
        hGenJetEtaPP->Draw("same");

        canv2->cd(i + 2 + 2*centLow.size() + 2);
        setPadStyle();
        hGenJetPhiPbPb[i]->Scale(1./hGenJetPhiPbPb[i]->Integral());
        hGenJetPhiPbPb[i]->Draw();
        hGenJetPhiPP->Draw("same");
    } // for (Int_t i{0}; i<centLow.size(); i++)
}

//________________
void drawJESvsCentrality(TFile *inFile) {

    // Jet histograms
    THnSparseD *hJESPtEtaPhiCent = (THnSparseD*)inFile->Get("hJESPtEtaPhiCent");
    THnSparseD *hJESPtEtaPhiCentWeighted = (THnSparseD*)inFile->Get("hJESPtEtaPhiCentWeighted");

    // Centrality projections. Centrality bins correspond 0: -10-0; 1: 0-10; 2: 10-20, ...
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {8, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    vector<Int_t> ptLowBin{9};
    vector<Int_t> ptHiBin{49};
    updateVectorLimits( ptLowBin );
    updateVectorLimits( ptHiBin );

    TH2D *hJESvsGenPt[ centLow.size() ];
    TH1D *hJESMean[ centLow.size() ];
    TH1D *hJESSigma[ centLow.size() ];
    TH2D *hJESvsGenPtW[ centLow.size() ];
    TH1D *hJESMeanW[ centLow.size() ];
    TH1D *hJESSigmaW[ centLow.size() ];

    TCanvas *canvJES = new TCanvas("canvJES","canvJES", 1400, 800);
    canvJES->Divide( centLow.size(), 3 );

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    Int_t rebinX = 2;
    Int_t rebinY = 4;

    TLegend *leg;
    // Centrality dependence
    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {

        //std::cout << "iCent: " << iCent << " padIndex: " << padIndex << " centHistoBin: " << centHistoBin << std::endl;
        
        // JES
        hJESPtEtaPhiCent->GetAxis(4)->SetRange( centLow.at(iCent), centHi.at(iCent) );
        hJESvsGenPt[iCent] = (TH2D*)hJESPtEtaPhiCent->Projection(0,1);
        hJESvsGenPt[iCent]->SetName(Form("hJESvsGenPt_%d",iCent));
        hJESvsGenPt[iCent]->RebinX(rebinX);
        hJESvsGenPt[iCent]->RebinY(rebinY);

        hJESvsGenPt[iCent]->FitSlicesY();
        hJESMean[iCent] = (TH1D*)gDirectory->Get(Form("hJESvsGenPt_%d_%d", iCent, 1));
        hJESMean[iCent]->SetName(Form("hJESMean_%d", iCent));
        hJESSigma[iCent] = (TH1D*)gDirectory->Get(Form("hJESvsGenPt_%d_%d", iCent, 2));
        hJESSigma[iCent]->SetName(Form("hJESSigma_%d",iCent));

        // JES ptHat weighted
        hJESPtEtaPhiCentWeighted->GetAxis(4)->SetRange( centLow.at(iCent), centHi.at(iCent) );
        hJESvsGenPtW[iCent] = (TH2D*)hJESPtEtaPhiCentWeighted->Projection(0,1);
        hJESvsGenPtW[iCent]->SetName(Form("hJESvsGenPtW_%d",iCent));
        hJESvsGenPtW[iCent]->RebinX(rebinX);
        hJESvsGenPtW[iCent]->RebinY(rebinY);

        hJESvsGenPtW[iCent]->FitSlicesY();
        hJESMeanW[iCent] = (TH1D*)gDirectory->Get(Form("hJESvsGenPtW_%d_%d", iCent, 1));
        hJESMeanW[iCent]->SetName(Form("hJESMeanW_%d", iCent));
        hJESSigmaW[iCent] = (TH1D*)gDirectory->Get(Form("hJESvsGenPtW_%d_%d", iCent, 2));
        hJESSigmaW[iCent]->SetName(Form("hJESSigmaW_%d",iCent));

        // JES (2D)
        canvJES->cd(iCent+1);
        setPadStyle();
        gPad->SetLeftMargin(0.17);
        set2DStyle(hJESvsGenPt[iCent]);
        hJESvsGenPt[iCent]->Draw("colz");
        hJESvsGenPt[iCent]->GetYaxis()->SetRangeUser(-1., 2.);
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.35, 0.9, Form("JES: %d-%d%% PYTHIA+HYDJET",(centLow.at(iCent)-2)*10, (centHi.at(iCent)-1)*10));

        // Jet energy scale
        canvJES->cd(iCent + 1 + centLow.size());
        setPadStyle();
        set1DStyle(hJESMean[iCent],0);
        hJESMean[iCent]->Draw();
        hJESMean[iCent]->GetYaxis()->SetRangeUser(0.9, 1.3);
        hJESMean[iCent]->GetYaxis()->SetTitle("#mu(JES)");
        set1DStyle(hJESMeanW[iCent],1);
        hJESMeanW[iCent]->Draw("same");
        if (iCent==0) {
            leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg->AddEntry(hJESMean[iCent],"w/o weight","p");
            leg->AddEntry(hJESMeanW[iCent]," w  weight","p");
            //leg->SetBorderColor(0);
            leg->SetBorderSize(0);
            leg->Draw();
        }

        // Jet energy resolution
        canvJES->cd(iCent + 1 + 2 * centLow.size());
        setPadStyle();
        set1DStyle(hJESSigma[iCent],0);
        hJESSigma[iCent]->Draw();
        hJESSigma[iCent]->GetYaxis()->SetRangeUser(0., 0.35);
        hJESSigma[iCent]->GetYaxis()->SetTitle("#sigma(JES)");
        set1DStyle(hJESSigmaW[iCent],1);
        hJESSigmaW[iCent]->Draw("same");
    } // for (Int_t iCent=1; iCent<=centBins; iCent++)

}

//________________
void drawJESvsEta(TFile* inFile) {

    // Histograms to read
    THnSparseD *hJESPtEtaPhiCent = (THnSparseD*)inFile->Get("hJESPtEtaPhiCent");
    THnSparseD *hJESPtEtaPhiCentWeighted = (THnSparseD*)inFile->Get("hJESPtEtaPhiCentWeighted");

    // Centrality projections. Centrality bins correspond 0: -10-0; 1: 0-10; 2: 10-20, ...
    vector<Int_t> centLow{ 6, 4, 1 };
    vector<Int_t> centHi { 8, 5, 3 };
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    Double_t etaLowVal = -2.5;
    Double_t etaStep = 0.1;
    
    vector<Int_t> etaLow{ 10, 15, 20, 25, 30, 35 };
    vector<Int_t> etaHi { 14, 19, 24, 29, 34, 39 };
    updateVectorLimits( etaLow );
    updateVectorLimits( etaHi );

    // Histograms to fill
    TH2D *hJESVsGenPtEtaDep[centLow.size()][ etaLow.size() ];
    TH1D *hJESMeanEtaDep[centLow.size()][etaLow.size()];
    TH1D *hJESSigmaEtaDep[centLow.size()][ etaLow.size() ];
    TH2D *hJESVsGenPtEtaDepW[centLow.size()][ etaLow.size() ];
    TH1D *hJESMeanEtaDepW[centLow.size()][ etaLow.size() ];
    TH1D *hJESSigmaEtaDepW[centLow.size()][ etaLow.size() ];

    // Canvases
    TCanvas *canvJESeta[centLow.size()];
    for(Int_t i{0}; i<centLow.size(); i++) {
        canvJESeta[i] = new TCanvas( Form("canvJESeta_%d",i),
                                     Form("canvJESeta_%d",i), 
                                     1400, 800);
        canvJESeta[i]->Divide( etaLow.size(), 3 );
    }

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TLegend *leg[centLow.size()];

    // Loop over centrality bins
    for( Int_t i{0}; i<centLow.size(); i++) {
        hJESPtEtaPhiCent->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hJESPtEtaPhiCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );

        Int_t rebinX{4};
        Int_t rebinY{4};

        // Loop over pseudorapidity bins
        for (Int_t iEta{0}; iEta<etaLow.size(); iEta++) {

            // JES

            // Unweighted
            hJESPtEtaPhiCent->GetAxis(2)->SetRange( etaLow.at(iEta), etaHi.at(iEta) );
            hJESVsGenPtEtaDep[i][iEta] = (TH2D*)hJESPtEtaPhiCent->Projection(0,1);
            hJESVsGenPtEtaDep[i][iEta]->SetName(Form("hJESVsGenPtEtaDep_%d_%d",i , iEta));
            hJESVsGenPtEtaDep[i][iEta]->RebinX(rebinX);
            hJESVsGenPtEtaDep[i][iEta]->RebinY(rebinY);

            // Weighted
            hJESPtEtaPhiCentWeighted->GetAxis(2)->SetRange( etaLow.at(iEta), etaHi.at(iEta) );
            hJESVsGenPtEtaDepW[i][iEta] = (TH2D*)hJESPtEtaPhiCentWeighted->Projection(0,1);
            hJESVsGenPtEtaDepW[i][iEta]->SetName(Form("hJESVsGenPtEtaDepW_%d_%d", i, iEta));
            hJESVsGenPtEtaDepW[i][iEta]->RebinX(rebinX);
            hJESVsGenPtEtaDepW[i][iEta]->RebinY(rebinY);

            // Unweighted fit slices
            hJESVsGenPtEtaDep[i][iEta]->FitSlicesY();
            hJESMeanEtaDep[i][iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDep_%d_%d_%d", i, iEta, 1));
            hJESMeanEtaDep[i][iEta]->SetName(Form("hJESMeanEtaDep_%d_%d", i, iEta));
            hJESSigmaEtaDep[i][iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDep_%d_%d_%d", i, iEta, 2));
            hJESSigmaEtaDep[i][iEta]->SetName(Form("hJESSigmaEtaDep_%d_%d", i, iEta));

            // Weighted fit slices
            hJESVsGenPtEtaDepW[i][iEta]->FitSlicesY();
            hJESMeanEtaDepW[i][iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDepW_%d_%d_%d", i, iEta, 1));
            hJESMeanEtaDepW[i][iEta]->SetName(Form("hJESMeanEtaDepW_%d_%d", i, iEta));
            hJESSigmaEtaDepW[i][iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDepW_%d_%d_%d", i, iEta, 2));
            hJESSigmaEtaDepW[i][iEta]->SetName(Form("hJESSigmaEtaDepW_%d_%d", i, iEta));


            //
            // For each centrality class draw JES as a funciton of eta
            //

            // JES (2D)
            canvJESeta[i]->cd(iEta + 1);
            setPadStyle();
            gPad->SetLeftMargin(0.17);
            set2DStyle(hJESVsGenPtEtaDep[i][iEta]);
            hJESVsGenPtEtaDep[i][iEta]->Draw("colz");
            hJESVsGenPtEtaDep[i][iEta]->GetYaxis()->SetRangeUser(-1., 2.);
            gPad->SetLogz(1);
            t.DrawLatexNDC(0.25, 0.9, Form("JES: %2.1f<#eta<%2.1f Pb+Pb #sqrt{s_{NN}}=5020 GeV", etaLowVal + ( etaLow.at(iEta) - 1 ) * etaStep, etaLowVal + ( etaHi.at(iEta) ) * etaStep));
            t.DrawLatexNDC(0.5, 0.8, Form("%d-%d%% centrality", ( centLow.at(i) - 2 ) * 10, ( centHi.at(i) - 1 ) * 10));

            // Jet energy scale
            canvJESeta[i]->cd(etaLow.size() + iEta + 1);
            setPadStyle();
            set1DStyle(hJESMeanEtaDep[i][iEta],0);
            hJESMeanEtaDep[i][iEta]->Draw();
            hJESMeanEtaDep[i][iEta]->GetYaxis()->SetRangeUser(0.9, 1.3);
            hJESMeanEtaDep[i][iEta]->GetYaxis()->SetTitle("#mu(JES)");
            set1DStyle(hJESMeanEtaDepW[i][iEta],1);
            hJESMeanEtaDepW[i][iEta]->Draw("same");
            if (iEta == 0) {
            leg[i] = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg[i]->AddEntry(hJESMeanEtaDep[i][iEta],"w/o weight","p");
            leg[i]->AddEntry(hJESMeanEtaDepW[i][iEta]," w  weight","p");
            leg[i]->SetBorderSize(0);
            leg[i]->Draw();
        }

            // Jet energy resolution
            canvJESeta[i]->cd( 2 * etaLow.size() + iEta + 1 );
            setPadStyle();
            set1DStyle(hJESSigmaEtaDep[i][iEta],0);
            hJESSigmaEtaDep[i][iEta]->Draw();
            hJESSigmaEtaDep[i][iEta]->GetYaxis()->SetRangeUser(0., 0.35);
            hJESSigmaEtaDep[i][iEta]->GetYaxis()->SetTitle("#sigma(JES)");
            set1DStyle(hJESSigmaEtaDepW[i][iEta],1);
            hJESSigmaEtaDepW[i][iEta]->Draw("same");
        } // for (Int_t iEta=0; iEta<etaBins; iEta++)

    } // for( Int_t i{0}; i<centLow.size(); i++)
}

//________________
void drawRecoUnmatched2Inclusive(TFile *inFile) {

    // Need to plot ptHat vs jet pt

    // Read N-dimensional histograms
    auto hRecoJetPtFlavPtHatCentInclusiveWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentInclusive");
    auto hRecoUnmatchedJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCent");
    // Matched jets
    auto hRecoJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCent");

    //auto hRecoJetPtFlavPtHatCentInclusiveWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentInclusiveWeighted");
    //auto hRecoUnmatchedJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCentWeighted");
    // Matched jets
    //auto hRecoJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentWeighted");

    // Setup centralities to read
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    // Histograms to fill
    TH1D *hPtInclusiveW[ centLow.size() ];
    TH1D *hPtUnmatchedW[ centLow.size() ];
    TH1D* hPtMatchedW[ centLow.size() ];

    TH2D *hPtHatVsInclusivePt[ centLow.size() ];
    TH2D *hPtHatVsUnmatchedPt[ centLow.size() ];

    // Create canvas
    TCanvas *canv = new TCanvas("canv","canv", 1200, 400);
    canv->Divide( centLow.size(), 1 );

    TCanvas *canv2 = new TCanvas("canv2", "canv2", 1200, 800);
    canv2->Divide( centLow.size(), 2 );

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TLegend *leg;

    // Loop over centralities
    for (Int_t i{0}; i<centLow.size(); i++) {
        // Fill inclusive reco jets
        hRecoJetPtFlavPtHatCentInclusiveWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hPtInclusiveW[i] = (TH1D*)hRecoJetPtFlavPtHatCentInclusiveWeighted->Projection(0);
        hPtInclusiveW[i]->SetName(Form("hPtInclusiveW_%d",i));
        set1DStyle(hPtInclusiveW[i], 1);
        hPtInclusiveW[i]->Scale(1./hPtInclusiveW[i]->Integral());
        hPtInclusiveW[i]->GetYaxis()->SetTitle("1/N_{jet} dN_{jet}/dp_{T}^{reco}");

        hPtHatVsInclusivePt[i] = (TH2D*)hRecoJetPtFlavPtHatCentInclusiveWeighted->Projection(2,0);
        hPtHatVsInclusivePt[i]->SetName(Form("hPtHatVsInclusivePt_%d",i));
        set2DStyle( hPtHatVsInclusivePt[i] );

        // Fill unmatched reco jets
        hRecoUnmatchedJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hPtUnmatchedW[i] = (TH1D*)hRecoUnmatchedJetPtFlavPtHatCentWeighted->Projection(0);
        hPtUnmatchedW[i]->SetName(Form("hPtUnmatchedW_%d",i));
        set1DStyle(hPtUnmatchedW[i], 0);
        hPtUnmatchedW[i]->Scale(1./hPtUnmatchedW[i]->Integral());

        hPtHatVsUnmatchedPt[i] = (TH2D*)hRecoUnmatchedJetPtFlavPtHatCentWeighted->Projection(2,0);
        hPtHatVsUnmatchedPt[i]->SetName( Form("hPtHatVsUnmatchedPt[i]_%d",i) );
        set2DStyle( hPtHatVsUnmatchedPt[i] );

        // Fill matched reco jets
        hRecoJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hPtMatchedW[i] = (TH1D*)hRecoJetPtFlavPtHatCentWeighted->Projection(0);
        hPtMatchedW[i]->SetName( Form("hPtMatchedW_%d",i) );
        setJetStyle(hPtMatchedW[i], 2);
        hPtMatchedW[i]->Scale(1./hPtMatchedW[i]->Integral());

        // Compare unmatched to inclusive
        canv->cd(i+1);
        setPadStyle();
        hPtInclusiveW[i]->Draw();
        hPtUnmatchedW[i]->Draw("same");
        hPtMatchedW[i]->Draw("same");
        gPad->SetLogy(1);
        t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        if ( i == 0) {
            leg = new TLegend(0.5, 0.75, 0.9, 0.85);
            leg->AddEntry(hPtInclusiveW[i],"Inclusive","p");
            leg->AddEntry(hPtUnmatchedW[i],"Unmatched","p");
            leg->AddEntry(hPtMatchedW[i],"Matched","p");
            leg->SetBorderSize(0);
            leg->SetTextSize(0.05);
            leg->Draw();
            //t.DrawLatexNDC(0.55, 0.78, "#hat{p_{T}} > 50 GeV/c");
        }

        // Plot ptHat vs pT of inclusive jet pT and unmatched jet pT
        canv2->cd(i+1);
        setPadStyle();
        hPtHatVsInclusivePt[i]->Draw("colz");
        //gPad->SetLogz(1);
        t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));

        canv2->cd( i+1+centLow.size() );
        setPadStyle();
        hPtHatVsUnmatchedPt[i]->Draw("colz");
        //gPad->SetLogz(1);
        t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
    }

    //
    // Compare inclusive jet pT spectra for different pThat
    //

    Int_t ptHatMin = 15;
    Int_t ptHatStep = 5;
    vector<Int_t> ptHatLow{0, 1, 2, 3, 4, 5, 6};
    vector<Int_t> ptHatHi{119, 119, 119, 119, 119, 119, 119};
    updateVectorLimits( ptHatLow );
    updateVectorLimits( ptHatHi );

    TH1D* hInclusivePt[ centLow.size() ][ ptHatLow.size() ];
    TH1D* hInclusivePtratio[ centLow.size() ][ ptHatLow.size() ];

    TCanvas *canv3 = new TCanvas("canv3", "canv3", 1200, 400);
    canv3->Divide( centLow.size(), 1 );

    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {
        hRecoJetPtFlavPtHatCentInclusiveWeighted->GetAxis(3)->SetRange( centLow.at(iCent), centHi.at(iCent) );
        leg = new TLegend(0.5, 0.75, 0.9, 0.85);
        for (Int_t iPtHat{0}; iPtHat<ptHatLow.size(); iPtHat++) {            
            hRecoJetPtFlavPtHatCentInclusiveWeighted->GetAxis(2)->SetRange( ptHatLow.at(iPtHat), ptHatHi.at(iPtHat) );
            hInclusivePt[iCent][iPtHat] = (TH1D*)hRecoJetPtFlavPtHatCentInclusiveWeighted->Projection(0);
            hInclusivePt[iCent][iPtHat]->SetName( Form("hInclusivePt_%d_%d",iCent,iPtHat) );
            setJetStyle( hInclusivePt[iCent][iPtHat], iPtHat );
            Int_t nBinsPt = hInclusivePt[iCent][iPtHat]->GetNbinsX();
            hInclusivePt[iCent][iPtHat]->Scale(1./hInclusivePt[iCent][iPtHat]->Integral( nBinsPt / 2, nBinsPt ) );
            hInclusivePt[iCent][iPtHat]->GetYaxis()->SetTitle("1/N_{jets} dN/dp_{T}");

            canv3->cd(iCent+1);
            setPadStyle();
            hInclusivePt[iCent][iPtHat]->Draw("same");
            gPad->SetLogy(1);

            
            if ( iCent == 0 ) {
                leg->AddEntry(hInclusivePt[0][iPtHat],Form("#hat{p_{T}}>%d GeV/c",ptHatMin + iPtHat*ptHatStep), "p");
                leg->SetBorderSize(0);
                //t.DrawLatexNDC(0.55, 0.78, "#hat{p_{T}} > 50 GeV/c");
            }
        } // for (Int_t iPtHat{0}; iPtHat<ptHatLow.size(); iPtHat++)
        t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(iCent)-2)*10, (centHi.at(iCent)-1)*10));
        if ( iCent == 0 ) {
            leg->SetTextSize(0.05);
            leg->Draw();
        }
    } // for (Int_t iCent{0}; iCent<centLow.size(); iCent++)
}

//________________
void drawRecoUnmatched(TFile *inFile) {

    auto hVzPtHatCent = (THnSparseD*)inFile->Get("hVzPtHatCent");
    auto hVzPtHatCentWeighted = (THnSparseD*)inFile->Get("hVzPtHatCentWeighted");
    auto hRecoUnmatchedJetPtFlavPtHatCent = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCent");
    auto hRecoUnmatchedJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCentWeighted");
    // Matched
    auto hRecoJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentWeighted");

    // Setup centralities to read
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    Int_t ptHatLowVal = 15;
    Int_t ptHatStep = 5;
    vector<Int_t> ptHatLow{10};
    vector<Int_t> ptHatHi{119};
    updateVectorLimits( ptHatLow );
    updateVectorLimits( ptHatHi );

    // Create histograms
    TH1D *hVz[ centLow.size() ];
    TH1D *hVzW[ centLow.size() ];
    TH1D *hPtUnmatched[ centLow.size() ];
    TH1D *hPtUnmatchedW[ centLow.size() ];
    TH1D *hPtMatchedW[ centLow.size() ];
    TH1D *hPtWeighted2Unweighted[ centLow.size() ];
    TH1D *hPtUnmatched2MatchedW[ centLow.size() ];

    auto canv = new TCanvas("canv", "canv", 1200, 800);
    canv->Divide( centLow.size(), 3 );

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TLegend *leg;

    hVzPtHatCent->GetAxis(1)->SetRange( ptHatLow.at(0), ptHatHi.at(0) );
    hVzPtHatCentWeighted->GetAxis(1)->SetRange( ptHatLow.at(0), ptHatHi.at(0) );

    hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(2)->SetRange( ptHatLow.at(0), ptHatHi.at(0) );
    hRecoUnmatchedJetPtFlavPtHatCentWeighted->GetAxis(2)->SetRange( ptHatLow.at(0), ptHatHi.at(0) );
    hRecoJetPtFlavPtHatCentWeighted->GetAxis(2)->SetRange( ptHatLow.at(0), ptHatHi.at(0) );

    // Loop over centralities
    for (Int_t i{0}; i<centLow.size(); i++) {
        
        // Select centrality range
        //hVzPtHatCent->GetAxis(1)->SetRange( ptHatLow.at(i), ptHatHi.at(i) );
        //hVzPtHatCentWeighted->GetAxis(1)->SetRange( ptHatLow.at(i), ptHatHi.at(i) );

        std::cout << "Here" << std::endl;
        hVzPtHatCent->GetAxis(2)->SetRange( centLow.at(i), centHi.at(i) );
        hVzPtHatCentWeighted->GetAxis(2)->SetRange( centLow.at(i), centHi.at(i) );

        std::cout << "Here2" << std::endl;
        //hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(2)->SetRange( ptHatLow.at(i), ptHatHi.at(i) );
        //hRecoUnmatchedJetPtFlavPtHatCentWeighted->GetAxis(2)->SetRange( ptHatLow.at(i), ptHatHi.at(i) );
        //hRecoJetPtFlavPtHatCentWeighted->GetAxis(2)->SetRange( ptHatLow.at(i), ptHatHi.at(i) );

        std::cout << "Here3" << std::endl;
        hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoUnmatchedJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        
        // Event-wise histograms and integrals
        hVz[i] = (TH1D*)hVzPtHatCent->Projection(0);
        hVz[i]->SetName(Form("hVz_%d",i));
        hVzW[i] = (TH1D*)hVzPtHatCentWeighted->Projection(0);
        hVzW[i]->SetName(Form("hVzW_%d",i));
        Double_t scaleUnweighted = hVz[i]->Integral();
        Double_t scaleWeighted = hVzW[i]->Integral();

        // Project spectra
        hPtUnmatched[i] = (TH1D*)hRecoUnmatchedJetPtFlavPtHatCent->Projection(0);
        hPtUnmatched[i]->SetName(Form("hPtUnmatched_%d",i));
        hPtUnmatchedW[i] = (TH1D*)hRecoUnmatchedJetPtFlavPtHatCentWeighted->Projection(0);
        hPtUnmatchedW[i]->SetName(Form("hPtUnmatchedW_%d",i));
        hPtMatchedW[i] = (TH1D*)hRecoJetPtFlavPtHatCentWeighted->Projection(0);
        hPtMatchedW[i]->SetName(Form("hPtMatchedW_%d",i));
        setJetStyle(hPtUnmatched[i], 0);
        setJetStyle(hPtUnmatchedW[i], 1);
        setJetStyle(hPtMatchedW[i], 3);

        hPtUnmatched[i]->Scale( 1./scaleUnweighted );
        hPtUnmatchedW[i]->Scale( 1./scaleWeighted );
        hPtMatchedW[i]->Scale( 1./scaleWeighted );

        hPtWeighted2Unweighted[i] = new TH1D( Form("hPtWeighted2Unweighted_%d",i),"Ratio of unweighted to weighted;p_{T}^{corr} (GeV/c); weighted / unweighted",
                                              hPtUnmatched[i]->GetNbinsX(),
                                              hPtUnmatched[i]->GetXaxis()->GetBinLowEdge(1),
                                              hPtUnmatched[i]->GetXaxis()->GetBinUpEdge( hPtUnmatched[i]->GetNbinsX() ) );
        hPtWeighted2Unweighted[i]->Sumw2();
        setJetStyle(hPtWeighted2Unweighted[i], 2);
        hPtWeighted2Unweighted[i]->Divide( hPtUnmatchedW[i], hPtUnmatched[i] );

        hPtUnmatched2MatchedW[i] = new TH1D( Form("hPtUnmatched2MatchedW_%d",i),"Unmatched over matched (both weighted);p_{T}^{corr} (GeV/c);unmatched / matched",
                                             hPtUnmatched[i]->GetNbinsX(),
                                             hPtUnmatched[i]->GetXaxis()->GetBinLowEdge(1),
                                             hPtUnmatched[i]->GetXaxis()->GetBinUpEdge( hPtUnmatched[i]->GetNbinsX() ) );
        setJetStyle(hPtUnmatched2MatchedW[i], 5);
        hPtUnmatched2MatchedW[i]->Divide( hPtUnmatchedW[i], hPtMatchedW[i] );


        // Plot spectra on the same plot
        canv->cd(i + 1);
        setPadStyle();
        hPtUnmatched[i]->Draw();
        hPtUnmatchedW[i]->Draw("same");
        hPtMatchedW[i]->Draw("same");
        gPad->SetLogy(1);
        t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        if ( i == 0) {
            leg = new TLegend(0.5, 0.75, 0.9, 0.85);
            leg->AddEntry(hPtUnmatched[i],"Unweighted","p");
            leg->AddEntry(hPtUnmatchedW[i],"Weighted","p");
            leg->AddEntry(hPtMatchedW[i],"Matched (w)", "p");
            leg->SetBorderSize(0);
            leg->SetTextSize(0.05);
            leg->Draw();
            //t.DrawLatexNDC(0.55, 0.78, "#hat{p_{T}} > 50 GeV/c");
        }

        // Plot ratio of unmatched weighted to unweighted
        canv->cd(i + 1 + centLow.size());
        setPadStyle();
        hPtWeighted2Unweighted[i]->Draw();
        hPtWeighted2Unweighted[i]->GetYaxis()->SetRangeUser(-0.5, 2.);

        // Plot ratio of matched to unmatched (both weighted)
        canv->cd(i + 1 + 2 * centLow.size());
        setPadStyle();
        hPtUnmatched2MatchedW[i]->Draw();
        hPtUnmatched2MatchedW[i]->GetYaxis()->SetRangeUser(0.01, 100.);
        gPad->SetLogy(1);
    } // for (Int_t i{0}; i<centLow.size(); i++)
}

//________________
void drawPtHatVsRecoJets(TFile* inFile, Bool_t isPbPb) {

    // Matched (centrality * ptHat weighted)
    auto hRecoJetPtFlavPtHatCent = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentWeighted");
    // To protect from possible crash due to the name inconsistency
    hRecoJetPtFlavPtHatCent->SetName("hRecoJetPtFlavPtHatCent");
    // Inclusive (centrality weighted)
    auto hRecoJetPtFlavPtHatCentInclusive = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentInclusiveWeighted");
    // Unmatched (centrality weighted)
    auto hRecoUnmatchedJetPtFlavPtHatCent = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCent");
    // Leading jet (centrality * ptHat weighted)
    auto hRecoLeadJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoLeadJetPtFlavPtHatCentWeighted");
    // Gen jet (centrality * ptHat weighted)
    auto hGenJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hGenJetPtFlavPtHatCentWeighted");

    // Event information
    auto hVzPtHatCent = (THnSparseD*)inFile->Get("hVzPtHatCent");

    vector<Int_t> centLow;
    vector<Int_t> centHi;
    // Setup centralities to read (PbPb)
    if (isPbPb) {
        centLow.push_back(6);
        centLow.push_back(4);
        centLow.push_back(2);
        centLow.push_back(1);

        centHi.push_back(8);
        centHi.push_back(5);
        centHi.push_back(3);
        centHi.push_back(1);
    }
    else {
        centLow.push_back(0);

        centHi.push_back(9);
    }
    updateVectorLimits( centLow );
    updateVectorLimits( centHi );

    Int_t ptLowVal = 20;
    Int_t ptStep = 10;
    // vector<Int_t> jetPtLow{0, 1, 2, 3, 5};
    // vector<Int_t> jetPtHi{1, 2, 3, 5, 49};

    // vector<Int_t> jetPtLow{4, 5, 6, 7, 8 };
    // vector<Int_t> jetPtHi{4, 5, 6, 7, 8  };

    vector<Int_t> jetPtLow{4, 6, 10 };
    vector<Int_t> jetPtHi {7, 9, 12  };
    updateVectorLimits( jetPtLow );
    updateVectorLimits( jetPtHi );

    Int_t ptHatLowVal = 15;
    Int_t ptHatStep = 5;
    vector<Int_t> ptHatLow{0, 5, 10, 20, 40};
    vector<Int_t> ptHatHi {4, 9, 19, 39, 119};
    updateVectorLimits( ptHatLow );
    updateVectorLimits( ptHatHi );

    TCanvas *canvVz;
    TCanvas *canvVzCentForPtHat[ centLow.size() ];
    TCanvas *canvUnmatchedJetPtForPtHat;
    TCanvas *canvMatchedJetPtHat;
    TCanvas *canvUnmatchedJetPtHat;
    TCanvas *canvInclusiveJetPtHat;
    TCanvas *canvLeadJetPtHat;
    TCanvas *canvGenJetPtHat;

    if ( isPbPb ) {
        canvVz = new TCanvas("canvVz","canvVz", 1200, 300);
        for (Int_t i{0}; i<centLow.size(); i++) {
            canvVzCentForPtHat[i] = new TCanvas(Form("canvVzCentForPtHat_%d",i), Form("canvVzCentForPtHat_%d",i), 1200, 300);
            canvVzCentForPtHat[i]->Divide( ptHatLow.size(), 1);
        }
        canvUnmatchedJetPtForPtHat = new TCanvas("canvUnmatchedJetPtForPtHat","canvUnmatchedJetPtForPtHat", 1200, 800);
        canvMatchedJetPtHat = new TCanvas("canvMatchedJetPtHat","canvMatchedJetPtHat", 1200, 800);
        canvUnmatchedJetPtHat = new TCanvas("canvUnmatchedJetPtHat","canvUnmatchedJetPtHat", 1200, 800);
        canvInclusiveJetPtHat = new TCanvas("canvInclusiveJetPtHat","canvInclusiveJetPtHat", 1200, 800);
        canvLeadJetPtHat = new TCanvas("canvLeadJetPtHat","canvLeadJetPtHat", 1200, 800);
        canvGenJetPtHat = new TCanvas("canvGenJetPtPtHat","canvGenJetPtPtHat", 1200, 800);
    }
    else {
        canvVz = new TCanvas("canvVz","canvVz", 300, 300);
        canvUnmatchedJetPtForPtHat = new TCanvas("canvUnmatchedJetPtForPtHat","canvUnmatchedJetPtForPtHat", 1200, 800);
        for (Int_t i{0}; i<centLow.size(); i++) {
            canvVzCentForPtHat[i] = new TCanvas(Form("canvVzCentForPtHat_%d",i), Form("canvVzCentForPtHat_%d",i), 1200, 300);
            canvVzCentForPtHat[i]->Divide( ptHatLow.size(), 1);
        }
        canvMatchedJetPtHat = new TCanvas("canvMatchedJetPtHat","canvMatchedJetPtHat", 300, 800);
        canvUnmatchedJetPtHat = new TCanvas("canvUnmatchedJetPtHat","canvUnmatchedJetPtHat", 300, 800);
        canvInclusiveJetPtHat = new TCanvas("canvInclusiveJetPtHat","canvInclusiveJetPtHat", 300, 800);
        canvLeadJetPtHat = new TCanvas("canvLeadJetPtHat","canvLeadJetPtHat", 300, 800);
        canvGenJetPtHat = new TCanvas("canvGenJetPtPtHat","canvGenJetPtPtHat", 300, 800);
    }
    canvMatchedJetPtHat->Divide( centLow.size(), 3 );
    canvUnmatchedJetPtHat->Divide( centLow.size(), 3 );
    canvInclusiveJetPtHat->Divide( centLow.size(), 3 );
    canvLeadJetPtHat->Divide( centLow.size(), 3 );
    canvVz->Divide( centLow.size(), 1 );
    canvUnmatchedJetPtForPtHat->Divide( centLow.size(), 1);
    canvGenJetPtHat->Divide( centLow.size(), 3 );

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TLegend *leg;   // Matched
    TLegend *leg2;  // Unmatched
    TLegend *leg3;  // pt for different ptHat
    TLegend *leg4;  // Leading jet
    TLegend *leg5;  // Inclusive
    TLegend *leg6;  // Gen jet

    TH2D* hMatchedPtHatVsPtReco[ centLow.size() ];
    TH2D* hUnmatchedPtHatVsPtReco[ centLow.size() ];
    TH2D* hInclusivePtHatVsPtReco[ centLow.size() ];
    TH2D* hMatchedLeadPtHatVsPtReco[ centLow.size() ];
    TH2D* hGenPtHatVsPtGen[ centLow.size() ];

    TH1D* hPtHatMatched[centLow.size()][jetPtLow.size()];
    TH1D* hPtHatUnmatched[centLow.size()][jetPtLow.size()];
    TH1D* hPtHatInclusive[centLow.size()][jetPtLow.size()];
    TH1D* hPtHatMatchedLead[centLow.size()][jetPtLow.size()];
    TH1D* hPtHatGen[centLow.size()][jetPtLow.size()];

    TGraph* hPtHatMatchedIntegral[centLow.size()][jetPtLow.size()];
    TGraph* hPtHatUnmatchedIntegral[centLow.size()][jetPtLow.size()];
    TGraph* hPtHatInclusiveIntegral[centLow.size()][jetPtLow.size()];
    TGraph* hPtHatMatchedLeadIntegral[centLow.size()][jetPtLow.size()];
    TGraph* hPtHatGenIntegral[centLow.size()][jetPtLow.size()];

    TMultiGraph *mgPtHatMatchedIntegral[centLow.size()];
    TMultiGraph *mgPtHatUnmatchedIntegral[centLow.size()];
    TMultiGraph *mgPtHatInclusiveIntegral[centLow.size()];
    TMultiGraph *mgPtHatMatchedLeadIntegral[centLow.size()];
    TMultiGraph *mgPtHatGenIntegral[centLow.size()];

    TH1D *hVz[ centLow.size() ];
    TH1D *hVzPtHatDep[ centLow.size() ][ ptHatLow.size() ]; 
    TH1D *hUnmatchedPtForPtHat[ centLow.size() ][ ptHatLow.size() ];

    // Loop over centralities
    for (Int_t i{0}; i<centLow.size(); i++) {
        
        // Select centrality range
        hRecoJetPtFlavPtHatCent->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoJetPtFlavPtHatCentInclusive->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoLeadJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hGenJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hVzPtHatCent->GetAxis(2)->SetRange( centLow.at(i), centHi.at(i) );

        // Vz
        hVz[i] = (TH1D*)hVzPtHatCent->Projection(0);
        hVz[i]->SetName( Form("hVz_%d",i) );
        setJetStyle( hVz[i], 0 );
        // Vz canvas for different centralities
        canvVz->cd( i + 1 );
        setPadStyle();
        hVz[i]->Draw();
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.55, 0.8, Form( "Integral: %.2f", hVz[i]->Integral() ) );

        //
        // Loop over ptHat bins and plot reco jet pt
        //
        for (Int_t j{0}; j<ptHatLow.size(); j++) {
            hVzPtHatCent->GetAxis(1)->SetRange( ptHatLow.at(j), ptHatHi.at(j) );
            hVzPtHatDep[i][j] = (TH1D*)hVzPtHatCent->Projection(0);
            hVzPtHatDep[i][j]->SetName( Form("hVzPtHatDep_%d_%d", i, j) );
            setJetStyle( hVzPtHatDep[i][j], 0 );

            // Plot Vz for each ptHat
            canvVzCentForPtHat[i]->cd( j + 1 );
            setPadStyle();
            hVzPtHatDep[i][j]->Draw();
            if ( isPbPb ) {
                t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
            }
            else {
                t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
            }
            t.DrawLatexNDC(0.32, 0.8, Form( "%d<#hat{p_{}} (GeV/c)<%d", 
                           ptHatLowVal + ptHatStep * (ptHatLow.at(j) - 1), 
                           ptHatLowVal + ptHatStep * (ptHatHi.at(j) ) ) );
            Double_t integral = hVzPtHatDep[i][j]->Integral();
            t.DrawLatexNDC(0.55, 0.7, Form( "Integral: %.2f", integral ) );

            // Set range of ptHat for unmatched jets
            hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(2)->SetRange( ptHatLow.at(j), ptHatHi.at(j) );
            hUnmatchedPtForPtHat[i][j] = (TH1D*)hRecoUnmatchedJetPtFlavPtHatCent->Projection(0);
            hUnmatchedPtForPtHat[i][j]->SetName( Form("hUnmatchedPtForPtHat_%d_%d", i, j) );
            setJetStyle( hUnmatchedPtForPtHat[i][j], j );
            hUnmatchedPtForPtHat[i][j]->Scale(1. / integral);

            // Plot jet pt
            canvUnmatchedJetPtForPtHat->cd( i + 1 );
            setPadStyle();
            hUnmatchedPtForPtHat[i][j]->Draw("same");
            hUnmatchedPtForPtHat[i][j]->GetYaxis()->SetRangeUser(1./10000, 10);
            gPad->SetLogy(1);

            if ( i == 0 && j == 0) {
                leg3 = new TLegend(0.45, 0.65, 0.9, 0.85);
            }
            if (i == 0) {
                leg3->AddEntry(hUnmatchedPtForPtHat[i][j], Form("%d<#hat{p_{T}} (GeV/c)<%d", 
                              ptHatLowVal + ptHatStep * (ptHatLow.at(j) - 1) , 
                              ptHatLowVal + ptHatStep * (ptHatHi.at(j) ) ), "p");
            }
            if (i == 0 && j == (ptHatLow.size() - 1) ) {
                leg3->SetTextSize(0.05);
                leg3->SetLineWidth(0);
                leg3->Draw();
            }

        } // for (Int_t j{0}; j<ptHatLow.size(); j++)

        // Restore ptHat projection back
        hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(2)->SetRange(0, 120);

        // Matched jets
        hMatchedPtHatVsPtReco[i] = (TH2D*)hRecoJetPtFlavPtHatCent->Projection(2, 0);
        hMatchedPtHatVsPtReco[i]->SetName(Form("hMatchedPtHatVsPtReco_%d",i));
        set2DStyle(hMatchedPtHatVsPtReco[i]);
        // Unmatched jets
        hUnmatchedPtHatVsPtReco[i] = (TH2D*)hRecoUnmatchedJetPtFlavPtHatCent->Projection(2, 0);
        hUnmatchedPtHatVsPtReco[i]->SetName(Form("hUnmatchedPtHatVsPtReco_%d",i));
        set2DStyle(hUnmatchedPtHatVsPtReco[i]);
        // Inclusive jets
        hInclusivePtHatVsPtReco[i] = (TH2D*)hRecoJetPtFlavPtHatCentInclusive->Projection(2, 0);
        hInclusivePtHatVsPtReco[i]->SetName(Form("hInclusivePtHatVsPtReco_%d",i));
        set2DStyle(hInclusivePtHatVsPtReco[i]);
        // Matched lead jets
        hMatchedLeadPtHatVsPtReco[i] = (TH2D*)hRecoLeadJetPtFlavPtHatCentWeighted->Projection(2, 0);
        hMatchedLeadPtHatVsPtReco[i]->SetName(Form("hMatchedLeadPtHatVsPtReco_%d",i));
        set2DStyle(hMatchedLeadPtHatVsPtReco[i]);
        // Generated jets
        hGenPtHatVsPtGen[i] = (TH2D*)hGenJetPtFlavPtHatCentWeighted->Projection(2, 0);
        hGenPtHatVsPtGen[i]->SetName(Form("hGenPtHatVsPtGen_%d",i));
        set2DStyle(hGenPtHatVsPtGen[i]);


        // Set multigraphs

        // Matched jets
        mgPtHatMatchedIntegral[i] = new TMultiGraph();
        mgPtHatMatchedIntegral[i]->SetNameTitle(Form("mgPtHatMatchedIntegral_%d",i),
                                                Form("mgPtHatMatchedIntegral_%d;#hat{p_{T}} (GeV/c);Fraction of total integral", i) );
        // Unmatched jets
        mgPtHatUnmatchedIntegral[i] = new TMultiGraph();
        mgPtHatUnmatchedIntegral[i]->SetNameTitle(Form("mgPtHatUnmatchedIntegral_%d",i),
                                                  Form("mgPtHatUnmatchedIntegral_%d;#hat{p_{T}} (GeV/c);Fraction of total integral", i) );
        // Inclusive jets
        mgPtHatInclusiveIntegral[i] = new TMultiGraph();
        mgPtHatInclusiveIntegral[i]->SetNameTitle(Form("mgPtHatInclusiveIntegral_%d",i),
                                                  Form("mgPtHatInclusiveIntegral_%d;#hat{p_{T}} (GeV/c);Fraction of total integral", i) );
        // Matched lead jets
        mgPtHatMatchedLeadIntegral[i] = new TMultiGraph();
        mgPtHatMatchedLeadIntegral[i]->SetNameTitle(Form("mgPtHatMatchedLeadIntegral_%d",i),
                                                    Form("mgPtHatMatchedLeadIntegral_%d;#hat{p_{T}} (GeV/c);Fraction of total integral",i));

        // Gen jets
        mgPtHatGenIntegral[i] = new TMultiGraph();
        mgPtHatGenIntegral[i]->SetNameTitle(Form("mgPtHatGenIntegral_%d",i),
                                            Form("mgPtHatGenIntegral_%d;#hat{p_{T}} (GeV/c);Fraction of total integral",i));


        //
        // Loop over reconstructed jet transverse momentum
        //
        for (Int_t j{0}; j<jetPtLow.size(); j++) {

            //
            // Matched jets
            //
            //Int_t ptBins[2] = {1, 50};
            hPtHatMatched[i][j] = (TH1D*)hMatchedPtHatVsPtReco[i]->ProjectionY(Form("hPtHatMatched_%d_%d",i,j), jetPtLow.at(j),jetPtHi.at(j));
            hPtHatMatched[i][j]->SetName( Form("hPtHatMatched_%d_%d", i, j) );
            setJetStyle(hPtHatMatched[i][j], j);
            //hPtHatMatched[i][j]->Scale(1./hPtHatMatched[i][j]->Integral(ptBins[0], ptBins[1]));
            hPtHatMatchedIntegral[i][j] = new TGraph();
            hPtHatMatchedIntegral[i][j]->SetNameTitle( Form("hPtHatMatchedIntegral_%d_%d", i, j),
                                                       Form("hPtHatMatchedIntegral_%d_%d;#hat{p_{T}} (GeV/c);Fraction of total integral",i,j) );
            calculateCumulativeIntegral(hPtHatMatchedIntegral[i][j], hPtHatMatched[i][j], j);
            mgPtHatMatchedIntegral[i]->Add(hPtHatMatchedIntegral[i][j]);

            //
            // Unmatched jets
            //
            hPtHatUnmatched[i][j] = (TH1D*)hUnmatchedPtHatVsPtReco[i]->ProjectionY(Form("hPtHatUnmatched_%d_%d",i,j), jetPtLow.at(j),jetPtHi.at(j));
            hPtHatUnmatched[i][j]->SetName( Form("hPtHatUnmatched_%d_%d", i, j) );
            setJetStyle(hPtHatUnmatched[i][j], j);
            //hPtHatUnmatched[i][j]->Scale(1./hPtHatUnmatched[i][j]->Integral(ptBins[0], ptBins[1]));
            hPtHatUnmatchedIntegral[i][j] = new TGraph();
            hPtHatUnmatchedIntegral[i][j]->SetNameTitle( Form("hPtHatUnmatchedIntegral_%d_%d",i,j),
                                                         Form("hPtHatUnmatchedIntegral_%d_%d;#hat{p_{T}} (GeV/c);Fraction of total integral",i,j) );
            calculateCumulativeIntegral(hPtHatUnmatchedIntegral[i][j], hPtHatUnmatched[i][j], j);
            mgPtHatUnmatchedIntegral[i]->Add(hPtHatUnmatchedIntegral[i][j]);

            //
            // Inclusive jets
            //
            hPtHatInclusive[i][j] = (TH1D*)hInclusivePtHatVsPtReco[i]->ProjectionY(Form("hPtHatInclusive_%d_%d",i,j), jetPtLow.at(j),jetPtHi.at(j));
            hPtHatInclusive[i][j]->SetName( Form("hPtHatInclusive_%d_%d",i,j) );
            setJetStyle(hPtHatInclusive[i][j], j);
            //hPtHatInclusive[i][j]->Scale(1./hPtHatInclusive[i][j]->Integral(ptBins[0], ptBins[1]));
            hPtHatInclusiveIntegral[i][j] = new TGraph();
            hPtHatInclusiveIntegral[i][j]->SetNameTitle( Form("hPtHatInclusiveIntegral_%d_%d", i, j),
                                                         Form("hPtHatInclusiveIntegral_%d_%d;#hat{p_{T}} (GeV/c);Fraction of total integral", i, j) );
            calculateCumulativeIntegral(hPtHatInclusiveIntegral[i][j], hPtHatInclusive[i][j], j);
            mgPtHatInclusiveIntegral[i]->Add(hPtHatInclusiveIntegral[i][j]);

            //
            // Matched leading jets
            //
            hPtHatMatchedLead[i][j] = (TH1D*)hMatchedLeadPtHatVsPtReco[i]->ProjectionY(Form("hPtHatMatchedLead_%d_%d",i,j), jetPtLow.at(j),jetPtHi.at(j));
            hPtHatMatchedLead[i][j]->SetName( Form("hPtHatMatchedLead_%d_%d",i,j) );
            setJetStyle(hPtHatMatchedLead[i][j], j);
            hPtHatMatchedLeadIntegral[i][j] = new TGraph();
            hPtHatMatchedLeadIntegral[i][j]->SetNameTitle( Form("hPtHatMatchedLeadIntegral_%d_%d",i,j),
                                                           Form("hPtHatMatchedLeadIntegral_%d_%d;#hat{p_{T}} (GeV/c);Fraction of total integral", i, j) );
            calculateCumulativeIntegral(hPtHatMatchedLeadIntegral[i][j], hPtHatMatchedLead[i][j], j);
            mgPtHatMatchedLeadIntegral[i]->Add( hPtHatMatchedLeadIntegral[i][j] );

            //
            // Gen jets
            //
            hPtHatGen[i][j] = (TH1D*)hGenPtHatVsPtGen[i]->ProjectionY(Form("hPtHatGen_%d_%d",i,j), jetPtLow.at(j), jetPtHi.at(j));
            hPtHatGen[i][j]->SetName( Form("hPtHatGen_%d_%d",i,j) );
            setJetStyle(hPtHatGen[i][j], j);
            hPtHatGenIntegral[i][j] = new TGraph();
            hPtHatGenIntegral[i][j]->SetNameTitle( Form("hPtHatGenIntegral_%d_%d",i,j),
                                                   Form("hPtHatGenIntegral_%d_%d;#hat{p_{T}} (GeV/c);Fraction of total integral", i, j) );
            calculateCumulativeIntegral(hPtHatGenIntegral[i][j], hPtHatGen[i][j], j);
            mgPtHatGenIntegral[i]->Add( hPtHatGenIntegral[i][j] );
        }


        //
        // Matched jets
        //

        // ptHat vs. matched jet pT
        canvMatchedJetPtHat->cd( i + 1 );
        setPadStyle();
        hMatchedPtHatVsPtReco[i]->Draw("colz");
        gPad->SetLogz(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Matched jets"));

        // ptHat projections for different jet pT
        canvMatchedJetPtHat->cd( centLow.size() + i + 1 );
        setPadStyle();
        if (i == 0) {
            leg = new TLegend(0.5, 0.65, 0.9, 0.85);
        }
        for (Int_t j{0}; j<jetPtLow.size(); j++) {
            hPtHatMatched[i][j]->Draw("same");
            if (i == 0) {
                leg->AddEntry(hPtHatMatched[i][j], Form("%d<p_{T} (GeV/c)<%d", 
                              ptLowVal + ptStep * (jetPtLow.at(j)-1) , 
                              ptLowVal + ptStep * ( jetPtHi.at(j) )), "p");
            }
        }
        gPad->SetLogy(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Matched jets"));
     
        if ( i == 0 ) {
            leg->SetLineWidth(0);
            leg->SetTextSize(0.05);
            leg->Draw();
        }

        // Cumulative integrals
        canvMatchedJetPtHat->cd( 2 * centLow.size() + i + 1 );
        setPadStyle();
        mgPtHatMatchedIntegral[i]->Draw("AP");
        mgPtHatMatchedIntegral[i]->GetYaxis()->SetTitleSize(0.06);
        mgPtHatMatchedIntegral[i]->GetYaxis()->SetLabelSize(0.06);
        mgPtHatMatchedIntegral[i]->GetXaxis()->SetTitleSize(0.06);
        mgPtHatMatchedIntegral[i]->GetXaxis()->SetLabelSize(0.06);
        mgPtHatMatchedIntegral[i]->GetXaxis()->SetNdivisions(208);
        mgPtHatMatchedIntegral[i]->GetYaxis()->SetNdivisions(208);    
        mgPtHatMatchedIntegral[i]->GetYaxis()->SetTitleOffset(1.0);

        // Zoom in
        // mgPtHatMatchedIntegral[i]->GetYaxis()->SetRangeUser(0.84, 1.01);
        // mgPtHatMatchedIntegral[i]->GetXaxis()->SetRangeUser(65., 155.);
        mgPtHatMatchedIntegral[i]->GetYaxis()->SetRangeUser(0., 0.35);
        mgPtHatMatchedIntegral[i]->GetXaxis()->SetRangeUser(15., 75.);
        gPad->SetGridx(1);
        gPad->SetGridy(1);

        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Matched jets"));


        //
        // Unmatched jets
        //

        // ptHat vs. unmatched jet pT
        canvUnmatchedJetPtHat->cd( i + 1 );
        setPadStyle();
        hUnmatchedPtHatVsPtReco[i]->Draw("colz");
        gPad->SetLogz(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Unmatched jets"));

        // ptHat projections for different jet pT
        canvUnmatchedJetPtHat->cd( centLow.size() + i + 1 );
        setPadStyle();
        for (Int_t j{0}; j<jetPtLow.size(); j++) {
            hPtHatUnmatched[i][j]->Draw("same");
        }
        gPad->SetLogy(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Unmatched jets"));

        if (i == 0) {
            leg2 = new TLegend(0.5, 0.4, 0.9, 0.65);
            for (Int_t j{0}; j<jetPtLow.size(); j++) {
                leg2->AddEntry(hPtHatUnmatched[i][j], Form("%d<p_{T} (GeV/c)<%d", 
                              ptLowVal + ptStep * (jetPtLow.at(j)-1) , 
                              ptLowVal + ptStep * (jetPtHi.at(j))), "p");
            }
            leg2->SetLineWidth(0);
            leg2->SetTextSize(0.05);
            leg2->Draw();
        }

        // Cumulative integrals
        canvUnmatchedJetPtHat->cd( 2 * centLow.size() + i + 1 );
        setPadStyle();
        mgPtHatUnmatchedIntegral[i]->Draw("AP");
        mgPtHatUnmatchedIntegral[i]->GetYaxis()->SetTitleSize(0.06);
        mgPtHatUnmatchedIntegral[i]->GetYaxis()->SetLabelSize(0.06);
        mgPtHatUnmatchedIntegral[i]->GetXaxis()->SetTitleSize(0.06);
        mgPtHatUnmatchedIntegral[i]->GetXaxis()->SetLabelSize(0.06);
        mgPtHatUnmatchedIntegral[i]->GetXaxis()->SetNdivisions(208);
        mgPtHatUnmatchedIntegral[i]->GetYaxis()->SetNdivisions(208);    
        mgPtHatUnmatchedIntegral[i]->GetYaxis()->SetTitleOffset(1.0);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.55, 0.3, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.55, 0.3, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.5, 0.2, Form("Unmatched jets"));


        //
        // Inclusive jets
        //

        // ptHat vs. inclusive jet pT
        canvInclusiveJetPtHat->cd( i + 1 );
        setPadStyle();
        hInclusivePtHatVsPtReco[i]->Draw("colz");
        gPad->SetLogz(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Inclusive jets"));

        // ptHat projections for different jet pT
        canvInclusiveJetPtHat->cd( centLow.size() + i + 1 );
        setPadStyle();
        if (i == 0) {
            leg5 = new TLegend(0.5, 0.65, 0.9, 0.85);
        }
        for (Int_t j{0}; j<jetPtLow.size(); j++) {
            hPtHatInclusive[i][j]->Draw("same");
            if (i == 0) {
                leg5->AddEntry(hPtHatInclusive[i][j], Form("%d<p_{T} (GeV/c)<%d", 
                               ptLowVal + ptStep * (jetPtLow.at(j)-1) , 
                               ptLowVal + ptStep * ( jetPtHi.at(j) )), "p");
            }
        }
        gPad->SetLogy(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Inclusive jets"));

        if ( i == 0 ) {
            leg5->SetLineWidth(0);
            leg5->SetTextSize(0.05);
            leg5->Draw();
        }

        // Cumulative integrals
        canvInclusiveJetPtHat->cd( 2 * centLow.size() + i + 1 );
        setPadStyle();
        mgPtHatInclusiveIntegral[i]->Draw("AP");
        mgPtHatInclusiveIntegral[i]->GetYaxis()->SetTitleSize(0.06);
        mgPtHatInclusiveIntegral[i]->GetYaxis()->SetLabelSize(0.06);
        mgPtHatInclusiveIntegral[i]->GetXaxis()->SetTitleSize(0.06);
        mgPtHatInclusiveIntegral[i]->GetXaxis()->SetLabelSize(0.06);
        mgPtHatInclusiveIntegral[i]->GetXaxis()->SetNdivisions(208);
        mgPtHatInclusiveIntegral[i]->GetYaxis()->SetNdivisions(208);    
        mgPtHatInclusiveIntegral[i]->GetYaxis()->SetTitleOffset(1.0);


        // Zoom in
        // mgPtHatInclusiveIntegral[i]->GetYaxis()->SetRangeUser(0.84, 1.01);
        // mgPtHatInclusiveIntegral[i]->GetXaxis()->SetRangeUser(65., 155.);

        mgPtHatInclusiveIntegral[i]->GetYaxis()->SetRangeUser(0., 0.35);
        mgPtHatInclusiveIntegral[i]->GetXaxis()->SetRangeUser(15., 75.);
        gPad->SetGridx(1);
        gPad->SetGridy(1);

        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Inclusive jets"));


        //
        // Leading jets 
        // 

        // ptHat vs. leading jet pT
        canvLeadJetPtHat->cd( i + 1 );
        setPadStyle();
        hMatchedLeadPtHatVsPtReco[i]->Draw("colz");
        gPad->SetLogz(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Leading jets (matched)"));

        // ptHat for different jet pT
        canvLeadJetPtHat->cd( centLow.size() + i + 1 );
        setPadStyle();
        if (i == 0) {
            leg4 = new TLegend(0.5, 0.65, 0.9, 0.85);
        }
        for (Int_t j{0}; j<jetPtLow.size(); j++) {
            hPtHatMatchedLead[i][j]->Draw("same");
            if (i == 0) {
                leg4->AddEntry(hPtHatMatchedLead[i][j], Form("%d<p_{T} (GeV/c)<%d", 
                              ptLowVal + ptStep * (jetPtLow.at(j)-1) , 
                              ptLowVal + ptStep * ( jetPtHi.at(j) )), "p");
            }
        }
        gPad->SetLogy(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Leading jets (matched)"));
     
        if ( i == 0 ) {
            leg4->SetLineWidth(0);
            leg4->SetTextSize(0.05);
            leg4->Draw();
        }

        // Leading jet pT cumulative integrals
        canvLeadJetPtHat->cd( 2 * centLow.size() + i + 1 );
        setPadStyle();
        mgPtHatMatchedLeadIntegral[i]->Draw("AP");
        mgPtHatMatchedLeadIntegral[i]->GetYaxis()->SetTitleSize(0.06);
        mgPtHatMatchedLeadIntegral[i]->GetYaxis()->SetLabelSize(0.06);
        mgPtHatMatchedLeadIntegral[i]->GetXaxis()->SetTitleSize(0.06);
        mgPtHatMatchedLeadIntegral[i]->GetXaxis()->SetLabelSize(0.06);
        mgPtHatMatchedLeadIntegral[i]->GetXaxis()->SetNdivisions(208);
        mgPtHatMatchedLeadIntegral[i]->GetYaxis()->SetNdivisions(208);    
        mgPtHatMatchedLeadIntegral[i]->GetYaxis()->SetTitleOffset(1.0);

        // Zoom in
        // mgPtHatMatchedLeadIntegral[i]->GetYaxis()->SetRangeUser(0.84, 1.01);
        // mgPtHatMatchedLeadIntegral[i]->GetXaxis()->SetRangeUser(65., 155.);

        mgPtHatMatchedLeadIntegral[i]->GetYaxis()->SetRangeUser(0., 0.35);
        mgPtHatMatchedLeadIntegral[i]->GetXaxis()->SetRangeUser(15., 75.);

        gPad->SetGridx(1);
        gPad->SetGridy(1);

        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Leading jets (matched)"));



        //
        // Gen jets 
        // 

        // ptHat vs. gen jet pT
        canvGenJetPtHat->cd( i + 1 );
        setPadStyle();
        hGenPtHatVsPtGen[i]->Draw("colz");
        gPad->SetLogz(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Gen jets"));

        // ptHat for different jet pT
        canvGenJetPtHat->cd( centLow.size() + i + 1 );
        setPadStyle();
        if (i == 0) {
            leg4 = new TLegend(0.5, 0.65, 0.9, 0.85);
        }
        for (Int_t j{0}; j<jetPtLow.size(); j++) {
            hPtHatGen[i][j]->Draw("same");
            if (i == 0) {
                leg4->AddEntry(hPtHatGen[i][j], Form("%d<p_{T} (GeV/c)<%d", 
                              ptLowVal + ptStep * (jetPtLow.at(j)-1) , 
                              ptLowVal + ptStep * ( jetPtHi.at(j) )), "p");
            }
        }
        gPad->SetLogy(1);
        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Gen jets"));
     
        if ( i == 0 ) {
            leg4->SetLineWidth(0);
            leg4->SetTextSize(0.05);
            leg4->Draw();
        }

        // Gen jet pT cumulative integrals
        canvGenJetPtHat->cd( 2 * centLow.size() + i + 1 );
        setPadStyle();
        mgPtHatGenIntegral[i]->Draw("AP");
        mgPtHatGenIntegral[i]->GetYaxis()->SetTitleSize(0.06);
        mgPtHatGenIntegral[i]->GetYaxis()->SetLabelSize(0.06);
        mgPtHatGenIntegral[i]->GetXaxis()->SetTitleSize(0.06);
        mgPtHatGenIntegral[i]->GetXaxis()->SetLabelSize(0.06);
        mgPtHatGenIntegral[i]->GetXaxis()->SetNdivisions(208);
        mgPtHatGenIntegral[i]->GetYaxis()->SetNdivisions(208);    
        mgPtHatGenIntegral[i]->GetYaxis()->SetTitleOffset(1.0);

        // Zoom in
        //mgPtHatGenIntegral[i]->GetYaxis()->SetRangeUser(0.84, 1.01);
        //mgPtHatGenIntegral[i]->GetXaxis()->SetRangeUser(65., 155.);

        mgPtHatGenIntegral[i]->GetYaxis()->SetRangeUser(0.0, 0.35);
        mgPtHatGenIntegral[i]->GetXaxis()->SetRangeUser(15., 75.);

        gPad->SetGridx(1);
        gPad->SetGridy(1);

        if ( isPbPb ) {
            t.DrawLatexNDC(0.32, 0.9, Form("%d-%d%% PYTHIA+HYDJET",(centLow.at(i)-2)*10, (centHi.at(i)-1)*10));
        }
        else {
            t.DrawLatexNDC(0.32, 0.9, Form("PYTHIA"));
        }
        t.DrawLatexNDC(0.32, 0.8, Form("Gen jets"));



    } // for (Int_t i{0}; i<centLow.size(); i++)

}

//________________
void drawJetAnalysisFigs(const Char_t *inFileNamePbPb = "../build/../build/oTestReadForest_PbPb_mcHiBinShift_centWeight_leadJetSelection.root",
                         const Char_t *inFileNamePP = "../build/oTestReadForest_pp_noPtHatCut.root", 
                         const Char_t *oFile = "oDrawJetAna.root") {

    // "../build/oTestReadForest_PbPb_ptHat50.root"
    // "../build/oTestReadForest_PbPb_ptHat50_jetpt50.root"
    // "../build/oTestReadForest_PbPb_noPtHatCut.root"
    // "../build/oTestReadForest_PbPb_mcHiBinShift_centWeight_leadJetSelection.root"

    // "../build/oTestReadForest_pp_ptHat50.root"
    // "../build/oTestReadForest_pp_noPtHatCut.root"

    //setStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Read ROOT file                             
    TFile *inFilePbPb = TFile::Open(inFileNamePbPb);
    if ( !inFilePbPb->IsOpen() ) {
        std::cout << "Input file for PbPb not opened. Terminating." << std::endl;
        exit(0);
    }

    // Read ROOT file                             
    TFile *inFilePP = TFile::Open(inFileNamePP);
    if ( !inFilePP->IsOpen() ) {
        std::cout << "Input file for pp not opened. Terminating." << std::endl;
        exit(0);
    }

    // Draw event histograms
    //drawEventQuantities(inFilePbPb);

    // Draw number of jets
    //drawNumberOfJets(inFilePbPb);

    //drawReco2GenSpectraPtHatComparison(inFilePbPb);

    // Draw gen pt distribution to find pthat
    //drawFindPtHatCut(inFilePP);

    // Draw correlation between pT corr vs. raw, corr vs. gen
    //drawPtRawCorrGen(inFilePbPb);

    // Draw correlation between ptHat and pt of leading jet with and without weight
    //drawPtHatVsLeadingJet(inFilePbPb);

    //
    // Draw JES and JER as a function of centrality
    //
    //drawJESvsCentrality(inFilePbPb);

    // Draw comparison of unmatched to inclusive reco jets
    //drawRecoUnmatched2Inclusive(inFilePbPb);

    // Draw scaled unmatched jets
    //drawRecoUnmatched(inFilePbPb);

    //
    // Compare ptHat to matched, unmatched and inclusive reco jets
    // Extract information about the pThat and xJets from HYDJET
    //
    drawPtHatVsRecoJets(inFilePbPb, kTRUE);

    //
    // For the events with ptHat>x compare gen jet pT spectra
    //
    //drawCompareGenPP2PbPb(inFilePbPb, inFilePP);

    //
    // Draw JES and JER as a function of pseudorapidity
    //
    //drawJESvsEta(inFilePbPb);
}