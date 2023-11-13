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

    Double_t ptHatStep = 20.;
    vector<Int_t> ptHatLow{0, 1, 2, 3, 4};
    vector<Int_t> ptHatHi {9, 9, 9, 9, 9};
    Int_t ptHatBin = 0;
    Double_t ptHatVal[2] = {15. + ptHatLow.at(ptHatBin), 215.};

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
    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {
        centLow.at(iCent) += firstHistBin;
        centHi.at(iCent) +=  firstHistBin;
    }

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
    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {
        centLow.at(iCent) += firstHistBin;
        centHi.at(iCent) +=  firstHistBin;
    }

    // PtHat projections starting 15 GeV with 20 GeV step
    Int_t ptHatMin = 15;
    Int_t pHatStep = 20;
    vector<Int_t> ptHatLow{0, 1, 2, 3, 4};
    vector<Int_t> ptHatHi{29, 29, 29, 29, 29};
    Int_t ptHatBins = ptHatLow.size();
    for (Int_t i{0}; i<ptHatLow.size(); i++) {
        ptHatLow.at(i) += firstHistBin;
        ptHatHi.at(i) += firstHistBin;
    }

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
    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {
        centLow.at(iCent) += firstHistBin;
        centHi.at(iCent) +=  firstHistBin;
    }

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
        hRecoJetRawPtCorrPtGenPtCent->GetAxis(3)->SetRange(centLow.at(i), centHi.at(i));

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

    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    Int_t ptBinLow{6};
    Int_t ptBinHi{50};
    Int_t etaBinLow{10};
    Int_t etaBinHi{40};
    for (Int_t i{0}; i<centLow.size(); i++) {
        centLow.at(i) += firstHistBin;
        centHi.at(i) += firstHistBin;
    }
    ptBinLow += firstHistBin;
    ptBinHi += firstHistBin;
    etaBinLow += firstHistBin;
    etaBinHi  += firstHistBin;

    // Set common pt and eta limits
    hGenJetPtEtaPhiCentPbPb->GetAxis(0)->SetRange(ptBinLow, ptBinHi);
    hGenJetPtEtaPhiCentPbPb->GetAxis(1)->SetRange(etaBinLow, etaBinHi);

    hGenJetPtEtaPhiCentPP->GetAxis(0)->SetRange(ptBinLow, ptBinHi);
    hGenJetPtEtaPhiCentPP->GetAxis(1)->SetRange(etaBinLow, etaBinHi);

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
    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {
        centLow.at(iCent) += firstHistBin;
        centHi.at(iCent) +=  firstHistBin;
    }
    Int_t ptLowBin{10};
    Int_t ptHiBin{50};

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
    Int_t rebinY = 2;

    TLegend *leg;
    // Centrality dependence
    for (Int_t iCent{0}; iCent<centLow.size(); iCent++) {

        //std::cout << "iCent: " << iCent << " padIndex: " << padIndex << " centHistoBin: " << centHistoBin << std::endl;
        
        // JES
        hJESPtEtaPhiCent->GetAxis(4)->SetRange(centLow.at(iCent), centHi.at(iCent));
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
        hJESPtEtaPhiCentWeighted->GetAxis(4)->SetRange(centLow.at(iCent), centHi.at(iCent));
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
    THnSparseD *hJESRecoEtaPhiPtHat = (THnSparseD*)inFile->Get("hJESRecoEtaPhiPtHat");
    THnSparseD *hJESRecoEtaPhiPtHatWeighted = (THnSparseD*)inFile->Get("hJESRecoEtaPhiPtHatWeighted");

    const Int_t etaBins = 5;

    // Histograms to fill
    TH2D *hJESVsGenPtEtaDep[etaBins];
    TH1D *hJESMeanEtaDep[etaBins];
    TH1D *hJESSigmaEtaDep[etaBins];
    TH2D *hJESVsGenPtEtaDepW[etaBins];
    TH1D *hJESMeanEtaDepW[etaBins];
    TH1D *hJESSigmaEtaDepW[etaBins];
    TCanvas *canvJESeta = new TCanvas("canvJESeta","canvJESeta", 1400, 800);
    canvJESeta->Divide(etaBins, 3);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);


    // Eta dependence
    for (Int_t iEta=0; iEta<etaBins; iEta++) {

        Int_t padIndex = iEta + 1;
        Double_t etaHi = 2.5;
        Double_t etaLow = -2.5;
        Double_t etaStep = (etaHi - etaLow) / etaBins;

        //std::cout << "iEta: " << iEta << " padIndex: " << padIndex << std::endl;

        // JES
        hJESRecoEtaPhiPtHat->GetAxis(2)->SetRange(iEta, iEta);
        hJESVsGenPtEtaDep[iEta] = (TH2D*)hJESRecoEtaPhiPtHat->Projection(0,1);
        hJESVsGenPtEtaDep[iEta]->SetName(Form("hJESVsGenPtEtaDep_%d",iEta));
        hJESVsGenPtEtaDep[iEta]->RebinX(6);
        hJESVsGenPtEtaDep[iEta]->RebinY(4);

        hJESRecoEtaPhiPtHatWeighted->GetAxis(2)->SetRange(iEta, iEta);
        hJESVsGenPtEtaDepW[iEta] = (TH2D*)hJESRecoEtaPhiPtHatWeighted->Projection(0,1);
        hJESVsGenPtEtaDepW[iEta]->SetName(Form("hJESVsGenPtEtaDepW_%d",iEta));
        hJESVsGenPtEtaDepW[iEta]->RebinX(6);
        hJESVsGenPtEtaDepW[iEta]->RebinY(4);

        
        hJESVsGenPtEtaDep[iEta]->FitSlicesY();
        hJESMeanEtaDep[iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDep_%d_%d", iEta, 1));
        hJESMeanEtaDep[iEta]->SetName(Form("hJESMeanEtaDep_%d", iEta));
        hJESSigmaEtaDep[iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDep_%d_%d", iEta, 2));
        hJESSigmaEtaDep[iEta]->SetName(Form("hJESSigmaEtaDep_%d",iEta));

        hJESVsGenPtEtaDepW[iEta]->FitSlicesY();
        hJESMeanEtaDepW[iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDepW_%d_%d", iEta, 1));
        hJESMeanEtaDepW[iEta]->SetName(Form("hJESMeanEtaDepW_%d", iEta));
        hJESSigmaEtaDepW[iEta] = (TH1D*)gDirectory->Get(Form("hJESVsGenPtEtaDepW_%d_%d", iEta, 2));
        hJESSigmaEtaDepW[iEta]->SetName(Form("hJESSigmaEtaDepW_%d",iEta));

        // JES (2D)
        canvJESeta->cd(padIndex);
        setPadStyle();
        gPad->SetLeftMargin(0.17);
        set2DStyle(hJESVsGenPtEtaDep[iEta]);
        hJESVsGenPtEtaDep[iEta]->Draw("colz");
        hJESVsGenPtEtaDep[iEta]->GetYaxis()->SetRangeUser(-1., 2.);
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.25, 0.9, Form("JES: %2.1f<#eta<%2.1f Pb+Pb #sqrt{s_{NN}}=5020 GeV",etaLow+iEta*etaStep, etaLow+(iEta+1)*etaStep));

        
        // Jet energy scale
        canvJESeta->cd(padIndex + etaBins);
        setPadStyle();
        set1DStyle(hJESMeanEtaDep[iEta],0);
        hJESMeanEtaDep[iEta]->Draw();
        hJESMeanEtaDep[iEta]->GetYaxis()->SetRangeUser(0.9, 1.3);
        hJESMeanEtaDep[iEta]->GetYaxis()->SetTitle("#mu(JES)");
        set1DStyle(hJESMeanEtaDepW[iEta],1);
        hJESMeanEtaDepW[iEta]->Draw("same");

        // Jet energy resolution
        canvJESeta->cd(padIndex + 2 * etaBins);
        setPadStyle();
        set1DStyle(hJESSigmaEtaDep[iEta],0);
        hJESSigmaEtaDep[iEta]->Draw();
        hJESSigmaEtaDep[iEta]->GetYaxis()->SetRangeUser(0., 0.35);
        hJESSigmaEtaDep[iEta]->GetYaxis()->SetTitle("#sigma(JES)");
        set1DStyle(hJESSigmaEtaDepW[iEta],1);
        hJESSigmaEtaDepW[iEta]->Draw("same");
    } // for (Int_t iEta=0; iEta<etaBins; iEta++)
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
    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    for (Int_t i{0}; i<centLow.size(); i++) {
        centLow.at(i) += firstHistBin;
        centHi.at(i) += firstHistBin;
    }

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
    Int_t ptHatStep = 20;
    vector<Int_t> ptHatLow{0, 1, 2, 3};
    vector<Int_t> ptHatHi{29, 29, 29, 29};
    for (Int_t i{0}; i<ptHatLow.size(); i++) {
        ptHatLow.at(i) += firstHistBin;
        ptHatHi.at(i) += firstHistBin;
    }

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

    auto hVzCent = (THnSparseD*)inFile->Get("hVzCent");
    auto hVzCentWeighted = (THnSparseD*)inFile->Get("hVzCentWeighted");
    auto hRecoUnmatchedJetPtFlavPtHatCent = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCent");
    auto hRecoUnmatchedJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoUnmatchedJetPtFlavPtHatCentWeighted");
    // Matched
    auto hRecoJetPtFlavPtHatCentWeighted = (THnSparseD*)inFile->Get("hRecoJetPtFlavPtHatCentWeighted");

    // Setup centralities to read
    Int_t firstHistBin = 1;
    vector<Int_t> centLow{6, 4, 2, 1};
    vector<Int_t> centHi {9, 5, 3, 1};
    for (Int_t i{0}; i<centLow.size(); i++) {
        centLow.at(i) += firstHistBin;
        centHi.at(i) += firstHistBin;
    }

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

    // Loop over centralities
    for (Int_t i{0}; i<centLow.size(); i++) {
        
        // Select centrality range
        hVzCent->GetAxis(1)->SetRange( centLow.at(i), centHi.at(i) );
        hVzCentWeighted->GetAxis(1)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoUnmatchedJetPtFlavPtHatCent->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoUnmatchedJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        hRecoJetPtFlavPtHatCentWeighted->GetAxis(3)->SetRange( centLow.at(i), centHi.at(i) );
        
        // Event-wise histograms and integrals
        hVz[i] = (TH1D*)hVzCent->Projection(0);
        hVz[i]->SetName(Form("hVz_%d",i));
        hVzW[i] = (TH1D*)hVzCentWeighted->Projection(0);
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
        setJetStyle(hPtUnmatched2MatchedW[i], 4);
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
void drawJetAnalysisFigs(const Char_t *inFileNamePbPb = "../build/oTestReadForest_PbPb_noPtHatCut.root",
                         const Char_t *inFileNamePP = "../build/oTestReadForest_pp_ptHat50.root", 
                         const Char_t *oFile = "oDrawJetAna.root") {

    // "../build/oTestReadForest_PbPb_ptHat50.root"
    // "../build/oTestReadForest_PbPb_noPtHatCut.root"

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

    // Draw JES and JER as a function of centrality
    //drawJESvsCentrality(inFilePbPb);

    // Draw comparison of unmatched to inclusive reco jets
    //drawRecoUnmatched2Inclusive(inFilePbPb);

    drawRecoUnmatched(inFilePbPb);

    // For the events with ptHat>x compare gen jet pT spectra
    //drawCompareGenPP2PbPb(inFilePbPb, inFilePP);

    // Draw JES and JER as a function of pseudorapidity
    //drawJESvsEta(inFile);
}