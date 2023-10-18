#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"

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
    gPad->SetRightMargin(0.);
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
    h->GetYaxis()->SetTitleOffset(0.8);
}

//________________
void set2DStyle(TH2* h) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(0.8);
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
void drawJESvsCentrality(TFile *inFile) {

    // Jet histograms
    THnSparseD *hJESReco = (THnSparseD*)inFile->Get("hJESReco");
    THnSparseD *hJESRecoWeighted = (THnSparseD*)inFile->Get("hJESRecoWeighted");

    const Int_t centBins = 9; // 0-10%, 10-20%, ...

    TH2D *hJESvsGenPt[centBins];
    TH1D *hJESMean[centBins];
    TH1D *hJESSigma[centBins];
    TH2D *hJESvsGenPtW[centBins];
    TH1D *hJESMeanW[centBins];
    TH1D *hJESSigmaW[centBins];
    TCanvas *canvJES = new TCanvas("canvJES","canvJES", 1400, 800);
    canvJES->Divide(centBins, 3);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Centrality dependence
    for (Int_t iCent=centBins; iCent>0; iCent--) {

        Int_t padIndex = centBins - iCent + 1;
        Int_t centHistoBin = centBins - iCent;

        //std::cout << "iCent: " << iCent << " padIndex: " << padIndex << " centHistoBin: " << centHistoBin << std::endl;
        
        // JES
        hJESReco->GetAxis(3)->SetRange(iCent-1, iCent-1);
        hJESvsGenPt[centHistoBin] = (TH2D*)hJESReco->Projection(0,1);
        hJESvsGenPt[centHistoBin]->SetName(Form("hJESvsGenPt_%d",centHistoBin));
        hJESvsGenPt[centHistoBin]->RebinX(6);
        hJESvsGenPt[centHistoBin]->RebinY(4);

        hJESvsGenPt[centHistoBin]->FitSlicesY();
        hJESMean[centHistoBin] = (TH1D*)gDirectory->Get(Form("hJESvsGenPt_%d_%d", centHistoBin, 1));
        hJESMean[centHistoBin]->SetName(Form("hJESMean_%d", centHistoBin));
        hJESSigma[centHistoBin] = (TH1D*)gDirectory->Get(Form("hJESvsGenPt_%d_%d", centHistoBin, 2));
        hJESSigma[centHistoBin]->SetName(Form("hJESSigma_%d",centHistoBin));

        // JES ptHat weighted
        hJESRecoWeighted->GetAxis(3)->SetRange(iCent-1, iCent-1);
        hJESvsGenPtW[centHistoBin] = (TH2D*)hJESRecoWeighted->Projection(0,1);
        hJESvsGenPtW[centHistoBin]->SetName(Form("hJESvsGenPtW_%d",centHistoBin));
        hJESvsGenPtW[centHistoBin]->RebinX(5);
        hJESvsGenPtW[centHistoBin]->RebinY(4);

        hJESvsGenPtW[centHistoBin]->FitSlicesY();
        hJESMeanW[centHistoBin] = (TH1D*)gDirectory->Get(Form("hJESvsGenPtW_%d_%d", centHistoBin, 1));
        hJESMeanW[centHistoBin]->SetName(Form("hJESMeanW_%d", centHistoBin));
        hJESSigmaW[centHistoBin] = (TH1D*)gDirectory->Get(Form("hJESvsGenPtW_%d_%d", centHistoBin, 2));
        hJESSigmaW[centHistoBin]->SetName(Form("hJESSigmaW_%d",centHistoBin));

        // JES (2D)
        canvJES->cd(padIndex);
        setPadStyle();
        gPad->SetLeftMargin(0.17);
        set2DStyle(hJESvsGenPt[centHistoBin]);
        hJESvsGenPt[centHistoBin]->Draw("colz");
        hJESvsGenPt[centHistoBin]->GetYaxis()->SetRangeUser(-1., 2.);
        gPad->SetLogz(1);
        t.DrawLatexNDC(0.25, 0.9, Form("JES: %d-%d%% Pb+Pb #sqrt{s_{NN}}=5020 GeV",(iCent-1)*10, (iCent)*10));

        // Jet energy scale
        canvJES->cd(padIndex + centBins);
        setPadStyle();
        set1DStyle(hJESMean[centHistoBin],0);
        hJESMean[centHistoBin]->Draw();
        hJESMean[centHistoBin]->GetYaxis()->SetRangeUser(0.9, 1.3);
        hJESMean[centHistoBin]->GetYaxis()->SetTitle("#mu(JES)");
        set1DStyle(hJESMeanW[centHistoBin],1);
        hJESMeanW[centHistoBin]->Draw("same");

        // Jet energy resolution
        canvJES->cd(padIndex + 2 * centBins);
        setPadStyle();
        set1DStyle(hJESSigma[centHistoBin],0);
        hJESSigma[centHistoBin]->Draw();
        hJESSigma[centHistoBin]->GetYaxis()->SetRangeUser(0., 0.35);
        hJESSigma[centHistoBin]->GetYaxis()->SetTitle("#sigma(JES)");
        set1DStyle(hJESSigmaW[centHistoBin],1);
        hJESSigmaW[centHistoBin]->Draw("same");

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
void drawJetAnalysisFigs(const Char_t *inFileName = "../build/oTestSimpleReadForest.root", 
                         const Char_t *oFile = "oDrawJetAna.root") {
    
    //setStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Read ROOT file                             
    TFile *inFile = TFile::Open(inFileName);
    if ( !inFile->IsOpen() ) {
        std::cout << "Input file not opened. Terminating." << std::endl;
        exit(0);
    }

    // Draw event histograms
    drawEventQuantities(inFile);

    // Draw JES and JER as a function of centrality
    drawJESvsCentrality(inFile);

    // Draw JES and JER as a function of pseudorapidity
    drawJESvsEta(inFile);
}