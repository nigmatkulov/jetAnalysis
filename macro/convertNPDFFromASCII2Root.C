#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"

#include <iostream>

std::vector<Int_t> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55 };
std::vector<Int_t> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94 };

const Int_t dijetEtaBins{30};
Double_t dijetEtaVals[dijetEtaBins+1] { -5.0, -4.0, -3.0, -2.4, -2.2, 
                                        -2.0, -1.8, -1.6, -1.4, -1.2, 
                                        -1.0, -0.8, -0.6, -0.4, -0.2,  
                                        0.0,  0.2,  0.4,  0.6,  0.8,  
                                        1.0,  1.2,  1.4,  1.6,  1.8,  
                                        2.0,  2.2,  2.4,  3.0,  4.0,  
                                        5.0 };

//________________
void fillDijetPtBins(std::vector<Int_t> &ptDijetLow, std::vector<Int_t> &ptDijetHi) {
    Int_t ptStep {5};
    Int_t ptLow {30};
    for (UInt_t i{0}; i<ptDijetBinLow.size(); i++) {
        ptDijetLow.push_back( ptLow + (ptDijetBinLow.at(i)-1) * ptStep );
        ptDijetHi.push_back( ptLow + ptDijetBinHi.at(i) * ptStep );
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
    } // for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++)
    h->Scale( 1. / h->Integral() );
}

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

//_________________
TH1D* makeRatio(TH1* hPb, TH1* hPP, Int_t ptLow, Int_t ptHi) {

    TH1D *hRat = dynamic_cast<TH1D*> ( hPb->Clone( Form("hRatio_%d_%d", ptLow, ptHi) ) );
    hRat->GetYaxis()->SetTitle("pPb/pp");
    hRat->Divide( hPb, hPP );
    set1DStyle( hRat, 2 );
    hRat->GetXaxis()->SetRangeUser(-3., 3.);
    hRat->GetYaxis()->SetRangeUser(0.8, 1.15);

    return hRat;
}

//_________________
TH1D* readFromASCII(Int_t ptLow, Int_t ptHi, Bool_t isPPb = kFALSE, Bool_t isX = kFALSE) {

    Int_t type = (isPPb) ? 0 : 1;
    TString system = (isPPb) ? "pPb" : "pp";
    TString xSystem = (isX) ? "x" : "";

    TString inFileName = "epps21/oEPPS21";
    inFileName += xSystem; inFileName += "_"; inFileName += system; inFileName += "_pt_"; 
    inFileName += ptLow; inFileName += "_"; inFileName += ptHi; inFileName += ".txt";

    TGraphErrors* grErr = new TGraphErrors(inFileName.Data(), "%lg %lg %lg");

    TString hName = ( isPPb ) ? Form("hPPb_%d_%d", ptLow, ptHi) : Form("hPP_%d_%d", ptLow, ptHi);

    TH1D* hTmp = new TH1D( hName.Data(), 
                           Form("EPPS21 %s %d-%d GeV;#eta_{dijet};1/N_{dijet} dN/d#eta_{dijet}", system.Data(), ptLow, ptHi), 
                           dijetEtaBins, dijetEtaVals);
    hTmp->Sumw2();
    ( isPPb ) ? set1DStyle(hTmp, 0) : set1DStyle(hTmp, 1);

    // Fill the histograms
    for (Int_t i{0}; i<dijetEtaBins; i++) {
        hTmp->SetBinContent( i+1, grErr->GetY()[i] );
        hTmp->SetBinError( i+1, grErr->GetEY()[i] );
    } // for (Int_t i{0}; i<dijetEtaBins; i++)

    rescaleEta( hTmp );

    hTmp->GetXaxis()->SetRangeUser(-3., 3.);
    hTmp->GetYaxis()->SetRangeUser(0.001, 0.12);

    return hTmp;
}

//_________________
void convertNPDFFromASCII2Root() {

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::vector<Int_t> ptDijetLow = {};
    std::vector<Int_t> ptDijetHi = {};
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    Bool_t isX = kFALSE;
    Int_t nBins2Read = 5;

    TString outFileName = "oEPPS21";
    (isX) ? outFileName += "x.root" : outFileName += ".root";

    TH1D *hPP[ ptDijetBinLow.size() ];
    TH1D *hPPb[ ptDijetBinLow.size() ];
    TH1D *hRat[ ptDijetBinLow.size() ];

    // Loop over bins
    for (Int_t i{0}; i<nBins2Read; i++) {
        hPP[i]  = readFromASCII( ptDijetLow.at(i), ptDijetHi.at(i), kFALSE, isX);
        hPPb[i] = readFromASCII( ptDijetLow.at(i), ptDijetHi.at(i), kTRUE, isX);
        hRat[i] = makeRatio( hPPb[i], hPP[i], ptDijetLow.at(i), ptDijetHi.at(i) );
    } // for (Int_t i{0}; i<nBins2Read; i++)

    TLatex t;
    t.SetTextSize(0.04);
    t.SetTextFont(42);

    TLegend *leg;

    // Plot the histograms on the canvas
    TCanvas *can = new TCanvas("can", "can", 1200, 800);
    can->Divide(nBins2Read, 2);

    // Loop over pT bins
    for (Int_t i{0}; i<nBins2Read; i++) {
        can->cd( i + 1 );
        setPadStyle();
        hPPb[i]->Draw();
        hPP[i]->Draw("same");
        t.DrawLatexNDC(0.35, 0.92, "50 < p_{T}^{ave} < 60 GeV");
        leg = new TLegend(0.75, 0.75, 0.85, 0.85);
        leg->AddEntry(hPPb[i], "pPb", "p");
        leg->AddEntry(hPP[i], "pp", "p");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->Draw();

        can->cd( nBins2Read + i + 1 );
        setPadStyle();
        hRat[i]->Draw();
    } // for (Int_t i{0}; i<nBins2Read; i++)
    

    can->SaveAs("epps21/oEPPS21.pdf");
    can->SaveAs("epps21/oEPPS21.png");

    // Create output file and save histograms
    TFile *outFile = TFile::Open(outFileName.Data(), "RECREATE");

    // Loop over bins
    for (Int_t i{0}; i<nBins2Read; i++) {
        hPP[i]->Write();
        hPPb[i]->Write();
        hRat[i]->Write();
    }
    outFile->Close();
}