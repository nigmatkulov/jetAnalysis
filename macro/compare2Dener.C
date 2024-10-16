#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THnSparse.h"

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
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
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


//________________
void make1DRatio(TH1D *hRat, TH1D *hDen, const Char_t *ratioName = "Ratio to default", Int_t style = 0, Bool_t isBinomial = kFALSE) {

    if ( !hDen ) {
        std::cout << "Denominator does not exist" << std::endl;
    }

    if ( isBinomial ) {
        hRat->Divide( hRat, hDen, 1., 1., "b" );
    }
    else {
        hRat->Divide( hRat, hDen, 1., 1. );
    }
    hRat->GetYaxis()->SetTitle( ratioName );
    hRat->GetYaxis()->SetRangeUser(0.8, 1.2);
    //hRat->GetXaxis()->SetRangeUser(-3., 3.);
    set1DStyle(hRat, style);
}

//________________
void plotComparison(TCanvas *c, TH1D *h1, TH1D *h2, TH1D *hRatio, 
                    const char *h1Type, const char *h2Type, 
                    int ptLow = 95, int ptHi = 115) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    set1DStyle(h1, 0);
    set1DStyle(h2, 1);
    set1DStyle(hRatio, 2);

    c->cd(1);
    setPadStyle();
    h1->Draw();
    h2->Draw("same");
    h1->GetYaxis()->SetRangeUser(0., 0.12);
    h1->GetXaxis()->SetRangeUser(-3., 3.);
    t.DrawLatexNDC(0.35, 0.92, Form("%d < p_{T}^{ave} GeV < %d", ptLow, ptHi ) );

    TLegend *leg = new TLegend(0.17, 0.75, 0.45, 0.85);
    leg->SetTextSize(0.04);
    leg->AddEntry(h1, Form("%s", h1Type), "p");
    leg->AddEntry(h2, Form("%s", h2Type), "p");
    leg->SetBorderSize(0);
    leg->Draw();

    c->cd(2);
    setPadStyle();
    hRatio->Draw();
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetXaxis()->SetRangeUser(-3., 3.);
}

//________________
void compare2Dener() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    std::vector<Int_t> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55};
    std::vector<Int_t> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94};

    const Int_t dijetEtaBins{30};
    Double_t dijetEtaVals[dijetEtaBins+1] { -5.0, -4.0, -3.0, -2.4, -2.2, 
                                            -2.0, -1.8, -1.6, -1.4, -1.2, 
                                            -1.0, -0.8, -0.6, -0.4, -0.2,  
                                            0.0,  0.2,  0.4,  0.6,  0.8,  
                                            1.0,  1.2,  1.4,  1.6,  1.8,  
                                            2.0,  2.2,  2.4,  3.0,  4.0,  
                                            5.0 };

    TFile *inDenerFile = TFile::Open("fromDener/output_Dener.root");
    TFile *inMyFile = TFile::Open("fromDener/oMB_pPb8160_Pbgoing_comp2Dener.root");

    TH3D *hMyPtEtaDphi = dynamic_cast<TH3D*>(inMyFile->Get("hRecoDijetPtEtaDphiWeighted"));
    TH1D *hMyEta = dynamic_cast<TH1D*>( hMyPtEtaDphi->ProjectionY("hMyEtaProj", ptDijetBinLow.at(0), ptDijetBinHi.at(0)) );
    cout << "hMyEta->Integral() = " << hMyEta->Integral() << endl;
    //rescaleEta(hMyEta);
    set1DStyle(hMyEta, 0);

    THnSparse *hDener = dynamic_cast<THnSparse*>( inDenerFile->Get("dijet_histograms/hist_etaDijet_reco") );
    hDener->GetAxis(4)->SetRange(21, 31);
    hDener->GetAxis(8)->SetRange(6, 6);
    TH1D *hEtaTmp = dynamic_cast<TH1D*>(hDener->Projection(0));
    TH1D *hDenerEta = new TH1D("hDenerEta", "hDenerEta", dijetEtaBins, dijetEtaVals);
    hDenerEta->Sumw2();
    for (Int_t i=1; i<=hEtaTmp->GetNbinsX(); i++) {
        Double_t binCenter = hEtaTmp->GetBinCenter(i);
        Int_t bin = hDenerEta->FindBin(binCenter);
        Int_t binContent = hEtaTmp->GetBinContent(i);
        binContent += hDenerEta->GetBinContent(bin);
        hDenerEta->SetBinContent(bin, binContent);
        hDenerEta->SetBinError(bin, TMath::Sqrt(binContent));
    }
    cout << "hDenerEta->Integral() = " << hDenerEta->Integral() << endl;
    //rescaleEta(hDenerEta);
    set1DStyle(hDenerEta, 1);

    TH1D *hMy2DenerRatio = dynamic_cast<TH1D*>( hMyEta->Clone("hMy2DenerRatio") );
    //hMy2DenerRatio->Reset();
    make1DRatio(hMy2DenerRatio, hDenerEta, "My / Dener", 2);

    // Plot comparison
    TCanvas *c = new TCanvas("c", "c", 800, 1200);
    c->Divide(1, 2);
    plotComparison(c, hMyEta, hDenerEta, hMy2DenerRatio, "My", "Dener", 50, 60);

    TFile *oFile = TFile::Open("fromDener/compare2Dener.root", "recreate");
    hMyEta->Write();
    hDenerEta->Write();
    oFile->Close();
}