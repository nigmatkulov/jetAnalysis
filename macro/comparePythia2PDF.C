#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"

#include <iostream>

std::vector<Int_t> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55 };
std::vector<Int_t> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94 };

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
    else if (type == 5) {
        color = 8;        // 
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
void makeRatio(TH1* hRat, TH1* hPb, TH1* hPP) {
    rescaleEta(hPb);
    rescaleEta(hPP);
    hPP->Scale( 1. / hPP->Integral() );
    hRat->Divide( hPb, hPP );
}

//________________
void comparePythia2PDF() {

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    Bool_t usePDFx = kFALSE;
    //Int_t nPtBins2Proc = 5;

    // std::vector<Int_t> ptDijetLow = {};
    // std::vector<Int_t> ptDijetHi = {};

    std::vector<Int_t> ptDijetLow =  { 50, 60, 100, 140};
    std::vector<Int_t> ptDijetHi =   { 60, 70, 110, 150};

    Int_t pdfSetNum = 0;     // 0 - EPPS21, 1 - nCTEQ15HQ
    TString npdfName;
    TString pdfName;
    Int_t nErrorSets = 100;
    if ( pdfSetNum == 0 ) {
        npdfName = "EPPS21";
        pdfName = "CT18";
    }
    else {
        npdfName = "nCTEQ15HQ";
    }
    TString frame = "lab";

    //fillDijetPtBins(ptDijetLow, ptDijetHi);
    Int_t ptBins = ptDijetLow.size();
    Int_t nPtBins2Proc = ptDijetLow.size();
    
    TString pdfInFileName = npdfName;
    pdfInFileName.ToLower();

    TFile *pythiaFile = TFile::Open( Form("20241008/oSystematics_%s.root", frame.Data()) );
    TFile *pdfFile = TFile::Open( Form("epps21/%s_pPb8160.root", pdfInFileName.Data()) );
    TFile *pdfxFile;
    if ( usePDFx ) {
        pdfxFile = TFile::Open("epps21/oEPPS21x.root");
    }

    TH1D *hPP_Pythia[ptDijetLow.size()];
    TH1D *hPP_PDF[ptDijetLow.size()];
    TH1D *hPb_data[ptDijetLow.size()];
    TH1D *hPb_nPDF[ptDijetLow.size()];
    TH1D *hRatio_pp[ptDijetLow.size()];
    TH1D *hRatio_pPb[ptDijetLow.size()];
    TH1D *hPb_fb_data[ptDijetLow.size()];
    TH1D *hPb_fb_nPDF[ptDijetLow.size()];
    TH1D *hPb_fb_ratio[ptDijetLow.size()];

    TH1D *hPP_PDF_x;
    TH1D *hPP_nPDF_x;
    if ( usePDFx ) { 
        hPP_PDF_x = dynamic_cast<TH1D*> ( pdfxFile->Get("hPP_50_60"));
        hPP_nPDF_x = dynamic_cast<TH1D*> ( pdfxFile->Get("hPPb_50_60"));
    }

    TH1D *hData2Pythia[ptDijetLow.size()];
    TH1D *hData2PDF[ptDijetLow.size()];
    TH1D *hNPDF2Pythia[ptDijetLow.size()];
    TH1D *hNPDF2PDF[ptDijetLow.size()];

    TLatex t;
    TLegend *leg;


    TCanvas *canv = new TCanvas("canv", "canv", 1800, 600);
    canv->Divide( nPtBins2Proc, 2 );

    TCanvas *canv2 = new TCanvas("canv2", "canv2", 1800, 600);
    canv2->Divide( nPtBins2Proc, 2 );

    // TCanvas *canv4 = new TCanvas("canv4", "canv4", 1800, 600);
    // canv4->Divide( nPtBins2Proc, 2 );

    TCanvas *canv3 = new TCanvas("canv3", "canv3", 1800, 600);
    canv3->Divide( nPtBins2Proc, 2 );

    // Loop over pt bins and read histograms
    for (Int_t i{0}; i<nPtBins2Proc; i++) {

        cout << "Processing pt bin " << i << " ptLow = " << ptDijetLow.at(i) << " ptHi = " << ptDijetHi.at(i) << endl;

        // proton-proton
        hPP_Pythia[i] = dynamic_cast<TH1D*> ( pythiaFile->Get( Form("Embed_pPb8160_etaDijet_gen_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        hPP_Pythia[i]->SetName( Form("pythia_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) );
        set1DStyle( hPP_Pythia[i], 0 );
        // rescaleEta( hPP_Pythia[i] );

        hPP_PDF[i] = dynamic_cast<TH1D*> ( pdfFile->Get( Form("ppLab_pt_%d_%d_0", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        hPP_PDF[i]->SetName( Form("pdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) );
        set1DStyle( hPP_PDF[i], 1 );
        // rescaleEta( hPP_PDF[i] ); 

        hRatio_pp[i] = dynamic_cast<TH1D*> ( hPP_PDF[i]->Clone( Form("pdf2pythia_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        hRatio_pp[i]->SetName( Form("pdf2pythia_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hRatio_pp[i], 2);
        hRatio_pp[i]->GetYaxis()->SetTitle( "CT18 / PYTHIA8" );
        makeRatio( hRatio_pp[i], hPP_PDF[i], hPP_Pythia[i] );

        // proton-lead
        TString trgName;
        if ( ptDijetHi.at(i) < 80 ) {
            trgName = "MB";
        }
        else if ( ptDijetHi.at(i) < 100 ) {
            trgName = "Jet60";
        }
        else if ( ptDijetHi.at(i) < 120 ) {
            trgName = "Jet80";
        }
        else {
            trgName = "Jet100";
        }

        hPb_data[i] = dynamic_cast<TH1D*> ( pythiaFile->Get( Form("%s_pPb8160_etaDijet_pt_%d_%d", trgName.Data(), ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        hPb_data[i]->SetName( Form("pPb_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hPb_data[i], 0);

        hPb_nPDF[i] = dynamic_cast<TH1D*> ( pdfFile->Get( Form("pPbLab_pt_%d_%d_0", ptDijetLow.at(i), ptDijetHi.at(i) ) ) );
        hPb_nPDF[i]->SetName( Form("npdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hPb_nPDF[i], 1);

        // hPb_fb_data[i] = dynamic_cast<TH1D*> ( pythiaFile->Get( Form("%s_pPb8160_etaDijet_fb_pt_%d_%d", trgName.Data(), ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        // hPb_fb_data[i]->SetName( Form("pPb_fb_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hPb_fb_data[i], 0);

        // hPb_fb_nPDF[i] = dynamic_cast<TH1D*> ( pdfFile->Get( Form("pPbCM_fb_pt_%d_%d_0", ptDijetLow.at(i), ptDijetHi.at(i) ) ) );
        // hPb_fb_nPDF[i]->SetName( Form("npdf_fb_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hPb_fb_nPDF[i], 1);

        // hPb_fb_ratio[i] = dynamic_cast<TH1D*> ( hPb_fb_data[i]->Clone( Form("fb2npdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        // hPb_fb_ratio[i]->SetName( Form("fb2npdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hPb_fb_ratio[i], 2);

        hRatio_pPb[i] = dynamic_cast<TH1D*> ( hPb_data[i]->Clone( Form("data2npdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        hRatio_pPb[i]->SetName( Form("data2npdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); set1DStyle( hRatio_pPb[i], 2);
        hRatio_pPb[i]->GetYaxis()->SetTitle( "pPb (CMS) / nPDF" );
        makeRatio( hRatio_pPb[i], hPb_data[i], hPb_nPDF[i] );

        // Plot various ratios pPb to pp
        hData2Pythia[i] = dynamic_cast<TH1D*> ( hPb_data[i]->Clone( Form("data2pythia_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        makeRatio( hData2Pythia[i], hPb_data[i], hPP_Pythia[i] );
        hData2Pythia[i]->GetYaxis()->SetTitle("pPb / pp");
        set1DStyle( hData2Pythia[i], 2 );

        hData2PDF[i] = dynamic_cast<TH1D*> ( hPb_data[i]->Clone( Form("data2npdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        makeRatio( hData2PDF[i], hPb_data[i], hPP_PDF[i] );
        hData2PDF[i]->GetYaxis()->SetTitle("pPb / pp");
        set1DStyle( hData2PDF[i], 5 );

        hNPDF2Pythia[i] = dynamic_cast<TH1D*> ( hPb_nPDF[i]->Clone( Form("npdf2pythia_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) );
        makeRatio( hNPDF2Pythia[i], hPb_nPDF[i], hPP_Pythia[i] );
        hNPDF2Pythia[i]->GetYaxis()->SetTitle("pPb / pp");
        set1DStyle( hNPDF2Pythia[i], 0 );

        hNPDF2PDF[i] = dynamic_cast<TH1D*> ( pdfFile->Get( Form("pPb2ppLabRatio_pt_%d_%d_0", ptDijetLow.at(i), ptDijetHi.at(i) ) ) );
        hNPDF2PDF[i]->SetName( Form("npdf2pdf_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ); 
        set1DStyle( hNPDF2PDF[i], 1);

        //
        // Plot comparisons
        // 

        cout << "There 1" << endl;

        // proton-proton
        canv->cd( i + 1 );
        setPadStyle();
        hPP_PDF[i]->Draw();
        hPP_Pythia[i]->Draw("same");
        hPP_PDF[i]->GetXaxis()->SetRangeUser( -3., 3. );
        hPP_PDF[i]->GetYaxis()->SetRangeUser( 0.001, 0.12 );
        t.DrawLatexNDC(0.35, 0.92, Form("%d < p_{T}^{ave} GeV < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
        leg = new TLegend(0.4, 0.75, 0.85, 0.85);
        leg->AddEntry(hPP_Pythia[i], "PYTHIA8", "p");
        if ( pdfSetNum == 0 ) {
            leg->AddEntry(hPP_PDF[i], "CT18", "p");
        }
        else {
            leg->AddEntry(hPP_PDF[i], "nCTEQ15HQ pp", "p");
        }
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->Draw();

        canv->cd( nPtBins2Proc + i + 1 );
        setPadStyle();
        hRatio_pp[i]->Draw();
        hRatio_pp[i]->GetXaxis()->SetRangeUser(-3., 3.);
        hRatio_pp[i]->GetYaxis()->SetRangeUser( 0.8, 1.15 );

        // proton-lead
        canv2->cd( i + 1 );
        setPadStyle();
        hPb_data[i]->Draw();
        hPb_nPDF[i]->Draw("same");
        hPb_data[i]->GetXaxis()->SetRangeUser( -3., 3. );
        hPb_data[i]->GetYaxis()->SetRangeUser( 0.001, 0.12 );
        t.DrawLatexNDC(0.35, 0.92, Form("%d < p_{T}^{ave} GeV < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
        leg = new TLegend(0.4, 0.75, 0.85, 0.85);
        leg->AddEntry(hPb_data[i], "pPb (CMS)", "p");
        if ( pdfSetNum == 0 ) {
            leg->AddEntry(hPb_nPDF[i], "EPPS21 x CT18", "p");
        }
        else {
            leg->AddEntry(hPb_nPDF[i], "nCTEQ15HQ", "p");
        }
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->Draw();

        canv2->cd( nPtBins2Proc + i + 1 );
        setPadStyle();
        hRatio_pPb[i]->Draw();
        hRatio_pPb[i]->GetXaxis()->SetRangeUser(-3., 3.);
        hRatio_pPb[i]->GetYaxis()->SetRangeUser( 0.8, 1.15 );

        // Compare different pPb to pp ratios
        canv3->cd( i + 1 );
        setPadStyle();
        hData2Pythia[i]->Draw();
        hData2PDF[i]->Draw("same");
        hData2Pythia[i]->GetXaxis()->SetRangeUser(-3., 3.);
        hData2Pythia[i]->GetYaxis()->SetRangeUser(0.75, 1.25);
        t.DrawLatexNDC(0.35, 0.92, Form("%d < p_{T}^{ave} GeV < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
        leg = new TLegend(0.25, 0.75, 0.65, 0.85);
        leg->AddEntry(hData2Pythia[i], "pPb (CMS) / PYTHIA8", "p");
        if ( pdfSetNum == 0 ) {
            leg->AddEntry(hData2PDF[i], "pPb (CMS) / CT18", "p");
        }
        else {
            leg->AddEntry(hData2PDF[i], "pPb (CMS) / nCTEQ15HQ pp", "p");
        }
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->Draw();

        canv3->cd( nPtBins2Proc + i + 1 );
        hNPDF2Pythia[i]->Draw();
        hNPDF2PDF[i]->Draw("same");
        hNPDF2Pythia[i]->GetXaxis()->SetRangeUser(-3., 3.);
        hNPDF2Pythia[i]->GetYaxis()->SetRangeUser(0.75, 1.25);
        leg = new TLegend(0.25, 0.75, 0.65, 0.85);
        if ( pdfSetNum == 0 ) {
            leg->AddEntry(hNPDF2Pythia[i], "EPPS21 x CT18 / PYTHIA", "p");
            leg->AddEntry(hNPDF2PDF[i], "EPPS21 x CT18 / CT18", "p");
        }
        else {
            leg->AddEntry(hNPDF2Pythia[i], "nCTEQ15HQ / PYTHIA", "p");
            leg->AddEntry(hNPDF2PDF[i], "nCTEQ15HQ / PDF", "p");
        }
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->Draw();

        // // Forward-backward ratios
        // canv4->cd( i + 1 );
        // setPadStyle();
        // hPb_fb_data[i]->Draw();
        // hPb_fb_nPDF[i]->Draw("same");
        // hPb_fb_data[i]->GetXaxis()->SetRangeUser( 0., 3. );
        // hPb_fb_data[i]->GetYaxis()->SetRangeUser( 0.65, 1.35 );
        // t.DrawLatexNDC(0.35, 0.92, Form("%d < p_{T}^{ave} GeV < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
        // leg = new TLegend(0.35, 0.7, 0.85, 0.85);
        // leg->AddEntry(hPb_fb_data[i], "pPb (CMS)", "p");
        // leg->AddEntry(hPb_fb_nPDF[i], "EPPS21 x CT18", "p");
        // leg->SetBorderSize(0);
        // leg->SetTextSize(0.05);
        // leg->Draw();

        // canv4->cd( nPtBins2Proc + i + 1 );
        // setPadStyle();
        // hPb_fb_ratio[i]->Draw();
        // hPb_fb_ratio[i]->GetXaxis()->SetRangeUser(0., 3.);
        // hPb_fb_ratio[i]->GetYaxis()->SetRangeUser( 0.85, 1.15 );

    } // 

    canv->SaveAs( Form("epps21/%s_pdf2pythiaComp_%s.pdf", npdfName.Data(), frame.Data()) );
    canv2->SaveAs( Form("epps21/%s_data2npdfComp_%s.pdf", npdfName.Data(), frame.Data()) );
    canv3->SaveAs( Form("epps21/%s_pPb2ppComp_%s.pdf", npdfName.Data(), frame.Data()) );

}