#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>
#include <vector>

std::vector<int> ptLow { 60, 100, 140 };
std::vector<int> ptHi  { 70, 110, 150 };

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
}

//________________
void setPadFBStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.25);
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
        color = 8;        // green
        markerStyle = 39; // radiation
    }
    else if ( type == 6 ) {
        color = 9;        // purple
        markerStyle = 20; // circle
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
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void make1DRatio(std::vector<TH1D*> &hNum, std::vector<TH1D*> &hDen, std::vector<TH1D*> &hRat,
                 const Char_t *name, Int_t style = 2) {
    for (Int_t i{0}; i<hNum.size(); i++) {
        hRat.push_back( (TH1D*)hNum[i]->Clone( Form("%s_pt_%d_%d", name, ptLow[i], ptHi[i]) ) );
        hRat.back()->Divide( hDen[i] );
        set1DStyle( hRat.back(), style );
    } // for (Int_t i{0}; i<hNum.size(); i++)
}

//________________
void readHistograms(TFile *f, 
                    std::vector<TH1D*> &ppLab, std::vector<TH1D*> &ppCM, 
                    std::vector<TH1D*> &ppCMForward, std::vector<TH1D*> &ppCMBackward, 
                    std::vector<TH1D*> &ppCMFBRatio, std::vector<TH1D*> &ppCMBFRatio, 
                    std::vector<TH1D*> &pPbLab, std::vector<TH1D*> &pPbCM, 
                    std::vector<TH1D*> &pPbCMForward, std::vector<TH1D*> &pPbCMBackward, 
                    std::vector<TH1D*> &pPbCMFBRatio, std::vector<TH1D*> &pPbCMBFRatio, 
                    std::vector<TH1D*> &RpPbLab, std::vector<TH1D*> &RpPbCM) {

    // Style codes
    // 0: red   - pPb
    // 1: blue  - pp
    // 2: black - RpPb
    Int_t pPbStyle{0}, ppStyle{1}, RpPbStyle{2};

    // Loop over pT bins
    for (Int_t i{0}; i<ptLow.size(); i++) {
        ppLab.push_back( (TH1D*)f->Get( Form("ppLab_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( ppLab.back(), ppStyle, kTRUE );
        ppCM.push_back( (TH1D*)f->Get( Form("ppCM_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( ppCM.back(), ppStyle, kTRUE );
        ppCMForward.push_back( (TH1D*)f->Get( Form("ppCMForward_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( ppCMForward.back(), ppStyle, kTRUE );
        ppCMBackward.push_back( (TH1D*)f->Get( Form("ppCMBackward_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( ppCMBackward.back(), ppStyle, kTRUE );

        pPbLab.push_back( (TH1D*)f->Get( Form("pPbLab_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( pPbLab.back(), pPbStyle, kTRUE );
        pPbCM.push_back( (TH1D*)f->Get( Form("pPbCM_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( pPbCM.back(), pPbStyle, kTRUE );
        pPbCMForward.push_back( (TH1D*)f->Get( Form("pPbCMForward_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( pPbCMForward.back(), pPbStyle, kTRUE );
        pPbCMBackward.push_back( (TH1D*)f->Get( Form("pPbCMBackward_pt_%d_%d_0", ptLow[i], ptHi[i]) ) ); set1DStyle( pPbCMBackward.back(), pPbStyle, kTRUE );
    } // for (Int_t i{0}; i<ptLow.size(); i++)

    // Make ratios
    make1DRatio( pPbLab, ppLab, RpPbLab, "RpPbLab", RpPbStyle );
    make1DRatio( pPbCM, ppCM, RpPbCM, "RpPbCM", RpPbStyle );
    
    make1DRatio( ppCMForward, ppCMBackward, ppCMFBRatio, "ppCMFBRatio", ppStyle );
    make1DRatio( ppCMBackward, ppCMForward, ppCMBFRatio, "ppCMBFRatio", ppStyle );

    make1DRatio( pPbCMForward, pPbCMBackward, pPbCMFBRatio, "pPbCMFBRatio", pPbStyle );
    make1DRatio( pPbCMBackward, pPbCMForward, pPbCMBFRatio, "pPbCMBFRatio", pPbStyle );
}

//________________
void plotHeader(TString npdfName) {
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    t.DrawLatexNDC(0.2, 0.92, Form("%s", npdfName.Data()) );
    t.DrawLatexNDC(0.65, 0.92, "#sqrt{s} = 8.16 TeV" );
}

//________________
void plotRpPb(std::vector<TH1D*> &RpPb, TString npdfName, Bool_t isCM = kFALSE) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.7, 1.2};
    TLine *line = new TLine( xRange[0], 1.0, xRange[1], 1.0 );
    line->SetLineStyle(2); // dashed
    line->SetLineColor(6); // magenta
    line->SetLineWidth(2);

    TString frame = ( isCM ) ? "CM" : "Lab";

    TCanvas *c = new TCanvas( "c", "c", 800, 600 );
    for (Int_t i{0}; i<RpPb.size(); i++) {
        c->cd();
        setPadStyle();
        RpPb[i]->Draw();
        RpPb[i]->GetXaxis()->SetRangeUser( xRange[0], xRange[1] );
        RpPb[i]->GetYaxis()->SetRangeUser( yRange[0], yRange[1] );
        RpPb[i]->GetYaxis()->SetTitle("R_{pPb}");
        line->Draw("same");
        t.DrawLatexNDC(0.35, 0.82, Form("%d < p_{T}^{dijet} (GeV/c) < %d", ptLow[i], ptHi[i]) );
        plotHeader( npdfName );

        //c->SaveAs( Form("npdf/figs/%s_RpPb%s_pt_%d_%d.pdf", npdfName.Data(), frame.Data(), ptLow[i], ptHi[i]) );
        c->SaveAs( Form("npdf/figs/%s_RpPb%s_pt_%d_%d.png", npdfName.Data(), frame.Data(), ptLow[i], ptHi[i]) );
        //c->SaveAs( Form("npdf/figs/%s_RpPb%s_pt_%d_%d.C", npdfName.Data(), frame.Data(), ptLow[i], ptHi[i]) );
    } // for (Int_t i{0}; i<RpPbLab.size(); i++)
}

//________________
void plotComparison(std::vector<TH1D*> &pPb, std::vector<TH1D*> &pp, 
                    Int_t type = 0, TString npdfName = "epps21", 
                    Bool_t isCM = kFALSE) {

    // type: 0 - full, 1 - forward, 2 - backward, 3 - fb, 4 - bf
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.0001, 0.12};
    Double_t xLegend[2] = {0.7, 0.8};
    Double_t yLegend[2] = {0.65, 0.75};
    TString direction = ( isCM ) ? "dN/d#eta^{dijet}_{CM}" : "dN/d#eta^{dijet}";
    TString dir = (isCM) ? "CM" : "Lab";
    if ( type == 1 || type == 2 ) {
        xRange[0] = 0.0; xRange[1] = 3.0;
        yRange[0] = 0.0001; yRange[1] = 0.18;
        if ( type == 1 ) {
            direction = "dN/d#eta^{dijet}_{CM} (#eta>0)";
            dir = "_forward";
        }
        else {
            direction = "dN/d#eta^{dijet}_{CM} (#eta<0)";
            dir = "_backward";
        }
    }
    else if ( type == 3 || type == 4 ) {
        xRange[0] = 0.0; xRange[1] = 3.0;
        yRange[0] = 0.6; yRange[1] = 1.4;
        if ( type == 3 ) {
            direction = "Forward / Backward";
            dir = "_fb";
        }
        else {
            direction = "Backward / Forward";
            dir = "_bf";
        }
        xLegend[0] = 0.2; xLegend[1] = 0.35;
        yLegend[0] = 0.2; yLegend[1] = 0.3;
    }

    TLine *line = new TLine( xRange[0], 1.0, xRange[1], 1.0 );
    line->SetLineStyle(2); // dashed
    line->SetLineColor(6); // magenta
    line->SetLineWidth(2);

    TCanvas *c = new TCanvas( "c", "c", 800, 600 );

    for (Int_t i{0}; i<pPb.size(); i++) {

        c->cd();
        setPadStyle();
        pPb[i]->Draw();
        pp[i]->Draw("same");
        pPb[i]->GetXaxis()->SetRangeUser( xRange[0], xRange[1] );
        pPb[i]->GetYaxis()->SetRangeUser( yRange[0], yRange[1] );
        pPb[i]->GetYaxis()->SetTitle( Form( "%s", direction.Data() ) );
        if ( type == 3 || type == 4 ) {
            line->Draw("same");
        }
        t.DrawLatexNDC(0.35, 0.82, Form("%d < p_{T}^{dijet} (GeV/c) < %d", ptLow[i], ptHi[i]) );
        plotHeader( npdfName );

        TLegend *leg = new TLegend(xLegend[0], yLegend[0], xLegend[1], yLegend[1]);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        leg->AddEntry( pPb[i], "pPb", "p" );
        leg->AddEntry( pp[i], "pp", "p" );
        leg->Draw("same");

        //c->SaveAs( Form("npdf/figs/%s_%sRatio_pt_%d_%d.pdf", npdfName.Data(), dir.Data(), ptLow[i], ptHi[i]) );
        c->SaveAs( Form("npdf/figs/%s_eta%s_pt_%d_%d.png", npdfName.Data(), dir.Data(), ptLow[i], ptHi[i]) );
        // c->SaveAs( Form("npdf/figs/%s_%sRatio_pt_%d_%d.C", npdfName.Data(), dir.Data(), ptLow[i], ptHi[i]) );

    } // for (Int_t i{0}; i<pPb.size(); i++)
}


//________________
void plotDirRatios(std::vector<TH1D*> &pPbDirRatio, std::vector<TH1D*> &ppDirRatio, 
                   Bool_t isFB = kTRUE, TString npdfName = "epps21") {
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    Double_t xRange[2] = {0., 3.};
    Double_t yRange[2] = {0.8, 1.2};
    TLine *line = new TLine( xRange[0], 1.0, xRange[1], 1.0 );
    line->SetLineStyle(2); // dashed
    line->SetLineColor(6); // magenta
    line->SetLineWidth(2);

    TString direction = ( isFB ) ? " Forward / Backward" : "Backward / Forward";
    TString dir = ( isFB ) ? "fb" : "bf";

    TCanvas *c = new TCanvas( "c", "c", 800, 600 );

    for (Int_t i{0}; i<pPbDirRatio.size(); i++) {
        c->cd();
        setPadStyle();
        pPbDirRatio[i]->Draw();
        ppDirRatio[i]->Draw("same");
        pPbDirRatio[i]->GetXaxis()->SetRangeUser( xRange[0], xRange[1] );
        pPbDirRatio[i]->GetYaxis()->SetRangeUser( yRange[0], yRange[1] );
        pPbDirRatio[i]->GetYaxis()->SetTitle( Form( "%s", direction.Data() ) );
        line->Draw("same");
        t.DrawLatexNDC(0.35, 0.82, Form("%d < p_{T}^{dijet} (GeV/c) < %d", ptLow[i], ptHi[i]) );
        plotHeader( npdfName );

        TLegend *leg = new TLegend(0.2, 0.2, 0.35, 0.3);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        leg->AddEntry( pPbDirRatio[i], "pPb", "p" );
        leg->AddEntry( ppDirRatio[i], "pp", "p" );
        leg->Draw("same");

        //c->SaveAs( Form("npdf/figs/%s_%sRatio_pt_%d_%d.pdf", npdfName.Data(), dir.Data(), ptLow[i], ptHi[i]) );
        c->SaveAs( Form("npdf/figs/%s_%sRatio_pt_%d_%d.png", npdfName.Data(), dir.Data(), ptLow[i], ptHi[i]) );
        // c->SaveAs( Form("npdf/figs/%s_%sRatio_pt_%d_%d.C", npdfName.Data(), dir.Data(), ptLow[i], ptHi[i]) );
    } // for (Int_t i{0}; i<pPbDirRatio.size(); i++)
}

//________________
void plotDistributions(std::vector<TH1D*> &ppLab, std::vector<TH1D*> &ppCM, 
                       std::vector<TH1D*> &ppCMForward, std::vector<TH1D*> &ppCMBackward, 
                       std::vector<TH1D*> &ppCMFBRatio, std::vector<TH1D*> &ppCMBFRatio, 
                       std::vector<TH1D*> &pPbLab, std::vector<TH1D*> &pPbCM, 
                       std::vector<TH1D*> &pPbCMForward, std::vector<TH1D*> &pPbCMBackward, 
                       std::vector<TH1D*> &pPbCMFBRatio, std::vector<TH1D*> &pPbCMBFRatio, 
                       std::vector<TH1D*> &RpPbLab, std::vector<TH1D*> &RpPbCM,
                       TString npdfName) {

    plotRpPb( RpPbLab, npdfName, kFALSE );
    plotRpPb( RpPbCM, npdfName, kTRUE );

    // plotDirRatios( pPbCMFBRatio, ppCMFBRatio, kTRUE, npdfName );
    // plotDirRatios( pPbCMBFRatio, ppCMBFRatio, kFALSE, npdfName );

    plotComparison( pPbLab, ppLab, 0, npdfName, kFALSE );
    plotComparison( pPbCM, ppCM, 0, npdfName, kTRUE );
    plotComparison( pPbCMForward, ppCMForward, 1, npdfName, kTRUE );
    plotComparison( pPbCMBackward, ppCMBackward, 2, npdfName, kTRUE );
    plotComparison( pPbCMFBRatio, ppCMFBRatio, 3, npdfName, kTRUE );
    plotComparison( pPbCMBFRatio, ppCMBFRatio, 4, npdfName, kTRUE );

}

//________________
void plotNPDFDistributions() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TString npdfName = "epps21";
    //npdfName = "ncteq15hq";

    // Open the ROOT file to read
    TFile *inFile = TFile::Open( Form( "npdf/%s_pPb8160.root", npdfName.Data() ), "READ");

    // pp histograms
    std::vector<TH1D*> ppLab;
    std::vector<TH1D*> ppCM;
    std::vector<TH1D*> ppCMForward;
    std::vector<TH1D*> ppCMBackward;
    std::vector<TH1D*> ppCMFBRatio;
    std::vector<TH1D*> ppCMBFRatio;

    // pPb histograms
    std::vector<TH1D*> pPbLab;
    std::vector<TH1D*> pPbCM;
    std::vector<TH1D*> pPbCMForward;
    std::vector<TH1D*> pPbCMBackward;
    std::vector<TH1D*> pPbCMFBRatio;
    std::vector<TH1D*> pPbCMBFRatio;

    // RpPb histograms
    std::vector<TH1D*> RpPbLab;
    std::vector<TH1D*> RpPbCM;

    // Read histograms
    readHistograms( inFile, ppLab, ppCM, ppCMForward, ppCMBackward, ppCMFBRatio, ppCMBFRatio, 
                    pPbLab, pPbCM, pPbCMForward, pPbCMBackward, pPbCMFBRatio, pPbCMBFRatio, 
                    RpPbLab, RpPbCM );


    // Plot distributions
    plotDistributions( ppLab, ppCM, ppCMForward, ppCMBackward, ppCMFBRatio, ppCMBFRatio, 
                       pPbLab, pPbCM, pPbCMForward, pPbCMBackward, pPbCMFBRatio, pPbCMBFRatio, 
                       RpPbLab, RpPbCM, npdfName );

}