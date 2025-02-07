#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TLine.h"
#include "TSystem.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TLatex.h"
#include "TRatioPlot.h"
#include "TPad.h"

#include <vector>

//________________
void plotCMSHeader(int collSystem = 0, double energy = 5.02) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Simulation}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, Form("%s #sqrt{s_{NN}} = %3.2f TeV", collSystemStr.Data(), energy) );
    t.SetTextSize(0.05);
}

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
}

//________________
void setRatioPadStyle(TPad *pad, bool isUpper = true) {
    if ( isUpper ) {
        pad->SetTopMargin(0.1);
        pad->SetBottomMargin(0.2);
        pad->SetRightMargin(0.1);
        pad->SetLeftMargin(0.2);
        pad->SetFrameLineWidth(2);
    }
    else {
        pad->SetTopMargin(0.1);
        pad->SetBottomMargin(0.2);
        pad->SetRightMargin(0.1);
        pad->SetLeftMargin(0.2);
        pad->SetFrameLineWidth(2);
    }
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

    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);    
    h->GetYaxis()->SetTitleOffset(1.2);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void setGraphStyle(TGraph *gr, Int_t type = 0, Bool_t doRenorm = kFALSE) {
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

    gr->SetLineWidth( lineWidth );
    gr->SetLineColor( color );
    
    gr->SetMarkerStyle( markerStyle );
    gr->SetMarkerColor( color );
    gr->SetMarkerSize( markerSize );

    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetNdivisions(205);
    gr->GetYaxis()->SetNdivisions(205);    
    gr->GetYaxis()->SetTitleOffset(1.1);
}

//________________
void set2DStyle(TH2* h) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(210);
    h->GetYaxis()->SetNdivisions(210);    
    h->GetYaxis()->SetTitleOffset(1.0);
}

//________________
void drawPlots(TFile *f, int collSystem = 0, double energy = 5.02) {

    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.4;
    double yTextPosition = 0.85;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );

    TH3D *hGenPtEtaPtHat = dynamic_cast<TH3D*>( f->Get("hGenInclusiveJetPtEtaPtHat") );
    TH3D *hRefPtEtaPtHat = dynamic_cast<TH3D*>( f->Get("hRefInclusiveJetPtEtaPtHat") );

    double ptRange[2] = {30., 400.};
    double etaRange[2] = {-3.6, 3.6};

    TH2D *hGenEtaPt = dynamic_cast<TH2D*>( hGenPtEtaPtHat->Project3D("yx") );
    hGenEtaPt->SetName("hGenEtaPt");
    set2DStyle( hGenEtaPt );
    TH2D *hRefEtaPt = dynamic_cast<TH2D*>( hRefPtEtaPtHat->Project3D("yx") );
    hRefEtaPt->SetName("hRefEtaPt");
    set2DStyle( hRefEtaPt );
    hRefEtaPt->Divide( hGenEtaPt );

    TCanvas *cEfficiency2D = new TCanvas( "cEfficiency2D", "cEfficiency2D", 800, 800 );
    setPadStyle();
    hRefEtaPt->Draw("colz");
    hRefEtaPt->GetXaxis()->SetRangeUser(etaRange[0], etaRange[1]);
    hRefEtaPt->GetYaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    hRefEtaPt->SetMaximum(1.0); 
    hRefEtaPt->SetMinimum(0.85);
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, "Efficiency");
    //cEfficiency2D->SaveAs("efficiency_plot.png");
}

//________________
void plotEfficiency() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    int collSystem = 0;
    double energy = 5.02;

    // Pythia for pp5020
    // TFile *f = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", uname.Data()) );
    // if ( !f ) {
    //     std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", uname.Data()) << std::endl;
    //     return;
    // }
    // collSystem = 0;
    // energy = 5.02;

    // Embedding for pPb8160
    // TFile *f = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", uname.Data()) );
    TFile *f = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pgoing_eta2.root", uname.Data()) );
    if ( !f ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pbgoing_eta2.root", uname.Data()) << std::endl;
        return;
    }
    collSystem = 1;
    energy = 8.16;

    // Plot efficiency
    drawPlots(f, collSystem, energy);
}