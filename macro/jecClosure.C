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
void plotCMSHeader() {
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Simulation}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, "PYTHIA #sqrt{s_{NN}} = 5.02 TeV");
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
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
    h->GetYaxis()->SetTitleOffset(1.0);
}

//________________
void plotComparison(TCanvas *c, TH1D* h1, TH1D* h2, 
                    int ptLow=20, int ptHatLow=15,
                    const char* h1Name="Reco", const char* h2Name="Gen",
                    bool isCM = false,
                    bool isJet = true) {

    // Text 
    TLatex t;
    TString frameT = (isCM) ? "CM" : "Lab";
    TString jetType = (isJet) ? "jet" : "dijet";
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );

    int maximumBin = h1->GetMaximumBin();
    double maximumVal = h1->GetBinContent( maximumBin );

    // Number for plotting position
    double xRange[2] = { -4.2, 4.2 };
    double yRange[2] = {0.0000001, maximumVal * 1.25 };
    double legX[2] = {0.4, 0.65};
    double legY[2] = {0.2, 0.35};
    double ratioYRange[2] = {0.6, 1.4};

    h1->GetXaxis()->SetTitle( Form("#eta^{%s}", jetType.Data() ));
    h1->GetYaxis()->SetTitle( Form("1/N_{%s} dN/d#eta^{%s}", jetType.Data(), jetType.Data() ) );
    h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);

    std::vector<double> gridLineValues { 0.95, 1.0, 1.05 };

    // Create ratio plot
    TRatioPlot *ratioPlot = new TRatioPlot( h1, h2 );
    ratioPlot->SetH1DrawOpt("E");
    ratioPlot->SetH2DrawOpt("E");
    ratioPlot->SetGridlines( gridLineValues );
    ratioPlot->Draw();

    // Set pad parameters
    ratioPlot->GetUpperPad()->SetFrameLineWidth(2);
    ratioPlot->GetLowerPad()->SetFrameLineWidth(2);
    ratioPlot->SetLeftMargin(0.15);
    ratioPlot->SetRightMargin(0.05);
    ratioPlot->SetSeparationMargin(0.03);

    // Lower plot style and titles
    ratioPlot->GetLowerRefXaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefXaxis()->SetTitleOffset(0.85);
    ratioPlot->GetLowerRefYaxis()->SetTitle( Form( "%s / %s", h1Name, h2Name ) );
    ratioPlot->GetLowerRefYaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefYaxis()->SetTitleOffset(1.2);
    

    ratioPlot->GetLowerRefGraph()->SetMarkerStyle(20);
    ratioPlot->GetLowerRefGraph()->SetMarkerSize(1.2);
    ratioPlot->GetLowerRefGraph()->SetMarkerColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineWidth(2);
    ratioPlot->GetLowYaxis()->SetNdivisions(205);

    ratioPlot->GetLowerRefGraph()->SetMinimum( ratioYRange[0] );
    ratioPlot->GetLowerRefGraph()->SetMaximum( ratioYRange[1] );
    ratioPlot->GetLowerRefXaxis()->SetRangeUser( xRange[0], xRange[1] );

    ratioPlot->GetUpperPad()->cd();
    t.DrawLatexNDC(0.5, 0.84, Form("p_{T}^{%s} > %d GeV #hat{p}_{T} > %d GeV", jetType.Data(), ptLow, ptHatLow));
    t.DrawLatexNDC(0.2, 0.84, Form("%s frame", frameT.Data() ) );
    plotCMSHeader();

    /*
    // Make ratios
    TH1D *hRatio = dynamic_cast<TH1D*>( h1->Clone("hRatio") );
    hRatio->Divide( h2 );
    set1DStyle( hRatio, 2 );

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot comparison
    //

    c->cd(1);

    setPadStyle();
    // Plot distributions
    h1->Draw();
    h2->Draw("same");
    // Set ranges
    h1->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    h1->GetYaxis()->SetTitle("dN/d#eta^{jet}");
    h1->GetXaxis()->SetTitle("#eta^{jet}");

    t.DrawLatexNDC(0.2, 0.84, Form("p_{T} > %d GeV #hat{p}_{T} > %d GeV", ptLow, ptHatLow));
    t.DrawLatexNDC(0.55, 0.75, Form("%s frame", frameT.Data() ) );
    plotCMSHeader();
    

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    leg->AddEntry(h1, Form("%s", h1Name), "p");
    leg->AddEntry(h2, Form("%s", h2Name), "p");
    leg->Draw();

    //
    // Plot ratio
    //

    c->cd(2);
    setPadStyle();
    hRatio->Draw();
    hRatio->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetYaxis()->SetTitle( Form( "%s / %s", h1Name, h2Name ) );
    hRatio->GetXaxis()->SetTitle("#eta^{jet}");
    t.DrawLatexNDC(0.2, 0.84, Form("p_{T} > %d GeV #hat{p}_{T} > %d GeV", ptLow, ptHatLow));
    t.DrawLatexNDC(0.55, 0.75, Form("%s frame", frameT.Data() ) );
    plotCMSHeader();

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    leg->Draw();

    // Line at unity
    line = new TLine(xRange[0], 1.0, xRange[1], 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kMagenta);
    line->Draw();
    */
}

//________________
void plot3DClosures(TFile *f) {

    // Retrieve histograms
    TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D*>( f->Get("hRecoInclusiveJetPtEtaPtHat") );
    TH3D *hGenPtEtaPtHat = dynamic_cast<TH3D*>( f->Get("hGenInclusiveJetPtEtaPtHat") );

    // Create vector of ptHat and jet pT bins for projecitons
    int ptHatStart = 0;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins { 1, 4, 5, 6, 7 };
    int jetPtStart = 5;
    int jetPtStep = 5;  // Starting from 5 GeV: jetPtStart + (jetPtBins(i) - 1) * jetPtStep
    int jetPtBinsMax = 150;
    std::vector<int> jetPtBins { 1, 6, 10, 14, 18 };

    // Declare canvases and histograms
    TCanvas *cClosure[ ptHatBins.size() ][ jetPtBins.size() ];
    TH1D *hRecoEta[ ptHatBins.size() ][ jetPtBins.size() ];
    TH1D *hGenEta[ ptHatBins.size() ][ jetPtBins.size() ];

    // Loop over ptHat and jet pT bins
    for (unsigned int i = 0; i < ptHatBins.size(); i++) {
        for (unsigned int j = 0; j < jetPtBins.size(); j++) {

            // Create canvas
            cClosure[i][j] = new TCanvas( Form("cClosure_%d_%d", i, j), Form("cClosure_%d_%d", i, j), 700, 800 );
            // Make projection of the 3D histograms
            hRecoEta[i][j] = dynamic_cast<TH1D*>( hRecoPtEtaPtHat->ProjectionX( Form("hRecoEta_%d_%d", i, j), 
                                                                               jetPtBins[j], jetPtBinsMax, 
                                                                               ptHatBins[i], ptHatBinsMax ) );
            set1DStyle( hRecoEta[i][j], 0, kTRUE );
            hGenEta[i][j] = dynamic_cast<TH1D*>( hGenPtEtaPtHat->ProjectionX( Form("hGenEta_%d_%d", i, j), 
                                                                              jetPtBins[j], jetPtBinsMax, 
                                                                              ptHatBins[i], ptHatBinsMax ) );
            set1DStyle( hGenEta[i][j], 1, kTRUE );

            // Plot distributions
            plotComparison(cClosure[i][j], hRecoEta[i][j], hGenEta[i][j], 
                           jetPtStart + (jetPtBins[j] - 1) * jetPtStep, 
                           ptHatStart + (ptHatBins[i] - 1) * ptHatStep, 
                           "Reco", "Gen", false, true);
                           
        } // for (unsigned int j = 0; j < jetPtBins.size(); j++)
    } // for (unsigned int i = 0; i < ptHatBins.size(); i++)

}

//________________
void plotSimple1DClosure(TFile *f) {

    // Lab frame
    TH1D *hGenEtaLab  = dynamic_cast<TH1D*>( f->Get("hGenGoodInclusiveJetEtaLabFrame") );
    TH1D *hRecoEtaLab = dynamic_cast<TH1D*>( f->Get("hRecoGoodInclusiveJetEtaLabFrame") );

    set1DStyle( hGenEtaLab, 1, kTRUE );
    set1DStyle( hRecoEtaLab, 0, kTRUE );

    TCanvas *cSimpleJESLab = new TCanvas( "cSimpleJESLab", "cSimpleJESLab", 700, 800 );
    plotComparison(cSimpleJESLab, hRecoEtaLab, hGenEtaLab, 30, 15, "Reco", "Gen", false);

    // CM frame
    TH1D *hGenEtaCM  = dynamic_cast<TH1D*>( f->Get("hGenGoodInclusiveJetEtaCMFrame") );
    TH1D *hRecoEtaCM = dynamic_cast<TH1D*>( f->Get("hRecoGoodInclusiveJetEtaCMFrame") );

    set1DStyle( hGenEtaCM, 1, kTRUE );
    set1DStyle( hRecoEtaCM, 0, kTRUE );

    TCanvas *cSimpleJESCM = new TCanvas( "cSimpleJESCM", "cSimpleJESCM", 700, 800 );
    plotComparison(cSimpleJESCM, hRecoEtaCM, hGenEtaCM, 30, 15, "Reco", "Gen", true);
}

//________________
void jecClosure() {
    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    // Pythia for pp5020
    TFile *pp5020PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", uname.Data()) );
    if ( !pp5020PythiaFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", uname.Data()) << std::endl;
        return;
    }

    // Embedding for pPb8160
    TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", uname.Data()) );
    if ( !pPb8160EmbedFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", uname.Data()) << std::endl;
        return;
    }

    // Plot simple JES
    // plotSimple1DClosure(pp5020PythiaFile);

    // Plot 3D JES closures
    //plot3DClosures(pp5020PythiaFile);

    plot3DClosures(pPb8160EmbedFile);
}