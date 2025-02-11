// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TRatioPlot.h"
#include "TPad.h"

// C++ headers
#include <iostream>
#include <vector>
#include <sys/stat.h>

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
void plotCMSHeader(int collSystem = 0, double energy = 5.02) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, Form("%s #sqrt{s_{NN}} = %3.2f TeV", collSystemStr.Data(), energy) );
    t.SetTextSize(0.05);
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
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
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
    }
    h->Scale( 1. / h->Integral() );
}

//________________
void plotComparison(TCanvas *c, TH1D* h1, TH1D* h2, 
                    int ptLow=50, int ptHi=60,
                    const char* h1Name="Reco", const char* h2Name="Gen",
                    int collSystem = 0, double energy = 5.02,
                    bool isRpPb = false) {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // Energy in TeV

    // Text 
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );

    int maximumBin = h1->GetMaximumBin();
    double maximumVal = h1->GetBinContent( maximumBin );

    // Number for plotting position
    double xRange[2] = { -4.2, 4.2 };
    // double yRange[2] = {0.0000001, maximumVal * 1.25 };
    double yRange[2] = {0.0000001, 0.14 };
    if ( isRpPb ) {
        yRange[0] = 0.7; yRange[1] = 1.3;
    }
    double legX[2] = {0.2, 0.4};
    double legY[2] = {0.55, 0.7};
    double ratioYRange[2] = {0.8, 1.2};

    // h1->GetXaxis()->SetTitle( Form("#eta^{dijet}_{%s}", frameT.Data() ));
    // h1->GetYaxis()->SetTitle( Form("1/N_{dijet} dN/d#eta^{dijet}_{%s}", frameT.Data() ) );
    if ( isRpPb ) {
        h1->GetYaxis()->SetTitle( Form("R_{pPb}}") );
    }
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
    t.DrawLatexNDC(0.4, 0.84, Form("%d< p_{T}^{ave} (GeV) < %d ", ptLow, ptHi));
    plotCMSHeader(collSystem, energy);
    TLegend *leg = new TLegend( legX[0], legY[0], legX[1], legY[1] );
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->AddEntry( h1, h1Name, "p" );
    leg->AddEntry( h2, h2Name, "p" );
    leg->Draw();
}

//________________
void compareReco2GenInclusiveJetPtEta(TFile *f) {

    // Bins to project 
    int etaBinsProj[2] = {20, 32};
    int ptBinsProj[2] = {1, 1};
    TH2D *hGenInclusiveJetPtEta = dynamic_cast<TH2D*>( f->Get("hGenInclusiveJetPtEta") );
    TH2D *hRecoInclusiveJetPtEta = dynamic_cast<TH2D*>( f->Get("hRecoInclusiveAllJetPtVsEta") );

    TH1D *hRecoInclusiveJetPt = dynamic_cast<TH1D*>( hRecoInclusiveJetPtEta->ProjectionY( "hRecoInclusiveJetPt", etaBinsProj[0], etaBinsProj[1] ) );
    TH1D *hGenInclusiveJetPt = dynamic_cast<TH1D*>( hGenInclusiveJetPtEta->ProjectionY( "hGenInclusiveJetPt", etaBinsProj[0], etaBinsProj[1] ) );

    // hRecoInclusiveJetPt->Scale( 1./hRecoInclusiveJetPt->Integral() );
    set1DStyle( hRecoInclusiveJetPt, 0 );
    //hGenInclusiveJetPt->Scale( 1./hGenInclusiveJetPt->Integral() );
    set1DStyle( hGenInclusiveJetPt, 1 );

    int canvX{500}, canvY{1000};
    TCanvas *cReco2GenComp = new TCanvas( "cReco2GenComp", "cReco2GenComp", canvX, canvY );
    cReco2GenComp->Divide(1, 2);
    plotComparison(cReco2GenComp, hRecoInclusiveJetPt, hGenInclusiveJetPt, 0, 1500, "Reco", "Gen");
}


//________________
void rescaleForwardBackward(TH1D *hForward, TH1D *hBackward) {
    // Rescale forward and backward histograms
    // to have the same integral
    double intForward = hForward->Integral();
    double intBackward = hBackward->Integral();
    double scaleFactor = intForward + intBackward;
    hForward->Scale( intForward / scaleFactor );
    hBackward->Scale( intBackward / scaleFactor );
}

//________________
void comparisons2gen(TFile *f, int collisionSystem = 0, double collisionEnergy = 5.02, TString date = "20250129") {
    
    // collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    // Dijet ptAve binning
    int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
                             100, 110,  120, 130, 140,
                             150, 160,  180, 200, 250, 
                             300, 500};
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);

    std::vector<int> dijetPtVals; 
    dijetPtVals.assign(dijetPtNewVals, dijetPtNewVals + sizeOfPtVals);

    std::cout << "Number of pT bins: " << dijetPtVals.size() << std::endl;

    // Bins for projections from 3D
    std::vector<int> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55 };
    std::vector<int> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94 };

    TH1D *h1DProj[ ptDijetBinLow.size() ];
    TH1D *h1DDirect[ ptDijetBinLow.size() ];
    TCanvas *cComp[ ptDijetBinLow.size() ];

    TCanvas *cDijetEtaLabComp[ ptDijetBinLow.size() ];
    TCanvas *cDijetEtaLabRat[ ptDijetBinLow.size() ];
    TCanvas *cDijetEtaCMComp[ ptDijetBinLow.size() ];
    TCanvas *cDijetEtaCMRat[ ptDijetBinLow.size() ];
    TCanvas *cDijetEtaFB[ ptDijetBinLow.size() ];

    //
    // Lab frame
    //
    TH1D *hRecoDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hReco2GenDijetEta1DLab[ ptDijetBinLow.size() ];
    TH1D *hRef2GenDijetEta1DLab[ ptDijetBinLow.size() ];

    //
    // CM frame
    //
    TH1D *hRecoDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hReco2GenDijetEta1DCM[ ptDijetBinLow.size() ];
    TH1D *hRef2GenDijetEta1DCM[ ptDijetBinLow.size() ];

    //
    // Forward/backward ratios
    //

    TH1D *hRecoDijetEtaCMForward1D[ ptDijetBinLow.size() ];
    TH1D *hRecoDijetEtaCMBackward1D[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEtaCMForward1D[ ptDijetBinLow.size() ];
    TH1D *hGenDijetEtaCMBackward1D[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEtaCMForward1D[ ptDijetBinLow.size() ];
    TH1D *hRefDijetEtaCMBackward1D[ ptDijetBinLow.size() ];

    TH1D *hRecoDijetFBEtaCM1D[ ptDijetBinLow.size() ];
    TH1D *hGenDijetFBEtaCM1D[ ptDijetBinLow.size() ];
    TH1D *hRefDijetFBEtaCM1D[ ptDijetBinLow.size() ];

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );
    TLegend *leg;

    // Loop over dijet ptAve bins
    for (unsigned int i = 0; i < ptDijetBinLow.size(); i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{1000}, canvY{1000};

        //
        // Lab frame
        //
        hRecoDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1DWeighted_%d", i) ) );
        hRecoDijetEta1DLab[i]->SetName( Form("hRecoDijetEta1DLab_%d", i) );
        set1DStyle( hRecoDijetEta1DLab[i], 0 );
        rescaleEta( hRecoDijetEta1DLab[i] );
        hGenDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1DWeighted_%d", i) ) );
        hGenDijetEta1DLab[i]->SetName( Form("hGenDijetEta1DLab_%d", i) );
        set1DStyle( hGenDijetEta1DLab[i], 2 );
        rescaleEta( hGenDijetEta1DLab[i] );
        hRefDijetEta1DLab[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1DWeighted_%d", i) ) );
        hRefDijetEta1DLab[i]->SetName( Form("hRefDijetEta1DLab_%d", i) );
        set1DStyle( hRefDijetEta1DLab[i], 1 );
        rescaleEta( hRefDijetEta1DLab[i] );

        hReco2GenDijetEta1DLab[i] = dynamic_cast<TH1D*>( hRecoDijetEta1DLab[i]->Clone( Form("hReco2GenDijetEta1DLab_%d", i) ) );
        hReco2GenDijetEta1DLab[i]->Divide( hGenDijetEta1DLab[i] );
        hRef2GenDijetEta1DLab[i] = dynamic_cast<TH1D*>( hRefDijetEta1DLab[i]->Clone( Form("hRef2GenDijetEta1DLab_%d", i) ) );
        hRef2GenDijetEta1DLab[i]->Divide( hGenDijetEta1DLab[i] );

        //
        // CM frame
        //
        hRecoDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEta1DCMWeighted_%d", i) ) );
        hRecoDijetEta1DCM[i]->SetName( Form("hRecoDijetEta1DCM_%d", i) );
        set1DStyle( hRecoDijetEta1DCM[i], 0 );
        rescaleEta( hRecoDijetEta1DCM[i] );
        hGenDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEta1DCMWeighted_%d", i) ) );
        hGenDijetEta1DCM[i]->SetName( Form("hGenDijetEta1DCM_%d", i) );
        set1DStyle( hGenDijetEta1DCM[i], 2 );
        rescaleEta( hGenDijetEta1DCM[i] );
        hRefDijetEta1DCM[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEta1DCMWeighted_%d", i) ) );
        hRefDijetEta1DCM[i]->SetName( Form("hRefDijetEta1DCM_%d", i) );
        set1DStyle( hRefDijetEta1DCM[i], 1 );
        rescaleEta( hRefDijetEta1DCM[i] );
        hReco2GenDijetEta1DCM[i] = dynamic_cast<TH1D*>( hRecoDijetEta1DCM[i]->Clone( Form("hReco2GenDijetEta1DCM_%d", i) ) );
        hReco2GenDijetEta1DCM[i]->Divide( hGenDijetEta1DCM[i] );
        hRef2GenDijetEta1DCM[i] = dynamic_cast<TH1D*>( hRefDijetEta1DCM[i]->Clone( Form("hRef2GenDijetEta1DCM_%d", i) ) );
        hRef2GenDijetEta1DCM[i]->Divide( hGenDijetEta1DCM[i] );

        //
        // Forward/backward ratios
        //
        hRecoDijetEtaCMForward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMForward1DWeighted_%d", i) ) );
        hRecoDijetEtaCMForward1D[i]->SetName( Form("hRecoDijetEtaCMForward1D_%d", i) );
        hRecoDijetEtaCMBackward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRecoDijetEtaCMBackward1DWeighted_%d", i) ) );
        hRecoDijetEtaCMBackward1D[i]->SetName( Form("hRecoDijetEtaCMBackward1D_%d", i) );
        hGenDijetEtaCMForward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMForward1DWeighted_%d", i) ) );
        hGenDijetEtaCMForward1D[i]->SetName( Form("hGenDijetEtaCMForward1D_%d", i) );
        hGenDijetEtaCMBackward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hGenDijetEtaCMBackward1DWeighted_%d", i) ) );
        hGenDijetEtaCMBackward1D[i]->SetName( Form("hGenDijetEtaCMBackward1D_%d", i) );
        hRefDijetEtaCMForward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMForward1DWeighted_%d", i) ) );
        hRefDijetEtaCMForward1D[i]->SetName( Form("hRefDijetEtaCMForward1D_%d", i) );
        hRefDijetEtaCMBackward1D[i] = dynamic_cast<TH1D*>( f->Get( Form("hRefDijetEtaCMBackward1DWeighted_%d", i) ) );
        hRefDijetEtaCMBackward1D[i]->SetName( Form("hRefDijetEtaCMBackward1D_%d", i) );

        rescaleForwardBackward( hRecoDijetEtaCMForward1D[i], hRecoDijetEtaCMBackward1D[i] );
        rescaleForwardBackward( hGenDijetEtaCMForward1D[i], hGenDijetEtaCMBackward1D[i] );
        rescaleForwardBackward( hRefDijetEtaCMForward1D[i], hRefDijetEtaCMBackward1D[i] );
        set1DStyle( hRecoDijetEtaCMForward1D[i], 0 );
        set1DStyle( hGenDijetEtaCMForward1D[i], 2 );
        set1DStyle( hRefDijetEtaCMForward1D[i], 1 );

        hRecoDijetFBEtaCM1D[i] = dynamic_cast<TH1D*>( hRecoDijetEtaCMForward1D[i]->Clone( Form("hRecoDijetFBEtaCM1D_%d", i) ) );
        hRecoDijetFBEtaCM1D[i]->Divide( hRecoDijetEtaCMBackward1D[i] );
        hRecoDijetFBEtaCM1D[i]->GetYaxis()->SetTitle("Forward/Backward");
        hGenDijetFBEtaCM1D[i] = dynamic_cast<TH1D*>( hGenDijetEtaCMForward1D[i]->Clone( Form("hGenDijetFBEtaCM1D_%d", i) ) );
        hGenDijetFBEtaCM1D[i]->Divide( hGenDijetEtaCMBackward1D[i] );
        hGenDijetFBEtaCM1D[i]->GetYaxis()->SetTitle("Forward/Backward");
        hRefDijetFBEtaCM1D[i] = dynamic_cast<TH1D*>( hRefDijetEtaCMForward1D[i]->Clone( Form("hRefDijetFBEtaCM1D_%d", i) ) );
        hRefDijetFBEtaCM1D[i]->Divide( hRefDijetEtaCMBackward1D[i] );
        hRefDijetFBEtaCM1D[i]->GetYaxis()->SetTitle("Forward/Backward");


        double xRange[2] = { -3.2, 3.2 };
        double yRange[2] = {0.0000001, 0.14 };

        //
        // Plot comparisons in the lab frame
        //
        cDijetEtaLabComp[i] = new TCanvas( Form("cDijetEtaLabComp_%d", i), 
                                           Form("cDijetEtaLabComp_%d", i), 
                                           canvX, canvY );
        setPadStyle();
        hRecoDijetEta1DLab[i]->Draw();
        hGenDijetEta1DLab[i]->Draw("same");
        hRefDijetEta1DLab[i]->Draw("same");
        hRecoDijetEta1DLab[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hRecoDijetEta1DLab[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        plotCMSHeader(collisionSystem, collisionEnergy);        
        t.DrawLatexNDC(0.4, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", dijetPtVals[i], dijetPtVals[i+1]) );
        leg = new TLegend(0.2, 0.55, 0.4, 0.7);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->AddEntry( hRecoDijetEta1DLab[i], "Reco", "p" );
        leg->AddEntry( hGenDijetEta1DLab[i], "Gen", "p" );
        leg->AddEntry( hRefDijetEta1DLab[i], "Ref", "p" );        
        leg->Draw();

        //
        // Plot ratios in the lab frame
        //
        cDijetEtaLabRat[i] = new TCanvas( Form("cDijetEtaLabRat_%d", i), 
                                          Form("cDijetEtaLabRat_%d", i), 
                                          canvX, canvY );
        setPadStyle();
        hReco2GenDijetEta1DLab[i]->Draw();
        hRef2GenDijetEta1DLab[i]->Draw("same");
        hReco2GenDijetEta1DLab[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hReco2GenDijetEta1DLab[i]->GetYaxis()->SetRangeUser(0.6, 1.4);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", dijetPtVals[i], dijetPtVals[i+1]) );
        leg = new TLegend(0.4, 0.2, 0.6, 0.4);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->AddEntry( hReco2GenDijetEta1DLab[i], "Reco/Gen", "p" );
        leg->AddEntry( hRef2GenDijetEta1DLab[i], "Ref/Gen", "p" );
        leg->Draw();

        //
        // Plot comparisons in the CM frame
        //

        xRange[0] = -2.5; xRange[1] = 2.5;

        cDijetEtaCMComp[i] = new TCanvas( Form("cDijetEtaCMComp_%d", i), 
                                           Form("cDijetEtaCMComp_%d", i), 
                                           canvX, canvY );
        setPadStyle();
        hRecoDijetEta1DCM[i]->Draw();
        hGenDijetEta1DCM[i]->Draw("same");
        hRefDijetEta1DCM[i]->Draw("same");
        hRecoDijetEta1DCM[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hRecoDijetEta1DCM[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", dijetPtVals[i], dijetPtVals[i+1]) );
        leg = new TLegend(0.2, 0.55, 0.4, 0.7);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->AddEntry( hRecoDijetEta1DCM[i], "Reco", "p" );
        leg->AddEntry( hGenDijetEta1DCM[i], "Gen", "p" );
        leg->AddEntry( hRefDijetEta1DCM[i], "Ref", "p" );
        leg->Draw();

        //
        // Plot ratios in the CM frame
        //
        cDijetEtaCMRat[i] = new TCanvas( Form("cDijetEtaCMRat_%d", i), 
                                          Form("cDijetEtaCMRat_%d", i), 
                                          canvX, canvY );

        setPadStyle();
        hReco2GenDijetEta1DCM[i]->Draw();
        hRef2GenDijetEta1DCM[i]->Draw("same");
        hReco2GenDijetEta1DCM[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hReco2GenDijetEta1DCM[i]->GetYaxis()->SetRangeUser(0.6, 1.4);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", dijetPtVals[i], dijetPtVals[i+1]) );
        leg = new TLegend(0.4, 0.2, 0.6, 0.4);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->AddEntry( hReco2GenDijetEta1DCM[i], "Reco/Gen", "p" );
        leg->AddEntry( hRef2GenDijetEta1DCM[i], "Ref/Gen", "p" );
        leg->Draw();

        //
        // Forward/backward ratios
        //
        cDijetEtaFB[i] = new TCanvas( Form("cDijetEtaFB_%d", i), 
                                      Form("cDijetEtaFB_%d", i), 
                                      canvX, canvY );
        setPadStyle();
        hRecoDijetFBEtaCM1D[i]->Draw();
        hGenDijetFBEtaCM1D[i]->Draw("same");
        hRefDijetFBEtaCM1D[i]->Draw("same");
        hRecoDijetFBEtaCM1D[i]->GetXaxis()->SetRangeUser(0, 2.5);
        hRecoDijetFBEtaCM1D[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        plotCMSHeader(collisionSystem, collisionEnergy);
        gPad->SetGrid();
        t.DrawLatexNDC(0.4, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", dijetPtVals[i], dijetPtVals[i+1]) );
        leg = new TLegend(0.4, 0.2, 0.6, 0.4);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->AddEntry( hRecoDijetFBEtaCM1D[i], "Reco", "p" );
        leg->AddEntry( hGenDijetFBEtaCM1D[i], "Gen", "p" );
        leg->AddEntry( hRefDijetFBEtaCM1D[i], "Ref", "p" );
        leg->Draw();        

    } // end loop over dijet pt bins
}



//________________
void plotMcClosures() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Date extraction
    TDatime dt; 
    TString date { Form( "%d",dt.GetDate() ) };

    if ( !directoryExists( date.Data() ) ) {
        createDirectory( date.Data() );
    } 


    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    int collisionSystem = 0;  // 0 - pp, 1 - pPb, 2 - pPb5020, 3 - pPb8160
    double collisionEnergy = 5.02;  // 5.02 TeV

    //
    // pPb8160
    //

    collisionSystem = 1;
    collisionEnergy = 8.16;

    // MC p-going direction new (coincides with the pPb5020)
    TFile *pPb8160EmbedNewFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_eta25.root", uname.Data()) );
    if ( !pPb8160EmbedNewFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pPb8160_pgoing_new.root", uname.Data()) << std::endl;
        return;
    }


    comparisons2gen( pPb8160EmbedNewFile, collisionSystem, collisionEnergy, date );



}
