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
void drawJecComparison(TCanvas *c, TH1D* h1, TH1D* h2, 
    int ptLow=20, int ptHatLow=15,
    const char* h1Name="Reco", const char* h2Name="Gen",
    int collSystem = 0, double energy = 5.02,
    bool isCM = false,
    bool isJet = true,
    bool isPt = false,
    double etaLow = -1.6,
    double etaHi = 1.6) {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // Energy in TeV

    // Text
    TLatex t;
    TString frameT = (isCM) ? "CM" : "Lab";
    TString jetType = (isJet) ? "jet" : "dijet";
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    int maximumBin = h1->GetMaximumBin();
    double maximumVal = h1->GetBinContent(maximumBin);

    // Number for plotting position
    double xRange[2] = {-3.2, 3.2};
    if (isPt)
    {
        xRange[0] = 0;
        xRange[1] = 500;
    }
    double yRange[2] = {0.0000001, maximumVal * 1.25};
    double legX[2] = {0.5, 0.7};
    double legY[2] = {0.2, 0.35};
    double ratioYRange[2] = {0.8, 1.2};
    if (isPt)
    {
        ratioYRange[0] = 0.9;
        ratioYRange[1] = 1.3;
    }

    if (!isPt)
    {
        h1->GetXaxis()->SetTitle(Form("#eta^{%s}", jetType.Data()));
        h1->GetYaxis()->SetTitle(Form("1/N_{%s} dN/d#eta^{%s}", jetType.Data(), jetType.Data()));
        h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    }
    else
    {
        h1->GetXaxis()->SetTitle(Form("p_{T}^{%s} (GeV)", jetType.Data()));
        h1->GetYaxis()->SetTitle(Form("1/N_{%s} dN/dp_{T}^{%s}", jetType.Data(), jetType.Data()));
    }

    std::vector<double> gridLineValues{0.95, 1.0, 1.05};

    // Create ratio plot
    TRatioPlot *ratioPlot = new TRatioPlot(h1, h2);
    ratioPlot->SetH1DrawOpt("E");
    ratioPlot->SetH2DrawOpt("E");
    ratioPlot->SetGridlines(gridLineValues);
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
    ratioPlot->GetLowerRefYaxis()->SetTitle(Form("%s / %s", h1Name, h2Name));
    ratioPlot->GetLowerRefYaxis()->SetTitleSize(0.05);
    ratioPlot->GetLowerRefYaxis()->SetTitleOffset(1.2);

    ratioPlot->GetLowerRefGraph()->SetMarkerStyle(20);
    ratioPlot->GetLowerRefGraph()->SetMarkerSize(1.2);
    ratioPlot->GetLowerRefGraph()->SetMarkerColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineColor(kBlack);
    ratioPlot->GetLowerRefGraph()->SetLineWidth(2);
    ratioPlot->GetLowYaxis()->SetNdivisions(205);

    ratioPlot->GetLowerRefGraph()->SetMinimum(ratioYRange[0]);
    ratioPlot->GetLowerRefGraph()->SetMaximum(ratioYRange[1]);
    ratioPlot->GetLowerRefXaxis()->SetRangeUser(xRange[0], xRange[1]);

    ratioPlot->GetUpperPad()->cd();
    if (!isPt)
    {
        t.DrawLatexNDC(0.5, 0.84, Form("p_{T}^{%s} > %d GeV #hat{p}_{T} > %d GeV", jetType.Data(), ptLow, ptHatLow));
    }
    else
    {
        t.DrawLatexNDC(0.5, 0.84, Form("%2.1f< #eta^{jet} < %2.1f #hat{p}_{T} > %d GeV", etaLow, etaHi, ptHatLow));
    }
    t.DrawLatexNDC(0.2, 0.84, Form("%s frame", frameT.Data()));
    plotCMSHeader(collSystem, energy);
    if (isPt)
    {
        ratioPlot->GetUpperPad()->SetLogy();
    }

    TLegend *leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h1, h1Name, "p");
    leg->AddEntry(h2, h2Name, "p");
    leg->Draw();
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
void drawToGenComparison(TCanvas *c, TH1D *hReco, TH1D *hRef, TH1D *hGen, 
                        int ptLow = 50, int ptHi = 60, 
                        bool isCM = false, bool isFB = false,
                        int collisionSystem = 1, double energy = 8.16) {

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.05 );
    TLegend *leg;

    double xRange[2] = { -3.2, 3.2 };
    double yRange[2] = {0.0000001, 0.08 };

    if ( isCM ) {
        if ( !isFB ) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] = 0.0000001; yRange[1] = 0.08;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.5;
            yRange[0] = 0.8; yRange[1] = 1.2;
        }
    }

    c->cd();
    setPadStyle();
    hReco->Draw();
    hGen->Draw("same");
    hRef->Draw("same");
    hReco->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hReco->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);        
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi) );
    if ( isFB ) {
        leg = new TLegend(0.2, 0.25, 0.4, 0.4);
    }
    else {
        leg = new TLegend(0.2, 0.55, 0.4, 0.7);
    }
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->AddEntry( hReco, "Reco", "p" );
    leg->AddEntry( hGen, "Gen", "p" );
    leg->AddEntry( hRef, "Ref", "p" );        
    leg->Draw();
}

//________________
void drawToGenRatio(TCanvas *c, TH1D *hReco2Gen, TH1D *hRef2Gen, 
                    int ptLow = 50, int ptHi = 60, 
                    bool isCM = false, bool isFB = false,
                    int collisionSystem = 1, double energy = 8.16) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TLegend *leg;

    double xRange[2] = {-3.2, 3.2};
    double yRange[2] = {0.7, 1.3};

    if (isCM) {
        if (!isFB) {
            xRange[0] = -2.5; xRange[1] = 2.5;
            yRange[0] =0.7; yRange[1] = 1.3;
        }
        else {
            xRange[0] = 0; xRange[1] = 2.5;
            yRange[0] = 0.7; yRange[1] = 1.3;
        }
    }

    c->cd();
    setPadStyle();
    hReco2Gen->Draw();
    hRef2Gen->Draw("same");
    hReco2Gen->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hReco2Gen->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    gPad->SetGrid();
    plotCMSHeader(collisionSystem, energy);
    t.DrawLatexNDC(0.35, 0.84, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi));
    leg = new TLegend(0.2, 0.25, 0.4, 0.4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->AddEntry(hReco2Gen, "Reco/Gen", "p");
    leg->AddEntry(hRef2Gen, "Ref/Gen", "p");
    leg->Draw();
}

//________________
void comparisons2gen(TFile *f, int collisionSystem = 0, double collisionEnergy = 5.02, TString date = "20250129") {
    
    // collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

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
    TCanvas *c = new TCanvas( "c", "c", 1000, 1000 );

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
        hReco2GenDijetEta1DLab[i]->Divide( hReco2GenDijetEta1DLab[i], hGenDijetEta1DLab[i], 1., 1., "b" );
        hRef2GenDijetEta1DLab[i] = dynamic_cast<TH1D*>( hRefDijetEta1DLab[i]->Clone( Form("hRef2GenDijetEta1DLab_%d", i) ) );
        hRef2GenDijetEta1DLab[i]->Divide( hRef2GenDijetEta1DLab[i], hGenDijetEta1DLab[i], 1., 1., "b" );

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
        hReco2GenDijetEta1DCM[i]->Divide( hReco2GenDijetEta1DCM[i], hGenDijetEta1DCM[i], 1., 1., "b" );
        hRef2GenDijetEta1DCM[i] = dynamic_cast<TH1D*>( hRefDijetEta1DCM[i]->Clone( Form("hRef2GenDijetEta1DCM_%d", i) ) );
        hRef2GenDijetEta1DCM[i]->Divide( hRef2GenDijetEta1DCM[i], hGenDijetEta1DCM[i], 1., 1., "b" );

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

        //
        // Plot comparisons in the lab frame
        //
        drawToGenComparison(c, hRecoDijetEta1DLab[i], hRefDijetEta1DLab[i], hGenDijetEta1DLab[i],
                            dijetPtVals[i], dijetPtVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaLab_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Plot ratios in the lab frame
        //
        drawToGenRatio(c, hReco2GenDijetEta1DLab[i], hRef2GenDijetEta1DLab[i],
                       dijetPtVals[i], dijetPtVals[i+1], false, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaLab_RecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Plot comparisons in the CM frame
        //

        drawToGenComparison(c, hRecoDijetEta1DCM[i], hRefDijetEta1DCM[i], hGenDijetEta1DCM[i],
                            dijetPtVals[i], dijetPtVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaCM_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Plot ratios in the CM frame
        //
        drawToGenRatio(c, hReco2GenDijetEta1DCM[i], hRef2GenDijetEta1DCM[i],
                       dijetPtVals[i], dijetPtVals[i+1], true, false, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaCM_RecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );

        //
        // Forward/backward ratios
        //
        drawToGenComparison(c, hRecoDijetFBEtaCM1D[i], hRefDijetFBEtaCM1D[i], hGenDijetFBEtaCM1D[i],
                            dijetPtVals[i], dijetPtVals[i+1], true, true, collisionSystem, collisionEnergy);
        c->SaveAs( Form("%s/%s_dijetEtaFB_RecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), dijetPtVals[i], dijetPtVals[i+1]) );
            
    } // end loop over dijet pt bins
}

//________________
void plotInclusiveJetJECClosures(TFile *f, int collSystem = 0, double energy = 5.02) {
    // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.6;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Retrieve histograms
    TH3D *hRecoPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRecoInclusiveJetPtEtaPtHat"));
    TH3D *hGenPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hGenInclusiveJetPtEtaPtHat"));
    TH3D *hRefPtEtaPtHat = dynamic_cast<TH3D *>(f->Get("hRefInclusiveJetPtEtaPtHat"));

    TH2D *hRecoPtVsPtHat = dynamic_cast<TH2D *>(hRecoPtEtaPtHat->Project3D("yz"));
    hRecoPtVsPtHat->SetName("hRecoPtVsPtHat");
    set2DStyle(hRecoPtVsPtHat);
    TH2D *hGenPtVsPtHat = dynamic_cast<TH2D *>(hGenPtEtaPtHat->Project3D("yz"));
    hGenPtVsPtHat->SetName("hGenPtVsPtHat");
    set2DStyle(hGenPtVsPtHat);

    double ptRange[2] = {0, 250};
    TCanvas *cJetPtVsPtHat = new TCanvas("cJetPtVsPtHat", "cJetPtVsPtHat", 1000, 500);
    cJetPtVsPtHat->Divide(2, 1);
    cJetPtVsPtHat->cd(1);
    setPadStyle();
    hRecoPtVsPtHat->Draw("colz");
    hRecoPtVsPtHat->GetXaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    hRecoPtVsPtHat->GetYaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    // gPad->SetLogz();
    hRecoPtVsPtHat->GetXaxis()->SetTitle("#hat{p}_{T} (GeV)");
    hRecoPtVsPtHat->GetYaxis()->SetTitle("Reco p_{T}^{jet} (GeV)");
    plotCMSHeader(collSystem, energy);
    cJetPtVsPtHat->cd(2);
    setPadStyle();
    hGenPtVsPtHat->Draw("colz");
    hGenPtVsPtHat->GetXaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    hGenPtVsPtHat->GetYaxis()->SetRangeUser(ptRange[0], ptRange[1]);
    // gPad->SetLogz();
    hGenPtVsPtHat->GetXaxis()->SetTitle("#hat{p}_{T} (GeV)");
    hGenPtVsPtHat->GetYaxis()->SetTitle("Gen p_{T}^{jet} (GeV)");
    plotCMSHeader(collSystem, energy);

    // Create vector of ptHat and jet pT bins for projections
    int ptHatStart = 0;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins{1}; // 20

    // Jet pT binning
    int jetPtStart = 5;
    int jetPtStep = 10; // Starting from 5 GeV: jetPtStart + (jetPtBins(i) - 1) * jetPtStep
    int jetPtBinsMax = 150;
    std::vector<int> jetPtBins{4, 7, 12}; // 35, 55, 105

    // Eta binning
    // 52 bins from (-5.2, 5.2)
    int nEtaBins = 52;
    double etaStep = 0.2;
    std::vector<int> jetEtaBinsLow{19, 12, 38};
    std::vector<int> jetEtaBinsHigh{35, 16, 42};

    // Declare canvases and histograms
    TCanvas *cPtVsEta[ptHatBins.size()];
    TCanvas *cClosureEta[ptHatBins.size()][jetPtBins.size()];
    TCanvas *cClosurePt[ptHatBins.size()][jetEtaBinsLow.size()];

    TH2D *hRecoPtVsEta[ptHatBins.size()];
    TH1D *hRecoEta[ptHatBins.size()][jetPtBins.size()];
    TH2D *hGenPtVsEta[ptHatBins.size()];
    TH1D *hGenEta[ptHatBins.size()][jetPtBins.size()];
    TH1D *hRecoPt[ptHatBins.size()][jetEtaBinsLow.size()];
    TH1D *hGenPt[ptHatBins.size()][jetEtaBinsLow.size()];

    // Loop over ptHat and jet pT bins
    for (unsigned int i = 0; i < ptHatBins.size(); i++) {

        hRecoPtEtaPtHat->GetZaxis()->SetRange(ptHatBins[i], ptHatBinsMax);
        hGenPtEtaPtHat->GetZaxis()->SetRange(ptHatBins[i], ptHatBinsMax);

        // Make 2D histograms
        hRecoPtVsEta[i] = dynamic_cast<TH2D *>(hRecoPtEtaPtHat->Project3D("yx"));
        hRecoPtVsEta[i]->SetName(Form("hRecoPtVsEta_%d", i));
        set2DStyle(hRecoPtVsEta[i]);

        hGenPtVsEta[i] = dynamic_cast<TH2D *>(hGenPtEtaPtHat->Project3D("yx"));
        hGenPtVsEta[i]->SetName(Form("hGenPtVsEta_%d", i));
        set2DStyle(hGenPtVsEta[i]);

        cPtVsEta[i] = new TCanvas(Form("cPtVsEta_%d", i), Form("cPtVsEta_%d", i), 1000, 500);
        cPtVsEta[i]->Divide(2, 1);
        cPtVsEta[i]->cd(1);
        setPadStyle();
        hRecoPtVsEta[i]->Draw("colz");
        hRecoPtVsEta[i]->GetXaxis()->SetRangeUser(-5.2, 5.2);
        hRecoPtVsEta[i]->GetYaxis()->SetRangeUser(0, 120);
        hRecoPtVsEta[i]->GetXaxis()->SetTitle("Reco #eta^{jet}");
        hRecoPtVsEta[i]->GetYaxis()->SetTitle("Reco p_{T}^{jet} (GeV)");
        plotCMSHeader(collSystem, energy);
        gPad->SetLogz();
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("#hat{p}_{T} > %d GeV", ptHatStart + (ptHatBins[i] - 1) * ptHatStep));
        cPtVsEta[i]->cd(2);
        setPadStyle();
        hGenPtVsEta[i]->Draw("colz");
        hGenPtVsEta[i]->GetXaxis()->SetRangeUser(-5.2, 5.2);
        hGenPtVsEta[i]->GetYaxis()->SetRangeUser(0, 120);
        hGenPtVsEta[i]->GetXaxis()->SetTitle("Gen #eta^{jet}");
        hGenPtVsEta[i]->GetYaxis()->SetTitle("Gen p_{T}^{jet} (GeV)");
        plotCMSHeader(collSystem, energy);
        gPad->SetLogz();
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("#hat{p}_{T} > %d GeV", ptHatStart + (ptHatBins[i] - 1) * ptHatStep));

        for (unsigned int j = 0; j < jetPtBins.size(); j++)
        {

            // Create canvas
            cClosureEta[i][j] = new TCanvas(Form("cClosureEta_%d_%d", i, j), Form("cClosureEta_%d_%d", i, j), 700, 800);
            // Make projection of the 3D histograms
            hRecoEta[i][j] = dynamic_cast<TH1D *>(hRecoPtEtaPtHat->ProjectionX(Form("hRecoEta_%d_%d", i, j),
                                                                               jetPtBins[j], jetPtBinsMax,
                                                                               ptHatBins[i], ptHatBinsMax));
            set1DStyle(hRecoEta[i][j], 0, kTRUE);
            hGenEta[i][j] = dynamic_cast<TH1D *>(hGenPtEtaPtHat->ProjectionX(Form("hGenEta_%d_%d", i, j),
                                                                             jetPtBins[j], jetPtBinsMax,
                                                                             ptHatBins[i], ptHatBinsMax));
            set1DStyle(hGenEta[i][j], 1, kTRUE);

            // Plot distributions
            drawJecComparison(cClosureEta[i][j], hRecoEta[i][j], hGenEta[i][j],
                           jetPtStart + (jetPtBins[j] - 1) * jetPtStep,
                           ptHatStart + (ptHatBins[i] - 1) * ptHatStep,
                           "Reco", "Gen", collSystem, energy, false, true);

        } // for (unsigned int j = 0; j < jetPtBins.size(); j++)

        for (unsigned int j = 0; j < jetEtaBinsLow.size(); j++)
        {

            // Create canvas
            cClosurePt[i][j] = new TCanvas(Form("cClosurePt_%d_%d", i, j), Form("cClosurePt_%d_%d", i, j), 700, 800);
            // Make projection of the 3D histograms
            hRecoPt[i][j] = dynamic_cast<TH1D *>(hRecoPtEtaPtHat->ProjectionY(Form("hRecoPt_%d_%d", i, j),
                                                                              jetEtaBinsLow[j], jetEtaBinsHigh[j],
                                                                              ptHatBins[i], ptHatBinsMax));
            set1DStyle(hRecoPt[i][j], 0, kFALSE);
            hGenPt[i][j] = dynamic_cast<TH1D *>(hGenPtEtaPtHat->ProjectionY(Form("hGenPt_%d_%d", i, j),
                                                                            jetEtaBinsLow[j], jetEtaBinsHigh[j],
                                                                            ptHatBins[i], ptHatBinsMax));
            set1DStyle(hGenPt[i][j], 1, kFALSE);

            // Plot distributions
            drawJecComparison(cClosurePt[i][j], hRecoPt[i][j], hGenPt[i][j],
                           jetPtStart + (jetPtBins[j] - 1) * jetPtStep,
                           ptHatStart + (ptHatBins[i] - 1) * ptHatStep,
                           "Reco", "Gen", collSystem, energy, false, true, true,
                           -5.2 + (jetEtaBinsLow[j] - 1) * etaStep, -5.2 + (jetEtaBinsHigh[j] - 1) * etaStep);

        } // for (unsigned int j = 0; j < jetEtaBinsLow.size(); j++)
    } // for (unsigned int i = 0; i < ptHatBins.size(); i++)
}

//________________
void plotSimpleInclusiveJetJECClosures(TFile *f, int collSystem = 0, double energy = 5.02) {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // Energy in TeV

    // Lab frame
    TH1D *hGenEtaLab  = dynamic_cast<TH1D*>( f->Get("hGenGoodInclusiveJetEtaLabFrame") );
    TH1D *hRecoEtaLab = dynamic_cast<TH1D*>( f->Get("hRecoGoodInclusiveJetEtaLabFrame") );

    set1DStyle( hGenEtaLab, 1, kTRUE );
    set1DStyle( hRecoEtaLab, 0, kTRUE );

    TCanvas *cSimpleJESLab = new TCanvas( "cSimpleJESLab", "cSimpleJESLab", 700, 800 );
    drawJecComparison(cSimpleJESLab, hRecoEtaLab, hGenEtaLab, 30, 15, "Reco", "Gen", collSystem, energy, false);

    // CM frame
    TH1D *hGenEtaCM  = dynamic_cast<TH1D*>( f->Get("hGenGoodInclusiveJetEtaCMFrame") );
    TH1D *hRecoEtaCM = dynamic_cast<TH1D*>( f->Get("hRecoGoodInclusiveJetEtaCMFrame") );

    set1DStyle( hGenEtaCM, 1, kTRUE );
    set1DStyle( hRecoEtaCM, 0, kTRUE );

    TCanvas *cSimpleJESCM = new TCanvas( "cSimpleJESCM", "cSimpleJESCM", 700, 800 );
    drawJecComparison(cSimpleJESCM, hRecoEtaCM, hGenEtaCM, 30, 15, "Reco", "Gen", collSystem, energy, true);
}

//________________
void calculateMeanAndErrorInRange(TH2D* h2D, TH1D *h1D) {


    double xValue = 0.0; double yValue = 0.0;
    double sum = 0.0;
    double weightedSum = 0.0;
    double sumOfWeights = 0.0;
    double sumOfWeightsSquared = 0.0;

    // Loop over x bins
    for (int i = 1; i <= h2D->GetNbinsX(); i++) {
        xValue = h2D->GetXaxis()->GetBinCenter(i);

        sum = 0.0;
        weightedSum = 0.0;

        int nonEmptyBins = 0;
        // Loop over y bins
        for (int j = 1; j <= h2D->GetNbinsY(); j++) {
            if (h2D->GetBinContent(i, j) <= 0) continue;
            nonEmptyBins++;
            yValue = h2D->GetYaxis()->GetBinCenter(j);
            weightedSum += yValue * h2D->GetBinContent(i, j);
        } // for (int j = 1; j <= h2D->GetNbinsY(); j++)

        if (nonEmptyBins != 0) {
            weightedSum /= nonEmptyBins;
        }
    }
    

    double mean = weightedSum / sum;
    double error = sqrt(sumOfWeightsSquared - pow(sumOfWeights, 2)) / sumOfWeights;

    

    h1D->SetBinContent(1, mean);
    h1D->SetBinError(1, error);
}


///________________
///
/// Plot JEC factor as a function of eta for the given pt bin ranges
///
void plotJECFactor(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250222" ) {

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    // Energy in TeV

    // Retrieve histogram
    TH3D *hJecFactor = dynamic_cast<TH3D *>(f->Get("hRecoInclusiveJetJECFactorVsPtEta"));
    if (!hJecFactor) {
        std::cerr << "Histogram 'hRecoInclusiveJetJECFactorVsPtEta' not found." << std::endl;
        return;
    }

    // Default histogram has 150 bins from 5 to 1505 GeV with 10 GeV step
    // int nPtBins = hJecFactor->GetNbinsY();
    // double ptMin = hJecFactor->GetYaxis()->GetBinLowEdge(1);
    // double ptMax = hJecFactor->GetYaxis()->GetBinUpEdge(nPtBins);
    // double ptStep = (ptMax - ptMin) / nPtBins;
    int nPtBins = 150;
    double ptMin = 5;
    double ptMax = 1505;
    double ptStep = 10;
    std::vector<int> ptBinsLow{1};
    std::vector<int> ptBinsHigh{150};
    // std::vector<int> ptBinsLow{2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 25};
    // std::vector<int> ptBinsHigh{2, 3, 4, 5, 6, 7, 8, 9, 11, 14, 24, 150};

    TH2D *hJecVsEta[ptBinsLow.size()];

    //
    // Loop over pt bins, make projection on JEC factor vs. eta, calculate
    // mean of JEC as a factor of eta and fill corresponding histograms
    //
    for (unsigned int iPt = 0; iPt < ptBinsLow.size(); ++iPt) {

        // Set pt range
        hJecFactor->GetYaxis()->SetRange( ptBinsLow[iPt], ptBinsHigh[iPt] );

        // Projection into JEC vs. eta axis
        hJecVsEta[iPt] = dynamic_cast<TH2D *>( hJecFactor->Project3D("zx") );
        hJecVsEta[iPt]->SetName( Form("hJecVsEta_%d", iPt) );
        set2DStyle( hJecVsEta[iPt] );

    } // for (unsigned int i = 0; i < ptBinsLow.size(); ++i)

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
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_80.root", uname.Data()) );
    TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_eta25.root", uname.Data()) );
    if ( !pPb8160EmbedFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_eta25.root", uname.Data()) << std::endl;
        return;
    }

    // Comparison of dijet reco and ref to gen distributions
    comparisons2gen( pPb8160EmbedFile, collisionSystem, collisionEnergy, date );

    // Plot simple inclusicve jet JEC closure (inclusive jets within |eta|<1.4)
    // plotSimpleInclusiveJetJECClosures(pPb8160EmbedFile, collisionSystem, collisionEnergy);
    
    // Plot for inclusive jets JEC closures (scan in eta and pT)
    // plotInclusiveJetJECClosures(pPb8160EmbedFile, collisionSystem, collisionEnergy);
}
