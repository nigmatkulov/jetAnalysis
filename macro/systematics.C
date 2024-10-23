// ROOT headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TSystemDirectory.h"
#include <sys/stat.h>
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraph.h"

// C++ headers
#include <fstream>
#include <iostream>
#include <vector>
                               // 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200, 250, 300, 500,  extra: 95 -- 115
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

const Int_t dijetEtaFBBins{15};
Double_t dijetEtaFBVals[dijetEtaFBBins+1] { 0.0,  0.2,  0.4,  0.6,  0.8,  
                                            1.0,  1.2,  1.4,  1.6,  1.8,  
                                            2.0,  2.2,  2.4,  3.0,  4.0,  
                                            5.0 };

Int_t nPads{4};

// nPDF pT bins
std::vector<Int_t> npdfPtLow = {50, 60, 100, 140, 300}; 
std::vector<Int_t> npdfPtHi =  {60, 70, 110, 150, 500};

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
bool directoryExists(const char* directoryPath) {
    TSystemFile file(directoryPath, "");
    return file.IsDirectory();
}

//________________
bool checkFileIsGood(TFile *inFile) {
    bool isGood{true};
    std::cout << Form("Check that file %s is good", inFile->GetName());
    if ( !inFile ) {
        std::cout << "[ERROR] Input file does not exist\n";
        isGood = false;
    }
    if ( inFile->IsZombie() ) {
        std::cout << "[ERROR] Beware of zombies!!!!\n";
        isGood = false;
    }
    std::cout << "\t[DONE]" << std::endl;
    return isGood;
}

//________________
void plotCMSHeader() {
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.51, 0.93, "pPb 174.6 nb^{-1} (8.16 TeV)");
    t.SetTextSize(0.05);
}

//________________
void createDirectory(const char* directoryPath) {
    // Create the directory with read, write, and execute permissions for owner, group, and others
    if (mkdir(directoryPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
        std::cout << "Directory created successfully." << std::endl;
    } 
    else {
        std::cerr << "Failed to create directory." << std::endl;
    }
    TString jeuPath = directoryPath;
    TString jeuOldPath = directoryPath;
    TString jerPath = directoryPath;
    TString pileupPath = directoryPath;
    TString pointingPath = directoryPath;
    TString dataPath = directoryPath;
    TString data2mcPath = directoryPath;
    TString systematicsPath = directoryPath;
    jeuPath += "/jes";
    jeuOldPath += "/jeu";
    jerPath += "/jer";
    pileupPath += "/pileup";
    pointingPath += "/pointing";
    dataPath += "/data";
    data2mcPath += "/data2mc";
    systematicsPath += "/systematics";

    if ( !directoryExists( jeuPath.Data() ) ) {
        if (mkdir(jeuPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "jeu directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( jeuOldPath.Data() ) ) {
        if (mkdir(jeuOldPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "jeu directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( jerPath.Data() ) ) {
        if (mkdir(jerPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "jer directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( pileupPath.Data() ) ) {
        if (mkdir(pileupPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "pileup directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( pointingPath.Data() ) ) {
        if (mkdir(pointingPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "pointing resolution directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( dataPath.Data() ) ) {
        if (mkdir(dataPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "data directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( data2mcPath.Data() ) ) {
        if (mkdir(data2mcPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "data2mc directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }

    if ( !directoryExists( systematicsPath.Data() ) ) {
        if (mkdir(systematicsPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
            std::cout << "systematics directory created successfully." << std::endl;
        } 
        else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    }
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
void rescaleEta(TH2* h) {
    for (Int_t iBin=1; iBin<=h->GetNbinsX(); iBin++) {
        for (Int_t jBin=1; jBin<=h->GetNbinsY(); jBin++) {
            Double_t val = h->GetBinContent( iBin, jBin );
            Double_t valErr = h->GetBinError( iBin, jBin );
            Double_t binWidthX = h->GetXaxis()->GetBinWidth( iBin );
            Double_t binWidthY = h->GetYaxis()->GetBinWidth( jBin );
            h->SetBinContent( iBin, jBin, val / (binWidthX * binWidthY) );
            h->SetBinError( iBin, jBin, valErr / (binWidthX * binWidthY) );
        } // for (Int_t jBin=1; jBin<=h->GetNbinsY(); jBin++)
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
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./h->Integral() );
    }
}

//________________
void setSystUncrtStyle(TH1* h, Int_t type = 0) {
    Int_t color = 29;        // Green-like
    Int_t lineWidth = 0;
    Int_t markerStyle = 20;
    if ( type == 0 ) {
        color = 29;          // Green-like
        markerStyle = 1;
    }
    else if ( type == 1) {
        color = 2;          // Red
        lineWidth = 2;
        markerStyle = 20;
    }
    else if ( type == 2 ) {
        color = 4;           // Blue
        lineWidth = 2;
        markerStyle = 24;
    }
    else {
        color = 46;          // Dark red-like
        lineWidth = 2;
        markerStyle = 22;
    }
    h->SetFillColorAlpha(color, 0.35);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(lineWidth);
    h->SetFillStyle(1001);           // Solid fill

    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);
}

//________________
void setGraphStyle(TGraph* gr, Bool_t isCM = "kFALSE") {
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetNdivisions(208);
    gr->GetYaxis()->SetNdivisions(208);    
    gr->GetYaxis()->SetTitleOffset(1.1);
    if ( isCM ) {
        gr->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
        gr->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}_{CM}");
    }
    else {
        gr->GetXaxis()->SetTitle("#eta^{dijet}");
        gr->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    }
}

//________________
void set2DStyle(TH2* h, Bool_t doRenorm = kFALSE) {
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.1);

    if ( doRenorm ) {
        h->Scale( 1./ h->Integral() );
    }
}

//________________
void makeProjectionsFrom2D(TH2D *h2D, TH1D *hProjX[], TH1D *hProjY[], 
                           Int_t nBinsX = 1, Int_t nBinsY = 1, 
                           const Char_t *hNameX = "hProjX", const Char_t *hNameY = "hProjY",
                           Int_t style = 1) {

    // Make projections on X axis
    for (Int_t i{0}; i<nBinsY; i++) {
        hProjX[i] = dynamic_cast<TH1D*>( h2D->ProjectionX( Form("%s_%d", hNameX, i), i+1, i+1) );
        hProjX[i]->SetNameTitle( Form("%s_%d", hNameX, i), ";#eta");
        set1DStyle( hProjX[i], style );
        hProjX[i]->SetMarkerSize(0.7);
    } // for (Int_t i{1}; i<=nBinsX; i++)

    // Make projections on Y axis
    for (Int_t i{0}; i<nBinsX; i++) {
        hProjY[i] = dynamic_cast<TH1D*>( h2D->ProjectionY( Form("%s_%d", hNameY, i), i+1, i+1) );
        hProjY[i]->SetNameTitle( Form("%s_%d", hNameY, i), ";p_{T} (GeV)");
        set1DStyle( hProjY[i], style );
        hProjY[i]->SetMarkerSize(0.7);
    } // for (Int_t i{1}; i<=nBinsY; i++)
}

//________________
void make1DRatio(TH1D *hRat, TH1D *hDen, const Char_t *ratioName = "Ratio to default", Int_t style = 0, 
                 Bool_t isBinomial = kFALSE) {

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
TH1D* projectEtaFrom3D(TH3D* h3D, const Char_t *name, 
                       Int_t ptDijetBinLow, Int_t ptDijetBinHi) {

    TH1D *tmp = dynamic_cast<TH1D*>( h3D->ProjectionY( name, ptDijetBinLow, ptDijetBinHi ) );
    
    tmp->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    return tmp;
}

//________________
TH1D* projectEtaFrom2D(TH2D* h2D, const Char_t *name, 
                       Int_t ptDijetBinLow, Int_t ptDijetBinHi) {

    TH1D *tmp = dynamic_cast<TH1D*>( h2D->ProjectionY( name, ptDijetBinLow, ptDijetBinHi ) );
    tmp->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    return tmp;
}

//________________
void dirSystematics(TF1 *fitf, TH1D *h, TH1D* hSyst, TGraph* syst) {

    std::cout << Form("Def file name: %s\n", h->GetName());
    // Loop over the bins and estimate systematics for those that have entries
    for (Int_t i{1}; i<=h->GetNbinsX(); i++ ) {
        //if ( h->GetBinContent(i) == 0 ) continue;
        Double_t xVal = h->GetXaxis()->GetBinCenter(i);
        Double_t yVal = TMath::Abs( fitf->Eval( xVal ) - 1.);
        if ( TMath::Abs( h->GetBinContent(i) ) < 0.0001 ) {
            yVal = 0.;
        }
        std::cout << Form("bin content: %f x: %3.2f y [perc.]: %.3f ", h->GetBinContent(i), xVal, yVal * 100.) << std::endl;;
        syst->SetPoint(i-1, xVal, yVal );

        hSyst->SetBinContent(i, yVal);
    } // for (Int_t i{1}; i<=h->GetNbinsX(); i++ )
}

//________________
void compareData2McDifferentDirections(TFile *expPbGoing, TFile *expPGoing, 
                                       TFile *mcPbGoing, TFile *mcPGoing, TString date,
                                       TFile *defFile) {

    TH2D *hExpPtEtaPbGoing = dynamic_cast<TH2D*>( expPbGoing->Get("hRecoDijetPtEta") );
    hExpPtEtaPbGoing->SetName("hExpPtEtaPbGoing");
    TH2D *hExpPtEtaPGoing = dynamic_cast<TH2D*>( expPGoing->Get("hRecoDijetPtEta") );
    hExpPtEtaPGoing->SetName("hExpPtEtaPGoing");

    TH2D *hMcPtEtaPbGoing = dynamic_cast<TH2D*>( mcPbGoing->Get("hRecoDijetPtEta") );
    hMcPtEtaPbGoing->SetName("hMcPtEtaPbGoing");
    TH2D *hMcPtEtaPGoing = dynamic_cast<TH2D*>( mcPGoing->Get("hRecoDijetPtEta") );
    hMcPtEtaPGoing->SetName("hMcPtEtaPGoing");

    TH3D *hPtEtaDphiDef = (TH3D*)defFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiDef->SetName("hPtEtaDphiDef");

    TString frame = "lab";

    // Do rebinning if needed
    Int_t rebinX{1}, rebinY{2};
    hExpPtEtaPbGoing->RebinX( rebinX );
    hExpPtEtaPbGoing->RebinY( rebinY );
    hExpPtEtaPGoing->RebinX( rebinX );
    hExpPtEtaPGoing->RebinY( rebinY );
    hMcPtEtaPbGoing->RebinX( rebinX );
    hMcPtEtaPbGoing->RebinY( rebinY );
    hMcPtEtaPGoing->RebinX( rebinX );
    hMcPtEtaPGoing->RebinY( rebinY );
    // hPtEtaDphiDef->RebinX( rebinX );
    // hPtEtaDphiDef->RebinY( rebinY );

    // Double_t dijetPtVals[dijetPtBins+1] {  40.,  50.,   60.,  70.,  80.,
    //                                        90., 100.,  110., 120., 130.,
    //                                       140., 150.,  160., 180., 200., 
    //                                       240., 300., 1000.};
    // fDijetPtBins{194}, fDijetPtRange{30., 1000.}

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};

    // Styles
    Int_t expPbGoingType{0};
    Int_t expPGoingType{3};
    Int_t mcPbGoingType{4};
    Int_t mcPGoingType{1};
    Int_t expOverMcType{2};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hExpEtaPbGoing[ ptDijetBinLow.size() ];
    TH1D *hExpEtaPGoing[ ptDijetBinLow.size() ];
    TH1D *hExpEtaPb2PGoingRatio[ ptDijetBinLow.size() ];
    TH1D *hMcEtaPbGoing[ ptDijetBinLow.size() ];
    TH1D *hMcEtaPGoing[ ptDijetBinLow.size() ];
    TH1D *hMcEtaPb2PGoingRatio[ ptDijetBinLow.size() ];
    TH1D *hExpOverMcPb2PGoingRatio[ ptDijetBinLow.size() ];
    TH1D *hDefEta[ ptDijetBinLow.size() ];

    TF1 *fitRatio[ ptDijetBinLow.size() ];
    TGraph *grSyst[ ptDijetBinLow.size() ];

    TH1D *hDirSyst[ ptDijetBinLow.size() ];
    TH1D *hJeuSyst[ ptDijetBinLow.size() ];
    TH1D *hJerSyst[ ptDijetBinLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetBinLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        // Exp Pb-going projection
        hExpEtaPbGoing[i] = dynamic_cast<TH1D*>( hExpPtEtaPbGoing->ProjectionY( Form("hExpEtaPbGoing_%d", i), 
                                                                                ptDijetBinLow.at(i), ptDijetBinHi.at(i) ) );
        rescaleEta( hExpEtaPbGoing[i] );
        set1DStyle(hExpEtaPbGoing[i], expPbGoingType);

        // Exp p-going projection
        hExpEtaPGoing[i] = dynamic_cast<TH1D*>( hExpPtEtaPGoing->ProjectionY( Form("hExpEtaPGoing_%d", i), 
                                                                              ptDijetBinLow.at(i), ptDijetBinHi.at(i) ));
        rescaleEta( hExpEtaPGoing[i] );
        set1DStyle(hExpEtaPGoing[i], expPGoingType);

        // Ratio of Pb-going to p-going
        hExpEtaPb2PGoingRatio[i] = dynamic_cast<TH1D*>( hExpEtaPbGoing[i]->Clone( Form("hExpEtaPb2PGoingRatio_%d", i) ) );
        make1DRatio(hExpEtaPb2PGoingRatio[i], hExpEtaPGoing[i], Form("hExpEtaPb2PGoingRatio_%d", i) );
        hExpEtaPb2PGoingRatio[i]->GetYaxis()->SetTitle("Pb-going / p-going");
        set1DStyle( hExpEtaPb2PGoingRatio[i], expPbGoingType);

        // Mc Pb-going projection
        hMcEtaPbGoing[i] = dynamic_cast<TH1D*>( hMcPtEtaPbGoing->ProjectionY( Form("hMcEtaPbGoing_%d", i), 
                                                                                ptDijetBinLow.at(i), ptDijetBinHi.at(i) ) );
        rescaleEta( hMcEtaPbGoing[i] );
        set1DStyle(hMcEtaPbGoing[i], mcPbGoingType);

        // Mc p-going projection
        hMcEtaPGoing[i] = dynamic_cast<TH1D*>( hMcPtEtaPGoing->ProjectionY( Form("hMcEtaPGoing_%d", i), 
                                                                              ptDijetBinLow.at(i), ptDijetBinHi.at(i) ));
        rescaleEta( hMcEtaPGoing[i] );
        set1DStyle(hMcEtaPGoing[i], mcPGoingType);

        // Default dijet eta distribution
        hDefEta[i] = projectEtaFrom3D(hPtEtaDphiDef, Form("hEtaDef_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hDefEta[i] ); 
        set1DStyle( hDefEta[i], expOverMcType );

        // Ratio of Pb-going to p-going
        hMcEtaPb2PGoingRatio[i] = dynamic_cast<TH1D*>( hMcEtaPbGoing[i]->Clone( Form("hMcEtaPb2PGoingRatio_%d", i) ) );
        make1DRatio(hMcEtaPb2PGoingRatio[i], hMcEtaPGoing[i], Form("hMcEtaPb2PGoingRatio_%d", i) );
        hMcEtaPb2PGoingRatio[i]->GetYaxis()->SetTitle("Pb-going / p-going");
        set1DStyle( hMcEtaPb2PGoingRatio[i], mcPbGoingType);

        // Make double ratio
        hExpOverMcPb2PGoingRatio[i] = dynamic_cast<TH1D*>( hExpEtaPb2PGoingRatio[i]->Clone( Form("hExpOverMcPb2PGoingRatio_%d", i) ) );
        make1DRatio( hExpOverMcPb2PGoingRatio[i], hMcEtaPb2PGoingRatio[i], Form("hExpOverMcPb2PGoingRatio_%d", i) );
        hExpOverMcPb2PGoingRatio[i]->GetYaxis()->SetTitle("Data / MC (Pb-going / p-going)");
        set1DStyle( hExpOverMcPb2PGoingRatio[i], expOverMcType );

        // Fit function
        // fitRatio[i] = new TF1(Form("fitRatio_%d",i), "[0] + [1] * x + [2] * x * x", -3.5, 3.5);
        fitRatio[i] = new TF1(Form("fitRatio_%d",i), "[0] + [1] * (x - [2]) * (x - [2])", -3.5, 3.5);
        fitRatio[i]->SetParameters(0., 0.001, 0.465);
        // fitRatio[i]->SetParLimits(2, 0., 1e6);
        fitRatio[i]->SetParLimits(1, 0., 1e6);
        fitRatio[i]->SetLineColor( kRed );
        fitRatio[i]->SetLineWidth( 2 );
        hExpOverMcPb2PGoingRatio[i]->Fit(Form("fitRatio_%d",i), "MRE0");
        hDirSyst[i] = dynamic_cast<TH1D*>( hDefEta[i]->Clone( Form("hDirSyst_%d",i) ) );
        grSyst[i] = new TGraph();
        dirSystematics(fitRatio[i], hDefEta[i], hDirSyst[i], grSyst[i]);
        grSyst[i]->SetName( Form("dirSyst_%d",i) );

        // Write systematic uncertainty to ASCII file
        std::ofstream outFile( Form("%s/pPb8160_etaDijet_dirSyst_%d_%d_%s.txt", date.Data(), 
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFile.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSyst[i]->GetN(); iPoint++) {
                std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSyst[i]->GetPoint(iPoint, x, y);
                std::cout << " x: " << x << " y: " << y << std::endl;
                outFile << x << " " << y << std::endl;
            }
            outFile.close();
        }
        else {
            std::cerr << "[DirSyst] Unable to open file for writing" << std::endl;
        }


        // Plot comparison
        canv->cd();
        setPadStyle();
        hExpEtaPb2PGoingRatio[i]->Draw();
        hMcEtaPb2PGoingRatio[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hExpEtaPb2PGoingRatio[i], Form("Data"), "p");
        leg->AddEntry(hMcEtaPb2PGoingRatio[i], Form("MC"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_exp2mc_dirComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios of Pb-going to p-going
        canv->cd();
        setPadStyle();
        hExpOverMcPb2PGoingRatio[i]->Draw();
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hExpOverMcPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_ exp2mc_dirRat_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hExpEtaPb2PGoingRatio[i]->Draw();
        hMcEtaPb2PGoingRatio[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hExpEtaPb2PGoingRatio[i], Form("Data"), "p");
        leg->AddEntry(hMcEtaPb2PGoingRatio[i], Form("MC"), "p");
        leg->Draw();

        // Plot all ratios on one canvas
        cRat->cd(i+1);
        setPadStyle();
        hExpOverMcPb2PGoingRatio[i]->Draw();
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hExpOverMcPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        fitRatio[i]->Draw("same");
    } // for (Int_t i=0; i<ptDijetBinLow.size(); i++)

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_exp2mc_dirComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_exp2mc_dirRat_all_%s.pdf", date.Data(), frame.Data() ) );

    TFile *oFile = new TFile( Form("%s/pPb8160_dirSyst_%s.root", date.Data(), frame.Data() ), "recreate" );
    for (Int_t i{0}; i<ptBins; i++) {
        hDirSyst[i]->GetYaxis()->SetRangeUser(0., 0.3);
        hDirSyst[i]->GetXaxis()->SetRangeUser(-4.8, 4.8);
        hDirSyst[i]->Write();
    } // for (Int_t i{0}; i<ptBins; i++)
    oFile->Close();
}

//________________
void makeDijetFB(TH1D* hOrig, TH1D* hBack, TH1D* hForw, TH1D* hBFRatio, TH1D* hFBRatio) {
    Int_t nBins = hOrig->GetNbinsX();
    Int_t middle = nBins / 2;
    for (Int_t i{1}; i <= middle; i++) {
        Double_t xBack = TMath::Abs( hOrig->GetBinCenter(middle-i+1) );
        Int_t binBack = hBack->FindBin(xBack);
        Double_t binBackVal = hBack->GetBinContent(binBack);
        binBackVal += hOrig->GetBinContent(middle-i+1);
        Double_t binBackErr = hBack->GetBinError(binBack);
        binBackErr = TMath::Sqrt( binBackErr * binBackErr + hOrig->GetBinError(middle-i+1) * hOrig->GetBinError(middle-i+1) );
        hBack->SetBinContent(binBack,  binBackVal );
        hBack->SetBinError(binBack, binBackErr );

        hBFRatio->SetBinContent(binBack,  binBackVal );
        hBFRatio->SetBinError(binBack, binBackErr );

        Double_t xForw = TMath::Abs( hOrig->GetBinCenter(middle+i) );
        Int_t binForw = hForw->FindBin(xForw);
        Double_t binForwVal = hForw->GetBinContent(binForw);
        binForwVal += hOrig->GetBinContent(middle+i);
        Double_t binForwErr = hForw->GetBinError(binForw);
        binForwErr = TMath::Sqrt( binForwErr * binForwErr + hOrig->GetBinError(middle+i) * hOrig->GetBinError(middle+i) );
        hForw->SetBinContent(binForw, binForwVal );
        hForw->SetBinError(binForw, binForwErr );

        hFBRatio->SetBinContent(binForw, binForwVal );
        hFBRatio->SetBinError(binForw, binForwErr );

        // std::cout << "i: " << i << " xBackOrig: " << hOrig->GetBinCenter(middle-i+1) << " yBackOrg: " << hOrig->GetBinContent(middle-i+1) << std::endl;
        // std::cout << "i: " << i << " xForwOrig: " << hOrig->GetBinCenter(middle+i) << " yForwOrg: " << hOrig->GetBinContent(middle+i) << std::endl;
        // std::cout << "i: " << i << " xBack: " << hBack->GetBinCenter(i) << " yBack: " << hBack->GetBinContent(i) << std::endl;
        // std::cout << "i: " << i << " xForw: " << hForw->GetBinCenter(i) << " yForw: " << hForw->GetBinContent(i) << std::endl;
        // std::cout << "i: " << i << " x: " << hBack->GetBinCenter(i) << " B/F: " << hBack->GetBinContent(i) / hForw->GetBinContent(i) << " F/B: " << hForw->GetBinContent(i) / hBack->GetBinContent(i) << std::endl;

    } // for (Int_t i{1}; i <= middle; i++)

    hBFRatio->Divide(hForw);
    hFBRatio->Divide(hBack);
}

//________________
std::vector<TH1D*> createFBRatios(std::vector<TH1D*> hForw, std::vector<TH1D*> hBack, const char* hName, 
                                  Int_t type = 2, Bool_t isCM = kTRUE, Bool_t isForward = kTRUE) {

    TString name = hName;
    if ( isForward ) {
        name += "_fb";
    }
    else {
        name += "_bf";
    }

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // std::cout << Form("[Create FB] hName: %s", name.Data() ) << std::endl;
    
    std::vector<TH1D*> hFBRatio;
    for (Int_t i{0}; i<hBack.size(); i++) {
        if ( isForward ) {
            hFBRatio.push_back( dynamic_cast<TH1D*>( hForw[i]->Clone( Form("%s_%d", name.Data(), i) ) ) );
            hFBRatio[i]->SetName( Form("%s_pt_%d_%d", name.Data(), ptDijetLow.at(i), ptDijetHi.at(i)) );
            hFBRatio[i]->Divide( hBack[i] );
            set1DStyle(hFBRatio[i], type);
            if ( isCM ) {
                hFBRatio[i]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
                hFBRatio[i]->GetYaxis()->SetTitle("#frac{ dN/d#eta^{dijet}_{CM} (forward) }{ dN/d#eta^{dijet}_{CM} (backward) }");
            }
            else {
                hFBRatio[i]->GetXaxis()->SetTitle("#eta^{dijet}");
                hFBRatio[i]->GetYaxis()->SetTitle("#frac{ dN/d#eta^{dijet} (forward) }{ dN/d#eta^{dijet} (backward) }");
            }
        }
        else {
            hFBRatio.push_back( dynamic_cast<TH1D*>( hBack[i]->Clone( Form("%s_%d", hName, i) ) ) );
            hFBRatio[i]->SetName( Form("%s_pt_%d_%d", name.Data(), ptDijetLow.at(i), ptDijetHi.at(i)) );
            hFBRatio[i]->Divide( hForw[i] );
            if ( isCM ) {
                hFBRatio[i]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
                hFBRatio[i]->GetYaxis()->SetTitle("#frac{ dN/d#eta^{dijet}_{CM} (backward) }{ dN/d#eta^{dijet}_{CM} (forward) }");
            }
            else {
                hFBRatio[i]->GetXaxis()->SetTitle("#eta^{dijet}");
                hFBRatio[i]->GetYaxis()->SetTitle("#frac{ dN/d#eta^{dijet} (backward) }{ dN/d#eta^{dijet} (forward) }");
            }
        }
    } // for (Int_t i{0}; i<hBack.size(); i++)
    return hFBRatio;
}

//________________
std::vector< std::vector<TH1D*> > createDijetDirectionEtaHistogramsFromLab(TFile *inFile, const char *inputHistoName,
                                                                           const char* hName, int type = 2) {

    // std::cout << Form("[Dir lab] Direction histogram input histo name : %s", inputHistoName) << std::endl;
    // std::cout << Form("[Dir lab] Direction histogram output histo name: %s", hName) << std::endl;

    // Retreive the histogram
    TH3D *hPtEtaDphi = dynamic_cast<TH3D*>( inFile->Get( inputHistoName ) );
    // Reset histogram name
    hPtEtaDphi->SetName( hName );

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Create vector of histograms to return
    std::vector< std::vector<TH1D*> > hEtaDir;
    std::vector<TH1D*> hEta;
    std::vector<TH1D*> hEtaForward;
    std::vector<TH1D*> hEtaBackward;

    Int_t upType{0};
    Int_t downType{1};

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {
        
        // Make projection on pseudorapidity axis
        hEta.push_back( projectEtaFrom3D(hPtEtaDphi, Form("%s_%d", hName, i), ptDijetBinLow.at(i), ptDijetBinHi.at(i)) );
        hEta[i]->SetName( Form("%s_pt_%d_%d", hName, ptDijetLow.at(i), ptDijetHi.at(i)) );
        rescaleEta(hEta[i]);
        hEtaForward.push_back( new TH1D( Form("%s_forward_%d", hName, i), Form("%s_forward_%d", hName, i), 50, 0., 5.) );
        hEtaForward[i]->SetName( Form("%s_forward_pt_%d_%d", hName, ptDijetLow.at(i), ptDijetHi.at(i)) );
        hEtaBackward.push_back( new TH1D( Form("%s_backward_pt_%d", hName, i), Form("%s_backward_%d", hName, i), 50, 0., 5.) );
        hEtaBackward[i]->SetName( Form("%s_backward_pt_%d_%d", hName, ptDijetLow.at(i), ptDijetHi.at(i)) );

        hEtaForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaForward[i]->Sumw2(); set1DStyle(hEtaForward[i], upType);
        hEtaBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaBackward[i]->Sumw2(); set1DStyle(hEtaBackward[i], downType);

        Int_t nBins = hEta[i]->GetNbinsX();
        Int_t middle = nBins / 2;
        for (Int_t j{1}; j <= middle; j++) {
            hEtaBackward[i]->SetBinContent(j,  hEta[i]->GetBinContent(middle-j+1) );
            hEtaBackward[i]->SetBinError(j, hEta[i]->GetBinError(middle-j+1) );
            hEtaForward[i]->SetBinContent(j, hEta[i]->GetBinContent(middle+j) );
            hEtaForward[i]->SetBinError(j, hEta[i]->GetBinError(middle+j) );
        } // for (Int_t j{1}; j <= middle; j++)
    } // for (Int_t i{0}; i<ptDijetBinLow.size(); i++)
    hEtaDir.push_back( hEtaForward );
    hEtaDir.push_back( hEtaBackward );
    return hEtaDir;
}

//________________
std::vector<TH1D*> createDijetDirectionEtaHistograms(TFile *inFile, const char *inputHistoName,
                                                     const char* hName, int type = 2, 
                                                     Bool_t isCM = kFALSE, Bool_t isForward = kTRUE) {

    TString name = hName;
    if ( isForward ) {
        name += "_forward";
    }
    else {
        name += "_backward";
    }
    // std::cout << Form("[Dir] Direction histogram input histo name : %s", inputHistoName) << std::endl;
    // std::cout << Form("[Dir] Direction histogram output histo name: %s", name.Data()) << std::endl;

    // Retreive the histogram
    TH2D* hPtEta = dynamic_cast<TH2D*>( inFile->Get( inputHistoName ) );
    // Reset histogram name
    hPtEta->SetName( hName );

    Int_t upType{0};
    Int_t downType{1};

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Create vector of histograms to return
    std::vector<TH1D*> hEtaDir;

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptDijetBinLow.size(); i++) {
            
        // Make projection on pseudorapidity axis
        hEtaDir.push_back( projectEtaFrom2D(hPtEta, Form("%s_%d", name.Data(), i), ptDijetBinLow.at(i), ptDijetBinHi.at(i)) );
        hEtaDir[i]->SetName( Form("%s_pt_%d_%d", name.Data(), ptDijetLow.at(i), ptDijetHi.at(i)) );
        rescaleEta(hEtaDir[i]); set1DStyle(hEtaDir[i], type);
        if ( isCM ) {
            hEtaDir[i]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
        }
        else {
            hEtaDir[i]->GetXaxis()->SetTitle("#eta^{dijet}");
        }
        
        if ( isForward ) {
            if ( isCM ) {
                hEtaDir[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet}_{CM} (forward)");
            }
            else {
                hEtaDir[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet} (forward)");
            }
        }
        else {
            if ( isCM ) {
                hEtaDir[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet}_{CM} (backward)");
            }
            else {
                hEtaDir[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet} (backward)");
            }
        }
    } // for (Int_t i{0}; i<ptDijetBinLow.size(); i++)

    return hEtaDir;
}

//________________
std::vector<TH1D*> createDijetEtaFullHistograms(TFile *inFile, const char *inputHistoName,
                                                const char* hName, int type = 2, 
                                                Bool_t isCM = kFALSE) {

    // Retreive the histogram
    // std::cout << Form("[Full] Input histo name: %s", inputHistoName) << std::endl;
    // std::cout << Form("[Full] Output histo name: %s", hName) << std::endl;
    TH3D* hPtEtaDphi = dynamic_cast<TH3D*>( inFile->Get( inputHistoName ) );
    // Reset histogram name
    hPtEtaDphi->SetName( hName );

    // std::cout << Form("\t[Full] hPtEtaDphi entries: %f", hPtEtaDphi->Integral() ) << std::endl;

    Int_t upType{0};
    Int_t downType{1};

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Create vector of histograms to return
    std::vector<TH1D*> hEta;

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptDijetBinLow.size(); i++) {
            
        // Make projection on pseudorapidity axis
        hEta.push_back( projectEtaFrom3D(hPtEtaDphi, Form("%s_%d", hName, i), ptDijetBinLow.at(i), ptDijetBinHi.at(i)) );
        hEta[i]->SetName( Form("%s_pt_%d_%d", hName, ptDijetLow.at(i), ptDijetHi.at(i)) );
        rescaleEta(hEta[i]); set1DStyle(hEta[i], 2);
        if ( isCM ) {
            hEta[i]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hEta[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet}_{CM}");
        }
        else {
            hEta[i]->GetXaxis()->SetTitle("#eta^{dijet}");
            hEta[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet}");
        }

        //std::cout << Form("\t[Full] %s integral: %f entries: %f", hEta[i]->GetName(), hEta[i]->Integral(), hEta[i]->GetEntries() ) << std::endl;
    } // for (Int_t i{0}; i<ptDijetBinLow.size(); i++)

    return hEta;
}

//________________
std::vector< std::vector<TH1D*> > createRatios2MC(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist,
                                                  std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist,
                                                  std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hGenEtaDist,
                                                  Bool_t isCM = kFALSE) {

    TString frame;
    if ( isCM ) {
        frame = "cms";
    }
    else {
        frame = "lab";
    }

    // Create vector of histograms to return
    std::vector< std::vector<TH1D*> > hRatios;
    
    std::vector<TH1D*> hMBEta = hMBEtaDist[0];
    std::vector<TH1D*> hMBEtaForward = hMBEtaDist[1];
    std::vector<TH1D*> hMBEtaBackward = hMBEtaDist[2];
    std::vector<TH1D*> hMBEtaFBRatio = hMBEtaDist[3];
    std::vector<TH1D*> hMBEtaBFRatio = hMBEtaDist[4];

    std::vector<TH1D*> hJet60Eta = hJet60EtaDist[0];
    std::vector<TH1D*> hJet60EtaForward = hJet60EtaDist[1];
    std::vector<TH1D*> hJet60EtaBackward = hJet60EtaDist[2];
    std::vector<TH1D*> hJet60EtaFBRatio = hJet60EtaDist[3];
    std::vector<TH1D*> hJet60EtaBFRatio = hJet60EtaDist[4];

    std::vector<TH1D*> hJet80Eta = hJet80EtaDist[0];
    std::vector<TH1D*> hJet80EtaForward = hJet80EtaDist[1];
    std::vector<TH1D*> hJet80EtaBackward = hJet80EtaDist[2];
    std::vector<TH1D*> hJet80EtaFBRatio = hJet80EtaDist[3];
    std::vector<TH1D*> hJet80EtaBFRatio = hJet80EtaDist[4];

    std::vector<TH1D*> hJet100Eta = hJet100EtaDist[0];
    std::vector<TH1D*> hJet100EtaForward = hJet100EtaDist[1];
    std::vector<TH1D*> hJet100EtaBackward = hJet100EtaDist[2];
    std::vector<TH1D*> hJet100EtaFBRatio = hJet100EtaDist[3];
    std::vector<TH1D*> hJet100EtaBFRatio = hJet100EtaDist[4];

    std::vector<TH1D*> hEmbeddingEta = hEmbeddingEtaDist[0];
    std::vector<TH1D*> hEmbeddingEtaForward = hEmbeddingEtaDist[1];
    std::vector<TH1D*> hEmbeddingEtaBackward = hEmbeddingEtaDist[2];
    std::vector<TH1D*> hEmbeddingEtaFBRatio = hEmbeddingEtaDist[3];
    std::vector<TH1D*> hEmbeddingEtaBFRatio = hEmbeddingEtaDist[4];

    std::vector<TH1D*> hGenEta = hGenEtaDist[0];
    std::vector<TH1D*> hGenEtaForward = hGenEtaDist[1];
    std::vector<TH1D*> hGenEtaBackward = hGenEtaDist[2];
    std::vector<TH1D*> hGenEtaFBRatio = hGenEtaDist[3];
    std::vector<TH1D*> hGenEtaBFRatio = hGenEtaDist[4];

    // Create vectors with ratios
    std::vector<TH1D*> hMBEta2GenEtaRatio;
    std::vector<TH1D*> hMBEtaForward2GenEtaForwardRatio;
    std::vector<TH1D*> hMBEtaBackward2GenEtaBackwardRatio;
    std::vector<TH1D*> hMBEtaFBRatio2GenEtaFBRatio;
    std::vector<TH1D*> hMBEtaBFRatio2GenEtaBFRatio;

    std::vector<TH1D*> hJet60Eta2GenEtaRatio;
    std::vector<TH1D*> hJet60EtaForward2GenEtaForwardRatio;
    std::vector<TH1D*> hJet60EtaBackward2GenEtaBackwardRatio;
    std::vector<TH1D*> hJet60EtaFBRatio2GenEtaFBRatio;
    std::vector<TH1D*> hJet60EtaBFRatio2GenEtaBFRatio;

    std::vector<TH1D*> hJet80Eta2GenEtaRatio;
    std::vector<TH1D*> hJet80EtaForward2GenEtaForwardRatio;
    std::vector<TH1D*> hJet80EtaBackward2GenEtaBackwardRatio;
    std::vector<TH1D*> hJet80EtaFBRatio2GenEtaFBRatio;
    std::vector<TH1D*> hJet80EtaBFRatio2GenEtaBFRatio;

    std::vector<TH1D*> hJet100Eta2GenEtaRatio;
    std::vector<TH1D*> hJet100EtaForward2GenEtaForwardRatio;
    std::vector<TH1D*> hJet100EtaBackward2GenEtaBackwardRatio;
    std::vector<TH1D*> hJet100EtaFBRatio2GenEtaFBRatio;
    std::vector<TH1D*> hJet100EtaBFRatio2GenEtaBFRatio;

    std::vector<TH1D*> hEmbeddingEta2GenEtaRatio;
    std::vector<TH1D*> hEmbeddingEtaForward2GenEtaForwardRatio;
    std::vector<TH1D*> hEmbeddingEtaBackward2GenEtaBackwardRatio;
    std::vector<TH1D*> hEmbeddingEtaFBRatio2GenEtaFBRatio;
    std::vector<TH1D*> hEmbeddingEtaBFRatio2GenEtaBFRatio;

    std::vector<TH1D*> hMBEta2EmbeddingEtaRatio;
    std::vector<TH1D*> hMBEtaForward2EmbeddingEtaForwardRatio;
    std::vector<TH1D*> hMBEtaBackward2EmbeddingEtaBackwardRatio;
    std::vector<TH1D*> hMBEtaFBRatio2EmbeddingEtaFBRatio;
    std::vector<TH1D*> hMBEtaBFRatio2EmbeddingEtaBFRatio;

    std::vector<TH1D*> hJet60Eta2EmbeddingEtaRatio;
    std::vector<TH1D*> hJet60EtaForward2EmbeddingEtaForwardRatio;
    std::vector<TH1D*> hJet60EtaBackward2EmbeddingEtaBackwardRatio;
    std::vector<TH1D*> hJet60EtaFBRatio2EmbeddingEtaFBRatio;
    std::vector<TH1D*> hJet60EtaBFRatio2EmbeddingEtaBFRatio;

    std::vector<TH1D*> hJet80Eta2EmbeddingEtaRatio;
    std::vector<TH1D*> hJet80EtaForward2EmbeddingEtaForwardRatio;
    std::vector<TH1D*> hJet80EtaBackward2EmbeddingEtaBackwardRatio;
    std::vector<TH1D*> hJet80EtaFBRatio2EmbeddingEtaFBRatio;
    std::vector<TH1D*> hJet80EtaBFRatio2EmbeddingEtaBFRatio;

    std::vector<TH1D*> hJet100Eta2EmbeddingEtaRatio;
    std::vector<TH1D*> hJet100EtaForward2EmbeddingEtaForwardRatio;
    std::vector<TH1D*> hJet100EtaBackward2EmbeddingEtaBackwardRatio;
    std::vector<TH1D*> hJet100EtaFBRatio2EmbeddingEtaFBRatio;
    std::vector<TH1D*> hJet100EtaBFRatio2EmbeddingEtaBFRatio;

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    Int_t dataStyle{2};
    Int_t mcStyle{0};

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        // Create histograms for ratios to gen
        hMBEta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hMBEta[i]->Clone( Form("%s_ratio2gen_%s", hMBEta[i]->GetName(), frame.Data()) ) ) );
        hMBEtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaForward[i]->Clone( Form("%s_ratio2gen_%s", hMBEtaForward[i]->GetName(), frame.Data()) ) ) );
        hMBEtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBackward[i]->Clone( Form("%s_ratio2gen_%s", hMBEtaBackward[i]->GetName(), frame.Data()) ) ) );
        hMBEtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio[i]->Clone( Form("%s_ratio2gen_%s", hMBEtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hMBEtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio[i]->Clone( Form("%s_ratio2gen_%s", hMBEtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hJet60Eta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hJet60Eta[i]->Clone( Form("%s_ratio2gen_%s", hJet60Eta[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaForward[i]->Clone( Form("%s_ratio2gen_%s", hJet60EtaForward[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward[i]->Clone( Form("%s_ratio2gen_%s", hJet60EtaBackward[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio[i]->Clone( Form("%s_ratio2gen_%s", hJet60EtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio[i]->Clone( Form("%s_ratio2gen_%s", hJet60EtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hJet80Eta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hJet80Eta[i]->Clone( Form("%s_ratio2gen_%s", hJet80Eta[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaForward[i]->Clone( Form("%s_ratio2gen_%s", hJet80EtaForward[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBackward[i]->Clone( Form("%s_ratio2gen_%s", hJet80EtaBackward[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaFBRatio[i]->Clone( Form("%s_ratio2gen_%s", hJet80EtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBFRatio[i]->Clone( Form("%s_ratio2gen_%s", hJet80EtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hJet100Eta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hJet100Eta[i]->Clone( Form("%s_ratio2gen_%s", hJet100Eta[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaForward[i]->Clone( Form("%s_ratio2gen_%s", hJet100EtaForward[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBackward[i]->Clone( Form("%s_ratio2gen_%s", hJet100EtaBackward[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaFBRatio[i]->Clone( Form("%s_ratio2gen_%s", hJet100EtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBFRatio[i]->Clone( Form("%s_ratio2gen_%s", hJet100EtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hEmbeddingEta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEta[i]->Clone( Form("%s_ratio2gen_%s", hEmbeddingEta[i]->GetName(), frame.Data()) ) ) );
        hEmbeddingEtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaForward[i]->Clone( Form("%s_ratio2gen_%s", hEmbeddingEtaForward[i]->GetName(), frame.Data()) ) ) );
        hEmbeddingEtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaBackward[i]->Clone( Form("%s_ratio2gen_%s", hEmbeddingEtaBackward[i]->GetName(), frame.Data()) ) ) );
        hEmbeddingEtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaFBRatio[i]->Clone( Form("%s_ratio2gen_%s", hEmbeddingEtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hEmbeddingEtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaBFRatio[i]->Clone( Form("%s_ratio2gen_%s", hEmbeddingEtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        // Create histograms for ratios to embedding
        hMBEta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hMBEta[i]->Clone( Form("%s_ratio2embedding_%s", hMBEta[i]->GetName(), frame.Data()) ) ) );
        hMBEtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaForward[i]->Clone( Form("%s_ratio2embedding_%s", hMBEtaForward[i]->GetName(), frame.Data()) ) ) );
        hMBEtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBackward[i]->Clone( Form("%s_ratio2embedding_%s", hMBEtaBackward[i]->GetName(), frame.Data()) ) ) );
        hMBEtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio[i]->Clone( Form("%s_ratio2embedding_%s", hMBEtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hMBEtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio[i]->Clone( Form("%s_ratio2embedding_%s", hMBEtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hJet60Eta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hJet60Eta[i]->Clone( Form("%s_ratio2embedding_%s", hJet60Eta[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaForward[i]->Clone( Form("%s_ratio2embedding_%s", hJet60EtaForward[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward[i]->Clone( Form("%s_ratio2embedding_p%s", hJet60EtaBackward[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio[i]->Clone( Form("%s_ratio2embedding_%s", hJet60EtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hJet60EtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio[i]->Clone( Form("%s_ratio2embedding_%s", hJet60EtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hJet80Eta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hJet80Eta[i]->Clone( Form("%s_ratio2embedding_%s", hJet80Eta[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaForward[i]->Clone( Form("%s_ratio2embedding_%s", hJet80EtaForward[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBackward[i]->Clone( Form("%s_ratio2embedding_%s", hJet80EtaBackward[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaFBRatio[i]->Clone( Form("%s_ratio2embedding_%s", hJet80EtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hJet80EtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBFRatio[i]->Clone( Form("%s_ratio2embedding_%s", hJet80EtaBFRatio[i]->GetName(), frame.Data()) ) ) );

        hJet100Eta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hJet100Eta[i]->Clone( Form("%s_ratio2embedding_%s", hJet100Eta[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaForward[i]->Clone( Form("%s_ratio2embedding_t_%s", hJet100EtaForward[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBackward[i]->Clone( Form("%s_ratio2embedding_%s", hJet100EtaBackward[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaFBRatio[i]->Clone( Form("%s_ratio2embedding_%s", hJet100EtaFBRatio[i]->GetName(), frame.Data()) ) ) );
        hJet100EtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBFRatio[i]->Clone( Form("%s_ratio2embedding_%s", hJet100EtaBFRatio[i]->GetName(), frame.Data()) ) ) );


        // Build the ratios
        hMBEta2GenEtaRatio[i]->Divide( hGenEta[i] );
        hMBEtaForward2GenEtaForwardRatio[i]->Divide( hGenEtaForward[i] );
        hMBEtaBackward2GenEtaBackwardRatio[i]->Divide( hGenEtaBackward[i] );
        hMBEtaFBRatio2GenEtaFBRatio[i]->Divide( hGenEtaFBRatio[i] );
        hMBEtaBFRatio2GenEtaBFRatio[i]->Divide( hGenEtaBFRatio[i] );

        hJet60Eta2GenEtaRatio[i]->Divide( hGenEta[i] );
        hJet60EtaForward2GenEtaForwardRatio[i]->Divide( hGenEtaForward[i] );
        hJet60EtaBackward2GenEtaBackwardRatio[i]->Divide( hGenEtaBackward[i] );
        hJet60EtaFBRatio2GenEtaFBRatio[i]->Divide( hGenEtaFBRatio[i] );
        hJet60EtaBFRatio2GenEtaBFRatio[i]->Divide( hGenEtaBFRatio[i] );

        hJet80Eta2GenEtaRatio[i]->Divide( hGenEta[i] );
        hJet80EtaForward2GenEtaForwardRatio[i]->Divide( hGenEtaForward[i] );
        hJet80EtaBackward2GenEtaBackwardRatio[i]->Divide( hGenEtaBackward[i] );
        hJet80EtaFBRatio2GenEtaFBRatio[i]->Divide( hGenEtaFBRatio[i] );
        hJet80EtaBFRatio2GenEtaBFRatio[i]->Divide( hGenEtaBFRatio[i] );

        hJet100Eta2GenEtaRatio[i]->Divide( hGenEta[i] );
        hJet100EtaForward2GenEtaForwardRatio[i]->Divide( hGenEtaForward[i] );
        hJet100EtaBackward2GenEtaBackwardRatio[i]->Divide( hGenEtaBackward[i] );
        hJet100EtaFBRatio2GenEtaFBRatio[i]->Divide( hGenEtaFBRatio[i] );
        hJet100EtaBFRatio2GenEtaBFRatio[i]->Divide( hGenEtaBFRatio[i] );

        hEmbeddingEta2GenEtaRatio[i]->Divide( hGenEta[i] );
        hEmbeddingEtaForward2GenEtaForwardRatio[i]->Divide( hGenEtaForward[i] );
        hEmbeddingEtaBackward2GenEtaBackwardRatio[i]->Divide( hGenEtaBackward[i] );
        hEmbeddingEtaFBRatio2GenEtaFBRatio[i]->Divide( hGenEtaFBRatio[i] );
        hEmbeddingEtaBFRatio2GenEtaBFRatio[i]->Divide( hGenEtaBFRatio[i] );

        hMBEta2EmbeddingEtaRatio[i]->Divide( hEmbeddingEta[i] );
        hMBEtaForward2EmbeddingEtaForwardRatio[i]->Divide( hEmbeddingEtaForward[i] );
        hMBEtaBackward2EmbeddingEtaBackwardRatio[i]->Divide( hEmbeddingEtaBackward[i] );
        hMBEtaFBRatio2EmbeddingEtaFBRatio[i]->Divide( hEmbeddingEtaFBRatio[i] );
        hMBEtaBFRatio2EmbeddingEtaBFRatio[i]->Divide( hEmbeddingEtaBFRatio[i] );

        hJet60Eta2EmbeddingEtaRatio[i]->Divide( hEmbeddingEta[i] );
        hJet60EtaForward2EmbeddingEtaForwardRatio[i]->Divide( hEmbeddingEtaForward[i] );
        hJet60EtaBackward2EmbeddingEtaBackwardRatio[i]->Divide( hEmbeddingEtaBackward[i] );
        hJet60EtaFBRatio2EmbeddingEtaFBRatio[i]->Divide( hEmbeddingEtaFBRatio[i] );
        hJet60EtaBFRatio2EmbeddingEtaBFRatio[i]->Divide( hEmbeddingEtaBFRatio[i] );

        hJet80Eta2EmbeddingEtaRatio[i]->Divide( hEmbeddingEta[i] );
        hJet80EtaForward2EmbeddingEtaForwardRatio[i]->Divide( hEmbeddingEtaForward[i] );
        hJet80EtaBackward2EmbeddingEtaBackwardRatio[i]->Divide( hEmbeddingEtaBackward[i] );
        hJet80EtaFBRatio2EmbeddingEtaFBRatio[i]->Divide( hEmbeddingEtaFBRatio[i] );
        hJet80EtaBFRatio2EmbeddingEtaBFRatio[i]->Divide( hEmbeddingEtaBFRatio[i] );

        hJet100Eta2EmbeddingEtaRatio[i]->Divide( hEmbeddingEta[i] );
        hJet100EtaForward2EmbeddingEtaForwardRatio[i]->Divide( hEmbeddingEtaForward[i] );
        hJet100EtaBackward2EmbeddingEtaBackwardRatio[i]->Divide( hEmbeddingEtaBackward[i] );
        hJet100EtaFBRatio2EmbeddingEtaFBRatio[i]->Divide( hEmbeddingEtaFBRatio[i] );
        hJet100EtaBFRatio2EmbeddingEtaBFRatio[i]->Divide( hEmbeddingEtaBFRatio[i] );

        // Set the style
        set1DStyle(hMBEta2GenEtaRatio[i], dataStyle);
        set1DStyle(hMBEtaForward2GenEtaForwardRatio[i], dataStyle);
        set1DStyle(hMBEtaBackward2GenEtaBackwardRatio[i], dataStyle);
        set1DStyle(hMBEtaFBRatio2GenEtaFBRatio[i], dataStyle);
        set1DStyle(hMBEtaBFRatio2GenEtaBFRatio[i], dataStyle);

        set1DStyle(hJet60Eta2GenEtaRatio[i], dataStyle);
        set1DStyle(hJet60EtaForward2GenEtaForwardRatio[i], dataStyle);
        set1DStyle(hJet60EtaBackward2GenEtaBackwardRatio[i], dataStyle);
        set1DStyle(hJet60EtaFBRatio2GenEtaFBRatio[i], dataStyle);
        set1DStyle(hJet60EtaBFRatio2GenEtaBFRatio[i], dataStyle);

        set1DStyle(hJet80Eta2GenEtaRatio[i], dataStyle);
        set1DStyle(hJet80EtaForward2GenEtaForwardRatio[i], dataStyle);
        set1DStyle(hJet80EtaBackward2GenEtaBackwardRatio[i], dataStyle);
        set1DStyle(hJet80EtaFBRatio2GenEtaFBRatio[i], dataStyle);
        set1DStyle(hJet80EtaBFRatio2GenEtaBFRatio[i], dataStyle);

        set1DStyle(hJet100Eta2GenEtaRatio[i], dataStyle);
        set1DStyle(hJet100EtaForward2GenEtaForwardRatio[i], dataStyle);
        set1DStyle(hJet100EtaBackward2GenEtaBackwardRatio[i], dataStyle);
        set1DStyle(hJet100EtaFBRatio2GenEtaFBRatio[i], dataStyle);
        set1DStyle(hJet100EtaBFRatio2GenEtaBFRatio[i], dataStyle);

        set1DStyle(hEmbeddingEta2GenEtaRatio[i], mcStyle);
        set1DStyle(hEmbeddingEtaForward2GenEtaForwardRatio[i], mcStyle);
        set1DStyle(hEmbeddingEtaBackward2GenEtaBackwardRatio[i], mcStyle);
        set1DStyle(hEmbeddingEtaFBRatio2GenEtaFBRatio[i], mcStyle);
        set1DStyle(hEmbeddingEtaBFRatio2GenEtaBFRatio[i], mcStyle);

        set1DStyle(hMBEta2EmbeddingEtaRatio[i], dataStyle);
        set1DStyle(hMBEtaForward2EmbeddingEtaForwardRatio[i], dataStyle);
        set1DStyle(hMBEtaBackward2EmbeddingEtaBackwardRatio[i], dataStyle);
        set1DStyle(hMBEtaFBRatio2EmbeddingEtaFBRatio[i], dataStyle);
        set1DStyle(hMBEtaBFRatio2EmbeddingEtaBFRatio[i], dataStyle);

        set1DStyle(hJet60Eta2EmbeddingEtaRatio[i], dataStyle);
        set1DStyle(hJet60EtaForward2EmbeddingEtaForwardRatio[i], dataStyle);
        set1DStyle(hJet60EtaBackward2EmbeddingEtaBackwardRatio[i], dataStyle);
        set1DStyle(hJet60EtaFBRatio2EmbeddingEtaFBRatio[i], dataStyle);
        set1DStyle(hJet60EtaBFRatio2EmbeddingEtaBFRatio[i], dataStyle);

        set1DStyle(hJet80Eta2EmbeddingEtaRatio[i], dataStyle);
        set1DStyle(hJet80EtaForward2EmbeddingEtaForwardRatio[i], dataStyle);
        set1DStyle(hJet80EtaBackward2EmbeddingEtaBackwardRatio[i], dataStyle);
        set1DStyle(hJet80EtaFBRatio2EmbeddingEtaFBRatio[i], dataStyle);
        set1DStyle(hJet80EtaBFRatio2EmbeddingEtaBFRatio[i], dataStyle);

        set1DStyle(hJet100Eta2EmbeddingEtaRatio[i], dataStyle);
        set1DStyle(hJet100EtaForward2EmbeddingEtaForwardRatio[i], dataStyle);
        set1DStyle(hJet100EtaBackward2EmbeddingEtaBackwardRatio[i], dataStyle);
        set1DStyle(hJet100EtaFBRatio2EmbeddingEtaFBRatio[i], dataStyle);
        set1DStyle(hJet100EtaBFRatio2EmbeddingEtaBFRatio[i], dataStyle);
    } // for (Int_t i{0}; i<ptBins; i++)

    hRatios.push_back( hMBEta2GenEtaRatio );
    hRatios.push_back( hMBEtaForward2GenEtaForwardRatio );
    hRatios.push_back( hMBEtaBackward2GenEtaBackwardRatio );
    hRatios.push_back( hMBEtaFBRatio2GenEtaFBRatio );
    hRatios.push_back( hMBEtaBFRatio2GenEtaBFRatio );

    hRatios.push_back( hJet60Eta2GenEtaRatio );
    hRatios.push_back( hJet60EtaForward2GenEtaForwardRatio );
    hRatios.push_back( hJet60EtaBackward2GenEtaBackwardRatio );
    hRatios.push_back( hJet60EtaFBRatio2GenEtaFBRatio );
    hRatios.push_back( hJet60EtaBFRatio2GenEtaBFRatio );

    hRatios.push_back( hJet80Eta2GenEtaRatio );
    hRatios.push_back( hJet80EtaForward2GenEtaForwardRatio );
    hRatios.push_back( hJet80EtaBackward2GenEtaBackwardRatio );
    hRatios.push_back( hJet80EtaFBRatio2GenEtaFBRatio );
    hRatios.push_back( hJet80EtaBFRatio2GenEtaBFRatio );

    hRatios.push_back( hJet100Eta2GenEtaRatio );
    hRatios.push_back( hJet100EtaForward2GenEtaForwardRatio );
    hRatios.push_back( hJet100EtaBackward2GenEtaBackwardRatio );
    hRatios.push_back( hJet100EtaFBRatio2GenEtaFBRatio );
    hRatios.push_back( hJet100EtaBFRatio2GenEtaBFRatio );

    hRatios.push_back( hEmbeddingEta2GenEtaRatio );
    hRatios.push_back( hEmbeddingEtaForward2GenEtaForwardRatio );
    hRatios.push_back( hEmbeddingEtaBackward2GenEtaBackwardRatio );
    hRatios.push_back( hEmbeddingEtaFBRatio2GenEtaFBRatio );
    hRatios.push_back( hEmbeddingEtaBFRatio2GenEtaBFRatio );

    hRatios.push_back( hMBEta2EmbeddingEtaRatio );
    hRatios.push_back( hMBEtaForward2EmbeddingEtaForwardRatio );
    hRatios.push_back( hMBEtaBackward2EmbeddingEtaBackwardRatio );
    hRatios.push_back( hMBEtaFBRatio2EmbeddingEtaFBRatio );
    hRatios.push_back( hMBEtaBFRatio2EmbeddingEtaBFRatio );

    hRatios.push_back( hJet60Eta2EmbeddingEtaRatio );
    hRatios.push_back( hJet60EtaForward2EmbeddingEtaForwardRatio );
    hRatios.push_back( hJet60EtaBackward2EmbeddingEtaBackwardRatio );
    hRatios.push_back( hJet60EtaFBRatio2EmbeddingEtaFBRatio );
    hRatios.push_back( hJet60EtaBFRatio2EmbeddingEtaBFRatio );

    hRatios.push_back( hJet80Eta2EmbeddingEtaRatio );
    hRatios.push_back( hJet80EtaForward2EmbeddingEtaForwardRatio );
    hRatios.push_back( hJet80EtaBackward2EmbeddingEtaBackwardRatio );
    hRatios.push_back( hJet80EtaFBRatio2EmbeddingEtaFBRatio );
    hRatios.push_back( hJet80EtaBFRatio2EmbeddingEtaBFRatio );

    hRatios.push_back( hJet100Eta2EmbeddingEtaRatio );
    hRatios.push_back( hJet100EtaForward2EmbeddingEtaForwardRatio );
    hRatios.push_back( hJet100EtaBackward2EmbeddingEtaBackwardRatio );
    hRatios.push_back( hJet100EtaFBRatio2EmbeddingEtaFBRatio );
    hRatios.push_back( hJet100EtaBFRatio2EmbeddingEtaBFRatio );

    return hRatios;
}

//________________
void appendVectorOfVectors(std::vector< std::vector<TH1D*> > &orig, std::vector< std::vector<TH1D*> > append) {
    for (Int_t i{0}; i<append.size(); i++) {
        orig.push_back( append.at(i) );
    }
}

//________________
std::vector< TF1* > fitSystRatio(std::vector< TH1D* > hSyst, Bool_t isUp = kTRUE, 
                                 Int_t fitType = 2, Int_t etaType = 0, Bool_t drawFits = kTRUE) {

    std::vector< TF1* > fitFunc( hSyst.size() );
    Double_t xMin{-3.}, xMax{3.};
    if (etaType != 0) {
        xMin = 0.0001;
        xMax = 3.;
    }

    // For all histograms in the vector
    for (Int_t i{0}; i<hSyst.size(); i++) {

        // Choose fit function
        if ( fitType == 0 ) {
            fitFunc.at(i) = new TF1(Form("%s_fitFunc", hSyst[i]->GetName()), "[0]", xMin, xMax);
            fitFunc.at(i)->SetParameter(0, 1.);
        }
        else if ( fitType == 1 ) {
            fitFunc.at(i) = new TF1(Form("%s_fitFunc", hSyst[i]->GetName()), "[0]+[1]*x", xMin, xMax);
            fitFunc.at(i)->SetParameters(1., 0.);
        }
        else if ( fitType == 2 ) {
            fitFunc.at(i) = new TF1(Form("%s_fitFunc", hSyst[i]->GetName()), "[0]+[1]*x+[2]*x*x", xMin, xMax);
            fitFunc.at(i)->SetParameters(1., 0., 0.);
        }
        else if ( fitType == 3 ) {
            fitFunc.at(i) = new TF1(Form("%s_fitFunc", hSyst[i]->GetName()), "[0]+[1]*x+[2]*x*x+[3]*x*x*x", xMin, xMax);
            fitFunc.at(i)->SetParameters(1., 0., 0., 0.);
        }
        else if ( fitType == 4 ) {
            fitFunc.at(i) = new TF1(Form("%s_fitFunc", hSyst[i]->GetName()), "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", xMin, xMax);
            fitFunc.at(i)->SetParameters(1., 0., 0., 0., 0.);
        }
        else {
            fitFunc.at(i) = new TF1(Form("%s_fitFunc", hSyst[i]->GetName()), "[0]+[1]*x", xMin, xMax);
            fitFunc.at(i)->SetParameters(1., 0.);
        }

        if ( isUp ) {
            fitFunc.at(i)->SetLineColor(kRed);
        }
        else {
            fitFunc.at(i)->SetLineColor(kBlue); 
        }
        fitFunc.at(i)->SetLineWidth(2);

        if ( drawFits ) {
            hSyst.at(i)->Fit( fitFunc.at(i), "MRE" );
        }
        else {
            hSyst.at(i)->Fit( fitFunc.at(i), "MRE0" );
        }

        // std::cout << Form("Fit function for %s: %s. chi2/ndf = %f", hSyst.at(i)->GetName(), fitFunc.at(i)->GetExpFormula().Data(), fitFunc.at(i)->GetChisquare() / fitFunc.at(i)->GetNDF() ) << std::endl;
    }

    return fitFunc;
}


//________________
std::vector< TH1D* > ratio2def(std::vector< TH1D* > hDef, std::vector< TH1D* > hSyst, 
                               Int_t systType = 0, Int_t upType = 0, Int_t etaType = 0, 
                               Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    // systType: 0 - jeu, 1 - jer, 2 - pointingRes, 3 - pileup                                
    // upType: 0 = up, 1 = down
    // etaType: 0 = eta, 1 = etaForward, 2 = etaBackward, 3 = etaFBRatio, 4 = etaBFRatio

    TString frame = (isCM) ? "cm" : "lab";
    std::vector< TH1D* > hRatio{};
    for (Int_t i{0}; i<hDef.size(); i++) {
        TH1D *h = dynamic_cast<TH1D*>( hSyst.at(i)->Clone( Form("%s_ratio2def_%s", hSyst.at(i)->GetName(), frame.Data()) ) );
        make1DRatio(h, hDef.at(i), Form("%s_ratio2def_%s", hSyst.at(i)->GetName(), frame.Data()) );
        set1DStyle(h, upType);
        hRatio.push_back(h);
    }

    Int_t fitType = 1;
    if (systType == 0) {
        if ( etaType == 0) {
            fitType = 4;
        }
        else {
            fitType = 2;
        }
    }
    else if (systType == 1) {
        if ( etaType == 0) {
            fitType = 2;
        }
        else {
            fitType = 1;
        }
        
    }
    else if (systType == 2) {
        if ( etaType == 0) {
            fitType = 2;
        }
        else {
            fitType = 1;
        }
    }
    else if (systType == 3) {
        if ( etaType == 0) {
            fitType = 1;
        }
        else {
            fitType = 1;
        }
    }

    std::vector< TF1* > fitFunc = fitSystRatio(hRatio, upType, fitType, etaType, drawFits);
    return hRatio;
}

//________________
std::vector< std::vector<TH1D*> > createVectorOfVariationRatio(std::vector< std::vector<TH1D*> > hDefDist, 
                                                               std::vector< std::vector<TH1D*> > hVarDist,
                                                               Int_t systType = 0,
                                                               Int_t upType = 0,
                                                               Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    // systType: 0 - jeu, 1 - jer, 2 - pointingRes, 3 - pileup
    // upType: 0 = up, 1 = down

    // Create vector of histograms to return
    std::vector< std::vector< TH1D* > > hRatios;

    // Vectors of default
    std::vector< TH1D* > hDefEta = hDefDist.at(0);
    std::vector< TH1D* > hDefEtaForward = hDefDist.at(1);
    std::vector< TH1D* > hDefEtaBackward = hDefDist.at(2);
    std::vector< TH1D* > hDefEtaFBRatio = hDefDist.at(3);
    std::vector< TH1D* > hDefEtaBFRatio = hDefDist.at(4);

    // Vectors of variation
    std::vector< TH1D* > hVarEta = hVarDist.at(0);
    std::vector< TH1D* > hVarEtaForward = hVarDist.at(1);
    std::vector< TH1D* > hVarEtaBackward = hVarDist.at(2);
    std::vector< TH1D* > hVarEtaFBRatio = hVarDist.at(3);
    std::vector< TH1D* > hVarEtaBFRatio = hVarDist.at(4);

    // Create histograms for ratios
    Int_t etaTypeDef = 0;
    Int_t etaTypeForward = 1;
    Int_t etaTypeBackward = 2;
    Int_t etaTypeFBRatio = 3;
    Int_t etaTypeBFRatio = 4;
    std::vector< TH1D* > hEtaVar2DefRatio = ratio2def(hDefEta, hVarEta, systType, upType, etaTypeDef, isCM, drawFits);
    std::vector< TH1D* > hEtaForwardVar2DefRatio = ratio2def(hDefEtaForward, hVarEtaForward, systType, upType, etaTypeForward, isCM, drawFits);
    std::vector< TH1D* > hEtaBackwardVar2DefRatio = ratio2def(hDefEtaBackward, hVarEtaBackward, systType, upType, etaTypeBackward, isCM, drawFits);
    std::vector< TH1D* > hEtaFBRatioVar2DefRatio = ratio2def(hDefEtaFBRatio, hVarEtaFBRatio, systType, upType, etaTypeFBRatio, isCM, drawFits);
    std::vector< TH1D* > hEtaBFRatioVar2DefRatio = ratio2def(hDefEtaBFRatio, hVarEtaBFRatio, systType, upType, etaTypeBFRatio, isCM, drawFits);

    hRatios.push_back( hEtaVar2DefRatio );
    hRatios.push_back( hEtaForwardVar2DefRatio );
    hRatios.push_back( hEtaBackwardVar2DefRatio );
    hRatios.push_back( hEtaFBRatioVar2DefRatio );
    hRatios.push_back( hEtaBFRatioVar2DefRatio );

    return hRatios;
}

//________________
std::vector< std::vector<TH1D*> > createVectorOfUpAndDownRatios(std::vector< std::vector<TH1D*> > hDefDist, 
                                                                std::vector< std::vector<TH1D*> > hUpDist, 
                                                                std::vector< std::vector<TH1D*> > hDownDist,
                                                                Int_t systType = 0,
                                                                Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {
                                                                  
    // Create vector of histograms to return
    std::vector< std::vector< TH1D* > > hRatios;

    // systType: 0 - jeu, 1 - jer, 2 - pointingRes, 3 - pileup
    
    Int_t upType = 0;
    Int_t downType = 1;

    // Vectors of default
    std::vector< std::vector<TH1D* > > hUp2DefRatios = createVectorOfVariationRatio(hDefDist, hUpDist, systType, upType, isCM, drawFits);
    std::vector< std::vector<TH1D* > > hDown2DefRatios = {};

    if ( !hDownDist.empty() ) {
        hDown2DefRatios = createVectorOfVariationRatio(hDefDist, hDownDist, systType, downType, isCM, drawFits);
    }
    
    // Append the vectors
    appendVectorOfVectors(hRatios, hUp2DefRatios);
    if ( !hDownDist.empty() ) {
        appendVectorOfVectors(hRatios, hDown2DefRatios);
    }

    return hRatios;
}


//________________
std::vector< std::vector<TH1D*> > createUpAndDown2Def(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hJeuMBUpEtaDist, std::vector< std::vector<TH1D*> > hJeuMBDownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hJet60EtaDist, std::vector< std::vector<TH1D*> > hJeuJet60UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet60DownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJeuJet80UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet80DownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hJet100EtaDist, std::vector< std::vector<TH1D*> > hJeuJet100UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet100DownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hEmbeddingRefEtaDist, std::vector< std::vector<TH1D*> > hJerUpEtaDist, std::vector< std::vector<TH1D*> > hJerDownEtaDist,
                                                      std::vector< std::vector<TH1D*> > hMBGplusEtaDist, std::vector< std::vector<TH1D*> > hMBVtx1EtaDist,
                                                      std::vector< std::vector<TH1D*> > hJet60GplusEtaDist, std::vector< std::vector<TH1D*> > hJet60Vtx1EtaDist,
                                                      std::vector< std::vector<TH1D*> > hJet80GplusEtaDist, std::vector< std::vector<TH1D*> > hJet80Vtx1EtaDist,
                                                      std::vector< std::vector<TH1D*> > hJet100GplusEtaDist, std::vector< std::vector<TH1D*> > hJet100Vtx1EtaDist,
                                                      Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    // Create vector of histograms to return
    std::vector< std::vector<TH1D*> > hRatios;
    std::vector< std::vector<TH1D*> > hEmpty = {};

    Int_t jeuSystType = 0;
    Int_t jerSystType = 1;
    Int_t pointingResSystType = 2;
    Int_t pileupSystType = 3;

    // Retrieve histograms to divide
    std::vector< std::vector<TH1D*> > hMBEtaJeuRatios = createVectorOfUpAndDownRatios(hMBEtaDist, hJeuMBUpEtaDist, hJeuMBDownEtaDist, jeuSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hJet60EtaJeuRatios = createVectorOfUpAndDownRatios(hJet60EtaDist, hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, jeuSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hJet80EtaJeuRatios = createVectorOfUpAndDownRatios(hJet80EtaDist, hJeuJet80UpEtaDist, hJeuJet80DownEtaDist, jeuSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hJet100EtaJeuRatios = createVectorOfUpAndDownRatios(hJet100EtaDist, hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, jeuSystType, isCM, drawFits);

    std::vector< std::vector<TH1D*> > hEmbeddingJerRatios = createVectorOfUpAndDownRatios(hEmbeddingEtaDist, hJerUpEtaDist, hJerDownEtaDist, jerSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hEmbeddingPointingResRatios = createVectorOfUpAndDownRatios(hEmbeddingEtaDist, hEmbeddingRefEtaDist, hEmpty, pointingResSystType, isCM, drawFits);

    std::vector< std::vector<TH1D*> > hMBPileup = createVectorOfUpAndDownRatios(hMBEtaDist, hMBGplusEtaDist, hMBVtx1EtaDist, pileupSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hJet60Pileup = createVectorOfUpAndDownRatios(hJet60EtaDist, hJet60GplusEtaDist, hJet60Vtx1EtaDist, pileupSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hJet80Pileup = createVectorOfUpAndDownRatios(hJet80EtaDist, hJet80GplusEtaDist, hJet80Vtx1EtaDist, pileupSystType, isCM, drawFits);
    std::vector< std::vector<TH1D*> > hJet100Pileup = createVectorOfUpAndDownRatios(hJet100EtaDist, hJet100GplusEtaDist, hJet100Vtx1EtaDist, pileupSystType, isCM, drawFits);

    // Append the vectors
    appendVectorOfVectors(hRatios, hMBEtaJeuRatios);      // 0
    appendVectorOfVectors(hRatios, hJet60EtaJeuRatios);   // 10
    appendVectorOfVectors(hRatios, hJet80EtaJeuRatios);   // 20
    appendVectorOfVectors(hRatios, hJet100EtaJeuRatios);  // 30

    appendVectorOfVectors(hRatios, hEmbeddingJerRatios);  // 40

    appendVectorOfVectors(hRatios, hMBPileup);            // 50
    appendVectorOfVectors(hRatios, hJet60Pileup);         // 60
    appendVectorOfVectors(hRatios, hJet80Pileup);         // 70
    appendVectorOfVectors(hRatios, hJet100Pileup);        // 80

    appendVectorOfVectors(hRatios, hEmbeddingPointingResRatios); // 90-95


    return hRatios;
}

//________________
/*
* Create a vector of histograms with dijet eta distributions from a 3D histogram for the given file
*/
std::vector< std::vector<TH1D*> >createDijetEtaHistograms(TFile* inFile, const char* inputHistoName, 
                                                          const char* hName, bool isCM = false, int defType = 2) {

    TString baseHistoName = inputHistoName;
    TString fullEtaHistName = baseHistoName;
    TString forwardEtaHistName = baseHistoName;
    TString backwardEtaHistName = baseHistoName;
    if ( isCM ) {
        fullEtaHistName += "DphiCMWeighted";
        forwardEtaHistName += "CMForwardWeighted";
        backwardEtaHistName += "CMBackwardWeighted";
    }
    else {
        fullEtaHistName += "DphiWeighted";
        forwardEtaHistName += "ForwardWeighted";
        backwardEtaHistName += "BackwardWeighted";
    }

    // std::cout << "[Main creator] Input histo name: " << inputHistoName << std::endl;
    // std::cout << "[Main creator] Output histo name: " << hName << std::endl;
    // std::cout << Form("\tFull: %s \n\tForward: %s \n\tBackward: %s", fullEtaHistName.Data(), forwardEtaHistName.Data(), backwardEtaHistName.Data() ) << std::endl;

    Int_t upType{0};
    Int_t downType{1};

    // Create vector of histograms to return
    std::vector< std::vector<TH1D*> > hDijetEtaHistos;
    std::vector<TH1D*> hEtaForward;
    std::vector<TH1D*> hEtaBackward;
    std::vector< std::vector<TH1D*> > hEtaDir;
    if ( isCM ) {
        hEtaForward = createDijetDirectionEtaHistograms(inFile, forwardEtaHistName.Data(), hName, upType, isCM, kTRUE);
        hEtaBackward = createDijetDirectionEtaHistograms(inFile, backwardEtaHistName.Data(), hName, downType, isCM, kFALSE);
    }
    else {
        hEtaDir = createDijetDirectionEtaHistogramsFromLab(inFile, fullEtaHistName.Data(), hName, defType);
        hEtaForward = hEtaDir[0];
        hEtaBackward = hEtaDir[1];
    }
    std::vector<TH1D*> hEta = createDijetEtaFullHistograms(inFile, fullEtaHistName.Data(), hName, defType, isCM);
    std::vector<TH1D*> hEtaFBRatio = createFBRatios(hEtaForward, hEtaBackward, hName, upType, isCM, kTRUE);
    std::vector<TH1D*> hEtaBFRatio = createFBRatios(hEtaForward, hEtaBackward, hName, downType, isCM, kFALSE);
   

    hDijetEtaHistos.push_back(hEta);
    hDijetEtaHistos.push_back(hEtaForward);
    hDijetEtaHistos.push_back(hEtaBackward);
    hDijetEtaHistos.push_back(hEtaFBRatio);
    hDijetEtaHistos.push_back(hEtaBFRatio);

    return hDijetEtaHistos;
}

//________________
void compareJetCollections(TFile *ak4, TFile *akCs4, TString date) {

    TString trigName = "MB";
    TString fName = ak4->GetName();
    if ( fName.Contains("MB") ) {
        trigName = "MB";
    }
    else if ( fName.Contains("Jet60") ) {
        trigName = "Jet60";
    }
    else if ( fName.Contains("Jet80") ) {
        trigName = "Jet80";
    }
    else if ( fName.Contains("Jet100") ) {
        trigName = "Jet100";
    }
    else if ( fName.Contains("Jet120") ) {
        trigName = "Jet120";
    }
    else {
        trigName = "MB";
    }

    TH3D *hPtEtaDphiAk4 = (TH3D*)ak4->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiAk4->SetName("hPtEtaDphiAk4");
    TH3D *hPtEtaDphiAkCs4 = (TH3D*)akCs4->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiAkCs4->SetName("hPtEtaDphiAkCs4");

    TString frame = "lab";

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};

    // Styles
    Int_t ak4Type{0};
    Int_t akCs4Type{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaAk4[ ptDijetBinLow.size() ];
    TH1D *hEtaAkCs4[ ptDijetBinLow.size() ];
    TH1D *hEtaRat[ ptDijetBinLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};

    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetBinLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // ak4
        hEtaAk4[i] = projectEtaFrom3D(hPtEtaDphiAk4, Form("hEtaAk4_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaAk4[i] );
        set1DStyle( hEtaAk4[i], ak4Type );

        // akCs4
        hEtaAkCs4[i] = projectEtaFrom3D(hPtEtaDphiAkCs4, Form("hEtaAkCs4_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaAkCs4[i] );
        set1DStyle( hEtaAkCs4[i], akCs4Type );

        // Ratio of akCs4 to ak4
        hEtaRat[i] = dynamic_cast<TH1D*>( hEtaAkCs4[i]->Clone( Form("hEtaRat_%d", i) ) );
        make1DRatio(hEtaRat[i], hEtaAk4[i], Form("hEtaRat_%d", i) );
        hEtaRat[i]->GetYaxis()->SetTitle("akCs4 / ak4");

        // Plot comparison
        canv->cd();
        setPadStyle();
        hEtaAk4[i]->Draw();
        hEtaAkCs4[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaAk4[i], Form("ak4"), "p");
        leg->AddEntry(hEtaAkCs4[i], Form("akCs4"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/%s_pPb8160_etaDijet_jetCollComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratio
        canv->cd();
        setPadStyle();
        hEtaRat[i]->Draw();
        hEtaRat[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaRat[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRat[i]->GetXaxis()->GetBinUpEdge(hEtaRat[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/%s_pPb8160_etaDijet_jetCollRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hEtaAk4[i]->Draw();
        hEtaAkCs4[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaAk4[i], Form("ak4"), "p");
        leg->AddEntry(hEtaAkCs4[i], Form("akCs4"), "p");
        leg->Draw();

        // Plot all ratios of one canvas
        cRat->cd(i+1);
        setPadStyle();
        hEtaRat[i]->Draw();
        hEtaRat[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaRat[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRat[i]->GetXaxis()->GetBinUpEdge(hEtaRat[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
    } // for (Int_t i{0}; i<ptBins; i++)

    cComp->SaveAs( Form("%s/%s_pPb8160_etaDijet_jetCollComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/%s_pPb8160_etaDijet_jetCollRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
}

//________________
TString returnTrigName(TFile *f) {
    TString trigName = "MB";
    TString fName = f->GetName();
    if ( fName.Contains("MB") ) {
        trigName = "MB";
    }
    else if ( fName.Contains("Jet60") ) {
        trigName = "Jet60";
    }
    else if ( fName.Contains("Jet80") ) {
        trigName = "Jet80";
    }
    else if ( fName.Contains("Jet100") ) {
        trigName = "Jet100";
    }
    else if ( fName.Contains("Jet120") ) {
        trigName = "Jet120";
    }
    else if ( fName.Contains("Embed") ) {
        trigName = "Embed";
    }
    else {
        trigName = "MB";
    }
    return trigName;
}

//________________
TString returnTrigName(const char *name) {
    TString trigName = "MB";
    TString hName = name;
    if ( hName.Contains("MB") ) {
        trigName = "MB";
    }
    else if ( hName.Contains("Jet60") ) {
        trigName = "Jet60";
    }
    else if ( hName.Contains("Jet80") ) {
        trigName = "Jet80";
    }
    else if ( hName.Contains("Jet100") ) {
        trigName = "Jet100";
    }
    else if ( hName.Contains("Jet120") ) {
        trigName = "Jet120";
    }
    else if ( hName.Contains("Embed") ) {
        trigName = "Embed";
    }
    else {
        trigName = "MB";
    }
    return trigName;
}

//________________
TGraphErrors* retrieveTotalAbsSystUncrt(TH1D *hDef, TH1D *hRelSyst, const Char_t *setName = "MBEta", Int_t iPt = 0) {
    TGraphErrors *gr = new TGraphErrors();
    TString name = Form("gr%s_totalAbsUncrt_%d", setName, iPt);
    gr->SetName( name.Data() );
    for (Int_t i{1}; i<=hDef->GetNbinsX(); i++) {
        Double_t x{0}, y{0};
        x = hDef->GetBinCenter(i);
        y = hDef->GetBinContent(i);
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, hDef->GetBinWidth(i) / 2, y * hRelSyst->GetBinContent(i));
    }
    //gr->SetFillColorAlpha(29, 0.55);
    gr->SetLineColor(2);
    gr->SetFillColor(2);
    gr->SetFillStyle(1001);
    return gr;
}

//________________
void plotDistributionWithUncrt(TCanvas *c, TH1D *h1, TH1D *h2, 
                              int ptLow = 50, int ptHi = 60,
                              Bool_t isCM = kFALSE, Int_t fbType = 0) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.001, 0.12};

    if ( fbType == 1 || fbType == 2 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.001;
        yRange[1] = 0.2;
    }
    else if ( fbType == 3 ) {               // Forward/backward
        xRange[0] = 0.;
        xRange[1] = 2.4;
        yRange[0] = 0.75;
        yRange[1] = 1.5;
    }
    else if ( fbType == 4 ) {               // Backward/forward
        xRange[0] = 0.;
        xRange[1] = 2.4;
        yRange[0] = 0.4;
        yRange[1] = 1.15;
    }

    // Set style for the data points
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.3);
    h1->SetMarkerColor(kBlack);
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(2);

    // Set style for systematic uncertainty
    setSystUncrtStyle(h2, 0);

    // Create pad
    TLegend *leg;

    // Set pad style
    setPadStyle();

    // Draw histograms
    h2->Draw("E2");
    h1->Draw("same");
    h2->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h2->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    t.DrawLatexNDC(0.35, 0.83, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi ) );
    
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
    t.DrawLatexNDC(0.2, 0.66, "p_{T}^{Subleading} > 40 GeV");
    if ( isCM ) {
        t.DrawLatexNDC(0.2, 0.57, "|#eta_{CM}| < 2.4");
    }
    else {
        t.DrawLatexNDC(0.2, 0.57, "|#eta| < 3");
    }
    t.DrawLatexNDC(0.2, 0.49, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
    t.SetTextSize(0.06);

    leg = new TLegend(0.6, 0.65, 0.8, 0.8);
    leg->SetTextSize(0.06);
    leg->SetLineWidth(0);
    leg->AddEntry(h1, "Data", "lep");
    leg->AddEntry(h2, "Syst. uncrt.", "f");
    leg->Draw();

    plotCMSHeader();
}

//________________
void plotManyDistributionsOnCanvas(TCanvas *c, std::vector< TH1D* > hData, std::vector< TH1D* > hSyst, 
                                   Bool_t isCM = kFALSE, Int_t fbType = 0) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.001, 0.12};

    if ( fbType == 1 || fbType == 2 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.001;
        yRange[1] = 0.2;
    }
    else if ( fbType == 3 ) {               // Forward/backward
        xRange[0] = 0.;
        xRange[1] = 2.4;
        yRange[0] = 0.75;
        yRange[1] = 1.5;
    }
    else if ( fbType == 4 ) {               // Backward/forward
        xRange[0] = 0.;
        xRange[1] = 2.4;
        yRange[0] = 0.4;
        yRange[1] = 1.15;
    }

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Set style for the data points
    for (Int_t i{0}; i<hData.size(); i++) {
        hData.at(i)->SetMarkerStyle(20);
        hData.at(i)->SetMarkerSize(1.3);
        hData.at(i)->SetMarkerColor(kBlack);
        hData.at(i)->SetLineColor(kBlack);
        hData.at(i)->SetLineWidth(2);

        setSystUncrtStyle( hSyst.at(i), 0 );
    }

    // Create legend
    TLegend *leg;

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        c->cd(i+1);
        setPadStyle();
        hSyst.at(i)->Draw("E2");
        hData.at(i)->Draw("same");
        // hData.at(i)->Draw();
        // hSyst.at(i)->Draw("E2 same");
        hSyst.at(i)->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hSyst.at(i)->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        t.DrawLatexNDC(0.33, 0.83, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetLow.at(i), ptDijetHi.at(i) ) );
        
        if ( i == 0 ) {
            t.SetTextSize(0.05);
            t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
            t.DrawLatexNDC(0.2, 0.66, "p_{T}^{Subleading} > 40 GeV");

            if ( isCM ) {
                t.DrawLatexNDC(0.2, 0.57, "|#eta_{CM}| < 2.4");
            }
            else {
                t.DrawLatexNDC(0.2, 0.57, "|#eta| < 3");
            }

            t.DrawLatexNDC(0.2, 0.49, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
            t.SetTextSize(0.06);
        }

        if ( i == 1 ) {
            leg = new TLegend(0.6, 0.65, 0.8, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hData.at(i), "Data", "lep");
            leg->AddEntry(hSyst.at(i), "Syst. uncrt.", "f");
            leg->Draw();
        }

        plotCMSHeader();
    } // for (Int_t i{0}; i<ptBins; i++)
}

//________________
void plotFinalEtaDistributions(std::vector< std::vector<TH1D*> > hFinalDist, 
                               std::vector< std::vector<TH1D*> > hFinalAbsSystUncrtDist,
                               TString date, Bool_t isCM = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Obtain vectors with distributions
    std::vector< TH1D* > hEtaDist = hFinalDist.at(0);
    std::vector< TH1D* > hEtaForwardDist = hFinalDist.at(1);
    std::vector< TH1D* > hEtaBackwardDist = hFinalDist.at(2);
    std::vector< TH1D* > hEtaFBRatioDist = hFinalDist.at(3);
    std::vector< TH1D* > hEtaBFRatioDist = hFinalDist.at(4);

    std::vector< TH1D* > hEtaAbsSystUncrtDist = hFinalAbsSystUncrtDist.at(0);
    std::vector< TH1D* > hEtaForwardAbsSystUncrtDist = hFinalAbsSystUncrtDist.at(1);
    std::vector< TH1D* > hEtaBackwardAbsSystUncrtDist = hFinalAbsSystUncrtDist.at(2);
    std::vector< TH1D* > hEtaFBRatioAbsSystUncrtDist = hFinalAbsSystUncrtDist.at(3);
    std::vector< TH1D* > hEtaBFRatioAbsSystUncrtDist = hFinalAbsSystUncrtDist.at(4);

    // Create canvases
    Int_t sizeX{1200}, sizeY{1200};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);

    TCanvas *cEta = new TCanvas("cEta", "cEta", sizeX, sizeY);
    cEta->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaForward = new TCanvas("cEtaForward", "cEtaForward", sizeX, sizeY);
    cEtaForward->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaBackward = new TCanvas("cEtaBackward", "cEtaBackward", sizeX, sizeY);
    cEtaBackward->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaFB = new TCanvas("cEtaFB", "cEtaFB", sizeX, sizeY);
    cEtaFB->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaBF = new TCanvas("cEtaBF", "cEtaBF", sizeX, sizeY);
    cEtaBF->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    // Plot distributions on canvases
    for (Int_t i=0; i<ptBins; i++) {

        // Full eta
        canv->cd();
        plotDistributionWithUncrt(canv, hEtaDist.at(i), hEtaAbsSystUncrtDist.at(i), ptDijetLow.at(i), ptDijetHi.at(i), isCM, 0);
        canv->SaveAs( Form("%s/data/pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ));

        // Forward eta
        canv->cd();
        plotDistributionWithUncrt(canv, hEtaForwardDist.at(i), hEtaForwardAbsSystUncrtDist.at(i), ptDijetLow.at(i), ptDijetHi.at(i), isCM, 1);
        canv->SaveAs( Form("%s/data/pPb8160_etaDijet_forward_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ));
        
        // Backward eta
        canv->cd();
        plotDistributionWithUncrt(canv, hEtaBackwardDist.at(i), hEtaBackwardAbsSystUncrtDist.at(i), ptDijetLow.at(i), ptDijetHi.at(i), isCM, 2);
        canv->SaveAs( Form("%s/data/pPb8160_etaDijet_backward_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ));

        // Forward/backward ratio
        canv->cd();
        plotDistributionWithUncrt(canv, hEtaFBRatioDist.at(i), hEtaFBRatioAbsSystUncrtDist.at(i), ptDijetLow.at(i), ptDijetHi.at(i), isCM, 3);
        canv->SaveAs( Form("%s/data/pPb8160_etaDijet_fbRatio_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ));

        // Backward/forward ratio
        canv->cd();
        plotDistributionWithUncrt(canv, hEtaBFRatioDist.at(i), hEtaBFRatioAbsSystUncrtDist.at(i), ptDijetLow.at(i), ptDijetHi.at(i), isCM, 4);
        canv->SaveAs( Form("%s/data/pPb8160_etaDijet_bfRatio_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ));
    } // for (Int_t i=0; i<ptBins; i++)

    plotManyDistributionsOnCanvas(cEta, hEtaDist, hEtaAbsSystUncrtDist, isCM, 0);
    plotManyDistributionsOnCanvas(cEtaForward, hEtaForwardDist, hEtaForwardAbsSystUncrtDist, isCM, 1);
    plotManyDistributionsOnCanvas(cEtaBackward, hEtaBackwardDist, hEtaBackwardAbsSystUncrtDist, isCM, 2);
    plotManyDistributionsOnCanvas(cEtaFB, hEtaFBRatioDist, hEtaFBRatioAbsSystUncrtDist, isCM, 3);
    plotManyDistributionsOnCanvas(cEtaBF, hEtaBFRatioDist, hEtaBFRatioAbsSystUncrtDist, isCM, 4);

    cEta->SaveAs( Form("%s/data/pPb8160_etaDijet_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaForward->SaveAs( Form("%s/data/pPb8160_etaDijet_forward_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaBackward->SaveAs( Form("%s/data/pPb8160_etaDijet_backward_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaFB->SaveAs( Form("%s/data/pPb8160_etaDijet_fb_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaBF->SaveAs( Form("%s/data/pPb8160_etaDijet_bf_all_%s.pdf", date.Data(), frame.Data() ) );
}

//________________
void plotIndividualUpDownDefComparison(TCanvas *c, TH1D *hDef, TH1D *hUp, TH1D *hDown,
                                       Int_t ptLow = 50, Int_t ptHi = 60, 
                                       Bool_t isCM = kFALSE, Int_t fbType = 0) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.001, 0.12};
    Double_t legX[2] = {0.27, 0.47};
    Double_t legY[2] = {0.65, 0.85};

    if ( fbType == 1 || fbType == 2 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.001;
        yRange[1] = 0.2;
    }
    else if ( fbType == 3 ) {               // Forward/backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.75;
        yRange[1] = 1.5;
    }
    else if ( fbType == 4 ) {               // Backward/forward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.4;
        yRange[1] = 1.15;
        
        legY[0] = 0.27;
        legY[1] = 0.47;
    }

    // Set style for the data points
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};
    set1DStyle( hDef, defType );
    set1DStyle( hUp, upType );
    if ( hDown ) {
        set1DStyle( hDown, downType );
    }

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot comparison
    //
    c->cd();

    // Set pad style
    if ( fbType == 0 ) {
        setPadStyle();
    }
    else {
        setPadFBStyle();
    }

    // Draw histograms
    hDef->Draw();
    hUp->Draw("same");
    if ( hDown) {
        hDown->Draw("same");
    }
    hDef->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hDef->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi ) );
    
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetTextSize(0.04);
    leg->SetLineWidth(0);
    leg->AddEntry(hDef, "Default", "p");
    leg->AddEntry(hUp, "Up", "p");
    if ( hDown ) {
        leg->AddEntry(hDown, "Down", "p");
    }
    leg->Draw();
}

//________________
void plotIndividualUpDownDefRatio(TCanvas *c, TH1D *hRatioUp, TH1D *hRatioDown,
                                  Int_t ptLow = 50, Int_t ptHi = 60, 
                                  Bool_t isCM = kFALSE, Int_t fbType = 0) {

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.75, 1.25};
    Double_t legX[2] = {0.27, 0.47};
    Double_t legY[2] = {0.65, 0.85};

    if ( fbType == 1 || fbType == 2 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.001;
        yRange[1] = 0.2;
    }
    else if ( fbType == 3 ) {               // Forward/backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.75;
        yRange[1] = 1.5;
    }
    else if ( fbType == 4 ) {               // Backward/forward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.4;
        yRange[1] = 1.15;
        
        legY[0] = 0.27;
        legY[1] = 0.47;
    }

    // Set style for the data points
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};
    set1DStyle( hRatioUp, upType );
    if ( hRatioDown ) {
        set1DStyle( hRatioDown, downType );
    }

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot ratios
    //
    c->cd();

    // Set pad style
    if ( fbType == 0 ) {
        setPadStyle();
    }
    else {
        setPadFBStyle();
    }

    // Draw histograms
    hRatioUp->Draw();
    if ( hRatioDown ) {
        hRatioDown->Draw("same");
    }
    hRatioUp->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hRatioUp->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi ) );
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetTextSize(0.04);
    leg->SetLineWidth(0);
    leg->AddEntry(hRatioUp, "Up/Def", "p");
    if ( hRatioDown ) {
        leg->AddEntry(hRatioDown, "Down/Def", "p");
    }
    leg->Draw();
    line = new TLine(xRange[0], 1., xRange[1], 1.);
    line->SetLineColor(kMagenta);
    line->SetLineWidth(3);
    line->SetLineStyle(3);
    line->Draw();
}

//________________
void plotSystComparison(TCanvas *c, std::vector<TH1D*> hDef, 
                        std::vector<TH1D*> hUp, 
                        std::vector<TH1D*> hDown,
                        Int_t fbType = 0) {

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.001, 0.12};
    Double_t legX[2] = {0.27, 0.47};
    Double_t legY[2] = {0.65, 0.85};

    if ( fbType == 1 || fbType == 2 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.001;
        yRange[1] = 0.2;
    }
    else if ( fbType == 3 ) {               // Forward/backward
        xRange[0] = 0.;
        xRange[1] = 2.4;
        yRange[0] = 0.75;
        yRange[1] = 1.5;
    }
    else if ( fbType == 4 ) {               // Backward/forward
        xRange[0] = 0.;
        xRange[1] = 2.4;
        yRange[0] = 0.4;
        yRange[1] = 1.15;

        legY[0] = 0.27;
        legY[1] = 0.47;
    }

    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    TLegend *leg;

    // std::cout << "hDef.size() = " << hDef.size() << std::endl;
    // std::cout << "hUp.size() = " << hUp.size() << std::endl;
    // if ( !hDown.empty() ) {
    //     std::cout << "hDown.size() = " << hDown.size() << std::endl;
    // }

    // Plot comparison
    for (Int_t i{0}; i<hDef.size(); i++) {

        // std::cout << "i = " << i << std::endl;
        // std::cout << "hDef.at(i) name = " << hDef.at(i)->GetName() << std::endl;
        // std::cout << "hUp.at(i) name = " << hUp.at(i)->GetName() << std::endl;
        // if ( !hDown.empty() ) {
        //     std::cout << "hDown.at(i) name = " << hDown.at(i)->GetName() << std::endl;
        // }

        c->cd(i+1);

        if ( fbType == 0 ) {
            setPadStyle();
        }
        else {
            setPadFBStyle();
        }

        set1DStyle( hDef.at(i), defType );
        set1DStyle( hUp.at(i), upType );
        if ( !hDown.empty() ) {
            set1DStyle( hDown.at(i), downType );
        }

        hDef.at(i)->Draw();
        hUp.at(i)->Draw("same");
        if ( !hDown.empty() ) {
            hDown.at(i)->Draw("same");
        }
        hDef.at(i)->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hDef.at(i)->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        t.DrawLatexNDC( 0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );

        if ( i == 0 ) {
            TLegend *leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
            leg->SetTextSize(0.04);
            leg->SetLineWidth(0);
            leg->AddEntry(hDef.at(i), "Default", "p");
            leg->AddEntry(hUp.at(i), "Up", "p");
            if ( !hDown.empty() ) {
                leg->AddEntry(hDown.at(i), "Down", "p");
            }
            leg->Draw();
        }
    } // for (Int_t i{0}; i<hDef.size(); i++)
}

//________________
void plotSystRatio(TCanvas *c, std::vector<TH1D*> hRatioUp, std::vector<TH1D*> hRatioDown, Int_t fbType = 0) {
    
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.7, 1.3};

    if ( fbType != 0 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.75;
        yRange[1] = 1.25;
    }

    Int_t upType{0};
    Int_t downType{1};

    TLegend *leg;
    TLine *line;

    // std::cout << "hRatioUp.size() = " << hRatioUp.size() << std::endl;
    // if ( !hRatioDown.empty() ) {
    //     std::cout << "hRatioDown.size() = " << hRatioDown.size() << std::endl;
    // }

    // Plot comparison
    for (Int_t i{0}; i<hRatioUp.size(); i++) {

        // std::cout << "i = " << i << std::endl;
        // std::cout << "hRatioUp.at(i) name = " << hRatioUp.at(i)->GetName() << std::endl;
        // if ( !hRatioDown.empty() ) {
        //     std::cout << "hRatioDown.at(i) name = " << hRatioDown.at(i)->GetName() << std::endl;
        // }

        // Switch to the pad
        c->cd(i+1);

        // Set pad style
        if ( fbType == 0 ) {
            setPadStyle();
        }
        else {
            setPadFBStyle();
        }

        // Set histogram style
        set1DStyle( hRatioUp.at(i), upType );
        if ( !hRatioDown.empty() ) {
            set1DStyle( hRatioDown.at(i), downType );
        }

        // Draw histograms
        hRatioUp.at(i)->Draw();
        TF1 *fitUp = hRatioUp.at(i)->GetFunction( Form("%s_fitFunc", hRatioUp.at(i)->GetName()) );
        fitUp->SetLineColor(kRed);
        fitUp->Draw("same");
        if ( !hRatioDown.empty() ) {
            hRatioDown.at(i)->Draw("same");
            TF1 *fitDown = hRatioDown.at(i)->GetFunction( Form("%s_fitFunc", hRatioDown.at(i)->GetName()) );
            fitDown->SetLineColor(kBlue);
            fitDown->Draw("same");
        }
        // Adjust ranges
        hRatioUp.at(i)->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        hRatioUp.at(i)->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hRatioUp.at(i)->GetYaxis()->SetTitle("Ratio to default");
        // Draw pT range
        t.DrawLatexNDC( 0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );

        // Draw line at y=1
        line = new TLine(xRange[0], 1., xRange[1], 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();

        // Draw legend in the first pad
        if ( i == 0 ) {
            TLegend *leg = new TLegend(0.27, 0.7, 0.47, 0.85);
            leg->SetTextSize(0.04);
            leg->SetLineWidth(0);
            leg->AddEntry(hRatioUp.at(i), "Up/Def", "p");
            if ( !hRatioDown.empty() ) {
                leg->AddEntry(hRatioDown.at(i), "Down/Def", "p");
            }
            leg->Draw();
        }
    } // for (Int_t i{0}; i<hRatioUp.size(); i++)
}

//________________
void plotUp2DownComparison(std::vector<TH1D*> hDef, std::vector<TH1D*> hUp, std::vector<TH1D*> hDown,
                           std::vector<TH1D*> hRatioUp, std::vector<TH1D*> hRatioDown,
                           Int_t systematics, Int_t etaType, Int_t fitType, TString date, 
                           Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Read systematics type
    TString systType;
    if ( systematics == 0 ) {
        systType = "jes";
    }
    else if ( systematics == 1 ) {
        systType = "jer";
    }
    else if ( systematics == 2 ) {
        systType = "pileup";
    }
    else if ( systematics == 3 ) {
        systType = "pointinResolution";
    }
    else {
        drawFits = kFALSE;
        systType = "data2mc";
    }

    // Direction name
    TString dirName;
    // 0 - symmetric, 1 - forward, backward, fb, bf
    if (etaType == 0) {
        dirName = "";
    }
    else if (etaType == 1) {
        dirName = "forward";
    }
    else if (etaType == 2) {
        dirName = "backward";
    }
    else if (etaType == 3) {
        dirName = "fb";
    }
    else if (etaType == 4) {
        dirName = "bf";
    }

    // Fits (Fit type: 0 -constant, 1 - linear, 2 - quadratic, 3 - cubic, 4 - quartic)
    // std::vector<TF1*> fitUp = fitSystRatio(hRatioUp, kTRUE, fitType, etaType, drawFits);
    // std::vector<TF1*> fitDown = {};
    // if ( !hRatioDown.empty() ) {
    //     fitDown = fitSystRatio(hRatioDown, kFALSE, fitType, etaType, drawFits);
    // }

    // Create canvas
    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", sizeX, sizeY);

    TCanvas *cComp = new TCanvas(Form("%s_%sComp_%s", returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data() ), 
                                 Form("%s_%sComp_%s", returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data() ), 
                                 sizeX, sizeY);
    cComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cRat = new TCanvas(Form("%s_%sRat_%s", returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data() ), 
                                Form("%s_%sRat_%s", returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data() ), 
                                sizeX, sizeY);
    cRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    // Fill canvases with many histograms
    plotSystComparison(cComp, hDef, hUp, hDown, etaType);
    plotSystRatio(cRat, hRatioUp, hRatioDown, etaType);

    // Save canvases
    cComp->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sComp_%sall_%s.pdf", date.Data(), systType.Data(),  
                   returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sRat_%sall_%s.pdf", date.Data(), systType.Data(), 
                   returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data(), frame.Data() ) );

    // Loop over pT average bins
    for (Int_t i{0}; i<ptDijetBinLow.size(); i++) {

        TH1D *hDownIndividual = nullptr;
        if ( !hDown.empty() ) {
            hDownIndividual = hDown.at(i);
        }

        // Plot and save individual comparisons
        plotIndividualUpDownDefComparison(canv, hDef.at(i), hUp.at(i), hDownIndividual, 
                                          ptDijetLow.at(i), ptDijetHi.at(i), kFALSE, etaType);
        canv->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sComp_%s_pt_%d_%d_%s.pdf", date.Data(), systType.Data(), 
                           returnTrigName( hDef.at(i)->GetName() ).Data(), systType.Data(), dirName.Data(),
                           ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

        hDownIndividual = nullptr;
        if ( !hRatioDown.empty() ) {
            hDownIndividual = hRatioDown.at(i);
        }

        // Plot and save individual ratios
        plotIndividualUpDownDefRatio(canv, hRatioUp.at(i), hDownIndividual, 
                                     ptDijetLow.at(i), ptDijetHi.at(i), kFALSE, etaType);
        canv->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sRat_%s_pt_%d_%d_%s.pdf", date.Data(), systType.Data(), 
                           returnTrigName( hDef.at(i)->GetName() ).Data(), systType.Data(), dirName.Data(),
                           ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );
    } // for (Int_t i{0}; i<ptDijetBinLow.size(); i++)
}

//________________
void plotData2McComparison(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist, 
                           std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist,
                           std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hGenEtaDist,
                           std::vector< std::vector<TH1D*> > hRatios2McDist,
                           TString date, Bool_t isCM = kFALSE) {

    Bool_t plotMcReco = {kFALSE};

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.001, 0.12};

    Double_t xRangeForward[2] = {0., 3.};
    Double_t yRangeForward[2] = {0.001, 0.2};

    Double_t xRangeFB[2] = {0., 2.4};
    Double_t yRangeFB[2] = {0.75, 1.5};
    Double_t yRangeBF[2] = {0.4, 1.15};

    Double_t yRangeRat[2] = {0.7, 1.25};

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Data
    std::vector< TH1D* > hMBEta = hMBEtaDist.at(0); std::vector< TH1D* > hMBEtaForward = hMBEtaDist.at(1); 
    std::vector< TH1D* > hMBEtaBackward = hMBEtaDist.at(2); std::vector< TH1D* > hMBEtaFBRatio = hMBEtaDist.at(3);
    std::vector< TH1D* > hMBEtaBFRatio = hMBEtaDist.at(4);

    std::vector< TH1D* > hJet60Eta = hJet60EtaDist.at(0); std::vector< TH1D* > hJet60EtaForward = hJet60EtaDist.at(1);
    std::vector< TH1D* > hJet60EtaBackward = hJet60EtaDist.at(2); std::vector< TH1D* > hJet60EtaFBRatio = hJet60EtaDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatio = hJet60EtaDist.at(4);

    std::vector< TH1D* > hJet80Eta = hJet80EtaDist.at(0); std::vector< TH1D* > hJet80EtaForward = hJet80EtaDist.at(1);
    std::vector< TH1D* > hJet80EtaBackward = hJet80EtaDist.at(2); std::vector< TH1D* > hJet80EtaFBRatio = hJet80EtaDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatio = hJet80EtaDist.at(4);

    std::vector< TH1D* > hJet100Eta = hJet100EtaDist.at(0); std::vector< TH1D* > hJet100EtaForward = hJet100EtaDist.at(1);
    std::vector< TH1D* > hJet100EtaBackward = hJet100EtaDist.at(2); std::vector< TH1D* > hJet100EtaFBRatio = hJet100EtaDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatio = hJet100EtaDist.at(4);

    // Embedding
    std::vector<TH1D*> hEmbeddingEta = hEmbeddingEtaDist.at(0); std::vector<TH1D*> hEmbeddingEtaForward = hEmbeddingEtaDist.at(1);
    std::vector<TH1D*> hEmbeddingEtaBackward = hEmbeddingEtaDist.at(2); std::vector<TH1D*> hEmbeddingEtaFBRatio = hEmbeddingEtaDist.at(3);
    std::vector<TH1D*> hEmbeddingEtaBFRatio = hEmbeddingEtaDist.at(4);

    // Gen
    std::vector<TH1D*> hGenEta = hGenEtaDist.at(0); std::vector<TH1D*> hGenEtaForward = hGenEtaDist.at(1);
    std::vector<TH1D*> hGenEtaBackward = hGenEtaDist.at(2); std::vector<TH1D*> hGenEtaFBRatio = hGenEtaDist.at(3);
    std::vector<TH1D*> hGenEtaBFRatio = hGenEtaDist.at(4);

    // Ratios
    std::vector<TH1D*> hMBEta2GenEtaRatio = hRatios2McDist.at(0);
    std::vector<TH1D*> hMBEtaForward2GenEtaRatio = hRatios2McDist.at(1);
    std::vector<TH1D*> hMBEtaBackward2GenEtaRatio = hRatios2McDist.at(2);
    std::vector<TH1D*> hMBEtaFBRatio2GenEtaRatio = hRatios2McDist.at(3);
    std::vector<TH1D*> hMBEtaBFRatio2GenEtaRatio = hRatios2McDist.at(4);

    std::vector<TH1D*> hJet60Eta2GenEtaRatio = hRatios2McDist.at(5);
    std::vector<TH1D*> hJet60EtaForward2GenEtaRatio = hRatios2McDist.at(6);
    std::vector<TH1D*> hJet60EtaBackward2GenEtaRatio = hRatios2McDist.at(7);
    std::vector<TH1D*> hJet60EtaFBRatio2GenEtaRatio = hRatios2McDist.at(8);
    std::vector<TH1D*> hJet60EtaBFRatio2GenEtaRatio = hRatios2McDist.at(9);

    std::vector<TH1D*> hJet80Eta2GenEtaRatio = hRatios2McDist.at(10);
    std::vector<TH1D*> hJet80EtaForward2GenEtaRatio = hRatios2McDist.at(11);
    std::vector<TH1D*> hJet80EtaBackward2GenEtaRatio = hRatios2McDist.at(12);
    std::vector<TH1D*> hJet80EtaFBRatio2GenEtaRatio = hRatios2McDist.at(13);
    std::vector<TH1D*> hJet80EtaBFRatio2GenEtaRatio = hRatios2McDist.at(14);

    std::vector<TH1D*> hJet100Eta2GenEtaRatio = hRatios2McDist.at(15);
    std::vector<TH1D*> hJet100EtaForward2GenEtaRatio = hRatios2McDist.at(16);
    std::vector<TH1D*> hJet100EtaBackward2GenEtaRatio = hRatios2McDist.at(17);
    std::vector<TH1D*> hJet100EtaFBRatio2GenEtaRatio = hRatios2McDist.at(18);
    std::vector<TH1D*> hJet100EtaBFRatio2GenEtaRatio = hRatios2McDist.at(19);

    std::vector<TH1D*> hEmbedding2GenEtaRatio = hRatios2McDist.at(20);
    std::vector<TH1D*> hEmbedding2GenEtaForwardRatio = hRatios2McDist.at(21);
    std::vector<TH1D*> hEmbedding2GenEtaBackwardRatio = hRatios2McDist.at(22);
    std::vector<TH1D*> hEmbedding2GenEtaFBRatio = hRatios2McDist.at(23);
    std::vector<TH1D*> hEmbedding2GenEtaBFRatio = hRatios2McDist.at(24);

    std::vector<TH1D*> hMBEta2EmbeddingEtaRatio = hRatios2McDist.at(25);
    std::vector<TH1D*> hMBEtaForward2EmbeddingEtaRatio = hRatios2McDist.at(26);
    std::vector<TH1D*> hMBEtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(27);
    std::vector<TH1D*> hMBEtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(28);
    std::vector<TH1D*> hMBEtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(29);

    std::vector<TH1D*> hJet60Eta2EmbeddingEtaRatio = hRatios2McDist.at(30);
    std::vector<TH1D*> hJet60EtaForward2EmbeddingEtaRatio = hRatios2McDist.at(31);
    std::vector<TH1D*> hJet60EtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(32);
    std::vector<TH1D*> hJet60EtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(33);
    std::vector<TH1D*> hJet60EtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(34);

    std::vector<TH1D*> hJet80Eta2EmbeddingEtaRatio = hRatios2McDist.at(35);
    std::vector<TH1D*> hJet80EtaForward2EmbeddingEtaRatio = hRatios2McDist.at(36);
    std::vector<TH1D*> hJet80EtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(37);
    std::vector<TH1D*> hJet80EtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(38);
    std::vector<TH1D*> hJet80EtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(39);

    std::vector<TH1D*> hJet100Eta2EmbeddingEtaRatio = hRatios2McDist.at(40);
    std::vector<TH1D*> hJet100EtaForward2EmbeddingEtaRatio = hRatios2McDist.at(41);
    std::vector<TH1D*> hJet100EtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(42);
    std::vector<TH1D*> hJet100EtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(43);
    std::vector<TH1D*> hJet100EtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(44);

    for (Int_t i{0}; i<ptBins; i++) {
        hMBEtaFBRatio[i]->SetLineColor(kBlack);
        hMBEtaFBRatio[i]->SetMarkerColor(kBlack);
        hMBEtaFBRatio[i]->SetMarkerStyle(20);
        hMBEtaBFRatio[i]->SetLineColor(kBlack);
        hMBEtaBFRatio[i]->SetMarkerColor(kBlack);
        hMBEtaBFRatio[i]->SetMarkerStyle(20);

        hJet60EtaFBRatio[i]->SetLineColor(kBlack);
        hJet60EtaFBRatio[i]->SetMarkerColor(kBlack);
        hJet60EtaFBRatio[i]->SetMarkerStyle(20);
        hJet60EtaBFRatio[i]->SetLineColor(kBlack);
        hJet60EtaBFRatio[i]->SetMarkerColor(kBlack);
        hJet60EtaBFRatio[i]->SetMarkerStyle(20);

        hJet80EtaFBRatio[i]->SetLineColor(kBlack);
        hJet80EtaFBRatio[i]->SetMarkerColor(kBlack);
        hJet80EtaFBRatio[i]->SetMarkerStyle(20);
        hJet80EtaBFRatio[i]->SetLineColor(kBlack);
        hJet80EtaBFRatio[i]->SetMarkerColor(kBlack);
        hJet80EtaBFRatio[i]->SetMarkerStyle(20);
        
        hJet100EtaFBRatio[i]->SetLineColor(kBlack);
        hJet100EtaFBRatio[i]->SetMarkerColor(kBlack);
        hJet100EtaFBRatio[i]->SetMarkerStyle(20);
        hJet100EtaBFRatio[i]->SetLineColor(kBlack);
        hJet100EtaBFRatio[i]->SetMarkerColor(kBlack);
        hJet100EtaBFRatio[i]->SetMarkerStyle(20);

        hEmbeddingEta[i]->SetLineColor(kBlue);
        hEmbeddingEta[i]->SetMarkerColor(kBlue);
        hEmbeddingEtaFBRatio[i]->SetLineColor(kBlue);
        hEmbeddingEtaFBRatio[i]->SetMarkerColor(kBlue);
        hEmbeddingEtaBFRatio[i]->SetLineColor(kBlue);
        hEmbeddingEtaBFRatio[i]->SetMarkerColor(kBlue);

        hGenEta[i]->SetLineColor(kMagenta);
        hGenEta[i]->SetMarkerColor(kMagenta);
        hGenEtaFBRatio[i]->SetLineColor(kMagenta);
        hGenEtaFBRatio[i]->SetMarkerColor(kMagenta);
        hGenEtaBFRatio[i]->SetLineColor(kMagenta);
        hGenEtaBFRatio[i]->SetMarkerColor(kMagenta);

        hEmbedding2GenEtaRatio[i]->SetLineColor(kBlue);
        hEmbedding2GenEtaRatio[i]->SetMarkerColor(kBlue);
    }

    Int_t sizeX{1200}, sizeY{1200};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);

    TCanvas *cEta = new TCanvas("cEta", "cEta", sizeX, sizeY);
    cEta->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaForwardAndBackward = new TCanvas("cEtaForwardAndBackward", "cEtaForwardAndBackward", sizeX, sizeY);
    cEtaForwardAndBackward->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaFB = new TCanvas("cEtaFB", "cEtaFB", sizeX, sizeY);
    cEtaFB->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEtaBF = new TCanvas("cEtaBF", "cEtaBF", sizeX, sizeY);
    cEtaBF->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TCanvas *cEta2Gen = new TCanvas("cEta2Gen", "cEta2Gen", sizeX, sizeY);
    cEta2Gen->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1), 0.001, 0.001);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.07);

    TLegend *leg;
    TLine *line;

    Double_t drawPt[2] = {0.35, 0.83};

    // Loop over pT average bins
    for (Int_t i{0}; i<ptBins; i++) {

        if (ptDijetHi.at(i) < 80) {

            canv->cd();
            setPadStyle();
            hMBEta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            hMBEta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hMBEta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            plotCMSHeader();
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hMBEta[i], "MinBias", "p");
            leg->AddEntry(hEmbeddingEta[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            canv->SaveAs(Form("%s/data2mc/MB_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data()));

            cEta->cd(i + 1);
            setPadStyle();
            hMBEta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            hMBEta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hMBEta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            if (i == 0) {
                leg = new TLegend(0.3, 0.6, 0.5, 0.8);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEta[i], "Data", "p");
                leg->AddEntry(hEmbeddingEta[i], "PYTHIA (no nPDF)", "p");
                leg->Draw();
            }
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hMBEtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hMBEtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hMBEtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);

            leg->AddEntry(hMBEtaFBRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaFBRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaFBRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs(Form("%s/data2mc/MB_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data()));

            cEtaFB->cd(i + 1);
            setPadFBStyle();
            hMBEtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hMBEtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hMBEtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            if (i == 0) {
                leg = new TLegend(0.3, 0.6, 0.5, 0.8);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEtaFBRatio[i], "Data", "p");
                if (plotMcReco)
                {
                    leg->AddEntry(hEmbeddingEtaFBRatio[i], "PYTHIA reco (no nPDF)", "p");
                }
                leg->AddEntry(hGenEtaFBRatio[i], "PYTHIA (no nPDF)", "p");
                leg->Draw();
            }
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hMBEtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hMBEtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hMBEtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hMBEtaBFRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaBFRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaBFRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs(Form("%s/data2mc/MB_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data()));

            cEtaBF->cd(i + 1);
            setPadFBStyle();
            hMBEtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hMBEtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hMBEtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            if (i == 0)
            {
                leg = new TLegend(0.3, 0.6, 0.5, 0.8);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEtaBFRatio[i], "Data", "p");
                if (plotMcReco)
                {
                    leg->AddEntry(hEmbeddingEtaBFRatio[i], "PYTHIA reco (no nPDF)", "p");
                }
                leg->AddEntry(hGenEtaBFRatio[i], "PYTHIA reco (no nPDF)", "p");
                leg->Draw();
            }
            plotCMSHeader();

            // Ratios to gen
            canv->cd();
            setPadStyle();
            hMBEta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hMBEta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hMBEta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            leg = new TLegend(0.4, 0.25, 0.6, 0.45);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hMBEta2GenEtaRatio[i], "Data/PYTHIA8", "p");
            leg->AddEntry(hEmbedding2GenEtaRatio[i], "PYTHIA8 reco/PYTHIA8", "p");
            leg->Draw();
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
            canv->SaveAs(Form("%s/data2mc/MB_pPb8160_etaDijet2Gen_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data()));

            cEta2Gen->cd(i + 1);
            setPadStyle();
            hMBEta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hMBEta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hMBEta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i)));
            if (i == 0) {
                leg = new TLegend(0.4, 0.25, 0.6, 0.45);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEta2GenEtaRatio[i], "Reco(data)/Gen", "p");
                leg->AddEntry(hEmbedding2GenEtaRatio[i], "Reco(MC)/Gen", "p");
                leg->Draw();
            }
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
        }
        else if (ptDijetHi.at(i) < 100) {
            canv->cd();
            setPadStyle();
            hJet60Eta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            // hGenEta[i]->Draw("same");
            hJet60Eta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet60Eta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet60Eta[i], "Data", "p");
            leg->AddEntry(hEmbeddingEta[i], "PYTHIA reco (no nPDF)", "p");
            // leg->AddEntry(hGenEta[i], "PYTHIA (no nPDF)", "p");
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet60_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            hJet60Eta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            // hGenEta[i]->Draw("same");
            hJet60Eta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet60Eta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet60EtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hJet60EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet60EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet60EtaFBRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaFBRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaFBRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet60_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            hJet60EtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hGenEtaFBRatio[i]->Draw("same");
            hJet60EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet60EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet60EtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hJet60EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet60EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet60EtaBFRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaBFRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaBFRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet60_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hJet60EtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hJet60EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet60EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            // Ratios
            canv->cd();
            setPadStyle();
            hJet60Eta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hJet60Eta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet60Eta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.4, 0.25, 0.6, 0.45);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet60Eta2GenEtaRatio[i], "Reco(data)/Gen", "p");
            leg->AddEntry(hEmbedding2GenEtaRatio[i], "Reco(MC)/Gen", "p");
            leg->Draw();
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet60_pPb8160_etaDijet2Gen_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta2Gen->cd(i+1);
            setPadStyle();
            hJet60Eta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hJet60Eta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet60Eta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();

        }
        else if (ptDijetHi.at(i) < 120) {
            canv->cd();
            setPadStyle();
            hJet80Eta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            //hGenEta[i]->Draw("same");
            hJet80Eta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet80Eta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet80Eta[i], "Data", "p");
            leg->AddEntry(hEmbeddingEta[i], "PYTHIA reco (no nPDF)", "p");
            // leg->AddEntry(hGenEta[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet80_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            hJet80Eta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            // hGenEta[i]->Draw("same");
            hJet80Eta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet80Eta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet80EtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hJet80EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet80EtaFBRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaFBRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaFBRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet80_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            hJet80EtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hJet80EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet80EtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hJet80EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet80EtaBFRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaBFRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaBFRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet80_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hJet80EtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hJet80EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            // Ratios
            canv->cd();
            setPadStyle();
            hJet80Eta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hJet80Eta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet80Eta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.4, 0.25, 0.6, 0.45);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet80Eta2GenEtaRatio[i], "Reco(data)/Gen", "p");
            leg->AddEntry(hEmbedding2GenEtaRatio[i], "Reco(MC)/Gen", "p");
            leg->Draw();
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet80_pPb8160_etaDijet2Gen_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta2Gen->cd(i+1);
            setPadStyle();
            hJet80Eta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hJet80Eta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet80Eta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
        }
        else {
            canv->cd();
            setPadStyle();
            hJet100Eta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            hJet100Eta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet100Eta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet100Eta[i], "Data", "p");
            leg->AddEntry(hEmbeddingEta[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet100_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            hJet100Eta[i]->Draw();
            hEmbeddingEta[i]->Draw("same");
            hJet100Eta[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet100Eta[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet100EtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hJet100EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet100EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet100EtaFBRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaFBRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaFBRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet100_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            hJet100EtaFBRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaFBRatio[i]->Draw("same");
            }
            hGenEtaFBRatio[i]->Draw("same");
            hJet100EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet100EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet100EtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hJet100EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet100EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.3, 0.6, 0.5, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet100EtaBFRatio[i], "Data", "p");
            if (plotMcReco) {
                leg->AddEntry(hEmbeddingEtaBFRatio[i], "PYTHIA reco (no nPDF)", "p");
            }
            leg->AddEntry(hGenEtaBFRatio[i], "PYTHIA (no nPDF)", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet100_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hJet100EtaBFRatio[i]->Draw();
            if (plotMcReco) {
                hEmbeddingEtaBFRatio[i]->Draw("same");
            }
            hGenEtaBFRatio[i]->Draw("same");
            hJet100EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet100EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            // Ratios
            canv->cd();
            setPadStyle();
            hJet100Eta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hJet100Eta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet100Eta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.4, 0.25, 0.6, 0.45);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet100Eta2GenEtaRatio[i], "Reco(data)/Gen", "p");
            leg->AddEntry(hEmbedding2GenEtaRatio[i], "Reco(MC)/Gen", "p");
            leg->Draw();
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data2mc/Jet100_pPb8160_etaDijet2Gen_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta2Gen->cd(i+1);
            setPadStyle();
            hJet100Eta2GenEtaRatio[i]->Draw();
            hEmbedding2GenEtaRatio[i]->Draw("same");
            hJet100Eta2GenEtaRatio[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            hJet100Eta2GenEtaRatio[i]->GetYaxis()->SetRangeUser(yRangeRat[0], yRangeRat[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            line = new TLine(xRange[0], 1., xRange[1], 1.);
            line->SetLineColor(kMagenta);
            line->SetLineStyle(2);
            line->Draw();
            plotCMSHeader();
        }
    } // for (Int_t i=0; i<ptBins; i++)

    cEta->SaveAs( Form("%s/data2mc/etaDijetComparison_%s.pdf", date.Data(), frame.Data() ) );
    cEtaFB->SaveAs( Form("%s/data2mc/etaDijetComparison_fb_%s.pdf", date.Data(), frame.Data() ) );
    cEtaBF->SaveAs( Form("%s/data2mc/etaDijetComparison_bf_%s.pdf", date.Data(), frame.Data() ) );
    cEta2Gen->SaveAs( Form("%s/data2mc/etaDijetComparison2Gen_%s.pdf", date.Data(), frame.Data() ) );
}


//________________
void plotJES(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist, 
             std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist, 
             std::vector< std::vector<TH1D*> > hJeuMBUpEtaDist, std::vector< std::vector<TH1D*> > hJeuMBDownEtaDist, 
             std::vector< std::vector<TH1D*> > hJeuJet60UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet60DownEtaDist, 
             std::vector< std::vector<TH1D*> > hJeuJet80UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet80DownEtaDist, 
             std::vector< std::vector<TH1D*> > hJeuJet100UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet100DownEtaDist,
             std::vector< std::vector<TH1D*> > hMBEtaJeuUp2DefDist, std::vector< std::vector<TH1D*> > hMBEtaJeuDown2DefDist,
             std::vector< std::vector<TH1D*> > hJet60EtaJeuUp2DefDist, std::vector< std::vector<TH1D*> > hJet60EtaJeuDown2DefDist,
             std::vector< std::vector<TH1D*> > hJet80EtaJeuUp2DefDist, std::vector< std::vector<TH1D*> > hJet80EtaJeuDown2DefDist, 
             std::vector< std::vector<TH1D*> > hJet100EtaJeuUp2DefDist, std::vector< std::vector<TH1D*> > hJet100EtaJeuDown2DefDist,
             TString date, Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Data
    std::vector< TH1D* > hMBEta = hMBEtaDist.at(0); std::vector< TH1D* > hMBEtaForward = hMBEtaDist.at(1); 
    std::vector< TH1D* > hMBEtaBackward = hMBEtaDist.at(2); std::vector< TH1D* > hMBEtaFBRatio = hMBEtaDist.at(3);
    std::vector< TH1D* > hMBEtaBFRatio = hMBEtaDist.at(4);

    std::vector< TH1D* > hJet60Eta = hJet60EtaDist.at(0); std::vector< TH1D* > hJet60EtaForward = hJet60EtaDist.at(1);
    std::vector< TH1D* > hJet60EtaBackward = hJet60EtaDist.at(2); std::vector< TH1D* > hJet60EtaFBRatio = hJet60EtaDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatio = hJet60EtaDist.at(4);

    std::vector< TH1D* > hJet80Eta = hJet80EtaDist.at(0); std::vector< TH1D* > hJet80EtaForward = hJet80EtaDist.at(1);
    std::vector< TH1D* > hJet80EtaBackward = hJet80EtaDist.at(2); std::vector< TH1D* > hJet80EtaFBRatio = hJet80EtaDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatio = hJet80EtaDist.at(4);

    std::vector< TH1D* > hJet100Eta = hJet100EtaDist.at(0); std::vector< TH1D* > hJet100EtaForward = hJet100EtaDist.at(1);
    std::vector< TH1D* > hJet100EtaBackward = hJet100EtaDist.at(2); std::vector< TH1D* > hJet100EtaFBRatio = hJet100EtaDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatio = hJet100EtaDist.at(4);

    // Jet energy scale
    std::vector< TH1D* > hJeuMBUpEta = hJeuMBUpEtaDist.at(0); std::vector< TH1D* > hJeuMBUpEtaForward = hJeuMBUpEtaDist.at(1);
    std::vector< TH1D* > hJeuMBUpEtaBackward = hJeuMBUpEtaDist.at(2); std::vector< TH1D* > hJeuMBUpEtaFBRatio = hJeuMBUpEtaDist.at(3);
    std::vector< TH1D* > hJeuMBUpEtaBFRatio = hJeuMBUpEtaDist.at(4);

    std::vector< TH1D* > hJeuMBDownEta = hJeuMBDownEtaDist.at(0); std::vector< TH1D* > hJeuMBDownEtaForward = hJeuMBDownEtaDist.at(1);
    std::vector< TH1D* > hJeuMBDownEtaBackward = hJeuMBDownEtaDist.at(2); std::vector< TH1D* > hJeuMBDownEtaFBRatio = hJeuMBDownEtaDist.at(3);
    std::vector< TH1D* > hJeuMBDownEtaBFRatio = hJeuMBDownEtaDist.at(4);

    std::vector< TH1D* > hJeuJet60UpEta = hJeuJet60UpEtaDist.at(0); std::vector< TH1D* > hJeuJet60UpEtaForward = hJeuJet60UpEtaDist.at(1);
    std::vector< TH1D* > hJeuJet60UpEtaBackward = hJeuJet60UpEtaDist.at(2); std::vector< TH1D* > hJeuJet60UpEtaFBRatio = hJeuJet60UpEtaDist.at(3);
    std::vector< TH1D* > hJeuJet60UpEtaBFRatio = hJeuJet60UpEtaDist.at(4);

    std::vector< TH1D* > hJeuJet60DownEta = hJeuJet60DownEtaDist.at(0); std::vector< TH1D* > hJeuJet60DownEtaForward = hJeuJet60DownEtaDist.at(1);
    std::vector< TH1D* > hJeuJet60DownEtaBackward = hJeuJet60DownEtaDist.at(2); std::vector< TH1D* > hJeuJet60DownEtaFBRatio = hJeuJet60DownEtaDist.at(3);
    std::vector< TH1D* > hJeuJet60DownEtaBFRatio = hJeuJet60DownEtaDist.at(4);

    std::vector< TH1D* > hJeuJet80UpEta = hJeuJet80UpEtaDist.at(0); std::vector< TH1D* > hJeuJet80UpEtaForward = hJeuJet80UpEtaDist.at(1);
    std::vector< TH1D* > hJeuJet80UpEtaBackward = hJeuJet80UpEtaDist.at(2); std::vector< TH1D* > hJeuJet80UpEtaFBRatio = hJeuJet80UpEtaDist.at(3);
    std::vector< TH1D* > hJeuJet80UpEtaBFRatio = hJeuJet80UpEtaDist.at(4);

    std::vector< TH1D* > hJeuJet80DownEta = hJeuJet80DownEtaDist.at(0); std::vector< TH1D* > hJeuJet80DownEtaForward = hJeuJet80DownEtaDist.at(1);
    std::vector< TH1D* > hJeuJet80DownEtaBackward = hJeuJet80DownEtaDist.at(2); std::vector< TH1D* > hJeuJet80DownEtaFBRatio = hJeuJet80DownEtaDist.at(3);
    std::vector< TH1D* > hJeuJet80DownEtaBFRatio = hJeuJet80DownEtaDist.at(4);

    std::vector< TH1D* > hJeuJet100UpEta = hJeuJet100UpEtaDist.at(0); std::vector< TH1D* > hJeuJet100UpEtaForward = hJeuJet100UpEtaDist.at(1);  
    std::vector< TH1D* > hJeuJet100UpEtaBackward = hJeuJet100UpEtaDist.at(2); std::vector< TH1D* > hJeuJet100UpEtaFBRatio = hJeuJet100UpEtaDist.at(3);
    std::vector< TH1D* > hJeuJet100UpEtaBFRatio = hJeuJet100UpEtaDist.at(4);

    std::vector< TH1D* > hJeuJet100DownEta = hJeuJet100DownEtaDist.at(0); std::vector< TH1D* > hJeuJet100DownEtaForward = hJeuJet100DownEtaDist.at(1);
    std::vector< TH1D* > hJeuJet100DownEtaBackward = hJeuJet100DownEtaDist.at(2); std::vector< TH1D* > hJeuJet100DownEtaFBRatio = hJeuJet100DownEtaDist.at(3);
    std::vector< TH1D* > hJeuJet100DownEtaBFRatio = hJeuJet100DownEtaDist.at(4);

    // Ratios
    std::vector< TH1D* > hMBJeuUpToDefEta = hMBEtaJeuUp2DefDist.at(0); std::vector< TH1D* > hMBJeuUpToDefEtaForward = hMBEtaJeuUp2DefDist.at(1);
    std::vector< TH1D* > hMBJeuUpToDefEtaBackward = hMBEtaJeuUp2DefDist.at(2); std::vector< TH1D* > hMBJeuUpToDefEtaFBRatio = hMBEtaJeuUp2DefDist.at(3);
    std::vector< TH1D* > hMBJeuUpToDefEtaBFRatio = hMBEtaJeuUp2DefDist.at(4);

    std::vector< TH1D* > hMBJeuDownToDefEta = hMBEtaJeuDown2DefDist.at(0); std::vector< TH1D* > hMBJeuDownToDefEtaForward = hMBEtaJeuDown2DefDist.at(1);
    std::vector< TH1D* > hMBJeuDownToDefEtaBackward = hMBEtaJeuDown2DefDist.at(2); std::vector< TH1D* > hMBJeuDownToDefEtaFBRatio = hMBEtaJeuDown2DefDist.at(3);
    std::vector< TH1D* > hMBJeuDownToDefEtaBFRatio = hMBEtaJeuDown2DefDist.at(4);

    std::vector< TH1D* > hJet60JeuUpToDefEta = hJet60EtaJeuUp2DefDist.at(0); std::vector< TH1D* > hJet60JeuUpToDefEtaForward = hJet60EtaJeuUp2DefDist.at(1);
    std::vector< TH1D* > hJet60JeuUpToDefEtaBackward = hJet60EtaJeuUp2DefDist.at(2); std::vector< TH1D* > hJet60JeuUpToDefEtaFBRatio = hJet60EtaJeuUp2DefDist.at(3);
    std::vector< TH1D* > hJet60JeuUpToDefEtaBFRatio = hJet60EtaJeuUp2DefDist.at(4);

    std::vector< TH1D* > hJet60JeuDownToDefEta = hJet60EtaJeuDown2DefDist.at(0); std::vector< TH1D* > hJet60JeuDownToDefEtaForward = hJet60EtaJeuDown2DefDist.at(1);
    std::vector< TH1D* > hJet60JeuDownToDefEtaBackward = hJet60EtaJeuDown2DefDist.at(2); std::vector< TH1D* > hJet60JeuDownToDefEtaFBRatio = hJet60EtaJeuDown2DefDist.at(3);
    std::vector< TH1D* > hJet60JeuDownToDefEtaBFRatio = hJet60EtaJeuDown2DefDist.at(4);

    std::vector< TH1D* > hJet80JeuUpToDefEta = hJet80EtaJeuUp2DefDist.at(0); std::vector< TH1D* > hJet80JeuUpToDefEtaForward = hJet80EtaJeuUp2DefDist.at(1);
    std::vector< TH1D* > hJet80JeuUpToDefEtaBackward = hJet80EtaJeuUp2DefDist.at(2); std::vector< TH1D* > hJet80JeuUpToDefEtaFBRatio = hJet80EtaJeuUp2DefDist.at(3);
    std::vector< TH1D* > hJet80JeuUpToDefEtaBFRatio = hJet80EtaJeuUp2DefDist.at(4);

    std::vector< TH1D* > hJet80JeuDownToDefEta = hJet80EtaJeuDown2DefDist.at(0); std::vector< TH1D* > hJet80JeuDownToDefEtaForward = hJet80EtaJeuDown2DefDist.at(1);
    std::vector< TH1D* > hJet80JeuDownToDefEtaBackward = hJet80EtaJeuDown2DefDist.at(2); std::vector< TH1D* > hJet80JeuDownToDefEtaFBRatio = hJet80EtaJeuDown2DefDist.at(3);
    std::vector< TH1D* > hJet80JeuDownToDefEtaBFRatio = hJet80EtaJeuDown2DefDist.at(4);

    std::vector< TH1D* > hJet100JeuUpToDefEta = hJet100EtaJeuUp2DefDist.at(0); std::vector< TH1D* > hJet100JeuUpToDefEtaForward = hJet100EtaJeuUp2DefDist.at(1);
    std::vector< TH1D* > hJet100JeuUpToDefEtaBackward = hJet100EtaJeuUp2DefDist.at(2); std::vector< TH1D* > hJet100JeuUpToDefEtaFBRatio = hJet100EtaJeuUp2DefDist.at(3);
    std::vector< TH1D* > hJet100JeuUpToDefEtaBFRatio = hJet100EtaJeuUp2DefDist.at(4);

    std::vector< TH1D* > hJet100JeuDownToDefEta = hJet100EtaJeuDown2DefDist.at(0); std::vector< TH1D* > hJet100JeuDownToDefEtaForward = hJet100EtaJeuDown2DefDist.at(1);
    std::vector< TH1D* > hJet100JeuDownToDefEtaBackward = hJet100EtaJeuDown2DefDist.at(2); std::vector< TH1D* > hJet100JeuDownToDefEtaFBRatio = hJet100EtaJeuDown2DefDist.at(3);
    std::vector< TH1D* > hJet100JeuDownToDefEtaBFRatio = hJet100EtaJeuDown2DefDist.at(4);

    Int_t etaTypeDef = 0;
    Int_t etaTypeForward = 1;
    Int_t etaTypeBackward = 2;
    Int_t etaTypeFB = 3;
    Int_t etaTypeBF = 4;

    plotUp2DownComparison(hMBEta, hJeuMBUpEta, hJeuMBDownEta, hMBJeuUpToDefEta, hMBJeuDownToDefEta, 0, etaTypeDef, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet60Eta, hJeuJet60UpEta, hJeuJet60DownEta, hJet60JeuUpToDefEta, hJet60JeuDownToDefEta, 0, etaTypeDef, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet80Eta, hJeuJet80UpEta, hJeuJet80DownEta, hJet80JeuUpToDefEta, hJet80JeuDownToDefEta, 0, etaTypeDef, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet100Eta, hJeuJet100UpEta, hJeuJet100DownEta, hJet100JeuUpToDefEta, hJet100JeuDownToDefEta, 0, etaTypeDef, 4, date, isCM, drawFits);

    plotUp2DownComparison(hMBEtaForward, hJeuMBUpEtaForward, hJeuMBDownEtaForward, hMBJeuUpToDefEtaForward, hMBJeuDownToDefEtaForward, 0, etaTypeForward, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet60EtaForward, hJeuJet60UpEtaForward, hJeuJet60DownEtaForward, hJet60JeuUpToDefEtaForward, hJet60JeuDownToDefEtaForward, 0, etaTypeForward, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet80EtaForward, hJeuJet80UpEtaForward, hJeuJet80DownEtaForward, hJet80JeuUpToDefEtaForward, hJet80JeuDownToDefEtaForward, 0, etaTypeForward, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet100EtaForward, hJeuJet100UpEtaForward, hJeuJet100DownEtaForward, hJet100JeuUpToDefEtaForward, hJet100JeuDownToDefEtaForward, 0, etaTypeForward, 4, date, isCM, drawFits);

    plotUp2DownComparison(hMBEtaBackward, hJeuMBUpEtaBackward, hJeuMBDownEtaBackward, hMBJeuUpToDefEtaBackward, hMBJeuDownToDefEtaBackward, 0, etaTypeBackward, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet60EtaBackward, hJeuJet60UpEtaBackward, hJeuJet60DownEtaBackward, hJet60JeuUpToDefEtaBackward, hJet60JeuDownToDefEtaBackward, 0, etaTypeBackward, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet80EtaBackward, hJeuJet80UpEtaBackward, hJeuJet80DownEtaBackward, hJet80JeuUpToDefEtaBackward, hJet80JeuDownToDefEtaBackward, 0, etaTypeBackward, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet100EtaBackward, hJeuJet100UpEtaBackward, hJeuJet100DownEtaBackward, hJet100JeuUpToDefEtaBackward, hJet100JeuDownToDefEtaBackward, 0, etaTypeBackward, 4, date, isCM, drawFits);

    plotUp2DownComparison(hMBEtaFBRatio, hJeuMBUpEtaFBRatio, hJeuMBDownEtaFBRatio, hMBJeuUpToDefEtaFBRatio, hMBJeuDownToDefEtaFBRatio, 0, etaTypeFB, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet60EtaFBRatio, hJeuJet60UpEtaFBRatio, hJeuJet60DownEtaFBRatio, hJet60JeuUpToDefEtaFBRatio, hJet60JeuDownToDefEtaFBRatio, 0, etaTypeFB, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet80EtaFBRatio, hJeuJet80UpEtaFBRatio, hJeuJet80DownEtaFBRatio, hJet80JeuUpToDefEtaFBRatio, hJet80JeuDownToDefEtaFBRatio, 0, etaTypeFB, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet100EtaFBRatio, hJeuJet100UpEtaFBRatio, hJeuJet100DownEtaFBRatio, hJet100JeuUpToDefEtaFBRatio, hJet100JeuDownToDefEtaFBRatio, 0, etaTypeFB, 4, date, isCM, drawFits);

    plotUp2DownComparison(hMBEtaBFRatio, hJeuMBUpEtaBFRatio, hJeuMBDownEtaBFRatio, hMBJeuUpToDefEtaBFRatio, hMBJeuDownToDefEtaBFRatio, 0, etaTypeBF, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet60EtaBFRatio, hJeuJet60UpEtaBFRatio, hJeuJet60DownEtaBFRatio, hJet60JeuUpToDefEtaBFRatio, hJet60JeuDownToDefEtaBFRatio, 0, etaTypeBF, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet80EtaBFRatio, hJeuJet80UpEtaBFRatio, hJeuJet80DownEtaBFRatio, hJet80JeuUpToDefEtaBFRatio, hJet80JeuDownToDefEtaBFRatio, 0, etaTypeBF, 4, date, isCM, drawFits);
    plotUp2DownComparison(hJet100EtaBFRatio, hJeuJet100UpEtaBFRatio, hJeuJet100DownEtaBFRatio, hJet100JeuUpToDefEtaBFRatio, hJet100JeuDownToDefEtaBFRatio, 0, etaTypeBF, 4, date, isCM, drawFits);
}

//________________
void plotRelSystUncrt(std::vector<TH1D*> hJeu, std::vector<TH1D*> hJer, std::vector<TH1D*> hPileup, 
                      std::vector<TH1D*> hPointingRes, std::vector<TH1D*> hTotal, 
                      Int_t etaType, TString date, Bool_t isCM) {

    // Frame name
    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Direction name
    TString dirName;
    // 0 - symmetric, 1 - forward, backward, fb, bf
    if (etaType == 0) {
        dirName = "";
    }
    else if (etaType == 1) {
        dirName = "forward";
    }
    else if (etaType == 2) {
        dirName = "backward";
    }
    else if (etaType == 3) {
        dirName = "fb";
    }
    else if (etaType == 4) {
        dirName = "bf";
    }

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // fbType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0., 30.};
    Double_t legX[2] = {0.35, 0.75};
    Double_t legY[2] = {0.65, 0.85};

    if ( etaType != 0 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
        yRange[0] = 0.0;
        yRange[1] = 30.;
    }

    // Create canvas
    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *cSyst = new TCanvas(Form("%s_%s_relSystUncrt", returnTrigName( hTotal.at(0)->GetName() ).Data(), dirName.Data() ), 
                                 Form("%s_%s_relSystUncrt", returnTrigName( hTotal.at(0)->GetName() ).Data(), dirName.Data() ), 
                                 sizeX, sizeY);
    cSyst->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    std::vector< TH1D* > hTotalInPercent{};
    std::vector< TH1D* > hJeuInPercent{};
    std::vector< TH1D* > hJerInPercent{};
    std::vector< TH1D* > hPileupInPercent{};
    std::vector< TH1D* > hPointingResInPercent{};

    // Loop over pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        // Switch to pad
        cSyst->cd(i+1);
        if ( etaType == 0 ) {
            setPadStyle();
        }
        else {
            setPadFBStyle();
        }

        hTotalInPercent.push_back( dynamic_cast<TH1D*>( hTotal.at(i)->Clone( Form("%s_inPerc", hTotal.at(i)->GetName() ) ) ) );
        hJeuInPercent.push_back( dynamic_cast<TH1D*>( hJeu.at(i)->Clone( Form("%s_inPerc", hJeu.at(i)->GetName() ) ) ) );
        hJerInPercent.push_back( dynamic_cast<TH1D*>( hJer.at(i)->Clone( Form("%s_inPerc", hJer.at(i)->GetName() ) ) ) );
        hPileupInPercent.push_back( dynamic_cast<TH1D*>( hPileup.at(i)->Clone( Form("%s_inPerc", hPileup.at(i)->GetName() ) ) ) );
        if ( !hPointingRes.empty() ) {
            hPointingResInPercent.push_back( dynamic_cast<TH1D*>( hPointingRes.at(i)->Clone( Form("%s_inPerc", hPointingRes.at(i)->GetName() ) ) ) );
        }

        set1DStyle( hTotalInPercent.at(i), 2 );
        set1DStyle( hJeuInPercent.at(i), 0 );
        set1DStyle( hJerInPercent.at(i), 1 );
        set1DStyle( hPileupInPercent.at(i), 5 );
        if ( !hPointingRes.empty() ) {
            set1DStyle( hPointingResInPercent.at(i), 6 );
        }

        hTotalInPercent.at(i)->Scale(100.0);
        hJeuInPercent.at(i)->Scale(100.0);
        hJerInPercent.at(i)->Scale(100.0);
        hPileupInPercent.at(i)->Scale(100.0);
        if ( !hPointingRes.empty() ) {
            hPointingResInPercent.at(i)->Scale(100.0);
        }

        // Draw histograms
        hTotalInPercent.at(i)->Draw();
        hJeuInPercent.at(i)->Draw("same");
        hJerInPercent.at(i)->Draw("same");
        hPileupInPercent.at(i)->Draw("same");
        if ( !hPointingResInPercent.empty() ) {
            hPointingRes.at(i)->Draw("same");
        }
        hTotalInPercent.at(i)->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        hTotalInPercent.at(i)->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        hTotalInPercent.at(i)->GetYaxis()->SetTitle("Rel. syst. uncrt. [%]");

        // Draw text
        t.DrawLatexNDC(0.35, 0.94, Form("%d < p_{T}^{jet} < %d GeV", ptDijetLow.at(i), ptDijetHi.at(i) ) );

        // Draw legend
        TLegend *leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hTotalInPercent.at(i), "Total", "l");
        leg->AddEntry(hJeuInPercent.at(i), "JES", "l");
        leg->AddEntry(hJerInPercent.at(i), "JER", "l");
        leg->AddEntry(hPileupInPercent.at(i), "Pileup", "l");
        if ( !hPointingResInPercent.empty() ) {
            leg->AddEntry(hPointingResInPercent.at(i), "Pointing res.", "l");
        }
        leg->Draw();
    } // for (Int_t i{0}; i<ptBins; i++)

    cSyst->SaveAs( Form("%s/systematics/%s_%s_relSystUncrt_%s.pdf", date.Data(), returnTrigName( hTotal.at(0)->GetName() ).Data(), dirName.Data(), frame.Data() ) );
}

//________________
void plotRelativeSystematicUncertainties(std::vector< std::vector<TH1D*> > hMBJeuRelSystDist, std::vector< std::vector<TH1D*> > hJet60JeuRelSystDist, 
                                         std::vector< std::vector<TH1D*> > hJet80JeuRelSystDist, std::vector< std::vector<TH1D*> > hJet100JeuRelSystDist,
                                         std::vector< std::vector<TH1D*> > hEmbeddingJerRelSystDist, std::vector< std::vector<TH1D*> > hEmbeddingPointResRelSystDist,
                                         std::vector< std::vector<TH1D*> > hMBPileupRelSystDist, std::vector< std::vector<TH1D*> > hJet60PileupRelSystDist, 
                                         std::vector< std::vector<TH1D*> > hJet80PileupRelSystDist, std::vector< std::vector<TH1D*> > hJet100PileupRelSystDist,
                                         std::vector< std::vector<TH1D*> > hMBTotalRelSystDist, std::vector< std::vector<TH1D*> > hJet60TotalRelSystDist, 
                                         std::vector< std::vector<TH1D*> > hJet80TotalRelSystDist, std::vector< std::vector<TH1D*> > hJet100TotalRelSystDist, 
                                         TString date, Bool_t isCM) {
    
    // Retrieve individual systematic uncertainties

    // JES
    std::vector< TH1D* > hMBEtaJeuRelSyst = hMBJeuRelSystDist.at(0); std::vector< TH1D* > hMBEtaForwardJeuRelSyst = hMBJeuRelSystDist.at(1);
    std::vector< TH1D* > hMBEtaBackwardJeuRelSyst = hMBJeuRelSystDist.at(2); std::vector< TH1D* > hMBEtaFBJeuRelSyst = hMBJeuRelSystDist.at(3);
    std::vector< TH1D* > hMBEtaBFJeuRelSyst = hMBJeuRelSystDist.at(4);

    std::vector< TH1D* > hJet60EtaJeuRelSyst = hJet60JeuRelSystDist.at(0); std::vector< TH1D* > hJet60EtaForwardJeuRelSyst = hJet60JeuRelSystDist.at(1);
    std::vector< TH1D* > hJet60EtaBackwardJeuRelSyst = hJet60JeuRelSystDist.at(2); std::vector< TH1D* > hJet60EtaFBJeuRelSyst = hJet60JeuRelSystDist.at(3);
    std::vector< TH1D* > hJet60EtaBFJeuRelSyst = hJet60JeuRelSystDist.at(4);

    std::vector< TH1D* > hJet80EtaJeuRelSyst = hJet80JeuRelSystDist.at(0); std::vector< TH1D* > hJet80EtaForwardJeuRelSyst = hJet80JeuRelSystDist.at(1);
    std::vector< TH1D* > hJet80EtaBackwardJeuRelSyst = hJet80JeuRelSystDist.at(2); std::vector< TH1D* > hJet80EtaFBJeuRelSyst = hJet80JeuRelSystDist.at(3);
    std::vector< TH1D* > hJet80EtaBFJeuRelSyst = hJet80JeuRelSystDist.at(4);

    std::vector< TH1D* > hJet100EtaJeuRelSyst = hJet100JeuRelSystDist.at(0); std::vector< TH1D* > hJet100EtaForwardJeuRelSyst = hJet100JeuRelSystDist.at(1);
    std::vector< TH1D* > hJet100EtaBackwardJeuRelSyst = hJet100JeuRelSystDist.at(2); std::vector< TH1D* > hJet100EtaFBJeuRelSyst = hJet100JeuRelSystDist.at(3);
    std::vector< TH1D* > hJet100EtaBFJeuRelSyst = hJet100JeuRelSystDist.at(4);

    // JER
    std::vector< TH1D* > hEmbeddingEtaJerRelSyst = hEmbeddingJerRelSystDist.at(0); std::vector< TH1D* > hEmbeddingEtaForwardJerRelSyst = hEmbeddingJerRelSystDist.at(1);
    std::vector< TH1D* > hEmbeddingEtaBackwardJerRelSyst = hEmbeddingJerRelSystDist.at(2); std::vector< TH1D* > hEmbeddingEtaFBJerRelSyst = hEmbeddingJerRelSystDist.at(3);
    std::vector< TH1D* > hEmbeddingEtaBFJerRelSyst = hEmbeddingJerRelSystDist.at(4);

    // Pointing resolution
    std::vector< TH1D* > hEmbeddingEtaPointingResRelSyst = hEmbeddingPointResRelSystDist.at(0); std::vector< TH1D* > hEmbeddingEtaForwardPointingResRelSyst = hEmbeddingPointResRelSystDist.at(1);
    std::vector< TH1D* > hEmbeddingEtaBackwardPointingResRelSyst = hEmbeddingPointResRelSystDist.at(2); std::vector< TH1D* > hEmbeddingEtaFBPointingResRelSyst = hEmbeddingPointResRelSystDist.at(3);
    std::vector< TH1D* > hEmbeddingEtaBFPointingResRelSyst = hEmbeddingPointResRelSystDist.at(4);

    // Pileup
    std::vector< TH1D* > hMBEtaPileupRelSyst = hMBPileupRelSystDist.at(0); std::vector< TH1D* > hMBEtaForwardPileupRelSyst = hMBPileupRelSystDist.at(1);
    std::vector< TH1D* > hMBEtaBackwardPileupRelSyst = hMBPileupRelSystDist.at(2); std::vector< TH1D* > hMBEtaFBPileupRelSyst = hMBPileupRelSystDist.at(3);
    std::vector< TH1D* > hMBEtaBFPileupRelSyst = hMBPileupRelSystDist.at(4);

    std::vector< TH1D* > hJet60EtaPileupRelSyst = hJet60PileupRelSystDist.at(0); std::vector< TH1D* > hJet60EtaForwardPileupRelSyst = hJet60PileupRelSystDist.at(1);
    std::vector< TH1D* > hJet60EtaBackwardPileupRelSyst = hJet60PileupRelSystDist.at(2); std::vector< TH1D* > hJet60EtaFBPileupRelSyst = hJet60PileupRelSystDist.at(3);
    std::vector< TH1D* > hJet60EtaBFPileupRelSyst = hJet60PileupRelSystDist.at(4);

    std::vector< TH1D* > hJet80EtaPileupRelSyst = hJet80PileupRelSystDist.at(0); std::vector< TH1D* > hJet80EtaForwardPileupRelSyst = hJet80PileupRelSystDist.at(1);
    std::vector< TH1D* > hJet80EtaBackwardPileupRelSyst = hJet80PileupRelSystDist.at(2); std::vector< TH1D* > hJet80EtaFBPileupRelSyst = hJet80PileupRelSystDist.at(3);
    std::vector< TH1D* > hJet80EtaBFPileupRelSyst = hJet80PileupRelSystDist.at(4);

    std::vector< TH1D* > hJet100EtaPileupRelSyst = hJet100PileupRelSystDist.at(0); std::vector< TH1D* > hJet100EtaForwardPileupRelSyst = hJet100PileupRelSystDist.at(1);
    std::vector< TH1D* > hJet100EtaBackwardPileupRelSyst = hJet100PileupRelSystDist.at(2); std::vector< TH1D* > hJet100EtaFBPileupRelSyst = hJet100PileupRelSystDist.at(3);
    std::vector< TH1D* > hJet100EtaBFPileupRelSyst = hJet100PileupRelSystDist.at(4);

    // Total
    std::vector< TH1D* > hMBTotalEtaRelSyst = hMBTotalRelSystDist.at(0); std::vector< TH1D* > hMBTotalEtaForwardRelSyst = hMBTotalRelSystDist.at(1);
    std::vector< TH1D* > hMBTotalEtaBackwardRelSyst = hMBTotalRelSystDist.at(2); std::vector< TH1D* > hMBTotalEtaFBRelSyst = hMBTotalRelSystDist.at(3);
    std::vector< TH1D* > hMBTotalEtaBFRelSyst = hMBTotalRelSystDist.at(4);

    std::vector< TH1D* > hJet60TotalEtaRelSyst = hJet60TotalRelSystDist.at(0); std::vector< TH1D* > hJet60TotalEtaForwardRelSyst = hJet60TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet60TotalEtaBackwardRelSyst = hJet60TotalRelSystDist.at(2); std::vector< TH1D* > hJet60TotalEtaFBRelSyst = hJet60TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet60TotalEtaBFRelSyst = hJet60TotalRelSystDist.at(4);

    std::vector< TH1D* > hJet80TotalEtaRelSyst = hJet80TotalRelSystDist.at(0); std::vector< TH1D* > hJet80TotalEtaForwardRelSyst = hJet80TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet80TotalEtaBackwardRelSyst = hJet80TotalRelSystDist.at(2); std::vector< TH1D* > hJet80TotalEtaFBRelSyst = hJet80TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet80TotalEtaBFRelSyst = hJet80TotalRelSystDist.at(4);

    std::vector< TH1D* > hJet100TotalEtaRelSyst = hJet100TotalRelSystDist.at(0); std::vector< TH1D* > hJet100TotalEtaForwardRelSyst = hJet100TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet100TotalEtaBackwardRelSyst = hJet100TotalRelSystDist.at(2); std::vector< TH1D* > hJet100TotalEtaFBRelSyst = hJet100TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet100TotalEtaBFRelSyst = hJet100TotalRelSystDist.at(4);


    Int_t etaTypeDef = 0;
    Int_t etaTypeForward = 1;
    Int_t etaTypeBackward = 2;
    Int_t etaTypeFB = 3;
    Int_t etaTypeBF = 4;

    // Plot relative systematic uncertainties
    plotRelSystUncrt(hMBEtaJeuRelSyst, hEmbeddingEtaJerRelSyst, hMBEtaPileupRelSyst, hEmbeddingEtaPointingResRelSyst, hMBTotalEtaRelSyst, etaTypeDef, date, isCM);
    plotRelSystUncrt(hJet60EtaJeuRelSyst, hEmbeddingEtaJerRelSyst, hJet60EtaPileupRelSyst, hEmbeddingEtaPointingResRelSyst, hJet60TotalEtaRelSyst, etaTypeDef, date, isCM);
    plotRelSystUncrt(hJet80EtaJeuRelSyst, hEmbeddingEtaJerRelSyst, hJet80EtaPileupRelSyst, hEmbeddingEtaPointingResRelSyst, hJet80TotalEtaRelSyst, etaTypeDef, date, isCM);
    plotRelSystUncrt(hJet100EtaJeuRelSyst, hEmbeddingEtaJerRelSyst, hJet100EtaPileupRelSyst, hEmbeddingEtaPointingResRelSyst, hJet100TotalEtaRelSyst, etaTypeDef, date, isCM);

}

//________________
void plotJER(std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hJerUpEtaDist, 
             std::vector< std::vector<TH1D*> > hJerDownEtaDist, 
             std::vector< std::vector<TH1D*> > hEmbeddingEtaJerUpToDefDist, std::vector< std::vector<TH1D*> > hEmbeddingEtaJerDownToDefDist,
             TString date, Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Data

    std::vector< TH1D* > hEmbeddingEta = hEmbeddingEtaDist.at(0); std::vector< TH1D* > hEmbeddingEtaForward = hEmbeddingEtaDist.at(1); 
    std::vector< TH1D* > hEmbeddingEtaBackward = hEmbeddingEtaDist.at(2); std::vector< TH1D* > hEmbeddingEtaFBRatio = hEmbeddingEtaDist.at(3);
    std::vector< TH1D* > hEmbeddingEtaBFRatio = hEmbeddingEtaDist.at(4);

    std::vector< TH1D* > hJerUpEta = hJerUpEtaDist.at(0); std::vector< TH1D* > hJerUpEtaForward = hJerUpEtaDist.at(1);
    std::vector< TH1D* > hJerUpEtaBackward = hJerUpEtaDist.at(2); std::vector< TH1D* > hJerUpEtaFBRatio = hJerUpEtaDist.at(3);
    std::vector< TH1D* > hJerUpEtaBFRatio = hJerUpEtaDist.at(4);

    std::vector< TH1D* > hJerDownEta = hJerDownEtaDist.at(0); std::vector< TH1D* > hJerDownEtaForward = hJerDownEtaDist.at(1);
    std::vector< TH1D* > hJerDownEtaBackward = hJerDownEtaDist.at(2); std::vector< TH1D* > hJerDownEtaFBRatio = hJerDownEtaDist.at(3);
    std::vector< TH1D* > hJerDownEtaBFRatio = hJerDownEtaDist.at(4);

    std::vector< TH1D* > hEmbeddingJerUpToDefEta = hEmbeddingEtaJerUpToDefDist.at(0); std::vector< TH1D* > hEmbeddingJerUpToDefEtaForward = hEmbeddingEtaJerUpToDefDist.at(1);
    std::vector< TH1D* > hEmbeddingJerUpToDefEtaBackward = hEmbeddingEtaJerUpToDefDist.at(2); std::vector< TH1D* > hEmbeddingJerUpToDefEtaFBRatio = hEmbeddingEtaJerUpToDefDist.at(3);
    std::vector< TH1D* > hEmbeddingJerUpToDefEtaBFRatio = hEmbeddingEtaJerUpToDefDist.at(4);

    std::vector< TH1D* > hEmbeddingJerDownToDefEta = hEmbeddingEtaJerDownToDefDist.at(0); std::vector< TH1D* > hEmbeddingJerDownToDefEtaForward = hEmbeddingEtaJerDownToDefDist.at(1);
    std::vector< TH1D* > hEmbeddingJerDownToDefEtaBackward = hEmbeddingEtaJerDownToDefDist.at(2); std::vector< TH1D* > hEmbeddingJerDownToDefEtaFBRatio = hEmbeddingEtaJerDownToDefDist.at(3);
    std::vector< TH1D* > hEmbeddingJerDownToDefEtaBFRatio = hEmbeddingEtaJerDownToDefDist.at(4);

    Int_t etaTypeDef = 0;
    Int_t etaTypeForward = 1;
    Int_t etaTypeBackward = 2;
    Int_t etaTypeFB = 3;
    Int_t etaTypeBF = 4;

    plotUp2DownComparison(hEmbeddingEta, hJerUpEta, hJerDownEta, hEmbeddingJerUpToDefEta, hEmbeddingJerDownToDefEta, 1, etaTypeDef, 2, date, frame, drawFits);
    plotUp2DownComparison(hEmbeddingEtaForward, hJerUpEtaForward, hJerDownEtaForward, hEmbeddingJerUpToDefEtaForward, hEmbeddingJerDownToDefEtaForward, 1, etaTypeForward, 2, date, frame, drawFits);
    plotUp2DownComparison(hEmbeddingEtaBackward, hJerUpEtaBackward, hJerDownEtaBackward, hEmbeddingJerUpToDefEtaBackward, hEmbeddingJerDownToDefEtaBackward, 1, etaTypeBackward, 2, date, frame, drawFits);
    plotUp2DownComparison(hEmbeddingEtaFBRatio, hJerUpEtaFBRatio, hJerDownEtaFBRatio, hEmbeddingJerUpToDefEtaFBRatio, hEmbeddingJerDownToDefEtaFBRatio, 1, etaTypeFB, 2, date, frame, drawFits);
    plotUp2DownComparison(hEmbeddingEtaBFRatio, hJerUpEtaBFRatio, hJerDownEtaBFRatio, hEmbeddingJerUpToDefEtaBFRatio, hEmbeddingJerDownToDefEtaBFRatio, 1, etaTypeBF, 2, date, frame, drawFits);
}

//________________
void plotPileup(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hMBGplusEtaDist, 
                std::vector< std::vector<TH1D*> > hMBVtx1EtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist, 
                std::vector< std::vector<TH1D*> > hJet60GplusEtaDist, std::vector< std::vector<TH1D*> > hJet60Vtx1EtaDist, 
                std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet80GplusEtaDist, 
                std::vector< std::vector<TH1D*> > hJet80Vtx1EtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist, 
                std::vector< std::vector<TH1D*> > hJet100GplusEtaDist, std::vector< std::vector<TH1D*> > hJet100Vtx1EtaDist,
                std::vector< std::vector<TH1D*> > hMBEtaGplus2DefDist, std::vector< std::vector<TH1D*> > hMBEtaVtxOne2DefDist,
                std::vector< std::vector<TH1D*> > hJet60EtaGplus2DefDist, std::vector< std::vector<TH1D*> > hJet60EtaVtxOne2DefDist,
                std::vector< std::vector<TH1D*> > hJet80EtaGplus2DefDist, std::vector< std::vector<TH1D*> > hJet80EtaVtxOne2DefDist,
                std::vector< std::vector<TH1D*> > hJet100EtaGplus2DefDist, std::vector< std::vector<TH1D*> > hJet100EtaVtxOne2DefDist,
                TString date, Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Data

    std::vector< TH1D* > hMBEta = hMBEtaDist.at(0); std::vector< TH1D* > hMBEtaForward = hMBEtaDist.at(1); 
    std::vector< TH1D* > hMBEtaBackward = hMBEtaDist.at(2); std::vector< TH1D* > hMBEtaFBRatio = hMBEtaDist.at(3);
    std::vector< TH1D* > hMBEtaBFRatio = hMBEtaDist.at(4);

    std::vector< TH1D* > hMBGplusEta = hMBGplusEtaDist.at(0); std::vector< TH1D* > hMBGplusEtaForward = hMBGplusEtaDist.at(1);
    std::vector< TH1D* > hMBGplusEtaBackward = hMBGplusEtaDist.at(2); std::vector< TH1D* > hMBGplusEtaFBRatio = hMBGplusEtaDist.at(3);
    std::vector< TH1D* > hMBGplusEtaBFRatio = hMBGplusEtaDist.at(4);

    std::vector< TH1D* > hMBVtx1Eta = {};
    std::vector< TH1D* > hMBVtx1EtaForward = {};
    std::vector< TH1D* > hMBVtx1EtaBackward = {};
    std::vector< TH1D* > hMBVtx1EtaFBRatio = {};
    std::vector< TH1D* > hMBVtx1EtaBFRatio = {};
    
    // std::vector< TH1D* > hMBVtx1Eta = hMBVtx1EtaDist.at(0); std::vector< TH1D* > hMBVtx1EtaForward = hMBVtx1EtaDist.at(1);
    // std::vector< TH1D* > hMBVtx1EtaBackward = hMBVtx1EtaDist.at(2); std::vector< TH1D* > hMBVtx1EtaFBRatio = hMBVtx1EtaDist.at(3);
    // std::vector< TH1D* > hMBVtx1EtaBFRatio = hMBVtx1EtaDist.at(4);

    std::vector< TH1D* > hJet60Eta = hJet60EtaDist.at(0); std::vector< TH1D* > hJet60EtaForward = hJet60EtaDist.at(1);
    std::vector< TH1D* > hJet60EtaBackward = hJet60EtaDist.at(2); std::vector< TH1D* > hJet60EtaFBRatio = hJet60EtaDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatio = hJet60EtaDist.at(4);

    std::vector< TH1D* > hJet60GplusEta = hJet60GplusEtaDist.at(0); std::vector< TH1D* > hJet60GplusEtaForward = hJet60GplusEtaDist.at(1);
    std::vector< TH1D* > hJet60GplusEtaBackward = hJet60GplusEtaDist.at(2); std::vector< TH1D* > hJet60GplusEtaFBRatio = hJet60GplusEtaDist.at(3);
    std::vector< TH1D* > hJet60GplusEtaBFRatio = hJet60GplusEtaDist.at(4);

    std::vector< TH1D* > hJet60Vtx1Eta = {};
    std::vector< TH1D* > hJet60Vtx1EtaForward = {};
    std::vector< TH1D* > hJet60Vtx1EtaBackward = {};
    std::vector< TH1D* > hJet60Vtx1EtaFBRatio = {};
    std::vector< TH1D* > hJet60Vtx1EtaBFRatio = {};

    // std::vector< TH1D* > hJet60Vtx1Eta = hJet60Vtx1EtaDist.at(0); std::vector< TH1D* > hJet60Vtx1EtaForward = hJet60Vtx1EtaDist.at(1);
    // std::vector< TH1D* > hJet60Vtx1EtaBackward = hJet60Vtx1EtaDist.at(2); std::vector< TH1D* > hJet60Vtx1EtaFBRatio = hJet60Vtx1EtaDist.at(3);
    // std::vector< TH1D* > hJet60Vtx1EtaBFRatio = hJet60Vtx1EtaDist.at(4);
    
    std::vector< TH1D* > hJet80Eta = hJet80EtaDist.at(0); std::vector< TH1D* > hJet80EtaForward = hJet80EtaDist.at(1);
    std::vector< TH1D* > hJet80EtaBackward = hJet80EtaDist.at(2); std::vector< TH1D* > hJet80EtaFBRatio = hJet80EtaDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatio = hJet80EtaDist.at(4);

    std::vector< TH1D* > hJet80GplusEta = hJet80GplusEtaDist.at(0); std::vector< TH1D* > hJet80GplusEtaForward = hJet80GplusEtaDist.at(1);
    std::vector< TH1D* > hJet80GplusEtaBackward = hJet80GplusEtaDist.at(2); std::vector< TH1D* > hJet80GplusEtaFBRatio = hJet80GplusEtaDist.at(3);
    std::vector< TH1D* > hJet80GplusEtaBFRatio = hJet80GplusEtaDist.at(4);

    std::vector< TH1D* > hJet80Vtx1Eta = {};
    std::vector< TH1D* > hJet80Vtx1EtaForward = {};
    std::vector< TH1D* > hJet80Vtx1EtaBackward = {};
    std::vector< TH1D* > hJet80Vtx1EtaFBRatio = {};
    std::vector< TH1D* > hJet80Vtx1EtaBFRatio = {};

    // std::vector< TH1D* > hJet80Vtx1Eta = hJet80Vtx1EtaDist.at(0); std::vector< TH1D* > hJet80Vtx1EtaForward = hJet80Vtx1EtaDist.at(1);
    // std::vector< TH1D* > hJet80Vtx1EtaBackward = hJet80Vtx1EtaDist.at(2); std::vector< TH1D* > hJet80Vtx1EtaFBRatio = hJet80Vtx1EtaDist.at(3);
    // std::vector< TH1D* > hJet80Vtx1EtaBFRatio = hJet80Vtx1EtaDist.at(4);

    std::vector< TH1D* > hJet100Eta = hJet100EtaDist.at(0); std::vector< TH1D* > hJet100EtaForward = hJet100EtaDist.at(1);
    std::vector< TH1D* > hJet100EtaBackward = hJet100EtaDist.at(2); std::vector< TH1D* > hJet100EtaFBRatio = hJet100EtaDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatio = hJet100EtaDist.at(4);

    std::vector< TH1D* > hJet100GplusEta = hJet100GplusEtaDist.at(0); std::vector< TH1D* > hJet100GplusEtaForward = hJet100GplusEtaDist.at(1);
    std::vector< TH1D* > hJet100GplusEtaBackward = hJet100GplusEtaDist.at(2); std::vector< TH1D* > hJet100GplusEtaFBRatio = hJet100GplusEtaDist.at(3);
    std::vector< TH1D* > hJet100GplusEtaBFRatio = hJet100GplusEtaDist.at(4);

    std::vector< TH1D* > hJet100Vtx1Eta = {};
    std::vector< TH1D* > hJet100Vtx1EtaForward = {};
    std::vector< TH1D* > hJet100Vtx1EtaBackward = {};
    std::vector< TH1D* > hJet100Vtx1EtaFBRatio = {};
    std::vector< TH1D* > hJet100Vtx1EtaBFRatio = {};

    // std::vector< TH1D* > hJet100Vtx1Eta = hJet100Vtx1EtaDist.at(0); std::vector< TH1D* > hJet100Vtx1EtaForward = hJet100Vtx1EtaDist.at(1);
    // std::vector< TH1D* > hJet100Vtx1EtaBackward = hJet100Vtx1EtaDist.at(2); std::vector< TH1D* > hJet100Vtx1EtaFBRatio = hJet100Vtx1EtaDist.at(3);
    // std::vector< TH1D* > hJet100Vtx1EtaBFRatio = hJet100Vtx1EtaDist.at(4);

    // Ratios
    std::vector< TH1D* > hMBGplusToDefEta = hMBEtaGplus2DefDist.at(0); std::vector< TH1D* > hMBGplusToDefEtaForward = hMBEtaGplus2DefDist.at(1);
    std::vector< TH1D* > hMBGplusToDefEtaBackward = hMBEtaGplus2DefDist.at(2); std::vector< TH1D* > hMBGplusToDefEtaFBRatio = hMBEtaGplus2DefDist.at(3);
    std::vector< TH1D* > hMBGplusToDefEtaBFRatio = hMBEtaGplus2DefDist.at(4);

    std::vector< TH1D* > hMBVtx1ToDefEta = {}; std::vector< TH1D* > hMBVtx1ToDefEtaForward = {};
    std::vector< TH1D* > hMBVtx1ToDefEtaBackward = {}; std::vector< TH1D* > hMBVtx1ToDefEtaFBRatio = {};
    std::vector< TH1D* > hMBVtx1ToDefEtaBFRatio = {};

    // std::vector< TH1D* > hMBVtx1ToDefEta = hMBEtaVtxOne2DefDist.at(0); std::vector< TH1D* > hMBVtx1ToDefEtaForward = hMBEtaVtxOne2DefDist.at(1);
    // std::vector< TH1D* > hMBVtx1ToDefEtaBackward = hMBEtaVtxOne2DefDist.at(2); std::vector< TH1D* > hMBVtx1ToDefEtaFBRatio = hMBEtaVtxOne2DefDist.at(3);
    // std::vector< TH1D* > hMBVtx1ToDefEtaBFRatio = hMBEtaVtxOne2DefDist.at(4);

    std::vector< TH1D* > hJet60GplusToDefEta = hJet60EtaGplus2DefDist.at(0); std::vector< TH1D* > hJet60GplusToDefEtaForward = hJet60EtaGplus2DefDist.at(1);
    std::vector< TH1D* > hJet60GplusToDefEtaBackward = hJet60EtaGplus2DefDist.at(2); std::vector< TH1D* > hJet60GplusToDefEtaFBRatio = hJet60EtaGplus2DefDist.at(3);
    std::vector< TH1D* > hJet60GplusToDefEtaBFRatio = hJet60EtaGplus2DefDist.at(4);

    std::vector< TH1D* > hJet60Vtx1ToDefEta = {}; std::vector< TH1D* > hJet60Vtx1ToDefEtaForward = {};
    std::vector< TH1D* > hJet60Vtx1ToDefEtaBackward = {}; std::vector< TH1D* > hJet60Vtx1ToDefEtaFBRatio = {};
    std::vector< TH1D* > hJet60Vtx1ToDefEtaBFRatio = {};

    // std::vector< TH1D* > hJet60Vtx1ToDefEta = hJet60EtaVtxOne2DefDist.at(0); std::vector< TH1D* > hJet60Vtx1ToDefEtaForward = hJet60EtaVtxOne2DefDist.at(1);
    // std::vector< TH1D* > hJet60Vtx1ToDefEtaBackward = hJet60EtaVtxOne2DefDist.at(2); std::vector< TH1D* > hJet60Vtx1ToDefEtaFBRatio = hJet60EtaVtxOne2DefDist.at(3);
    // std::vector< TH1D* > hJet60Vtx1ToDefEtaBFRatio = hJet60EtaVtxOne2DefDist.at(4);

    std::vector< TH1D* > hJet80GplusToDefEta = hJet80EtaGplus2DefDist.at(0); std::vector< TH1D* > hJet80GplusToDefEtaForward = hJet80EtaGplus2DefDist.at(1);
    std::vector< TH1D* > hJet80GplusToDefEtaBackward = hJet80EtaGplus2DefDist.at(2); std::vector< TH1D* > hJet80GplusToDefEtaFBRatio = hJet80EtaGplus2DefDist.at(3);
    std::vector< TH1D* > hJet80GplusToDefEtaBFRatio = hJet80EtaGplus2DefDist.at(4);

    std::vector< TH1D* > hJet80Vtx1ToDefEta = {}; std::vector< TH1D* > hJet80Vtx1ToDefEtaForward = {};
    std::vector< TH1D* > hJet80Vtx1ToDefEtaBackward = {}; std::vector< TH1D* > hJet80Vtx1ToDefEtaFBRatio = {};
    std::vector< TH1D* > hJet80Vtx1ToDefEtaBFRatio = {};

    // std::vector< TH1D* > hJet80Vtx1ToDefEta = hJet80EtaVtxOne2DefDist.at(0); std::vector< TH1D* > hJet80Vtx1ToDefEtaForward = hJet80EtaVtxOne2DefDist.at(1);
    // std::vector< TH1D* > hJet80Vtx1ToDefEtaBackward = hJet80EtaVtxOne2DefDist.at(2); std::vector< TH1D* > hJet80Vtx1ToDefEtaFBRatio = hJet80EtaVtxOne2DefDist.at(3);
    // std::vector< TH1D* > hJet80Vtx1ToDefEtaBFRatio = hJet80EtaVtxOne2DefDist.at(4);

    std::vector< TH1D* > hJet100GplusToDefEta = hJet100EtaGplus2DefDist.at(0); std::vector< TH1D* > hJet100GplusToDefEtaForward = hJet100EtaGplus2DefDist.at(1);
    std::vector< TH1D* > hJet100GplusToDefEtaBackward = hJet100EtaGplus2DefDist.at(2); std::vector< TH1D* > hJet100GplusToDefEtaFBRatio = hJet100EtaGplus2DefDist.at(3);
    std::vector< TH1D* > hJet100GplusToDefEtaBFRatio = hJet100EtaGplus2DefDist.at(4);

    std::vector< TH1D* > hJet100Vtx1ToDefEta = {}; std::vector< TH1D* > hJet100Vtx1ToDefEtaForward = {};
    std::vector< TH1D* > hJet100Vtx1ToDefEtaBackward = {}; std::vector< TH1D* > hJet100Vtx1ToDefEtaFBRatio = {};
    std::vector< TH1D* > hJet100Vtx1ToDefEtaBFRatio = {};

    // std::vector< TH1D* > hJet100Vtx1ToDefEta = hJet100EtaVtxOne2DefDist.at(0); std::vector< TH1D* > hJet100Vtx1ToDefEtaForward = hJet100EtaVtxOne2DefDist.at(1);
    // std::vector< TH1D* > hJet100Vtx1ToDefEtaBackward = hJet100EtaVtxOne2DefDist.at(2); std::vector< TH1D* > hJet100Vtx1ToDefEtaFBRatio = hJet100EtaVtxOne2DefDist.at(3);
    // std::vector< TH1D* > hJet100Vtx1ToDefEtaBFRatio = hJet100EtaVtxOne2DefDist.at(4);

    Int_t etaTypeDef = 0;
    Int_t etaTypeForward = 1;
    Int_t etaTypeBackward = 2;
    Int_t etaTypeFB = 3;
    Int_t etaTypeBF = 4;

    // Plotting
    plotUp2DownComparison(hMBEta, hMBGplusEta, hMBVtx1Eta, hMBGplusToDefEta, hMBVtx1ToDefEta, 2, etaTypeDef, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet60Eta, hJet60GplusEta, hJet60Vtx1Eta, hJet60GplusToDefEta, hJet60Vtx1ToDefEta, 2, etaTypeDef, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet80Eta, hJet80GplusEta, hJet80Vtx1Eta, hJet80GplusToDefEta, hJet80Vtx1ToDefEta, 2, etaTypeDef, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet100Eta, hJet100GplusEta, hJet100Vtx1Eta, hJet100GplusToDefEta, hJet100Vtx1ToDefEta, 2, etaTypeDef, 1, date, frame, drawFits);

    plotUp2DownComparison(hMBEtaForward, hMBGplusEtaForward, hMBVtx1EtaForward, hMBGplusToDefEtaForward, hMBVtx1ToDefEtaForward, 2, etaTypeForward, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet60EtaForward, hJet60GplusEtaForward, hJet60Vtx1EtaForward, hJet60GplusToDefEtaForward, hJet60Vtx1ToDefEtaForward, 2, etaTypeForward, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet80EtaForward, hJet80GplusEtaForward, hJet80Vtx1EtaForward, hJet80GplusToDefEtaForward, hJet80Vtx1ToDefEtaForward, 2, etaTypeForward, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet100EtaForward, hJet100GplusEtaForward, hJet100Vtx1EtaForward, hJet100GplusToDefEtaForward, hJet100Vtx1ToDefEtaForward, 2, etaTypeForward, 1, date, frame, drawFits);

    plotUp2DownComparison(hMBEtaBackward, hMBGplusEtaBackward, hMBVtx1EtaBackward, hMBGplusToDefEtaBackward, hMBVtx1ToDefEtaBackward, 2, etaTypeBackward, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet60EtaBackward, hJet60GplusEtaBackward, hJet60Vtx1EtaBackward, hJet60GplusToDefEtaBackward, hJet60Vtx1ToDefEtaBackward, 2, etaTypeBackward, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet80EtaBackward, hJet80GplusEtaBackward, hJet80Vtx1EtaBackward, hJet80GplusToDefEtaBackward, hJet80Vtx1ToDefEtaBackward, 2, etaTypeBackward, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet100EtaBackward, hJet100GplusEtaBackward, hJet100Vtx1EtaBackward, hJet100GplusToDefEtaBackward, hJet100Vtx1ToDefEtaBackward, 2, etaTypeBackward, 1, date, frame, drawFits);

    plotUp2DownComparison(hMBEtaFBRatio, hMBGplusEtaFBRatio, hMBVtx1EtaFBRatio, hMBGplusToDefEtaFBRatio, hMBVtx1ToDefEtaFBRatio, 2, etaTypeFB, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet60EtaFBRatio, hJet60GplusEtaFBRatio, hJet60Vtx1EtaFBRatio, hJet60GplusToDefEtaFBRatio, hJet60Vtx1ToDefEtaFBRatio, 2, etaTypeFB, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet80EtaFBRatio, hJet80GplusEtaFBRatio, hJet80Vtx1EtaFBRatio, hJet80GplusToDefEtaFBRatio, hJet80Vtx1ToDefEtaFBRatio, 2, etaTypeFB, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet100EtaFBRatio, hJet100GplusEtaFBRatio, hJet100Vtx1EtaFBRatio, hJet100GplusToDefEtaFBRatio, hJet100Vtx1ToDefEtaFBRatio, 2, etaTypeFB, 1, date, frame, drawFits);

    plotUp2DownComparison(hMBEtaBFRatio, hMBGplusEtaBFRatio, hMBVtx1EtaBFRatio, hMBGplusToDefEtaBFRatio, hMBVtx1ToDefEtaBFRatio, 2, etaTypeBF, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet60EtaBFRatio, hJet60GplusEtaBFRatio, hJet60Vtx1EtaBFRatio, hJet60GplusToDefEtaBFRatio, hJet60Vtx1ToDefEtaBFRatio, 2, etaTypeBF, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet80EtaBFRatio, hJet80GplusEtaBFRatio, hJet80Vtx1EtaBFRatio, hJet80GplusToDefEtaBFRatio, hJet80Vtx1ToDefEtaBFRatio, 2, etaTypeBF, 1, date, frame, drawFits);
    plotUp2DownComparison(hJet100EtaBFRatio, hJet100GplusEtaBFRatio, hJet100Vtx1EtaBFRatio, hJet100GplusToDefEtaBFRatio, hJet100Vtx1ToDefEtaBFRatio, 2, etaTypeBF, 1, date, frame, drawFits);
}

//________________
void writeVectors2RootFile(std::vector< std::vector<TH1D*> > vectorOfVectors) {
    for (Int_t i = 0; i < vectorOfVectors.size(); i++) {
        for (Int_t j = 0; j < vectorOfVectors.at(i).size(); j++) {
            vectorOfVectors.at(i).at(j)->Write();
        }
    }
}

//________________
void writeDistributions2RootFile(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hMBPbGoingEtaDist, 
                                 std::vector< std::vector<TH1D*> > hMBPGoingEtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet60PbGoingEtaDist, std::vector< std::vector<TH1D*> > hJet60PGoingEtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet80PbGoingEtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet80PGoingEtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet100PbGoingEtaDist, std::vector< std::vector<TH1D*> > hJet100PGoingEtaDist, 
                                 std::vector< std::vector<TH1D*> > hJeuMBUpEtaDist, std::vector< std::vector<TH1D*> > hJeuMBDownEtaDist, 
                                 std::vector< std::vector<TH1D*> > hJeuJet60UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet60DownEtaDist, 
                                 std::vector< std::vector<TH1D*> > hJeuJet80UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet80DownEtaDist, 
                                 std::vector< std::vector<TH1D*> > hJeuJet100UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet100DownEtaDist,
                                 std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hEmbeddingRefEtaDist, std::vector< std::vector<TH1D*> > hEmbeddingGenEtaDist,
                                 std::vector< std::vector<TH1D*> > hJerUpEtaDist, std::vector< std::vector<TH1D*> > hJerDownEtaDist, 
                                 std::vector< std::vector<TH1D*> > hPbGoingEmbeddingRecoEtaDist, std::vector< std::vector<TH1D*> > hPGoingEmbeddingRecoEtaDist, 
                                 std::vector< std::vector<TH1D*> > hPbGoingEmbeddingGenEtaDist, std::vector< std::vector<TH1D*> > hPGoingEmbeddingGenEtaDist,
                                 //std::vector< std::vector<TH1D*> > hEmbeddingPointResDist,
                                 std::vector< std::vector<TH1D*> > hMBGplusEtaDist, std::vector< std::vector<TH1D*> > hMBVtx1EtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet60GplusEtaDist, std::vector< std::vector<TH1D*> > hJet60Vtx1EtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet80GplusEtaDist, std::vector< std::vector<TH1D*> > hJet80Vtx1EtaDist, 
                                 std::vector< std::vector<TH1D*> > hJet100GplusEtaDist, std::vector< std::vector<TH1D*> > hJet100Vtx1EtaDist,
                                 std::vector< std::vector<TH1D*> > hRatios2McDist, std::vector< std::vector<TH1D*> > hUpAndDown2DefDist,
                                 std::vector< std::vector<TH1D*> > hMBJeuRelSystDist, std::vector< std::vector<TH1D*> > hJet60JeuRelSystDist,
                                 std::vector< std::vector<TH1D*> > hJet80JeuRelSystDist, std::vector< std::vector<TH1D*> > hJet100JeuRelSystDist,
                                 std::vector< std::vector<TH1D*> > hEmbeddingJerRelSystDist,
                                 std::vector< std::vector<TH1D*> > hMBPileupRelSystDist, std::vector< std::vector<TH1D*> > hJet60PileupRelSystDist,
                                 std::vector< std::vector<TH1D*> > hJet80PileupRelSystDist, std::vector< std::vector<TH1D*> > hJet100PileupRelSystDist,
                                 std::vector< std::vector<TH1D*> > hMBTotalRelSystDist, std::vector< std::vector<TH1D*> > hJet60TotalRelSystDist,
                                 std::vector< std::vector<TH1D*> > hJet80TotalRelSystDist, std::vector< std::vector<TH1D*> > hJet100TotalRelSystDist,
                                 std::vector< std::vector<TH1D*> > hFinalDist, std::vector< std::vector<TH1D*> > hFinalRelSystUncrtDist, 
                                 std::vector< std::vector<TH1D*> > hFinalAbsSystUncrtDist,
                                 std::vector< std::vector<TH1D*> > hEpps21Dist, std::vector< std::vector<TH1D*> > hNcteq15hqDist,
                                 TString date, Bool_t isCM = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    // Write all histograms from all vectors to the output file
    TFile *oFile = new TFile( Form("%s/oSystematics_%s.root", date.Data(), frame.Data()), "recreate");
    writeVectors2RootFile(hMBEtaDist);
    writeVectors2RootFile(hMBPbGoingEtaDist);
    writeVectors2RootFile(hMBPGoingEtaDist);
    writeVectors2RootFile(hJet60EtaDist);
    writeVectors2RootFile(hJet60PbGoingEtaDist);
    writeVectors2RootFile(hJet60PGoingEtaDist);
    writeVectors2RootFile(hJet80EtaDist);
    writeVectors2RootFile(hJet80PbGoingEtaDist);
    writeVectors2RootFile(hJet80PGoingEtaDist);
    writeVectors2RootFile(hJet100EtaDist);
    writeVectors2RootFile(hJet100PbGoingEtaDist);
    writeVectors2RootFile(hJet100PGoingEtaDist);

    writeVectors2RootFile(hJeuMBUpEtaDist);
    writeVectors2RootFile(hJeuMBDownEtaDist);
    writeVectors2RootFile(hJeuJet60UpEtaDist);
    writeVectors2RootFile(hJeuJet60DownEtaDist);
    writeVectors2RootFile(hJeuJet80UpEtaDist);
    writeVectors2RootFile(hJeuJet80DownEtaDist);
    writeVectors2RootFile(hJeuJet100UpEtaDist);
    writeVectors2RootFile(hJeuJet100DownEtaDist);

    writeVectors2RootFile(hEmbeddingEtaDist);
    writeVectors2RootFile(hEmbeddingRefEtaDist);
    writeVectors2RootFile(hEmbeddingGenEtaDist);
    writeVectors2RootFile(hJerUpEtaDist);
    writeVectors2RootFile(hJerDownEtaDist);
    writeVectors2RootFile(hPbGoingEmbeddingRecoEtaDist);
    writeVectors2RootFile(hPGoingEmbeddingRecoEtaDist);
    //writeVectors2RootFile(hEmbeddingPointResDist);

    writeVectors2RootFile(hMBGplusEtaDist);
    writeVectors2RootFile(hMBVtx1EtaDist);
    writeVectors2RootFile(hJet60GplusEtaDist);
    writeVectors2RootFile(hJet60Vtx1EtaDist);
    writeVectors2RootFile(hJet80GplusEtaDist);
    writeVectors2RootFile(hJet80Vtx1EtaDist);
    writeVectors2RootFile(hJet100GplusEtaDist);
    writeVectors2RootFile(hJet100Vtx1EtaDist);

    writeVectors2RootFile(hRatios2McDist);
    writeVectors2RootFile(hUpAndDown2DefDist);

    writeVectors2RootFile(hMBJeuRelSystDist);
    writeVectors2RootFile(hJet60JeuRelSystDist);
    writeVectors2RootFile(hJet80JeuRelSystDist);
    writeVectors2RootFile(hJet100JeuRelSystDist);

    writeVectors2RootFile(hEmbeddingJerRelSystDist);

    writeVectors2RootFile(hMBPileupRelSystDist);
    writeVectors2RootFile(hJet60PileupRelSystDist);
    writeVectors2RootFile(hJet80PileupRelSystDist);
    writeVectors2RootFile(hJet100PileupRelSystDist);

    writeVectors2RootFile(hMBTotalRelSystDist);
    writeVectors2RootFile(hJet60TotalRelSystDist);
    writeVectors2RootFile(hJet80TotalRelSystDist);
    writeVectors2RootFile(hJet100TotalRelSystDist);

    writeVectors2RootFile(hFinalDist);
    writeVectors2RootFile(hFinalRelSystUncrtDist);
    writeVectors2RootFile(hFinalAbsSystUncrtDist);

    writeVectors2RootFile(hEpps21Dist);
    writeVectors2RootFile(hNcteq15hqDist);

    oFile->Write();
    oFile->Close();
}

//________________
void applySmoothing(TH1D* h, Int_t systType = 0) {
    TString name = h->GetName();
    name.ToLower();
    if ( ( ( systType == 0 || systType == 2 ) && 
           ( !name.Contains("fb") && !name.Contains("forward") && 
             !name.Contains("backward") && !name.Contains("bf") ) ) ||
          ( ( systType == 0 || systType == 1 ) && name.Contains("fb") ) ) {
        Int_t binLow = ( h->GetXaxis()->GetBinCenter(1) < 0 ) ? h->FindBin(-1.2) : 1;
        Int_t binHi = h->FindBin(1.2);
        Double_t maxVal{-1};
        for (Int_t i = binLow; i<=binHi; i++) {
            if (maxVal < h->GetBinContent(i) ) {
                maxVal = h->GetBinContent(i);
            }
        } // for (Int_t i = binLow; i<=binHi; i++)

        for (Int_t i = 1; i<=h->GetNbinsX(); i++) {
            if ( h->GetBinContent(i) < maxVal) {
                h->SetBinContent( i, maxVal );
            }
        } // for (Int_t i = binLow; i<=binHi; i++)
    } // if ( name.Contains("jeuSyst") || name.Contains("pointingSyst") )
}

//________________
TH1D* calculateSystUncrtBinByBin(TH1D *hUp, TH1D *hDown = nullptr, 
                                const Char_t *setName = "hMB", 
                                Int_t systType = 0, Int_t bin = 0,
                                Bool_t useFit = kFALSE) {

    TString hName = setName;
    TString systSuffix;
    Int_t color;
    if (systType == 0) {
        systSuffix = "jeu";
        color = 4; // blue
    } 
    else if (systType == 1) {
        systSuffix = "jer";
        color = 6; // magenta
    }
    else if (systType == 2) {
        systSuffix = "pointing";
        color = 8; // green

    } 
    else if (systType == 3) {
        systSuffix = "pileup";
        color = 9; // red
    }

    hName += "_";
    hName += systSuffix;
    hName += "_";
    hName += "relSystUncrt";
    hName += "_";
    hName += bin;
    // Create histogram to store systematic uncertainties
    TH1D *hSystUncrt = dynamic_cast<TH1D*>( hUp->Clone( hName.Data() ) );
    hSystUncrt->Reset();
    TF1 *fitUp = nullptr;
    TF1 *fitDown = nullptr;
    if ( useFit ) {
        fitUp = hUp->GetFunction( Form( "%s_fitFunc", hUp->GetName() ) );
        if ( hDown ) {
            fitDown = hDown->GetFunction( Form( "%s_fitFunc", hDown->GetName() ) );
        }
    }

    // Calculate systematic uncertainties
    for (Int_t i = 1; i <= hUp->GetNbinsX(); i++) {
        Double_t up = hUp->GetBinContent(i);
        if ( useFit ) {
            up = fitUp->Eval( hUp->GetBinCenter(i) );
        }
        Double_t down = up; 
        if ( hDown ) {
            down = hDown->GetBinContent(i);
            if ( useFit ) {
                down = fitDown->Eval( hDown->GetBinCenter(i) );
            }
        }
        Double_t systUncrt = ( TMath::Abs(up - 1) + TMath::Abs(down - 1) ) / 2.0;
        if ( up < 0.000001 && down < 0.000001 ) {
            systUncrt = 0;
        }
        //std::cout << "bin: " << i << " up: " << up << " down: " << down << " systUncrt: " << systUncrt << std::endl;
        hSystUncrt->SetBinContent(i, systUncrt);
    }
    applySmoothing(hSystUncrt, systType);
    hSystUncrt->SetLineColor(color);
    hSystUncrt->SetMarkerColor(color);
    return hSystUncrt;
}

//________________
std::vector<TH1D*> calculateSystUncrt(std::vector<TH1D*> hUp, std::vector<TH1D*> hDown = {},
                                      const Char_t *setName = "hMB", Int_t systType = 0,
                                      Bool_t useFit = kFALSE) {
 
    // systType: 0 - JEU, 1 - JER, 2 - Pointing, 3 - Pileup

    // Create vector to store systematic uncertainties
    std::vector<TH1D*> hSystUncrt{};

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Calculate systematic uncertainties
    for (Int_t i = 0; i < hUp.size(); i++) {
        if ( systType != 3 ) {
            hSystUncrt.push_back( calculateSystUncrtBinByBin(hUp.at(i), ( !hDown.empty() ) ? hDown.at(i) : nullptr, setName, systType, i, useFit) );
        }
        else {
            if ( i <=9 ) {
                hSystUncrt.push_back( calculateSystUncrtBinByBin(hUp.at(i), ( !hDown.empty() ) ? hDown.at(i) : nullptr, setName, systType, i, useFit) );
            }
            else {
                hSystUncrt.push_back( calculateSystUncrtBinByBin(hUp.at(9), ( !hDown.empty() ) ? hDown.at(9) : nullptr, setName, systType, i, useFit) );
            }
        }
    } // for (Int_t i = 0; i < hUp.size(); i++)

    return hSystUncrt;
}

//________________
std::vector< std::vector<TH1D*> > createRelSystUncrt(std::vector< std::vector<TH1D*> > hUp, 
                                                     std::vector< std::vector<TH1D*> > hDown,
                                                     const Char_t* setName = "hMB", Int_t systType = 0,
                                                     Bool_t useFit = kFALSE) {

    // systType: 0 - JEU, 1 - JER, 2 - Pointing, 3 - Pileup

    TString datasetName = setName;
    // Create vector to store relative systematic uncertainties
    std::vector< std::vector<TH1D*> > hRelSystUncrt{};

    // Vectors to process
    std::vector< TH1D* > hUpEta = hUp.at(0); std::vector< TH1D* > hUpEtaForward = hUp.at(1);
    std::vector< TH1D* > hUpEtaBackward = hUp.at(2); std::vector< TH1D* > hUpEtaFBRatio = hUp.at(3);
    std::vector< TH1D* > hUpEtaBFRatio = hUp.at(4);

    std::vector< TH1D* > hDownEta{};
    std::vector< TH1D* > hDownEtaForward{};
    std::vector< TH1D* > hDownEtaBackward{};
    std::vector< TH1D* > hDownEtaFBRatio{};
    std::vector< TH1D* > hDownEtaBFRatio{};
    if ( !hDown.empty() ) {
        hDownEta = hDown.at(0); hDownEtaForward = hDown.at(1);
        hDownEtaBackward = hDown.at(2); hDownEtaFBRatio = hDown.at(3);
        hDownEtaBFRatio = hDown.at(4);
    }

    // Create relative systematic uncertainties
    datasetName = setName; datasetName += "Eta";         hRelSystUncrt.push_back( calculateSystUncrt(hUpEta, hDownEta, datasetName.Data(), systType, useFit) );
    datasetName = setName; datasetName += "EtaForward";  hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaForward, hDownEtaForward, datasetName.Data(), systType, useFit) );
    datasetName = setName; datasetName += "EtaBackward"; hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaBackward, hDownEtaBackward, datasetName.Data(), systType, useFit) );
    datasetName = setName; datasetName += "EtaFB";       hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaFBRatio, hDownEtaFBRatio, datasetName.Data(), systType, useFit) );
    datasetName = setName; datasetName += "EtaBF";       hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaBFRatio, hDownEtaBFRatio, datasetName.Data(), systType, useFit) );

    return hRelSystUncrt;
}

//________________
std::vector<TH1D*> combineRelSystUncrt(std::vector<TH1D*> hJeuRelSystUncrt, std::vector<TH1D*> hJerRelSystUncrt, 
                                       std::vector<TH1D*> hPointingResRelSystUncert, std::vector<TH1D*> hPileupRelSystUncrt,
                                       const Char_t *setName = "hMB") {
    std::vector<TH1D*> hTotalRelSystUncrt{};
    TString hName = setName;

    // pT bins
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    for (Int_t i = 0; i < ptBins; i++) {
        hName = setName;
        hName += "_";
        hName += "totalSystUncrt";
        hName += "_";
        hName += i;
        TH1D *hTotal = (TH1D*)hJeuRelSystUncrt.at(i)->Clone( hName.Data() );
        hTotal->Reset();
        std::cout << "pT bin: " << i << std::endl;
        for (Int_t j = 1; j <= hJeuRelSystUncrt.at(i)->GetNbinsX(); j++) {
            Double_t jeu = hJeuRelSystUncrt.at(i)->GetBinContent(j);
            Double_t jer = hJerRelSystUncrt.at(i)->GetBinContent(j);
            Double_t pointingRes = hPointingResRelSystUncert.at(i)->GetBinContent(j);
            Double_t pileup = hPileupRelSystUncrt.at(i)->GetBinContent(j);
            Double_t total = TMath::Sqrt( jeu * jeu + jer * jer + pointingRes * pointingRes + pileup * pileup );
            hTotal->SetBinContent(j, total);
            std::cout << "eta bin: " << j << " jeu: " << jeu << " jer: " << jer << " pointingRes: " << pointingRes << " pileup: " << pileup << " total: " << total << std::endl;
            if ( total < jeu ) std::cout << "total < jeu" << std::endl;
            if ( total < jer ) std::cout << "total < jer" << std::endl;
            if ( total < pointingRes ) std::cout << "total < pointingRes" << std::endl;
            if ( total < pileup ) std::cout << "total < pileup" << std::endl;
        }
        hTotal->SetLineColor(1);
        hTotal->SetMarkerColor(1);
        hTotalRelSystUncrt.push_back(hTotal);
    }

    return hTotalRelSystUncrt;
}

//________________
std::vector< std::vector<TH1D*> > createTotalRelSystUncrt(std::vector< std::vector<TH1D*> > hJeuRelSystDist, 
                                                          std::vector< std::vector<TH1D*> > hJerRelSystDist, 
                                                          std::vector< std::vector<TH1D*> > hPointingResRelSystDist,
                                                          std::vector< std::vector<TH1D*> > hPileupRelSystDist, 
                                                          const Char_t *setName = "hMB") {

    // Create vector to store total relative systematic uncertainties
    std::vector< std::vector<TH1D*> > hTotalRelSystUncrt{};

    // Retrieve vectors
    std::vector< TH1D* > hEtaJeuRelSyst = hJeuRelSystDist.at(0); 
    std::vector< TH1D* > hEtaForwardJeuRelSyst = hJeuRelSystDist.at(1);
    std::vector< TH1D* > hEtaBackwardJeuRelSyst = hJeuRelSystDist.at(2);
    std::vector< TH1D* > hEtaFBJeuRelSyst = hJeuRelSystDist.at(3);
    std::vector< TH1D* > hEtaBFJeuRelSyst = hJeuRelSystDist.at(4);

    std::vector< TH1D* > hEtaJerRelSyst = hJerRelSystDist.at(0);
    std::vector< TH1D* > hEtaForwardJerRelSyst = hJerRelSystDist.at(1);
    std::vector< TH1D* > hEtaBackwardJerRelSyst = hJerRelSystDist.at(2);
    std::vector< TH1D* > hEtaFBJerRelSyst = hJerRelSystDist.at(3);
    std::vector< TH1D* > hEtaBFJerRelSyst = hJerRelSystDist.at(4);

    std::vector< TH1D* > hEtaPointingResRelSyst = hPointingResRelSystDist.at(0);
    std::vector< TH1D* > hEtaForwardPointingResRelSyst = hPointingResRelSystDist.at(1);
    std::vector< TH1D* > hEtaBackwardPointingResRelSyst = hPointingResRelSystDist.at(2);
    std::vector< TH1D* > hEtaFBPointingResRelSyst = hPointingResRelSystDist.at(3);
    std::vector< TH1D* > hEtaBFPointingResRelSyst = hPointingResRelSystDist.at(4);

    std::vector< TH1D* > hEtaPileupRelSyst = hPileupRelSystDist.at(0);
    std::vector< TH1D* > hEtaForwardPileupRelSyst = hPileupRelSystDist.at(1);
    std::vector< TH1D* > hEtaBackwardPileupRelSyst = hPileupRelSystDist.at(2);
    std::vector< TH1D* > hEtaFBPileupRelSyst = hPileupRelSystDist.at(3);
    std::vector< TH1D* > hEtaBFPileupRelSyst = hPileupRelSystDist.at(4);

    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaJeuRelSyst, hEtaJerRelSyst, hEtaPointingResRelSyst, hEtaPileupRelSyst, Form("%sEta", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaForwardJeuRelSyst, hEtaForwardJerRelSyst, hEtaForwardPointingResRelSyst, hEtaForwardPileupRelSyst, Form("%sEtaForward", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaBackwardJeuRelSyst, hEtaBackwardJerRelSyst, hEtaBackwardPointingResRelSyst, hEtaBackwardPileupRelSyst, Form("%sEtaBackward", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaFBJeuRelSyst, hEtaFBJerRelSyst, hEtaFBPointingResRelSyst, hEtaFBPileupRelSyst, Form("%sEtaFB", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaBFJeuRelSyst, hEtaBFJerRelSyst, hEtaBFPointingResRelSyst, hEtaBFPileupRelSyst, Form("%sEtaBF", setName) ) );

    return hTotalRelSystUncrt;
}

//________________
std::vector< std::vector<TH1D*> > makeFinalDistributions(std::vector< std::vector<TH1D*> > hMBEtaDist, 
                                                         std::vector< std::vector<TH1D*> > hJet60EtaDist,
                                                         std::vector< std::vector<TH1D*> > hJet80EtaDist, 
                                                         std::vector< std::vector<TH1D*> > hJet100EtaDist,
                                                         std::vector< std::vector<TH1D*> > hGenEtaDist,
                                                         std::vector< std::vector<TH1D*> > hEmbeddingEtaDist,
                                                         std::vector< std::vector<TH1D*> > hRatios2McDist) {

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    //
    // Vectors to fill
    //
    std::vector< std::vector<TH1D*> > hFinalDistributions{};

    // Data to fill
    std::vector< TH1D* > hEtaFinal;                     // 0
    std::vector< TH1D* > hEtaForwardFinal;              // 1
    std::vector< TH1D* > hEtaBackwardFinal;             // 2
    std::vector< TH1D* > hEtaFBRatioFinal;              // 3
    std::vector< TH1D* > hEtaBFRatioFinal;              // 4
    std::vector< TH1D* > hEtaFinalGen;                  // 5
    std::vector< TH1D* > hEtaForwardFinalGen;           // 6
    std::vector< TH1D* > hEtaBackwardFinalGen;          // 7
    std::vector< TH1D* > hEtaFBRatioFinalGen;           // 8
    std::vector< TH1D* > hEtaBFRatioFinalGen;           // 9
    std::vector< TH1D* > hEtaFinalEmbedding;            // 10
    std::vector< TH1D* > hEtaForwardFinalEmbedding;     // 11
    std::vector< TH1D* > hEtaBackwardFinalEmbedding;    // 12
    std::vector< TH1D* > hEtaFBRatioFinalEmbedding;     // 13
    std::vector< TH1D* > hEtaBFRatioFinalEmbedding;     // 14

    // Ratios to Monte Carlo to fill
    std::vector< TH1D* > hEtaData2GenRatioFinal;               // 15
    std::vector< TH1D* > hEtaForwardData2GenRatioFinal;        // 16
    std::vector< TH1D* > hEtaBackwardData2GenRatioFinal;       // 17
    std::vector< TH1D* > hEtaFBRatioData2GenRatioFinal;        // 18
    std::vector< TH1D* > hEtaBFRatioData2GenRatioFinal;        // 19

    std::vector< TH1D* > hEtaEmedding2GenRatioFinal;           // 20
    std::vector< TH1D* > hEtaForwardEmbedding2GenRatioFinal;   // 21
    std::vector< TH1D* > hEtaBackwardEmbedding2GenRatioFinal;  // 22
    std::vector< TH1D* > hEtaFBRatioEmbedding2GenRatioFinal;   // 23
    std::vector< TH1D* > hEtaBFRatioEmbedding2GenRatioFinal;   // 24

    std::vector< TH1D* > hEtaData2EmbeddingRatioFinal;         // 25
    std::vector< TH1D* > hEtaForwardData2EmbeddingRatioFinal;  // 26
    std::vector< TH1D* > hEtaBackwardData2EmbeddingRatioFinal; // 27
    std::vector< TH1D* > hEtaFBRatioData2EmbeddingRatioFinal;  // 28
    std::vector< TH1D* > hEtaBFRatioData2EmbeddingRatioFinal;  // 29
    
    //
    // Read input vectors
    //

    // Data
    std::vector< TH1D* > hMBEta = hMBEtaDist.at(0); std::vector< TH1D* > hMBEtaForward = hMBEtaDist.at(1); 
    std::vector< TH1D* > hMBEtaBackward = hMBEtaDist.at(2); std::vector< TH1D* > hMBEtaFBRatio = hMBEtaDist.at(3);
    std::vector< TH1D* > hMBEtaBFRatio = hMBEtaDist.at(4);

    std::vector< TH1D* > hJet60Eta = hJet60EtaDist.at(0); std::vector< TH1D* > hJet60EtaForward = hJet60EtaDist.at(1);
    std::vector< TH1D* > hJet60EtaBackward = hJet60EtaDist.at(2); std::vector< TH1D* > hJet60EtaFBRatio = hJet60EtaDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatio = hJet60EtaDist.at(4);

    std::vector< TH1D* > hJet80Eta = hJet80EtaDist.at(0); std::vector< TH1D* > hJet80EtaForward = hJet80EtaDist.at(1);
    std::vector< TH1D* > hJet80EtaBackward = hJet80EtaDist.at(2); std::vector< TH1D* > hJet80EtaFBRatio = hJet80EtaDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatio = hJet80EtaDist.at(4);

    std::vector< TH1D* > hJet100Eta = hJet100EtaDist.at(0); std::vector< TH1D* > hJet100EtaForward = hJet100EtaDist.at(1);
    std::vector< TH1D* > hJet100EtaBackward = hJet100EtaDist.at(2); std::vector< TH1D* > hJet100EtaFBRatio = hJet100EtaDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatio = hJet100EtaDist.at(4);

    // Gen
    std::vector< TH1D* > hGenEta = hGenEtaDist.at(0); std::vector< TH1D* > hGenEtaForward = hGenEtaDist.at(1);
    std::vector< TH1D* > hGenEtaBackward = hGenEtaDist.at(2); std::vector< TH1D* > hGenEtaFBRatio = hGenEtaDist.at(3);
    std::vector< TH1D* > hGenEtaBFRatio = hGenEtaDist.at(4);

    // Embedding
    std::vector< TH1D* > hEmbeddingEta = hEmbeddingEtaDist.at(0); std::vector< TH1D* > hEmbeddingEtaForward = hEmbeddingEtaDist.at(1);
    std::vector< TH1D* > hEmbeddingEtaBackward = hEmbeddingEtaDist.at(2); std::vector< TH1D* > hEmbeddingEtaFBRatio = hEmbeddingEtaDist.at(3);
    std::vector< TH1D* > hEmbeddingEtaBFRatio = hEmbeddingEtaDist.at(4);

    // Ratios to Monte Carlo
    std::vector<TH1D*> hMBEta2GenEtaRatio = hRatios2McDist.at(0);
    std::vector<TH1D*> hMBEtaForward2GenEtaRatio = hRatios2McDist.at(1);
    std::vector<TH1D*> hMBEtaBackward2GenEtaRatio = hRatios2McDist.at(2);
    std::vector<TH1D*> hMBEtaFBRatio2GenEtaRatio = hRatios2McDist.at(3);
    std::vector<TH1D*> hMBEtaBFRatio2GenEtaRatio = hRatios2McDist.at(4);

    std::vector<TH1D*> hJet60Eta2GenEtaRatio = hRatios2McDist.at(5);
    std::vector<TH1D*> hJet60EtaForward2GenEtaRatio = hRatios2McDist.at(6);
    std::vector<TH1D*> hJet60EtaBackward2GenEtaRatio = hRatios2McDist.at(7);
    std::vector<TH1D*> hJet60EtaFBRatio2GenEtaRatio = hRatios2McDist.at(8);
    std::vector<TH1D*> hJet60EtaBFRatio2GenEtaRatio = hRatios2McDist.at(9);

    std::vector<TH1D*> hJet80Eta2GenEtaRatio = hRatios2McDist.at(10);
    std::vector<TH1D*> hJet80EtaForward2GenEtaRatio = hRatios2McDist.at(11);
    std::vector<TH1D*> hJet80EtaBackward2GenEtaRatio = hRatios2McDist.at(12);
    std::vector<TH1D*> hJet80EtaFBRatio2GenEtaRatio = hRatios2McDist.at(13);
    std::vector<TH1D*> hJet80EtaBFRatio2GenEtaRatio = hRatios2McDist.at(14);

    std::vector<TH1D*> hJet100Eta2GenEtaRatio = hRatios2McDist.at(15);
    std::vector<TH1D*> hJet100EtaForward2GenEtaRatio = hRatios2McDist.at(16);
    std::vector<TH1D*> hJet100EtaBackward2GenEtaRatio = hRatios2McDist.at(17);
    std::vector<TH1D*> hJet100EtaFBRatio2GenEtaRatio = hRatios2McDist.at(18);
    std::vector<TH1D*> hJet100EtaBFRatio2GenEtaRatio = hRatios2McDist.at(19);

    std::vector<TH1D*> hEmbedding2GenEtaRatio = hRatios2McDist.at(20);
    std::vector<TH1D*> hEmbedding2GenEtaForwardRatio = hRatios2McDist.at(21);
    std::vector<TH1D*> hEmbedding2GenEtaBackwardRatio = hRatios2McDist.at(22);
    std::vector<TH1D*> hEmbedding2GenEtaFBRatio = hRatios2McDist.at(23);
    std::vector<TH1D*> hEmbedding2GenEtaBFRatio = hRatios2McDist.at(24);

    std::vector<TH1D*> hMBEta2EmbeddingEtaRatio = hRatios2McDist.at(25);
    std::vector<TH1D*> hMBEtaForward2EmbeddingEtaRatio = hRatios2McDist.at(26);
    std::vector<TH1D*> hMBEtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(27);
    std::vector<TH1D*> hMBEtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(28);
    std::vector<TH1D*> hMBEtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(29);

    std::vector<TH1D*> hJet60Eta2EmbeddingEtaRatio = hRatios2McDist.at(30);
    std::vector<TH1D*> hJet60EtaForward2EmbeddingEtaRatio = hRatios2McDist.at(31);
    std::vector<TH1D*> hJet60EtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(32);
    std::vector<TH1D*> hJet60EtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(33);
    std::vector<TH1D*> hJet60EtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(34);

    std::vector<TH1D*> hJet80Eta2EmbeddingEtaRatio = hRatios2McDist.at(35);
    std::vector<TH1D*> hJet80EtaForward2EmbeddingEtaRatio = hRatios2McDist.at(36);
    std::vector<TH1D*> hJet80EtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(37);
    std::vector<TH1D*> hJet80EtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(38);
    std::vector<TH1D*> hJet80EtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(39);

    std::vector<TH1D*> hJet100Eta2EmbeddingEtaRatio = hRatios2McDist.at(40);
    std::vector<TH1D*> hJet100EtaForward2EmbeddingEtaRatio = hRatios2McDist.at(41);
    std::vector<TH1D*> hJet100EtaBackward2EmbeddingEtaRatio = hRatios2McDist.at(42);
    std::vector<TH1D*> hJet100EtaFBRatio2EmbeddingEtaRatio = hRatios2McDist.at(43);
    std::vector<TH1D*> hJet100EtaBFRatio2EmbeddingEtaRatio = hRatios2McDist.at(44);


    // Create final distributions
    for (Int_t i = 0; i < ptBins; i++) {
        if ( ptDijetLow.at(i) < 80 ) {
            hEtaFinal.push_back( dynamic_cast<TH1D*>( hMBEta.at(i)->Clone( Form("hData_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardFinal.push_back( dynamic_cast<TH1D*>( hMBEtaForward.at(i)->Clone( Form("hData_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardFinal.push_back( dynamic_cast<TH1D*>( hMBEtaBackward.at(i)->Clone( Form("hData_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio.at(i)->Clone( Form("hData_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio.at(i)->Clone( Form("hData_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );

            hEtaData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEta2GenEtaRatio.at(i)->Clone( Form("hData2Gen_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaForward2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaBackward2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );

            hEtaEmedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaRatio.at(i)->Clone( Form("hEmbedding2Gen_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaForwardRatio.at(i)->Clone( Form("hEmbedding2Gen_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaBackwardRatio.at(i)->Clone( Form("hEmbedding2Gen_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaFBRatio.at(i)->Clone( Form("hEmbedding2Gen_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaBFRatio.at(i)->Clone( Form("hEmbedding2Gen_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );

            hEtaData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEta2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaForward2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaBackward2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        }
        else if ( ptDijetLow.at(i) < 100 ) {
            hEtaFinal.push_back( dynamic_cast<TH1D*>( hJet60Eta.at(i)->Clone( Form("hData_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaForward.at(i)->Clone( Form("hData_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward.at(i)->Clone( Form("hData_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio.at(i)->Clone( Form("hData_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio.at(i)->Clone( Form("hData_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );

            hEtaData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60Eta2GenEtaRatio.at(i)->Clone( Form("hData2Gen_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaForward2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioData2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio2GenEtaRatio.at(i)->Clone( Form("hData2Gen_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        
            hEtaEmedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaRatio.at(i)->Clone( Form("hEmbedding2Gen_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaForwardRatio.at(i)->Clone( Form("hEmbedding2Gen_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaBackwardRatio.at(i)->Clone( Form("hEmbedding2Gen_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaFBRatio.at(i)->Clone( Form("hEmbedding2Gen_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioEmbedding2GenRatioFinal.push_back( dynamic_cast<TH1D*>( hEmbedding2GenEtaBFRatio.at(i)->Clone( Form("hEmbedding2Gen_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );

            hEtaData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60Eta2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaForward2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioData2EmbeddingRatioFinal.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio2EmbeddingEtaRatio.at(i)->Clone( Form("hData2Embedding_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        }
        else if ( ptDijetLow.at(i) < 120) {
            hEtaFinal.push_back( dynamic_cast<TH1D*>( hJet80Eta.at(i)->Clone( Form("hData_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardFinal.push_back( dynamic_cast<TH1D*>( hJet80EtaForward.at(i)->Clone( Form("hData_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardFinal.push_back( dynamic_cast<TH1D*>( hJet80EtaBackward.at(i)->Clone( Form("hData_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioFinal.push_back( dynamic_cast<TH1D*>( hJet80EtaFBRatio.at(i)->Clone( Form("hData_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioFinal.push_back( dynamic_cast<TH1D*>( hJet80EtaBFRatio.at(i)->Clone( Form("hData_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        }
        else {
            hEtaFinal.push_back( dynamic_cast<TH1D*>( hJet100Eta.at(i)->Clone( Form("hData_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaForwardFinal.push_back( dynamic_cast<TH1D*>( hJet100EtaForward.at(i)->Clone( Form("hData_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBackwardFinal.push_back( dynamic_cast<TH1D*>( hJet100EtaBackward.at(i)->Clone( Form("hData_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaFBRatioFinal.push_back( dynamic_cast<TH1D*>( hJet100EtaFBRatio.at(i)->Clone( Form("hData_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
            hEtaBFRatioFinal.push_back( dynamic_cast<TH1D*>( hJet100EtaBFRatio.at(i)->Clone( Form("hData_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        }

        // Gen
        hEtaFinalGen.push_back( dynamic_cast<TH1D*>( hGenEta.at(i)->Clone( Form("hPythia_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaForwardFinalGen.push_back( dynamic_cast<TH1D*>( hGenEtaForward.at(i)->Clone( Form("hPythia_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaBackwardFinalGen.push_back( dynamic_cast<TH1D*>( hGenEtaBackward.at(i)->Clone( Form("hPythia_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaFBRatioFinalGen.push_back( dynamic_cast<TH1D*>( hGenEtaFBRatio.at(i)->Clone( Form("hPythia_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaBFRatioFinalGen.push_back( dynamic_cast<TH1D*>( hGenEtaBFRatio.at(i)->Clone( Form("hPythia_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );

        // Embedding
        hEtaFinalEmbedding.push_back( dynamic_cast<TH1D*>( hEmbeddingEta.at(i)->Clone( Form("hEmbedding_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaForwardFinalEmbedding.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaForward.at(i)->Clone( Form("hEmbedding_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaBackwardFinalEmbedding.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaBackward.at(i)->Clone( Form("hEmbedding_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaFBRatioFinalEmbedding.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaFBRatio.at(i)->Clone( Form("hEmbedding_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
        hEtaBFRatioFinalEmbedding.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaBFRatio.at(i)->Clone( Form("hEmbedding_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) )  ) );
    } // for (Int_t i = 0; i < ptBins; i++)

    hFinalDistributions.push_back(hEtaFinal);                       // 0
    hFinalDistributions.push_back(hEtaForwardFinal);                // 1
    hFinalDistributions.push_back(hEtaBackwardFinal);               // 2
    hFinalDistributions.push_back(hEtaFBRatioFinal);                // 3
    hFinalDistributions.push_back(hEtaBFRatioFinal);                // 4

    hFinalDistributions.push_back(hEtaFinalGen);                    // 5
    hFinalDistributions.push_back(hEtaForwardFinalGen);             // 6
    hFinalDistributions.push_back(hEtaBackwardFinalGen);            // 7
    hFinalDistributions.push_back(hEtaFBRatioFinalGen);             // 8
    hFinalDistributions.push_back(hEtaBFRatioFinalGen);             // 9

    hFinalDistributions.push_back(hEtaFinalEmbedding);              // 10
    hFinalDistributions.push_back(hEtaForwardFinalEmbedding);       // 11
    hFinalDistributions.push_back(hEtaBackwardFinalEmbedding);      // 12
    hFinalDistributions.push_back(hEtaFBRatioFinalEmbedding);       // 13
    hFinalDistributions.push_back(hEtaBFRatioFinalEmbedding);       // 14

    hFinalDistributions.push_back(hEtaData2GenRatioFinal);          // 15
    hFinalDistributions.push_back(hEtaForwardData2GenRatioFinal);   // 16
    hFinalDistributions.push_back(hEtaBackwardData2GenRatioFinal);  // 17
    hFinalDistributions.push_back(hEtaFBRatioData2GenRatioFinal);   // 18
    hFinalDistributions.push_back(hEtaBFRatioData2GenRatioFinal);   // 19

    hFinalDistributions.push_back(hEtaEmedding2GenRatioFinal);      // 20
    hFinalDistributions.push_back(hEtaForwardEmbedding2GenRatioFinal);   // 21
    hFinalDistributions.push_back(hEtaBackwardEmbedding2GenRatioFinal);  // 22
    hFinalDistributions.push_back(hEtaFBRatioEmbedding2GenRatioFinal);   // 23
    hFinalDistributions.push_back(hEtaBFRatioEmbedding2GenRatioFinal);   // 24

    hFinalDistributions.push_back(hEtaData2EmbeddingRatioFinal);    // 25
    hFinalDistributions.push_back(hEtaForwardData2EmbeddingRatioFinal); // 26
    hFinalDistributions.push_back(hEtaBackwardData2EmbeddingRatioFinal); // 27
    hFinalDistributions.push_back(hEtaFBRatioData2EmbeddingRatioFinal); // 28
    hFinalDistributions.push_back(hEtaBFRatioData2EmbeddingRatioFinal); // 29

    return hFinalDistributions;
}

//________________
std::vector< TH1D* > makeRelSystUncrtAtUnity(std::vector<TH1D*> hUncrt) {
    std::vector< TH1D* > hRelSystUncrt{};
    for (Int_t i{0}; i<hUncrt.size(); i++) {
        hRelSystUncrt.push_back( dynamic_cast<TH1D*>( hUncrt.at(i)->Clone( Form("%s_rescaledAtUnity", hUncrt.at(i)->GetName()) ) ) );
        // hRelSystUncrt.at(i)->Scale(100.);
        for (Int_t j{1}; j<=hRelSystUncrt.at(i)->GetNbinsX(); j++) {
            Double_t val = hRelSystUncrt.at(i)->GetBinContent(j);
            hRelSystUncrt.at(i)->SetBinContent(j, 1.0);
            hRelSystUncrt.at(i)->SetBinError(j, val);
        }
    }
    return hRelSystUncrt;
}

//________________
std::vector< std::vector<TH1D*> > makeFinalRelSystUncrt(std::vector< std::vector<TH1D*> > hMBTotalRelSystDist, 
                                                        std::vector< std::vector<TH1D*> > hJet60TotalRelSystDist, 
                                                        std::vector< std::vector<TH1D*> > hJet80TotalRelSystDist, 
                                                        std::vector< std::vector<TH1D*> > hJet100TotalRelSystDist) {
    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    //
    // Vectors to fill
    //
    std::vector< std::vector<TH1D*> > hFinalDistributions{};

    // Vectors for uncertainties
    std::vector<TH1D*> hEtaFinalRelSystUncrt{};
    std::vector<TH1D*> hEtaForwardFinalRelSystUncrt{};
    std::vector<TH1D*> hEtaBackwardFinalRelSystUncrt{};
    std::vector<TH1D*> hEtaFBRatioFinalRelSystUncrt{};
    std::vector<TH1D*> hEtaBFRatioFinalRelSystUncrt{};

    // Vectors with incoming uncertainties
    std::vector< TH1D* > hMBEtaRelSystUncrt = hMBTotalRelSystDist.at(0);
    std::vector< TH1D* > hMBEtaForwardRelSystUncrt = hMBTotalRelSystDist.at(1);
    std::vector< TH1D* > hMBEtaBackwardRelSystUncrt = hMBTotalRelSystDist.at(2);
    std::vector< TH1D* > hMBEtaFBRatioRelSystUncrt = hMBTotalRelSystDist.at(3);
    std::vector< TH1D* > hMBEtaBFRatioRelSystUncrt = hMBTotalRelSystDist.at(4);

    std::vector< TH1D* > hJet60EtaRelSystUncrt = hJet60TotalRelSystDist.at(0);
    std::vector< TH1D* > hJet60EtaForwardRelSystUncrt = hJet60TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet60EtaBackwardRelSystUncrt = hJet60TotalRelSystDist.at(2);
    std::vector< TH1D* > hJet60EtaFBRatioRelSystUncrt = hJet60TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatioRelSystUncrt = hJet60TotalRelSystDist.at(4);

    std::vector< TH1D* > hJet80EtaRelSystUncrt = hJet80TotalRelSystDist.at(0);
    std::vector< TH1D* > hJet80EtaForwardRelSystUncrt = hJet80TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet80EtaBackwardRelSystUncrt = hJet80TotalRelSystDist.at(2);
    std::vector< TH1D* > hJet80EtaFBRatioRelSystUncrt = hJet80TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatioRelSystUncrt = hJet80TotalRelSystDist.at(4);

    std::vector< TH1D* > hJet100EtaRelSystUncrt = hJet100TotalRelSystDist.at(0);
    std::vector< TH1D* > hJet100EtaForwardRelSystUncrt = hJet100TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet100EtaBackwardRelSystUncrt = hJet100TotalRelSystDist.at(2);
    std::vector< TH1D* > hJet100EtaFBRatioRelSystUncrt = hJet100TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatioRelSystUncrt = hJet100TotalRelSystDist.at(4);

    // Loop over pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        if ( ptDijetLow.at(i) < 80 ) {
            hEtaFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hMBEtaRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaForwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hMBEtaForwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBackwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hMBEtaBackwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaFBRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBFRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
        }
        else if ( ptDijetLow.at(i) < 100 ) {
            hEtaFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet60EtaRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaForwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet60EtaForwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBackwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet60EtaBackwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaFBRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBFRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
        }
        else if ( ptDijetLow.at(i) < 120 ) {
            hEtaFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet80EtaRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaForwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet80EtaForwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBackwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet80EtaBackwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaFBRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet80EtaFBRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBFRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet80EtaBFRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
        }
        else {
            hEtaFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet100EtaRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaForwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet100EtaForwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBackwardFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet100EtaBackwardRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaFBRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet100EtaFBRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBFRatioFinalRelSystUncrt.push_back( dynamic_cast<TH1D*>( hJet100EtaBFRatioRelSystUncrt.at(i)->Clone( Form("hDataRelSystUncrt_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
        }
    } // for (Int_t i{0}; i<ptBins; i++)

    hFinalDistributions.push_back(hEtaFinalRelSystUncrt);                       // 0
    hFinalDistributions.push_back(hEtaForwardFinalRelSystUncrt);                // 1
    hFinalDistributions.push_back(hEtaBackwardFinalRelSystUncrt);               // 2
    hFinalDistributions.push_back(hEtaFBRatioFinalRelSystUncrt);                // 3
    hFinalDistributions.push_back(hEtaBFRatioFinalRelSystUncrt);                // 4

    std::vector< TH1D* > hEtaFinalRelSystUncrtAtUnity = makeRelSystUncrtAtUnity(hEtaFinalRelSystUncrt);
    std::vector< TH1D* > hEtaForwardFinalRelSystUncrtAtUnity = makeRelSystUncrtAtUnity(hEtaForwardFinalRelSystUncrt);
    std::vector< TH1D* > hEtaBackwardFinalRelSystUncrtAtUnity = makeRelSystUncrtAtUnity(hEtaBackwardFinalRelSystUncrt);
    std::vector< TH1D* > hEtaFBRatioFinalRelSystUncrtAtUnity = makeRelSystUncrtAtUnity(hEtaFBRatioFinalRelSystUncrt);
    std::vector< TH1D* > hEtaBFRatioFinalRelSystUncrtAtUnity = makeRelSystUncrtAtUnity(hEtaBFRatioFinalRelSystUncrt);

    hFinalDistributions.push_back(hEtaFinalRelSystUncrtAtUnity);                // 5
    hFinalDistributions.push_back(hEtaForwardFinalRelSystUncrtAtUnity);         // 6
    hFinalDistributions.push_back(hEtaBackwardFinalRelSystUncrtAtUnity);        // 7
    hFinalDistributions.push_back(hEtaFBRatioFinalRelSystUncrtAtUnity);         // 8
    hFinalDistributions.push_back(hEtaBFRatioFinalRelSystUncrtAtUnity);         // 9

    return hFinalDistributions;
}

//________________
std::vector< std::vector< TH1D* > > retrieveNLOJetDistributions(TFile *inFile, Int_t npdfSet = 0) {

    // npdfSet: 0 - epps21; 1 - nCTEQ15HQ

    // pT bin edges
    std::vector<int> ptLow = {50, 60, 100, 140, 300}; 
    std::vector<int> ptHi =  {60, 70, 110, 150, 500};

    const char* npdfName = "epps21";
    Int_t nPDFSets = 107;
    if ( npdfSet == 0 ) {
        npdfName = "epps21";
        nPDFSets = 107;
    }
    else if ( npdfSet == 1 ) {
        npdfName = "ncteq15hq";
        nPDFSets = 39;
    }

    // std::vector< std::vector< TH1D* > > hNpdfDistirbutions;
    std::vector< std::vector< TH1D* > > hNpdfEtaLab{};
    std::vector< TH1D* > hNpdfEtaLabPt{};
    // std::vector< TH1D* > hNpdfEtaCMForward{};
    // std::vector< TH1D* > hNpdfEtaCMBackward{};
    // std::vector< TH1D* > hNpdfEtaCMFBRatio{};
    // std::vector< TH1D* > hNpdfEtaCMBFRatio{};

    // Loop over pT bins
    for (Int_t i{0}; i<ptLow.size(); i++) {

        hNpdfEtaLabPt.clear();
        for (Int_t j{0}; j<nPDFSets; j++ ) {

            // Full eta in lab frame
            hNpdfEtaLabPt.push_back( dynamic_cast<TH1D*>( inFile->Get( Form("pPbLab_pt_%d_%d_%d", ptLow.at(i), ptHi.at(i), j) ) ) );
            hNpdfEtaLabPt.at(j)->SetName( Form("%s_pPb8160_etaLab_pt_%d_%d_%d", npdfName, ptLow.at(i), ptHi.at(i), j) );
            //rescaleEta(hNpdfEtaLabPt.at(j));
            hNpdfEtaLabPt.at(j)->Scale( 1./hNpdfEtaLabPt.at(j)->Integral() );

            // // Forward eta in CM frame
            // hNpdfEtaCMForward.push_back( dynamic_cast<TH1D*>( inFile->Get( Form("pPbCMForward_pt_%d_%d_%d", ptLow.at(i), ptHi.at(i), j) ) ) );
            // hNpdfEtaCMForward.at(j)->SetName( Form("%s_pPb8160_etaCM_forward_pt_%d_%d_%d", npdfName, ptLow.at(i), ptHi.at(i), j) );

            // // Backward eta in CM frame
            // hNpdfEtaCMBackward.push_back( dynamic_cast<TH1D*>( inFile->Get( Form("pPbCMBackward_pt_%d_%d_%d", ptLow.at(i), ptHi.at(i), j) ) ) );
            // hNpdfEtaCMBackward.at(j)->SetName( Form("%s_pPb8160_etaCM_backward_pt_%d_%d_%d", npdfName, ptLow.at(i), ptHi.at(i), j  );

            // // Forward/Backward ratio in CM frame
            // hNpdfEtaCMFBRatio.push_back( dynamic_cast<TH1D*>( hNpdfEtaCMForward.at(i)->Clone( Form("%s_pPb8160_etaCM_fb_pt_%d_%d_%d", npdfName, ptLow.at(i), ptHi.at(i), j) ) ) );
            // hNpdfEtaCMFBRatio.at(j)->Divide( hNpdfEtaCMBackward.at(i) );

            // // Backward/Forward ratio in CM frame
            // hNpdfEtaCMBFRatio.push_back( dynamic_cast<TH1D*>( hNpdfEtaCMBackward.at(i)->Clone( Form("%s_pPb8160_etaCM_bf_pt_%d_%d_%d", npdfName, ptLow.at(i), ptHi.at(i), j) ) ) );
            // hNpdfEtaCMBFRatio.at(j)->Divide( hNpdfEtaCMForward.at(i) );
        } // for (Int_t j{0}; j<nPDFSets; j++ )
        hNpdfEtaLab.push_back(hNpdfEtaLabPt);
    } // for (Int_t i{0}; i<ptLow.size(); i++)

    // hNpdfDistirbutions.push_back(hNpdfEtaLab);
    // hNpdfDistirbutions.push_back(hNpdfEtaCMForward);
    // hNpdfDistirbutions.push_back(hNpdfEtaCMBackward);
    // hNpdfDistirbutions.push_back(hNpdfEtaCMFBRatio);
    // hNpdfDistirbutions.push_back(hNpdfEtaCMBFRatio);

    return hNpdfEtaLab;
}

//________________
TH1D* calculateSystUncrtForNPDF(std::vector<TH1D*> hRatios, Int_t npdfSet = 0, Int_t ptLow = 50, Int_t ptHi = 60) {

    const char* npdfName = "data2epps21";
    if ( npdfSet == 1 ) npdfName = "data2ncteq15hq";

    TH1D* hSystUncrt = dynamic_cast<TH1D*>( hRatios.at(0)->Clone( Form("%s_systUncrt_pt_%d_%d", npdfName, ptLow, ptHi)) );
    std::vector<TH1D*> hSystSidePositive{};
    std::vector<TH1D*> hSystSideNegative{};
    TH1D* hErrors = dynamic_cast<TH1D*>( hRatios.at(0)->Clone( Form("%s_errors", hRatios.at(0)->GetName())) );
    hErrors->Reset();


    // Try to estimate simmetric uncertainties accroding
    // the Eq. 10 from https://arxiv.org/pdf/1903.09832
    Int_t iter = 0;
    for (Int_t i{1}; i<hRatios.size(); ) {
        // Negative side
        hSystSideNegative.push_back( dynamic_cast<TH1D*>( hRatios.at(i)->Clone( Form("%s_systNegative_%d", hRatios.at(i)->GetName(), iter) ) ) );
        // Positive side
        hSystSidePositive.push_back( dynamic_cast<TH1D*>( hRatios.at(i+1)->Clone( Form("%s_systPositive_%d", hRatios.at(i+1)->GetName(), iter) ) ) );

        // Square the numbers and sum with the previous set
        for (Int_t j{1}; j<=hErrors->GetNbinsX(); j++) {
            Double_t oldVal = hErrors->GetBinContent(j);
            Double_t newVal = hSystSidePositive.at(iter)->GetBinContent(j) - hSystSideNegative.at(iter)->GetBinContent(j);
            newVal *= newVal;
            newVal += oldVal;
            hErrors->SetBinContent(j, newVal);
        }
        i+=2;
        iter++;
    } // for (Int_t i{1}; i<hRatios.size(); )

    // Take a half of the root of each bin and assign it as a final uncertainty
    // to the main nPDF set
    for (Int_t i{1}; i<=hSystUncrt->GetNbinsX(); i++) {
        Double_t val = hSystUncrt->GetBinContent(i);
        Double_t valErr = hErrors->GetBinContent(i);
        valErr = 0.5 * TMath::Sqrt(valErr);
        hSystUncrt->SetBinContent(i, val);
        hSystUncrt->SetBinError(i, valErr);

        //std::cout << "Bin: " << i << " eta: " << hSystUncrt->GetBinCenter(i) << " value: " << val << " error: " << valErr << std::endl;
    }

    return hSystUncrt;
}

//________________
std::vector< std::vector< TH1D* > > makeRatiosOfData2nPDF(std::vector< std::vector< TH1D* > > hDataDist, 
                                                          std::vector< std::vector< TH1D* > > hNpdfDist,
                                                          Int_t npdfSet = 0) {

    // npdfSet: 0 - epps21; 1 - ncteq15hq
    const char* npdfName = "data2epps21";
    if ( npdfSet == 1 ) npdfName = "data2ncteq15hq";

    // Dijet pT standard binning
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Vector that corresponds to all pT bins
    std::vector< TH1D* > hDataEta = hDataDist.at(0);
    std::vector< std::vector< TH1D* > > hReturnVector{};
    std::vector< TH1D* > hData2NpdfStat{};
    std::vector< TH1D* > hData2NpdfSyst{};

    // For all existing nPDF calculated pT bins
    for (Int_t i{0}; i<npdfPtLow.size(); i++) {

        // Skip bins that are not matching pT binning
        Int_t ptBin = -1;
        for (Int_t j{0}; j<ptDijetLow.size(); j++) {
            if ( npdfPtLow.at(i) == ptDijetLow.at(j) && npdfPtHi.at(i) == ptDijetHi.at(j) ) {
                ptBin = j;
                std::cout << "Found matching pT bin: " << ptBin << " ptLow: " << npdfPtLow.at(i) << " ptHi: " << npdfPtHi.at(i) << std::endl;
                break;
            }
        }

        // Retrieve data distribution for the given MC pT bin
        TH1D* hData = hDataEta.at(ptBin);
        // Retrieve nPDF distributions for the given MC pT bin
        std::vector<TH1D*> nloSets = hNpdfDist.at(i);
        std::vector<TH1D*> hData2NpdfRatio{};

        // Loop over nPDF sets
        for (Int_t j{0}; j<nloSets.size(); j++) {
            if ( j == 0 ) {
                hData2NpdfRatio.push_back( dynamic_cast<TH1D*>( hData->Clone( Form("%s_statUcrt_eta_pt_%d_%d", npdfName, npdfPtLow.at(i), npdfPtHi.at(i)) ) ) );
            }
            else {
                hData2NpdfRatio.push_back( dynamic_cast<TH1D*>( hData->Clone( Form("%s_statUcrt_eta_pt_%d_%d_%d", npdfName, npdfPtLow.at(i), npdfPtHi.at(i), j) ) ) );
            }
            hData2NpdfRatio.at(j)->Divide( nloSets.at(j) );
        } // for (Int_t j{0}; j<nloSets.size(); j++)

        // Calculate systematic uncertainties for the given nPDF set
        hData2NpdfStat.push_back( hData2NpdfRatio.at(0) );
        hData2NpdfSyst.push_back( calculateSystUncrtForNPDF(hData2NpdfRatio, npdfSet, npdfPtLow.at(i), npdfPtHi.at(i)) );
    } // for (Int_t i{0}; i<ptLow.size(); i++)

    // Ratio with statistical uncertainties
    hReturnVector.push_back(hData2NpdfStat);
    // Ratio with systematic uncertainties
    hReturnVector.push_back(hData2NpdfSyst);

    return hReturnVector;
}

//________________
std::vector< std::vector<TH1D*> > makeFinalAbsSystUncrt(std::vector< std::vector<TH1D*> > hFinalDist, 
                                                        std::vector< std::vector<TH1D*> > hFinalRelSystUncrtDist) {
    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    //
    // Vectors to fill
    //
    std::vector< std::vector<TH1D*> > hFinalDistributions{};

    // Vectors for uncertainties
    std::vector<TH1D*> hEtaFinalAbsSystUncrt{};
    std::vector<TH1D*> hEtaForwardFinalAbsSystUncrt{};
    std::vector<TH1D*> hEtaBackwardFinalAbsSystUncrt{};
    std::vector<TH1D*> hEtaFBRatioFinalAbsSystUncrt{};
    std::vector<TH1D*> hEtaBFRatioFinalAbsSystUncrt{};

    // Vectors with incoming uncertainties
    std::vector< TH1D* > hEtaFinal = hFinalDist.at(0);
    std::vector< TH1D* > hEtaForwardFinal = hFinalDist.at(1);
    std::vector< TH1D* > hEtaBackwardFinal = hFinalDist.at(2);
    std::vector< TH1D* > hEtaFBRatioFinal = hFinalDist.at(3);
    std::vector< TH1D* > hEtaBFRatioFinal = hFinalDist.at(4);

    std::vector< TH1D* > hEtaFinalRelSystUncrt = hFinalRelSystUncrtDist.at(0);
    std::vector< TH1D* > hEtaForwardFinalRelSystUncrt = hFinalRelSystUncrtDist.at(1);
    std::vector< TH1D* > hEtaBackwardFinalRelSystUncrt = hFinalRelSystUncrtDist.at(2);
    std::vector< TH1D* > hEtaFBRatioFinalRelSystUncrt = hFinalRelSystUncrtDist.at(3);
    std::vector< TH1D* > hEtaBFRatioFinalRelSystUncrt = hFinalRelSystUncrtDist.at(4);

    // Loop over pT bins
    for (Int_t i{0}; i<ptBins; i++) {
            
            hEtaFinalAbsSystUncrt.push_back( dynamic_cast<TH1D*>( hEtaFinal.at(i)->Clone( Form("hDataAbsSystUncrt_eta_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaForwardFinalAbsSystUncrt.push_back( dynamic_cast<TH1D*>( hEtaForwardFinal.at(i)->Clone( Form("hDataAbsSystUncrt_etaForward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBackwardFinalAbsSystUncrt.push_back( dynamic_cast<TH1D*>( hEtaBackwardFinal.at(i)->Clone( Form("hDataAbsSystUncrt_etaBackward_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaFBRatioFinalAbsSystUncrt.push_back( dynamic_cast<TH1D*>( hEtaFBRatioFinal.at(i)->Clone( Form("hDataAbsSystUncrt_etaFB_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
            hEtaBFRatioFinalAbsSystUncrt.push_back( dynamic_cast<TH1D*>( hEtaBFRatioFinal.at(i)->Clone( Form("hDataAbsSystUncrt_etaBF_pt_%d_%d", ptDijetLow.at(i), ptDijetHi.at(i)) ) ) );
    
            // Absolute uncertainties for full eta
            for (Int_t j{1}; j<hEtaFinalAbsSystUncrt.at(i)->GetNbinsX()+1; j++) {
                Double_t xVal{-999.};
                Double_t yVal{-999.};
                Double_t relUncrt{0.};
                Double_t absUncrt{0.};

                xVal = hEtaFinalAbsSystUncrt.at(i)->GetBinCenter(j);
                yVal = hEtaFinal.at(i)->GetBinContent(j);
                relUncrt = hEtaFinalRelSystUncrt.at(i)->GetBinContent(j); // Check if it is in percent or absolute value
                absUncrt = hEtaFinal.at(i)->GetBinContent(j) * relUncrt;
                hEtaFinalAbsSystUncrt.at(i)->SetBinContent(j, yVal);
                hEtaFinalAbsSystUncrt.at(i)->SetBinError(j, absUncrt);
            } //  for (Int_t j{1}; j< hEtaFinalAbsSystUncrt.at(i)->GetNbinsX()+1; j++)

            for (Int_t j{1}; j<hEtaForwardFinalAbsSystUncrt.at(i)->GetNbinsX()+1; j++) {
                Double_t xVal{-999.};
                Double_t yVal{-999.};
                Double_t relUncrt{0.};
                Double_t absUncrt{0.};

                xVal = hEtaForwardFinalAbsSystUncrt.at(i)->GetBinCenter(j);
                yVal = hEtaForwardFinal.at(i)->GetBinContent(j);
                relUncrt = hEtaForwardFinalRelSystUncrt.at(i)->GetBinContent(j); 
                absUncrt = hEtaForwardFinal.at(i)->GetBinContent(j) * relUncrt;
                hEtaForwardFinalAbsSystUncrt.at(i)->SetBinContent(j, yVal);
                hEtaForwardFinalAbsSystUncrt.at(i)->SetBinError(j, absUncrt);

                xVal = hEtaBackwardFinalAbsSystUncrt.at(i)->GetBinCenter(j);
                yVal = hEtaBackwardFinal.at(i)->GetBinContent(j);
                relUncrt = hEtaBackwardFinalRelSystUncrt.at(i)->GetBinContent(j);
                absUncrt = hEtaBackwardFinal.at(i)->GetBinContent(j) * relUncrt;
                hEtaBackwardFinalAbsSystUncrt.at(i)->SetBinContent(j, yVal);
                hEtaBackwardFinalAbsSystUncrt.at(i)->SetBinError(j, absUncrt);

                xVal = hEtaFBRatioFinalAbsSystUncrt.at(i)->GetBinCenter(j);
                yVal = hEtaFBRatioFinal.at(i)->GetBinContent(j);
                relUncrt = hEtaFBRatioFinalRelSystUncrt.at(i)->GetBinContent(j);
                absUncrt = hEtaFBRatioFinal.at(i)->GetBinContent(j) * relUncrt;
                hEtaFBRatioFinalAbsSystUncrt.at(i)->SetBinContent(j, yVal);
                hEtaFBRatioFinalAbsSystUncrt.at(i)->SetBinError(j, absUncrt);

                xVal = hEtaBFRatioFinalAbsSystUncrt.at(i)->GetBinCenter(j);
                yVal = hEtaBFRatioFinal.at(i)->GetBinContent(j);
                relUncrt = hEtaBFRatioFinalRelSystUncrt.at(i)->GetBinContent(j);
                absUncrt = hEtaBFRatioFinal.at(i)->GetBinContent(j) * relUncrt;
                hEtaBFRatioFinalAbsSystUncrt.at(i)->SetBinContent(j, yVal);
                hEtaBFRatioFinalAbsSystUncrt.at(i)->SetBinError(j, absUncrt);
            } //  for (Int_t j{1}; j< hEtaForwardFinalAbsSystUncrt.at(i)->GetNbinsX()+1; j++)

    } // for (Int_t i{0}; i<ptBins; i++)

    hFinalDistributions.push_back(hEtaFinalAbsSystUncrt);                       // 0
    hFinalDistributions.push_back(hEtaForwardFinalAbsSystUncrt);                // 1
    hFinalDistributions.push_back(hEtaBackwardFinalAbsSystUncrt);               // 2
    hFinalDistributions.push_back(hEtaFBRatioFinalAbsSystUncrt);                // 3
    hFinalDistributions.push_back(hEtaBFRatioFinalAbsSystUncrt);                // 4

    return hFinalDistributions;
}

//________________
void plotIndividualDataOverNpdf(TCanvas *c, TH1D* data, TH1D* epps21Stat, TH1D* epps21Syst, 
                                TH1D* ncteq15hqStat, TH1D* ncteq15hqSyst,
                                Int_t etaType=0, Int_t ptLow=50, Int_t ptHi=60) {

    // Text 
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Number for plotting position
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.7, 1.3};
    Double_t legX[2] = {0.4, 0.65};
    Double_t legY[2] = {0.2, 0.35};

    // etaType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    if ( etaType > 0 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
    }

    // Set style for the data points
    Int_t dataStyle{0};
    Int_t epps21Style{1};
    Int_t ncteq15hqStyle{2};
    setSystUncrtStyle(data, dataStyle);
    setSystUncrtStyle(epps21Syst, epps21Style);
    setSystUncrtStyle(ncteq15hqSyst, ncteq15hqStyle);
    set1DStyle(epps21Stat, 0);
    set1DStyle(ncteq15hqStat, 1);

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot comparison
    //
    c->cd();

    // Set pad style
    if ( etaType == 0 ) {
        setPadStyle();
    }
    else {
        setPadFBStyle();
    }

    // Plot distributions
    data->Draw("E2");
    epps21Syst->Draw("E2 same");
    ncteq15hqSyst->Draw("E2 same");
    epps21Stat->Draw("same");
    ncteq15hqStat->Draw("same");

    // Set ranges
    data->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    data->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    data->GetYaxis()->SetTitle("Data / nPDF");

    t.DrawLatexNDC(0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
    plotCMSHeader();

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(data, "Data uncrt.", "f");
    leg->AddEntry(epps21Syst, "EPPS21 x CT18", "pf");
    leg->AddEntry(ncteq15hqSyst, "nCTEQ15HQ", "pf");
    leg->Draw();

    // Line at unity
    line = new TLine(xRange[0], 1.0, xRange[1], 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kMagenta);
    line->Draw();
}

//________________
void plotDataOverNpdf(std::vector< std::vector<TH1D*> > epps21Dist, 
                      std::vector< std::vector<TH1D*> > ncteq15hqDist, 
                      std::vector< std::vector<TH1D*> > dataRelSystUncrt,
                      Int_t etaType, TString date, Bool_t isCM) {

    // Dijet pT standard binning
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    std::vector< TH1D* > dataUncrtAtUnity = dataRelSystUncrt.at(5);
    std::vector< TH1D* > epps21StatUncrt = epps21Dist.at(0);
    std::vector< TH1D* > epps21SystUncrt = epps21Dist.at(1);
    std::vector< TH1D* > ncteq15hqStatUncrt = ncteq15hqDist.at(0);
    std::vector< TH1D* > ncteq15hqSystUncrt = ncteq15hqDist.at(1);

    Int_t ptBin = -1;

    // Text 
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Number for plotting position
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.7, 1.3};
    Double_t legX[2] = {0.4, 0.65};
    Double_t legY[2] = {0.2, 0.35};

    // etaType: 0 - all, 1 - forward, 2 - backward, 3 - forward/backward, 4 - backward/forward
    if ( etaType > 0 ) {     // Forward or backward
        xRange[0] = 0.;
        xRange[1] = 3.0;
    }

    // Set style for the data points
    Int_t dataStyle{0};
    Int_t epps21Style{1};
    Int_t ncteq15hqStyle{2};

    // Canvases
    Int_t sizeX{1200}, sizeY{800};
    TCanvas *c = new TCanvas("c", "c", sizeX, sizeY);

    TCanvas *comp = new TCanvas("comp", "comp", 1500, 300 );
    comp->Divide(npdfPtLow.size(), 1);

    TLegend *leg; 

    for (Int_t i{0}; i<epps21StatUncrt.size(); i++) {

        // Find the proper pT bin
        for (Int_t j{0}; j<ptDijetLow.size(); j++) {
            if ( npdfPtLow.at(i) == ptDijetLow.at(j) && npdfPtHi.at(i) == ptDijetHi.at(j) ) {
                ptBin = j;
                break;
            }
        }

        // Switch to canvas
        plotIndividualDataOverNpdf(c, dataUncrtAtUnity.at(ptBin), epps21StatUncrt.at(i), epps21SystUncrt.at(i), 
                                   ncteq15hqStatUncrt.at(i), ncteq15hqSystUncrt.at(i),
                                   etaType, ptDijetLow.at(ptBin), ptDijetHi.at(ptBin));
        c->SaveAs( Form("%s/data2mc/comp2npdf_eta_pt_%d_%d_%s.png", date.Data(), ptDijetLow.at(ptBin), ptDijetHi.at(ptBin), date.Data()) );

        // Comparison on a single canvas
        comp->cd(i+1);
        if ( etaType == 0 ) {
            setPadStyle();
        }
        else {
            setPadFBStyle();
        }

        // Set histogram styles
        setSystUncrtStyle(dataUncrtAtUnity.at(ptBin), dataStyle);
        setSystUncrtStyle(epps21SystUncrt.at(i), epps21Style);
        setSystUncrtStyle(ncteq15hqSystUncrt.at(i), ncteq15hqStyle);
        set1DStyle(epps21StatUncrt.at(i), 0);
        set1DStyle(ncteq15hqStatUncrt.at(i), 1);

        // Draw histograms
        dataUncrtAtUnity.at(ptBin)->Draw("E2");
        epps21SystUncrt.at(i)->Draw("E2 same");
        ncteq15hqSystUncrt.at(i)->Draw("E2 same");
        epps21StatUncrt.at(i)->Draw("same");
        ncteq15hqStatUncrt.at(i)->Draw("same");
        dataUncrtAtUnity.at(ptBin)->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
        dataUncrtAtUnity.at(ptBin)->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
        dataUncrtAtUnity.at(ptBin)->GetYaxis()->SetTitle("Data / nPDF");

        // Plot pT label
        t.DrawLatexNDC(0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", ptDijetLow.at(ptBin), ptDijetHi.at(ptBin)));
        plotCMSHeader();

        // Legend
        leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        leg->AddEntry(dataUncrtAtUnity.at(ptBin), "Data uncrt.", "f");
        leg->AddEntry(epps21SystUncrt.at(i), "EPPS21 x CT18", "pf");
        leg->AddEntry(ncteq15hqSystUncrt.at(i), "nCTEQ15HQ", "pf");
        leg->Draw();

        // Line at unity
        TLine *line = new TLine(xRange[0], 1.0, xRange[1], 1.0);
        line->SetLineStyle(2);
        line->SetLineColor(kMagenta);
        line->Draw();
    } // for (Int_t i{0}; i<epps21Dist.size(); i++)
    comp->SaveAs( Form("%s/data2mc/comp2npdf_eta_%s.png", date.Data(), date.Data()) );
}

//________________
void retrieveDistributions(TFile *mbFile, TFile *mbPbGoingFile, TFile *mbPGoingFile, 
                           TFile *jet60File, TFile *jet60PbGoingFile, TFile *jet60PGoingFile,
                           TFile *jet80File, TFile *jet80PbGoingFile, TFile *jet80PGoingFile,
                           TFile *jet100File, TFile *jet100PbGoingFile, TFile *jet100PGoingFile,
                           TFile *jeuMBUpFile, TFile *jeuMBDownFile,
                           TFile *jeuJet60UpFile, TFile *jeuJet60DownFile,
                           TFile *jeuJet80UpFile, TFile *jeuJet80DownFile,
                           TFile *jeuJet100UpFile, TFile *jeuJet100DownFile,
                           TFile *embeddingFile, TFile *jerUpFile, TFile *jerDownFile,
                           TFile *pbGoingEmbeddingFile, TFile *pGoingEmbeddingFile,
                           TFile *mbGplusFile, TFile *mbVtx1File,
                           TFile *jet60GplusFile, TFile *jet60Vtx1File,
                           TFile *jet80GplusFile, TFile *jet80Vtx1File,
                           TFile *jet100GplusFile, TFile *jet100Vtx1File,
                           TFile *epps21File, TFile *ncteq15hqFile,
                           TString date, Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";
    TString histoName = "hRecoDijetPtEta";
    // if ( isCM ) {
    //     histoName = "hRecoDijetPtEta";
    // }

    TString genHistoName = "hGenDijetPtEta";
    // if ( isCM ) {
    //     genHistoName = "hGenDijetPtEta";
    // }

    TString refHistoName = "hRefDijetPtEta";
    // if ( isCM ) {
    //     refHistoName = "hRefDijetPtEta";
    // }

    // Output histogram name
    TString oHistoName;

    // Dijet pT selection
    // Int_t ptBins = ptDijetBinLow.size();
    // std::vector<Int_t> ptDijetLow{};
    // std::vector<Int_t> ptDijetHi{};
    // fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};
    Int_t genType{4};

    // TLatex t;
    // t.SetTextFont(42);
    // t.SetTextSize(0.06);

    //
    // Retrieve dijet distributions
    //

    // Data

    // MB
    oHistoName = returnTrigName(mbFile) + "_pPb8160_etaDijet";
    std::vector< std::vector<TH1D*> > hMBEtaDist = createDijetEtaHistograms(mbFile, histoName.Data(), oHistoName.Data(), isCM, defType);
    // MB Pb-going
    oHistoName = returnTrigName(mbPbGoingFile) + "_pPb8160_pbGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hMBPbGoingEtaDist = createDijetEtaHistograms(mbPbGoingFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // MB P-going
    oHistoName = returnTrigName(mbPGoingFile) + "_pPb8160_pGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hMBPGoingEtaDist = createDijetEtaHistograms(mbPGoingFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Jet60
    oHistoName = returnTrigName(jet60File) + "_pPb8160_etaDijet";
    std::vector< std::vector<TH1D*> > hJet60EtaDist = createDijetEtaHistograms(jet60File, histoName.Data(), oHistoName.Data(), isCM, defType);
    // Jet60 Pb-going
    oHistoName = returnTrigName(jet60PbGoingFile) + "_pPb8160_pbGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hJet60PbGoingEtaDist = createDijetEtaHistograms(jet60PbGoingFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Jet60 P-going
    oHistoName = returnTrigName(jet60PGoingFile) + "_pPb8160_pGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hJet60PGoingEtaDist = createDijetEtaHistograms(jet60PGoingFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Jet80
    oHistoName = returnTrigName(jet80File) + "_pPb8160_etaDijet";
    std::vector< std::vector<TH1D*> > hJet80EtaDist = createDijetEtaHistograms(jet80File, histoName.Data(), oHistoName.Data(), isCM, defType);
    // Jet80 Pb-going
    oHistoName = returnTrigName(jet80PbGoingFile) + "_pPb8160_pbGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hJet80PbGoingEtaDist = createDijetEtaHistograms(jet80PbGoingFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Jet80 P-going
    oHistoName = returnTrigName(jet80PGoingFile) + "_pPb8160_pGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hJet80PGoingEtaDist = createDijetEtaHistograms(jet80PGoingFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Jet100
    oHistoName = returnTrigName(jet100File) + "_pPb8160_etaDijet";
    std::vector< std::vector<TH1D*> > hJet100EtaDist = createDijetEtaHistograms(jet100File, histoName.Data(), oHistoName.Data(), isCM, defType);
    // Jet100 Pb-going
    oHistoName = returnTrigName(jet100PbGoingFile) + "_pPb8160_pbGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hJet100PbGoingEtaDist = createDijetEtaHistograms(jet100PbGoingFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Jet100 P-going
    oHistoName = returnTrigName(jet100PGoingFile) + "_pPb8160_pGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hJet100PGoingEtaDist = createDijetEtaHistograms(jet100PGoingFile, histoName.Data(), oHistoName.Data(), isCM, downType);

    // JEU

    // JEU MB Up
    oHistoName = returnTrigName(jeuMBUpFile) + "_pPb8160_etaDijet_jeuUp";
    std::vector< std::vector<TH1D*> > hJeuMBUpEtaDist = createDijetEtaHistograms(jeuMBUpFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // JEU MB Down
    oHistoName = returnTrigName(jeuMBDownFile) + "_pPb8160_etaDijet_jeuDown";
    std::vector< std::vector<TH1D*> > hJeuMBDownEtaDist = createDijetEtaHistograms(jeuMBDownFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // JEU Jet60 Up
    oHistoName = returnTrigName(jeuJet60UpFile) + "_pPb8160_etaDijet_jeuUp";
    std::vector< std::vector<TH1D*> > hJeuJet60UpEtaDist = createDijetEtaHistograms(jeuJet60UpFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // JEU Jet60 Down
    oHistoName = returnTrigName(jeuJet60DownFile) + "_pPb8160_etaDijet_jeuDown";
    std::vector< std::vector<TH1D*> > hJeuJet60DownEtaDist = createDijetEtaHistograms(jeuJet60DownFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // JEU Jet80 Up
    oHistoName = returnTrigName(jeuJet80UpFile) + "_pPb8160_etaDijet_jeuUp";
    std::vector< std::vector<TH1D*> > hJeuJet80UpEtaDist = createDijetEtaHistograms(jeuJet80UpFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // JEU Jet80 Down
    oHistoName = returnTrigName(jeuJet80DownFile) + "_pPb8160_etaDijet_jeuDown";
    std::vector< std::vector<TH1D*> > hJeuJet80DownEtaDist = createDijetEtaHistograms(jeuJet80DownFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // JEU Jet100 Up
    oHistoName = returnTrigName(jeuJet100UpFile) + "_pPb8160_etaDijet_jeuUp";
    std::vector< std::vector<TH1D*> > hJeuJet100UpEtaDist = createDijetEtaHistograms(jeuJet100UpFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // JEU Jet100 Down
    oHistoName = returnTrigName(jeuJet100DownFile) + "_pPb8160_etaDijet_jeuDown";
    std::vector< std::vector<TH1D*> > hJeuJet100DownEtaDist = createDijetEtaHistograms(jeuJet100DownFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    
    // Monte Carlo

    // Embedding
    oHistoName = returnTrigName(embeddingFile) + "_pPb8160_etaDijet";
    std::vector< std::vector<TH1D*> > hEmbeddingEtaDist = createDijetEtaHistograms(embeddingFile, histoName.Data(), oHistoName.Data(), isCM, defType);
    oHistoName = returnTrigName(embeddingFile) + "_pPb8160_etaDijet_ref";
    std::vector< std::vector<TH1D*> > hEmbeddingRefEtaDist = createDijetEtaHistograms(embeddingFile, refHistoName.Data(), oHistoName.Data(), isCM, genType);
    oHistoName = returnTrigName(embeddingFile) + "_pPb8160_etaDijet_gen";
    std::vector< std::vector<TH1D*> > hEmbeddingGenEtaDist = createDijetEtaHistograms(embeddingFile, genHistoName.Data(), oHistoName.Data(), isCM, genType);
    // JER Up
    oHistoName = returnTrigName(jerUpFile) + "_pPb8160_etaDijet_jerUp";
    std::vector< std::vector<TH1D*> > hJerUpEtaDist = createDijetEtaHistograms(jerUpFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // JER Down
    oHistoName = returnTrigName(jerDownFile) + "_pPb8160_etaDijet_jerDown";
    std::vector< std::vector<TH1D*> > hJerDownEtaDist = createDijetEtaHistograms(jerDownFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Embedding Pb-going
    oHistoName = returnTrigName(pbGoingEmbeddingFile) + "_pPb8160_pbGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hPbGoingEmbeddingRecoEtaDist = createDijetEtaHistograms(pbGoingEmbeddingFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Embedding P-going
    oHistoName = returnTrigName(pGoingEmbeddingFile) + "_pPb8160_pGoing_etaDijet";
    std::vector< std::vector<TH1D*> > hPGoingEmbeddingRecoEtaDist = createDijetEtaHistograms(pGoingEmbeddingFile, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Embedding Pb-going
    oHistoName = returnTrigName(pbGoingEmbeddingFile) + "_pPb8160_pbGoing_etaDijet_gen";
    std::vector< std::vector<TH1D*> > hPbGoingEmbeddingGenEtaDist = createDijetEtaHistograms(pbGoingEmbeddingFile, genHistoName.Data(), oHistoName.Data(), isCM, upType);
    // Embedding P-going
    oHistoName = returnTrigName(pGoingEmbeddingFile) + "_pPb8160_pGoing_etaDijet_gen";
    std::vector< std::vector<TH1D*> > hPGoingEmbeddingGenEtaDist = createDijetEtaHistograms(pGoingEmbeddingFile, genHistoName.Data(), oHistoName.Data(), isCM, downType);

    // Pileup

    // MB Gplus
    oHistoName = returnTrigName(mbGplusFile) + "_pPb8160_etaDijet_gplus";
    std::vector< std::vector<TH1D*> > hMBGplusEtaDist = createDijetEtaHistograms(mbGplusFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // MB vtx1
    oHistoName = returnTrigName(mbVtx1File) + "_pPb8160_etaDijet_vtx1";
    std::vector< std::vector<TH1D*> > hMBVtx1EtaDist = createDijetEtaHistograms(mbVtx1File, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Jet60 Gplus
    oHistoName = returnTrigName(jet60GplusFile) + "_pPb8160_etaDijet_gplus";
    std::vector< std::vector<TH1D*> > hJet60GplusEtaDist = createDijetEtaHistograms(jet60GplusFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Jet60 vtx1
    oHistoName = returnTrigName(jet60Vtx1File) + "_pPb8160_etaDijet_vtx1";
    std::vector< std::vector<TH1D*> > hJet60Vtx1EtaDist = createDijetEtaHistograms(jet60Vtx1File, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Jet80 Gplus
    oHistoName = returnTrigName(jet80GplusFile) + "_pPb8160_etaDijet_gplus";
    std::vector< std::vector<TH1D*> > hJet80GplusEtaDist = createDijetEtaHistograms(jet80GplusFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Jet80 vtx1
    oHistoName = returnTrigName(jet80Vtx1File) + "_pPb8160_etaDijet_vtx1";
    std::vector< std::vector<TH1D*> > hJet80Vtx1EtaDist = createDijetEtaHistograms(jet80Vtx1File, histoName.Data(), oHistoName.Data(), isCM, downType);
    // Jet100 Gplus
    oHistoName = returnTrigName(jet100GplusFile) + "_pPb8160_etaDijet_gplus";
    std::vector< std::vector<TH1D*> > hJet100GplusEtaDist = createDijetEtaHistograms(jet100GplusFile, histoName.Data(), oHistoName.Data(), isCM, upType);
    // Jet100 vtx1
    oHistoName = returnTrigName(jet100Vtx1File) + "_pPb8160_etaDijet_vtx1";
    std::vector< std::vector<TH1D*> > hJet100Vtx1EtaDist = createDijetEtaHistograms(jet100Vtx1File, histoName.Data(), oHistoName.Data(), isCM, downType);

    // nPDF
    std::vector< std::vector<TH1D*> > hEpps21Dist = retrieveNLOJetDistributions(epps21File, 0);
    std::vector< std::vector<TH1D*> > hNcteq15hqDist = retrieveNLOJetDistributions(ncteq15hqFile, 1);

    // Create ratios to MC
    std::vector< std::vector<TH1D*> > hRatios2McDist = createRatios2MC(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist, hEmbeddingEtaDist, hEmbeddingGenEtaDist, isCM);

    // Create up and down variation to default ratios
    std::vector< std::vector<TH1D*> > hUpAndDown2DefDist = createUpAndDown2Def(hMBEtaDist, hJeuMBUpEtaDist, hJeuMBDownEtaDist, 
                                                                               hJet60EtaDist, hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, 
                                                                               hJet80EtaDist, hJeuJet80UpEtaDist, hJeuJet80DownEtaDist, 
                                                                               hJet100EtaDist, hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, 
                                                                               hEmbeddingEtaDist, hEmbeddingRefEtaDist, hJerUpEtaDist, hJerDownEtaDist,
                                                                               hMBGplusEtaDist, hMBVtx1EtaDist,
                                                                               hJet60GplusEtaDist, hJet60Vtx1EtaDist,
                                                                               hJet80GplusEtaDist, hJet80Vtx1EtaDist,
                                                                               hJet100GplusEtaDist, hJet100Vtx1EtaDist, 
                                                                               isCM, drawFits);
    
    
    // Create JEU up and down histograms
    std::vector< std::vector< TH1D* > > hMBEtaJeuUp2DefDist; hMBEtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(0) ); hMBEtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(1) ); hMBEtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(2) ); hMBEtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(3) ); hMBEtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(4) );
    std::vector< std::vector< TH1D* > > hMBEtaJeuDown2DefDist; hMBEtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(5) ); hMBEtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(6) ); hMBEtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(7) ); hMBEtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(8) ); hMBEtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(9) );
    std::vector< std::vector< TH1D* > > hJet60EtaJeuUp2DefDist; hJet60EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(10) ); hJet60EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(11) ); hJet60EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(12) ); hJet60EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(13) ); hJet60EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(14) );
    std::vector< std::vector< TH1D* > > hJet60EtaJeuDown2DefDist; hJet60EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(15) ); hJet60EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(16) ); hJet60EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(17) ); hJet60EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(18) ); hJet60EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(19) );
    std::vector< std::vector< TH1D* > > hJet80EtaJeuUp2DefDist; hJet80EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(20) ); hJet80EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(21) ); hJet80EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(22) ); hJet80EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(23) ); hJet80EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(24) );
    std::vector< std::vector< TH1D* > > hJet80EtaJeuDown2DefDist; hJet80EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(25) ); hJet80EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(26) ); hJet80EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(27) ); hJet80EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(28) ); hJet80EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(29) );
    std::vector< std::vector< TH1D* > > hJet100EtaJeuUp2DefDist; hJet100EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(30) ); hJet100EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(31) ); hJet100EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(32) ); hJet100EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(33) ); hJet100EtaJeuUp2DefDist.push_back( hUpAndDown2DefDist.at(34) );
    std::vector< std::vector< TH1D* > > hJet100EtaJeuDown2DefDist; hJet100EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(35) ); hJet100EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(36) ); hJet100EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(37) ); hJet100EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(38) ); hJet100EtaJeuDown2DefDist.push_back( hUpAndDown2DefDist.at(39) );

    // Create JER up and down histograms
    std::vector< std::vector<TH1D*> > hEmbeddingJerUp2DefDist; hEmbeddingJerUp2DefDist.push_back( hUpAndDown2DefDist.at(40) ); hEmbeddingJerUp2DefDist.push_back( hUpAndDown2DefDist.at(41) ); hEmbeddingJerUp2DefDist.push_back( hUpAndDown2DefDist.at(42) ); hEmbeddingJerUp2DefDist.push_back( hUpAndDown2DefDist.at(43) ); hEmbeddingJerUp2DefDist.push_back( hUpAndDown2DefDist.at(44) );
    std::vector< std::vector<TH1D*> > hEmbeddingJerDown2DefDist; hEmbeddingJerDown2DefDist.push_back( hUpAndDown2DefDist.at(45) ); hEmbeddingJerDown2DefDist.push_back( hUpAndDown2DefDist.at(46) ); hEmbeddingJerDown2DefDist.push_back( hUpAndDown2DefDist.at(47) ); hEmbeddingJerDown2DefDist.push_back( hUpAndDown2DefDist.at(48) ); hEmbeddingJerDown2DefDist.push_back( hUpAndDown2DefDist.at(49) );

    // Create pileup up and down histograms
    std::vector< std::vector<TH1D*> > hMBPileupGplus2DefDist; hMBPileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(50) ); hMBPileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(51) ); hMBPileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(52) ); hMBPileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(53) ); hMBPileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(54) );
    std::vector< std::vector<TH1D*> > hMBPileupVtxOne2DefDist; hMBPileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(55) ); hMBPileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(56) ); hMBPileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(57) ); hMBPileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(58) ); hMBPileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(59) );
    std::vector< std::vector<TH1D*> > hJet60PileupGplus2DefDist; hJet60PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(60) ); hJet60PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(61) ); hJet60PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(62) ); hJet60PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(63) ); hJet60PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(64) );
    std::vector< std::vector<TH1D*> > hJet60PileupVtxOne2DefDist; hJet60PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(65) ); hJet60PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(66) ); hJet60PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(67) ); hJet60PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(68) ); hJet60PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(69) );
    std::vector< std::vector<TH1D*> > hJet80PileupGplus2DefDist; hJet80PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(70) ); hJet80PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(71) ); hJet80PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(72) ); hJet80PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(73) ); hJet80PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(74) );
    std::vector< std::vector<TH1D*> > hJet80PileupVtxOne2DefDist; hJet80PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(75) ); hJet80PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(76) ); hJet80PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(77) ); hJet80PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(78) ); hJet80PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(79) );
    std::vector< std::vector<TH1D*> > hJet100PileupGplus2DefDist; hJet100PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(80) ); hJet100PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(81) ); hJet100PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(82) ); hJet100PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(83) ); hJet100PileupGplus2DefDist.push_back( hUpAndDown2DefDist.at(84) );
    std::vector< std::vector<TH1D*> > hJet100PileupVtxOne2DefDist; hJet100PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(85) ); hJet100PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(86) ); hJet100PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(87) ); hJet100PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(88) ); hJet100PileupVtxOne2DefDist.push_back( hUpAndDown2DefDist.at(89) );

    // Create pointing resolution ratios
    std::vector< std::vector<TH1D*> > hPointingRes2DefDist; hPointingRes2DefDist.push_back( hUpAndDown2DefDist.at(90) ); hPointingRes2DefDist.push_back( hUpAndDown2DefDist.at(91) ); hPointingRes2DefDist.push_back( hUpAndDown2DefDist.at(92) ); hPointingRes2DefDist.push_back( hUpAndDown2DefDist.at(93) ); hPointingRes2DefDist.push_back( hUpAndDown2DefDist.at(94) );

    // Systematics type
    Int_t jeuType{0};
    Int_t jerType{1};
    Int_t pointingType{2};
    Int_t pileupType{3};

    Bool_t useFitJeu = kTRUE;
    Bool_t useFitJer = kFALSE;
    Bool_t useFitPointing = kTRUE;
    Bool_t useFitPileup = kTRUE;

    // Create JEU systematic tables
    std::vector< std::vector<TH1D*> > hMBJeuRelSystDist = createRelSystUncrt(hMBEtaJeuUp2DefDist, hMBEtaJeuDown2DefDist, "hMB", jeuType, useFitJeu);
    std::vector< std::vector<TH1D*> > hJet60JeuRelSystDist = createRelSystUncrt(hJet60EtaJeuUp2DefDist, hJet60EtaJeuDown2DefDist, "hJet60", jeuType, useFitJeu);
    std::vector< std::vector<TH1D*> > hJet80JeuRelSystDist = createRelSystUncrt(hJet80EtaJeuUp2DefDist, hJet80EtaJeuDown2DefDist, "hJet80", jeuType, useFitJeu);
    std::vector< std::vector<TH1D*> > hJet100JeuRelSystDist = createRelSystUncrt(hJet100EtaJeuUp2DefDist, hJet100EtaJeuDown2DefDist, "hJet100", jeuType, useFitJeu);

    // Create JER systematic tables
    std::vector< std::vector<TH1D*> > hEmbeddingJerRelSystDist = createRelSystUncrt(hEmbeddingJerUp2DefDist, hEmbeddingJerDown2DefDist, "hEmbedding", jerType, useFitJer);
    // Create pointing resolution systematic tables
    std::vector< std::vector<TH1D*> > hEmpty{};
    std::vector< std::vector<TH1D*> > hPointingResRelSystDist = createRelSystUncrt(hPointingRes2DefDist, hEmpty, "hPointingRes", pointingType, useFitPointing);

    // Create pileup systematic tables
    std::vector< std::vector<TH1D*> > hMBPileupRelSystDist = createRelSystUncrt(hMBPileupGplus2DefDist, hEmpty, "hMB", pileupType, useFitPileup);
    std::vector< std::vector<TH1D*> > hJet60PileupRelSystDist = createRelSystUncrt(hJet60PileupGplus2DefDist, hEmpty, "hJet60", pileupType, useFitPileup);
    std::vector< std::vector<TH1D*> > hJet80PileupRelSystDist = createRelSystUncrt(hJet80PileupGplus2DefDist, hEmpty, "hJet80", pileupType, useFitPileup);
    std::vector< std::vector<TH1D*> > hJet100PileupRelSystDist = createRelSystUncrt(hJet100PileupGplus2DefDist, hEmpty, "hJet100", pileupType, useFitPileup);

    // std::vector< std::vector<TH1D*> > hMBPileupRelSystDist = createRelSystUncrt(hMBPileupGplus2DefDist, hMBPileupVtxOne2DefDist, "hMB", pileupType);
    // std::vector< std::vector<TH1D*> > hJet60PileupRelSystDist = createRelSystUncrt(hJet60PileupGplus2DefDist, hJet60PileupVtxOne2DefDist, "hJet60", pileupType);
    // std::vector< std::vector<TH1D*> > hJet80PileupRelSystDist = createRelSystUncrt(hJet80PileupGplus2DefDist, hJet80PileupVtxOne2DefDist, "hJet80", pileupType);
    // std::vector< std::vector<TH1D*> > hJet100PileupRelSystDist = createRelSystUncrt(hJet100PileupGplus2DefDist, hJet100PileupVtxOne2DefDist, "hJet100", pileupType);

    // Calculate total systematic uncertainty
    std::vector< std::vector<TH1D*> > hMBTotalRelSystDist = createTotalRelSystUncrt(hMBJeuRelSystDist, hEmbeddingJerRelSystDist, hPointingResRelSystDist, hMBPileupRelSystDist, "hMB");
    std::vector< std::vector<TH1D*> > hJet60TotalRelSystDist = createTotalRelSystUncrt(hJet60JeuRelSystDist, hEmbeddingJerRelSystDist, hPointingResRelSystDist, hJet60PileupRelSystDist, "hJet60");
    std::vector< std::vector<TH1D*> > hJet80TotalRelSystDist = createTotalRelSystUncrt(hJet80JeuRelSystDist, hEmbeddingJerRelSystDist, hPointingResRelSystDist, hJet80PileupRelSystDist, "hJet80");
    std::vector< std::vector<TH1D*> > hJet100TotalRelSystDist = createTotalRelSystUncrt(hJet100JeuRelSystDist, hEmbeddingJerRelSystDist, hPointingResRelSystDist, hJet100PileupRelSystDist, "hJet100");

    // Make final distributions
    std::vector< std::vector<TH1D*> > hFinalDist = makeFinalDistributions(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist, 
                                                                          hEmbeddingGenEtaDist, hEmbeddingEtaDist, hRatios2McDist);

    Int_t epps21type{0};
    Int_t ncteq15hqtype{1};
    std::vector< std::vector< TH1D* > > hFinalOverEpps21Dist = makeRatiosOfData2nPDF(hFinalDist, hEpps21Dist, epps21type);
    std::vector< std::vector< TH1D* > > hFinalOverNcteq15hqDist = makeRatiosOfData2nPDF(hFinalDist, hNcteq15hqDist, ncteq15hqtype);

    // Make final systematic uncertainty
    std::vector< std::vector<TH1D*> > hFinalRelSystUncrtDist = makeFinalRelSystUncrt(hMBTotalRelSystDist, hJet60TotalRelSystDist, hJet80TotalRelSystDist, hJet100TotalRelSystDist); 


    // Make final absolute systematic uncertainty
    std::vector< std::vector<TH1D*> > hFinalAbsSystUncrtDist = makeFinalAbsSystUncrt(hFinalDist, hFinalRelSystUncrtDist);

    // Plot final distributions
    plotFinalEtaDistributions(hFinalDist, hFinalAbsSystUncrtDist, date, isCM);

    // Plot comparison of data and nPDF
    Int_t etaType{0};
    // Int_t etaForwardType{1};
    // plotDataOverNpdf(hFinalOverEpps21Dist, hFinalOverNcteq15hqDist, hFinalRelSystUncrtDist, etaType, date, isCM);

    // Plot comparison of data and Monte Carlo
    // plotData2McComparison(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist, hEmbeddingEtaDist, hEmbeddingGenEtaDist, hRatios2McDist, date, isCM);

    plotRelativeSystematicUncertainties(hMBJeuRelSystDist, hJet60JeuRelSystDist, hJet80JeuRelSystDist, hJet100JeuRelSystDist,
                                        hEmbeddingJerRelSystDist, hPointingResRelSystDist,
                                        hMBPileupRelSystDist, hJet60PileupRelSystDist, hJet80PileupRelSystDist, hJet100PileupRelSystDist,
                                        hMBTotalRelSystDist, hJet60TotalRelSystDist, hJet80TotalRelSystDist, hJet100TotalRelSystDist, 
                                        date, isCM);

    // plotJES(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist,
    //         hJeuMBUpEtaDist, hJeuMBDownEtaDist,
    //         hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, 
    //         hJeuJet80UpEtaDist, hJeuJet80DownEtaDist,
    //         hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, 
    //         hMBEtaJeuUp2DefDist, hMBEtaJeuDown2DefDist,
    //         hJet60EtaJeuUp2DefDist, hJet60EtaJeuDown2DefDist,
    //         hJet80EtaJeuUp2DefDist, hJet80EtaJeuDown2DefDist,
    //         hJet100EtaJeuUp2DefDist, hJet100EtaJeuDown2DefDist,
    //         date, isCM, drawFits);

    // plotJER(hEmbeddingEtaDist, hJerUpEtaDist, hJerDownEtaDist, hEmbeddingJerUp2DefDist, hEmbeddingJerDown2DefDist,  date, isCM, drawFits);

    // plotPileup(hMBEtaDist, hMBGplusEtaDist, hMBVtx1EtaDist,
    //            hJet60EtaDist, hJet60GplusEtaDist, hJet60Vtx1EtaDist,
    //            hJet80EtaDist, hJet80GplusEtaDist, hJet80Vtx1EtaDist,
    //            hJet100EtaDist, hJet100GplusEtaDist, hJet100Vtx1EtaDist,
    //            hMBPileupGplus2DefDist, hMBPileupVtxOne2DefDist,
    //            hJet60PileupGplus2DefDist, hJet60PileupVtxOne2DefDist,
    //            hJet80PileupGplus2DefDist, hJet80PileupVtxOne2DefDist,
    //            hJet100PileupGplus2DefDist, hJet100PileupVtxOne2DefDist,
    //            date, isCM, drawFits);

    // plotData2McComparison(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist,
    //                       hEmbeddingEtaDist, date, isCM);

    // Write default distributions to output file
    writeDistributions2RootFile(hMBEtaDist, hMBPbGoingEtaDist, hMBPGoingEtaDist, 
                                hJet60EtaDist, hJet60PbGoingEtaDist, hJet60PGoingEtaDist,
                                hJet80EtaDist, hJet80PbGoingEtaDist, hJet80PGoingEtaDist, 
                                hJet100EtaDist, hJet100PbGoingEtaDist, hJet100PGoingEtaDist,
                                hJeuMBUpEtaDist, hJeuMBDownEtaDist, 
                                hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, 
                                hJeuJet80UpEtaDist, hJeuJet80DownEtaDist,
                                hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, 
                                hEmbeddingEtaDist, hEmbeddingRefEtaDist, hEmbeddingGenEtaDist, 
                                hJerUpEtaDist, hJerDownEtaDist,
                                hPbGoingEmbeddingRecoEtaDist, hPGoingEmbeddingRecoEtaDist,
                                hPbGoingEmbeddingGenEtaDist, hPGoingEmbeddingGenEtaDist,
                                //hEmbeddingPointResDist,
                                hMBGplusEtaDist, hMBVtx1EtaDist, 
                                hJet60GplusEtaDist, hJet60Vtx1EtaDist, 
                                hJet80GplusEtaDist, hJet80Vtx1EtaDist, 
                                hJet100GplusEtaDist, hJet100Vtx1EtaDist, 
                                hRatios2McDist, hUpAndDown2DefDist,
                                hMBJeuRelSystDist, hJet60JeuRelSystDist, hJet80JeuRelSystDist, hJet100JeuRelSystDist,
                                hEmbeddingJerRelSystDist,
                                hMBPileupRelSystDist, hJet60PileupRelSystDist, hJet80PileupRelSystDist, hJet100PileupRelSystDist,
                                hMBTotalRelSystDist, hJet60TotalRelSystDist, hJet80TotalRelSystDist, hJet100TotalRelSystDist,
                                hFinalDist, hFinalRelSystUncrtDist, hFinalAbsSystUncrtDist,
                                hEpps21Dist, hNcteq15hqDist,
                                date, isCM);

}

//________________
void compareJER(TFile *embFile, TFile *jerDefFile, TFile *jerUpFile, TFile *jerDownFile, TString date) {

    // Retrieve JES histogram
    THnSparseD *hJESPtEtaPhiPure = (THnSparseD*)embFile->Get("hJESInclusiveJetPtEtaPhiWeighted");
    hJESPtEtaPhiPure->SetName("hJESPtEtaPhiPure");
    THnSparseD *hJESPtEtaPhiJerDef = (THnSparseD*)jerDefFile->Get("hJESInclusiveJetPtEtaPhiWeighted");
    hJESPtEtaPhiJerDef->SetName("hJESPtEtaPhiJerDef");
    THnSparseD *hJESPtEtaPhiJerUp = (THnSparseD*)jerUpFile->Get("hJESInclusiveJetPtEtaPhiWeighted");
    hJESPtEtaPhiJerUp->SetName("hJESPtEtaPhiJerUp");
    THnSparseD *hJESPtEtaPhiJerDown = (THnSparseD*)jerDownFile->Get("hJESInclusiveJetPtEtaPhiWeighted");
    hJESPtEtaPhiJerDown->SetName("hJESPtEtaPhiJerDown");

    Int_t rebinX{1}, rebinY{1};

    // Styles
    Int_t pureType{2};
    Int_t upType{0};
    Int_t downType{1};
    Int_t defType{5};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    TLegend *leg;

    // Get JES vs ref pT
    hJESPtEtaPhiPure->GetAxis(2)->SetRange(20, 32); // -1.2 < eta < 1.2
    TH2D* hJESvsPtPure = (TH2D*)hJESPtEtaPhiPure->Projection(0, 1);
    hJESvsPtPure->SetName("hJESvsPtPure");
    hJESvsPtPure->RebinX( rebinX );
    hJESvsPtPure->RebinY( rebinY );
    set2DStyle( hJESvsPtPure );

    hJESPtEtaPhiJerDef->GetAxis(2)->SetRange(20, 32); // -1.2 < eta < 1.2
    TH2D* hJESvsPtJerDef = (TH2D*)hJESPtEtaPhiJerDef->Projection(0, 1);
    hJESvsPtJerDef->SetName("hJESvsPtJerDef");
    hJESvsPtJerDef->RebinX( rebinX );
    hJESvsPtJerDef->RebinY( rebinY );
    set2DStyle( hJESvsPtJerDef );

    hJESPtEtaPhiJerUp->GetAxis(2)->SetRange(20, 32); // -1.2 < eta < 1.2
    TH2D* hJESvsPtJerUp = (TH2D*)hJESPtEtaPhiJerUp->Projection(0, 1);
    hJESvsPtJerUp->SetName("hJESvsPtJerUp");
    hJESvsPtJerUp->RebinX( rebinX );
    hJESvsPtJerUp->RebinY( rebinY );
    set2DStyle( hJESvsPtJerUp );

    hJESPtEtaPhiJerDown->GetAxis(2)->SetRange(20, 32); // -1.2 < eta < 1.2
    TH2D* hJESvsPtJerDown = (TH2D*)hJESPtEtaPhiJerDown->Projection(0, 1);
    hJESvsPtJerDown->SetName("hJESvsPtJerDown");
    hJESvsPtJerDown->RebinX( rebinX );
    hJESvsPtJerDown->RebinY( rebinY );
    set2DStyle( hJESvsPtJerDown );

    // Prepare 1D histograms
    TH1D *hJERVsPtPure;
    TH1D *hJERVsPtJerUp;
    TH1D *hJERVsPtJerDown;
    TH1D *hJERVsPtJerDef;

    TH1D *hJESVsPtPure;
    TH1D *hJESVsPtJerUp;
    TH1D *hJESVsPtJerDown;
    TH1D *hJESVsPtJerDef;

    TF1 *fitPure = new TF1("fitPure", "TMath::Sqrt([0] + [1]/x)", 30., 800.);
    fitPure->SetParameters(0.0415552, 0.960013);
    fitPure->SetLineColor(pureType);
    fitPure->SetLineWidth(2);
    TF1 *fitJerDef = new TF1("fitJerDef", "TMath::Sqrt([0] + [1]/x)", 30., 800.);
    fitJerDef->SetParameters(0.0415552, 0.960013);
    fitJerDef->SetLineColor(defType);
    fitJerDef->SetLineWidth(2);
    TF1 *fitJerUp = new TF1("fitJerUp", "TMath::Sqrt([0] + [1]/x)", 30., 800.);
    fitJerUp->SetParameters(0.0415552, 0.960013);
    fitJerUp->SetLineColor(upType);
    fitJerUp->SetLineWidth(2);
    TF1 *fitJerDown = new TF1("fitJerDown", "TMath::Sqrt([0] + [1]/x)", 30., 800.);
    fitJerDown->SetParameters(0.0415552, 0.960013);
    fitJerDown->SetLineColor(downType);
    fitJerDown->SetLineWidth(2);

    // Make slices
    hJESvsPtPure->FitSlicesY();
    hJESvsPtJerDef->FitSlicesY();
    hJESvsPtJerUp->FitSlicesY();
    hJESvsPtJerDown->FitSlicesY();

    // Retrieve JER
    hJERVsPtPure = (TH1D*)gDirectory->Get("hJESvsPtPure_2");
    hJERVsPtPure->SetName("hJERVsPtPure");
    hJERVsPtPure->GetYaxis()->SetTitle("JER");
    set1DStyle(hJERVsPtPure, pureType);

    hJERVsPtJerUp = (TH1D*)gDirectory->Get("hJESvsPtJerUp_2");
    hJERVsPtJerUp->SetName("hJERVsPtJerUp");
    hJERVsPtJerUp->GetYaxis()->SetTitle("JER");
    set1DStyle(hJERVsPtJerUp, upType);

    hJERVsPtJerDown = (TH1D*)gDirectory->Get("hJESvsPtJerDown_2");
    hJERVsPtJerDown->SetName("hJERVsPtJerDown");
    hJERVsPtJerDown->GetYaxis()->SetTitle("JER");
    set1DStyle(hJERVsPtJerDown, downType);

    hJERVsPtJerDef = (TH1D*)gDirectory->Get("hJESvsPtJerDef_2");
    hJERVsPtJerDef->SetName("hJERVsPtJerDef");
    hJERVsPtJerDef->GetYaxis()->SetTitle("JER");
    set1DStyle(hJERVsPtJerDef, defType);

    // Retrieve JES
    hJESVsPtPure = (TH1D*)gDirectory->Get("hJESvsPtPure_1");
    hJESVsPtPure->SetName("hJESVsPtPure");
    hJESVsPtPure->GetYaxis()->SetTitle("JES");
    set1DStyle(hJESVsPtPure, pureType);

    hJESVsPtJerUp = (TH1D*)gDirectory->Get("hJESvsPtJerUp_1");
    hJESVsPtJerUp->SetName("hJESVsPtJerUp");
    hJESVsPtJerUp->GetYaxis()->SetTitle("JES");
    set1DStyle(hJESVsPtJerUp, upType);

    hJESVsPtJerDown = (TH1D*)gDirectory->Get("hJESvsPtJerDown_1");
    hJESVsPtJerDown->SetName("hJESVsPtJerDown");
    hJESVsPtJerDown->GetYaxis()->SetTitle("JES");
    set1DStyle(hJESVsPtJerDown, downType);

    hJESVsPtJerDef = (TH1D*)gDirectory->Get("hJESvsPtJerDef_1");
    hJESVsPtJerDef->SetName("hJESVsPtJerDef");
    hJESVsPtJerDef->GetYaxis()->SetTitle("JES");
    set1DStyle(hJESVsPtJerDef, defType);

    // Set fit line properties
    hJERVsPtPure->Fit("fitPure","MRE0");
    hJERVsPtPure->SetLineColor(kBlack);
    hJERVsPtJerDef->Fit("fitJerDef","MRE0");
    hJERVsPtJerDef->SetLineColor(kMagenta);
    hJERVsPtJerUp->Fit("fitJerUp","MRE0");
    hJERVsPtJerUp->SetLineColor(kRed);
    hJERVsPtJerDown->Fit("fitJerDown","MRE0");
    hJERVsPtJerDown->SetLineColor(kBlue);


    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    canv->cd();
    setPadStyle();

    hJERVsPtPure->Draw();
    fitPure->Draw("same");
    hJERVsPtJerDef->Draw("same");
    fitJerDef->Draw("same");
    hJERVsPtJerUp->Draw("same");
    fitJerUp->Draw("same");
    hJERVsPtJerDown->Draw("same");
    fitJerDown->Draw("same");

    leg = new TLegend(0.4, 0.6, 0.65, 0.8);
    leg->SetTextSize(0.04);
    leg->SetLineWidth(0);
    leg->AddEntry(hJERVsPtPure, Form("JER pure"), "p");
    leg->AddEntry(hJERVsPtJerDef, Form("JER def"), "p");
    leg->AddEntry(hJERVsPtJerUp, Form("JER+"), "p");
    leg->AddEntry(hJERVsPtJerDown, Form("JER-"), "p");
    leg->Draw();

    t.DrawLatexNDC(0.4, 0.82, Form("Fit func: #sqrt{a + b/x}") );
    t.DrawLatexNDC(0.6, 0.75, Form("Pure a: %5.4f b: %4.3f", fitPure->GetParameter(0), fitPure->GetParameter(1)) );
    t.DrawLatexNDC(0.6, 0.7,  Form("Def  a: %5.4f b: %4.3f", fitJerDef->GetParameter(0), fitJerDef->GetParameter(1)) );
    t.DrawLatexNDC(0.6, 0.65, Form("JER+ a: %5.4f b: %4.3f", fitJerUp->GetParameter(0), fitJerUp->GetParameter(1)) );
    t.DrawLatexNDC(0.6, 0.6,  Form("JER- a: %5.4f b: %4.3f", fitJerDown->GetParameter(0), fitJerDown->GetParameter(1)) );

    TCanvas *canv2 = new TCanvas("canv2", "canv2", 1200, 800);
    canv2->cd();
    setPadStyle();

    hJESVsPtPure->Draw();
    hJESVsPtJerUp->Draw("same");
    hJESVsPtJerDown->Draw("same");
    hJESVsPtJerDef->Draw("same");

    leg = new TLegend(0.4, 0.6, 0.65, 0.8);
    leg->SetTextSize(0.04);
    leg->SetLineWidth(0);
    leg->AddEntry(hJESVsPtPure, Form("JES pure"), "p");
    leg->AddEntry(hJESVsPtJerDef, Form("JES def"), "p");
    leg->AddEntry(hJESVsPtJerUp, Form("JES+"), "p");
    leg->AddEntry(hJESVsPtJerDown, Form("JES-"), "p");
    leg->Draw();

    canv->SaveAs( Form("%s/pPb8160_JER.pdf", date.Data() ) );
    canv2->SaveAs( Form("%s/pPb8160_JES.pdf", date.Data() ) );
}

//________________
void runNewCalculations(Bool_t isCM = kFALSE, Bool_t drawFits = kFALSE) {

    TString cmStr;
    if ( isCM ) {
        cmStr = "_cmSel";
    } 
    else {
        cmStr = "";
    }

    // Check the directory for storing figures exists
    TDatime dt;
    TString date { Form( "%d",dt.GetDate() ) };
    if ( directoryExists( date.Data() ) ) {
        //std::cout << "Directory exists." << std::endl;
    } 
    else {
        createDirectory( date.Data() );
    }

    //
    // File names
    //

    TString pathPrefix = "/Users/gnigmat/cernbox";

    // Data file names
    TString mbFileName( Form("%s/ana/pPb8160/exp/MB_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString mbPbGoingFileName( Form("%s/ana/pPb8160/exp/Pbgoing/def/MB_Pbgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString mbPGoingFileName( Form("%s/ana/pPb8160/exp/pgoing/def/MB_pgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60FileName( Form("%s/ana/pPb8160/exp/Jet60_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60PbGoingFileName( Form("%s/ana/pPb8160/exp/Pbgoing/def/Jet60_Pbgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60PGoingFileName( Form("%s/ana/pPb8160/exp/pgoing/def/Jet60_pgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80FileName( Form("%s/ana/pPb8160/exp/Jet80_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80PbGoingFileName( Form("%s/ana/pPb8160/exp/Pbgoing/def/Jet80_Pbgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80PGoingFileName( Form("%s/ana/pPb8160/exp/pgoing/def/Jet80_pgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100FileName( Form("%s/ana/pPb8160/exp/Jet100_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100PbGoingFileName( Form("%s/ana/pPb8160/exp/Pbgoing/def/Jet100_Pbgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100PGoingFileName( Form("%s/ana/pPb8160/exp/pgoing/def/Jet100_pgoing_pPb8160_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );

    // JEU up/down file names
    TString mbJeuUpFileName( Form("%s/ana/pPb8160/exp/MB_pPb8160_jeu_up_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString mbJeuDownFileName( Form("%s/ana/pPb8160/exp/MB_pPb8160_jeu_down_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60JeuUpFileName( Form("%s/ana/pPb8160/exp/Jet60_pPb8160_jeu_up_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60JeuDownFileName( Form("%s/ana/pPb8160/exp/Jet60_pPb8160_jeu_down_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80JeuUpFileName( Form("%s/ana/pPb8160/exp/Jet80_pPb8160_jeu_up_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80JeuDownFileName( Form("%s/ana/pPb8160/exp/Jet80_pPb8160_jeu_down_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100JeuUpFileName( Form("%s/ana/pPb8160/exp/Jet100_pPb8160_jeu_up_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100JeuDownFileName( Form("%s/ana/pPb8160/exp/Jet100_pPb8160_jeu_down_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );

    // Embedding file names
    TString embeddingFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    // TString embeddingFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_weightMB_ak4.root", pathPrefix.Data()) );
    TString jerUpFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerUp_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jerDownFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDown_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString pbGoingEmbeddingFileName( Form("%s/ana/pPb8160/embedding/Pbgoing/jer/oEmbedding_pPb8160_Pbgoing_jerDef_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString pGoingEmbeddingFileName( Form("%s/ana/pPb8160/embedding/pgoing/jer/oEmbedding_pPb8160_pgoing_jerDef_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );

    // Pileup file names
    TString mbGplusFileName( Form("%s/ana/pPb8160/exp/MB_pPb8160_gplus_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString mbVtx1FileName( Form("%s/ana/pPb8160/exp/MB_pPb8160_vtx1_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60GplusFileName( Form("%s/ana/pPb8160/exp/Jet60_pPb8160_gplus_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet60Vtx1FileName( Form("%s/ana/pPb8160/exp/Jet60_pPb8160_vtx1_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80GplusFileName( Form("%s/ana/pPb8160/exp/Jet80_pPb8160_gplus_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet80Vtx1FileName( Form("%s/ana/pPb8160/exp/Jet80_pPb8160_vtx1_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100GplusFileName( Form("%s/ana/pPb8160/exp/Jet100_pPb8160_gplus_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jet100Vtx1FileName( Form("%s/ana/pPb8160/exp/Jet100_pPb8160_vtx1_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );

    // NLOJet++ file names
    TString epps21FileName( Form("npdf/epps21_pPb8160.root") );
    TString ncteq15hqFileName( Form("npdf/ncteq15hq_pPb8160.root") );

    //
    // Open files
    //

    // Data
    TFile* mbFile = TFile::Open( mbFileName.Data() );
    TFile* mbPbGoingFile = TFile::Open( mbPbGoingFileName.Data() );
    TFile* mbPGoingFile = TFile::Open( mbPGoingFileName.Data() );
    TFile* jet60File = TFile::Open( jet60FileName.Data() );
    TFile* jet60PbGoingFile = TFile::Open( jet60PbGoingFileName.Data() );
    TFile* jet60PGoingFile = TFile::Open( jet60PGoingFileName.Data() );
    TFile* jet80File = TFile::Open( jet80FileName.Data() );
    TFile* jet80PbGoingFile = TFile::Open( jet80PbGoingFileName.Data() );
    TFile* jet80PGoingFile = TFile::Open( jet80PGoingFileName.Data() );
    TFile* jet100File = TFile::Open( jet100FileName.Data() );
    TFile* jet100PbGoingFile = TFile::Open( jet100PbGoingFileName.Data() );
    TFile* jet100PGoingFile = TFile::Open( jet100PGoingFileName.Data() );

    // JEU up/down
    TFile* mbJeuUpFile = TFile::Open( mbJeuUpFileName.Data() );
    TFile* mbJeuDownFile = TFile::Open( mbJeuDownFileName.Data() );
    TFile* jet60JeuUpFile = TFile::Open( jet60JeuUpFileName.Data() );
    TFile* jet60JeuDownFile = TFile::Open( jet60JeuDownFileName.Data() );
    TFile* jet80JeuUpFile = TFile::Open( jet80JeuUpFileName.Data() );
    TFile* jet80JeuDownFile = TFile::Open( jet80JeuDownFileName.Data() );
    TFile* jet100JeuUpFile = TFile::Open( jet100JeuUpFileName.Data() );
    TFile* jet100JeuDownFile = TFile::Open( jet100JeuDownFileName.Data() );

    // Embedding
    TFile* embeddingFile = TFile::Open( embeddingFileName.Data() );
    TFile* jerUpFile = TFile::Open( jerUpFileName.Data() );
    TFile* jerDownFile = TFile::Open( jerDownFileName.Data() );
    TFile* pbGoingEmbeddingFile = TFile::Open( pbGoingEmbeddingFileName.Data() );
    TFile* pGoingEmbeddingFile = TFile::Open( pGoingEmbeddingFileName.Data() );

    // Pileup
    TFile* mbGplusFile = TFile::Open( mbGplusFileName.Data() );
    TFile* mbVtx1File = TFile::Open( mbVtx1FileName.Data() );
    TFile* jet60GplusFile = TFile::Open( jet60GplusFileName.Data() );
    TFile* jet60Vtx1File = TFile::Open( jet60Vtx1FileName.Data() );
    TFile* jet80GplusFile = TFile::Open( jet80GplusFileName.Data() );
    TFile* jet80Vtx1File = TFile::Open( jet80Vtx1FileName.Data() );
    TFile* jet100GplusFile = TFile::Open( jet100GplusFileName.Data() );
    TFile* jet100Vtx1File = TFile::Open( jet100Vtx1FileName.Data() );

    // NLOJet++
    TFile* epps21File = TFile::Open( epps21FileName.Data() );
    TFile* ncteq15hqFile = TFile::Open( ncteq15hqFileName.Data() );

    
    // Check that files exist
    

    // Data
    if ( !checkFileIsGood( mbFile ) ) return;
    if ( !checkFileIsGood( mbPbGoingFile ) ) return;
    if ( !checkFileIsGood( mbPGoingFile ) ) return;
    if ( !checkFileIsGood( jet60File ) ) return;
    if ( !checkFileIsGood( jet60PbGoingFile ) ) return;
    if ( !checkFileIsGood( jet60PGoingFile ) ) return;
    if ( !checkFileIsGood( jet80File ) ) return;
    if ( !checkFileIsGood( jet80PbGoingFile ) ) return;
    if ( !checkFileIsGood( jet80PGoingFile ) ) return;
    if ( !checkFileIsGood( jet100File ) ) return;
    if ( !checkFileIsGood( jet100PbGoingFile ) ) return;
    if ( !checkFileIsGood( jet100PGoingFile ) ) return;

    // JEU
    if ( !checkFileIsGood( mbJeuUpFile ) ) return;
    if ( !checkFileIsGood( mbJeuDownFile ) ) return;
    if ( !checkFileIsGood( jet60JeuUpFile ) ) return;
    if ( !checkFileIsGood( jet60JeuDownFile ) ) return;
    if ( !checkFileIsGood( jet80JeuUpFile ) ) return;
    if ( !checkFileIsGood( jet80JeuDownFile ) ) return;
    if ( !checkFileIsGood( jet100JeuUpFile ) ) return;
    if ( !checkFileIsGood( jet100JeuDownFile ) ) return;

    // Embedding
    if ( !checkFileIsGood( embeddingFile ) ) return;
    if ( !checkFileIsGood( jerUpFile ) ) return;
    if ( !checkFileIsGood( jerDownFile ) ) return;
    if ( !checkFileIsGood( pbGoingEmbeddingFile ) ) return;
    if ( !checkFileIsGood( pGoingEmbeddingFile ) ) return;

    // Pileup
    if ( !checkFileIsGood( mbGplusFile ) ) return;
    if ( !checkFileIsGood( mbVtx1File ) ) return;
    if ( !checkFileIsGood( jet60GplusFile ) ) return;
    if ( !checkFileIsGood( jet60Vtx1File ) ) return;
    if ( !checkFileIsGood( jet80GplusFile ) ) return;
    if ( !checkFileIsGood( jet80Vtx1File ) ) return;
    if ( !checkFileIsGood( jet100GplusFile ) ) return;
    if ( !checkFileIsGood( jet100Vtx1File ) ) return;

    // NLOJet++
    if ( !checkFileIsGood( epps21File ) ) return;
    if ( !checkFileIsGood( ncteq15hqFile ) ) return;

    retrieveDistributions(mbFile, mbPbGoingFile, mbPGoingFile, 
                          jet60File, jet60PbGoingFile, jet60PGoingFile,
                          jet80File, jet80PbGoingFile, jet80PGoingFile,
                          jet100File, jet100PbGoingFile, jet100PGoingFile,
                          mbJeuUpFile, mbJeuDownFile,
                          jet60JeuUpFile, jet60JeuDownFile,
                          jet80JeuUpFile, jet80JeuDownFile,
                          jet100JeuUpFile, jet100JeuDownFile,
                          embeddingFile, jerUpFile, jerDownFile,
                          pbGoingEmbeddingFile, pGoingEmbeddingFile,
                          mbGplusFile, mbVtx1File,
                          jet60GplusFile, jet60Vtx1File,
                          jet80GplusFile, jet80Vtx1File,
                          jet100GplusFile, jet100Vtx1File,
                          epps21File, ncteq15hqFile,
                          date, isCM, drawFits);
}

//________________
void systematics() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // CMS or Lab frame
    // Bool_t isCM{kTRUE};
    Bool_t isCM{kFALSE};

    // Bool_t drawFits{kTRUE};
    Bool_t drawFits{kFALSE};

    // Run old calculations
    // runOldCalculations(isCM);

    // Run new calculations
    runNewCalculations(isCM, drawFits);
}