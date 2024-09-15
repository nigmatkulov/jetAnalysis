// ROOT headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystemDirectory.h"
#include <sys/stat.h>
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

// C++ headers
#include <fstream>
#include <iostream>
#include <vector>

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
}

//________________
void setMultiGraphStyle(TMultiGraph *mg) {

    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.06);
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetNdivisions(208);
    mg->GetYaxis()->SetNdivisions(205);    
    mg->GetYaxis()->SetTitleOffset(1.1);

    mg->GetXaxis()->SetTitle("#eta^{dijet}");
    mg->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    // std::cout << "Set style to histo: " << h->GetName() << " type: " << type << std::endl;
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.1;
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
        markerStyle = 22; // filled triangle
    }
    else if (type == 3) {
        color = 2;        // red
        markerStyle = 24; // open circle
    }
    else if (type == 4) {
        color = 4;        // blue
        markerStyle = 20; // filled circle
    }
    else if (type == 5 ) {
        color = 8;        // green
        markerStyle = 43; // diamond
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
void set1DStyle(TGraphErrors *h) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.1;
    Int_t lineWidth = 2;
    Int_t color = 1;
    markerStyle = 20; // filled circle

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
}

//________________
void plotCMSHeader() {
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.57, 0.93, "pPb 174.56 nb^{-1} (8.16 TeV)");
    t.SetTextSize(0.05);
}

//________________
void relSyst(TGraph *gr, TH1D* h, Int_t style = 0) {
    // Create histogram
    set1DStyle(h, style);
    TString name = gr->GetName();
    std::cout << "Source: " << gr->GetName() << std::endl;
    for (Int_t i{1}; i <= h->GetNbinsX(); i++) {
        h->SetBinContent(i, 0.);
    }
    
    // Loop over points
    for (Int_t i{0}; i<gr->GetN(); i++) {

        Double_t xVal{0.}, yVal{0.};
        gr->GetPoint(i, xVal, yVal);

        //std::cout << "i: " << i << " x: " << xVal << " y: " << yVal << std::endl;
        if ( TMath::Abs( yVal ) < 0.0001 ) {
            h->SetBinContent(i+1, 0);
        }
        else {
            h->SetBinContent(i+1, TMath::Abs(yVal) );
        }
        //std::cout << "i: " << i << " x: " << h->GetBinCenter(i+1) << " bin width: " << h->GetXaxis()->GetBinWidth(i+1) <<" uncrt [%]: " << h->GetBinContent(i+1) * 100 << std::endl;
    } // for (Int_t i{0}; i<gr->GetN(); i++)

    //
    // Temporary smoothing procedure
    //
    // h->Smooth();
    if ( name.Contains("JeuSyst") || name.Contains("PointingSyst") ) {
        Int_t binLow = h->FindBin(-1.2);
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
                if ( gr->GetPointY(i-1) < 0 ) {
                    gr->SetPointY(i-1, -1. * maxVal);
                }
            }
        } // for (Int_t i = binLow; i<=binHi; i++)
    } // if ( name.Contains("jeuSyst") || name.Contains("pointingSyst") )
}

//________________
void addSystSource(TH1D* hTot, TH1D* hJeu, TH1D* hJer, TH1D* hPointing, TH1D* hPileup) {
    // std::cout << "Add systematic source" << std::endl;
    Double_t yTot{0}, yJeu{0}, yJer{0}, yPointing{0}, yPileup{0};
    for (Int_t i{1}; i <= hTot->GetNbinsX(); i++) {
        yJeu = hJeu->GetBinContent(i);
        yJer = hJer->GetBinContent(i);
        yPointing = hPointing->GetBinContent(i);
        yPileup = hPileup->GetBinContent(i);
        yTot = TMath::Sqrt(yJeu*yJeu + yJer*yJer + yPointing*yPointing + yPileup*yPileup);
        std::cout << "i: " << i << " eta: " << hTot->GetBinCenter(i) << " yJES: " << yJeu << " yJER: " << yJer << " yPoint: " << yPointing
                  << " yPileup: " << yPileup
                  << " yTot: " << yTot << std::endl;
        hTot->SetBinContent(i, yTot);
    }
}

//________________
TGraphErrors *totalSyst(TGraph* jeuSyst = nullptr, TGraph *jerSyst = nullptr,
                        TGraph* pointingSyst = nullptr, TGraph* pileupSyst = nullptr,
                        TGraphErrors *data = nullptr, 
                        TGraph* dirSyst = nullptr) {
    Double_t dirVal{0}, jeuVal{0}, jerVal{0}, pointingVal{0}, pileupVal{0}, totalUncrt{0};
    Double_t xVal{0}, yVal{0}, tmp{0};
    TGraphErrors *gr = new TGraphErrors();

    std::cout << "Creating final systematic uncertainty ";
    // Loop over data points
    for (Int_t i=0; i<data->GetN(); i++) {
        data->GetPoint(i, xVal, yVal);

        if (pointingSyst) {
            pointingSyst->GetPoint(i, tmp, pointingVal);
        }
        if (dirSyst) {
            dirSyst->GetPoint(i, tmp, dirVal);
        }
        if (jeuSyst) {
            jeuSyst->GetPoint(i, tmp, jeuVal);
        }
        if (jerSyst) {
            jerSyst->GetPoint(i, tmp, jerVal);
        }
        if (pileupSyst) {
            pileupSyst->GetPoint(i, tmp, pileupVal);
        }

        totalUncrt = TMath::Sqrt( dirVal*dirVal + jeuVal*jeuVal + jerVal*jerVal + pointingVal*pointingVal + pileupVal*pileupVal);

        // std::cout << "Total uncert [%]: " << totalUncrt * 100 << " dirVal: " << dirVal * 100 << " jeuVal: " 
        //           << jeuVal * 100 << " jerVal: " << jerVal * 100  << " pointingVal: " << pointingVal * 100
        //           << " pileupVal: " << pileupVal*100
        //           << std::endl;
        // std::cout << "Data y: " << yVal << " ey: " << data->GetErrorY(i) << std::endl;


        gr->SetPoint(i, xVal, yVal);
        gr->SetPointError(i, data->GetErrorX(i), totalUncrt * yVal);
    }

    gr->GetXaxis()->SetTitle("#eta^{dijet}");
    gr->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");

    std::cout << "\t[DONE]" << std::endl;

    return gr;
}

//________________
void writeCovPoints(TGraphErrors* grData, TGraph* grSyst, TH1* hSyst, Int_t systSource, Int_t iPt) {

    //
    // systSource: 
    // 0 - Data
    // 1 - JES
    // 2 - JER
    // 3 - Pointing
    // 4 - Pileup
    // 5 - Total syst
    //

    //
    // iPt - corresponds to the pTave bin in the array 
    //

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetBinLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    std::vector<Int_t> ptDijetBinHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194 };

    std::vector<Int_t> ptDijetPtLow{};
    std::vector<Int_t> ptDijetPtHi{};

    for (UInt_t i{0}; i<ptDijetBinLow.size(); i++) {
        ptDijetPtLow.push_back( ptLow + (ptDijetBinLow.at(i)-1) * ptStep );
        ptDijetPtHi.push_back( ptLow + ptDijetBinHi.at(i) * ptStep );
    }

    // Define variables
    Double_t xDataVal{0}, yDataVal{0}, xDataValErr{0}, yDataValErr{0}, ySystSign{1};


    TGraphErrors *gr2write = new TGraphErrors;
    Int_t iPoint{0};
    Int_t nPoints2Skip{2};

    TString systText;
    if ( systSource == 0 ) {
        systText = "data";
    }
    else if ( systSource == 1 ) {
        systText = "jesSyst";
    }
    else if ( systSource == 2 ) {
        systText = "jerSyst";
    }
    else if ( systSource == 3 ) {
        systText = "pointingSyst";
    }
    else if ( systSource == 4 ) {
        systText = "pileupSyst";
    }
    else if ( systSource == 5 ) {
        systText = "totalSyst";
    }
    else {
        systText = "unknownSyst";
    }

    // Create file to write points and errors to ASCII file
    std::ofstream outFile( Form("freezeSyst/covMtx/pPb8160_dijet_%s_pt_%d_%d.txt", systText.Data(), 
                                ptDijetPtLow.at( iPt ), ptDijetPtHi.at( iPt ) ) );

    if ( outFile.is_open() ) { // Check if file is open

        // Loop over points in TGraph
        for (Int_t i=nPoints2Skip; i<(grData->GetN()-nPoints2Skip); i++) {

            // Initialize values
            xDataVal = {0}; yDataVal = {0}; xDataValErr = {0}; yDataValErr = {0}; ySystSign = {1};

            // Retrieve point for the data
            grData->GetPoint(i, xDataVal, yDataVal);

            //
            // Retrieve uncertainties
            //
            xDataValErr = grData->GetErrorX(i);
            if ( systSource == 0 ) { // Data
                yDataValErr = grData->GetErrorY(i);
            }
            else { // Other sources
                xDataValErr /= 2;
                yDataValErr = hSyst->GetBinContent( i+1 );
                yDataValErr *= 0.01;
                yDataValErr *= yDataVal;
                // std::cout << "Data x: " << xDataVal << " data y: " << yDataVal << " systVal (%): " << hSyst->GetBinContent( i+1 ) << " systUncr: " << yDataValErr << std::endl;
                if ( systSource < 5 ) { // For all sources, except total
                    Double_t tmp = grSyst->GetPointY(i);
                    //std::cout << "tmp: " << tmp << " yDataValErr: " << yDataValErr << std::endl; 
                    yDataValErr = ( tmp >=0 ) ? yDataValErr : (-1.) * yDataValErr;
                }  
            }

            gr2write->SetPoint(iPoint, xDataVal, yDataVal);
            gr2write->SetPointError(iPoint, xDataValErr, yDataValErr);

            outFile << xDataVal << " " << yDataVal << " " << xDataValErr << " " << yDataValErr << std::endl;

            iPoint++;
        } // for (Int_t i=0; i<grSyst->GetN(); i++)

        outFile.close();
    } // if ( outFile.is_open() )
    else {
        std::cerr << "[writeCovPoints] Unable to open file for writing" << std::endl;
    }

}

//________________
void createValues() {

    // Bool_t useDirSyst{kTRUE};
    Bool_t useDirSyst{kFALSE};

    Bool_t useJeuSyst{kTRUE};
    // Bool_t useJeuSyst{kFALSE};

    Bool_t useJerSyst{kTRUE};
    // Bool_t useJerSyst{kFALSE};

    Bool_t usePointingSyst{kTRUE};
    // Bool_t usePointingSyst{kFALSE};

    Bool_t usePileupSyst{kTRUE};
    // Bool_t usePileupSyst{kFALSE};

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetBinLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    std::vector<Int_t> ptDijetBinHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194 };

    std::vector<Int_t> ptDijetPtLow{};
    std::vector<Int_t> ptDijetPtHi{};

    for (UInt_t i{0}; i<ptDijetBinLow.size(); i++) {
        ptDijetPtLow.push_back( ptLow + (ptDijetBinLow.at(i)-1) * ptStep );
        ptDijetPtHi.push_back( ptLow + ptDijetBinHi.at(i) * ptStep );
    }

    const Int_t dijetEtaBins{30};
    Double_t dijetEtaVals[dijetEtaBins+1] { -5.0, -4.0, -3.0, -2.4, -2.2, 
                                            -2.0, -1.8, -1.6, -1.4, -1.2, 
                                            -1.0, -0.8, -0.6, -0.4, -0.2,  
                                             0.0,  0.2,  0.4,  0.6,  0.8,  
                                             1.0,  1.2,  1.4,  1.6,  1.8,  
                                             2.0,  2.2,  2.4,  3.0,  4.0,  
                                             5.0 };

    Int_t totalType{2};
    Int_t dirType{0};
    Int_t jeuType{1};
    Int_t pointingType{5};
    Int_t jerType{6};
    Int_t pileupType{0};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TLegend *leg;

    // Create histograms
    TH1D *hDirSyst[ ptDijetPtLow.size() ];
    TH1D *hJeuSyst[ ptDijetPtLow.size() ];
    TH1D *hJerSyst[ ptDijetPtLow.size() ];
    TH1D *hPointingSyst[ ptDijetPtLow.size() ];
    TH1D *hPileupSyst[ ptDijetPtLow.size() ];
    TH1D *hTotalSyst[ ptDijetPtLow.size() ];
    TH1D *hTotalSystRel[ ptDijetPtLow.size() ];

    // Create graphs
    TGraphErrors *grDijetEta[ ptDijetPtLow.size() ];
    TGraph *grDijetDirRelSyst[ ptDijetPtLow.size() ];
    TGraph *grDijetJeuRelSyst[ ptDijetPtLow.size() ];
    TGraph *grDijetJerRelSyst[ ptDijetPtLow.size() ];
    TGraph *grDijetPointingRelSyst[ ptDijetPtLow.size() ];
    TGraph *grDijetPileupRelSyst[ ptDijetPtLow.size() ];
    for (int i{0}; i<ptDijetPtLow.size(); i++ ) {
        grDijetEta[i] = nullptr;
        grDijetDirRelSyst[i] = nullptr;
        grDijetJeuRelSyst[i] = nullptr;
        grDijetJerRelSyst[i] = nullptr;
        grDijetPointingRelSyst[i] = nullptr;
        grDijetPileupRelSyst[i] = nullptr;
    }
    TGraphErrors *grErrTotalSystUncrt[ ptDijetPtLow.size() ];
    TMultiGraph *mgDijetEta[ ptDijetPtLow.size() ];

    // Create canvases
    Int_t sizeX = 1000;
    Int_t sizeY = 1200;

    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetPtLow.size();

    TCanvas *cSyst = new TCanvas("cSyst", "cSyst", sizeX, sizeY);
    cSyst->Divide(4, ( (ptBins % 4) == 0 ) ? (ptBins / 4) : (ptBins / 4 + 1) );

    TString triggerName = "MB";

    // Loop over pT bins
    for (UInt_t i{0}; i<ptDijetPtLow.size(); i++) {

        if (ptDijetPtLow.at(i) < 80) {
            triggerName = "MB";
        }
        else if (ptDijetPtLow.at(i) < 100) {
            triggerName = "Jet60";
        }
        else if (ptDijetPtLow.at(i) < 120) {
            triggerName = "Jet80";
        }
        else {
            triggerName = "Jet100";
        }

        // std::cout << "ptDijetPtLow.at(i): " << ptDijetPtLow.at(i) << " trgName: " << triggerName.Data() << std::endl;
        // std::cout << "Going to open: " <<  Form("freezeSyst/%s_pPb8160_etaDijet_data_%d_%d_lab.txt", triggerName.Data(), ptDijetPtLow.at(i), ptDijetPtHi.at(i)) << std::endl;

        // Retrieve data for the give pT average bin
        grDijetEta[i] = new TGraphErrors( Form("freezeSyst/systematics/%s_pPb8160_etaDijet_data_%d_%d_lab.txt", triggerName.Data(), ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");

        grDijetEta[i]->SetName( Form("grDijetEta_%d", i) );
        set1DStyle(grDijetEta[i]);
        // for (Int_t j{0}; j<grDijetEta[i]->GetN(); j++) {
        //     std::cout << "Data j: " << j << " x: " << grDijetEta[i]->GetPointX(j) << " y: " << grDijetEta[i]->GetPointY(j) << " ey: " << grDijetEta[i]->GetErrorY(j) << std::endl;
        // }

        writeCovPoints(grDijetEta[i], nullptr, nullptr, 0, i);

        // Create histogram for total realtive systematic uncertainty
        hTotalSyst[i] = new TH1D( Form("hTotalSyst_%d",i), "Total relative systematic uncertainty;#eta^{dijet};Relative uncertainty [%]", 50, -5., 5. );
        hTotalSyst[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hTotalSyst[i]->Sumw2();
        set1DStyle(hTotalSyst[i], 2);
        //std::cout << "i: " << i << " Nbins: " << hTotalSyst[i]->GetNbinsX() << std::endl;

        // Retrieve p-going vs Pb-going systematics
        if ( useDirSyst ) {
            grDijetDirRelSyst[i] = new TGraph(Form("freezeSyst/systematics/pPb8160_etaDijet_dirSyst_%d_%d_lab.txt", 40, 1000), "%lg %lg", " ");
            grDijetDirRelSyst[i]->SetName( Form("grDirSyst_%d", i) );
            hDirSyst[i] = new TH1D(  Form("hDirSyst_%d",i), "Relative systematic uncertainty [%]", 30, -5., 5.);
            hDirSyst[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hDirSyst[i]->Sumw2();
            relSyst(grDijetDirRelSyst[i], hDirSyst[i], dirType);
            hDirSyst[i]->Scale(100.);
            hDirSyst[i]->SetName( Form("hDirSyst_%d",i) );
        }

        // Retrieve JEU systematics
        if ( useJeuSyst ) {

            std::cout << "JES systematics" << std::endl;
            grDijetJeuRelSyst[i] = new TGraph(Form("freezeSyst/systematics/%s_pPb8160_etaDijet_jeuSyst_%d_%d_lab.txt", triggerName.Data(), ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg", " ");
            grDijetJeuRelSyst[i]->SetName( Form("grJeuSyst_%d", i) );
            hJeuSyst[i] = new TH1D(  Form("hJeuSyst_%d",i), "Relative systematic uncertainty [%]", 30, -5., 5.);
            hJeuSyst[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hJeuSyst[i]->Sumw2();
            relSyst(grDijetJeuRelSyst[i], hJeuSyst[i], jeuType);
            hJeuSyst[i]->Scale(100.);
            hJeuSyst[i]->SetName( Form("hJeuSyst_%d",i) );

            writeCovPoints(grDijetEta[i], grDijetJeuRelSyst[i], hJeuSyst[i], 1, i);
        }

        // Retrieve JER systematics
        if ( useJerSyst ) {
            std::cout << "JER systematics" << std::endl;
            grDijetJerRelSyst[i] = new TGraph(Form("freezeSyst/systematics/MB_pPb8160_etaDijet_jerSyst_%d_%d_lab.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg", " ");
            grDijetJerRelSyst[i]->SetName( Form("grJerSyst_%d", i) );
            hJerSyst[i] = new TH1D(  Form("hJerSyst_%d",i), "Relative systematic uncertainty [%]", 30, -5., 5.);
            hJerSyst[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hJerSyst[i]->Sumw2();
            relSyst(grDijetJerRelSyst[i], hJerSyst[i], jerType);
            hJerSyst[i]->Scale(100.);
            hJerSyst[i]->SetName( Form("hJerSyst_%d",i) );

            writeCovPoints(grDijetEta[i], grDijetJerRelSyst[i], hJerSyst[i], 2, i);
        }

        // Retrieve pointing resolution systematics
        if ( usePointingSyst ) {
            std::cout << "Pointing resolution systematics" << std::endl;
            grDijetPointingRelSyst[i] = new TGraph(Form("freezeSyst/systematics/MB_pPb8160_etaDijet_pointingSyst_%d_%d_lab.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg", " ");
            grDijetPointingRelSyst[i]->SetName( Form("grPointingSyst_%d", i) );
            hPointingSyst[i] = new TH1D(  Form("hPointingSyst_%d",i), "Relative systematic uncertainty [%]", 30, -5., 5.);
            hPointingSyst[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hPointingSyst[i]->Sumw2();
            relSyst(grDijetPointingRelSyst[i], hPointingSyst[i], pointingType);
            hPointingSyst[i]->Scale(100.);
            hPointingSyst[i]->SetName( Form("hPointingSyst_%d",i) );

            writeCovPoints(grDijetEta[i], grDijetPointingRelSyst[i], hPointingSyst[i], 3, i);
        }

        // Retrieve systematics due to pileup effect
        if ( usePileupSyst ) {
            std::cout << "Pileup systematics" << std::endl;
            grDijetPileupRelSyst[i] = new TGraph(Form("freezeSyst/systematics/%s_pPb8160_etaDijet_pileupSyst_%d_%d_lab.txt", triggerName.Data(), ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg", " ");
            grDijetPileupRelSyst[i]->SetName( Form("grPileupSyst_%d", i) );
            hPileupSyst[i] = new TH1D(  Form("hPileupSyst_%d",i), "Relative systematic uncertainty [%]", 30, -5., 5.);
            hPileupSyst[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hPileupSyst[i]->Sumw2();
            relSyst(grDijetPileupRelSyst[i], hPileupSyst[i], pileupType);
            hPileupSyst[i]->Scale(100.);
            hPileupSyst[i]->SetName( Form("hPileupSyst_%d",i) );
            writeCovPoints(grDijetEta[i], grDijetPileupRelSyst[i], hPileupSyst[i], 4, i);
        }

        //
        // Calculate total systematic uncertainty
        //
        addSystSource(hTotalSyst[i], hJeuSyst[i], hJerSyst[i], hPointingSyst[i], hPileupSyst[i]);
        writeCovPoints(grDijetEta[i], nullptr, hTotalSyst[i], 5, i);
        

        //
        // Create graph with total systematic uncertainties
        //
        //std::cout << "Address of dir syst: " << grDijetDirRelSyst[i] << " name: " << grDijetDirRelSyst[i]->GetName() << std::endl;
        grErrTotalSystUncrt[i] = totalSyst( grDijetJeuRelSyst[i], grDijetJerRelSyst[i], grDijetPointingRelSyst[i], grDijetPileupRelSyst[i], grDijetEta[i] /*, grDijetDirRelSyst[i] */);
        grErrTotalSystUncrt[i]->SetName( Form("grErrTotalSystUncrt_%d", i));
        grErrTotalSystUncrt[i]->SetFillColor(29);
        grErrTotalSystUncrt[i]->SetFillStyle(1001); // Solid fill

        mgDijetEta[i] = new TMultiGraph();
        mgDijetEta[i]->SetName( Form("mgDijetEta_%d",i) );
        mgDijetEta[i]->Add(grErrTotalSystUncrt[i], "A5");
        mgDijetEta[i]->Add(grDijetEta[i], "AP");

        //
        // Plot individual uncertainties
        //
        canv->cd();
        setPadStyle();
        hTotalSyst[i]->Draw();
        if ( useDirSyst ) {
            hDirSyst[i]->Draw("same");
            // hDirSyst[i]->Draw();
        }
        if ( useJeuSyst ) {
            hJeuSyst[i]->Draw("same");
            // hJeuSyst[i]->Draw();
        }
        if ( useJerSyst ) {
            hJerSyst[i]->Draw("same");
        }
        if ( usePointingSyst ) {
            hPointingSyst[i]->Draw("same");
        }
        if ( usePileupSyst ) {
            hPileupSyst[i]->Draw("same");
        }
        hTotalSyst[i]->GetYaxis()->SetRangeUser(0., 20.);
        hTotalSyst[i]->GetXaxis()->SetRangeUser(-3., 3.);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        leg = new TLegend(0.4, 0.65, 0.65, 0.85);
        leg->SetTextSize(0.05);
        leg->SetLineWidth(0);
        leg->AddEntry(hTotalSyst[i], "Total", "l");
        if (useDirSyst) {
            leg->AddEntry(hDirSyst[i], Form("Beam dir."), "l");  
        }
        if (useJeuSyst) {
            leg->AddEntry(hJeuSyst[i], Form("JES"), "l");  
        }
        if (useJerSyst) {
            leg->AddEntry(hJerSyst[i], Form("JER"), "l");  
        }
        if ( usePointingSyst ) {
            leg->AddEntry(hPointingSyst[i], Form("Pointing res."), "l");  
        }
        if ( usePileupSyst ) {
            leg->AddEntry(hPileupSyst[i], Form("Pileup"), "l");
        }
        leg->Draw();
        canv->SaveAs( Form("freezeSyst/systematics/pPb8160_syst_pt_%d_%d_lab.pdf", ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        canv->SaveAs( Form("freezeSyst/systematics/pPb8160_syst_pt_%d_%d_lab.png", ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        
        //
        // Plot systematic uncertainties for all pT bins
        //
        cSyst->cd(i+1);
        setPadStyle();
        hTotalSyst[i]->Draw();
        if ( useDirSyst ) {
            hDirSyst[i]->Draw("same");
            // hDirSyst[i]->Draw();
        }
        if ( useJeuSyst ) {
            hJeuSyst[i]->Draw("same");
            // hJeuSyst[i]->Draw();
        }
        if ( useJerSyst ) {
            hJerSyst[i]->Draw("same");
        }
        if ( usePointingSyst ) {
            hPointingSyst[i]->Draw("same");
        }
        if ( usePileupSyst ) {
            hPileupSyst[i]->Draw("same");
        }
        hTotalSyst[i]->GetYaxis()->SetRangeUser(0., 20.);
        hTotalSyst[i]->GetXaxis()->SetRangeUser(-3., 3.);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        leg = new TLegend(0.4, 0.65, 0.65, 0.85);
        leg->SetTextSize(0.05);
        leg->SetLineWidth(0);
        leg->AddEntry(hTotalSyst[i], "Total", "l");
        if (useDirSyst) {
            leg->AddEntry(hDirSyst[i], Form("Beam dir."), "l");  
        }
        if (useJeuSyst) {
            leg->AddEntry(hJeuSyst[i], Form("JES"), "l");  
        }
        if (useJerSyst) {
            leg->AddEntry(hJerSyst[i], Form("JER"), "l");  
        }
        if ( usePointingSyst ) {
            leg->AddEntry(hPointingSyst[i], Form("Pointing res."), "l");  
        }
        if ( usePileupSyst ) {
            leg->AddEntry(hPileupSyst[i], Form("Pileup"), "l");  
        }

        leg->Draw();

    } // for (UInt_t i{0}; i<ptDijetBinLow.size(); i++)
    cSyst->SaveAs( Form("freezeSyst/systematics/pPb8160_syst_pt_all_lab.pdf") );
    cSyst->SaveAs( Form("freezeSyst/systematics/pPb8160_syst_pt_all_lab.png") );

    TFile *oFile = new TFile("freezeSyst/oSystematics.root", "recreate");
    for (UInt_t i{0}; i<ptDijetBinLow.size(); i++) {
        hTotalSystRel[i] = dynamic_cast<TH1D*> ( hTotalSyst[i]->Clone( Form("hTotalSystRel_%d", i) ) );
        //hTotalSystRel[i]->Sumw2();
        hTotalSystRel[i]->Scale(1./100);
        hTotalSystRel[i]->SetMarkerStyle(1);
        hTotalSystRel[i]->SetFillColor(41);
        hTotalSystRel[i]->SetMarkerColor(41);
        hTotalSystRel[i]->SetLineColor(41);
        for (Int_t j=1; j<hTotalSystRel[i]->GetNbinsX(); j++) {
            Double_t uncrt = hTotalSystRel[i]->GetBinContent(j);
            hTotalSystRel[i]->SetBinContent(j, 1.);
            hTotalSystRel[i]->SetBinError(j, uncrt);
        }
        hTotalSystRel[i]->Write();
        hTotalSyst[i]->Write();
        hJeuSyst[i]->Write();
        hJerSyst[i]->Write();
        hPointingSyst[i]->Write();
        hPileupSyst[i]->Write();
        grDijetEta[i]->Write();
        grErrTotalSystUncrt[i]->Write();
    } // for (UInt_t i{0}; i<ptDijetBinLow.size(); i++)
    oFile->Close();
}

//________________
void plotFinalDistributions() {

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetBinLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    std::vector<Int_t> ptDijetBinHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194 };

    std::vector<Int_t> ptDijetPtLow{};
    std::vector<Int_t> ptDijetPtHi{};

    for (UInt_t i{0}; i<ptDijetBinLow.size(); i++) {
        ptDijetPtLow.push_back( ptLow + (ptDijetBinLow.at(i)-1) * ptStep );
        ptDijetPtHi.push_back( ptLow + ptDijetBinHi.at(i) * ptStep );
    }

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    TLegend *leg;

    // Create graphs
    TGraphErrors *grDijetEta[ ptDijetPtLow.size() ];
    TGraphErrors *grErrTotalSystUncrt[ ptDijetPtLow.size() ];
    TMultiGraph *mgDijetEta[ ptDijetPtLow.size() ];

    // Create canvases
    Int_t sizeX = 1000;
    Int_t sizeY = 1200;
    Int_t ptBins = ptDijetPtLow.size();
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    TCanvas *cTot = new TCanvas("cTot", "cTot", sizeX, sizeY);
    cTot->Divide(4, ( (ptBins % 4) == 0 ) ? (ptBins / 4) : (ptBins / 4 + 1) );

    // Loop over pT bins
    for (UInt_t i{0}; i<ptDijetPtLow.size(); i++) {

        // Retrieve data for the given pT average bin
        grDijetEta[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_data_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grDijetEta[i]->SetName( Form("grDijetEta_%d", i) );
        set1DStyle(grDijetEta[i]);

        // Retrieve total systematic uncertainty for the given pT average bin
        grErrTotalSystUncrt[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_totalSyst_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grErrTotalSystUncrt[i]->SetName( Form("grDijetEta_%d", i) );
        set1DStyle(grErrTotalSystUncrt[i]);
        grErrTotalSystUncrt[i]->SetFillColor(29);
        grErrTotalSystUncrt[i]->SetFillStyle(1001); // Solid fill

        // Create multigraph with the data and systematic uncertainty
        mgDijetEta[i] = new TMultiGraph();
        mgDijetEta[i]->SetName( Form("mgDijetEta_%d",i) );
        mgDijetEta[i]->Add(grErrTotalSystUncrt[i], "A5");
        mgDijetEta[i]->Add(grDijetEta[i], "AP");

        //
        // Plot individual pTave final dijet eta distribution
        //
        canv->cd();
        setPadStyle();
        setMultiGraphStyle(mgDijetEta[i]);
        mgDijetEta[i]->Draw("AP");
        mgDijetEta[i]->GetYaxis()->SetRangeUser(0.001, 0.15);
        //mgDijetEta[i]->GetXaxis()->SetRangeUser(-3., 3.);
        t.DrawLatexNDC(0.2, 0.83, Form("%d < p_{T}^{ave} GeV < %d", 
                       ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
        t.DrawLatexNDC(0.2, 0.67, "p_{T}^{Subleading} > 40 GeV");
        t.DrawLatexNDC(0.2, 0.59, "|#eta| < 3");
        t.DrawLatexNDC(0.2, 0.51, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
        plotCMSHeader();

        leg = new TLegend(0.7, 0.7, 0.85, 0.85);
        leg->SetTextSize(0.05);
        leg->SetLineWidth(0);
        leg->AddEntry(grDijetEta[i], "Data", "pl");
        leg->AddEntry(grErrTotalSystUncrt[i], Form("Syst. ucrt."), "f");
        leg->Draw();
        canv->SaveAs( Form("freezeSyst/finalPlots/pPb8160_dijet_pt_%d_%d_lab.pdf", ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        canv->SaveAs( Form("freezeSyst/finalPlots/pPb8160_dijet_pt_%d_%d_lab.png", ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );

        //
        // Plot all final plots
        //
        cTot->cd(i+1);
        setPadStyle();
        mgDijetEta[i]->Draw("AP");
        mgDijetEta[i]->GetYaxis()->SetRangeUser(0.005, 0.15);
        mgDijetEta[i]->GetXaxis()->SetRangeUser(-3., 3.);
        // gPad->SetLogy(1);
        t.DrawLatexNDC(0.2, 0.83, Form("%d < p_{T}^{ave} GeV < %d", 
                       ptDijetPtLow.at(i), ptDijetPtHi.at(i)) );
        t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
        t.DrawLatexNDC(0.2, 0.67, "p_{T}^{Subleading} > 40 GeV");
        t.DrawLatexNDC(0.2, 0.59, "|#eta| < 3");
        t.DrawLatexNDC(0.2, 0.51, "#Delta#phi^{dijet} > #frac{5#pi}{6}");

        plotCMSHeader();

        leg = new TLegend(0.7, 0.7, 0.85, 0.85);
        leg->SetTextSize(0.05);
        leg->SetLineWidth(0);
        leg->AddEntry(grDijetEta[i], "Data", "pl");
        leg->AddEntry(grErrTotalSystUncrt[i], Form("Syst. ucrt."), "f");
        leg->Draw();
    } // for (UInt_t i{0}; i<ptDijetPtLow.size(); i++)

    cTot->SaveAs( Form("freezeSyst/finalPlots/pPb8160_dijetEta_pt_all_lab.pdf") );
    cTot->SaveAs( Form("freezeSyst/finalPlots/pPb8160_dijetEta_pt_all_lab.png") );
}

//________________
void calculateCovMatrix() {

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetBinLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55, 75, 95,  55  };
    std::vector<Int_t> ptDijetBinHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 74, 94, 194, 194 };

    std::vector<Int_t> ptDijetPtLow{};
    std::vector<Int_t> ptDijetPtHi{};

    for (UInt_t i{0}; i<ptDijetBinLow.size(); i++) {
        ptDijetPtLow.push_back( ptLow + (ptDijetBinLow.at(i)-1) * ptStep );
        ptDijetPtHi.push_back( ptLow + ptDijetBinHi.at(i) * ptStep );
    }

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    const Int_t nPtBins = ptDijetPtLow.size();

    // Create graphs
    TGraphErrors *grData[ ptDijetPtLow.size() ];
    TGraphErrors *grJES[ ptDijetPtLow.size() ];
    TGraphErrors *grJER[ ptDijetPtLow.size() ];
    TGraphErrors *grPileup[ ptDijetPtLow.size() ];
    TGraphErrors *grPointing[ ptDijetPtLow.size() ];

    // Loop over pT bins to read 
    for (UInt_t i{0}; i<ptDijetPtLow.size(); i++) {
        grData[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_data_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grData[i]->SetName( Form("grData_%d", i) );
        grJES[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_jesSyst_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grJES[i]->SetName( Form("grJES_%d", i) );
        grJER[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_jerSyst_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grJER[i]->SetName( Form("grJER_%d", i) );
        grPileup[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_pileupSyst_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grPileup[i]->SetName( Form("grPileup_%d", i) );
        grPointing[i] = new TGraphErrors( Form("freezeSyst/covMtx/pPb8160_dijet_pointingSyst_pt_%d_%d.txt", ptDijetPtLow.at(i), ptDijetPtHi.at(i)), "%lg %lg %lg %lg", "");
        grPointing[i]->SetName( Form("grPointing_%d", i) );
    } // for (UInt_t i{0}; i<ptDijetPtLow.size(); i++)

    const Int_t nPoints = grJES[0]->GetN();

    TH2D *hCovMtxData = new TH2D("hCovMtxData", "Covariance matrix for data;N_{p_{T}^{ave}} x N_{#eta^{dijet}};N_{p_{T}^{ave}} x N_{#eta^{dijet}}", 
                                nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1,
                                nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1);
    hCovMtxData->Sumw2();
    TH2D *hCovMtxJES = new TH2D("hCovMtxJES", "Covariance matrix for the JES source;N_{p_{T}^{ave}} x N_{#eta^{dijet}};N_{p_{T}^{ave}} x N_{#eta^{dijet}}", 
                                nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1,
                                nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1);
    hCovMtxJES->Sumw2();
    TH2D *hCovMtxJER = new TH2D("hCovMtxJER", "Covariance matrix for the JER source;N_{p_{T}^{ave}} x N_{#eta^{dijet}};N_{p_{T}^{ave}} x N_{#eta^{dijet}}", 
                                nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1,
                                nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1);
    hCovMtxJER->Sumw2();
    TH2D *hCovMtxPileup = new TH2D("hCovMtxPileup", "Covariance matrix for the pileup source;N_{p_{T}^{ave}} x N_{#eta^{dijet}};N_{p_{T}^{ave}} x N_{#eta^{dijet}}", 
                                   nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1,
                                   nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1);
    hCovMtxPileup->Sumw2();
    TH2D *hCovMtxPointing = new TH2D("hCovMtxPointing", "Covariance matrix for the pointing resolution source;N_{p_{T}^{ave}} x N_{#eta^{dijet}};N_{p_{T}^{ave}} x N_{#eta^{dijet}}", 
                                   nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1,
                                   nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1);
    hCovMtxPointing->Sumw2();
    TH2D *hCovMtxFull = new TH2D("hCovMtxFull", "Full covariance matrix;N_{p_{T}^{ave}} x N_{#eta^{dijet}};N_{p_{T}^{ave}} x N_{#eta^{dijet}}", 
                                 nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1,
                                 nPoints * ptDijetPtLow.size() + 2, -1, nPoints * ptDijetPtLow.size() + 1);
    hCovMtxFull->Sumw2();

    Double_t value{0.};
    // First loop over pT bins to construct covariance matrix
    for (UInt_t iPt{0}; iPt<ptDijetPtLow.size(); iPt++) {
        for (Int_t i{0}; i<nPoints; i++) {
            Double_t iErrJES = grJES[iPt]->GetErrorY(i);
            Double_t iErrJER = grJER[iPt]->GetErrorY(i);
            Double_t iErrPileup = grPileup[iPt]->GetErrorY(i);
            Double_t iErrPointing = grPointing[iPt]->GetErrorY(i);

            // Second loop over pT bins
            for (Int_t jPt{0}; jPt<ptDijetPtLow.size(); jPt++) {
                for (Int_t j{0}; j<nPoints; j++) {
                    Double_t jErrJES = grJES[jPt]->GetErrorY(j);
                    Double_t jErrJER = grJER[jPt]->GetErrorY(j);
                    Double_t jErrPileup = grPileup[jPt]->GetErrorY(j);
                    Double_t jErrPointing = grPointing[jPt]->GetErrorY(j);
                    hCovMtxJES->Fill(iPt * nPoints + i, jPt * nPoints + j, iErrJES * jErrJES );
                    hCovMtxJER->Fill(iPt * nPoints + i, jPt * nPoints + j, iErrJER * jErrJER );
                    hCovMtxPileup->Fill(iPt * nPoints + i, jPt * nPoints + j, iErrPileup * jErrPileup );
                    hCovMtxPointing->Fill(iPt * nPoints + i, jPt * nPoints + j, iErrPointing * jErrPointing );

                    if ( iPt == jPt && i == j ) {
                        Double_t iErrData = grData[iPt]->GetErrorY(i);
                        hCovMtxData->Fill(iPt * nPoints + i, jPt * nPoints + j, iErrData * iErrData );
                    }
                } // for (Int_t j{0}; j<nPoints; j++)
            } // for (Int_t jPt{0}; jPt<ptDijetPtLow.size(); jPt++)
        } // for (Int_t i{0}; i<nPoints; i++)
    } // for (UInt_t iPt{0}; iPt<ptDijetPtLow.size(); iPt++)

    hCovMtxFull->Add(hCovMtxJES);
    hCovMtxFull->Add(hCovMtxJER);
    hCovMtxFull->Add(hCovMtxPileup);
    hCovMtxFull->Add(hCovMtxPointing);

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    hCovMtxJES->Draw("colz");
    t.DrawLatexNDC(0.25, 0.93, "Covariance matrix (JES)" );
    //hCovMtxJES->GetZaxis()->SetRangeUser(-0.001, 0.001);
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_cov_jes_lab.pdf") );

    c->cd();
    setPadStyle();
    hCovMtxJER->Draw("colz");
    t.DrawLatexNDC(0.25, 0.93, "Covariance matrix (JER)" );
    //hCovMtxJER->GetZaxis()->SetRangeUser(-0.001, 0.001);
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_cov_jer_lab.pdf") );

    c->cd();
    setPadStyle();
    hCovMtxPileup->Draw("colz");
    t.DrawLatexNDC(0.25, 0.93, "Covariance matrix (pileup)" );
    //hCovMtxPileup->GetZaxis()->SetRangeUser(-0.001, 0.001);
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_cov_pileup_lab.pdf") );

    c->cd();
    setPadStyle();
    hCovMtxPointing->Draw("colz");
    t.DrawLatexNDC(0.25, 0.93, "Covariance matrix (pointing resolution)" );
    //hCovMtxPointing->GetZaxis()->SetRangeUser(-0.001, 0.001);
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_cov_pointing_lab.pdf") );

    c->cd();
    setPadStyle();
    hCovMtxFull->Draw("colz");
    t.DrawLatexNDC(0.25, 0.93, "Covariance matrix" );
    //hCovMtxFull->GetZaxis()->SetRangeUser(-0.001, 0.001);    
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_cov_full_lab.pdf") );

    c->cd();
    setPadStyle();
    hCovMtxData->Draw("colz");
    t.DrawLatexNDC(0.25, 0.93, "Covariance matrix (data)" );
    //hCovMtxData->GetZaxis()->SetRangeUser(-0.001, 0.001);
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_cov_diag_lab.pdf") );

    //
    // Build correlation matrix
    //
    TH2D *hCorrMtx = dynamic_cast<TH2D*>( hCovMtxFull->Clone() );
    hCorrMtx->SetNameTitle("hCorrMtx", "Correlation matrix");

    hCorrMtx->Add(hCovMtxData); // ?

    for (Int_t i=1; i<=hCovMtxFull->GetNbinsX(); i++) {
        for (Int_t j=1; j<=hCovMtxFull->GetNbinsY(); j++) {
            Double_t covVal = hCorrMtx->GetBinContent(i, j);
            Double_t iDiag = hCorrMtx->GetBinContent(i, i);
            Double_t jDiag = hCorrMtx->GetBinContent(j, j);
            hCorrMtx->SetBinContent(i, j, covVal / TMath::Sqrt(iDiag * jDiag) );
        } // for (Int_t j=1; j<=hCovMtxData->GetNbinsX(); j++)
    } // for (Int_t i=1; i<=hCovMtxData->GetNbinsX(); i++)

    c->cd();
    setPadStyle();
    hCorrMtx->Draw("colz");
    hCorrMtx->GetZaxis()->SetRangeUser(-1., 1.);
    t.DrawLatexNDC(0.25, 0.93, "Correlation matrix" );
    c->SaveAs( Form("freezeSyst/covMtx/pPb8160_dijet_corr_lab.pdf") );

}

//________________
void pPbFinalDijetEtaDistributions() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    //createValues();

    calculateCovMatrix();

    //plotFinalDistributions();
}