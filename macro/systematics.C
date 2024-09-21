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
    t.SetTextSize(0.06);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.5, 0.93, "pPb 174.56 nb^{-1} (8.16 TeV)");
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
    jeuPath += "/jes";
    jeuOldPath += "/jeu";
    jerPath += "/jer";
    pileupPath += "/pileup";
    pointingPath += "/pointing";
    dataPath += "/data";
    data2mcPath += "/data2mc";

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
void plotDifferentDirections(TFile *pbGoingFile, TFile *pGoingFile, TString date) {

    TString trigName = "MB";
    TString fName = pbGoingFile->GetName();
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

    TH3D *hPtEtaDphiPbGoing = (TH3D*)pbGoingFile->Get("hRecoDijetPtEtaDphiWeighted");
    //TH3D *hPtEtaDphiPbGoing = (TH3D*)pbGoingFile->Get("hRecoDijetPtEtaDphi");
    hPtEtaDphiPbGoing->SetName("hPtEtaDphiPbGoing");
    TH3D *hPtEtaDphiPGoing = (TH3D*)pGoingFile->Get("hRecoDijetPtEtaDphiWeighted");
    //TH3D *hPtEtaDphiPGoing = (TH3D*)pGoingFile->Get("hRecoDijetPtEtaDphi");
    hPtEtaDphiPGoing->SetName("hPtEtaDphiPGoing");

    TString frame = "lab";

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};

    // Styles
    Int_t pbGoingType{0};
    Int_t pGoingType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaPbGoing[ ptDijetBinLow.size() ];
    TH1D *hEtaPGoing[ ptDijetBinLow.size() ];
    TH1D *hEtaPb2PGoingRatio[ ptDijetBinLow.size() ];

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

    std::cout << "ptBins: " << ptBins << std::endl;

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        // Pb-going projection
        hEtaPbGoing[i] = projectEtaFrom3D(hPtEtaDphiPbGoing, Form("hEtaPbGoing_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaPbGoing[i] );
        set1DStyle(hEtaPbGoing[i], pbGoingType);

        // p-going projection
        hEtaPGoing[i] = projectEtaFrom3D(hPtEtaDphiPGoing, Form("hEtaPGoing_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaPGoing[i] );
        set1DStyle(hEtaPGoing[i], pGoingType);

        // Ratio of Pb-going to p-going
        hEtaPb2PGoingRatio[i] = (TH1D*)hEtaPbGoing[i]->Clone( Form("hEtaPb2PGoingRatio_%d", i) );
        make1DRatio(hEtaPb2PGoingRatio[i], hEtaPGoing[i], Form("hEtaPb2PGoingRatio_%d", i) );
        hEtaPb2PGoingRatio[i]->GetYaxis()->SetTitle("Pb-going / p-going");

        // Plot comparison
        canv->cd();
        setPadStyle();
        hEtaPbGoing[i]->Draw();
        hEtaPGoing[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaPbGoing[i], Form("Pb-going"), "p");
        leg->AddEntry(hEtaPGoing[i], Form("p-going"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/%s_pPb8160_etaDijet_dirComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios of Pb-going to p-going
        canv->cd();
        setPadStyle();
        hEtaPb2PGoingRatio[i]->Draw();
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hEtaPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/%s_pPb8160_etaDijet_dirRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hEtaPbGoing[i]->Draw();
        hEtaPGoing[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaPbGoing[i], Form("Pb-going"), "p");
        leg->AddEntry(hEtaPGoing[i], Form("p-going"), "p");
        leg->Draw();

        // Plot all ratios on one canvas
        cRat->cd(i+1);
        setPadStyle();
        hEtaPb2PGoingRatio[i]->Draw();
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hEtaPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
    } // for (Int_t i=0; i<ptDijetBinLow.size(); i++)

    cComp->SaveAs( Form("%s/%s_pPb8160_etaDijet_dirComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/%s_pPb8160_etaDijet_dirRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
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
        hBack->SetBinContent(i,  hOrig->GetBinContent(middle-i+1) );
        hBack->SetBinError(i, hOrig->GetBinError(middle-i+1) );

        hBFRatio->SetBinContent(i,  hOrig->GetBinContent(middle-i+1) );
        hBFRatio->SetBinError(i, hOrig->GetBinError(middle-i+1) );

        hForw->SetBinContent(i, hOrig->GetBinContent(middle+i) );
        hForw->SetBinError(i, hOrig->GetBinError(middle+i) );

        hFBRatio->SetBinContent(i, hOrig->GetBinContent(middle+i) );
        hFBRatio->SetBinError(i, hOrig->GetBinError(middle+i) );

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
        name += "_forward";
    }
    else {
        name += "_backward";
    }
    std::cout << Form("[Create FB] hName: %s", name.Data() ) << std::endl;
    
    std::vector<TH1D*> hFBRatio;
    for (Int_t i{0}; i<hBack.size(); i++) {
        if ( isForward ) {
            hFBRatio.push_back( dynamic_cast<TH1D*>( hForw[i]->Clone( Form("%s_%d", name.Data(), i) ) ) );
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

    std::cout << Form("[Dir lab] Direction histogram input histo name : %s", inputHistoName) << std::endl;
    std::cout << Form("[Dir lab] Direction histogram output histo name: %s", hName) << std::endl;
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
        rescaleEta(hEta[i]);
        hEtaForward.push_back( new TH1D( Form("%s_forward_%d", hName, i), Form("%s_forward_%d", hName, i), 50, 0., 5.) );
        hEtaBackward.push_back( new TH1D( Form("%s_backward_%d", hName, i), Form("%s_backward_%d", hName, i), 50, 0., 5.) );

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
    std::cout << Form("[Dir] Direction histogram input histo name : %s", inputHistoName) << std::endl;
    std::cout << Form("[Dir] Direction histogram output histo name: %s", name.Data()) << std::endl;
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
    std::cout << Form("[Full] Input histo name: %s", inputHistoName) << std::endl;
    std::cout << Form("[Full] Output histo name: %s", hName) << std::endl;
    TH3D* hPtEtaDphi = dynamic_cast<TH3D*>( inFile->Get( inputHistoName ) );
    // Reset histogram name
    hPtEtaDphi->SetName( hName );

    std::cout << Form("\t[Full] hPtEtaDphi entries: %f", hPtEtaDphi->Integral() ) << std::endl;

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
        rescaleEta(hEta[i]); set1DStyle(hEta[i], 2);
        if ( isCM ) {
            hEta[i]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hEta[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet}_{CM}");
        }
        else {
            hEta[i]->GetXaxis()->SetTitle("#eta^{dijet}");
            hEta[i]->GetYaxis()->SetTitle("dN/d#eta^{dijet}");
        }

        std::cout << Form("\t[Full] %s integral: %f entries: %f", hEta[i]->GetName(), hEta[i]->Integral(), hEta[i]->GetEntries() ) << std::endl;
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
        hMBEta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hMBEta[i]->Clone( Form("%s_mb2gen_%d_%s", hMBEta[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaForward[i]->Clone( Form("%s_mb2gen_forward_%d_%s", hMBEtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBackward[i]->Clone( Form("%s_mb2gen_backward_%d_%s", hMBEtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio[i]->Clone( Form("%s_mb2gen_fb_%d_%s", hMBEtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio[i]->Clone( Form("%s_mb2gen_bf_%d_%s", hMBEtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hJet60Eta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hJet60Eta[i]->Clone( Form("%s_jet60_2gen_%d_%s", hJet60Eta[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaForward[i]->Clone( Form("%s_jet60_2gen_forward_%d_%s", hJet60EtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward[i]->Clone( Form("%s_jet60_2gen_backward_%d_%s", hJet60EtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio[i]->Clone( Form("%s_jet60_2gen_fb_%d_%s", hJet60EtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio[i]->Clone( Form("%s_jet60_2gen_bf_%d_%s", hJet60EtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hJet80Eta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hJet80Eta[i]->Clone( Form("%s_jet80_2gen_%d_%s", hJet80Eta[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaForward[i]->Clone( Form("%s_jet80_2gen_forward_%d_%s", hJet80EtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBackward[i]->Clone( Form("%s_jet80_2gen_backward_%d_%s", hJet80EtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaFBRatio[i]->Clone( Form("%s_jet80_2gen_fb_%d_%s", hJet80EtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBFRatio[i]->Clone( Form("%s_jet80_2gen_bf_%d_%s", hJet80EtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hJet100Eta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hJet100Eta[i]->Clone( Form("%s_jet100_2gen_%d_%s", hJet100Eta[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaForward[i]->Clone( Form("%s_jet100_2gen_forward_%d_%s", hJet100EtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBackward[i]->Clone( Form("%s_jet100_2gen_backward_%d_%s", hJet100EtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaFBRatio[i]->Clone( Form("%s_jet100_2gen_fb_%d_%s", hJet100EtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBFRatio[i]->Clone( Form("%s_jet100_2gen_bf_%d_%s", hJet100EtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hEmbeddingEta2GenEtaRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEta[i]->Clone( Form("%s_embedding_2gen_%d_%s", hEmbeddingEta[i]->GetName(), i, frame.Data()) ) ) );
        hEmbeddingEtaForward2GenEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaForward[i]->Clone( Form("%s_embedding_2gen_forward_%d_%s", hEmbeddingEtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hEmbeddingEtaBackward2GenEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaBackward[i]->Clone( Form("%s_embedding_2gen_backward_%d_%s", hEmbeddingEtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hEmbeddingEtaFBRatio2GenEtaFBRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaFBRatio[i]->Clone( Form("%s_embedding_2gen_fb_%d_%s", hEmbeddingEtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hEmbeddingEtaBFRatio2GenEtaBFRatio.push_back( dynamic_cast<TH1D*>( hEmbeddingEtaBFRatio[i]->Clone( Form("%s_embedding_2gen_bf_%d_%s", hEmbeddingEtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        // Create histograms for ratios to embedding
        hMBEta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hMBEta[i]->Clone( Form("%s_mb2embedding_%d_%s", hMBEta[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaForward[i]->Clone( Form("%s_mb2embedding_forward_%d_%s", hMBEtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBackward[i]->Clone( Form("%s_mb2embedding_backward_%d_%s", hMBEtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hMBEtaFBRatio[i]->Clone( Form("%s_mb2embedding_fb_%d_%s", hMBEtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hMBEtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hMBEtaBFRatio[i]->Clone( Form("%s_mb2embedding_bf_%d_%s", hMBEtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hJet60Eta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hJet60Eta[i]->Clone( Form("%s_jet60_2embedding_%d_%s", hJet60Eta[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaForward[i]->Clone( Form("%s_jet60_2embedding_forward_%d_%s", hJet60EtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBackward[i]->Clone( Form("%s_jet60_2embedding_backward_%d_%s", hJet60EtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaFBRatio[i]->Clone( Form("%s_jet60_2embedding_fb_%d_%s", hJet60EtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hJet60EtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet60EtaBFRatio[i]->Clone( Form("%s_jet60_2embedding_bf_%d_%s", hJet60EtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hJet80Eta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hJet80Eta[i]->Clone( Form("%s_jet80_2embedding_%d_%s", hJet80Eta[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaForward[i]->Clone( Form("%s_jet80_2embedding_forward_%d_%s", hJet80EtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBackward[i]->Clone( Form("%s_jet80_2embedding_backward_%d_%s", hJet80EtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaFBRatio[i]->Clone( Form("%s_jet80_2embedding_fb_%d_%s", hJet80EtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hJet80EtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet80EtaBFRatio[i]->Clone( Form("%s_jet80_2embedding_bf_%d_%s", hJet80EtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );

        hJet100Eta2EmbeddingEtaRatio.push_back( dynamic_cast<TH1D*>( hJet100Eta[i]->Clone( Form("%s_jet100_2embedding_%d_%s", hJet100Eta[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaForward2EmbeddingEtaForwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaForward[i]->Clone( Form("%s_jet100_2embedding_forward_%d_%s", hJet100EtaForward[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaBackward2EmbeddingEtaBackwardRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBackward[i]->Clone( Form("%s_jet100_2embedding_backward_%d_%s", hJet100EtaBackward[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaFBRatio2EmbeddingEtaFBRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaFBRatio[i]->Clone( Form("%s_jet100_2embedding_fb_%d_%s", hJet100EtaFBRatio[i]->GetName(), i, frame.Data()) ) ) );
        hJet100EtaBFRatio2EmbeddingEtaBFRatio.push_back( dynamic_cast<TH1D*>( hJet100EtaBFRatio[i]->Clone( Form("%s_jet100_2embedding_bf_%d_%s", hJet100EtaBFRatio[i]->GetName(), i, frame.Data()) ) ) );


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
std::vector< std::vector<TH1D*> > createVectorOfVariationRatio(std::vector< std::vector<TH1D*> > hDefDist, 
                                                               std::vector< std::vector<TH1D*> > hVarDist,
                                                               Bool_t isCM = kFALSE, Int_t varType = 0) {

    TString frame = "lab";
    if (isCM) frame = "cm";

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
    std::vector< TH1D* > hEtaVar2DefRatio;
    std::vector< TH1D* > hEtaForwardVar2DefRatio;
    std::vector< TH1D* > hEtaBackwardVar2DefRatio;
    std::vector< TH1D* > hEtaFBRatioVar2DefRatio;
    std::vector< TH1D* > hEtaBFRatioVar2DefRatio;

    // Dijet pT selection
    Int_t ptBins = ptDijetBinLow.size();
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    TString varStr = "up2def_ratio";
    if (varType == 1) varStr = "down2def_ratio";

    Int_t style = 2;
    if (varType == 0) {
        style = 0;
    }
    else {
        style = 1;
    }

    // Loop over pT bins
    for (Int_t i{0}; i<ptBins; i++) {
            
        // Create histograms for ratios
        hEtaVar2DefRatio.push_back( dynamic_cast<TH1D*>( hVarEta[i]->Clone( Form("%s_%s_%d_%s", hVarEta[i]->GetName(), varStr.Data(), i, frame.Data()) ) ) );
        hEtaForwardVar2DefRatio.push_back( dynamic_cast<TH1D*>( hVarEtaForward[i]->Clone( Form("%s_%s_forward_%d_%s", hVarEtaForward[i]->GetName(), varStr.Data(), i, frame.Data()) ) ) );
        hEtaBackwardVar2DefRatio.push_back( dynamic_cast<TH1D*>( hVarEtaBackward[i]->Clone( Form("%s_%s_backward_%d_%s", hVarEtaBackward[i]->GetName(), varStr.Data(), i, frame.Data()) ) ) );
        hEtaFBRatioVar2DefRatio.push_back( dynamic_cast<TH1D*>( hVarEtaFBRatio[i]->Clone( Form("%s_%s_fb_%d_%s", hVarEtaFBRatio[i]->GetName(), varStr.Data(), i, frame.Data()) ) ) );
        hEtaBFRatioVar2DefRatio.push_back( dynamic_cast<TH1D*>( hVarEtaBFRatio[i]->Clone( Form("%s_%s_bf_%d_%s", hVarEtaBFRatio[i]->GetName(), varStr.Data(), i, frame.Data()) ) ) );

        // Build the ratios
        hEtaVar2DefRatio[i]->Divide( hDefEta[i] );
        hEtaForwardVar2DefRatio[i]->Divide( hDefEtaForward[i] );
        hEtaBackwardVar2DefRatio[i]->Divide( hDefEtaBackward[i] );
        hEtaFBRatioVar2DefRatio[i]->Divide( hDefEtaFBRatio[i] );
        hEtaBFRatioVar2DefRatio[i]->Divide( hDefEtaBFRatio[i] );

        // Set the style
        set1DStyle(hEtaVar2DefRatio[i], style);
        set1DStyle(hEtaForwardVar2DefRatio[i], style);
        set1DStyle(hEtaBackwardVar2DefRatio[i], style);
        set1DStyle(hEtaFBRatioVar2DefRatio[i], style);
        set1DStyle(hEtaBFRatioVar2DefRatio[i], style);
    } // for (Int_t i{0}; i<ptBins; i++)

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
                                                                Bool_t isCM = kFALSE) {

    TString frame = "lab";
    if (isCM) frame = "cm";                                                                    

    // Create vector of histograms to return
    std::vector< std::vector< TH1D* > > hRatios;

    // Vectors of default
    std::vector< std::vector<TH1D* > > hUp2DefRatios = createVectorOfVariationRatio(hDefDist, hUpDist, isCM, 0);
    std::vector< std::vector<TH1D* > > hDown2DefRatios = createVectorOfVariationRatio(hDefDist, hDownDist, isCM, 1);
    
    // Append the vectors
    appendVectorOfVectors(hRatios, hUp2DefRatios);
    appendVectorOfVectors(hRatios, hDown2DefRatios);

    return hRatios;
}



//________________
std::vector< std::vector<TH1D*> > createUpAndDown2Def(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hJeuMBUpEtaDist, std::vector< std::vector<TH1D*> > hJeuMBDownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hJet60EtaDist, std::vector< std::vector<TH1D*> > hJeuJet60UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet60DownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJeuJet80UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet80DownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hJet100EtaDist, std::vector< std::vector<TH1D*> > hJeuJet100UpEtaDist, std::vector< std::vector<TH1D*> > hJeuJet100DownEtaDist, 
                                                      std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hJerUpEtaDist, std::vector< std::vector<TH1D*> > hJerDownEtaDist,
                                                      std::vector< std::vector<TH1D*> > hMBGplusEtaDist, std::vector< std::vector<TH1D*> > hMBVtx1EtaDist,
                                                      std::vector< std::vector<TH1D*> > hJet60GplusEtaDist, std::vector< std::vector<TH1D*> > hJet60Vtx1EtaDist,
                                                      std::vector< std::vector<TH1D*> > hJet80GplusEtaDist, std::vector< std::vector<TH1D*> > hJet80Vtx1EtaDist,
                                                      std::vector< std::vector<TH1D*> > hJet100GplusEtaDist, std::vector< std::vector<TH1D*> > hJet100Vtx1EtaDist,
                                                      Bool_t isCM = kFALSE) {

    // Create vector of histograms to return
    std::vector< std::vector<TH1D*> > hRatios;

    // Retrieve histograms to divide
    std::vector< std::vector<TH1D*> > hMBEtaJeuRatios = createVectorOfUpAndDownRatios(hMBEtaDist, hJeuMBUpEtaDist, hJeuMBDownEtaDist, isCM);
    std::vector< std::vector<TH1D*> > hJet60EtaJeuRatios = createVectorOfUpAndDownRatios(hJet60EtaDist, hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, isCM);
    std::vector< std::vector<TH1D*> > hJet80EtaJeuRatios = createVectorOfUpAndDownRatios(hJet80EtaDist, hJeuJet80UpEtaDist, hJeuJet80DownEtaDist, isCM);
    std::vector< std::vector<TH1D*> > hJet100EtaJeuRatios = createVectorOfUpAndDownRatios(hJet100EtaDist, hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, isCM);

    std::vector< std::vector<TH1D*> > hEmbeddingJerRatios = createVectorOfUpAndDownRatios(hEmbeddingEtaDist, hJerUpEtaDist, hJerDownEtaDist, isCM);

    std::vector< std::vector<TH1D*> > hMBPileup = createVectorOfUpAndDownRatios(hMBEtaDist, hMBGplusEtaDist, hMBVtx1EtaDist, isCM);
    std::vector< std::vector<TH1D*> > hJet60Pileup = createVectorOfUpAndDownRatios(hJet60EtaDist, hJet60GplusEtaDist, hJet60Vtx1EtaDist, isCM);
    std::vector< std::vector<TH1D*> > hJet80Pileup = createVectorOfUpAndDownRatios(hJet80EtaDist, hJet80GplusEtaDist, hJet80Vtx1EtaDist, isCM);
    std::vector< std::vector<TH1D*> > hJet100Pileup = createVectorOfUpAndDownRatios(hJet100EtaDist, hJet100GplusEtaDist, hJet100Vtx1EtaDist, isCM);

    // Append the vectors
    appendVectorOfVectors(hRatios, hMBEtaJeuRatios);
    appendVectorOfVectors(hRatios, hJet60EtaJeuRatios);
    appendVectorOfVectors(hRatios, hJet80EtaJeuRatios);
    appendVectorOfVectors(hRatios, hJet100EtaJeuRatios);

    appendVectorOfVectors(hRatios, hEmbeddingJerRatios);

    appendVectorOfVectors(hRatios, hMBPileup);
    appendVectorOfVectors(hRatios, hJet60Pileup);
    appendVectorOfVectors(hRatios, hJet80Pileup);
    appendVectorOfVectors(hRatios, hJet100Pileup);

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

    std::cout << "[Main creator] Input histo name: " << inputHistoName << std::endl;
    std::cout << "[Main creator] Output histo name: " << hName << std::endl;
    std::cout << Form("\tFull: %s \n\tForward: %s \n\tBackward: %s", fullEtaHistName.Data(), forwardEtaHistName.Data(), backwardEtaHistName.Data() ) << std::endl;

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
void jeuSystematics(TF1 *jeuUp, TF1 *jeuDown, TH1D *hDef, TGraph* syst) {

    // Loop over the bins and estimate systematics for those that have entries
    for (Int_t i{1}; i<=hDef->GetNbinsX(); i++ ) {
        // if ( hDef->GetBinContent(i) == 0 ) continue;
        Double_t xVal = hDef->GetBinCenter(i);

        Double_t yUp = jeuUp->Eval( xVal );
        Double_t sign = ( (yUp - 1.) >= 0) ? 1 : -1;
        Double_t yDown = jeuDown->Eval( xVal );
        std::cout << Form("x: %3.2f up: %.3f down: %.3f ", xVal, yUp, yDown);
        Double_t sysYrel = ( TMath::Abs(yUp - 1.) + TMath::Abs(yDown - 1.) ) / 2;
        Double_t sysYabs = sysYrel * hDef->GetBinContent(i);

        // if ( TMath::Abs(hDef->GetBinContent(i)) < 0.00001 ) {
        //     sysYrel = 0.;
        // }

        std::cout << Form("bin content: %f syst [rel. uncrt.]: %.3f syst [abs. uncrt.]: %.8f\n", hDef->GetBinContent(i), sysYrel * 100., sysYabs );

        syst->SetPoint(i-1, xVal, sign * sysYrel );
    } // for (Int_t i{1}; i<=h->GetNbinsX(); i++ )
}

//________________
void jeuFBSystematics(TH1D *jeuUp, TH1D* jeuDown, TH1D *hDef, TGraph* syst) {
    // Loop over the bins and estimate systematics
    for (Int_t i{1}; i<=hDef->GetNbinsX(); i++ ) {
        Double_t xVal = hDef->GetBinCenter(i);
        Double_t yUp = jeuUp->GetBinContent(i);
        Double_t yDown{1};
        if (jeuDown != nullptr) {
            yDown = jeuDown->GetBinContent(i);
        }
        else {
            yDown = yUp;
        }
        Double_t sysYrel = ( TMath::Abs(yUp - 1.) + TMath::Abs(yDown - 1.) ) / 2;
        Double_t sysYabs = sysYrel * hDef->GetBinContent(i);
        syst->SetPoint(i-1, xVal, sysYrel );
    } // for (Int_t i{1}; i<=h->GetNbinsX(); i++ )
}

//________________
void plotJEU(TFile *defaultFile, TFile *jeuUpFile, TFile *jeuDownFile, TFile *dataFile, TString date, 
             Bool_t isCM = kFALSE, Bool_t drawFits = kTRUE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";
    TString histoName = "hRecoDijetPtEtaDphiWeighted";
    if ( isCM ) {
        histoName = "hRecoDijetPtEtaDphiCMWeighted";
    }
    TString trigName = "MB";
    TString fName = defaultFile->GetName();
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

    TH3D *hPtEtaDphiDef = (TH3D*)defaultFile->Get( histoName.Data() );
    hPtEtaDphiDef->SetName("hPtEtaDphiDef");
    TH3D *hPtEtaDphiUp = (TH3D*)jeuUpFile->Get( histoName.Data() );
    hPtEtaDphiUp->SetName("hPtEtaDphiUp");
    TH3D *hPtEtaDphiDown = (TH3D*)jeuDownFile->Get( histoName.Data() );
    hPtEtaDphiDown->SetName("hPtEtaDphiDown");
    
    TH3D *hPtEtaDphiData = (TH3D*)dataFile->Get( histoName.Data() );
    hPtEtaDphiData->SetName("hPtEtaDphiData");

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaDef[ ptDijetBinLow.size() ];
    TH1D *hEtaDefForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDefBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDefBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDefFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaUp[ ptDijetBinLow.size() ];
    TH1D *hEtaUpForward[ ptDijetBinLow.size() ];
    TH1D *hEtaUpBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaUpBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaUpFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaDown[ ptDijetBinLow.size() ];
    TH1D *hEtaDownForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDownBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDownBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDownFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaData[ ptDijetBinLow.size() ];
    TH1D *hEtaDataForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDataBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDataBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDataFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaRatioUp[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFBUp[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBFUp[ ptDijetBinLow.size() ];

    TH1D *hEtaRatioDown[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFBDown[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBFDown[ ptDijetBinLow.size() ];

    TF1 *fitRatioUp[ ptDijetBinLow.size() ];
    TF1 *fitRatioDown[ ptDijetBinLow.size() ];

    TGraph *grSyst[ ptDijetBinLow.size() ];
    TGraph *grSystFB[ ptDijetBinLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetBinLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBComp = new TCanvas("cFBComp", "cFBComp", sizeX, sizeY);
    cFBComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFComp = new TCanvas("cBFComp", "cBFComp", sizeX, sizeY);
    cBFComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBRat = new TCanvas("cFBRat", "cFBRat", sizeX, sizeY);
    cFBRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFRat = new TCanvas("cBFRat", "cBFRat", sizeX, sizeY);
    cBFRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TString fbTitle = "#eta^{dijet}; #frac{ dN/d#eta^{forward} }{ dN/d#eta^{backward} }";
    TString bfTitle = "#eta^{dijet}; #frac{ dN/d#eta^{backward} }{ dN/d#eta^{forward} }";
    if ( isCM ) {
        fbTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (forward) }{ dN/d#eta^{dijet}_{CM} (backward) }";
        bfTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (backward) }{ dN/d#eta^{dijet}_{CM} (forward) }";
    }

    //
    // Loop over dijet pT bins
    //
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        //
        // Read histograms for the real DATA
        //
        hEtaData[i] = projectEtaFrom3D(hPtEtaDphiData, Form("hEtaData_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaData[i] );
        set1DStyle(hEtaData[i], defType);
        hEtaDataForward[i] = new TH1D( Form("hEtaDataForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDataBackward[i] = new TH1D( Form("hEtaDataBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDataBFRatio[i] = new TH1D( Form("hEtaDataBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDataFBRatio[i] = new TH1D( Form("hEtaDataFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDataForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDataForward[i]->Sumw2();  set1DStyle(hEtaDataForward[i], upType);
        hEtaDataBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDataBackward[i]->Sumw2(); set1DStyle(hEtaDataBackward[i], downType);
        hEtaDataBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDataBFRatio[i]->Sumw2();  set1DStyle(hEtaDataBFRatio[i], defType);
        hEtaDataFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDataFBRatio[i]->Sumw2();  set1DStyle(hEtaDataFBRatio[i], defType);
        makeDijetFB(hEtaData[i], hEtaDataBackward[i], hEtaDataForward[i], hEtaDataBFRatio[i], hEtaDataFBRatio[i]);

        //
        // Nominal (default) JEU selection
        //
        hEtaDef[i] = projectEtaFrom3D(hPtEtaDphiDef, Form("hEtaDef_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaDef[i] );
        set1DStyle(hEtaDef[i], defType);
        hEtaDefForward[i] = new TH1D( Form("hEtaDefForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDefBackward[i] = new TH1D( Form("hEtaDefBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDefBFRatio[i] = new TH1D( Form("hEtaDefBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDefFBRatio[i] = new TH1D( Form("hEtaDefFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDefForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefForward[i]->Sumw2();  set1DStyle(hEtaDefForward[i], upType);
        hEtaDefBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDefBackward[i]->Sumw2(); set1DStyle(hEtaDefBackward[i], downType);
        hEtaDefBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefBFRatio[i]->Sumw2();  set1DStyle(hEtaDefBFRatio[i], defType);
        hEtaDefFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefFBRatio[i]->Sumw2();  set1DStyle(hEtaDefFBRatio[i], defType);
        makeDijetFB(hEtaDef[i], hEtaDefBackward[i], hEtaDefForward[i], hEtaDefBFRatio[i], hEtaDefFBRatio[i]);

        //
        // JEU up selection
        //
        hEtaUp[i] = projectEtaFrom3D(hPtEtaDphiUp, Form("hEtaUp_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaUp[i] );
        set1DStyle(hEtaUp[i], upType);
        hEtaUpForward[i] = new TH1D( Form("hEtaUpForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaUpBackward[i] = new TH1D( Form("hEtaUpBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaUpBFRatio[i] = new TH1D( Form("hEtaUpBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaUpFBRatio[i] = new TH1D( Form("hEtaUpFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaUpForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpForward[i]->Sumw2();  set1DStyle(hEtaUpForward[i], upType);
        hEtaUpBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaUpBackward[i]->Sumw2(); set1DStyle(hEtaUpBackward[i], downType);
        hEtaUpBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpBFRatio[i]->Sumw2();  set1DStyle(hEtaUpBFRatio[i], upType);
        hEtaUpFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpFBRatio[i]->Sumw2();  set1DStyle(hEtaUpFBRatio[i], upType);
        makeDijetFB(hEtaUp[i], hEtaUpBackward[i], hEtaUpForward[i], hEtaUpBFRatio[i], hEtaUpFBRatio[i]);   

        //
        // JEU down selection
        //
        hEtaDown[i] = projectEtaFrom3D(hPtEtaDphiDown, Form("hEtaDown_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaDown[i] );
        set1DStyle(hEtaDown[i], downType);
        hEtaDownForward[i] = new TH1D( Form("hEtaDownForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDownBackward[i] = new TH1D( Form("hEtaDownBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDownBFRatio[i] = new TH1D( Form("hEtaDownBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDownFBRatio[i] = new TH1D( Form("hEtaDownFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDownForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownForward[i]->Sumw2();  set1DStyle(hEtaDownForward[i], downType);
        hEtaDownBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDownBackward[i]->Sumw2(); set1DStyle(hEtaDownBackward[i], downType);
        hEtaDownBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownBFRatio[i]->Sumw2();  set1DStyle(hEtaDownBFRatio[i], downType);
        hEtaDownFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownFBRatio[i]->Sumw2();  set1DStyle(hEtaDownFBRatio[i], downType);
        makeDijetFB(hEtaDown[i], hEtaDownBackward[i], hEtaDownForward[i], hEtaDownBFRatio[i], hEtaDownFBRatio[i]);

        //
        // Ratio of JEU up to default ratio
        //
        hEtaRatioUp[i] = dynamic_cast<TH1D*> ( hEtaUp[i]->Clone( Form("hEtaRatioUp_%d", i) ) );
        make1DRatio(hEtaRatioUp[i], hEtaDef[i], Form("hEtaRatioUp_%d", i), upType );
        hEtaRatioUp[i]->GetYaxis()->SetTitle("JES / Default");

        hEtaRatioBFUp[i] = dynamic_cast<TH1D*>( hEtaUpBFRatio[i]->Clone( Form("hEtaRatioBFUp_%d", i) ) );
        hEtaRatioBFUp[i]->Divide( hEtaRatioBFUp[i], hEtaDefBFRatio[i], 1., 1. );
        set1DStyle(hEtaRatioBFUp[i], upType);
        hEtaRatioFBUp[i] = dynamic_cast<TH1D*>( hEtaUpFBRatio[i]->Clone( Form("hEtaRatioFBUp_%d", i) ) );
        hEtaRatioFBUp[i]->Divide( hEtaRatioFBUp[i], hEtaDefFBRatio[i], 1., 1. );
        set1DStyle(hEtaRatioFBUp[i], upType);    

        fitRatioUp[i] = new TF1(Form("fitRatioUp_%d", i), "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", -5., 5.);
        // fitRatioUp[i] = new TF1(Form("fitRatioUp_%d", i), "[0]+[1]*x+[2]*x*x", -5., 5.);
        // fitRatioUp[i] = new TF1(Form("fitRatioUp_%d", i), "[0]+[1]*(x-0.465)*(x-0.465)", -5., 5.);
        fitRatioUp[i]->SetParameters(0.986138, -0.0026037, 0.0119519);
        fitRatioUp[i]->SetParLimits(2, 0., 1e6);
        fitRatioUp[i]->SetLineColor(kRed);
        fitRatioUp[i]->SetLineWidth(2);
        hEtaRatioUp[i]->Fit(Form("fitRatioUp_%d", i), "MRE0");

        // Ratio of JEU down to default ratio
        hEtaRatioDown[i] = dynamic_cast<TH1D*>( hEtaDown[i]->Clone( Form("hEtaRatioDown_%d", i) ) );
        make1DRatio(hEtaRatioDown[i], hEtaDef[i], Form("hEtaRatioDown_%d", i), downType );
        hEtaRatioDown[i]->GetYaxis()->SetTitle("JES / Default");

        hEtaRatioBFDown[i] = dynamic_cast<TH1D*>( hEtaDownBFRatio[i]->Clone( Form("hEtaRatioBFDown_%d", i) ) );
        hEtaRatioBFDown[i]->Divide( hEtaRatioBFDown[i], hEtaDefBFRatio[i], 1., 1. );
        set1DStyle(hEtaRatioBFDown[i], downType);
        hEtaRatioFBDown[i] = dynamic_cast<TH1D*>( hEtaDownFBRatio[i]->Clone( Form("hEtaRatioFBDown_%d", i) ) );
        hEtaRatioFBDown[i]->Divide( hEtaRatioFBDown[i], hEtaDefFBRatio[i], 1., 1. );     
        set1DStyle(hEtaRatioFBDown[i], downType);

        fitRatioDown[i] = new TF1(Form("fitRatioDown_%d", i), "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", -5., 5.);
        // fitRatioDown[i] = new TF1(Form("fitRatioDown_%d", i), "[0]+[1]*x+[2]*x*x", -5., 5.);
        // fitRatioDown[i] = new TF1(Form("fitRatioDown_%d", i), "[0]+[1]*(x-0.465)*(x-0.465)", -5., 5.);
        fitRatioDown[i]->SetParameters(1.01234, 0.00637344, -0.0104644);
        fitRatioDown[i]->SetParLimits(2, -1e6, 0);
        fitRatioDown[i]->SetLineColor(kBlue);
        fitRatioDown[i]->SetLineWidth(2);
        hEtaRatioDown[i]->Fit(Form("fitRatioDown_%d", i), "MRE0");


        // Calculate systematic uncertainty
        grSyst[i] = new TGraph();
        jeuSystematics(fitRatioUp[i], fitRatioDown[i], hEtaData[i], grSyst[i]);
        grSyst[i]->SetName( Form("grSystJEU_%d", i) );

        grSystFB[i] = new TGraph();
        jeuFBSystematics(hEtaUpFBRatio[i], hEtaDownFBRatio[i], hEtaDataFBRatio[i], grSystFB[i]);
        grSystFB[i]->SetName( Form("grSystJEU_FB_%d", i) );

        // std::cout << "Number of points in graph: " << grSyst[i]->GetN() << std::endl;
        // grSyst[i]->Print();

        // Write systematic uncertainty to ASCII file
        std::ofstream outFile( Form("%s/jeu/%s_pPb8160_etaDijet_jeuSyst_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFile.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSyst[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSyst[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFile << x << " " << y << std::endl;
            }
            outFile.close();
        }
        else {
            std::cerr << "[JeuSyst] Unable to open file for writing" << std::endl;
        }

        // Write FB systematic uncertainty to ASCII file
        std::ofstream outFileFB( Form("%s/jeu/%s_pPb8160_etaDijet_jeuSyst_FB_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFileFB.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSystFB[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSystFB[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFileFB << x << " " << y << std::endl;
            }
            outFileFB.close();
        }
        else {
            std::cerr << "[JeuSyst] Unable to open file for writing" << std::endl;
        }

        // Write data points to ASCII file
        std::ofstream outFile2( Form("%s/data/%s_pPb8160_etaDijet_data_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFile2.is_open() ) {
            for (Int_t iPoint = 1; iPoint <= hEtaDef[i]->GetNbinsX(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0}, xErr{0}, yErr{0}; 
                x = hEtaDef[i]->GetXaxis()->GetBinCenter(iPoint);
                y = hEtaDef[i]->GetBinContent(iPoint);
                xErr = hEtaDef[i]->GetXaxis()->GetBinWidth(iPoint) / 2;
                yErr = hEtaDef[i]->GetBinError(iPoint);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFile2 << x << " " << y << " " << xErr << " " << yErr << std::endl;
            }
            outFile2.close();
        }
        else {
            std::cerr << "[DATA] Unable to open file for writing" << std::endl;
        }

        std::ofstream outFile3( Form("%s/data/%s_pPb8160_etaDijet_data_fb_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if (outFile3.is_open()) {
            for (Int_t iPoint = 1; iPoint <= hEtaDefFBRatio[i]->GetNbinsX(); iPoint++) {
                Double_t x{0}, y{0}, xErr{0}, yErr{0};
                x = hEtaDefFBRatio[i]->GetXaxis()->GetBinCenter(iPoint);
                y = hEtaDefFBRatio[i]->GetBinContent(iPoint);
                xErr = hEtaDefFBRatio[i]->GetXaxis()->GetBinWidth(iPoint) / 2;
                yErr = hEtaDefFBRatio[i]->GetBinError(iPoint);
                outFile3 << x << " " << y << " " << xErr << " " << yErr << std::endl;
            }
            outFile3.close();
        }
        else {
            std::cerr << "[DATA] Unable to open file for writing" << std::endl;
        }

        //
        // Plot comparison eta distributions
        //
        canv->cd();
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("Default"), "p");
        leg->AddEntry(hEtaUp[i], Form("JES+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JES-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot BF comparison
        //
        canv->cd();
        setPadFBStyle();
        hEtaDefBFRatio[i]->Draw();
        hEtaUpBFRatio[i]->Draw("same");
        hEtaDownBFRatio[i]->Draw("same");
        hEtaDefBFRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefBFRatio[i]->GetYaxis()->SetRangeUser(0., 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.7, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefBFRatio[i], Form("Backward/Forward"), "p");
        leg->AddEntry(hEtaUpBFRatio[i], Form("JES+"), "p");
        leg->AddEntry(hEtaDownBFRatio[i], Form("JES-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot FB comparison
        //
        canv->cd();
        setPadFBStyle();
        hEtaDefFBRatio[i]->Draw();
        hEtaUpFBRatio[i]->Draw("same");
        hEtaDownFBRatio[i]->Draw("same");
        hEtaDefFBRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefFBRatio[i]->GetYaxis()->SetRangeUser(0.8, 1.8);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.55, 0.65, 0.65, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefBFRatio[i], Form("Forward/Backward"), "p");
        leg->AddEntry(hEtaUpBFRatio[i], Form("JES+"), "p");
        leg->AddEntry(hEtaDownBFRatio[i], Form("JES-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Ratio of FB up to default
        //
        canv->cd(i+1);
        setPadFBStyle();
        hEtaRatioFBUp[i]->Draw();
        hEtaRatioFBDown[i]->Draw("same");
        hEtaRatioFBUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioFBUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioFBUp[i]->GetYaxis()->SetTitle("JES/Default forward/backward");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.3, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioFBUp[i], Form("JES+/def"), "p");
        leg->AddEntry(hEtaRatioFBDown[i], Form("JES-/def"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Ratio of BF up to default
        //
        canv->cd(i+1);
        setPadFBStyle();
        hEtaRatioBFUp[i]->Draw();
        hEtaRatioBFDown[i]->Draw("same");
        hEtaRatioBFUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioBFUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioBFUp[i]->GetYaxis()->SetTitle("JES/Default backward/forward");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.3, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioBFUp[i], Form("JES+/def"), "p");
        leg->AddEntry(hEtaRatioBFDown[i], Form("JES-/def"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot ratios of full eta distributions
        //
        canv->cd();
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        if (drawFits) {
            fitRatioUp[i]->SetLineColor(kRed);
            fitRatioDown[i]->SetLineColor(kBlue);
            fitRatioUp[i]->Draw("same");
            fitRatioDown[i]->Draw("same");
        }
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("JES+ / Default"), "p");
        leg->AddEntry(hEtaRatioDown[i], Form("JES- / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot comparisons of full dijet pseudorapidity on one canvas
        //
        cComp->cd(i+1);
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("Default"), "p");
        leg->AddEntry(hEtaUp[i], Form("JES+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JES-"), "p");
        leg->Draw();

        //
        // Plot comparisons of backward to forward ratios on one canvas
        //
        cBFComp->cd(i+1);
        setPadFBStyle();
        hEtaDefBFRatio[i]->Draw();
        hEtaUpBFRatio[i]->Draw("same");
        hEtaDownBFRatio[i]->Draw("same");
        hEtaDefBFRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefBFRatio[i]->GetYaxis()->SetRangeUser(0., 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.35, 0.65, 0.45, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefBFRatio[i], Form("Backward/Forward"), "p");
        leg->AddEntry(hEtaUpBFRatio[i], Form("JES+"), "p");
        leg->AddEntry(hEtaDownBFRatio[i], Form("JES-"), "p");
        leg->Draw();

        //
        // Plot comparisons of forward to backward ratios on one canvas
        //
        cFBComp->cd(i+1);
        setPadFBStyle();
        hEtaDefFBRatio[i]->Draw();
        hEtaUpFBRatio[i]->Draw("same");
        hEtaDownFBRatio[i]->Draw("same");
        hEtaDefFBRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefFBRatio[i]->GetYaxis()->SetRangeUser(0.8, 1.8);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.6, 0.65, 0.7, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefFBRatio[i], Form("Forward/Backward"), "p");
        leg->AddEntry(hEtaUpFBRatio[i], Form("JES+"), "p");
        leg->AddEntry(hEtaDownFBRatio[i], Form("JES-"), "p");
        leg->Draw();

        //
        // Plot ratios of full dijet pseudorapidity on one canvas
        //
        cRat->cd(i+1);
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        if (drawFits) {
            fitRatioUp[i]->Draw("same");
            fitRatioDown[i]->Draw("same");
        }
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("JES+ / Default"), "p");
        leg->AddEntry(hEtaRatioDown[i], Form("JES- / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();

        //
        // Ratio of FB up to default on one canvas
        //
        cFBRat->cd(i+1);
        setPadFBStyle();
        hEtaRatioFBUp[i]->Draw();
        hEtaRatioFBDown[i]->Draw("same");
        hEtaRatioFBUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioFBUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioFBUp[i]->GetYaxis()->SetTitle("JES/Default forward/backward");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.3, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioFBUp[i], Form("JES+/def"), "p");
        leg->AddEntry(hEtaRatioFBDown[i], Form("JES-/def"), "p");
        leg->Draw();

        //
        // Ratio of BF up to default on one canvas
        //
        cBFRat->cd(i+1);
        setPadFBStyle();
        hEtaRatioBFUp[i]->Draw();
        hEtaRatioBFDown[i]->Draw("same");
        hEtaRatioBFUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioBFUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioBFUp[i]->GetYaxis()->SetTitle("JES/Default backward/FBforward");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.3, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioBFUp[i], Form("JES+/def"), "p");
        leg->AddEntry(hEtaRatioBFDown[i], Form("JES-/def"), "p");
        leg->Draw();
    } // for (Int_t i=0; i<ptDijetBinLow.size(); i++)

    cComp->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cComp->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_jeuRat_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cBFComp->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFComp->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cFBComp->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cFBComp->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_bf_jeuRat_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cFBRat->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cFBRat->SaveAs( Form("%s/jeu/%s_pPb8160_etaDijet_fb_jeuRat_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
}

//________________
void jerSystematics(TF1 *jerUp, TF1 *jerDown, TH1D *hDef, TGraph* syst) {

    // Loop over the bins and estimate systematics for those that have entries
    for (Int_t i{1}; i<=hDef->GetNbinsX(); i++ ) {
        // if ( hDef->GetBinContent(i) == 0 ) continue;
        Double_t xVal = hDef->GetBinCenter(i);

        Double_t yUp = jerUp->Eval( xVal );
        Double_t sign = ( (yUp - 1.) >= 0 ) ? 1: -1;
        Double_t yDown = jerDown->Eval( xVal );
        std::cout << Form("x: %3.2f up: %.3f down: %.3f ", xVal, yUp, yDown);
        Double_t sysYrel = ( TMath::Abs(yUp - 1.) + TMath::Abs(yDown - 1.) ) / 2;
        Double_t sysYabs = sysYrel * hDef->GetBinContent(i);
        // if ( TMath::Abs(hDef->GetBinContent(i)) < 0.00001 ) {
        //     sysYrel = 0.;
        // }

        std::cout << Form("bin content: %f syst [perc.]: %.3f syst [abs. val.]: %.8f\n", hDef->GetBinContent(i), sysYrel * 100., sysYabs );

        syst->SetPoint(i-1, xVal, sign * sysYrel );
    } // for (Int_t i{1}; i<=h->GetNbinsX(); i++ )
    std::cout << "JER up chi2/ndf: " << jerUp->GetChisquare() / jerUp->GetNDF() 
              << " JER down chi2/ndf: " << jerDown->GetChisquare() / jerDown->GetNDF()
              << std::endl;
}

//________________
void updateJERRatio(TH1D *h) {
    for (Int_t i{1}; i<=h->GetNbinsX(); i++) {
        if (h->GetBinContent(i) != 0 ) {
            Double_t val = h->GetBinContent(i);
            val = TMath::Abs(val - 1.) + 1.;
            h->SetBinContent(i, val);
        }
    } // for (Int_t i{1}; i<=h->GetNBinsX(); i++)
}

//________________
void plotJER(TFile *defaultFile, TFile *jerUpFile, TFile *jerDownFile, TString date, 
             Bool_t isCM = kFALSE, Bool_t drawFits = kTRUE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";
    TString histoName = "hGenDijetPtEtaDphiWeighted";
    if ( isCM ) {
        histoName = "hGenDijetPtEtaDphiCMWeighted";
    }
    TString trigName = "MB";
    TString fName = defaultFile->GetName();
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

    TH3D *hPtEtaDphiDef = (TH3D*)defaultFile->Get( histoName.Data() );
    hPtEtaDphiDef->SetName("hPtEtaDphiDef");
    TH3D *hPtEtaDphiUp = (TH3D*)jerUpFile->Get( histoName.Data() );
    hPtEtaDphiUp->SetName("hPtEtaDphiUp");
    TH3D *hPtEtaDphiDown = (TH3D*)jerDownFile->Get( histoName.Data() );
    hPtEtaDphiDown->SetName("hPtEtaDphiDown");

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    // std::vector<Int_t> ptDijetBinLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55 , 3};
    // std::vector<Int_t> ptDijetBinHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194, 194};

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaDef[ ptDijetBinLow.size() ];
    TH1D *hEtaDefForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDefBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDefBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDefFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaUp[ ptDijetBinLow.size() ];
    TH1D *hEtaUpForward[ ptDijetBinLow.size() ];
    TH1D *hEtaUpBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaUpBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaUpFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaDown[ ptDijetBinLow.size() ];
    TH1D *hEtaDownForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDownBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDownBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDownFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaRatioUp[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFBUp[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBFUp[ ptDijetBinLow.size() ];

    TH1D *hEtaRatioDown[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFBDown[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBFDown[ ptDijetBinLow.size() ];

    TF1 *fitRatioUp[ ptDijetBinLow.size() ];
    TF1 *fitRatioDown[ ptDijetBinLow.size() ];

    TGraph *grSyst[ ptDijetBinLow.size() ];
    TGraph *grSystFB[ ptDijetBinLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetBinLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBComp = new TCanvas("cFBComp", "cFBComp", sizeX, sizeY);
    cFBComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFComp = new TCanvas("cBFComp", "cBFComp", sizeX, sizeY);
    cBFComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBRat = new TCanvas("cFBRat", "cFBRat", sizeX, sizeY);
    cFBRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFRat = new TCanvas("cBFRat", "cBFRat", sizeX, sizeY);
    cBFRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TString fbTitle = "#eta^{dijet}; #frac{ dN/d#eta^{forward} }{ dN/d#eta^{backward} }";
    TString bfTitle = "#eta^{dijet}; #frac{ dN/d#eta^{backward} }{ dN/d#eta^{forward} }";
    if ( isCM ) {
        fbTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (forward) }{ dN/d#eta^{dijet}_{CM} (backward) }";
        bfTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (backward) }{ dN/d#eta^{dijet}_{CM} (forward) }";
    }

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        //
        // JER default selection
        //
        hEtaDef[i] = projectEtaFrom3D(hPtEtaDphiDef, Form("hEtaDef_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaDef[i] );
        set1DStyle(hEtaDef[i], defType);
        hEtaDefForward[i] = new TH1D( Form("hEtaDefForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDefBackward[i] = new TH1D( Form("hEtaDefBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDefBFRatio[i] = new TH1D( Form("hEtaDefBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDefFBRatio[i] = new TH1D( Form("hEtaDefFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDefForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefForward[i]->Sumw2();  set1DStyle(hEtaDefForward[i], upType);
        hEtaDefBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDefBackward[i]->Sumw2(); set1DStyle(hEtaDefBackward[i], downType);
        hEtaDefBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefBFRatio[i]->Sumw2();  set1DStyle(hEtaDefBFRatio[i], defType);
        hEtaDefFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefFBRatio[i]->Sumw2();  set1DStyle(hEtaDefFBRatio[i], defType);
        makeDijetFB(hEtaDef[i], hEtaDefBackward[i], hEtaDefForward[i], hEtaDefBFRatio[i], hEtaDefFBRatio[i]);

        //
        // JER up selection
        //
        hEtaUp[i] = projectEtaFrom3D(hPtEtaDphiUp, Form("hEtaUp_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaUp[i] );
        set1DStyle(hEtaUp[i], upType);
        hEtaUpForward[i] = new TH1D( Form("hEtaUpForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaUpBackward[i] = new TH1D( Form("hEtaUpBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaUpBFRatio[i] = new TH1D( Form("hEtaUpBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaUpFBRatio[i] = new TH1D( Form("hEtaUpFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaUpForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpForward[i]->Sumw2();  set1DStyle(hEtaUpForward[i], upType);
        hEtaUpBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaUpBackward[i]->Sumw2(); set1DStyle(hEtaUpBackward[i], downType);
        hEtaUpBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpBFRatio[i]->Sumw2();  set1DStyle(hEtaUpBFRatio[i], upType);
        hEtaUpFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpFBRatio[i]->Sumw2();  set1DStyle(hEtaUpFBRatio[i], upType);
        makeDijetFB(hEtaUp[i], hEtaUpBackward[i], hEtaUpForward[i], hEtaUpBFRatio[i], hEtaUpFBRatio[i]);

        //
        // JER down selection
        //
        hEtaDown[i] = projectEtaFrom3D(hPtEtaDphiDown, Form("hEtaDown_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaDown[i] );
        set1DStyle(hEtaDown[i], downType);
        hEtaDownForward[i] = new TH1D( Form("hEtaDownForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDownBackward[i] = new TH1D( Form("hEtaDownBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDownBFRatio[i] = new TH1D( Form("hEtaDownBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDownFBRatio[i] = new TH1D( Form("hEtaDownFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDownForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownForward[i]->Sumw2();  set1DStyle(hEtaDownForward[i], downType);
        hEtaDownBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDownBackward[i]->Sumw2(); set1DStyle(hEtaDownBackward[i], downType);
        hEtaDownBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownBFRatio[i]->Sumw2();  set1DStyle(hEtaDownBFRatio[i], downType);
        hEtaDownFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownFBRatio[i]->Sumw2();  set1DStyle(hEtaDownFBRatio[i], downType);
        makeDijetFB(hEtaDown[i], hEtaDownBackward[i], hEtaDownForward[i], hEtaDownBFRatio[i], hEtaDownFBRatio[i]);

        // Ratio of JER up to default ratio
        hEtaRatioUp[i] = dynamic_cast<TH1D*>( hEtaUp[i]->Clone( Form("hEtaRatioUp_%d", i) ) );
        make1DRatio(hEtaRatioUp[i], hEtaDef[i], Form("hEtaRatioUp_%d", i), upType );
        hEtaRatioUp[i]->GetYaxis()->SetTitle("JER / Default");
        updateJERRatio( hEtaRatioUp[i] );
        hEtaRatioFBUp[i] = dynamic_cast<TH1D*>( hEtaUpFBRatio[i]->Clone( Form("hEtaRatioFBUp_%d", i) ) );
        // make1DRatio(hEtaRatioFBUp[i], hEtaDefFBRatio[i], Form("hEtaRatioFBUp_%d", i), upType);
        hEtaRatioBFUp[i] = dynamic_cast<TH1D*>( hEtaUpBFRatio[i]->Clone( Form("hEtaRatioBFUp_%d", i) ) );
        hEtaRatioFBUp[i]->Divide(hEtaRatioFBUp[i], hEtaDefFBRatio[i], 1., 1.);
        set1DStyle(hEtaRatioBFUp[i], upType);

        // make1DRatio(hEtaRatioBFUp[i], hEtaDefBFRatio[i], Form("hEtaRatioBFUp_%d", i), upType);
        hEtaRatioBFUp[i]->Divide(hEtaRatioBFUp[i], hEtaDefBFRatio[i], 1., 1.);
        set1DStyle(hEtaRatioBFUp[i], upType);
        hEtaRatioFBUp[i]->GetYaxis()->SetTitle("JER / Default");
        hEtaRatioBFUp[i]->GetYaxis()->SetTitle("JER / Default");

        // Ratio of JER down to default ratio
        hEtaRatioDown[i] = dynamic_cast<TH1D*>( hEtaDown[i]->Clone( Form("hEtaRatioDown_%d", i) ) );
        make1DRatio(hEtaRatioDown[i], hEtaDef[i], Form("hEtaRatioDown_%d", i), downType );
        hEtaRatioDown[i]->GetYaxis()->SetTitle("JER / Default");
        updateJERRatio( hEtaRatioDown[i] );
        hEtaRatioFBDown[i] = dynamic_cast<TH1D*>( hEtaDownFBRatio[i]->Clone( Form("hEtaRatioFBDown_%d", i) ) );
        // make1DRatio(hEtaRatioFBDown[i], hEtaDefFBRatio[i], Form("hEtaRatioFBDown_%d", i), downType);
        hEtaRatioFBDown[i]->Divide(hEtaRatioFBDown[i], hEtaDefFBRatio[i], 1., 1.);
        set1DStyle(hEtaRatioFBDown[i], downType);

        hEtaRatioBFDown[i] = dynamic_cast<TH1D*>( hEtaDownBFRatio[i]->Clone( Form("hEtaRatioBFDown_%d", i) ) );
        // make1DRatio(hEtaRatioBFDown[i], hEtaDefBFRatio[i], Form("hEtaRatioBFDown_%d", i), downType);
        hEtaRatioBFDown[i]->Divide(hEtaRatioBFDown[i], hEtaDefBFRatio[i], 1., 1.);
        set1DStyle(hEtaRatioBFDown[i], downType);
        hEtaRatioFBDown[i]->GetYaxis()->SetTitle("JER / Default");
        hEtaRatioBFDown[i]->GetYaxis()->SetTitle("JER / Default");

        // Make fit functions
        fitRatioUp[i] = new TF1(Form("fitRatioUp_%d",i), "[0]+[1]*x + [2]*x*x", -3.5, 3.5);
        // fitRatioUp[i] = new TF1(Form("fitRatioUp_%d",i), "[0]", -5., 5.);
        // fitRatioUp[i]->SetParameters(0., 0.0001);
        fitRatioUp[i]->SetParLimits(2, 0.000001, 1e6);
        fitRatioUp[i]->SetParameters(0.);
        fitRatioUp[i]->SetLineColor(kRed);
        fitRatioUp[i]->SetLineWidth(2);

        fitRatioDown[i] = new TF1(Form("fitRatioDown_%d",i), "[0]+[1]*x +[2]*x*x", -3.5, 3.5);
        // fitRatioDown[i] = new TF1(Form("fitRatioDown_%d",i), "[0]", -5., 5.);
        // fitRatioDown[i]->SetParameters(0., 0.0001);
        fitRatioDown[i]->SetParLimits(2, 0.000001, 1e6);
        fitRatioDown[i]->SetParameters(0.);
        fitRatioDown[i]->SetLineColor(kBlue);
        fitRatioDown[i]->SetLineWidth(2);

        // Perform fits
        hEtaRatioUp[i]->Fit(Form("fitRatioUp_%d",i), "MRE0");
        hEtaRatioDown[i]->Fit(Form("fitRatioDown_%d",i), "MRE0");
        std::cout << "chi2/ndf up: " << fitRatioUp[i]->GetChisquare() / fitRatioUp[i]->GetNDF() << std::endl;
        std::cout << "chi2/ndf down: " << fitRatioDown[i]->GetChisquare() / fitRatioDown[i]->GetNDF() << std::endl;

        // Calculate systematics
        grSyst[i] = new TGraph();
        jerSystematics(fitRatioUp[i], fitRatioDown[i], hEtaDef[i], grSyst[i]);
        grSyst[i]->SetName( Form("grSyst_%d", i) );

        grSystFB[i] = new TGraph();
        jeuFBSystematics(hEtaRatioFBUp[i], hEtaRatioFBDown[i], hEtaDefFBRatio[i], grSystFB[i]);
        grSystFB[i]->SetName( Form("grSystFB_%d", i) );

        // Write systematic uncertainty to ASCII file
        std::ofstream outFile( Form("%s/jer/%s_pPb8160_etaDijet_jerSyst_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFile.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSyst[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSyst[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFile << x << " " << y << std::endl;
            }
            outFile.close();
        }
        else {
            std::cerr << "[JerSyst] Unable to open file for writing" << std::endl;
        }

        std::ofstream outFileFB( Form("%s/jer/%s_pPb8160_etaDijet_jerSyst_FB_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFileFB.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSystFB[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSystFB[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFileFB << x << " " << y << std::endl;
            }
            outFileFB.close();
        }
        else {
            std::cerr << "[JerSyst] Unable to open file for writing" << std::endl;
        }

        //
        // Plot comparison of dijet pseudorapidity distributions
        //
        canv->cd();
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("JER def."), "p");
        leg->AddEntry(hEtaUp[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JER-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot BF comparison
        //
        canv->cd();
        setPadFBStyle();
        hEtaDefBFRatio[i]->Draw();
        hEtaUpBFRatio[i]->Draw("same");
        hEtaDownBFRatio[i]->Draw("same");
        hEtaDefBFRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefBFRatio[i]->GetYaxis()->SetRangeUser(0., 2.);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefBFRatio[i], Form("Backward/Forward"), "p");
        leg->AddEntry(hEtaUpBFRatio[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDownBFRatio[i], Form("JER-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot FB comparison
        //
        canv->cd();
        setPadFBStyle();
        hEtaDefFBRatio[i]->Draw();
        hEtaUpFBRatio[i]->Draw("same");
        hEtaDownFBRatio[i]->Draw("same");
        hEtaDefFBRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefFBRatio[i]->GetYaxis()->SetRangeUser(0., 2.);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.55, 0.65, 0.65, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefBFRatio[i], Form("Forward/Backward"), "p");
        leg->AddEntry(hEtaUpBFRatio[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDownBFRatio[i], Form("JER-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot ratios of pseudorapidity distributions to default
        //
        canv->cd();
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        if (drawFits) {
            fitRatioUp[i]->SetLineColor(kRed);
            fitRatioDown[i]->SetLineColor(kBlue);
            fitRatioUp[i]->Draw("same");
            fitRatioDown[i]->Draw("same");
        }
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("JER+ / Default"), "p");
        leg->AddEntry(hEtaRatioDown[i], Form("JER- / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot ratios of BF to default
        //
        canv->cd();
        setPadFBStyle();
        hEtaRatioBFUp[i]->Draw();
        hEtaRatioBFDown[i]->Draw("same");
        hEtaRatioBFUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioBFUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.26, 0.65, 0.46, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioBFUp[i], Form("JER+/Def"), "p");
        leg->AddEntry(hEtaRatioBFDown[i], Form("JER-/Def"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot ratios of FB to default
        //
        canv->cd();
        setPadFBStyle();
        hEtaRatioFBUp[i]->Draw();
        hEtaRatioFBDown[i]->Draw("same");
        hEtaRatioFBUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioFBUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.26, 0.65, 0.46, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioFBUp[i], Form("JER+/Def"), "p");
        leg->AddEntry(hEtaRatioFBDown[i], Form("JER-/Def"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot all dijet pseudorapidity comparisons on one canvas
        //
        cComp->cd(i+1);
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("JER def."), "p");
        leg->AddEntry(hEtaUp[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JER-"), "p");
        leg->Draw();

        //
        // Plot dijet pseudorapidity ratios on one canvas
        //
        cRat->cd(i+1);
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("JER+ / Default"), "p");
        leg->AddEntry(hEtaRatioDown[i], Form("JER- / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        if (drawFits) {
            fitRatioUp[i]->Draw("same");
            fitRatioDown[i]->Draw("same");
        }

        //
        // Plot comparisons of backward to forward ratios on one canvas
        //
        cBFComp->cd(i+1);
        setPadFBStyle();
        hEtaDefBFRatio[i]->Draw();
        hEtaUpBFRatio[i]->Draw("same");
        hEtaDownBFRatio[i]->Draw("same");
        hEtaDefBFRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaDefBFRatio[i]->GetYaxis()->SetRangeUser(0., 2.);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.26, 0.65, 0.46, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefBFRatio[i], Form("Backward/Forward"), "p");
        leg->AddEntry(hEtaUpBFRatio[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDownBFRatio[i], Form("JER-"), "p");
        leg->Draw();

        //
        // Plot BF ratios to default on one canvas
        //
        cBFRat->cd(i+1);
        setPadFBStyle();
        hEtaRatioBFUp[i]->Draw();
        hEtaRatioBFDown[i]->Draw("same");
        hEtaRatioBFUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioBFUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.26, 0.65, 0.46, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioBFUp[i], Form("JER+/Def"), "p");
        leg->AddEntry(hEtaRatioBFDown[i], Form("JER-/Def"), "p");
        leg->Draw();

        //
        // Plot comparisons of forward to backward ratios on one canvas
        //
        cFBComp->cd(i+1);
        setPadFBStyle();
        hEtaDefFBRatio[i]->Draw();
        hEtaUpFBRatio[i]->Draw("same");
        hEtaDownFBRatio[i]->Draw("same");
        hEtaDefFBRatio[i]->GetYaxis()->SetRangeUser(0., 3.);
        hEtaDefFBRatio[i]->GetYaxis()->SetRangeUser(0., 2.);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.55, 0.65, 0.65, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDefFBRatio[i], Form("Forward/Backward"), "p");
        leg->AddEntry(hEtaUpFBRatio[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDownFBRatio[i], Form("JER-"), "p");
        leg->Draw();

        //
        // Plot FB ratios to default on one canvas
        //
        cFBRat->cd(i+1);
        setPadFBStyle();
        hEtaRatioFBUp[i]->Draw();
        hEtaRatioFBDown[i]->Draw("same");
        hEtaRatioFBUp[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioFBUp[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        // t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.26, 0.65, 0.46, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioFBUp[i], Form("JER+/Def"), "p");
        leg->AddEntry(hEtaRatioFBDown[i], Form("JER-/Def"), "p");
        leg->Draw();
    } // for (Int_t i{0}; i<ptBins; i++)

    cComp->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cComp->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_jerRat_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cBFComp->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFComp->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cFBComp->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cFBComp->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_bf_jerRat_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/jer/%s_pPb8160_etaDijet_fb_jerRat_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
}

//________________
void pointingSystematics(TF1 *fitf, TH1D *h, TGraph* syst) {

    std::cout << Form("Def file name: %s\n", h->GetName());
    // Loop over the bins and estimate systematics for those that have entries
    for (Int_t i{1}; i<=h->GetNbinsX(); i++ ) {
        //if ( h->GetBinContent(i) == 0 ) continue;
        Double_t xVal = h->GetXaxis()->GetBinCenter(i);
        Double_t yVal = fitf->Eval( xVal );
        Double_t sign = ( (yVal - 1.) >= 0 ) ? 1 : -1;
        yVal = TMath::Abs( yVal - 1.);
        // if ( TMath::Abs( h->GetBinContent(i) ) < 0.0001 ) {
        //     yVal = 0.;
        // }
        std::cout << Form("bin content: %f x: %3.2f y [perc.]: %.3f ", h->GetBinContent(i), xVal, yVal * 100.) << std::endl;;
        syst->SetPoint(i-1, xVal, sign * yVal );
    } // for (Int_t i{1}; i<=h->GetNbinsX(); i++ )
}

//________________
void plotPointingResolution(TFile *embeddingFile, TString date, Bool_t isCM = kFALSE, 
                            Bool_t drawFits = kTRUE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";
    TString histoName = "hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted";
    if ( isCM ) {
        histoName = "hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted";
    }

    TString trigName = "MB";
    TString fName = embeddingFile->GetName();
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

    TH3D *hRecoEtaRefEtaRecoPt = (TH3D*)embeddingFile->Get( histoName.Data() );
    hRecoEtaRefEtaRecoPt->SetName("hRecoEtaRefEtaRecoPt");

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    // std::vector<Int_t> ptDijetBinLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    // std::vector<Int_t> ptDijetBinHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // Styles
    Int_t recoType{0};
    Int_t refType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaReco[ ptDijetBinLow.size() ];
    TH1D *hEtaRecoForward[ ptDijetBinLow.size() ];
    TH1D *hEtaRecoBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaRecoBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaRecoFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaRef[ ptDijetBinLow.size() ];
    TH1D *hEtaRefForward[ ptDijetBinLow.size() ];
    TH1D *hEtaRefBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaRefBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaRefFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBF[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFB[ ptDijetBinLow.size() ];
    TH2D *hEtaRefVsEtaReco[ ptDijetBinLow.size() ];

    TF1 *fitRatio[ ptDijetBinLow.size() ];

    TGraph *grSyst[ ptDijetBinLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetBinLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *c2D = new TCanvas("c2D", "c2D", sizeX, sizeY);
    c2D->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBComp = new TCanvas("cFBComp", "cFBComp", sizeX, sizeY);
    cFBComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFComp = new TCanvas("cBFComp", "cBFComp", sizeX, sizeY);
    cBFComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBRat = new TCanvas("cFBRat", "cFBRat", sizeX, sizeY);
    cFBRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFRat = new TCanvas("cBFRat", "cBFRat", sizeX, sizeY);
    cBFRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );


    TString fbTitle = "#eta^{dijet}; #frac{ dN/d#eta^{forward} }{ dN/d#eta^{backward} }";
    TString bfTitle = "#eta^{dijet}; #frac{ dN/d#eta^{backward} }{ dN/d#eta^{forward} }";
    if ( isCM ) {
        fbTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (forward) }{ dN/d#eta^{dijet}_{CM} (backward) }";
        bfTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (backward) }{ dN/d#eta^{dijet}_{CM} (forward) }";
    }

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // Project to 2D (ref vs reco)
        hRecoEtaRefEtaRecoPt->GetZaxis()->SetRange( ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        hEtaRefVsEtaReco[i] = dynamic_cast<TH2D*>( hRecoEtaRefEtaRecoPt->Project3D("yx") );
        hEtaRefVsEtaReco[i]->SetNameTitle( Form("hEtaRefVsEtaReco_%d", i),";#eta^{dijet}_{reco};#eta^{dijet}_{ref}" );
        rescaleEta( hEtaRefVsEtaReco[i] );
        set2DStyle( hEtaRefVsEtaReco[i] );

        // Project on eta reco
        hEtaReco[i] = dynamic_cast<TH1D*>(hRecoEtaRefEtaRecoPt->ProjectionX( Form("hEtaReco_%d", i), 1, -1, ptDijetBinLow.at(i), ptDijetBinHi.at(i) ) );
        rescaleEta( hEtaReco[i] );
        set1DStyle( hEtaReco[i], recoType );
        ( isCM ) ? hEtaReco[i]->GetXaxis()->SetTitle("#eta_{CM}") : hEtaReco[i]->GetXaxis()->SetTitle("#eta");
        hEtaRecoForward[i] = new TH1D( Form("hEtaRecoForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaRecoBackward[i] = new TH1D( Form("hEtaRecoBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaRecoBFRatio[i] = new TH1D( Form("hEtaRecoBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaRecoFBRatio[i] = new TH1D( Form("hEtaRecoFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaRecoForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaRecoForward[i]->Sumw2();  set1DStyle(hEtaRecoForward[i], recoType);
        hEtaRecoBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaRecoBackward[i]->Sumw2(); set1DStyle(hEtaRecoBackward[i], refType);
        hEtaRecoBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaRecoBFRatio[i]->Sumw2();  set1DStyle(hEtaRecoBFRatio[i], recoType);
        hEtaRecoFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaRecoFBRatio[i]->Sumw2();  set1DStyle(hEtaRecoFBRatio[i], recoType);
        makeDijetFB(hEtaReco[i], hEtaRecoBackward[i], hEtaRecoForward[i], hEtaRecoBFRatio[i], hEtaRecoFBRatio[i]);

        // Project on eta ref
        hEtaRef[i] = dynamic_cast<TH1D*>(hRecoEtaRefEtaRecoPt->ProjectionY( Form("hEtaRef_%d", i), 1, -1, ptDijetBinLow.at(i), ptDijetBinHi.at(i) ) );
        rescaleEta( hEtaRef[i] );
        set1DStyle( hEtaRef[i], refType );
        ( isCM ) ? hEtaRef[i]->GetXaxis()->SetTitle("#eta_{CM}") : hEtaRef[i]->GetXaxis()->SetTitle("#eta");
        hEtaRefForward[i] = new TH1D( Form("hEtaRefForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaRefBackward[i] = new TH1D( Form("hEtaRefBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaRefBFRatio[i] = new TH1D( Form("hEtaRefBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaRefFBRatio[i] = new TH1D( Form("hEtaRefFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaRefForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaRefForward[i]->Sumw2();  set1DStyle(hEtaRefForward[i], recoType);
        hEtaRefBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaRefBackward[i]->Sumw2(); set1DStyle(hEtaRefBackward[i], refType);
        hEtaRefBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaRefBFRatio[i]->Sumw2();  set1DStyle(hEtaRefBFRatio[i], refType);
        hEtaRefFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaRefFBRatio[i]->Sumw2();  set1DStyle(hEtaRefFBRatio[i], refType);
        makeDijetFB(hEtaRef[i], hEtaRefBackward[i], hEtaRefForward[i], hEtaRefBFRatio[i], hEtaRefFBRatio[i]);

        hEtaRatioBF[i] = dynamic_cast<TH1D*>( hEtaRecoBFRatio[i]->Clone( Form("hEtaRatioBF_%d", i) ) );
        hEtaRatioBF[i]->Divide( hEtaRatioBF[i], hEtaRefBFRatio[i], 1., 1.);
        set1DStyle(hEtaRatioBF[i], recoType);

        hEtaRatioFB[i] = dynamic_cast<TH1D*>( hEtaRecoFBRatio[i]->Clone( Form("hEtaRatioFB_%d", i) ) );
        hEtaRatioFB[i]->Divide( hEtaRatioBF[i], hEtaRefFBRatio[i], 1., 1.);
        set1DStyle(hEtaRatioFB[i], recoType);

        // Make ratios of reco to ref
        hEtaRatio[i] = dynamic_cast<TH1D*>( hEtaReco[i]->Clone( Form("hEtaRatio_%d", i) ) );
        make1DRatio(hEtaRatio[i], hEtaRef[i], Form("make1DRatio_%d", i), recoType);
        hEtaRatio[i]->GetYaxis()->SetTitle("reco / ref");

        // Create fit function
        fitRatio[i] = new TF1( Form("fitRatio_%d", i), "[0]+[1]*x+[2]*x*x", -4.8, 4.8 );
        fitRatio[i]->SetLineColor(kBlack);
        fitRatio[i]->SetLineWidth(2);
        fitRatio[i]->SetParameters(0., 0., 0.);
        fitRatio[i]->SetParLimits(2, -1e6, 0.);
        hEtaRatio[i]->Fit( Form("fitRatio_%d", i), "MRE0" );
        std::cout << "chi2/ndf down: " << fitRatio[i]->GetChisquare() / fitRatio[i]->GetNDF() << std::endl;

        // Create graph to store the data
        grSyst[i] = new TGraph();
        pointingSystematics( fitRatio[i], hEtaRatio[i], grSyst[i] );

        // Write systematic uncertainty to ASCII file
        std::ofstream outFile( Form("%s/pointing/%s_pPb8160_etaDijet_pointingSyst_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFile.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSyst[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSyst[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFile << x << " " << y << std::endl;
            }
            outFile.close();
        }
        else {
            std::cerr << "[PointingSyst] Unable to open file for writing" << std::endl;
        }

        //
        // Plot 2D matrix
        //
        canv->cd();
        setPadStyle();
        hEtaRefVsEtaReco[i]->Draw("colz");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefComp_2D_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefComp_2D_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot comparison of reco and ref dijet pseudorapidity distributions
        //
        canv->cd();
        setPadStyle();
        hEtaReco[i]->Draw();
        hEtaRef[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaReco[i], Form("Reco"), "p");
        leg->AddEntry(hEtaRef[i], Form("Ref"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot comparison of reco and ref BF distributions
        //
        canv->cd();
        setPadFBStyle();
        hEtaRecoBFRatio[i]->Draw();
        hEtaRefBFRatio[i]->Draw("same");
        hEtaRecoBFRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRecoBFRatio[i]->GetYaxis()->SetRangeUser(0.5, 4.5);
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRecoBFRatio[i], Form("Reco backward/Forward"), "p");
        leg->AddEntry(hEtaRefBFRatio[i], Form("Gen backward/Forward"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot ratio of reco and ref BF distributions
        //
        canv->cd();
        setPadFBStyle();
        hEtaRatioBF[i]->Draw();
        hEtaRatioBF[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioBF[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioBF[i]->GetYaxis()->SetTitle("Reco/Gen");
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioBF[i], Form("Reco/Gen"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );


        //
        // Plot comparison of reco and ref FB distributions
        //
        canv->cd();
        setPadFBStyle();
        hEtaRecoFBRatio[i]->Draw();
        hEtaRefFBRatio[i]->Draw("same");
        hEtaRecoFBRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRecoFBRatio[i]->GetYaxis()->SetRangeUser(0., 4.0);
        leg = new TLegend(0.7, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRecoFBRatio[i], Form("Reco forward/backward"), "p");
        leg->AddEntry(hEtaRefFBRatio[i], Form("Gen forward/backward"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefComp_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot ratio of reco and ref FB distributions
        //
        canv->cd();
        setPadFBStyle();
        hEtaRatioFB[i]->Draw();
        hEtaRatioFB[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioFB[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioFB[i]->GetYaxis()->SetTitle("Data/MC");
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioFB[i], Form("Reco/Gen"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefRat_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );


        //
        // Plot reco to ref ratio of dijet pseudorapidity distributions
        //
        canv->cd();
        setPadStyle();
        hEtaRatio[i]->Draw();
        if (drawFits) {
            fitRatio[i]->Draw("same");
        }
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefRatio_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefRatio_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot reco to ref ratio of dijet pseudorapidity distributions
        //
        canv->cd();
        setPadStyle();
        hEtaRatio[i]->Draw();
        if (drawFits) {
            fitRatio[i]->Draw("same");
        }
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefRatio_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        canv->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefRatio_pt_%d_%d_%s.png", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        //
        // Plot all 2D distributions dijet pseudorapidity distributions on one canvas
        //
        c2D->cd(i+1);
        setPadStyle();
        hEtaRefVsEtaReco[i]->Draw("colz");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );

        //
        // Plot all comparisons of dijet pseudorapidity distributions on one canvas
        //
        cComp->cd(i+1);
        setPadStyle();
        hEtaReco[i]->Draw();
        hEtaRef[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaReco[i], Form("Reco"), "p");
        leg->AddEntry(hEtaRef[i], Form("Gen"), "p");
        leg->Draw();

        //
        // Plot FB comparisons on one canvas
        //
        cFBComp->cd(i+1);
        setPadFBStyle();
        hEtaRecoFBRatio[i]->Draw();
        hEtaRefFBRatio[i]->Draw("same");
        hEtaRecoFBRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRecoFBRatio[i]->GetYaxis()->SetRangeUser(0., 1.1);
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRecoFBRatio[i], Form("Reco forward/backward"), "p");
        leg->AddEntry(hEtaRefFBRatio[i], Form("Gen forward/backward"), "p");
        leg->Draw();

        //
        // Plot BF comparisons on one canvas
        //
        cBFComp->cd(i+1);
        setPadFBStyle();  
        hEtaRecoBFRatio[i]->Draw();
        hEtaRefBFRatio[i]->Draw("same");
        hEtaRecoBFRatio[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRecoBFRatio[i]->GetYaxis()->SetRangeUser(0., 4.0);
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRecoBFRatio[i], Form("Reco backward/forward"), "p");
        leg->AddEntry(hEtaRefBFRatio[i], Form("Gen backward/forward"), "p");
        leg->Draw();

        //
        // Plot ratio of reco and ref BF distributions on one canvas
        //
        cBFRat->cd(i+1);
        setPadFBStyle();
        hEtaRatioBF[i]->Draw();
        hEtaRatioBF[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioBF[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioBF[i]->GetYaxis()->SetTitle("Reco/Gen");
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioBF[i], Form("Reco/Gen"), "p");
        leg->Draw();

        //
        // Plot ratio of reco and ref FB distributions
        //
        cFBRat->cd(i+1);
        setPadFBStyle();
        hEtaRatioFB[i]->Draw();
        hEtaRatioFB[i]->GetXaxis()->SetRangeUser(0., 3.);
        hEtaRatioFB[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hEtaRatioFB[i]->GetYaxis()->SetTitle("Data/MC");
        leg = new TLegend(0.65, 0.45, 0.85, 0.65);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioFB[i], Form("Data/MC"), "p");
        leg->Draw();

        //
        // Plot ratios on one canvas
        //
        cRat->cd(i+1);
        setPadStyle();
        hEtaRatio[i]->Draw();
        if (drawFits) {
            fitRatio[i]->Draw("same");
        }
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    } // for (Int_t i{0}; i<ptBins; i++) 

    cComp->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cComp->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    c2D->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefCorr_2D_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    c2D->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefCorr_2D_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefRatio_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_RecoRefRatio_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cBFComp->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFComp->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cFBComp->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cFBComp->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefComp_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );    
    cBFRat->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefRatio_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cBFRat->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_bf_RecoRefRatio_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
    cFBRat->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefRatio_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cFBRat->SaveAs( Form("%s/pointing/%s_pPb8160_etaDijet_fb_RecoRefRatio_all_%s.png", date.Data(), trigName.Data(), frame.Data() ) );
}

//________________
void pileupSystematics(TF1 *pileupUp, TF1 *pileupDown, TH1D *hDef, TGraph* syst) {

    // Loop over the bins and estimate systematics for those that have entries
    for (Int_t i{1}; i<=hDef->GetNbinsX(); i++ ) {
        // if ( hDef->GetBinContent(i) == 0 ) continue;
        Double_t xVal = hDef->GetBinCenter(i);

        Double_t yUp = pileupUp->Eval( xVal );
        Double_t yDown = pileupDown->Eval( xVal );
        std::cout << Form("x: %3.2f up: %.3f down: %.3f ", xVal, yUp, yDown);
        // Will use only gplus as a systematic variation
        Double_t sysYrel = yUp - 1.;
        // Combination of gplus and vtx1
        // Double_t sysYrel = ( TMath::Abs(yUp - 1.) + TMath::Abs(yDown - 1.) ) / 2;
        Double_t sysYabs = sysYrel * hDef->GetBinContent(i);

        std::cout << Form("bin content: %f syst [perc.]: %.4f syst [abs. val.]: %.8f\n", hDef->GetBinContent(i), sysYrel * 100., sysYabs );

        syst->SetPoint(i-1, xVal, sysYrel );
    } // for (Int_t i{1}; i<=h->GetNbinsX(); i++ )
    std::cout << "pileup Gplus chi2/ndf: " << pileupUp->GetChisquare() / pileupUp->GetNDF() 
              << " pileup Vtx1 chi2/ndf: " << pileupDown->GetChisquare() / pileupDown->GetNDF()
              << std::endl;
}

//________________
void plotPileup(TFile *defaultFile, TFile *gplusFile, TFile *vtx1File, TString date, 
                Bool_t isCM = kFALSE, Bool_t drawFits = kTRUE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";
    TString histoName = "hRecoDijetPtEtaDphiWeighted";
    if ( isCM ) {
        histoName = "hRecoDijetPtEtaDphiCMWeighted";
    }

    TString trigName = "MB";
    TString fName = defaultFile->GetName();
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

    TH3D *hPtEtaDphiDef = (TH3D*)defaultFile->Get( histoName.Data() );
    hPtEtaDphiDef->SetName("hPtEtaDphiDef");
    TH3D *hPtEtaDphiUp = (TH3D*)gplusFile->Get( histoName.Data() );
    hPtEtaDphiUp->SetName("hPtEtaDphiUp");
    TH3D *hPtEtaDphiDown = (TH3D*)vtx1File->Get( histoName.Data() );
    hPtEtaDphiDown->SetName("hPtEtaDphiDown");

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaDef[ ptDijetBinLow.size() ];
    TH1D *hEtaDefForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDefBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDefBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDefFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaUp[ ptDijetBinLow.size() ];
    TH1D *hEtaUpForward[ ptDijetBinLow.size() ];
    TH1D *hEtaUpBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaUpBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaUpFBRatio[ ptDijetBinLow.size() ];

    TH1D *hEtaDown[ ptDijetBinLow.size() ];
    TH1D *hEtaDownForward[ ptDijetBinLow.size() ];
    TH1D *hEtaDownBackward[ ptDijetBinLow.size() ];
    TH1D *hEtaDownBFRatio[ ptDijetBinLow.size() ];
    TH1D *hEtaDownFBRatio[ ptDijetBinLow.size() ];    

    TH1D *hEtaRatioUp[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFBUp[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBFUp[ ptDijetBinLow.size() ];

    TH1D *hEtaRatioDown[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioFBDown[ ptDijetBinLow.size() ];
    TH1D *hEtaRatioBFDown[ ptDijetBinLow.size() ];

    TF1 *fitRatioUp[ ptDijetBinLow.size() ];
    TF1 *fitRatioDown[ ptDijetBinLow.size() ];

    TGraph *grSyst[ ptDijetBinLow.size() ];
    TGraph *grSystFB[ ptDijetBinLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetBinLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBComp = new TCanvas("cFBComp", "cFBComp", sizeX, sizeY);
    cFBComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFComp = new TCanvas("cBFComp", "cBFComp", sizeX, sizeY);
    cBFComp->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cFBRat = new TCanvas("cFBRat", "cFBRat", sizeX, sizeY);
    cFBRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TCanvas *cBFRat = new TCanvas("cBFRat", "cBFRat", sizeX, sizeY);
    cBFRat->Divide(nPads, ( (ptBins % nPads) == 0 ) ? (ptBins / nPads) : (ptBins / nPads + 1) );

    TString fbTitle = "#eta^{dijet}; #frac{ dN/d#eta^{forward} }{ dN/d#eta^{backward} }";
    TString bfTitle = "#eta^{dijet}; #frac{ dN/d#eta^{backward} }{ dN/d#eta^{forward} }";
    if ( isCM ) {
        fbTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (forward) }{ dN/d#eta^{dijet}_{CM} (backward) }";
        bfTitle = "#eta^{dijet}_{CM}; #frac{ dN/d#eta^{dijet}_{CM} (backward) }{ dN/d#eta^{dijet}_{CM} (forward) }";
    }

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // Pileup default selection
        hEtaDef[i] = projectEtaFrom3D(hPtEtaDphiDef, Form("hEtaDef_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaDef[i] );
        set1DStyle(hEtaDef[i], defType);
        hEtaDefForward[i] = new TH1D( Form("hEtaDefForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDefBackward[i] = new TH1D( Form("hEtaDefBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDefBFRatio[i] = new TH1D( Form("hEtaDefBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDefFBRatio[i] = new TH1D( Form("hEtaDefFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDefForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefForward[i]->Sumw2();  set1DStyle(hEtaDefForward[i], upType);
        hEtaDefBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDefBackward[i]->Sumw2(); set1DStyle(hEtaDefBackward[i], downType);
        hEtaDefBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefBFRatio[i]->Sumw2();  set1DStyle(hEtaDefBFRatio[i], defType);
        hEtaDefFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDefFBRatio[i]->Sumw2();  set1DStyle(hEtaDefFBRatio[i], defType);
        makeDijetFB(hEtaDef[i], hEtaDefBackward[i], hEtaDefForward[i], hEtaDefBFRatio[i], hEtaDefFBRatio[i]);

        // Gplus (up) selection
        hEtaUp[i] = projectEtaFrom3D(hPtEtaDphiUp, Form("hEtaUp_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaUp[i] );
        set1DStyle(hEtaUp[i], upType);
        hEtaUpForward[i] = new TH1D( Form("hEtaUpForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaUpBackward[i] = new TH1D( Form("hEtaUpBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaUpBFRatio[i] = new TH1D( Form("hEtaUpBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaUpFBRatio[i] = new TH1D( Form("hEtaUpFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaUpForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpForward[i]->Sumw2();  set1DStyle(hEtaUpForward[i], upType);
        hEtaUpBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaUpBackward[i]->Sumw2(); set1DStyle(hEtaUpBackward[i], downType);
        hEtaUpBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpBFRatio[i]->Sumw2();  set1DStyle(hEtaUpBFRatio[i], upType);
        hEtaUpFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaUpFBRatio[i]->Sumw2();  set1DStyle(hEtaUpFBRatio[i], upType);
        makeDijetFB(hEtaUp[i], hEtaUpBackward[i], hEtaUpForward[i], hEtaUpBFRatio[i], hEtaUpFBRatio[i]);


        // Vtx1 (down) selection
        hEtaDown[i] = projectEtaFrom3D(hPtEtaDphiDown, Form("hEtaDown_%d", i), ptDijetBinLow.at(i), ptDijetBinHi.at(i) );
        rescaleEta( hEtaDown[i] );
        set1DStyle(hEtaDown[i], downType);
        hEtaDownForward[i] = new TH1D( Form("hEtaDownForward_%d",i), Form("#eta for forward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDownBackward[i] = new TH1D( Form("hEtaDownBackward_%d",i), Form("#eta for backward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDownBFRatio[i] = new TH1D( Form("hEtaDownBFRatio_%d",i), Form("#eta for backward/forward direction;%s", bfTitle.Data()), 15, 0., 5.);
        hEtaDownFBRatio[i] = new TH1D( Form("hEtaDownFBRatio_%d",i), Form("#eta for forward/backward direction;%s", fbTitle.Data()), 15, 0., 5.);
        hEtaDownForward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownForward[i]->Sumw2();  set1DStyle(hEtaDownForward[i], downType);
        hEtaDownBackward[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals); hEtaDownBackward[i]->Sumw2(); set1DStyle(hEtaDownBackward[i], downType);
        hEtaDownBFRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownBFRatio[i]->Sumw2();  set1DStyle(hEtaDownBFRatio[i], downType);
        hEtaDownFBRatio[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);  hEtaDownFBRatio[i]->Sumw2();  set1DStyle(hEtaDownFBRatio[i], downType);
        makeDijetFB(hEtaDown[i], hEtaDownBackward[i], hEtaDownForward[i], hEtaDownBFRatio[i], hEtaDownFBRatio[i]);

        // Ratio of pileup up to default ratio
        hEtaRatioUp[i] = dynamic_cast<TH1D*>( hEtaUp[i]->Clone( Form("hEtaRatioUp_%d", i) ) );
        make1DRatio(hEtaRatioUp[i], hEtaDef[i], Form("hEtaRatioUp_%d", i), upType );
        hEtaRatioUp[i]->GetYaxis()->SetTitle("vtx.flt. / dz1p0");
        //updateJERRatio( hEtaRatioUp[i] );

        // Ratio of pileup down to default ratio
        hEtaRatioDown[i] = dynamic_cast<TH1D*>( hEtaDown[i]->Clone( Form("hEtaRatioDown_%d", i) ) );
        make1DRatio(hEtaRatioDown[i], hEtaDef[i], Form("hEtaRatioDown_%d", i), downType );
        hEtaRatioDown[i]->GetYaxis()->SetTitle("vtx.flt. / dz1p0");
        //updateJERRatio( hEtaRatioDown[i] );

        // Make fit functions
        fitRatioUp[i] = new TF1(Form("fitRatioUp_%d",i), "[0]+[1]*x", -3.5, 3.5);
        // fitRatioUp[i] = new TF1(Form("fitRatioUp_%d",i), "[0]", -5., 5.);
        // fitRatioUp[i]->SetParameters(0., 0.0001);
        // fitRatioUp[i]->SetParLimits(2, 0.000001, 1e6);
        fitRatioUp[i]->SetParameters(0.);
        fitRatioUp[i]->SetLineColor(kRed);
        fitRatioUp[i]->SetLineWidth(2);

        fitRatioDown[i] = new TF1(Form("fitRatioDown_%d",i), "[0]+[1]*x", -3.5, 3.5);
        // fitRatioDown[i] = new TF1(Form("fitRatioDown_%d",i), "[0]", -5., 5.);
        // fitRatioDown[i]->SetParameters(0., 0.0001);
        // fitRatioDown[i]->SetParLimits(2, 0.000001, 1e6);
        fitRatioDown[i]->SetParameters(0.);
        fitRatioDown[i]->SetLineColor(kBlue);
        fitRatioDown[i]->SetLineWidth(2);

        // Perform fits
        hEtaRatioUp[i]->Fit(Form("fitRatioUp_%d",i), "MRE0");
        hEtaRatioDown[i]->Fit(Form("fitRatioDown_%d",i), "MRE0");
        std::cout << "chi2/ndf up: " << fitRatioUp[i]->GetChisquare() / fitRatioUp[i]->GetNDF() << std::endl;
        std::cout << "chi2/ndf down: " << fitRatioDown[i]->GetChisquare() / fitRatioDown[i]->GetNDF() << std::endl;

        // Calculate systematics
        grSyst[i] = new TGraph();
        pileupSystematics(fitRatioUp[i], fitRatioDown[i], hEtaDef[i], grSyst[i]);
        grSyst[i]->SetName( Form("grSyst_%d", i) );

        grSystFB[i] = new TGraph();
        jeuFBSystematics(hEtaUpFBRatio[i], nullptr, hEtaDefFBRatio[i], grSystFB[i]);
        grSystFB[i]->SetName( Form("grSystFB_%d", i) );

        // Write systematic uncertainty to ASCII file
        std::ofstream outFile( Form("%s/pileup/%s_pPb8160_etaDijet_pileupSyst_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFile.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSyst[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSyst[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFile << x << " " << y << std::endl;
            }
            outFile.close();
        }
        else {
            std::cerr << "[PileupSyst] Unable to open file for writing" << std::endl;
        }

        std::ofstream outFileFB( Form("%s/pileup/%s_pPb8160_etaDijet_pileupSystFB_%d_%d_%s.txt", date.Data(), trigName.Data(),
                                    ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                                    ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );
        if ( outFileFB.is_open() ) {
            for (Int_t iPoint = 0; iPoint < grSystFB[i]->GetN(); iPoint++) {
                // std::cout << "i: " << iPoint;
                Double_t x{0}, y{0};
                grSystFB[i]->GetPoint(iPoint, x, y);
                // std::cout << " x: " << x << " y: " << y << std::endl;
                outFileFB << x << " " << y << std::endl;
            }
            outFileFB.close();
        }
        else {
            std::cerr << "[PileupSyst] Unable to open file for writing" << std::endl;
        }

        // Plot comparison
        canv->cd();
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        //hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("dz1p0 def."), "p");
        leg->AddEntry(hEtaUp[i], Form("Gplus"), "p");
        //leg->AddEntry(hEtaDown[i], Form("Vtx1"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pileup/%s_pPb8160_etaDijet_pileupComp_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios
        canv->cd();
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        //hEtaRatioDown[i]->Draw("same");
        if (drawFits) {
            fitRatioUp[i]->SetLineColor(kBlack);
            fitRatioUp[i]->Draw("same");
            //fitRatioDown[i]->Draw("same");
        }
        hEtaRatioUp[i]->GetYaxis()->SetRangeUser(0.95, 1.05);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("Gplus / dz1p0"), "p");
        //leg->AddEntry(hEtaRatioDown[i], Form("Vtx1 / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pileup/%s_pPb8160_etaDijet_pileupRat_pt_%d_%d_%s.pdf", date.Data(), trigName.Data(),
                      ptLow + (ptDijetBinLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetBinHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        //hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("dz1p0"), "p");
        leg->AddEntry(hEtaUp[i], Form("Gplus"), "p");
        //leg->AddEntry(hEtaDown[i], Form("Vtx1"), "p");
        leg->Draw();

        // Plot all ratios on one canvas
        cRat->cd(i+1);
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        //hEtaRatioDown[i]->Draw("same");
        hEtaRatioUp[i]->GetYaxis()->SetRangeUser(0.95, 1.05);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptLow + (ptDijetBinLow.at(i) - 1) * ptStep, ptLow + ptDijetBinHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("Gplus / dz1p0"), "p");
        //leg->AddEntry(hEtaRatioDown[i], Form("Vtx1 / dz1p0"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        if (drawFits) {
            fitRatioUp[i]->Draw("same");
            //fitRatioDown[i]->Draw("same");
        }
    } // for (Int_t i{0}; i<ptBins; i++)

    cComp->SaveAs( Form("%s/pileup/%s_pPb8160_etaDijet_pileupComp_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pileup/%s_pPb8160_etaDijet_pileupRat_all_%s.pdf", date.Data(), trigName.Data(), frame.Data() ) );

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
void plotFinalEtaDistributions(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist, 
                               std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist,
                               std::vector< std::vector<TH1D*> > hMBTotalRelSystDist, std::vector< std::vector<TH1D*> > hJet60TotalRelSystDist,
                               std::vector< std::vector<TH1D*> > hJet80TotalRelSystDist, std::vector< std::vector<TH1D*> > hJet100TotalRelSystDist,
                               TString date, Bool_t isCM = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";

    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0.001, 0.12};

    Double_t xRangeForward[2] = {0., 3.};
    Double_t yRangeForward[2] = {0.001, 0.2};

    Double_t xRangeFB[2] = {0., 2.4};
    Double_t yRangeFB[2] = {0.75, 1.5};
    Double_t yRangeBF[2] = {0.4, 1.15};

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

    // Systematics (relative systematic uncertainty)
    std::vector< TH1D* > hMBEtaTotalRelSyst = hMBTotalRelSystDist.at(0);
    std::vector< TH1D* > hMBEtaForwardTotalRelSyst = hMBTotalRelSystDist.at(1);
    std::vector< TH1D* > hMBEtaBackwardTotalRelSyst = hMBTotalRelSystDist.at(2);
    std::vector< TH1D* > hMBEtaFBRatioTotalRelSyst = hMBTotalRelSystDist.at(3);
    std::vector< TH1D* > hMBEtaBFRatioTotalRelSyst = hMBTotalRelSystDist.at(4);

    std::vector< TH1D* > hJet60EtaTotalRelSyst = hJet60TotalRelSystDist.at(0);
    std::vector< TH1D* > hJet60EtaForwardTotalRelSyst = hJet60TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet60EtaBackwardTotalRelSyst = hJet60TotalRelSystDist.at(2);
    std::vector< TH1D* > hJet60EtaFBRatioTotalRelSyst = hJet60TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatioTotalRelSyst = hJet60TotalRelSystDist.at(4);

    std::vector< TH1D* > hJet80EtaTotalRelSyst = hJet80TotalRelSystDist.at(0);
    std::vector< TH1D* > hJet80EtaForwardTotalRelSyst = hJet80TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet80EtaBackwardTotalRelSyst = hJet80TotalRelSystDist.at(2);
    std::vector< TH1D* > hJet80EtaFBRatioTotalRelSyst = hJet80TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatioTotalRelSyst = hJet80TotalRelSystDist.at(4);

    std::vector< TH1D* > hJet100EtaTotalRelSyst = hJet100TotalRelSystDist.at(0);
    std::vector< TH1D* > hJet100EtaForwardTotalRelSyst = hJet100TotalRelSystDist.at(1);
    std::vector< TH1D* > hJet100EtaBackwardTotalRelSyst = hJet100TotalRelSystDist.at(2);
    std::vector< TH1D* > hJet100EtaFBRatioTotalRelSyst = hJet100TotalRelSystDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatioTotalRelSyst = hJet100TotalRelSystDist.at(4);

    // Create TGraphErrors with total absolute systematic uncertainty
    std::vector< TGraphErrors* > grMBEtaTotalAbsSyst = {};
    std::vector< TGraphErrors* > grMBEtaForwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grMBEtaBackwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grMBEtaFBRatioTotalAbsSyst = {};
    std::vector< TGraphErrors* > grMBEtaBFRatioTotalAbsSyst = {};

    std::vector< TGraphErrors* > grJet60EtaTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet60EtaForwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet60EtaBackwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet60EtaFBRatioTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet60EtaBFRatioTotalAbsSyst = {};

    std::vector< TGraphErrors* > grJet80EtaTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet80EtaForwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet80EtaBackwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet80EtaFBRatioTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet80EtaBFRatioTotalAbsSyst = {};

    std::vector< TGraphErrors* > grJet100EtaTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet100EtaForwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet100EtaBackwardTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet100EtaFBRatioTotalAbsSyst = {};
    std::vector< TGraphErrors* > grJet100EtaBFRatioTotalAbsSyst = {};

    // Loop over pT bins
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

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.07);

    TLegend *leg;

    Double_t drawPt[2] = {0.35, 0.83};

    // Plot distributions on canvases
    for (Int_t i=0; i<ptBins; i++) {


        grMBEtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hMBEta[i], hMBEtaTotalRelSyst[i], "MBEta", i) );
        grMBEtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hMBEtaForward[i], hMBEtaForwardTotalRelSyst[i], "MBEtaForward", i) );
        grMBEtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hMBEtaBackward[i], hMBEtaBackwardTotalRelSyst[i], "MBEtaBackward", i) );
        grMBEtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hMBEtaFBRatio[i], hMBEtaFBRatioTotalRelSyst[i], "MBEtaFBRatio", i) );
        grMBEtaBFRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hMBEtaBFRatio[i], hMBEtaBFRatioTotalRelSyst[i], "MBEtaBFRatio", i) );

        grJet60EtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60Eta[i], hJet60EtaTotalRelSyst[i], "Jet60Eta", i) );
        grJet60EtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaForward[i], hJet60EtaForwardTotalRelSyst[i], "Jet60EtaForward", i) );
        grJet60EtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaBackward[i], hJet60EtaBackwardTotalRelSyst[i], "Jet60EtaBackward", i) );
        grJet60EtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaFBRatio[i], hJet60EtaFBRatioTotalRelSyst[i], "Jet60EtaFBRatio", i) );

        grJet80EtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80Eta[i], hJet80EtaTotalRelSyst[i], "Jet80Eta", i) );
        grJet80EtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaForward[i], hJet80EtaForwardTotalRelSyst[i], "Jet80EtaForward", i) );
        grJet80EtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaBackward[i], hJet80EtaBackwardTotalRelSyst[i], "Jet80EtaBackward", i) );
        grJet80EtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaFBRatio[i], hJet80EtaFBRatioTotalRelSyst[i], "Jet80EtaFBRatio", i) );
        grJet80EtaBFRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaBFRatio[i], hJet80EtaBFRatioTotalRelSyst[i], "Jet80EtaBFRatio", i) );

        grJet100EtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100Eta[i], hJet100EtaTotalRelSyst[i], "Jet100Eta", i) );
        grJet100EtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaForward[i], hJet100EtaForwardTotalRelSyst[i], "Jet100EtaForward", i) );
        grJet100EtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaBackward[i], hJet100EtaBackwardTotalRelSyst[i], "Jet100EtaBackward", i) );
        grJet100EtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaFBRatio[i], hJet100EtaFBRatioTotalRelSyst[i], "Jet100EtaFBRatio", i) );
        grJet100EtaBFRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaBFRatio[i], hJet100EtaBFRatioTotalRelSyst[i], "Jet100EtaBFRatio", i) );

        if (ptDijetHi.at(i) < 80) {



            canv->cd();
            setPadStyle();
            setGraphStyle(grMBEtaTotalAbsSyst[i], isCM);
            grMBEtaTotalAbsSyst[i]->Draw("AE5");
            hMBEta[i]->Draw("same");
            grMBEtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grMBEtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
            t.DrawLatexNDC(0.2, 0.66, "p_{T}^{Subleading} > 40 GeV");
            if ( isCM ) {
                t.DrawLatexNDC(0.2, 0.57, "|#eta_{CM}| < 2.4");
            }
            else {
                t.DrawLatexNDC(0.2, 0.57, "|#eta| < 3");
            }
            t.DrawLatexNDC(0.2, 0.49, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
            leg = new TLegend(0.6, 0.65, 0.8, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hMBEta[i], "Data", "p");
            leg->AddEntry(grMBEtaTotalAbsSyst[i], "Syst. uncrt.", "f");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/MB_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            setGraphStyle(grMBEtaTotalAbsSyst[i], isCM);
            grMBEtaTotalAbsSyst[i]->Draw("AE5");
            hMBEta[i]->Draw("same");
            grMBEtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grMBEtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            if (i == 0) {
                t.DrawLatexNDC(0.2, 0.75, "p_{T}^{Leading} > 50 GeV");
                t.DrawLatexNDC(0.2, 0.65, "p_{T}^{Subleading} > 40 GeV");
                if ( isCM ) {
                    t.DrawLatexNDC(0.2, 0.57, "|#eta_{CM}| < 2.4");
                }
                else {
                    t.DrawLatexNDC(0.2, 0.57, "|#eta| < 3");
                }
                t.DrawLatexNDC(0.2, 0.49, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
            }
            if ( i == 1 ) {
                leg = new TLegend(0.6, 0.65, 0.8, 0.8);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEta[i], "Data", "p");
                leg->AddEntry(grMBEtaTotalAbsSyst[i], "Syst. uncrt.", "f");
                leg->Draw();
            }
            plotCMSHeader();

            canv->cd();
            setPadStyle();
            hMBEtaForward[i]->Draw();
            hMBEtaBackward[i]->Draw("same");
            hMBEtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hMBEtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            leg = new TLegend(0.6, 0.65, 0.8, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hMBEtaForward[i], "Forward", "p");
            leg->AddEntry(hMBEtaBackward[i], "Backward", "p");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/MB_pPb8160_etaDijet_forward_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaForwardAndBackward->cd(i+1);
            setPadStyle();
            hMBEtaForward[i]->Draw();
            hMBEtaBackward[i]->Draw("same");
            hMBEtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hMBEtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            if (i==0) {
                leg = new TLegend(0.6, 0.65, 0.8, 0.8);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEtaForward[i], "Forward", "p");
                leg->AddEntry(hMBEtaBackward[i], "Backward", "p");
                leg->Draw();
            }
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            setGraphStyle(grMBEtaFBRatioTotalAbsSyst[i], isCM);
            grMBEtaFBRatioTotalAbsSyst[i]->Draw("AE5");
            hMBEtaFBRatio[i]->Draw("same");
            grMBEtaFBRatioTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            grMBEtaFBRatioTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            t.DrawLatexNDC(0.3, 0.75, "p_{T}^{Leading} > 50 GeV");
            t.DrawLatexNDC(0.3, 0.65, "p_{T}^{Subleading} > 40 GeV");
            if ( isCM ) {
                t.DrawLatexNDC(0.3, 0.57, "|#eta_{CM}| < 2.4");
            }
            else {
                t.DrawLatexNDC(0.3, 0.57, "|#eta| < 3");
            }
            t.DrawLatexNDC(0.3, 0.49, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
            leg = new TLegend(0.3, 0.2, 0.5, 0.4);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hMBEtaFBRatio[i], "Data", "p");
            leg->AddEntry(grMBEtaFBRatioTotalAbsSyst[i], "Syst. uncrt.", "f");
            leg->Draw();
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/MB_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            setGraphStyle(grMBEtaFBRatioTotalAbsSyst[i], isCM);
            grMBEtaFBRatioTotalAbsSyst[i]->Draw("AE5");
            hMBEtaFBRatio[i]->Draw("same");
            grMBEtaFBRatioTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            grMBEtaFBRatioTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            if (i == 0) {
                t.DrawLatexNDC(0.3, 0.75, "p_{T}^{Leading} > 50 GeV");
                t.DrawLatexNDC(0.3, 0.65, "p_{T}^{Subleading} > 40 GeV");
                if ( isCM ) {
                    t.DrawLatexNDC(0.3, 0.57, "|#eta_{CM}| < 2.4");
                }
                else {
                    t.DrawLatexNDC(0.3, 0.57, "|#eta| < 3");
                }
                t.DrawLatexNDC(0.3, 0.49, "#Delta#phi^{dijet} > #frac{5#pi}{6}");
            }
            if ( i == 1 ) {
                leg = new TLegend(0.3, 0.2, 0.5, 0.4);
                leg->SetTextSize(0.06);
                leg->SetLineWidth(0);
                leg->AddEntry(hMBEtaFBRatio[i], "Data", "p");
                leg->AddEntry(grMBEtaFBRatioTotalAbsSyst[i], "Syst. uncrt.", "f");
                leg->Draw();
            }
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hMBEtaBFRatio[i]->Draw();
            hMBEtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hMBEtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/MB_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hMBEtaBFRatio[i]->Draw();
            hMBEtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hMBEtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
        }
        else if (ptDijetHi.at(i) < 100) {

            std::cout << "Ahhhtung 0" << std::endl;
            grJet60EtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60Eta[i], hJet60EtaTotalRelSyst[i], "Jet60Eta", i) );
            grJet60EtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaForward[i], hJet60EtaForwardTotalRelSyst[i], "Jet60EtaForward", i) );
            grJet60EtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaBackward[i], hJet60EtaBackwardTotalRelSyst[i], "Jet60EtaBackward", i) );
            grJet60EtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaFBRatio[i], hJet60EtaFBRatioTotalRelSyst[i], "Jet60EtaFBRatio", i) );
            grJet60EtaBFRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet60EtaBFRatio[i], hJet60EtaBFRatioTotalRelSyst[i], "Jet60EtaBFRatio", i) );

            canv->cd();std::cout << "Ahhhtung  1" << std::endl;

            setPadStyle();
            std::cout << Form("Name: %s", grJet60EtaTotalAbsSyst[i]->GetName()) << std::endl;
            setGraphStyle(grJet60EtaTotalAbsSyst[i], isCM);

            std::cout << "Ahhhtung" << std::endl;
            grJet60EtaTotalAbsSyst[i]->Draw("AE5");
            hJet60Eta[i]->Draw("same");
            grJet60EtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grJet60EtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet60_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            setGraphStyle(grJet60EtaTotalAbsSyst[i], isCM);
            grJet60EtaTotalAbsSyst[i]->Draw("AE5");
            hJet60Eta[i]->Draw("same");
            grJet60EtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grJet60EtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadStyle();
            hJet60EtaForward[i]->Draw();
            hJet60EtaForward[i]->Draw("same");
            hJet60EtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hJet60EtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            leg = new TLegend(0.6, 0.65, 0.8, 0.8);
            leg->SetTextSize(0.06);
            leg->SetLineWidth(0);
            leg->AddEntry(hJet60EtaForward[i], "Forward", "p");
            leg->AddEntry(hJet60EtaBackward[i], "Backward", "p");
            leg->Draw();
            canv->SaveAs( Form("%s/data/Jet60_pPb8160_etaDijet_forward_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaForwardAndBackward->cd(i+1);
            setPadStyle();
            hJet60EtaForward[i]->Draw();
            hJet60EtaForward[i]->Draw("same");
            hJet60EtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hJet60EtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            setGraphStyle(grJet60EtaFBRatioTotalAbsSyst[i], isCM);
            grJet60EtaFBRatioTotalAbsSyst[i]->Draw("AE5");
            hJet60EtaFBRatio[i]->Draw("same");
            grJet60EtaFBRatioTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            grJet60EtaFBRatioTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet60_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            setGraphStyle(grJet60EtaFBRatioTotalAbsSyst[i], isCM);
            grJet60EtaFBRatioTotalAbsSyst[i]->Draw("AE5");
            hJet60EtaFBRatio[i]->Draw("same");
            grJet60EtaFBRatioTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            grJet60EtaFBRatioTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet60EtaBFRatio[i]->Draw();
            hJet60EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet60EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet60_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hJet60EtaBFRatio[i]->Draw();
            hJet60EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet60EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
        }
        else if (ptDijetHi.at(i) < 120) {

            grJet80EtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80Eta[i], hJet80EtaTotalRelSyst[i], "Jet80Eta", i) );
            grJet80EtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaForward[i], hJet80EtaForwardTotalRelSyst[i], "Jet80EtaForward", i) );
            grJet80EtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaBackward[i], hJet80EtaBackwardTotalRelSyst[i], "Jet80EtaBackward", i) );
            grJet80EtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaFBRatio[i], hJet80EtaFBRatioTotalRelSyst[i], "Jet80EtaFBRatio", i) );
            grJet80EtaBFRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet80EtaBFRatio[i], hJet80EtaBFRatioTotalRelSyst[i], "Jet80EtaBFRatio", i) );

            canv->cd();
            setPadStyle();
            setGraphStyle(grJet80EtaTotalAbsSyst[i], isCM);
            grJet80EtaTotalAbsSyst[i]->Draw("AE5");
            hJet80Eta[i]->Draw("same");
            grJet80EtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grJet80EtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet80_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            setGraphStyle(grJet80EtaTotalAbsSyst[i], isCM);
            grJet80EtaTotalAbsSyst[i]->Draw("AE5");
            hJet80Eta[i]->Draw("same");
            grJet80EtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grJet80EtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadStyle();
            hJet80EtaForward[i]->Draw();
            hJet80EtaBackward[i]->Draw("same");
            hJet80EtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hJet80EtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet80_pPb8160_etaDijet_forward_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );
            
            cEtaForwardAndBackward->cd(i+1);
            setPadStyle();
            hJet80EtaForward[i]->Draw();
            hJet80EtaBackward[i]->Draw("same");
            hJet80EtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hJet80EtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet80EtaFBRatio[i]->Draw();
            hJet80EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet80_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            hJet80EtaFBRatio[i]->Draw();
            hJet80EtaFBRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaFBRatio[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet80EtaBFRatio[i]->Draw();
            hJet80EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet80_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hJet80EtaBFRatio[i]->Draw();
            hJet80EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet80EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
        }
        else {

            grJet100EtaTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100Eta[i], hJet100EtaTotalRelSyst[i], "Jet100Eta", i) );
            grJet100EtaForwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaForward[i], hJet100EtaForwardTotalRelSyst[i], "Jet100EtaForward", i) );
            grJet100EtaBackwardTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaBackward[i], hJet100EtaBackwardTotalRelSyst[i], "Jet100EtaBackward", i) );
            grJet100EtaFBRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaFBRatio[i], hJet100EtaFBRatioTotalRelSyst[i], "Jet100EtaFBRatio", i) );
            grJet100EtaBFRatioTotalAbsSyst.push_back( retrieveTotalAbsSystUncrt(hJet100EtaBFRatio[i], hJet100EtaBFRatioTotalRelSyst[i], "Jet100EtaBFRatio", i) );

            canv->cd();
            setPadStyle();
            setGraphStyle(grJet100EtaTotalAbsSyst[i], isCM);
            grJet100EtaTotalAbsSyst[i]->Draw("AE5");
            hJet100Eta[i]->Draw("same");
            grJet100EtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grJet100EtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet100_pPb8160_etaDijet_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEta->cd(i+1);
            setPadStyle();
            setGraphStyle(grJet100EtaTotalAbsSyst[i], isCM);
            grJet100EtaTotalAbsSyst[i]->Draw("AE5");
            hJet100Eta[i]->Draw("same");
            grJet100EtaTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
            grJet100EtaTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadStyle();
            hJet100EtaForward[i]->Draw();
            hJet100EtaBackward[i]->Draw("same");
            hJet100EtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hJet100EtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet100_pPb8160_etaDijet_forward_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaForwardAndBackward->cd(i+1);
            setPadStyle();
            hJet100EtaForward[i]->Draw();
            hJet100EtaBackward[i]->Draw("same");
            hJet100EtaForward[i]->GetXaxis()->SetRangeUser(xRangeForward[0], xRangeForward[1]);
            hJet100EtaForward[i]->GetYaxis()->SetRangeUser(yRangeForward[0], yRangeForward[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            setGraphStyle(grJet100EtaFBRatioTotalAbsSyst[i], isCM);
            grJet100EtaFBRatioTotalAbsSyst[i]->Draw("AE5");
            hJet100EtaFBRatio[i]->Draw("same");
            grJet100EtaFBRatioTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            grJet100EtaFBRatioTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet100_pPb8160_etaDijet_fb_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaFB->cd(i+1);
            setPadFBStyle();
            setGraphStyle(grJet100EtaFBRatioTotalAbsSyst[i], isCM);
            grJet100EtaFBRatioTotalAbsSyst[i]->Draw("AE5");
            hJet100EtaFBRatio[i]->Draw("same");
            grJet100EtaFBRatioTotalAbsSyst[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            grJet100EtaFBRatioTotalAbsSyst[i]->GetYaxis()->SetRangeUser(yRangeFB[0], yRangeFB[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();

            canv->cd();
            setPadFBStyle();
            hJet100EtaBFRatio[i]->Draw();
            hJet100EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet100EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
            canv->SaveAs( Form("%s/data/Jet100_pPb8160_etaDijet_bf_pt_%d_%d_%s.pdf", date.Data(), ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

            cEtaBF->cd(i+1);
            setPadFBStyle();
            hJet100EtaBFRatio[i]->Draw();
            hJet100EtaBFRatio[i]->GetXaxis()->SetRangeUser(xRangeFB[0], xRangeFB[1]);
            hJet100EtaBFRatio[i]->GetYaxis()->SetRangeUser(yRangeBF[0], yRangeBF[1]);
            t.DrawLatexNDC(drawPt[0], drawPt[1], Form("%d < p_{T}^{ave} (GeV) < %d", ptDijetLow.at(i), ptDijetHi.at(i) ) );
            plotCMSHeader();
        }

    } // for (Int_t i=0; i<ptBins; i++)
    cEta->SaveAs( Form("%s/data/pPb8160_etaDijet_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaForwardAndBackward->SaveAs( Form("%s/data/pPb8160_etaDijet_forward_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaFB->SaveAs( Form("%s/data/pPb8160_etaDijet_fb_all_%s.pdf", date.Data(), frame.Data() ) );
    cEtaBF->SaveAs( Form("%s/data/pPb8160_etaDijet_bf_all_%s.pdf", date.Data(), frame.Data() ) );
}

//________________
void plotUp2DownComparison(std::vector<TH1D*> hDef, std::vector<TH1D*> hUp, std::vector<TH1D*> hDown,
                           std::vector<TH1D*> *hRatioUp, std::vector<TH1D*> *hRatioDown,
                           Int_t systematics, Int_t fitType, TString date, TString frame) {

    // Dijet pT selection
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);
    TLegend *leg;
    TLine *line;

    Bool_t drawFits{kFALSE};
    // 0 - symmetric, 1 - forward, backward, fb, bf
    Int_t etaType{0}; 
    // Ranges for plotting
    Double_t xMin{-3.}, xMax{3.};
    Double_t yMin{0.001}, yMax{0.15};
    Double_t legXmin{0.27}, legXmax{0.42}, legYmin{0.65}, legYmax{0.85};

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
    else {
        drawFits = kFALSE;
        systType = "data2mc";
    }

    TString hDefName = hDef.at(0)->GetName();
    TString dirName;
    if ( hDefName.Contains("Forward") || hDefName.Contains("Backward") || hDefName.Contains("FB") || hDefName.Contains("BF") ) {
        etaType = 1;
        xMin = 0.;
        drawFits = kFALSE;
        if ( hDefName.Contains("Forward") ) {
            dirName = "forward_";
        }
        if ( hDefName.Contains("Backward") ) {
            dirName = "backward_";
        }
        if ( hDefName.Contains("FBRatio") ) {
            dirName = "fb_";
            yMin = 0.001; yMax = 1.2;
            legYmin = 0.3; legYmax = 0.5;
        }
        if ( hDefName.Contains("BFRatio") ) {
            dirName = "bf_";
            yMin = 0.9; yMax = 5.1;
            // legXmin = 0.6; legXmax = 0.8; legYmin = 0.3; legYmax = 0.5;
        }
    }

    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    // Fits (Fit type: 0 -constant, 1 - linear, 2 - quadratic, 3 - cubic, 4 - quartic)
    std::vector<TF1*> fitUp( ptDijetBinLow.size() );
    std::vector<TF1*> fitDown( ptDijetBinLow.size() );

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

    // Loop over pT average bins
    for (Int_t i{0}; i<ptDijetBinLow.size(); i++) {

        set1DStyle(hDef.at(i), defType);
        set1DStyle(hUp.at(i), upType);
        if ( !hDown.empty() ) {
            set1DStyle(hDown.at(i), downType);
        }

        // Make ratio histograms
        hRatioUp->push_back( dynamic_cast<TH1D*>( hUp.at(i)->Clone( Form("%s_%s_ratioUp", hUp.at(i)->GetName(), systType.Data() ) ) ) );
        make1DRatio(hRatioUp->at(i), hDef.at(i), Form("%s_%s_ratioUp", hUp.at(i)->GetName(), systType.Data() ) );
        hRatioUp->at(i)->GetYaxis()->SetTitle( Form("%s variation / def", systType.Data() ) );
        set1DStyle(hRatioUp->at(i), upType);

        if ( !hDown.empty() ) {
            hRatioDown->push_back( dynamic_cast<TH1D*>( hDown.at(i)->Clone( Form("%s_%s_ratioDown", hDown.at(i)->GetName(), systType.Data() ) ) ) );
            make1DRatio(hRatioDown->at(i), hDef.at(i), Form("%s_%s_ratioDown", hDown.at(i)->GetName(), systType.Data() ) );
            hRatioDown->at(i)->GetYaxis()->SetTitle( Form("%s variation / def", systType.Data() ) );
            set1DStyle(hRatioDown->at(i), downType);
        }

        if ( fitType == 0 ) {
            fitUp.at(i) = new TF1(Form("fitUp_%d",i), "[0]", xMin, xMax);
            fitDown.at(i) = new TF1(Form("fitDown_%d",i), "[0]", xMin, xMax);
        }
        else if ( fitType == 1 ) {
            fitUp.at(i) = new TF1(Form("fitUp_%d",i), "[0]+[1]*x", xMin, xMax);
            fitDown.at(i) = new TF1(Form("fitDown_%d",i), "[0]+[1]*x", xMin, xMax);
        }
        else if ( fitType == 2 ) {
            fitUp.at(i) = new TF1(Form("fitUp_%d",i), "[0]+[1]*x+[2]*x*x", xMin, xMax);
            fitDown.at(i) = new TF1(Form("fitDown_%d",i), "[0]+[1]*x+[2]*x*x", xMin, xMax);
        }
        else if ( fitType == 3 ) {
            fitUp.at(i) = new TF1(Form("fitUp_%d",i), "[0]+[1]*x+[2]*x*x+[3]*x*x*x", xMin, xMax);
            fitDown.at(i) = new TF1(Form("fitDown_%d",i), "[0]+[1]*x+[2]*x*x+[3]*x*x*x", xMin, xMax);
        }
        else if ( fitType == 4 ) {
            fitUp.at(i) = new TF1(Form("fitUp_%d",i), "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", xMin, xMax);
            fitDown.at(i) = new TF1(Form("fitDown_%d",i), "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", xMin, xMax);
        }
        else {
            fitUp.at(i) = new TF1(Form("fitUp_%d",i), "[0]+[1]*x", xMin, xMax);
            fitDown.at(i) = new TF1(Form("fitDown_%d",i), "[0]+[1]*x", xMin, xMax);
        }

        fitUp.at(i)->SetLineColor(kRed); fitUp.at(i)->SetLineWidth(2);
        fitDown.at(i)->SetLineColor(kBlue); fitDown.at(i)->SetLineWidth(2);

        // Perform fits
        if ( drawFits ) {
            hRatioUp->at(i)->Fit(Form("fitUp_%d",i), "MRE0");
            if ( !hDown.empty() ) {
                hRatioDown->at(i)->Fit(Form("fitDown_%d",i), "MRE0");
            }
        }

        // std::cout << "chi2/ndf up: " << fitUp.at(i)->GetChisquare() / fitUp.at(i)->GetNDF() << std::endl;
        // std::cout << "chi2/ndf down: " << fitDown.at(i)->GetChisquare() / fitDown.at(i)->GetNDF() << std::endl;

        // Plot comparison
        canv->cd();
        setPadFBStyle();
        hDef.at(i)->Draw();
        hUp.at(i)->Draw("same");
        if ( !hDown.empty() ) {
            hDown.at(i)->Draw("same");
        }
        hDef.at(i)->GetYaxis()->SetRangeUser(yMin, yMax);
        hDef.at(i)->GetXaxis()->SetRangeUser(xMin, xMax);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetLow.at(i), ptDijetHi.at(i) ) );
        //t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );

        leg = new TLegend(legXmin, legYmin, legXmax, legYmax);

        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hDef.at(i), "Default", "p");
        leg->AddEntry(hUp.at(i), "Up", "p");
        if ( !hDown.empty() ) {
            leg->AddEntry(hDown.at(i), "Down", "p");
        }
        leg->Draw();
        // canv->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sComp_%spt_%d_%d_%s.pdf", date.Data(), systType.Data(), 
        //                   returnTrigName( hDef.at(i)->GetName() ).Data(), systType.Data(), dirName.Data(),
        //                   ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );
        canv->Print( Form("%s/%s/%s_pPb8160_etaDijet_%sComp_%spt_%d_%d_%s.pdf", date.Data(), systType.Data(), 
                          returnTrigName( hDef.at(i)->GetName() ).Data(), systType.Data(), dirName.Data(),
                          ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

        // Plot ratios
        canv->cd();
        setPadFBStyle();
        hRatioUp->at(i)->Draw();
        if ( drawFits ) {
            fitUp.at(i)->Draw("same");
        }
        if ( !hDown.empty() ) {
            hRatioDown->at(i)->Draw("same");
            if ( drawFits ) {
                fitDown.at(i)->Draw("same");
            }
        }
        hRatioUp->at(i)->GetXaxis()->SetRangeUser(xMin, xMax);
        hRatioUp->at(i)->GetYaxis()->SetRangeUser(0.85, 1.15);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetLow.at(i), ptDijetHi.at(i) ) );
        //t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hRatioUp->at(i), "Up / Default", "p");
        if ( !hDown.empty() ) {
            leg->AddEntry(hRatioDown->at(i), "Down / Default", "p");
        }
        leg->Draw();
        if ( etaType == 0 ) {
            line = new TLine(-3., 1., 3., 1.);
        }
        else {
            line = new TLine(0., 1., 3., 1.);
        }
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        // canv->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sRat_%spt_%d_%d_%s.pdf", date.Data(), systType.Data(), 
        //                   returnTrigName( hDef.at(i)->GetName() ).Data(), systType.Data(), dirName.Data(),
        //                   ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );
        canv->Print( Form("%s/%s/%s_pPb8160_etaDijet_%sRat_%spt_%d_%d_%s.pdf", date.Data(), systType.Data(), 
                          returnTrigName( hDef.at(i)->GetName() ).Data(), systType.Data(), dirName.Data(),
                          ptDijetLow.at(i), ptDijetHi.at(i), frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadFBStyle();
        hDef.at(i)->Draw();
        hUp.at(i)->Draw("same");
        if ( !hDown.empty() ) {
            hDown.at(i)->Draw("same");
        }
        hDef.at(i)->GetYaxis()->SetRangeUser(yMin, yMax);
        hDef.at(i)->GetXaxis()->SetRangeUser(xMin, xMax);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetLow.at(i), ptDijetHi.at(i) ) );
        leg = new TLegend(legXmin, legYmin, legXmax, legYmax);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hDef.at(i), "Default", "p");
        leg->AddEntry(hUp.at(i), "Up", "p");
        if ( !hDown.empty() ) {
            leg->AddEntry(hDown.at(i), "Down", "p");
        }
        leg->Draw();

        // Plot all ratios on one canvas
        cRat->cd(i+1);
        setPadFBStyle();
        hRatioUp->at(i)->Draw();
        if ( drawFits ) {
            fitUp.at(i)->Draw("same");
        }
        if ( !hDown.empty() ) {
            hRatioDown->at(i)->Draw("same");
            if ( drawFits ) {
                fitDown.at(i)->Draw("same");
            }
        }
        hRatioUp->at(i)->GetXaxis()->SetRangeUser(xMin, xMax);
        hRatioUp->at(i)->GetYaxis()->SetRangeUser(0.85, 1.15);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{ave} (GeV) < %d", 
                       ptDijetLow.at(i), ptDijetHi.at(i) ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hRatioUp->at(i), "Up / Default", "p");
        if ( !hDown.empty() ) {
            leg->AddEntry(hRatioDown->at(i), "Down / Default", "p");
        }
        leg->Draw();
        if ( etaType == 0 ) {
            line = new TLine(-3., 1., 3., 1.);
        }
        else {
            line = new TLine(0., 1., 3., 1.);
        }
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();

    } // for (Int_t i{0}; i<ptDijetBinLow.size(); i++)

    cComp->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sComp_%sall_%s.pdf", date.Data(), systType.Data(),  
                   returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/%s/%s_pPb8160_etaDijet_%sRat_%sall_%s.pdf", date.Data(), systType.Data(), 
                   returnTrigName( hDef.at(0)->GetName() ).Data(), systType.Data(), dirName.Data(), frame.Data() ) );
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
             TString date, Bool_t isCM = kFALSE) {

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

    std::vector< TH1D* > hMBJeuUpToDefEta;
    std::vector< TH1D* > hMBJeuDownToDefEta;
    std::vector< TH1D* > hJet60JeuUpToDefEta;
    std::vector< TH1D* > hJet60JeuDownToDefEta;
    std::vector< TH1D* > hJet80JeuUpToDefEta;
    std::vector< TH1D* > hJet80JeuDownToDefEta;
    std::vector< TH1D* > hJet100JeuUpToDefEta;
    std::vector< TH1D* > hJet100JeuDownToDefEta;

    std::vector< TH1D* > hMBJeuUpToDefEtaForward;
    std::vector< TH1D* > hMBJeuDownToDefEtaForward;
    std::vector< TH1D* > hJet60JeuUpToDefEtaForward;
    std::vector< TH1D* > hJet60JeuDownToDefEtaForward;
    std::vector< TH1D* > hJet80JeuUpToDefEtaForward;
    std::vector< TH1D* > hJet80JeuDownToDefEtaForward;
    std::vector< TH1D* > hJet100JeuUpToDefEtaForward;
    std::vector< TH1D* > hJet100JeuDownToDefEtaForward;

    std::vector< TH1D* > hMBJeuUpToDefEtaBackward;
    std::vector< TH1D* > hMBJeuDownToDefEtaBackward;
    std::vector< TH1D* > hJet60JeuUpToDefEtaBackward;
    std::vector< TH1D* > hJet60JeuDownToDefEtaBackward;
    std::vector< TH1D* > hJet80JeuUpToDefEtaBackward;
    std::vector< TH1D* > hJet80JeuDownToDefEtaBackward;
    std::vector< TH1D* > hJet100JeuUpToDefEtaBackward;
    std::vector< TH1D* > hJet100JeuDownToDefEtaBackward;

    std::vector< TH1D* > hMBJeuUpToDefEtaFBRatio;
    std::vector< TH1D* > hMBJeuDownToDefEtaFBRatio;
    std::vector< TH1D* > hJet60JeuUpToDefEtaFBRatio;
    std::vector< TH1D* > hJet60JeuDownToDefEtaFBRatio;
    std::vector< TH1D* > hJet80JeuUpToDefEtaFBRatio;
    std::vector< TH1D* > hJet80JeuDownToDefEtaFBRatio;
    std::vector< TH1D* > hJet100JeuUpToDefEtaFBRatio;
    std::vector< TH1D* > hJet100JeuDownToDefEtaFBRatio;

    std::vector< TH1D* > hMBJeuUpToDefEtaBFRatio;
    std::vector< TH1D* > hMBJeuDownToDefEtaBFRatio;
    std::vector< TH1D* > hJet60JeuUpToDefEtaBFRatio;
    std::vector< TH1D* > hJet60JeuDownToDefEtaBFRatio;
    std::vector< TH1D* > hJet80JeuUpToDefEtaBFRatio;
    std::vector< TH1D* > hJet80JeuDownToDefEtaBFRatio;
    std::vector< TH1D* > hJet100JeuUpToDefEtaBFRatio;
    std::vector< TH1D* > hJet100JeuDownToDefEtaBFRatio;

    plotUp2DownComparison(hMBEta, hJeuMBUpEta, hJeuMBDownEta, &hMBJeuUpToDefEta, &hMBJeuDownToDefEta, 0, 4, date, frame);
    plotUp2DownComparison(hJet60Eta, hJeuJet60UpEta, hJeuJet60DownEta, &hJet60JeuUpToDefEta, &hJet60JeuDownToDefEta, 0, 4, date, frame);
    plotUp2DownComparison(hJet80Eta, hJeuJet80UpEta, hJeuJet80DownEta, &hJet80JeuUpToDefEta, &hJet80JeuDownToDefEta, 0, 4, date, frame);
    plotUp2DownComparison(hJet100Eta, hJeuJet100UpEta, hJeuJet100DownEta, &hJet100JeuUpToDefEta, &hJet100JeuDownToDefEta, 0, 4, date, frame);

    plotUp2DownComparison(hMBEtaForward, hJeuMBUpEtaForward, hJeuMBDownEtaForward, &hMBJeuUpToDefEtaForward, &hMBJeuDownToDefEtaForward, 0, 4, date, frame);
    plotUp2DownComparison(hJet60EtaForward, hJeuJet60UpEtaForward, hJeuJet60DownEtaForward, &hJet60JeuUpToDefEtaForward, &hJet60JeuDownToDefEtaForward, 0, 4, date, frame);
    plotUp2DownComparison(hJet80EtaForward, hJeuJet80UpEtaForward, hJeuJet80DownEtaForward, &hJet80JeuUpToDefEtaForward, &hJet80JeuDownToDefEtaForward, 0, 4, date, frame);
    plotUp2DownComparison(hJet100EtaForward, hJeuJet100UpEtaForward, hJeuJet100DownEtaForward, &hJet100JeuUpToDefEtaForward, &hJet100JeuDownToDefEtaForward, 0, 4, date, frame);

    plotUp2DownComparison(hMBEtaBackward, hJeuMBUpEtaBackward, hJeuMBDownEtaBackward, &hMBJeuUpToDefEtaBackward, &hMBJeuDownToDefEtaBackward, 0, 4, date, frame);
    plotUp2DownComparison(hJet60EtaBackward, hJeuJet60UpEtaBackward, hJeuJet60DownEtaBackward, &hJet60JeuUpToDefEtaBackward, &hJet60JeuDownToDefEtaBackward, 0, 4, date, frame);
    plotUp2DownComparison(hJet80EtaBackward, hJeuJet80UpEtaBackward, hJeuJet80DownEtaBackward, &hJet80JeuUpToDefEtaBackward, &hJet80JeuDownToDefEtaBackward, 0, 4, date, frame);
    plotUp2DownComparison(hJet100EtaBackward, hJeuJet100UpEtaBackward, hJeuJet100DownEtaBackward, &hJet100JeuUpToDefEtaBackward, &hJet100JeuDownToDefEtaBackward, 0, 4, date, frame);

    plotUp2DownComparison(hMBEtaFBRatio, hJeuMBUpEtaFBRatio, hJeuMBDownEtaFBRatio, &hMBJeuUpToDefEtaFBRatio, &hMBJeuDownToDefEtaFBRatio, 0, 4, date, frame);
    plotUp2DownComparison(hJet60EtaFBRatio, hJeuJet60UpEtaFBRatio, hJeuJet60DownEtaFBRatio, &hJet60JeuUpToDefEtaFBRatio, &hJet60JeuDownToDefEtaFBRatio, 0, 4, date, frame);
    plotUp2DownComparison(hJet80EtaFBRatio, hJeuJet80UpEtaFBRatio, hJeuJet80DownEtaFBRatio, &hJet80JeuUpToDefEtaFBRatio, &hJet80JeuDownToDefEtaFBRatio, 0, 4, date, frame);
    plotUp2DownComparison(hJet100EtaFBRatio, hJeuJet100UpEtaFBRatio, hJeuJet100DownEtaFBRatio, &hJet100JeuUpToDefEtaFBRatio, &hJet100JeuDownToDefEtaFBRatio, 0, 4, date, frame);

    plotUp2DownComparison(hMBEtaBFRatio, hJeuMBUpEtaBFRatio, hJeuMBDownEtaBFRatio, &hMBJeuUpToDefEtaBFRatio, &hMBJeuDownToDefEtaBFRatio, 0, 4, date, frame);
    plotUp2DownComparison(hJet60EtaBFRatio, hJeuJet60UpEtaBFRatio, hJeuJet60DownEtaBFRatio, &hJet60JeuUpToDefEtaBFRatio, &hJet60JeuDownToDefEtaBFRatio, 0, 4, date, frame);
    plotUp2DownComparison(hJet80EtaBFRatio, hJeuJet80UpEtaBFRatio, hJeuJet80DownEtaBFRatio, &hJet80JeuUpToDefEtaBFRatio, &hJet80JeuDownToDefEtaBFRatio, 0, 4, date, frame);
    plotUp2DownComparison(hJet100EtaBFRatio, hJeuJet100UpEtaBFRatio, hJeuJet100DownEtaBFRatio, &hJet100JeuUpToDefEtaBFRatio, &hJet100JeuDownToDefEtaBFRatio, 0, 4, date, frame);   
}

//________________
void plotJER(std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hJerUpEtaDist, 
             std::vector< std::vector<TH1D*> > hJerDownEtaDist, TString date, Bool_t isCM = kFALSE) {

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

    std::vector< TH1D* > hEmbeddingJerUpToDefEta;
    std::vector< TH1D* > hEmbeddingJerDownToDefEta;
    std::vector< TH1D* > hEmbeddingJerUpToDefEtaForward;
    std::vector< TH1D* > hEmbeddingJerDownToDefEtaForward;
    std::vector< TH1D* > hEmbeddingJerUpToDefEtaBackward;
    std::vector< TH1D* > hEmbeddingJerDownToDefEtaBackward;
    std::vector< TH1D* > hEmbeddingJerUpToDefEtaFBRatio;
    std::vector< TH1D* > hEmbeddingJerDownToDefEtaFBRatio;
    std::vector< TH1D* > hEmbeddingJerUpToDefEtaBFRatio;
    std::vector< TH1D* > hEmbeddingJerDownToDefEtaBFRatio;

    plotUp2DownComparison(hEmbeddingEta, hJerUpEta, hJerDownEta, &hEmbeddingJerUpToDefEta, &hEmbeddingJerDownToDefEta, 1, 2, date, frame);
    plotUp2DownComparison(hEmbeddingEtaForward, hJerUpEtaForward, hJerDownEtaForward, &hEmbeddingJerUpToDefEtaForward, &hEmbeddingJerDownToDefEtaForward, 1, 2, date, frame);
    plotUp2DownComparison(hEmbeddingEtaBackward, hJerUpEtaBackward, hJerDownEtaBackward, &hEmbeddingJerUpToDefEtaBackward, &hEmbeddingJerDownToDefEtaBackward, 1, 2, date, frame);
    plotUp2DownComparison(hEmbeddingEtaFBRatio, hJerUpEtaFBRatio, hJerDownEtaFBRatio, &hEmbeddingJerUpToDefEtaFBRatio, &hEmbeddingJerDownToDefEtaFBRatio, 1, 2, date, frame);
    plotUp2DownComparison(hEmbeddingEtaBFRatio, hJerUpEtaBFRatio, hJerDownEtaBFRatio, &hEmbeddingJerUpToDefEtaBFRatio, &hEmbeddingJerDownToDefEtaBFRatio, 1, 2, date, frame);
}

//________________
void plotPileup(std::vector< std::vector<TH1D*> > hMBEtaDist, std::vector< std::vector<TH1D*> > hMBGplusEtaDist, 
                std::vector< std::vector<TH1D*> > hMBVtx1EtaDist, std::vector< std::vector<TH1D*> > hJet60EtaDist, 
                std::vector< std::vector<TH1D*> > hJet60GplusEtaDist, std::vector< std::vector<TH1D*> > hJet60Vtx1EtaDist, 
                std::vector< std::vector<TH1D*> > hJet80EtaDist, std::vector< std::vector<TH1D*> > hJet80GplusEtaDist, 
                std::vector< std::vector<TH1D*> > hJet80Vtx1EtaDist, std::vector< std::vector<TH1D*> > hJet100EtaDist, 
                std::vector< std::vector<TH1D*> > hJet100GplusEtaDist, std::vector< std::vector<TH1D*> > hJet100Vtx1EtaDist, 
                TString date, Bool_t isCM = kFALSE) {

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

    std::vector< TH1D* > hMBVtx1Eta = hMBVtx1EtaDist.at(0); std::vector< TH1D* > hMBVtx1EtaForward = hMBVtx1EtaDist.at(1);
    std::vector< TH1D* > hMBVtx1EtaBackward = hMBVtx1EtaDist.at(2); std::vector< TH1D* > hMBVtx1EtaFBRatio = hMBVtx1EtaDist.at(3);
    std::vector< TH1D* > hMBVtx1EtaBFRatio = hMBVtx1EtaDist.at(4);

    std::vector< TH1D* > hJet60Eta = hJet60EtaDist.at(0); std::vector< TH1D* > hJet60EtaForward = hJet60EtaDist.at(1);
    std::vector< TH1D* > hJet60EtaBackward = hJet60EtaDist.at(2); std::vector< TH1D* > hJet60EtaFBRatio = hJet60EtaDist.at(3);
    std::vector< TH1D* > hJet60EtaBFRatio = hJet60EtaDist.at(4);

    std::vector< TH1D* > hJet60GplusEta = hJet60GplusEtaDist.at(0); std::vector< TH1D* > hJet60GplusEtaForward = hJet60GplusEtaDist.at(1);
    std::vector< TH1D* > hJet60GplusEtaBackward = hJet60GplusEtaDist.at(2); std::vector< TH1D* > hJet60GplusEtaFBRatio = hJet60GplusEtaDist.at(3);
    std::vector< TH1D* > hJet60GplusEtaBFRatio = hJet60GplusEtaDist.at(4);

    std::vector< TH1D* > hJet60Vtx1Eta = hJet60Vtx1EtaDist.at(0); std::vector< TH1D* > hJet60Vtx1EtaForward = hJet60Vtx1EtaDist.at(1);
    std::vector< TH1D* > hJet60Vtx1EtaBackward = hJet60Vtx1EtaDist.at(2); std::vector< TH1D* > hJet60Vtx1EtaFBRatio = hJet60Vtx1EtaDist.at(3);
    std::vector< TH1D* > hJet60Vtx1EtaBFRatio = hJet60Vtx1EtaDist.at(4);
    
    std::vector< TH1D* > hJet80Eta = hJet80EtaDist.at(0); std::vector< TH1D* > hJet80EtaForward = hJet80EtaDist.at(1);
    std::vector< TH1D* > hJet80EtaBackward = hJet80EtaDist.at(2); std::vector< TH1D* > hJet80EtaFBRatio = hJet80EtaDist.at(3);
    std::vector< TH1D* > hJet80EtaBFRatio = hJet80EtaDist.at(4);

    std::vector< TH1D* > hJet80GplusEta = hJet80GplusEtaDist.at(0); std::vector< TH1D* > hJet80GplusEtaForward = hJet80GplusEtaDist.at(1);
    std::vector< TH1D* > hJet80GplusEtaBackward = hJet80GplusEtaDist.at(2); std::vector< TH1D* > hJet80GplusEtaFBRatio = hJet80GplusEtaDist.at(3);
    std::vector< TH1D* > hJet80GplusEtaBFRatio = hJet80GplusEtaDist.at(4);

    std::vector< TH1D* > hJet80Vtx1Eta = hJet80Vtx1EtaDist.at(0); std::vector< TH1D* > hJet80Vtx1EtaForward = hJet80Vtx1EtaDist.at(1);
    std::vector< TH1D* > hJet80Vtx1EtaBackward = hJet80Vtx1EtaDist.at(2); std::vector< TH1D* > hJet80Vtx1EtaFBRatio = hJet80Vtx1EtaDist.at(3);
    std::vector< TH1D* > hJet80Vtx1EtaBFRatio = hJet80Vtx1EtaDist.at(4);

    std::vector< TH1D* > hJet100Eta = hJet100EtaDist.at(0); std::vector< TH1D* > hJet100EtaForward = hJet100EtaDist.at(1);
    std::vector< TH1D* > hJet100EtaBackward = hJet100EtaDist.at(2); std::vector< TH1D* > hJet100EtaFBRatio = hJet100EtaDist.at(3);
    std::vector< TH1D* > hJet100EtaBFRatio = hJet100EtaDist.at(4);

    std::vector< TH1D* > hJet100GplusEta = hJet100GplusEtaDist.at(0); std::vector< TH1D* > hJet100GplusEtaForward = hJet100GplusEtaDist.at(1);
    std::vector< TH1D* > hJet100GplusEtaBackward = hJet100GplusEtaDist.at(2); std::vector< TH1D* > hJet100GplusEtaFBRatio = hJet100GplusEtaDist.at(3);
    std::vector< TH1D* > hJet100GplusEtaBFRatio = hJet100GplusEtaDist.at(4);

    std::vector< TH1D* > hJet100Vtx1Eta = hJet100Vtx1EtaDist.at(0); std::vector< TH1D* > hJet100Vtx1EtaForward = hJet100Vtx1EtaDist.at(1);
    std::vector< TH1D* > hJet100Vtx1EtaBackward = hJet100Vtx1EtaDist.at(2); std::vector< TH1D* > hJet100Vtx1EtaFBRatio = hJet100Vtx1EtaDist.at(3);
    std::vector< TH1D* > hJet100Vtx1EtaBFRatio = hJet100Vtx1EtaDist.at(4);

    std::vector< TH1D* > hMBGplusToDefEta;
    std::vector< TH1D* > hMBVtx1ToDefEta;
    std::vector< TH1D* > hJet60GplusToDefEta;
    std::vector< TH1D* > hJet60Vtx1ToDefEta;
    std::vector< TH1D* > hJet80GplusToDefEta;
    std::vector< TH1D* > hJet80Vtx1ToDefEta;
    std::vector< TH1D* > hJet100GplusToDefEta;
    std::vector< TH1D* > hJet100Vtx1ToDefEta;
    
    std::vector< TH1D* > hMBGplusToDefEtaForward;
    std::vector< TH1D* > hMBVtx1ToDefEtaForward;
    std::vector< TH1D* > hJet60GplusToDefEtaForward;
    std::vector< TH1D* > hJet60Vtx1ToDefEtaForward;
    std::vector< TH1D* > hJet80GplusToDefEtaForward;
    std::vector< TH1D* > hJet80Vtx1ToDefEtaForward;
    std::vector< TH1D* > hJet100GplusToDefEtaForward;
    std::vector< TH1D* > hJet100Vtx1ToDefEtaForward;

    std::vector< TH1D* > hMBGplusToDefEtaBackward;
    std::vector< TH1D* > hMBVtx1ToDefEtaBackward;
    std::vector< TH1D* > hJet60GplusToDefEtaBackward;
    std::vector< TH1D* > hJet60Vtx1ToDefEtaBackward;
    std::vector< TH1D* > hJet80GplusToDefEtaBackward;
    std::vector< TH1D* > hJet80Vtx1ToDefEtaBackward;
    std::vector< TH1D* > hJet100GplusToDefEtaBackward;
    std::vector< TH1D* > hJet100Vtx1ToDefEtaBackward;

    std::vector< TH1D* > hMBGplusToDefEtaFBRatio;
    std::vector< TH1D* > hMBVtx1ToDefEtaFBRatio;
    std::vector< TH1D* > hJet60GplusToDefEtaFBRatio;
    std::vector< TH1D* > hJet60Vtx1ToDefEtaFBRatio;
    std::vector< TH1D* > hJet80GplusToDefEtaFBRatio;
    std::vector< TH1D* > hJet80Vtx1ToDefEtaFBRatio;
    std::vector< TH1D* > hJet100GplusToDefEtaFBRatio;
    std::vector< TH1D* > hJet100Vtx1ToDefEtaFBRatio;

    std::vector< TH1D* > hMBGplusToDefEtaBFRatio;
    std::vector< TH1D* > hMBVtx1ToDefEtaBFRatio;
    std::vector< TH1D* > hJet60GplusToDefEtaBFRatio;
    std::vector< TH1D* > hJet60Vtx1ToDefEtaBFRatio;
    std::vector< TH1D* > hJet80GplusToDefEtaBFRatio;
    std::vector< TH1D* > hJet80Vtx1ToDefEtaBFRatio;
    std::vector< TH1D* > hJet100GplusToDefEtaBFRatio;
    std::vector< TH1D* > hJet100Vtx1ToDefEtaBFRatio;

    plotUp2DownComparison(hMBEta, hMBGplusEta, hMBVtx1Eta, &hMBGplusToDefEta, &hMBVtx1ToDefEta, 2, 1, date, frame);
    plotUp2DownComparison(hJet60Eta, hJet60GplusEta, hJet60Vtx1Eta, &hJet60GplusToDefEta, &hJet60Vtx1ToDefEta, 2, 1, date, frame);
    plotUp2DownComparison(hJet80Eta, hJet80GplusEta, hJet80Vtx1Eta, &hJet80GplusToDefEta, &hJet80Vtx1ToDefEta, 2, 1, date, frame);
    plotUp2DownComparison(hJet100Eta, hJet100GplusEta, hJet100Vtx1Eta, &hJet100GplusToDefEta, &hJet100Vtx1ToDefEta, 2, 1, date, frame);
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
                                 std::vector< std::vector<TH1D*> > hEmbeddingEtaDist, std::vector< std::vector<TH1D*> > hEmbeddingGenEtaDist,
                                 std::vector< std::vector<TH1D*> > hJerUpEtaDist, std::vector< std::vector<TH1D*> > hJerDownEtaDist, 
                                 std::vector< std::vector<TH1D*> > hPbGoingEmbeddingRecoEtaDist, std::vector< std::vector<TH1D*> > hPGoingEmbeddingRecoEtaDist, 
                                 std::vector< std::vector<TH1D*> > hPbGoingEmbeddingGenEtaDist, std::vector< std::vector<TH1D*> > hPGoingEmbeddingGenEtaDist,
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
    writeVectors2RootFile(hEmbeddingGenEtaDist);
    writeVectors2RootFile(hJerUpEtaDist);
    writeVectors2RootFile(hJerDownEtaDist);
    writeVectors2RootFile(hPbGoingEmbeddingRecoEtaDist);
    writeVectors2RootFile(hPGoingEmbeddingRecoEtaDist);

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

    oFile->Write();
    oFile->Close();
}

//________________
void applySmoothing(TH1D* h) {
    TString name = h->GetName();
    name.ToLower();
    if ( ( name.Contains("jeu") || name.Contains("pointing") ) && 
           ( !name.Contains("fb") && !name.Contains("forward") && 
             !name.Contains("backward") && !name.Contains("bf") ) ) {
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
            }
        } // for (Int_t i = binLow; i<=binHi; i++)
    } // if ( name.Contains("jeuSyst") || name.Contains("pointingSyst") )
}

//________________
TH1D* calculateSystUncrtBinByBin(TH1D *hUp, TH1D *hDown = nullptr, 
                                const Char_t *setName = "hMB", Int_t systType = 0, Int_t bin = 0) {

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
    TH1D *hSystUncrt = (TH1D*)hUp->Clone( hName.Data() );
    hSystUncrt->Reset();

    // Calculate systematic uncertainties
    for (Int_t i = 1; i <= hUp->GetNbinsX(); i++) {
        Double_t up = hUp->GetBinContent(i);
        Double_t down = up; 
        if ( hDown ) {
            down = hDown->GetBinContent(i);
        }
        Double_t systUncrt = ( TMath::Abs(up - 1) + TMath::Abs(down - 1) ) / 2.0;
        hSystUncrt->SetBinContent(i, systUncrt);
    }
    applySmoothing(hSystUncrt);
    hSystUncrt->SetLineColor(color);
    hSystUncrt->SetMarkerColor(color);
    return hSystUncrt;
}

//________________
std::vector<TH1D*> calculateSystUncrt(std::vector<TH1D*> hUp, std::vector<TH1D*> hDown = {},
                                      const Char_t *setName = "hMB", Int_t systType = 0) {
 
    // Create vector to store systematic uncertainties
    std::vector<TH1D*> hSystUncrt{};

    // Create vectors with pT values
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    Int_t ptBins = ptDijetBinLow.size();
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Calculate systematic uncertainties
    for (Int_t i = 0; i < hUp.size(); i++) {
        hSystUncrt.push_back( calculateSystUncrtBinByBin(hUp.at(i), ( !hDown.empty() ) ? hDown.at(i) : nullptr, setName, systType, i) );
    }

    return hSystUncrt;
}

//________________
std::vector< std::vector<TH1D*> > createRelSystUncrt(std::vector< std::vector<TH1D*> > hUp, std::vector< std::vector<TH1D*> > hDown,
                                                     const Char_t* setName = "hMB", Int_t systType = 0) {

    TString datasetName = setName;
    // Create vector to store relative systematic uncertainties
    std::vector< std::vector<TH1D*> > hRelSystUncrt{};

    // Vectors to process
    std::vector< TH1D* > hUpEta = hUp.at(0); std::vector< TH1D* > hUpEtaForward = hUp.at(1);
    std::vector< TH1D* > hUpEtaBackward = hUp.at(2); std::vector< TH1D* > hUpEtaFBRatio = hUp.at(3);
    std::vector< TH1D* > hUpEtaBFRatio = hUp.at(4);

    std::vector< TH1D* > hDownEta = hDown.at(0); std::vector< TH1D* > hDownEtaForward = hDown.at(1);
    std::vector< TH1D* > hDownEtaBackward = hDown.at(2); std::vector< TH1D* > hDownEtaFBRatio = hDown.at(3);
    std::vector< TH1D* > hDownEtaBFRatio = hDown.at(4);

    // Create relative systematic uncertainties
    datasetName = setName; datasetName += "Eta";         hRelSystUncrt.push_back( calculateSystUncrt(hUpEta, hDownEta, datasetName.Data(), systType) );
    datasetName = setName; datasetName += "EtaForward";  hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaForward, hDownEtaForward, datasetName.Data(), systType) );
    datasetName = setName; datasetName += "EtaBackward"; hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaBackward, hDownEtaBackward, datasetName.Data(), systType) );
    datasetName = setName; datasetName += "EtaFB";       hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaFBRatio, hDownEtaFBRatio, datasetName.Data(), systType) );
    if (systType == 3) {
        std::vector<TH1D*> empty{};
        datasetName = setName; datasetName += "EtaBF";       hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaBFRatio, empty, datasetName.Data(), systType) );
    }
    else {
        datasetName = setName; datasetName += "EtaBF";       hRelSystUncrt.push_back( calculateSystUncrt(hUpEtaBFRatio, hDownEtaBFRatio, datasetName.Data(), systType) );
    }

    return hRelSystUncrt;
}

//________________
std::vector<TH1D*> combineRelSystUncrt(std::vector<TH1D*> hJeuRelSystUncrt, std::vector<TH1D*> hJerRelSystUncrt, 
                                       std::vector<TH1D*> hPileupRelSystUncrt,
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
        for (Int_t j = 1; j <= hJeuRelSystUncrt.at(i)->GetNbinsX(); j++) {
            Double_t jeu = hJeuRelSystUncrt.at(i)->GetBinContent(j);
            Double_t jer = hJerRelSystUncrt.at(i)->GetBinContent(j);
            Double_t pileup = hPileupRelSystUncrt.at(i)->GetBinContent(j);
            Double_t total = TMath::Sqrt( TMath::Power(jeu, 2) + TMath::Power(jer, 2) + TMath::Power(pileup, 2) );
            hTotal->SetBinContent(j, total);
        }
        hTotal->SetLineColor(1);
        hTotal->SetMarkerColor(1);
        hTotalRelSystUncrt.push_back(hTotal);
    }

    return hTotalRelSystUncrt;
}

//________________
std::vector< std::vector<TH1D*> > createTotalRelSystUncrt(std::vector< std::vector<TH1D*> > hJeuRelSystDist, 
                                                          std::vector< std::vector<TH1D*> >hJerRelSystDist, 
                                                          std::vector< std::vector<TH1D*> >hPileupRelSystDist, 
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

    std::vector< TH1D* > hEtaPileupRelSyst = hPileupRelSystDist.at(0);
    std::vector< TH1D* > hEtaForwardPileupRelSyst = hPileupRelSystDist.at(1);
    std::vector< TH1D* > hEtaBackwardPileupRelSyst = hPileupRelSystDist.at(2);
    std::vector< TH1D* > hEtaFBPileupRelSyst = hPileupRelSystDist.at(3);
    std::vector< TH1D* > hEtaBFPileupRelSyst = hPileupRelSystDist.at(4);

    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaJeuRelSyst, hEtaJerRelSyst, hEtaPileupRelSyst, Form("%sEta", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaForwardJeuRelSyst, hEtaForwardJerRelSyst, hEtaForwardPileupRelSyst, Form("%sEtaForward", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaBackwardJeuRelSyst, hEtaBackwardJerRelSyst, hEtaBackwardPileupRelSyst, Form("%sEtaBackward", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaFBJeuRelSyst, hEtaFBJerRelSyst, hEtaFBPileupRelSyst, Form("%sEtaFB", setName) ) );
    hTotalRelSystUncrt.push_back( combineRelSystUncrt(hEtaBFJeuRelSyst, hEtaBFJerRelSyst, hEtaBFPileupRelSyst, Form("%sEtaBF", setName) ) );

    return hTotalRelSystUncrt;
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
                           TString date, Bool_t isCM = kFALSE) {

    TString frame;
    frame = ( isCM ) ? "cms" : "lab";
    TString histoName = "hRecoDijetPtEta";
    if ( isCM ) {
        histoName = "hRecoDijetPtEta";
    }

    TString genHistoName = "hGenDijetPtEta";
    if ( isCM ) {
        genHistoName = "hGenDijetPtEta";
    }

    // Output histogram name
    TString oHistoName;

    // Dijet pT selection
    Int_t ptBins = ptDijetBinLow.size();
    std::vector<Int_t> ptDijetLow{};
    std::vector<Int_t> ptDijetHi{};
    fillDijetPtBins(ptDijetLow, ptDijetHi);

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};
    Int_t genType{4};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

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

    // Create ratios to MC
    std::vector< std::vector<TH1D*> > hRatios2McDist = createRatios2MC(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist, hEmbeddingEtaDist, hEmbeddingGenEtaDist, isCM);
    // Create up and down variation to default ratios
    std::vector< std::vector<TH1D*> > hUpAndDown2DefDist = createUpAndDown2Def(hMBEtaDist, hJeuMBUpEtaDist, hJeuMBDownEtaDist, 
                                                                               hJet60EtaDist, hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, 
                                                                               hJet80EtaDist, hJeuJet80UpEtaDist, hJeuJet80DownEtaDist, 
                                                                               hJet100EtaDist, hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, 
                                                                               hEmbeddingEtaDist, hJerUpEtaDist, hJerDownEtaDist,
                                                                               hMBGplusEtaDist, hMBVtx1EtaDist,
                                                                               hJet60GplusEtaDist, hJet60Vtx1EtaDist,
                                                                               hJet80GplusEtaDist, hJet80Vtx1EtaDist,
                                                                               hJet100GplusEtaDist, hJet100Vtx1EtaDist);
    
    
    // Create JEU up and down histograms
    std::vector< std::vector< TH1D* > > hMBEtaJeuUpDist;
    hMBEtaJeuUpDist.push_back( hUpAndDown2DefDist.at(0) );
    hMBEtaJeuUpDist.push_back( hUpAndDown2DefDist.at(1) );
    hMBEtaJeuUpDist.push_back( hUpAndDown2DefDist.at(2) );
    hMBEtaJeuUpDist.push_back( hUpAndDown2DefDist.at(3) );
    hMBEtaJeuUpDist.push_back( hUpAndDown2DefDist.at(4) );

    std::vector< std::vector< TH1D* > > hMBEtaJeuDownDist;
    hMBEtaJeuDownDist.push_back( hUpAndDown2DefDist.at(5) );
    hMBEtaJeuDownDist.push_back( hUpAndDown2DefDist.at(6) );
    hMBEtaJeuDownDist.push_back( hUpAndDown2DefDist.at(7) );
    hMBEtaJeuDownDist.push_back( hUpAndDown2DefDist.at(8) );
    hMBEtaJeuDownDist.push_back( hUpAndDown2DefDist.at(9) );

    std::vector< std::vector< TH1D* > > hJet60EtaJeuUpDist;
    hJet60EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(10) );
    hJet60EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(11) );
    hJet60EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(12) );
    hJet60EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(13) );
    hJet60EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(14) );

    std::vector< std::vector< TH1D* > > hJet60EtaJeuDownDist;
    hJet60EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(15) );
    hJet60EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(16) );
    hJet60EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(17) );
    hJet60EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(18) );
    hJet60EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(19) );

    std::vector< std::vector< TH1D* > > hJet80EtaJeuUpDist;
    hJet80EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(20) );
    hJet80EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(21) );
    hJet80EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(22) );
    hJet80EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(23) );
    hJet80EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(24) );

    std::vector< std::vector< TH1D* > > hJet80EtaJeuDownDist;
    hJet80EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(25) );
    hJet80EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(26) );
    hJet80EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(27) );
    hJet80EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(28) );
    hJet80EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(29) );

    std::vector< std::vector< TH1D* > > hJet100EtaJeuUpDist;
    hJet100EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(30) );
    hJet100EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(31) );
    hJet100EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(32) );
    hJet100EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(33) );
    hJet100EtaJeuUpDist.push_back( hUpAndDown2DefDist.at(34) );

    std::vector< std::vector< TH1D* > > hJet100EtaJeuDownDist;
    hJet100EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(35) );
    hJet100EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(36) );
    hJet100EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(37) );
    hJet100EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(38) );
    hJet100EtaJeuDownDist.push_back( hUpAndDown2DefDist.at(39) );

    // Create JER up and down histograms
    std::vector< std::vector<TH1D*> > hEmbeddingJerUpDist;
    hEmbeddingJerUpDist.push_back( hUpAndDown2DefDist.at(40) );
    hEmbeddingJerUpDist.push_back( hUpAndDown2DefDist.at(41) );
    hEmbeddingJerUpDist.push_back( hUpAndDown2DefDist.at(42) );
    hEmbeddingJerUpDist.push_back( hUpAndDown2DefDist.at(43) );
    hEmbeddingJerUpDist.push_back( hUpAndDown2DefDist.at(44) );

    std::vector< std::vector<TH1D*> > hEmbeddingJerDownDist;
    hEmbeddingJerDownDist.push_back( hUpAndDown2DefDist.at(45) );
    hEmbeddingJerDownDist.push_back( hUpAndDown2DefDist.at(46) );
    hEmbeddingJerDownDist.push_back( hUpAndDown2DefDist.at(47) );
    hEmbeddingJerDownDist.push_back( hUpAndDown2DefDist.at(48) );
    hEmbeddingJerDownDist.push_back( hUpAndDown2DefDist.at(49) );

    // Create pileup up and down histograms
    std::vector< std::vector<TH1D*> > hMBPileupGplusDist;
    hMBPileupGplusDist.push_back( hUpAndDown2DefDist.at(50) );
    hMBPileupGplusDist.push_back( hUpAndDown2DefDist.at(51) );
    hMBPileupGplusDist.push_back( hUpAndDown2DefDist.at(52) );
    hMBPileupGplusDist.push_back( hUpAndDown2DefDist.at(53) );
    hMBPileupGplusDist.push_back( hUpAndDown2DefDist.at(54) );

    std::vector< std::vector<TH1D*> > hMBPileupVtx1Dist;
    hMBPileupVtx1Dist.push_back( hUpAndDown2DefDist.at(55) );
    hMBPileupVtx1Dist.push_back( hUpAndDown2DefDist.at(56) );
    hMBPileupVtx1Dist.push_back( hUpAndDown2DefDist.at(57) );
    hMBPileupVtx1Dist.push_back( hUpAndDown2DefDist.at(58) );
    hMBPileupVtx1Dist.push_back( hUpAndDown2DefDist.at(59) );

    std::vector< std::vector<TH1D*> > hJet60PileupGplusDist;
    hJet60PileupGplusDist.push_back( hUpAndDown2DefDist.at(60) );
    hJet60PileupGplusDist.push_back( hUpAndDown2DefDist.at(61) );
    hJet60PileupGplusDist.push_back( hUpAndDown2DefDist.at(62) );
    hJet60PileupGplusDist.push_back( hUpAndDown2DefDist.at(63) );
    hJet60PileupGplusDist.push_back( hUpAndDown2DefDist.at(64) );

    std::vector< std::vector<TH1D*> > hJet60PileupVtx1Dist;
    hJet60PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(65) );
    hJet60PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(66) );
    hJet60PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(67) );
    hJet60PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(68) );
    hJet60PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(69) );

    std::vector< std::vector<TH1D*> > hJet80PileupGplusDist;
    hJet80PileupGplusDist.push_back( hUpAndDown2DefDist.at(70) );
    hJet80PileupGplusDist.push_back( hUpAndDown2DefDist.at(71) );
    hJet80PileupGplusDist.push_back( hUpAndDown2DefDist.at(72) );
    hJet80PileupGplusDist.push_back( hUpAndDown2DefDist.at(73) );
    hJet80PileupGplusDist.push_back( hUpAndDown2DefDist.at(74) );

    std::vector< std::vector<TH1D*> > hJet80PileupVtx1Dist;
    hJet80PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(75) );
    hJet80PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(76) );
    hJet80PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(77) );
    hJet80PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(78) );
    hJet80PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(79) );

    std::vector< std::vector<TH1D*> > hJet100PileupGplusDist;
    hJet100PileupGplusDist.push_back( hUpAndDown2DefDist.at(80) );
    hJet100PileupGplusDist.push_back( hUpAndDown2DefDist.at(81) );
    hJet100PileupGplusDist.push_back( hUpAndDown2DefDist.at(82) );
    hJet100PileupGplusDist.push_back( hUpAndDown2DefDist.at(83) );
    hJet100PileupGplusDist.push_back( hUpAndDown2DefDist.at(84) );

    std::vector< std::vector<TH1D*> > hJet100PileupVtx1Dist;
    hJet100PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(85) );
    hJet100PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(86) );
    hJet100PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(87) );
    hJet100PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(88) );
    hJet100PileupVtx1Dist.push_back( hUpAndDown2DefDist.at(89) );


    // Create JEU systematic tables
    std::vector< std::vector<TH1D*> > hMBJeuRelSystDist = createRelSystUncrt(hMBEtaJeuUpDist, hMBEtaJeuDownDist, "hMB", 0);
    std::vector< std::vector<TH1D*> > hJet60JeuRelSystDist = createRelSystUncrt(hJet60EtaJeuUpDist, hJet60EtaJeuDownDist, "hJet60", 0);
    std::vector< std::vector<TH1D*> > hJet80JeuRelSystDist = createRelSystUncrt(hJet80EtaJeuUpDist, hJet80EtaJeuDownDist, "hJet80", 0);
    std::vector< std::vector<TH1D*> > hJet100JeuRelSystDist = createRelSystUncrt(hJet100EtaJeuUpDist, hJet100EtaJeuDownDist, "hJet100", 0);

    // Create JER systematic tables
    std::vector< std::vector<TH1D*> > hEmbeddingJerRelSystDist = createRelSystUncrt(hEmbeddingJerUpDist, hEmbeddingJerDownDist, "hEmbedding", 1);

    // Create pileup systematic tables
    std::vector< std::vector<TH1D*> > hMBPileupRelSystDist = createRelSystUncrt(hMBPileupGplusDist, hMBPileupVtx1Dist, "hMB", 3);
    std::vector< std::vector<TH1D*> > hJet60PileupRelSystDist = createRelSystUncrt(hJet60PileupGplusDist, hJet60PileupVtx1Dist, "hJet60", 3);
    std::vector< std::vector<TH1D*> > hJet80PileupRelSystDist = createRelSystUncrt(hJet80PileupGplusDist, hJet80PileupVtx1Dist, "hJet80", 3);
    std::vector< std::vector<TH1D*> > hJet100PileupRelSystDist = createRelSystUncrt(hJet100PileupGplusDist, hJet100PileupVtx1Dist, "hJet100", 3);

    // Calculate totaly systematic uncertainty
    std::vector< std::vector<TH1D*> > hMBTotalRelSystDist = createTotalRelSystUncrt(hMBJeuRelSystDist, hEmbeddingJerRelSystDist, hMBPileupRelSystDist, "hMB");
    std::vector< std::vector<TH1D*> > hJet60TotalRelSystDist = createTotalRelSystUncrt(hJet60JeuRelSystDist, hEmbeddingJerRelSystDist, hJet60PileupRelSystDist, "hJet60");
    std::vector< std::vector<TH1D*> > hJet80TotalRelSystDist = createTotalRelSystUncrt(hJet80JeuRelSystDist, hEmbeddingJerRelSystDist, hJet80PileupRelSystDist, "hJet80");
    std::vector< std::vector<TH1D*> > hJet100TotalRelSystDist = createTotalRelSystUncrt(hJet100JeuRelSystDist, hEmbeddingJerRelSystDist, hJet100PileupRelSystDist, "hJet100");

    // Plot final distributions
    plotFinalEtaDistributions(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist, 
                              hMBTotalRelSystDist, hJet60TotalRelSystDist, 
                              hJet80TotalRelSystDist, hJet100TotalRelSystDist, 
                              date, isCM);

    // Plot comparison of data and Monte Carlo
    // plotData2McComparison(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist, hEmbeddingEtaDist, hEmbeddingGenEtaDist, hRatios2McDist, date, isCM);

    // plotJES(hMBEtaDist, hJet60EtaDist, hJet80EtaDist, hJet100EtaDist,
    //         hJeuMBUpEtaDist, hJeuMBDownEtaDist,
    //         hJeuJet60UpEtaDist, hJeuJet60DownEtaDist, 
    //         hJeuJet80UpEtaDist, hJeuJet80DownEtaDist,
    //         hJeuJet100UpEtaDist, hJeuJet100DownEtaDist, 
    //         date, isCM);

    // plotJER(hEmbeddingEtaDist, hJerUpEtaDist, hJerDownEtaDist, date, isCM);

    // plotPileup(hMBEtaDist, hMBGplusEtaDist, hMBVtx1EtaDist,
    //            hJet60EtaDist, hJet60GplusEtaDist, hJet60Vtx1EtaDist,
    //            hJet80EtaDist, hJet80GplusEtaDist, hJet80Vtx1EtaDist,
    //            hJet100EtaDist, hJet100GplusEtaDist, hJet100Vtx1EtaDist,
    //            date, isCM);

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
                                hEmbeddingEtaDist, hEmbeddingGenEtaDist, 
                                hJerUpEtaDist, hJerDownEtaDist, 
                                hPbGoingEmbeddingRecoEtaDist, hPGoingEmbeddingRecoEtaDist,
                                hPbGoingEmbeddingGenEtaDist, hPGoingEmbeddingGenEtaDist,
                                hMBGplusEtaDist, hMBVtx1EtaDist, 
                                hJet60GplusEtaDist, hJet60Vtx1EtaDist, 
                                hJet80GplusEtaDist, hJet80Vtx1EtaDist, 
                                hJet100GplusEtaDist, hJet100Vtx1EtaDist, 
                                hRatios2McDist, hUpAndDown2DefDist,
                                hMBJeuRelSystDist, hJet60JeuRelSystDist, hJet80JeuRelSystDist, hJet100JeuRelSystDist,
                                hEmbeddingJerRelSystDist,
                                hMBPileupRelSystDist, hJet60PileupRelSystDist, hJet80PileupRelSystDist, hJet100PileupRelSystDist,
                                hMBTotalRelSystDist, hJet60TotalRelSystDist, hJet80TotalRelSystDist, hJet100TotalRelSystDist,
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
void runOldCalculations(Bool_t isCM = kFALSE) {

    Bool_t drawFits{kTRUE};
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

    // Retrieve the base branch name to work on
    // Int_t branchId{0};
    // if ( defaultFileName.Contains("akCs4") ) {
    //     branchId = {0};
    // }
    // else {
    //     branchId = {1};
    // }

    Int_t trigVal{0}; // 0-MB, 1-Jet60, 2-Jet80, 3-Jet100

    TString trigName;
    if ( trigVal == 0 ) {
        trigName = "MB";
    }
    else if ( trigVal == 1 ) {
        trigName = "Jet60";
    }
    else if ( trigVal == 2 ) {
        trigName = "Jet80";
    }
    else {
        trigName = "Jet100";
    }

    TString pathPrefix = "/Users/nigmatkulov/cernbox";

    // List of file names
    TString defaultFileName( Form("%s/ana/pPb8160/exp/%s_pPb8160_ak4%s.root", pathPrefix.Data(), trigName.Data(), cmStr.Data()) );
    TString pbGoingFileName( Form("%s/ana/pPb8160/exp/Pbgoing/def/%s_Pbgoing_pPb8160_ak4%s.root", pathPrefix.Data(), trigName.Data(), cmStr.Data()) );
    TString pGoingFileName( Form("%s/ana/pPb8160/exp/pgoing/def/%s_pgoing_pPb8160_ak4%s.root", pathPrefix.Data(), trigName.Data(), cmStr.Data()) );

    TString pbGoingEmbeddingFileName( Form("%s/ana/pPb8160/embedding/Pbgoing/jer/oEmbedding_pPb8160_Pbgoing_jerDef_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString pGoingEmbeddingFileName( Form("%s/ana/pPb8160/embedding/pgoing/jer/oEmbedding_pPb8160_pgoing_jerDef_ak4%s.root", pathPrefix.Data() , cmStr.Data()) );

    TString akcs4FileName( Form("%s/ana/pPb8160/exp/MB_pPb8160_akCs4.root", pathPrefix.Data(), trigName.Data()) );

    TString jeuUpFileName( Form("%s/ana/pPb8160/exp/%s_pPb8160_jeu_up_ak4%s.root", pathPrefix.Data(), trigName.Data(), cmStr.Data()) );
    TString jeuDownFileName( Form("%s/ana/pPb8160/exp/%s_pPb8160_jeu_down_ak4%s.root", pathPrefix.Data(), trigName.Data(), cmStr.Data()) );
    TString embeddingFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );

    // TString jerDefFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_weight%s_ak4.root", pathPrefix.Data(), trigName.Data() ) );
    // TString jerUpFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerUp_weight%s_ak4.root", pathPrefix.Data(), trigName.Data() ) );
    // TString jerDownFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDown_weight%s_ak4.root", pathPrefix.Data(), trigName.Data() ) );
    TString jerDefFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jerUpFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerUp_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );
    TString jerDownFileName( Form("%s/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDown_ak4%s.root", pathPrefix.Data(), cmStr.Data()) );

    TString pileupGplusFileName( Form("%s/ana/pPb8160/exp/%s_pPb8160_gplus_ak4%s.root", pathPrefix.Data(), trigName.Data(), cmStr.Data()) );
    TString pileupVtx1FileName( Form("%s/ana/pPb8160/exp/%s_pPb8160_vtx1_ak4.root", pathPrefix.Data(), trigName.Data()) );

    std::cout << "Default file: " << defaultFileName.Data() << std::endl;
    std::cout << "PbGoing file: " << pbGoingFileName.Data() << std::endl;
    std::cout << "PGoing file: " << pGoingFileName.Data() << std::endl;
    std::cout << "JeuUp file: " << jeuUpFileName.Data() << std::endl;
    std::cout << "JeuDown file: " << jeuDownFileName.Data() << std::endl;


    // Files to work on
    TFile *defaultFile = TFile::Open( defaultFileName.Data() );
    TFile *pbGoingFile = TFile::Open( pbGoingFileName.Data() );
    TFile *pGoingFile = TFile::Open( pGoingFileName.Data() );
    TFile *akcs4File = TFile::Open( akcs4FileName.Data() );
    TFile *jeuUpFile = TFile::Open( jeuUpFileName.Data() );
    TFile *jeuDownFile = TFile::Open( jeuDownFileName.Data() );
    TFile *embeddingFile = TFile::Open( embeddingFileName.Data() );
    TFile *jerDefFile = TFile::Open( jerDefFileName.Data() );
    TFile *jerUpFile = TFile::Open( jerUpFileName.Data() );
    TFile *jerDownFile = TFile::Open( jerDownFileName.Data() );
    TFile *pbGoingEmbeddingFile = TFile::Open( pbGoingEmbeddingFileName.Data() );
    TFile *pGoingEmbeddingFile = TFile::Open( pGoingEmbeddingFileName.Data() );
    TFile *gplusFile = TFile::Open( pileupGplusFileName.Data() );
    TFile *vtx1File = TFile::Open( pileupVtx1FileName.Data() );

    // Checks files are okay
    if ( !checkFileIsGood( defaultFile ) ) return;
    if ( !checkFileIsGood( pbGoingFile ) ) return;
    if ( !checkFileIsGood( pGoingFile ) ) return;
    if ( !checkFileIsGood( akcs4File ) ) return;
    if ( !checkFileIsGood( jeuUpFile ) ) return;
    if ( !checkFileIsGood( jeuDownFile ) ) return;
    if ( !checkFileIsGood( embeddingFile ) ) return;
    if ( !checkFileIsGood( jerDefFile ) ) return;
    if ( !checkFileIsGood( jerUpFile ) ) return;
    if ( !checkFileIsGood( jerDownFile ) ) return;
    if ( !checkFileIsGood( gplusFile ) ) return;
    if ( !checkFileIsGood( vtx1File ) ) return;


    // plotDifferentDirections( pbGoingFile, pGoingFile, isCM, date );
    // plotDifferentDirections( pbGoingEmbeddingFile, pGoingEmbeddingFile, date );
    // compareData2McDifferentDirections(pbGoingFile, pGoingFile, pbGoingEmbeddingFile, pGoingEmbeddingFile, isCM, date, defaultFile);
    // plotJEU( defaultFile, jeuUpFile, jeuDownFile, defaultFile, date, isCM, drawFits );
    // plotJER(jerDefFile, jerUpFile, jerDownFile, date, isCM, drawFits);
    // plotPointingResolution( jerDefFile, date, isCM, drawFits );
    // plotPileup( defaultFile, gplusFile, vtx1File, date, isCM, drawFits );
    // compareJetCollections( defaultFile, akcs4File, date );
    // compareJER( embeddingFile, jerDefFile, jerUpFile, jerDownFile, date );
}

//________________
void runNewCalculations(Bool_t isCM = kFALSE) {

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

    TString pathPrefix = "/Users/nigmatkulov/cernbox";

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
                          date, isCM);
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

    // Run old calculations
    // runOldCalculations(isCM);

    // Run new calculations
    runNewCalculations(isCM);
}