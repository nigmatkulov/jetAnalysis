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

// C++ headers
#include <iostream>
#include <vector>


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
void createDirectory(const char* directoryPath) {
    // Create the directory with read, write, and execute permissions for owner, group, and others
    if (mkdir(directoryPath, S_IRWXU | S_IRWXG | S_IRWXO) == 0) {
        std::cout << "Directory created successfully." << std::endl;
    } else {
        std::cerr << "Failed to create directory." << std::endl;
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
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
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
        hProjY[i]->SetNameTitle( Form("%s_%d", hNameY, i), ";p_{T} (GeV/c)");
        set1DStyle( hProjY[i], style );
        hProjY[i]->SetMarkerSize(0.7);
    } // for (Int_t i{1}; i<=nBinsY; i++)
}

//________________
void make1DRatio(TH1D *hRat, TH1D *hDen, const Char_t *ratioName = "Ratio to default", Int_t style = 0) {

    if ( !hDen ) {
        std::cout << "Denominator does not exist" << std::endl;
    }

    hRat->Divide( hRat, hDen, 1., 1., "b" );
    hRat->GetYaxis()->SetTitle( ratioName );
    hRat->GetYaxis()->SetRangeUser(0.8, 1.2);
    set1DStyle(hRat, style);
}

//________________
TH1D* projectEtaFrom3D(TH3D* h3D, const Char_t *name, 
                       Int_t ptDijetLow, Int_t ptDijetHi) {

    TH1D *tmp = dynamic_cast<TH1D*>( h3D->ProjectionY( name, ptDijetLow, ptDijetHi ) );
    
    tmp->GetYaxis()->SetTitle("1/N dN/d#eta^{dijet}");
    return tmp;
}

//________________
void plotDifferentDirections(TFile *pbGoingFile, TFile *pGoingFile, TString date) {

    TH3D *hPtEtaDphiPbGoing = (TH3D*)pbGoingFile->Get("hRecoDijetPtEtaDphiWeighted");
    //TH3D *hPtEtaDphiPbGoing = (TH3D*)pbGoingFile->Get("hRecoDijetPtEtaDphi");
    hPtEtaDphiPbGoing->SetName("hPtEtaDphiPbGoing");
    TH3D *hPtEtaDphiPGoing = (TH3D*)pGoingFile->Get("hRecoDijetPtEtaDphiWeighted");
    //TH3D *hPtEtaDphiPGoing = (TH3D*)pGoingFile->Get("hRecoDijetPtEtaDphi");
    hPtEtaDphiPGoing->SetName("hPtEtaDphiPGoing");

    TString frame = "lab";

    // Double_t dijetPtVals[dijetPtBins+1] {  40.,  50.,   60.,  70.,  80.,
    //                                        90., 100.,  110., 120., 130.,
    //                                       140., 150.,  160., 180., 200., 
    //                                       240., 300., 1000.};
    // fDijetPtBins{194}, fDijetPtRange{30., 1000.}

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // Styles
    Int_t pbGoingType{0};
    Int_t pGoingType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaPbGoing[ ptDijetLow.size() ];
    TH1D *hEtaPGoing[ ptDijetLow.size() ];
    TH1D *hEtaPb2PGoingRatio[ ptDijetLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    std::cout << "ptBins: " << ptBins << std::endl;

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        // Pb-going projection
        hEtaPbGoing[i] = projectEtaFrom3D(hPtEtaDphiPbGoing, Form("hEtaPbGoing_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaPbGoing[i] );
        set1DStyle(hEtaPbGoing[i], pbGoingType);

        // p-going projection
        hEtaPGoing[i] = projectEtaFrom3D(hPtEtaDphiPGoing, Form("hEtaPGoing_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
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
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaPbGoing[i], Form("Pb-going"), "p");
        leg->AddEntry(hEtaPGoing[i], Form("p-going"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_dirComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios of Pb-going to p-going
        canv->cd();
        setPadStyle();
        hEtaPb2PGoingRatio[i]->Draw();
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hEtaPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_dirRat_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hEtaPbGoing[i]->Draw();
        hEtaPGoing[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
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
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hEtaPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
    } // for (Int_t i=0; i<ptDijetLow.size(); i++)

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_dirComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_dirRat_all_%s.pdf", date.Data(), frame.Data() ) );
}


//________________
void compareData2McDifferentDirections(TFile *expPbGoing, TFile *expPGoing, 
                                       TFile *mcPbGoing, TFile *mcPGoing, TString date) {

    TH2D *hExpPtEtaPbGoing = dynamic_cast<TH2D*>( expPbGoing->Get("hRecoDijetPtEta") );
    hExpPtEtaPbGoing->SetName("hExpPtEtaPbGoing");
    TH2D *hExpPtEtaPGoing = dynamic_cast<TH2D*>( expPGoing->Get("hRecoDijetPtEta") );
    hExpPtEtaPGoing->SetName("hExpPtEtaPGoing");

    TH2D *hMcPtEtaPbGoing = dynamic_cast<TH2D*>( mcPbGoing->Get("hRecoDijetPtEta") );
    hMcPtEtaPbGoing->SetName("hMcPtEtaPbGoing");
    TH2D *hMcPtEtaPGoing = dynamic_cast<TH2D*>( mcPGoing->Get("hRecoDijetPtEta") );
    hMcPtEtaPGoing->SetName("hMcPtEtaPGoing");

    TString frame = "lab";

    // Do rebinning if needed
    Int_t rebinX{1}, rebinY{4};
    hExpPtEtaPbGoing->RebinX( rebinX );
    hExpPtEtaPbGoing->RebinY( rebinY );
    hExpPtEtaPGoing->RebinX( rebinX );
    hExpPtEtaPGoing->RebinY( rebinY );
    hMcPtEtaPbGoing->RebinX( rebinX );
    hMcPtEtaPbGoing->RebinY( rebinY );
    hMcPtEtaPGoing->RebinX( rebinX );
    hMcPtEtaPGoing->RebinY( rebinY );

    // Double_t dijetPtVals[dijetPtBins+1] {  40.,  50.,   60.,  70.,  80.,
    //                                        90., 100.,  110., 120., 130.,
    //                                       140., 150.,  160., 180., 200., 
    //                                       240., 300., 1000.};
    // fDijetPtBins{194}, fDijetPtRange{30., 1000.}

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

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
    TH1D *hExpEtaPbGoing[ ptDijetLow.size() ];
    TH1D *hExpEtaPGoing[ ptDijetLow.size() ];
    TH1D *hExpEtaPb2PGoingRatio[ ptDijetLow.size() ];
    TH1D *hMcEtaPbGoing[ ptDijetLow.size() ];
    TH1D *hMcEtaPGoing[ ptDijetLow.size() ];
    TH1D *hMcEtaPb2PGoingRatio[ ptDijetLow.size() ];
    TH1D *hExpOverMcPb2PGoingRatio[ ptDijetLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        // Exp Pb-going projection
        hExpEtaPbGoing[i] = dynamic_cast<TH1D*>( hExpPtEtaPbGoing->ProjectionY( Form("hExpEtaPbGoing_%d", i), 
                                                                                ptDijetLow.at(i), ptDijetHi.at(i) ) );
        rescaleEta( hExpEtaPbGoing[i] );
        set1DStyle(hExpEtaPbGoing[i], expPbGoingType);

        // Exp p-going projection
        hExpEtaPGoing[i] = dynamic_cast<TH1D*>( hExpPtEtaPGoing->ProjectionY( Form("hExpEtaPGoing_%d", i), 
                                                                              ptDijetLow.at(i), ptDijetHi.at(i) ));
        rescaleEta( hExpEtaPGoing[i] );
        set1DStyle(hExpEtaPGoing[i], expPGoingType);

        // Ratio of Pb-going to p-going
        hExpEtaPb2PGoingRatio[i] = dynamic_cast<TH1D*>( hExpEtaPbGoing[i]->Clone( Form("hExpEtaPb2PGoingRatio_%d", i) ) );
        make1DRatio(hExpEtaPb2PGoingRatio[i], hExpEtaPGoing[i], Form("hExpEtaPb2PGoingRatio_%d", i) );
        hExpEtaPb2PGoingRatio[i]->GetYaxis()->SetTitle("Pb-going / p-going");
        set1DStyle( hExpEtaPb2PGoingRatio[i], expPbGoingType);

        // Mc Pb-going projection
        hMcEtaPbGoing[i] = dynamic_cast<TH1D*>( hMcPtEtaPbGoing->ProjectionY( Form("hMcEtaPbGoing_%d", i), 
                                                                                ptDijetLow.at(i), ptDijetHi.at(i) ) );
        rescaleEta( hMcEtaPbGoing[i] );
        set1DStyle(hMcEtaPbGoing[i], mcPbGoingType);

        // Mc p-going projection
        hMcEtaPGoing[i] = dynamic_cast<TH1D*>( hMcPtEtaPGoing->ProjectionY( Form("hMcEtaPGoing_%d", i), 
                                                                              ptDijetLow.at(i), ptDijetHi.at(i) ));
        rescaleEta( hMcEtaPGoing[i] );
        set1DStyle(hMcEtaPGoing[i], mcPGoingType);

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

        // Plot comparison
        canv->cd();
        setPadStyle();
        hExpEtaPb2PGoingRatio[i]->Draw();
        hMcEtaPb2PGoingRatio[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hExpEtaPb2PGoingRatio[i], Form("Data"), "p");
        leg->AddEntry(hMcEtaPb2PGoingRatio[i], Form("MC"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_exp2mc_dirComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios of Pb-going to p-going
        canv->cd();
        setPadStyle();
        hExpOverMcPb2PGoingRatio[i]->Draw();
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hExpOverMcPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_ exp2mc_dirRat_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hExpEtaPb2PGoingRatio[i]->Draw();
        hMcEtaPb2PGoingRatio[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
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
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hExpOverMcPb2PGoingRatio[i]->GetXaxis()->GetBinUpEdge(hExpOverMcPb2PGoingRatio[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
    } // for (Int_t i=0; i<ptDijetLow.size(); i++)

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_exp2mc_dirComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_exp2mc_dirRat_all_%s.pdf", date.Data(), frame.Data() ) );
}

//________________
void plotJEU(TFile *defaultFile, TFile *jeuUpFile, TFile *jeuDownFile, TString date) {

    TH3D *hPtEtaDphiDef = (TH3D*)defaultFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiDef->SetName("hPtEtaDphiDef");
    TH3D *hPtEtaDphiUp = (TH3D*)jeuUpFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiUp->SetName("hPtEtaDphiUp");
    TH3D *hPtEtaDphiDown = (TH3D*)jeuDownFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiDown->SetName("hPtEtaDphiDown");

    TString frame = "lab";

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaDef[ ptDijetLow.size() ];
    TH1D *hEtaUp[ ptDijetLow.size() ];
    TH1D *hEtaDown[ ptDijetLow.size() ];

    TH1D *hEtaRatioUp[ ptDijetLow.size() ];
    TH1D *hEtaRatioDown[ ptDijetLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // Nominal selection
        hEtaDef[i] = projectEtaFrom3D(hPtEtaDphiDef, Form("hEtaDef_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaDef[i] );
        set1DStyle(hEtaDef[i], defType);

        // JEU up selection
        hEtaUp[i] = projectEtaFrom3D(hPtEtaDphiUp, Form("hEtaUp_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaUp[i] );
        set1DStyle(hEtaUp[i], upType);

        // JEU down selection
        hEtaDown[i] = projectEtaFrom3D(hPtEtaDphiDown, Form("hEtaDown_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaDown[i] );
        set1DStyle(hEtaDown[i], downType);

        // Ratio of JEU up to default ratio
        hEtaRatioUp[i] = (TH1D*)hEtaUp[i]->Clone( Form("hEtaRatioUp_%d", i) );
        make1DRatio(hEtaRatioUp[i], hEtaDef[i], Form("hEtaRatioUp_%d", i), upType );
        hEtaRatioUp[i]->GetYaxis()->SetTitle("JEU / Default");

        // Ratio of JEU down to default ratio
        hEtaRatioDown[i] = (TH1D*)hEtaDown[i]->Clone( Form("hEtaRatioDown_%d", i) );
        make1DRatio(hEtaRatioDown[i], hEtaDef[i], Form("hEtaRatioDown_%d", i), downType );
        hEtaRatioDown[i]->GetYaxis()->SetTitle("JEU / Default");

        // Plot comparison
        canv->cd();
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("Default"), "p");
        leg->AddEntry(hEtaUp[i], Form("JEU+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JEU-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_jeuComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios
        canv->cd();
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("JEU+ / Default"), "p");
        leg->AddEntry(hEtaRatioDown[i], Form("JEU- / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_jeuRat_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
                hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("Default"), "p");
        leg->AddEntry(hEtaUp[i], Form("JEU+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JEU-"), "p");
        leg->Draw();

        // Plot all rations on one canvas
        cRat->cd(i+1);
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.45, 0.2, 0.65, 0.3);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaRatioUp[i], Form("JEU+ / Default"), "p");
        leg->AddEntry(hEtaRatioDown[i], Form("JEU- / Default"), "p");
        leg->Draw();
        line = new TLine(hEtaRatioUp[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRatioUp[i]->GetXaxis()->GetBinUpEdge(hEtaRatioUp[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
    } // for (Int_t i=0; i<ptDijetLow.size(); i++)

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_jeuComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_jeuRat_all_%s.pdf", date.Data(), frame.Data() ) );
}

//________________
void plotJER(TFile *defaultFile, TFile *jerUpFile, TFile *jerDownFile, TString date) {

    TH3D *hPtEtaDphiDef = (TH3D*)defaultFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiDef->SetName("hPtEtaDphiDef");
    TH3D *hPtEtaDphiUp = (TH3D*)jerUpFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiUp->SetName("hPtEtaDphiUp");
    TH3D *hPtEtaDphiDown = (TH3D*)jerDownFile->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiDown->SetName("hPtEtaDphiDown");

    TString frame = "lab";

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // Styles
    Int_t defType{2};
    Int_t upType{0};
    Int_t downType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaDef[ ptDijetLow.size() ];
    TH1D *hEtaUp[ ptDijetLow.size() ];
    TH1D *hEtaDown[ ptDijetLow.size() ];

    TH1D *hEtaRatioUp[ ptDijetLow.size() ];
    TH1D *hEtaRatioDown[ ptDijetLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // JER default selection
        hEtaDef[i] = projectEtaFrom3D(hPtEtaDphiDef, Form("hEtaDef_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaDef[i] );
        set1DStyle(hEtaDef[i], defType);

        // JER up selection
        hEtaUp[i] = projectEtaFrom3D(hPtEtaDphiUp, Form("hEtaUp_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaUp[i] );
        set1DStyle(hEtaUp[i], upType);

        // JER down selection
        hEtaDown[i] = projectEtaFrom3D(hPtEtaDphiDown, Form("hEtaDown_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaDown[i] );
        set1DStyle(hEtaDown[i], downType);

        // Ratio of JER up to default ratio
        hEtaRatioUp[i] = dynamic_cast<TH1D*>( hEtaUp[i]->Clone( Form("hEtaRatioUp_%d", i) ) );
        make1DRatio(hEtaRatioUp[i], hEtaDef[i], Form("hEtaRatioUp_%d", i), upType );
        hEtaRatioUp[i]->GetYaxis()->SetTitle("JER / Default");

        // Ratio of JER down to default ratio
        hEtaRatioDown[i] = dynamic_cast<TH1D*>( hEtaDown[i]->Clone( Form("hEtaRatioDown_%d", i) ) );
        make1DRatio(hEtaRatioDown[i], hEtaDef[i], Form("hEtaRatioDown_%d", i), downType );
        hEtaRatioDown[i]->GetYaxis()->SetTitle("JER / Default");

        // Plot comparison
        canv->cd();
        setPadStyle();
        hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("JER def."), "p");
        leg->AddEntry(hEtaUp[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JER-"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_jerComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratios
        canv->cd();
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
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
        canv->SaveAs( Form("%s/pPb8160_etaDijet_jerRat_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
                hEtaDef[i]->Draw();
        hEtaUp[i]->Draw("same");
        hEtaDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaDef[i], Form("JER def."), "p");
        leg->AddEntry(hEtaUp[i], Form("JER+"), "p");
        leg->AddEntry(hEtaDown[i], Form("JER-"), "p");
        leg->Draw();

        // Plot all rations on one canvas
        cRat->cd(i+1);
        setPadStyle();
        hEtaRatioUp[i]->Draw();
        hEtaRatioDown[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
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
    } // for (Int_t i{0}; i<ptBins; i++)

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_jerComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_jerRat_all_%s.pdf", date.Data(), frame.Data() ) );
}

//________________
void plotPointingResolution(TFile *embeddingFile, TString date) {

    TH3D *hRecoEtaRefEtaRecoPt = (TH3D*)embeddingFile->Get("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted");
    hRecoEtaRefEtaRecoPt->SetName("hRecoEtaRefEtaRecoPt");

    TString frame = "lab";

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // Styles
    Int_t recoType{0};
    Int_t refType{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaReco[ ptDijetLow.size() ];
    TH1D *hEtaRef[ ptDijetLow.size() ];
    TH1D *hEtaDiff[ ptDijetLow.size() ];
    TH2D *hEtaRefVsEtaReco[ ptDijetLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};
    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *c2D = new TCanvas("c2D", "c2D", sizeX, sizeY);
    c2D->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // Project to 2D (ref vs reco)
        hRecoEtaRefEtaRecoPt->GetZaxis()->SetRange( ptDijetLow.at(i), ptDijetHi.at(i) );
        hEtaRefVsEtaReco[i] = dynamic_cast<TH2D*>( hRecoEtaRefEtaRecoPt->Project3D("yx") );
        hEtaRefVsEtaReco[i]->SetNameTitle( Form("hEtaRefVsEtaReco_%d", i),";#eta^{dijet}_{reco};#eta^{dijet}_{ref}" );
        rescaleEta( hEtaRefVsEtaReco[i] );
        set2DStyle( hEtaRefVsEtaReco[i] );

        // Project on eta reco
        hEtaReco[i] = dynamic_cast<TH1D*>(hRecoEtaRefEtaRecoPt->ProjectionX( Form("hEtaReco_%d", i), 1, -1, ptDijetLow.at(i), ptDijetHi.at(i) ) );
        rescaleEta( hEtaReco[i] );
        set1DStyle( hEtaReco[i], recoType );
        hEtaReco[i]->GetXaxis()->SetTitle("#eta");

        // Project on eta ref
        hEtaRef[i] = dynamic_cast<TH1D*>(hRecoEtaRefEtaRecoPt->ProjectionY( Form("hEtaRef_%d", i), 1, -1, ptDijetLow.at(i), ptDijetHi.at(i) ) );
        rescaleEta( hEtaRef[i] );
        set1DStyle( hEtaRef[i], refType );
        hEtaRef[i]->GetXaxis()->SetTitle("#eta");

        // Take difference of reco and ref
        hEtaDiff[i] = dynamic_cast<TH1D*>( hEtaReco[i]->Clone( Form( "hEtaDiff_%d", i ) ) );
        hEtaDiff[i]->Add(hEtaRef[i], -1.);
        hEtaDiff[i]->GetYaxis()->SetTitle("reco - ref");

        // Plot 2D matrix
        canv->cd();
        setPadStyle();
        hEtaRefVsEtaReco[i]->Draw("colz");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/pPb8160_etaDijet_RecoRefComp_2D_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot comparison
        canv->cd();
        setPadStyle();
        hEtaReco[i]->Draw();
        hEtaRef[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaReco[i], Form("Reco"), "p");
        leg->AddEntry(hEtaRef[i], Form("Ref"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_RecoRefComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot difference
        canv->cd();
        setPadStyle();
        hEtaDiff[i]->Draw();
        hEtaDiff[i]->GetYaxis()->SetRangeUser(-0.002, 0.002);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        canv->SaveAs( Form("%s/pPb8160_etaDijet_RecoRefDiff_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot all 2D distributions on one canvas
        c2D->cd(i+1);
        setPadStyle();
        hEtaRefVsEtaReco[i]->Draw("colz");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hEtaReco[i]->Draw();
        hEtaRef[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.65, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaReco[i], Form("Reco"), "p");
        leg->AddEntry(hEtaRef[i], Form("Ref"), "p");
        leg->Draw();

        // Plot all differences on one canvas
        cRat->cd(i+1);
        setPadStyle();
        hEtaDiff[i]->Draw();
        hEtaDiff[i]->GetYaxis()->SetRangeUser(-0.002, 0.002);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
    } // for (Int_t i{0}; i<ptBins; i++) 

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_RecoRefComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_RecoRefDiff_all_%s.pdf", date.Data(), frame.Data() ) );
    c2D->SaveAs( Form("%s/pPb8160_etaDijet_RecoRefCorr_2D_all_%s.pdf", date.Data(), frame.Data() ) );
}

//________________
void compareJetCollections(TFile *ak4, TFile *akCs4, TString date) {

    TH3D *hPtEtaDphiAk4 = (TH3D*)ak4->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiAk4->SetName("hPtEtaDphiAk4");
    TH3D *hPtEtaDphiAkCs4 = (TH3D*)akCs4->Get("hRecoDijetPtEtaDphiWeighted");
    hPtEtaDphiAkCs4->SetName("hPtEtaDphiAkCs4");

    TString frame = "lab";

    // Dijet pT selection
    Int_t ptStep {5};
    Int_t ptLow {30};
    std::vector<Int_t> ptDijetLow {3, 5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 43, 55};
    std::vector<Int_t> ptDijetHi  {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 42, 54, 194};

    // Styles
    Int_t ak4Type{0};
    Int_t akCs4Type{1};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.06);

    // Dijet eta distributions
    TH1D *hEtaAk4[ ptDijetLow.size() ];
    TH1D *hEtaAkCs4[ ptDijetLow.size() ];
    TH1D *hEtaRat[ ptDijetLow.size() ];

    TLine *line;
    TLegend *leg;

    Int_t sizeX{1200};
    Int_t sizeY{800};

    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    Int_t ptBins = ptDijetLow.size();

    TCanvas *cComp = new TCanvas("cComp", "cComp", sizeX, sizeY);
    cComp->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    TCanvas *cRat = new TCanvas("cRat", "cRat", sizeX, sizeY);
    cRat->Divide(5, ( (ptBins % 5) == 0 ) ? (ptBins / 5) : (ptBins / 5 + 1) );

    // Loop over dijet pT bins
    for (Int_t i{0}; i<ptBins; i++) {

        std::cout << "pT bin: " << i << std::endl;

        // ak4
        hEtaAk4[i] = projectEtaFrom3D(hPtEtaDphiAk4, Form("hEtaAk4_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
        rescaleEta( hEtaAk4[i] );
        set1DStyle( hEtaAk4[i], ak4Type );

        // akCs4
        hEtaAkCs4[i] = projectEtaFrom3D(hPtEtaDphiAkCs4, Form("hEtaAkCs4_%d", i), ptDijetLow.at(i), ptDijetHi.at(i) );
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
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        leg = new TLegend(0.2, 0.75, 0.4, 0.85);
        leg->SetTextSize(0.04);
        leg->SetLineWidth(0);
        leg->AddEntry(hEtaAk4[i], Form("ak4"), "p");
        leg->AddEntry(hEtaAkCs4[i], Form("akCs4"), "p");
        leg->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_jetCollComp_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot ratio
        canv->cd();
        setPadStyle();
        hEtaRat[i]->Draw();
        hEtaRat[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaRat[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRat[i]->GetXaxis()->GetBinUpEdge(hEtaRat[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
        canv->SaveAs( Form("%s/pPb8160_etaDijet_jetCollRat_pt_%d_%d_%s.pdf", date.Data(), 
                      ptLow + (ptDijetLow.at(i)-1) * ptStep, 
                      ptLow + ptDijetHi.at(i) * ptStep, frame.Data() ) );

        // Plot all comparisons on one canvas
        cComp->cd(i+1);
        setPadStyle();
        hEtaAk4[i]->Draw();
        hEtaAkCs4[i]->Draw("same");
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
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
        t.DrawLatexNDC(0.25, 0.93, Form("%d < p_{T}^{dijet} (GeV/c) < %d", 
                       ptLow + (ptDijetLow.at(i) - 1) * ptStep, ptLow + ptDijetHi.at(i) * ptStep) );
        t.DrawLatexNDC( 0.65, 0.8, Form("%s frame", frame.Data() ) );
        line = new TLine(hEtaRat[i]->GetXaxis()->GetBinLowEdge(1), 1., 
                         hEtaRat[i]->GetXaxis()->GetBinUpEdge(hEtaRat[i]->GetNbinsX()), 1.);
        line->SetLineColor(kMagenta);
        line->SetLineWidth(3);
        line->SetLineStyle(3);
        line->Draw();
    } // for (Int_t i{0}; i<ptBins; i++)

    cComp->SaveAs( Form("%s/pPb8160_etaDijet_jetCollComp_all_%s.pdf", date.Data(), frame.Data() ) );
    cRat->SaveAs( Form("%s/pPb8160_etaDijet_jetCollRat_all_%s.pdf", date.Data(), frame.Data() ) );
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

}

//________________
void systematics() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Date
    TDatime dt;
    TString date { Form( "%d",dt.GetDate() ) };

    // List of file names
    TString defaultFileName("../build/MB_pPb8160_ak4.root");
    TString pbGoingFileName("../build/MB_Pbgoing_pPb8160_ak4.root");
    TString pGoingFileName("../build/MB_pgoing_pPb8160_ak4.root");

    TString pbGoingEmbeddingFileName("../build/oEmbedding_pPb8160_Pbgoing_jerDef_ak4.root");
    TString pGoingEmbeddingFileName("../build/oEmbedding_pPb8160_pgoing_jerDef_ak4.root");

    TString akcs4FileName("../build/MB_pPb8160_akCs4.root");
    TString jeuUpFileName("../build/MB_pPb8160_jeu_up_ak4.root");
    TString jeuDownFileName("../build/MB_pPb8160_jeu_down_ak4.root");
    TString embeddingFileName("../build/oEmbedding_pPb8160_ak4.root");

    TString jerDefFileName("../build/oEmbedding_pPb8160_jerDef_ak4.root");
    TString jerUpFileName("../build/oEmbedding_pPb8160_jerUp_ak4.root");
    TString jerDownFileName("../build/oEmbedding_pPb8160_jerDown_ak4.root");

    // Check the directory for storing figures exists
    if ( directoryExists( date.Data() ) ) {
        //std::cout << "Directory exists." << std::endl;
    } 
    else {
        createDirectory( date.Data() );
    }

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
    TFile *pbGoingEmbeddingFile = TFile::Open( pbGoingEmbeddingFileName );
    TFile *pGoingEmbeddingFile = TFile::Open( pGoingEmbeddingFileName );

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

    // Retrieve the base branch name to work on
    Int_t branchId{0};
    if ( defaultFileName.Contains("akCs4") ) {
        branchId = {0};
    }
    else {
        branchId = {1};
    }

    // plotDifferentDirections( pbGoingFile, pGoingFile, date );
    // plotDifferentDirections( pbGoingEmbeddingFile, pGoingEmbeddingFile, date );

    // compareData2McDifferentDirections(pbGoingFile, pGoingFile, pbGoingEmbeddingFile, pGoingEmbeddingFile, date);

    plotJEU( defaultFile, jeuUpFile, jeuDownFile, date );

    // plotJER(jerDefFile, jerUpFile, jerDownFile, date);

    // plotPointingResolution( embeddingFile, date );

    // compareJetCollections( defaultFile, akcs4File, date );

    // compareJER( embeddingFile, jerDefFile, jerUpFile, jerDownFile, date );

}