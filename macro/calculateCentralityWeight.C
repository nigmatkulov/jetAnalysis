// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"

#include <iostream>

using namespace std;

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.15);
}

//________________
Double_t fitFunction(Double_t *x, Double_t *par) {
    return par[0] + par[1] * (*x) + TMath::Power(*x, 2) * par[2] + TMath::Power(*x, 3) * par[3] +
           TMath::Power(*x, 4) * par[4] + TMath::Power(*x, 5) * par[5];
}

//________________
void calculateCentralityWeight(const Char_t* expFName = "../build/oTestReadForest_PbPb_exp.root",
                               const Char_t* mcFName = "../build/oTestReadForest_PbPb_mcHiBinShift.root") {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Read ROOT file                             
    TFile *expFile = TFile::Open(expFName);
    if ( !expFile->IsOpen() ) {
        std::cout << "Input file for data not opened. Terminating." << std::endl;
        exit(0);
    }

    // Read ROOT file                             
    TFile *mcFile = TFile::Open(mcFName);
    if ( !mcFile->IsOpen() ) {
        std::cout << "Input file for MC not opened. Terminating." << std::endl;
        exit(0);
    }

    TH1D* hExpHiBin = (TH1D*)expFile->Get("hHiBin");
    hExpHiBin->SetName("hExpHiBin");
    hExpHiBin->SetLineWidth(2);
    hExpHiBin->SetLineColor(kRed);
    hExpHiBin->SetMarkerColor(kRed);
    hExpHiBin->SetMarkerStyle(20);
    hExpHiBin->Scale( 1. / hExpHiBin->Integral() );

    TH1D* hMcHiBin = (TH1D*)mcFile->Get("hHiBin");
    hMcHiBin->SetName("hMcHiBin");
    hMcHiBin->SetLineWidth(2);
    hMcHiBin->SetLineColor(kBlue);
    hMcHiBin->SetMarkerColor(kBlue);
    hMcHiBin->SetMarkerStyle(24);
    hMcHiBin->Scale( 1. / hMcHiBin->Integral() );

    TH1D *hCentralityWeight = new TH1D( "hCentralityWeight",
                                        "hCentralityWeight;HiBin;Data/MC",
                                        hExpHiBin->GetNbinsX(),
                                        hExpHiBin->GetXaxis()->GetBinLowEdge(1),
                                        hExpHiBin->GetXaxis()->GetBinUpEdge( hExpHiBin->GetNbinsX() ) );
    hCentralityWeight->Sumw2();
    hCentralityWeight->SetLineWidth(2);
    hCentralityWeight->SetLineColor(kBlack);
    hCentralityWeight->SetMarkerColor(kBlack);
    hCentralityWeight->SetMarkerStyle(20);

    hCentralityWeight->Divide( hExpHiBin, hMcHiBin );

    TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
    canv->Divide(2, 1);

    // Compare data to MC
    canv->cd(1);
    setPadStyle();
    hExpHiBin->Draw();
    hMcHiBin->Draw("same");
    auto leg = new TLegend(0.45, 0.75, 0.9, 0.9);
    leg->AddEntry(hExpHiBin,"Experiment","p");
    leg->AddEntry(hMcHiBin,"PYTHIA+HYDJET","p");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->Draw();

    // Ratio of data over MC
    canv->cd(2);
    setPadStyle();
    hCentralityWeight->Draw();
    hCentralityWeight->GetYaxis()->SetRangeUser( -0.2, 6.2 );

    TF1 *f = new TF1("f", fitFunction, -0.5, 181., 6);
    f->SetParameters( 4.67807, -0.131387, 0.00227593, 0., 0., 0.);
    //f->FixParameter(3, 0.);
    //f->FixParameter(4, 0.);
    f->FixParameter(5, 0.);
    f->SetLineWidth(3);
    hCentralityWeight->Fit("f", "MRE");
    std::cout << "Fit function: " << f->GetExpFormula() << std::endl;

    // Extracting fit parameters to an array
    const int numParameters = f->GetNpar();
    double fitParameters[numParameters];

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    t.DrawLatexNDC( 0.2, 0.9, Form("[0]+[1]*x+[2]*x^{2}+[3]*x^{3}+[4]*x^{5}+[5]*x^{6}") );
    for (int i = 0; i < numParameters; ++i) {
        fitParameters[i] = f->GetParameter(i);
        std::cout << Form("par[%d] = %e", i, fitParameters[i] ) << std::endl;
        t.DrawLatexNDC( 0.5, 0.85 - 0.05*i, Form("[%d]=%e", i, fitParameters[i]) );
    }
}