#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLine.h"
#include "TSystem.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TLatex.h"

#include <vector>

//________________
void plotCMSHeader() {
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.55, 0.93, "pp #sqrt{s_{NN}} = 5.02 TeV");
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
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);    
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
void plotJESvsPtHat(TFile *f) {

    double xTextPosition = 0.7;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);

    // Retrieve JES
    THnSparseD *hJESPars = dynamic_cast<THnSparseD*>( f->Get("hInclusiveJetJESGenPtGenEtaPtHatWeighted") );
    if ( !hJESPars ) {
        std::cerr << "[ERROR] Could not retrieve THnSparseD from the file" << std::endl;
        return;
    }

    // Create vector of ptHat and jet pT bins for projecitons
    int ptHatStart = 0;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins { 1, 4, 5, 6, 7 };

    int etaBinsProj[2] = {20, 32};   // eta range from (-1.2, 1.2):  52 bins from -5.2, 5.2

    TH2D *hJES[ ptHatBins.size() ];
    TH1D *hJESMean[ ptHatBins.size() ];
    TH1D *hJESSigma[ ptHatBins.size() ];;
    TCanvas *cJES[ ptHatBins.size() ];

    // Loop over ptHat bins
    for (unsigned int i{0}; i<ptHatBins.size(); i++) {

        //
        // JES vs. jet pT
        //

        // Set axis limits for THnSparse
        hJESPars->GetAxis(2)->SetRange( etaBinsProj[0], etaBinsProj[1] );
        hJESPars->GetAxis(3)->SetRange( ptHatBins[i], ptHatBinsMax );

        // Create 2D histogram
        hJES[i] = dynamic_cast<TH2D*>( hJESPars->Projection(1, 0) );
        hJES[i]->SetName( Form("hJES_%d", i) );
        set2DStyle(hJES[i]);

        // Retrieve JES and JER
        hJES[i]->FitSlicesY();

        // Jet energy scale
        hJESMean[i] = (TH1D*)gDirectory->Get( Form("hJES_%d_1", i) );
        hJESMean[i]->SetName( Form("hJESMean_%d", i) );
        set1DStyle(hJESMean[i], 2);
        hJESMean[i]->GetYaxis()->SetTitle("JES");

        // Jet energy resolution
        hJESSigma[i] = (TH1D*)gDirectory->Get( Form("hJES_%d_2", i) );
        hJESSigma[i]->SetName( Form("hJESSigma_%d", i) );
        set1DStyle(hJESSigma[i], 2);
        hJESSigma[i]->GetYaxis()->SetTitle("JER");

        // Create canvas
        cJES[i] = new TCanvas( Form("cJES_%d", i), Form("cJES_%d", i), 1500, 500 );
        cJES[i]->Divide(3, 1);
        
        cJES[i]->cd(1);
        setPadStyle();
        hJES[i]->Draw("colz");
        hJES[i]->GetXaxis()->SetRangeUser(0, 600);
        hJES[i]->GetYaxis()->SetRangeUser(0., 2.);
        plotCMSHeader();
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("#hat{p}_{T} > %d GeV",  ptHatStart + (ptHatBins.at(i) - 1) * ptHatStep) );

        cJES[i]->cd(2);
        setPadStyle();
        hJESMean[i]->Draw();
        hJESMean[i]->GetXaxis()->SetRangeUser(0, 600);
        hJESMean[i]->GetYaxis()->SetRangeUser(0.98, 1.3);
        plotCMSHeader();
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%2.1f < #eta^{jet} < %2.1f", -5.2 + etaBinsProj[0] * 0.2, -5.2 + etaBinsProj[1] * 0.2) );

        cJES[i]->cd(3);
        setPadStyle();
        hJESSigma[i]->Draw();
        hJESSigma[i]->GetXaxis()->SetRangeUser(0, 600);
        hJESSigma[i]->GetYaxis()->SetRangeUser(0.005, 0.17);
    } // for (unsigned int i{0}; i<ptHatBins.size(); i++)

}

//________________
void plotSimpleJES(TFile *f) {
    // The function reads the JES (reco pt/gen pt vs. gen pt for |eta|<1.4) from the file

    double xTextPosition = 0.7;
    double yTextPosition = 0.8;
    TLatex t;

    // Retrieve JES
    TH2D *h2D = dynamic_cast<TH2D*>( f->Get("hInclusiveJetJESVsPtGen") );
    set2DStyle(h2D);
    h2D->FitSlicesY();
    // Jet energy scale
    TH1D *hJESMean = (TH1D*)gDirectory->Get("hInclusiveJetJESVsPtGen_1");
    hJESMean->SetName("hJESMean");
    set1DStyle(hJESMean, 2);
    hJESMean->GetYaxis()->SetTitle("JES");

    // Jet energy resolution
    TH1D *hJESSigma = (TH1D*)gDirectory->Get("hInclusiveJetJESVsPtGen_2");
    hJESSigma->SetName("hJESSigma");
    set1DStyle(hJESSigma, 2);
    hJESSigma->GetYaxis()->SetTitle("JER");

    //
    // Create canvas
    //
    TCanvas *cSimpleJES = new TCanvas( "cSimpleJES", "cSimpleJES", 1500, 500 );
    cSimpleJES->Divide(3, 1);

    cSimpleJES->cd(1);
    setPadStyle();
    h2D->Draw("colz");
    //gPad->SetLogz();
    h2D->GetXaxis()->SetRangeUser(0, 800);

    plotCMSHeader();
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");

    cSimpleJES->cd(2);
    setPadStyle();
    hJESMean->Draw();
    plotCMSHeader();
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");
    hJESMean->GetXaxis()->SetRangeUser(0, 800);
    hJESMean->GetYaxis()->SetRangeUser(0.98, 1.3);

    cSimpleJES->cd(3);
    setPadStyle();
    hJESSigma->Draw();
    plotCMSHeader();
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");
    hJESSigma->GetXaxis()->SetRangeUser(0, 800);
    hJESSigma->GetYaxis()->SetRangeUser(-0.005, 0.17);
}

//________________
void plotJES() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    // pp5020 PYTHIA
    // TFile *f = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", uname.Data()) );
    // if ( !f ) {
    //     std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", uname.Data()) << std::endl;
    //     return;
    // }

    // pPb8160 embedding
    TFile *f = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", uname.Data()) );
    if ( !f ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", uname.Data()) << std::endl;
        return;
    }

    // Plot simple JES
    // plotSimpleJES(f);

    // Plot JES for different ptHat selections
    plotJESvsPtHat(f);
}