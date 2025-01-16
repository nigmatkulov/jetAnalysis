#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TLine.h"
#include "TSystem.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TLatex.h"

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
    gPad->SetLogz();
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

    // Pythia 
    // TFile *pp5020PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_wExtraJEC.root", uname.Data()) );
    TFile *pp5020PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC.root", uname.Data()) );
    if ( !pp5020PythiaFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC.root", uname.Data()) << std::endl;
        return;
    }

    // Plot simple JES
    plotSimpleJES(pp5020PythiaFile);
}