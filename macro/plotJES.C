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
#include "THnSparse.h"

#include <vector>

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
void plotCMSHeader(int collSystem = 1, double energy = 8.16) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Simulation}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, Form("%s #sqrt{s_{NN}} = %3.2f TeV", collSystemStr.Data(), energy) );
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
void plotJESandJER(TCanvas *c, TH2D *h2D, 
                   double lowVal = -1.6, double hiVal = 1.6,
                   int collSystem = 1, double energy = 8.16,
                   bool isPt = true, bool performFit = false) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.35;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);

    // Set xAxis range
    double xRange[2] = {10., 800.};
    double yMuRange[2] = {0.95, 1.05};
    double ySigmaRange[2] = {0., 0.25};
    
    if ( !isPt ) {
        xRange[0] = -3.5;
        xRange[1] = 3.5;
    }
    else {
        h2D->RebinX(2);
    }

    // Retrieve JES and JER
    h2D->FitSlicesY();

    // Jet energy scale
    TH1D *hJES = (TH1D*)gDirectory->Get( Form("%s_1", h2D->GetName()) );
    set1DStyle(hJES, 2);
    hJES->GetYaxis()->SetTitle("JES");

    // Jet energy resolution
    TH1D *hJER = (TH1D*)gDirectory->Get( Form("%s_2", h2D->GetName()) );
    set1DStyle(hJER, 2);
    hJER->GetYaxis()->SetTitle("JER");

    // Plot 2D distribution
    c->cd(1);
    setPadStyle();
    h2D->Draw("colz");
    h2D->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h2D->GetYaxis()->SetRangeUser(0., 2.);
    gPad->SetLogz();
    plotCMSHeader(collSystem, energy);

    // Plot JES
    c->cd(2);
    setPadStyle();
    hJES->Draw();
    hJES->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hJES->GetYaxis()->SetRangeUser(yMuRange[0], yMuRange[1]);
    if ( isPt ) {
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%2.1f < #eta^{jet} < %2.1f ", lowVal, hiVal) );
    }
    else {
        t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%4.0f < p_{T}^{jet} (GeV) < %4.0f ", lowVal, hiVal) );
    }
    plotCMSHeader(collSystem, energy);

    // Plot JER
    c->cd(3);
    setPadStyle();
    hJER->Draw();
    hJER->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hJER->GetYaxis()->SetRangeUser(ySigmaRange[0], ySigmaRange[1]);
    plotCMSHeader(collSystem, energy);

    TF1 *jerFit;
    if ( performFit ) {
        jerFit = new TF1( Form("%s_fit", hJER->GetName()), "sqrt([0] + [1] / x)", 30., 600.);
        jerFit->SetParameters(0.002, 1.0);
        jerFit->SetLineColor(kRed);
        jerFit->SetLineWidth(2);
        hJER->Fit(jerFit, "MRE");
        jerFit->Draw("same");

        t.DrawLatexNDC(0.35, 0.8, Form("Fit: #sqrt{a+b/x}")); 
        t.DrawLatexNDC(0.35, 0.7, Form("a = %5.4f, b = %4.3f", jerFit->GetParameter(0), jerFit->GetParameter(1)) );
    }
}

//________________
void plotJESandJER(TFile *f, int collSystem = 0, double energy = 5.02, TString date = "20250129") {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    double xTextPosition = 0.5;
    double yTextPosition = 0.8;
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);

    // Retrieve JES hInclusiveJetJESGenPtGenEtaPtHatWeighted hInclusiveJetJESRecoPtRecoEtaPtHatWeighted
    THnSparseD *hJESPars = (THnSparseD*)f->Get("hInclusiveJetJESGenPtGenEtaPtHatWeighted");
    hJESPars->SetName("hJESPars");
    if ( hJESPars == nullptr ) {
        std::cerr << Form( "[ERROR] Could not retrieve hInclusiveJetJESGenPtGenEtaPtHatWeighted from the file: %s", f->GetName() ) << std::endl;
        return;
    }

    // Create vector of ptHat and jet pT bins for projecitons
    int ptHatStart = 0;
    int ptHatStep = 10; // Starting from 10 GeV: ptHatStart + (ptHatBins(i) - 1) * ptHatStep
    int ptHatBinsMax = 100;
    std::vector<int> ptHatBins { 3 };

    // Event ptHat binning
    double ptHatVals[] = { 15. };
    int sizeOfPtHatBins = sizeof(ptHatVals) / sizeof(ptHatVals[0]);

    // Jet eta binning
    double etaVals[] = { -5.2, -3.6, -3.0, -2.8, -2.6, -2.4, -2.0, -1.6, 0., 1.6, 2.0, 2.4, 2.6, 2.8, 3.0, 3.6, 5.2 };
    int sizeOfEtaBins = sizeof(etaVals) / sizeof(etaVals[0]);

    // Jet pT binning
    double ptVals[] = { 25., 55., 95., 125., 200., 1000. };
    int sizeOfPtBins = sizeof(ptVals) / sizeof(ptVals[0]);

    TH2D *h2D {nullptr};
    TCanvas *c = new TCanvas( "c", "cJESandJER", 1500, 500 );
    c->Divide(3, 1);

    // Loop over ptHat bins
    for (unsigned int i{0}; i<sizeOfPtHatBins; i++) {

        double ptHatLow = ptHatVals[i];

        // Set ptHat binning
        hJESPars->GetAxis(3)->SetRange( hJESPars->GetAxis(3)->GetBinLowEdge(ptHatLow), 
                                        hJESPars->GetAxis(3)->GetNbins() );

        //
        // JES vs. jet pT for different eta bins
        //
        for (unsigned int j{0}; j<(sizeOfEtaBins-1); j++) {

            double etaLow = etaVals[j];
            double etaHi = etaVals[j+1];

            // Set axis limits for THnSparse
            hJESPars->GetAxis(2)->SetRange( hJESPars->GetAxis(2)->GetBinLowEdge( etaLow ), 
                                            hJESPars->GetAxis(2)->GetBinUpEdge( etaHi ) - 1 );

            // Create 2D histogram
            h2D = dynamic_cast<TH2D*>( hJESPars->Projection(0, 1) );
            h2D->SetName( Form("hJES_pt_%d_%d", i, j) );
            set2DStyle(h2D);


            plotJESandJER(c, h2D, etaLow, etaHi, collSystem, energy, true, true);

            c->SaveAs( Form("%s/JES_vs_pt_%d_%d.pdf", date.Data(), i, j) );

        } // for (unsigned int j{0}; j<(sizeOfEtaBins-1); j++)

        hJESPars->GetAxis(2)->SetRange( 1, hJESPars->GetAxis(2)->GetNbins() ); // Restore eta binning

        //
        // JES vs. jet eta for different pt bins
        //
        for (unsigned int j{0}; j<(sizeOfPtBins-1); j++) {

            double ptLow = ptVals[j];
            double ptHi = ptVals[j+1];

            // Set axis limits for THnSparse
            hJESPars->GetAxis(1)->SetRange( hJESPars->GetAxis(1)->GetBinLowEdge( ptLow ), 
                                            hJESPars->GetAxis(1)->GetBinUpEdge( ptHi ) - 1 );

            // Create 2D histogram
            h2D = dynamic_cast<TH2D*>( hJESPars->Projection(0, 2) );
            h2D->SetName( Form("hJES_eta_%d_%d", i, j) );
            set2DStyle(h2D);

            plotJESandJER(c, h2D, lowVal, hiVal, collSystem, energy, false, false);
            c->SaveAs( Form("%s/JES_vs_eta_%d_%d.pdf", date.Data(), i, j) );
        } // for (unsigned int j{0}; j<(sizeOfPtBins-1); j++)

    } // for (unsigned int i{0}; i<sizeOfPtHatBins; i++)

}

//________________
void plotSimpleJES(TFile *f, int collSystem = 0, double energy = 5.02) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

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
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");

    cSimpleJES->cd(2);
    setPadStyle();
    hJESMean->Draw();
    plotCMSHeader(collSystem, energy);
    t.DrawLatexNDC(xTextPosition, yTextPosition, "|#eta|<1.4");
    hJESMean->GetXaxis()->SetRangeUser(0, 800);
    hJESMean->GetYaxis()->SetRangeUser(0.98, 1.3);

    cSimpleJES->cd(3);
    setPadStyle();
    hJESSigma->Draw();
    plotCMSHeader(collSystem, energy);
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
    TString username = gSystem->GetFromPipe("whoami");

    // Date extraction
    TDatime dt; 
    TString date { Form( "%d",dt.GetDate() ) };

    if ( !directoryExists( date.Data() ) ) {
        createDirectory( date.Data() );
    } 

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    int collSystem = 0;
    double energy = 5.02;

    // // pp5020 PYTHIA
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", username.Data()) );
    // if ( !inputFile ) {
    //     std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC_3020.root", username.Data()) << std::endl;
    //     return;
    // }
    // collSystem = 0;
    // energy = 5.02;

    // pPb8160 embedding
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", username.Data()) );
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pgoing_def_ak4_eta25.root", username.Data()) );
    TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_eta20.root", username.Data()) );
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oPythia_Pbgoing_def_ak4_eta20.root", username.Data()) );
    // TFile *inputFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/build/oTest_pPb8160_dijet_ptHat_50_80_noTrkMax.root", username.Data()) );
    if ( !inputFile || inputFile->IsZombie() ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oPythia_Pbgoing_def_ak4_eta20.root", username.Data()) << std::endl;
        return;
    }
    collSystem = 1;
    energy = 8.16;

    // Plot simple JES
    // plotSimpleJES( inputFile, collSystem, energy );

    // Plot JES for different ptHat selections
    plotJESandJER( inputFile, collSystem, energy, date );
}
