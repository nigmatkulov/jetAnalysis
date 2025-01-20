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

// C++ headers
#include <iostream>
#include <vector>

//________________
void plotCMSHeader() {
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.65, 0.93, "#sqrt{s_{NN}} = 5.02 TeV");
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
                    int ptLow=55, int ptHi=75,
                    const char* h1Name="h1", const char* h2Name="h2",
                    bool isRpPb = false) {

    // Text 
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.06);

    // Number for plotting position
    Double_t xRange[2] = {-3., 3.};
    Double_t yRange[2] = {0., 0.14};
    if (isRpPb) {
        yRange[0] = 0.8; yRange[1] = 1.2;
    }
    Double_t legX[2] = {0.4, 0.65};
    Double_t legY[2] = {0.2, 0.35};

    // Make ratios
    TH1D *hRatio = dynamic_cast<TH1D*>( h2->Clone("hRatio") );
    hRatio->Divide( h1 );

    // Create pad
    TLegend *leg;
    TLine *line;

    //
    // Plot comparison
    //

    c->cd(1);

    setPadStyle();
    // Plot distributions
    h1->Draw();
    h2->Draw("same");
    // Set ranges
    h1->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    h1->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
    h1->GetYaxis()->SetTitle("dN/d#eta^{dijet}");
    h1->GetXaxis()->SetTitle("#eta^{dijet}");

    t.DrawLatexNDC(0.3, 0.85, Form("%d < p_{T}^{ave} (GeV) < %d", ptLow, ptHi));
    plotCMSHeader();
    

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    leg->AddEntry(h1, Form("%s", h1Name), "p");
    leg->AddEntry(h2, Form("%s", h2Name), "p");
    leg->Draw();

    //
    // Plot ratio
    //

    c->cd(2);
    setPadStyle();
    hRatio->Draw();
    hRatio->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetYaxis()->SetTitle( Form( "%s / %s", h2Name, h1Name ) );
    hRatio->GetXaxis()->SetTitle("#eta^{dijet}");

    // Legend
    leg = new TLegend(legX[0], legY[0], legX[1], legY[1]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont( 42 );
    //leg->AddEntry(hRatio, Form( "%s / %s", h2Name, h1Name ), "p");
    leg->Draw();

    // Line at unity
    line = new TLine(xRange[0], 1.0, xRange[1], 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kMagenta);
    line->Draw();
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
void compare_pp5020(TFile *pubFile, TFile *dataFile, TFile *pythiaFile) {

    // Dijet ptAve binning
    int dijetPtOldVals[7] {25, 55, 75, 95, 115, 150, 400};
    std::vector<int> dijetPtVals(dijetPtOldVals, dijetPtOldVals + sizeof(dijetPtOldVals) / sizeof(int));

    // Histograms 
    TH1D *hPubDijetEta[6];
    TH1D *hDataDijetEta[6];
    TH1D *hPythiaRecoDijetEta[6];
    TH1D *hPythiaGenDijetEta[6];

    // Loop over dijet pT bins
    for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++) {
        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];

        // Initialize histograms with nullptr first
        hPubDijetEta[i] = nullptr;
        hDataDijetEta[i] = nullptr;
        hPythiaRecoDijetEta[i] = nullptr;
        hPythiaGenDijetEta[i] = nullptr;

        // Published data
        if ( i != 0 ) {
            hPubDijetEta[i] = dynamic_cast<TH1D*>( pubFile->Get( Form("ppEta_pt_%d_%d", ptLow, ptHi) ) );
            hPubDijetEta[i]->SetName( Form("hPubDijetEta_%d", i) );
            set1DStyle( hPubDijetEta[i], 2 );
        }

        hDataDijetEta[i] = dynamic_cast<TH1D*>( dataFile->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hDataDijetEta[i]->SetName( Form("hDataDijetEta_%d", i) );
        rescaleEta( hDataDijetEta[i] );
        set1DStyle( hDataDijetEta[i], 0 );

        hPythiaRecoDijetEta[i] = dynamic_cast<TH1D*>( pythiaFile->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hPythiaRecoDijetEta[i]->SetName( Form("hPythiaRecoDijetEta_%d", i) );
        rescaleEta( hPythiaRecoDijetEta[i] );
        set1DStyle( hPythiaRecoDijetEta[i], 1 );

        hPythiaGenDijetEta[i] = dynamic_cast<TH1D*>( pythiaFile->Get( Form("hGenDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hPythiaGenDijetEta[i]->SetName( Form("hPythiaGenDijetEta_%d", i) );
        rescaleEta( hPythiaGenDijetEta[i] );
        set1DStyle( hPythiaGenDijetEta[i], 3 );
    } // for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++)

    //
    // Plot comparisons
    //

    TCanvas *cPub2DataComparison[5];
    TCanvas *cData2PythiaRecoComparison[5];
    TCanvas *cPythiaGen2PythiaRecoComparison[5];
    TCanvas *cPub2PythiaGenComparison[5];

    // Loop over dijet ptAve bins
    for (int i{1}; i < dijetPtVals.size() - 1; i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{500}, canvY{1000};

        //
        // Published vs. data
        //

        // cPub2DataComparison[i-1] = new TCanvas( Form("ppPub2DataComparison_%d", i-1), 
        //                                         Form("ppPub2DataComparison_%d", i-1), 
        //                                          canvX, canvY );
        // cPub2DataComparison[i-1]->Divide(1, 2);
        // plotComparison(cPub2DataComparison[i-1], hPubDijetEta[i], hDataDijetEta[i], 
        //                ptLow, ptHi, "pp5020 pub.", "pp5020 my");

        //
        // Data vs. pythia reco
        //

        // cData2PythiaRecoComparison[i-1] = new TCanvas( Form("ppData2PythiaRecoComparison_%d", i-1), 
        //                                                Form("ppData2PythiaRecoComparison_%d", i-1), 
        //                                                canvX, canvY );
        // cData2PythiaRecoComparison[i-1]->Divide(1, 2);
        // plotComparison( cData2PythiaRecoComparison[i-1], hDataDijetEta[i], hPythiaRecoDijetEta[i], 
        //                ptLow, ptHi, "pp5020 my", "pp5020 pythia reco");

        //
        // Pythia gen vs. pythia reco
        //

        cPythiaGen2PythiaRecoComparison[i-1] = new TCanvas( Form("ppPythiaGen2PythiaRecoComparison_%d", i-1), 
                                                            Form("ppPythiaGen2PythiaRecoComparison_%d", i-1), 
                                                            canvX, canvY );
        cPythiaGen2PythiaRecoComparison[i-1]->Divide(1, 2);
        plotComparison( cPythiaGen2PythiaRecoComparison[i-1], hPythiaGenDijetEta[i], hPythiaRecoDijetEta[i], 
                       ptLow, ptHi, "pp5020 pythia gen", "pp5020 pythia reco");

        
        //
        // Published vs. pythia gen
        //

        // cPub2PythiaGenComparison[i-1] = new TCanvas( Form("ppPub2PythiaGenComparison_%d", i-1), 
        //                                              Form("ppPub2PythiaGenComparison_%d", i-1), 
        //                                              canvX, canvY );
        // cPub2PythiaGenComparison[i-1]->Divide(1, 2);
        // plotComparison( cPub2PythiaGenComparison[i-1], hPubDijetEta[i], hPythiaGenDijetEta[i], 
        //                ptLow, ptHi, "pp5020 pub.", "pp5020 pythia gen");

    } // for (int i{1}; i < dijetPtVals.size() - 1; i++)

}

//________________
void compare_pPb5020(TFile *pub, TFile *runB, TFile *runD) {

    // Dijet ptAve binning
    int dijetPtOldVals[7] {25, 55, 75, 95, 115, 150, 400};
    std::vector<int> dijetPtVals(dijetPtOldVals, dijetPtOldVals + sizeof(dijetPtOldVals) / sizeof(int));

    // Histograms 
    TH1D *hPubDijetEta[6];
    TH1D *hRunBDijetEta[6];
    TH1D *hRunDDijetEta[6];

    // Loop over dijet ptAve bins
    for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++) {
        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];

        // Initialize histograms with nullptr first
        hPubDijetEta[i] = nullptr;
        hRunBDijetEta[i] = nullptr;
        hRunDDijetEta[i] = nullptr;

        // Published data
        if ( i != 0 ) {
            hPubDijetEta[i] = dynamic_cast<TH1D*>( pub->Get( Form("pPbEta_pt_%d_%d", ptLow, ptHi) ) );
            hPubDijetEta[i]->SetName( Form("hPubDijetEta_%d", i) );
            set1DStyle( hPubDijetEta[i], 2 );
        }

        hRunBDijetEta[i] = dynamic_cast<TH1D*>( runB->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hRunBDijetEta[i]->SetName( Form("hRunBDijetEta_%d", i) );
        rescaleEta( hRunBDijetEta[i] );
        set1DStyle( hRunBDijetEta[i], 0 );

        hRunDDijetEta[i] = dynamic_cast<TH1D*>( runD->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hRunDDijetEta[i]->SetName( Form("hRunDDijetEta_%d", i) );
        rescaleEta( hRunDDijetEta[i] );
        set1DStyle( hRunDDijetEta[i], 1 );
    } // for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++)

    //
    // Plot comparisons
    //

    TCanvas *cPub2RunBComparison[5];
    TCanvas *cPub2RunDComparison[5];
    TCanvas *cRunB2RunDComparison[5];
    
    // Loop over dijet ptAve bins
    for (int i{1}; i < dijetPtVals.size() - 1; i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{500}, canvY{1000};

        //
        // Published vs. RunB
        //

        // cPub2RunBComparison[i-1] = new TCanvas( Form("pPbPub2RunBComparison_%d", i-1), 
        //                                         Form("pPbPub2RunBComparison_%d", i-1), 
        //                                         canvX, canvY );
        // cPub2RunBComparison[i-1]->Divide(1, 2);
        // plotComparison(cPub2RunBComparison[i-1], hPubDijetEta[i], hRunBDijetEta[i], 
        //                ptLow, ptHi, "pPb5020 pub.", "pPb5020 RunB");

        //
        // Published vs. RunD
        //

        // cPub2RunDComparison[i-1] = new TCanvas( Form("pPbPub2RunDComparison_%d", i-1), 
        //                                         Form("pPbPub2RunDComparison_%d", i-1), 
        //                                         canvX, canvY );
        // cPub2RunDComparison[i-1]->Divide(1, 2);
        // plotComparison(cPub2RunDComparison[i-1], hPubDijetEta[i], hRunDDijetEta[i], 
        //                ptLow, ptHi, "pPb5020 pub.", "pPb5020 RunD");

        //
        // RunB vs. RunD
        //

        // cRunB2RunDComparison[i-1] = new TCanvas( Form("pPbRunB2RunDComparison_%d", i-1), 
        //                                          Form("pPbRunB2RunDComparison_%d", i-1), 
        //                                          canvX, canvY );
        // cRunB2RunDComparison[i-1]->Divide(1, 2);
        // plotComparison(cRunB2RunDComparison[i-1], hRunBDijetEta[i], hRunDDijetEta[i], 
        //                ptLow, ptHi, "pPb5020 RunB", "pPb5020 RunD");

    } // for (int i{1}; i < dijetPtVals.size() - 1; i++)
}

//________________
void compare_RpPb(TFile *pub, TFile *pPb, TFile *pp) {

    // Dijet ptAve binning
    int dijetPtOldVals[7] {25, 55, 75, 95, 115, 150, 400};
    std::vector<int> dijetPtVals(dijetPtOldVals, dijetPtOldVals + sizeof(dijetPtOldVals) / sizeof(int));

    // Histograms
    TH1D *hPubRpPb[ dijetPtVals.size() - 1 ];
    TH1D *hRpPbMy[ dijetPtVals.size() - 1 ];
    TH1D *hPPbDijetEta[ dijetPtVals.size() - 1 ];
    TH1D *hPPDijetEta[ dijetPtVals.size() - 1 ];

    // Loop over dijet ptAve bins
    for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++) {
        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];

        // Initialize histograms with nullptr first
        hPubRpPb[i] = nullptr;
        hRpPbMy[i] = nullptr;
        hPPbDijetEta[i] = nullptr;
        hPPDijetEta[i] = nullptr;

        // Published RpPb
        if ( i != 0 ) {
            hPubRpPb[i] = dynamic_cast<TH1D*>( pub->Get( Form("RpPb_pt_%d_%d", ptLow, ptHi) ) );
            hPubRpPb[i]->SetName( Form("hPubRpPb_%d", i) );
            set1DStyle( hPubRpPb[i], 2 );
        }

        // pPb data (my)
        hPPbDijetEta[i] = dynamic_cast<TH1D*>( pPb->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hPPbDijetEta[i]->SetName( Form("hPPbDijetEta_%d", i) );
        rescaleEta( hPPbDijetEta[i] );
        set1DStyle( hPPbDijetEta[i], 0 );

        // pp data (my)
        hPPDijetEta[i] = dynamic_cast<TH1D*>( pp->Get( Form("hRecoDijetEta1DOldPtBinningWeighted_%d", i) ) );
        hPPDijetEta[i]->SetName( Form("hPPDijetEta_%d", i) );
        rescaleEta( hPPDijetEta[i] );
        set1DStyle( hPPDijetEta[i], 1 );

        // Create RpPb from pPb and pp
        hRpPbMy[i] = dynamic_cast<TH1D*>( hPPbDijetEta[i]->Clone( Form("hRpPbMy_%d", i) ) );
        hRpPbMy[i]->Divide( hPPDijetEta[i] );
    } // for (unsigned int i = 0; i < dijetPtVals.size() - 1; i++)


    //
    // Plot comparisons
    //

    TCanvas *cPub2MyComparison[5];

    for (int i{1}; i < dijetPtVals.size() - 1; i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{500}, canvY{1000};

        //
        // Published vs. my
        //

        cPub2MyComparison[i-1] = new TCanvas( Form("RpPbPub2MyComparison_%d", i-1), 
                                              Form("RpPbPub2MyComparison_%d", i-1), 
                                              canvX, canvY );
        cPub2MyComparison[i-1]->Divide(1, 2);
        plotComparison(cPub2MyComparison[i-1], hPubRpPb[i], hRpPbMy[i], 
                       ptLow, ptHi, "RpPb pub.", "RpPb my", true);

    } // for (int i{1}; i < dijetPtVals.size() - 1; i++)

}

//________________
void crossCheckProjections(TFile *f, const char *collisionSystem = "pp5020",int iCase = 0) {

    // iCase: 0 - reco, 1 - gen, 2 - ref
    if ( iCase < 0 || iCase > 2 ) {
        std::cerr << "[ERROR] Invalid case: " << iCase << std::endl;
        return;
    }
    TString jetType = (iCase == 0) ? "Reco" : (iCase == 1) ? "Gen" : "Ref";

    // Dijet ptAve binning
    int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
                             100, 110,  120, 130, 140,
                             150, 160,  180, 200, 250, 
                             300, 500};
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);

    std::vector<int> dijetPtVals; 
    dijetPtVals.assign(dijetPtNewVals, dijetPtNewVals + sizeOfPtVals);

    // Bins for projections from 3D
    std::vector<int> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55 };
    std::vector<int> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94 };

    // Retrieve 3D histogram
    TH3D *h3D = dynamic_cast<TH3D*>( f->Get( Form( "h%sDijetPtEtaDphiWeighted", jetType.Data() ) ) );
    TH1D *h1DProj[ ptDijetBinLow.size() ];
    TH1D *h1DDirect[ ptDijetBinLow.size() ];
    TCanvas *cComp[ ptDijetBinLow.size() ];

    // Loop over dijet ptAve bins
    for (unsigned int i = 0; i < ptDijetBinLow.size(); i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{500}, canvY{1000};

        // Make projection
        h1DProj[i] = dynamic_cast<TH1D*>( h3D->ProjectionY( Form("h%sDijetEta_%d_%d", jetType.Data(), dijetPtVals[i], dijetPtVals[i+1]), 
                                                            ptDijetBinLow[i], ptDijetBinHi[i] ) );
        h1DProj[i]->SetName( Form("h%sDijetEtaProj_%d_%d", jetType.Data(), dijetPtVals[i], dijetPtVals[i+1]) );
        rescaleEta( h1DProj[i] );
        set1DStyle( h1DProj[i], 0 );

        // Retrieve a directly filled 1D histogram
        h1DDirect[i] = dynamic_cast<TH1D*>( f->Get( Form("h%sDijetEta1DWeighted_%d", jetType.Data(), i) ) );
        h1DDirect[i]->SetName( Form("h%sDijetEtaDirect_%d", jetType.Data(), i) );
        rescaleEta( h1DDirect[i] );
        set1DStyle( h1DDirect[i], 1 );

        // Create canvas
        cComp[i] = new TCanvas( Form("c%sDijetEtaComp_%d", jetType.Data(), i), 
                                Form("c%sDijetEtaComp_%d", jetType.Data(), i), 
                                canvX, canvY );
        cComp[i]->Divide(1, 2);
        plotComparison( cComp[i], h1DProj[i], h1DDirect[i], 
                        ptLow, ptHi, Form("%s proj.", jetType.Data()), Form("%s direct", jetType.Data()) );

    }
}

//________________
void compareNew2Old(TFile *newFile, TFile *oldFile, int iCase = 0) {
    
    // iCase: 0 - reco, 1 - gen, 2 - ref
    if ( iCase < 0 || iCase > 2 ) {
        std::cerr << "[ERROR] Invalid case: " << iCase << std::endl;
        return;
    }
    TString jetType = (iCase == 0) ? "Reco" : (iCase == 1) ? "Gen" : "Ref";

    // Dijet ptAve binning
    int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
                             100, 110,  120, 130, 140,
                             150, 160,  180, 200, 250, 
                             300, 500};
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);

    std::vector<int> dijetPtVals; 
    dijetPtVals.assign(dijetPtNewVals, dijetPtNewVals + sizeOfPtVals);

    // Bins for projections from 3D
    std::vector<int> ptDijetBinLow {5, 7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 35, 45, 55 };
    std::vector<int> ptDijetBinHi  {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 34, 44, 54, 94 };

    // Retrieve 3D histogram
    TH3D *h3DNew = dynamic_cast<TH3D*>( newFile->Get( Form( "h%sDijetPtEtaDphiWeighted", jetType.Data() ) ) );
    TH3D *h3DOld = dynamic_cast<TH3D*>( oldFile->Get( Form( "h%sDijetPtEtaDphiWeighted", jetType.Data() ) ) );
    TH1D *h1DProjNew[ ptDijetBinLow.size() ];
    TH1D *h1DProjOld[ ptDijetBinLow.size() ];

    TCanvas *cComp[ ptDijetBinLow.size() ];

    // Loop over dijet ptAve bins
    for (unsigned int i = 0; i < ptDijetBinLow.size(); i++) {

        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i+1];
        int canvX{500}, canvY{1000};

        // Make projection
        h1DProjNew[i] = dynamic_cast<TH1D*>( h3DNew->ProjectionY( Form("h%sDijetEta_%d_%d", jetType.Data(), dijetPtVals[i], dijetPtVals[i+1]), 
                                                                 ptDijetBinLow[i], ptDijetBinHi[i] ) );
        h1DProjNew[i]->SetName( Form("h%sDijetEtaProjNew_%d_%d", jetType.Data(), dijetPtVals[i], dijetPtVals[i+1]) );
        rescaleEta( h1DProjNew[i] );
        set1DStyle( h1DProjNew[i], 0 );

        // Make projection
        h1DProjOld[i] = dynamic_cast<TH1D*>( h3DOld->ProjectionY( Form("h%sDijetEta_%d_%d", jetType.Data(), dijetPtVals[i], dijetPtVals[i+1]), 
                                                                 ptDijetBinLow[i], ptDijetBinHi[i] ) );
        h1DProjOld[i]->SetName( Form("h%sDijetEtaProjOld_%d_%d", jetType.Data(), dijetPtVals[i], dijetPtVals[i+1]) );
        rescaleEta( h1DProjOld[i] );
        set1DStyle( h1DProjOld[i], 1 );

        // Create canvas
        cComp[i] = new TCanvas( Form("c%sDijetEtaComp_%d", jetType.Data(), i), 
                                Form("c%sDijetEtaComp_%d", jetType.Data(), i), 
                                canvX, canvY );
        cComp[i]->Divide(1, 2);
        plotComparison( cComp[i], h1DProjNew[i], h1DProjOld[i], 
                        ptLow, ptHi, Form("%s proj. new", jetType.Data()), Form("%s proj. old", jetType.Data()) );

    }
}

//________________
void comparisons_5TeV() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);


    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    //
    // Published data
    //

    TFile *pubFile = TFile::Open("pPb5020/cms_dijet_eta_5TeV_pub.root");
    if ( !pubFile ) {
        std::cerr << "File not found: pPb5020/cms_dijet_eta_5TeV_pub.root" << std::endl;
        return;
    }

    //
    // pp5020
    //

    // Processed data
    TFile *pp5020DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/exp/pp5020_2017_woExtraJEC.root", uname.Data()) );
    // TFile *pp5020DataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/exp/pp5020_2017_woExtraJEC.root", uname.Data()) );
    if ( !pp5020DataFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/exp/pp5020_2017_wExtraJEC.root", uname.Data()) << std::endl;
        return;
    }

    // Pythia 
    // TFile *pp5020PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_wExtraJEC.root", uname.Data()) );
    TFile *pp5020PythiaFile = TFile::Open( Form("/Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC.root", uname.Data()) );
    if ( !pp5020PythiaFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pp5020/pythia/pp5020_pythia8_woExtraJEC.root", uname.Data()) << std::endl;
        return;
    }

    //
    // pPb5020
    //

    // RunB
    TFile *pPb5020RunBDataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb5020/exp/RunB/PAEGJet_RunB_pPb5020_ak4PF.root", uname.Data() ) );
    if ( !pPb5020RunBDataFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb5020/exp/RunB/PAEGJet_RunB_pPb5020_ak4PF.root", uname.Data()) << std::endl;
        return;
    }

    // RunD
    TFile *pPb5020RunDDataFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb5020/exp/RunD/PAEGJet_RunD_pPb5020_ak4PF.root", uname.Data() ) );
    if ( !pPb5020RunDDataFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb5020/exp/RunD/PAEGJet_RunD_pPb5020_ak4PF.root", uname.Data()) << std::endl;
        return;
    }

    //
    // pPb8160
    //

    // MC p-going direction new (coincides with the pPb5020)
    TFile *pPb8160EmbedNewFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pPb8160_pgoing_new.root", uname.Data()) );
    if ( !pPb8160EmbedNewFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pPb8160_pgoing_new.root", uname.Data()) << std::endl;
        return;
    }

    TFile *pPb8160EmbedOldFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/pgoing/jer/oEmbedding_pPb8160_pgoing_jerDef_ak4.root", uname.Data()) );
    if ( !pPb8160EmbedOldFile ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/pgoing/jer/oEmbedding_pPb8160_pgoing_jerDef_ak4.root", uname.Data()) << std::endl;
        return;
    }


    //
    // Compare distributions for pp5020
    //
    compare_pp5020(pubFile, pp5020DataFile, pp5020PythiaFile);

    // compareReco2GenInclusiveJetPt(pp5020PythiaFile);

    //
    // Compare distributions for pPb5020
    //
    // compare_pPb5020(pubFile, pPb5020RunBDataFile, pPb5020RunDDataFile);

    //
    // Compare distributions for RpPb
    //
    // compare_RpPb(pubFile, pPb5020RunBDataFile, pp5020DataFile);

    //
    // Cross check projections
    //
    // crossCheckProjections(pp5020PythiaFile, "pp5020", 0);  // Reco
    // crossCheckProjections(pp5020PythiaFile, "pp5020", 1);  // Gen
    // crossCheckProjections(pp5020PythiaFile, "pp5020", 2);  // Ref

    //
    // Compare new and old for pPb8160
    //
    // compareNew2Old(pPb8160EmbedNewFile, pPb8160EmbedOldFile, 0);  // Reco
    // compareNew2Old(pPb8160EmbedNewFile, pPb8160EmbedOldFile, 1);  // Gen
    // compareNew2Old(pPb8160EmbedNewFile, pPb8160EmbedOldFile, 2);  // Ref

}
