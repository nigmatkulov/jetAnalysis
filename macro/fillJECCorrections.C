#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "../JetCorrector.h"
#include "../JetCorrector.cc"

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
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
}

//________________
void fillJECCorrections() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);
    
    JetCorrector *jec = new JetCorrector("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt");
    //JetCorrector *jec = new JetCorrector("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt");

    // TFile *f = new TFile("pPb8160_JEC_L2Relative.root", "RECREATE");

    std::vector<double> ptValues { 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 180., 200., 250., 300., 350., 500., 700., 1000., 1500., 6500.}; 
    // std::vector<double> ptValues { 20., 30., 40., 50., 60., 70., 80., 90., 100., 110. };
    // std::vector<double> ptValues { 120., 130., 140., 150., 160., 180., 200., 250., 300., 350. };
    // std::vector<double> ptValues { 500., 700., 1000., 1500., 6500. }; 
    TH1D *hJECvsEta[ptValues.size()];
    TH1D *hJECvsPt[104];

    double ptVals[ptValues.size()];
    for (int i=0; i<ptValues.size(); i++) {
        ptVals[i] = ptValues[i];
    }
    int ptBins = ptValues.size() - 1;

    //
    // eta dependence
    //

    TLegend *leg;
    TCanvas *cEta = new TCanvas("cEta", "cEta", 800, 800);
    cEta->cd();
    setPadStyle();

    leg = new TLegend(0.45, 0.45, 0.75, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
    
    // pT bins
    for (int i=0; i<(ptValues.size()-1); i++) {

        double pt = (ptValues[i] + ptValues[i+1]) / 2.;
        hJECvsEta[i] = new TH1D( Form("hJEC_eta_%d", i), Form("JEC vs #eta for pt=%.0f;#eta;JEC", pt), 104, -5.2, 5.2);
        set1DStyle(hJECvsEta[i], 2);
        hJECvsEta[i]->Sumw2();
        hJECvsEta[i]->SetLineColor(i+1);
        hJECvsEta[i]->SetMarkerColor(i+1);
        hJECvsEta[i]->SetLineWidth(2);
        
        // Eta bins
        for (int j=1; j<104; j++) {

            double eta = -5.2 + j * 0.1;
            double phi = 0.;
            jec->SetJetPT(pt);
            jec->SetJetEta(eta);
            jec->SetJetPhi(phi);
            double ptCorr = jec->GetCorrectedPT();
            double jecFactor = ptCorr / pt;
            hJECvsEta[i]->SetBinContent(j, jecFactor);
            hJECvsEta[i]->SetBinError(j, 0.001);
            // std::cout << Form("pT: %.0f eta: %2.1f jec: %3.2f", pt, eta, jecFactor) << std::endl;
        } // for (int j=1; j<=104; j++)

        cEta->cd();
        // hJECvsEta[i]->SetLineWidth(3);
        if ( i == 0 ) {
            hJECvsEta[i]->Draw();
            hJECvsEta[i]->GetYaxis()->SetRangeUser(0.8, 2.2);
            hJECvsEta[i]->GetXaxis()->SetRangeUser(-3.1, 3.1);
        }
        else {
            hJECvsEta[i]->Draw("same");
        }
        leg->AddEntry(hJECvsEta[i], Form("p_{T}=%.0f GeV", pt), "lp");
    } // for (int i=0; i<(ptValues.size()-1); i++)

    leg->Draw();

    //
    // pT dependence
    //

    TCanvas *cPt = new TCanvas("cPt", "cPt", 800, 800);
    cPt->cd();
    setPadStyle();
    TLegend *leg2;
    leg2 = new TLegend(0.45, 0.45, 0.75, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->SetTextFont(42);

    // Loop over eta bins
    for (int i=1; i<104; i++) {

        double eta = -5.2 + i * 0.1;
        hJECvsPt[i] = new TH1D( Form("hJEC_pt_%d", i), Form("JEC vs p_{T} for #eta=%3.2f;p_{T} (GeV);JEC", eta), ptBins, ptVals);
        set1DStyle(hJECvsPt[i], 2);
        hJECvsPt[i]->Sumw2();
        hJECvsPt[i]->SetLineColor(i+1);
        hJECvsPt[i]->SetMarkerColor(i+1);
        hJECvsPt[i]->SetLineWidth(2);
        
        // Loop over pT bins
        for (int j=0; j<(ptValues.size()-1); j++) {

            double pt = (ptValues[j] + ptValues[j+1]) / 2.;
            double phi = 0.;
            jec->SetJetPT(pt);
            jec->SetJetEta(eta);
            jec->SetJetPhi(phi);
            double ptCorr = jec->GetCorrectedPT();
            double jecFactor = ptCorr / pt;
            hJECvsPt[i]->SetBinContent(j+1, jecFactor);
            hJECvsPt[i]->SetBinError(j+1, 0.001);
            // std::cout << Form("pT: %.0f eta: %2.1f jec: %3.2f", pt, eta, jecFactor) << std::endl;
        } // for (int j=0; j<(ptValues.size()-1); j++)

        cPt->cd();
        if ( i >= 53 && i < 71) {
            hJECvsPt[i]->Draw("same");
            hJECvsPt[i]->SetLineColor(i - 53 + 1);
            hJECvsPt[i]->SetMarkerColor(i - 53 + 1);
            hJECvsPt[i]->GetYaxis()->SetRangeUser(0.95, 1.35);
            hJECvsPt[i]->GetXaxis()->SetRangeUser(20., 300.);
            // hJECvsPt[i]->Draw("same");
            leg2->AddEntry(hJECvsPt[i], Form("#eta=%3.2f", eta), "lp");
        }
    } // for (int i=1; i<104; i++)

    leg2->Draw();

    // f->Write();
    // f->Close();
}