// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TRatioPlot.h"
#include "TPad.h"
#include "TF1.h"
#include "TGraphErrors.h"

// C++ headers
#include <iostream>
#include <vector>
#include <sys/stat.h>

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
void plotCMSHeader(int collSystem = 0, double energy = 5.02) {
    // collSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV
    TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize(0.05);
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}");
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.6, 0.93, Form("%s #sqrt{s_{NN}} = %3.2f TeV", collSystemStr.Data(), energy) );
    t.SetTextSize(0.05);
}

//________________
void set1DStyle(TH1 *h, Int_t type = 0, Bool_t doRenorm = kFALSE) {
    Int_t markerStyle = 20; // Full circle
    Double_t markerSize = 1.3;
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
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetNdivisions(208);
    h->GetYaxis()->SetNdivisions(208);    
    h->GetYaxis()->SetTitleOffset(1.5);

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
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.15);
}

//________________
void recalculateStatisticalUncertainty(TH1D *h) {
    for (int iBin{1}; iBin<=h->GetNbinsX(); iBin++) {
        double content = h->GetBinContent(iBin);
        double error = h->GetBinError(iBin);
        if (content > 0) {
            h->SetBinError(iBin, sqrt(content));
        } else {
            h->SetBinError(iBin, 0.);
        }
    }
}

//________________
void makeFullEtaFromForwardAndBackward(TH1D *hFull, TH1D *hForward, TH1D *hBackward) {

    int nBinsForward = hForward->GetNbinsX();
    int nBinsBackward = hBackward->GetNbinsX();

    for (int iBin{1}; iBin<=nBinsBackward; iBin++) {
        double binCenter = hBackward->GetBinCenter(iBin);
        // Flip the sign of the bin center for backward distribution since it's in negative eta
        binCenter *= -1;
        double binContent = hBackward->GetBinContent(iBin);
        double binError = hBackward->GetBinError(iBin);
        
        int fullBin = hFull->GetXaxis()->FindBin(binCenter);
        hFull->SetBinContent(fullBin, binContent);
        hFull->SetBinError(fullBin, binError);
    }

    for (int iBin{1}; iBin<=nBinsForward; iBin++) {
        double binCenter = hForward->GetBinCenter(iBin);
        double binContent = hForward->GetBinContent(iBin);
        double binError = hForward->GetBinError(iBin);
        
        int fullBin = hFull->GetXaxis()->FindBin(binCenter);
        hFull->SetBinContent(fullBin, binContent);
        hFull->SetBinError(fullBin, binError);
    }
}

//________________
void makeForwardBackwardFromFullEta(TH1D *hRatio, TH1D *hForward, TH1D *hBackward) {
    hRatio->Divide(hForward, hBackward);
}

//________________
void makeForwardBackwardDoubleRatio(TH1D *hDR, TH1D *hNpdfForward, TH1D *hNpdfBackward, TH1D *hPdfForward, TH1D *hPdfBackward) {
    TH1D *hRatioForward = (TH1D*)hNpdfForward->Clone("hRatioForward");
    hRatioForward->Divide(hPdfForward);
    TH1D *hRatioBackward = (TH1D*)hNpdfBackward->Clone("hRatioBackward");
    hRatioBackward->Divide(hPdfBackward);
    hDR->Divide(hRatioForward, hRatioBackward);
}

//________________
// Rescale a 1D histogram by its area-weighted integral
// This function rescales a 1D histogram by its area-weighted integral.
// It calculates the total area-weighted sum of the histogram and then normalizes each bin content and error accordingly.
// The function assumes that the histogram has non-equidistant axes.
void rescaleHisto1D(TH1* h1) {
    if (!h1) {
        std::cerr << "Null histogram pointer!" << std::endl;
        return;
    }

    const int nBins = h1->GetNbinsX();
    TAxis* xAxis = h1->GetXaxis();

    double total = 0.0;

    // First pass: calculate total sum (content × bin width)
    for (int i = 1; i <= nBins; ++i) {
        double width = xAxis->GetBinWidth(i);
        double content = h1->GetBinContent(i);
        total += content * width;
    }

    // Second pass: normalize content and error
    if (total > 0) {
        for (int i = 1; i <= nBins; ++i) {
            double content = h1->GetBinContent(i);
            double error   = h1->GetBinError(i);

            h1->SetBinContent(i, content / total);
            h1->SetBinError(i, error / total);
        }
    } 
    else {
        std::cerr << "Warning: total area-normalized sum is zero!" << std::endl;
    }
    // h1->Scale(1.0 / h1->Integral()); // Normalize to unity
}

//________________
void rescaleForwardBackward(TH1D *hForward, TH1D *hBackward) {
    // Rescale forward and backward histograms
    // to have the same integral
    double intForward = hForward->Integral();
    double intBackward = hBackward->Integral();
    double scaleFactor = intForward + intBackward;
    hForward->Scale( intForward / scaleFactor );
    hBackward->Scale( intBackward / scaleFactor );
}

//________________
void recalculateUncertaintyOfNpdf(TH1D *h) {

    TString fitName = Form("%s_fit", h->GetName());
    int lastNonZeroBin = h->FindLastBinAbove(0);
    TF1 *f = new TF1(fitName.Data(), "pol3", 
                     h->GetXaxis()->GetXmin(), 
                     h->GetXaxis()->GetBinUpEdge(lastNonZeroBin) );
    f->SetLineColor( h->GetLineColor() );
    f->SetLineWidth(3);
    h->Fit(fitName.Data(), "MRE0");

    for (int iBin{1}; iBin<=lastNonZeroBin; iBin++) {
        double x = h->GetXaxis()->GetBinCenter(iBin);
        double y = h->GetBinContent(iBin);
        if ( abs(y) < 1e-5 ) continue; // Skip bins with zero content to avoid issues with fitting
        double yFit = f->Eval(x);
        double newUnc = std::abs(y - yFit);
        h->SetBinError(iBin, newUnc);
    }
}


//________________
// Function calculates forward, backward and ratio of forward/backward distribution from
// the full 1D histograms in the CM frame
void recalculateFBRatioFromFullDistribution(TH1* hRatio, TH1* hForward, TH1* hBackward, TH1* hFull) {
    if (!hRatio || !hForward || !hBackward || !hFull) {
        std::cerr << "Error: Null histogram pointer passed." << std::endl;
        return;
    }

    if (hRatio->GetNbinsX() != hForward->GetNbinsX() || 
        hForward->GetNbinsX() != hBackward->GetNbinsX()) {
        std::cerr << "Error: Histograms have different number of bins." << std::endl;
        return;
    }

    // std::cout << Form("Number of bins Full: %d Forward: %d Backward: %d Ratio: %d\n",
    //                   hFull->GetNbinsX(), hForward->GetNbinsX(), hBackward->GetNbinsX(), hRatio->GetNbinsX());

    int middleBin = hFull->FindBin(0.0001);

    // Make forward
    for (int i=middleBin; i<=hFull->GetNbinsX(); i++) {
        double binVal = hFull->GetBinContent(i);
        double binError = hFull->GetBinError(i);
        hForward->SetBinContent(i-middleBin+1, binVal);
        hForward->SetBinError(i-middleBin+1, binError);
    } // for (int i=middleBin; i<=hFull->GetNbinsX(); i++)

    // Make backward
    for (int i=middleBin-1; i>0; i--) {
        double binVal = hFull->GetBinContent(i);
        double binError = hFull->GetBinError(i);
        hBackward->SetBinContent(middleBin-i, binVal);
        hBackward->SetBinError(middleBin-i, binError);
    } // for (int i=1; i<middleBin; i++)

    hRatio->Divide(hForward, hBackward, 1., 1.);
    for (int i=1; i<=hRatio->GetNbinsX(); i++) {
        double binVal = hRatio->GetBinContent(i);
        double binError = hRatio->GetBinError(i);
    }
}

//________________
// Plot dijet eta comparison of reco, gen, ref and refSel
void dijetClosures(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    
    // collisionSystem: 0 = pp, 1 = pPb, 2 = PbPb
    // energy in TeV

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // Dijet ptAve binning
    // double dijetPtVals[] {  50.,  60., 70.,  80.,  90.,
    //                         100., 110.,  120., 130., 140.,
    //                         150., 160.,  180., 200., 250., 
    //                         300., 500.};
    double dijetPtVals[] {  50.,  60., 80., 100., 120., 140.,
                            160.,  180., 200., 250., 300., 500.};
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    // Dijet selection ranges |eta_CM|<1.4, 1.5, 1.6, 1.7, 1.8, 1.9
    int dijetEtaCMIntVals[] { 14, 15, 16, 17, 18, 19 };
    float dijetEtaCMVals[] { 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 };
    int sizeOfDijetEtaCMVals = sizeof(dijetEtaCMVals)/sizeof(dijetEtaCMVals[0]);
    int rebinPseudoRapidity = 2; // Rebin factor for the pseudorapidity axis of the 2D histograms

    // Colors for plotting
    static constexpr std::array<const char*, 7> p8Colors {"kBlue", "kRed", "kMagenta", "kP8Orange", "kP8Green", "kP8Azure", "kBlack"};

    // Histograms for dijet pt vs eta CM frame
    TH2D *hRecoDijetPtEtaCMForward[6]{nullptr};
    TH2D *hRecoDijetPtEtaCMBackward[6]{nullptr};
    TH2D *hGenDijetPtEtaCMForward[6]{nullptr};
    TH2D *hGenDijetPtEtaCMBackward[6]{nullptr};
    TH2D *hRefDijetPtEtaCMForward[6]{nullptr};
    TH2D *hRefDijetPtEtaCMBackward[6]{nullptr};
    TH2D *hRefSelDijetPtEtaCMForward[6]{nullptr};
    TH2D *hRefSelDijetPtEtaCMBackward[6]{nullptr};

    // Loop over dijet etaCM selections and retrieve 2D distributions
    for (int iEtaCM = 0; iEtaCM < sizeOfDijetEtaCMVals; iEtaCM++) {

        //
        // Retrieve reco histograms
        //
        hRecoDijetPtEtaCMForward[iEtaCM] = dynamic_cast<TH2D *>( f->Get( Form("hRecoDijetPtEtaForwardArr_%d", iEtaCM) ) );
        if ( !hRecoDijetPtEtaCMForward[iEtaCM] ) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << iEtaCM << " not found in data file." << std::endl;
            return;
        }
        hRecoDijetPtEtaCMForward[iEtaCM]->SetDirectory(0);
        hRecoDijetPtEtaCMForward[iEtaCM]->RebinY(rebinPseudoRapidity);
        hRecoDijetPtEtaCMForward[iEtaCM]->SetName( Form("hRecoDijetPtEtaCMForward_%d", iEtaCM) );

        hRecoDijetPtEtaCMBackward[iEtaCM] = dynamic_cast<TH2D *>( f->Get( Form("hRecoDijetPtEtaBackwardArr_%d", iEtaCM) ) );
        if ( !hRecoDijetPtEtaCMBackward[iEtaCM] ) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << iEtaCM << " not found in data file." << std::endl;
            return;
        }
        hRecoDijetPtEtaCMBackward[iEtaCM]->SetDirectory(0);
        hRecoDijetPtEtaCMBackward[iEtaCM]->RebinY(rebinPseudoRapidity);
        hRecoDijetPtEtaCMBackward[iEtaCM]->SetName( Form("hRecoDijetPtEtaCMBackward_%d", iEtaCM) );

        //
        // Retrieve gen histograms
        //
        hGenDijetPtEtaCMForward[iEtaCM] = dynamic_cast<TH2D *>( f->Get( Form("hGenDijetPtEtaForwardArr_%d", iEtaCM) ) );
        if ( !hGenDijetPtEtaCMForward[iEtaCM] ) {
            std::cerr << "Error: hGenDijetPtEtaForwardArr_" << iEtaCM << " not found in data file." << std::endl;
            return;
        }
        hGenDijetPtEtaCMForward[iEtaCM]->SetDirectory(0);
        hGenDijetPtEtaCMForward[iEtaCM]->RebinY(rebinPseudoRapidity);
        hGenDijetPtEtaCMForward[iEtaCM]->SetName( Form("hGenDijetPtEtaCMForward_%d", iEtaCM) );

        hGenDijetPtEtaCMBackward[iEtaCM] = dynamic_cast<TH2D *>( f->Get( Form("hGenDijetPtEtaBackwardArr_%d", iEtaCM) ) );
        if ( !hGenDijetPtEtaCMBackward[iEtaCM] ) {
            std::cerr << "Error: hGenDijetPtEtaBackwardArr_" << iEtaCM << " not found in data file." << std::endl;
            return;
        }
        hGenDijetPtEtaCMBackward[iEtaCM]->SetDirectory(0);
        hGenDijetPtEtaCMBackward[iEtaCM]->RebinY(rebinPseudoRapidity);
        hGenDijetPtEtaCMBackward[iEtaCM]->SetName( Form("hGenDijetPtEtaCMBackward_%d", iEtaCM) );

    } // for (int iEtaCM = 0; iEtaCM < sizeOfDijetEtaCMVals; iEtaCM++)

    // Close the input file as we have retrieved all necessary histograms
    f->Close();

    //
    // 1D histograms for forward, backward and ratio of forward/backward distributions
    //

    TH1D *hRecoDijetEtaCMForward[6]{nullptr};
    TH1D *hRecoDijetEtaCMBackward[6]{nullptr};
    TH1D *hRecoDijetEtaCMForward2Backward[6]{nullptr};
    TH1D *hRecoDijetEtaCMFull[6]{nullptr};

    TH1D *hGenDijetEtaCMForward[6]{nullptr};
    TH1D *hGenDijetEtaCMBackward[6]{nullptr};
    TH1D *hGenDijetEtaCMForward2Backward[6]{nullptr};
    TH1D *hGenDijetEtaCMFull[6]{nullptr};

    TLatex t;
    t.SetTextFont( 42 );
    t.SetTextSize( 0.04 );

    TLegend *leg;
    TCanvas *c = new TCanvas( "c", "c", 1000, 1000 );

    // Loop over dijet ptAve bins
    for (int iPt = 0; iPt < sizeOfPtVals-1; iPt++) {

        double ptLow = dijetPtVals[iPt];
        double ptHi = dijetPtVals[iPt+1];

        // Select the appropriate bin range for the ptAve selection in the 2D histograms
        int ptLowBin = hRecoDijetPtEtaCMForward[0]->GetXaxis()->FindBin(ptLow);
        int ptHiBin = hRecoDijetPtEtaCMForward[0]->GetXaxis()->FindBin(ptHi) - 1; // Subtract 1 to include the upper edge of the last bin
        if (ptHiBin < ptLowBin) {
            std::cerr << "Error: Invalid pt bin range for ptLow = " << ptLow << " and ptHi = " << ptHi << std::endl;
            continue;
        }

        // Loop over dijet etaCM selections and create 1D projections for forward, backward and make their ratios
        for (int iEtaCM = 0; iEtaCM < sizeOfDijetEtaCMVals; iEtaCM++) {

            // Create 1D histograms for forward, backward and ratio of forward/backward distributions
            hRecoDijetEtaCMForward[iEtaCM] = hRecoDijetPtEtaCMForward[iEtaCM]->ProjectionY( Form("hRecoDijetEtaCMForward_pt_%d_%d_etaCM_%d", (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), ptLowBin, ptHiBin );
            set1DStyle(hRecoDijetEtaCMForward[iEtaCM], 0);

            hRecoDijetEtaCMBackward[iEtaCM] = hRecoDijetPtEtaCMBackward[iEtaCM]->ProjectionY( Form("hRecoDijetEtaCMBackward_pt_%d_%d_etaCM_%d", (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), ptLowBin, ptHiBin );
            set1DStyle(hRecoDijetEtaCMBackward[iEtaCM], 0);

            hRecoDijetEtaCMForward2Backward[iEtaCM] = (TH1D*)hRecoDijetEtaCMForward[iEtaCM]->Clone( Form("hRecoDijetEtaCMForward2Backward_pt_%d_%d_etaCM_%d", (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]) );
            hRecoDijetEtaCMForward2Backward[iEtaCM]->Divide(hRecoDijetEtaCMBackward[iEtaCM]);
            hRecoDijetEtaCMForward2Backward[iEtaCM]->SetLineColor( p8Colors[iEtaCM] );
            hRecoDijetEtaCMForward2Backward[iEtaCM]->SetMarkerColor( p8Colors[iEtaCM] );

            if (hRecoDijetEtaCMFull[iEtaCM]) { hRecoDijetEtaCMFull[iEtaCM]->Delete(); }
            hRecoDijetEtaCMFull[iEtaCM] = new TH1D( Form("hRecoDijetEtaCMFull_pt_%d_%d_etaCM_%d", 
                                                        (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), 
                                                    Form("hRecoDijetEtaCMFull_pt_%d_%d_etaCM_%d;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", 
                                                        (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), 
                                                    20, -2., 2. );
            hRecoDijetEtaCMFull[iEtaCM]->Sumw2();
            set1DStyle(hRecoDijetEtaCMFull[iEtaCM], 0);
            hRecoDijetEtaCMFull[iEtaCM]->SetLineColor( p8Colors[iEtaCM] );
            hRecoDijetEtaCMFull[iEtaCM]->SetMarkerColor( p8Colors[iEtaCM] );
            makeFullEtaFromForwardAndBackward(hRecoDijetEtaCMFull[iEtaCM], hRecoDijetEtaCMForward[iEtaCM], hRecoDijetEtaCMBackward[iEtaCM]);
            hRecoDijetEtaCMFull[iEtaCM]->Scale( 1./hRecoDijetEtaCMFull[iEtaCM]->Integral() );
            
            // Create gen projections for forward, backward and ratio of forward/backward distributions
            hGenDijetEtaCMForward[iEtaCM] = hGenDijetPtEtaCMForward[iEtaCM]->ProjectionY( Form("hGenDijetEtaCMForward_pt_%d_%d_etaCM_%d", (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), ptLowBin, ptHiBin );
            set1DStyle(hGenDijetEtaCMForward[iEtaCM], 1);

            hGenDijetEtaCMBackward[iEtaCM] = hGenDijetPtEtaCMBackward[iEtaCM]->ProjectionY( Form("hGenDijetEtaCMBackward_pt_%d_%d_etaCM_%d", (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), ptLowBin, ptHiBin );
            set1DStyle(hGenDijetEtaCMBackward[iEtaCM], 1);

            hGenDijetEtaCMForward2Backward[iEtaCM] = (TH1D*)hGenDijetEtaCMForward[iEtaCM]->Clone( Form("hGenDijetEtaCMForward2Backward_pt_%d_%d_etaCM_%d", (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]) );
            hGenDijetEtaCMForward2Backward[iEtaCM]->Divide(hGenDijetEtaCMBackward[iEtaCM]);
            hGenDijetEtaCMForward2Backward[iEtaCM]->SetLineColor( p8Colors[iEtaCM] );
            hGenDijetEtaCMForward2Backward[iEtaCM]->SetMarkerColor( p8Colors[iEtaCM] );

            if (hGenDijetEtaCMFull[iEtaCM]) { hGenDijetEtaCMFull[iEtaCM]->Delete(); }
            hGenDijetEtaCMFull[iEtaCM] = new TH1D( Form("hGenDijetEtaCMFull_pt_%d_%d_etaCM_%d", 
                                                        (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), 
                                                    Form("hGenDijetEtaCMFull_pt_%d_%d_etaCM_%d;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", 
                                                        (int)ptLow, (int)ptHi, dijetEtaCMIntVals[iEtaCM]), 
                                                    20, -2., 2. );
            hGenDijetEtaCMFull[iEtaCM]->Sumw2();
            set1DStyle(hGenDijetEtaCMFull[iEtaCM], 1);
            hGenDijetEtaCMFull[iEtaCM]->SetLineColor( p8Colors[iEtaCM] );
            hGenDijetEtaCMFull[iEtaCM]->SetMarkerColor( p8Colors[iEtaCM] );
            makeFullEtaFromForwardAndBackward(hGenDijetEtaCMFull[iEtaCM], hGenDijetEtaCMForward[iEtaCM], hGenDijetEtaCMBackward[iEtaCM]);
            hGenDijetEtaCMFull[iEtaCM]->Scale( 1./hGenDijetEtaCMFull[iEtaCM]->Integral() );

        } // for (int iEtaCM = 0; iEtaCM < sizeOfDijetEtaCMVals; iEtaCM++)

    } // for (int iPt = 0; iPt < sizeOfPtVals-1; iPt++)
}


//________________
//
// Plot comparison of inclusive jet distributions for different pT and eta bins
// Compare both reco data to reco MC and reco data to gen MC
//
// void data2mcInclusiveJetComparison(TFile *fData, TFile *fMc, int collSystem = 1, double energy = 8.16, TString date = "20250224") {
    
//     // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
//     // energy in TeV

//     TString prefix = "mb_";
//     TString dataFName = fData->GetName();
//     if ( dataFName.Contains("jet60", TString::kIgnoreCase) ) {
//         prefix = "jet60";
//     }
//     else if ( dataFName.Contains("jet80", TString::kIgnoreCase) ) {
//         prefix = "jet80";
//     }
//     else if ( dataFName.Contains("jet100", TString::kIgnoreCase) ) {
//         prefix = "jet100";
//     }

//     // Direction 
//     TString directionStr;
//     TString filename = fMc->GetName();
//     if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
//         directionStr = "Pbgoing";
//     } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
//         directionStr = "pgoing";
//     } else {
//         directionStr = "combined";
//     }

//     TString collSystemStr = (collSystem == 0) ? "pp" : (collSystem == 1) ? "pPb" : "PbPb";
//     collSystemStr += Form("%d", (int)(1000 * energy));
//     collSystemStr += Form("_%s", prefix.Data());

//     TH2D *hRecoDataPtEta = dynamic_cast<TH2D *>(fData->Get("hRecoInclusiveAllJetPtEtaCM"));
//     if (!hRecoDataPtEta) {
//         std::cerr << "Error: hRecoDataPtEta not found in data file!" << std::endl;
//         return;
//     }
//     hRecoDataPtEta->SetName("hRecoDataPtEta");
//     TH2D *hRecoMcPtEta = dynamic_cast<TH2D *>(fMc->Get("hRecoInclusiveAllJetPtEtaCM"));
//     if (!hRecoMcPtEta) {
//         std::cerr << "Error: hRecoMcPtEta not found in MC file!" << std::endl;
//         return;
//     }
//     hRecoMcPtEta->SetName("hRecoMcPtEta");
//     // TH3D *hRecoMcPtEtaPtHat = dynamic_cast<TH3D *>(fMc->Get("hRecoMatchedJetPtEtaPtHat"));
//     // hRecoMcPtEtaPtHat->SetName("hRecoMcPtEtaPtHat");
//     TH2D *hGenMcPtEta = dynamic_cast<TH2D *>(fMc->Get("hGenInclusiveJetPtEtaCM"));
//     if (!hGenMcPtEta) {
//         std::cerr << "Error: hGenMcPtEta not found in MC file!" << std::endl;
//         return;
//     }
//     hGenMcPtEta->SetName("hGenMcPtEta");

//     // jet pT and eta binning
//     double jetPtVals[] = {40., 80., 120., 180., 250., 500.};
//     int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

//     double jetEtaVals[] = {-3.6, -2.4, -1.6, 1.6, 2.4, 3.6};
//     int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);


//     TH1D *hRecoDataEta{nullptr};
//     TH1D *hRecoMcEta{nullptr};
//     TH1D *hGenMcEta{nullptr};

//     TCanvas *c = new TCanvas("c", "c", 1000, 1000);

//     // Jet pT binning
//     for ( int i = 0; i < sizeOfJetPtVals-1; i++) {

//         double jetPtLow = jetPtVals[i];
//         double jetPtHigh = jetPtVals[i+1];
//         int jetPtBinLow = hRecoDataPtEta->GetYaxis()->FindBin(jetPtLow);
//         int jetPtBinHigh = hRecoDataPtEta->GetYaxis()->FindBin(jetPtHigh)-1;

//         //
//         // Reco data 2 Reco MC comparison
//         //
//         // Make projection of the 3D histograms
//         hRecoDataEta = dynamic_cast<TH1D *>(hRecoDataPtEta->ProjectionX(Form("hRecoDataEta_%d", i),
//                                             jetPtBinLow, jetPtBinHigh));
//         if (!hRecoDataEta) {
//             std::cerr << Form("Error: hRecoDataEta_%d not found!", i) << std::endl;
//             return;
//         }
//         rescaleHisto1D(hRecoDataEta);
//         set1DStyle(hRecoDataEta, 0);

//         hRecoMcEta = dynamic_cast<TH1D *>(hRecoMcPtEta->ProjectionX(Form("hRecoMcEta_%d", i),
//                                           jetPtBinLow, jetPtBinHigh));
//         if (!hRecoMcEta) {
//             std::cerr << Form("Error: hRecoMcEta_%d not found!", i) << std::endl;
//             return;
//         }
//         rescaleHisto1D(hRecoMcEta);
//         set1DStyle(hRecoMcEta, 1);

//         // Plot distributions
//         drawJecComparison(c, hRecoDataEta, hRecoMcEta,
//             jetPtLow, jetPtHigh, 15,
//             "Data", "Reco MC", collSystem, energy, false, true,
//             hRecoDataPtEta->GetXaxis()->GetBinLowEdge(1),
//             hRecoDataPtEta->GetXaxis()->GetBinUpEdge(hRecoDataPtEta->GetXaxis()->GetNbins()));

//         c->SaveAs( Form("%s/%s_closureEta_Reco2Reco_pt_%d_%d.pdf", 
//             date.Data(), collSystemStr.Data(), (int)jetPtLow, (int)jetPtHigh ) );


//         //
//         // Reco data 2 Gen MC comparison
//         //

//         c->cd(); c->Clear(); setPadStyle();
//         // Make projection of the 3D histograms
//         hGenMcEta = dynamic_cast<TH1D *>(hGenMcPtEta->ProjectionX(Form("hGenMcEta_%d", i),
//             jetPtBinLow, jetPtBinHigh));
//         if (!hGenMcEta) {
//             std::cerr << Form("Error: hGenMcEta_%d not found!", i) << std::endl;
//             return;
//         }
//         rescaleHisto1D(hGenMcEta);
//         set1DStyle(hGenMcEta, 5);

//         // Plot distributions
//         drawJecComparison(c, hRecoDataEta, hGenMcEta,
//                           jetPtLow, jetPtHigh, 15,
//                           "Data", "Gen MC", collSystem, energy, false, true,
//                           hGenMcPtEta->GetXaxis()->GetBinLowEdge(1),
//                           hGenMcPtEta->GetXaxis()->GetBinUpEdge(hGenMcPtEta->GetXaxis()->GetNbins())
//                         );

//         // Save canvas
//         c->SaveAs( Form("%s/%s_closureEta_Reco2Gen_pt_%d_%d.pdf", 
//                     date.Data(), collSystemStr.Data(), (int)jetPtLow, (int)jetPtHigh) );
//     } // for (unsigned int i = 0; i < jetPtBinsLow.size(); i++)

//     if (c) { delete c; c = nullptr; }
// }

//________________
//
// Plot comparison of dijet eta distributions for different pTave bins
// Compare both reco data to reco MC and reco data to gen MC
//
// void data2mcDijetComparison(TFile *fData, TFile *fMc, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250224") {

//     // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
//     // energy in TeV
//     TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
//     collSystemStr += Form("%d", int(collisionEnergy * 1000) );

//     TString prefix = "mb";
//     TString dataFName = fData->GetName();
//     if ( dataFName.Contains("jet60", TString::kIgnoreCase) ) {
//         prefix = "jet60";
//     }
//     else if ( dataFName.Contains("jet80", TString::kIgnoreCase) ) {
//         prefix = "jet80";
//     }
//     else if ( dataFName.Contains("jet100", TString::kIgnoreCase) ) {
//         prefix = "jet100";
//     }

//     // Direction
//     TString directionStr;
//     TString filename = fMc->GetName();
//     if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
//         directionStr = "Pbgoing";
//     } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
//         directionStr = "pgoing";
//     } else {
//         directionStr = "combined";
//     }

//     // double dijetPtVals[] {  50.,  60.,   70.,  80.,  90.,
//     //                         100., 110.,  120., 130., 140.,
//     //                         150., 160.,  180., 200., 250., 
//     //                         300., 500.};
//     double dijetPtVals[] {  50., 80., 100., 120., 180., 200., 250., 300., 350., 500.};
//     int sizeOfDijetPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

//     //
//     // Declare histograms
//     //

//     // Data
//     TH3D *hRecoDataDijetPtEtaLab = dynamic_cast<TH3D *>(fData->Get("hRecoDijetPtEtaWeighted"));
//     if ( !hRecoDataDijetPtEtaLab ) {
//         std::cerr << "Data histogram hRecoDijetPtEtaWeighted not found in file." << std::endl; return;
//     }
//     hRecoDataDijetPtEtaLab->SetName("hRecoDataDijetPtEtaLab");
//     TH3D *hRecoDataDijetPtEtaCM = dynamic_cast<TH3D *>(fData->Get("hRecoDijetPtEtaCMWeighted"));
//     if ( !hRecoDataDijetPtEtaCM ) {
//         std::cerr << "Data histogram hRecoDijetPtEtaCMWeighted not found in file." << std::endl; return;
//     }
//     hRecoDataDijetPtEtaCM->SetName("hRecoDataDijetPtEtaCM");
//     TH2D *hRecoDataDijetPtEtaCMForward = dynamic_cast<TH2D *>(fData->Get("hRecoDijetPtEtaCMForwardWeighted"));
//     if ( !hRecoDataDijetPtEtaCMForward ) {
//         std::cerr << "Data histogram hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
//     }
//     hRecoDataDijetPtEtaCMForward->SetName("hRecoDataDijetPtEtaCMForward");
//     TH2D *hRecoDataDijetPtEtaCMBackward = dynamic_cast<TH2D *>(fData->Get("hRecoDijetPtEtaCMBackwardWeighted"));
//     if ( !hRecoDataDijetPtEtaCMBackward ) {
//         std::cerr << "Data histogram hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
//     }
//     hRecoDataDijetPtEtaCMBackward->SetName("hRecoDataDijetPtEtaCMBackward");

//     TH1D *hRecoDataDijetEtaLab{nullptr};
//     TH1D *hRecoDataDijetEtaCM{nullptr};
//     TH1D *hRecoDataDijetEtaForward{nullptr};
//     TH1D *hRecoDataDijetEtaBackward{nullptr};
//     TH1D *hRecoDataDijetFBRatio{nullptr};

//     // MC

//     // Reco
//     TH3D *hRecoMcDijetPtEtaLab = dynamic_cast<TH3D *>(fMc->Get("hRecoDijetPtEtaWeighted"));
//     if ( !hRecoMcDijetPtEtaLab ) {
//         std::cerr << "MC histogram hRecoDijetPtEtaWeighted not found in file." << std::endl; return;
//     }
//     hRecoMcDijetPtEtaLab->SetName("hRecoMcDijetPtEtaLab");
//     TH3D *hRecoMcDijetPtEtaCM = dynamic_cast<TH3D *>(fMc->Get("hRecoDijetPtEtaCMWeighted"));
//     if ( !hRecoMcDijetPtEtaCM ) {
//         std::cerr << "MC histogram hRecoDijetPtEtaCMWeighted not found in file." << std::endl; return;
//     }
//     hRecoMcDijetPtEtaCM->SetName("hRecoMcDijetPtEtaCM");
//     TH2D *hRecoMcDijetPtEtaCMForward = dynamic_cast<TH2D *>(fMc->Get("hRecoDijetPtEtaCMForwardWeighted"));
//     if ( !hRecoMcDijetPtEtaCMForward ) {
//         std::cerr << "MC histogram hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
//     }
//     hRecoMcDijetPtEtaCMForward->SetName("hRecoMcDijetPtEtaCMForward");
//     TH2D *hRecoMcDijetPtEtaCMBackward = dynamic_cast<TH2D *>(fMc->Get("hRecoDijetPtEtaCMBackwardWeighted"));
//     if ( !hRecoMcDijetPtEtaCMBackward ) {
//         std::cerr << "MC histogram hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
//     }
//     hRecoMcDijetPtEtaCMBackward->SetName("hRecoMcDijetPtEtaCMBackward");

//     TH1D *hRecoMcDijetEtaLab{nullptr};
//     TH1D *hRecoMcDijetEtaCM{nullptr};
//     TH1D *hRecoMcDijetEtaForward{nullptr};
//     TH1D *hRecoMcDijetEtaBackward{nullptr};
//     TH1D *hRecoMcDijetFBRatio{nullptr};

//     // Gen
//     TH3D *hGenMcDijetPtEtaLab = dynamic_cast<TH3D *>(fMc->Get("hGenDijetPtEtaWeighted"));
//     if ( !hGenMcDijetPtEtaLab ) {
//         std::cerr << "MC histogram hGenDijetPtEtaWeighted not found in file." << std::endl; return;
//     }
//     hGenMcDijetPtEtaLab->SetName("hGenMcDijetPtEtaLab");
//     TH3D *hGenMcDijetPtEtaCM = dynamic_cast<TH3D *>(fMc->Get("hGenDijetPtEtaCMWeighted"));
//     if ( !hGenMcDijetPtEtaCM ) {
//         std::cerr << "MC histogram hGenDijetPtEtaCMWeighted not found in file." << std::endl; return;
//     }
//     hGenMcDijetPtEtaCM->SetName("hGenMcDijetPtEtaCM");
//     TH2D *hGenMcDijetPtEtaCMForward = dynamic_cast<TH2D *>(fMc->Get("hGenDijetPtEtaCMForwardWeighted"));
//     if ( !hGenMcDijetPtEtaCMForward ) {
//         std::cerr << "MC histogram hGenDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
//     }
//     hGenMcDijetPtEtaCMForward->SetName("hGenMcDijetPtEtaCMForward");
//     TH2D *hGenMcDijetPtEtaCMBackward =  dynamic_cast<TH2D *>(fMc->Get("hGenDijetPtEtaCMBackwardWeighted"));
//     if ( !hGenMcDijetPtEtaCMBackward ) {
//         std::cerr << "MC histogram hGenDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
//     }
//     hGenMcDijetPtEtaCMBackward->SetName("hGenMcDijetPtEtaCMBackward");

//     TH1D *hGenMcDijetEtaLab{nullptr};
//     TH1D *hGenMcDijetEtaCM{nullptr};
//     TH1D *hGenMcDijetEtaForward{nullptr};
//     TH1D *hGenMcDijetEtaBackward{nullptr};
//     TH1D *hGenMcDijetFBRatio{nullptr};

//     // Ref
//     TH3D *hRefMcDijetPtEtaLab = dynamic_cast<TH3D *>(fMc->Get("hRefDijetPtEtaWeighted"));
//     if ( !hRefMcDijetPtEtaLab ) {
//         std::cerr << "MC histogram hRefDijetPtEtaWeighted not found in file." << std::endl; return;
//     }
//     hRefMcDijetPtEtaLab->SetName("hRefMcDijetPtEtaLab");
//     TH3D *hRefMcDijetPtEtaCM = dynamic_cast<TH3D *>(fMc->Get("hRefDijetPtEtaCMWeighted"));
//     if ( !hRefMcDijetPtEtaCM ) {
//         std::cerr << "MC histogram hRefDijetPtEtaCMWeighted not found in file." << std::endl; return;
//     }
//     hRefMcDijetPtEtaCM->SetName("hRefMcDijetPtEtaCM");
//     TH2D *hRefMcDijetPtEtaCMForward = dynamic_cast<TH2D *>(fMc->Get("hRefDijetPtEtaCMForwardWeighted"));
//     if ( !hRefMcDijetPtEtaCMForward ) {
//         std::cerr << "MC histogram hRefDijetPtEtaCMForwardWeighted not found in file." << std::endl; return;
//     }
//     hRefMcDijetPtEtaCMForward->SetName("hRefMcDijetPtEtaCMForward");
//     TH2D *hRefMcDijetPtEtaCMBackward = dynamic_cast<TH2D *>(fMc->Get("hRefDijetPtEtaCMBackwardWeighted"));
//     if ( !hRefMcDijetPtEtaCMBackward ) {
//         std::cerr << "MC histogram hRefDijetPtEtaCMBackwardWeighted not found in file." << std::endl; return;
//     }
//     hRefMcDijetPtEtaCMBackward->SetName("hRefMcDijetPtEtaCMBackward");

//     TH1D *hRefMcDijetEtaLab{nullptr};
//     TH1D *hRefMcDijetEtaCM{nullptr};
//     TH1D *hRefMcDijetEtaForward{nullptr};
//     TH1D *hRefMcDijetEtaBackward{nullptr};
//     TH1D *hRefMcDijetFBRatio{nullptr};

//     // Ratio of full distributions to gen
//     TH1D *hRecoData2GenDijetEtaLab{nullptr};
//     TH1D *hRecoMc2GenDijetEtaLab{nullptr};
//     TH1D *hRef2GenDijetEtaLab{nullptr};

//     TH1D *hRecoData2GenDijetEtaCM{nullptr};
//     TH1D *hRecoMc2GenDijetEtaCM{nullptr};
//     TH1D *hRef2GenDijetEtaCM{nullptr};

//     //
//     // Canvas 
//     //
//     TLatex t;
//     t.SetTextFont(42);
//     t.SetTextSize(0.05);
//     TCanvas *c = new TCanvas("c", "c", 1200, 1200);
//     setPadStyle();

//     //
//     // Loop over the dijet pTave bins, make projections and plot distributions
//     //
//     for (int i = 0; i < sizeOfDijetPtVals - 1; ++i) {

//         double ptLow = dijetPtVals[i];
//         double ptHi = dijetPtVals[i + 1];
//         int ptBinLow = hRecoDataDijetPtEtaLab->GetXaxis()->FindBin(ptLow);
//         int ptBinHigh = hRecoDataDijetPtEtaLab->GetXaxis()->FindBin(ptHi) - 1;

//         //
//         // Data
//         //

//         hRecoDataDijetEtaLab = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaLab->ProjectionY( Form("hRecoDataDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoDataDijetEtaLab ) {
//             std::cerr << Form("Data histogram not found: hRecoDataDijetEtaLab_%d", i) << std::endl; return;
//         }
//         hRecoDataDijetEtaLab->SetName( Form("hRecoDataDijetEtaLab_%d", i) );
//         set1DStyle( hRecoDataDijetEtaLab, 2 );
//         rescaleHisto1D( hRecoDataDijetEtaLab );

//         hRecoDataDijetEtaCM = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaCM->ProjectionY( Form("hRecoDataDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoDataDijetEtaCM ) {
//             std::cerr << Form("Data histogram not found: hRecoDataDijetEtaCM_%d", i) << std::endl; return;
//         }
//         hRecoDataDijetEtaCM->SetName( Form("hRecoDataDijetEtaCM_%d", i) );
//         set1DStyle( hRecoDataDijetEtaCM, 2 );
//         rescaleHisto1D( hRecoDataDijetEtaCM );

//         hRecoDataDijetEtaForward = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaCMForward->ProjectionY( Form("hRecoDataDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoDataDijetEtaForward ) {
//             std::cerr << Form("Data histogram not found: hRecoDataDijetEtaForward_%d", i) << std::endl; return;
//         }
//         hRecoDataDijetEtaForward->SetName( Form("hRecoDataDijetEtaForward_%d", i) );
//         set1DStyle( hRecoDataDijetEtaForward, 2 );
        
//         hRecoDataDijetEtaBackward = dynamic_cast<TH1D*>( hRecoDataDijetPtEtaCMBackward->ProjectionY( Form("hRecoDataDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoDataDijetEtaBackward ) {
//             std::cerr << Form("Data histogram not found: hRecoDataDijetEtaBackward_%d", i) << std::endl; return;
//         }
//         hRecoDataDijetEtaBackward->SetName( Form("hRecoDataDijetEtaBackward_%d", i) );
//         set1DStyle( hRecoDataDijetEtaBackward, 2 );

//         // Normalize forward and backward distributions
//         // double integralForward = hRecoDataDijetEtaForward->Integral();
//         // double integralBackward = hRecoDataDijetEtaBackward->Integral();
//         // hRecoDataDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
//         // hRecoDataDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

//         // Ratios
//         hRecoDataDijetFBRatio = dynamic_cast<TH1D*>( hRecoDataDijetEtaForward->Clone( Form("hRecoDataDijetFBRatio_%d", i) ) );
//         hRecoDataDijetFBRatio->SetName( Form("hRecoDataDijetFBRatio_%d", i) );
//         hRecoDataDijetFBRatio->Divide( hRecoDataDijetFBRatio, hRecoDataDijetEtaBackward, 1., 1. );

//         c->cd(); c->Clear(); setPadStyle();
//         gPad->SetGrid(1, 1);
//         hRecoDataDijetEtaLab->Draw();
//         hRecoDataDijetEtaLab->GetXaxis()->SetRangeUser(-3.0, 3.0);
//         hRecoDataDijetEtaLab->GetYaxis()->SetRangeUser(0.00001, 0.15);
//         t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
//         c->SaveAs( Form("%s/%s_%s_dijetEtaLab_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)dijetPtVals[i], (int)dijetPtVals[i+1]) );

//         c->cd(); c->Clear(); setPadStyle();
//         gPad->SetGrid(1, 1);
//         hRecoDataDijetEtaCM->Draw();
//         hRecoDataDijetEtaCM->GetXaxis()->SetRangeUser(-2.5, 2.5);
//         hRecoDataDijetEtaCM->GetYaxis()->SetRangeUser(0.00001, 0.15);
//         t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
//         c->SaveAs( Form("%s/%s_%s_dijetEtaCM_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

//         c->cd(); c->Clear(); setPadStyle();
//         gPad->SetGrid(1, 1);
//         hRecoDataDijetFBRatio->Draw();
//         hRecoDataDijetFBRatio->GetXaxis()->SetRangeUser(-2.5, 2.5);
//         hRecoDataDijetFBRatio->GetYaxis()->SetRangeUser(0.75, 1.25);
//         t.DrawLatexNDC( 0.35, 0.85, Form("%d < p_{T}^{ave} < %d GeV", (int)ptLow, (int)ptHi) );
//         c->SaveAs( Form("%s/%s_%s_dijetEtaFBRatio_data_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

//         //
//         // MC (reco)
//         //

//         // Eta full (lab frame)
//         hRecoMcDijetEtaLab = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaLab->ProjectionY( Form("hRecoMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoMcDijetEtaLab ) {
//             std::cerr << Form("MC histogram not found: hRecoDijetEta1DLab_%d", i) << std::endl; return;
//         }
//         hRecoMcDijetEtaLab->SetName( Form("hRecoMcDijetEtaLab_%d", i) );
//         set1DStyle( hRecoMcDijetEtaLab, 0 );
//         rescaleHisto1D( hRecoMcDijetEtaLab );

//         // Eta full (CM frame)
//         hRecoMcDijetEtaCM = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaCM->ProjectionY( Form("hRecoMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoMcDijetEtaCM ) {
//             std::cerr << Form("MC histogram not found: hRecoMcDijetEtaCM_%d", i) << std::endl; return;
//         }
//         hRecoMcDijetEtaCM->SetName( Form("hRecoMcDijetEtaCM_%d", i) );
//         set1DStyle( hRecoMcDijetEtaCM, 0 );
//         rescaleHisto1D( hRecoMcDijetEtaCM );

//         // Eta forward (CM frame)
//         hRecoMcDijetEtaForward = dynamic_cast<TH1D*>(hRecoMcDijetPtEtaCMForward->ProjectionY( Form("hRecoMcDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoMcDijetEtaForward ) {
//             std::cerr << Form("MC histogram not found: hRecoDijetEtaCMForward1D_%d", i) << std::endl; return;
//         }
//         hRecoMcDijetEtaForward->SetName( Form("hRecoMcDijetEtaForward_%d", i) );
//         set1DStyle( hRecoMcDijetEtaForward, 0 );

//         // Eta backward (CM frame)
//         hRecoMcDijetEtaBackward = dynamic_cast<TH1D*>( hRecoMcDijetPtEtaCMBackward->ProjectionY( Form("hRecoMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRecoMcDijetEtaBackward ) {
//             std::cerr << Form("MC histogram not found: hRecoMcDijetEtaBackward_%d", i) << std::endl; return;
//         }
//         hRecoMcDijetEtaBackward->SetName( Form("hRecoMcDijetEtaBackward_%d", i) );
//         set1DStyle( hRecoMcDijetEtaBackward, 0 );

//         // Normalize the forward and backward distributions
//         // integralForward = hRecoMcDijetEtaForward->Integral();
//         // integralBackward = hRecoMcDijetEtaBackward->Integral();
//         // hRecoMcDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
//         // hRecoMcDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

//         // Ratios
//         hRecoMcDijetFBRatio = dynamic_cast<TH1D*>( hRecoMcDijetEtaForward->Clone( Form("hRecoMcDijetFBRatio_%d", i) ) );
//         hRecoMcDijetFBRatio->SetName( Form("hRecoMcDijetFBRatio_%d", i) );
//         hRecoMcDijetFBRatio->Divide( hRecoMcDijetFBRatio, hRecoMcDijetEtaBackward, 1., 1. );

//         //
//         // MC (gen)
//         //

//         // Eta full (lab frame)
//         hGenMcDijetEtaLab = dynamic_cast<TH1D*>( hGenMcDijetPtEtaLab->ProjectionY( Form("hGenMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hGenMcDijetEtaLab ) {
//             std::cerr << Form("MC histogram not found: hGenMcDijetEtaLab_%d", i) << std::endl; return;
//         }
//         hGenMcDijetEtaLab->SetName( Form("hGenMcDijetEtaLab_%d", i) );
//         set1DStyle( hGenMcDijetEtaLab, 5 );
//         rescaleHisto1D( hGenMcDijetEtaLab );    

//         // Eta full (CM frame)
//         hGenMcDijetEtaCM = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCM->ProjectionY( Form("hGenMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hGenMcDijetEtaCM ) {
//             std::cerr << Form("MC histogram not found: hGenMcDijetEtaCM_%d", i) << std::endl; return;
//         }
//         hGenMcDijetEtaCM->SetName( Form("hGenMcDijetEtaCM_%d", i) );
//         set1DStyle( hGenMcDijetEtaCM, 5 );
//         rescaleHisto1D( hGenMcDijetEtaCM );
        
//         // Eta forward (CM frame)
//         hGenMcDijetEtaForward = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCMForward->ProjectionY( Form("hGenMcDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hGenMcDijetEtaForward ) {
//             std::cerr << Form("MC histogram not found: hGenMcDijetEtaForward_%d", i) << std::endl; return;
//         }
//         hGenMcDijetEtaForward->SetName( Form("hGenMcDijetEtaForward_%d", i) );
//         set1DStyle( hGenMcDijetEtaForward, 5 );

//         // Eta backward (CM frame)
//         hGenMcDijetEtaBackward = dynamic_cast<TH1D*>( hGenMcDijetPtEtaCMBackward->ProjectionY( Form("hGenMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hGenMcDijetEtaBackward ) {
//             std::cerr << Form("MC histogram not found: hGenMcDijetEtaBackward_%d", i) << std::endl; return;
//         }
//         hGenMcDijetEtaBackward->SetName( Form("hGenMcDijetEtaBackward_%d", i) );
//         set1DStyle( hGenMcDijetEtaBackward, 5 );
        
//         // Normalize the forward and backward distributions
//         // integralForward = hGenMcDijetEtaForward->Integral();
//         // integralBackward = hGenMcDijetEtaBackward->Integral();
//         // hGenMcDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
//         // hGenMcDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

//         // Ratios
//         hGenMcDijetFBRatio = dynamic_cast<TH1D*>( hGenMcDijetEtaForward->Clone( Form("hGenMcDijetFBRatio_%d", i) ) );
//         hGenMcDijetFBRatio->SetName( Form("hGenMcDijetFBRatio_%d", i) );
//         hGenMcDijetFBRatio->Divide( hGenMcDijetFBRatio, hGenMcDijetEtaBackward, 1., 1. );

//         //
//         // MC (ref)
//         //

//         // Eta full (lab frame)
//         hRefMcDijetEtaLab = dynamic_cast<TH1D*>( hRefMcDijetPtEtaLab->ProjectionY( Form("hRefMcDijetEtaLab_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRefMcDijetEtaLab ) {
//             std::cerr << Form("MC histogram not found: hRefMcDijetEtaLab_%d", i) << std::endl; return;
//         }
//         hRefMcDijetEtaLab->SetName( Form("hRefMcDijetEtaLab_%d", i) );
//         set1DStyle( hRefMcDijetEtaLab, 1 );
//         rescaleHisto1D( hRefMcDijetEtaLab );

//         // Eta full (CM frame)
//         hRefMcDijetEtaCM = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCM->ProjectionY( Form("hRefMcDijetEtaCM_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRefMcDijetEtaCM ) {
//             std::cerr << Form("MC histogram not found: hRefMcDijetEtaCM_%d", i) << std::endl; return;
//         }
//         hRefMcDijetEtaCM->SetName( Form("hRefMcDijetEtaCM_%d", i) );
//         set1DStyle( hRefMcDijetEtaCM, 1 );
//         rescaleHisto1D( hRefMcDijetEtaCM );

//         // Eta forward (CM frame)
//         hRefMcDijetEtaForward = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCMForward->ProjectionY( Form("hRefMcDijetEtaForward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRefMcDijetEtaForward ) {
//             std::cerr << Form("MC histogram not found: hRefMcDijetEtaForward_%d", i) << std::endl; return;
//         }
//         hRefMcDijetEtaForward->SetName( Form("hRefMcDijetEtaForward_%d", i) );
//         set1DStyle( hRefMcDijetEtaForward, 1 );
        

//         // Eta backward (CM frame)
//         hRefMcDijetEtaBackward = dynamic_cast<TH1D*>( hRefMcDijetPtEtaCMBackward->ProjectionY( Form("hRefMcDijetEtaBackward_%d", i), ptBinLow, ptBinHigh ) );
//         if ( !hRefMcDijetEtaBackward ) {
//             std::cerr << Form("MC histogram not found: hRefMcDijetEtaBackward_%d", i) << std::endl; return;
//         }
//         hRefMcDijetEtaBackward->SetName( Form("hRefMcDijetEtaBackward_%d", i) );
//         set1DStyle( hRefMcDijetEtaBackward, 1 );

//         // Normalize the forward and backward distributions
//         // integralForward = hRefMcDijetEtaForward->Integral();
//         // integralBackward = hRefMcDijetEtaBackward->Integral();
//         // hRefMcDijetEtaForward->Scale( integralForward / (integralForward + integralBackward) );
//         // hRefMcDijetEtaBackward->Scale( integralBackward / (integralForward + integralBackward) );

//         // Ratios
//         hRefMcDijetFBRatio = dynamic_cast<TH1D*>( hRefMcDijetEtaForward->Clone( Form("hRefMcDijetFBRatio_%d", i) ) );
//         hRefMcDijetFBRatio->SetName( Form("hRefMcDijetFBRatio_%d", i) );
//         hRefMcDijetFBRatio->Divide( hRefMcDijetFBRatio, hRefMcDijetEtaBackward, 1., 1. );

//         //
//         // Ratios of full distributions to gen in the laboratory 
//         //
//         hRecoData2GenDijetEtaLab = dynamic_cast<TH1D*>( hRecoDataDijetEtaLab->Clone( Form("hRecoData2GenDijetEtaLab_%d", i) ) );
//         hRecoData2GenDijetEtaLab->SetName( Form("hRecoData2GenDijetEtaLab_%d", i) );
//         hRecoData2GenDijetEtaLab->GetYaxis()->SetTitle( "Ratio to Gen" );
//         hRecoData2GenDijetEtaLab->Divide( hRecoData2GenDijetEtaLab, hGenMcDijetEtaLab, 1., 1. );

//         hRecoMc2GenDijetEtaLab = dynamic_cast<TH1D*>( hRecoMcDijetEtaLab->Clone( Form("hRecoMc2GenDijetEtaLab_%d", i) ) );
//         hRecoMc2GenDijetEtaLab->SetName( Form("hRecoMc2GenDijetEtaLab_%d", i) );
//         hRecoMc2GenDijetEtaLab->GetYaxis()->SetTitle( "Ratio to Gen" );
//         hRecoMc2GenDijetEtaLab->Divide( hRecoMc2GenDijetEtaLab, hGenMcDijetEtaLab, 1., 1., "b" );

//         hRef2GenDijetEtaLab = dynamic_cast<TH1D*>( hRefMcDijetEtaLab->Clone( Form("hRef2GenDijetEtaLab_%d", i) ) );
//         hRef2GenDijetEtaLab->SetName( Form("hRef2GenDijetEtaLab_%d", i) );
//         hRef2GenDijetEtaLab->GetYaxis()->SetTitle( "Ratio to Gen" );
//         hRef2GenDijetEtaLab->Divide( hRef2GenDijetEtaLab, hGenMcDijetEtaLab, 1., 1., "b" );

//         //
//         // Ratios of full distributions to gen in the center-of-mass frame
//         //
//         hRecoData2GenDijetEtaCM = dynamic_cast<TH1D*>( hRecoDataDijetEtaCM->Clone( Form("hRecoData2GenDijetEtaCM_%d", i) ) );
//         hRecoData2GenDijetEtaCM->SetName( Form("hRecoData2GenDijetEtaCM_%d", i) );
//         hRecoData2GenDijetEtaCM->GetYaxis()->SetTitle( "Ratio to Gen" );
//         hRecoData2GenDijetEtaCM->Divide( hRecoData2GenDijetEtaCM, hGenMcDijetEtaCM, 1., 1. );

//         hRecoMc2GenDijetEtaCM = dynamic_cast<TH1D*>( hRecoMcDijetEtaCM->Clone( Form("hRecoMc2GenDijetEtaCM_%d", i) ) );
//         hRecoMc2GenDijetEtaCM->SetName( Form("hRecoMc2GenDijetEtaCM_%d", i) );
//         hRecoMc2GenDijetEtaCM->GetYaxis()->SetTitle( "Ratio to Gen" );
//         hRecoMc2GenDijetEtaCM->Divide( hRecoMc2GenDijetEtaCM, hGenMcDijetEtaCM, 1., 1., "b" );

//         hRef2GenDijetEtaCM = dynamic_cast<TH1D*>( hRefMcDijetEtaCM->Clone( Form("hRef2GenDijetEtaCM_%d", i) ) );
//         hRef2GenDijetEtaCM->SetName( Form("hRef2GenDijetEtaCM_%d", i) );
//         hRef2GenDijetEtaCM->GetYaxis()->SetTitle( "Ratio to Gen" );
//         hRef2GenDijetEtaCM->Divide( hRef2GenDijetEtaCM, hGenMcDijetEtaCM, 1., 1., "b" );

//         //
//         // Plot comparisons
//         //


//         // Lab frame
//         c->cd();
//         drawDijetToGenComparison(c, hRecoMcDijetEtaLab, hGenMcDijetEtaLab, hRefMcDijetEtaLab, nullptr, hRecoDataDijetEtaLab,
//                             (int)ptLow, (int)ptHi, false, false, collisionSystem, collisionEnergy);
//         c->SaveAs( Form("%s/%s_%s_dijetEtaLab_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );


//         // CM frame
//         drawDijetToGenComparison(c, hRecoMcDijetEtaCM, hGenMcDijetEtaCM, hRefMcDijetEtaCM, nullptr, hRecoDataDijetEtaCM,
//                             (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
//         c->SaveAs( Form("%s/%s_%s_dijetEtaCM_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );



//         // Forward-backward ratio
//         drawDijetToGenComparison(c, hRecoMc2GenDijetEtaCM, nullptr, hRef2GenDijetEtaCM, nullptr, hRecoData2GenDijetEtaCM, 
//                             (int)ptLow, (int)ptHi, true, true, collisionSystem, collisionEnergy);
//         c->SaveAs( Form("%s/%s_%s_dijetEtaFBRatio_DataRecoRefGenComp_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

//         //
//         // Plot ratios
//         //

//         // Lab frame
//         drawDijetToGenRatio(c, hRecoMc2GenDijetEtaLab, hRef2GenDijetEtaLab, nullptr, hRecoData2GenDijetEtaLab,
//                        (int)ptLow, (int)ptHi, false, false, collisionSystem, collisionEnergy);
//         c->SaveAs( Form("%s/%s_%s_dijetEtaLab_DataRecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );

//         // CM frame
//         drawDijetToGenRatio(c, hRecoMc2GenDijetEtaCM, hRef2GenDijetEtaCM, nullptr, hRecoData2GenDijetEtaCM,
//                        (int)ptLow, (int)ptHi, true, false, collisionSystem, collisionEnergy);
//         c->SaveAs( Form("%s/%s_%s_dijetEtaCM_DataRecoRef2GenRatio_ptAve_%d_%d.pdf", date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi) );
      
//     } // for (int i = 0; i < sizeOfDijetPtVals - 1; ++i)

//     if (c) { delete c; c = nullptr; }
// }

//________________
// Plot eta and pT distributions of inclusive jets from Pythia and embedding.
// Compare dijet distributions.
// void pythia2embeddingSingleJet(TFile *fEmbedding, TFile *fPythia, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250305") {
//     // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
//     // energy in TeV

//     TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
//     collSystemStr += Form("%d", (int)(1000 * collisionEnergy));
//     collSystemStr += Form("_emb2pythia");

//     // Event ptHat binning
//     double ptHatVals[] { 15, 45, 90 };
//     int sizeOfPtHatVals = sizeof(ptHatVals)/sizeof(ptHatVals[0]);

//     // Jet pT binning
//     double jetPtVals[] = {30., 50., 90., 120., 200., 500., 1500.};
//     int sizeOfJetPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

//     // Jet eta binning
//     // double jetEtaVals[] = {-3.0, -2.4, -1.6, 1.6, 2.4, 3.0};
//     // int sizeOfJetEtaVals = sizeof(jetEtaVals)/sizeof(jetEtaVals[0]);

//     //
//     // Retrieve inclusive jet histograms
//     //

//     // Embedding
//     TH3D *hRecoPtEtaPtHatEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRecoInclusiveJetPtEtaPtHat"));
//     if (!hRecoPtEtaPtHatEmb) {
//         std::cerr << "Error: hRecoPtEtaPtHatEmb not found!" << std::endl;
//         return;
//     }
//     hRecoPtEtaPtHatEmb->SetName("hRecoPtEtaPtHatEmb");
//     TH3D *hGenPtEtaPtHatEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hGenInclusiveJetPtEtaPtHat"));
//     if (!hGenPtEtaPtHatEmb) {
//         std::cerr << "Error: hGenPtEtaPtHatEmb not found!" << std::endl;
//         return;
//     }
//     hGenPtEtaPtHatEmb->SetName("hGenPtEtaPtHatEmb");

//     TH1D *hRecoEtaEmb{nullptr};
//     TH1D *hGenEtaEmb{nullptr};
//     TH1D *hRecoPtEmb{nullptr};
//     TH1D *hGenPtEmb{nullptr};

//     // Pythia
//     TH3D *hRecoPtEtaPtHatPythia = dynamic_cast<TH3D *>(fPythia->Get("hRecoInclusiveJetPtEtaPtHat"));
//     if (!hRecoPtEtaPtHatPythia) {
//         std::cerr << "Error: hRecoPtEtaPtHatPythia not found!" << std::endl;
//         return;
//     }
//     hRecoPtEtaPtHatPythia->SetName("hRecoPtEtaPtHatPythia");
//     TH3D *hGenPtEtaPtHatPythia = dynamic_cast<TH3D *>(fPythia->Get("hGenInclusiveJetPtEtaPtHat"));
//     if (!hGenPtEtaPtHatPythia) {
//         std::cerr << "Error: hGenPtEtaPtHatPythia not found!" << std::endl;
//         return;
//     }
//     hGenPtEtaPtHatPythia->SetName("hGenPtEtaPtHatPythia");

//     TH1D *hRecoEtaPythia{nullptr};
//     TH1D *hGenEtaPythia{nullptr};
//     TH1D *hRecoPtPythia{nullptr};
//     TH1D *hGenPtPythia{nullptr};

//     TCanvas *c = new TCanvas("c", "c", 1000, 1000);
//     setPadStyle();
//     gPad->SetGrid();

//     // Loop over ptHat bins
//     for (size_t iPtHat = 0; iPtHat < sizeOfPtHatVals-1; ++iPtHat) {

//         double ptHatLow = ptHatVals[iPtHat];
//         double ptHatHigh = hRecoPtEtaPtHatEmb->GetZaxis()->GetBinUpEdge( hRecoPtEtaPtHatEmb->GetNbinsZ() );
//         int ptHatBinLow = hRecoPtEtaPtHatEmb->GetZaxis()->FindBin(ptHatLow);
//         int ptHatBinHigh = hRecoPtEtaPtHatEmb->GetZaxis()->FindBin(ptHatHigh)-1;

//         // Loop over ptHat bins and make projections of eta
//         for (size_t iJetPt = 0; iJetPt < sizeOfJetPtVals-1; ++iJetPt) {

//             double jetPtLow = jetPtVals[iJetPt];
//             double jetPtHigh = jetPtVals[iJetPt+1];
//             int jetPtBinLow = hRecoPtEtaPtHatEmb->GetYaxis()->FindBin(jetPtLow);
//             int jetPtBinHigh = hRecoPtEtaPtHatEmb->GetYaxis()->FindBin(jetPtHigh)-1;

//             // Embedding
//             hRecoEtaEmb = dynamic_cast<TH1D *>(hRecoPtEtaPtHatEmb->ProjectionX( Form("hRecoEtaEmb_%zu_%zu", iPtHat, iJetPt), 
//                                                                                 jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
//             rescaleHisto1D(hRecoEtaEmb);
//             set1DStyle(hRecoEtaEmb, 0);

//             hGenEtaEmb = dynamic_cast<TH1D *>(hGenPtEtaPtHatEmb->ProjectionX( Form("hGenEtaEmb_%zu_%zu", iPtHat, iJetPt), 
//                                                                                     jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
//             rescaleHisto1D(hGenEtaEmb);
//             set1DStyle(hGenEtaEmb, 0);

//             // Pythia
//             hRecoEtaPythia = dynamic_cast<TH1D *>(hRecoPtEtaPtHatPythia->ProjectionX(Form("hRecoEtaPythia_%zu_%zu", iPtHat, iJetPt), 
//                                                                                           jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
//             rescaleHisto1D(hRecoEtaPythia);
//             set1DStyle(hRecoEtaPythia, 1);

//             hGenEtaPythia = dynamic_cast<TH1D *>(hGenPtEtaPtHatPythia->ProjectionX(Form("hGenEtaPythia_%zu_%zu", iPtHat, iJetPt), 
//                                                                                                  jetPtBinLow, jetPtBinHigh, ptHatBinLow, ptHatBinHigh));
//             rescaleHisto1D(hGenEtaPythia);
//             set1DStyle(hGenEtaPythia, 1);

//             // Plot comparisons
//             c->cd(); c->Clear(); setPadStyle();
//             drawEtaPtComparison(c, hRecoEtaEmb, hRecoEtaPythia, 
//                                 jetPtLow, jetPtHigh, -5.2, 5.2, ptHatLow, 
//                                 "Reco (embed)", "Reco (pythia)", collisionSystem, collisionEnergy, 
//                                 false, true, false);
//             c->SaveAs(Form("%s/%s_reco_eta_ptHat_ptHat%d_jetPt_%d_%d.pdf", 
//                       date.Data(), collSystemStr.Data(), (int)ptHatLow, (int)jetPtLow, (int)jetPtHigh));

//             c->cd(); c->Clear(); setPadStyle();
//             drawEtaPtComparison(c, hGenEtaEmb, hGenEtaPythia, 
//                                 jetPtLow, jetPtHigh, -5.2, 5.2, ptHatLow, 
//                                 "Gen (embed)", "Gen (pythia)", collisionSystem, collisionEnergy, 
//                                 false, true, false);
//             c->SaveAs(Form("%s/%s_gen_eta_ptHat_ptHat%d_jetPt_%d_%d.pdf", 
//                       date.Data(), collSystemStr.Data(), (int)ptHatLow, (int)jetPtLow, (int)jetPtHigh));

//         } // for (size_t iJetPt = 0; iJetPt < sizeOfJetPtVals-1; ++iJetPt)
//     } // for (size_t iPtHat = 0; iPtHat < sizeOfPtHatVals-1; ++iPtHat)

//     if (c) { delete c; c = nullptr; }
// }

//________________
// Compare dijet distributions from embedding and Pythia.
// Compare reco, ref, and gen dijet eta distributions.
// void pythia2embeddingDijetComparison(TFile *fEmbedding, TFile *fPythia, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250305") {
//     // Collisions system: 0 = pp, 1 = pPb, 2 = PbPb
//     // energy in TeV

//     TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
//     collSystemStr += Form("%d", (int)(1000 * collisionEnergy));
//     collSystemStr += Form("_emb2pythia");

//     double xTextPosition = 0.35;
//     double yTextPosition = 0.8;
//     TLatex t;
//     t.SetTextFont(42);
//     t.SetTextSize(0.05);

//     // Determine the direction based on the filename
//     TString directionStr;
//     TString fileName = fEmbedding->GetName();
//     if (fileName.Contains("pbgoing", TString::kIgnoreCase)) {
//         directionStr = "Pbgoing";
//     } else if (fileName.Contains("pgoing", TString::kIgnoreCase)) {
//         directionStr = "pgoing";
//     } else {
//         directionStr = "combined";
//     }

//     // Dijet ptAve binning
//     // double dijetPtVals[] { 50.,  60.,   70.,  80.,  90.,
//     //                           100., 110.,  120., 130., 140.,
//     //                           150., 160.,  180., 200., 250., 
//     //                           300., 500. };
//     int dijetPtVals[] { 50,  90,   120,  180, 200, 250, 300, 350, 500 };
//     int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

//     // Dijet dijstributions for embedding
//     TH3D *hRecoDijetPtEtaCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRecoDijetPtEtaCMWeighted"));
//     if (!hRecoDijetPtEtaCMEmb) {
//         std::cerr << "Error: hRecoDijetPtEtaCMWeighted for embedding not found in file." << std::endl;
//         return;
//     }
//     TH3D *hGenDijetPtEtaCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hGenDijetPtEtaCMWeighted"));
//     if (!hGenDijetPtEtaCMEmb) {
//         std::cerr << "Error: hGenDijetPtEtaCMWeighted for embedding not found in file." << std::endl;
//         return;
//     }
//     TH3D *hRefDijetPtEtaCMEmb = dynamic_cast<TH3D *>(fEmbedding->Get("hRefDijetPtEtaCMWeighted"));
//     if (!hRefDijetPtEtaCMEmb) {
//         std::cerr << "Error: hRefDijetPtEtaCMWeighted for embedding not found in file." << std::endl;
//         return;
//     }

//     // Dijet dijstributions for pythia
//     TH3D *hRecoDijetPtEtaCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hRecoDijetPtEtaCMWeighted"));
//     if (!hRecoDijetPtEtaCMPythia) {
//         std::cerr << "Error: hRecoDijetPtEtaCMWeighted for pythia not found in file." << std::endl;
//         return;
//     }
//     TH3D *hGenDijetPtEtaCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hGenDijetPtEtaCMWeighted"));
//     if (!hGenDijetPtEtaCMPythia) {
//         std::cerr << "Error: hGenDijetPtEtaCMWeighted for pythia not found in file." << std::endl;
//         return;
//     }
//     TH3D *hRefDijetPtEtaCMPythia = dynamic_cast<TH3D *>(fPythia->Get("hRefDijetPtEtaCMWeighted"));
//     if (!hRefDijetPtEtaCMPythia) {
//         std::cerr << "Error: hRefDijetPtEtaCMWeighted for pythia not found in file." << std::endl;
//         return;
//     }

//     // Create 1D histograms for dijet eta distributions
//     TH1D *hRecoDijetEtaCMEmb{nullptr};
//     TH1D *hGenDijetEtaCMEmb{nullptr};
//     TH1D *hRefDijetEtaCMEmb{nullptr};

//     TH1D *hRecoDijetEtaCMPythia{nullptr};
//     TH1D *hGenDijetEtaCMPythia{nullptr};
//     TH1D *hRefDijetEtaCMPythia{nullptr};

//     TH1D *hRecoDijetEtaCM_emb2pythia{nullptr};
//     TH1D *hGenDijetEtaCM_emb2pythia{nullptr};
//     TH1D *hRefDijetEtaCM_emb2pythia{nullptr};

//     TCanvas *c = new TCanvas("c", "c", 1000, 1000);
//     setPadStyle();
//     TLegend *leg{nullptr};

//     // Loop over dijet ptAve bins
//     for (int i = 0; i < sizeOfPtVals - 1; ++i) {
//         double ptLow = dijetPtVals[i];
//         double ptHi = dijetPtVals[i + 1];
//         int ptBinLow = hRecoDijetPtEtaCMEmb->GetXaxis()->FindBin(ptLow);
//         int ptBinHi = hRecoDijetPtEtaCMEmb->GetXaxis()->FindBin(ptHi)-1;

//         // Create 1D histograms for embedding
//         hRecoDijetEtaCMEmb = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMEmb->ProjectionY(Form("hRecoDijetEtaCMEmb_%d", i), 
//                                                   ptBinLow, ptBinHi));
//         if (!hRecoDijetEtaCMEmb) {
//             std::cerr << "Error: hRecoDijetEtaCMEmb_" << i << " not found." << std::endl;
//             return;
//         }
//         set1DStyle(hRecoDijetEtaCMEmb, 0);
//         rescaleHisto1D(hRecoDijetEtaCMEmb);
//         hGenDijetEtaCMEmb = dynamic_cast<TH1D *>(hGenDijetPtEtaCMEmb->ProjectionY(Form("hGenDijetEtaCMEmb_%d", i), 
//                                                  ptBinLow, ptBinHi));
//         if (!hGenDijetEtaCMEmb) {
//             std::cerr << "Error: hGenDijetEtaCMEmb_" << i << " not found." << std::endl;
//             return;
//         }
//         set1DStyle(hGenDijetEtaCMEmb, 0);
//         rescaleHisto1D(hGenDijetEtaCMEmb);
//         hRefDijetEtaCMEmb = dynamic_cast<TH1D *>(hRefDijetPtEtaCMEmb->ProjectionY(Form("hRefDijetEtaCMEmb_%d", i),
//                                                   ptBinLow, ptBinHi));
//         if (!hRefDijetEtaCMEmb) {
//             std::cerr << "Error: hRefDijetEtaCMEmb_" << i << " not found." << std::endl;
//             return;
//         }
//         set1DStyle(hRefDijetEtaCMEmb, 0);
//         rescaleHisto1D(hRefDijetEtaCMEmb);

//         // Create 1D histograms for pythia
//         hRecoDijetEtaCMPythia = dynamic_cast<TH1D *>(hRecoDijetPtEtaCMPythia->ProjectionY(Form("hRecoDijetEtaCMPythia_%d", i), 
//                                                      ptBinLow, ptBinHi));
//         if (!hRecoDijetEtaCMPythia) {
//             std::cerr << "Error: hRecoDijetEtaCMPythia_" << i << " not found." << std::endl;
//             return;
//         }
//         set1DStyle(hRecoDijetEtaCMPythia, 1);
//         rescaleHisto1D(hRecoDijetEtaCMPythia);
//         hGenDijetEtaCMPythia = dynamic_cast<TH1D *>(hGenDijetPtEtaCMPythia->ProjectionY(Form("hGenDijetEtaCMPythia_%d", i), 
//                                                   ptBinLow, ptBinHi));
//         if (!hGenDijetEtaCMPythia) {
//             std::cerr << "Error: hGenDijetEtaCMPythia_" << i << " not found." << std::endl;
//             return;
//         }
//         set1DStyle(hGenDijetEtaCMPythia, 1);
//         rescaleHisto1D(hGenDijetEtaCMPythia);
//         hRefDijetEtaCMPythia = dynamic_cast<TH1D *>(hRefDijetPtEtaCMPythia->ProjectionY(Form("hRefDijetEtaCMPythia_%d", i), 
//                                                   ptBinLow, ptBinHi));
//         if (!hRefDijetEtaCMPythia) {
//             std::cerr << "Error: hRefDijetEtaCMPythia_" << i << " not found." << std::endl;
//             return;
//         }
//         set1DStyle(hRefDijetEtaCMPythia, 1);
//         rescaleHisto1D(hRefDijetEtaCMPythia);

//         // Draw and save the comparisons
//         c->cd(); c->Clear(); setPadStyle();
//         drawDijetToGenComparison(c, hRecoDijetEtaCMEmb, hRecoDijetEtaCMPythia, nullptr, nullptr, nullptr,
//                                  (int)ptLow, (int)ptHi, 
//                                  true, false, collisionSystem, collisionEnergy);
//         c->SaveAs(Form("%s/%s_%s_dijetEtaCM_recoComp_ptAve_%d_%d.pdf", date.Data(), 
//                        collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));

//         c->cd(); c->Clear(); setPadStyle();
//         drawDijetToGenComparison(c, hGenDijetEtaCMEmb, hGenDijetEtaCMPythia, nullptr, nullptr, nullptr,
//                                  (int)ptLow, (int)ptHi,
//                                  true, false, collisionSystem, collisionEnergy);
//         c->SaveAs(Form("%s/%s_%s_dijetEtaCM_genComp_ptAve_%d_%d.pdf", date.Data(), 
//                        collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));

//         c->cd(); c->Clear(); setPadStyle();
//         drawDijetToGenComparison(c, hRefDijetEtaCMEmb, hRefDijetEtaCMPythia, nullptr, nullptr, nullptr,
//                                  (int)ptLow, (int)ptHi,
//                                  true, false, collisionSystem, collisionEnergy);
//         c->SaveAs(Form("%s/%s_%s_dijetEtaCM_refComp_ptAve_%d_%d.pdf", date.Data(), 
//                        collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));

//         // Create 1D histograms for embedding to pythia ratios
//         hRecoDijetEtaCM_emb2pythia = dynamic_cast<TH1D *>(hRecoDijetEtaCMEmb->Clone(Form("hRecoDijetEtaCM_emb2pythia_%d", i)));
//         if (!hRecoDijetEtaCM_emb2pythia) {
//             std::cerr << "Error: hRecoDijetEtaCM_emb2pythia_" << i << " not found." << std::endl;
//             return;
//         }
//         hRecoDijetEtaCM_emb2pythia->Divide(hRecoDijetEtaCMEmb, hRecoDijetEtaCMPythia, 1., 1., "b");
//         set1DStyle(hRecoDijetEtaCM_emb2pythia, 0);
        
//         hGenDijetEtaCM_emb2pythia = dynamic_cast<TH1D *>(hGenDijetEtaCMEmb->Clone(Form("hGenDijetEtaCM_emb2pythia_%d", i)));
//         if (!hGenDijetEtaCM_emb2pythia) {
//             std::cerr << "Error: hGenDijetEtaCM_emb2pythia_" << i << " not found." << std::endl;
//             return;
//         }
//         hGenDijetEtaCM_emb2pythia->Divide(hGenDijetEtaCMEmb, hGenDijetEtaCMPythia, 1., 1., "b");
//         set1DStyle(hGenDijetEtaCM_emb2pythia, 4);
        
//         hRefDijetEtaCM_emb2pythia = dynamic_cast<TH1D *>(hRefDijetEtaCMEmb->Clone(Form("hRefDijetEtaCM_emb2pythia_%d", i)));
//         if (!hRefDijetEtaCM_emb2pythia) {
//             std::cerr << "Error: hRefDijetEtaCM_emb2pythia_" << i << " not found." << std::endl;
//             return;


//     } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
//         hRefDijetEtaCM_emb2pythia->Divide(hRefDijetEtaCMEmb, hRefDijetEtaCMPythia, 1., 1., "b");
//         set1DStyle(hRefDijetEtaCM_emb2pythia, 2);

//         // Draw and save the ratios
//         c->cd(); c->Clear(); setPadStyle();
//         hRecoDijetEtaCM_emb2pythia->Draw();
//         hGenDijetEtaCM_emb2pythia->Draw("same");
//         hRefDijetEtaCM_emb2pythia->Draw("same");
//         hRecoDijetEtaCM_emb2pythia->GetXaxis()->SetRangeUser(-2.4, 2.4);
//         hRecoDijetEtaCM_emb2pythia->GetYaxis()->SetRangeUser(0.8, 1.2);
//         hRecoDijetEtaCM_emb2pythia->GetYaxis()->SetTitle("Embedding/Pythia");
//         t.DrawLatexNDC(xTextPosition, yTextPosition, Form("%d < p_{T}^{ave} < %d (GeV)", (int)ptLow, (int)ptHi));
//         leg = new TLegend(0.35, 0.25, 0.65, 0.4);
//         leg->SetTextSize(0.04);
//         leg->SetFillColor(0);
//         leg->SetBorderSize(0);
//         leg->AddEntry(hRecoDijetEtaCM_emb2pythia, "Reco", "p");
//         leg->AddEntry(hGenDijetEtaCM_emb2pythia, "Gen", "p");
//         leg->AddEntry(hRefDijetEtaCM_emb2pythia, "Ref", "p");
//         leg->Draw();
//         c->SaveAs(Form("%s/%s_%s_dijetEtaCM_ratio_ptAve_%d_%d.pdf", 
//                        date.Data(), collSystemStr.Data(), directionStr.Data(), (int)ptLow, (int)ptHi));
//     } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
//     if (c) { delete c; c = nullptr; }
// }


//________________
// Plot dijet distributions for data by projecting those from 2D (ptAve, eta) histogram
void dijetDistributionsForData(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    TString triggerName;
    if (filename.Contains("Jet100", TString::kIgnoreCase)) {
        triggerName = "Jet100";
    } else if (filename.Contains("Jet80", TString::kIgnoreCase)) {
        triggerName = "Jet80";
    } else if (filename.Contains("Jet60", TString::kIgnoreCase)) {
        triggerName = "Jet60";
    } else {
        triggerName = "MB";
    }

    // Dijet ptAve binning
    // int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
    //                          100, 110,  120, 130, 140,
    //                          150, 160,  180, 200, 250, 
    //                          300, 500 };
    double dijetPtNewVals[] {  50., 80., 120., 180., 200., 250., 300., 350., 500. };
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);



    // Retrieve 2D distributions from data
    TH2D *hRecoDijetPtEtaWeighted = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaWeighted"));
    if (!hRecoDijetPtEtaWeighted) {
        std::cerr << "Error: hRecoDijetPtEtaWeighted not found in file." << std::endl;
        return;
    }
    TH2D *hRecoDijetPtEtaCMWeighted = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMWeighted"));
    if (!hRecoDijetPtEtaCMWeighted) {
        std::cerr << "Error: hRecoDijetPtEtaCMWeighted not found in file." << std::endl;
        return;
    }

}

//________________
void usefullDistributions(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    // Dijet ptAve binning
    // int dijetPtNewVals[17] {  50,  60,   70,  80,  90,
    //                          100, 110,  120, 130, 140,
    //                          150, 160,  180, 200, 250, 
    //                          300, 500 };
    int dijetPtNewVals[9] {  50,  90,   120,  180, 200, 250, 300, 350, 500 };
    int sizeOfPtVals = sizeof(dijetPtNewVals)/sizeof(dijetPtNewVals[0]);


}

//________________
void plotSingleJetData2McComparison(TFile *fExp, TFile *fMc, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    // TString directionStr;
    // TString filename = f->GetName();
    // if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
    //     directionStr = "Pbgoing";
    // } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
    //     directionStr = "pgoing";
    // } else {
    //     directionStr = "combined";
    // }
}

//_________________
void plotForwardBackwardRatios(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // Determine the direction based on the filename
    TString directionStr;
    TString filename = f->GetName();
    if (filename.Contains("pbgoing", TString::kIgnoreCase)) {
        directionStr = "Pbgoing";
    } else if (filename.Contains("pgoing", TString::kIgnoreCase)) {
        directionStr = "pgoing";
    } else {
        directionStr = "combined";
    }

    TH2D *hForward2D = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMForwardWeighted"));
    if (!hForward2D) {
        std::cerr << "Error: hRecoDijetPtEtaCMForwardWeighted not found in file." << std::endl;
        return;
    }
    hForward2D->SetName("hForward2D");
    hForward2D->SetDirectory(0);
    TH2D *hBackward2D = dynamic_cast<TH2D *>(f->Get("hRecoDijetPtEtaCMBackwardWeighted"));
    if (!hBackward2D) {
        std::cerr << "Error: hRecoDijetPtEtaCMBackwardWeighted not found in file." << std::endl;
        return;
    }
    hBackward2D->SetName("hBackward2D");
    hBackward2D->SetDirectory(0);


    // int dijetPtVals[] { 50, 80, 100, 140 };
    int dijetPtVals[] { 120, 150, 180, 250, 500 };
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);
    for (int i = 0; i < sizeOfPtVals - 1; ++i) {
        int ptLow = dijetPtVals[i];
        int ptHi = dijetPtVals[i + 1];
        int ptBinLow = hForward2D->GetXaxis()->FindBin(ptLow);
        int ptBinHi = hForward2D->GetXaxis()->FindBin(ptHi) - 1;


        // Project to 1D
        TH1D *hForward1D = dynamic_cast<TH1D *>( hForward2D->ProjectionY(Form("hForward1D_%d_%d", ptLow, ptHi)) );
        set1DStyle(hForward1D, 0, true);
        TH1D *hBackward1D = dynamic_cast<TH1D *>( hBackward2D->ProjectionY(Form("hBackward1D_%d_%d", ptLow, ptHi)) );
        set1DStyle(hBackward1D, 1, true);

        // Create ratio
        TH1D *hFwdBwdRatio = dynamic_cast<TH1D *>(hForward1D->Clone(Form("hFwdBwdRatio_%d_%d", ptLow, ptHi)));
        hFwdBwdRatio->Divide(hForward1D, hBackward1D, 1., 1., "b");
        set1DStyle(hFwdBwdRatio, 2);

        TCanvas *c = new TCanvas("c", "c", 800, 800);
        setPadStyle();
        gPad->SetGrid();
        hFwdBwdRatio->Draw();
        hFwdBwdRatio->GetXaxis()->SetRangeUser(0., 2.);
        hFwdBwdRatio->GetYaxis()->SetRangeUser(0.75, 1.25);
        hFwdBwdRatio->GetYaxis()->SetTitle("Forward / Backward");
        TLatex t;
        t.SetTextFont(42);
        t.SetTextSize(0.05);
        t.DrawLatexNDC(0.3, 0.8, Form("%d < p_{T}^{ave} < %d GeV", ptLow, ptHi));
        plotCMSHeader(collisionSystem, collisionEnergy);
        c->SaveAs(Form("%s/%s_%s_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), directionStr.Data(), ptLow, ptHi));
        delete c;   
    } // for (int i = 0; i < sizeOfPtVals - 1; ++i)
}

//________________
void plotDijetFBRatiosForEta(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TFile *fExpPbgoing = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/Pbgoing/MB_Pbgoing_ak4_jetId_eta19.root") );
    if (!fExpPbgoing || fExpPbgoing->IsZombie()) {
        std::cerr << "Error: Could not open experimental Pb-going file." << std::endl;
        return;
    }
    TFile *fExpPgoing = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/pgoing/MB_pgoing_ak4_jetId_eta19.root") );
    if (!fExpPgoing || fExpPgoing->IsZombie()) {
        std::cerr << "Error: Could not open experimental p-going file." << std::endl;
        return;
    }
    TFile *fMcPbgoing = TFile::Open( "~/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_jetId_eta19.root" );
    if (!fMcPbgoing || fMcPbgoing->IsZombie()) {
        std::cerr << "Error: Could not open MC Pb-going file." << std::endl;
        return;
    }
    TFile *fMcPgoing = TFile::Open( "~/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pgoing_def_ak4_jetId_eta19.root" );
    if (!fMcPgoing || fMcPgoing->IsZombie()) {
        std::cerr << "Error: Could not open MC p-going file." << std::endl;
        return;
    }
    

    TH2D *hExpPbgoingForward2D[6];
    TH2D *hExpPbgoingBackward2D[6];
    TH2D *hExpPgoingForward2D[6];
    TH2D *hExpPgoingBackward2D[6];
    TH2D *hMcPbgoingForward2D[6];
    TH2D *hMcPbgoingBackward2D[6];
    TH2D *hMcPgoingForward2D[6];
    TH2D *hMcPgoingBackward2D[6];

    // Loop over dijet eta_CM cuts |eta_CM|<1.4, 1.5, 1.6, 1.7, 1.8, 1.9
    for (int iEta{0}; iEta<6; iEta++) {

        hExpPbgoingForward2D[iEta] = dynamic_cast<TH2D *>(fExpPbgoing->Get( Form("hRecoDijetPtEtaForwardArr_%d", iEta) ));
        if (!hExpPbgoingForward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << iEta << " not found in experimental Pb-going file." << std::endl;
            return;
        }
        hExpPbgoingForward2D[iEta]->SetName(Form("hExpPbgoingForward2D_%d", iEta));
        hExpPbgoingForward2D[iEta]->SetDirectory(0);
        hExpPbgoingForward2D[iEta]->RebinY(2);
        hExpPbgoingBackward2D[iEta] = dynamic_cast<TH2D *>(fExpPbgoing->Get( Form("hRecoDijetPtEtaBackwardArr_%d", iEta) ));
        if (!hExpPbgoingBackward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << iEta << " not found in experimental Pb-going file." << std::endl;
            return;
        }
        hExpPbgoingBackward2D[iEta]->SetName(Form("hExpPbgoingBackward2D_%d", iEta));
        hExpPbgoingBackward2D[iEta]->SetDirectory(0);
        hExpPbgoingBackward2D[iEta]->RebinY(2);

        hExpPgoingForward2D[iEta] = dynamic_cast<TH2D *>(fExpPgoing->Get( Form("hRecoDijetPtEtaForwardArr_%d", iEta) ));
        if (!hExpPgoingForward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << iEta << " not found in experimental p-going file." << std::endl;
            return;
        }
        hExpPgoingForward2D[iEta]->SetName(Form("hExpPgoingForward2D_%d", iEta));
        hExpPgoingForward2D[iEta]->SetDirectory(0);
        hExpPgoingForward2D[iEta]->RebinY(2);
        hExpPgoingBackward2D[iEta] = dynamic_cast<TH2D *>(fExpPgoing->Get( Form("hRecoDijetPtEtaBackwardArr_%d", iEta) ));
        if (!hExpPgoingBackward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << iEta << " not found in experimental p-going file." << std::endl;
            return;
        }
        hExpPgoingBackward2D[iEta]->SetName(Form("hExpPgoingBackward2D_%d", iEta));
        hExpPgoingBackward2D[iEta]->SetDirectory(0);
        hExpPgoingBackward2D[iEta]->RebinY(2);

        hMcPbgoingForward2D[iEta] = dynamic_cast<TH2D *>(fMcPbgoing->Get( Form("hRecoDijetPtEtaForwardArr_%d", iEta) ));
        if (!hMcPbgoingForward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << iEta << " not found in MC Pb-going file." << std::endl;
            return;
        }
        hMcPbgoingForward2D[iEta]->SetName(Form("hMcPbgoingForward2D_%d", iEta));
        hMcPbgoingForward2D[iEta]->SetDirectory(0);
        hMcPbgoingForward2D[iEta]->RebinY(2);
        hMcPbgoingBackward2D[iEta] = dynamic_cast<TH2D *>(fMcPbgoing->Get( Form("hRecoDijetPtEtaBackwardArr_%d", iEta) ));
        if (!hMcPbgoingBackward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << iEta << " not found in MC Pb-going file." << std::endl;
            return;
        }
        hMcPbgoingBackward2D[iEta]->SetName(Form("hMcPbgoingBackward2D_%d", iEta));
        hMcPbgoingBackward2D[iEta]->SetDirectory(0);
        hMcPbgoingBackward2D[iEta]->RebinY(2);

        hMcPgoingForward2D[iEta] = dynamic_cast<TH2D *>(fMcPgoing->Get( Form("hRecoDijetPtEtaForwardArr_%d", iEta) ));
        if (!hMcPgoingForward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << iEta << " not found in MC p-going file." << std::endl;
            return;
        }
        hMcPgoingForward2D[iEta]->SetName(Form("hMcPgoingForward2D_%d", iEta));
        hMcPgoingForward2D[iEta]->SetDirectory(0);
        hMcPgoingForward2D[iEta]->RebinY(2);
        hMcPgoingBackward2D[iEta] = dynamic_cast<TH2D *>(fMcPgoing->Get( Form("hRecoDijetPtEtaBackwardArr_%d", iEta) ));
        if (!hMcPgoingBackward2D[iEta]) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << iEta << " not found in MC p-going file." << std::endl;
            return;
        }
        hMcPgoingBackward2D[iEta]->SetName(Form("hMcPgoingBackward2D_%d", iEta));
        hMcPgoingBackward2D[iEta]->SetDirectory(0);
        hMcPgoingBackward2D[iEta]->RebinY(2);
    }

    fExpPbgoing->Close();
    fExpPgoing->Close();
    fMcPbgoing->Close();
    fMcPgoing->Close();

    TH1D *hExpPbgoingForward1D[6];
    TH1D *hExpPbgoingBackward1D[6];
    TH1D *hExpPgoingForward1D[6];
    TH1D *hExpPgoingBackward1D[6];
    TH1D *hMcPbgoingForward1D[6];
    TH1D *hMcPbgoingBackward1D[6];
    TH1D *hMcPgoingForward1D[6];
    TH1D *hMcPgoingBackward1D[6];

    TH1D *hExpPbgoingFwdBwdRatio1D[6];
    TH1D *hExpPgoingFwdBwdRatio1D[6];
    TH1D *hMcPbgoingFwdBwdRatio1D[6];
    TH1D *hMcPgoingFwdBwdRatio1D[6];

    TH1D *hExpPbgoingOverPgoingFwdBwdRatio1D[6];
    TH1D *hMcPbgoingOverPgoingFwdBwdRatio1D[6];

    double dijetPtVals[] { 60., 80. };
    int sizeOfPtVals = sizeof( dijetPtVals )/sizeof(dijetPtVals[0]);
    static constexpr std::array<const char*, 7> p8Colors {"kBlue", "kRed", "kMagenta", "kP8Orange", "kP8Green", "kP8Azure", "kBlack"};

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    TLegend *leg{nullptr};
    TString directionStr;

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    for (int iPt{0}; iPt<sizeOfPtVals-1; iPt++) {

        double ptLow = dijetPtVals[iPt];
        double ptHi = dijetPtVals[iPt+1];
        int ptBinLow = hExpPbgoingForward2D[0]->GetXaxis()->FindBin(ptLow);
        int ptBinHi = hExpPbgoingForward2D[0]->GetXaxis()->FindBin(ptHi) - 1;
        
        for (int iEta{0}; iEta<6; iEta++) {
            // Project to 1D
            hExpPbgoingForward1D[iEta] = dynamic_cast<TH1D *>( hExpPbgoingForward2D[iEta]->ProjectionY(Form("hExpPbgoingForward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hExpPbgoingForward1D[iEta], 2);
            hExpPbgoingBackward1D[iEta] = dynamic_cast<TH1D *>( hExpPbgoingBackward2D[iEta]->ProjectionY(Form("hExpPbgoingBackward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hExpPbgoingBackward1D[iEta], 2);
            hExpPgoingForward1D[iEta] = dynamic_cast<TH1D *>( hExpPgoingForward2D[iEta]->ProjectionY(Form("hExpPgoingForward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hExpPgoingForward1D[iEta], 2);
            hExpPgoingBackward1D[iEta] = dynamic_cast<TH1D *>( hExpPgoingBackward2D[iEta]->ProjectionY(Form("hExpPgoingBackward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hExpPgoingBackward1D[iEta], 2);
            hMcPbgoingForward1D[iEta] = dynamic_cast<TH1D *>( hMcPbgoingForward2D[iEta]->ProjectionY(Form("hMcPbgoingForward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hMcPbgoingForward1D[iEta], 1);
            hMcPbgoingBackward1D[iEta] = dynamic_cast<TH1D *>( hMcPbgoingBackward2D[iEta]->ProjectionY(Form("hMcPbgoingBackward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hMcPbgoingBackward1D[iEta], 1);
            hMcPgoingForward1D[iEta] = dynamic_cast<TH1D *>( hMcPgoingForward2D[iEta]->ProjectionY(Form("hMcPgoingForward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hMcPgoingForward1D[iEta], 1);
            hMcPgoingBackward1D[iEta] = dynamic_cast<TH1D *>( hMcPgoingBackward2D[iEta]->ProjectionY(Form("hMcPgoingBackward1D_%d", iEta), ptBinLow, ptBinHi) );
            set1DStyle(hMcPgoingBackward1D[iEta], 1);

            hExpPbgoingFwdBwdRatio1D[iEta] = dynamic_cast<TH1D *>(hExpPbgoingForward1D[iEta]->Clone(Form("hExpPbgoingFwdBwdRatio1D_%d", iEta)));
            hExpPbgoingFwdBwdRatio1D[iEta]->Divide(hExpPbgoingForward1D[iEta], hExpPbgoingBackward1D[iEta], 1., 1.);
            hExpPbgoingFwdBwdRatio1D[iEta]->SetLineColor( p8Colors[iEta] );
            hExpPbgoingFwdBwdRatio1D[iEta]->SetMarkerColor( p8Colors[iEta] );
            
            hExpPgoingFwdBwdRatio1D[iEta] = dynamic_cast<TH1D *>(hExpPgoingForward1D[iEta]->Clone(Form("hExpPgoingFwdBwdRatio1D_%d", iEta)));
            hExpPgoingFwdBwdRatio1D[iEta]->Divide(hExpPgoingForward1D[iEta], hExpPgoingBackward1D[iEta], 1., 1.);
            hExpPgoingFwdBwdRatio1D[iEta]->SetLineColor( p8Colors[iEta] );
            hExpPgoingFwdBwdRatio1D[iEta]->SetMarkerColor( p8Colors[iEta] );

            hMcPbgoingFwdBwdRatio1D[iEta] = dynamic_cast<TH1D *>(hMcPbgoingForward1D[iEta]->Clone(Form("hMcPbgoingFwdBwdRatio1D_%d", iEta)));
            hMcPbgoingFwdBwdRatio1D[iEta]->Divide(hMcPbgoingForward1D[iEta], hMcPbgoingBackward1D[iEta], 1., 1.);
            hMcPbgoingFwdBwdRatio1D[iEta]->SetLineColor( p8Colors[iEta] );
            hMcPbgoingFwdBwdRatio1D[iEta]->SetMarkerColor( p8Colors[iEta] );

            hMcPgoingFwdBwdRatio1D[iEta] = dynamic_cast<TH1D *>(hMcPgoingForward1D[iEta]->Clone(Form("hMcPgoingFwdBwdRatio1D_%d", iEta)));
            hMcPgoingFwdBwdRatio1D[iEta]->Divide(hMcPgoingForward1D[iEta], hMcPgoingBackward1D[iEta], 1., 1.);
            hMcPgoingFwdBwdRatio1D[iEta]->SetLineColor( p8Colors[iEta] );
            hMcPgoingFwdBwdRatio1D[iEta]->SetMarkerColor( p8Colors[iEta] );

            // Compare Pb-going over p-going for the ratios
            hExpPbgoingOverPgoingFwdBwdRatio1D[iEta] = dynamic_cast<TH1D *>(hExpPbgoingFwdBwdRatio1D[iEta]->Clone(Form("hExpPbgoingOverPgoingFwdBwdRatio1D_%d", iEta)));
            hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->Divide(hExpPbgoingFwdBwdRatio1D[iEta], hExpPgoingFwdBwdRatio1D[iEta], 1., 1.);

            hMcPbgoingOverPgoingFwdBwdRatio1D[iEta] = dynamic_cast<TH1D *>(hMcPbgoingFwdBwdRatio1D[iEta]->Clone(Form("hMcPbgoingOverPgoingFwdBwdRatio1D_%d", iEta)));
            hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->Divide(hMcPbgoingFwdBwdRatio1D[iEta], hMcPgoingFwdBwdRatio1D[iEta], 1., 1.);
        } // for (int iEta{0}; iEta<6; iEta++)



        //
        // Draw the ratios for different |eta_CM|<X cuts on the same canvas for Pb-going data
        //
        c->cd(); 
        setPadStyle();
        gPad->SetGrid();
        for (int iEta{0}; iEta<6; iEta++) {

            hExpPbgoingFwdBwdRatio1D[iEta]->Draw( (iEta == 0) ? "" : "same" );
            hExpPbgoingFwdBwdRatio1D[iEta]->SetLineColor(p8Colors[iEta]);
            hExpPbgoingFwdBwdRatio1D[iEta]->SetMarkerColor(p8Colors[iEta]);

            if (iEta == 0) {
                hExpPbgoingFwdBwdRatio1D[iEta]->GetXaxis()->SetRangeUser(0., 2.);
                hExpPbgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetRangeUser(0.7, 1.3);
                hExpPbgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetTitle("Forward / Backward");
                leg = new TLegend(0.16, 0.16, 0.35, 0.45);
                leg->SetTextSize(0.04);
                leg->SetFillColor(0);
                leg->SetBorderSize(0);
                leg->Draw();
            }
            leg->AddEntry(hExpPbgoingFwdBwdRatio1D[iEta], Form("|#eta_{CM}| < %.1f", 1.4 + iEta * 0.1), "p");
        } // for (int iEta{0}; iEta<6; iEta++)

        leg->Draw();
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.85, Form("%d < p_{T}^{ave} < %d (GeV)", (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]) );
        t.DrawLatexNDC(0.4, 0.8, "Pb-going MB data");
        c->SaveAs(Form("%s/%s_exp_Pbgoing_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]));

        
        //
        // Draw the ratios for different |eta_CM|<X cuts on the same canvas for p-going data
        //
        c->cd();
        for (int iEta{0}; iEta<6; iEta++) {

            hExpPgoingFwdBwdRatio1D[iEta]->Draw( (iEta == 0) ? "" : "same" );
            hExpPgoingFwdBwdRatio1D[iEta]->SetLineColor(p8Colors[iEta]);
            hExpPgoingFwdBwdRatio1D[iEta]->SetMarkerColor(p8Colors[iEta]);

            if (iEta == 0) {
                hExpPgoingFwdBwdRatio1D[iEta]->GetXaxis()->SetRangeUser(0., 2.);
                hExpPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetRangeUser(0.7, 1.3);
                hExpPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetTitle("Forward / Backward");
                leg = new TLegend(0.16, 0.16, 0.35, 0.45);
                leg->SetTextSize(0.04);
                leg->SetFillColor(0);
                leg->SetBorderSize(0);
                leg->Draw();
            }
            leg->AddEntry(hExpPgoingFwdBwdRatio1D[iEta], Form("|#eta_{CM}| < %.1f", 1.4 + iEta * 0.1), "p");
        } // for (int iEta{0}; iEta<6; iEta++)

        leg->Draw();
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.85, Form("%d < p_{T}^{ave} < %d (GeV)", (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]) );
        t.DrawLatexNDC(0.4, 0.8, "p-going MB data");
        c->SaveAs(Form("%s/%s_exp_pgoing_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]));

        //
        // Ratios Pb-going/p-going data for |eta_CM|<X cuts on the same canvas
        //
        c->cd();
        for (int iEta{0}; iEta<6; iEta++) {

            hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->Draw( (iEta == 0) ? "" : "same" );
            hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->SetLineColor(p8Colors[iEta]);
            hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->SetMarkerColor(p8Colors[iEta]);

            if (iEta == 0) {
                hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->GetXaxis()->SetRangeUser(0., 2.);
                hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetRangeUser(0.8, 1.2);
                hExpPbgoingOverPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetTitle("Pb-going / p-going");
                leg = new TLegend(0.16, 0.16, 0.35, 0.45);
                leg->SetTextSize(0.04);
                leg->SetFillColor(0);
                leg->SetBorderSize(0);
                leg->Draw();
            }
            leg->AddEntry(hExpPbgoingOverPgoingFwdBwdRatio1D[iEta], Form("|#eta_{CM}| < %.1f", 1.4 + iEta * 0.1), "p");
        } // for (int iEta{0}; iEta<6; iEta++)
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.85, Form("%d < p_{T}^{ave} < %d (GeV)", (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]) );
        t.DrawLatexNDC(0.4, 0.8, "Pb-going / p-going MB data");
        c->SaveAs(Form("%s/%s_exp_PbgoingOverPgoing_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]));

        //
        // Draw the ratios for different |eta_CM|<X cuts on the same canvas for Pb-going embedding (MC)
        //
        c->cd(); 
        setPadStyle();
        gPad->SetGrid();
        for (int iEta{0}; iEta<6; iEta++) {

            hMcPbgoingFwdBwdRatio1D[iEta]->Draw( (iEta == 0) ? "" : "same" );
            hMcPbgoingFwdBwdRatio1D[iEta]->SetLineColor(p8Colors[iEta]);
            hMcPbgoingFwdBwdRatio1D[iEta]->SetMarkerColor(p8Colors[iEta]);

            if (iEta == 0) {
                hMcPbgoingFwdBwdRatio1D[iEta]->GetXaxis()->SetRangeUser(0., 2.);
                hMcPbgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetRangeUser(0.7, 1.3);
                hMcPbgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetTitle("Forward / Backward");
                leg = new TLegend(0.16, 0.16, 0.35, 0.45);
                leg->SetTextSize(0.04);
                leg->SetFillColor(0);
                leg->SetBorderSize(0);
                leg->Draw();
            }
            leg->AddEntry(hMcPbgoingFwdBwdRatio1D[iEta], Form("|#eta_{CM}| < %.1f", 1.4 + iEta * 0.1), "p");
        } // for (int iEta{0}; iEta<6; iEta++)

        leg->Draw();
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.85, Form("%d < p_{T}^{ave} < %d (GeV)", (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]) );
        t.DrawLatexNDC(0.4, 0.8, "Pb-going embedding");
        c->SaveAs(Form("%s/%s_mc_Pbgoing_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]));


        //
        // Draw the ratios for different |eta_CM|<X cuts on the same canvas for p-going embedding (MC)
        //
        c->cd(); 
        setPadStyle();
        gPad->SetGrid();
        for (int iEta{0}; iEta<6; iEta++) {

            hMcPgoingFwdBwdRatio1D[iEta]->Draw( (iEta == 0) ? "" : "same" );
            hMcPgoingFwdBwdRatio1D[iEta]->SetLineColor(p8Colors[iEta]);
            hMcPgoingFwdBwdRatio1D[iEta]->SetMarkerColor(p8Colors[iEta]);

            if (iEta == 0) {
                hMcPgoingFwdBwdRatio1D[iEta]->GetXaxis()->SetRangeUser(0., 2.);
                hMcPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetRangeUser(0.7, 1.3);
                hMcPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetTitle("Forward / Backward");
                leg = new TLegend(0.16, 0.16, 0.35, 0.45);
                leg->SetTextSize(0.04);
                leg->SetFillColor(0);
                leg->SetBorderSize(0);
                leg->Draw();
            }
            leg->AddEntry(hMcPgoingFwdBwdRatio1D[iEta], Form("|#eta_{CM}| < %.1f", 1.4 + iEta * 0.1), "p");
        } // for (int iEta{0}; iEta<6; iEta++)

        leg->Draw();
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.85, Form("%d < p_{T}^{ave} < %d (GeV)", (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]) );
        t.DrawLatexNDC(0.4, 0.8, "p-going embedding");
        c->SaveAs(Form("%s/%s_mc_pgoing_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                       date.Data(), collSystemStr.Data(), (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]));

        //
        // Ratios Pb-going/p-going embedding for |eta_CM|<X cuts on the same canvas
        //
        c->cd();
        for (int iEta{0}; iEta<6; iEta++) {

            hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->Draw( (iEta == 0) ? "" : "same" );
            hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->SetLineColor(p8Colors[iEta]);
            hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->SetMarkerColor(p8Colors[iEta]);

            if (iEta == 0) {
                hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->GetXaxis()->SetRangeUser(0., 2.);
                hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetRangeUser(0.8, 1.2);
                hMcPbgoingOverPgoingFwdBwdRatio1D[iEta]->GetYaxis()->SetTitle("Pb-going / p-going");
                leg = new TLegend(0.16, 0.16, 0.35, 0.45);
                leg->SetTextSize(0.04);
                leg->SetFillColor(0);
                leg->SetBorderSize(0);
                leg->Draw();
            }
            leg->AddEntry(hMcPbgoingOverPgoingFwdBwdRatio1D[iEta], Form("|#eta_{CM}| < %.1f", 1.4 + iEta * 0.1), "p");
        } // for (int iEta{0}; iEta<6; iEta++)
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.4, 0.85, Form("%d < p_{T}^{ave} < %d (GeV)", (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]) );
        t.DrawLatexNDC(0.4, 0.8, "Pb-going / p-going MC embedding");
        c->SaveAs(Form("%s/%s_mc_PbgoingOverPgoing_dijetEtaCM_forwardBackwardRatio_ptAve_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), (int)dijetPtVals[iPt], (int)dijetPtVals[iPt+1]));
    } // for (int iPt{0}; iPt<sizeOfPtVals-1; iPt++)

    // Loop over dijet eta_CM cuts |eta_CM|<1.4, 1.5, 1.6, 1.7, 1.8, 1.9

}

//________________
void plotXjComparison(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TString residualStr = "";
    // residualStr = "_noResidual";

    TFile *fExpPbgoing = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/Pbgoing/MB_Pbgoing_ak4_jetId_eta19%s.root", residualStr.Data()) );
    if (!fExpPbgoing || fExpPbgoing->IsZombie()) {
        std::cerr << "Error: Could not open experimental Pb-going file." << std::endl;
        return;
    }
    // TFile *fExpPbgoing = TFile::Open( "~/cernbox/ana/pPb8160/exp/Pbgoing/MB_Pbgoing_ak4_jetId_eta19_noResidual.root" );
    // if (!fExpPbgoing || fExpPbgoing->IsZombie()) {
    //     std::cerr << "Error: Could not open experimental Pb-going noResidual file." << std::endl;
    //     return;
    // }
    TFile *fExppgoing = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/pgoing/MB_pgoing_ak4_jetId_eta19%s.root", residualStr.Data()) );
    if (!fExppgoing || fExppgoing->IsZombie()) {
        std::cerr << "Error: Could not open experimental p-going file." << std::endl;
        return;
    }
    // TFile *fExppgoing = TFile::Open( "~/cernbox/ana/pPb8160/exp/pgoing/MB_pgoing_ak4_jetId_eta19_noResidual.root" );
    // if (!fExppgoing || fExppgoing->IsZombie()) {
    //     std::cerr << "Error: Could not open experimental p-going noResidual file." << std::endl;
    //     return;
    // }

    TFile *fMcPbgoing = TFile::Open( "~/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_jetId_eta19.root" );
    if (!fMcPbgoing || fMcPbgoing->IsZombie()) {
        std::cerr << "Error: Could not open MC Pb-going file." << std::endl;
        return;
    }
    TFile *fMcpgoing = TFile::Open( "~/cernbox/ana/pPb8160/embedding/pgoing/oEmbedding_pgoing_def_ak4_jetId_eta19.root" );
    if (!fMcpgoing || fMcpgoing->IsZombie()) {
        std::cerr << "Error: Could not open MC p-going file." << std::endl;
        return;
    }


    TH1D *hExpPbgoingXj[3];
    TH1D *hExppgoingXj[3];
    TH1D *hMcPbgoingXj[3];
    TH1D *hMcpgoingXj[3];

    // TH1D *hExpPbgoingXjCM[3];
    // TH1D *hExppgoingXjCM[3];
    // TH1D *hMcPbgoingXjCM[3];
    // TH1D *hMcpgoingXjCM[3];

    TH1D *hExpForwardXjRatio;
    TH1D *hMcForwardXjRatio;
    TH1D *hExpForwardXjCMRatio;
    TH1D *hMcForwardXjCMRatio;

    // Retrieve xJ distributions
    for (int i{0}; i<3; i++) {

        // Laboratory frame

        hExpPbgoingXj[i] = dynamic_cast<TH1D *>(fExpPbgoing->Get( Form("hRecoDijetXj_%d", i) ));
        if ( !hExpPbgoingXj[i] ) {
            std::cerr << "Error: hRecoDijetXj_" << i << " not found in experimental Pb-going file." << std::endl;
            return;
        }
        set1DStyle( hExpPbgoingXj[i], 0 );
        hExpPbgoingXj[i]->SetName( Form("hExpPbgoingXj_%d", i) );
        hExpPbgoingXj[i]->SetDirectory(0);
        hExpPbgoingXj[i]->Scale( 1.0 / hExpPbgoingXj[i]->Integral() );

        hExppgoingXj[i] = dynamic_cast<TH1D *>(fExppgoing->Get( Form("hRecoDijetXj_%d", i) ));
        if ( !hExppgoingXj[i] ) {
            std::cerr << "Error: hRecoDijetXj_" << i << " not found in experimental p-going file." << std::endl;
            return;
        }
        set1DStyle( hExppgoingXj[i], 1 );
        hExppgoingXj[i]->SetName( Form("hExppgoingXj_%d", i) );
        hExppgoingXj[i]->SetDirectory(0);
        hExppgoingXj[i]->Scale( 1.0 / hExppgoingXj[i]->Integral() );

        hMcPbgoingXj[i] = dynamic_cast<TH1D *>(fMcPbgoing->Get( Form("hRecoDijetXj_%d", i) ));
        if ( !hMcPbgoingXj[i] ) {
            std::cerr << "Error: hRecoDijetXj_" << i << " not found in MC Pb-going file." << std::endl;
            return;
        }
        set1DStyle( hMcPbgoingXj[i], 1 );
        hMcPbgoingXj[i]->SetName( Form("hMcPbgoingXj_%d", i) );
        hMcPbgoingXj[i]->SetDirectory(0);
        hMcPbgoingXj[i]->Scale( 1.0 / hMcPbgoingXj[i]->Integral() );

        hMcpgoingXj[i] = dynamic_cast<TH1D *>(fMcpgoing->Get( Form("hRecoDijetXj_%d", i) ));
        if ( !hMcpgoingXj[i] ) {
            std::cerr << "Error: hRecoDijetXj_" << i << " not found in MC p-going file." << std::endl;
            return;
        }
        set1DStyle( hMcpgoingXj[i], 0 );
        hMcpgoingXj[i]->SetName( Form("hMcpgoingXj_%d", i) );
        hMcpgoingXj[i]->SetDirectory(0);
        hMcpgoingXj[i]->Scale( 1.0 / hMcpgoingXj[i]->Integral() );

        // Center-of-mass frame

        // hExpPbgoingXjCM[i] = dynamic_cast<TH1D *>(fExpPbgoing->Get( Form("hRecoDijetXjCM_%d", i) ));
        // if ( !hExpPbgoingXjCM[i] ) {
        //     std::cerr << "Error: hRecoDijetXjCM_" << i << " not found in experimental Pb-going file." << std::endl;
        //     return;
        // }
        // set1DStyle( hExpPbgoingXjCM[i], 0 );
        // hExpPbgoingXjCM[i]->SetName( Form("hExpPbgoingXjCM_%d", i) );
        // hExpPbgoingXjCM[i]->SetDirectory(0);
        // hExpPbgoingXjCM[i]->Scale( 1.0 / hExpPbgoingXjCM[i]->Integral() );

        // hExppgoingXjCM[i] = dynamic_cast<TH1D *>(fExppgoing->Get( Form("hRecoDijetXjCM_%d", i) ));
        // if ( !hExppgoingXjCM[i] ) {
        //     std::cerr << "Error: hRecoDijetXjCM_" << i << " not found in experimental p-going file." << std::endl;
        //     return;
        // }
        // set1DStyle( hExppgoingXjCM[i], 1 );
        // hExppgoingXjCM[i]->SetName( Form("hExppgoingXjCM_%d", i) );
        // hExppgoingXjCM[i]->SetDirectory(0);
        // hExppgoingXjCM[i]->Scale( 1.0   / hExppgoingXjCM[i]->Integral() );
        // hMcPbgoingXjCM[i] = dynamic_cast<TH1D *>(fMcPbgoing->Get( Form("hRecoDijetXjCM_%d", i) ));
    } // for (int i{0}; i<3; i++)

    // Make sure that histograms do not belong to files anymore
    fExpPbgoing->Close();
    fExppgoing->Close();
    fMcPbgoing->Close();
    fMcpgoing->Close();

    // Make ratios of forward (Pb-going) to backward (p-going) in data (Pb-> +eta, p-> -eta)
    hExpForwardXjRatio = dynamic_cast<TH1D *>( hExpPbgoingXj[2]->Clone("hExpForwardXjRatio") );
    hExpForwardXjRatio->Divide( hExpForwardXjRatio, hExppgoingXj[2], 1.0, 1.0 );
    // set1DStyle( hExpForwardXjRatio, 2 );
    // Make ratios of forward (Pb-going) to backward (p-going) in MC (Pb-> -eta, p-> +eta)
    hMcForwardXjRatio = dynamic_cast<TH1D *>( hMcPbgoingXj[0]->Clone("hMcForwardXjRatio") );
    hMcForwardXjRatio->Divide( hMcForwardXjRatio, hMcpgoingXj[0], 1.0, 1.0 );
    // set1DStyle( hMcForwardXjRatio, 2 );

    TLegend *leg{nullptr};
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.03);

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    c->SetGrid();


    // Plot comparisons of xJ distributions for data
    // hExpPbgoingXj[2]->Draw();
    // hExppgoingXj[2]->Draw("same");
    // hExpPbgoingXj[2]->GetXaxis()->SetRangeUser(0.3, 1.);
    // hExpPbgoingXj[2]->GetYaxis()->SetRangeUser(0.01, 0.14);
    // hExpPbgoingXj[2]->GetYaxis()->SetTitle("dN/dx_{J}");
    // plotCMSHeader(collisionSystem, collisionEnergy);
    // t.DrawLatexNDC(0.2, 0.85, "Data");
    // t.DrawLatexNDC(0.2, 0.8, "p_{T}^{Lead} > 50 GeV");
    // t.DrawLatexNDC(0.2, 0.75, "p_{T}^{SubLead} > 40 GeV");
    // t.DrawLatexNDC(0.2, 0.7, "50 < p_{T}^{ave} < 90 GeV");
    // t.DrawLatexNDC(0.2, 0.65, "1.8 < #eta^{Lead} < 2.4");
    // t.DrawLatexNDC(0.2, 0.6, "|#eta^{SubLead}| < 1.2");
    // leg = new TLegend(0.2, 0.45, 0.35, 0.55);
    // leg->SetTextSize(0.03);
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->AddEntry(hExpPbgoingXj[2], "Pb-going (1.8 < #eta^{Lead} < 2.4)", "p");
    // leg->AddEntry(hExppgoingXj[2], "p-going (1.8 < #eta^{Lead} < 2.4)", "p");
    // leg->Draw();
    // c->SaveAs(Form("%s/%s_dijetXj_exp_Pbgoing_vs_pgoing_lab_comparison.pdf", 
    //                date.Data(), collSystemStr.Data()) );
    
    // Plot comparisons of xJ distributions for MC
    // hMcPbgoingXj[0]->Draw();
    // hMcpgoingXj[0]->Draw("same");
    // hMcPbgoingXj[0]->GetXaxis()->SetRangeUser(0.3, 1.);
    // hMcPbgoingXj[0]->GetYaxis()->SetRangeUser(0.01, 0.14);
    // hMcPbgoingXj[0]->GetYaxis()->SetTitle("dN/dx_{J}");
    // plotCMSHeader(collisionSystem, collisionEnergy);
    // t.DrawLatexNDC(0.2, 0.85, "MC (Embedding)");
    // t.DrawLatexNDC(0.2, 0.8, "p_{T}^{Lead} > 50 GeV");
    // t.DrawLatexNDC(0.2, 0.75, "p_{T}^{SubLead} > 40 GeV");
    // t.DrawLatexNDC(0.2, 0.7, "50 < p_{T}^{ave} < 90 GeV");
    // t.DrawLatexNDC(0.2, 0.65, "1.8 < #eta^{Lead} < 2.4");
    // t.DrawLatexNDC(0.2, 0.6, "|#eta^{SubLead}| < 1.2");
    // leg = new TLegend(0.2, 0.45, 0.35, 0.55);
    // leg->SetTextSize(0.03);
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->AddEntry(hMcPbgoingXj[0], "Pb-going (1.8 < #eta^{Lead} < 2.4)", "p");
    // leg->AddEntry(hMcpgoingXj[0], "p-going (1.8 < #eta^{Lead} < 2.4)", "p");
    // leg->Draw();
    // c->SaveAs(Form("%s/%s_dijetXj_mc_Pbgoing_vs_pgoing_lab_comparison.pdf", 
    //                date.Data(), collSystemStr.Data()) );

    // Plot ratios of forward to backward xJ 
    hExpForwardXjRatio->Draw();
    hExpForwardXjRatio->GetXaxis()->SetRangeUser(0.3, 1.);
    hExpForwardXjRatio->GetYaxis()->SetRangeUser(0.75, 1.25);
    hExpForwardXjRatio->GetYaxis()->SetTitle("Pb-going / p-going");
    hMcForwardXjRatio->Draw("same");
    plotCMSHeader(collisionSystem, collisionEnergy);
    // t.DrawLatexNDC(0.3, 0.8, "Pb-going / p-going (in forward)");
    leg = new TLegend(0.5, 0.75, 0.6, 0.85);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hExpForwardXjRatio, "Data", "p");
    leg->AddEntry(hMcForwardXjRatio, "Embedding", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_dijetXj_forwardBackwardRatio_lab_comparison.pdf", 
                   date.Data(), collSystemStr.Data()) );

    // Plot comparisons of xJ distributions for data in CM
    // hExpPbgoingXjCM[2]->Draw();
    // hExppgoingXjCM[2]->Draw("same");
    // hExpPbgoingXjCM[2]->GetXaxis()->SetRangeUser(0.3, 1.);
    // hExpPbgoingXjCM[2]->GetYaxis()->SetRangeUser(0.01, 0.14);
    // hExpPbgoingXjCM[2]->GetYaxis()->SetTitle("dN/dx_{J}");
    // plotCMSHeader(collisionSystem, collisionEnergy);
    // t.SetTextFont(42);
    // t.SetTextSize(0.03);
    // //t.DrawLatexNDC(0.3, 0.8, "Reco dijets, Pb-going vs p-going");
    // t.DrawLatexNDC(0.2, 0.85, "Data");
    // t.DrawLatexNDC(0.2, 0.8, "p_{T}^{Lead} > 50 GeV");
    // t.DrawLatexNDC(0.2, 0.75, "p_{T}^{SubLead} > 40 GeV");
    // t.DrawLatexNDC(0.2, 0.7, "50 < p_{T}^{ave} < 90 GeV");
    // t.DrawLatexNDC(0.2, 0.65, "1.6 < #eta^{Lead}_{CM} < 1.9");
    // t.DrawLatexNDC(0.2, 0.6, "|#eta^{SubLead}_{CM}| < 1.2");
    // leg = new TLegend(0.2, 0.45, 0.35, 0.55);
    // leg->SetTextSize(0.03);
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->AddEntry(hExpPbgoingXjCM[2], "Pb-going (1.6 < #eta^{Lead}_{CM} < 1.9)", "p");
    // leg->AddEntry(hExppgoingXjCM[0], "p-going (-1.9 < #eta^{Lead}_{CM} < -1.6)", "p");
    // leg->Draw();
    // c->SaveAs(Form("%s/%s_dijetXj_exp_Pbgoing_vs_pgoing_cm_comparison.pdf", 
    //                date.Data(), collSystemStr.Data()) );

    // Plot comparisons of xJ distributions for MC in CM
    // hMcPbgoingXj[0]->Draw();
    // hMcpgoingXj[0]->Draw("same");
    // hMcPbgoingXj[0]->GetXaxis()->SetRangeUser(0.3, 1.);
    // hMcPbgoingXj[0]->GetYaxis()->SetRangeUser(0.01, 0.14);
    // hMcPbgoingXj[0]->GetYaxis()->SetTitle("dN/dx_{J}");
    // plotCMSHeader(collisionSystem, collisionEnergy);
    // t.DrawLatexNDC(0.2, 0.85, "MC (Embedding)");
    // t.DrawLatexNDC(0.2, 0.8, "p_{T}^{Lead} > 50 GeV");
    // t.DrawLatexNDC(0.2, 0.75, "p_{T}^{SubLead} > 40 GeV");
    // t.DrawLatexNDC(0.2, 0.7, "50 < p_{T}^{ave} < 90 GeV");
    // t.DrawLatexNDC(0.2, 0.65, "1.6 < #eta^{Lead}_{CM} < 1.9");
    // t.DrawLatexNDC(0.2, 0.6, "|#eta^{SubLead}_{CM}| < 1.2");
    // leg = new TLegend(0.2, 0.45, 0.35, 0.55);
    // leg->SetTextSize(0.03);
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->AddEntry(hMcPbgoingXj[0], "Pb-going (-1.9 < #eta^{Lead}_{CM} < -1.6)", "p");
    // leg->AddEntry(hMcpgoingXj[2], "p-going (1.6 < #eta^{Lead}_{CM} < 1.9)", "p");
    // leg->Draw();
    // c->SaveAs(Form("%s/%s_dijetXj_mc_Pbgoing_vs_pgoing_cm_comparison.pdf", 
    //                date.Data(), collSystemStr.Data()) );

}

//________________
void compareSingleJetEta(int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TFile *fExpNoResidual = TFile::Open( "~/cernbox/ana/pPb8160/exp/Pbgoing/MB_Pbgoing_ak4_jetId_eta19_noResidual.root" );
    if (!fExpNoResidual || fExpNoResidual->IsZombie()) {
        std::cerr << "Error: Could not open experimental no residual." << std::endl;
        return;
    }

    TFile *fExpAll = TFile::Open( "~/cernbox/ana/pPb8160/exp/Pbgoing/MB_pgoing_ak4_jetId_eta19.root" );
    if (!fExpAll || fExpAll->IsZombie()) {
        std::cerr << "Error: Could not open experimental all corrections." << std::endl;
        return;
    }


    TFile *fMcAll = TFile::Open( "~/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_jetId_eta19.root" );
    if (!fMcAll || fMcAll->IsZombie()) {
        std::cerr << "Error: Could not open MC all corrections." << std::endl;
        return;
    }

    TH2D *hExpPbgoingJetPtEtaRaw = dynamic_cast<TH2D *>(fExpNoResidual->Get("hRecoInclusiveAllJetPtRawEtaStdBins"));
    TH2D *hExpPbgoingJetPtEtaNoResidual = dynamic_cast<TH2D *>(fExpNoResidual->Get("hRecoInclusiveAllJetPtEtaStdBins"));
    TH2D *hExpPbgoingJetPtEtaAll = dynamic_cast<TH2D *>(fExpAll->Get("hRecoInclusiveAllJetPtEtaStdBins"));

    hExpPbgoingJetPtEtaRaw->SetDirectory(0);
    hExpPbgoingJetPtEtaNoResidual->SetDirectory(0);
    hExpPbgoingJetPtEtaAll->SetDirectory(0);

    TH2D *hMcPbgoingJetPtEtaRaw = dynamic_cast<TH2D *>(fMcAll->Get("hRecoInclusiveAllJetPtRawEtaStdBins"));
    TH2D *hMcPbgoingJetPtEtaAll = dynamic_cast<TH2D *>(fMcAll->Get("hRecoInclusiveAllJetPtEtaStdBins"));
    hMcPbgoingJetPtEtaRaw->SetDirectory(0);
    hMcPbgoingJetPtEtaAll->SetDirectory(0);

    fExpNoResidual->Close();
    fExpAll->Close();
    fMcAll->Close();

    int jetPtVals[] { 50, 90 };
    int sizeOfPtVals = sizeof(jetPtVals)/sizeof(jetPtVals[0]);

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2,1);
    TLegend *leg{nullptr};
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    // Loop over pT bins
    for (int i{0}; i<sizeOfPtVals; i++) {

        int ptLow = jetPtVals[i];
        int ptHi = jetPtVals[i + 1];
        int ptBinLow = hExpPbgoingJetPtEtaRaw->GetYaxis()->FindBin(ptLow);
        int ptBinHi = hExpPbgoingJetPtEtaRaw->GetYaxis()->FindBin(ptHi) - 1;

        TH1D *hExpJetEtaRaw = dynamic_cast<TH1D *>( hExpPbgoingJetPtEtaRaw->ProjectionX( Form("hExpJetEtaRaw_%d_%d", ptLow, ptHi), ptBinLow, ptBinHi ) );
        set1DStyle( hExpJetEtaRaw, 0 );
        TH1D *hExpJetEtaNoResidual = dynamic_cast<TH1D *>( hExpPbgoingJetPtEtaNoResidual->ProjectionX( Form("hExpJetEtaNoResidual_%d_%d", ptLow, ptHi), ptBinLow, ptBinHi ) );
        set1DStyle( hExpJetEtaNoResidual, 1 );
        TH1D *hExpJetEtaAll = dynamic_cast<TH1D *>( hExpPbgoingJetPtEtaAll->ProjectionX( Form("hExpJetEtaAll_%d_%d", ptLow, ptHi), ptBinLow, ptBinHi ) );
        set1DStyle( hExpJetEtaAll, 2 );

        TH1D *hExpNoResidualRatio = dynamic_cast<TH1D *>( hExpJetEtaNoResidual->Clone( Form("hExpNoResidualRatio_%d_%d", ptLow, ptHi) ) );
        hExpNoResidualRatio->Divide( hExpJetEtaNoResidual, hExpJetEtaRaw, 1., 1.);
        TH1D *hExpAllRatio = dynamic_cast<TH1D *>( hExpJetEtaAll->Clone( Form("hExpAllRatio_%d_%d", ptLow, ptHi) ) );
        hExpAllRatio->Divide( hExpJetEtaAll, hExpJetEtaRaw, 1., 1. );

        TH1D *hMcJetEtaRaw = dynamic_cast<TH1D *>( hMcPbgoingJetPtEtaRaw->ProjectionX( Form("hMcJetEtaRaw_%d_%d", ptLow, ptHi), ptBinLow, ptBinHi ) );
        set1DStyle( hMcJetEtaRaw, 0 );
        TH1D *hMcJetEtaAll = dynamic_cast<TH1D *>( hMcPbgoingJetPtEtaAll->ProjectionX( Form("hMcJetEtaAll_%d_%d", ptLow, ptHi), ptBinLow, ptBinHi ) );
        set1DStyle( hMcJetEtaAll, 1 );
        
        TH1D *hMcAllRatio = dynamic_cast<TH1D *>( hMcJetEtaAll->Clone( Form("hMcAllRatio_%d_%d", ptLow, ptHi) ) );
        hMcAllRatio->Divide( hMcJetEtaAll, hMcJetEtaRaw, 1., 1. );


    }

}

//________________
void plotEtaDistributionsForRunId(int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TFile *fMBFile = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4_jetId_eta19.root") );

    TH1D *hInclusiveJetEta[6];
    TH1D *hLeadJetEta[6];
    TH1D *hSubLeadJetEta[6];
    TH1D *hDijetEtaCM[6];

    TH1D *hInclusiveJetEtaRatioToAll[5];
    TH1D *hLeadJetEtaRatioToAll[5];
    TH1D *hSubLeadJetEtaRatioToAll[5];
    TH1D *hDijetEtaCMRatioToAll[5];

    TString colors[6] {"kBlack", "kRed", "kBlue", "kMagenta", "kGreen", "kCyan"};
    TString runId[6] {"All", "285480", "285505", "285517", "285832", "285993"};

    for (int iRun{0}; iRun<6; iRun++) {
        hInclusiveJetEta[iRun] = dynamic_cast<TH1D *>( fMBFile->Get( Form("hRecoInclusiveJetEtaRun_%d", iRun) ) );
        if ( !hInclusiveJetEta[iRun] ) {
            std::cerr << "Error: hRecoInclusiveJetEtaRun_" << iRun << " not found in file." << std::endl;
            return;
        }
        set1DStyle( hInclusiveJetEta[iRun], 2 );
        hInclusiveJetEta[iRun]->Scale( 1./hInclusiveJetEta[iRun]->Integral() );
        hInclusiveJetEta[iRun]->SetLineColor( colors[iRun].Data() );
        hInclusiveJetEta[iRun]->SetMarkerColor( colors[iRun].Data() );
        hInclusiveJetEta[iRun]->SetName( Form("hInclusiveJetEta_%d", iRun) );
        hInclusiveJetEta[iRun]->SetDirectory(0);

        if (iRun != 0) {
            hInclusiveJetEtaRatioToAll[iRun - 1] = dynamic_cast<TH1D *>( hInclusiveJetEta[iRun]->Clone( Form("hInclusiveJetEtaRatioToAll_%d", iRun-1) ) );
            hInclusiveJetEtaRatioToAll[iRun - 1]->Divide( hInclusiveJetEtaRatioToAll[iRun - 1], hInclusiveJetEta[0], 1., 1. );
            set1DStyle( hInclusiveJetEtaRatioToAll[iRun - 1], 2 );
            hInclusiveJetEtaRatioToAll[iRun - 1]->SetLineColor( colors[iRun].Data() );
            hInclusiveJetEtaRatioToAll[iRun - 1]->SetMarkerColor( colors[iRun].Data() );
        }

        hLeadJetEta[iRun] = dynamic_cast<TH1D *>( fMBFile->Get( Form("hRecoLeadJetEtaRun_%d", iRun) ) );
        if ( !hLeadJetEta[iRun] ) {
            std::cerr << "Error: hRecoLeadJetEtaRun_" << iRun << " not found in file." << std::endl;
            return;
        }
        set1DStyle( hLeadJetEta[iRun], 2, true );
        hLeadJetEta[iRun]->SetLineColor( colors[iRun].Data() );
        hLeadJetEta[iRun]->SetMarkerColor( colors[iRun].Data() );
        hLeadJetEta[iRun]->SetName( Form("hLeadJetEta_%d", iRun) );
        hLeadJetEta[iRun]->SetDirectory(0);

        hSubLeadJetEta[iRun] = dynamic_cast<TH1D *>( fMBFile->Get( Form("hRecoSubLeadJetEtaRun_%d", iRun) ) );
        if ( !hSubLeadJetEta[iRun] ) {
            std::cerr << "Error: hRecoSubLeadJetEtaRun_" << iRun << " not found in file." << std::endl;
            return;     
        }
        set1DStyle( hSubLeadJetEta[iRun], 2, true );
        hSubLeadJetEta[iRun]->SetLineColor( colors[iRun].Data() );
        hSubLeadJetEta[iRun]->SetMarkerColor( colors[iRun].Data() );
        hSubLeadJetEta[iRun]->SetName( Form("hSubLeadJetEta_%d", iRun) );
        hSubLeadJetEta[iRun]->SetDirectory(0);

        hDijetEtaCM[iRun] = dynamic_cast<TH1D *>( fMBFile->Get( Form("hRecoDijetEtaCMRun_%d", iRun) ) );
        if ( !hDijetEtaCM[iRun] ) {
            std::cerr << "Error: hRecoDijetEtaCMRun_" << iRun << " not found in file." << std::endl;
            return;
        }
        set1DStyle( hDijetEtaCM[iRun], 2, true );
        hDijetEtaCM[iRun]->SetLineColor( colors[iRun].Data() );
        hDijetEtaCM[iRun]->SetMarkerColor( colors[iRun].Data() );
        hDijetEtaCM[iRun]->SetName( Form("hDijetEtaCM_%d", iRun) );
        hDijetEtaCM[iRun]->SetDirectory(0);
    }

    fMBFile->Close();

    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    c->Divide(2,1);
    c->cd(1);
    setPadStyle();
    gPad->SetGrid();
    c->cd(2);
    setPadStyle();
    gPad->SetGrid();

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);
    
    TLegend *leg{nullptr};

    // Loop over runs
    c->cd(1);
    for (int iRun{0}; iRun<6; iRun++) {
        hInclusiveJetEta[iRun]->Draw( (iRun==0) ? "" : "same" );
        if (iRun == 0) {
            hInclusiveJetEta[iRun]->GetXaxis()->SetTitle("#eta^{Inclusive}");
            hInclusiveJetEta[iRun]->GetYaxis()->SetTitle("dN/d#eta");
            hInclusiveJetEta[iRun]->GetYaxis()->SetRangeUser(0., 0.25);
            leg = new TLegend(0.7, 0.7, 0.85, 0.85);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hInclusiveJetEta[iRun], Form("Run %s", runId[iRun].Data()), "p");
    } // for (int i{0}; i<6; i++)
    t.DrawLatexNDC(0.2, 0.85, "Inclusive jets");
    t.DrawLatexNDC(0.2, 0.8, "40 < p_{T} < 90 GeV");
    leg->Draw();

    c->cd(2);
    for (int iRun{0}; iRun<5; iRun++) {
        hLeadJetEta[iRun]->Draw( (iRun==0) ? "" : "same" );
        if (iRun == 0) {
            hLeadJetEta[iRun]->GetXaxis()->SetTitle("#eta^{Lead}");
            hLeadJetEta[iRun]->GetYaxis()->SetTitle("dN/d#eta");
            hLeadJetEta[iRun]->GetYaxis()->SetRangeUser(0.5, 1.5);
            leg = new TLegend(0.7, 0.7, 0.85, 0.85);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hLeadJetEta[iRun], Form("Run %s", runId[iRun].Data()), "p");
    } // for (int i{0}; i<6; i++)
}

//________________
void plotNpdfToDataComparison(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    int rebinFactor = 2;
    static constexpr std::array<int, 6> etaCMIntCut = {19, 18, 17, 16, 15, 14};
    static constexpr std::array<double, 6> etaCMDoubleCut = {1.9, 1.8, 1.7, 1.6, 1.5, 1.4};
    static constexpr std::array<const char*, 6> p8Colors {"kBlack", "kBlue", "kRed", "kMagenta", "kP8Orange", "kP8Green"};
    const int etaCMCutSize = etaCMIntCut.size();
    const double ptLow = 60.;
    const double ptHigh = 80.;
    const double dijetEtaCMMax = 3.2;
    const int dijetEtaCMBins = 64;

    // nPDF and PDF histograms
    TH1D *hNpdfEtaCMForward[etaCMCutSize];
    TH1D *hNpdfEtaCMBackward[etaCMCutSize];
    TH1D *hNpdfEtaCMForwardBackwardRatio[etaCMCutSize];
    TH1D *hNpdfEtaCMFull[etaCMCutSize];

    TH1D *hPdfEtaCMForward[etaCMCutSize];
    TH1D *hPdfEtaCMBackward[etaCMCutSize];
    TH1D *hPdfEtaCMForwardBackwardRatio[etaCMCutSize];
    TH1D *hPdfEtaCMFull[etaCMCutSize];

    TH1D *hNpdfToPdfEtaCMFullRatio[etaCMCutSize];
    TH1D *hNpdfToPdfEtaCMForwardBackwardRatio[etaCMCutSize];

    // Data histograms
    TH2D *hData2DEtaCMForward[etaCMCutSize];
    TH2D *hData2DEtaCMBackward[etaCMCutSize];

    TH1D *hDataEtaCMForward[etaCMCutSize];
    TH1D *hDataEtaCMBackward[etaCMCutSize];
    TH1D *hDataEtaCMForwardBackwardRatio[etaCMCutSize];
    TH1D *hDataEtaCMFull[etaCMCutSize];

    // Data to nPDF ratios
    TH1D *hDataToNpdfForwardBackwardRatio[etaCMCutSize];
    TH1D *hDataToNpdfFullEtaRatio[etaCMCutSize];

    TFile *fNpdf = TFile::Open( Form("./npdf/epps21_pPb8160_dijet_eta_dset.root" ) );
    if ( !fNpdf || fNpdf->IsZombie() ) {
        std::cerr << "Error: Could not open npdf file." << std::endl;
        return;
    }

    // TFile *fData = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4_jetId_eta19.root") );
    TFile *fData = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4_jetId_eta19_newResidual.root") );
    if ( !fData || fData->IsZombie() ) {
        std::cerr << "Error: Could not open data file." << std::endl;
        return;
    }

    // Read histograms from files
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {

        //
        // nPDF
        //

        // Full
        hNpdfEtaCMFull[iEta] = dynamic_cast<TH1D *>( fNpdf->Get( Form("pPb_hDijetEtaCMFull_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hNpdfEtaCMFull[iEta] ) {
            std::cerr << "Error: pPb_hDijetEtaCMFull_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hNpdfEtaCMFull[iEta]->SetDirectory(0);
        set1DStyle( hNpdfEtaCMFull[iEta], 0, true );
        hNpdfEtaCMFull[iEta]->SetName( Form("hNpdfEtaCMFull_%d", iEta) );
        hNpdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->Rebin(rebinFactor);

        // Forward
        hNpdfEtaCMForward[iEta] = dynamic_cast<TH1D *>( fNpdf->Get( Form("pPb_hDijetEtaCMForward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hNpdfEtaCMForward[iEta] ) {
            std::cerr << "Error: pPb_hDijetEtaCMForward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hNpdfEtaCMForward[iEta]->SetDirectory(0);
        set1DStyle( hNpdfEtaCMForward[iEta], 0 );
        hNpdfEtaCMForward[iEta]->SetName( Form("hNpdfEtaCMForward_%d", iEta) );
        hNpdfEtaCMForward[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMForward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMForward[iEta]->Rebin(rebinFactor);

        // Backward
        hNpdfEtaCMBackward[iEta] = dynamic_cast<TH1D *>( fNpdf->Get( Form("pPb_hDijetEtaCMBackward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hNpdfEtaCMBackward[iEta] ) {
            std::cerr << "Error: pPb_hDijetEtaCMBackward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hNpdfEtaCMBackward[iEta]->SetDirectory(0);
        set1DStyle( hNpdfEtaCMBackward[iEta], 0 );
        hNpdfEtaCMBackward[iEta]->SetName( Form("hNpdfEtaCMBackward_%d", iEta) );
        hNpdfEtaCMBackward[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMBackward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMBackward[iEta]->Rebin(rebinFactor);

        // Forward/backward ratio
        hNpdfEtaCMForwardBackwardRatio[iEta] = (TH1D*)hNpdfEtaCMForward[iEta]->Clone(Form("hNpdfEtaCMForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hNpdfEtaCMForwardBackwardRatio[iEta]->Divide(hNpdfEtaCMBackward[iEta]);
        hNpdfEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

        //
        // PDF
        //

        // Full
        hPdfEtaCMFull[iEta] = dynamic_cast<TH1D *>( fNpdf->Get( Form("pp_hDijetEtaCMFull_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfEtaCMFull[iEta] ) {
            std::cerr << "Error: pp_hDijetEtaCMFull_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfEtaCMFull[iEta]->SetDirectory(0);
        set1DStyle( hPdfEtaCMFull[iEta], 1, true );
        hPdfEtaCMFull[iEta]->SetName( Form("hPdfEtaCMFull_%d", iEta) );
        hPdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMFull[iEta]->Rebin(rebinFactor);

        // Forward
        hPdfEtaCMForward[iEta] = dynamic_cast<TH1D *>( fNpdf->Get( Form("pp_hDijetEtaCMForward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfEtaCMForward[iEta] ) {
            std::cerr << "Error: pp_hDijetEtaCMForward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfEtaCMForward[iEta]->SetDirectory(0);
        set1DStyle( hPdfEtaCMForward[iEta], 1 );
        hPdfEtaCMForward[iEta]->SetName( Form("hPdfEtaCMForward_%d", iEta) );
        hPdfEtaCMForward[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMForward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMForward[iEta]->Rebin(rebinFactor);

        // Backward
        hPdfEtaCMBackward[iEta] = dynamic_cast<TH1D *>( fNpdf->Get( Form("pp_hDijetEtaCMBackward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfEtaCMBackward[iEta] ) {
            std::cerr << "Error: pp_hDijetEtaCMBackward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfEtaCMBackward[iEta]->SetDirectory(0);
        set1DStyle( hPdfEtaCMBackward[iEta], 1 );
        hPdfEtaCMBackward[iEta]->SetName( Form("hPdfEtaCMBackward_%d", iEta) );
        hPdfEtaCMBackward[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMBackward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMBackward[iEta]->Rebin(rebinFactor);

        // Forward/backward ratio
        hPdfEtaCMForwardBackwardRatio[iEta] = (TH1D*)hPdfEtaCMForward[iEta]->Clone(Form("hPdfEtaCMForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hPdfEtaCMForwardBackwardRatio[iEta]->Divide(hPdfEtaCMBackward[iEta]);
        hPdfEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

        // Ratio of nPDF to PDF
        hNpdfToPdfEtaCMFullRatio[iEta] = (TH1D*)hNpdfEtaCMFull[iEta]->Clone(Form("hNpdfToPdfEtaCMFullRatio_%d_0", etaCMIntCut[iEta]));
        hNpdfToPdfEtaCMFullRatio[iEta]->Divide(hPdfEtaCMFull[iEta]);
        hNpdfToPdfEtaCMFullRatio[iEta]->SetDirectory(0);

        hNpdfToPdfEtaCMForwardBackwardRatio[iEta] = (TH1D*)hNpdfEtaCMForwardBackwardRatio[iEta]->Clone(Form("hNpdfToPdfEtaCMForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->Divide(hPdfEtaCMForwardBackwardRatio[iEta]);
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

        recalculateUncertaintyOfNpdf( hNpdfToPdfEtaCMForwardBackwardRatio[iEta] );

        //
        // Data
        //

        // Forward
        hData2DEtaCMForward[iEta] = dynamic_cast<TH2D *>( fData->Get( Form("hRecoDijetPtEtaForwardArr_%d", etaCMCutSize - iEta - 1) ) );
        if ( !hData2DEtaCMForward[iEta] ) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << etaCMCutSize - iEta - 1 << " not found in data file." << std::endl;
            return;
        }
        hData2DEtaCMForward[iEta]->SetDirectory(0);
        hData2DEtaCMForward[iEta]->RebinY(rebinFactor);
        int ptBinLow = hData2DEtaCMForward[iEta]->GetXaxis()->FindBin(ptLow);
        int ptBinHigh = hData2DEtaCMForward[iEta]->GetXaxis()->FindBin(ptHigh) - 1;
        hDataEtaCMForward[iEta] = dynamic_cast<TH1D *>( hData2DEtaCMForward[iEta]->ProjectionY( Form("hDataEtaCMForward_%d", iEta), ptBinLow, ptBinHigh ) );
        set1DStyle( hDataEtaCMForward[iEta], 0 );
        hDataEtaCMForward[iEta]->SetName( Form("hDataEtaCMForward_%d", iEta) );
        hDataEtaCMForward[iEta]->SetLineColor( p8Colors[iEta] );
        hDataEtaCMForward[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Backward
        hData2DEtaCMBackward[iEta] = dynamic_cast<TH2D *>( fData->Get( Form("hRecoDijetPtEtaBackwardArr_%d", etaCMCutSize - iEta - 1) ) );
        if ( !hData2DEtaCMBackward[iEta] ) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << etaCMCutSize - iEta - 1 << " not found in data file." << std::endl;
            return;
        }
        hData2DEtaCMBackward[iEta]->SetDirectory(0);
        hData2DEtaCMBackward[iEta]->RebinY(rebinFactor);
        hDataEtaCMBackward[iEta] = dynamic_cast<TH1D *>( hData2DEtaCMBackward[iEta]->ProjectionY( Form("hDataEtaCMBackward_%d", iEta), ptBinLow, ptBinHigh ) );
        set1DStyle( hDataEtaCMBackward[iEta], 0 );
        hDataEtaCMBackward[iEta]->SetName( Form("hDataEtaCMBackward_%d", iEta) );
        hDataEtaCMBackward[iEta]->SetLineColor( p8Colors[iEta] );
        hDataEtaCMBackward[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Make full eta distribution by combining forward and backward
        hDataEtaCMFull[iEta] = new TH1D( Form("hDataEtaCMFull_%d", iEta), 
                                         Form("hDataEtaCMFull_%d;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", iEta), 
                                         20, -2.0, 2.0 );
        hDataEtaCMFull[iEta]->Sumw2();
        hDataEtaCMFull[iEta]->SetDirectory(0);
        makeFullEtaFromForwardAndBackward( hDataEtaCMFull[iEta], hDataEtaCMForward[iEta], hDataEtaCMBackward[iEta] );
        set1DStyle( hDataEtaCMFull[iEta], 0, true );
        hDataEtaCMFull[iEta]->SetName( Form("hDataEtaCMFull_%d", iEta) );
        hDataEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hDataEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Ratio
        hDataEtaCMBackward[iEta]->SetName( Form("hDataEtaCMBackward_%d", iEta) );
        hDataEtaCMForwardBackwardRatio[iEta] = (TH1D*)hDataEtaCMForward[iEta]->Clone(Form("hDataEtaCMForwardBackwardRatio_%d", iEta));
        hDataEtaCMForwardBackwardRatio[iEta]->Divide(hDataEtaCMBackward[iEta]);
        hDataEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

        //
        // Data to nPDF ratios
        //

        // Forward/backward ratio
        hDataToNpdfForwardBackwardRatio[iEta] = (TH1D*)hDataEtaCMForwardBackwardRatio[iEta]->Clone(Form("hDataToNpdfForwardBackwardRatio_%d", iEta));
        hDataToNpdfForwardBackwardRatio[iEta]->Divide(hNpdfToPdfEtaCMForwardBackwardRatio[iEta]);
        hDataToNpdfForwardBackwardRatio[iEta]->SetDirectory(0);

        // Full eta distribution
        hDataToNpdfFullEtaRatio[iEta] = (TH1D*)hDataEtaCMFull[iEta]->Clone(Form("hDataToNpdfFullEtaRatio_%d", iEta));
        hDataToNpdfFullEtaRatio[iEta]->Divide(hNpdfEtaCMFull[iEta]);
        hDataToNpdfFullEtaRatio[iEta]->SetDirectory(0);
    } // for (int i{0}; i<etaCMCutSize; i++)

    fData->Close();
    fNpdf->Close();

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TLegend *leg{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();

    // Plot forward/backward ratios for data on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataEtaCMForwardBackwardRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hDataEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
            hDataEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.3);
            hDataEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
            leg = new TLegend(0.16, 0.5, 0.34, 0.8);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hDataEtaCMForwardBackwardRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int i{0}; i<etaCMCutSize; i++)
     plotCMSHeader(collisionSystem, collisionEnergy);
     t.DrawLatexNDC(0.4, 0.85, "MB data");
     t.DrawLatexNDC(0.4, 0.8, "60 < p_{T}^{ave} < 80 GeV");
     leg->Draw();
     c->SaveAs(Form("%s/MB_%s_forwardBackwardRatio_etaCM_comparison.pdf", 
                    date.Data(), collSystemStr.Data()) );

    // Plot full eta comparison of data to nPDF for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {

        hDataEtaCMFull[iEta]->SetLineColor( kRed );
        hDataEtaCMFull[iEta]->SetMarkerColor( kRed );
        hNpdfEtaCMFull[iEta]->SetLineColor( kBlue );
        hNpdfEtaCMFull[iEta]->SetMarkerColor( kBlue );

        hDataEtaCMFull[iEta]->Draw();
        hNpdfEtaCMFull[iEta]->Draw("same");
        hDataEtaCMFull[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
        hDataEtaCMFull[iEta]->GetYaxis()->SetRangeUser(0., 0.16);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "MB to nPDF comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.6, 0.34, 0.7);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hDataEtaCMFull[iEta], "MB data", "p");
        leg->AddEntry(hNpdfEtaCMFull[iEta], "nPDF (EPPS21)", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_data2nPDF_fullEtaRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );

        hDataEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hDataEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int i{0}; i<etaCMCutSize; i++)


    TGraphErrors *grNpdfToPdfForwardBackwardRatio[etaCMCutSize];
    for (int iEta{0}; iEta<etaCMCutSize; ++iEta) {
        int nPoints = hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetNbinsX();
        grNpdfToPdfForwardBackwardRatio[iEta] = new TGraphErrors(nPoints);
        for (int iPoint{0}; iPoint<nPoints; ++iPoint) {
            double x = hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetBinCenter(iPoint + 1);
            double y = hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetBinContent(iPoint + 1);
            double ex = 0.;
            double ey = hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetBinError(iPoint + 1);
            grNpdfToPdfForwardBackwardRatio[iEta]->SetPoint(iPoint, x, y);
            grNpdfToPdfForwardBackwardRatio[iEta]->SetPointError(iPoint, ex, ey);
        }
        grNpdfToPdfForwardBackwardRatio[iEta]->SetName( Form("grNpdfToPdfForwardBackwardRatio_%d", iEta) );
        grNpdfToPdfForwardBackwardRatio[iEta]->SetLineColor( kBlue );
        grNpdfToPdfForwardBackwardRatio[iEta]->SetMarkerColor( kBlue );
        grNpdfToPdfForwardBackwardRatio[iEta]->SetFillColorAlpha( kBlue, 0.35 );
        grNpdfToPdfForwardBackwardRatio[iEta]->SetFillStyle(1001);
    } //    for (int iEta{0}; iEta<etaCMCutSize; ++iEta) {

    // Plot forward/backward ratio comparison of data to nPDF for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataEtaCMForwardBackwardRatio[iEta]->SetLineColor( kRed );
        hDataEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( kRed );
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->SetLineColor( kBlue );
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( kBlue );

        grNpdfToPdfForwardBackwardRatio[iEta]->Draw("AL2");
        hDataEtaCMForwardBackwardRatio[iEta]->Draw("same");
        grNpdfToPdfForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
        grNpdfToPdfForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.1);
        grNpdfToPdfForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
        // hDataEtaCMForwardBackwardRatio[iEta]->Draw();
        // hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->Draw("E1 SAME");
        hDataEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
        hDataEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.1);
        hDataEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "MB to nPDF comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.25, 0.34, 0.35);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hDataEtaCMForwardBackwardRatio[iEta], "MB data", "p");
        leg->AddEntry(hNpdfToPdfEtaCMForwardBackwardRatio[iEta], "nPDF (EPPS21)", "A");
        leg->Draw();
        c->SaveAs(Form("%s/%s_data2nPDF_forwardBackwardRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );

        hDataEtaCMForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hDataEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)

    // Compilation of ratios of data to nPDF for different eta cuts on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataToNpdfFullEtaRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hDataToNpdfFullEtaRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hDataToNpdfFullEtaRatio[iEta]->GetYaxis()->SetTitle("Data / nPDF");
            hDataToNpdfFullEtaRatio[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
            hDataToNpdfFullEtaRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.4);
            leg = new TLegend(0.35, 0.5, 0.65, 0.75);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hDataToNpdfFullEtaRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "MB to nPDF comparison");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_data2nPDF_fullEtaRatio_etaCM_ratios.pdf", 
                    date.Data(), collSystemStr.Data()) );
}

//________________
void plotNpdfToPdfComparison(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TFile *f = TFile::Open( Form("./npdf/epps21_pPb8160_dijet_eta_dset.root" ) );
    if ( !f || f->IsZombie() ) {
        std::cerr << "Error: Could not open npdf file." << std::endl;
        return;
    }

    int rebinFactor = 2;
    static constexpr std::array<int, 7> etaCMIntCut = {25, 19, 18, 17, 16, 15, 14};
    static constexpr std::array<double, 7> etaCMDoubleCut = {2.5, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4};
    static constexpr std::array<const char*, 7> p8Colors {"kBlack", "kBlue", "kRed", "kMagenta", "kP8Orange", "kP8Green", "kP8Azure"};
    const int etaCMCutSize = etaCMIntCut.size();
    const double ptLow = 60.;
    const double ptHigh = 80.;

    TH1D *hNpdfEtaCMForward[etaCMCutSize];
    TH1D *hNpdfEtaCMBackward[etaCMCutSize];
    TH1D *hNpdfEtaCMForwardBackwardRatio[etaCMCutSize];
    TH1D *hNpdfEtaCMFull[etaCMCutSize];

    TH1D *hPdfEtaCMForward[etaCMCutSize];
    TH1D *hPdfEtaCMBackward[etaCMCutSize];
    TH1D *hPdfEtaCMForwardBackwardRatio[etaCMCutSize];
    TH1D *hPdfEtaCMFull[etaCMCutSize];

    TH1D *hNpdfToPdfEtaCMFullRatio[etaCMCutSize];
    TH1D *hNpdfToPdfEtaCMForwardBackwardRatio[etaCMCutSize];
    TH1D *hNpdfToPdfEtaCMForwardBackwardRatioTo25[etaCMCutSize - 1];

    // Read histograms from files and make comparisons
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        //
        // nPDF
        //

        // Full
        hNpdfEtaCMFull[iEta] = dynamic_cast<TH1D *>( f->Get( Form("pPb_hDijetEtaCMFull_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hNpdfEtaCMFull[iEta] ) {
            std::cerr << "Error: pPb_hDijetEtaCMFull_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hNpdfEtaCMFull[iEta]->SetDirectory(0);
        set1DStyle( hNpdfEtaCMFull[iEta], 0, true );
        hNpdfEtaCMFull[iEta]->SetName( Form("hNpdfEtaCMFull_%d", iEta) );
        hNpdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->Rebin(rebinFactor);

        // Forward
        hNpdfEtaCMForward[iEta] = dynamic_cast<TH1D *>( f->Get( Form("pPb_hDijetEtaCMForward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hNpdfEtaCMForward[iEta] ) {
            std::cerr << "Error: pPb_hDijetEtaCMForward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hNpdfEtaCMForward[iEta]->SetDirectory(0);
        set1DStyle( hNpdfEtaCMForward[iEta], 0 );
        hNpdfEtaCMForward[iEta]->SetName( Form("hNpdfEtaCMForward_%d", iEta) );
        hNpdfEtaCMForward[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMForward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMForward[iEta]->Rebin(rebinFactor);

        // Backward
        hNpdfEtaCMBackward[iEta] = dynamic_cast<TH1D *>( f->Get( Form("pPb_hDijetEtaCMBackward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hNpdfEtaCMBackward[iEta] ) {
            std::cerr << "Error: pPb_hDijetEtaCMBackward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hNpdfEtaCMBackward[iEta]->SetDirectory(0);
        set1DStyle( hNpdfEtaCMBackward[iEta], 0 );
        hNpdfEtaCMBackward[iEta]->SetName( Form("hNpdfEtaCMBackward_%d", iEta) );
        hNpdfEtaCMBackward[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMBackward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hNpdfEtaCMBackward[iEta]->Rebin(rebinFactor);

        // Forward/backward ratio
        hNpdfEtaCMForwardBackwardRatio[iEta] = (TH1D*)hNpdfEtaCMForward[iEta]->Clone(Form("hNpdfEtaCMForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hNpdfEtaCMForwardBackwardRatio[iEta]->Divide(hNpdfEtaCMBackward[iEta]);
        hNpdfEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

        //
        // PDF
        //

        // Full
        hPdfEtaCMFull[iEta] = dynamic_cast<TH1D *>( f->Get( Form("pp_hDijetEtaCMFull_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfEtaCMFull[iEta] ) {
            std::cerr << "Error: pp_hDijetEtaCMFull_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfEtaCMFull[iEta]->SetDirectory(0);
        set1DStyle( hPdfEtaCMFull[iEta], 1, true );
        hPdfEtaCMFull[iEta]->SetName( Form("hPdfEtaCMFull_%d", iEta) );
        hPdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMFull[iEta]->Rebin(rebinFactor);

        // Forward
        hPdfEtaCMForward[iEta] = dynamic_cast<TH1D *>( f->Get( Form("pp_hDijetEtaCMForward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfEtaCMForward[iEta] ) {
            std::cerr << "Error: pp_hDijetEtaCMForward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfEtaCMForward[iEta]->SetDirectory(0);
        set1DStyle( hPdfEtaCMForward[iEta], 1 );
        hPdfEtaCMForward[iEta]->SetName( Form("hPdfEtaCMForward_%d", iEta) );
        hPdfEtaCMForward[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMForward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMForward[iEta]->Rebin(rebinFactor);

        // Backward
        hPdfEtaCMBackward[iEta] = dynamic_cast<TH1D *>( f->Get( Form("pp_hDijetEtaCMBackward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfEtaCMBackward[iEta] ) {
            std::cerr << "Error: pp_hDijetEtaCMBackward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfEtaCMBackward[iEta]->SetDirectory(0);
        set1DStyle( hPdfEtaCMBackward[iEta], 1 );
        hPdfEtaCMBackward[iEta]->SetName( Form("hPdfEtaCMBackward_%d", iEta) );
        hPdfEtaCMBackward[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMBackward[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMBackward[iEta]->Rebin(rebinFactor);

        // Forward/backward ratio
        hPdfEtaCMForwardBackwardRatio[iEta] = (TH1D*)hPdfEtaCMForward[iEta]->Clone(Form("hPdfEtaCMForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hPdfEtaCMForwardBackwardRatio[iEta]->Divide(hPdfEtaCMBackward[iEta]);
        hPdfEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

        // Ratio of nPDF to PDF
        hNpdfToPdfEtaCMFullRatio[iEta] = (TH1D*)hNpdfEtaCMFull[iEta]->Clone(Form("hNpdfToPdfEtaCMFullRatio_%d_0", etaCMIntCut[iEta]));
        hNpdfToPdfEtaCMFullRatio[iEta]->Divide(hPdfEtaCMFull[iEta]);
        hNpdfToPdfEtaCMFullRatio[iEta]->SetDirectory(0);

        hNpdfToPdfEtaCMForwardBackwardRatio[iEta] = (TH1D*)hNpdfEtaCMForwardBackwardRatio[iEta]->Clone(Form("hNpdfToPdfEtaCMForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->Divide(hPdfEtaCMForwardBackwardRatio[iEta]);
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);

    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)

    f->Close();

    // Plotting components
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TLegend *leg{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();

    // Plot forward/backward ratio comparison of nPDF and PDF for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hNpdfEtaCMForwardBackwardRatio[iEta]->SetLineColor( kRed );
        hNpdfEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( kRed );
        hPdfEtaCMForwardBackwardRatio[iEta]->SetLineColor( kBlue );
        hPdfEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( kBlue );
 
        hNpdfEtaCMForwardBackwardRatio[iEta]->Draw();
        hPdfEtaCMForwardBackwardRatio[iEta]->Draw("same");
        hNpdfEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.6);
        hNpdfEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.1);
        hNpdfEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "nPDF to PDF comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.25, 0.34, 0.35);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hNpdfEtaCMForwardBackwardRatio[iEta], "EPPS21xCT18", "p");
        leg->AddEntry(hPdfEtaCMForwardBackwardRatio[iEta], "CT18", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_npdf2pdf_forwardBackwardRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );
        hNpdfEtaCMForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)

    // Plot full eta comparison of nPDF to PDF for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hNpdfEtaCMFull[iEta]->SetLineColor( kRed );
        hNpdfEtaCMFull[iEta]->SetMarkerColor( kRed );
        hPdfEtaCMFull[iEta]->SetLineColor( kBlue );
        hPdfEtaCMFull[iEta]->SetMarkerColor( kBlue );
        hNpdfEtaCMFull[iEta]->Draw();
        hPdfEtaCMFull[iEta]->Draw("same");
        hNpdfEtaCMFull[iEta]->GetXaxis()->SetRangeUser(-2.6, 2.6);
        hNpdfEtaCMFull[iEta]->GetYaxis()->SetRangeUser(0., 0.16);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "nPDF to PDF comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.6, 0.34, 0.7);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hNpdfEtaCMFull[iEta], "EPPS21xCT18", "p");
        leg->AddEntry(hPdfEtaCMFull[iEta], "CT18", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_npdf2pdf_fullEtaRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );
        hNpdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hNpdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfEtaCMFull[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfEtaCMFull[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)

    // Plot ratio of nPDF to PDF for different eta cuts on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hNpdfToPdfEtaCMFullRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hNpdfToPdfEtaCMFullRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hNpdfToPdfEtaCMFullRatio[iEta]->GetYaxis()->SetTitle("nPDF / PDF");
            hNpdfToPdfEtaCMFullRatio[iEta]->GetXaxis()->SetRangeUser(-2.6, 2.6);
            hNpdfToPdfEtaCMFullRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.15);
            leg = new TLegend(0.35, 0.2, 0.55, 0.55);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hNpdfToPdfEtaCMFullRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "nPDF/PDF");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_npdf2pdf_fullEtaRatio_etaCM_ratios.pdf", 
                    date.Data(), collSystemStr.Data()) );

    // Plot double ratio of nPDF to PDF forward/backward ratio for different eta cuts on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward/Backward (nPDF/PDF)");
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.6);
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.05);
            leg = new TLegend(0.65, 0.6, 0.85, 0.88);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hNpdfToPdfEtaCMForwardBackwardRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "nPDF/PDF");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_npdf2pdf_forwardBackwardRatio_doubleRatio_etaCM_ratios.pdf", 
                    date.Data(), collSystemStr.Data()) );

    //
    // Perform uncertainty correction for forward/backward double ratio and plot comparisons
    // on the same canvas again to see the effect of the correction
    //
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        recalculateUncertaintyOfNpdf( hNpdfToPdfEtaCMForwardBackwardRatio[iEta] );
        if (iEta == 0) {
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward/Backward (nPDF/PDF)");
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.6);
            hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.05);
            leg = new TLegend(0.65, 0.6, 0.85, 0.88);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hNpdfToPdfEtaCMForwardBackwardRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "nPDF/PDF");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_npdf2pdf_forwardBackwardRatio_doubleRatio_etaCM_ratios_uncertaintyCorrected.pdf", 
                    date.Data(), collSystemStr.Data()) );

    // Look at the ratio to the most inclusive cut and plot them on the same canvas
    for (int iEta{1}; iEta<etaCMCutSize; iEta++) {
        hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1] = (TH1D*)hNpdfToPdfEtaCMForwardBackwardRatio[iEta]->Clone(Form("hNpdfToPdfEtaCMForwardBackwardRatioTo25_%d", etaCMIntCut[iEta]));
        hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->Divide(hNpdfToPdfEtaCMForwardBackwardRatio[0]);
        hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->SetDirectory(0);

        hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->Draw( (iEta==1) ? "" : "same" );
        if (iEta == 1) {
            hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->GetYaxis()->SetTitle("Ratios to |#eta_{CM}|<2.5");
            hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->GetXaxis()->SetRangeUser(0., 2.0);
            hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1]->GetYaxis()->SetRangeUser(0.96, 1.01);
            leg = new TLegend(0.2, 0.2, 0.45, 0.45);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hNpdfToPdfEtaCMForwardBackwardRatioTo25[iEta - 1], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p"); 
    } // for (int iEta{1}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "F/B ratios to the most inclusive one (|#eta_{CM}|<2.5)");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_npdf2pdf_forwardBackwardRatio_doubleRatio_ratios_to25.pdf", 
                    date.Data(), collSystemStr.Data()) );
}

//________________
void plotPdfAndPythiaToDataComparison(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TFile *fPdf = TFile::Open( Form("./npdf/ct18_pp8160_dijet_eta.root" ) );
    if ( !fPdf || fPdf->IsZombie() ) {
        std::cerr << "Error: Could not open npdf file." << std::endl;
        return;
    }

    TFile *fData = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4_jetId_eta19_newResidual.root") );
    // TFile *fData = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/MB_pPb8160_ak4_jetId_eta19_newResidual.root") );
    if ( !fData || fData->IsZombie() ) {
        std::cerr << "Error: Could not open data file." << std::endl;
        return;
    }

    TFile *fPythia = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_jetId_eta19.root") );
    // TFile *fPythia = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_jetId_eta19_ptHat_gt_30.root") );
    if ( !fPythia || fPythia->IsZombie() ) {
        std::cerr << "Error: Could not open pythia file." << std::endl;
        return;
    }

    bool isRecalculateStatUnc = true;
    static constexpr std::array<int, 6> etaCMIntCut = {19, 18, 17, 16, 15, 14};
    static constexpr std::array<double, 6> etaCMDoubleCut = {1.9, 1.8, 1.7, 1.6, 1.5, 1.4};
    static constexpr std::array<const char*, 6> p8Colors {"kBlack", "kBlue", "kRed", "kMagenta", "kP8Orange", "kP8Green"};
    const int etaCMCutSize = etaCMIntCut.size();
    const double ptLow = 60.;
    const double ptHigh = 80.;

    // PDF histograms
    TH1D *hPdfForwardCMEta[etaCMCutSize];
    TH1D *hPdfBackwardCMEta[etaCMCutSize];
    TH1D *hPdfForwardBackwardRatio[etaCMCutSize];
    TH1D *hPdfFullCMEta[etaCMCutSize];

    // Data histograms
    TH2D *hData2DForwardCMEta[etaCMCutSize];
    TH2D *hData2DBackwardCMEta[etaCMCutSize];

    TH1D *hDataForwardCMEta[etaCMCutSize];
    TH1D *hDataBackwardCMEta[etaCMCutSize];
    TH1D *hDataForwardBackwardRatio[etaCMCutSize];
    TH1D *hDataFullCMEta[etaCMCutSize];

    TH1D *hDataToPdfForwardBackwardRatio[etaCMCutSize];
    TH1D *hDataToPdfFullEtaRatio[etaCMCutSize];

    // Pythia histograms
    TH2D *hPythia2DForwardCMEta[etaCMCutSize];
    TH2D *hPythia2DBackwardCMEta[etaCMCutSize];

    TH1D *hPythiaForwardCMEta[etaCMCutSize];
    TH1D *hPythiaBackwardCMEta[etaCMCutSize];
    TH1D *hPythiaForwardBackwardRatio[etaCMCutSize];
    TH1D *hPythiaFullCMEta[etaCMCutSize];

    TH1D *hDataToPythiaForwardBackwardRatio[etaCMCutSize];
    TH1D *hDataToPythiaFullEtaRatio[etaCMCutSize];

    // Pythia to PDF ratios
    TH1D *hPDFToPythiaFullEtaRatio[etaCMCutSize];

    // Read histograms from files
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        //
        // PDF
        //

        // Forward
        hPdfForwardCMEta[iEta] = dynamic_cast<TH1D *>( fPdf->Get( Form("hDijetEtaCMForward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfForwardCMEta[iEta] ) {
            std::cerr << "Error: hDijetEtaCMForward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfForwardCMEta[iEta]->SetDirectory(0);
        set1DStyle( hPdfForwardCMEta[iEta], 1);
        hPdfForwardCMEta[iEta]->SetName( Form("hPdfForwardCMEta_%d", iEta) );
        hPdfForwardCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfForwardCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Backward
        hPdfBackwardCMEta[iEta] = dynamic_cast<TH1D *>( fPdf->Get( Form("hDijetEtaCMBackward_%d_0", etaCMIntCut[iEta]) ) );
        if ( !hPdfBackwardCMEta[iEta] ) {
            std::cerr << "Error: hDijetEtaCMBackward_" << etaCMIntCut[iEta] << " not found in npdf file." << std::endl;
            return;
        }
        hPdfBackwardCMEta[iEta]->SetDirectory(0);
        set1DStyle( hPdfBackwardCMEta[iEta], 1);
        hPdfBackwardCMEta[iEta]->SetName( Form("hPdfBackwardCMEta_%d", iEta) );
        hPdfBackwardCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfBackwardCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Optionally recalculate statistical uncertainties 
        if (isRecalculateStatUnc) {
            recalculateStatisticalUncertainty( hPdfForwardCMEta[iEta] );
            recalculateStatisticalUncertainty( hPdfBackwardCMEta[iEta] );
        }

        // Make full eta distribution by combining forward and backward
        hPdfFullCMEta[iEta] = new TH1D( Form("hPdfFullCMEta_%d", iEta), 
                                         Form("hPdfFullCMEta_%d;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", iEta), 
                                         20, -2.0, 2.0 );
        hPdfFullCMEta[iEta]->Sumw2();
        hPdfFullCMEta[iEta]->SetDirectory(0);
        makeFullEtaFromForwardAndBackward( hPdfFullCMEta[iEta], hPdfForwardCMEta[iEta], hPdfBackwardCMEta[iEta] );
        set1DStyle( hPdfFullCMEta[iEta], 1, true );
        hPdfFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Ratio
        hPdfForwardBackwardRatio[iEta] = (TH1D*)hPdfForwardCMEta[iEta]->Clone(Form("hPdfForwardBackwardRatio_%d_0", etaCMIntCut[iEta]));
        hPdfForwardBackwardRatio[iEta]->Divide(hPdfBackwardCMEta[iEta]);
        hPdfForwardBackwardRatio[iEta]->SetDirectory(0);

        //
        // Data
        //

        // Forward
        hData2DForwardCMEta[iEta] = dynamic_cast<TH2D *>( fData->Get( Form("hRecoDijetPtEtaForwardArr_%d", etaCMCutSize - iEta - 1) ) );
        if ( !hData2DForwardCMEta[iEta] ) {
            std::cerr << "Error: hRecoDijetPtEtaForwardArr_" << etaCMCutSize - iEta - 1 << " not found in data file." << std::endl;
            return;
        }
        hData2DForwardCMEta[iEta]->SetDirectory(0);
        hData2DForwardCMEta[iEta]->RebinY(2);
        int ptBinLow = hData2DForwardCMEta[iEta]->GetXaxis()->FindBin(ptLow);
        int ptBinHigh = hData2DForwardCMEta[iEta]->GetXaxis()->FindBin(ptHigh) - 1;
        hDataForwardCMEta[iEta] = dynamic_cast<TH1D *>( hData2DForwardCMEta[iEta]->ProjectionY( Form("hDataForwardCMEta_%d", iEta), ptBinLow, ptBinHigh ) );
        set1DStyle( hDataForwardCMEta[iEta], 0 );
        hDataForwardCMEta[iEta]->SetName( Form("hDataForwardCMEta_%d", iEta) );
        hDataForwardCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hDataForwardCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Backward
        hData2DBackwardCMEta[iEta] = dynamic_cast<TH2D *>( fData->Get( Form("hRecoDijetPtEtaBackwardArr_%d", etaCMCutSize - iEta - 1) ) );
        if ( !hData2DBackwardCMEta[iEta] ) {
            std::cerr << "Error: hRecoDijetPtEtaBackwardArr_" << etaCMCutSize - iEta - 1 << " not found in data file." << std::endl;
            return;
        }
        hData2DBackwardCMEta[iEta]->SetDirectory(0);
        hData2DBackwardCMEta[iEta]->RebinY(2);
        hDataBackwardCMEta[iEta] = dynamic_cast<TH1D *>( hData2DBackwardCMEta[iEta]->ProjectionY( Form("hDataBackwardCMEta_%d", iEta), ptBinLow, ptBinHigh ) );
        set1DStyle( hDataBackwardCMEta[iEta], 0 );
        hDataBackwardCMEta[iEta]->SetName( Form("hDataBackwardCMEta_%d", iEta) );
        hDataBackwardCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hDataBackwardCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Make full eta distribution by combining forward and backward
        hDataFullCMEta[iEta] = new TH1D( Form("hDataFullCMEta_%d", iEta), 
                                         Form("hDataFullCMEta_%d;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", iEta), 
                                         20, -2.0, 2.0 );
        hDataFullCMEta[iEta]->Sumw2();
        hDataFullCMEta[iEta]->SetDirectory(0);
        makeFullEtaFromForwardAndBackward( hDataFullCMEta[iEta], hDataForwardCMEta[iEta], hDataBackwardCMEta[iEta] );
        set1DStyle( hDataFullCMEta[iEta], 0, true );
        hDataFullCMEta[iEta]->SetName( Form("hDataFullCMEta_%d", iEta) );
        hDataFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hDataFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Ratio
        hDataForwardBackwardRatio[iEta] = (TH1D*)hDataForwardCMEta[iEta]->Clone(Form("hDataForwardBackwardRatio_%d", iEta));
        hDataForwardBackwardRatio[iEta]->Divide(hDataBackwardCMEta[iEta]);
        hDataForwardBackwardRatio[iEta]->SetDirectory(0);

        //
        // PYTHIA
        //

        // Forward
        hPythia2DForwardCMEta[iEta] = dynamic_cast<TH2D *>( fPythia->Get( Form("hGenDijetPtEtaForwardArr_%d", etaCMCutSize - iEta - 1) ) );
        if ( !hPythia2DForwardCMEta[iEta] ) {
            std::cerr << "Error: hGenDijetPtEtaForwardArr_" << etaCMCutSize - iEta - 1 << " not found in pythia file." << std::endl;
            return;
        }
        hPythia2DForwardCMEta[iEta]->SetDirectory(0);
        hPythia2DForwardCMEta[iEta]->RebinY(2);
        hPythiaForwardCMEta[iEta] = dynamic_cast<TH1D *>( hPythia2DForwardCMEta[iEta]->ProjectionY( Form("hPythiaForwardCMEta_%d", iEta), ptBinLow, ptBinHigh ) );
        set1DStyle( hPythiaForwardCMEta[iEta], 1 );
        hPythiaForwardCMEta[iEta]->SetName( Form("hPythiaForwardCMEta_%d", iEta) );
        hPythiaForwardCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPythiaForwardCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Backward
        hPythia2DBackwardCMEta[iEta] = dynamic_cast<TH2D *>( fPythia->Get( Form("hGenDijetPtEtaBackwardArr_%d", etaCMCutSize - iEta - 1) ) );
        if ( !hPythia2DBackwardCMEta[iEta] ) {
            std::cerr << "Error: hGenDijetPtEtaBackwardArr_" << etaCMCutSize - iEta - 1 << " not found in pythia file." << std::endl;
            return;
        }
        hPythia2DBackwardCMEta[iEta]->SetDirectory(0);
        hPythia2DBackwardCMEta[iEta]->RebinY(2);
        hPythiaBackwardCMEta[iEta] = dynamic_cast<TH1D *>( hPythia2DBackwardCMEta[iEta]->ProjectionY( Form("hPythiaBackwardCMEta_%d", iEta), ptBinLow, ptBinHigh ) );
        set1DStyle( hPythiaBackwardCMEta[iEta], 1 );
        hPythiaBackwardCMEta[iEta]->SetName( Form("hPythiaBackwardCMEta_%d", iEta) );
        hPythiaBackwardCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPythiaBackwardCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Make full eta distribution by combining forward and backward
        hPythiaFullCMEta[iEta] = new TH1D( Form("hPythiaFullCMEta_%d", iEta), 
                                         Form("hPythiaFullCMEta_%d;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", iEta), 
                                         20, -2.0, 2.0 );
        hPythiaFullCMEta[iEta]->Sumw2();
        hPythiaFullCMEta[iEta]->SetDirectory(0);
        makeFullEtaFromForwardAndBackward( hPythiaFullCMEta[iEta], hPythiaForwardCMEta[iEta], hPythiaBackwardCMEta[iEta] );
        set1DStyle( hPythiaFullCMEta[iEta], 1, true );
        hPythiaFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPythiaFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );

        // Ratio
        hPythiaForwardBackwardRatio[iEta] = (TH1D*)hPythiaForwardCMEta[iEta]->Clone(Form("hPythiaForwardBackwardRatio_%d", iEta));
        hPythiaForwardBackwardRatio[iEta]->Divide(hPythiaBackwardCMEta[iEta]);
        hPythiaForwardBackwardRatio[iEta]->SetDirectory(0);

        //
        // Data to PDF ratios
        //

        // Forward/backward ratio
        hDataToPdfForwardBackwardRatio[iEta] = (TH1D*)hDataForwardBackwardRatio[iEta]->Clone(Form("hDataToPdfForwardBackwardRatio_%d", iEta));
        hDataToPdfForwardBackwardRatio[iEta]->Divide(hPdfForwardBackwardRatio[iEta]);
        hDataToPdfForwardBackwardRatio[iEta]->SetDirectory(0);

        // Full eta distribution
        hDataToPdfFullEtaRatio[iEta] = (TH1D*)hDataFullCMEta[iEta]->Clone(Form("hDataToPdfFullEtaRatio_%d", iEta));
        hDataToPdfFullEtaRatio[iEta]->Divide(hPdfFullCMEta[iEta]);
        hDataToPdfFullEtaRatio[iEta]->SetDirectory(0);

        //
        // Data to PYTHIA ratios
        //

        // Forward/backward ratio
        hDataToPythiaForwardBackwardRatio[iEta] = (TH1D*)hDataForwardBackwardRatio[iEta]->Clone(Form("hDataToPythiaForwardBackwardRatio_%d", iEta));
        hDataToPythiaForwardBackwardRatio[iEta]->Divide(hPythiaForwardBackwardRatio[iEta]);
        hDataToPythiaForwardBackwardRatio[iEta]->SetDirectory(0);

        // Full eta distribution
        hDataToPythiaFullEtaRatio[iEta] = (TH1D*)hDataFullCMEta[iEta]->Clone(Form("hDataToPythiaFullEtaRatio_%d", iEta));
        hDataToPythiaFullEtaRatio[iEta]->Divide(hPythiaFullCMEta[iEta]);
        hDataToPythiaFullEtaRatio[iEta]->SetDirectory(0);

        // PDF to PYTHIA ratio for full eta distribution
        hPDFToPythiaFullEtaRatio[iEta] = (TH1D*)hPdfFullCMEta[iEta]->Clone(Form("hPDFToPythiaFullEtaRatio_%d", iEta));
        hPDFToPythiaFullEtaRatio[iEta]->Divide(hPythiaFullCMEta[iEta]);
        hPDFToPythiaFullEtaRatio[iEta]->SetDirectory(0);
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)

    // Close files
    fPdf->Close();
    fData->Close();
    fPythia->Close();

    // Plotting components
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);
    TLegend *leg{nullptr};

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();

    // Plot forward/backward ratios for data on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataForwardBackwardRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hDataForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
            hDataForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.3);
            hDataForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
            leg = new TLegend(0.16, 0.5, 0.34, 0.8);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hDataForwardBackwardRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int i{0}; i<etaCMCutSize; i++)
     plotCMSHeader(collisionSystem, collisionEnergy);
     t.DrawLatexNDC(0.4, 0.85, "MB data");
     t.DrawLatexNDC(0.4, 0.8, "60 < p_{T}^{ave} < 80 GeV");
     leg->Draw();
     c->SaveAs(Form("%s/MB_%s_forwardBackwardRatio_etaCM_comparison.pdf", 
                    date.Data(), collSystemStr.Data()) );

    // Plot full eta comparison of data to PDF for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {

        hDataFullCMEta[iEta]->SetLineColor( kRed );
        hDataFullCMEta[iEta]->SetMarkerColor( kRed );
        hPdfFullCMEta[iEta]->SetLineColor( kBlue );
        hPdfFullCMEta[iEta]->SetMarkerColor( kBlue );

        hDataFullCMEta[iEta]->Draw();
        hPdfFullCMEta[iEta]->Draw("same");
        hDataFullCMEta[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
        hDataFullCMEta[iEta]->GetYaxis()->SetRangeUser(0., 0.16);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "MB to PDF comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.6, 0.34, 0.7);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hDataFullCMEta[iEta], "MB data", "p");
        leg->AddEntry(hPdfFullCMEta[iEta], "PDF (CT18)", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_data2PDF_fullEtaRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );

        hDataFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hDataFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int i{0}; i<etaCMCutSize; i++)

    // Plot full eta comparison of data to PYTHIA for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataFullCMEta[iEta]->SetLineColor( kRed );
        hDataFullCMEta[iEta]->SetMarkerColor( kRed );
        hPythiaFullCMEta[iEta]->SetLineColor( kBlue );
        hPythiaFullCMEta[iEta]->SetMarkerColor( kBlue );

        hDataFullCMEta[iEta]->Draw();
        hPythiaFullCMEta[iEta]->Draw("same");
        hDataFullCMEta[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
        hDataFullCMEta[iEta]->GetYaxis()->SetRangeUser(0., 0.16);
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "MB to PYTHIA comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.6, 0.34, 0.7);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hDataFullCMEta[iEta], "MB data", "p");
        leg->AddEntry(hPythiaFullCMEta[iEta], "PYTHIA", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_data2PYTHIA_fullEtaRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );

        hDataFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hDataFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPythiaFullCMEta[iEta]->SetLineColor( p8Colors[iEta] );
        hPythiaFullCMEta[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int i{0}; i<etaCMCutSize; i++)

    // Plot forward/backward ratio comparison of data to PDF for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataForwardBackwardRatio[iEta]->SetLineColor( kRed );
        hDataForwardBackwardRatio[iEta]->SetMarkerColor( kRed );
        hPdfForwardBackwardRatio[iEta]->SetLineColor( kBlue );
        hPdfForwardBackwardRatio[iEta]->SetMarkerColor( kBlue );

        hDataForwardBackwardRatio[iEta]->Draw();
        hPdfForwardBackwardRatio[iEta]->Draw("same");
        hDataForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
        hDataForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.1);
        hDataForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "MB to PDF comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.25, 0.34, 0.35);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hDataForwardBackwardRatio[iEta], "MB data", "p");
        leg->AddEntry(hPdfForwardBackwardRatio[iEta], "PDF (CT18)", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_data2PDF_forwardBackwardRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );
        hDataForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hDataForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPdfForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hPdfForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int i{0}; i<etaCMCutSize; i++)

    // Plot forward/backward ratio comparison of data to PYTHIA for different eta cuts individually
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataForwardBackwardRatio[iEta]->SetLineColor( kRed );
        hDataForwardBackwardRatio[iEta]->SetMarkerColor( kRed );
        hPythiaForwardBackwardRatio[iEta]->SetLineColor( kBlue );
        hPythiaForwardBackwardRatio[iEta]->SetMarkerColor( kBlue );

        hDataForwardBackwardRatio[iEta]->Draw();
        hPythiaForwardBackwardRatio[iEta]->Draw("same");
        hDataForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
        hDataForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.1);
        hDataForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.SetTextSize(0.04);
        t.DrawLatexNDC(0.16, 0.85, "MB to PYTHIA comparison");
        t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
        t.DrawLatexNDC(0.16, 0.75, Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]) );
        leg = new TLegend(0.16, 0.25, 0.34, 0.35);
        leg->SetTextSize(0.03);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hDataForwardBackwardRatio[iEta], "MB data", "p");
        leg->AddEntry(hPythiaForwardBackwardRatio[iEta], "PYTHIA", "p");
        leg->Draw();
        c->SaveAs(Form("%s/%s_data2PYTHIA_forwardBackwardRatio_etaCM_comparison_%d.pdf", 
                    date.Data(), collSystemStr.Data(), etaCMIntCut[iEta]) );
        hDataForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hDataForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
        hPythiaForwardBackwardRatio[iEta]->SetLineColor( p8Colors[iEta] );
        hPythiaForwardBackwardRatio[iEta]->SetMarkerColor( p8Colors[iEta] );
    } // for (int i{0}; i<etaCMCutSize; i++)

    // Compilation of ratios of data to PDF for different eta cuts on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataToPdfFullEtaRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hDataToPdfFullEtaRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hDataToPdfFullEtaRatio[iEta]->GetYaxis()->SetTitle("Data / PDF");
            hDataToPdfFullEtaRatio[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
            hDataToPdfFullEtaRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.4);
            leg = new TLegend(0.35, 0.5, 0.65, 0.75);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hDataToPdfFullEtaRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "MB to PDF comparison");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_data2PDF_fullEtaRatio_etaCM_ratios.pdf", 
                    date.Data(), collSystemStr.Data()) );

    // Compilation of ratios of data to PYTHIA for different eta cuts on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hDataToPythiaFullEtaRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hDataToPythiaFullEtaRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hDataToPythiaFullEtaRatio[iEta]->GetYaxis()->SetTitle("Data / PYTHIA");
            hDataToPythiaFullEtaRatio[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
            hDataToPythiaFullEtaRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.4);
            leg = new TLegend(0.35, 0.5, 0.65, 0.75);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hDataToPythiaFullEtaRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "MB to PYTHIA comparison");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_data2PYTHIA_fullEtaRatio_etaCM_ratios.pdf", 
                    date.Data(), collSystemStr.Data()) );

    // Compilation of ratios of PDF to PYTHIA for different eta cuts on the same canvas
    for (int iEta{0}; iEta<etaCMCutSize; iEta++) {
        hPDFToPythiaFullEtaRatio[iEta]->SetMarkerStyle(kFullCircle);
        hPDFToPythiaFullEtaRatio[iEta]->Draw( (iEta==0) ? "" : "same" );
        if (iEta == 0) {
            hPDFToPythiaFullEtaRatio[iEta]->GetXaxis()->SetTitle("#eta^{dijet}_{CM}");
            hPDFToPythiaFullEtaRatio[iEta]->GetYaxis()->SetTitle("PDF / PYTHIA");
            hPDFToPythiaFullEtaRatio[iEta]->GetXaxis()->SetRangeUser(-2.0, 2.0);
            hPDFToPythiaFullEtaRatio[iEta]->GetYaxis()->SetRangeUser(0.8, 1.1);
            leg = new TLegend(0.45, 0.2, 0.65, 0.45);
            leg->SetTextSize(0.03);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
        }
        leg->AddEntry(hPDFToPythiaFullEtaRatio[iEta], Form("|#eta^{jet}_{CM}| < %.1f", etaCMDoubleCut[iEta]), "p");
    } // for (int iEta{0}; iEta<etaCMCutSize; iEta++)
    plotCMSHeader(collisionSystem, collisionEnergy);
    t.SetTextSize(0.04);
    t.DrawLatexNDC(0.16, 0.85, "PDF to PYTHIA comparison");
    t.DrawLatexNDC(0.16, 0.8, "60 < p_{T}^{ave} < 80 GeV");
    leg->Draw();
    c->SaveAs(Form("%s/%s_PDF2PYTHIA_fullEtaRatio_etaCM_ratios.pdf", 
                    date.Data(), collSystemStr.Data()) ); 
}

//________________
void setHistogramStyle(TH1D *h, int color, int lineWidth = 2, double markerSize = 1.3, int markerStyle = 20) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(markerSize);
    h->SetLineWidth(lineWidth);
}

//________________
void plotSingleJetClosures(int collisionSystem, double collisionEnergy, TString date) {

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // pT bin selection
    double deltaPt = 0.001;
    double ptLow = 80. + deltaPt;
    double ptHigh = 120. - deltaPt;
    int ptLowInt = std::round(ptLow);
    int ptHighInt = std::round(ptHigh);
    double etaLabMax = 2.4;
    double etaCMMax = 2.0;
    double yAxisRatio2Gen[2] = {0.85, 1.15};

    std::vector<double> ptVals = {40., 80., 120., 180., 250., 500.};
    int nPtBins = ptVals.size() - 1;


    TString mcType = "PYTHIA"; // "pythia" or "embedding"
    TString mcTypeLower = mcType;
    mcTypeLower.ToLower();

    TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/pgoing/oPythia_pgoing_def_ak4_jetId_eta19.root") );
    // TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_jetId_eta19.root") );
    // TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_jetId_eta19.root") );

    // TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/Pbgoing/oPythia_pPb8160_Pbgoing_30.root") );
    // TFile *fPythia = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_jetId_eta19_ptHat_gt_30.root") );
    if ( !fMC || fMC->IsZombie() ) {
        std::cerr << "Error: Could not open pythia file." << std::endl;
        return;
    }

    // Gen histograms
    TH2D *hGenInclusiveJetPtEta;
    TH2D *hGenInclusiveJetPtEtaStdBins;
    TH2D *hGenLeadJetPtEta;
    TH2D *hGenLeadJetPtEtaStdBins;
    TH2D *hGenSubLeadJetPtEta;
    TH2D *hGenSubLeadJetPtEtaStdBins;
    TH2D *hGenInclusiveJetPtEtaCM;
    TH2D *hGenLeadJetPtEtaCM;
    TH2D *hGenSubLeadJetPtEtaCM;

    // From selected dijets
    TH2D *hGenDijetLeadPtEtaCM;
    TH2D *hGenDijetLeadPtEtaStdBins;
    TH2D *hGenDijetSubLeadPtEtaCM;
    TH2D *hGenDijetSubLeadPtEtaStdBins;

    // Reco histograms
    TH2D *hRecoInclusiveJetPtEta;
    TH2D *hRecoInclusiveJetPtEtaStdBins;
    TH2D *hRecoLeadJetPtEta;
    TH2D *hRecoSubLeadJetPtEta;
    TH2D *hRecoInclusiveJetPtEtaCM;
    TH2D *hRecoLeadJetPtEtaCM;
    TH2D *hRecoSubLeadJetPtEtaCM;

    // From selected dijets
    TH2D *hRecoDijetLeadPtEtaCM;
    TH2D *hRecoDijetLeadPtEtaStdBins;
    TH2D *hRecoDijetSubLeadPtEtaCM;
    TH2D *hRecoDijetSubLeadPtEtaStdBins;

    // Ref histograms
    TH2D *hRefInclusiveJetPtEta;
    TH2D *hRefLeadJetPtEta;
    TH2D *hRefSubLeadJetPtEta;
    TH2D *hRefInclusiveJetPtEtaCM;
    TH2D *hRefLeadJetPtEtaCM;
    TH2D *hRefSubLeadJetPtEtaCM;

    // From selected dijets
    TH2D *hRefDijetLeadPtEtaCM;
    TH2D *hRefDijetLeadPtEtaStdBins;
    TH2D *hRefDijetSubLeadPtEtaCM;
    TH2D *hRefDijetSubLeadPtEtaStdBins;

    //////////////////////////////////////////
    //    Retrieve histograms from MC file  //
    //////////////////////////////////////////

    // Gen level
    
    // Lab frame
    hGenInclusiveJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hGenInclusiveJetPtEta") );
    if ( !hGenInclusiveJetPtEta ) {
        std::cerr << "Error: hGenInclusiveJetPtEta not found in MC file." << std::endl;
        return;
    }
    hGenInclusiveJetPtEta->SetDirectory(0);
    hGenInclusiveJetPtEta->SetName("hGenInclusiveJetPtEta");
    hGenInclusiveJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");
    hGenInclusiveJetPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hGenInclusiveJetPtEtaStdBins") );
    if ( !hGenInclusiveJetPtEtaStdBins ) {
        std::cerr << "Error: hGenInclusiveJetPtEtaStdBins not found in MC file." << std::endl;
        return;
    }
    hGenInclusiveJetPtEtaStdBins->SetDirectory(0);
    hGenInclusiveJetPtEtaStdBins->SetName("hGenInclusiveJetPtEtaStdBins");
    hGenInclusiveJetPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{std bins}");
    hGenLeadJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hGenLeadJetPtEta") );
    if ( !hGenLeadJetPtEta ) {
        std::cerr << "Error: hGenLeadJetPtEta not found in MC file." << std::endl;
        return;
    }
    hGenLeadJetPtEta->SetDirectory(0);
    hGenLeadJetPtEta->SetName("hGenLeadJetPtEta");
    hGenLeadJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");
    hGenSubLeadJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hGenSubLeadJetPtEta") );
    if ( !hGenSubLeadJetPtEta ) {
        std::cerr << "Error: hGenSubLeadJetPtEta not found in MC file." << std::endl;
        return;
    }
    hGenSubLeadJetPtEta->SetDirectory(0);
    hGenSubLeadJetPtEta->SetName("hGenSubLeadJetPtEta");
    hGenSubLeadJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");

    // CM frame
    hGenInclusiveJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hGenInclusiveJetPtEtaCM") );
    if ( !hGenInclusiveJetPtEtaCM ) {
        std::cerr << "Error: hGenInclusiveJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hGenInclusiveJetPtEtaCM->SetDirectory(0);
    hGenInclusiveJetPtEtaCM->SetName("hGenInclusiveJetPtEtaCM");
    hGenInclusiveJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hGenLeadJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hGenLeadJetPtEtaCM") );
    if ( !hGenLeadJetPtEtaCM ) {
        std::cerr << "Error: hGenLeadJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hGenLeadJetPtEtaCM->SetDirectory(0);
    hGenLeadJetPtEtaCM->SetName("hGenLeadJetPtEtaCM");
    hGenLeadJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hGenSubLeadJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hGenSubLeadJetPtEtaCM") );
    if ( !hGenSubLeadJetPtEtaCM ) {
        std::cerr << "Error: hGenSubLeadJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hGenSubLeadJetPtEtaCM->SetDirectory(0);
    hGenSubLeadJetPtEtaCM->SetName("hGenSubLeadJetPtEtaCM");
    hGenSubLeadJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");

    // From selected dijets
    hGenDijetLeadPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hGenDijetLeadPtEtaCM") );
    if ( !hGenDijetLeadPtEtaCM ) {
        std::cerr << "Error: hGenDijetLeadPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hGenDijetLeadPtEtaCM->SetDirectory(0);
    hGenDijetLeadPtEtaCM->SetName("hGenDijetLeadPtEtaCM");
    hGenDijetLeadPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hGenDijetLeadPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hGenDijetLeadPtEtaLabStdBins") );
    if ( !hGenDijetLeadPtEtaStdBins ) {
        std::cerr << "Error: hGenDijetLeadPtEtaLabStdBins not found in MC file." << std::endl;
        return;
    }
    hGenDijetLeadPtEtaStdBins->SetDirectory(0);
    hGenDijetLeadPtEtaStdBins->SetName("hGenDijetLeadPtEtaStdBins");
    hGenDijetLeadPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{CM, std bins}");
    hGenDijetSubLeadPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hGenDijetSubLeadPtEtaCM") );
    if ( !hGenDijetSubLeadPtEtaCM ) {
        std::cerr << "Error: hGenDijetSubLeadPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hGenDijetSubLeadPtEtaCM->SetDirectory(0);
    hGenDijetSubLeadPtEtaCM->SetName("hGenDijetSubLeadPtEtaCM");
    hGenDijetSubLeadPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hGenDijetSubLeadPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hGenDijetSubLeadPtEtaLabStdBins") );
    if ( !hGenDijetSubLeadPtEtaStdBins ) {
        std::cerr << "Error: hGenDijetSubLeadPtEtaLabStdBins not found in MC file." << std::endl;
        return;
    }
    hGenDijetSubLeadPtEtaStdBins->SetDirectory(0);
    hGenDijetSubLeadPtEtaStdBins->SetName("hGenDijetSubLeadPtEtaStdBins");
    hGenDijetSubLeadPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{CM, std bins}");


    // Reco level

    // Lab frame
    hRecoInclusiveJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hRecoInclusiveAllJetPtEta") );
    if ( !hRecoInclusiveJetPtEta ) {
        std::cerr << "Error: hRecoInclusiveAllJetPtEta not found in MC file." << std::endl;
        return;
    }
    hRecoInclusiveJetPtEta->SetDirectory(0);
    hRecoInclusiveJetPtEta->SetName("hRecoInclusiveJetPtEta");
    hRecoInclusiveJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");
    hRecoInclusiveJetPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hRecoInclusiveAllJetPtEtaStdBins") );
    if ( !hRecoInclusiveJetPtEtaStdBins ) {
        std::cerr << "Error: hRecoInclusiveAllJetPtEtaStdBins not found in MC file." << std::endl;
        return;
    }
    hRecoInclusiveJetPtEtaStdBins->SetDirectory(0);
    hRecoInclusiveJetPtEtaStdBins->SetName("hRecoInclusiveJetPtEtaStdBins");
    hRecoInclusiveJetPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{std bins}");
    hRecoLeadJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hRecoLeadAllJetPtEta") );
    if ( !hRecoLeadJetPtEta ) {
        std::cerr << "Error: hRecoLeadAllJetPtEta not found in MC file." << std::endl;
        return;
    }
    hRecoLeadJetPtEta->SetDirectory(0);
    hRecoLeadJetPtEta->SetName("hRecoLeadJetPtEta");
    hRecoLeadJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");
    hRecoSubLeadJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hRecoSubLeadAllJetPtEta") );
    if ( !hRecoSubLeadJetPtEta ) {
        std::cerr << "Error: hRecoSubLeadAllJetPtEta not found in MC file." << std::endl;
        return;
    }
    hRecoSubLeadJetPtEta->SetDirectory(0);
    hRecoSubLeadJetPtEta->SetName("hRecoSubLeadJetPtEta");
    hRecoSubLeadJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");

    // CM frame
    hRecoInclusiveJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRecoInclusiveAllJetPtEtaCM") );
    if ( !hRecoInclusiveJetPtEtaCM ) {
        std::cerr << "Error: hRecoInclusiveAllJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRecoInclusiveJetPtEtaCM->SetDirectory(0);
    hRecoInclusiveJetPtEtaCM->SetName("hRecoInclusiveJetPtEtaCM");
    hRecoInclusiveJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRecoLeadJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRecoLeadAllJetPtEtaCM") );
    if ( !hRecoLeadJetPtEtaCM ) {
        std::cerr << "Error: hRecoLeadAllJetPtEtaCM not found in MC file." << std::endl;
        return; 
    }
    hRecoLeadJetPtEtaCM->SetDirectory(0);
    hRecoLeadJetPtEtaCM->SetName("hRecoLeadJetPtEtaCM");
    hRecoLeadJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRecoSubLeadJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRecoSubLeadAllJetPtEtaCM") );
    if ( !hRecoSubLeadJetPtEtaCM ) {
        std::cerr << "Error: hRecoSubLeadAllJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRecoSubLeadJetPtEtaCM->SetDirectory(0);
    hRecoSubLeadJetPtEtaCM->SetName("hRecoSubLeadJetPtEtaCM");
    hRecoSubLeadJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");

    // From selected dijets
    hRecoDijetLeadPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRecoDijetLeadPtEtaCM") );
    if ( !hRecoDijetLeadPtEtaCM ) {
        std::cerr << "Error: hRecoDijetLeadPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRecoDijetLeadPtEtaCM->SetDirectory(0);
    hRecoDijetLeadPtEtaCM->SetName("hRecoDijetLeadPtEtaCM");
    hRecoDijetLeadPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRecoDijetLeadPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hRecoDijetLeadPtEtaLabStdBins") );
    if ( !hRecoDijetLeadPtEtaStdBins ) {
        std::cerr << "Error: hRecoDijetLeadPtEtaLabStdBins not found in MC file." << std::endl;
        return;
    }
    hRecoDijetLeadPtEtaStdBins->SetDirectory(0);
    hRecoDijetLeadPtEtaStdBins->SetName("hRecoDijetLeadPtEtaStdBins");
    hRecoDijetLeadPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{CM, std bins}");
    hRecoDijetSubLeadPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRecoDijetSubLeadPtEtaCM") );
    if ( !hRecoDijetSubLeadPtEtaCM ) {
        std::cerr << "Error: hRecoDijetSubLeadPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRecoDijetSubLeadPtEtaCM->SetDirectory(0);
    hRecoDijetSubLeadPtEtaCM->SetName("hRecoDijetSubLeadPtEtaCM");
    hRecoDijetSubLeadPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRecoDijetSubLeadPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hRecoDijetSubLeadPtEtaLabStdBins") );
    if ( !hRecoDijetSubLeadPtEtaStdBins ) {
        std::cerr << "Error: hRecoDijetSubLeadPtEtaLabStdBins not found in MC file." << std::endl;
        return;
    }
    hRecoDijetSubLeadPtEtaStdBins->SetDirectory(0);
    hRecoDijetSubLeadPtEtaStdBins->SetName("hRecoDijetSubLeadPtEtaStdBins");
    hRecoDijetSubLeadPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{CM, std bins}");

    // Ref level

    // Lab frame
    hRefInclusiveJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hRefInclusiveJetPtEta") );
    if ( !hRefInclusiveJetPtEta ) {
        std::cerr << "Error: hRefInclusiveJetPtEta not found in MC file." << std::endl;
        return;
    }
    hRefInclusiveJetPtEta->SetDirectory(0);
    hRefInclusiveJetPtEta->SetName("hRefInclusiveJetPtEta");
    hRefInclusiveJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");
    hRefLeadJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hRefLeadJetPtEta") );
    if ( !hRefLeadJetPtEta ) {
        std::cerr << "Error: hRefLeadJetPtEta not found in MC file." << std::endl;
        return;
    }
    hRefLeadJetPtEta->SetDirectory(0);
    hRefLeadJetPtEta->SetName("hRefLeadJetPtEta");
    hRefLeadJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");
    hRefSubLeadJetPtEta = dynamic_cast<TH2D *>( fMC->Get("hRefSubLeadJetPtEta") );
    if ( !hRefSubLeadJetPtEta ) {
        std::cerr << "Error: hRefSubLeadJetPtEta not found in MC file." << std::endl;
        return;
    }
    hRefSubLeadJetPtEta->SetDirectory(0);
    hRefSubLeadJetPtEta->SetName("hRefSubLeadJetPtEta");
    hRefSubLeadJetPtEta->GetXaxis()->SetTitle("#eta^{jet}");

    // CM frame
    hRefInclusiveJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRefInclusiveJetPtEtaCM") );
    if ( !hRefInclusiveJetPtEtaCM ) {
        std::cerr << "Error: hRefInclusiveJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRefInclusiveJetPtEtaCM->SetDirectory(0);
    hRefInclusiveJetPtEtaCM->SetName("hRefInclusiveJetPtEtaCM");
    hRefInclusiveJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRefLeadJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRefLeadJetPtEtaCM") );
    if ( !hRefLeadJetPtEtaCM ) {
        std::cerr << "Error: hRefLeadJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRefLeadJetPtEtaCM->SetDirectory(0);
    hRefLeadJetPtEtaCM->SetName("hRefLeadJetPtEtaCM");
    hRefLeadJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRefSubLeadJetPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRefSubLeadJetPtEtaCM") );
    if ( !hRefSubLeadJetPtEtaCM ) {
        std::cerr << "Error: hRefSubLeadJetPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRefSubLeadJetPtEtaCM->SetDirectory(0);
    hRefSubLeadJetPtEtaCM->SetName("hRefSubLeadJetPtEtaCM");
    hRefSubLeadJetPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");

    // From selected dijets
    hRefDijetLeadPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRefDijetLeadPtEtaCM") );
    if ( !hRefDijetLeadPtEtaCM ) {
        std::cerr << "Error: hRefDijetLeadPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRefDijetLeadPtEtaCM->SetDirectory(0);
    hRefDijetLeadPtEtaCM->SetName("hRefDijetLeadPtEtaCM");
    hRefDijetLeadPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRefDijetLeadPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hRefDijetLeadPtEtaLabStdBins") );
    if ( !hRefDijetLeadPtEtaStdBins ) {
        std::cerr << "Error: hRefDijetLeadPtEtaLabStdBins not found in MC file." << std::endl;
        return;
    }
    hRefDijetLeadPtEtaStdBins->SetDirectory(0);
    hRefDijetLeadPtEtaStdBins->SetName("hRefDijetLeadPtEtaStdBins");
    hRefDijetLeadPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{CM, std bins}");   
    hRefDijetSubLeadPtEtaCM = dynamic_cast<TH2D *>( fMC->Get("hRefDijetSubLeadPtEtaCM") );
    if ( !hRefDijetSubLeadPtEtaCM ) {
        std::cerr << "Error: hRefDijetSubLeadPtEtaCM not found in MC file." << std::endl;
        return;
    }
    hRefDijetSubLeadPtEtaCM->SetDirectory(0);
    hRefDijetSubLeadPtEtaCM->SetName("hRefDijetSubLeadPtEtaCM");
    hRefDijetSubLeadPtEtaCM->GetXaxis()->SetTitle("#eta^{jet}_{CM}");
    hRefDijetSubLeadPtEtaStdBins = dynamic_cast<TH2D *>( fMC->Get("hRefDijetSubLeadPtEtaLabStdBins") );
    if ( !hRefDijetSubLeadPtEtaStdBins ) {
        std::cerr << "Error: hRefDijetSubLeadPtEtaLabStdBins not found in MC file." << std::endl;
        return;
    }
    hRefDijetSubLeadPtEtaStdBins->SetDirectory(0);
    hRefDijetSubLeadPtEtaStdBins->SetName("hRefDijetSubLeadPtEtaStdBins");
    hRefDijetSubLeadPtEtaStdBins->GetXaxis()->SetTitle("#eta^{jet}_{CM, std bins}");

    fMC->Close();

    /////////////////////////////////////////////
    //  Histograms for projections and ratios  //
    /////////////////////////////////////////////

    // Gen level

    // Lab frame
    TH1D *hGenInclusiveJetEta;
    TH1D *hGenInclusiveJetEtaStdBins;
    TH1D *hGenInclusiveJetEtaStdBinsPtBins[nPtBins];
    TH1D *hGenLeadJetEta;
    TH1D *hGenSubLeadJetEta;

    // CM frame
    TH1D *hGenInclusiveJetEtaCM;
    TH1D *hGenInclusiveJetEtaCMForward;
    TH1D *hGenInclusiveJetEtaCMBackward;
    TH1D *hGenInclusiveJetForwardBackwardRatioCM;
    TH1D *hGenInclusiveJetEtaCMPtBins[nPtBins];

    TH1D *hGenLeadJetEtaCM;
    TH1D *hGenLeadJetEtaCMForward;
    TH1D *hGenLeadJetEtaCMBackward;
    TH1D *hGenLeadJetForwardBackwardRatioCM;

    TH1D *hGenSubLeadJetEtaCM;
    TH1D *hGenSubLeadJetEtaCMForward;
    TH1D *hGenSubLeadJetEtaCMBackward;
    TH1D *hGenSubLeadJetForwardBackwardRatioCM;

    // From selected dijets
    TH1D *hGenDijetLeadJetEtaStdBins;
    TH1D *hGenDijetLeadJetEtaCM;
    TH1D *hGenDijetLeadJetEtaCMForward;
    TH1D *hGenDijetLeadJetEtaCMBackward;
    TH1D *hGenDijetLeadJetForwardBackwardRatioCM;

    TH1D *hGenDijetSubLeadJetEtaStdBins;
    TH1D *hGenDijetSubLeadJetEtaCM;
    TH1D *hGenDijetSubLeadJetEtaCMForward;
    TH1D *hGenDijetSubLeadJetEtaCMBackward;
    TH1D *hGenDijetSubLeadJetForwardBackwardRatioCM;

    // Reco level

    // Lab frame
    TH1D *hRecoInclusiveJetEta;
    TH1D *hRecoInclusiveJetEtaStdBins;
    TH1D *hRecoInclusiveJetEtaStdBinsPtBins[nPtBins];
    TH1D *hRecoLeadJetEta;
    TH1D *hRecoSubLeadJetEta;

    TH1D *hInclusiveJetReco2GenEta;
    TH1D *hInclusiveJetReco2GenEtaStdBins;
    TH1D *hInclusiveJetReco2GenEtaStdBinsPtBins[nPtBins];
    TH1D *hLeadJetReco2GenEta;
    TH1D *hSubLeadJetReco2GenEta;

    // CM frame
    TH1D *hRecoInclusiveJetEtaCM;
    TH1D *hRecoInclusiveJetEtaCMForward;
    TH1D *hRecoInclusiveJetEtaCMBackward;
    TH1D *hRecoInclusiveJetForwardBackwardRatioCM;
    TH1D *hRecoInclusiveJetEtaCMPtBins[nPtBins];

    TH1D *hRecoLeadJetEtaCM;
    TH1D *hRecoLeadJetEtaCMForward;
    TH1D *hRecoLeadJetEtaCMBackward;
    TH1D *hRecoLeadJetForwardBackwardRatioCM;

    TH1D *hRecoSubLeadJetEtaCM;
    TH1D *hRecoSubLeadJetEtaCMForward;
    TH1D *hRecoSubLeadJetEtaCMBackward;
    TH1D *hRecoSubLeadJetForwardBackwardRatioCM;

    // From selected dijets
    TH1D *hRecoDijetLeadJetEtaStdBins;
    TH1D *hRecoDijetLeadJetEtaCM;
    TH1D *hRecoDijetLeadJetEtaCMForward;
    TH1D *hRecoDijetLeadJetEtaCMBackward;
    TH1D *hRecoDijetLeadJetForwardBackwardRatioCM;

    TH1D *hRecoDijetSubLeadJetEtaStdBins;
    TH1D *hRecoDijetSubLeadJetEtaCM;
    TH1D *hRecoDijetSubLeadJetEtaCMForward;
    TH1D *hRecoDijetSubLeadJetEtaCMBackward;
    TH1D *hRecoDijetSubLeadJetForwardBackwardRatioCM;

    // Ratios to gen level
    TH1D *hInclusiveJetReco2GenEtaCM;
    TH1D *hLeadJetReco2GenEtaCM;
    TH1D *hSubLeadJetReco2GenEtaCM;
    TH1D *hDijetLeadJetReco2GenEtaCM;
    TH1D *hDijetSubLeadJetReco2GenEtaCM;
    TH1D *hDijetLeadJetReco2GenEtaStdBins;
    TH1D *hDijetSubLeadJetReco2GenEtaStdBins;

    // Ref level

    // Lab frame
    TH1D *hRefInclusiveJetEta;
    TH1D *hRefLeadJetEta;
    TH1D *hRefSubLeadJetEta;

    TH1D *hInclusiveJetRef2GenEta;
    TH1D *hLeadJetRef2GenEta;
    TH1D *hSubLeadJetRef2GenEta;

    // CM frame
    TH1D *hRefInclusiveJetEtaCM;
    TH1D *hRefInclusiveJetEtaCMForward;
    TH1D *hRefInclusiveJetEtaCMBackward;
    TH1D *hRefInclusiveJetForwardBackwardRatioCM;

    TH1D *hRefLeadJetEtaCM;
    TH1D *hRefLeadJetEtaCMForward;
    TH1D *hRefLeadJetEtaCMBackward;
    TH1D *hRefLeadJetForwardBackwardRatioCM;

    TH1D *hRefSubLeadJetEtaCM;
    TH1D *hRefSubLeadJetEtaCMForward;
    TH1D *hRefSubLeadJetEtaCMBackward;
    TH1D *hRefSubLeadJetForwardBackwardRatioCM;

    // From selected dijets
    TH1D *hRefDijetLeadJetEtaStdBins;
    TH1D *hRefDijetLeadJetEtaCM;
    TH1D *hRefDijetLeadJetEtaCMForward;
    TH1D *hRefDijetLeadJetEtaCMBackward;
    TH1D *hRefDijetLeadJetForwardBackwardRatioCM;

    TH1D *hRefDijetSubLeadJetEtaStdBins;
    TH1D *hRefDijetSubLeadJetEtaCM;
    TH1D *hRefDijetSubLeadJetEtaCMForward;
    TH1D *hRefDijetSubLeadJetEtaCMBackward;
    TH1D *hRefDijetSubLeadJetForwardBackwardRatioCM;

    // Ratios to gen level
    TH1D *hInclusiveJetRef2GenEtaCM;
    TH1D *hLeadJetRef2GenEtaCM;
    TH1D *hSubLeadJetRef2GenEtaCM;
    TH1D *hDijetLeadJetRef2GenEtaCM;
    TH1D *hDijetSubLeadJetRef2GenEtaCM;
    TH1D *hDijetLeadJetRef2GenEtaStdBins;
    TH1D *hDijetSubLeadJetRef2GenEtaStdBins;


    ///////////////////////////
    //    Make projections   //
    ///////////////////////////

    // Gen level

    // Lab frame
    hGenInclusiveJetEta = dynamic_cast<TH1D *>( hGenInclusiveJetPtEta->ProjectionX("hGenInclusiveJetEta", 
                                        hGenInclusiveJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hGenInclusiveJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    hGenInclusiveJetEtaStdBins = dynamic_cast<TH1D *>( hGenInclusiveJetPtEtaStdBins->ProjectionX("hGenInclusiveJetEtaStdBins", 
                                        hGenDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                        hGenDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    hGenLeadJetEta = dynamic_cast<TH1D *>( hGenLeadJetPtEta->ProjectionX("hGenLeadJetEta", 
                                        hGenLeadJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hGenLeadJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    hGenSubLeadJetEta = dynamic_cast<TH1D *>( hGenSubLeadJetPtEta->ProjectionX("hGenSubLeadJetEta", 
                                        hGenSubLeadJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hGenSubLeadJetPtEta->GetYaxis()->FindBin(ptHigh)) );

    set1DStyle(hGenInclusiveJetEta, 1, true);
    set1DStyle(hGenInclusiveJetEtaStdBins, 1);
    set1DStyle(hGenLeadJetEta, 1, true);
    set1DStyle(hGenSubLeadJetEta, 1, true);
    hGenInclusiveJetEta->Scale(1./hGenInclusiveJetEta->Integral( hGenInclusiveJetEta->GetXaxis()->FindBin(-1.0), 
                                                                 hGenInclusiveJetEta->GetXaxis()->FindBin(1.0) ) );
    // hGenInclusiveJetEtaStdBins->Scale(1./hGenInclusiveJetEtaStdBins->Integral("width"));
    rescaleHisto1D(hGenInclusiveJetEtaStdBins);
    hGenInclusiveJetEtaStdBins->Scale(1./hGenInclusiveJetEtaStdBins->Integral( hGenInclusiveJetEtaStdBins->GetXaxis()->FindBin(-1.0), 
                                                                               hGenInclusiveJetEtaStdBins->GetXaxis()->FindBin(1.0) ) );
    for (int iPtBin = 0; iPtBin < nPtBins; ++iPtBin) {
        hGenInclusiveJetEtaStdBinsPtBins[iPtBin] = dynamic_cast<TH1D *>( hGenInclusiveJetPtEtaStdBins->ProjectionX(Form("hGenInclusiveJetEtaStdBinsPtBin%d", iPtBin), 
                                                        hGenInclusiveJetPtEtaStdBins->GetYaxis()->FindBin(ptVals[iPtBin]+deltaPt), 
                                                        hGenInclusiveJetPtEtaStdBins->GetYaxis()->FindBin(ptVals[iPtBin+1]-deltaPt) ) );
        set1DStyle(hGenInclusiveJetEtaStdBinsPtBins[iPtBin], 1, true);
        hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->SetName(Form("hGenInclusiveJetEtaStdBinsPtBin_%d", iPtBin));
        // hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->Scale(1./hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->Integral("width"));
        rescaleHisto1D(hGenInclusiveJetEtaStdBinsPtBins[iPtBin]);
        hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->Scale(1./hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->Integral( hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->GetXaxis()->FindBin(-1.0), 
                                                                                                                 hGenInclusiveJetEtaStdBinsPtBins[iPtBin]->GetXaxis()->FindBin(1.0) ) );
    }
    hGenLeadJetEta->Scale(1./hGenLeadJetEta->Integral( hGenLeadJetEta->GetXaxis()->FindBin(-etaLabMax), 
                                                       hGenLeadJetEta->GetXaxis()->FindBin(etaLabMax) ) );
    hGenSubLeadJetEta->Scale(1./hGenSubLeadJetEta->Integral( hGenSubLeadJetEta->GetXaxis()->FindBin(-etaLabMax), 
                                                             hGenSubLeadJetEta->GetXaxis()->FindBin(etaLabMax) ) );

    // CM frame
    hGenInclusiveJetEtaCM = dynamic_cast<TH1D *>( hGenInclusiveJetPtEtaCM->ProjectionX("hGenInclusiveJetEtaCM", 
                                                  hGenInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                  hGenInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenInclusiveJetEtaCM, 1, true);
    hGenInclusiveJetEtaCM->Scale(1./hGenInclusiveJetEtaCM->Integral( hGenInclusiveJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                   hGenInclusiveJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hGenInclusiveJetEtaCMForward = new TH1D("hGenInclusiveJetEtaCMForward", "hGenInclusiveJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hGenInclusiveJetPtEtaCM->GetNbinsX() / 2, 0., hGenInclusiveJetPtEtaCM->GetXaxis()->GetXmax() );
    hGenInclusiveJetEtaCMBackward = new TH1D("hGenInclusiveJetEtaCMBackward", "hGenInclusiveJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hGenInclusiveJetPtEtaCM->GetNbinsX() / 2, 0., hGenInclusiveJetPtEtaCM->GetXaxis()->GetXmax() );
    hGenInclusiveJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hGenInclusiveJetEtaCMForward->Clone("hGenInclusiveJetForwardBackwardRatioCM") );
    hGenInclusiveJetForwardBackwardRatioCM->Divide(hGenInclusiveJetEtaCMBackward);
    set1DStyle(hGenInclusiveJetForwardBackwardRatioCM, 1);
    for ( int iPtBin = 0; iPtBin < nPtBins; ++iPtBin ) {
        hGenInclusiveJetEtaCMPtBins[iPtBin] = dynamic_cast<TH1D *>( hGenInclusiveJetPtEtaCM->ProjectionX(Form("hGenInclusiveJetEtaCMPtBin%d", iPtBin), 
                                                        hGenInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptVals[iPtBin]+deltaPt), 
                                                        hGenInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptVals[iPtBin+1]-deltaPt) ) );
        set1DStyle(hGenInclusiveJetEtaCMPtBins[iPtBin], 1, true);
        hGenInclusiveJetEtaCMPtBins[iPtBin]->SetName(Form("hGenInclusiveJetEtaCMPtBin_%d", iPtBin));
        hGenInclusiveJetEtaCMPtBins[iPtBin]->Scale(1./hGenInclusiveJetEtaCMPtBins[iPtBin]->Integral( hGenInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->FindBin(-etaCMMax), 
                                                                                                     hGenInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->FindBin(etaCMMax) ) );
    }

    hGenLeadJetEtaCM = dynamic_cast<TH1D *>( hGenLeadJetPtEtaCM->ProjectionX("hGenLeadJetEtaCM", 
                                             hGenLeadJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                             hGenLeadJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenLeadJetEtaCM, 1, true);
    hGenLeadJetEtaCM->Scale(1./hGenLeadJetEtaCM->Integral( hGenLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                   hGenLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hGenLeadJetEtaCMForward = new TH1D("hGenLeadJetEtaCMForward", "hGenLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hGenLeadJetPtEtaCM->GetNbinsX() / 2, 0., hGenLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hGenLeadJetEtaCMBackward = new TH1D("hGenLeadJetEtaCMBackward", "hGenLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hGenLeadJetPtEtaCM->GetNbinsX() / 2, 0., hGenLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hGenLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hGenLeadJetEtaCMForward->Clone("hGenLeadJetForwardBackwardRatioCM") );
    hGenLeadJetForwardBackwardRatioCM->Divide(hGenLeadJetEtaCMBackward);
    set1DStyle(hGenLeadJetForwardBackwardRatioCM, 1);

    hGenSubLeadJetEtaCM = dynamic_cast<TH1D *>( hGenSubLeadJetPtEtaCM->ProjectionX("hGenSubLeadJetEtaCM", 
                                                hGenSubLeadJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hGenSubLeadJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenSubLeadJetEtaCM, 1, true);
    hGenSubLeadJetEtaCM->Scale(1./hGenSubLeadJetEtaCM->Integral( hGenSubLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hGenSubLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hGenSubLeadJetEtaCMForward = new TH1D("hGenSubLeadJetEtaCMForward", "hGenSubLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}",
                                            hGenSubLeadJetPtEtaCM->GetNbinsX() / 2, 0., hGenSubLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hGenSubLeadJetEtaCMBackward = new TH1D("hGenSubLeadJetEtaCMBackward", "hGenSubLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hGenSubLeadJetPtEtaCM->GetNbinsX() / 2, 0., hGenSubLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hGenSubLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hGenSubLeadJetEtaCMForward->Clone("hGenSubLeadJetForwardBackwardRatioCM") );
    hGenSubLeadJetForwardBackwardRatioCM->Divide(hGenSubLeadJetEtaCMBackward);
    set1DStyle(hGenSubLeadJetForwardBackwardRatioCM, 1);

    // From selected dijets
    hGenDijetLeadJetEtaStdBins = dynamic_cast<TH1D *>( hGenDijetLeadPtEtaStdBins->ProjectionX("hGenDijetLeadJetEtaStdBins", 
                                                hGenDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                                hGenDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenDijetLeadJetEtaStdBins, 1);
    rescaleHisto1D(hGenDijetLeadJetEtaStdBins);
    hGenDijetLeadJetEtaStdBins->Scale(1./hGenDijetLeadJetEtaStdBins->Integral( hGenDijetLeadJetEtaStdBins->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hGenDijetLeadJetEtaStdBins->GetXaxis()->FindBin(etaCMMax) ) );
    hGenDijetLeadJetEtaCM = dynamic_cast<TH1D *>( hGenDijetLeadPtEtaCM->ProjectionX("hGenDijetLeadJetEtaCM", 
                                                hGenDijetLeadPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hGenDijetLeadPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenDijetLeadJetEtaCM, 1, true);
    hGenDijetLeadJetEtaCM->Scale(1./hGenDijetLeadJetEtaCM->Integral( hGenDijetLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hGenDijetLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hGenDijetLeadJetEtaCMForward = new TH1D("hGenDijetLeadJetEtaCMForward", "hGenDijetLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hGenDijetLeadPtEtaCM->GetNbinsX() / 2, 0., hGenDijetLeadPtEtaCM->GetXaxis()->GetXmax() );
    hGenDijetLeadJetEtaCMBackward = new TH1D("hGenDijetLeadJetEtaCMBackward", "hGenDijetLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hGenDijetLeadPtEtaCM->GetNbinsX() / 2, 0., hGenDijetLeadPtEtaCM->GetXaxis()->GetXmax() );
    hGenDijetLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hGenDijetLeadJetEtaCMForward->Clone("hGenDijetLeadJetForwardBackwardRatioCM") );
    hGenDijetLeadJetForwardBackwardRatioCM->Divide(hGenDijetLeadJetEtaCMBackward);
    set1DStyle(hGenDijetLeadJetForwardBackwardRatioCM, 1);

    hGenDijetSubLeadJetEtaStdBins = dynamic_cast<TH1D *>( hGenDijetSubLeadPtEtaStdBins->ProjectionX("hGenDijetSubLeadJetEtaStdBins", 
                                                hGenDijetSubLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                                hGenDijetSubLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenDijetSubLeadJetEtaStdBins, 1);
    rescaleHisto1D(hGenDijetSubLeadJetEtaStdBins);
    hGenDijetSubLeadJetEtaStdBins->Scale(1./hGenDijetSubLeadJetEtaStdBins->Integral( hGenDijetSubLeadJetEtaStdBins->GetXaxis()->FindBin(-etaCMMax), 
                                                                                       hGenDijetSubLeadJetEtaStdBins->GetXaxis()->FindBin(etaCMMax) ) );
    hGenDijetSubLeadJetEtaCM = dynamic_cast<TH1D *>( hGenDijetSubLeadPtEtaCM->ProjectionX("hGenDijetSubLeadJetEtaCM", 
                                                hGenDijetSubLeadPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hGenDijetSubLeadPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hGenDijetSubLeadJetEtaCM, 1, true);
    hGenDijetSubLeadJetEtaCM->Scale(1./hGenDijetSubLeadJetEtaCM->Integral( hGenDijetSubLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hGenDijetSubLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hGenDijetSubLeadJetEtaCMForward = new TH1D("hGenDijetSubLeadJetEtaCMForward", "hGenDijetSubLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hGenDijetSubLeadPtEtaCM->GetNbinsX() / 2, 0., hGenDijetSubLeadPtEtaCM->GetXaxis()->GetXmax() );
    hGenDijetSubLeadJetEtaCMBackward = new TH1D("hGenDijetSubLeadJetEtaCMBackward", "hGenDijetSubLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hGenDijetSubLeadPtEtaCM->GetNbinsX() / 2, 0., hGenDijetSubLeadPtEtaCM->GetXaxis()->GetXmax() );
    hGenDijetSubLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hGenDijetSubLeadJetEtaCMForward->Clone("hGenDijetSubLeadJetForwardBackwardRatioCM") );
    hGenDijetSubLeadJetForwardBackwardRatioCM->Divide(hGenDijetSubLeadJetEtaCMBackward);
    set1DStyle(hGenDijetSubLeadJetForwardBackwardRatioCM, 1);

    // Make forward/backward ratios
    recalculateFBRatioFromFullDistribution(hGenInclusiveJetForwardBackwardRatioCM, hGenInclusiveJetEtaCMForward, hGenInclusiveJetEtaCMBackward, hGenInclusiveJetEtaCM);
    recalculateFBRatioFromFullDistribution(hGenLeadJetForwardBackwardRatioCM, hGenLeadJetEtaCMForward, hGenLeadJetEtaCMBackward, hGenLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hGenSubLeadJetForwardBackwardRatioCM, hGenSubLeadJetEtaCMForward, hGenSubLeadJetEtaCMBackward, hGenSubLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hGenDijetLeadJetForwardBackwardRatioCM, hGenDijetLeadJetEtaCMForward, hGenDijetLeadJetEtaCMBackward, hGenDijetLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hGenDijetSubLeadJetForwardBackwardRatioCM, hGenDijetSubLeadJetEtaCMForward, hGenDijetSubLeadJetEtaCMBackward, hGenDijetSubLeadJetEtaCM);

    // Reco level

    // Lab frame
    hRecoInclusiveJetEta = dynamic_cast<TH1D *>( hRecoInclusiveJetPtEta->ProjectionX("hRecoInclusiveJetEta", 
                                        hRecoInclusiveJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hRecoInclusiveJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    hRecoInclusiveJetEtaStdBins = dynamic_cast<TH1D *>( hRecoInclusiveJetPtEtaStdBins->ProjectionX("hRecoInclusiveJetEtaStdBins", 
                                        hRecoInclusiveJetPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                        hRecoInclusiveJetPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    hRecoLeadJetEta = dynamic_cast<TH1D *>( hRecoLeadJetPtEta->ProjectionX("hRecoLeadJetEta", 
                                        hRecoLeadJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hRecoLeadJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    hRecoSubLeadJetEta = dynamic_cast<TH1D *>( hRecoSubLeadJetPtEta->ProjectionX("hRecoSubLeadJetEta", 
                                        hRecoSubLeadJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hRecoSubLeadJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoInclusiveJetEta, 0, true);
    set1DStyle(hRecoInclusiveJetEtaStdBins, 0);
    set1DStyle(hRecoLeadJetEta, 0, true);
    set1DStyle(hRecoSubLeadJetEta, 0, true);
    hRecoInclusiveJetEta->Scale(1./hRecoInclusiveJetEta->Integral( hRecoInclusiveJetEta->GetXaxis()->FindBin(-1.0), 
                                                                 hRecoInclusiveJetEta->GetXaxis()->FindBin(1.0) ) );
    // hRecoInclusiveJetEtaStdBins->Scale(1./hRecoInclusiveJetEtaStdBins->Integral("width"));
    rescaleHisto1D(hRecoInclusiveJetEtaStdBins);
    hRecoInclusiveJetEtaStdBins->Scale(1./hRecoInclusiveJetEtaStdBins->Integral( hRecoInclusiveJetEtaStdBins->GetXaxis()->FindBin(-1.0), 
                                                                               hRecoInclusiveJetEtaStdBins->GetXaxis()->FindBin(1.0) ) );
    for (int iPtBin = 0; iPtBin < nPtBins; ++iPtBin) {
        hRecoInclusiveJetEtaStdBinsPtBins[iPtBin] = dynamic_cast<TH1D *>( hRecoInclusiveJetPtEtaStdBins->ProjectionX(Form("hRecoInclusiveJetEtaStdBinsPtBin%d", iPtBin), 
                                                        hRecoInclusiveJetPtEtaStdBins->GetYaxis()->FindBin(ptVals[iPtBin]+deltaPt), 
                                                        hRecoInclusiveJetPtEtaStdBins->GetYaxis()->FindBin(ptVals[iPtBin+1]-deltaPt) ) );
        set1DStyle(hRecoInclusiveJetEtaStdBinsPtBins[iPtBin], 0, true);
        hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->SetName(Form("hRecoInclusiveJetEtaStdBinsPtBin_%d", iPtBin));
        // hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->Scale(1./hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->Integral("width"));
        rescaleHisto1D(hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]);
        hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->Scale(1./hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->Integral( hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->GetXaxis()->FindBin(-1.0), 
                                                                                                                 hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->GetXaxis()->FindBin(1.0) ) );
    }
    hRecoLeadJetEta->Scale(1./hRecoLeadJetEta->Integral( hRecoLeadJetEta->GetXaxis()->FindBin(-etaLabMax), 
                                                       hRecoLeadJetEta->GetXaxis()->FindBin(etaLabMax) ) );
    hRecoSubLeadJetEta->Scale(1./hRecoSubLeadJetEta->Integral( hRecoSubLeadJetEta->GetXaxis()->FindBin(-etaLabMax), 
                                                             hRecoSubLeadJetEta->GetXaxis()->FindBin(etaLabMax) ) );

    hInclusiveJetReco2GenEta = dynamic_cast<TH1D *>( hRecoInclusiveJetEta->Clone("hInclusiveJetReco2GenEta") );
    hInclusiveJetReco2GenEta->Divide(hGenInclusiveJetEta);
    set1DStyle(hInclusiveJetReco2GenEta, 0);
    hInclusiveJetReco2GenEtaStdBins = dynamic_cast<TH1D *>( hRecoInclusiveJetEtaStdBins->Clone("hInclusiveJetReco2GenEtaStdBins") );
    hInclusiveJetReco2GenEtaStdBins->Divide(hGenInclusiveJetEtaStdBins);
    set1DStyle(hInclusiveJetReco2GenEtaStdBins, 0); 
    for (int iPtBin = 0; iPtBin < nPtBins; ++iPtBin) {
        hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin] = dynamic_cast<TH1D *>( hRecoInclusiveJetEtaStdBinsPtBins[iPtBin]->Clone(Form("hInclusiveJetReco2GenEtaStdBinsPtBin%d", iPtBin)) );
        hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->Divide(hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin], hGenInclusiveJetEtaStdBinsPtBins[iPtBin], 1., 1., "B");
        set1DStyle(hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin], 0);
    }
    hLeadJetReco2GenEta = dynamic_cast<TH1D *>( hRecoLeadJetEta->Clone("hLeadJetReco2GenEta") );
    hLeadJetReco2GenEta->Divide(hGenLeadJetEta);
    set1DStyle(hLeadJetReco2GenEta, 0);
    hSubLeadJetReco2GenEta = dynamic_cast<TH1D *>( hRecoSubLeadJetEta->Clone("hSubLeadJetReco2GenEta") );
    hSubLeadJetReco2GenEta->Divide(hGenSubLeadJetEta);
    set1DStyle(hSubLeadJetReco2GenEta, 0);

    // CM frame
    hRecoInclusiveJetEtaCM = dynamic_cast<TH1D *>( hRecoInclusiveJetPtEtaCM->ProjectionX("hRecoInclusiveJetEtaCM", 
                                                   hRecoInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                   hRecoInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoInclusiveJetEtaCM, 0, true);
    hRecoInclusiveJetEtaCM->Scale(1./hRecoInclusiveJetEtaCM->Integral( hRecoInclusiveJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                   hRecoInclusiveJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoInclusiveJetEtaCMForward = new TH1D("hRecoInclusiveJetEtaCMForward", "hRecoInclusiveJetEtaCMForward;#eta_{CM};dN/d#eta_{CM}", 
                                            hRecoInclusiveJetPtEtaCM->GetNbinsX() / 2, 0., hRecoInclusiveJetPtEtaCM->GetXaxis()->GetXmax() );
    hRecoInclusiveJetEtaCMBackward = new TH1D("hRecoInclusiveJetEtaCMBackward", "hRecoInclusiveJetEtaCMBackward;#eta_{CM};dN/d#eta_{CM}", 
                                             hRecoInclusiveJetPtEtaCM->GetNbinsX() / 2, 0., hRecoInclusiveJetPtEtaCM->GetXaxis()->GetXmax() );
    hRecoInclusiveJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRecoInclusiveJetEtaCMForward->Clone("hRecoInclusiveJetForwardBackwardRatioCM") );
    hRecoInclusiveJetForwardBackwardRatioCM->Divide(hRecoInclusiveJetEtaCMBackward);
    set1DStyle(hRecoInclusiveJetForwardBackwardRatioCM, 0, true);
    for ( int iPtBin = 0; iPtBin < nPtBins; ++iPtBin ) {
        hRecoInclusiveJetEtaCMPtBins[iPtBin] = dynamic_cast<TH1D *>( hRecoInclusiveJetPtEtaCM->ProjectionX(Form("hRecoInclusiveJetEtaCMPtBin%d", iPtBin), 
                                                        hRecoInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptVals[iPtBin]+deltaPt), 
                                                        hRecoInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptVals[iPtBin+1]-deltaPt) ) );
        set1DStyle(hRecoInclusiveJetEtaCMPtBins[iPtBin], 0, true);
        hRecoInclusiveJetEtaCMPtBins[iPtBin]->SetName(Form("hRecoInclusiveJetEtaCMPtBin_%d", iPtBin));
        hRecoInclusiveJetEtaCMPtBins[iPtBin]->Scale(1./hRecoInclusiveJetEtaCMPtBins[iPtBin]->Integral( hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->FindBin(-etaCMMax), 
                                                                                                     hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->FindBin(etaCMMax) ) );
    }

    hRecoLeadJetEtaCM = dynamic_cast<TH1D *>( hRecoLeadJetPtEtaCM->ProjectionX("hRecoLeadJetEtaCM", 
                                             hRecoLeadJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                             hRecoLeadJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoLeadJetEtaCM, 0, true);
    hRecoLeadJetEtaCM->Scale(1./hRecoLeadJetEtaCM->Integral( hRecoLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                   hRecoLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoLeadJetEtaCMForward = new TH1D("hRecoLeadJetEtaCMForward", "hRecoLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}",
                                            hRecoLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRecoLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRecoLeadJetEtaCMBackward = new TH1D("hRecoLeadJetEtaCMBackward", "hRecoLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRecoLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRecoLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRecoLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRecoLeadJetEtaCMForward->Clone("hRecoLeadJetForwardBackwardRatioCM") );
    hRecoLeadJetForwardBackwardRatioCM->Divide(hRecoLeadJetEtaCMBackward);
    set1DStyle(hRecoLeadJetForwardBackwardRatioCM, 0, true);

    hRecoSubLeadJetEtaCM = dynamic_cast<TH1D *>( hRecoSubLeadJetPtEtaCM->ProjectionX("hRecoSubLeadJetEtaCM", 
                                                hRecoSubLeadJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hRecoSubLeadJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoSubLeadJetEtaCM, 0, true);
    hRecoSubLeadJetEtaCM->Scale(1./hRecoSubLeadJetEtaCM->Integral( hRecoSubLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hRecoSubLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoSubLeadJetEtaCMForward = new TH1D("hRecoSubLeadJetEtaCMForward", "hRecoSubLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}",
                                            hRecoSubLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRecoSubLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRecoSubLeadJetEtaCMBackward = new TH1D("hRecoSubLeadJetEtaCMBackward", "hRecoSubLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRecoSubLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRecoSubLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRecoSubLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRecoSubLeadJetEtaCMForward->Clone("hRecoSubLeadJetForwardBackwardRatioCM") );
    hRecoSubLeadJetForwardBackwardRatioCM->Divide(hRecoSubLeadJetEtaCMBackward);
    set1DStyle(hRecoSubLeadJetForwardBackwardRatioCM, 0, true);

    // From selected dijets
    hRecoDijetLeadJetEtaStdBins = dynamic_cast<TH1D *>( hRecoDijetLeadPtEtaStdBins->ProjectionX("hRecoDijetLeadJetEtaStdBins", 
                                                hRecoDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                                hRecoDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoDijetLeadJetEtaStdBins, 0);
    rescaleHisto1D(hRecoDijetLeadJetEtaStdBins);
    hRecoDijetLeadJetEtaStdBins->Scale(1./hRecoDijetLeadJetEtaStdBins->Integral( hRecoDijetLeadJetEtaStdBins->GetXaxis()->FindBin(-etaCMMax), 
                                                                                       hRecoDijetLeadJetEtaStdBins->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoDijetLeadJetEtaCM = dynamic_cast<TH1D *>( hRecoDijetLeadPtEtaCM->ProjectionX("hRecoDijetLeadJetEtaCM", 
                                                hRecoDijetLeadPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hRecoDijetLeadPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoDijetLeadJetEtaCM, 0, true);
    hRecoDijetLeadJetEtaCM->Scale(1./hRecoDijetLeadJetEtaCM->Integral( hRecoDijetLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hRecoDijetLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoDijetLeadJetEtaCMForward = new TH1D("hRecoDijetLeadJetEtaCMForward", "hRecoDijetLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRecoDijetLeadPtEtaCM->GetNbinsX() / 2, 0., hRecoDijetLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRecoDijetLeadJetEtaCMBackward = new TH1D("hRecoDijetLeadJetEtaCMBackward", "hRecoDijetLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRecoDijetLeadPtEtaCM->GetNbinsX() / 2, 0., hRecoDijetLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRecoDijetLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRecoDijetLeadJetEtaCMForward->Clone("hRecoDijetLeadJetForwardBackwardRatioCM") );
    hRecoDijetLeadJetForwardBackwardRatioCM->Divide(hRecoDijetLeadJetEtaCMBackward);
    set1DStyle(hRecoDijetLeadJetForwardBackwardRatioCM, 0, true);

    hRecoDijetSubLeadJetEtaStdBins = dynamic_cast<TH1D *>( hRecoDijetSubLeadPtEtaStdBins->ProjectionX("hRecoDijetSubLeadJetEtaStdBins", 
                                                hRecoDijetSubLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                                hRecoDijetSubLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoDijetSubLeadJetEtaStdBins, 0);
    rescaleHisto1D(hRecoDijetSubLeadJetEtaStdBins);
    hRecoDijetSubLeadJetEtaStdBins->Scale(1./hRecoDijetSubLeadJetEtaStdBins->Integral( hRecoDijetSubLeadJetEtaStdBins->GetXaxis()->FindBin(-etaCMMax), 
                                                                                       hRecoDijetSubLeadJetEtaStdBins->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoDijetSubLeadJetEtaCM = dynamic_cast<TH1D *>( hRecoDijetSubLeadPtEtaCM->ProjectionX("hRecoDijetSubLeadJetEtaCM", 
                                                hRecoDijetSubLeadPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hRecoDijetSubLeadPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRecoDijetSubLeadJetEtaCM, 0, true);
    hRecoDijetSubLeadJetEtaCM->Scale(1./hRecoDijetSubLeadJetEtaCM->Integral( hRecoDijetSubLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hRecoDijetSubLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRecoDijetSubLeadJetEtaCMForward = new TH1D("hRecoDijetSubLeadJetEtaCMForward", "hRecoDijetSubLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRecoDijetSubLeadPtEtaCM->GetNbinsX() / 2, 0., hRecoDijetSubLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRecoDijetSubLeadJetEtaCMBackward = new TH1D("hRecoDijetSubLeadJetEtaCMBackward", "hRecoDijetSubLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRecoDijetSubLeadPtEtaCM->GetNbinsX() / 2, 0., hRecoDijetSubLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRecoDijetSubLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRecoDijetSubLeadJetEtaCMForward->Clone("hRecoDijetSubLeadJetForwardBackwardRatioCM") );
    hRecoDijetSubLeadJetForwardBackwardRatioCM->Divide(hRecoDijetSubLeadJetEtaCMBackward);
    set1DStyle(hRecoDijetSubLeadJetForwardBackwardRatioCM, 0, true);

    // Make forward/backward ratios
    recalculateFBRatioFromFullDistribution(hRecoInclusiveJetForwardBackwardRatioCM, hRecoInclusiveJetEtaCMForward, hRecoInclusiveJetEtaCMBackward, hRecoInclusiveJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRecoLeadJetForwardBackwardRatioCM, hRecoLeadJetEtaCMForward, hRecoLeadJetEtaCMBackward, hRecoLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRecoSubLeadJetForwardBackwardRatioCM, hRecoSubLeadJetEtaCMForward, hRecoSubLeadJetEtaCMBackward, hRecoSubLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRecoDijetLeadJetForwardBackwardRatioCM, hRecoDijetLeadJetEtaCMForward, hRecoDijetLeadJetEtaCMBackward, hRecoDijetLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRecoDijetSubLeadJetForwardBackwardRatioCM, hRecoDijetSubLeadJetEtaCMForward, hRecoDijetSubLeadJetEtaCMBackward, hRecoDijetSubLeadJetEtaCM);

    hInclusiveJetReco2GenEtaCM = dynamic_cast<TH1D *>( hRecoInclusiveJetEtaCM->Clone("hInclusiveJetReco2GenEtaCM") );
    hInclusiveJetReco2GenEtaCM->Divide(hGenInclusiveJetEtaCM);
    set1DStyle(hInclusiveJetReco2GenEtaCM, 0);
    hLeadJetReco2GenEtaCM = dynamic_cast<TH1D *>( hRecoLeadJetEtaCM->Clone("hLeadJetReco2GenEtaCM") );
    hLeadJetReco2GenEtaCM->Divide(hGenLeadJetEtaCM);
    set1DStyle(hLeadJetReco2GenEtaCM, 0);
    hSubLeadJetReco2GenEtaCM = dynamic_cast<TH1D *>( hRecoSubLeadJetEtaCM->Clone("hSubLeadJetReco2GenEtaCM") );
    hSubLeadJetReco2GenEtaCM->Divide(hGenSubLeadJetEtaCM);
    set1DStyle(hSubLeadJetReco2GenEtaCM, 0);
    hDijetLeadJetReco2GenEtaCM = dynamic_cast<TH1D *>( hRecoDijetLeadJetEtaCM->Clone("hDijetLeadJetReco2GenEtaCM") );
    hDijetLeadJetReco2GenEtaCM->Divide(hGenDijetLeadJetEtaCM);
    set1DStyle(hDijetLeadJetReco2GenEtaCM, 0);
    hDijetSubLeadJetReco2GenEtaCM = dynamic_cast<TH1D *>( hRecoDijetSubLeadJetEtaCM->Clone("hDijetSubLeadJetReco2GenEtaCM") );
    hDijetSubLeadJetReco2GenEtaCM->Divide(hGenDijetSubLeadJetEtaCM);
    set1DStyle(hDijetSubLeadJetReco2GenEtaCM, 0);
    hDijetLeadJetReco2GenEtaStdBins = dynamic_cast<TH1D *>( hRecoDijetLeadJetEtaStdBins->Clone("hDijetLeadJetReco2GenEtaStdBins") );
    hDijetLeadJetReco2GenEtaStdBins->Divide(hGenDijetLeadJetEtaStdBins);
    set1DStyle(hDijetLeadJetReco2GenEtaStdBins, 0);
    hDijetSubLeadJetReco2GenEtaStdBins = dynamic_cast<TH1D *>( hRecoDijetSubLeadJetEtaStdBins->Clone("hDijetSubLeadJetReco2GenEtaStdBins") );
    hDijetSubLeadJetReco2GenEtaStdBins->Divide(hGenDijetSubLeadJetEtaStdBins);
    set1DStyle(hDijetSubLeadJetReco2GenEtaStdBins, 0);

    // Ref level

    // Lab frame
    hRefInclusiveJetEta = dynamic_cast<TH1D *>( hRefInclusiveJetPtEta->ProjectionX("hRefInclusiveJetEta", 
                                        hRefInclusiveJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hRefInclusiveJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    hRefLeadJetEta = dynamic_cast<TH1D *>( hRefLeadJetPtEta->ProjectionX("hRefLeadJetEta", 
                                        hRefLeadJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hRefLeadJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    hRefSubLeadJetEta = dynamic_cast<TH1D *>( hRefSubLeadJetPtEta->ProjectionX("hRefSubLeadJetEta", 
                                        hRefSubLeadJetPtEta->GetYaxis()->FindBin(ptLow), 
                                        hRefSubLeadJetPtEta->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefInclusiveJetEta, 2, true);
    set1DStyle(hRefLeadJetEta, 2, true);
    set1DStyle(hRefSubLeadJetEta, 2, true);
    hRefInclusiveJetEta->Scale(1./hRefInclusiveJetEta->Integral( hRefInclusiveJetEta->GetXaxis()->FindBin(-1.0), 
                                                                 hRefInclusiveJetEta->GetXaxis()->FindBin(1.0) ) );
    hRefLeadJetEta->Scale(1./hRefLeadJetEta->Integral( hRefLeadJetEta->GetXaxis()->FindBin(-etaLabMax), 
                                                       hRefLeadJetEta->GetXaxis()->FindBin(etaLabMax) ) );
    hRefSubLeadJetEta->Scale(1./hRefSubLeadJetEta->Integral( hRefSubLeadJetEta->GetXaxis()->FindBin(-etaLabMax), 
                                                             hRefSubLeadJetEta->GetXaxis()->FindBin(etaLabMax) ) );

    hInclusiveJetRef2GenEta = dynamic_cast<TH1D *>( hRefInclusiveJetEta->Clone("hInclusiveJetRef2GenEta") );
    hInclusiveJetRef2GenEta->Divide(hGenInclusiveJetEta);
    set1DStyle(hInclusiveJetRef2GenEta, 2);
    hLeadJetRef2GenEta = dynamic_cast<TH1D *>( hRefLeadJetEta->Clone("hLeadJetRef2GenEta") );
    hLeadJetRef2GenEta->Divide(hGenLeadJetEta);
    set1DStyle(hLeadJetRef2GenEta, 2);
    hSubLeadJetRef2GenEta = dynamic_cast<TH1D *>( hRefSubLeadJetEta->Clone("hSubLeadJetRef2GenEta") );
    hSubLeadJetRef2GenEta->Divide(hGenSubLeadJetEta);
    set1DStyle(hSubLeadJetRef2GenEta, 2);

    // CM frame
    hRefInclusiveJetEtaCM = dynamic_cast<TH1D *>( hRefInclusiveJetPtEtaCM->ProjectionX("hRefInclusiveJetEtaCM", 
                                                  hRefInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                  hRefInclusiveJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefInclusiveJetEtaCM, 2, true);
    hRefInclusiveJetEtaCM->Scale(1./hRefInclusiveJetEtaCM->Integral( hRefInclusiveJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                   hRefInclusiveJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRefInclusiveJetEtaCMForward = new TH1D("hRefInclusiveJetEtaCMForward", "hRefInclusiveJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRefInclusiveJetPtEtaCM->GetNbinsX() / 2, 0., hRefInclusiveJetPtEtaCM->GetXaxis()->GetXmax() );
    hRefInclusiveJetEtaCMBackward = new TH1D("hRefInclusiveJetEtaCMBackward", "hRefInclusiveJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRefInclusiveJetPtEtaCM->GetNbinsX() / 2, 0., hRefInclusiveJetPtEtaCM->GetXaxis()->GetXmax() );
    hRefInclusiveJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRefInclusiveJetEtaCMForward->Clone("hRefInclusiveJetForwardBackwardRatioCM") );
    hRefInclusiveJetForwardBackwardRatioCM->Divide(hRefInclusiveJetEtaCMBackward);
    set1DStyle(hRefInclusiveJetForwardBackwardRatioCM, 3, true);

    hRefLeadJetEtaCM = dynamic_cast<TH1D *>( hRefLeadJetPtEtaCM->ProjectionX("hRefLeadJetEtaCM", 
                                             hRefLeadJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                             hRefLeadJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefLeadJetEtaCM, 2, true);
    hRefLeadJetEtaCM->Scale(1./hRefLeadJetEtaCM->Integral( hRefLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                   hRefLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRefLeadJetEtaCMForward = new TH1D("hRefLeadJetEtaCMForward", "hRefLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRefLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRefLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRefLeadJetEtaCMBackward = new TH1D("hRefLeadJetEtaCMBackward", "hRefLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRefLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRefLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRefLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRefLeadJetEtaCMForward->Clone("hRefLeadJetForwardBackwardRatioCM") );
    hRefLeadJetForwardBackwardRatioCM->Divide(hRefLeadJetEtaCMBackward);
    set1DStyle(hRefLeadJetForwardBackwardRatioCM, 2, true);

    hRefSubLeadJetEtaCM = dynamic_cast<TH1D *>( hRefSubLeadJetPtEtaCM->ProjectionX("hRefSubLeadJetEtaCM", 
                                                hRefSubLeadJetPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hRefSubLeadJetPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefSubLeadJetEtaCM, 2, true);
    hRefSubLeadJetEtaCM->Scale(1./hRefSubLeadJetEtaCM->Integral( hRefSubLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hRefSubLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRefSubLeadJetEtaCMForward = new TH1D("hRefSubLeadJetEtaCMForward", "hRefSubLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRefSubLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRefSubLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRefSubLeadJetEtaCMBackward = new TH1D("hRefSubLeadJetEtaCMBackward", "hRefSubLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRefSubLeadJetPtEtaCM->GetNbinsX() / 2, 0., hRefSubLeadJetPtEtaCM->GetXaxis()->GetXmax() );
    hRefSubLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRefSubLeadJetEtaCMForward->Clone("hRefSubLeadJetForwardBackwardRatioCM") );
    hRefSubLeadJetForwardBackwardRatioCM->Divide(hRefSubLeadJetEtaCMBackward);
    set1DStyle(hRefSubLeadJetForwardBackwardRatioCM, 2, true);

    // From selected dijets
    hRefDijetLeadJetEtaStdBins = dynamic_cast<TH1D *>( hRefDijetLeadPtEtaStdBins->ProjectionX("hRefDijetLeadJetEtaStdBins", 
                                                hRefDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                                hRefDijetLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefDijetLeadJetEtaStdBins, 2);
    rescaleHisto1D(hRefDijetLeadJetEtaStdBins);
    hRefDijetLeadJetEtaStdBins->Scale(1./hRefDijetLeadJetEtaStdBins->Integral( hRefDijetLeadJetEtaStdBins->GetXaxis()->FindBin(-etaCMMax), 
                                                                                       hRefDijetLeadJetEtaStdBins->GetXaxis()->FindBin(etaCMMax) ) );
    hRefDijetLeadJetEtaCM = dynamic_cast<TH1D *>( hRefDijetLeadPtEtaCM->ProjectionX("hRefDijetLeadJetEtaCM", 
                                                hRefDijetLeadPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hRefDijetLeadPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefDijetLeadJetEtaCM, 2, true);
    hRefDijetLeadJetEtaCM->Scale(1./hRefDijetLeadJetEtaCM->Integral( hRefDijetLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hRefDijetLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRefDijetLeadJetEtaCMForward = new TH1D("hRefDijetLeadJetEtaCMForward", "hRefDijetLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRefDijetLeadPtEtaCM->GetNbinsX() / 2, 0., hRefDijetLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRefDijetLeadJetEtaCMBackward = new TH1D("hRefDijetLeadJetEtaCMBackward", "hRefDijetLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRefDijetLeadPtEtaCM->GetNbinsX() / 2, 0., hRefDijetLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRefDijetLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRefDijetLeadJetEtaCMForward->Clone("hRefDijetLeadJetForwardBackwardRatioCM") );
    hRefDijetLeadJetForwardBackwardRatioCM->Divide(hRefDijetLeadJetEtaCMBackward);
    set1DStyle(hRefDijetLeadJetForwardBackwardRatioCM, 2, true);

    hRefDijetSubLeadJetEtaStdBins = dynamic_cast<TH1D *>( hRefDijetSubLeadPtEtaStdBins->ProjectionX("hRefDijetSubLeadJetEtaStdBins", 
                                                hRefDijetSubLeadPtEtaStdBins->GetYaxis()->FindBin(ptLow), 
                                                hRefDijetSubLeadPtEtaStdBins->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefDijetSubLeadJetEtaStdBins, 2);
    rescaleHisto1D(hRefDijetSubLeadJetEtaStdBins);
    hRefDijetSubLeadJetEtaStdBins->Scale(1./hRefDijetSubLeadJetEtaStdBins->Integral( hRefDijetSubLeadJetEtaStdBins->GetXaxis()->FindBin(-etaCMMax), 
                                                                                       hRefDijetSubLeadJetEtaStdBins->GetXaxis()->FindBin(etaCMMax) ) );
    hRefDijetSubLeadJetEtaCM = dynamic_cast<TH1D *>( hRefDijetSubLeadPtEtaCM->ProjectionX("hRefDijetSubLeadJetEtaCM", 
                                                hRefDijetSubLeadPtEtaCM->GetYaxis()->FindBin(ptLow), 
                                                hRefDijetSubLeadPtEtaCM->GetYaxis()->FindBin(ptHigh)) );
    set1DStyle(hRefDijetSubLeadJetEtaCM, 2, true);
    hRefDijetSubLeadJetEtaCM->Scale(1./hRefDijetSubLeadJetEtaCM->Integral( hRefDijetSubLeadJetEtaCM->GetXaxis()->FindBin(-etaCMMax), 
                                                                 hRefDijetSubLeadJetEtaCM->GetXaxis()->FindBin(etaCMMax) ) );
    hRefDijetSubLeadJetEtaCMForward = new TH1D("hRefDijetSubLeadJetEtaCMForward", "hRefDijetSubLeadJetEtaCMForward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                            hRefDijetSubLeadPtEtaCM->GetNbinsX() / 2, 0., hRefDijetSubLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRefDijetSubLeadJetEtaCMBackward = new TH1D("hRefDijetSubLeadJetEtaCMBackward", "hRefDijetSubLeadJetEtaCMBackward;#eta^{jet}_{CM};dN/d#eta^{jet}_{CM}", 
                                             hRefDijetSubLeadPtEtaCM->GetNbinsX() / 2, 0., hRefDijetSubLeadPtEtaCM->GetXaxis()->GetXmax() );
    hRefDijetSubLeadJetForwardBackwardRatioCM = dynamic_cast<TH1D *>( hRefDijetSubLeadJetEtaCMForward->Clone("hRefDijetSubLeadJetForwardBackwardRatioCM") );
    hRefDijetSubLeadJetForwardBackwardRatioCM->Divide(hRefDijetSubLeadJetEtaCMBackward);
    set1DStyle(hRefDijetSubLeadJetForwardBackwardRatioCM, 2, true);

    // Make forward/backward ratios
    recalculateFBRatioFromFullDistribution(hRefInclusiveJetForwardBackwardRatioCM, hRefInclusiveJetEtaCMForward, hRefInclusiveJetEtaCMBackward, hRefInclusiveJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRefLeadJetForwardBackwardRatioCM, hRefLeadJetEtaCMForward, hRefLeadJetEtaCMBackward, hRefLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRefSubLeadJetForwardBackwardRatioCM, hRefSubLeadJetEtaCMForward, hRefSubLeadJetEtaCMBackward, hRefSubLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRefDijetLeadJetForwardBackwardRatioCM, hRefDijetLeadJetEtaCMForward, hRefDijetLeadJetEtaCMBackward, hRefDijetLeadJetEtaCM);
    recalculateFBRatioFromFullDistribution(hRefDijetSubLeadJetForwardBackwardRatioCM, hRefDijetSubLeadJetEtaCMForward, hRefDijetSubLeadJetEtaCMBackward, hRefDijetSubLeadJetEtaCM);

    hInclusiveJetRef2GenEtaCM = dynamic_cast<TH1D *>( hRefInclusiveJetEtaCM->Clone("hInclusiveJetRef2GenEtaCM") );
    hInclusiveJetRef2GenEtaCM->Divide(hGenInclusiveJetEtaCM);
    set1DStyle(hInclusiveJetRef2GenEtaCM, 2);
    hLeadJetRef2GenEtaCM = dynamic_cast<TH1D *>( hRefLeadJetEtaCM->Clone("hLeadJetRef2GenEtaCM") );
    hLeadJetRef2GenEtaCM->Divide(hGenLeadJetEtaCM);
    set1DStyle(hLeadJetRef2GenEtaCM, 2);
    hSubLeadJetRef2GenEtaCM = dynamic_cast<TH1D *>( hRefSubLeadJetEtaCM->Clone("hSubLeadJetRef2GenEtaCM") );
    hSubLeadJetRef2GenEtaCM->Divide(hGenSubLeadJetEtaCM);
    set1DStyle(hSubLeadJetRef2GenEtaCM, 2);

    hDijetLeadJetRef2GenEtaCM = dynamic_cast<TH1D *>( hRefDijetLeadJetEtaCM->Clone("hDijetLeadJetRef2GenEtaCM") );
    hDijetLeadJetRef2GenEtaCM->Divide(hGenDijetLeadJetEtaCM);
    set1DStyle(hDijetLeadJetRef2GenEtaCM, 2);
    hDijetSubLeadJetRef2GenEtaCM = dynamic_cast<TH1D *>( hRefDijetSubLeadJetEtaCM->Clone("hDijetSubLeadJetRef2GenEtaCM") );
    hDijetSubLeadJetRef2GenEtaCM->Divide(hGenDijetSubLeadJetEtaCM);
    set1DStyle(hDijetSubLeadJetRef2GenEtaCM, 2);
    hDijetLeadJetRef2GenEtaStdBins = dynamic_cast<TH1D *>( hRefDijetLeadJetEtaStdBins->Clone("hDijetLeadJetRef2GenEtaStdBins") );
    hDijetLeadJetRef2GenEtaStdBins->Divide(hGenDijetLeadJetEtaStdBins);
    set1DStyle(hDijetLeadJetRef2GenEtaStdBins, 2);
    hDijetSubLeadJetRef2GenEtaStdBins = dynamic_cast<TH1D *>( hRefDijetSubLeadJetEtaStdBins->Clone("hDijetSubLeadJetRef2GenEtaStdBins") );
    hDijetSubLeadJetRef2GenEtaStdBins->Divide(hGenDijetSubLeadJetEtaStdBins);
    set1DStyle(hDijetSubLeadJetRef2GenEtaStdBins, 2);


    /////////////////////////////////////
    //             Plotting            //
    /////////////////////////////////////

    double meanRecoEtaCM = 0.;
    double meanGenEtaCM = 0.;
    double meanRefEtaCM = 0.;

    TLegend *leg;
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.04);

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();

    // Overlay inclusive reco, ref and gen distributions in lab frame
    hGenInclusiveJetEta->Draw();
    hRecoInclusiveJetEta->Draw("same");
    hRefInclusiveJetEta->Draw("same");
    hGenInclusiveJetEta->GetXaxis()->SetTitle("#eta^{Inclusive}");
    hGenInclusiveJetEta->GetYaxis()->SetTitle("dN/d#eta^{Inclusive}");
    hGenInclusiveJetEta->GetYaxis()->SetRangeUser(0., 0.1);
    hGenInclusiveJetEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenInclusiveJetEta, "Gen", "p");
    leg->AddEntry(hRecoInclusiveJetEta, "Reco", "p");
    leg->AddEntry(hRefInclusiveJetEta, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_etaLab_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );

    if (leg) { delete leg; leg = nullptr; }
    
    // Overlay lead jet reco, ref and gen distributions in lab frame
    hGenLeadJetEta->Draw();
    hRecoLeadJetEta->Draw("same");
    hRefLeadJetEta->Draw("same");
    hGenLeadJetEta->GetXaxis()->SetTitle("#eta^{Lead}");
    hGenLeadJetEta->GetYaxis()->SetTitle("dN/d#eta^{Lead}");
    hGenLeadJetEta->GetYaxis()->SetRangeUser(0., 0.1);
    hGenLeadJetEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenLeadJetEta, "Gen", "p");
    leg->AddEntry(hRecoLeadJetEta, "Reco", "p");
    leg->AddEntry(hRefLeadJetEta, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_leadJet_etaLab_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }
    
    // Overlay sublead jet reco, ref and gen distributions in lab frame
    hGenSubLeadJetEta->Draw();
    hRecoSubLeadJetEta->Draw("same");
    hRefSubLeadJetEta->Draw("same");
    hGenSubLeadJetEta->GetXaxis()->SetTitle("#eta^{SubLead}");
    hGenSubLeadJetEta->GetYaxis()->SetTitle("dN/d#eta^{SubLead}");
    hGenSubLeadJetEta->GetYaxis()->SetRangeUser(0., 0.1);
    hGenSubLeadJetEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenSubLeadJetEta, "Gen", "p");
    leg->AddEntry(hRecoSubLeadJetEta, "Reco", "p");
    leg->AddEntry(hRefSubLeadJetEta, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_subLeadJet_etaLab_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of inclusive reco and ref to gen ratio in lab frame
    hInclusiveJetReco2GenEta->Draw();
    hInclusiveJetRef2GenEta->Draw("same");
    hInclusiveJetReco2GenEta->GetXaxis()->SetTitle("#eta^{Inclusive}");
    hInclusiveJetReco2GenEta->GetYaxis()->SetTitle("Reco / Gen");
    hInclusiveJetReco2GenEta->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hInclusiveJetReco2GenEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hInclusiveJetReco2GenEta, "Reco / Gen", "p");
    leg->AddEntry(hInclusiveJetRef2GenEta, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_etaLab_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay inclusive reco and gen eta distributions in lab frame with the standard eta JEC bins
    hGenInclusiveJetEtaStdBins->Draw();
    hRecoInclusiveJetEtaStdBins->Draw("same");
    hGenInclusiveJetEtaStdBins->GetXaxis()->SetTitle("#eta^{Inclusive}_{std bins}");
    hGenInclusiveJetEtaStdBins->GetYaxis()->SetTitle("dN/d#eta^{Inclusive}_{std bins}");
    hGenInclusiveJetEtaStdBins->GetYaxis()->SetRangeUser(0.00001, 0.05);
    hGenInclusiveJetEtaStdBins->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s std JEC binning", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenInclusiveJetEtaStdBins, "Gen", "p");
    leg->AddEntry(hRecoInclusiveJetEtaStdBins, "Reco", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_etaLab_stdBins_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Plot inclusive reco to gen ratio in lab frame with the standard eta JEC bins
    hInclusiveJetReco2GenEtaStdBins->Draw();
    hInclusiveJetReco2GenEtaStdBins->GetXaxis()->SetTitle("#eta^{Inclusive}_{std bins}");
    hInclusiveJetReco2GenEtaStdBins->GetYaxis()->SetTitle("Reco / Gen");
    hInclusiveJetReco2GenEtaStdBins->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hInclusiveJetReco2GenEtaStdBins->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hInclusiveJetReco2GenEtaStdBins, "Reco / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_etaLab_stdBins_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay leading jet (from dijet) reco, ref and gen eta distributions in lab frame with the standard eta JEC bins
    hGenDijetLeadJetEtaStdBins->Draw();
    hRecoDijetLeadJetEtaStdBins->Draw("same");
    hRefDijetLeadJetEtaStdBins->Draw("same");
    hGenDijetLeadJetEtaStdBins->GetXaxis()->SetTitle("#eta^{Lead}_{std bins}");
    hGenDijetLeadJetEtaStdBins->GetYaxis()->SetTitle("dN/d#eta^{Lead}_{std bins}");
    hGenDijetLeadJetEtaStdBins->GetYaxis()->SetRangeUser(0.00001, 0.05);
    hGenDijetLeadJetEtaStdBins->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s std JEC binning", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetLeadJetEtaStdBins, "Gen", "p");
    leg->AddEntry(hRecoDijetLeadJetEtaStdBins, "Reco", "p");
    leg->AddEntry(hRefDijetLeadJetEtaStdBins, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_dijetLeadJet_etaLab_stdBins_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay ratios of leading jet (from dijet) reco and ref to gen eta distributions in lab frame with the standard eta JEC bins
    hDijetLeadJetReco2GenEtaStdBins->Draw();
    hDijetLeadJetRef2GenEtaStdBins->Draw("same");
    hDijetLeadJetReco2GenEtaStdBins->GetXaxis()->SetTitle("#eta^{Lead}_{std bins}");
    hDijetLeadJetReco2GenEtaStdBins->GetYaxis()->SetTitle("Reco / Gen");
    hDijetLeadJetReco2GenEtaStdBins->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hDijetLeadJetReco2GenEtaStdBins->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hDijetLeadJetReco2GenEtaStdBins, "Reco / Gen", "p");
    leg->AddEntry(hDijetLeadJetRef2GenEtaStdBins, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_dijetLeadJet_etaLab_stdBins_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay subleading jet (from dijet) reco, ref and gen eta distributions in lab frame with the standard eta JEC bins
    hGenDijetSubLeadJetEtaStdBins->Draw();
    hRecoDijetSubLeadJetEtaStdBins->Draw("same");
    hRefDijetSubLeadJetEtaStdBins->Draw("same");
    hGenDijetSubLeadJetEtaStdBins->GetXaxis()->SetTitle("#eta^{SubLead}_{std bins}");
    hGenDijetSubLeadJetEtaStdBins->GetYaxis()->SetTitle("dN/d#eta^{SubLead}_{std bins}");
    hGenDijetSubLeadJetEtaStdBins->GetYaxis()->SetRangeUser(0.00001, 0.05);
    hGenDijetSubLeadJetEtaStdBins->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s std JEC binning", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetSubLeadJetEtaStdBins, "Gen", "p");
    leg->AddEntry(hRecoDijetSubLeadJetEtaStdBins, "Reco", "p");
    leg->AddEntry(hRefDijetSubLeadJetEtaStdBins, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_dijetSubLeadJet_etaLab_stdBins_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of ratios of subleading jet (from dijet) reco and ref to gen eta distributions in lab frame with the standard eta JEC bins
    hDijetSubLeadJetReco2GenEtaStdBins->Draw();
    hDijetSubLeadJetRef2GenEtaStdBins->Draw("same");
    hDijetSubLeadJetReco2GenEtaStdBins->GetXaxis()->SetTitle("#eta^{SubLead}_{std bins}");
    hDijetSubLeadJetReco2GenEtaStdBins->GetYaxis()->SetTitle("Reco / Gen");
    hDijetSubLeadJetReco2GenEtaStdBins->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hDijetSubLeadJetReco2GenEtaStdBins->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hDijetSubLeadJetReco2GenEtaStdBins, "Reco / Gen", "p");
    leg->AddEntry(hDijetSubLeadJetRef2GenEtaStdBins, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_dijetSubLeadJet_etaLab_stdBins_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of lead jet reco and ref to gen ratio in lab frame
    hLeadJetReco2GenEta->Draw();
    hLeadJetRef2GenEta->Draw("same");
    hLeadJetReco2GenEta->GetXaxis()->SetTitle("#eta^{Lead}");
    hLeadJetReco2GenEta->GetYaxis()->SetTitle("Reco / Gen");
    hLeadJetReco2GenEta->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hLeadJetReco2GenEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hLeadJetReco2GenEta, "Reco / Gen", "p");
    leg->AddEntry(hLeadJetRef2GenEta, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_leadJet_etaLab_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of sublead jet reco and ref to gen ratio in lab frame
    hSubLeadJetReco2GenEta->Draw();
    hSubLeadJetRef2GenEta->Draw("same");
    hSubLeadJetReco2GenEta->GetXaxis()->SetTitle("#eta^{SubLead}");
    hSubLeadJetReco2GenEta->GetYaxis()->SetTitle("Reco / Gen");
    hSubLeadJetReco2GenEta->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hSubLeadJetReco2GenEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hSubLeadJetReco2GenEta, "Reco / Gen", "p");
    leg->AddEntry(hSubLeadJetRef2GenEta, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_subLeadJet_etaLab_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of inclusive gen, reco and ref eta distributions in CM frame
    hGenInclusiveJetEtaCM->Draw();
    hRecoInclusiveJetEtaCM->Draw("same");
    hRefInclusiveJetEtaCM->Draw("same");
    meanGenEtaCM = hGenInclusiveJetEtaCM->GetMean();
    meanRecoEtaCM = hRecoInclusiveJetEtaCM->GetMean();
    meanRefEtaCM = hRefInclusiveJetEtaCM->GetMean();
    hGenInclusiveJetEtaCM->GetXaxis()->SetTitle("#eta^{Inclusive}_{CM}");
    hGenInclusiveJetEtaCM->GetYaxis()->SetTitle("dN/d#eta^{Inclusive}_{CM}");
    hGenInclusiveJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenInclusiveJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenInclusiveJetEtaCM, "Gen", "p");
    leg->AddEntry(hRecoInclusiveJetEtaCM, "Reco", "p");
    leg->AddEntry(hRefInclusiveJetEtaCM, "Ref", "p");
    leg->Draw();
    t.SetTextSize(0.03);
    t.DrawLatexNDC(0.7, 0.85, Form("Gen <#eta_{CM}> = %.4f", meanGenEtaCM) );
    t.DrawLatexNDC(0.7, 0.8, Form("Reco <#eta_{CM}> = %.4f", meanRecoEtaCM) );
    t.DrawLatexNDC(0.7, 0.75, Form("Ref <#eta_{CM}> = %.4f", meanRefEtaCM) );
    t.SetTextSize(0.04);
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_etaCM_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of lead jet gen, reco and ref eta distributions in CM frame
    hGenLeadJetEtaCM->Draw();
    hRecoLeadJetEtaCM->Draw("same");
    hRefLeadJetEtaCM->Draw("same");
    meanGenEtaCM = hGenLeadJetEtaCM->GetMean();
    meanRecoEtaCM = hRecoLeadJetEtaCM->GetMean();
    meanRefEtaCM = hRefLeadJetEtaCM->GetMean();
    hGenLeadJetEtaCM->GetXaxis()->SetTitle("#eta^{Lead}_{CM}");
    hGenLeadJetEtaCM->GetYaxis()->SetTitle("dN/d#eta^{Lead}_{CM}");
    hGenLeadJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenLeadJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenLeadJetEtaCM, "Gen", "p");
    leg->AddEntry(hRecoLeadJetEtaCM, "Reco", "p");
    leg->AddEntry(hRefLeadJetEtaCM, "Ref", "p");
    leg->Draw();
    t.SetTextSize(0.03);
    t.DrawLatexNDC(0.7, 0.85, Form("Gen <#eta_{CM}> = %.4f", meanGenEtaCM) );
    t.DrawLatexNDC(0.7, 0.8, Form("Reco <#eta_{CM}> = %.4f", meanRecoEtaCM) );
    t.DrawLatexNDC(0.7, 0.75, Form("Ref <#eta_{CM}> = %.4f", meanRefEtaCM) );
    t.SetTextSize(0.04);
    c->SaveAs(Form("%s/%s_%s_leadJet_etaCM_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), int(ptLow), int(ptHigh) ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of sublead jet gen, reco and ref eta distributions in CM frame
    hGenSubLeadJetEtaCM->Draw();
    hRecoSubLeadJetEtaCM->Draw("same");
    hRefSubLeadJetEtaCM->Draw("same");
    meanGenEtaCM = hGenSubLeadJetEtaCM->GetMean();
    meanRecoEtaCM = hRecoSubLeadJetEtaCM->GetMean();
    meanRefEtaCM = hRefSubLeadJetEtaCM->GetMean();
    hGenSubLeadJetEtaCM->GetXaxis()->SetTitle("#eta^{SubLead}_{CM}");
    hGenSubLeadJetEtaCM->GetYaxis()->SetTitle("dN/d#eta^{SubLead}_{CM}");
    hGenSubLeadJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenSubLeadJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenSubLeadJetEtaCM, "Gen", "p");
    leg->AddEntry(hRecoSubLeadJetEtaCM, "Reco", "p");
    leg->AddEntry(hRefSubLeadJetEtaCM, "Ref", "p");
    leg->Draw();
    t.SetTextSize(0.03);
    t.DrawLatexNDC(0.7, 0.85, Form("Gen <#eta_{CM}> = %.4f", meanGenEtaCM) );
    t.DrawLatexNDC(0.7, 0.8, Form("Reco <#eta_{CM}> = %.4f", meanRecoEtaCM) );
    t.DrawLatexNDC(0.7, 0.75, Form("Ref <#eta_{CM}> = %.4f", meanRefEtaCM) );
    t.SetTextSize(0.04);
    c->SaveAs(Form("%s/%s_%s_subLeadJet_etaCM_comparison_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of inclusive reco and ref to gen ratio in CM frame
    hInclusiveJetReco2GenEtaCM->Draw();
    hInclusiveJetRef2GenEtaCM->Draw("same");
    hInclusiveJetReco2GenEtaCM->GetXaxis()->SetTitle("#eta^{Inclusive}_{CM}");
    hInclusiveJetReco2GenEtaCM->GetYaxis()->SetTitle("Reco / Gen");
    hInclusiveJetReco2GenEtaCM->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hInclusiveJetReco2GenEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hInclusiveJetReco2GenEtaCM, "Reco / Gen", "p");
    leg->AddEntry(hInclusiveJetRef2GenEtaCM, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_etaCM_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of lead jet reco and ref to gen ratio in CM frame
    hLeadJetReco2GenEtaCM->Draw();
    hLeadJetRef2GenEtaCM->Draw("same");
    hLeadJetReco2GenEtaCM->GetXaxis()->SetTitle("#eta^{Lead}_{CM}");
    hLeadJetReco2GenEtaCM->GetYaxis()->SetTitle("Reco / Gen");
    hLeadJetReco2GenEtaCM->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hLeadJetReco2GenEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hLeadJetReco2GenEtaCM, "Reco / Gen", "p");
    leg->AddEntry(hLeadJetRef2GenEtaCM, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_leadJet_etaCM_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of sublead jet reco and ref to gen ratio in CM frame
    hSubLeadJetReco2GenEtaCM->Draw();
    hSubLeadJetRef2GenEtaCM->Draw("same");
    hSubLeadJetReco2GenEtaCM->GetXaxis()->SetTitle("#eta^{SubLead}_{CM}");
    hSubLeadJetReco2GenEtaCM->GetYaxis()->SetTitle("Reco / Gen");
    hSubLeadJetReco2GenEtaCM->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hSubLeadJetReco2GenEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hSubLeadJetReco2GenEtaCM, "Reco / Gen", "p");
    leg->AddEntry(hSubLeadJetRef2GenEtaCM, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_subLeadJet_etaCM_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of gen forward/backward ratio for inclusive, lead and sublead jets in CM frame
    set1DStyle(hGenInclusiveJetForwardBackwardRatioCM, 2);
    set1DStyle(hGenLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hGenSubLeadJetForwardBackwardRatioCM, 1);
    hGenInclusiveJetForwardBackwardRatioCM->Draw();
    hGenLeadJetForwardBackwardRatioCM->Draw("same");
    hGenSubLeadJetForwardBackwardRatioCM->Draw("same");
    hGenInclusiveJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hGenInclusiveJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hGenInclusiveJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenInclusiveJetForwardBackwardRatioCM, "Inclusive", "p");
    leg->AddEntry(hGenLeadJetForwardBackwardRatioCM, "Lead", "p");
    leg->AddEntry(hGenSubLeadJetForwardBackwardRatioCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_singleJet_forwardBackwardRatio_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    
    // Overlay of reco forward/backward ratio for inclusive, lead and sublead jets in CM frame
    set1DStyle(hRecoInclusiveJetForwardBackwardRatioCM, 2);
    set1DStyle(hRecoLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hRecoSubLeadJetForwardBackwardRatioCM, 1);
    hRecoInclusiveJetForwardBackwardRatioCM->Draw();
    hRecoLeadJetForwardBackwardRatioCM->Draw("same");
    hRecoSubLeadJetForwardBackwardRatioCM->Draw("same");
    hRecoInclusiveJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hRecoInclusiveJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hRecoInclusiveJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s reco", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRecoInclusiveJetForwardBackwardRatioCM, "Inclusive", "p");
    leg->AddEntry(hRecoLeadJetForwardBackwardRatioCM, "Lead", "p");
    leg->AddEntry(hRecoSubLeadJetForwardBackwardRatioCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_reco_singleJet_forwardBackwardRatio_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of ref forward/backward ratio for inclusive, lead and sublead jets in CM frame
    set1DStyle(hRefInclusiveJetForwardBackwardRatioCM, 2);
    set1DStyle(hRefLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hRefSubLeadJetForwardBackwardRatioCM, 1);
    hRefInclusiveJetForwardBackwardRatioCM->Draw();
    hRefLeadJetForwardBackwardRatioCM->Draw("same");
    hRefSubLeadJetForwardBackwardRatioCM->Draw("same");
    hRefInclusiveJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hRefInclusiveJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hRefInclusiveJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s ref", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRefInclusiveJetForwardBackwardRatioCM, "Inclusive", "p");
    leg->AddEntry(hRefLeadJetForwardBackwardRatioCM, "Lead", "p");
    leg->AddEntry(hRefSubLeadJetForwardBackwardRatioCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_ref_singleJet_forwardBackwardRatio_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of gen eta distributions for inclusive, lead and sublead jets in CM frame
    set1DStyle(hGenInclusiveJetEtaCM, 2, true);
    set1DStyle(hGenLeadJetEtaCM, 0, true);
    set1DStyle(hGenSubLeadJetEtaCM, 1, true);
    hGenInclusiveJetEtaCM->Draw();
    hGenLeadJetEtaCM->Draw("same");
    hGenSubLeadJetEtaCM->Draw("same");
    hGenInclusiveJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hGenInclusiveJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenInclusiveJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s gen", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenInclusiveJetEtaCM, "Inclusive", "p");
    leg->AddEntry(hGenLeadJetEtaCM, "Lead", "p");
    leg->AddEntry(hGenSubLeadJetEtaCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_singleJet_etaCM_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of reco eta distributions for inclusive, lead and sublead jets in CM frame
    set1DStyle(hRecoInclusiveJetEtaCM, 2, true);
    set1DStyle(hRecoLeadJetEtaCM, 0, true);
    set1DStyle(hRecoSubLeadJetEtaCM, 1, true);
    hRecoInclusiveJetEtaCM->Draw();
    hRecoLeadJetEtaCM->Draw("same");
    hRecoSubLeadJetEtaCM->Draw("same");
    hRecoInclusiveJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hRecoInclusiveJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hRecoInclusiveJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s reco", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRecoInclusiveJetEtaCM, "Inclusive", "p");
    leg->AddEntry(hRecoLeadJetEtaCM, "Lead", "p");
    leg->AddEntry(hRecoSubLeadJetEtaCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_reco_singleJet_etaCM_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }
    
    // Overlay of ref eta distributions for inclusive, lead and sublead jets in CM frame
    set1DStyle(hRefInclusiveJetEtaCM, 2, true);
    set1DStyle(hRefLeadJetEtaCM, 0, true);
    set1DStyle(hRefSubLeadJetEtaCM, 1, true);
    hRefInclusiveJetEtaCM->Draw();
    hRefLeadJetEtaCM->Draw("same");
    hRefSubLeadJetEtaCM->Draw("same");
    hRefInclusiveJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hRefInclusiveJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hRefInclusiveJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s ref", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRefInclusiveJetEtaCM, "Inclusive", "p");
    leg->AddEntry(hRefLeadJetEtaCM, "Lead", "p");
    leg->AddEntry(hRefSubLeadJetEtaCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_ref_singleJet_etaCM_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of gen eta distributions for inclusive, lead and sublead jets in lab frame
    set1DStyle(hGenInclusiveJetEta, 2, true);
    set1DStyle(hGenLeadJetEta, 0, true);
    set1DStyle(hGenSubLeadJetEta, 1, true);
    hGenInclusiveJetEta->Draw();
    hGenLeadJetEta->Draw("same");
    hGenSubLeadJetEta->Draw("same");
    hGenInclusiveJetEta->GetYaxis()->SetTitle("dN/d#eta");
    hGenInclusiveJetEta->GetYaxis()->SetRangeUser(0., 0.1);
    hGenInclusiveJetEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s gen", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenInclusiveJetEta, "Inclusive", "p");
    leg->AddEntry(hGenLeadJetEta, "Lead", "p");
    leg->AddEntry(hGenSubLeadJetEta, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_singleJet_etaLab_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of reco eta distributions for inclusive, lead and sublead jets in lab frame
    set1DStyle(hRecoInclusiveJetEta, 2, true);
    set1DStyle(hRecoLeadJetEta, 0, true);
    set1DStyle(hRecoSubLeadJetEta, 1, true);
    hRecoInclusiveJetEta->Draw();
    hRecoLeadJetEta->Draw("same");
    hRecoSubLeadJetEta->Draw("same");
    hRecoInclusiveJetEta->GetYaxis()->SetTitle("dN/d#eta");
    hRecoInclusiveJetEta->GetYaxis()->SetRangeUser(0., 0.1);
    hRecoInclusiveJetEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s reco", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRecoInclusiveJetEta, "Inclusive", "p");
    leg->AddEntry(hRecoLeadJetEta, "Lead", "p");
    leg->AddEntry(hRecoSubLeadJetEta, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_reco_singleJet_etaLab_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of ref eta distributions for inclusive, lead and sublead jets in lab frame
    set1DStyle(hRefInclusiveJetEta, 2, true);
    set1DStyle(hRefLeadJetEta, 0, true);
    set1DStyle(hRefSubLeadJetEta, 1, true);
    hRefInclusiveJetEta->Draw();
    hRefLeadJetEta->Draw("same");
    hRefSubLeadJetEta->Draw("same");
    hRefInclusiveJetEta->GetYaxis()->SetTitle("dN/d#eta");
    hRefInclusiveJetEta->GetYaxis()->SetRangeUser(0., 0.1);
    hRefInclusiveJetEta->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s ref", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRefInclusiveJetEta, "Inclusive", "p");
    leg->AddEntry(hRefLeadJetEta, "Lead", "p");
    leg->AddEntry(hRefSubLeadJetEta, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_ref_singleJet_etaLab_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of gen lead and sublead full eta distributions from selected dijet events in CM frame
    set1DStyle(hGenDijetLeadJetEtaCM, 0, true);
    set1DStyle(hGenDijetSubLeadJetEtaCM, 1, true);
    hGenDijetLeadJetEtaCM->Draw();
    hGenDijetSubLeadJetEtaCM->Draw("same");
    hGenDijetLeadJetEtaCM->GetXaxis()->SetTitle("#eta_{CM}");
    hGenDijetLeadJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hGenDijetLeadJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenDijetLeadJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s gen (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetLeadJetEtaCM, "Lead", "p");
    leg->AddEntry(hGenDijetSubLeadJetEtaCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_selectedDijet_leadSubleadJet_etaCM_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of gen lead and sublead forward/bacward ratios from selected dijet events
    set1DStyle(hGenDijetLeadJetForwardBackwardRatioCM, 1);
    set1DStyle(hGenDijetSubLeadJetForwardBackwardRatioCM, 0);
    hGenDijetLeadJetForwardBackwardRatioCM->Draw();
    hGenDijetSubLeadJetForwardBackwardRatioCM->Draw("same");
    hGenDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hGenDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hGenDijetLeadJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s gen (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetLeadJetForwardBackwardRatioCM, "Lead", "p");
    leg->AddEntry(hGenDijetSubLeadJetForwardBackwardRatioCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_selectedDijet_leadSubleadJet_forwardBackwardRatio_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of ref lead and sublead forward/bacward ratios from selected dijet events
    set1DStyle(hRefDijetLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hRefDijetSubLeadJetForwardBackwardRatioCM, 1);
    hRefDijetLeadJetForwardBackwardRatioCM->Draw();
    hRefDijetSubLeadJetForwardBackwardRatioCM->Draw("same");
    hRefDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hRefDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hRefDijetLeadJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s ref (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRefDijetLeadJetForwardBackwardRatioCM, "Lead", "p");
    leg->AddEntry(hRefDijetSubLeadJetForwardBackwardRatioCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_ref_selectedDijet_leadSubleadJet_forwardBackwardRatio_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of reco lead and sublead full eta distributions from selected dijet events in CM frame
    set1DStyle(hRecoDijetLeadJetEtaCM, 0);
    set1DStyle(hRecoDijetSubLeadJetEtaCM, 1);
    hRecoDijetLeadJetEtaCM->Draw();
    hRecoDijetSubLeadJetEtaCM->Draw("same");
    hRecoDijetLeadJetEtaCM->GetXaxis()->SetTitle("#eta_{CM}");
    hRecoDijetLeadJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hRecoDijetLeadJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hRecoDijetLeadJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s reco (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRecoDijetLeadJetEtaCM, "Lead", "p");
    leg->AddEntry(hRecoDijetSubLeadJetEtaCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_reco_selectedDijet_leadSubleadJet_etaCM_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of reco lead and sublead forward/bacward ratios from selected dijet events
    set1DStyle(hRecoDijetLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hRecoDijetSubLeadJetForwardBackwardRatioCM, 1);
    hRecoDijetLeadJetForwardBackwardRatioCM->Draw();
    hRecoDijetSubLeadJetForwardBackwardRatioCM->Draw("same");
    hRecoDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hRecoDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hRecoDijetLeadJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s reco (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.16, 0.65, 0.34, 0.75);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hRecoDijetLeadJetForwardBackwardRatioCM, "Lead", "p");
    leg->AddEntry(hRecoDijetSubLeadJetForwardBackwardRatioCM, "SubLead", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_reco_selectedDijet_leadSubleadJet_forwardBackwardRatio_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of lead jet reco and ref to gen ratio from selected dijet events in CM frame
    set1DStyle(hDijetLeadJetReco2GenEtaCM, 0);
    set1DStyle(hDijetLeadJetRef2GenEtaCM, 1);
    hDijetLeadJetReco2GenEtaCM->Draw();
    hDijetLeadJetRef2GenEtaCM->Draw("same");
    hDijetLeadJetReco2GenEtaCM->GetXaxis()->SetTitle("#eta^{Lead}_{CM}");
    hDijetLeadJetReco2GenEtaCM->GetYaxis()->SetTitle("Ratio to Gen");
    hDijetLeadJetReco2GenEtaCM->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hDijetLeadJetReco2GenEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s (selected dijet events)", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hDijetLeadJetReco2GenEtaCM, "Reco / Gen", "p");
    leg->AddEntry(hDijetLeadJetRef2GenEtaCM, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_selectedDijet_leadJet_etaCM_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of sublead jet reco and ref to gen ratio from selected dijet events in CM frame
    set1DStyle(hDijetSubLeadJetReco2GenEtaCM, 0);
    set1DStyle(hDijetSubLeadJetRef2GenEtaCM, 1);
    hDijetSubLeadJetReco2GenEtaCM->Draw();
    hDijetSubLeadJetRef2GenEtaCM->Draw("same");
    hDijetSubLeadJetReco2GenEtaCM->GetXaxis()->SetTitle("#eta^{SubLead}_{CM}");
    hDijetSubLeadJetReco2GenEtaCM->GetYaxis()->SetTitle("Ratio to Gen");
    hDijetSubLeadJetReco2GenEtaCM->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
    hDijetSubLeadJetReco2GenEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s (selected dijet events)", mcType.Data()) );
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hDijetSubLeadJetReco2GenEtaCM, "Reco / Gen", "p");
    leg->AddEntry(hDijetSubLeadJetRef2GenEtaCM, "Ref / Gen", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_selectedDijet_subleadJet_etaCM_ratio2gen_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of lead jet reco, ref and gen full eta distributions from selected dijet events in CM frame
    set1DStyle(hGenDijetLeadJetEtaCM, 1);
    set1DStyle(hRecoDijetLeadJetEtaCM, 0);
    set1DStyle(hRefDijetLeadJetEtaCM, 2);
    hGenDijetLeadJetEtaCM->Draw();
    hRecoDijetLeadJetEtaCM->Draw("same");
    hRefDijetLeadJetEtaCM->Draw("same");
    meanGenEtaCM = hGenDijetLeadJetEtaCM->GetMean();
    meanRecoEtaCM = hRecoDijetLeadJetEtaCM->GetMean();
    meanRefEtaCM = hRefDijetLeadJetEtaCM->GetMean();
    hGenDijetLeadJetEtaCM->GetXaxis()->SetTitle("#eta^{Lead}_{CM}");
    hGenDijetLeadJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hGenDijetLeadJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenDijetLeadJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetLeadJetEtaCM, "Gen", "p");
    leg->AddEntry(hRecoDijetLeadJetEtaCM, "Reco", "p");
    leg->AddEntry(hRefDijetLeadJetEtaCM, "Ref", "p");
    leg->Draw();
    t.SetTextSize(0.03);
    t.DrawLatexNDC(0.7, 0.85, Form("Gen <#eta_{CM}> = %.4f", meanGenEtaCM) );
    t.DrawLatexNDC(0.7, 0.8, Form("Reco <#eta_{CM}> = %.4f", meanRecoEtaCM) );
    t.DrawLatexNDC(0.7, 0.75, Form("Ref <#eta_{CM}> = %.4f", meanRefEtaCM) );
    t.SetTextSize(0.04);
    c->SaveAs(Form("%s/%s_%s_selectedDijet_leadJet_etaCM_genRecoRef_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of sublead jet reco, ref and gen full eta distributions from selected dijet events in CM frame
    set1DStyle(hGenDijetSubLeadJetEtaCM, 1);
    set1DStyle(hRecoDijetSubLeadJetEtaCM, 0);
    set1DStyle(hRefDijetSubLeadJetEtaCM, 2);
    hGenDijetSubLeadJetEtaCM->Draw();
    hRecoDijetSubLeadJetEtaCM->Draw("same");
    hRefDijetSubLeadJetEtaCM->Draw("same");
    meanGenEtaCM = hGenDijetSubLeadJetEtaCM->GetMean();
    meanRecoEtaCM = hRecoDijetSubLeadJetEtaCM->GetMean();
    meanRefEtaCM = hRefDijetSubLeadJetEtaCM->GetMean();
    hGenDijetSubLeadJetEtaCM->GetXaxis()->SetTitle("#eta^{SubLead}_{CM}");
    hGenDijetSubLeadJetEtaCM->GetYaxis()->SetTitle("dN/d#eta_{CM}");
    hGenDijetSubLeadJetEtaCM->GetYaxis()->SetRangeUser(0., 0.1);
    hGenDijetSubLeadJetEtaCM->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetSubLeadJetEtaCM, "Gen", "p");
    leg->AddEntry(hRecoDijetSubLeadJetEtaCM, "Reco", "p");
    leg->AddEntry(hRefDijetSubLeadJetEtaCM, "Ref", "p");
    leg->Draw();
    t.SetTextSize(0.03);
    t.DrawLatexNDC(0.7, 0.85, Form("Gen <#eta_{CM}> = %.4f", meanGenEtaCM) );
    t.DrawLatexNDC(0.7, 0.8, Form("Reco <#eta_{CM}> = %.4f", meanRecoEtaCM) );
    t.DrawLatexNDC(0.7, 0.75, Form("Ref <#eta_{CM}> = %.4f", meanRefEtaCM) );
    t.SetTextSize(0.04);
    c->SaveAs(Form("%s/%s_%s_selectedDijet_subleadJet_etaCM_genRecoRef_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of reco, ref and gen lead forward/bacward ratios from selected dijet events
    set1DStyle(hGenDijetLeadJetForwardBackwardRatioCM, 1);
    set1DStyle(hRecoDijetLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hRefDijetLeadJetForwardBackwardRatioCM, 2);
    hGenDijetLeadJetForwardBackwardRatioCM->Draw();
    hRecoDijetLeadJetForwardBackwardRatioCM->Draw("same");
    hRefDijetLeadJetForwardBackwardRatioCM->Draw("same");
    hGenDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hGenDijetLeadJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hGenDijetLeadJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s leading jet (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetLeadJetForwardBackwardRatioCM, "Gen", "p");
    leg->AddEntry(hRecoDijetLeadJetForwardBackwardRatioCM, "Reco", "p");
    leg->AddEntry(hRefDijetLeadJetForwardBackwardRatioCM, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_selectedDijet_leadJet_forwardBackwardRatio_genRecoRef_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of reco, ref and gen sublead forward/bacward ratios from selected dijet events
    set1DStyle(hGenDijetSubLeadJetForwardBackwardRatioCM, 1);
    set1DStyle(hRecoDijetSubLeadJetForwardBackwardRatioCM, 0);
    set1DStyle(hRefDijetSubLeadJetForwardBackwardRatioCM, 2);
    hGenDijetSubLeadJetForwardBackwardRatioCM->Draw();
    hRecoDijetSubLeadJetForwardBackwardRatioCM->Draw("same");
    hRefDijetSubLeadJetForwardBackwardRatioCM->Draw("same");
    hGenDijetSubLeadJetForwardBackwardRatioCM->GetYaxis()->SetTitle("Forward / Backward");
    hGenDijetSubLeadJetForwardBackwardRatioCM->GetYaxis()->SetRangeUser(0.85, 1.15);
    hGenDijetSubLeadJetForwardBackwardRatioCM->GetXaxis()->SetRangeUser(0., etaCMMax);
    t.DrawLatexNDC(0.16, 0.85, Form("%s subleading jet (selected dijet events)", mcType.Data()));
    t.DrawLatexNDC(0.16, 0.8, Form("%d < p_{T}^{jet} < %d GeV", ptLowInt, ptHighInt) );
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg = new TLegend(0.45, 0.2, 0.65, 0.3);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenDijetSubLeadJetForwardBackwardRatioCM, "Gen", "p");
    leg->AddEntry(hRecoDijetSubLeadJetForwardBackwardRatioCM, "Reco", "p");
    leg->AddEntry(hRefDijetSubLeadJetForwardBackwardRatioCM, "Ref", "p");
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_selectedDijet_subleadJet_forwardBackwardRatio_genRecoRef_pt_%d_%d.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data(), ptLowInt, ptHighInt ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay various gen eta cm distributions for different pT bins
    static constexpr std::array<const char*, 7> p8Colors {"kBlue", "kRed", "kMagenta", "kP8Orange", "kP8Green", "kP8Azure", "kBlack"};
    static constexpr std::array<int, 7> markerStyles {20, 21, 22, 23, 24, 25, 26};
    leg = new TLegend(0.16, 0.65, 0.34, 0.83);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    for ( int iPtBin = 0; iPtBin < nPtBins; ++iPtBin ) {
        hGenInclusiveJetEtaCMPtBins[iPtBin]->SetLineColor( p8Colors[iPtBin] );
        hGenInclusiveJetEtaCMPtBins[iPtBin]->SetMarkerColor( p8Colors[iPtBin] );
        hGenInclusiveJetEtaCMPtBins[iPtBin]->SetMarkerStyle( markerStyles[iPtBin] );
        meanGenEtaCM = hGenInclusiveJetEtaCMPtBins[iPtBin]->GetMean();

        hGenInclusiveJetEtaCMPtBins[iPtBin]->Draw( (iPtBin == 0) ? "" : "same" );
        if (iPtBin == 0) {
            hGenInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->SetTitle("#eta_{CM}");
            hGenInclusiveJetEtaCMPtBins[iPtBin]->GetYaxis()->SetTitle("dN/d#eta_{CM}");
            hGenInclusiveJetEtaCMPtBins[iPtBin]->GetYaxis()->SetRangeUser(0., 0.1);
            hGenInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->SetRangeUser(-etaCMMax-1., etaCMMax+1.);
        }
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.16, 0.85, Form("%s gen (inclusive jets)", mcType.Data()));
        leg->AddEntry(hGenInclusiveJetEtaCMPtBins[iPtBin], Form("%3.f < p_{T}^{jet} < %3.f GeV; <#eta_{CM}> = %.4f", 
                                                                ptVals[iPtBin], ptVals[iPtBin+1], meanGenEtaCM), "p");
    }
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_inclusiveJet_etaCM_ptBin_comparison.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data() ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay various reco eta cm distributions for different pT bins
    leg = new TLegend(0.16, 0.65, 0.34, 0.83);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    for ( int iPtBin = 0; iPtBin < nPtBins; ++iPtBin ) {
        hRecoInclusiveJetEtaCMPtBins[iPtBin]->SetLineColor( p8Colors[iPtBin] );
        hRecoInclusiveJetEtaCMPtBins[iPtBin]->SetMarkerColor( p8Colors[iPtBin] );
        hRecoInclusiveJetEtaCMPtBins[iPtBin]->SetMarkerStyle( markerStyles[iPtBin] );
        meanRecoEtaCM = hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetMean();
        hRecoInclusiveJetEtaCMPtBins[iPtBin]->Draw( (iPtBin == 0) ? "" : "same" );
        if (iPtBin == 0) {
            hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->SetTitle("#eta_{CM}");
            hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetYaxis()->SetTitle("dN/d#eta_{CM}");
            hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetYaxis()->SetRangeUser(0., 0.1);
            hRecoInclusiveJetEtaCMPtBins[iPtBin]->GetXaxis()->SetRangeUser(-etaCMMax, etaCMMax);
        }
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.16, 0.85, Form("%s reco (inclusive jets)", mcType.Data()));
        leg->AddEntry(hRecoInclusiveJetEtaCMPtBins[iPtBin], Form("%3.f < p_{T}^{jet} < %3.f GeV; <#eta_{CM}> = %.4f", 
                                                                ptVals[iPtBin], ptVals[iPtBin+1], meanRecoEtaCM), "p");
    }
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_reco_inclusiveJet_etaCM_ptBin_comparison.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data() ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay inclusive jet reco to gen ratios in standard JEC bins for different pT bins
    leg = new TLegend(0.16, 0.18, 0.34, 0.4);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    for ( int iPtBin = 0; iPtBin < nPtBins; ++iPtBin ) {
        hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->SetLineColor( p8Colors[iPtBin] );
        hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->SetMarkerColor( p8Colors[iPtBin] );
        hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->SetMarkerStyle( markerStyles[iPtBin] );
        hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->Draw( (iPtBin == 0) ? "" : "same" );
        if (iPtBin == 0) {
            hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->GetXaxis()->SetRangeUser(-etaLabMax, etaLabMax);
            hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->GetXaxis()->SetTitle("#eta");
            hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->GetYaxis()->SetTitle("Reco / Gen");
            hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin]->GetYaxis()->SetRangeUser(yAxisRatio2Gen[0], yAxisRatio2Gen[1]);
        }
        t.DrawLatexNDC(0.16, 0.85, Form("%s inclusive jets (standard JEC bins)", mcType.Data()));
        leg->AddEntry(hInclusiveJetReco2GenEtaStdBinsPtBins[iPtBin], Form("%3.f < p_{T}^{jet} < %3.f GeV", ptVals[iPtBin], ptVals[iPtBin+1]), "p");
    }
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_inclusiveJet_eta_reco2gen_stdJECbins_ptBin_comparison.pdf", 
                    date.Data(), collSystemStr.Data(), mcTypeLower.Data() ) );
    if (leg) { delete leg; leg = nullptr; }

    if ( c ) { delete c; c = nullptr; }
    if ( leg ) { delete leg; leg = nullptr; }
}


//________________
void plotDiJetClosures(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    // pT bin selection
    double deltaPt = 0.001;
    double ptLow = 60. + deltaPt;
    double ptHigh = 119.999;
    int ptLowInt = std::round(ptLow);
    int ptHighInt = std::round(ptHigh);
    double etaLabMax = 2.4;
    double etaCMMax = 2.0;
    double yAxisRatio2Gen[2] = {0.85, 1.15};

    TString mcType = "PYTHIA"; // "pythia" or "embedding"
    TString mcTypeLower = mcType;
    mcTypeLower.ToLower();

    TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_jetId_eta19.root") );
    // TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_jetId_eta19.root") );

    // TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/Pbgoing/oPythia_pPb8160_Pbgoing_30.root") );
    // TFile *fPythia = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_jetId_eta19_ptHat_gt_30.root") );
    if ( !fMC || fMC->IsZombie() ) {
        std::cerr << "Error: Could not open pythia file." << std::endl;
        return;
    }
}

//________________
void plotOverweightProtection(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TFile *fMC = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_jerDef_ak4_jetId_eta19.root") );

    if ( !fMC || fMC->IsZombie() ) {
        std::cerr << "Error: Could not open pythia file." << std::endl;
        return;
    }

    TH2D *hRecoDijetPtAveOverPtHatVsPtHatWeighted = dynamic_cast<TH2D*>( fMC->Get("hRecoDijetPtAveOverPtHatVsPtHatWeighted") );
    hRecoDijetPtAveOverPtHatVsPtHatWeighted->SetDirectory(0);

    fMC->Close();
    
    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();
    hRecoDijetPtAveOverPtHatVsPtHatWeighted->Draw("colz");
}

//________________
void plotGenEtaShiftCheck(int collisionSystem, double collisionEnergy, TString date) {
    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    bool isPythia = false; // Set to false for embedding
    TString generatorType = isPythia ? "pythia" : "embedding";
    TString direction = "Pbgoing"; // "pgoing", "Pbgoing", "combined"
    // TString inputFileName = Form("eta_shift/%s_%s.root", generatorType.Data(), direction.Data() );
    TString inputFileName = Form("eta_shift/%s_%s_ptHat30.root", generatorType.Data(), direction.Data() );
    const int nEtaShifts = 13;
    static constexpr std::array<float, nEtaShifts> etaShift{0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495 };
    static constexpr std::array<int, nEtaShifts> markerStyles{20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
    static constexpr std::array<const char*, nEtaShifts> p8Colors {"kBlue", "kRed", "kBlack", "kMagenta", "kP8Orange", "kP8Green", "kP8Azure", "kPink", "kCyan", "kTeal", "kGray", "kSpring", "kViolet"};
    int lineWidth = 3;
    double markerSize = 1.3;


    TFile *f = TFile::Open( inputFileName.Data() );
    if ( !f || f->IsZombie() ) {
        std::cerr << Form("Error: Could not open %s", inputFileName.Data()) << std::endl;
        return;
    }

    int rebin = 4;

    // Histograms for different eta shifts
    TH1D *hGenEtaCM[nEtaShifts]{nullptr};
    TH1D *hGenEtaCMForward[nEtaShifts]{nullptr};
    TH1D *hGenEtaCMBackward[nEtaShifts]{nullptr};
    TH1D *hGenEtaCMForwardBackwardRatio[nEtaShifts]{nullptr};

    for ( int iEta = 0; iEta < nEtaShifts; ++iEta ) {
        hGenEtaCM[iEta] = dynamic_cast<TH1D*>( f->Get( Form("hGenJetEtaCMShifted_%d", iEta) ) );
        if ( !hGenEtaCM[iEta] ) {
            std::cerr << Form("Error: Could not retrieve hGenJetEtaCMShifted_%d from file.", iEta) << std::endl;
            return;
        }
        hGenEtaCM[iEta]->SetDirectory(0);
        set1DStyle(hGenEtaCM[iEta], 0);
        hGenEtaCM[iEta]->SetName( Form("hGenJetEtaCM_%d", iEta) );
        hGenEtaCM[iEta]->Rebin(rebin);
        hGenEtaCM[iEta]->SetMarkerStyle(markerStyles[iEta]);
        hGenEtaCM[iEta]->SetMarkerColor(p8Colors[iEta]);
        hGenEtaCM[iEta]->SetLineColor(p8Colors[iEta]);
        hGenEtaCM[iEta]->SetLineWidth(lineWidth);
        hGenEtaCM[iEta]->SetMarkerSize(markerSize);
        hGenEtaCM[iEta]->Scale(1./hGenEtaCM[iEta]->Integral(hGenEtaCM[iEta]->GetXaxis()->FindBin(-1.899), hGenEtaCM[iEta]->GetXaxis()->FindBin(1.8999)) );

        hGenEtaCMForward[iEta] = new TH1D(Form("hGenEtaCMForward_%d", iEta), 
                                          Form("hGenEtaCMForward_%d;#eta_{CM};dN/d#eta_{CM}", iEta), 
                                          hGenEtaCM[iEta]->GetNbinsX()/2, hGenEtaCM[iEta]->GetXaxis()->GetXmin(), hGenEtaCM[iEta]->GetXaxis()->GetXmax() );
        hGenEtaCMForward[iEta]->SetDirectory(0);
        hGenEtaCMForward[iEta]->Sumw2();
        set1DStyle(hGenEtaCMForward[iEta], 0);
        hGenEtaCMForward[iEta]->SetMarkerStyle(markerStyles[iEta]);
        hGenEtaCMForward[iEta]->SetMarkerColor(p8Colors[iEta]);
        hGenEtaCMForward[iEta]->SetLineColor(p8Colors[iEta]);
        hGenEtaCMForward[iEta]->SetLineWidth(lineWidth);
        hGenEtaCMForward[iEta]->SetMarkerSize(markerSize);

        hGenEtaCMBackward[iEta] = new TH1D(Form("hGenEtaCMBackward_%d", iEta), 
                                           Form("hGenEtaCMBackward_%d;#eta_{CM};dN/d#eta_{CM}", iEta), 
                                           hGenEtaCM[iEta]->GetNbinsX()/2, hGenEtaCM[iEta]->GetXaxis()->GetXmin(), hGenEtaCM[iEta]->GetXaxis()->GetXmax() );
        hGenEtaCMBackward[iEta]->SetDirectory(0);
        hGenEtaCMBackward[iEta]->Sumw2();
        set1DStyle(hGenEtaCMBackward[iEta], 0);
        hGenEtaCMBackward[iEta]->SetMarkerStyle(markerStyles[iEta]);
        hGenEtaCMBackward[iEta]->SetMarkerColor(p8Colors[iEta]);
        hGenEtaCMBackward[iEta]->SetLineColor(p8Colors[iEta]);
        hGenEtaCMBackward[iEta]->SetLineWidth(lineWidth);
        hGenEtaCMBackward[iEta]->SetMarkerSize(markerSize);

        hGenEtaCMForwardBackwardRatio[iEta] = new TH1D(Form("hGenEtaCMForwardBackwardRatio_%d", iEta), 
                                                        Form("hGenEtaCMForwardBackwardRatio_%d;#eta_{CM};Forward / Backward", iEta), 
                                                        hGenEtaCM[iEta]->GetNbinsX()/2, 0., hGenEtaCM[iEta]->GetXaxis()->GetXmax() );
        hGenEtaCMForwardBackwardRatio[iEta]->SetDirectory(0);
        hGenEtaCMForwardBackwardRatio[iEta]->Sumw2();
        set1DStyle(hGenEtaCMForwardBackwardRatio[iEta], 0);
        hGenEtaCMForwardBackwardRatio[iEta]->SetMarkerStyle(markerStyles[iEta]);
        hGenEtaCMForwardBackwardRatio[iEta]->SetMarkerColor(p8Colors[iEta]);
        hGenEtaCMForwardBackwardRatio[iEta]->SetLineColor(p8Colors[iEta]);
        hGenEtaCMForwardBackwardRatio[iEta]->SetLineWidth(lineWidth);
        hGenEtaCMForwardBackwardRatio[iEta]->SetMarkerSize(markerSize);

        recalculateFBRatioFromFullDistribution(hGenEtaCMForwardBackwardRatio[iEta], hGenEtaCMForward[iEta], 
                                               hGenEtaCMBackward[iEta], hGenEtaCM[iEta]);
    } // for ( int iEta = 0; iEta < nEtaShifts; ++iEta )

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    setPadStyle();
    gPad->SetGrid();

    TLegend *leg{nullptr};
    TLatex t;
    t.SetTextSize(0.04);
    t.SetTextFont(42);   

    // Overlay of gen eta cm distributions for different eta shifts
    leg = new TLegend(0.16, 0.55, 0.34, 0.83);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    for ( int iEta = 0; iEta < nEtaShifts; ++iEta ) {
        hGenEtaCM[iEta]->Draw( (iEta == 0) ? "" : "same" );
        if (iEta == 0) {
            hGenEtaCM[iEta]->GetXaxis()->SetTitle("#eta_{CM}");
            hGenEtaCM[iEta]->GetYaxis()->SetTitle("dN/d#eta_{CM}");
            hGenEtaCM[iEta]->GetYaxis()->SetRangeUser(0., 0.06);
            hGenEtaCM[iEta]->GetXaxis()->SetRangeUser(-4.8, 4.8);
        }
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.16, 0.85, Form("%s gen (inclusive jets)", generatorType.Data()));
        leg->AddEntry(hGenEtaCM[iEta], Form("#Delta#eta = %.3f; #LT#eta_{CM}#GT = %.4f; #gamma_{1} = %.4f", etaShift[iEta], hGenEtaCM[iEta]->GetMean(), hGenEtaCM[iEta]->GetSkewness()), "p");
    }
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_inclusiveJet_etaCM_etaShift_comparison.pdf", 
                    date.Data(), generatorType.Data(), direction.Data() ) );
    if (leg) { delete leg; leg = nullptr; }

    // Overlay of gen forward/backward cm distributions for different eta shifts
    leg = new TLegend(0.16, 0.55, 0.34, 0.83);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    TF1 *fLine[nEtaShifts]{nullptr};

    for ( int iEta = 0; iEta < nEtaShifts; ++iEta ) {
        hGenEtaCMForwardBackwardRatio[iEta]->Draw( (iEta == 0) ? "" : "same" );
        fLine[iEta] = new TF1(Form("fLine_%d", iEta), "pol0", 0., 1.9);
        fLine[iEta]->SetLineColor(p8Colors[iEta]);
        fLine[iEta]->SetLineWidth(lineWidth);
        fLine[iEta]->SetLineStyle(2);
        hGenEtaCMForwardBackwardRatio[iEta]->Fit( fLine[iEta], "MRQE", "", 0., 1.9 );
        if (iEta == 0) {
            hGenEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetTitle("#eta_{CM}");
            hGenEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetTitle("Forward / Backward");
            hGenEtaCMForwardBackwardRatio[iEta]->GetYaxis()->SetRangeUser(0.85, 1.25);
            hGenEtaCMForwardBackwardRatio[iEta]->GetXaxis()->SetRangeUser(0., 2.0);
        }
        plotCMSHeader(collisionSystem, collisionEnergy);
        t.DrawLatexNDC(0.16, 0.85, Form("%s gen (inclusive jets)", generatorType.Data()));
        leg->AddEntry(hGenEtaCMForwardBackwardRatio[iEta], Form("#Delta#eta = %.3f; Offset = %.3f #pm %.3f", etaShift[iEta], fLine[iEta]->GetParameter(0), fLine[iEta]->GetParError(0)), "p");
    }
    plotCMSHeader(collisionSystem, collisionEnergy);
    leg->Draw();
    c->SaveAs(Form("%s/%s_%s_gen_inclusiveJet_forwardBackwardRatio_etaCM_etaShift_comparison.pdf", 
                    date.Data(), generatorType.Data(), direction.Data() ) );
    if (leg) { delete leg; leg = nullptr; }
    for ( int iEta = 0; iEta < nEtaShifts; ++iEta ) {
        if (fLine[iEta]) { delete fLine[iEta]; fLine[iEta] = nullptr; }
    }
    if ( c ) { delete c; c = nullptr; }
}

//________________
void plotMcClosures() {

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Date extraction
    TDatime dt; 
    TString date { Form( "%d",dt.GetDate() ) };

    if ( !directoryExists( date.Data() ) ) {
        createDirectory( date.Data() );
    } 

    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    //
    // pPb8160
    //

    int collisionSystem = 1;         // 0 - pp, 1 - pPb, 2 - pPb5020, 3 - pPb8160
    double collisionEnergy = 8.16;   // 8.16 TeV
    int direction = 2;               // 0-p-going, 1-Pb-going, 2 - combined
    TString directionStr = (direction == 0) ? "pgoing" : ((direction == 1) ? "Pbgoing" : "");
    int dataTrigger = 3;               // 0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100
    TString dataStr = (dataTrigger == 0) ? "MB" : ((dataTrigger == 1) ? "Jet60" : ((dataTrigger == 2) ? "Jet80" : ((dataTrigger == 3) ? "Jet100" : "unknownData")));
    TString dataDirectionStr = (direction == 0) ? "Pbgoing" : ((direction == 1) ? "pgoing" : "");
    int jetType = 0; // 0 - inclusive, 1 - lead, 2 - sublead
    int matchType = 0; // 0 - inclusive, 1 - matched, 2 - unmatched
    int cutType = 2; // 0 - no cuts, 1 - trkMax, 2 - jetId
    TString cutTypeStr = (cutType == 0) ? "noCut_" : ((cutType == 1) ? "trkMax_" : ((cutType == 2) ? "jetId_" : "unknownCut"));
    int etaRange = 19;

    //
    // Embedding
    //
    TFile *pPb8160EmbedFile = nullptr;
    if ( direction < 2 ) {
        pPb8160EmbedFile = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_%seta%d.root", 
                                        directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160EmbedFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/embedding/%s/oEmbedding_%s_def_ak4_%seta%d.root", 
                               directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    else {
        pPb8160EmbedFile = TFile::Open( Form("~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_%seta%d.root", 
                                        cutTypeStr.Data(), etaRange) );
        if ( !pPb8160EmbedFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/embedding/oEmbedding_pPb8160_def_ak4_%seta%d.root", 
                               cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    // TFile *pPb8160EmbedFile = TFile::Open( Form("/Users/%s/work/cms/soft/jetAnalysis/macro/pPb8160_Pbgoing_noTrkMax/oEmbedding_pPb8160_Pbgoing_80.root", uname.Data()) );

    //
    // Pythia
    //
    TFile *pPb8160PythiaFile = nullptr;
    if ( direction < 2 ) {
        pPb8160PythiaFile = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_%seta%d.root", 
                                         directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160PythiaFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/pythia/%s/oPythia_%s_def_ak4_%seta%d.root", 
                              directionStr.Data(), directionStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    else {
        pPb8160PythiaFile = TFile::Open( Form("~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_%seta%d.root", 
                                         cutTypeStr.Data(), etaRange) );
        // if ( !pPb8160PythiaFile ) {
        //     std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/pythia/oPythia_pPb8160_def_ak4_%seta%d.root", 
        //                       cutTypeStr.Data(), etaRange) << std::endl;
        //     return;
        // }
    }

    //
    // Data
    //
    TFile *pPb8160DataFile = nullptr;
    if ( direction < 2 ) {
        pPb8160DataFile = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/%s/%s_%s_ak4_%seta%d.root", 
            dataDirectionStr.Data(), dataStr.Data(), dataDirectionStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160DataFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/exp/%s/%s_%s_ak4_%seta%d.root", 
                dataDirectionStr.Data(), dataStr.Data(), dataDirectionStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }
    else {
        pPb8160DataFile = TFile::Open( Form("~/cernbox/ana/pPb8160/exp/%s_pPb8160_ak4_%seta%d.root", 
            dataStr.Data(), cutTypeStr.Data(), etaRange) );
        if ( !pPb8160DataFile ) {
            std::cerr << Form("File not found: ~/cernbox/ana/pPb8160/exp/%s_pPb8160_ak4_%seta%d.root", 
                dataStr.Data(), cutTypeStr.Data(), etaRange) << std::endl;
            return;
        }
    }


    //
    // Comparison of dijet reco and ref to gen distributions
    //
    // dijetClosures( pPb8160EmbedFile, collisionSystem, collisionEnergy, date );

    //
    // Plot comparison of inclusive jet eta distributions to check/validate the JEC
    //
    // data2mcInclusiveJetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of dijet reco and ref to gen distributions
    //
    // data2mcDijetComparison(pPb8160DataFile, pPb8160EmbedFile, collisionSystem, collisionEnergy, date);
    // data2mcDijetComparison(pPb8160DataFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of inclusive jets and dijets for embedding and PYTHIA
    //
    // pythia2embeddingSingleJet(pPb8160EmbedFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);

    //
    // Plot comparison of dijet distributions from embedding and Pythia
    //
    // pythia2embeddingDijetComparison(pPb8160EmbedFile, pPb8160PythiaFile, collisionSystem, collisionEnergy, date);


    //
    // Plot distributions for the reco dijets
    //
    // recoDijetDistributions(pPb8160EmbedFile, collisionSystem, collisionEnergy, date, true);

    // plotForwardBackwardRatios(pPb8160DataFile, collisionSystem, collisionEnergy, date);

    // plotXjComparison(collisionSystem, collisionEnergy, date);

    // Plot forward/backward ratios for dijets in different |eta_CM|<X cuts
    // plotDijetFBRatiosForEta(collisionSystem, collisionEnergy, date);

    // Plot eta distributions for inclusive jets, leading jets, subleading jets, and dijets in CM frame for different run periods
    // plotEtaDistributionsForRunId(collisionSystem, collisionEnergy, date);

    // Plot comparison of forward/backward ratios and full eta distributions between nPDF calculations and data
    //plotNpdfToDataComparison(collisionSystem, collisionEnergy, date);

    // Plot comparison of forward/backward ratios and full eta distributions between PDF calculations and data
    // plotPdfAndPythiaToDataComparison(collisionSystem, collisionEnergy, date);

    // Plot comparison of forward/backward ratios and full eta distributions between nPDF and PDF calculations
    // plotNpdfToPdfComparison(collisionSystem, collisionEnergy, date);

    // Plot single jet closures (gen, reco, ref, forward/backward ratios, eta distributions) 
    // for inclusive, lead and sublead jets in CM frame
    // plotSingleJetClosures(collisionSystem, collisionEnergy, date);

    // plotOverweightProtection(collisionSystem, collisionEnergy, date);

    // Plot check of gen eta shift in pPb collisions
    plotGenEtaShiftCheck(collisionSystem, collisionEnergy, date);
}
