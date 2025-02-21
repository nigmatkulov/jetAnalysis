// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TRatioPlot.h"
#include "TGraphErrors.h"
#include "TPad.h"

// C++ headers
#include <iostream>
#include <vector>
#include <cmath>
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
std::vector<int> findFivePercBins(TH1D *h) {
    // Find 5% bins
    std::vector<int> bins;
    double cumulativeIntegral = 0.;
    double totalIntegral = h->Integral();
    double fivePercMark = 0.05 * totalIntegral;
    double previousBinIntegral = 0.;
    int bin = 1;
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        cumulativeIntegral += h->GetBinContent(i);
        if (cumulativeIntegral > fivePercMark) {
            bin = ( abs(fivePercMark - previousBinIntegral) < abs(fivePercMark - cumulativeIntegral) ) ? i - 1 : i;
            bins.push_back(bin);
            cumulativeIntegral = 0.;
        }
        previousBinIntegral = cumulativeIntegral;
    } // for (int i = 1; i <= h->GetNbinsX(); i++)

    for (int j = 0; j < bins.size(); j++) {
        std::cout << "5% bin: " << bins[j] << std::endl;
    }

    return bins;
}

//________________
void calculateMeanAndErrorInRange(TH1D* h, int x1, int x2, double &mean, double &meanError) {
    double sum = 0.0;
    double weightedSum = 0.0;
    double sumOfWeights = 0.0;
    double sumOfWeightsSquared = 0.0;

    // std::cout << Form("hName: %s, x1 = %d, x2 = %d", h->GetName(), x1, x2) << std::endl;

    // In case entries exist in the specified range
    if (h->Integral(x1, x2) > 0) {
        for (int i = x1; i < x2; ++i) {
            double binContent = h->GetBinContent(i);
            double binCenter = h->GetBinCenter(i);
            double binError = h->GetBinError(i);

            sum += binContent;
            weightedSum += binContent * binCenter;
            sumOfWeights += binContent;
            sumOfWeightsSquared += binError * binError;
        }

        mean = weightedSum / sum;
        meanError = std::sqrt(sumOfWeightsSquared) / sum;
    }
    else { // In case of empty bins just assign the mean of the interval
        mean = ( h->GetXaxis()->GetBinLowEdge(x1) + h->GetXaxis()->GetBinUpEdge(x2) ) / 2;
        meanError = 0.000001;
    }

    // std::cout << Form("Mean = %.4f, Mean error = %.4f", mean, meanError) << std::endl;

    
}

//________________
double fitL2L3Corrections(double *x, double *par) {
    double xx = x[0];
    // double value = par[0] * exp( par[1] / pow(xx, par[2]) ) ; // Simple-efficiency
    // double value = par[0] * exp( (par[1] + par[6]*xx) / pow(xx, par[2]) ) + par[3] * TMath::Gaus(xx, par[4], par[5]); // Simple-efficiency + Gaussian
    double value = par[0] + par[1] / ( par[2] * xx + par[3] ) + par[4] * exp( par[5] / pow(xx, par[6]) ) + par[7] * exp( par[8] * log10(xx) + par[9]); // Simple-efficiency + Gaussian + Gaussian
    // double value = par[0] + par[1] / (pow(log10(xx),2) + par[2]) + par[3]*exp(-par[4]*pow((log10(xx)-par[5]),2)) - par[6]*exp(-par[7]*(pow(log10(xx)+par[8],2))); // From pp 2017
    return value;
}

//________________
void print1DArray(double *array, int size) {
    std::cout << "Array with size: " << size << std::endl;
    std::cout << "{";
    for (int i{0}; i<size; i++) {
        if (i == size - 1) {
            std::cout << array[i];
        }
        else {
            std::cout << array[i] << ", ";
        }
    }
    std::cout << "};" << std::endl;
}

//________________
void printCorrFactorArray(double *corrFactor, int jetNEtaL2L3StdBins, int jetNPtBins) {
    std::cout << "Correction factor array: " << std::endl;
    std::cout << "nEtaBins = " << jetNEtaL2L3StdBins << ", nPtBins = " << jetNPtBins << std::endl;
    std::cout << "{" << std::endl;
    for (int iEta{0}; iEta<jetNEtaL2L3StdBins; iEta++) {
        std::cout << "    {";
        for (int jPt{0}; jPt<jetNPtBins; jPt++) {
            if (jPt == jetNPtBins - 1) {
                std::cout << corrFactor[iEta*jetNPtBins + jPt];
            }
            else {
                std::cout << corrFactor[iEta*jetNPtBins + jPt] << ", ";
            }
        }
        if (iEta == jetNEtaL2L3StdBins - 1) {
            std::cout << "}" << std::endl;  
        }
        else {
            std::cout << "}," << std::endl;
        }
    }
    std::cout << "};\n" << std::endl;
}

//________________
bool histogramWithAtLeastThreeNonEmptyBins(TH1D *h) {
    int nBins = h->GetNbinsX();
    int nNonEmptyBins = 0;
    for (int i{1}; i<=nBins; i++) {
        if ( h->GetBinContent(i) > 0 ) {
            nNonEmptyBins++;
        }
    }
    return (nNonEmptyBins >= 3);
}

//________________
void findCorrections(TFile *f, int collisionSystem = 1, double collisionEnergy = 8.16, TString date = "20250129") {

    TString collSystemStr = (collisionSystem == 0) ? "pp" : (collisionSystem == 1) ? "pPb" : "PbPb";
    collSystemStr += Form("%d", int(collisionEnergy * 1000) );

    TH3D *h3D = (TH3D*)f->Get("hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning");
    if ( !h3D ) {
        std::cerr << "Histogram hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning not found" << std::endl;
        return;
    }

    // Standard binning

    // pt  : 1300 bins from 10. to 6510. GeV
    // eta : 52 bins from -5.2 to 5.2
    int jetPtBinNumbers[] = {
        1,   // 10.0
        2,   // 15.0
        3,   // 20.0
        4,   // 25.0
        6,   // 35.0
        8,   // 45.0
        10,  // 55.0
        12,  // 65.0
        14,  // 75.0
        16,  // 85.0
        18,  // 95.0
        21,  // 110.0
        25,  // 130.0
        29,  // 150.0
        35,  // 180.0
        41,  // 210.0
        49,  // 250.0
        59,  // 300.0
        71,  // 360.0
        83,  // 420.0
        99,  // 500.0
        119, // 600.0
        143, // 720.0
        171, // 860.0
        219, // 1100.0
        1300 // 6510.0
    };
    int jetNPtBins = sizeof(jetPtBinNumbers) / sizeof(int);
    std::cout << "pt bins: " << jetNPtBins << std::endl;

    double jetEtaL2L3StdVals[] = { -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, 
                                   -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, 
                                   -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
                                   -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, 
                                   -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, 
                                   -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
                                   -0.435, -0.348, -0.261, -0.174, -0.087,  0.000, 
                                    0.087,  0.174,  0.261,  0.348,  0.435,  0.522, 
                                    0.609,  0.696,  0.783,  0.879,  0.957,  1.044, 
                                    1.131,  1.218,  1.305,  1.392,  1.479,  1.566, 
                                    1.653,  1.740,  1.830,  1.930,  2.043,  2.172, 
                                    2.322,  2.500,  2.650,  2.853,  2.964,  3.139, 
                                    3.314,  3.489,  3.664,  3.839,  4.013,  4.191, 
                                    4.363,  4.538,  4.716,  4.889,  5.191 };
    int jetNEtaL2L3StdBins = sizeof(jetEtaL2L3StdVals)/sizeof(double);
    std::cout << "eta values: " << std::endl;
    print1DArray(jetEtaL2L3StdVals, jetNEtaL2L3StdBins);

    // Create array with correction factors
    double jetPtVals[jetNPtBins];
    for (int i{0}; i<jetNPtBins; i++) {
        if ( i < (jetNPtBins - 1 ) ) {
            jetPtVals[i] = 10. + (jetPtBinNumbers[i] - 1) * 5.;
        }
        else {
            jetPtVals[i] = 10. + jetPtBinNumbers[i] * 5.;
        }
    }
    std::cout << "pt values: " << std::endl;
    print1DArray(jetPtVals, jetNPtBins);

    double corrFactor[jetNEtaL2L3StdBins-1][jetNPtBins-1];

    TCanvas *cRawOverGenPt[jetNEtaL2L3StdBins];
    TH1D *hRawOverGenPt[jetNEtaL2L3StdBins][jetNPtBins];
    TH1D *hPtRaw[jetNEtaL2L3StdBins];

    TF1 *ratioFit[jetNEtaL2L3StdBins][jetNPtBins];

    TGraphErrors *gCorrFactorVsPt[jetNEtaL2L3StdBins];

    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Loop over eta bins
    int nColumns = 7;
    int nRawsPt = ( jetNPtBins % nColumns ) ? (jetNPtBins / nColumns) + 1 : jetNPtBins / nColumns;
    int nRawsEta = ( jetNEtaL2L3StdBins % nColumns ) ? (jetNEtaL2L3StdBins / nColumns) + 1 : jetNEtaL2L3StdBins / nColumns;

    // TCanvas *cMeanValues = new TCanvas("cMeanValues", "cMeanValues", 1200, 800);
    // cMeanValues->Divide(nColumns, nRawsEta);

    TCanvas *c = new TCanvas("c", "c", 1000, 800);

    TF1 *l2l3CorrectionsFit[jetNEtaL2L3StdBins];

    for (int iEta{0}; iEta<(jetNEtaL2L3StdBins-1); iEta++) {

        // Create canvas (to draw histograms)
        // cRawOverGenPt[iEta] = new TCanvas(Form("cRawOverGenPt_%d", iEta), Form("cRawOverGenPt_%d", iEta), 1200, 800);        
        // cRawOverGenPt[iEta]->Divide(nColumns, nRawsPt);
        hPtRaw[iEta] = (TH1D*)h3D->ProjectionY(Form("hPtRaw_%d", iEta), 1, -1, iEta+1, iEta+1);
        set1DStyle(hPtRaw[iEta], 2);
        hPtRaw[iEta]->Scale( 1./hPtRaw[iEta]->Integral() );

        // std::cout << "Eta bin: " << iEta << std::endl;
        // std::cout << "Expected projection range: " << jetEtaL2L3StdVals[iEta] << " - " << jetEtaL2L3StdVals[iEta+1] << std::endl;
        // std::cout << Form("Eta bin: %d, range: %.3f - %.3f", iEta, h3D->GetZaxis()->GetBinLowEdge(iEta + 1), h3D->GetZaxis()->GetBinUpEdge(iEta + 1)) << std::endl;
        // std::cout << "Eta bin center: " << h3D->GetZaxis()->GetBinCenter(iEta + 1) << std::endl;

        gCorrFactorVsPt[iEta] = new TGraphErrors();
        gCorrFactorVsPt[iEta]->SetName(Form("gCorrFactorVsPt_%d", iEta));


        // Draw and save ptRaw spectra
        // c->cd();
        // setPadStyle();
        // hPtRaw[iEta]->Draw();
        // hPtRaw[iEta]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV)");
        // hPtRaw[iEta]->GetYaxis()->SetTitle("1/N dN/dp_{T}^{jet}");
        // hPtRaw[iEta]->GetXaxis()->SetRangeUser(10., 400.);
        // c->SetLogy(1);
        // t.DrawLatexNDC(0.45, 0.8, Form("%.2f < #eta < %.2f", jetEtaL2L3StdVals[iEta], jetEtaL2L3StdVals[iEta+1]));
        // c->SaveAs(Form("%s/%s_ptRaw_%d.pdf", date.Data(), collSystemStr.Data(), iEta));
        

        // Loop over pT bins
        for (int jPt{0}; jPt<(jetNPtBins-1); jPt++) {

            corrFactor[iEta][jPt] = 0.;
            double meanX{0.}, meanXError{0.};
            double meanY{0.}, meanYError{0.};

            // Calculate mean X (ptRaw) and its error
            calculateMeanAndErrorInRange(hPtRaw[iEta], jetPtBinNumbers[jPt], jetPtBinNumbers[jPt+1], meanX, meanXError);

            // std::cout << "Pt bin: " << jPt << std::endl;
            // std::cout << "Expected projection range: " << jetPtVals[jPt] << " - " << jetPtVals[jPt+1] << std::endl;
            // std::cout << Form("Pt bin: %d, range: %.3f - %.3f", jPt, 
            //                   h3D->GetYaxis()->GetBinLowEdge(jetPtBinNumbers[jPt]), 
            //                   h3D->GetYaxis()->GetBinUpEdge(jetPtBinNumbers[jPt+1]-1)) << std::endl;

            // Make projection
            hRawOverGenPt[iEta][jPt] = (TH1D*)h3D->ProjectionX(Form("hRawOverGenPt_%d_%d", iEta, jPt), 
                                                               jetPtBinNumbers[jPt], jetPtBinNumbers[jPt+1]-1, iEta+1, iEta+1);
            set1DStyle(hRawOverGenPt[iEta][jPt], 2);
            hRawOverGenPt[iEta][jPt]->Scale( 1./hRawOverGenPt[iEta][jPt]->Integral() );

            //
            // If empty histogram (no entries for the given eta and pt bin)
            //
            if ( hRawOverGenPt[iEta][jPt]->GetEntries() <= 0 ) {
                // If not the first point, copy values from the previous point. Otherwise, set to 0.
                if ( gCorrFactorVsPt[iEta]->GetN() != 0 ) {
                    double tmpX{0.};
                    gCorrFactorVsPt[iEta]->GetPoint(gCorrFactorVsPt[iEta]->GetN()-1, tmpX, meanY);
                    meanYError = gCorrFactorVsPt[iEta]->GetErrorY(gCorrFactorVsPt[iEta]->GetN()-1);
                }
            } // if ( hRawOverGenPt[iEta][jPt]->GetEntries() <= 0 )
            //
            // Non-empty histogram
            //
            else { 

                //
                // Analyze only histograms with at least three non-empty bins
                // (otherwise large fluctuations)
                //
                if ( histogramWithAtLeastThreeNonEmptyBins( hRawOverGenPt[iEta][jPt] ) ) {
                    // Mean y and its error from histogram
                    meanY = hRawOverGenPt[iEta][jPt]->GetMean();
                    meanYError = hRawOverGenPt[iEta][jPt]->GetMeanError();
                    // std::cout << "Mean: " << meanY << std::endl;
                    // std::cout << "Mean error: " << meanYError << std::endl;
                
                    // Draw and save raw over ref pt and perform fit
                    c->cd();
                    setPadStyle();
                    hRawOverGenPt[iEta][jPt]->Draw();
                    ratioFit[iEta][jPt] = new TF1(Form("ratioFit_%d_%d", iEta, jPt), "gaus", 0., 1.7);
                    ratioFit[iEta][jPt]->SetLineColor(kRed);
                    ratioFit[iEta][jPt]->SetParameters(0.6, 1., 0.1);
                    hRawOverGenPt[iEta][jPt]->Fit(ratioFit[iEta][jPt], "MRQE");
                    double fitMeanY = ratioFit[iEta][jPt]->GetParameter(1);
                    double fitSigmaY = ratioFit[iEta][jPt]->GetParameter(2);
                    t.SetTextSize(0.04);
                    t.DrawLatexNDC(0.18, 0.8, Form("%.2f < #eta < %.2f", jetEtaL2L3StdVals[iEta], jetEtaL2L3StdVals[iEta+1]));
                    t.DrawLatexNDC(0.18, 0.7, Form("%.0f < p_{T}^{jet} < %.0f", 10.+(jetPtBinNumbers[jPt]-1) * 5., 10. + (jetPtBinNumbers[jPt+1]-1) * 5.));
                    t.DrawLatexNDC(0.6, 0.8, Form("#mu = %.3f #pm %.3f", ratioFit[iEta][jPt]->GetParameter(1), ratioFit[iEta][jPt]->GetParError(1)));
                    t.DrawLatexNDC(0.6, 0.7, Form("#sigma = %.3f #pm %.3f", ratioFit[iEta][jPt]->GetParameter(2), ratioFit[iEta][jPt]->GetParError(2)));
                    t.DrawLatexNDC(0.6, 0.6, Form("#chi^{2}/ndf = %.2f", ratioFit[iEta][jPt]->GetChisquare() / ratioFit[iEta][jPt]->GetNDF()));
                    c->SetLogy(0);
                    c->SaveAs( Form("%s/%s_rawOverGenPt_%d_pt_%d_%d.pdf", date.Data(), collSystemStr.Data(), iEta, 
                                (int)h3D->GetYaxis()->GetBinLowEdge(jetPtBinNumbers[jPt]), (int)h3D->GetYaxis()->GetBinUpEdge(jetPtBinNumbers[jPt+1]-1) ) );

                    // Use fitted values    
                    // meanY = fitMeanY;
                    // meanYError = fitSigmaY;

                    if ( !std::isnan(meanY) ) {
                        meanY = 1. / meanY;
                        meanYError = meanYError / (meanY * meanY);
                    }
                    else {
                        meanY = 0.;
                        meanYError = 0.;
                    }

                    // std::cout << Form("Filling values to the graph: %.3f, %.3f", meanX, meanY) << std::endl;
                    gCorrFactorVsPt[iEta]->AddPoint(meanX, meanY);
                    gCorrFactorVsPt[iEta]->SetPointError(gCorrFactorVsPt[iEta]->GetN()-1, meanXError, meanYError);
                } // if ( histogramWithAtLeastThreeNonEmptyBins ( hRawOverGenPt[iEta][jPt] ) )
                else {
                    // If not the first point, copy values from the previous point. Otherwise, set to 0.
                    if ( gCorrFactorVsPt[iEta]->GetN() != 0 ) {
                        double tmpX{0.};
                        gCorrFactorVsPt[iEta]->GetPoint(gCorrFactorVsPt[iEta]->GetN()-1, tmpX, meanY);
                        meanYError = gCorrFactorVsPt[iEta]->GetErrorY(gCorrFactorVsPt[iEta]->GetN()-1);
                    }
                }

            } // else { // Non-empty histogram
    
            // Draw
            // cRawOverGenPt[iEta]->cd(jPt+1);
            // setPadStyle();
            // hRawOverGenPt[iEta][jPt]->Draw();
            // t.DrawLatexNDC(0.45, 0.8, Form("%.2f < #eta < %.2f", jetEtaL2L3StdVals[iEta], jetEtaL2L3StdVals[iEta+1]));
            // t.DrawLatexNDC(0.45, 0.7, Form("%.0f < p_{T}^{jet} < %.0f", 10.+(jetPtBinNumbers[jPt]-1) * 5., 10. + (jetPtBinNumbers[jPt+1]-1) * 5.));

            // Remember, the correction factor is already inverted!!
            corrFactor[iEta][jPt] = meanY; 

        } // for (int jPt{0}; jPt<jetNPtBins; jPt++)

        // cMeanValues->cd(iEta+1);
        // setPadStyle();
        gCorrFactorVsPt[iEta]->SetMarkerStyle(20);
        gCorrFactorVsPt[iEta]->SetMarkerSize(1.2);
        gCorrFactorVsPt[iEta]->SetMarkerColor(kBlack);
        gCorrFactorVsPt[iEta]->SetLineColor(kBlack);
        gCorrFactorVsPt[iEta]->SetLineWidth(2);
        gCorrFactorVsPt[iEta]->Draw("AP");
        gCorrFactorVsPt[iEta]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV)");
        gCorrFactorVsPt[iEta]->GetXaxis()->SetRangeUser(0., 600.);
        gCorrFactorVsPt[iEta]->GetYaxis()->SetTitle("L2L3 Correction Factor");
        gCorrFactorVsPt[iEta]->GetYaxis()->SetRangeUser(0.8, 2.2);
        // t.DrawLatexNDC(0.45, 0.8, Form("%.2f < #eta < %.2f", jetEtaL2L3StdVals[iEta], jetEtaL2L3StdVals[iEta+1]));

        c->cd();
        setPadStyle();
        c->SetLogy(0);
        gCorrFactorVsPt[iEta]->Draw("AP");
        t.DrawLatexNDC(0.4, 0.8, Form("%.2f < #eta < %.2f", jetEtaL2L3StdVals[iEta], jetEtaL2L3StdVals[iEta+1]));

        // pp2017 fit
        // l2l3CorrectionsFit[iEta] = new TF1( Form("l2l3CorrectionsFit_%d", iEta), fitL2L3Corrections, 10., 400., 9);
        // l2l3CorrectionsFit[iEta]->SetParameters(12., 775., 0.5, 8.7, 4.3, -0.98, 0.35, 0.6, 6.9);

        // Current fit
        // l2l3CorrectionsFit[iEta] = new TF1( Form("l2l3CorrectionsFit_%d", iEta), fitL2L3Corrections, 10., 400., 10);
        // l2l3CorrectionsFit[iEta]->SetParameters(1.86529, -1.61543e+03, 1.62368e+01, 3.57560e+03, -4.32881e+05, -5.24556e+01, -1.93128e+15, -8.56230e-02, 3.99420e-01, 9.76207e-01 );
        // gCorrFactorVsPt[iEta]->Fit(l2l3CorrectionsFit[iEta], "MRE");

        c->SaveAs(Form("%s/%s_L2L3Correction_%d.pdf", date.Data(), collSystemStr.Data(), iEta));
    } // for (int i{0}; i<jetNEtaL2L3StdBins; i++)


    // Print correction factor array
    printCorrFactorArray(&corrFactor[0][0], jetNEtaL2L3StdBins-1, jetNPtBins-1);
}

//________________
/// Macro derives L2L3 corrections from the ratio of reco raw pt / ref pt.
/// Raw pt is the reconstructed jet pt, ref pt is the generated jet pt.
/// The ratio is the L2L3 correction factor.
/// The correction factors are plotted for each eta bin as a function of pt,
/// then total integrals is splitted into several pT intervals with 5-10% width of the total area.
/// For each eta bin and pT interval, the mean value of the correction factor is calculated.
/// The mean values is plotted as a function of pT, and fitted with a function.
/// The correction values are stored as a 2D array C++ array, and written to a file.
void deriveL2L3Corrections() {

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

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    int collSystem = 0;
    double energy = 5.02;

    // pPb8160 embedding
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", uname.Data()) );
    TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_eta25.root", uname.Data()) );
    if ( !inputFile || inputFile->IsZombie() ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_def_ak4_eta25.root", uname.Data()) << std::endl;
        return;
    }
    collSystem = 1;
    energy = 8.16;

    findCorrections(inputFile, collSystem, energy, date.Data());
}