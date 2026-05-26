#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom3.h"

//________________
void setPadStyle() {
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.15);
}

void testSmearing() {
    TRandom3 rand(12345);

    auto inFile = TFile::Open("oFilePtSpectra.root", "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: could not open input file oFilePtSpectra.root" << std::endl;
        return;
    }
    TH1D *hPtGen = dynamic_cast<TH1D*>( inFile->Get("hPt") );
    if (!hPtGen) {
        std::cerr << "Error: histogram hPt not found in oFilePtSpectra.root" << std::endl;
        inFile->Close();
        return;
    }
    hPtGen->SetDirectory(0); // Detach from file to prevent deletion when file is closed
    hPtGen->SetName("hPtGen");
    inFile->Close();
    TH1D *hPtClone = dynamic_cast<TH1D*>( hPtGen->Clone("hPt_clone") );
    if (!hPtClone) {
        std::cerr << "Error: failed to clone hPt" << std::endl;
        return;
    }
    auto hPt = std::unique_ptr<TH1D>( hPtClone );
    if (!hPt->GetSumw2()) hPt->Sumw2();
    auto hPtSmeared = std::make_unique<TH1D>("hPtSmeared", "Smeared Pt", hPt->GetNbinsX(), hPt->GetXaxis()->GetXmin(), hPt->GetXaxis()->GetXmax()); 
    hPtSmeared->Sumw2();

    auto fSmear = std::make_unique<TF1>("fSmear", "sqrt(max(0., [0] + [1]/x))", 20., 800.);
    fSmear->SetParameters(0.00183, 0.755);
    fSmear->SetLineColor(kBlue);
    fSmear->SetLineWidth(2);

    auto hSmeared2D = std::make_unique<TH2D>("hSmeared2D", "Smeared 2D;p_{T}^{gen} (GeV);p_{T}^{smeared}/p_{T}^{gen}", 100, 0., 1000., 200, 0., 2.);
    hSmeared2D->Sumw2();

    long int n2gen = 1000000000;
    for (long int i = 0; i < n2gen; ++i) {
        const double randomPt = hPt->GetRandom();
        const double evalX = std::max(1.0, randomPt);
        const double smearFactor = fSmear->Eval(evalX);
        const double smearedPt = randomPt * rand.Gaus(1., smearFactor);
        hSmeared2D->Fill(randomPt, smearedPt/randomPt);
        hPtSmeared->Fill(smearedPt);
    }

    auto c = std::make_unique<TCanvas>("c", "c", 1000, 1000);
    setPadStyle();
    hSmeared2D->Draw("colz");
    c->SaveAs("smeared2D.pdf");

    hSmeared2D->FitSlicesY();
    auto hJER = dynamic_cast<TH1D*>( gDirectory->Get( Form("%s_2", hSmeared2D->GetName()) ) );
    auto fSmearFit = std::make_unique<TF1>("fSmearFit", "sqrt(max(0., [0] + [1]/x))", 30., 600.);
    if (!hJER) {
        std::cerr << "Error: failed to retrieve JER histogram from FitSlicesY." << std::endl;
    } else {
        hJER->SetDirectory(0);
        hJER->SetLineColor(kBlack);
        hJER->SetLineWidth(2);
        hJER->SetMarkerStyle(20);
        hJER->SetMarkerSize(1.3);
        hJER->Fit(fSmearFit.get(), "MRE0");

        hJER->Draw();
        hJER->GetXaxis()->SetRangeUser(10., 300.);
        hJER->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
        hJER->GetYaxis()->SetRangeUser(0., 0.3);
        hJER->GetYaxis()->SetTitle("JER (p_{T}^{smeared}/p_{T}^{gen})");
        fSmearFit->Draw("same");
        fSmear->Draw("same");
        TLatex t;
        t.SetTextSize(0.03);
        t.SetTextFont(42);
        t.DrawLatexNDC(0.16, 0.85, "Testing smearing procedure");
        t.DrawLatexNDC(0.16, 0.80, Form("Initial smearing function: #sqrt{%.5f + %.5f/x}", 0.00183, 0.755) );
        t.DrawLatexNDC(0.16, 0.75, Form("JER fit: #sqrt{%.5f + %.5f/x}", fSmearFit->GetParameter(0), fSmearFit->GetParameter(1)) );
        c->SaveAs("smeared2D_JER.pdf");
    }

        // Write histograms to output ROOT file
        auto outFile = std::unique_ptr<TFile>( TFile::Open("testSmearing_out.root", "RECREATE") );
        if ( outFile && !outFile->IsZombie() ) {
            outFile->cd();
            // Write histograms into the file without transferring ownership
            hSmeared2D->Write("hSmeared2D");
            if (hJER) hJER->Write("hJER");
            hPt->Write("hPt");
            hPtSmeared->Write("hPtSmeared");
            outFile->Close();
        } else {
            std::cerr << "Error: could not create output file testSmearing_out.root" << std::endl;
        }
}
