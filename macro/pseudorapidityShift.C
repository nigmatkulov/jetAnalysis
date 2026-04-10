#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"

#include <iostream>
#include <array>
#include <limits>
#include <cmath>

float ptHat;
float vz;

int ngen;
float genpt[5000];
float geneta[5000];

//________________
// Event weight calculation for pPb8160 MC samples
double eventWeight(float ptHat, float vz, TF1 *fVzWeight, Long64_t eventsInSample) {

    double genWeight = 1.0;
    double vzWeight = 1./ fVzWeight->Eval( vz );
    // Magic numbers are (cross section x Nevents generated)
    if (ptHat > 15.0 && ptHat <= 30.)       { genWeight = 1.0404701e-06 * 961104 ; }
    else if (ptHat > 30. && ptHat <= 50.)   { genWeight = 7.7966624e-08 * 952110 ; }
    else if (ptHat > 50. && ptHat <= 80.)   { genWeight = 1.0016052e-08 * 952554 ; }
    else if (ptHat > 80. && ptHat <= 120.)  { genWeight = 1.3018269e-09 * 996844 ; }
    else if (ptHat > 120.&& ptHat <= 170.)  { genWeight = 2.2648493e-10 * 964681 ; }
    else if (ptHat > 170. && ptHat <= 220.) { genWeight = 4.0879112e-11 * 999260 ; }
    else if (ptHat > 220. && ptHat <= 280.) { genWeight = 1.1898939e-11 * 964336 ; }
    else if (ptHat > 280. && ptHat <= 370.) { genWeight = 3.3364433e-12 * 995036 ; }
    else if (ptHat > 370. && ptHat <= 460.) { genWeight = 7.6612402e-13 * 958160 ; }
    else if (ptHat > 460. && ptHat <= 540.) { genWeight = 2.1341026e-13 * 981427 ; }
    else if (ptHat > 540.)                  { genWeight = 7.9191586e-14 * 1000000; }
    genWeight /= eventsInSample;

    return genWeight * vzWeight;
} 

//________________
void setupChain(const char *inputFileName, TChain *mainTree, TChain *eventTree, TChain *jetTree) {
    // Merge trees into a single TChain
    mainTree->Add(inputFileName);
    eventTree->Add(inputFileName);
    jetTree->Add(inputFileName);

    mainTree->AddFriend(eventTree);
    mainTree->AddFriend(jetTree);

    // Disable all branches first
    mainTree->SetBranchStatus("*", 0);

    // Enable only the branches needed for this analysis
    mainTree->SetBranchStatus("vz", 1);
    mainTree->SetBranchStatus("pthat", 1);

    mainTree->SetBranchStatus("ngen", 1);
    mainTree->SetBranchStatus("genpt", 1);
    mainTree->SetBranchStatus("geneta", 1);

    // Set branch addresses
    mainTree->SetBranchAddress("vz", &vz);
    mainTree->SetBranchAddress("pthat", &ptHat);

    mainTree->SetBranchAddress("ngen", &ngen);
    mainTree->SetBranchAddress("genpt", &genpt);
    mainTree->SetBranchAddress("geneta", &geneta);
}

//________________
float etaCM(float etaLab, float etaShift, bool isPbGoing) {
    float etaCM;
    if (isPbGoing) {
        etaCM = -1.0 * (etaLab + etaShift);
    }
    else {
        etaCM = etaLab - etaShift;
    }
    return etaCM;
}

//________________
void pseudorapidityShift() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Define the standard eta bins for L2L3 corrections
    bool isPbGoing = false;
    const int ptHatSample = 15;
    bool isPythia = true;
    TString generator = isPythia ? "pythia" : "embedding";
    TString tag = isPythia ? "unembedded" : "embedded";
    TString direction = isPbGoing ? "Pbgoing" : "pgoing";

    const int nEtaShifts = 8;
    static constexpr std::array<float, nEtaShifts> etaShift{0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470};

    float ptHatRange[2];
    if (ptHatSample == 15)       { ptHatRange[0] = 15.; ptHatRange[1] = 30.; }
    else if (ptHatSample == 30)  { ptHatRange[0] = 30.; ptHatRange[1] = 50.; }
    else if (ptHatSample == 50)  { ptHatRange[0] = 50.; ptHatRange[1] = 80.; }
    else if (ptHatSample == 80)  { ptHatRange[0] = 80.; ptHatRange[1] = 120.; }
    else if (ptHatSample == 120) { ptHatRange[0] = 120.; ptHatRange[1] = 170.; }
    else if (ptHatSample == 170) { ptHatRange[0] = 170.; ptHatRange[1] = 220.; }
    else if (ptHatSample == 220) { ptHatRange[0] = 220.; ptHatRange[1] = 280.; }
    else if (ptHatSample == 280) { ptHatRange[0] = 280.; ptHatRange[1] = 370.; }
    else if (ptHatSample == 370) { ptHatRange[0] = 370.; ptHatRange[1] = 460.; }
    else if (ptHatSample == 460) { ptHatRange[0] = 460.; ptHatRange[1] = 540.; }
    else if (ptHatSample == 540) { ptHatRange[0] = 540.; ptHatRange[1] = std::numeric_limits<float>::max(); }

    TString fname = Form("~/cernbox/ana/pPb8160/%s/%s/forest/HiForestSkim_pPb_MC_pthat%d_%s_%s.root", 
                         generator.Data(), direction.Data(), ptHatSample, direction.Data(), tag.Data());


    TFile *f = TFile::Open( fname );
    TChain *mainTree = new TChain("hltanalysis/HltTree");
    TChain *eventTree = new TChain("hiEvtAnalyzer/HiTree");
    TChain *jetTree = new TChain("ak4PFJetAnalyzer/t");
    setupChain(fname.Data(), mainTree, eventTree, jetTree);

    TF1 *fVzWeight = new TF1("fVzWeight", "pol8", -15.1, 15.1);
    fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);

    //////////////////////////
    //       Histograms     //
    //////////////////////////

    TH1D *hPtHatUnweighted = new TH1D("hPtHatUnweighted", "#hat{p}_{T} unweighted;#hat{p}_{T} (GeV);dN/d#hat{p}_{T}", 100, 0., 1000.);
    hPtHatUnweighted->Sumw2();
    TH1D *hPtHat = new TH1D("hPtHat", "#hat{p}_{T};#hat{p}_{T} (GeV);dN/d#hat{p}_{T}", 100, 0., 1000.);
    hPtHat->Sumw2();
    TH1D *hGenJetEtaCMShiftedUnweighted[nEtaShifts];
    TH1D *hGenJetEtaCMShifted[nEtaShifts];
    for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
        hGenJetEtaCMShiftedUnweighted[iShift] = new TH1D(Form("hGenJetEtaCMShiftedUnweighted_%d", iShift), 
                                               Form("Gen jet #eta_{CM} with #Delta#eta = %.4f (unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                               480, -6., 6.); // bin width 0.05
        hGenJetEtaCMShiftedUnweighted[iShift]->Sumw2();
        hGenJetEtaCMShifted[iShift] = new TH1D(Form("hGenJetEtaCMShifted_%d", iShift), 
                                               Form("Gen jet #eta_{CM} with #Delta#eta = %.4f;#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                               480, -6., 6.); // bin width 0.05
        hGenJetEtaCMShifted[iShift]->Sumw2();
    }


    double weight{0.};

    //////////////////////////////
    //        Processing        //
    //////////////////////////////

    Long64_t nEntries = mainTree->GetEntries();
    std::cout << "Number of entries: " << nEntries << std::endl;
    
    int nEventsProcessed{0};
    int nGoodEvents{0};

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {

        mainTree->GetEntry(iEntry);
        nEventsProcessed++;

        // std::cout << "\n========================================\n";
        // std::cout << Form("Event: %d ptHat = %.1f GeV, vz = %.2f cm", iEntry, ptHat, vz) << std::endl;
        // std::cout << Form("%.1f <= ptHat <= %.1f ? %s", ptHatRange[0], ptHatRange[1], 
        //                  (ptHat >= ptHatRange[0] && ptHat <= ptHatRange[1]) ? "true" : "false") << std::endl;
        // std::cout << Form("|vz| <= 15 cm ? %s", (std::abs(vz) <= 15.) ? "true" : "false") << std::endl;

        // Check the event is satisfies basic event selection
        if (ptHat < ptHatRange[0] || ptHat > ptHatRange[1]) continue; 
        if ( std::abs(vz) > 15.) continue;

        nGoodEvents++;

        if (isPbGoing) vz = -vz; 
        weight = eventWeight(ptHat, vz, fVzWeight, nEntries);

        hPtHatUnweighted->Fill(ptHat);
        hPtHat->Fill(ptHat, weight);

        // Loop over generated jets in the event
        for (int iGenJet = 0; iGenJet < ngen; ++iGenJet) {
            // Skip low-pT jets
            if (genpt[iGenJet] < 40.) continue;
            // std::cout << "Gen jet " << iGenJet << ": pT = " << genpt[iGenJet] << " GeV, eta = " << geneta[iGenJet] << std::endl;

            // Loop over different eta shifts and fill histograms
            for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
                float etaCMShifted = etaCM(geneta[iGenJet], etaShift[iShift], isPbGoing);
                hGenJetEtaCMShiftedUnweighted[iShift]->Fill(etaCMShifted);
                hGenJetEtaCMShifted[iShift]->Fill(etaCMShifted, weight);
            } // for (int iShift = 0; iShift < nEtaShifts; ++iShift)
        } // for (int iGenJet = 0; iGenJet < ngen; ++iGenJet)
    } // for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)

    std::cout << Form("Processed %d events, of which %d passed selection.", nEventsProcessed, nGoodEvents) << std::endl;
    
    // Output file
    int compressionSetting = 208; // LZMA compression
    TFile *fOut = TFile::Open( Form("eta_shift/%s_%s_ptHat%d.root", 
                                    generator.Data(), direction.Data(), ptHatSample), 
                               "RECREATE", "", compressionSetting);
    hPtHatUnweighted->Write();
    hPtHat->Write();
    for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
        hGenJetEtaCMShiftedUnweighted[iShift]->Write();
        hGenJetEtaCMShifted[iShift]->Write();
    }
    fOut->Close();
    f->Close();
}
