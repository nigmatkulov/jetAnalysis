#include "call_libraries.h"
#include "JetCorrector.h"   // reader for JEC
#include "JetUncertainty.h" // reader for JEU

// Dummy global variables
// event quantities
float vertexz; // event z vertex
int hiBin; // event centrality (used if use_centrality = true in input_variables.h)

// trigger quantities
int jet_trigger_bit; // jet HLT path trigger used for analysis (jet_trigger variable in input_variables.h)

// reco jets
int nref;         // number of jets
float jteta[9999]; // jet eta
float jtphi[9999]; // jet phi
float rawpt[9999]; // jet pT without JEC
float trackMax[9999]; // track maximum pT in a jet
// reco tracks
int ntrk;                 // number of track
float trkpt[9999];       // track pT
float trketa[9999];      // track eta
float trkphi[9999];      // track phi
float trkpterr[9999];    // track pT error (uncertainty)
float trkdcaxy[9999];    // track dxy impact parameter (transverse distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaz[9999];     // track dz impact parameter (longitudinal distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaxyerr[9999]; // track dxy error (uncertainty)
float trkdcazerr[9999];  // track dxy error (uncertainty)
float trkchi2[9999];     // track reconstruction chi2 of the fitting
float pfEcal[9999];      // particle flow energy deposit in ECAL
float pfHcal[9999];      // particle flow energy deposit in HCAL
float trkmva[9999];      // track mva for each step
unsigned char trkalgo[9999];       // track algorithm/step
unsigned char trkndof[9999];       // track number of degrees of freedom in the fitting 
int trkcharge[9999];     // track charge
unsigned char trknhits[9999];      // number of hits in the tracker
unsigned char trknlayer[9999];     // number of layers with measurement in the tracker
bool highpur[9999];      // tracker steps MVA selection

// events quantities from gen
float weight; // event weight --> pthat weight
float pthat;  // pthat (initial parton pT)

// gen jets
int ngen;             // number of gen jets
float gen_jtpt[9999];  // gen jet pT
float gen_jteta[9999]; // gen jet eta
float gen_jtphi[9999]; // gen jet phi
float gen_jtetaX[9999]; // gen jet eta
float gen_jtphiX[9999]; // gen jet phi

// matched jets
float refpt[9999]; // jet pT matched with Gen pT
float refeta[9999]; // jet eta matched with Gen eta
float refphi[9999]; // jet phi matched with Gen phi
int refparton_flavor[9999]; // jet phi matched with Gen phi
int refparton_flavorForB[9999]; // jet phi matched with Gen phi

// gen tracks
std::vector<float> *gen_trkpt = 0;  // gen particle pT
std::vector<float> *gen_trketa = 0; // gen particle eta
std::vector<float> *gen_trkphi = 0; // gen particle phi
std::vector<int> *gen_trkchg = 0;   // gen particle charge
std::vector<int> *gen_trkpid = 0;   // gen particle pid
std::vector<int> *gen_trksube = 0;   // gen particle pid

//________________
Bool_t isGoodEvent(Float_t vz, Float_t ptHat) {
    return TMath::Abs(vz) < 15.;
}

//________________
Float_t get_event_weight(Float_t vz, int mult, Float_t weighttree, Float_t leadjetpt){
    Float_t vzweight = 1.0;
    Float_t multweight = 1.0;
    Float_t evtweight = 1.0; // = weighttree
    Float_t multefficiency = 1.0;
    Float_t jetefficiency = 1.0;		
    Float_t totalweight = 1.0;

    return evtweight * multweight * vzweight * multefficiency * jetefficiency;
}

//________________
int get_Ntrkoff(int size, float *eta, 
                float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, 
                float *dcaxyerr, float *dcaz, float *dcazerr, float* chi2, 
                unsigned char* ndof, unsigned char* nlayer, unsigned char* nhits, 
                unsigned char* algo, float* mva) {

    int Ntrk_off{0};

    // Loop over tracks
    //std::cout << "Total number of tracks: " << size;
    for (int iTrk=0; iTrk<size; iTrk++) { 
        if ( TMath::Abs( eta[iTrk] ) > 2.4) continue; 
        if (TMath::Abs( charge[iTrk] ) == 0)continue;
        if (TMath::Abs( charge[iTrk] ) == 0)continue;
        if (hp[iTrk] != 1) continue;
        if (TMath::Abs( pterr[iTrk]/pt[iTrk] ) >= 0.1) continue;
        if (TMath::Abs( dcaxy[iTrk]/dcaxyerr[iTrk] ) >= 3.0) continue;
        if (TMath::Abs( dcaz[iTrk]/dcazerr[iTrk] ) >= 3.0) continue;
        double calomatching = ((pfEcal[iTrk]+pfHcal[iTrk])/cosh(eta[iTrk]))/pt[iTrk];

        if( pt[iTrk] <= 0.5)  continue; 
        if( (chi2[iTrk]/ndof[iTrk])/Double_t(nlayer[iTrk]) >= 0.18 ) continue;
        if(nhits[iTrk] < 11) continue;
        if(pt[iTrk] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
        if(algo[iTrk]==6 && mva[iTrk]<0.98) continue;

        Ntrk_off=Ntrk_off+1;
    }
    //std::cout << "\tNtracks selected: " << Ntrk_off << std::endl;
    return Ntrk_off;
}

//________________
void setupChains(TString input, TChain *hltChain, TChain *eveChain, TChain *partFlowChain, 
                 TChain *trkChain, Bool_t useMC, TChain *genTrkChain) {

    // Check input exists
    if ( input.Length() <= 0 ) {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        exit(0);
    }
    // Normail input
    else {
        // If input is a single ROOT file
        if ( input.Index(".root") > 0 ) {
            std::cout << Form( "Adding %s file to chains\n", input.Data() );
            hltChain->Add( input.Data() );
            eveChain->Add( input.Data() );
            partFlowChain->Add( input.Data() );
            trkChain->Add( input.Data() );
            if ( useMC ) {
                genTrkChain->Add( input.Data() );
            }
        }
        // Assuming that list of files is provided instead of a single file
        else {
            std::ifstream inputStream( input.Data() );

            if ( !inputStream ) std::cout << Form( "ERROR: Cannot open file list: %s\n", input.Data() );
            Int_t nFiles = 0;
            std::string file;
            size_t pos;
            while ( getline( inputStream, file ) ) {
                // NOTE: our external formatters may pass "file NumEvents"
                //       Take only the first part
                //cout << "DEBUG found " <<  file << endl;
                pos = file.find_first_of(" ");
                if ( pos != std::string::npos ) file.erase( pos, file.length() - pos );
                //cout << "DEBUG found [" <<  file << "]" << endl;

                // Check that file is of a correct name
                if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos &&
                     file.find("AOD") != std::string::npos ) {
                    
                    // Open file
                    TFile* ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if ( ftmp && !ftmp->IsZombie() && ftmp->GetNkeys() ) {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str() );
                        // Adding file to chains
                        hltChain->Add( file.c_str() );
                        eveChain->Add( file.c_str() );
                        partFlowChain->Add( file.c_str() );
                        trkChain->Add( file.c_str() );
                        if ( useMC ) {
                            genTrkChain->Add( file.c_str() );
                        }
                        ++nFiles;
                    } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp) {
                        ftmp->Close();
                    } //if (ftmp)
                } //if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos && file.find("AOD") != std::string::npos )
            } //while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
        } // else {   if file list
    } // else {   if normal input

    // Connect chains between each other
    hltChain->AddFriend(eveChain);
    hltChain->AddFriend(partFlowChain);
    hltChain->AddFriend(trkChain);
}

//________________
void setupBranches(TChain *tree, Bool_t useMC) {

    // Disable all branches - this is important while reading big files
    tree->SetBranchStatus("*", 0); 

    // enable branches of interest -> see definition of each variables above

    // Event quantities
    tree->SetBranchStatus("vz", 1);
    tree->SetBranchStatus("hiBin", 1); //centrality only for PbPb and XeXe
    tree->SetBranchAddress("vz", &vertexz);
    tree->SetBranchAddress("hiBin", &hiBin); 

    if( useMC ) {
        tree->SetBranchStatus("weight", 1);
        tree->SetBranchStatus("pthat", 1); 
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("pthat", &pthat);
    }

    // Jet quantities
    tree->SetBranchStatus("nref", 1);
    tree->SetBranchStatus("rawpt", 1);
    tree->SetBranchStatus("trackMax", 1);
    tree->SetBranchStatus("jteta", 1);
    tree->SetBranchStatus("jtphi", 1);

    tree->SetBranchAddress("nref", &nref);
    tree->SetBranchAddress("rawpt", &rawpt);
    tree->SetBranchAddress("trackMax", &trackMax);
    tree->SetBranchAddress("jteta", &jteta);
    tree->SetBranchAddress("jtphi", &jtphi);

    // Gen jet quantities
    if ( useMC ) {
        tree->SetBranchStatus("ngen", 1);
        tree->SetBranchStatus("genpt", 1);

        tree->SetBranchStatus("geneta", 1);
        tree->SetBranchStatus("genphi", 1);
        tree->SetBranchStatus("WTAgeneta", 1);
        tree->SetBranchStatus("WTAgenphi", 1);

        tree->SetBranchAddress("ngen", &ngen);
        tree->SetBranchAddress("genpt", &gen_jtpt);

        tree->SetBranchAddress("geneta", &gen_jteta);
        tree->SetBranchAddress("genphi", &gen_jtphi);
        tree->SetBranchAddress("WTAgeneta", &gen_jtetaX);
        tree->SetBranchAddress("WTAgeneta", &gen_jtphiX); 

    }
    // Jet-matching quantities
    if ( useMC ) {
        tree->SetBranchStatus("refpt", 1);
        tree->SetBranchAddress("refpt", &refpt);
        tree->SetBranchStatus("refeta", 1);
        tree->SetBranchAddress("refeta", &refeta);
        tree->SetBranchStatus("refphi", 1);
        tree->SetBranchAddress("refphi", &refphi);
        tree->SetBranchStatus("refparton_flavor", 1);
        tree->SetBranchAddress("refparton_flavor", &refparton_flavor);
        tree->SetBranchStatus("refparton_flavorForB", 1);
        tree->SetBranchAddress("refparton_flavorForB", &refparton_flavorForB);
    }

    // Track quantities
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchStatus("trkPtError", 1);
    tree->SetBranchStatus("trkDxy1", 1);
    tree->SetBranchStatus("trkDxyError1", 1);
    tree->SetBranchStatus("trkDz1", 1);
    tree->SetBranchStatus("trkDzError1", 1);
    tree->SetBranchStatus("trkChi2", 1);
    tree->SetBranchStatus("trkNdof", 1);
    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchStatus("trkNHit", 1);
    tree->SetBranchStatus("trkNlayer", 1);
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchStatus("pfEcal", 1);
    tree->SetBranchStatus("pfHcal", 1);

    tree->SetBranchAddress("nTrk", &ntrk);
    tree->SetBranchAddress("trkPt", &trkpt);
    tree->SetBranchAddress("trkEta", &trketa);
    tree->SetBranchAddress("trkPhi", &trkphi);
    tree->SetBranchAddress("trkPtError", &trkpterr);
    tree->SetBranchAddress("trkDxy1", &trkdcaxy);
    tree->SetBranchAddress("trkDxyError1", &trkdcaxyerr);
    tree->SetBranchAddress("trkDz1", &trkdcaz);
    tree->SetBranchAddress("trkDzError1", &trkdcazerr);
    tree->SetBranchAddress("trkChi2", &trkchi2);
    tree->SetBranchAddress("trkNdof", &trkndof);
    tree->SetBranchAddress("trkCharge", &trkcharge);
    tree->SetBranchAddress("trkNHit", &trknhits);
    tree->SetBranchAddress("trkNlayer", &trknlayer);
    tree->SetBranchAddress("highPurity", &highpur);
    tree->SetBranchAddress("pfEcal", &pfEcal);
    tree->SetBranchAddress("pfHcal", &pfHcal);


    tree->SetBranchStatus("trkMVA", 1);
    tree->SetBranchStatus("trkAlgo", 1);
    tree->SetBranchAddress("trkMVA", &trkmva);
    tree->SetBranchAddress("trkAlgo", &trkalgo);

    // Gen particle quantities
    if( useMC ){
        tree->SetBranchStatus("pt", 1);
        tree->SetBranchStatus("eta", 1);
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchStatus("chg", 1);
        tree->SetBranchStatus("pdg", 1);
        tree->SetBranchStatus("sube", 1);

        tree->SetBranchAddress("pt", &gen_trkpt);
        tree->SetBranchAddress("eta", &gen_trketa);
        tree->SetBranchAddress("phi", &gen_trkphi);
        tree->SetBranchAddress("chg", &gen_trkchg);
        tree->SetBranchAddress("pdg", &gen_trkpid);
        tree->SetBranchAddress("sube", &gen_trksube);
    }
}

//________________
void setupInput(TString input, TChain *hltChain, TChain *eveChain, TChain *partFlowChain, 
                TChain *trkChain, Bool_t useMC, TChain *genTrkChain) {
    setupChains(input, hltChain, eveChain, partFlowChain, trkChain, useMC, genTrkChain);
    setupBranches(hltChain, useMC);
}

/**
 * @brief A macro that provides simle reading of the CMS Forest AOD
 * 
 * @param argc Provides number of arguments (0 - first argument always corresponds to the program name)
 * @param argv Array of arguments
 * 
 * @return Zero if finished normally
 */
int main(int argc, char* argv[]) {
    std::cout << "Executing simpleReadForest macro" << std::endl;

    // Set default values for arguments
    TString inFileName = "../../../data/HiForestAOD_PbPbMC2018skim_10.root";
    TString oFileName = "oTestSimpleReadForest.root";
    Long64_t nEventsToRead = 500;

    // Read input argument list 
    if (argc > 1) inFileName = argv[1];
    if (argc > 2) oFileName = argv[2];

    TString jet_collection = "akCs4PFJetAnalyzer";

    // Create chains to fill
    TChain *hltTree = new TChain("hltanalysis/HltTree");
    TChain *eventTree = new TChain("hiEvtAnalyzer/HiTree");
    TChain *partFlowTree = new TChain(Form("%s/t",jet_collection.Data()));
    TChain *trkTree = new TChain("ppTrack/trackTree");
    Bool_t useMC = kTRUE;
    TChain *genTrkTree{nullptr};
    if ( useMC ) genTrkTree = new TChain("HiGenParticleAna/hi");

    // Setup readout
    setupInput(inFileName, hltTree, eventTree, partFlowTree, trkTree, useMC, genTrkTree);

    // Read JEC file
    // std::vector< std::string > correctionFiles;
    // correctionFiles.push_back(Form("aux_correctionFiles/PbPb_5020/JEC/%s",JEC_file.Data()));
    // JetCorrector JEC(correctionFiles);

    // Histograms
    Int_t multBins = 1800;
    Double_t multRange[2] = {-0.5, 1799.5};
    Int_t hiBinBins = 203;
    Double_t hiBinRange[2] = {-1.5, 201.5};
    Int_t centralityBins = 101;
    Double_t centralityRange[2] = {-0.5, 100.5};
    Int_t weightBins = 220;
    Double_t weightRange[2] = {-1.1, 1.1};
    Int_t ptHatBins = 100;
    Double_t ptHatRange[2] = {0., 1000.};
    TH1D *hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", 320, -40., 40.);
    TH1D *hMult = new TH1D("hMult","Charged particle multiplicity;Multiplicity;Entries", 
                            multBins, multRange[0], multRange[1]);
    TH1D *hHiBin = new TH1D("hHiBin","HiBin a.k.a. centrality;HiBin;Entries", 
                            hiBinBins, hiBinRange[0], hiBinRange[1]);
    TH1D *hWeight = new TH1D("hWeight","Event weight;weight;Entries", 
                             weightBins, weightRange[0], weightRange[1]);
    TH1D *hPtHat = new TH1D("hPtHat","#hat{p_{T}};#hat{p_{T}} (GeV/c);Entries", 
                            ptHatBins, ptHatRange[0], ptHatRange[1]);
    TH1D *hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries",
                                 centralityBins, centralityRange[0], centralityRange[1]);
    TH2F *hCentralityVsMultiplicity = new TH2F("hCentralityVsMultiplicity","Centrality vs. multiplicity;Multiplicity;Centrality (%);Entries", 
                                                multBins, multRange[0], multRange[1], 
                                                centralityBins, centralityRange[0], centralityRange[1]);
    //TH1D *hTrackPt = new TH1D("hTrackPt","Track transverse momentum;p_{T} (GeV/c)")
    
    TH2F *hCentVsMultPartFlowJet[3];

    for (Int_t i=0; i<3; i++) {
        Int_t ptCutLow = 50 + 50 * i;
        hCentVsMultPartFlowJet[i] = new TH2F(Form("hCentVsMultPartFlowJet_%d",i),
                                             Form("Centrality vs. multiplicity for PartFlowJets with p_{T}>%d;Multiplicity;Centrality;Entries", ptCutLow),
                                             multBins, multRange[0], multRange[1], 
                                             centralityBins, centralityRange[0], centralityRange[1]);
    }
    TH1D *hRecoJetPt = new TH1D("hRecoJetPt","Reconstructed jet p_{T};p_{T}^{reco} (GeV/c);Entries", 75, 15., 315.);

    // Retrieve number of events to read
    nEventsToRead = hltTree->GetEntries(); 
	std::cout << "Total number of events in files: "<< nEventsToRead << std::endl;

    // Start event loop
    for (Long64_t iEvent{0}; iEvent < nEventsToRead; iEvent++) {

        // Load event
        hltTree->GetEntry(iEvent);
        hVz->Fill( vertexz );
        Int_t mult = get_Ntrkoff(ntrk, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, 
                                 trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, 
                                 trknhits, trkalgo, trkmva);
        hMult->Fill( mult );
        hHiBin->Fill( hiBin );
        hWeight->Fill( weight );
        hPtHat->Fill( pthat );

        Int_t centrality = 100 - Int_t(Double_t(200 - hiBin) / 2.);
        hCentrality->Fill( centrality );
        hCentralityVsMultiplicity->Fill( mult, centrality );

        // Get the event weight
        Double_t event_weight = get_event_weight(vertexz, mult, weight, pthat);

        // Loop over reconstructed jets
        Bool_t isGoodJetLow = kFALSE;
        Bool_t isGoodJetMid = kFALSE;
        Bool_t isGoodJetHi = kFALSE;
        for (Int_t iJet{0}; iJet<nref; iJet++) {
            hRecoJetPt->Fill( rawpt[iJet] );
            if ( rawpt[iJet] > 50 ) {
                isGoodJetLow = kTRUE;
                if ( rawpt[iJet] > 100 ) {
                    isGoodJetMid = kTRUE;
                    if ( rawpt[iJet] > 150 ) {
                        isGoodJetHi = kTRUE;
                    } // if ( rawpt[iJet] > 150 ) {
                } // if ( rawpt[iJet] > 100 ) {
            } // if ( rawpt[iJet] > 50 ) {

            // // Cut for jets with only very low pT particles
            // if (trackMax[j] / rawpt[j] < 0.05) continue; 
            // // Cut for jets where all the pT is taken by one track
            // if (trackMax[j] / rawpt[j] > 0.95) continue;

            // Define jet kinematics
            // float jet_rawpt = rawpt[j];
            // float jet_eta = jteta[j];
            // float jet_phi = jtphi[j];

            // // Apply JEC
            // JEC.SetJetPT(jet_rawpt);
            // JEC.SetJetEta(jet_eta);
            // JEC.SetJetPhi(jet_phi);
            // float jet_pt_corr = JEC.GetCorrectedPT();
            // // resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
            // double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr);                                      // Jet weight (specially for MC)
            // jet_weight = jet_weight * get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr, do_jet_smearing, 0.663); // Jet smearing (For systematics)

            // if (jet_eta > jet_eta_min_cut && jet_eta < jet_eta_max_cut)
            // {                                                                                                                                                               // Jet eta cut
            //     find_leading_subleading(jet_pt_corr, jet_eta, jet_phi, leadrecojet_pt, leadrecojet_eta, leadrecojet_phi, sublrecojet_pt, sublrecojet_eta, sublrecojet_phi); // Find leading and subleading jets
            //     hist_reco_jet_weighted_nocut->Fill(jet_pt_corr, event_weight * jet_weight);                                                                                 // Fill histogram without any cut
            //     if (jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
            //     {   // Jet pT cut
            //         // Fill reco jet QA histograms
            //         double x4D_reco_jet[4] = {jet_rawpt, jet_eta, jet_phi, (double)multcentbin};
            //         hist_reco_jet->Fill(x4D_reco_jet);
            //         hist_reco_jet_weighted->Fill(x4D_reco_jet, event_weight * jet_weight);
            //         double x4D_reco_jet_corr[4] = {jet_pt_corr, jet_eta, jet_phi, (double)multcentbin};
            //         hist_reco_jet_corr->Fill(x4D_reco_jet_corr);
            //         hist_reco_jet_corr_weighted->Fill(x4D_reco_jet_corr, event_weight * jet_weight);
            //         // Fill reco jet vectors
            //         TVector3 GoodJets;
            //         GoodJets.SetPtEtaPhi(jet_pt_corr, jet_eta, jet_phi);
            //         jets_reco.push_back(GoodJets);
            //         jet_w_reco.push_back(jet_weight);
            //     }
            // }
        } // for (Int_t iJet{0}; iJet<nref; iJet++) {

        if (isGoodJetLow) hCentVsMultPartFlowJet[0]->Fill(mult, centrality);
        if (isGoodJetMid) hCentVsMultPartFlowJet[1]->Fill(mult, centrality);
        if (isGoodJetHi) hCentVsMultPartFlowJet[2]->Fill(mult, centrality);
        
    } // for (Long64_t iEvent{0}; iEvent < nEventsToRead; iEvent++)
    

    // Write to output file
    TFile *oFile = new TFile(oFileName, "recreate");
    hVz->Write();
    hMult->Write();
    hHiBin->Write();
    hWeight->Write();
    hPtHat->Write();
    hCentrality->Write();
    hCentralityVsMultiplicity->Write();
    for (Int_t i{0}; i<3; i++) {
        hCentVsMultPartFlowJet[i]->Write();
    }
    hRecoJetPt->Write();

    oFile->Close();
    
    std::cout << "Finish" << std::endl;

    return 0;
}
