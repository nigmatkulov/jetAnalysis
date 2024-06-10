/**
 * @file ForestAODReader.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief CMS ForestAOD reader
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "ForestAODReader.h"

// ROOT headers
#include "TFile.h"

// C++ headers
#include <cstddef>
#include <cstring>
#include <fstream>

ClassImp(ForestAODReader)

//_________________
ForestAODReader::ForestAODReader() : fEvent{nullptr}, fInFileName{nullptr}, fEvents2Read{0}, fEventsProcessed{0},
    fIsMc{kFALSE}, fCorrectCentMC{kFALSE}, fUseHltBranch{kTRUE}, fUseSkimmingBranch{kTRUE}, 
    fUsePartFlowJetBranch{kTRUE}, fUseCaloJetBranch{kFALSE}, 
    fUseTrackBranch{kFALSE}, fUseGenTrackBranch{kFALSE},
    fHltTree{nullptr}, fSkimTree{nullptr}, fEventTree{nullptr}, fCaloJetTree{nullptr},
    fPartFlowJetTree{nullptr}, fTrkTree{nullptr}, fGenTrkTree{nullptr},
    fPFTreeName{"akCs4PFJetAnalyzer"}, fCaloTreeName{"fCaloTreeName"},
    fJEC{nullptr}, fJECFiles{}, fJECPath{}, fJEU{nullptr}, fJEUInputFileName{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{kFALSE}, 
    fFixJetArrays{kFALSE}, fEventCut{nullptr}, fJetCut{nullptr},
    fRecoPFJet2GenJetId{}, fGenJet2RecoPFJet{}, 
    fRecoCaloJet2GenJetId{}, fGenJet2RecoCaloJet{},
    fUseExtraJEC{kFALSE}, fJECScaleCorr{nullptr}, fUseJEU{0},
    fUseJERSystematics{0}, fAlphaJER{0.0415552}, fBetaJER{0.960013},
    fJERSmearFunc{nullptr}, fRndm{nullptr},
    fVerbose{kFALSE} {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::ForestAODReader()" << std::endl;
    }
    // Initialize many variables
    fRndm = new TRandom3(0);
    clearVariables();
    setJERSystParams();
}

//_________________
ForestAODReader::ForestAODReader(const Char_t* inputStream,  
                                 const Bool_t& useHltBranch, const Bool_t& useSkimmingBranch, 
                                 const Bool_t& usePFJetBranch, const Bool_t& useCaloJetBranch, 
                                 const Bool_t& useTrackBranch, const Bool_t& useGenTrackBranch, 
                                 const Bool_t& isMc) : 
    fEvent{nullptr}, fInFileName{inputStream}, fEvents2Read{0}, 
    fEventsProcessed{0}, fIsMc{isMc}, fCorrectCentMC{kFALSE},
    fUseHltBranch{useHltBranch}, fUseSkimmingBranch{useSkimmingBranch}, 
    fUsePartFlowJetBranch{usePFJetBranch}, fUseCaloJetBranch{useCaloJetBranch}, 
    fUseTrackBranch{useTrackBranch}, fUseGenTrackBranch{useGenTrackBranch},
    fPFTreeName{"akCs4PFJetAnalyzer"}, fCaloTreeName{"fCaloTreeName"},
    fJEC{nullptr}, fJECFiles{}, fJECPath{}, fJEU{nullptr}, fJEUInputFileName{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{kFALSE}, 
    fFixJetArrays{kFALSE}, fEventCut{nullptr}, fJetCut{nullptr},
    fJECScaleCorr{nullptr}, fUseJEU{0}, fUseJERSystematics{0}, 
    fAlphaJER{0.0415552}, fBetaJER{0.960013}, fJERSmearFunc{nullptr}, 
    fVerbose{kFALSE} {
    // Initialize many variables
    fRndm = new TRandom3(0);
    clearVariables();

    if ( fVerbose ) {
        std::cout << "ForestAODReader::ForestAODReader()" << std::endl;
    }
    setJERSystParams();
}

//_________________
ForestAODReader::~ForestAODReader() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::~ForestAODReader()";
    }
    if (fEvent) delete fEvent;
    if (fHltTree) delete fHltTree;
    if (fSkimTree) delete fSkimTree;
    if (fEventTree) delete fEventTree;
    if (fCaloJetTree) delete fCaloJetTree;
    if (fPartFlowJetTree) delete fPartFlowJetTree;
    if (fTrkTree) delete fTrkTree;
    if (fGenTrkTree) delete fGenTrkTree;
    if (fJEC) delete fJEC;
    if (fJEU) delete fJEU;
    if (fEventCut) delete fEventCut;
    if (fJetCut) delete fJetCut;
    if (fJECScaleCorr) delete fJECScaleCorr;
    if (fJERSmearFunc) delete fJERSmearFunc;
}

//________________
void ForestAODReader::setJERSystParams() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::setJERSystParams";
    }
    fJerEtaLow.clear(); fJerEtaHi.clear(); fJerDef.clear(); fJerLow.clear(); fJerHi.clear();
    fJerEtaLow.push_back(-5.191); fJerEtaHi.push_back(-3.139); fJerDef.push_back(1.1922); fJerLow.push_back(1.0434); fJerHi.push_back(1.341);
    fJerEtaLow.push_back(-3.139); fJerEtaHi.push_back(-2.964); fJerDef.push_back(1.1869); fJerLow.push_back(1.0626); fJerHi.push_back(1.3112);
    fJerEtaLow.push_back(-2.964); fJerEtaHi.push_back(-2.853); fJerDef.push_back(1.7788); fJerLow.push_back(1.578); fJerHi.push_back(1.9796);
    fJerEtaLow.push_back(-2.853); fJerEtaHi.push_back(-2.500); fJerDef.push_back(1.3418); fJerLow.push_back(1.1327); fJerHi.push_back(1.5509);
    fJerEtaLow.push_back(-2.500); fJerEtaHi.push_back(-2.322); fJerDef.push_back(1.2963); fJerLow.push_back(1.0592); fJerHi.push_back(1.5334);
    fJerEtaLow.push_back(-2.322); fJerEtaHi.push_back(-2.043); fJerDef.push_back(1.1512); fJerLow.push_back(1.0372); fJerHi.push_back(1.2652);
    fJerEtaLow.push_back(-2.043); fJerEtaHi.push_back(-1.930); fJerDef.push_back(1.1426); fJerLow.push_back(1.0212); fJerHi.push_back(1.264);
    fJerEtaLow.push_back(-1.930); fJerEtaHi.push_back(-1.740); fJerDef.push_back(1.1000); fJerLow.push_back(0.9921); fJerHi.push_back(1.2079);
    fJerEtaLow.push_back(-1.740); fJerEtaHi.push_back(-1.305); fJerDef.push_back(1.1278); fJerLow.push_back(1.0292); fJerHi.push_back(1.2264);
    fJerEtaLow.push_back(-1.305); fJerEtaHi.push_back(-1.131); fJerDef.push_back(1.1609); fJerLow.push_back(1.0584); fJerHi.push_back(1.2634);
    fJerEtaLow.push_back(-1.131); fJerEtaHi.push_back(-0.783); fJerDef.push_back(1.1464); fJerLow.push_back(1.0832); fJerHi.push_back(1.2096);
    fJerEtaLow.push_back(-0.783); fJerEtaHi.push_back(-0.522); fJerDef.push_back(1.1948); fJerLow.push_back(1.1296); fJerHi.push_back(1.26);
    fJerEtaLow.push_back(-0.522); fJerEtaHi.push_back(0.000); fJerDef.push_back(1.15958); fJerLow.push_back(1.095); fJerHi.push_back(1.224);
    fJerEtaLow.push_back(0.000); fJerEtaHi.push_back(0.522); fJerDef.push_back(1.15958); fJerLow.push_back(1.095); fJerHi.push_back(1.224);
    fJerEtaLow.push_back(0.522); fJerEtaHi.push_back(0.783); fJerDef.push_back(1.1948); fJerLow.push_back(1.1296); fJerHi.push_back(1.26);
    fJerEtaLow.push_back(0.783); fJerEtaHi.push_back(1.131); fJerDef.push_back(1.1464); fJerLow.push_back(1.0832); fJerHi.push_back(1.2096);
    fJerEtaLow.push_back(1.131); fJerEtaHi.push_back(1.305); fJerDef.push_back(1.1609); fJerLow.push_back(1.0584); fJerHi.push_back(1.2634);
    fJerEtaLow.push_back(1.305); fJerEtaHi.push_back(1.740); fJerDef.push_back(1.1278); fJerLow.push_back(1.0292); fJerHi.push_back(1.2264);
    fJerEtaLow.push_back(1.740); fJerEtaHi.push_back(1.930); fJerDef.push_back(1.1000); fJerLow.push_back(0.9921); fJerHi.push_back(1.2079);
    fJerEtaLow.push_back(1.930); fJerEtaHi.push_back(2.043); fJerDef.push_back(1.1426); fJerLow.push_back(1.0212); fJerHi.push_back(1.264);
    fJerEtaLow.push_back(2.043); fJerEtaHi.push_back(2.322); fJerDef.push_back(1.1512); fJerLow.push_back(1.0372); fJerHi.push_back(1.2652);
    fJerEtaLow.push_back(2.322); fJerEtaHi.push_back(2.500); fJerDef.push_back(1.2963); fJerLow.push_back(1.0592); fJerHi.push_back(1.5334);
    fJerEtaLow.push_back(2.500); fJerEtaHi.push_back(2.853); fJerDef.push_back(1.3418); fJerLow.push_back(1.1327); fJerHi.push_back(1.5509);
    fJerEtaLow.push_back(2.853); fJerEtaHi.push_back(2.964); fJerDef.push_back(1.7788); fJerLow.push_back(1.578); fJerHi.push_back(1.9796);
    fJerEtaLow.push_back(2.964); fJerEtaHi.push_back(3.139); fJerDef.push_back(1.1869); fJerLow.push_back(1.0626); fJerHi.push_back(1.3112);
    fJerEtaLow.push_back(3.139); fJerEtaHi.push_back(5.191); fJerDef.push_back(1.1922); fJerLow.push_back(1.0434); fJerHi.push_back(1.341);

    fJERSmearFunc = new TF1("fJERSmearFunc","sqrt( [0] * [0] + [1] * [1] / x )", 30., 800.);
    fJERSmearFunc->SetParameter(0, fAlphaJER);
    fJERSmearFunc->SetParameter(1, fBetaJER);

    if ( fVerbose ) {
        std::cout << "\t[DONE]\n";
    }
}

//________________
Double_t ForestAODReader::retrieveResolutionFactor(const Double_t& eta) {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::retrieveResolutionFactor\n";
    }

    Double_t val{1.};
    Double_t res{0.};

    // Search for the bin index
    for (Int_t i{0}; i<fJerEtaLow.size(); i++) {
        if ( eta>=fJerEtaLow.at(i) && eta<fJerEtaHi.at(i) ) {
            if ( fUseJERSystematics == -1 ) {
                val = fJerLow.at(i);
            }
            else if ( fUseJERSystematics == 0 ) {
                val = fJerDef.at(i);
            }
            else if ( fUseJERSystematics == 1 ) {
                val = fJerHi.at(i);
            }
            else {
                val = {1.};
            }
            break;
        }
    } // for (Int_t i{0}; i<fJerEtaLow.size(); i++)

    res = TMath::Sqrt( TMath::Max(val * val - 1., 0.) );

    if ( fVerbose ) {
        std::cout << "eta: " << eta << "JER val: " << val << " Resolution factor: " << res << std::endl;
        std::cout << "\t[DONE]\n";
    }
    return res;
}

//________________
Double_t ForestAODReader::extraJERCorr(const Double_t &ptCorr, const Double_t &eta) {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::extraJERCorr\n";
    }

    Double_t res = retrieveResolutionFactor(eta);
    Double_t sigmaSmear{0.};
    if ( ptCorr <= 30.) {
        sigmaSmear = res * fJERSmearFunc->Eval( 31. );
    }
    else if ( ptCorr >= 800 ) {
        sigmaSmear = res * fJERSmearFunc->Eval( 799. );
    }
    else {
        sigmaSmear = res * fJERSmearFunc->Eval( ptCorr );
    }

    Double_t extraCorr = fRndm->Gaus( 1., sigmaSmear );

    if ( fVerbose ) {
        std::cout << "Resolution factor: " << res << " sigma: " << sigmaSmear 
                  << " correction factor: " << extraCorr << std::endl;
        std::cout << "\t[DONE]\n";
    }
    return extraCorr;
}

//_________________
void ForestAODReader::clearVariables() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::clearVariables " << std::endl;
    }
    fRunId = {0};
    fEventId = {0};
    fLumi = {0};
    fVertexZ = {-999.f};
    fHiBin = {-1};
    fPtHatWeight = {-1.f};
    fPtHat = {-1.f};

    // bad jets and multiplicity to be added

    fNPFRecoJets = {0};
    fNPFGenJets = {0};
    fNCaloRecoJets = {0};
    fNCaloGenJets = {0};
    fNTracks = {0};

    fHLT_PAAK4CaloJet60_Eta5p1_v3 = {0};
    fHLT_PAAK4CaloJet80_Eta5p1_v3 = {0};
    fHLT_PAAK4CaloJet100_Eta5p1_v3 = {0};
    fHLT_PAAK4PFJet60_Eta5p1_v4 = {0};
    fHLT_PAAK4PFJet80_Eta5p1_v3 = {0};
    fHLT_PAAK4PFJet100_Eta5p1_v3 = {0};
    fHLT_PAAK4PFJet120_Eta5p1_v2 = {0};

    fHLT_HIAK4PFJet15_v1 = {0};
    fHLT_HIAK4PFJet15_v1_Prescl = {0};
    fHLT_HIAK4PFJet30_v1 = {0};
    fHLT_HIAK4PFJet30_v1_Prescl = {0};
    fHLT_HIAK4PFJet40_v1 = {0};
    fHLT_HIAK4PFJet40_v1_Prescl = {0};
    fHLT_HIAK4PFJet60_v1 = {0};
    fHLT_HIAK4PFJet60_v1_Prescl = {0};
    fHLT_HIAK4PFJet80_v1 = {0};
    fHLT_HIAK4PFJet80_v1_Prescl = {0};
    fHLT_HIAK4PFJet120_v1 = {0};
    fHLT_HIAK4PFJet120_v1_Prescl = {0};

    fHLT_HIAK8PFJet15_v1 = {0};
    fHLT_HIAK8PFJet15_v1_Prescl = {0};
    fHLT_HIAK8PFJet25_v1 = {0};
    fHLT_HIAK8PFJet25_v1_Prescl = {0};
    fHLT_HIAK8PFJet40_v1 = {0};
    fHLT_HIAK8PFJet40_v1_Prescl = {0};
    fHLT_HIAK8PFJet60_v1 = {0};
    fHLT_HIAK8PFJet60_v1_Prescl = {0};
    fHLT_HIAK8PFJet80_v1 = {0};
    fHLT_HIAK8PFJet80_v1_Prescl = {0};
    fHLT_HIAK8PFJet140_v1 = {0};
    fHLT_HIAK8PFJet140_v1_Prescl = {0};

    fHLT_HIPFJet25_v1 = {0};
    fHLT_HIPFJet25_v1_Prescl = {0};
    fHLT_HIPFJet140_v1 = {0};
    fHLT_HIPFJet140_v1_Prescl = {0};

    fHLT_HIPuAK4CaloJet80Eta5p1_v1 = {0};
    fHLT_HIPuAK4CaloJet100Eta5p1_v1 = {0};

    fHBHENoiseFilterResultRun2Loose = {0};
    fHBHENoiseFilterResultRun2Tight = {0};
    fHBHEIsoNoiseFilterResult = {0};
    fCollisionEventSelectionAODv2 = {0};
    fPhfCoincFilter2Th4 = {0};
    fPPAprimaryVertexFilter = {0};
    fPBeamScrapingFilter = {0};
    fPprimaryVertexFilter = {0};
    fPVertexFilterCutG = {0};
    fPVertexFilterCutGloose = {0};
    fPVertexFilterCutGtight = {0};
    fPVertexFilterCutE = {0};
    fPVertexFilterCutEandG = {0};
    fPClusterCompatibilityFilter = {0};

    fPhfCoincFilter = {0};
    fPVertexFilterCutdz1p0 = {0};
    fPVertexFilterCutGplus = {0};
    fPVertexFilterCutVtx1 = {0};

    // Loop over jets and tracks
    for (Short_t i{0}; i<9999; i++) {

        // Jet variables
        if (i<100) {
            fPFRecoJetPt[i] = {0.f};
            fPFRecoJetEta[i] = {0.f};
            fPFRecoJetPhi[i] = {0.f};
            fPFRecoJetWTAEta[i] = {0.f};
            fPFRecoJetWTAPhi[i] = {0.f};
            fPFRecoJetTrackMax[i] = {0.f};
            fPFRefJetPt[i] = {0.f};
            fPFRefJetEta[i] = {0.f};
            fPFRefJetPhi[i] = {0.f};
            fPFRefJetWTAEta[i] = {0.f};
            fPFRefJetWTAPhi[i] = {0.f};
            fPFRefJetPartonFlavor[i] = {-999};
            fPFRefJetPartonFlavorForB[i] = {-99};
            fPFGenJetPt[i] = {0.f};
            fPFGenJetEta[i] = {0.f};
            fPFGenJetPhi[i] = {0.f};
            fPFGenJetWTAEta[i] = {0.f};
            fPFGenJetWTAPhi[i] = {0.f};

            fCaloRecoJetPt[i] = {0.f};
            fCaloRecoJetEta[i] = {0.f};
            fCaloRecoJetPhi[i] = {0.f};
            fCaloRecoJetWTAEta[i] = {0.f};
            fCaloRecoJetWTAPhi[i] = {0.f};
            fCaloRecoJetTrackMax[i] = {0.f};
            fCaloRefJetPt[i] = {0.f};
            fCaloRefJetEta[i] = {0.f};
            fCaloRefJetPhi[i] = {0.f};
            fCaloRefJetWTAEta[i] = {0.f};
            fCaloRefJetWTAPhi[i] = {0.f};
            fCaloRefJetPartonFlavor[i] = {-999};
            fCaloRefJetPartonFlavorForB[i] = {-99};
            fCaloGenJetPt[i] = {0.f};
            fCaloGenJetEta[i] = {0.f};
            fCaloGenJetPhi[i] = {0.f};
            fCaloGenJetWTAEta[i] = {0.f};
            fCaloGenJetWTAPhi[i] = {0.f};
        } // if (i<100)

        // Track variables
        fTrackPt[i] = {0.f};
        fTrackEta[i] = {0.f};
        fTrackPhi[i] = {0.f};
        fTrackPtErr[i] = {0.f};
        fTrackDcaXY[i] = {0.f};
        fTrackDcaZ[i] = {0.f};
        fTrackDcaXYErr[i] = {0.f};
        fTrackDcaZErr[i] = {0.f};
        fTrackChi2[i] = {0.f};
        fTrackNDOF[i] = {0};
        fTrackPartFlowEcal[i] = {0.f};
        fTrackPartFlowHcal[i] = {0.f};
        fTrackMVA[i] = {0.f};
        fTrackAlgo[i] = {0};
        fTrackCharge[i] = {0};
        fTrackNHits[i] = {0};
        fTrackNLayers[i] = {0};
        fTrackHighPurity[i] = {kFALSE};
    } // for (Short_t i{0}; i<9999; i++)

    fGenTrackPt.clear();
    fGenTrackEta.clear();
    fGenTrackPhi.clear();
    fGenTrackCharge.clear();
    fGenTrackPid.clear();
    fGenTrackSube.clear();

    if (fIsMc) {
        fRecoPFJet2GenJetId.clear();
        fGenJet2RecoPFJet.clear();

        fRecoCaloJet2GenJetId.clear();
        fGenJet2RecoCaloJet.clear();
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//_________________
Int_t ForestAODReader::init() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::init()" << std::endl;
    }
    Int_t status = 0;
    // Setup chains to read
    status = setupChains();
    // Setup branches to read
    setupBranches();
    // Setup jet energy correction files and pointer
    setupJEC();
    // Setup jet energy uncertainty files and pointer
    setupJEU();
    if ( fIsMc ) {
        if ( TMath::Abs(fUseJERSystematics)<=1 ) {
            setJERSystParams();
        }
    }
    if ( fVerbose ) {
        std::cout << "ForestAODReader::init() is finished " << std::endl;
    }
    return status;
}

//________________
void ForestAODReader::setupJEC() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupJEC()" << std::endl;
    }

    // If no path to the aux_file
    if ( fJECPath.Length() <= 0 ) {
        // Set default values
        std::cout << "[WARNING] Default path to JEC files will be used" << std::endl;
        setPath2JetAnalysis();
    }

    if ( fJECFiles.empty() ) {
        std::cout << "[WARNING] Default JEC file with parameters will be used" << std::endl;
        addJECFile();
    }

    std::vector< std::string > tmp;
    for (UInt_t i{0}; i<fJECFiles.size(); i++) {
        tmp.push_back( Form( "%s/aux_files/%s_%i/JEC/%s", 
                             fJECPath.Data(), fCollidingSystem.Data(),
                             fCollidingEnergyGeV, fJECFiles.at(i).c_str() ) );
    }
        
    fJECFiles.clear();
    fJECFiles = tmp;

    std::cout << "JEC files added: " << std::endl;
    for (UInt_t i{0}; i<fJECFiles.size(); i++) {
        std::cout << i << fJECFiles.at(i) << std::endl;
    }
	
	fJEC = new JetCorrector( fJECFiles );

    if ( fUseExtraJEC ) {
        createExtraJECScaleCorrFunction();
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//________________
void ForestAODReader::setupJEU() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupJEU()" << std::endl;
    }

    // Next part is needed only if JEU correction is applied
    if ( fUseJEU == 0) return;

    // If no path to the aux_file
    if ( fJECPath.Length() <= 0 ) {
        // Set default values
        std::cout << "[WARNING] Default path to JEU files will be used" << std::endl;
        setPath2JetAnalysis();
    }

    // If no correction file is specified
    if ( fJEUInputFileName.Length() <= 0 ) {
        std::cout << "[WARNING] Default JEU file with parameters will be used" << std::endl;
        setJEUFileName();
    }

    TString tmp = Form( "%s/aux_files/%s_%i/JEC/%s", 
                        fJECPath.Data(), fCollidingSystem.Data(),
                        fCollidingEnergyGeV, fJEUInputFileName.Data() );
    fJEUInputFileName = tmp;

    fJEU = new JetUncertainty( fJEUInputFileName.Data() );
    std::cout << "JEU file: " << fJEUInputFileName.Data() << std::endl;

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }   
}

//________________
void ForestAODReader::createExtraJECScaleCorrFunction() {
    fJECScaleCorr = new TF1("JetScaleCorrection","[3] + ([0]-[3]) / ( 1.0 + pow( x/[2],[1] ) )", 30, 800);
    fJECScaleCorr->SetParameters(2.27008e+00, 9.18625e-01, 1.43067e+00, 1.00002e+00);
}

//________________
Float_t ForestAODReader::jetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, 
                                     const Int_t &energy, float jetpt) const {
    Float_t jetptweight = 1.0;

    // JetPtWeightFunction is derived from MC vs data jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *JetPtWeightFunction = new TF1("JetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120 GeV and JECv6
        JetPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09);
        jetptweight = JetPtWeightFunction->Eval(jetpt);
    }
    */
    return jetptweight;
}

//________________
Float_t ForestAODReader::leadJetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t& year, 
                                         const Int_t& energy, const Float_t& leadjetpt) const {
    Float_t leadjetptweight = 1.0;

    // LeadJetPtWeightFunction is derived from MC vs data leading jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *LeadJetPtWeightFunction = new TF1("LeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading jets above 120 GeV and JECv6
        LeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
        leadjetptweight = LeadJetPtWeightFunction->Eval(leadjetpt);
    }
    */
    return leadjetptweight;
}

//________________
Float_t ForestAODReader::subleadJetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, 
                                            const Int_t &energy, const Float_t &subleadjetpt) {
    Float_t subleadjetptweight = 1.0;

    // SubLeadJetPtWeightFunction is derived from MC vs data subleading jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *SubLeadJetPtWeightFunction = new TF1("SubLeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading jets above 120 GeV and JECv6
        SubLeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
        subleadjetptweight = SubLeadJetPtWeightFunction->Eval(subleadjetpt);
    }
    */
    return subleadjetptweight;
}

//________________
Float_t ForestAODReader::jetPtSmeringWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, 
                                            const Int_t &energy, const Float_t &jetpt, 
                                            const Bool_t &dosmearing, const Float_t resolutionfactor) const {
    Float_t jetptsmearweight = 1.0;
    if(!dosmearing) return jetptsmearweight;

    // JetPtSmearingWeightFunction is derived from MC vs data jet pT spectra.
    if(!isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *JetPtSmearingWeightFunction = new TF1("JetPtSmearingWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120 GeV and JECv6
        JetPtSmearingWeightFunction->SetParameters(0.174881, -0.00091979, 3.50064e-06, -6.52541e-09, 4.64199e-12);
        jetptsmearweight = JetPtSmearingWeightFunction->Eval(jetpt);
        jetptsmearweight = jetptsmearweight*resolutionfactor;
    }

    return jetptsmearweight;
}

//________________
Float_t ForestAODReader::trkEtaMixWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, const Int_t &energy, 
                                         const Float_t &trketa, const Bool_t &reco) const {
    float trketamixweight = 1.0;

    // TrkEtaMixWeightFunction is derived from trk eta from signal over trk eta from mixing
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017 && !reco){
        TF1 *TrkEtaMixWeightFunction = new TF1("TrkEtaMixWeightFunction", "pol3", 0.0, 500.0); 
        TrkEtaMixWeightFunction->SetParameters(0.174881, -0.00091979, 3.50064e-06, -6.52541e-09, 4.64199e-12);
        trketamixweight = TrkEtaMixWeightFunction->Eval(jetpt);
    }
    */
    return trketamixweight;
}

//________________
Float_t ForestAODReader::eventWeight(const Bool_t &isMC, const Bool_t &use_centrality, 
                                     const std::string& system, const Int_t &year, const Int_t &energy, 
                                     const Float_t &vz, const Int_t mult, const Float_t &weighttree, 
                                     const Float_t &leadjetpt) const {

    Float_t vzweight = 1.0;
    Float_t multweight = 1.0;
    Float_t evtweight = 1.0;
    Float_t multefficiency = 1.0;
    Float_t jetefficiency = 1.0;		
    Float_t totalweight = 1.0;

    // VzWeightFunction is derived from MC vs data event Vz --> MC only --> vzweight
    // MultCentWeightFunction is derived from MC vs data event multiplicity or centrality --> MC only --> multweight
    // MultTriggerWeightFunction is derived from the turn on plots as function of multiplicity --> RECO only
    // JetTriggerWeightFunction is derived from the turn on plots as function of leading jet pT --> RECO only
    // weighttree is the pthat weight --> MC only 

	if (isMC && !use_centrality && system == "pp" && 
        energy == 5020 && year == 2017) {

		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol6", -15.0, 15.0);
		VzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
		vzweight = VzWeightFunction->Eval(vz);

		TF1 *MultCentWeightFunction = new TF1("MultCentWeightFunction", "pol0", 0.0, 500.0);
		MultCentWeightFunction->SetParameter(0,1.0);
		multweight = MultCentWeightFunction->Eval(mult);

		TF1 *MultTriggerWeightFunction = new TF1("MultTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
		MultTriggerWeightFunction->SetParameter(0,1.0);
		Float_t multtrigweight = 1.0;
		multtrigweight = MultTriggerWeightFunction->Eval(mult);
		multefficiency = 1./multtrigweight;

		TF1 *JetTriggerWeightFunction = new TF1("JetTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
		JetTriggerWeightFunction->SetParameter(0,1.0);
		Float_t jettrigweight = 1.0;
		jettrigweight = JetTriggerWeightFunction->Eval(leadjetpt);
		jetefficiency = 1./jettrigweight;

		evtweight = weighttree;
	}

	totalweight = evtweight * multweight * vzweight * multefficiency * jetefficiency;
	return totalweight;
}

//________________
Double_t ForestAODReader::evalCentralityWeight(const Double_t& x) {
    Double_t weight{1.};
    Double_t p0{4.363352};
    Double_t p1{-8.957467e-02};
    Double_t p2{7.301890e-04};
    Double_t p3{-2.885492e-06};
    Double_t p4{4.741175e-09};
    Double_t p5{0.};
    //Double_t p5{-1.407975e-09};

    weight = p0 + p1 * x + p2 * TMath::Power(x, 2) +
             p3 * TMath::Power(x, 3) + p4 * TMath::Power(x, 4) +
             p5 * TMath::Power(x, 5);

    return weight;
}

//_________________
void ForestAODReader::finish() {

}

//_________________
Int_t ForestAODReader::setupChains() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupChains()";
    }

    // Setup chains (0-good, 1-bad)
    Int_t returnStatus = 1;

    // Setup chains to read

    std::cout << "Setting chains... ";

    // Use event branch
    fEventTree = new TChain("hiEvtAnalyzer/HiTree");

    // Use HLT branch
    if ( fUseHltBranch ) {
        fHltTree = new TChain("hltanalysis/HltTree");
    }
    // Use skimming branch
    if ( fUseSkimmingBranch ) {
        fSkimTree = new TChain("skimanalysis/HltTree");
    }
    // Use calo jet branch
    if ( fUseCaloJetBranch ) {
        fCaloJetTree = new TChain( Form( "%s/t", fCaloTreeName.Data() ) );
    }
    // Use particle flow jet branch
    if ( fUsePartFlowJetBranch ) {
        fPartFlowJetTree = new TChain( Form( "%s/t", fPFTreeName.Data() ) );
    }
    // Use reconstructed track branch
    if ( fUseTrackBranch ) {
        fTrkTree = new TChain("ppTrack/trackTree");
    }
    // Use generated track branch
    if ( fIsMc && fUseGenTrackBranch ) {
        fGenTrkTree = new TChain("HiGenParticleAna/hi");
    }
    std::cout << "\t[DONE]\n";

    // Initialize input file name (should switch to const char* processing later)
    TString input(fInFileName);

    // Check input exists
    if (  input.Length()<= 0 ) {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        returnStatus = 1;
        exit(0);
    }
    // Normail input
    else {
        // If input is a single ROOT file
        if ( input.Index(".root") > 0 ) {
            std::cout << Form( "Adding %s file to chains\n", input.Data() );
            fEventTree->Add( input.Data() );
            if ( fUseHltBranch ) fHltTree->Add( input.Data() );
            if ( fUseSkimmingBranch ) fSkimTree->Add( input.Data() );
            if ( fUseCaloJetBranch ) fCaloJetTree->Add(input.Data() );
            if ( fUsePartFlowJetBranch ) fPartFlowJetTree->Add( input.Data() );
            if ( fUseTrackBranch ) fTrkTree->Add( input.Data() );
            if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( input.Data() );

            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
            Long64_t fEvents2Read2 = fPartFlowJetTree->GetEntries();
            std::cout << Form("Total number of events to read2: %lld\n", fEvents2Read2 );
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
                if ( file.find(".root") != std::string::npos /* && file.find("Forest") != std::string::npos &&
                     file.find("AOD") != std::string::npos */ ) {
                    
                    // Open file
                    TFile* ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if ( ftmp && !ftmp->IsZombie() && ftmp->GetNkeys() ) {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str() );
                        // Adding file to chains
                        fEventTree->Add( file.c_str() );
                        if ( fUseHltBranch ) fHltTree->Add( file.c_str() );
                        if ( fUseSkimmingBranch ) fSkimTree->Add( file.c_str() );
                        if ( fUseCaloJetBranch ) fCaloJetTree->Add( file.c_str() );
                        if ( fUsePartFlowJetBranch ) fPartFlowJetTree->Add( file.c_str() );
                        if ( fUseTrackBranch ) fTrkTree->Add( file.c_str() );
                        if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( file.c_str() );
                        ++nFiles;
                    } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp) {
                        ftmp->Close();
                    } //if (ftmp)
                } //if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos && file.find("AOD") != std::string::npos )
            } //while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
        } // else {   if file list
        returnStatus = 0;
    } // else {   if normal input

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
    return returnStatus;
}

//_________________
void ForestAODReader::setupBranches() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupChains()";
    }

    // Disable all branches - this is important while reading big files
    fEventTree->SetBranchStatus("*", 0);
    if ( fUseHltBranch ) fHltTree->SetBranchStatus("*", 0);
    if ( fUseSkimmingBranch ) fSkimTree->SetBranchStatus("*", 0);
    if ( fUseCaloJetBranch) fCaloJetTree->SetBranchStatus("*", 0);
    if ( fUsePartFlowJetBranch ) fPartFlowJetTree->SetBranchStatus("*", 0);
    if ( fUseTrackBranch ) fTrkTree->SetBranchStatus("*", 0);
    if ( fUseGenTrackBranch && fIsMc ) fGenTrkTree->SetBranchStatus("*", 0);


    // enable branches of interest -> see definition of each variables above

    // Event quantities
    fEventTree->SetBranchStatus("run", 1);
    fEventTree->SetBranchStatus("evt", 1);
    fEventTree->SetBranchStatus("lumi", 1);
    fEventTree->SetBranchStatus("vz", 1);
    fEventTree->SetBranchStatus("hiBin", 1); //centrality only for PbPb and XeXe
    fEventTree->SetBranchAddress("run", &fRunId);
    fEventTree->SetBranchAddress("evt", &fEventId);
    fEventTree->SetBranchAddress("lumi", &fLumi);
    fEventTree->SetBranchAddress("vz", &fVertexZ);
    fEventTree->SetBranchAddress("hiBin", &fHiBin); 

    if( fIsMc ) {
        fEventTree->SetBranchStatus("weight", 1);
        fEventTree->SetBranchStatus("pthat", 1); 
        fEventTree->SetBranchAddress("weight", &fPtHatWeight);
        fEventTree->SetBranchAddress("pthat", &fPtHat);
    }

    if ( fUseHltBranch ) {

        // Status
        fHltTree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet120_Eta5p1_v2", 1);

        fHltTree->SetBranchStatus("HLT_HIAK4PFJet15_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet15_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet30_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet30_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet40_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet40_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet60_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet80_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet120_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet120_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIAK8PFJet15_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet15_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet25_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet25_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet40_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet40_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet60_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet80_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet140_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet140_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIPFJet25_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet25_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet140_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet140_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1", 1);

        // Address
        fHltTree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3", &fHLT_PAAK4CaloJet60_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3", &fHLT_PAAK4CaloJet80_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3", &fHLT_PAAK4CaloJet100_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4", &fHLT_PAAK4PFJet60_Eta5p1_v4);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3", &fHLT_PAAK4PFJet80_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3", &fHLT_PAAK4PFJet100_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet120_Eta5p1_v2", &fHLT_PAAK4PFJet120_Eta5p1_v2);

        fHltTree->SetBranchAddress("HLT_HIAK4PFJet15_v1", &fHLT_HIAK4PFJet15_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet15_v1_Prescl", &fHLT_HIAK4PFJet15_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet30_v1", &fHLT_HIAK4PFJet30_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet30_v1_Prescl", &fHLT_HIAK4PFJet30_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet40_v1", &fHLT_HIAK4PFJet40_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet40_v1_Prescl", &fHLT_HIAK4PFJet40_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet60_v1", &fHLT_HIAK4PFJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet60_v1_Prescl", &fHLT_HIAK4PFJet60_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet80_v1", &fHLT_HIAK4PFJet80_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet80_v1_Prescl", &fHLT_HIAK4PFJet80_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet120_v1", &fHLT_HIAK4PFJet120_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet120_v1_Prescl", &fHLT_HIAK4PFJet120_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIAK8PFJet15_v1", &fHLT_HIAK8PFJet15_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet15_v1_Prescl", &fHLT_HIAK8PFJet15_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet25_v1", &fHLT_HIAK8PFJet25_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet25_v1_Prescl", &fHLT_HIAK8PFJet25_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet40_v1", &fHLT_HIAK8PFJet40_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet40_v1_Prescl", &fHLT_HIAK8PFJet40_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet60_v1", &fHLT_HIAK8PFJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet60_v1_Prescl", &fHLT_HIAK8PFJet60_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet80_v1", &fHLT_HIAK8PFJet80_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet80_v1_Prescl", &fHLT_HIAK8PFJet80_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet140_v1", &fHLT_HIAK8PFJet140_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet140_v1_Prescl", &fHLT_HIAK8PFJet140_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIPFJet25_v1", &fHLT_HIPFJet25_v1);
        fHltTree->SetBranchAddress("HLT_HIPFJet25_v1_Prescl", &fHLT_HIPFJet25_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIPFJet140_v1", &fHLT_HIPFJet140_v1);
        fHltTree->SetBranchAddress("HLT_HIPFJet140_v1_Prescl", &fHLT_HIPFJet140_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1", &fHLT_HIPuAK4CaloJet80Eta5p1_v1);
        fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1", &fHLT_HIPuAK4CaloJet100Eta5p1_v1);
    } // if ( fUseHltBranch )

    // Skimming quantities
    if ( fUseSkimmingBranch ) {
        // Status
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Tight", 1);
        fSkimTree->SetBranchStatus("HBHEIsoNoiseFilterResult", 1);
        fSkimTree->SetBranchStatus("collisionEventSelectionAODv2", 1);
        fSkimTree->SetBranchStatus("phfCoincFilter2Th4", 1);
        fSkimTree->SetBranchStatus("pPAprimaryVertexFilter", 1);
        fSkimTree->SetBranchStatus("pBeamScrapingFilter", 1);
        fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutG", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGloose", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGtight", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutE", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutEandG", 1);
        fSkimTree->SetBranchStatus("pclusterCompatibilityFilter", 1);

        fSkimTree->SetBranchStatus("phfCoincFilter", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutdz1p0", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGplus", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutVtx1", 1);


        // Address
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &fHBHENoiseFilterResultRun2Loose);
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &fHBHENoiseFilterResultRun2Tight);
        fSkimTree->SetBranchAddress("HBHEIsoNoiseFilterResult", &fHBHEIsoNoiseFilterResult);
        fSkimTree->SetBranchAddress("collisionEventSelectionAODv2", &fCollisionEventSelectionAODv2);
        fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fPhfCoincFilter2Th4);
        fSkimTree->SetBranchAddress("pPAprimaryVertexFilter", &fPPAprimaryVertexFilter);
        fSkimTree->SetBranchAddress("pBeamScrapingFilter", &fPBeamScrapingFilter);
        fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPprimaryVertexFilter);
        fSkimTree->SetBranchAddress("pVertexFilterCutG", &fPVertexFilterCutG);
        fSkimTree->SetBranchAddress("pVertexFilterCutGloose", &fPVertexFilterCutGloose);
        fSkimTree->SetBranchAddress("pVertexFilterCutGtight", &fPVertexFilterCutGtight);
        fSkimTree->SetBranchAddress("pVertexFilterCutE", &fPVertexFilterCutE);
        fSkimTree->SetBranchAddress("pVertexFilterCutEandG", &fPVertexFilterCutEandG);
        fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fPClusterCompatibilityFilter);

        fSkimTree->SetBranchAddress("phfCoincFilter", &fPhfCoincFilter);
        fSkimTree->SetBranchAddress("pVertexFilterCutdz1p0", &fPVertexFilterCutdz1p0);
        fSkimTree->SetBranchAddress("pVertexFilterCutGplus", &fPVertexFilterCutGplus);
        fSkimTree->SetBranchAddress("pVertexFilterCutVtx1", &fPVertexFilterCutVtx1);
    } // if ( fUseSkimmingBranch )

    // Jet quantities
    if ( fUsePartFlowJetBranch ) {

        fPartFlowJetTree->SetBranchStatus("nref", 1);
        fPartFlowJetTree->SetBranchStatus("rawpt", 1);
        fPartFlowJetTree->SetBranchStatus("trackMax", 1);
        fPartFlowJetTree->SetBranchStatus("jteta", 1);
        fPartFlowJetTree->SetBranchStatus("jtphi", 1);
        fPartFlowJetTree->SetBranchStatus("WTAeta", 1);
        fPartFlowJetTree->SetBranchStatus("WTAphi", 1);

        fPartFlowJetTree->SetBranchAddress("nref", &fNPFRecoJets);
        fPartFlowJetTree->SetBranchAddress("rawpt", &fPFRecoJetPt);
        fPartFlowJetTree->SetBranchAddress("trackMax", &fPFRecoJetTrackMax);
        fPartFlowJetTree->SetBranchAddress("jteta", &fPFRecoJetEta);
        fPartFlowJetTree->SetBranchAddress("jtphi", &fPFRecoJetPhi);
        fPartFlowJetTree->SetBranchAddress("WTAeta", &fPFRecoJetWTAEta);
        fPartFlowJetTree->SetBranchAddress("WTAphi", &fPFRecoJetWTAPhi);

        fPartFlowJetTree->SetBranchStatus("jtPfNHF", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfNEF", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfCHF", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfMUF", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfCEF", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfCHM", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfCEM", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfNHM", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfNEM", 1);
        fPartFlowJetTree->SetBranchStatus("jtPfMUM", 1);

      	fPartFlowJetTree->SetBranchAddress("jtPfNHF", &fPFRecoJtPfNHF);    
 	    fPartFlowJetTree->SetBranchAddress("jtPfNEF", &fPFRecoJtPfNEF);
      	fPartFlowJetTree->SetBranchAddress("jtPfCHF", &fPFRecoJtPfCHF);
      	fPartFlowJetTree->SetBranchAddress("jtPfMUF", &fPFRecoJtPfMUF);
      	fPartFlowJetTree->SetBranchAddress("jtPfCEF", &fPFRecoJtPfCEF);
      	fPartFlowJetTree->SetBranchAddress("jtPfCHM", &fPFRecoJtPfCHM);
      	fPartFlowJetTree->SetBranchAddress("jtPfCEM", &fPFRecoJtPfCEM);
      	fPartFlowJetTree->SetBranchAddress("jtPfNHM", &fPFRecoJtPfNHM);
      	fPartFlowJetTree->SetBranchAddress("jtPfNEM", &fPFRecoJtPfNEM);
  	    fPartFlowJetTree->SetBranchAddress("jtPfMUM", &fPFRecoJtPfMUM);

        // Gen jet quantities
        if ( fIsMc ) {
            fPartFlowJetTree->SetBranchStatus("ngen", 1);
            fPartFlowJetTree->SetBranchStatus("genpt", 1);
            fPartFlowJetTree->SetBranchStatus("geneta", 1);
            fPartFlowJetTree->SetBranchStatus("genphi", 1);
            fPartFlowJetTree->SetBranchStatus("WTAgeneta", 1);
            fPartFlowJetTree->SetBranchStatus("WTAgenphi", 1);
            fPartFlowJetTree->SetBranchAddress("ngen", &fNPFGenJets);
            fPartFlowJetTree->SetBranchAddress("genpt", &fPFGenJetPt);
            fPartFlowJetTree->SetBranchAddress("geneta", &fPFGenJetEta);
            fPartFlowJetTree->SetBranchAddress("genphi", &fPFGenJetPhi);
            fPartFlowJetTree->SetBranchAddress("WTAgeneta", &fPFGenJetWTAEta);
            fPartFlowJetTree->SetBranchAddress("WTAgenphi", &fPFGenJetWTAPhi);
        }

        // Jet-matching quantities
        if ( fIsMc ) {
            fPartFlowJetTree->SetBranchStatus("refpt", 1);
            fPartFlowJetTree->SetBranchStatus("refeta", 1);
            fPartFlowJetTree->SetBranchStatus("refphi", 1);
            fPartFlowJetTree->SetBranchStatus("refWTAeta", 1);
            fPartFlowJetTree->SetBranchStatus("refWTAphi", 1);
            fPartFlowJetTree->SetBranchStatus("refparton_flavor", 1);
            fPartFlowJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fPartFlowJetTree->SetBranchAddress("refpt", &fPFRefJetPt);
            fPartFlowJetTree->SetBranchAddress("refeta", &fPFRefJetEta);
            fPartFlowJetTree->SetBranchAddress("refphi", &fPFRefJetPhi);
            fPartFlowJetTree->SetBranchAddress("refWTAeta", &fPFRefJetWTAEta);
            fPartFlowJetTree->SetBranchAddress("refWTAphi", &fPFRefJetWTAPhi);
            fPartFlowJetTree->SetBranchAddress("refparton_flavor", &fPFRefJetPartonFlavor);
            fPartFlowJetTree->SetBranchAddress("refparton_flavorForB", &fPFRefJetPartonFlavorForB);
        }

    } // if ( fUsePartFlowJetBranch )

    if ( fUseCaloJetBranch ) {
        fCaloJetTree->SetBranchStatus("nref", 1);
        fCaloJetTree->SetBranchStatus("rawpt", 1);
        fCaloJetTree->SetBranchStatus("trackMax", 1);
        fCaloJetTree->SetBranchStatus("jteta", 1);
        fCaloJetTree->SetBranchStatus("jtphi", 1);
        fCaloJetTree->SetBranchStatus("WTAeta", 1);
        fCaloJetTree->SetBranchStatus("WTAphi", 1);

        fCaloJetTree->SetBranchAddress("nref", &fNCaloRecoJets);
        fCaloJetTree->SetBranchAddress("rawpt", &fCaloRecoJetPt);
        fCaloJetTree->SetBranchAddress("trackMax", &fCaloRecoJetTrackMax);
        fCaloJetTree->SetBranchAddress("jteta", &fCaloRecoJetEta);
        fCaloJetTree->SetBranchAddress("jtphi", &fCaloRecoJetPhi);
        fCaloJetTree->SetBranchAddress("WTAeta", &fCaloRecoJetWTAEta);
        fCaloJetTree->SetBranchAddress("WTAphi", &fCaloRecoJetWTAPhi);

        // Gen jet quantities
        if ( fIsMc ) {
            fCaloJetTree->SetBranchStatus("ngen", 1);
            fCaloJetTree->SetBranchStatus("genpt", 1);
            fCaloJetTree->SetBranchStatus("geneta", 1);
            fCaloJetTree->SetBranchStatus("genphi", 1);
            fCaloJetTree->SetBranchStatus("WTAgeneta", 1);
            fCaloJetTree->SetBranchStatus("WTAgenphi", 1);

            fCaloJetTree->SetBranchAddress("ngen", &fNCaloGenJets);
            fCaloJetTree->SetBranchAddress("genpt", &fCaloGenJetPt);
            fCaloJetTree->SetBranchAddress("geneta", &fCaloGenJetEta);
            fCaloJetTree->SetBranchAddress("genphi", &fCaloGenJetPhi);
            fCaloJetTree->SetBranchAddress("WTAgeneta", &fCaloGenJetWTAEta);
            fCaloJetTree->SetBranchAddress("WTAgenphi", &fCaloGenJetWTAPhi);
        }

        // Jet-matching quantities
        if ( fIsMc ) {
            fCaloJetTree->SetBranchStatus("refpt", 1);
            fCaloJetTree->SetBranchStatus("refeta", 1);
            fCaloJetTree->SetBranchStatus("refphi", 1);
            fCaloJetTree->SetBranchAddress("refWTAeta", &fCaloRefJetWTAEta);
            fCaloJetTree->SetBranchAddress("refWTAphi", &fCaloRefJetWTAPhi);
            fCaloJetTree->SetBranchStatus("refparton_flavor", 1);
            fCaloJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fCaloJetTree->SetBranchAddress("refpt", &fCaloRefJetPt);
            fCaloJetTree->SetBranchAddress("refeta", &fCaloRefJetEta);
            fCaloJetTree->SetBranchAddress("refphi", &fCaloRefJetPhi);
            fCaloJetTree->SetBranchAddress("refWTAeta", &fCaloRefJetWTAEta);
            fCaloJetTree->SetBranchAddress("refWTAphi", &fCaloRefJetWTAPhi);
            fCaloJetTree->SetBranchAddress("refparton_flavor", &fCaloRefJetPartonFlavor);
            fCaloJetTree->SetBranchAddress("refparton_flavorForB", &fCaloRefJetPartonFlavorForB);
        }
    } // if ( fUsePartFlowJetBranch )

    // Track quantities
    if ( fUseTrackBranch ) {
        fTrkTree->SetBranchStatus("nTrk", 1);
        fTrkTree->SetBranchStatus("trkPt", 1);
        fTrkTree->SetBranchStatus("trkEta", 1);
        fTrkTree->SetBranchStatus("trkPhi", 1);
        fTrkTree->SetBranchStatus("trkPtError", 1);
        fTrkTree->SetBranchStatus("trkDxy1", 1);
        fTrkTree->SetBranchStatus("trkDxyError1", 1);
        fTrkTree->SetBranchStatus("trkDz1", 1);
        fTrkTree->SetBranchStatus("trkDzError1", 1);
        fTrkTree->SetBranchStatus("trkChi2", 1);
        fTrkTree->SetBranchStatus("trkNdof", 1);
        fTrkTree->SetBranchStatus("trkCharge", 1);
        fTrkTree->SetBranchStatus("trkNHit", 1);
        fTrkTree->SetBranchStatus("trkNlayer", 1);
        fTrkTree->SetBranchStatus("highPurity", 1);
        fTrkTree->SetBranchStatus("pfEcal", 1);
        fTrkTree->SetBranchStatus("pfHcal", 1);
        fTrkTree->SetBranchStatus("trkMVA", 1);
        fTrkTree->SetBranchStatus("trkAlgo", 1);

        fTrkTree->SetBranchAddress("nTrk", &fNTracks);
        fTrkTree->SetBranchAddress("trkPt", &fTrackPt);
        fTrkTree->SetBranchAddress("trkEta", &fTrackEta);
        fTrkTree->SetBranchAddress("trkPhi", &fTrackPhi);
        fTrkTree->SetBranchAddress("trkPtError", &fTrackPtErr);
        fTrkTree->SetBranchAddress("trkDxy1", &fTrackDcaXY);
        fTrkTree->SetBranchAddress("trkDxyError1", &fTrackDcaXYErr);
        fTrkTree->SetBranchAddress("trkDz1", &fTrackDcaZ);
        fTrkTree->SetBranchAddress("trkDzError1", &fTrackDcaZErr);
        fTrkTree->SetBranchAddress("trkChi2", &fTrackChi2);
        fTrkTree->SetBranchAddress("trkNdof", &fTrackNDOF);
        fTrkTree->SetBranchAddress("trkCharge", &fTrackCharge);
        fTrkTree->SetBranchAddress("trkNHit", &fTrackNHits);
        fTrkTree->SetBranchAddress("trkNlayer", &fTrackNLayers);
        fTrkTree->SetBranchAddress("highPurity", &fTrackHighPurity);
        fTrkTree->SetBranchAddress("pfEcal", &fTrackPartFlowEcal);
        fTrkTree->SetBranchAddress("pfHcal", &fTrackPartFlowHcal);
        fTrkTree->SetBranchAddress("trkMVA", &fTrackMVA);
        fTrkTree->SetBranchAddress("trkAlgo", &fTrackAlgo);
    } // if ( fUseTrackBranch ) 

    // Gen particle quantities
    if( fIsMc && fUseGenTrackBranch ) {
        fGenTrkTree->SetBranchStatus("pt", 1);
        fGenTrkTree->SetBranchStatus("eta", 1);
        fGenTrkTree->SetBranchStatus("phi", 1);
        fGenTrkTree->SetBranchStatus("chg", 1);
        fGenTrkTree->SetBranchStatus("pdg", 1);
        fGenTrkTree->SetBranchStatus("sube", 1);

        fGenTrkTree->SetBranchAddress("pt", &fGenTrackPt);
        fGenTrkTree->SetBranchAddress("eta", &fGenTrackEta);
        fGenTrkTree->SetBranchAddress("phi", &fGenTrackPhi);
        fGenTrkTree->SetBranchAddress("chg", &fGenTrackCharge);
        fGenTrkTree->SetBranchAddress("pdg", &fGenTrackPid);
        fGenTrkTree->SetBranchAddress("sube", &fGenTrackSube);
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//_________________
void ForestAODReader::report() {
    std::cout << "ForestAODReader::reporting" << std::endl;
}

//_________________
void ForestAODReader::readEvent() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::readEvent()\n";
    }

    if (fIsMc) {
        fRecoPFJet2GenJetId.clear();
        fGenJet2RecoPFJet.clear();

        fRecoCaloJet2GenJetId.clear();
        fGenJet2RecoCaloJet.clear();
    }

    // Or one can call the clearVariables() function (will take more time)
    if (fUseGenTrackBranch && fIsMc) {
        fGenTrackPt.clear();
        fGenTrackEta.clear();
        fGenTrackPhi.clear();
        fGenTrackCharge.clear();
        fGenTrackPid.clear();
        fGenTrackSube.clear();
    }

    if ( fEventsProcessed >= fEvents2Read ) { 
        std::cerr << "ForestAODReader::readEvent() out of entry numbers\n"; 
        fReaderStatus = 2; // End of input stream
    }
    fEventTree->GetEntry( fEventsProcessed );
    if (fUseHltBranch) fHltTree->GetEntry(fEventsProcessed);
    if (fUseSkimmingBranch) fSkimTree->GetEntry(fEventsProcessed);
    if (fUseCaloJetBranch) fCaloJetTree->GetEntry(fEventsProcessed);
    if (fUsePartFlowJetBranch) fPartFlowJetTree->GetEntry(fEventsProcessed);
    if (fUseTrackBranch) fTrkTree->GetEntry(fEventsProcessed);
    if (fUseGenTrackBranch) fGenTrkTree->GetEntry(fEventsProcessed);
    fEventsProcessed++;

    if ( fVerbose ) {
        std::cout << "Events processed: " << fEventsProcessed << std::endl;
        std::cout << "ForestAODReader::readEvent() \t[DONE]" << std::endl;
    }
}

//________________
void ForestAODReader::fixIndices() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::fixIndices()";
    }

    if (fUsePartFlowJetBranch) {

        // Loop over reconstructed jets
        for (Int_t iRecoJet{0}; iRecoJet<fNPFRecoJets; iRecoJet++) {

            if ( fNPFGenJets <= 0 ) {
                fRecoPFJet2GenJetId.push_back(-1);
                continue;
            }

            // Must have a gen-matched jet
            if ( fPFRefJetPt[iRecoJet] < 0) {
                fRecoPFJet2GenJetId.push_back(-1);
                continue;
            }

            for (Int_t iGenJet{0}; iGenJet<fNPFGenJets; iGenJet++) {
                // Skip Ref and Gen jets that do not match on pT within computational precision
                //std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fPFGenJetPt[iGenJet] - fPFRefJetPt[iRecoJet]) << std::endl;

                //std::cout << "iGen: " << iGenJet << " genPt: " << fPFGenJetPt[iGenJet] << std::endl;
                if ( TMath::Abs(fPFGenJetPt[iGenJet] - fPFRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON ) {
                    // If it is the last one
                    if ( iGenJet == (fNPFGenJets - 1) ) fRecoPFJet2GenJetId.push_back(-1);
                    continue;
                }
                else {
                    fPFRefJetEta[iRecoJet] = fPFGenJetEta[iGenJet];
                    fPFRefJetPhi[iRecoJet] = fPFGenJetPhi[iGenJet];
                    fPFRefJetWTAEta[iRecoJet] = fPFGenJetWTAEta[iGenJet];
                    fPFRefJetWTAPhi[iRecoJet] = fPFGenJetWTAPhi[iGenJet];
                    fRecoPFJet2GenJetId.push_back(iGenJet);
                    break;
                }

            }
        } //for (Int_t iRecoJet=0; iRecoJet<fNPFRecoJets; iRecoJet++)

        // Fill the corresponding index in the reco vector and fill the gen
        for (Int_t iGenJet{0}; iGenJet<fNPFGenJets; iGenJet++) {
            std::vector<Int_t>::iterator it=std::find(fRecoPFJet2GenJetId.begin(), fRecoPFJet2GenJetId.end(), iGenJet);
            if (it != fRecoPFJet2GenJetId.end()) {
                fGenJet2RecoPFJet.push_back( std::distance(fRecoPFJet2GenJetId.begin(), it) );
            }
            else {
                fGenJet2RecoPFJet.push_back(-1);
            }
        }
    } // if (fUsePartFlowJetBranch)

    if (fUseCaloJetBranch) {
        //std::cout << "Will fix Calo jets. nRecoJets: " << fNCaloRecoJets << " nGenJet: " << fNCaloGenJets << std::endl;

        // Loop over reconstructed jets
        for (Int_t iRecoJet{0}; iRecoJet<fNCaloRecoJets; iRecoJet++) {

            if ( fNCaloGenJets <= 0 ) {
                fRecoCaloJet2GenJetId.push_back(-1);
                continue;    
            }

            // Must have a gen-matched jet
            if ( fCaloRefJetPt[iRecoJet] < 0) {
                fRecoCaloJet2GenJetId.push_back(-1);
                continue;
            }
            for (Int_t iGenJet{0}; iGenJet<fNCaloGenJets; iGenJet++) {
                // Skip Ref and Gen jets that do not match on pT within computational precision
                //std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fCaloGenJetPt[iGenJet] - fCaloRefJetPt[iRecoJet]) << std::endl;
                if ( TMath::Abs(fCaloGenJetPt[iGenJet] - fCaloRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON ) {
                    // In case of the last generated jet
                    if ( iGenJet == (fNCaloGenJets-1) ) fRecoCaloJet2GenJetId.push_back(-1);;
                    continue;
                }
                else {
                    fCaloRefJetEta[iRecoJet] = fCaloGenJetEta[iGenJet];
                    fCaloRefJetPhi[iRecoJet] = fCaloGenJetPhi[iGenJet];
                    fCaloRefJetWTAEta[iRecoJet] = fCaloGenJetWTAEta[iGenJet];
                    fCaloRefJetWTAPhi[iRecoJet] = fCaloGenJetWTAPhi[iGenJet];
                    fRecoCaloJet2GenJetId.push_back(iGenJet);
                    break;
                }
                
            }
        } //for (Int_t iRecoJet=0; iRecoJet<fNCaloRecoJets; iRecoJet++)

        // Fill the corresponding index in the reco vector and fill the gen
        for (Int_t iGenJet{0}; iGenJet<fNCaloGenJets; iGenJet++) {
            std::vector<Int_t>::iterator it=std::find(fRecoCaloJet2GenJetId.begin(), fRecoCaloJet2GenJetId.end(), iGenJet);
            if (it != fRecoCaloJet2GenJetId.end()) {
                fGenJet2RecoCaloJet.push_back( std::distance(fRecoCaloJet2GenJetId.begin(), it) );
            }
            else {
                fGenJet2RecoCaloJet.push_back(-1);
            }
        }
    } // if (fUseCaloJetBranch)

    if ( fVerbose ) {
        std::cout << "\t[DONE]\n";
    }
}

//_________________
Event* ForestAODReader::returnEvent() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::returnEvent() \n";
    }

    //std::cout << "ForestAODReader::returnEvent" << std::endl;
    readEvent();

    Int_t nBadPFJets{0};
    Int_t nBadCaloJets{0};

    if (fFixJetArrays && fIsMc) {
        fixIndices();
    }

    fEvent = new Event();

    // Remove UPC bins
    if ( fIsMc && fCorrectCentMC && fHiBin<10) {
        delete fEvent;
        fEvent = nullptr;
        return fEvent;
    }

    fEvent->setRunId( fRunId );
    fEvent->setEventId( fEventId );
    fEvent->setLumi( fLumi );
    fEvent->setVz( fVertexZ );
    Float_t centW{1.f};
    if ( fIsMc && fCorrectCentMC) {
        // To handle centrality weight
        fEvent->setHiBin( fHiBin - 10 );
        centW = evalCentralityWeight( fHiBin - 10 );
    }
    else if ( fIsMc ) {
        fEvent->setHiBin( fHiBin );
    }
    else {
        fEvent->setHiBin( fHiBin );
    }
    fEvent->setCentralityWeight( centW );
    
    if ( fIsMc ) {
        fEvent->setPtHat( fPtHat );
        fEvent->setPtHatWeight( fPtHatWeight );
    }
    else {
        fEvent->setPtHat( 1. );
        fEvent->setPtHatWeight( 1. );        
    }

    // Fill HLT branch
    if ( fUseHltBranch ) {

        fEvent->trigAndSkim()->setHLT_PAAK4CaloJet60_Eta5p1_v3(fHLT_PAAK4CaloJet60_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4CaloJet80_Eta5p1_v3(fHLT_PAAK4CaloJet80_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4CaloJet100_Eta5p1_v3(fHLT_PAAK4CaloJet100_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet60_Eta5p1_v4(fHLT_PAAK4PFJet60_Eta5p1_v4);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet80_Eta5p1_v3(fHLT_PAAK4PFJet80_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet100_Eta5p1_v3(fHLT_PAAK4PFJet100_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet120_Eta5p1_v2(fHLT_PAAK4PFJet120_Eta5p1_v2);

        fEvent->trigAndSkim()->setHLT_HIAK4PFJet15_v1(fHLT_HIAK4PFJet15_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet15_v1_Prescl(fHLT_HIAK4PFJet15_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet30_v1(fHLT_HIAK4PFJet30_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet30_v1_Prescl(fHLT_HIAK4PFJet30_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet40_v1(fHLT_HIAK4PFJet40_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet40_v1_Prescl(fHLT_HIAK4PFJet40_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet60_v1(fHLT_HIAK4PFJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet60_v1_Prescl(fHLT_HIAK4PFJet60_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet80_v1(fHLT_HIAK4PFJet80_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet80_v1_Prescl(fHLT_HIAK4PFJet80_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet120_v1(fHLT_HIAK4PFJet120_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet120_v1_Prescl(fHLT_HIAK4PFJet120_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIAK8PFJet15_v1(fHLT_HIAK8PFJet15_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet15_v1_Prescl(fHLT_HIAK8PFJet15_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet25_v1(fHLT_HIAK8PFJet25_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet25_v1_Prescl(fHLT_HIAK8PFJet25_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet40_v1(fHLT_HIAK8PFJet40_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet40_v1_Prescl(fHLT_HIAK8PFJet40_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet60_v1(fHLT_HIAK8PFJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet60_v1_Prescl(fHLT_HIAK8PFJet60_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet80_v1(fHLT_HIAK8PFJet80_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet80_v1_Prescl(fHLT_HIAK8PFJet80_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet140_v1(fHLT_HIAK8PFJet140_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet140_v1_Prescl(fHLT_HIAK8PFJet140_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIPFJet25_v1(fHLT_HIPFJet25_v1);
        fEvent->trigAndSkim()->setHLT_HIPFJet25_v1_Prescl(fHLT_HIPFJet25_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIPFJet140_v1(fHLT_HIPFJet140_v1);
        fEvent->trigAndSkim()->setHLT_HIPFJet140_v1_Prescl(fHLT_HIPFJet140_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIPuAK4CaloJet80Eta5p1_v1(fHLT_HIPuAK4CaloJet80Eta5p1_v1);
        fEvent->trigAndSkim()->setHLT_HIPuAK4CaloJet100Eta5p1_v1(fHLT_HIPuAK4CaloJet100Eta5p1_v1);
    }

    // Fill skimming branch
    if ( fUseSkimmingBranch ) {
        fEvent->trigAndSkim()->setHBHENoiseFilterResultRun2Loose(fHBHENoiseFilterResultRun2Loose);
        fEvent->trigAndSkim()->setHBHENoiseFilterResultRun2Tight(fHBHENoiseFilterResultRun2Tight);
        fEvent->trigAndSkim()->setHBHEIsoNoiseFilterResult(fHBHEIsoNoiseFilterResult);
        fEvent->trigAndSkim()->setCollisionEventSelectionAODv2(fCollisionEventSelectionAODv2);
        fEvent->trigAndSkim()->setPhfCoincFilter2Th4(fPhfCoincFilter2Th4);
        fEvent->trigAndSkim()->setPPAprimaryVertexFilter(fPPAprimaryVertexFilter);
        fEvent->trigAndSkim()->setPBeamScrapingFilter(fPBeamScrapingFilter);
        fEvent->trigAndSkim()->setPprimaryVertexFilter(fPprimaryVertexFilter);
        fEvent->trigAndSkim()->setPVertexFilterCutG(fPVertexFilterCutG);
        fEvent->trigAndSkim()->setPVertexFilterCutGloose(fPVertexFilterCutGloose);
        fEvent->trigAndSkim()->setPVertexFilterCutGtight(fPVertexFilterCutGtight);
        fEvent->trigAndSkim()->setPVertexFilterCutE(fPVertexFilterCutE);
        fEvent->trigAndSkim()->setPVertexFilterCutEandG(fPVertexFilterCutEandG);
        fEvent->trigAndSkim()->setPClusterCompatibilityFilter(fPClusterCompatibilityFilter);

        fEvent->trigAndSkim()->setPhfCoincFilter(fPhfCoincFilter);
        fEvent->trigAndSkim()->setPVertexFilterCutdz1p0(fPVertexFilterCutdz1p0);
        fEvent->trigAndSkim()->setPVertexFilterCutGplus(fPVertexFilterCutGplus);
        fEvent->trigAndSkim()->setPVertexFilterCutVtx1(fPVertexFilterCutVtx1);
    }

    //fEvent->print();
    
    // Create particle flow jet instances
    if ( fUsePartFlowJetBranch ) {

        if ( fVerbose ) {
            std::cout << "Use PF branch \n";
        }

        // Loop over generated jets
        if ( fIsMc && !fEvent->isGenJetCollectionFilled() ) {

            if ( fVerbose ) {
                std::cout << "nGenPFJets: " << fNPFGenJets << std::endl;
            }

            for (Int_t iGenJet{0}; iGenJet<fNPFGenJets; iGenJet++) {
                GenJet *jet = new GenJet{};
                jet->setPt( fPFGenJetPt[iGenJet] );
                jet->setEta( fPFGenJetEta[iGenJet] );
                jet->setPhi( fPFGenJetPhi[iGenJet] );
                jet->setWTAEta( fPFGenJetWTAEta[iGenJet] );
                jet->setWTAPhi( fPFGenJetWTAPhi[iGenJet] );
                jet->setFlavor( fPFRefJetPartonFlavor[fGenJet2RecoPFJet.at(iGenJet)] );
                jet->setFlavorForB( fPFRefJetPartonFlavorForB[fGenJet2RecoPFJet.at(iGenJet)] );
                jet->setPtWeight( jetPtWeight(fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, 
                                              fCollidingEnergyGeV, fPFGenJetPt[iGenJet]) );
                fEvent->genJetCollection()->push_back( jet );
            } // for (Int_t iGenJet{0}; iGenJet<fNPFGenJets; iGenJet++)

            // Projection from filling the collection several times
            fEvent->setGenJetCollectionIsFilled();
        } // if ( fIsMc )
        

        // Loop over reconstructed jets
        
        if ( fVerbose ) {
            std::cout << "nRecoPFJets: " << fNPFRecoJets << std::endl;
        }

        for (Int_t iJet{0}; iJet<fNPFRecoJets; iJet++) {

            // Create a new jet instance
            RecoJet *jet = new RecoJet{};

            if ( fIsMc ) {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if ( fPFRecoJetPt[iJet] > fPtHat ) {
                    nBadPFJets++;
                }

                // Add index of the matched GenJet
                jet->setGenJetId( fRecoPFJet2GenJetId.at(iJet) );
            } // if ( fIsMc )

            // Reco
            jet->setPt( fPFRecoJetPt[iJet] );
            jet->setEta( fPFRecoJetEta[iJet] );
            jet->setPhi( fPFRecoJetPhi[iJet] );
            jet->setWTAEta( fPFRecoJetWTAEta[iJet] );
            jet->setWTAPhi( fPFRecoJetWTAPhi[iJet] );
            jet->setRawPt( fPFRecoJetPt[iJet] );
            jet->setTrackMaxPt( fPFRecoJetTrackMax[iJet] );
            if ( fJEC ) {
                fJEC->SetJetPT( fPFRecoJetPt[iJet] );
                fJEC->SetJetEta( fPFRecoJetEta[iJet] );
                fJEC->SetJetPhi( fPFRecoJetPhi[iJet] );
                double pTcorr = fJEC->GetCorrectedPT();
                if ( fVerbose ) {
                    std::cout << "pTCorr: " << pTcorr << std::endl; 
                }
                if ( fUseExtraJEC ) {
                    if ( fJECScaleCorr ){
                        pTcorr *= fJECScaleCorr->Eval( pTcorr );
                    }
                    else {
                        std::cerr << "No extra correction exists!" << std::endl;
                    }
                }

                if ( fVerbose && fUseExtraJEC ) {
                    std::cout << "pTCorr after axtra correction for the jetType: " << pTcorr << std::endl; 
                }

                if ( fIsMc && ( TMath::Abs( fUseJERSystematics ) <= 1 ) ) {
                    // pTcorr *= extraJERCorr( pTcorr, fPFRecoJetEta[iJet]);

                    pTcorr *= extraJERCorr( fEvent->genJetCollection()->at( fRecoPFJet2GenJetId.at(iJet) )->pt(), 
                                            fPFRecoJetEta[iJet]);
                }

                // JEU correction for the real data for systematic uncertainty calculation
                if ( fUseJEU !=0 && fJEU && !fIsMc ) {
                    fJEU->SetJetPT( pTcorr );
			        fJEU->SetJetEta( fPFRecoJetEta[iJet] );
			        fJEU->SetJetPhi( fPFRecoJetPhi[iJet] );
                    if ( fUseJEU > 0 ) {
                        pTcorr *= (1. + fJEU->GetUncertainty().first);
                    }
                    else {
                        pTcorr *= (1. - fJEU->GetUncertainty().second);
                    }

                    if ( fVerbose ) {
                        std::cout << "pTCorr after JEU: " << pTcorr << std::endl; 
                    }
                }
                jet->setPtJECCorr( pTcorr );
            }
            else { // If no JEC available
                jet->setPtJECCorr( -999.f );
            }
            jet->setJtPfNHF( fPFRecoJtPfNHF[iJet] );
            jet->setJtPfNEF( fPFRecoJtPfNEF[iJet] );
            jet->setJtPfCHF( fPFRecoJtPfCHF[iJet] );
            jet->setJtPfMUF( fPFRecoJtPfMUF[iJet] );
            jet->setJtPfCEF( fPFRecoJtPfCEF[iJet] );
            jet->setJtPfCHM( fPFRecoJtPfCHM[iJet] );
            jet->setJtPfCEM( fPFRecoJtPfCEM[iJet] );
            jet->setJtPfNHM( fPFRecoJtPfNHM[iJet] );
            jet->setJtPfNEM( fPFRecoJtPfNEM[iJet] );
            jet->setJtPfMUM( fPFRecoJtPfMUM[iJet] );

            if ( fVerbose ) {
                jet->print();
            }

            // Check fronе-loaded cut
            if ( fJetCut && !fJetCut->pass(jet) ) {
                delete jet;
                continue;
            }

            fEvent->pfJetCollection()->push_back( jet );
        } // for (Int_t iJet{0}; iJet<fNPFRecoJets; iJet++)
    } // if ( fUsePartFlowJetBranch )


    // Create jet instances
    if ( fUseCaloJetBranch ) {

        // Loop over reconstructed jets
        for (Int_t iJet{0}; iJet<fNCaloRecoJets; iJet++) {

            // Create a new jet instance
            RecoJet *jet = new RecoJet{};

                if ( fIsMc ) {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if ( fCaloRecoJetPt[iJet] > fPtHat ) {
                    nBadPFJets++;
                }

                // Add index of the matched GenJet
                jet->setGenJetId( fRecoCaloJet2GenJetId.at(iJet) );
            } // if ( fIsMc )


            // Reco
            jet->setPt( fCaloRecoJetPt[iJet] );
            jet->setEta( fCaloRecoJetEta[iJet] );
            jet->setPhi( fCaloRecoJetPhi[iJet] );
            jet->setWTAEta( fCaloRecoJetWTAEta[iJet] );
            jet->setWTAPhi( fCaloRecoJetWTAPhi[iJet] );
            if ( fJEC ) {
                fJEC->SetJetPT( fCaloRecoJetPt[iJet] );
                fJEC->SetJetEta( fCaloRecoJetEta[iJet] );
                fJEC->SetJetPhi( fCaloRecoJetPhi[iJet] );
                double pTcorr = fJEC->GetCorrectedPT();
                //std::cout << "pTCorr: " << pTcorr << std::endl;
                if ( fUseExtraJEC ) {
                    pTcorr *= fJECScaleCorr->Eval( pTcorr );
                }
                jet->setPtJECCorr( pTcorr );
                //std::cout << "pTCorr: " << jet->recoJetPtJECCorr() << std::endl; 
            }
            else { // If no JEC available
                jet->setPtJECCorr( -999.f );
            } 

            // Check fron-loaded cut
            if ( fJetCut && !fJetCut->pass(jet) ) {
                delete jet;
                continue;
            }

            fEvent->caloJetCollection()->push_back( jet );
        } // for (Int_t iJet{0}; iJet<fNCaloRecoJets; iJet++)
    } // if ( fUseCaloJetBranch )

    if ( fEventCut && !fEventCut->pass(fEvent) ) {
        delete fEvent;
        fEvent = nullptr;
    }

    return fEvent;
}