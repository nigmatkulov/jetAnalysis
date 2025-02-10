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

//_________________
ForestAODReader::ForestAODReader() : fEvent{nullptr}, fInFileName{nullptr}, fEvents2Read{0}, fEventsProcessed{0},
    fIsMc{false}, fCorrectCentMC{false}, fUseHltBranch{kTRUE}, fUseSkimmingBranch{kTRUE}, 
    fUseRecoJetBranch{kTRUE}, 
    fUseTrackBranch{false}, fUseGenTrackBranch{false},
    fHltTree{nullptr}, fSkimTree{nullptr}, fEventTree{nullptr},
    fRecoJetTree{nullptr}, fTrkTree{nullptr}, fGenTrkTree{nullptr},
    fRecoJetTreeName{"akCs4PFJetAnalyzer"},
    fJEC{nullptr}, fJECFiles{}, fJECPath{}, fJEU{nullptr}, fJEUInputFileName{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{false}, 
    fFixJetArrays{false}, fEventCut{nullptr}, fJetCut{nullptr},
    fRecoJet2GenJetId{}, fGenJet2RecoJet{}, 
    fUseExtraJECforAk4Cs{false}, fJECScaleCorr{nullptr}, fUseJEU{0},
    fUseJERSystematics{0}, fAlphaJER{0.0415552}, fBetaJER{0.960013},
    fJERSmearFunc{nullptr}, fRndm{nullptr},
    fVerbose{false} {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::ForestAODReader()" << std::endl;
    }
    // Initialize many variables
    fRndm = new TRandom3(0);
    clearVariables();
    setJERSystParams();
}

//_________________
ForestAODReader::ForestAODReader(const char* inputStream,  
                                 const bool& useHltBranch, const bool& useSkimmingBranch, 
                                 const bool& useRecoJetBranch, 
                                 const bool& useTrackBranch, const bool& useGenTrackBranch, 
                                 const bool& isMc) : 
    fEvent{nullptr}, fInFileName{inputStream}, fEvents2Read{0}, 
    fEventsProcessed{0}, fIsMc{isMc}, fCorrectCentMC{false},
    fUseHltBranch{useHltBranch}, fUseSkimmingBranch{useSkimmingBranch}, 
    fUseRecoJetBranch{useRecoJetBranch}, 
    fUseTrackBranch{useTrackBranch}, fUseGenTrackBranch{useGenTrackBranch},
    fRecoJetTreeName{"akCs4PFJetAnalyzer"},
    fJEC{nullptr}, fJECFiles{}, fJECPath{}, fJEU{nullptr}, fJEUInputFileName{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{false}, 
    fFixJetArrays{false}, fEventCut{nullptr}, fJetCut{nullptr},
    fJECScaleCorr{nullptr}, fUseJEU{0}, fUseJERSystematics{0}, 
    fAlphaJER{0.0415552}, fBetaJER{0.960013}, fJERSmearFunc{nullptr}, 
    fVerbose{false} {
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
    if (fRecoJetTree) delete fRecoJetTree;
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
double ForestAODReader::retrieveResolutionFactor(const double& eta) {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::retrieveResolutionFactor\n";
    }

    double val{1.};
    double res{0.};

    // Search for the bin index
    for (int i{0}; i<fJerEtaLow.size(); i++) {
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
    } // for (int i{0}; i<fJerEtaLow.size(); i++)

    res = TMath::Sqrt( TMath::Max(val * val - 1., 0.) );

    if ( fVerbose ) {
        std::cout << "eta: " << eta << "JER val: " << val << " Resolution factor: " << res << std::endl;
        std::cout << "\t[DONE]\n";
    }
    return res;
}

//________________
double ForestAODReader::extraJERCorr(const double &ptCorr, const double &eta) {
    // This factor applied to the reco jet pT after JEC to match data.
    // By default JEC is not fully cover data/MC JEC difference.
    // Extra correction should be applied to MC only.
    if ( fVerbose ) {
        std::cout << "ForestAODReader::extraJERCorr\n";
    }

    double res = retrieveResolutionFactor(eta);
    double sigmaSmear{0.};
    if ( ptCorr <= 30.) {
        sigmaSmear = res * fJERSmearFunc->Eval( 31. );
    }
    else if ( ptCorr >= 800 ) {
        sigmaSmear = res * fJERSmearFunc->Eval( 799. );
    }
    else {
        sigmaSmear = res * fJERSmearFunc->Eval( ptCorr );
    }

    double extraCorr = fRndm->Gaus( 1., sigmaSmear );

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

    fNRecoJets = {0};
    fNGenJets = {0};
    fNTracks = {0};

    fHLT_HIAK4CaloJet60_v1 = {0};
    fHLT_HIAK4CaloJet80_v1 = {0};
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
    for (short i{0}; i<20000; i++) {

        // Jet variables
        if (i<10000) {
            fRecoJetPt[i] = {0.f};
            fRecoJetEta[i] = {0.f};
            fRecoJetPhi[i] = {0.f};
            fRecoJetWTAEta[i] = {0.f};
            fRecoJetWTAPhi[i] = {0.f};
            fRecoJetTrackMax[i] = {0.f};
            fRefJetPt[i] = {0.f};
            fRefJetEta[i] = {0.f};
            fRefJetPhi[i] = {0.f};
            fRefJetWTAEta[i] = {0.f};
            fRefJetWTAPhi[i] = {0.f};
            fRefJetPartonFlavor[i] = {-999};
            fRefJetPartonFlavorForB[i] = {-99};
            fGenJetPt[i] = {0.f};
            fGenJetEta[i] = {0.f};
            fGenJetPhi[i] = {0.f};
            fGenJetWTAEta[i] = {0.f};
            fGenJetWTAPhi[i] = {0.f};
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
        fTrackHighPurity[i] = {false};
    } // for (short i{0}; i<9999; i++)

    fGenTrackPt.clear();
    fGenTrackEta.clear();
    fGenTrackPhi.clear();
    fGenTrackCharge.clear();
    fGenTrackPid.clear();
    fGenTrackSube.clear();

    if (fIsMc) {
        fRecoJet2GenJetId.clear();
        fGenJet2RecoJet.clear();
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//_________________
int ForestAODReader::init() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::init()" << std::endl;
    }
    int status = 0;
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

    // Need to add path in front of the file name
    std::vector< std::string > tmp;
    for (unsigned int i{0}; i<fJECFiles.size(); i++) {
        tmp.push_back( Form( "%s/aux_files/%s_%i/JEC/%s", 
                             fJECPath.Data(), fCollidingSystem.Data(),
                             fCollidingEnergyGeV, fJECFiles.at(i).c_str() ) );
    }
        
    fJECFiles.clear();
    fJECFiles = tmp;

    std::cout << "JEC files added: " << std::endl;
    for (unsigned int i{0}; i<fJECFiles.size(); i++) {
        std::cout << i << fJECFiles.at(i) << std::endl;
    }
	
	fJEC = new JetCorrector( fJECFiles );

    if ( fUseExtraJECforAk4Cs ) {
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
double ForestAODReader::evalCentralityWeight(const double& x) {
    double weight{1.};
    double p0{4.363352};
    double p1{-8.957467e-02};
    double p2{7.301890e-04};
    double p3{-2.885492e-06};
    double p4{4.741175e-09};
    double p5{0.};
    //double p5{-1.407975e-09};

    weight = p0 + p1 * x + p2 * TMath::Power(x, 2) +
             p3 * TMath::Power(x, 3) + p4 * TMath::Power(x, 4) +
             p5 * TMath::Power(x, 5);

    return weight;
}

//_________________
void ForestAODReader::finish() {
    // Nothing to do here
}

//_________________
int ForestAODReader::setupChains() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupChains()";
    }

    // Setup chains (0-good, 1-bad)
    int returnStatus = 1;

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
    // Use particle flow jet branch
    if ( fUseRecoJetBranch ) {
        fRecoJetTree = new TChain( Form( "%s/t", fRecoJetTreeName.Data() ) );
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
            if ( fUseRecoJetBranch ) fRecoJetTree->Add( input.Data() );
            if ( fUseTrackBranch ) fTrkTree->Add( input.Data() );
            if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( input.Data() );

            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
            Long64_t fEvents2Read2 = fRecoJetTree->GetEntries();
            std::cout << Form("Total number of events to read2: %lld\n", fEvents2Read2 );
        }
        // Assuming that list of files is provided instead of a single file
        else {
            std::ifstream inputStream( input.Data() );

            if ( !inputStream ) std::cout << Form( "ERROR: Cannot open file list: %s\n", input.Data() );
            int nFiles = 0;
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
                        if ( fUseRecoJetBranch ) fRecoJetTree->Add( file.c_str() );
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
    if ( fUseRecoJetBranch ) fRecoJetTree->SetBranchStatus("*", 0);
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
        fHltTree->SetBranchStatus("HLT_HIAK4CaloJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1", 1);
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
        fHltTree->SetBranchAddress("HLT_HIAK4CaloJet60_v1", &fHLT_HIAK4CaloJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1", &fHLT_HIAK4CaloJet80_v1);
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
    if ( fUseRecoJetBranch ) {

        fRecoJetTree->SetBranchStatus("nref", 1);
        fRecoJetTree->SetBranchStatus("rawpt", 1);
        fRecoJetTree->SetBranchStatus("trackMax", 1);
        fRecoJetTree->SetBranchStatus("jteta", 1);
        fRecoJetTree->SetBranchStatus("jtphi", 1);
        fRecoJetTree->SetBranchStatus("WTAeta", 1);
        fRecoJetTree->SetBranchStatus("WTAphi", 1);

        fRecoJetTree->SetBranchAddress("nref", &fNRecoJets);
        fRecoJetTree->SetBranchAddress("rawpt", &fRecoJetPt);
        fRecoJetTree->SetBranchAddress("trackMax", &fRecoJetTrackMax);
        fRecoJetTree->SetBranchAddress("jteta", &fRecoJetEta);
        fRecoJetTree->SetBranchAddress("jtphi", &fRecoJetPhi);
        fRecoJetTree->SetBranchAddress("WTAeta", &fRecoJetWTAEta);
        fRecoJetTree->SetBranchAddress("WTAphi", &fRecoJetWTAPhi);

        fRecoJetTree->SetBranchStatus("jtPfNHF", 1);
        fRecoJetTree->SetBranchStatus("jtPfNEF", 1);
        fRecoJetTree->SetBranchStatus("jtPfCHF", 1);
        fRecoJetTree->SetBranchStatus("jtPfMUF", 1);
        fRecoJetTree->SetBranchStatus("jtPfCEF", 1);
        fRecoJetTree->SetBranchStatus("jtPfCHM", 1);
        fRecoJetTree->SetBranchStatus("jtPfCEM", 1);
        fRecoJetTree->SetBranchStatus("jtPfNHM", 1);
        fRecoJetTree->SetBranchStatus("jtPfNEM", 1);
        fRecoJetTree->SetBranchStatus("jtPfMUM", 1);

      	fRecoJetTree->SetBranchAddress("jtPfNHF", &fRecoJtPfNHF);    
 	    fRecoJetTree->SetBranchAddress("jtPfNEF", &fRecoJtPfNEF);
      	fRecoJetTree->SetBranchAddress("jtPfCHF", &fRecoJtPfCHF);
      	fRecoJetTree->SetBranchAddress("jtPfMUF", &fRecoJtPfMUF);
      	fRecoJetTree->SetBranchAddress("jtPfCEF", &fRecoJtPfCEF);
      	fRecoJetTree->SetBranchAddress("jtPfCHM", &fRecoJtPfCHM);
      	fRecoJetTree->SetBranchAddress("jtPfCEM", &fRecoJtPfCEM);
      	fRecoJetTree->SetBranchAddress("jtPfNHM", &fRecoJtPfNHM);
      	fRecoJetTree->SetBranchAddress("jtPfNEM", &fRecoJtPfNEM);
  	    fRecoJetTree->SetBranchAddress("jtPfMUM", &fRecoJtPfMUM);

        // Gen jet quantities
        if ( fIsMc ) {
            fRecoJetTree->SetBranchStatus("ngen", 1);
            fRecoJetTree->SetBranchStatus("genpt", 1);
            fRecoJetTree->SetBranchStatus("geneta", 1);
            fRecoJetTree->SetBranchStatus("genphi", 1);
            fRecoJetTree->SetBranchStatus("WTAgeneta", 1);
            fRecoJetTree->SetBranchStatus("WTAgenphi", 1);
            fRecoJetTree->SetBranchAddress("ngen", &fNGenJets);
            fRecoJetTree->SetBranchAddress("genpt", &fGenJetPt);
            fRecoJetTree->SetBranchAddress("geneta", &fGenJetEta);
            fRecoJetTree->SetBranchAddress("genphi", &fGenJetPhi);
            fRecoJetTree->SetBranchAddress("WTAgeneta", &fGenJetWTAEta);
            fRecoJetTree->SetBranchAddress("WTAgenphi", &fGenJetWTAPhi);
        }

        // Jet-matching quantities
        if ( fIsMc ) {
            fRecoJetTree->SetBranchStatus("refpt", 1);
            fRecoJetTree->SetBranchStatus("refeta", 1);
            fRecoJetTree->SetBranchStatus("refphi", 1);
            fRecoJetTree->SetBranchStatus("refWTAeta", 1);
            fRecoJetTree->SetBranchStatus("refWTAphi", 1);
            fRecoJetTree->SetBranchStatus("refparton_flavor", 1);
            fRecoJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fRecoJetTree->SetBranchAddress("refpt", &fRefJetPt);
            fRecoJetTree->SetBranchAddress("refeta", &fRefJetEta);
            fRecoJetTree->SetBranchAddress("refphi", &fRefJetPhi);
            fRecoJetTree->SetBranchAddress("refWTAeta", &fRefJetWTAEta);
            fRecoJetTree->SetBranchAddress("refWTAphi", &fRefJetWTAPhi);
            fRecoJetTree->SetBranchAddress("refparton_flavor", &fRefJetPartonFlavor);
            fRecoJetTree->SetBranchAddress("refparton_flavorForB", &fRefJetPartonFlavorForB);
        }

    } // if ( fUseRecoJetBranch )

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
        fRecoJet2GenJetId.clear();
        fGenJet2RecoJet.clear();
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
    if (fUseRecoJetBranch) fRecoJetTree->GetEntry(fEventsProcessed);
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

    if (fUseRecoJetBranch) {

        // Loop over reconstructed jets
        for (int iRecoJet{0}; iRecoJet<fNRecoJets; iRecoJet++) {

            if ( fNGenJets <= 0 ) {
                fRecoJet2GenJetId.push_back(-1);
                continue;
            }

            // Must have a gen-matched jet
            if ( fRefJetPt[iRecoJet] < 0) {
                fRecoJet2GenJetId.push_back(-1);
                continue;
            }

            for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++) {
                // Skip Ref and Gen jets that do not match on pT within computational precision
                //std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fGenJetPt[iGenJet] - fRefJetPt[iRecoJet]) << std::endl;

                //std::cout << "iGen: " << iGenJet << " genPt: " << fGenJetPt[iGenJet] << std::endl;
                if ( TMath::Abs(fGenJetPt[iGenJet] - fRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON ) {
                    // If it is the last one
                    if ( iGenJet == (fNGenJets - 1) ) fRecoJet2GenJetId.push_back(-1);
                    continue;
                }
                else {
                    fRefJetEta[iRecoJet] = fGenJetEta[iGenJet];
                    fRefJetPhi[iRecoJet] = fGenJetPhi[iGenJet];
                    fRefJetWTAEta[iRecoJet] = fGenJetWTAEta[iGenJet];
                    fRefJetWTAPhi[iRecoJet] = fGenJetWTAPhi[iGenJet];
                    fRecoJet2GenJetId.push_back(iGenJet);
                    break;
                }

            }
        } //for (int iRecoJet=0; iRecoJet<fNRecoJets; iRecoJet++)

        // Fill the corresponding index in the reco vector and fill the gen
        for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++) {
            std::vector<int>::iterator it=std::find(fRecoJet2GenJetId.begin(), fRecoJet2GenJetId.end(), iGenJet);
            if (it != fRecoJet2GenJetId.end()) {
                fGenJet2RecoJet.push_back( std::distance(fRecoJet2GenJetId.begin(), it) );
            }
            else {
                fGenJet2RecoJet.push_back(-1);
            }
        }
    } // if (fUseRecoJetBranch)

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

    int nBadRecoJets{0};

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
    float centW{1.f};
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

        fEvent->trigAndSkim()->setHLT_HIAK4CaloJet60_v1(fHLT_HIAK4CaloJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4CaloJet80_v1(fHLT_HIAK4CaloJet80_v1);
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
    if ( fUseRecoJetBranch ) {

        if ( fVerbose ) {
            std::cout << "Use PF branch \n";
        }

        // Loop over generated jets
        if ( fIsMc && !fEvent->isGenJetCollectionFilled() ) {

            if ( fVerbose ) {
                std::cout << "nGenRecoJets: " << fNGenJets << std::endl;
            }

            if ( fVerbose ) {
                std::cout << "Filling GenJets: " << fNGenJets << std::endl;
            }
            for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++) {
                GenJet *jet = new GenJet{};
                jet->setPt( fGenJetPt[iGenJet] );
                jet->setEta( fGenJetEta[iGenJet] );
                jet->setPhi( fGenJetPhi[iGenJet] );
                jet->setWTAEta( fGenJetWTAEta[iGenJet] );
                jet->setWTAPhi( fGenJetWTAPhi[iGenJet] );
                jet->setFlavor( fRefJetPartonFlavor[fGenJet2RecoJet.at(iGenJet)] );
                jet->setFlavorForB( fRefJetPartonFlavorForB[fGenJet2RecoJet.at(iGenJet)] );
                jet->setPtWeight( 1. );
                if ( fVerbose ) {
                    jet->print();
                }
                
                fEvent->genJetCollection()->push_back( jet );
            } // for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++)
            if ( fVerbose) {
                std::cout << "GenJetCollection filled" << std::endl;
            }

            // Projection from filling the collection several times
            fEvent->setGenJetCollectionIsFilled();
        } // if ( fIsMc )
        

        // Loop over reconstructed jets
        
        if ( fVerbose ) {
            std::cout << "nRecoJets: " << fNRecoJets << std::endl;
        }

        for (int iJet{0}; iJet<fNRecoJets; iJet++) {

            // Create a new jet instance
            RecoJet *jet = new RecoJet{};

            if ( fIsMc ) {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if ( fRecoJetPt[iJet] > fPtHat ) {
                    nBadRecoJets++;
                }

                // Add index of the matched GenJet
                jet->setGenJetId( fRecoJet2GenJetId.at(iJet) );
            } // if ( fIsMc )

            // Reco
            jet->setPt( fRecoJetPt[iJet] );
            jet->setEta( fRecoJetEta[iJet] );
            jet->setPhi( fRecoJetPhi[iJet] );
            jet->setWTAEta( fRecoJetWTAEta[iJet] );
            jet->setWTAPhi( fRecoJetWTAPhi[iJet] );
            jet->setRawPt( fRecoJetPt[iJet] );
            jet->setTrackMaxPt( fRecoJetTrackMax[iJet] );
            if ( fJEC ) {
                fJEC->SetJetPT( fRecoJetPt[iJet] );
                fJEC->SetJetEta( fRecoJetEta[iJet] );
                fJEC->SetJetPhi( fRecoJetPhi[iJet] );
                double pTcorr = fJEC->GetCorrectedPT();
                if ( fVerbose ) {
                    std::cout << "pTCorr: " << pTcorr << std::endl; 
                }
                if ( fUseExtraJECforAk4Cs ) {
                    if ( fJECScaleCorr ){
                        pTcorr *= fJECScaleCorr->Eval( pTcorr );
                    }
                    else {
                        std::cerr << "No extra correction for ak4cs exists!" << std::endl;
                    }
                }

                if ( fVerbose && fUseExtraJECforAk4Cs ) {
                    std::cout << "pTCorr after extra ak4cs correction for the jetType: " << pTcorr << std::endl; 
                }

                // To check JEC fUseJECSystematics should be outside [-1, 1] range
                if ( fIsMc && ( TMath::Abs( fUseJERSystematics ) <= 1 ) ) {
                    // pTcorr *= extraJERCorr( pTcorr, fRecoJetEta[iJet]);
                    if ( jet->hasMatching() ) {
                        pTcorr *= extraJERCorr( fEvent->genJetCollection()->at( fRecoJet2GenJetId.at(iJet) )->pt(), 
                                                fRecoJetEta[iJet]);
                    }
                }

                // JEU correction for the real data for systematic uncertainty calculation
                if ( fUseJEU !=0 && fJEU && !fIsMc ) {
                    fJEU->SetJetPT( pTcorr );
			        fJEU->SetJetEta( fRecoJetEta[iJet] );
			        fJEU->SetJetPhi( fRecoJetPhi[iJet] );
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
                if ( fVerbose ) {
                    std::cout << "No JEC available" << std::endl;
                }
                jet->setPtJECCorr( -999.f );
            }
            jet->setJtPfNHF( fRecoJtPfNHF[iJet] );
            jet->setJtPfNEF( fRecoJtPfNEF[iJet] );
            jet->setJtPfCHF( fRecoJtPfCHF[iJet] );
            jet->setJtPfMUF( fRecoJtPfMUF[iJet] );
            jet->setJtPfCEF( fRecoJtPfCEF[iJet] );
            jet->setJtPfCHM( fRecoJtPfCHM[iJet] );
            jet->setJtPfCEM( fRecoJtPfCEM[iJet] );
            jet->setJtPfNHM( fRecoJtPfNHM[iJet] );
            jet->setJtPfNEM( fRecoJtPfNEM[iJet] );
            jet->setJtPfMUM( fRecoJtPfMUM[iJet] );

            if ( fVerbose ) {
                jet->print();
            }

            // Check fronе-loaded cut
            if ( fJetCut && !fJetCut->pass(jet) ) {
                delete jet;
                continue;
            }

            fEvent->recoJetCollection()->push_back( jet );
        } // for (int iJet{0}; iJet<fNRecoJets; iJet++)
    } // if ( fUseRecoJetBranch )

    if ( fEventCut && !fEventCut->pass(fEvent) ) {
        delete fEvent;
        fEvent = nullptr;
    }

    return fEvent;
}