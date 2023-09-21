#ifndef JETANALYSIS_H
#define JETANALYSIS_H

// Load ROOT libraries
#include "TObject.h"
#include "TString.h"
#include "TChain.h"
#include "Rtypes.h"

#include "Jet.h"

/*
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TDatime.h"
#include "THnSparse.h"


// Load C++ libraries
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
*/

//________________
class JetAnalysis : public TObject {
  public:
    /// @brief Constructor
    /// @param inFileName Input filename (name.root) or file with list of ROOT files 
    /// @param oFileName Output ROOT file name
    /// @param nEventsToRead Number of events to process
    JetAnalysis(TString inFileName = "../../../data/HiForestAOD_PbPbMC2018skim_10.root",
                TString oFileName = "oTestSimpleReadForest.root",
                Long64_t nEventsToRead = 500);
    /// @brief Destructor
    virtual ~JetAnalysis();

    /// @brief Initialize variables and functions
    void init();
    /// @brief Process data
    void processData();
    /// @brief Finish analysis
    void finish();

    /// Turn-on event branch to be read
    void useEventBranch()       { fUseEventBranch = kTRUE; }
    /// Turn-on HLT branch to be read
    void useHltBranch()         { fUseHltBranch = kTRUE; }
    /// Turn-on particle flow branch to be read
    void usePartFlowJetBranch() { fUsePartFlowBranch = kTRUE; }
    /// Turn-on calorimeter jet branch to be read
    void useCaloJetBranch()     { fUseCaloJetBranch = kTRUE; }
    /// Turn-on calorimeter jet branch to be read
    void useTrackBranch()       { fUseTrackBranch = kTRUE; }

  private:

    /// Setup input all input
    void setupInput(TString input, TChain *hltChain, TChain *eveChain, TChain *partFlowChain, 
                    TChain *trkChain, Bool_t useMC, TChain *genTrkChain);
    /// Setup chains to be filled
    void setupChains();
    /// Setup branches
    void setupBranches();
    /// Read event and fill objects
    void readEvent();


    /// @brief Input filename (name.root) or file with list of ROOT files 
    TString fInFileName;
    /// @brief Output ROOT file name
    const Char_t *fOutFileName;
    /// @brief Number of events to process from input file(s)
    Long64_t fEvent2Read;

    /// @brief Is file with MC information
    Bool_t fIsMc;

  
    Bool_t fUseEventBranch;
    Bool_t fUseHltBranch;
    Bool_t fUsePartFlowJetBranch;
    Bool_t fUseCaloJetBranch;
    Bool_t fUseTrackBranch;

    /// @brief Chain conaining HLT information (used to friend other trees)
    TChain *fHltTree;
    /// @brief Chain containing event information
    TChain *fEventTree;
    /// @brief Chain containing jets
    TChain *fJetTree;
    /// @brief Chain containing tracks
    TChain *fTrkTree;

  ClassDef(JetAnalysis, 0)
};

#endif // #define JETANALYSIS_H