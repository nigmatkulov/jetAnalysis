/**
 * @file TriggerAndSkim.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Trigger and skim properties of the event
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef TriggerAndSkim_h
#define TriggerAndSkim_h

// ROOT headers
#include "TObject.h"
#include "Rtypes.h"

//________________
class TriggerAndSkim : public TObject {
  public:
    /// @brief Constructor
    TriggerAndSkim();
    /// @brief Destructor
    virtual ~TriggerAndSkim() {/* empty */ }

    //
    // Setters
    //

    // Set HLT flags

    void setHLT_PAAK4CaloJet60_Eta5p1_v3(const Int_t& b) { fHLT_PAAK4CaloJet60_Eta5p1_v3 = (Char_t)b; }
    void setHLT_PAAK4CaloJet80_Eta5p1_v3(const Int_t& b) { fHLT_PAAK4CaloJet80_Eta5p1_v3 = (Char_t)b; }
    void setHLT_PAAK4CaloJet100_Eta5p1_v3(const Int_t& b) { fHLT_PAAK4CaloJet100_Eta5p1_v3 = (Char_t)b; }
    void setHLT_PAAK4PFJet60_Eta5p1_v4(const Int_t& b) { fHLT_PAAK4PFJet60_Eta5p1_v4 = (Char_t)b; }
    void setHLT_PAAK4PFJet80_Eta5p1_v3(const Int_t& b) { fHLT_PAAK4PFJet80_Eta5p1_v3 = (Char_t)b; }
    void setHLT_PAAK4PFJet100_Eta5p1_v3(const Int_t& b) { fHLT_PAAK4PFJet100_Eta5p1_v3 = (Char_t)b; }
    void setHLT_PAAK4PFJet120_Eta5p1_v2(const Int_t& b) { fHLT_PAAK4PFJet120_Eta5p1_v2 = (Char_t)b; }

    void setHLT_HIAK4PFJet15_v1(const Int_t& b)        { fHLT_HIAK4PFJet15_v1 = (Char_t)b; }
    void setHLT_HIAK4PFJet15_v1_Prescl(const Int_t& b) { fHLT_HIAK4PFJet15_v1_Prescl = (Char_t)b; }  
    void setHLT_HIAK4PFJet30_v1(const Int_t& b)        { fHLT_HIAK4PFJet30_v1 = (Char_t)b; }
    void setHLT_HIAK4PFJet30_v1_Prescl(const Int_t& b) { fHLT_HIAK4PFJet30_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK4PFJet40_v1(const Int_t& b)        { fHLT_HIAK4PFJet40_v1 = (Char_t)b; }
    void setHLT_HIAK4PFJet40_v1_Prescl(const Int_t& b) { fHLT_HIAK4PFJet40_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK4PFJet60_v1(const Int_t& b)        { fHLT_HIAK4PFJet60_v1 = (Char_t)b; }
    void setHLT_HIAK4PFJet60_v1_Prescl(const Int_t& b) { fHLT_HIAK4PFJet60_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK4PFJet80_v1(const Int_t& b)        { fHLT_HIAK4PFJet80_v1 = (Char_t)b; }
    void setHLT_HIAK4PFJet80_v1_Prescl(const Int_t& b) { fHLT_HIAK4PFJet80_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK4PFJet120_v1(const Int_t& b)       { fHLT_HIAK4PFJet120_v1 = (Char_t)b; }
    void setHLT_HIAK4PFJet120_v1_Prescl(const Int_t& b){ fHLT_HIAK4PFJet120_v1_Prescl = (Char_t)b; }

    void setHLT_HIAK8PFJet15_v1(const Int_t& b)         { fHLT_HIAK8PFJet15_v1 = (Char_t)b; }
    void setHLT_HIAK8PFJet15_v1_Prescl(const Int_t& b)  { fHLT_HIAK8PFJet15_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK8PFJet25_v1(const Int_t& b)         { fHLT_HIAK8PFJet25_v1 = (Char_t)b; }
    void setHLT_HIAK8PFJet25_v1_Prescl(const Int_t& b)  { fHLT_HIAK8PFJet25_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK8PFJet40_v1(const Int_t& b)         { fHLT_HIAK8PFJet40_v1 = (Char_t)b; }
    void setHLT_HIAK8PFJet40_v1_Prescl(const Int_t& b)  { fHLT_HIAK8PFJet40_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK8PFJet60_v1(const Int_t& b)         { fHLT_HIAK8PFJet60_v1 = (Char_t)b; }
    void setHLT_HIAK8PFJet60_v1_Prescl(const Int_t& b)  { fHLT_HIAK8PFJet60_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK8PFJet80_v1(const Int_t& b)         { fHLT_HIAK8PFJet80_v1 = (Char_t)b; }
    void setHLT_HIAK8PFJet80_v1_Prescl(const Int_t& b)  { fHLT_HIAK8PFJet80_v1_Prescl = (Char_t)b; }
    void setHLT_HIAK8PFJet140_v1(const Int_t& b)        { fHLT_HIAK8PFJet140_v1 = (Char_t)b; }
    void setHLT_HIAK8PFJet140_v1_Prescl(const Int_t& b) { fHLT_HIAK8PFJet140_v1_Prescl = (Char_t)b; }

    void setHLT_HIPFJet25_v1(const Int_t& b)            { fHLT_HIPFJet25_v1 = (Char_t)b; }
    void setHLT_HIPFJet25_v1_Prescl(const Int_t& b)     { fHLT_HIPFJet25_v1_Prescl = (Char_t)b; }
    void setHLT_HIPFJet140_v1(const Int_t& b)           { fHLT_HIPFJet140_v1 = (Char_t)b; }
    void setHLT_HIPFJet140_v1_Prescl(const Int_t& b)    { fHLT_HIPFJet140_v1_Prescl = (Char_t)b; }

    void setHLT_HIPuAK4CaloJet80Eta5p1_v1(const Int_t& b)  { fHLT_HIPuAK4CaloJet80Eta5p1_v1 = (Char_t)b; }
    void setHLT_HIPuAK4CaloJet100Eta5p1_v1(const Int_t& b) { fHLT_HIPuAK4CaloJet100Eta5p1_v1 = (Char_t)b; }

    // Set skim flags
    void setHBHENoiseFilterResultRun2Loose(const Int_t& b) { fHBHENoiseFilterResultRun2Loose = (Char_t)b; }
    void setHBHENoiseFilterResultRun2Tight(const Int_t& b) { fHBHENoiseFilterResultRun2Tight = (Char_t)b; }
    void setHBHEIsoNoiseFilterResult(const Int_t& b)       { fHBHEIsoNoiseFilterResult = (Char_t)b; }
    void setCollisionEventSelectionAODv2(const Int_t& b)   { fCollisionEventSelectionAODv2 = (Char_t)b; }
    void setPhfCoincFilter2Th4(const Int_t& b)             { fPhfCoincFilter2Th4 = (Char_t)b; }
    void setPPAprimaryVertexFilter(const Int_t& b)         { fPPAprimaryVertexFilter = (Char_t)b; }
    void setPBeamScrapingFilter(const Int_t& b)            { fPBeamScrapingFilter = (Char_t)b; }
    void setPprimaryVertexFilter(const Int_t& b)           { fPprimaryVertexFilter = (Char_t)b; }
    void setPVertexFilterCutG(const Int_t& b)              { fPVertexFilterCutG = (Char_t)b; }
    void setPVertexFilterCutGloose(const Int_t& b)         { fPVertexFilterCutGloose = (Char_t)b; }
    void setPVertexFilterCutGtight(const Int_t& b)         { fPVertexFilterCutGtight = (Char_t)b; }
    void setPVertexFilterCutE(const Int_t& b)              { fPVertexFilterCutE = (Char_t)b; }
    void setPVertexFilterCutEandG(const Int_t& b)          { fPVertexFilterCutEandG = (Char_t)b; }
    void setPClusterCompatibilityFilter(const Int_t& b)    { fPClusterCompatibilityFilter = (Char_t)b; }

    void setPhfCoincFilter(const Int_t& b) { fPhfCoincFilter = (Char_t)b; }
    void setPVertexFilterCutdz1p0(const Int_t& b) { fPVertexFilterCutdz1p0 = (Char_t)b; }
    void setPVertexFilterCutGplus(const Int_t& b) { fPVertexFilterCutGplus = (Char_t)b; }
    void setPVertexFilterCutVtx1(const Int_t& b) { fPVertexFilterCutVtx1 = (Char_t)b; }

    //
    // Getters
    //

    // HLT flags
    Int_t HLT_PAAK4CaloJet60_Eta5p1_v3() const { return (Int_t)fHLT_PAAK4CaloJet60_Eta5p1_v3; }
    Int_t HLT_PAAK4CaloJet80_Eta5p1_v3() const { return (Int_t)fHLT_PAAK4CaloJet80_Eta5p1_v3; }
    Int_t HLT_PAAK4CaloJet100_Eta5p1_v3() const { return (Int_t)fHLT_PAAK4CaloJet100_Eta5p1_v3; }
    Int_t HLT_PAAK4PFJet60_Eta5p1_v4() const { return (Int_t)fHLT_PAAK4PFJet60_Eta5p1_v4; }
    Int_t HLT_PAAK4PFJet80_Eta5p1_v3() const { return (Int_t)fHLT_PAAK4PFJet80_Eta5p1_v3; }
    Int_t HLT_PAAK4PFJet100_Eta5p1_v3() const { return (Int_t)fHLT_PAAK4PFJet100_Eta5p1_v3; }
    Int_t HLT_PAAK4PFJet120_Eta5p1_v2() const { return (Int_t)fHLT_PAAK4PFJet120_Eta5p1_v2; }

    Int_t HLT_HIAK4PFJet15_v1() const { return (Int_t)fHLT_HIAK4PFJet15_v1; }
    Int_t HLT_HIAK4PFJet15_v1_Prescl() const { return (Int_t)fHLT_HIAK4PFJet15_v1_Prescl; }
    Int_t HLT_HIAK4PFJet30_v1() const { return (Int_t)fHLT_HIAK4PFJet30_v1; }
    Int_t HLT_HIAK4PFJet30_v1_Prescl() const { return (Int_t)fHLT_HIAK4PFJet30_v1_Prescl; }
    Int_t HLT_HIAK4PFJet40_v1() const { return (Int_t)fHLT_HIAK4PFJet40_v1; }
    Int_t HLT_HIAK4PFJet40_v1_Prescl() const { return (Int_t)fHLT_HIAK4PFJet40_v1_Prescl; }
    Int_t HLT_HIAK4PFJet60_v1() const { return (Int_t)fHLT_HIAK4PFJet60_v1; }
    Int_t HLT_HIAK4PFJet60_v1_Prescl() const { return (Int_t)fHLT_HIAK4PFJet60_v1_Prescl; }
    Int_t HLT_HIAK4PFJet80_v1() const { return (Int_t)fHLT_HIAK4PFJet80_v1; }
    Int_t HLT_HIAK4PFJet80_v1_Prescl() const { return (Int_t)fHLT_HIAK4PFJet80_v1_Prescl; }
    Int_t HLT_HIAK4PFJet120_v1() const { return (Int_t)fHLT_HIAK4PFJet120_v1; }
    Int_t HLT_HIAK4PFJet120_v1_Prescl() const { return (Int_t)fHLT_HIAK4PFJet120_v1_Prescl; }

    Int_t HLT_HIAK8PFJet15_v1() const { return (Int_t)fHLT_HIAK8PFJet15_v1; }
    Int_t HLT_HIAK8PFJet15_v1_Prescl() const { return (Int_t)fHLT_HIAK8PFJet15_v1_Prescl; }
    Int_t HLT_HIAK8PFJet25_v1() const { return (Int_t)fHLT_HIAK8PFJet25_v1; }
    Int_t HLT_HIAK8PFJet25_v1_Prescl() const { return (Int_t)fHLT_HIAK8PFJet25_v1_Prescl; }
    Int_t HLT_HIAK8PFJet40_v1() const { return (Int_t)fHLT_HIAK8PFJet40_v1; }
    Int_t HLT_HIAK8PFJet40_v1_Prescl() const { return (Int_t)fHLT_HIAK8PFJet40_v1_Prescl; }
    Int_t HLT_HIAK8PFJet60_v1() const { return (Int_t)fHLT_HIAK8PFJet60_v1; }
    Int_t HLT_HIAK8PFJet60_v1_Prescl() const { return (Int_t)fHLT_HIAK8PFJet60_v1_Prescl; }
    Int_t HLT_HIAK8PFJet80_v1() const { return (Int_t)fHLT_HIAK8PFJet80_v1; }
    Int_t HLT_HIAK8PFJet80_v1_Prescl() const { return (Int_t)fHLT_HIAK8PFJet80_v1_Prescl; }
    Int_t HLT_HIAK8PFJet140_v1() const { return (Int_t)fHLT_HIAK8PFJet140_v1; }
    Int_t HLT_HIAK8PFJet140_v1_Prescl() const { return (Int_t)fHLT_HIAK8PFJet140_v1_Prescl; }

    Int_t HLT_HIPFJet25_v1() const { return (Int_t)fHLT_HIPFJet25_v1; }
    Int_t HLT_HIPFJet25_v1_Prescl() const { return (Int_t)fHLT_HIPFJet25_v1_Prescl; }
    Int_t HLT_HIPFJet140_v1() const { return (Int_t)fHLT_HIPFJet25_v1; }
    Int_t HLT_HIPFJet140_v1_Prescl() const { return (Int_t)fHLT_HIPFJet25_v1_Prescl; }

    Int_t HLT_HIPuAK4CaloJet80Eta5p1_v1() const { return (Int_t)fHLT_HIPuAK4CaloJet80Eta5p1_v1; }
    Int_t HLT_HIPuAK4CaloJet100Eta5p1_v1() const { return (Int_t)fHLT_HIPuAK4CaloJet100Eta5p1_v1; }

    // Skim flags
    Int_t HBHENoiseFilterResultRun2Loose() const { return (Int_t)fHBHENoiseFilterResultRun2Loose; }
    Int_t HBHENoiseFilterResultRun2Tight() const { return (Int_t)fHBHENoiseFilterResultRun2Tight; }
    Int_t HBHEIsoNoiseFilterResult() const { return (Int_t)fHBHEIsoNoiseFilterResult; }
    Int_t collisionEventSelectionAODv2() const { return (Int_t)fCollisionEventSelectionAODv2; }
    Int_t phfCoincFilter2Th4() const { return (Int_t)fPhfCoincFilter2Th4; }
    Int_t pPAprimaryVertexFilter() const { return (Int_t)fPPAprimaryVertexFilter; }
    Int_t pBeamScrapingFilter() const { return (Int_t)fPBeamScrapingFilter; }
    Int_t pprimaryVertexFilter() const { return (Int_t)fPprimaryVertexFilter; }
    Int_t pVertexFilterCutG() const { return (Int_t)fPVertexFilterCutG; }

    Int_t pVertexFilterCutGloose() const { return (Int_t)fPVertexFilterCutGloose; }
    Int_t pVertexFilterCutGtight() const { return (Int_t)fPVertexFilterCutGtight; }
    Int_t pVertexFilterCutE() const { return (Int_t)fPVertexFilterCutE; }
    Int_t pVertexFilterCutEandG() const { return (Int_t)fPVertexFilterCutEandG; }
    Int_t pClusterCompatibilityFilter() const {return (Int_t)fPClusterCompatibilityFilter; }

    Int_t phfCoincFilter() const        { return (Int_t)fPhfCoincFilter; }
    Int_t pVertexFilterCutdz1p0() const { return (Int_t)fPVertexFilterCutdz1p0; }
    Int_t pVertexFilterCutGplus() const { return (Int_t)fPVertexFilterCutGplus; }
    Int_t pVertexFilterCutVtx1() const  { return (Int_t)fPVertexFilterCutVtx1; }

  private:

    // HLT part
    Char_t fHLT_PAAK4CaloJet60_Eta5p1_v3;
    Char_t fHLT_PAAK4CaloJet80_Eta5p1_v3;
    Char_t fHLT_PAAK4CaloJet100_Eta5p1_v3;
    Char_t fHLT_PAAK4PFJet60_Eta5p1_v4;
    Char_t fHLT_PAAK4PFJet80_Eta5p1_v3;
    Char_t fHLT_PAAK4PFJet100_Eta5p1_v3;
    Char_t fHLT_PAAK4PFJet120_Eta5p1_v2;

    Char_t fHLT_HIAK4PFJet15_v1;
    Char_t fHLT_HIAK4PFJet15_v1_Prescl;
    Char_t fHLT_HIAK4PFJet30_v1;
    Char_t fHLT_HIAK4PFJet30_v1_Prescl;
    Char_t fHLT_HIAK4PFJet40_v1;
    Char_t fHLT_HIAK4PFJet40_v1_Prescl;
    Char_t fHLT_HIAK4PFJet60_v1;
    Char_t fHLT_HIAK4PFJet60_v1_Prescl;
    Char_t fHLT_HIAK4PFJet80_v1;
    Char_t fHLT_HIAK4PFJet80_v1_Prescl;
    Char_t fHLT_HIAK4PFJet120_v1;
    Char_t fHLT_HIAK4PFJet120_v1_Prescl;

    Char_t fHLT_HIAK8PFJet15_v1;
    Char_t fHLT_HIAK8PFJet15_v1_Prescl;
    Char_t fHLT_HIAK8PFJet25_v1;
    Char_t fHLT_HIAK8PFJet25_v1_Prescl;
    Char_t fHLT_HIAK8PFJet40_v1;
    Char_t fHLT_HIAK8PFJet40_v1_Prescl;
    Char_t fHLT_HIAK8PFJet60_v1;
    Char_t fHLT_HIAK8PFJet60_v1_Prescl;
    Char_t fHLT_HIAK8PFJet80_v1;
    Char_t fHLT_HIAK8PFJet80_v1_Prescl;
    Char_t fHLT_HIAK8PFJet140_v1;
    Char_t fHLT_HIAK8PFJet140_v1_Prescl;

    Char_t fHLT_HIPFJet25_v1;
    Char_t fHLT_HIPFJet25_v1_Prescl;
    Char_t fHLT_HIPFJet140_v1;
    Char_t fHLT_HIPFJet140_v1_Prescl;

    Char_t fHLT_HIPuAK4CaloJet80Eta5p1_v1;
    Char_t fHLT_HIPuAK4CaloJet100Eta5p1_v1;


    // Skimanalysis part
    Char_t fHBHENoiseFilterResultRun2Loose;
    Char_t fHBHENoiseFilterResultRun2Tight;
    Char_t fHBHEIsoNoiseFilterResult;
    Char_t fCollisionEventSelectionAODv2;
    Char_t fPhfCoincFilter2Th4;
    Char_t fPPAprimaryVertexFilter;
    Char_t fPBeamScrapingFilter;
    Char_t fPprimaryVertexFilter;
    Char_t fPVertexFilterCutG;
    Char_t fPVertexFilterCutGloose;
    Char_t fPVertexFilterCutGtight;
    Char_t fPVertexFilterCutE;
    Char_t fPVertexFilterCutEandG;
    Char_t fPClusterCompatibilityFilter;

    Char_t fPhfCoincFilter;
    Char_t fPVertexFilterCutdz1p0;
    Char_t fPVertexFilterCutGplus;
    Char_t fPVertexFilterCutVtx1;

    ClassDef(TriggerAndSkim, 0)
};

#endif // #define TriggerAndSkim_h