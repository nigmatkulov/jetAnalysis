/**
 * @file EventCut.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Example of the event cut
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef EventCut_h
#define EventCut_h

// Jet analysis headers
#include "Event.h"

// ROOT headers
#include "Rtypes.h"

// C++ headers
#include <limits>

//________________
class EventCut {
  public:
    /// @brief Base constructor
    EventCut();
    /// @brief Destructor
    virtual ~EventCut();

    //
    // Setters
    //

    /// @brief Set limits on Vx
    void setVx(const double& lo = -1e9, const double& hi = 1e9) { fVx[0]=lo; fVx[1]=hi; }
    /// @brief Set limits on Vy
    void setVy(const double& lo = -1e9, const double& hi = 1e9) { fVy[0]=lo; fVy[1]=hi; }
    /// @brief Set limits on Vz
    void setVz(const double& lo = -1e9, const double& hi = 1e9) { fVz[0]=lo; fVz[1]=hi; }
    /// @brief Shift of the nominal beam position in x direction
    void setShiftVx(const double& x = 0.) { fShiftVx = x; }
    /// @brief Shift of the nominal beam position in y direction
    void setShiftVy(const double& y = 0.) { fShiftVy = y; }
    /// @brief Set cut on maximal radial vertex displacement
    void setVertexR(const double& r = 1e9) { fVR = r; }
    /// @brief Set limits on HiBin
    void setHiBin(const short& lo = -10000, const short& hi = 10000) { fHiBin[0]=lo; fHiBin[1]=hi; }
    /// @brief Set limits on centrality
    void setCentrality(const double& lo=-1000, const double& hi = 1000) { fCentVal[0]=lo; fCentVal[1]=hi;}
    /// @brief Set limits on luminosity
    void setLumi(const unsigned int& lo = 0, const unsigned int& hi = std::numeric_limits<unsigned int>::max()) { fLumi[0]=lo; fLumi[1]=hi; }
    /// @brief Set limits on pT hat
    void setPtHat(const double& lo=-1e9, const double& hi=1e9) { fPtHat[0]=lo; fPtHat[1]=hi; }
    /// @brief Set limits on event weight
    void setPtHatWeight(const double& lo=-1e9, const double& hi=1e9) { fPtHatWeight[0]=lo; fPtHatWeight[1]=hi; }
    /// @brief Print information each event
    void setVerbose() { fVerbose = {true}; }
    // Skim selection criteria
    void usePPrimaryVertexFilter()           { fPPrimaryVertexFilter = {true}; }
    void useHBHENoiseFilterResultRun2Loose() { fHBHENoiseFilterResultRun2Loose = {true}; }
    void useCollisionEventSelectionAODv2()   { fCollisionEventSelectionAODc2 = {true}; }
    void usePhfCoincFilter2Th4()             { fPhfCoincFilter2Th4 = {true}; }
    void usePPAprimaryVertexFilter()         { fPPAprimaryVertexFilter = {true}; }
    void usePBeamScrapingFilter()            { fPBeamScrapingFilter = {true}; }
    void usePClusterCompatibilityFilter()    { fPClusterCompatibilityFilter = {true}; }
    void useHLT_HIPuAK4CaloJet80Eta5p1_v1()  { fHLT_HIPuAK4CaloJet80Eta5p1_v1 = {true}; }
    void useHLT_PAAK4PFJet60_Eta5p1_v4()     { fHLT_PAAK4PFJet60_Eta5p1_v4 = {true}; }
    void useHLT_PAAK4PFJet80_Eta5p1_v3()     { fHLT_PAAK4PFJet80_Eta5p1_v3 = {true}; }
    void useHLT_PAAK4PFJet100_Eta5p1_v3()    { fHLT_PAAK4PFJet100_Eta5p1_v3 = {true}; }
    void useHLT_PAAK4PFJet120_Eta5p1_v2()    { fHLT_PAAK4PFJet120_Eta5p1_v2 = {true}; }
    void useHLT_HIAK4CaloJet60_v1()          { fHLT_HIAK4CaloJet60_v1 = {true}; }
    void useHLT_HIAK4CaloJet80_v1()          { fHLT_HIAK4CaloJet80_v1 = {true}; }
    void useHLT_HIAK4PFJet60_v1()            { fHLT_HIAK4PFJet60_v1 = {true}; }
    void useHLT_HIAK4PFJet80_v1()            { fHLT_HIAK4PFJet80_v1 = {true}; }

    void usePhfCoincFilter()                { fPhfCoincFilter = {true}; }
    void usePVertexFilterCutdz1p0()         { fPVertexFilterCutdz1p0 = {true}; }
    void usePVertexFilterCutGplus()         { fPVertexFilterCutGplus = {true}; }
    void usePVertexFilterCutVtx1()          { fPVertexFilterCutVtx1 = {true}; }
    /// @brief Report information about
    void report();
    /// @brief Check if evn 
    virtual bool pass(const Event* ev);

  private:
    /// @brief Vertex x interval
    double fVx[2];
    /// @brief Vertex y interval
    double fVy[2];
    /// @brief Vertex z interval
    double fVz[2];

    /// @brief Beam nominal shift in X direction
    double fShiftVx;
    /// @brief Beam nominal shift in Y direction
    double fShiftVy;
    /// @brief Cut on radial displacement
    double fVR;
    /// @brief Cut on HiBin
    short  fHiBin[2];
    /// @brief Cut on centrality values (assuming 0.5% centrality step from 0-100%)
    double  fCentVal[2];
    /// @brief Cut on luminosity value
    unsigned int fLumi[2];
    /// @brief PtHat cut
    double fPtHat[2];
    /// @brief Event weight cut
    double fPtHatWeight[2];
    /// @brief  Print information each time
    bool fVerbose;

    // Skim flags 
    bool fPPrimaryVertexFilter;
    bool fHBHENoiseFilterResultRun2Loose;
    bool fCollisionEventSelectionAODc2;
    bool fPhfCoincFilter2Th4;

    bool fPPAprimaryVertexFilter;
    bool fPBeamScrapingFilter;
    bool fPClusterCompatibilityFilter;
    bool fPhfCoincFilter;
    bool fPVertexFilterCutdz1p0;
    bool fPVertexFilterCutGplus;
    bool fPVertexFilterCutVtx1;

    // Triggers
    bool fHLT_HIAK4CaloJet60_v1; // pp5020
    bool fHLT_HIAK4CaloJet80_v1; // pp5020

    bool fHLT_HIPuAK4CaloJet80Eta5p1_v1;
    bool fHLT_PAAK4PFJet60_Eta5p1_v4;  // pPb 8160
    bool fHLT_PAAK4PFJet80_Eta5p1_v3;  // pPb 8160
    bool fHLT_PAAK4PFJet100_Eta5p1_v3; // pPb 8160
    bool fHLT_PAAK4PFJet120_Eta5p1_v2; // pPb 8160

    bool fHLT_HIAK4PFJet60_v1;         // pp 5020
    bool fHLT_HIAK4PFJet80_v1;         // pp 5020

    /// @brief Number of events passed cut
    Long64_t fEventsPassed;
    /// @brief Number of events failed cut
    Long64_t fEventsFailed;

    ClassDef(EventCut, 0)
};

#endif // #define EventCut_h