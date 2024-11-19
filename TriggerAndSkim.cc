/**
 * @file TriggerAndSkim.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Trigger and skim properties of the event
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "TriggerAndSkim.h"

// C++ headers
#include <limits>

ClassImp(TriggerAndSkim)

//________________
TriggerAndSkim::TriggerAndSkim() :
  fHLT_HIAK4CaloJet60_v1{0},
  fHLT_HIAK4CaloJet80_v1{0},
  fHLT_PAAK4CaloJet60_Eta5p1_v3{0},
  fHLT_PAAK4CaloJet80_Eta5p1_v3{0},
  fHLT_PAAK4CaloJet100_Eta5p1_v3{0},
  fHLT_PAAK4PFJet60_Eta5p1_v4{0},
  fHLT_PAAK4PFJet80_Eta5p1_v3{0},
  fHLT_PAAK4PFJet100_Eta5p1_v3{0},
  fHLT_PAAK4PFJet120_Eta5p1_v2{0},
  fHLT_HIAK4PFJet15_v1{0},
  fHLT_HIAK4PFJet15_v1_Prescl{0},
  fHLT_HIAK4PFJet30_v1{0},
  fHLT_HIAK4PFJet30_v1_Prescl{0},
  fHLT_HIAK4PFJet40_v1{0},
  fHLT_HIAK4PFJet40_v1_Prescl{0},
  fHLT_HIAK4PFJet60_v1{0},
  fHLT_HIAK4PFJet60_v1_Prescl{0},
  fHLT_HIAK4PFJet80_v1{0},
  fHLT_HIAK4PFJet80_v1_Prescl{0},
  fHLT_HIAK4PFJet120_v1{0},
  fHLT_HIAK4PFJet120_v1_Prescl{0},

  fHLT_HIAK8PFJet15_v1{0},
  fHLT_HIAK8PFJet15_v1_Prescl{0},
  fHLT_HIAK8PFJet25_v1{0},
  fHLT_HIAK8PFJet25_v1_Prescl{0},
  fHLT_HIAK8PFJet40_v1{0},
  fHLT_HIAK8PFJet40_v1_Prescl{0},
  fHLT_HIAK8PFJet60_v1{0},
  fHLT_HIAK8PFJet60_v1_Prescl{0},
  fHLT_HIAK8PFJet80_v1{0},
  fHLT_HIAK8PFJet80_v1_Prescl{0},
  fHLT_HIAK8PFJet140_v1{0},
  fHLT_HIAK8PFJet140_v1_Prescl{0},

  fHLT_HIPFJet25_v1{0},
  fHLT_HIPFJet25_v1_Prescl{0},
  fHLT_HIPFJet140_v1{0},
  fHLT_HIPFJet140_v1_Prescl{0},

  fHLT_HIPuAK4CaloJet80Eta5p1_v1{0},
  fHLT_HIPuAK4CaloJet100Eta5p1_v1{0},
  
  fHBHENoiseFilterResultRun2Loose{0},
  fHBHENoiseFilterResultRun2Tight{0},
  fHBHEIsoNoiseFilterResult{0},
  fCollisionEventSelectionAODv2{0},
  fPhfCoincFilter2Th4{0},
  fPPAprimaryVertexFilter{0},
  fPBeamScrapingFilter{0},
  fPprimaryVertexFilter{0},
  fPVertexFilterCutG{0},
  fPVertexFilterCutGloose{0},
  fPVertexFilterCutGtight{0},
  fPVertexFilterCutE{0},
  fPVertexFilterCutEandG{0},
  fPClusterCompatibilityFilter{0},
  
  fPhfCoincFilter{0},
  fPVertexFilterCutdz1p0{0},
  fPVertexFilterCutGplus{0},
  fPVertexFilterCutVtx1{0} {
    /* empty */
}