/**
 * @file Collections.h
 * @author Grigory Nigmatkulov (you@domain.com)
 * @brief Header that contains definition of collection wrappers
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef Collections_h
#define Collections_h

// Jeta analysis headers
#include "BaseAnalysis.h"
#include "RecoJet.h"
#include "GenJet.h"
#include "Track.h"
#include "GenTrack.h"

// C++ headers
#include <vector>

// Analysis collection
typedef std::vector<BaseAnalysis*> AnalysisCollection;
typedef std::vector<BaseAnalysis*>::iterator AnalysisIterator;

// Reconstructed jet collection
typedef std::vector<RecoJet*> RecoJetCollection;
typedef std::vector<RecoJet*>::iterator RecoJetIterator;

// Generated jet collection
typedef std::vector<GenJet*> GenJetCollection;
typedef std::vector<GenJet*>::iterator GenJetIterator;

// Reconstructed track collection
typedef std::vector<Track*> TrackCollection;
typedef std::vector<Track*>::iterator TrackIterator;

// Generated track collection
typedef std::vector<GenTrack*> GenTrackCollection;
typedef std::vector<GenTrack*>::iterator GenTrackIterator;

#endif // #define Collections_h