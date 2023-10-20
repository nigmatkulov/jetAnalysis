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
#include "Jet.h"
#include "Track.h"
#include "GenTrack.h"

// C++ headers
#include <vector>

// Analysis collection
typedef std::vector<BaseAnalysis*> AnalysisCollection;
typedef std::vector<BaseAnalysis*>::iterator AnalysisIterator;

// Particle flow jet collection
typedef std::vector<Jet*> PartFlowJetCollection;
typedef std::vector<Jet*>::iterator PartFlowJetIterator;

// Calorimeter jet collection
typedef std::vector<Jet*> CaloJetCollection;
typedef std::vector<Jet*>::iterator CaloJetIterator;

// Reconstructed track collection
typedef std::vector<Track*> TrackCollection;
typedef std::vector<Track*>::iterator TrackIterator;

// Generated track collection
typedef std::vector<GenTrack*> GenTrackCollection;
typedef std::vector<GenTrack*>::iterator GenTrackIterator;

#endif // #define Collections_h