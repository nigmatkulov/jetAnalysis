#ifndef Collections_h
#define Collections_h

// Jeta analysis headers
#include "BaseAnalysis.h"
#include "Jet.h"
#include "Track.h"
#include "GenTrack.h"

// C++ headers
#include <list>

// Analysis collection
typedef std::list<BaseAnalysis*> AnalysisCollection;
typedef std::list<BaseAnalysis*>::iterator AnalysisIterator;

// Particle flow jet collection
typedef std::list<Jet*> PartFlowJetCollection;
typedef std::list<Jet*>::iterator PartFlowJetIterator;

// Calorimeter jet collection
typedef std::list<Jet*> CaloJetCollection;
typedef std::list<Jet*>::iterator CaloJetIterator;

// Reconstructed track collection
typedef std::list<Track*> TrackCollection;
typedef std::list<Track*>::iterator TrackIterator;

// Generated track collection
typedef std::list<GenTrack*> GenTrackCollection;
typedef std::list<GenTrack*>::iterator GenTrackIterator;

#endif // #define Collections_h