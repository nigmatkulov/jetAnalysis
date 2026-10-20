#pragma once
#include <array>
#include <cstdint>
#include <vector>
namespace jetanalysis {
struct GenJet {
    float pt;
    float eta;
    float phi;
    float smearedPt;
};

//________________
struct RecoJet {
    float recoPtRaw;
    float recoPt;
    float recoEta;
    float recoPhi;
    float recoPfNHF;
    float recoPfNEF;
    float recoPfCHF;
    float recoPfMUF;
    float recoPfCEF;
    std::int32_t recoPfCHM;
    std::int32_t recoPfCEM;
    std::int32_t recoPfNHM;
    std::int32_t recoPfNEM;
    std::int32_t recoPfMUM;
    float recoTrackMaxPt;
    float refPt;
    float refEta;
    float refPhi;

    float recoPtJeuUp;
    float recoPtJeuDown;

    // Historical name: in Hybrid/Stochastic modes this is SMEARED RECO pT;
    // only TruthToPseudoReco constructs it around the matched reference pT.
    float refPtSmeared; // unit-SF JER with the base resolution
    float refPtSmearedEta = -1.f; // unit-SF JER with eta-dependent resolution
    float refPtSmearedJerDef = -1.f;
    float refPtSmearedJerUp = -1.f;
    float refPtSmearedJerDown = -1.f;
    // All existing output branches: unit, x2, x0.9, base Down/Def/Up,
    // eta-dependent Down/Def/Up. Cache travels with the jet through sorting.
    std::array<float, 9> jerMomenta{};
    std::int32_t smearingIndex = -1; // index in initial corrected-pT ordering, before response sorts
};

struct TriggerBits {
    std::int32_t HLT_PAAK4PFJet60_Eta5p1_v4=0, HLT_PAAK4PFJet80_Eta5p1_v3=0, HLT_PAAK4PFJet100_Eta5p1_v3=0;
};
struct FilterBits {
    std::int32_t pBeamScrapingFilter=0, pPAprimaryVertexFilter=0, HBHENoiseFilterResultRun2Loose=0,
        phfCoincFilter=0, pVertexFilterCutdz1p0=0, pVertexFilterCutGplus=0, pVertexFilterCutVtx1=0;
};
struct EventInfo : TriggerBits, FilterBits {
    std::uint32_t run=0, lumi=0;
    std::uint64_t evt=0;
    float ptHat=0, inputEventWeight=0, vz=0, eventRho=0;
};
struct Track { float pt=0, eta=0, phi=0; std::int32_t charge=0; bool highPurity=false; };
struct GenParticle { float pt=0, eta=0, phi=0; std::int32_t charge=0, pdg=0, subevent=0; };
enum class SampleType { Data, Embedding, Pythia };
struct Event {
    EventInfo info;
    SampleType sample=SampleType::Data;
    bool hasRecoJets=false, hasGenJets=false, hasTracks=false, hasParticles=false;
    std::vector<RecoJet> recoJets; // raw input, before analysis corrections
    std::vector<GenJet> genJets;
    std::vector<Track> tracks;
    std::vector<GenParticle> particles;
};
} // namespace jetanalysis
