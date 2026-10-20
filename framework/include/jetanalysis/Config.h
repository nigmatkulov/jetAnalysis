#pragma once
#include <string>
#include <vector>
namespace jetanalysis {
struct ReaderConfig {
    std::string input, jetTree="ak4PFJetAnalyzer/t", trackTree="ppTrack/trackTree", particleTree="HiGenParticleAna/hi";
    int mcType=2;
    bool pp=false, readRecoJets=true, readGenJets=true, readTracks=false, readParticles=false;
};
struct AnalysisConfig {
    ReaderConfig reader;
    bool profileStages=false;
    std::string output, auxDirectory="./aux_files/pPb_8160";
    int isPbGoingDir=1, ptHatSample=50, triggerId=0, recoJetSelMethod=2,
        vertexFilterSelection=0, jerSmearingMode=0, useEtaDependentJER=0;
    std::vector<std::string> jecFiles, excludeOutputPrefixes;
    std::string jeuFile;
};
} // namespace jetanalysis
