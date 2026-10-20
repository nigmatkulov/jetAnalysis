#include "jetanalysis/RootBranchContract.h"
#include <TFile.h>
#include <iostream>
#include <memory>

int main(int argc, char** argv) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: check_forest_schema file.root mcType pPb|pp [jetTree=ak4PFJetAnalyzer/t]\n";
        return 2;
    }
    try {
        const std::string mcType = argv[2], system = argv[3];
        if ((mcType != "0" && mcType != "1" && mcType != "2") ||
            (system != "pPb" && system != "pp"))
            throw std::runtime_error("Invalid mcType or collision system");
        const bool mc = mcType != "0", pp = system == "pp";
        std::unique_ptr<TFile> file(TFile::Open(argv[1], "READ"));
        if (!file || file->IsZombie()) throw std::runtime_error("Cannot open input file");
        auto tree = [&](const char* path) -> TTree& {
            auto* result = file->Get<TTree>(path);
            if (!result) throw std::runtime_error(std::string("Missing TTree: ") + path);
            return *result;
        };
        using namespace jetanalysis;
        auto& events = tree("hiEvtAnalyzer/HiTree");
        requireLeaf<Float_t>(events, "vz");
        requireLeaf<UInt_t>(events, "run");
        requireLeaf<UInt_t>(events, "lumi");
        requireLeaf<ULong64_t>(events, "evt");
        if (mc) requireLeaf<Float_t>(events, "pthat");
        if (mc && pp) requireLeaf<Float_t>(events, "weight");
        if (events.GetBranch("rho")) requireLeaf<Float_t>(events, "rho");
        else if (events.GetBranch("fixedGridRhoFastjetAll"))
            requireLeaf<Float_t>(events, "fixedGridRhoFastjetAll");
        if (!pp) {
            auto& skim = tree("skimanalysis/HltTree");
            for (const char* name : {"pBeamScrapingFilter", "pPAprimaryVertexFilter",
                 "HBHENoiseFilterResultRun2Loose", "phfCoincFilter",
                 "pVertexFilterCutdz1p0", "pVertexFilterCutGplus", "pVertexFilterCutVtx1"})
                requireLeaf<Int_t>(skim, name);
        }
        if (!mc) {
            auto& hlt = tree("hltanalysis/HltTree");
            for (const char* name : {"HLT_PAAK4PFJet60_Eta5p1_v4",
                 "HLT_PAAK4PFJet80_Eta5p1_v3", "HLT_PAAK4PFJet100_Eta5p1_v3"})
                requireLeaf<Int_t>(hlt, name);
        }
        auto& jets = tree(argc == 5 ? argv[4] : "ak4PFJetAnalyzer/t");
        requireLeaf<Int_t>(jets, "nref");
        for (const char* name : {"rawpt", "trackMax", "jteta", "jtphi", "jtPfNHF",
             "jtPfNEF", "jtPfCHF", "jtPfMUF", "jtPfCEF"})
            requireLeaf<Float_t>(jets, name, "nref");
        for (const char* name : {"jtPfCHM", "jtPfCEM", "jtPfNHM", "jtPfNEM", "jtPfMUM"})
            requireLeaf<Int_t>(jets, name, "nref");
        if (mc) {
            requireLeaf<Int_t>(jets, "ngen");
            for (const char* name : {"genpt", "geneta", "genphi"})
                requireLeaf<Float_t>(jets, name, "ngen");
            for (const char* name : {"refpt", "refeta", "refphi"})
                requireLeaf<Float_t>(jets, name, "nref");
        }
        std::cout << "All branches used by the reference analysis have exact expected types and shapes.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
