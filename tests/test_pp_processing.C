// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C before this macro.
#include <stdexcept>

// Minimal forests exercise actual friend-chain branch activation and file changes.
void writePpBranchFixture(const TString &path, float scale) {
    TFile file(path, "RECREATE");
    Int_t zero = 0;
    Float_t value = 0.f, jetPtHat = 9999.f;
    file.mkdir("hltanalysis")->cd();
    TTree hlt("HltTree", "HLT"); hlt.Branch("dummy", &zero, "dummy/I");
    file.mkdir("skimanalysis")->cd();
    TTree skim("HltTree", "skim"); skim.Branch("dummy", &zero, "dummy/I");
    for (auto name : {"pBeamScrapingFilter", "pPAprimaryVertexFilter",
                      "HBHENoiseFilterResultRun2Loose", "phfCoincFilter",
                      "pVertexFilterCutdz1p0", "pVertexFilterCutGplus", "pVertexFilterCutVtx1"})
        skim.Branch(name, &zero, Form("%s/I", name));
    file.mkdir("hiEvtAnalyzer")->cd();
    TTree events("HiTree", "events");
    Float_t eventPtHat, eventVz, eventWeightValue, rhoValue;
    events.Branch("rho", &rhoValue, "rho/F");
    UInt_t eventRun = 1, eventLumi = 2;
    ULong64_t eventId;
    events.Branch("run", &eventRun, "run/i");
    events.Branch("lumi", &eventLumi, "lumi/i");
    events.Branch("evt", &eventId, "evt/l");
    events.Branch("pthat", &eventPtHat, "pthat/F");
    events.Branch("vz", &eventVz, "vz/F");
    events.Branch("weight", &eventWeightValue, "weight/F");
    file.mkdir("ak4PFJetAnalyzer")->cd();
    TTree jets("t", "jets");
    jets.Branch("nref", &zero, "nref/I"); jets.Branch("ngen", &zero, "ngen/I");
    jets.Branch("pthat", &jetPtHat, "pthat/F");
    for (auto name : {"genpt", "geneta", "genphi"})
        jets.Branch(name, &value, Form("%s[ngen]/F", name));
    for (auto name : {"rawpt", "trackMax", "jteta", "jtphi", "jtPfNHF", "jtPfNEF",
                      "jtPfCHF", "jtPfMUF", "jtPfCEF", "refpt", "refeta", "refphi"})
        jets.Branch(name, &value, Form("%s[nref]/F", name));
    for (auto name : {"jtPfCHM", "jtPfCEM", "jtPfNHM", "jtPfNEM", "jtPfMUM"})
        jets.Branch(name, &zero, Form("%s[nref]/I", name));
    for (int i = 0; i < 2; ++i) {
        eventId = static_cast<ULong64_t>(scale) * 100 + i;
        rhoValue = scale * 10 + i;
        eventPtHat = i == 0 ? 1.f : 1000.f;
        eventVz = i == 0 ? 0.f : 14.f;
        eventWeightValue = (i == 0 ? 2.f : -0.5f) * scale;
        hlt.Fill(); skim.Fill(); events.Fill(); jets.Fill();
    }
    file.Write();
}

void test_pp_processing() {
    auto check = [](bool condition, const char *message) {
        if (!condition) throw std::runtime_error(message);
    };
    TF1 zeroVzWeight("testZeroVzWeight", "0", -20., 20.);
    for (float stored : {-2.5f, 0.f, 0.125f, 3.f}) {
        // Neither the zero vz function nor event-count normalization may affect pp.
        check(eventWeight(1.f, 0.f, zeroVzWeight, 0, true, stored) == stored,
              "pp must return the exact stored weight");
        check(eventWeight(1000.f, 14.f, zeroVzWeight, 12345, true, stored) == stored,
              "pp weight must be independent of pthat, vz and sample size");
    }
    bool rejected = false;
    try {
        eventWeight(30.f, 0.f, zeroVzWeight, 1, true,
                    std::numeric_limits<float>::quiet_NaN());
    } catch (const std::runtime_error &) { rejected = true; }
    check(rejected, "non-finite pp weights must fail");
    TF1 constantVzWeight("testConstantVzWeight", "2", -20., 20.);
    const double legacyWeight = 7.7966624e-08 * 952110 / 100. / 2.;
    check(eventWeight(40.f, 0.f, constantVzWeight, 100) == legacyWeight,
          "pPb normalization and vz weighting changed");

    ptHatRange[0] = 30.f; ptHatRange[1] = 50.f;
    pBeamScrapingFilter = pPAprimaryVertexFilter = HBHENoiseFilterResultRun2Loose = 0;
    phfCoincFilter = pVertexFilterCutdz1p0 = 0;
    for (float value : {1.f, 40.f, 1000.f}) {
        ptHat = value;
        for (float z : {-15.f, 0.f, 15.f}) {
            vz = z;
            check(isGoodEvent(true, true, 0, 0, true),
                  "pp must ignore pthat range and pPb filter flags");
        }
    }
    vz = 15.01f;
    check(!isGoodEvent(true, true, 0, 0, true), "pp must retain the vz cut");
    vz = 0.f; ptHat = 40.f;
    check(!isGoodEvent(true, true, 0, 0), "pPb must retain event filters");
    pBeamScrapingFilter = pPAprimaryVertexFilter = HBHENoiseFilterResultRun2Loose = 1;
    phfCoincFilter = pVertexFilterCutdz1p0 = 1;
    check(isGoodEvent(true, true, 0, 0), "pPb should accept the in-range event");
    ptHat = 1000.f;
    check(!isGoodEvent(true, true, 0, 0), "pPb must retain the pthat range cut");
    const TString prefix = Form("%s/HiForest_pp_branches_%d", gSystem->TempDirectory(),
                                gSystem->GetPid());
    const TString first = prefix + "_1.root", second = prefix + "_2.root";
    const TString list = prefix + ".list";
    writePpBranchFixture(first, 1.f);
    writePpBranchFixture(second, 2.f);
    {
        std::ofstream files(list.Data()); files << first << '\n' << second << '\n';
    }
    for (bool isPp : {true, false}) {
        TChain hlt("hltanalysis/HltTree"), events("hiEvtAnalyzer/HiTree");
        TChain skim("skimanalysis/HltTree"), jets("ak4PFJetAnalyzer/t");
        setupInput(list, hlt, events, skim, jets, true, isPp);
        hlt.GetMaximum("nref"); // Same pre-read used by processEvents.
        hlt.GetMaximum("ngen");
        for (int i = 0; i < 4; ++i) {
            check(hlt.GetEntry(i) > 0, "failed to read fixture entry");
            check(run == 1 && lumi == 2 && evt == static_cast<ULong64_t>((i / 2 + 1) * 100 + i % 2),
                  "EventID must update across entries and input files");
            check(eventRho == (i / 2 + 1) * 10 + i % 2,
                  "rho must update across entries and input files");
            check(ptHat == (i % 2 == 0 ? 1.f : 1000.f),
                  "pthat must be read from the event tree, not the jet tree");
            if (isPp) {
                check(inputEventWeight == (i % 2 == 0 ? 2.f : -0.5f) * (i / 2 + 1),
                      "weight must update for every event and across files");
            }
        }
    }
    gSystem->Unlink(first); gSystem->Unlink(second); gSystem->Unlink(list);
    std::cout << "pp/pPb weight and event-selection checks passed\n";
}
