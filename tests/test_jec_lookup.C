// Load JetCorrector.cc first, then run test_jec_lookup().
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <stdexcept>

struct JetCorrectorLookupTest {
    static void compare(SingleJetCorrector &corrector, double eta, double pt) {
        corrector.SetJetEta(eta);
        corrector.SetJetPT(pt);
        corrector.SetRho(17.);
        corrector.SetJetPhi(.3);
        corrector.SetJetArea(.5);
        const double indexed = corrector.GetCorrection();
        // An empty index selects the original full scan, with the exact same
        // payload, dependency clamping, TF1 evaluation, and first-match rule.
        std::vector<std::vector<int>> saved;
        saved.swap(corrector.LookupRows);
        const double scanned = corrector.GetCorrection();
        saved.swap(corrector.LookupRows);
        if (indexed != scanned && !(std::isnan(indexed) && std::isnan(scanned)))
            throw std::runtime_error("JEC indexed/linear mismatch at eta=" +
                std::to_string(eta) + " pt=" + std::to_string(pt));
    }
    static void payload(const std::string &file) {
        SingleJetCorrector corrector(file);
        if (corrector.LookupRows.empty()) throw std::runtime_error("Real payload was not indexed: " + file);
        size_t checks = 0;
        for (size_t row = 0; row < corrector.BinTypes.size(); ++row) {
            std::vector<double> etas{.3}, pts{75.};
            for (size_t axis = 0; axis < corrector.BinTypes[row].size(); ++axis) {
                auto *values = corrector.BinTypes[row][axis] == SingleJetCorrector::TypeJetEta ? &etas
                    : corrector.BinTypes[row][axis] == SingleJetCorrector::TypeJetPT ? &pts : nullptr;
                if (!values) continue;
                values->clear();
                const double low = corrector.BinRanges[row][axis * 2];
                const double high = corrector.BinRanges[row][axis * 2 + 1];
                for (double edge : {low, high}) {
                    values->push_back(std::nextafter(edge, -std::numeric_limits<double>::infinity()));
                    values->push_back(edge);
                    values->push_back(std::nextafter(edge, std::numeric_limits<double>::infinity()));
                }
                values->push_back((low + high) / 2.);
            }
            for (double eta : etas) for (double pt : pts) { compare(corrector, eta, pt); ++checks; }
        }
        std::mt19937 random(12345);
        std::uniform_real_distribution<double> eta(-6., 6.), pt(-100., 7000.);
        for (int i = 0; i < 2000; ++i) { compare(corrector, eta(random), pt(random)); ++checks; }
        compare(corrector, std::numeric_limits<double>::quiet_NaN(), 75.);
        compare(corrector, .3, std::numeric_limits<double>::infinity());
        std::cout << file << ": " << checks + 2 << " exact comparisons passed\n";
    }
};

void test_jec_lookup() {
    const std::string base = "aux_files/pPb_8160/JEC/";
    for (const char *name : {
        "Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt",
        "Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt",
        "Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF.txt",
        "Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PF.txt"})
        JetCorrectorLookupTest::payload(base + name);
    // Overlaps, reversed bounds, gaps, and differently valued adjacent rows
    // make first-match ordering and boundary errors observable.
    const std::string file = std::string(gSystem->TempDirectory()) +
        "/jec-lookup-" + std::to_string(gSystem->GetPid()) + ".txt";
    {
        std::ofstream out(file);
        out << "{2 JetEta JetPt 1 JetPt [0]+[1]*x Correction Test}\n"
            << "-1 0 10 20 4 12 18 1 .01\n"
            << "0 1 10 20 4 12 18 2 .02\n"
            << "1 -1 30 15 4 12 18 3 .03\n"
            << "2 3 50 60 4 12 18 4 .04\n";
    }
    JetCorrectorLookupTest::payload(file);
    // Reinitialization appends in the original implementation. Its index must
    // be rebuilt and retain the first payload's precedence.
    SingleJetCorrector repeated(file);
    repeated.Initialize(file);
    JetCorrectorLookupTest::compare(repeated, 0., 20.);
    {
        std::ofstream out(file);
        out << "{1 JetEta 1 JetPt [0] Correction Test}\n"
            << "-1 0 3 10 20 2\n"
            << "{1 JetPt 1 JetPt [0] Correction Test}\n"
            << "10 20 3 10 20 3\n";
    }
    SingleJetCorrector heterogeneous(file);
    JetCorrectorLookupTest::compare(heterogeneous, 0., 15.);
    JetCorrectorLookupTest::compare(heterogeneous, .5, 15.);
    // Every level must see pT after the preceding correction, not raw pT.
    {
        std::ofstream out(file);
        out << "{1 JetEta 1 JetPt [0]*x Correction Test}\n"
            << "-5 5 3 0 1000 .01\n";
    }
    JetCorrector levels(std::vector<std::string>{file, file});
    levels.SetJetPT(50.);
    levels.SetJetEta(.3);
    if (levels.GetCorrectedPT() != 6.25)
        throw std::runtime_error("Sequential correction levels changed");
    std::remove(file.c_str());
    std::cout << "JEC lookup regression passed\n";
}
