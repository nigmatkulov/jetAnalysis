// Load JetCorrector.cc and JetUncertainty.cc first.
#include <fstream>
#include <stdexcept>
#include <cstdio>

void test_calibration_input_validation() {
    const std::string path = std::string(gSystem->TempDirectory()) +
        "/jet-calibration-validation-" + std::to_string(gSystem->GetPid()) + ".txt";
    auto mustReject = [&](auto make, const std::string &file) {
        try { make(file); }
        catch (const std::runtime_error &error) {
            if (std::string(error.what()).find(file) == std::string::npos)
                throw std::runtime_error("calibration error omitted the filename");
            return;
        }
        throw std::runtime_error("invalid calibration was accepted");
    };
    auto jec = [](const std::string &file) { SingleJetCorrector c(file); };
    auto jeu = [](const std::string &file) { JetUncertainty c(file); };
    mustReject(jec, path + ".missing");
    mustReject(jeu, path + ".missing");
    // MC intentionally has no uncertainty payload; preserve that sentinel.
    JetUncertainty unused("");
    if (unused.GetUncertainty() != std::make_pair(-1., -1.))
        throw std::runtime_error("empty MC uncertainty sentinel changed");
    {
        std::ofstream out(path);
        out << std::string(1048576, '1') << '\n';
    }
    mustReject(jec, path); mustReject(jeu, path);
    {
        std::ofstream out(path);
        out << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}\n";
        out << "-5 5 6 10 0.1 0.2\n"; // six values declared, only three provided
    }
    mustReject(jeu, path);
    std::remove(path.c_str());
    JetCorrector valid("aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_Unembedded_MC_L2Relative_AK4PF.txt");
    valid.SetJetPT(100.); valid.SetJetEta(0.); valid.SetJetPhi(0.);
    if (!(valid.GetCorrectedPT() > 100.)) throw std::runtime_error("valid JEC failed");
    std::cout << "calibration input: missing files, oversized lines, truncated triples and valid payload passed\n";
}
