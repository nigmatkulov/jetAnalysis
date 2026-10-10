// Load JetCorrector.cc first. This checks real payload loading/evaluation;
// executable logs separately check which payload the main program selects.
#include <stdexcept>
void test_jec_payloads() {
    const std::string base = "aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_";
    JetCorrector embedded(base + "Embedded_MC_L2Relative_AK4PF.txt");
    JetCorrector unembedded(base + "Unembedded_MC_L2Relative_AK4PF.txt");
    for (double eta : {-1.5, 0., 1.5}) {
        embedded.SetJetPT(100.); embedded.SetJetEta(eta); embedded.SetJetPhi(0.);
        unembedded.SetJetPT(100.); unembedded.SetJetEta(eta); unembedded.SetJetPhi(0.);
        const double a = embedded.GetCorrectedPT(), b = unembedded.GetCorrectedPT();
        if (!std::isfinite(a) || !std::isfinite(b) || a <= 0. || b <= 0. || a == b)
            throw std::runtime_error("Pb-going JEC payloads did not load as distinct valid corrections");
        std::cout << "Pb-going rawPt=100 eta=" << eta << " Embedded=" << a
                  << " Unembedded=" << b << std::endl;
    }
}
