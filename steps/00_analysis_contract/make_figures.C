#include "../../AnalysisContract.h"

#include "TArrow.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"

void make_figures() {
    gSystem->mkdir("steps/00_analysis_contract/figures", true);

    constexpr int n = 41;
    double stored[n], mcPb[n], mcP[n], dataPb[n], dataP[n];
    for (int i = 0; i < n; ++i) {
        stored[i] = -2. + 4. * i / (n - 1.);
        const double shift = analysis_contract::kNominalEtaShift;
        mcPb[i] = -(stored[i] + shift);
        mcP[i] = stored[i] - shift;
        dataPb[i] = stored[i] - shift;
        dataP[i] = -(stored[i] + shift);
    }

    TCanvas overlay("overlay", "orientation transformations", 900, 700);
    TGraph graphs[4] = {TGraph(n, stored, mcPb), TGraph(n, stored, mcP),
                        TGraph(n, stored, dataPb), TGraph(n, stored, dataP)};
    const int colors[4] = {2, 4, 6, 8};
    const char *labels[4] = {"Pb-going MC", "p-going MC", "Pb-going data", "p-going data"};
    graphs[0].SetTitle("Step 00 orientation contract;stored jet #eta;transformed #eta_{CM}");
    TLegend legend(0.14, 0.68, 0.43, 0.88);
    for (int i = 0; i < 4; ++i) {
        graphs[i].SetLineColor(colors[i]);
        graphs[i].SetLineWidth(3);
        graphs[i].Draw(i == 0 ? "AL" : "L SAME");
        legend.AddEntry(&graphs[i], labels[i], "l");
    }
    legend.Draw();
    overlay.SaveAs("steps/00_analysis_contract/figures/direction_overlay.png");
    overlay.SaveAs("steps/00_analysis_contract/figures/direction_overlay.pdf");

    TCanvas event("event", "signed eta example", 900, 450);
    event.DrawFrame(-2.5, -1., 2.5, 1., "Signed-#eta example;#eta;beam / jet annotation");
    TLine axis(-2.3, 0., 2.3, 0.); axis.SetLineWidth(2); axis.Draw();
    TArrow proton(0., 0.35, 2.0, 0.35, 0.03, "|>"); proton.SetLineColor(4); proton.Draw();
    TArrow lead(0., -0.35, 1.2, -0.35, 0.03, "|>"); lead.SetLineColor(2); lead.Draw();
    TLatex text; text.SetTextSize(0.045);
    text.DrawLatex(1.45, 0.5, "proton-going: +#eta_{CM}");
    text.DrawLatex(0.65, -0.55, "example selected jet");
    text.DrawLatex(-2.25, 0.78, "All four storage branches are transformed to this common orientation");
    event.SaveAs("steps/00_analysis_contract/figures/signed_eta_event.png");
    event.SaveAs("steps/00_analysis_contract/figures/signed_eta_event.pdf");
}
