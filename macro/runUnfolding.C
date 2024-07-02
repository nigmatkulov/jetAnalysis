#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"

//________________
void runUnfolding() {
    gSystem->Load("~/work/RooUnfold/build/libRooUnfold.dylib");
    gROOT->LoadMacro("./unfoldDistributions.C");
    gROOT->ProcessLine( Form("unfoldDistributions()") );
    // gROOT->LoadMacro("./testUnfolding.C");
    // gROOT->ProcessLine("testUnfolding()");
}