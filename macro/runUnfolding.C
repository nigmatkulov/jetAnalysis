#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"

//________________
void runUnfolding() {
    gSystem->Load("~/work/RooUnfold/build/libRooUnfold.dylib");
    gROOT->LoadMacro("./unfoldDistributions.C");
    const Char_t *date = "20240419";
    gROOT->ProcessLine( Form("unfoldDistributions(\"%s\")", date) );
}