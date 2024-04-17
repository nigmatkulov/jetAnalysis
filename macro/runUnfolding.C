#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"

//________________
void runUnfolding() {
    gSystem->Load("~/work/RooUnfold/build/libRooUnfold.dylib");
    gROOT->LoadMacro("./unfoldDistributions.C");
    gROOT->ProcessLine("unfoldDistributions()");
}