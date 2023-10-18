/**
  Macro reads input file reads and saves all TH1 and TH2 derived objects
  as .png and .C in output directory
**/

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TKey.h"
#include "TList.h"

#include <iostream>

//________________
// Set common style
void setStyle() {
  TStyle* myStyle = new TStyle("myStyle","My style");
  myStyle->SetOptStat(0);
  myStyle->SetOptTitle(1);
  myStyle->SetOptDate(0);
  myStyle->SetLabelSize(0.03,"xyz");   // size of axis value font
  myStyle->SetTitleSize(0.035,"xyz");  // size of axis title font
  myStyle->SetTitleFont(22,"xyz");     // font option
  myStyle->SetLabelFont(22, "xyz");
  myStyle->SetTitleOffset(1.2, "y");
  //myStyle->SetFillColor(0);            // Does some weird things with 2D figures
  myStyle->SetTitleColor(1);
  myStyle->SetLineWidth(2);
  myStyle->SetMarkerStyle(20);
  myStyle->SetMarkerSize(1.4);
  myStyle->SetFrameLineWidth(2);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleColor(1);
  myStyle->SetStatColor(0);
  myStyle->SetTitleTextColor(1);

  gROOT->SetStyle("myStyle");
  gROOT->ForceStyle();
}

//________________
void extractFigures(const Char_t *inFile, const Char_t *outputPath = "figs/",
                    bool savePNG = true, bool saveMacro = false) {
   setStyle();
   TFile *f1 = TFile::Open(inFile);
   TList *list = f1->GetListOfKeys();
   TKey *key;
   TCanvas c1;
   TString hname, tmpName, className;
   TString textForLogo("CMS Pb+Pb #sqrt{s_{NN}} = 5020 GeV");
   TIter iter(list);
   static TString h1dClassName("TH1D");
   static TString h2fClassName("TH2F");
   static TString h2dClassName("TH2D");
   while ((key = (TKey*)iter())) {
      std::cout << "Class name: " << key->GetClassName() << std::endl;
      if ( key->GetClassName() == h1dClassName ) {
        TH1D *h = (TH1D*)key->ReadObj();
        hname = outputPath;
        hname += h->GetName();
        tmpName = h->GetName();
        tmpName.ToLower();
        h->Draw();
        gPad->SetLogy(0);
        if ( tmpName.Contains("refmult") ||
             tmpName.Contains("chemistry") ) {
          gPad->SetLogy(1);
        }
        TLatex logo;
        logo.DrawLatexNDC(0.15, 0.85, textForLogo.Data());
        c1.Update();
        if (savePNG) {
          c1.Print(Form("%s.png",hname.Data()));
        }
        if (saveMacro) {
          c1.Print(Form("%s.C",hname.Data()));
        }
      }
      else if ( key->GetClassName() == h2fClassName ) {
        TH2F *h = (TH2F*)key->ReadObj();
        hname = outputPath;
        hname += h->GetName();
        h->Draw("colz");
        gPad->SetLogy(0);
        gPad->SetLogz(1);
        TLatex logo;
        logo.DrawLatexNDC(0.15, 0.85, textForLogo.Data());
        c1.Update();
        if (savePNG) {
          c1.Print(Form("%s.png",hname.Data()));
        }
        if (saveMacro) {
          c1.Print(Form("%s.C",hname.Data()));
        }
      }
      else if ( key->GetClassName() == h2dClassName ) {
        TH2D *h = (TH2D*)key->ReadObj();
        hname = outputPath;
        hname += h->GetName();
        h->Draw("colz");
        gPad->SetLogy(0);
        gPad->SetLogz(1);
        TLatex logo;
        logo.DrawLatexNDC(0.15, 0.85, textForLogo.Data());
        c1.Update();
        if (savePNG) {
          c1.Print(Form("%s.png",hname.Data()));
        }
        if (saveMacro) {
          c1.Print(Form("%s.C",hname.Data()));
        }
      }
   }
}