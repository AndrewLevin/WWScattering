#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlotVBS.C"

 //.x finalPlotVBS.C+(0,1,"m_{jj} (GeV)","GeV","histo_nice8TeV.root","wwss_mjj",0,19.4)

void finalPlotVBS(int nsel = 0, int ReBin = 1, TString XTitle = "N_{jets}", TString units = "", 
                  TString plotName = "histo_nice.root", TString outputName = "njets",
                  bool isLogY = false, double lumi = 19.4, bool blindedData = false) {
  // nsel == 0 --> WW EWK not overlaid with SM background
  // nsel == 1 --> WW EWK overlaid with SM background

  gInterpreter->ExecuteMacro("GoodStyle.C");
  TFile* file = new TFile(plotName.Data(), "read");

  StandardPlot myPlot;
  myPlot.setLabel(XTitle.Data());
  if     (lumi ==    4.9) myPlot.setLumiLabel("4.9 fb^{-1} (7 TeV)");
  else if(lumi ==   19.4) myPlot.setLumiLabel("19.4 fb^{-1} (8 TeV)");
  else if(lumi ==   24.4) myPlot.setLumiLabel("4.9 fb^{-1} (7 TeV) + 19.4 fb^{-1} (8 TeV)");
  else                    myPlot.setLumiLabel(""); 
  myPlot.setUnits(units.Data());

  TH1D* hWWEWK   = (TH1D*)file->Get("histo0");
  TH1D* hVV      = (TH1D*)file->Get("histo1");
  TH1D* hWWQCD   = (TH1D*)file->Get("histo2");
  TH1D* hWJets   = (TH1D*)file->Get("histo3");
  TH1D* hWW      = (TH1D*)file->Get("histo4");
  TH1D* hVVV     = (TH1D*)file->Get("histo5");
  TH1D* hHiggs   = (TH1D*)file->Get("histo7");
  if(!hHiggs) {hHiggs   = (TH1D*)file->Get("histo0"); hHiggs->Scale(0);}
  TH1D *hData    = (TH1D*)file->Get("histo6");
  if(blindedData == true) hData->Scale(0);

  printf("%f + %f + %f + %f + %f + %f = %f - %f - %f\n",
          hWWEWK->GetSumOfWeights(),hVV->GetSumOfWeights(),hWWQCD->GetSumOfWeights(),
  	  hWJets->GetSumOfWeights(),hWW->GetSumOfWeights(),hVVV->GetSumOfWeights(),
	  hWWEWK->GetSumOfWeights()+hVV->GetSumOfWeights()+hWWQCD->GetSumOfWeights()+
	  hWJets->GetSumOfWeights()+hWW->GetSumOfWeights()+hVVV->GetSumOfWeights(),
	  hData->GetSumOfWeights(),hHiggs->GetSumOfWeights());

  double SFHiggs = 1.0;
  if(hHiggs->GetSumOfWeights() > 0) {
    SFHiggs = (hWWEWK->GetSumOfWeights()+hVV->GetSumOfWeights()+hWWQCD->GetSumOfWeights()+
               hWJets->GetSumOfWeights()+hWW->GetSumOfWeights()+hVVV  ->GetSumOfWeights())/hHiggs->GetSumOfWeights();
  }

  double scale = 1;
  hWWEWK  ->Scale(scale);
  hVV	  ->Scale(scale);
  hWWQCD  ->Scale(scale);
  hWJets  ->Scale(scale);
  hWW	  ->Scale(scale);
  hVVV	  ->Scale(scale);
  hHiggs  ->Scale(scale*SFHiggs);

  if(nsel == -1 || nsel == 0 || nsel == 1 || nsel == 2 || nsel == 3){
    myPlot.setMCHist(iWWEWK,(TH1D*)hWWEWK->Clone("hWWEWK"));
    myPlot._mass = 0;
    if(nsel >= 1 || nsel == -1 || hHiggs->GetSumOfWeights() > 0) myPlot.setHWWOverlaid(true);
    myPlot.setUnits(units);
    myPlot.setBreakdown(true);
    if(nsel == 1) myPlot.setMass(200);
    else          myPlot.setMass(800);
    
    if(nsel == 3) myPlot.setAlternativeOption(3);
    
  } else assert(0);

  TH1D* hOther = (TH1D*) hWWQCD->Clone("hOther");
  hOther->Add(hVVV);
  hOther->Add(hWW);

  myPlot.setMCHist(iVV,     (TH1D*)hVV     ->Clone("hVV"));
  if(nsel == -1) {
    myPlot.setMCHist(iOther,  (TH1D*)hOther  ->Clone("hOther"));
  } else {
    myPlot.setMCHist(iWWQCD,  (TH1D*)hWWQCD  ->Clone("hWWQCD"));
    myPlot.setMCHist(iWW,     (TH1D*)hWW     ->Clone("hWW"));
    myPlot.setMCHist(iVVV,    (TH1D*)hVVV    ->Clone("hVVV"));
  }
  myPlot.setMCHist(iWJets,  (TH1D*)hWJets  ->Clone("hWJets")); 
  myPlot.setMCHist(iHiggs,  (TH1D*)hHiggs  ->Clone("hHiggs"));
  myPlot.setDataHist((TH1D*)hData->Clone("data"));

  TCanvas* c1 = new TCanvas("c1", "c1",5,50,500,500);

  if(isLogY == true) c1->SetLogy();
  myPlot.Draw(ReBin);  // Can pass a rebin 
  c1->GetFrame()->DrawClone();

  char CommandToExec[300];
  sprintf(CommandToExec,"mkdir -p plots");
  gSystem->Exec(CommandToExec);  

  char myOutputFile[300];
  sprintf(myOutputFile,"plots/%s.eps",outputName.Data());
  c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"plots/%s.png",outputName.Data());
  c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"plots/%s.pdf",outputName.Data());
  c1->SaveAs(myOutputFile);

  if(nsel == 0 || nsel == 10){
    TCanvas* c2 = new TCanvas("c2", "c2",700,50,500,500);
    c2->cd(1);
    TH1D* hSignal = (TH1D*)file->Get("histo0");
    TH1D* hSumBck = (TH1D*)file->Get("histo1");
    hSumBck->Add(hWWQCD  );
    hSumBck->Add(hWJets  );
    hSumBck->Add(hWW     );
    hSumBck->Add(hVVV    );
    hSignal->Rebin(ReBin);
    hSumBck->Rebin(ReBin);
    printf("S/B(%f/%f) = %f\n",hSignal->GetSumOfWeights(),hSumBck->GetSumOfWeights(),hSignal->GetSumOfWeights()/hSumBck->GetSumOfWeights());
    hSignal->Divide(hSumBck);
    hSignal->Draw("e");
  }
}
