#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlotVBS.C"

 //.x finalPlotVBS.C+(0,1,"m_{jj}","GeV","histo_nice8TeV.root","wwss_mjj",0,19.4)

void finalPlotVBS(int nsel = 0, int ReBin = 1, TString XTitle = "N_{jets}", TString units = "", 
                  TString plotName = "histo_nice.root", TString outputName = "njets",
                  bool isLogY = false, double lumi = 19.4) {
  // nsel == 0 --> WW EWK not overlaid with SM background
  // nsel == 1 --> WW EWK overlaid with SM background

  gInterpreter->ExecuteMacro("GoodStyle.C");
  TFile* file = new TFile(plotName.Data(), "read");

  StandardPlot myPlot;
  myPlot.setLumi(lumi);
  myPlot.setLabel(XTitle.Data());
  if     (lumi ==    4.9)  myPlot.addLabel("#sqrt{s} = 7 TeV");
  else if(lumi ==   19.4)  myPlot.addLabel("#sqrt{s} = 8 TeV");
  else if(lumi ==   24.4)  myPlot.addLabel("#sqrt{s} = 7+8 TeV");
  else                    myPlot.addLabel(""); 
  myPlot.setUnits(units.Data());

  TH1F* hWWEWK   = (TH1F*)file->Get("histo0");
  TH1F* hVV      = (TH1F*)file->Get("histo1");
  TH1F* hWWQCD   = (TH1F*)file->Get("histo2");
  TH1F* hWJets   = (TH1F*)file->Get("histo3");
  TH1F* hWW      = (TH1F*)file->Get("histo4");
  TH1F *hData    = (TH1F*)file->Get("histo5");

  double scale = 1;
  hWWEWK  ->Scale(scale);
  hVV	  ->Scale(scale);
  hWWQCD  ->Scale(scale);
  hWJets  ->Scale(scale);
  hWW	  ->Scale(scale);

  if(nsel == 0 || nsel == 1){
    myPlot.setMCHist(iWWEWK,(TH1F*)hWWEWK->Clone("hWWEWK"));
    myPlot._mass = 0;
    if(nsel == 1) myPlot.setHWWOverlaid(true);
    myPlot.setUnits(units);
    myPlot.setBreakdown(true);
  } else assert(0);

  myPlot.setMCHist(iVV,     (TH1F*)hVV     ->Clone("hVV"));
  myPlot.setMCHist(iWWQCD,  (TH1F*)hWWQCD  ->Clone("hWWQCD"));
  myPlot.setMCHist(iWJets,  (TH1F*)hWJets  ->Clone("hWJets")); 
  myPlot.setMCHist(iWW,     (TH1F*)hWW	   ->Clone("hWW"));
  myPlot.setDataHist((TH1F*)hData->Clone("data"));

  printf("%f + %f + %f + %f + %f = %f - %f\n",
          hWWEWK->GetSumOfWeights(),hVV->GetSumOfWeights(),hWWQCD->GetSumOfWeights(),
  	  hWJets->GetSumOfWeights(),hWW->GetSumOfWeights(),
	  hWWEWK->GetSumOfWeights()+hVV->GetSumOfWeights()+hWWQCD->GetSumOfWeights()+
	  hWJets->GetSumOfWeights()+hWW->GetSumOfWeights(),
	  hData->GetSumOfWeights());

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

  TCanvas* c2 = new TCanvas("c2", "c2",700,50,500,500);
  c2->cd(1);
  if(nsel == 0 || nsel == 10){
    TH1F* hSignal = (TH1F*)file->Get("histo0");
    TH1F* hSumBck = (TH1F*)file->Get("histo1");
    hSumBck->Add(hWWQCD  );
    hSumBck->Add(hWJets  );
    hSumBck->Add(hWW     );
    printf("S/B(%f/%f) = %f\n",hSignal->GetSumOfWeights(),hSumBck->GetSumOfWeights(),hSignal->GetSumOfWeights()/hSumBck->GetSumOfWeights());
    //hSignal->Divide(hSumBck);
    hSignal->Draw("e");
  }
}
