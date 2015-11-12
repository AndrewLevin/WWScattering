#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/makeSystematicEffects3l.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/LeptonEfficiencyZH.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLorentzVector.h"

const int verboseLevel =   1;
bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC
SmurfTree systEvent;
const unsigned int nSelTypes = 2;
const unsigned int nSelTypesSyst = 9;
const bool showSignalOnly = false;

enum selType {BTAGSEL, WZSEL};
TString selTypeName[nSelTypes*2] = {"BTAGSEL-OS", "WZSEL-OS",
                                    "BTAGSEL-SS", "WZSEL-SS"};
enum selTypeSyst {JESUP=0, JESDOWN, LEPP, LEPM, MET, JERUP, JERDOWN, EFFP, EFFM};
TString selTypeNameSyst[nSelTypesSyst*2] = {"JESUP-OS", "JESDOWN-OS", "LEPP-OS", "LEPM-OS", "MET-OS", "JERUP-OS", "JERDOWN-OS", "EFFP-OS", "EFFM-OS",
                                            "JESUP-SS", "JESDOWN-SS", "LEPP-SS", "LEPM-SS", "MET-SS", "JERUP-SS", "JERDOWN-SS", "EFFP-SS", "EFFM-SS"};

bool use_fake_rate_method = true;
bool run_over_data = true;

void scaleFactor_WS(LorentzVector l,int q, int ld, int mcld, double val[2], int opt);

// thePlot == 0 (mjj), 2 (ptlmax), 9 (mll), 19(mt), anything else (mlljj)

void wz_ana
(
 int thePlot = 0,
 int lSel = 4,
 TString bgdInputFile    = "ntuples_53x/backgroundEWK_skim14_sm.root",
 TString dataInputFile   = "ntuples_53x/data_skim14.root",
 int period = 3
 )
{

  double genLevelNorm[8] = {0.,0.,0.,0.,0.,0.,0.,0.};

  double frCorr = 0.78;
  double lumi = 1.0;
  double ptJetMin = 30.0;

  bool fCheckProblem = true;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  TString ECMsb  = "";
  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  double lumiE = 1.099; int year = 1;
  if	 (period == 3){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365;minRun =      0;maxRun = 999999;ECMsb="8TeV";lumiE = 1.026; year = 2012;
  }
  else if(period == 4){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/puWeights_Fall11_42x_True.root";
    lumi     = 4.924;minRun =	 0;maxRun = 999999;ECMsb="7TeV"; lumiE = 1.022; year = 2011;
    UseDyttDataDriven = false;fCheckProblem = false;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  //----------------------------------------------------------------------------
  // radio photon to electron
  //----------------------------------------------------------------------------
  TFile *fRatioPhotonElectron = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/ratio_photon_electron.root");
  TH1D *fhDRatioPhotonElectron = (TH1D*)(fRatioPhotonElectron->Get("hDRatioPhotonElectron"));
  assert(fhDRatioPhotonElectron);
  fhDRatioPhotonElectron->SetDirectory(0);
  fRatioPhotonElectron->Close();
  delete fRatioPhotonElectron;

  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open(fakePath.Data());
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open(fakePath.Data());
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  //Fake rate systematics
  TFile *fLeptonFRFileSyst = TFile::Open(fakePath.Data());
  TH2D *fhDFRMuSyst = (TH2D*)(fLeptonFRFileSyst->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  TH2D *fhDFRElSyst = (TH2D*)(fLeptonFRFileSyst->Get("ElectronFakeRate_V4_ptThreshold20_PtEta"));
  assert(fhDFRMuSyst);
  assert(fhDFRElSyst);
  fhDFRMuSyst->SetDirectory(0);
  fhDFRElSyst->SetDirectory(0);
  fLeptonFRFileSyst->Close();
  delete fLeptonFRFileSyst;
 
  LeptonScaleLookup trigLookup(effPath.Data());

  // useful if using ZZ lepton selection
  LeptonEfficiencyZH theLeptonEfficiencyZH(year);

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  double JERVal[5]  = {1.052,1.057,1.096,1.134,1.288}; // 0.0.0.5, 0.5-1.1, 1.1-1.7, 1.7-2.3, 2.3-5.0
  double JERDown[5] = {0.990,1.001,1.032,1.042,1.089}; 
  double JERUp[5]   = {1.115,1.114,1.161,1.228,1.488};

  const int nBin = 4;
  Float_t xbins[nBin+1] = {700, 1100, 1500, 2000, 3000};
  if     (thePlot == 0) {xbins[0] = 500; xbins[1] = 700; xbins[2] = 1100; xbins[3] = 1600; xbins[4] = 2000;}
  else if(thePlot == 2) {xbins[0] =   0; xbins[1] = 100; xbins[2] =  200; xbins[3] =  300; xbins[4] =  500;}
  else if(thePlot == 9) {xbins[0] =  50; xbins[1] = 100; xbins[2] =  200; xbins[3] =  300; xbins[4] =  500;}
  else if(thePlot ==19) {xbins[0] = 0.0; xbins[1] = 150; xbins[2] =  250; xbins[3] =  350; xbins[4] = 1000;}
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBin, xbins);
  histoMVA->Sumw2();
  TH1D *histo_Data      = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_WZ        = (TH1D*) histoMVA->Clone("histo_WZ");
  TH1D *histo_ZZ        = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_VVV       = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_Wjets     = (TH1D*) histoMVA->Clone("histo_Wjets");
  TH1D *histo_Higgs     = (TH1D*) histoMVA->Clone("histo_Higgs");

  char finalStateName[2],effName[10],momName[10];sprintf(effName,"CMS_eff_l");sprintf(momName,"CMS_p_scale_l");
  if     (lSel == 0) {sprintf(finalStateName,"mm");}
  else if(lSel == 1) {sprintf(finalStateName,"me");}
  else if(lSel == 2) {sprintf(finalStateName,"em");}
  else if(lSel == 3) {sprintf(finalStateName,"ee");}
  else if(lSel == 4) {sprintf(finalStateName,"ll");}
  else if(lSel == 5) {sprintf(finalStateName,"sf");}
  else if(lSel == 6) {sprintf(finalStateName,"of");}
  else {printf("Wrong lSel: %d\n",lSel); assert(0);}

  int nBinPlot      = 400;
  double xminPlot   = 0.0;
  double xmaxPlot   = 400.0;

  if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 8; xminPlot = 0.0; xmaxPlot =  2000.0;} // mjj
  else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 20; xminPlot = 0.0; xmaxPlot = 2000.0;} // mlljj
  else if(thePlot >=  2 && thePlot <=  8) {}
  else if(thePlot >=  9 && thePlot <=  9) {nBinPlot = 50;  xminPlot =  0.0; xmaxPlot =  500.0;} // mll
  else if(thePlot >= 10 && thePlot <= 10) {nBinPlot = 10;  xminPlot =  -0.5; xmaxPlot =  9.5;}
  else if(thePlot >= 11 && thePlot <= 11) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot >= 12 && thePlot <= 12) {nBinPlot = 36; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 13 && thePlot <= 14) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 15 && thePlot <= 15) {nBinPlot = 7; xminPlot = 0.0; xmaxPlot =  8.75;} // detajjs
  else if(thePlot >= 16 && thePlot <= 16) {nBinPlot = 4; xminPlot = -0.5; xmaxPlot = 3.5;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 44; xminPlot = 0.0; xmaxPlot = 4.4;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 40; xminPlot = 0.0; xmaxPlot = 4.0;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 50; xminPlot =  0.0; xmaxPlot = 1000.0;}
  else if(thePlot >= 20 && thePlot <= 20) {nBinPlot =100; xminPlot =  0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 21 && thePlot <= 21) {nBinPlot =40; xminPlot =  0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 22 && thePlot <= 22) {nBinPlot =50; xminPlot =  0.0; xmaxPlot = 5.0;}
  else assert(0);

  TH1D* histo0;
  if(thePlot != 0 && thePlot != 1 && thePlot != 2 && thePlot != 9 && thePlot != 19) histo0 = new TH1D("histo0", "histo0", nBinPlot, xminPlot, xmaxPlot);
  else                                                                              histo0 = new TH1D("histo0", "histo0", nBin, xbins);  
  histo0->Sumw2();
  TH1D* histo1 = (TH1D*) histo0->Clone("histo1");
  TH1D* histo2 = (TH1D*) histo0->Clone("histo2");
  TH1D* histo3 = (TH1D*) histo0->Clone("histo3");
  TH1D* histo4 = (TH1D*) histo0->Clone("histo4");
  TH1D* histo5 = (TH1D*) histo0->Clone("histo5");
  TH1D* histo6 = (TH1D*) histo0->Clone("histo6");
  TH1D* histo7 = (TH1D*) histo0->Clone("histo7");
  histo0->Scale(0.0);
  histo1->Scale(0.0);
  histo2->Scale(0.0);
  histo3->Scale(0.0);
  histo4->Scale(0.0);
  histo5->Scale(0.0);
  histo6->Scale(0.0);
  histo7->Scale(0.0);

  TH1D* histo_ZZ_ZZStatUp         = new TH1D( Form("histo_ZZ_CMS_qqwz%s_MVAZZStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_qqwz%s_MVAZZStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_ZZ_ZZStatUp  ->Sumw2();
  TH1D* histo_ZZ_ZZStatDown       = new TH1D( Form("histo_ZZ_CMS_qqwz%s_MVAZZStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_qqwz%s_MVAZZStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_ZZ_ZZStatDown->Sumw2();
  TH1D* histo_WZ_WZStatUp         = new TH1D( Form("histo_WZ_CMS_qqwz%s_MVAWZStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_qqwz%s_MVAWZStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_WZ_WZStatUp  ->Sumw2();
  TH1D* histo_WZ_WZStatDown       = new TH1D( Form("histo_WZ_CMS_qqwz%s_MVAWZStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_qqwz%s_MVAWZStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_WZ_WZStatDown->Sumw2();
  TH1D* histo_VVV_VVVStatUp       = new TH1D( Form("histo_VVV_CMS_qqwz%s_MVAVVVStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_qqwz%s_MVAVVVStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_VVV_VVVStatUp  ->Sumw2();
  TH1D* histo_VVV_VVVStatDown     = new TH1D( Form("histo_VVV_CMS_qqwz%s_MVAVVVStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_qqwz%s_MVAVVVStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_VVV_VVVStatDown->Sumw2();
  TH1D* histo_Wjets_WjetsStatUp   = new TH1D( Form("histo_Wjets_CMS_qqwz%s_MVAWjetsStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_qqwz%s_MVAWjetsStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_Wjets_WjetsStatUp  ->Sumw2();
  TH1D* histo_Wjets_WjetsStatDown = new TH1D( Form("histo_Wjets_CMS_qqwz%s_MVAWjetsStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_qqwz%s_MVAWjetsStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_Wjets_WjetsStatDown->Sumw2();
  TH1D* histo_Higgs_HiggsStatUp   = new TH1D( Form("histo_Higgs_CMS_qqwz%s__MVAHiggsStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_qqwz%s__MVAHiggsStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_Higgs_HiggsStatUp  ->Sumw2();
  TH1D* histo_Higgs_HiggsStatDown = new TH1D( Form("histo_Higgs_CMS_qqwz%s_MVAHiggsStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_qqwz%s__MVAHiggsStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_Higgs_HiggsStatDown->Sumw2();

  TH1D* histo_WZ_LepEffUp         = new TH1D( Form("histo_WZ_%sUp",effName)  , Form("histo_WZ_%sUp",effName)  , nBin, xbins); histo_WZ_LepEffUp  ->Sumw2();
  TH1D* histo_WZ_LepEffDown       = new TH1D( Form("histo_WZ_%sDown",effName), Form("histo_WZ_%sDown",effName), nBin, xbins); histo_WZ_LepEffDown->Sumw2();
  TH1D* histo_ZZ_LepEffUp         = new TH1D( Form("histo_ZZ_%sUp",effName)  , Form("histo_ZZ_%sUp",effName)  , nBin, xbins); histo_ZZ_LepEffUp  ->Sumw2();
  TH1D* histo_ZZ_LepEffDown       = new TH1D( Form("histo_ZZ_%sDown",effName), Form("histo_ZZ_%sDown",effName), nBin, xbins); histo_ZZ_LepEffDown->Sumw2();
  TH1D* histo_VVV_LepEffUp        = new TH1D( Form("histo_VVV_%sUp",effName)  , Form("histo_VVV_%sUp",effName)  , nBin, xbins); histo_VVV_LepEffUp  ->Sumw2();
  TH1D* histo_VVV_LepEffDown      = new TH1D( Form("histo_VVV_%sDown",effName), Form("histo_VVV_%sDown",effName), nBin, xbins); histo_VVV_LepEffDown->Sumw2();
  TH1D* histo_Higgs_LepEffUp      = new TH1D( Form("histo_Higgs_%sUp",effName)  , Form("histo_Higgs_%sUp",effName)  , nBin, xbins); histo_Higgs_LepEffUp  ->Sumw2();
  TH1D* histo_Higgs_LepEffDown    = new TH1D( Form("histo_Higgs_%sDown",effName), Form("histo_Higgs_%sDown",effName), nBin, xbins); histo_Higgs_LepEffDown->Sumw2();

  TH1D* histo_WZ_LepResUp         = new TH1D( Form("histo_WZ_%sUp",momName)  , Form("histo_WZ_%sUp",momName)  , nBin, xbins); histo_WZ_LepResUp  ->Sumw2();
  TH1D* histo_WZ_LepResDown       = new TH1D( Form("histo_WZ_%sDown",momName), Form("histo_WZ_%sDown",momName), nBin, xbins); histo_WZ_LepResDown->Sumw2();
  TH1D* histo_ZZ_LepResUp         = new TH1D( Form("histo_ZZ_%sUp",momName)  , Form("histo_ZZ_%sUp",momName)  , nBin, xbins); histo_ZZ_LepResUp  ->Sumw2();
  TH1D* histo_ZZ_LepResDown       = new TH1D( Form("histo_ZZ_%sDown",momName), Form("histo_ZZ_%sDown",momName), nBin, xbins); histo_ZZ_LepResDown->Sumw2();
  TH1D* histo_VVV_LepResUp        = new TH1D( Form("histo_VVV_%sUp",momName)  , Form("histo_VVV_%sUp",momName)  , nBin, xbins); histo_VVV_LepResUp  ->Sumw2();
  TH1D* histo_VVV_LepResDown      = new TH1D( Form("histo_VVV_%sDown",momName), Form("histo_VVV_%sDown",momName), nBin, xbins); histo_VVV_LepResDown->Sumw2();
  TH1D* histo_Higgs_LepResUp      = new TH1D( Form("histo_Higgs_%sUp",momName)  , Form("histo_Higgs_%sUp",momName)  , nBin, xbins); histo_Higgs_LepResUp  ->Sumw2();
  TH1D* histo_Higgs_LepResDown    = new TH1D( Form("histo_Higgs_%sDown",momName), Form("histo_Higgs_%sDown",momName), nBin, xbins); histo_Higgs_LepResDown->Sumw2();

  TH1D* histo_WZ_METResUp         = new TH1D( Form("histo_WZ_%sUp","CMS_scale_met")  , Form("histo_WZ_%sUp","CMS_scale_met")  , nBin, xbins); histo_WZ_METResUp  ->Sumw2();
  TH1D* histo_WZ_METResDown       = new TH1D( Form("histo_WZ_%sDown","CMS_scale_met"), Form("histo_WZ_%sDown","CMS_scale_met"), nBin, xbins); histo_WZ_METResDown->Sumw2();
  TH1D* histo_ZZ_METResUp         = new TH1D( Form("histo_ZZ_%sUp","CMS_scale_met")  , Form("histo_ZZ_%sUp","CMS_scale_met")  , nBin, xbins); histo_ZZ_METResUp  ->Sumw2();
  TH1D* histo_ZZ_METResDown       = new TH1D( Form("histo_ZZ_%sDown","CMS_scale_met"), Form("histo_ZZ_%sDown","CMS_scale_met"), nBin, xbins); histo_ZZ_METResDown->Sumw2();
  TH1D* histo_VVV_METResUp        = new TH1D( Form("histo_VVV_%sUp","CMS_scale_met")  , Form("histo_VVV_%sUp","CMS_scale_met")  , nBin, xbins); histo_VVV_METResUp  ->Sumw2();
  TH1D* histo_VVV_METResDown      = new TH1D( Form("histo_VVV_%sDown","CMS_scale_met"), Form("histo_VVV_%sDown","CMS_scale_met"), nBin, xbins); histo_VVV_METResDown->Sumw2();
  TH1D* histo_Higgs_METResUp      = new TH1D( Form("histo_Higgs_%sUp","CMS_scale_met")  , Form("histo_Higgs_%sUp","CMS_scale_met")  , nBin, xbins); histo_Higgs_METResUp  ->Sumw2();
  TH1D* histo_Higgs_METResDown    = new TH1D( Form("histo_Higgs_%sDown","CMS_scale_met"), Form("histo_Higgs_%sDown","CMS_scale_met"), nBin, xbins); histo_Higgs_METResDown->Sumw2();

  TH1D* histo_WZ_JESUp            = new TH1D( Form("histo_WZ_%sUp","CMS_scale_j")  , Form("histo_WZ_%sUp","CMS_scale_j")  , nBin, xbins); histo_WZ_JESUp  ->Sumw2();
  TH1D* histo_WZ_JESDown          = new TH1D( Form("histo_WZ_%sDown","CMS_scale_j"), Form("histo_WZ_%sDown","CMS_scale_j"), nBin, xbins); histo_WZ_JESDown->Sumw2();
  TH1D* histo_ZZ_JESUp            = new TH1D( Form("histo_ZZ_%sUp","CMS_scale_j")  , Form("histo_ZZ_%sUp","CMS_scale_j")  , nBin, xbins); histo_ZZ_JESUp  ->Sumw2();
  TH1D* histo_ZZ_JESDown          = new TH1D( Form("histo_ZZ_%sDown","CMS_scale_j"), Form("histo_ZZ_%sDown","CMS_scale_j"), nBin, xbins); histo_ZZ_JESDown->Sumw2();
  TH1D* histo_VVV_JESUp           = new TH1D( Form("histo_VVV_%sUp","CMS_scale_j")  , Form("histo_VVV_%sUp","CMS_scale_j")  , nBin, xbins); histo_VVV_JESUp  ->Sumw2();
  TH1D* histo_VVV_JESDown         = new TH1D( Form("histo_VVV_%sDown","CMS_scale_j"), Form("histo_VVV_%sDown","CMS_scale_j"), nBin, xbins); histo_VVV_JESDown->Sumw2();
  TH1D* histo_Higgs_JESUp         = new TH1D( Form("histo_Higgs_%sUp","CMS_scale_j")  , Form("histo_Higgs_%sUp","CMS_scale_j")  , nBin, xbins); histo_Higgs_JESUp  ->Sumw2();
  TH1D* histo_Higgs_JESDown       = new TH1D( Form("histo_Higgs_%sDown","CMS_scale_j"), Form("histo_Higgs_%sDown","CMS_scale_j"), nBin, xbins); histo_Higgs_JESDown->Sumw2();

  TH1D* histo_WZ_JERUp            = new TH1D( Form("histo_WZ_%sUp","CMS_res_j")  , Form("histo_WZ_%sUp","CMS_res_j")  , nBin, xbins); histo_WZ_JERUp  ->Sumw2();
  TH1D* histo_WZ_JERDown          = new TH1D( Form("histo_WZ_%sDown","CMS_res_j"), Form("histo_WZ_%sDown","CMS_res_j"), nBin, xbins); histo_WZ_JERDown->Sumw2();
  TH1D* histo_ZZ_JERUp            = new TH1D( Form("histo_ZZ_%sUp","CMS_res_j")  , Form("histo_ZZ_%sUp","CMS_res_j")  , nBin, xbins); histo_ZZ_JERUp  ->Sumw2();
  TH1D* histo_ZZ_JERDown          = new TH1D( Form("histo_ZZ_%sDown","CMS_res_j"), Form("histo_ZZ_%sDown","CMS_res_j"), nBin, xbins); histo_ZZ_JERDown->Sumw2();
  TH1D* histo_VVV_JERUp           = new TH1D( Form("histo_VVV_%sUp","CMS_res_j")  , Form("histo_VVV_%sUp","CMS_res_j")  , nBin, xbins); histo_VVV_JERUp  ->Sumw2();
  TH1D* histo_VVV_JERDown         = new TH1D( Form("histo_VVV_%sDown","CMS_res_j"), Form("histo_VVV_%sDown","CMS_res_j"), nBin, xbins); histo_VVV_JERDown->Sumw2();
  TH1D* histo_Higgs_JERUp         = new TH1D( Form("histo_Higgs_%sUp","CMS_res_j")  , Form("histo_Higgs_%sUp","CMS_res_j")  , nBin, xbins); histo_Higgs_JERUp  ->Sumw2();
  TH1D* histo_Higgs_JERDown       = new TH1D( Form("histo_Higgs_%sDown","CMS_res_j"), Form("histo_Higgs_%sDown","CMS_res_j"), nBin, xbins); histo_Higgs_JERDown->Sumw2();

  TH1D* histo_Wjets_WUp           = new TH1D( Form("histo_Wjets_CMS_qqwz_MVAWUp"),   Form("histo_Wjets_CMS_qqwz_MVAWUp"),   nBin, xbins); histo_Wjets_WUp  ->Sumw2();
  TH1D* histo_Wjets_WDown         = new TH1D( Form("histo_Wjets_CMS_qqwz_MVAWDown"), Form("histo_Wjets_CMS_qqwz_MVAWDown"), nBin, xbins); histo_Wjets_WDown->Sumw2();

  double nSelectedData[nSelTypes*2];
  double bgdDecay[nSelTypes*2][50],weiDecay[nSelTypes*2][50];
  double bgdDecaySyst[nSelTypesSyst*2][50],weiDecaySyst[nSelTypesSyst*2][50];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    nSelectedData[i] = 0.0; 
    for(int j=0; j<50; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    for(int j=0; j<50; j++) {
      bgdDecaySyst[i][j] = 0.0; weiDecaySyst[i][j] = 0.0; 
    }
  }

  unsigned int patternTopVeto = SmurfTree::TopVeto;
  float ewkMVA = -999.;
  bgdEvent.tree_->SetBranchAddress("ewkMVA", &ewkMVA );

  int nBgd=bgdEvent.tree_->GetEntries();

  for (int evt=0; evt<nBgd; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nBgd);
    bgdEvent.tree_->GetEntry(evt);

    // generator level selection
    bool minGenCuts = !(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
                        ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
			((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection && bgdEvent.lid3_ != 0));
    bool genLevelSel[2] = {false, false};
    if(minGenCuts == true) {
      genLevelNorm[0]++;
      if(bgdEvent.genlep1_.Pt() > 0 
      && bgdEvent.genlep2_.Pt() > 0
      //&& bgdEvent.genlep3_.Pt() > 0
        ) {
        genLevelNorm[1]++;
        bool passJetsCuts[3]; 
        passJetsCuts[0] = bgdEvent.genjet1_.Pt() > 20 && TMath::Abs(bgdEvent.genjet1_.Eta()) < 5.0 && bgdEvent.genjet2_.Pt() > 20 && TMath::Abs(bgdEvent.genjet2_.Eta()) < 5.0 && TMath::Abs(bgdEvent.genjet1_.Eta()-bgdEvent.genjet2_.Eta()) > 2.5 && (bgdEvent.genjet1_+bgdEvent.genjet2_).M() > 300.;
        passJetsCuts[1] = bgdEvent.genjet1_.Pt() > 20 && TMath::Abs(bgdEvent.genjet1_.Eta()) < 5.0 && bgdEvent.genjet3_.Pt() > 20 && TMath::Abs(bgdEvent.genjet3_.Eta()) < 5.0 && TMath::Abs(bgdEvent.genjet1_.Eta()-bgdEvent.genjet3_.Eta()) > 2.5 && (bgdEvent.genjet1_+bgdEvent.genjet3_).M() > 300.;
        passJetsCuts[2] = bgdEvent.genjet2_.Pt() > 20 && TMath::Abs(bgdEvent.genjet2_.Eta()) < 5.0 && bgdEvent.genjet3_.Pt() > 20 && TMath::Abs(bgdEvent.genjet3_.Eta()) < 5.0 && TMath::Abs(bgdEvent.genjet2_.Eta()-bgdEvent.genjet3_.Eta()) > 2.5 && (bgdEvent.genjet2_+bgdEvent.genjet3_).M() > 300.;
        if(bgdEvent.genlep1_.Pt() > 10 && TMath::Abs(bgdEvent.genlep1_.Eta()) < 2.5 && 
           bgdEvent.genlep2_.Pt() > 10 && TMath::Abs(bgdEvent.genlep2_.Eta()) < 2.5 &&
           bgdEvent.genlep3_.Pt() > 10 && TMath::Abs(bgdEvent.genlep3_.Eta()) < 2.5 &&
	   (passJetsCuts[0] || passJetsCuts[1] || passJetsCuts[2])) {
          genLevelNorm[2]++;
	  genLevelSel[0] = true;
          passJetsCuts[0] = bgdEvent.genjet1_.Pt() > 30 && TMath::Abs(bgdEvent.genjet1_.Eta()) < 4.7 && bgdEvent.genjet2_.Pt() > 30 && TMath::Abs(bgdEvent.genjet2_.Eta()) < 4.7 && TMath::Abs(bgdEvent.genjet1_.Eta()-bgdEvent.genjet2_.Eta()) > 2.5 && (bgdEvent.genjet1_+bgdEvent.genjet2_).M() > 500.;
          passJetsCuts[1] = bgdEvent.genjet1_.Pt() > 30 && TMath::Abs(bgdEvent.genjet1_.Eta()) < 4.7 && bgdEvent.genjet3_.Pt() > 30 && TMath::Abs(bgdEvent.genjet3_.Eta()) < 4.7 && TMath::Abs(bgdEvent.genjet1_.Eta()-bgdEvent.genjet3_.Eta()) > 2.5 && (bgdEvent.genjet1_+bgdEvent.genjet3_).M() > 500.;
          passJetsCuts[2] = bgdEvent.genjet2_.Pt() > 30 && TMath::Abs(bgdEvent.genjet2_.Eta()) < 4.7 && bgdEvent.genjet3_.Pt() > 30 && TMath::Abs(bgdEvent.genjet3_.Eta()) < 4.7 && TMath::Abs(bgdEvent.genjet2_.Eta()-bgdEvent.genjet3_.Eta()) > 2.5 && (bgdEvent.genjet2_+bgdEvent.genjet3_).M() > 500.;
          if(bgdEvent.genlep1_.Pt() > 20 && TMath::Abs(bgdEvent.genlep1_.Eta()) < 2.5 && 
             bgdEvent.genlep2_.Pt() > 20 && TMath::Abs(bgdEvent.genlep2_.Eta()) < 2.5 &&
             bgdEvent.genlep3_.Pt() > 20 && TMath::Abs(bgdEvent.genlep3_.Eta()) < 2.5 &&
	     (passJetsCuts[0] || passJetsCuts[1] || passJetsCuts[2])) {
            genLevelNorm[3]++;
	    genLevelSel[1] = true;
	  }
        }
      }
    }

    if(bgdEvent.lep1_.Pt() < 1.0) continue;

    if(!(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
         ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
	  bgdEvent.dstype_ != SmurfTree::data)) continue;
    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::wjets 	   ) fDecay =  3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay =  5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  	   ) fDecay =  9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  	   ) fDecay =  9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqww2j  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wwewk  	   ) fDecay = 31;
    else if(bgdEvent.dstype_ == SmurfTree::wzewk  	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(bgdEvent.processId_==121 ||
            bgdEvent.processId_==122)   fDecay = 41;
    else if(bgdEvent.processId_==24)    fDecay = 42;
    else if(bgdEvent.processId_==26)    fDecay = 43;
    else if(bgdEvent.processId_==10001) fDecay = 44;
    else if(bgdEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;}

    if (!use_fake_rate_method && fDecay == 1)
      continue;

    int correctQ = 0;
    if(bgdEvent.lep1McId_ ==  11 && bgdEvent.lq1_ < 0) correctQ++;
    if(bgdEvent.lep1McId_ == -11 && bgdEvent.lq1_ > 0) correctQ++;
    if(bgdEvent.lep1McId_ ==  13 && bgdEvent.lq1_ < 0) correctQ++;
    if(bgdEvent.lep1McId_ == -13 && bgdEvent.lq1_ > 0) correctQ++;
    if(bgdEvent.lep2McId_ ==  11 && bgdEvent.lq2_ < 0) correctQ++;
    if(bgdEvent.lep2McId_ == -11 && bgdEvent.lq2_ > 0) correctQ++;
    if(bgdEvent.lep2McId_ ==  13 && bgdEvent.lq2_ < 0) correctQ++;
    if(bgdEvent.lep2McId_ == -13 && bgdEvent.lq2_ > 0) correctQ++;

    if(fDecay == 29 && correctQ != 2) fDecay = 30;
    if(fDecay == 27 && correctQ  < 2) fDecay = 30;

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false, false, false},
			                     {false, false, false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false, false},
                                   {false, false}};
    bool isRealLepton = false;
    if(bgdEvent.lid3_ == 0) {
      if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
         (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;
    }
    else {
      if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
         (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13) &&
         (TMath::Abs(bgdEvent.lep3McId_) == 11 || TMath::Abs(bgdEvent.lep3McId_) == 13)) isRealLepton = true;
    }

    double theMET = bgdEvent.met_; double theMETPHI = bgdEvent.metPhi_; 
    
    double massZMin = trilepton_info(0,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                       bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                       bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				       bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);

    int lType = 1;
    //if     (bgdEvent.lq1_ * bgdEvent.lq2_ < 0) lType = 0;

    double zeppenfeld = TMath::Min(TMath::Max(TMath::Abs(bgdEvent.lep1_.Eta()-(bgdEvent.jet1_.Eta()+bgdEvent.jet2_.Eta())/2.),
                                              TMath::Abs(bgdEvent.lep2_.Eta()-(bgdEvent.jet1_.Eta()+bgdEvent.jet2_.Eta())/2.)),3.999);
    //int centrality = 0;
    //if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta()+0.1 > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
    //    (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta()+0.1 > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
    //   ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta()+0.1 > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
    //    (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta()+0.1 > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 40.0; if(bgdEvent.type_ != SmurfTree::ee) metMin = 40.0;

    // jet1McId and jet2McId actually provide lepton rejection information
    // hasZCand = reject events with |mll-mZ|<15
    // trackSel[0] == reject events with |mll-mZ|<10 for pf candidates with |eta|<3.0
    // trackSel[1] == reject events with |mll-mZ|<10 for pf candidates with |eta|<4.7
    // trackSel[2] == reject events with isolated reconstructed leptons with pt>10 and iso/pt<0.1
    // trackSel[3] == reject events with isolated tracks with pt>10 and iso/pt<0.1 (not used by default)
    int newId=int(bgdEvent.jet1McId_);
    //int tauId=int((bgdEvent.jet1McId_%100-bgdEvent.jet1McId_%10)/10);
    int qDisAgree=int((newId%1000-newId%100)/100);
    //int hasZCand=int(newId/1000);
    //int trackSel[4] = {int((bgdEvent.jet2McId_%100-bgdEvent.jet2McId_%10)/10),int((bgdEvent.jet2McId_%1000-bgdEvent.jet2McId_%100)/100),int((bgdEvent.jet2McId_%10000-bgdEvent.jet2McId_%1000)/1000),int(bgdEvent.jet2McId_/10000)};

    bool passNjets    = bgdEvent.njets_ >= 2;
    bool passMET      = bgdEvent.met_ > metMin;
    bool preselCuts   = bgdEvent.lep1_.Pt() > 20. && bgdEvent.lep2_.Pt() > 20.;
    bool passBtagVeto = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool passVBFSel   = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 2.5;
    // quick way to accept W+W+ (1) or W-W- (2)
    bool passNsignSel = true;
    if(TMath::Abs(bgdEvent.lq1_+bgdEvent.lq2_+bgdEvent.lq3_) != 1) passNsignSel = false;

    // JER
    LorentzVector genjet_jer1;
    LorentzVector genjet_jer2;
    int jerbin[2] = {-1, -1};
    if     (DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet1_.Phi(),bgdEvent.genjet1_.Eta()) < 0.5) {
      genjet_jer1 = bgdEvent.genjet1_;
    }
    else if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet2_.Phi(),bgdEvent.genjet2_.Eta()) < 0.5) {
      genjet_jer1 = bgdEvent.genjet2_;
    }
    else if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet3_.Phi(),bgdEvent.genjet3_.Eta()) < 0.5) {
      genjet_jer1 = bgdEvent.genjet3_;
    }
    else {
      genjet_jer1 = bgdEvent.jet1_;
    }
    if     (DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.genjet1_.Phi(),bgdEvent.genjet1_.Eta()) < 0.5) {
      genjet_jer2 = bgdEvent.genjet1_;
    }
    else if(DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.genjet2_.Phi(),bgdEvent.genjet2_.Eta()) < 0.5) {
      genjet_jer2 = bgdEvent.genjet2_;
    }
    else if(DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.genjet3_.Phi(),bgdEvent.genjet3_.Eta()) < 0.5) {
      genjet_jer2 = bgdEvent.genjet3_;
    }
    else {
      genjet_jer2 = bgdEvent.jet2_;
    }
    if     (TMath::Abs(genjet_jer1.Eta()) < 0.5) jerbin[0] = 0;
    else if(TMath::Abs(genjet_jer1.Eta()) < 1.1) jerbin[0] = 1;
    else if(TMath::Abs(genjet_jer1.Eta()) < 1.7) jerbin[0] = 2;
    else if(TMath::Abs(genjet_jer1.Eta()) < 2.3) jerbin[0] = 3;
    else                                         jerbin[0] = 4;
    if     (TMath::Abs(genjet_jer2.Eta()) < 0.5) jerbin[1] = 0;
    else if(TMath::Abs(genjet_jer2.Eta()) < 1.1) jerbin[1] = 1;
    else if(TMath::Abs(genjet_jer2.Eta()) < 1.7) jerbin[1] = 2;
    else if(TMath::Abs(genjet_jer2.Eta()) < 2.3) jerbin[1] = 3;
    else                                         jerbin[1] = 4;

    LorentzVector genjetU_jer1 = genjet_jer1; 
    LorentzVector genjetD_jer1 = genjet_jer1; 
    LorentzVector genjetU_jer2 = genjet_jer2;
    LorentzVector genjetD_jer2 = genjet_jer2;
    double factJER[3]; // order is crucial
    if(genjet_jer1.Pt()>0){
      factJER[0] = TMath::Max(0.0,genjet_jer1.Pt()+JERVal [jerbin[0]]*(bgdEvent.jet1_.Pt()-genjet_jer1.Pt()))/genjet_jer1.Pt();
      factJER[1] = TMath::Max(0.0,genjet_jer1.Pt()+JERUp  [jerbin[0]]*(bgdEvent.jet1_.Pt()-genjet_jer1.Pt()))/genjet_jer1.Pt();
      factJER[2] = TMath::Max(0.0,genjet_jer1.Pt()+JERDown[jerbin[0]]*(bgdEvent.jet1_.Pt()-genjet_jer1.Pt()))/genjet_jer1.Pt();
      genjetU_jer1.SetPx(genjet_jer1.Px()*factJER[1]);
      genjetU_jer1.SetPy(genjet_jer1.Py()*factJER[1]);
      genjetU_jer1.SetPz(genjet_jer1.Pz()*factJER[1]);
      genjetU_jer1.SetE (genjet_jer1.E( )*factJER[1]);
      genjetD_jer1.SetPx(genjet_jer1.Px()*factJER[2]);
      genjetD_jer1.SetPy(genjet_jer1.Py()*factJER[2]);
      genjetD_jer1.SetPz(genjet_jer1.Pz()*factJER[2]);
      genjetD_jer1.SetE (genjet_jer1.E( )*factJER[2]);
      genjet_jer1.SetPx (genjet_jer1.Px()*factJER[0]);
      genjet_jer1.SetPy (genjet_jer1.Py()*factJER[0]);
      genjet_jer1.SetPz (genjet_jer1.Pz()*factJER[0]);
      genjet_jer1.SetE  (genjet_jer1.E( )*factJER[0]);
    }
    if(genjet_jer2.Pt()>0){
      factJER[0] = TMath::Max(0.0,genjet_jer2.Pt()+JERVal [jerbin[1]]*(bgdEvent.jet2_.Pt()-genjet_jer2.Pt()))/genjet_jer2.Pt();
      factJER[1] = TMath::Max(0.0,genjet_jer2.Pt()+JERUp  [jerbin[1]]*(bgdEvent.jet2_.Pt()-genjet_jer2.Pt()))/genjet_jer2.Pt();
      factJER[2] = TMath::Max(0.0,genjet_jer2.Pt()+JERDown[jerbin[1]]*(bgdEvent.jet2_.Pt()-genjet_jer2.Pt()))/genjet_jer2.Pt();
      genjetU_jer2.SetPx(genjet_jer2.Px()*factJER[1]);
      genjetU_jer2.SetPy(genjet_jer2.Py()*factJER[1]);
      genjetU_jer2.SetPz(genjet_jer2.Pz()*factJER[1]);
      genjetU_jer2.SetE (genjet_jer2.E( )*factJER[1]);
      genjetD_jer2.SetPx(genjet_jer2.Px()*factJER[2]);
      genjetD_jer2.SetPy(genjet_jer2.Py()*factJER[2]);
      genjetD_jer2.SetPz(genjet_jer2.Pz()*factJER[2]);
      genjetD_jer2.SetE (genjet_jer2.E( )*factJER[2]);
      genjet_jer2.SetPx (genjet_jer2.Px()*factJER[0]);
      genjet_jer2.SetPy (genjet_jer2.Py()*factJER[0]);
      genjet_jer2.SetPz (genjet_jer2.Pz()*factJER[0]);
      genjet_jer2.SetE  (genjet_jer2.E( )*factJER[0]);
    }

    double deltaRlJMin = 999.0;
    if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta());
    if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta());
    if(DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta());
    if(DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta());

    // 0      1      2       3     4   5      6        7           8  9            10            11     12  13    14  15
    // lep1pt,lep2pt,dilmass,dilpt,met,metPhi,trackMet,trackMetPhi,mt,dPhiDiLepMET,dPhiMETTrkMET,pTFrac,mtH,mlljj,mjj,lep3pt;
    double outputVarLepP[16];
    makeSystematicEffects3l(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lid3_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.lep3_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 0, outputVarLepP);
    double outputVarLepM[16];
    makeSystematicEffects3l(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lid3_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.lep3_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_,
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 1, outputVarLepM);
    double outputVarMET[16];
    makeSystematicEffects3l(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lid3_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.lep3_, bgdEvent.dilep_,
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 2, outputVarMET);
    double outputVar[16];
    makeSystematicEffects3l(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lid3_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.lep3_, bgdEvent.dilep_,
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 3, outputVar);
    double outputVarJESP[16];
    makeSystematicEffects3l(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lid3_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.lep3_, bgdEvent.dilep_,
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 4, outputVarJESP);
    double outputVarJESM[16];
    makeSystematicEffects3l(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lid3_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.lep3_, bgdEvent.dilep_,
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 5, outputVarJESM);
    double MVAVar[8] = {outputVar[13],outputVarJESP[13],outputVarJESM[13],outputVarLepP[13],outputVarLepM[13],outputVarMET[13],(bgdEvent.lep1_+bgdEvent.lep2_+genjetU_jer1+genjetU_jer2).M(),(bgdEvent.lep1_+bgdEvent.lep2_+genjetD_jer1+genjetD_jer2).M()};
    if     (thePlot == 0) {MVAVar[0]=outputVar[14];MVAVar[1]=outputVarJESP[14];MVAVar[2]=outputVarJESM[14];MVAVar[3]=outputVarLepP[14];MVAVar[4]=outputVarLepM[14];MVAVar[5]=outputVarMET[14];MVAVar[6]=(genjetU_jer1+genjetU_jer2).M();MVAVar[7]=(genjetD_jer1+genjetD_jer2).M();}
    else if(thePlot == 2) {MVAVar[0]=outputVar[ 0];MVAVar[1]=outputVarJESP[ 0];MVAVar[2]=outputVarJESM[ 0];MVAVar[3]=outputVarLepP[ 0];MVAVar[4]=outputVarLepM[ 0];MVAVar[5]=outputVarMET[ 0];MVAVar[6]=outputVar[ 0];MVAVar[7]=outputVar[ 0];}
    else if(thePlot == 9) {MVAVar[0]=outputVar[ 2];MVAVar[1]=outputVarJESP[ 2];MVAVar[2]=outputVarJESM[ 2];MVAVar[3]=outputVarLepP[ 2];MVAVar[4]=outputVarLepM[ 2];MVAVar[5]=outputVarMET[ 2];MVAVar[6]=outputVar[ 2];MVAVar[7]=outputVar[ 2];}
    else if(thePlot ==19) {MVAVar[0]=outputVar[12];MVAVar[1]=outputVarJESP[12];MVAVar[2]=outputVarJESM[12];MVAVar[3]=outputVarLepP[12];MVAVar[4]=outputVarLepM[12];MVAVar[5]=outputVarMET[12];MVAVar[6]=outputVar[12];MVAVar[7]=outputVar[12];}
    for(int nv=0; nv<8; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBin]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff  = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 0)*
    		 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 1);
      addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1);
    } else {addLepEff = 1.0;}

    double NjetSyst[4] = {0., 0., 0., 0.};
    if(bgdEvent.jet1_.Pt()*1.04 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet2_.Pt()*1.04 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet3_.Pt()*1.04 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet4_.Pt()*1.04 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet1_.Pt()*0.96 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet2_.Pt()*0.96 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet3_.Pt()*0.96 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet4_.Pt()*0.96 > ptJetMin) NjetSyst[1]++;
    if(genjetU_jer1.Pt()        > ptJetMin) NjetSyst[2]++;
    if(genjetU_jer2.Pt()        > ptJetMin) NjetSyst[2]++;
    if(bgdEvent.jet3_.Pt()      > ptJetMin) NjetSyst[2]++;
    if(bgdEvent.jet4_.Pt()      > ptJetMin) NjetSyst[2]++;
    if(genjetD_jer1.Pt()        > ptJetMin) NjetSyst[3]++;
    if(genjetD_jer2.Pt()        > ptJetMin) NjetSyst[3]++;
    if(bgdEvent.jet3_.Pt()      > ptJetMin) NjetSyst[3]++;
    if(bgdEvent.jet4_.Pt()      > ptJetMin) NjetSyst[3]++;

    bool passLSel = false;
    if     (lSel == 0 && bgdEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && bgdEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && bgdEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && bgdEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (bgdEvent.type_ == SmurfTree::me || bgdEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passNsignSel && qDisAgree == 0 && NjetSyst[0] >= 2	        && outputVarJESP[4] > metMin && outputVar[0]	 > 20.0 && outputVar[1]     > 20.0 && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && bgdEvent.lep3_.Pt() > 10. && outputVarJESP[14] > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 2.5) passSystCuts[lType][JESUP] = true;
    if(passNsignSel && qDisAgree == 0 && NjetSyst[1] >= 2	        && outputVarJESM[4] > metMin && outputVar[0]	 > 20.0 && outputVar[1]     > 20.0 && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && bgdEvent.lep3_.Pt() > 10. && outputVarJESM[14] > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 2.5) passSystCuts[lType][JESDOWN] = true;
    if(passNsignSel && qDisAgree == 0 && bgdEvent.jet2_.Pt() > ptJetMin && outputVarLepP[4] > metMin && outputVarLepP[0] > 20.0 && outputVarLepP[1] > 20.0 && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && outputVarLepP[15]   > 10. && outputVarLepP[14] > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 2.5) passSystCuts[lType][LEPP] = true;
    if(passNsignSel && qDisAgree == 0 && bgdEvent.jet2_.Pt() > ptJetMin && outputVarLepM[4] > metMin && outputVarLepM[0] > 20.0 && outputVarLepM[1] > 20.0 && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && outputVarLepM[15]   > 10. && outputVarLepM[14] > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 2.5) passSystCuts[lType][LEPM] = true;
    if(passNsignSel && qDisAgree == 0 && bgdEvent.jet2_.Pt() > ptJetMin && outputVarMET[4]  > metMin && outputVarMET[0]  > 20.0 && outputVarMET[1]  > 20.0 && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && bgdEvent.lep3_.Pt() > 10. && outputVarMET[14]  > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 2.5) passSystCuts[lType][MET] = true;

    if(passNsignSel && qDisAgree == 0 && NjetSyst[2] >= 2 && passMET && preselCuts && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && bgdEvent.lep3_.Pt() > 10. && (genjetU_jer1+genjetU_jer2).M() > 500 && TMath::Abs(genjetU_jer1.Eta()-genjetU_jer2.Eta()) > 2.5) passSystCuts[lType][JERUP]   = true;
    if(passNsignSel && qDisAgree == 0 && NjetSyst[3] >= 2 && passMET && preselCuts && passBtagVeto && bgdEvent.lid3_ != 0 && massZMin < 15.0 && bgdEvent.lep3_.Pt() > 10. && (genjetD_jer1+genjetD_jer2).M() > 500 && TMath::Abs(genjetD_jer1.Eta()-genjetD_jer2.Eta()) > 2.5) passSystCuts[lType][JERDOWN] = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && bgdEvent.dilep_.M() > 15.0 && qDisAgree == 0) {
       
       if(passNsignSel && !passBtagVeto && passVBFSel == true && bgdEvent.lid3_ != 0 && 		   bgdEvent.lep3_.Pt() > 10.) passCuts[lType][BTAGSEL] = true;
       if(passNsignSel &&  passBtagVeto && passVBFSel == true && bgdEvent.lid3_ != 0 && massZMin < 15.0 && bgdEvent.lep3_.Pt() > 10.) passCuts[lType][WZSEL] = true;

      if(use_fake_rate_method && isRealLepton == false &&
         (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
          bgdEvent.dstype_ == SmurfTree::qqww	|| bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
          bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) 
        {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false;
	 for(unsigned int i=0; i<nSelTypesSyst-2; i++) passSystCuts[lType][i] = false;
	}
    }

    //otherwise the ttbar fake background from mc would get counted as wrong sign
    if (!use_fake_rate_method && isRealLepton == false && fDecay == 5)
      fDecay = 45;
    

    if(1){
      double theWeight = 0.0;
      double add       = 1.0;
      int nFake = 0;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(nFake < 0) assert(0);
 
      if(nFake > 1){

	if (!use_fake_rate_method)
	  continue;

	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	fDecay = 22;
	theWeight	       = -1.0*add*frCorr;
      }
      else if(nFake == 1){

	if (!use_fake_rate_method)
	  continue;

        if(bgdEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
          if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  // new category, W+jetsM
	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0*frCorr;
	}
	else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);

	  add = add*nPUScaleFactor2012(fhDPU ,bgdEvent.npu_);

          add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add = add*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt(), 
								   fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						   TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));

          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
            double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	 							      fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        						             TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      	    double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	  							      fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
      	  							     TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      	    double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     	  							      fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
      	  							     TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      	    trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
         }
	  
	  add = add*trigEff;
	  if(fCheckProblem == true && fDecay != 44 && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMB: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  fDecay                 = 1;

	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * bgdEvent.scale1fb_*lumi*add*frCorr;
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	        					        fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
							        TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        double sf_eff = 1.0;
	sf_eff = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_)*
        	 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);

        theWeight = ZttScaleFactor(period,bgdEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(bgdEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,bgdEvent.npu_);
        double add2 = 1.0;
	add2 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	add2 = add2*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add2 = add2*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
								 fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						 TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));

        if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
           double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     								     fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
                						    TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      	   double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     								     fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
      								    TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      	   double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     								     fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
      								    TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      	   trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
        }
        add = add1*add2*trigEff;

	 if(fCheckProblem == true && fDecay != 44 && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) && add != 0 && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	 printf("PROBLEMCB(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",bgdEvent.event_,add1,add2,trigEff,add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);


	if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,theMET);
        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (bgdEvent.dstype_ == SmurfTree::dymm || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dytt) &&
          (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) add = 0.0;

	theWeight              = bgdEvent.scale1fb_*lumi*add;
      }

      // uncertainty related to wrong-sign leptons
      double weightWS[2] = {theWeight,theWeight};
      scaleFactor_WS(bgdEvent.lep1_,bgdEvent.lq1_,bgdEvent.lid1_,bgdEvent.lep1McId_,weightWS,0);
      scaleFactor_WS(bgdEvent.lep2_,bgdEvent.lq2_,bgdEvent.lid2_,bgdEvent.lep2McId_,weightWS,0);
      if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
        scaleFactor_WS(bgdEvent.lep3_,bgdEvent.lq3_,bgdEvent.lid3_,bgdEvent.lep3McId_,weightWS,0);
      }
      theWeight = weightWS[0];

      if(passBtagVeto) {
        if(bgdEvent.dstype_ == SmurfTree::ttbar) theWeight = theWeight * 1.08767;
        if(bgdEvent.dstype_ == SmurfTree::tw)    theWeight = theWeight * 1.08767;
      } else {
        if(bgdEvent.dstype_ == SmurfTree::ttbar) theWeight = theWeight / 1.08767;
        if(bgdEvent.dstype_ == SmurfTree::tw)    theWeight = theWeight / 1.08767;
      }

      if(minGenCuts == true && (passCuts[0][WZSEL]||passCuts[1][WZSEL]) && genLevelSel[0] == false) genLevelNorm[4]++;
      if(minGenCuts == true && (passCuts[0][WZSEL]||passCuts[1][WZSEL]) && genLevelSel[0] == true)  genLevelNorm[5]++;
      if(minGenCuts == true && (passCuts[0][WZSEL]||passCuts[1][WZSEL]) && genLevelSel[1] == false) genLevelNorm[6]++;
      if(minGenCuts == true && (passCuts[0][WZSEL]||passCuts[1][WZSEL]) && genLevelSel[1] == true)  genLevelNorm[7]++;

      double MT3lTotx = bgdEvent.met_*cos(bgdEvent.metPhi_)-bgdEvent.lep1_.Px()-bgdEvent.lep2_.Px()-bgdEvent.lep3_.Px();
      double MT3lToty = bgdEvent.met_*sin(bgdEvent.metPhi_)-bgdEvent.lep1_.Py()-bgdEvent.lep2_.Py()-bgdEvent.lep3_.Py();
      double MT3lTot = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);

      if(passCuts[1][WZSEL]){ // begin making plots
	double myVar = -1.0;
	if     (thePlot == 0) myVar = TMath::Max(TMath::Min((bgdEvent.jet1_+bgdEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot == 1) myVar = TMath::Max(TMath::Min((bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot == 2) myVar = TMath::Min(bgdEvent.lep1_.Pt(),499.999);
	else if(thePlot == 3) myVar = TMath::Min(bgdEvent.lep2_.Pt(),399.999);
	else if(thePlot == 4) myVar = bgdEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = bgdEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = bgdEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = bgdEvent.mt_;
	else if(thePlot == 8) myVar = bgdEvent.dilep_.Pt();
	else if(thePlot == 9) myVar = TMath::Min(bgdEvent.dilep_.M(),499.999);
	else if(thePlot ==10) myVar = bgdEvent.njets_;
	else if(thePlot ==11) myVar = bgdEvent.nvtx_;
	else if(thePlot ==12) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==13) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==14) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==15) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
	else if(thePlot ==16) myVar = bgdEvent.type_;
	else if(thePlot ==17) myVar = bgdEvent.dR_;
	else if(thePlot ==18) myVar = zeppenfeld;
	else if(thePlot ==19) myVar = TMath::Min(MT3lTot,999.999);
	else if(thePlot ==20) myVar = TMath::Min(massZMin,99.999);
	else if(thePlot ==21) myVar = TMath::Min((double)bgdEvent.met_,199.999);
	else if(thePlot ==22) myVar = TMath::Min(deltaRlJMin,4.999);
	else assert(0);

      	if(fDecay == 27){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 1 || fDecay == 23 || fDecay == 3 || fDecay == 45){
      	  histo3->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 30 || fDecay == 28 || fDecay == 29 || fDecay == 31 ||
                fDecay ==  5 || fDecay == 13 || fDecay == 20 || 
		fDecay == 10 || fDecay ==  9 || fDecay == 19){
      	  histo4->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 21){
      	  histo5->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 41 || fDecay == 42 || fDecay == 43 || fDecay == 44){
      	  histo7->Fill(myVar,theWeight);
      	}
      	else {
      	  printf("NOOOOOOOOOOOOOOOOOOOO %d\n",fDecay);
      	}
      } // end making plots
      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            bgdDecay[i+j*nSelTypes][(int)fDecay] += theWeight;
            weiDecay[i+j*nSelTypes][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      for(unsigned int i=0; i<7; i++) {
        for(int j=0; j<2; j++){
          if(passSystCuts[j][i]) {
            bgdDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight;
            weiDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][WZSEL]) {
        bgdDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffUp  /addLepEff;
        weiDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        bgdDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffDown/addLepEff;
        weiDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }

      if     (fDecay == 21){
        if(passCuts[1][WZSEL])  	     histo_VVV           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WZSEL])  	     histo_VVV_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WZSEL])  	     histo_VVV_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_VVV_JESUp	 ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_VVV_JESDown	 ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_VVV_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_VVV_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_VVV_METResUp  ->Fill(MVAVar[5], theWeight);;
        if(passSystCuts[1][JERUP  ] == true) histo_VVV_JERUp	 ->Fill(MVAVar[6], theWeight);
        if(passSystCuts[1][JERDOWN] == true) histo_VVV_JERDown	 ->Fill(MVAVar[7], theWeight);
      }
      else if(fDecay == 44){
        if(passCuts[1][WZSEL])  	     histo_Higgs           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WZSEL])  	     histo_Higgs_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WZSEL])  	     histo_Higgs_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Higgs_JESUp	   ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Higgs_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Higgs_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Higgs_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Higgs_METResUp  ->Fill(MVAVar[5], theWeight);;
        if(passSystCuts[1][JERUP  ] == true) histo_Higgs_JERUp	   ->Fill(MVAVar[6], theWeight);
        if(passSystCuts[1][JERDOWN] == true) histo_Higgs_JERDown   ->Fill(MVAVar[7], theWeight);
      }
      else if(fDecay == 27){
        if(passCuts[1][WZSEL])  	     histo_WZ           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WZSEL])  	     histo_WZ_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WZSEL])  	     histo_WZ_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WZ_JESUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WZ_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WZ_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WZ_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WZ_METResUp  ->Fill(MVAVar[5], theWeight);;
        if(passSystCuts[1][JERUP  ] == true) histo_WZ_JERUp     ->Fill(MVAVar[6], theWeight);
        if(passSystCuts[1][JERDOWN] == true) histo_WZ_JERDown   ->Fill(MVAVar[7], theWeight);
      }
      else if(fDecay == 30 || fDecay == 28 || fDecay == 29 || fDecay == 31 ||
              fDecay ==  5 || fDecay == 13 || fDecay == 20 || 
	      fDecay == 10 || fDecay ==  9 || fDecay == 19){
        if(passCuts[1][WZSEL])  	     histo_ZZ           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WZSEL])  	     histo_ZZ_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WZSEL])  	     histo_ZZ_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_ZZ_JESUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_ZZ_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_ZZ_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_ZZ_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_ZZ_METResUp  ->Fill(MVAVar[5], theWeight);;
        if(passSystCuts[1][JERUP  ] == true) histo_ZZ_JERUp     ->Fill(MVAVar[6], theWeight);
        if(passSystCuts[1][JERDOWN] == true) histo_ZZ_JERDown   ->Fill(MVAVar[7], theWeight);
      }
      else if(fDecay == 1 || fDecay == 23 || fDecay == 3 || fDecay == 45 ){
	if (use_fake_rate_method){
	  double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
					  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
	  addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
				   (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
					  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
	  addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
				   (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  if(passCuts[1][WZSEL]) 	     histo_Wjets		       ->Fill(MVAVar[0], theWeight);
	  if(passCuts[1][WZSEL]) 	     histo_Wjets_WUp    ->Fill(MVAVar[0], theWeight*addFRS/addFR);
	}
	else {
	  if(passCuts[1][WZSEL]) 	     histo_Wjets		       ->Fill(MVAVar[0], theWeight);
	  if(passCuts[1][WZSEL]) 	     histo_Wjets_WUp    ->Fill(MVAVar[0], theWeight);
	}

      }
      else {
        printf("Forbidden fDecay: %d\n",fDecay);
	assert(0);
      }
    } // if passCuts
  } // end background loop
  
  if(run_over_data){
  dataEvent.tree_->SetBranchAddress("ewkMVA", &ewkMVA );
  int nData=dataEvent.tree_->GetEntries();
  for (int evt=0; evt<nData; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nData);
    dataEvent.tree_->GetEntry(evt);

    if(dataEvent.lep1_.Pt() < 1.0) continue;

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;

    bool passCuts[3][nSelTypes] = {{false, false},
                                   {false, false}};

    double theMET = dataEvent.met_; double theMETPHI = dataEvent.metPhi_; 

    
    double massZMin = trilepton_info(0,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                       dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                       dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				       dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
    int lType = 1;
    //if     (dataEvent.lq1_ * dataEvent.lq2_ < 0) lType = 0;

    double deltaRlJMin = 999.0;
    if(DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta());
    if(DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta());
    if(DeltaR(dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta(),dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta(),dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta());
    if(DeltaR(dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta()) < deltaRlJMin) deltaRlJMin = DeltaR(dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta());

    double zeppenfeld = TMath::Min(TMath::Max(TMath::Abs(dataEvent.lep1_.Eta()-(dataEvent.jet1_.Eta()+dataEvent.jet2_.Eta())/2.),
                                              TMath::Abs(dataEvent.lep2_.Eta()-(dataEvent.jet1_.Eta()+dataEvent.jet2_.Eta())/2.)),3.999);
    //int centrality = 0;
    //if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta()+0.1 > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
    //    (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta()+0.1 > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
    //   ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta()+0.1 > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
    //    (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta()+0.1 > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 40.0; if(dataEvent.type_ != SmurfTree::ee) metMin = 40.0;

    int newId=int(dataEvent.jet1McId_);
    //int tauId=int((dataEvent.jet1McId_%100-dataEvent.jet1McId_%10)/10);
    int qDisAgree=int((newId%1000-newId%100)/100);
    //int hasZCand=int(newId/1000);
    //int trackSel[4] = {int((dataEvent.jet2McId_%100-dataEvent.jet2McId_%10)/10),int((dataEvent.jet2McId_%1000-dataEvent.jet2McId_%100)/100),int((dataEvent.jet2McId_%10000-dataEvent.jet2McId_%1000)/1000),int(dataEvent.jet2McId_/10000)};

    bool passNjets    = dataEvent.njets_ >= 2;
    bool passMET      = dataEvent.met_ > metMin;
    bool preselCuts   = dataEvent.lep1_.Pt() > 20. && dataEvent.lep2_.Pt() > 20.;
    bool passBtagVeto = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool passVBFSel   = (dataEvent.jet1_+dataEvent.jet2_).M() > 500 && TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 2.5;
    // quick way to accept W+W+ (1) or W-W- (2)
    bool passNsignSel = true;
    if(TMath::Abs(dataEvent.lq1_+dataEvent.lq2_+dataEvent.lq3_) != 1) passNsignSel = false;

    bool passLSel = false;
    if     (lSel == 0 && dataEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && dataEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && dataEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && dataEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                     passLSel = true;
    else if(lSel == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (dataEvent.type_ == SmurfTree::me || dataEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && dataEvent.dilep_.M() > 15.0 && qDisAgree == 0) {
       
       if(passNsignSel && !passBtagVeto && passVBFSel == true &&  dataEvent.lid3_ != 0 && (dataEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection &&                    dataEvent.lep3_.Pt() > 10.) passCuts[lType][BTAGSEL] = true;
       if(passNsignSel &&  passBtagVeto && passVBFSel == true &&  dataEvent.lid3_ != 0 && (dataEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection && massZMin < 15.0 && dataEvent.lep3_.Pt() > 10.) passCuts[lType][WZSEL] = true;

    }

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && dataEvent.dilep_.M() > 15.0) {

      double outputVar[16];
      makeSystematicEffects3l(dataEvent.lid1_, dataEvent.lid2_, dataEvent.lid3_, dataEvent.lep1_, dataEvent.lep2_, dataEvent.lep3_, dataEvent.dilep_, 
                              dataEvent.mt_, theMET, theMETPHI, 
                              dataEvent.trackMet_, dataEvent.trackMetPhi_, 
			      dataEvent.njets_, dataEvent.jet1_, dataEvent.jet2_, 
			      year, 3, outputVar);
      double MVAVar[1] = {outputVar[13]};
      if     (thePlot == 0) {MVAVar[0]=outputVar[14];}
      else if(thePlot == 2) {MVAVar[0]=outputVar[ 0];}
      else if(thePlot == 9) {MVAVar[0]=outputVar[ 2];}
      else if(thePlot ==19) {MVAVar[0]=outputVar[12];}

      for(int nv=0; nv<1; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBin]-0.001);

      double MT3lTotx = dataEvent.met_*cos(dataEvent.metPhi_)-dataEvent.lep1_.Px()-dataEvent.lep2_.Px()-dataEvent.lep3_.Px();
      double MT3lToty = dataEvent.met_*sin(dataEvent.metPhi_)-dataEvent.lep1_.Py()-dataEvent.lep2_.Py()-dataEvent.lep3_.Py();
      double MT3lTot = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);

      if(passCuts[1][WZSEL]){ // begin making plots
	double myVar = -1.0;
	if     (thePlot == 0) myVar = TMath::Max(TMath::Min((dataEvent.jet1_+dataEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot == 1) myVar = TMath::Max(TMath::Min((dataEvent.lep1_+dataEvent.lep2_+dataEvent.jet1_+dataEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot == 2) myVar = TMath::Min(dataEvent.lep1_.Pt(),499.999);
	else if(thePlot == 3) myVar = TMath::Min(dataEvent.lep2_.Pt(),399.999);
	else if(thePlot == 4) myVar = dataEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = dataEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = dataEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = dataEvent.mt_;
	else if(thePlot == 8) myVar = dataEvent.dilep_.Pt();
	else if(thePlot == 9) myVar = TMath::Min(dataEvent.dilep_.M(),499.999);
	else if(thePlot ==10) myVar = dataEvent.njets_;
	else if(thePlot ==11) myVar = dataEvent.nvtx_;
	else if(thePlot ==12) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==13) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==14) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==15) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
	else if(thePlot ==16) myVar = dataEvent.type_;
	else if(thePlot ==17) myVar = dataEvent.dR_;
	else if(thePlot ==18) myVar = zeppenfeld;
	else if(thePlot ==19) myVar = TMath::Min(MT3lTot,999.999);
	else if(thePlot ==20) myVar = TMath::Min(massZMin,99.999);
	else if(thePlot ==21) myVar = TMath::Min((double)dataEvent.met_,199.999);
	else if(thePlot ==22) myVar = TMath::Min(deltaRlJMin,4.999);
	else assert(0);
      	histo6->Fill(myVar,1.0);
      } // end making plots

      if(passCuts[1][WZSEL]){
	histo_Data->Fill(MVAVar[0], 1.0);
      }

      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            nSelectedData[i+j*nSelTypes]  += 1.0;
          }
        }
      }

    } // if passCuts
  } // End loop data
  } // true if we want to look at the data

  char output[200];
  sprintf(output,Form("histo_nice%s.root",ECMsb.Data()));	 
  TFile* outFilePlotsNote = new TFile(output,"recreate");

  printf("gen_eff: %f / %f / %f / %f | rec_eff: %f / %f = %f | %f / %f = %f\n",genLevelNorm[0],genLevelNorm[1],genLevelNorm[2],genLevelNorm[3],
         genLevelNorm[4],genLevelNorm[5],genLevelNorm[4]/genLevelNorm[5],genLevelNorm[6],genLevelNorm[7],genLevelNorm[6]/genLevelNorm[7]);

  outFilePlotsNote->cd();
    double nOldH[6] = {histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights(),histo5->GetSumOfWeights()};
    for(int i=1; i<=histo0->GetNbinsX(); i++){
      if(histo0->GetBinContent(i) < 0) {histo0->SetBinContent(i,0.000001);histo0->SetBinError(i,0.000001);}
      if(histo1->GetBinContent(i) < 0) {histo1->SetBinContent(i,0.000001);histo1->SetBinError(i,0.000001);}
      if(histo2->GetBinContent(i) < 0) {histo2->SetBinContent(i,0.000001);histo2->SetBinError(i,0.000001);}
      if(histo3->GetBinContent(i) < 0) {histo3->SetBinContent(i,0.000001);histo3->SetBinError(i,0.000001);}
      if(histo4->GetBinContent(i) < 0) {histo4->SetBinContent(i,0.000001);histo4->SetBinError(i,0.000001);}
      if(histo5->GetBinContent(i) < 0) {histo5->SetBinContent(i,0.000001);histo5->SetBinError(i,0.000001);}
    }
    if(nOldH[0] > 0) histo0->Scale(nOldH[0]/histo0->GetSumOfWeights());
    if(nOldH[1] > 0) histo1->Scale(nOldH[1]/histo1->GetSumOfWeights());
    if(nOldH[2] > 0) histo2->Scale(nOldH[2]/histo2->GetSumOfWeights());
    if(nOldH[3] > 0) histo3->Scale(nOldH[3]/histo3->GetSumOfWeights());
    if(nOldH[4] > 0) histo4->Scale(nOldH[4]/histo4->GetSumOfWeights());
    if(nOldH[5] > 0) histo5->Scale(nOldH[5]/histo5->GetSumOfWeights());

    printf("histo -> d: %8.2f b: %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f | h: %8.2f\n",histo6->GetSumOfWeights(),
    histo0->GetSumOfWeights()+histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights()+histo5->GetSumOfWeights(),
    histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights(),histo5->GetSumOfWeights(),
    histo7->GetSumOfWeights());

    histo0->Write();
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();
    histo6->Write();
    histo7->Write();

  outFilePlotsNote->Close();
  
  double nOldWjets = histo_Wjets->GetSumOfWeights();
  for(int i=1; i<=histo_Wjets->GetNbinsX(); i++){
    if(histo_Wjets->GetBinContent(i)     < 0) {histo_Wjets    ->SetBinContent(i,0.0);histo_Wjets    ->SetBinError(i,0.0);}
    if(histo_Wjets_WUp->GetBinContent(i) < 0) {histo_Wjets_WUp->SetBinContent(i,0.0);histo_Wjets_WUp->SetBinError(i,0.0);}
  }
  if(nOldWjets > 0){
    histo_Wjets    ->Scale(nOldWjets/histo_Wjets   ->GetSumOfWeights());
    histo_Wjets_WUp->Scale(nOldWjets/histo_Wjets_WUp->GetSumOfWeights());
  }
  else {
    histo_Wjets    ->Scale(0.0);
    histo_Wjets_WUp->Scale(0.0);
  }
  if(showSignalOnly == false) printf("WjetsNorm ini/end: %f / %f\n",nOldWjets,histo_Wjets->GetSumOfWeights());
  double WjetsSyst = 1.36;

  const unsigned int nBkg = 4;
  double nTot[nSelTypes*2]; double nETot[nSelTypes*2];
  double bgdCombined[nSelTypes*2][nBkg],bgdCombinedE[nSelTypes*2][nBkg];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) {bgdCombined[i][j] = 0.0; bgdCombinedE[i][j] = 0.0;}
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("selection: %s\n",selTypeName[i].Data());
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("data(%2d): %f\n",i,nSelectedData[i]);
    nTot[i] = 0.0; nETot[i] = 0.0;
    for(int j=0; j<50; j++){
      // WWdps treatment
      if(j == 29 && bgdDecay[i][j] < 0) {printf("negative(29,%d) = %f +/- %f\n",i,bgdDecay[i][j],sqrt(weiDecay[i][j]));bgdDecay[i][j] = 0; weiDecay[i][j] = 0;}

      if(showSignalOnly == false || i%nSelTypes == WZSEL) if(bgdDecay[i][j] != 0) printf("bdg(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecay[i][j],sqrt(weiDecay[i][j]));

      if(j == 44) continue;

      nTot[i]  += bgdDecay[i][j];
      nETot[i] += weiDecay[i][j];

      if     (j == 27)			       {bgdCombined[i][0] += bgdDecay[i][j]; bgdCombinedE[i][0] += weiDecay[i][j];}
      else if(j == 30 || j == 28 || j == 29 || j == 31 ||
              j ==  5 || j == 13 || j == 20 || 
	      j == 10 || j ==  9 || j == 19)   {bgdCombined[i][1] += bgdDecay[i][j]; bgdCombinedE[i][1] += weiDecay[i][j];}
      else if(j == 21)			       {bgdCombined[i][2] += bgdDecay[i][j]; bgdCombinedE[i][2] += weiDecay[i][j];}
      else if(j == 1 || j == 23 || j == 45)    {bgdCombined[i][3] += bgdDecay[i][j]; bgdCombinedE[i][3] += weiDecay[i][j];}
    }
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTot[i],sqrt(nETot[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]-nTot[i])/nTot[i] > 0.00001) 
                    {printf("%f\n",bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]);assert(0);}
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("------\n");
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("bgd(xWZ) = %11.3f +/- %8.3f\n",bgdCombined[i][0],sqrt(bgdCombinedE[i][0]));
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("bgd(xZZ) = %11.3f +/- %8.3f\n",bgdCombined[i][1],sqrt(bgdCombinedE[i][1]));
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("bgd(VVV) = %11.3f +/- %8.3f\n",bgdCombined[i][2],sqrt(bgdCombinedE[i][2]));
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("bgd(xWj) = %11.3f +/- %8.3f\n",bgdCombined[i][3],sqrt(bgdCombinedE[i][3]));
    if(showSignalOnly == false || i%nSelTypes == WZSEL) printf("*******************************\n");
  }

  if(showSignalOnly == false) printf("WjetsNorm: %f --> %f\n",bgdCombined[WZSEL+nSelTypes][3],histo_Wjets->GetSumOfWeights());

  if(showSignalOnly == false) printf("+++++++++++++++++++++++++++++++\n");
  double nTotSyst[nSelTypesSyst*2]; double nETotSyst[nSelTypesSyst*2];
  double bgdCombinedSyst[nSelTypesSyst*2][nBkg],bgdCombinedESyst[nSelTypesSyst*2][nBkg];
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    if(showSignalOnly == false) printf("selectionSyst: %s\n",selTypeNameSyst[i].Data());
    for(unsigned int j=0; j<nBkg; j++) {bgdCombinedSyst[i][j] = 0.0; bgdCombinedESyst[i][j] = 0.0;}
    nTotSyst[i] = 0.0; nETotSyst[i] = 0.0;
    for(int j=0; j<50; j++){
      if(showSignalOnly == false) if(bgdDecaySyst[i][j] != 0) printf("bdgSyst(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecaySyst[i][j],sqrt(weiDecaySyst[i][j]));
      
      if(j == 44) continue;

      nTotSyst[i]  += bgdDecaySyst[i][j];
      nETotSyst[i] += weiDecaySyst[i][j];

     if     (j == 27)			      {bgdCombinedSyst[i][0] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][0] += weiDecaySyst[i][j];}
     else if(j == 30 || j == 28 || j == 29 || j == 31 ||		      
     	     j ==  5 || j == 13 || j == 20 || 
     	     j == 10 || j ==  9 || j == 19)   {bgdCombinedSyst[i][1] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][1] += weiDecaySyst[i][j];}
     else if(j == 21)			      {bgdCombinedSyst[i][2] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][2] += weiDecaySyst[i][j];}
     else if(j == 1 || j == 23 || j == 45)    {bgdCombinedSyst[i][3] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][3] += weiDecaySyst[i][j];}
    }
    if(showSignalOnly == false) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTotSyst[i],sqrt(nETotSyst[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]-nTotSyst[i])/nTotSyst[i] > 0.00001) 
                    {printf("%f\n",bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]);assert(0);}
    if(showSignalOnly == false) printf("------\n");
    if(showSignalOnly == false) printf("bgdSyst(xWZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][0],sqrt(bgdCombinedESyst[i][0]));
    if(showSignalOnly == false) printf("bgdSyst(xZZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][1],sqrt(bgdCombinedESyst[i][1]));
    if(showSignalOnly == false) printf("bgdSyst(VVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][2],sqrt(bgdCombinedESyst[i][2]));
    if(showSignalOnly == false) printf("bgdSyst(xWj) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][3],sqrt(bgdCombinedESyst[i][3]));
    if(showSignalOnly == false) printf("*******************************\n");
  }

  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) if(bgdCombined[i][j] == 0) {bgdCombined[i][j] = 0.0000000001; bgdCombinedE[i][j] = 0.0;}
  }

  double systEffect[nSelTypesSyst][nBkg];
  for(unsigned int i=0 ; i<nSelTypesSyst; i++){
    for(unsigned int j=0 ; j<nBkg; j++){
      if(bgdCombinedE[WZSEL+nSelTypes][j] > 0){
        systEffect[i][j] = bgdCombinedSyst[i+nSelTypesSyst][j]/bgdCombined[WZSEL+nSelTypes][j];
        if(systEffect[i][j] < 1) systEffect[i][j] = 1.0/systEffect[i][j];
      } else {systEffect[i][j] = 1.0;}
    }
  }
  if(showSignalOnly == false) printf("Syst(xWZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][0]-1,systEffect[JESDOWN][0]-1,systEffect[LEPP][0]-1,systEffect[LEPM][0]-1,systEffect[MET][0]-1,systEffect[EFFP][0]-1,systEffect[EFFM][0]-1);
  if(showSignalOnly == false) printf("Syst(xZZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][1]-1,systEffect[JESDOWN][1]-1,systEffect[LEPP][1]-1,systEffect[LEPM][1]-1,systEffect[MET][1]-1,systEffect[EFFP][1]-1,systEffect[EFFM][1]-1);
  if(showSignalOnly == false) printf("Syst(VVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][2]-1,systEffect[JESDOWN][2]-1,systEffect[LEPP][2]-1,systEffect[LEPM][2]-1,systEffect[MET][2]-1,systEffect[EFFP][2]-1,systEffect[EFFM][2]-1);
  if(showSignalOnly == false) printf("Syst(xWj) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][3]-1,systEffect[JESDOWN][3]-1,systEffect[LEPP][3]-1,systEffect[LEPM][3]-1,systEffect[MET][3]-1,systEffect[EFFP][3]-1,systEffect[EFFM][3]-1);

  double pdf_qqbar[2] = {1.071,1.065};
  for(int i=1; i<=histo_WZ->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_WZ_WZStatUp	      	    ->SetBinContent(i,TMath::Max(histo_WZ    	->GetBinContent(i)+factorUp  *histo_WZ   	 ->GetBinError(i),0.000001));
    histo_WZ_WZStatDown        	    ->SetBinContent(i,TMath::Max(histo_WZ    	->GetBinContent(i)+factorDown*histo_WZ   	 ->GetBinError(i),0.000001));
    histo_ZZ_ZZStatUp	      	    ->SetBinContent(i,TMath::Max(histo_ZZ    	->GetBinContent(i)+factorUp  *histo_ZZ   	 ->GetBinError(i),0.000001));
    histo_ZZ_ZZStatDown        	    ->SetBinContent(i,TMath::Max(histo_ZZ    	->GetBinContent(i)+factorDown*histo_ZZ   	 ->GetBinError(i),0.000001));
    histo_VVV_VVVStatUp        	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorUp  *histo_VVV  	 ->GetBinError(i),0.000001));
    histo_VVV_VVVStatDown      	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorDown*histo_VVV  	 ->GetBinError(i),0.000001));
    histo_Wjets_WjetsStatUp    	    ->SetBinContent(i,TMath::Max(histo_Wjets    ->GetBinContent(i)+factorUp  *histo_Wjets        ->GetBinError(i),0.000001));
    histo_Wjets_WjetsStatDown  	    ->SetBinContent(i,TMath::Max(histo_Wjets    ->GetBinContent(i)+factorDown*histo_Wjets        ->GetBinError(i),0.000001));
    histo_Higgs_HiggsStatUp	    ->SetBinContent(i,TMath::Max(histo_Higgs    ->GetBinContent(i)+factorUp  *histo_Higgs	 ->GetBinError(i),0.000001));
    histo_Higgs_HiggsStatDown       ->SetBinContent(i,TMath::Max(histo_Higgs    ->GetBinContent(i)+factorDown*histo_Higgs	 ->GetBinError(i),0.000001));
  }
  double mean,up,diff;

  if(showSignalOnly == false) {
    printf("nuisance Wj: %f/%f\n",histo_Wjets->GetSumOfWeights(),histo_Wjets_WUp->GetSumOfWeights());
  }
  histo_Wjets_WUp->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_WUp->GetSumOfWeights());

  for(int i=1; i<=histo_WZ->GetNbinsX(); i++){
    // METRes
    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_ZZ			   ->GetBinContent(i);
    up   = histo_ZZ_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZZ_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZZ_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Higgs			      ->GetBinContent(i);
    up   = histo_Higgs_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Higgs_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Higgs_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    // LepRes
    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_ZZ			   ->GetBinContent(i);
    up   = histo_ZZ_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZZ_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZZ_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Higgs			      ->GetBinContent(i);
    up   = histo_Higgs_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Higgs_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Higgs_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    // JES
    mean = histo_WZ			->GetBinContent(i);
    up   = histo_WZ_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_ZZ			->GetBinContent(i);
    up   = histo_ZZ_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZZ_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZZ_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			 ->GetBinContent(i);
    up   = histo_VVV_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Higgs			   ->GetBinContent(i);
    up   = histo_Higgs_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Higgs_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Higgs_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    // GEN
    mean = histo_Wjets 		         ->GetBinContent(i);
    up   = histo_Wjets_WUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Wjets_WDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Wjets_WDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
  }
  histo_Wjets_WDown->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_WDown->GetSumOfWeights());

  //----------------------------------------------------------------------------
  // Produce output cards for shape-based analyses
  //----------------------------------------------------------------------------
  if(showSignalOnly == false){
  char outputLimits[200];
  sprintf(outputLimits,"qqwz%2s.input_%4s.root",finalStateName,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data	 ->Write();
  histo_WZ	 ->Write();
  histo_ZZ	 ->Write();
  histo_VVV	 ->Write();
  histo_Wjets	 ->Write();
  histo_Higgs	 ->Write();

  cout << histo_Data	 ->GetSumOfWeights() << " ";
  cout << histo_WZ	 ->GetSumOfWeights() << " ";
  cout << histo_ZZ	 ->GetSumOfWeights() << " ";
  cout << histo_VVV	 ->GetSumOfWeights() << " ";
  cout << histo_Wjets	 ->GetSumOfWeights() << " ";
  cout << histo_Higgs	 ->GetSumOfWeights() << " ";
  cout << endl;


  printf("uncertainties Stat\n");
  histo_WZ_WZStatUp	  	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_WZStatUp	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_WZStatDown	 	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_WZStatDown	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_ZZStatUp	  	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_ZZStatUp	->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_ZZStatDown	 	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_ZZStatDown	->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_VVVStatUp	  	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_VVVStatUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_VVVStatDown   	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_VVVStatDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_WjetsStatUp  	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WjetsStatUp  ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_WjetsStatDown	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WjetsStatDown->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_HiggsStatUp	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_HiggsStatUp  ->GetBinContent(i)/histo_Higgs->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_HiggsStatDown	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_HiggsStatDown->GetBinContent(i)/histo_Higgs->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEff\n");
  histo_WZ_LepEffUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepEffUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_LepEffDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepEffDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_LepEffUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_LepEffUp     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_LepEffDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_LepEffDown   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepEffUp          ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepEffUp  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepEffDown        ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepEffDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_LepEffUp        ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_LepEffUp        ->GetBinContent(i)/histo_Higgs	     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_LepEffDown      ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_LepEffDown   ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LetRes\n");
  histo_WZ_LepResUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepResUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_LepResDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepResDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_LepResUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_LepResUp     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_LepResDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_LepResDown   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepResUp          ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepResUp  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepResDown        ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepResDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_LepResUp        ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_LepResUp        ->GetBinContent(i)/histo_Higgs	     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_LepResDown      ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_LepResDown   ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties METRes\n");
  histo_WZ_METResUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_METResUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_METResDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_METResDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_METResUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_METResUp     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_METResDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_METResDown   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_METResUp          ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_METResUp  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_METResDown        ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_METResDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_METResUp        ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs	 ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_METResUp	  ->GetBinContent(i)/histo_Higgs	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_METResDown      ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs	 ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_METResDown ->GetBinContent(i)/histo_Higgs	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_WZ_JESUp              ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_JESUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_JESDown            ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_JESDown	   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_JESUp              ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_JESUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_JESDown            ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_JESDown	   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_JESUp             ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_JESUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_JESDown           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_JESDown	  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_JESUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_JESUp   ->GetBinContent(i)/histo_Higgs	     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_JESDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_JESDown   ->GetBinContent(i)/histo_Higgs      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JER\n");
  histo_WZ_JERUp              ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_JERUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_JERDown            ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_JERDown	   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_JERUp              ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_JERUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_JERDown            ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_JERDown	   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_JERUp             ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_JERUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_JERDown           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_JERDown	  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_JERUp           ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_JERUp   ->GetBinContent(i)/histo_Higgs	     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_JERDown         ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_JERDown   ->GetBinContent(i)/histo_Higgs      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties GEN\n");
  histo_Wjets_WUp	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WUp       ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_WDown	  ->Write(); for(int i=1; i<=histo_WZ->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WDown     ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  for(int nb=1; nb<=nBin; nb++){

    double systNLO[3] = {1.10,1.10,1.0}; // WZ, ZZ, Wjets
    if     (histo_Wjets->GetBinContent(nb) > 0 && histo_Wjets_WUp  ->GetBinContent(nb) > 0) systNLO[2] = histo_Wjets_WUp  ->GetBinContent(nb)/histo_Wjets->GetBinContent(nb);
    else if(histo_Wjets->GetBinContent(nb) > 0 && histo_Wjets_WDown->GetBinContent(nb) > 0) systNLO[2] = histo_Wjets->GetBinContent(nb)/histo_Wjets_WDown->GetBinContent(nb);

    double systEff[4] = {1.0,1.0,1.0,1.0};
    if(histo_WZ   ->GetBinContent(nb) > 0 && histo_WZ_LepEffUp	->GetBinContent(nb) > 0)  systEff[0] =  histo_WZ_LepEffUp   ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    if(histo_ZZ   ->GetBinContent(nb) > 0 && histo_ZZ_LepEffUp	->GetBinContent(nb) > 0)  systEff[1] =  histo_ZZ_LepEffUp   ->GetBinContent(nb)/histo_ZZ   ->GetBinContent(nb);
    if(histo_VVV  ->GetBinContent(nb) > 0 && histo_VVV_LepEffUp	->GetBinContent(nb) > 0)  systEff[2] =  histo_VVV_LepEffUp  ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_LepEffUp->GetBinContent(nb) > 0) systEff[3] =  histo_Higgs_LepEffUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);

    double systLep[4] = {1.0,1.0,1.0,1.0};
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_LepResUp     ->GetBinContent(nb) > 0) systLep[0] =  histo_WZ_LepResUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_LepResDown   ->GetBinContent(nb) > 0) systLep[0] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_LepResDown   ->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_LepResUp     ->GetBinContent(nb) > 0) systLep[1] =  histo_ZZ_LepResUp     ->GetBinContent(nb)/histo_ZZ   ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_LepResDown   ->GetBinContent(nb) > 0) systLep[1] =  histo_ZZ   ->GetBinContent(nb)/histo_ZZ_LepResDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_LepResUp    ->GetBinContent(nb) > 0) systLep[2] =  histo_VVV_LepResUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_LepResDown  ->GetBinContent(nb) > 0) systLep[2] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_LepResDown  ->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_LepResUp  ->GetBinContent(nb) > 0) systLep[3] =  histo_Higgs_LepResUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_LepResDown->GetBinContent(nb) > 0) systLep[3] =  histo_Higgs->GetBinContent(nb)/histo_Higgs_LepResDown->GetBinContent(nb);

    double systMet[4] = {1.0,1.0,1.0,1.0};
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_METResUp     ->GetBinContent(nb) > 0) systMet[0] =  histo_WZ_METResUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_METResDown   ->GetBinContent(nb) > 0) systMet[0] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_METResDown   ->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_METResUp     ->GetBinContent(nb) > 0) systMet[1] =  histo_ZZ_METResUp     ->GetBinContent(nb)/histo_ZZ   ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_METResDown   ->GetBinContent(nb) > 0) systMet[1] =  histo_ZZ   ->GetBinContent(nb)/histo_ZZ_METResDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_METResUp    ->GetBinContent(nb) > 0) systMet[2] =  histo_VVV_METResUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_METResDown  ->GetBinContent(nb) > 0) systMet[2] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_METResDown  ->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_METResUp  ->GetBinContent(nb) > 0) systMet[3] =  histo_Higgs_METResUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_METResDown->GetBinContent(nb) > 0) systMet[3] =  histo_Higgs->GetBinContent(nb)/histo_Higgs_METResDown->GetBinContent(nb);

    double systJes[4] = {1.0,1.0,1.0,1.0};
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_JESUp     ->GetBinContent(nb) > 0) systJes[0] =  histo_WZ_JESUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_JESDown   ->GetBinContent(nb) > 0) systJes[0] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_JESDown   ->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_JESUp     ->GetBinContent(nb) > 0) systJes[1] =  histo_ZZ_JESUp     ->GetBinContent(nb)/histo_ZZ   ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_JESDown   ->GetBinContent(nb) > 0) systJes[1] =  histo_ZZ   ->GetBinContent(nb)/histo_ZZ_JESDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_JESUp    ->GetBinContent(nb) > 0) systJes[2] =  histo_VVV_JESUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_JESDown  ->GetBinContent(nb) > 0) systJes[2] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_JESDown  ->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_JESUp  ->GetBinContent(nb) > 0) systJes[3] =  histo_Higgs_JESUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_JESDown->GetBinContent(nb) > 0) systJes[3] =  histo_Higgs->GetBinContent(nb)/histo_Higgs_JESDown->GetBinContent(nb);

    double systJER[4] = {1.0,1.0,1.0,1.0};
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_JERUp     ->GetBinContent(nb) > 0) systJER[0] =  histo_WZ_JERUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_JERDown   ->GetBinContent(nb) > 0) systJER[0] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_JERDown   ->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_JERUp     ->GetBinContent(nb) > 0) systJER[1] =  histo_ZZ_JERUp     ->GetBinContent(nb)/histo_ZZ   ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)    > 0 && histo_ZZ_JERDown   ->GetBinContent(nb) > 0) systJER[1] =  histo_ZZ   ->GetBinContent(nb)/histo_ZZ_JERDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_JERUp    ->GetBinContent(nb) > 0) systJER[2] =  histo_VVV_JERUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_JERDown  ->GetBinContent(nb) > 0) systJER[2] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_JERDown  ->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_JERUp  ->GetBinContent(nb) > 0) systJER[3] =  histo_Higgs_JERUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_JERDown->GetBinContent(nb) > 0) systJER[3] =  histo_Higgs->GetBinContent(nb)/histo_Higgs_JERDown->GetBinContent(nb);

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_qqwz%2s_shape_%4s_Bin%d.txt",finalStateName,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin qqwz%2s%4s%d qqwz%2s%4s%d qqwz%2s%4s%d qqwz%2s%4s%d qqwz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process WZ ZZ VVV Wjets Higgs\n");
    if(histo_Higgs->GetSumOfWeights() <=0)
    newcardShape << Form("process 0 1 2 3 4\n");
    else
    newcardShape << Form("process 1 2 3 4 0\n");
    newcardShape << Form("rate %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_VVV->GetBinContent(nb),histo_Wjets->GetBinContent(nb),histo_Higgs->GetBinContent(nb));
    newcardShape << Form("lumi_%4s                                 lnN %5.3f %5.3f %5.3f  -   %5.3f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE);		 
    newcardShape << Form("%s                                       lnN %5.3f %5.3f %5.3f  -   %5.3f\n",effName,systEff[0],systEff[1],systEff[2],systEff[3]);
    newcardShape << Form("%s                                       lnN %5.3f %5.3f %5.3f  -   %5.3f\n",momName,systLep[0],systLep[1],systLep[2],systLep[3]);
    newcardShape << Form("CMS_scale_met                            lnN %5.3f %5.3f %5.3f  -   %5.3f\n",systMet[0],systMet[1],systMet[2],systMet[3]);
    newcardShape << Form("CMS_scale_j                              lnN %5.3f %5.3f %5.3f  -   %5.3f\n",systJes[0],systJes[1],systJes[2],systJes[3]);		    
    newcardShape << Form("CMS_res_j                                lnN %5.3f %5.3f %5.3f  -   %5.3f\n",systJER[0],systJER[1],systJER[2],systJER[3]);		    
    newcardShape << Form("CMS_eff_b                                lnN %5.3f %5.3f %5.3f  -   %5.3f\n",1.02,1.02,1.02,1.02);
    newcardShape << Form("pdf_qqbar                                lnN %5.3f   -     -    -   %5.3f\n",pdf_qqbar[0],pdf_qqbar[1]);
    newcardShape << Form("QCDscale_Higgs                           lnN   -     -     -    -    1.03\n");    
    newcardShape << Form("CMS_qqwz_QCDWZ                           lnN %5.3f   -     -    -	-  \n",systNLO[0]);
    newcardShape << Form("CMS_qqwz_QCDZZ                           lnN  -    %5.3f   -    -	-  \n",systNLO[1]);	    
    newcardShape << Form("QCDscale_VVV		                   lnN  -      -   1.500  -	-  \n");       
    newcardShape << Form("CMS_FakeRate                             lnN  -      -     -   %5.3f  -  \n",WjetsSyst);  
    newcardShape << Form("CMS_qqwz_MVAW                            lnN  -      -     -   %5.3f  -  \n",systNLO[2]);
    if(histo_WZ->GetBinContent(nb) > 0.0001)
    newcardShape << Form("CMS_qqwz%s_MVAWZStat_%s_Bin%d            lnN  %5.3f	-     -    -	 -  \n",finalStateName,ECMsb.Data(),nb-1,histo_WZ_WZStatUp    ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb));
    if(histo_ZZ->GetBinContent(nb) > 0.0001)
    newcardShape << Form("CMS_qqwz%s_MVAZZStat_%s_Bin%d            lnN   -    %5.3f   -    -	 -  \n",finalStateName,ECMsb.Data(),nb-1,histo_ZZ_ZZStatUp    ->GetBinContent(nb)/histo_ZZ   ->GetBinContent(nb));
    if(histo_VVV->GetBinContent(nb) > 0.0001)
    newcardShape << Form("CMS_qqwz%s_MVAVVVStat_%s_Bin%d           lnN   -	-   %5.3f  -	 -  \n",finalStateName,ECMsb.Data(),nb-1,histo_VVV_VVVStatUp	->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb));
    if(histo_Wjets->GetBinContent(nb) > 0.0001)
    newcardShape << Form("CMS_qqwz%s_MVAWjetsStat_%s_Bin%d         lnN   -	-     -   %5.3f  -  \n",finalStateName,ECMsb.Data(),nb-1,histo_Wjets_WjetsStatUp->GetBinContent(nb)/histo_Wjets->GetBinContent(nb));
    if(histo_Higgs->GetBinContent(nb) > 0.0001)
    newcardShape << Form("CMS_qqwz%s_MVAHiggsStat_%s_Bin%d         lnN   -	-     -    -   %5.3f\n",finalStateName,ECMsb.Data(),nb-1,histo_Higgs_HiggsStatUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb));
    newcardShape.close();
  }
  } // if showSignalOnly == true

  return;
}

void scaleFactor_WS(LorentzVector l,int lq, int ld, int mcld, double val[2], int opt){
//---------------------------------------------------------------------
// |eta|        data                  mc                    factor
//---------------------------------------------------------------------
//0.0-0.5 0.000109 +/- 0.000017 | 0.000106 +/- 0.000017 ==> 1.028 +/- 0.160
//0.5-1.0 0.000210 +/- 0.000025 | 0.000188 +/- 0.000025 ==> 1.117 +/- 0.133
//1.0-1.5 0.001302 +/- 0.000087 | 0.001022 +/- 0.000087 ==> 1.274 +/- 0.085
//1.5-2.0 0.003437 +/- 0.000193 | 0.002562 +/- 0.000193 ==> 1.342 +/- 0.075
//2.0-2.5 0.003270 +/- 0.000241 | 0.002245 +/- 0.000240 ==> 1.456 +/- 0.107
// additional 10% uncertainty for the overall normalization
  if(opt == 0) {
    double factor[5]  = {1.028,1.117,1.274,1.342,1.456};
    double factorE[5] = {0.160,0.133,0.085,0.075,0.107};

    if(abs(ld) == 11){
      if((mcld ==  11 && lq > 0) || 
	 (mcld == -11 && lq < 0)){ // wrong charge
	if     (abs(l.Eta()) >= 0.0 && abs(l.Eta()) < 0.5) {val[0] = val[0]*factor[0]; val[1] = val[1]*(factor[0]+sqrt(factorE[0]*factorE[0]+0.10*0.10));}
	else if(abs(l.Eta()) >= 0.5 && abs(l.Eta()) < 1.0) {val[0] = val[0]*factor[1]; val[1] = val[1]*(factor[1]+sqrt(factorE[1]*factorE[1]+0.10*0.10));}
	else if(abs(l.Eta()) >= 1.0 && abs(l.Eta()) < 1.5) {val[0] = val[0]*factor[2]; val[1] = val[1]*(factor[2]+sqrt(factorE[2]*factorE[2]+0.10*0.10));}
	else if(abs(l.Eta()) >= 1.5 && abs(l.Eta()) < 2.0) {val[0] = val[0]*factor[3]; val[1] = val[1]*(factor[3]+sqrt(factorE[3]*factorE[3]+0.10*0.10));}
	else if(abs(l.Eta()) >= 2.0 && abs(l.Eta()) < 2.5) {val[0] = val[0]*factor[4]; val[1] = val[1]*(factor[4]+sqrt(factorE[4]*factorE[4]+0.10*0.10));}
      }
    }
  } else {
    double factor[5]  = {0.000109,0.000210,0.001302,0.003437,0.003270};

    if(abs(ld) == 11){
	if     (abs(l.Eta()) >= 0.0 && abs(l.Eta()) < 0.5) {val[0] = val[0]*factor[0]; val[1] = 0.0;}
	else if(abs(l.Eta()) >= 0.5 && abs(l.Eta()) < 1.0) {val[0] = val[0]*factor[1]; val[1] = 0.0;}
	else if(abs(l.Eta()) >= 1.0 && abs(l.Eta()) < 1.5) {val[0] = val[0]*factor[2]; val[1] = 0.0;}
	else if(abs(l.Eta()) >= 1.5 && abs(l.Eta()) < 2.0) {val[0] = val[0]*factor[3]; val[1] = 0.0;}
	else if(abs(l.Eta()) >= 2.0 && abs(l.Eta()) < 2.5) {val[0] = val[0]*factor[4]; val[1] = 0.0;}
    }
  }
}
