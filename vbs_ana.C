#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/makeSystematicEffects.h"
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
//root -l -q -b vbs_ana.C+'(0,"ntuples_53x/backgroundA_skim8_lt012.root","ntuples_53x/data_skim8.root","ntuples_53x/hww_syst_skim8.root",3,4)';
//root -l -q -b vbs_ana.C+'(0,"ntuples_53x/backgroundA_skim8_lt012.root","ntuples_53x/data_skim8.root","ntuples_53x/hww_syst_skim8.root",3,14)';
//root -l -q -b vbs_ana.C+'(0,"ntuples_53x/backgroundA_skim8_lt012.root","ntuples_53x/data_skim8.root","ntuples_53x/hww_syst_skim8.root",3,24)'

//std::string file_for_grid="/afs/cern.ch/work/a/anlevin/data/lhe/qed_4_qcd_99_lt012_grid.lhe";
std::string file_for_grid="/afs/cern.ch/user/a/anlevin/public/forGuillelmo04Feb2014/unweighted_events_9.lhe";
int x_param_number = 11;
int y_param_number = 12;
std::vector<std::pair<float,float> > grid_points;
std::vector<float> histo_grid;
std::vector<int> lhe_weight_index;

const int verboseLevel =   1;
bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC
SmurfTree systEvent;
const unsigned int nSelTypes = 3;
const unsigned int nSelTypesSyst = 7;
const bool showSignalOnly = false;

enum selType {WWSEL, BTAGSEL, WZSEL};
TString selTypeName[nSelTypes*2] = {"WWSEL-OS", "BTAGSEL-OS", "WZSEL-OS",
                                    "WWSEL-SS", "BTAGSEL-SS", "WZSEL-SS"};
enum selTypeSyst {JESUP=0, JESDOWN, LEPP, LEPM, MET, EFFP, EFFM};
TString selTypeNameSyst[nSelTypesSyst*2] = {"JESUP-OS", "JESDOWN-OS", "LEPP-OS", "LEPM-OS", "MET-OS", "EFFP-OS", "EFFM-OS",
                                            "JESUP-SS", "JESDOWN-SS", "LEPP-SS", "LEPM-SS", "MET-SS", "EFFP-SS", "EFFM-SS"};

bool run_over_data = false;
bool doAQGCsAna = false;
bool use_anom_sample = true;
int which_lhe_weight = 60;


void scaleFactor_WS(LorentzVector l,int q, int ld, int mcld, double val[2]);

void parse_grid(string lhe_filename){
  grid_points.push_back(pair<float,float>(0,0));
  histo_grid.push_back(0);
  lhe_weight_index.push_back(0);

  ifstream infile(lhe_filename.c_str());
  assert(infile.is_open());

  while(!infile.eof()){
    std::string line;
    getline(infile,line);

    if(line=="<initrwgt>\0"){
      getline(infile,line);
      assert(line=="<weightgroup type='mg_reweighting'>");

      int i = 1;

      while(true){
	getline(infile,line);

	if(line=="</initrwgt>\0")
	  return;

	if (line == "</weight>\0" || line=="</weightgroup>\0")
	  continue;

	int param_number1 = -1;
	int param_number2 = -1;
	float param1 = 0;
	float param2 = 0;

	assert(line.find("set param_card anoinputs") != string::npos);
	std::string paraminfo1=line.substr(line.find("set param_card anoinputs ")+std::string("set param_card anoinputs ").size(),line.find("#")-line.find("set param_card anoinputs ")-std::string("set param_card anoinputs ").size());
	stringstream ss1;
	ss1 << paraminfo1;
	ss1 >> param_number1;
	if(param_number1 == x_param_number)
	  ss1 >> param1;
	else if (param_number1==y_param_number)
	  ss1 >> param2;
	//else
	//  assert(0);

	getline(infile,line);

	if (line != "</weight>\0"){

	  assert(line.find("set param_card anoinputs") != string::npos);
	  std::string paraminfo2=line.substr(line.find("set param_card anoinputs ")+std::string("set param_card anoinputs ").size(),line.find("#")-line.find("set param_card anoinputs ")-std::string("set param_card anoinputs ").size());
	  stringstream ss2;
	  ss2 << paraminfo2;
	  ss2 >> param_number2;
	  if(param_number2 == x_param_number)
	    ss2 >> param1;
	  else if (param_number2==y_param_number)
	    ss2 >> param2;
	  //else
	  //  assert(0);

	  assert(param_number1 != param_number2);

	}
	if((param_number1 == x_param_number && param_number2 == y_param_number) || (param_number2 == x_param_number && param_number1 == y_param_number)|| ( param_number1 == x_param_number && param_number2 == -1) || (param_number1 == y_param_number && param_number2 == -1)) {
	  
	  //the same grid point may happen multiple times
	  //make sure to only add each grid point once
	  bool found =false;
	  for(unsigned int j = 0; j < grid_points.size(); j++){
	    if (grid_points[j] == pair<float,float>(param1,param2))
	      found = true;
	  }
	  
	  if(!found){
	    grid_points.push_back(pair<float,float>(param1,param2));
	    histo_grid.push_back(0);
	    lhe_weight_index.push_back(i);
	  }
	}

	i++;

      }
    }
  }

  std::cout << "reweight block not found, exiting" << std::endl;
  exit(1);
}

void vbs_ana
(
 int thePlot = 0,
 TString bgdInputFile    = "ntuples_53x/backgroundA_skim8_lt012.root",
 TString dataInputFile   = "ntuples_53x/data_skim8.root",
 TString systInputFile   = "ntuples_53x/hww_syst_skim8.root",
 int period = 3,
 int lSel = 4
 )
{

  if(doAQGCsAna == true){
    parse_grid(file_for_grid);
    std::cout << "grid_points.size() = " << grid_points.size() << std::endl;
    for(unsigned int i = 0; i < grid_points.size(); i++){
      std::cout << grid_points[i].first << ", " << grid_points[i].second << std::endl;
    }
    //change to more convenient units  
    for(unsigned int i = 0; i < grid_points.size(); i++){
      grid_points[i].first = grid_points[i].first*pow(10.,11);
      grid_points[i].second = grid_points[i].second*pow(10.,11);
    }
    for(unsigned int i = 0; i < grid_points.size(); i++){
      std::cout << grid_points[i].first << ", " << grid_points[i].second << std::endl;
    }
    for(unsigned int i = 0; i < lhe_weight_index.size(); i++){
      std::cout << "lhe_weight_index[i] = " << lhe_weight_index[i] << std::endl;
    }
  }

  double frCorr = 0.78;
  double lumi = 1.0;
  double ptJetMin = 30.0;

  bool fCheckProblem = true;

  int signSel = 0; // both W+W+ and W-W-
  if     (              lSel < 10) {lSel = lSel -  0; signSel = 0;} // both W+W+ and W-W-
  else if(lSel >= 10 && lSel < 20) {lSel = lSel - 10; signSel = 1;} // only W+W+
  else if(lSel >= 20 && lSel < 30) {lSel = lSel - 20; signSel = 2;} // only W-W-
  else assert(0);

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  if(systInputFile != ""){
    systEvent.LoadTree(systInputFile,-1);
    systEvent.InitTree(0);
  }

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

  const int nBin = 4;
  Float_t xbins[nBin+1] = {700, 1100, 1500, 2000, 3000};
  if(thePlot == 0) {xbins[0] = 500; xbins[1] = 700; xbins[2] = 1100; xbins[3] = 1600; xbins[4] = 2000;}
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBin, xbins);
  histoMVA->Sumw2();
  TH1D *histo_Data      = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_WWewk     = (TH1D*) histoMVA->Clone("histo_WWewk");
  TH1D *histo_WWqcd     = (TH1D*) histoMVA->Clone("histo_WWqcd");
  TH1D *histo_WZ        = (TH1D*) histoMVA->Clone("histo_WZ");
  TH1D *histo_WS        = (TH1D*) histoMVA->Clone("histo_WS");
  TH1D *histo_VVV       = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_Wjets     = (TH1D*) histoMVA->Clone("histo_Wjets");

  std::vector<TH1D *> histo_WWewk_anom;

  if(doAQGCsAna == true){
    for(unsigned int a = 0; a < grid_points.size(); a++){
      stringstream ss;
      ss << a;
      histo_WWewk_anom.push_back((TH1D*) histoMVA->Clone("histo_WWewk_anom"+a));
    }
  }

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

  if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 8; xminPlot = 0.0; xmaxPlot =  2000;} // mjj
  else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 60; xminPlot = 0.0; xmaxPlot = 4000.0;} // mlljj
  else if(thePlot >=  2 && thePlot <=  9) {nBinPlot = 60; xminPlot = 0.0; xmaxPlot = 1500;} //mll
  else if(thePlot >= 10 && thePlot <= 10) {nBinPlot = 10;  xminPlot =  -0.5; xmaxPlot =  9.5;}
  else if(thePlot >= 11 && thePlot <= 11) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot >= 12 && thePlot <= 12) {nBinPlot = 36; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 13 && thePlot <= 14) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 15 && thePlot <= 15) {nBinPlot = 9; xminPlot = 0.0; xmaxPlot =  8.75;} // detajjs
  else if(thePlot >= 16 && thePlot <= 16) {nBinPlot = 4; xminPlot = -0.5; xmaxPlot = 3.5;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 44; xminPlot = 0.0; xmaxPlot = 4.4;}
  else assert(0);

  TH1D* histo0;
  if(thePlot != 0)   histo0 = new TH1D("histo0", "histo0", nBinPlot, xminPlot, xmaxPlot);
  else                              histo0 = new TH1D("histo0", "histo0", nBin, xbins);  
  histo0->Sumw2();
  TH1D* histo1 = (TH1D*) histo0->Clone("histo1");
  TH1D* histo2 = (TH1D*) histo0->Clone("histo2");
  TH1D* histo3 = (TH1D*) histo0->Clone("histo3");
  TH1D* histo4 = (TH1D*) histo0->Clone("histo4");
  TH1D* histo5 = (TH1D*) histo0->Clone("histo5");
  histo0->Scale(0.0);
  histo1->Scale(0.0);
  histo2->Scale(0.0);
  histo3->Scale(0.0);
  histo4->Scale(0.0);
  histo5->Scale(0.0);

  TH1D* histo_WWewk_WWewkStatUp   = new TH1D( Form("histo_WWewk_CMS_wwss%s__MVAWWewkStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WWewk_CMS_wwss%s__MVAWWewkStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_WWewk_WWewkStatUp  ->Sumw2();
  TH1D* histo_WWewk_WWewkStatDown = new TH1D( Form("histo_WWewk_CMS_wwss%s_MVAWWewkStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_WWewk_CMS_wwss%s__MVAWWewkStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_WWewk_WWewkStatDown->Sumw2();
  TH1D* histo_WWqcd_WWqcdStatUp   = new TH1D( Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_WWqcd_WWqcdStatUp  ->Sumw2();
  TH1D* histo_WWqcd_WWqcdStatDown = new TH1D( Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_WWqcd_WWqcdStatDown->Sumw2();
  TH1D* histo_WS_WSStatUp         = new TH1D( Form("histo_WS_CMS_wwss%s_MVAWSStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WS_CMS_wwss%s_MVAWSStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_WS_WSStatUp  ->Sumw2();
  TH1D* histo_WS_WSStatDown       = new TH1D( Form("histo_WS_CMS_wwss%s_MVAWSStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_WS_CMS_wwss%s_MVAWSStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_WS_WSStatDown->Sumw2();
  TH1D* histo_WZ_WZStatUp         = new TH1D( Form("histo_WZ_CMS_wwss%s_MVAWZStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wwss%s_MVAWZStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_WZ_WZStatUp  ->Sumw2();
  TH1D* histo_WZ_WZStatDown       = new TH1D( Form("histo_WZ_CMS_wwss%s_MVAWZStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wwss%s_MVAWZStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_WZ_WZStatDown->Sumw2();
  TH1D* histo_VVV_VVVStatUp       = new TH1D( Form("histo_VVV_CMS_wwss%s_MVAVVVStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwss%s_MVAVVVStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_VVV_VVVStatUp  ->Sumw2();
  TH1D* histo_VVV_VVVStatDown     = new TH1D( Form("histo_VVV_CMS_wwss%s_MVAVVVStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwss%s_MVAVVVStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_VVV_VVVStatDown->Sumw2();
  TH1D* histo_Wjets_WjetsStatUp   = new TH1D( Form("histo_Wjets_CMS_wwss%s_MVAWjetsStat_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_wwss%s_MVAWjetsStat_%sUp"  ,finalStateName,ECMsb.Data()), nBin, xbins); histo_Wjets_WjetsStatUp  ->Sumw2();
  TH1D* histo_Wjets_WjetsStatDown = new TH1D( Form("histo_Wjets_CMS_wwss%s_MVAWjetsStat_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_wwss%s_MVAWjetsStat_%sDown",finalStateName,ECMsb.Data()), nBin, xbins); histo_Wjets_WjetsStatDown->Sumw2();

  TH1D* histo_WWewk_LepEffUp      = new TH1D( Form("histo_WWewk_%sUp",effName)  , Form("histo_WWewk_%sUp",effName)  , nBin, xbins); histo_WWewk_LepEffUp  ->Sumw2();
  TH1D* histo_WWewk_LepEffDown    = new TH1D( Form("histo_WWewk_%sDown",effName), Form("histo_WWewk_%sDown",effName), nBin, xbins); histo_WWewk_LepEffDown->Sumw2();
  TH1D* histo_WWqcd_LepEffUp      = new TH1D( Form("histo_WWqcd_%sUp",effName)  , Form("histo_WWqcd_%sUp",effName)  , nBin, xbins); histo_WWqcd_LepEffUp  ->Sumw2();
  TH1D* histo_WWqcd_LepEffDown    = new TH1D( Form("histo_WWqcd_%sDown",effName), Form("histo_WWqcd_%sDown",effName), nBin, xbins); histo_WWqcd_LepEffDown->Sumw2();
  TH1D* histo_WZ_LepEffUp         = new TH1D( Form("histo_WZ_%sUp",effName)  , Form("histo_WZ_%sUp",effName)  , nBin, xbins); histo_WZ_LepEffUp  ->Sumw2();
  TH1D* histo_WZ_LepEffDown       = new TH1D( Form("histo_WZ_%sDown",effName), Form("histo_WZ_%sDown",effName), nBin, xbins); histo_WZ_LepEffDown->Sumw2();
  TH1D* histo_WS_LepEffUp         = new TH1D( Form("histo_WS_%sUp",effName)  , Form("histo_WS_%sUp",effName)  , nBin, xbins); histo_WS_LepEffUp  ->Sumw2();
  TH1D* histo_WS_LepEffDown       = new TH1D( Form("histo_WS_%sDown",effName), Form("histo_WS_%sDown",effName), nBin, xbins); histo_WS_LepEffDown->Sumw2();
  TH1D* histo_VVV_LepEffUp        = new TH1D( Form("histo_VVV_%sUp",effName)  , Form("histo_VVV_%sUp",effName)  , nBin, xbins); histo_VVV_LepEffUp  ->Sumw2();
  TH1D* histo_VVV_LepEffDown      = new TH1D( Form("histo_VVV_%sDown",effName), Form("histo_VVV_%sDown",effName), nBin, xbins); histo_VVV_LepEffDown->Sumw2();

  TH1D* histo_WWewk_LepResUp      = new TH1D( Form("histo_WWewk_%sUp",momName)  , Form("histo_WWewk_%sUp",momName)  , nBin, xbins); histo_WWewk_LepResUp  ->Sumw2();
  TH1D* histo_WWewk_LepResDown    = new TH1D( Form("histo_WWewk_%sDown",momName), Form("histo_WWewk_%sDown",momName), nBin, xbins); histo_WWewk_LepResDown->Sumw2();
  TH1D* histo_WWqcd_LepResUp      = new TH1D( Form("histo_WWqcd_%sUp",momName)  , Form("histo_WWqcd_%sUp",momName)  , nBin, xbins); histo_WWqcd_LepResUp  ->Sumw2();
  TH1D* histo_WWqcd_LepResDown    = new TH1D( Form("histo_WWqcd_%sDown",momName), Form("histo_WWqcd_%sDown",momName), nBin, xbins); histo_WWqcd_LepResDown->Sumw2();
  TH1D* histo_WZ_LepResUp         = new TH1D( Form("histo_WZ_%sUp",momName)  , Form("histo_WZ_%sUp",momName)  , nBin, xbins); histo_WZ_LepResUp  ->Sumw2();
  TH1D* histo_WZ_LepResDown       = new TH1D( Form("histo_WZ_%sDown",momName), Form("histo_WZ_%sDown",momName), nBin, xbins); histo_WZ_LepResDown->Sumw2();
  TH1D* histo_WS_LepResUp         = new TH1D( Form("histo_WS_%sUp",momName)  , Form("histo_WS_%sUp",momName)  , nBin, xbins); histo_WS_LepResUp  ->Sumw2();
  TH1D* histo_WS_LepResDown       = new TH1D( Form("histo_WS_%sDown",momName), Form("histo_WS_%sDown",momName), nBin, xbins); histo_WS_LepResDown->Sumw2();
  TH1D* histo_VVV_LepResUp        = new TH1D( Form("histo_VVV_%sUp",momName)  , Form("histo_VVV_%sUp",momName)  , nBin, xbins); histo_VVV_LepResUp  ->Sumw2();
  TH1D* histo_VVV_LepResDown      = new TH1D( Form("histo_VVV_%sDown",momName), Form("histo_VVV_%sDown",momName), nBin, xbins); histo_VVV_LepResDown->Sumw2();

  TH1D* histo_WWewk_METResUp      = new TH1D( Form("histo_WWewk_%sUp","CMS_scale_met")  , Form("histo_WWewk_%sUp","CMS_scale_met")  , nBin, xbins); histo_WWewk_METResUp  ->Sumw2();
  TH1D* histo_WWewk_METResDown    = new TH1D( Form("histo_WWewk_%sDown","CMS_scale_met"), Form("histo_WWewk_%sDown","CMS_scale_met"), nBin, xbins); histo_WWewk_METResDown->Sumw2();
  TH1D* histo_WWqcd_METResUp      = new TH1D( Form("histo_WWqcd_%sUp","CMS_scale_met")  , Form("histo_WWqcd_%sUp","CMS_scale_met")  , nBin, xbins); histo_WWqcd_METResUp  ->Sumw2();
  TH1D* histo_WWqcd_METResDown    = new TH1D( Form("histo_WWqcd_%sDown","CMS_scale_met"), Form("histo_WWqcd_%sDown","CMS_scale_met"), nBin, xbins); histo_WWqcd_METResDown->Sumw2();
  TH1D* histo_WZ_METResUp         = new TH1D( Form("histo_WZ_%sUp","CMS_scale_met")  , Form("histo_WZ_%sUp","CMS_scale_met")  , nBin, xbins); histo_WZ_METResUp  ->Sumw2();
  TH1D* histo_WZ_METResDown       = new TH1D( Form("histo_WZ_%sDown","CMS_scale_met"), Form("histo_WZ_%sDown","CMS_scale_met"), nBin, xbins); histo_WZ_METResDown->Sumw2();
  TH1D* histo_WS_METResUp         = new TH1D( Form("histo_WS_%sUp","CMS_scale_met")  , Form("histo_WS_%sUp","CMS_scale_met")  , nBin, xbins); histo_WS_METResUp  ->Sumw2();
  TH1D* histo_WS_METResDown       = new TH1D( Form("histo_WS_%sDown","CMS_scale_met"), Form("histo_WS_%sDown","CMS_scale_met"), nBin, xbins); histo_WS_METResDown->Sumw2();
  TH1D* histo_VVV_METResUp        = new TH1D( Form("histo_VVV_%sUp","CMS_scale_met")  , Form("histo_VVV_%sUp","CMS_scale_met")  , nBin, xbins); histo_VVV_METResUp  ->Sumw2();
  TH1D* histo_VVV_METResDown      = new TH1D( Form("histo_VVV_%sDown","CMS_scale_met"), Form("histo_VVV_%sDown","CMS_scale_met"), nBin, xbins); histo_VVV_METResDown->Sumw2();

  TH1D* histo_WWewk_JESUp         = new TH1D( Form("histo_WWewk_%sUp","CMS_scale_j")  , Form("histo_WWewk_%sUp","CMS_scale_j")  , nBin, xbins); histo_WWewk_JESUp  ->Sumw2();
  TH1D* histo_WWewk_JESDown       = new TH1D( Form("histo_WWewk_%sDown","CMS_scale_j"), Form("histo_WWewk_%sDown","CMS_scale_j"), nBin, xbins); histo_WWewk_JESDown->Sumw2();
  TH1D* histo_WWqcd_JESUp         = new TH1D( Form("histo_WWqcd_%sUp","CMS_scale_j")  , Form("histo_WWqcd_%sUp","CMS_scale_j")  , nBin, xbins); histo_WWqcd_JESUp  ->Sumw2();
  TH1D* histo_WWqcd_JESDown       = new TH1D( Form("histo_WWqcd_%sDown","CMS_scale_j"), Form("histo_WWqcd_%sDown","CMS_scale_j"), nBin, xbins); histo_WWqcd_JESDown->Sumw2();
  TH1D* histo_WZ_JESUp            = new TH1D( Form("histo_WZ_%sUp","CMS_scale_j")  , Form("histo_WZ_%sUp","CMS_scale_j")  , nBin, xbins); histo_WZ_JESUp  ->Sumw2();
  TH1D* histo_WZ_JESDown          = new TH1D( Form("histo_WZ_%sDown","CMS_scale_j"), Form("histo_WZ_%sDown","CMS_scale_j"), nBin, xbins); histo_WZ_JESDown->Sumw2();
  TH1D* histo_WS_JESUp            = new TH1D( Form("histo_WS_%sUp","CMS_scale_j")  , Form("histo_WS_%sUp","CMS_scale_j")  , nBin, xbins); histo_WS_JESUp  ->Sumw2();
  TH1D* histo_WS_JESDown          = new TH1D( Form("histo_WS_%sDown","CMS_scale_j"), Form("histo_WS_%sDown","CMS_scale_j"), nBin, xbins); histo_WS_JESDown->Sumw2();
  TH1D* histo_VVV_JESUp           = new TH1D( Form("histo_VVV_%sUp","CMS_scale_j")  , Form("histo_VVV_%sUp","CMS_scale_j")  , nBin, xbins); histo_VVV_JESUp  ->Sumw2();
  TH1D* histo_VVV_JESDown         = new TH1D( Form("histo_VVV_%sDown","CMS_scale_j"), Form("histo_VVV_%sDown","CMS_scale_j"), nBin, xbins); histo_VVV_JESDown->Sumw2();

  TH1D* histo_WZ_CMS_WZNLOUp      = new TH1D( Form("histo_WZ_CMS_wwss_WZNLOUp"),   Form("histo_WZ_CMS_wwss_WZNLOUp"),	nBin, xbins); histo_WZ_CMS_WZNLOUp  ->Sumw2();
  TH1D* histo_WZ_CMS_WZNLODown    = new TH1D( Form("histo_WZ_CMS_wwss_WZNLODown"), Form("histo_WZ_CMS_wwss_WZNLODown"), nBin, xbins); histo_WZ_CMS_WZNLODown->Sumw2();

  TH1D* histo_WS_WSUp             = new TH1D( Form("histo_WS_CMS_wwss_MVAWSUp"),   Form("histo_WS_CMS_wwss_MVAWSUp"),   nBin, xbins); histo_WS_WSUp  ->Sumw2();
  TH1D* histo_WS_WSDown           = new TH1D( Form("histo_WS_CMS_wwss_MVAWSDown"), Form("histo_WS_CMS_wwss_MVAWSDown"), nBin, xbins); histo_WS_WSDown->Sumw2();

  TH1D* histo_Wjets_WUp           = new TH1D( Form("histo_Wjets_CMS_wwss_MVAWUp"),   Form("histo_Wjets_CMS_wwss_MVAWUp"),   nBin, xbins); histo_Wjets_WUp  ->Sumw2();
  TH1D* histo_Wjets_WDown         = new TH1D( Form("histo_Wjets_CMS_wwss_MVAWDown"), Form("histo_Wjets_CMS_wwss_MVAWDown"), nBin, xbins); histo_Wjets_WDown->Sumw2();

  double nSelectedData[nSelTypes*2];
  double bgdDecay[nSelTypes*2][45],weiDecay[nSelTypes*2][45];
  double bgdDecaySyst[nSelTypesSyst*2][45],weiDecaySyst[nSelTypesSyst*2][45];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    nSelectedData[i] = 0.0; 
    for(int j=0; j<45; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    for(int j=0; j<45; j++) {
      bgdDecaySyst[i][j] = 0.0; weiDecaySyst[i][j] = 0.0; 
    }
  }

  unsigned int patternTopVeto = SmurfTree::TopVeto;

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int evt=0; evt<nBgd; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nBgd);
    bgdEvent.tree_->GetEntry(evt);

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
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(bgdEvent.processId_==121 ||
            bgdEvent.processId_==122)   fDecay = 41;
    else if(bgdEvent.processId_==24)    fDecay = 42;
    else if(bgdEvent.processId_==26)    fDecay = 43;
    else if(bgdEvent.processId_==10001) fDecay = 44;
    else if(bgdEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;}

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

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false, false, false},
                                   {false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
       (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = bgdEvent.met_; double theMETPHI = bgdEvent.metPhi_; 
    
    int lType = 1;
    if     (bgdEvent.lq1_ * bgdEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
       ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(bgdEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(bgdEvent.type_ == SmurfTree::mm) metMin = 40.0;

    // jet1McId and jet2McId actually provide lepton rejection information
    // hasZCand = reject events with |mll-mZ|<15
    // trackSel[0] == reject events with |mll-mZ|<10 for pf candidates with |eta|<3.0
    // trackSel[1] == reject events with |mll-mZ|<10 for pf candidates with |eta|<4.7
    // trackSel[2] == reject events with isolated reconstructed leptons with pt>10 and iso/pt<0.1
    // trackSel[3] == reject events with isolated tracks with pt>10 and iso/pt<0.1 (not used by default)
    int newId=int(bgdEvent.jet1McId_);
    //double wzId=bgdEvent.jet1McId_%10;
    //int tauId=int((bgdEvent.jet1McId_%100-bgdEvent.jet1McId_%10)/10)+int((bgdEvent.jet2McId_%100-bgdEvent.jet2McId_%10)/10)+int((bgdEvent.jet3McId_%100-bgdEvent.jet3McId_%10)/10)+int((bgdEvent.jet4McId_%100-bgdEvent.jet4McId_%10)/10);
    int qDisAgree=int((newId%1000-newId%100)/100);
    int hasZCand=int(newId/1000);
    int trackSel[4] = {int((bgdEvent.jet2McId_%1000-bgdEvent.jet2McId_%100)/10),int((bgdEvent.jet2McId_%10000-bgdEvent.jet2McId_%1000)/10),int((bgdEvent.jet2McId_%100000-bgdEvent.jet2McId_%10000)/10),int(bgdEvent.jet2McId_/100000)};

    bool passNjets         = bgdEvent.njets_ >= 2;
    bool passMET           = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > metMin;
    bool preselCuts        = bgdEvent.lep1_.Pt() > 20. && bgdEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && hasZCand == 0 && trackSel[0]+trackSel[1]+trackSel[2] == 0;
    bool passMass          = bgdEvent.dilep_.M() > 50.0 && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 && centrality == 1;
    // quick way to accept W+W+ (1) or W-W- (2)
    bool passNsignSel = true;
    if(signSel == 1 && bgdEvent.lq1_+bgdEvent.lq2_ == -2) passNsignSel = false;
    if(signSel == 2 && bgdEvent.lq1_+bgdEvent.lq2_ ==  2) passNsignSel = false;

    if(lType == 0) passMass = passMass && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ != SmurfTree::mm);

    // 0      1      2       3     4   5      6        7           8  9            10            11     12  13    14
    // lep1pt,lep2pt,dilmass,dilpt,met,metPhi,trackMet,trackMetPhi,mt,dPhiDiLepMET,dPhiMETTrkMET,pTFrac,mtZ,mlljj,mjj;
    double outputVarLepP[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 0, outputVarLepP);
    double outputVarLepM[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 1, outputVarLepM);
    double outputVarMET[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 2, outputVarMET);
    double outputVar[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 3, outputVar);
    double outputVarJESP[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 4, outputVarJESP);
    double outputVarJESM[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 5, outputVarJESM);
    double MVAVar[6] = {outputVar[13],outputVarJESP[13],outputVarJESM[13],outputVarLepP[13],outputVarLepM[13],outputVarMET[13]};
    if(thePlot == 0) {MVAVar[0]=outputVar[14];MVAVar[1]=outputVarJESP[14];MVAVar[2]=outputVarJESM[14];MVAVar[3]=outputVarLepP[14];MVAVar[4]=outputVarLepM[14];MVAVar[5]=outputVarMET[14];}
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBin]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff  = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 0)*
    		 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 1);
      addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1);
    } else {addLepEff = 1.0;}

    double NjetSyst[2] = {0., 0.};
    if(bgdEvent.jet1_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet2_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet3_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet4_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet1_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet2_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet3_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet4_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;

    bool passLSel = false;
    if     (lSel == 0 && bgdEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && bgdEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && bgdEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && bgdEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (bgdEvent.type_ == SmurfTree::me || bgdEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passNsignSel && qDisAgree == 0 && NjetSyst[0] >= 2	  && TMath::Min(outputVarJESP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_)       > metMin && outputVar[0]	  > 20.0 && outputVar[1]     > 20.0 && passBtagVeto && pass3rLVeto && outputVar[2]     > 50.0 && passVBFSel) passSystCuts[lType][JESUP] = true;
    if(passNsignSel && qDisAgree == 0 && NjetSyst[1] >= 2	  && TMath::Min(outputVarJESM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_)       > metMin && outputVar[0]	  > 20.0 && outputVar[1]     > 20.0 && passBtagVeto && pass3rLVeto && outputVar[2]     > 50.0 && passVBFSel) passSystCuts[lType][JESDOWN] = true;
    if(passNsignSel && qDisAgree == 0 && bgdEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarLepP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarLepP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVarLepP[0] > 20.0 && outputVarLepP[1] > 20.0 && passBtagVeto && pass3rLVeto && outputVarLepP[2] > 50.0 && passVBFSel) passSystCuts[lType][LEPP] = true;
    if(passNsignSel && qDisAgree == 0 && bgdEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarLepM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarLepM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVarLepM[0] > 20.0 && outputVarLepM[1] > 20.0 && passBtagVeto && pass3rLVeto && outputVarLepM[2] > 50.0 && passVBFSel) passSystCuts[lType][LEPM] = true;
    if(passNsignSel && qDisAgree == 0 && bgdEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarMET[4] /bgdEvent.met_*bgdEvent.pmet_,outputVarMET[6] /bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVarMET[0]  > 20.0 && outputVarMET[1]  > 20.0 && passBtagVeto && pass3rLVeto && outputVarMET[2]  > 50.0 && passVBFSel) passSystCuts[lType][MET] = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && bgdEvent.dilep_.M() > 15.0 && qDisAgree == 0) {
       
       if(passNsignSel &&  passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;
       if(passNsignSel && !passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][BTAGSEL] = true;
       if(passNsignSel &&  passBtagVeto && passVBFSel == true                      && !pass3rLVeto && bgdEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL] = true;

      if(isRealLepton == false &&
         (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
          bgdEvent.dstype_ == SmurfTree::qqww	|| bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
          bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) 
        {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false;}
    }

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
	  if(fCheckProblem == true && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)+add)/add>0.0001)
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

        if(fCheckProblem == true && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) && add != 0 && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
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
      scaleFactor_WS(bgdEvent.lep1_,bgdEvent.lq1_,bgdEvent.lid1_,bgdEvent.lep1McId_,weightWS);
      scaleFactor_WS(bgdEvent.lep2_,bgdEvent.lq2_,bgdEvent.lid2_,bgdEvent.lep2McId_,weightWS);
      if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
        scaleFactor_WS(bgdEvent.lep3_,bgdEvent.lq3_,bgdEvent.lid3_,bgdEvent.lep3McId_,weightWS);
      }
      theWeight = weightWS[0];
      
      if(bgdEvent.dstype_ == SmurfTree::ttbar) theWeight = theWeight * 1.07841;
      if(bgdEvent.dstype_ == SmurfTree::tw)    theWeight = theWeight * 1.07841;

      if(passCuts[1][WWSEL]){ // begin making plots
	double myVar = -1.0;
	if     (thePlot == 0) myVar = TMath::Max(TMath::Min((bgdEvent.jet1_+bgdEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot == 1) myVar = TMath::Max(TMath::Min((bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot == 2) myVar = bgdEvent.lep1_.Pt();
	else if(thePlot == 3) myVar = bgdEvent.lep2_.Pt();
	else if(thePlot == 4) myVar = bgdEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = bgdEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = bgdEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = bgdEvent.mt_;
	else if(thePlot == 8) myVar = bgdEvent.dilep_.M();
	else if(thePlot == 9) myVar = bgdEvent.dilep_.Pt();
	else if(thePlot ==10) myVar = bgdEvent.njets_;
	else if(thePlot ==11) myVar = bgdEvent.nvtx_;
	else if(thePlot ==12) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==13) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==14) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==15) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
	else if(thePlot ==16) myVar = bgdEvent.type_;
	else if(thePlot ==17) myVar = bgdEvent.dR_;
	else assert(0);
      	if     (fDecay == 31){
	  if(use_anom_sample){
	    if(bgdEvent.scale1fb_ > 0){
	      theWeight=theWeight*bgdEvent.lheWeights_[which_lhe_weight]/bgdEvent.lheWeights_[0];
	    }
	  }
      	  histo0->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 21){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 27){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 29){
      	  histo2->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 1 || fDecay == 23 || fDecay == 3){
      	  histo3->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 30 || fDecay == 28 ||
                fDecay ==  5 || fDecay == 13 || fDecay == 20 || 
		fDecay == 10 || fDecay ==  9 || fDecay == 19){
      	  histo4->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 41 || fDecay == 42 || fDecay == 43){
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
      for(unsigned int i=0; i<5; i++) {
        for(int j=0; j<2; j++){
          if(passSystCuts[j][i]) {
            bgdDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight;
            weiDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][WWSEL]) {
        bgdDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffUp  /addLepEff;
        weiDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        bgdDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffDown/addLepEff;
        weiDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }

      if     (fDecay == 21){
        if(passCuts[1][WWSEL])  	     histo_VVV           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_VVV_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_VVV_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_VVV_JESUp	 ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_VVV_JESDown	 ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_VVV_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_VVV_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_VVV_METResUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 31){
        if(passCuts[1][WWSEL])  	     histo_WWewk           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WWewk_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WWewk_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WWewk_JESUp	   ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WWewk_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WWewk_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WWewk_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WWewk_METResUp  ->Fill(MVAVar[5], theWeight);;

        if(passCuts[1][WWSEL] && doAQGCsAna == true){
	  //the qcd WW that we subtract from the signal is not reweighted
	  if(bgdEvent.scale1fb_ > 0){
	    assert(bgdEvent.lheWeights_.size() >= grid_points.size());
	    assert(grid_points.size() == histo_grid.size());
	    assert(lhe_weight_index.size() == histo_grid.size());
	  }
	  for(unsigned int a = 0; a < grid_points.size(); a++){
	    if(bgdEvent.scale1fb_ > 0)
	      histo_WWewk_anom[a]->Fill(MVAVar[0],theWeight*bgdEvent.lheWeights_[lhe_weight_index[a]]/bgdEvent.lheWeights_[0]);
	    else
	      histo_WWewk_anom[a]->Fill(MVAVar[0],theWeight);
	  }
	}
      }
      else if(fDecay == 27){
        if(passCuts[1][WWSEL])  	     histo_WZ           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WZ_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WZ_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WZ_JESUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WZ_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WZ_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WZ_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WZ_METResUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 29){
        if(passCuts[1][WWSEL])  	     histo_WWqcd           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WWqcd_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WWqcd_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WWqcd_JESUp	   ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WWqcd_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WWqcd_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WWqcd_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WWqcd_METResUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 30 || fDecay == 28 ||
              fDecay ==  5 || fDecay == 13 || fDecay == 20 || 
	      fDecay == 10 || fDecay ==  9 || fDecay == 19){
        if(passCuts[1][WWSEL])  	     histo_WS           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WS_LepEffUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WS_LepEffDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WS_JESUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WS_JESDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WS_LepResUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WS_LepResDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WS_METResUp  ->Fill(MVAVar[5], theWeight);;
        if(passCuts[1][WWSEL])  	     histo_WS_WSUp      ->Fill(MVAVar[0], weightWS[1]);
      }
      else if(fDecay == 1 || fDecay == 23 || fDecay == 3){
        double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        												     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        												     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        												     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        												     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        if(passCuts[1][WWSEL]) 	     histo_Wjets		       ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL]) 	     histo_Wjets_WUp    ->Fill(MVAVar[0], theWeight*addFRS/addFR);
      }
      else {
        printf("Forbidden fDecay: %d\n",fDecay);
	assert(0);
      }
    } // if passCuts
  } // end background loop

  if(systInputFile != ""){
  int nSyst=systEvent.tree_->GetEntries();
  for (int evt=0; evt<nSyst; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nSyst);
    systEvent.tree_->GetEntry(evt);

    if(systEvent.dstype_ == SmurfTree::data &&
      (systEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ <  minRun) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (systEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(systEvent.dstype_ == SmurfTree::wjets 	   ) fDecay =  3;
    else if(systEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay =  5;
    else if(systEvent.dstype_ == SmurfTree::dyee  	   ) fDecay =  9;
    else if(systEvent.dstype_ == SmurfTree::dymm  	   ) fDecay =  9;
    else if(systEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::dyttDataDriven ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(systEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(systEvent.dstype_ == SmurfTree::wgstar         ) fDecay = 20;
    else if(systEvent.dstype_ == SmurfTree::www            ) fDecay = 21;
    else if(systEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(systEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(systEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(systEvent.dstype_ == SmurfTree::other          ) fDecay = 40;
    else if(systEvent.processId_==121 ||
            systEvent.processId_==122)   fDecay = 41;
    else if(systEvent.processId_==24)    fDecay = 42;
    else if(systEvent.processId_==26)    fDecay = 43;
    else if(systEvent.processId_==10001) fDecay = 44;
    else if(systEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << systEvent.dstype_ << std::endl;}

    bool passCuts[3][nSelTypes] = {{false, false, false},
                                   {false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13) &&
       (TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = systEvent.met_; double theMETPHI = systEvent.metPhi_; 
    
    int lType = 1;
    if     (systEvent.lq1_ * systEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((systEvent.jet1_.Eta()-systEvent.lep1_.Eta() > 0 && systEvent.jet2_.Eta()-systEvent.lep1_.Eta() < 0) ||
        (systEvent.jet2_.Eta()-systEvent.lep1_.Eta() > 0 && systEvent.jet1_.Eta()-systEvent.lep1_.Eta() < 0)) &&
       ((systEvent.jet1_.Eta()-systEvent.lep2_.Eta() > 0 && systEvent.jet2_.Eta()-systEvent.lep2_.Eta() < 0) ||
        (systEvent.jet2_.Eta()-systEvent.lep2_.Eta() > 0 && systEvent.jet1_.Eta()-systEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(systEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(systEvent.type_ == SmurfTree::mm) metMin = 40.0;

    int newId=int(systEvent.jet1McId_);
    //double wzId=systEvent.jet1McId_%10;
    //int tauId=int((systEvent.jet1McId_%100-systEvent.jet1McId_%10)/10)+int((systEvent.jet2McId_%100-systEvent.jet2McId_%10)/10)+int((systEvent.jet3McId_%100-systEvent.jet3McId_%10)/10)+int((systEvent.jet4McId_%100-systEvent.jet4McId_%10)/10);
    int qDisAgree=int((newId%1000-newId%100)/100);
    int hasZCand=int(newId/1000);
    int trackSel[4] = {int((systEvent.jet2McId_%1000-systEvent.jet2McId_%100)/10),int((systEvent.jet2McId_%10000-systEvent.jet2McId_%1000)/10),int((systEvent.jet2McId_%100000-systEvent.jet2McId_%10000)/10),int(systEvent.jet2McId_/100000)};

    bool passNjets         = systEvent.njets_ >= 2;
    bool passMET           = TMath::Min(systEvent.pmet_,systEvent.pTrackMet_) > metMin;
    bool preselCuts        = systEvent.lep1_.Pt() > 20. && systEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (systEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && hasZCand == 0 && trackSel[0]+trackSel[1]+trackSel[2] == 0;
    bool passMass          = systEvent.dilep_.M() > 50.0 && (TMath::Abs(systEvent.dilep_.M()-91.1876) > 15 || systEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (systEvent.jet1_+systEvent.jet2_).M() > 500 && TMath::Abs(systEvent.jet1_.Eta()-systEvent.jet2_.Eta()) > 3.5 && centrality == 1;
    // quick way to accept W+W+ (1) or W-W- (2)
    bool passNsignSel = true;
    if(signSel == 1 && systEvent.lq1_+systEvent.lq2_ == -2) passNsignSel = false;
    if(signSel == 2 && systEvent.lq1_+systEvent.lq2_ ==  2) passNsignSel = false;

    if(lType == 0) passMass = passMass && (TMath::Abs(systEvent.dilep_.M()-91.1876) > 15 || systEvent.type_ != SmurfTree::mm);

    bool passLSel = false;
    if     (lSel == 0 && systEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && systEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && systEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && systEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                     passLSel = true;
    else if(lSel == 5 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (systEvent.type_ == SmurfTree::me || systEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && systEvent.dilep_.M() > 15.0 && qDisAgree == 0) {
      if(passNsignSel &&  passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;

      if(isRealLepton == false &&
         (systEvent.dstype_ == SmurfTree::ttbar  || systEvent.dstype_ == SmurfTree::tw   || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dymm ||
          systEvent.dstype_ == SmurfTree::qqww   || systEvent.dstype_ == SmurfTree::ggww || systEvent.dstype_ == SmurfTree::wz   || systEvent.dstype_ == SmurfTree::zz   ||
          systEvent.dstype_ == SmurfTree::wgstar || systEvent.dstype_ == SmurfTree::dytt || systEvent.dstype_ == SmurfTree::www)) 
        {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false;}

    }

    if(passCuts[lType][WWSEL]){
      double theWeight = 0.0;
      double add       = 1.0;
      int nFake = 0;
      if(((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(nFake < 0) assert(0);
 
      if(nFake > 1){
	add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											  (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	fDecay = 22;

	theWeight	       = -1.0*add;
      }
      else if(nFake == 1){
        if(systEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          if(fCheckProblem == true && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  // new category, W+jetsM
	  if((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0;
	}
	else if(isRealLepton == true || systEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  add = add*nPUScaleFactor2012(fhDPU ,systEvent.npu_);
          add = add*leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	  add = add*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
          if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add = add*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt(), 
								   fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						   TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
          add = add*trigEff;
	  if(fCheckProblem == true && (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMBSyst: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  fDecay                 = 1;

	  if((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * systEvent.scale1fb_*lumi*add;
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(systEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
	        					        fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
							        TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        double sf_eff = 1.0;
	sf_eff = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_)*
        	 leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);

        theWeight = ZttScaleFactor(period,systEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(systEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,systEvent.npu_);
        double add2 = 1.0;
	add2 = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	add2 = add2*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
        if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add2 = add2*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
								 fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						 TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        add = add1*add2*trigEff;

        if(fCheckProblem == true && add != 0 && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	printf("PROBLEMCSy(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",systEvent.event_,add1,add2,trigEff,add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);

	if(systEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(systEvent.type_,theMET);

        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (systEvent.dstype_ == SmurfTree::dymm || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dytt) &&
          (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me)) add = 0.0;

        //----------------------------------------------------------------------------      
        // Apply weighting factor to wgamma (gamma->electron ratio)
        //----------------------------------------------------------------------------
        if(systEvent.dstype_ == SmurfTree::wgamma) {
          if(!(TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13)) add = add * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(systEvent.lep1_.Eta()));
          if(!(TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) add = add * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(systEvent.lep2_.Eta()));      
        }

	theWeight              = systEvent.scale1fb_*lumi*add;
      }

      double outputVar[15];
      makeSystematicEffects(systEvent.lid1_, systEvent.lid2_, systEvent.lep1_, systEvent.lep2_, systEvent.dilep_, 
                            systEvent.mt_, theMET, theMETPHI, 
                            systEvent.trackMet_, systEvent.trackMetPhi_, 
			    systEvent.njets_, systEvent.jet1_, systEvent.jet2_,
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[13],0,0,0,0,0};
      if(thePlot == 0) {MVAVar[0]=outputVar[14];}
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBin]-0.001);
      if(passCuts[1][WWSEL]){
	if     (fDecay == 27){
	  histo_WZ_CMS_WZNLOUp->Fill(MVAVar[0], theWeight);
        }
      }

    } // if passCuts
  } // end syst loop
  } // if want to use it at all
  
  if(run_over_data){
  int nData=dataEvent.tree_->GetEntries();
  for (int evt=0; evt<nData; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nData);
    dataEvent.tree_->GetEntry(evt);

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;

    bool passCuts[3][nSelTypes] = {{false, false, false},
                                   {false, false, false}};

    double theMET = dataEvent.met_; double theMETPHI = dataEvent.metPhi_; 

    int lType = 1;
    if     (dataEvent.lq1_ * dataEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
       ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(dataEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(dataEvent.type_ == SmurfTree::mm) metMin = 40.0;

    int newId=int(dataEvent.jet1McId_);
    //double wzId=dataEvent.jet1McId_%10;
    //int tauId=int((dataEvent.jet1McId_%100-dataEvent.jet1McId_%10)/10)+int((dataEvent.jet2McId_%100-dataEvent.jet2McId_%10)/10)+int((dataEvent.jet3McId_%100-dataEvent.jet3McId_%10)/10)+int((dataEvent.jet4McId_%100-dataEvent.jet4McId_%10)/10);
    int qDisAgree=int((newId%1000-newId%100)/100);
    int hasZCand=int(newId/1000);
    int trackSel[4] = {int((dataEvent.jet2McId_%1000-dataEvent.jet2McId_%100)/10),int((dataEvent.jet2McId_%10000-dataEvent.jet2McId_%1000)/10),int((dataEvent.jet2McId_%100000-dataEvent.jet2McId_%10000)/10),int(dataEvent.jet2McId_/100000)};

    bool passNjets         = dataEvent.njets_ >= 2;
    bool passMET           = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_) > metMin;
    bool preselCuts        = dataEvent.lep1_.Pt() > 20. && dataEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && hasZCand == 0 && trackSel[0]+trackSel[1]+trackSel[2] == 0;
    bool passMass          = dataEvent.dilep_.M() > 50.0 && (TMath::Abs(dataEvent.dilep_.M()-91.1876) > 15 || dataEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (dataEvent.jet1_+dataEvent.jet2_).M() > 500 && TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 3.5 && centrality == 1;
    // quick way to accept W+W+ (1) or W-W- (2)
    bool passNsignSel = true;
    if(signSel == 1 && dataEvent.lq1_+dataEvent.lq2_ == -2) passNsignSel = false;
    if(signSel == 2 && dataEvent.lq1_+dataEvent.lq2_ ==  2) passNsignSel = false;

    if(lType == 0) passMass = passMass && (TMath::Abs(dataEvent.dilep_.M()-91.1876) > 15 || dataEvent.type_ != SmurfTree::mm);

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
       
       if(passNsignSel &&  passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;
       if(passNsignSel && !passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][BTAGSEL] = true;
       if(passNsignSel &&  passBtagVeto && passVBFSel == true                      && !pass3rLVeto && dataEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL] = true;

    }

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && dataEvent.dilep_.M() > 15.0) {

      if(passCuts[1][WWSEL]){ // begin making plots
	double myVar = -1.0;
	if     (thePlot == 0) myVar = TMath::Max(TMath::Min((dataEvent.jet1_+dataEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot == 1) myVar = TMath::Max(TMath::Min((dataEvent.lep1_+dataEvent.lep2_+dataEvent.jet1_+dataEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot == 3) myVar = dataEvent.lep2_.Pt();
	else if(thePlot == 4) myVar = dataEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = dataEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = dataEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = dataEvent.mt_;
	else if(thePlot == 8) myVar = dataEvent.dilep_.M();
	else if(thePlot == 9) myVar = dataEvent.dilep_.Pt();
	else if(thePlot ==10) myVar = dataEvent.njets_;
	else if(thePlot ==11) myVar = dataEvent.nvtx_;
	else if(thePlot ==12) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==13) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==14) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==15) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
	else if(thePlot ==16) myVar = dataEvent.type_;
	else if(thePlot ==17) myVar = dataEvent.dR_;
	else assert(0);
      	histo5->Fill(myVar,1.0);
      } // end making plots

      double outputVar[15];
      makeSystematicEffects(dataEvent.lid1_, dataEvent.lid2_, dataEvent.lep1_, dataEvent.lep2_, dataEvent.dilep_, 
                            dataEvent.mt_, theMET, theMETPHI, 
                            dataEvent.trackMet_, dataEvent.trackMetPhi_, 
			    dataEvent.njets_, dataEvent.jet1_, dataEvent.jet2_, 
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[13],0,0,0,0,0};
      if(thePlot == 0) {MVAVar[0]=outputVar[14];}
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBin]-0.001);
      if(passCuts[1][WWSEL]){
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

  TFile *th2d_outfile; 
  if(doAQGCsAna == true){
    th2d_outfile = new TFile("aQGC_grids.root","recreate");
  }

  outFilePlotsNote->cd();
    double nOldH[5] = {histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights()};
    for(int i=1; i<=histo0->GetNbinsX(); i++){
      if(histo0->GetBinContent(i) < 0) {histo0->SetBinContent(i,0.000001);histo0->SetBinError(i,0.000001);}
      if(histo1->GetBinContent(i) < 0) {histo1->SetBinContent(i,0.000001);histo1->SetBinError(i,0.000001);}
      if(histo2->GetBinContent(i) < 0) {histo2->SetBinContent(i,0.000001);histo2->SetBinError(i,0.000001);}
      if(histo3->GetBinContent(i) < 0) {histo3->SetBinContent(i,0.000001);histo3->SetBinError(i,0.000001);}
      if(histo4->GetBinContent(i) < 0) {histo4->SetBinContent(i,0.000001);histo4->SetBinError(i,0.000001);}
    }
    if(nOldH[0] > 0) histo0->Scale(nOldH[0]/histo0->GetSumOfWeights());
    if(nOldH[1] > 0) histo1->Scale(nOldH[1]/histo1->GetSumOfWeights());
    if(nOldH[2] > 0) histo2->Scale(nOldH[2]/histo2->GetSumOfWeights());
    if(nOldH[3] > 0) histo3->Scale(nOldH[3]/histo3->GetSumOfWeights());
    if(nOldH[4] > 0) histo4->Scale(nOldH[4]/histo4->GetSumOfWeights());

    printf("histo -> d: %8.2f b: %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f\n",histo5->GetSumOfWeights(),
    histo0->GetSumOfWeights()+histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights(),
    histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights());

    histo0->Write();
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();

  outFilePlotsNote->Close();
  
  const unsigned int nBkg = 6;
  double nTot[nSelTypes*2]; double nETot[nSelTypes*2];
  double bgdCombined[nSelTypes*2][nBkg],bgdCombinedE[nSelTypes*2][nBkg];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) {bgdCombined[i][j] = 0.0; bgdCombinedE[i][j] = 0.0;}
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("selection: %s\n",selTypeName[i].Data());
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("data(%2d): %f\n",i,nSelectedData[i]);
    nTot[i] = 0.0; nETot[i] = 0.0;
    for(int j=0; j<45; j++){
      // WWqcd treatment
      if(j == 29 && bgdDecay[i][j] < 0) {printf("negative(29,%d) = %f +/- %f\n",i,bgdDecay[i][j],sqrt(weiDecay[i][j]));bgdDecay[i][j] = 0; weiDecay[i][j] = 0;}

      if(showSignalOnly == false || i%nSelTypes == WWSEL) if(bgdDecay[i][j] != 0) printf("bdg(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecay[i][j],sqrt(weiDecay[i][j]));

      nTot[i]  += bgdDecay[i][j];
      nETot[i] += weiDecay[i][j];

      if     (j == 31)                         {bgdCombined[i][0] += bgdDecay[i][j]; bgdCombinedE[i][0] += weiDecay[i][j];}
      else if(j == 29)                         {bgdCombined[i][1] += bgdDecay[i][j]; bgdCombinedE[i][1] += weiDecay[i][j];}
      else if(j == 27)			       {bgdCombined[i][2] += bgdDecay[i][j]; bgdCombinedE[i][2] += weiDecay[i][j];}
      else if(j == 30 || j == 28 ||
              j ==  5 || j == 13 || j == 20 || 
	      j == 10 || j ==  9 || j == 19)   {bgdCombined[i][3] += bgdDecay[i][j]; bgdCombinedE[i][3] += weiDecay[i][j];}
      else if(j == 21)			       {bgdCombined[i][4] += bgdDecay[i][j]; bgdCombinedE[i][4] += weiDecay[i][j];}
      else if(j == 1 || j == 23)               {bgdCombined[i][5] += bgdDecay[i][j]; bgdCombinedE[i][5] += weiDecay[i][j];}
    }
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTot[i],sqrt(nETot[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]-nTot[i])/nTot[i] > 0.00001) 
                    {printf("%f\n",bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]);assert(0);}
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("------\n");
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(WWe) = %11.3f +/- %8.3f\n",bgdCombined[i][0],sqrt(bgdCombinedE[i][0]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(WWq) = %11.3f +/- %8.3f\n",bgdCombined[i][1],sqrt(bgdCombinedE[i][1]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWZ) = %11.3f +/- %8.3f\n",bgdCombined[i][2],sqrt(bgdCombinedE[i][2]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWS) = %11.3f +/- %8.3f\n",bgdCombined[i][3],sqrt(bgdCombinedE[i][3]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(VVV) = %11.3f +/- %8.3f\n",bgdCombined[i][4],sqrt(bgdCombinedE[i][4]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWj) = %11.3f +/- %8.3f\n",bgdCombined[i][5],sqrt(bgdCombinedE[i][5]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("*******************************\n");
  }

  if(showSignalOnly == false) printf("+++++++++++++++++++++++++++++++\n");
  double nTotSyst[nSelTypesSyst*2]; double nETotSyst[nSelTypesSyst*2];
  double bgdCombinedSyst[nSelTypesSyst*2][nBkg],bgdCombinedESyst[nSelTypesSyst*2][nBkg];
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    if(showSignalOnly == false) printf("selectionSyst: %s\n",selTypeNameSyst[i].Data());
    for(unsigned int j=0; j<nBkg; j++) {bgdCombinedSyst[i][j] = 0.0; bgdCombinedESyst[i][j] = 0.0;}
    nTotSyst[i] = 0.0; nETotSyst[i] = 0.0;
    for(int j=0; j<45; j++){
      if(showSignalOnly == false) if(bgdDecaySyst[i][j] != 0) printf("bdgSyst(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecaySyst[i][j],sqrt(weiDecaySyst[i][j]));
      nTotSyst[i]  += bgdDecaySyst[i][j];
      nETotSyst[i] += weiDecaySyst[i][j];

     if     (j == 31)			      {bgdCombinedSyst[i][0] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][0] += weiDecaySyst[i][j];}
     else if(j == 29)			      {bgdCombinedSyst[i][1] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][1] += weiDecaySyst[i][j];}
     else if(j == 27)			      {bgdCombinedSyst[i][2] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][2] += weiDecaySyst[i][j];}
     else if(j == 30 || j == 28 ||		      
     	     j ==  5 || j == 13 || j == 20 || 
     	     j == 10 || j ==  9 || j == 19)   {bgdCombinedSyst[i][3] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][3] += weiDecaySyst[i][j];}
     else if(j == 21)			      {bgdCombinedSyst[i][4] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][4] += weiDecaySyst[i][j];}
     else if(j == 1 || j == 23) 	      {bgdCombinedSyst[i][5] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][5] += weiDecaySyst[i][j];}
    }
    if(showSignalOnly == false) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTotSyst[i],sqrt(nETotSyst[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]-nTotSyst[i])/nTotSyst[i] > 0.00001) 
                    {printf("%f\n",bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]);assert(0);}
    if(showSignalOnly == false) printf("------\n");
    if(showSignalOnly == false) printf("bgdSyst(WWe) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][0],sqrt(bgdCombinedESyst[i][0]));
    if(showSignalOnly == false) printf("bgdSyst(WWq) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][1],sqrt(bgdCombinedESyst[i][1]));
    if(showSignalOnly == false) printf("bgdSyst(xWZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][2],sqrt(bgdCombinedESyst[i][2]));
    if(showSignalOnly == false) printf("bgdSyst(xWS) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][3],sqrt(bgdCombinedESyst[i][3]));
    if(showSignalOnly == false) printf("bgdSyst(VVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][4],sqrt(bgdCombinedESyst[i][4]));
    if(showSignalOnly == false) printf("bgdSyst(xWj) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][5],sqrt(bgdCombinedESyst[i][5]));
    if(showSignalOnly == false) printf("*******************************\n");
  }

  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) if(bgdCombined[i][j] == 0) {bgdCombined[i][j] = 0.0000000001; bgdCombinedE[i][j] = 0.0;}
  }

  double QCDscale_WWewk = 1.05;
  double QCDscale_WWqcd = 1.16;

  double systEffect[nSelTypesSyst][nBkg];
  for(unsigned int i=0 ; i<nSelTypesSyst; i++){
    for(unsigned int j=0 ; j<nBkg; j++){
      if(bgdCombinedE[WWSEL+nSelTypes][j] > 0){
        systEffect[i][j] = bgdCombinedSyst[i+nSelTypesSyst][j]/bgdCombined[WWSEL+nSelTypes][j];
        if(systEffect[i][j] < 1) systEffect[i][j] = 1.0/systEffect[i][j];
      } else {systEffect[i][j] = 1.0;}
    }
  }
  if(showSignalOnly == false) printf("Syst(WWe) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][0]-1,systEffect[JESDOWN][0]-1,systEffect[LEPP][0]-1,systEffect[LEPM][0]-1,systEffect[MET][0]-1,systEffect[EFFP][0]-1,systEffect[EFFM][0]-1);
  if(showSignalOnly == false) printf("Syst(WWq) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][1]-1,systEffect[JESDOWN][1]-1,systEffect[LEPP][1]-1,systEffect[LEPM][1]-1,systEffect[MET][1]-1,systEffect[EFFP][1]-1,systEffect[EFFM][1]-1);
  if(showSignalOnly == false) printf("Syst(xWZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][2]-1,systEffect[JESDOWN][2]-1,systEffect[LEPP][2]-1,systEffect[LEPM][2]-1,systEffect[MET][2]-1,systEffect[EFFP][2]-1,systEffect[EFFM][2]-1);
  if(showSignalOnly == false) printf("Syst(xWS) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][3]-1,systEffect[JESDOWN][3]-1,systEffect[LEPP][3]-1,systEffect[LEPM][3]-1,systEffect[MET][3]-1,systEffect[EFFP][3]-1,systEffect[EFFM][3]-1);
  if(showSignalOnly == false) printf("Syst(VVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][4]-1,systEffect[JESDOWN][4]-1,systEffect[LEPP][4]-1,systEffect[LEPM][4]-1,systEffect[MET][4]-1,systEffect[EFFP][4]-1,systEffect[EFFM][4]-1);
  if(showSignalOnly == false) printf("Syst(xWj) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][5]-1,systEffect[JESDOWN][5]-1,systEffect[LEPP][5]-1,systEffect[LEPM][5]-1,systEffect[MET][5]-1,systEffect[EFFP][5]-1,systEffect[EFFM][5]-1);

  double WjetsSyst = 1.0;
  if(bgdCombined[WWSEL+nSelTypes][5] > 0){
    WjetsSyst = bgdCombined[WWSEL+nSelTypes][5]*0.36;
    WjetsSyst = 1.0+WjetsSyst/(bgdCombined[WWSEL+nSelTypes][5]);
  }
  if(showSignalOnly == false) printf("WjetsSyst: %f --> %f\n",bgdCombined[WWSEL+nSelTypes][5],WjetsSyst);
  double pdf_qqbar[3] = {1.073,1.068,1.069};
  double syst_WZ3l = 1.010;

  double nOldWjets = TMath::Max(histo_Wjets->GetSumOfWeights(),0.000001);
  for(int i=1; i<=histo_Wjets->GetNbinsX(); i++){
    if(histo_Wjets->GetBinContent(i)                    < 0) {histo_Wjets		    ->SetBinContent(i,0.000001);histo_Wjets		      ->SetBinError(i,0.000001);}
    if(histo_Wjets_WUp->GetBinContent(i) < 0) {histo_Wjets_WUp->SetBinContent(i,0.000001);histo_Wjets_WUp->SetBinError(i,0.000001);}
  }
  histo_Wjets                   ->Scale(nOldWjets/histo_Wjets                   ->GetSumOfWeights());
  histo_Wjets_WUp->Scale(nOldWjets/histo_Wjets_WUp->GetSumOfWeights());

  // WWqcd treatment
  if(histo_WWqcd->GetSumOfWeights() > 0){
    double nOldWWQCD[8] = {TMath::Max(histo_WWqcd->GetSumOfWeights(),0.000001),
                           TMath::Max(histo_WWqcd_LepEffUp->GetSumOfWeights(),0.000001),TMath::Max(histo_WWqcd_LepEffDown->GetSumOfWeights(),0.000001),
                           TMath::Max(histo_WWqcd_JESUp->GetSumOfWeights(),0.000001),   TMath::Max(histo_WWqcd_JESDown->GetSumOfWeights(),0.000001),
			   TMath::Max(histo_WWqcd_LepResUp->GetSumOfWeights(),0.000001),TMath::Max(histo_WWqcd_LepResDown->GetSumOfWeights(),0.000001),
			   TMath::Max(histo_WWqcd_METResUp->GetSumOfWeights(),0.000001)};
    for(int i=1; i<=histo_WWqcd->GetNbinsX(); i++){
      if(histo_WWqcd			    ->GetBinContent(i) < 0) {histo_WWqcd			  ->SetBinContent(i,0.000001);histo_WWqcd			   ->SetBinError(i,0.000001);}
      if(histo_WWqcd_LepEffUp  ->GetBinContent(i) < 0) {histo_WWqcd_LepEffUp  ->SetBinContent(i,0.000001);histo_WWqcd_LepEffUp  ->SetBinError(i,0.000001);}
      if(histo_WWqcd_LepEffDown->GetBinContent(i) < 0) {histo_WWqcd_LepEffDown->SetBinContent(i,0.000001);histo_WWqcd_LepEffDown->SetBinError(i,0.000001);}
      if(histo_WWqcd_JESUp     ->GetBinContent(i) < 0) {histo_WWqcd_JESUp	  ->SetBinContent(i,0.000001);histo_WWqcd_JESUp     ->SetBinError(i,0.000001);}
      if(histo_WWqcd_JESDown   ->GetBinContent(i) < 0) {histo_WWqcd_JESDown   ->SetBinContent(i,0.000001);histo_WWqcd_JESDown   ->SetBinError(i,0.000001);}
      if(histo_WWqcd_LepResUp  ->GetBinContent(i) < 0) {histo_WWqcd_LepResUp  ->SetBinContent(i,0.000001);histo_WWqcd_LepResUp  ->SetBinError(i,0.000001);}
      if(histo_WWqcd_LepResDown->GetBinContent(i) < 0) {histo_WWqcd_LepResDown->SetBinContent(i,0.000001);histo_WWqcd_LepResDown->SetBinError(i,0.000001);}
      if(histo_WWqcd_METResUp  ->GetBinContent(i) < 0) {histo_WWqcd_METResUp  ->SetBinContent(i,0.000001);histo_WWqcd_METResUp  ->SetBinError(i,0.000001);}
    }
    histo_WWqcd			         ->Scale(nOldWWQCD[0]/histo_WWqcd			         ->GetSumOfWeights());
    histo_WWqcd_LepEffUp  ->Scale(nOldWWQCD[1]/histo_WWqcd_LepEffUp  ->GetSumOfWeights());
    histo_WWqcd_LepEffDown->Scale(nOldWWQCD[2]/histo_WWqcd_LepEffDown->GetSumOfWeights());
    histo_WWqcd_JESUp     ->Scale(nOldWWQCD[3]/histo_WWqcd_JESUp     ->GetSumOfWeights());
    histo_WWqcd_JESDown   ->Scale(nOldWWQCD[4]/histo_WWqcd_JESDown   ->GetSumOfWeights());
    histo_WWqcd_LepResUp  ->Scale(nOldWWQCD[5]/histo_WWqcd_LepResUp  ->GetSumOfWeights());
    histo_WWqcd_LepResDown->Scale(nOldWWQCD[6]/histo_WWqcd_LepResDown->GetSumOfWeights());
    histo_WWqcd_METResUp  ->Scale(nOldWWQCD[7]/histo_WWqcd_METResUp  ->GetSumOfWeights());
  } else {
    histo_WWqcd			         ->Scale(0.0);
  }

  for(int i=1; i<=histo_WWewk->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_WWewk_WWewkStatUp	    ->SetBinContent(i,TMath::Max(histo_WWewk    ->GetBinContent(i)+factorUp  *histo_WWewk	 ->GetBinError(i),0.000001));
    histo_WWewk_WWewkStatDown        ->SetBinContent(i,TMath::Max(histo_WWewk    ->GetBinContent(i)+factorDown*histo_WWewk	 ->GetBinError(i),0.000001));
    histo_WWqcd_WWqcdStatUp	    ->SetBinContent(i,TMath::Max(histo_WWqcd    ->GetBinContent(i)+factorUp  *histo_WWqcd	 ->GetBinError(i),0.000001));
    histo_WWqcd_WWqcdStatDown        ->SetBinContent(i,TMath::Max(histo_WWqcd    ->GetBinContent(i)+factorDown*histo_WWqcd	 ->GetBinError(i),0.000001));
    histo_WZ_WZStatUp	      	    ->SetBinContent(i,TMath::Max(histo_WZ    	->GetBinContent(i)+factorUp  *histo_WZ   	 ->GetBinError(i),0.000001));
    histo_WZ_WZStatDown        	    ->SetBinContent(i,TMath::Max(histo_WZ    	->GetBinContent(i)+factorDown*histo_WZ   	 ->GetBinError(i),0.000001));
    histo_WS_WSStatUp	      	    ->SetBinContent(i,TMath::Max(histo_WS    	->GetBinContent(i)+factorUp  *histo_WS   	 ->GetBinError(i),0.000001));
    histo_WS_WSStatDown        	    ->SetBinContent(i,TMath::Max(histo_WS    	->GetBinContent(i)+factorDown*histo_WS   	 ->GetBinError(i),0.000001));
    histo_VVV_VVVStatUp        	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorUp  *histo_VVV  	 ->GetBinError(i),0.000001));
    histo_VVV_VVVStatDown      	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorDown*histo_VVV  	 ->GetBinError(i),0.000001));
    histo_Wjets_WjetsStatUp    	    ->SetBinContent(i,TMath::Max(histo_Wjets    ->GetBinContent(i)+factorUp  *histo_Wjets        ->GetBinError(i),0.000001));
    histo_Wjets_WjetsStatDown  	    ->SetBinContent(i,TMath::Max(histo_Wjets    ->GetBinContent(i)+factorDown*histo_Wjets        ->GetBinError(i),0.000001));
  }
  double mean,up,diff;

  if(showSignalOnly == false) {
    printf("nuisance WZ | Wj: %f/%f | %f/%f\n",histo_WZ   ->GetSumOfWeights(),histo_WZ_CMS_WZNLOUp  ->GetSumOfWeights(),
                                               histo_Wjets->GetSumOfWeights(),histo_Wjets_WUp->GetSumOfWeights());
  }
  histo_Wjets_WUp->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_WUp->GetSumOfWeights());
  histo_WZ_CMS_WZNLOUp  ->Scale(histo_WZ   ->GetSumOfWeights()/histo_WZ_CMS_WZNLOUp  ->GetSumOfWeights());

  for(int i=1; i<=histo_WWewk->GetNbinsX(); i++){
    // METRes
    mean = histo_WWewk			      ->GetBinContent(i);
    up   = histo_WWewk_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWewk_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWewk_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WWqcd			      ->GetBinContent(i);
    up   = histo_WWqcd_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWqcd_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWqcd_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WZ		       ->GetBinContent(i);
    up   = histo_WZ_CMS_WZNLOUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_WZNLODown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_WZNLODown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WS			   ->GetBinContent(i);
    up   = histo_WS_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WS_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WS_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_METResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_METResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_METResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    // LepRes
    mean = histo_WWewk			      ->GetBinContent(i);
    up   = histo_WWewk_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWewk_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWewk_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WWqcd			      ->GetBinContent(i);
    up   = histo_WWqcd_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWqcd_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWqcd_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WZ		       ->GetBinContent(i);
    up   = histo_WZ_CMS_WZNLOUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_WZNLODown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_WZNLODown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WS			   ->GetBinContent(i);
    up   = histo_WS_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WS_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WS_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_LepResUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_LepResDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_LepResDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    // JES
    mean = histo_WWewk			   ->GetBinContent(i);
    up   = histo_WWewk_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWewk_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWewk_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WWqcd			   ->GetBinContent(i);
    up   = histo_WWqcd_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWqcd_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWqcd_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WZ			->GetBinContent(i);
    up   = histo_WZ_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WZ		       ->GetBinContent(i);
    up   = histo_WZ_CMS_WZNLOUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_WZNLODown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_WZNLODown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WS			->GetBinContent(i);
    up   = histo_WS_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WS_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WS_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			 ->GetBinContent(i);
    up   = histo_VVV_JESUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_JESDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_JESDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    // GEN
    mean = histo_Wjets 		         ->GetBinContent(i);
    up   = histo_Wjets_WUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Wjets_WDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Wjets_WDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WS 		       ->GetBinContent(i);
    up   = histo_WS_WSUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WS_WSDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WS_WSDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  }
  histo_Wjets_WDown->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_WDown->GetSumOfWeights());
  histo_WZ_CMS_WZNLODown  ->Scale(histo_WZ   ->GetSumOfWeights()/histo_WZ_CMS_WZNLODown  ->GetSumOfWeights());

  //----------------------------------------------------------------------------
  // Produce output cards for shape-based analyses
  //----------------------------------------------------------------------------
  if(showSignalOnly == false){
  char outputLimits[200];
  sprintf(outputLimits,"wwss%2s_nsign%1d.input_%4s.root",finalStateName,signSel,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data	 ->Write();
  histo_WWewk	 ->Write();
  histo_WWqcd	 ->Write();
  histo_WZ	 ->Write();
  histo_WS	 ->Write();
  histo_VVV	 ->Write();
  histo_Wjets	 ->Write();

  cout << histo_Data	 ->GetSumOfWeights() << " ";
  cout << histo_WWewk	 ->GetSumOfWeights() << " ";
  cout << histo_WWqcd	 ->GetSumOfWeights() << " ";
  cout << histo_WZ	 ->GetSumOfWeights() << " ";
  cout << histo_WS	 ->GetSumOfWeights() << " ";
  cout << histo_VVV	 ->GetSumOfWeights() << " ";
  cout << histo_Wjets	 ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_WWewk_WWewkStatUp	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk	->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_WWewkStatUp	     ->GetBinContent(i)/histo_WWewk   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_WWewkStatDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk	->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_WWewkStatDown      ->GetBinContent(i)/histo_WWewk   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_WWqcdStatUp	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_WWqcdStatUp  	->GetBinContent(i)/histo_WWqcd   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_WWqcdStatDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_WWqcdStatDown      ->GetBinContent(i)/histo_WWqcd   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_WZStatUp	  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_WZStatUp	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_WZStatDown	 	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_WZStatDown	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_WSStatUp	  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_WSStatUp	->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_WSStatDown	 	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_WSStatDown	->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_VVVStatUp	  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_VVVStatUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_VVVStatDown   	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_VVVStatDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_WjetsStatUp  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WjetsStatUp  ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_WjetsStatDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WjetsStatDown->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEff\n");
  histo_WWewk_LepEffUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_LepEffUp	  ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_LepEffDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_LepEffDown   ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_LepEffUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_LepEffUp	     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_LepEffDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_LepEffDown	->GetBinContent(i)/histo_WWqcd     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_LepEffUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepEffUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_LepEffDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepEffDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_LepEffUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_LepEffUp     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_LepEffDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_LepEffDown   ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepEffUp          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepEffUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepEffDown        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepEffDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LetRes\n");
  histo_WWewk_LepResUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_LepResUp	  ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_LepResDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_LepResDown   ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_LepResUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_LepResUp	     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_LepResDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_LepResDown	->GetBinContent(i)/histo_WWqcd     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_LepResUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepResUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_LepResDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_LepResDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_LepResUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_LepResUp     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_LepResDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_LepResDown   ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepResUp          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepResUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_LepResDown        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_LepResDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties METRes\n");
  histo_WWewk_METResUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_METResUp	     ->GetBinContent(i)/histo_WWewk	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_METResDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_METResDown	->GetBinContent(i)/histo_WWewk     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_METResUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_METResUp	     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_METResDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_METResDown	->GetBinContent(i)/histo_WWqcd     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_METResUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_METResUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_METResDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_METResDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_METResUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_METResUp     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_METResDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_METResDown   ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_METResUp          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_METResUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_METResDown        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_METResDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_WWewk_JESUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_JESUp	  ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_JESDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_JESDown   ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_JESUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_JESUp     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_JESDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_JESDown   ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_JESUp              ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_JESUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_JESDown            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_JESDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_JESUp              ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_JESUp	     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_JESDown            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_JESDown      ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_JESUp             ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_JESUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_JESDown           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_JESDown	     ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties GEN\n");
  histo_Wjets_WUp	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WUp       ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_WDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_WDown     ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_WZNLOUp            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_WZNLOUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_WZNLODown          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_WZNLODown       ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_WSUp            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_WSUp         ->GetBinContent(i)/histo_WS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_WSDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_WSDown       ->GetBinContent(i)/histo_WS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  for(int nb=1; nb<=nBin; nb++){
    if(doAQGCsAna == true){
      stringstream ss;
      ss << nb;
    
      TH2D *th2d  = new TH2D(string("aQGC_scaling"+ss.str()).c_str(),string("aQGC_scaling"+ss.str()).c_str(),11,-11,11,11,-11,11);

      for(unsigned int a = 0; a < grid_points.size(); a++){
        //histo_grid[nb][a] = histo_WWewk_anom[a]->GetBinContent(nb);
        assert(histo_WWewk_anom[0]->GetBinContent(nb) > 0);
        //assert(histo_grid[nb][0]>0);
        th2d->SetBinContent(th2d->GetXaxis()->FindFixBin(grid_points[a].first), th2d->GetYaxis()->FindFixBin(grid_points[a].second), histo_WWewk_anom[a]->GetBinContent(nb)/histo_WWewk_anom[0]->GetBinContent(nb));
      }
      th2d_outfile->cd();
      th2d->Write();
    }

    double systNLO[3] = {1.0,1.0,1.0}; // WZ, WS, Wjets
    if     (histo_WZ   ->GetBinContent(nb) > 0 && histo_WZ_CMS_WZNLOUp    ->GetBinContent(nb) > 0) systNLO[0] = histo_WZ_CMS_WZNLOUp    ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ   ->GetBinContent(nb) > 0 && histo_WZ_CMS_WZNLODown  ->GetBinContent(nb) > 0) systNLO[0] = histo_WZ   ->GetBinContent(nb)/histo_WZ_CMS_WZNLODown  ->GetBinContent(nb);
    if     (histo_WS   ->GetBinContent(nb) > 0 && histo_WS_WSUp    ->GetBinContent(nb) > 0) systNLO[1] = histo_WS_WSUp    ->GetBinContent(nb)/histo_WS   ->GetBinContent(nb);
    else if(histo_WS   ->GetBinContent(nb) > 0 && histo_WS_WSDown  ->GetBinContent(nb) > 0) systNLO[1] = histo_WS   ->GetBinContent(nb)/histo_WS_WSDown  ->GetBinContent(nb);
    if     (histo_Wjets->GetBinContent(nb) > 0 && histo_Wjets_WUp  ->GetBinContent(nb) > 0) systNLO[2] = histo_Wjets_WUp  ->GetBinContent(nb)/histo_Wjets->GetBinContent(nb);
    else if(histo_Wjets->GetBinContent(nb) > 0 && histo_Wjets_WDown->GetBinContent(nb) > 0) systNLO[2] = histo_Wjets->GetBinContent(nb)/histo_Wjets_WDown->GetBinContent(nb);

    double systEff[5] = {1.0,1.0,1.0,1.0,1.0};
    if(histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_LepEffUp->GetBinContent(nb) > 0) systEff[0] =  histo_WWewk_LepEffUp->GetBinContent(nb)/histo_WWewk->GetBinContent(nb);
    if(histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_LepEffUp->GetBinContent(nb) > 0) systEff[1] =  histo_WWqcd_LepEffUp->GetBinContent(nb)/histo_WWqcd->GetBinContent(nb);
    if(histo_WZ   ->GetBinContent(nb) > 0 && histo_WZ_LepEffUp	->GetBinContent(nb) > 0) systEff[2] =  histo_WZ_LepEffUp   ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    if(histo_WS   ->GetBinContent(nb) > 0 && histo_WS_LepEffUp	->GetBinContent(nb) > 0) systEff[3] =  histo_WS_LepEffUp   ->GetBinContent(nb)/histo_WS   ->GetBinContent(nb);
    if(histo_VVV  ->GetBinContent(nb) > 0 && histo_VVV_LepEffUp	->GetBinContent(nb) > 0) systEff[4] =  histo_VVV_LepEffUp  ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);

    double systLep[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_LepResUp  ->GetBinContent(nb) > 0) systLep[0] =  histo_WWewk_LepResUp  ->GetBinContent(nb)/histo_WWewk->GetBinContent(nb);
    else if(histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_LepResDown->GetBinContent(nb) > 0) systLep[0] =  histo_WWewk->GetBinContent(nb)/histo_WWewk_LepResDown->GetBinContent(nb);
    if     (histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_LepResUp  ->GetBinContent(nb) > 0) systLep[1] =  histo_WWqcd_LepResUp  ->GetBinContent(nb)/histo_WWqcd->GetBinContent(nb);
    else if(histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_LepResDown->GetBinContent(nb) > 0) systLep[1] =  histo_WWqcd->GetBinContent(nb)/histo_WWqcd_LepResDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_LepResUp     ->GetBinContent(nb) > 0) systLep[2] =  histo_WZ_LepResUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_LepResDown   ->GetBinContent(nb) > 0) systLep[2] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_LepResDown   ->GetBinContent(nb);
    if     (histo_WS->GetBinContent(nb)    > 0 && histo_WS_LepResUp     ->GetBinContent(nb) > 0) systLep[3] =  histo_WS_LepResUp     ->GetBinContent(nb)/histo_WS   ->GetBinContent(nb);
    else if(histo_WS->GetBinContent(nb)    > 0 && histo_WS_LepResDown   ->GetBinContent(nb) > 0) systLep[3] =  histo_WS   ->GetBinContent(nb)/histo_WS_LepResDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_LepResUp    ->GetBinContent(nb) > 0) systLep[4] =  histo_VVV_LepResUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_LepResDown  ->GetBinContent(nb) > 0) systLep[4] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_LepResDown  ->GetBinContent(nb);

    double systMet[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_METResUp  ->GetBinContent(nb) > 0) systMet[0] =  histo_WWewk_METResUp  ->GetBinContent(nb)/histo_WWewk->GetBinContent(nb);
    else if(histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_METResDown->GetBinContent(nb) > 0) systMet[0] =  histo_WWewk->GetBinContent(nb)/histo_WWewk_METResDown->GetBinContent(nb);
    if     (histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_METResUp  ->GetBinContent(nb) > 0) systMet[1] =  histo_WWqcd_METResUp  ->GetBinContent(nb)/histo_WWqcd->GetBinContent(nb);
    else if(histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_METResDown->GetBinContent(nb) > 0) systMet[1] =  histo_WWqcd->GetBinContent(nb)/histo_WWqcd_METResDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_METResUp     ->GetBinContent(nb) > 0) systMet[2] =  histo_WZ_METResUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_METResDown   ->GetBinContent(nb) > 0) systMet[2] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_METResDown   ->GetBinContent(nb);
    if     (histo_WS->GetBinContent(nb)    > 0 && histo_WS_METResUp     ->GetBinContent(nb) > 0) systMet[3] =  histo_WS_METResUp     ->GetBinContent(nb)/histo_WS   ->GetBinContent(nb);
    else if(histo_WS->GetBinContent(nb)    > 0 && histo_WS_METResDown   ->GetBinContent(nb) > 0) systMet[3] =  histo_WS   ->GetBinContent(nb)/histo_WS_METResDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_METResUp    ->GetBinContent(nb) > 0) systMet[4] =  histo_VVV_METResUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_METResDown  ->GetBinContent(nb) > 0) systMet[4] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_METResDown  ->GetBinContent(nb);

    double systJes[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_JESUp  ->GetBinContent(nb) > 0) systJes[0] =  histo_WWewk_JESUp  ->GetBinContent(nb)/histo_WWewk->GetBinContent(nb);
    else if(histo_WWewk->GetBinContent(nb) > 0 && histo_WWewk_JESDown->GetBinContent(nb) > 0) systJes[0] =  histo_WWewk->GetBinContent(nb)/histo_WWewk_JESDown->GetBinContent(nb);
    if     (histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_JESUp  ->GetBinContent(nb) > 0) systJes[1] =  histo_WWqcd_JESUp  ->GetBinContent(nb)/histo_WWqcd->GetBinContent(nb);
    else if(histo_WWqcd->GetBinContent(nb) > 0 && histo_WWqcd_JESDown->GetBinContent(nb) > 0) systJes[1] =  histo_WWqcd->GetBinContent(nb)/histo_WWqcd_JESDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_JESUp     ->GetBinContent(nb) > 0) systJes[2] =  histo_WZ_JESUp     ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)    > 0 && histo_WZ_JESDown   ->GetBinContent(nb) > 0) systJes[2] =  histo_WZ   ->GetBinContent(nb)/histo_WZ_JESDown   ->GetBinContent(nb);
    if     (histo_WS->GetBinContent(nb)    > 0 && histo_WS_JESUp     ->GetBinContent(nb) > 0) systJes[3] =  histo_WS_JESUp     ->GetBinContent(nb)/histo_WS   ->GetBinContent(nb);
    else if(histo_WS->GetBinContent(nb)    > 0 && histo_WS_JESDown   ->GetBinContent(nb) > 0) systJes[3] =  histo_WS   ->GetBinContent(nb)/histo_WS_JESDown   ->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_JESUp    ->GetBinContent(nb) > 0) systJes[4] =  histo_VVV_JESUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)   > 0 && histo_VVV_JESDown  ->GetBinContent(nb) > 0) systJes[4] =  histo_VVV  ->GetBinContent(nb)/histo_VVV_JESDown  ->GetBinContent(nb);

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_wwss%2s_nsign%1d_shape_%4s_Bin%d.txt",finalStateName,signSel,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin wwss%2s%4s%d wwss%2s%4s%d wwss%2s%4s%d wwss%2s%4s%d wwss%2s%4s%d wwss%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process WWewk WWqcd WZ WS VVV Wjets\n");
    newcardShape << Form("process 0 1 2 3 4 5\n");
    newcardShape << Form("rate %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",histo_WWewk->GetBinContent(nb),histo_WWqcd->GetBinContent(nb),histo_WZ->GetBinContent(nb),histo_WS->GetBinContent(nb),histo_VVV->GetBinContent(nb),histo_Wjets->GetBinContent(nb));
    newcardShape << Form("lumi_%4s                                 lnN %5.3f %5.3f %5.3f %5.3f %5.3f  -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE);		     
    newcardShape << Form("%s                                       lnN %5.3f %5.3f %5.3f %5.3f %5.3f  -  \n",effName,systEff[0],systEff[1],systEff[2],systEff[3],systEff[4]);
    newcardShape << Form("%s                                       lnN %5.3f %5.3f %5.3f %5.3f %5.3f  -  \n",momName,systLep[0],systLep[1],systLep[2],systLep[3],systLep[4]);
    newcardShape << Form("CMS_scale_met                            lnN %5.3f %5.3f %5.3f %5.3f %5.3f  -  \n",systMet[0],systMet[1],systMet[2],systMet[3],systMet[4]);
    newcardShape << Form("CMS_scale_j                              lnN %5.3f %5.3f %5.3f %5.3f %5.3f  -  \n",systJes[0],systJes[1],systJes[2],systJes[3],systJes[4]);		      
    newcardShape << Form("pdf_qqbar                                lnN %5.3f %5.3f %5.3f   -     -    -  \n",pdf_qqbar[0],pdf_qqbar[1],pdf_qqbar[2]);
    newcardShape << Form("QCDscale_WWewk		           lnN %5.3f   -     -	   -     -    -  \n",QCDscale_WWewk);	  
    newcardShape << Form("QCDscale_WWqcd		           lnN   -   %5.3f   -	   -     -    -  \n",QCDscale_WWqcd);	  
    newcardShape << Form("QCDscale_VV		                   lnN   -     -   1.100   -     -    -  \n");  	
    newcardShape << Form("CMS_wwss_WZ3l                            lnN   -     -   %5.3f   -     -    -  \n",syst_WZ3l);		
    newcardShape << Form("CMS_wwss_WZNLO                           lnN   -     -   %5.3f   -     -    -  \n",systNLO[0]);
    newcardShape << Form("CMS_wwss_MVAWS                           lnN   -     -    -    %5.3f   -    -  \n",systNLO[1]);		
    newcardShape << Form("QCDscale_VVV		                   lnN   -     -    -	  -   1.500   -  \n");  	   
    newcardShape << Form("CMS_FakeRate                             lnN   -     -    -	  -	-   %5.3f\n",WjetsSyst);  
    newcardShape << Form("CMS_wwss_MVAW                            lnN   -     -    -	  -	-   %5.3f\n",systNLO[2]);
    if(histo_WWewk->GetBinContent(nb) > 0)
    newcardShape << Form("CMS_wwss%s_MVAWWewkStat_%s_Bin%d         lnN  %5.3f  -    -     -     -      - \n",finalStateName,ECMsb.Data(),nb-1,histo_WWewk_WWewkStatUp->GetBinContent(nb)/histo_WWewk->GetBinContent(nb));
    if(histo_WWqcd->GetBinContent(nb) > 0)
    newcardShape << Form("CMS_wwss%s_MVAWWqcdStat_%s_Bin%d         lnN    -   %5.3f -	  -     -      - \n",finalStateName,ECMsb.Data(),nb-1,histo_WWqcd_WWqcdStatUp->GetBinContent(nb)/histo_WWqcd->GetBinContent(nb));
    if(histo_WZ->GetBinContent(nb) > 0)
    newcardShape << Form("CMS_wwss%s_MVAWZStat_%s_Bin%d            lnN    -    -   %5.3f  -     -      - \n",finalStateName,ECMsb.Data(),nb-1,histo_WZ_WZStatUp    ->GetBinContent(nb)/histo_WZ   ->GetBinContent(nb));
    if(histo_WS->GetBinContent(nb) > 0)
    newcardShape << Form("CMS_wwss%s_MVAWSStat_%s_Bin%d            lnN    -    -    -	%5.3f   -      - \n",finalStateName,ECMsb.Data(),nb-1,histo_WS_WSStatUp    ->GetBinContent(nb)/histo_WS   ->GetBinContent(nb));
    if(histo_VVV->GetBinContent(nb) > 0)
    newcardShape << Form("CMS_wwss%s_MVAVVVStat_%s_Bin%d           lnN    -    -    -	  -   %5.3f    - \n",finalStateName,ECMsb.Data(),nb-1,histo_VVV_VVVStatUp    ->GetBinContent(nb)/histo_VVV  ->GetBinContent(nb));
    if(histo_Wjets->GetBinContent(nb) > 0)
    newcardShape << Form("CMS_wwss%s_MVAWjetsStat_%s_Bin%d         lnN    -    -    -	  -     -   %5.3f\n",finalStateName,ECMsb.Data(),nb-1,histo_Wjets_WjetsStatUp->GetBinContent(nb)/histo_Wjets->GetBinContent(nb));
    newcardShape.close();
  }
  } // if showSignalOnly == true

  return;
}

void scaleFactor_WS(LorentzVector l,int lq, int ld, int mcld, double val[2]){
//---------------------------------------------------------------------
// |eta|        data                  mc                    factor
//---------------------------------------------------------------------
//0.0-0.5 0.000109 +/- 0.000017 | 0.000106 +/- 0.000017 ==> 1.028 +/- 0.160
//0.5-1.0 0.000210 +/- 0.000025 | 0.000188 +/- 0.000025 ==> 1.117 +/- 0.133
//1.0-1.5 0.001302 +/- 0.000087 | 0.001022 +/- 0.000087 ==> 1.274 +/- 0.085
//1.5-2.0 0.003437 +/- 0.000193 | 0.002562 +/- 0.000193 ==> 1.342 +/- 0.075
//2.0-2.5 0.003270 +/- 0.000241 | 0.002245 +/- 0.000240 ==> 1.456 +/- 0.107
// additional 10% uncertainty for the overall normalization
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
}
