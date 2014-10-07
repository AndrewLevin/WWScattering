#include <TAxis.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "tdrstyle.C"


Bool_t verbose = false;

void DrawTLatex(Double_t x, Double_t y, Double_t tsize, const char* text);


//------------------------------------------------------------------------------
//
// ratio
//    0 = don't take ratio
//    1 = simple ratio
//    2 = ratio in units of sigma
//    3 = sigma
//
//------------------------------------------------------------------------------
void PlotLimitForVBS(string  limitFiles   = "inputs/ana_ICHEP_limits_nj_shape7teV_cut8TeV.txt",
	       string  outputPrefix = "combined",
	       string  luminosity   = "5.1 fb^{-1} (8 TeV) + 4.9 fb^{-1} (7 TeV)",
	       Float_t mhmin        = 110,
	       Float_t mhmax        = 160,
	       Int_t   setLogx      = 0,
	       Int_t   setLogy      = 1,
	       string  title        = "H #rightarrow WW #rightarrow 2l2#nu",
	       Bool_t  drawObserved = 1,
	       Int_t   ratio        = 0,
	       string  format       = "pdf")
{
  setTDRStyle();

  gSystem->mkdir(format.c_str(), kTRUE);

  // Get input files
  //----------------------------------------------------------------------------
  stringstream fnames(limitFiles);
  string fname;
  vector<string> LimitFile;

  while (getline(fnames, fname, ',')) {
    if (verbose) printf(" Using file(s) %s\n", fname.c_str());
    LimitFile.push_back(fname);
  }


  // Read in the nominal limits
  //----------------------------------------------------------------------------
  vector<float> vMass;
  vector<float> vObsLimit; 
  vector<float> vMeanExpLimit; 
  vector<float> vMedianExpLimit; 
  vector<float> vExpLim68Down; 
  vector<float> vExpLim68Up; 
  vector<float> vExpLim95Down; 
  vector<float> vExpLim95Up;

  float Mass;
  float ObsLimit; 
  float MeanExpLimit; 
  float MedianExpLimit; 
  float ExpLim68Down; 
  float ExpLim68Up; 
  float ExpLim95Down; 
  float ExpLim95Up;

  ifstream indata;
  indata.open(LimitFile[0].c_str());
  if (!indata) {
    cerr << " Error: file could not be opened" << endl;
    return;
  }

  while (indata
	 >> Mass
	 >> ObsLimit
	 >> MeanExpLimit
	 >> MedianExpLimit
	 >> ExpLim95Down
	 >> ExpLim68Down
	 >> ExpLim68Up
	 >> ExpLim95Up) {

    vMass          .push_back(Mass);
    vObsLimit      .push_back(ObsLimit); 
    vMeanExpLimit  .push_back(MeanExpLimit); 
    vMedianExpLimit.push_back(MedianExpLimit); 
    vExpLim68Down  .push_back(ExpLim68Down); 
    vExpLim68Up    .push_back(ExpLim68Up); 
    vExpLim95Down  .push_back(ExpLim95Down); 
    vExpLim95Up    .push_back(ExpLim95Up);
  }

  UInt_t npoints = vMass.size();


  // Draw
  //----------------------------------------------------------------------------
  TString ctitle = Form("%s_from%.0fto%.0f_logx%d_logy%d",
			outputPrefix.c_str(),
			mhmin,
			mhmax,
			setLogx,
			setLogy);

  //TCanvas* canvas = new TCanvas(ctitle, ctitle);
  TCanvas* canvas = new TCanvas(ctitle, ctitle, 600, 600, 600, 600);

  canvas->SetLeftMargin  (1.30 * canvas->GetLeftMargin());
  canvas->SetRightMargin (2.10 * canvas->GetRightMargin());
  canvas->SetBottomMargin(1.35 * canvas->GetBottomMargin());

  canvas->SetLogx(setLogx);
  canvas->SetLogy(setLogy);

  float min = +9999;
  float max = -9999;


  // Expected Limit
  //----------------------------------------------------------------------------
  float x   [npoints];
  float ex  [npoints];
  float y   [npoints];
  float yu68[npoints];
  float yd68[npoints];
  float yu95[npoints];
  float yd95[npoints]; 
  
  for (UInt_t i=0; i<npoints; ++i) {

    x [i] = vMass.at(i);
    ex[i] = 0; 

    y   [i] = vMedianExpLimit.at(i);
    yu68[i] = vExpLim68Up.at(i) - y[i];
    yu95[i] = vExpLim95Up.at(i) - y[i];   
    yd68[i] = y[i] - vExpLim68Down.at(i);  
    yd95[i] = y[i] - vExpLim95Down.at(i);  

    if (ratio == 1) {
      y   [i] /= vMedianExpLimit.at(i);
      yu68[i] /= vMedianExpLimit.at(i);
      yu95[i] /= vMedianExpLimit.at(i);
      yd68[i] /= vMedianExpLimit.at(i);
      yd95[i] /= vMedianExpLimit.at(i);
    }
    else if (ratio == 2) {
      y   [i] = 0;
      yu68[i] = 1;
      yu95[i] = 2;
      yd68[i] = 1;
      yd95[i] = 2;
    }

    if (y[i] + yu95[i] > max) max = y[i] + yu95[i];
    if (y[i] - yd95[i] < min) min = y[i] - yd95[i]; // 0.01;
  }
  if(setLogy == 0) min =0.01;


  TGraph*            ExpLim    = new TGraph           (npoints, x, y);
  TGraphAsymmErrors* ExpBand95 = new TGraphAsymmErrors(npoints, x, y, ex, ex, yd95, yu95);
  TGraphAsymmErrors* ExpBand68 = new TGraphAsymmErrors(npoints, x, y, ex, ex, yd68, yu68);

  ExpBand95->GetXaxis()->SetRangeUser(mhmin, mhmax);

  // Cosmetics
  //----------------------------------------------------------------------------
  ExpBand95->SetTitle("");

  //ExpBand95->GetXaxis()->SetTitle("Higgs boson mass (GeV)");
  ExpBand95->GetXaxis()->SetTitle("m_{H^{#pm#pm}} (GeV)");
    
  if (ratio == 0) ExpBand95->GetYaxis()->SetTitle("95% CL limit on #sigma/#sigma_{SM}");
  if (ratio == 1) ExpBand95->GetYaxis()->SetTitle("ratio to expected");
  if (ratio == 2) ExpBand95->GetYaxis()->SetTitle("(observed - expected) / 1#sigma");
  if (ratio == 3) ExpBand95->GetYaxis()->SetTitle("#sigma X BR(H #rightarrow WW #rightarrow 2l2#nu) (pb)");
  if (ratio == 4) ExpBand95->GetYaxis()->SetTitle("#sigma X BR(H #rightarrow invisible)/#sigma_{SM}");
  if (ratio == 5) ExpBand95->GetYaxis()->SetTitle("#sigma X #Beta(H^{#pm#pm} #rightarrow W^{#pm}W^{#pm}) (fb)");

  ExpBand95->GetXaxis()->SetTitleOffset(1.05);
  ExpBand95->GetYaxis()->SetTitleOffset(1.05);
  ExpBand95->GetYaxis()->SetNdivisions(505);

  ExpBand68->SetFillColor(211); 
  ExpBand68->SetLineColor( 10);

  ExpBand95->SetFillColor(90); 
  ExpBand95->SetLineColor(10);

  ExpLim->SetLineStyle(2);
  ExpLim->SetLineWidth(2);

  ExpBand95->Draw("a3");
  ExpBand68->Draw("3");
  ExpLim   ->Draw("l");


  // Observed limit
  //----------------------------------------------------------------------------
  TGraph* ObsLim = NULL;

  if (drawObserved) {

    float yObs[npoints];    

    for (UInt_t i=0; i<npoints; ++i) {

      yObs[i] = vObsLimit.at(i);

      if (ratio == 1) {
	yObs[i] /= vMedianExpLimit.at(i);
      }
      else if (ratio == 2) {
        float onesigma =
	  (fabs(log(vExpLim68Down.at(i)) - log(vMedianExpLimit.at(i))) +
	   fabs(log(vExpLim68Up.at(i))   - log(vMedianExpLimit.at(i))) +
	   fabs(log(vExpLim95Down.at(i)) - log(vMedianExpLimit.at(i))) +
	   fabs(log(vExpLim95Up.at(i))   - log(vMedianExpLimit.at(i)))) / 6;
	
        yObs[i] = (log(vObsLimit.at(i)) - log(vMedianExpLimit.at(i))) / onesigma;
      }

      if (yObs[i] > max) max = yObs[i];
      if (yObs[i] < min) min = yObs[i];
    }

    ObsLim = new TGraph(npoints, x, yObs);
    ObsLim->SetLineWidth(4);
    ObsLim->SetMarkerStyle(kFullCircle);
    ObsLim->Draw("l");
  }


  // y-axis
  //----------------------------------------------------------------------------
  if (canvas->GetLogy())
    ExpBand95->GetYaxis()->SetRangeUser(0.03, max+400);
  else
    ExpBand95->GetYaxis()->SetRangeUser(min-0.0, max+1);


  // canvas dimensions
  //----------------------------------------------------------------------------
  canvas->Update();

  Float_t uxmin = canvas->GetUxmin();
  Float_t uxmax = canvas->GetUxmax();
  Float_t uymin = canvas->GetUymin();


  // x-axis ticks
  //----------------------------------------------------------------------------
  if (canvas->GetLogx()) {

    ExpBand95->GetXaxis()->SetNdivisions(0);

    TLine tick;

    tick.SetLineWidth(1);
    tick.SetLineColor(1);

    for (Int_t i=100; i<=600; i+=10) {
	
      if (i < mhmin || i > mhmax) continue;

      Float_t xx = i;

      if (canvas->GetLogy())
	tick.DrawLine(xx, pow(10,uymin), xx, pow(10,uymin) + (i%100 == 0 ? 0.025 : 0.005));
      else
	tick.DrawLine(xx, uymin, xx, uymin + (i%100 == 0 ? 0.75 : 0.3));
    }

    // TLatex
    //--------------------------------------------------------------------------
    Float_t ylatex = (canvas->GetLogy()) ? pow(10,uymin) - 0.01 : uymin - 0.75;

    Float_t xbins[6] = {100, 200, 300, 400, 500, 600};

    while (mhmin > xbins[0]) xbins[0] += 10;

    for (Int_t i=0; i<6; i++) {
      
      if (xbins[i] < mhmin || xbins[i] > mhmax) continue;
    
      TLatex* latex = new TLatex(xbins[i], ylatex, Form("%.0f", xbins[i]));
      
      latex->SetTextAlign(  22);
      latex->SetTextFont (  42);
      latex->SetTextSize (0.040);
    
      latex->Draw("same");
    }
  }

  // Cosmetics
  //----------------------------------------------------------------------------
  //DrawTLatex(0.94, 0.850, 0.050, "CMS preliminary");
  //DrawTLatex(0.94, 0.850, 0.032, "CMS");
  //DrawTLatex(0.94, 0.795, 0.032, title.c_str());
  //DrawTLatex(0.94, 0.740, 0.032, TString("L = "+ luminosity).Data());
  DrawTLatex(0.21, 0.94, 0.04, "#bf{CMS}");
  //DrawTLatex(0.43, 0.94, 0.04, "#bf{CMS (preliminary)}");
  DrawTLatex(0.54, 0.85, 0.04, title.c_str());
  DrawTLatex(0.95, 0.94, 0.04, TString(luminosity).Data());

  TLegend* leg = new TLegend(0.65, 0.66, 0.85, 0.88, "");

  leg->SetBorderSize(    0);
  leg->SetFillColor (    0);
  leg->SetFillStyle (    0);
  leg->SetTextFont  (   42);
  leg->SetTextSize  (0.032);

  if(ObsLim != NULL)
  leg->AddEntry(ObsLim,    " Observed",             "l");
  leg->AddEntry(ExpLim,    " Median expected",      "l");
  leg->AddEntry(ExpBand68, " Expected #pm 1#sigma", "f");
  leg->AddEntry(ExpBand95, " Expected #pm 2#sigma", "f");
  //leg->AddEntry(ExpBand68, " Expected (68%)", "f");
  //leg->AddEntry(ExpBand95, " Expected (95%)", "f");

  leg->Draw("same");


  // Unit line
  //----------------------------------------------------------------------------
  canvas->Update();

  if (!ratio || ratio == 4) {

    TLine *line = NULL;

    if (canvas->GetLogx())
      line = new TLine(pow(10,uxmin), 1, pow(10,uxmax), 1);
    else
      line = new TLine(uxmin, 1, uxmax, 1);

    line->SetLineColor(kRed+1);
    line->SetLineWidth(2);
    line->Draw("same");

    if(ObsLim != NULL)
    ObsLim->Draw("same");
  }
    

  // Save
  //----------------------------------------------------------------------------
  canvas->GetFrame()->DrawClone();
  canvas->RedrawAxis();
  canvas->Update();

  canvas->SaveAs(Form("%s/%s.%s", format.c_str(), ctitle.Data(), format.c_str()));
}


//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t x, Double_t y, Double_t tsize, const char* text)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC();
  tl->SetTextAlign(   32);
  tl->SetTextFont (   42);
  tl->SetTextSize (tsize);

  tl->Draw("same");
}
