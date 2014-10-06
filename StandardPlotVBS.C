#include<vector>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THStack.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TExec.h"
#include <iostream>
#include "TPaveText.h"
#endif

enum samp { iVV,iWWQCD,iWJets,iWW,iVVV,iOther,iWWEWK,iHiggs,nSamples};

float xPos[nSamples+1] = {0.65,0.65,0.65,0.65,0.65,0.40,0.40,0.40,0.40}; 
float yOff[nSamples+1] = {0,1,2,3,4,0,1,2,3};

const Float_t _tsize   = 0.040;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;

//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1D* h)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
        TString coordinate,
        TString title)
{
    axis->SetLabelFont  (   42);
    axis->SetLabelOffset(0.015);
    axis->SetLabelSize  (0.050);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (   42);
    axis->SetTitleOffset(  1.5);
    axis->SetTitleSize  (0.050);

    if (coordinate == "y") axis->SetTitleOffset(1.6);

    axis->SetTitle(title);
}

//------------------------------------------------------------------------------
// THStackAxisFonts
//------------------------------------------------------------------------------
void THStackAxisFonts(THStack* h,
        TString  coordinate,
        TString  title)
{
    TAxis* axis = NULL;

    if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
    if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();

    AxisFonts(axis, coordinate, title);
}

//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
        Float_t y1,
        TH1D*   hist,
        TString label,
        TString option)
{
    TLegend* legend = new TLegend(x1,
            y1,
            x1 + _xoffset,
            y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (_tsize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}


class StandardPlot {

    public: 
        StandardPlot() { _hist.resize(nSamples,0); _data = 0; _breakdown = false; _mass = 0;_isHWWOverlaid = false;}
        void setMCHist   (const samp &s, TH1D * h)  { _hist[s]       = h;  } 
        void setDataHist  (TH1D * h)                { _data          = h;  } 
        void setHWWOverlaid(bool b)                 { _isHWWOverlaid = b;  }

  TH1D* getDataHist() { return _data; }

        void setMass(const int &m) {_mass=m;}

        TH1* DrawAndRebinTo(const int &rebinTo) {

            if(rebinTo == 0) return Draw();
            int rebin = 0, nbins = 0;
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;

                nbins = _hist[i]->GetNbinsX();
            }
            if (nbins == 0) return Draw();

            rebin = nbins / rebinTo;
            while(nbins % rebin != 0) rebin--;
            return Draw(rebin);

        }

        TH1* Draw(const int &rebin=1) {

            Color_t _sampleColor[nSamples];
            _sampleColor[iWWEWK  ] = kBlue+1;
            _sampleColor[iVV 	 ] = kAzure-2;
            _sampleColor[iWWQCD  ] = kGreen+2;
            _sampleColor[iWJets  ] = kGray+1;
            _sampleColor[iWW 	 ] = kAzure-9;
            _sampleColor[iVVV 	 ] = kBlue-1;
	    _sampleColor[iOther  ] = kAzure-9;
            _sampleColor[iHiggs  ] = kMagenta;

            THStack* hstack = new THStack();
	    TH1D* hSum = (TH1D*)_data->Clone();
	    hSum->Rebin(rebin);
	    hSum->Scale(0.0);
	    TAxis *xa;
	    
	    if(_hist[iHiggs] && _hist[iHiggs]->GetSumOfWeights() > 0) _isHWWOverlaid = true;
 
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;
  	    	bool modifyXAxis = false;
		if(modifyXAxis == true){
		  xa =_hist[i]->GetXaxis();
   	          xa->SetLabelSize(0.001);
	          xa->SetLabelColor(4);
  	    	  for(Int_t k=1;k<=_hist[i]->GetNbinsX();++k){
  	    	    xa->SetBinLabel(1 ,"m_{ll}< 250&m_{jj}< 750");
  	    	    xa->SetBinLabel(2 ,"m_{ll}>=250&m_{jj}< 750");
  	    	    xa->SetBinLabel(3 ,"m_{ll}< 250&m_{jj}>=750");
  	    	    xa->SetBinLabel(4 ,"m_{ll}>=250&m_{jj}>=750");
  	    	    xa->SetRangeUser(1,4);
  	    	  }
		}
                _hist[i]->Rebin(rebin);
                _hist[i]->SetLineColor(_sampleColor[i]);

                // signal gets overlaid
                if (i == iWWEWK && _isHWWOverlaid == false) continue;

                if (i == iHiggs) continue;

                _hist[i]->SetFillColor(_sampleColor[i]);
                _hist[i]->SetFillStyle(1001);

                hstack->Add(_hist[i]);
		hSum->Add(_hist[i]);
            }

            if(_hist[iWWEWK]) _hist[iWWEWK]->SetLineWidth(3);
            if(_hist[iHiggs]) _hist[iHiggs]->SetLineWidth(3);
            if(_data) _data->Rebin(rebin);
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);
	    hstack->Draw("hist");

	    bool plotSystErrorBars = true;
	    if(plotSystErrorBars == true) {
  	      TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
              for (int i = 0; i < gsyst->GetN(); ++i) {
                double systBck = 0;
		if(_hist[iWWEWK]) systBck = systBck + 0.107*0.107*_hist[iWWEWK]->GetBinContent(i+1)*_hist[iWWEWK]->GetBinContent(i+1);
		if(_hist[iWWQCD]) systBck = systBck + 0.089*0.089*_hist[iWWQCD]->GetBinContent(i+1)*_hist[iWWQCD]->GetBinContent(i+1);
		if(_hist[   iVV]) systBck = systBck + 0.380*0.380*_hist[   iVV]->GetBinContent(i+1)*_hist[   iVV]->GetBinContent(i+1);
		if(_hist[   iWW]) systBck = systBck + 0.114*0.114*_hist[   iWW]->GetBinContent(i+1)*_hist[   iWW]->GetBinContent(i+1);
		if(_hist[  iVVV]) systBck = systBck + 0.503*0.503*_hist[  iVVV]->GetBinContent(i+1)*_hist[  iVVV]->GetBinContent(i+1);
		if(_hist[iOther]) systBck = systBck + 0.400*0.400*_hist[iOther]->GetBinContent(i+1)*_hist[iOther]->GetBinContent(i+1);
		if(_hist[iWJets]) systBck = systBck + 0.360*0.360*_hist[iWJets]->GetBinContent(i+1)*_hist[iWJets]->GetBinContent(i+1);

                double total = 0;
		if(_hist[iWWEWK]) total = total + _hist[iWWEWK]->GetBinContent(i+1);
		if(_hist[iWWQCD]) total = total + _hist[iWWQCD]->GetBinContent(i+1);
		if(_hist[   iVV]) total = total + _hist[   iVV]->GetBinContent(i+1);
		if(_hist[   iWW]) total = total + _hist[   iWW]->GetBinContent(i+1);
		if(_hist[  iVVV]) total = total + _hist[  iVVV]->GetBinContent(i+1);
		if(_hist[iOther]) total = total + _hist[iOther]->GetBinContent(i+1);
		if(_hist[iWJets]) total = total + _hist[iWJets]->GetBinContent(i+1);

                if(total > 0) systBck = sqrt(systBck)/total;
                gsyst->SetPointEYlow (i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
                gsyst->SetPointEYhigh(i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
	      }
              gsyst->SetFillColor(12);
              gsyst->SetFillStyle(3345);
              gsyst->SetMarkerSize(0);
              gsyst->SetLineWidth(0);
              gsyst->SetLineColor(kWhite);
	      gsyst->Draw("E2same");
              //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
              //setex1->Draw();
	    }
	    else { // plot statistical errors
  	      TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
              for (int i = 0; i < gsyst->GetN(); ++i) {
                gsyst->SetPointEYlow (i, hSum->GetBinError(i+1));
                gsyst->SetPointEYhigh(i, hSum->GetBinError(i+1));
	      }
              gsyst->SetFillColor(12);
              gsyst->SetFillStyle(3345);
              gsyst->SetMarkerSize(0);
              gsyst->SetLineWidth(0);
              gsyst->SetLineColor(kWhite);
	      gsyst->Draw("E2same");
	    }

            if     (_hist[iWWEWK] && _isHWWOverlaid == false)              _hist[iWWEWK]->Draw("hist,same");
	    else if(_hist[iHiggs] && _hist[iHiggs]->GetSumOfWeights() > 0) _hist[iHiggs]->Draw("hist,same");

            if(_data && _data->GetSumOfWeights()) {
 	      bool plotCorrectErrorBars = true;
 	      if(plotCorrectErrorBars == true) {
 	        TGraphAsymmErrors * g = new TGraphAsymmErrors(_data);
                for (int i = 0; i < g->GetN(); ++i) {
 	          double N = g->GetY()[i];
 	          double alpha=(1-0.6827);
 	          double L = (N==0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
 	          double U = ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
 	          g->SetPointEYlow(i,double(N)-L);
 	          if(N >= 0) g->SetPointEYhigh(i, U-double(N));
 	          else       g->SetPointEYhigh(i, 0.0);
 	        }
 	        g->Draw("P");
 	      }
 	      else {
 	        _data->Draw("ep,same");
 	      }
	      //_data->SetBinErrorOption(TH1::kPoisson);
	      //_data->Sumw2(kFALSE);
	      //_data->Draw("E0,same");
            }

            Float_t theMax = hstack->GetMaximum();
            Float_t theMin = hstack->GetMinimum();

            if (_hist[iWWEWK]) {
                if (_hist[iWWEWK]->GetMaximum() > theMax) theMax = _hist[iWWEWK]->GetMaximum();
                if (_hist[iWWEWK]->GetMinimum() < theMin) theMin = _hist[iWWEWK]->GetMinimum();
            }

            if (_data) {

                Float_t dataMax = GetMaximumIncludingErrors(_data);

                if (dataMax > theMax) theMax = dataMax;
            }

            if (gPad->GetLogy()) {
            	hstack->SetMaximum(18 * theMax);
            	hstack->SetMinimum(0.10);
            } else {
              hstack->SetMaximum(1.55 * theMax);
            }

            if(_breakdown) {
                THStackAxisFonts(hstack, "x", _xLabel.Data());
                THStackAxisFonts(hstack, "y", "Events / bin");
                hstack->GetHistogram()->LabelsOption("v");
            } else {
                THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                if(_units.Sizeof() == 1) {
                    THStackAxisFonts(hstack, "x", _xLabel.Data());
                    THStackAxisFonts(hstack, "y", "Events / bin");
                } else {
                    THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                    if(_data->GetBinWidth(0) < 1) THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f %s", _data->GetBinWidth(0),_units.Data()));
		    else                          THStackAxisFonts(hstack, "y", TString::Format("Events / %.0f %s", _data->GetBinWidth(0),_units.Data()));
                }
            }

            // total mess to get it nice, should be redone
            size_t j=0;
            TString higgsLabel = " HWW";
            higgsLabel.Form(" W^{#pm}W^{#pm}jj");

            if(_data->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data,          " Data",    "lp"); j++; }

            if(_hist[iHiggs ] &&_hist[iHiggs ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHiggs ], Form("H^{++}(%3d) #rightarrow W^{#pm}W^{#pm}",_mass)  ,       "f" ); j++; }

            if     (_hist[iWWEWK] && _hist[iWWEWK]->GetSumOfWeights() > 0 && _isHWWOverlaid) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWWEWK], higgsLabel, "f" ); j++; }
            else if(_hist[iWWEWK] && _hist[iWWEWK]->GetSumOfWeights() > 0)		     { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWWEWK], higgsLabel, "l" ); j++; }

	    if(_hist[iVVV   ] &&_hist[iVVV   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVVV   ], " VVV"  ,       "f" ); j++; }

	    if(_hist[iOther   ] &&_hist[iOther   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iOther   ], " Other Bkgs."  ,       "f" ); j++; }

            if(_alternativeOption != 3){
	    if(_hist[iWW    ] &&_hist[iWW    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWW    ], " Wrong sign"  ,"f" ); j++; }
            }
	    else{
	    if(_hist[iWW    ] &&_hist[iWW    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWW    ], " Top-quark+WW"  ,"f" ); j++; }
            }

            if(_hist[iWJets ] &&_hist[iWJets ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWJets ], " Nonprompt",  "f" ); j++; }

            if(_hist[iWWQCD ] &&_hist[iWWQCD ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWWQCD ], " WW DPS",      "f" ); j++; }

            if(_hist[iVV    ] &&_hist[iVV    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVV    ], " WZ",          "f" ); j++; }

            TLatex * CMSLabel = new TLatex (0.18, 0.93, "#bf{CMS}");
            //TLatex * CMSLabel = new TLatex (0.18, 0.93, "#bf{CMS (preliminary)}");
            CMSLabel->SetNDC ();
            CMSLabel->SetTextAlign (10);
            CMSLabel->SetTextFont (42);
            CMSLabel->SetTextSize (_tsize);
            CMSLabel->Draw ("same") ;

           _lumiLabel->Draw ("same") ;

            return hstack->GetHistogram();
        }
        void setLabel(const TString &s) { _xLabel = s; }
        void setUnits(const TString &s) { _units = s; }
        void setBreakdown(const bool &b = true) { _breakdown = b; }
        void setAlternativeOption(const int sel) { _alternativeOption = sel; }
        void setLumiLabel (const std::string &s) {
            _lumiLabel = new TLatex (0.95, 0.93, TString (s));
            _lumiLabel->SetNDC ();
            _lumiLabel->SetTextAlign (30);
            _lumiLabel->SetTextFont (42);
            _lumiLabel->SetTextSize (_tsize);
        }

    private: 
        std::vector<TH1D*> _hist;
        TH1D* _data;

        //MWL
        TLatex * _lumiLabel;     //PG label with the centre of mass energy and lumi info
        TString  _xLabel;
        TString  _units;
        bool     _breakdown;
        int      _mass;
	int      _alternativeOption;
        bool     _isHWWOverlaid;

};
