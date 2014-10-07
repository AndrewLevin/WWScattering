# poor man macro to remind how to produce the plots

root -l -q -b finalPlotVBS.C+'(-1,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mjj.root","wwss_mjj",0,19.4)'
root -l -q -b finalPlotVBS.C+'(-1,1,"m_{ll} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mll.root","wwss_mll",0,19.4)'
root -l -q -b finalPlotVBS.C+'(-1,1,"p_{T}^{l max} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_ptlmax.root","wwss_ptlmax",0,19.4)'

# standard
root -l -q -b finalPlotVBS.C+'(1,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wz.root","wz",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,1,"N_{jets}","","/data/smurf/ceballos/distributions/wwss/wwss_btag_njets.root","wwss_btag_njets",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,10,"p_{T}^{j max} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetpt1.root","wwss_presel_jetpt1",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,10,"p_{T}^{j min} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetpt2.root","wwss_presel_jetpt2",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,5,"|#eta^{j}|^{max}","","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetetamax.root","wwss_presel_jetetamax",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,5,"|#eta^{j}|^{min}","","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetetamin.root","wwss_presel_jetetamin",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,1,"m_{ll} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_presel_mll.root","wwss_presel_mll",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,1,"|#Delta #eta_{jj}|","","/data/smurf/ceballos/distributions/wwss/wwss_presel_detajj.root","wwss_presel_detajj",0,19.4)'

# N-1
root -l -q -b finalPlotVBS.C+'(1,1,"|#Delta #eta_{jj}|","","/data/smurf/ceballos/distributions/wwss/wwss_np1_detajj.root","wwss_np1_detajj",0,19.4)'

# change h legend
root -l -q -b finalPlotVBS.C+'(1,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mjj_hpp200.root","wwss_mjj_hpp200",0,19.4)'
root -l -q -b finalPlotVBS.C+'(1,1,"m_{ll} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mll_hpp200.root","wwss_mll_hpp200",0,19.4)'
root -l -q -b finalPlotVBS.C+'(2,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mjj_hpp800.root","wwss_mjj_hpp800",0,19.4)'
root -l -q -b finalPlotVBS.C+'(2,1,"m_{ll} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mll_hpp800.root","wwss_mll_hpp800",0,19.4)'

# change wrong-sign legend
root -l -q -b finalPlotVBS.C+'(3,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/toptag.root","toptag",0,19.4)'

root -l -q -b finalPlotVBS.C+'(3,8,"p_{T}^{l max} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/toptag_ptlmax.root","toptag_ptlmax",0,19.4)'
root -l -q -b finalPlotVBS.C+'(3,8,"p_{T}^{l min} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/toptag_ptlmin.root","toptag_ptlmin",0,19.4)'

# Limit plot
root -l -b -q PlotLimitForVBS.C+'("inputs/ana_hpp_samesign.txt","ana_hpp_samesign","19.4 fb^{-1} (8 TeV)",160,800,0,0,"VBF H^{#pm#pm} #rightarrow W^{#pm}W^{#pm}",1,5,"pdf")';
root -l -b -q PlotLimitForVBS.C+'("inputs/ana_hpp_samesign.txt","ana_hpp_samesign","19.4 fb^{-1} (8 TeV)",160,800,0,0,"VBF H^{#pm#pm} #rightarrow W^{#pm}W^{#pm}",1,5,"png")';
root -l -b -q PlotLimitForVBS.C+'("inputs/ana_hpp_samesign.txt","ana_hpp_samesign","19.4 fb^{-1} (8 TeV)",160,800,0,0,"VBF H^{#pm#pm} #rightarrow W^{#pm}W^{#pm}",1,5,"eps")';
