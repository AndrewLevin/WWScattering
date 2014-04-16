// poor man macro to remind how to produce the plots

// remove data while blinded
.x finalPlotVBS.C+(1,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mjj.root","wwss_mjj",0,19.4,1)
.x finalPlotVBS.C+(1,1,"m_{ll} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_mll.root","wwss_mll",0,19.4,1)
.x finalPlotVBS.C+(1,1,"p_{T}^{l max} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_ptlmax.root","wwss_ptlmax",0,19.4,1)

// change wrong-sign legend
.x finalPlotVBS.C+(1,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/toptag.root","toptag",0,19.4)

// standard
.x finalPlotVBS.C+(1,1,"m_{jj} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wz.root","wz",0,19.4)
.x finalPlotVBS.C+(1,1,"N_{jets}","","/data/smurf/ceballos/distributions/wwss/wwss_btag_njets.root","wwss_btag_njets",0,19.4)
.x finalPlotVBS.C+(1,10,"p_{T}^{j max} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetpt1.root","wwss_presel_jetpt1",0,19.4)
.x finalPlotVBS.C+(1,10,"p_{T}^{j min} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetpt2.root","wwss_presel_jetpt2",0,19.4)
.x finalPlotVBS.C+(1,5,"|#eta^{j}|^{max}","","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetetamax.root","wwss_presel_jetetamax",0,19.4)
.x finalPlotVBS.C+(1,5,"|#eta^{j}|^{min}","","/data/smurf/ceballos/distributions/wwss/wwss_presel_jetetamin.root","wwss_presel_jetetamin",0,19.4)
.x finalPlotVBS.C+(1,1,"m_{ll} (GeV)","GeV","/data/smurf/ceballos/distributions/wwss/wwss_presel_mll.root","wwss_presel_mll",0,19.4)
.x finalPlotVBS.C+(1,1,"|#Delta #eta_{jj}|","","/data/smurf/ceballos/distributions/wwss/wwss_presel_detajj.root","wwss_presel_detajj",0,19.4)

// N-1
.x finalPlotVBS.C+(1,1,"|#Delta #eta_{jj}|","","/data/smurf/ceballos/distributions/wwss/wwss_np1_detajj.root","wwss_np1_detajj",0,19.4)
