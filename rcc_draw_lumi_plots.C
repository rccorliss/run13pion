//code to read in bunch by bunch scalers and provide relative luminosity per bunch
//#include <TFileCollection>
#include <TCut.h>

void rcc_draw_lumi_plots(){
  //printf("tried to run rcc_gen_lumi directly instead of using the proper function calls.\n");

  //plot k_S and k_N vs BBC rate per-crossing, populating this plot with one entry for every crossing in every fill we are considering.

  //assume we're able to fill this vector:
  std::vector<int> run; // run numbers
  //here the prefix r_ denotes the vector contains run-by-run data
  std::vector<int> r_pattern; //spin pattern.
  //here the prefix b_ denotes the vector is 120x longer, and contains bunch-by-bunch data
  std::vector<double> b_bbc_n_cnt; //raw counts
  std::vector<double> b_bbc_s_cnt; //raw counts
  std::vector<double> b_bbc_ns_cnt; //raw counts
  std::vector<double> b_zdc_n_cnt; //raw counts
  std::vector<double> b_zdc_s_cnt; //raw counts
  std::vector<double> b_zdc_ns_cnt; //raw counts
  std::vector<double> b_clk;  // clock

  //even better, make a tchain that loads all the runs from pedro
  TChain *t=new TChain("t");


  //TFileCollection fc("files"); // The name is irrelevant
  //  fc.AddFromDirectory("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/*/rlstar.root");
  //can't do this this way:  t->Add("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/*/rlstar.root"); // can't glob in this fashion.

  /*
  const char *dirname="/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/";
  auto dir = gSystem->OpenDirectory(dirname);
  char *f;//[500];
  while (f = gSystem->GetDirEntry(dir)) { 
    if (!strcmp(f,".") || !strcmp(f,"..")) continue;
    t->Add(TString(dirname) + f + "/*.root");
  }
  gSystem->FreeDirectory(dir);
  printf("t has %d entries\n",t->GetEntries());
  */

  t->Add("unique_db.root");
   t->Add("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/unique_db.root");//same number of entries as previous, but syncs with gl1 and starscalers both.
  //one annoying catch:  all variables must be prefixed with "star_", "gl1_", or "gl1p_"
  // t->AddFileInfoList(fc.GetList());
  //delete fc;

  TCanvas *c;
  int nc=0;


  TCut rcc_cross_qa="star_bbcwide >0.05 && star_fill!=17443 && star_fill>17211"; //various cuts on basic bunch and rates
  TCut rcc_clip_loud_runs="1";
  int loud_runlist[]={391372,391581,391869,393906,398132,398026,398027,398019,398017,398018,398011,398010,398009,397176,389437,390177,391868,391870,398143,398030,398031,398029,398028,398020,398021,398014,398015,398013,398012,398005,398007,397177,397178};
  for (int i=0;i<33;i++){
    rcc_clip_loud_runs=rcc_clip_loud_runs+Form("star_run!=%d",loud_runlist[i]);
  }
    
  TCut minclocks="star_clk>1e6"; //require a bunch to have at least 1e6 live clocks.
  TCut allcuts=minclocks;
  TCut live70="star_clk/star_rclk>0.70";//require a bunch to be live at least 70% of the time.
  allcuts=allcuts && live70;
  TCut minbbcrate="star_bbcwide>0.05";//reject bunches where the bbcwide trigger isn't acting right.
  allcuts=allcuts && minbbcrate;
  TCut bbnslope="abs((star_bbncnt/star_bbcwidecnt-1.1)/(bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  allcuts=allcuts && bbnslope;
  TCut bbsslope="abs((star_bbscnt/star_bbcwidecnt-1.1)/(star_bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  allcuts=allcuts && bbsslope;
  TCut abortgap="star_cross<111";
  allcuts=allcuts && abortgap;
  TCut crossing1="star_cross!=1";
  allcuts=allcuts+crossing1;
  TCut stableratio="star_zdcwidecnt/star_bbcwidecnt<0.20";
  allcuts=allcuts && stableratio;
  TCut minzdccnt="star_zdcwidecnt/star_rclk>0.001";//reject bunches where the zdc trigger isn't acting right.
  allcuts=allcuts && minzdccnt;

  TCut rccCutSetA=rcc_cross_qa+rcc_clip_loud_runs+abortgap+crossing1;


  //show the SetB parameters we're cutting on.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw("star_bbcwiderun/gl1_bbcwidelive",rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(2)->SetLogy();
    t->Draw("star_bbc30run/gl1_bbc30live",rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(3)->SetLogy();
    t->Draw("star_zdcwiderun/gl1_zdcwidelive",rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
   c->cd(4)->SetLogy();
   t->Draw("star_bbcwiderun/gl1_bbcwidelive","star_clkrun"*(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(5)->SetLogy();
    t->Draw("star_bbc30run/gl1_bbc30live","star_clkrun"*(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(6)->SetLogy();
    t->Draw("star_zdcwiderun/gl1_zdcwidelive","star_clkrun"*(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //check if those values are correlated with run length:
  TCut bbcwidescalerratio="(star_bbcwiderun/gl1_bbcwidelive>0.99 && star_bbcwiderun/gl1_bbcwidelive<1.015)";
  TCut bbc30scalerratio="(star_bbc30run/gl1_bbc30live>0.99 && star_bbc30run/gl1_bbc30live<1.015)";
  TCut zdcwidescalerratio="(star_zdcwiderun/gl1_zdcwidelive>0.99 && star_zdcwiderun/gl1_zdcwidelive<1.015)";
  TCut starvsgl1ratios=bbcwidescalerratio+bbc30scalerratio+zdcwidescalerratio;
  TCut rccCutSetB=starvsgl1ratios;
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw(Form("%s>>bbcwidecut(2,-0.5,1.5)",bbcwidescalerratio.GetTitle()),"star_clkrun"*rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(2)->SetLogy();
    t->Draw(Form("%s>>bbc30cut(2,-0.5,1.5)",bbc30scalerratio.GetTitle()),rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(3)->SetLogy();
    t->Draw(Form("%s>>zdcwidecut(2,-0.5,1.5)",zdcwidescalerratio.GetTitle()),rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(4);//->SetLogy();
   t->Draw("star_bbcwiderun/gl1_bbcwidelive:star_clkrun",(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
   c->cd(5);//->SetLogy();
    t->Draw("star_bbc30run/gl1_bbc30live:star_clkrun",(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(6);//->SetLogy();
    t->Draw("star_zdcwiderun/gl1_zdcwidelive:star_clkrun",(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //some checks trying to understand what gl1p_ scalers are doing.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1);//->SetLogy();
   t->Draw("star_zdc30:log10(gl1p_zdc_narrow)",(rccCutSetB+rccCutSetA),"colz");
   c->cd(2);//->SetLogy();
     t->Draw("star_zdcwide:log10(gl1p_zdc_wide)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(3);//->SetLogy();
   t->Draw("star_bbc30:log10(gl1p_bbc_30)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(4);//->SetLogy();
   t->Draw("log10(star_zdc30/gl1p_zdc_narrow):log10(star_clkrun)",(rccCutSetB+rccCutSetA),"colz");
   c->cd(5);//->SetLogy();
    t->Draw("log10(star_zdcwide/gl1p_zdc_wide):log10(star_clkrun)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(6);//->SetLogy();
    t->Draw("log10(star_bbc30/gl1p_bbc_30):log10(star_clkrun)",(rccCutSetB+rccCutSetA),"colz");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //plot that shows the way the star and gl1p scalers track each other
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw("log10(star_zdc30cnt/gl1p_zdc_narrow)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
   c->cd(2)->SetLogy();
     t->Draw("log10(star_zdcwidecnt/gl1p_zdc_wide)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(3)->SetLogy();
   t->Draw("log10(star_bbc30cnt/gl1p_bbc_30)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(4);//->SetLogy();
   t->Draw("log10(star_zdc30cnt):log10(gl1p_zdc_narrow)",(rccCutSetB+rccCutSetA),"colz");
   c->cd(5);//->SetLogy();
    t->Draw("log10(star_zdcwidecnt):log10(gl1p_zdc_wide)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(6);//->SetLogy();
    t->Draw("log10(star_bbc30cnt):log10(gl1p_bbc_30)",(rccCutSetB+rccCutSetA),"colz");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  TCut starvsgl1pzdcwideratio="abs(star_zdcwidecnt/gl1p_zdc_wide-1)<0.001";
  TCut starvsgl1pzdc30ratio="abs(star_zdc30cnt/gl1p_zdc_narrow-1)<0.001";
  TCut starvsgl1pbbcratio="abs(star_bbc30cnt/gl1p_bbc_30-1)<0.001";
  TCut starvsgl1pratios=starvsgl1pzdcwideratio+starvsgl1pzdc30ratio+starvsgl1pbbcratio;
  rccCutSetB+=starvsgl1pratios;
  //plots showing the cut applied
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw("(star_zdc30cnt/gl1p_zdc_narrow)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
   c->cd(2)->SetLogy();
     t->Draw("(star_zdcwidecnt/gl1p_zdc_wide)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(3)->SetLogy();
   t->Draw("(star_bbc30cnt/gl1p_bbc_30)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(4);//->SetLogy();
   t->Draw("log10(star_zdc30cnt):log10(gl1p_zdc_narrow)",(rccCutSetB+rccCutSetA),"colz");
   c->cd(5);//->SetLogy();
    t->Draw("log10(star_zdcwidecnt):log10(gl1p_zdc_wide)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(6);//->SetLogy();
    t->Draw("log10(star_bbc30cnt):log10(gl1p_bbc_30)",(rccCutSetB+rccCutSetA),"colz");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //to ease redundant cuts, let's define a TEntryList and set it:
  t->Draw(">>elist", rccCutSetA+rccCutSetB, "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  t->SetEntryList(elist);

  TF1 *test=new TF1("testpoly","[0]+[1]*x+[2]*x^2",0.1,0.6);
  TGraph *g;
  TLatex tex;
  float texpos=0.75;
  float texshift=0;
  tex.SetTextAlign(12);
  tex.SetTextSize(0.03);
  int fiti=0;
  double pars[4][3];
  double parerrs[4][3];

  if (1){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(2,2);
    c->cd(1);//->SetLogy();
    t->Draw("star_zdscnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,0.25,Form("y=(%1.3E #pm %1.3E) + (%1.3E #pm %1.3E)x + (%1.3E #pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;

    c->cd(2);//->SetLogy();
    t->Draw("star_zdncnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,texpos,Form("y=(%1.3E#pm %1.3E) + (%1.3E#pm %1.3E)x + (%1.3E#pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;    c->cd(3);//->SetLogy();
    t->Draw("star_bbscnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,texpos,Form("y=(%1.3E#pm %1.3E) + (%1.3E#pm %1.3E)x + (%1.3E#pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;    c->cd(4);//->SetLogy();
    t->Draw("star_bbncnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,texpos,Form("y=(%1.3E#pm %1.3E) + (%1.3E#pm %1.3E)x + (%1.3E#pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;
    nc++;
  }


  
  return;
  printf("Caution!  Below this point are old Draw commands that assume we're running with the old starscaler-only files.  The names of variables will not match the combined db file.\n");
  return;
  
  //display rates vs bunch and fill
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:fill>>hnew1(120,17180,17620,120,-0.5,119.5)","log10(rclk)","colz");
    c->cd(2);
    t->Draw("cross:fill>>hnew2(120,17180,17620,120,-0.5,119.5)","log10(clk)","colz");
    c->SetLogz();
    c->cd(3);
    t->Draw("cross:fill>>hnew3(120,17180,17620,120,-0.5,119.5)","bbcwide","colz");
    c->SetLogz();
    c->cd(4);
    t->Draw("cross:fill","bbncnt/bbcwidecnt","colz");
    c->SetLogz();
    nc++;
  }

   //check some bunch crossing-related anomalies
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,900);
    c->Divide(2,3);
    c->cd(1);
    t->Draw("cross:fill","bbcwide *(cross>20 && cross<40)","colz");
    c->cd(2);
    t->Draw("cross:fill","bbcwide *(cross>60 && cross<80)","colz");
    c->cd(3);
    t->Draw("cross","bbcwide * (fill<17210 && cross<10)","colz");
    c->cd(4);
    t->Draw("cross","bbcwide * (fill>17210 && cross<10)","colz");
    c->cd(5);
    t->Draw("cross:fill","bbcwide *(cross>100 && fill>17430 && fill<17450)","colz");
    c->cd(6);
    t->Draw("(cross>110):bbcwide>>hnew(10,-0.05,0.95,2,-0.5,1.5)","1","colz");
   nc++;
  }

    //detailed check of those anomalies
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,900);
    c->Divide(2,3);
    c->cd(1);
    t->Draw("bbcwide","((cross==29 || cross==30 || cross==69 || cross==70) && fill>17380 && fill<17420)","colz");
    c->cd(2);
    t->Draw("bbcwide>0.1:run-386700","((cross==29 || cross==30 || cross==69 || cross==70) && fill>17380 && fill<17420)","colz");
    c->cd(3);
    t->Draw("cross:run-386700","bbcwide * (fill<17212 && cross<5)","colz");
    c->cd(4);
    t->Draw("fill:run-386700","(fill<17212 && cross==0)","colz");
    c->cd(5);
    t->Draw("cross:run-386700","bbcwide *(cross>100 && fill==17443)","colz");
    c->cd(6);
    t->Draw("bbcwide:run-386700","(cross>110 && fill>17435 && fill<17450)","colz");
   nc++;
  }

 //display rates vs bunch and fill after cuts designed to clean them.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:fill>>hnew1(120,17180,17620,120,-0.5,119.5)","log10(rclk)"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->cd(2);
    t->Draw("cross:fill>>hnew2(120,17180,17620,120,-0.5,119.5)","log10(clk)"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->SetLogz();
    c->cd(3);
    t->Draw("cross:fill>>hnew3(120,17180,17620,120,-0.5,119.5)","bbcwide"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->SetLogz();
    c->cd(4);
    t->Draw("cross:fill","bbncnt/bbcwidecnt"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->SetLogz();
    nc++;
  }

  //a marginally useful look at the abort gap spectrum after the rcc_cross_qa cuts
   if (0){
     c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("run-386700","cross>110"+rcc_cross_qa,"colz");
    c->cd(2);
    t->Draw("bbcwide","cross>110"+rcc_cross_qa,"colz");
    c->cd(3);
    t->Draw("run-386700","cross==29 || cross==30 || cross==69 || cross==70"+rcc_cross_qa,"colz");
    c->cd(4);
    t->Draw("bbcwide","cross==29 || cross==30 || cross==69 || cross==70"+rcc_cross_qa,"colz");
    nc++;
  }

   //a better look at the abort gap spectrum after the rcc_cross_qa cuts
   if (0){
     c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("run-386700","cross>110"+rcc_cross_qa,"colz");
    c->cd(2);
    t->Draw("bbcwide","cross>110"+rcc_cross_qa,"colz");
    c->cd(3);
    t->Draw("run-386700:fill>>hnew(350,17299.5,17649.5,100,2000,12000)","cross>110"+rcc_cross_qa,"colz");
    c->cd(4);
    t->Draw("run-386700:fill>>hnew2(350,17299.5,17649.5,100,2000,12000)","cross>110","colz");
    nc++;
  }
  
 //show how many live clocks we have in each bunch:
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("log10(clk)");
    t->SetLineColor(kRed);
    t->Draw("log10(clk)",minclocks,"same");
    t->SetLineColor(kBlack);
    nc++;
  }
   
  //show how the livetime goes with BBC rate
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("clk/rclk:bbcwide",minclocks);
    t->SetMarkerColor(kRed);
    t->Draw("clk/rclk:bbcwide",minclocks && live70,"same");
    t->SetMarkerColor(kBlack);
    nc++;
  }

 //show the uncorrected(?) BBC rate for the survivors:
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("bbcwide",minclocks && live70);
    t->SetLineColor(kRed);
    t->Draw("bbcwide",minclocks&& live70 &&minbbcrate,"same");
    t->SetLineColor(kBlack);
    nc++;
  }
  
  //show how bbcs singles to doubles goes with bbc rate, with livetime, rate, and runlength cuts
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("bbscnt/bbcwidecnt:bbcwide",minclocks && live70 && minbbcrate,"colz");
    nc++;
  }

  //assuming these all intercept at y=1.1, x=0.7,
  //we construct their slope
  //show how ~k slope is distributed after cuts
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("(bbscnt/bbcwidecnt-1.1)/(bbcwide-0.7)",minclocks && live70 && minbbcrate);
    t->SetLineColor(kRed);
    t->Draw("(bbscnt/bbcwidecnt-1.1)/(bbcwide-0.7)",minclocks&& live70 &&minbbcrate &&bbsslope,"same");
    t->SetLineColor(kBlack);
    nc++;
  }
  
   //show how bbcs singles to doubles goes with bbc rate, with all cuts
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("bbscnt/bbcwidecnt:bbcwide",minclocks && live70 && minbbcrate);
    t->SetMarkerColor(kRed);
    t->Draw("bbscnt/bbcwidecnt:bbcwide",minclocks && live70 && minbbcrate && bbsslope,"same");
    t->SetMarkerColor(kBlack);   
    nc++;
  }

  //show how bbc north singles to doubles goes with bbc rate, with all cuts -- do the south cuts apply to the north?
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("bbncnt/bbcwidecnt:bbcwide",minclocks && live70 && minbbcrate && bbsslope);
    t->SetMarkerColor(kRed);
    t->Draw("bbncnt/bbcwidecnt:bbcwide",minclocks && live70 && minbbcrate && bbsslope && bbnslope,"same");
    t->SetMarkerColor(kBlack);   
    nc++;
  }

  //show how zdc south singles to doubles goes with bbc rate, with all cuts -- do the bbc cuts apply to the zdc?
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->cd(1);
    t->Draw("zdscnt:zdcwidecnt",minclocks && live70 && minbbcrate);
    t->SetMarkerColor(kRed);
    t->Draw("zdscnt:zdcwidecnt",minclocks && live70 && minbbcrate && bbsslope,"same");
    t->SetMarkerColor(kBlack);   
    nc++;
  }

  //show zdc vs bunch xing, to make the abort gaps visible.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:zdcwidecnt","1","colz");
    t->SetMarkerColor(kRed);
    c->cd(2);
    t->Draw("cross:zdcwidecnt",minclocks && live70 && minbbcrate && bbsslope,"colz");
    t->SetMarkerColor(kBlack);   
    c->cd(3);
    t->Draw("cross:zdcwidecnt",abortgap,"colz");
    t->SetMarkerColor(kBlack);   
    c->cd(4);
    t->Draw("cross:zdcwidecnt",abortgap&&minclocks && live70 && minbbcrate && bbsslope,"colz");
    t->SetMarkerColor(kBlack);   
    nc++;
  }

    //after all cuts, zdc to bbc ratio, to look for excursions
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(3,3);

    //ratio vs run
    c->cd(1);
    t->Draw("zdcwidecnt/bbcwidecnt:run","1","colz");
    t->SetMarkerColor(kRed);
    c->cd(2);
    t->Draw("zdcwidecnt/bbcwidecnt:run","zdcwidecnt/bbcwidecnt<0.45","colz");
    t->SetMarkerColor(kRed);
    c->cd(3);
    t->Draw("zdcwidecnt/bbcwidecnt:run",allcuts,"colz");
    t->SetMarkerColor(kBlack);

    //ratios vs one part
    c->cd(4);
    t->Draw("zdcwidecnt/bbcwidecnt:zdcwidecnt",allcuts,"colz");
     c->cd(5);
     t->Draw("zdcwidecnt:bbcwidecnt",allcuts,"colz");
     c->cd(6);
    t->Draw("zdcwidecnt/bbcwidecnt:bbcwidecnt",allcuts,"colz");  
  //ratios vs one part, normalized by live clocks
    c->cd(7);
    t->Draw("zdcwidecnt/bbcwidecnt:zdcwidecnt/rclk",allcuts,"colz");
     c->cd(8);
     t->Draw("zdcwidecnt/rclk:bbcwidecnt/rclk",allcuts,"colz");
     c->cd(9);
    t->Draw("zdcwidecnt/bbcwidecnt:bbcwidecnt/rclk",allcuts,"colz");  
    nc++;
  }

  //after all cuts, inspect the corrected ratios and see how smooth they're getting...
  if (1){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(1,3);

  //ratios vs one part, normalized by live clocks
    c->cd(1);
    t->Draw("zdcwidecnt/bbcwidecnt:zdcwidecnt/clk",allcuts);
     c->cd(2);
     t->Draw("zdcwidecnt/rclk:bbcwidecnt/clk",allcuts);
     c->cd(3);
    t->Draw("zdcwidecnt/bbcwidecnt:bbcwidecnt/clk",allcuts);  
    nc++;
  }
  

  return;
   //show how bbcn singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("bbncnt/bbcwidecnt:bbcwide","bbcwide>0.01");
  nc++;

    //show how zdcs singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("zdscnt/zdcwidecnt:bbcwide","bbcwide>0.01","colz");
  nc++;
  
   //show how zdcn singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("zdncnt/zdcwidecnt:bbcwide","bbcwide>0.01","colz");
  nc++;
  
  /*
  
  
  const int nSamples=120*500;
  double kN[nSamples];//ratio of north singles to coin for a particular xing of a particular run
  double kS[nSamples];//ratio of south singles to coin for a particular xing of a particular run
  double pkl_A[nSamples];//inferred fraction of clocks that contain a trigger (we flip between these as we iterate)
  double pkl_B[nSamples];//inferred fraction of clocks that contain a trigger (we flip between these as we iterate)
  double *pkl=pkl_A; //start with pkl_A;
  double *pklNext=pkl_B; //and fill in pkl_B;
  //these k_'s are fixed truth values -- the singles scaler and the doubles scaler.  It is the x-axis, the true BBC rate per crossing, that we are adjusting.  In iteration zero we use Pkl=BBCLL1/CLOCK.
  */
  
  //For the next iteration, we fit k(Pkl) with a polynomial and extrapolate to zero,( the constant term) which gives us k_ prime.

  
  //We use these new values of k_ , and Pkl to solve 1.26 (pasted above this line) to solve for Pkl'=mu*eps_NS .  ALthough they don't show the equation in the form we would like, that's what is being done.  Pkl, is put into the left hand side and the equation is numerically solved for Pkl'=mu*eps_NS.
  
  //Once this is done, Pkl' is used to plot the same k_ values as before, where we now fit (k(Pkl') and iterate to get Pkl'', and so on, until this converges.  They state that this always converges by ~ 4 iterations.


  
  return;
}
