//code to read in bunch by bunch scalers and provide relative luminosity per bunch
#include <TFileCollection>

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
  TFileCollection fc("files"); // The name is irrelevant
  //  fc.AddFromDirectory("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/*/rlstar.root");
  //can't do this this way:  t->Add("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/*/rlstar.root"); // can't glob in this fashion.

  const char *dirname="/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/";
  auto dir = gSystem->OpenDirectory(dirname);
  char *f;//[500];
  while (f = gSystem->GetDirEntry(dir)) { 
    if (!strcmp(f,".") || !strcmp(f,"..")) continue;
    t->Add(TString(dirname) + f + "/*.root");
  }
  gSystem->FreeDirectory(dir);
  printf("t has %d entries\n",t->GetEntries());
  // t->AddFileInfoList(fc.GetList());
  //delete fc;

  TCanvas *c;
  int nc=0;


  TCut rcc_cross_qa="bbcwide >0.05 && fill!=17443 && fill>17211"; //various cuts on basic bunch and rates
    
  TCut minclocks="clk>1e6"; //require a bunch to have at least 1e6 live clocks.
  TCut allcuts=minclocks;
  TCut live70="clk/rclk>0.70";//require a bunch to be live at least 70% of the time.
  allcuts=allcuts && live70;
  TCut minbbcrate="bbcwide>0.05";//reject bunches where the bbcwide trigger isn't acting right.
  allcuts=allcuts && minbbcrate;
  TCut bbnslope="abs((bbncnt/bbcwidecnt-1.1)/(bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  allcuts=allcuts && bbnslope;
  TCut bbsslope="abs((bbscnt/bbcwidecnt-1.1)/(bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  allcuts=allcuts && bbsslope;
  TCut abortgap="cross<111";
  allcuts=allcuts && abortgap;
  TCut stableratio="zdcwidecnt/bbcwidecnt<0.20";
  allcuts=allcuts && stableratio;
  TCut minzdccnt="zdcwidecnt/rclk>0.001";//reject bunches where the zdc trigger isn't acting right.
  allcuts=allcuts && minzdccnt;

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
  if (1){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:fill>>hnew1(120,17180,17620,120,-0.5,119.5)","log10(rclk)"+rcc_cross_qa,"colz");
    c->cd(2);
    t->Draw("cross:fill>>hnew2(120,17180,17620,120,-0.5,119.5)","log10(clk)"+rcc_cross_qa,"colz");
    c->SetLogz();
    c->cd(3);
    t->Draw("cross:fill>>hnew3(120,17180,17620,120,-0.5,119.5)","bbcwide"+rcc_cross_qa,"colz");
    c->SetLogz();
    c->cd(4);
    t->Draw("cross:fill","bbncnt/bbcwidecnt","colz");
    c->SetLogz();
    nc++;
  }
  
  return;
  
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
