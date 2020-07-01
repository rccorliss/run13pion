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
  
  TCut minclocks="clk>1e6"; //require a bunch to have at least 1e6 live clocks.
  TCut live70="clk/rclk>0.70";//require a bunch to be live at least 70% of the time.
  TCut minbbcrate="bbcwide>0.05";//reject bunches where the bbcwide trigger isn't acting right.
  TCut bbnslope="abs((bbncnt/bbcwidecnt-1.1)/(bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  TCut bbsslope="abs((bbscnt/bbcwidecnt-1.1)/(bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  
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
  if (1){
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
  if (1){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,1);
    c->cd(1);
    t->Draw("cross:zdcwidecnt","1","colz");
    t->SetMarkerColor(kRed);
    c->cd(2);
    t->Draw("cross:zdcwidecnt",minclocks && live70 && minbbcrate && bbsslope,"colz");
    t->SetMarkerColor(kBlack);   
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
