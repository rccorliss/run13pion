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
  
  //c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  //c->cd(1);
  //t->Draw("rclk/clk:bbcwide","1","colz");
  //show how the livetime goes with BBC rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("rclk/clk:bbcwide","rclk/clk<2","colz");
  nc++;


  //show how bbcs singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("bbscnt/bbcwidecnt:bbcwide","rclk/clk<2","colz");
  nc++;
  
   //show how bbcn singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("bbncnt/bbcwidecnt:bbcwide","rclk/clk<2","colz");
  nc++;

    //show how zdcs singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("zdscnt/zdcwidecnt:bbcwide","rclk/clk<2","colz");
  nc++;
  
   //show how zdcn singles to doubles goes with bbc rate
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
  c->cd(1);
  t->Draw("zdncnt/zdcwidecnt:bbcwide","rclk/clk<2","colz");
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
