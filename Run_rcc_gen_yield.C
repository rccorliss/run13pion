//This macro runs rcc_gen_yield over all runs in the listfile.
// for future versions:  it really ought to define inputs and outputs directories.
// This produces the hist files that are consumed by rcc_calc_all.C
// In this version, it also calls rcc_calc_all to produce the final asymmetries.


void Run_rcc_gen_yield()
{
  int minbias = 0;
  gSystem->AddIncludePath("-I${OFFLINE_MAIN}/include");
  gSystem->Load("libmpc.so");
  gSystem->Load("libuspin.so");
  gROOT->ProcessLine(".L rcc_gen_yield.C+");
  gROOT->ProcessLine(".L rcc_calc_all.C+");
  //ifstream listfile("Run11_singlefile.txt");
  //ifstream listfile("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/macros/final_run_list.txt");
  ifstream listfile("rcc_runlist_jan2021.txt");
  int filename;
  //GenALL(runnum,minbias);
  while (listfile.good()){
    listfile >> filename;
    printf("loading %d\n", filename);
    //uncomment the following line to regenerate the yield files:
    //rcc_gen_yield(filename,"/phenix/spin/phnxsp01/rosscorliss/taxi/Run13pp510MPC/15944/data/","./yields/");
    rcc_gen_yield(filename,"/phenix/spin/phnxsp01/rosscorliss/taxi/Run13pp510MPC/16882/data/","./yields2021/");
    //rcc_calc_all(filename,"./yields2021/","./asyms/");
      
   
  }

  
}
