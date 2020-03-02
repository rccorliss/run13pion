//This macro runs rcc_gen_yield over all runs in the listfile.
// for future versions:  it really ought to define inputs and outputs directories.
// This produces the hist files that are consumed by rcc_calc_all.C


void Run_rcc_gen_yield()
{
  int minbias = 0;
  gSystem->AddIncludePath("-I${OFFLINE_MAIN}/include");
  gSystem->Load("libmpc.so");
  gSystem->Load("libuspin.so");
  gROOT->ProcessLine(".L rcc_gen_all.C+");
  //ifstream listfile("Run11_singlefile.txt");
  ifstream listfile("/phenix/spin2/pmontu/offline/analysis/pmontu/"
                    "relative_luminosity/macros/final_run_list.txt");
  int filename;
  //GenALL(runnum,minbias);
  while (listfile.good()){
    listfile >> filename;
    // if (filename != 398149)
    //   continue;
    if (listfile.good()){
      rcc_gen_yield(filename,minbias);
    }
  }
}
