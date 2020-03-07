void rcc_draw_all_plots()
{
  gSystem->AddIncludePath("-I${OFFLINE_MAIN}/include");
  ifstream listfile("/phenix/spin2/pmontu/offline/analysis/pmontu/"
                    "relative_luminosity/macros/final_run_list.txt");
  int runnumber=0;

  const int nptbins=10;
  double pt_limits[nptbins+1] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};


  
  TString filename;
  TFile *histfile=NULL;


  int ngoodruns=0;
  const int maxruns=1000;
  const int maxptbins=10;
  double rawAsym[maxruns*maxptbins];
  double *asymByBin[maxptbins];
  TH1F *hWeightedAsym=new TH1F("hWeightedAsym","Weighted average asymmetry vs pt",nptbins,pt_limits);
  TH1F *hWeightSum=new TH1F("hWeightSum","Sum of Weights vs pt",nptbins,pt_limits);
  hWeightedAsym->Sumw2();
  
  for (int i=0;i<maxptbins;i++)
    asymByBin[i]=&(rawAsym[i*maxruns]);
  
  while (listfile.good()){
    listfile >> runnumber;
    // if (filename != 398149)
    //   continue;

    //must be matched to the output of rcc_calc_all.C:

    filename= Form("./asyms/%d.MPC.ALL.rcc.hist.root",runnumber);
    histfile=NULL;
    histfile=new TFile(filename,"READONLY");
    if ((histfile==NULL) || histfile->IsZombie() || !histfile->GetNkeys()) continue;

    
    TH1F *temp=0;
    temp=((TH1F*)(histfile->Get("hAllByPt")));


    //two things I can do:
    //1) plot vs pt and run
    //2) plot vs pt with sumw2.
    for (int i=0;i<nptbins;i++){
      double ptmid=(pt_limits[i]+pt_limits[i+1])/2;
      int sourcebin=temp->FindBin(ptmid);
      double asym=temp->GetBinContent(sourcebin);
      double err=temp->GetBinError(sourcebin);
      if (err<1e-6) err=asym; //if no error, assign 100% error for play.
      if (err<1e-6) continue; //if still no error, the asym was 0.  skip it.
      double err2=err*err;
      double w=1/err2;
      hWeightedAsym->Fill(ptmid,asym*w);
      hWeightSum->Fill(ptmid,w);
      //hErrSum->Fill(ptmid,err2/(asym*asym));
    //eventually, I should move this all to a simple ttree that puts entries in more granular way so I can divvy however I want.
    //Note that if input hist has Sumw2 set, Sumw2 is automatically called for this if not already set.
    }
  }

  
  hWeightedAsym->Divide(hWeightSum);
  
  for (int i=0;i<nptbins;i++){
    double ptmid=(pt_limits[i]+pt_limits[i+1])/2;
    int sourcebin=hWeightedAsym->FindBin(ptmid);
    double sumasym=hWeightedAsym->GetBinContent(sourcebin);
    double sumweight=hWeightSum->GetBinContent(sourcebin);
    if (sumweight==0) continue; //skip if we have no weight in this bin.
    double asym=sumasym/sumweight;
    //err_on_average=sqrt(1/sum of weights) as long as we pick weight=1/err2;
    double err=sqrt(1/hWeightSum->GetBinContent(sourcebin));
    //if (err>1) continue; //if large error,s omething went wrong.  skip it.

    hWeightedAsym->SetBinContent(sourcebin,asym);
    hWeightedAsym->SetBinError(sourcebin,err);
    }

hWeightedAsym->Draw();
TFile *outfile=TFile::Open("rcc_draw_all_plots.hist.root","RECREATE");
outfile->cd();
hWeightedAsym->Write();
outfile->Close();
return;


  
}
