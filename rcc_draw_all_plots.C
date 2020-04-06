void rcc_draw_all_plots()
{
  gSystem->AddIncludePath("-I${OFFLINE_MAIN}/include");
  ifstream listfile("/phenix/spin2/pmontu/offline/analysis/pmontu/"
                    "relative_luminosity/macros/final_run_list.txt");
  int runnumber=0;

  const int nptbins=10;
  double pt_limits[nptbins+1] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};
  const int nzdcbins=3;
  double zdc_limits[nzdcbins+1]={0,2e9,4e9,1e10};
  const int nratiobins=3;
  double ratio_limits[nratiobins+1]={0.8,0.95,1.03,1.2};
  const int nspinpats=10;
  const int spinpat_lower_limit=20;
  float spinpat_bound[]={spinpat_lower_limit-0.5,spinpat_lower_limit+nspinpats-0.5};

  
  TString filename;
  TFile *histfile=NULL;


  int ngoodruns=0;
  const int maxruns=1000;
  const int maxptbins=10;
  double rawAsym[maxruns*maxptbins];
  double *asymByBin[maxptbins];

  //some sanity check plots:
  TH2F *hPolarizationRMS=new TH2F("hPolarizationRMS","RMS of bpol and ypol per run (checks if they vary over the run);bpol RMS;ypol RMS",20,0,0.01,20,0,0.01);
  TH2F *hPolarizationMean=new TH2F("hPolarizationMean","Mean bpol and ypol for each run;bpol;ypol",100,0.2,1,100,0.2,1);


  //some other monitor plots:
  TH1F *hZdcNarrowSum=new TH1F("hZdcNarrowSum","ZdcNarrow sum for all runs",200,0,1e10);
  TH1F *hZdcNarrowRatio=new TH1F("hZdcNarrowRatio","ZdcNarrow Ratio (Unlike over Like spin) for all runs",200,0.5,1.5);

  TTree *mtree=new TTree("mTree","Run-wise Metadata Tree");
  float m_zdcsum,m_zdcratio;
  int m_spinpat,m_run;
  float m_bpol,m_ypol;
  int m_i;//increments by one every tiem we fill.  Unique run ID, but not run number.
  mtree->Branch("zdcsum",&m_zdcsum);
  mtree->Branch("zdcratio",&m_zdcratio);
  mtree->Branch("spinpat",&m_spinpat);
  mtree->Branch("run",&m_run);
  mtree->Branch("bpol",&m_bpol);
  mtree->Branch("ypol",&m_ypol);
  mtree->Branch("index",&m_i);
  

  
  TH1F *hWeightedAsym=new TH1F("hWeightedAsym","Weighted average asymmetry vs pt;pt [GeV];asym",nptbins,pt_limits);
  TH1F *hWeightSum=new TH1F("hWeightSum","Sum of Weights vs pt;pt [GeV];weight",nptbins,pt_limits);
  hWeightedAsym->Sumw2();

  //accumulators for asymmetry divided into zdc bins
  TH2F *hWeightedAsymByLumi=new TH2F("hWeightedAsymByLumi","Weighted average asymmetry vs pt",nzdcbins,zdc_limits,nptbins,pt_limits);
  TH2F *hWeightSumByLumi=new TH2F("hWeightSumByLumi","Sum of Weights vs pt",nratiobins,ratio_limits,nptbins,pt_limits);

   //accumulators for asymmetry divided into zdc bins
  TH2F *hWeightedAsymByRatio=new TH2F("hWeightedAsymByRatio","Weighted average asymmetry vs pt and zdc ratio",nzdcbins,zdc_limits,nptbins,pt_limits);
  TH2F *hWeightSumByRatio=new TH2F("hWeightSumByRatio","Sum of Weights vs pt and zdc ratio",nratiobins,ratio_limits,nptbins,pt_limits);

  //accumulators for asymmetry divided into spin patterns bins
  TH2F *hWeightedAsymBySpinpat=new TH2F("hWeightedAsymBySpinpat","Weighted average asymmetry vs pt;spinpat;pt [GeV]",nspinpats,spinpat_bound[0],spinpat_bound[1],nptbins,pt_limits);
  TH2F *hWeightSumBySpinpat=new TH2F("hWeightSumBySpinpat","Sum of Weights vs pt",nspinpats,spinpat_bound[0],spinpat_bound[1],nptbins,pt_limits);

  m_i=0;
  for (int i=0;i<maxptbins;i++)
    asymByBin[i]=&(rawAsym[i*maxruns]);


  //loop over all runs.
  while (listfile.good()){
    listfile >> runnumber;
    // if (filename != 398149)
    //   continue;

    //must be matched to the output of rcc_calc_all.C:

    filename= Form("./asyms/%d.MPC.ALL.rcc.hist.root",runnumber);
    histfile=NULL;
    histfile=new TFile(filename,"READONLY");
    if ((histfile==NULL) || histfile->IsZombie() || !histfile->GetNkeys()) continue;

    TH1F* htags=((TH1F*)(histfile->Get("hTags")));
    int spinpat=htags->GetBinContent(htags->GetBin(htags->GetXaxis()->FindBin("spinpattern")));
  
    TH1F *temp=0;
    temp=((TH1F*)(histfile->Get("hAllByPt")));
    
    TH1F *hLumi=((TH1F*)(histfile->Get("hTotalLumi")));
    float lumi=hLumi->Integral();//does this work?  Tired.
    float lumiLike=hLumi->GetBinContent(2);
    float lumiUnlike=hLumi->GetBinContent(1);
    hZdcNarrowSum->Fill(lumi);
    hZdcNarrowRatio->Fill(lumiUnlike/lumiLike);
 
    
    TH2F *hPol=0;
    for (int i=0;i<2;i++){
      char b=(i?'P':'N');
      for (int j=0;j<2;j++){
	char y=(j?'P':'N');
	hPol=((TH2F*)(histfile->Get(Form("hPolarizationBySpin%c%c",b,y))));
	hPolarizationRMS->Fill(hPol->GetRMS(1),hPol->GetRMS(2));
	hPolarizationMean->Fill(hPol->GetMean(1),hPol->GetMean(2));
      }
    }

    m_zdcsum=lumi;
    m_zdcratio=lumiUnlike/lumiLike;
    m_spinpat=spinpat;
    m_run=runnumber;
    m_bpol=hPol->GetMean(1);
    m_ypol=hPol->GetMean(2);
    mtree->Fill();
    m_i++;
    
    //two things I can do:
    //1) plot vs pt and runnumber
    
    //1.1) plot vs pt for 'low' and 'high' luminosity:
    
      
    //2) plot vs pt with sumw2.
    for (int i=0;i<nptbins;i++){
      double ptmid=(pt_limits[i]+pt_limits[i+1])/2;
      int sourcebin=temp->FindBin(ptmid);
      float asym=temp->GetBinContent(sourcebin);
      float  err=temp->GetBinError(sourcebin);
      if (err<1e-9) err=asym; //if no error, assign 100% error for play.
      if (err<1e-9) {
	err=1e-9; //originally I skipped if this happened, but that might've artificially inflated if my error was very small.
      }
      double err2=err*err;
      double w=1/err2;
      hWeightedAsym->Fill(ptmid,asym*w);
      hWeightSum->Fill(ptmid,w); //so every bin gets the same weight, for later division.

      //add this to a weighted sum for the zdc bin of choice:
      hWeightedAsymByLumi->Fill(lumi,ptmid,asym*w);
      hWeightSumByLumi->Fill(lumi,ptmid,w);
      hWeightedAsymByRatio->Fill(m_zdcratio,ptmid,asym*w);
      hWeightSumByRatio->Fill(m_zdcratio,ptmid,w);
      hWeightedAsymBySpinpat->Fill(spinpat,ptmid,asym*w);
      hWeightSumBySpinpat->Fill(spinpat,ptmid,w);
      
      //hErrSum->Fill(ptmid,err2/(asym*asym));
    //eventually, I should move this all to a simple ttree that puts entries in more granular way so I can divvy however I want.
    //Note that if input hist has Sumw2 set, Sumw2 is automatically called for this if not already set.
    }
  }

  //at this point we have two histograms:
  //weightedAsym has the sum of asym/asym_err^2 over all runs
  //weightSum has the sum of 1/asym_err^2 over all runs, and for convenience has the same bin structure
  //the average Asym per bin would be the ratio of these two, but that might not respect the error propagation.
  //as long as our weighting was weight=1/err2, the err on the average is sqrt(1/sum of weights).
  
  hWeightedAsym->Divide(hWeightSum);
  hWeightedAsymByLumi->Divide(hWeightSumByLumi);
  hWeightedAsymBySpinpat->Divide(hWeightSumBySpinpat);
  hWeightedAsymByRatio->Divide(hWeightSumByRatio);

  //fix the errorbars:
  for (int i=0;i<nptbins;i++){
    float sumweight=hWeightSum->GetBinContent(i+1);
    if (sumweight==0) continue; //skip if we have no weight in this bin.
    //err_on_average=sqrt(1/sum of weights) as long as we pick weight=1/err2;
    float err=sqrt(1/sumweight);
    hWeightedAsym->SetBinError(i+1,err);
  }

  //fix the error bars on the 2D plots and extract their 1D components, so they can be more easily visualized:
  TH1F *hFinalAsymByLumi[nzdcbins];
  for (int j=0;j<nzdcbins;j++){
    hFinalAsymByLumi[j]=new TH1F(Form("hFinalAsymByLumi%d",j),
				 Form("Weighted asymmetry vs pt for runs with %f<zdc_narrow<=%f;pt [GeV];asym",
				      zdc_limits[j],zdc_limits[j+1]),
				 nptbins,pt_limits);
    for (int i=0;i<nptbins;i++){
      float sumweight=hWeightSumByLumi->GetBinContent(j+1,i+1);
      float asym=hWeightedAsymByLumi->GetBinContent(j+1,i+1);
      if (sumweight==0) continue; //skip if we have no weight in this bin.
      //err_on_average=sqrt(1/sum of weights) as long as we pick weight=1/err2;
      float err=sqrt(1/sumweight);
      hFinalAsymByLumi[j]->SetBinContent(i+1,asym);
      hFinalAsymByLumi[j]->SetBinError(i+1,err);
    }
  }

  TH1F *hFinalAsymByRatio[nratiobins];
  for (int j=0;j<nratiobins;j++){
    hFinalAsymByRatio[j]=new TH1F(Form("hFinalAsymByRatio%d",j),
				 Form("Weighted asymmetry vs pt for runs with %f<unlike/like<=%f;pt [GeV];asym",
				      ratio_limits[j],ratio_limits[j+1]),
				 nptbins,pt_limits);
    for (int i=0;i<nptbins;i++){
      float sumweight=hWeightSumByRatio->GetBinContent(j+1,i+1);
      float asym=hWeightedAsymByRatio->GetBinContent(j+1,i+1);
      if (sumweight==0) continue; //skip if we have no weight in this bin.
      //err_on_average=sqrt(1/sum of weights) as long as we pick weight=1/err2;
      float err=sqrt(1/sumweight);
      hFinalAsymByRatio[j]->SetBinContent(i+1,asym);
      hFinalAsymByRatio[j]->SetBinError(i+1,err);
    }
  }

  TH1F *hFinalAsymBySpinpat[nspinpats];
  for (int j=0;j<nspinpats;j++){
    hFinalAsymBySpinpat[j]=new TH1F(Form("hFinalAsymBySpinpat%d",j),
				 Form("Weighted asymmetry vs pt for runs with spinpat=%d;pt [GeV];asym",j+spinpat_lower_limit),
				 nptbins,pt_limits);
    for (int i=0;i<nptbins;i++){
      float sumweight=hWeightSumBySpinpat->GetBinContent(j+1,i+1);
      float asym=hWeightedAsymBySpinpat->GetBinContent(j+1,i+1);
      if (sumweight==0) continue; //skip if we have no weight in this bin.
      //err_on_average=sqrt(1/sum of weights) as long as we pick weight=1/err2;
      float err=sqrt(1/sumweight);
      hFinalAsymBySpinpat[j]->SetBinContent(i+1,asym);
      hFinalAsymBySpinpat[j]->SetBinError(i+1,err);
    }
  }
    
  
  hWeightedAsym->Draw();
  TFile *outfile=TFile::Open("rcc_draw_all_plots.hist.root","RECREATE");
  outfile->cd();
  hWeightedAsym->Write();
  hWeightedAsymByLumi->Write();
  hWeightedAsymBySpinpat->Write();
  hPolarizationRMS->Write();
  hPolarizationMean->Write();
  hZdcNarrowSum->Write();
  hZdcNarrowRatio->Write();

  //only draw the ones that had some entries.
  for (int j=0;j<nspinpats;j++){
    //if (hWeightSumBySpinpat[j]->GetBinContent(1,1)>0)
      hFinalAsymBySpinpat[j]->Write();
  }
  for (int j=0;j<nzdcbins;j++){
    //if (hWeightSumByLumi[j]->GetBinContent(1,1)>0)
      hFinalAsymByLumi[j]->Write();
  }
  for (int j=0;j<nratiobins;j++){
      hFinalAsymByRatio[j]->Write();
  }
  mtree->Write();
  outfile->Close();
  return;


  
  }
