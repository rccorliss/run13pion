

TTree *uLumi, *uLumiXL;
TTree *uBunch;
const int nPats=16;
const int spinpat[]={1,2,3,4,5,6,7,8,21,22,23,24,25,26,27,28};

  const int nptbins=10;
const double pt_limits[] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};//one more than nptbins, to cover lower and upper bounds

TCanvas *c; int nc=0; //canvas and a running counter of how many canvases I've made
TGraph *gt;
TF1 *ft;
TLegend *leg;

void DrawPionALL();
TGraphErrors * GeneratePionAllTGraph(TString uBunchCut, TString uLumiCut, int arm=0);
TGraphErrors * GenerateBunchwisePionAllTGraph(TString uBunchCut, TString uLumiCut, int arm=0);
void GenerateRelativeLumiZDC(double *rellumi,double *rellumi_err, TString uBunchCut, TString uLumiCut);
void CompareEarlyAndLateSpinConfigs();



void DrawPionALLwithArms();
void DrawPionALLwithEvenVsOdd();
void DrawPionALLwithPatternSets();
void DrawPionALLwithPatternSetsAndEvenVsOdd();
void DrawPionALLwithEarlyVsLate();

void DrawRelativeLumiALL();
void DrawRelativeLumiZDC();


void CalcAsymAndErr(double *asym, double *asym_err,
		    double bpol, double bpol_err,
		    double ypol, double ypol_err,
		    double unlike_signal, double unlike_signal_err,
		    double like_signal, double like_signal_err,
		    double unlike_lumi, double unlike_lumi_err,
		    double like_lumi, double like_lumi_err);



void rcc_draw_lumi_asym(){
  //a little booster program to plot some stuff generated from lots of calcs in draw_lumi_plots that don't need to be repeaated every time:
  TFile *uLumiFile=TFile::Open("uLumi.ttree.root","READ");
  TFile *uLumiXLfile=TFile::Open("uLumiXL.ttree.root","READ");
  TFile *uBunchFile=TFile::Open("uBunch.ttree.root","READ");

  
  uLumi=(TTree*)uLumiFile->Get("uLumi");
  uLumiXL=(TTree*)uLumiXLfile->Get("uLumiXL");
  int u_fill; uLumi->SetBranchAddress("fill",&u_fill);
  int u_run; uLumi->SetBranchAddress("run",&u_run);
  int u_pat; uLumi->SetBranchAddress("pat",&u_pat);
  double u_bpol; uLumi->SetBranchAddress("bpol",&u_bpol);
  double u_ypol; uLumi->SetBranchAddress("ypol",&u_ypol);
  double u_bpol_err; uLumi->SetBranchAddress("bpol_err",&u_bpol_err);
  double u_ypol_err; uLumi->SetBranchAddress("ypol_err",&u_ypol_err);

  double u_likemuz0; uLumi->SetBranchAddress("likemuz0",&u_likemuz0);
  double u_likemuz1; uLumi->SetBranchAddress("likemuz1",&u_likemuz1);
  double u_likemub0; uLumi->SetBranchAddress("likemub0",&u_likemub0);
  double u_likemub1; uLumi->SetBranchAddress("likemub1",&u_likemub1);
  double u_unlikemuz0; uLumi->SetBranchAddress("unlikemuz0",&u_unlikemuz0);
  double u_unlikemuz1; uLumi->SetBranchAddress("unlikemuz1",&u_unlikemuz1);
  double u_unlikemub0; uLumi->SetBranchAddress("unlikemub0",&u_unlikemub0);
  double u_unlikemub1; uLumi->SetBranchAddress("unlikemub1",&u_unlikemub1);
  double u_rellumizdc; uLumi->SetBranchAddress("rellumizdc",&u_rellumizdc);
  double u_rellumizdc_err; uLumi->SetBranchAddress("rellumizdc_err",&u_rellumizdc_err);
  double u_rellumibbc; uLumi->SetBranchAddress("rellumibbc",&u_rellumibbc);
  double u_rellumibbc_err; uLumi->SetBranchAddress("rellumibbc_err",&u_rellumibbc_err);
  double u_allzdcbbc; uLumi->SetBranchAddress("allzdcbbc",&u_allzdcbbc);
  double u_allzdcbbc_err; uLumi->SetBranchAddress("allzdcbbc_err",&u_allzdcbbc_err);

  uBunch=(TTree*)uBunchFile->Get("uBunch");
  int ub_fill; uBunch->SetBranchAddress("fill",&ub_fill);
  int ub_run; uBunch->SetBranchAddress("run",&ub_run);
  int ub_pat; uBunch->SetBranchAddress("pat",&ub_pat);
  int ub_bunch; uBunch->SetBranchAddress("bunch",&ub_bunch);
  int ub_bspin; uBunch->SetBranchAddress("bspin",&ub_bspin);
  int ub_yspin; uBunch->SetBranchAddress("yspin",&ub_yspin);
  DrawPionALLwithEvenVsOdd();

  //CompareEarlyAndLateSpinConfigs();
  return;
  
  //DrawRelativeLumiZDC();
  //DrawRelativeLumiALL();
  //DrawPionALLwithArms();
  //DrawPionALLwithEvenVsOdd();
  //DrawPionALLwithPatternSetsAndEvenVsOdd();
  TString patset[]={"pat==1||pat==4||pat==21||pat==24",
		     "pat==2||pat==3||pat==22||pat==23",
 		     "pat==5||pat==8||pat==25||pat==28",
 		     "pat==6||pat==7||pat==26||pat==27"};
  TString patname[]={"(1,4,21,21)",
		      "(2,3,22,23)",
		      "(5,8,25,28)",
		      "(6,7,26,27)"};
  if(0){
    DrawPionALLwithPatternSets(); //draws lots of evens+odds together.
    //then manually adds them separated:
  c->cd(1);
  for (int i=0;i<4;i++){
    gt=GeneratePionAllTGraph("bunch%2==0",patset[i]);
    gt->SetMarkerColor(kBlack+i+1);
   gt->SetLineColor(kBlue);
    gt->SetTitle(patname[i]+ " evens");
    gt->Draw("L*");
    gt=GeneratePionAllTGraph("bunch%2==1",patset[i]);
    gt->SetTitle(patname[i]+" odds");
    gt->SetMarkerColor(kBlack+i+1);
    gt->SetLineColor(kGreen);
    gt->Draw("L*");
  }
  }

  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,700);nc++;
  c->Divide(1,2);
  c->cd(1);
  gt=GeneratePionAllTGraph("fill<17400","fill<17400");
  gt->SetTitle("MPC Mean Asymmetry, fill<17400;pT(GeV/c);asym (#)");
  gt->GetHistogram()->SetMaximum(0.2);
  gt->GetHistogram()->SetMinimum(-0.2);
  gt->Draw("A*");

   for (int i=0;i<4;i++){
     gt=GeneratePionAllTGraph("bunch%2==0 && fill<17400",Form("%s && fill<17400",patset[i].Data()));
    gt->SetMarkerColor(kBlack+i+1);
   gt->SetLineColor(kBlue);
    gt->SetTitle(patname[i]+ " evens");
    gt->Draw("L*");
    gt=GeneratePionAllTGraph("bunch%2==1&& fill<17400",Form("%s && fill<17400",patset[i].Data()));
    gt->SetTitle(patname[i]+" odds");
    gt->SetMarkerColor(kBlack+i+1);
    gt->SetLineColor(kGreen);
    gt->Draw("L*");
  }
   leg=(c->cd(1))->BuildLegend(0.80,0.25,0.99,0.75);
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel("Full set");

   c->cd(2);
  gt=GeneratePionAllTGraph("fill>17400","fill>17400");
  gt->SetTitle("MPC Mean Asymmetry, fill>17400;pT(GeV/c);asym (#)");
  gt->GetHistogram()->SetMaximum(0.2);
  gt->GetHistogram()->SetMinimum(-0.2);
  gt->Draw("A*");
  
  for (int i=0;i<4;i++){
    gt=GeneratePionAllTGraph("bunch%2==0&& fill>17400",Form("%s && fill>17400",patset[i].Data()));
    gt->SetMarkerColor(kBlack+i+1);
   gt->SetLineColor(kBlue);
    gt->SetTitle(patname[i]+ " evens");
    gt->Draw("L*");
    gt=GeneratePionAllTGraph("bunch%2==1&& fill>17400",Form("%s && fill>17400",patset[i].Data()));
    gt->SetTitle(patname[i]+" odds");
    gt->SetMarkerColor(kBlack+i+1);
    gt->SetLineColor(kGreen);
    gt->Draw("L*");
  }
   leg=(c->cd(2))->BuildLegend(0.80,0.25,0.99,0.75);
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel("Full set");
  return;
}



void CompareEarlyAndLateSpinConfigs(){
    TString patset[]={"pat==21||pat==24",
		     "pat==22||pat==23",
 		     "pat==25||pat==28",
 		     "pat==26||pat==27"};
  TString patname[]={"(21,24)",
		      "(22,23)",
		      "(25,28)",
		      "(26,27)"};
     TString earlypatset[]={"pat==1||pat==4",
		     "pat==2||pat==3",
 		     "pat==5||pat==8",
 		     "pat==6||pat==7"};
  TString earlypatname[]={"(1,4)",
		      "(2,3)",
		      "(5,8)",
		      "(6,7)"};
  TGraph *egt,*lgt;
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1600,500);nc++;
  c->Divide(4,1);

 gt=egt=GeneratePionAllTGraph("pat<20","pat<20");
  gt->SetTitle("MPC Asymmetry;pT(GeV/c);asym (#)");
  gt->GetHistogram()->SetMaximum(0.5);
  gt->GetHistogram()->SetMinimum(-0.2);
  gt->Draw("A*"); return;
  gt=lgt=GeneratePionAllTGraph("pat>20","pat>20");
  gt->SetMarkerColor(kRed);


   for (int i=0;i<4;i++){
      c->cd(1+i);
      egt->Draw("A*");
      c->cd(2);
      lgt->Draw("L*");
      return;
     gt=GeneratePionAllTGraph("bunch%2==0",patset[i]);
    gt->SetMarkerColor(kBlack);
   gt->SetLineColor(kBlue);
    gt->SetTitle(patname[i]+ " evens");
    gt->Draw("L*");
     gt=GeneratePionAllTGraph("bunch%2==1",patset[i]);
    gt->SetTitle(patname[i]+" odds");
    gt->SetMarkerColor(kBlack);
    gt->SetLineColor(kGreen);
    gt->Draw("L*");
      gt=GeneratePionAllTGraph("bunch%2==0",earlypatset[i]);
    gt->SetMarkerColor(kRed);
   gt->SetLineColor(kOrange);
    gt->SetTitle(earlypatname[i]+ " evens");
    gt->Draw("L*");
     gt=GeneratePionAllTGraph("bunch%2==1",earlypatset[i]);
    gt->SetTitle(earlypatname[i]+" odds");
    gt->SetMarkerColor(kRed);
    gt->SetLineColor(kRed);
    gt->Draw("L*");   
  
   leg=(c->cd(1+i))->BuildLegend(0.80,0.25,0.99,0.75);
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel("Full set<20");
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(1))->SetLabel("Full set>20");
   }
  return;
}

void DrawPionALL(){

  gt=GeneratePionAllTGraph("1","1");
  gt->SetTitle("MPC North+South Asymmetry;pT;A_LL");
  //gt->SetMarkerColor(kRed);
  gt->Draw("A*");

  return;
}


void DrawPionALLwithArms(){

  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run;
  vector<double>run_err;
  vector<double>all[3][nptbins],all_err[3][nptbins];

    //dividing by patterns later, maybe:
  vector<double>pat_lumi[nPats];
  vector<double>pat_lumi_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
  
  double bpol,ypol,bpol_err,ypol_err;
  TH2F* hYield[3];

  double yield[3][nptbins];
  for (int arm=0;arm<3;arm++){
    for (int j=0;j<nptbins;j++){//zero out our yield
      yield[arm][j]=0;
    }
  }
  

  uLumi->Draw("run:rellumizdc:rellumizdc_err:bpol:ypol:bpol_err:ypol_err","1","goff");
  int nRuns=uLumi->GetSelectedRows();
  for (int i=0;i<nRuns;i++){
    int thisrun=uLumi->GetVal(0)[i];
    TFile *yieldfile=NULL;
    yieldfile=TFile::Open(Form("./yields/%d.MPC.yields.rcc.hist.root",thisrun),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",thisrun);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");
    
    run.push_back(thisrun);
    run_err.push_back(0);
    double thisrel=uLumi->GetVal(1)[i];
    double thisrel_err=uLumi->GetVal(2)[i];
    lumi.push_back(thisrel);
    lumi_err.push_back(thisrel_err);

    bpol=uLumi->GetVal(3)[i];
    ypol=uLumi->GetVal(4)[i];
    bpol_err=uLumi->GetVal(5)[i];
    ypol_err=uLumi->GetVal(6)[i];
    

    //get the bunch data for this run:
    uBunch->Draw("bunch:bspin:yspin",Form("run==%d",thisrun),"goff");
    int nBunches=uBunch->GetSelectedRows();
    int nPi[2][nptbins];

    for (int arm=0;arm<3;arm++){
    for (int j=0;j<nptbins;j++){//zero out our counters for this set.
      nPi[0][j]=0;//unlike
      nPi[1][j]=0;//like
    }

    //sum the pions in each bin for each spin config:
    for (int j=0;j<nBunches;j++){
      int bun=uBunch->GetVal(0)[j];
      int bspin=uBunch->GetVal(1)[j];
      int yspin=uBunch->GetVal(2)[j];
      for (int k=0;k<nptbins;k++){
	int lowbin=hYield[arm]->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield[arm]->GetXaxis()->FindBin(pt_limits[k+1]-0.0001);
	int bunchbin=hYield[arm]->GetYaxis()->FindBin(bun);
	nPi[bspin==yspin][k]+=hYield[arm]->Integral(lowbin,highbin,bunchbin,bunchbin);
	yield[arm][k]+=hYield[arm]->Integral(lowbin,highbin,bunchbin,bunchbin);
      }
    }//j nBunches

    //calc the asymmetry for this run:
    for (int k=0;k<nptbins;k++){
      double asym, asym_err;
      CalcAsymAndErr(&asym, &asym_err,
		     bpol,  bpol_err,
		     ypol,  ypol_err,
		     nPi[0][k], sqrt( nPi[0][k]),
		     nPi[1][k], sqrt( nPi[1][k]),
		     1, 0,
		     thisrel, thisrel_err);
      if (nPi[0][k]+nPi[1][k]>0){
	all[arm][k].push_back(asym);
	all_err[arm][k].push_back(asym_err);
      } else { //if there are no events in this file, that's worrisome, but we can turn off the uncertainty for it:
	all[arm][k].push_back(0);
	all_err[arm][k].push_back(1e9);//huge uncertainty so this doesn't impact our result
	printf("No events in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
      }
      if (asym_err<0){
	printf("Negative all_err in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
	return;
      }

    }//k ptbins
    }//arm arms
    yieldfile->Close();
  }//i nRuns
  //create the weighted average of the ALLs for each ptbin:
  vector<double>final_pt[3];
  vector<double>final_pt_err;
  vector<double>final_all[3];
  vector<double>final_all_err[3];
  for (int k=0;k<nptbins;k++){
    final_pt[0].push_back(0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin
    final_pt[1].push_back(-0.1+0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin
    final_pt[2].push_back(0.1+0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin
    final_pt_err.push_back(0.5*(pt_limits[k+1]-pt_limits[k]));//and let the error bars span the bin
    //ALL_ave=Sum(ALL_i/sig_i^2)/Sum(1/sig_i^2)
    //ALL_ave_err^2=Sum( (1/sig_i^2)^2*sig_i^2)/(Sum(1/sig_i^2))^2
    //... and simplify:
    //ALL_ave_err^2=1/(Sum(1/sig_i^2))
    //1/ALL_ave_err^2=Sum(1/sig_i^2)
    for (int arm=0;arm<3;arm++){
      double weightedsum=0;
      double weight2sum=0;
      for (int i=0;i<run.size();i++){
	weightedsum+=all[arm][k][i]/(all_err[arm][k][i]*all_err[arm][k][i]);
	weight2sum+=1/(all_err[arm][k][i]*all_err[arm][k][i]);
	//  if (run[i]==389752)       printf("bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);

	printf("bin%d: r%d: sumsofar=%f, errsum=%f\n",k,run[i],weightedsum,weight2sum);
      }// i nRuns
      final_all[arm].push_back(weightedsum/weight2sum);
      final_all_err[arm].push_back(sqrt(1/weight2sum));
    }//arm arms
  }//k ptbins

  //plot our resulting asymmetry!
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,700);nc++;
  c->Divide(1,2);
  c->cd(1);
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(final_all[0][0]),&(final_pt_err[0]),&(final_all_err[0][0]));
  gt->SetTitle("MPC North+South Asymmetry;pT;A_LL");
  gt->GetHistogram()->SetMaximum(0.01);//Yaxis()->SetLimits(1,1.5e8);
  gt->GetHistogram()->SetMinimum(-0.01);//Yaxis()->SetLimits(1,1.5e8);

  gt->Draw("A*");
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[1][0]),&(final_all[1][0]),&(run_err[0]),&(final_all_err[1][0]));
  gt->SetTitle("MPC North Asymmetry;pT;A_LL");
  gt->SetLineColor(kRed);
  gt->SetMarkerColor(kRed);
  gt->Draw("*");
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[2][0]),&(final_all[2][0]),&(run_err[0]),&(final_all_err[2][0]));
  gt->SetTitle("MPC South Asymmetry;pT;A_LL");
  gt->SetLineColor(kBlue);
  gt->SetMarkerColor(kBlue);
  gt->Draw("*");
  c->cd(2)->SetLogy();
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[0][0]),&(final_pt_err[0]),&(run_err[0]));
  gt->SetTitle("MPC North+South Yield;pT;#");
  gt->GetHistogram()->SetMinimum(1);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");
  //c->cd(2)->Update();

  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[1][0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC North Yield;pT;#");
  gt->SetLineColor(kRed);
  gt->SetMarkerColor(kRed);
  gt->Draw("*");
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[2][0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC South Yield;pT;#");
  gt->SetLineColor(kBlue);
  gt->SetMarkerColor(kBlue);
  gt->Draw("*");

  double sumyield[nptbins];
  for (int j=0;j<nptbins;j++){
    sumyield[j]=yield[1][j]+yield[2][j];
  }
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(sumyield[0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC summed Yield;pT;#");
  gt->SetLineColor(kGreen);
  gt->SetMarkerColor(kGreen);
  gt->Draw("*");

  return;
}


void DrawPionALLwithEvenVsOdd(){

  int nDivisions=3;//number of ways we're splitting the data, +1 for the mainline
  
  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run;
  vector<double>run_err;
  vector<double>all[nDivisions][nptbins],all_err[nDivisions][nptbins];

    //dividing by patterns later, maybe:
  vector<double>pat_lumi[nPats];
  vector<double>pat_lumi_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
  
  double bpol,ypol,bpol_err,ypol_err;
  TH2F* hYield[nDivisions];

  double yield[nDivisions][nptbins];
  for (int div=0;div<nDivisions;div++){
    for (int j=0;j<nptbins;j++){//zero out our yield
      yield[div][j]=0;
    }
  }
  

  uLumi->Draw("run:rellumizdc:rellumizdc_err:bpol:ypol:bpol_err:ypol_err","1","goff");
  int nRuns=uLumi->GetSelectedRows();
  for (int i=0;i<nRuns;i++){
    int thisrun=uLumi->GetVal(0)[i];
    TFile *yieldfile=NULL;
    yieldfile=TFile::Open(Form("./yields/%d.MPC.yields.rcc.hist.root",thisrun),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",thisrun);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");
    
    run.push_back(thisrun);
    run_err.push_back(0);
    //get the bunch data for this run:
    uBunch->Draw("bunch:bspin:yspin",Form("run==%d",thisrun),"goff");
    int nBunches=uBunch->GetSelectedRows();
    int nPi[2][nptbins];

    for (int div=0;div<nDivisions;div++){
    //sum the pions in each bin for each spin config:
    for (int j=0;j<nBunches;j++){
      int bun=uBunch->GetVal(0)[j];
      if (div>0 && bun%2!=div%2) continue; //skip ones with the wrong bunch parity.
      for (int k=0;k<nptbins;k++){
	int lowbin=hYield[0]->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield[0]->GetXaxis()->FindBin(pt_limits[k+1]-0.0001);
	int bunchbin=hYield[0]->GetYaxis()->FindBin(bun);
	yield[div][k]+=hYield[0]->Integral(lowbin,highbin,bunchbin,bunchbin);
      }
    }//j nBunches
    }//div divs
    yieldfile->Close();
  }//i nRuns
  //create the weighted average of the ALLs for each ptbin:
  vector<double>final_pt[nDivisions];
  vector<double>final_pt_err;
  for (int k=0;k<nptbins;k++){
    final_pt[0].push_back(0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin
    final_pt_err.push_back(0.5*(pt_limits[k+1]-pt_limits[k]));//and let the error bars span the bin
  }//k ptbins

  //plot our resulting asymmetry!
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,700);nc++;
  c->Divide(1,2);
  c->cd(1);
  gt=GeneratePionAllTGraph("1","1");
  gt->SetTitle("MPC Even+Odd Asymmetry;pT;A_LL");
  gt->GetHistogram()->SetMaximum(0.01);//Yaxis()->SetLimits(1,1.5e8);
  gt->GetHistogram()->SetMinimum(-0.01);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");
  
  gt=GeneratePionAllTGraph("bunch%2==1","1");
  gt->SetTitle("MPC Odd Asymmetry;pT;A_LL");
  gt->SetLineColor(kRed);
  gt->SetMarkerColor(kRed);
  gt->Draw("*");
  
  gt=GeneratePionAllTGraph("bunch%2==0","1");
 gt->SetTitle("MPC Even Asymmetry;pT;A_LL");
  gt->SetLineColor(kBlue);
  gt->SetMarkerColor(kBlue);
  gt->Draw("*");

    gt=GenerateBunchwisePionAllTGraph("bunch%2==1","1");
  gt->SetTitle("Bunchwise Odd;pT;A_LL");
  gt->SetLineColor(kOrange);
  gt->SetMarkerColor(kOrange);
  gt->Draw("*");
  
  gt=GenerateBunchwisePionAllTGraph("bunch%2==0","1");
 gt->SetTitle("Bunchwise Even;pT;A_LL");
  gt->SetLineColor(kGreen);
  gt->SetMarkerColor(kGreen);
  gt->Draw("*");
    leg=(c->cd(1))->BuildLegend(0.80,0.25,0.99,0.75);
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel("Full set");

  c->cd(2)->SetLogy();
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[0][0]),&(final_pt_err[0]),&(run_err[0]));
  gt->SetTitle("MPC Even+Odd Yield;pT;#");
  gt->GetHistogram()->SetMinimum(1);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");
  //c->cd(2)->Update();

  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[1][0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC Odd Yield;pT;#");
  gt->SetLineColor(kRed);
  gt->SetMarkerColor(kRed);
  gt->Draw("*");
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[2][0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC Even Yield;pT;#");
  gt->SetLineColor(kBlue);
  gt->SetMarkerColor(kBlue);
  gt->Draw("*");

  double sumyield[nptbins];
  for (int j=0;j<nptbins;j++){
    sumyield[j]=yield[1][j]+yield[2][j];
  }
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(sumyield[0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC summed Yield;pT;#");
  gt->SetLineColor(kGreen);
  gt->SetMarkerColor(kGreen);
  gt->Draw("*");

  return;
}


void DrawPionALLwithPatternSets(){

  int nDivisions=5;//number of ways we're splitting the data, +1 for the mainline
  char * divName[]={"Full","(1,4,21,24)","(2,3,22,23)","(5,8,25,28)","(6,7,26,27)"};
  vector<int> divset[nDivisions];
  divset[1].push_back(1);divset[1].push_back(4);divset[1].push_back(21);divset[1].push_back(24);
  divset[2].push_back(2);divset[2].push_back(3);divset[2].push_back(22);divset[2].push_back(23);
  divset[3].push_back(5);divset[3].push_back(8);divset[3].push_back(25);divset[3].push_back(28);
  divset[4].push_back(6);divset[4].push_back(7);divset[4].push_back(26);divset[4].push_back(27);

  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run;
  vector<double>run_err;
  vector<double>all[nDivisions][nptbins],all_err[nDivisions][nptbins];

    //dividing by patterns later, maybe:
  vector<double>pat_lumi[nPats];
  vector<double>pat_lumi_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
  
  double bpol,ypol,bpol_err,ypol_err;
  TH2F* hYield[nDivisions];

  double yield[nDivisions][nptbins];
  for (int div=0;div<nDivisions;div++){
    for (int j=0;j<nptbins;j++){//zero out our yield
      yield[div][j]=0;
    }
  }
  

  uLumi->Draw("run:rellumizdc:rellumizdc_err:bpol:ypol:bpol_err:ypol_err:pat","1","goff");
  int nRuns=uLumi->GetSelectedRows();
  for (int i=0;i<nRuns;i++){
    int thisrun=uLumi->GetVal(0)[i];
    TFile *yieldfile=NULL;
    yieldfile=TFile::Open(Form("./yields/%d.MPC.yields.rcc.hist.root",thisrun),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",thisrun);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");
    
    run.push_back(thisrun);
    run_err.push_back(0);
    double thisrel=uLumi->GetVal(1)[i];
    double thisrel_err=uLumi->GetVal(2)[i];

    lumi.push_back(thisrel);
    lumi_err.push_back(thisrel_err);

    bpol=uLumi->GetVal(3)[i];
    ypol=uLumi->GetVal(4)[i];
    bpol_err=uLumi->GetVal(5)[i];
    ypol_err=uLumi->GetVal(6)[i];
    int spinpat=uLumi->GetVal(7)[i];
    

    //get the bunch data for this run:
    uBunch->Draw("bunch:bspin:yspin",Form("run==%d",thisrun),"goff");
    int nBunches=uBunch->GetSelectedRows();
    int nPi[2][nptbins];

    for (int div=0;div<nDivisions;div++){
      thisrel=lumi[lumi.size()-1];
      //if (div==2 || div==4) thisrel=1/lumi[lumi.size()-1]; //rcc ad-hoc flip.

      bool valid=(div==0);//guilty at first, unless we're doing the overall one.
      if (div>0){
	for (int q=0;q<divset[div].size();q++){
	  if (spinpat==divset[div][q]) {
	    valid=true;
	    break;
	  }
	}
      }
      if (!valid) continue; //skip ones with the wrong spinpat

      for (int j=0;j<nptbins;j++){//zero out our counters for this set.
	nPi[0][j]=0;//unlike
	nPi[1][j]=0;//like
      }

    //sum the pions in each bin for each spin config:
    for (int j=0;j<nBunches;j++){
      int bun=uBunch->GetVal(0)[j];
      int bspin=uBunch->GetVal(1)[j];
      int yspin=uBunch->GetVal(2)[j];
      int index=(bspin==yspin);
      //if (div==2 || div==4) index=(1-index); //rcc ad-hoc flip
      for (int k=0;k<nptbins;k++){
	int lowbin=hYield[0]->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield[0]->GetXaxis()->FindBin(pt_limits[k+1]-0.0001);
	int bunchbin=hYield[0]->GetYaxis()->FindBin(bun);
	nPi[index][k]+=hYield[0]->Integral(lowbin,highbin,bunchbin,bunchbin);
	yield[div][k]+=hYield[0]->Integral(lowbin,highbin,bunchbin,bunchbin);
      }
    }//j nBunches

    //calc the asymmetry for this run:
    for (int k=0;k<nptbins;k++){
      double asym, asym_err;
      CalcAsymAndErr(&asym, &asym_err,
		     bpol,  bpol_err,
		     ypol,  ypol_err,
		     nPi[0][k], sqrt( nPi[0][k]),
		     nPi[1][k], sqrt( nPi[1][k]),
		     1, 0,
		     thisrel, thisrel_err);
      if (nPi[0][k]+nPi[1][k]>0){
	all[div][k].push_back(asym);
	all_err[div][k].push_back(asym_err);
      } else { //if there are no events in this file, that's worrisome, but we can turn off the uncertainty for it:
	all[div][k].push_back(0);
	all_err[div][k].push_back(1e9);//huge uncertainty so this doesn't impact our result
	printf("No events in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
      }
      if (asym_err<0){
	printf("Negative all_err in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
	return;
      }

    }//k ptbins
    }//div divs
    yieldfile->Close();
  }//i nRuns
  //create the weighted average of the ALLs for each ptbin:
  vector<double>final_pt[nDivisions];
  vector<double>final_pt_err;
  vector<double>final_all[nDivisions];
  vector<double>final_all_err[nDivisions];
  for (int k=0;k<nptbins;k++){
    final_pt[0].push_back(0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin

    for (int div=1;div<nDivisions;div++){
      final_pt[div].push_back(-0.08*(nDivisions/2+div-1)+0.5*(pt_limits[k]+pt_limits[k+1]));//shift the datapoint in the bin
    }
    final_pt_err.push_back(0.5*(pt_limits[k+1]-pt_limits[k]));//and let the error bars span the bin
    //ALL_ave=Sum(ALL_i/sig_i^2)/Sum(1/sig_i^2)
    //ALL_ave_err^2=Sum( (1/sig_i^2)^2*sig_i^2)/(Sum(1/sig_i^2))^2
    //... and simplify:
    //ALL_ave_err^2=1/(Sum(1/sig_i^2))
    //1/ALL_ave_err^2=Sum(1/sig_i^2)
    for (int div=0;div<nDivisions;div++){
      double weightedsum=0;
      double weight2sum=0;
      for (int i=0;i<all[div][k].size();i++){
	weightedsum+=all[div][k][i]/(all_err[div][k][i]*all_err[div][k][i]);
	weight2sum+=1/(all_err[div][k][i]*all_err[div][k][i]);
	//  if (run[i]==389752)       printf("bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);

	printf("div %d bin%d: r%d: sumsofar=%f, errsum=%f\n",div,k,run[i],weightedsum,weight2sum);
      }// i nRuns
      final_all[div].push_back(weightedsum/weight2sum);
      final_all_err[div].push_back(sqrt(1/weight2sum));
    }//div divs
  }//k ptbins

  //plot our resulting asymmetry!
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,700);nc++;
  c->Divide(1,2);
  c->cd(1);
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(final_all[0][0]),&(final_pt_err[0]),&(final_all_err[0][0]));
  gt->SetTitle("MPC Asymmetry;pT;A_LL");
  gt->GetHistogram()->SetMaximum(0.03);//Yaxis()->SetLimits(1,1.5e8);
  gt->GetHistogram()->SetMinimum(-0.03);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");

  for (int div=1;div<nDivisions;div++){
    gt=new TGraphErrors(final_pt[0].size(),
			&(final_pt[div][0]),&(final_all[div][0]),
			&(run_err[0]),&(final_all_err[div][0]));
    gt->SetTitle(divName[div]);
    gt->SetLineColor(kBlack+div);
    gt->SetMarkerColor(kBlack+div);
    gt->Draw("*");
  }
  c->cd(2)->SetLogy();
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[0][0]),&(final_pt_err[0]),&(run_err[0]));
  gt->SetTitle("MPC Total Yield;pT;#");
  gt->GetHistogram()->SetMinimum(1);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");
  //c->cd(2)->Update();
 for (int div=1;div<nDivisions;div++){
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[div][0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle(divName[div]);
  gt->SetLineColor(kBlack+div);
  gt->SetMarkerColor(kBlack+div);
  gt->Draw("*");
  }

  double sumyield[nptbins];
  for (int j=0;j<nptbins;j++){
    sumyield[j]=0;
    for (int div=1;div<nDivisions;div++){
      sumyield[j]+=yield[div][j];
    }
  }
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(sumyield[0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC summed Yield;pT;#");
  gt->SetLineColor(kGreen);
  gt->SetMarkerColor(kGreen);
  gt->Draw("*");

  return;
}



void DrawPionALLwithPatternSetsAndEvenVsOdd(){

  int nDivisions=5;//number of ways we're splitting the data, +1 for the mainline
  char * divName[]={"Full","(1,4,21,24)","(2,3,22,23)","(5,8,25,28)","(6,7,26,27)"};
  vector<int> divset[nDivisions];
  divset[1].push_back(1);divset[1].push_back(4);divset[1].push_back(21);divset[1].push_back(24);
  divset[2].push_back(2);divset[2].push_back(3);divset[2].push_back(22);divset[2].push_back(23);
  divset[3].push_back(5);divset[3].push_back(8);divset[3].push_back(25);divset[3].push_back(28);
  divset[4].push_back(6);divset[4].push_back(7);divset[4].push_back(26);divset[4].push_back(27);

  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run;
  vector<double>run_err;
  vector<double>all[nDivisions][nptbins],all_err[nDivisions][nptbins];

    //dividing by patterns later, maybe:
  vector<double>pat_lumi[nPats];
  vector<double>pat_lumi_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
  
  double bpol,ypol,bpol_err,ypol_err;
  TH2F* hYield[nDivisions];

  double yield[nDivisions][nptbins];
  for (int div=0;div<nDivisions;div++){
    for (int j=0;j<nptbins;j++){//zero out our yield
      yield[div][j]=0;
    }
  }
  

  uLumi->Draw("run:rellumizdc:rellumizdc_err:bpol:ypol:bpol_err:ypol_err:pat","1","goff");
  int nRuns=uLumi->GetSelectedRows();
  for (int i=0;i<nRuns;i++){
    int thisrun=uLumi->GetVal(0)[i];
    TFile *yieldfile=NULL;
    yieldfile=TFile::Open(Form("./yields/%d.MPC.yields.rcc.hist.root",thisrun),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",thisrun);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");
    
    run.push_back(thisrun);
    run_err.push_back(0);
    double thisrel=uLumi->GetVal(1)[i];
    double thisrel_err=uLumi->GetVal(2)[i];

    lumi.push_back(thisrel);
    lumi_err.push_back(thisrel_err);

    bpol=uLumi->GetVal(3)[i];
    ypol=uLumi->GetVal(4)[i];
    bpol_err=uLumi->GetVal(5)[i];
    ypol_err=uLumi->GetVal(6)[i];
    int spinpat=uLumi->GetVal(7)[i];
    

    //get the bunch data for this run:
    uBunch->Draw("bunch:bspin:yspin",Form("run==%d",thisrun),"goff");
    int nBunches=uBunch->GetSelectedRows();
    int nPi[2][nptbins];

    for (int div=0;div<nDivisions;div++){
      thisrel=lumi[lumi.size()-1];
      if (div==2 || div==4) thisrel=1/lumi[lumi.size()-1];//rcc this is the ad-hoc lumi flip

      bool valid=(div==0);//guilty at first, unless we're doing the overall one.
      if (div>0){
	for (int q=0;q<divset[div].size();q++){
	  if (spinpat==divset[div][q]) {
	    valid=true;
	    break;
	  }
	}
      }
      if (!valid) continue; //skip ones with the wrong spinpat

      for (int j=0;j<nptbins;j++){//zero out our counters for this set.
	nPi[0][j]=0;//unlike
	nPi[1][j]=0;//like
      }

    //sum the pions in each bin for each spin config:
    for (int j=0;j<nBunches;j++){
      int bun=uBunch->GetVal(0)[j];
      int bspin=uBunch->GetVal(1)[j];
      int yspin=uBunch->GetVal(2)[j];
      int index=(bspin==yspin);
      if (div==2 || div==4) index=(1-index);
      for (int k=0;k<nptbins;k++){
	int lowbin=hYield[0]->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield[0]->GetXaxis()->FindBin(pt_limits[k+1]-0.0001);
	int bunchbin=hYield[0]->GetYaxis()->FindBin(bun);
	nPi[index][k]+=hYield[0]->Integral(lowbin,highbin,bunchbin,bunchbin);
	yield[div][k]+=hYield[0]->Integral(lowbin,highbin,bunchbin,bunchbin);
      }
    }//j nBunches

    //calc the asymmetry for this run:
    for (int k=0;k<nptbins;k++){
      double asym, asym_err;
      CalcAsymAndErr(&asym, &asym_err,
		     bpol,  bpol_err,
		     ypol,  ypol_err,
		     nPi[0][k], sqrt( nPi[0][k]),
		     nPi[1][k], sqrt( nPi[1][k]),
		     1, 0,
		     thisrel, thisrel_err);
      if (nPi[0][k]+nPi[1][k]>0){
	all[div][k].push_back(asym);
	all_err[div][k].push_back(asym_err);
      } else { //if there are no events in this file, that's worrisome, but we can turn off the uncertainty for it:
	all[div][k].push_back(0);
	all_err[div][k].push_back(1e9);//huge uncertainty so this doesn't impact our result
	printf("No events in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
      }
      if (asym_err<0){
	printf("Negative all_err in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
	return;
      }

    }//k ptbins
    }//div divs
    yieldfile->Close();
  }//i nRuns
  //create the weighted average of the ALLs for each ptbin:
  vector<double>final_pt[nDivisions];
  vector<double>final_pt_err;
  vector<double>final_all[nDivisions];
  vector<double>final_all_err[nDivisions];
  for (int k=0;k<nptbins;k++){
    final_pt[0].push_back(0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin

    for (int div=1;div<nDivisions;div++){
      final_pt[div].push_back(-0.08*(nDivisions/2+div-1)+0.5*(pt_limits[k]+pt_limits[k+1]));//shift the datapoint in the bin
    }
    final_pt_err.push_back(0.5*(pt_limits[k+1]-pt_limits[k]));//and let the error bars span the bin
    //ALL_ave=Sum(ALL_i/sig_i^2)/Sum(1/sig_i^2)
    //ALL_ave_err^2=Sum( (1/sig_i^2)^2*sig_i^2)/(Sum(1/sig_i^2))^2
    //... and simplify:
    //ALL_ave_err^2=1/(Sum(1/sig_i^2))
    //1/ALL_ave_err^2=Sum(1/sig_i^2)
    for (int div=0;div<nDivisions;div++){
      double weightedsum=0;
      double weight2sum=0;
      for (int i=0;i<all[div][k].size();i++){
	weightedsum+=all[div][k][i]/(all_err[div][k][i]*all_err[div][k][i]);
	weight2sum+=1/(all_err[div][k][i]*all_err[div][k][i]);
	//  if (run[i]==389752)       printf("bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);

	printf("div %d bin%d: r%d: sumsofar=%f, errsum=%f\n",div,k,run[i],weightedsum,weight2sum);
      }// i nRuns
      final_all[div].push_back(weightedsum/weight2sum);
      final_all_err[div].push_back(sqrt(1/weight2sum));
    }//div divs
  }//k ptbins

  //plot our resulting asymmetry!
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,700);nc++;
  c->Divide(1,2);
  c->cd(1);
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(final_all[0][0]),&(final_pt_err[0]),&(final_all_err[0][0]));
  gt->SetTitle("MPC Asymmetry;pT;A_LL");
  gt->GetHistogram()->SetMaximum(0.03);//Yaxis()->SetLimits(1,1.5e8);
  gt->GetHistogram()->SetMinimum(-0.03);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");

  for (int div=1;div<nDivisions;div++){
    gt=new TGraphErrors(final_pt[0].size(),
			&(final_pt[div][0]),&(final_all[div][0]),
			&(run_err[0]),&(final_all_err[div][0]));
    gt->SetTitle(divName[div]);
    gt->SetLineColor(kBlack+div);
    gt->SetMarkerColor(kBlack+div);
    gt->Draw("*");
  }
  c->cd(2)->SetLogy();
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[0][0]),&(final_pt_err[0]),&(run_err[0]));
  gt->SetTitle("MPC Total Yield;pT;#");
  gt->GetHistogram()->SetMinimum(1);//Yaxis()->SetLimits(1,1.5e8);
  gt->Draw("A*");
  //c->cd(2)->Update();
 for (int div=1;div<nDivisions;div++){
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(yield[div][0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle(divName[div]);
  gt->SetLineColor(kBlack+div);
  gt->SetMarkerColor(kBlack+div);
  gt->Draw("*");
  }

  double sumyield[nptbins];
  for (int j=0;j<nptbins;j++){
    sumyield[j]=0;
    for (int div=1;div<nDivisions;div++){
      sumyield[j]+=yield[div][j];
    }
  }
  gt=new TGraphErrors(final_pt[0].size(),&(final_pt[0][0]),&(sumyield[0]),&(run_err[0]),&(run_err[0]));
  gt->SetTitle("MPC summed Yield;pT;#");
  gt->SetLineColor(kGreen);
  gt->SetMarkerColor(kGreen);
  gt->Draw("*");

  return;
}



void DrawRelativeLumiALL(){
  //A_LL of the relative luminosity
  vector<double>lumiall;
  vector<double>lumiall_err;
  vector<double>run;
  vector<double>run_err;
  vector<double>pat_lumiall[nPats];
  vector<double>pat_lumiall_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
 
  uLumi->Draw("run:allzdcbbc:allzdcbbc_err","1","goff");
  int length=uLumi->GetSelectedRows();
  for (int i=0;i<length;i++){
    run.push_back(uLumi->GetVal(0)[i]);
    run_err.push_back(0);
    lumiall.push_back(uLumi->GetVal(1)[i]);
    lumiall_err.push_back(uLumi->GetVal(2)[i]);
  }

  for (int p=0;p<nPats;p++){
    uLumi->Draw("run:allzdcbbc:allzdcbbc_err",Form("pat==%d",spinpat[p]),"goff");
    int length=uLumi->GetSelectedRows();
    for (int i=0;i<length;i++){
      pat_run[p].push_back(uLumi->GetVal(0)[i]);
      pat_run_err[p].push_back(0);
      pat_lumiall[p].push_back(uLumi->GetVal(1)[i]);
      pat_lumiall_err[p].push_back(uLumi->GetVal(2)[i]);
    }
  }
       
  c=new TCanvas(Form("c%d",nc++),Form("c%d",nc),800,600);
  gt=new TGraphErrors(run.size(),&(run[0]),&(lumiall[0]),&(run_err[0]),&(lumiall_err[0]));
  gt->SetTitle("ZDC/BBC Asymmetries vs run;run;asym");
  gt->Draw("A*");
  ft=new TF1("flat","[0]",0,1e6);
  ft->SetLineColor(kBlack);
  gt->Fit(ft);
  for (int p=0;p<nPats;p++){
    gt=new TGraphErrors(pat_run[p].size(),&(pat_run[p][0]),&(pat_lumiall[p][0]),
			&(pat_run_err[p][0]),&(pat_lumiall_err[p][0]));
    ft->SetLineColor(p%5+2);
    gt->Fit(ft);

    gt->SetTitle(Form("P%d: ALL=%1.2E+%1.2E, chi2/ndf=%1.2f",spinpat[p],ft->GetParameter(0),ft->GetParError(0),ft->GetChisquare()/ft->GetNDF()));
    gt->SetMarkerColor(p%5+2);
    gt->SetLineColor(p%5+2);
    gt->Draw("*");
  }
  return;
}


void DrawRelativeLumiZDC(){
  //A_LL of the relative luminosity
  vector<double>lumi;
  vector<double>lumi_err;
  vector<double>run;
  vector<double>run_err;
  vector<double>pat_lumi[nPats];
  vector<double>pat_lumi_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
 
  uLumi->Draw("run:rellumizdc:rellumizdc_err","1","goff");
  int length=uLumi->GetSelectedRows();
  for (int i=0;i<length;i++){
    run.push_back(uLumi->GetVal(0)[i]);
    run_err.push_back(0);
    lumi.push_back(uLumi->GetVal(1)[i]);
    lumi_err.push_back(uLumi->GetVal(2)[i]);
  }

  for (int p=0;p<nPats;p++){
    uLumi->Draw("run:rellumizdc:rellumizdc_err",Form("pat==%d",spinpat[p]),"goff");
    int length=uLumi->GetSelectedRows();
    for (int i=0;i<length;i++){
      pat_run[p].push_back(uLumi->GetVal(0)[i]);
      pat_run_err[p].push_back(0);
      pat_lumi[p].push_back(uLumi->GetVal(1)[i]);
      pat_lumi_err[p].push_back(uLumi->GetVal(2)[i]);
    }
  }
       
  c=new TCanvas(Form("c%d",nc++),Form("c%d",nc),800,600);
  gt=new TGraphErrors(run.size(),&(run[0]),&(lumi[0]),&(run_err[0]),&(lumi_err[0]));
  gt->SetTitle("ZDC Relative Luminosity vs run;run;ZDC++/ZDC+-");
  gt->Draw("A*");
  ft=new TF1("flat","[0]",0,1e6);
  ft->SetLineColor(kBlack);

  for (int p=0;p<nPats;p++){
    gt=new TGraphErrors(pat_run[p].size(),&(pat_run[p][0]),&(pat_lumi[p][0]),
			&(pat_run_err[p][0]),&(pat_lumi_err[p][0]));
    ft->SetLineColor(p%5+2);
    gt->Fit(ft);

    gt->SetTitle(Form("P%d: Mean=%1.2E",spinpat[p],ft->GetParameter(0)));
    gt->SetMarkerColor(p%5+2);
    gt->SetLineColor(p%5+2);
    gt->Draw("*");
  }
  return;
}


void CalcAsymAndErr(double *asym, double *asym_err,
		    double bpol, double bpol_err,
		    double ypol, double ypol_err,
		    double unlike_signal, double unlike_signal_err,
		    double like_signal, double like_signal_err,
		    double unlike_lumi, double unlike_lumi_err,
		    double like_lumi, double like_lumi_err){

  //assumes each term in the signature is uncorrelated.
		    
  //translated from rcc_calc_all:
  double tunlike=unlike_signal;
  double tlike=like_signal;
  double trellumi=like_lumi/unlike_lumi;
  double trelunlike=trellumi*tunlike;
  double tsum=tlike+trelunlike;
  double tdiff=tlike-trelunlike;
  double tpolfactor=1/(bpol*ypol);
  double tasym=tpolfactor*(tdiff/tsum);

  double trellumi_err=sqrt(trellumi*trellumi*
			   ((like_lumi_err*like_lumi_err)/(like_lumi*like_lumi)
			    +(unlike_lumi_err*unlike_lumi_err)/(unlike_lumi*unlike_lumi)
			    ));
       
  double tbpolerrterm=-tasym/bpol * bpol_err;
  double typolerrterm=-tasym/ypol * ypol_err;

  double tlikesumerrterm=tpolfactor*(2*trelunlike)/(tsum*tsum)*like_signal_err;
  double tunlikesumerrterm=-tpolfactor*(2*trelunlike*trellumi)/(tsum*tsum)*unlike_signal_err;
  double trellumierrterm=-tpolfactor*(2*trelunlike*tunlike)/(tsum*tsum)*trellumi_err;

  double err2=(tbpolerrterm*tbpolerrterm
	       +typolerrterm*typolerrterm
	       +tlikesumerrterm*tlikesumerrterm
	       +tunlikesumerrterm*tunlikesumerrterm
	       +trellumierrterm*trellumierrterm);

  *asym=tasym;
  *asym_err=sqrt(err2);
  return;
}


TGraphErrors *GeneratePionAllTGraph(TString uBunchCut, TString uLumiCut, int arm)
{
  //given a cut we can apply to uBunch (eg a cut on bunch number, fill pattern, etc),
  //generate a tgrapherrors containing the ALL and statistical errors.
  //arm=0 for both arms, 1=north, 2=south.
  if (arm<0 || arm>2) {
    printf("arm (%d) out of bounds.\n",arm);
    return 0;
  }

  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run;
  vector<double>run_err;
  vector<double>all[nptbins],all_err[nptbins];
  
  double bpol,ypol,bpol_err,ypol_err;
  TH2F* hYield[3];

  double yield[nptbins];
  for (int j=0;j<nptbins;j++){//zero out our yield
    yield[j]=0;
  }
  

  uLumi->Draw("run:rellumizdc:rellumizdc_err:bpol:ypol:bpol_err:ypol_err:pat",uLumiCut,"goff");
  int nRuns=uLumi->GetSelectedRows();
  printf("lumiCut \"%s\" produces %d runs\n",uLumiCut.Data(),nRuns);
  for (int i=0;i<nRuns;i++){ //loop over all runs that pass our lumi-level cut:
    int thisrun=uLumi->GetVal(0)[i];
    TFile *yieldfile=NULL;
    yieldfile=TFile::Open(Form("./yields/%d.MPC.yields.rcc.hist.root",thisrun),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",thisrun);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");
    
    run.push_back(thisrun);
    run_err.push_back(0);
    double thisrel=uLumi->GetVal(1)[i];
    double thisrel_err=uLumi->GetVal(2)[i];

    lumi.push_back(thisrel);
    lumi_err.push_back(thisrel_err);

    bpol=uLumi->GetVal(3)[i];
    ypol=uLumi->GetVal(4)[i];
    bpol_err=uLumi->GetVal(5)[i];
    ypol_err=uLumi->GetVal(6)[i];
    int spinpat=uLumi->GetVal(7)[i];
    

    //get the bunch data for this run:
    uBunch->Draw("bunch:bspin:yspin",Form("(run==%d)&&(%s)",thisrun,uBunchCut.Data()),"goff");
    int nBunches=uBunch->GetSelectedRows();
    int nPi[2][nptbins];
    printf("run %d has %d valid bunches\n",thisrun,nBunches);

    // for (int div=0;div<nDivisions;div++){
    // thisrel=lumi[lumi.size()-1];
      //this is the line rcc wants to correct where it appears elsewhere:  if (div==2 || div==4) thisrel=1/lumi[lumi.size()-1];

    for (int j=0;j<nptbins;j++){//zero out our counters for this set.
      nPi[0][j]=0;//unlike
      nPi[1][j]=0;//like
    }

    //sum the pions in each valid bunch, divvying into like and unlike configurations:
    for (int j=0;j<nBunches;j++){
      int bun=uBunch->GetVal(0)[j];
      int bspin=uBunch->GetVal(1)[j];
      int yspin=uBunch->GetVal(2)[j];
      int index=(bspin==yspin);
      for (int k=0;k<nptbins;k++){
	int lowbin=hYield[arm]->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield[arm]->GetXaxis()->FindBin(pt_limits[k+1]-0.0001);
	int bunchbin=hYield[arm]->GetYaxis()->FindBin(bun);
	nPi[index][k]+=hYield[arm]->Integral(lowbin,highbin,bunchbin,bunchbin);
	yield[k]+=hYield[arm]->Integral(lowbin,highbin,bunchbin,bunchbin);
      }
    }//j nBunches

    //calc the asymmetry for this run:
    for (int k=0;k<nptbins;k++){
      double asym, asym_err;
      CalcAsymAndErr(&asym, &asym_err,
		     bpol,  bpol_err,
		     ypol,  ypol_err,
		     nPi[0][k], sqrt( nPi[0][k]),
		     nPi[1][k], sqrt( nPi[1][k]),
		     1, 0,
		     thisrel, thisrel_err);
      if (nPi[0][k]+nPi[1][k]>0){
	all[k].push_back(asym);
	all_err[k].push_back(asym_err);
      } else { //if there are no events in this file, that's worrisome, but we can turn off the uncertainty for it:
	all[k].push_back(0);
	all_err[k].push_back(1e9);//huge uncertainty so this doesn't impact our result
	printf("No events in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
      }
      if (asym_err<0){
	printf("Negative all_err in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
	return 0;
      }
    }//k ptbins
    yieldfile->Close();
  }//i nRuns
  //create the weighted average of the ALLs for each ptbin:
  vector<double>final_pt;
  vector<double>final_pt_err;
  vector<double>final_all;
  vector<double>final_all_err;
  for (int k=0;k<nptbins;k++){
    final_pt.push_back(0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin
    final_pt_err.push_back(0.5*(pt_limits[k+1]-pt_limits[k]));//and let the error bars span the bin
    //ALL_ave=Sum(ALL_i/sig_i^2)/Sum(1/sig_i^2)
    //ALL_ave_err^2=Sum( (1/sig_i^2)^2*sig_i^2)/(Sum(1/sig_i^2))^2
    //... and simplify:
    //ALL_ave_err^2=1/(Sum(1/sig_i^2))
    //1/ALL_ave_err^2=Sum(1/sig_i^2)
    double weightedsum=0;
    double weight2sum=0;
    for (int i=0;i<all[k].size();i++){
      weightedsum+=all[k][i]/(all_err[k][i]*all_err[k][i]);
      weight2sum+=1/(all_err[k][i]*all_err[k][i]);
      //  if (run[i]==389752)       printf("bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);

      printf("bin%d: r%d: sumsofar=%f, errsum=%f\n",k,run[i],weightedsum,weight2sum);
    }// i nRuns
    final_all.push_back(weightedsum/weight2sum);
    final_all_err.push_back(sqrt(1/weight2sum));
  }//k ptbins


  return new TGraphErrors(final_pt.size(),&(final_pt[0]),&(final_all[0]),&(final_pt_err[0]),&(final_all_err[0]));
}



TGraphErrors *GenerateBunchwisePionAllTGraph(TString uBunchCut, TString uLumiCut, int arm)
{
  //given a cut we can apply to uBunch (eg a cut on bunch number, fill pattern, etc),
  //generate a tgrapherrors containing the ALL and statistical errors.
  //arm=0 for both arms, 1=north, 2=south.
  if (arm<0 || arm>2) {
    printf("arm (%d) out of bounds.\n",arm);
    return 0;
  }

  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run;
  vector<double>run_err;
  vector<double>all[nptbins],all_err[nptbins];
  
  double bpol,ypol,bpol_err,ypol_err;
  TH2F* hYield[3];

  double yield[nptbins];
  for (int j=0;j<nptbins;j++){//zero out our yield
    yield[j]=0;
  }
  

  uLumi->Draw("run:rellumizdc:rellumizdc_err:bpol:ypol:bpol_err:ypol_err:pat",uLumiCut,"goff");
  int nRuns=uLumi->GetSelectedRows();
  printf("lumiCut \"%s\" produces %d runs\n",uLumiCut.Data(),nRuns);
  for (int i=0;i<nRuns;i++){ //loop over all runs that pass our lumi-level cut:
    int thisrun=uLumi->GetVal(0)[i];
    TFile *yieldfile=NULL;
    yieldfile=TFile::Open(Form("./yields/%d.MPC.yields.rcc.hist.root",thisrun),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",thisrun);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");
    
    run.push_back(thisrun);
    run_err.push_back(0);
    double thisrel;
    double thisrel_err;

    double lumilike=0;
    double lumiunlike=0;
    uLumiXL->Draw("likemuz1:unlikemuz1",Form("(run==%d)&&(%s)",thisrun,uBunchCut.Data()),"goff");
    int nXLbunches=uLumiXL->GetSelectedRows();
    for (int j=0;j<nXLbunches;j++){
      lumilike+=uLumiXL->GetVal(0)[j];
      lumiunlike+=uLumiXL->GetVal(1)[j];
    }
    thisrel=lumilike/lumiunlike;
    thisrel_err=thisrel*sqrt(1/lumilike+1/lumiunlike);
    lumi.push_back(thisrel);
    lumi_err.push_back(thisrel_err);

    bpol=uLumi->GetVal(3)[i];
    ypol=uLumi->GetVal(4)[i];
    bpol_err=uLumi->GetVal(5)[i];
    ypol_err=uLumi->GetVal(6)[i];
    int spinpat=uLumi->GetVal(7)[i];
    

    //get the bunch data for this run:
    uBunch->Draw("bunch:bspin:yspin",Form("(run==%d)&&(%s)",thisrun,uBunchCut.Data()),"goff");
    int nBunches=uBunch->GetSelectedRows();
    int nPi[2][nptbins];
    printf("run %d has %d valid bunches\n",thisrun,nBunches);

    // for (int div=0;div<nDivisions;div++){
    // thisrel=lumi[lumi.size()-1];
      //this is the line rcc wants to correct where it appears elsewhere:  if (div==2 || div==4) thisrel=1/lumi[lumi.size()-1];

    for (int j=0;j<nptbins;j++){//zero out our counters for this set.
      nPi[0][j]=0;//unlike
      nPi[1][j]=0;//like
    }

    //sum the pions in each valid bunch, divvying into like and unlike configurations:
    for (int j=0;j<nBunches;j++){
      int bun=uBunch->GetVal(0)[j];
      int bspin=uBunch->GetVal(1)[j];
      int yspin=uBunch->GetVal(2)[j];
      int index=(bspin==yspin);
      for (int k=0;k<nptbins;k++){
	int lowbin=hYield[arm]->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield[arm]->GetXaxis()->FindBin(pt_limits[k+1]-0.0001);
	int bunchbin=hYield[arm]->GetYaxis()->FindBin(bun);
	nPi[index][k]+=hYield[arm]->Integral(lowbin,highbin,bunchbin,bunchbin);
	yield[k]+=hYield[arm]->Integral(lowbin,highbin,bunchbin,bunchbin);
      }
    }//j nBunches

    //calc the asymmetry for this run:
    for (int k=0;k<nptbins;k++){
      double asym, asym_err;
      CalcAsymAndErr(&asym, &asym_err,
		     bpol,  bpol_err,
		     ypol,  ypol_err,
		     nPi[0][k], sqrt( nPi[0][k]),
		     nPi[1][k], sqrt( nPi[1][k]),
		     1, 0,
		     thisrel, thisrel_err);
      if (nPi[0][k]+nPi[1][k]>0){
	all[k].push_back(asym);
	all_err[k].push_back(asym_err);
      } else { //if there are no events in this file, that's worrisome, but we can turn off the uncertainty for it:
	all[k].push_back(0);
	all_err[k].push_back(1e9);//huge uncertainty so this doesn't impact our result
	printf("No events in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
      }
      if (asym_err<0){
	printf("Negative all_err in bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);
	return 0;
      }
    }//k ptbins
    yieldfile->Close();
  }//i nRuns
  //create the weighted average of the ALLs for each ptbin:
  vector<double>final_pt;
  vector<double>final_pt_err;
  vector<double>final_all;
  vector<double>final_all_err;
  for (int k=0;k<nptbins;k++){
    final_pt.push_back(0.5*(pt_limits[k]+pt_limits[k+1]));//center data point in the bin
    final_pt_err.push_back(0.5*(pt_limits[k+1]-pt_limits[k]));//and let the error bars span the bin
    //ALL_ave=Sum(ALL_i/sig_i^2)/Sum(1/sig_i^2)
    //ALL_ave_err^2=Sum( (1/sig_i^2)^2*sig_i^2)/(Sum(1/sig_i^2))^2
    //... and simplify:
    //ALL_ave_err^2=1/(Sum(1/sig_i^2))
    //1/ALL_ave_err^2=Sum(1/sig_i^2)
    double weightedsum=0;
    double weight2sum=0;
    for (int i=0;i<all[k].size();i++){
      weightedsum+=all[k][i]/(all_err[k][i]*all_err[k][i]);
      weight2sum+=1/(all_err[k][i]*all_err[k][i]);
      //  if (run[i]==389752)       printf("bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);

      printf("bin%d: r%d: sumsofar=%f, errsum=%f\n",k,run[i],weightedsum,weight2sum);
    }// i nRuns
    final_all.push_back(weightedsum/weight2sum);
    final_all_err.push_back(sqrt(1/weight2sum));
  }//k ptbins


  return new TGraphErrors(final_pt.size(),&(final_pt[0]),&(final_all[0]),&(final_pt_err[0]),&(final_all_err[0]));
}
