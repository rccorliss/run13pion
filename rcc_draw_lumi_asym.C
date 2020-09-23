

TTree *uLumi;
TTree *uBunch;
const int nPats=16;
const int spinpat[]={1,2,3,4,5,6,7,8,21,22,23,24,25,26,27,28};

  const int nptbins=10;
const double pt_limits[] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};//one more than nptbins, to cover lower and upper bounds

TCanvas *c; int nc=0; //canvas and a running counter of how many canvases I've made
TGraph *gt;
TF1 *ft;

void DrawPionALL();

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
  TFile *uBunchFile=TFile::Open("uBunch.ttree.root","READ");

  
  uLumi=(TTree*)uLumiFile->Get("uLumi");
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

  //DrawRelativeLumiZDC();
  //DrawRelativeLumiALL();
  DrawPionALL();
  
 
  return;
}


void DrawPionALL(){

  //get the run-by-run luminosity numbers:
  vector<double>lumi,lumi_err;
  vector<int>run,run_err;
  vector<double>all[nptbins],all_err[nptbins];

  //dividing by patterns later, maybe:
  vector<double>pat_lumi[nPats];
  vector<double>pat_lumi_err[nPats];
  vector<double>pat_run[nPats];
  vector<double>pat_run_err[nPats];
  
  double bpol,ypol,bpol_err,ypol_err;
 

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
    TH2F* hYield=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    //TH2F* hYieldN=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    //TH2F* hYieldS=(TH2F*)yieldfile->Get("hYieldByBunchAndPtNorth");
    
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
	int lowbin=hYield->GetXaxis()->FindBin(pt_limits[k]);
	int highbin=hYield->GetXaxis()->FindBin(pt_limits[k+1]);
	int bunchbin=hYield->GetYaxis()->FindBin(bun);
	nPi[bspin==yspin][k]+=hYield->Integral(lowbin,highbin,bunchbin,bunchbin);
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
	return;
      }

    }
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
    for (int i=0;i<run.size();i++){
      weightedsum+=all[k][i]/(all_err[k][i]*all_err[k][i]);
      weight2sum+=1/(all_err[k][i]*all_err[k][i]);
      //  if (run[i]==389752)       printf("bin%d: r%d: same=%d diff=%d, all=%f, err=%f\n",k,thisrun, nPi[1][k],nPi[0][k],asym,asym_err);

      printf("bin%d: r%d: sumsofar=%f, errsum=%f\n",k,run[i],weightedsum,weight2sum);
    }// i nRuns
    final_all.push_back(weightedsum/weight2sum);
    final_all_err.push_back(sqrt(1/weight2sum));
  }

  //plot our resulting asymmetry!
  c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);nc++;
  gt=new TGraphErrors(final_pt.size(),&(final_pt[0]),&(final_all[0]),&(final_pt_err[0]),&(final_all_err[0]));
  gt->SetTitle("MPC North+South Asymmetry;pT;A_LL");
  gt->Draw("A*");
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
