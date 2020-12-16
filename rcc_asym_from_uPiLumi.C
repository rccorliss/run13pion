
void CalcAsymAndErr(double *asym, double *asym_err,
		    double bpol, double bpol_err,
		    double ypol, double ypol_err,
		    double unlike_signal, double unlike_signal_err,
		    double like_signal, double like_signal_err,
		    double unlike_lumi, double unlike_lumi_err,
		    double like_lumi, double like_lumi_err);

float SumArray(double *arr, int n){
  float ret=0;
  for (int i=0;i<n;i++) ret+=arr[i];
  return ret;
}

void rcc_asym_from_uPiLumi(){
  TFile *uPiLumiFile=TFile::Open("uPiLumi.ttree.root","READ");
  TTree *uPiLumi=(TTree*)uPiLumiFile->Get("uPiLumi");
  TTree *uPiLumiBinning=(TTree*)uPiLumiFile->Get("uPiLumiBinning");

  
  //get our binning:
  uPiLumiBinning->Draw("ptlow:pthigh","1","goff");
  int nBins=uPiLumiBinning->GetSelectedRows();
  double *ptlow=uPiLumiBinning->GetVal(0);
  double *pthigh=uPiLumiBinning->GetVal(1);
  vector<double> ptmid;
  for (int i=0;i<nBins;i++){
    ptmid.push_back(0.5*(ptlow[i]+pthigh[i]));
  }

  int nSets=4;
  vector<double> likeyield[4];
  vector<double> unlikeyield[4];
  vector<double> asym[4];
  double rel[4];
  TString cut[4]={"(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"!(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)",
		"!(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)"};
  TString cutname[4]={"SSOO odd bunch",
		      "SSOO even bunch",
		      "OOSS odd bunch",
		      "OOSS even bunch"};

  
  //get the average polarizations:
  double bpol[4], bpolerr[4], ypol[4],ypolerr[4];
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("bpol:bpol_err:ypol:ypol_err",Form("(zdc>0)*(%s)",cut[i].Data()),"goff");
    bpol[i]=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    bpolerr[i]=SumArray(uPiLumi->GetVal(1),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    ypol[i]=SumArray(uPiLumi->GetVal(2),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    ypolerr[i]=SumArray(uPiLumi->GetVal(2),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
  }
  
  //calc the rel lumi
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("zdc",Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    int nLike=uPiLumi->GetSelectedRows();
    double zdcLike=SumArray(uPiLumi->GetVal(0),nLike);
    uPiLumi->Draw("zdc",Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    int nUnlike=uPiLumi->GetSelectedRows();
    double zdcUnlike=SumArray(uPiLumi->GetVal(0),nUnlike);
    rel[i]=(zdcLike/zdcUnlike);
  }

  //calc the yields
  for (int i=0;i<nSets;i++){
    for (int j=0;j<nBins;j++){
      uPiLumi->Draw(Form("yield%d",j),Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
      likeyield[i].push_back(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      
      uPiLumi->Draw(Form("yield%d",j),Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
      unlikeyield[i].push_back(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    }
  }

  //calc the asyms:
  for (int i=0;i<nSets;i++){
    for (int j=0;j<nBins;j++){
      double tasym, tasymerr;
      CalcAsymAndErr(&tasym,&tasymerr,
		     bpol[i],bpolerr[i],
		     ypol[i],ypolerr[i],
		     (unlikeyield[i])[j],sqrt((unlikeyield[i])[j]),
		     (likeyield[i])[j],sqrt((likeyield[i])[j]),
		     1,0,
		     rel[i],0);
		     
      asym[i].push_back(tasym);
    }
  }

  TGraph *g;
  for (int i=0;i<nSets;i++){
    g=new TGraph(nBins,&(ptmid[0]),&((asym[i])[0]));
    g->SetLineColor(i+1);
    g->SetTitle(cutname[i].Data());
    if (i==0){
      g->Draw("AC*");
    } else {
      g->Draw("C*");
    }
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
