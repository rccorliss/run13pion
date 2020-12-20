
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

void PlotFullAverageAsym();



//useful things to get once and get out of the way:
TTree *uPiLumi, *uPiLumiBinning; //data and binning data.

int nBins; //number of ptbins
vector<double> ptmid,pterr;//bin centers and distances from center to edge


void rcc_asym_from_uPiLumi(){
  TFile *uPiLumiFile=TFile::Open("uPiLumi.ttree.root","READ");
  uPiLumi=(TTree*)uPiLumiFile->Get("uPiLumi");
  uPiLumiBinning=(TTree*)uPiLumiFile->Get("uPiLumiBinning");

  
  //get our binning:
  uPiLumiBinning->Draw("ptlow:pthigh","1","goff");
  nBins=uPiLumiBinning->GetSelectedRows();
  double *ptlow=uPiLumiBinning->GetVal(0);
  double *pthigh=uPiLumiBinning->GetVal(1);
  for (int i=0;i<nBins;i++){
    ptmid.push_back(0.5*(ptlow[i]+pthigh[i]));
    pterr.push_back(0.5*(ptlow[i]-pthigh[i]));
  }

  PlotFullAverageAsym();
  return;
}

void   PlotFullAverageAsym(){
  //sums yields over all bunches in all runs and produces a single asym from that.
  int nSets=4;
  vector<double> likeYield[4];
  vector<double> unlikeYield[4];
  vector<double> asym[4],asymErr[4];
  double rel[4],relErr[4];
  TString cut[4]={"(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"!(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)",
		"!(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)"};
  TString cutName[4]={"SSOO odd bunch",
		      "SSOO even bunch",
		      "OOSS odd bunch",
		      "OOSS even bunch"};

  
  //get the average polarizations:
  double bpol[4], bpolErr[4], ypol[4],ypolErr[4];
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("bpol:bpol_err:ypol:ypol_err",Form("(zdc>0)*(%s)",cut[i].Data()),"goff");
    bpol[i]=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    bpolErr[i]=SumArray(uPiLumi->GetVal(1),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    ypol[i]=SumArray(uPiLumi->GetVal(2),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    ypolErr[i]=SumArray(uPiLumi->GetVal(3),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
  }
  
  //calc the rel lumi
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("zdc",Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    int nLike=uPiLumi->GetSelectedRows();
    double zdcLike=SumArray(uPiLumi->GetVal(0),nLike);
    uPiLumi->Draw("zdc",Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    int nUnlike=uPiLumi->GetSelectedRows();
    double zdcUnlike=SumArray(uPiLumi->GetVal(0),nUnlike);
    //naively, the uncertainties in like an unlike are independent (though really they share the global scaler)
    uPiLumi->Draw("zdc",Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");//square of error term
    double zdcLikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    uPiLumi->Draw("zdc",Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double zdcUnlikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    rel[i]=(zdcLike/zdcUnlike);
    relErr[i]=rel[i]*sqrt((zdcLikeErr/zdcLike)*(zdcLikeErr/zdcLike)
			  +(zdcUnlikeErr/zdcUnlike)*(zdcUnlikeErr/zdcUnlike));
  }

  //calc the yields
  for (int i=0;i<nSets;i++){
    for (int j=0;j<nBins;j++){
      uPiLumi->Draw(Form("yield%d",j),Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
      likeYield[i].push_back(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      
      uPiLumi->Draw(Form("yield%d",j),Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
      unlikeYield[i].push_back(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    }
  }

  //calc the asyms:
  vector<double>asymcomp,asymerrcomp;
  for (int i=0;i<nSets;i++){
    for (int j=0;j<nBins;j++){
      double tasym, tasymerr;
      CalcAsymAndErr(&tasym,&tasymerr,
		     bpol[i],bpolErr[i],
		     ypol[i],ypolErr[i],
		     (unlikeYield[i])[j],sqrt((unlikeYield[i])[j]),
		     (likeYield[i])[j],sqrt((likeYield[i])[j]),
		     1,0,
		     rel[i],relErr[i]);
		     
      asym[i].push_back(tasym);
      asymErr[i].push_back(tasymerr);
      printf("%d,%d: asym=%1.2E + %1.2E\n",i,j,tasym,tasymerr);
      asymcomp.push_back(abs(tasym));
      asymerrcomp.push_back(tasymerr);
    }
  }

  TCanvas *c;
  TGraphErrors *g;
  TLegend *leg;
 c=new TCanvas("cComp","cComp",800,400);
 TGraph *g0=new TGraph(asymcomp.size(),&(asymcomp[0]),&(asymerrcomp[0]));
 g0->SetTitle("Asymmetries vs errors;abs(asym);err");
 g0->Draw("A*");
 
  c=new TCanvas("cAsym","cAsym",1600,400);
  c->Divide(2,1);
  c->cd(1);
  for (int i=0;i<nSets;i++){
    g=new TGraphErrors(nBins,&(ptmid[0]),&((asym[i])[0]),&(pterr[0]),&((asymErr[i])[0]));
    g->SetLineColor(i+1);
    g->SetMarkerColor(i+1);
    //g->SetTitle(cutName[i].Data());
    if (i==0){
      g->SetTitle("MPC asym by spin group");
      g->GetHistogram()->SetMaximum(0.15);
      g->GetHistogram()->SetMinimum(-0.15);
      g->Draw("AC*");
    } else{
          g->SetTitle(cutName[i].Data());
      g->Draw("C*");
    }
  }
  c->cd(1)->SetGridy();
  leg=c->cd(1)->BuildLegend();
  ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());

  
 //calc the bbc asyms:
  vector<double> setindex,bbcasym,bbcasymerr;
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("bbc",Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double bbcLike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
    uPiLumi->Draw("bbc",Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double bbcUnlike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
    double tasym, tasymerr;
    CalcAsymAndErr(&tasym,&tasymerr,
		   bpol[i],bpolErr[i],
		   ypol[i],ypolErr[i],
		   bbcUnlike,sqrt(bbcUnlike),
		   bbcLike,sqrt(bbcLike),
		   1,0,
		   rel[i],relErr[i]);
    setindex.push_back(i*1.0);
    bbcasym.push_back(tasym);
    bbcasymerr.push_back(tasymerr);
    printf("set %1.1f has bbc:%E,%E ==> asym=%E\n",setindex[i],bbcLike,bbcUnlike,tasym);

    }
  
  c->cd(2);
  for (int i=0;i<nSets;i++){
    g=new TGraphErrors(1,&(setindex[i]),&(bbcasym[i]),&(bbcasymerr[i]),&(bbcasymerr[i]));
    g->SetMarkerColor(i+1);
    if (i==0){
    g->SetTitle("BBC asym by spin group");
      g->GetXaxis()->SetLimits(-0.5,3.5);
      g->GetHistogram()->SetMaximum(0.002);
      g->GetHistogram()->SetMinimum(-0.002);
      g->Draw("AC*");
      //g->SetTitle(cutName[i].Data());
    } else {
      g->SetTitle(cutName[i].Data());
      g->Draw("C*");
    }
  }
  c->cd(2)->SetGridy();
  leg=c->cd(2)->BuildLegend();
  ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());


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
  printf("asym=%1.2E\terror terms: bpol=%1.2E\t ypol=%1.2E\t like=%1.2E\t unlike=%1.2E\trel=%1.2E\n",
	 tasym,tbpolerrterm,typolerrterm,tlikesumerrterm,tunlikesumerrterm,trellumierrterm);

  *asym=tasym;
  *asym_err=sqrt(err2);
  return;
}
