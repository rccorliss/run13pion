//code to read in bunch by bunch scalers and provide relative luminosity per bunch
//#include <TFileCollection>
#include <TCut.h>
#include "Math/IFunction.h"

class TF1wrapper: public ROOT::Math::IBaseFunctionOneDim
{
public:
  TF1 *f,*d;
  TF1wrapper(){};
  void SetF(TF1* fin){f=(TF1*)fin->Clone(); return;};
  void SetD(TF1* din){d=(TF1*)din->Clone(); return;};
  double DoEval(double x) const{if (x<0) return f->Eval(0)+x;if(x>1) return 1+x;return f->Eval(x);};//was previously -2+x
  double DoDerive(double x) const{if (x<0) return d->Eval(0)-x; if(x>1) return 1+x; return d->Eval(x);};//was previously 10-x
  double operator() (double x) const{return f->Eval(x);};
  ROOT::Math::IBaseFunctionOneDim* Clone() const{
    TF1wrapper *n=new TF1wrapper();
    n->SetF(f);
    n->SetD(d);
    return n; };
};

TF1wrapper *globwrap;
double global_eval(double x){return globwrap->DoEval(x);}
double global_eval_gsl(double x, void *){return global_eval(x);}
double global_derive(double x){return globwrap->DoDerive(x);}
double global_derive_gsl(double x, void *){return global_derive(x);}


bool IterateRateCorrection(int length, double *rate_arr, double kn0, double ks0, double *kn_arr, double *ks_arr, double *new_kn0, double *new_ks0, double *new_mu_arr);
//iteratively correct the true rate calculation using linear fits to the set of kn(s) for each fill

int DetermineSpinPattern(vector<int> bpat, vector<int> ypat);
//match the spin orientation vectors to the spin pattern label.

int RetrieveSpinPattern(int run, TTree *t=0);

void CalcAsymAndErr(double *asym, double *asym_err,
		    double bpol, double bpol_err,
		    double ypol, double ypol_err,
		    double unlike_signal, double unlike_signal_err,
		    double like_signal, double like_signal_err,
		    double unlike_lumi, double unlike_lumi_err,
		    double like_lumi, double like_lumi_err);


void rcc_draw_lumi_plots(){
  globwrap=new TF1wrapper();
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


  //TFileCollection fc("files"); // The name is irrelevant
  //  fc.AddFromDirectory("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/*/rlstar.root");
  //can't do this this way:  t->Add("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/*/rlstar.root"); // can't glob in this fashion.

  /*
    const char *dirname="/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/star/run13pp510/";
    auto dir = gSystem->OpenDirectory(dirname);
    char *f;//[500];
    while (f = gSystem->GetDirEntry(dir)) { 
    if (!strcmp(f,".") || !strcmp(f,"..")) continue;
    t->Add(TString(dirname) + f + "/*.root");
    }
    gSystem->FreeDirectory(dir);
    printf("t has %d entries\n",t->GetEntries());
  */

  t->Add("unique_db.root");
  t->Add("/phenix/spin2/pmontu/offline/analysis/pmontu/relative_luminosity/SpinDB/unique_db.root");//same number of entries as previous, but syncs with gl1 and starscalers both.
  //one annoying catch:  all variables must be prefixed with "star_", "gl1_", or "gl1p_"
  // t->AddFileInfoList(fc.GetList());
  //delete fc;

  TCanvas *c;
  int nc=0;
  TPad *pad; //placeholder pointer to a pad


  //plot the spin pattern as a function of run.  If disabled, The rest of this will misbehave.
  int test_of_spin_pattern=RetrieveSpinPattern(398134,t);
  printf("spin pattern for 398134=%d\n",test_of_spin_pattern);
  test_of_spin_pattern=RetrieveSpinPattern(398029,t);
  printf("spin pattern for 398029=%d\n",test_of_spin_pattern);




  
  //old code to determine spin patterns etc, now folded into 'RetreiveSpinPattern' messy, I know...
  if (0){
    vector<int> spinpat;
    vector<int> runvec; 
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,300);
    //c->Divide(3,1);
    pad=(TPad*)c->cd(1);
    t->Draw("crossing:star_run:bpat:ypat","1","goff");
    int nbunches=t->GetSelectedRows();
    //find the starts to each run:
    vector<int> runstart;
    int thisrun, lastrun=0;
    for (int i=0;i<nbunches;i++){
      thisrun=t->GetVal(1)[i];
      if (thisrun!=lastrun){
	runstart.push_back(i);
	runvec.push_back(thisrun);
	lastrun=thisrun;
      }
    }
    int nruns=runstart.size();
    vector<int> bpattern[nruns], ypattern[nruns];
    int bp,yp;
    for (int i=0;i<nruns;i++){
      for (int j=0;j<24;j++){
	bp=t->GetVal(2)[runstart[i]+j];
	yp=t->GetVal(3)[runstart[i]+j];
	bpattern[i].push_back(bp);
	ypattern[i].push_back(yp);
      }

      int patternID=DetermineSpinPattern(bpattern[i],ypattern[i]);
      if (patternID<0){
	printf("Pattern not found:\n");
	if(1){//printf blue pattern
	  printf("run %d blue pat: ",runvec[i]);
	  for (int j=0;j<16;j++){
	    printf("%c",bpattern[i][j]==1?'+':bpattern[i][j]==-1?'-':'X');
	  }
	  printf("\n");
	}
	if(1){//printf yellow pattern
	  printf("run %d yell pat: ",runvec[i]);
	  for (int j=0;j<16;j++){
	    printf("%c",ypattern[i][j]==1?'+':ypattern[i][j]==-1?'-':'X');
	  }
	  printf("\n");
	}
      }
      spinpat.push_back(patternID);
    }
    if (1){
      TGraph *g=new TGraph(spinpat.size(),&(runvec[0]),&(spinpat[0]));
      g->SetTitle("Reconstructed SpinPattern ID per run;run;pattern #");
      g->Draw("A*");
    }
    
    nc++;
  }


  

  TCut rcc_cross_qa="star_bbcwide >0.05 && star_fill!=17443 && star_fill>17211"; //various cuts on basic bunch and rates
  TCut rcc_clip_loud_runs="1";
  int loud_runlist[]={391372,391581,391869,393906,398132,398026,398027,398019,398017,398018,398011,398010,398009,397176,389437,390177,391868,391870,398143,398030,398031,398029,398028,398020,398021,398014,398015,398013,398012,398005,398007,397177,397178};
  for (int i=0;i<33;i++){
    rcc_clip_loud_runs=rcc_clip_loud_runs+Form("star_run!=%d",loud_runlist[i]);
  }
    
  TCut minclocks="star_clk>1e6"; //require a bunch to have at least 1e6 live clocks.
  TCut allcuts=minclocks;
  TCut live70="star_clk/star_rclk>0.70";//require a bunch to be live at least 70% of the time.
  allcuts=allcuts && live70;
  TCut minbbcrate="star_bbcwide>0.05";//reject bunches where the bbcwide trigger isn't acting right.
  allcuts=allcuts && minbbcrate;
  TCut bbnslope="abs((star_bbncnt/star_bbcwidecnt-1.1)/(bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  allcuts=allcuts && bbnslope;
  TCut bbsslope="abs((star_bbscnt/star_bbcwidecnt-1.1)/(star_bbcwide-0.7)+0.21)<0.2"; //reject bunches where the singles to doubles rates are far from expected -- could be SBB or other detector issue.
  allcuts=allcuts && bbsslope;
  TCut abortgap="star_cross<111";
  allcuts=allcuts && abortgap;
  TCut crossing1="star_cross!=1";
  allcuts=allcuts+crossing1;
  TCut stableratio="star_zdcwidecnt/star_bbcwidecnt<0.20";
  allcuts=allcuts && stableratio;
  TCut minzdccnt="star_zdcwidecnt/star_rclk>0.001";//reject bunches where the zdc trigger isn't acting right.
  allcuts=allcuts && minzdccnt;

  TCut rccCutSetA=rcc_cross_qa+rcc_clip_loud_runs+abortgap+crossing1;


  //show the SetB parameters we're cutting on.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw("star_bbcwiderun/gl1_bbcwidelive",rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(2)->SetLogy();
    t->Draw("star_bbc30run/gl1_bbc30live",rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(3)->SetLogy();
    t->Draw("star_zdcwiderun/gl1_zdcwidelive",rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(4)->SetLogy();
    t->Draw("star_bbcwiderun/gl1_bbcwidelive","star_clkrun"*(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(5)->SetLogy();
    t->Draw("star_bbc30run/gl1_bbc30live","star_clkrun"*(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(6)->SetLogy();
    t->Draw("star_zdcwiderun/gl1_zdcwidelive","star_clkrun"*(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //check if those values are correlated with run length:
  TCut bbcwidescalerratio="(star_bbcwiderun/gl1_bbcwidelive>0.99 && star_bbcwiderun/gl1_bbcwidelive<1.015)";
  TCut bbc30scalerratio="(star_bbc30run/gl1_bbc30live>0.99 && star_bbc30run/gl1_bbc30live<1.015)";
  TCut zdcwidescalerratio="(star_zdcwiderun/gl1_zdcwidelive>0.99 && star_zdcwiderun/gl1_zdcwidelive<1.015)";
  TCut starvsgl1ratios=bbcwidescalerratio+bbc30scalerratio+zdcwidescalerratio;
  TCut rccCutSetB=starvsgl1ratios;
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw(Form("%s>>bbcwidecut(2,-0.5,1.5)",bbcwidescalerratio.GetTitle()),"star_clkrun"*rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(2)->SetLogy();
    t->Draw(Form("%s>>bbc30cut(2,-0.5,1.5)",bbc30scalerratio.GetTitle()),rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(3)->SetLogy();
    t->Draw(Form("%s>>zdcwidecut(2,-0.5,1.5)",zdcwidescalerratio.GetTitle()),rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0","colz");
    c->cd(4);//->SetLogy();
    t->Draw("star_bbcwiderun/gl1_bbcwidelive:star_clkrun",(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(5);//->SetLogy();
    t->Draw("star_bbc30run/gl1_bbc30live:star_clkrun",(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    c->cd(6);//->SetLogy();
    t->Draw("star_zdcwiderun/gl1_zdcwidelive:star_clkrun",(rcc_cross_qa+rcc_clip_loud_runs+"star_cross==0"),"");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //some checks trying to understand what gl1p_ scalers are doing.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1);//->SetLogy();
    t->Draw("star_zdc30:log10(gl1p_zdc_narrow)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(2);//->SetLogy();
    t->Draw("star_zdcwide:log10(gl1p_zdc_wide)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(3);//->SetLogy();
    t->Draw("star_bbc30:log10(gl1p_bbc_30)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(4);//->SetLogy();
    t->Draw("log10(star_zdc30/gl1p_zdc_narrow):log10(star_clkrun)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(5);//->SetLogy();
    t->Draw("log10(star_zdcwide/gl1p_zdc_wide):log10(star_clkrun)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(6);//->SetLogy();
    t->Draw("log10(star_bbc30/gl1p_bbc_30):log10(star_clkrun)",(rccCutSetB+rccCutSetA),"colz");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //plot that shows the way the star and gl1p scalers track each other
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw("log10(star_zdc30cnt/gl1p_zdc_narrow)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(2)->SetLogy();
    t->Draw("log10(star_zdcwidecnt/gl1p_zdc_wide)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(3)->SetLogy();
    t->Draw("log10(star_bbc30cnt/gl1p_bbc_30)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(4);//->SetLogy();
    t->Draw("log10(star_zdc30cnt):log10(gl1p_zdc_narrow)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(5);//->SetLogy();
    t->Draw("log10(star_zdcwidecnt):log10(gl1p_zdc_wide)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(6);//->SetLogy();
    t->Draw("log10(star_bbc30cnt):log10(gl1p_bbc_30)",(rccCutSetB+rccCutSetA),"colz");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  TCut starvsgl1pzdcwideratio="abs(star_zdcwidecnt/gl1p_zdc_wide-1)<0.001";
  TCut starvsgl1pzdc30ratio="abs(star_zdc30cnt/gl1p_zdc_narrow-1)<0.001";
  TCut starvsgl1pbbcratio="abs(star_bbc30cnt/gl1p_bbc_30-1)<0.001";
  TCut starvsgl1pratios=starvsgl1pzdcwideratio+starvsgl1pzdc30ratio+starvsgl1pbbcratio;
  rccCutSetB+=starvsgl1pratios;
  //plots showing the cut applied
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(3,2);
    c->cd(1)->SetLogy();
    t->Draw("(star_zdc30cnt/gl1p_zdc_narrow)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(2)->SetLogy();
    t->Draw("(star_zdcwidecnt/gl1p_zdc_wide)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(3)->SetLogy();
    t->Draw("(star_bbc30cnt/gl1p_bbc_30)","star_clkrun"*(rccCutSetB+rccCutSetA),"colz");
    c->cd(4);//->SetLogy();
    t->Draw("log10(star_zdc30cnt):log10(gl1p_zdc_narrow)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(5);//->SetLogy();
    t->Draw("log10(star_zdcwidecnt):log10(gl1p_zdc_wide)",(rccCutSetB+rccCutSetA),"colz");
    c->cd(6);//->SetLogy();
    t->Draw("log10(star_bbc30cnt):log10(gl1p_bbc_30)",(rccCutSetB+rccCutSetA),"colz");
    // c->cd(4)->SetLogy();
    //t->Draw("star_zdc30cnt/gl1_zdc30live",rcc_cross_qa+rcc_clip_loud_runs,"colz");
    nc++;
  }

  //to ease redundant cuts, let's define a TEntryList and set it:
  t->Draw(">>elistAB", rccCutSetA+rccCutSetB, "entrylist");
  t->SetEntryList((TEntryList*)gDirectory->Get("elistAB"));

 

  //deprecated way to show the trendlines
  if (0){
    TF1 *test=new TF1("testpoly","[0]+[1]*x+[2]*x^2",0.1,0.6);
    TGraph *g;
    TLatex tex;
    float texpos=0.75;
    float texshift=0;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.03);
    int fiti=0;
    double pars[4][3];
    double parerrs[4][3];
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    c->Divide(2,2);
    c->cd(1);//->SetLogy();
    t->Draw("star_zdscnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,0.25,Form("y=(%1.3E #pm %1.3E) + (%1.3E #pm %1.3E)x + (%1.3E #pm %1.3E)x^2",
				 pars[fiti][0],parerrs[fiti][0],
				 pars[fiti][1],parerrs[fiti][1],
				 pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;

    c->cd(2);//->SetLogy();
    t->Draw("star_zdncnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,texpos,Form("y=(%1.3E#pm %1.3E) + (%1.3E#pm %1.3E)x + (%1.3E#pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;    c->cd(3);//->SetLogy();
    t->Draw("star_bbscnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,texpos,Form("y=(%1.3E#pm %1.3E) + (%1.3E#pm %1.3E)x + (%1.3E#pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;    c->cd(4);//->SetLogy();
    t->Draw("star_bbncnt/star_clk:star_bbcwide");
    g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    g->Fit(test);
    for (int i=0;i<3;i++){
      pars[fiti][i]=test->GetParameter(i);
      parerrs[fiti][i]=test->GetParError(i);
    }
    test->Draw("same");
    tex.DrawLatex(0.05,texpos,Form("y=(%1.3E#pm %1.3E) + (%1.3E#pm %1.3E)x + (%1.3E#pm %1.3E)x^2",
				   pars[fiti][0],parerrs[fiti][0],
				   pars[fiti][1],parerrs[fiti][1],
				   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
    fiti++;
    nc++;
  }



 
  if (0){
    TF1 *test=new TF1("testpoly","[0]+[1]*x+[2]*x^2",0.1,0.6);
    TGraph *g;
    TLatex tex;
    float texpos=0.75;
    float texshift=0;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.03);
    int fiti=0;
    double pars[4][3];
    double parerrs[4][3];
    //char *cuttemplate="%%f+%%f*%%s+%%f*pow(%%s,2)-%%s";//abs([0]+[1]*x+[2]*x^2-y)<thresh
    char xtemp[100];
    char ytemp[100];
    char *temptrimrange;
    fiti=0;
    char ych[4][100]={"star_zdscnt/star_clk","star_zdncnt/star_clk","star_bbscnt/star_clk","star_bbncnt/star_clk"};
    char xch[100]="star_bbcwide";
    TCut tempcut[4];
    float thresh[3]={0.1,0.015,0.015};
    TLine *line;
    TH1 *h;
    for (int pass=0;pass<3;pass++){
      c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,600);
      c->Divide(4,2);
      for (int fiti=0;fiti<4;fiti++){
	if (pass==0) tempcut[fiti]="1";
	pad=(TPad*)c->cd(fiti+1);//->SetLogy();
	t->Draw(Form("%s:%s",ych[fiti],xch),tempcut[fiti],"colz");
	g=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
	g->Fit(test);
	for (int i=0;i<3;i++){
	  pars[fiti][i]=test->GetParameter(i);
	  parerrs[fiti][i]=test->GetParError(i);
	}
	test->DrawCopy("same");
	h=(TH1*)pad->GetPrimitive("htemp");
	float ymid=(0.9*h->GetYaxis()->GetXmax()+(1.0-0.9)*h->GetYaxis()->GetXmin());
	float xmid=(0.05*h->GetXaxis()->GetXmax()+(1.0-0.05)*h->GetXaxis()->GetXmin());
	tex.DrawLatex(xmid,ymid,
		      Form("y=(%1.3E #pm %1.3E) + (%1.3E #pm %1.3E)x + (%1.3E #pm %1.3E)x^2",
			   pars[fiti][0],parerrs[fiti][0],
			   pars[fiti][1],parerrs[fiti][1],
			   pars[fiti][2],parerrs[fiti][2]));texpos-=texshift;
	temptrimrange=Form("abs(%f+%f*%s+%f*pow(%s,2)-%s)",
			   pars[fiti][0],pars[fiti][1],xch,pars[fiti][2],xch,ych[fiti]);
	printf("trying %s\n",temptrimrange);
	pad=(TPad*)c->cd(fiti+1+4);
	pad->SetLogy();
	t->Draw(temptrimrange,tempcut[fiti]);
	h=(TH1*)pad->GetPrimitive("htemp");
	line=new TLine(thresh[pass],1,thresh[pass],1000);
	line->SetLineColor(kRed);
	line->Draw();
	tempcut[fiti]=Form("%s<%f",temptrimrange,thresh[pass]);
      }
      nc++;
    }
    for (int i=0;i<4;i++){
      tempcut[i].Print();
    }
  }


  TCut zdspoly="abs(0.004996+0.318485*star_bbcwide+0.179259*pow(star_bbcwide,2)-star_zdscnt/star_clk)<0.015000";
  TCut zdnpoly="abs(0.006605+0.297460*star_bbcwide+0.184102*pow(star_bbcwide,2)-star_zdncnt/star_clk)<0.015000";
  TCut bbspoly="abs(-0.003347+1.252117*star_bbcwide+-0.248236*pow(star_bbcwide,2)-star_bbscnt/star_clk)<0.015000";
  TCut bbnpoly="abs(-0.002720+1.243970*star_bbcwide+-0.231012*pow(star_bbcwide,2)-star_bbncnt/star_clk)<0.015000";
  TCut rccCutSetC=zdspoly+zdnpoly+bbspoly+bbnpoly;

  //to ease redundant cuts, let's define a TEntryList and set it:
  t->Draw(">>elistABC", rccCutSetA+rccCutSetB+rccCutSetC, "entrylist");
  t->SetEntryList((TEntryList*)gDirectory->Get("elistABC"));


  TLine newline;
  newline.SetLineColor(kRed);
  //look at livetime vs run for the remaining xing+runs:
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("star_clk/star_rclk:star_run","1","colz");
    newline.DrawLine(386,0.7,399,0.7);
    c->cd(2);
    t->Draw("star_clk/star_rclk:star_run","star_bbcwidecnt","colz");
    newline.DrawLine(386,0.7,399,0.7);
    c->cd(3);
    t->Draw("star_clk/star_rclk:crossing","1","colz");
    newline.DrawLine(0,0.7,111,0.7);
    c->cd(4);
    t->Draw("star_clk/star_rclk","star_bbcwidecnt","colz");
    newline.DrawLine(0.7,0,0.7,1e12);

    nc++;
  }
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),600,400);
    c->Divide(2,1);
    c->cd(1);
    t->Draw("star_clk/star_rclk:gl1_clocklive/gl1_clockraw","1","colz");
    newline.DrawLine(0,0.7,1.0,0.7);
    newline.SetLineColor(kGreen);
    newline.DrawLine(0,0.0,1.0,1.0);
    c->cd(2);
    t->Draw("(star_clk/star_rclk)/(gl1_clocklive/gl1_clockraw)","1");
    nc++;
  }

  TCut rccCutSetD=live70;

  //to ease redundant cuts, let's define a TEntryList and set it:
  t->Draw(">>elistABCD", rccCutSetA+rccCutSetB+rccCutSetC+rccCutSetD, "entrylist");
  t->SetEntryList((TEntryList*)gDirectory->Get("elistABCD"));


  //look at how many crossings survive per run
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,300);
    c->Divide(3,1);
    pad=(TPad*)c->cd(1);
    t->Draw("crossing:star_run");
    int nbunches=t->GetSelectedRows();
    int runlist[2000];
    int buncheshere=0;
    int rt;
    int rtold=t->GetV2()[0];
    TH2F *hXingPerRun=new TH2F("hXingPerRun","Remaining Xings Per Run;run;crossings",100,386e3,399e3,120,-0.5,119.5);
    TH1F *hXingsLeft=new TH1F("hXingsLeft","Remaining Xings Per Run;crossings",120,-0.5,119.5);
    
    for (int i=0;i<nbunches;i++){
      rt=t->GetV2()[i];
      if (rt!=rtold){
	hXingPerRun->Fill(rtold,buncheshere);
	hXingsLeft->Fill(buncheshere);
	if (buncheshere<100) printf("star_run!=%d && ",rtold);
	rtold=rt;
	buncheshere=0;
      }
      buncheshere++;
    }
    if (buncheshere<100) printf("star_run!=%d && ",rtold);
    printf("\n");

    hXingPerRun->Fill(rtold,buncheshere);
    c->cd(2);
    hXingPerRun->Draw("colz");
    c->cd(3);
    hXingsLeft->Draw();
    newline.SetLineColor(kRed);
    newline.DrawLine(100,0,100,400);
  }

  TCut xing100="star_run!=387788 && star_run!=392820 && star_run!=392842 && star_run!=393531 && star_run!=394060 && star_run!=394067 && star_run!=394526 && star_run!=395734 && star_run!=397205 && star_run!=397531";


  //look at the BBc to ZDC ratios in the surviving bunches vs the run average.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),900,700);
    c->Divide(3,2);
    pad=(TPad*)c->cd(1);
    t->Draw("star_zdcwide/star_bbcwide:star_run");
    t->Draw("star_zdcwide/star_bbcwide:star_run:crossing","1","goff");
    int nel=t->GetSelectedRows();
    int runlist[2000];
    float ratiosumhere=0;
    int buncheshere=0;
    int rt;
    int rtold=t->GetV2()[0];
    TH2F *hRatioAvePerRun=new TH2F("hRatioAvePerRun","Average ZDC/BBCwide ratio in remaining Xings Per Run;run;ave",100,386e3,399e3,100,0.05,0.2);
    TH2F *hRatioXingVsRun=new TH2F("hRatioXingVsRun","Bunch ZDCwide/BBCwide Ratio vs Run Ave ratio;run ave;bunch ave",100,0.05,.2,100,0.05,.2);
    TH2F *hRatioDiffPerXing=new TH2F("hRatioDiffPerXing","Bunch Ratio/ Run Ave BBcwide;run,ratio",100,386e3,399e3,100,0.5,1.5);
    TH2F *hImpactedPairs=new TH2F("hImpactedPairs","Bunches With >20% Deviation from run ave zdc/bbc ratio;run;bunch",80,391e3,397e3,120,-0.5,119.5);
    TH1F *hImpactedBunches=new TH1F("hImpactedBunches","Bunches With >20% Deviation from run ave zdc/bbc ratio;bunch",120,-0.5,119.5);
    for (int i=0;i<nel;i++){
      rt=t->GetV2()[i];
      if (rt!=rtold){
	float averatio=ratiosumhere/buncheshere;
	hRatioAvePerRun->Fill(rtold,averatio);
	for (int j=i-buncheshere;j<i;j++){
	  //go back and fill in the bunches with this value.
	  hRatioXingVsRun->Fill(averatio,t->GetV1()[j]);
	  hRatioDiffPerXing->Fill(rtold,t->GetV1()[j]/averatio);
	  if (abs(t->GetV1()[j]/averatio-1)>0.2){
	    printf("!(star_run==%d && crossing==%1.0f) && ",rtold,t->GetV3()[j]);
	    hImpactedPairs->Fill(rtold,t->GetV3()[j]);
	    hImpactedBunches->Fill(t->GetV3()[j]);
	  }
	}
	rtold=rt;
	buncheshere=0;
	ratiosumhere=0;
      }
      ratiosumhere+=t->GetV1()[i];
      buncheshere++;
    }
    printf("\n");
    c->cd(2);
    hRatioAvePerRun->Draw("colz");
    c->cd(3);
    hRatioXingVsRun->Draw("colz");
    c->cd(4);
    hRatioDiffPerXing->Draw("colz");
    newline.SetLineColor(kRed);
    newline.DrawLine(386e3,0.8,399e3,0.8);
    newline.DrawLine(386e3,1.2,399e3,1.2);
    c->cd(5);
    hImpactedPairs->Draw("colz");
    c->cd(6);
    hImpactedBunches->Draw();
    
  }

  TCut badbunchratio[4];
  badbunchratio[0]="!(star_run==391566 && crossing==18) && !(star_run==391567 && crossing==18) && !(star_run==391569 && crossing==18) && !(star_run==391966 && crossing==36) && !(star_run==391966 && crossing==39) && !(star_run==391967 && crossing==31) && !(star_run==391967 && crossing==33) && !(star_run==391967 && crossing==34) && !(star_run==391967 && crossing==35) && !(star_run==391967 && crossing==36) && !(star_run==391967 && crossing==37) && !(star_run==391967 && crossing==39) && !(star_run==391968 && crossing==31) && !(star_run==391968 && crossing==33) && !(star_run==391968 && crossing==34) && !(star_run==391968 && crossing==35) && !(star_run==391968 && crossing==36) && !(star_run==391968 && crossing==37) && !(star_run==391968 && crossing==38) && !(star_run==391968 && crossing==39) && !(star_run==391969 && crossing==31) && !(star_run==391969 && crossing==33) && !(star_run==391969 && crossing==34) && !(star_run==391969 && crossing==35) && !(star_run==391969 && crossing==36) && !(star_run==391969 && crossing==38) && !(star_run==391969 && crossing==39) && !(star_run==391970 && crossing==33) && !(star_run==391970 && crossing==34) && !(star_run==391970 && crossing==35) && !(star_run==391970 && crossing==36) && !(star_run==391970 && crossing==39) && !(star_run==393530 && crossing==75) && !(star_run==393530 && crossing==78) && !(star_run==393534 && crossing==72) && !(star_run==393534 && crossing==75)";
 
  badbunchratio[1]="!(star_run==393534 && crossing==77) && !(star_run==393534 && crossing==78) && !(star_run==394002 && crossing==75) && !(star_run==394002 && crossing==78) && !(star_run==394003 && crossing==75) && !(star_run==394003 && crossing==78) && !(star_run==394004 && crossing==75) && !(star_run==394004 && crossing==78) && !(star_run==394005 && crossing==75) && !(star_run==394005 && crossing==78) && !(star_run==394368 && crossing==71) && !(star_run==394368 && crossing==72) && !(star_run==394368 && crossing==78) && !(star_run==394388 && crossing==71) && !(star_run==394388 && crossing==73) && !(star_run==394388 && crossing==74) && !(star_run==394388 && crossing==75) && !(star_run==394388 && crossing==76) && !(star_run==394388 && crossing==77) && !(star_run==394388 && crossing==78) && !(star_run==394388 && crossing==79) && !(star_run==394389 && crossing==71) && !(star_run==394389 && crossing==72) && !(star_run==394389 && crossing==73) && !(star_run==394389 && crossing==74) && !(star_run==394389 && crossing==75) && !(star_run==394389 && crossing==76) && !(star_run==394389 && crossing==77) && !(star_run==394389 && crossing==78) && !(star_run==394389 && crossing==79) && !(star_run==394390 && crossing==71) && !(star_run==394390 && crossing==72) && !(star_run==394390 && crossing==73) && !(star_run==394390 && crossing==74) && !(star_run==394390 && crossing==75) && !(star_run==394390 && crossing==76)";
 
  badbunchratio[2]="!(star_run==394390 && crossing==77) && !(star_run==394390 && crossing==78) && !(star_run==394390 && crossing==79) && !(star_run==394391 && crossing==71) && !(star_run==394391 && crossing==72) && !(star_run==394391 && crossing==73) && !(star_run==394391 && crossing==74) && !(star_run==394391 && crossing==75) && !(star_run==394391 && crossing==76) && !(star_run==394391 && crossing==77) && !(star_run==394391 && crossing==78) && !(star_run==394391 && crossing==79) && !(star_run==395228 && crossing==0) && !(star_run==395229 && crossing==0) && !(star_run==395230 && crossing==0) && !(star_run==395231 && crossing==0) && !(star_run==395233 && crossing==0) && !(star_run==395239 && crossing==0) && !(star_run==395242 && crossing==0) && !(star_run==395244 && crossing==0) && !(star_run==395390 && crossing==0) && !(star_run==395397 && crossing==0) && !(star_run==395402 && crossing==0) && !(star_run==395405 && crossing==0) && !(star_run==395407 && crossing==0) && !(star_run==395408 && crossing==0) && !(star_run==395411 && crossing==0) && !(star_run==395413 && crossing==0) && !(star_run==395419 && crossing==0) && !(star_run==395420 && crossing==0) && !(star_run==395421 && crossing==0) && !(star_run==395429 && crossing==0) && !(star_run==395430 && crossing==0) && !(star_run==395432 && crossing==0) && !(star_run==395526 && crossing==0) && !(star_run==395527 && crossing==0)";

  badbunchratio[3]="!(star_run==395544 && crossing==0) && !(star_run==395544 && crossing==54) && !(star_run==395544 && crossing==68) && !(star_run==395545 && crossing==0) && !(star_run==395545 && crossing==54) && !(star_run==395545 && crossing==68) && !(star_run==395549 && crossing==0) && !(star_run==395549 && crossing==54) && !(star_run==395549 && crossing==68) && !(star_run==395550 && crossing==0) && !(star_run==395550 && crossing==54) && !(star_run==395550 && crossing==68) && !(star_run==395551 && crossing==0) && !(star_run==395551 && crossing==68) && !(star_run==395552 && crossing==0) && !(star_run==395552 && crossing==68) && !(star_run==395553 && crossing==0) && !(star_run==395553 && crossing==68) && !(star_run==395587 && crossing==0) && !(star_run==395588 && crossing==0) && !(star_run==395589 && crossing==0) && !(star_run==395590 && crossing==0) && !(star_run==395591 && crossing==0) && !(star_run==395883 && crossing==30) && !(star_run==395884 && crossing==30)";
 
  TCut rccCutSetE=xing100;

  //to ease redundant cuts, let's define a TEntryList and set it:
  for (int i=0;i<4;i++){
    t->Draw(Form(">>elistABCDEpart%d",i), rccCutSetE+badbunchratio[i], "entrylist");
    t->SetEntryList((TEntryList*)gDirectory->Get(Form("elistABCDEpart%d",i)));
  }


  //before we get into the kn and ks stuff, let's make sure we know what's in the ttree, by checking whether certain variables are indeed related:
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("star_bbscnt/star_clk:star_bbs","1");
    c->cd(2);
    t->Draw("star_bbncnt/star_clk:star_bbn","1");
    c->cd(3);
    t->Draw("star_zdncnt/star_clk:star_zdn","1");
    c->cd(4);
    t->Draw("star_zdncnt/star_clk:star_zdn","1");
    nc++;
  }

  const char* ch_kn_bbc="(log(1-star_bbs)-log(1-star_bbn-star_bbs+star_bbcwide))/(log(1-star_bbn-star_bbs+star_bbcwide)-log(1-star_bbs)-log(1-star_bbn))";
  const char* ch_ks_bbc="(log(1-star_bbn)-log(1-star_bbn-star_bbs+star_bbcwide))/(log(1-star_bbn-star_bbs+star_bbcwide)-log(1-star_bbs)-log(1-star_bbn))";
  const char* ch_kn_zdc="(log(1-star_zds)-log(1-star_zdn-star_zds+star_zdcwide))/(log(1-star_zdn-star_zds+star_zdcwide)-log(1-star_zds)-log(1-star_zdn))";
  const char* ch_ks_zdc="(log(1-star_zdn)-log(1-star_zdn-star_zds+star_zdcwide))/(log(1-star_zdn-star_zds+star_zdcwide)-log(1-star_zds)-log(1-star_zdn))";
  const char* ch_mu_bbc="(log(1-star_bbn-star_bbs+star_bbcwide)-log(1-star_bbs)-log(1-star_bbn))";
  const char* ch_mu_zdc="(log(1-star_zdn-star_zds+star_zdcwide)-log(1-star_zds)-log(1-star_zdn))";
  //draw our four exclusive singles to doubles ratios first pass:
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(4,2);
    pad=(TPad*)c->cd(1);
    t->Draw(Form("%s:star_bbcwide",ch_kn_bbc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("kn_bbc vs bbcwide rate;bbcwide;kn_bbc");
    pad=(TPad*)c->cd(2);
    t->Draw(Form("%s:star_bbcwide",ch_ks_bbc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("ks_bbc vs bbcwide rate;bbcwide;ks_bbc");
    pad=(TPad*)c->cd(3);
    t->Draw(Form("%s:star_zdcwide",ch_kn_zdc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("kn_zdc vs zdcwide rate;zdcwide;kn_zdc");
    pad=(TPad*)c->cd(4);
    t->Draw(Form("%s:star_zdcwide",ch_ks_zdc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("ks_zdc vs zdcwide rate;zdcwide;ks_zdc");
    pad=(TPad*)c->cd(5);
    t->Draw(Form("%s:%s",ch_kn_bbc,ch_mu_bbc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("kn_bbc vs rate-corrected bbc rate;bbc corrected;kn_bbc");
    pad=(TPad*)c->cd(6);
    t->Draw(Form("%s:%s",ch_ks_bbc,ch_mu_bbc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("ks_bbc vs rate-corrected bbc rate;bbc corrected;ks_bbc");
    pad=(TPad*)c->cd(7);
    t->Draw(Form("%s:%s",ch_kn_zdc,ch_mu_zdc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("kn_zdc vs rate-corrected zdc rate;zdc corrected;kn_zdc");
    pad=(TPad*)c->cd(8);
    t->Draw(Form("%s:%s",ch_ks_zdc,ch_mu_zdc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("ks_zdc vs rate-corrected zdc rate;zdc corrected;ks_zdc");
    nc++;
  }

  //check if the 'corr' variables are the first-pass corrections (they're not.)
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    pad=(TPad*)c->cd(1);
    t->Draw(Form("star_bbcwidecorr:%s",ch_mu_bbc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("bbcwidecorr vs rate-corrected bbc rate;bbc math;bbcwidecorr");
    pad=(TPad*)c->cd(2);
    t->Draw(Form("star_zdcwidecorr:%s",ch_mu_zdc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("zdcwidecorr vs rate-corrected zdc rate;zdc math;zdcwidecorr");
    pad=(TPad*)c->cd(3);
    t->Draw(Form("star_bbcwidecorr/(%s)",ch_mu_bbc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("bbcwidecorr / rate-corrected bbc rate;ratio");
    pad=(TPad*)c->cd(4);
    t->Draw(Form("star_zdcwidecorr/(%s)",ch_mu_zdc),"1");
    ((TH1*)pad->GetPrimitive("htemp"))->SetTitle("zdcwidecorr / rate-corrected zdc rate;ratio");
    nc++;
  }


  if (0){
    TF1 *rateobs=new TF1("rateobs","-[0]+1-exp(-([1]+1)*x)-exp(-([2]+1)*x)+exp(-([1]+[2]+1)*x)",0,1);
    TF1 *ratederivative=new TF1("ratederivative",
				"([1]+1)*exp(-([1]+1)*x)+([2]+1)*exp(-([2]+1)*x)-([1]+[2]+1)*exp(-([1]+[2]+1)*x)",0,1);
    //have to use global.  sigh.TF1wrapper *ratewrap=new TF1wrapper();
    rateobs->SetTitle("Random Subset of 0=Robs-f(kn,ks,#mu);#mu;'0'");
    globwrap->SetF(rateobs);
    globwrap->SetD(ratederivative);
    ROOT::Math::RootFinder *rf4 = new ROOT::Math::RootFinder(ROOT::Math::RootFinder::kGSL_STEFFENSON);
    ROOT::Math::GradFunctor1D gfunc( &global_eval, &global_derive);

    TH2F *hkns[2];
    TH1F *hkn[2];
    TH1F *hks[2];
    TH1F *hrns[2];
    TH1F *hmuns[2];
    TH2F *hmurns[2];
    TH2F *hmumu[2];
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    int nplots=5;
    c->Divide(nplots,2);
    c->cd(1);
    char *trigname[]={"bbc","zdc"};
    for (int j=0;j<2;j++){
      float l=!j?0.15:3.2;
      float h=!j?0.3:4.5;
      float fl=!j?0.0:0.0;
      float fh=!j?1.0:0.12;
      float rl=!j?0.8:0.14;
      float rh=!j?1.8:0.30;
      char *trig=trigname[j];
      hkns[j]=new TH2F(Form("hkns%s",trig),Form("%s kn and ks before correction;%s kn;%s ks",trig,trig,trig),50,l,h,50,l,h);
      hkn[j]=new TH1F(Form("hkn%s",trig),Form("%s kn before correction;%s kn",trig,trig),50,l,h);
      hks[j]=new TH1F(Form("hks%s",trig),Form("%s ks before correction;%s ks",trig,trig),50,l,h);
      hrns[j]=new TH1F(Form("hrns%s",trig),Form("%s coin rate before correction;%s R_{NS}",trig,trig),50,fl,fh);
      hmuns[j]=new TH1F(Form("hmuns%s",trig),Form("%s underlying collision rate (0th order);%s #mu_{NS}",trig,trig),50,rl,rh);
      hmurns[j]=new TH2F(Form("hmurns%s",trig),Form("%s underlying collision rate (0th order) vs coin rate;%s R_{NS};%s #mu_{NS}",trig,trig,trig),50,fl,fh,50,fl,fh);
      hmumu[j]=new TH2F(Form("hmumu%s",trig),Form("%s underlying collision rate (0th order) vs precalc mu;calc'd %s mu_{NS};%s #mu_{NS}",trig,trig,trig),50,fl,fh,50,fl,fh);

      // if (j==0)t->Draw(Form("%s:%s:%s",ch_kn_bbc,ch_ks_bbc,ch_mu_bbc),"1","goff");
      if (j==0)t->Draw(Form("%s:%s:star_bbcwide:%s",ch_kn_bbc,ch_ks_bbc,ch_mu_bbc),"1","goff");
      if (j==1)t->Draw(Form("%s:%s:star_zdcwide:%s",ch_kn_zdc,ch_ks_zdc,ch_mu_zdc),"1","goff");

      c->cd(1+j*nplots);
      for (int i=0;i<t->GetSelectedRows();i++){
	double trns=t->GetV3()[i];
	double tkn=t->GetV1()[i];
	double tks=t->GetV2()[i];
	double tmuns=t->GetV4()[i];
	//printf("pars:%f,%f,%f\n",t->GetV3()[i],t->GetV2()[i],t->GetV1()[i]);
	globwrap->f->SetParameters(trns,tkn,tks);//t->GetV3()[i],t->GetV2()[i],t->GetV1()[i]);
	globwrap->d->SetParameters(trns,tkn,tks);//t->GetV3()[i],t->GetV2()[i],t->GetV1()[i]);
	if (i==0) {globwrap->f->DrawCopy();newline.SetLineColor(kBlack);newline.DrawLine(0.0,0.0,1.0,0.0);
	} else {if (i%49==0) globwrap->f->DrawCopy("same");}
	rf4->SetFunction(gfunc,1);
	bool ret=rf4->Solve();
	//printf("return code=%d\n",ret);
	float root=rf4->Root();
	if (root<0.01 || root>1.9) printf("root=%f\n",root);
	hkns[j]->Fill(tkn,tks);
	hkn[j]->Fill(tkn);
	hks[j]->Fill(tks);
	hrns[j]->Fill(trns);
	hmuns[j]->Fill(root);
	hmurns[j]->Fill(trns,root);
	hmumu[j]->Fill(tmuns,root);
      }
      c->cd(2+j*nplots);
      hkns[j]->Draw("colz");
      c->cd(3+j*nplots);
      hrns[j]->Draw();
      c->cd(4+j*nplots);
      hmurns[j]->Draw("colz");
      c->cd(5+j*nplots);
      hmumu[j]->Draw("colz");
      newline.DrawLine(fl,fl,fh,fh);
    }

    //t->Draw(Form("%s:%s:star_zdcwide",ch_kn_zdc),"1");
    nc++;
  }





  //iteratively compute the fill-by-fill extrapolations to zero rate in order to remove rate-dependent effects from the true kn and ks values:
  const char* ch_kn_bbc_corr="(log(1-star_bscorrs)-log(1-star_bncorrs-star_bscorrs+star_bbcwidecorrs))/(log(1-star_bncorrs-star_bscorrs+star_bbcwidecorrs)-log(1-star_bscorrs)-log(1-star_bncorrs))";
  const char* ch_ks_bbc_corr="(log(1-star_bncorrs)-log(1-star_bncorrs-star_bscorrs+star_bbcwidecorrs))/(log(1-star_bncorrs-star_bscorrs+star_bbcwidecorrs)-log(1-star_bscorrs)-log(1-star_bncorrs))";
  const char* ch_kn_zdc_corr="(log(1-star_zscorrs)-log(1-star_zncorrs-star_zscorrs+star_zdcwidecorrs))/(log(1-star_zncorrs-star_zscorrs+star_zdcwidecorrs)-log(1-star_zscorrs)-log(1-star_zncorrs))";
  const char* ch_ks_zdc_corr="(log(1-star_zncorrs)-log(1-star_zncorrs-star_zscorrs+star_zdcwidecorrs))/(log(1-star_zncorrs-star_zscorrs+star_zdcwidecorrs)-log(1-star_zscorrs)-log(1-star_zncorrs))";
  const char* ch_mu_bbc_corr="(log(1-star_bncorrs-star_bscorrs+star_bbcwidecorrs)-log(1-star_bscorrs)-log(1-star_bncorrs))";
  const char* ch_mu_zdc_corr="(log(1-star_zncorrs-star_zscorrs+star_zdcwidecorrs)-log(1-star_zscorrs)-log(1-star_zncorrs))";

   
  if (0){//do the iterated fits and save the uLumi tree
    TF1 *rateobs=new TF1("rateobs","-[0]+1-exp(-([1]+1)*x)-exp(-([2]+1)*x)+exp(-([1]+[2]+1)*x)",-0.01,1);
    TF1 *ratederivative=new TF1("ratederivative",
				"([1]+1)*exp(-([1]+1)*x)+([2]+1)*exp(-([2]+1)*x)-([1]+[2]+1)*exp(-([1]+[2]+1)*x)",-0.01,1);
    rateobs->SetTitle("Random Subset of 0=Robs-f(kn,ks,#mu);#mu;'0'");
    globwrap->SetF(rateobs);
    globwrap->SetD(ratederivative);
    //ROOT::Math::RootFinder *rf4 = new ROOT::Math::RootFinder(ROOT::Math::RootFinder::kGSL_STEFFENSON);
    //ROOT::Math::GradFunctor1D gfunc( &global_eval, &global_derive);
    const int nIterations=10;
    vector<double>zmu[nIterations],bmu[nIterations],zkn,zks,bkn,bks;
    vector<int>fill,run,bpat,ypat;
    vector<double>zdcr,bbcr;
    vector<int>fillstart,fillend;
    vector<int>runstart,runend;

    vector<double> bpol, bpol_err, ypol, ypol_err, clk;//, pedzdcks, pedbbckn, pedbbcks;//pedro's corrected values, I think?
    vector<double> fillpedzkn,fillpedzks,fillpedbkn,fillpedbks;//the fill-by-fill values of kn and ks.
   
    int length;
    int tfill=0,trun=0,lastfill=0,lastrun=0;

    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1200,800);
    int nplots=4;
    c->Divide(4,2);


    //load the arrays we will need for the iterative kn procedure:
    t->Draw(Form("star_fill:star_run:%s:%s:%s:%s:%s:%s:star_zdcwide:star_bbcwide:bpat:ypat:gl1p_bpol:gl1p_bpolerr:gl1p_ypol:gl1p_ypolerr:star_clk",
		 ch_kn_bbc,ch_ks_bbc,ch_mu_bbc,
		 ch_kn_zdc,ch_ks_zdc,ch_mu_zdc),"1","goff");
    length=t->GetSelectedRows();
    for (int i=0;i<length;i++){
       
      fill.push_back(t->GetVal(0)[i]);
      run.push_back(t->GetVal(1)[i]);
      bkn.push_back(t->GetVal(2)[i]);
      bks.push_back(t->GetVal(3)[i]);
      bmu[0].push_back(t->GetVal(4)[i]);

      zkn.push_back(t->GetVal(5)[i]);
      zks.push_back(t->GetVal(6)[i]);
      zmu[0].push_back(t->GetVal(7)[i]);
      zdcr.push_back(t->GetVal(8)[i]);
      bbcr.push_back(t->GetVal(9)[i]);
      bpat.push_back(t->GetVal(10)[i]);
      ypat.push_back(t->GetVal(11)[i]);
      bpol.push_back(t->GetVal(12)[i]);
      bpol_err.push_back(t->GetVal(13)[i]);
      ypol.push_back(t->GetVal(14)[i]);
      ypol_err.push_back(t->GetVal(15)[i]);
      clk.push_back(t->GetVal(16)[i]);


      for (int j=1;j<nIterations;j++){
	//the Iteration for kn wants to interact with the corrected mu vectors like arrays, so have to pre-load them with enough elements that they dont' explode:
	bmu[j].push_back(0);
	zmu[j].push_back(0);
      }
       

      //printf("fill %d run=%d\n",fill[i],run[i]);
      if (i==0 ||fill[i]!=fill[i-1]){
	fillstart.push_back(i);
      }
      if (i!=0 &&fill[i]!=fill[i-1]){
	fillend.push_back(i);
	//printf("fill %d %d<=i<%d\n",fill[i-1],fillstart[fillstart.size()-2],fillend[fillend.size()-1]);
	//printf("zmu=%f, zkn=%f\n",zmu[0][i],zkn[i]);
      }
      if (i==0 ||run[i]!=run[i-1]){
	runstart.push_back(i);
      }
      if (i!=0 &&run[i]!=run[i-1]){
	runend.push_back(i);
      }
    }
    fillend.push_back(length);
    runend.push_back(length);

    const int nRuns=runstart.size();
    const int nFills=fillstart.size();
    vector<double>zkn0[nFills],zks0[nFills],bkn0[nFills],bks0[nFills];//the extrapolations to zero, one per run (should be one per /fill/, but we can check run by run...)

    //so we can plot evolution of parameters against iteration number:
    vector<double>dummyindex;
    for (int i=0;i<nIterations;i++){
      dummyindex.push_back(i*1.0);
    }

     

    TGraph *gt;
    TF1 *linear=new TF1("lineark","[0]+[1]*x",0,1);

    c->cd(1);
    gt=new TGraph(fillend[0],&(zmu[0][0]),&(zkn[0]));
    //gt=new TGraph(3,test,test1);
    gt->SetTitle("ZDC north exclusive 1-2 ratio vs corrected rate (0th order);ZDC #mu;ZDC kn");
    gt->Draw("*A");

    c->cd(2);
    gt=new TGraph(fillend[0],&(zmu[0][0]),&(zks[0]));
    gt->SetTitle("ZDC south exclusive 1-2 ratio vs corrected rate (0th order);ZDC #mu;ZDC ks");
    gt->Draw("*A");
    /*
      c->cd(3);
      gt=new TGraph(length,&(bmu[0][0]),&(bkn[0]));
      gt->SetTitle("BBC north exclusive 1-2 ratio vs corrected rate (0th order);BBC #mu;BBC kn");
      gt->Draw("*A");
     
      c->cd(4);
      gt=new TGraph(length,&(bmu[0][0]),&(bks[0]));
      gt->SetTitle("BBC south exclusive 1-2 ratio vs corrected rate (0th order);BBC #mu;BBC ks");
      gt->Draw("*A");
    */
    int testfilli=0;
    for (int i=0;i<fillstart.size();i++){
      //just compare to Pedro.
      //if (fill[fillstart[i]]!=17318) continue;
      printf("found fill 17318 @ i=%d\n",i);
      testfilli=i;
       
      //prime the iteration with a 0th order guess of the extrapolated-to-zero exclusive ratios:
      gt=new TGraph(fillend[i]-fillstart[i],&(zmu[0][fillstart[i]]),&(zkn[fillstart[i]]));
      gt->SetMarkerColor(i%6+2);
      gt->Fit(linear);
      zkn0[i].push_back(linear->GetParameter(0));
      c->cd(1);
      if (!(i%10))gt->Draw("*");
   
      gt=new TGraph(fillend[i]-fillstart[i],&(zmu[0][fillstart[i]]),&(zks[fillstart[i]]));
      gt->SetMarkerColor(i%6+2);
      gt->Fit(linear);
      zks0[i].push_back(linear->GetParameter(0));
      c->cd(2);
      if (!(i%10))gt->Draw("*");
      c->cd(3);

      //and repeat for the BBC:
      gt=new TGraph(fillend[i]-fillstart[i],&(bmu[0][fillstart[i]]),&(bkn[fillstart[i]]));
      gt->Fit(linear);
      bkn0[i].push_back(linear->GetParameter(0));
      gt=new TGraph(fillend[i]-fillstart[i],&(bmu[0][fillstart[i]]),&(bks[fillstart[i]]));
      gt->Fit(linear);
      bks0[i].push_back(linear->GetParameter(0));

       

      //now iterate:
      double newks,newkn;//temporary holders to get the new kn0 and ks0 out of the iteration function.

      //iteratively correct the BBC:
      for (int j=1;j<nIterations;j++){
	//bool IterateRateCorrection(int length, double *rate_arr,
	//                           double kn0, double ks0, double *kn_arr, double *ks_arr,
	//                           double *new_kn0, double *new_ks0, double *new_mu_arr)

	bool result=IterateRateCorrection(fillend[i]-fillstart[i], &(bbcr[fillstart[i]]),
					  bkn0[i][j-1], bks0[i][j-1], &(bkn[fillstart[i]]), &(bks[fillstart[i]]),
					  &newkn,&newks,&(bmu[j][fillstart[i]]));
	bkn0[i].push_back(newkn);
	bks0[i].push_back(newks);
      }

      //repeat the correction for zdc
      for (int j=1;j<nIterations;j++){
	//bool IterateRateCorrection(int length, double *rate_arr,
	//                           double kn0, double ks0, double *kn_arr, double *ks_arr,
	//                           double *new_kn0, double *new_ks0, double *new_mu_arr)

	bool result=IterateRateCorrection(fillend[i]-fillstart[i], &(zdcr[fillstart[i]]),
					  zkn0[i][j-1], zks0[i][j-1], &(zkn[fillstart[i]]), &(zks[fillstart[i]]),
					  &newkn,&newks,&(zmu[j][fillstart[i]]));
	zkn0[i].push_back(newkn);
	zks0[i].push_back(newks);
	if (zmu[j][fillstart[i]]<0.0002){
	  printf("run %d has a problem in iteration %d\n",run[fillstart[i]],j);
	  c->cd(3);	   
	  globwrap->f->Draw();
	  return;
	     
	  assert(false);
	}
      }
       
    }
    c->cd(4);
    gt=new TGraph(fillend[testfilli]-fillstart[testfilli],
		  &(zmu[0][fillstart[testfilli]]),&(zkn[fillstart[testfilli]]));
    //gt=new TGraph(3,test,test1);
    gt->SetTitle("ZDC north exclusive 1-2 ratio vs corrected rate (0th order);ZDC #mu;ZDC kn");
    gt->Draw("*A");
    gt=new TGraph(fillend[testfilli]-fillstart[testfilli],&(zmu[1][fillstart[testfilli]]),&(zkn[fillstart[testfilli]]));
    gt->SetMarkerColor(kGreen);
    gt->Draw("*");


    c->cd(5);
    gt=new TGraph(nIterations,&(dummyindex[0]),&(zkn0[testfilli][0]));
    gt->SetTitle("ZDC north exclusive ratio extrapolated to zero;iteration;ZDC kn(0)");
    gt->Draw("*A");
    printf("final order intercept: ZDC kn0=%f\n",zkn0[testfilli][nIterations-1]);

    c->cd(6);
    gt=new TGraph(nIterations,&(dummyindex[0]),&(zks0[testfilli][0]));
    gt->SetTitle("ZDC south exclusive ratio extrapolated to zero;iteration;ZDC kn(0)");
    gt->Draw("*A");
    printf("final order intercept: ZDC ks0=%f\n",zks0[testfilli][nIterations-1]);

    //compare to pedro's numbers:?
    //compare to each other:
    c->cd(6);
    gt=new TGraph(zmu[nIterations-1].size(),&(zmu[0][0]),&(zmu[nIterations-1][0]));
    gt->SetTitle("ZDC underlying rate before and after coorection;before;after");
    gt->Draw("*A");
    c->cd(7);
    gt=new TGraph(bmu[nIterations-1].size(),&(bmu[0][0]),&(bmu[nIterations-1][0]));
    gt->SetTitle("BBC underlying rate before and after coorection;before;after");
    gt->Draw("*A");
    c->cd(8);
    gt=new TGraph(zmu[nIterations-1].size(),&(zmu[0][0]),&(bmu[0][0]));
    gt->SetTitle("uncorrected underlying rates ZDC vs BBC;ZDC mu;BBC mu");
    gt->Draw("*A");
    c->cd(8);
    gt=new TGraph(zmu[nIterations-1].size(),&(zmu[nIterations-1][0]),&(bmu[nIterations-1][0]));
    gt->SetTitle("corrected underlying rates ZDC vs BBC;ZDC mu;BBC mu");
    gt->Draw("*A");

    vector<double>ratio[2];
    for (int i=0;i<bmu[nIterations-1].size();i++){
      ratio[0].push_back(zmu[nIterations-1][i]/zmu[0][i]);
      ratio[1].push_back(bmu[nIterations-1][i]/bmu[0][i]);
    }
    c->cd(6);
    gt=new TGraph(bmu[nIterations-1].size(),&(zmu[0][0]),&(ratio[0][0]));
    gt->SetTitle("ratio of corrected to uncorrected ZDC mu;ZDC mu;ratio");
    gt->Draw("*A");
    c->cd(7);
    gt=new TGraph(bmu[nIterations-1].size(),&(bmu[0][0]),&(ratio[1][0]));
    gt->SetTitle("ratio of corrected to uncorrected BBC mu;BBC mu;ratio");
    gt->Draw("*A");

    nc++;

    //plot the distribution before and after the corrections

    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),1000,400);
    c->Divide(4,1);
     
    TH1F *hbeforecor[4];
    TH1F *haftercor[4];
    for (int i=0;i<4;i++){
      hbeforecor[i]=new TH1F(Form("hbeforecorr%d",i),Form("k%d Distribution before correction",i),200,0,5);
      haftercor[i]=new TH1F(Form("haftercor%d",i),Form("k%d Distribution after correction",i),200,0,5);
    }
    for (int i=0;i<fillstart.size();i++){
      hbeforecor[0]->Fill(bks0[i][0]);
      haftercor[0]->Fill(bks0[i][nIterations-1]);
      hbeforecor[1]->Fill(bkn0[i][0]);
      haftercor[1]->Fill(bkn0[i][nIterations-1]);
      hbeforecor[2]->Fill(zks0[i][0]);
      haftercor[2]->Fill(zks0[i][nIterations-1]);
      hbeforecor[3]->Fill(zkn0[i][0]);
      haftercor[3]->Fill(zkn0[i][nIterations-1]);
    }
    for (int i=0;i<4;i++){
      c->cd(i+1);
      hbeforecor[i]->Draw();
      haftercor[i]->SetLineColor(kRed);
      haftercor[i]->Draw("same");
    }
       
    //now we can plot the ZDC/BBC asymmetry:
    //ALL=1/pp*(zbpp-zbpm)/(zbpp+zbpm)
    TTree *uLumi=new TTree("uLumi","micro luminosity tree");
    int u_fill; uLumi->Branch("fill",&u_fill);
    int u_run; uLumi->Branch("run",&u_run);
    int u_pat; uLumi->Branch("pat",&u_pat);
    double u_bpol; uLumi->Branch("bpol",&u_bpol);
    double u_ypol; uLumi->Branch("ypol",&u_ypol);
    double u_bpol_err; uLumi->Branch("bpol_err",&u_bpol_err);
    double u_ypol_err; uLumi->Branch("ypol_err",&u_ypol_err);

    double u_likemuz0; uLumi->Branch("likemuz0",&u_likemuz0);
    double u_likemuz1; uLumi->Branch("likemuz1",&u_likemuz1);
    double u_likemub0; uLumi->Branch("likemub0",&u_likemub0);
    double u_likemub1; uLumi->Branch("likemub1",&u_likemub1);
    double u_unlikemuz0; uLumi->Branch("unlikemuz0",&u_unlikemuz0);
    double u_unlikemuz1; uLumi->Branch("unlikemuz1",&u_unlikemuz1);
    double u_unlikemub0; uLumi->Branch("unlikemub0",&u_unlikemub0);
    double u_unlikemub1; uLumi->Branch("unlikemub1",&u_unlikemub1);
    double u_rellumizdc; uLumi->Branch("rellumizdc",&u_rellumizdc);
    double u_rellumizdc_err; uLumi->Branch("rellumizdc_err",&u_rellumizdc_err);
    double u_rellumibbc; uLumi->Branch("rellumibbc",&u_rellumibbc);
    double u_rellumibbc_err; uLumi->Branch("rellumibbc_err",&u_rellumibbc_err);
    double u_allzdcbbc; uLumi->Branch("allzdcbbc",&u_allzdcbbc);
    double u_allzdcbbc_err; uLumi->Branch("allzdcbbc_err",&u_allzdcbbc_err);

    vector<double> uLumi_mon_run;
    vector<double> uLumi_mon_all;
    vector<double> uLumi_mon_all_check;
    vector<double> uLumi_mon_all_check2;

    vector<double> uLumi_mon_run_bar;
    vector<double> uLumi_mon_all_err;

    //loop over the arrays run by run to compose the u_(un)like* stuff:
    for (int i=0;i<nRuns;i++){
      //find the start of our fill:
      u_fill=fill[runstart[i]];
      int filli=0;
      for (filli=0;filli<fillstart.size();filli++){
	if (fill[fillstart[filli]]==u_fill) {
	  printf("fill %d, run%d",u_fill,u_run);
	  break;
	}
      }
       
      //run-by-run info:
      u_run=run[runstart[i]];
       
      u_pat=RetrieveSpinPattern(u_run);
      u_bpol=bpol[runstart[i]];
      u_ypol=ypol[runstart[i]];
      u_bpol_err=bpol_err[runstart[i]];
      u_ypol_err=ypol_err[runstart[i]];
      if (u_bpol==0 || u_ypol==0) {
	printf("run %d has bpol=%f, ypol=%f, skipping\n",u_run,u_bpol,u_ypol);
	continue;
      }
      u_likemuz0=0;
      u_likemuz1=0;
      u_likemub0=0;
      u_likemub1=0;
      u_unlikemuz0=0;
      u_unlikemuz1=0;
      u_unlikemub0=0;
      u_unlikemub1=0;
      for (int j=runstart[i];j<runend[i];j++){
	//bunch-by-bunch summation:
	if (bpat[j]==ypat[j]){
	  u_likemuz0+=zmu[0][j]*clk[j];//[filli];//*clk[j];
	  u_likemuz1+=zmu[nIterations-1][j]*clk[j];
	  u_likemub0+=bmu[0][j]*clk[j];
	  u_likemub1+=bmu[nIterations-1][j]*clk[j];
	} else {
	  u_unlikemuz0+=zmu[0][j]*clk[j];
	  u_unlikemuz1+=zmu[nIterations-1][j]*clk[j];
	  u_unlikemub0+=bmu[0][j]*clk[j];
	  u_unlikemub1+=bmu[nIterations-1][j]*clk[j];
	}
      }
      printf("(zl=%2.2E zu=%2.2E diff=%2.2E) (bl=%2.2E bu=%2.2E diff=%2.2E) ",
	     u_likemuz1,u_unlikemuz1,u_likemuz1-u_unlikemuz1,
	     u_likemub1,u_unlikemub1, u_likemub1-u_unlikemub1);
      // printf("(zl=%2.2E zu=%2.2E) (bl=%2.2E bu=%2.2E) ",u_likemuz1,u_unlikemuz1,u_likemub1,u_unlikemub1);

 
      uLumi_mon_run.push_back(u_run);
      /*
	uLumi_mon_all.push_back(1/(u_bpol*u_ypol)*
	(u_likemuz1/u_likemub1-u_unlikemuz1/u_unlikemub1)
	/(u_likemuz1/u_likemub1+u_unlikemuz1/u_unlikemub1));
      */
      //float qqmuz1,qqmub1,qqmuz0,qqmub0,qpmuz1,qpmub1,qpmuz0,qpmub0;
	 
      uLumi_mon_all.push_back(1/((float)u_bpol*(float)u_ypol)*
			      ((float)u_likemuz1/(float)u_likemub1-(float)u_unlikemuz1/(float)u_unlikemub1)
			      /((float)u_likemuz1/(float)u_likemub1+(float)u_unlikemuz1/(float)u_unlikemub1));

      double asym,asym_err;
      CalcAsymAndErr(&asym, &asym_err,
		     u_bpol,u_bpol_err,
		     u_ypol,u_ypol_err,
		     u_unlikemuz1,sqrt(u_unlikemuz1),
		     u_likemuz1,sqrt(u_likemuz1),
		     u_unlikemub1,sqrt(u_unlikemub1),
		     u_likemub1,sqrt(u_likemub1));
      uLumi_mon_all_check2.push_back(asym);
   
      u_rellumizdc=u_likemuz1/u_unlikemuz1;
      u_rellumizdc_err=u_rellumizdc*sqrt(1/u_likemuz1+1/u_unlikemuz1);
      u_rellumibbc=u_likemub1/u_unlikemub1;
      u_rellumibbc_err=u_rellumibbc*sqrt(1/u_likemub1+1/u_unlikemub1);
      u_allzdcbbc=asym;
      u_allzdcbbc_err=asym_err;

      uLumi->Fill();
      //translated from rcc_calc_all:
      double tunlike=u_unlikemuz1;
      double tlike=u_likemuz1;
      double trellumi=u_likemub1/u_unlikemub1;
      double trellumi_err=sqrt(trellumi*trellumi*(1/u_likemub1+1/u_unlikemub1));//assumes these are counting.
      double trelunlike=trellumi*tunlike;
      double tsum=tlike+trelunlike;
      double tdiff=tlike-trelunlike;
      double tpolfactor=1/(u_bpol*u_ypol);
      double tasym=tpolfactor*(tdiff/tsum);

      double tbpolerrterm=-tasym/u_bpol * u_bpol_err;
      double typolerrterm=-tasym/u_ypol * u_ypol_err;

      double tlikesumerrterm=tpolfactor*(2*trelunlike)/(tsum*tsum)*sqrt(tlike);//assumes tlike is a counting exercise: sqrt(counts)
      double tunlikesumerrterm=-tpolfactor*(2*trelunlike*trellumi)/(tsum*tsum)*sqrt(tunlike);//assumes tlike is a counting exercise
      double trellumierrterm=-tpolfactor*(2*trelunlike*tunlike)/(tsum*tsum)*trellumi_err;

      double err2=(tbpolerrterm*tbpolerrterm
		   +typolerrterm*typolerrterm
		   +tlikesumerrterm*tlikesumerrterm
		   +tunlikesumerrterm*tunlikesumerrterm
		   +trellumierrterm*trellumierrterm);

      uLumi_mon_all_check.push_back(tasym);
      uLumi_mon_all_err.push_back(sqrt(err2));
      uLumi_mon_run_bar.push_back(0);


      printf(" all=%2.2E vs tasym=%2.2E vs Get()=%2.2E (zl=%2.2E zu=%2.2E) (bl=%2.2E bu=%2.2E)\n",uLumi_mon_all[uLumi_mon_all.size()-1],tasym, asym,u_likemuz1,u_unlikemuz1,u_likemub1,u_unlikemub1);
      //printf(" all=%2.2E \n",uLumi_mon_all[uLumi_mon_all.size()-1]);

       
    }
    uLumi->SaveAs("uLumi.ttree.root");
    nc++;

    //plot the distribution before and after the corrections
     
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,400);
    c->Divide(1,3);
    c->cd(1);
    gt=new TGraphErrors(uLumi_mon_run.size(),
			&(uLumi_mon_run[0]),&(uLumi_mon_all[0]),
			&(uLumi_mon_run_bar[0]),&(uLumi_mon_all_err[0]));
    gt->SetTitle("Rel Lumi ZDC/BBC asymmetry vs run;run;asym");
    gt->Draw("A*");
    c->cd(2);
    gt=new TGraphErrors(uLumi_mon_run.size(),
			&(uLumi_mon_run[0]),&(uLumi_mon_all_check[0]),
			&(uLumi_mon_run_bar[0]),&(uLumi_mon_all_err[0]));
    gt->SetTitle("Rel Lumi ZDC/BBC asymmetry vs run;run;asym");
    gt->Draw("A*");
    c->cd(3);
    gt=new TGraph(uLumi_mon_run.size(),
		  &(uLumi_mon_all_check[0]),&(uLumi_mon_all[0]));
    gt->SetTitle("Cross check of asym calculation;asym (chk);asym");
    gt->Draw("A*");
  }

  if (1){//bunch by bunch status and spin bits

    t->Draw(Form("star_fill:star_run:crossing:bpat:ypat"),"1","goff");
    int length=t->GetSelectedRows();
    TTree *uBunch=new TTree("uBunch","micro bunch bits tree");
    int u_fill; uBunch->Branch("fill",&u_fill);
    int u_run; uBunch->Branch("run",&u_run);
    int u_pat; uBunch->Branch("pat",&u_pat);
    int u_bunch; uBunch->Branch("bunch",&u_bunch);
    int u_bspin; uBunch->Branch("bspin",&u_bspin);
    int u_yspin; uBunch->Branch("yspin",&u_yspin);
    for (int i=0;i<length;i++){
      u_fill=t->GetVal(0)[i];
      u_run=t->GetVal(1)[i];
      u_pat=RetrieveSpinPattern(u_run);//redundant since it doesn't change with bunchID, but convenient.
      u_bunch=t->GetVal(2)[i];
      u_bspin=t->GetVal(3)[i];
      u_yspin=t->GetVal(4)[i];
      if (abs(u_bspin*u_yspin)!=1) continue; //remove anything without spin info still.
      uBunch->Fill();
    }
    uBunch->SaveAs("uBunch.ttree.root");

  

    
  
    //eventually, I need to add bunch-by-bunch luminosity counters, but skip for now.
  }





   


  
  return;
  printf("Caution!  Below this point are old Draw commands that assume we're running with the old starscaler-only files.  The names of variables will not match the combined db file.\n");
  return;
  
  //display rates vs bunch and fill
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:fill>>hnew1(120,17180,17620,120,-0.5,119.5)","log10(rclk)","colz");
    c->cd(2);
    t->Draw("cross:fill>>hnew2(120,17180,17620,120,-0.5,119.5)","log10(clk)","colz");
    c->SetLogz();
    c->cd(3);
    t->Draw("cross:fill>>hnew3(120,17180,17620,120,-0.5,119.5)","bbcwide","colz");
    c->SetLogz();
    c->cd(4);
    t->Draw("cross:fill","bbncnt/bbcwidecnt","colz");
    c->SetLogz();
    nc++;
  }

  //check some bunch crossing-related anomalies
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,900);
    c->Divide(2,3);
    c->cd(1);
    t->Draw("cross:fill","bbcwide *(cross>20 && cross<40)","colz");
    c->cd(2);
    t->Draw("cross:fill","bbcwide *(cross>60 && cross<80)","colz");
    c->cd(3);
    t->Draw("cross","bbcwide * (fill<17210 && cross<10)","colz");
    c->cd(4);
    t->Draw("cross","bbcwide * (fill>17210 && cross<10)","colz");
    c->cd(5);
    t->Draw("cross:fill","bbcwide *(cross>100 && fill>17430 && fill<17450)","colz");
    c->cd(6);
    t->Draw("(cross>110):bbcwide>>hnew(10,-0.05,0.95,2,-0.5,1.5)","1","colz");
    nc++;
  }

  //detailed check of those anomalies
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,900);
    c->Divide(2,3);
    c->cd(1);
    t->Draw("bbcwide","((cross==29 || cross==30 || cross==69 || cross==70) && fill>17380 && fill<17420)","colz");
    c->cd(2);
    t->Draw("bbcwide>0.1:run-386700","((cross==29 || cross==30 || cross==69 || cross==70) && fill>17380 && fill<17420)","colz");
    c->cd(3);
    t->Draw("cross:run-386700","bbcwide * (fill<17212 && cross<5)","colz");
    c->cd(4);
    t->Draw("fill:run-386700","(fill<17212 && cross==0)","colz");
    c->cd(5);
    t->Draw("cross:run-386700","bbcwide *(cross>100 && fill==17443)","colz");
    c->cd(6);
    t->Draw("bbcwide:run-386700","(cross>110 && fill>17435 && fill<17450)","colz");
    nc++;
  }

  //display rates vs bunch and fill after cuts designed to clean them.
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:fill>>hnew1(120,17180,17620,120,-0.5,119.5)","log10(rclk)"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->cd(2);
    t->Draw("cross:fill>>hnew2(120,17180,17620,120,-0.5,119.5)","log10(clk)"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->SetLogz();
    c->cd(3);
    t->Draw("cross:fill>>hnew3(120,17180,17620,120,-0.5,119.5)","bbcwide"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->SetLogz();
    c->cd(4);
    t->Draw("cross:fill","bbncnt/bbcwidecnt"*rcc_cross_qa+rcc_clip_loud_runs,"colz");
    c->SetLogz();
    nc++;
  }

  //a marginally useful look at the abort gap spectrum after the rcc_cross_qa cuts
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("run-386700","cross>110"+rcc_cross_qa,"colz");
    c->cd(2);
    t->Draw("bbcwide","cross>110"+rcc_cross_qa,"colz");
    c->cd(3);
    t->Draw("run-386700","cross==29 || cross==30 || cross==69 || cross==70"+rcc_cross_qa,"colz");
    c->cd(4);
    t->Draw("bbcwide","cross==29 || cross==30 || cross==69 || cross==70"+rcc_cross_qa,"colz");
    nc++;
  }

  //a better look at the abort gap spectrum after the rcc_cross_qa cuts
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("run-386700","cross>110"+rcc_cross_qa,"colz");
    c->cd(2);
    t->Draw("bbcwide","cross>110"+rcc_cross_qa,"colz");
    c->cd(3);
    t->Draw("run-386700:fill>>hnew(350,17299.5,17649.5,100,2000,12000)","cross>110"+rcc_cross_qa,"colz");
    c->cd(4);
    t->Draw("run-386700:fill>>hnew2(350,17299.5,17649.5,100,2000,12000)","cross>110","colz");
    nc++;
  }
  
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
  if (0){
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
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(2,2);
    c->cd(1);
    t->Draw("cross:zdcwidecnt","1","colz");
    t->SetMarkerColor(kRed);
    c->cd(2);
    t->Draw("cross:zdcwidecnt",minclocks && live70 && minbbcrate && bbsslope,"colz");
    t->SetMarkerColor(kBlack);   
    c->cd(3);
    t->Draw("cross:zdcwidecnt",abortgap,"colz");
    t->SetMarkerColor(kBlack);   
    c->cd(4);
    t->Draw("cross:zdcwidecnt",abortgap&&minclocks && live70 && minbbcrate && bbsslope,"colz");
    t->SetMarkerColor(kBlack);   
    nc++;
  }

  //after all cuts, zdc to bbc ratio, to look for excursions
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(3,3);

    //ratio vs run
    c->cd(1);
    t->Draw("zdcwidecnt/bbcwidecnt:run","1","colz");
    t->SetMarkerColor(kRed);
    c->cd(2);
    t->Draw("zdcwidecnt/bbcwidecnt:run","zdcwidecnt/bbcwidecnt<0.45","colz");
    t->SetMarkerColor(kRed);
    c->cd(3);
    t->Draw("zdcwidecnt/bbcwidecnt:run",allcuts,"colz");
    t->SetMarkerColor(kBlack);

    //ratios vs one part
    c->cd(4);
    t->Draw("zdcwidecnt/bbcwidecnt:zdcwidecnt",allcuts,"colz");
    c->cd(5);
    t->Draw("zdcwidecnt:bbcwidecnt",allcuts,"colz");
    c->cd(6);
    t->Draw("zdcwidecnt/bbcwidecnt:bbcwidecnt",allcuts,"colz");  
    //ratios vs one part, normalized by live clocks
    c->cd(7);
    t->Draw("zdcwidecnt/bbcwidecnt:zdcwidecnt/rclk",allcuts,"colz");
    c->cd(8);
    t->Draw("zdcwidecnt/rclk:bbcwidecnt/rclk",allcuts,"colz");
    c->cd(9);
    t->Draw("zdcwidecnt/bbcwidecnt:bbcwidecnt/rclk",allcuts,"colz");  
    nc++;
  }

  //after all cuts, inspect the corrected ratios and see how smooth they're getting...
  if (0){
    c=new TCanvas(Form("c%d",nc),Form("c%d",nc),800,600);
    c->Divide(1,3);

    //ratios vs one part, normalized by live clocks
    c->cd(1);
    t->Draw("zdcwidecnt/bbcwidecnt:zdcwidecnt/clk",allcuts);
    c->cd(2);
    t->Draw("zdcwidecnt/rclk:bbcwidecnt/clk",allcuts);
    c->cd(3);
    t->Draw("zdcwidecnt/bbcwidecnt:bbcwidecnt/clk",allcuts);  
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

unsigned long ToBitPattern(vector<int> vec, int nbits){
  //convert a vector to an unsigned long/
  unsigned long pattern;
  for (int bit=0;bit<nbits;bit+=2){
    if ( vec[bit]!=vec[bit+1]) {
      printf("parity wrong in spin pattern\n");
      return 0;
    }
    if (abs(vec[bit])!=1){
      printf("spin pattern bit isn't +/-1.\n");
      return 0;
    }
    //make the 'even bits only' object:
    //note that to make this human-readable, the leftmost bit is the first bit, not the last one
    pattern+=(vec[bit]==1)<<((nbits-bit)/2-1);
  }
  return pattern;
}

int DetermineSpinPattern(vector<int> bpat, vector<int> ypat){
  //convert the vector of patterns into the bit string representation:
  unsigned long bshort=ToBitPattern(bpat,16);
  unsigned long yshort=ToBitPattern(ypat,16);

  //note that in these printouts, we write the leftmost bit first, not last.
  if (0){
    printf("blue pattern=%lu =0b",bshort);
    for (int i=0;i<8;i++){
      printf("%d",(bshort&(1<<(7-i)))?1:0);
    }
    printf("\n");
  }
  if (0){
    printf("yell pattern=%lu =0b",yshort);
    for (int i=0;i<8;i++){
      printf("%d",(yshort&(1<<(7-i)))?1:0);
    }
    printf("\n");
  }

  //define the possible valid patterns
  //see tables 5 and 6 of AN1125 rev01:
  //these read left to right, so the lowest bit is actually the largest crossing number.
  unsigned long base[6];
  base[0]=0b10101010;//as in P1 or P21 blue
  base[1]=0b01010101;//as in P2 or P22 blue
  base[2]=0b11001101;//as in P1 or P2 yellow  
  base[3]=0b00110010;//as in P3 or P4 yellow
  base[4]=0b01100110;//as in P21 or P22 yellow, with the initial bin offset by one to match Ross's observations
  base[5]=0b10011001;//as in P23 or P24 yellow, per the above.


  if (0){
    printf("base[2] pattern=%lu =0b",base[2]);
    for (int i=0;i<8;i++){
      printf("%d",(base[2]&(1<<(7-i)))?1:0);
    }
    printf("\n");
  }
  //define the combinations of patterns that are valid spin patterns:
  //see tables 5 and 6 of AN1125 rev01:

  //initialize all lookup patterns to 'invalid' (=0):
  //patlookup[a][b] resolves to the spin pattern ID# for blue=base[a] && yellow=base[b]
  int patlookup[6][6];
  for (int i=0;i<6;i++){
    for (int j=0;j<6;j++){
      patlookup[i][j]=0;
    }
  }
  //set the valid ones to their proper settings.
  patlookup[0][2]=1;
  patlookup[1][2]=2;
  patlookup[0][3]=3;
  patlookup[1][3]=4;
  patlookup[2][0]=5;
  patlookup[2][1]=6;
  patlookup[3][0]=7;
  patlookup[3][1]=8;

  patlookup[0][4]=21;
  patlookup[1][4]=22;
  patlookup[0][5]=23;
  patlookup[1][5]=24;
  patlookup[4][0]=25;
  patlookup[4][1]=26;
  patlookup[5][0]=27;
  patlookup[5][1]=28;
  

  //find which pattern we have:
  int blabel=-1;
  int ylabel=-1;
  for (int i=0;i<6;i++){
    //printf("trying blue base=%lu\n",base[i]);
    if (bshort==base[i]) {
      blabel=i;
      break;
    }
  }
  for (int i=0;i<6;i++){
    //printf("trying yell base=%lu\n",base[i]);
    if (yshort==base[i]) {
      ylabel=i;
      break;
    }
  }
  if (blabel<0){ 
    printf("couldn't find a blue spin pattern.  parity is good, but not balanced.\n");
    return -1;
  }
  if (ylabel<0) {
    printf("couldn't find a yellow spin pattern.  parity is good, but not balanced.\n");
    return -1;
  }
  return patlookup[blabel][ylabel];
}
  




bool IterateRateCorrection(int length, double *rate_arr, double kn0, double ks0, double *kn_arr, double *ks_arr, double *new_kn0, double *new_ks0, double *new_mu_arr){
  //NOTE:  Assumes new_mu_arr has enough memory to hold length*size_of(float)!
  /* for now, assume the globwrap has been set up
     TF1 *rateobs=new TF1("rateobs","-[0]+1-exp(-([1]+1)*x)-exp(-([2]+1)*x)+exp(-([1]+[2]+1)*x)",0,1);
     TF1 *ratederivative=new TF1("ratederivative",
     "([1]+1)*exp(-([1]+1)*x)+([2]+1)*exp(-([2]+1)*x)-([1]+[2]+1)*exp(-([1]+[2]+1)*x)",0,1);
     rateobs->SetTitle("Random Subset of 0=Robs-f(kn,ks,#mu);#mu;'0'");
     globwrap->SetF(rateobs);
     globwrap->SetD(ratederivative);
  */

  
  //for each measured rate in the sample, compute the associated underlying collision rate:     
  //ROOT::Math::RootFinder *rf4 = new ROOT::Math::RootFinder(ROOT::Math::RootFinder::kGSL_STEFFENSON);
  ROOT::Math::RootFinder *rf4 = new ROOT::Math::RootFinder(ROOT::Math::RootFinder::kGSL_BISECTION);
  ROOT::Math::GradFunctor1D gfunc( &global_eval, &global_derive);
  for (int i=0;i<length;i++){
    //printf("i=%d/%d, zmu=%f, zknc0=%f, zksc0=%f\n",i, length, zmu[0][k],zknc0,zksc0);
    globwrap->f->SetParameters(rate_arr[i],kn0,ks0);
    //if (i==0) {globwrap->f->DrawCopy();newline.SetLineColor(kBlack);newline.DrawLine(0.0,0.0,1.0,0.0);
    //} else {if (i%1==0) globwrap->f->DrawCopy("same");}
    // rf4->SetFunction(gfunc,0.5);
    rf4->SetFunction(gfunc,0,1);//for not using the derivative
    bool ret=rf4->Solve(300,1E-6);
    //printf("return code=%d\n",ret);
    if (ret!=1) {
      printf("EEP. return code=%d.  Last guess was %f\n",ret,rf4->Root());
      return ret;
    }
    float root=rf4->Root();
    new_mu_arr[i]=root;
  }

  
  //for each set of measured exclusives, fit to a linear slope with the new corrected mu:
  TGraph *gt;
  TF1 *linear=new TF1("lineark","[0]+[1]*x",0,1);

  gt=new TGraph(length,new_mu_arr,kn_arr);
  gt->Fit(linear);
  *new_kn0=linear->GetParameter(0);
  gt=new TGraph(length,new_mu_arr,ks_arr);
  gt->Fit(linear);
  *new_ks0=linear->GetParameter(0);

  return true;;
}


int RetrieveSpinPattern(int run, TTree *t){
  static bool firstRun=true;
  static vector<int> runnum;
  static vector<int> pattern;

  if (firstRun){
    firstRun=false;
    //build the vectors.
    t->Draw("crossing:star_run:bpat:ypat","1","goff");
    int nbunches=t->GetSelectedRows();
    //find the starts to each run:
    vector<int> runstart;
    int thisrun, lastrun=0;
    for (int i=0;i<nbunches;i++){
      thisrun=t->GetVal(1)[i];
      if (thisrun!=lastrun){
	runstart.push_back(i);
	runnum.push_back(thisrun);
	lastrun=thisrun;
      }
    }
    int nruns=runstart.size();
    vector<int> bpattern[nruns], ypattern[nruns];
    int bp,yp;
    for (int i=0;i<nruns;i++){
      for (int j=0;j<24;j++){
	bp=t->GetVal(2)[runstart[i]+j];
	yp=t->GetVal(3)[runstart[i]+j];
	bpattern[i].push_back(bp);
	ypattern[i].push_back(yp);
      }

      int patternID=DetermineSpinPattern(bpattern[i],ypattern[i]);
      if (patternID<0){
	printf("Pattern not found:\n");
	if(1){//printf blue pattern
	  printf("run %d blue pat: ",runnum[i]);
	  for (int j=0;j<16;j++){
	    printf("%c",bpattern[i][j]==1?'+':bpattern[i][j]==-1?'-':'X');
	  }
	  printf("\n");
	}
	if(1){//printf yellow pattern
	  printf("run %d yell pat: ",runnum[i]);
	  for (int j=0;j<16;j++){
	    printf("%c",ypattern[i][j]==1?'+':ypattern[i][j]==-1?'-':'X');
	  }
	  printf("\n");
	}
      }
      pattern.push_back(patternID);
    }



     
  }

  if (run<0) return -1;

  for (int i=0;i<runnum.size();i++){
    if (run==runnum[i]) return pattern[i];
  }
  printf("asked for pattern from run %d, but no matching run found.\n",run);
  return -1;
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
