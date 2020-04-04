#include <TString.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <SpinDBContent.hh>
#include <SpinDBOutput.hh>

void rcc_calc_all(const int runnumber = 398149,
			 const char * inputdir="/",
                         const char * outputdir="/phenix/spin/spin1/phnxsp01/rosscorliss/trees/") {

  //define our bins and divisions:
    // spin patterns are in order: ++,+-,--,-+


  const int nptbins=10;
  double pt_limits[nptbins+1] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};

 

  //get run cluster yield
  TString yieldfname = inputdir;//"/direct/phenix+u/rosscorliss/pion_ana/output/"
  yieldfname += runnumber;
  yieldfname += ".MPC.yields.rcc.hist.root";
  TFile *yieldfile = NULL;
  yieldfile=TFile::Open(yieldfname);
  if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys())
    return;

  TH2F *hRunYieldByBunchAndPt=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
   
  //load run scalers from the spindb, do some error checks on them
  //handled in the macro that runs this:  gSystem->Load("libuspin.so");
  SpinDBContent spin_cont;
  SpinDBOutput spin_out("phnxrc");
  spin_out.StoreDBContent(runnumber, runnumber);
  spin_out.GetDBContentStore(spin_cont, runnumber);
  int error = spin_cont.GetErrorValue();
  int qa_level = spin_cont.GetQALevel();
  int fill = spin_cont.GetFillNumber();
  int badrun = spin_cont.GetBadRunFlag();
  int cross_shift = spin_cont.GetCrossingShift(); //spindb stores data in the C-AD convention, which doesn't match the PHENIX numbering.  This is the relative offset.


  //for future use, accumulate the luminosity in a 1x2 histogram:
  TH1F *hTotalLumi=new TH1F("hTotalLumi","Unlike (like) zdc-narrow lumi scaler sums;unlike=0,like=1;sum",2,-0.5,1.5);
  TH2F *hPolarizationBySpin[2][2];
  hPolarizationBySpin[0][0]=new TH2F("hPolarizationBySpinNN","Polarization for B-,Y-;B;Y",100,-1,1,100,-1,1);
  hPolarizationBySpin[1][0]=new TH2F("hPolarizationBySpinPN","Polarization for B+,Y-;B;Y",100,-1,1,100,-1,1);
  hPolarizationBySpin[0][1]=new TH2F("hPolarizationBySpinNP","Polarization for B-,Y+;B;Y",100,-1,1,100,-1,1);
  hPolarizationBySpin[1][1]=new TH2F("hPolarizationBySpinPP","Polarization for B+,Y+;B;Y",100,-1,1,100,-1,1);

  // for each bunch, accumulate yields and bunch info into the appropriate spinpattern grouping:
  TH1F *hYieldByPtAndSpin[2][2];
  hYieldByPtAndSpin[0][0]=new TH1F("hYieldByPtNN","Pion Yield by ptbin for B-,Y-",nptbins,pt_limits);
  hYieldByPtAndSpin[1][0]=new TH1F("hYieldByPtPN","Pion Yield by ptbin for B+,Y-",nptbins,pt_limits);
  hYieldByPtAndSpin[0][1]=new TH1F("hYieldByPtNP","Pion Yield by ptbin for B-,Y+",nptbins,pt_limits);
  hYieldByPtAndSpin[1][1]=new TH1F("hYieldByPtPP","Pion Yield by ptbin for B+,Y+",nptbins,pt_limits);

  //
  double weighted_bpol_sum[2][2];
  double weighted_ypol_sum[2][2];
  double bpolerr2_zdc2_sum[2][2]; //square of zdc counts times square of pol error
  double ypolerr2_zdc2_sum[2][2]; //square of zdc counts times square of pol error
  double bpol2_zdc_sum[2][2]; //square of polarization times square of zdc error (which is just zdc counts)
  double ypol2_zdc_sum[2][2]; //square of polarization times square of zdc error (which is just zdc counts)

  //do I have a scale separation between sum and contribution?  only 120 bins, so no.
  double zdc_narrow_sum[2][2];
  double bbc_nocut_sum[2][2];
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      zdc_narrow_sum[i][j]=0;
      bbc_nocut_sum[i][j]=0;
      weighted_bpol_sum[i][j]=0;
      weighted_ypol_sum[i][j]=0;
      bpolerr2_zdc2_sum[i][j]=0;
      ypolerr2_zdc2_sum[i][j]=0;
      bpol2_zdc_sum[i][j]=0;
      ypol2_zdc_sum[i][j]=0;
    }
  }
  
  for (int phenix_i=0;phenix_i<120;phenix_i++){
    int cad_i=(phenix_i+cross_shift)%120; //compute the C-AD bunch number
    double bpol,bpolerr,bpolsys;
    double ypol,ypolerr,ypolsys;
    spin_cont.GetPolarizationBlue(cad_i, bpol, bpolerr, bpolsys);
    spin_cont.GetPolarizationYellow(cad_i, ypol, ypolerr, ypolsys);
    int bspin = spin_cont.GetSpinPatternBlue(cad_i); //helicity of blue bunch
    int yspin = spin_cont.GetSpinPatternYellow(cad_i); //helicity of yellow bunch
    if (bspin>1 || bspin==0 || yspin > 1 || yspin==0) continue; // skip the empty and unpolarized crossings.
    int bspinbin=(bspin>0);//0 for neg. helicity, 1 for pos. helicity.
    int yspinbin=(yspin>0);
    
    
    long long  scaler_bbc_vtxcut =  spin_cont.GetScalerBbcVertexCut(cad_i);
    long long scaler_bbc_nocut  =  spin_cont.GetScalerBbcNoCut(cad_i);
    long long  scaler_zdc_wide   =  spin_cont.GetScalerZdcWide(cad_i);
    long long scaler_zdc_narrow = spin_cont.GetScalerZdcNarrow(cad_i);

    //accumulate total zdc narrow counts:
    hTotalLumi->Fill((bspinbin==yspinbin),scaler_zdc_narrow);

    //fill bare polarization plot:
    hPolarizationBySpin[bspinbin][yspinbin]->Fill(bpol,ypol);

    //accumulate average polarizations by spin configuration:
    double bpol_times_zdc=bpol*scaler_zdc_narrow;
    double ypol_times_zdc=ypol*scaler_zdc_narrow;
    weighted_bpol_sum[bspinbin][yspinbin]+=bpol_times_zdc;
    weighted_ypol_sum[bspinbin][yspinbin]+=ypol_times_zdc;

    //accumulate various moments of the weighted sums for error propagation:
    //sum of (square of zdc counts times square of pol error)
    double bpolerr_times_zdc=bpolerr*scaler_zdc_narrow;
    double ypolerr_times_zdc=ypolerr*scaler_zdc_narrow;
    bpolerr2_zdc2_sum[bspinbin][yspinbin]+=bpolerr_times_zdc*bpolerr_times_zdc;
    ypolerr2_zdc2_sum[bspinbin][yspinbin]+=ypolerr_times_zdc*ypolerr_times_zdc;

    //sum of (square of polarization times square of zdc error):
    bpol2_zdc_sum[bspinbin][yspinbin]+=bpol*bpol_times_zdc;
    ypol2_zdc_sum[bspinbin][yspinbin]+=ypol*ypol_times_zdc;

    
    //todo:  use corrected relative luminosity following pedro!
    zdc_narrow_sum[bspinbin][yspinbin]+=scaler_zdc_narrow;
    bbc_nocut_sum[bspinbin][yspinbin]+=scaler_bbc_nocut;

    //accumulate this run's data by spin pattern:
    for (int i=0;i<nptbins;i++){
      float bincenter=(pt_limits[i]+pt_limits[i]+1)/2;
      int sourcebin=hRunYieldByBunchAndPt->FindBin(bincenter,phenix_i);
      hYieldByPtAndSpin[bspinbin][yspinbin]->Fill(bincenter, hRunYieldByBunchAndPt->GetBinContent(sourcebin));
    }
  }


  //get the average polarizations for each spin state:
  double average_bpol[2][2]; //weighted average polarization for particular spin pairing
  double average_ypol[2][2]; //weightedaverage polarization for particular spin pairing
  double average_bpol_err[2][2]; //err on average polarization for that spin pairing.
  double average_ypol_err[2][2];//err on average polarization for that spin pairing.
  double bpol_sum=0;
  double ypol_sum=0;
  double zdc_sum=0;
  double bpol_error_numerator=0;
  double ypol_error_numerator=0;

  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      average_bpol[i][j]=weighted_bpol_sum[i][j]/zdc_narrow_sum[i][j];
      average_ypol[i][j]=weighted_ypol_sum[i][j]/zdc_narrow_sum[i][j];
      average_bpol_err[i][j]=sqrt(bpolerr2_zdc2_sum[i][j]
				  +bpol2_zdc_sum[i][j]
				  -average_bpol[i][j]*average_bpol[i][j]*zdc_narrow_sum[i][j]
				  )/zdc_narrow_sum[i][j];
      average_ypol_err[i][j]=sqrt(ypolerr2_zdc2_sum[i][j]
				  +ypol2_zdc_sum[i][j]
				  -average_ypol[i][j]*average_ypol[i][j]*zdc_narrow_sum[i][j]
				  )/zdc_narrow_sum[i][j];
      
      bpol_sum+=weighted_bpol_sum[i][j];
      ypol_sum+=weighted_ypol_sum[i][j];
      bpol_error_numerator+=(bpolerr2_zdc2_sum[i][j]
			     +bpol2_zdc_sum[i][j]
			     -average_bpol[i][j]*average_bpol[i][j]*zdc_narrow_sum[i][j]);
      ypol_error_numerator+=(ypolerr2_zdc2_sum[i][j]
			     +ypol2_zdc_sum[i][j]
			     -average_ypol[i][j]*average_ypol[i][j]*zdc_narrow_sum[i][j]);
      zdc_sum+=zdc_narrow_sum[i][j];

      //uncorrected lumi scalers, for like and unlike:
      hTotalLumi->Fill((i==j),zdc_narrow_sum[i][j]);      
    }
  }
  //polarization averaged over all spin states:
  double bpol_ALL=bpol_sum/zdc_sum;
  double ypol_ALL=ypol_sum/zdc_sum;
  double bpol_ALL_err=sqrt(bpol_error_numerator)/zdc_sum;
  double ypol_ALL_err=sqrt(ypol_error_numerator)/zdc_sum;



  //this should use the corrected luminosity, not the straightforward one:
  double rellumi=(zdc_narrow_sum[0][0]+zdc_narrow_sum[1][1])/(zdc_narrow_sum[0][1]+zdc_narrow_sum[1][0]);
  double rellumi_err=rellumi*sqrt(1/(zdc_narrow_sum[0][0]+zdc_narrow_sum[1][1])
				  +1/(zdc_narrow_sum[0][1]+zdc_narrow_sum[1][0]));

  TH1F* hAllByPt=new TH1F("hAllByPt","A_LL by pT;pT;A_LL",nptbins,pt_limits);
  TH1F* hDenom=new TH1F("hDenom","Denominator of ALL",nptbins,pt_limits);
  TH1F* hNumer=new TH1F("hNumer","Numerator of ALL",nptbins,pt_limits);
  
  TH1F* hRelLumi=new TH1F("hRelLumi","relative luminosity for multiplyin'",nptbins,pt_limits);

  //sum the numerator and denominator for the ALL:
  TH1F *hLikeSum=new TH1F("hLikeSum","Sum of ++ and -- bins of pion yield",nptbins,pt_limits);
  hLikeSum->Add(hYieldByPtAndSpin[0][0]);
  hLikeSum->Add(hYieldByPtAndSpin[1][1]);

  TH1F *hUnlikeSum=new TH1F("hUnlikeSum","Sum of +- and -+ bins of pion yield",nptbins,pt_limits);
  hUnlikeSum->Add(hYieldByPtAndSpin[1][0]);
  hUnlikeSum->Add(hYieldByPtAndSpin[0][1]);
  for (int i=0;i<nptbins;i++){
    double rawyield=hUnlikeSum->GetBinContent(i+1);
    double err=sqrt(rawyield*rawyield*rellumi_err*rellumi_err+rellumi*rellumi*rawyield);
    hUnlikeSum->SetBinError(i+1,err);
  }
  
  hNumer->Add(hLikeSum);
  hNumer->Add(hUnlikeSum,-rellumi);

  hDenom->Add(hLikeSum);
  hDenom->Add(hUnlikeSum,rellumi);

  
  hAllByPt->Sumw2();
  hAllByPt->Add(hNumer);
  hAllByPt->Divide(hDenom);
  hAllByPt->Scale(1/(bpol_ALL*ypol_ALL));

  //because numerator and denominator have correlated errors, sumw2 isn't sufficient to get the right uncertainty
  //so we do the math elsewhere and implement that here:
  for (int i=0;i<nptbins;i++){
    //'t' just to avoid repeating a variable I've named elsewhere in this mess.
    double tunlike=hUnlikeSum->GetBinContent(i+1);
    double trelunlike=tunlike*rellumi;
    double tlike=hLikeSum->GetBinContent(i+1);
    double tsum=tlike+trelunlike;
    double tdiff=tlike-trelunlike;
    double tasym=hAllByPt->GetBinContent(i+1);
    
    double tbpolerrterm=-tasym/bpol_ALL * bpol_ALL_err;
    double typolerrterm=-tasym/ypol_ALL * ypol_ALL_err;

    //note that this particular way of expressing it adds some 0/0 poles that root may not resolve correctly
    double tlikesumcoeff=1/(bpol_ALL*ypol_ALL)*(2*trelunlike)/(tsum*tsum);
    double tunlikesumcoeff=-1/(bpol_ALL*ypol_ALL)*(2*trelunlike*rellumi)/(tsum*tsum);
    double trellumierrterm=-1/(bpol_ALL*ypol_ALL)*(2*trelunlike*tunlike)/(tsum*tsum)*rellumi_err;


    double err2=(tbpolerrterm*tbpolerrterm
		 +typolerrterm*typolerrterm
		 +tlikesumcoeff*tlikesumcoeff*tlike
		 +tunlikesumcoeff*tunlikesumcoeff*tunlike
		 +trellumierrterm*trellumierrterm);
    
    double err=sqrt(err2);
    hAllByPt->SetBinError(i+1,err);
    //hAllByPt->SetBinContent(i+1,abs(tasym));//uncomment this if you want a wrong answer that is definitely not zero
    if (isnan(err2))
	hAllByPt->SetBinError(i+1,0.00001);
    if (tlike+trelunlike==0){
      printf("bin %d has sum=%f+%f=0\n",i,tlike,trelunlike);
	hAllByPt->SetBinContent(i+1,0);//temporary, I swear.
    }
  }
  

  TString allfname = outputdir;//"/direct/phenix+u/rosscorliss/pion_ana/output/"
  allfname += runnumber;
  allfname += ".MPC.ALL.rcc.hist.root";
  TFile *allfile = TFile::Open(allfname,"RECREATE");
  allfile->cd();
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      hYieldByPtAndSpin[i][j]->Write();
      hPolarizationBySpin[i][j]->Write();
    }
  }
  hTotalLumi->Write();
  hAllByPt->Write();

  
  hAllByPt->Draw();
  
  return;
}
