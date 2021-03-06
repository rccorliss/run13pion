#include <TString.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <SpinDBContent.hh>
#include <SpinDBOutput.hh>
#include <assert.h>

int lookupSpinPattern(char bpat, char ypat);


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

  //this assumes yields vs pt_limit bins have been put into the following histogram:
  TH2F *hRunYieldByBunchAndPt=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
  //the rest of the data needed to generate asymmetries are loaded from the spin db.

  
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


  //accumulate a few numbers in a hTags holder:
  TH1F *hTags=new TH1F("hTags","Various named numbers that can be stored safely as floats",100,0,100);
  //firstly, the run number, just in case:
  hTags->Fill("runnumber",runnumber);
  hTags->Fill("fill",fill);

  //accumulate the luminosity in a 1x2 histogram:
  TH1F *hTotalLumi=new TH1F("hTotalLumi","Unlike (like) zdc-narrow lumi scaler sums;unlike=0,like=1;sum",2,-0.5,1.5);
  TH1F *hTotalZdcWide=new TH1F("hTotalZdcWide","Unlike (like) zdc-wide lumi scaler sums;unlike=0,like=1;sum",2,-0.5,1.5);
  TH1F *hTotalBbc=new TH1F("hTotalBbc","Unlike (like) bbc-narrow lumi scaler sums;unlike=0,like=1;sum",2,-0.5,1.5);
  TH1F *hTotalBbcWide=new TH1F("hTotalBbcWide","Unlike (like) bbc-wide lumi scaler sums;unlike=0,like=1;sum",2,-0.5,1.5);

  //polarizations are actually single-valued across a run, but these accumulate them bunch by bunch
  //which allows that to be cross-checked for consistency:
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

  //polarization variables
  //NOTE:  detailed error handling is not needed.  See comment in for loop.
  //these were previously per-bunch, but now just overall values:
  double bpol=0, bpolerr=0, bpolsys=0;
  double ypol=0, ypolerr=0, ypolsys=0;
  
  //double weighted_bpol_sum[2][2];
  //double weighted_ypol_sum[2][2];
  //double bpolerr2_zdc2_sum[2][2]; //square of zdc counts times square of pol error
  //double ypolerr2_zdc2_sum[2][2]; //square of zdc counts times square of pol error
  //double bpol2_zdc_sum[2][2]; //square of polarization times square of zdc error (which is just zdc counts)
  //double ypol2_zdc_sum[2][2]; //square of polarization times square of zdc error (which is just zdc counts)

  //do I have a scale separation between sum and contribution?  only 120 bins, so no.
  double zdc_narrow_sum[2][2];
  double zdc_wide_sum[2][2];
  double bbc_vtxcut_sum[2][2];
  double bbc_nocut_sum[2][2];
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      zdc_narrow_sum[i][j]=0;
      zdc_wide_sum[i][j]=0;
      bbc_vtxcut_sum[i][j]=0;
      bbc_nocut_sum[i][j]=0;
      //weighted_bpol_sum[i][j]=0;
      //weighted_ypol_sum[i][j]=0;
      //bpolerr2_zdc2_sum[i][j]=0;
      //ypolerr2_zdc2_sum[i][j]=0;
      //bpol2_zdc_sum[i][j]=0;
      //ypol2_zdc_sum[i][j]=0;
    }
  }

  //get the spin pattern:
  char bspinpat=0;
  char yspinpat=0;
  for (int cad_i=0;cad_i<8;cad_i++){
    int bspin = spin_cont.GetSpinPatternBlue(cad_i); //helicity of blue bunch
    int yspin = spin_cont.GetSpinPatternYellow(cad_i); //helicity of yellow bunch
    if (bspin>1 || bspin==0 || yspin > 1 || yspin==0) {
      printf("an empty bunch in the first eight in cad numbering!\n");
      assert (1==2);//panic.  this shouldn't happen.
    }
    char bbit=(bspin>0);//0 for neg. helicity, 1 for pos. helicity.
    char ybit=(yspin>0);
    bspinpat|=(bbit<<(7-cad_i));
    yspinpat|=(ybit<<(7-cad_i));
    if(0){
      //checks to make sure I'm loading the spin pattern correctly.
      printf("blue spinpat:%d%d%d%d%d%d%d%d\n",
	     (bspinpat&128)>0,(bspinpat&64)>0,(bspinpat&32)>0,(bspinpat&16)>0,
	     (bspinpat&8)>0,(bspinpat&4)>0,(bspinpat&2)>0,(bspinpat&1)>0);
      printf("yellow spinpat:%d%d%d%d%d%d%d%d\n",
	     (yspinpat&128)>0,(yspinpat&64)>0,(yspinpat&32)>0,(yspinpat&16)>0,
	     (yspinpat&8)>0,(yspinpat&4)>0,(yspinpat&2)>0,(yspinpat&1)>0);
    }
  }
  int spinpattern=lookupSpinPattern(bspinpat,  yspinpat);
  hTags->Fill("spinpattern",spinpattern);
 

    
  for (int phenix_i=0;phenix_i<120;phenix_i++){
    int cad_i=(phenix_i+cross_shift)%120; //compute the C-AD bunch number
    //double bpol,bpolerr,bpolsys;
    //double ypol,ypolerr,ypolsys;
    if (phenix_i==0){
    spin_cont.GetPolarizationBlue(cad_i, bpol, bpolerr, bpolsys);
    spin_cont.GetPolarizationYellow(cad_i, ypol, ypolerr, ypolsys);

    hTags->Fill("ypolerr",ypolerr);
    hTags->Fill("bpolerr",bpolerr);
   hTags->Fill("ypol",ypol);
    hTags->Fill("bpol",bpol);
 
    
    }
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
    hTotalZdcWide->Fill((bspinbin==yspinbin),scaler_zdc_wide);
    hTotalBbcWide->Fill((bspinbin==yspinbin),scaler_bbc_nocut);
    hTotalBbc->Fill((bspinbin==yspinbin),scaler_bbc_vtxcut);


    //manage polarization monitoring and averaging:
    //NOTE:  the polarization is ~always constant across all bunches, so it is inappropriate to treat the polarization as the weighted average of all the measurements -- they're completely correlated, and this only serves to artificially reduce the polarization error.
    
    //fill bare polarization plot:
    hPolarizationBySpin[bspinbin][yspinbin]->Fill(bpol,ypol);
      
    //accumulate average polarizations by spin configuration:
    //double bpol_times_zdc=bpol*scaler_zdc_narrow;
    //double ypol_times_zdc=ypol*scaler_zdc_narrow;
    //weighted_bpol_sum[bspinbin][yspinbin]+=bpol_times_zdc;
    //weighted_ypol_sum[bspinbin][yspinbin]+=ypol_times_zdc;

    //accumulate various moments of the weighted sums for error propagation:
    //sum of (square of zdc counts times square of pol error)
    //double bpolerr_times_zdc=bpolerr*scaler_zdc_narrow;
    //double ypolerr_times_zdc=ypolerr*scaler_zdc_narrow;
    //bpolerr2_zdc2_sum[bspinbin][yspinbin]+=bpolerr_times_zdc*bpolerr_times_zdc;
    //ypolerr2_zdc2_sum[bspinbin][yspinbin]+=ypolerr_times_zdc*ypolerr_times_zdc;

    //sum of (square of polarization times square of zdc error):
    //bpol2_zdc_sum[bspinbin][yspinbin]+=bpol*bpol_times_zdc;
    //ypol2_zdc_sum[bspinbin][yspinbin]+=ypol*ypol_times_zdc;


    
    
    //todo:  use corrected relative luminosity following pedro!
    zdc_narrow_sum[bspinbin][yspinbin]+=scaler_zdc_narrow;
    bbc_nocut_sum[bspinbin][yspinbin]+=scaler_bbc_nocut;
    zdc_wide_sum[bspinbin][yspinbin]+=scaler_zdc_wide;
    bbc_vtxcut_sum[bspinbin][yspinbin]+=scaler_bbc_vtxcut;

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
      /*
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
      */
      zdc_sum+=zdc_narrow_sum[i][j];

      //uncorrected lumi scalers, for like and unlike:
      //this is done in the loop above. hTotalLumi->Fill((i==j),zdc_narrow_sum[i][j]);      
    }
  }
  //polarization averaged over all spin states:
  //per the above, this is single-valued, so we just copy the data from phenix-numbered bunch 0:
  double bpol_ALL=bpol;//bpol_sum/zdc_sum;
  double ypol_ALL=ypol;//ypol_sum/zdc_sum;
  double bpol_ALL_err=bpolerr;//sqrt(bpol_error_numerator)/zdc_sum;
  double ypol_ALL_err=ypolerr;//sqrt(ypol_error_numerator)/zdc_sum;



  //this should use the corrected luminosity, not the straightforward one:
  //options are zdc_narrow_sum , zdc_wide_sum, bbc_vtxcut_sum, bbc_nocut_sum
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

  
  hAllByPt->Add(hNumer);
  hAllByPt->Divide(hDenom);
  hAllByPt->Scale(1/(bpol_ALL*ypol_ALL));

  //because numerator and denominator have correlated errors, sumw2 isn't sufficient to get the right uncertainty
  //so we do the math elsewhere and implement that here:
  for (int i=0;i<nptbins;i++){
    //prepend 't' to avoid repeating a variable I've named elsewhere in this mess.
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
	hAllByPt->SetBinError(i+1,1.0);
    if (tsum==0){
      printf("bin %d has sum=%f+%f=0\n",i,tlike,trelunlike);
	hAllByPt->SetBinContent(i+1,0);//temporary, I swear.
    }
  }
  

  TString allfname = outputdir;//"/direct/phenix+u/rosscorliss/pion_ana/output/"
  allfname += runnumber;
  allfname += ".MPC.ALL.rcc.hist.root";
  TFile *allfile = TFile::Open(allfname,"RECREATE");
  allfile->cd();

  hTags->Write();
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      hYieldByPtAndSpin[i][j]->Write();
      hPolarizationBySpin[i][j]->Write();
    }
  }
  hTotalLumi->Write();
  hTotalZdcWide->Write();
  hTotalBbc->Write();
  hTotalBbcWide->Write();
  hAllByPt->Write();

  
  hAllByPt->Draw();
  
  return;
}






int lookupSpinPattern(char bpat, char ypat){
  //take eight consecutive spinbits, where the lowest bit is the 7th bunch, and the highest bit is the 0th
  //look these up against the labeled spin patterns and return the correct spin pattern number.
  //assumes '1' = spin up, and '0' = spin down;

  //note that all valid patterns have pairs of bunches in the same config.
  //if we check this parity first, we only have to check four bits instead of 8
  char bshort=0;
  char yshort=0;
  for (int bit=0;bit<8;bit+=2){
    char bbit=(bpat>>bit)&1;//just the bit'th bit, shifted to the ones place.
    char bnextbit=(bpat>>(bit+1))&1;
    //printf("comparing bbits: %d%d\n", bbit,bnextbit);
    if ( bbit!=bnextbit) {
      printf("parity wrong in spin pattern\n");
      return -2;
    }
    bshort+=bbit<<(bit/2);//make the 'even bits only' object.
    
    char ybit=(ypat>>bit)&1;//just the bit'th bit, shifted to the ones place.
     char ynextbit=(ypat>>(bit+1))&1;
   if ( ybit!=ynextbit) {
      printf("parity wrong in spin pattern\n");
      return -2;
    }
    yshort+=ybit<<(bit/2);//make the 'even bits only' object.
  }
  
  char base[6];
  base[0]=0b1010;//base1 and 1a
  base[1]=0b0101;//base2 and 2a
  base[2]=0b1100;//base3  
  base[3]=0b0011;//base4
  base[4]=0b0110;//base3a
  base[5]=0b1001;//base4a

  //define the valid patterns:
  int patlookup[6][6];//blue then yellow
  for (int i=0;i<6;i++){
    for (int j=0;j<6;j++){
      patlookup[i][j]=0;
    }
  }
  //defined by AN1125:
  patlookup[0][2]=1;//B=1, Y=3
  patlookup[1][2]=2;//B=2, Y=4
  patlookup[0][3]=3;//add one both to indices, etc
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
    if (bshort==base[i]) {
      blabel=i;
      break;
    }
  }
  for (int i=0;i<6;i++){
    if (yshort==base[i]) {
      ylabel=i;
      break;
    }
  }
  if (blabel<0 || ylabel<0){
    printf("couldn't find a spin pattern.  parity is good, but not balanced.\n");
    return -1;
  }
  return patlookup[blabel][ylabel];
}
  
