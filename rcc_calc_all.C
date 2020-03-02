#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"

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
  TFile *yieldfile = TFile::Open(yieldfname);
  TH2F *hRunYieldByBunchAndPt=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
   
  //load run scalers from the spindb, do some error checks on them
  gSystem->Load("libuspin.so");
  SpinDBContent spin_cont;
  SpinDBOutput spin_out("phnxrc");
  spin_out.StoreDBContent(runnumber, runnumber);
  spin_out.GetDBContentStore(spin_cont, runnumber);
  error = spin_cont.GetErrorValue();
  qa_level = spin_cont.GetQALevel();
  fill = spin_cont.GetFillNumber();
  badrun = spin_cont.GetBadRunFlag();
  cross_shift = spin_cont.GetCrossingShift(); //spindb stores data in the C-AD convention, which doesn't match the PHENIX numbering.  This is the relative offset.

 
  // for each bunch, accumulate yields and bunch info into the appropriate spinpattern grouping:
  TH1F *hYieldByPtAndSpin[2][2];
  hYieldByPtAndSpin[0][0]=new TH1F("hYieldByPtNN","Pion Yield by ptbin for B-,Y-",ptbins,pt_limits);
  hYieldByPtAndSpin[1][0]=new TH1F("hYieldByPtPN","Pion Yield by ptbin for B+,Y-",ptbins,pt_limits);
  hYieldByPtAndSpin[0][1]=new TH1F("hYieldByPtNP","Pion Yield by ptbin for B-,Y+",ptbins,pt_limits);
  hYieldByPtAndSpin[1][1]=new TH1F("hYieldByPtPP","Pion Yield by ptbin for B+,Y+",ptbins,pt_limits);

  //
  double weighted_bpol_sum[2][2];
  double weighted_ypol_sum[2][2];
 

  //do I have a scale separation between sum and contribution?  only 120 bins, so no.
  double zdc_narrow_sum[2][2];
  double bbc_nocut_sum[2][2];
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      zdc_narrow_sum[i][j]=0;
      bbc_nocut_sum[i][j]=0;
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
    if (bspin>1 || bspin==0 || yspin > 1 || yspin=0) continue; // skip the empty and unpolarized crossings.
    int bspinbin=(bspin>0);//0 for neg. helicity, 1 for pos. helicity.
    int yspinbin=(yspin>0);
    
    
    long long  scaler_bbc_vtxcut =  spin_cont.GetScalerBbcVertexCut(cad_i);
    long long scaler_bbc_nocut  =  spin_cont.GetScalerBbcNoCut(cad_i);
    long long  scaler_zdc_wide   =  spin_cont.GetScalerZdcWide(cad_i);
    long long scaler_zdc_narrow = spin_cont.GetScalerZdcNarrow(cad_i);

    //accumulate average polarizations by spin configuration:
    weighted_bpol_sum[bspinbin][yspinbin]+=bpol*scaler_zdc_narrow;
    weighted_ypol_sum[bspinbin][yspinbin]+=ypol*scaler_zdc_narrow;
    
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
  double average_bpol[2][2];
  double average_ypol[2][2];
  double bpol_sum=0;
  double ypol_sum=0;
  double zdc_sum=0;
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      average_bpol[i][j]=weighted_bpol_sum[i][j]/zdc_narrow_sum[i][j];
      average_ypol[i][j]=weighted_ypol_sum[i][j]/zdc_narrow_sum[i][j];
      bpol_sum+=weighted_bpol_sum[i][j];
      ypol_sum+=weighted_ypol_sum[i][j];
      zdc_sum+=zdc_narrow_sum;
    }
  }
  //averaged over all spin states:
  double bpol_ALL=bpol_sum/zdc_sum;
  double ypol_ALL=ypol_sum/zdc_sum;


  //this should use the corrected luminosity, not the straightforward one:
  double rellumi=(zdc_narrow_sum[0][0]+zdc_narrow_sum[1][1])/(zdc_narrow_sum[0][1]+zdc_narrow_sum[1][0]);

  TH1F* hAllByPt=new TH1F("hAllByPt","A_LL by pT;pT;A_LL",nptbins,pt_limits);
  TH1F* hDenom=new TH1F("hDenom","Denominator of ALL",nptbins,pt_limits);
  TH1F* hNumer=new TH1F("hNumer","Numerator of ALL",nptbins,pt_limits);

  //sum the numerator and denominator for the ALL:
  TH1F *hLikeSum= new TH1F(hYieldByPtAndSpin[0][0]);
  hLikeSum.Add(hYieldByPtAndSpin[1][1]);

  TH1F *hUnlikeSum= new TH1F(hYieldByPtAndSpin[0][1]);
  hUnlikeSum.Add(hYieldByPtAndSpin[1][0]);
	       
  hNumer.Add(hLikeSum);
  hNumer.Add(hUnlikeSum,-rellumi);

  hDenom.Add(hLikeSum);
  hDenom.Add(hUnlikeSum,rellumi);

  hAllByPt->Add(hNumer);
  hAllByPt->Divide(hDenom);
  hAllByPt->Scale(1/(bpol_ALL*ypol_ALL));

  TString allfname = outputdir;//"/direct/phenix+u/rosscorliss/pion_ana/output/"
  allfname += runnumber;
  allfname += ".MPC.ALL.rcc.hist.root";
  TFile *allfile = TFile::Open(allfname,"RECREATE");
  allfile->cd();
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      hYieldByPtAndSpin[i][i]->Write();
    }
  }
  hAllByPt->Write();

  
  hAllByPt->Draw();
  
  return;
}
