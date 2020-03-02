#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <fstream>
#include <string>
#include <TString.h>
#include <TTree.h>
#include <recoConsts.h>
#include <MpcCalib.h>
#include <MpcMap.h>
#include <TSystem.h>
#include <cmath>
#include <algorithm>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <SpinDBContent.hh>
#include <SpinDBOutput.hh>

using namespace std;

const static int MAXCLUSTERS = 100;
const static int MAXTOWERS = 1000;
const static int NBUNCHES = 120;
const static int NPTBINS = 10;
TTree *ttree;
TTree *t;
// Input tree event-wide variables
int event;
int nclus;
int ntow;
float zvtx;
short bunch;
short trig;

// Input tree cluster variables
float ecore[MAXCLUSTERS];
short feecore[MAXCLUSTERS]; // fee576ch of central tower
float pt[MAXCLUSTERS];
float x[MAXCLUSTERS];
float y[MAXCLUSTERS];
float z[MAXCLUSTERS];
short mult[MAXCLUSTERS];         // num. of towers in cluster
float disp[MAXCLUSTERS];       // max(dispx,dispy) of cluster
float chi2core[MAXCLUSTERS];   // measure of goodness of shower fit
float e8e9[MAXCLUSTERS];
short int tdc_core[MAXCLUSTERS];//for checking TDC overflow
short int lg_post_core[MAXCLUSTERS];//for checking ADC overflow
int num_adcovers;
int num_tdcovers;

// Input tree tower variables
// float etow[MAXTOWERS];
// int chtow[MAXTOWERS];

// QA objects
bool isWarn[576];
int tdcover[576];
int adcover[576];

TH1D *vtx[2][2];
TH2D *toweryields[2];
// Outfile file objects
TFile *histfile;
TFile *treefile;

TH2F *hYieldByBunchAndPt;
TH1F *ptyield[NPTBINS][2];
TH1F *checkpt;
TH1F *ptspectrum_raw[2][2];
TH1F *ptspectrum_ecut[2][2];
TH1F *ptspectrum[2][2];
TH1F *espectrum_raw[2][2];
TH1F *espectrum_ecut[2][2];
TH1F *espectrum[2][2];
TH2F *crystalcheck;
TH1F *rspectrum_pre;
TH1F *rspectrum_clustcuts;
TH1F *rspectrum_candcuts;
TH1F *chi2spec;
TH1F *e8e9spec;

TH1F *MPCA[2];
TH1F *MPCB[2];

TH1F *r_ch_066;
TH1F *r_ch_468;
TH1F *r_ch_066_post;
TH1F *r_ch_468_post;

double pT_arr[2][10][120];
double ptYields[2][10];

bool isBunchBad[NBUNCHES];
void InitWarn(int runnum);
void InitInTree();
void InitOutput(int runnum, const char *outputdir);
void InitDB(int n_runnum);
void get_entry(int ientry);
bool PassesClusterCuts(int iclus);
void End();

SpinDBContent spin_cont;
MpcMap *mpcmap;
recoConsts *rc;
void rcc_gen_yield(int runnum,
		   const char * inputdir="/",
		   const char * outputdir="/phenix/spin/spin1/phnxsp01/rosscorliss/trees/") {
  //MB is at:
  //        "/phenix/spin/phnxsp01/pmontu/taxi/Run13pp510MinBias/10699/data/";
  TString fullfile = inputdir;//"/phenix/spin/phnxsp01/rosscorliss/taxi/Run13pp510MPC/15944/data/";
  fullfile += runnum;
  fullfile += ".root";

  TFile *rootin = new TFile(fullfile, "READONLY");
  if ((rootin->IsZombie()) || !rootin->GetNkeys())
    return;

  ttree = (TTree *)rootin->Get("T");

  InitInTree();
  InitOutput(runnum, outputdir);
  InitDB(runnum);
  InitWarn(runnum);

  int is_north, even_or_odd, spin_pattern;
  double cluster_r;
  int nentries = ttree->GetEntries();
  cout << "Number of entries: " << nentries << endl;
  for (int ievent = 0; ievent < nentries; ievent++) {
    if (ievent % 10000 == 0)
      cout << "event: " << ievent << endl;
    get_entry(ievent);

    int corrbunch = (bunch + spin_cont.GetCrossingShift()) %
                    120;
    
    even_or_odd = (corrbunch % 2); // 0 for even, 1 for odd

    for (int iclus = 0; iclus < nclus; iclus++) {
      is_north = (feecore[iclus] < 288) ? 0 : 1;

      cluster_r=sqrt(x[iclus] * x[iclus] + y[iclus] * y[iclus]);

      crystalcheck->Fill(feecore[iclus],mpcmap->isCrystal(feecore[iclus]));
      rspectrum_pre->Fill(cluster_r);
      chi2spec->Fill(chi2core[iclus]);
      e8e9spec->Fill(e8e9[iclus]);

      if (feecore[iclus] == 66)
	r_ch_066->Fill(cluster_r);
      if (feecore[iclus] == 468)
	r_ch_468->Fill(cluster_r);
      //rspectrum->Fill(cluster_r);
      
      if (trig & (1 << 5))
	MPCA[is_north]->Fill(pt[iclus]);
      if (trig & (1 << 4))
	MPCB[is_north]->Fill(pt[iclus]);

      if (!PassesClusterCuts(iclus))
        continue;
      rspectrum_clustcuts->Fill(cluster_r);

      ptspectrum_raw[is_north][corrbunch % 2]->Fill(pt[iclus]);
      espectrum_raw[is_north][corrbunch % 2]->Fill(ecore[iclus]);
      
      if (ecore[iclus] < 15)
	continue; // looking for merged clusters

      ptspectrum_ecut[is_north][corrbunch % 2]->Fill(pt[iclus]);
      espectrum_ecut[is_north][corrbunch % 2]->Fill(ecore[iclus]);

      if (mult[iclus] <= 2)
	continue;
      if (disp[iclus] < 0.0005)
      	continue;
      if (chi2core[iclus] > 30.)
      	continue;
      if (e8e9[iclus] < 0.2)
	continue;
      
      if (feecore[iclus] == 66)
	r_ch_066_post->Fill(cluster_r);
      if (feecore[iclus] == 468)
	r_ch_468_post->Fill(cluster_r);
      rspectrum_candcuts->Fill(cluster_r);


      vtx[is_north][corrbunch % 2]->Fill(zvtx);

      int ptbin = checkpt->Fill(pt[iclus]);
 
       if (ptbin<0){
	 std::cout << "under/overflow.  not including event." << std::endl;
	 continue;
       }
       
      pT_arr[is_north][ptbin-1][corrbunch]++;
      ptyield[ptbin - 1][is_north]->Fill(pt[iclus]);
      ptspectrum[is_north][corrbunch % 2]->Fill(pt[iclus]);
      espectrum[is_north][corrbunch % 2]->Fill(ecore[iclus]);

      int ix = mpcmap->getGridX(feecore[iclus]);
      int iy = mpcmap->getGridY(feecore[iclus]);
      toweryields[is_north]->Fill(ix, iy);
      hYieldByBunchAndPt->Fill(pt[iclus],bunch);
    }
  }
  rootin->Close();
  delete rootin;
  End();
}

void get_entry(int ientry) {
  ttree->GetEntry(ientry);
  return;
}

void InitInTree() {
  // Set Branches to point to variables in this macro
  ttree->SetBranchAddress("event", &event);
  ttree->SetBranchAddress("nclus", &nclus);
  ttree->SetBranchAddress("bbcvtx", &zvtx);
  ttree->SetBranchAddress("bunch", &bunch);
  ttree->SetBranchAddress("trig", &trig);
  ttree->SetBranchAddress("ecore", ecore);
  ttree->SetBranchAddress("feecore", feecore);
  ttree->SetBranchAddress("pt", pt);
  ttree->SetBranchAddress("x", x);
  ttree->SetBranchAddress("y", y);
  ttree->SetBranchAddress("mult", mult);
  ttree->SetBranchAddress("disp", disp);
  ttree->SetBranchAddress("chi2core", chi2core);
  ttree->SetBranchAddress("e8e9", e8e9);
  // ttree->SetBranchAddress("tdc_max", tdc_core);
  // ttree->SetBranchAddress("lg_post_core", lg_post_core);
  return;
}

void InitOutput(int runnum, const char* outputdir){
  TString yieldfname = outputdir;//"/direct/phenix+u/rosscorliss/pion_ana/output/"
  yieldfname += runnum;
  yieldfname += ".MPC.yields.rcc.hist.root";
  histfile = new TFile(yieldfname, "RECREATE");
  printf("opening SpinDb histfile at %s\n",yieldfname.Data());
  Int_t sparsebins[4] = {2, 2, 120, 10}; // N/S, Even/Odd,crossing num., NPTBINS
  // spin patterns are in order: ++,+-,--,-+
  Double_t sparsebinsmin[4] = {-.5, -.5, -.5, 1.0};
  Double_t sparsebinsmax[4] = {1.5, 1.5, 119.5, 12};
  Double_t pt_limits[11] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};

  hYieldByBunchAndPt=new TH2F("hYieldByBunchAndPt","yield by bunch and pt",10,pt_limits,120,-0.5,119.5);




  
  // TString sparsename = "clus_yields";
  // outsparse = new THnSparseD(sparsename, sparsename, 4, sparsebins,
  //                            sparsebinsmin, sparsebinsmax);
  // outsparse->GetAxis(3)->Set(9, pt_limits);
  vtx[0][0] = new TH1D("evenvtxS", "evenvtxS", 600, -300, 300);
  vtx[0][1] = new TH1D("oddvtxS", "oddvtxS", 600, -300, 300);
  vtx[1][0] = new TH1D("evenvtxN", "evenvtxN", 600, -300, 300);
  vtx[1][1] = new TH1D("oddvtxN", "oddvtxN", 600, -300, 300);
  toweryields[0] =
      new TH2D("toweryieldsS", "toweryieldsS", 18, -.5, 17.5, 18, -.5, 17.5);
  toweryields[1] =
      new TH2D("toweryieldsN", "toweryieldsN", 18, -.5, 17.5, 18, -.5, 17.5);
  TString yieldname;
  for (int ipt = 0; ipt < NPTBINS; ipt++) {
    for (int iarm = 0; iarm <= 1; iarm++) {
      yieldname = "hpt_ptbin_";
      yieldname += ipt;
      yieldname += "_arm_";
      yieldname += iarm;
      ptyield[ipt][iarm] = new TH1F(yieldname, yieldname, 100, pt_limits[ipt],
                                    pt_limits[ipt + 1]);
    }
  }

  for (std::size_t iarm = 0; iarm < 2; iarm++) {
    for (std::size_t ipt = 0; ipt < 10; ipt++) {
      for (std::size_t ibunch = 0; ibunch < 120; ibunch++) {
	pT_arr[iarm][ipt][ibunch] = 0;
      }
    }
  }
  
  checkpt = new TH1F("checkpt", "checkpt", 10, 0, 10);
  checkpt->GetXaxis()->Set(10, pt_limits);

  ptspectrum_raw[0][0] =
      new TH1F("ptspectrum_rawarm0even", "ptspectrum_rawarm0even", 300, 0, 60);
  ptspectrum_raw[0][1] =
      new TH1F("ptspectrum_rawarm0odd", "ptspectrum_rawarm0odd", 300, 0, 60);
  ptspectrum_raw[1][0] =
      new TH1F("ptspectrum_rawarm1even", "ptspectrum_rawarm1even", 300, 0, 60);
  ptspectrum_raw[1][1] =
      new TH1F("ptspectrum_rawarm1odd", "ptspectrum_rawarm1odd", 300, 0, 60);

  ptspectrum_ecut[0][0] =
      new TH1F("ptspectrum_ecutarm0even", "ptspectrum_ecutarm0even", 300, 0, 60);
  ptspectrum_ecut[0][1] =
      new TH1F("ptspectrum_ecutarm0odd", "ptspectrum_ecutarm0odd", 300, 0, 60);
  ptspectrum_ecut[1][0] =
      new TH1F("ptspectrum_ecutarm1even", "ptspectrum_ecutarm1even", 300, 0, 60);
  ptspectrum_ecut[1][1] =
      new TH1F("ptspectrum_ecutarm1odd", "ptspectrum_ecutarm1odd", 300, 0, 60);

  ptspectrum[0][0] =
      new TH1F("ptspectrumarm0even", "ptspectrumarm0even", 300, 0, 60);
  ptspectrum[0][1] =
      new TH1F("ptspectrumarm0odd", "ptspectrumarm0odd", 300, 0, 60);
  ptspectrum[1][0] =
      new TH1F("ptspectrumarm1even", "ptspectrumarm1even", 300, 0, 60);
  ptspectrum[1][1] =
      new TH1F("ptspectrumarm1odd", "ptspectrumarm1odd", 300, 0, 60);

  MPCA[0] = new TH1F("MPC_A_S", "MPC_A only trigger, South Arm", 300, 0, 60);
  MPCA[1] = new TH1F("MPC_A_N", "MPC_A only trigger, North Arm", 300, 0, 60);

  MPCB[0] = new TH1F("MPC_B_S", "MPC_B trigger, South Arm", 300, 0, 60);
  MPCB[1] = new TH1F("MPC_B_N", "MPC_B trigger, North Arm", 300, 0, 60);
  
  espectrum_raw[0][0] =
      new TH1F("espectrum_rawarm0even", "espectrum_rawarm0even", 300, 0, 600);
  espectrum_raw[0][1] =
      new TH1F("espectrum_rawarm0odd", "espectrum_rawarm0odd", 300, 0, 600);
  espectrum_raw[1][0] =
      new TH1F("espectrum_rawarm1even", "espectrum_rawarm1even", 300, 0, 600);
  espectrum_raw[1][1] =
      new TH1F("espectrum_rawarm1odd", "espectrum_rawarm1odd", 300, 0, 600);

  espectrum_ecut[0][0] =
      new TH1F("espectrum_ecutarm0even", "espectrum_ecutarm0even", 300, 0, 600);
  espectrum_ecut[0][1] =
      new TH1F("espectrum_ecutarm0odd", "espectrum_ecutarm0odd", 300, 0, 600);
  espectrum_ecut[1][0] =
      new TH1F("espectrum_ecutarm1even", "espectrum_ecutarm1even", 300, 0, 600);
  espectrum_ecut[1][1] =
      new TH1F("espectrum_ecutarm1odd", "espectrum_ecutarm1odd", 300, 0, 600);

  espectrum[0][0] =
      new TH1F("espectrumarm0even", "espectrumarm0even", 300, 0, 600);
  espectrum[0][1] =
      new TH1F("espectrumarm0odd", "espectrumarm0odd", 300, 0, 600);
  espectrum[1][0] =
      new TH1F("espectrumarm1even", "espectrumarm1even", 300, 0, 600);
  espectrum[1][1] =
      new TH1F("espectrumarm1odd", "espectrumarm1odd", 300, 0, 600);

  crystalcheck=new TH2F("crystalcheck","isCrystal() flag vs feecore",500,-0.5,499.5,2,-0.5,1.5);

  rspectrum_clustcuts = new TH1F("rspectrum_clustcuts", "r distribution after clust cuts", 100, 10, 20);
  rspectrum_candcuts = new TH1F("rspectrum_candcuts", "r distribution after clust and candidate cuts", 100, 10, 20);
  rspectrum_pre = new TH1F("rspectrum_pre", "r distribution before cuts", 100, 10, 20);

  chi2spec = new TH1F("chi2spec", "Chi2 distribution - before cuts", 100, 0, 100);
  e8e9spec = new TH1F("e8e9spec", "E8/E9 Ratio", 100, 0, 1);

  r_ch_066 = new TH1F("r_ch_066", "r distribution for channel 066", 200, 10, 20);
  r_ch_468 = new TH1F("r_ch_468", "r distribution for channel 468", 200, 10, 20);
  r_ch_066_post = new TH1F("r_ch_066_post", "r distribution for channel 066", 200, 10, 20);
  r_ch_468_post = new TH1F("r_ch_468_post", "r distribution for channel 468", 200, 10, 20);

  TString ttreefname = "/direct/phenix+u/rosscorliss/pion_ana/output/"
    "relative_luminosity.SpinDB.clusters.";
  if (minbias)
    ttreefname += "MinBias.trees.";
  else
    ttreefname += "MPC.trees.";
  ttreefname += runnum;
  ttreefname += ".root";
  treefile = new TFile(ttreefname, "RECREATE");
  
  t = new TTree("t", "t");

  t->Branch("ptbin0arm0", &ptYields[0][0], "ptbin0arm0/D");
  t->Branch("ptbin0arm1", &ptYields[1][0], "ptbin0arm1/D");
  t->Branch("ptbin1arm0", &ptYields[0][1], "ptbin1arm0/D");
  t->Branch("ptbin1arm1", &ptYields[1][1], "ptbin1arm1/D");
  t->Branch("ptbin2arm0", &ptYields[0][2], "ptbin2arm0/D");
  t->Branch("ptbin2arm1", &ptYields[1][2], "ptbin2arm1/D");
  t->Branch("ptbin3arm0", &ptYields[0][3], "ptbin3arm0/D");
  t->Branch("ptbin3arm1", &ptYields[1][3], "ptbin3arm1/D");
  t->Branch("ptbin4arm0", &ptYields[0][4], "ptbin4arm0/D");
  t->Branch("ptbin4arm1", &ptYields[1][4], "ptbin4arm1/D");
  t->Branch("ptbin5arm0", &ptYields[0][5], "ptbin5arm0/D");
  t->Branch("ptbin5arm1", &ptYields[1][5], "ptbin5arm1/D");
  t->Branch("ptbin6arm0", &ptYields[0][6], "ptbin6arm0/D");
  t->Branch("ptbin6arm1", &ptYields[1][6], "ptbin6arm1/D");
  t->Branch("ptbin7arm0", &ptYields[0][7], "ptbin7arm0/D");
  t->Branch("ptbin7arm1", &ptYields[1][7], "ptbin7arm1/D");
  t->Branch("ptbin8arm0", &ptYields[0][8], "ptbin8arm0/D");
  t->Branch("ptbin8arm1", &ptYields[1][8], "ptbin8arm1/D");
  t->Branch("ptbin9arm0", &ptYields[0][9], "ptbin9arm0/D");
  t->Branch("ptbin9arm1", &ptYields[1][9], "ptbin9arm1/D");

  return;
}

void InitDB(int n_runnum) {
  SpinDBOutput spin_out("phnxrc");
  spin_out.StoreDBContent(n_runnum, n_runnum);
  if (spin_out.CheckRunRowStore(n_runnum) != 1) {
    cout << "Error connecting to DB" << endl;
    return;
  }
  spin_out.GetDBContentStore(spin_cont, n_runnum);
  int fillnum = spin_cont.GetFillNumber();

  rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", n_runnum);
  mpcmap = MpcMap::instance();

  // Populate array that tells which bunches to exclude because they fail QA
  // checks
  char sname[200];

  sprintf(sname, "/phenix/spin2/pmontu/offline/analysis/pmontu/"
                 "relative_luminosity/SpinDB/good_bunches/%d.txt",
          n_runnum);
  printf("loading good bunches list from %s\b",sname);
  ifstream bstatusfile(sname);
  if (!bstatusfile.good())
    cout << "Error! Can't open bunch status file." << endl;
  int ibunch = 0;
  int ibunchstatus = 0;
  while (bstatusfile.good()) {
    bstatusfile >> ibunch >> ibunchstatus;
    if (!ibunchstatus)
      isBunchBad[ibunch] = true; // status in file is 1 for good, 0 for bad
  }
  printf("good bunches loaded.\n");
  return;
}

bool PassesClusterCuts(int iclus) {
  if (!mpcmap->isCrystal(feecore[iclus]))
    return false;
  if (isWarn[feecore[iclus]])
    return false;
  float r = sqrt(x[iclus] * x[iclus] + y[iclus] * y[iclus]);
  if ((r < 11) || (r > 19))
    return false;
  if (pt[iclus] < 1.)
    return false; // lower than low edge of lowest bin
  // Clusterness cuts

  // Overflow cuts
  if (disp[iclus] > 4)
    return false; // no chi2 cut since we're looking for overlapping clusters

  // if (tdc_core[iclus] > tdcover[feecore[iclus]]) {
  //   num_tdcovers++;
  //   return false;
  // }
  // if (lg_post_core[iclus] > adcover[feecore[iclus]]) {
  //   num_adcovers++;
  //   return false;
  // }

  return true;
}

void InitWarn(int runnum) {
  for (int ich = 0; ich < 576; ich++) {
    isWarn[ich] = false;
  }
  printf("seeking warnmap\n");
  TString warnfile = "/phenix/spin2/pmontu/offline/packages/mpc/calibrations/"
                     "pi0cal_fast/macros/iterative/db_update/final_warnmap.txt";
  printf("loading warnmap from %s\n",warnfile.Data());
  ifstream warnmap(warnfile.Data());
  int wch;
  while (warnmap.good()) {
    warnmap >> wch;
    if (wch == -1)
      break;
    isWarn[wch] = true;
  }
  warnmap.close();
  return;
}

void InitOverflows() {
  // ifstream
  // tdcoverfile("/direct/phenix+u/cmckinn4/run11/analysis/Calibrations/TdcOverflow_final_m15.txt");
  TString tdcfilename="/direct/phenix+u/cmckinn4/run11/analysis/Calibrations/"
                       "MpcCal_Run11.overflow";
  printf("loading overflow mpc calibration from %s\n",tdcfilename.Data());
  ifstream tdcoverfile(tdcfilename.Data());
  int och;
  float oflow, error;
  while (tdcoverfile.good()) {
    tdcoverfile >> och >> oflow >> error;
    if ((och >= 0) && (och < 576)) {
      tdcover[och] = oflow;
    }
  }
  tdcoverfile.close();
  TString adcfilename="/direct/phenix+u/cmckinn4/run11/analysis/Calibrations/"
                       "MpcCal_Run11_lopostoverflow.txt";
  printf("loading adcoverfile from %s\n",adcfilename.Data());
  ifstream adcoverfile(adcfilename.Data());
  while (adcoverfile.good()) {
    adcoverfile >> och >> oflow;
    if ((och >= 0) && (och < 576))
      adcover[och] = oflow;
  }
  adcoverfile.close();

  return;
}

void End() {
  cout << "Wrapping up." << endl;
  // ofstream rawovers("adc_tdc_overs.txt", ios::app);
  // rawovers << num_adcovers << "\t" << num_tdcovers << endl;
  // rawovers.close();
  ofstream rawovers("adc_tdc_overs.txt",ios::app);
  rawovers << num_adcovers << "\t" << num_tdcovers << endl;
  rawovers.close();
  histfile->cd();

  hYieldByBunchAndPt->Write();




  // outsparse->Write();
  // delete outsparse;
  vtx[0][0]->Write();
  vtx[0][1]->Write();
  vtx[1][0]->Write();
  vtx[1][1]->Write();
  toweryields[0]->Write();
  toweryields[1]->Write();
  for (int ipt = 0; ipt < NPTBINS; ipt++) {
    for (int iarm = 0; iarm <= 1; iarm++) {
      ptyield[ipt][iarm]->Write();
    }
  }


  

  ptspectrum_raw[0][0]->Write();
  ptspectrum_raw[0][1]->Write();
  ptspectrum_raw[1][0]->Write();
  ptspectrum_raw[1][1]->Write();

  ptspectrum_ecut[0][0]->Write();
  ptspectrum_ecut[0][1]->Write();
  ptspectrum_ecut[1][0]->Write();
  ptspectrum_ecut[1][1]->Write();

  ptspectrum[0][0]->Write();
  ptspectrum[0][1]->Write();
  ptspectrum[1][0]->Write();
  ptspectrum[1][1]->Write();

  espectrum_raw[0][0]->Write();
  espectrum_raw[0][1]->Write();
  espectrum_raw[1][0]->Write();
  espectrum_raw[1][1]->Write();

  espectrum_ecut[0][0]->Write();
  espectrum_ecut[0][1]->Write();
  espectrum_ecut[1][0]->Write();
  espectrum_ecut[1][1]->Write();

  espectrum[0][0]->Write();
  espectrum[0][1]->Write();
  espectrum[1][0]->Write();
  espectrum[1][1]->Write();

  crystalcheck->Write();

  rspectrum_clustcuts->Write();
  rspectrum_candcuts->Write();
  rspectrum_pre->Write();
  chi2spec->Write();
  e8e9spec->Write();

  MPCA[0]->SetMarkerStyle(20);
  MPCA[0]->SetMarkerColor(kBlue);
  MPCA[0]->SetMarkerSize(1.3);
  MPCA[1]->SetMarkerStyle(22);
  MPCA[1]->SetMarkerColor(kBlue);
  MPCA[1]->SetMarkerSize(1.3);
  MPCB[0]->SetMarkerStyle(20);
  MPCB[0]->SetMarkerColor(kRed);
  MPCB[0]->SetMarkerSize(1.3);
  MPCB[1]->SetMarkerStyle(22);
  MPCB[1]->SetMarkerColor(kRed);
  MPCB[1]->SetMarkerSize(1.3);

  MPCA[0]->Write();
  MPCA[1]->Write();
  MPCB[0]->Write();
  MPCB[1]->Write();

  r_ch_066->Write();
  r_ch_468->Write();
  r_ch_066_post->Write();
  r_ch_468_post->Write();

  histfile->Close();
  delete histfile;

  for (std::size_t ibunch = 0; ibunch < 120; ibunch++) {
    for (int ipt = 0; ipt < NPTBINS; ipt++) {
      for (int iarm = 0; iarm <= 1; iarm++) {
	ptYields[iarm][ipt] = pT_arr[iarm][ipt][ibunch];
      }
    }
    t->Fill();
  }
  treefile->cd();
  t->Write();
  t->ResetBranchAddresses();
  treefile->Close();
  delete treefile;
  
  return;
}
