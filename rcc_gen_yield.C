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
const static float zpos=220;//cm
TTree *ttree;
TTree *t;
// Input tree event-wide variables
int event;
int nclus;
int ntow;
float zvtx;
short bunch;
short trig;

// Input tree cluster variables used for cluster cuts:
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

TTree *rccRunTree;
//contains: run,fill,neve (not clusters!),nbins,nbounds(=nbins+1)
int rccRun, rccFill, rccNeve, rccTotRawClust, rccTotGoodClust, rccNbins, rccNbounds;
const Float_t rccBounds[]={1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};
TTree *rccBunchTree;
TTree *rccClusterTree;
TTree *splitClusterTree;
int rccBunch, rccIx, rccIy, rccFeecore,rccMult;
float rccX, rccY, rccVtx, rccEcore,rccE8, rccE9, rccDisp,rccChi;
bool rccNorth;

//contains nEvenBunchEve,nOddBunchEve,nBunchEve
//(iDet=0 ==> north, 1==> south, 2==>both
//to make sure I have the implicit ordering correctly, doing this out long-hand:
int rccRawClust[NBUNCHES*NPTBINS];//raw events in this bunch and ptbin, index by iBun*NPT+iPt)
int rccGoodClust[NBUNCHES*NPTBINS];//good clusters in this bunch and ptbin
int rccRawClustN[NBUNCHES*NPTBINS];//raw events in this bunch and ptbin, index by iBun*NPT+iPt)
int rccGoodClustN[NBUNCHES*NPTBINS];//good clusters in this bunch and ptbin
int rccRawClustS[NBUNCHES*NPTBINS];//raw events in this bunch and ptbin, index by iBun*NPT+iPt)
int rccGoodClustS[NBUNCHES*NPTBINS];//good clusters in this bunch and ptbin

int *rccRawClustPtr, *rccGoodClustPtr;
int *rccRawClustPtrN, *rccGoodClustPtrN;
int *rccRawClustPtrS, *rccGoodClustPtrS;



//like:  rccRawClust[iBun*NPTBINS+iPt]++;
//and: rccRawClust[(isNorth+1)*NBUNCHES*NPTBINS+iBun*NPTBINS+iPt]++;



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
TH2F *hYieldByBunchAndPtNorth;
TH2F *hYieldByBunchAndPtSouth;
TH2F *hTightYieldByBunchAndPt;
TH2F *hTightYieldByBunchAndPtNorth;
TH2F *hTightYieldByBunchAndPtSouth;

TH1F *hRegionMassSpectrum[4][2];//[region][raw/summedenergy]
TH1F *hRegionClusts[4][3];//[region][raw/after loose/after tight]
TH1F *hRegionClustEcore[4][3];//[region][raw/after loose/after tight]
TH1F *hRegionClustMult[4][3];//[region][raw/after loose/after tight]
TH1F *hRegionClustDisp[4][3];//[region][raw/after loose/after tight]
TH1F *hRegionClustChi2[4][3];//[region][raw/after loose/after tight]
TH1F *hRegionClustE8e9[4][3];//[region][raw/after loose/after tight]



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
int GetRegion(int bx);
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

  TFile *rootin = NULL;
  rootin=TFile::Open(fullfile, "READONLY");
  if (rootin==NULL || rootin->IsZombie()|| !rootin->GetNkeys())
    return;

  ttree = (TTree *)rootin->Get("T");

  InitInTree();
  /*
  fillSets[1]="(";
  fillSets[1]+="(fill>=17217 && fill<=17232) ||";
  fillSets[1]+="(fill>=17238 && fill<=17240) ||";
  fillSets[1]+="(fill>=17247 && fill<=17256) ||";
  fillSets[1]+="(fill>=17263 && fill<=17276) ||";
  fillSets[1]+="(fill>=17284 && fill<=17297) ||";
  fillSets[1]+="(fill>=17302 && fill<=17305) ||";
  fillSets[1]+="(fill>=17308 && fill<=17317) ||";
  fillSets[1]+="(fill>=17328 && fill<=17333) ||";
  fillSets[1]+="(fill>=17338 && fill<=17382) ||";
  fillSets[1]+="(fill>=17391 && fill<=17396) ||";
  fillSets[1]+="(fill>=17403 && fill<=17407)";
  fillSets[1]+=")";
  fillNames[1]="29+69 Fills";
  */
  InitOutput(runnum, outputdir);
  InitDB(runnum);
  InitWarn(runnum);

  rccRun=runnum;

  //set rcc tree data:
  rccRun=runnum;
  rccNeve=ttree->GetEntries();
  //rccFill=;//this is set in InitDB
  rccTotRawClust=0;//accumulate these over the file
  rccTotGoodClust=0;//accumulate these over the file
  rccNbins=NPTBINS;
  rccNbounds=NPTBINS+1;
//to make sure I have the implicit ordering correctly, doing this out long-hand:
  for (int i=0;i<NBUNCHES;i++){
    for (int j=0;j<NPTBINS;j++){
      rccRawClust[i*NPTBINS+j]=0;
      rccGoodClust[i*NPTBINS+j]=0;
      rccRawClustN[i*NPTBINS+j]=0;
      rccGoodClustN[i*NPTBINS+j]=0;
      rccRawClustS[i*NPTBINS+j]=0;
      rccGoodClustS[i*NPTBINS+j]=0;
    }
  }
//these are filled into the tree by incrementing a pointer through them at fixed intervals.
  

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
    if (corrbunch<0) printf("corrected bunch less than zero (b=%d, shift=%d, corr=%d)!\n",bunch,spin_cont.GetCrossingShift(),corrbunch);
    //the bunch correction is wrong for some reason...
    
    even_or_odd = (corrbunch % 2); // 0 for even, 1 for odd

    
    int region=GetRegion(corrbunch);
    if (region!=-1){
	hRegionClusts[region][0]->Fill(nclus);
    }
    int nNominalClusters=0;
    int nTightClusters=0;
    for (int iclus = 0; iclus < nclus; iclus++) {

      int rccPtBin=0;
      for (rccPtBin=0;rccPtBin<NPTBINS;rccPtBin++){
	if (rccBounds[rccPtBin]<=pt[iclus] && pt[iclus]<rccBounds[rccPtBin+1]) break; //stop searching when we find the match;
      }
      int rccBinID=corrbunch*NPTBINS+rccPtBin;

      
      is_north = (feecore[iclus] < 288) ? 0 : 1;

 

	
      rccTotRawClust++;
      rccRawClust[rccBinID]++;
      if (is_north){
	rccRawClustN[rccBinID]++;
      }else{
	rccRawClustS[rccBinID]++;
      }
      

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

      //check cluster QA cuts:
      if (!PassesClusterCuts(iclus))
        continue;
      rspectrum_clustcuts->Fill(cluster_r);

      ptspectrum_raw[is_north][corrbunch % 2]->Fill(pt[iclus]);
      espectrum_raw[is_north][corrbunch % 2]->Fill(ecore[iclus]);


     //look for all possible pions in this arm, using the 0907.4832 paper cut definitions:
      float clusterE=ecore[iclus]/(1-e8e9[iclus]);
      TVector3 clusterVec(x[iclus],y[iclus],z[iclus]);
      TLorentzVector cluster4(clusterVec,clusterVec.Mag());//e=p, hence massless

      for (int pairclus = iclus+1; pairclus < nclus; pairclus++) {
	//printf("trying pair %d + %d\n",iclus,pairclus);
	bool pair_is_north = (feecore[pairclus] < 288) ? 0 : 1;
	if (pair_is_north!=is_north) continue; //skip if they're in different arms;
	if (!PassesClusterCuts(pairclus)) continue; //skip if it's not a good cluster;
      TVector3 pairVec(x[pairclus],y[pairclus],z[pairclus]);
      TLorentzVector pair4(pairVec,pairVec.Mag());//e=p, hence massless

      TLorentzVector sum4=pair4+cluster4;
	
	float pairE=ecore[pairclus]/(1-e8e9[pairclus]);
	float Egg=pairE+clusterE;
	float Eggcore=ecore[pairclus]+ecore[iclus];
	if (Egg<7 || Egg>17) continue; //skip if the energy is low or merged;
	float xrel=x[iclus]-x[pairclus];
	float yrel=y[iclus]-y[pairclus];
	float delr=sqrt(xrel*xrel+yrel*yrel);
	if (delr<3.5) continue;
	float alpha=abs(pairE-clusterE)/(pairE+clusterE);
	if (alpha>0.6) continue; //skip if the energy is too asymmetric
	//sin of half the opening angle:
	float sinth2=delr/2./zpos;//zpos should be the vertex, but I don't have that yet.  strictly speaking, the angle is:
	//sin(atan((delr/2)/zpos)), but even with very large separations of 40cm, this is nearly correct.

	float Mgg=sqrt(4*pairE*clusterE)*sinth2;
	float Mggcore=sqrt(4*ecore[pairclus]*ecore[iclus])*sinth2;
	splitClusterMgg=Mgg;
	splitClusterMggcore=Mggcore;
	splitClusterMvec=sum4.M();
	splitClusterPt=sum4.Pt();
	splitClusterTree->Fill();
	hRegionMassSpectrum[region][0]->Fill(Mggcore);//[region]
	hRegionMassSpectrum[region][1]->Fill(Mgg);//[region]
      }


      
      
      if (ecore[iclus] > 15.){
	//looking for merged clusters
	ptspectrum_ecut[is_north][corrbunch % 2]->Fill(pt[iclus]);
	espectrum_ecut[is_north][corrbunch % 2]->Fill(ecore[iclus]);
      }
      //old cut definition:
      // if (mult[iclus] <= 2)
      //	continue;
      //if (disp[iclus] < 0.0005)
      // 	continue;
      //if (chi2core[iclus] > 30.)
      //	continue;
      //if (e8e9[iclus] < 0.2)
      //	continue;

      bool nominalCut=false;
      bool tightCut=false;
      
      if (ecore[iclus]>15. &&
	  mult[iclus]>2 &&
	  disp[iclus]>0.0005 &&
	  chi2core[iclus] < 30. &&
	  e8e9[iclus]>0.2){
	nNominalClusters++;
	nominalCut=true;
      }

      if (ecore[iclus]>30 &&
	  mult[iclus]>2 &&
	  disp[iclus]>=0.0005 &&
	  chi2core[iclus] < 30. &&
	  e8e9[iclus]>0.25){
	nTightClusters++;
	tightCut=true;
      }

      if (region!=-1){
	for (int j=0;j<1+nominalCut+tightCut;j++){//assumes tightcut implies nominalcut.
	  //this is a per-event, not per-cluster variable:  hRegionClusts[region][j]->Fill();
	  hRegionClustEcore[region][j]->Fill(ecore[iclus]);
	  hRegionClustMult[region][j]->Fill(mult[iclus]);
	  hRegionClustDisp[region][j]->Fill(disp[iclus]);
	  hRegionClustChi2[region][j]->Fill(chi2core[iclus]);
	  hRegionClustE8e9[region][j]->Fill(e8e9[iclus]);
	}
      }


      int ix = mpcmap->getGridX(feecore[iclus]);
      int iy = mpcmap->getGridY(feecore[iclus]);

      //assign bunch tree variables:
       rccBunch=corrbunch;
      rccIx=ix;
      rccIy=iy;
      rccFeecore=feecore[iclus];
      rccMult=mult[iclus];
      rccX=x[iclus];
      rccY=y[iclus];
      rccVtx=zvtx;
      rccEcore=ecore[iclus];
      rccE9=ecore[iclus]/(1-e8e9[iclus]);
      rccE8=rccE9*e8e9[iclus];
      //rccE8e9=e8e9[iclus];
      rccDisp=disp[iclus];
      rccChi=chi2core[iclus];
      rccNorth=is_north;

      if (ecore[iclus]>10.){
	rccClusterTree->Fill();
      }
      
      
      if (!nominalCut && !tightCut)
	continue;

      //get bin and core coordinates:
      int ptbin = checkpt->Fill(pt[iclus]);
      if (ptbin<0){
	std::cout << "under/overflow.  not including event." << std::endl;
	continue;
      }
  
      
      if (nominalCut){
	if (feecore[iclus] == 66)
	  r_ch_066_post->Fill(cluster_r);
	if (feecore[iclus] == 468)
	  r_ch_468_post->Fill(cluster_r);
	rspectrum_candcuts->Fill(cluster_r);
	
	
	vtx[is_north][corrbunch % 2]->Fill(zvtx);
       
	pT_arr[is_north][ptbin-1][corrbunch]++;
	ptyield[ptbin - 1][is_north]->Fill(pt[iclus]);
	ptspectrum[is_north][corrbunch % 2]->Fill(pt[iclus]);
	espectrum[is_north][corrbunch % 2]->Fill(ecore[iclus]);

	toweryields[is_north]->Fill(ix, iy);
	hYieldByBunchAndPt->Fill(pt[iclus],corrbunch);
	if (is_north){
	  hYieldByBunchAndPtNorth->Fill(pt[iclus],corrbunch);
	}else{
	  hYieldByBunchAndPtSouth->Fill(pt[iclus],corrbunch);
	}

	rccTotGoodClust++;
	rccGoodClust[rccBinID]++;
	if (is_north){
	  rccGoodClustN[rccBinID]++;
	}else{
	  rccGoodClustS[rccBinID]++;
	}
      }
      if (tightCut){
	hTightYieldByBunchAndPt->Fill(pt[iclus],corrbunch);
	if (is_north){
	  hTightYieldByBunchAndPtNorth->Fill(pt[iclus],corrbunch);
	}else{
	  hTightYieldByBunchAndPtSouth->Fill(pt[iclus],corrbunch);
	}
      }
      
    }
    if (region!=-1){
	hRegionClusts[region][1]->Fill(nNominalClusters);
	hRegionClusts[region][2]->Fill(nTightClusters);
    }
  }
  rootin->Close();
  delete rootin;

  //fill our one-liner run variables:
  rccRunTree->Fill();
  //rccClusterTree->Fill();
  
  //fill our 120-line bunch variables:
  for (int b=0;b<120;b++){
    rccRawClustPtr=rccRawClust+(b*NPTBINS);
    rccGoodClustPtr=rccGoodClust+(b*NPTBINS);
    rccRawClustPtrN=rccRawClustN+(b*NPTBINS);
    rccGoodClustPtrN=rccGoodClustN+(b*NPTBINS);
    rccRawClustPtrS=rccRawClustS+(b*NPTBINS);
    rccGoodClustPtrS=rccGoodClustS+(b*NPTBINS);
    rccBunchTree->Fill();
  }
  
  End();//write and close output files.
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
  printf("creating yield histfile at %s\n",yieldfname.Data());

  rccRunTree=new TTree("runTree","per-run variables and totals");
  rccRunTree->Branch("run",&rccRun);
  rccRunTree->Branch("fill",&rccFill);
  rccRunTree->Branch("neve",&rccNeve);
  rccRunTree->Branch("nraw",&rccTotRawClust);
  rccRunTree->Branch("npi0",&rccTotGoodClust);
  rccRunTree->Branch("nptbins",&rccNbins);
  rccRunTree->Branch("nbounds",&rccNbounds);
  rccRunTree->Branch("ptbound",(void *)(&rccBounds[0]),"ptbound[nbounds]/F",32000);//convoluted way to force it to recognize the right version of the ambiguous signature.

  rccBunchTree=new TTree("bunchTree","per-bunch total raw and good clusters per det and sum");
  rccBunchTree->Branch("nRawByPt",rccRawClustPtr,Form("nRawByPt[%d]/F",NPTBINS));
  rccBunchTree->Branch("nPiByPt",rccGoodClustPtr,Form("nPiByPt[%d]/F",NPTBINS));
  rccBunchTree->Branch("nRawByPtN",rccRawClustPtrN,Form("nRawByPtN[%d]/F",NPTBINS));
  rccBunchTree->Branch("nPiByPtN",rccGoodClustPtrN,Form("nPiByPtN[%d]/F",NPTBINS));
  rccBunchTree->Branch("nRawByPtS",rccRawClustPtrS,Form("nRawByPtS[%d]/F",NPTBINS));
  rccBunchTree->Branch("nPiByPtS",rccGoodClustPtrS,Form("nPiByPtS[%d]/F",NPTBINS));

  rccClusterTree=new TTree("cTree","per-cluster info");
  rccClusterTree->Branch("fill",&rccFill);
  rccClusterTree->Branch("run",&rccRun);
  rccClusterTree->Branch("bunch",&rccBunch);
  rccClusterTree->Branch("ix",&rccIx);
  rccClusterTree->Branch("iy",&rccIy);
  rccClusterTree->Branch("x",&rccX);
  rccClusterTree->Branch("y",&rccY);
  rccClusterTree->Branch("vtx",&rccVtx);
  rccClusterTree->Branch("ecore",&rccEcore);
  rccClusterTree->Branch("feecore",&rccFeecore);
  rccClusterTree->Branch("e8",&rccE8);
  rccClusterTree->Branch("e9",&rccE9);
  rccClusterTree->Branch("mult",&rccMult);
  rccClusterTree->Branch("disp",&rccDisp);
  rccClusterTree->Branch("chi2core",&rccChi);
  rccClusterTree->Branch("north",&rccNorth);

  splitClusterTree=new TTree("piTree","pion clusters");
  rccClusterTree->Branch("M9",&splitClusterMgg);
  rccClusterTree->Branch("Mcore",&splitClusterMggcore);
  rccClusterTree->Branch("Mvec",&splitClusterMvec);
  rccClusterTree->Branch("pT",&splitClusterPt);

  Int_t sparsebins[4] = {2, 2, 120, 10}; // N/S, Even/Odd,crossing num., NPTBINS
  // spin patterns are in order: ++,+-,--,-+
  Double_t sparsebinsmin[4] = {-.5, -.5, -.5, 1.0};
  Double_t sparsebinsmax[4] = {1.5, 1.5, 119.5, 12};
  Double_t pt_limits[11] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};

  hYieldByBunchAndPt=new TH2F("hYieldByBunchAndPt","yield by bunch and pt",10,pt_limits,120,-0.5,119.5);
  hYieldByBunchAndPtNorth=new TH2F("hYieldByBunchAndPtNorth","yield by bunch and pt",10,pt_limits,120,-0.5,119.5);
  hYieldByBunchAndPtSouth=new TH2F("hYieldByBunchAndPtSouth","yield by bunch and pt",10,pt_limits,120,-0.5,119.5);

  hTightYieldByBunchAndPt=new TH2F("hTightYieldByBunchAndPt","Tight yield by bunch and pt",10,pt_limits,120,-0.5,119.5);
  hTightYieldByBunchAndPtNorth=new TH2F("hTightYieldByBunchAndPtNorth","Tight yield by bunch and pt",10,pt_limits,120,-0.5,119.5);
  hTightYieldByBunchAndPtSouth=new TH2F("hTightYieldByBunchAndPtSouth","Tight yield by bunch and pt",10,pt_limits,120,-0.5,119.5);


  TString regionname[]={"0<=bx<11","(29<=bx<40)||(69<=bx<80)","stable bxings","abort gap"};
  TString cutname[]={"raw","after loose cut","after tight cut"};
  for (int i=0;i<4;i++){
    hRegionMassSpectrum[i][0]=new TH1F(Form("hRegionMassSpectrum%d_0",i),
				   Form("Split Clusters (ecore) Mass Spectrum in %s;GeV/c^2",regionname[i].Data()),
				   300,0,3);//[region]
    hRegionMassSpectrum[i][1]=new TH1F(Form("hRegionMassSpectrum%d_1",i),
				   Form("Split Clusters (e9) Mass Spectrum in %sGeV/c^2",regionname[i].Data()),
				   300,0,3);//[region]
    for (int j=0;j<3;j++){
      hRegionClusts[i][j]=new TH1F(Form("hRegionClusts%d_%d",i,j),
				   Form("nClusters (%s) in %s",cutname[j].Data(),regionname[i].Data()),
				   10,-0.5,9.5);//[region][raw/after loose/after tight]
      hRegionClustEcore[i][j]=new TH1F(Form("hRegionClustEcore%d_%d",i,j),
				   Form("Cluster Core E (%s) in %s",cutname[j].Data(),regionname[i].Data()),
				   100,0,300);//[region][raw/after loose/after tight]
      hRegionClustMult[i][j]=new TH1F(Form("hRegionClustMult%d_%d",i,j),
				   Form("Cluster Multiplicity (%s) in %s",cutname[j].Data(),regionname[i].Data()),
				   10,-0.5,9.5);//[region][raw/after loose/after tight]
      hRegionClustDisp[i][j]=new TH1F(Form("hRegionClustDisp%d_%d",i,j),
				   Form("Cluster Dispersion (%s) in %s",cutname[j].Data(),regionname[i].Data()),
				       100,0.0,0.0001);//[region][raw/after loose/after tight]
      hRegionClustChi2[i][j]=new TH1F(Form("hRegionClustChi2%d_%d",i,j),
				   Form("Cluster Chi2 (%s) in %s",cutname[j].Data(),regionname[i].Data()),
				   100,0,200);//[region][raw/after loose/after tight]
      hRegionClustE8e9[i][j]=new TH1F(Form("hRegionClustE8e9%d_%d",i,j),
				   Form("e8e9 ratio (%s) in %s",cutname[j].Data(),regionname[i].Data()),
				   50,0,1.0);//[region][raw/after loose/after tight]
    }
  }



  
  
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
  rccFill=fillnum;

  rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", n_runnum);
  mpcmap = MpcMap::instance();


  printf("Not loading bunch status from pmontu.  Assuming all bunches good -- will be excluded later.\n");
  
  for (int i=0;i<120;i++){
    isBunchBad[i]=false;
  }

  return;

  int this_is_not_ever_reached=false;
  assert(this_is_not_ever_reached);
  
  // Populate array that tells which bunches to exclude because they fail QA
  // checks
  char sname[200];

  sprintf(sname, "/phenix/spin2/pmontu/offline/analysis/pmontu/"
                 "relative_luminosity/SpinDB/good_bunches/%d.txt",
          n_runnum);
  printf("loading good bunches list from %s\b",sname);
  ifstream bstatusfile(sname);
  if (!bstatusfile.good())
    cout << "Error! Can't open bunch status file. All bunches defaulted to GOOD" << endl;

  
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

int GetRegion(int bx){
  if (bx<11) return 0;
  if (bx>=29 && bx<40) return 1;
  if (bx>=69 && bx<80) return 1;
  if (bx>110) return 3;
  return 2;
}

bool PassesClusterCuts(int iclus) {
  //checks fiducial cuts on a cluster candidate.

  
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
  ifstream tdcoverfile(tdcfilename.Data());//is this still applicable in run13?  the electronics changed!
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
  hYieldByBunchAndPtNorth->Write();
  hYieldByBunchAndPtSouth->Write();
  hTightYieldByBunchAndPt->Write();
  hTightYieldByBunchAndPtNorth->Write();
  hTightYieldByBunchAndPtSouth->Write();
  rccBunchTree->Write();
  rccRunTree->Write();
  rccClusterTree->Write();
  splitClusterTree->Write();
  for (int i=0;i<4;i++){
    hRegionMassSpectrum[i][0]->Write();
    hRegionMassSpectrum[i][1]->Write();
    for (int j=0;j<3;j++){
      hRegionClusts[i][j]->Write();
      hRegionClustEcore[i][j]->Write();
      hRegionClustMult[i][j]->Write();
      hRegionClustDisp[i][j]->Write();
       hRegionClustChi2[i][j]->Write();
      hRegionClustE8e9[i][j]->Write();
      			
    }
  }


  

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

 

  for (std::size_t ibunch = 0; ibunch < 120; ibunch++) {
    for (int ipt = 0; ipt < NPTBINS; ipt++) {
      for (int iarm = 0; iarm <= 1; iarm++) {
	ptYields[iarm][ipt] = pT_arr[iarm][ipt][ibunch];
      }
    }
    t->Fill();
  }
  t->Write();
  t->ResetBranchAddresses();

   histfile->Close();
  delete histfile;
  return;
}
