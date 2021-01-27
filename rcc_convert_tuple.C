
//adjusting globally mis-aligned bunch IDs:
const int MASTER_BUNCH_OFFSET=0; //for the old march2020 set it is -5.  for new dec2020 set it is 0.

//this value needs to be added to the lumi bunch ID to get the corresponding MPC histogram bunch ID

void rcc_convert_tuple(){
  //this code produces a single master ntuple for basic asymmetry analysis using two sources:
  //1) histograms of yield vs pt and bin
  //2) the uLumiXL tuple from other code.

  TFile *uLumiFile=TFile::Open("uLumi.ttree.root","READ");
  TFile *uLumiXLfile=TFile::Open("uLumiXL.ttree.root","READ");
  TFile *uBunchFile=TFile::Open("uBunch.ttree.root","READ");
  TString yieldDir="yields2021";
  TString outputFileBase="uPiLumi2021";//"uPiLumi";

  
  // uLumi=(TTree*)uLumiFile->Get("uLumi");
  TTree *uLumiXL=(TTree*)uLumiXLfile->Get("uLumiXL");
  TTree *uBunch=(TTree*)uBunchFile->Get("uBunch");
  
  //prepare our output tree:
  TFile *uPiLumiFile=TFile::Open(Form("%s.ttree.root",outputFileBase.Data()),"RECREATE");
  TTree *uPiLumiBinning=new TTree("uPiLumiBinning","bin boundary info for yields in uPiLumi");
  int ubin_bin; uPiLumiBinning->Branch("bin",&ubin_bin);
  float ubin_ptlow; uPiLumiBinning->Branch("ptlow",&ubin_ptlow);
  float ubin_pthigh; uPiLumiBinning->Branch("pthigh",&ubin_pthigh);
  const int nptbins=10;
  const double pt_limits[] = {1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 12};
  double ptcenter[nptbins];
  for (int i=0;i<nptbins;i++){
    ubin_bin=i;
    ubin_ptlow=pt_limits[i];
    ubin_pthigh=pt_limits[i+1];
    uPiLumiBinning->Fill();
    ptcenter[i]=0.5*(ubin_ptlow+ubin_pthigh);
  }
 
  
  TTree *uPiLumi=new TTree("uPiLumi","per-bunch pion yields and other data needed to generate asymmetry");
  int u_fill; uPiLumi->Branch("fill",&u_fill);
  int u_run; uPiLumi->Branch("run",&u_run);
  int u_pat; uPiLumi->Branch("pat",&u_pat);
  double u_bpol; uPiLumi->Branch("bpol",&u_bpol);
  double u_ypol; uPiLumi->Branch("ypol",&u_ypol);
  double u_bpol_err; uPiLumi->Branch("bpol_err",&u_bpol_err);
  double u_ypol_err; uPiLumi->Branch("ypol_err",&u_ypol_err);
  
  int u_bunch; uPiLumi->Branch("bunch",&u_bunch);
  int u_bspin; uPiLumi->Branch("bspin",&u_bspin);
  int u_yspin; uPiLumi->Branch("yspin",&u_yspin);
  float u_yield[nptbins];
  float u_tightyield[nptbins];
  for (int i=0;i<nptbins;i++){
    uPiLumi->Branch(Form("yield%d",i),&(u_yield[i]));
    uPiLumi->Branch(Form("tightyield%d",i),&(u_tightyield[i]));
  }

  double u_zdc;uPiLumi->Branch("zdc",&u_zdc);
  double u_zdc_err;uPiLumi->Branch("zdc_err",&u_zdc_err);
  double u_bbc;uPiLumi->Branch("bbc",&u_bbc);
  double u_bbc_err;uPiLumi->Branch("bbc_err",&u_bbc_err);

  //check that the uBunch and uLumiXL data matches in terms of spin orientations.
  uLumiXL->Draw("run","bunch==0","goff");
  int nRuns=uLumiXL->GetSelectedRows();
  vector<int> runlist;
  printf("uLumiXL has %d runs\n",nRuns);
  for (int i=0;i<nRuns;i++){ //loop over all runs that pass our lumi-level cut:
    runlist.push_back(uLumiXL->GetVal(0)[i]);
  }




  
  TH2F * hYield[3];//pointers to be used in the loop.


  
  for (int i=0;i<nRuns;i++){
    uLumiXL->Draw("(likemuz1>0):fill:run:pat:bpol:ypol:bpol_err:ypol_err:bunch:likemuz1+unlikemuz1:likemub1+unlikemub1:likemuz1_err+unlikemuz1_err:likemub1_err+unlikemub1_err",Form("run==%d",runlist[i]),"goff");
    uBunch->Draw("(bspin==yspin):bspin:yspin",Form("run==%d",runlist[i]),"goff");
    int nBunches_L=uLumiXL->GetSelectedRows();
    int nBunches_B=uBunch->GetSelectedRows();
    if (nBunches_L!=nBunches_B){
      printf("number of bunches differs for run %d, (L=%d, B=%d)\n",runlist[i],nBunches_L,nBunches_B);
    }

    //load variables that won't change across the run:
    u_fill=uLumiXL->GetVal(1)[0];
    u_run=uLumiXL->GetVal(2)[0];
    u_pat=uLumiXL->GetVal(3)[0];
    u_bpol=uLumiXL->GetVal(4)[0];
    u_ypol=uLumiXL->GetVal(5)[0];
    u_bpol_err=uLumiXL->GetVal(6)[0];
    u_ypol_err=uLumiXL->GetVal(7)[0];

    
    //load the appropriate histogram, if available:
   TFile *yieldfile=NULL;
   yieldfile=TFile::Open(Form("./%s/%d.MPC.yields.rcc.hist.root",yieldDir.Data(),runlist[i]),"READ");
    if (yieldfile==NULL || yieldfile->IsZombie()|| !yieldfile->GetNkeys()){
      printf("couldn't find yields for run %d. skipping.\n",runlist[i]);
      continue;
    }
    hYield[0]=(TH2F*)yieldfile->Get("hYieldByBunchAndPt");
    hYield[1]=(TH2F*)yieldfile->Get("hTightYieldByBunchAndPt");
    hYield[2]=(TH2F*)yieldfile->Get("hYieldByBunchAndPtSouth");


    
    for (int j=0;j<nBunches_L;j++){
      int l_parity=uLumiXL->GetVal(0)[j];
      int b_parity=uBunch->GetVal(0)[j];
      u_bunch=uLumiXL->GetVal(8)[j];
      u_bspin=uBunch->GetVal(1)[j];
      u_yspin=uBunch->GetVal(2)[j];
      u_zdc=uLumiXL->GetVal(9)[j];
      u_bbc=uLumiXL->GetVal(10)[j];
      u_zdc_err=uLumiXL->GetVal(11)[j];
      u_bbc_err=uLumiXL->GetVal(12)[j];
      for (int k=0;k<nptbins;k++){
	int shifted_bunch=(120+u_bunch+MASTER_BUNCH_OFFSET)%120;
	int bin=hYield[0]->FindBin(ptcenter[k],shifted_bunch);
	u_yield[k]=hYield[0]->GetBinContent(bin);
	u_tightyield[k]=hYield[1]->GetBinContent(bin);
      }
      uPiLumi->Fill();
      
      if (l_parity!=b_parity){
	printf("run %d bunch %d parity differs (L=%d, B=%d)\n",runlist[i],j,l_parity,b_parity);
      }
    }
    yieldfile->Close();

  }
  uPiLumiFile->cd();
  printf("all values match unless you saw a complaint.  Phew.\n");
  uPiLumiBinning->Write();
  uPiLumi->Write();
  uPiLumiFile->Close();

  return;
}
