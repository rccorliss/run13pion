
void rcc_convert_tuple(){
  //this code produces a single master ntuple for basic asymmetry analysis using two sources:
  //1) histograms of yield vs pt and bin
  //2) the uLumiXL tuple from other code.

  TFile *uLumiFile=TFile::Open("uLumi.ttree.root","READ");
  TFile *uLumiXLfile=TFile::Open("uLumiXL.ttree.root","READ");
  TFile *uBunchFile=TFile::Open("uBunch.ttree.root","READ");

  
  // uLumi=(TTree*)uLumiFile->Get("uLumi");
  TTree *uLumiXL=(TTree*)uLumiXLfile->Get("uLumiXL");
  double u_bpol; uLumiXL->SetBranchAddress("bpol",&u_bpol);
  double u_ypol; uLumiXL->SetBranchAddress("ypol",&u_ypol);
  double u_bpol_err; uLumiXL->SetBranchAddress("bpol_err",&u_bpol_err);
  double u_ypol_err; uLumiXL->SetBranchAddress("ypol_err",&u_ypol_err);

  double u_zdclumi_like; uLumiXL->SetBranchAddress("likemuz1",&u_zdclumi_like);
  double u_zdclumi_unlike; uLumiXL->SetBranchAddress("unlikemuz1",&u_zdclumi_unlike);


  
  TTree *uBunch=(TTree*)uBunchFile->Get("uBunch");
  int ub_fill; uBunch->SetBranchAddress("fill",&ub_fill);
  int ub_run; uBunch->SetBranchAddress("run",&ub_run);
  int ub_pat; uBunch->SetBranchAddress("pat",&ub_pat);
  int ub_bunch; uBunch->SetBranchAddress("bunch",&ub_bunch);
  int ub_bspin; uBunch->SetBranchAddress("bspin",&ub_bspin);
  int ub_yspin; uBunch->SetBranchAddress("yspin",&ub_yspin);
  //blah.



  //check that the uBunch and uLumiXL data matches in terms of spin orientations.
  uLumiXL->Draw("run","bunch==0","goff");
  int nRuns=uLumiXL->GetSelectedRows();
  vector<int> runlist;
  printf("uLumiXL has %d runs\n",nRuns);
  for (int i=0;i<nRuns;i++){ //loop over all runs that pass our lumi-level cut:
    runlist.push_back(uLumiXL->GetVal(0)[i]);
  }

  for (int i=0;i<nRuns;i++){
    uLumiXL->Draw("(likemuz1>0)",Form("run==%d",runlist[i]),"goff");
    uBunch->Draw("(bspin==yspin)",Form("run==%d",runlist[i]),"goff");
    int nBunches_L=uLumiXL->GetSelectedRows();
    int nBunches_B=uBunch->GetSelectedRows();
    if (nBunches_L!=nBunches_B){
      printf("number of bunches differs for run %d, (L=%d, B=%d)\n",runlist[i],nBunches_L,nBunches_B);
    }
    for (int j=0;j<nBunches_L;j++){
      int parity_L=uLumiXL->GetVal(0)[j];
      int parity_B=uBunch->GetVal(0)[j];
      if (parity_L!=parity_B){
	printf("run %d bunch %d parity differs (L=%d, B=%d)\n",runlist[i],j,parity_L,parity_B);
      }
    }
  }
        printf("all values match unless you saw a complaint.  Phew.\n");


  return;
}
