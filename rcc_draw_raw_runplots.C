//this code is meant to do a basic qa to look at low mass pions and other things, per run.

void rcc_draw_raw_runplots(){
  //a little booster program to plot raw spectra in the runs/fills, grouped conveniently by defineable sets of fills
  
  char outputfilename[100];
  char indexName[100];//variable to use to iterate over dataset, eg bunch or fill
  sprintf(outputfilename,"spectra2021_byrun.pdf");
  sprintf(indexName,"fill");//"run");//generally only fill or run are useful terms here


  char yieldDirectory[100];
  sprintf(yieldDirectory,"yields2021_2021.02.17");

  TFile *scratch=TFile::Open("scratch.hist.root","RECREATE");

  //define what plots you want to see in stacks which we will draw:
  int nHistSets=6;
  
  vector<TString> histName[nHistSets];
  vector<TH1F*> histSet[nHistSets];
  vector<TString> regionName;
  regionName.push_back("0<=bx<11");
  regionName.push_back("(29<=bx<40)||(69<=bx<80)");
  regionName.push_back("stable bxings");
  regionName.push_back("abort gap");
  int nRegions=4;
  for (int i=0;i<nRegions;i++){
    //for now, all regions, just the raw plots.
    histName[0].push_back(Form("hRegionClusts%d_0",i));
    histName[1].push_back(Form("hRegionClustEcore%d_0",i));
    histName[2].push_back(Form("hRegionClustMult%d_0",i));
    histName[3].push_back(Form("hRegionClustDisp%d_0",i));
    histName[4].push_back(Form("hRegionClustChi2%d_0",i));
    histName[5].push_back(Form("hRegionClustE8e9%d_0",i));
 
  }
  
  
  TFile *uLumiFile=TFile::Open("uLumi.ttree.root","READ");
  TFile *uLumiXLfile=TFile::Open("uLumiXL.ttree.root","READ");
  TFile *uBunchFile=TFile::Open("uBunch.ttree.root","READ");

  
  TTree *uLumi=(TTree*)uLumiFile->Get("uLumi");
  TTree *uLumiXL=(TTree*)uLumiXLfile->Get("uLumiXL");
  //uBunch=(TTree*)uBunchFile->Get("uBunch");

  //generate the list of unique indices (assuming they're sequential and not shuffled randomly):
  vector<int> indexList;
  vector<TString> indexCut;
  vector<TString> indexDisplayName;
  int nIndices=0;
  /*
  uLumi->Draw(indexName,"1","goff");
  uLumi->SetLineColor(kRed);
  //int nIndices=uLumi->GetSelectedRows();
  for (int i=0;i<nIndices;i++){
    int newIndex=uLumi->GetVal(0)[i];
    if (i==0){//always push the first value to the list.
      indexList.push_back(newIndex);
    } else if(indexList.back()!=newIndex){//only push other values if they're not already in there (assume sorted)
      indexList.push_back(newIndex);
indexDisplayName.push_back(Form(%d,newIndex));
    }
  }
  */



  int nFillSets=2;
  TString fillSets[nFillSets];
  TString fillNames[nFillSets];
  fillSets[0]="(1 ||";
  //fillSets[0]+="(fill>=17410 && fill<=17415) ||";
  fillSets[0]+="fill==17417 ||";
  fillSets[0]+="fill==17429 ||";
  // fillSets[0]+="(fill>=17431 && fill<=17434) ||";
  // fillSets[0]+="(fill>=17439 && fill<=17451) ||";
  fillSets[0]+="fill==17455 ||";
  fillSets[0]+="(fill>=17474 && fill<=17479) ||";
  fillSets[0]+="(fill>=17486 && fill<=17488) ||";
  fillSets[0]+="(fill>=17492 && fill<=17514) ||";
  fillSets[0]+="fill==17518 ||";
  fillSets[0]+="(fill>=17520 && fill<=17524) ||";
  fillSets[0]+="(fill>=17530 && fill<=17533) ||";
  fillSets[0]+="(fill>=17536 && fill<=17538) ||";
  fillSets[0]+="(fill>=17544 && fill<=17545) ||";
  fillSets[0]+="fill==17550 ||";
  fillSets[0]+="(fill>=17558 && fill<=17561) ||";
  fillSets[0]+="(fill>=17568 && fill<=17573) ||";
  fillSets[0]+="(fill>=17579 && fill<=17601)";
  fillSets[0]+=")";
  fillNames[0]="No-Gap Fills";

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

  
  indexDisplayName.push_back(fillNames[0]);
  indexCut.push_back(fillSets[0]);
  indexList.push_back(0);
    indexDisplayName.push_back(fillNames[1]);
  indexCut.push_back(fillSets[1]);
  indexList.push_back(0);

  nIndices=indexList.size();
  // for (int i=0;i<nIndices;i++){
  //  printf("i=%d, index=%d\n",i,indexList[i]);
  //}



  
  bool isFirstFile=true;
  bool isFirstPage=true;
  int indicesDrawn=0;
  TCanvas *cbase=new TCanvas("cme","cme",850,1100);
  TLegend *leg;
  TLatex tex;
  TH1F* htemp;
  //TCanvas *cjunk=new TCanvas("cj","cj",80,80);
  TPad *ctitle=new TPad("ctitle","ctitle",0.0,0.9,1.0,1.0);
  TPad *c=new TPad("cdata","cdata",0.0,0.0,1.0,0.9);
  c->Divide(4,6);
  cbase->cd();
  ctitle->Draw();
  c->Draw();

  for (int i=0;i<nIndices && i<1;i++){
    //int ipad=indicesDrawn%6+1;

    int thisIndex=indexList[i];
    TFile *runfile=NULL;
    
    //get all runs that share this index:
    uLumi->Draw("run:fill",indexCut[i],"goff");
    // uLumi->Draw("run",Form("%s==%d",indexName,thisIndex),"goff");
    int nRuns=uLumi->GetSelectedRows();
    for (int j=0;j<nRuns && j<4;j++){
      int thisRun=uLumi->GetVal(0)[j];
      int thisFill=uLumi->GetVal(1)[j];
      
      runfile=TFile::Open(Form("./%s/%d.MPC.yields.rcc.hist.root",yieldDirectory,thisRun),"READ");
      if (runfile==NULL || runfile->IsZombie()|| !runfile->GetNkeys()){
	printf("couldn't find yields for run %d. skipping.\n", thisRun);
	continue;
      }
      printf("found yields for run %d\n",thisRun);
      TTree *cTree=(TTree*)runfile->Get("cTree");
      TTree *piTree=(TTree*)runfile->Get("piTree");
      ctitle>cd();
      ctitle->Clear();
      tex.SetTextSize(0.25);
      tex.DrawLatex(0.1,0.5,Form("Run=%d Fill=%d",thisRun,thisFill));
      
      c->cd(1);
      cTree->Draw("ix:iy","north","colz");
      c->cd(2);
      cTree->Draw("ix:iy","!north","colz");
       c->cd(3);
      cTree->Draw("x:y","north","colz");
      c->cd(4);
      cTree->Draw("x:y","!north","colz");
      c->cd(5);
      cTree->Draw("ix:iy","north && ecore>1000","colz");
      c->cd(6);
      cTree->Draw("ix:iy","!north && ecore>1000","colz");
       c->cd(7);
      cTree->Draw("x:y","north && ecore>1000","colz");
      c->cd(8);
      cTree->Draw("x:y","!north && ecore>1000","colz");
     c->cd(9);
      cTree->Draw("feecore:ecore","1","colz");
     c->cd(10);
      cTree->Draw("feecore");
      htemp=(TH1F*)(c->cd(10)->GetPrimitive("htemp"));
      tex.SetTextSize(0.08);
     tex.DrawLatex(0.1,htemp->GetMaximum()*0.9,Form("MaxBin=%1.1f",htemp->GetBinLowEdge(htemp->GetMaximumBin())));
      
     c->cd(11);
     cTree->Draw("feecore","ecore>1000");
     htemp=(TH1F*)(c->cd(11)->GetPrimitive("htemp"));
      tex.SetTextSize(0.08);
     tex.DrawLatex(0.1,htemp->GetMaximum()*0.9,Form("MaxBin=%1.1f",htemp->GetBinLowEdge(htemp->GetMaximumBin())));
        runfile->Close();
      //           histSet[0][0]->DrawNormalized("hist"); return;

      if (isFirstPage){
	cbase->Print(Form("%s(",outputfilename),"pdf");
	isFirstPage=false;
      }
      else{
	cbase->Print(outputfilename,"pdf");
      }
    }

  }
  cbase->Clear();
   	cbase->Print(Form("%s)",outputfilename),"pdf");
	return;
}
