//this code is meant to do a basic qa to confirm that the crossing shift is correctly applied.
//it generates a pdf with one page for every few runs

void rcc_draw_raw_fillplots(){
  //a little booster program to plot raw spectra in the runs/fills, grouped conveniently by defineable sets of fills
  
  char outputfilename[100];
  char indexName[100];//variable to use to iterate over dataset, eg bunch or fill
  sprintf(outputfilename,"spectra2021_byfill_0good.pdf");
  sprintf(indexName,"fill");//"run");//generally only fill or run are useful terms here

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
  fillSets[0]="(";
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
  TCanvas *c=new TCanvas("cme","cme",800,1200);
  TLegend *leg;
  //TCanvas *cjunk=new TCanvas("cj","cj",80,80);
  c->Divide(2,3);

  for (int i=0;i<nIndices;i++){
    //int ipad=indicesDrawn%6+1;

    int thisIndex=indexList[i];
    TFile *runfile=NULL;

    //clear all our hist sets:
      if (!isFirstFile){
	//if we have already made the hists, clear them and assume the dimensions won't change file to file.
       for (int h=0;h<nHistSets;h++){
	 for (int n=0;n<histName[h].size();n++){
	   (histSet[h])[n]->Reset();
	 }
       }
      }
    
    //get all runs that share this index:
    uLumi->Draw("run",indexCut[i],"goff");
    // uLumi->Draw("run",Form("%s==%d",indexName,thisIndex),"goff");
    int nRuns=uLumi->GetSelectedRows();
    for (int j=0;j<nRuns;j++){
      int thisRun=uLumi->GetVal(0)[j];
      
      runfile=TFile::Open(Form("./yields_a2021/%d.MPC.yields.rcc.hist.root",thisRun),"READ");
      if (runfile==NULL || runfile->IsZombie()|| !runfile->GetNkeys()){
	printf("couldn't find yields for run %d. skipping.\n", thisRun);
	continue;
      }

      if (isFirstFile){
	//if this is our first successful open, load all our hists from what we find in the file.
       for (int h=0;h<nHistSets;h++){
	 for (int n=0;n<histName[h].size();n++){
	   scratch->cd();
	   TH1F *hTemp=NULL;
	   hTemp=(TH1F*)(runfile->Get(histName[h][n].Data()));
	   if (hTemp==NULL){
	     printf("couldn't find %s in file %s\n",histName[h][n].Data(), Form("./yields_a2021/%d.MPC.yields.rcc.hist.root",thisRun));
	     return;
	   }
	   histSet[h].push_back( (TH1F*)(hTemp->Clone(Form("hLocal%d_%d",h,n))) );
	   (histSet[h])[n]->SetLineColor(kBlack+n);
	   (histSet[h])[n]->SetTitle(Form("%d:%s",thisIndex, (histSet[h])[n]->GetTitle()));
	 }
       }
       isFirstFile=false;
      } else {
	//if this is not our first successful open, just Add these contents to the hists.  But do rename:
	for (int h=0;h<nHistSets;h++){
	  for (int n=0;n<histName[h].size();n++){
	    TH1F *hTemp=NULL;
	    hTemp=(TH1F*)(runfile->Get(histName[h][n].Data()));
	    if (hTemp==NULL){
	     printf("couldn't find %s in file %s\n",histName[h][n].Data(), Form("./yields_a2021/%d.MPC.yields.rcc.hist.root",thisRun));
	     continue;
	   }
	    (histSet[h])[n]->SetTitle(Form("%s:%s",indexDisplayName[i].Data(), hTemp->GetTitle()));

	    (histSet[h])[n]->Add(hTemp);
	 }
       }
      }

      runfile->Close();
      //           histSet[0][0]->DrawNormalized("hist"); return;


    }

    printf("drawing plots from %s %d to canvas\n", indexName,thisIndex);
    //hYield->SetLineColor(kBlack);
    bool hasGoodPlots=false;
    for (int h=0;h<nHistSets;h++){
      c->cd(h+1)->Clear();
      if ((histSet[h])[0]->Integral()>0){
	(histSet[h])[0]->DrawNormalized("hist");
	hasGoodPlots=true;
      }
      for (int n=1;n<histName[h].size();n++){
	if ((histSet[h])[n]->Integral()>0){
	  (histSet[h])[n]->DrawNormalized("same,hist");
	  hasGoodPlots=true;
	}
      }
   
    }
    
    if (hasGoodPlots){
      c->cd(1)->SetLogy();
      c->cd(1)->BuildLegend();
      if (isFirstPage){
	c->Print(Form("%s(",outputfilename),"pdf");
	isFirstPage=false;
      }
      else{
	c->Print(outputfilename,"pdf");
      }
    }
      //clear off the canvas so we don't accidentally save stale data later:
      for (int iclear=1;iclear<7;iclear++){
	c->cd(iclear)->Clear();
      }
  }

  	c->Print(Form("%s)",outputfilename),"pdf");


  return;
}
