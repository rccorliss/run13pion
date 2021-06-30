
void CalcAsymAndErr(double *asym, double *asym_err,
		    double bpol, double bpol_err,
		    double ypol, double ypol_err,
		    double unlike_signal, double unlike_signal_err,
		    double like_signal, double like_signal_err,
		    double unlike_lumi, double unlike_lumi_err,
		    double like_lumi, double like_lumi_err);

float SumArray(double *arr, int n){
  float ret=0;
  for (int i=0;i<n;i++) ret+=arr[i];
  return ret;
}

void PlotFullAverageAsym();
void   PlotAsymByRun();
void   PlotAsymByFill();


const int MAXCUTS=100;
int nSets=0;
TString cut[MAXCUTS],cutName[MAXCUTS];

void LoadEvenlySpacedSegments(int nSegments){
  int nPatSets=2;
  TString patSets[]={"pat == 21 || pat == 24 || pat == 25 || pat == 28",
		     "pat == 22 || pat == 23 || pat == 26 || pat == 27"};
  TString patNames[]={"SSOO","OOSS"};
  int nDivSets=nSegments;
  TString divSets[nDivSets];//={"abs(bunch-10)<10","abs(bunch-35)<15","abs(bunch-75)<15","abs(bunch-105)<15"};
  TString divNames[nDivSets];//={"0<bun<30","30<bun<60","60<bun<90","90<bun<120"};
  int divSpacing=120/nDivSets;
  for (int i=0;i<nDivSets;i++){
    divSets[i]=Form("abs(bunch-%d)<%d",divSpacing/2+divSpacing*i,divSpacing);
    divNames[i]=Form("%d<bunch<%d",divSpacing*i,divSpacing*(i+1));
  }
  int cuti=0;
  for (int i=0;i<nPatSets;i++){
    for (int j=0;j<nDivSets;j++){
      //don't dead reckon.  It goes wrong if I fuss with the limits (ie skip the first one): cuti=i*nDivSets+j;
      if (cuti>=MAXCUTS){
	printf("tried to make too many divisions!  What do you need more than %d for?\n",MAXCUTS);
      }
      cut[cuti]=Form("(%s)&&(%s)",patSets[i].Data(),divSets[j].Data());
      cutName[cuti]=Form("%s %s",patNames[i].Data(),divNames[j].Data());
      printf("Cut %d:  \"%s\": %s\n",cuti,cutName[cuti].Data(),cut[cuti].Data());
      cuti++;
    }
  }
  nSets=cuti;
  printf("total cuts=%d\n",nSets);

  return;
}

void LoadTailoredSegments(){
  int nPatSets=2;
  TString patSets[]={"pat == 21 || pat == 24 || pat == 25 || pat == 28",
		     "pat == 22 || pat == 23 || pat == 26 || pat == 27"};
  TString patNames[]={"SSOO","OOSS"};
  int nDivSets=7;
  TString divSets[nDivSets];
  TString divNames[nDivSets];
  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  for (int i=0;i<nDivSets;i++){
    divSets[i]=Form("%d<=bunch && bunch<%d",divBounds[i],divBounds[i+1]);
    divNames[i]=Form("%d<=bunch<%d",divBounds[i],divBounds[i+1]);
  }

  int cuti=0;
  for (int i=1;i<nPatSets;i++){
    for (int j=0;j<nDivSets;j++){
      //don't dead reckon.  It goes wrong if I fuss with the limits (ie skip the first one): cuti=i*nDivSets+j;
      if (cuti>=MAXCUTS){
	printf("tried to make too many divisions!  What do you need more than %d for?\n",MAXCUTS);
      }
      cut[cuti]=Form("(%s)&&(%s)",patSets[i].Data(),divSets[j].Data());
      cutName[cuti]=Form("%s %s",patNames[i].Data(),divNames[j].Data());
      printf("Cut %d:  \"%s\": %s\n",cuti,cutName[cuti].Data(),cut[cuti].Data());
      cuti++;
    }
  }

    nSets=cuti;
    printf("total cuts=%d\n",nSets);
  return;
}



void LoadFillTailoredSegments(){
  int nPatSets=1;
  TString patSets[]={"pat == 21 || pat == 24 || pat == 25 || pat == 28",
		     "pat == 22 || pat == 23 || pat == 26 || pat == 27"};
  TString patNames[]={"SSOO","OOSS"};

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
  
  vector<TString> divSets[nFillSets];
  vector<TString> divNames[nFillSets];
  /*
  divSets[0].push_back("bunch<11"); divNames[0].push_back("0<=bx<11");
  divSets[0].push_back("bunch>10 && bunch<28"); divNames[0].push_back("11<=bx<28");
  divSets[0].push_back("bunch<11"); divNames[0].push_back("0<=bx<11");
  divSets[0].push_back("bunch<11"); divNames[0].push_back("0<=bx<11");
  divSets[0].push_back("bunch<11"); divNames[0].push_back("0<=bx<11");
  TString divSets[nDivSets];
  TString divNames[nDivSets];
  */
  int nDivs=6;
  //old bounds:  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  int divBounds[]={0,14,28,41,68,81,110,120};//lower bound is included, upper bound is excluded.
  for (int i=0;i<nDivs;i++){
    divSets[0].push_back(Form("%d<=bunch && bunch<%d",divBounds[i],divBounds[i+1]));
    divNames[0].push_back(Form("%d<=bunch<%d",divBounds[i],divBounds[i+1]));
  }
  if (0){
    int nAutoDivs=7;
    for (int i=0;i<nAutoDivs;i++){
      int lowbound=11+(100/nAutoDivs*i);
      int highbound=lowbound+(100/nAutoDivs);
      divSets[0].push_back(Form("bunch>=%d && bunch<%d",lowbound,highbound));
      divNames[0].push_back(Form("%d<=bx<%d",lowbound,highbound));
    }
  }
  //divSets[0].push_back("bunch>10 && bunch<21");
  // divNames[0].push_back("11<=bx<21");
  //divSets[0].push_back("bunch>20 && bunch<110 && !(bunch%2)");
  //divNames[0].push_back("11<=bx<110 even");ecore[iclus]>10. && 
  //divSets[0].push_back("bunch>20 && bunch<110 && (bunch%2)");
  //divNames[0].push_back("11<=bx<110 odd");
  /*
  int nDivSets=7;
  TString divSets[nDivSets];
  TString divNames[nDivSets];
  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  for (int i=0;i<nDivSets;i++){
    divSets[i]=Form("%d<=bunch && bunch<%d",divBounds[i],divBounds[i+1]);
    divNames[i]=Form("%d<=bunch<%d",divBounds[i],divBounds[i+1]);
  }
  */
   nFillSets=1;
  int cuti=0;
  for (int i=0;i<nPatSets;i++){
    for (int j=0;j<nFillSets;j++){
      for (int k=0;k<divSets[0].size();k++){
      //don't dead reckon.  It goes wrong if I fuss with the limits (ie skip the first one): cuti=i*nDivSets+j;
      if (cuti>=MAXCUTS){
	printf("tried to make too many divisions!  What do you need more than %d for?\n",MAXCUTS);
      }
      cut[cuti]=Form("(%s)&&(%s)&&(%s)",patSets[i].Data(),fillSets[j].Data(),divSets[0][k].Data());
      cutName[cuti]=Form("%s %s %s",patNames[i].Data(),fillNames[j].Data(), divNames[0][k].Data());
      printf("Cut %d:  \"%s\": %s \n",cuti,cutName[cuti].Data(),cut[cuti].Data());
      cuti++;
      }
    }
  }

    nSets=cuti;
    printf("total cuts=%d\n",nSets);
  return;
}


void LoadBunchSafeTailoredSegments(int selectPat=0, int selectFill=0){
  int nPatSets=2;
  TString patSets[]={"pat == 21 || pat == 24 || pat == 25 || pat == 28",
		     "pat == 22 || pat == 23 || pat == 26 || pat == 27"};
  TString patNames[]={"SSOO","OOSS"};

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
  
  vector<TString> divSets[nFillSets];
  vector<TString> divNames[nFillSets];

  /* as of may 2021 presentation:
  int nDivs=3;
  //old bounds:  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  int divBounds[]={0,40,80};//lower bound is included, upper bound is excluded.
  int divUpperBounds[]={11,65,111};//lower bound is included, upper bound is excluded.
  for (int i=0;i<nDivs;i++){
    divSets[0].push_back(Form("%d<=bunch && bunch<%d && (bunch%%2)",divBounds[i],divUpperBounds[i]));
    divNames[0].push_back(Form("%d<=bx<%d odd",divBounds[i],divUpperBounds[i]));
    divSets[0].push_back(Form("%d<=bunch && bunch<%d&& !(bunch%%2)",divBounds[i],divUpperBounds[i]));
    divNames[0].push_back(Form("%d<=bx<%d even",divBounds[i],divUpperBounds[i]));
  }
  */
  int nDivs=11;
  //old bounds:  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  int divBounds[]={0,40,80};//lower bound is included, upper bound is excluded.
  int divUpperBounds[]={11,65,111};//lower bound is included, upper bound is excluded.
  for (int i=0;i<nDivs;i++){
    divSets[0].push_back(Form("%d<=bunch && bunch<%d && (bunch%%2)",i*10,i*10+10));
    divNames[0].push_back(Form("%d<=bx<%d odd",i*10,i*10+10));
    //divSets[0].push_back(Form("%d<=bunch && bunch<%d&& !(bunch%%2)",i*10,i*10+10));
    //divNames[0].push_back(Form("%d<=bx<%d even",i*10,i*10+10));
  }
 
   nFillSets=2;
   int cuti=0;
   for (int i=0;i<nPatSets;i++){
     if (selectPat!=-1){
       i=selectPat;
     }
     for (int j=0;j<nFillSets;j++){
       if (selectFill!=-1){
	 j=selectFill;
       }
       for (int k=0;k<divSets[0].size();k++){
	 //don't dead reckon.  It goes wrong if I fuss with the limits (ie skip the first one): cuti=i*nDivSets+j;
	 if (cuti>=MAXCUTS){
	   printf("tried to make too many divisions!  What do you need more than %d for?\n",MAXCUTS);
	 }
	 cut[cuti]=Form("(%s)&&(%s)&&(%s)",patSets[i].Data(),fillSets[j].Data(),divSets[0][k].Data());
	 cutName[cuti]=Form("%s %s %s",patNames[i].Data(),fillNames[j].Data(), divNames[0][k].Data());
	 printf("Cut %d:  \"%s\": %s \n",cuti,cutName[cuti].Data(),cut[cuti].Data());
	 cuti++;
       }
     if (selectFill!=-1) break; //don't loop if we picked a particular set.
     }
     if (selectPat!=-1) break; //don't loop if we picked a particular set.

   }

    nSets=cuti;
    printf("total cuts=%d\n",nSets);
  return;
}



void LoadJune14Segments(int selectPat=0, int selectFill=0){
  int nPatSets=2;
  TString patSets[]={"pat == 21 || pat == 24 || pat == 25 || pat == 28",
		     "pat == 22 || pat == 23 || pat == 26 || pat == 27"};
  TString patNames[]={"SSOO","OOSS"};

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
  
  vector<TString> divSets[nFillSets];
  vector<TString> divNames[nFillSets];

  /*
  int nDivs=5;
  //old bounds:  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  int divEvens[]={0,0,35,68,75,110};
  int divOdds[]={0,28,35,68,75,110};

  for (int i=0;i<nDivs;i++){
    divSets[0].push_back(Form("%d<=bunch && bunch<%d && (bunch%%2)",divOdds[i],divOdds[i+1]));
    divNames[0].push_back(Form("%d<=bx<%d odd",divOdds[i],divOdds[i+1]));
    if (i>0){
    divSets[0].push_back(Form("%d<=bunch && bunch<%d&& !(bunch%%2)",divEvens[i],divEvens[i+1]));
    divNames[0].push_back(Form("%d<=bx<%d even",divEvens[i],divEvens[i+1]));
    }
  }
  */

      int nDivs=5;
  //old bounds:  int divBounds[]={0,11,29,40,69,80,111,120};//lower bound is included, upper bound is excluded.
  int divEvens[]={35,75};
  int divEvensHigh[]={68,110};
  int divOdds[]={0,35,75};
  int divOddsHigh[]={28,68,110};


  for (int i=0;i<2;i++){

    divSets[0].push_back(Form("%d<=bunch && bunch<%d&& !(bunch%%2)",divEvens[i],divEvensHigh[i]));
    divNames[0].push_back(Form("%d<=bx<%d even",divEvens[i],divEvensHigh[i]));
  }
  for (int i=0;i<3;i++){
    divSets[0].push_back(Form("%d<=bunch && bunch<%d && (bunch%%2)",divOdds[i],divOddsHigh[i]));
    divNames[0].push_back(Form("%d<=bx<%d odd",divOdds[i],divOddsHigh[i]));
  }
  TString fullset=Form("(%s)",divSets[0][0].Data());
  for (int i=1;i<divSets[0].size();i++){
    fullset=fullset+Form("||(%s)",divSets[0][i].Data());
  }
  divSets[0].push_back(fullset);
  divNames[0].push_back("all good bunches");
 
   nFillSets=2;
   int cuti=0;
   for (int i=0;i<nPatSets;i++){
     if (selectPat!=-1){
       i=selectPat;
     }
     for (int j=0;j<nFillSets;j++){
       if (selectFill!=-1){
	 j=selectFill;
       }
       for (int k=0;k<divSets[0].size();k++){
	 //don't dead reckon.  It goes wrong if I fuss with the limits (ie skip the first one): cuti=i*nDivSets+j;
	 if (cuti>=MAXCUTS){
	   printf("tried to make too many divisions!  What do you need more than %d for?\n",MAXCUTS);
	 }
	 cut[cuti]=Form("(%s)&&(%s)&&(%s)",patSets[i].Data(),fillSets[j].Data(),divSets[0][k].Data());
	 cutName[cuti]=Form("%s %s %s",patNames[i].Data(),fillNames[j].Data(), divNames[0][k].Data());
	 printf("Cut %d:  \"%s\": %s \n",cuti,cutName[cuti].Data(),cut[cuti].Data());
	 cuti++;
       }
     if (selectFill!=-1) break; //don't loop if we picked a particular set.
     }
     if (selectPat!=-1) break; //don't loop if we picked a particular set.

   }

    nSets=cuti;
    printf("total cuts=%d\n",nSets);
  return;
}



//useful things to get once and get out of the way:
TTree *uPiLumi, *uPiLumiBinning; //data and binning data.
TString yieldName="yields";//yield, yields, yieldn, tightyield.

int nBins; //number of ptbins
vector<double> ptmid,pterr;//bin centers and distances from center to edge


void rcc_asym_from_uPiLumi(){
  TFile *uPiLumiFile=TFile::Open("uPiLumi2021.06.ttree.root","READ");
  uPiLumi=(TTree*)uPiLumiFile->Get("uPiLumi");
  uPiLumiBinning=(TTree*)uPiLumiFile->Get("uPiLumiBinning");

  
  //get our binning:
  uPiLumiBinning->Draw("ptlow:pthigh","1","goff");
  nBins=uPiLumiBinning->GetSelectedRows();
  double *ptlow=uPiLumiBinning->GetVal(0);
  double *pthigh=uPiLumiBinning->GetVal(1);
  for (int i=0;i<nBins;i++){
    ptmid.push_back(0.5*(ptlow[i]+pthigh[i]));
    pterr.push_back(0.5*(ptlow[i]-pthigh[i]));
  }
  
  //PlotAsymByFill();
  PlotFullAverageAsym();
  return;
}

void   PlotFullAverageAsym(){
  //sums yields over all bunches in all runs and produces a single asym from that.


  //fill in the cut and cutname variables (defined globally, because this is messy) with the chosen way to divde the data into sets.
  //LoadFillTailoredSegments();
  //LoadBunchSafeTailoredSegments(-1,-1);
  //LoadEvenlySpacedSegments(10);
  LoadJune14Segments(-1,-1);
  
  //int nDivSets=10;
  //int nSets=8;//1*nDivSets;
  vector<double> likeYield[nSets];
  vector<double> unlikeYield[nSets];
  vector<double> asym[nSets],asymErr[nSets];
  vector<double> asym_byptbin[nBins],asymErr_byptbin[nBins];
  double rel[nSets],relErr[nSets];
  
  TString lumiMon="zdc";
  TString lumiMonErr2=Form("%s_err*%s_err",lumiMon.Data(),lumiMon.Data());

  //get the average polarizations:
  double bpol[nSets], bpolErr[nSets], ypol[nSets],ypolErr[nSets];
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("bpol:bpol_err:ypol:ypol_err",Form("(zdc>0)*(%s)",cut[i].Data()),"goff");
    bpol[i]=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    bpolErr[i]=SumArray(uPiLumi->GetVal(1),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    ypol[i]=SumArray(uPiLumi->GetVal(2),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
    ypolErr[i]=SumArray(uPiLumi->GetVal(3),uPiLumi->GetSelectedRows())/uPiLumi->GetSelectedRows();
  }
  
  //calc the rel lumi
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw(lumiMon.Data(),Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    int nLike=uPiLumi->GetSelectedRows();
    double zdcLike=SumArray(uPiLumi->GetVal(0),nLike);
    uPiLumi->Draw(lumiMon.Data(),Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    int nUnlike=uPiLumi->GetSelectedRows();
    double zdcUnlike=SumArray(uPiLumi->GetVal(0),nUnlike);
    //naively, the uncertainties in like an unlike are independent (though really they share the global scaler)
    uPiLumi->Draw(lumiMonErr2.Data(),Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");//square of error term
    double zdcLikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    uPiLumi->Draw(lumiMonErr2.Data(),Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double zdcUnlikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    rel[i]=(zdcLike/zdcUnlike);
    relErr[i]=rel[i]*sqrt((zdcLikeErr/zdcLike)*(zdcLikeErr/zdcLike)
			  +(zdcUnlikeErr/zdcUnlike)*(zdcUnlikeErr/zdcUnlike));
    //rel[i]*=1.05;
    //rel[i]*=1.01;
  }

  //calc the yields
  for (int i=0;i<nSets;i++){
    for (int j=0;j<nBins;j++){
      uPiLumi->Draw(Form("%s%d",yieldName.Data(),j),Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
      likeYield[i].push_back(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      
      uPiLumi->Draw(Form("%s%d",yieldName.Data(),j),Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
      unlikeYield[i].push_back(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    }
  }

  //calc the asyms:
  vector<double>asymcomp,asymerrcomp;
  
  for (int i=0;i<nSets;i++){
    for (int j=0;j<nBins;j++){
      double tasym, tasymerr;
      CalcAsymAndErr(&tasym,&tasymerr,
		     bpol[i],bpolErr[i],
		     ypol[i],ypolErr[i],
		     (unlikeYield[i])[j],sqrt((unlikeYield[i])[j]),
		     (likeYield[i])[j],sqrt((likeYield[i])[j]),
		     1,0,
		     rel[i],relErr[i]);
		     
      asym[i].push_back(tasym);
      asymErr[i].push_back(tasymerr);
      asym_byptbin[j].push_back(tasym);
      asymErr_byptbin[j].push_back(tasymerr);
      printf("%d,%d: asym=%1.2E + %1.2E\n",i,j,tasym,tasymerr);
      asymcomp.push_back(abs(tasym));
      asymerrcomp.push_back(tasymerr);
    }
  }

  TCanvas *c;
  TPad *cMinor;
  TGraphErrors *g;
  TLegend *leg;
  TH1F *h;
  if (0){//plot a check of the correlation between asymmetries and errors
    c=new TCanvas("cComp","cComp",800,400);
    TGraph *g0=new TGraph(asymcomp.size(),&(asymcomp[0]),&(asymerrcomp[0]));
    g0->SetTitle("Asymmetries vs errors;abs(asym);err");
    g0->Draw("A*");
  }
 
  
  
  
  c=new TCanvas("cAsym","cAsym",1600,800);
  c->Divide(1,2);
  cMinor=(TPad*)c->cd(1);
  cMinor->Divide(2,1);
  cMinor->cd(1);
  for (int i=0;i<nSets;i++){
    g=new TGraphErrors(nBins,&(ptmid[0]),&((asym[i])[0]),&(pterr[0]),&((asymErr[i])[0]));
    g->SetLineColor(i+1);
    g->SetMarkerColor(i+1);
    //g->SetTitle(cutName[i].Data());
    if (i==0){
      g->SetTitle("MPC asym by spin group;pT(GeV);asym");
      g->GetHistogram()->SetMaximum(0.15);
      g->GetHistogram()->SetMinimum(-0.15);
      g->Draw("AC*");
    } else{
          g->SetTitle(cutName[i].Data());
      g->Draw("C*");
    }
  }
  cMinor->cd(1)->SetGridy();
  leg=cMinor->cd(1)->BuildLegend();
  ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());

  
 //calc the bbc asyms:
  vector<double> setindex,setindexerr,bbcasym,bbcasymerr;
  for (int i=0;i<nSets;i++){
    uPiLumi->Draw("bbc",Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double bbcLike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
    uPiLumi->Draw("bbc",Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double bbcUnlike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
    uPiLumi->Draw("bbc_err*bbc_err",Form("(bspin==yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double bbcLike_err=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
    uPiLumi->Draw("bbc_err*bbc_err",Form("(bspin!=yspin && zdc>0)*(%s)",cut[i].Data()),"goff");
    double bbcUnlike_err=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
     double tasym, tasymerr;
    CalcAsymAndErr(&tasym,&tasymerr,
		   bpol[i],bpolErr[i],
		   ypol[i],ypolErr[i],
		   bbcUnlike,bbcUnlike_err,
		   bbcLike,bbcLike_err,
		   1,0,
		   rel[i],relErr[i]);
    setindex.push_back(i*1.0);
    setindexerr.push_back(0.1);
    bbcasym.push_back(tasym);
    bbcasymerr.push_back(tasymerr);
    printf("set %1.1f has bbc:%E,%E ==> asym=%E\n",setindex[i],bbcLike,bbcUnlike,tasym);

    }

  if (1){ //plot the BBC asymmetry by spin group:
    TH1F *hBbcAsym=new TH1F("hBbcAsym","BBC/ZDC Asym",40,-0.003,0.003);
    cMinor->cd(2);
    for (int i=0;i<nSets;i++){
      hBbcAsym->Fill(bbcasym[i]);
      g=new TGraphErrors(1,&(bbcasym[i]),&(setindex[i]),&(bbcasymerr[i]),&(bbcasymerr[i]));
      g->SetMarkerColor(i+1);
      if (i==0){
	g->SetTitle("BBC asym by spin group;group;asym");
	g->GetXaxis()->SetLimits(-0.003,0.003);
	g->GetHistogram()->SetMaximum(nSets*1.0-0.5);
	g->GetHistogram()->SetMinimum(-0.5);

      //g->GetXaxis()->SetLimits(-0.5,nSets*1.0-0.5);
      //g->GetHistogram()->SetMaximum(0.002);
      //g->GetHistogram()->SetMinimum(-0.002);
	g->Draw("AC*");
	//g->SetTitle(cutName[i].Data());
      } else {
	g->SetTitle(cutName[i].Data());
	g->Draw("C*");
      }
    }
    cMinor->cd(2)->SetGridy();
    leg=cMinor->cd(2)->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());
    hBbcAsym->Draw("same");
  }

  if (1){//plot the distributions in each pt bin:
    cMinor=(TPad*)c->cd(2);
    cMinor->Divide(5,2);
    for (int i=0;i<nBins;i++){
      cMinor->cd(i+1);
 
      g=new TGraphErrors(nSets,&(setindex[0]),&((asym_byptbin[i])[0]),&(setindexerr[0]),&((asymErr_byptbin[i])[0]));
      g->SetTitle(Form("Asym in bin %d",i));
      /*
      g->Draw("A*");
      for (int j=0;j<setindex.size();j++){
	printf("setindex=%1.1f, err=%1.1f, asym[%d][%d]=%1.1E, err=%1.1E\n",
	       setindex[j],setindexerr[j],i,j,asym_byptbin[i][j],asymErr_byptbin[i][j]);
	g=new TGraphErrors(1,&(setindex[j]),&((asym_byptbin[i])[j]),&(setindexerr[j]),&((asymErr_byptbin[i])[j]));
	    g->SetMarkerColor(j+1);
	g->Draw("*");
      }
      printf("test\n");
      */
      float minh=g->GetMinimum();
      float maxh=g->GetMaximum();
      float span=maxh-minh;
      float margin=0.25;
      //h=new TH1F(Form("h%d",i),Form("Asym in bin %d;asym",i),50,-0.05,0.05);//20,minh-span*0.25,minh+span*0.25);
      h=new TH1F(Form("h%d",i),Form("Asym in bin %d;asym",i),60,-0.1,0.1);//20,minh-span*0.25,minh+span*0.25);
      for (int j=0;j<setindex.size();j++){
	h->Fill(asym_byptbin[i][j]);
	if (fabs(asym_byptbin[i][j])>0.1){
	  printf("%d %d out of bounds: %s (bin %d)\n",i,j,cutName[j].Data(),i);
	}
      }
      h->Draw();
    }
  }
  

  if (0){//plot a check of the relative luminosities per selection
    c=new TCanvas("cRelLumi","cRelLumi",800,400);
    for (int i=0;i<nSets;i++){
      g=new TGraphErrors(1,&(setindex[i]),&(rel[i]),&(relErr[i]),&(relErr[i]));
      g->SetMarkerColor(i+1);
      if (i==0){
	g->SetTitle("Relative Lumi spin group;group;asym");
	g->GetXaxis()->SetLimits(-0.5,nSets*1.0-0.5);
	g->GetHistogram()->SetMaximum(1.5);
	g->GetHistogram()->SetMinimum(0.5);
	g->Draw("AC*");
	//g->SetTitle(cutName[i].Data());
      } else {
	g->SetTitle(cutName[i].Data());
	g->Draw("C*");
      }
    }
    c->cd(1)->SetGridy();
    leg=c->cd(1)->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());
  }

  return;
}




void   PlotAsymByRun(){
  //sums bbc scalers over each run and produces the per-run asym per group
  int nSets=4;
  vector<double> likeYield[4];
  vector<double> unlikeYield[4];
  vector<double> asym[4],asymErr[4];
  double rel[4],relErr[4];
  TString cut[4]={"(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"!(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)",
		"!(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)"};
  TString cutName[4]={"SSOO odd bunch",
		      "SSOO even bunch",
		      "OOSS odd bunch",
		      "OOSS even bunch"};

  //for each set we will have a runlist that may differ (definitely, since runs can't be more than one pattern)
  vector<double> runlist[nSets],runlisterr;
  //we will end up with a bbc asymmetry for each run for each group
    vector<double> bbc_asym[nSets],bbc_asym_err[nSets];
  //as well as an asymmetry for each pt bin for the same:
    vector<double> ptbin_asym[nSets][nBins],ptbin_asym_err[nSets][nBins];

  
  //get the runlist for each group:
  for (int i=0;i<nSets;i++){
    int ci=i;
    if (i==0) ci=1; //to force using the 'evens' definition.
    if (i==2) ci=3;
    uPiLumi->Draw("run",Form("(bunch==0 )&&(%s)",cut[ci].Data()),"goff"); 
    int nRuns=uPiLumi->GetSelectedRows();
    for (int r=0;r<nRuns;r++){
      runlist[i].push_back(uPiLumi->GetVal(0)[r]);
      printf("adding to set %d run %f (%1.2f)\n",i,uPiLumi->GetVal(0)[r],runlist[i][runlist[i].size()-1]);
      runlisterr.push_back(0.5);
    }
  }


  //for each set, iterate through its runlist:
  TString runcut;
  for (int s=0;s<nSets;s++){
    for(int r=0;r<runlist[s].size();r++){
      int run=(runlist[s])[r];

      //build the cut for this run.
      runcut=Form("(zdc>0)*(run==%d)*(%s)",run,cut[s].Data());
    
    //get the polarization for this run:
      double bpol, bpolErr, ypol,ypolErr;
      uPiLumi->Draw("bpol:bpol_err:ypol:ypol_err",Form("(%s)",runcut.Data()),"goff");
      if (uPiLumi->GetSelectedRows()==0) {
	printf("skipping run=%d : no data.\n",run);
	continue; //skip if there was no valid bunch for this run!
      }
      bpol=(uPiLumi->GetVal(0))[0];
      bpolErr=(uPiLumi->GetVal(1))[0];
      ypol=(uPiLumi->GetVal(2))[0];
      ypolErr=(uPiLumi->GetVal(3))[0];

      //assumption is that each bunch in the run has the same polarization, se we can read the first element of the return.

  
      //calc the rel lumi
      uPiLumi->Draw("zdc",Form("(bspin==yspin )*(%s)",runcut.Data()),"goff");
      int nLike=uPiLumi->GetSelectedRows();
      double zdcLike=SumArray(uPiLumi->GetVal(0),nLike);
      uPiLumi->Draw("zdc",Form("(bspin!=yspin)*(%s)",runcut.Data()),"goff");
      int nUnlike=uPiLumi->GetSelectedRows();
      double zdcUnlike=SumArray(uPiLumi->GetVal(0),nUnlike);
      //naively, the uncertainties in like an unlike are independent (though really they share the global scaler)
      uPiLumi->Draw("zdc_err*zdc_err",Form("(bspin==yspin)*(%s)",runcut.Data()),"goff");//square of error term
      double zdcLikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      uPiLumi->Draw("zdc_err*zdc_err",Form("(bspin!=yspin)*(%s)",runcut.Data()),"goff");
      double zdcUnlikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      double rel=(zdcLike/zdcUnlike);
      double relErr=rel*sqrt((zdcLikeErr/zdcLike)*(zdcLikeErr/zdcLike)
			    +(zdcUnlikeErr/zdcUnlike)*(zdcUnlikeErr/zdcUnlike));

      //get the yields and calc the asym per ptbin:
      for (int j=0;j<nBins;j++){
	uPiLumi->Draw(Form("%s%d",yieldName.Data(),j),Form("(bspin==yspin)*(%s)",runcut.Data()),"goff");
	double likeYield=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());	
	uPiLumi->Draw(Form("%s%d",yieldName.Data(),j),Form("(bspin!=yspin)*(%s)",runcut.Data()),"goff");
	double unlikeYield=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
	//calc the asyms:
	double tasym, tasymerr;
	if (likeYield+unlikeYield>0){
	  CalcAsymAndErr(&tasym,&tasymerr,
			 bpol,bpolErr,
			 ypol,ypolErr,
			 unlikeYield,sqrt(unlikeYield),
			 likeYield,sqrt(likeYield),
			 1,0,
			 rel,relErr);
	}
	if (likeYield+unlikeYield<1){
	  printf("run=%d,set=%d,bin=%d has yields %1.2E,%1.2E. Setting large errors.\n",run,s,j,likeYield,unlikeYield);
	  tasym=0;
	  tasymerr=1;
	}
	ptbin_asym[s][j].push_back(tasym);
	ptbin_asym_err[s][j].push_back(tasymerr);
	//printf("set=%d bin=%d: asym=%1.2E + %1.2E\n",s,j,tasym,tasymerr);
      }

      //get the bbc counts and calc the asym for the bbc:
      uPiLumi->Draw("bbc",Form("(bspin==yspin )*(%s)",runcut.Data()),"goff");
      double bbcLike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
      uPiLumi->Draw("bbc",Form("(bspin!=yspin)*(%s)",runcut.Data()),"goff");
      double bbcUnlike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
      uPiLumi->Draw("bbc_err*bbc_err",Form("(bspin==yspin )*(%s)",runcut.Data()),"goff");
      double bbcLike_err=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      uPiLumi->Draw("bbc_err*bbc_err",Form("(bspin!=yspin)*(%s)",runcut.Data()),"goff");
      double bbcUnlike_err=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      double tasym, tasymerr;
      CalcAsymAndErr(&tasym,&tasymerr,
		     bpol,bpolErr,
		     ypol,ypolErr,
		     bbcUnlike,bbcUnlike_err,
		     bbcLike,bbcLike_err,
		     1,0,
		     rel,relErr);
      bbc_asym[s].push_back(tasym);
      bbc_asym_err[s].push_back(tasymerr);
      printf("set %d has bbc:l=%E,u=%E ==> asym=%E\n",s,bbcLike,bbcUnlike,tasym);
     }//loop over runs
  }//loop over sets

  TCanvas *c;
  TGraphErrors *g;
  TLegend *leg;
 c=new TCanvas("cByRun","cByRun",1600,400);
  int nPtBinsToDraw=2;
  int ptBinToDraw[]={1,4};
  c->Divide(nPtBinsToDraw+1,1);

  //plot the MPC asyms for the selected bins, separated by group
  for (int i=0;i<nPtBinsToDraw;i++){
    c->cd(i+1);
    for (int s=0;s<nSets;s++){
      g=new TGraphErrors(runlist[s].size(),&((runlist[s])[0]),&((ptbin_asym[s][ptBinToDraw[i]])[0]),&(runlisterr[0]),&((ptbin_asym_err[s][ptBinToDraw[i]])[0]));
      g->SetLineColor(s+1);
      g->SetMarkerColor(s+1);
      //g->SetTitle(cutName[i].Data());
      if (s==0){
	g->SetTitle(Form("MPC bin%d asym by spin group;run;asym",ptBinToDraw[i]));
	g->GetHistogram()->SetMaximum(0.25);
	g->GetHistogram()->SetMinimum(-0.25);
	g->Draw("A*");
      } else{
	g->SetTitle(cutName[s].Data());
	g->Draw("*");
      }
    }
    c->cd(i+1)->SetGridy();
    leg=c->cd(i+1)->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());
  }
    //plot the BBC asyms for the selected bins, separated by group
  c->cd(nPtBinsToDraw+1);
  for (int s=0;s<nSets;s++){
    g=new TGraphErrors(runlist[s].size(),&((runlist[s])[0]),&((bbc_asym[s])[0]),&(runlisterr[0]),&((bbc_asym_err[s])[0]));

    //g=new TGraphErrors(1,&(setindex[i]),&(bbcasym[i]),&(bbcasymerr[i]),&(bbcasymerr[i]));
    g->SetMarkerColor(s+1);
    if (s==0){
      g->SetTitle("BBC asym by spin group;run;asym");
      g->GetHistogram()->SetMaximum(0.02);
      g->GetHistogram()->SetMinimum(-0.02);
      g->Draw("A*");
    } else {
      g->SetTitle(cutName[s].Data());
      g->Draw("*");
    }
  }
  c->cd(nPtBinsToDraw+1)->SetGridy();
  leg=c->cd(nPtBinsToDraw+1)->BuildLegend();
  ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());


  return;
  }



void   PlotAsymByFill(){
  //sums bbc scalers over each run and produces the per-run asym per group
  int nSets=4;
  vector<double> likeYield[nSets];
  vector<double> unlikeYield[nSets];
  vector<double> asym[nSets],asymErr[nSets];
  vector<double> rel_check[nSets],rel_err_check[nSets];
  double rel[nSets],relErr[nSets];
  TString cut[4]={"(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"!(bunch%2)&&(pat == 21 || pat == 24 || pat == 25 || pat == 28)",
		"(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)",
		"!(bunch%2)&&(pat == 22 || pat == 23 || pat == 26 || pat == 27)"};
  TString cutName[4]={"SSOO odd bunch",
		      "SSOO even bunch",
		      "OOSS odd bunch",
		      "OOSS even bunch"};

  //for each set we will have a runlist that may differ (definitely, since runs can't be more than one pattern)
  vector<double> filllist[nSets],filllisterr;
  //we will end up with a bbc asymmetry for each run for each group
    vector<double> bbc_asym[nSets],bbc_asym_err[nSets];
  //as well as an asymmetry for each pt bin for the same:
    vector<double> ptbin_asym[nSets][nBins],ptbin_asym_err[nSets][nBins];

  
  //get the runlist for each group:
  for (int i=0;i<nSets;i++){
    int ci=i;
    if (i==0) ci=1; //to force using the 'evens' definition.
    if (i==2) ci=3;
    uPiLumi->Draw("fill",Form("(bunch==0 )&&(%s)",cut[ci].Data()),"goff"); 
    int nFills=uPiLumi->GetSelectedRows();
    for (int r=0;r<nFills;r++){
      double newfill=uPiLumi->GetVal(0)[r];
      if (filllist[i].size()>0){
	if (filllist[i][filllist[i].size()-1]==newfill)
	  continue; //skip if we've already got this fill in our list.
      }
      filllist[i].push_back(uPiLumi->GetVal(0)[r]);
      printf("adding to set %d fill %f (%1.2f)\n",i,uPiLumi->GetVal(0)[r],filllist[i][filllist[i].size()-1]);
      filllisterr.push_back(0.5);
    }
  }


  //for each set, iterate through its filllist:
  TString fillcut;
  for (int s=0;s<nSets;s++){
    for(int r=0;r<filllist[s].size();r++){
      int fill=(filllist[s])[r];

      //build the cut for this fill.
      fillcut=Form("(zdc>0)*(fill==%d)*(%s)",fill,cut[s].Data());
    
    //get the polarization for this fill:
      double bpol, bpolErr, ypol,ypolErr;
      uPiLumi->Draw("bpol:bpol_err:ypol:ypol_err",Form("(%s)",fillcut.Data()),"goff");
      if (uPiLumi->GetSelectedRows()==0) {
	printf("skipping fill=%d : no data.\n",fill);
	continue; //skip if there was no valid bunch for this fill!
      }
      bpol=(uPiLumi->GetVal(0))[0];
      bpolErr=(uPiLumi->GetVal(1))[0];
      ypol=(uPiLumi->GetVal(2))[0];
      ypolErr=(uPiLumi->GetVal(3))[0];

      //assumption is that each bunch in the fill has the same polarization, se we can read the first element of the return.

  
      //calc the rel lumi
      uPiLumi->Draw("bbc:bbc_err*bbc_err",Form("(bspin==yspin )*(%s)",fillcut.Data()),"goff");
      int nLike=uPiLumi->GetSelectedRows();
      double zdcLike=SumArray(uPiLumi->GetVal(0),nLike);
      double zdcLikeErr=sqrt(SumArray(uPiLumi->GetVal(1),nLike));
      uPiLumi->Draw("bbc:bbc_err*bbc_err",Form("(bspin!=yspin )*(%s)",fillcut.Data()),"goff");
      int nUnlike=uPiLumi->GetSelectedRows();
      double zdcUnlike=SumArray(uPiLumi->GetVal(0),nUnlike);
      double zdcUnlikeErr=sqrt(SumArray(uPiLumi->GetVal(1),nUnlike));

      /*
      uPiLumi->Draw("zdc",Form("(bspin!=yspin)*(%s)",fillcut.Data()),"goff");
      int nUnlike=uPiLumi->GetSelectedRows();
      double zdcUnlike=SumArray(uPiLumi->GetVal(0),nUnlike);
      //naively, the uncertainties in like an unlike are independent (though really they share the global scaler)
      uPiLumi->Draw("zdc_err*zdc_err",Form("(bspin==yspin)*(%s)",fillcut.Data()),"goff");//square of error term
      double zdcLikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      uPiLumi->Draw("zdc_err*zdc_err",Form("(bspin!=yspin)*(%s)",fillcut.Data()),"goff");
      double zdcUnlikeErr=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      */
      double rel=(zdcLike/zdcUnlike);
      double relErr=rel*sqrt((zdcLikeErr/zdcLike)*(zdcLikeErr/zdcLike)
			    +(zdcUnlikeErr/zdcUnlike)*(zdcUnlikeErr/zdcUnlike));

      //get the yields and calc the asym per ptbin:
      for (int j=0;j<nBins;j++){
	uPiLumi->Draw(Form("%s%d",yieldName.Data(),j),Form("(bspin==yspin)*(%s)",fillcut.Data()),"goff");
	double likeYield=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());	
	uPiLumi->Draw(Form("%s%d",yieldName.Data(),j),Form("(bspin!=yspin)*(%s)",fillcut.Data()),"goff");
	double unlikeYield=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
	//calc the asyms:
	double tasym, tasymerr;
	if (likeYield+unlikeYield>0){
	  CalcAsymAndErr(&tasym,&tasymerr,
			 bpol,bpolErr,
			 ypol,ypolErr,
			 unlikeYield,sqrt(unlikeYield),
			 likeYield,sqrt(likeYield),
			 zdcLike,zdcLikeErr,
			 zdcUnlike,zdcUnlikeErr);
	  //			 1,0,
	  //			 rel,relErr);
	}
	if (likeYield+unlikeYield<1){
	  printf("fill=%d,set=%d,bin=%d has yields %1.2E,%1.2E. Setting large errors.\n",fill,s,j,likeYield,unlikeYield);
	  tasym=0;
	  tasymerr=1;
	}
	if (abs(tasym)>1.0){
	  printf("======> fill=%d,set=%d,bin=%d has asym=%1.2E. yields u=%1.2E (buns=%d), l=%1.2E (buns=%d), pols b=%1.2E, y=%1.2E.\n",
		 fill,s,j,tasym,likeYield,nLike,unlikeYield,nUnlike,bpol,ypol);
	}
	ptbin_asym[s][j].push_back(tasym);
	ptbin_asym_err[s][j].push_back(tasymerr);
	//printf("set=%d bin=%d: asym=%1.2E + %1.2E\n",s,j,tasym,tasymerr);
      }

      //get the bbc counts and calc the asym for the bbc:
      uPiLumi->Draw("bbc",Form("(bspin==yspin )*(%s)",fillcut.Data()),"goff");
      double bbcLike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
      uPiLumi->Draw("bbc",Form("(bspin!=yspin)*(%s)",fillcut.Data()),"goff");
      double bbcUnlike=SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows());
      uPiLumi->Draw("bbc_err*bbc_err",Form("(bspin==yspin )*(%s)",fillcut.Data()),"goff");
      double bbcLike_err=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      uPiLumi->Draw("bbc_err*bbc_err",Form("(bspin!=yspin)*(%s)",fillcut.Data()),"goff");
      double bbcUnlike_err=sqrt(SumArray(uPiLumi->GetVal(0),uPiLumi->GetSelectedRows()));
      double tasym, tasymerr;
      CalcAsymAndErr(&tasym,&tasymerr,
		     bpol,bpolErr,
		     ypol,ypolErr,
		     bbcUnlike,bbcUnlike_err,
		     bbcLike,bbcLike_err,
		     1,0,
		     rel,relErr);
      bbc_asym[s].push_back(tasym);
      bbc_asym_err[s].push_back(tasymerr);
      rel_check[s].push_back(rel);
      rel_err_check[s].push_back(relErr);
      printf("fill=%d,set=%d has bbc:l=%E,u=%E ==> asym=%E\n",fill,s,bbcLike,bbcUnlike,tasym);
     }//loop over fills
  }//loop over sets

  TCanvas *c;
  TGraphErrors *g;
  TLegend *leg;
 c=new TCanvas("cByFill","cByFill",1600,400);
  int nPtBinsToDraw=2;
  int ptBinToDraw[]={1,4};
  c->Divide(nPtBinsToDraw+1,1);

  //plot the relative lumi on its own, as a temporary curiosity:
  if (1){
    c->cd(1);
    for (int s=0;s<nSets;s++){
      g=new TGraphErrors(filllist[s].size(),&((filllist[s])[0]),&((rel_check[s])[0]),&(filllisterr[0]),&((rel_err_check[s])[0]));
      g->SetLineColor(s+1);
      g->SetMarkerColor(s+1);
      //g->SetTitle(cutName[i].Data());
      if (s==0){
	g->SetTitle("Relative Lumi (ZDC) vs fill;fill;Like/unlike");
	g->GetHistogram()->SetMaximum(0.25);
	g->GetHistogram()->SetMinimum(-0.25);
	g->Draw("A*");
      } else{
	g->SetTitle(cutName[s].Data());
	g->Draw("*");
      }
    }
    c->cd(1)->SetGridy();
    leg=c->cd(1)->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());
  }
  
  //plot the MPC asyms for the selected bins, separated by group
  for (int i=1;i<nPtBinsToDraw;i++){
    c->cd(i+1);
    for (int s=0;s<nSets;s++){
      g=new TGraphErrors(filllist[s].size(),&((filllist[s])[0]),&((ptbin_asym[s][ptBinToDraw[i]])[0]),&(filllisterr[0]),&((ptbin_asym_err[s][ptBinToDraw[i]])[0]));
      g->SetLineColor(s+1);
      g->SetMarkerColor(s+1);
      //g->SetTitle(cutName[i].Data());
      if (s==0){
	g->SetTitle(Form("MPC bin%d asym by spin group;fill;asym",ptBinToDraw[i]));
	g->GetHistogram()->SetMaximum(0.25);
	g->GetHistogram()->SetMinimum(-0.25);
	g->Draw("A*");
      } else{
	g->SetTitle(cutName[s].Data());
	g->Draw("*");
      }
    }
    c->cd(i+1)->SetGridy();
    leg=c->cd(i+1)->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());
  }
    //plot the BBC asyms for the selected bins, separated by group
  c->cd(nPtBinsToDraw+1);
  for (int s=0;s<nSets;s++){
    g=new TGraphErrors(filllist[s].size(),&((filllist[s])[0]),&((bbc_asym[s])[0]),&(filllisterr[0]),&((bbc_asym_err[s])[0]));

    //g=new TGraphErrors(1,&(setindex[i]),&(bbcasym[i]),&(bbcasymerr[i]),&(bbcasymerr[i]));
    g->SetMarkerColor(s+1);
    if (s==0){
      g->SetTitle("BBC asym by spin group;fill;asym");
      g->GetHistogram()->SetMaximum(0.02);
      g->GetHistogram()->SetMinimum(-0.02);
      g->Draw("A*");
    } else {
      g->SetTitle(cutName[s].Data());
      g->Draw("*");
    }
  }
  c->cd(nPtBinsToDraw+1)->SetGridy();
  leg=c->cd(nPtBinsToDraw+1)->BuildLegend();
  ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel(cutName[0].Data());


  return;
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
  if (0){
    //debugging info
    printf("asym=%1.2E\terror terms: bpol=%1.2E\t ypol=%1.2E\t like=%1.2E\t unlike=%1.2E\trel=%1.2E\n",
	 tasym,tbpolerrterm,typolerrterm,tlikesumerrterm,tunlikesumerrterm,trellumierrterm);
  }
  *asym=tasym;
  *asym_err=sqrt(err2);
  return;
}
