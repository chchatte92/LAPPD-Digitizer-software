#include <sys/wait.h>

void timebin(Int_t ncpu, std::string datafiles, Int_t evmax=0, std::string cond="1") {
  printf("Here");
  Int_t nch = 8; //32
  Int_t fCh= 24;   // first channel of group`	
  Int_t chg = 31;
  Int_t useGuardCh = 1; //1

  Int_t nch_mask[32] = { 1, 1, 1, 1, 1, 1, 1, 1,   //set back all to 1
                         1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1 };

   Double_t fixmax[32] = { -1,  -1, -1, -1, -1, -1, -1, -1,    // if neg check max by itself
                           -1,  -1, -1, -1, -1, -1, -1, -1,	// o/w put it here by hand
                           -1,  -1, -1, -1, -1, -1, -1, -1,
                           -1,  -1, -1, -1, -1, -1, -1, -1 };

  Int_t polarity[32] = { 1, 1, 1, 1, 1, 1, 1, 1,         // really polarity of the signal
                         1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1 };

  Double_t frmax[32] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,     //only leading edge at 50% if set say 0.1 and 0.9 fits 10-90% 
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

  Double_t frmin[32] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

  Double_t flatrange[32] = { 40, 40, 40, 40, 40, 40, 40, 40,     //baseline fit
                             40, 40, 40, 40, 40, 40, 40, 40,
                             40, 40, 40, 40, 40, 40, 40, 40,
                             40, 40, 40, 40, 40, 40, 40, 40 };

  Int_t maxfind = 900;				// upto 900 cells it will check max
  Int_t minfind[32] = {   5,   5,   5,   5,   5,   5,   5,   5,   //starting from 5/6th cell to be checked
                          5,   5,   5,   5,   5,   5,   5,   5, 
                          0,   0,   0,   0,   0,   0,   0,   0, 
                          0,   0,   0,   0,   0,   0,   0,   0 }; 
  
  Int_t doPlot = 0;
  if(ncpu==0) doPlot = 1;
///////////////////////

  //datafiles="ch_0_1/"+datafiles;
  //datafiles="/storage/gpfs_data/local/lhcb/users/vagnoni/data/"+datafiles;
	
  Int_t nchild=0;

  Int_t frk=0;
  for(Int_t i=0; i<ncpu; ++i) {
    frk=fork();
    if(frk<0) {
       cout << "Error: cannot fork" << endl;
       exit(1);
    }
    if(!frk) break;
    ++nchild;
  }

  if(frk) {
    int status;
    //printf("Are we here?\n");
    while(waitpid(-1,&status,0)>0);
	//printf("We are here\n");
    TChain chain("T","T");
    for(Int_t i=0; i<ncpu; ++i) {
      chain.Add(Form("buffer/tobemerged_%d.root/tt",i));
    }

    cout << "Childs done. Merging files" << endl;

    chain.Merge("merged.root");
    chain.Merge(Form("%s/merged.root",datafiles.c_str()));

    cout << "Done" << endl;

    exit(0);
  }

  if(ncpu==0) ncpu=1;

  TCanvas *c[32];

  struct {
    Double_t t[32];
    Double_t max[32];
    Double_t area[32];
    Double_t base[32];
    Int_t counter[32][51];
    Int_t n;
  } data;

  std::string branch="t[32]/D:max[32]/D:area[32]/D:base[32]/D:counter[32][51]/I:n/I";

  TTree *ttree = new TTree("tt","");
  ttree->Branch("data", &data.t[0], branch.c_str());

  TTree *ttreehelp = new TTree("thelp","");
  ttreehelp->Branch("data", &data.t[0], branch.c_str());

  std::string file[32];
  std::string fcal[32];
  std::string fcals[32];
  std::string ftimes[32];

  for(Int_t j=fCh+0; j<fCh+nch; j++) {
    if(!nch_mask[j]) continue;
    file[j] = Form("%s/wave_%d.dat",datafiles.c_str(),j);
    fcal[j] = Form("../calib_Bol/calib_cell_%d.txt",j);
    fcals[j] = Form("../calib_Bol/calib_sample_%d.txt",j);
    ftimes[j] = Form("../calib_Bol/calib_time_%d.txt",j); 
    //fcal[j] = Form("calib/calib_cell_%d.txt",j);
    //fcals[j] = Form("calib/calib_sample_%d.txt",j);
    //ftimes[j] = Form("calib/calib_time_%d.txt",j);
  }

  int fd[32];
  FILE *fdcal[32];
  FILE *fdcals[32];
  FILE *fdtimes[32];
  for(Int_t j=fCh+0; j<fCh+nch; j++) {
    //printf("Looop\n");
    if(!nch_mask[j]) continue;
    fd[j]=open(file[j].c_str(),O_RDONLY);
    if(fd[j]<0) nch_mask[j]=0;
    fdcal[j]=fopen(fcal[j].c_str(),"r");
    fdcals[j]=fopen(fcals[j].c_str(),"r");
    fdtimes[j]=fopen(ftimes[j].c_str(),"r");
    //if(fdcal[j] && fdcals[j] && fdtimes[j]) printf("ALL FILE FIUND");	
  }

  if(doPlot) for(Int_t i=fCh+0; i<fCh+nch; ++i) if(nch_mask[i]) c[i] = new TCanvas();

  TGraph *hg[32];

  for(Int_t i=fCh+0; i<fCh+nch; ++i) {
    if(!nch_mask[i]) continue;
    hg[i] = new TGraph(1024);
    hg[i]->SetName(Form("hg%d",i));
    hg[i]->SetMarkerStyle(20);
    hg[i]->SetMarkerSize(0.5); 
  }

  Double_t p0[32][1024];
  Double_t p1[32][1024];
  Double_t p2[32][1024];

  Double_t pa0[32][1024];
  Double_t pa1[32][1024];
  Double_t pa2[32][1024];

  Double_t timev[32][1024];
  for(Int_t j=fCh+0; j<fCh+nch; ++j) {
    if(!nch_mask[j]) continue;
    for(Int_t i=0; i<1024; ++i) {
      fscanf(fdcal[j],"%lg %lg %lg", &p0[j][i], &p1[j][i], &p2[j][i]);
      fscanf(fdcals[j],"%lg %lg %lg", &pa0[j][i], &pa1[j][i], &pa2[j][i]);
      fscanf(fdtimes[j],"%lg", &timev[j][i]);
    }
  }

  Int_t n=0;
  while(++n) {
    if(evmax && n>evmax){printf("persing out\n"); goto out;};
    uint16_t startcell[32];
    Float_t fbuf[32][1024];
    for(Int_t j=fCh+0; j<fCh+nch; ++j) {
      if(!nch_mask[j]) continue;
      ULong64_t offset = (ULong64_t)((nchild+(n-1)*ncpu))*(sizeof(uint16_t)+1024*sizeof(Float_t));
      ULong64_t lsk=lseek64(fd[j], offset, SEEK_SET);
      if(lsk!=offset) {printf("persing out");goto out;}
      if(!read(fd[j],&startcell[j],sizeof(uint16_t))) {printf("persing out");goto out;}
      read(fd[j],&fbuf[j][0],sizeof(Float_t)*1024);

      //if(nchild==0) cout << lsk << " " << offset << " " << startcell << endl;

      Double_t tinc=0;
      for(Int_t i=0; i<1024; ++i) {
        Int_t index = (i+startcell[j])%1024;
        Double_t v = p0[j][index]+p1[j][index]*(fbuf[j][i]-p2[j][index]);
        v = pa0[j][i]+pa1[j][i]*(v-pa2[j][i]);
        hg[j]->SetPoint(i,tinc,v*polarity[j]);
//        if(j==30) cout << n << " " << j << " " << i << " " << index << " " << tinc << " " << v << endl;
        tinc += timev[j][index];
      }
    }

    for(Int_t j=fCh+0; j<fCh+nch; ++j) {
      if(!nch_mask[j]) continue;
      Double_t t0;
      Int_t index = (1024-startcell[j])%1024;
      if(j<8) t0 = hg[0]->GetPointX(index);
      else if(j<16) t0 = hg[8]->GetPointX(index);
      else if(j<24) t0 = hg[16]->GetPointX(index);
      else t0 = hg[24]->GetPointX(index);
      //else t0 = hg[25]->GetPointX(index);

      Double_t tch0 = hg[j]->GetPointX(index);
      for(Int_t i=0; i<1024; ++i) {
        hg[j]->SetPointX(i,hg[j]->GetPointX(i)-(tch0-t0));
      }
    }

    if(useGuardCh) {

      for(Int_t j=fCh+0; j<fCh+nch; ++j) {
        if(j==chg) continue;
        if(!nch_mask[j]) continue;

        Int_t probs[1024];
        memset(probs,0,1024*sizeof(Int_t));
        hg[chg]->Fit("pol0", "0Q", "", 5, 970);
        Double_t zerog=hg[chg]->GetFunction("pol0")->GetParameter(0);
        for(Int_t i=0; i<1024; ++i) {
          if(fabs(hg[chg]->GetPointY(i)-zerog)>0.005) probs[i]=1;   
          //cout << j << " " << i << " " << probs[i] << " " << hg[chg]->GetPointY(i) << " " << zerog << endl;
        }

        Int_t pr[1024];
        memcpy(pr,probs,1024*sizeof(Int_t));

        if(pr[0]) {
          Int_t m;
          for(m=1; m<1024; ++m) if(!pr[m]) break;
          if(m==1024) { cout << "Impossible... exiting" << endl; exit(1); }
          hg[j]->SetPointY(0,hg[j]->GetPointY(m));
          pr[0]=0;
        }

        if(pr[1023]) {
          Int_t m;
          for(m=1022; m>=0; --m) if(!pr[m]) break;
          if(m==-1) { cout << "Impossible... exiting" << endl; exit(1); }
          hg[j]->SetPointY(1023,hg[j]->GetPointY(m));
          pr[1023]=0;
        }

        for(Int_t i=1; i<1023; ++i) {
          if(pr[i]) {
            Int_t mu, md;
            for(mu=i+1; mu<1024; ++mu) if(!pr[mu]) break;
            for(md=i-1; md<1024; --md) if(!pr[md]) break;
            Double_t diff = hg[j]->GetPointY(mu)-hg[j]->GetPointY(md);

            Double_t diffX = hg[j]->GetPointX(mu)-hg[j]->GetPointX(md);
            Double_t diffY = hg[j]->GetPointY(mu)-hg[j]->GetPointY(md);

            hg[j]->SetPointY(i,diffY/diffX*(hg[j]->GetPointX(i)-hg[j]->GetPointX(md))+hg[j]->GetPointY(md));
          //cout << "point " << i << " channel " << j << " set to " << hg[j]->GetPointY(i) << endl;
            pr[i]=0;
          }
        }
      }
    }


    Double_t max[32];
    Double_t xMin[32];
    Double_t xMax[32];
    Double_t zero[32];

    Double_t area[32];
    memset(area,0,32*sizeof(Double_t));

    for(Int_t j=fCh+0; j<fCh+nch; ++j) {
      if(!nch_mask[j]) continue;
      hg[j]->Fit("pol0", "0Q", "", 5, flatrange[j]);
      TF1 *fit = hg[j]->GetFunction("pol0");
      zero[j] = fit->GetParameter(0);

      Int_t i;
      //max[j] = TMath::MaxElement(maxfind,&(hg[j]->GetY()[1]))-zero[j];
      max[j] = TMath::MaxElement(maxfind-minfind[j],&(hg[j]->GetY()[minfind[j]+1]))-zero[j];
      if(fixmax[j]>0) max[j]=fixmax[j];

      //for(xMax[j] = 0, i = 1; i < 1024; ++i)
      for(xMax[j] = 0, i = minfind[j]+1; i < maxfind; ++i)
        if(hg[j]->GetPointY(i)-zero[j] > frmax[j]*max[j]) {
          xMax[j] = hg[j]->GetPointX(i);
          break;
        }
      Int_t imem = i;
      Int_t izero = 0;
      for(; i < 1023; ++i)
        if(hg[j]->GetPointY(i)-zero[j] < 0) {
          izero = i;
          break;
        }
      i = imem;
      for(xMin[j] = 200; i > 0; --i)
        if(hg[j]->GetPointY(i)-zero[j] < frmin[j]*max[j]) {
          xMin[j] = hg[j]->GetPointX(i);
          break;
        }
      //for(Int_t i=5; i<maxfind; ++i) {
      for(Int_t i=minfind[j]; i<izero; ++i) {
        //area[j]+=(hg[j]->GetPointY(i)-zero[j])*(hg[j]->GetPointX(i+1)-hg[j]->GetPointX(i));
         area[j]+=0.5*(hg[j]->GetPointY(i)+hg[j]->GetPointY(i+1)-2*zero[j])*(hg[j]->GetPointX(i+1)-hg[j]->GetPointX(i));
      }

    }
    Double_t t[32];

    for(Int_t j=fCh+0; j<fCh+nch; ++j) {
      if(!nch_mask[j]) continue;

      hg[j]->Fit("pol1", "Q", "", xMin[j], xMax[j]);
      TF1 *fit = hg[j]->GetFunction("pol1");
      Double_t p0 = fit->GetParameter(0);
      Double_t p1 = fit->GetParameter(1);
      t[j] = (0.5*(frmax[j]+frmin[j])*max[j]+zero[j]-p0)/p1;
      //printf("CHECK %d %g %g %g %g %g\n",j,xMin[j],xMax[j],t[j],max[j]+zero[j],zero[j]);

    }

    if(!(n%1000)) cout << "Child " << nchild << " Event " << n << endl;

    data.n=nchild+(n-1)*ncpu;
    data.max[0]=max[0];

    memset(&data.counter[0][0], 0, sizeof(data.counter));

    for(Int_t k=fCh+0; k<fCh+nch; ++k) {
      if(!nch_mask[k]) continue;
      data.t[k] = t[k];
      data.max[k] = max[k];
      data.area[k] = area[k];
      data.base[k] = zero[k];
      for(Int_t jj = 0; jj<=50; ++jj) {
        for(Int_t kk = minfind[k]; kk<maxfind; ++kk) {
          if(hg[k]->GetPointY(kk)-zero[k] > jj*1e-3) {
            ++(data.counter[k][jj]);
            kk+=19;
          }
        }
      }
    }    

    if(doPlot) {

      ttreehelp->Fill();
      if(ttreehelp->GetEntries(Form("%s%s",Form("Entry$==%d&&",n-1),cond.c_str()))==0) continue;

      //cout << data.t[1] - data.t[5] << endl;
      //cout << data.t[5] << endl;

      //if(n!=30228) continue;
      //if(n!=30 && n!=2108 && n!=9809 && n!=11299 && n!=11676 && n!=12978  && n!=16806  && n!=19733 && n!=19767 && n!=25557) continue;

      Int_t lastk;
       	
      for(Int_t k=fCh+0; k<fCh+nch; ++k) {
        if(!nch_mask[k]) continue;
        lastk=k;
        c[k]->cd();
	//for(int cBol=0; cBol<1024; cBol++){
	//	printf("%lf\n",hg[k]->GetPointY(cBol));
	//}
        hg[k]->Draw("AP");
        c[k]->Update();
        cout << k << " " << zero[k] << " " << data.counter[k][5] << endl;
      }
      //for(int kkk=0; kkk<1024; ++kkk) cout << hg[0]->GetPointX(kkk) << " " << hg[0]->GetPointY(kkk) << endl;
      //return;
      c[lastk]->WaitPrimitive();
    }
    //printf("Filling tree\n");	
    ttree->Fill();
    //printf("Tree Filled\n");	
  }

  out:

  TFile *outfile=TFile::Open(Form("buffer/tobemerged_%d.root",nchild),"RECREATE");
  outfile->WriteTObject(ttree,"","Overwrite");
  outfile->Close(); 

  exit(0);

}

