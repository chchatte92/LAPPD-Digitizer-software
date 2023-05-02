void calib_3points() {

  // Total number of channels
  //Int_t nch = 19;
  Int_t nch = 1;
 	
  // Should differ from zero if one does not want to read the full calibration sample (e.g. for small stat tests)
  Int_t evmax_cell=0;
  Int_t evmax_sample=0;
  Int_t evmax_tc=0;
  Int_t evmax_tgc=20000;

  // DC voltage points taken for voltage calibration (9 points for each channel, but could be done with less)
  Double_t v[19][3] = { {-0.400,0,0.400},   //0.3
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300},
                        {-0.300,0,0.300} };

  std::string file[3][19];
  std::string filet[19];
  std::string fileg[19];
  std::string filetgl[19];
  std::string fileglg[19];

  // Calibration data for each channel (3 files x 16 channels)
  file[0][0] = "data/2021_11_05/calibMarco_amp_ch0_m300mV/wave_0.dat";
  file[1][0] = "data/2021_11_05/calibMarco_amp_ch0_0mV/wave_0.dat";
  file[2][0] = "data/2021_11_05/calibMarco_amp_ch0_300mV/wave_0.dat";

  file[0][1] = "data/2021_11_05/calibMarco_amp_ch1_m300mV/wave_1.dat";
  file[1][1] = "data/2021_11_05/calibMarco_amp_ch1_0mV/wave_1.dat";
  file[2][1] = "data/2021_11_05/calibMarco_amp_ch1_300mV/wave_1.dat";

  file[0][2] = "data/2021_11_05/calibMarco_amp_ch2_m300mV/wave_2.dat";
  file[1][2] = "data/2021_11_05/calibMarco_amp_ch2_0mV/wave_2.dat";
  file[2][2] = "data/2021_11_05/calibMarco_amp_ch2_300mV/wave_2.dat";

  file[0][3] = "data/2021_11_05/calibMarco_amp_ch3_m300mV/wave_3.dat";
  file[1][3] = "data/2021_11_05/calibMarco_amp_ch3_0mV/wave_3.dat";
  file[2][3] = "data/2021_11_05/calibMarco_amp_ch3_300mV/wave_3.dat";

  file[0][4] = "data/2021_11_05/calibMarco_amp_ch4_m300mV/wave_4.dat";
  file[1][4] = "data/2021_11_05/calibMarco_amp_ch4_0mV/wave_4.dat";
  file[2][4] = "data/2021_11_05/calibMarco_amp_ch4_300mV/wave_4.dat";

  file[0][5] = "data/2021_11_05/calibMarco_amp_ch5_m300mV/wave_5.dat";
  file[1][5] = "data/2021_11_05/calibMarco_amp_ch5_0mV/wave_5.dat";
  file[2][5] = "data/2021_11_05/calibMarco_amp_ch5_300mV/wave_5.dat";

  file[0][6] = "data/2021_11_05/calibMarco_amp_ch6_m300mV/wave_6.dat";
  file[1][6] = "data/2021_11_05/calibMarco_amp_ch6_0mV/wave_6.dat";
  file[2][6] = "data/2021_11_05/calibMarco_amp_ch6_300mV/wave_6.dat";

  file[0][7] = "data/2021_11_05/calibMarco_amp_ch7_m300mV/wave_7.dat";
  file[1][7] = "data/2021_11_05/calibMarco_amp_ch7_0mV/wave_7.dat";
  file[2][7] = "data/2021_11_05/calibMarco_amp_ch7_300mV/wave_7.dat";

  file[0][8] = "data/2021_11_05/calibMarco_amp_ch8_m300mV/wave_8.dat";
  file[1][8] = "data/2021_11_05/calibMarco_amp_ch8_0mV/wave_8.dat";
  file[2][8] = "data/2021_11_05/calibMarco_amp_ch8_300mV/wave_8.dat";

  file[0][9] = "data/2021_11_05/calibMarco_amp_ch9_m300mV/wave_9.dat";
  file[1][9] = "data/2021_11_05/calibMarco_amp_ch9_0mV/wave_9.dat";
  file[2][9] = "data/2021_11_05/calibMarco_amp_ch9_300mV/wave_9.dat";

  file[0][10] = "data/2021_11_05/calibMarco_amp_ch10_m300mV/wave_10.dat";
  file[1][10] = "data/2021_11_05/calibMarco_amp_ch10_0mV/wave_10.dat";
  file[2][10] = "data/2021_11_05/calibMarco_amp_ch10_300mV/wave_10.dat";

  file[0][11] = "data/2021_11_05/calibMarco_amp_ch11_m300mV/wave_11.dat";
  file[1][11] = "data/2021_11_05/calibMarco_amp_ch11_0mV/wave_11.dat";
  file[2][11] = "data/2021_11_05/calibMarco_amp_ch11_300mV/wave_11.dat";

  file[0][12] = "data/2021_11_05/calibMarco_amp_ch12_m300mV/wave_12.dat";
  file[1][12] = "data/2021_11_05/calibMarco_amp_ch12_0mV/wave_12.dat";
  file[2][12] = "data/2021_11_05/calibMarco_amp_ch12_300mV/wave_12.dat";

  file[0][13] = "data/2021_11_05/calibMarco_amp_ch13_m300mV/wave_13.dat";
  file[1][13] = "data/2021_11_05/calibMarco_amp_ch13_0mV/wave_13.dat";
  file[2][13] = "data/2021_11_05/calibMarco_amp_ch13_300mV/wave_13.dat";

  file[0][14] = "data/2021_11_05/calibMarco_amp_ch14_m300mV/wave_14.dat";
  file[1][14] = "data/2021_11_05/calibMarco_amp_ch14_0mV/wave_14.dat";
  file[2][14] = "data/2021_11_05/calibMarco_amp_ch14_300mV/wave_14.dat";

  file[0][15] = "data/2021_11_05/calibMarco_amp_ch15_m300mV/wave_15.dat";
  file[1][15] = "data/2021_11_05/calibMarco_amp_ch15_0mV/wave_15.dat";
  file[2][15] = "data/2021_11_05/calibMarco_amp_ch15_300mV/wave_15.dat";

  file[0][16] = "data/2021_11_05/calibMarco_amp_ch16_m300mV/wave_16.dat";
  file[1][16] = "data/2021_11_05/calibMarco_amp_ch16_0mV/wave_16.dat";
  file[2][16] = "data/2021_11_05/calibMarco_amp_ch16_300mV/wave_16.dat";

  file[0][17] = "data/2021_11_05/calibMarco_amp_ch17_m300mV/wave_17.dat";
  file[1][17] = "data/2021_11_05/calibMarco_amp_ch17_0mV/wave_17.dat";
  file[2][17] = "data/2021_11_05/calibMarco_amp_ch17_300mV/wave_17.dat";

  file[0][18] = "data/2021_11_05/calibMarco_amp_ch23_m300mV/wave_23.dat";
  file[1][18] = "data/2021_11_05/calibMarco_amp_ch23_0mV/wave_23.dat";
  file[2][18] = "data/2021_11_05/calibMarco_amp_ch23_300mV/wave_23.dat";

  // Sine waveforms for local time calibration (50 MHz waveform, were actually acquired at 1 kHz)
  filet[0] = "data/2021_11_05/calibMarco_0_sine_50MHz/wave_0.dat";
  filet[1] = "data/2021_11_05/calibMarco_1_sine_50MHz/wave_1.dat";
  filet[2] = "data/2021_11_05/calibMarco_2_sine_50MHz/wave_2.dat";
  filet[3] = "data/2021_11_05/calibMarco_3_sine_50MHz/wave_3.dat";
  filet[4] = "data/2021_11_05/calibMarco_4_sine_50MHz/wave_4.dat";
  filet[5] = "data/2021_11_05/calibMarco_5_sine_50MHz/wave_5.dat";
  filet[6] = "data/2021_11_05/calibMarco_6_sine_50MHz/wave_6.dat";
  filet[7] = "data/2021_11_05/calibMarco_7_sine_50MHz/wave_7.dat";
  filet[8] = "data/2021_11_05/calibMarco_8_sine_50MHz/wave_8.dat";
  filet[9] = "data/2021_11_05/calibMarco_9_sine_50MHz/wave_9.dat";
  filet[10] = "data/2021_11_05/calibMarco_10_sine_50MHz/wave_10.dat";
  filet[11] = "data/2021_11_05/calibMarco_11_sine_50MHz/wave_11.dat";
  filet[12] = "data/2021_11_05/calibMarco_12_sine_50MHz/wave_12.dat";
  filet[13] = "data/2021_11_05/calibMarco_13_sine_50MHz/wave_13.dat";
  filet[14] = "data/2021_11_05/calibMarco_14_sine_50MHz/wave_14.dat";
  filet[15] = "data/2021_11_05/calibMarco_15_sine_50MHz/wave_15.dat";
  filet[16] = "data/2021_11_05/calibMarco_16_sine_50MHz/wave_16.dat";
  filet[17] = "data/2021_11_05/calibMarco_17_sine_50MHz/wave_17.dat";
  filet[18] = "data/2021_11_05/calibMarco_23_sine_50MHz/wave_23.dat";

  // Files of guard channels for local time calibration. This is for correcting random spikes.
  // The idea is that we acquire an open channel (which should be constant)
  // along with the signal channel, and look for spikes there
  fileg[0] = "data/2021_11_05/calibMarco_0_sine_50MHz/wave_7.dat";
  fileg[1] = "data/2021_11_05/calibMarco_1_sine_50MHz/wave_7.dat";
  fileg[2] = "data/2021_11_05/calibMarco_2_sine_50MHz/wave_7.dat";
  fileg[3] = "data/2021_11_05/calibMarco_3_sine_50MHz/wave_7.dat";
  fileg[4] = "data/2021_11_05/calibMarco_4_sine_50MHz/wave_7.dat";
  fileg[5] = "data/2021_11_05/calibMarco_5_sine_50MHz/wave_7.dat";
  fileg[6] = "data/2021_11_05/calibMarco_6_sine_50MHz/wave_7.dat";
  fileg[7] = "data/2021_11_05/calibMarco_7_sine_50MHz/wave_0.dat";
  fileg[8] = "data/2021_11_05/calibMarco_8_sine_50MHz/wave_15.dat";
  fileg[9] = "data/2021_11_05/calibMarco_9_sine_50MHz/wave_15.dat";
  fileg[10] = "data/2021_11_05/calibMarco_10_sine_50MHz/wave_15.dat";
  fileg[11] = "data/2021_11_05/calibMarco_11_sine_50MHz/wave_15.dat";
  fileg[12] = "data/2021_11_05/calibMarco_12_sine_50MHz/wave_15.dat";
  fileg[13] = "data/2021_11_05/calibMarco_13_sine_50MHz/wave_15.dat";
  fileg[14] = "data/2021_11_05/calibMarco_14_sine_50MHz/wave_15.dat";
  fileg[15] = "data/2021_11_05/calibMarco_15_sine_50MHz/wave_8.dat";
  fileg[16] = "data/2021_11_05/calibMarco_16_sine_50MHz/wave_23.dat";
  fileg[17] = "data/2021_11_05/calibMarco_17_sine_50MHz/wave_23.dat";
  fileg[18] = "data/2021_11_05/calibMarco_23_sine_50MHz/wave_16.dat";

  // Calibration files for global time calibration (sine wave at 100 MHz, were actually acquired at 1000 Hz)
  Double_t frequency=100;
  filetgl[0] = "data/2021_11_05/calibMarco_0_sine_100MHz/wave_0.dat";
  filetgl[1] = "data/2021_11_05/calibMarco_1_sine_100MHz/wave_1.dat";
  filetgl[2] = "data/2021_11_05/calibMarco_2_sine_100MHz/wave_2.dat";
  filetgl[3] = "data/2021_11_05/calibMarco_3_sine_100MHz/wave_3.dat";
  filetgl[4] = "data/2021_11_05/calibMarco_4_sine_100MHz/wave_4.dat";
  filetgl[5] = "data/2021_11_05/calibMarco_5_sine_100MHz/wave_5.dat";
  filetgl[6] = "data/2021_11_05/calibMarco_6_sine_100MHz/wave_6.dat";
  filetgl[7] = "data/2021_11_05/calibMarco_7_sine_100MHz/wave_7.dat";
  filetgl[8] = "data/2021_11_05/calibMarco_8_sine_100MHz/wave_8.dat";
  filetgl[9] = "data/2021_11_05/calibMarco_9_sine_100MHz/wave_9.dat";
  filetgl[10] = "data/2021_11_05/calibMarco_10_sine_100MHz/wave_10.dat";
  filetgl[11] = "data/2021_11_05/calibMarco_11_sine_100MHz/wave_11.dat";
  filetgl[12] = "data/2021_11_05/calibMarco_12_sine_100MHz/wave_12.dat";
  filetgl[13] = "data/2021_11_05/calibMarco_13_sine_100MHz/wave_13.dat";
  filetgl[14] = "data/2021_11_05/calibMarco_14_sine_100MHz/wave_14.dat";
  filetgl[15] = "data/2021_11_05/calibMarco_15_sine_100MHz/wave_15.dat";
  filetgl[16] = "data/2021_11_05/calibMarco_16_sine_100MHz/wave_16.dat";
  filetgl[17] = "data/2021_11_05/calibMarco_17_sine_100MHz/wave_17.dat";
  filetgl[18] = "data/2021_11_05/calibMarco_23_sine_100MHz/wave_23.dat";

  // Files of guard channels for global time calibration
  fileglg[0] = "data/2021_11_05/calibMarco_0_sine_100MHz/wave_7.dat";
  fileglg[1] = "data/2021_11_05/calibMarco_1_sine_100MHz/wave_7.dat";
  fileglg[2] = "data/2021_11_05/calibMarco_2_sine_100MHz/wave_7.dat";
  fileglg[3] = "data/2021_11_05/calibMarco_3_sine_100MHz/wave_7.dat";
  fileglg[4] = "data/2021_11_05/calibMarco_4_sine_100MHz/wave_7.dat";
  fileglg[5] = "data/2021_11_05/calibMarco_5_sine_100MHz/wave_7.dat";
  fileglg[6] = "data/2021_11_05/calibMarco_6_sine_100MHz/wave_7.dat";
  fileglg[7] = "data/2021_11_05/calibMarco_7_sine_100MHz/wave_0.dat";
  fileglg[8] = "data/2021_11_05/calibMarco_8_sine_100MHz/wave_15.dat";
  fileglg[9] = "data/2021_11_05/calibMarco_9_sine_100MHz/wave_15.dat";
  fileglg[10] = "data/2021_11_05/calibMarco_10_sine_100MHz/wave_15.dat";
  fileglg[11] = "data/2021_11_05/calibMarco_11_sine_100MHz/wave_15.dat";
  fileglg[12] = "data/2021_11_05/calibMarco_12_sine_100MHz/wave_15.dat";
  fileglg[13] = "data/2021_11_05/calibMarco_13_sine_100MHz/wave_15.dat";
  fileglg[14] = "data/2021_11_05/calibMarco_14_sine_100MHz/wave_15.dat";
  fileglg[15] = "data/2021_11_05/calibMarco_15_sine_100MHz/wave_8.dat";
  fileglg[16] = "data/2021_11_05/calibMarco_16_sine_100MHz/wave_23.dat";
  fileglg[17] = "data/2021_11_05/calibMarco_17_sine_100MHz/wave_23.dat";
  fileglg[18] = "data/2021_11_05/calibMarco_23_sine_100MHz/wave_16.dat";

  // Corresponding guard channel for each signal channel. We use always the last one of a group of 8 channels
  // except when calibrating the last channel, in which case we use one of the others
  Int_t guardch[19] = { 7,  7,  7,  7,  7,  7,  7,  0, 
                        15, 15, 15, 15, 15, 15, 15, 8,
                        18, 18, 16 };

  Int_t tbins = 1024;
  Double_t tmin = 0-0.5;
  Double_t tmax = 1023+0.5;

  // Open all input files. No check for failures... to be fixed
  int fd[3][19];
  TH1D *h[3][19];
  TH1D *ha[3][19];
  for(Int_t j=0; j<nch; ++j) {
    for(Int_t k=0; k<3; ++k) {
      fd[k][j]=open(file[k][j].c_str(),O_RDONLY);
      if(fd[k][j]<0) cout << "Missing file " << file[k][j].c_str() << " for CH " << j << endl;
      h[k][j] = new TH1D(Form("h%d%d",k,j),Form("h%d%d",k,j),tbins,tmin,tmax);
      ha[k][j] = new TH1D(Form("ha%d%d",k,j),Form("ha%d%d",k,j),tbins,tmin,tmax);
    }
  }

  int fdt[19];
  int fdg[19];
  int fdtgl[19];
  int fdglg[19];
  for(Int_t j=0; j<nch; ++j) {
    fdt[j]=open(filet[j].c_str(),O_RDONLY);
    if(fdt[j]<0) cout << "Missing file " << filet[j].c_str() << " for CH " << j << endl;
    fdg[j]=open(fileg[j].c_str(),O_RDONLY);
    if(fdg[j]<0) cout << "Missing file " << fileg[j].c_str() << " for CH " << j << endl;
    fdtgl[j]=open(filetgl[j].c_str(),O_RDONLY);
    if(fdtgl[j]<0) cout << "Missing file " << filetgl[j].c_str() << " for CH " << j << endl;
    fdglg[j]=open(fileglg[j].c_str(),O_RDONLY);
    if(fdglg[j]<0) cout << "Missing file " << fileglg[j].c_str() << " for CH " << j << endl;
  }


  TH1D *ht[19];
  TH1D *hg[19];

  for(Int_t j=0; j<nch; ++j) {
    ht[j] = new TH1D(Form("ht%d",j),Form("ht%d",j),tbins,tmin,tmax);
    hg[j] = new TH1D(Form("hg%d",j),Form("hg%d",j),tbins,tmin,tmax);
  }

  TGraph *gr=new TGraph(3);
  gr->SetName("gr");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);

  // This is the first estimate of one readout period (quite accurate BTW). In any case it will be improved
  // by the global time calibration, which works with real time differences
  Double_t timerange=204.8;
  //Double_t timerange=409.6;

  TF1 *func = new TF1("func","[0]+[1]*(x-[2])",-5000,5000);

  cout << "Cell calibration" << endl;

  Int_t n[3][19][1024];
  memset(n, 0, 3*19*1024*sizeof(Int_t));
  Int_t flag=1;
  Int_t evnt=0;
  while(flag) {
    ++evnt;
    if(evmax_cell && evnt==evmax_cell+1) break;
    if(!(evnt%1000)) cout << "Event " << evnt << endl;
    flag=0;
    uint16_t startcell[3];
    Float_t fbuf[3][1024];
    for(Int_t j=0; j<nch; ++j) {
      for(Int_t k=0; k<3; ++k) {
        // In our record, the first is the starting cells, and then the list of 1024 ADC counts
        // It's a float just because CAEN wavedump stores floats (as they could be online corrected, although
        // we take raw data which are uncorrected
        if(!read(fd[k][j],&startcell[k],sizeof(uint16_t))) continue;
        read(fd[k][j],&fbuf[k][0],sizeof(Float_t)*1024);
        for(Int_t i=0; i<1024; ++i) {
           // Discard the first sample and about the last 50 samples, don't want to waste time
           // with problematic regions that we are not going to use anyway
           if(i==0||i>970) continue;
           Int_t index = (i+startcell[k])%1024;
           ++n[k][j][index];
           h[k][j]->SetBinContent(index+1,(Double_t)fbuf[k][i]+h[k][j]->GetBinContent(index+1));
        }
        ++flag;
      }
    }
  }

  for(Int_t j=0; j<nch; ++j) {
    for(Int_t k=0; k<3; ++k) {
      for(Int_t i=0; i<1024; ++i) {
        h[k][j]->SetBinContent(i+1,h[k][j]->GetBinContent(i+1)/n[k][j][i]);
      }
    }
  }

  // Now perform fits for all cells of each channel, and save voltage calibration constants

  Double_t p0[19][1024];
  Double_t p1[19][1024];
  Double_t p2[19][1024];

  for(Int_t j=0; j<nch; ++j) {
    TString filename=Form("calib_Marco/calib_cell_%d.txt",j);
    FILE *file=fopen(filename.Data(),"w");

    for(Int_t i=0; i<1024; ++i) {
      Double_t a[3];
      for(Int_t k=0; k<3; ++k) {
        a[k]=h[k][j]->GetBinContent(i+1);
        gr->SetPoint(k,a[k],v[j][k]);
      }
      func->FixParameter(2,a[1]);

      gr->Fit("func","0Q","");
      p0[j][i] = func->GetParameter(0);
      p1[j][i] = func->GetParameter(1);
      p2[j][i] = func->GetParameter(2);

      fprintf(file,"%g %g %g\n", p0[j][i], p1[j][i], p2[j][i]); 

    }
    fclose(file);
  }


  for(Int_t j=0; j<nch; ++j) {
    for(Int_t k=0; k<3; ++k) {
      lseek(fd[k][j], 0, SEEK_SET);
    }
  }

  // Besides correcting each capacitor cell, we also apply a correction to ADC channels using the same data
  // Probably not as fundamental, but it does not hurt

  cout << "Sample calibration" << endl;

  memset(n, 0, 3*19*1024*sizeof(Int_t));
  flag=1;
  evnt=0;
  while(flag) {
    ++evnt;
    if(evmax_sample && evnt==evmax_sample+1) break;
    if(!(evnt%1000)) cout << "Event " << evnt << endl;
    flag=0;
    uint16_t startcell[3];
    Float_t fbuf[3][1024];
    for(Int_t j=0; j<nch; ++j) {
      for(Int_t k=0; k<3; ++k) {
        if(!read(fd[k][j],&startcell[k],sizeof(uint16_t))) continue;
        read(fd[k][j],&fbuf[k][0],sizeof(Float_t)*1024);
        for(Int_t i=0; i<1024; ++i) {
          Int_t index = (i+startcell[k])%1024;
          ++n[k][j][i];
          Double_t y=p0[j][index]+p1[j][index]*(fbuf[k][i]-p2[j][index]);
          ha[k][j]->SetBinContent(i+1,y+ha[k][j]->GetBinContent(i+1));
        }
        ++flag;
      }
    }
  }  

  for(Int_t j=0; j<nch; ++j) {
    for(Int_t k=0; k<3; ++k) {
      for(Int_t i=0; i<1024; ++i) {
        ha[k][j]->SetBinContent(i+1,ha[k][j]->GetBinContent(i+1)/n[k][j][i]);
      }
    }
  }

  // Make same fits, and save calibration constants

  Double_t pa0[19][1024];
  Double_t pa1[19][1024];
  Double_t pa2[19][1024];

  for(Int_t j=0; j<nch; ++j) {
    TString filename=Form("calib_Marco/calib_sample_%d.txt",j);
    FILE *file=fopen(filename.Data(),"w");

    for(Int_t i=0; i<1024; ++i) {
      Double_t a[3];
      for(Int_t k=0; k<3; ++k) {
        a[k]=ha[k][j]->GetBinContent(i+1);
        gr->SetPoint(k,a[k],v[j][k]);
      }
      func->FixParameter(2,a[1]);
      gr->Fit("func","0Q");

      pa0[j][i] = func->GetParameter(0);
      pa1[j][i] = func->GetParameter(1);
      pa2[j][i] = func->GetParameter(2);
     
      fprintf(file,"%g %g %g\n", pa0[j][i], pa1[j][i], pa2[j][i]); 
    }
    fclose(file);
  } 

  // Now moving to time calibration, first with local time calibration

  Double_t dvr[19][1024];
  Int_t ndvr[19][1024];
  memset(dvr,0,19*1024*sizeof(Double_t));
  memset(ndvr,0,19*1024*sizeof(Int_t));
  Double_t dvf[19][1024];
  Int_t ndvf[19][1024];
  memset(dvf,0,19*1024*sizeof(Double_t));
  memset(ndvf,0,19*1024*sizeof(Int_t));

  cout << "Time local calibration" << endl;

  flag=1;
  evnt=0;
  Int_t cont[19];
  memset(cont,0,19*sizeof(Int_t));
  while(flag) {
    ++evnt;
    if(evmax_tc && evnt==evmax_tc+1) break;
    if(!(evnt%1000)) cout << "Event " << evnt << endl;
    flag=0;
    Float_t fbuf[19][1024];
    Float_t fbufg[19][1024];
    uint16_t startcell[19];
    uint16_t startcellg[19];
    for(Int_t j=0; j<nch; ++j) {
      // Read signal channel file and guard channel file (the latter used to remove the spikes)
      if(!read(fdt[j],&startcell[j],sizeof(uint16_t))) {cont[j]=1; continue;}
      read(fdt[j],&fbuf[j][0],sizeof(Float_t)*1024);
      if(!read(fdg[j],&startcellg[j],sizeof(uint16_t))) {cont[j]=1; continue;}
      read(fdg[j],&fbufg[j][0],sizeof(Float_t)*1024);
      for(Int_t i=0; i<1024; ++i) {
        Int_t index = (i+startcell[j])%1024;
        // Apply first voltage cell calibration
        Double_t y=p0[j][index]+p1[j][index]*(fbuf[j][i]-p2[j][index]);
        y = pa0[j][i]+pa1[j][i]*(y-pa2[j][i]);
        ht[j]->SetBinContent(i+1,y);

        Int_t indexg = (i+startcellg[j])%1024;
        // and second apply voltage ADC calibration in sequence
        Double_t yg=p0[guardch[j]][indexg]+p1[guardch[j]][indexg]*(fbufg[j][i]-p2[guardch[j]][indexg]);
        yg = pa0[guardch[j]][i]+pa1[guardch[j]][i]*(yg-pa2[guardch[j]][i]);
        hg[j]->SetBinContent(i+1,yg);

      }
      ++flag;
    }

    for(Int_t j=0; j<nch; ++j) {
      if(cont[j]) continue;
      // We don't care about the first few samples and about the last 50 samples
      for(Int_t i=5; i<970; ++i) {
        // First look for zero crossings on the rising front
        if( ( (ht[j]->GetBinContent(i)>0 && ht[j]->GetBinContent(i-1)<0) ||
               (ht[j]->GetBinContent(i)<0 && ht[j]->GetBinContent(i-1)>0) ) &&
             ht[j]->GetBinContent(i+1)>0 && ht[j]->GetBinContent(i-2)<0 ) {
          for(Int_t m=0; i+m<970; ++m) {
            // Skip if there's a spike on the guard channel
            if(fabs(hg[j]->GetBinContent(i+m)-hg[j]->GetBinContent(i+m-1))>0.005
 ||
               fabs(hg[j]->GetBinContent(i+m+1)-hg[j]->GetBinContent(i+m))>0.005
)
               continue;

            // Use only events around the zero up to 0.3 V
            if(ht[j]->GetBinContent(i+m)>0.3) break;
            dvr[j][(i+m-2+startcell[j])%1024] += ht[j]->GetBinContent(i+m)-ht[j]->GetBinContent(i+m-1);
            ++ndvr[j][(i+m-2+startcell[j])%1024];
          }
          for(Int_t m=-1; i+m>5; --m) {
            if(fabs(hg[j]->GetBinContent(i+m)-hg[j]->GetBinContent(i+m-1))>0.005 ||
               fabs(hg[j]->GetBinContent(i+m+1)-hg[j]->GetBinContent(i+m))>0.005)
               continue;
            // And around the zero down to -0.3 V
            if(ht[j]->GetBinContent(i+m)<-0.3) break;
            dvr[j][(i+m-2+startcell[j])%1024] += ht[j]->GetBinContent(i+m)-ht[j]->GetBinContent(i+m-1);
            ++ndvr[j][(i+m-2+startcell[j])%1024];
          }

        // Now repeat everything for the falling front
        } else if( ((ht[j]->GetBinContent(i)>0 && ht[j]->GetBinContent(i-1)<0) ||
                     (ht[j]->GetBinContent(i)<0 && ht[j]->GetBinContent(i-1)>0)) &&
                     ht[j]->GetBinContent(i+1)<0 && ht[j]->GetBinContent(i-2)>0 ) {

          for(Int_t m=0; i+m<970; ++m) {
            if(fabs(hg[j]->GetBinContent(i+m)-hg[j]->GetBinContent(i+m-1))>0.005 ||
               fabs(hg[j]->GetBinContent(i+m+1)-hg[j]->GetBinContent(i+m))>0.005)
               continue;
            if(ht[j]->GetBinContent(i+m)<-0.3) break;
            dvf[j][(i+m-2+startcell[j])%1024] += -ht[j]->GetBinContent(i+m)+ht[j]->GetBinContent(i+m-1);
            ++ndvf[j][(i+m-2+startcell[j])%1024];
          }
          for(Int_t m=-1; i+m>5; --m) {
            if(fabs(hg[j]->GetBinContent(i+m)-hg[j]->GetBinContent(i+m-1))>0.005 ||
               fabs(hg[j]->GetBinContent(i+m+1)-hg[j]->GetBinContent(i+m))>0.005)
               continue;
            if(ht[j]->GetBinContent(i+m)>0.3) break;
            dvf[j][(i+m-2+startcell[j])%1024] += -ht[j]->GetBinContent(i+m)+ht[j]->GetBinContent(i+m-1);
            ++ndvf[j][(i+m-2+startcell[j])%1024];
          }
        }
      }
    }
  

  }


  Double_t dtu[19][1024];
  Double_t dtd[19][1024];
  for(Int_t j=0; j<nch; ++j) {
    Double_t tvr=0;
    Double_t tvf=0;
    for(Int_t i=0; i<1024; ++i) {
      dvr[j][i] /= ndvr[j][i];
      dvf[j][i] /= ndvf[j][i];

      tvr += dvr[j][i];
      tvf += dvf[j][i];
    }
    for(Int_t i=0; i<1024; ++i) {
      dtu[j][i]=dvr[j][i]/tvr*timerange;
      dtd[j][i]=dvf[j][i]/tvf*timerange;
    }
  }

  // Determine cell widths by averaging rising and falling fronts (this cancels some residual voltage effect)
  Double_t dt[19][1024];
  for(Int_t j=0; j<nch; ++j) {
    for(Int_t i=0; i<1024; ++i) {
      dt[j][i]=0.5*(dtu[j][i]+dtd[j][i]);
    }
  }


  // Now moving to the global time calibration. Start with the rising front
  cout << "Time global calibration rising" << endl;

  Double_t dtr[19][1024]; 
  memcpy(dtr,dt,19*1024*sizeof(Double_t));

  Double_t ucorr[19][1024];
  Int_t nucorr[19][1024];
  memset(ucorr,0,19*1024*sizeof(Double_t));
  memset(nucorr,0,19*1024*sizeof(Int_t));

  flag=1;
  evnt=0;
  memset(cont,0,19*sizeof(Int_t));
  while(flag) {
    ++evnt;
    if(!(evnt%1000)) cout << "Event " << evnt << endl;
    if(evmax_tgc && evnt==evmax_tgc+1) break;
    flag=0;
    Float_t fbuf[19][1024];
    Float_t fbufg[19][1024];
    uint16_t startcell[19];
    uint16_t startcellg[19];

    for(Int_t j=0; j<nch; ++j) {
      // Read signal and guard data
      if(!read(fdtgl[j],&startcell[j],sizeof(uint16_t))) {cont[j]=1; continue;}
      read(fdtgl[j],&fbuf[j][0],sizeof(Float_t)*1024);
      if(!read(fdglg[j],&startcellg[j],sizeof(uint16_t))) {cont[j]=1; continue;}
      read(fdglg[j],&fbufg[j][0],sizeof(Float_t)*1024);

      // Apply voltage calibration for cells and ADCs 
      for(Int_t i=0; i<1024; ++i) {
        Int_t index = (i+startcell[j])%1024;
        Double_t y=p0[j][index]+p1[j][index]*(fbuf[j][i]-p2[j][index]);
        y = pa0[j][i]+pa1[j][i]*(y-pa2[j][i]);
        ht[j]->SetBinContent(i+1,y);

        Int_t indexg = (i+startcellg[j])%1024;
        Double_t yg=p0[guardch[j]][indexg]+p1[guardch[j]][indexg]*(fbufg[j][i]-p2[guardch[j]][indexg]);
        yg = pa0[guardch[j]][i]+pa1[guardch[j]][i]*(yg-pa2[guardch[j]][i]);
        hg[j]->SetBinContent(i+1,yg);

      }
      ++flag;
    }

    for(Int_t j=0; j<nch; ++j) {
      if(cont[j]) continue;
      Int_t zeroc[100];
      memset(zeroc,0,100*sizeof(Int_t));
      Int_t nzeroc=0;

      // Determine the zero crossings of the waveform
      for(Int_t i=5; i<970; ++i) {
        if( ( (ht[j]->GetBinContent(i)>0 && ht[j]->GetBinContent(i-1)<0) ||
              (ht[j]->GetBinContent(i)<0 && ht[j]->GetBinContent(i-1)>0) ) &&
             ht[j]->GetBinContent(i+1)>0 && ht[j]->GetBinContent(i-2)<0 ) {
           zeroc[nzeroc]=i;
           ++nzeroc;
        }
      }

      // Loop on all zero crossings (they are much less than 100...)
      for(Int_t m=0; m<100; m++) {
        for(Int_t i=0; ; ++i) {
          if(i+m+1 >= nzeroc) break;
          // Discard spikes if any
          if(fabs(hg[j]->GetBinContent(zeroc[i])-hg[j]->GetBinContent(zeroc[i]-1))>0.005 ||
             fabs(hg[j]->GetBinContent(zeroc[i])-hg[j]->GetBinContent(zeroc[i]+1))>0.005)
            continue;

          Double_t Uk = ht[j]->GetBinContent(zeroc[i]);
          Double_t dUkm1 = Uk-ht[j]->GetBinContent(zeroc[i]-1);
          Double_t dtkm1 = dtr[j][(zeroc[i]-2+startcell[j])%1024];

          if(fabs(hg[j]->GetBinContent(zeroc[i+1+m]-1)-hg[j]->GetBinContent(zeroc[i+1+m]))>0.005 ||
             fabs(hg[j]->GetBinContent(zeroc[i+1+m]-1)-hg[j]->GetBinContent(zeroc[i+1+m]-2))>0.005)
            continue;

          Double_t Uq = ht[j]->GetBinContent(zeroc[i+1+m]-1);
          Double_t dUq = ht[j]->GetBinContent(zeroc[i+1+m])-Uq;
          Double_t dtq = dtr[j][(zeroc[i+1+m]-2+startcell[j])%1024];

          Double_t tcor = Uk/dUkm1*dtkm1 - Uq/dUq*dtq;

          Double_t tkq=0;
          for(Int_t k=zeroc[i]-1; k<zeroc[i+m+1]-2; ++k) tkq += dtr[j][(k+startcell[j])%1024];

          // Determine correction for global calibration (see relevant publication by Stricker-Shaver et al.
          Double_t uCor = (m+1)/(frequency*1e-3)/(tkq+tcor);

          // determine average correction (this is different from Stricker-Shaver et al., which apply an
          // interative procedure
          for(Int_t k=zeroc[i]-1; k<zeroc[i+m+1]-2; ++k) { ucorr[j][(k+startcell[j])%1024] += uCor; ++nucorr[j][(k+startcell[j])%1024]; }
        }
      }
    }
  }

  // Redo the same for the falling fronts
  cout << "Time global calibration falling" << endl;

  for(Int_t j=0; j<nch; ++j) {
    lseek(fdtgl[j], 0, SEEK_SET);
    lseek(fdglg[j], 0, SEEK_SET);
  }

  Double_t dtf[19][1024];
  memcpy(dtf,dt,19*1024*sizeof(Double_t));

  Double_t ucorf[19][1024];
  Int_t nucorf[19][1024];
  memset(ucorf,0,19*1024*sizeof(Double_t));
  memset(nucorf,0,19*1024*sizeof(Int_t));

  flag=1;
  evnt=0;
  memset(cont,0,19*sizeof(Int_t));
  while(flag) {
    ++evnt;
    if(!(evnt%1000)) cout << "Event " << evnt << endl;
    if(evmax_tgc && evnt==evmax_tgc+1) break;
    flag=0;
    Float_t fbuf[19][1024];
    Float_t fbufg[19][1024];
    uint16_t startcell[19];
    uint16_t startcellg[19];
    for(Int_t j=0; j<nch; ++j) {
      if(!read(fdtgl[j],&startcell[j],sizeof(uint16_t))) {cont[j]=1; continue;}
      read(fdtgl[j],&fbuf[j][0],sizeof(Float_t)*1024);
      if(!read(fdglg[j],&startcellg[j],sizeof(uint16_t))) {cont[j]=1; continue;}
      read(fdglg[j],&fbufg[j][0],sizeof(Float_t)*1024);

      for(Int_t i=0; i<1024; ++i) {
        Int_t index = (i+startcell[j])%1024;
        Double_t y=p0[j][index]+p1[j][index]*(fbuf[j][i]-p2[j][index]);
        y = pa0[j][i]+pa1[j][i]*(y-pa2[j][i]);
        ht[j]->SetBinContent(i+1,y);

        Int_t indexg = (i+startcellg[j])%1024;
        Double_t yg=p0[guardch[j]][indexg]+p1[guardch[j]][indexg]*(fbufg[j][i]-p2[guardch[j]][indexg]);
        yg = pa0[guardch[j]][i]+pa1[guardch[j]][i]*(yg-pa2[guardch[j]][i]);
        hg[j]->SetBinContent(i+1,yg);

      }
      ++flag;

    }

    for(Int_t j=0; j<nch; ++j) {
      if(cont[j]) continue;
      Int_t zeroc[100];
      memset(zeroc,0,100*sizeof(Int_t));
      Int_t nzeroc=0;

      for(Int_t i=5; i<970; ++i) {
        if( ((ht[j]->GetBinContent(i)>0 && ht[j]->GetBinContent(i-1)<0) ||
             (ht[j]->GetBinContent(i)<0 && ht[j]->GetBinContent(i-1)>0)) &&
              ht[j]->GetBinContent(i+1)<0 && ht[j]->GetBinContent(i-2)>0 ) {
           zeroc[nzeroc]=i;
           ++nzeroc;
        }
      }

      for(Int_t m=0; m<100; m++) {
        for(Int_t i=0; ; ++i) {
          if(i+m+1 >= nzeroc) break;

          if(fabs(hg[j]->GetBinContent(zeroc[i])-hg[j]->GetBinContent(zeroc[i]-1))>0.005 ||
             fabs(hg[j]->GetBinContent(zeroc[i])-hg[j]->GetBinContent(zeroc[i]+1))>0.005)
            continue;

          Double_t Uk = -ht[j]->GetBinContent(zeroc[i]);
          Double_t dUkm1 = Uk+ht[j]->GetBinContent(zeroc[i]-1);
          Double_t dtkm1 = dtf[j][(zeroc[i]-2+startcell[j])%1024];

          if(fabs(hg[j]->GetBinContent(zeroc[i+1+m]-1)-hg[j]->GetBinContent(zeroc[i+1+m]))>0.005 ||
             fabs(hg[j]->GetBinContent(zeroc[i+1+m]-1)-hg[j]->GetBinContent(zeroc[i+1+m]-2))>0.005)
            continue;

          Double_t Uq = -ht[j]->GetBinContent(zeroc[i+m+1]-1);
          Double_t dUq = -ht[j]->GetBinContent(zeroc[i+m+1])-Uq;
          Double_t dtq = dtf[j][(zeroc[i+m+1]-2+startcell[j])%1024];

          Double_t tcor = Uk/dUkm1*dtkm1 - Uq/dUq*dtq;

          Double_t tkq=0;
          for(Int_t k=zeroc[i]-1; k<zeroc[i+m+1]-2; ++k) tkq += dtf[j][(k+startcell[j])%1024];

          Double_t uCor = (m+1)/(frequency*1e-3)/(tkq+tcor);

          for(Int_t k=zeroc[i]-1; k<zeroc[i+m+1]-2; ++k) { ucorf[j][(k+startcell[j])%1024] += uCor; ++nucorf[j][(k+startcell[j])%1024]; }
        }   
      }
    }
  }

  // Apply global time calibration to the cell widths determined by the local time calibration
  for(Int_t j=0; j<nch; ++j) {
    for(Int_t i=0; i<1024; ++i) {
      dt[j][i] *= 0.5*(ucorr[j][i]/nucorr[j][i] + ucorf[j][i]/nucorf[j][i]);
    }
  }

  // Finally store the cell widths
  for(Int_t j=0; j<nch; ++j) {
    TString filename=Form("calib_Marco/calib_time_%d.txt",j);
    FILE *filetout=fopen(filename.Data(),"w");
    for(Int_t i=0; i<1024; ++i) {
      cout << j << " " << i << " " << dt[j][i] <<  " " << " " << dtu[j][i] << " " << dtd[j][i] << " " << ndvr[j][i] << " " << ndvf[j][i] << " " << endl;
      fprintf(filetout,"%g\n", dt[j][i]);
    }
    fclose(filetout);
  }

}

