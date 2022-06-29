// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include <unordered_map>
#include "TROOT.h"

// Standard includes
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>

// Lorenzo includes
#include "Analysis_Digitizer.h"

using namespace std;

int main()
{

  bool enable_MT = true;
  if(enable_MT){
   ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(8);
    cout << "Multithread is enabled"<< endl;
  }
  //bool showFFT=false;  
  bool showFFT=true;
  bool doFit = false;
  
  //  int samplesfact=1;  //divider to sample rate for waveform memorization

  TFile *OutputFile = new TFile("RunXX.root","recreate");

  ifstream InputCARD("Input_Folder_Digitizer_wfm.txt");
  
  if(!InputCARD)
    {
      cout << "Error: could not find the InputCARD" << endl;
      return (0);
    }



  // read the active numnber of channels
  int ntmp;
  int  Nth = 0;
  int np_Max = 300; // 800;
  int np_acq = 0; // points in acquisizion, set to == size  
  int np_offset = 100; //100; //starting point for waveform analysis

 
  //  TRandom *xi = new TRandom();
//samples in t
  InputCARD >> pip >> ntmp;
  nchro = ntmp;
  
  bool FWF2 = false; //file weightfield

  //Waveform variables



  //  const Int_t max_samples=20000;
  // const Int_t max_samplesrec=int(max_samples/samplesfact)+100:

  // SetRootFile(doFit,showFFT,"RunXX.root");
  
  // Tree definition
  // TTree *w = new TTree("Wave","Wave");
  // w->SetDirectory(OutputFile);
  
  TTree *OutTree = new TTree("Analysis","Analysis");
  OutTree->SetDirectory(OutputFile);
 

  // Analysis branches
  OutTree->Branch("ntrig",&ntrig,"ntrig/I");
  OutTree->Branch("event",&event,"event/I");
  OutTree->Branch("nrun",&nrun,"nrun/I");
  OutTree->Branch("nchro",&nchro,"nchro/I");
  OutTree->Branch("npos",&npos,"npos/I");
  
  //  OutTree->Branch("WF2gain",&WF2gain,"WF2gain/D");
  // OutTree->Branch("WF2xpos",&WF2xpos,"WF2xpos/D");
  // OutTree->Branch("WF2angle",&WF2angle,"WF2angle/D");

  OutTree->Branch("XPos",&XPos,"XPos/D");
  OutTree->Branch("YPos",&YPos,"YPos/D");
  OutTree->Branch("XPad",XPad,"XPad[nchro]/D");
  OutTree->Branch("YPad",YPad,"YPad[nchro]/D");
  OutTree->Branch("Dist",Dist,"Dist[nchro]/D"); //disatanza dei punti dai pad
  

  OutTree->Branch("t_bck",t_bck,"t_bck[nchro]/D"); //beginning od the signal
  OutTree->Branch("t_pul",t_pul,"t_pul[nchro]/D"); //end of sìgnal

  OutTree->Branch("maxD",maxD,"maxD[nchro]/D"); 
  OutTree->Branch("area",area,"area[nchro]/D");
  OutTree->Branch("ampl",ampl,"ampl[nchro]/D");
  OutTree->Branch("Toffset",Toffset,"Toffset[nchro]/D"); //forse vuoto
  OutTree->Branch("ampl_chi2",ampl_chi2,"ampl_chi2[nchro]/D"); //chi2 of the ampl fit
  OutTree->Branch("dVdt3070",dVdt3070,"dVdt3070[nchro]/D"); 
  OutTree->Branch("dVdt1030",dVdt1030,"dVdt1030[nchro]/D");
  OutTree->Branch("dVdt2080",dVdt2080,"dVdt2080[nchro]/D");

  OutTree->Branch("bck",bck,"bck[nchro]/D"); //baseline?
  OutTree->Branch("max_bck_before",max_bck_before,"max_bck_before[nchro]/D"); 
  OutTree->Branch("rms_bck_before",rms_bck_before,"rms_bck_before[nchro]/D");
  OutTree->Branch("max_bck_after",max_bck_after,"max_bck_after[nchro]/D");
  OutTree->Branch("rms_bck_after",rms_bck_after,"rms_bck_after[nchro]/D");


  OutTree->Branch("t_cent",t_cent,"t_cent[nchro]/D"); //time at centroid of the signal
  
  OutTree->Branch("t_level10",t_level10,"t_level10[nchro]/D"); //probably does not make sense for the digitizer 
  OutTree->Branch("t_level15",t_level15,"t_level15[nchro]/D"); 
  OutTree->Branch("t_level20",t_level20,"t_level20[nchro]/D");
  OutTree->Branch("t_level30",t_level30,"t_level30[nchro]/D");
  OutTree->Branch("t_level40",t_level40,"t_level40[nchro]/D");
  OutTree->Branch("t_level50",t_level50,"t_level50[nchro]/D");
  OutTree->Branch("t_level60",t_level60,"t_level60[nchro]/D");
  OutTree->Branch("t_level80",t_level80,"t_level80[nchro]/D");
  OutTree->Branch("t_level100",t_level100,"t_level100[nchro]/D");
  OutTree->Branch("t_level200",t_level200,"t_level200[nchro]/D");
  OutTree->Branch("t_level300",t_level300,"t_level300[nchro]/D");

  OutTree->Branch("trail_t_level10",trail_t_level10,"trail_t_level10[nchro]/D"); //on the other part of the signal
  OutTree->Branch("trail_t_level20",trail_t_level20,"trail_t_level20[nchro]/D");
  OutTree->Branch("trail_t_level30",trail_t_level30,"trail_t_level30[nchro]/D");
  OutTree->Branch("trail_t_level40",trail_t_level40,"trail_t_level40[nchro]/D");
  OutTree->Branch("trail_t_level50",trail_t_level50,"trail_t_level50[nchro]/D");
  OutTree->Branch("trail_t_level60",trail_t_level60,"trail_t_level60[nchro]/D");
  OutTree->Branch("trail_t_level80",trail_t_level80,"trail_t_level80[nchro]/D"); 
  OutTree->Branch("trail_t_level100",trail_t_level100,"trail_t_level100[nchro]/D");
  OutTree->Branch("trail_t_level120",trail_t_level120,"trail_t_level120[nchro]/D");
  OutTree->Branch("trail_t_level140",trail_t_level140,"trail_t_level140[nchro]/D");
  OutTree->Branch("trail_t_level160",trail_t_level160,"trail_t_level160[nchro]/D");
  OutTree->Branch("trail_t_level180",trail_t_level180,"trail_t_level180[nchro]/D");


  OutTree->Branch("cfd05",cfd05,"cfd05[nchro]/D");
  OutTree->Branch("cfd10",cfd10,"cfd10[nchro]/D");
  OutTree->Branch("cfd15",cfd15,"cfd15[nchro]/D");
  OutTree->Branch("cfd20",cfd20,"cfd20[nchro]/D");
  OutTree->Branch("cfd25",cfd25,"cfd25[nchro]/D");
  OutTree->Branch("cfd30",cfd30,"cfd30[nchro]/D");
  OutTree->Branch("cfd40",cfd40,"cfd40[nchro]/D");
  OutTree->Branch("cfd50",cfd50,"cfd50[nchro]/D");
  OutTree->Branch("cfd60",cfd60,"cfd60[nchro]/D");
  OutTree->Branch("cfd70",cfd70,"cfd70[nchro]/D");
  OutTree->Branch("cfd80",cfd80,"cfd80[nchro]/D");
  OutTree->Branch("cfd90",cfd90,"cfd90[nchro]/D");

  OutTree->Branch("Ncfd05",Ncfd05,"Ncfd05[nchro]/I"); //N at which the signal reaches 5 pc
  OutTree->Branch("Ncfd10",Ncfd10,"Ncfd10[nchro]/I");
  OutTree->Branch("Ncfd15",Ncfd15,"Ncfd15[nchro]/I");
  OutTree->Branch("Ncfd20",Ncfd20,"Ncfd20[nchro]/I");
  OutTree->Branch("Ncfd25",Ncfd25,"Ncfd25[nchro]/I");
  OutTree->Branch("Ncfd30",Ncfd30,"Ncfd30[nchro]/I");
  OutTree->Branch("Ncfd40",Ncfd40,"Ncfd40[nchro]/I");
  OutTree->Branch("Ncfd50",Ncfd50,"Ncfd50[nchro]/I");
  OutTree->Branch("Ncfd60",Ncfd60,"Ncfd60[nchro]/I");
  OutTree->Branch("Ncfd70",Ncfd70,"Ncfd70[nchro]/I");
  OutTree->Branch("Ncfd80",Ncfd80,"Ncfd80[nchro]/I");
  OutTree->Branch("Ncfd90",Ncfd90,"Ncfd90[nchro]/I");

  OutTree->Branch("trail_Ncfd10",trail_Ncfd10,"trail_Ncfd10[nchro]/I");
  OutTree->Branch("trail_Ncfd30",trail_Ncfd30,"trail_Ncfd30[nchro]/I");
  OutTree->Branch("trail_Ncfd50",trail_Ncfd50,"trail_Ncfd50[nchro]/I");
  OutTree->Branch("trail_Ncfd70",trail_Ncfd70,"trail_Ncfd70[nchro]/I");
  OutTree->Branch("trail_Ncfd90",trail_Ncfd90,"trail_Ncfd90[nchro]/I");

  
  OutTree->Branch("t_rms3",t_rms3,"t_rms3[nchro]/D"); //time at 3xRMS
  OutTree->Branch("t_rms5",t_rms5,"t_rms5[nchro]/D");
  OutTree->Branch("trail_cfd10",trail_cfd30,"trail_cfd10[nchro]/D");
  OutTree->Branch("trail_cfd30",trail_cfd30,"trail_cfd30[nchro]/D");
  OutTree->Branch("trail_cfd50",trail_cfd50,"trail_cfd50[nchro]/D");
  OutTree->Branch("trail_cfd70",trail_cfd70,"trail_cfd70[nchro]/D");
  OutTree->Branch("trail_cfd90",trail_cfd90,"trail_cfd90[nchro]/D");

  OutTree->Branch("samples",samples,"samples[nchro]/I"); //n of samples of the waveform
  //  OutTree->Branch("samplesrec",samples,"samples[nchro]/I");
  OutTree->Branch("t_max",t_max,"t_max[nchro]/D"); //time at maximum
 
  
  if(doFit)
    {
      OutTree->Branch("t_zero",t_zero,"t_zero[nchro]/D");
      OutTree->Branch("rise_lin0",rise_lin0,"rise_lin0[nchro]/D");
      OutTree->Branch("rise_lin1",rise_lin1,"rise_lin1[nchro]/D");
      OutTree->Branch("rise_lin_chi2",&rise_lin_chi2,"rise_lin_chi2[nchro]/D");
      OutTree->Branch("rise_exp0",rise_exp0,"rise_exp0[nchro]/D");
      OutTree->Branch("rise_exp1",rise_exp1,"rise_exp1[nchro]/D");
      OutTree->Branch("rise_exp_chi2",rise_exp_chi2,"rise_exp_chi2[nchro]/D");
    }

  cout << "nchro = " << nchro << endl;
  for(i=0;i<nchro;i++)
    {
      
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "samples_" << i;
      leafl << "samples_" << i <<"/I";
      OutTree->Branch(leaf.str().c_str(),&samples[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "samplesrec_" << i;
      leafl << "samplesrec_" << i <<"/I";
      OutTree->Branch(leaf.str().c_str(),&samplesrec[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
    
      leaf << "amp" << i; //wfm as it comes from digi
      leafl << "amp" << i <<"[samplesrec_" << i << "]/D";
      OutTree->Branch(leaf.str().c_str(),&amprec[i][0],leafl.str().c_str());
      
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "m_amp" << i; //positive and baseline suptracted
      leafl << "m_amp" << i <<"[samplesrec_" << i << "]/D";
	    //leafl << "m_amp" << i <<"[samples[" << i << "]]/D";
      OutTree->Branch(leaf.str().c_str(),&m_amprec[i][0],leafl.str().c_str());
      cout << " =  " << leafl.str().c_str() << endl;

      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "der_amp" << i; //derivative
      leafl << "der_amp" << i <<"[samplesrec_" << i << "]/D";
      OutTree->Branch(leaf.str().c_str(),&d_amprec[i][0],leafl.str().c_str());


  //OutTree->Branch("gum_amp",gum_amp,"gum_amp[samples_1]/D");

      if(showFFT)
	{
	  
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "FFT_abs" << i;
	  leafl << "FFT_abs" << i <<"[samplesrec_" << i << "]/D";
	  OutTree->Branch(leaf.str().c_str(),&FFT_abs[i][0],leafl.str().c_str());
	  
	  
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "FFT_real" << i;
	  leafl << "FFT_real" << i <<"[samplesrec_" << i << "]/D";
	  OutTree->Branch(leaf.str().c_str(),&FFT_real[i][0],leafl.str().c_str());
	  
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "FFT_comp" << i;
	  leafl << "FFT_comp" << i <<"[samplesrec_" << i << "]/D";
	  OutTree->Branch(leaf.str().c_str(),&FFT_comp[i][0],leafl.str().c_str());
	}


      
    }
						       
 OutTree->Branch("time",timerec,"time[samples_0]/D"); 
  //  OutTree->Branch("freq",freq,"freq[samples_0]/D");
  
  
  
 cout << "Branch settings done" << endl;
  // float ie, ih, ieg, ihg;
 float Totie, Totih, Totieg, Totihg; 
 bool FToffeeCh[4] = {0,0,0,0};
 

  bool Interpolation = true;
  //bool Interpolation = false;
 
  event = 0;


  
  // Levels array
  const int NArrayL = 11;
  double ThresholdL[NArrayL] = {10, 15,20,30,40,50,60,80,100,200,300};
  double TValuesL[NArrayL];
  int NTValuesL[NArrayL];
  int NArr = 0;
  
  // Levels trail array
  const int TrailNArrayL = 12;
  double TrailThresholdL[TrailNArrayL] = {10,20,30,40,50,60,80,100,120,140,160,180};
  double TrailTValuesL[TrailNArrayL];
  int TrailNTValuesL[TrailNArrayL];
  int TrailNArr = 0;
  
  // fraction array
  const int NArrayF = 12;
  double ThresholdF[NArrayF] = {0.05, 0.1,0.15,0.20,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  double TValuesF[NArrayF];
  int NTValuesF[NArrayF];
  
  // fraction trail array
  const int TrailNArrayF = 5;
  double TrailThresholdF[TrailNArrayF] = {0.1,0.3,0.5,0.7,0.9};
  double TrailTValuesF[TrailNArrayF];
  int TrailNTValuesF[TrailNArrayF];
  
  npos= 0;
  
  //old
  //  int XPosArray[18] = {0,240,0,240,140,440,140,340,540,640,0,0,440,640,0,0,0,0};
  // int YPosArray[18] = {0,609,0,263,436,263,90,90,90,609,0,0,609,263,0,0,0,0};
  
  //new
 // int YPosArray[18] = {0,240,0,240,140,440,140,340,540,640,0,0,440,640,0,0,0,0};
 // int XPosArray[18] = {0,609,0,263,436,263,90,90,90,609,0,0,609,263,0,0,0,0};
  

   int XPosArray[18] = {1680,1680,1230,0,780,780,330,330,330,0,780,0, 0, 1230, 0, 1680, 0,0};
   int YPosArray[18] = {451,901,451,0, 901,451,901,451,1,  0,  1,0, 0,1,0,1,0,0};

  for(i=0;i<nchro;i++)
    {
      XPad[i] = XPosArray[i];
      YPad[i] = YPosArray[i];
      cout << "Pad " << i << " is in " << XPad[i]<<","<<YPad[i] << endl;
    }
	
	//Loop on different runs
    while(1)
	  {
	    pip="";
	    toffee="";
	    if(InputCARD.eof())
	      break;
	    
	    //Read run variables
	    InputCARD >> pip >> nrun >> pip >> ntrig >>  pip >> Max_eventperpoint >>pip >> Max_skippos >> pip >> nchro >> pip >> MaxEvt;
	    if (nrun == -1 ) break;
	    cout << endl << endl << "Run " << nrun << "\t " << nchro << " channels" << " Events to be analized = " << MaxEvt - ntrig << " Event per point = " << Max_eventperpoint << " Skipping points # " << Max_skippos << endl;
      //Note: ntrig is the number of the FIRST trigger of the run

      //     if ( (nrun >= 23000 && nrun <= 24000 ) || (nrun >= 33000 && nrun <= 34000 ) )
      //	if (np_offset = 600;
   
      //Inizialization of data run variables
      ntrig=ntrig-1;
     
      ostringstream convert;
      convert << nrun ;
      String_nrun = convert.str();
      

      //Inizialization of program run variables
      running=1;
      char title[200];
      //Inizialization of channel variables for this run
 
      InputCARD >> toffee >> InputNA >> pip >>nme;  //NOTE: InputNAME is the complete path + file name without trigger number and final .txt
      sprintf(title,"%s",InputNA.c_str());

      cout << endl << " File name:" << endl << title << " averaging " << nme << " points " <<endl;

      
      //      FreqCut[1] = 0;
      if (nme<0)
	{
	  FreqCut[0] = -nme*1.e6;
	  //FreqCut[1] = 1.e9;
	  // FreqCut[1] = 0;
	  nmedia[1] = 3;
	  nmedia[2] = 3;
	  cout << endl << "Channel 0  with Freq. Cut = " << FreqCut[0]*1e-6 << " MHz" <<endl;
	  cout << "Channel 1  with  Freq. Cut  = " <<  FreqCut[1]*1e-6 << " MHz"  <<endl;
	}
	
      

      TFile *f = new TFile(title);
      TTree *tree = (TTree*)f->Get("wfm");
  

   //  tree->SetBranchAddress("i_timestamp", &i_timestamp, &b_i_timestamp);
   // tree->SetBranchAddress("i_current", &i_current, &b_i_current);
      //variables from input root files
      tree->SetBranchAddress("w0", &w0, &b_w0);
      tree->SetBranchAddress("w1", &w1, &b_w1);
      tree->SetBranchAddress("w2", &w2, &b_w2);
      tree->SetBranchAddress("w3", &w3, &b_w3);
      tree->SetBranchAddress("w4", &w4, &b_w4);
      tree->SetBranchAddress("w5", &w5, &b_w5);
      tree->SetBranchAddress("w6", &w6, &b_w6);
      tree->SetBranchAddress("w7", &w7, &b_w7);
      tree->SetBranchAddress("w8", &w8, &b_w8);
      tree->SetBranchAddress("w9", &w9, &b_w9);
      tree->SetBranchAddress("w10", &w10, &b_w10);
      tree->SetBranchAddress("w11", &w11, &b_w11);
      tree->SetBranchAddress("w12", &w12, &b_w12);
      tree->SetBranchAddress("w13", &w13, &b_w13);
      tree->SetBranchAddress("w14", &w14, &b_w14);
      tree->SetBranchAddress("w15", &w15, &b_w15);
      tree->SetBranchAddress("trg0", &trg0, &b_trg0);
      tree->SetBranchAddress("trg1", &trg1, &b_trg1);
      tree->SetBranchAddress("pos", &pos, &b_pos);
      tree->SetBranchAddress("size", &size, &b_size);
      

      
      
      Long64_t nentries = tree->GetEntries();
      
      cout << "Number of events = " << nentries << endl;

      
      if (MaxEvt == -1) MaxEvt = nentries;
      string mystr;
      

      string RunBit;
      string RunRoot;
      RunRoot = String_nrun[0];


     
      emptycount=0;

      skippos =0;

      // here we are looping on each event
      while(running)  //running is set on 1 as long as all the channels have data for a trigger event
	{
	 
	  ntrig++;
	  event++;

	  //	   cout << "Event = " << event << " at trigger = " << ntrig << endl;
	  
	  if (MaxEvt > 0 && ntrig >MaxEvt) goto Write;
	  if(event==1)
	    cout << "Event = " << event << " at trigger = " << ntrig << endl;
	  if(ntrig%200==0)
	    cout << "Event = " << event << " at trigger = " << ntrig << endl;
	  tree->GetEntry(ntrig);
	  


	  //  if(pos[0][0]>1000 || pos[0][1]>1000) continue; 
 

	  if ( pos[0][0] != XPos || pos[0][1] != YPos)
	      
	    {
	      eventperpoint = 0;
	      XPos = pos[0][0];
	      YPos = pos[0][1];
	      // if(npos%50==0)
	      //	cout << "Position = " << npos << " at " << XPos <<","<<YPos << endl;
	      if ( skippos>=Max_skippos)
		{
		  skippos = 0;
		  npos++;
		  if(npos%50==0)
		    cout << "Position = " << npos << " at " << XPos <<","<<YPos << " " << skippos << " " << Max_skippos << endl;
		}
	      else if(skippos<Max_skippos)
		{
		  skippos++;
		}
	    }
	  if (skippos !=0) continue;
	      
	  eventperpoint++;

	  if (eventperpoint>Max_eventperpoint && Max_eventperpoint>-1) continue;



	  
	  for(i=0;i<nchro;i++)
	    {

	      FreqCut[i] = 0;

	      if (nme<0)
		  FreqCut[i] = -nme*1.e6;
	      else
		nmedia[i] = nme;		

	      samples[i]=0;
	      ampl_chi2[i]=0;
	      AMax[i]=0;
	      NMax[i]=0;
	      RunBit = String_nrun[i+1];
	      

	      samples[i]=0;
	      tempsum=0;
	      np=0; 
	      nprec=0;
	      Totie = 0;
	      Totih = 0;
	      Totieg = 0;
	      Totihg = 0;
	      WF2gain = 0;
	    	      

	      // reading channel i in the wfm root tree
	     
	      if (i ==0 && ntrig ==0)
		{
		  cout << "Number of samples in the acquisition = " << size << endl;
		  //	  cout << "Position:  x =  " << pos[0][0] << " y =  " << pos[0][1] <<endl;
		  
		  np_acq = size;

		  //		  if (np_acq>np_Max) np_Max = np_acq-1;// np_offset= 600;
		  cout << "Number of samples in the analysis = " << np_Max << " starting at samples " << np_offset << endl;
		  if (np_Max+np_offset>np_acq)
		    {
		      cout << "Number of samples in analysis larger than in the acquisition = " << np_Max+np_offset  << " > " << np_acq << endl;
		      cout << "The program STOPS" << endl;
		      return 1;
		    }		  				
		    
		}

	      Dist[i] = sqrt( pow(XPad[i]-XPos,2)+pow(YPad[i]-YPos,2));

	      //	      cout << "check 0" << endl;
	      
	      for (int npp =np_offset;npp<np_Max+np_offset;npp++)
		{
		  timeS[np] =  0.2*(npp-np_offset);
		  if (np ==0 ) Toffset[i] = 0;
		  
		  if (i== 0)		     
		      amp[i][np] = w0[0][npp];		    
		  else if (i ==1 )		    
		      amp[i][np] = w1[0][npp];		    
		  else if (i ==2 )		    
		      amp[i][np] = w2[0][npp];		    
		  else if (i ==3 )		    
		      amp[i][np] = w3[0][npp];
		  else if (i ==4 )		    
		      amp[i][np] = w4[0][npp];		    
		  else if (i ==5 )		    
		      amp[i][np] = w5[0][npp];		    
		  else if (i ==6 )		    
		      amp[i][np] = w6[0][npp];
		  else if (i ==7 )		    
		      amp[i][np] = w7[0][npp];		    
		  else if (i ==8 )		    
		      amp[i][np] = w8[0][npp];		    
		  else if (i ==9 )		    
		      amp[i][np] = w9[0][npp];
		  else if (i ==10 )		    
		      amp[i][np] = w10[0][npp];		    
		  else if (i ==11)		    
		      amp[i][np] = w11[0][npp];		    
		  else if (i == 12)		    
		      amp[i][np] = w12[0][npp];
		  else if (i ==13 )		    
		    amp[i][np] = w13[0][npp];		    
		  else if (i ==14)		    
		    amp[i][np] = w14[0][npp];		    
		  else if (i ==15)		    
		    amp[i][np] = w15[0][npp];
		  else if (i ==16)		    
		    amp[i][np] = trg0[0][npp];
		  else if (i ==17)		    
		    amp[i][np] = trg1[0][npp];		  
		  np++;
		}



	      if ( np < np_Max-10 && event<20) cout << "Warning: file reading ended too soon: read point np = " << np << " instead of np_max = " << np_Max << endl;
	      
	      DT = (timeS[10]-timeS[9]); // delta T in nanosecond
	      
	      samples[i]=np-1;
	      	      
	      AMax[i]=0;
			     
	      for (np = 0;np<samples[i];np++)
		// for (np = 20/DT;np<25/DT;np++)
		{
		  
		  if (fabs(amp[i][np]-amp[i][0])>AMax[i]) 
		    {
		      AMax[i] = fabs(amp[i][np]-amp[i][0]); 
		      TMax[i] = timeS[np];
		      polarity[i] = (amp[i][np]-amp[i][0])/fabs(amp[i][np]-amp[i][0]);
		      NMax[i] = np;		      		     

		    }

		}
	      

	      if (ntrig == 2  && i == 0) cout << "Number of points: " << np  << " with a time step of " << DT << " [ns] " << endl;

	      if (ntrig <20 ) cout << "Event: " << event << " channel " << i << " has maximum value of " << AMax[i] << " [mV] at  " <<  TMax[i] << " [ns] " << " on sample " << NMax[i] << " polarity = " << polarity[i] << " at position x-y " << pos[0][0] <<", " << pos[0][1] <<  endl;
	      if ( TMax[i] ==0 ||  AMax[i]>10000)
		{
		  //  cout << "Maximum at 0 ns, or signal too large ==> something must be wrong. The program skips this event" << endl;
		  //  continue;
		}
		  
	      //  if (ampl[0]<10) goto NoRoot;
	      if (!FWF2)
		{
		  t_bck[i] = 3; // (TMax[i]-10>0) ? TMax[i]-10 : 0. ;
		  t_pul[i] = timeS[np-2]-5; // (TMax[i]+10 < timeS[np-2]) ? TMax[i]+10 : timeS[np-2] ;
		}
	      else
		{
		  t_bck[i] = 5. ;
		  t_pul[i] = 30.;
		}

	      //	      t_bck[i] = 8. ;
	      // t_pul[i] = 12.;


	      //	      cout << "check 1" << endl;
	      
	      if (FreqCut[i] == 0)
		{
		  //	  cout << " Check " << i << endl;
		  mobileAVG(samples[i],&amp[i][0],nmedia[i],&m_amp[i][0]);
		  // cout << " Check " << i <<  " samples " << samples[i] << " " << m_amp[i][100] << endl; 
		}
	      else
		{
		  if (i==0)
		    for(int l=0;l<samples[i];l++)		      
		      freq[l]=(1.e9/DT)*(1./samples[i])*l; //http://it.mathworks.com/help/matlab/math/fast-fourier-transform-fft.html			  
		  
		  FFTrealtocomplex(samples[i],&amp[i][0],&FFT_real[i][0],&FFT_comp[i][0],&FFT_abs[i][0]);
		  for(int l=0;l<samples[i];l++)
		    { 
		      if (freq[l]< FreqCut[i])
			{
			  FFT_real_cut[i][l] = FFT_real[i][l];
			  FFT_comp_cut[i][l] = FFT_comp[i][l];
			}
		      FFTcomplextoreal(samples[i],&FFT_real_cut[i][0],&FFT_comp_cut[i][0],&m_amp[i][0]);
		    }
		}


	
	      for(j=0;j<samples[i];j++)
		{
		  if(j%CAMPFACT==0) // option to undersample the data by CAMPFACT
		    {
		      timerec[nprec]=timeS[j]; // fills the time axis
		      amprec[i][nprec]=amp[i][j]; // fills the amplitude
		      m_amprec[i][nprec]=m_amp[i][j]; //fill the smoothed amplitude
		      if(polarity[i]<0)
			{
			  amprec[i][nprec]=-amprec[i][nprec];
			  m_amprec[i][nprec]=-m_amprec[i][nprec];
			}
		      // cout << nprec << "\t" << timerec[nprec] << "\t" << amprec[i][nprec] << endl;
		      nprec++;
		    }
		}
	      
	      samplesrec[i]=nprec;
		
	      // cout << "check 2 " << i << endl;

	      // compute the derivative on the signal shape
	      Derivative( DT, samples[i],&m_amprec[i][0],&der_amp[i][0]);

	      //   cout << " Check " << i <<  " samples " << samples[i] << " " << m_amprec[i][100] << endl; 
	      //calculate the baseline rms and amplitude after the signal, the from 10 to 15 ns after the max
	     
	      //	      if (!FWF2) Background(&m_amprec[i][0], samples[i]-2./DT,samples[i]-1./DT,&bck[i],&max_bck_after[i],&rms_bck_after[i]);
	      Background(&m_amprec[i][0], np_Max-5./DT,np_Max-5,&bck[i],&max_bck_after[i],&rms_bck_after[i]);

	      //calculate the baseline rms and amplitude before the signal, from  0 to 5 ns 
	      //   cout << " Check " << i <<  "bck " << bck[i] << " " <<  max_bck_after[i]<< endl; 
	      //	      Background(&m_amprec[i][0], NMax[i]-10./DT, NMax[i]-5./DT,&bck[i],&max_bck[i],&rms_bck[i]);
	     Background(&m_amprec[i][0], 0./DT, 5./DT,&bck[i],&max_bck_before[i],&rms_bck_before[i]);

	     //	      std::cout << t_bck[i]-5 << " " << t_bck[i] << std::endl;
	     //  cout << " Check " << i <<  "bck " << bck[i] << " " <<  max_bck_before[i]<< endl; 
	     
	      
	     Amplitudes(samples[i],&m_amprec[i][0],timeS,bck[i],t_bck[i],t_pul[i],&ampl[i],&t_max[i],&ampl_chi2[i]);
	      


	      //  cout << " Check " << i <<  "bck " << bck[i] << "  ampl " <<  ampl[i]<< endl; 
	      //Amplitudes_NF(&m_amprec[i][0],timeS,bck[i],NMax[i],&ampl[i],&t_max[i]);

	      
	      //  if (i==0 &&( ampl[0] < 200 || rms_bck_before[0]>10)  && toffee =="dig") goto NoRoot; // for FAST2_Dig
	      // if (i==1 && rms_bck_after[1]>8  && toffee =="dig") goto NoRoot; // for FAST2_Dig
	      //     if (i==0 && ampl[0] < 10) goto NoRoot;
	      //  if (i==0 && rms_bck_after[0]>10) goto NoRoot;
	      
	      	     
	      //	      Charges(&m_amprec[i][0],DT,NMax[i]-5./DT, NMax[i]+5./DT,&area[i],&totcha[i],50);

	      // Time of signal at 3 and 5 sigma noise
	      
	      Levtime_N(&m_amprec[i][0],DT,NMax[i],bck[i],3.*rms_bck_before[i],&t_rms3[i],Interpolation, &Nth);
	      Levtime_N(&m_amprec[i][0],DT,NMax[i],bck[i],5.*rms_bck_before[i],&t_rms5[i],Interpolation, &Nth);		      

	      // Count the levels below the signal amplitude
	      
	      NArr = 0;		  
	      for (int kl = 0; kl<NArrayL; kl++)
		{
		  NTValuesL[kl] = 0;
		  TValuesL[kl] = 0;
		  //		  if ( ThresholdL[kl]<(m_amprec[i][ (int) NMax[i]]) ) NArr++;
		  if ( ThresholdL[kl]<ampl[i] ) NArr++;
		}
	      
	      TrailNArr = 0;		  
	      for (int kl = 0; kl<TrailNArrayL; kl++)
		{
		  TrailTValuesL[kl] = 0;
		  //  if ( TrailThresholdL[kl]<(m_amprec[i][ (int) NMax[i]]) ) TrailNArr++;
		  if ( TrailThresholdL[kl]<ampl[i] )  TrailNArr++;
		}

	      for (int kl = 0;  kl<NArrayL; kl++)
		{
		  TValuesL[kl]=0;
		  NTValuesL[kl]=0;
		}

	      for (int kl = 0; kl<TrailNArrayL; kl++)
		{
		  TrailTValuesL[kl]=0;
		  TrailNTValuesL[kl]=0;
		}

	      
	       if (NArr>0)      Levtime_N_array( &NTValuesL[0], &TValuesL[0],  NArr, &ThresholdL[0], &m_amprec[i][0],DT,NMax[i],bck[i],Interpolation);
	       if (TrailNArr>0) TrailLevtime_N_array(&TrailTValuesL[0],  TrailNArr, &TrailThresholdL[0], samples[i], &m_amprec[i][0],DT,NMax[i],bck[i],Interpolation);

	      
	      // Calculate Time of Constant Fraction values

	       if (ampl[i]>10) ConstFractime_N_array( &NTValuesF[0], &TValuesF[0],  NArrayF, &ThresholdF[0],  &m_amprec[i][0],DT,bck[i],ampl[i],NMax[i],Interpolation);
	       if (ampl[i]>10) TrailConstFractime_N_array(&TrailNTValuesF[0], &TrailTValuesF[0],  TrailNArrayF, &TrailThresholdF[0], samples[i],  &m_amprec[i][0],DT,bck[i],ampl[i],NMax[i],Interpolation);
        
		   // Fill Level Values
		   
	       t_level10[i] = TValuesL[0];
	       t_level15[i] = TValuesL[1];
	       t_level20[i] = TValuesL[2];
	       t_level30[i] = TValuesL[3];
	       t_level40[i] = TValuesL[4];
	       t_level50[i] = TValuesL[5];
	       t_level60[i] = TValuesL[6];
	       t_level80[i] = TValuesL[7];
	       t_level100[i] = TValuesL[8];
	       t_level200[i] = TValuesL[9];
	       t_level300[i] = TValuesL[10];

	       //	       if (i==0 && (t_level100[0]< 8 || t_level100[0]> 15) && toffee =="dig") goto NoRoot; // for FAST2_Dig
	       
	       trail_t_level10[i] = TrailTValuesL[0];
	       trail_t_level20[i] = TrailTValuesL[1];
	       trail_t_level30[i] = TrailTValuesL[2];
	       trail_t_level40[i] = TrailTValuesL[3];
	       trail_t_level50[i] = TrailTValuesL[4];
	       trail_t_level60[i] = TrailTValuesL[5];
	       trail_t_level80[i] = TrailTValuesL[6];
	       trail_t_level100[i] = TrailTValuesL[7];
	       trail_t_level120[i] = TrailTValuesL[8];
	       trail_t_level140[i] = TrailTValuesL[9];
	       trail_t_level160[i] = TrailTValuesL[10];
	       trail_t_level180[i] = TrailTValuesL[11];
	       
	       //	      cout << t_level30[i] << " " <<  trail_t_level30[i] << endl;
	       
	       // Fill Constant Fraction Values
	       
	       Ncfd05[i] = NTValuesF[0];
	       Ncfd10[i] = NTValuesF[1];
	       Ncfd15[i] = NTValuesF[2];
	       Ncfd20[i] = NTValuesF[3];
	       Ncfd25[i] = NTValuesF[4];
	       Ncfd30[i] = NTValuesF[5];
	       Ncfd40[i] = NTValuesF[6];
	       Ncfd50[i] = NTValuesF[7];
	       Ncfd60[i] = NTValuesF[8];
	       Ncfd70[i] = NTValuesF[9];
	       Ncfd80[i] = NTValuesF[10];
	       Ncfd90[i] = NTValuesF[11];
	       
	       
	       cfd05[i] = TValuesF[0];
	       cfd10[i] = TValuesF[1];
	       cfd15[i] = TValuesF[2];
	       cfd20[i] = TValuesF[3];
	       cfd25[i] = TValuesF[4];
	       cfd30[i] = TValuesF[5];
	       cfd40[i] = TValuesF[6];
	       cfd50[i] = TValuesF[7];
	       cfd60[i] = TValuesF[8];
	       cfd70[i] = TValuesF[9];
	       cfd80[i] = TValuesF[10];
	       cfd90[i] = TValuesF[11];
	       
	       //	      trail_cfd10[i] = TrailTValuesF[0];
	       trail_cfd30[i] = TrailTValuesF[1];
	       trail_cfd50[i] = TrailTValuesF[2];
	       trail_cfd70[i] = TrailTValuesF[3];
	       trail_cfd90[i] = TrailTValuesF[4];
	       
	       trail_Ncfd10[i] = TrailNTValuesF[0];
	       trail_Ncfd30[i] = TrailNTValuesF[1];
	       trail_Ncfd50[i] = TrailNTValuesF[2];
	       trail_Ncfd70[i] = TrailNTValuesF[3];
	       trail_Ncfd90[i] = TrailNTValuesF[4];

	       Centroid_N( DT, Ncfd10[i], trail_Ncfd30[i] , &m_amprec[i][0], &t_cent[i]);
	       
	       // Compute area between the starting 10% and falling 30% of the signal
	       
	       Charges(&m_amprec[i][0],DT,Ncfd10[i],trail_Ncfd30[i],&area[i],&totcha[i],50);
	       
	       //	      cout << "check 3 " << i << endl;
	       // 
	       if (( (cfd10[i]/DT) <samples[i] &&   (cfd20[i]/DT) < samples[i] && (cfd30[i]/DT)<samples[i]
		     &&  (cfd70[i]/DT) <samples[i]&&   (cfd80[i]/DT) <samples[i]) &&
		   ( (cfd10[i]/DT) >0 &&   (cfd20[i]/DT) >0 && (cfd30[i]/DT)>0
		     &&  (cfd70[i]/DT) >0 &&   (cfd80[i]/DT) >0))
		 {
		   
		   if ((cfd30[i]- cfd10[i]) != 0)
		     dVdt1030[i] =( m_amprec[i][ int (cfd30[i]/DT)]- m_amprec[i][int (cfd10[i]/DT)])/(cfd30[i]- cfd10[i]);
		   
		   if  ((cfd70[i]- cfd30[i]) != 0)
		     dVdt3070[i] =( m_amprec[i][ int (cfd70[i]/DT)]- m_amprec[i][int (cfd30[i]/DT)])/(cfd70[i]- cfd30[i]);
		   
		   if  ((cfd80[i]- cfd20[i]) != 0)
		     dVdt2080[i] =( m_amprec[i][ int (cfd80[i]/DT)]- m_amprec[i][int (cfd20[i]/DT)])/(cfd80[i]- cfd20[i]);
		 }
	       // cout << dVdt3070[i] << endl;
	       
	       if(doFit)
		 {
		   
		   riselinfit(samples[i],&m_amprec[i][0],timeS, bck[i], t_level20[i],t_level40[i], &rise_lin0[i], &rise_lin1[i],&rise_lin_chi2[i]) ;
		   if (rise_lin1[i] != 0) t_zero[i]= (-rise_lin0[i])/rise_lin1[i];
		   riseexpfit(samples[i],&m_amprec[i][0],timeS, bck[i], cfd20[i]-5,cfd20[i], &rise_exp0[i], &rise_exp1[i],&rise_exp_chi2[i]) ;
		 }
	       
	       // } // toffee=dig
	       
	       nprec=0;
	      //
	      maxD[i] = 0;
	      for(int l=0;l<samples[i];l++)
		{
		  m_amprec[i][l] = m_amprec[i][l]-bck[i];
		  if(l%CAMPFACT==0)
		    {
		      d_amprec[i][nprec]=der_amp[i][l];
		      if (fabs(d_amprec[i][nprec])>fabs(maxD[i])) 		       
			maxD[i] = fabs(d_amprec[i][nprec]);
			  
		      nprec++;
		    }
		}
	      //	      cout << " # of chr = " << i << endl;  
	    } // end loop on channels
	  //	  f->Close();

	  
	  if(running==0 && emptycount<100)
	    {
	      running=1;
	      continue;
	    }
	  if(running==0 && emptycount>=100)
	    {
	      break;
	    }
	  
	  
	  
	  //  if (ampl[0] > 10 || ampl[1] > 10) 
	  // {

	  
	  //	  w->Fill();
	  OutTree->Fill();

	NoRoot:
	  continue;
	  // }
	  
	}  // while{running}

 Write:

      continue;
    

    }

   OutputFile->cd();
   // w->Write();
   OutTree->Write();
   OutputFile->Close();

   return(EXIT_SUCCESS);
}


