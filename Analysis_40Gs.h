#ifndef WAVE_ANALYZER_H
#define WAVE_ANALYZER_H

#define NCHRO 6
#define MAX_CAMP 20000
#define CAMPFACT 1
#define MAX_CAMPREC MAX_CAMP/CAMPFACT+100

Int_t osc[NCHRO];
Int_t chan[NCHRO];
Double_t bck[NCHRO];
Double_t t_bck[NCHRO];
Double_t t_pul[NCHRO];
Int_t camp[NCHRO];
Int_t camprec[NCHRO];
Int_t samples[NCHRO];
Int_t samplesrec[NCHRO];
Double_t FreqCut[NCHRO];
Double_t Toffset[NCHRO];


Double_t DT;
Double_t timeS[MAX_CAMP]; 
Double_t Itot[NCHRO][MAX_CAMP];
Double_t amp[NCHRO][MAX_CAMP];
Double_t m_amp[NCHRO][MAX_CAMP];
Double_t der_amp[NCHRO][MAX_CAMP];
Double_t amprec[NCHRO][MAX_CAMPREC];
Double_t m_amprec[NCHRO][MAX_CAMPREC];
Double_t d_amprec[NCHRO][MAX_CAMP];
Double_t FFT_real[NCHRO][MAX_CAMP];
Double_t FFT_comp[NCHRO][MAX_CAMP];
Double_t FFT_abs[NCHRO][MAX_CAMP];

Double_t FFT_real_cut[NCHRO][MAX_CAMP];
Double_t FFT_comp_cut[NCHRO][MAX_CAMP];
Double_t FFT_abs_cut[NCHRO][MAX_CAMP];

Double_t freq[MAX_CAMPREC];
Double_t timerec[MAX_CAMPREC];

Double_t XCross;
Double_t YCross;

  //Analysis variables

Double_t day;
Double_t month;
Double_t year;
Double_t hour;
Double_t minute;

Double_t MaxDeriv;
Double_t sec[NCHRO];
Double_t ampl[NCHRO];
Double_t area[NCHRO];
Double_t dVdt1030[NCHRO];
Double_t dVdt3070[NCHRO];
Double_t dVdt2080[NCHRO];

  Double_t totcha[NCHRO];
  Double_t maxS[NCHRO];
  Double_t maxD[NCHRO];
  Double_t max_bck_before[NCHRO];
  Double_t rms_bck_before[NCHRO];
  Double_t max_bck_after[NCHRO];
  Double_t rms_bck_after[NCHRO];
  Double_t charge[NCHRO];
  Int_t strip[NCHRO];
  Double_t FEcharge[NCHRO];
  Double_t FEchargeerr[NCHRO];
  Double_t totcharge[NCHRO];
  Double_t p_charge;
  Double_t p_chargeerr;
  Int_t sat[NCHRO];
  Double_t maxamp;
  Int_t maxstrip;
  Double_t eff_dis08[NCHRO];
  Double_t eff_dis10[NCHRO];
  Double_t eff_dis15[NCHRO];
  Double_t eff_dis20[NCHRO];
  Double_t eff_dis30[NCHRO];
  Double_t eff_dis40[NCHRO];
  Double_t eff_dis50[NCHRO];
  Double_t eff_dis60[NCHRO];
  Double_t eff_dis70[NCHRO];
Double_t eff_dis100[NCHRO];
  Double_t eff_dis150[NCHRO];
  Double_t eff_dis200[NCHRO];
  Double_t eff_dis250[NCHRO];
  Double_t eff_dis300[NCHRO];
  Double_t eff_dis400[NCHRO];
  Double_t eff_dis500[NCHRO];
  Double_t eff_dis800[NCHRO];
  Double_t eff_dis1000[NCHRO];
  Double_t eff_dis2000[NCHRO];
  Double_t eff_dis5000[NCHRO];
  Double_t eff_dis10000[NCHRO];
  Double_t t_level10[NCHRO];
  Double_t t_level15[NCHRO];
  Double_t t_level18[NCHRO];
  Double_t t_level20[NCHRO];
  Double_t t_level30[NCHRO];
  Double_t t_level40[NCHRO];
  Double_t t_level50[NCHRO];
  Double_t t_level60[NCHRO];
  Double_t t_level80[NCHRO];
  Double_t t_level100[NCHRO];
  Double_t t_level200[NCHRO];
  Double_t t_level300[NCHRO];

 
  Double_t trail_t_level10[NCHRO];
  Double_t trail_t_level15[NCHRO];
  Double_t trail_t_level18[NCHRO];
  Double_t trail_t_level20[NCHRO];
  Double_t trail_t_level30[NCHRO];
  Double_t trail_t_level40[NCHRO];
  Double_t trail_t_level50[NCHRO];
  Double_t trail_t_level60[NCHRO];
  Double_t trail_t_level80[NCHRO];
  Double_t trail_t_level100[NCHRO];
  Double_t trail_t_level120[NCHRO];
  Double_t trail_t_level140[NCHRO];
 Double_t trail_t_level160[NCHRO];
  Double_t trail_t_level180[NCHRO];
  Double_t trail_t_level200[NCHRO];
  Double_t trail_t_level300[NCHRO];

  Double_t cfd05[NCHRO];
  Double_t cfd10[NCHRO];
  Double_t cfd15[NCHRO];
  Double_t cfd20[NCHRO];
  Double_t cfd25[NCHRO];
  Double_t cfd30[NCHRO];
  Double_t cfd40[NCHRO];
  Double_t cfd50[NCHRO];
  Double_t cfd60[NCHRO];
  Double_t cfd70[NCHRO];
  Double_t cfd80[NCHRO];
  Double_t cfd90[NCHRO];

  Int_t Ncfd05[NCHRO];
  Int_t Ncfd10[NCHRO];
  Int_t Ncfd15[NCHRO];
  Int_t Ncfd20[NCHRO];
  Int_t Ncfd25[NCHRO];
  Int_t Ncfd30[NCHRO];
  Int_t Ncfd40[NCHRO];
  Int_t Ncfd50[NCHRO];
  Int_t Ncfd60[NCHRO];
  Int_t Ncfd70[NCHRO];
  Int_t Ncfd80[NCHRO];
  Int_t Ncfd90[NCHRO];

  Int_t trail_Ncfd10[NCHRO];
  Int_t trail_Ncfd30[NCHRO];
  Int_t trail_Ncfd50[NCHRO];
  Int_t trail_Ncfd70[NCHRO];
  Int_t trail_Ncfd90[NCHRO];

  Double_t dist[NCHRO];
  Int_t ind[NCHRO];
  Int_t XPad[NCHRO];
  Int_t YPad[NCHRO];
  Double_t XDist[NCHRO];
  Double_t YDist[NCHRO];

  Double_t trail_cfd30[NCHRO];
  Double_t trail_cfd50[NCHRO];
  Double_t trail_cfd70[NCHRO];
  Double_t trail_cfd90[NCHRO];
 


  Int_t Nt_level30[NCHRO];
  Int_t Ntrail_t_level30[NCHRO];
  Double_t t_cent30[NCHRO];
  Double_t t_centMax[NCHRO];
  Double_t t_cent[NCHRO];

  Double_t t_rms3[NCHRO];
  Double_t t_rms5[NCHRO];
  Double_t t_max[NCHRO];
Double_t t_zero[NCHRO];
Double_t fitchi[NCHRO];
Double_t ampl_chi2[NCHRO];
Double_t rise_lin0[NCHRO];
Double_t rise_lin1[NCHRO];
Double_t rise_lin_chi2[NCHRO];
  Double_t rise_exp0[NCHRO];
  Double_t rise_exp1[NCHRO];
Double_t rise_exp_chi2[NCHRO];
Double_t AMax[NCHRO];
Double_t TMax[NCHRO];
Double_t NMax[NCHRO];
int polarity[NCHRO];


  // Analysis branches
Double_t WF2gain;
Int_t Wafer;
Int_t LGain;
Int_t FType;
Int_t DMetal;
Int_t DPitch;

Double_t WF2xpos;
Double_t WF2angle;
Double_t TotCharge;
Int_t MaxEvt;
Double_t XCord;
Double_t YCord;
Double_t XAmpl;
Double_t YAmpl;
Double_t XTime;
Double_t YTime;
Double_t XCirc;
Double_t YCirc;
Double_t XArea;
Double_t YArea;


Int_t CS;
Int_t eff;
Double_t noise=3.;
Double_t bar0;
Double_t bar1;
Double_t A0;
Double_t A1;
Double_t delta0;
Double_t delta1;
Double_t chinorm0;
Double_t chinorm1;


  // Double_t gum_amp[MAX_CAMP];


Int_t nchro; // = NCHRO;  //Maximum number of Read out channels to be present in the rootple
//Int_t nchActive = NCHRO;  //Number of active Read out channels, to be filled written in the InputFolder
Int_t nchro_active;
bool doFFT = false;  
bool showFFT = false;

std::string pip;
std::string toffee;
std::string InputNAME[100];
std::string InputNA;
int i;
int j;
int k;
int np;
int nprec;


Double_t hv;
Double_t hvcorr;
Double_t temp;
Double_t pres;
Double_t TTCfd30;
Double_t TTCfd50;
Double_t TCfd20;
Double_t TCfd30;
Double_t TCfd40;
Double_t TCfd50;
Double_t TTmax;

Int_t nrun=0;
std::string String_nrun;
Int_t ntrig=0;
Int_t ds=0;
Int_t event=0;
Int_t CAmax=0;


TFile *OutputFile;
TTree *OutTree;

std::ostringstream leaf(std::ostringstream::out);
std::ostringstream leafl(std::ostringstream::out);

  // Program variables
std::ostringstream datafile;

  float Tpeak = 0;
  float BallDef = 0;

int nmedia[NCHRO];
int nme;

  double tempsum;
  int running=1;
  int emptycount=0;
  double time0temp;
  int segsize;
  char number[5];
  char date[12];
  char timeofday[8];
  double FFT_ratio;



   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        i_timestamp;
   Double_t        i_current;
   std::vector<double>  *w1;
   std::vector<double>  *t1;
   std::vector<double>  *w2;
std::vector<double>  *t2;
std::vector<double>  *w3;
std::vector<double>  *t3;

std::vector<double>  *w4;
std::vector<double>  *t4;

TBranch        *b_i_timestamp;   //!
TBranch        *b_i_current;   //!
TBranch        *b_w1;   //!
TBranch        *b_t1;   //!
TBranch        *b_w2;   //!
TBranch        *b_t2;   //!
TBranch        *b_w3;   //!
TBranch        *b_t3;   //!
// TBranch        *b_w4;   //!
// TBranch        *b_t4;   //!


//Convert little endian int to big endian int
unsigned int intToBigEndianess(const unsigned int x)
{
  return  ( x >> 24 ) |  // Move first byte to the end,
    ( ( x << 8 ) & 0x00FF0000 ) | // move 2nd byte to 3rd,
    ( ( x >> 8 ) & 0x0000FF00 ) | // move 3rd byte to 2nd,
    ( x << 24 ); // move last byte to start.
}


//returns the float value corresponding to a given bit represention of a scalar int value or vector of int values
float intBitsToFloat(const int x)
{
  union 
  {
    float f;  // assuming 32-bit IEEE 754 single-precision
    int i;    // assuming 32-bit 2's complement int
  } 
  u;
  
  u.i = x;
  return u.f;
}


//Always run background function before doing further analyses

//Compute average amplitude for background
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, max background amplitude)
void Background(Double_t amp[],Int_t Nt0, Int_t Nt1, Double_t *bck, Double_t *maxbck, Double_t *rmsbck)
{

  //  double min=1000000;
  double max=-1000000;
  float sumamp=0;
  float sumampsq=0;
  int pointsavg=0;
  float meanbck;

  // TH1D hlite("hlite","hlite",camp,min,max);
  *bck=0;
  //  std::cout << Nt0 << " " << Nt1 << std::endl;
  for(int j=Nt0;j<Nt1;j++)
    {
      pointsavg++;
      sumamp += amp[j];
      sumampsq += amp[j]*amp[j];
    }
  meanbck=sumamp/float(pointsavg);
  *rmsbck=sqrt((1./float(pointsavg))*(sumampsq-float(pointsavg)*meanbck*meanbck));
  
  //*bck=hlite.GetBinCenter(hlite.GetMaximumBin());
  *bck=meanbck;

  for(int j=Nt0;j<Nt1;j++)
    {
      if(std::fabs(amp[j]-meanbck)>max)
	max=amp[j];
    }

  *maxbck=max;

  return;
}


//Inversione degli impulsi negativi
void Inversion(Int_t camp, Double_t neg_amp[])
{
  for(int j=0;j<camp;j++)
	{
	  neg_amp[j]=-neg_amp[j];
        }
}


//DEVELOPEMENT
// Compute maximum amplitude and peak time for waveform using a gaussian fit to the peak
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse maximum amplitude (averaged), pulse maximum amplitude)
void Amplitudes(Int_t camp, Double_t amp[], Double_t timeS[],  Double_t bck, Double_t t_bck,Double_t t_pul, Double_t *pamp, Double_t *tamp, Double_t *chi)
{
  
  Int_t np=0;
  //  int cont=0;
  int j_max= 0;
  int j_start = 0;
  int n_j = 0;
  // Limit of parabola fitting around j_max
  float AmpFrac = 0.8;
  // double tempmean;
  // double tempsig;
  double_t   max=-100000;
  double_t   tmax=-100000;

  //Maximum calculation
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  np++;
	  if(amp[j]-bck>max)
	    {
	      max=amp[j]-bck;
	      tmax=(timeS[j+1]+timeS[j])/2.;
	      j_max=j;
	    }
	}
    }

  //Fit amplitude calculation
  // Find limits at AmpFrac
  
  //from the maximum backwards 
  for(int j=j_max;j>0;j--)
    {
      if(amp[j]-bck<AmpFrac*(max))
	{
	  j_start=j;
	  break;
	}
    }

  //from the minimum limit, up the maximum, and then down till the limit on the falling edge  
  for(int j=j_start+1;j<camp-1;j++)
    {
      n_j++;
      if(amp[j]-bck<AmpFrac*(max))
	{
	  break;
	}
    }

  
  //  std::cout << j_start << " " <<j_max << " Number of points in the fit= " << n_j << std::endl; 
  
  double temp0;
  double temp1;
  double temp2;
  TF1 *f1 = new TF1("gausfit","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",timeS[j_start-1],timeS[j_start+n_j+1]);
  // TGraph g(n_j);
  TGraph *g = new TGraph();
  g->Set(n_j);
  for(int j=0;j<n_j;j++)
    {
      g->SetPoint(j,timeS[j+j_start],amp[j+j_start]-bck);
    }

  temp0=max;
  temp1=tmax;
  temp2=(timeS[j_start]-timeS[j_start+n_j])/2.;
  f1->SetParameters(temp0,temp1,temp2);
  g->Fit("gausfit","QN","",timeS[j_start]-1,timeS[j_start+n_j-1]+1);
  //  g.Fit("gausfit","QN","",timeS[j_start]-1,timeS[j_start+n_j-1]+1);
  *pamp=f1->GetParameter(0);
  *tamp=f1->GetParameter(1);

  *chi=f1->GetChisquare()/f1->GetNDF();
  g->Delete();
  f1->Delete();

  return;
}

//Amplitude without the Fit
void Amplitudes_NF( Double_t amp[], Double_t timeS[],  Double_t bck, Int_t NMax, Double_t *max,Double_t *tmax)
{
  

  *max=-100000;
  *tmax=-100000;
  *tmax=(timeS[NMax-1]*amp[NMax-1]+timeS[NMax]*amp[NMax]+timeS[NMax+1]*amp[NMax+1])/(amp[NMax-1]+amp[NMax]+amp[NMax+1]);
  *max=(amp[NMax-1]+amp[NMax]+amp[NMax+1]-3*bck)/3;
 
  return;
}



//DEVELOPEMENT
// Compute maximum, amplitude and peak time for waveform
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse maximum amplitude (averaged), pulse maximum amplitude)
void Neg_Amplitudes(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck, Double_t *pamp, Double_t *pmax, Double_t *tmax)
{
  *pamp = 0;
  Int_t np=0;
  // int cont=0;
  double maxtemp=-100000;
  //  double tempmean;
  //  double tempsig;
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  np++;
	  if(-(amp[j]-bck)>maxtemp)
	    {
	      maxtemp=-(amp[j]-bck);
	      *tmax=(timeS[j+1]+timeS[j])/2.;
	    }
	}
    }
  *pmax=maxtemp;

  /*
  TGraph *g = new TGraph(np);
  TF1 *f= new TF1("gauss","[0]*exp(-((x-[1])*(x-[1]))/(2.*[2]*[2]))",t_bck,t_pul);
  f->SetParameters(*pmax,(t_pul+t_bck)/2.,(t_pul-t_bck)/2.);
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  cont++;
	  g->SetPoint(cont,timeS[j],amp[j]-bck);
	}
    }
  g->Fit("gauss","QN","",t_bck,t_pul);
  
  tempmean=f->GetParameter(1);
  tempsig=f->GetParameter(2);
  g->Fit("gauss","QN","",tempmean-tempsig,tempmean+tempsig);
  
  *pamp=f->GetParameter(0);
  */
  return;
}

// Calibration using curve [0](1-e^(x/[1]))
// Set the function with the proper front end each time
// (Parameter zero, parameter 1, amplitude of amplified pulse, corresponding charge/amplitude of input pulse, error on input charge/amplitude, error on output amplitude - must be known)
void Calibration(Double_t cal_low_par0,Double_t cal_low_par1,Double_t cal_high_par0,Double_t cal_high_par1, Double_t max,Double_t *charge, Double_t *chargeerr, Double_t amperr)
{
  TF1 f_low("cal_low","[0]*log(1.-x/[1])",0.,cal_low_par1-0.1);
  TF1 f_high("cal_high","[0]*log(1.-pow((x/[1]),2))",0.,cal_high_par1-0.1);
  f_low.SetParameters(cal_low_par0,cal_low_par1);
  f_high.SetParameters(cal_high_par0,cal_high_par1);
  if(max<444)
    {
      *charge=f_low(max);
      *chargeerr=(-cal_low_par0/(cal_low_par1-(max)))*amperr;
    }
  
  if(max>=444 && max<cal_high_par1)
    {
      *charge=f_high(max);
      *chargeerr=(-cal_high_par0)*(1./(1.-pow(1.-max/cal_high_par1,2)))*2.*(1.-max/cal_high_par1)*amperr;
    }
  if(max>=cal_high_par1)
    {
      *charge=1200.;
      *chargeerr=1000.;
    }
  
  return;
}

// Calibration using curve [0](1-e^(x/[1]))
// Set the function with the proper front end each time
// (Parameter zero, parameter 1, amplitudesof amplified pulse, corresponding charge/amplitude of input pulse, error on input charge/amplitude, error on output amplitude - must be known)
void Calibration_ip(Double_t cal_A0,Double_t cal_x0,Double_t cal_B0,Double_t max,Double_t *charge, Double_t *chargeerr, Double_t amperr)
{
  TF1 f("cal_ip","[1]+pow([0]/(x-[2]),0.25)",-10.,900);
  f.SetParameters(cal_A0,cal_x0,cal_B0);
  
  *charge=f(max);
  *chargeerr=-0.25*(pow(cal_A0/(max-cal_B0),1.25))*amperr/cal_A0;
    
   
  if(max>=cal_B0)
    {
      *charge=1500.;
      *chargeerr=1000.;
    }
  
  return;
}


// Integrates the pulse to give charge - readout impedance must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse charge, total charge, readout resistance)
void Charges(Double_t amp[], Double_t Dt, Int_t t_start, Int_t t_stop, Double_t *area, Double_t *charge, double R)
{
  double tempminch=0;
  double bck=0;
  for(int j=0;j<10;j++)
    bck += amp[j]/10;
    
  
  //  std::cout << "start - stop " << t_start << " - " << t_stop << std::endl;
  
  for(int j=t_start;j<t_stop;j++)
    {
      // std::cout << " j = " << j <<   " total= " <<  tempminch <<  " amp " << amp[j] << std::endl;
      tempminch += amp[j]-bck;
    }
      *charge=tempminch*Dt/R;
      *area=tempminch*Dt;
  return;
}


//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, threshold to comupute efficiency, efficiency)
void Efficiencies(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck,Double_t threshold, Double_t *eff)
{
  *eff=0;
  for(int j=0;j<camp;j++)
    {
      if(timeS[j]<t_bck)
	continue;
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  if((amp[j]-bck-threshold)>0 && threshold>0)
	    {
	      *eff=1;
	      break;
	    }
	  if((amp[j]-bck-threshold)<0 && threshold<0)
	    {
	      *eff=1;
	      break;
	    }
	}
    }
  return;
}

// Time at which a fixed threshold is passed. Makes a linear fit on 6 points and inverts the equation
// time  = (threshold-q)/q
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, threshold to measure time, time at which the threshold is passed)
void Thrtime_old(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck,Double_t threshold, Double_t *thrtime)
{
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=0;j<camp;j++)
    {
      if(timeS[j]<t_bck)
	continue;
      if(timeS[j]>t_pul)
	continue;
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  //	  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
	  if((amp[j]-bck-threshold)>0 && threshold>0)
	    {
	      //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	      //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	      //std::cout << "Der = " << Der << std::endl;
	      //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //   std::cout << "th1 = " << *thrtime << std::endl;
	      for ( int kk = 0;kk<npoints;kk++)
		{
		  g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
		  //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
		  //  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
		}

	      break;

	    }
	  //	  if((amp[j]-bck-threshold)<0 && threshold<0)
	  // {
	  //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
	}
    }
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
    //    Der= lin->GetParameter(1);
  if (lin->GetParameter(1)>0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
  // std::cout << "th1bis = " << *thrtime << std::endl;
  g->Delete();
  lin->Delete();
  return;
}

// Time at which a threshold relative to the amplitude of the pulse is passed. Amplitude of the pulse must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, relative threshold to measure time, time at which the threshold is passed)


void ConstFractime(Double_t amp[], Double_t timeS[], Double_t bck,Double_t threshold, Double_t max, Int_t NMax,Double_t *thrtime)
{
  
  
  *thrtime=0;

  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  //  *thrtime=0;
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  //  for(int j=0;j<camp;j++)
    for(int j=NMax;j>1;j--)
    {
      //      if(timeS[j]>t_max || timeS[j]<t_bck )
      //	continue;
      // std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  //	  for ( int kk = 0;kk<npoints;kk++)
	  // {
	  //   g->SetPoint(kk,(amp[j-npoints/2+kk]-bck)/max,timeS[j-npoints/2+kk]);
	      //   std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << (amp[j-npoints/2+kk]-bck)/max << " " << threshold << std::endl;
	  // }
	  C2 = (amp[j+1]-bck)/max;
	  C1 = (amp[j]-bck)/max;
	  //  std::cout << "time = " << timeS[j]  <<  " " << timeS[j+1]<<  " C1 = " << C1 << " C2 = " << C2 <<std::endl;
	  if ( (C2-C1) != 0) TCrossing =  timeS[j]+(timeS[j+1]-timeS[j])*(threshold-C1)/(C2-C1);
	  break;
	  
	}
	
    }
  //  std::cout << "th1 = " << *thrtime << std::endl;
    //  TF1 *lin = new TF1("lin","[0]+[1]*x");
    // g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  
    //  if (lin->GetParameter(1) !=0)
    //  *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
    //  *thrtime= lin->GetParameter(0) +threshold*lin->GetParameter(1);
    *thrtime= TCrossing;
    //else
    // *thrtime = 0;
    
  //  std::cout << "th fit = " << *thrtime <<  " Th frac = " << TCrossing << " "  << " C1 = " << C1 << " C2 = " << C2 << std::endl;
    // g->Delete();
    // lin->Delete();
  return;
}


// Time at which a threshold relative to the amplitude of the pulse is passed. Amplitude of the pulse must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, relative threshold to measure time, time at which the threshold is passed)
void Zerofittime(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck,Double_t threshold1, Double_t threshold2, Double_t t_max,Double_t *thrtime)
{
  t_pul = 0;
  *thrtime=0;
  double N=0;
  double sumV=0;
  double sumT=0;
  double sumVT=0;
  double sumTT=0;
  for(int j=camp-1;j>=0;j--)
    {
      if(timeS[j]>t_max)
	continue;
      if(timeS[j]<t_bck)
	continue;
      if(amp[j]-bck>threshold1 && amp[j]-bck<threshold2 && timeS[j]<t_max && N<3)
	{
	  N=N+1;
	  sumV=sumV+amp[j]-bck;
	  sumT=sumT+timeS[j];
	  sumVT=sumVT+(amp[j]-bck)*timeS[j];
	  sumTT=sumTT+timeS[j]*timeS[j];
	}
    }
  if(N>0)
    {
      double Vm=(N*sumVT-sumV*sumT)/(sumTT-sumT*sumT);
      double V0=(sumV-Vm*sumT)/N;
      *thrtime=-V0/Vm;
    }
  if(N<2)
    *thrtime=-99;
  return;
}

// Time at which a fixed threshold is passed (descending)
//
void Trailtime(Int_t camp, Double_t amp[], Double_t timeS[], Int_t NMax, Double_t bck,Double_t threshold, Double_t *thrtime)
{
  //  std::cout << "NMAX = " << NMax << std::endl;
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j<camp;j++)
    {
     
      //  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	{
	  //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	  //std::cout << "Der = " << Der << std::endl;
	  //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //	       std::cout << "th1 = " << timeS[j] << std::endl;
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
	      //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
	      //	  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
	    }
	  
	  break;
	  
	}
      //	  if((amp[j]-bck-threshold)<0 && threshold<0)
      // {
      //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
    }
    
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  //    Der= lin->GetParameter(1);
  if (lin->GetParameter(1)<0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
  //  std::cout << "th1bis = " << *thrtime << std::endl;
  g->Delete();
  lin->Delete();
  return;
}


// Time at which a fixed threshold is passed (going to the maximum)
// The search starts from the maximum and goes backwards
void Thrtime( Double_t amp[], Double_t timeS[], Int_t NMax, Double_t bck,Double_t threshold, Double_t *thrtime)
{
  //  std::cout << "NMAX = " << NMax << std::endl;
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j>0;j--)
    {
     
      //std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	{
	  //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	  //std::cout << "Der = " << Der << std::endl;
	  //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //	       std::cout << "th1 = " << timeS[j] << std::endl;
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
	      //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
	      //	  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
	    }
	  
	  break;
	  
	}
      //	  if((amp[j]-bck-threshold)<0 && threshold<0)
      // {
      //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
    }
    
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  //    Der= lin->GetParameter(1);
  // if (lin->GetParameter(1)<0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
    //else
    // *thrtime = 0;
    //    std::cout << "th1bis = " << *thrtime << std::endl;
    g->Delete();
    lin->Delete();
  return;
}


// Time at which a constant fraction threshold is passed (descending)
//(record lenght, amplitude vector, time vector, average background amplitude, relative threshold to measure time, starting time for the search (example: time of the matimum), time of the trailing edge, isteresys with respect to the leading edge)
void TrailConstFractime(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t threshold, Double_t max, Int_t NMax,Double_t *trailtime)
{
  *trailtime=0;

  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  //  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j<camp;j++)
    {
      //	  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],(amp[j-npoints/2+kk]-bck)/max);
	      //		  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << (amp[j-npoints/2+kk]-bck)/max << " " << threshold << std::endl;
	    }
	  
	  break;
	  
	}
      
    }
    // std::cout << "th1 = " << *trailtime << std::endl;

    // TF1 *quad = new TF1("quad","[0]+[1]*x+[2]*x*x");
    TF1 *lin = new TF1("lin","[0]+[1]*x");
    g->Fit("lin","QN","goff");
    //g->Fit("lin");
    // g->Fit("quad");
  if (lin->GetParameter(1) !=0)
    *trailtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *trailtime = 0;
    
  //  std::cout << "th2 = " << *trailtime << std::endl;
  g->Delete();
  lin->Delete();

  
  return;
}


// Says if the event suffered from saturation of the waveform readout or not. Must be calibrated each time to evaluate the sensitivity of the routine
//(record lenght, amplitude vector, time vector,saturation flag)
void Saturation(Int_t camp, Double_t amp[], Double_t timeS[], Double_t bck, Int_t *satur)
{
  *satur=0;
  timeS = 0;
  int cont=0;
  int tempcont=0;
  double value=100000;
  for(int j=1;j<camp;j++)
    {
      if(fabs(amp[j]-amp[j-1])>0.0001)
	{
	  tempcont=0;
	}
      if(fabs(amp[j]-amp[j-1])<0.0001)
	{
	  tempcont++;
	  if(tempcont>cont)
	    {
	      cont=tempcont;
	      value=amp[j];
	    }
	}
    }
  if(fabs(value-bck)>1. && cont>int(camp/500)+5)
    *satur=1;
  return;
}

// Computed the charge centroid using the A0/cosh((x-x0)/delta0) function.
//(Charge, Error of charge, number of RO channels, distribution width, distribution width 2nd method, Distribution const, Distribution const 2nd method, Centroid , Centroid 2nd method , ChiSquare, ChiSquare 2nd method,pitch, Cluster size, threshold)
void Centroid(Double_t ch[], Double_t cherr[],Int_t nchro, Double_t *delta0, Double_t *delta1, Double_t *A0, Double_t *A1, Double_t *x0, Double_t *x1, Double_t *chinorm0, Double_t *chinorm1, Double_t pitch, Int_t CS, Double_t threshold)
{
  TF1 cint("cint",Form("[0]*atan(exp((x+%f-[1])/[2]))-[0]*atan(exp((x-%f-[1])/[2]))",pitch/2.,pitch/2.),-100.,100.);
  TF1 cdis("cdis","[0]*(1./cosh((x-[1])/[2]))",-100.,100.);
  TF1 cintf("cintf",Form("[0]*atan(exp((x+%f-[1])/[2]))-[0]*atan(exp((x-%f-[1])/[2]))",pitch/2.,pitch/2.),-100.,100.);
  TGraphErrors gbar(CS);
  double tmpmax=-100;
  double tmpbar=0;
  int k=0;
  *delta0=-1;
  *delta1=-1;
  *A0=-1;
  *A1=-1;
  *x0=-10;
  *x1=-10;
  *chinorm0=-1;
  *chinorm1=-1;
  cintf.FixParameter(2,3.9);
  if(CS>2)
    {
      for(int i=0;i<nchro;i++)
	{
	  if(ch[i]>tmpmax)
	    {
	      tmpmax=ch[i];
	      tmpbar=double(i)*pitch+pitch/2.;
	    }
	}
      
      cint.SetParameters(tmpmax,tmpbar,3.5);
      cint.SetParLimits(1,-pitch*5.,double(nchro)*pitch+pitch*5.);
      cintf.SetParameters(tmpmax,tmpbar,3.9);
      cintf.FixParameter(2,3.9);
      cintf.SetParLimits(1,-pitch*5.,double(nchro)*pitch+pitch*5.);

      for(int i=0;i<nchro;i++)
	{
	  if(ch[i]>threshold)
	    {
	      k++;
	      gbar.SetPoint(k,double(i)*pitch+pitch/2.,ch[i]);
	      gbar.SetPointError(k,0.,cherr[i]);
	    }
	}
      
      gbar.Fit("cint","QN","",0.,nchro*pitch+pitch/2.);
      *delta0=cint.GetParameter(2);
      *x0=cint.GetParameter(1);
      *A0=cint.GetParameter(0);
      *chinorm0=cint.GetChisquare()/cint.GetNDF();
      gbar.Fit("cintf","QN","",0.,nchro*pitch+pitch/2.);
      *delta1=cintf.GetParameter(2);
      *x1=cintf.GetParameter(1);
      *A1=cintf.GetParameter(0);
      *chinorm1=cintf.GetChisquare()/cintf.GetNDF();
    }
  if(CS==2)
    {
      *x0=0;
      float tempnum=0;
      float tempden=0;
      float bar0temp=0;
      float tempor=0;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      bar0temp=bar0temp+(float(j)*pitch+pitch/2.);
	    }
	}
      bar0temp=bar0temp/2.;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      tempnum=tempnum+ch[j]*(float(j)*pitch+pitch/2.-bar0temp);
	      tempden=tempden+ch[j];
	      //cout << FEcharge[j] << "\t";
	    }
	}
      tempor=tempnum/tempden;
      *x0=bar0temp+7.5*tempor;  //7.5 conversion factor from charge centroid to mm
    }
  if(CS==1)
    {
      *x0=0;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      *x0=*x0+(float(j)*pitch+pitch/2.);
	    }
	}
    }
  return;
}


//  Computes the slope of the rise of a signal as an exponential
void riseexpfit(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t t_low,  Double_t t_high, Double_t *rise_exp0, Double_t *rise_exp1, Double_t *chi)
{
  TGraph *g = new TGraph();
  int npoints=0;
  int nnpoints=0;
  int Nlow=0;

  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)  continue;
      if(timeS[i]>t_high) break;
      npoints++;
      if (npoints == 1) 
	{
	  Nlow = i;
	}
    }

  g->Set(npoints);
  nnpoints=0;
  for(int i=Nlow;i<Nlow+npoints-1;i++)
    {
      g->SetPoint(nnpoints,timeS[i]-t_low,amp[i]-bck);
	    //      g->SetPoint(nnpoints,timeS[i],amp[i]-bck);
	    //  std::cout << timeS[i]<< " " << amp[i]-bck << std::endl;
      nnpoints++;
    }
  
  TF1 *exp = new TF1("exp","exp([0])*exp([1]*x) +[2]");
  

  g->Draw();


  exp->SetParLimits(1, 0.,5.);
  //  exp->SetParLimits(2, -5.,5.);
  exp->SetParameters(-1,1.,bck);
  
  g->Fit("exp","QN","goff",0.,t_high-t_low);
  // g->Fit("exp"," "," ",0.,t_high-t_low);
  *rise_exp0=exp->GetParameter(0);
  *rise_exp1=exp->GetParameter(1);
 
  *chi=exp->GetChisquare()/exp->GetNDF();

  //  std::cout << " rise = " << *rise_exp1 << " chi2=  " << *chi << std::endl;
  exp->Delete();
  g->Delete();

  return;
}


void riselinfit(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t t_low,  Double_t t_high, Double_t *rise_lin0, Double_t *rise_lin1, Double_t *chi)
{
  TGraph *g = new TGraph();
  int npoints=0;
  int nnpoints=0;
  int Nlow=0;
  double_t amp_low =0;
  double_t amp_high =0;
  //  t_low = t_high-0.2;
  // std::cout << " t low " << t_low << " t_high " << t_high << std::endl;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)  continue;
      if(timeS[i]>t_high) break;
      npoints++;
      if (npoints == 1) 
	{
	  Nlow = i;
	}
    }
  amp_low = amp[Nlow]-bck;
  amp_high = amp[npoints]-bck;
  g->Set(npoints);
  nnpoints=0;
 
  for(int i=Nlow;i<Nlow+npoints-1;i++)
    {
      g->SetPoint(nnpoints,timeS[i],amp[i]-bck);
      //  std::cout << nnpoints << " \t" << timeS[i]<< " \t" <<amp[i]-bck << " \n";
      nnpoints++;
    }
  
  TF1 *lin = new TF1("lin","[0]+[1]*x",t_low,t_high);
  // TF1 *exp = new TF1("exp","[0]*exp([1]*x)",t_low,t_high);

  lin->SetParameters(amp_low-t_low*(amp_high-amp_low)/(t_high-t_low),(amp_high-amp_low)/(t_high-t_low));
  g->Draw();
  g->Fit("lin","QN","goff",t_low,t_high);
  // g->Fit("lin"," ","g",t_low,t_high);
  *rise_lin0=lin->GetParameter(0);
  *rise_lin1=lin->GetParameter(1);
  *chi=lin->GetChisquare()/lin->GetNDF();
  lin->Delete();
  g->Delete();

  return;
}
// Computes the mobile average of the signal
// (record lenght,       amplitude vector, time vector, number of points for the mobile average (odd number), averaged amplitude vector)
void mobileAVG(Int_t camp, Double_t amp[], Int_t navg, Double_t m_amp[])
{
  float tempsum=0;
  for(int j=int((navg-1)/2);j<camp-int((navg-1)/2)-1;j++)
    {
      tempsum=0;
      for(int k=-int((navg-1)/2);k<int((navg-1)/2)+1;k++)
	{
	  tempsum=tempsum+amp[j+k];
	}
      m_amp[j]=tempsum/float(navg);
    }
  for(int j=0;j<int((navg-1)/2);j++)
    {
      m_amp[j]=m_amp[int((navg-1)/2)];
    }
  for(int j=camp-int((navg-1)/2)-1;j<camp;j++)
    {
      m_amp[j]=m_amp[camp-int((navg-1)/2)-2];
    }

  return;
}

void Centroid_N( Double_t TIMEUNIT, Int_t Nmin, Int_t Nmax, Double_t Itot[], double *t_centr)
{
  float Num =  0;// BB RC in seconds
  float Den = 0; // BB Tau in seconds

  //tmp = new double[2*camp];//output voltage
  for (int i=Nmin;i<Nmax;i++)
    {
      Num += Itot[i]*TIMEUNIT*i;
      Den += Itot[i];
      //      tmp[i]=0;
    }
  * t_centr =  Num/Den;
  // std::cout  <<  t_centr  << " " <<  Num <<" " <<Den << "\n";
  return;
}


void CSA_NA62( Double_t TIMEUNIT, Int_t camp, Double_t Itot[],Double_t ShaperOut_Q[], Int_t *IMaxSh)
{

  *IMaxSh = 7./TIMEUNIT; //TimeUnit in [ns]
  double tr = 5.6E-9 ; 
  double til = 1.8E-9 ;
  double tac = 47.E-9 ;
  double tau = 0;
  double tt = 0;
  double qtot = 0; 
  double  Qdif_Shaper = 0;
  double sh_max = 0;
  TIMEUNIT *=1e-9; // in seconds

  for (int i=0;i< *IMaxSh;i++)  ShaperOut_Q[i] = 0;
  for (int i=0;i<camp;i++)
    {
      qtot +=  Itot[i]* TIMEUNIT*1e-6;

      if (fabs(Itot[i]) == 0) continue;
      //   std::cout <<  "i = " << i <<  " Itot[i] = " << Itot[i] << " \n";
      Qdif_Shaper  = Itot[i]*1e-6*(TIMEUNIT);// This is the charge on the Cap      

      for (int ll = 0; ll<*IMaxSh-i;ll++)
	{

  // NA62
  
	  tt = ll*TIMEUNIT;
	  tau = 0;
	  ShaperOut_Q[i+ll] += Qdif_Shaper*(-1.*tac * tr)
	    * ( ((-1. * exp(-(tt-tau)/tac) * tac) / (pow((tac-til),2.) * (tac -tr)))+((exp(-(tt-tau)/tr) * tr) / (pow((til-tr),2.) * (tac -tr))) +
		((exp(-(tt-tau)/til) * (pow(til,3.) + (tt-tau)*(tac-til)*(til-tr) -tac*til*tr)) / (pow((tac - til),2.) * til*pow((til - tr),2.))) ) ;

	  if (fabs(ShaperOut_Q[i+ll]) > fabs(sh_max))  sh_max  = ShaperOut_Q[i+ll]; 
	}

    }
  for (int ll = 0; ll<*IMaxSh;ll++)
    {
      ShaperOut_Q[ll] *=  75*1e15*qtot/sh_max; // [mV/fQ *Q/Q]
      // cout << ShaperOut_V[i] << endl;
    }
      

  // std::cout << qtot << " " << sh_max << " " << 75.*qtot/sh_max << "\n";
  // for (int ll = 0; ll<5*camp-i;ll++) std::cout  << " ll = " << ll << " time =  " << ll*TIMEUNIT*1e9 << " " << ShaperOut_Q[ll] << "\n";
  
  return;
}
void BBA( Double_t TIMEUNIT, Double_t Cdet, Double_t BBImp, Double_t BBGain, double_t BBABW, Int_t camp, Double_t Itot[],Double_t BBGraph[],Int_t *IMaxSh)
{

  *IMaxSh = 2./TIMEUNIT; //Number of time steps to be computed 
  double tau_BB_RC = 1.0e-12*BBImp*Cdet; // BB RC in seconds
  double tau_BB_BW = 0.35/(1.0e9*BBABW)/2.2; // BB Tau in seconds
  TIMEUNIT *=1e-9; // in seconds
  double Qi = 0;
  double   *VonR_BB_RC = 0;

  
  VonR_BB_RC=new double[*IMaxSh];//output voltage
  
  for (int i=0;i<*IMaxSh;i++)
    {
      BBGraph[i]=0;
      VonR_BB_RC[i] = 0;

    }
  
  for (int i=0;i<camp;i++)
    {
            
      if (i == 0) Qi = Itot[0]*1e-6*(TIMEUNIT);
      else      Qi  = (Itot[i]-Itot[i-1])*1e-6*(TIMEUNIT);// This is the charge on the Cap added in step i
      
      //     pippo += Idif/(1e-6*(TIMEUNIT));
      //	  std::cout << " i-1 = " << i-1 << " " << Itot[i-1] <<  " " << pippo << std::endl;
      //   std::cout << " i = " << i << " " <<  Qi << " " << Itot[i] << "\n";
      for (int ll = 0; ll<*IMaxSh-i;ll++)
	{
	  
 //Voltage on Rinput, this is the discharge of Qi on Rinput
	  
	  VonR_BB_RC[i+ll] += Qi/(1e-12*Cdet)*TMath::Exp(-ll*TIMEUNIT/tau_BB_RC);	  
	  BBGraph[i+ll] += 1e+3*BBGain*VonR_BB_RC[i]*(1.-TMath::Exp( (float) -ll*TIMEUNIT/tau_BB_BW));
	  
	}
    }



  delete[] VonR_BB_RC;
  return;
}
void CSA( Double_t TimeUnit, Double_t TransImp, Double_t TauFall, double_t TauRise, Int_t camp, Double_t Itot[],Double_t CSAout[], Int_t *IMaxSh)
{

  *IMaxSh = 12./TimeUnit; //TimeUnit in [ns]
  int IMax;
  float  PreAmp_Q[*IMaxSh];
  float CSAMax = 0;
  float Qtot = 0;

  for(int k=0;k<*IMaxSh;k++)
    {
      PreAmp_Q[k]=0.0;
      CSAout[k] = 0.0;
    }


  IMax = camp; // In case Itot never goes to zero  
  for(int j=0;j<camp;j++)
    {
      if ( j>2 && Itot[j] == 0)
	{
	  IMax = j;
	  break;
	}
    }
  
  for(int i=1;i<IMax;i++)
    {
      
      PreAmp_Q[i]= Itot[i]*TimeUnit;
      for (int ll = 0; ll<*IMaxSh-i;ll++)  // valid only up to IMaxSh 
	{
	  CSAout[i+ll] += TauFall/(TauFall+TauRise)*(PreAmp_Q[i])*
	    (TMath::Exp(-ll*TimeUnit/TauFall)-TMath::Exp(-ll*TimeUnit/TauRise)); // [Q] HS eq 4.3 This is the shaper	  
	  if (fabs(CSAout[i+ll]) > fabs(CSAMax))  CSAMax = CSAout[i+ll]; 
      	}
	  

      Qtot += Itot[i]*TimeUnit;
      
     }

  for (int ll = 0; ll<*IMaxSh;ll++)
    {
      CSAout[ll] *= fabs(TransImp*Qtot/CSAMax);
      //std::cout  <<  "CSA = " << CSAout[ll] << std::endl;
    }
  //   std::cout  << " Qtot = " << Qtot <<  " CSAMax  = " << CSAMax << std::endl;
  // for (int ll = 0; ll<IMaxSh;ll++)  // valid only up to IMaxSh 
  // std::cout <<  " Total Charge = " << Qtot << " Tpeak = " << Tpeak << " Ball. Deficit = " << BallDef << " CSAMax = "<< CSAMax << std::endl;       

  
  return;
}


//
//Calcola la derivata a quattro punti
// (record lenght, amplitude vector,reshaping vector)
void Derivative( Double_t TimeUnit,Int_t camp, Double_t amp[],Double_t resh[])
{

  TGraph *g = new TGraph();
  int npoints=3;
  int nnpoints=0;
  // int Nlow=0;
  g->Set(2*npoints);
 
  for(int j=npoints;j<camp-npoints;j++)
	{
	   nnpoints=0;
	  for (int kk = j-npoints; kk< j+npoints; kk++)
	    {
	      g->SetPoint(nnpoints,nnpoints*TimeUnit,amp[kk]);
	      //	      std::cout << "Fit = " <<  nnpoints << " \t" << nnpoints*TimeUnit << " \t" <<amp[kk] << " \n";
	      nnpoints++;
	    }	    
	  // resh[j]=(amp[j+10]-amp[j-10])/(2*npoints*TimeUnit);
	  TF1 *lin = new TF1("lin","[0]+[1]*x",0,2*nnpoints*TimeUnit);
	  g->Fit("lin","QN","goff",0,2*nnpoints*TimeUnit);
	  //*rise_lin1=lin->GetParameter(1);
	  resh[j]=lin->GetParameter(1);
	  //	  std::cout << "j = " << j << " "  << resh[j] << " " << *rise_lin1 << std::endl;
	  lin->Delete();
	}
   g->Delete();

return;
}



// Generates FFT of a real vector (inreal) of lenght camp. Real part of output FFT is outre, immaginary part is outim. Module output is outmod.
void FFTrealtocomplex(Int_t camp, Double_t inreal[], Double_t outre[], Double_t outim[], Double_t outmod[])
{
  Double_t re;
  Double_t im;
  Double_t *in = new Double_t[camp];
  for(int i=0;i<camp;i++)
    {
      *(&in[i])=inreal[i];
    }
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &camp, "R2C");
  fftr2c->SetPoints(in);
  fftr2c->Transform();
  for (Int_t i=0; i<camp; i++)
    {
      fftr2c->GetPointComplex(i, re, im);
      outre[i]=re;
      outim[i]=im;
      outmod[i]=sqrt(re*re+im*im);
    }
  fftr2c->Delete();
  delete[] in;

  return;
}


// Generates the inverse FFT from a complex vector of lenght camp (two vectors: inre real part, inim immaginary part). The output vector is outreal.
void FFTcomplextoreal(Int_t camp, Double_t inre[], Double_t inim[], Double_t outreal[])
{
  Double_t *imma = new Double_t[camp];
  Double_t *real = new Double_t[camp];
  
  for(int i=0;i<camp;i++)
    {
      *(&real[i])=inre[i];
      *(&imma[i])=inim[i];
    }
  TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &camp, "C2R");
  fftc2r->SetPointsComplex(real,imma);
  fftc2r->Transform();
  for (Int_t i=0; i<camp; i++)
    {
      outreal[i]=fftc2r->GetPointReal(i)/double(camp);
    }
  fftc2r->Delete();
  delete[] imma;
  delete[] real;

  return;
}


// Redines coordinates correction distorsion of centroid. polinomial correction of coordinates (expansion at center of strips and compressio at borders)
//(centroid variable address, curve parameter)
void positionbias(Double_t *bar,Double_t sca)
{
  sca = 0;
  float barfz=*bar;
  float xloop=(1/100000.)*(int((100000.)*(barfz-4))%800000);
  
  *bar=((1/(70.+16.))*((70.)*(xloop-4)+pow((xloop-4),3))+(barfz-xloop-4)+8);

  return;
}


//Function to fit the waveform with Gumbel function
void Gumbel_fit(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t t_low,Double_t t_up,Double_t t_amp, Double_t t_width, Double_t max, Double_t *bck_fit, Double_t *t_amp_fit, Double_t *t_width_fit, Double_t *amp_fit)
{
  TF1 *fgum = new TF1("gumbel","[0]+[3]*(exp(-(x-[1])/[2]+exp(-(x-[1])/[2])))",t_low,t_up);
  fgum->SetParameters(bck,t_amp,t_width,max);

  TGraph *g = new TGraph();
  int npoints=0;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)
	continue;
      if(timeS[i]>t_up)
	break;
      npoints++;
    }
  g->Set(npoints);
  npoints=0;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)
	continue;
      if(timeS[i]>t_up)
	break;
      g->SetPoint(npoints,timeS[i],amp[i]);
      npoints++;
    }

  g->Draw();
  g->Fit("gumbel","QN","goff",t_low,t_up);
  *bck_fit=fgum->GetParameter(0);
  *t_amp_fit=fgum->GetParameter(1);
  *t_width_fit=fgum->GetParameter(2);
  *amp_fit=fgum->GetParameter(3);

  g->Delete();
  fgum->Delete();

  return;
}

double string_to_double( const std::string& s )
 {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 } 





void ConstFractime_N(Double_t amp[], Double_t DT, Double_t bck,Double_t threshold, Double_t max, Int_t NMax,Double_t *thrtime, bool Interpolation, Int_t  *Nth)
{    
  *thrtime=0;

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  //  for(int j=0;j<camp;j++)
    for(int j=NMax;j>1;j--)
    {
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j+1]-bck)/max; // Above Vth; too late
	      C1 = (amp[j]-bck)/max;  // Below Vth; too early
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
	    }
	  else   TCrossing =  timeS[j];

	  break;	  
	}
	
    }
    *thrtime= TCrossing;

  return;
}




void Trailtime_N(Int_t camp, Double_t amp[], Double_t DT, Int_t NMax, Double_t bck,Double_t threshold, Double_t *thrtime, bool Interpolation, Int_t  *Nth)
{
  
  *thrtime=0;
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  for(int j=NMax;j<camp;j++)
    {
      if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	{
	    if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j-1]-bck); //above Vth, too early
	      C1 = (amp[j]-bck);  //below Vth, too late	      
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1);	  
	    }
	  else   TCrossing =  timeS[j];

	  break;
	  
	}
    }
    
  *thrtime= TCrossing;
  return;
}

void TrailConstFractime_N(Int_t camp, Double_t amp[], Double_t DT,Double_t bck,Double_t threshold, Double_t max, Int_t NMax,Double_t *trailtime, bool Interpolation, Int_t  *Nth)
{
  *trailtime=0;
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  for(int j=NMax;j<camp;j++)
    {

      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j-1]-bck)/max; //above Vth, too early
	      C1 = (amp[j]-bck)/max;  // below vth, too late
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1); //timeS[j] is too late, I need to subtract  a bit of DT
	      //	      std:: cout << " C1 = " << C1 << " C2 = " << C2 << " threshold = " << threshold << std:: endl;
	      //
	      // std:: cout <<  timeS[j-1] << "  threshold - C1 = " << threshold - C1 << " C2 - C1 = " << C2-C1 << " TCrossing = " << TCrossing << std:: endl; 
	    }
	  else   TCrossing =  timeS[j];
	  
	  break;	  
	}
      
    }
  *trailtime= TCrossing;  
  return;
}


#endif


void Levtime_N( Double_t amp[], Double_t DT, Int_t NMax, Double_t bck,Double_t threshold, Double_t *thrtime, bool Interpolation, Int_t  *Nth)
{

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  for(int j=NMax;j>0;j--)
    {
 
      if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	{
	  if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j+1]-bck); // above Vth
	      C1 = (amp[j]-bck); //Below Vth
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1);  //time[j] is too early, so I need to add a fraction of time bin  
	    }
	  else   TCrossing =  timeS[j];

	  break;
	  
	}
    }
    *thrtime= TCrossing;    

  return;
}


void Levtime_N_array(int NTValues[], double TValues[], Int_t NArray, Double_t Threshold[], Double_t amp[], Double_t DT, Int_t NMax, Double_t bck, bool Interpolation)
{

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;

  NArray--;
  NTValues[NArray] = -100;
  TValues[NArray] = -100;
  int Nth;
  for(int j=NMax;j>0;j--)
    {
      if (NArray < 0) break;
      threshold = Threshold[NArray]; 
      if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	{
	  if (Interpolation)
	    {
	      Nth = j;
	      C2 = (amp[j+1]-bck); // above Vth
	      C1 = (amp[j]-bck); //Below Vth
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1);  //time[j] is too early, so I need to add a fraction of time bin

	      // std::cout << "Crossing 2-points: = "<< TCrossing << std::endl;
	      // std::cout << "Level = "<<  threshold << std::endl;
	      
	      //	      TF1 *f1 = new TF1("f1","[0]+[1]*x");
	      if (1>2) // dofit
		{
		  TF1 *f1 = new TF1("f1","pol2");
		  TGraph *g = new TGraph();
		  g->Set(7);
		  for(int b=0;b<7;b++)
		    {
		      g->SetPoint(b,amp[j-2+b]-bck,timeS[j-2+b]);
		      //  std::cout <<amp[j-2+b]-bck << " " << timeS[j-2+b]<<std::endl;
		    }
		  g->Fit("f1","QN");
		  TCrossing = f1->GetParameter(0)+ threshold*f1->GetParameter(1)+ threshold*threshold*f1->GetParameter(2);
		  //  std::cout << "Parameter Crossing fits: = "<<f1->GetParameter(0) << " " << f1->GetParameter(1) <<  " " << f1->GetParameter(2)<< std::endl;
		  // std::cout << "Crossing fits: = "<< TCrossing << std::endl;
		  g->Delete();
		  f1->Delete();
		}
	    }
	  
	  else   TCrossing =  timeS[j];
	  
	  NTValues[NArray] = Nth;
	  TValues[NArray] = TCrossing ;
	  
	  //	  std::cout << " Level th = " <<   Threshold[NArray]  << " at time " <<  TValues[NArray] <<  " at index = " << NTValues[NArray] << std::endl;
	  NArray--;
	  while(amp[j]-bck<Threshold[NArray]) j++;
	    

	  
	}
    }

  return;
}

void TrailLevtime_N_array(double TValues[], Int_t NArray, Double_t Threshold[], Int_t camp, Double_t amp[], Double_t DT, Int_t NMax, Double_t bck, bool Interpolation)
{

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;

  NArray--;
  TValues[NArray] = -100;
  int Nth;
  for(int j=NMax;j<camp;j++)
    {
       if (NArray < 0) break;
      threshold = Threshold[NArray]; 
      if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	{
	  if (Interpolation)
	    {
	      Nth = j;
	      C2 = (amp[j-1]-bck); // above Vth
	      C1 = (amp[j]-bck); //Below Vth
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1);  //time[j] is too early, so I need to add a fraction of time bin  
	    }
	  else   TCrossing =  timeS[j];
	  TValues[NArray] = TCrossing ;

	  NArray--;
	  while(amp[j]-bck<Threshold[NArray]) j--;

	  
	}
    }

  return;
}

void ConstFractime_N_array(int NTValues[], double TValues[], Int_t NArray, Double_t Threshold[], Double_t amp[], Double_t DT, Double_t bck, Double_t max, Int_t NMax, bool Interpolation)
{    
 
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;
  int Nth;
  // double thrtime;

  //  for(int j=0;j<camp;j++)
  
  NArray--;
  NTValues[NArray] = -100;
  TValues[NArray] = -100;
    for(int j=NMax;j>1;j--)
    {
      if (NArray < 0) break;
      threshold = Threshold[NArray]; 
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  if (Interpolation)
	    {
	      Nth = j;
	      C2 = (amp[j+1]-bck)/max; // Above Vth; too late
	      C1 = (amp[j]-bck)/max;  // Below Vth; too early
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
	    }
	  else   TCrossing =  timeS[j];
	  NTValues[NArray] = Nth;
	  TValues[NArray] = TCrossing ;
	  
	  //	  std::cout << " Frac th = " <<   Threshold[NArray]  << " at time " <<  TValues[NArray] <<  " at index = " << NTValues[NArray] << std::endl;
	  NArray--;	  
	}
	
    }
  
    //    *ANTValues = NNTValues[0];
    // std::cout << " &NNTValues[0]  " <<  &NNTValues[0] << std::endl;
    // std::cout << " * ANTValues " <<   *ANTValues  << std::endl;
    // std::cout << " ANTValues " <<   ANTValues  << std::endl;
    // std::cout << " &ANTValues " <<   &ANTValues  << std::endl;
    return;
}


void TrailConstFractime_N_array(int NTValues[], double TValues[], Int_t NArray, Double_t Threshold[], Int_t camp, Double_t amp[], Double_t DT, Double_t bck, Double_t max, Int_t NMax, bool Interpolation)
{    
 
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;
  int Nth;
  // double thrtime;

  //  for(int j=0;j<camp;j++)
  
  NArray--;
  TValues[NArray] = -100;
  NTValues[NArray] = -100;
  
    for(int j=NMax;j<camp;j++)
    {
      if (NArray < 0) break;
      threshold = Threshold[NArray]; 
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  if (Interpolation)
	    {
	      Nth = j;
	      C2 = (amp[j-1]-bck)/max; // Above Vth; too late
	      C1 = (amp[j]-bck)/max;  // Below Vth; too early
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
	    }
	  else   TCrossing =  timeS[j];
	  TValues[NArray] = TCrossing ;
	  NTValues[NArray] = Nth;
	  
	  //	  std::cout << " Frac th = " <<   Threshold[NArray]  << " at time " <<  TValues[NArray] <<  " at index = " << NTValues[NArray] << std::endl;
	  NArray--;	  
	}
	
    }
  
    return;
}



void calculateThreeCircleIntersectionN(Int_t xc[], Int_t yc[], double rc[], double am[])
			  
{
  double a, dx, dy, d, h, rx, ry;
  double point2_x, point2_y;
  double EPSILON = 20;
  int i,j;
  XCross =-100;
  YCross = -100;
  for (i = 0; i < 3 ; i++)
    {
       for (j = i; j < 4; j++)
         if (am[j] > am[i])
	   {
	     //  cout << i << " " << j << endl;
	     temp = am[i];
	     am[i] = am[j];
	     am[j] = temp;

	     
	     temp = xc[i];
	     xc[i] = xc[j];
	     xc[j] = temp;


	     
	     temp = rc[i];
	     rc[i] = rc[j];
	     rc[j] = temp;
	     
	     	     
	     temp = yc[i];
	     yc[i] = yc[j];
	     yc[j] = temp;
	     /**/
	     
	     
	   }
    }
  //  std::cout << "pippo " << am[0] << " " << am[1] << " " << am[2] << " " << am[3] <<std::endl;
 
  

    /* dx and dy are the vertical and horizontal distances between
    * the circle centers.
    */
    dx = xc[1] - xc[0];
    dy = yc[1] - yc[0];

    /* Determine the straight-line distance between the centers. */
    d = sqrt((dy*dy) + (dx*dx));

    //    std::cout << "distance between centers: " << d << std::endl;
    /* Check for solvability. */
    if (d-EPSILON > (rc[0] + rc[1]))
    {
        /* no solution. circles do not intersect. */
        return;
    }
    if (d+EPSILON < fabs(rc[0] - rc[1]))
    {
        /* no solution. one circle is contained in the other */
        return;
    }

    /* 'point 2' is the point where the line through the circle
    * intersection points crosses the line between the circle
    * centers.
    */

    /* Determine the distance from point 0 to point 2. */
    a = ((rc[0]*rc[0]) - (rc[1]*rc[1]) + (d*d)) / (2.0 * d) ;

    //  std::cout << "point 2  " <<  a << "," << r0 <<  "," << r1 << "," << d << std::endl;
    /* Determine the coordinates of point 2. */
    point2_x = xc[0] + (dx * a/d);
    point2_y = yc[0] + (dy * a/d);
    
    //    std::cout << "point 2  " <<  point2_x << x0 <<  "," << dx*a/d <<  std::endl;
    /* Determine the distance from point 2 to either of the
    * intersection points.
    */
    h = sqrt((rc[0]*rc[0]) - (a*a));

    /* Now determine the offsets of the intersection points from
    * point 2.
    */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    /* Determine the absolute intersection points. */
    double intersectionPoint1_x = point2_x + rx;
    double intersectionPoint2_x = point2_x - rx;
    double intersectionPoint1_y = point2_y + ry;
    double intersectionPoint2_y = point2_y - ry;

    //    std::cout << "INTERSECTION Circle1 AND Circle2: (" <<  intersectionPoint1_x << "," << intersectionPoint1_y<< ") AND ("<< intersectionPoint2_x << "," << intersectionPoint2_y <<")" << std::endl;
    //    XCross =(intersectionPoint2_x+intersectionPoint1_x)/2;
    // YCross = (intersectionPoint2_y+intersectionPoint1_y)/2;      
    
    /* Lets determine if circle 3 intersects at either of the above intersection points. */
    dx = intersectionPoint1_x - xc[2];
    dy = intersectionPoint1_y - yc[2];
    double d1 = sqrt((dy*dy) + (dx*dx));

    dx = intersectionPoint2_x - xc[2];
    dy = intersectionPoint2_y - yc[2];
    double d2 = sqrt((dy*dy) + (dx*dx));

    if(fabs(d1 - rc[2]) < EPSILON) {
      XCross =intersectionPoint1_x;
      YCross = intersectionPoint1_y;
     // std::cout << "INTERSECTION Circle1 AND Circle2 AND Circle3: (" <<  XCross << "," << YCross << ")"<< std::endl;


    }
    else if(fabs(d2 - rc[2]) < EPSILON) {
      XCross =intersectionPoint2_x;
      YCross = intersectionPoint2_y;      
      //  std::cout << "INTERSECTION Circle1 AND Circle2 AND Circle3: (" <<  XCross << "," << YCross << ")"<< std::endl;
    }
    else {
      //        Log.d("INTERSECTION Circle1 AND Circle2 AND Circle3:", "NONE");
    }
    return;
}

void SetRootFile(bool doFit, bool showFFT, std::string FileN)
//,TFile *OutputFile,TTree *OutTree)
{

  OutputFile = new TFile(FileN.c_str(),"recreate"); 
  OutTree = new TTree("Analysis","Analysis");
  OutTree->SetDirectory(OutputFile);

  // Analysis branches
  OutTree->Branch("Wafer",&Wafer,"Wafer/I");
  OutTree->Branch("FType",&FType,"FType/I");
  OutTree->Branch("LGain",&LGain,"LGain/I");
  OutTree->Branch("DMetal",&DMetal,"DMetal/I");
  OutTree->Branch("DPitch",&DPitch,"DPitch/I");
  OutTree->Branch("ntrig",&ntrig,"ntrig/I");
  OutTree->Branch("ds",&ds,"ds/I");
  OutTree->Branch("event",&event,"event/I");
  OutTree->Branch("nrun",&nrun,"nrun/I");
  OutTree->Branch("nchro",&nchro,"nchro/I");
  OutTree->Branch("CAmax",&CAmax,"CAmax/I");
  OutTree->Branch("XTime",&XTime,"XTime/D");
  OutTree->Branch("YTime",&YTime,"YTime/D");
  OutTree->Branch("XCord",&XCord,"XCord/D");
  OutTree->Branch("YCord",&YCord,"YCord/D");
  OutTree->Branch("XAmpl",&XAmpl,"XAmpl/D");
  OutTree->Branch("YAmpl",&YAmpl,"YAmpl/D");
  OutTree->Branch("XArea",&XArea,"XArea/D");
  OutTree->Branch("YArea",&YArea,"YArea/D");
  OutTree->Branch("XCirc",&XCirc,"XCirc/D");
  OutTree->Branch("YCirc",&YCirc,"YCirc/D");
  OutTree->Branch("TTCfd30",&TTCfd30,"TTCfd30/D");
  OutTree->Branch("TTCfd50",&TTCfd50,"TTCfd50/D");
  OutTree->Branch("TCfd20",&TCfd20,"TCfd20/D");
  OutTree->Branch("TCfd30",&TCfd30,"TCfd30/D");
  OutTree->Branch("TCfd40",&TCfd40,"TCfd40/D");
  OutTree->Branch("TCfd50",&TCfd50,"TCfd50/D");
  OutTree->Branch("TTmax",&TTmax,"TTmax/D");
  
  

  OutTree->Branch("t_bck",t_bck,"t_bck[nchro]/D");
  OutTree->Branch("t_pul",t_pul,"t_pul[nchro]/D");

  OutTree->Branch("maxD",maxD,"maxD[nchro]/D");
  OutTree->Branch("area",area,"area[nchro]/D");
  OutTree->Branch("ampl",ampl,"ampl[nchro]/D");
  OutTree->Branch("ampl_chi2",ampl_chi2,"ampl_chi2[nchro]/D");
  OutTree->Branch("dVdt3070",dVdt3070,"dVdt3070[nchro]/D");
  OutTree->Branch("dVdt1030",dVdt1030,"dVdt1030[nchro]/D");
  OutTree->Branch("dVdt2080",dVdt2080,"dVdt2080[nchro]/D");

  OutTree->Branch("bck",bck,"bck[nchro]/D");
  OutTree->Branch("max_bck_before",max_bck_before,"max_bck_before[nchro]/D");
  OutTree->Branch("rms_bck_before",rms_bck_before,"rms_bck_before[nchro]/D");
  OutTree->Branch("max_bck_after",max_bck_after,"max_bck_after[nchro]/D");
  OutTree->Branch("rms_bck_after",rms_bck_after,"rms_bck_after[nchro]/D");


  OutTree->Branch("t_centMax",t_centMax,"t_centMax[nchro]/D");
  
  OutTree->Branch("t_level10",t_level10,"t_level10[nchro]/D");
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

  OutTree->Branch("trail_t_level30",trail_t_level30,"trail_t_level30[nchro]/D");
  OutTree->Branch("trail_t_level40",trail_t_level40,"trail_t_level40[nchro]/D");
  OutTree->Branch("trail_t_level80",trail_t_level80,"trail_t_level80[nchro]/D");
  OutTree->Branch("trail_t_level100",trail_t_level100,"trail_t_level100[nchro]/D");
  OutTree->Branch("trail_t_level200",trail_t_level200,"trail_t_level200[nchro]/D");
  OutTree->Branch("trail_t_level300",trail_t_level300,"trail_t_level300[nchro]/D");

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

  OutTree->Branch("Ncfd05",Ncfd05,"Ncfd05[nchro]/I");
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

  
  OutTree->Branch("XDist",XDist,"XDist[nchro]/D");
  OutTree->Branch("YDist",YDist,"YDist[nchro]/D");
  OutTree->Branch("dist",dist,"dist[nchro]/D");
  OutTree->Branch("ind",ind,"ind[nchro]/I");
  OutTree->Branch("XPad",XPad,"XPad[nchro]/I");
  OutTree->Branch("YPad",YPad,"YPad[nchro]/I");
  
  OutTree->Branch("trail_cfd30",trail_cfd30,"trail_cfd30[nchro]/D");
  OutTree->Branch("trail_cfd50",trail_cfd50,"trail_cfd50[nchro]/D");
  OutTree->Branch("trail_cfd70",trail_cfd70,"trail_cfd70[nchro]/D");
  OutTree->Branch("trail_cfd90",trail_cfd90,"trail_cfd90[nchro]/D");

  OutTree->Branch("samples",samples,"samples[nchro]/I");
  //  OutTree->Branch("samplesrec",samples,"samples[nchro]/I");
  OutTree->Branch("t_max",t_max,"t_max[nchro]/D");
  
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
  // std::cout << "nchro = "<< nchro << std::endl;
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
    
      leaf << "amp" << i;
      leafl << "amp" << i <<"[samplesrec_" << i << "]/D";
      OutTree->Branch(leaf.str().c_str(),&amprec[i][0],leafl.str().c_str());
      
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "m_amp" << i;
      leafl << "m_amp" << i <<"[samplesrec_" << i << "]/D";
	    //leafl << "m_amp" << i <<"[samples[" << i << "]]/D";
      OutTree->Branch(leaf.str().c_str(),&m_amprec[i][0],leafl.str().c_str());
      //      cout << " =  " << leafl.str().c_str() << endl;

      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "der_amp" << i;
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
  std::cout << " sam " << samplesrec << std::endl;
  OutTree->Branch("time",timerec,"time[samplesrec_0]/D"); 
  OutTree->Branch("freq",freq,"freq[samples_0]/D");
  

  
  std::cout << "Branch settings done" << std::endl;
  return; 
}
