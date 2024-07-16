//// KinematicDistributions.h //////////////////////////////////////////////////////////////////////////
//
// collection of momentum spectra of various particles
//
// currently implemented:
//
//   TF1 Hagedorn (name, mass, upperlim, lowerlim, A=504.5, a=0.5169, b=0.1626, p0=0.7366, n=-8.274)
//   TF1 HagedornMt (name, mass, upperlim, lowerlim, A=504.5, a=0.5169, b=0.1626, p0=0.7366, n=-8.274)
//   TF1 HagedornYield (name, mass, upperlim, lowerlim, A=504.5, a=0.5169, b=0.1626, p0=0.7366, n=-8.274)
//       default parameters for Au+Au 200 (see Alan Dion thesis for details)
//       pp 200 parameters: 377., 0.356, 0.068, 0.7, -8.25
//
//   TF1 PromptPhotonYield(name, upperlim, lowerlim, a=0.0066, b=6.4, c=0.4, n=17.6, s=200.)
//       Prompt photon function fit from PPG140 fit to pp 200
//
//   double  Weight_GPR_pion_pp200(pt, opt=0)
//           - uses GPR for pi0 from ppg202, and pi+/pi-  ppg030 & ppg101  - Roli Esha 11/10/2021  
//           - return value+opt*error
//   double  Weight_GPR_etapi_uni(pt, opt=0)
//           - universal eta/pi ratio fro Yuanjie Ren thesis  
//
//   TF1 SPSptYield(name,mass,ptmin=0,ptmax=4)
//   TF1 SPSrapidity(name,mass,ymin,ymax,ebeam=200)
//       pt and rapidity functions for SPS fix target collisions (taken from EXODUS)
//
//   TF1 Flat(const Char_t *name, Double_t min=0, Double_t max = 2*3.1415926)  
//       generates generic flat distribution between min,max default assumes 0-2pi
//
//   double Pt_Mt_Scaled(pt,m,mref)
//       mt scales pt from particle of mass m based on pt of particle with mass mref
//   double RapidtyToEta(y,pt,m)
//       calculates eta from y,pt,m
//
// Axel Drees  11/11/2021 
// Axel Drees  3/18/2022    add SPS kinematic distributions
// Axel Drees  6/9/2022     add conversion from rapidity to eta, and mt scaling from m to m* 
// Axel Drees  11/17/2022   debug and fix SPS rapidity distributions 
//

#ifndef KinDist_h
#define KinDist_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/GPRFiles/GPRFileNames.h"

//*******************************************************************
//
// Hagedorn Function for invariant yield mesons, 1/pt dN/dpt
//
// see thesis of Alan Dion for papameters for pi0
// all other particle distributions are derived using mt scaling
// assumes mt scaling 
//
Double_t Hagedorn_Func(Double_t *x, Double_t *p) {
//
// Hagedorn function in 1/pt dN/pt from ppg088 
//
  Double_t pt   = x[0];
  Double_t mass = p[0];
  Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.13498*0.13498);
  Double_t A    = p[1];
  Double_t a    = p[2];
  Double_t b    = p[3];
  Double_t p0   = p[4];
  Double_t n    = p[5];
 
  Double_t value = A* pow( exp(-a*mt-b*mt*mt)+mt/p0 ,n);
  return value;
}
//*******************************************************************
//
// Hagedorn Function for mesons yield, dN/dpt
//
// see thesis of Alan Dion for papameters for pi0
// all other particle distributions are derived using mt scaling
// assumes mt scaling 
//
Double_t Hagedorn_Yield_Func(Double_t *x, Double_t *p) {
//
// Hagedorn function in 1/pt dN/pt from ppg088 
//
  Double_t pt   = x[0];
  Double_t mass = p[0];
  Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.13498*0.13498);
  Double_t A    = p[1];
  Double_t a    = p[2];
  Double_t b    = p[3];
  Double_t p0   = p[4];
  Double_t n    = p[5];
 
  Double_t value = pt*A/2/3.141* pow( exp(-a*mt-b*mt*mt)+mt/p0 ,n);
  return value;
}

Double_t Hagedorn_Func1(const Double_t *x, const Double_t *p) {
//
// Hagedorn function in 1/pt dN/pt from ppg088 
//
  Double_t pt   = x[0];
  Double_t mass = p[0];
  Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.13498*0.13498);
  Double_t A    = p[1];
  Double_t a    = p[2];
  Double_t b    = p[3];
  Double_t p0   = p[4];
  Double_t n    = p[5];
 
  Double_t value = A*TMath::Erf(mt)* pow( exp(-a*mt-b*mt*mt)+mt/p0 ,n);
  return value;
}




///////////////////////////////////////////////////////////////////////////////
TF1 *HagedornYield(const Char_t *name, Double_t mass, Double_t upperlim, Double_t lowerlim, Double_t A=504.5, Double_t a=0.5169, Double_t b=0.1626, Double_t p0=0.7366, Double_t n=-8.274){
//
// modified Hagedorn function used in ppg088, default parameters for min bias AuAu as documented in Alan Dion's PhD thesis
//
//
  if (lowerlim <= 0.) lowerlim = 0.001;
  TF1 *fHagedorn = new TF1(name, Hagedorn_Yield_Func, lowerlim, upperlim, 6);
  fHagedorn->FixParameter(0,mass); 
  fHagedorn->SetParameters(mass, A, a, b, p0, n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
    }

///////////////////////////////////////////////////////////////////////////////
TF1 *Hagedorn(const Char_t *name, Double_t mass, Double_t upperlim, Double_t lowerlim, Double_t A=504.5, Double_t a=0.5169, Double_t b=0.1626, Double_t p0=0.7366, Double_t n=-8.274){
//
// modified Hagedorn function used in ppg088, default parameters for min bias AuAu as documented in Alan Dion's PhD thesis
//
//
  if (lowerlim <= 0.) lowerlim = 0.001;
  TF1 *fHagedorn = new TF1(name, Hagedorn_Func, lowerlim, upperlim, 6);
  fHagedorn->FixParameter(0,mass); 
  fHagedorn->SetParameters(mass, A, a, b, p0, n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
    }


#include "TMath.h"

//*******************************************************************

TF1 *HagedornErf(const Char_t *name, Double_t mass, Double_t upperlim, Double_t lowerlim, Double_t A=504.5, Double_t a=0.5169, Double_t b=0.1626, Double_t p0=0.7366, Double_t n=-8.274){
//
// modified Hagedorn function used in ppg088, default parameters for min bias AuAu as documented in Alan Dion's PhD thesis
//
//
  TF1 *fHagedorn = new TF1(name, Hagedorn_Func1, lowerlim, upperlim, 6);
  fHagedorn->FixParameter(0,mass); 
  fHagedorn->SetParameters(mass, A, a, b, p0,n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
    }
//


//*******************************************************************
//
Double_t PromptPhoton_Func(Double_t *x, Double_t *p) {
//
// Prompt photon function fit from PPG140 
//
// the units are in mb/GeV^2
//
  Double_t pt    = x[0];
  Double_t a     = p[0];
  Double_t b     = p[1];
  Double_t c     = p[2];
  Double_t n     = p[3];
  Double_t s     = p[4];
  Double_t xt    = 2.*pt/s;
  Double_t value = pt* a / pow(pt, b+c*log(xt)) * pow(1-xt*xt,n);
//  std::cout << "... fit ... pt " << pt << " xt " << xt << " a " << a << " value " << value << std::endl;
  return value;
}

TF1 *PromptPhotonYield(const Char_t *name, Double_t upperlim, Double_t lowerlim, Double_t a=0.0066, Double_t b=6.4, Double_t c=0.4, Double_t n=17.6, Double_t s=200.){
//
//
  if (lowerlim <= 0.) lowerlim = 0.001;
  TF1 *fHagedorn = new TF1(name, PromptPhoton_Func, lowerlim, upperlim, 5);
  fHagedorn->FixParameter(5,s); 
  fHagedorn->SetParameters(a, b, c, n, s); 
  fHagedorn->SetParNames("a", "b", "c", "n", "s");
  return fHagedorn;
    }

////////////////////////////////////////////////////////////////////////////////////////////
//
// distributions created with Gaussian Process Regression
//    
// txt files located in /GPRfiles/... 
//
//    GPR_pion_pp200.txt - uses pi0 from ppg202, and pi+/pi-  ppg030 & ppg101  - Roli Esha 11/10/2021
//    GPR_etapi_universal.txt - for pp  and pA collisons at all energies - from Yuanjie Ren MA thesis 2020
//
const int GPR_pi_pp200_n=1000;
const double GPR_pi_pp200_ptmax = 25.;   // GeV/c
const double GPR_pi_pp200_ptmin = 0.3;   // GeV/c
const double GPR_pi_pp200_gap = (GPR_pi_pp200_ptmax-GPR_pi_pp200_ptmin)/GPR_pi_pp200_n;

const int GPR_etapi_uni_n=500;
const double GPR_etapi_uni_ptmax = 1.477;  // log10(pt in GeV/c)
const double GPR_etapi_uni_ptmin = -1.8;    // log10(pt in GeV/c)
const double GPR_etapi_uni_gap = (GPR_etapi_uni_ptmax-GPR_etapi_uni_ptmin)/GPR_etapi_uni_n;

double GPR_pi_pp200_pt[GPR_pi_pp200_n], GPR_pi_pp200_yield[GPR_pi_pp200_n], GPR_pi_pp200_yield_error[GPR_pi_pp200_n];
double GPR_etapi_uni_pt[GPR_etapi_uni_n], GPR_etapi_uni[GPR_etapi_uni_n], GPR_etapi_uni_error[GPR_etapi_uni_n];

//////// pp 200 pion pt distribution ////////////////////////////////////////////////////////////////
// read textfile

void Read_GPR_pion_pp200(){
  std::ifstream infile;
  infile.open(GPR_pion_pp200);
  for(int i = 0; i < GPR_pi_pp200_n; i++){
    infile >> GPR_pi_pp200_pt[i] >> GPR_pi_pp200_yield[i] >> GPR_pi_pp200_yield_error[i];
  }
  infile.close();
}
// calculate weight 
double Weight_GPR_pion_pp200(double pt, double opt=0){
  
  double gap, value=0, error=0, target_index, alpha;
  int index;
  Read_GPR_pion_pp200();

  if (pt<=GPR_pi_pp200_ptmin or pt>=GPR_pi_pp200_ptmax) return 0;
  gap = GPR_pi_pp200_gap;
  target_index = (pt-GPR_pi_pp200_ptmin)/gap;
  index = floor(target_index);
  alpha = target_index - index;

  if (index < GPR_pi_pp200_n-1) {
     value = GPR_pi_pp200_yield[index]*(1-alpha) + GPR_pi_pp200_yield[index+1]*alpha; 
     error = GPR_pi_pp200_yield_error[index]*(1-alpha) + GPR_pi_pp200_yield_error[index+1]*alpha;
  }
  
  return value+opt*error;
}
////////// eta to pi0 ration univesal function valid for all pp and pA collisions //////////////////////
void Read_GPR_etapi_uni(){
  std::string trashline;
  std::ifstream infile;
  infile.open(GPR_etapi_universal);
  getline(infile,trashline); 

  getline(infile,trashline); 

  for(int i = 0; i < GPR_etapi_uni_n; i++){
    infile >> GPR_etapi_uni_pt[i] >> GPR_etapi_uni[i] >> GPR_etapi_uni_error[i];

  }
  infile.close();
}
double Weight_GPR_etapi_uni(double pt, double opt=0){
  double gap, lpt, value=0, error=0, target_index, alpha;
  int index;
  Read_GPR_etapi_uni();
  lpt = log10(pt);

  if (lpt>=GPR_etapi_uni_ptmax) return 0;
  if (lpt<=GPR_etapi_uni_ptmin) return GPR_etapi_uni[0]+opt*GPR_etapi_uni_error[0];
  gap = GPR_etapi_uni_gap; 
  target_index = (lpt-GPR_etapi_uni_ptmin)/gap;
  index = floor(target_index);
  alpha = target_index - index;

  if (index < GPR_etapi_uni_n-1) {
     value = GPR_etapi_uni[index]*(1-alpha) + GPR_etapi_uni[index+1]*alpha; 
     error = GPR_etapi_uni_error[index]*(1-alpha) + GPR_etapi_uni_error[index+1]*alpha;
  }  
  return value+opt*error;
}

double Pt_Mt_Scaled(double pt, double m, double m_ref){

  double pt_ref = sqrt(pt*pt + m*m - m_ref*m_ref);
  return pt_ref;
}


//-----------------------------------------------------------------------------
//
// Generate transverse-momentum distribution of pions for top SPS energies 
// based on parameterization from CERES (included in EXODUS) 
//
//-----------------------------------------------------------------------------

double SPSpt_Func(double *x, double *p)
{
  double pt = x[0];
  double mass = p[0];
  double a1 = p[1];
  double a2 = p[2];
  double a3 = p[3];
  double T1 = p[4];
  double T2 = p[5];
  double T3 = p[6];
  double mt = sqrt(pt*pt+mass*mass);
  double weight;

  weight = pt*(a1*exp(-1.0*mt/T1)+a2*exp(-1.0*mt/T2)+a3*exp(-1.0*mt/T3));

  return weight;
}

TF1 *SPSptYield(const Char_t *name, Double_t mass, Double_t ptmin=0, Double_t ptmax=4.0){
//
//
  const double a1 = 1.0;
  const double a2 = 0.139;
  const double a3 = 0.107;
  const double T1 = 0.1;
  const double T2 = 0.23;
  const double T3 = 0.102;
  if (ptmin <= 0.) ptmin = 0.001;
  TF1 *fSPSptYield = new TF1(name, SPSpt_Func, ptmin, ptmax, 7);
  fSPSptYield->SetParameters(mass, a1, a2, a3, T1, T2, T3); 
  fSPSptYield->FixParameter(0,mass); 
  fSPSptYield->SetParNames("mass","a1", "a2", "a3", "T1", "T2", "T3");
  return fSPSptYield;
}


//-----------------------------------------------------------------------------
//
// Generate rapidity distributions for the CERES generator
//
//-----------------------------------------------------------------------------
double SPSrapidity_Func(double *x, double *p)
{
  double s2pi = 0.39894245; // 1/sqrt(2*pi)
  double y = x[0];
  double sigma = p[1];
  double y0    = p[0];

  double weight = exp(-pow((y-y0)/sigma,2)/2);

  return weight;
}

// double y_shift(const double, const double, const double, const double);


// void SPSrapidity(const Char_t *name, Double_t mass, Double_t ymin, Double_t ymax, Double_t ebeam=200) {

TF1 *SPSrapidity(const Char_t *name, Double_t mass, Double_t ymin=-999, Double_t ymax = -999, Double_t ebeam=200) {

  const double mnuc  = NucleonMass;
  const double m_pi0 = pi0Mass;
  
  double sqrts        = sqrt(2.0*mnuc*ebeam+2.0*mnuc*mnuc);
  double y0           = log(2.0*ebeam/mnuc)/2.0;
  double gamma        = sqrts/(2.0*mnuc);
  double sigma_landau = sqrt(log(gamma));
  double ymax_pi0     = log(sqrts/m_pi0);
  double ymax_part    = log(sqrts/mass);
  double sig          = sigma_landau*(ymax_part/ymax_pi0);
  if (ymin == -999 or ymin < y0-ymax_part) ymin = y0-ymax_part;
  if (ymax == -999 or ymax > y0+ymax_part) ymax = y0+ymax_part;
  
  std::cout << "name         " << name         << std::endl;
  std::cout << "sqrt(s)      " << sqrts        << std::endl;
  std::cout << "y_cms        " << y0           << std::endl;
  std::cout << "gamma        " << gamma        << std::endl;
  std::cout << "sigma_Landau " << sigma_landau << std::endl;
  std::cout << "ymax_pi0     " << ymax_pi0     << std::endl;
  std::cout << "ymax_part    " << ymax_part    << std::endl;
  std::cout << "sig          " << sig          << std::endl;
  std::cout << "ymin         " << ymin         << std::endl;
  std::cout << "ymax         " << ymax         << std::endl;
  std::cout << std::endl;
 
  TF1 *fSPSrapidity = new TF1(name, SPSrapidity_Func, ymin, ymax, 2);
  fSPSrapidity->SetParameters(y0, sig);
  return fSPSrapidity;
}

  double RapidityToEta(double rapidity, double pt, double mass){

    double eta;
    float pz = sqrt(mass*mass+pt*pt) * sinh(rapidity);
    float p = sqrt(pt*pt+pz*pz);
    eta = log((p+pz)/(p-pz))/2;
    return eta;
  }

 // legacy code from EXODUS for assymetric systems
 //  or ( int ibin=1; ibin<=nbins; ibin++ )
 // {
 //   const double y  = ymin+(ibin-1)*binwidth+binwidth/2.0;
 //   double sigma;
 //   if ( y<=ymean )
 //   {
 //     sigma = sig*(1.0-((y0-ymean)*(y0-ymean)/(ymax_part-ymin_part)));
 //   }
 //   else
 //   {
 //     sigma = sig*(1.0+((y0-ymean)*(y0-ymean)/(ymax_part-ymin_part)));
 //   }
 //   double weight = std::exp(0.75*(y-ymean)/sigma)+std::exp(-0.75*(y-ymean)/sigma);
 //   weight = 1.0/(weight*weight);
 //   yhistogram->AddBinContent(ibin,weight);
 // }
 //  return yhistogram;
 

// double y_shift(const double Aproj, const double Atarg, const double bpar)
// {
//  const double density = 0.168;

//  if ( Aproj<2.0 )
//  {
//    if ( Atarg<12.0 ) return 0.0;
//    return 1.0;
//  }
//  else
//  {
//    if ( sqrt((Atarg/Aproj-1.0)*(Atarg/Aproj-1.0))<0.1 ) return 0.0;
//  }

//  double rmin = std::exp((1.0/3.0)*std::log(0.75*Aproj/TMath::Pi()/density));
//  double rmax = std::exp((1.0/3.0)*std::log(0.75*Atarg/TMath::Pi()/density));
//  if ( bpar>(rmin+rmax) ) return 0.0;

//  double hmax = 2.0*std::sqrt(rmax*rmax-rmin*rmin);

//  double t_part, p_part;
//  if ( bpar<(rmax-rmin) )
//  {
//    double rcyl   = rmin;
//    double vmax   = 2.0*TMath::Pi()*rcyl*rcyl*hmax;
//    double hmic   = rmax-std::sqrt(rmax*rmax-rmin*rmin);
//    double rapb   = (rmax-hmic)/rmax;
//      double vmin   = 2.0/3.0*TMath::Pi()*rmax*rmax*rmax*(1.0+0.5*rapb*rapb*rapb-1.5*rapb);
//      t_part = density*(2.0*vmin+vmax);
//      p_part = Aproj;
//    }
//    else
//    {
//      double rcyl   = 0.5*(rmax+rmin-bpar);
//      double vmax   = 2.0*TMath::Pi()*rcyl*rcyl*hmax;
//      double hmic   = rmax-std::sqrt(rmax*rmax-rmin*rmin);
//      double rapc   = (rmin-2.0*rcyl)/rmin;
//      t_part = density*(2.0/3.0*TMath::Pi()*rcyl*rcyl*hmic+vmax);
//      p_part = density*2.0/3.0*TMath::Pi()*rmin*rmin*rmin*
//        (1.0+0.5*rapc*rapc*rapc-1.5*rapc);
//    }
 
//    double shift = 0.0055*(t_part-p_part);
//    return shift;
//  }

//-----------------------------------------------------------------------------
//
// Generate flat distributions for any purpose
//
// default min max set to 0 - 2 pi
//-----------------------------------------------------------------------------
double Flat_Func(double *x, double *p)
{
  double phi = x[0];
  double c = p[0];

  double weight = c;

  return weight;
}

TF1 *Flat(const Char_t *name, Double_t min=0, Double_t max = 2*3.1415926) {

  double cc[0];
  cc[0]=1.;
  TF1 *fFlat = new TF1(name, Flat_Func, min, max, 1);
  fFlat->SetParameters(cc);
  return fFlat;
}


#endif
