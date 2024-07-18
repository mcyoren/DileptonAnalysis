////// DecayAuxilliaryFunctions.h ///////////////////////////////////////////////////////////////////////////////////////
//
// Auxiliary functions used in Class Decays  
//
// mass distributions functions use to define TF1's for Dalitz Decays: 
//
// f_pi0Dalitz           ee Dalitz
// f_etaDalitz 
// f_etapDalitz 
// f_omegaDalitz 
// f_etaDalitz2          mumu Dalitz
// f_etapDalitz2 
// f_omegaDalitz2 
//
// mass distributions functions used to define TF1's for 2 body decays with finite width:
//
// f_omegaee
// f_rho0ee
// f_phiee
// f_omegamm
// f_rho0mm
// f_phimm
//
// Support functions:
// 
// KrollWada:        is a generalized Kroll Wada Function for p -> o e+e-, where p is the parent 
//                   and o is the "other" particle. For o being a photon this will be the well 
//                   known Kroll Wada distribution. 
//
// BreitWigner:      Breit-Wigner parameterisation of resonances.
//
// GounarisSakurai:  Generalized Breit-Wigner function, commonly used for all 2 body decays 
//                   with finite width. 
//
// EMTFormFactor:    Lepton G parameterization of Electro Magnetic Tranition Form Factor in 
//                   Dalitz decays.
//
// Axel Drees 8/30/2021
// expanded   9/9/2021    
// muon decays added 3/18/2022
//  
#ifndef Aux_h
#define Aux_h

#include "PDG.h"

// generalized Kroll Wada function 
double KrollWada(const double mll, const double mp, const double ml, const double mo=0)
//
// mll  = lepton pair mass
// mp   = parent mass
// ml   = lepton mass
// mo   = other decay particle mass - default 0 for photon
{
  Double_t val = 0;
  Double_t c1 = 2*alpha/3/Pi;

  Double_t q = pow(mll/mp,2);
  Double_t eps = pow(ml/mp,2);
  Double_t delta = pow(mo/mp,2);

  Double_t f = (1+q/(1-delta))*(1+q/(1-delta)) - 4*q/(1-delta)/(1-delta); 
  val = 2*c1/mll *  pow(f,3/2) * sqrt(1-4*eps/q) * (1+2*eps/q);
  return val;
}
// short lived resonances
double GounarisSakurai(const double mll, const double mp, const double Gp, const double ml)
{ 
  Double_t value = 0;
  Double_t eps = pow(ml/mll,2);
  Double_t delta = pow(pi0Mass/mll,2);

  Double_t Gamma = Gp*mp/mll * pow((mll*mll/4-ml*ml)/(mp*mp/4-ml),3/2); 

  value = sqrt(1-4.*eps)*(1+2.*eps)*pow((1-4.*delta),3/2);
  value = value / (pow((mp*mp-mll*mll),2) + mp*mp*Gamma*Gamma);
  return value;

  }
double BreitWigner(const double mll, const double mp, const double Gp, const double ml)
{ 
  Double_t value = 0;
  value = Gp / ( pow(mll-mp ,2) + Gp*Gp/4.);
//  std::cout << mll << "   " << Gp << "  "  << mp << "   " << pow(mll-mp,2) << "  " << Gp*Gp/4 << "   " << value << std::endl;
  return value;
}

// Electro magnetic Transition Form Factor
double EMTFormFactor(const double mee, const double b)
{
  Double_t value = pow(1.0/(1.0-b*mee*mee),2);
  return value;

}
// Mass distributions used in Decays to define TF1's  
//
// Dalitz Decay functions for individual decays including formfactors
double f_pi0Dalitz(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,pi0Mass,eMass) * EMTFormFactor(mll,b_pi0_Dalitz);
  return y; 
}
//
double f_etaDalitz(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,etaMass,eMass) * EMTFormFactor(mll,b_eta_Dalitz);
  return y; 
}
double f_etaDalitz2(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,etaMass,muMass) * EMTFormFactor(mll,b_eta_Dalitz);
  return y; 
}
//
// taken from EXODUS code
double f_etapDalitz(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,etapMass,eMass) * pow(0.764, 4) 
             / (pow( pow(0.764, 2) - mll*mll, 2) + pow(0.1020*0.764, 2));
  return y; 
}
double f_etapDalitz2(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,etapMass,muMass) * pow(0.764, 4) 
             / (pow( pow(0.764, 2) - mll*mll, 2) + pow(0.1020*0.764, 2));
  return y; 
}
//
double f_omegaDalitz(Double_t *x, Double_t *p) {
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,omegaMass,eMass,pi0Mass) *  EMTFormFactor(mll,b_omega_pi0ll);
  return y; 
}
double f_omegaDalitz2(Double_t *x, Double_t *p) {
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,omegaMass,muMass,pi0Mass) *  EMTFormFactor(mll,b_omega_pi0ll);
  return y; 
}
double f_omegaee(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//  Double_t y = BreitWigner(x[0],omegaMass,omegaWidth,eMass);
  Double_t y = GounarisSakurai(mll,omegaMass,omegaWidth,eMass);
  return y; 
}
double f_omegamm(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//  Double_t y = BreitWigner(x[0],omegaMass,omegaWidth,eMass);
  Double_t y = GounarisSakurai(mll,omegaMass,omegaWidth,muMass);
  return y; 
}
double f_rho0ee(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//  Double_t y = BreitWigner(x[0],omegaMass,omegaWidth,eMass);
  Double_t y = GounarisSakurai(mll,rho0Mass,rho0Width,eMass) 
               * pow(mll*T_rho0_ll,3/2) * exp(-mll/T_rho0_ll);
  return y; 
}
double f_rho0mm(Double_t *x, Double_t *p) {
  Double_t mll=x[0];
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//  Double_t y = BreitWigner(x[0],omegaMass,omegaWidth,eMass);
  Double_t y = GounarisSakurai(mll,rho0Mass,rho0Width,muMass) 
               * pow(mll*T_rho0_ll,3/2) * exp(-mll/T_rho0_ll);
  return y; 
}
double f_phiee(Double_t *x, Double_t *p) {
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//  Double_t y = BreitWigner(x[0],omegaMass,omegaWidth,eMass);
  Double_t y = GounarisSakurai(x[0],phiMass,phiWidth,eMass);
  return y; 
}
double f_phimm(Double_t *x, Double_t *p) {
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//  Double_t y = BreitWigner(x[0],omegaMass,omegaWidth,eMass);
  Double_t y = GounarisSakurai(x[0],phiMass,phiWidth,muMass);
  return y; 
}
double f_etapRho0g(Double_t *x, Double_t *p) {
  Double_t y = 0;
  return y;
}

double f_phietaee(Double_t *x, Double_t *p) {
//  std::cout << " ----------- " << x[0] << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  Double_t mll=x[0];
  Double_t y = KrollWada(mll,phiMass,eMass,etaMass);
  return y; 
}


#endif 
