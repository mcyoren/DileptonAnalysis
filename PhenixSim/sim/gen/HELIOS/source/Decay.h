////// Decay.h //////////////////////////////////////////////////////////////////////////////////////////////////
//
// Header for Class Decay that handels all decay branches defined for particles in HELIOS
//
// currently defined decay names (selfexplanatory) used to identify the correct parent specific decay branch 
//
//     "pi0->gg"       "pi0->gee"
//     "eta->gg"       "eta->gee"         "eta->gmm"      "eta->mm"
//     "etap->gg"      "etap->gee"        "etap->rho0g"   "etap->gmm"        "etap->wg"       "etap->2pi0eta"
//     "omega->pi0g"   "omega->pi0ee"     "omega->ee"     "omega->pi0mm"     "omega->mm"
//     "rho0->ee"      "rho0->mm"
//     "phi->ee"       "phi->mm" 
//     "jpsi->ee"      "jpsi->mm" 
//     "psip->ee"      "psip->mm" 
//     "K0s->2pi0"
//     "Delta->Ng"
// 
// currently defined decay types of decay branches used internaly to select correct decay function                                   
//     "TwoBody"    two body decay 
//     "Dalitz"     Dalitz decay
//     "ThreeBody"  three body decay
//
// Particle class defines each decay branch for non stable particle with Decay constructor and sets 
// properties using:
//   SetName("name")                   as defined above
//   SetType("type")                   as defined above
//   SetBR(BR)                         branching ratio of decay 
//   SetNumberOfDecayParticles(n)      # of decay particles
//   SetMassDistributions()            of lepton pairs from Dalitz or vector meson decays 
//   DefineDaughters(*4vector)         preset array of TLorentz vectors of decay particles  
//
// random decay particle generation for defined branch:
//   GenerateDecay()                   produces random decay particles for this branch boosted in to parent frame
//   SetParent(pt,eta,phi,m)           used to update parent particle prior to decay simulation
//
// decay particles and properties can be accessed via:
//  GetNumberDecayParticles()    returns number of decay particles for this decay branch
//  GetDaughter(i)               returns TLorentzVector of daughter i    
//  GetBR()                      returns branching ratio
//  GetName()                    returns name of decay branch
//  GetType()                    returns type of decay branch
//
// defined operators:    " = "   set one decay branch equal to another
//
//
//  8/23/2021     originally developed Axel Drees 
//  9/21/2021     addpted to HELIOS   
//  3/18/2022     muon decay channels added 
//  11/13/2022    3 body decays           
//

#ifndef Decay_h
#define Decay_h

#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TF1.h>
#include "PDG.h"
#include "DecayAuxilliaryFunctions.h"

class Decay {
 public:
                                                           
  Decay();                                                  // constructor for decay                 
  virtual ~Decay();                                         // destructor

// access to member variables 
  TLorentzVector GetDaughter(Int_t index){                  // returns daughter 4 vector from list
    return Daughter[index];
  }
  Int_t GetNumberDecayParticles() {                         // returns number of decay particles generated
    return NumberDecayParticles;
  }
  TString GetType() {                                       // return type of decay
    return DecayType;
  }
  TString GetName() {                                       // return name of decay
    return DecayName;
  }
  Double_t GetBR() {                                        // returns branching ratio of decay
    return BR;
  }

// functions allowing to set some member variables          
  void SetBR(Double_t br){                                  // set branching ratios
    BR = br;
  }
  void SetType(TString t){                                  // set decay time  
    if (t == "TwoBody" or t == "Dalitz" or t == "ThreeBody") {
      DecayType = t;
    } else {
      std::cout << " decay type " << t << " not known" << std::endl;
      DecayType = "unknown";
    }
  }
  void SetName(TString n){
    if (   n == "pi0->gg" or n == "pi0->gee" 
        or n == "eta->gg" or n == "eta->gee" or "eta->gmm" or "eta->mm" or "eta->3pi0"
        or n == "omega->pi0g"  or n == "omega->pi0ee" or n == "omega->ee" or n == "omega->gmm"
        or n == "rho0->ee" or n == "rho0->mm"
        or n == "etap->gg" or n == "etap->gee" or n == "etap->rho0g" or n == "etap->wg" or n == "etap->gmm" or n == "etap->2pieta"
        or n == "phi->ee" or  n == "phi->mm" or n == "phi->etaee"
        or n == "jpsi->ee" or  n == "jspi->mm"
        or n == "psip->ee" or  n == "psip->mm"
        or n == "Delta->Ng") {
      DecayName = n;
    } else {
      std::cout << " decay type " << n << " not known" << std::endl;
      DecayName = "unknown";
    }
  }
  void SetNumberOfDecayParticles(Int_t n){                   // set number of decay particles
    NumberDecayParticles = n;
  }
  void DefineDaughters(TLorentzVector *daughter){            // set daughter 4 vectors
    for (int i=0; i<NumberDecayParticles; i++){
      Daughter[i] = daughter[i];
    }
  }                                                         // set parent kinematics
  void SetParent(Double_t pt, Double_t phi, Double_t eta, Double_t m){
    Parent.SetPtEtaPhiM(pt,eta,phi,m);
  }
  void SetMassDistributions();         // used for Dalitz decays and decays with finite width 


// define = operator for this class
  void operator=(Decay &D){
   DecayType = D.DecayType;
   DecayName = D.DecayName;
   NumberDecayParticles = D.NumberDecayParticles;
   for (Int_t i; i<10;i++){
    Daughter[i] = D.Daughter[i];
   }
   BR = D.BR;
//   std::cout << " created Decay through = operator " << std::endl;
  }

// generator of decay  
  void GenerateDecay();

 private:
// core member variables for Deacy class
   TString DecayType;
   TString DecayName;
   Int_t NumberDecayParticles;
   TLorentzVector Daughter[10];
   TLorentzVector Parent;
   Double_t BR;
   Bool_t debug = false;

  Double_t pi = 3.14159;                              // define pi
// defined decay member functions
   void TwoBodyDecay(TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2);
   void DalitzDecay (TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2, TLorentzVector &daughter3);
   void ThreeBodyDecay(TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2, TLorentzVector &daughter3);

}; //end of class


#endif 
