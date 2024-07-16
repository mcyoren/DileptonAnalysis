//// Particle.h //////////////////////////////////////////////////////////////////////////////////////////////
//
// HELIOS Particle Class used for stand alone simulation in ROOT 6
//
// inherits from TLorentzVector
//
// Particle(name)  - Defines particles as Lorentz vectors, with parameters define in PDG.h. 
//                 - Defines decay branches as particle property if particle is not stable. 
//                 - Decays are handeled in Class Decay.C 
//
// currently implemented particles and decay branches by name used in HELIOS
// 
//     mesons        decays 
//     pi0           pi0->gg, pi0->gee
//     pi+           stable
//     pi-           stable
//     eta           eta->gg, eta->gee, eta->gmm, eta->mm, eta->3pi0
//     etap          etap->gg, etap->gee, etap-gmm, etap->wg, etap->2pi0eta
//     rho0          rho0->ee, rho0->mm
//     omega         omega->ee, omega->mm  omega->pi0g, omega->pi0ee, omega->pi0mm     
//     phi           phi->ee, phi->mm
//     K0s           K0s->2pi0
//     jpsi          jpsi->ee, jpsi->mm
//     psip          psip->ee, psip->mm
//
//     leptons     
//     photon        stable
//     electron      stable
//     positron      stable
//     mu+           stable
//     mu-           stable
//
//     baryons       
//     Delta         Delta->Ng
//     Nucleon       generic nucleon
//     proton        stable
//     neutron       stable
//     Lambda        Lambda->ppi-
// 
// generate random 3 vector 
// GenerateP()     - there are three options implemented through over loading 
//                 - flat, from pt spectrum, or pt, eta or y and phi spectra (see below)
//
// decays particle, three different options 
// Decay()         - generate one random decay with probability equal to branching ratio
// DecaySingleBranch(name)   - forces a particular decay
// DecayMultiBranch(name[])  - forces a slected subset of decays 
// DeacyFlat        - generates one random decay with equal probability for all defined branches
//
// member variables
//
// Double_t weight - allows to assign a weight to a given particle, e.g. dN/dpt for p generated flat
//                 - Weight() returns value
// Int_t PID       - PDG particle id, ID() returns value 
// TString Name    - name used by HELIOS - defined in PDG, Name() returns value 
// Int_t Charge    - charge of particle - Charge() returns value
//
// member functions:
//
// SetWeight(w)       - sets value of weight
// AddWeight(w)       - multiplies w to existing weight
// ResetP()           - set Lorentz vector to (0,0,0,m)
// UpdateParticle(p)  - sets TLorentz vector of particle to p
// 
// access to decay particle information:
//
//  Int_t GetNumberOfDaughters()    - returns number of daugthers generated for this particle
//  Int_t GetDaughterID(i)          - returns PDG ID of daughter i 
//  Double_t GetDaughterWeight(i)   - returns get weight for daughter, depending on method used for decay generator  
//  TLorentzVector GetDecayDaughter(i) - returns 4 vector of daughter
//  TString GetDecayName()          - returns decay name (eg. pi0->gg)
//  Int_t GetDecayNumber()          - returns decay number (eg. 0 for pi0->gg) [from DefineDecays() in Particle.C]
// 
// 11/19/2019   started by                  Axel Drees 
// 9/21/2021    integrated to HELIOS        Axel Drees
// 3/18/2022    muon decay channels added   Axel Drees
// 5/13/2022    retrieve more decay info    Roli Esha
// 6/9/2022     generate in y or eta        Axel Drees
// 11/11/2022   3 body decays and updates   Axel Drees
// 3/21/2022    add more particles & decays Axel Drees
//
#ifndef Particle_h
#define Particle_h

#include <TLorentzVector.h>                                    // MyParticle inherits from TLorentz Vector
#include <TRandom3.h>
#include "Decay.C"                                             // class of member function used for decay branches
#include "PDG.h"                                               // 

class Particle : public TLorentzVector{

public:
  Particle();                                                  // default constructor
  Particle(TString name);                                      // sets particle properties for known particles 
  virtual ~Particle();                                         // destructor

// public random generator functions for particle 4 vector 
  void GenerateP(Double_t ptmin, Double_t ptmax, Bool_t rap=true); 
                                                               // generate flat pt between min and max, -pi<phi<pi, -0.5<y<0.5
                                                               // rap=false uses pseudorapidity, if true use rapidity
  void GenerateP(Double_t ptmin, Double_t ptmax, TF1* YSpectrum, Bool_t rap=true); 
                                                               // generate flat pt between min and max, -pi<phi<pi
                                                               // uses TF1 for rapidity distribution 
                                                               // rap=false uses pseudorapidity, if true use rapidity
  void GenerateP(TF1* PtSpectrum, Bool_t rap=true);           // generates pt from ptSpectrum, -pi<phi<pi, -0.5<y<0.5
                                                               // rap=false uses pseudorapidity, if true use rapidity
  void GenerateP(TF1* Pt,TF1* PhiSpectrum, TF1* YSpectrum, Bool_t rap=true);    
                                                               // generates 4 vector from TF1 ptSpectrum, PhiSpectrum, and EtaSpectrum
                                                               // rap=false uses pseudorapidity, if true use rapidity
  void GenerateVtx();                                          // generates decay vertex location


// public random decay generators
  void Decay();                                                // generate random decay from know decay branches
  void DecaySingleBranch(TString branch);                      // generate decay for specified branch only
  void DecayMultiBranch(TString *branch, Int_t n);             // generates random decay with equal probability for
                                                               // specified branches only and sets weight corresponding weight
  void DecayFlat();                                            // generates random decay with equal probability for each 
                                                               // branch and sets weight corresponding weight
// public access to decay particles
  TString GetDecayName();                                      // get decay name (eg. pi0->gg)
  Int_t GetDecayNumber();                                      // get decay numer (eg. 0 for pi0->gg)
  Int_t GetNumberOfDaughters();                                // get number of decay particles
  Int_t GetDaughterID(Int_t index);                            // get PDG ID of decay particle index
  Double_t GetDaughterWeight(Int_t index);                     // get weight of decay 
                                                               // = 1 if random decay branch is generated with Decay() 
  TLorentzVector GetDecayDaughter(Int_t index);                // return 4 vector for decay particle index
  TVector3 GetDecayVtx();                                      // returns position of decay

  Double_t Weight();                                           // returns weight - not part of TLorentzVector
  Int_t   Charge();                                           // returns charge - not part of TLorentzVector
  TString Name();                                              // returns name   - not part of TLorentzVector
  Int_t ID();                                                  // returns ID     - not part of TLorentzVector

// manipulate particle characteristics
  void ResetP();                                               // resets momentum to 0 
  void SetWeight(Double_t w);                                  // sets particle weight
  void AddWeight(Double_t w);                                  // updates particle weight

  void UpdateParticle(TLorentzVector &p);                      // update 4 vector of particle

// Define = operator for particle class: Particle = Particle
void operator = (Particle const &otherParticle)
{
    name   = otherParticle.name;                      // set name 
    id     = otherParticle.id;                        // set id  
    charge = otherParticle.charge;                    // set charge  
    weight = otherParticle.weight;                    // set weight  
    mass   = otherParticle.mass;                      // set mass  
    stable = otherParticle.stable;                    // set stable
    double t_pt     = otherParticle.Pt();             // set 4-momentum 
    double t_eta    = otherParticle.Eta();
    double t_phi    = otherParticle.Phi();
    double t_mass   = otherParticle.M();
    SetPtEtaPhiM(t_pt,t_eta,t_phi,t_mass);
    if (!stable) DefineDecays();                      // define decays id not stable
}

// Define = operator for particle class:  Particle = TLorentzVector
void operator = ( TLorentzVector const &otherParticle)
{
    double t_pt     = otherParticle.Pt();             // set 4-momentum 
    double t_eta    = otherParticle.Eta();
    double t_phi    = otherParticle.Phi();
    double t_mass   = otherParticle.M();
    SetPtEtaPhiM(t_pt,t_eta,t_phi,t_mass);
}

 private:  
// internal variables
  TRandom3 randy = TRandom3(0);                                // Random Generator
  Bool_t debug = false;                                         // debug flag

// characteristics of particle set in constructor
  TString name;                                                // unique name 
  Double_t charge;                                             // charge 
  Double_t weight;                                             // weight 
  Double_t mass;                                               // mass in GeV
  Double_t ct;                                                 // decay length for weak decays
  Int_t id;                                                    // iD using PDG scheme
  Bool_t stable;                                               // if false decay branches will be defined in DefineDecays()
  TVector3 vtx;                                                // vertex or origin of particle, default 0,0,0 unless from weak decay

// characteristics of decay Branches and decay particles 
  void DefineDecays();                                         // defines decay of parent particle if not stable 
  void ExecuteDecay(Int_t branch, Double_t ww, Int_t i=0);     // manages random decay after branch is selected
  class Decay DecayBranch[10];                                 // possible decay branches - maximum of 10  
  Int_t NumberOfBranches = 0;                                  // actual number of decay branches 
  TString DecayName = "null";                                  // name of the decay channel
  Int_t DecayNumber = -999;                                    // branch number of the decay channel
  Int_t NumberOfDaughters = 0;
  TLorentzVector daughter[10][10];                             // 4 vector for daughters [branch][index], max 10 per branch
  Int_t daughterID[10][10];                                    // ID of daughters [branch][index]
  TLorentzVector DecayDaughter[100];                           // 4 vector for daughters [branch][index], max 10 per branch
  Int_t DecayDaughterID[100];                                  // PDG particle ID of daughter 
  Double_t DecayDaughterWeight[100];                           // weigth of decay branch that produced daughter 
  
  Double_t pi = 3.14159;                              // define pi

}; // end of Particle class
//
// End of Header File
//

#endif 


