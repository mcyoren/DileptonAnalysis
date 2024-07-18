//////   PHENIXDetector.h //////////////////////////////////////////////////////////////////
//
// HELIOS PHENIXDetector Class used to approximate acceptance and detector response 
// of the main PHENIX central arm detectors, EMCal and DC/PC1
// 
//
// PHENIXDetector()     constructor, sets private member variables not specific to a particle   
//
// high level member functions to generated PHENIX "default" response 
// - ReconstructParticle(Particle)            One stop shop ... works with Particle Class
// - CharacterizeTrack(particle, charge, id)  sets private member variables for particle
// - ReconstructTrack(particle, id)           returns TLorentzVector with reconstructed charged track
// - ReconstructTrack(particle, id, r)        returns TLorentzVector with reconstructed charged track
//                                                    assumes off vertex 
// - ReconstructShower(particle, id)          returns TLorentzVector with reconstrcuted EMCal shower 
// 
// more specific public member functions that can be used for systematic studies by 
// variing default response, more detailed description can be found in PHENIXDetector.C
// 
// general:
// - InAcceptance(particle,charge)            Ideal PHENIX acceptance - checks if particle is in acceptance returns arm
// - Acceptance()                             used in InAcceptance() 
// EMCal
// - EMCalSector(phi)                         returns sector number 
// - EMCalSectorCoordinates(phi,theta, y, z)  returns x,z in local sector coordinates 
// - EMCalLive(sector)                        returns true or false, implemented statistically without dead
// - SmearEnergy(energy, opt, c1, c2)         returns smeared energy, default parameters opt 0=PbSc, 1=PbGl,  
//                                            all other values use provided c1, c2
// - SmearEMCalPosition(energy, x, opt)       returns smeared position, 
//                                            x is distance from sector center, opt 0=PbSc, 1=PbGl
// - ShowerEnergyCorrection(energy, sinT)     returns energy correction based on impact angle 
// - ShowerEnergyCorrection(energy, y, z)     returns energy correction based on impact angle (for photons only)
// - EMCalImpactAngle(particle, id)           returns sinT impact angle on calorimeter 
// - EMCalPhi(particle, q)                    returns phi angle of charged particle at calorimeter  
// - NonLinearEnergy(energy, c0, c1, c2)      optional nonlinearty with parameters c0=1.003, c1=0.05, c2=1.77
//                                            not used in default shower reconstruction
// DC/PC1
// - SmearPt(pt, q)                           returns smeared pt, charge may change at very high momenta
// - SmearDCphi0(pt, phi0)                    returns smeared phi             
// - SmearDCeta(pt, eta0)                     returns smeared eta
// - DCPhi(particle, q)                       returns phi angle at DC 
// - RICHPhi(particle, q)                     dito at RICH 
// other functions
// - ElectronMomentum(E, p)                   calculates weighted average of electron momentum from reconstructed 
//                                            E and p   
// - VTXConversion()                          returns 1-4 if a conversion in VTX layer will occure, else 0 
//
// member functions that return private member variables, can be accessed after creating defaul response to particle 
// - Arm()
// - Sector()
// - SectorY()
// - SectorZ()
// - SectorSinT()
// - Phi_EMCal()
// - Phi_DC()
// - Phi_DC()
// - Phi_RICH()
//
// private members
//
// internal variables
//   TRandom3 randy = TRandom3(0);                                // Random Generator
//   Double_t SectorX0[8],SectorXmax[8],SectorXmin[8];
//   Double_t SectorY0[8],SectorYmax[8],SectorYmin[8];
//   Double_t P_R_EMCal_max[8];
//   Bool_t Intersect(Double_t &Px,Double_t &Py,Double_t *x, Double_t *y);
//   Double_t  phi_EMCal;
//   Double_t  theta_EMCal;
//   Double_t  phi_RICH;
//   Double_t  phi_DC;
//   Double_t  pt_DC;
//   Double_t  eta_DC;
//   Double_t  q_charge;
//   Double_t  alpha_DC;
//   Double_t  alpha_EMCal;
//   Double_t  sectorZ;
//   Double_t  sectorY;
//   Double_t  sectorSinT;
//   Int_t     sector;
//   Int_t     arm; 
//
// Axel Drees
// updated 11/22/2022
//
///////////////////////////////////////////////////////////////////////////////////

#ifndef PHENIX_h
#define PHENIX_h

#include <TLorentzVector.h>                                    // MyParticle inherits from TLorentz Vector
#include <TRandom3.h>
#include <TMath.h>
#include "PDG.h"                                             // 
#include "PHENIXSetup.h"
#include "InteractionWithMaterial.h"

class PHENIXDetector
{

public:
  PHENIXDetector();                                                // default constructor
  virtual ~PHENIXDetector();                                       // destructor

// high level interface to PHENIX stand alone fast simulation
  TLorentzVector ReconstructParticleShower(Int_t ID, TLorentzVector particle, TVector3 vtx2); 

// EMCal acceptance and local variables
  Bool_t InEMCalAcceptance(Int_t ID, TLorentzVector particle, TVector3 vtx2);

// apply energy resolution etc .. 
  TLorentzVector ReconstructShower(Int_t id);

// LEGACY stuff
// member functions - see definitions below for details 
// general central arm functions  
  Int_t InAcceptance(TLorentzVector particle, Int_t q);            // standalone acceptance calculation
  Int_t Acceptance();                                              // same as InAcceptance() but expects 
                                                                   // private member variables preset
// calculate inportant variables for given particle
  void CharacterizeTrack(TLorentzVector VT, Int_t q, Int_t id);

// EMCal 
  TLorentzVector ReconstructShower(TLorentzVector VT, Int_t id);
  Int_t    EMCalSector(Double_t phi);
  Int_t    EMCalSectorCoordinates(Double_t phi, Double_t theta, Double_t& y, Double_t& z);
  Bool_t   EMCalLive(Int_t sector);
  Bool_t   EMCalEff(Int_t sector);
  Double_t NonLinearEnergy(Double_t energy, Double_t c0=1.003, Double_t c1=0.05, Double_t c2=1.77);
  Double_t SmearEnergy(Double_t energy, Int_t opt=0, Double_t c1=0.081, Double_t c2=0.021 );
  Double_t SmearEMCalPosition(Double_t energy, Double_t x, Int_t opt=0 );
  Double_t ShowerEnergyCorrection(Double_t energy, Double_t sinT);
  Double_t ShowerEnergyCorrection(Double_t energy, Double_t y, Double_t z);
  Double_t EMCalImpactAngle(TLorentzVector VT, Int_t id);
  Double_t EMCalPhi(TLorentzVector VT, Double_t q); 
  Double_t EMCalPhi(TLorentzVector VT, TVector3 vtx); 

// Tracking
  TLorentzVector ReconstructTrack(TLorentzVector VT, Int_t id);
//  TLorentzVector ReconstructTrack(TLorentzVector VT, Int_t id, Double_t r=0);
  Double_t SmearPt(Double_t pt, Int_t &q);
  Double_t SmearDCphi0(Double_t pt, Double_t phi0);
  Double_t SmearDCeta(Double_t pt, Double_t eta0);
  Double_t DCPhi(TLorentzVector VT, Double_t q); 
  Double_t RICHPhi(TLorentzVector VT, Double_t q); 

  Double_t ElectronMomentum(Double_t E, Double_t p);

// physics processes 
  Bool_t VTXConversion(TVector3 &vtx2);

// get Private variables
  Int_t Arm(){
    return arm;
  }
  Int_t Sector(){
    return sector;
  }
  Double_t SectorY(){
    return sectorY;
  }
  Double_t SectorZ(){
    return sectorZ;
  }
  Double_t SectorSinT(){
//   if (sector != 0) std::cout << " --+-- " << sectorSinT << std::endl;
    return sectorSinT;
  }
  Double_t Phi_EMCal(){
    return phi_EMCal;
  }
  Double_t Phi_DC(){
    return phi_DC;
  }
  Double_t Phi_RICH(){
    return phi_RICH;
  }

private:

// internal variables
  TRandom3 randy = TRandom3(0);                                // Random Generator
  Bool_t debug = true;

  Double_t SectorX0[8],SectorXmax[8],SectorXmin[8];
  Double_t SectorY0[8],SectorYmax[8],SectorYmin[8];
  Double_t P_R_EMCal_max[8];
  Bool_t Intersect(Double_t &Px,Double_t &Py,Double_t *x, Double_t *y);
  Double_t  E_EMCal;
  Double_t  phi_EMCal;
  Double_t  theta_EMCal;
  Double_t  phi_RICH;
  Double_t  phi_DC;
  Double_t  pt_DC;
  Double_t  eta_DC;
  Double_t  q_charge;
  Double_t  alpha_DC;
  Double_t  alpha_EMCal;
  Double_t  sectorZ;
  Double_t  sectorY;
  Double_t  sectorSinT;
  Int_t     sector;
  Int_t     arm; 
  
  Double_t pi = 3.14159;                              // define pi

};

#endif

