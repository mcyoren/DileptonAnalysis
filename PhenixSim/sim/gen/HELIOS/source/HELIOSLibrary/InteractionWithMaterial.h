//// InteractionWithMaterial.h //////////////////////////////////////////////////////////////////////////
//
// collection of functions to simulate interaction of particles with material
//
// currently implemented:
//
//   PhotonConversion
//   
//
// Axel Drees  9/22/2021 
//

#ifndef IntMat_h
#define IntMat_h

#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

// used for photon conversions 
Double_t EnergySplit(const Double_t *x, const Double_t *p){
   Double_t value;
   value = 1-4./3.*x[0]*(1-x[0]);
   return value;
}

// energy distribution of brems strahlungs photon in deltaE/E (relative to electron energy) 
// see PDG 2002 chapter 26 for this approximation
Double_t BremsEnergy(const Double_t *x, const Double_t *p){
  Double_t value;
  value = 1./x[0]*(4./3.*(1-x[0])+x[0]*x[0]);
  return value;
}

// average number of bremsstrahlungs photons
Double_t BremsPhotons(Double_t xX0, Double_t kmin, Double_t kmax) {
  Double_t value;
//  xX0 = fraction of radiation length 
//  kmin = lower energy cutoff
//  kmax = upper energy cutoff
  value = xX0*(4./3.*log(kmax/kmin)-4./3.*(kmax-kmin) + (kmax*kmax-kmin*kmin)/2);
  return value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Bethe Bloch equation for average energy loss
//
// Double_t BetheBloch(Double_t x){
// }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PhotonConvert(TLorentzVector &photon, TLorentzVector &electron, TLorentzVector &positron){
//
// Converts photon to e+e- in high energy limit
//
// input:  TF1 ConversionEnergySplit      TF1 used to generate electron energy set outside using EnergySplit defined below
//         photon                         4 vector
//         decay particle 1&2             4 vector does not need to be set
//         
// 
// output: parent and decay    4 vectors
//
// Axel Drees 10/19/2018 
//
  Double_t pi = 3.14159;                              // define pi
  Double_t px,py,pz,E;                                // define generic 4 vector component 
  Double_t me = 0.000511;                             // electron mass in GeV
  Double_t e0 = me/photon.E();                        // minimum energy is leectron mass 
//
//
  if (gROOT->FindObject("ConvSplit") == NULL){
    TF1 *ConversionEnergySplit = new TF1("ConvSplit",EnergySplit,0.,1.,0.);                                                      // function that gives energy split for conversion
  }
  TF1 *h = (TF1 *) gROOT->FindObject("ConvSplit");    // get eepair mass
  Double_t phi   = photon.Phi();                      // phi of e+e- same as photon
  Double_t theta = photon.Theta();                    // theta same as photon 
  Double_t e = h->GetRandom();                        // energy of electron - assume E=p 

  E = e*photon.E();
  pz = E*cos(theta);                                  // longitudinal momentum
  px = E*sin(theta)*cos(phi);                         // transverse momentum x component
  py = E*sin(theta)*sin(phi);                         // transverse momentum y component
  electron.SetPxPyPzE(px,py,pz,E);                    // set 4 vector of electron
  E = (1-e)*photon.E();
  pz = E*cos(theta);                                  // longitudinal momentum
  px = E*sin(theta)*cos(phi);                         // transverse momentum x component
  py = E*sin(theta)*sin(phi);                         // transverse momentum y component
  positron.SetPxPyPzE(px,py,pz,E);                    // set 4 vector of electron

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t BremsEnergyLoss(Double_t xX0) {
//
//
//  
  TRandom3 randy = TRandom3(0);                                // Random Generator
  const Double_t kmin = 0.005;
  const Double_t kmax = 0.99; 
  Double_t BE, BEtot=0; 

  BE = 0;
  BEtot = 0;
  if (gROOT->FindObject("BremsE") == NULL){
    TF1 *BremsE = new TF1("BremsE",BremsEnergy,kmin,kmax,0.);                                                      // function that gives energy split for conversion
  }
  TF1 *h = (TF1 *) gROOT->FindObject("BremsE");    // get bremsstrahlungs distribution

  Double_t AvgBrems = BremsPhotons(xX0, kmin, kmax);
  Int_t Nbrems = randy.Poisson(AvgBrems);
  //  std::cout << "avg number of brems photons  " << AvgBrems << " number of brems photons " << Nbrems << std::endl;

  for (int i=0; i<Nbrems; i++){
     BE = h->GetRandom();
      // std::cout << "energy of photon " << i << "  " << BE << std::endl;
     BEtot = +BE;
  } 
  return BEtot; 
}


#endif
