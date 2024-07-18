////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Particle Properties used in HELIOS from PDG: https://pdg.lbl.gov/2021
//
//
// Axel Drees 8/30/2021
//            9/21/2021  updated
//            3/13/2021  muon decays added
// Roli Esha  5/13/2022  added photonMass = 0
// Axel Drees 11/18/2022 added several decays 
// Axel Drees 3/21/2023  added J/psi and psi(2s), and previously Lamda and K0s
//  
#ifndef PDG_h
#define PDG_h



// constants
  Double_t Pi = 3.14159;                              // define pi
  Double_t alpha = 1/137.;                            // fine structure constant

// masses of defined particles
  Double_t photonMass = 0;                             // photon
  Double_t pi0Mass  = 0.134977;                        // pi0
  Double_t piMass   = 0.1395702;                       // pi+/pi-
  Double_t etaMass  = 0.54786;                         // eta
  Double_t eMass    = 0.000511;                        // e+/e-
  Double_t muMass   = 0.10566;                         // mu+/mu-
  Double_t rho0Mass = 0.77526;                         // rho meson
  Double_t omegaMass= 0.78266;                         // omega meson
  Double_t etapMass = 0.95778;                         // eta' meson
  Double_t phiMass  = 1.0195;                          // phi meson
  Double_t K0sMass  = 0.497611;                        // K0 short 
  Double_t jpsiMass = 3.0969;                          // J/Psi
  Double_t psipMass = 3.6861;                          // psi(2s)


// baryons
  Double_t DeltaMass   = 1.232;                          // Delta baryon average
  Double_t NucleonMass = 0.939;                        // Nucleon (p+n)/2
  Double_t protonMass  = 0.93827;
  Double_t neutronMass = 1.00866;
  Double_t LambdaMass  = 1.11568;                        // lightest starnge baryon

// width of particles
  Double_t omegaWidth = 0.00849;  
  Double_t rho0Width  = 0.1491;   
  Double_t phiWidth   = 0.0042;

//  Double_t DeltaWidth = 0.117;
  
// PDG group particle IDs
  Int_t pi0ID      =  111;
  Int_t pipID      =  211;
  Int_t pimID      = -211;
  Int_t etaID      =  221;
  Int_t photonID   =  22;
  Int_t virtualg   = -22;
  Int_t electronID =  11;
  Int_t positronID = -11;
  Int_t mupID      =  13; 
  Int_t mumID      = -13; 
  Int_t omegaID    = 223;
  Int_t rho0ID     = 113;
  Int_t etapID     = 331;
  Int_t phiID      = 333;
  Int_t K0sID      = 310;
  Int_t DeltaID    = 2114;                            // this is the ID of the Delta_0
  Int_t NucleonID  = 2112;                            // neutron ID  
  Int_t protonID   = 2212;
  Int_t neutronID  = 2112;
  Int_t LambdaID   = 3122;
  Int_t jpsiID     = 443;
  Int_t psipID     = 100443;

// PDG decay lengths for weak decays
  Double_t K0s_ct   =  2.6844;                        // ct in cm
  Double_t Lambda_ct = 7.89; 


// Branching ratios of defined decays
  Double_t BR_pi0_gg        = 0.98823; 
  Double_t BR_pi0_Dalitz    = 1-0.98823;
  Double_t BR_eta_gg        = 0.3941;
  Double_t BR_eta_Dalitz    = 0.0069;
  Double_t BR_eta_Dalitz2   = 0.0031;                 // g mu mu
  Double_t BR_eta_mm        = 0.000058;
  Double_t BR_eta_3pi0      = 0.3257;
  Double_t BR_omega_pi0g    = 0.0834;
  Double_t BR_omega_pi0ee   = 0.00077;
  Double_t BR_omega_ee      = 0.0000739;
  Double_t BR_omega_pi0mm   = 0.000134;
  Double_t BR_omega_mm      = 0.000074;
  Double_t BR_etap_gg       = 0.0231;
  Double_t BR_etap_Dalitz   = 0.000491;
  Double_t BR_etap_Dalitz2  = 0.000113;
  Double_t BR_etap_wg       = 0.0252;
  Double_t BR_etap_rho0g    = 0.295;
  Double_t BR_etap_2pi0eta  = 0.224;
  Double_t BR_rho0_ee       = 0.0000472;
  Double_t BR_rho0_mm       = 0.000045;
  Double_t BR_phi_ee        = 0.000297;
  Double_t BR_phi_etaee        = 0.000108;
  Double_t BR_phi_mm        = 0.000286;
  Double_t BR_K0s_2pi0      = 0.314 ;
  Double_t BR_Delta_Ng      = 0.006;
  Double_t BR_Lambda_ppim   = 0.639;
  Double_t BR_jpsi_ee       = 0.05971;
  Double_t BR_jpsi_mm       = 0.05961;
  Double_t BR_psip_ee       = 0.00793; 
  Double_t BR_psip_mm       = 0.008; 

// transition formfactor slope values -> NA60 and Lepton G
  Double_t b_pi0_Dalitz = 5.5;
  Double_t b_eta_Dalitz = 1.934;
  Double_t b_etap_Dalitz;
  Double_t b_omega_pi0ll = 2.223; 

// effective temperature defining the line shape of rho -> ll (NA60)
  Double_t T_rho0_ll = 0.161; 

// particle ratios at pt->inft. (can be used as weights for mt-scaled pt distributions)
  Double_t Eta_to_Pi0 = 0.487; // 0.487±0.024 from Yuanjee Ren Thesis
  Double_t Etap_to_Pi0 =  0.25; // 0.25±0.075 from Wenqing Fan Thesis table 4.4
  Double_t Omega_to_Pi0 =  0.81; // 0.81+/-0.02(stat)+/-0.09(sys) from ppg118 
  Double_t Rho_to_Pi0 = 1.15; 
  Double_t Phi_to_Pi0 = 0.25;
  Double_t Jpsi_to_Pi0 = 0.25;
  Double_t Psip_to_Pi0 = 0.025;


  Int_t PDG_Charge(Int_t ID){

  Int_t q = 0; 

  if (ID == pipID) q = 1;
  if (ID == pimID) q = -1;
  if (ID == electronID) q = -1;
  if (ID == positronID) q = 1;
  if (ID == mupID) q = 1;
  if (ID == mumID) q = -1;
  if (ID == protonID) q = 1;
  
  return q;
}


#endif 
