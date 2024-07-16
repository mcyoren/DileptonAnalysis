////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// See Particle.h for detailed discription 
//
// Particle properties are set with constructor, only particle know to HELIOS can be defined
// 
//
// Axel Drees 11/19/2019
// addaped to HELIOS 9/21/2021
// updated to write ROOT and Oscar output 05/13/2022 

#include "Particle.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Particle Constructors
//
// Particle has 
//                  4Vector intialized with momntum = 0
//                  mass, ID, name, default weight = 1.
//                  
// Known particles:  uses PDG convention for ID
//                   see Particle.h for details
//
// Axel Drees      11/19/2019
// last updated    3/21/2023 see Partile.h 
// 
//

  Particle::Particle(){                 // default constructor makes pi0
    // name   = "pi0";
    // id     = pi0ID;
    // charge = 0;
    // weight = 1.;
    // mass   = pi0Mass;
    // stable = true;
    // SetPxPyPzE(0.,0.,0.,mass);
  }

  Particle::Particle(TString n){        // sets particle properties for known particles    
    name = n;

    ct     = 0;                             // ct is zero unless explicitly set
    if (name == "pi0"){                     // neutral pion
      id     = pi0ID;
      charge = 0;
      weight = 1.;
      mass   = pi0Mass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name=="pi+"){                // pi +
      id     = pipID;
      charge = 1;
      weight = 1.;
      mass   = piMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name == "pi-"){              // pi -
      id     = pimID;
      charge = -1;
      weight = 1.;
      mass   = piMass; 
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name =="eta") {              // eta 
      id = etaID;
      charge = 0;
      weight = 1.;
      mass = etaMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name == "photon") {          // gamma
      id = photonID;
      charge = 0;
      weight = 1.;
      mass = 0;
      stable = true;
      SetPxPyPzE(0.,0.,0.,0.);
    } else if (name == "electron") {        // electron
      id = electronID;
      charge = -1;
      weight = 1.;
      mass = eMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name =="positron") {         // positron
      id = positronID;
      charge = +1;
      weight = 1.;
      mass = eMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name =="mu+") {             // mu+
      id = mupID;
      charge = +1;
      weight = 1.;
      mass = muMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name =="mu-") {             // mu+
      id = mumID;
      charge = +1;
      weight = 1.;
      mass = muMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);
    } else if (name == "omega") {
      id = omegaID;
      charge = 0;
      weight = 1.;
      mass = omegaMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if (name == "rho0") {
      id = rho0ID;
      charge = 0;
      weight = 1.;
      mass = rho0Mass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if (name == "etap") {
      id = etapID;
      charge = 0;
      weight = 1.;
      mass = etapMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if (name == "phi") {
      id = phiID;
      charge = 0;
      weight = 1.;
      mass = phiMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if (name == "K0s") {
      id = K0sID;
      charge = 0;
      weight = 1.;
      mass = K0sMass;
      ct   = K0s_ct;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "Delta"){        // this is an "average" Delta baryon 
      id = DeltaID;
      charge = 0;
      weight = 1.;
      mass = DeltaMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "Nucleon"){      // this is an "average" Nucleon baryon 
      id = NucleonID;
      charge = 0;
      weight = 1.;
      mass = NucleonMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "proton"){      
      id = protonID;
      charge = 1;
      weight = 1.;
      mass = protonMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "neutron"){      
      id = neutronID;
      charge = 0;
      weight = 1.;
      mass = neutronMass;
      stable = true;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "Lambda"){      
      id = LambdaID;
      charge = 0;
      weight = 1.;
      mass = LambdaMass;
      stable = false;
      ct   = Lambda_ct;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "jpsi"){      
      id = jpsiID;
      charge = 0;
      weight = 1.;
      mass = jpsiMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    
    } else if ( name == "psip"){      
      id = psipID;
      charge = 0;
      weight = 1.;
      mass = psipMass;
      stable = false;
      SetPxPyPzE(0.,0.,0.,mass);    


    } else {
     std::cout << name << " particle is not defined" << std::endl;
    }

    if (debug) std::cout << " particle " << name << " is being created with mass " << mass << " ID " << id << std::endl;

    if (!stable){
       DefineDecays();
    }


  }
  Particle::~Particle(void){}


  // public access to decay particles
  TString Particle::GetDecayName(){                                // get name of the decay channel
    return DecayName;
  }

  Int_t Particle::GetDecayNumber(){                                // get branch_number of the decay channel
    return DecayNumber;
  }

  Int_t Particle::GetNumberOfDaughters(){                                // get number of decay particles
    return NumberOfDaughters;
  }
  Int_t Particle::GetDaughterID(Int_t index){                            // get PDG ID of decay particle index
    return DecayDaughterID[index];
  }  
  Double_t Particle::GetDaughterWeight(Int_t index){                        // get weight of decay 
    return DecayDaughterWeight[index];                         // = 1 if random decay branch is generated with Decay()
  }                                                            // 
  TLorentzVector Particle::GetDecayDaughter(Int_t index){                // return 4 vector for decay particle index
    return DecayDaughter[index];
  }
  TVector3 Particle::GetDecayVtx(){
    return vtx;
  }
 
  Double_t Particle::Weight(){                                        // returns weight - not part of TLorentzVector
    return weight;
  }
  Int_t Particle::Charge(){                                           // returns charge - not part of TLorentzVector
    return charge;
  }
  TString Particle::Name(){                                           // returns name   - not part of TLorentzVector
    return name;
  }
  Int_t Particle::ID(){                                               // returns ID     - not part of TLorentzVector
    return id;
  }

// manipulate particle characteristics
  void Particle::ResetP(){                                               // resets momentum to 0 
    SetPxPyPzE(0.,0.,0.,mass);
  }
  void Particle::SetWeight(Double_t w){                                  // sets particle weight
    weight = w;
  };  
  void Particle::AddWeight(Double_t w){                                  // updates particle weight
    weight = w*weight;
  };

  void Particle::UpdateParticle(TLorentzVector &p){                      // update 4 vector of particle
    SetPtEtaPhiM(p.Pt(),p.Eta(),p.Phi(),p.M());
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate particle with flat pt from pt_low to pt_high
//                               phi over 2pi
//                               rapidity from -0.5 to 0.5
//
// input  Double_t pt_low      - lower pt bound 
//        Double_t pt_high     - upper pt bound
//
// sets 4 vector of particle
//
// Axel Drees    10/18/2018 
// updates       8/28/2020
//
void Particle::GenerateP(Double_t pt_low, Double_t pt_high, Bool_t rap) {
   Double_t pt,phi,rapidity,eta; 
   pt        = randy.Uniform(pt_low,pt_high);
   phi       = randy.Uniform(0.,2*pi);
   rapidity  = randy.Uniform(-.5,.5);

   eta       = rapidity;
//   generate flat in rapidity rather than pseudorapidity
   if (rap) eta = RapidityToEta(rapidity,pt,mass);
   SetPtEtaPhiM(pt,eta,phi,mass);
   GenerateVtx();
  }
 
////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate particle with flat pt from pt_low to pt_high
//                               phi over 2pi
//                               rapidity from TF1 
//
// input  Double_t pt_low      - lower pt bound 
//        Double_t pt_high     - upper pt bound
//
// sets 4 vector of particle
//
// Axel Drees    11/17/2022 
//
void Particle::GenerateP(Double_t pt_low, Double_t pt_high, TF1*RapiditySpectrum, Bool_t rap) {
   Double_t pt,phi,rapidity,eta; 
   pt        = randy.Uniform(pt_low,pt_high);
   phi       = randy.Uniform(0.,2*pi);
   rapidity  = RapiditySpectrum->GetRandom(); 

   eta       = rapidity;
//   generate flat in rapidity rather than pseudorapidity
   if (rap) eta = RapidityToEta(rapidity,pt,mass);
   SetPtEtaPhiM(pt,eta,phi,mass);
   GenerateVtx();
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate particle from input TF1 
//                               
// input  TF1 PtSpectrum       - TF1 with pt function (must be yield not invariant yield!)
//        TF1 PhiSpectrum      - TF1 with phi spectrum to generate from
//        TF1 RapiditySpectrum - TF1 with Rapidity spectrum to generate from          
// 
// sets 4 vector of particle
//
// Axel Drees    10/21/2018
// updated       8/28/2021
//
void Particle::GenerateP(TF1* PtSpectrum, TF1* PhiSpectrum, TF1* RapiditySpectrum, Bool_t rap) {
   Double_t pt,phi,rapidity,eta;

   pt        = PtSpectrum->GetRandom(); 
   phi       = PhiSpectrum->GetRandom();
   rapidity  = RapiditySpectrum->GetRandom(); 

   eta       = rapidity;
//   generate flat in rapidity rather than pseudorapidity
   if (rap) eta = RapidityToEta(rapidity,pt,mass);
   SetPtEtaPhiM(pt,eta,phi,mass);
   GenerateVtx();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate particle pt from input TF1, with uniform phi and eta 
//                               
// input  TF1 PtSpectrum       - TF1 with pt function (must be yield not invariant yield!)
// 
// sets 4 vector of particle
//
// Axel Drees    10/18/2018
// updated       8/28/2021
//
void Particle::GenerateP(TF1* PtSpectrum, Bool_t rap) {
   Double_t pt,phi,rapidity,eta;
   
   pt        = PtSpectrum->GetRandom(); 
   phi       = randy.Uniform(-pi,pi);
   rapidity  = randy.Uniform(-0.5,0.5);
   eta       = rapidity;
//   generate flat in rapidity rather than pseudorapidity
   if (rap) eta = RapidityToEta(rapidity,pt,mass);
   SetPtEtaPhiM(pt,eta,phi,mass);
   GenerateVtx();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate decay vertex for weak decays using ct
// 
// Axel Drees    11/20/2022 
//
  void Particle::GenerateVtx(){
    Double_t L;
     
    if (debug) std::cout << std::endl;
    if (debug) std::cout << " GenerateVTX ct " << ct << std::endl;
    if (debug) std::cout << " GenerateVTX theta " << Theta() << " phi " << Phi() << std::endl;
    if (debug) std::cout << " GenerateVTX energy " << Energy() << " mass " << mass << std::endl;

    if (ct != 0){
      Double_t ct1   = randy.Exp(ct);
      Double_t betagamma = Pt()/mass;
      L = ct1*betagamma;      
      if (debug) std::cout << " GenerateVTX ct " << ct1 << " beta gamma " << betagamma << std::endl;
      if (debug) std::cout << " GenerateVTX L " << L << std::endl;
      vtx.SetMagThetaPhi(L,Theta(),Phi());
    } else {
      vtx.SetXYZ(0.,0.,0.);
    }
 
  }


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define possible decay channel for !stable particle
//
// currently implemented decays
//
//    pi0 -> gamma gamma
//
// 
// Axel Drees 8/28/2021
//

void Particle::DefineDecays(){

  if (stable) return;                                    // this particle is stable no decays defined

// preset lorentz vector of possible decay particles
  TLorentzVector pi0(0.,0.,0.,pi0Mass);                  // define pi0 as possible decay particle
  TLorentzVector pip(0.,0.,0.,piMass);                  // define pi0 as possible decay particle
  TLorentzVector pim(0.,0.,0.,piMass);                  // define pi0 as possible decay particle
  TLorentzVector photon(0.,0.,0.,0.);                    // define photon as possible decay particle
  TLorentzVector electron(0.,0.,0.,eMass);               // define electron as possible decay particle
  TLorentzVector positron(0.,0.,0.,eMass);               // define positron as possible decay particle
  TLorentzVector mup(0.,0.,0.,muMass);                   // define mu+ as possible decay particle
  TLorentzVector mum(0.,0.,0.,eMass);                    // define mu- as possible decay particle
  TLorentzVector rho0(0.,0.,0.,rho0Mass);                // define rho0 as possible decay particle
  TLorentzVector omega(0.,0.,0.,omegaMass);              // define omega as possible decay particle
  TLorentzVector eta(0.,0.,0.,etaMass);                  // define eta as possible decay particle
  TLorentzVector Nucleon(0.,0.,0.,NucleonMass);
  TLorentzVector proton(0.,0.,0.,protonMass);            // define eta as possible decay particle
 
  if (debug) std::cout << "#### DefineDecay ####################################" << std::endl;

  NumberOfBranches = 0;

  if (name == "pi0") {                             
     if (debug) std::cout << " definging decay of pi0" << std::endl;
                                                         // pi0 -> gamma gamma
     DecayBranch[0].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[0].SetName("pi0->gg");
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_pi0_gg);
     daughter[0][0] = photon;                            // set daughter 1
     daughterID[0][0] = photonID;
     daughter[0][1] = photon;                            // set daughter 2 
     daughterID[0][1] = photonID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;
                                                         // pi0 -> e+ e- gamma
     DecayBranch[1].SetType("Dalitz");                   // this is Dalitz decay 
     DecayBranch[1].SetName("pi0->gee");                 // define name of decay
     DecayBranch[1].SetMassDistributions();              // set virtual photon mass distribution for "pi0-Dalitz"
     DecayBranch[1].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[1].SetBR(BR_pi0_Dalitz);
     daughter[1][0] = electron;                          // set daughter 1
     daughterID[1][0] = electronID;
     daughter[1][1] = positron;                          // set daughter 2 
     daughterID[1][1] = positronID;
     daughter[1][2] = photon;                            // set daughter 3 
     daughterID[1][2] = photonID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

  } else if (name == "eta") {
     if (debug) std::cout << " definging decay of eta" << std::endl;
                                                         // eta -> gamma gamma
     DecayBranch[0].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[0].SetName("eta->gg");
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_eta_gg);
     daughter[0][0] = photon;                            // set daughter 1
     daughterID[0][0] = photonID;
     daughter[0][1] = photon;                            // set daughter 2 
     daughterID[0][1] = photonID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;
                                                         // eta -> e+ e- gamma
     DecayBranch[1].SetType("Dalitz");                   // this is Dalitz decay 
     DecayBranch[1].SetName("eta->gee");                 // define name of decay
     DecayBranch[1].SetMassDistributions();              // set virtual photon mass distribution for "eta-Dalitz"
     DecayBranch[1].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[1].SetBR(BR_eta_Dalitz);
     daughter[1][0] = electron;                          // set daughter 1
     daughterID[1][0] = electronID;
     daughter[1][1] = positron;                          // set daughter 2 
     daughterID[1][1] = positronID;
     daughter[1][2] = photon;                            // set daughter 3 
     daughterID[1][2] = photonID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;
                                                         // eta -> mu+mu- gamma  
     DecayBranch[2].SetType("Dalitz");                   // this is Dalitz decay 
     DecayBranch[2].SetName("eta->gmm");                 // define name of decay
     DecayBranch[2].SetMassDistributions();              // set virtual photon mass distribution for "eta-Dalitz"
     DecayBranch[2].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[2].SetBR(BR_eta_Dalitz2);
     daughter[2][0] = mum;                               // set daughter 1
     daughterID[2][0] = mumID;
     daughter[2][1] = mup;                               // set daughter 2 
     daughterID[2][1] = mupID;
     daughter[2][2] = photon;                            // set daughter 3 
     daughterID[2][2] = photonID;
     DecayBranch[2].DefineDaughters(daughter[2]);           
     NumberOfBranches++;
                                                         // eta -> mu+mu- 
     DecayBranch[3].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[3].SetName("eta->mm");
     DecayBranch[3].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[3].SetBR(BR_eta_mm);
     daughter[3][0] = mum;                               // set daughter 1
     daughterID[3][0] = mumID;
     daughter[3][1] = mup;                               // set daughter 2 
     daughterID[3][1] = mupID;
     DecayBranch[3].DefineDaughters(daughter[3]);           
     NumberOfBranches++;
                                                         // eta -> pi0 pi0 pi0  
     DecayBranch[4].SetType("ThreeBody");                // this is a 3 body  decay 
     DecayBranch[4].SetName("eta->3pi0");                // define name of decay
     DecayBranch[4].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[4].SetBR(BR_eta_3pi0);
     daughter[4][0] = pi0;                               // set daughter 1
     daughterID[4][0] = pi0ID;
     daughter[4][1] = pi0;                               // set daughter 2 
     daughterID[4][1] = pi0ID;
     daughter[4][2] = pi0;                               // set daughter 3 
     daughterID[4][2] = pi0ID;
     DecayBranch[4].DefineDaughters(daughter[4]);           
     NumberOfBranches++;

  } else if (name == "omega") {
     if (debug) std::cout << " definging decay of omega" << std::endl;
                                                         // omega -> pi0 gamma
     DecayBranch[0].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[0].SetName("omega->pi0g");
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_omega_pi0g);
     daughter[0][0] = pi0;                               // set daughter 1
     daughterID[0][0] = pi0ID;
     daughter[0][1] = photon;                            // set daughter 2 
     daughterID[0][1] = photonID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;

     DecayBranch[1].SetType("Dalitz");                   // this is Dalitz decay 
     DecayBranch[1].SetName("omega->pi0ee");             // define name of decay
     DecayBranch[1].SetMassDistributions();              // set virtual photon mass distribution for "omega-Dalitz"
     DecayBranch[1].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[1].SetBR(BR_omega_pi0ee);
     daughter[1][0] = electron;                          // set daughter 1
     daughterID[1][0] = electronID;
     daughter[1][1] = positron;                          // set daughter 2 
     daughterID[1][1] = positronID;
     daughter[1][2] = pi0;                               // set daughter 3 
     daughterID[1][2] = pi0ID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

     DecayBranch[2].SetType("TwoBody");                   // this is two body decay  
     DecayBranch[2].SetName("omega->ee");                // define name of decay
     DecayBranch[2].SetMassDistributions();              // set virtual photon mass distribution for "omega-ee"
     DecayBranch[2].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[2].SetBR(BR_omega_ee);
     daughter[2][0] = electron;                          // set daughter 1
     daughterID[2][0] = electronID;
     daughter[2][1] = positron;                          // set daughter 2 
     daughterID[2][1] = positronID;
     DecayBranch[2].DefineDaughters(daughter[2]);           
     NumberOfBranches++;

     DecayBranch[3].SetType("Dalitz");                   // this is Dalitz decay to muons 
     DecayBranch[3].SetName("omega->pi0mm");             // define name of decay
     DecayBranch[3].SetMassDistributions();              // set virtual photon mass distribution for "omega-Dalitz"
     DecayBranch[3].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[3].SetBR(BR_omega_pi0mm);
     daughter[3][0] = mum;                               // set daughter 1
     daughterID[3][0] = mumID;
     daughter[3][1] = mup;                               // set daughter 2 
     daughterID[3][1] = mupID;
     daughter[3][2] = pi0;                               // set daughter 3 
     daughterID[3][2] = pi0ID;
     DecayBranch[3].DefineDaughters(daughter[3]);           
     NumberOfBranches++;

     DecayBranch[4].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[4].SetName("omega->mm");                // define name of decay
     DecayBranch[4].SetMassDistributions();              // set virtual photon mass distribution for "omega-ee"
     DecayBranch[4].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[4].SetBR(BR_omega_mm);
     daughter[4][0] = mum;                               // set daughter 1
     daughterID[4][0] = mumID;
     daughter[4][1] = mup;                               // set daughter 2 
     daughterID[4][1] = mupID;
     DecayBranch[4].DefineDaughters(daughter[4]);           
     NumberOfBranches++;

  } else if (name == "rho0") {
     if (debug) std::cout << " definging decay of rho0'" << std::endl;

     DecayBranch[0].SetType("TwoBody");                  // rho->ee 
     DecayBranch[0].SetName("rho0->ee");                 // define name of decay
     DecayBranch[0].SetMassDistributions();              // set virtual photon mass distribution for "rho-ee"
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_rho0_ee);
     daughter[0][0] = electron;                          // set daughter 1
     daughterID[0][0] = electronID;
     daughter[0][1] = positron;                          // set daughter 2 
     daughterID[0][1] = positronID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;

     DecayBranch[1].SetType("TwoBody");                  // rho->mu+mu- 
     DecayBranch[1].SetName("rho0->mm");                 // define name of decay
     DecayBranch[1].SetMassDistributions();              // set virtual photon mass distribution for "rho-ee"
     DecayBranch[1].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[1].SetBR(BR_rho0_mm);
     daughter[1][0] = mum;                               // set daughter 1
     daughterID[1][0] = mumID;
     daughter[1][1] = mup;                               // set daughter 2 
     daughterID[1][1] = mupID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

  } else if (name == "etap") {
     if (debug) std::cout << " definging decay of eta'" << std::endl;
                                                         // eta -> gamma gamma
     DecayBranch[0].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[0].SetName("etap->gg");
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_etap_gg);
     daughter[0][0] = photon;                            // set daughter 1
     daughterID[0][0] = photonID;
     daughter[0][1] = photon;                            // set daughter 2 
     daughterID[0][1] = photonID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;
  
     DecayBranch[1].SetType("Dalitz");                   // this is Dalitz decay 
     DecayBranch[1].SetName("etap->gee");                // define name of decay
     DecayBranch[1].SetMassDistributions();              // set virtual photon mass distribution for "eta'-Dalitz"
     DecayBranch[1].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[1].SetBR(BR_etap_Dalitz);
     daughter[1][0] = electron;                          // set daughter 1
     daughterID[1][0] = electronID;
     daughter[1][1] = positron;                          // set daughter 2 
     daughterID[1][1] = positronID;
     daughter[1][2] = photon;                            // set daughter 3 
     daughterID[1][2] = photonID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

     DecayBranch[2].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[2].SetName("etap->rho0g");
     DecayBranch[2].SetMassDistributions();              // set rho0 mass distribution for decay
     DecayBranch[2].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[2].SetBR(BR_etap_rho0g);
     daughter[2][0] = rho0;                              // set daughter 1
     daughterID[2][0] = rho0ID;
     daughter[2][1] = photon;                            // set daughter 2 
     daughterID[2][1] = photonID;
     DecayBranch[2].DefineDaughters(daughter[2]);           
     NumberOfBranches++;

     DecayBranch[3].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[3].SetName("etap->wg");
     DecayBranch[3].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[3].SetBR(BR_etap_wg);
     daughter[3][0] = omega;                             // set daughter 1
     daughterID[3][0] = omegaID;
     daughter[3][1] = photon;                            // set daughter 2 
     daughterID[3][1] = photonID;
     DecayBranch[3].DefineDaughters(daughter[3]);           
     NumberOfBranches++;

     DecayBranch[4].SetType("Dalitz");                   // this is Dalitz decay 
     DecayBranch[4].SetName("etap->gmm");                // define name of decay
     DecayBranch[4].SetMassDistributions();              // set virtual photon mass distribution for "eta'-Dalitz"
     DecayBranch[4].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[4].SetBR(BR_etap_Dalitz2);
     daughter[4][0] = mum;                               // set daughter 1
     daughterID[4][0] = mumID;
     daughter[4][1] = mup;                               // set daughter 2 
     daughterID[4][1] = mupID;
     daughter[4][2] = photon;                            // set daughter 3 
     daughterID[4][2] = photonID;
     DecayBranch[4].DefineDaughters(daughter[4]);           
     NumberOfBranches++;

     DecayBranch[5].SetType("ThreeBody");                // this is a 3 body decay 
     DecayBranch[5].SetName("etap->2pi0eta");            // define name of decay
     DecayBranch[5].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[5].SetBR(BR_etap_2pi0eta);
     daughter[5][0] = eta;                               // set daughter 1
     daughterID[5][0] = etaID;
     daughter[5][1] = pi0;                               // set daughter 2 
     daughterID[5][1] = pi0ID;
     daughter[5][2] = pi0;                               // set daughter 3 
     daughterID[5][2] = pi0ID;
     DecayBranch[5].DefineDaughters(daughter[5]);           
     NumberOfBranches++;

  } else if (name == "phi") {
     if (debug) std::cout << " definging decay of phi" << std::endl;

     DecayBranch[0].SetType("TwoBody");                  // this is 2 body decay 
     DecayBranch[0].SetName("phi->ee");                  // define name of decay
     DecayBranch[0].SetMassDistributions();              // set virtual photon mass distribution for "phi-ee"
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_phi_ee);
     daughter[0][0] = electron;                          // set daughter 1
     daughterID[0][0] = electronID;
     daughter[0][1] = positron;                          // set daughter 2 
     daughterID[0][1] = positronID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;

     DecayBranch[1].SetType("TwoBody");                  // this is a 2 body decay 
     DecayBranch[1].SetName("phi->mm");                  // define name of decay
     DecayBranch[1].SetMassDistributions();              // set virtual photon mass distribution for "phi->mm"
     DecayBranch[1].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[1].SetBR(BR_phi_mm);
     daughter[1][0] = mup;                               // set daughter 1
     daughterID[1][0] = mupID;
     daughter[1][1] = mum;                               // set daughter 2 
     daughterID[1][1] = mumID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

     DecayBranch[2].SetType("Dalitz");                  // this is a 2 body decay 
     DecayBranch[2].SetName("phi->etaee");                  // define name of decay
     DecayBranch[2].SetMassDistributions();              // set virtual photon mass distribution for "phi->mm"
     DecayBranch[2].SetNumberOfDecayParticles(3);        // number of daughters  
     DecayBranch[2].SetBR(BR_phi_etaee);
     daughter[2][0] = electron;                               // set daughter 1
     daughterID[2][0] = electronID;
     daughter[2][1] = positron;                               // set daughter 2 
     daughterID[2][1] = positronID;
     daughter[2][2] = eta;                               // set daughter 2 
     daughterID[2][2] = etaID;
     DecayBranch[2].DefineDaughters(daughter[2]);           
     NumberOfBranches++;

  } else if (name == "jpsi") {
     if (debug) std::cout << " definging decay of J/Psi" << std::endl;


     DecayBranch[0].SetType("TwoBody");                  // this is 2 body decay 
     DecayBranch[0].SetName("jpsi->ee");                  // define name of decay
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_jpsi_ee);
     daughter[0][0] = electron;                          // set daughter 1
     daughterID[0][0] = electronID;
     daughter[0][1] = positron;                          // set daughter 2 
     daughterID[0][1] = positronID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;

     DecayBranch[1].SetType("TwoBody");                  // this is a 2 body decay 
     DecayBranch[1].SetName("jpsi->mm");                  // define name of decay
     DecayBranch[1].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[1].SetBR(BR_jpsi_mm);
     daughter[1][0] = mup;                               // set daughter 1
     daughterID[1][0] = mupID;
     daughter[1][1] = mum;                               // set daughter 2 
     daughterID[1][1] = mumID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

  } else if (name == "psip") {
     if (debug) std::cout << " definging decay of psi(2s)" << std::endl;


     DecayBranch[0].SetType("TwoBody");                  // this is 2 body decay 
     DecayBranch[0].SetName("psip->ee");                  // define name of decay
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_psip_ee);
     daughter[0][0] = electron;                          // set daughter 1
     daughterID[0][0] = electronID;
     daughter[0][1] = positron;                          // set daughter 2 
     daughterID[0][1] = positronID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;

     DecayBranch[1].SetType("TwoBody");                  // this is a 2 body decay 
     DecayBranch[1].SetName("psip->mm");                  // define name of decay
     DecayBranch[1].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[1].SetBR(BR_psip_mm);
     daughter[1][0] = mup;                               // set daughter 1
     daughterID[1][0] = mupID;
     daughter[1][1] = mum;                               // set daughter 2 
     daughterID[1][1] = mumID;
     DecayBranch[1].DefineDaughters(daughter[1]);           
     NumberOfBranches++;

  } else if (name == "K0s") {
     if (debug) std::cout << " definging decay of K0s" << std::endl;

     DecayBranch[0].SetType("TwoBody");                  // this is two body decay 
     DecayBranch[0].SetName("K0s->2pi0");                // define name of decay
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_K0s_2pi0);
     daughter[0][0] = pi0;                               // set daughter 1
     daughterID[0][0] = pi0ID;     
     daughter[0][1] = pi0;                               // set daughter 2 
     daughterID[0][1] = pi0ID;
     DecayBranch[0].DefineDaughters(daughter[0]);           
     NumberOfBranches++;

  } else if (name == "Delta") {
     if (debug) std::cout << " definging decay of Delta" << std::endl;
 
     DecayBranch[0].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[0].SetName("Delta->Ng");
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_Delta_Ng);
     daughter[0][0] = Nucleon;                            // set daughter 1
     daughterID[0][0] = NucleonID;
     daughter[0][1] = photon;                            // set daughter 2 
     daughterID[0][1] = photonID;
     DecayBranch[0].DefineDaughters(daughter[2]);           
     NumberOfBranches++;

  } else if (name == "Lambda") {
     if (debug) std::cout << " definging decay of Lambda" << std::endl;
 
     DecayBranch[0].SetType("TwoBody");                  // this is a two body decay 
     DecayBranch[0].SetName("Lambda->ppi-");
     DecayBranch[0].SetNumberOfDecayParticles(2);        // number of daughters  
     DecayBranch[0].SetBR(BR_Lambda_ppim);
     daughter[0][0] = proton;                            // set daughter 1
     daughterID[0][0] = protonID;
     daughter[0][1] = pim;                            // set daughter 2 
     daughterID[0][1] = pimID;
     DecayBranch[0].DefineDaughters(daughter[2]);           
     NumberOfBranches++;

  } else {
     std::cout << name << " particle has no defined decay channels" << std::endl;
  }

  if (debug) std::cout << "####################################################" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Decay particle if not stable according to defined decay branches
//
// If multiple decay branches exist, choose one randomly according to Branching Ratios.
//
// If less than 100% of the branching ratio is defined, decay according to total defined branching ration 
// and update weigth of decay branch.  
//
// Axel Drees 8/29/2021 
//

void Particle::Decay(){
// define variables
  Int_t nbr=0;                                               
  Double_t br=0,brtot=0; 

  if (debug) std::cout << "---- Decay ------------------------------------------------" << std::endl;
  if (debug) std::cout << " number of decay branches defined: " << NumberOfBranches << std::endl;

// calculate total defined BR
  for (int i=0; i<NumberOfBranches; i++) {                   // calculate total defined branching ratio
     brtot = brtot + DecayBranch[i].GetBR();
  }
  if (debug) std::cout << " total defined branching ratio: " << brtot << std::endl;

// get random variale  
  Double_t r = randy.Uniform(0.,brtot);                      // get random variable 0 to total defined BR

// identify which decay branch to use
  for (int i=0; i<NumberOfBranches; i++) {                   // loop over decay branches to decide which decay to generate
    br = DecayBranch[i].GetBR();                             // get branching ratio of this branch
    if (r <= br) {                                           // use this branch if smaller than random variable  
      nbr = i;                                               // set branch to this index
      if (debug) std::cout << " branch CAR$#@W^*@&#*(T    " << i << " " << r << " " << br << " " << nbr << std::endl;
      r = 1.;                                                // set r = 1 to end loop
    } else {
      r = r-br;                                              // reduce r by tested BR
    } 
  }
  if (nbr>=0) ExecuteDecay(nbr,brtot);                       // execute selected decay
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Decay particle if not stable according to requested decay branches
//
// Axel Drees 8/31/2021 
//
void Particle::DecaySingleBranch(TString branch){

  Int_t nbr=-1;  
  Double_t br=0;                                             

  if (debug) std::cout << "---- Particle::DecaySingleBranch---------------------" << std::endl;
  if (debug) std::cout << " number of decay branches defined: " << NumberOfBranches << std::endl;
  if (debug) std::cout << " looking for decay branch: " << branch << std::endl;

// calculate total defined BR
  for (int i=0; i<NumberOfBranches; i++) {                   // find decay branch
    if (debug) std::cout << " decay branch " << DecayBranch[i].GetName() << " found: " << i << std::endl;
    if (branch == DecayBranch[i].GetName()) {
       nbr = i;
       br  = DecayBranch[i].GetBR();
       if (debug) std::cout << " desired decay branch  found: " << i << std::endl;
     }
   }
   if (nbr>=0) ExecuteDecay(nbr,br);                         // execute selected decay

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Decay particle if not stable with equal probability for all requested decays
//
// Axel Drees 11/13/2022 
//
void Particle::DecayMultiBranch(TString *ActiveBranch, Int_t nActive){
 
  Int_t nbr;
  Int_t activeBR[10];
  Double_t br=0 , brtot=0;                                             
  if (debug) std::cout << "---- DecayMultiBranch ----------------------------------------" << std::endl;
  if (debug) std::cout << " number of decay branches defined: " << NumberOfBranches << std::endl;
  if (debug) std::cout << " number of decay branches requested: " << nActive << std::endl;

// calculate total defined BR and identify decay branches to be generated
  int n=0;  
  for (int i=0; i<NumberOfBranches; i++) {       // loop over all defined branches 
    if (debug) std::cout << " decay branch " << DecayBranch[i].GetName() << " found: " << i << std::endl;
    for (int j=0; j<nActive; j++) {              // loop over active branches
      if (ActiveBranch[j] == DecayBranch[i].GetName()) {  // find requested branch
        brtot = brtot + DecayBranch[i].GetBR();  // accumulate total BR 
        activeBR[n] = i;                         // store branch #
        if (debug) std::cout << " active decay branch " << DecayBranch[i].GetName() << " found index: " << n << std::endl;
        n++;
      }
    }
  }
  if (debug) std::cout << " total defined branching ratio: " << brtot << std::endl;

// randomly select one decay branch from active branches

  n   = randy.Uniform(0.,nActive);               // get random variable 0 to total active BR
  nbr = activeBR[n];                             // get actual branch number 
  br = DecayBranch[nbr].GetBR();                 // get branching ratio of this branch
  br = br*NumberOfBranches;              
  ExecuteDecay(nbr,br);                          // execute selected decay
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Decay particle if not stable with equal probability for all decays
//
// Axel Drees 9/1/2021 
//
void Particle::DecayFlat(){
 
  Int_t nbr;
  Double_t br=0 , brtot=0;                                             
  if (debug) std::cout << "---- DecayFlat ----------------------------------------" << std::endl;
  if (debug) std::cout << " number of decay branches defined: " << NumberOfBranches << std::endl;

// calculate total defined BR
  for (int i=0; i<NumberOfBranches; i++) {            // calculate total defined branching ratio
     brtot = brtot + DecayBranch[i].GetBR();
  }
  if (debug) std::cout << " total defined branching ratio: " << brtot << std::endl;

// randomly select one decay branch  
  nbr = randy.Uniform(0.,NumberOfBranches);           // get random variable 0 to total defined BR
  br = DecayBranch[nbr].GetBR();                      // get branching ratio of this branch
  br = br*NumberOfBranches;              
  ExecuteDecay(nbr,br);                               // execute selected decay
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Decay particle through preset decay branch using MyDecay class
//
// input: Int_t nbr          index of selected decay branch
//        Double_t ww        weight assigned to decay
//        Int_t index=0      0 unless multiple decay branches are generated
// 
// private variables generated 
//
//         Int_t NumberOfDaughters   counter for number of daughter particles set
//         DecayDaughter[i]          4 vector of daughter particle i  
//         DecayDaughterID[i]        PDG ID of particle
//         DecayDaughterWeight[i]    weight assigned to decay
//         TString DecayName         name of the decay channel
//         Int_t DecayNumber         branch number of the decay channel in DefineDecays() [above]
//
// these variables are accessed through member functions 
//
// Axel Drees 9/1/2021
// Roli Esha  5/13/2022
//
void Particle::ExecuteDecay(Int_t nbr, Double_t ww, Int_t index){
  Int_t n=0;
  Double_t pt,phi,eta,m;

// set kinematics of parent particle
  pt = Pt();                                                  // get parent kinematics
  phi = Phi();                                                // pt, phi, and eta
  eta = Eta();
  m = M();                                                    // get mass
  DecayBranch[nbr].SetParent(pt,phi,eta,m);                   // set parent 4 vector for decay

  if (debug) std::cout << " generate decay partiles " << std::endl;
// generate decay
  DecayBranch[nbr].GenerateDecay();                           // generate decay particles 
  DecayName = DecayBranch[nbr].GetName();                     // get name of the decay branch
  DecayNumber = nbr;                                          // get branch number of the decay
  n = DecayBranch[nbr].GetNumberDecayParticles();             // get number of decay particle generated
  if (debug) std::cout << " number of decay partiles " << n << std::endl;

// update list of decay particles
  index = 0;                                                  // fill list  
  for (int i=0; i<n; i++){                                    // loop over generated decay particles
     DecayDaughter[index] = DecayBranch[nbr].GetDaughter(index);  // set 4 vector of daughter particle 
     if (debug) DecayDaughter[index].Print();                     
     DecayDaughterID[index] = daughterID[nbr][index];          // set PDG particle ID   
     DecayDaughterWeight[index] = ww;                          // set weight to total defined BR
     index++;                             
  }
  NumberOfDaughters = index;                                   // log number of decay particles in list
}



