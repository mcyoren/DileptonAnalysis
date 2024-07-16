//// Decay.C /////////////////////////////////////////////////////////////////////////////////////////////////
//
// see Decay.h for detailed description
//
// handels decay branches defined as particle properties in HELIOS 
// 
//
// Axel Drees 11/19/2019
// addaped to HELIOS 9/21/2021
//  
//
//

#include "Decay.h"

Decay::Decay(){             // Default constructor
  BR = 0;
  DecayType = "unkown";
  DecayName = "unkown";
  NumberDecayParticles = 0; 
}
Decay::~Decay(void){}      // destructor

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::GenerateDecay(){
//
// generate random decay particles depending on the defined decay type
//
  TLorentzVector parent;
  TLorentzVector daughter1;
  TLorentzVector daughter2;
  TLorentzVector daughter3;
  Double_t pt, eta, phi, pairMass; 

  if (debug)   std::cout << "Decay name " << DecayName << std::endl; 
  if (DecayType == "TwoBody") {
    if (DecayName == "omega->ee" or DecayName == "omega->mm" or DecayName == "rho0->ee" or DecayName == "rho0->mm"
        or DecayName == "phi->ee" or DecayName == "phi->mm" ) {
       TF1 *h = (TF1 *)gROOT->FindObject(DecayName);   // get leptonPair mass
       if (debug) std::cout << " root object found" << DecayName << std::endl;
       pairMass = h->GetRandom(); 
       pt  = Parent.Pt();
       eta = Parent.Eta();
       phi = Parent.Phi();
       if (debug)  std::cout << " ----- " << pairMass  << std:: endl; 
       parent.SetPtEtaPhiM(pt,eta,phi,pairMass);
    } else {
      parent = Parent;
    }
    daughter1.SetPtEtaPhiM(0.,0.,0.,Daughter[0].M());
    daughter2.SetPtEtaPhiM(0.,0.,0.,Daughter[1].M());
    TwoBodyDecay(parent,daughter1,daughter2);
    Daughter[0] = daughter1;
    Daughter[1] = daughter2;
  } else if (DecayType == "Dalitz" ){
    parent = Parent;
    daughter1.SetPtEtaPhiM(0.,0.,0.,Daughter[0].M());
    daughter2.SetPtEtaPhiM(0.,0.,0.,Daughter[1].M());
    daughter3.SetPtEtaPhiM(0.,0.,0.,Daughter[2].M());
    DalitzDecay(parent,daughter1,daughter2,daughter3);
    Daughter[0] = daughter1;
    Daughter[1] = daughter2;
    Daughter[2] = daughter3;
  } else if (DecayType == "ThreeBody" ){
    parent = Parent;
    daughter1.SetPtEtaPhiM(0.,0.,0.,Daughter[0].M());
    daughter2.SetPtEtaPhiM(0.,0.,0.,Daughter[1].M());
    daughter3.SetPtEtaPhiM(0.,0.,0.,Daughter[2].M());
    ThreeBodyDecay(parent,daughter1,daughter2,daughter3);
    Daughter[0] = daughter1;
    Daughter[1] = daughter2;
    Daughter[2] = daughter3;
  } else {
    BR = 0;
    NumberDecayParticles=0;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::TwoBodyDecay(TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2){
//
// Isotropic two body decay of parent into two particles
//
// input:  parent 4 vector
// output: daughter1 &2 4 vector 
//
// Axel Drees 10/17/2018 - checked to work for 2 photon decays
//            8/29/2021 - adapted from standalone function to this class
//
  TRandom3 randy = TRandom3(0);                // Random Generator
  Double_t pi = 3.14159;                              // define pi

  Double_t px,py,pz,E;                                // define generic 4 vector component 
  Double_t mass = parent.M();                         // mass of parent particle
  Double_t m1 = daughter1.M();                      // mass of decay particle 1 
  Double_t m2 = daughter2.M();                      // mass of decay particle 2
  TVector3 boost = parent.BoostVector();              // boost vector of parent particle

  // std::cout << "-------------------------------------------------- " << std::endl;
  // float beta = parent.P()/parent.E(); 
  // float betax = parent.Px()/parent.E(); 
  // float betay = parent.Py()/parent.E(); 
  // float betaz = parent.Pz()/parent.E(); 
  // float gamma = parent.E()/parent.M() ;
  // std::cout << " parent  beta=" << beta 
  //           << " gamma=" << gamma
  //           << " betax=" << betax
  //           << " betay=" << betay
  //           << " betaz=" << betaz << std::endl;
  // parent.Print();
  // boost.Print();
 
//
// decay in rest frame 
//
  Double_t p     = sqrt((mass*mass-(m1+m2)*(m1+m2))*(mass*mass-(m1-m2)*(m1-m2)))/2/mass; 
//  std::cout << "----- p --------- " << p << std::endl;
  Double_t phi = randy.Uniform(0.,2*pi);            // phi random from 0. to 2pi
  Double_t z   = randy.Uniform(-1.,1.);             // z used to calculate theta 
//  std::cout << "----- phi,z --------- " << phi << " " << z << std::endl;
  Double_t theta = acos(z);                           //     so that each angle in 4pi is equaly likely
  pz = p*cos(theta);                                  // longitudinal momentum
  px = p*sin(theta)*cos(phi);                         // transverse momentum x component
  py = p*sin(theta)*sin(phi);                         // transverse momentum y component
  E = sqrt(px*px+py*py+pz*pz+m1*m1);                  // energy of decay particle 1
  daughter1.SetPxPyPzE(px,py,pz,E);                 // set 4 vector of 1 decay particle
  E = sqrt(px*px+py*py+pz*pz+m2*m2);                  // energy of decay particle 1
  daughter2.SetPxPyPzE(-px,-py,-pz,E);              // set 4 vector of 2 decay particle in opposit direction
//  std::cout << "-----------------" << E << std::endl;
//
// boost to parent frame
//
  // std::cout << "first daughter" << std::endl;
  // daughter1.Print();
  // E = gamma*daughter1.E()+betax*gamma*daughter1.Px()+betay*gamma*daughter1.Py()+betaz*gamma*daughter1.Pz();
  // px = -gamma*betax*daughter1.E()-(1+(gamma-1)*betax*betax/beta/beta)*daughter1.Px()
  //      -((gamma-1)*betax*betay/beta/beta)*daughter1.Py()-((gamma-1)*betax*betaz/beta/beta)*daughter1.Pz();
  // py = -gamma*betay*daughter1.E()-((gamma-1)*betay*betax/beta/beta)*daughter1.Px()
  //      -(1+(gamma-1)*betay*betay/beta/beta)*daughter1.Py()-((gamma-1)*betay*betaz/beta/beta)*daughter1.Pz();
  // pz = -gamma*betaz*daughter1.E()-((gamma-1)*betaz*betax/beta/beta)*daughter1.Px()
  //      -((gamma-1)*betaz*betay/beta/beta)*daughter1.Py()-(1+(gamma-1)*betaz*betaz/beta/beta)*daughter1.Pz();
  // std::cout << px << "  " << py << "  " << pz << "  " << E << std::endl;
  daughter1.Boost(boost);                           // boost with parent momentum
  // daughter1.Print();
  // std::cout << "second daughter" << std::endl;
  // daughter2.Print();
  // E = gamma*daughter2.E()+betax*gamma*daughter2.Px()+betay*gamma*daughter2.Py()+betaz*gamma*daughter2.Pz();
  // px = -gamma*betax*daughter2.E()-(1+(gamma-1)*betax*betax/beta/beta)*daughter2.Px()
  //      -((gamma-1)*betax*betay/beta/beta)*daughter2.Py()-((gamma-1)*betax*betaz/beta/beta)*daughter2.Pz();
  // py = -gamma*betay*daughter2.E()-((gamma-1)*betay*betax/beta/beta)*daughter2.Px()
  //      -(1+(gamma-1)*betay*betay/beta/beta)*daughter2.Py()-((gamma-1)*betay*betaz/beta/beta)*daughter2.Pz();
  // pz = -gamma*betaz*daughter2.E()-((gamma-1)*betaz*betax/beta/beta)*daughter2.Px()
  //      -((gamma-1)*betaz*betay/beta/beta)*daughter2.Py()-(1+(gamma-1)*betaz*betaz/beta/beta)*daughter2.Pz();
  // std::cout << px << "  " << py << "  " << pz << "  " << E << std::endl;
  daughter2.Boost(boost);                           // boost with parent momentum
  // daughter2.Print();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::DalitzDecay(TLorentzVector &parent, TLorentzVector &decay1, TLorentzVector &decay2, TLorentzVector &decay3){
//
// Isotropic two body decay of parent into photon and virtual photon, followed by two body decay of virtual photon 
//
// input:  parent              4 vector
// output: decay particle 3    4 vector must have E = mass p=0
//         decay particle 1&2  4 vector of electron and positron, must have E = mass electron and p=0
//         MassDistribution    - function with mass distribution for Dalitz ee pair
// 
// output: parent and decay    4 vectors
//
// The decay is calculated in several seps
//
//   1. determine mass of ee pair using MassDistribution
//
// Axel Drees 11/16/2019 - 
//
  Double_t px,py,pz,E;                                // define generic 4 vector component 
  Double_t pairMass=0;
  Double_t mass = parent.M();                         // mass of parent particle
  Double_t m1 = decay1.M();                           // mass of decay particle 1 
  Double_t m2 = decay2.M();                           // mass of decay particle 2
  Double_t m3 = decay3.M();                           // mass of decay particle 2
  TLorentzVector leptonPair;
  TVector3 parentBoost = parent.BoostVector();        // boost vector of parent particle
//
// mass of ee pair
  TF1 *h = (TF1 *)gROOT->FindObject(DecayName);           // get leptonPair mass
  pairMass = h->GetRandom(); 
  leptonPair.SetPxPyPzE(0.,0.,0.,pairMass);                 // set leptonPair lorentz vector

  // std::cout << std::endl;
  if(debug) std::cout << "pair mass: " << pairMass << std::endl;      
//
// first 2 body decay parent->decay3 leptonPair
  TwoBodyDecay(parent,decay3,leptonPair);
//
// second 2 body decay of leptonPair to decay1 and decay2
  TwoBodyDecay(leptonPair,decay1,decay2);

  
  // parent.Print();
  // decay3.Print();                   
  // decay1.Print();
  // decay2.Print();
  // leptonPair.Print();                   
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::ThreeBodyDecay(TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2, TLorentzVector &daughter3){
//
// Isotropic three body decay of parent into three particles
//
// input:  parent 4 vector
// output: daughter1&2&3 4 vector 
//
// Axel Drees 11/11/2022 - based on EXODUS code, also check PDG chaper on kinematics
//
  TRandom3 randy = TRandom3(0);                       // Random Generator
  Double_t mp = parent.M();                           // parent mass
  Double_t md[3];                                     // 3 daughter masses
  md[0] = daughter1.M();
  md[1] = daughter2.M();
  md[2] = daughter3.M();

  TVector3 boost = parent.BoostVector();              // boost vector of parent particle

  Double_t m12,m13;                                   // enery of daugther 1&2 1&3 pair in parent rest frame
  Double_t m12lo = pow(md[0]+md[1], 2);               // min. energy for daughters 1&2
  Double_t m12hi = pow(mp-md[2], 2);                  // max. energy for daugthers 1&2
  Double_t m13lo = pow(md[0]+md[2], 2);               // dito daughters 1&3 
  Double_t m13hi = pow(mp-md[1], 2);
  
  Double_t t1, t2, t3, e1s, e3s, m13min, m13max;      // temp variables
  while (true) {                                      // populate Dalitz plot randomly in allowed region
                                                      // see PDG chapter on Kinematics
    m12 = sqrt(randy.Uniform(m12lo,m12hi));           // start with random value m12lo<m12<m12hi
    m13 = sqrt(randy.Uniform(m13lo,m13hi));           // dito m13
    e1s = (m12*m12+md[0]*md[0]-md[1]*md[1])/(2.0*m12);// used to calculate boundary of dalitz plot
    e3s = (mp*mp-m12*m12-md[2]*md[2])/(2.0*m12);  
    t1 = e1s*e1s-md[0]*md[0];                         // used to calculate min/max m13
    t2 = e3s*e3s-md[2]*md[2];
    t3 = pow(e1s+e3s,2);
    m13max = sqrt(t3 - pow(sqrt(t1)-sqrt(t2),2));     // maximum energy for 1&3 given m12
    m13min = sqrt(t3 - pow(sqrt(t1)+sqrt(t2),2));     // minimum energy for 1&3 given m12
    if( m13<=m13max && m13>=m13min ) break;           // found a valid m12 m13 combimation
  }

                                                    // generate daughter 3 in restframe of parent
  Double_t Ed3 = (mp*mp + md[2]*md[2] - m12*m12)/(2.0*mp);   // calculate energy of daughter 3
  Double_t pd3 = sqrt((Ed3+md[2])*(Ed3-md[2]));              // calculate momentum of daughter 3
  
  Double_t phi = randy.Uniform(0.,2*pi);              // phi random from 0. to 2pi
  Double_t z   = randy.Uniform(-1.,1.);               // z used to calculate theta 
  Double_t theta = acos(z);                           // so that each angle in 4pi is equaly likely
  Double_t pz = pd3*cos(theta);                       // longitudinal momentum
  Double_t px = pd3*sin(theta)*cos(phi);              // transverse momentum x component
  Double_t py = pd3*sin(theta)*sin(phi);              // transverse momentum y component
  daughter3.SetPxPyPzE(px,py,pz,Ed3);                 // set 4 vector of 3 decay particle
  daughter3.Boost(boost);                             // boost with parent momentum

  TLorentzVector temp = -(parent-daughter3);

//
// generate daughter 1 & 2 as two body decay of -(parent-daughter3)
  TwoBodyDecay(temp,daughter1,daughter2);

  if(debug) {
    std::cout << std::endl;
    std::cout << " 3 body: " << " daughter3 "; daughter3.Print(); 
    std::cout << " 3 body: " << " m12 " << m12 << " " << temp.M() << std::endl;
    std::cout << " 3 body: " << " daughter1 "; daughter1.Print(); 
    std::cout << " 3 body: " << " daughter1 "; daughter1.Print(); 
    std::cout << " 3 body: " << " daughter2 "; daughter2.Print(); 
    TLorentzVector mom = daughter1+daughter2+daughter3;
    std::cout << " 3 body: " << " parent mass " << parent.M() << "    " << mom.M() << std::endl;
    std::cout << " 3 body: " << " parent      "; parent.Print(); 
    std::cout << " 3 body: " << " mom         "; mom.Print(); 
    std::cout << std::endl;
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::SetMassDistributions(){
//
// define TF1's used to sample mass for various decays including
//
//     - virtual photon mass from Dalitz Decays
//     - short lived resonances omega, phi J/psi etc
//     - longer lived rho meson
//
//
//
//   
  Double_t c1 = 2*alpha/3/pi;                         // constant factor
  
  if (DecayName == "pi0->gee") {
      TF1 *h_Mass = new TF1(DecayName,f_pi0Dalitz,2.*eMass,pi0Mass);

  } else if (DecayName == "eta->gee") {
      TF1 *h_Mass = new TF1(DecayName,f_etaDalitz,2.*eMass,etaMass);
  } else if (DecayName == "eta->gmm") {
      TF1 *h_Mass = new TF1(DecayName,f_etaDalitz2,2.*muMass,etaMass);

  } else if (DecayName == "etap->gee") {
      TF1 *h_Mass = new TF1(DecayName,f_etapDalitz,2.*eMass,etapMass);
  } else if (DecayName == "etap->gmm") {
      TF1 *h_Mass = new TF1(DecayName,f_etapDalitz2,2.*muMass,etapMass);

  } else if (DecayName == "etap->rho0g") {
      TF1 *h_Mass = new TF1(DecayName,f_etapRho0g,2.*eMass,etapMass);

  } else if (DecayName == "omega->pi0ee") {
      TF1 *h_Mass = new TF1(DecayName,f_omegaDalitz,2.*eMass,omegaMass-pi0Mass);

  } else if (DecayName == "omega->ee") {
      TF1 *h_Mass = new TF1(DecayName,f_omegaee,2.*pi0Mass,2*omegaMass); 

  } else if (DecayName == "omega->pi0mm") {
      TF1 *h_Mass = new TF1(DecayName,f_omegaDalitz2,2.*muMass,omegaMass-pi0Mass);
  } else if (DecayName == "omega->mm") {
      TF1 *h_Mass = new TF1(DecayName,f_omegamm,2.*pi0Mass,2*omegaMass); 


  } else if (DecayName == "rho0->ee") {
      TF1 *h_Mass = new TF1(DecayName,f_rho0ee,2.*pi0Mass,2*rho0Mass); 
  } else if (DecayName == "rho0->mm") {
      TF1 *h_Mass = new TF1(DecayName,f_rho0ee,2.*pi0Mass,2*rho0Mass); 

  } else if (DecayName == "phi->ee") {
      TF1 *h_Mass = new TF1(DecayName,f_phiee,2.*pi0Mass,2*phiMass); 
  } else if (DecayName == "phi->mm") {
      TF1 *h_Mass = new TF1(DecayName,f_phimm,2.*pi0Mass,2*phiMass); 
  }
  else if (DecayName == "phi->etaee") {
      TF1 *h_Mass = new TF1(DecayName,f_phietaee,2.*eMass,phiMass-etaMass); 
  }
}

