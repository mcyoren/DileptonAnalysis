HH  HH EEEEEE LL     IIII   OOO    SSSSS
HH  HH EE     LL      II   O   O  SS
HHHHHH EEEEE  LL      II  O     O  SSSS
HH  HH EE     LL      II   O   O      SS
HH  HH EEEEEE LLLLLL IIII   OOO   SSSSS         

FAST SIMULATION          Axel Drees 2021
last updated 11/20/2022

Standalone fast simulation package that run in ROOT framework (version 6). 
Designed to aid systematic studies of measurements of photons and lepton 
pairs in high energy p+p, p+A and A+A collisions, in particular with PHENIX. 

All necessary files of the HELIOS package are in /HELIOSLibrary/. Example code that
uses HELIOS package can be found in the top directory. After updateing the download 
locations in /HELIOSLibrary/HELIOSLibrary.h and /HELIOSLibrary/GPRfiles/GPRFileNames.h 
the example code can be executed at the ROOT prompt in ROOT 6.  

content of /HELIOSLibrary/

main include file for user program
- HELIOSLibrary.h             use #include "...path../HELIOSLibrary.h" to include all necessary 
                              components of HELIOS in your code
Particle generator
- PDG.h                       Collection of particle properties and constants
- KinematicDistributions.h    Collection of kinematic distributions of various particles
- Particle.C, Particle.h      Particle Class - inherets from TLorentzVector, 
                              defines, particles and decays, includes 3-momentum generator 
                              and generates decay particles using decay class 
- Decay.C, Decay.h            Decay Class, used by Particle to generate decays
- DecayAuxilliaryFunctions.h  collection of functions used in Decay Class 

Interaction with matter
- InteractionWithMaterial.h   simulation of physics processes like Photon conversions, bremmsstrahlung etc. 

Experiment specific simulation
- PHENIXSetup.h               collection of constants defining PHENIXDetector
- PHENIXDetector.C, PHENIXDetector.h                  
                              PHENIX detector class used to approximate acceptance and detector response 
                              of the main PHENIX central arm detectors, EMCal and DC/PC1

Example codes running HELIOS fast simulation, can be executed at ROOT prompt (needs MyPlot class see below)

- TestDecay.C                 Example of use of Particle class without PHENIX detector simulation
                              tests Particle and Decay class: generates all implemented decays 
                              and plots photons and e+e- pairs properties

- TestEMCalReso.C             Combines use of Particle class and EMCal fast simulation to simulate 
                              pi0 and EMCal response. Calculate energy resolution and shift of scale from pion mass (gamma+gamma) and its RMS as function of energy.      

- TestElectron.C              Simulates electrons and determines E/p distribution.  

- TestDCReso.C                Simulates pi+ and pi- to determine momentum and position resolution 
                              of tracking as function of pt.

- DirectGamma.C               Example code to detemine photons from hadron decays
                              generates direct photons, and photons from hadron decays, 
                              true and reconstructed and histograms various properties like g_hadron/g_pion
 
- efFastSim.C                 simulates pi0 decay to gamma+gamma, converts one photon to e+e- and 
                              determines the conditional acceptance efficiency <ef> for the second 
                              decay photon to be detected in the EMCal if the e+e- pair was 
                              reconstructed. Also simulates systematic change of <ef> with variations of energy scale uncertainty, energy resolution, non linearity of energy scale, and different input pt distributions.                            

the example codes heavily use the MyPlot class for displaying graphs, this class is not needed 
to run HELIOS fast simulations. Code is located in /MyPlotting/

TCanvas *Canvas (name,  width, height, x-position, y-position, logy=0, logx=0 )
•   Creates default TCanvas of “name” with width and height and upper left corner at x-, x-position     
•   Log x or y scales optional
•   Sets  LeftMargin, RightMargin, TopMargin, and BottomMargin to default values, which can be changed (see below)

TFrame *Frame (name,  xAxisLable,  yAxisLable,  xmin,  xmax, ymin, ymax, centerTitle=0);
•   Creates a default TFrame with axis labels and axis ranges 
•   Center axis labels optional  
•   Axis title offset, size and font set to default values, which can be changed (see below) 

TLegend *Legend (name, xmin, ymin, xmax, ymax, ScaleDefaultSize=1.);
•   Creates a default legend at with 0 < xmin,xmax,ymin,ymax < 1
•   Sets default size and font, size can be scaled optimally

StyleMe (TObject, marker=20, color=1, msize=1.2,  lstyle=1 , lwidth=2);
•   Sets marker and line properties for TH1D, TF1, TGraph, TGraphErrors

Change default parameters before creating canvas:
•   SetLeftMargin(t)        default 0.18
•   SetRightMargin(t)       default 0.02
•   SetBottomMargin(t)      default 0.15
•   SetTopMargin(t)         default 0.02

Change default parameters before creating frame:
•   SetFont(i)              default   42 
•   SetxTitleOffset(t)      default    1.1  
•   SetyTitleOffset(t)      default   1.1
•   SetTitleSize(t)         default   0.06
•   void SetLabelSize(t)    default   0.05 
•   SetLegendSize(t)        default   0.05
•   SetLegendColor(t)       default  kBlack

Restore default parameters for next canvas/frame with Reset()

 

Description of HELIOS package from header files

////// HELIOSLibrary.h /////////////////////////////////////////////////////////////////////
//
//  includes all necessary code for the HELIOS packackage and 
//  defines string with full pathlength to local directory of HELIOS package
//  used for some data files
//  
//  this string need to be updated for installation on new computer 
//  e.g.  C:/root_v6.22.06/macros/HELIOS/
//
//  Axel Drees 11/15/2021
//
//
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
//  
//// KinematicDistributions.h //////////////////////////////////////////////////////////////////////////
//
// collection of momentum spectra of various particles
//
// currently implemented:
//
//   TF1 Hagedorn (name, mass, upperlim, lowerlim, A=504.5, a=0.5169, b=0.1626, p0=0.7366, n=-8.274)
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
// 11/11/2022   3 body dacys and updates    Axel Drees
//
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
/////// PHENIXSetup.h //////////////////////////////////////////////////////////////////////////////////////////
//
// Collection of constants defining PHENIX Setup 
// used in Class PHENIX 
//
// Axel Drees 9/14/2021
//  
//////   PHENIXDetector.h //////////////////////////////////////////////////////////////////
//
// HELIOS PHENIXDetector Class used to approximate acceptance and detector response 
// of the main PHENIX central arm detectors, EMCal and DC/PC1
// 
//
// PHENIXDetector()     constructor, sets private member variables not specific to a particle   
//
// high level member functions to generated PHENIX "default" response 
// - CharacterizeTrack(particle, charge, id)  sets private member variables for particle
// - ReconstructTrack(particle, id)           returns TLorentzVector with reconstructed charged track
// - ReconstructShower(particle, id)          returns TLorentzVector with reconstrcuted EMCal shower 
// 
// more specific public member functions that can be used for systematic studies by 
// variing default response, more detailed description can be found in PHENIXDetector.C
// 
// general:
// - InAcceptance(particle,charge)            Ideal PHENIX acceptance - checks if particle is in acceptance 
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
// - Phi_EMCal(){
// - Phi_DC()
// - Phi_DC()
// - Phi_RICH()
//
// Axel Drees
// updated 11/16/2021
//