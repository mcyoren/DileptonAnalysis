#include <iostream>
#include <fstream>
#include <stdio.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TGraph.h"
#include "TCanvas.h"
//gSystem->Load("lib/MyEvent_C.so");
//#include "../../ee_QA/AnaTrain/November/MyEvent.h"
#include "../AnaTrain/Run14AuAuLeptonComby/MyEvent.h"
#include "MyEvent_OLD.h"

const int in_id = 2; //pi+,pi-,e+,e- = 8,9,3,2
const TString inFile[] = {"/home/yoren/bnl/PHENIX/ee/AuAu/my-10M_photon_embed_single_v0.root", //2M_ccbar_embed_pythia8_v0  2M_phi_embed_helios_v0 
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_ccbar_embed_pythia8_v0.root", 
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_jpsi_embed_helios_v0.root",
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_phi_embed_helios_v0.root",
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_photon_embed_helios_v0.root"};
const char field_file_path[200]="../AnaTrain/November/field_map.root";
//const double dNdy_pp[] = {43.5,0.41,0.759/1000.,0.12,0.01,4.3,10,0.133/1000};//pi0, phi, jpsi, ccbar, bbbar, omega, qgp, psi2s
const double dNdy_pp[13] = {42.2*95/257,42.2*11/257,8.6*42/257,4.3,4.3,0.92,0.93,0.93,1.77e-5,0.133/1000/42.2,0.12,0.01,10};//pi0, eta, rho->ee, omega->pi0ee, omega, eta'->gamma ee, phi, jpsi, psi2s, ccbar, bbbar, qgp
const double T_pp[13] = {0.113,0.113,0.139,0.139,0.110,0.110,0.139,0.139,0.149,0.164,0.3,0.3,0.3};
const double n_pp[13] = {9.57, 9.57, 10.8, 9.7, 9.78, 9.78, 10.8, 10.8, 12.3, 14, 0, 0, 0};
const double m_pp[13] = {0.1349766, 0.547862, 0.77526, 0.78265, 0.78265, 0.95778, 1.019461, 1.019461, 3.096916, 3.68609, -1.019461, -9.46030, 1};
const double Ncolls[5] = {771,282.4,82.6,12.1,3};
const double Br_to_ee[13] = {1.174e-2, 6.9e-3, 4.73e-5, 7.38e-5, 7.7e-4, 4.73e-4, 2.98e-4, 1.08e-4, 5.94e-2, 7.9e-3, 1,1,1};
///////////////////////////     0        1        2        3        4      5        6        7       8        9    10   11
const double hith_pt_ratio[13] = {1.0,0.48,    1.0,     0.78,    0.78,   0.25,     0.4,     0.4,   0.01,      0.5,  1,1,1};
const double pi0_high_pt_ratio[5] = {0.00623224,0.00228273,0.000667683,9.78082e-05,2.425e-05};
const TString part_names[13] = {"pi0","eta","rho->ee","omega","omega->pi0ee","eta'->gamma ee","phi->ee","phi->etaee","jpsi","psi2s","ccbar","bbbar","qgp"};
/////////////////////////////      0     1      2        3             4               5              6      7           8     9       10      11     12
const double RAA_jpsi[5] = {0.2,0.3,0.7,0.8,0.8};
// #include "AnaTrain/MyEvent.C"
bool ifcout = false;
//const float me2 = 0.000510998918*0.000510998918;788
const double RAA_pi0[5] = {0.25,0.3,0.5,0.7,0.8};
const double pi0_multiplicity_AuAu[6] = {
  240.75, // 0–20%
  117.15, // 20–40%
   47.10, // 40–60%
   13.35, // 60–80%
    5.00,  // 80–93%
    95.7, //0-93%

};

const double hag_params_try[9][5] = {
  //         c      ,     a     ,     b     ,    p0     ,     n
  {1505.2,  0.4347,  0.2610,  0.7205,  8.319}, // 0–20%   
  {735.9 ,  0.4103,  0.1884,  0.7316,  8.272}, // 20–40%
  {296.2 ,  0.4038,  0.1290,  0.7313,  8.207}, // 40–60%
  {118.1 ,  0.4416 , 0.0559 , 0.7230 , 8.163}, // 60–80%
  {51.1  ,  0.2470,  0.0619,  0.7101,  8.459}  // 80–93%
};
const double hag_params_try1[9][5] = {
  //         c      ,     a     ,     b     ,    p0     ,     n
  {1166.0-60, 0.5260 , 0.1628 , 0.7511 , 8.348}, // 0–20%
  {643.0    , 0.4534 , 0.1325 , 0.7525 , 8.333}, // 20–40%
  {297.85   , 0.4220 , 0.1027 , 0.7258 , 8.220}, // 40–60%
  {118.1    , 0.4416 , 0.0559 , 0.7230 , 8.163}, // 60–80%
  {51.1     , 0.2470 , 0.0619 , 0.7101 , 8.459}  // 80–93%
};
const double hag_params[5][5] = {
  //        c       ,     a     ,     b     ,    p0     ,     n
  {1166.0   , 0.5457 , 0.17865 , 0.7470 , 8.311 }, // 0–20% (avg of 0–10 and 10–20)
  {643.0    , 0.4717 , 0.14155 , 0.7502 , 8.316 }, // 20–40%
  {297.85   , 0.42765, 0.1124  , 0.73215, 8.241 }, // 40–60%
  {93.65    , 0.3633 , 0.0453  , 0.7501 , 8.348 }, // 60–80%
  {51.1     , 0.2470 , 0.0619  , 0.7101 , 8.459 }  // 80–93% (same as 80–92%)
};
const double hag_params_full[9][5] = {
  //         c      ,     a     ,     b     ,    p0     ,     n
  {1331.0   , 0.5654 , 0.1945 , 0.7429 , 8.274}, // 0–10%
  {1001.0   , 0.5260 , 0.1628 , 0.7511 , 8.348}, // 10–20%
  {750.7    , 0.4900 , 0.1506 , 0.7478 , 8.299}, // 20–30%
  {535.3    , 0.4534 , 0.1325 , 0.7525 , 8.333}, // 30–40%
  {364.5    , 0.4333 , 0.1221 , 0.7385 , 8.261}, // 40–50%
  {231.2    , 0.4220 , 0.1027 , 0.7258 , 8.220}, // 50–60%
  {118.1    , 0.4416 , 0.0559 , 0.7230 , 8.163}, // 60–70%
  {69.2     , 0.2850 , 0.0347 , 0.7787 , 8.532}, // 70–80%
  {51.1     , 0.2470 , 0.0619 , 0.7101 , 8.459}  // 80–92%
};


// Scale factors for sectors 0–7
const double scale[8] = {
    0.988, 0.990, 0.985, 0.981,
    0.980, 0.985, 1.025, 1.023
};

// Smear c1 values for sectors 0–7
const double smear_c1[8] = {
    0.055, 0.066, 0.055, 0.059,
    0.057, 0.056, 0.072, 0.072
};

// Smear c2 values for sectors 0–7
const double smear_c2[8] = {
    0.011, 0.017, 0.011, 0.013,
    0.012, 0.012, 0.063, 0.062
};
