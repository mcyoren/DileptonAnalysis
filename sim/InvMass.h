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
const double dNdy_pp[13] = {42.2*95/257,42.2*11/257,8.6*42/257,4.3,4.3,0.92,0.93,0.93,0.788/1000.,0.133/1000,0.12,0.01,10};//pi0, eta, rho->ee, omega->pi0ee, omega, eta'->gamma ee, phi, jpsi, psi2s, ccbar, bbbar, qgp
const double T_pp[13] = {0.113,0.113,0.139,0.139,0.110,0.110,0.139,0.139,0.149,0.164,0.3,0.3,0.3};
const double n_pp[13] = {9.57, 9.57, 10.8, 9.7, 9.78, 9.78, 10.8, 10.8, 12.3, 14, 0, 0, 0};
const double m_pp[13] = {0.1349766, 0.547862, 0.77526, 0.78265, 0.78265, 0.95778, 1.019461, 1.019461, 3.096916, 3.68609, -1.019461, -9.46030, 1};
const double Ncolls[5] = {771,282.4,82.6,12.1,3};
const double Br_to_ee[13] = {1.174e-2, 6.9e-3, 4.73e-5, 7.38e-5, 7.7e-4, 4.73e-4, 2.98e-4, 1.08e-4, 5.94e-2, 7.9e-3, 0.01, 0.01, 1e-2};
///////////////////////////     0        1        2        3        4      5        6        7       8        9    10   11
const double hith_pt_ratio[13] = {1.0,0.48,    1.0,     0.9,    0.9,   0.25,     0.4,     0.4,   0.01,      0.5,  0.5,0.5,0.5};
const double pi0_high_pt_ratio[5] = {0.00623224,0.00228273,0.000667683,9.78082e-05,2.425e-05};
const TString part_names[13] = {"pi0","eta","rho->ee","omega","omega->pi0ee","eta'->gamma ee","phi->ee","phi->etaee","jpsi","psi2s","ccbar","bbbar","qgp"};
/////////////////////////////      0     1      2        3             4               5              6      7           8     9       10      11     12
// #include "AnaTrain/MyEvent.C"
bool ifcout = false;
//const float me2 = 0.000510998918*0.000510998918;788