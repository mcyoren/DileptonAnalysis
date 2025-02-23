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
const double dNdy_pp[] = {42.2*95/257,0.93,0.788/1000.,0.12,0.01,4.3,10,0.133/1000};//pi0, phi, jpsi, ccbar, bbbar, omega, qgp, psi2s
const double T_pp[] = {0.113,0.139,0.149,0.3,0.3,0.110,0.3,0.164};
const double n_pp[] = {9.57, 10.8, 12.3, 0, 0, 9.78, 0, 14};
const double m_pp[] = {0.1349766, 1.019461, 3.096916, -1.019461, -9.46030, 0.78265, 0.938272, 3.68609};
const double Ncolls[] = {771,282.4,82.6,12.1,3};
const double Br_to_ee[] = {1.174e-2, 2.97e-4, 5.94e-2, 0.01, 0.01, 7.28e-5, 1e-2, 7.9e-3};
const double hith_pt_ratio[] = {1.0,0.4,1.0,0.5,0.5,0.9,0.5,0.5};
const double pi0_high_pt_ratio[] = {0.00623224,0.00228273,0.000667683,9.78082e-05,2.425e-05};
const TString part_names[] = {"pi0","phi","jpsi","ccbar","bbbar","omega","qgp","psi2s"};
// #include "AnaTrain/MyEvent.C"
bool ifcout = false;
//const float me2 = 0.000510998918*0.000510998918;

using namespace std;


