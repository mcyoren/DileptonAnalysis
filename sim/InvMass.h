#include <iostream>
#include <fstream>
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
const char inFile[][200] = {"/home/yoren/bnl/PHENIX/ee/AuAu/my-10M_photon_embed_single_v0.root", //2M_ccbar_embed_pythia8_v0  2M_phi_embed_helios_v0 
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_ccbar_embed_pythia8_v0.root", 
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_jpsi_embed_helios_v0.root",
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_phi_embed_helios_v0.root",
                            "/home/yoren/bnl/PHENIX/ee/AuAu/2M_photon_embed_helios_v0.root"};
const char field_file_path[200]="../AnaTrain/November/field_map.root";

// #include "AnaTrain/MyEvent.C"
bool ifcout = false;
//const float me2 = 0.000510998918*0.000510998918;

using namespace std;


