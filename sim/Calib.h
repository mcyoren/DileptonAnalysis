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

const char inFile[][200] = {"/home/yoren/bnl/PHENIX/ee/event_display/input/my-2M_piminus_embed_v1.root", 
                            "/home/yoren/bnl/PHENIX/ee/event_display/input/analysis_output_ccbar_pthat0.root", 
                            "/home/yoren/bnl/PHENIX/ee/event_display/input/merged_output_photon.root",
                            "/home/yoren/bnl/PHENIX/ee/event_display/input/jpsi_sim.root",
                            "/home/yoren/bnl/PHENIX/ee/event_display/input/phi_ee_sim.root"};
const char field_file_path[200]="../AnaTrain/November/field_map.root";

// #include "AnaTrain/MyEvent.C"
bool ifcout = false;
//const float me2 = 0.000510998918*0.000510998918;

using namespace std;


