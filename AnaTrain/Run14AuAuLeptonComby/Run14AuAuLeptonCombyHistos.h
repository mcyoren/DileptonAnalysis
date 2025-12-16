#ifndef __RUN14AUAULEPTONCOMBYHISTOS_H__
#define __RUN14AUAULEPTONCOMBYHISTOS_H__

#include "Run14AuAuLeptonCombyEnum.h"
#include "CabanaBoy/cbMasterHistos.h"
#include "UltraLight/UltraLight.h"
#include "UltraLight/UltraLightTrack.h"

#include <string>
#include <deque>
#include <vector>
#include <iostream>
#include <cmath>
#include <TH2.h>
#include <TH3.h>
#include "TMath.h"
#include "TVector3.h"

class Run14AuAuLeptonCombyHistos: public cbMasterHistos {

public:

    Run14AuAuLeptonCombyHistos();

      static float get_pt_V0 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V1 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V2 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V3 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V4 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V5 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V6 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V7 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V8 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V9 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V10(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V11(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V12(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V13(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V14(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V15(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V16(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V17(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V18(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V19(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V20(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V21(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V22(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V23(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V24(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V25(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V26(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V27(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V28(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V29(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V30(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V31(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_pt_V32(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      

      static float get_DCA(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_DCA1(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_DCA2(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_DCA_BG(PHParticle *, const unsigned int,PHParticle *, const unsigned int);

      static float get_mass_ee(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_phi_ee(PHParticle *, const unsigned int,PHParticle *, const unsigned int);

      static float get_phi_V_ee(const UltraLightTrack *p1, const UltraLightTrack *p2);
      static float set_weight_1 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float set_weight_2 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);

};

const int N_pt = 16;  
const float pt_bins[N_pt+1] = {0.,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,4.0,6.0,9.0}; //4
const int N_inv_m = 35;  
const float inv_m_bins[N_inv_m+1] = {0.,0.05,0.1,0.15,0.25,0.35,0.45,0.55,0.65,0.725,0.8,0.875,0.95,1.0,1.05,1.15,1.3,1.5,1.75,2.0,2.25,2.5,2.75,2.9,3.0,3.05,3.1,3.15,3.2,3.3,3.45,3.55,3.65,3.75,3.85,4.5};
//const int N_DCA = 42;  
//const float DCA_bins[N_DCA+1] = {0,2,5,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,300,350,400,450,500,600,700,800,900,1000,1200,1500,2000,3000,4000,5000,10000,100000};
//const float DCA_bins[N_DCA+1] = {
//  -100000, -5000, -3000, -1500, -1000, -800, -600, -450, -350, -260, -220,
//  -180, -150, -130, -110, -90, -70, -50, -30, -15, -5,
//   0,
//   5, 15, 30, 50, 70, 90, 110, 130, 150, 180, 220, 260, 350, 450, 600, 800, 1000, 1500, 3000, 5000, 100000
//};
const int N_DCA = 40;
const float DCA_bins[N_DCA + 1] = {
  -3100, -2700, -2300, -1900,
  -1500, -1300, -1100, -900,                              
  -700, -600, -500, -400,                               
  -300, -250, -200, -150,                               
  -100, -75, -50, -25, 0, 25, 50, 75, 100,
  150, 200, 250, 300,
  400, 500, 600, 700,                                   
  900, 1100, 1300, 1500,
  1900, 2300, 2700, 3100      
};

const int N_phi_bins = 100;
const double phi_bins[N_phi_bins+1]= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
const int N_dca_bins = 17;
const double new_DCA_bins[N_dca_bins+1] = {0.,25,50,75,100,125,150,175,200,225,250,275,300,400,500,750,1000,2000};
const int N_pt_bins = 10;
const double new_pt_bins[N_pt_bins+1] = {0.0,0.3,0.6,0.9,1.2,1.6,2.0,2.5,3.0,5.0,10.0};
// --- your edges (Double_t) ---
const float mass_edges[] = {
  0.,0.025,0.05,0.075,0.1,0.15,0.25,0.35,0.45,0.55,0.65,0.7,0.75,0.8,0.9,0.95,1.0,1.05,
  1.15,1.3,1.5,1.75,2.0,2.25,2.5,2.75,2.9,3.0,3.05,3.1,3.15,3.2,3.3,3.45,
  3.55,3.65,3.75,3.85,4.5
};
const int nMassBins = sizeof(mass_edges)/sizeof(mass_edges[0]) - 1;

// --- example: 3D THnF: mass, DCA, pT ---
const int myndim = 4;
//Int_t    nbins[ndim] = { nMassBins, 40, 12, 6 };
//Double_t xmin [ndim] = { mass_edges[0], 0.0, 0.0, 0.0 };
//Double_t xmax [ndim] = { mass_edges[nMassBins], 1000.0, 4.8, 3.1415926535 };
const int    mynbins[myndim] = { nMassBins, N_DCA, N_pt_bins, N_phi_bins };
const double myxmin [myndim] = { mass_edges[0], new_DCA_bins[0], new_pt_bins[0], phi_bins[0] };
const double myxmax [myndim] = { mass_edges[nMassBins], new_DCA_bins[N_dca_bins], new_pt_bins[N_pt_bins], phi_bins[N_phi_bins] };

#endif

