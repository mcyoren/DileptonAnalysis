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
#include <TH2D.h>
#include <TH3D.h>
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

      static float get_mass_ee(PHParticle *, const unsigned int,PHParticle *, const unsigned int);
      static float get_phi_ee(PHParticle *, const unsigned int,PHParticle *, const unsigned int);

      static float get_phi_V_ee(const UltraLightTrack *p1, const UltraLightTrack *p2);
      static float set_weight_1 (PHParticle *, const unsigned int,PHParticle *, const unsigned int);

};

const int N_pt = 16;  
const float pt_bins[N_pt+1] = {0.,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,4.0,6.0,9.0}; //4
const int N_inv_m = 35;  
const float inv_m_bins[N_inv_m+1] = {0.,0.05,0.1,0.15,0.25,0.35,0.45,0.55,0.65,0.725,0.8,0.875,0.95,1.0,1.05,1.15,1.3,1.5,1.75,2.0,2.25,2.5,2.75,2.9,3.0,3.05,3.1,3.15,3.2,3.3,3.45,3.55,3.65,3.75,3.85,4.5};
const int N_DCA = 42;  
//const float DCA_bins[N_DCA+1] = {0,2,5,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,300,350,400,450,500,600,700,800,900,1000,1200,1500,2000,3000,4000,5000,10000,100000};
const float DCA_bins[N_DCA+1] = {
  -100000, -5000, -3000, -1500, -1000, -800, -600, -450, -350, -260, -220,
  -180, -150, -130, -110, -90, -70, -50, -30, -15, -5,
   0,
   5, 15, 30, 50, 70, 90, 110, 130, 150, 180, 220, 260, 350, 450, 600, 800, 1000, 1500, 3000, 5000, 100000
};
const int N_phi_bins = 32;
const float phi_bins[N_phi_bins+1]= {0.0,0.126,0.251,0.377,0.503,0.628,0.754,0.88,1.005,1.131,1.194,1.257,1.319,1.382,1.445,1.508,1.571,1.634,1.696,1.759,1.822,1.885,1.948,2.011,2.136,2.262,2.388,2.513,2.639,2.765,2.89,3.016,3.142};

#endif

