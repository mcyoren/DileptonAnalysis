#include "Run14AuAuLeptonCombyHistos.h"

using namespace std;

Run14AuAuLeptonCombyHistos::Run14AuAuLeptonCombyHistos()
{
  // 1D hists
  Master1D.clear();
  fcn_1Dx.clear();

  // 2D hists
  Master2D.clear();
  fcn_2Dx.clear();
  fcn_2Dy.clear();

  // 3D hists // DCA, sigma DCA, 2*with 3 hits, 2*phi_v and pc3
  Master3D.clear();

  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V0",  "inv_mass_ee_DCA_V0",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V1",  "inv_mass_ee_DCA_V1",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V2",  "inv_mass_ee_DCA_V2",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V3",  "inv_mass_ee_DCA_V3",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V4",  "inv_mass_ee_DCA_V4",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V5",  "inv_mass_ee_DCA_V5",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V6",  "inv_mass_ee_DCA_V6",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V7",  "inv_mass_ee_DCA_V7",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V8",  "inv_mass_ee_DCA_V8",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V9",  "inv_mass_ee_DCA_V9",  N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V10", "inv_mass_ee_DCA_V10", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V11", "inv_mass_ee_DCA_V11", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V12", "inv_mass_ee_DCA_V12", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V13", "inv_mass_ee_DCA_V13", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V14", "inv_mass_ee_DCA_V14", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V15", "inv_mass_ee_DCA_V15", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V16", "inv_mass_ee_DCA_V16", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V17", "inv_mass_ee_DCA_V17", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V18", "inv_mass_ee_DCA_V18", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V19", "inv_mass_ee_DCA_V19", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V20", "inv_mass_ee_DCA_V20", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V21", "inv_mass_ee_DCA_V21", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V22", "inv_mass_ee_DCA_V22", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V23", "inv_mass_ee_DCA_V23", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V24", "inv_mass_ee_DCA_V24", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V25", "inv_mass_ee_DCA_V25", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V26", "inv_mass_ee_DCA_V26", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V27", "inv_mass_ee_DCA_V27", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V28", "inv_mass_ee_DCA_V28", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V29", "inv_mass_ee_DCA_V29", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V30", "inv_mass_ee_DCA_V30", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V31", "inv_mass_ee_DCA_V31", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V32", "inv_mass_ee_DCA_V32", N_inv_m,    inv_m_bins, N_DCA,      DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V0",  "delt_phi_ee_DCA_V0",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V1",  "delt_phi_ee_DCA_V1",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V2",  "delt_phi_ee_DCA_V2",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V3",  "delt_phi_ee_DCA_V3",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V4",  "delt_phi_ee_DCA_V4",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V5",  "delt_phi_ee_DCA_V5",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V6",  "delt_phi_ee_DCA_V6",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V7",  "delt_phi_ee_DCA_V7",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V8",  "delt_phi_ee_DCA_V8",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V9",  "delt_phi_ee_DCA_V9",  N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V10", "delt_phi_ee_DCA_V10", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V11", "delt_phi_ee_DCA_V11", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V12", "delt_phi_ee_DCA_V12", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V13", "delt_phi_ee_DCA_V13", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V14", "delt_phi_ee_DCA_V14", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V15", "delt_phi_ee_DCA_V15", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V16", "delt_phi_ee_DCA_V16", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V17", "delt_phi_ee_DCA_V17", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V18", "delt_phi_ee_DCA_V18", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V19", "delt_phi_ee_DCA_V19", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V20", "delt_phi_ee_DCA_V20", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V21", "delt_phi_ee_DCA_V21", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V22", "delt_phi_ee_DCA_V22", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V23", "delt_phi_ee_DCA_V23", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V24", "delt_phi_ee_DCA_V24", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V25", "delt_phi_ee_DCA_V25", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V26", "delt_phi_ee_DCA_V26", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V27", "delt_phi_ee_DCA_V27", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V28", "delt_phi_ee_DCA_V28", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V29", "delt_phi_ee_DCA_V29", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V30", "delt_phi_ee_DCA_V30", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V31", "delt_phi_ee_DCA_V31", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V32", "delt_phi_ee_DCA_V32", N_phi_bins,   phi_bins, N_inv_m,  inv_m_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("mass_dca_pt_upsilon", "mass_dca_pt_upsilon",  100, 8, 12, 10, 0, 1000, 10, 0, 10));


  fcn_3Dx.clear();
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_mass_ee);

  fcn_3Dy.clear();
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dy.push_back(get_DCA);
  
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);

  fcn_3Dz.clear();
  fcn_3Dz.push_back(get_pt_V0);
  fcn_3Dz.push_back(get_pt_V1);
  fcn_3Dz.push_back(get_pt_V2);
  fcn_3Dz.push_back(get_pt_V3);
  fcn_3Dz.push_back(get_pt_V4);
  fcn_3Dz.push_back(get_pt_V5);
  fcn_3Dz.push_back(get_pt_V6);
  fcn_3Dz.push_back(get_pt_V7);
  fcn_3Dz.push_back(get_pt_V8);
  fcn_3Dz.push_back(get_pt_V9);
  fcn_3Dz.push_back(get_pt_V10);
  fcn_3Dz.push_back(get_pt_V11);
  fcn_3Dz.push_back(get_pt_V12);
  fcn_3Dz.push_back(get_pt_V13);
  fcn_3Dz.push_back(get_pt_V14);
  fcn_3Dz.push_back(get_pt_V15);
  fcn_3Dz.push_back(get_pt_V16);
  fcn_3Dz.push_back(get_pt_V17);
  fcn_3Dz.push_back(get_pt_V18);
  fcn_3Dz.push_back(get_pt_V19);
  fcn_3Dz.push_back(get_pt_V20);
  fcn_3Dz.push_back(get_pt_V21);
  fcn_3Dz.push_back(get_pt_V22);
  fcn_3Dz.push_back(get_pt_V23);
  fcn_3Dz.push_back(get_pt_V24);
  fcn_3Dz.push_back(get_pt_V25);
  fcn_3Dz.push_back(get_pt_V26);
  fcn_3Dz.push_back(get_pt_V27);
  fcn_3Dz.push_back(get_pt_V28);
  fcn_3Dz.push_back(get_pt_V29);
  fcn_3Dz.push_back(get_pt_V30);
  fcn_3Dz.push_back(get_pt_V31);
  fcn_3Dz.push_back(get_pt_V32);
  fcn_3Dz.push_back(get_pt_V0);
  fcn_3Dz.push_back(get_pt_V1);
  fcn_3Dz.push_back(get_pt_V2);
  fcn_3Dz.push_back(get_pt_V3);
  fcn_3Dz.push_back(get_pt_V4);
  fcn_3Dz.push_back(get_pt_V5);
  fcn_3Dz.push_back(get_pt_V6);
  fcn_3Dz.push_back(get_pt_V7);
  fcn_3Dz.push_back(get_pt_V8);
  fcn_3Dz.push_back(get_pt_V9);
  fcn_3Dz.push_back(get_pt_V10);
  fcn_3Dz.push_back(get_pt_V11);
  fcn_3Dz.push_back(get_pt_V12);
  fcn_3Dz.push_back(get_pt_V13);
  fcn_3Dz.push_back(get_pt_V14);
  fcn_3Dz.push_back(get_pt_V15);
  fcn_3Dz.push_back(get_pt_V16);
  fcn_3Dz.push_back(get_pt_V17);
  fcn_3Dz.push_back(get_pt_V18);
  fcn_3Dz.push_back(get_pt_V19);
  fcn_3Dz.push_back(get_pt_V20);
  fcn_3Dz.push_back(get_pt_V21);
  fcn_3Dz.push_back(get_pt_V22);
  fcn_3Dz.push_back(get_pt_V23);
  fcn_3Dz.push_back(get_pt_V24);
  fcn_3Dz.push_back(get_pt_V25);
  fcn_3Dz.push_back(get_pt_V26);
  fcn_3Dz.push_back(get_pt_V27);
  fcn_3Dz.push_back(get_pt_V28);
  fcn_3Dz.push_back(get_pt_V29);
  fcn_3Dz.push_back(get_pt_V30);
  fcn_3Dz.push_back(get_pt_V31);
  fcn_3Dz.push_back(get_pt_V32);
  fcn_3Dz.push_back(get_pt_V0);

  fcn_3D_weight.clear();
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
  fcn_3D_weight.push_back(set_weight_1);
}

// ee

float Run14AuAuLeptonCombyHistos::get_DCA(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  //const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  //const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);
  //const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  //const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);
  //const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  const double DCA1 = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA2 = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA = TMath::Sqrt( TMath::Abs ( DCA1 * DCA1 - DCA2 * DCA2 ) );
  if (DCA > 3000) return 3000;
  //if (DCA < -3000) return -3000;
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_mass_ee(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const double me = 0.000510998918;

  const double px1 = p1->get_px();
  const double py1 = p1->get_py();
  const double pz1 = p1->get_pz();

  const double px2 = p2->get_px();
  const double py2 = p2->get_py();
  const double pz2 = p2->get_pz();

  const double pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
  const double pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
  const double es = TMath::Sqrt(pm1 + me * me) + TMath::Sqrt(pm2 + me * me);
  const double px = px1 + px2;
  const double py = py1 + py2;
  const double pz = pz1 + pz2;

  return TMath::Sqrt(es * es - px * px - py * py - pz * pz);
}

float Run14AuAuLeptonCombyHistos::get_pt_V0(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  //if ( ghost1 > 0 || ghost2 > 0 ) std::cout<<"ev1 = "<<p1->get_double(Run14AuAuLeptonCombyEnum::ZVTX)<<" ev2 = "<<p2->get_double(Run14AuAuLeptonCombyEnum::ZVTX)<<std::endl;
  //if ( ghost1 > 0 || ghost2 > 0 ) std::cout<<"ghost1 = "<<p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT)<<" "<<ghost1<<" ghost2 = "<<p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT)<<" "<<ghost2<<std::endl;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 0 || hit_assoc2 < 0 ) return -999;
  
  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V1(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 0 || hit_assoc2 < 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V2(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 0 || hit_assoc2 < 0 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V3(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V4(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_pt_V5(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V6(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V7(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V8(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V9(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V10(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V11(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V12(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;
  
  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);

  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V13(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V14(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V15(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V16(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V17(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);  

  if( hadron_reject1 < 10000 || hadron_reject2 < 10000 ) return -999;
  if( hadron_reject1%10 <6 || hadron_reject2%10 <6 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V18(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);  

  if( hadron_reject1 < 10060 || hadron_reject2 < 10060 ) return -999;
  if( hadron_reject1%10 <6 || hadron_reject2%10 <6 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;
  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V19(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);  

  if( hadron_reject1 < 10050 || hadron_reject2 < 10050 ) return -999;
  if( hadron_reject1%10 <6 || hadron_reject2%10 <6 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V20(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);  

  if( hadron_reject1 < 10060 || hadron_reject2 < 10060 ) return -999;
  if( hadron_reject1%10 <6 || hadron_reject2%10 <6 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const float phi11 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi12 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi13 = p1->get_double(Run14AuAuLeptonCombyEnum::PHI3);
  const float phi21 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI1);
  const float phi22 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI2);
  const float phi23 = p2->get_double(Run14AuAuLeptonCombyEnum::PHI3);

  const float the11 = p1->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the12 = p1->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the13 = p1->get_double(Run14AuAuLeptonCombyEnum::THE3);
  const float the21 = p2->get_double(Run14AuAuLeptonCombyEnum::THE1);
  const float the22 = p2->get_double(Run14AuAuLeptonCombyEnum::THE2);
  const float the23 = p2->get_double(Run14AuAuLeptonCombyEnum::THE3);

  if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
      return -999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V21(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V22(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V23(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V24(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V25(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if (! ( (conv_reject1 < 0 && conv_reject2 > 999) || (conv_reject2 < 0 && conv_reject1 > 999) ) ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}


float Run14AuAuLeptonCombyHistos::get_pt_V26(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if (! ( (conv_reject1 == -10 && conv_reject2 > 999) || (conv_reject2 == -10 && conv_reject1 > 999) ) ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V27(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 > -10 || conv_reject2 > -10 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}


float Run14AuAuLeptonCombyHistos::get_pt_V28(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 > -1 || conv_reject2 > -1 ) return -999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V29(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if (! ( (conv_reject1 < 0 && conv_reject2 > 999) || (conv_reject2 < 0 && conv_reject1 > 999) ) ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V30(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if (! ( (conv_reject1 == -10 && conv_reject2 > 999) || (conv_reject2 == -10 && conv_reject1 > 999) ) ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V31(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 > -10 || conv_reject2 > -10 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_pt_V32(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 > -1 || conv_reject2 > -1 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_phi_ee(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const TVector3 ee1(p1->get_px(), p1->get_py(), 0); //p1->get_pz());
  const TVector3 ee2(p2->get_px(), p2->get_py(), 0); //p2->get_pz());

  return ee1.Angle(ee2);
}

float Run14AuAuLeptonCombyHistos::get_phi_V_ee(const UltraLightTrack *p1, const UltraLightTrack *p2)
{
  const TVector3 ee1(p1->get_px(), p1->get_py(), p1->get_pz());
  const TVector3 ee2(p2->get_px(), p2->get_py(), p2->get_pz());

  const TVector3 z(0, 0, 1);
  const TVector3 u = (ee1 + ee2).Unit();
  const TVector3 v = ee1.Cross(ee2);
  const TVector3 w = u.Cross(v);
  const TVector3 ua = (u.Cross(z)).Unit();

  return w.Angle(ua);
}


float Run14AuAuLeptonCombyHistos::set_weight_1(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  return 1;
}
