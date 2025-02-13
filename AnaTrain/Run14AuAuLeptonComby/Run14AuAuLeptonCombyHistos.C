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

  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V0",  "delt_phi_ee_DCA_V0",  N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V0",  "inv_mass_ee_DCA_V0",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V1",  "inv_mass_ee_DCA_V1",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V2",  "inv_mass_ee_DCA_V2",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V3",  "inv_mass_ee_DCA_V3",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V4",  "inv_mass_ee_DCA_V4",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V5",  "inv_mass_ee_DCA_V5",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V6",  "inv_mass_ee_DCA_V6",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V7",  "inv_mass_ee_DCA_V7",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V8",  "inv_mass_ee_DCA_V8",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V9",  "inv_mass_ee_DCA_V9",  N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V10", "inv_mass_ee_DCA_V10", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V11", "inv_mass_ee_DCA_V11", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V12", "inv_mass_ee_DCA_V12", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V13", "inv_mass_ee_DCA_V13", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V14", "inv_mass_ee_DCA_V14", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V15", "inv_mass_ee_DCA_V15", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V16", "inv_mass_ee_DCA_V16", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V17", "inv_mass_ee_DCA_V17", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V18", "inv_mass_ee_DCA_V18", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V19", "inv_mass_ee_DCA_V19", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V20", "inv_mass_ee_DCA_V20", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V21", "inv_mass_ee_DCA_V21", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V22", "inv_mass_ee_DCA_V22", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V23", "inv_mass_ee_DCA_V23", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));

  fcn_3Dx.clear();
  fcn_3Dx.push_back(get_phi_ee);
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

  fcn_3Dy.clear();
  fcn_3Dy.push_back(get_DCA_V11);
  fcn_3Dy.push_back(get_DCA_V0);
  fcn_3Dy.push_back(get_DCA_V1);
  fcn_3Dy.push_back(get_DCA_V2);
  fcn_3Dy.push_back(get_DCA_V3);
  fcn_3Dy.push_back(get_DCA_V4);
  fcn_3Dy.push_back(get_DCA_V5);
  fcn_3Dy.push_back(get_DCA_V6);
  fcn_3Dy.push_back(get_DCA_V7);
  fcn_3Dy.push_back(get_DCA_V8);
  fcn_3Dy.push_back(get_DCA_V9);
  fcn_3Dy.push_back(get_DCA_V10);
  fcn_3Dy.push_back(get_DCA_V11);
  fcn_3Dy.push_back(get_DCA_V12);
  fcn_3Dy.push_back(get_DCA_V13);
  fcn_3Dy.push_back(get_DCA_V14);
  fcn_3Dy.push_back(get_DCA_V15);
  fcn_3Dy.push_back(get_DCA_V16);
  fcn_3Dy.push_back(get_DCA_V17);
  fcn_3Dy.push_back(get_DCA_V18);
  fcn_3Dy.push_back(get_DCA_V19);
  fcn_3Dy.push_back(get_DCA_V20);
  fcn_3Dy.push_back(get_DCA_V21);
  fcn_3Dy.push_back(get_DCA_V22);
  fcn_3Dy.push_back(get_DCA_V23);

  fcn_3Dz.clear();
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);
  fcn_3Dz.push_back(get_pt);

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
}

// ee

float Run14AuAuLeptonCombyHistos::get_pt(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  if(pt>2.9) return 2.9;
  return pt;
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

float Run14AuAuLeptonCombyHistos::get_DCA_V0(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1%10 < 1 || hit_assoc2%10 < 1 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V1(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1%10 < 5 || hit_assoc2%10 < 5 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V2(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10 || hit_assoc2 < 10 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V3(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V4(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V5(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10 || conv_reject2 < 10 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V6(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10 || hit_assoc2 < 10 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V7(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V8(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V9(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V10(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10 || hit_assoc2 < 10 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V11(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V12(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V13(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V14(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10 || hit_assoc2 < 10 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V15(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V16(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V17(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 10000 || conv_reject2 < 10000 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V18(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;
  
  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

  if ( hadron_reject1/10%10 < 1 || hadron_reject2/10%10 < 1 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V19(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 100 || conv_reject2 < 100 ) return -999;
  
  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

  if ( hadron_reject1/10%10 < 5 || hadron_reject2/10%10 < 5 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V20(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

  if ( hadron_reject1/10%10 < 1 || hadron_reject2/10%10 < 1 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V21(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 1000 || hit_assoc2 < 1000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

  if ( hadron_reject1/10%10 < 5 || hadron_reject2/10%10 < 5 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V22(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

  if ( hadron_reject1/10%10 < 1 || hadron_reject2/10%10 < 1 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_DCA_V23(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -999;
  
  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

  if ( hadron_reject1/10%10 < 5 || hadron_reject2/10%10 < 5 ) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = TMath::Sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return DCA;
}

float Run14AuAuLeptonCombyHistos::get_phi_ee(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const TVector3 ee1(p1->get_px(), p1->get_py(), p1->get_pz());
  const TVector3 ee2(p2->get_px(), p2->get_py(), p2->get_pz());

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
