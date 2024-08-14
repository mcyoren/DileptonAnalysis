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

  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V0", "inv_mass_ee_DCA_V0", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V0", "delt_phi_ee_DCA_V0", N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V1", "inv_mass_ee_DCA_V1", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V1", "delt_phi_ee_DCA_V1", N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V2", "inv_mass_ee_DCA_V2", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V2", "delt_phi_ee_DCA_V2", N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V3", "inv_mass_ee_DCA_V3", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V3", "delt_phi_ee_DCA_V3", N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V4", "inv_mass_ee_DCA_V4", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V4", "delt_phi_ee_DCA_V4", N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("inv_mass_ee_DCA_V5", "inv_mass_ee_DCA_V5", N_inv_m,    inv_m_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));
  Master3D.push_back(new TH3D("delt_phi_ee_DCA_V5", "delt_phi_ee_DCA_V5", N_phi_bins,   phi_bins, N_DCA,  DCA_bins,  N_pt, pt_bins));

  fcn_3Dx.clear();
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dx.push_back(get_phi_ee);

  fcn_3Dy.clear();
  fcn_3Dy.push_back(get_DCA_V0);
  fcn_3Dy.push_back(get_DCA_V0);
  fcn_3Dy.push_back(get_DCA_V1);
  fcn_3Dy.push_back(get_DCA_V1);
  fcn_3Dy.push_back(get_DCA_V2);
  fcn_3Dy.push_back(get_DCA_V2);
  fcn_3Dy.push_back(get_DCA_V3);
  fcn_3Dy.push_back(get_DCA_V3);
  fcn_3Dy.push_back(get_DCA_V4);
  fcn_3Dy.push_back(get_DCA_V4);
  fcn_3Dy.push_back(get_DCA_V5);
  fcn_3Dy.push_back(get_DCA_V5);

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

  const double pt = sqrt(px*px + py*py);
  
  if(pt>2.9) return 2.9;
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_recon_pt(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const double KEFF_pip = 1./0.98;
  const double KEFF_pim = 1./0.98;

  const double px = p1->get_px()*KEFF_pip + p2->get_px()*KEFF_pim;
  const double py = p1->get_py()*KEFF_pip + p2->get_py()*KEFF_pim;

  const double pt = sqrt(px*px + py*py);
  
  if(pt>2.9 && pt<5.0) return 2.9;
  return pt;
}

float Run14AuAuLeptonCombyHistos::get_mass_ee_recon(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const double me = 0.000510998918;

  const double KEFF_pip = 1.0/0.98;
  const double KEFF_pim = 1.0/0.98;

  const double px1 = p1->get_px()*KEFF_pip;
  const double py1 = p1->get_py()*KEFF_pip;
  const double pz1 = p1->get_pz()*KEFF_pip;

  const double px2 = p2->get_px()*KEFF_pim;
  const double py2 = p2->get_py()*KEFF_pim;
  const double pz2 = p2->get_pz()*KEFF_pim;

  const double pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
  const double pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
  const double es = sqrt(pm1 + me * me) + sqrt(pm2 + me * me);
  const double px = px1 + px2;
  const double py = py1 + py2;
  const double pz = pz1 + pz2;

  return sqrt(es * es - px * px - py * py - pz * pz);
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
  const double es = sqrt(pm1 + me * me) + sqrt(pm2 + me * me);
  const double px = px1 + px2;
  const double py = py1 + py2;
  const double pz = pz1 + pz2;

  return sqrt(es * es - px * px - py * py - pz * pz);
}

float Run14AuAuLeptonCombyHistos::get_DCA_V0(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  if(id11==id21 || id12==id22 || id13==id23) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return fabs(DCA);
}

float Run14AuAuLeptonCombyHistos::get_DCA_V1(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int match1 = p1->get_integer(Run14AuAuLeptonCombyEnum::MATCH);
  const int match2 = p2->get_integer(Run14AuAuLeptonCombyEnum::MATCH);

  if ( match1 == 0 || match2 == 0 ) return -999;

  const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  if(id11==id21 || id12==id22 || id13==id23) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return fabs(DCA);
}

float Run14AuAuLeptonCombyHistos::get_DCA_V2(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int match1 = p1->get_integer(Run14AuAuLeptonCombyEnum::MATCH);
  const int match2 = p2->get_integer(Run14AuAuLeptonCombyEnum::MATCH);

  if ( match1 < 11 || match2 < 11 ) return -999;

  const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  if(id11==id21 || id12==id22 || id13==id23) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return fabs(DCA);
}

float Run14AuAuLeptonCombyHistos::get_DCA_V3(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 < 111 || ghost2 < 111 ) return -999;

  const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  if(id11==id21 || id12==id22 || id13==id23) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return fabs(DCA);
}

float Run14AuAuLeptonCombyHistos::get_DCA_V4(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 < 110 || ghost2 < 110 ) return -999;

  const int match1 = p1->get_integer(Run14AuAuLeptonCombyEnum::MATCH);
  const int match2 = p2->get_integer(Run14AuAuLeptonCombyEnum::MATCH);

  if ( match1 == 0 || match2 == 0 ) return -999;

  const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  if(id11==id21 || id12==id22 || id13==id23) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return fabs(DCA);
}

float Run14AuAuLeptonCombyHistos::get_DCA_V5(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 < 111 || ghost2 < 111 ) return -999;

  const int match1 = p1->get_integer(Run14AuAuLeptonCombyEnum::MATCH);
  const int match2 = p2->get_integer(Run14AuAuLeptonCombyEnum::MATCH);

  if ( match1 == 0 || match2 == 0 ) return -999;

  const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
  const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
  const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

  if(id11==id21 || id12==id22 || id13==id23) return -999;

  const double DCA_X_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pip = p1->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA_X_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  const double DCA_Y_pim = p2->get_double(Run14AuAuLeptonCombyEnum::DCAY);

  const double DCA = sqrt( (DCA_X_pip-DCA_X_pim)*(DCA_X_pip-DCA_X_pim) + (DCA_Y_pip-DCA_Y_pim)*(DCA_Y_pip-DCA_Y_pim) );
  
  return fabs(DCA);
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
