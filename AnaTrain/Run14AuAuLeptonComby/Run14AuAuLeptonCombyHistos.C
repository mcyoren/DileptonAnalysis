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
  fcn_3Dx.clear();
  fcn_3Dy.clear();
  fcn_3Dz.clear();
  fcn_3D_weight.clear();

  for (int i = 0; i < 24; ++i) {
    Master3D.push_back(new TH3D(Form("inv_mass_ee_DCA_V%d", i), Form("inv_mass_ee_DCA_V%d;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]", i), 180, 0, 4.5, 40, 0, 1000, 12, 0, 4.8));
    fcn_3Dx.push_back(get_mass_ee);
    if (i<20 && i!=3) fcn_3Dy.push_back(get_DCA);
    else fcn_3Dy.push_back(get_DCA_BG);
    switch(i) {
      case 0:  fcn_3Dz.push_back(get_pt_V0); break;
      case 1:  fcn_3Dz.push_back(get_pt_V1); break;
      case 2:  fcn_3Dz.push_back(get_pt_V2); break;
      case 3:  fcn_3Dz.push_back(get_pt_V3); break;
      case 4:  fcn_3Dz.push_back(get_pt_V4); break;
      case 5:  fcn_3Dz.push_back(get_pt_V5); break;
      case 6:  fcn_3Dz.push_back(get_pt_V6); break;
      case 7:  fcn_3Dz.push_back(get_pt_V7); break;
      case 8:  fcn_3Dz.push_back(get_pt_V8); break;
      case 9:  fcn_3Dz.push_back(get_pt_V9); break;
      case 10: fcn_3Dz.push_back(get_pt_V10); break;
      case 11: fcn_3Dz.push_back(get_pt_V11); break;
      case 12: fcn_3Dz.push_back(get_pt_V12); break;
      case 13: fcn_3Dz.push_back(get_pt_V13); break;
      case 14: fcn_3Dz.push_back(get_pt_V14); break;
      case 15: fcn_3Dz.push_back(get_pt_V15); break;
      case 16: fcn_3Dz.push_back(get_pt_V16); break;
      case 17: fcn_3Dz.push_back(get_pt_V17); break;
      case 18: fcn_3Dz.push_back(get_pt_V18); break;
      case 19: fcn_3Dz.push_back(get_pt_V19); break;
      case 20: fcn_3Dz.push_back(get_pt_V20); break;
      case 21: fcn_3Dz.push_back(get_pt_V21); break;
      case 22: fcn_3Dz.push_back(get_pt_V22); break;
      case 23: fcn_3Dz.push_back(get_pt_V23); break;
      case 24: fcn_3Dz.push_back(get_pt_V24); break;
      case 25: fcn_3Dz.push_back(get_pt_V25); break;
      case 26: fcn_3Dz.push_back(get_pt_V26); break;
      case 27: fcn_3Dz.push_back(get_pt_V27); break;
      case 28: fcn_3Dz.push_back(get_pt_V28); break;
      case 29: fcn_3Dz.push_back(get_pt_V29); break;
      case 30: fcn_3Dz.push_back(get_pt_V30); break;
      case 31: fcn_3Dz.push_back(get_pt_V31); break;
      default: fcn_3Dz.push_back(get_pt_V0); break;
    }
    fcn_3D_weight.push_back(set_weight_1);
  }
  for (int i = 0; i < 24; ++i) {
    Master3D.push_back(new TH3D(Form("delt_phi_ee_DCA_V%d", i), Form("delt_phi_ee_DCA_V%d;#Delta#phi_{ee} [rad];m_{ee} [GeV];p_{T} [GeV]", i), 32, 0, 3.2, 45, 0, 4.5,  12, 0, 4.8));
    fcn_3Dx.push_back(get_phi_ee);
    fcn_3Dy.push_back(get_mass_ee);
    switch(i) {
      case 0:  fcn_3Dz.push_back(get_pt_V0); break;
      case 1:  fcn_3Dz.push_back(get_pt_V1); break;
      case 2:  fcn_3Dz.push_back(get_pt_V2); break;
      case 3:  fcn_3Dz.push_back(get_pt_V3); break;
      case 4:  fcn_3Dz.push_back(get_pt_V4); break;
      case 5:  fcn_3Dz.push_back(get_pt_V5); break;
      case 6:  fcn_3Dz.push_back(get_pt_V6); break;
      case 7:  fcn_3Dz.push_back(get_pt_V7); break;
      case 8:  fcn_3Dz.push_back(get_pt_V8); break;
      case 9:  fcn_3Dz.push_back(get_pt_V9); break;
      case 10: fcn_3Dz.push_back(get_pt_V10); break;
      case 11: fcn_3Dz.push_back(get_pt_V11); break;
      case 12: fcn_3Dz.push_back(get_pt_V12); break;
      case 13: fcn_3Dz.push_back(get_pt_V13); break;
      case 14: fcn_3Dz.push_back(get_pt_V14); break;
      case 15: fcn_3Dz.push_back(get_pt_V15); break;
      case 16: fcn_3Dz.push_back(get_pt_V16); break;
      case 17: fcn_3Dz.push_back(get_pt_V17); break;
      case 18: fcn_3Dz.push_back(get_pt_V18); break;
      case 19: fcn_3Dz.push_back(get_pt_V19); break;
      case 20: fcn_3Dz.push_back(get_pt_V20); break;
      case 21: fcn_3Dz.push_back(get_pt_V21); break;
      case 22: fcn_3Dz.push_back(get_pt_V22); break;
      case 23: fcn_3Dz.push_back(get_pt_V23); break;
      case 24: fcn_3Dz.push_back(get_pt_V24); break;
      case 25: fcn_3Dz.push_back(get_pt_V25); break;
      case 26: fcn_3Dz.push_back(get_pt_V26); break;
      case 27: fcn_3Dz.push_back(get_pt_V27); break;
      case 28: fcn_3Dz.push_back(get_pt_V28); break;
      case 29: fcn_3Dz.push_back(get_pt_V29); break;
      case 30: fcn_3Dz.push_back(get_pt_V30); break;
      case 31: fcn_3Dz.push_back(get_pt_V31); break;
      default: fcn_3Dz.push_back(get_pt_V0); break;
    }
    fcn_3D_weight.push_back(set_weight_1);
  }

  Master3D.push_back(new TH3D("mass_dca_pt_upsilon", "mass_dca_pt_upsilon;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]",  100, 8, 12, 10, 0, 1000, 10, 0, 10));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V32);
  fcn_3D_weight.push_back(set_weight_1);
  
  Master3D.push_back(new TH3D("mass_DCA1_DCA2_v0" , "mass_DCA1_DCA2_v0;m_{ee} [GeV];DCA_{e} [#mum];DCA_{e} [#mum]" , 90, 0, 4.5, 20, -500, 500, 20, -500, 500));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA1);
  fcn_3Dz.push_back(get_DCA2);
  fcn_3D_weight.push_back(set_weight_1);

  Master3D.push_back(new TH3D("mass_dalitz_v0" , "mass_dalitz_v0;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V11);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v1" , "mass_dalitz_v1;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V15);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v2" , "mass_dalitz_v2;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V7);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v3" , "mass_dalitz_v3;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V9);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v4" , "mass_dalitz_v4;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V10);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v5" , "mass_dalitz_v5;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V12);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v6" , "mass_dalitz_v6;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V17);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v7" , "mass_dalitz_v7;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V2);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v8" , "mass_dalitz_v8;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA_BG);
  fcn_3Dz.push_back(get_pt_V20);
  fcn_3D_weight.push_back(set_weight_1);
  Master3D.push_back(new TH3D("mass_dalitz_v9" , "mass_dalitz_v9;m_{ee} [GeV];DCA_{ee} [#mum];p_{T} [GeV]" , 200, 0, 0.5, 40, 0, 1000, 10, 0, 5));
  fcn_3Dx.push_back(get_mass_ee);
  fcn_3Dy.push_back(get_DCA);
  fcn_3Dz.push_back(get_pt_V18);
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
  //const double DCA = TMath::Sqrt( TMath::Abs ( DCA1 * DCA1 + DCA2 * DCA2 ) );
  //if (DCA > 1000) return 995;
  //if (DCA < -3000) return -3000;
  //double sign = (DCA1 + DCA2 >= 0 ? 1.0 : -1.0);
  const double DCA = TMath::Sqrt( /*0.5**/  TMath::Abs ( DCA1 * DCA1 + DCA2 * DCA2 ) );
  if (DCA > 990 ) return 990;
  //if (DCA < -490 ) return -490;
  
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

float Run14AuAuLeptonCombyHistos::get_DCA_BG(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  double DCA1 = p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  double DCA2 = p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
  
  //const double DCA = TMath::Sqrt( TMath::Abs ( 2*DCA1 * DCA1 + 2*DCA2 * DCA2 ) );
  //if (DCA > 3000) return 3000;
  //if (DCA < -3000) return -3000;
  //double sign = (DCA1 + DCA2 >= 0 ? 1.0 : -1.0);
  const double DCA = TMath::Sqrt( /*0.5**/  TMath::Abs ( DCA1 * DCA1 + DCA2 * DCA2 ) );
  if (DCA > 990 ) return 990;
  //if (DCA < -490 ) return -490;
  
  return DCA;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7; // Limit pt to 5 GeV/c
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V3(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)//26
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -999;

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
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt < 4.7 ? pt : 4.7;
}
float Run14AuAuLeptonCombyHistos::get_pt_V4(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 || hit_assoc2 < 10000 ) return -999;

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
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -999;

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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

  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V17(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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
  
  return pt < 4.7 ? pt : 4.7;
}


float Run14AuAuLeptonCombyHistos::get_pt_V18(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
  const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);  

  if( hadron_reject1 < 10000 || hadron_reject2 < 10000 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V19(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int centrality1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CENTRALITY);
  const int centrality2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CENTRALITY);

  if ( centrality1 < 10 || centrality2 < 10 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V20(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)//27
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  
  if (! ( (conv_reject1 == -10 && conv_reject2 > 999) || (conv_reject2 == -10 && conv_reject1 > 999) ) ) return -999;

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 && hit_assoc2 < 100 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}


float Run14AuAuLeptonCombyHistos::get_pt_V21(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)//29
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  
  if (! ( (conv_reject1 == -10 && conv_reject2 > 999) || (conv_reject2 == -10 && conv_reject1 > 999) ) ) return -999;

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 10000 && hit_assoc2 < 10000 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V22(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)//30
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V23(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)//32
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if  ( conv_reject1 > -1 || conv_reject2 > -1 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V24(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V25(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);

  const int centrality1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CENTRALITY);
  const int centrality2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CENTRALITY);

  if ( centrality1 < 10 || centrality2 < 10 ) return -999;

  const int bbcq = p1->get_integer(Run14AuAuLeptonCombyEnum::BBCQ);
  const int bbcq2 = p2->get_integer(Run14AuAuLeptonCombyEnum::BBCQ);
  if ( centrality1 > 50 && bbcq > 200 ) return -999;
  if ( centrality2 > 50 && bbcq2 > 200 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V26(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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
  
  return pt < 4.7 ? pt : 4.7;
}
float Run14AuAuLeptonCombyHistos::get_pt_V27(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V28(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V29(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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

  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V30(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V31(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  
  if (! ( (conv_reject1 < 0 && conv_reject2 > 999) || (conv_reject2 < 0 && conv_reject1 > 999) ) ) return -999;

  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 && hit_assoc2 < 100 ) return -999;
  
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
  
  return pt < 4.7 ? pt : 4.7;
}

float Run14AuAuLeptonCombyHistos::get_pt_V32(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const double px = p1->get_px() + p2->get_px();
  const double py = p1->get_py() + p2->get_py();

  const double pt = TMath::Sqrt(px*px + py*py);
  
  return pt < 4.7 ? pt : 4.7;
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

float Run14AuAuLeptonCombyHistos::get_DCA1(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
  const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

  if ( conv_reject1 < 1000 || conv_reject2 < 1000 ) return -9999;
  
  const int hit_assoc1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);
  const int hit_assoc2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC);

  if ( hit_assoc1 < 100 || hit_assoc2 < 100 ) return -9999;
  
  const int ghost1 = p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST);
  const int ghost2 = p2->get_integer(Run14AuAuLeptonCombyEnum::GHOST);

  if ( ghost1 > 0 || ghost2 > 0 ) return -9999;

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
      return -9999;

  if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
      return -9999;

  if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
      return -9999;

  return p1->get_double(Run14AuAuLeptonCombyEnum::DCAX);
}

float Run14AuAuLeptonCombyHistos::get_DCA2(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p2 = ct2->GetTrack(i2);

  return p2->get_double(Run14AuAuLeptonCombyEnum::DCAX);
}


float Run14AuAuLeptonCombyHistos::set_weight_1(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
  UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
  UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

  UltraLightTrack *p1 = ct1->GetTrack(i1);
  UltraLightTrack *p2 = ct2->GetTrack(i2);
  
  const double weigth1 = p1->get_double(Run14AuAuLeptonCombyEnum::WEIGHT);
  const double weigth2 = p2->get_double(Run14AuAuLeptonCombyEnum::WEIGHT);

  return weigth1*weigth2;
}
