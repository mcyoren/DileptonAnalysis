#include "Run14AuAuLeptonCombyCutter.h"

using namespace std;
using namespace findNode;

Run14AuAuLeptonCombyCutter::Run14AuAuLeptonCombyCutter(const char *name) : cbMasterCutter(name) {}

bool Run14AuAuLeptonCombyCutter::isEventOK(PHCompositeNode *topNode)
{
    PHGlobal *global = getClass<PHGlobal>(topNode, "PHGlobal");

    if (!global)
    {
        std::cout << "ERROR: No global  node!" << std::endl;
        return false;
    }

    const float bbc_vertex = global->getBbcZVertex();
    if (TMath::Abs(bbc_vertex) > 10)
        return false;

    return true;
}

bool Run14AuAuLeptonCombyCutter::isParticleType1(
    PHParticle *part, const unsigned int iTrack)
{
    UltraLight *ul = dynamic_cast<UltraLight *>(part);
    UltraLightTrack *ult = ul->GetTrack(iTrack);

    const int particleType =
        ult->get_integer(Run14AuAuLeptonCombyEnum::PTYPE);

    return (particleType == 1);
}

bool Run14AuAuLeptonCombyCutter::isParticleType2(
    PHParticle *part, const unsigned int iTrack)
{
    UltraLight *ul = dynamic_cast<UltraLight *>(part);
    UltraLightTrack *ult = ul->GetTrack(iTrack);

    const int particleType =
        ult->get_integer(Run14AuAuLeptonCombyEnum::PTYPE);

    return (particleType == 2);
}

bool Run14AuAuLeptonCombyCutter::isCollectionOK(
    PHParticle *type1, PHParticle *type2)
{
    return true;
}

bool Run14AuAuLeptonCombyCutter::isPairOK(
    PHParticle *type1, const unsigned int iTrack1,
    PHParticle *type2, const unsigned int iTrack2)
{
    UltraLight *ct1 = dynamic_cast<UltraLight *>(type1);
    UltraLight *ct2 = dynamic_cast<UltraLight *>(type2);

    UltraLightTrack *p1 = ct1->GetTrack(iTrack1);
    UltraLightTrack *p2 = ct2->GetTrack(iTrack2);

    const float phi_pip = p1->get_double(Run14AuAuLeptonCombyEnum::PHI);
    const float phi_pim = p2->get_double(Run14AuAuLeptonCombyEnum::PHI);

    // if (phi_pip < phi_pim)
    //     return false; //-- kak Run5pp

    // Crosses by Sarah
    const float alpha_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ALPHA);
    const float alpha_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ALPHA);

    const float zvtx_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ZVTX);
    const float zvtx_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ZVTX);

    const float zed_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ZED) - zvtx_pip;
    const float zed_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ZED) - zvtx_pim;

    const float dalpha = alpha_pip - alpha_pim;
    const float dphi = phi_pip - phi_pim;
    const float dzed = zed_pip - zed_pim;

    if (GetRichGhost(p1->get_double(Run14AuAuLeptonCombyEnum::CRKPHI), p2->get_double(Run14AuAuLeptonCombyEnum::CRKZED),
                     p2->get_double(Run14AuAuLeptonCombyEnum::CRKPHI), p2->get_double(Run14AuAuLeptonCombyEnum::CRKZED)) < 5. )
        return false;

    //if( TMath::Abs(zvtx_pip - zvtx_pim) > 1.0 )
    //    return false;

    // pc1
    if (TMath::Abs(dzed) < 6.0 && TMath::Abs(dphi - (0.13 * dalpha)) < 0.015)
        return false;
    // x1x2_1
    if (TMath::Abs(dphi - (0.04 * dalpha)) < 0.015)
        return false;
    // x1x2_2
    if (TMath::Abs(dphi - (-0.065 * dalpha)) < 0.015)
        return false;

    const float psi_pip = p1->get_double(Run14AuAuLeptonCombyEnum::PSI);
    const float psi_pim = p2->get_double(Run14AuAuLeptonCombyEnum::PSI);
    //if ( TMath::Abs(psi_pip-psi_pim)>TMath::Pi()/16. && TMath::Abs(psi_pip-psi_pim)<15.*TMath::Pi()/16. )
    //    return false;

    
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


    const int id11 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
    const int id12 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
    const int id13 = p1->get_integer(Run14AuAuLeptonCombyEnum::ID3); 
    const int id21 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID1); 
    const int id22 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID2); 
    const int id23 = p2->get_integer(Run14AuAuLeptonCombyEnum::ID3); 

    if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
       std::cout<<id11<< " "<<id21<<" "<<psi_pip<<" "<<psi_pim<<std::endl;

    if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
       std::cout<<id12<< " "<<id22<<" "<<psi_pip<<" "<<psi_pim<<std::endl;

    if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
       std::cout<<id13<< " "<<id23<<" "<<psi_pip<<" "<<psi_pim<<std::endl;


    if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.017 )
        return false;

    if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.0085 )
        return false;

    if ( TMath::Abs(phi13-phi23)<0.0008 && TMath::Abs(the13-the23)<0.01 )
        return false;
    return true;
}

int Run14AuAuLeptonCombyCutter::Cleaner(PHParticle *Part)
{
    return 0;
}

void Run14AuAuLeptonCombyCutter::CrossClean(
    PHParticle *type1, PHParticle *type2)
{
    return;
}

float Run14AuAuLeptonCombyCutter::GetRichGhost(float phi1, float z1, float phi2, float z2)
{
    const float dcenter_phi_sigma = 0.013;
    const float dcenter_z_sigma = 5.0;

    const float dcenter_z = (z1 - z2) / dcenter_z_sigma;
    const float dcenter_phi = (phi1 - phi2) / dcenter_phi_sigma;

    return sqrt(dcenter_phi * dcenter_phi + dcenter_z * dcenter_z);
}