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

    const int centrality_pip = p1->get_integer(Run14AuAuLeptonCombyEnum::CENTR);
    const int centrality_pim = p2->get_integer(Run14AuAuLeptonCombyEnum::CENTR);
    const float zvtx_pip =     p1->get_double(Run14AuAuLeptonCombyEnum::ZVTX);
    const float zvtx_pim =     p2->get_double(Run14AuAuLeptonCombyEnum::ZVTX);
    const float psi_pip =      p1->get_double(Run14AuAuLeptonCombyEnum::PSI);
    const float psi_pim =      p2->get_double(Run14AuAuLeptonCombyEnum::PSI);

    if(centrality_pip<20||centrality_pim<20)
    {
        if( TMath::Abs(centrality_pip-centrality_pim)>5 ) return false;
        if( TMath::Abs(zvtx_pip-zvtx_pim)>2.5 ) return false;
        if( TMath::Abs(psi_pip-psi_pim)>TMath::Pi()/8 && TMath::Abs(psi_pip-psi_pim)<7*TMath::Pi()/8) return false;
    }
    else if(centrality_pip<60||centrality_pim<60)
    {
        if( TMath::Abs(centrality_pip-centrality_pim)>10 ) return false;
        if( TMath::Abs(zvtx_pip-zvtx_pim)>5 ) return false;
        if( TMath::Abs(psi_pip-psi_pim)>TMath::Pi()/4 && TMath::Abs(psi_pip-psi_pim)<3*TMath::Pi()/4) return false;
    }
    else{
        if( TMath::Abs(centrality_pip-centrality_pim)>20 ) return false;
    }

    const float phi_pip = p1->get_double(Run14AuAuLeptonCombyEnum::PHI);
    const float phi_pim = p2->get_double(Run14AuAuLeptonCombyEnum::PHI);

    // if (phi_pip < phi_pim)
    //     return false; //-- kak Run5pp

    // Crosses by Sarah
    const float alpha_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ALPHA);
    const float alpha_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ALPHA);

    const float zed_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ZED);
    const float zed_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ZED);

    const float dalpha = alpha_pip - alpha_pim;
    const float dphi = phi_pip - phi_pim;
    const float dzed = zed_pip - zed_pim;

    if (GetRichGhost(p1->get_double(Run14AuAuLeptonCombyEnum::CRKPHI), p2->get_double(Run14AuAuLeptonCombyEnum::CRKZED),
                     p2->get_double(Run14AuAuLeptonCombyEnum::CRKPHI), p2->get_double(Run14AuAuLeptonCombyEnum::CRKZED)) < 3 )
        return false;
    if (GetRichGhost(p1->get_double(Run14AuAuLeptonCombyEnum::CRKPHI), p1->get_double(Run14AuAuLeptonCombyEnum::CRKZED),
                     p2->get_double(Run14AuAuLeptonCombyEnum::CRKPHI), p2->get_double(Run14AuAuLeptonCombyEnum::CRKZED)) < 5.0 )
        return false;

    // pc1
    if (fabs(dzed) < 6.0 && fabs(dphi - (0.13 * dalpha)) < 0.015)
        return false;
    // x1x2_1
    if (fabs(dphi - (0.04 * dalpha)) < 0.015)
        return false;
    // x1x2_2
    if (fabs(dphi - (-0.065 * dalpha)) < 0.015)
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