#include "Run14AuAuLeptonCombyCutter.h"

using namespace std;
using namespace findNode;

Run14AuAuLeptonCombyCutter::Run14AuAuLeptonCombyCutter(const char *name) : cbMasterCutter(name) {}

bool Run14AuAuLeptonCombyCutter::isEventOK(PHCompositeNode *topNode)
{
    PHGlobal *global = getClass<PHGlobal>(topNode, "PHGlobal");

    if (!global)
    {
        std::cout << "ERROR: No global node!" << std::endl;
        return false;
    }

    const float bbc_vertex = global->getBbcZVertex();
    if (TMath::Abs(bbc_vertex) > 10)
        return false;

    return true;
}

bool Run14AuAuLeptonCombyCutter::isParticleType1(PHParticle *part, const unsigned int iTrack)
{
    UltraLight *ul = dynamic_cast<UltraLight *>(part);
    UltraLightTrack *ult = ul->GetTrack(iTrack);

    const int particleType = ult->get_integer(Run14AuAuLeptonCombyEnum::PTYPE);

    return (particleType == 1);
}

bool Run14AuAuLeptonCombyCutter::isParticleType2(PHParticle *part, const unsigned int iTrack)
{
    UltraLight *ul = dynamic_cast<UltraLight *>(part);
    UltraLightTrack *ult = ul->GetTrack(iTrack);

    const int particleType = ult->get_integer(Run14AuAuLeptonCombyEnum::PTYPE);

    return (particleType == 2);
}

void Run14AuAuLeptonCombyCutter::CrossClean(PHParticle *type1, PHParticle *type2)
{
    UltraLight *php1 = dynamic_cast<UltraLight *>(type1);
    UltraLight *php2 = dynamic_cast<UltraLight *>(type2);

    unsigned int nT1 = php1->get_npart();
    unsigned int nT2 = php2->get_npart();

    // Loop over ALL the possible pairs in an event and check to see that EACH
    // pair is OK. Even if only one pair isn't, then the whole collection fails.
    //_______________________Check the Type1-Type2 Combinations_________________
    for (unsigned int it1 = 0; it1 < nT1; it1++)
    {
        NulifyGhost(php1, it1);
    }
    for (unsigned int it2 = 0; it2 < nT2; it2++)
    {
        NulifyGhost(php2, it2);
    }
    if (nT1 > 0 && nT2 > 0)
    {
        for (unsigned int it1 = 0; it1 < nT1; it1++)
        {
            for (unsigned int it2 = 0; it2 < nT2; it2++)
            {
                // Ring sharing cut
                int dcenter = getDcenter_phi(php1, it1, php2, it2);
                int EMCdistance = getEMCdistance(php1, it1, php2, it2);
                int dPC1 = getDC_ghost(php1, it1, php2, it2);

                if (dcenter || EMCdistance || dPC1)
                {
                    int del_type = ChooseBest(php1, it1, php2, it2);
                    if (del_type == 1)
                    {
                        nT1--;
                        it1--;
                        break;
                    }
                    else if (del_type == 2)
                    {
                        nT2--;
                        it2--;
                    }
                }
            }
        }
    }

    //_______________________Check the Type1-Type1 Combinations_________________
    if (nT1 > 0)
    {
        for (unsigned int it1 = 0; it1 < nT1 - 1; it1++)
        {
            for (unsigned int it2 = it1 + 1; it2 < nT1; it2++)
            {
                int dcenter = getDcenter_phi(php1, it1, php1, it2);
                int EMCdistance = getEMCdistance(php1, it1, php1, it2);
                int dPC1 = getDC_ghost(php1, it1, php1, it2);

                if (dcenter || EMCdistance || dPC1)
                {
                    int del_type = ChooseBest(php1, it1, php1, it2);
                    if (del_type == 1)
                    {
                        nT1--;
                        it1--;
                        break;
                    }
                    else if (del_type == 2)
                    {
                        nT1--;
                        it2--;
                    }
                }
            }
        }
    }

    //_______________________Check the Type2-Type2 Combinations_________________
    if (nT2 > 0)
    {
        for (unsigned int it1 = 0; it1 < nT2 - 1; it1++)
        {
            for (unsigned int it2 = it1 + 1; it2 < nT2; it2++)
            {
                int dcenter = getDcenter_phi(php2, it1, php2, it2);
                int EMCdistance = getEMCdistance(php2, it1, php2, it2);
                int dPC1 = getDC_ghost(php2, it1, php2, it2);

                if (dcenter || EMCdistance || dPC1)
                {
                    int del_type = ChooseBest(php2, it1, php2, it2);
                    if (del_type == 1)
                    {
                        nT2--;
                        it1--;
                        break;
                    }
                    else if (del_type == 2)
                    {
                        nT2--;
                        it2--;
                    }
                }
            }
        }
    }
    return;
}

bool Run14AuAuLeptonCombyCutter::isCollectionOK(PHParticle *type1, PHParticle *type2)
{
    return true;
}

bool Run14AuAuLeptonCombyCutter::isPairOK(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
    int dcenter = getDcenter_phi(Type1, i1, Type2, i2);
    int EMCdistance = getEMCdistance(Type1, i1, Type2, i2);
    int dPC1 = getDC_ghost(Type1, i1, Type2, i2);

    if (dcenter || EMCdistance || dPC1) return false;

    return true;
}

int Run14AuAuLeptonCombyCutter::Cleaner(PHParticle *Part)
{
    return 0;
}

float Run14AuAuLeptonCombyCutter::GetRichGhost(float phi1, float z1, float phi2, float z2)
{
    const float dcenter_phi_sigma = 0.013;
    const float dcenter_z_sigma = 5.0;

    const float dcenter_z = (z1 - z2) / dcenter_z_sigma;
    const float dcenter_phi = (phi1 - phi2) / dcenter_phi_sigma;

    return sqrt(dcenter_phi * dcenter_phi + dcenter_z * dcenter_z);
}

int Run14AuAuLeptonCombyCutter::getDcenter_phi(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
    UltraLight *php1 = dynamic_cast<UltraLight *>(Type1);
    UltraLight *php2 = dynamic_cast<UltraLight *>(Type2);

    float phi1 = php1->GetTrack(i1)->get_double(Run14AuAuLeptonCombyEnum::CRKPHI);
    float phi2 = php2->GetTrack(i2)->get_double(Run14AuAuLeptonCombyEnum::CRKPHI);

    return TMath::Abs(phi1 - phi2) / 0.013 < 5.0;
}

int Run14AuAuLeptonCombyCutter::getEMCdistance(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
    UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
    UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

    UltraLightTrack *p1 = ct1->GetTrack(i1);
    UltraLightTrack *p2 = ct2->GetTrack(i2);

    const int Nsect1 = p1->get_integer(Run14AuAuLeptonCombyEnum::SECTOR);
    const int Ysect1 = p1->get_integer(Run14AuAuLeptonCombyEnum::YSECT);
    const int Zsect1 = p1->get_integer(Run14AuAuLeptonCombyEnum::ZSECT);
    const int Nsect2 = p2->get_integer(Run14AuAuLeptonCombyEnum::SECTOR);
    const int Ysect2 = p2->get_integer(Run14AuAuLeptonCombyEnum::YSECT);
    const int Zsect2 = p2->get_integer(Run14AuAuLeptonCombyEnum::ZSECT);

    if (Nsect1 == Nsect2 && TMath::Abs(Ysect1 - Ysect2) < 3 && TMath::Abs(Zsect1 - Zsect2) < 3)
        return true;

    return false;
}

int Run14AuAuLeptonCombyCutter::getDC_ghost(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
    UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
    UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

    UltraLightTrack *p1 = ct1->GetTrack(i1);
    UltraLightTrack *p2 = ct2->GetTrack(i2);

    const float phi_pip = p1->get_double(Run14AuAuLeptonCombyEnum::PHI);
    const float phi_pim = p2->get_double(Run14AuAuLeptonCombyEnum::PHI);

    const float alpha_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ALPHA);
    const float alpha_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ALPHA);

    const float zvtx_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ZVTX);
    const float zvtx_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ZVTX);

    const float zed_pip = p1->get_double(Run14AuAuLeptonCombyEnum::ZED) - zvtx_pip;
    const float zed_pim = p2->get_double(Run14AuAuLeptonCombyEnum::ZED) - zvtx_pim;

    const float dalpha = alpha_pip - alpha_pim;
    const float dphi = phi_pip - phi_pim;
    const float dzed = zed_pip - zed_pim;

    // pc1
    if (TMath::Abs(dzed) < 6.0 && TMath::Abs(dphi - (0.13 * dalpha)) < 0.015)
        return true;
    // x1x2_1
    if (TMath::Abs(dphi - (0.04 * dalpha)) < 0.015)
        return true;
    // x1x2_2
    if (TMath::Abs(dphi - (-0.065 * dalpha)) < 0.015)
        return true;

    return false;
}

int Run14AuAuLeptonCombyCutter::ChooseBest(PHParticle *Type1, const unsigned int i1, PHParticle *Type2, const unsigned int i2)
{
    UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
    UltraLight *ct2 = dynamic_cast<UltraLight *>(Type2);

    UltraLightTrack *p1 = ct1->GetTrack(i1);
    UltraLightTrack *p2 = ct2->GetTrack(i2);

    const int conv_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);
    const int conv_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT);

    if (conv_reject1 < 0 && conv_reject2 < 0)
    {   
        const int hadron_reject1 = p1->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);
        const int hadron_reject2 = p2->get_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT);

        if (p1->get_double(Run14AuAuLeptonCombyEnum::ZVTX) == p2->get_double(Run14AuAuLeptonCombyEnum::ZVTX))
        {
            if (hadron_reject1 < hadron_reject2)
            {
                ct1->RemoveTrack(i1);
                return 1;
            }
            else
            {
                ct2->RemoveTrack(i2);
                return 2;
            }
        }
        if (hadron_reject1 < hadron_reject2)
            p1->set_integer(Run14AuAuLeptonCombyEnum::GHOST, 1);
        else
            p2->set_integer(Run14AuAuLeptonCombyEnum::GHOST, 1);

    }

    if (conv_reject1 >= 0 && conv_reject2 >= 0)
    {   
        if (p1->get_double(Run14AuAuLeptonCombyEnum::ZVTX) == p2->get_double(Run14AuAuLeptonCombyEnum::ZVTX))
        {
            if (conv_reject1 < conv_reject2)
            {
                ct1->RemoveTrack(i1);
                return 1;
            }
            else
            {
                ct2->RemoveTrack(i2);
                return 2;
            }
        }
        if (conv_reject1 < conv_reject2)
            p1->set_integer(Run14AuAuLeptonCombyEnum::GHOST, 1);
        else
            p2->set_integer(Run14AuAuLeptonCombyEnum::GHOST, 1);
    }


    return 0;
}

int Run14AuAuLeptonCombyCutter::NulifyGhost(PHParticle *Type1, const unsigned int i1)
{
    UltraLight *ct1 = dynamic_cast<UltraLight *>(Type1);
    UltraLightTrack *p1 = ct1->GetTrack(i1);
    if (p1->get_integer(Run14AuAuLeptonCombyEnum::GHOST))
        p1->set_integer(Run14AuAuLeptonCombyEnum::GHOST, 0);
    return 0;
}