#ifndef __RUN14AUAULEPTONCOMBYCUTTER_H__
#define __RUN14AUAULEPTONCOMBYCUTTER_H__

#include "CabanaBoy/cbMasterCutter.h"
#include "PHCompositeNode.h"
#include "PHParticle.h"
#include "UltraLight/UltraLight.h"
#include "UltraLight/UltraLightTrack.h"
#include "PHGlobal.h"
#include "getClass.h"
#include "Run14AuAuLeptonCombyEnum.h"
#include "Run14AuAuLeptonCombyConstants.h"

#include "TVector3.h"
#include "TMath.h"

#include <cmath>

class PHCompositeNode;
class PHParticle;

class Run14AuAuLeptonCombyCutter: public cbMasterCutter {

public:
    Run14AuAuLeptonCombyCutter(const char* name = "NONAME");
    virtual ~Run14AuAuLeptonCombyCutter() {};

    bool isEventOK(PHCompositeNode *topNode);
    bool isParticleType1(PHParticle *part, const unsigned int iTrack);
    bool isParticleType2(PHParticle *part, const unsigned int iTrack);
    bool isPairOK(PHParticle *Type1, const unsigned int TrackNumber1, PHParticle *Type2, const unsigned int TrackNumber2);
    bool isCollectionOK(PHParticle *Type1, PHParticle *Type2);
    int Cleaner(PHParticle *Part1);
    void CrossClean(PHParticle *Type1, PHParticle *Type2);
    static float GetRichGhost(float phi1, float z1, float phi2, float z2);
    //float Phi_V(const TVector3& p1, const TVector3& p2);

private:

};

#endif