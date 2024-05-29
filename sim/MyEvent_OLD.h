#ifndef __MYEVENT_OLD_H__
#define __MYEVENT_OLD_H__

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MyEvent                                                              //
//                                                                      //
// Description of the event and track parameters                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include <climits>

// using namespace std
//   UShort_t  USHORT_MAX = 65535;
//   Byte_t    UCHAR_MAX   = 255;
//   float Max_Float  = -99999;

class TDirectory;

namespace DileptonAnalysis
{

   class MyTrack : public TObject
   {

   private:

      int    trkid;         // DC track index 
      int    trkinfo;                           // arm: 0 east, 1 west 
                                 // side 0 south, 1 north
                                 // sector 0-7 with 0 being lowest west sector

// basic parameters from tracking
      int      q;                // charge of the track
      float   pt;           // transverse momentum
      float   phi_DC;        // phi of track at 220 cm
      float   z_DC;          // z poition of track at 220 cm (measured by PC1)
      float   phi0;          // phi at interaction point - assuming track is primary 
      float   the0;          // dito but theta
      float   alpha;         // alpha inclination at DC (total field bent = alpha + (phi_DC-phi0))
      int     arm;
      int     dcside;
      int     sect;
      float   chi2;
      float   qpt;
      int     n_hits;

// Updated track parameters
      
      int   q_prime;
      float phi0_prime;
      float the0_prime;
      float pt_prime;
      float alpha_prime;

// Track Hit Association

      int hit_index_l1;
      int hit_index_l2;
      int hit_index_l3;
      int hit_index_l4;
      int hit_index_l5;
      int hit_index_l6;
      int hit_index_l7;
      int hit_index_l8;
      int hit_counter_l1;
      int hit_counter_l2;
      int hit_counter_l3;
      int hit_counter_l4;
      int hit_counter_l5;
      int hit_counter_l6;
      int hit_counter_l7;
      int hit_counter_l8; 

// RICH variables
      float   crkphi;        // phi position of ring
      float   crkz;          // z position of ring
      int      n0;            // number of photo tubes in mask

// EMCal variables
      int    emcid;         // EMC cluster 
      int    emctower;      // emctower = 1000*iy + 10*iz + sector
                                 // EMC sector 0-4 from bottom of arm
                                 // iy tower coordinates in rphi
                                 // iz tower coordinates in z
      float   ecore;         // energy of cluster - assuning its a photon cluster
      float   prob;          // similar to chi2 (is possibly redundant)
      float   emcdz;         // delta z of cluster and track projection
      float   emcdphi;       // delta phi of cluster and track projection
      float   emctof;        // TOF of EMCal (PbSc only)
      float   dep;           // sigmalized E/p

// TOF variables (should include TOF east and TOF west)
      float   tofe;          // TOF west??
// ERT information
      int      isERT;         // indicates which ERT trigger fired 

// MC variables
      int    mcid;         // associated MC track 


   public:
      MyTrack()
      {  
         trkid    = -99999;  // max value of unsigned short
         trkinfo  = -99999;
         arm      = -99999;
         dcside   = -99999;
         sect     = -99999;
         chi2     = -99999;
         emcid    = -99999; 
         emctower = -99999; 
         n0       = -99999;    // max value for byte
         isERT    = -99999;
         mcid     = -99999;
         q        = -99999;
	      qpt      = -99999;
         q_prime  = -99999;
         n_hits   = -99999;
         hit_index_l1 = -99999;
         hit_index_l2 = -99999;
         hit_index_l3 = -99999;
         hit_index_l4 = -99999;
         hit_index_l5 = -99999;
         hit_index_l6 = -99999;
         hit_index_l7 = -99999;
         hit_index_l8 = -99999;
         hit_counter_l1  = -99999;
         hit_counter_l2  = -99999;
         hit_counter_l3  = -99999;
         hit_counter_l4  = -99999;
         hit_counter_l5  = -99999;
         hit_counter_l6  = -99999;
         hit_counter_l7  = -99999;
         hit_counter_l8  = -99999;
         pt       = -99999;
         pt_prime = -99999; 
         phi_DC   = -99999;
         z_DC     = -99999;
         phi0     = -99999;
         phi0_prime = -99999;
         the0     = -99999;
         the0_prime = -99999;
         alpha    = -99999;
         alpha_prime = -99999;
         ecore    = -99999;
         dep      = -99999;
         prob     = -99999; 
         emcdz    = -99999;
         emcdphi  = -99999;
         emctof   = -99999;
         crkphi   = -99999;
         crkz     = -99999;
         tofe     = -99999;
      };

      virtual ~MyTrack() {  };

      int      GetTrkId() const { return trkid; };
      int      GetTrkInfo() const { return trkinfo; };       
      int      GetArm() const { return arm;};
      int      GetDCside() const { return dcside;};
      int      GetSect() const { return sect; };
      float    GetChi2() const { return chi2; };
      float    GetPx() const { if (pt==-99999) { return pt;} else {return pt*cos(phi0);} };
      float    GetPy() const { if (pt==-99999) { return pt;} else {return pt*sin(phi0);} };
      float    GetPz() const { if (pt==-99999) { return pt;} else {return pt*cos(the0);} };
      float    GetPt() const { return pt; };
      float    GetPtPrime() const { return pt_prime; };
      float    GetPhiDC() const { return phi_DC; };
      float    GetZDC() const { return z_DC; };
      float    GetPhi0() const { return phi0; };
      int      GetNHits() const {return n_hits; };
      float    GetPhi0Prime() const { return phi0_prime; };
      float    GetThe0() const { return the0; };
      float    GetThe0Prime() const { return the0_prime; };
      float    GetAlpha() const { return alpha; };
      float    GetAlphaPrime() const { return alpha_prime; };
      int      GetCharge() const { return q; };
      int      GetChargePrime() const { return q_prime; };
      int      GetEmcId() const { return emcid; };
      float    GetEcore() const { return ecore; };
      int      GetN0() const { return n0; };
      float    GetDep() const { return dep; };
      float    GetProb() const { return prob; };
      float    GetEmcdz() const { return emcdz; };
      float    GetEmcdphi() const { return emcdphi; };
      float    GetCrkphi() const { return crkphi; };
      float    GetCrkz() const { return crkz; };
      int      GetEsect() const { return int(emctower-10*((emctower-1000*(emctower/1000))/10)-1000*(emctower/1000)); };
      int      GetYsect() const { return int(emctower/1000); };
      int      GetZsect() const { return int((emctower-1000*(emctower/1000))/10); };
      float    GetEmcTOF() const { return emctof; };
      float    GetTOFE() const { return tofe; };
      int      GetisERT() const { return isERT; };
      int      GetMcId() const { return mcid; };
      int      GetHitIndexL1() const { return hit_index_l1; };
      int      GetHitIndexL2() const { return hit_index_l2; };
      int      GetHitIndexL3() const { return hit_index_l3; };
      int      GetHitIndexL4() const { return hit_index_l4; };
      int      GetHitIndexL5() const { return hit_index_l5; };
      int      GetHitIndexL6() const { return hit_index_l6; };
      int      GetHitIndexL7() const { return hit_index_l7; };
      int      GetHitIndexL8() const { return hit_index_l8; };
      int     GetHitCounterL1() { return hit_counter_l1;};
      int     GetHitCounterL2() { return hit_counter_l2;};
      int     GetHitCounterL3() { return hit_counter_l3;};
      int     GetHitCounterL4() { return hit_counter_l4;};
      int     GetHitCounterL5() { return hit_counter_l5;};
      int     GetHitCounterL6() { return hit_counter_l6;};
      int     GetHitCounterL7() { return hit_counter_l7;};
      int     GetHitCounterL8() { return hit_counter_l8;};
      

      void     SetTrkId(int strkid) { trkid = strkid; };
      void     SetArm(int sarm) { arm = sarm;};
      void     SetTrkInfo(int sarm, int sdcside, int ssect) { trkinfo = short(ssect+10*sdcside+100*sarm);};
      void     SetQpt(float spx, float spy, float scharge) { qpt = scharge*sqrt(spx*spx+spy*spy); };
      void     SetDCSide(int sdcside) { dcside = sdcside;};
      void     SetSect(int ssect) { sect = ssect;};
      void     SetChi2(float schi2) { chi2 = schi2;};
      void     SetPt(float spt) { pt = spt; };
      void     SetPtPrime(float spt_prime) {pt_prime = spt_prime; };
      void     SetQ(int sq) { q = sq; };
      void     SetQPrime(int sq_prime) { q_prime = sq_prime; };
      void     SetPhiDC(float sphi_DC) { phi_DC = sphi_DC; };
      void     SetZDC(float sz_DC) { z_DC = sz_DC; };
      void     SetPhi0(float sphi0) { phi0 = sphi0; };
      void     SetPhi0Prime(float sphi0_prime) { phi0_prime = sphi0_prime; };
      void     SetThe0(float sthe0) { the0 = sthe0; };
      void     SetThe0Prime(float sthe0_prime) { the0_prime = sthe0_prime; };
      void     SetAlpha(float salpha) { alpha = salpha; };
      void     SetAlphaPrime(float salpha_prime) { alpha_prime = salpha_prime; };
      void     SetEmcId(int sid) { emcid = sid; };
      void     SetEcore(float secore) { ecore = secore; };
      void     SetN0(int sn0) { n0 = sn0; };
      void     SetDep(float sdep) { dep = sdep; };
      void     SetProb(float sprob) { prob = sprob; };
      void     SetEmcdz(float semcdz) { emcdz = semcdz; };
      void     SetEmcdphi(float semcdphi) { emcdphi = semcdphi; }
      void     SetCrkphi(float scrkphi) { crkphi = scrkphi; };
      void     SetCrkz(float scrkz) { crkz = scrkz; };
      void     SetEmcTower(float sesect, float sysect, float szsect) { emctower = 1000*sysect + 10*szsect + sesect; };
      void     SetEmcTOF(float semctof) { emctof = semctof; };
      void     SetTOFE(float stofe) { tofe = stofe; };
      void     SetisERT(int sisERT) { isERT = sisERT; };
      void     SetMcId(int smcid) { mcid = smcid; };
      void     SetNHits(int sn_hits) { n_hits = sn_hits; };
      void     SetHitIndexL1(int shit_index_l1) { hit_index_l1 = shit_index_l1;};
      void     SetHitIndexL2(int shit_index_l2) { hit_index_l2 = shit_index_l2;};
      void     SetHitIndexL3(int shit_index_l3) { hit_index_l3 = shit_index_l3;};
      void     SetHitIndexL4(int shit_index_l4) { hit_index_l4 = shit_index_l4;};
      void     SetHitIndexL5(int shit_index_l5) { hit_index_l5 = shit_index_l5;};
      void     SetHitIndexL6(int shit_index_l6) { hit_index_l6 = shit_index_l6;};
      void     SetHitIndexL7(int shit_index_l7) { hit_index_l7 = shit_index_l7;};
      void     SetHitIndexL8(int shit_index_l8) { hit_index_l8 = shit_index_l8;};
      void     SetHitCounterL1(int shit_counter_l1) { hit_counter_l1 = shit_counter_l1;};
      void     SetHitCounterL2(int shit_counter_l2) { hit_counter_l2 = shit_counter_l2;};
      void     SetHitCounterL3(int shit_counter_l3) { hit_counter_l3 = shit_counter_l3;};
      void     SetHitCounterL4(int shit_counter_l4) { hit_counter_l4 = shit_counter_l4;};
      void     SetHitCounterL5(int shit_counter_l5) { hit_counter_l5 = shit_counter_l5;};
      void     SetHitCounterL6(int shit_counter_l6) { hit_counter_l6 = shit_counter_l6;};
      void     SetHitCounterL7(int shit_counter_l7) { hit_counter_l7 = shit_counter_l7;};
      void     SetHitCounterL8(int shit_counter_l8) { hit_counter_l8 = shit_counter_l8;};

      ClassDef(MyTrack,1)  
   };

   class MyVTXHit : public TObject
   {

   private:

      int      clustid;
      int      layer;
      int      ladder;
      int      sensor;

      float      xhit;           // xyz position of cluster [cm]
      float      yhit;           
      float      zhit;

   public:

      MyVTXHit()
      {  
         clustid = -99999;
         layer   = -99999;
         ladder  = -99999;
         sensor  = -99999;

         xhit = -99999;
         yhit = -99999;
         zhit = -99999;

      };
      virtual ~MyVTXHit() {  };

      int      GetClustId() const { return clustid; };
      int      GetLayer() const { return layer; };
      int      GetLadder() const { return ladder; };
      int      GetSensor() const {return sensor; };
      
      float    GetXHit() const { return xhit; };
      float    GetYHit() const { return yhit; };
      float    GetZHit() const { return zhit; };

      float    GetPhi()  const { float phi = atan2(yhit,xhit); if (phi<-TMath::ACos(-1)/2) phi+=2*TMath::ACos(-1); return phi; }; 

      void     SetClustId(int sclustid) { clustid = sclustid; };
      void     SetLayer(int slayer) { layer = slayer; };
      void     SetLadder(int sladder) { ladder = sladder; };
      void     SetSensor(int ssensor) {sensor = ssensor; };


      void     SetXHit(float sx) { xhit = sx; };
      void     SetYHit(float sy) { yhit = sy; };
      void     SetZHit(float sz) { zhit = sz; };

      ClassDef(MyVTXHit,1)  
   };


   class MyCluster
   {
   private:

      UShort_t       id;          // cluster ID
      UShort_t       tower;       // tower = 1000*iy + 10*iz + sector
                                  // EMC sector 0-4 from bottom of arm
                                  // iy tower coordinates in rphi
                                  // iz tower coordinates in z 
      UShort_t       info;        // info = 10* ID of associated track + isERT  
                                  // isERT is the ERT trigger number - if shower triggered an ERT 
      float      x;           // x position in cm
      float      y;           // y position in cm           
      float      z;           // z position in cm           
      float      ecore;       // corrected energy assuming photon cluster
      float      prob;        // prob varible - encodes probability that this is EM shower
      float      tof;         // time of flight - needs to be checked!
      float      chi2;        // match with EM shower (alternative to prob)

      UShort_t       mcid;          // MC cluster ID

   public:
      MyCluster()
      {
         id = USHRT_MAX;
         x = -99999;
         y = -99999;
         z = -99999;
         ecore = -99999; 
         prob = -99999;
         tof = -99999;
         chi2 = -99999;
         info = USHRT_MAX;
         tower = USHRT_MAX;
         mcid = USHRT_MAX;
      };
 
      virtual ~MyCluster() {  };

      int      GetArm() const { int arm; if(x<-600){arm = int(x);} else if(x<0){arm=0;}else{arm=1;} return arm;};
      int      GetID() const { return id; };
      float    GetX() const { return x; };
      float    GetY() const { return y; };
      float    GetZ() const { return z; };
      float    GetEcore() const { return ecore; };
      float    GetProb() const { return prob; };
      float    GetTof() const { return tof; };
      float    GetChi2() const { return chi2; };
      int      GetisERT() const { return info - 10*(info/10); };
      int      GetAssoTrack() const {return info/10;};
      int      GetSect() const { int sect; if(tower==USHRT_MAX) {sect=tower;} else {sect=tower-10*((tower-1000*(tower/1000))/10)-1000*(tower/1000);} return sect; };
      int      GetIY() const { int iy; if(tower==USHRT_MAX) {iy=tower;} else { iy = tower/1000; } return iy;};
      int      GetIZ() const { int iz; if(tower==USHRT_MAX) {iz=tower;} else {iz = (tower-1000*(tower/1000))/10;} return iz; };
      int      GetTower() const { return int(tower); };
      int      GetInfo() const { return (info); };
      int      GetMcId() const { return mcid; };

      void     SetID(int sid) { id = sid; };
      void     SetX(float sx) { x = sx; };
      void     SetY(float sy) { y = sy; };
      void     SetZ(float sz) { z = sz; };
      void     SetEcore(float secore) { ecore = secore; }; 
      void     SetProb(float sprob) { prob = sprob; };
      void     SetTof(float stof) { tof = stof; };
      void     SetChi2(float schi2) { chi2 = schi2; };
      void     SetTower(float sesect, float sysect, float szsect) { tower = (1000*sysect + 10*szsect + sesect); };
      void     SetInfo(int sisERT, int sassoTrack) { info = sisERT + 10*sassoTrack; };
      void     SetMcId(int smcid) { mcid = smcid; };

      ClassDef(MyCluster,1)  // A e+e- Cluster
   };


   class McTrack : public TObject
   {

   private:

      int        eventId;
      int        mcIndex;      
      int        particleid;
      int        parentid;
      int        primaryid;
      int        gen;

      float      px;           
      float      py;           
      float      pz;

      float      parentpx;           
      float      parentpy;           
      float      parentpz;

      float      primarypx;           
      float      primarypy;           
      float      primarypz;

      float      vertexx;
      float      vertexy;
      float      vertexz;

      float      parentvertexx;
      float      parentvertexy;
      float      parentvertexz;

      float      primaryvertexx;
      float      primaryvertexy;
      float      primaryvertexz;

      float      phi_DC;
      float      z_DC;
      float      phi0;
      float      the0;
      float      alpha;

   public:
      
      McTrack()
      {  
         eventId = -99999;
         mcIndex = -99999;      
         particleid = -99999;
         parentid = -99999;
         primaryid = -99999;
         gen = -99999;

         px = -99999;           
         py = -99999;           
         pz = -99999;

         parentpx = -99999;           
         parentpy = -99999;           
         parentpz = -99999;

         primarypx = -99999;           
         primarypy = -99999;           
         primarypz = -99999;

         vertexx = -99999;
         vertexy = -99999;
         vertexz = -99999;

         parentvertexx = -99999;
         parentvertexy = -99999;
         parentvertexz = -99999;

         primaryvertexx = -99999;
         primaryvertexy = -99999;
         primaryvertexz = -99999;

         phi_DC = -99999;
         z_DC = -99999;
         phi0 = -99999;
         the0 = -99999;
         alpha = -99999;
      };

      virtual ~McTrack() {  };

      int      GetEventId() const { return eventId; };
      int      GetMcIndex() const { return mcIndex; };
      int      GetParticleID() const { return particleid; };
      int      GetParentID() const { return parentid; };
      int      GetPrimaryID() const { return primaryid; };
      int      GetGen() const { return gen; };

      float    GetPx() const { return px; };
      float    GetPy() const { return py; };
      float    GetPz() const { return pz; };

      float    GetParentPx() const { return parentpx; };
      float    GetParentPy() const { return parentpy; };
      float    GetParentPz() const { return parentpz; };

      float    GetPrimaryPx() const { return primarypx; };
      float    GetPrimaryPy() const { return primarypy; };
      float    GetPrimaryPz() const { return primarypz; };

      float    GetVertexx() const { return vertexx; };
      float    GetVertexy() const { return vertexy; };
      float    GetVertexz() const { return vertexz; };

      float    GetParentVertexx() const { return parentvertexx; };
      float    GetParentVertexy() const { return parentvertexy; };
      float    GetParentVertexz() const { return parentvertexz; };

      float    GetPrimaryVertexx() const { return primaryvertexx; };
      float    GetPrimaryVertexy() const { return primaryvertexy; };
      float    GetPrimaryVertexz() const { return primaryvertexz; };

      float    GetPhiDC() const { return phi_DC; };
      float    GetZDC() const { return z_DC; };
      float    GetPhi0() const { return phi0; };
      float    GetThe0() const { return the0; };
      float    GetAlpha() const { return alpha; };


      void     SetEventId(int seventId) { eventId = seventId; };
      void     SetMcIndex(int smcIndex) { mcIndex = smcIndex; };
      void     SetParticleID(int sparticleid) { particleid = sparticleid; };
      void     SetParentID(int sparentid) { parentid = sparentid; };
      void     SetPrimaryID(int sprimaryid) { primaryid = sprimaryid; };
      void     SetGen(int sgen) { gen = sgen; };

      void     SetPx(float spx) { px = spx; };
      void     SetPy(float spy) { py = spy; };
      void     SetPz(float spz) { pz = spz; };

      void     SetParentPx(float sparentpx) { parentpx = sparentpx; };
      void     SetParentPy(float sparentpy) { parentpy = sparentpy; };
      void     SetParentPz(float sparentpz) { parentpz = sparentpz; };

      void     SetPrimaryPx(float sprimarypx) { primarypx = sprimarypx; };
      void     SetPrimaryPy(float sprimarypy) { primarypy = sprimarypy; };
      void     SetPrimaryPz(float sprimarypz) { primarypz = sprimarypz; };

      void     SetVertexx(float svertexx) { vertexx = svertexx; };
      void     SetVertexy(float svertexy) { vertexy = svertexy; };
      void     SetVertexz(float svertexz) { vertexz = svertexz; };

      void     SetParentVertexx(float sparentvertexx) { parentvertexx = sparentvertexx; };
      void     SetParentVertexy(float sparentvertexy) { parentvertexy = sparentvertexy; };
      void     SetParentVertexz(float sparentvertexz) { parentvertexz = sparentvertexz; };

      void     SetPrimaryVertexx(float sprimaryvertexx) { primaryvertexx = sprimaryvertexx; };
      void     SetPrimaryVertexy(float sprimaryvertexy) { primaryvertexy = sprimaryvertexy; };
      void     SetPrimaryVertexz(float sprimaryvertexz) { primaryvertexz = sprimaryvertexz; };

      void     SetPhiDC(float sphi_DC) { phi_DC = sphi_DC; };
      void     SetZDC(float sz_DC) { z_DC = sz_DC; };
      void     SetPhi0(float sphi0) { phi0 = sphi0; };
      void     SetThe0(float sthe0) { the0 = sthe0; };
      void     SetAlpha(float salpha) { alpha = salpha; };

      ClassDef(McTrack,1)  
   };


   class McCluster : public TObject
   {

   private:

      UShort_t        eventId;
      UShort_t        mcIndex;      
      UShort_t        particleid;
      UShort_t        parentid;
      UShort_t        gen;

      float      px;           
      float      py;           
      float      pz;

      float      impx;
      float      impy;
      float      impz;
      float      ekin;

      float      vertexx;
      float      vertexy;
      float      vertexz;

      // associated tower (largest energy deposit)
      UShort_t        towerid;

      // associated cluster (largest energy deposit)
      UShort_t       sect;
      UShort_t       emcid;
      UShort_t       iy;
      UShort_t       iz;

      float      ecore;
      float      x;           // xyz position of associated cluster [cm]
      float      y;           
      float      z;           
      float      chi2;

   public:
      
      McCluster()
      {  
         eventId = USHRT_MAX;
         mcIndex = USHRT_MAX;      
         particleid = USHRT_MAX;
         parentid = USHRT_MAX;
         gen = USHRT_MAX;

         px = -99999;           
         py = -99999;           
         pz = -99999;

         vertexx = -99999;
         vertexy = -99999;
         vertexz = -99999;

         impx = -99999;
         impy = -99999;
         impz = -99999;
         ekin = -99999;

         towerid = USHRT_MAX;

         sect  = USHRT_MAX;
         emcid = USHRT_MAX;
         iy = USHRT_MAX;
         iz = USHRT_MAX;

         ecore = -99999;
         x = -99999;
         y = -99999;
         z = -99999;
         chi2 = -99999;
         
      };

      virtual ~McCluster() {  };

      int      GetEventId() const { return eventId; };
      int      GetMcIndex() const { return mcIndex; };
      int      GetParticleID() const { return particleid; };
      int      GetParentID() const { return parentid; };
      int      GetGen() const { return gen; };

      float    GetPx() const { return px; };
      float    GetPy() const { return py; };
      float    GetPz() const { return pz; };

      float    GetImpx() const { return impx; };
      float    GetImpy() const { return impy; };
      float    GetImpz() const { return impz; };
      float    GetEkin() const { return ekin; };

      float    GetVertexx() const { return vertexx; };
      float    GetVertexy() const { return vertexy; };
      float    GetVertexz() const { return vertexz; };

      int      GetTowerID() const { return towerid; };

      float    GetEcore() const { return ecore; };
      int      GetSect() const { return sect; };
      int      GetClusterID() const { return emcid; };
      int      GetIY() const { return iy; };
      int      GetIZ() const { return iz; };
      float    GetX() const { return x; };
      float    GetY() const { return y; };
      float    GetZ() const { return z; };
      float    GetChi2() const { return chi2; };

      void     SetEventId(int seventId) { eventId = seventId; };
      void     SetMcIndex(int smcIndex) { mcIndex = smcIndex; };
      void     SetParticleID(int sparticleid) { particleid = sparticleid; };
      void     SetParentID(int sparentid) { parentid = sparentid; };
      void     SetGen(int sgen) { gen = sgen; };

      void     SetPx(float spx) { px = spx; };
      void     SetPy(float spy) { py = spy; };
      void     SetPz(float spz) { pz = spz; };

      void     SetImpx(float simpx) { impx = simpx; };
      void     SetImpy(float simpy) { impy = simpy; };
      void     SetImpz(float simpz) { impz = simpz; };
      void     SetEkin(float sekin) { ekin = sekin; };      

      void     SetVertexx(float svertexx) { vertexx = svertexx; };
      void     SetVertexy(float svertexy) { vertexy = svertexy; };
      void     SetVertexz(float svertexz) { vertexz = svertexz; };

      void     SetTowerID(int stowerid) { towerid = stowerid; };
      
      void     SetEcore(float secore) { ecore = secore; };
      void     SetSect(int ssect) { sect = ssect; };
      void     SetClusterID(int semcid) { emcid = semcid; };
      void     SetIY(int siy) { iy = siy; }; 
      void     SetIZ(int siz) { iz = siz; };
      void     SetX(float sx) { x = sx; };
      void     SetY(float sy) { y = sy; };
      void     SetZ(float sz) { z = sz; };
      void     SetChi2(float schi2) { chi2 = schi2; };

      ClassDef(McCluster,1)  
   };


   class MyPair
   {
   private:
      float phi_e;
      float phi_p;
      float theta_e;
      float theta_p;
      float r_pair;

      UShort_t id_e;
      UShort_t id_p;

   public:
      MyPair()
      {
         phi_e = -99999;
         phi_p = -99999;
         theta_e = -99999;
         theta_p = -99999;
         r_pair = -99999;
         id_e = USHRT_MAX;
         id_p = USHRT_MAX;
      };
      virtual ~MyPair() {  };

      float GetPhiElectron() const { return phi_e; }
      float GetPhiPositron() const { return phi_p; }
      float GetThetaElectron() const { return theta_e; }
      float GetThetaPositron() const { return theta_p; }
      float GetRPair() const { return r_pair; }
      float GetIDElectron() const { return id_e; }
      float GetIDPositron() const { return id_p; }

      void SetPhiElectron(float sphi_e) { phi_e = sphi_e; }
      void SetPhiPositron(float sphi_p) { phi_p = sphi_p; }
      void SetThetaElectron(float stheta_e) { theta_e = stheta_e; }
      void SetThetaPositron(float stheta_p) { theta_p = stheta_p; }
      void SetRPair(float sr_pair) { r_pair = sr_pair; }
      void SetIDElectron(int sid_e) { id_e = sid_e; }
      void SetIDPositron(int sid_p ) { id_p = sid_p; }

      ClassDef(MyPair,1)  
   };

   class MyEvent : public TObject
   {
 
   private:
      int              evtno;
      float        BBCcharge;
      float        BBCchargeN;
      float        BBCchargeS;
      float        BBCtimeN;
      float        BBCtimeS;
      float        zvertex;
      float        centrality;
      float        preciseX;
      float        preciseY;
      float        preciseZ;
      float        run_number;

      float        psi2_BBC; // BBC (2nd and 3rd)
      float        psi2_BBCn;
      float        psi2_BBCs;
      float        psi3_BBC;
      float        psi3_BBCn;
      float        psi3_BBCs;
      float        psi2_FVTXA0; // FVTXA0
      float        psi2_FVTXA0n;
      float        psi2_FVTXA0s;
      float        psi3_FVTXA0;
      float        psi3_FVTXA0n;
      float        psi3_FVTXA0s;

      std::vector<DileptonAnalysis::MyTrack> TrackList;
      std::vector<DileptonAnalysis::MyVTXHit> VTXHitList; 
      std::vector<DileptonAnalysis::MyCluster> ClusterList;
      std::vector<DileptonAnalysis::McTrack> McTrackList;
      std::vector<DileptonAnalysis::McCluster> McClusterList;
      std::vector<DileptonAnalysis::MyPair> PairList;
      
   public:
      MyEvent()
      {
         evtno      = -999;
         BBCcharge  = -999;
         BBCchargeN = -999;
         BBCchargeS = -999;
         BBCtimeN   = -999;
         BBCtimeS   = -999;
         zvertex    = -999;
         centrality = -999;
         
         preciseX = -999;
         preciseY = -999;
         preciseZ = -999;
         run_number = -999;

         psi2_BBC    = -999;
         psi2_BBCn   = -999;
         psi2_BBCs   = -999;
         psi3_BBC    = -999;
         psi3_BBCn   = -999;
         psi3_BBCs   = -999;
         psi2_FVTXA0    = -999;
         psi2_FVTXA0n   = -999;
         psi2_FVTXA0s   = -999;
         psi3_FVTXA0    = -999;
         psi3_FVTXA0n   = -999;
         psi3_FVTXA0s   = -999;

      };
      virtual ~MyEvent() {};

      void       ClearEvent();

      void       SetEvtNo(int sevtno) { evtno = sevtno; };
      int        GetEvtNo() { return evtno; };

      void       SetRunNumber(int srun_number) { run_number = srun_number; };
      int        GetRunNumber() { return run_number; };
      
      void       SetBBCcharge(float sBBCcharge) { BBCcharge = sBBCcharge; };
      float      GetBBCcharge(){ return BBCcharge; };

      void       SetBBCchargeN(float sBBCchargeN) { BBCchargeN = sBBCchargeN; };
      float      GetBBCchargeN(){ return BBCchargeN; };

      void       SetBBCchargeS(float sBBCchargeS) { BBCchargeS = sBBCchargeS; };
      float      GetBBCchargeS(){ return BBCchargeS; };

      void       SetBBCtimeN(float sBBCtimeN) { BBCtimeN = sBBCtimeN; };
      float      GetBBCtimeN(){ return BBCtimeN; };

      void       SetBBCtimeS(float sBBCtimeS) { BBCtimeS = sBBCtimeS; };
      float      GetBBCtimeS(){ return BBCtimeS; };

      void       SetVtxZ(float szvertex) { zvertex = szvertex; };
      float      GetVtxZ() { return zvertex; };

      void       SetCentrality(float scentrality) { centrality = scentrality; };
      float      GetCentrality(){ return centrality; };

      void       SetPsi2BBC(float spsi2_BBC) { psi2_BBC = spsi2_BBC; };
      float      GetPsi2BBC(){ return psi2_BBC; };

      void       SetPsi2BBCN(float spsi2_BBCn) { psi2_BBCn = spsi2_BBCn; };
      float      GetPsi2BBCN(){ return psi2_BBCn; };

      void       SetPsi2BBCS(float spsi2_BBCs) { psi2_BBCs = spsi2_BBCs; };
      float      GetPsi2BBCS(){ return psi2_BBCs; };

      void       SetPsi3BBC(float spsi3_BBC) { psi3_BBC = spsi3_BBC; };
      float      GetPsi3BBC(){ return psi3_BBC; };

      void       SetPsi3BBCN(float spsi3_BBCn) { psi3_BBCn = spsi3_BBCn; };
      float      GetPsi3BBCN(){ return psi3_BBCn; };

      void       SetPsi3BBCS(float spsi3_BBCs) { psi3_BBCs = spsi3_BBCs; };
      float      GetPsi3BBCS(){ return psi3_BBCs; };

      void       SetPsi2FVTXA0(float spsi2_FVTXA0) { psi2_FVTXA0 = spsi2_FVTXA0; };
      float      GetPsi2FVTXA0(){ return psi2_FVTXA0; };

      void       SetPsi2FVTXA0N(float spsi2_FVTXA0n) { psi2_FVTXA0n = spsi2_FVTXA0n; };
      float      GetPsi2FVTXA0N(){ return psi2_FVTXA0n; };

      void       SetPsi2FVTXA0S(float spsi2_FVTXA0s) { psi2_FVTXA0s = spsi2_FVTXA0s; };
      float      GetPsi2FVTXA0S(){ return psi2_FVTXA0s; };

      void       SetPsi3FVTXA0(float spsi3_FVTXA0) { psi3_FVTXA0 = spsi3_FVTXA0; };
      float      GetPsi3FVTXA0(){ return psi3_FVTXA0; };

      void       SetPsi3FVTXA0N(float spsi3_FVTXA0n) { psi3_FVTXA0n = spsi3_FVTXA0n; };
      float      GetPsi3FVTXA0N(){ return psi3_FVTXA0n; };

      void       SetPsi3FVTXA0S(float spsi3_FVTXA0s) { psi3_FVTXA0s = spsi3_FVTXA0s; };
      float      GetPsi3FVTXA0S(){ return psi3_FVTXA0s; };

      void       SetPreciseX(float sPreciseX) { preciseX = sPreciseX; };
      float      GetPreciseX(){ return preciseX; };

      void       SetPreciseY(float sPreciseY) { preciseY = sPreciseY; };
      float      GetPreciseY(){ return preciseY; };
 
      void       SetPreciseZ(float sPreciseZ) { preciseZ = sPreciseZ; };
      float      GetPreciseZ(){ return preciseZ; };

      void       AddTrack(MyTrack newTrack) { TrackList.push_back(newTrack); };
      void       AddVTXHit(MyVTXHit newVTXHit) { VTXHitList.push_back(newVTXHit); };
      void       AddCluster(MyCluster newCluster) { ClusterList.push_back(newCluster); };
      void       AddMcTrack(McTrack newMcTrack) { McTrackList.push_back(newMcTrack); };
      void       AddMcCluster(McCluster newMcCluster) { McClusterList.push_back(newMcCluster); };
      void       AddPair(MyPair newPair) { PairList.push_back(newPair); };
      
      Long64_t        GetNtrack() { return TrackList.size(); };
      Long64_t        GetNVTXhit() { return VTXHitList.size(); };
      Long64_t        GetNcluster() { return ClusterList.size(); };
      Long64_t        GetNMcTrack() { return McTrackList.size(); };
      Long64_t        GetNMcCluster() { return McClusterList.size(); };
      Long64_t        GetNpair() { return PairList.size(); };
      
      MyTrack&      GetEntry(int i) { return TrackList[i]; };
      MyVTXHit&      GetVTXHitEntry(int i) { return VTXHitList[i]; };
      MyCluster&    GetClusterEntry(int i) { return ClusterList[i]; };
      McTrack&      GetMcEntry(int i) { return McTrackList[i]; };
      McCluster&      GetMcClusterEntry(int i) { return McClusterList[i]; };
      MyPair&    GetPairEntry(int i) { return PairList[i]; };

      std::vector<DileptonAnalysis::MyTrack> GetTracks() { return TrackList; };
      std::vector<DileptonAnalysis::MyVTXHit> GetVTXHits() { return VTXHitList; };
      std::vector<DileptonAnalysis::MyCluster> GetClusters() { return ClusterList; };
      std::vector<DileptonAnalysis::McTrack> GetMcTracks() { return McTrackList; };
      std::vector<DileptonAnalysis::McCluster> GetMcClusters() { return McClusterList; };
      std::vector<DileptonAnalysis::MyPair> GetPairs() { return PairList; };
      
      ClassDef(MyEvent,1)  //MyEvent structure
   };
}

#endif /*__MYEVENT_OLD_H__*/

