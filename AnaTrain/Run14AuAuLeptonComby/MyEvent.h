#ifndef __MYEVENT_H__
#define __MYEVENT_H__

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
#include <map>
#include <vector>
#include <iostream>
#include "MyEventConstants.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "stdio.h"
#include "TGraph.h"
#include "TF1.h"

class TDirectory;

namespace MyDileptonAnalysis
{

      class MyTrack : public TObject
      {

      private:
            int trkid;       // DC track index
            int trkinfo;     // arm: 0 east, 1 west
            int trk_quality; // side 0 south, 1 north
                             // sector 0-7 with 0 being lowest west sector

            // basic parameters from tracking
            int q;        // charge of the track
            float pt;     // transverse momentum
            float phi_DC; // phi of track at 220 cm
            float z_DC;   // z poition of track at 220 cm (measured by PC1)
            float phi0;   // phi at interaction point - assuming track is primary
            float the0;   // dito but theta
            float alpha;  // alpha inclination at DC (total field bent = alpha + (phi_DC-phi0))
            int arm;
            int dcside;
            int sect;
            float qpt;
            int n_hits;

            // Updated track parameters

            int q_prime;
            float phi0_prime;
            float the0_prime;
            float pt_prime;
            float alpha_prime;

            // EMCal variables
            int emcid;     // EMC cluster
            int emctower;  // emctower = 1000*iy + 10*iz + sector
                           // EMC sector 0-4 from bottom of arm
                           // iy tower coordinates in rphi
                           // iz tower coordinates in z
            float ecore;   // energy of cluster - assuning its a photon cluster
            float prob;    // similar to chi2 (is possibly redundant)
            float emcdz;   // delta z of cluster and track projection
            float emcdphi; // delta phi of cluster and track projection
            float emctof;  // TOF of EMCal (PbSc only)
            float dep;     // sigmalized E/p

            // TOF variables (should include TOF east and TOF west)
            float tofe; // TOF west??
            float tofdphi;
            float tofdz;
            // pc3 info
            float pc3sdphi;
            float pc3sdz;
            // ERT information
            int isERT; // indicates which ERT trigger fired

            float crkphi; // phi position of ring
            float crkz;   // z position of ring

            // MC variables
            int mcid; // associated MC track
            
            std::vector<float> dangle0,dangle1,dangle2,dangle3;

      public:
            MyTrack()
            {
                  trkid = -99999; // max value of unsigned short
                  trkinfo = -99999;
                  trk_quality = -99999;
                  arm = -99999;
                  dcside = -99999;
                  sect = -99999;
                  emcid = -99999;
                  emctower = -99999;
                  isERT = -99999;
                  mcid = -99999;
                  q = -99999;
                  qpt = -99999;
                  q_prime = -99999;
                  n_hits = -99999;
                  pt = -99999;
                  pt_prime = -99999;
                  phi_DC = -99999;
                  z_DC = -99999;
                  phi0 = -99999;
                  phi0_prime = -99999;
                  the0 = -99999;
                  the0_prime = -99999;
                  alpha = -99999;
                  alpha_prime = -99999;
                  ecore = -99999;
                  dep = -99999;
                  prob = -99999;
                  emcdz = -99999;
                  emcdphi = -99999;
                  emctof = -99999;
                  tofe = -99999;
                  tofdphi = -99990;
                  tofdz = -99999;
                  pc3sdphi = -99999;
                  pc3sdz = -99999;
                  crkphi = -99999;
                  crkz = -99999;
            };

            virtual ~MyTrack(){};

            int GetTrkId() const { return trkid; };
            int GetTrkInfo() const { return trkinfo; };
            int GetTrkQuality() const { return trk_quality; };
            int GetArm() const { return arm; };
            int GetDCside() const { return dcside; };
            int GetSect() const { return sect; };
            float GetPx() const
            {
                  if (pt == -99999)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime * cos(phi0);
                  }
            };
            float GetPy() const
            {
                  if (pt == -99999)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime * sin(phi0);
                  }
            };
            float GetPz() const
            {
                  if (pt == -99999)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime / tan(the0);
                  }
            };
            float GetPtot() const
            {
                  if (pt == -99999)
                  {
                        return pt;
                  }
                  else
                  {
                        return pt_prime * sqrt(1.+1./tan(the0)/tan(the0));
                  }
            };
            float GetPt() const { return pt; };
            float GetPtPrime() const { return pt_prime; };
            float GetPhiDC() const { return phi_DC; };
            float GetZDC() const { return z_DC; };
            float GetPhi0() const { return phi0; };
            int GetNHits() const { return n_hits; };
            float GetPhi0Prime() const { return phi0_prime; };
            float GetThe0() const { return the0; };
            float GetThe0Prime() const { return the0_prime; };
            float GetAlpha() const { return alpha; };
            float GetAlphaPrime() const { return alpha_prime; };
            int GetCharge() const { return q; };
            int GetChargePrime() const { return q_prime; };
            int GetEmcId() const { return emcid; };
            float GetEcore() const { return ecore; };
            float GetDep() const { return dep; };
            float GetProb() const { return prob; };
            float GetEmcdz() const { return emcdz; };
            float GetEmcdphi() const { return emcdphi; };
            int GetEsect() const { return int(emctower - 10 * ((emctower - 1000 * (emctower / 1000)) / 10) - 1000 * (emctower / 1000)); };
            int GetYsect() const { return int(emctower / 1000); };
            int GetZsect() const { return int((emctower - 1000 * (emctower / 1000)) / 10); };
            float GetEmcTOF() const { return emctof; };
            float GetTOFE() const { return tofe; };
            float GetTOFDPHI() const { return tofdphi; };
            float GetTOFDZ() const { return tofdz; };
            float GetPC3SDPHI() const { return pc3sdphi; };
            float GetPC3SDZ() const { return pc3sdz; };
            int GetisERT() const { return isERT; };
            int GetMcId() const { return mcid; };
            float GetCrkphi() const { return crkphi; };
            float GetCrkz() const { return crkz; };

            void SetTrkId(int strkid) { trkid = strkid; };
            void SetTrkQuality(int strk_quality) { trk_quality = strk_quality; };
            void SetArm(int sarm) { arm = sarm; };
            void SetTrkInfo(int sarm, int sdcside, int ssect) { trkinfo = short(ssect + 10 * sdcside + 100 * sarm); };
            void SetQpt(float spx, float spy, float scharge) { qpt = scharge * sqrt(spx * spx + spy * spy); };
            void SetDCSide(int sdcside) { dcside = sdcside; };
            void SetSect(int ssect) { sect = ssect; };
            void SetPt(float spt) { pt = spt; pt_prime = spt; };
            void SetPtPrime(float spt_prime) { pt_prime = spt_prime; };
            void SetQ(int sq) { q = sq; };
            void SetQPrime(int sq_prime) { q_prime = sq_prime; };
            void SetPhiDC(float sphi_DC) { phi_DC = sphi_DC; };
            void SetZDC(float sz_DC) { z_DC = sz_DC; };
            void SetPhi0(float sphi0) { phi0 = sphi0; };
            void SetPhi0Prime(float sphi0_prime) { phi0_prime = sphi0_prime; };
            void SetThe0(float sthe0) { the0 = sthe0; };
            void SetThe0Prime(float sthe0_prime) { the0_prime = sthe0_prime; };
            void SetAlpha(float salpha) { alpha = salpha; };
            void SetAlphaPrime(float salpha_prime) { alpha_prime = salpha_prime; };
            void SetEmcId(int sid) { emcid = sid; };
            void SetEcore(float secore) { ecore = secore; };
            void SetDep(float sdep) { dep = sdep; };
            void SetProb(float sprob) { prob = sprob; };
            void SetEmcdz(float semcdz) { emcdz = semcdz; };
            void SetEmcdphi(float semcdphi) { emcdphi = semcdphi; }
            void SetEmcTower(float sesect, float sysect, float szsect) { emctower = 1000 * sysect + 10 * szsect + sesect; };
            void SetEmcTOF(float semctof) { emctof = semctof; };
            void SetTOFE(float stofe) { tofe = stofe; };
            void SetTOFDPHI(float stofdphi) { tofdphi = stofdphi; };
            void SetTOFDZ(float stofdz) { tofdz = stofdz; };
            void SetPC3SDPHI(float spc3sdphi) { pc3sdphi = spc3sdphi; };
            void SetPC3SDZ(float spc3sdz) { pc3sdz = spc3sdz; };
            void SetisERT(int sisERT) { isERT = sisERT; };
            void SetMcId(int smcid) { mcid = smcid; };
            void SetCrkphi(float scrkphi) { crkphi = scrkphi; };
            void SetCrkz(float scrkz) { crkz = scrkz; };

            void SetNHits(int sn_hits) { n_hits = sn_hits; };

            void SetPrimes(const float bbcz = 0, const float svxz = 0, const int rg_beamoffset = 0, const int rungroup = 0);

            float get_mean_phi_data(int rungroup, int centr_bin, int layer) { return mean_phi_pars[rungroup][centr_bin][layer]; };
            float get_mean_theta_data(int rungroup, int centr_bin, int layer) { return mean_theta_pars[rungroup][centr_bin][layer]; };
            float get_sigma_phi_data(int rungroup, int centr_bin, int layer)
            {
                  return sigma_phi_pars[rungroup][centr_bin][layer][0] + sigma_phi_pars[rungroup][centr_bin][layer][1] * exp(sigma_phi_pars[rungroup][centr_bin][layer][2] * pt_prime);
            };
            float get_sigma_theta_data(int rungroup, int centr_bin, int layer)
            {
                  return sigma_theta_pars[rungroup][centr_bin][layer][0] + sigma_theta_pars[rungroup][centr_bin][layer][1] * exp(sigma_theta_pars[rungroup][centr_bin][layer][2] * pt_prime);
            };

            float get_dynamic_mean_phi_data(int rungroup, int ilay, float phi_prev) 
            { 
                  const int arg0 = 2*ilay + (1-q_prime)/2;
                  return phi_mean_phi_params[rungroup][arg0][0]+phi_mean_phi_params[rungroup][arg0][1]*phi_prev; 
            };
            float get_dynamic_mean_theta_data(int rungroup, int ilay, float phi_prev)
            { 
                  const int arg0 = 2*ilay + (1-q_prime)/2;
                  return the_mean_the_params[rungroup][arg0][0]+the_mean_the_params[rungroup][arg0][1]*phi_prev; 
            };
            float get_dynamic_sigma_phi_data(int rungroup, int ilay, float phi_prev)
            {
                  const int arg0 = 2*ilay + (1-q_prime)/2;
                  const float sigma_pt = phi_sigma_pt_params[rungroup][arg0][0] + phi_sigma_pt_params[rungroup][arg0][1] * exp(phi_sigma_pt_params[rungroup][arg0][2] * pt_prime);
                  return sigma_pt*(phi_sigma_phi_params[rungroup][arg0][0]+phi_sigma_phi_params[rungroup][arg0][1]* phi_prev+phi_sigma_phi_params[rungroup][arg0][2]* phi_prev* phi_prev);
            };
            float get_dynamic_sigma_theta_data(int rungroup, int ilay, float phi_prev)
            {
                  const int arg0 = 2*ilay + (1-q_prime)/2;
                  return the_sigma_pt_params[rungroup][arg0][0] + the_sigma_pt_params[rungroup][arg0][1] * exp(the_sigma_pt_params[rungroup][arg0][2] * pt_prime);
            };
            float get_dynamic_smean_phi_data(int rungroup, int ilay, float phi_prev)
            { 
                  const int arg0 = 2*ilay + (1-q_prime)/2;
                  return phi_sMean_pt_params[rungroup][arg0][0] + phi_sMean_pt_params[rungroup][arg0][1] * exp(phi_sMean_pt_params[rungroup][arg0][2] * pt_prime); 
            };


            void SetdPhidThe(int layer, float val1, float val2, float val3, float val4, float val5, float val6) 
            {
                  if(layer==0) {dangle0.push_back(val1);dangle0.push_back(val2);dangle0.push_back(val3);dangle0.push_back(val4);dangle0.push_back(val5);dangle0.push_back(val6);}
                  if(layer==1) {dangle1.push_back(val1);dangle1.push_back(val2);dangle1.push_back(val3);dangle1.push_back(val4);dangle1.push_back(val5);dangle1.push_back(val6);}
                  if(layer==2) {dangle2.push_back(val1);dangle2.push_back(val2);dangle2.push_back(val3);dangle2.push_back(val4);dangle2.push_back(val5);dangle2.push_back(val6);}
                  if(layer==3) {dangle3.push_back(val1);dangle3.push_back(val2);dangle3.push_back(val3);dangle3.push_back(val4);dangle3.push_back(val5);dangle3.push_back(val6);}
            };

            float GetdPhi(int layer, unsigned int iter) 
            {
                  if     (layer==0&& 6*iter<dangle0.size()) return dangle0[iter*6];
                  else if(layer==1&& 6*iter<dangle1.size()) return dangle1[iter*6];
                  else if(layer==2&& 6*iter<dangle2.size()) return dangle2[iter*6];
                  else if(layer==3&& 6*iter<dangle3.size()) return dangle3[iter*6];
                  std::cout<<"GetdPhi "<<layer<<" "<<6*iter+0<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetdThe(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+1<dangle0.size()) return dangle0[iter*6+1];
                  else if(layer==1 && 6*iter+1<dangle1.size()) return dangle1[iter*6+1];
                  else if(layer==2 && 6*iter+1<dangle2.size()) return dangle2[iter*6+1];
                  else if(layer==3 && 6*iter+1<dangle3.size()) return dangle3[iter*6+1];
                  std::cout<<"GetdThe "<<layer<<" "<<6*iter+1<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetsdPhi(int layer, unsigned int iter) 
            {
                  if     (layer==0&& 6*iter+2<dangle0.size()) return dangle0[iter*6+2];
                  else if(layer==1&& 6*iter+2<dangle1.size()) return dangle1[iter*6+2];
                  else if(layer==2&& 6*iter+2<dangle2.size()) return dangle2[iter*6+2];
                  else if(layer==3&& 6*iter+2<dangle3.size()) return dangle3[iter*6+2];
                  std::cout<<"GetsdPhi "<<layer<<" "<<6*iter+2<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetsdThe(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+3<dangle0.size()) return dangle0[iter*6+3];
                  else if(layer==1 && 6*iter+3<dangle1.size()) return dangle1[iter*6+3];
                  else if(layer==2 && 6*iter+3<dangle2.size()) return dangle2[iter*6+3];
                  else if(layer==3 && 6*iter+3<dangle3.size()) return dangle3[iter*6+3];
                  std::cout<<"GetsdThe "<<layer<<" "<<6*iter+3<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            float GetDist(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+4<dangle0.size()) return dangle0[iter*6+4];
                  else if(layer==1 && 6*iter+4<dangle1.size()) return dangle1[iter*6+4];
                  else if(layer==2 && 6*iter+4<dangle2.size()) return dangle2[iter*6+4];
                  else if(layer==3 && 6*iter+4<dangle3.size()) return dangle3[iter*6+4];
                  std::cout<<"GetDist "<<layer<<" "<<6*iter+4<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };
            int GetHits(int layer, unsigned int iter) 
            {
                  if     (layer==0 && 6*iter+5<dangle0.size()) return (int) dangle0[iter*6+5];
                  else if(layer==1 && 6*iter+5<dangle1.size()) return (int) dangle1[iter*6+5];
                  else if(layer==2 && 6*iter+5<dangle2.size()) return (int) dangle2[iter*6+5];
                  else if(layer==3 && 6*iter+5<dangle3.size()) return (int) dangle3[iter*6+5];
                  std::cout<<"GetHits "<<layer<<" "<<6*iter+5<<" "<<dangle0.size()<<std::endl;
                  return -999;
            };

            /// below defined virtual functions for ide to suggest based on inherentence in electron and hadron classes

            virtual float GetChi2() const { return 0; };
            virtual int GetN0() const { return 0; };
            virtual int GetNpe0() const { return 0; };
            virtual int GetDisp() const { return 0; };
            virtual float GetEmcdz_e() const { return 0; };
            virtual float GetEmcdphi_e() const { return 0; };

            virtual void SetChi2(float schi2) { };
            virtual void SetN0(int sn0) { };
            virtual void SetNPE0(int snpe0) { };
            virtual void SetDISP(int sdisp) { };
            virtual void SetEmcdz_e(float semcdz) { };
            virtual void SetEmcdphi_e(float semcdphi) { }

            virtual void SetHitIndex(int iindex, int iilayer) { };
            virtual void AddHitCounter(int iilayer) {  };
            virtual void RemoveHitCounter(int iilayer) {  };
            virtual void SetMinDist(float imindist, int iilayer) { };
            virtual void SetGhost(int value) { };

            virtual int GetHitIndex(int iilayer) { return 0; };
            virtual int GetHitCounter(int iilayer) { return 0; };
            virtual float GetMinDist(int iilayer) { return 0; };
            virtual int GetGhost() { return 0; } 

            ClassDef(MyTrack, 1)
      };

      class MyElectron : public MyTrack
      {

      private:
            int n0;       // number of photo tubes in mask
            float disp;
            float npe0;
            float chi2;

            float emcdz_e;   // delta z of cluster and track projection
            float emcdphi_e; // delta phi of cluster and track projection

            int hit_index[4];
            int hit_counter[4];
            int isGhost;
            float min_dist[4];
            float min_dphi[4];
            float min_dthe[4];
            float min_sdphi[4];
            float min_sdthe[4];
            float DCA, DCA2;
            float sDCA;
            float DCA_X,DCA_X2;
            float DCA_Y,DCA_Y2;
            float true_pt;

      public:
            MyElectron() : MyTrack()
            {
                  n0 = -99999; // max value for byte
                  disp = -99999;
                  chi2 = -99999;
                  npe0 = -99999;
                  emcdz_e = -99999;
                  emcdphi_e = -99999;
                  DCA = -99999;
                  sDCA = -99999;
                  DCA_X = -99999;
                  DCA_Y = -99999;
                  DCA2 = -99999;
                  DCA_X2 = -99999;
                  DCA_Y2 = -99999;
                  true_pt = -99999;

                  for (int iteri = 0; iteri < 4; iteri++)
                  {
                        hit_index[iteri] = -99999;
                        hit_counter[iteri] = 0;
                        min_dist[iteri] = 100.;
                        min_dphi[iteri] = 100;
                        min_dthe[iteri] = 100;
                        min_sdphi[iteri] = 100;
                        min_sdthe[iteri] = 100;
                  }
                  isGhost = 0;
            };

            virtual ~MyElectron(){};

            float GetChi2() const { return chi2; };
            int GetN0() const { return n0; };
            int GetNpe0() const { return npe0; };
            int GetDisp() const { return disp; };
            float GetEmcdz_e() const { return emcdz_e; };
            float GetEmcdphi_e() const { return emcdphi_e; };

            void SetChi2(float schi2) { chi2 = schi2; };
            void SetN0(int sn0) { n0 = sn0; };
            void SetNPE0(int snpe0) { npe0 = snpe0; };
            void SetDISP(int sdisp) { disp = sdisp; };
            void SetEmcdz_e(float semcdz) { emcdz_e = semcdz; };
            void SetEmcdphi_e(float semcdphi) { emcdphi_e = semcdphi; }

            void SetHitIndex(int iindex, int iilayer) { hit_index[iilayer] = iindex; };
            void AddHitCounter(int iilayer) { hit_counter[iilayer]++; };
            void SetHitCounter(int iilayer, int value) { hit_counter[iilayer] = value; };
            void RemoveHitCounter(int iilayer) { hit_counter[iilayer]--; };
            void ZeroHitCounters() { for (int i = 0; i < 4; i++) hit_counter[i]=0; };
            void SetMinDist(float imindist, int iilayer) { min_dist[iilayer] = imindist; };
            void SetMinDphi(float imindist, int iilayer) { min_dphi[iilayer] = imindist; };
            void SetMinDthe(float imindist, int iilayer) { min_dthe[iilayer] = imindist; };
            void SetGhost(int value) { isGhost=value; };
            void SetDCA(float value) {DCA = value;};
            void SetsDCA(float value) {sDCA = value;};
            void SetDCAX(float value) {DCA_X = value;};
            void SetDCAY(float value) {DCA_Y = value;};
            void SetDCA2(float value) {DCA2 = value;};
            void SetDCAX2(float value) {DCA_X2 = value;};
            void SetDCAY2(float value) {DCA_Y2 = value;};
            void SetReconPT(float value) {true_pt = value;};
            void SetMinsDphi(float imindist, int iilayer) { min_sdphi[iilayer] = imindist; };
            void SetMinsDthe(float imindist, int iilayer) { min_sdthe[iilayer] = imindist; };

            int GetHitIndex(int iilayer) { return hit_index[iilayer]; };
            int GetHitCounter(int iilayer) { return hit_counter[iilayer]; };
            float GetMinDist(int iilayer) { return min_dist[iilayer]; };
            float GetMinDphi(int iilayer) { return min_dphi[iilayer]; };
            float GetMinDthe(int iilayer) { return min_dthe[iilayer]; };
            int GetGhost() { return isGhost; }
            float GetDCA() {return DCA;};
            float GetsDCA() {return sDCA;};
            float GetDCAX() {return DCA_X;};
            float GetDCAY() {return DCA_Y;};
            float GetDCA2() {return DCA2;};
            float GetDCAX2() {return DCA_X2;};
            float GetDCAY2() {return DCA_Y2;};
            float GetReconPT() {return true_pt;};
            float GetMinsDphi(int iilayer) { return min_sdphi[iilayer]; };
            float GetMinsDthe(int iilayer) { return min_sdthe[iilayer]; };

            ClassDef(MyElectron, 1)
      };

      class MyHadron : public MyTrack
      {
      private:
            int hit_index[4];
            int hit_counter[4];
            int isGhost;
            float min_dist[4];

      public:
            MyHadron() : MyTrack()
            {
                  for (int iteri = 0; iteri < 4; iteri++)
                  {
                        hit_index[iteri] = -99999;
                        hit_counter[iteri] = 0;
                        min_dist[iteri] = 100.;
                  }
                  isGhost = 0;
            };

            virtual ~MyHadron(){};

            void SetHitIndex(int iindex, int iilayer) { hit_index[iilayer] = iindex; };
            void AddHitCounter(int iilayer) { hit_counter[iilayer]++; };
            void RemoveHitCounter(int iilayer) { hit_counter[iilayer]--; };
            void SetMinDist(float imindist, int iilayer) { min_dist[iilayer] = imindist; };
            void SetGhost(int value) { isGhost=value; };

            int GetHitIndex(int iilayer) { return hit_index[iilayer]; };
            int GetHitCounter(int iilayer) { return hit_counter[iilayer]; };
            float GetMinDist(int iilayer) { return min_dist[iilayer]; };
            int GetGhost() { return isGhost; }

            ClassDef(MyHadron, 1)
      };

      class MyVTXHit : public TObject
      {

      private:
            int clustid;
            int layer;
            int ladder;
            int sensor;

            float xhit; // xyz position of cluster [cm]
            float yhit;
            float zhit;

            float rhit;
            float phihit;
            float thetahit;

            short ilayer;

            float adc[4];

            std::map<int, float> associated_tracks;
            //std::vector<double> dangle1;

      public:
            MyVTXHit()
            {
                  clustid = -99999;
                  layer = -99999;
                  ladder = -8888;
                  sensor = -99999;

                  xhit = -99999;
                  yhit = -99999;
                  zhit = -99999;

                  rhit = -8888;
                  phihit = -8888;
                  thetahit = -8888;

                  ilayer = -8888;
                  adc[0]=-9999;
                  adc[1]=-9999;
                  adc[2]=-9999;
                  adc[3]=-9999;
            };
            virtual ~MyVTXHit(){};

            int GetClustId() const { return clustid; };
            int GetLayer() const { return layer; };
            int GetLadder() const { return ladder; };
            int GetSensor() const { return sensor; };

            float GetXHit() const { return xhit; };
            float GetYHit() const { return yhit; };
            float GetZHit() const { return zhit; };

            float GetRHit() const { return rhit; };
            float GetPhiHit() const { return phihit; };
            float GetTheHit() const { return thetahit; };

            short GetiLayer() const { return ilayer; };

            float GetPhi() const
            {
                  float phi = atan2(yhit, xhit);
                  if (phi < -TMath::ACos(-1) / 2)
                        phi += 2 * TMath::ACos(-1);
                  return phi;
            };

            void SetClustId(int sclustid) { clustid = sclustid; };
            void SetLayer(int slayer) { layer = slayer; };
            void SetLadder(int sladder) { ladder = sladder; };
            void SetSensor(int ssensor) { sensor = ssensor; };

            void SetXHit(float sx) { xhit = sx; };
            void SetYHit(float sy) { yhit = sy; };
            void SetZHit(float sz) { zhit = sz; };

            void SetPolars(const float xvtx = 0, const float yvtx = 0, const float zvtx = 0);
            void SetiLayerFromR();

            void AddAssociatedTrack(int trackID, float dist) { associated_tracks[trackID] = dist; };

            int N_AssociatedTracks() { return associated_tracks.size(); };

            int GetAssociatedTrack(int itrack)
            {
                  std::map<int, float>::const_iterator iterator = associated_tracks.begin();
                  for (int iass = 0; iass < itrack; iass++)
                        iterator++;
                  return iterator->first;
            };
            float GetAssociatedTrackDistance(int itrackid) { return associated_tracks[itrackid]; };
            void DeleteAssociatedTrack(int itrackid) { associated_tracks.erase(itrackid); };

            void SetAdc(int layer, float adc1, float adc2){if(layer>1) {adc[layer*2-4]=adc1;adc[layer*2-3]=adc2;}};
            float GetAdc1(int layer){if (layer>1) return adc[layer*2-4]; else return -9999;}
            float GetAdc2(int layer){if (layer>1) return adc[layer*2-3]; else return -9999;}

            ClassDef(MyVTXHit, 1)
      };

      class MyEvent : public TObject
      {

      private:
            int evtno;
            float BBCcharge;
            float BBCchargeN;
            float BBCchargeS;
            float BBCtimeN;
            float BBCtimeS;
            float zvertex;
            float centrality;
            float preciseX;
            float preciseY;
            float preciseZ;
            float run_number;

            std::vector<MyDileptonAnalysis::MyElectron> TrackList;
            std::vector<MyDileptonAnalysis::MyHadron> HadronList;
            std::vector<MyDileptonAnalysis::MyVTXHit> VTXHitList;
            std::vector<MyDileptonAnalysis::MyElectron> ElecCandList;

      public:
            MyEvent()
            {
                  evtno = -999;
                  BBCcharge = -999;
                  BBCchargeN = -999;
                  BBCchargeS = -999;
                  BBCtimeN = -999;
                  BBCtimeS = -999;
                  zvertex = -999;
                  centrality = -999;

                  preciseX = -999;
                  preciseY = -999;
                  preciseZ = -999;
                  run_number = -999;
            };

            virtual ~MyEvent(){};

            void ClearEvent();

            void SetEvtNo(int sevtno) { evtno = sevtno; };
            int GetEvtNo() { return evtno; };

            void SetRunNumber(int srun_number) { run_number = srun_number; };
            int GetRunNumber() { return run_number; };

            void SetBBCcharge(float sBBCcharge) { BBCcharge = sBBCcharge; };
            float GetBBCcharge() { return BBCcharge; };

            void SetBBCchargeN(float sBBCchargeN) { BBCchargeN = sBBCchargeN; };
            float GetBBCchargeN() { return BBCchargeN; };

            void SetBBCchargeS(float sBBCchargeS) { BBCchargeS = sBBCchargeS; };
            float GetBBCchargeS() { return BBCchargeS; };

            void SetBBCtimeN(float sBBCtimeN) { BBCtimeN = sBBCtimeN; };
            float GetBBCtimeN() { return BBCtimeN; };

            void SetBBCtimeS(float sBBCtimeS) { BBCtimeS = sBBCtimeS; };
            float GetBBCtimeS() { return BBCtimeS; };

            void SetVtxZ(float szvertex) { zvertex = szvertex; };
            float GetVtxZ() { return zvertex; };

            void SetCentrality(float scentrality) { centrality = scentrality; };
            float GetCentrality() { return centrality; };

            void SetPreciseX(float sPreciseX) { preciseX = sPreciseX; };
            float GetPreciseX() { return preciseX; };

            void SetPreciseY(float sPreciseY) { preciseY = sPreciseY; };
            float GetPreciseY() { return preciseY; };

            void SetPreciseZ(float sPreciseZ) { preciseZ = sPreciseZ; };
            float GetPreciseZ() { return preciseZ; };

            int GetRunGroup(int in_run_number);

            void AddTrack(const MyElectron *newTrack) { TrackList.push_back(*newTrack); };
            void AddHadron(const MyHadron *newTrack) { HadronList.push_back(*newTrack); };
            void AddVTXHit(const MyVTXHit *newVTXHit) { VTXHitList.push_back(*newVTXHit); };
            void AddElecCand(const MyElectron *newTrack) { ElecCandList.push_back(*newTrack); };

            Long64_t GetNtrack() { return TrackList.size(); };
            Long64_t GetNhadron() { return HadronList.size(); };
            Long64_t GetNVTXhit() { return VTXHitList.size(); };
            Long64_t GetNeleccand() { return ElecCandList.size(); };

            MyElectron *GetEntry(unsigned int i) const { return const_cast<MyElectron *>(&(TrackList[i])); };
            MyHadron *GetHadronEntry(unsigned int i) const { return const_cast<MyHadron *>(&(HadronList[i])); };
            MyVTXHit *GetVTXHitEntry(unsigned int i) const { return const_cast<MyVTXHit *>(&(VTXHitList[i])); };
            MyElectron *GetElecCand(unsigned int i) const { return const_cast<MyElectron *>(&(ElecCandList[i])); };

            void RemoveTrackEntry(const unsigned int i){ TrackList.erase(TrackList.begin() + i); };

            std::vector<MyDileptonAnalysis::MyElectron> GetTracks() { return TrackList; };
            std::vector<MyDileptonAnalysis::MyHadron> GetHadrons() { return HadronList; };
            std::vector<MyDileptonAnalysis::MyVTXHit> GetVTXHits() { return VTXHitList; };

            void SetDCA(const unsigned int itrk = 0, const int layer2 = 1);
            void SetDCA2(const unsigned int itrk = 0, const int layer3 = 2);
            void ReshuffleElectrons();

            ClassDef(MyEvent, 1) // MyEvent structure
      };

      class MyEventContainer : public TObject
      {

      private:
            MyEvent *event;
            std::vector<MyDileptonAnalysis::MyEvent> EventList;
            TFile *infile, *outfile;
            TH2D *hist_br, *hist_bz;
            TTree *tree;
            TH3D *dphi_hist[N_centr], *dthe_hist[N_centr], *sdphi_hist[N_centr], *sdthe_hist[N_centr];
            TH3D *chi2_ndf[N_centr];
            TH3D *dphi_hist_el[N_centr], *dthe_hist_el[N_centr], *sdphi_hist_el[N_centr], *sdthe_hist_el[N_centr];
            TH3D *dphi_hist_el_dynamic[N_dynamic], *dthe_hist_el_dynamic[N_dynamic], *sdphi_hist_el_dynamic[N_dynamic], *sdthe_hist_el_dynamic[N_dynamic];
            TH3D *d_dphi_hist[N_centr], *d_dthe_hist[N_centr], *DCA_hist[N_centr];
            TH3D *sd_dphi_hist[N_centr], *sd_dthe_hist[N_centr], *sDCA_hist[N_centr];
            TH3D *temc, *ttof, *n0_hist, *ep_hist, *prob_hist, *disp_hist, *chi2npe0_hist;
            TH3D *DCA_2D_hist[N_centr], *sDCA_2D_hist[N_centr];
            TH3D *emc_dphi_el, *emc_dz_el, *n0_hist_el, *ep_hist_el, *prob_hist_el, *disp_hist_el, *chi2npe0_hist_el;
            TH3D *el_had_dphi, *el_had_dz, *el_had_dr;
            TH3D *DCPT_ReconPT, *sDCPT_ReconPT, *DCA2_hist[N_centr], *sDCA2_hist[N_centr], *DCA2_2D_hist[N_centr], *sDCA2_2D_hist[N_centr], *DCA12_hist[N_centr], *charge_hist;
            TH3D *veto_hist[N_centr], *veto_hist_the[N_centr], *sveto_hist[N_centr];
            TH2D *couter_veto_hist;
            TH3D *adc_hist;
            int is_fill_hsits, is_fill_hadron_hsits, is_fill_tree, is_fill_dphi_hist, is_fill_DCA_hist, is_fill_track_QA, 
            is_fill_reveal, is_fill_DCA2_hist, is_check_veto;
           
      public:
            MyEventContainer()
            {
                  event = nullptr;
                  infile = nullptr;
                  outfile = nullptr;
                  hist_br = nullptr;
                  hist_bz = nullptr;
                  tree = nullptr;
                  temc = nullptr; ttof = nullptr; n0_hist = nullptr; ep_hist=nullptr; prob_hist=nullptr; disp_hist=nullptr; chi2npe0_hist=nullptr;
                  emc_dphi_el = nullptr; emc_dz_el = nullptr; n0_hist_el = nullptr; ep_hist_el=nullptr; prob_hist_el=nullptr; disp_hist_el=nullptr; chi2npe0_hist_el=nullptr;
                  el_had_dphi = nullptr, el_had_dz = nullptr, el_had_dr = nullptr, DCPT_ReconPT = nullptr, sDCPT_ReconPT = nullptr, charge_hist = nullptr;
                  couter_veto_hist = nullptr;
                  adc_hist = nullptr;
                  for (int i = 0; i < N_dynamic; i++)
                  {
                        dphi_hist_el_dynamic[i] = nullptr;
                        dthe_hist_el_dynamic[i] = nullptr;
                        sdphi_hist_el_dynamic[i] = nullptr;
                        sdthe_hist_el_dynamic[i] = nullptr;
                  }
                  for (int i = 0; i < N_centr; i++)
                  {
                        dphi_hist[i] = nullptr;
                        dthe_hist[i] = nullptr;
                        sdphi_hist[i] = nullptr;
                        sdthe_hist[i] = nullptr;
                        dphi_hist_el[i] = nullptr;
                        dthe_hist_el[i] = nullptr;
                        sdphi_hist_el[i] = nullptr;
                        sdthe_hist_el[i] = nullptr;
                        chi2_ndf[i] = nullptr;
                        d_dphi_hist[i] = nullptr;
                        d_dthe_hist[i] = nullptr;
                        DCA_hist[i] = nullptr;
                        sd_dphi_hist[i] = nullptr;
                        sd_dthe_hist[i] = nullptr;
                        sDCA_hist[i] = nullptr;
                        DCA_2D_hist[i] = nullptr;
                        sDCA_2D_hist[i] = nullptr;
                        DCA2_hist[i] = nullptr;
                        sDCA2_hist[i] = nullptr;
                        DCA2_2D_hist[i] = nullptr;
                        sDCA2_2D_hist[i] = nullptr;
                        DCA12_hist[i]  = nullptr;
                        veto_hist[i] = nullptr, veto_hist_the[i] = nullptr, sveto_hist[i] = nullptr;
                  }
                  is_fill_hsits = 0;
                  is_fill_hadron_hsits = 0;
                  is_fill_tree = 0;
                  is_fill_dphi_hist = 0;
                  is_fill_DCA_hist = 0;
                  is_fill_track_QA = 0;
                  is_fill_reveal = 0;
                  is_fill_DCA2_hist = 0;
                  is_check_veto = 0;
            };

            virtual ~MyEventContainer()
            {
                  std::cout << "Starting to Delete Event in my Container" << std::endl;
                  if(is_fill_tree||is_fill_hadron_hsits||is_fill_hsits||is_fill_dphi_hist||is_fill_DCA_hist||is_fill_track_QA
                  ||is_fill_reveal||is_fill_DCA2_hist||is_check_veto) {delete outfile;}
                  delete event;
                  std::cout << "Event in the container was deleted" << std::endl;
            };

            void InitEvent() { event = new MyEvent; };
            void ClearEvent() { event->ClearEvent(); };
            void GetHistsFromFile(const std::string loc);
            void CreateOutFileAndInitHists(std::string outfilename, const int fill_el = 0, const int fill_had = 0, const int fill_tree = 0, const int fill_dphi = 0, 
            const int fill_DCA = 0, const int fill_track_QA = 0, const int fill_reveal = 0, const int fill_true_DCA= 0, const int check_veto= 0);
            void ResetTree() {tree->Reset();};
            void FillTree() {tree->Fill();};
            void WriteOutFile();
            void Associate_Hits_to_Leptons(bool test = false);
            void Associate_Hits_to_Hadrons();

            void Reveal_Hadron();

            void FillDphiHists();
            void FillTrueDCA();

            MyEvent *GetEvent() { return event; };
            void SetEvent(MyDileptonAnalysis::MyEvent *ev) { event = ev; };

            void CheckVeto();

            ClassDef(MyEventContainer, 1) // MyEvent structure
      };
}

#endif /*__MYEVENT_H__*/