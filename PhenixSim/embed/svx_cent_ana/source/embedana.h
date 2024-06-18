#ifndef __EMBEDANA_H__
#define __EMBEDANA_H__

#include <SubsysReco.h>
#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class TFile;
class TNtuple;
//class TTree;
//class TH1F;
//class TH2F;
//class TH3F;
//class TH3D;

//class SvxRawhitClusterList;
//class SvxRawhitList;
//class SvxClusterList;
//class SvxCluster;

class PHCentralTrack;
class SvxCentralTrackList;
class SvxCentralTrack;

//class BbcOut;
//class VtxOut;
//class EventHeader;
//class SvxClusterContainer;
//class SvxClusterInfo;
//class RpSumXYObject;
//class McEvalSingleList;
//class SvxGhitList;
//class SvxGhit;
//class SvxGhitRawhitList;
//class SvxHitMap;
//class fkinWrapper;
//class svxAddress;
//class svxDetectorGeo;


////////////

class embedana: public SubsysReco {
public:
  embedana(std::string filename="embedana.root");
  virtual ~embedana();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
  inline int ResetEvent(PHCompositeNode *topNode) 	{return 0;}
  int InitRun(PHCompositeNode *topNode);
  inline void Print(const char *what) const 		{return;}

private:
  void calcDCA_BCbyCircleProjection
            (
              float pt, float phi, float the,  // pt in xy-plane, phi, theta at inner most layer
              int charge,                      // charge of the track
              float hx, float hy, float hz,    // hit position at inner most layer
              float vx, float vy,              // beam center
              float fieldPolarity,             //
              float* dx, float* dy, float* dz, // dca position
              float* d2dca_bc                  // return
            );

  void remapCntSvx(PHCentralTrack* cnt, SvxCentralTrackList* svxcnt);

private:
  int EventNumber;

  std::string m_outFileName;
  TFile*      m_outfile;

  TNtuple*    m_ntp_embed;

  std::vector<SvxCentralTrack*> m_vsvxcnt;
};

#endif

