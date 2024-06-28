#ifndef __EMBEDANA_H__
#define __EMBEDANA_H__

#include <SubsysReco.h>
#include <string>
#include <vector>
#include <map>


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <stdlib.h>
#include <cassert>
#include <stdio.h>

#include "TMath.h"
#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "recoConsts.h"


#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"
#include "EventHeader.h"
#include "VtxOut.h"
#include "PHGlobal.h"

#include <DchHitLineTable.hh>
#include <PadCluster.h>
#include <CrkHit.h>
#include <PHTrackOut.h>
#include <PHDchTrackOut.h>
#include <CglTrack.h>
#include "SvxClusterList.h"
#include "SvxGhitClusterList.h"
#include "SvxCluster.h"
#include "SvxGhitCluster.h"
#include <PHPointerList.h>
#include <PHEmbedMcRecoTrack.hh>

#include "TFile.h"
#include "TNtuple.h"

#include "getClass.h"

using namespace std;
using namespace findNode;

#include "MyEvent.h"
#include "TOAD.h"
#include "TMath.h"

class PHCompositeNode;
class TFile;
class TNtuple;

class PHCentralTrack;

struct InData {
    double px, py, pz, vx, vy, vz;
    int id;
};

////////////

class embedana: public SubsysReco {
public:
  embedana(std::string filename="embedana.root", std::string filepath="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
  std::string oscarpath = "/phenix/plhf/mitran/Simul/Dileptons/output_single/single/00002.oscar.particles.dat");
  virtual ~embedana();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);

  inline int Reset(PHCompositeNode *topNode) 		{return 0;}
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

  //void remapCntSvx(PHCentralTrack* cnt, SvxCentralTrackList* svxcnt);

private:
  int EventNumber;

  std::string m_outFileName;
  TFile*      m_outfile;
  int fill_TTree, remove_hadron_hits;
  TNtuple*    m_ntp_embed;
  std::string local_filepath;
  std::string local_oscarpath;
  std::vector<double> vertexes;

  int ReadOrigPartMoms();
  InData InData_read[10000];


  MyDileptonAnalysis::MyEventContainer *event_container;

  //std::vector<SvxCentralTrack*> m_vsvxcnt;
};

#endif

