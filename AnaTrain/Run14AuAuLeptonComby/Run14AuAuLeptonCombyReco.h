#ifndef __RUN14AUAULEPTONCOMBYRECO_H__
#define __RUN14AUAULEPTONCOMBYRECO_H__

#include "SubsysReco.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "phool.h"
#include "PHCentralTrack.h"
#include "PHGlobal.h"
#include "TrigLvl1.h"
#include "RunHeader.h"
#include "TriggerHelper.h"
#include "emcClusterContainer.h"
#include "emcClusterContent.h"
#include "ErtOutv1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TFile.h"
#include "getClass.h"
#include "TOAD.h"
#include "UltraLight/UltraLight.h"
#include "UltraLight/UltraLightTrack.h"
#include "Run14AuAuLeptonCombyEnum.h"
#include "Run14AuAuLeptonCombyConstants.h"
#include "TTree.h"
#include "Declarations.h"

#include "recoConsts.h"

#include "MyEvent.h"
#include "SvxClusterList.h"
#include "SvxCluster.h"
#include <TMath.h>
#include "VtxOut.h"
#include "TVector3.h"
#include "PHPoint.h"

#include <iostream>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <cmath>


class Run14AuAuLeptonCombyReco : public SubsysReco
{
public:
    Run14AuAuLeptonCombyReco(const char *outfile = "test.root");
    virtual ~Run14AuAuLeptonCombyReco();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int ResetEvent(PHCompositeNode *topNode);
    int Reset(PHCompositeNode *topNode) { return 0; }
    int End(PHCompositeNode *topNode);
    void Print(const char *what) const { return; }

private:
    void InitParams();
    int InitUltraLight(PHCompositeNode *topNode);
    const std::string GetFilePath();
    float get_board(float phi_dc);
    void init_rungroup();
    int get_rungroup(int run_num);
    void init_dead_area();
    bool left_bound(float x, float y, float xx1, float yy1, float xx2, float yy2);
    bool right_bound(float x, float y, float xx1, float yy1, float xx2, float yy2);
    bool dead_region(float x, float y, float xx1, float yy1, float xx2, float yy2, float xx3, float yy3, float xx4, float yy4);
    bool IsCentralSupportCut(const float theta0, const float bbcVertex);
    int applySingleTrackCut(const PHCentralTrack *d_trk, const int itrk, const float vertex, const float centrality, const int run_number);
    template <typename track>
    void set_track(track *newTrack, const PHCentralTrack *trk, const int itrk_reco, const float bbcz = 0, const float svxz = 0, const int rg_beamoffset = 0);
    void fill_SVXHits_to_myevent(const SvxClusterList *svxhitlist, MyDileptonAnalysis::MyEvent *event);


protected:
    UltraLight *ul;
    MyDileptonAnalysis::MyEventContainer *event_container;

    int remove_hadron_hits, fill_QA_hadron_hists, fill_QA_lepton_hists, fill_TTree, fill_d_dphi_hists, 
        fill_DCA_hists, use_d_dphi_DCA, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto;
    int dcmap_runs[N_RUN_GRP][MAX_RUN];
    float dcmap_xx1[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_yy1[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_xx2[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_yy2[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_xx3[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_yy3[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_xx4[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];
    float dcmap_yy4[N_RUN_GRP][N_SIDE][N_ARM][MAX_DEAD_AREA];

    int TRIGGERBIT;
    float BBC_VERTEX_CUT;

    // single cut
    int QUALITY[3];
    float Z_GLOBAL; // <|cut value|
    int DC_DEADMAP; // 0--w/o dc dead map, 1--w/ dc dead map

    // eid cut
    float E_PT;
    float MAX_PT;
    float N0;        // >cut value
    float DISP;      // <cut value
    float CHI2_NPE0; //< cut value

    float EOVERP;
    float DEP[2];
    float PROB; // >cut value
    float TEMC;
    float EMCDPHI;
    float EMCDZ;
    float EMCSDPHI;
    float EMCSDZ;
    float RICH_GHOST; // RICH ghost cut

    int ncalls, npassed;
    std::string outfilename;
};

#endif /* __DATAANALYZER3DRECO_H__ */