#include <iostream>
//#include <iomanip>
//#include <vector>
//#include <map>
//#include <algorithm>
//#include <cstdlib>
//#include "gsl/gsl_rng.h"
//#include <math.h>

#include "TMath.h"
#include "phool.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"
#include "recoConsts.h"

#include "embedana.h"

#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"
#include "EventHeader.h"
#include "VtxOut.h"
#include "PHGlobal.h"
//#include "TriggerHelper.h"
//#include "RunHeader.h"
//#include <PreviousEvent.h>
//#include "TrigLvl1.h"
//#include "TrigRunLvl1.h"
//#include "Bbc.hh"
//#include "BbcOut.h"

#include <DchHitLineTable.hh>
#include <PadCluster.h>
#include <CrkHit.h>
#include <PHTrackOut.h>
#include <PHDchTrackOut.h>
#include <CglTrack.h>

//#include "fkinWrapper.h"

//#include "McEvalSingleList.h"

//#include "SvxRawhitList.h"
//#include "SvxRawhit.h"
//#include "SvxSegmentList.h"
//#include "SvxSegment.h"
//#include "SvxTracker.h"
//#include "SvxRawhitClusterList.h"
//#include "SvxRawhitCluster.h"
#include "SvxClusterList.h"
#include "SvxGhitClusterList.h"
#include "SvxCluster.h"
#include "SvxGhitCluster.h"
//#include "SvxCentralTrackList.h"
//#include "SvxCentralTrack.h"
//#include "SvxClusterInfo.h"
//#include "SvxResidualInfo.h"
//#include "SvxBeamCenterPar.h"
//#include "SvxEventInfo.h"

//#include "SvxGhitList.h"
//#include "SvxGhit.h"
//#include "SvxGhitRawhitList.h"
//#include "SvxGhitRawhit.h"

//#include "SvxClusterContainer.h"

//#include "svxAddress.hh"
//#include "svxDetectorGeo.hh"
//#include "SvxSensor.h"
//
//#include "compactCNT/SvxHitMap.h"
//#include "compactCNT/SvxTrackMap.h"
//#include <SvxCentralTrackRecalList.h>

//#include <TofwHit.h>
//#include <HbdBlobList.h>
#include <PHPointerList.h>
#include <PHEmbedMcRecoTrack.hh>

#include "TFile.h"
#include "TNtuple.h"
//#include "TTree.h"
//#include "TF1.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TH3.h"
//#include "TLorentzVector.h"
//#include "TVector3.h"
//#include "TRandom2.h"
//#include "TMath.h"

#include "getClass.h"


using namespace std;
using namespace findNode;

//==============================================================

embedana::embedana(string filename) : m_outFileName(filename)
{
  ThisName   = "embedana";
  EventNumber=0;
  event_container = nullptr;
  fill_TTree = 0;
}

//==============================================================

embedana::~embedana() {
  delete event_container;
    std::cout<<"Event Container was deleted"<<std::endl;
}

//==============================================================

int embedana::Init(PHCompositeNode *topNode) {

    cout << "embedana::Init started..." << endl;
    //OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
    //cout << "embedana::Init: output file " << OutputFileName << " opened." << endl;
    recoConsts* rc  = recoConsts::instance(); 
    const int remove_hadron_hits = rc->get_IntFlag("Remove_hadron_hits", 0);
    const int fill_QA_hadron_hists = rc->get_IntFlag("Fill_QA_hadron_hists", 0);
    const int fill_QA_lepton_hists = rc->get_IntFlag("Fill_QA_lepton_hists", 0);
    fill_TTree = rc->get_IntFlag("Fill_TTree", 0);
    const int fill_d_dphi_hists = rc->get_IntFlag("Fill_d_dphi_hists", 0);
    const int fill_DCA_hists = rc->get_IntFlag("Fill_DCA_hists", 0);
    const int use_iden = rc->get_IntFlag("Use_ident", 0);
    const int do_track_QA = rc->get_IntFlag("Do_track_QA", 0);
    const int do_reveal_hadron = rc->get_IntFlag("Do_reveal_hadron", 0);
    const int fill_true_DCA = rc->get_IntFlag("Fill_true_DCA", 0);
    const int check_veto = rc->get_IntFlag("Check_Veto", 0);


    std::cout<<"Remove_hadron_hits: "<<remove_hadron_hits<<std::endl;
    std::cout<<"fill_QA_hadron_hists: "<<fill_QA_hadron_hists<<std::endl;
    std::cout<<"fill_QA_lepton_hists: "<<fill_QA_lepton_hists<<std::endl;
    std::cout<<"fill_TTree: "<<fill_TTree<<std::endl;
    std::cout<<"fill_d_dphi_hists: "<<fill_d_dphi_hists<<std::endl;
    std::cout<<"fill_DCA_hists: "<<fill_DCA_hists<<std::endl;
    std::cout<<"use_iden: "<<use_iden<<std::endl;
    std::cout<<"Do_track_QA: "<<do_track_QA<<std::endl;
    std::cout<<"do_reveal_hadron: "<<do_reveal_hadron<<std::endl;
    std::cout<<"fill_true_DCA: "<<fill_true_DCA<<std::endl;
    std::cout<<"check_veto: "<<check_veto<<std::endl;

    TOAD toad("Run14AuAuLeptonComby");

    const std::string loc = toad.location("field_map.root");
    
    event_container = new MyDileptonAnalysis::MyEventContainer();
    event_container->InitEvent();
    event_container->GetHistsFromFile(loc);
    event_container->CreateOutFileAndInitHists(m_outFileName.c_str(),fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,
                                               fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto);
    
    if(fill_TTree) event_container->ResetTree();
    
  cout << "embedana::Init ended." << endl;
  return 0;
}

//==============================================================
  
int embedana::InitRun(PHCompositeNode *topNode) {
  cout << "embedana::InitRun started..." << endl;

/*
  // check magnet current 
  RunHeader* runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runheader==NULL) {
    cout << PHWHERE<< "Can't find runheader. " << endl;
    return ABORTRUN;
  }
  if(runheader->get_currentCentral()>0){
   m_fieldScale = 1.0;
  } else {
   m_fieldScale = -1.0;
  }
  cout<<"SvxCentralTrackReco::InitRun  fieldScale="<<m_fieldScale<<endl;
*/

  //TofwHit        *tofw = findNode::getClass<TofwHit>(       topNode, "TofwHit");
  //SvxClusterList *svx  = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
  //HbdBlobList    *hbd  = findNode::getClass<HbdBlobList>(   topNode, "HbdBlobList");
  //cout<<"tofw "<<(tofw==NULL?"doesnot":"")<<" exist in merged topNode"<<endl;
  //cout<<"svx  "<<(svx ==NULL?"doesnot":"")<<" exist in merged topNode"<<endl;
  //cout<<"hbd  "<<(hbd ==NULL?"doesnot":"")<<" exist in merged topNode"<<endl;

  m_outfile = new TFile(m_outFileName.c_str(), "RECREATE");

  m_ntp_embed = new TNtuple("ntp_embed","embedding",
      "zvtx"
      ":pt:charge:mom:phi0:the0:dcqual"
      ":emcdz:emcdphi:pc3dz:pc3dphi"
      ":ecore:prob:n0:npe0:ch2npe0:disp"
      ":dca2d:dcaz:chi2ndf:hitptn"
      ":simpt:simmom:simphi0:simthe0"
      ":simecore:simprob:simn0:simnpe0:simch2npe0:simdisp"
      ":genpt:genmom:genphi0:genthe0:genvx:genvy:genvz:genpid:gendca2d"
    );

  cout << "embedana::InitRun ended." << endl;
  return 0;
}

//==============================================================
int embedana::ResetEvent(PHCompositeNode *topNode) {
    event_container->ClearEvent();
    return 0;
}
//==============================================================

int embedana::process_event(PHCompositeNode *topNode) {

  static int first_event = 0;
  if(first_event==0){
    cout<<"embedana::process_event firstevent"<<endl;
    topNode->print();
    first_event++;
  }

  // event info
  //--EventHeader*  evthdr  = getClass<EventHeader>(topNode,"EventHeader");
  VtxOut*       vtxout  = getClass<VtxOut>(     topNode,"VtxOut");

  // track info
  //--DchHitLineTable     *dchit   = getClass<DchHitLineTable>(topNode,"DchHitLineTable");
  //PadCluster          *pc1     = getClass<PadCluster>(     topNode,"Pc1Cluster");
  //PadCluster          *pc2     = getClass<PadCluster>(     topNode,"Pc2Cluster");
  //PadCluster          *pc3     = getClass<PadCluster>(     topNode,"Pc3Cluster");
  //CrkHit              *crkhit  = getClass<CrkHit>(         topNode,"CrkHit");
  PHCentralTrack      *trk     = getClass<PHCentralTrack>( topNode,"PHCentralTrack");
  //PHTrackOut          *phtrk   = getClass<PHTrackOut>(     topNode,"PHTrackOut");
  PHDchTrackOut       *phdchtrk= getClass<PHDchTrackOut>(  topNode,"PHDchTrackOut");
  //SvxCentralTrackList *svxcnt  = getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");
  
  
  recoConsts*    rc = recoConsts::instance();
  Fun4AllServer* se = Fun4AllServer::instance();
  PHCompositeNode* mcnode   = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCentralTrack  *trk_mc    = getClass<PHCentralTrack>(mcnode,"PHCentralTrack");

  SvxClusterList      *svx     = getClass<SvxClusterList>( topNode,"SvxClusterList");
  SvxClusterList      *svxsim     = getClass<SvxClusterList>( mcnode,"SvxClusterList");
  SvxGhitClusterList  *svxembed  = getClass<SvxGhitClusterList>( topNode,"SvxGhitClusterList");
  
  std::cout << "real embed and sim Nhits and Ntraks: " << svx->get_nClusters() << " " << svxsim->get_nClusters() << " " <<
  svxembed->get_nGhitClusters() << " " << trk->get_npart() << " " << trk_mc->get_npart() << std::endl;
  //cout<<"event : "<<EventNumber<<"  "<<((evthdr!=NULL) ? evthdr->get_EvtSequence() : -1 ) <<endl;
  
  // topnode
  //--recoConsts*    rc = recoConsts::instance();
  //--Fun4AllServer* se = Fun4AllServer::instance();
  //--PHCompositeNode* mcnode   = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  //--PHCompositeNode* realnode = se->topNode(rc->get_CharFlag("EMBED_REAL_TOPNODE"));
  if(EventNumber>=10) {
    CglTrack   *cgl   = getClass<CglTrack>(  topNode,"CglTrack");   cgl->ShutUp(1);
    PHTrackOut *phtrk = getClass<PHTrackOut>(topNode,"PHTrackOut"); phtrk->ShutUp(1);
  }

    
  //EventHeader    *evthdr_mc = getClass<EventHeader>(   mcnode,"EventHeader");
  //--DchHitLineTable *dchit_mc  = getClass<DchHitLineTable>(mcnode,"DchHitLineTable");
  //--PadCluster      *pc1_mc    = getClass<PadCluster>(    mcnode,"Pc1Cluster");
  //--PadCluster      *pc2_mc    = getClass<PadCluster>(    mcnode,"Pc2Cluster");
  //--PadCluster      *pc3_mc    = getClass<PadCluster>(    mcnode,"Pc3Cluster");
  //--CrkHit          *crkhit_mc = getClass<CrkHit>(        mcnode,"CrkHit");
  //--PHCentralTrack  *trk_mc    = getClass<PHCentralTrack>(mcnode,"PHCentralTrack");
  //--PHTrackOut      *phtrk_mc  = getClass<PHTrackOut>(    mcnode,"PHTrackOut");
    
  //EventHeader    *evthdr_real = getClass<EventHeader>(   realnode,"EventHeader");
  //--DchHitLineTable *dchit_real  = getClass<DchHitLineTable>(realnode,"DchHitLineTable");
  //--PadCluster      *pc1_real    = getClass<PadCluster>(    realnode,"Pc1Cluster");
  //--PadCluster      *pc2_real    = getClass<PadCluster>(    realnode,"Pc2Cluster");
  //--PadCluster      *pc3_real    = getClass<PadCluster>(    realnode,"Pc3Cluster");
  //--CrkHit          *crkhit_real = getClass<CrkHit>(        realnode,"CrkHit");
  //--PHCentralTrack  *trk_real    = getClass<PHCentralTrack>(realnode,"PHCentralTrack");
  //--PHTrackOut      *phtrk_real  = getClass<PHTrackOut>(    realnode,"PHTrackOut");
    
    
/*
  if(dchit_mc   !=NULL) cout<<"dchit nhit : "<<dchit_mc->Entries()<<flush;
  if(dchit_real !=NULL) cout<<          " : "<<dchit_real->Entries()<<flush;
  if(dchit      !=NULL) cout<<          " : "<<dchit->Entries()<<flush;
  cout<<endl;
  if(pc1_mc   !=NULL) cout<<"pc1 nhit : "<<pc1_mc->get_PadNCluster()<<flush;
  if(pc1_real !=NULL) cout<<        " : "<<pc1_real->get_PadNCluster()<<flush;
  if(pc1      !=NULL) cout<<        " : "<<pc1->get_PadNCluster()<<flush;
  cout<<endl;
  if(pc2_mc   !=NULL) cout<<"pc2 nhit : "<<pc2_mc->get_PadNCluster()<<flush;
  if(pc2_real !=NULL) cout<<        " : "<<pc2_real->get_PadNCluster()<<flush;
  if(pc2      !=NULL) cout<<        " : "<<pc2->get_PadNCluster()<<flush;
  cout<<endl;
  if(pc3_mc   !=NULL) cout<<"pc3 nhit : "<<pc3_mc->get_PadNCluster()<<flush;
  if(pc3_real !=NULL) cout<<        " : "<<pc3_real->get_PadNCluster()<<flush;
  if(pc3      !=NULL) cout<<        " : "<<pc3->get_PadNCluster()<<flush;
  cout<<endl;
  if(crkhit_mc  !=NULL) cout<<"crk nhit : "<<crkhit_mc->get_CrkNHit()<<flush;
  if(crkhit_real!=NULL) cout<<        " : "<<crkhit_real->get_CrkNHit()<<flush;
  if(crkhit     !=NULL) cout<<        " : "<<crkhit->get_CrkNHit()<<flush;
  cout<<endl;
  cout<<"cnt ntrk : "<<((trk_mc  !=NULL) ? trk_mc->get_npart()   : 999999)<<flush;
  cout<<        " : "<<((trk_real!=NULL) ? trk_real->get_npart() : 999999)<<flush;
  cout<<        " : "<<((trk     !=NULL) ? trk->get_npart()      : 999999)<<flush;
  cout<<endl;
  if(phtrk_mc  !=NULL) cout<<"phtrk ntrk : "<<phtrk_mc->get_PHNTrack()<<flush;
  if(phtrk_real!=NULL) cout<<          " : "<<phtrk_real->get_PHNTrack()<<flush;
  if(phtrk     !=NULL) cout<<          " : "<<phtrk->get_PHNTrack()<<flush;
  cout<<endl;
  cout<<        "svxcnttrk : "<<((svxcnt  !=NULL) ? svxcnt->get_nCentralTracks()      : 999999)<<endl;
  cout<<        "svxcls    : "<<((svx     !=NULL) ? svx->get_nClusters()      : 999999)<<endl;
*/

  //////////////////////////////////////////////////////////////////////
  PHPointerList<PHEmbedMcRecoTrack>  *embedtrk = getClass< PHPointerList<PHEmbedMcRecoTrack> >(topNode, "PHEmbedMcRecoTrack");
  MyDileptonAnalysis::MyEvent *event = event_container->GetEvent();
  event->SetPreciseX((vtxout->get_Vertex()).getX());
  event->SetPreciseY((vtxout->get_Vertex()).getY());
  event->SetPreciseZ((vtxout->get_Vertex()).getZ());
  //std::cout<<(vtxout->get_Vertex()).getX()<<" "<<(vtxout->get_Vertex()).getY()<<" "<<(vtxout->get_Vertex()).getZ()<<std::endl;
  //event->ClearEvent();

  // mapping between PHCentralTrack and SvxCentralTrack
  //remapCntSvx(trk, svxcnt);

  if(embedtrk==nullptr) { cout<<" no PHEmbedMcRecoTrack"<<endl;}
  else {
    //--cout<<"PHEmbedMcRecoTrack: length = "<<embedtrk->length()<<endl;
    for(unsigned int i=0; i<embedtrk->length(); i++){
      PHEmbedMcRecoTrack *embed = (*embedtrk)[i];
      //--cout<<"   id_g, _s, _r = "<<embed->get_dctrkidG()<<" "<<embed->get_dctrkidS()<<" "<<embed->get_dctrkidR()<<" ";
      //--cout<<" mom="<<embed->get_momG()<<" "<<embed->get_momR()<<endl;

      float momr = -9999, momr1=-9999;
      int dcidr = embed->get_dctrkidR();
      if(0<=dcidr&&dcidr<(int)trk->get_npart()){
         PHSnglCentralTrack* sngl = trk->get_track(dcidr);
         momr = sngl->get_mom();
         momr1= phdchtrk->get_momentum(dcidr);
         //--cout<<"  reco mom= "<<momr<<" "<<momr1<<endl;
      }
      else {
        cout<<"out of range : "<<dcidr<<momr<<momr1<<endl;
      }
    }

    int Ncnt = trk->get_npart();

    // ntp fill
    unsigned int Nembed = (unsigned int)(0.5*embedtrk->length());

    for(unsigned int itrk=0; itrk<Nembed; itrk++){
      PHEmbedMcRecoTrack *embed = (*embedtrk)[itrk];
      //int idG = embed->get_dctrkidG();
      int idR = embed->get_dctrkidR();
      //int idS = embed->get_dctrkidS();
      std::cout<<embed->get_bbccent()<<std::endl;

      if(idR<0||Ncnt<=idR) {
        cout<<" RealtrackID is out of range = "<<idR<<" : "<<Ncnt<<endl;
        continue;
      }

      PHSnglCentralTrack* sngl       = trk->get_track(idR);
      ///SvxCentralTrack*    svxcntsngl = m_vsvxcnt[idR];
      if(false)
      {
        PHSnglCentralTrack* mcsngl       = trk_mc->get_track(0);
        std::cout << "mom of embed, mc, embed sim, geant " << sngl->get_mom() << " " << mcsngl->get_mom() << " " << embed->get_momS() << " " << embed->get_momG() << std::endl;
      }

      float ntp[100];
      for(int i=0; i<100; i++) ntp[i] = -9999.;

      float charge = sngl->get_charge();

      ntp[ 0] = vtxout->get_ZVertex(); //zvtx
      ntp[ 1] = sngl->get_mom()*sin(sngl->get_the0());
      ntp[ 2] = charge;
      ntp[ 3] = sngl->get_mom();
      ntp[ 4] = sngl->get_phi0();
      ntp[ 5] = sngl->get_the0();
      ntp[ 6] = sngl->get_quality();
      ntp[ 7] = sngl->get_emcdz();
      ntp[ 8] = sngl->get_emcdphi();
      ntp[ 9] = sngl->get_pc3dz();
      ntp[10] = sngl->get_pc3dphi();
      ntp[11] = sngl->get_ecore();
      ntp[12] = sngl->get_prob();
      ntp[13] = sngl->get_n0();
      ntp[14] = sngl->get_npe0();
      ntp[15] = sngl->get_chi2()/sngl->get_npe0();
      ntp[16] = sngl->get_disp();
      
      //if(svxcntsngl!=NULL){
      //  float ndf        = svxcntsngl->getNDF();
      //  float chisqr     = svxcntsngl->getChiSquare();
      //  ntp[17] = svxcntsngl->getDCA2D(); // dca2d
      //  ntp[18] = svxcntsngl->getDCAZ();  // dcaz
      //  ntp[19] = (ndf>0) ? chisqr/ndf : -9999.; // chi2ndf
      //  ntp[20] = svxcntsngl->getHitPattern(); // hitptn
      //}

      ntp[21] = embed->get_ptS();
      ntp[22] = embed->get_momS();
      ntp[23] = embed->get_phi0S();
      ntp[24] = embed->get_the0S();
      ntp[25] = embed->get_emcecoreS();
      ntp[26] = embed->get_emcprobphotS();
      ntp[27] = embed->get_crknpmt0S();
      ntp[28] = embed->get_crknpe0S();
      ntp[29] = embed->get_crkchi2S()/embed->get_crknpe0S();
      ntp[30] = embed->get_crkdispS();

      // use gen info
      float genpt  = embed->get_ptG();
      float genphi = embed->get_phi0G();
      float genthe = embed->get_the0G();
      float genvx  = embed->get_xvtxG();
      float genvy  = embed->get_xvtxG()*TMath::Tan(embed->get_phi0G());
      float genvz  = embed->get_zvtxG();

      float beamxy[2] = {0., 0.};
      float fieldPol  = 1.0;
      float gendca[3], gendca2d;

      calcDCA_BCbyCircleProjection
            (
              genpt, genphi, genthe,  // pt in xy-plane, phi, theta at inner most layer
              charge,                 // charge of the track
              genvx, genvy, genvz,    // hit position at inner most layer
              beamxy[0], beamxy[1],   // beam center
              fieldPol,               //
              &gendca[0], &gendca[1], &gendca[2], // dca position
              &gendca2d               // return
            );

      ntp[31] = embed->get_ptG();
      ntp[32] = embed->get_momG();
      ntp[33] = embed->get_phi0G();
      ntp[34] = embed->get_the0G();
      ntp[35] = embed->get_xvtxG();
      ntp[36] = embed->get_xvtxG()*TMath::Tan(embed->get_phi0G());
      ntp[37] = embed->get_zvtxG();
      ntp[38] = embed->get_partidG();
      ntp[39] = gendca2d;

      m_ntp_embed->Fill(ntp);

      MyDileptonAnalysis::MyElectron newElectron;

      newElectron.SetPt(embed->get_momS()*sin(embed->get_the0S()));
      newElectron.SetPtPrime(sngl->get_mom()*sin(sngl->get_the0()));
      newElectron.SetReconPT(embed->get_momG()*sin(embed->get_the0G()));
      newElectron.SetTrkId(itrk);
      newElectron.SetTrkQuality(trk->get_quality(itrk));
      newElectron.SetArm(trk->get_dcarm(itrk));
      newElectron.SetDCSide(trk->get_dcside(itrk));
      newElectron.SetSect(trk->get_sect(itrk));
      newElectron.SetQ(trk->get_charge(itrk));
      newElectron.SetPhiDC(trk->get_phi(itrk));
      newElectron.SetPhi0(trk->get_phi0(itrk));
      newElectron.SetThe0(trk->get_the0(itrk));
      newElectron.SetPhi0Prime(trk->get_phi0(itrk));
      newElectron.SetThe0Prime(trk->get_the0(itrk));
      newElectron.SetZDC(trk->get_zed(itrk));
      newElectron.SetAlpha(trk->get_alpha(itrk));
      newElectron.SetAlphaPrime(trk->get_alpha(itrk));
      newElectron.SetEmcId(trk->get_emcid(itrk));
      newElectron.SetEcore(trk->get_ecore(itrk));
      newElectron.SetDep(trk->get_dep(itrk));
      newElectron.SetProb(trk->get_prob(itrk));
      newElectron.SetEmcdz(trk->get_emcdz(itrk));
      newElectron.SetEmcdphi(trk->get_emcdphi(itrk));
      newElectron.SetEmcTower(trk->get_sect(itrk), trk->get_ysect(itrk), trk->get_zsect(itrk));
      newElectron.SetTOFDPHI(trk->get_n0(itrk));
      newElectron.SetTOFDZ(trk->get_plemc(itrk));
      newElectron.SetPC3SDPHI(trk->get_pc3sdphi(itrk));
      newElectron.SetPC3SDZ(trk->get_pc3sdz(itrk));
      newElectron.SetCrkphi(trk->get_center_phi(itrk));
      newElectron.SetCrkz(trk->get_center_z(itrk));
      newElectron.SetTOFE((trk->get_m2tof(itrk)));
      newElectron.SetEmcTOF(trk->get_temc(itrk));
      newElectron.SetChi2(trk->get_chi2(itrk));
      newElectron.SetN0(trk->get_n0(itrk));
      newElectron.SetNPE0(trk->get_npe0(itrk));
      newElectron.SetDISP(trk->get_disp(itrk));
      newElectron.SetEmcdz_e(trk->get_emcsdz_e(itrk));
      newElectron.SetEmcdphi_e(trk->get_emcsdphi_e(itrk));
      newElectron.SetMcId(embed->get_partidG());

      event->AddTrack(&newElectron);
    }
    
  }
  for (int ihit = 0; ihit < svxsim->get_nClusters(); ihit++)
    {

        SvxCluster *svxhit = svxsim->get_Cluster(ihit);

        if (svxhit == nullptr)
        {
            std::cout << "cluster NULL : " << ihit << std::endl;
            continue;
        }
        
        MyDileptonAnalysis::MyVTXHit newHit;
        
        newHit.SetClustId(ihit);
        newHit.SetLayer(svxhit->get_layer());
        newHit.SetLadder(svxhit->get_ladder());
        newHit.SetSensor(1);
        newHit.SetXHit(svxhit->get_xyz_global(0));
        newHit.SetYHit(svxhit->get_xyz_global(1));
        newHit.SetZHit(svxhit->get_xyz_global(2));
        if(false)std::cout<<newHit.GetXHit()<<" "<<newHit.GetYHit()<<" "<<newHit.GetZHit()<<std::endl;
        newHit.SetiLayerFromR();
        if( svxhit->get_layer()!=newHit.GetLayer()||svxhit->get_ladder()!=newHit.GetLadder()||
        newHit.GetSensor()!=1) 
        {
            std::cout<<" smth is wrong "<<std::endl;
            continue;
        }
        event->AddVTXHit(&newHit);
    }
    for (int ihit = 0; ihit < svxembed->get_nGhitClusters(); ihit++)
    {

        SvxGhitCluster *svxembeshit = svxembed->get_GhitCluster(ihit);
        SvxCluster *svxhit = svx->get_Cluster(svxembeshit->get_clusterID());

        if (svxhit == nullptr)
        {
            std::cout << "cluster NULL : " << ihit << std::endl;
            continue;
        }
        
        MyDileptonAnalysis::MyVTXHit newHit;
        
        newHit.SetClustId(ihit);
        newHit.SetLayer(svxhit->get_layer());
        newHit.SetLadder(svxhit->get_ladder());
        newHit.SetSensor(0);
        newHit.SetXHit(svxhit->get_xyz_global(0));
        newHit.SetYHit(svxhit->get_xyz_global(1));
        newHit.SetZHit(svxhit->get_xyz_global(2));
        if(false)std::cout<<newHit.GetXHit()<<" "<<newHit.GetYHit()<<" "<<newHit.GetZHit()<<std::endl;
        newHit.SetiLayerFromR();
        if( svxhit->get_layer()!=newHit.GetLayer()||svxhit->get_ladder()!=newHit.GetLadder()||
        newHit.GetSensor()!=0) 
        {
            std::cout<<" smth is wrong "<<std::endl;
            continue;
        }
        event->AddVTXHit(&newHit);
    }

    if(fill_TTree) event_container->FillTree();
  // dchit
/*
  for(int ihit=0; ihit<dchit->Entries(); ihit++){
    PHPoint pos = dchit->getXYZ(ihit);
    cout<<"pos : "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<endl;
  }
*/
/*
  //  cout<<"aa"<<endl;
  //  PHTrackOut         *proj = getClass<PHTrackOut>(topNode,"PHTrackOut");
  //  CrkRing             *crk = getClass<CrkRing>(topNode,"CrkRing");
  //  emcClusterContainer *emc = getClass<emcClusterContainer>(topNode,"emcClusterContainer");
  //  DchTrack            *dch = getClass<DchTrack>(topNode,"DchTrack");
  //  CglTrack            *cgl = getClass<CglTrack>(topNode,"CglTrack");
  //TrigRunLvl1 *lvl1    = getClass<TrigRunLvl1>(topNode,"TrigRunLvl1");
  PHGlobal    *global       = getClass<PHGlobal>(topNode,"PHGlobal");
  RunHeader   *run          = getClass<RunHeader>(topNode,"RunHeader");
  TrigLvl1    *trglvl1      = getClass<TrigLvl1>(topNode,"TrigLvl1");
  BbcOut *bbc               = getClass<BbcOut>(topNode,"BbcOut");

  PreviousEvent *peve = getClass<PreviousEvent>(topNode, "PreviousEvent");
  if(EventNumber==0)
    if(peve==NULL){ cout<<"no PreviousEvent"<<endl; }


  int   RunNumber = (run!=NULL)          ? run->get_RunNumber()            : -9999;
  int   EventNum  = (event_header!=NULL) ? event_header->get_EvtSequence() : -9999;
  float cent      = (global!=NULL)       ? global->getCentrality()         : -9999;

  if(EventNumber==0) {
    cout << "run = "<<RunNumber<<" event = "<<EventNum<<endl;
  }
*/


/*

  SvxRawhitList  *svxrawhit = getClass<SvxRawhitList>(topNode,"SvxRawhitList");
  SvxClusterList *svxsel    = getClass<SvxClusterList>(topNode,"SvxSelectedClusterList");
  SvxSegmentList *svxtracks = getClass<SvxSegmentList>(topNode,"SvxSegmentList");
  SvxRawhitClusterList 
             *svxrawhitclus = getClass<SvxRawhitClusterList>(topNode,"SvxRawhitClusterList");

  SvxCentralTrackList
             *svxcnttrklist = getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");


  McEvalSingleList *mctrk   = getClass<McEvalSingleList>(topNode,"McSingle");
  SvxGhitList      *svxghit = getClass<SvxGhitList>(topNode,"SvxGhitList");
  fkinWrapper      *fkinW   = getClass<fkinWrapper>(topNode, "fkin");
  SvxGhitRawhitList *svxghitrawlist = getClass<SvxGhitRawhitList>(topNode,"SvxGhitRawhitList");
  if(EventNumber==0){
    if(m_simmode){
      if(mctrk==NULL){
        cout<<"No McSingle Info in Nodetree"<<endl;
      }
      if(svxghit==NULL){
        cout<<"No SvxGhitList Info in Nodetree"<<endl;
      }
      if(fkinW==NULL){
        cout<<"No fkinWrapper Info in Nodetree"<<endl;
      }
      if(svxghitrawlist==NULL){
        cout<<"No SvxGhitRawhitList Info in Nodetree"<<endl;
      }
    }
  }

  int ntrk       = (trk!=NULL) ? trk->get_npart() : 0;
  int nsvxtracks = (svxtracks!=NULL) ? svxtracks->get_nSegments() : -9999;
  int nsvx       = (svx!=NULL)? svx->get_nClusters() : -9999;
  int nsvxrawhit = (svxrawhit!=NULL) ? svxrawhit->get_nRawhits() : -9999;
  //int nsvxrawhitclus = 0;
  int nkalfit = 0;
  float zvtx = (bbc!=NULL) ? bbc->get_VertexPoint() : -999.0;
  if(zvtx<-900){
    zvtx = (vtxout!=NULL) ? vtxout->get_ZVertex("BBC") : -999.0;
  }
  float bbcq = 50;
  if(bbc!=NULL) {
    bbcq = bbc->get_ChargeSum(Bbc::North) + bbc->get_ChargeSum(Bbc::South);
  }
  else if(global!=NULL) {
    bbcq = global->getBbcChargeN() + global->getBbcChargeS();
  }
  //if(bbcq>200 || bbcq<10)return 1; 

  float vtx_z = (vtxout!=NULL) ? vtxout->get_ZVertex() : -999.;
  if(m_goodOnlyFlag&&(fabs(vtx_z)>GOODZVTX)){ // skip this event
    if(verbosity>0) cout << "embedana  skip this event : run = "<<RunNumber<<" event = "<<EventNum<<endl;
    return EVENT_OK;
  }

  

  if(init_ana==0) {
    cout<<"init ana"<<endl;
    topNode->print(); 
    init_ana=1;
 
    if(svx!=NULL)          svx->identify();          else cout<<"no SvxCluster object"<<endl; 
    if(svxtracks!=NULL)    svxtracks->identify();    else cout<<"no SvxSegment object"<<endl;
    if(svxrawhitclus!=NULL)svxrawhitclus->identify();else cout<<"no SvxRawhitClus object"<<endl;
    if(svxrawhit!=NULL)    svxrawhit->identify();    else cout<<"no SvxRawhit object"<<endl;
    if(trglvl1!=NULL)      trglvl1->identify();      else cout<<"no TrigLvl1 object"<<endl;

    cout<<endl<<endl;
    cout<<"init ana ends"<<endl;
  }

  if(verbosity>1){
    if(EventNumber%1==0) {
      cout << "------- Event # " << EventNumber << " nsvxrawhit= " << nsvxrawhit << " nsvx= " << nsvx << " ntrack= " << ntrk << " nkalfit= " << nkalfit << " nsvxtracks=" << nsvxtracks << endl;
    }
  }

  if((EventNumber%1000)==0){
    cout << "------- Event # " << EventNumber << endl;
  }
*/

  EventNumber++;
  return 0;
}

int embedana::End(PHCompositeNode *topNode)
{
  cout << "Writing out..." << endl;
  m_outfile->Write();
  cout << "Closing output file..." << endl;
  m_outfile->Close();
  delete m_outfile;
  event_container->WriteOutFile();
  
  return 0;
}

void embedana::calcDCA_BCbyCircleProjection
(
 float pt, float phi, float the,  // pt in xy-plane, phi, theta at inner most layer
 int charge,                      // charge of the track
 float hx, float hy, float hz,    // hit position at inner most layer
 float vx, float vy,              // beam center
 float fieldPolarity,             //
 float* dx, float* dy, float* dz, // dca position
 float* d2dca_bc)                 // return
{
  static const float B = 0.90;
  static const float b = 0.003*B;

  // vector to rotation center at hit1 
  float R  = pt/b;
  float pz = pt/tan(the);

  float b_sign = (fieldPolarity>0) ? -1.0 : 1.0;
  float dir = ( (b_sign)*charge>0. ) ? -1.0 : 1.0;
  float cx = hx - dir*R*sin(phi);
  float cy = hy + dir*R*cos(phi);

  // L is a distance btw the rotation center and primary vtx  
  // psi is a angle of the vector starting from the rotation center to primary vertex
  float L = sqrt((vx-cx)*(vx-cx) + (vy-cy)*(vy-cy));
  float psi = atan2((vy-cy), (vx-cx));

  // DCA point
  *dx = cx + R*cos(psi);
  *dy = cy + R*sin(psi);

  float dphi = phi - dir*0.5*M_PI - psi;
  if ( dphi>M_PI*1.5 ) dphi -= M_PI*2.;
  else if ( dphi<-M_PI*1.5 ) dphi += M_PI*2.;
  float dzdphi = dir*pz*R/pt;
  *dz = hz + dphi*dzdphi;

  // DCA value
  *d2dca_bc = b_sign*(charge*(R - L));
}
  
