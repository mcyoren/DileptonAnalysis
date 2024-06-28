#include "embedana.h"
//==============================================================

embedana::embedana(string filename, string filepath, string oscarpath) : m_outFileName(filename)
{
  ThisName = "embedana";
  EventNumber = 0;
  event_container = nullptr;
  fill_TTree = 0;
  remove_hadron_hits = 0;
  local_filepath = filepath;
  local_oscarpath = oscarpath;
  vertexes.clear();

  memset(InData_read, 0, 10000 * sizeof(InData));
}

//==============================================================

embedana::~embedana()
{
  delete event_container;
  std::cout << "Event Container was deleted" << std::endl;
}

//==============================================================

int embedana::Init(PHCompositeNode *topNode)
{

  cout << "embedana::Init started..." << endl;
  // OutputNtupleFile = new TFile(OutputFileName.c_str(),"RECREATE");
  // cout << "embedana::Init: output file " << OutputFileName << " opened." << endl;
  recoConsts *rc = recoConsts::instance();
  remove_hadron_hits = rc->get_IntFlag("Remove_hadron_hits", 0);
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

  std::cout << "Remove_hadron_hits: " << remove_hadron_hits << std::endl;
  std::cout << "fill_QA_hadron_hists: " << fill_QA_hadron_hists << std::endl;
  std::cout << "fill_QA_lepton_hists: " << fill_QA_lepton_hists << std::endl;
  std::cout << "fill_TTree: " << fill_TTree << std::endl;
  std::cout << "fill_d_dphi_hists: " << fill_d_dphi_hists << std::endl;
  std::cout << "fill_DCA_hists: " << fill_DCA_hists << std::endl;
  std::cout << "use_iden: " << use_iden << std::endl;
  std::cout << "Do_track_QA: " << do_track_QA << std::endl;
  std::cout << "do_reveal_hadron: " << do_reveal_hadron << std::endl;
  std::cout << "fill_true_DCA: " << fill_true_DCA << std::endl;
  std::cout << "check_veto: " << check_veto << std::endl;

  TOAD toad("Run14AuAuLeptonComby");

  const std::string loc = toad.location("field_map.root");

  event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile(loc);
  event_container->CreateOutFileAndInitHists(m_outFileName.c_str(), fill_QA_lepton_hists, fill_QA_hadron_hists, fill_TTree, fill_d_dphi_hists,
                                             fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto);

  if (fill_TTree)
    event_container->ResetTree();

  cout << "embedana::Init ended." << endl;
  return 0;
}

//==============================================================

int embedana::InitRun(PHCompositeNode *topNode)
{
  cout << "embedana::InitRun started..." << endl;

  m_outfile = new TFile(m_outFileName.c_str(), "RECREATE");

  m_ntp_embed = new TNtuple("ntp_embed", "embedding",
                            "zvtx"
                            ":pt:charge:mom:phi0:the0:dcqual"
                            ":emcdz:emcdphi:pc3dz:pc3dphi"
                            ":ecore:prob:n0:npe0:ch2npe0:disp"
                            ":dca2d:dcaz:chi2ndf:hitptn"
                            ":simpt:simmom:simphi0:simthe0"
                            ":simecore:simprob:simn0:simnpe0:simch2npe0:simdisp"
                            ":genpt:genmom:genphi0:genthe0:genvx:genvy:genvz:genpid:gendca2d");

  const int error = ReadOrigPartMoms();
  if (error)
  {
    printf("Something wrong!!!!\n");
    return 1;
  }

  std::ifstream myfile(local_filepath.c_str());
  string line;
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      string s;
      std::stringstream ss(line);
      while (getline(ss, s, ' '))
      {
        vertexes.push_back(atof(s.c_str()));
      }
    }
    myfile.close();
  }

  cout << "embedana::InitRun ended." << endl;
  return 0;
}

//==============================================================
int embedana::ResetEvent(PHCompositeNode *topNode)
{
  event_container->ClearEvent();
  return 0;
}
//==============================================================

int embedana::process_event(PHCompositeNode *topNode)
{

  static int first_event = 0;
  if (first_event == 0)
  {
    cout << "embedana::process_event firstevent" << endl;
    topNode->print();
    first_event++;
  }

  if (EventNumber > (int)vertexes.size() / 4)
    return 0;

  // event info
  //--EventHeader*  evthdr  = getClass<EventHeader>(topNode,"EventHeader");
  VtxOut *vtxout = getClass<VtxOut>(topNode, "VtxOut");

  PHCentralTrack *trk = getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  PHDchTrackOut *phdchtrk = getClass<PHDchTrackOut>(topNode, "PHDchTrackOut");

  recoConsts *rc = recoConsts::instance();
  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *mcnode = se->topNode(rc->get_CharFlag("EMBED_MC_TOPNODE"));
  PHCentralTrack *trk_mc = getClass<PHCentralTrack>(mcnode, "PHCentralTrack");

  SvxClusterList *svx = getClass<SvxClusterList>(topNode, "SvxClusterList");
  SvxClusterList *svxsim = getClass<SvxClusterList>(mcnode, "SvxClusterList");
  SvxGhitClusterList *svxembed = getClass<SvxGhitClusterList>(topNode, "SvxGhitClusterList");

  std::cout << "real embed and sim Nhits and Ntraks: " << svx->get_nClusters() << " " << svxsim->get_nClusters() << " " << svxembed->get_nGhitClusters() << " " << trk->get_npart() << " " << trk_mc->get_npart() << std::endl;
  // cout<<"event : "<<EventNumber<<"  "<<((evthdr!=NULL) ? evthdr->get_EvtSequence() : -1 ) <<endl;

  if (EventNumber >= 10)
  {
    CglTrack *cgl = getClass<CglTrack>(topNode, "CglTrack");
    cgl->ShutUp(1);
    PHTrackOut *phtrk = getClass<PHTrackOut>(topNode, "PHTrackOut");
    phtrk->ShutUp(1);
  }

  //////////////////////////////////////////////////////////////////////

  PHPointerList<PHEmbedMcRecoTrack> *embedtrk = getClass<PHPointerList<PHEmbedMcRecoTrack> >(topNode, "PHEmbedMcRecoTrack");
  MyDileptonAnalysis::MyEvent *event = event_container->GetEvent();
  event->SetPreciseX(vertexes[4 * EventNumber + 0]);
  event->SetPreciseY(vertexes[4 * EventNumber + 1]);
  event->SetPreciseZ(vertexes[4 * EventNumber + 2]);
  event->SetCentrality(vertexes[4 * EventNumber + 3]);

  MyDileptonAnalysis::MyElectron initial_particle;
  initial_particle.SetMcId(InData_read[EventNumber].id);
  const double init_pt = sqrt(InData_read[EventNumber].px * InData_read[EventNumber].px + InData_read[EventNumber].py * InData_read[EventNumber].py);
  initial_particle.SetPt(init_pt);
  initial_particle.SetPhi0(atan2(InData_read[EventNumber].py, InData_read[EventNumber].px));
  initial_particle.SetThe0(atan2(init_pt, InData_read[EventNumber].pz));
  if (false)
    std::cout << "x,y,z,cent: " << vertexes[4 * EventNumber] << " " << vertexes[4 * EventNumber + 1] << " " << vertexes[4 * EventNumber + 2] << " " << vertexes[4 * EventNumber + 3] << " " << std::endl;
  if (false)
    std::cout << "px,py,pz,id : " << InData_read[EventNumber].px << " " << InData_read[EventNumber].py << " " << InData_read[EventNumber].pz << " " << InData_read[EventNumber].id << " " << std::endl;
  if (false)
    std::cout << "px,py,pz,id : " << initial_particle.GetPx() << " " << initial_particle.GetPy() << " " << initial_particle.GetPz() << " " << initial_particle.GetMcId() << " " << std::endl;
  event->AddElecCand(&initial_particle);
  // std::cout<<(vtxout->get_Vertex()).getX()<<" "<<(vtxout->get_Vertex()).getY()<<" "<<(vtxout->get_Vertex()).getZ()<<std::endl;
  // event->ClearEvent();

  if (embedtrk == nullptr)
  {
    cout << " no PHEmbedMcRecoTrack" << endl;
  }
  else
  {
    //--cout<<"PHEmbedMcRecoTrack: length = "<<embedtrk->length()<<endl;
    for (unsigned int i = 0; i < embedtrk->length(); i++)
    {
      PHEmbedMcRecoTrack *embed = (*embedtrk)[i];
      //--cout<<"   id_g, _s, _r = "<<embed->get_dctrkidG()<<" "<<embed->get_dctrkidS()<<" "<<embed->get_dctrkidR()<<" ";
      //--cout<<" mom="<<embed->get_momG()<<" "<<embed->get_momR()<<endl;

      float momr = -9999, momr1 = -9999;
      int dcidr = embed->get_dctrkidR();
      if (0 <= dcidr && dcidr < (int)trk->get_npart())
      {
        PHSnglCentralTrack *sngl = trk->get_track(dcidr);
        momr = sngl->get_mom();
        momr1 = phdchtrk->get_momentum(dcidr);
        //--cout<<"  reco mom= "<<momr<<" "<<momr1<<endl;
      }
      else
      {
        cout << "out of range : " << dcidr << momr << momr1 << endl;
      }
    }

    int Ncnt = trk->get_npart();

    // ntp fill
    unsigned int Nembed = (unsigned int)(0.5 * embedtrk->length());

    for (unsigned int itrk = 0; itrk < Nembed; itrk++)
    {
      PHEmbedMcRecoTrack *embed = (*embedtrk)[itrk];
      // int idG = embed->get_dctrkidG();
      int idR = embed->get_dctrkidR();
      // int idS = embed->get_dctrkidS();
      if (false)
        std::cout << embed->get_bbccent() << " " << embed->get_ptG() * TMath::Cos(embed->get_phi0G()) << " " << embed->get_ptG() * TMath::Sin(embed->get_phi0G()) << " " << idR << " " << Nembed << std::endl;

      if (idR < 0 || Ncnt <= idR)
      {
        cout << " RealtrackID is out of range = " << idR << " : " << Ncnt << endl;
        continue;
      }

      PHSnglCentralTrack *sngl = trk->get_track(idR);
      if (!sngl)
        continue;
      /// SvxCentralTrack*    svxcntsngl = m_vsvxcnt[idR];
      if (false)
      {
        PHSnglCentralTrack *mcsngl = trk_mc->get_track(0);
        std::cout << "mom of embed, mc, embed sim, geant " << sngl->get_mom() << " " << mcsngl->get_mom() << " " << embed->get_momS() << " " << embed->get_momG() << std::endl;
      }

      float ntp[100];
      for (int i = 0; i < 100; i++)
        ntp[i] = -9999.;

      float charge = sngl->get_charge();

      ntp[0] = vtxout->get_ZVertex(); // zvtx
      ntp[1] = sngl->get_mom() * sin(sngl->get_the0());
      ntp[2] = charge;
      ntp[3] = sngl->get_mom();
      ntp[4] = sngl->get_phi0();
      ntp[5] = sngl->get_the0();
      ntp[6] = sngl->get_quality();
      ntp[7] = sngl->get_emcdz();
      ntp[8] = sngl->get_emcdphi();
      ntp[9] = sngl->get_pc3dz();
      ntp[10] = sngl->get_pc3dphi();
      ntp[11] = sngl->get_ecore();
      ntp[12] = sngl->get_prob();
      ntp[13] = sngl->get_n0();
      ntp[14] = sngl->get_npe0();
      ntp[15] = sngl->get_chi2() / sngl->get_npe0();
      ntp[16] = sngl->get_disp();

      // if(svxcntsngl!=NULL){
      //   float ndf        = svxcntsngl->getNDF();
      //   float chisqr     = svxcntsngl->getChiSquare();
      //   ntp[17] = svxcntsngl->getDCA2D(); // dca2d
      //   ntp[18] = svxcntsngl->getDCAZ();  // dcaz
      //   ntp[19] = (ndf>0) ? chisqr/ndf : -9999.; // chi2ndf
      //   ntp[20] = svxcntsngl->getHitPattern(); // hitptn
      // }

      ntp[21] = embed->get_ptS();
      ntp[22] = embed->get_momS();
      ntp[23] = embed->get_phi0S();
      ntp[24] = embed->get_the0S();
      ntp[25] = embed->get_emcecoreS();
      ntp[26] = embed->get_emcprobphotS();
      ntp[27] = embed->get_crknpmt0S();
      ntp[28] = embed->get_crknpe0S();
      ntp[29] = embed->get_crkchi2S() / embed->get_crknpe0S();
      ntp[30] = embed->get_crkdispS();

      // use gen info
      float genpt = embed->get_ptG();
      float genphi = embed->get_phi0G();
      float genthe = embed->get_the0G();
      float genvx = embed->get_xvtxG();
      float genvy = embed->get_xvtxG() * TMath::Tan(embed->get_phi0G());
      float genvz = embed->get_zvtxG();

      float beamxy[2] = {0., 0.};
      float fieldPol = 1.0;
      float gendca[3], gendca2d;

      calcDCA_BCbyCircleProjection(
          genpt, genphi, genthe,              // pt in xy-plane, phi, theta at inner most layer
          charge,                             // charge of the track
          genvx, genvy, genvz,                // hit position at inner most layer
          beamxy[0], beamxy[1],               // beam center
          fieldPol,                           //
          &gendca[0], &gendca[1], &gendca[2], // dca position
          &gendca2d                           // return
      );

      ntp[31] = embed->get_ptG();
      ntp[32] = embed->get_momG();
      ntp[33] = embed->get_phi0G();
      ntp[34] = embed->get_the0G();
      ntp[35] = embed->get_xvtxG();
      ntp[36] = embed->get_xvtxG() * TMath::Tan(embed->get_phi0G());
      ntp[37] = embed->get_zvtxG();
      ntp[38] = embed->get_partidG();
      ntp[39] = gendca2d;

      m_ntp_embed->Fill(ntp);

      MyDileptonAnalysis::MyElectron newElectron;

      newElectron.SetPt(embed->get_momS() * sin(embed->get_the0S()));
      newElectron.SetPtPrime(sngl->get_mom() * sin(sngl->get_the0()));
      newElectron.SetReconPT(embed->get_momG() * sin(embed->get_the0G()));
      newElectron.SetTrkId(itrk);
      newElectron.SetTrkQuality(trk->get_quality(itrk));
      newElectron.SetArm(trk->get_dcarm(itrk));
      newElectron.SetDCSide(trk->get_dcside(itrk));
      newElectron.SetSect(trk->get_sect(itrk));
      newElectron.SetQ(trk->get_charge(itrk));
      newElectron.SetQPrime(trk->get_charge(itrk));
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
      // newElectron.SetDep(trk->get_dep(itrk));
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
  /////////////////////////single sim hits///////////////////////
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
    if (false)
      std::cout << newHit.GetXHit() << " " << newHit.GetYHit() << " " << newHit.GetZHit() << std::endl;
    newHit.SetiLayerFromR();
    if (svxhit->get_layer() != newHit.GetLayer() || svxhit->get_ladder() != newHit.GetLadder() ||
        newHit.GetSensor() != 1)
    {
      std::cout << " smth is wrong " << std::endl;
      continue;
    }
    event->AddVTXHit(&newHit);
  }
  ///////////////embeded+sim hits////////////////
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
    newHit.SetZHit(svxhit->get_xyz_global(2) - event->GetPreciseZ());
    if (false)
      std::cout << newHit.GetXHit() << " " << newHit.GetYHit() << " " << newHit.GetZHit() << std::endl;
    newHit.SetiLayerFromR();
    if (svxhit->get_layer() != newHit.GetLayer() || svxhit->get_ladder() != newHit.GetLadder() ||
        newHit.GetSensor() != 0)
    {
      std::cout << " smth is wrong " << std::endl;
      continue;
    }
    event->AddVTXHit(&newHit);
  }
  ////////////////all hits////////////////////
  if (remove_hadron_hits)
  {
    for (int ihit = 0; ihit < svx->get_nClusters(); ihit++)
    {

      SvxCluster *svxhit = svx->get_Cluster(ihit);

      if (svxhit == nullptr)
      {
        std::cout << "cluster NULL : " << ihit << std::endl;
        continue;
      }

      MyDileptonAnalysis::MyVTXHit newHit;

      newHit.SetClustId(ihit);
      newHit.SetLayer(svxhit->get_layer());
      newHit.SetLadder(svxhit->get_ladder() + 50);
      newHit.SetSensor(2);
      newHit.SetXHit(svxhit->get_xyz_global(0));
      newHit.SetYHit(svxhit->get_xyz_global(1));
      newHit.SetZHit(svxhit->get_xyz_global(2));
      if (false)
        std::cout << newHit.GetXHit() << " " << newHit.GetYHit() << " " << newHit.GetZHit() << std::endl;
      newHit.SetiLayerFromR();
      if (svxhit->get_layer() != newHit.GetLayer() || svxhit->get_ladder() != newHit.GetLadder() - 50 ||
          newHit.GetSensor() != 2)
      {
        std::cout << " smth is wrong " << std::endl;
        continue;
      }
      event->AddVTXHit(&newHit);
    }

    event_container->Associate_Hits_to_Leptons();
    event_container->CleanUpHitList();
  }

  if (fill_TTree)
    event_container->FillTree();

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

int embedana::ReadOrigPartMoms()
{

  int i;

  printf("Small oscar file is %s\n", local_oscarpath.c_str());

  fstream fin(local_oscarpath.c_str());
  string line;

  if (!fin)
  {
    printf("No such small file!!\n");
    return -1;
  }

  i = 0;
  while (getline(fin, line))
  {
    if (line.find("#") == 0)
      continue;
    int i1, i2, i3;
    double px, py, pz, vx, vy, vz, dummy;
    stringstream s(line);
    s >> i1 >> i2;
    if (TMath::Abs(i2) > 1)
    {
      s >> i3 >> px >> py >> pz >> dummy >> dummy >> vx >> vy >> vz;
      // cout<<px<<" "<<py<<" "<<pz<<endl;
      if (i >= 10000)
      {
        printf("Too much particles in small file!\n");
        return -1;
      }

      InData_read[i].id = i2;
      InData_read[i].px = px;
      InData_read[i].py = py;
      InData_read[i].pz = pz;
      InData_read[i].vx = vx;
      InData_read[i].vy = vy;
      InData_read[i].vz = vz;

      ++i;
    }
  }
  cout << "Nlines: " << i << endl;
  if (i != 10000)
  {
    printf("Wrong number of tracks");
    return -1;
  }

  return 0;
}

void embedana::calcDCA_BCbyCircleProjection(
    float pt, float phi, float the,  // pt in xy-plane, phi, theta at inner most layer
    int charge,                      // charge of the track
    float hx, float hy, float hz,    // hit position at inner most layer
    float vx, float vy,              // beam center
    float fieldPolarity,             //
    float *dx, float *dy, float *dz, // dca position
    float *d2dca_bc)                 // return
{
  static const float B = 0.90;
  static const float b = 0.003 * B;

  // vector to rotation center at hit1
  float R = pt / b;
  float pz = pt / tan(the);

  float b_sign = (fieldPolarity > 0) ? -1.0 : 1.0;
  float dir = ((b_sign)*charge > 0.) ? -1.0 : 1.0;
  float cx = hx - dir * R * sin(phi);
  float cy = hy + dir * R * cos(phi);

  // L is a distance btw the rotation center and primary vtx
  // psi is a angle of the vector starting from the rotation center to primary vertex
  float L = sqrt((vx - cx) * (vx - cx) + (vy - cy) * (vy - cy));
  float psi = atan2((vy - cy), (vx - cx));

  // DCA point
  *dx = cx + R * cos(psi);
  *dy = cy + R * sin(psi);

  float dphi = phi - dir * 0.5 * M_PI - psi;
  if (dphi > M_PI * 1.5)
    dphi -= M_PI * 2.;
  else if (dphi < -M_PI * 1.5)
    dphi += M_PI * 2.;
  float dzdphi = dir * pz * R / pt;
  *dz = hz + dphi * dzdphi;

  // DCA value
  *d2dca_bc = b_sign * (charge * (R - L));
}
