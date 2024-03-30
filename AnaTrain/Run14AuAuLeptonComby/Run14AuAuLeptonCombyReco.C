#include <stdint.h>
#include "Run14AuAuLeptonCombyReco.h"


Run14AuAuLeptonCombyReco::Run14AuAuLeptonCombyReco(const char *outfile)
{
    InitParams();

    outfilename = outfile;
    ThisName = "Run14 Au+Au 200 GeV Run14AuAuLeptonComby";
    ncalls = 0;
    npassed = 0;
    verbosity = 0;

    remove_hadron_hits = 0;
    fill_QA_hadron_hists = 0; 
    fill_QA_lepton_hists = 0; 
    fill_TTree = 0;
    fill_d_dphi_hists = 0;
    fill_DCA_hists = 0;
    use_d_dphi_DCA = 0;
    do_track_QA = 0;
    do_reveal_hadron = 0;
    fill_true_DCA = 0;
    check_veto = 0;

    ul = nullptr;
    event_container = nullptr;

    return;
}

Run14AuAuLeptonCombyReco::~Run14AuAuLeptonCombyReco()
{
    delete event_container;
    std::cout<<"Event Container was deleted"<<std::endl;
    return;
}

int Run14AuAuLeptonCombyReco::Init(PHCompositeNode *topNode)
{
    const int ultraLightInitErrorCode = InitUltraLight(topNode);
    if (ultraLightInitErrorCode < 0)
        return -1;

    recoConsts* rc  = recoConsts::instance(); 
    remove_hadron_hits = rc->get_IntFlag("Remove_hadron_hits", 0);
    fill_QA_hadron_hists = rc->get_IntFlag("Fill_QA_hadron_hists", 0);
    fill_QA_lepton_hists = rc->get_IntFlag("Fill_QA_lepton_hists", 0);
    fill_TTree = rc->get_IntFlag("Fill_TTree", 0);
    fill_d_dphi_hists = rc->get_IntFlag("Fill_d_dphi_hists", 0);
    fill_DCA_hists = rc->get_IntFlag("Fill_DCA_hists", 0);
    use_d_dphi_DCA = rc->get_IntFlag("Use_d_dphi_DCA", 0);
    do_track_QA = rc->get_IntFlag("Do_track_QA", 0);
    do_reveal_hadron = rc->get_IntFlag("Do_reveal_hadron", 0);
    fill_true_DCA = rc->get_IntFlag("Fill_true_DCA", 0);
    check_veto = rc->get_IntFlag("Check_Veto", 0);


    std::cout<<"Remove_hadron_hits: "<<remove_hadron_hits<<std::endl;
    std::cout<<"fill_QA_hadron_hists: "<<fill_QA_hadron_hists<<std::endl;
    std::cout<<"fill_QA_lepton_hists: "<<fill_QA_lepton_hists<<std::endl;
    std::cout<<"fill_TTree: "<<fill_TTree<<std::endl;
    std::cout<<"fill_d_dphi_hists: "<<fill_d_dphi_hists<<std::endl;
    std::cout<<"fill_DCA_hists: "<<fill_DCA_hists<<std::endl;
    std::cout<<"use_d_dphi_DCA: "<<use_d_dphi_DCA<<std::endl;
    std::cout<<"Do_track_QA: "<<do_track_QA<<std::endl;
    std::cout<<"do_reveal_hadron: "<<do_reveal_hadron<<std::endl;
    std::cout<<"fill_true_DCA: "<<fill_true_DCA<<std::endl;
    std::cout<<"check_veto: "<<check_veto<<std::endl;
    
    event_container = new MyDileptonAnalysis::MyEventContainer();
    event_container->InitEvent();
    event_container->GetHistsFromFile(GetFilePath());
    event_container->CreateOutFileAndInitHists(outfilename,fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,
                                               fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto);

    return 0;
}

int Run14AuAuLeptonCombyReco::InitRun(PHCompositeNode *TopNode)
{
    if(fill_TTree) event_container->ResetTree();
    return 0;
}

int Run14AuAuLeptonCombyReco::ResetEvent(PHCompositeNode *topNode)
{
    ul->Reset();
    event_container->ClearEvent();
    return 0;
}

int Run14AuAuLeptonCombyReco::process_event(PHCompositeNode *TopNode)
{

    if (ncalls++ % 10000 == 0)
    {
        printf("Ncalls = %dk\n", ncalls / 1000);
    }

    const PHGlobal *globalCNT =
        findNode::getClass<PHGlobal>(TopNode, "PHGlobal");
    const PHCentralTrack *particleCNT =
        findNode::getClass<PHCentralTrack>(TopNode, "PHCentralTrack");
    const RunHeader *runHDR =
        findNode::getClass<RunHeader>(TopNode, "RunHeader");
    const TrigLvl1 *Trig =
        findNode::getClass<TrigLvl1>(TopNode, "TrigLvl1");
    const SvxClusterList *svxhitlist =
        findNode::getClass<SvxClusterList>(TopNode, "SvxClusterList");
    const VtxOut *vtxout =
        findNode::getClass<VtxOut>(TopNode, "VtxOut");

    if (!globalCNT)
        std::cout << "NO GLOBAL!!!!!!!!!!!!!!!\n";
    if (!particleCNT)
        std::cout << "NO CENTRAL TRACK!!!!!!!!!!!!!!!\n";
    if (!runHDR)
        std::cout << "NO RUN HEADER!!!!!!!!!!!!!!!\n";
    if (!Trig)
        std::cout << "NO TRIG!!!!!!!!!!!!!!!\n";
    if (!svxhitlist)
        std::cout << "NO SvxClusterList!!!!!!!!!!!!!!!\n";
    if (!vtxout)
        std::cout << "NO vtxout!!!!!!!!!!!!!!!\n";

    const int run_number = runHDR->get_RunNumber();
    int run_group = get_rungroup(run_number) - 1;
    if(run_group<0||run_group>6) run_group = 0;
    const float bbc_vertex = globalCNT->getBbcZVertex();
    const float centrality = globalCNT->getCentrality();

    const int trigscaled_on = Trig->get_lvl1_trigscaled_bit(TRIGGERBIT);
    if (!trigscaled_on)
        return false;

    // ZDC coincidence
    const float zdc1 = globalCNT->getZdcEnergyN();
    const float zdc2 = globalCNT->getZdcEnergyS();
    const float zdcz = globalCNT->getZdcZVertex();
    const bool isZDCOK = (zdc1 > 0 && zdc2 > 0 && zdcz > -9000);
    if (!isZDCOK)
        return false;

    if (run_number < 0)
        return 0;
    if (fabs(bbc_vertex) > BBC_VERTEX_CUT)
        return 0;

    if (centrality < 0 || centrality > 93)
        return 0;

    npassed++;
    
    const float bbcq = globalCNT->getBbcChargeN() + globalCNT->getBbcChargeS();
    const float bbcT0 = globalCNT->getBbcTimeZero();

    PHPoint precisevtx = vtxout->get_Vertex();

    const float precise_x = precisevtx.getX();
    const float precise_y = precisevtx.getY();
    const float precise_z = precisevtx.getZ();

    MyDileptonAnalysis::MyEvent *event = event_container->GetEvent();

    event->SetCentrality(centrality);
    event->SetPreciseX(precise_x);
    event->SetPreciseY(precise_y);
    event->SetPreciseZ(precise_z);
    event->SetEvtNo(ncalls);
    event->SetRunNumber(run_group);
    event->SetBBCcharge(bbcq);
    event->SetVtxZ(bbc_vertex);
    event->SetBBCtimeN(bbcT0);

    const int run_group_beamoffset = event->GetRunGroup(run_number);
    const int n_tracks = particleCNT->get_npart();
    for (int itrk_reco = 0; itrk_reco < n_tracks; ++itrk_reco)
    {
        const int itype = applySingleTrackCut(particleCNT, itrk_reco, precise_z, centrality, run_number);
        MyDileptonAnalysis::MyElectron newElectron;// = new MyDileptonAnalysis::MyElectron();
        MyDileptonAnalysis::MyHadron newHadron;// = new MyDileptonAnalysis::MyHadron();
        switch (itype)
        {
        case 0:
            if(remove_hadron_hits||do_track_QA){
                set_track(&newHadron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset);
                event->AddHadron(&newHadron);
            }
            break;
        case 1:
            if(remove_hadron_hits||do_track_QA){
                set_track(&newHadron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset);
                event->AddHadron(&newHadron);
            }
            break;
        case 3:
            if(do_reveal_hadron){
                set_track(&newElectron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset);
                event->AddElecCand(&newElectron);
            }
            if(remove_hadron_hits||do_track_QA){
                set_track(&newHadron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset);
                event->AddHadron(&newHadron);
            }
            break;
        case 2:
            set_track(&newElectron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset);
            newElectron.SetChi2(particleCNT->get_chi2(itrk_reco));
            newElectron.SetN0(particleCNT->get_n0(itrk_reco));
            newElectron.SetNPE0(particleCNT->get_npe0(itrk_reco));
            newElectron.SetDISP(particleCNT->get_disp(itrk_reco));
            newElectron.SetEmcdz_e(particleCNT->get_emcsdz_e(itrk_reco));
            newElectron.SetEmcdphi_e(particleCNT->get_emcsdphi_e(itrk_reco));
            if(newElectron.GetPtPrime()>E_PT && newElectron.GetPtPrime()<MAX_PT) event->AddTrack(&newElectron);
            else event->AddElecCand(&newElectron);
            break;
        default:
            continue;
        }
    }

    if(event->GetNtrack()<1) return 0;

    fill_SVXHits_to_myevent(svxhitlist, event);

    event_container->Associate_Hits_to_Leptons(remove_hadron_hits);
    
    int n_electrons = event->GetNtrack()*remove_hadron_hits;
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
        MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
        if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
           (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 ))
        {
            event->RemoveTrackEntry(itrk);
            n_electrons--;
            itrk--;
            continue;
        }
        mytrk->ZeroHitCounters();
    }

    if(event->GetNtrack()<1 || (centrality < 20 && event->GetNtrack() < 1 ) ) return 0;
    
    if(remove_hadron_hits) 
    {
        event_container->Associate_Hits_to_Hadrons();
        event_container->Associate_Hits_to_Leptons();
    }
    
    if(check_veto) event_container->CheckVeto();
    if(use_d_dphi_DCA)  event_container->FillDphiHists();
    if(do_reveal_hadron) event_container->Reveal_Hadron();
    if(fill_TTree) event_container->FillTree();

    for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
    {
        MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
        if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
           (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 )) continue;
           
        int addit_reject = 0;
        if( ((mytrk->GetMinsDphi(0) > 0) || mytrk->GetPtPrime()>0.6) && (mytrk->GetMinsDphi(0) > -1 || mytrk->GetPtPrime()>0.9) ) addit_reject = 1;
        if ( mytrk->GetMinsDphi(0) > 0)
             addit_reject += 10;
        
        int hadron_reject = 0;
        if ( mytrk->GetPtPrime() > 0.4 ) hadron_reject=10;
        if ( fabs(mytrk->GetEmcdphi_e())<2 && fabs(mytrk->GetEmcdz_e())<2 && 
        ((mytrk->GetEcore()/mytrk->GetPtot()>0.8 && mytrk->GetDep()<2 ) ||  mytrk->GetPtPrime()>0.9) )  hadron_reject+=1;
        
        const int ptype = 1 + (1 - mytrk->GetChargePrime()) / 2; //temporary changed to GetGharge cuase in fact its prime

        UltraLightTrack ult("Run14AuAuLeptonComby Track",
                            Run14AuAuLeptonCombyEnum::LAST_DOUBLE,
                            Run14AuAuLeptonCombyEnum::LAST_INTEGER);

        ult.set_px(mytrk->GetPx());
        ult.set_py(mytrk->GetPy());
        ult.set_pz(mytrk->GetPz());

        ult.set_double(Run14AuAuLeptonCombyEnum::ALPHA, mytrk->GetAlphaPrime());
        ult.set_double(Run14AuAuLeptonCombyEnum::PHI, mytrk->GetPhiDC());
        ult.set_double(Run14AuAuLeptonCombyEnum::ZED, mytrk->GetZDC());
        ult.set_double(Run14AuAuLeptonCombyEnum::CRKPHI, mytrk->GetCrkphi());
        ult.set_double(Run14AuAuLeptonCombyEnum::CRKZED, mytrk->GetCrkz());
        //ult.set_double(Run14AuAuLeptonCombyEnum::DCA,  mytrk->GetDCA());
        //ult.set_double(Run14AuAuLeptonCombyEnum::SDCA, mytrk->GetsDCA());
        ult.set_double(Run14AuAuLeptonCombyEnum::PHI0, mytrk->GetPhi0Prime());
        ult.set_double(Run14AuAuLeptonCombyEnum::KEFF,  mytrk->GetReconPT()/mytrk->GetPtPrime());
        
        ult.set_double(Run14AuAuLeptonCombyEnum::DCAX,  mytrk->GetDCAX2());
        ult.set_double(Run14AuAuLeptonCombyEnum::DCAY,  mytrk->GetDCAY2());

        ult.set_integer(Run14AuAuLeptonCombyEnum::PTYPE, ptype);
        ult.set_integer(Run14AuAuLeptonCombyEnum::MATCH, addit_reject);
        ult.set_integer(Run14AuAuLeptonCombyEnum::GHOST, hadron_reject);

        ul->AddTrack(&ult);
    }

    return 0;
}

int Run14AuAuLeptonCombyReco::End(PHCompositeNode *topNode)
{ 
    event_container->WriteOutFile();
    std::cout << "The routine was called " << ncalls << " times." << std::endl;
    std::cout << "There were " << npassed << " events that passed triggers" << std::endl;
    return 0;
}

int Run14AuAuLeptonCombyReco::applySingleTrackCut(const PHCentralTrack *d_trk, const int itrk, const float vertex, const float centrality, const int run_number)
{
    const float px = d_trk->get_px(itrk);
    const float py = d_trk->get_py(itrk);
    const float pz = d_trk->get_pz(itrk);
    const float ecore = d_trk->get_ecore(itrk);
    const float p = sqrt(px * px + py * py + pz * pz);
    const float pT = sqrt(px * px + py * py);
    const float rich_phi = d_trk->get_center_phi(itrk);
    
    const float dep = d_trk->get_dep(itrk);

    if (Z_GLOBAL > -99 && fabs(d_trk->get_zed(itrk)) >= Z_GLOBAL)
        return -1;

    if( ( pT<0.09 || d_trk->get_quality(itrk)<7 )) return -1;
    //if(IsCentralSupportCut(d_trk->get_the0(itrk), vertex)) 
    //    return -1;

    if (DC_DEADMAP)
    {
        const int run_grp = get_rungroup(run_number);
        const int side = d_trk->get_dcside(itrk);
        const int arm = d_trk->get_dcarm(itrk);
        const float board = get_board(d_trk->get_phi(itrk));
        const float alpha = d_trk->get_alpha(itrk);
        for (int iarea = 0; iarea < MAX_DEAD_AREA; ++iarea)
        {
            if (dead_region(board, alpha, dcmap_xx1[run_grp][side][arm][iarea], dcmap_yy1[run_grp][side][arm][iarea], dcmap_xx2[run_grp][side][arm][iarea], dcmap_yy2[run_grp][side][arm][iarea], dcmap_xx3[run_grp][side][arm][iarea], dcmap_yy3[run_grp][side][arm][iarea], dcmap_xx4[run_grp][side][arm][iarea], dcmap_yy4[run_grp][side][arm][iarea]))
            {
                return -1;
            }
        }
    }
    
    bool pass_quality0 = true, pass_quality1 = true, pass_quality2 = true;
    if (QUALITY[0] > -99 && d_trk->get_quality(itrk) != QUALITY[0])
        pass_quality0 = false;
    if (QUALITY[1] > -99 && d_trk->get_quality(itrk) != QUALITY[1])
        pass_quality1 = false;
    if (QUALITY[2] > -99 && d_trk->get_quality(itrk) != QUALITY[2])
        pass_quality2 = false;
    if (!pass_quality0 && !pass_quality1 && !pass_quality2 && rich_phi<-99)
        return 0;
    
    if (E_PT > -99 && pT <= E_PT*0.9 &&  rich_phi<-99)
        return 0;

    if (MAX_PT > -99 && pT >= MAX_PT*1.05 &&  rich_phi<-99)
        return 0;
        
    if (dep < -2 && ecore < 0.3 && d_trk->get_n0(itrk)<0 && rich_phi<-99)
        return 1;

    if(rich_phi>-99 && ((E_PT > -99 && pT <= E_PT*0.8) || (!pass_quality0 && !pass_quality1 && !pass_quality2) ||
        (MAX_PT > -99 && pT >= MAX_PT*1.2) )) return 3;
    //if (TEMC>-99 && d_trk->get_temc(itrk) > TEMC)
    //    return 1;
    if (N0 > -99 && d_trk->get_n0(itrk) <= N0)
        return 3;
    if (DISP > -99 && d_trk->get_disp(itrk) >= DISP)
        return 3;
    if (CHI2_NPE0 > -99 && (d_trk->get_chi2(itrk) / d_trk->get_npe0(itrk)) >= CHI2_NPE0)
        return 3;
    if (EMCDPHI > -99 && fabs(d_trk->get_emcdphi(itrk)) >= EMCDPHI)
        return 3;
    if (EMCDZ > -99 && fabs(d_trk->get_emcdz(itrk)) >= EMCDZ)
        return 3;
    if (EOVERP > -99 && ecore / p <= EOVERP)
        return 3;
    if (DEP[0] > -99 && dep <= DEP[0])
        return 3;
    if (DEP[1] > -99 && dep >= DEP[1])
        return 3;
    if(centrality<60)
    {
        if(pT > 0.40)
        {
            if (EMCSDPHI > -99 && fabs(d_trk->get_emcsdphi_e(itrk)) >= EMCSDPHI)
                return 3;
            if (EMCSDZ > -99 && fabs(d_trk->get_emcsdz_e(itrk)) >= EMCSDZ)
                return 3;
        }else{
            if (fabs(d_trk->get_emcdphi(itrk)) >= 0.02)
                return 3;
            if (fabs(d_trk->get_emcdz(itrk)-0.8) >= 8)
                return 3;
        }
    }
    if (PROB > -99 && d_trk->get_prob(itrk) <= PROB)
        return 3;
    return 2;
}

template <typename track>
void Run14AuAuLeptonCombyReco::set_track(track *newTrack, const PHCentralTrack *trk, int itrk_reco, const float bbcz, const float svxz, const int rg_beamoffset)
{
    newTrack->SetTrkId(itrk_reco);
    newTrack->SetTrkQuality(trk->get_quality(itrk_reco));
    newTrack->SetArm(trk->get_dcarm(itrk_reco));
    newTrack->SetDCSide(trk->get_dcside(itrk_reco));
    newTrack->SetSect(trk->get_sect(itrk_reco));

    newTrack->SetPt(sqrt(trk->get_px(itrk_reco) * trk->get_px(itrk_reco) + trk->get_py(itrk_reco) * trk->get_py(itrk_reco)));
    newTrack->SetQ(trk->get_charge(itrk_reco));

    newTrack->SetPhiDC(trk->get_phi(itrk_reco));
    newTrack->SetPhi0(trk->get_phi0(itrk_reco));
    newTrack->SetThe0(trk->get_the0(itrk_reco));
    newTrack->SetZDC(trk->get_zed(itrk_reco));
    newTrack->SetAlpha(trk->get_alpha(itrk_reco));
    newTrack->SetEmcId(trk->get_emcid(itrk_reco));
    newTrack->SetEcore(trk->get_ecore(itrk_reco));
    newTrack->SetDep(trk->get_dep(itrk_reco));
    newTrack->SetProb(trk->get_prob(itrk_reco));
    newTrack->SetEmcdz(trk->get_emcdz(itrk_reco));
    newTrack->SetEmcdphi(trk->get_emcdphi(itrk_reco));
    newTrack->SetEmcTower(trk->get_sect(itrk_reco), trk->get_ysect(itrk_reco), trk->get_zsect(itrk_reco));
    newTrack->SetTOFE(trk->get_ttof(itrk_reco));
    newTrack->SetEmcTOF(trk->get_temc(itrk_reco));
    
    newTrack->SetTOFDPHI(trk->get_n0(itrk_reco));
    newTrack->SetTOFDZ(trk->get_disp(itrk_reco));
    newTrack->SetPC3SDPHI(trk->get_pc3sdphi(itrk_reco));
    newTrack->SetPC3SDZ(trk->get_pc3sdz(itrk_reco));

    newTrack->SetCrkphi(trk->get_center_phi(itrk_reco));
    newTrack->SetCrkz(trk->get_center_z(itrk_reco));

    newTrack->SetPrimes(bbcz, svxz, rg_beamoffset);
}

void Run14AuAuLeptonCombyReco::fill_SVXHits_to_myevent(const SvxClusterList *svxhitlist, MyDileptonAnalysis::MyEvent *event)
{

    for (int ihit = 0; ihit < svxhitlist->get_nClusters(); ihit++)
    {

        SvxCluster *svxhit = svxhitlist->get_Cluster(ihit);

        if (svxhit == nullptr)
        {
            std::cout << "cluster NULL : " << ihit << std::endl;
            continue;
        }

        MyDileptonAnalysis::MyVTXHit newHit;

        newHit.SetClustId(ihit);
        newHit.SetLayer(svxhit->get_layer());
        newHit.SetLadder(svxhit->get_ladder());
        newHit.SetSensor(svxhit->get_sensor());
        newHit.SetXHit(svxhit->get_xyz_global(0));
        newHit.SetYHit(svxhit->get_xyz_global(1));
        newHit.SetZHit(svxhit->get_xyz_global(2));
        newHit.SetPolars(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
        newHit.SetiLayerFromR();
        event->AddVTXHit(&newHit);
    }
}

int Run14AuAuLeptonCombyReco::InitUltraLight(PHCompositeNode *topNode)
{
    PHNodeIterator iter(topNode);

    ul = new UltraLight("Ultra Light Run14AuAuLeptonComby");

    PHCompositeNode *dstNode =
        dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
        std::cout << PHWHERE << " Could not find DST node" << std::endl;
        return -1;
    }

    PHIODataNode<UltraLight> *ulNode =
        new PHIODataNode<UltraLight>(ul, "Run14AuAuLeptonComby", "PHObject");
    dstNode->addNode(ulNode);

    return 0;
};


const std::string Run14AuAuLeptonCombyReco::GetFilePath()
{
    TOAD toad("Run14AuAuLeptonComby");

    const std::string loc = toad.location("field_map.root");

    return loc;
}