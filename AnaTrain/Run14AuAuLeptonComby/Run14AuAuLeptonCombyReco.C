#include <stdint.h>
#include "Run14AuAuLeptonCombyReco.h"


Run14AuAuLeptonCombyReco::Run14AuAuLeptonCombyReco(const char *outfile, const char *lookup_file) : 
    SubsysReco("Run14AuAuLeptonCombyReco"), reco(lookup_file)
{
    InitParams();
    MoonWalk();

    outfilename = outfile;
    ThisName = "Run14 Au+Au 200 GeV Run14AuAuLeptonComby";
    ncalls = 0;
    npassed = 0;
    verbosity = 0;

    remove_hadron_hits = 0;
    fill_QA_hadron_hists = 0; 
    fill_QA_lepton_hists = 0;
    fill_ddphi_hadron = 0; 
    fill_TTree = 0;
    fill_d_dphi_hists = 0;
    fill_DCA_hists = 0;
    use_iden = 0;
    do_track_QA = 0;
    do_electron_QA = 0;
    do_reveal_hadron = 0;
    fill_flow_hists = 0;
    fill_true_DCA = 0;
    check_veto = 0;
    fill_inv_mass = 0;

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
    fill_ddphi_hadron = rc->get_IntFlag("Fill_ddphi_hadron",0);
    fill_TTree = rc->get_IntFlag("Fill_TTree", 0);
    fill_d_dphi_hists = rc->get_IntFlag("Fill_d_dphi_hists", 0);
    fill_DCA_hists = rc->get_IntFlag("Fill_DCA_hists", 0);
    use_iden = rc->get_IntFlag("Use_ident", 0);
    do_track_QA = rc->get_IntFlag("Do_track_QA", 0);
    do_electron_QA = rc->get_IntFlag("Do_electron_QA", 0);
    do_reveal_hadron = rc->get_IntFlag("Do_reveal_hadron", 0);
    fill_flow_hists = rc->get_IntFlag("Fill_flow", 0);
    fill_true_DCA = rc->get_IntFlag("Fill_true_DCA", 0);
    check_veto = rc->get_IntFlag("Check_Veto", 0);
    fill_inv_mass = rc->get_IntFlag("fill_inv_mass", 0);


    std::cout<<"Remove_hadron_hits:   "<<remove_hadron_hits<<std::endl;
    std::cout<<"fill_QA_hadron_hists: "<<fill_QA_hadron_hists<<std::endl;
    std::cout<<"fill_QA_lepton_hists: "<<fill_QA_lepton_hists<<std::endl;
    std::cout<<"fill_ddphi_hadron:    "<<fill_ddphi_hadron<<std::endl;
    std::cout<<"fill_TTree:           "<<fill_TTree<<std::endl;
    std::cout<<"fill_d_dphi_hists:    "<<fill_d_dphi_hists<<std::endl;
    std::cout<<"fill_DCA_hists:       "<<fill_DCA_hists<<std::endl;
    std::cout<<"use_iden:             "<<use_iden<<std::endl;
    std::cout<<"Do_track_QA:          "<<do_track_QA<<std::endl;
    std::cout<<"Do_electron_QA:       "<<do_electron_QA<<std::endl;
    std::cout<<"do_reveal_hadron:     "<<do_reveal_hadron<<std::endl;
    std::cout<<"fill_flow_hists:      "<<fill_flow_hists<<std::endl;
    std::cout<<"fill_true_DCA:        "<<fill_true_DCA<<std::endl;
    std::cout<<"check_veto:           "<<check_veto<<std::endl;
    std::cout<<"fill_inv_mass:        "<<fill_inv_mass<<std::endl;
    
    event_container = new MyDileptonAnalysis::MyEventContainer();
    event_container->InitEvent();
    event_container->GetHistsFromFile(GetFilePath());
    event_container->CreateOutFileAndInitHists(outfilename,fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,///temporary to be removed
                                               fill_DCA_hists, do_track_QA+do_electron_QA, fill_flow_hists, fill_true_DCA, check_veto,
                                               fill_flow_hists>0?0:fill_inv_mass);

    return 0;
}

int Run14AuAuLeptonCombyReco::InitRun(PHCompositeNode *TopNode)
{
    if(fill_TTree) event_container->ResetTree();
    InitWalk(TopNode);
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

    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(0);
    bool do_event_selection = false;
    
    PHGlobal *globalCNT =
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
    const emcClusterContainer* emccont =
        findNode::getClass<emcClusterContainer>(TopNode, "emcClusterContainer");
    const ReactionPlaneObject* rpobject = 
        findNode::getClass<ReactionPlaneObject>(TopNode, "ReactionPlaneObject");

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
    if (!emccont)
        std::cout << "NO emcClusterContainer!!!!!!!!!!!!!!!\n";
    if (!rpobject)
        std::cout << "NO ReactionPlaneObject!!!!!!!!!!!!!!!\n";

    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(1);

    const int run_number = runHDR->get_RunNumber();
    const float bbc_vertex = globalCNT->getBbcZVertex();
    const float centrality = globalCNT->getCentrality();
    globalCNT->setBbcZVertex(100.);
    
    if(fill_TTree||fill_true_DCA) 
        if(TMath::Abs(vtxout->get_Vertex().getZ())<10)
            event_container->FillEventHist(8);

    const int trigscaled_on = Trig->get_lvl1_trigscaled_bit(TRIGGERBIT);
    if (!trigscaled_on && do_event_selection)
        return false;

    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(2);

    // ZDC coincidence
    const float zdc1 = globalCNT->getZdcEnergyN();
    const float zdc2 = globalCNT->getZdcEnergyS();
    const float zdcz = globalCNT->getZdcZVertex();
    const bool isZDCOK = (zdc1 > 0 && zdc2 > 0 && zdcz > -9000);
    if (!isZDCOK&&do_event_selection)
        return false;
    if (run_number < 0)
        return 0;
    if (fabs(bbc_vertex) > BBC_VERTEX_CUT && do_event_selection)
        return 0;

    if (centrality < 0 || centrality > 93)
        return 0;

    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(3);

    if(fill_TTree||fill_true_DCA) event_container->FillCentrHist(centrality);

    npassed++;
    
    const float bbcq = globalCNT->getBbcChargeN() + globalCNT->getBbcChargeS();
    const float bbcT0 = globalCNT->getBbcTimeZero();

    PHPoint precisevtx = vtxout->get_Vertex();

    const float precise_x = precisevtx.getX();
    const float precise_y = precisevtx.getY();
    const float precise_z = precisevtx.getZ();

    if (TMath::Abs(precise_z) > BBC_VERTEX_CUT)
        return 0;

    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(4);

    if (precise_x < 0.27 || precise_x > 0.40 || precise_y < -0.02 || precise_y > 0.10)
        return 0;   
        
    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(5);
    
    globalCNT->setBbcZVertex(precise_z);

    MyDileptonAnalysis::MyEvent *event = event_container->GetEvent();

    event->SetCentrality(centrality);
    event->SetPreciseX(precise_x);
    event->SetPreciseY(precise_y);
    event->SetPreciseZ(precise_z);
    event->SetEvtNo(ncalls);
    event->SetRunNumber(run_number);/////neeeds a doctor
    event->SetBBCcharge(bbcq);
    event->SetVtxZ(bbc_vertex);
    event->SetBBCtimeN(bbcT0);

     // BBC sum
    ReactionPlaneSngl *rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 1));
    float psi2_BBC = (rpsngl) ? rpsngl->GetPsi() : -9999;

    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 2));
    float psi3_BBC = (rpsngl) ? rpsngl->GetPsi() : -9999;

    event->SetPsi2BBC(psi2_BBC);
    event->SetPsi3BBC(psi3_BBC);

    // FVTX, all sectors w/ eta>1.0
    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 42, 1));
    float psi2_FVTXA0 = (rpsngl) ? rpsngl->GetPsi() : -9999;

    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 42, 2));
    float psi3_FVTXA0 = (rpsngl) ? rpsngl->GetPsi() : -9999;

    event->SetPsi2FVTXA0(psi2_FVTXA0);
    event->SetPsi3FVTXA0(psi3_FVTXA0);

    rpobject->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 1))->SetPsi(psi2_FVTXA0);

    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 0, 1));
    float psi2_BBCS = (rpsngl) ? rpsngl->GetPsi() : -9999;
    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 1, 1));
    float psi2_BBCN = (rpsngl) ? rpsngl->GetPsi() : -9999;
    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 40, 1));
    float psi2_FVTXS = (rpsngl) ? rpsngl->GetPsi() : -9999;
    rpsngl = rpobject->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 41, 1));
    float psi2_FVTXN = (rpsngl) ? rpsngl->GetPsi() : -9999;

    if((fill_TTree||fill_true_DCA) && (psi2_BBC>-9000 || psi2_FVTXA0 >-9000)) event_container->FillEventHist(6);
    if((fill_TTree||fill_true_DCA)) event_container->FillMyVTXHist((int) (centrality/20),precise_x,precise_y,precise_z);

    const int run_group_beamoffset = event->GetRunGroup(run_number);// what is this???? <8?event->GetRunGroup(run_number):7;
    const int n_tracks = particleCNT->get_npart();
    
    for (int itrk_reco = 0; itrk_reco < n_tracks; ++itrk_reco)
    {
        const int itype = applySingleTrackCut(particleCNT, itrk_reco, precise_z, centrality, run_number);
        MyDileptonAnalysis::MyElectron newElectron;// = new MyDileptonAnalysis::MyElectron();
        MyDileptonAnalysis::MyHadron newHadron;// = new MyDileptonAnalysis::MyHadron();
        switch (itype)
        {
        case 0:
            if(do_track_QA){
                set_track(&newHadron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset, emccont);
                event->AddHadron(&newHadron);
            }
            break;
        case 1:
            if(remove_hadron_hits||do_track_QA||fill_QA_hadron_hists){
                set_track(&newHadron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset, emccont);
                if(remove_hadron_hits||do_track_QA)event->AddHadron(&newHadron);
                if(fill_QA_hadron_hists&&newHadron.GetPtPrime()>1.5&&fabs(newHadron.GetPC3SDPHI())<2&&fabs(newHadron.GetPC3SDZ())<2)
                    event->AddHadron(&newHadron);
                if(fill_ddphi_hadron&&fabs(newHadron.GetPC3SDPHI())<2&&fabs(newHadron.GetPC3SDZ())<2)
                {
                    set_track(&newElectron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset, emccont);
                    newElectron.SetCrkphi(-999);
                    event->AddTrack(&newElectron);
                }
            }
            break;
        case 3:
            if(do_reveal_hadron){
                set_track(&newElectron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset, emccont);
                event->AddElecCand(&newElectron);
            }
            if(do_track_QA){
                set_track(&newHadron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset, emccont);
                event->AddHadron(&newHadron);
            }
            if(True)
            {
                set_track(&newElectron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset,  emccont);
                newElectron.SetChi2(particleCNT->get_chi2(itrk_reco));
                newElectron.SetN0(particleCNT->get_n0(itrk_reco));
                newElectron.SetNPE0(particleCNT->get_npe0(itrk_reco));
                newElectron.SetDISP(particleCNT->get_disp(itrk_reco));
                newElectron.SetEmcdz_e(particleCNT->get_emcsdz_e(itrk_reco));
                newElectron.SetEmcdphi_e(particleCNT->get_emcsdphi_e(itrk_reco));
                if(newElectron.GetPtPrime()>E_PT && newElectron.GetPtPrime()<MAX_PT) event->AddTrack(&newElectron);
                else if(do_reveal_hadron) event->AddElecCand(&newElectron);       
            }
            break;
        case 2:
            set_track(&newElectron, particleCNT, itrk_reco, bbc_vertex, precise_z, run_group_beamoffset,  emccont);
            newElectron.SetChi2(particleCNT->get_chi2(itrk_reco));
            newElectron.SetN0(particleCNT->get_n0(itrk_reco));
            newElectron.SetNPE0(particleCNT->get_npe0(itrk_reco));
            newElectron.SetDISP(particleCNT->get_disp(itrk_reco));
            newElectron.SetEmcdz_e(particleCNT->get_emcsdz_e(itrk_reco));
            newElectron.SetEmcdphi_e(particleCNT->get_emcsdphi_e(itrk_reco));
            if(newElectron.GetPtPrime()>E_PT && newElectron.GetPtPrime()<MAX_PT) event->AddTrack(&newElectron);
            else if(do_reveal_hadron) event->AddElecCand(&newElectron);
            break;
        default:
            continue;
        }
    }

    if(event->GetNtrack()<1) return 0;

    if(use_iden)
    {
        Walking(TopNode);
        for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            if(mytrk->GetEmcId()>=0) 
            {
                emcClusterContent* emc = emccont->getCluster(mytrk->GetEmcId());
                if(!emc) continue;
                
	            //mytrk->SetEmcTOF(emc->tofcorr());  
                const int icnttrack = mytrk->GetTrkId();
                mytrk->SetEmcTOF( particleCNT->get_mom(icnttrack)*particleCNT->get_mom(icnttrack)*(emc->tofcorr()*emc->tofcorr()*900./particleCNT->get_plemc(icnttrack)/particleCNT->get_plemc(icnttrack)-1));
            } /// 100*mytrk->GetPtot()*mytrk->GetPtot()*(emc->tofcorr()*emc->tofcorr()*900/mytrk->GetTOFDZ()/mytrk->GetTOFDZ()-1)
        }
        int n_hadrons = event->GetNhadron();
        for (int itrk = 0; itrk < n_hadrons; itrk++)
        {
            MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
            if(mytrk->GetEmcId()>=0) 
            {
                emcClusterContent* emc = emccont->getCluster(mytrk->GetEmcId());
                if(!emc) continue;
	            mytrk->SetEmcTOF(emc->tofcorr());  
            } 
        }
    }
    
    if(fill_ddphi_hadron)
    {
        fill_SVXHits_to_myevent(svxhitlist, event);
        event_container->Associate_Hits_to_Leptons(4.,2.,4.,-fill_ddphi_hadron);
    }

    int n_electrons = event->GetNtrack();
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron mytrk = *event->GetEntry(itrk);
      
      bool skip = false; 
        if (mytrk.GetPtPrime()>MAX_PT || mytrk.GetPtPrime() < E_PT)
          skip = true;
        if (mytrk.GetCrkphi()<-99)
          skip = true;
        if (mytrk.GetEcore()/mytrk.GetPtot() < 0.6)
          skip = true;
        if (mytrk.GetN0() < 0 )
          skip = true;
        if (fabs(mytrk.GetEmcdphi())>0.05 || fabs(mytrk.GetEmcdz())>25 )
          skip = true;

        if( skip ){
          event->RemoveTrackEntry(itrk);
          //event->AddElecCand(&mytrk);
          n_electrons--;
          itrk--;
          continue;
        }
    }
    
    event_container->IdenElectrons();

    n_electrons = event->GetNtrack();
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron mytrk = *event->GetEntry(itrk);
     
      if ( mytrk.GetMcId()<1000 || mytrk.GetPtPrime() < 0.4 || (event->GetCentrality()<20 && mytrk.GetMcId()<10000) ) //adding regualr electron cuts
      {
          event->RemoveTrackEntry(itrk);
          //event->AddElecCand(&mytrk);
          n_electrons--;
          itrk--;
          continue;
      }
    }

    if(event->GetNtrack()<1) return 0;
    if(do_electron_QA) event_container->FillQAHistPreAssoc();

    if(!fill_ddphi_hadron)
    {
        fill_SVXHits_to_myevent(svxhitlist, event);
        event_container->Associate_Hits_to_Leptons(5.,5.,5,0,2,3.0);
        event_container->Associate_Hits_to_Leptons(5.,5.,5,0,1,3.0);
        event_container->Associate_Hits_to_Leptons(5.,5.,5,0,1,3.0);
        if(event_container->GetNGoodElectrons()>1) event_container->Associate_Hits_to_Leptons(5.,5.,5,0,0,3.);
        
        if(false)
        {
            for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
            {
                MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
                if(event_container->GetNGoodElectrons()>0&&mytrk->GetHitCounter(0)>0) std::cout<< "1: " <<mytrk->GetNHits()<<" "<<mytrk->GetTOFDPHI()<<" "<<mytrk->GetGhost()<<std::endl;
            }
            for (int ibdtrack = 0; ibdtrack < (int) event_container->GetNBDThit(); ibdtrack++)
            {       
                if(event_container->GetNGoodElectrons()<1) continue;
                for (int ihit = 0; ihit < event->GetNVTXhit(); ihit++)
                {
                  MyDileptonAnalysis::MyVTXHit myhit = *event->GetVTXHitEntry(ihit);
                if(myhit.GetLadder()>24&&myhit.GetLadder()<50)std::cout<<myhit.GetLadder()<<" "<<myhit.GetClustId()<<std::endl;
                }
                std::cout<<event_container->GetBDTHitEntry(ibdtrack)->GetPt()<<" "<<event_container->GetBDTHitEntry(ibdtrack)->GetEcore()<<" "<<event_container->GetBDTHitEntry(ibdtrack)->GetNBDThit()<<std::endl;
                if(event_container->GetBDTHitEntry(ibdtrack)->GetPt()<0.4) continue;
                MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(ibdtrack);
                std::cout<<mytrk->GetHitIndex(0)<<" "<<mytrk->GetHitIndex(1)<<" "<<mytrk->GetHitIndex(2)<<" "<<mytrk->GetHitIndex(3)<<std::endl;
                for (int jbdthit = 0; jbdthit < (int) event_container->GetBDTHitEntry(ibdtrack)->GetNBDThit(); jbdthit++)
                {
                    std::cout<<event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(0,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(1,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(2,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(3,1)<<" "<<std::endl;
                }
            }
        }
    }
    const int n_good_el = event_container->GetNGoodElectrons();
    if( n_good_el<1 && fill_inv_mass ) return 0;

    if(fill_TTree||fill_true_DCA) event_container->FillEventHist(7);

    if(n_good_el>=0 && (remove_hadron_hits||fill_QA_hadron_hists))
    {
        event_container->Associate_Hits_to_Hadrons();
    }

    n_electrons = event->GetNtrack();
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron mytrk = *event->GetEntry(itrk);
      bool do_reshuf = false;

      if (mytrk.GetHitCounter(0) < 1 || mytrk.GetHitCounter(1) < 1 ||
          (mytrk.GetHitCounter(2) < 1 && mytrk.GetHitCounter(3) < 1))
          do_reshuf = true;

      if (do_reshuf&&fill_inv_mass)
      {
          event->RemoveTrackEntry(itrk);
          //event->AddElecCand(&mytrk);
          n_electrons--;
          itrk--;
          continue;
      }
    }
    
    if(do_electron_QA) event_container->FillQAHist();
    
    n_electrons = event->GetNtrack()*0;
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron mytrk = *event->GetEntry(itrk);
      bool do_reshuf = false;

      if (mytrk.GetTOFDPHI()<100 || mytrk.GetNHits()<1000)
          do_reshuf = true;

      if (do_reshuf&&fill_inv_mass)
      {
          event->RemoveTrackEntry(itrk);
          //event->AddElecCand(&mytrk);
          n_electrons--;
          itrk--;
          continue;
      }
    }

    if(use_iden)
    {
        int n_electrons = event->GetNtrack();
        for (int itrk = 0; itrk < n_electrons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            if(mytrk->GetChargePrime()<0) continue;
            for (int jtrk = 0; jtrk < n_electrons; jtrk++)
            {
                MyDileptonAnalysis::MyElectron *mytrk2 = event->GetEntry(jtrk);
                if(mytrk2->GetChargePrime()>0) continue;
                const int solut =  Solution(mytrk,mytrk2,&reco,precise_z);
                if(solut>0)
                {
                    mytrk ->SetIsConv(solut);
                    mytrk2->SetIsConv(solut);
                }
            }
        }
    }

    //const int n_ghosts = event_container->isGhostEvent();
    //if(n_ghosts>0) std::cout<<"GHOST EVENT with centrality and n_ghosts: " <<event->GetCentrality()<<" "<<n_ghosts<<std::endl;  //removing ghost
    
    if(event_container->GetNGoodElectrons()>1) event_container->ResetRecoverFGVars();//needs check for hadrons
    if(check_veto) event_container->CheckVeto();
    
    //event->ReshuffleElectrons();
    if(fill_TTree) event_container->CleanUpHitList();
    if(fill_d_dphi_hists)  event_container->FillDphiHists();
    if(fill_true_DCA) event_container->FillTrueDCA();
    if(fill_flow_hists) event_container->FillFlow(psi2_BBCS, psi2_BBCN, psi2_FVTXS, psi2_FVTXN);
    if(do_reveal_hadron) event_container->Reveal_Hadron();
    if(fill_TTree) event_container->FillTree();
    if(fill_inv_mass&&!fill_QA_lepton_hists&&!fill_flow_hists) event_container->fill_inv_mass();

    for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
    {
        MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
        
        if ( mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
           ( mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 )) continue;
           
        if ( mytrk->GetPtPrime() < 0.4) continue;///add wenquing cut
        if ( mytrk->GetPtPrime() > 4.4) continue;///temporary cut
        //if (mytrk->GetMcId() < 1000) continue;// mytrk->GetMcId()%10<5
        //if (mytrk->GetNHits() < 10 || mytrk->GetTOFDPHI()<10) continue; 
        //if ( event->GetCentrality() < 20 && mytrk->GetN0() < 2 + SQR(mytrk->GetDisp()) / 8. ) continue;
        //if (event->GetCentrality() < 20  && !(mytrk->GetN0()>1 && mytrk->GetEcore()/mytrk->GetPtot()>0.8 &&
        //    mytrk->GetDisp()<4 && mytrk->GetChi2()/mytrk->GetNpe0()<10 && mytrk->GetProb()>0.01 )) continue;

        int hadron_reject = mytrk->GetMcId();
        if ( (mytrk->GetEmcTOF() > - 1 && mytrk->GetEmcTOF() < 0.4 && mytrk->GetTOFE() < -100) 
            || mytrk->GetMcId()%10 > 5 || TMath::Abs(mytrk->GetTOFE()*0.01) < 0.6 ) hadron_reject+=10;
        if ( mytrk->GetProb()>0.1  ) hadron_reject+=50;
        int hit_assocaition = mytrk->GetNHits();
        if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<3 && TMath::Abs(mytrk->GetMinsDthe(3))<3) ||
               (TMath::Abs(mytrk->GetMinsDphi(2))<3 && TMath::Abs(mytrk->GetMinsDthe(2))<3) ) && 
               (TMath::Abs(mytrk->GetMinsDthe(1))<3 && TMath::Abs(mytrk->GetMinsDphi(1))<3) && 
               (TMath::Abs(mytrk->GetMinsDthe(0))<3 && TMath::Abs(mytrk->GetMinsDphi(0))<3) && mytrk->GetMinsDphi(0)>-2) ) hit_assocaition+=1;
        if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<3 && TMath::Abs(mytrk->GetMinsDthe(3))<3) ||
               (TMath::Abs(mytrk->GetMinsDphi(2))<3 && TMath::Abs(mytrk->GetMinsDthe(2))<3) ) && 
               (TMath::Abs(mytrk->GetMinsDthe(1))<3 && TMath::Abs(mytrk->GetMinsDphi(1))<3) && 
               (TMath::Abs(mytrk->GetMinsDthe(0))<3 && TMath::Abs(mytrk->GetMinsDphi(0))<3)  && mytrk->GetMinsDphi(0)> 0 )) hit_assocaition+=5;
        //if(mytrk->GetNHits()>90) addit_reject = 1;
        //if(mytrk->GetTOFDPHI()>900 && mytrk->GetGhost()<1) addit_reject = 1;
        ////if(mytrk->GetGhost()!=-9 && mytrk->GetNHits()>900 && mytrk->GetTOFDPHI()>90) addit_reject = 10;
        //if(mytrk->GetNHits()>900 && mytrk->GetTOFDPHI()>900 && mytrk->GetGhost()<1) addit_reject = 10; 

        MyDileptonAnalysis::MyVTXHit *vtxhit0 = event->GetVTXHitEntry(mytrk->GetHitIndex(0));
        MyDileptonAnalysis::MyVTXHit *vtxhit1 = event->GetVTXHitEntry(mytrk->GetHitIndex(1));
        MyDileptonAnalysis::MyVTXHit *vtxhit2 = nullptr; 
        if (mytrk->GetHitCounter(2)>0) vtxhit2 = event->GetVTXHitEntry(mytrk->GetHitIndex(2));
        else                           vtxhit2 = event->GetVTXHitEntry(mytrk->GetHitIndex(3));
        
        int conv_reject = mytrk->GetTOFDPHI();//mytrk->GetMcId();
        //if ( mytrk->GetMcId() > 900  && ( ( mytrk->GetEmcTOF() > - 1 && mytrk->GetEmcTOF() < 0.4 && mytrk->GetTOFE() < -100) 
        //    || mytrk->GetMcId()%10 > 5 || TMath::Abs(mytrk->GetTOFE()*0.01) < 0.6) ) hadron_reject = 1000; //// not sure it works after pt=0.7
        //if (hadron_reject<100) continue;
        //if ( (  mytrk->GetMcId()%10 > 5 && mytrk->GetProb()>0.1  )
        //     || TMath::Abs(mytrk->GetTOFE()*0.01) < 0.4 ) hadron_reject = 1000;
        //std::cout<<event->GetCentrality()<<" "<<mytrk->GetPtPrime()*mytrk->GetChargePrime()<<" "<<mytrk->GetIsConv()<<" "<<hadron_reject<<" "<<hit_assocaition<<" "<<conv_reject<<std::endl;
        //const int local_charge = mytrk->GetPhiDC() -  vtxhit0->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()) > 0 ? 1 : -1;
        //if( mytrk->GetChargePrime() != mytrk->GetCharge() || mytrk->GetChargePrime() != local_charge || mytrk->GetCharge() != local_charge )
        //    hadron_reject = 100;///to be removed
        //std::cout<<event->GetCentrality()<<" "<<mytrk->GetPtPrime()*mytrk->GetChargePrime()<<" "<<mytrk->GetPhi0()<<" "<<mytrk->GetPhi0Prime()<< " " << mytrk->GetThe0()<< " " << mytrk->GetThe0Prime()<<std::endl;

        const int ptype = 1 + (1 - mytrk->GetChargePrime()) / 2;

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
        
        ult.set_double(Run14AuAuLeptonCombyEnum::DCAX,  mytrk->GetDCAX2());
        ult.set_double(Run14AuAuLeptonCombyEnum::DCAY,  mytrk->GetDCAY2());

        ult.set_double(Run14AuAuLeptonCombyEnum::ZVTX, event->GetPreciseZ());
        ult.set_double(Run14AuAuLeptonCombyEnum::PSI,  event->GetPsi2FVTXA0());

        ult.set_double(Run14AuAuLeptonCombyEnum::PHI1, vtxhit0->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()));
        ult.set_double(Run14AuAuLeptonCombyEnum::PHI2, vtxhit1->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()));
        ult.set_double(Run14AuAuLeptonCombyEnum::PHI3, vtxhit2->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()));
        ult.set_double(Run14AuAuLeptonCombyEnum::THE1, vtxhit0->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()));
        ult.set_double(Run14AuAuLeptonCombyEnum::THE2, vtxhit1->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()));
        ult.set_double(Run14AuAuLeptonCombyEnum::THE3, vtxhit2->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ()));

        ult.set_integer(Run14AuAuLeptonCombyEnum::PTYPE, ptype);
        ult.set_integer(Run14AuAuLeptonCombyEnum::HADRON_REJECT, hadron_reject);
        ult.set_integer(Run14AuAuLeptonCombyEnum::HIT_ASSOC, hit_assocaition);
        ult.set_integer(Run14AuAuLeptonCombyEnum::CONV_REJECT, conv_reject);
        ult.set_integer(Run14AuAuLeptonCombyEnum::ID1, mytrk->GetHitIndex(0));
        ult.set_integer(Run14AuAuLeptonCombyEnum::ID2, mytrk->GetHitIndex(1));
        ult.set_integer(Run14AuAuLeptonCombyEnum::ID3, mytrk->GetHitCounter(2)>0? mytrk->GetHitIndex(2):mytrk->GetHitIndex(3));
        ul->AddTrack(&ult);
    }

    return 0;
}

int Run14AuAuLeptonCombyReco::End(PHCompositeNode *topNode)
{ 
    event_container->WriteOutFile();
    StopWalking();
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
    
    if (E_PT > -99 && pT <= E_PT*0.5 &&  rich_phi<-99)
        return 0;

    if (MAX_PT > -99 && pT >= MAX_PT*2 &&  rich_phi<-99)
        return 0;
        
    if (ecore < 0.2 && ecore/p < 0.5 && d_trk->get_n0(itrk)<0 /*&& rich_phi<-99*/)
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
void Run14AuAuLeptonCombyReco::set_track(track *newTrack, const PHCentralTrack *trk, int itrk_reco, const float bbcz, const float svxz, const int rg_beamoffset, const emcClusterContainer* emccont)
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
    
    newTrack->SetTOFDPHI(trk->get_tofdphi(itrk_reco));
    newTrack->SetTOFDZ(trk->get_tofdz(itrk_reco));
    newTrack->SetPC3SDPHI(trk->get_pc3sdphi(itrk_reco));
    newTrack->SetPC3SDZ(trk->get_pc3sdz(itrk_reco));

    newTrack->SetCrkphi(trk->get_center_phi(itrk_reco));
    newTrack->SetCrkz(trk->get_center_z(itrk_reco));

    newTrack->SetPrimes(bbcz, svxz, rg_beamoffset);
    newTrack->ResetPrimes(bbcz, svxz, rg_beamoffset);
    newTrack->SetTOFE((trk->get_m2tof(itrk_reco))*100);
    newTrack->SetEmcTOF(trk->get_temc(itrk_reco));
	if(trk->get_emcid(itrk_reco) >= 0)
    {
	    emcClusterContent* emc = emccont->getCluster(trk->get_emcid(itrk_reco));
        if(!emc) return;
	    newTrack->SetEmcTOF( emc->tofcorr() ); 
        //newTrack->SetEmcTOF( 100*trk->get_mom(itrk_reco)*trk->get_mom(itrk_reco)*(emc->tofcorr()*emc->tofcorr()*900/trk->get_plemc(itrk_reco)/trk->get_plemc(itrk_reco)-1));
        //newTrack->SetEmcTOF( 1);
	    //newTrack.SetE( emc->e());
	    //newTrack->SetEcore(emc->ecore());
    //std::cout<<trk->get_temc(itrk_reco)<<" "<<emc->tof()<<" "<<emc->tofcorr()<<" "<<trk->get_ecore(itrk_reco)<<" "<<emc->e()<<std::endl;
	}
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
        newHit.SetLadder(svxhit->get_ladder()+50);
        newHit.SetSensor(svxhit->get_sensor());
        newHit.SetXHit(svxhit->get_xyz_global(0));
        newHit.SetYHit(svxhit->get_xyz_global(1));
        newHit.SetZHit(svxhit->get_xyz_global(2));
        newHit.SetiLayerFromR();
        if( svxhit->get_layer()!=newHit.GetLayer()||svxhit->get_ladder()!=newHit.GetLadder()-50||
        svxhit->get_sensor()!=newHit.GetSensor()) 
        {
            std::cout<<" smth is wrong "<<std::endl;
            continue;
        }
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

static double getPhiv(double px_e, double py_e, double pz_e, double px_p, double py_p, double pz_p)
{
    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector pair;

    p1.SetX(px_e);
    p1.SetY(py_e);
    p1.SetZ(pz_e);
    p1.SetE(sqrt(pow(p1.P(), 2) + SQR(0.0005)));
    p2.SetX(px_p);
    p2.SetY(py_p);
    p2.SetZ(pz_p);
    p2.SetE(sqrt(pow(p2.P(), 2) + SQR(0.0005)));
    pair = p1 + p2; // pair corresponds to photon if the pair matches

    TVector3 P1, P2, Photon, z, v, u, w, wc;

    z.SetX(0);
    z.SetY(0);
    z.SetZ(1); // unit vector along z

    P1.SetX(px_e);
    P1.SetY(py_e);
    P1.SetZ(pz_e);

    P2.SetX(px_p);
    P2.SetY(py_p);
    P2.SetZ(pz_p);

    Photon.SetX(pair.Px());
    Photon.SetY(pair.Py());
    Photon.SetZ(pair.Pz());

    v = (P1.Cross(P2)).Unit(); // unit vector v corresponds to the unit vector of the cross product of ep pair
    u = Photon.Unit();
    w = (u.Cross(v)).Unit();
    wc = (u.Cross(z)).Unit();

    double phiv = acos(-w.Dot(wc));
    return phiv;
}


int Run14AuAuLeptonCombyReco::Solution(MyDileptonAnalysis::MyTrack *mytrk1, MyDileptonAnalysis::MyTrack *mytrk2, MyDileptonAnalysis::Reconstruction *reco, float zVtx)
{ // only solution cut

    const float DPHI = 0.05, R_LO = 1, R_HI = 29, R_HI2 = 23, R_HI3 = 20, DTHETA=0.01, ZEDCUT=4, PHIVCUT=0.1, ZEDCUT2=2, DTHETA2=0.006, DPHI2 = 0.005;

    MyDileptonAnalysis::MyPair mypair;

    reco->findIntersection(mytrk1, mytrk2, &mypair, zVtx);
    const float radius_r = mypair.GetRPair();

    r_conv_hist->Fill(radius_r,0.5*(mytrk1->GetPt()+mytrk2->GetPt()));
    
    if ((radius_r <= R_LO) || (radius_r >= R_HI))
        return 0;

    float dphi_r = mypair.GetPhiElectron() - mypair.GetPhiPositron();
    if (dphi_r > TMath::Pi())
        dphi_r = 2 * TMath::Pi() - dphi_r; // 1.3pi->0.7pi
    if (dphi_r < -TMath::Pi())
        dphi_r = -2 * TMath::Pi() - dphi_r; // -1.3pi->-0.7pi

    const float phiv = getPhiv(mytrk2->GetPx(), mytrk2->GetPy(), mytrk2->GetPz(), mytrk1->GetPx(), mytrk1->GetPy(), mytrk1->GetPz());
    const float dzed =  mytrk2->GetZDC() - mytrk1->GetZDC();
    const float dtheta_r = mypair.GetThetaElectron() - mypair.GetThetaPositron();

    mytrk1->SetRConv(radius_r); mytrk2->SetRConv(radius_r);
    mytrk1->SetPhiConv(mypair.GetPhiPositron()); mytrk2->SetPhiConv(mypair.GetPhiElectron());
    mytrk1->SetTheConv(mypair.GetThetaPositron()); mytrk2->SetTheConv(mypair.GetThetaElectron());
    mytrk1->SetdZedConv(dzed); mytrk2->SetdZedConv(dzed);
    mytrk1->SetPhiVConv(phiv); mytrk2->SetPhiVConv(phiv);
    
    if (TMath::Abs(dzed) < ZEDCUT) phi_conv_hist->Fill(dphi_r,0.5*(mytrk1->GetPt()+mytrk2->GetPt()));
    if (TMath::Abs(dphi_r) < DPHI) dzed_conv_hist->Fill(dzed,0.5*(mytrk1->GetPt()+mytrk2->GetPt()));    

    if ( mytrk1->GetArm() != mytrk2->GetArm() ) return 1;
    if (TMath::Abs(dphi_r) >= DPHI) return 2;
    if(TMath::Abs(dzed) >= ZEDCUT)  return 3;

    the_conv_hist->Fill(dtheta_r,0.5*(mytrk1->GetPt()+mytrk2->GetPt()));
    phiv_conv_hist->Fill(phiv,0.5*(mytrk1->GetPt()+mytrk2->GetPt()));

    if (dtheta_r >= DTHETA) return 4;
    if (radius_r >= R_HI2)  return 4;

    if ( TMath::Abs(phiv)>=PHIVCUT ) return 5;

    if (radius_r >= R_HI3)   return 6;
    if (dtheta_r >= DTHETA2) return 6;

    if (TMath::Abs(dzed) >= ZEDCUT2) return 7;
    if (TMath::Abs(dphi_r) >= DPHI2) return 7;

    return 8;
}