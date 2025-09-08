void Run_Run14AuAuLeptonComby(char *outFile = "Phi_om.root") {
    gSystem->Load("libRun14AuAuLeptonEvent.so");
    gSystem->Load("libRun14AuAuLeptonConvReco.so");
    gSystem->Load("libRun14AuAuLeptonComby.so");  
    gSystem->Load("libTOAD");  
    //gSystem->Load("libRPCalibRun.so");
    gSystem->Load("libUltraLight");  
    gSystem->Load("libCabanaBoy");
    
    TOAD toad_loader("PhotonConversionAnalysis");
    toad_loader.SetVerbosity(4);
    string lookupfile_location = toad_loader.location("lookup_3D_one_phi.root");
    std::cout<<"lookup_3D_one_phi.root is located at "<<lookupfile_location.c_str()<<std::endl;

    recoConsts *reco_consts =  recoConsts::instance();
    reco_consts->set_IntFlag("Remove_hadron_hits", 0);
    reco_consts->set_IntFlag("Fill_QA_hadron_hists", 0);
    reco_consts->set_IntFlag("Fill_QA_lepton_hists", 0);
    reco_consts->set_IntFlag("Fill_ddphi_hadron", 0);
    reco_consts->set_IntFlag("Fill_TTree", 0);
    reco_consts->set_IntFlag("Fill_d_dphi_hists", 0);
    reco_consts->set_IntFlag("Fill_DCA_hists", 0);
    reco_consts->set_IntFlag("Use_ident", 2);
    reco_consts->set_IntFlag("Do_track_QA", 0);
    reco_consts->set_IntFlag("Do_electron_QA", 1);
    reco_consts->set_IntFlag("Do_reveal_hadron", 0);
    reco_consts->set_IntFlag("Fill_flow", 0);
    reco_consts->set_IntFlag("Fill_true_DCA", 1);
    reco_consts->set_IntFlag("Check_Veto", 0);
    reco_consts->set_IntFlag("fill_inv_mass", 1);
    reco_consts->set_IntFlag("do_reco_vertex", 1);
    reco_consts->set_IntFlag("do_conv_dalitz_finder", 2);
       
    SubsysReco *reco = new Run14AuAuLeptonCombyReco(outFile, lookupfile_location.c_str());

    cbMasterCutter *mc = new Run14AuAuLeptonCombyCutter();
    cbMasterHistos *mh = new Run14AuAuLeptonCombyHistos();

    CabanaBoy *cb = new CabanaBoy(5,1,1, "Run14AuAuLeptonComby");
    //CabanaBoy *cb = new CabanaBoy(5,1,1, "Run14AuAuLeptonComby");
	
	cb->SetHistoFileName(outFile);
	cb->setZVertexMax(10);
	cb->setReactionPlaneSelectionType(CabanaBoy::ReactionPlaneNotUsed);//
	//cb->setReactionPlaneSelectionType(CabanaBoy::ReactionPlaneBBCSNPsi2);//ReactionPlaneNotUsed
	cb->setCentralitySelectionType(CabanaBoy::CentralityTypeRun12CuAu);
	cb->setCuts(mc);
	cb->setHistos(mh);
	cb->setPoolType(CabanaBoy::MultiAkibaPools);
	cb->setNSubPools(4, 20, 1);
	cb->setFastMom(false);
	cb->setPoolDepth(500);
	cb->setMixingType11(true);  
	cb->setMixingType12(true);
	cb->setMixingType22(true);
	cb->Verbosity(0);

    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);

    SubsysReco *mstr = se->getSubsysReco("MASTERRECALIBRATORMANAGER");
    se->unregisterSubsystem(mstr);
   
    MasterRecalibrator *mr = new MasterRecalibrator();
    mr->Unlock(0);
    //mr->RemoveRecalibrator("MomChangeRecalReco");
    //mr->RemoveRecalibrator("Run14AuAu200PC2MatchRecal");
    //mr->RemoveRecalibrator("Run14AuAu200EMCMatchRecal");
    mr->RemoveRecalibrator("SvxCentralTrackRecalReco");
    mr->RemoveRecalibrator("SvxCentralTrackReFit");
    mr->RemoveRecalibrator("SvxCentralTrackFitRecal");

    //mr->RemoveRecalibrator("EmcTofWalkRecalReco");//EmcGenericEScaleRecalReco
    //mr->RemoveRecalibrator("EmctofrecalReco");//EmcPidrecalReco
    //mr->RemoveRecalibrator("EmcGenericDeadRecalReco");

    se->registerSubsystem(mr);


    //RPReadCalibTree *readT = new RPReadCalibTree();
    //readT->setTreeFileRecent("RP_recent_run14pro106_newcent_merge.root");
    //readT->setTreeFileFlat("RP_flat_run14pro106_newcent_merge.root");
    //readT->setTOADname("hachiya/15.08");	
    //se->registerSubsystem(readT);

    se->registerSubsystem(reco);
    se->registerSubsystem(cb);
}


void InputData(vector<string> &indata) {
  //indata.push_back("EWG");
  //indata.push_back("PWG");
  indata.push_back("CNT");
  indata.push_back("DST_SVX");
  //indata.push_back("DST_EVE");
  return;
}

