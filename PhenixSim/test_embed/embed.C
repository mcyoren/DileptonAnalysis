#include "Riostream.h"

void embed(
	const int nevent = 100,
	const char *rdDST = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/berd/sim/output/embedding/CNTmerge_MB-0000231429-0002.root", 
	const char *mcDST = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/berd/sim/output/pisaToDST/Run7AuAu200_kaon/DST/dst_kaon_0_0.root", 
	const char *dstout      = "kaon.root", 
	const char *ntname      = "embed.root",
	const int   runnum      = 231429)
{
  
  cout << "Input RD merged dst      : " << rdDST << endl;
  cout << "Input MC dst             : " << mcDST << endl;
  cout << "Output evaluation ntuple : " << ntname << endl;
  cout << "Output dst name          : " << dstout << endl;
  cout << "Number of events to analyze: " << nevent << endl;
  
  cout << "You are running with: " << endl;
  gSystem->Exec("echo $OFFLINE_MAIN"); 
  
  float vtxmatch = 1.0; // match the vertex between real and sim
  float vtxmax = 30.0;

  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libcgl.so");
  gSystem->Load("libembed.so");
  gSystem->Load("libembedreco.so");
  gSystem->Load("libCrkPID.so");

  //gROOT->ProcessLine(".L embed_IOManager.C");

  gSystem->ListLibraries();

  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUN7AUAU200GEV",1);    // Run 7 200 GeV Au+Au
  rc->set_IntFlag("RUNNUMBER", runnum);  // for 200GEV RUN7
  rc->set_IntFlag("EVALUATIONFLAG", 1);  // Requested by EMCal
  rc->set_IntFlag("DCHREQFLAG", 0);
  rc->set_IntFlag("VERBOSITY", 0); 
  
  Fun4AllServer *se = Fun4AllServer::instance(); 
  se->Verbosity(0);
  
  SubsysReco *mixrec       = new MixEmbedreco("MIX");
  SubsysReco *bbcrec       = new BbcEmbedreco("BBC");
  SubsysReco *vtxrec       = new VtxReco("VTX");
  SubsysReco *padrec       = new PadEmbedreco("PAD");
  SubsysReco *dchrec       = new DchEmbedreco("DCH");
  SubsysReco *tofrec       = new TofEmbedreco("TOF");
  SubsysReco *tofwrec      = new TofEmbedreco("TOFW");
  SubsysReco *crkrec       = new CrkEmbedreco("CRK");
  SubsysReco *emcrec       = new EmcEmbedreco("EMC");
  SubsysReco *accrec       = new AccEmbedreco("ACC");

  //HBD and ACC geometries are not in.  Expect to see hundreds of
  //"length of one vector = 0" errors coming from PHDetectorGeometry
  //trying to instantiate PHPanel objects from null vectors.
  SubsysReco *cglrec       = new CglEmbedreco("CGL");
  //SubsysReco *cglrec       = new CglReco("CGL");
  SubsysReco *ringrec      = new RingEmbedreco("RING");  
  EvaEmbedreco *evarec     = new EvaEmbedreco("ChargedEVA");

  mixrec->Verbosity(1);
  dchrec->Verbosity(0);
  evarec->Verbosity(10);

  // EmbedVertexSelect is now registered with the fun4all server, not
  // the input manager, because the server reads all inputs for an
  // event before deciding whether to ditch it. The manager only
  // ditches the event from one input stream (the one you register
  // with), which is undesirable if the DSTs were pre-engineered to be
  // synchronized. This selector needs a node to work on, and in the
  // current implementation, it shouldn't matter whether whether it's
  // REAL or SINGLE since the event from both inputs will get rejected
  // upon z-vertex mismatch.
  EmbedVertexSelect *vtxsel = new EmbedVertexSelect("VTXSEL","REAL");
  vtxsel->SetVertexRange(vtxmatch); //match vertex with in vtxmatch cm
  vtxsel->Verbosity(1); // verbosity, was set to 1

  // MC TopNode name
  rc->set_CharFlag("EMBED_MC_TOPNODE","SINGLE");
  // real event TopNode name
  rc->set_CharFlag("EMBED_REAL_TOPNODE","REAL");
  
  // If one arm has no MC hits, then kick out the hits from real DST,
  // the reconstruction will be much faster.
  rc->set_IntFlag("EMBED_KickOutHitsToSpeedupReco",0);
  //T0 information for DC East and DC West
  rc->set_FloatFlag("EMBED_DCEASTT0",40);
  rc->set_FloatFlag("EMBED_DCWESTT0",39);
  //rc->set_IntFlag("VERBOSITY", 0);
  
  // The output evaluation ntuples for charged tracks.
  // you can add your own evaluation modules to EvaEmbedreco class
  rc->set_CharFlag("EMBED_CHARGED_EVAOUT",ntname);

  // Or you can set the output file name directly: 
  //PHEmbedHistogrammer::instance()->setFileName(ntname);

  //se->registerSubsystem(vtxsel);
  se->registerSubsystem(mixrec);  
  se->registerSubsystem(bbcrec);
  se->registerSubsystem(vtxrec);
  se->registerSubsystem(padrec);
  se->registerSubsystem(dchrec);
  se->registerSubsystem(tofrec);
  se->registerSubsystem(tofwrec);
  se->registerSubsystem(crkrec);
  se->registerSubsystem(emcrec);
  se->registerSubsystem(accrec);
  // cglrec currently gives errors...see comment above.
  se->registerSubsystem(cglrec); 
  se->registerSubsystem(ringrec);
  
  // Evaluation module should always be the last one.
  // There could be multiple evaluation modules.
  se->registerSubsystem(evarec);
  
  // If an instance of CentralTrackReco was made:
  // se->registerSubsystem(central);
  Fun4AllInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1","DST","SINGLE"); 
  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTin2","DST","REAL");
  
  in2->registerSubsystem(vtxsel);

  // Add the DSTs directly instead of in a list.
  in1->AddFile(mcDST);
  in2->AddFile(rdDST);

  se->registerInputManager(in1);
  se->registerInputManager(in2);

  //if( dstout ) DST_IOManager(dstout, se);
  
  cout << "running ..." << endl;
  se->run(nevent);
  se->End();
  
  cout<<"finished"<<endl;
  gSystem->Exit(0);
  
}
