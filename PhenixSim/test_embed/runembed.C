#include "Riostream.h"
//This macro merges simulated DST defined by "mcdst" with real dst. 
//The output evaluation ntuple can be specified by variable "ntname"
void runembed(
	      const int nevent = 0, 
	      const char* rdDST =
	      //"/phenix/hl/tujuba/embedding/prod/dsts/DST_MinBias_run7AuAu_Central_200GeV_pro79-0000235613-0000.root",
	      "DST_EMBED_run11auau_200GeV_pro101-0000347381-0005_good.root",
	      const char* mcDST = 
	      "simDST_piplus_0_Run11_SL6_530.root", 
	      const char* dstout = "dstout.root",
	      const char* ntname = "ntout_test.root",
	      int runYear = 11
	      ) {
  if(!runYear) break;
  /*
  char realdstlist[500];
  char mcdstlist[500];
  char dstout[500];
  char ntname[200];
  
  sprintf( mcdstlist,"simDSTkplus.list");
  sprintf( realdstlist,"realDST.list");
  sprintf( dstout, "test_embedDst.root");
  sprintf( ntname,  "dummy.root" );

  cout << "Input RD merged dst list : " << realdstlist << endl;
  cout << "Input MC merged dst list : " << mcdstlist << endl;
  cout << "Output evaluation ntuple : " << ntname << endl;
  cout << "Output dst name          : " << dstout << endl;
  */

  cout << "Input RD merged dst      : " << rdDST << endl;
  cout << "Input MC dst             : " << mcDST << endl;
  cout << "Output evaluation ntuple : " << ntname << endl;
  cout << "Output dst name          : " << dstout << endl;
  
  cout << "You are running with: " << endl;
  gSystem->Exec("echo $OFFLINE_MAIN"); 
  
  float vtxmatch = 1.0;
  float vtxmax = 30.0;
  
  Int_t magField = 5; //no field integral.dat past run 10

  string base = "/afs/rhic/phenix/software/calibration/";
  string targ = " fieldIntegral.dat";
  string cmd1 = "ln -fs " + base + "run2001/fieldIntegral.dat" + targ;
  string cmd2 = "ln -fs " + base + "run2003/fieldIntegral.dat" + targ;
  string cmd3 = "ln -fs " + base + "run2004/fieldIntegral++.dat" + targ;
  string cmd4 = "ln -fs " + base + "run2004/fieldIntegral+-.dat" + targ;
  string cmd5 = "ln -fs " + base + "run2007/fieldIntegral+-.dat.run07" + targ;

  string cmd6 =  "ln -fs " + base + "run2010/fieldIntegral+-.dat.run10" + targ;
  string cmd7 =  "ln -fs " + base + "run2010/fieldIntegral+-.dat.run10" + targ;
  string cmd8 =  "ln -fs " + base + "run2010/fieldIntegral+-.dat.run10" + targ;


  /*
  switch(magField) {
    case 1:  // Run1/Run2 (3D01)
      gSystem->Exec(cmd1.c_str());
      cout << "\nMagnetic field map fieldIntegral.dat file set for Run1/Run2 (3D01)"
	   << endl;
      break;
    case 2:  // Run3 (3D03)
      gSystem->Exec(cmd2.c_str());
      cout << "\nMagnetic field map fieldIntegral.dat file set for Run3 (3D03), "
	   << "same as for Run4 (3D+0)"
	   << endl;
      break;
    case 3:  // Run4 (3D++)
      gSystem->Exec(cmd3.c_str());
      cout << "\nMagnetic field map fieldIntegral.dat file set for Run4 (3D++), "
	   << "same polarity in both coils" 
	   << endl;
      break;
    case 4: // Run4 (3D+-) not yet implemented
      gSystem->Exec(cmd4.c_str());
      cout << "\nMagnetic field map fieldIntegral.dat file set for 3D+-, "
	   << "fully reversed polarity in two coils" 
	   << endl;
      break;
    case 5: // Run7 (3D+-)
      gSystem->Exec(cmd5.c_str());
      cout << "\nMagnetic field map fieldIntegral.dat file set for 3D+-, "
	   << "fully reversed polarity in two coils as in run7" 
	   << endl;
      break;
  case 6:
    gSystem->Exec(cmd6.c_str());
    cout << "\nMagnetic field map fieldIntegral.dat file set for 3D+-, "
	 << "fully reversed polarity in two coils as in run10, maybe?" 
	 << endl;
    break;
  case 7:
    gSystem->Exec(cmd7.c_str());
    cout << "\nMagnetic field map fieldIntegral.dat file set for 3D+-, "
	 << "fully reversed polarity in two coils as in run11, maybe?" 
	 << endl;
    break;
  case 8:
    gSystem->Exec(cmd8.c_str());
    cout << "\nMagnetic field map fieldIntegral.dat file set for 3D+-, "
	 << "fully reversed polarity in two coils as in run14, maybe?" 
	 << endl;
    break;
    default:
      cout << "\nmagField value " << magField 
	   << " is not recognized; job is aborting" << endl;
      return;
    }
  */
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");

  gSystem->Load("libcgl.so");
  gSystem->Load("libembed.so");
  gSystem->Load("libembedreco.so");

  gSystem->Load("libCrkPID.so");

  gSystem -> Load("libpad.so");
  gSystem -> Load("libdch.so");
  gSystem -> Load("libtof.so");
  gSystem -> Load("libtofw.so");
  gSystem -> Load("libemc.so");
  gSystem -> Load("libemcEmbed4all.so");
  gSystem -> Load("libcompactCNT.so");
  gSystem -> Load("libcteval.so");
  gSystem -> Load("libMatchrecalRecoTest.so");
 
  gROOT->ProcessLine(".L embed_IOManager.C");

  gSystem->ListLibraries();

  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUN7AUAU200GEV",1);    // Run 7 200 GeV Au+Au
  //rc->set_IntFlag("RUN11AUAU200GEV",1);    // Run 11 200 GeV Au+Au
  //  rc->set_IntFlag("RUNNUMBER",+231428);  // for 200GEV RUN7
  rc->set_IntFlag("EVALUATIONFLAG", 1);  // Requested by EMCal
  rc->set_IntFlag("DCHREQFLAG", 0);
  rc->set_IntFlag("VERBOSITY", 0); 
  //   rc->set_IntFlag("AFSABSENT",1);
  
  Fun4AllServer *se = Fun4AllServer::instance(); 
  se->Verbosity(0);
  if(runYear == 10 || runYear == 11)
    {
      //SubsysReco *matching = new MatchrecalRecoTest();
      //se -> registerSubsystem(matching);
    }
  SubsysReco *mixrec       = new MixEmbedreco("MIX");
  SubsysReco *bbcrec       = new BbcEmbedreco("BBC");
  SubsysReco *vtxrec       = new VtxReco("VTX");
  SubsysReco *padrec       = new PadEmbedreco("PAD");
  SubsysReco *dchrec       = new DchEmbedreco("DCH");
  SubsysReco *tecrec       = new TecEmbedreco("TEC");
  SubsysReco *tofrec       = new TofEmbedreco("TOF");
  SubsysReco *crkrec       = new CrkEmbedreco("CRK");
  SubsysReco *emcrec       = new EmcEmbedreco("EMC");
  SubsysReco *accrec       = new AccEmbedreco("ACC");
 
  
  // HBD and ACC geometries are not in.  Expect to see hundreds of
  // "length of one vector = 0" errors coming from PHDetectorGeometry
  // trying to instantiate PHPanel objects from null vectors.
  SubsysReco *cglrec       = new CglEmbedreco("CGL");    
  SubsysReco *ringrec      = new RingEmbedreco("RING");  
  EvaEmbedreco *evarec     = new EvaEmbedreco("ChargedEVA");

  mixrec->Verbosity(1);
  dchrec->Verbosity(0);
  evarec->Verbosity(0);


  // This is the class which makes the Central Tracks nanoDST output
  // 22 corresponds to the version used in pro.78 for Run7 Au+Au
  //  SubsysReco *central      = new CentraltrackReco( 22 );

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
  vtxsel->Verbosity(0);

  // MC TopNode name
  rc->set_CharFlag("EMBED_MC_TOPNODE","SINGLE");
  // real event TopNode name
  rc->set_CharFlag("EMBED_REAL_TOPNODE","REAL");
  
  // If one arm has no MC hits, then kick out the hits from real DST,
  // the reconstruction will be much faster.
  rc->set_IntFlag ("EMBED_KickOutHitsToSpeedupReco",0);
  // T0 information for DC East and DC West
  rc->set_FloatFlag("EMBED_DCEASTT0",40);
  rc->set_FloatFlag("EMBED_DCWESTT0",39);
  rc->set_IntFlag("VERBOSITY", 0);
  
  // The output evaluation ntuples for charged tracks.
  // you can add your own evaluation modules to EvaEmbedreco class
  rc->set_CharFlag("EMBED_CHARGED_EVAOUT",ntname);

  // Or you can set the output file name directly: 
  // PHEmbedHistogrammer::instance()->setFileName(ntname);

  //  se->registerSubsystem(vtxsel);
  se->registerSubsystem(mixrec);  
  se->registerSubsystem(bbcrec);
  se->registerSubsystem(vtxrec);
  se->registerSubsystem(padrec);
  se->registerSubsystem(dchrec);
  se->registerSubsystem(tecrec);
  se->registerSubsystem(tofrec);
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

  Fun4AllInputManager *in1 = 
    new Fun4AllNoSyncDstInputManager("DSTin1","DST","SINGLE"); 
  Fun4AllInputManager *in2 = 
    new Fun4AllDstInputManager("DSTin2","DST","REAL");//real data tree
  
  in2->registerSubsystem(vtxsel);

  /*
  in1->AddListFile(mcdstlist);   //read into "SINGLE" Node  
  in2->AddListFile(realdstlist); //read into "REAL" Node
  */

  // Add the DSTs directly instead of in a list.
  in1->AddFile(mcDST);
  in2->AddFile(rdDST);

  se->registerInputManager(in1);
  se->registerInputManager(in2);
  
  // OutputManager setup functionsm originally commented out by Andrew
  //if( dstout ) DST_IOManager(dstout, se);
  
  cout << "running ..." << endl;
  se->run(nevent);
  se->End();
  
  cout<<"finished"<<endl;
  gSystem->Exit(0);
  
}

