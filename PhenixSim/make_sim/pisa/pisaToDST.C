// $Id: pisaToDST.C,v 1.52 2007/11/28 14:04:09 hpereira Exp $
/*!
   \file pisaToDST.C
   \brief pisa to DST reconstruction chain
   \author <a href="mailto:pereira@hep.saclay.cea.fr">Hugo Pereira</a>
   \version $Revision: 1.52 $
   \date $Date: 2007/11/28 14:04:09 $
*/

void pisaToDST(

  Int_t nEvents = 100, 
  char *filein="PISAEvent.root",
  char *dstout = "simDST.root", 
  // this is an arbitrary 200 GeV run number. 
  // It was selected using basic QA criterion from the DAQ database.
  int run_number = 289000 
  
) {
 
  // print output
  cout << "pisaToDST - nEvents: " << nEvents << endl;
  cout << "pisaToDST - filein: " << filein << endl;
  if( dstout ) cout << "pisaToDST - dstout: " << dstout << endl;
  cout << "pisaToDST - run_number: " << run_number << endl;
  cout << endl;
   
  // allow to disable muon simulations
  // they are enabled by default
  bool do_muon_arms = true;
  
  // load libraries
  gSystem->Load("libfun4all");
  gSystem->Load("libmutoo_subsysreco");
  gSystem->Load("libfun4allfuncs"); 
  gSystem->Load("libsimreco");    

  gROOT->ProcessLine(".L pisaToDST_IOManager.C");
    
  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
  recoConsts *rc = recoConsts::instance();

  // 2 means PISA-To-DST
  rc->set_IntFlag("SIMULATIONFLAG",2); 
  
  // disable embedding
  rc->set_IntFlag("EMBEDFLAG",0); 
  
  // Reference run number used in 2007 Au+Au 200 GeV 
  rc->set_IntFlag("RUNNUMBER",run_number); 
  
  // Requested by EMCal
  rc->set_IntFlag("EVALUATIONFLAG", 1); 
  
  // Run flag
  // rc->set_IntFlag("RUN9PP500GEV",1);        
  rc->set_IntFlag("RUN9PP200GEV",1);        

  // this should be moved to the Init method of TofSimreco
  rc->set_FloatFlag("TOFTIMINGRESOLUTION", 0.100);  
  
  /*  
  Flags to abort event if required number of GEANT hits is not present in the subsystem
  Defaults are all 0 except for the Drift Chamber
  default setting is 3 Drift Chamber wire plane hits  
  */
  rc->set_IntFlag("DCHREQFLAG", 3);   
  rc->set_IntFlag("PC1REQFLAG", 0);
  rc->set_IntFlag("PC2REQFLAG", 0);
  rc->set_IntFlag("PC3REQFLAG", 0);
  rc->set_IntFlag("TOFREQFLAG", 0);
  
  // not yet operational
  rc->set_IntFlag("EMCREQFLAG", 0);
  
  // assume AFS is present as at RCF
  rc->set_IntFlag("AFSABSENT", 0);                                
    
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
    
  ///////////////////////////////////////////
  // Activate the subsystems
  //////////////////////////////////////////
  
  // run header and trigger setting
  se->registerSubsystem( new HeadSimreco() );
  se->registerSubsystem( new TrigSimreco() );
 
  // event counter
  if( true ) 
  {
    MuonCounter* counter = new MuonCounter();
    se->registerSubsystem( counter );
  }
  
  // BBC simReco
  se->registerSubsystem(new BbcSimreco("BBC"));
  
  // pisa is used as an input vertex.
  // it overwrites the contents of the BBC out node.
  VtxSimreco* vtx_sim = new VtxSimreco();
  vtx_sim->UseVtx( VtxSimreco::PISA );
  vtx_sim->SmearZ( true );        // this is the default
  vtx_sim->UseXY( false );        // this is the default
  vtx_sim->OverwriteBBC( true );  // this is the default
  vtx_sim->ZVertexSigma( 2 );     // default error on the simulated vertex
  se->registerSubsystem( vtx_sim );

  // t0
  T0Simreco* t0_sim = new T0Simreco();
  t0_sim->T0Sigma(0.04);
  se->registerSubsystem( t0_sim );
  
  // pad chambers
  se->registerSubsystem(new PadSimreco("PAD"));
                                      
  // chiu Pad Vertexing Code, maybe needed for multiple vertices
  se->registerSubsystem(new PadVtxReco("PADVTX"));
                                      
  // The VtxReco works unchanged for both real and simulation events
  se->registerSubsystem(new VtxReco("VTX"));
                                      
  // The T0Reco works unchanged for both real and simulation events
  se->registerSubsystem(new T0Reco());
                                      
  // As of January 2, 2004 the Dch has uninitialized variable warnings from Valgrind
  // There are also log file output warning messages
  se->registerSubsystem( new DchSimreco("DCH") );
  
  // Time expansion chamber
  se->registerSubsystem( new TecSimreco("TEC"));
  
  // Time of flight detector
  se->registerSubsystem(new TofSimreco("TOF"));
  
  // Tof west
  se->registerSubsystem(new TfwSimreco("TFW"));
  
  // HBD
  se->registerSubsystem(new HbdSimreco("HBD"));
  
  // RICH
  se->registerSubsystem(new CrkSimreco("CRK"));
  
  // Aerogel subsystem as per e-mail from Narumi Kurihara on May 13, 2005
  se->registerSubsystem(new AccSimreco("ACC"));
  se->registerSubsystem(new AccReco());
  
  // EMCal uses the real data class
  rc->set_FloatFlag("EMCTOWERLOWGAIN", 0.0015625);
  rc->set_FloatFlag("EMCTOWERHIGHGAIN", 0.0125);
  se->registerSubsystem( new EmcReco3() );
  
  // The CglReco works unchanged for both real and simulation events
  se->registerSubsystem(new CglReco("CGL"));
  
  //Aerogel cluster  (Needs to be after cglRec)
  se->registerSubsystem(new AccclusterReco());
  
  //  This is the class which makes the RICH Ring data structure
  se->registerSubsystem( new RingReco() );
  
  // This is the class which makes the Central Tracks nanoDST output
  // 22 corresponds to the version used in pro.78 for Run7 Au+Au
  se->registerSubsystem(new CentraltrackReco( 22 ));
  
  //  This is the class which makes the GlobalEvent data on the nanoDST output
  se->registerSubsystem(new GlobalReco());
  
  // This is the class which checks for charged particles going into EMCal
  se->registerSubsystem(new ChargedvetoReco());
  
  //added the DC based global evaluation module
  se->registerSubsystem( new McEvalSimreco() );

  // muon arm reconstruction
  if( do_muon_arms )
  {
    
    // unfortunately the muon arm need the sign of the magnetic field
    // and does not have yet the logic to retrieve it from the pisa event header
    // it is hard-coded here
    mMfmMT::setMapFileScale( 1.0 );
    MuonUtil::set_check_mapfile_scale( false );
    
    // pisa unpacking
    // The muon reconstruction is not performed
    // and only the simulated objects are created. 
    // the reconstruction itself is performed in the Fun4All_RecoDST_sim afterburner macro
    se->registerSubsystem( new MuonUnpackPisa() );
    
  }
  
  // dumper
  if( false ) 
  {
    gSystem->Load( "libphnodedump" );
    se->registerSubsystem( new Dumper() );
  }
  
  ///////////////////////////////////////////
  // InputManager 
  ///////////////////////////////////////////
  Fun4AllInputManager *input_manager = new Fun4AllPisaInputManager("PisaIn","TOP");
  se->registerInputManager(input_manager);
  
  ///////////////////////////////////////////
  // OutputManagers Set up functions  
  ///////////////////////////////////////////
  if( dstout ) DST_IOManager(dstout, se);     
  
  ///////////////////////////////////////////
  // open input file
  se->fileopen(input_manager->Name(),filein);
  
  // process input events
  gBenchmark->Start("eventLoop");     
  se->run(nEvents);             
  se->End();
  gBenchmark->Show("eventLoop");      
            
  // If you do not see this message, the job failed
  cout << "Completed reconstruction." << endl;
}
