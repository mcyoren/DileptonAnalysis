////////////////////////////////////////////////////////////////
//
//  includes all necessary code for the HELIOS packackage and 
//  defines string with full pathlength to local directory of HELIOS package
//  used for some data files
//  
// this string need to be updated for installation on new computer 
//  e.g.  C:/root_v6.22.06/macros/HELIOS/
//
//  Axel Drees 11/15/2021
//
#include "TROOT.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#ifndef HELIOSLibrary_h
#define HELIOSLibrary_h

#include <string>
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/PDG.h"
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/InteractionWithMaterial.h"
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/KinematicDistributions.h"                       
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/Particle.C"
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/PHENIXSetup.h"
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/PHENIXDetector.C"
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/WriteEvent.h"

#endif
