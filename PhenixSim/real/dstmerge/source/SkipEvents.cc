#include <iostream>
#include "Fun4AllReturnCodes.h"
#include "phool.h"
#include "recoConsts.h"
#include <getClass.h>
#include "VtxOut.h"
#include "PHPoint.h"
#include "SkipEvents.h"
#include "TrigLvl1.h"
#include "PHGlobal.h"

//==============================================================
SkipEvents::SkipEvents(int nevents) 
{
  ThisName = "SkipEvents";
  std::cout << ThisName << " constructor called " << std::endl;
  EventNumber=0;
  Events2Skip = nevents;
  minZVertex=-99999.;
  maxZVertex= 99999.;
}
//==============================================================
int SkipEvents::InitRun(PHCompositeNode *topNode)
{
  myfile.open("vertexes.txt");
  return 0;
}
//==============================================================
int SkipEvents::process_event(PHCompositeNode *topNode) 
{
  EventNumber++;
  if(EventNumber<=Events2Skip) {return ABORTEVENT;} 

  VtxOut *d_vtx = findNode::getClass<VtxOut>(topNode,"VtxOut");
  if(!d_vtx) { std::cerr << PHWHERE << "VtxOut node not found." << std::endl; return DISCARDEVENT;}
  float xvtx0 =  (d_vtx->get_Vertex()).getX();
  float yvtx0 =  (d_vtx->get_Vertex()).getY();
  float zvtx0 =  (d_vtx->get_Vertex()).getZ();
  float xvtx  =  (d_vtx->get_Vertex("SIM")).getX();
  float yvtx  =  (d_vtx->get_Vertex("SIM")).getY();
  float zvtx  =  (d_vtx->get_Vertex("SIM")).getZ();
  float xvtx1 =  (d_vtx->get_Vertex("BBC")).getX();
  float yvtx1 =  (d_vtx->get_Vertex("BBC")).getY();
  float zvtx1 =  (d_vtx->get_Vertex("BBC")).getZ();
  std::string s_vtx = d_vtx->which_Vtx();
  if(verbosity>0) {
    std::cout << "SIM Event vertex = " << xvtx << " " << yvtx << " " << zvtx << std::endl;
    std::cout << "BBC Event vertex = " << xvtx1 << " " << yvtx1 << " " << zvtx1 << std::endl;
    std::cout << "Default Event vertex = " << xvtx0 << " " << yvtx0 << " " << zvtx0 << " " << s_vtx << std::endl;
  }
  
  const TrigLvl1 *Trig =
    findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  const int trigscaled_on = Trig->get_lvl1_trigscaled_bit(4);
  
  const PHGlobal *globalCNT =
    findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  
  const float centrality = globalCNT->getCentrality();
  const float zdc1 = globalCNT->getZdcEnergyN();
  const float zdc2 = globalCNT->getZdcEnergyS();
  const float zdcz = globalCNT->getZdcZVertex();
  const bool isZDCOK = (zdc1 > 0 && zdc2 > 0 && zdcz > -9000);
  if(verbosity>0) std::cout<<centrality<<" "<<isZDCOK<<std::endl;


    if (!Trig || !trigscaled_on || !globalCNT || !isZDCOK || centrality < 0 )
        return ABORTEVENT;
  if(zvtx0>minZVertex && zvtx0<maxZVertex) {
    myfile << xvtx0 << " " << yvtx0 << " " << zvtx0 << " " << centrality << "\n"; 
    // this is done to properly handle NAN
  }
  else {
    if(verbosity>0) {std::cout << "!!! Event ABORTED !!!" << std::endl;}
    return ABORTEVENT;
  }

  return EVENT_OK;
}
//==============================================================
int SkipEvents::End(PHCompositeNode *topNode)
{
  myfile.close();
  return 0;
}
