#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"

#include "MyEvent_OLD.h"


ClassImp(DileptonAnalysis::MyEvent)
ClassImp(DileptonAnalysis::MyTrack)
ClassImp(DileptonAnalysis::MyVTXHit)
ClassImp(DileptonAnalysis::MyCluster)
ClassImp(DileptonAnalysis::McTrack)
ClassImp(DileptonAnalysis::McCluster)
ClassImp(DileptonAnalysis::MyPair)

using namespace std;

namespace DileptonAnalysis
{
	void DileptonAnalysis::MyEvent::ClearEvent()
	{
         evtno      = -99999;
         BBCcharge  = -99999;
         zvertex    = -99999;
         centrality = -99999;
         psi2_BBC    = -99999;
         psi2_FVTXA0    = -99999;
   
         preciseX = -99999;
         preciseY = -99999;
         preciseZ = -99999;
                  
         run_number = -99999;

         TrackList.clear();
         VTXHitList.clear();
         ClusterList.clear();
         McTrackList.clear();
         McClusterList.clear();
         PairList.clear();
	}

 
}
