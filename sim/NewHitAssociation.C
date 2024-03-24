#include "NewHitAssociation.h"
void NewHitAssociation(int par = 0)
{

  TFile *input = new TFile(inFile[par], "READ");
  if (!(input))
  {
    cout << "no input file" << endl;
    return;
  }
  const int remove_hadron_hits = 0;
  const int fill_QA_hadron_hists = 0;
  const int fill_QA_lepton_hists = 1;
  const int fill_TTree = 0;
  const int fill_d_dphi_hists = 1;
  const int fill_DCA_hists = 0;
  const int use_d_dphi_DCA = 1;
  const int do_track_QA = 0;
  const int do_reveal_hadron = 0;
  const int fill_true_DCA = 1;
  const int check_veto = 1;


  TTree *T = (TTree *)input->Get("T");
  TBranch *br = T->GetBranch("MyEvent");
  MyDileptonAnalysis::MyEventContainer *event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile("../../ee_QA/AnaTrain/Run14AuAuLeptonComby/field_map.root");
  event_container->CreateOutFileAndInitHists("kek.root",fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,
                                               fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto);
  DileptonAnalysis::MyEvent *event = 0;                                             
  //event = 0;
  br->SetAddress(&event);

  cout << "Trees read!" << endl;

  int nevt = T->GetEntries();

  
 
  for (int ievent = 0; ievent < nevt; ievent++)
  {
    if (ievent % 1000 == 0)
      cout << "Event: " << ievent << " / " << nevt << endl;
    br->GetEntry(ievent);
    if (ievent > 100000)
      break;


    if(event->GetNtrack()<1) continue;
    MyDileptonAnalysis::MyEvent *myevent = new MyDileptonAnalysis::MyEvent;

    myevent->SetCentrality(1);
    myevent->SetPreciseX(event->GetPreciseX());
    myevent->SetPreciseY(event->GetPreciseY());
    myevent->SetPreciseZ(event->GetPreciseZ());
    myevent->SetEvtNo(event->GetEvtNo());
    myevent->SetRunNumber(0);
    myevent->SetBBCcharge(event->GetBBCcharge());
    myevent->SetVtxZ(event->GetVtxZ());
    myevent->SetBBCtimeN(event->GetBBCtimeN());

    for (int i = 0; i < event->GetNtrack(); i++)
    {
      DileptonAnalysis::MyTrack trk = event->GetEntry(i);
      MyDileptonAnalysis::MyElectron *newTrack = new MyDileptonAnalysis::MyElectron;
      newTrack->SetTrkId(i);
      newTrack->SetArm(trk.GetArm());
      newTrack->SetSect(trk.GetSect());
      newTrack->SetTrkQuality(63);
      newTrack->SetPt(trk.GetPt());
      newTrack->SetPtPrime(trk.GetPtPrime());
      newTrack->SetQ(trk.GetCharge());
      newTrack->SetQPrime(trk.GetChargePrime());
      newTrack->SetPhiDC(trk.GetPhiDC());
      newTrack->SetPhi0(trk.GetPhi0());
      newTrack->SetThe0(trk.GetThe0());
      newTrack->SetPhi0Prime(trk.GetPhi0Prime());
      newTrack->SetThe0Prime(trk.GetThe0Prime());
      newTrack->SetZDC(trk.GetZDC());
      newTrack->SetAlpha(trk.GetAlpha());
      newTrack->SetAlphaPrime(trk.GetAlphaPrime());
      newTrack->SetEmcId(trk.GetEmcId());
      newTrack->SetEcore(trk.GetEcore());
      newTrack->SetDep(trk.GetDep());
      newTrack->SetProb(trk.GetProb());
      newTrack->SetEmcdz(trk.GetEmcdz());
      newTrack->SetEmcdphi(trk.GetEmcdphi());
      newTrack->SetTOFE(trk.GetTOFE());
      newTrack->SetEmcTOF(trk.GetEmcTOF());
      newTrack->SetCrkphi(trk.GetCrkphi());
      newTrack->SetCrkz(trk.GetCrkz());
      newTrack->SetChi2(trk.GetChi2());
      newTrack->SetN0(trk.GetN0());
      newTrack->SetEmcdz_e(trk.GetEmcdz());
      newTrack->SetEmcdphi_e(trk.GetEmcdphi());
      myevent->AddTrack(newTrack);
    }
    for (int i = 0; i < event->GetNVTXhit(); i++)
    {
      DileptonAnalysis::MyVTXHit oldhit = event->GetVTXHitEntry(i);
      MyDileptonAnalysis::MyVTXHit *newHit = new MyDileptonAnalysis::MyVTXHit;
      newHit->SetClustId(i);
      newHit->SetLayer(oldhit.GetLayer());
      newHit->SetLadder(oldhit.GetLadder());
      newHit->SetSensor(oldhit.GetSensor());
      newHit->SetXHit(oldhit.GetXHit());
      newHit->SetYHit(oldhit.GetYHit());
      newHit->SetZHit(oldhit.GetZHit());
      newHit->SetPolars(myevent->GetPreciseX(), myevent->GetPreciseY(), myevent->GetPreciseZ());
      newHit->SetiLayerFromR();
      myevent->AddVTXHit(newHit);
    }
    

    
    event_container->SetEvent(myevent);

    event_container->Associate_Hits_to_Leptons(remove_hadron_hits);
    int n_electrons = myevent->GetNtrack()*remove_hadron_hits;
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
        MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
        if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
           (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 ))
        {
            myevent->RemoveTrackEntry(itrk);
            n_electrons--;
            itrk--;
            continue;
        }
        mytrk->ZeroHitCounters();
    }

    if(myevent->GetNtrack()<1  ) continue;
    
    if(remove_hadron_hits) 
    {
        event_container->Associate_Hits_to_Hadrons();
        event_container->Associate_Hits_to_Leptons();
    }
    if(use_d_dphi_DCA)  event_container->FillDphiHists();
    if(do_reveal_hadron) event_container->Reveal_Hadron();
    if(check_veto) event_container->CheckVeto();
    if(fill_TTree) event_container->FillTree();
    for (int i = 0; i < myevent->GetNtrack(); i++)
    {
      MyDileptonAnalysis::MyElectron *newTrack1 = myevent->GetEntry(i);
      if(newTrack1->GetChargePrime()<0) continue;
      for (int i = 0; i < myevent->GetNtrack(); i++)
    {
        MyDileptonAnalysis::MyElectron *newTrack2 = myevent->GetEntry(i);
        if(newTrack2->GetChargePrime()>0) continue;
        //std::cout<<newTrack1->GetDCA2()<<" "<<newTrack2->GetDCA2()<<std::endl;
        //std::cout<<"grad: "<<newTrack1->GetPhiDC()*180/3.14159<<" "<<newTrack2->GetPhiDC()*180/3.14159<<std::endl;
        //const float ptpair = newTrack1->GetPy()+newTrack2->GetPy();
        //if(newTrack1->GetPhiDC()<0.8) std::cout<<newTrack2->GetCrkphi()*180/3.14159-newTrack1->GetCrkphi()*180/3.14159<<" "<<ptpair<<std::endl;
        //else std::cout<<360-newTrack2->GetCrkphi()*180/3.14159-newTrack1->GetCrkphi()*180/3.14159<<" kek "<<ptpair<<std::endl;
    }
    }
  
  }
  event_container->WriteOutFile();
}
