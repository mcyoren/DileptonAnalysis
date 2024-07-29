#include "InvMass.h"
void InvMass(const TString inname = inFile[0],  int itread = 0, int ntreads = 1)
{
  std::cout<<"start"<<std::endl;
  TFile *input = new TFile(inname, "READ");
  if (!(input))
  {
    cout << "no input file" << endl;
    return;
  }
  const int associate_hits = 1;
  const int remove_hadron_hits = 0;
  const int fill_QA_hadron_hists = 0;
  const int fill_QA_lepton_hists = 0;
  const int fill_TTree = 0;
  const int fill_d_dphi_hists = 0;
  const int fill_DCA_hists = 0;
  const int use_d_dphi_DCA = 0;
  const int do_track_QA = 0;
  const int do_reveal_hadron = 0;
  const int Use_ident = 0;
  const int fill_true_DCA = 0;
  const int check_veto = 0;
  const int fill_inv_mass = 1;

  char outname[200];
  sprintf(outname,"kek_%d.root",itread);
  TTree *T = (TTree *)input->Get("tree");
  TBranch *br = T->GetBranch("MyEvent");
  MyDileptonAnalysis::MyEventContainer *event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile("../../ee_QA/AnaTrain/Run14AuAuLeptonComby/field_map.root");
  event_container->CreateOutFileAndInitHists(outname,fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,
                                               fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto, fill_inv_mass);
  MyDileptonAnalysis::MyEvent *myevent = 0;                                            
  //event = 0;
  br->SetAddress(&myevent);

  cout << "Trees read!" << endl;

  const int nevt = T->GetEntries();
  const int beggin = nevt * (itread - 1) / ntreads;
  const int endish = nevt * itread / ntreads;

  TH3D *myhist0 = new TH3D("myhist0","myhist0",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist1 = new TH3D("myhist1","myhist1",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist2 = new TH3D("myhist2","myhist2",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist3 = new TH3D("myhist3","myhist3",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist4 = new TH3D("myhist4","myhist4",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist5 = new TH3D("myhist5","myhist5",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist6 = new TH3D("myhist6","myhist6",500,-2000,2000,500,0,5,10,0,10);
  TH3D *myhist7 = new TH3D("myhist7","myhist7",500,-2000,2000,500,0,5,10,0,10);
 
  for (int ievent = beggin; ievent < endish; ievent++)
  {
    if ((ievent -beggin) % 50000 == 0)
      cout << "ithread, iEvent, N_events: " << itread<< ",  " << ievent -beggin<< " / " << nevt/ntreads << endl;
    myevent->ClearEvent();
    br->GetEntry(ievent);
    if (ievent - beggin > 3500000)
      break;


    if(myevent->GetNtrack()<1 && !do_track_QA) continue;


    myevent->SetVtxZ(myevent->GetPreciseZ());
    /*myevent->SetCentrality(event->GetCentrality());
    myevent->SetPreciseX(event->GetPreciseX());
    myevent->SetPreciseY(event->GetPreciseY());
    myevent->SetPreciseZ(event->GetPreciseZ());
    myevent->SetEvtNo(event->GetEvtNo());
    myevent->SetRunNumber(0);
    myevent->SetBBCcharge(event->GetBBCcharge());
    myevent->SetBBCtimeN(event->GetBBCtimeN());

    for (int i = 0; i < event->GetNtrack(); i++)
    {
      MyDileptonAnalysis::MyElectron trk = *event->GetEntry(i);
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
      newTrack->SetTOFDPHI(trk.GetTOFDPHI());
      newTrack->SetTOFDZ(trk.GetTOFDZ());
      newTrack->SetTOFE(trk.GetTOFE());
      newTrack->SetEmcTOF(trk.GetEmcTOF());
      newTrack->SetCrkphi(trk.GetCrkphi());
      newTrack->SetCrkz(trk.GetCrkz());
      newTrack->SetChi2(trk.GetChi2());
      newTrack->SetN0(trk.GetN0());
      newTrack->SetDISP(trk.GetDisp());
      newTrack->SetNPE0(trk.GetNpe0());
      newTrack->SetEmcdz_e(trk.GetEmcdz_e());
      newTrack->SetEmcdphi_e(trk.GetEmcdphi_e());
      newTrack->SetSect(trk.GetSect());
      newTrack->SetMcId(trk.GetMcId());
      myevent->AddTrack(newTrack);
    }

    for (int i = 0; i < event->GetNtrack()*fill_QA_hadron_hists; i++)
    {
      MyDileptonAnalysis::MyElectron trk = *event->GetEntry(i);
      if (trk.GetPtPrime()<-1.5||trk.GetPtPrime()>10.75) continue;
      MyDileptonAnalysis::MyHadron *newTrack = new MyDileptonAnalysis::MyHadron;
      newTrack->SetTrkId(i);
      newTrack->SetArm(trk.GetArm());
      newTrack->SetSect(trk.GetSect());
      newTrack->SetTrkQuality(63);
      newTrack->SetPt(trk.GetPt());
      newTrack->SetPtPrime(trk.GetPtPrime());
      newTrack->SetQ(trk.GetCharge());
      newTrack->SetQPrime(trk.GetCharge());
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
      newTrack->SetDISP(trk.GetDisp());
      newTrack->SetNPE0(trk.GetNpe0());
      newTrack->SetEmcdz_e(trk.GetEmcdz_e());
      newTrack->SetEmcdphi_e(trk.GetEmcdphi_e());
      newTrack->SetSect(trk.GetSect());
      newTrack->SetMcId(trk.GetMcId());
      myevent->AddHadron(newTrack);
    }
    if(myevent->GetNhadron()<1 && fill_QA_hadron_hists ) continue;


    for (int i = 0; i < event->GetNVTXhit()*associate_hits; i++)
    {
      MyDileptonAnalysis::MyVTXHit oldhit =*event->GetVTXHitEntry(i);
      MyDileptonAnalysis::MyVTXHit *newHit = new MyDileptonAnalysis::MyVTXHit;
      newHit->SetClustId(i);
      if(oldhit.GetSensor() == 0) continue;
      newHit->SetLayer(oldhit.GetLayer());
      newHit->SetLadder(oldhit.GetLadder());
      newHit->SetSensor(oldhit.GetSensor());
      newHit->SetXHit(oldhit.GetXHit());
      newHit->SetYHit(oldhit.GetYHit());
      newHit->SetZHit(oldhit.GetZHit());
      //newHit->SetPolars(myevent->GetPreciseX(), myevent->GetPreciseY(), myevent->GetPreciseZ());
      newHit->SetiLayerFromR();
      myevent->AddVTXHit(newHit);
    }
    */

    int n_hits = myevent->GetNVTXhit();
    for (int ihit = 0; ihit < n_hits; ihit++)
    {
      MyDileptonAnalysis::MyVTXHit myhit = *myevent->GetVTXHitEntry(ihit);

     if (myhit.GetSensor() == 1)
     {
         myevent->RemoveVTXHitEntry(ihit);
         n_hits--;
         ihit--;
         continue;
     }
    }
    
    event_container->SetEvent(myevent);

    event_container->correct_beam_offset();
    if(fill_QA_hadron_hists) event_container->Associate_Hits_to_Hadrons(400);
    if(do_track_QA) event_container->FillQAHist(in_id);
    if(associate_hits)event_container->Associate_Hits_to_Leptons(2,2,5);
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
    event_container->CheckVeto();
    if(fill_true_DCA) event_container->FillTrueDCA();
    if(fill_TTree) event_container->FillTree();
    myevent->ReshuffleElectrons();
    if(fill_inv_mass) event_container->fill_inv_mass();
    myevent->ClearEvent();
    //event->ClearEvent();
  }
  event_container->WriteOutFile();
}
