#include "Calib.h"
void Calib(const TString inname = inFile[0],  int itread = 0, int ntreads = 1)
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
  const int check_veto = 1;
  const int fill_inv_mass = 0;
  const int istruehitsigmacounter = 0;


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
    
    for (int i = 0; i < myevent->GetNtrack()*fill_QA_hadron_hists; i++)
    {
      MyDileptonAnalysis::MyElectron trk = *myevent->GetEntry(i);
      if (trk.GetPtPrime()<-1.5||trk.GetPtPrime()>10.75) continue;
      MyDileptonAnalysis::MyHadron *newTrack = new MyDileptonAnalysis::MyHadron;
      newTrack->SetTrkId(i);
      newTrack->SetArm(trk.GetArm());
      newTrack->SetSect(trk.GetSect());
      newTrack->SetTrkQuality(63);
      newTrack->SetPt(trk.GetPtPrime());
      newTrack->SetPtPrime(trk.GetPtPrime());
      newTrack->SetQ(trk.GetChargePrime());
      newTrack->SetQPrime(trk.GetChargePrime());
      newTrack->SetPhiDC(trk.GetPhiDC());
      newTrack->SetPhi0(trk.GetPhi0Prime());
      newTrack->SetThe0(trk.GetThe0Prime());
      newTrack->SetPhi0Prime(trk.GetPhi0Prime());
      newTrack->SetThe0Prime(trk.GetThe0Prime());
      newTrack->SetZDC(trk.GetZDC());
      newTrack->SetAlpha(trk.GetAlphaPrime());
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
    if(fill_QA_hadron_hists) event_container->correct_beam_offset();
    if(fill_QA_hadron_hists) event_container->Associate_Hits_to_Hadrons(5);
    if(do_track_QA) event_container->FillQAHist(in_id);
    if(associate_hits)event_container->Associate_Hits_to_Leptons(5,5,5);
    if(istruehitsigmacounter) 
    {
      event_container->Associate_Hits_to_Leptons(2,5);
      event_container->Associate_Hits_to_Leptons(3,5);
      event_container->Associate_Hits_to_Leptons(4,5);
      event_container->Associate_Hits_to_Leptons(5,5);
      event_container->Associate_Hits_to_Leptons(6,5);
      event_container->Associate_Hits_to_Leptons(8,5);
      event_container->Associate_Hits_to_Leptons(10,5);
    }
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
    if(fill_true_DCA) event_container->FillTrueDCA();
    if(fill_TTree) event_container->FillTree();
    myevent->ReshuffleElectrons();
    if(fill_inv_mass) event_container->fill_inv_mass();
    
  
  }

  event_container->WriteOutFile();

}
