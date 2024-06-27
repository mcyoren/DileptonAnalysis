#include "NewHitAssociation.h"
void NewHitAssociation(const char* inFile0, const char* outFile, int par = 0)
{

  TFile *input = new TFile(inFile0, "READ");
  if (!(input))
  {
    cout << "no input file" << endl;
    return;
  }
  const int remove_hadron_hits = 1;
  const int fill_QA_hadron_hists = 1;
  const int fill_QA_lepton_hists = 1;
  const int fill_TTree = 0;
  const int fill_d_dphi_hists = 0;
  const int fill_DCA_hists = 0;
  const int use_d_dphi_DCA = 0;
  const int do_track_QA = 0;
  const int do_reveal_hadron = 0;
  const int Use_ident = 0;
  const int fill_true_DCA = 1;
  const int check_veto = 1;
  const int is_only_conv = 0;
  const int fill_inv_mass = 0;
  const int use_only_vertex_hadrons = 1;

  TTree *T = (TTree *)input->Get("tree");
  TBranch *br = T->GetBranch("MyEvent");
  MyDileptonAnalysis::MyEventContainer *event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile(field_file_path);
  event_container->CreateOutFileAndInitHists(outFile, fill_QA_lepton_hists, fill_QA_hadron_hists, fill_TTree, fill_d_dphi_hists,
                                             fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto, fill_inv_mass);
  MyDileptonAnalysis::MyEvent *event = event_container->GetEvent();
  // event = 0;
  br->SetAddress(&event);

  cout << "Trees read!" << endl;

  int nevt = T->GetEntries();


  TH3D *vtxhist = new TH3D("vtxhist", "vtxhist", 250, -1, 1, 250, -1, 1, 10, 0, 10);

  for (int ievent = 0; ievent < nevt; ievent++)
  {
    if (ievent % 1000 == 0)
      cout << "Event: " << ievent << " / " << nevt << endl;
    event->ClearEvent();
    event_container->ClearEvent();

    br->GetEntry(ievent);
    if (ievent > 8000000)
      break;

    int n_hadrons = event->GetNhadron();
    for (int itrk = 0; itrk < n_hadrons*remove_hadron_hits; itrk++)
    {
      MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
      if (mytrk->GetPtPrime() < 0.5 || TMath::Abs(mytrk->GetPC3SDPHI())>2 || TMath::Abs(mytrk->GetPC3SDZ())>2) 
      {
        event->RemoveHadronEntry(itrk);
        n_hadrons--;
        itrk--;
      }
    }

    int n_electrons = event->GetNtrack();
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
      mytrk->ResetPrimes(event->GetVtxZ(),event->GetPreciseZ(),event->GetRunNumber());
      if (use_only_vertex_hadrons || ((mytrk->GetIsConv() < 8 || mytrk->GetRConv()>3) && is_only_conv))
      {
        event->RemoveTrackEntry(itrk);
        n_electrons--;
        itrk--;
      }
    }

     for (int i = 0; i < event->GetNhadron()*use_only_vertex_hadrons; i++)
    {
      MyDileptonAnalysis::MyHadron trk = *event->GetHadronEntry(i);
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
      event->AddTrack(newTrack);
    }
    
    if (event->GetNtrack() < 1) continue;

    //event_container->SetEvent(event);
    
    
    event_container->Associate_Hits_to_Leptons();

    if (check_veto)
      event_container->CheckVeto();
    if (fill_true_DCA)
      event_container->FillTrueDCA();


    n_electrons = event->GetNtrack();
        for (int itrk = 0; itrk < n_electrons; itrk++)
        {
            MyDileptonAnalysis::MyElectron mytrk = *event->GetEntry(itrk);
            bool do_reshuf = false;

            if (mytrk.GetHitCounter(0) < 1 || mytrk.GetHitCounter(1) < 1 ||
                (mytrk.GetHitCounter(2) < 1 && mytrk.GetHitCounter(3) < 1))
                do_reshuf = true;
            if(TMath::Abs(mytrk.GetDCA2())>50||mytrk.GetGhost()>0) do_reshuf = true;

            if (do_reshuf)
            {
                event->RemoveHadronEntry(itrk);
                event->RemoveTrackEntry(itrk);
                //event->AddElecCand(&mytrk);
                n_electrons--;
                itrk--;
                continue;
            }
        }

    vtxhist->Fill(event->GetPreciseX(),event->GetPreciseY(),event->GetRunNumber());
    
    if(use_only_vertex_hadrons==0)event->ReshuffleElectrons();

    if (event->GetNtrack() < 2-use_only_vertex_hadrons)
      continue;

    for (int itrk = 0; itrk <  event->GetNhadron(); itrk++)
    {
      MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
      if (remove_hadron_hits && (mytrk->GetPtPrime() < 0.5 || TMath::Abs(mytrk->GetPC3SDPHI())>2 || TMath::Abs(mytrk->GetPC3SDZ())>2)) 
      {
        std::cout<<"WTF"<<std::endl;
      }
      //std::cout<<"before: "<<mytrk->GetPhi0()<<" "<<mytrk->GetPhi0Prime()<<" "<<mytrk->GetThe0()<<" "<<mytrk->GetThe0Prime()<<" "<<std::endl;
      mytrk->ResetPrimes(event->GetVtxZ(),event->GetPreciseZ(),event->GetRunNumber());
      //std::cout<<"after: "<<mytrk->GetPhi0()<<" "<<mytrk->GetPhi0Prime()<<" "<<mytrk->GetThe0()<<" "<<mytrk->GetThe0Prime()<<" "<<std::endl;
    }

    if (remove_hadron_hits)
    {
      event_container->Associate_Hits_to_Hadrons();
    }
    
    if (use_d_dphi_DCA)
      event_container->FillDphiHists();
    if (do_reveal_hadron)
      event_container->Reveal_Hadron();
    if (fill_TTree)
      event_container->FillTree();
    if(fill_inv_mass)
      event_container->fill_inv_mass();

  }
  event_container->WriteOutFile();
}
