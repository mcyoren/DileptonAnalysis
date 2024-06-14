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
  const int fill_QA_lepton_hists = 0;
  const int fill_TTree = 0;
  const int fill_d_dphi_hists = 0;
  const int fill_DCA_hists = 0;
  const int use_d_dphi_DCA = 0;
  const int do_track_QA = 0;
  const int do_reveal_hadron = 0;
  const int Use_ident = 1;
  const int fill_true_DCA = 1;
  const int check_veto = 1;
  const int is_only_conv = 0;
  const int fill_inv_mass = 1;

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

    int n_electrons = event->GetNtrack();
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
      mytrk->ResetPrimes(event->GetVtxZ(),event->GetPreciseZ(),event->GetRunNumber());
      if ((mytrk->GetIsConv() < 8 || mytrk->GetRConv()>3) && is_only_conv) 
      {
        event->RemoveTrackEntry(itrk);
        n_electrons--;
        itrk--;
      }
    }
    if (n_electrons < 1)
      continue;

    //event_container->SetEvent(event);
    
    
    event_container->Associate_Hits_to_Leptons();

    event->ReshuffleElectrons();

    if (event->GetNtrack() < 2)
      continue;

    vtxhist->Fill(event->GetPreciseX(),event->GetPreciseY(),event->GetRunNumber());

    int n_hadrons = event->GetNhadron();
    for (int itrk = 0; itrk < n_hadrons*remove_hadron_hits; itrk++)
    {
      MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
      if (mytrk->GetPtPrime() < 1.5 || TMath::Abs(mytrk->GetPC3SDPHI())>2 || TMath::Abs(mytrk->GetPC3SDZ())>2) 
      {
        event->RemoveHadronEntry(itrk);
        n_hadrons--;
        itrk--;
      }
    }
    for (int itrk = 0; itrk <  event->GetNhadron(); itrk++)
    {
      MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
      if (remove_hadron_hits && (mytrk->GetPtPrime() < 1.5 || TMath::Abs(mytrk->GetPC3SDPHI())>2 || TMath::Abs(mytrk->GetPC3SDZ())>2)) 
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
    if (check_veto)
      event_container->CheckVeto();
    if (fill_true_DCA)
      event_container->FillTrueDCA();
    if (fill_TTree)
      event_container->FillTree();
    if(fill_inv_mass)
      event_container->fill_inv_mass();

  }
  event_container->WriteOutFile();
}
