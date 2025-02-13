#include "Calib.h"
void Testm2(const TString inname = inFile[0],  int itread = 0, int ntreads = 1, int N_max = 1000000)
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
  const int fill_true_DCA = 1;
  const int check_veto = 0;
  const int fill_inv_mass = 0;
  const int istruehitsigmacounter = 0;


  TH3D *hist_m2 = new TH3D("hist_m2","hist_m2",2000,0,2,50,0,5,5,0,100);
  TH3D *hist_m2_emb = new TH3D("hist_m2_emb","hist_m2_emb",2000,-10,40,50,0,5,5,0,100);
  TH3D *hist_m2_ass = new TH3D("hist_m2_ass","hist_m2_ass",2000,-10,40,50,0,5,5,0,100);
  TH2D *hist_pt_orig = new TH2D("hist_pt_orig","hist_pt_orig",50,0,5,5,0,100);
  TH2D *hist_pt_calc = new TH2D("hist_pt_calc","hist_pt_calc",50,0,5,5,0,100);
  TH2D *hist_pt_embed = new TH2D("hist_pt_embed","hist_pt_embed",50,0,5,5,0,100);
  TH2D *hist_pt_vtx = new TH2D("hist_pt_vtx","hist_pt_vtx",50,0,5,5,0,100);


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
  TF1 f("f","gaus",-5.*200./10000,5.*200./10000);
  f.SetParameter(0,1);
  f.SetParameter(1,0);
  f.SetParameter(2,1000./10000);
  const float deltaR[3]={0.080*0,0.016*0,0.200};
 
  for (int ievent = beggin; ievent < endish; ievent++)
  {
    if ((ievent -beggin) % 50000 == 0)
      cout << "ithread, iEvent, N_events: " << itread<< ",  " << ievent -beggin<< " / " << nevt/ntreads << endl;
    myevent->ClearEvent();
    br->GetEntry(ievent);
    if (ievent - beggin > N_max)
      break;

    for (int itrk = 0; itrk < myevent->GetNgentrack(); itrk++)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk = myevent->GetGenTrack(itrk);
      hist_pt_orig->Fill(mygentrk->GetPt(),myevent->GetCentrality());
    }
    myevent->SetVtxZ(myevent->GetPreciseZ());
    myevent->SetPreciseX(myevent->GetPreciseX()+f.GetRandom()*deltaR[0]);
    myevent->SetPreciseY(myevent->GetPreciseY()+f.GetRandom()*deltaR[1]);
    myevent->SetPreciseZ(myevent->GetPreciseZ()+f.GetRandom()*deltaR[2]);



    int n_hits = myevent->GetNVTXhit();
    for (int ihit = 0; ihit < n_hits; ihit++)
    {
      MyDileptonAnalysis::MyVTXHit myhit = *myevent->GetVTXHitEntry(ihit);

     if (myhit.GetSensor() == 1)//no sim hits
     {
         myevent->RemoveVTXHitEntry(ihit);
         n_hits--;
         ihit--;
         continue;
     }
    }
    
  
    event_container->SetEvent(myevent);
    
    if(fill_QA_hadron_hists) event_container->correct_beam_offset();
    for (int itrk = 0; itrk < myevent->GetNeleccand(); itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = myevent->GetElecCand(itrk);
      hist_m2->Fill(mytrk->GetTOFE(),mytrk->GetPt(),myevent->GetCentrality());
      if (mytrk->GetTOFE()>0.4 && mytrk->GetTOFE()<1.2)
        hist_pt_calc->Fill(mytrk->GetPt(),myevent->GetCentrality());
      if (mytrk->GetTOFE()>0.4 && mytrk->GetTOFE()<1.2 && mytrk->GetEmcId()>=0)
        hist_pt_embed->Fill(mytrk->GetPt(),myevent->GetCentrality());
    }
    event_container->Associate_Hits_to_Leptons(5,5,5,0,2);
    event_container->Associate_Hits_to_Leptons(5,5,5,0,1);
    for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
      hist_m2_emb->Fill(mytrk->GetTOFE(),mytrk->GetPtPrime(),myevent->GetCentrality());
      //if (mytrk->GetTOFE()>0.4 && mytrk->GetTOFE()<1.2)
      //  hist_pt_embed->Fill(mytrk->GetPtPrime(),myevent->GetCentrality());
      if (mytrk->GetNHits()>1)
      {
        hist_m2_ass->Fill(mytrk->GetTOFE(),mytrk->GetPt(),myevent->GetCentrality());
        if (mytrk->GetTOFE()>0.4 && mytrk->GetTOFE()<1.2)
          hist_pt_vtx->Fill(mytrk->GetPtPrime(),myevent->GetCentrality());
      }
    }
    continue;
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
    if(check_veto||fill_inv_mass) event_container->CheckVeto();
    if(fill_true_DCA) event_container->FillTrueDCA();
    if(fill_TTree) event_container->FillTree();
    myevent->ReshuffleElectrons();
    if(fill_inv_mass) event_container->fill_inv_mass();
    
  
  }
  hist_m2->Write();
  hist_m2_emb->Write();
  hist_m2_ass->Write();
  hist_pt_orig->Write();
  hist_pt_calc->Write();
  hist_pt_embed->Write();
  hist_pt_vtx->Write();
  event_container->WriteOutFile();

}
