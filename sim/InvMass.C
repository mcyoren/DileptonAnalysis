#include "InvMass.h"
void InvMass(const TString inname = inFile[0],  int itread = 0, int ntreads = 1, int N_max = 1000000)
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
  const int do_electron_QA = 1;
  const int do_ident_electrons = 1;
  const int do_reveal_hadron = 0;
  const int Use_ident = 0;
  const int fill_true_DCA = 1;
  const int check_veto = 0;
  const int fill_inv_mass = 1;
  const int fill_inv_mass_sim = 1;

  char outname[200];
  sprintf(outname,"kek_%d.root",itread);
  TTree *T = (TTree *)input->Get("tree");
  TBranch *br = T->GetBranch("MyEvent");
  MyDileptonAnalysis::MyEventContainer *event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile("../../ee_QA/AnaTrain/Run14AuAuLeptonComby/field_map.root");
  event_container->CreateOutFileAndInitHists(outname,fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,
                                               fill_DCA_hists, do_track_QA+do_electron_QA, do_reveal_hadron, fill_true_DCA, check_veto, fill_inv_mass);
  MyDileptonAnalysis::MyEvent *myevent = 0;                                            
  //event = 0;
  br->SetAddress(&myevent);

  TH3D *hist_pt_orig = new TH3D("hist_pt_orig","hist_pt_orig",50,0,5,5,0,100,2,0,2);

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
  TF1 f("f","gaus",-5.*250./10000,5.*250./10000);
  f.SetParameter(0,1);
  f.SetParameter(1,0);
  f.SetParameter(2,250./10000);
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
      if (mygentrk->GetID()==22) continue;
      if (mygentrk->GetID()>0)
        hist_pt_orig->Fill(mygentrk->GetPt(),myevent->GetCentrality(),0);
      else
        hist_pt_orig->Fill(mygentrk->GetPt(),myevent->GetCentrality(),1);
    }
    myevent->SetPreciseX(myevent->GetPreciseX()+f.GetRandom()*0);
    myevent->SetPreciseY(myevent->GetPreciseY()+f.GetRandom()*0);
    myevent->SetVtxZ(myevent->GetPreciseZ());

    if(fill_inv_mass_sim) event_container->fill_inv_mass_sim();

    if(myevent->GetNtrack()<1 && !do_track_QA) continue;

    int n_electrons = myevent->GetNtrack();
    for (int itrk = 0; itrk < n_electrons; itrk++)
    {
      MyDileptonAnalysis::MyElectron mytrk = *myevent->GetEntry(itrk);
      
      bool skip = false; 
        if (mytrk.GetPtPrime()>4.4 || mytrk.GetPtPrime() < 0.4)
          skip = true;
        if (mytrk.GetCrkphi()<-99)
          skip = true;
        if (mytrk.GetEcore()/mytrk.GetPtot() < 0.6)
          skip = true;
        if (mytrk.GetN0() < 0 )
          skip = true;
        if (fabs(mytrk.GetEmcdphi())>0.05 || fabs(mytrk.GetEmcdz())>25 )
          skip = true;

        if( skip ){
          myevent->RemoveTrackEntry(itrk);
          //event->AddElecCand(&mytrk);
          n_electrons--;
          itrk--;
          continue;
        }
    }
    
    int n_hits = myevent->GetNVTXhit();
    for (int ihit = 0; ihit < n_hits; ihit++)
    {
      MyDileptonAnalysis::MyVTXHit myhit = *myevent->GetVTXHitEntry(ihit);
    if(myhit.GetLadder()>24)std::cout<<myhit.GetLadder()<<std::endl;
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
    for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    {
        MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
        mytrk->ZeroHitCounters();
    }
    if(fill_QA_hadron_hists) event_container->Associate_Hits_to_Hadrons(400);
    if(do_track_QA) event_container->FillQAHist(in_id);

    if(do_ident_electrons) event_container->IdenElectrons();
    if(do_electron_QA)event_container->FillQAHistPreAssoc();

    if(associate_hits) event_container->Associate_Hits_to_Leptons(5,5,5,0,2);
    if(associate_hits) event_container->Associate_Hits_to_Leptons(5,5,5,0,1);
    if(associate_hits) event_container->Associate_Hits_to_Leptons(5,5,5,0,1);
    if(false) event_container->Associate_Hits_to_Leptons_OLD(20,20,20);

    if(associate_hits && event_container->GetNGoodElectrons()<1  ) continue;

    //int n_electrons = myevent->GetNtrack()*0;
    //for (int itrk = 0; itrk < n_electrons; itrk++)
    //{
    //    MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
    //    if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
    //       (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 ))
    //    {
    //        myevent->RemoveTrackEntry(itrk);
    //        n_electrons--;
    //        itrk--;
    //        continue;
    //    }
    //    mytrk->ZeroHitCounters();
    //}

    if(false)
        {

            for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
            {
                MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
                if(event_container->GetNGoodElectrons()>1) std::cout<< "000: " <<mytrk->GetNHits()<<" "<<mytrk->GetTOFDPHI()<<" "<<mytrk->GetGhost()<<" "<<mytrk->GetHitCounter(0)<<std::endl;
            }
            for (int ibdtrack = 0; ibdtrack < (int) event_container->GetNBDThit(); ibdtrack++)
            {       
                if(event_container->GetNGoodElectrons()<1) continue;
                std::cout<<event_container->GetBDTHitEntry(ibdtrack)->GetPt()<<" "<<event_container->GetBDTHitEntry(ibdtrack)->GetEcore()<<" "<<event_container->GetBDTHitEntry(ibdtrack)->GetNBDThit()<<std::endl;
                if(event_container->GetBDTHitEntry(ibdtrack)->GetPt()<0.4) continue;
                MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(ibdtrack);
                std::cout<<mytrk->GetHitIndex(0)<<" "<<mytrk->GetHitIndex(1)<<" "<<mytrk->GetHitIndex(2)<<" "<<mytrk->GetHitIndex(3)<<std::endl;
                for (int jbdthit = 0; jbdthit < (int) event_container->GetBDTHitEntry(ibdtrack)->GetNBDThit(); jbdthit++)
                {
                    std::cout<<event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(0,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(1,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(2,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(3,1)<<" "<<std::endl;
                }
            }
        }

        if(associate_hits && event_container->GetNGoodElectrons()>1) event_container->Associate_Hits_to_Leptons(5,5,5);
    
    if(false)
        {

            for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
            {
                MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
                if(event_container->GetNGoodElectrons()>1) std::cout<< "1: " <<mytrk->GetNHits()<<" "<<mytrk->GetTOFDPHI()<<" "<<mytrk->GetGhost()<<" "<<mytrk->GetHitCounter(0)<<std::endl;
            }
            for (int ibdtrack = 0; ibdtrack < (int) event_container->GetNBDThit(); ibdtrack++)
            {       
                if(event_container->GetNGoodElectrons()<1) continue;
                std::cout<<event_container->GetBDTHitEntry(ibdtrack)->GetPt()<<" "<<event_container->GetBDTHitEntry(ibdtrack)->GetEcore()<<" "<<event_container->GetBDTHitEntry(ibdtrack)->GetNBDThit()<<std::endl;
                if(event_container->GetBDTHitEntry(ibdtrack)->GetPt()<0.4) continue;
                MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(ibdtrack);
                std::cout<<mytrk->GetHitIndex(0)<<" "<<mytrk->GetHitIndex(1)<<" "<<mytrk->GetHitIndex(2)<<" "<<mytrk->GetHitIndex(3)<<std::endl;
                for (int jbdthit = 0; jbdthit < (int) event_container->GetBDTHitEntry(ibdtrack)->GetNBDThit(); jbdthit++)
                {
                    std::cout<<event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetIsTrue(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiL(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(0)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(2)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(3)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(0,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(1,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(2,1)<<" "<<
                    event_container->GetBDTHitEntry(ibdtrack)->GetBDTHitEntry(jbdthit)->GetSecondHitPhiR(3,1)<<" "<<std::endl;
                }
            }
        }
    
    if(remove_hadron_hits) 
    {
        event_container->Associate_Hits_to_Hadrons();
        event_container->Associate_Hits_to_Leptons();
    }

    n_electrons = myevent->GetNtrack();
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
    }
    if(associate_hits && event_container->GetNGoodElectrons()>1) event_container->ResetRecoverFGVars(); 

    if(use_d_dphi_DCA)  event_container->FillDphiHists();
    if(do_reveal_hadron) event_container->Reveal_Hadron();
    //event_container->CheckVeto();
    if(fill_true_DCA) event_container->FillTrueDCA();
    if(fill_TTree) event_container->FillTree();
    if(do_electron_QA)event_container->FillQAHist(-99);
    //myevent->ReshuffleElectrons();
    if(fill_inv_mass && event_container->GetNGoodElectrons()<2  ) continue;
    if(fill_inv_mass) event_container->fill_inv_mass();
    myevent->ClearEvent();
    //event->ClearEvent();
  }
  hist_pt_orig->Write();
  event_container->WriteOutFile();
}
