#include "InvMass.h"
#include "TMath.h"
#include "TF1.h"

double Tsallis(const double pt, const double mass = 0.135, const double n = 8, const double T = 0.15, const double dNdy = 1)
{
  if (n < 1) return 1;
  const double mT = sqrt(pt*pt+mass*mass);
  return  dNdy * pt * (n - 1) * (n - 2) / (n*T + mass*(n - 1)) / (n*T + mass) * pow ( (n*T + mT) / (n*T + mass) , -n) ; //1./(2*TMath::Pi())
}

double Hagedorn(const double pt, const double dNdy = 95.7, const double mass = 0.135, const double A = 504.5, const double a = 0.5169, const double b = 0.1626, const double p0 = 0.7366, const double n = 8.274)
{
  if (mass < 0) return 1;
  const double mT = sqrt(pt*pt + mass*mass - 0.018218683); ///mpi0^2
  return  dNdy / dNdy_pp[0] * pt * A * pow( exp(-a * mT - b *mT * mT ) + mT / p0, -n) ;///2 * TMath::Pi() *
}
Double_t Hagedorn_Yield_Func(Double_t *x, Double_t *p) {
  //
  // Hagedorn function in 1/pt dN/pt from ppg088 
  //
    Double_t pt   = x[0];
    Double_t mass = p[0];
    Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.13498*0.13498);
    Double_t A    = p[1];
    Double_t a    = p[2];
    Double_t b    = p[3];
    Double_t p0   = p[4];
    Double_t n    = p[5];
   
    Double_t value = 2 * TMath::Pi() * pt*A* pow( exp(-a*mt-b*mt*mt)+mt/p0 , -n);
    return value;
  }
  
TF1 *fHagedorn(const Char_t *name, Double_t lowerlim,  Double_t upperlim,  const double dNdy = 95.7, const double mass = 0.135, const double A = 504.5, const double a = 0.5169, const double b = 0.1626, const double p0 = 0.7366, const double n = 8.274)
{
  TF1 *fHagedorn = new TF1(name, Hagedorn_Yield_Func, lowerlim, upperlim, 6);
  fHagedorn->SetParameters(mass, A, a, b, p0, n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
}

void InvMass(const TString inname = inFile[0],  int itread = 0, int ntreads = 1, int N_max = 1000000, const int part = 0)
{
  std::cout<<"start"<<std::endl;

  TF1* hadron_yield[5];
  for (int i = 0; i < 5; i++)
  {
    TString name = Form("hagedorr_yield_%d", i);
    hadron_yield[i] = fHagedorn(name, 0.01, 20, dNdy_pp[part], m_pp[part]);
    std::cout<<" integral " << i << " " << hadron_yield[i]->Integral(0.01,20) <<std::endl;
    //hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) / hadron_yield[i]->Integral(0.4,10) * dNdy_pp[part]/42.2*Ncolls[i] * Br_to_ee[part] );
    hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) * Ncolls[i] / 257. * Br_to_ee[part] );
    if (part == 8 || part == 9 ) hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) / hadron_yield[i]->Integral(0.01,20) * dNdy_pp[part] / 42.2 * (Ncolls[i] / 257. * Br_to_ee[part]));
    if (part < 8 && part >0)
      hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) * hith_pt_ratio[part] / (hadron_yield[i]->Eval(5.0) / Br_to_ee[part] / pi0_high_pt_ratio[i] ) );
    std::cout<<" integral " << i << " " << hadron_yield[i]->Integral(0.01,20)/Br_to_ee[part] << " " << dNdy_pp[part]/42.2*Ncolls[i] <<" "<<Br_to_ee[part] << " " << hadron_yield[i]->Eval(5.0)/Br_to_ee[part]<<std::endl;
  }

  std::cout<<"following particle was chosen: "<< part_names[part] <<std::endl;
  TFile *input = new TFile(inname, "READ");
  if (!(input))
  {
    std::cout << "no input file" << std::endl;
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
  TH2D *hist_pt_mother = new TH2D("hist_pt_mother","hist_pt_mother",500,0,10,5,0,100);
  TH2D *hist_pt_mother_weight = new TH2D("hist_pt_mother_weight","hist_pt_mother_weight",500,0,10,5,0,100);
  TH2D *n_hits_hist = new TH2D("n_hits_hist","n_hits_hist",1000,0,1000,100,0,100); 

  std::cout << "Trees read!" << std::endl;

  const int nevt = T->GetEntries();
  const int beggin = nevt * (itread - 1) / ntreads;
  const int endish = nevt * itread / ntreads;
  
  TF1 f("f","gaus",-5,5);
  f.SetParameter(0,1);
  f.SetParameter(1,0);
  f.SetParameter(2,1);
  const float event_vertex_sigma[3] = {89./10000,76./10000,100./10000};

  for (int ievent = beggin; ievent < endish; ievent++)
  {
    if ((ievent -beggin) % 50000 == 0)
    std::cout << "ithread, iEvent, N_events: " << itread<< ",  " << ievent -beggin<< " / " << nevt/ntreads << std::endl;
    myevent->ClearEvent();
    event_container->ClearEvent();
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
    double weight = Ncolls[(int)myevent->GetCentrality()/20];
    if (myevent->GetNgentrack()==2)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk1 = myevent->GetGenTrack(0);
      MyDileptonAnalysis::MyGenTrack *mygentrk2 = myevent->GetGenTrack(1);
      const double px = mygentrk1->GetPx()+mygentrk2->GetPx();
      const double py = mygentrk1->GetPy()+mygentrk2->GetPy();
      const double pt = sqrt(px*px+py*py);
      //weight = Tsallis(pt, m_pp[part], n_pp[part], T_pp[part], dNdy_pp[part]/42.2*Ncolls[(int)myevent->GetCentrality()/20]);
      //weight = Hagedorn(pt, dNdy_pp[part]*Br_to_ee[part], m_pp[part]);
      if ( part < 10) weight = hadron_yield[(int)myevent->GetCentrality()/20]->Eval(pt);
      hist_pt_mother->Fill(pt,myevent->GetCentrality());
      hist_pt_mother_weight->Fill(pt,myevent->GetCentrality(),weight);
    }
    if (myevent->GetNgentrack()==3)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk1 = myevent->GetGenTrack(0);
      MyDileptonAnalysis::MyGenTrack *mygentrk2 = myevent->GetGenTrack(1);
      MyDileptonAnalysis::MyGenTrack *mygentrk3 = myevent->GetGenTrack(2);
      const double px = mygentrk1->GetPx()+mygentrk2->GetPx()+mygentrk3->GetPx();
      const double py = mygentrk1->GetPy()+mygentrk2->GetPy()+mygentrk3->GetPy();
      const double pt = sqrt(px*px+py*py);
      //weight = Hagedorn(pt, dNdy_pp[part]*Br_to_ee[part], m_pp[part]);
      //weight = Tsallis(pt, m_pp[part], n_pp[part], T_pp[part], dNdy_pp[part]/42.2*Ncolls[(int)myevent->GetCentrality()/20]);
      if ( part < 10) weight = hadron_yield[(int)myevent->GetCentrality()/20]->Eval(pt);
      hist_pt_mother->Fill(pt,myevent->GetCentrality());
      hist_pt_mother_weight->Fill(pt,myevent->GetCentrality(),weight);
    }

    myevent->SetPreciseX(myevent->GetPreciseX()+f.GetRandom()*event_vertex_sigma[0]);
    myevent->SetPreciseY(myevent->GetPreciseY()+f.GetRandom()*event_vertex_sigma[1]);
    myevent->SetPreciseZ(myevent->GetPreciseZ()+f.GetRandom()*event_vertex_sigma[2]);
    myevent->SetVtxZ(myevent->GetPreciseZ());
    //if(myevent->GetCentrality()>20) continue;

    if(fill_inv_mass_sim) event_container->fill_inv_mass_sim(weight);

    if(myevent->GetNtrack()<1 && !do_track_QA) continue;
    
    //for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    //{
    //  MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
    //  mytrk->SetPtPrime(mytrk->GetPt());
    //}

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
        //std::cout<<mytrk.GetPt()<<" "<<mytrk.GetPtPrime()<<" "<<mytrk.GetReconPT()<<std::endl;
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
    n_hits_hist->Fill(n_hits,myevent->GetCentrality());
    
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
    if(fill_true_DCA) event_container->FillTrueDCA(weight);
    if(fill_TTree) event_container->FillTree();
    if(do_electron_QA)event_container->FillQAHist(-99);
    //myevent->ReshuffleElectrons();
    if(fill_inv_mass && event_container->GetNGoodElectrons()<2  ) continue;
    if(fill_inv_mass) event_container->fill_inv_mass(weight);
    myevent->ClearEvent();
    event_container->ClearEvent();
    //event->ClearEvent();
  }
  hist_pt_orig->Write();
  hist_pt_mother->Write();
  hist_pt_mother_weight->Write();
  n_hits_hist->Write();
  event_container->WriteOutFile();
}
