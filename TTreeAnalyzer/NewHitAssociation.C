#include "NewHitAssociation.h"
void NewHitAssociation(int par = 0)
{

  TFile *input = new TFile(inFile[par], "READ");
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
  const int Use_ident = 1;
  const int fill_true_DCA = 1;
  const int check_veto = 1;
  const int is_only_conv = 0;

  TTree *T = (TTree *)input->Get("tree");
  TBranch *br = T->GetBranch("MyEvent");
  MyDileptonAnalysis::MyEventContainer *event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile("../AnaTrain/Run14AuAuLeptonComby/field_map.root");
  event_container->CreateOutFileAndInitHists("kek.root", fill_QA_lepton_hists, fill_QA_hadron_hists, fill_TTree, fill_d_dphi_hists,
                                             fill_DCA_hists, do_track_QA, do_reveal_hadron, fill_true_DCA, check_veto);
  MyDileptonAnalysis::MyEvent *event = event_container->GetEvent();
  // event = 0;
  br->SetAddress(&event);

  cout << "Trees read!" << endl;

  int nevt = T->GetEntries();

  TH3D *myhist0 = new TH3D("myhist0", "myhist0", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist1 = new TH3D("myhist1", "myhist1", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist2 = new TH3D("myhist2", "myhist2", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist3 = new TH3D("myhist3", "myhist3", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist4 = new TH3D("myhist4", "myhist4", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist5 = new TH3D("myhist5", "myhist5", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist6 = new TH3D("myhist6", "myhist6", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);
  TH3D *myhist7 = new TH3D("myhist7", "myhist7", 500, -2000, 2000, 500, 0, 5, 10, 0, 10);

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
    for (int itrk = 0; itrk < n_electrons*is_only_conv; itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
      if (mytrk->GetIsConv() < 4) 
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

    int n_hadrons = event->GetNhadron();
    for (int itrk = 0; itrk < n_hadrons*remove_hadron_hits; itrk++)
    {
      MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
      if (mytrk->GetPtPrime() < 1.0 || TMath::Abs(mytrk->GetPC3SDPHI())>2 || TMath::Abs(mytrk->GetPC3SDZ())>2) 
      {
        event->RemoveHadronEntry(itrk);
        n_hadrons--;
        itrk--;
      }
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

    float x0 = event->GetPreciseX();
    const float y0 = event->GetPreciseY();
    for (int i = 0; i < event->GetNtrack(); i++)
    {
      MyDileptonAnalysis::MyElectron *newTrack1 = event->GetEntry(i);
      if (newTrack1->GetGhost() > 0)
        continue;
      if (newTrack1->GetChargePrime() < 0)
        continue;
      const float a1 = newTrack1->GetMinDist(0);
      const float b1 = newTrack1->GetMinDist(1);
      const float c1 = newTrack1->GetMinDist(2);
      for (int j = 0; j < event->GetNtrack(); j++)
      {
        MyDileptonAnalysis::MyElectron *newTrack2 = event->GetEntry(j);
        if (newTrack2->GetChargePrime() > 0)
          continue;
        if (newTrack2->GetGhost() > 0)
          continue;
        const float a2 = newTrack2->GetMinDist(0);
        const float b2 = newTrack2->GetMinDist(1);
        const float c2 = newTrack2->GetMinDist(2);
        const float arm1 = (0.5 - newTrack1->GetArm()) * 2;
        const float arm2 = (0.5 - newTrack2->GetArm()) * 2;
        float dca0 = ((a1 * x0 * x0 + b1 * x0 + c1) - (a2 * x0 * x0 + b2 * x0 + c2)) * 10000;
        // x0 = (b1*arm1 + b2*arm2)/2/((a1*arm1 + a2*arm2));
        /*dca0 = 99999.;
        if (dca0>0)
        {
            if((b1-b2)*(b1-b2)>4*(a1-a2)*(c1-c2))
            {
              float x00 = (-(b1-b2)+sqrt((b1-b2)*(b1-b2)-4*(a1-a2)*(c1-c2)))/2/(a1-a2);
              float x01 = (-(b1-b2)-sqrt((b1-b2)*(b1-b2)-4*(a1-a2)*(c1-c2)))/2/(a1-a2);
              if(SQR(x00-x0) + SQR(a1*x00*x00+b1*x00+c1-y0)<SQR(x01-x0) + SQR(a1*x01*x01+b1*x01+c1-y0))
              {
                    dca0 = sqrt( SQR(x00-x0) + SQR(a1*x00*x00+b1*x00+c1-y0))*10000;
              }else dca0 = sqrt( SQR(x01-x0) + SQR(a1*x01*x01+b1*x01+c1-y0))*10000;
            }else std::cout<<"WTF"<<std::endl;
            const float xmin = abs(-(b1-b2)/2/(a1-a2))*10000;
            if (xmin<dca0) dca0=xmin;
        }else dca0 = -99999;
        dca0=99999;
        for (int iphi = 0; iphi < 10000; iphi++)
        {
            float phi1 = newTrack1->GetPhi0Prime()+(0.00001*iphi-0.05);
            if(newTrack1->GetPx()<0) phi1-=pi;
            float x1 =(tan(phi1)-b1)/2/a1;
            float phi2 = newTrack2->GetPhi0Prime()+(0.00001*iphi-0.05);
            if(newTrack1->GetPx()<0) phi2-=pi;
            float x2 =(tan(phi2)-b2)/2/a2;
            float dca00 = ((a1*x1*x1+b1*x1+c1) - (a2*x2*x2+b2*x2+c2))*10000;
            if(abs(dca00)<abs(dca0))dca0 = dca00;
        }
        */
        dca0 = abs(newTrack1->GetDCAX2() + newTrack2->GetDCAX2());
        const float dca1 = abs(newTrack1->GetDCAY2() + newTrack2->GetDCAY2());
        const float dca2 = abs(newTrack1->GetDCAX2() - newTrack2->GetDCAX2());
        const float dca3 = abs(newTrack1->GetDCAY2() - newTrack2->GetDCAY2());
        const float dca4 = abs(newTrack1->GetDCAX2() / abs(newTrack1->GetDCAX2()) * abs(newTrack1->GetDCA2()) + newTrack2->GetDCAX2() / abs(newTrack2->GetDCAX2()) * abs(newTrack2->GetDCA2()));
        const float dca5 = abs(newTrack1->GetDCAY2() / abs(newTrack1->GetDCAY2()) * abs(newTrack1->GetDCA2()) + newTrack2->GetDCAY2() / abs(newTrack2->GetDCAY2()) * abs(newTrack2->GetDCA2()));
        const float dca6 = abs(newTrack1->GetDCAX2() / abs(newTrack1->GetDCAX2()) * abs(newTrack1->GetDCA2()) - newTrack2->GetDCAX2() / abs(newTrack2->GetDCAX2()) * abs(newTrack2->GetDCA2()));
        const float dca7 = abs(newTrack1->GetDCAY2() / abs(newTrack1->GetDCAY2()) * abs(newTrack1->GetDCA2()) - newTrack2->GetDCAY2() / abs(newTrack2->GetDCAY2()) * abs(newTrack2->GetDCA2()));

        const float pair_pt = sqrt((newTrack1->GetPx() + newTrack2->GetPx()) * (newTrack1->GetPx() + newTrack2->GetPx()) + (newTrack1->GetPy() + newTrack2->GetPy()) * (newTrack1->GetPy() + newTrack2->GetPy()));

        const float px1 = newTrack1->GetPx();
        const float py1 = newTrack1->GetPy();
        const float pz1 = newTrack1->GetPz();
        const float px2 = newTrack2->GetPx();
        const float py2 = newTrack2->GetPy();
        const float pz2 = newTrack2->GetPz();
        const float pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
        const float pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
        const float es = sqrt(pm1 + me2) + sqrt(pm2 + me2);
        const float px = px1 + px2;
        const float py = py1 + py2;
        const float pz = pz1 + pz2;

        const float invm = sqrt(es * es - px * px - py * py - pz * pz);
        // std::cout<<dca*10000<<" "<<event->GetPreciseX()<<std::endl;
        myhist0->Fill(dca0, invm, pair_pt);
        myhist1->Fill(dca1, invm, pair_pt);
        myhist2->Fill(dca2, invm, pair_pt);
        myhist3->Fill(dca3, invm, pair_pt);
        myhist4->Fill(dca4, invm, pair_pt);
        myhist5->Fill(dca5, invm, pair_pt);
        myhist6->Fill(dca6, invm, pair_pt);
        myhist7->Fill(dca7, invm, pair_pt);

        // std::cout<<newTrack1->GetDCA2()<<" "<<newTrack2->GetDCA2()<<std::endl;
        // std::cout<<"grad: "<<newTrack1->GetPhiDC()*180/3.14159<<" "<<newTrack2->GetPhiDC()*180/3.14159<<std::endl;
        // const float ptpair = newTrack1->GetPy()+newTrack2->GetPy();
        // if(newTrack1->GetPhiDC()<0.8) std::cout<<newTrack2->GetCrkphi()*180/3.14159-newTrack1->GetCrkphi()*180/3.14159<<" "<<ptpair<<std::endl;
        // else std::cout<<360-newTrack2->GetCrkphi()*180/3.14159-newTrack1->GetCrkphi()*180/3.14159<<" kek "<<ptpair<<std::endl;
      }
    }
  }
  TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
  c1->Draw();
  myhist0->Draw();
  c1->SaveAs("output/pics/dca.png");
  event_container->WriteOutFile();
}
