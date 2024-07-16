#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
// 
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include "HELIOSLibrary/HELIOSLibrary.h"
#include "MyPlotting/MyPlot.C"


using namespace std;
Double_t f_EcalReso(Double_t *x,Double_t *p);

///////////////////////////////////////////////////////////////////////////////
void TestEMCalReso(){
//
  TRandom3 MyRandy = TRandom3(0);                // Random Generator

  cout << endl;

  PHENIXDetector PHENIX;
  Particle pi0("pi0");

  TLorentzVector gamma1, gamma2; 
  TLorentzVector Reco0, Reco1, Reco2;

  Double_t pt0,phi0,eta0,ww0;
  Double_t pt1,phi1,eta1,ww1;
  Double_t pt2,phi2,eta2,ww2;
  Int_t q1,q2;
  Int_t is1,is2;
  double y,z;
   
  Double_t pt_min = 0.1;
  Double_t pt_max = 20.;
  Double_t pt_low = 0.1;
  Double_t pt_high = 15.;
  Double_t AssyCut = 0.2;                      // Assymetry cut for calibration
  Double_t Ecut = 0.2;
  Int_t nevt = 100000;

  Bool_t g1 = false;
  Bool_t g2 = false;
  Bool_t isPbGl1 = false;
  Bool_t isPbGl2 = false;

  Int_t ndecay; 
  Double_t ww;

//  TF1 *piHagedorn      = Hagedorn("piHagedorn", pi0Mass, pt_max, pt_min);
  TF1 *piHagedornRE    = HagedornYield("piHagedornRE", pi0Mass, pt_max, pt_min, 58, 0.661, 0.015, 0.745, -9.167);
// ppg088
  TF1 *piHagedorn     = HagedornYield("piHagedorn", pi0Mass, pt_max, pt_min, 377., 0.356, 0.068, 0.7, -8.25);

  TH1D *h_ptpi             = new TH1D("h_ptpi0","pt",100,pt_low,pt_high);
  TH1D *h_ptpireco         = new TH1D("h_ptpireco","ptreco",100,pt_low,pt_high);
  TH2D *h_butsyk           = new TH2D("h_butsyk","butsyk",100,-2*pi,2*pi,100,-6.0,6.0);
  TH2D *h_Ereso            = new TH2D("h_Ereso","Ereso",150,pt_low,pt_max,250,-.3,.3);
  TH2D *h_piMass           = new TH2D("h_piMass","piMass",150,pt_low,pt_max,250,0.,1.3);
  TH2D *h_piMassPos        = new TH2D("h_piMassPos","piMassPos",150,pt_low,pt_max,250,0.,1.3);
  TH2D *h_Ereso1           = new TH2D("h_Ereso1","Ereso1",150,pt_low,pt_max,250,-.3,.3);
  TH2D *h_piMass1          = new TH2D("h_piMass1","piMass1",150,pt_low,pt_max,250,0.,1.3);
  TH2D *h_piMassPos1       = new TH2D("h_piMassPos1","piMassPos1",150,pt_low,pt_max,250,0.,1.3);


  TF1 *EMCalReso = new TF1("EMCalReso",f_EcalReso,0.,10.,2);
  EMCalReso->SetParameters(0.081,0.021); 
  TF1 *EMCalReso1 = new TF1("EMCalReso1",f_EcalReso,0.,10.,2);
  EMCalReso1->SetParameters(0.059,0.008); 


  for (int i=0; i<nevt; i++){                                     // event loop
    if(i%100000==0) cout << "event loop at " << i << endl;         // event counter printed 
/////////////////////////////////////////////// pi decays ////////////////////////////////////////////////
    g1 = false;
    g2 = false;

    pi0.ResetP();
    if  (i < nevt) {
      pi0.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi with flat distribution
      pi0.SetWeight(piHagedorn->Eval(pi0.Pt())/float(nevt));
    } 
    if (PHENIX.InAcceptance(pi0,pi0.Charge())>0) {                                // monitor PEHNIX acceptance 
      h_butsyk->Fill(pi0.Phi(),1./pi0.Pt(),pi0.Weight());                   //      for charged particles
      h_ptpi->Fill(pi0.Pt(),pi0.Weight());
      
      pi0.DecaySingleBranch("pi0->gg");
      ndecay = pi0.GetNumberOfDaughters();
      if (i<50)  cout << "----- Pi0Deacy --------------------------------------------------------" << endl; 

      for (Int_t j=0; j< ndecay; j++) {
        if (j==0) {
          ww     = pi0.GetDaughterWeight(j)*pi0.Weight();
          gamma1 = pi0.GetDecayDaughter(j);
          PHENIX.CharacterizeTrack(gamma1,0,photonID);
          Reco1  = PHENIX.ReconstructShower(gamma1,pi0.GetDaughterID(j));
          if (i<50) {
                cout << " recostructed   ";
                Reco1.Print();
                cout << endl;
          }
          if (Reco1.E()>Ecut) {
              g1 = true;
              is1 =  PHENIX.EMCalSectorCoordinates(gamma1.Phi(),gamma1.Eta(),y,z);
              isPbGl1 = (is1 == 7 or is1 == 8);
              if (isPbGl1) {
                h_Ereso1->Fill(Reco1.E(),(Reco1.E()-gamma1.E())/gamma1.E());
              } else {
                h_Ereso->Fill(Reco1.E(),(Reco1.E()-gamma1.E())/gamma1.E());               
              } 
          }
        } 
        if (j==1) {

          gamma2 = pi0.GetDecayDaughter(j);
          PHENIX.CharacterizeTrack(gamma2,0,photonID);
          Reco2  = PHENIX.ReconstructShower(gamma2,pi0.GetDaughterID(j));
          if (Reco2.E()>Ecut) {
              g2 = true;
              is2 =  PHENIX.EMCalSectorCoordinates(gamma2.Phi(),gamma2.Eta(),y,z);
              isPbGl2 = (is2 == 7 or is2 == 8);
              if (isPbGl2) {
                h_Ereso1->Fill(Reco2.E(),(Reco2.E()-gamma2.E())/gamma2.E());   
              } else {
                h_Ereso->Fill(Reco2.E(),(Reco2.E()-gamma2.E())/gamma2.E());         
              }
          }
          if (i<50) {
              cout << " recostructed   ";
              Reco2.Print();
              cout << endl;
          }
        }
      }

      if (g1 && g2) {
        if (Reco1.E() > Ecut && Reco2.E()>Ecut &&
             abs((Reco1.E()-Reco2.E())/(Reco1.E()+Reco2.E()))<=AssyCut ){      // both photons in same EMCal sector and within assymetry cut
          Reco0 = Reco1+Reco2;
          h_ptpireco->Fill(Reco0.Pt(),ww);
          if (is1 == is2) {
            if (isPbGl1) {
              h_piMass1->Fill((Reco1.E()+Reco2.E())/2,Reco0.M()/pi0Mass,ww);   
            } else {
              h_piMass->Fill((Reco1.E()+Reco2.E())/2,Reco0.M()/pi0Mass,ww);               
            }
          }
        }
      }
//
// isolate position resolution part only
//      
      if (g1 && g2) {
        phi1 = Reco1.Phi();
        eta1 = Reco1.Eta();
        pt1  = gamma1.E()*cos(eta1);                         // recalculate pt     
        Reco1.SetPtEtaPhiM(pt1,eta1,phi1,0);
        phi2 = Reco2.Phi();
        eta2 = Reco2.Eta();
        pt2  = gamma2.E()*cos(eta2);                         // recalculate pt     
        Reco2.SetPtEtaPhiM(pt2,eta2,phi2,0);
        Reco0 = Reco1+Reco2;
//        if (is1 == is2) {
            if (isPbGl1) {
              h_piMassPos1->Fill((Reco1.E()+Reco2.E())/2,Reco0.M()/pi0Mass,ww);   
            } else {
              h_piMassPos->Fill((Reco1.E()+Reco2.E())/2,Reco0.M()/pi0Mass,ww);               
            }
//          }
      }
    }  
  } // end of event loop

 ////////////////////////////////////////////////////////////////////////////

  MyPlot plot;

  h_Ereso->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_Ereso_mean = (TH1D*)gDirectory->Get("h_Ereso_1");
  h_Ereso_mean->SetName("h_Ereso_mean");
  plot.StyleMe(h_Ereso_mean,20,kTeal,.8);
  TH1D * h_Ereso_sigma = (TH1D*)gDirectory->Get("h_Ereso_2");
  h_Ereso_sigma->SetName("h_Ereso_sigma");
  plot.StyleMe(h_Ereso_sigma,20,kTeal,.8);
  plot.StyleMe(EMCalReso, 20,kTeal,.8);
  h_Ereso1->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt

  TH1D * h_Ereso1_mean = (TH1D*)gDirectory->Get("h_Ereso1_1");
  h_Ereso1_mean->SetName("h_Ereso1_mean");
  plot.StyleMe(h_Ereso1_mean,20,kBlue,.8);
  TH1D * h_Ereso1_sigma = (TH1D*)gDirectory->Get("h_Ereso1_2");
  h_Ereso1_sigma->SetName("h_Ereso1_sigma");
  plot.StyleMe(h_Ereso1_sigma,20,kBlue,.8);
  plot.StyleMe(EMCalReso1, 20,kBlue,.8);


  h_piMass->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_piMass_mean = (TH1D*)gDirectory->Get("h_piMass_1");
  h_piMass_mean->SetName("h_piMass_mean");
  plot.StyleMe(h_piMass_mean,20,kOrange+2,.8);
  TH1D * h_piMass_sigma = (TH1D*)gDirectory->Get("h_piMass_2");
  h_piMass_sigma->SetName("h_piMass_sigma");
  plot.StyleMe(h_piMass_sigma,20,kOrange+2,.8);

  h_piMass1->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_piMass1_mean = (TH1D*)gDirectory->Get("h_piMass1_1");
  h_piMass1_mean->SetName("h_piMass1_mean");
  plot.StyleMe(h_piMass1_mean,20,kRed,.8);
  TH1D * h_piMass1_sigma = (TH1D*)gDirectory->Get("h_piMass1_2");
  h_piMass1_sigma->SetName("h_piMass1_sigma");
  plot.StyleMe(h_piMass1_sigma,20,kRed,.8);


  h_piMassPos->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_piMassPos_mean = (TH1D*)gDirectory->Get("h_piMassPos_1");
  h_piMassPos_mean->SetName("h_piMassPos_mean");
  plot.StyleMe(h_piMassPos_mean,20,kGreen+2,.8);
  TH1D * h_piMassPos_sigma = (TH1D*)gDirectory->Get("h_piMassPos_2");
  h_piMassPos_sigma->SetName("h_piMassPos_sigma");
  plot.StyleMe(h_piMassPos_sigma,20,kGreen+2,.8);

  h_piMassPos1->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_piMassPos1_mean = (TH1D*)gDirectory->Get("h_piMassPos1_1");
  h_piMassPos1_mean->SetName("h_piMassPos1_mean");
  plot.StyleMe(h_piMassPos1_mean,20,kGreen,.8);
  TH1D * h_piMassPos1_sigma = (TH1D*)gDirectory->Get("h_piMassPos1_2");
  h_piMassPos1_sigma->SetName("h_piMassPos1_sigma");
  plot.StyleMe(h_piMassPos1_sigma,20,kGreen,.8);

  plot.StyleMe(h_ptpi,   20, kBlue, 1., 1, 2); 
  plot.StyleMe(h_ptpireco, 20, kTeal, 1., 1, 2); 

  plot.SetLeftMargin(.2);  
  plot.SetyTitleOffset(1.4);
  TCanvas *c1 = plot.Canvas ("c1",400,500,10,10,1);
  TH1D *f1    = plot.Frame("f1","p_{T} [GeV/c^{2}]","dN/dp_{T} [a.u.]",pt_low,pt_high*.999,1e-12,1.);
  f1->Draw();
  h_ptpi->Draw("sameL Chist");
  h_ptpireco->Draw("sameL Chist");
  TLegend *L1 = plot.Legend(" " ,0.5,.79,.7,.95);
  L1->AddEntry(h_ptpi,"#pi^{0}","l");
  L1->AddEntry(h_ptpireco,"reconstructed #pi^{0}","l");
  L1->Draw("same");
  plot.Reset();
 
  TCanvas *c2 = plot.Canvas ("c2",400,400,10,510);
  TH1D *f2 = plot.Frame("c2","E_{true} [GeV]","(E_{reco}-E_{true})/E_{true}",0.,15.,-0.3,.3);
  f2->Draw();
  h_Ereso->Draw("Col same");
  h_Ereso_mean->Draw("same");
  TLegend *L2 = plot.Legend(" PbSc " ,0.5,.79,.7,.95);
  L2->Draw("same");

  TCanvas *c3 = plot.Canvas ("c3",400,400,260,510);
  TH1D *f3 = plot.Frame("f3","E_{true} [GeV]","(E_{reco}-E_{true})/E_{true}",0.,15.,-0.3,.3);
  f3->Draw();
  h_Ereso1->Draw("Col same");
  h_Ereso1_mean->Draw("same");
  TLegend *L3 = plot.Legend(" PbGl " ,0.5,.79,.7,.95);
  L3->Draw("same");

  TCanvas *c4 = plot.Canvas ("c4",400,400,510,510);
  TH1D *f4 = plot.Frame("f4","<E_{1}+E_{2}> [GeV]","(m_{reco}-m_{true})/m_{true}",0.,5.,0.,1.3);
  f4->Draw();
  h_piMass->Draw("Col same");
  h_piMass_mean->Draw("same");
  TLegend *L4 = plot.Legend(" PbSc " ,0.5,.29,.7,.35);
  L4->Draw("same");

  TCanvas *c5 = plot.Canvas ("c5",400,400,760,510);
  TH1D *f5 = plot.Frame("f5","<E_{1}+E_{2}> [GeV]","(m_{reco}-m_{true})/m_{true}",0.,5.,0.,1.3);
  f5->Draw();
  h_piMass1->Draw("Col same");
  h_piMass1_mean->Draw("same");
  TLegend *L5 = plot.Legend(" PbGl " ,0.5,.29,.7,.35);
  L5->Draw("same");

  TCanvas *c6 = plot.Canvas ("c6",400,400,410,10);
  TH1D *f6 = plot.Frame("f6","<E_{1}+E_{2}> [GeV]","(m_{reco}-m_{true})/m_{true}",0.,5.,0.,1.3);
  f6->Draw();
  h_piMassPos->Draw("Col same");
  h_piMassPos_mean->Draw("same");
  TLegend *L6 = plot.Legend(" PbGl - position resolution only " ,0.3,.29,.7,.35);
  L6->Draw("same");

  TCanvas *c10 = plot.Canvas ("c10",500,400,1010,10);
  TH1D *frame10 = plot.Frame("frame10","E_{reco} [GeV]","#sigma_{E} [%]",0.,8.,0.,.25);
  frame10->Draw();
  h_Ereso_sigma->Draw("same");
  h_piMass_sigma->Scale(sqrt(2));
  h_piMass_sigma->Draw("same");
  h_piMassPos_sigma->Scale(sqrt(2));
  h_piMassPos_sigma->Draw("sameL");
  EMCalReso->Draw("sameL");

  TLegend *L10 = plot.Legend("PbSc fast simulation ",0.55,0.6,0.88,0.95,1.);
  L10->AddEntry(h_Ereso_sigma,"energy resolution","p");
  L10->AddEntry(h_piMass_sigma,"#sqrt{2} x #sigma ( m_{ #pi^{0} }) ","p");
  L10->AddEntry(h_piMassPos_sigma,"#pi^{0} only position resolution ","p");
  L10->AddEntry(EMCalReso,"8.1%/#sqrt{E_{true}} + 2.1%","l");
  L10->Draw("same");

  TCanvas *c11 = plot.Canvas ("c11",500,400,1010,410);
  TH1D *frame11 = plot.Frame("frame11","E_{reco} [GeV]","#sigma_{E_{true}} [%]",0.,8.,0.,.25);
  frame11->Draw();
  h_Ereso1_sigma->Draw("same");
  h_piMass1_sigma->Scale(sqrt(2));
  h_piMass1_sigma->Draw("same");
  h_piMassPos1_sigma->Scale(sqrt(2));
  h_piMassPos1_sigma->Draw("sameL");
  EMCalReso1->Draw("sameL");

  TLegend *L11 = plot.Legend("PbGl fast simulation ",0.55,0.6,0.88,0.95,1.);
  L11->AddEntry(h_Ereso1_sigma,"energy resolution","p");
  L11->AddEntry(h_piMass1_sigma,"#sqrt{2} x #sigma ( m_{ #pi^{0} })","p");
  L11->AddEntry(h_piMassPos1_sigma,"#pi^{0} only position resolution","p");
  L11->AddEntry(EMCalReso1,"5.9%/#sqrt{E_{true}} + 0.8%","l");
  L11->Draw("same");



}

Double_t f_EcalReso(Double_t *x, Double_t *p) {
//
  Double_t energy = x[0];
  Double_t c1    = p[0];
  Double_t c2    = p[1];
// parameters are square of resolution values!
  Double_t sigma = sqrt(c1*c1/energy + c2*c2); 
//  cout << energy << "  " << sigma << " c1=" << c1 << " c2=" << c2 << endl;  
  return sigma; 
}