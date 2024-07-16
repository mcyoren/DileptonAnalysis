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
Double_t f_DCReso(Double_t *x,Double_t *p);
Double_t f_EoPReso(Double_t *x,Double_t *p);


///////////////////////////////////////////////////////////////////////////////
void TestElectron(){
//
  TRandom3 MyRandy = TRandom3(0);                // Random Generator

  cout << endl;

  PHENIXDetector PHENIX;
  MyPlot plot;

  Particle electron("electron");
  Particle positron("positron");

  TLorentzVector eTrack, eEMCal; 
  TLorentzVector pTrack, pEMCal; 

  Double_t E1,E1nb,p1,pt1,phi1,eta1,ww1;
  Double_t E2,E2nb,p2,pt2,phi2,eta2,ww2;
  Double_t E1_true,E2_true,p1_true,p2_true;
  Int_t q1,q2;
  Double_t P1,P2;

  Double_t BE1=0, BE2=0;
   
  Double_t pt_min = 0.1;
  Double_t pt_max = 15.;
  Double_t pt_low = 0.0;
  Double_t pt_high = 15.;
  Double_t ptcut = 0.2;
  Int_t nevt = 1000000;
  Int_t arm1, arm2; 

  Bool_t bp = false;
  Bool_t bm = false;

  Int_t is1,is2;
  Double_t y1,z1,y2,z2;
  Double_t sinT1,sinT2;
  Double_t phi_EMCal1, phi_EMCal2; 
  Bool_t isPbSc1 = false;
  Bool_t isPbSc2 = false;



// this is not the righ tdistribution, need to make an electron distribution
  TF1 *piHagedorn      = HagedornYield("piHagedorn", pi0Mass, pt_max, pt_min);

  TF1 *EMCalReso = new TF1("EMCalReso",f_EcalReso,0.,10.,2);
  EMCalReso->SetParameters(sigma_E_PbSc_c1,sigma_E_PbSc_c2); 
  TF1 *MomReso = new TF1("MomReso",f_DCReso,0.,10.,2);
  MomReso->SetParameters(sigma_p_c1,sigma_p_c2); 
  TF1 *EoPReso = new TF1("EoPReso",f_EoPReso,0.,10.,4);
  EoPReso->SetParameters(sigma_p_c1,sigma_p_c2,sigma_E_PbSc_c1,sigma_E_PbSc_c2); 

  plot.StyleMe(EMCalReso, 0, kRed, 0, 2, 1.5);
  plot.StyleMe(MomReso, 0, kRed, 0, 2, 1.5);
  plot.StyleMe(EoPReso, 0, kBlue, 0, 1, 2.5);

  TH1D *h_pt           = new TH1D("h_pt","pt",100,pt_low,pt_high);
  plot.StyleMe(h_pt, 20, kBlue, 1., 1, 2); 
  TH1D *h_ptreco       = new TH1D("h_ptreco","ptreco",100,pt_low,pt_high);
  plot.StyleMe(h_ptreco, 20, kTeal, 1., 1, 2); 
  TH1D *h_Ereco        = new TH1D("h_Ereco","Ereco",100,pt_low,pt_high);
  plot.StyleMe(h_Ereco, 20, kOrange+2, 1., 1, 2); 

  TH2D *h_butsyk           = new TH2D("h_butsyk","butsyk",100,-2*pi,2*pi,100,-6.0,6.0);
  TH2D *h_EoP              = new TH2D("h_EoP","EoP",250,pt_low,pt_max,250,0.4,1.6);
  TH2D *h_EoP_Ereso        = new TH2D("h_EoP_Ereso","EoP_E",250,pt_low,pt_max,250,0.4,1.6);
  TH2D *h_EoP_ptreso       = new TH2D("h_EoP_ptreso","EoP_p",250,pt_low,pt_max,250,0.4,1.6);
  TH2D *h_EoP_corr         = new TH2D("h_EoP_corr","EoP_corr",250,pt_low,pt_max,250,0.9,1.1);
  TH2D *h_EoP_brems        = new TH2D("h_EoP_brems","EoP_brems",250,pt_low,pt_max,250,0.0,1.5);
  

  for (int i=0; i<nevt; i++){                                     // event loop
  
    if(i%100000==0) cout << "event loop at " << i << endl;         // event counter printed 
/////////////////////////////////////////////// pi decays ////////////////////////////////////////////////
    bp = false;
    bm = false;

//    cout << endl;    

    electron.ResetP();
    electron.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi with flat distribution
    electron.SetWeight(piHagedorn->Eval(electron.Pt())/float(nevt));
    positron.ResetP();
    positron.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi with flat distribution
    positron.SetWeight(piHagedorn->Eval(positron.Pt())/float(nevt));

    h_pt->Fill(positron.Pt(),positron.Weight());
    h_pt->Fill(electron.Pt(),electron.Weight());
 
    q1 = electron.Charge();
    q2 = positron.Charge();
    p1_true = electron.P();
    E1_true = electron.E();
    p2_true = positron.P();
    E2_true = positron.E();

    PHENIX.CharacterizeTrack(electron, q1, electron.ID());
    arm1 = PHENIX.Arm();
    if (arm1 > 0) {          // monitor PEHNIX acceptance 
      bm = true;
      h_butsyk->Fill(electron.Phi(),q1/electron.Pt(),electron.Weight());  
      is1 = PHENIX.Sector(); 
      phi_EMCal1 = PHENIX.Phi_EMCal(); 
      y1 = PHENIX.SectorY();
      z1 = PHENIX.SectorZ();
      sinT1 = PHENIX.SectorSinT();
//
      eTrack = PHENIX.ReconstructTrack(electron, q1);
      eEMCal = PHENIX.ReconstructShower(electron, electron.ID()); 
      E1nb   = eEMCal.E();
      BE1  = BremsEnergyLoss(.05)*E1_true;
      electron.SetE(E1_true-BE1);
      eEMCal = PHENIX.ReconstructShower(electron, electron.ID()); 
      E1   = eEMCal.E();
      p1   = eTrack.P();
      pt1  = eTrack.Pt();
      phi1 = eTrack.Phi();
      eta1 = eTrack.Eta();
      P1 = p1;
    }

    PHENIX.CharacterizeTrack(positron, q2, positron.ID());
    arm2 = PHENIX.Arm();
    if (arm2 > 0) {          // monitor PEHNIX acceptance 
      bp = true;
      h_butsyk->Fill(positron.Phi(),q2/positron.Pt(),positron.Weight());                   //      for charged particles
//
      is2 = PHENIX.Sector(); 
      phi_EMCal2 = PHENIX.Phi_EMCal(); 
      y2 = PHENIX.SectorY();
      z2 = PHENIX.SectorZ();
      sinT2 = PHENIX.SectorSinT();
//
      pTrack = PHENIX.ReconstructTrack(positron, q2);
      pEMCal = PHENIX.ReconstructShower(positron, positron.ID()); 
      E2nb   = pEMCal.E();
      BE2  = BremsEnergyLoss(.05)*E2_true;
      positron.SetE(E2_true-BE2);
      pEMCal = PHENIX.ReconstructShower(positron, positron.ID()); 
      E2   = pEMCal.E();
      p2   = pTrack.P();
      pt2  = pTrack.Pt();
      phi2 = pTrack.Phi();
      eta2 = pTrack.Eta();
      P2   = p2;
    }
    
    isPbSc1 = false; 
    if (bm){
      if (i<10) { 
        cout << endl;
        cout << i << endl;
        cout << " --- energy of electron " << electron.E() << endl; 
        if (is1>0) {
           cout << " --- total Energy lost from bremsstrahlung " << BE1*electron.E() << endl;
        }
      }
      h_ptreco->Fill(pt1,electron.Weight());
      if (E1>0 && pt1>ptcut) {
        isPbSc1 = (is1 > 0 && is1 <7 );
        if (isPbSc1) {
          h_Ereco->Fill(E1,electron.Weight());
          float Ecorr = PHENIX.ShowerEnergyCorrection(electron.E(),sinT1)
                        /PHENIX.ShowerEnergyCorrection(electron.E(),y1,z1);
          float E = E1*Ecorr;
          P1 = E;
//          P1 = PHENIX.ElectronMomentum(E,p1);
          h_EoP->Fill(P1,E/p1,electron.Weight());
//
          h_EoP_Ereso->Fill(P1,E1nb/p1_true,electron.Weight());
          h_EoP_ptreso->Fill(P1,E1_true/p1,electron.Weight());
          h_EoP_corr->Fill(P1,Ecorr,electron.Weight());
          h_EoP_brems->Fill(P1,E1/p1_true,electron.Weight());
         }
      }
    }
  
    isPbSc2 = false; 
    if (bp) {
      h_ptreco->Fill(pt2,positron.Weight());
      if (E2>0 && pt2>ptcut) {
        isPbSc2 = (is2 > 0 && is2 <7 );
        if (isPbSc2) {
          h_Ereco->Fill(E2,positron.Weight());
          float Ecorr = PHENIX.ShowerEnergyCorrection(positron.E(),sinT2)
                        /PHENIX.ShowerEnergyCorrection(positron.E(),y2,z2);
          float E = E2*Ecorr;
          P2 = E;
//          P2 = PHENIX.ElectronMomentum(E,p2);
          h_EoP->Fill(P2,E/p2,positron.Weight());
//
          h_EoP_Ereso->Fill(P2,E2nb/p2_true,positron.Weight());
          h_EoP_ptreso->Fill(P2,E2_true/p2,positron.Weight());
          h_EoP_corr->Fill(P2,Ecorr,positron.Weight());
          h_EoP_brems->Fill(P2,E2/p2_true,positron.Weight());
         }
      }
    }
  }

  if (nevt>1000){
  h_EoP->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_EoP_mean = (TH1D*)gDirectory->Get("h_EoP_1");
  h_EoP_mean->SetName("h_EoP_mean");
  plot.StyleMe(h_EoP_mean,20,kTeal,.8);
  TH1D * h_EoP_sigma = (TH1D*)gDirectory->Get("h_EoP_2");
  h_EoP_sigma->SetName("h_EoP_sigma");
  plot.StyleMe(h_EoP_sigma,20,kTeal,.8);

  h_EoP_Ereso->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_EoP_Emean = (TH1D*)gDirectory->Get("h_EoP_Ereso_1");
  h_EoP_Emean->SetName("h_EoP_Emean");
  plot.StyleMe(h_EoP_Emean,20,kOrange+2,.8);
  TH1D * h_EoP_Esigma = (TH1D*)gDirectory->Get("h_EoP_Ereso_2");
  h_EoP_Esigma->SetName("h_EoP_Esigma");
  plot.StyleMe(h_EoP_Esigma,20,kOrange+2,.8);

  h_EoP_ptreso->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_EoP_ptmean = (TH1D*)gDirectory->Get("h_EoP_ptreso_1");
  h_EoP_ptmean->SetName("h_EoP_ptmean");
  plot.StyleMe(h_EoP_ptmean,20,kGreen+2,.8);
  TH1D * h_EoP_ptsigma = (TH1D*)gDirectory->Get("h_EoP_ptreso_2");
  h_EoP_ptsigma->SetName("h_EoP_ptsigma");
  plot.StyleMe(h_EoP_ptsigma,20,kGreen+2,.8);

  h_EoP_corr->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_EoP_cmean = (TH1D*)gDirectory->Get("h_EoP_corr_1");
  h_EoP_cmean->SetName("h_EoP_cmean");
  plot.StyleMe(h_EoP_cmean,20,kGray+2,.8);
  TH1D * h_EoP_csigma = (TH1D*)gDirectory->Get("h_EoP_corr_2");
  h_EoP_csigma->SetName("h_EoP_csigma");
  plot.StyleMe(h_EoP_csigma,20,kGray+2,.8);

  h_EoP_brems->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_EoP_bmean = (TH1D*)gDirectory->Get("h_EoP_brems_1");
  h_EoP_bmean->SetName("h_EoP_bmean");
  plot.StyleMe(h_EoP_bmean,20,kMagenta-1,.8);
  TH1D * h_EoP_bsigma = (TH1D*)gDirectory->Get("h_EoP_brems_2");
  h_EoP_bsigma->SetName("h_EoP_bsigma");
  plot.StyleMe(h_EoP_bsigma,20,kMagenta-1,.8);
  
 ////////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = plot.Canvas ("c1",400,500,10,10,1);
  plot.SetyTitleOffset(1.3);
  TH1D *frame1 = plot.Frame("frame1","pt","counts",pt_low,pt_high*.999,1e-12,1.);
  frame1->Draw();
  h_pt->Draw("sameL Chist");
  h_ptreco->Draw("sameL Chist");
  h_Ereco->Draw("sameL Chist");

  plot.SetyTitleOffset(0.6);
  plot.SetLeftMargin(0.1);
  TCanvas *c2 = plot.Canvas ("c2",600,300,410,10);
  TH1D *frame2 = plot.Frame("frame2","phi","q/pt",-2.*pi,2*pi,-6.,6.);
  frame2->Draw();
  h_butsyk->Draw("Col same");
  plot.Reset();
 
  TCanvas *c3 = plot.Canvas ("c3",400,400,10,510);
  TH1D *frame3 = plot.Frame("frame3","p_{reco} [Gev/c]","E/p reconstrcuted",0.,10.,0.4,1.6);
  frame3->Draw();
  h_EoP->Draw("Col same");
  h_EoP_mean->Draw("same");

  TCanvas *c4 = plot.Canvas ("c4",400,400,210,510);
  TH1D *frame4 = plot.Frame("frame4","p_{reco} [Gev/c]","E/p reconstrcuted",0.,10.,0.4,1.6);
  frame4->Draw();
  h_EoP_Ereso->Draw("Col same");
  h_EoP_Emean->Draw("same");
  
  TCanvas *c5 = plot.Canvas ("c5",400,400,410,510);
  TH1D *frame5 = plot.Frame("frame5","p_{reco} [Gev/c]","E/p reconstrcuted",0.,10.,0.4,1.6);
  frame5->Draw();
  h_EoP_ptreso->Draw("Col same");
  h_EoP_ptmean->Draw("same");

  TCanvas *c6 = plot.Canvas ("c6",400,400,610,510);
  TH1D *frame6 = plot.Frame("frame6","p_{reco} [Gev/c]","E/p reconstrcuted",0.,10.,0.9,1.1);
  frame6->Draw();
  h_EoP_corr->Draw("Col same");
  h_EoP_cmean->Draw("same");

  TCanvas *c7 = plot.Canvas ("c7",400,400,810,510);
  TH1D *frame7 = plot.Frame("frame7","p_{reco} [Gev/c]","E/p reconstrcuted",0.,10.,0.0,1.5);
  frame7->Draw();
  h_EoP_brems->Draw("Col same");
  h_EoP_bmean->Draw("same");

  TCanvas *c10 = plot.Canvas ("c10",500,300,1010,10);
  TH1D *frame10 = plot.Frame("frame10","E_{reco}","#sigma_{E/p} (%)",0.,10.,0.,.2);
  frame10->Draw();
  h_EoP_sigma->Draw("same");
  h_EoP_Esigma->Draw("same");
  h_EoP_ptsigma->Draw("same");
  h_EoP_csigma->Draw("same");
  h_EoP_bsigma->Draw("same");
  EMCalReso->Draw("sameL");
  MomReso->Draw("sameL");
  EoPReso->Draw("sameL");

  TLegend *L10 = plot.Legend("E/p Fast Simulation Results",0.3,.6,.46,.92);
  L10->AddEntry(h_EoP_mean," all effects","p");
  L10->AddEntry(h_EoP_Emean," E_{reco}/p_{true} ","p");
  L10->AddEntry(h_EoP_bmean," E_{reco}/p_{true} + bremsstrahlung","p");
  L10->AddEntry(h_EoP_ptmean," E_{true}/p_{reco}","p");
  L10->AddEntry(h_EoP_cmean," impact angle contribution","p");
  L10->AddEntry(EoPReso, " energy and momentum resolution ");
  L10->Draw("same");


  TCanvas *c11 = plot.Canvas ("c11",500,300,1010,310);
  TH1D *frame11 = plot.Frame("frame11","E_{reco}","<E/p> ",0.,10.,0.92,1.02);
  frame11->Draw();
  h_EoP_mean->Draw("same");
  h_EoP_Emean->Draw("same");
  h_EoP_ptmean->Draw("same");
  h_EoP_cmean->Draw("same");
  h_EoP_bmean->Draw("same");
//  EMCalReso->Draw("sameL");
//  MomReso->Draw("sameL");
//  EoPReso->Draw("sameL");

  TLegend *L11 = plot.Legend("E/p Fast Simulation Results",0.3,.2,.5,.52);
  L11->AddEntry(h_EoP_mean," all effects","p");
  L11->AddEntry(h_EoP_Emean," E_{reco}/p_{true} ","p");
  L11->AddEntry(h_EoP_bmean," E_{reco}/p_{true} + bremsstrahlung","p");
  L11->AddEntry(h_EoP_ptmean," E_{true}/p_{reco}","p");
  L11->AddEntry(h_EoP_cmean," impact angle contribution","p");
  L11->Draw("same");


  }
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

Double_t f_DCReso(Double_t *x, Double_t *p) {
//
  Double_t pt    = x[0];
  Double_t c1    = p[0];
  Double_t c2    = p[1];
// parameters are square of resolution values!
  Double_t sigma = sqrt(c1*c1 + c2*c2*pt*pt); 
  return sigma; 
}

Double_t f_EoPReso(Double_t *x, Double_t *p) {
//
  Double_t E       = x[0];
  Double_t DCc1    = p[0];
  Double_t DCc2    = p[1];
  Double_t c1      = p[2];
  Double_t c2      = p[3];
// parameters are square of resolution values!
  Double_t sigma = sqrt(DCc1*DCc1 + DCc2*DCc2*E*E + c1*c1/E + c2*c2); 
  return sigma; 
}

