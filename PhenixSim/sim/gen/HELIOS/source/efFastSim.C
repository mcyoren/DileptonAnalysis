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
void efFastSim(){
//
  PHENIXDetector PHENIX;                       // define PHENIX detector 
  Particle pi0("pi0");                   // define particle - used for random generation and decay 

  TLorentzVector gamma1, gamma2;         // 4 vectors for pi0->gg
  TLorentzVector electron,positron;      // 4 vector for  g->e+e-

  TLorentzVector electron_reco, positron_reco;  // dito but with recostructed 4 vector
  TLorentzVector gamma1_reco, gamma2_reco;      //  

  Double_t E_up,E_down;
  Double_t k0,k1,k2;
   
  Double_t ww, ww2;                      // weight of parent particle

  Double_t pt_min = 0.4;                 // range of pi0 momentum
  Double_t pt_max = 12.;                 // 
  Double_t pt_low = 0.8;                 // range of histograms
  Double_t pt_high = 10.8;
  Double_t ptcut_t = 0.18;
  Double_t ptcut   = 0.2;
  Int_t nevt = 100000000;

// flags if particle is in PHENIX Acceptance
  Bool_t iA_pi0 = false; 
  Bool_t iA_g1  = false;
  Bool_t iA_g2  = false;
  Bool_t iA_ep  = false;
  Bool_t iA_e   = false;
  Bool_t iA_p   = false;


  Int_t bins = 200;


  TH1::SetDefaultSumw2(true);   
  TF1 *piHagedorn2      = HagedornYield("piHagedorn2", pi0Mass, pt_max, pt_min);
// Norbert
  TF1 *piHagedornNN   = HagedornYield("piHagedornNN", pi0Mass, pt_max, pt_min, 68, 0.31, 0.13, 0.66, -8.14);
  // Roli
  TF1 *piHagedornRE     = HagedornYield("piHagedornRE", pi0Mass, pt_max, pt_min, 58, 0.661, 0.015, 0.745, -9.167);
// ppg088
  TF1 *piHagedorn     = HagedornYield("piHagedorn", pi0Mass, pt_max, pt_min, 377., 0.356, 0.068, 0.7, -8.25);
  TH1D *h_ptpi0_t      = new TH1D("h_ptpi0t","pt",bins,pt_low,pt_high);
  TH1D *h_ptpi0        = new TH1D("h_ptpi0","pt",bins,pt_low,pt_high);
  TH1D *h_ptpi02        = new TH1D("h_ptpi02","pt",bins,pt_low,pt_high);
  TH1D *h_ptg          = new TH1D("h_ptg","ptg",bins,pt_low,pt_high);
  TH1D *h_pte          = new TH1D("h_ptg","pte",bins,pt_low,pt_high);

// step by step <ef>
  TH1D *h_Nincl_t      = new TH1D("h_Nincl_t","Nincl_t",bins,pt_low,pt_high);
  TH1D *h_Nincl        = new TH1D("h_Nincl","Nincl",bins,pt_low,pt_high);
  TH1D *h_Nincl2       = new TH1D("h_Nincl2","Nincl2",bins,pt_low,pt_high);
  TH1D *h_Ntag_ideal   = new TH1D("h_Ntag_ideal","Ntag_ideal",bins,pt_low,pt_high);
  TH1D *h_Ntag_fiducial= new TH1D("h_Ntag_fiducial","Ntag_fiducial",bins,pt_low,pt_high);
  TH1D *h_Ntag_deadarea= new TH1D("h_Ntag_deadarea","Ntag_deadarea",bins,pt_low,pt_high);
  TH1D *h_Ntag_VTXconv = new TH1D("h_Ntag_VTXconv","Ntag_VTXconv",bins,pt_low,pt_high);
  TH1D *h_Ntag_300     = new TH1D("h_Ntag_300","Ntag_300",bins,pt_low,pt_high);
  TH1D *h_Ntag_400     = new TH1D("h_Ntag_400","Ntag_400",bins,pt_low,pt_high);
  TH1D *h_Ntag_500     = new TH1D("h_Ntag_500","Ntag_500",bins,pt_low,pt_high);
  TH1D *h_Ntag_600     = new TH1D("h_Ntag_600","Ntag_600",bins,pt_low,pt_high);
 
 // energy scale systematic +/- 2%
  TH1D *h_Ntag_300_d   = new TH1D("h_Ntag_300_d","Ntag_300d",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_d   = new TH1D("h_Ntag_400_d","Ntag_400d",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_d   = new TH1D("h_Ntag_500_d","Ntag_500d",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_d   = new TH1D("h_Ntag_600_d","Ntag_600d",bins,pt_low,pt_high);
  TH1D *h_Ntag_300_u   = new TH1D("h_Ntag_300_u","Ntag_300u",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_u   = new TH1D("h_Ntag_400_u","Ntag_400u",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_u   = new TH1D("h_Ntag_500_u","Ntag_500u",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_u   = new TH1D("h_Ntag_600_u","Ntag_600u",bins,pt_low,pt_high);
 
// energy resolution change low energy part by +/-10%
  TH1D *h_Ntag_300_rd   = new TH1D("h_Ntag_300_rd","Ntag_300rd",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_rd   = new TH1D("h_Ntag_400_rd","Ntag_400rd",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_rd   = new TH1D("h_Ntag_500_rd","Ntag_500rd",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_rd   = new TH1D("h_Ntag_600_rd","Ntag_600rd",bins,pt_low,pt_high);
  TH1D *h_Ntag_300_ru   = new TH1D("h_Ntag_300_ru","Ntag_300ru",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_ru   = new TH1D("h_Ntag_400_ru","Ntag_400ru",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_ru   = new TH1D("h_Ntag_500_ru","Ntag_500ru",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_ru   = new TH1D("h_Ntag_600_ru","Ntag_600ru",bins,pt_low,pt_high);
 
// energy non linearity
  TH1D *h_Ntag_300_nld   = new TH1D("h_Ntag_300_nld","Ntag_300nld",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_nld   = new TH1D("h_Ntag_400_nld","Ntag_400nld",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_nld   = new TH1D("h_Ntag_500_nld","Ntag_500nld",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_nld   = new TH1D("h_Ntag_600_nld","Ntag_600nld",bins,pt_low,pt_high);
  TH1D *h_Ntag_300_nlu   = new TH1D("h_Ntag_300_nlu","Ntag_300nlu",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_nlu   = new TH1D("h_Ntag_400_nlu","Ntag_400nlu",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_nlu   = new TH1D("h_Ntag_500_nlu","Ntag_500nlu",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_nlu   = new TH1D("h_Ntag_600_nlu","Ntag_600nlu",bins,pt_low,pt_high);

// different pt distribution
  TH1D *h_Ntag_300_pt   = new TH1D("h_Ntag_300_pt","Ntag_300pt",bins,pt_low,pt_high);
  TH1D *h_Ntag_400_pt   = new TH1D("h_Ntag_400_pt","Ntag_400pt",bins,pt_low,pt_high);
  TH1D *h_Ntag_500_pt   = new TH1D("h_Ntag_500_pt","Ntag_500pt",bins,pt_low,pt_high);
  TH1D *h_Ntag_600_pt   = new TH1D("h_Ntag_600_pt","Ntag_600pt",bins,pt_low,pt_high);



  MyPlot plot;

  plot.StyleMe(h_ptpi0,   20, kBlue, 1., 1, 2); 
  plot.StyleMe(h_ptpi02,   20, kRed, 1., 2, 2); 
  plot.StyleMe(h_ptg  , 20, kTeal, 1., 1, 2); 
  plot.StyleMe(h_pte  , 20, kOrange+2, 1., 1, 2); 
  plot.StyleMe(h_Nincl_t , 20, kGray+2, 1., 2, 2); 
  plot.StyleMe(h_Nincl , 20, kBlack, 1., 1., 2); 
  plot.StyleMe(h_Ntag_ideal, 20, kBlack, 1., 1., 2); 
  plot.StyleMe(h_Ntag_fiducial, 20, kBlue, 1., 1., 2); 
  plot.StyleMe(h_Ntag_deadarea, 20, kGreen+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_VTXconv, 20, kRed, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300, 20, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400, 20, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500, 20, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600, 20, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_pt, 21, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_pt, 21, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_pt, 21, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_pt, 21, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_d, 24, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_d, 24, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_d, 24, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_d, 24, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_u, 24, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_u, 24, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_u, 24, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_u, 24, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_rd, 25, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_rd, 25, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_rd, 25, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_rd, 25, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_ru, 25, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_ru, 25, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_ru, 25, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_ru, 25, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_nld, 33, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_nld, 33, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_nld, 33, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_nld, 33, kMagenta+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_300_nlu, 33, kTeal-1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_400_nlu, 33, kOrange+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_500_nlu, 33, kAzure+1, 1., 1., 2); 
  plot.StyleMe(h_Ntag_600_nlu, 33, kMagenta+1, 1., 1., 2); 

// counters to monitor overall <ef>
  float cc1=0, cc2=0, cc3=0, cc4=0,  cc5=0,  cc6=0, 
        cc7=0, cc8=0, cc9=0, cc10=0, cc11=0, cc12=0,
        cc13=0, cc14=0, cc15=0;
                             
///////////////////////////////////////////////////////////////////////////////
// start of event loop
cout << "hello world" << endl;
  for (int i=0; i<nevt; i++)
  {
    if(i%100000==0) cout << "event loop at " << i << endl;      
// generate pi0
    pi0.ResetP();                                              
    pi0.GenerateP(pt_min,pt_max);                              
    pi0.SetWeight(piHagedorn->Eval(pi0.Pt())/float(nevt));  
    ww = pi0.Weight();
    ww2 = piHagedorn2->Eval(pi0.Pt())/float(nevt)/2.;

// cout << ww << "  "  << ww2 << "  " << ww2/ww << endl; 


    cc1=cc1+ww;
    h_ptpi0_t->Fill(pi0.Pt(),ww);
    iA_pi0 = true;
    while (iA_pi0) 
    {                          
// check if pi0 in central arm acceptance
      iA_pi0 = (PHENIX.InAcceptance(pi0,pi0.Charge())>0);        
      if (!iA_pi0) break;
      cc2=cc2+ww;
      h_ptpi0->Fill(pi0.Pt(),ww);
      h_ptpi02->Fill(pi0.Pt(),ww2);

// decay to gamma gamma
      pi0.DecaySingleBranch("pi0->gg");
      gamma1 = pi0.GetDecayDaughter(0);
      gamma2 = pi0.GetDecayDaughter(1);
      h_ptg->Fill(gamma1.Pt(),ww);
      h_ptg->Fill(gamma2.Pt(),ww); 

// ceck if gamma1 is in PHENIX acceptance
      iA_g1 = (PHENIX.InAcceptance(gamma1,0)>0);
      if (!iA_g1) break;
      cc3=cc3+ww;
// convert gamma1
      PhotonConvert(gamma1, electron, positron);
      h_pte->Fill(electron.Pt(),ww);
      h_pte->Fill(positron.Pt(),ww);

// check that electron and positron are in acceptance of same arm
      if (electron.Pt()<ptcut_t) break;
      if (positron.Pt()<ptcut_t) break;  
      cc4=cc4+ww;
      iA_ep =  (PHENIX.InAcceptance(electron,-1) > 0)
              && (PHENIX.InAcceptance(positron,1) > 0)
              && (PHENIX.InAcceptance(electron,-1) == PHENIX.InAcceptance(positron,+1))
              ;
      if (!iA_ep) break;
      cc5=cc5+ww;
      h_Nincl_t->Fill(gamma1.Pt(),ww);  

// reconstruct electron and positron
      PHENIX.CharacterizeTrack(electron,-1,electronID);
      electron_reco = PHENIX.ReconstructTrack(electron, -1);
      PHENIX.CharacterizeTrack(positron,1,positronID);
      positron_reco = PHENIX.ReconstructTrack(positron,  1);
    
//electron.Print(); 
//electron_reco.Print(); 

      if (electron_reco.Pt()<ptcut) break;
      if (positron_reco.Pt()<ptcut) break;  
      iA_ep =  (PHENIX.InAcceptance(electron_reco,-1) > 0)
               && (PHENIX.InAcceptance(positron_reco,1) > 0)
               && (PHENIX.InAcceptance(electron_reco, -1) == PHENIX.InAcceptance(positron_reco,+1))
               ;
      if (!iA_ep) break;
      cc6=cc6+ww;
      gamma1_reco = electron_reco+positron_reco;
      h_Nincl->Fill(gamma1_reco.Pt(),ww);
      h_Nincl2->Fill(gamma1_reco.Pt(),ww2);
//
// from here on out calculate varios conditional acceptances 
//
// ceck if gamma2 is in PHENIX acceptance
      iA_g2 = (PHENIX.InAcceptance(gamma2,0)>0);
      if (!iA_g2) break;
      cc7=cc7+ww;

// // reconstruct gamma2
      gamma2_reco = PHENIX.ReconstructShower(gamma2,photonID);
      iA_g2 = (PHENIX.InAcceptance(gamma2_reco,0)>0);
      if (!iA_g2) break;
      cc8=cc8+ww;
      h_Ntag_ideal->Fill(gamma1_reco.Pt(),ww);

// is in fiducial volume of EMCal
      int sector = PHENIX.EMCalSector(gamma2.Phi());
      if ( sector == 0) break;
      cc9=cc9+ww;
      h_Ntag_fiducial->Fill(gamma1_reco.Pt(),ww);

// is not in dead area
      if(!PHENIX.EMCalLive(sector)) break;
      cc10=cc10+ww;
      h_Ntag_deadarea->Fill(gamma1_reco.Pt(),ww);

// does not convert in VTX
      if (PHENIX.VTXConversion() != 0) break;
      cc11=cc11+ww;
      h_Ntag_VTXconv->Fill(gamma1_reco.Pt(),ww);

// central value for energy
      Bool_t e1 = true; 
      while (e1) {
    // energy cut 300 MeV
          if (gamma2_reco.E() < 0.3) break;
          cc12=cc12+ww;
          h_Ntag_300->Fill(gamma1_reco.Pt(),ww);
          h_Ntag_300_pt->Fill(gamma1_reco.Pt(),ww2);
           
    // energy cut 400 MeV
          if (gamma2_reco.E() < 0.4) break;
          cc13=cc13+ww;
          h_Ntag_400->Fill(gamma1_reco.Pt(),ww);
          h_Ntag_400_pt->Fill(gamma1_reco.Pt(),ww2);

    // energy cut 500 MeV
          if (gamma2_reco.E() < 0.5) break;
          cc14=cc14+ww;
          h_Ntag_500->Fill(gamma1_reco.Pt(),ww);
          h_Ntag_500_pt->Fill(gamma1_reco.Pt(),ww2);

    // energy cut 600 MeV
          if (gamma2_reco.E() < 0.6) break;
          cc15=cc15+ww;
          h_Ntag_600->Fill(gamma1_reco.Pt(),ww);
          h_Ntag_600_pt->Fill(gamma1_reco.Pt(),ww2);
          e1 = false;
      }

// lower value for energy
      float escale = 0.99;
      Bool_t e2 = true; 
      while (e2) {
    // energy cut 300 MeV
          if (gamma2_reco.E()*escale < 0.3) break;
          h_Ntag_300_d->Fill(gamma1_reco.Pt(),ww);
           
    // energy cut 400 MeV
          if (gamma2_reco.E()*escale < 0.4) break;
          h_Ntag_400_d->Fill(gamma1_reco.Pt(),ww);

    // energy cut 500 MeV
          if (gamma2_reco.E()*escale < 0.5) break;
          h_Ntag_500_d->Fill(gamma1_reco.Pt(),ww);

    // energy cut 600 MeV
          if (gamma2_reco.E()*escale < 0.6) break;
          h_Ntag_600_d->Fill(gamma1_reco.Pt(),ww);
          e2 = false;
      }

// upper value for energy
      escale = 1/escale;
      Bool_t e3 = true; 
      while (e3) {
    // energy cut 300 MeV
          if (gamma2_reco.E()*escale < 0.3) break;
          h_Ntag_300_u->Fill(gamma1_reco.Pt(),ww);
           
    // energy cut 400 MeV
          if (gamma2_reco.E()*escale < 0.4) break;
          h_Ntag_400_u->Fill(gamma1_reco.Pt(),ww);

    // energy cut 500 MeV
          if (gamma2_reco.E()*escale < 0.5) break;
          h_Ntag_500_u->Fill(gamma1_reco.Pt(),ww);

    // energy cut 600 MeV
          if (gamma2_reco.E()*escale < 0.6) break;
          h_Ntag_600_u->Fill(gamma1_reco.Pt(),ww);
          e3 = false;
      }

// smear energy up or down by 10% 
      if (sector < 7){
        k1 = sigma_E_PbSc_c1;
        k2 = sigma_E_PbSc_c2;
      } else {
        k1 = sigma_E_PbGl_c1;
        k2 = sigma_E_PbGl_c2;      
      }
      E_up   = PHENIX.SmearEnergy(gamma2.E(),2,k1*1.1,k2*1.1); 
      E_down = PHENIX.SmearEnergy(gamma2.E(),2,k1/1.1,k2/1.1);

//      cout << gamma1.E() << " " << E_up << " " << E_down << endl;

      Bool_t e4 = true; 
      while (e4) {
    // energy cut 300 MeV
          if (E_down < 0.3) break;
          h_Ntag_300_rd->Fill(gamma1_reco.Pt(),ww);
           
    // energy cut 400 MeV
          if (E_down < 0.4) break;
          h_Ntag_400_rd->Fill(gamma1_reco.Pt(),ww);

    // energy cut 500 MeV
          if (E_down < 0.5) break;
          h_Ntag_500_rd->Fill(gamma1_reco.Pt(),ww);

    // energy cut 600 MeV
          if (E_down < 0.6) break;
          h_Ntag_600_rd->Fill(gamma1_reco.Pt(),ww);
          e4 = false;
      }

      Bool_t e5 = true; 
      while (e5) {
    // energy cut 300 MeV
          if (E_up < 0.3) break;
          h_Ntag_300_ru->Fill(gamma1_reco.Pt(),ww);
           
    // energy cut 400 MeV
          if (E_up < 0.4) break;
          h_Ntag_400_ru->Fill(gamma1_reco.Pt(),ww);

    // energy cut 500 MeV
          if (E_up < 0.5) break;
          h_Ntag_500_ru->Fill(gamma1_reco.Pt(),ww);

    // energy cut 600 MeV
          if (E_up < 0.6) break;
          h_Ntag_600_ru->Fill(gamma1_reco.Pt(),ww);
          e5 = false;
      }

// non linearity from Wenqings thesis 
      if (sector < 7){
        k0 = 1.003;
        k1 = 0.049;
        k2 = 1.766;
      } else {
        k0 = 1.022;
        k1 = 0.083;
        k2 = 0.649;     
      }
      E_up   = PHENIX.NonLinearEnergy(gamma2.E(),1/k0,-k1,k2);
      E_down = PHENIX.NonLinearEnergy(gamma2.E(),k0,k1,k2);

      Bool_t e6 = true; 
      while (e6) {
    // energy cut 300 MeV
          if (E_down < 0.3) break;
          h_Ntag_300_nld->Fill(gamma1_reco.Pt(),ww);
           
    // energy cut 400 MeV
          if (E_down < 0.4) break;
          h_Ntag_400_nld->Fill(gamma1_reco.Pt(),ww);

    // energy cut 500 MeV
          if (E_down < 0.5) break;
          h_Ntag_500_nld->Fill(gamma1_reco.Pt(),ww);

    // energy cut 600 MeV
          if (E_down < 0.6) break;
          h_Ntag_600_nld->Fill(gamma1_reco.Pt(),ww);
          e6 = false;
      }

      Bool_t e7 = true; 
      while (e7) {
    // energy cut 300 MeV
          if (E_up < 0.3) break;
          h_Ntag_300_nlu->Fill(gamma1_reco.Pt(),ww);
           
    // energy cut 400 MeV
          if (E_up < 0.4) break;
          h_Ntag_400_nlu->Fill(gamma1_reco.Pt(),ww);

    // energy cut 500 MeV
          if (E_up < 0.5) break;
          h_Ntag_500_nlu->Fill(gamma1_reco.Pt(),ww);

    // energy cut 600 MeV
          if (E_up < 0.6) break;
          h_Ntag_600_nlu->Fill(gamma1_reco.Pt(),ww);
          e7 = false;
      }




// thats it      
      iA_pi0 = false;
    }
  } 
//  end of event loop
  cout << endl;
  cout << "------------ 1st photon -----------------------------------------"   << endl;
  cout << " PHENIX central arm acceptance           " << cc2/cc1   << endl;
  cout << " 1st photon in acceptance                " << cc3/cc2   << endl;
  cout << " conversion e+e- pass ptcut              " << cc4/cc3   << endl;
  cout << " conversion in acceptance                " << cc5/cc4   << endl;
  cout << " dito using reconstructed momentum       " << cc6/cc4   << endl;
  cout << "------------ 2nd photon -----------------------------------------"   << endl;
  cout << " 2nd photon in ideal acceptance          " << cc7/cc6   << endl;
  cout << " dito using reconstructed momentum       " << cc8/cc6   << endl;
  cout << " in fiducial area                        " << cc9/cc8   << "   " << cc9/cc6 << endl;
  cout << " in live area                            " << cc10/cc9  << "   " << cc10/cc6<< endl;
  cout << " VTX conversion loss                     " << cc11/cc10 << "   " << cc11/cc6<< endl;
  cout << " Ecut > 300 MeV                          " << cc12/cc10 << "   " << cc12/cc6<< endl;
  cout << " Ecut > 400 MeV                          " << cc13/cc12 << "   " << cc13/cc6<< endl;
  cout << " Ecut > 500 MeV                          " << cc14/cc13 << "   " << cc14/cc6<< endl;
  cout << " Ecut > 600 MeV                          " << cc15/cc14 << "   " << cc15/cc6<< endl;
  cout << endl;

/////////////////////////////////////////////////////////////////////////////
//  calculate <ef>

  TH1D *h_ef_nincl  = (TH1D*) h_Nincl_t->Clone("h_ef_nincl");
  plot.RatioBinomial(h_ef_nincl,h_Nincl,h_ef_nincl);
 
  TH1D *h_ef_ideal  = (TH1D*) h_Ntag_ideal->Clone("h_ef_ideal");
  plot.RatioBinomial(h_ef_ideal, h_Nincl, h_ef_ideal);
  TH1D *h_ef_fiducial  = (TH1D*) h_Ntag_fiducial->Clone("h_ef_fiducial");
  plot.RatioBinomial(h_ef_fiducial, h_Nincl, h_ef_fiducial);
  TH1D *h_ef_deadarea  = (TH1D*) h_Ntag_deadarea->Clone("h_ef_deadarea");
  plot.RatioBinomial(h_ef_deadarea, h_Nincl, h_ef_deadarea);
  TH1D *h_ef_VTXconv  = (TH1D*) h_Ntag_VTXconv->Clone("h_ef_VTXconv");
  plot.RatioBinomial(h_ef_VTXconv, h_Nincl, h_ef_VTXconv);
  TH1D *h_ef_300  = (TH1D*) h_Ntag_300->Clone("h_ef_300");
  plot.RatioBinomial(h_Ntag_300, h_Nincl, h_ef_300);
  TH1D *h_ef_400  = (TH1D*) h_Ntag_400->Clone("h_ef_400");
  plot.RatioBinomial(h_ef_400, h_Nincl, h_ef_400);
  TH1D *h_ef_500  = (TH1D*) h_Ntag_500->Clone("h_ef_500");
  plot.RatioBinomial(h_ef_500, h_Nincl, h_ef_500);
  TH1D *h_ef_600  = (TH1D*) h_Ntag_600->Clone("h_ef_600");
  plot.RatioBinomial(h_ef_600, h_Nincl, h_ef_600);


cout<< "------------------------------------------------------" << endl;
// different pt distribution
  TH1D *h_ef_300_pt = (TH1D*) h_Ntag_300_pt->Clone("h_ef_300_pt");
  plot.RatioBinomial(h_Ntag_300_pt, h_Nincl2, h_ef_300_pt);
  TH1D *h_ef_400_pt = (TH1D*) h_Ntag_400_pt->Clone("h_ef_400_pt");
  plot.RatioBinomial(h_Ntag_400_pt, h_Nincl2, h_ef_400_pt);
  TH1D *h_ef_500_pt = (TH1D*) h_Ntag_500_pt->Clone("h_ef_500_pt");
  plot.RatioBinomial(h_Ntag_500_pt, h_Nincl2, h_ef_500_pt);
  TH1D *h_ef_600_pt  = (TH1D*) h_Ntag_600_pt->Clone("h_ef_600_pt");
  plot.RatioBinomial(h_Ntag_600_pt, h_Nincl2, h_ef_600_pt);

  TH1D *R_300_pt  = (TH1D*) h_Ntag_300_pt->Clone("R_300_pt");
  plot.RatioBinomial(h_ef_300_pt, h_ef_300, R_300_pt);
  TH1D *R_400_pt  = (TH1D*) h_Ntag_400_pt->Clone("R_400_pt");
  plot.RatioBinomial(h_ef_400_pt, h_ef_400, R_400_pt);
  TH1D *R_500_pt  = (TH1D*) h_Ntag_500_pt->Clone("R_500_pt");
  plot.RatioBinomial(h_ef_500_pt, h_ef_500, R_500_pt);
  TH1D *R_600_pt  = (TH1D*) h_Ntag_600_pt->Clone("R_600_pt");
  plot.RatioBinomial(h_ef_600_pt, h_ef_600, R_600_pt);
cout<< "------------------------------------------------------" << endl;



// systematics of energy scale
  TH1D *R_300_d  = (TH1D*) h_Ntag_300_d->Clone("R_300_d");
  plot.RatioBinomial(h_Ntag_300_d, h_Ntag_300, R_300_d);
  TH1D *R_400_d  = (TH1D*) h_Ntag_400_d->Clone("R_400_d");
  plot.RatioBinomial(h_Ntag_400_d, h_Ntag_400, R_400_d);
  TH1D *R_500_d  = (TH1D*) h_Ntag_500_d->Clone("R_500_d");
  plot.RatioBinomial(h_Ntag_500_d, h_Ntag_500, R_500_d);
  TH1D *R_600_d  = (TH1D*) h_Ntag_600_d->Clone("R_600_d");
  plot.RatioBinomial(h_Ntag_600_d, h_Ntag_600, R_600_d);
  TH1D *R_300_u  = (TH1D*) h_Ntag_300_u->Clone("R_300_u");
  plot.RatioBinomial(h_Ntag_300_u, h_Ntag_300, R_300_u);
  TH1D *R_400_u  = (TH1D*) h_Ntag_400_u->Clone("R_400_u");
  plot.RatioBinomial(h_Ntag_400_u, h_Ntag_400, R_400_u);
  TH1D *R_500_u  = (TH1D*) h_Ntag_500_u->Clone("R_500_u");
  plot.RatioBinomial(h_Ntag_500_u, h_Ntag_500, R_500_u);
  TH1D *R_600_u  = (TH1D*) h_Ntag_600_u->Clone("R_600_u");
  plot.RatioBinomial(h_Ntag_600_u, h_Ntag_600, R_600_u);

// systematics of energy resolution 
  TH1D *R_300_rd  = (TH1D*) h_Ntag_300_rd->Clone("R_300_rd");
  plot.RatioBinomial(h_Ntag_300_rd, h_Ntag_300, R_300_rd);
  TH1D *R_400_rd  = (TH1D*) h_Ntag_400_rd->Clone("R_400_rd");
  plot.RatioBinomial(h_Ntag_400_rd, h_Ntag_400, R_400_rd);
  TH1D *R_500_rd  = (TH1D*) h_Ntag_500_rd->Clone("R_500_rd");
  plot.RatioBinomial(h_Ntag_500_rd, h_Ntag_500, R_500_rd);
  TH1D *R_600_rd  = (TH1D*) h_Ntag_600_rd->Clone("R_600_rd");
  plot.RatioBinomial(h_Ntag_600_rd, h_Ntag_600, R_600_rd);
  TH1D *R_300_ru  = (TH1D*) h_Ntag_300_ru->Clone("R_300_ru");
  plot.RatioBinomial(h_Ntag_300_ru, h_Ntag_300, R_300_ru);
  TH1D *R_400_ru  = (TH1D*) h_Ntag_400_ru->Clone("R_400_ru");
  plot.RatioBinomial(h_Ntag_400_ru, h_Ntag_400, R_400_ru);
  TH1D *R_500_ru  = (TH1D*) h_Ntag_500_ru->Clone("R_500_ru");
  plot.RatioBinomial(h_Ntag_500_ru, h_Ntag_500, R_500_ru);
  TH1D *R_600_ru  = (TH1D*) h_Ntag_600_ru->Clone("R_600_ru");
  plot.RatioBinomial(h_Ntag_600_ru, h_Ntag_600, R_600_ru);

// systematics of possible nonlinearity 
  TH1D *R_300_nld  = (TH1D*) h_Ntag_300_nld->Clone("R_300_nld");
  plot.RatioBinomial(h_Ntag_300_nld, h_Ntag_300, R_300_nld);
  TH1D *R_400_nld  = (TH1D*) h_Ntag_400_nld->Clone("R_400_nld");
  plot.RatioBinomial(h_Ntag_400_nld, h_Ntag_400, R_400_nld);
  TH1D *R_500_nld  = (TH1D*) h_Ntag_500_nld->Clone("R_500_nld");
  plot.RatioBinomial(h_Ntag_500_nld, h_Ntag_500, R_500_nld);
  TH1D *R_600_nld  = (TH1D*) h_Ntag_600_nld->Clone("R_600_nld");
  plot.RatioBinomial(h_Ntag_600_nld, h_Ntag_600, R_600_nld);
  TH1D *R_300_nlu  = (TH1D*) h_Ntag_300_nlu->Clone("R_300_nlu");
  plot.RatioBinomial(h_Ntag_300_nlu, h_Ntag_300, R_300_nlu);
  TH1D *R_400_nlu  = (TH1D*) h_Ntag_400_nlu->Clone("R_400_nlu");
  plot.RatioBinomial(h_Ntag_400_nlu, h_Ntag_400, R_400_nlu);
  TH1D *R_500_nlu  = (TH1D*) h_Ntag_500_nlu->Clone("R_500_nlu");
  plot.RatioBinomial(h_Ntag_500_nlu, h_Ntag_500, R_500_nlu);
  TH1D *R_600_nlu  = (TH1D*) h_Ntag_600_nlu->Clone("R_600_nlu");
  plot.RatioBinomial(h_Ntag_600_nlu, h_Ntag_600, R_600_nlu);



//
////////////////////////////////////////////////////////////////////////////
//  plot results

  pt_low = 0;
// simulated pi0, gamma, and electron sectra
  TCanvas *c1 = plot.Canvas ("c1",400,500,10,10,1,0);
  plot.SetyTitleOffset(1.4);
  TH1D *frame1 = plot.Frame("frame1","p_{T} [GeV/c]","counts/#pi^{0}",pt_low,pt_high*.999,1e-12,1.);
  frame1->Draw();
//  h_ptpi0_t->Draw("sameL Chist");
  h_ptpi0->Draw("sameL Chist");
  h_ptpi02->Draw("sameL Chist");
  h_ptg->Draw("sameL Chist");
  h_pte->Draw("sameL Chist");

  TLegend *L1 = plot.Legend("#pi^{0} decay simulation ",0.5,0.75,0.88,0.95,1.);
  L1->AddEntry(h_ptpi0,"parent #pi^{0}","l");
  L1->AddEntry(h_ptg,  "decay #gamma " ,"l");
  L1->AddEntry(h_pte,  "conversion electrons","l");
  L1->Draw("same");

// various Nincl and Ntag
  TCanvas *c2 = plot.Canvas ("c2",400,500,410,10,1,0);
  TH1D *frame2 = plot.Frame("frame2","p_{T} [GeV/c]","counts/#pi^{0}",pt_low,pt_high*.999,1e-11,1e-3);
  frame2->Draw();
  h_Nincl_t->Draw("sameL Chist");
  h_Nincl->Draw("sameL Chist");
  h_Ntag_ideal->Draw("sameL Chist");
  h_Ntag_fiducial->Draw("sameL Chist");
  h_Ntag_deadarea->Draw("sameL Chist");
  h_Ntag_VTXconv->Draw("sameL Chist");

  TLegend *L2 = plot.Legend(" ",0.5,0.75,0.88,0.95,1.);
  L2->AddEntry(h_Nincl_t,"N_{incl} true momentum","l");
  L2->AddEntry(h_Nincl,  "N_{incl}" ,"p");
  L2->Draw("same");


  pt_high = 6.;
// <ef>
  plot.SetLeftMargin(0.1);  
  TCanvas *c3 = plot.Canvas ("c3",700,450,810,10);
  plot.SetyTitleOffset(.8);
  TH1D *frame3 = plot.Frame("frame3","p_{T} [GeV/c]","<#varepsilonf>",pt_low,pt_high*.999,0.,1.099);
  frame3->Draw();
  h_ef_nincl->Draw("sameP");
  h_ef_ideal->Draw("sameP");
  h_ef_fiducial->Draw("sameP");
  h_ef_deadarea->Draw("sameP");
  h_ef_VTXconv->Draw("sameP");

  TLegend *L3a = plot.Legend("",0.12,.75,.42,.86);
  L3a->AddEntry(h_ef_nincl,"N_{incl} true/reco","p");
  L3a->Draw("same");
  TLegend *L3 = plot.Legend("Conditional Acceptace for #pi^{0} #rightarrow #gamma#gamma",0.52,.18,.95,.46);
  L3->AddEntry(h_ef_ideal,"ideal PHENIX acceptance","p");
  L3->AddEntry(h_ef_fiducial,"+ fiducial cut","p");
  L3->AddEntry(h_ef_deadarea,"+ dead area","p");
  L3->AddEntry(h_ef_VTXconv,"+ VTX conversion loss","p");
  L3->Draw("same");

// <ef>
  TCanvas *c4 = plot.Canvas ("c4",700,450,810,460);
  TH1D *frame4 = plot.Frame("frame4","p_{T} [GeV/c]","<#varepsilonf>",pt_low,pt_high*.999,0.,0.599);
  frame4->Draw();
  h_ef_300->Draw("sameP");
  h_ef_400->Draw("sameP");
  h_ef_500->Draw("sameP");
  h_ef_600->Draw("sameP");
  h_ef_300_pt->Draw("sameP");
  h_ef_400_pt->Draw("sameP");
  h_ef_500_pt->Draw("sameP");
  h_ef_600_pt->Draw("sameP");

  TLegend *L4 = plot.Legend("Conditional Acceptace for #pi^{0} #rightarrow #gamma#gamma with E_{cut}",0.2,.67,.6,.95);
  L4->AddEntry(h_ef_300,"300 MeV","p");
  L4->AddEntry(h_ef_400,"400 MeV","p");
  L4->AddEntry(h_ef_500,"500 MeV","p");
  L4->AddEntry(h_ef_600,"600 MeV","p");
  L4->Draw("same");

  TCanvas *c5 = plot.Canvas ("c5",700,450,110,460);
  TH1D *frame5 = plot.Frame("frame5","pt [GeV/c]","#sigma_{<#varepsilonf>}",0.,4.99*.999,0.901,1.099);
  frame5->Draw();

  R_300_d->Draw("sameP");
  R_400_d->Draw("sameP");
  R_500_d->Draw("sameP");
  R_600_d->Draw("sameP");
  R_300_u->Draw("sameP");
  R_400_u->Draw("sameP");
  R_500_u->Draw("sameP");
  R_600_u->Draw("sameP");

  TLegend *L5 = plot.Legend("#pm 2% energy scale",0.7,.67,.9,.95);
  L5->AddEntry(R_300_d,"300 MeV","p");
  L5->AddEntry(R_400_d,"400 MeV","p");
  L5->AddEntry(R_500_d,"500 MeV","p");
  L5->AddEntry(R_600_d,"600 MeV","p");
  L5->Draw("same");

 
  TCanvas *c6 = plot.Canvas ("c6",700,450,110,10);
  TH1D *frame6 = plot.Frame("frame6","pt [GeV/c]","#sigma_{<#varepsilonf>}",0.,4.99*.999,0.96,1.0399);
  frame6->Draw();

  R_300_rd->Draw("sameP");
  R_400_rd->Draw("sameP");
  R_500_rd->Draw("sameP");
  R_600_rd->Draw("sameP");
  R_300_ru->Draw("sameP");
  R_400_ru->Draw("sameP");
  R_500_ru->Draw("sameP");
  R_600_ru->Draw("sameP");

  TLegend *L6 = plot.Legend("#pm 10% energy resolution",0.7,.67,.9,.95);
  L6->AddEntry(R_300_rd,"300 MeV","p");
  L6->AddEntry(R_400_rd,"400 MeV","p");
  L6->AddEntry(R_500_rd,"500 MeV","p");
  L6->AddEntry(R_600_rd,"600 MeV","p");
  L6->Draw("same");

  TCanvas *c7 = plot.Canvas ("c7",700,450,360,10);
  TH1D *frame7 = plot.Frame("frame6","pt [GeV/c]","#sigma_{<#varepsilonf>}",0.,4.99*.999,0.9,1.1);
  frame7->Draw();

  R_300_nld->Draw("sameP");
  R_400_nld->Draw("sameP");
  R_500_nld->Draw("sameP");
  R_600_nld->Draw("sameP");
  R_300_nlu->Draw("sameP");
  R_400_nlu->Draw("sameP");
  R_500_nlu->Draw("sameP");
  R_600_nlu->Draw("sameP");

  TLegend *L7 = plot.Legend("using non linearity correction",0.7,.67,.9,.95);
  L7->AddEntry(R_300_nld,"300 MeV","p");
  L7->AddEntry(R_400_nld,"400 MeV","p");
  L7->AddEntry(R_500_nld,"500 MeV","p");
  L7->AddEntry(R_600_nld,"600 MeV","p");
  L7->Draw("same");

  TCanvas *c8 = plot.Canvas ("c8",700,450,360,460);
  TH1D *frame8 = plot.Frame("frame8","pt [GeV/c]","#sigma_{<#varepsilonf>}",0.,4.99*.999,0.8,1.2);
  frame8->Draw();

  R_300_pt->Draw("sameP");
  R_400_pt->Draw("sameP");
  R_500_pt->Draw("sameP");
  R_600_pt->Draw("sameP");

  TLegend *L8 = plot.Legend("using differnt pt weights",0.7,.67,.9,.95);
  L8->AddEntry(R_300_pt,"300 MeV","p");
  L8->AddEntry(R_400_pt,"400 MeV","p");
  L8->AddEntry(R_500_pt,"500 MeV","p");
  L8->AddEntry(R_600_pt,"600 MeV","p");
  L8->Draw("same");



}