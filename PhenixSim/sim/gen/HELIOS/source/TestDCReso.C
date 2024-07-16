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

///////////////////////////////////////////////////////////////////////////////
void TestDCReso(){
//
  TRandom3 MyRandy = TRandom3(0);                // Random Generator

  cout << endl;

  PHENIXDetector PHENIX;
  MyPlot plot;

  Particle pip("pi+");
  Particle pim("pi-");

  TLorentzVector Reco1, Reco2;

  Double_t pt1,phi1,eta1,ww1;
  Double_t pt2,phi2,eta2,ww2;
  Int_t q1,q2;
   
  Double_t pt_min = 0.1;
  Double_t pt_max = 20.;
  Double_t pt_low = 0.1;
  Double_t pt_high = 15.;
  Double_t ptcut = 0.2;
  Int_t nevt = 1000000;

  Bool_t bp = false;
  Bool_t bm = false;


  TF1 *piHagedorn      = Hagedorn("piHagedorn", piMass, pt_max, pt_min);

  TH1D *h_ptpi             = new TH1D("h_ptpi","pt",100,pt_low,pt_high);
  plot.StyleMe(h_ptpi,   20, kBlue, 1., 1, 2); 
  TH1D *h_ptpireco         = new TH1D("h_ptpireco","ptreco",100,pt_low,pt_high);
  plot.StyleMe(h_ptpireco, 20, kTeal, 1., 1, 2); 

  TH2D *h_butsyk           = new TH2D("h_butsyk","butsyk",100,-2*pi,2*pi,100,-6.0,6.0);
  TH2D *h_areso            = new TH2D("h_areso","areso",250,-6.,6.0,250,-6.0,6.0);
  TH2D *h_ptreso            = new TH2D("h_ptreso","ptreso",250,pt_low,pt_max,250,pt_low,2*pt_max);
  TH2D *h_ptreso2           = new TH2D("h_ptreso2","ptreso2",250,pt_low,pt_max,250,-0.3,.3);
  TH2D *h_phireso           = new TH2D("h_phireso","phireso",250,pt_low,pt_max,250,-0.1,.1);
  TH2D *h_etareso           = new TH2D("h_etareso","etareso",250,pt_low,pt_max,250,-0.1,.1);
  

  for (int i=0; i<nevt; i++){                                     // event loop
    if(i%100000==0) cout << "event loop at " << i << endl;         // event counter printed 
/////////////////////////////////////////////// pi decays ////////////////////////////////////////////////
    bp = false;
    bm = false;

    pip.ResetP();
    pip.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi with flat distribution
    pip.SetWeight(piHagedorn->Eval(pip.Pt())/float(nevt));
    pim.ResetP();
    pim.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi with flat distribution
    pim.SetWeight(piHagedorn->Eval(pim.Pt())/float(nevt));

    if (i<10) cout << " PHENIX Acceptance " << PHENIX.InAcceptance(pip,pip.Charge()) << " " << pip.Phi() << endl;
    h_ptpi->Fill(pip.Pt(),pip.Weight());
    h_ptpi->Fill(pim.Pt(),pim.Weight());
 
    
    if (PHENIX.InAcceptance(pip,pip.Charge())>0) {                                // monitor PEHNIX acceptance 
              bp = true;
              h_butsyk->Fill(pip.Phi(),1./pip.Pt(),pip.Weight());                   //      for charged particles
    }

    if (PHENIX.InAcceptance(pim,pim.Charge())>0) {                                // monitor PEHNIX acceptance 
              bm = true;
              h_butsyk->Fill(pim.Phi(),-1./pim.Pt(),pim.Weight());                   //      for charged particles
    }

    if (i<50) cout << pip.Pt() << "  " << pt1 << endl;

   
    if (bp){
      q1 = pip.Charge();
      PHENIX.CharacterizeTrack(pip,pip.Charge(),pipID);
      Reco1 = PHENIX.ReconstructTrack(pip, pip.Charge());
      if (i<100) { 
        cout << endl;
        cout << pip.ID() << endl;
        pip.Print();
        Reco1.Print();
      }
      pt1 = Reco1.Pt();
      phi1 = Reco1.Phi();
      eta1 = Reco1.Eta();
      h_ptpireco->Fill(pt1,pip.Weight());
      h_areso->Fill(1./pip.Pt(),1./pt1);
      h_ptreso->Fill(pip.Pt(),pt1);
      h_ptreso2->Fill(pip.Pt(),(pt1-pip.Pt())/pip.Pt(),pip.Weight());
      h_phireso->Fill(pip.Pt(),(phi1-pip.Phi()),pip.Weight());
      h_etareso->Fill(pip.Pt(),(eta1-pip.Eta()),pip.Weight());

      // if (abs(phi1-pip.Phi())>0.2){
      //   cout << pip.Pt() << "   " << pip.Phi() << "   " << phi1 
      //   << " PHENIX Acceptance " << PHENIX.InAcceptance(pip,pip.Charge()) << endl; 
      // }
    }
    if (bm) {
      q2 = pim.Charge();
      PHENIX.CharacterizeTrack(pim,pim.Charge(), pimID);
      Reco2 = PHENIX.ReconstructTrack(pim, pim.Charge());
      pt2 = Reco2.Pt();
      phi2 = Reco2.Phi();
      eta2 = Reco2.Eta();
      h_areso->Fill(pim.Charge()/pim.Pt(),-1./pt2);
      h_ptreso->Fill(pim.Pt(),pt2);
      h_ptreso2->Fill(pim.Pt(),(pt2-pim.Pt())/pim.Pt(),pim.Weight());
      h_ptpireco->Fill(pt2,pim.Weight());
      h_phireso->Fill(pim.Pt(),(phi2-pim.Phi()),pim.Weight());
      h_etareso->Fill(pim.Pt(),(eta2-pim.Eta()),pim.Weight());
    }
  }

  h_ptreso2->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_reso_mean = (TH1D*)gDirectory->Get("h_ptreso2_1");
  h_reso_mean->SetName("h_reso_mean");
  plot.StyleMe(h_reso_mean,20,kTeal,.8);
  TH1D * h_reso_sigma = (TH1D*)gDirectory->Get("h_ptreso2_2");
  h_reso_sigma->SetName("h_reso_sigma");
  plot.StyleMe(h_reso_sigma,20,kTeal,.8);

  h_phireso->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_phi_mean = (TH1D*)gDirectory->Get("h_phireso_1");
  h_phi_mean->SetName("h_phi_mean");
  plot.StyleMe(h_phi_mean,20,kOrange+2l,.8);
  TH1D * h_phi_sigma = (TH1D*)gDirectory->Get("h_phireso_2");
  h_phi_sigma->SetName("h_phi_sigma");
  plot.StyleMe(h_phi_sigma,20,kOrange+2,.8);

  h_etareso->FitSlicesY(0,0,-1,0,"Q");                            // extract mean and sigma of pi0 mass vs pt
  TH1D * h_eta_mean = (TH1D*)gDirectory->Get("h_etareso_1");
  h_eta_mean->SetName("h_eta_mean");
  plot.StyleMe(h_eta_mean,20,kGreen+2l,.8);
  TH1D * h_eta_sigma = (TH1D*)gDirectory->Get("h_etareso_2");
  h_eta_sigma->SetName("h_eta_sigma");
  plot.StyleMe(h_eta_sigma,20,kGreen+2,.8);

 ////////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = plot.Canvas ("c1",400,500,10,10,1);
  plot.SetyTitleOffset(1.3);
  TH1D *frame1 = plot.Frame("frame1","pt","counts",pt_low,pt_high*.999,1e-12,1.);
  frame1->Draw();
  h_ptpi->Draw("sameL Chist");
  h_ptpireco->Draw("sameL Chist");

  plot.SetyTitleOffset(0.6);
  plot.SetLeftMargin(0.1);
  TCanvas *c2 = plot.Canvas ("c2",600,300,410,10);
  TH1D *frame2 = plot.Frame("frame2","phi","q/pt",-2.*pi,2*pi,-6.,6.);
  frame2->Draw();
  h_butsyk->Draw("Col same");
  plot.Reset();
 
  TCanvas *c3 = plot.Canvas ("c3",400,400,10,510);
  TH1D *frame3 = plot.Frame("frame3","q/pt true","q/pt reco",-6.,6.,-6.,6.);
  frame3->Draw();
  h_areso->Draw("Col same");

  TCanvas *c4 = plot.Canvas ("c4",400,400,110,510);
  TH1D *frame4 = plot.Frame("frame4","pt true","pt reco",0.,15.,0.,30.);
  frame4->Draw();
  h_ptreso->Draw("Col same");

  TCanvas *c5 = plot.Canvas ("c5",400,400,210,510);
  TH1D *frame5 = plot.Frame("frame5","pt true","pt reco-true",0.,15.,-0.3,.3);
  frame5->Draw();
  h_ptreso2->Draw("Col same");
  h_reso_mean->Draw("same");

  TCanvas *c6 = plot.Canvas ("c6",400,400,310,510);
  plot.SetyTitleOffset(1.2);
  TH1D *frame6 = plot.Frame("frame6","pt true","Dphi ",0.,5.,-0.1,.1);
  frame6->Draw();
  h_phireso->Draw("Col same");
  h_phi_mean->Draw("same");

  TCanvas *c7 = plot.Canvas ("c7",400,400,410,510);
  TH1D *frame7 = plot.Frame("frame7","pt true","Deta ",0.,5.,-0.1,.1);
  frame7->Draw();
  h_etareso->Draw("Col same");
  h_phi_mean->Draw("same");
  plot.Reset();

  TCanvas *c10 = plot.Canvas ("c10",500,300,1010,10);
  TH1D *frame10 = plot.Frame("frame10","pt true","#sigma_{pt} (%)",0.,5.,0.,.06);
  frame10->Draw();
  h_reso_sigma->Draw("same");

  TCanvas *c11 = plot.Canvas ("c11",500,300,1010,310);
  TH1D *frame11 = plot.Frame("frame11","pt true","#sigma_{#phi} (rad)",0.,5.,0.,.03499);
  frame11->Draw();
  h_phi_sigma->Draw("same");
 
  TCanvas *c12 = plot.Canvas ("c12",500,300,1010,610);
  TH1D *frame12 = plot.Frame("frame12","pt true","#sigma_{#eta} (rad)",0.,5.,0.,.03499);
  frame12->Draw();
  h_eta_sigma->Draw("same");
 
}




