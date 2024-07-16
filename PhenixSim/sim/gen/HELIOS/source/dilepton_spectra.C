#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include "HELIOSLibrary/HELIOSLibrary.h"
#include "MyPlotting/MyPlot.C"


using namespace std;

double TsallisFit(double *x, double *par){

  double pt = x[0];
  double m0 = par[0];
  double yield = par[1];
  double T = par[2];
  double n = par[3];

  double mt = sqrt(pt*pt + m0*m0);
  return (pt/42.)*(yield)*(((n-1)*(n-2))/((n*T+m0*(n-1))*(n*T+m0)))*pow(((n*T+mt)/(n*T+m0)), -n);

}

TF1 *TsallisYield(const char *name, double mass, double a, double b, double c){

  TF1 *fTsallis = new TF1(name, TsallisFit, 0, 10, 4);
  fTsallis->SetParameters(mass, a, b, c);
  fTsallis->SetParNames("mass", "ds/dy", "T", "n");
  return fTsallis;
}


void dilepton_spectra(){

  TRandom3 MyRandy = TRandom3(0);                // Random Generator
  TLorentzVector temp;
  Double_t ww = 0;
  Int_t id = 0;

  Particle pi0("pi0");
  Particle eta("eta");
  Particle omega("omega");
  Particle etap("etap");
  Particle phi_meson("phi");
  Particle jpsi_meson("jpsi");
  Particle psip_meson("psip");
  Particle electron("electron");
  Particle positron("positron");
  TLorentzVector eepair;

  TH1D *h_masspi = new TH1D("hmasspi","mass",300,0.,4.);
  TH1D *h_masseta = new TH1D("hmasseta","mass",300,0.,4.);
  TH1D *h_massetap = new TH1D("hmassetap","mass",300,0.,4.0);
  TH1D *h_massomega = new TH1D("hmassomega","mass",300,0.,4.0);
  TH1D *h_massphi = new TH1D("hmassphi","mass",300,0.,4.0);
  TH1D *h_massjpsi = new TH1D("hmassjpsi","mass",300,0.,4.0);
  TH1D *h_masspsip = new TH1D("hmasspsip","mass",300,0.,4.0);

  TH1D *h_ptpi = new TH1D("hptpi","pt",100,0.,10.);
  TH1D *h_pteta = new TH1D("hpteta","pt",100,0.,10.);
  TH1D *h_ptetap = new TH1D("hptetap","pt",100,0.,10.0);
  TH1D *h_ptomega = new TH1D("hptomega","pt",100,0.,10.0);
  TH1D *h_ptphi = new TH1D("hptphi","pt",100,0.,10.0);
  TH1D *h_ptjpsi = new TH1D("hptjpsi","pt",100,0.,10.0);
  TH1D *h_ptpsip = new TH1D("hptpsip","pt",100,0.,10.0);

  TH1D *h_mass = new TH1D("hmass","mass",300,0.,4.0);

  Bool_t elecFound = false;
  Bool_t posiFound = false;

  Int_t nevt = 100000;
  int ndecay = 0;

  TF1 *pi0Tsallis     = TsallisYield("pi0Tsallis", pi0Mass, 40.5, 0.1142, 9.57);
  TF1 *etaTsallis     = TsallisYield("etaTsallis", etaMass, 4.46, 0.119, 9.84);
  TF1 *omegaTsallis   = TsallisYield("omegaTsallis", omegaMass, 3.64, 0.1218, 10.);
  TF1 *etapTsallis    = TsallisYield("etapTsallis", etapMass, 0.62, 0.1238, 10.11);
  TF1 *phiTsallis     = TsallisYield("phiTsallis", phiMass, 0.421, 0.1245, 10.15);
  TF1 *jpsiTsallis    = TsallisYield("jpsiTsallis", jpsiMass, 0.000761, 0.149, 11.5);
  TF1 *psipTsallis    = TsallisYield("psipTsallis", psipMass, 0.000133, 0.156, 11.9);

  for (int i=0; i<nevt; i++){

    if(i%100000==0) cout << "event loop at " << i << endl;
    pi0.ResetP();
    pi0.GenerateP(0.2,10.);
    pi0.SetWeight(pi0Tsallis->Eval(pi0.Pt()));
    h_ptpi->Fill(pi0.Pt(),pi0.Weight()); 
    pi0.DecaySingleBranch("pi0->gee");

    ndecay = pi0.GetNumberOfDaughters();

    elecFound = false;
    posiFound = false;

    for (Int_t j=0; j< ndecay; j++) {

        temp = pi0.GetDecayDaughter(j);
        ww   = pi0.GetDaughterWeight(j)*pi0.Weight();
        id   = pi0.GetDaughterID(j);
        if (id == electronID) {
            electron.UpdateParticle(temp);
            electron.SetWeight(ww);
            elecFound = true;
        } 
        else if (id == positronID) {
            positron.UpdateParticle(temp);
            positron.SetWeight(ww);
            posiFound = true;
        }
    }

    if (elecFound && posiFound) {

      eepair = electron+positron;
      h_masspi->Fill(eepair.M(),ww);

    }

    eta.ResetP();
    electron.ResetP();
    positron.ResetP();
    eta.GenerateP(0.2,10);                                // generate 4 vector for pi0 with flat distribution
    eta.SetWeight(etaTsallis->Eval(eta.Pt()));
    h_pteta->Fill(eta.Pt(),eta.Weight());
    eta.DecaySingleBranch("eta->gee");

    ndecay = eta.GetNumberOfDaughters();

    elecFound = false;
    posiFound = false;    

    for (Int_t j=0; j< ndecay; j++) {

      temp = eta.GetDecayDaughter(j);
      ww   = eta.GetDaughterWeight(j)*eta.Weight();
      id   = eta.GetDaughterID(j);
      
      if (id == electronID) {

            electron.UpdateParticle(temp);
            electron.SetWeight(ww);
            elecFound = true;
        } 
      
      else if (id == positronID) {

            positron.UpdateParticle(temp);
            positron.SetWeight(ww);
            posiFound = true;
        }
    }

    if (elecFound && posiFound) {
      eepair = electron+positron;
      h_masseta->Fill(eepair.M(),ww);     
    }

    omega.ResetP();
    electron.ResetP();
    positron.ResetP();
    omega.GenerateP(0.2,10);
    omega.SetWeight(omegaTsallis->Eval(omega.Pt()));
    h_ptomega->Fill(omega.Pt(),omega.Weight());
    omega.DecaySingleBranch("omega->ee");

    ndecay = omega.GetNumberOfDaughters();

    elecFound = false;
    posiFound = false;    
      
    for (Int_t j=0; j< ndecay; j++) {

       temp = omega.GetDecayDaughter(j);
       ww   = omega.GetDaughterWeight(j)*omega.Weight();
       id   = omega.GetDaughterID(j);
       
       if (id == electronID) {

            electron.UpdateParticle(temp);
            electron.SetWeight(ww);
            elecFound = true;
       } 
       else if (id == positronID) {
            positron.UpdateParticle(temp);
            positron.SetWeight(ww);
            posiFound = true;
        }
    } 

    if (elecFound && posiFound) {
       eepair = electron+positron;
       h_massomega->Fill(eepair.M(),ww);
    }

    etap.ResetP();
    electron.ResetP();
    positron.ResetP();
    etap.GenerateP(0.2, 10);                                // generate 4 vector for pi0 with flat distribution
    etap.SetWeight(etapTsallis->Eval(etap.Pt()));
    h_ptetap->Fill(etap.Pt(),etap.Weight());
    etap.DecaySingleBranch("etap->gee");

    ndecay = etap.GetNumberOfDaughters();

    elecFound = false;
    posiFound = false;    

    for (Int_t j=0; j< ndecay; j++) {
      temp = etap.GetDecayDaughter(j);
      ww   = etap.GetDaughterWeight(j)*etap.Weight();
      id   = etap.GetDaughterID(j);
        if (id == electronID) {
            electron.UpdateParticle(temp);
            electron.SetWeight(ww);
            elecFound = true;
        } 
        else if (id == positronID) {
            positron.UpdateParticle(temp);
            positron.SetWeight(ww);
            posiFound = true;
        }
    }

    if (elecFound && posiFound) {
      eepair = electron+positron;
      h_massetap->Fill(eepair.M(),ww);     
    }

    phi_meson.ResetP();
    electron.ResetP();
    positron.ResetP();
    phi_meson.GenerateP(0.2, 10.);
    phi_meson.SetWeight(phiTsallis->Eval(phi_meson.Pt()));
    h_ptphi->Fill(phi_meson.Pt(),phi_meson.Weight());
    phi_meson.DecaySingleBranch("phi->ee");

    ndecay = phi_meson.GetNumberOfDaughters();

    elecFound = false;
    posiFound = false;    

   for (Int_t j=0; j< ndecay; j++) {

     temp = phi_meson.GetDecayDaughter(j);
     ww   = phi_meson.GetDaughterWeight(j)*phi_meson.Weight();
     id   = phi_meson.GetDaughterID(j);
        if (id == electronID) {
            electron.UpdateParticle(temp);
            electron.SetWeight(ww);
            elecFound = true;
        } 
        else if (id == positronID) {
            positron.UpdateParticle(temp);
            positron.SetWeight(ww);
            posiFound = true;
       
        }
    }

  if (elecFound && posiFound) {
    eepair = electron+positron;
    h_massphi->Fill(eepair.M(),ww);     
  }

    jpsi_meson.ResetP();
    electron.ResetP();
    positron.ResetP();
    jpsi_meson.GenerateP(0.2,10);
    jpsi_meson.SetWeight(jpsiTsallis->Eval(jpsi_meson.Pt()));
    h_ptjpsi->Fill(jpsi_meson.Pt(),jpsi_meson.Weight());
    jpsi_meson.DecaySingleBranch("jpsi->ee");

    ndecay = jpsi_meson.GetNumberOfDaughters();

    elecFound = false;
    posiFound = false;    

   for (Int_t j=0; j< ndecay; j++) {

     temp = jpsi_meson.GetDecayDaughter(j);
     ww   = jpsi_meson.GetDaughterWeight(j)*jpsi_meson.Weight();
     id   = jpsi_meson.GetDaughterID(j);

     if (id == electronID) {
       electron.UpdateParticle(temp);
       electron.SetWeight(ww);
       elecFound = true;
     } 
     else if (id == positronID) {
       positron.UpdateParticle(temp);
       positron.SetWeight(ww);
       posiFound = true;
     }
    }

  if (elecFound && posiFound) {
    eepair = electron+positron;
    h_massjpsi->Fill(eepair.M(),ww);     
  }

  psip_meson.ResetP();
  electron.ResetP();
  positron.ResetP();
  psip_meson.GenerateP(0.2, 10);
  psip_meson.SetWeight(psipTsallis->Eval(psip_meson.Pt()));
  h_ptpsip->Fill(psip_meson.Pt(),psip_meson.Weight());
  psip_meson.DecaySingleBranch("psip->ee");

  ndecay = psip_meson.GetNumberOfDaughters();
    
  elecFound = false;
  posiFound = false;    

  for (Int_t j=0; j< ndecay; j++) {
    temp = psip_meson.GetDecayDaughter(j);
    
    ww   = psip_meson.GetDaughterWeight(j)*psip_meson.Weight();
    id   = psip_meson.GetDaughterID(j);
    
    if (id == electronID) {
      electron.UpdateParticle(temp);
      electron.SetWeight(ww);
      elecFound = true;
    } 
    else if (id == positronID) {
      positron.UpdateParticle(temp);
      positron.SetWeight(ww);
      posiFound = true;
  
    }
  }
  
  if (elecFound && posiFound) {
    eepair = electron+positron;
    h_masspsip->Fill(eepair.M(),ww);     
  }
}

  h_mass->Add(h_masspi);
  h_mass->Add(h_masseta);
  h_mass->Add(h_massomega);
  h_mass->Add(h_massetap);
  h_mass->Add(h_massjpsi);
  h_mass->Add(h_masspsip);
  h_mass->Add(h_massphi);

  TFile* fout = new TFile("mass_spectrum_wo_weighting.root","RECREATE");
  fout->cd();
  h_mass->Write();
  h_masspi->Write();
  h_masseta->Write();
  h_massomega->Write();
  h_massetap->Write();
  h_massphi->Write();
  h_massjpsi->Write();
  h_masspsip->Write();
  
  h_ptpi->Write();
  h_pteta->Write();
  h_ptomega->Write();
  h_ptetap->Write();
  h_ptphi->Write();
  h_ptjpsi->Write();
  h_ptpsip->Write();
  fout->Close();
}






