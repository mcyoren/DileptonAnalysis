
#include <TFile.h>
#include <TH1D.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TTree.h>
#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/WriteEvent.h"

void rapidity(){

  TH1D* hist = new TH1D("hist","eta_dist",100,-1.0,1.0);

  TFile* fhelios = new TFile("helios_pi0_gee_25GeV.root", "READ");
  TTree* t_helios = (TTree*)fhelios->Get("T");
  WriteEvent* helios_event = 0;
  TBranch* br_helios  = t_helios->GetBranch("MyEvent");
  br_helios->SetAddress(&helios_event);
  int nevts = t_helios->GetEntries();

  for(int ievent = 0; ievent < 1000000; ievent++){ 

      if(ievent%10000==0) cout << "At Event = " << ievent << endl;
      br_helios->GetEntry(ievent);
      int ntrks = helios_event->GetNStable(); 
      
      for(int itrk = 0; itrk < ntrks; itrk++){ 

	    WriteTrack daughter_track = helios_event->GetWriteTrack(itrk); 
	    if(!  (daughter_track.GetFinal() == 1 && fabs(daughter_track.GetID() == 11)) ) continue; 
	    TLorentzVector particle; 
	    particle.SetPxPyPzE(daughter_track.GetPx(), daughter_track.GetPy(), daughter_track.GetPz(), daughter_track.GetEnergy()); 
	    hist->Fill(particle.Eta());

    } 

  }

  TFile* fout = new TFile("eta_distribution_pi0_gee_25GeV.root","RECREATE");
  fout->cd();
  hist->Write();
  fout->Close();

}
