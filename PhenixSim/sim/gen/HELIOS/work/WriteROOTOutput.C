// Writes out ROOT output file 
// This example generates single pi0 events with the pi0 decay handled by HELIOS


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TString.h>

using namespace std;

#include "../source/HELIOSLibrary/HELIOSLibrary.h"
PHENIXDetector MyPHENIX; 

void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist, float px, float py, float pz, float en, float mass, float xpos, float ypos, float zpos, double br, int parent_index, long double weight);

TF1 *pi0Hagedorn     = HagedornYield("pi0Hagedorn",   pi0Mass,  0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *etaHagedorn     = HagedornYield("etaHagedorn",   etaMass,  0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *rho0Hagedorn    = HagedornYield("rhoHagedorn",   rho0Mass, 0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *omegaHagedorn   = HagedornYield("omegaHagedorn", omegaMass,0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *etapHagedorn    = HagedornYield("etapHagedorn",  etapMass, 0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *phiHagedorn     = HagedornYield("phiHagedorn",   phiMass,  0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *jpsiHagedorn    = HagedornYield("jpsiHagedorn",  jpsiMass, 0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);
TF1 *psipHagedorn    = HagedornYield("psipHagedorn",  psipMass, 0.2, 25, 377., 0.356, 0.068, 0.7, -8.25);

TF1 *InHagedorn[8] = {pi0Hagedorn, etaHagedorn, rho0Hagedorn, omegaHagedorn, etapHagedorn, phiHagedorn, jpsiHagedorn, psipHagedorn};
TString allpartnames[8]={"pi0","eta","rho0","omega","etap","phi","jpsi","psip"};

void WriteROOTOutput(TString OutputName="kek.root", const int Nev = 5000, const TString part_name = "phi", const TString decay_name = "phi->ee", const double pt_min = 0.2, const double pt_max = 10.){

	TRandom3 MyRandy = TRandom3(0);

	TFile* f = new TFile(OutputName,"RECREATE");

	WriteEvent MyEvent;
	WriteTrack MyTrack;

	TH2F* h_input = new TH2F("h_input", "PID vs decay", 2500, -0.5, 2499.5, 10, -0.5, 9.5);
	TH1F* h_input_pt = new TH1F("h_input_pt", "pt profile of input", 100, 0., 10.0);
	TH1F* h_input_eta = new TH1F("h_input_eta", "eta profile of input", 200, -1.0, 1.0);
	TH1F* h_input_phi = new TH1F("h_input_phi", "phi profile of input", 200, -10.0, 10.0);

	TH2D* h2d_opening_angle = new TH2D("h2d_opening_angle", "Opening Angle Distribution", 100, 0, 10, 100, 0, pi);
	TProfile* hprof_opening_angle = new TProfile("hprof_opening_angle", "Opening Angle Distribution", 100, 0, 10);
	
	TTree* T = new TTree("T","Tree with my events");
	T->Branch("MyEvent",&MyEvent);

	TLorentzVector temp, temp1, temp2;
	Int_t ndecay = 0;
	TString DecayBranch;
	Int_t DecayNumber, mother_index, daughter_index;
	Bool_t VTXconv = false;

	TF1 *func = new TF1("func", "gaus", -10, 10);
	func->SetParameter(0, 2.76901e+03);
	func->SetParameter(1, 7.84214e-02);
	func->SetParameter(2, 1.22454e1);

	Particle pi0(part_name);

	int ipart = 0;
	for (int i = 0; i < 8; i++)
	{
		if (allpartnames[i]==part_name) ipart = i;
	}
	std::cout<<"algortm suggest that the following paricle was chosen: "<<allpartnames[ipart]<<std::endl;
	
	Int_t nevt = Nev; 
	Int_t nentry = 0;

	for (int i = 0; i < nevt; i++){   									// event loop

	  	if(i%1000000 == 0) cout << "event loop at " << i << endl;     	// event counter printed 

	  	MyEvent.ClearEvent();
	  	
	  	pi0.ResetP();                                                 	// reset
	 	pi0.GenerateP(pt_min,pt_max);                                 	// generate
	 	double trk_wt = 1;

		pi0.DecaySingleBranch(decay_name);
		pi0.SetWeight((InHagedorn[ipart]->Eval(pi0.Pt())));
		ndecay = pi0.GetNumberOfDaughters();                          	// 2 for pi0->gg and 3 for pi0-e+e-g

		DecayNumber = pi0.GetDecayNumber();								// 0 for pi0->gg and 1 for pi0-e+e-g

		double zvertex = 0;//func->GetRandom() * pow(10,13); //in fm

		double phi0_A1 = -999;
		double phi0_A2 = -999;
		double parent_pt = pi0.Pt();

		setTrack(MyTrack, 0, nentry, pi0.ID(), 0, pi0.Px(), pi0.Py(), pi0.Pz(), pi0.E(), pi0.M(), 0, 0, zvertex, DecayNumber, -999, trk_wt);
		MyEvent.AddEntry(MyTrack); 
		mother_index = nentry;
		nentry++;

		h_input->Fill(pi0.ID(), DecayNumber);

		h_input_pt->Fill(pi0.Pt(), pi0.Weight()); //pi0
		h_input_eta->Fill(pi0.PseudoRapidity()); //pi0
		h_input_phi->Fill(pi0.Phi()); //pi0
    	
		for (Int_t j=0; j< ndecay; j++) {                             // loop over decay particles

    	  	temp = pi0.GetDecayDaughter(j);                             
    	  	int id = pi0.GetDaughterID(j);
    	  	trk_wt = pi0.GetDaughterWeight(j)*pi0.Weight();

    		setTrack(MyTrack, 1, nentry, id, 0, temp.Px(), temp.Py(), temp.Pz(), temp.E(), temp.M(), 0, 0, zvertex, pi0.GetDaughterWeight(j), mother_index, trk_wt);
    		MyEvent.AddEntry(MyTrack);
    		nentry++;

		if(j==0) phi0_A1 = temp.Phi();
		else phi0_A2 = temp.Phi();

		}

		h2d_opening_angle->Fill(parent_pt, fabs(phi0_A1 - phi0_A2));
		hprof_opening_angle->Fill(parent_pt, fabs(phi0_A1 - phi0_A2));

    	MyEvent.SetNStable(ndecay);
    	T->Fill();

    	MyEvent.ClearEvent();
    	nentry = 0;
    }

    f->cd();

    h_input->Write();
    h_input_pt->Write();
    h_input_eta->Write();
    h_input_phi->Write();
	h2d_opening_angle->Write();
	hprof_opening_angle->Write();

    T->Write();
    f->Close();
}


void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist, float px, float py, float pz, float en, float mass, float xpos, float ypos, float zpos, double br, int parent_index, long double weight){

	newTrack.SetFinal(isFinal);
	newTrack.SetNum(num);
	newTrack.SetID(id);
	newTrack.SetIst(ist);
	newTrack.SetPx(px);
	newTrack.SetPy(py);
	newTrack.SetPz(pz);
	newTrack.SetEnergy(en);
	newTrack.SetMass(mass);
	newTrack.SetXpos(xpos);
	newTrack.SetYpos(ypos);
	newTrack.SetZpos(zpos);
	newTrack.SetBranch(br);
	newTrack.SetParentIndex(parent_index);
	newTrack.SetWeight(weight);

}

