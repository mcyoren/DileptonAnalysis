// Writes out OSCAR text output from HELIOS
// This example writes out single pi0 events and also includes photon conversions

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

#include "HELIOSLibrary/HELIOSLibrary.h"

PHENIXDetector MyPHENIX;

void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist, float px, float py, float pz, float en, float mass, float xpos, float ypos, float zpos, int br);

void WriteOscarOutput(){

	char name[60];
	sprintf(name,"output.dat");
	ofstream* oscar_file = OpenOutputFile(name);

	TRandom3 MyRandy = TRandom3(0);                // Random Generator

	WriteEvent MyEvent;
	WriteTrack MyTrack;

	Particle pi0("pi0");
	TLorentzVector temp, temp1, temp2;
	Int_t ndecay = 0;
	TString DecayBranch;
	Int_t DecayNumber;
	Bool_t VTXconv = false;
	                              
	Double_t pt_min = 0.3;     // used for all simulation
	Double_t pt_max = 25.;

	Int_t nevt = 10;
	Int_t nentry = 0;

	for (int i = 0; i < nevt; i++){   									 // event loop

	  	if(i%10 == 0) cout << "event loop at " << i << endl;         // event counter printed 

	  	pi0.ResetP();                                                 // reset
	 	pi0.GenerateP(pt_min,pt_max);                                 // generate

	  	pi0.DecayFlat();                                              // generate decay photon

    	ndecay = pi0.GetNumberOfDaughters();                          // 2 for pi0->gg and 3 for pi0-e+e-g

    	DecayNumber = pi0.GetDecayNumber();
    	
    	setTrack(MyTrack, 0, nentry, pi0.ID(), 0, pi0.Px(), pi0.Py(), pi0.Pz(), pi0.E(), pi0.M(), 0, 0, 0, DecayNumber);
    	
    	MyEvent.AddEntry(MyTrack); 
    	nentry++;

    	for (Int_t j=0; j< ndecay; j++) {                             // loop over decay particles

    	  	temp = pi0.GetDecayDaughter(j);                             
    	  	int id   = pi0.GetDaughterID(j);

    	  	VTXconv = false; 
	  	    if (id == photonID) {                                       // if its a photon 
	  	        VTXconv = (MyPHENIX.VTXConversion()>0);                 // check if photon will convert 
    	  	
	    	  	if (VTXconv){                                           
	    	  	    
	    	  		setTrack(MyTrack, 0, nentry, photonID, 0, temp.Px(), temp.Py(), temp.Pz(), temp.E(), photonMass, 0, 0, 0, 1);
	    	  		MyEvent.AddEntry(MyTrack);
	    	  		nentry++;

	    	  	    PhotonConvert(temp, temp1, temp2);                   // if photon converts reconstruct e+e- 
	    	  	     
	    	  		setTrack(MyTrack, 1, nentry, electronID, 0, temp1.Px(), temp1.Py(), temp1.Pz(), temp1.E(), eMass, 0, 0, 0, 0);
	    	  		MyEvent.AddEntry(MyTrack);
	    	  		nentry++;

	    	  		setTrack(MyTrack, 1, nentry, positronID, 0, temp2.Px(), temp2.Py(), temp2.Pz(), temp2.E(), eMass, 0, 0, 0, 0);
	    	  		MyEvent.AddEntry(MyTrack);
	    	  		nentry++;

	    		} else {

	    			setTrack(MyTrack, 1, nentry, photonID, 0, temp.Px(), temp.Py(), temp.Pz(), temp.E(), photonMass, 0, 0, 0, 0);
	    			MyEvent.AddEntry(MyTrack);
	    			nentry++;
	    		}
	    	}

	    	else{

	    		setTrack(MyTrack, 1, nentry, id, 0, temp.Px(), temp.Py(), temp.Pz(), temp.E(), temp.M(), 0, 0, 0, 0);
	    		MyEvent.AddEntry(MyTrack);
	    		nentry++;

	    	}
    	}

    	WriteOscarEvent(oscar_file, MyEvent);

    	MyEvent.ClearEvent();
    	nentry = 0;
    }
}


void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist, float px, float py, float pz, float en, float mass, float xpos, float ypos, float zpos, int br){

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

}
