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
#include <TVector.h>

using namespace std;

void WriteROOT2OscarPythia(TString infile = "single_pi0_HELIOS_1B.root", TString output = "oscar.txt"){

	TFile* input = new TFile(infile,"READ");
	if(!(input))
	{
	  cout << "no input file" << endl;
	  exit(1);
	}

	ofstream file(output);

	//Read in the TTrees 

	TTree* T = (TTree*)input->Get("T");
  	const int max_trks = 200;
  	int ntracks = -999;
	int pid[max_trks] = {-999};
  	double px[max_trks] = {-999};
	double py[max_trks] = {-999};
	double pz[max_trks] = {-999};
	double mass[max_trks] = {-999};
	double energy[max_trks] = {-999};
	double vx[max_trks] = {-999};
	double vy[max_trks] = {-999};
	double vz[max_trks] = {-999};

	T->SetBranchAddress("ntracks",&ntracks);
	T->SetBranchAddress("pid",pid);
	T->SetBranchAddress("px",px);
	T->SetBranchAddress("py",py);
	T->SetBranchAddress("pz",pz);
	T->SetBranchAddress("mass",mass);
	T->SetBranchAddress("energy",energy);
	T->SetBranchAddress("vx",vx);
	T->SetBranchAddress("vy",vy);
	T->SetBranchAddress("vz",vz);

	file << "# OSC1999A" << endl;
	file << "# final_id_p_x" << endl;
	file << "# SimName 1.0" << endl;
	file << "#" << endl;
	file << "# Some comments..." << endl;
	file << endl;

	for(int ievt = 0; ievt < (int)T->GetEntries(); ievt++){

	  T->GetEntry(ievt);
       
	  file << 0 << "\t" << ntracks+1 << endl;

	  for(int i = -1; i < ntracks; i++){
		  
	    if(i == -1) file << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
	    else file << i+1 << "\t" << pid[i] << "\t" << 0 << "\t" << px[i] << "\t" << py[i] << "\t" << pz[i] << "\t" << energy[i] << "\t" << mass[i] << "\t" << vx[i]*pow(10,12) << "\t" << vy[i]*pow(10,12)  << "\t" << vz[i]*pow(10,12) << "\t" << 0 << endl;
		}

		file << 0 << "\t" << 0 << endl;
	}

}
