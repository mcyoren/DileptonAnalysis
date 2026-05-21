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
#include "TTree.h"
#include "TBranch.h"

using namespace std;

struct MyEvent {
    int ntracks;
    std::vector<int> *pid;
    std::vector<double> *mass;
    std::vector<double> *energy;
    std::vector<double> *px;
    std::vector<double> *py;
    std::vector<double> *pz;
    std::vector<double> *vx;
    std::vector<double> *vy;
    std::vector<double> *vz;

    void set_to_null() {
        ntracks = 0;
        pid = 0;
        mass = 0;
        energy = 0;
        px = 0;
        py = 0;
        pz = 0;
        vx = 0;
        vy = 0;
        vz = 0;
        return;
    };
};

void WriteROOT2OscarPythia(const TString filepath = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
						   const TString infile = "single_pi0_HELIOS_1B.root", 
						   const TString output = "oscar.txt"){

	TFile* input = new TFile(infile,"READ");
	if(!(input))
	{
	  cout << "no input file" << endl;
	  exit(1);
	}

	//vertex staff
	const double scale = 1e13; //cm to fm conversion
 	std::ifstream myfile (filepath);
 	vector<double> vertexes;
 	string line;
 	if (myfile.is_open())
 	{
 	  while ( getline (myfile,line) )
 	  {
 	    string s;
 	    std::stringstream ss(line);
 	    while(getline(ss, s, ' '))
 	    {
 	      vertexes.push_back(atof(s.c_str())*scale);
 	    }
 	  }
 	  myfile.close();
 	}
	if(false)
	{
		for (int i = 0; i < (int) vertexes.size()/4; i++)
  		{
  		  std::cout<<vertexes[4*i] << "     " <<vertexes[4*i+1] << "     " <<vertexes[4*i+2] << " "<<vertexes[4*i+3] << " " <<std::endl;
  		}
	}
	// ------------end for vertexes---------------


	ofstream file(output);

	//Read in the TTrees 

	TTree* T = (TTree*)input->Get("T");
  	const int max_trks = 200;
  	int ntracks = -999;
	MyEvent myevent;
	myevent.set_to_null();

	TBranch *branch[8];
	T->SetBranchAddress("ntracks",&ntracks);
	T->SetBranchAddress("pid",&myevent.pid, &branch[0]);
	T->SetBranchAddress("px",&myevent.px, &branch[1]);
	T->SetBranchAddress("py",&myevent.py, &branch[2]);
	T->SetBranchAddress("pz",&myevent.pz, &branch[3]);
	T->SetBranchAddress("energy",&myevent.energy, &branch[4]);
	T->SetBranchAddress("mass",&myevent.mass, &branch[5]);
	T->SetBranchAddress("vx",&myevent.vx, &branch[6]);
	T->SetBranchAddress("vy",&myevent.vy, &branch[7]);
	T->SetBranchAddress("vz",&myevent.vz, &branch[8]);

	file << "# OSC1999A" << endl;
	file << "# final_id_p_x" << endl;
	file << "# SimName 1.0" << endl;
	file << "#" << endl;
	file << "# Some comments..." << endl;
	file << endl;

	for(int ievt = 0; ievt < (int)T->GetEntries(); ievt++){
		myevent.set_to_null();
	  T->GetEntry(ievt);
       
	  file << 0 << "\t" << ntracks << endl;

	  for(int i = 0; i < ntracks; i++){
		  
	    if(i == -1) file << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
	    else file << i+1 << "\t" << myevent.pid->at(i) << "\t" << 0 << "\t" << myevent.px->at(i) << "\t" << myevent.py->at(i) << "\t" <<
		 							myevent.pz->at(i) << "\t" << myevent.energy->at(i) << "\t" << myevent.mass->at(i) << "\t" << 
									myevent.vx->at(i)*pow(10,12)+vertexes[4*ievt] << "\t" << myevent.vy->at(i)*pow(10,12)+vertexes[4*ievt+1]  << "\t" << 
									myevent.vz->at(i)*pow(10,12) +vertexes[4*ievt+2] << "\t" << 0 << endl;
		}

		file << 0 << "\t" << 0 << endl;
	}

}
