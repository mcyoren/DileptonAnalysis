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

#include "/phenix/plhf/mitran/Simul/Dileptons/sim/gen/HELIOS/source/HELIOSLibrary/WriteEvent.h"

using namespace std;

void WriteROOT2Oscar(const TString filepath = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
					 const int nfolder = 0, 
					 const int Nev = 10000, 
					 const TString infile = "helios_jpsi.root", 
					 const TString output = "oscar.particles.dat"
					 )
{

	TFile* input = new TFile(infile,"READ");
	if(!(input))
	{
	  cout << "no input file" << endl;
	  exit(1);
	}

	ofstream file(output);

	const int nstart = nfolder * Nev;
	const int nend = nstart + Nev;

	//Read in the TTrees 

	TTree* T = (TTree*)input->Get("T");
	TBranch* br = T->GetBranch("MyEvent");

	WriteEvent* event = 0;
	br->SetAddress(&event);

	file << "# OSC1999A" << endl;
	file << "# final_id_p_x" << endl;
	file << "# SimName 1.0" << endl;
	file << "#" << endl;
	file << "# Some comments..." << endl;
	file << endl;

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

	for(int ievt = nstart; ievt < nend; ievt++){

	  
		br->GetEntry(ievt);

		file << 0 << "\t" << event->GetNStable() << endl;
		const int ivertexevent = ievt - nstart;
		for(int i = 0; i < event->GetNEntries(); i++){

			WriteTrack Track = event->GetWriteTrack(i);

			if(Track.GetFinal() == 1){

				file << Track.GetNum() << "\t" << Track.GetID() << "\t" << 0 << "\t" << Track.GetPx() << "\t" << Track.GetPy() << "\t" << Track.GetPz() 
				<< "\t" << Track.GetEnergy() << "\t" << Track.GetMass() << "\t" 
				<< Track.GetXpos()+vertexes[4*ivertexevent] << "\t" << Track.GetYpos()+vertexes[4*ivertexevent+1] << "\t" << Track.GetZpos()+vertexes[4*ivertexevent+2] << "\t" << 0 << endl;
			}
		}

		file << 0 << "\t" << 0 << endl;
	}

}
