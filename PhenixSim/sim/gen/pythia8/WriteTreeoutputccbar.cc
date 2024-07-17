#include <iostream>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include "TVector3.h"
#include "Pythia8/Pythia.h"
#include <chrono>
#include <ctime>

struct MyEvent {
    int ntracks;
    std::vector<int> pid;
    std::vector<double> mass;
    std::vector<double> energy;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;

    void set_to_null() {
        ntracks = 0;
        pid.clear();
        mass.clear();
        energy.clear();
        px.clear();
        py.clear();
        pz.clear();
        vx.clear();
        vy.clear();
        vz.clear();
        return;
    };
};

using namespace Pythia8;

const int MaxTracks = 10;
const double pi = TMath::ACos(-1);

int main(int argc, char* argv[]){
	
  std::cout<<"lets begin"<<std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  std::string str_seed = argv[1];
  std::string InnEv = argv[2];
  int nevents = std::stoi(InnEv);


  bool IsWriteOscar = true;
  int nevt = 20000;
  Pythia pythia;
  std::cout<<"pythia was initialized"<<std::endl;

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");
  
  pythia.readString("HardQCD:hardccbar = on");  

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + str_seed);
  pythia.readString("Next:numberCount = 100000");
 
  //Pythia8 tune for STAR (2110.09447 )
  pythia.readString("PDF:pSet = 17");
  pythia.readString("MultipartonInteractions:ecmRef = 200");
  pythia.readString("MultipartonInteractions:bprofile = 2");

  pythia.readString("MultipartonInteractions:pT0Ref = 1.40"); //Gaussian kT term
  pythia.readString("MultipartonInteractions:ecmPow = 0.135");
  pythia.readString("MultipartonInteractions:coreRadius = 0.56");
  pythia.readString("MultipartonInteractions:coreFraction = 0.78");
  pythia.readString("ColourReconnection:range = 5.4");

  pythia.init();
  std::cout<<"pythia params were set"<<std::endl;
  //Intialize the tree 

  TTree* tree = new TTree("T","RECREATE");
  MyEvent myevent;
  int ntracks = 0;
  int pid[MaxTracks] = {-999};
  double mass[MaxTracks] = {-999};
  double px[MaxTracks] = {-999};
  double py[MaxTracks] = {-999};
  double pz[MaxTracks] = {-999};
  double energy[MaxTracks] = {-999};
  double vx[MaxTracks] = {-999};
  double vy[MaxTracks] = {-999};
  double vz[MaxTracks] = {-999};

  tree->Branch("ntracks", &myevent.ntracks,"ntracks/I");
  tree->Branch("pid",&myevent.pid);
  tree->Branch("mass",&myevent.mass);
  tree->Branch("energy",&myevent.energy);
  tree->Branch("px",&myevent.px);
  tree->Branch("py",&myevent.py);
  tree->Branch("pz",&myevent.pz);
  tree->Branch("vx",&myevent.vx);
  tree->Branch("vy",&myevent.vy);
  tree->Branch("vz",&myevent.vz);

  ofstream file("oscar.particles.dat");
	if(IsWriteOscar)
  {
    file << "# OSC1999A" << endl;
	  file << "# final_id_p_x" << endl;
	  file << "# SimName 1.0" << endl;
	  file << "#" << endl;
	  file << "# Some comments..." << endl;
  }

  int evt_count = 0;

  TH1D* hist = new TH1D("hist", "pT Distribution", 100, 0, 10);
  TH1D* pTHat = new TH1D("pTHat", "pTHat Distribution", 100, 0, 10);
  TH1D* norm = new TH1D("norm","EVENT COUNTER", 2,-0.5,1.5);

  std::cout<<"tree and histos were created"<<std::endl;

  for(int ievt = 0; ievt < nevt; ievt++){
	//std::cout<<ievt<<std::endl;
    if(evt_count == nevents) break;

    if(!pythia.next()){ nevt++; continue; }
    int nparticles = pythia.event.size();
    int npart = 0;

    norm->Fill(0);

    for(int j = 0; j < nparticles; j++){

      if (!( (pythia.event[j].daughter1() == 0 || fabs(pythia.event[pythia.event[j].daughter1()].id()) == 11 ) &&  
             (pythia.event[j].daughter2() == 0 || fabs(pythia.event[pythia.event[j].daughter2()].id()) == 11 ) && 
             fabs(pythia.event[j].eta()) < 0.5 ) ) continue;
      
      int pdg = pythia.event[j].id();
      if (! (pdg == 11 ||  pdg == -11) ) continue;
       
      double Px = pythia.event[j].px();
      double Py = pythia.event[j].py();
      double Pt = sqrt(Px*Px + Py*Py);
      if(Pt < 0.2) continue;

      int daughter_index = j;
      int mother_index = pythia.event[daughter_index].mother1();
      int mother_id = fabs(pythia.event[mother_index].id());

      bool mother_D_meson = mother_id == 411 || mother_id == 421 || mother_id == 10411 || mother_id == 10421 ||
                            mother_id == 413 || mother_id == 423 || mother_id == 10413 || mother_id == 10423 ||
                            mother_id == 20413 || mother_id == 20423 || mother_id == 415 || mother_id == 425 ||
                            mother_id == 431 || mother_id == 10431 || mother_id == 433 || mother_id == 10433 ||
	                          mother_id == 20433 || mother_id == 435;

      if(!(mother_D_meson)) continue;
      
      bool parentisB = false;

      int mother_index1 = pythia.event[mother_index].mother1();
      int mother_index2 = pythia.event[mother_index].mother2();

      std::vector<int> motherIndices;
      if(mother_index1 > 0) motherIndices.push_back(mother_index1);
      if(mother_index2 > 0) motherIndices.push_back(mother_index2);

      while( !motherIndices.empty() ) {

        int motherIndex = motherIndices.front();
        motherIndices.erase(motherIndices.begin());

        if (abs(pythia.event[motherIndex].id()) == 5) { parentisB = true; break; }

        if (pythia.event[motherIndex].mother1() > 0) motherIndices.push_back(pythia.event[motherIndex].mother1());
        if (pythia.event[motherIndex].mother2() > 0) motherIndices.push_back(pythia.event[motherIndex].mother2());

      }

      if (!parentisB) npart++;

    }

    if(npart<2){ nevt++; continue; }

    int count = 0;
    myevent.set_to_null();

    for(int j = 0; j < nparticles; j++){

      if (!( (pythia.event[j].daughter1() == 0 || fabs(pythia.event[pythia.event[j].daughter1()].id()) == 11 ) &&  
             (pythia.event[j].daughter2() == 0 || fabs(pythia.event[pythia.event[j].daughter2()].id()) == 11 ) && 
             fabs(pythia.event[j].eta()) < 0.5 ) ) continue;

      int pdg = pythia.event[j].id();
      if (! (pdg == 11 ||  pdg == -11) ) continue;
      if ( fabs(pythia.event[pythia.event[j].daughter1()].id()) == 11  || fabs(pythia.event[pythia.event[j].daughter2()].id()) == 11  ) std::cout<<"keks"<<std::endl;

      double Px = pythia.event[j].px();
      double Py = pythia.event[j].py();
      double Pt = sqrt(Px*Px + Py*Py);
      if(Pt < 0.2) continue;
      
      int daughter_index = j;
      int mother_index = pythia.event[daughter_index].mother1();
      int mother_id = fabs(pythia.event[mother_index].id());

      bool mother_D_meson = mother_id == 411 || mother_id == 421 || mother_id == 10411 || mother_id == 10421 ||
                            mother_id == 413 || mother_id == 423 || mother_id == 10413 || mother_id == 10423 ||
                            mother_id == 20413 || mother_id == 20423 || mother_id == 415 || mother_id == 425 ||
                            mother_id == 431 || mother_id == 10431 || mother_id == 433 || mother_id == 10433 ||
	                          mother_id == 20433 || mother_id == 435;

      if(!(mother_D_meson)) continue;

      bool parentisB = false;

      int mother_index1 = pythia.event[mother_index].mother1();
      int mother_index2 = pythia.event[mother_index].mother2();

      std::vector<int> motherIndices;
      if(mother_index1 > 0) motherIndices.push_back(mother_index1);
      if(mother_index2 > 0) motherIndices.push_back(mother_index2);

      while( !motherIndices.empty() ) {

        int motherIndex = motherIndices.front();
        motherIndices.erase(motherIndices.begin());

        if (abs(pythia.event[motherIndex].id()) == 5) { parentisB = true; break; }

        if (pythia.event[motherIndex].mother1() > 0) motherIndices.push_back(pythia.event[motherIndex].mother1());
        if (pythia.event[motherIndex].mother2() > 0) motherIndices.push_back(pythia.event[motherIndex].mother2());

      }

      if (parentisB) continue;
      
      hist->Fill(Pt);

      double Pz = pythia.event[j].pz();
      double Energy = pythia.event[j].e();
      double Vx = pythia.event[j].xProd();
      double Vy = pythia.event[j].yProd();
      double Vz = pythia.event[j].zProd();
      double Mass = -999;

      if(fabs(pdg) == 11) Mass = 0.000511;

      mass[count] = Mass; px[count] = Px; py[count] = Py; pz[count] = Pz;
      energy[count] = Energy; vx[count] = Vx; vy[count] = Vy; vz[count] = Vz;
      pid[count] = pdg; 

      myevent.mass.push_back(Mass); myevent.px.push_back(Px); myevent.py.push_back(Py); myevent.pz.push_back(Pz); myevent.energy.push_back(Energy); 
      myevent.pid.push_back(pdg);   myevent.vx.push_back(Vx); myevent.vy.push_back(Vy); myevent.vz.push_back(Vz);
      
      count++;

    }

    pTHat->Fill(pythia.info.pTHat());
    ntracks = count;

    myevent.ntracks = ntracks;
    tree->Fill();
    
    if(IsWriteOscar && ntracks>0)
    {
	    file << 0 << "\t" << ntracks << endl;
	    for(int i = 0; i < ntracks; i++){
	      if(i == -1) file << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
	      else file << i+1 << "\t" << pid[i] << "\t" << 0 << "\t" << px[i] << "\t" << py[i] << "\t" << pz[i] << "\t" << energy[i] << "\t" << mass[i] << "\t" << vx[i]*pow(10,12) << "\t" << vy[i]*pow(10,12)  << "\t" << vz[i]*pow(10,12) << "\t" << 0 << endl;
		  }
		  file << 0 << "\t" << 0 << endl;  
    }

    evt_count++;
    if(evt_count%1000==0) cout << evt_count << "\t" << "Completed" << endl;
    nevt++;
    
  }

  cout << std::setprecision(6) << pythia.info.sigmaGen() << endl;

  int nbinsX = hist->GetNbinsX();

  TH1D* pt_spectra = (TH1D*)hist->Clone();
  pt_spectra->Reset("ICESM");
  pt_spectra->SetName("pt_spectra");

  for(int ibin = 1; ibin < nbinsX + 1; ibin++){

    double bin_content = hist->GetBinContent(ibin);
    double bin_width = hist->GetBinWidth(ibin);
    double bin_center = hist->GetBinCenter(ibin);
    double bin_error = hist->GetBinError(ibin);

    pt_spectra->SetBinContent(ibin, bin_content/(2*pi*bin_width*bin_center));
    pt_spectra->SetBinError(ibin, bin_error/(2*pi*bin_width*bin_center));


  }

  TFile* fout = new TFile("tree_out.root", "RECREATE");
  fout->cd();
  pt_spectra->Write();
  tree->Write();
  pTHat->Write();
  norm->Write();
  fout->Close();

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

  std::cout << duration.count() << std::endl;

  return 0;

}
