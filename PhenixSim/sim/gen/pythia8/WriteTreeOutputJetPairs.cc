#include <iostream>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

const int MaxTracks = 100;
const double pi = TMath::ACos(-1);
const double me2 = 0.000511*0.000511;

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
  TH1D* hsect = new TH1D("hsect","Sector Distribution", 9, -0.5, 8.5);

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");

  pythia.readString("HardQCD:all = on");

  pythia.readString("111:oneChannel = 1 1.0 11 22 11 -11"); // Turning ON the Dalitz decay of pion only

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + str_seed);
  pythia.readString("Next:numberCount = 100000");
  pythia.readString("PhaseSpace:pTHatMin = 1.0");
 
  //Pythia8 tune for STAR (2110.09447)
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
  TH3D* m2hist = new TH3D("m2hist","m2hist",90,0,4.5,50,0,10,3,0,3);

  std::cout<<"tree and histos were created"<<std::endl;

  TF1* ftrigger = new TF1("ftrigger","[3] + ([0]-[3])/(1 + (x/[2])^[1])", 0, 10);

  //Intialize the tree

  for(int ievt = 0; ievt < nevt; ievt++){
	//std::cout<<ievt<<std::endl;
    if(evt_count == nevents) break;

    if(!pythia.next()){ nevt++; continue; }
    int nparticles = pythia.event.size();
    
    int npart = 0;

    bool BBC_North = false; bool BBC_South = false;

    for(int ii = 0; ii < nparticles; ii++){

      if (!( pythia.event[ii].daughter1() == 0 &&  pythia.event[ii].daughter2() == 0 ) ) continue;
      if (!pythia.event[ii].isCharged()) continue;

      double eta = pythia.event[ii].eta();

      if(eta > 3.0 && eta < 3.9) BBC_North = true;
      if(eta > -3.9 && eta < -3.0) BBC_South = true;

    }

    if(!(BBC_North && BBC_South)){ nevt++;  continue; }


    int count = 0;
    myevent.set_to_null();
    std::vector<int> already_used_pi0,already_used_eta, already_used_etap, already_used_omege;

    for(int j = 0; j < nparticles; j++){

      // It should be a final state particle. It can be an electron, positron

      if (!( pythia.event[j].daughter1() == 0 &&  pythia.event[j].daughter2() == 0 ) ) continue;
      
      int pdg = pythia.event[j].id();
      if (!(pdg == 11 ||  pdg == -11) ) continue;


      int mother_index = pythia.event[j].mother1();
      int pdg_mother = pythia.event[mother_index].id();

      int isCharm = 0;
      int isBottom = 0;

      bool mother_D_meson = pdg_mother == 411 || pdg_mother == 421 || pdg_mother == 10411 || pdg_mother == 10421 ||
                            pdg_mother == 413 || pdg_mother == 423 || pdg_mother == 10413 || pdg_mother == 10423 ||
                            pdg_mother == 20413 || pdg_mother == 20423 || pdg_mother == 415 || pdg_mother == 425 ||
                            pdg_mother == 431 || pdg_mother == 10431 || pdg_mother == 433 || pdg_mother == 10433 ||
	                          pdg_mother == 20433 || pdg_mother == 435;

      bool parentisB = false;

      int mother_index1 = pythia.event[mother_index].mother1();
      int mother_index2 = pythia.event[mother_index].mother2();

      std::vector<int> motherIndices;
      if(mother_index1 > 0 && mother_index1 < mother_index) motherIndices.push_back(mother_index1);
      if(mother_index2 > 0 && mother_index2 < mother_index) motherIndices.push_back(mother_index2);
      bool skip_track = false;

      while( !motherIndices.empty() ) {

        int motherIndex = motherIndices.front();
        if(pythia.event[motherIndex].id()==111)
        {
          for (auto aready_used_mom_id : already_used_pi0)
          {
            if(motherIndex==aready_used_mom_id) skip_track = true;
          }
          already_used_pi0.push_back(motherIndex);
        }
        if(pythia.event[motherIndex].id()==221)
        {
          for (auto aready_used_mom_id : already_used_eta)
          {
            if(motherIndex==aready_used_mom_id) skip_track = true;
          }
          already_used_eta.push_back(motherIndex);
        }
        if(pythia.event[motherIndex].id()==331)
        {
          for (auto aready_used_mom_id : already_used_etap)
          {
            if(motherIndex==aready_used_mom_id) skip_track = true;
          }
          already_used_etap.push_back(motherIndex);
        }
        if(pythia.event[motherIndex].id()==3334)
        {
          for (auto aready_used_mom_id : already_used_omege)
          {
            if(motherIndex==aready_used_mom_id) skip_track = true;
          }
          already_used_omege.push_back(motherIndex);
        }
        motherIndices.erase(motherIndices.begin());

        if (abs(pythia.event[motherIndex].id()) == 5) { parentisB = true; break; }

        if (pythia.event[motherIndex].mother1() > 0 && pythia.event[motherIndex].mother1() < motherIndex) motherIndices.push_back(pythia.event[motherIndex].mother1());
        if (pythia.event[motherIndex].mother2() > 0 && pythia.event[motherIndex].mother2() < motherIndex) motherIndices.push_back(pythia.event[motherIndex].mother2());

      }

      if(mother_D_meson && !parentisB) isCharm = 1;
      if(parentisB) isBottom = 1;
      
      double Px = pythia.event[j].px();
      double Py = pythia.event[j].py();
      double eta = TMath::Abs(pythia.event[j].eta());
      double Pt = sqrt(Px*Px + Py*Py);
      if(Pt < 0.4||eta>0.5||skip_track) continue;
      //if(count>0)
      //  std::cout<<mother_index1<<" "<<mother_index2<<" "<<pythia.event[mother_index1].id()<<" "<<pythia.event[mother_index2].id()<<std::endl;
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
    if(count<2||count>9) {nevt++; continue;}
    for (int ipart = 0; ipart < count; ipart++)
    {
      for (int jpart = ipart+1; jpart < count; jpart++)
      {
        const float pair_pt = sqrt( ( px[ipart] + px[jpart] )*( px[ipart] + px[jpart] ) + ( py[ipart] + py[jpart] )*( py[ipart] + py[jpart] ) );

        const float px1 = px[ipart];
        const float py1 = py[ipart];
        const float pz1 = pz[ipart];
        const float px2 = px[jpart];
        const float py2 = py[jpart];
        const float pz2 = pz[jpart];
        const float pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
        const float pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
        const float es = sqrt(pm1 + me2) + sqrt(pm2 + me2);
        const float px = px1 + px2;
        const float py = py1 + py2;
        const float pz = pz1 + pz2;

        const float invm = sqrt(es * es - px * px - py * py - pz * pz);
        const float qbin = (pid[ipart]>0?1:0)+(pid[jpart]>0?1:0);
        m2hist->Fill(invm,pair_pt,qbin);

      }
      
    }
    
  
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

  TFile* fout = new TFile("tree_out.root", "RECREATE");
  fout->cd();
  tree->Write();
  hsect->Write();
  m2hist->Write();
  fout->Close();

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

  std::cout << duration.count() << std::endl;

  return 0;

}
