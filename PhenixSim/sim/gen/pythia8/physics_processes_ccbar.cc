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

using namespace Pythia8;

const int MaxTracks = 10;
const double pi = TMath::ACos(-1);

int main(){

  int nevt = 1;
  int nevents = 20000;
  Pythia pythia;

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");
  
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 1.0");

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.readString("Next:numberCount = 100000");
 
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

  //Intialize the tree

  TTree* tree = new TTree("T","RECREATE");
  int ntracks = -999;
  int pid[MaxTracks] = {-999};
  double mass[MaxTracks] = {-999};
  double px[MaxTracks] = {-999};
  double py[MaxTracks] = {-999};
  double pz[MaxTracks] = {-999};
  double energy[MaxTracks] = {-999};
  double vx[MaxTracks] = {-999};
  double vy[MaxTracks] = {-999};
  double vz[MaxTracks] = {-999};

  tree->Branch("ntracks", &ntracks,"ntracks/I");
  tree->Branch("pid",pid,"pid[ntracks]/I");
  tree->Branch("mass",mass,"mass[ntracks]/D");
  tree->Branch("px",px,"px[ntracks]/D");
  tree->Branch("py",py,"py[ntracks]/D");
  tree->Branch("pz",pz,"pz[ntracks]/D");
  tree->Branch("energy",energy,"energy[ntracks]/D");
  tree->Branch("vx",vx,"vx[ntracks]/D");
  tree->Branch("vy",vy,"vy[ntracks]/D");
  tree->Branch("vz",vz,"vz[ntracks]/D");

  int evt_count = 0;

  TH1D* hist = new TH1D("hist", "pT Distribution", 100, 0, 10);
  TH1D* pTHat = new TH1D("pTHat", "pTHat Distribution", 100, 0, 10);
  TH1D* norm = new TH1D("norm","EVENT COUNTER", 2,-0.5,1.5);

  for(int ievt = 0; ievt < nevt; ievt++){

    if(evt_count == nevents) break;

    if(!pythia.next()){ nevt++; continue; }
    int nparticles = pythia.event.size();
    int npart = 0;

    norm->Fill(0);

    for(int j = 0; j < nparticles; j++){

      if (!( pythia.event[j].daughter1() == 0 &&  pythia.event[j].daughter2() == 0 && fabs(pythia.event[j].eta()) < 0.5 ) ) continue;

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

    if(npart<1){ nevt++; continue; }

    int count = 0;

    for(int j = 0; j < nparticles; j++){

      if (!( pythia.event[j].daughter1() == 0 &&  pythia.event[j].daughter2() == 0 && fabs(pythia.event[j].eta()) < 0.5 ) ) continue;

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

      count++;

    }

    pTHat->Fill(pythia.info.pTHat());
    ntracks = count;

    tree->Fill();
    evt_count++;
    if(evt_count%10==0) cout << evt_count << "\t" << "Completed" << endl;
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

  return 0;

}
