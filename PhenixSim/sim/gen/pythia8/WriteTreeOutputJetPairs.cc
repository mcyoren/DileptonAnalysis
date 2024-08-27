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

using namespace std;
using namespace Pythia8;

int main(){

  int nevt = 1;

  int nevents = 100000;
  Pythia pythia;

  TH1D* hsect = new TH1D("hsect","Sector Distribution", 9, -0.5, 8.5);

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");

  pythia.readString("HardQCD:all = on");

  pythia.readString("111:oneChannel = 1 1.0 11 22 11 -11"); // Turning ON the Dalitz decay of pion only

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
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

  TF1* ftrigger = new TF1("ftrigger","[3] + ([0]-[3])/(1 + (x/[2])^[1])", 0, 10);

  //Intialize the tree

  TTree* tree = new TTree("T","RECREATE");
  MyEvent event;
  tree->Branch("MyEvent", &event);

  int event_counter = 0;

  for(int ievt = 0; ievt < nevt; ievt++){

    if(event_counter == nevents) break;

    event.ClearEvent();

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

    event.SetEvtNo(ievt);

    for(int j = 0; j < nparticles; j++){

      // It should be a final state particle. It can be an electron, positron

      if (!( pythia.event[j].daughter1() == 0 &&  pythia.event[j].daughter2() == 0 ) ) continue;
      
      int pdg = pythia.event[j].id();
      if (!(pdg == 11 ||  pdg == -11) ) continue;

      double px = pythia.event[j].px();
      double py = pythia.event[j].py();
      double pz = pythia.event[j].pz();
      double vx = pythia.event[j].xProd();
      double vy = pythia.event[j].yProd();
      double vz = pythia.event[j].zProd();
      double ecore = pythia.event[j].e();
      double charge = pythia.event[j].charge();
      int mother_index = pythia.event[j].mother1();
      int pdg_mother = pythia.event[mother_index].id();

      int isCharm = 0;
      int isBottom = 0;

      double phi_EMCal = -999;
      double theta_EMCal = -999;

      TVector3 vertex(0.0, 0.0, 0.0);  // As pythia generate events at (0, 0, 0). And the purpose of this simulation is only primary or near primary particles.

      TLorentzVector particle(px, py, pz, ecore);

      int arm = InAcceptance(particle, charge, vertex, phi_EMCal, theta_EMCal);

      if( arm == 0 ) continue;

      int sector = EMCalSector(phi_EMCal);
      hsect->Fill(sector);

      //At this time, phi_EMcal and theta_EMCal will be other than -999 so the sector should not be zero. But we still put a safety condition.

      if( sector == 0 ) continue;

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

      while( !motherIndices.empty() ) {

        int motherIndex = motherIndices.front();
        motherIndices.erase(motherIndices.begin());

        if (abs(pythia.event[motherIndex].id()) == 5) { parentisB = true; break; }

        if (pythia.event[motherIndex].mother1() > 0 && pythia.event[motherIndex].mother1() < motherIndex) motherIndices.push_back(pythia.event[motherIndex].mother1());
        if (pythia.event[motherIndex].mother2() > 0 && pythia.event[motherIndex].mother2() < motherIndex) motherIndices.push_back(pythia.event[motherIndex].mother2());

      }

      if(mother_D_meson && !parentisB) isCharm = 1;
      if(parentisB) isBottom = 1;

      MyTrack track;

      // ParticleA1 includes electrons that are in the PHENIX acceptance
      // Now Calculate the Bremstruhlung Energy Loss

      double eloss = BremsEnergyLoss(0.1322)*ecore;
      double ecore_after = ecore - eloss;
      double scale = ecore_after/ecore;

      TLorentzVector particle_brem(px*scale, py*scale, pz*scale, ecore_after);

      // Now inlcude the DC resolution effects

      TLorentzVector particle_reco;
      particle_reco = ReconstructTrack(particle_brem, charge, arm);

      if(particle_reco.Pt() < 0.4) continue;

      TLorentzVector reco_shower;
      reco_shower = ReconstructShower(particle_brem, arm, sector, phi_EMCal, theta_EMCal, pdg);

      int isERT = 0;
      double rand = gRandom->Uniform(0,1);
      double energy_threshold = -999;
                                                                                                                                          
      if (sector == 1) { ftrigger->SetParameters(0.0, 7.74596, 1.80702, 1.0); energy_threshold = 1.652; } // 1.652                                                        
      else if (sector == 2) { ftrigger->SetParameters(0.0, 8.19638, 1.64472, 1.0); energy_threshold = 1.511; } // 1.511                                                   
      else if (sector == 3) { ftrigger->SetParameters(0.0, 7.56855, 1.68859, 1.0); energy_threshold = 1.541; } // 1.541                                                   
      else if (sector == 4) { ftrigger->SetParameters(0.0, 10.7913, 1.59629, 1.0); energy_threshold = 1.497; } // 1.497                                                                                
      else if (sector == 5){ ftrigger->SetParameters(0.0, 7.05411, 1.61651, 1.0); energy_threshold = 1.465; } // 1.465
      else if (sector == 6) { ftrigger->SetParameters(0.0, 8.078, 1.82598, 1.0); energy_threshold = 1.676; } // 1.676
      else if (sector == 7) { ftrigger->SetParameters(0.0, 6.87068, 2.83264, 1.0); energy_threshold = 2.561; } // 2.561
      else if (sector == 8) { ftrigger->SetParameters(0.0, 6.07937, 2.50881, 1.0); energy_threshold = 2.239; } // 2.239

      if(rand < ftrigger->Eval(reco_shower.E()) && particle_reco.E() > energy_threshold) isERT = 1;

      npart++;

      track.SetTrkID(j);
      track.SetCharge(charge);
      track.SetPID(pdg);
      track.SetMotherPID(pdg_mother);
      track.SetMotherIndex(mother_index);
      track.SetPx(particle_reco.Px());
      track.SetPy(particle_reco.Py());
      track.SetPz(particle_reco.Pz());
      track.SetVx(vx);
      track.SetVy(vy);
      track.SetVz(vz);
      track.SetPt(sqrt(particle_reco.Px()*particle_reco.Px() + particle_reco.Py()*particle_reco.Py()));
      track.SetEnergy(particle_reco.E());
      track.SetisCharm(isCharm);
      track.SetisBottom(isBottom);
      track.SetisERT(isERT);

      event.AddTrack(track);
      
    }

    nevt++;
    if(npart > 0){
      if(event_counter%1000 == 0) cout << "At Event Entry = " << event_counter << endl;
      event_counter++;
      tree->Fill();
    }

  }

  TFile* fout = new TFile("output.root", "RECREATE");
  fout->cd();
  tree->Write();
  hsect->Write();
  fout->Close();

  return 0;

}
