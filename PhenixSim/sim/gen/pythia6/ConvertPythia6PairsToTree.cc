#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TParameter.h"

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

  void clear() {
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
  }
};

struct PairInfo {
  int pair_id;
  int gen_event;
  int isub;
  double weight_br;
  double weight_rel;
  int parent1;
  int parent2;
  int pdg1;
  int pdg2;
};

struct TrackInfo {
  int pair_id;
  int itrack;
  int pid;
  int parent;
  double mass;
  double energy;
  double px;
  double py;
  double pz;
  double vx;
  double vy;
  double vz;
};

static int smartRound(double w, TRandom3& rng)
{
  if (w <= 0.0) return 0;
  int n = (int)std::floor(w);
  double frac = w - (double)n;
  if (rng.Uniform() < frac) n++;
  return n;
}

static bool readPairs(const char* filename, std::map<int, PairInfo>& pairs)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);

    PairInfo p;

    double dummy;
    ss >> p.pair_id
       >> p.gen_event
       >> p.isub
       >> p.weight_br
       >> p.weight_rel
       >> p.parent1
       >> p.parent2
       >> p.pdg1
       >> p.pdg2;

    // Skip the remaining kinematic columns.
    while (ss >> dummy) {}

    if (!ss.bad()) {
      pairs[p.pair_id] = p;
    }
  }

  return true;
}

static bool readTracks(const char* filename,
                       std::map<int, std::vector<TrackInfo> >& tracks)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);

    TrackInfo t;
    ss >> t.pair_id
       >> t.itrack
       >> t.pid
       >> t.parent
       >> t.mass
       >> t.energy
       >> t.px
       >> t.py
       >> t.pz
       >> t.vx
       >> t.vy
       >> t.vz;

    if (!ss.fail()) {
      tracks[t.pair_id].push_back(t);
    }
  }

  return true;
}

int main(int argc, char** argv)
{
  const char* pairFile = "pythia6_pairs.dat";
  const char* trackFile = "pythia6_tracks.dat";
  const char* outFile = "tree_out.root";

  if (argc > 1) outFile = argv[1];

  std::map<int, PairInfo> pairs;
  std::map<int, std::vector<TrackInfo> > tracks;

  if (!readPairs(pairFile, pairs)) return 1;
  if (!readTracks(trackFile, tracks)) return 1;

  TFile* fout = new TFile(outFile, "RECREATE");

  TTree* tree = new TTree("T", "PYTHIA6 forced-e open-HF dielectron tree");

  MyEvent ev;

  tree->Branch("ntracks", &ev.ntracks, "ntracks/I");
  tree->Branch("pid",    &ev.pid);
  tree->Branch("mass",   &ev.mass);
  tree->Branch("energy", &ev.energy);
  tree->Branch("px",     &ev.px);
  tree->Branch("py",     &ev.py);
  tree->Branch("pz",     &ev.pz);
  tree->Branch("vx",     &ev.vx);
  tree->Branch("vy",     &ev.vy);
  tree->Branch("vz",     &ev.vz);

  TRandom3 rng(12345);

  long long nPairsRead = 0;
  long long nTreeEntries = 0;
  double sumWeightRel = 0.0;

  for (std::map<int, PairInfo>::const_iterator it = pairs.begin();
       it != pairs.end(); ++it) {

    const int pair_id = it->first;
    const PairInfo& p = it->second;

    std::map<int, std::vector<TrackInfo> >::const_iterator jt =
      tracks.find(pair_id);

    if (jt == tracks.end()) continue;
    if (jt->second.size() < 2) continue;

    nPairsRead++;
    sumWeightRel += p.weight_rel;

    int nrep = smartRound(p.weight_rel, rng);
    if (nrep <= 0) continue;

    for (int ir = 0; ir < nrep; ++ir) {
      ev.clear();

      const std::vector<TrackInfo>& tv = jt->second;

      for (size_t k = 0; k < tv.size(); ++k) {
        ev.pid.push_back(tv[k].pid);
        ev.mass.push_back(tv[k].mass);
        ev.energy.push_back(tv[k].energy);
        ev.px.push_back(tv[k].px);
        ev.py.push_back(tv[k].py);
        ev.pz.push_back(tv[k].pz);
        ev.vx.push_back(tv[k].vx);
        ev.vy.push_back(tv[k].vy);
        ev.vz.push_back(tv[k].vz);
      }

      ev.ntracks = (int)ev.pid.size();

      tree->Fill();
      nTreeEntries++;
    }
  }

  TParameter<long long>* pPairsRead =
    new TParameter<long long>("n_pairs_read", nPairsRead);

  TParameter<long long>* pTreeEntries =
    new TParameter<long long>("n_tree_entries", nTreeEntries);

  TParameter<double>* pSumWeightRel =
    new TParameter<double>("sum_weight_rel", sumWeightRel);

  fout->cd();
  tree->Write();
  pPairsRead->Write();
  pTreeEntries->Write();
  pSumWeightRel->Write();
  fout->Close();

  std::cout << "Read pairs:      " << nPairsRead << std::endl;
  std::cout << "Tree entries:    " << nTreeEntries << std::endl;
  std::cout << "Sum weight_rel:  " << sumWeightRel << std::endl;
  std::cout << "Wrote:           " << outFile << std::endl;

  return 0;
}