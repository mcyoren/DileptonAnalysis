#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TParameter.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"

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
  int srcbin;
  double weight_br;
  double weight_rel;
  int parent1;
  int parent2;
  double ptmom1;
  double ptmom2;
  double ymom1;
  double ymom2;
  int pdg1;
  int pdg2;
  double pt1;
  double pt2;
  double y1;
  double y2;
  double eta1;
  double eta2;
  double phi1;
  double phi2;
  double pair_mass;
  double pair_pt;
  double pair_y;
  double pair_eta;
};

struct QAPairInfo {
  PairInfo p;
  int pass_phenix;
  int pass_star;
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

struct CharmHadronInfo {
  int gen_event;
  int isub;
  int srcbin;
  int pdg;
  int species;
  double br_e;
  double pt;
  double y;
  double eta;
  double phi;
  double px;
  double py;
  double pz;
  double e;
};

static int smartRound(double w, TRandom3& rng)
{
  if (w <= 0.0) return 0;
  int n = (int)std::floor(w);
  double frac = w - (double)n;
  if (rng.Uniform() < frac) n++;
  return n;
}

static int parentBin(int pdg)
{
  int a = std::abs(pdg);
  if (a == 421)  return 1; // D0
  if (a == 411)  return 2; // D+
  if (a == 431)  return 3; // Ds
  if (a == 4122) return 4; // Lambda_c
  return 5;                // other weak charm
}

static const char* parentLabel(int b)
{
  if (b == 1) return "D0";
  if (b == 2) return "D+";
  if (b == 3) return "Ds";
  if (b == 4) return "Lambda_c";
  if (b == 5) return "other";
  return "unknown";
}

static const char* srcLabel(int b)
{
  if (b == 1) return "qqbar->ccbar";
  if (b == 2) return "gg->ccbar";
  if (b == 3) return "qg->qg";
  if (b == 4) return "qq/qqbar";
  if (b == 5) return "gg->qqbar";
  if (b == 6) return "gg->gg";
  if (b == 7) return "semihard/MPI";
  if (b == 8) return "other";
  return "unknown";
}

static void labelProcessAxis(TH2D* h)
{
  for (int b = 1; b <= 8; ++b) h->GetYaxis()->SetBinLabel(b, srcLabel(b));
}

static void labelParentAxis(TH2D* h)
{
  for (int b = 1; b <= 5; ++b) h->GetYaxis()->SetBinLabel(b, parentLabel(b));
}

static bool readProductionPairs(const char* filename,
                                std::map<int, PairInfo>& pairs)
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

    ss >> p.pair_id
       >> p.gen_event
       >> p.isub
       >> p.srcbin
       >> p.weight_br
       >> p.weight_rel
       >> p.parent1
       >> p.parent2
       >> p.ptmom1
       >> p.ptmom2
       >> p.ymom1
       >> p.ymom2
       >> p.pdg1
       >> p.pdg2
       >> p.pt1
       >> p.pt2
       >> p.y1
       >> p.y2
       >> p.eta1
       >> p.eta2
       >> p.phi1
       >> p.phi2
       >> p.pair_mass
       >> p.pair_pt
       >> p.pair_y
       >> p.pair_eta;

    if (!ss.fail()) pairs[p.pair_id] = p;
  }

  return true;
}

static bool readQAPairs(const char* filename, std::vector<QAPairInfo>& pairs)
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
    QAPairInfo q;
    PairInfo& p = q.p;

    ss >> p.pair_id
       >> p.gen_event
       >> p.isub
       >> p.srcbin
       >> p.weight_br
       >> p.weight_rel
       >> p.parent1
       >> p.parent2
       >> p.ptmom1
       >> p.ptmom2
       >> p.ymom1
       >> p.ymom2
       >> p.pdg1
       >> p.pdg2
       >> p.pt1
       >> p.pt2
       >> p.y1
       >> p.y2
       >> p.eta1
       >> p.eta2
       >> p.phi1
       >> p.phi2
       >> p.pair_mass
       >> p.pair_pt
       >> p.pair_y
       >> p.pair_eta
       >> q.pass_phenix
       >> q.pass_star;

    if (!ss.fail()) pairs.push_back(q);
  }

  return true;
}

static bool readQACharmHadrons(const char* filename,
                               std::vector<CharmHadronInfo>& hadrons)
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
    CharmHadronInfo h;

    ss >> h.gen_event
       >> h.isub
       >> h.srcbin
       >> h.pdg
       >> h.species
       >> h.br_e
       >> h.pt
       >> h.y
       >> h.eta
       >> h.phi
       >> h.px
       >> h.py
       >> h.pz
       >> h.e;

    if (!ss.fail()) hadrons.push_back(h);
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

    if (!ss.fail()) tracks[t.pair_id].push_back(t);
  }

  return true;
}

static void saveSummaryText(const char* filename)
{
  std::ifstream in(filename);
  if (!in) return;

  std::string all;
  std::string line;
  while (std::getline(in, line)) {
    all += line;
    all += "\n";
  }

  TNamed summary("pythia6_summary_text", all.c_str());
  summary.Write();
}

static void saveSummaryNumbers(const char* filename)
{
  std::ifstream in(filename);
  if (!in) return;

  std::string key;
  double val;

  while (in >> key >> val) {
    if (key.size() == 0) continue;
    if (key[0] == '#') {
      std::string dummy;
      std::getline(in, dummy);
      continue;
    }
    TParameter<double> p(key.c_str(), val);
    p.Write();
  }
}

int main(int argc, char** argv)
{
  const char* outFile = "tree_qa.root";
  if (argc > 1) outFile = argv[1];

  const char* prodPairFile = "pythia6_pairs.dat";
  const char* trackFile    = "pythia6_tracks.dat";
  const char* qaPairFile   = "pythia6_qa_pairs.dat";
  const char* qaCharmFile  = "pythia6_qa_charmhadrons.dat";
  const char* summaryFile  = "pythia6_summary.dat";

  std::map<int, PairInfo> prodPairs;
  std::map<int, std::vector<TrackInfo> > tracks;
  std::vector<QAPairInfo> qaPairs;
  std::vector<CharmHadronInfo> qaCharm;

  if (!readProductionPairs(prodPairFile, prodPairs)) return 1;
  if (!readTracks(trackFile, tracks)) return 1;
  if (!readQAPairs(qaPairFile, qaPairs)) return 1;
  if (!readQACharmHadrons(qaCharmFile, qaCharm)) return 1;

  TFile* fout = new TFile(outFile, "RECREATE");

  // ------------------------------------------------------------------
  // Final no-weight replicated tree from production PHENIX pairs.
  // ------------------------------------------------------------------
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

  // ------------------------------------------------------------------
  // QA histograms from all-pair and all-charm-hadron QA files.
  // ------------------------------------------------------------------
  TH1D* hMee_all = new TH1D("hMee_all", "all open-HF e^{+}e^{-};m_{ee} [GeV];BR-weighted pairs", 300, 0.0, 6.0);
  TH1D* hMee_PHENIX = new TH1D("hMee_PHENIX", "PHENIX-perfect open-HF e^{+}e^{-};m_{ee} [GeV];BR-weighted pairs", 300, 0.0, 6.0);
  TH1D* hMee_STAR = new TH1D("hMee_STAR", "STAR-like open-HF e^{+}e^{-};m_{ee} [GeV];BR-weighted pairs", 300, 0.0, 6.0);

  TH2D* hMeeVsProcess_all = new TH2D("hMeeVsProcess_all", "m_{ee} vs PYTHIA6 process, all;m_{ee} [GeV];process", 300, 0.0, 6.0, 8, 0.5, 8.5);
  TH2D* hMeeVsProcess_PHENIX = new TH2D("hMeeVsProcess_PHENIX", "m_{ee} vs PYTHIA6 process, PHENIX;m_{ee} [GeV];process", 300, 0.0, 6.0, 8, 0.5, 8.5);
  TH2D* hMeeVsProcess_STAR = new TH2D("hMeeVsProcess_STAR", "m_{ee} vs PYTHIA6 process, STAR;m_{ee} [GeV];process", 300, 0.0, 6.0, 8, 0.5, 8.5);

  TH2D* hPtMother_all = new TH2D("hPtMother_all", "pair mothers, all;p_{T}^{mother} [GeV];mother species", 200, 0.0, 20.0, 5, 0.5, 5.5);
  TH2D* hPtMother_PHENIX = new TH2D("hPtMother_PHENIX", "pair mothers, PHENIX;p_{T}^{mother} [GeV];mother species", 200, 0.0, 20.0, 5, 0.5, 5.5);
  TH2D* hPtMother_STAR = new TH2D("hPtMother_STAR", "pair mothers, STAR;p_{T}^{mother} [GeV];mother species", 200, 0.0, 20.0, 5, 0.5, 5.5);

  TH2D* hPtMother_BRweight_all = new TH2D("hPtMother_BRweight_all", "pair mothers BR-weighted, all;p_{T}^{mother} [GeV];mother species", 200, 0.0, 20.0, 5, 0.5, 5.5);
  TH2D* hPtMother_BRweight_PHENIX = new TH2D("hPtMother_BRweight_PHENIX", "pair mothers BR-weighted, PHENIX;p_{T}^{mother} [GeV];mother species", 200, 0.0, 20.0, 5, 0.5, 5.5);
  TH2D* hPtMother_BRweight_STAR = new TH2D("hPtMother_BRweight_STAR", "pair mothers BR-weighted, STAR;p_{T}^{mother} [GeV];mother species", 200, 0.0, 20.0, 5, 0.5, 5.5);

  TH2D* hPtCharm_all = new TH2D("hPtCharm_all", "all weak open-charm hadrons;p_{T}^{charm hadron} [GeV];species", 200, 0.0, 20.0, 5, 0.5, 5.5);
  TH2D* hPtCharm_BRweight = new TH2D("hPtCharm_BRweight", "weak open-charm hadrons weighted by BR_{e};p_{T}^{charm hadron} [GeV];species", 200, 0.0, 20.0, 5, 0.5, 5.5);
  TH2D* hPtCharmVsProcess = new TH2D("hPtCharmVsProcess", "weak open-charm hadron p_{T} vs process;p_{T}^{charm hadron} [GeV];process", 200, 0.0, 20.0, 8, 0.5, 8.5);

  labelProcessAxis(hMeeVsProcess_all);
  labelProcessAxis(hMeeVsProcess_PHENIX);
  labelProcessAxis(hMeeVsProcess_STAR);
  labelProcessAxis(hPtCharmVsProcess);

  labelParentAxis(hPtMother_all);
  labelParentAxis(hPtMother_PHENIX);
  labelParentAxis(hPtMother_STAR);
  labelParentAxis(hPtMother_BRweight_all);
  labelParentAxis(hPtMother_BRweight_PHENIX);
  labelParentAxis(hPtMother_BRweight_STAR);
  labelParentAxis(hPtCharm_all);
  labelParentAxis(hPtCharm_BRweight);

  // Fill pair QA.
  for (size_t i = 0; i < qaPairs.size(); ++i) {
    const PairInfo& p = qaPairs[i].p;
    const int b1 = parentBin(p.parent1);
    const int b2 = parentBin(p.parent2);

    hMee_all->Fill(p.pair_mass, p.weight_br);
    hMeeVsProcess_all->Fill(p.pair_mass, p.srcbin, p.weight_br);

    hPtMother_all->Fill(p.ptmom1, b1, 1.0);
    hPtMother_all->Fill(p.ptmom2, b2, 1.0);
    hPtMother_BRweight_all->Fill(p.ptmom1, b1, p.weight_br);
    hPtMother_BRweight_all->Fill(p.ptmom2, b2, p.weight_br);

    if (qaPairs[i].pass_phenix) {
      hMee_PHENIX->Fill(p.pair_mass, p.weight_br);
      hMeeVsProcess_PHENIX->Fill(p.pair_mass, p.srcbin, p.weight_br);
      hPtMother_PHENIX->Fill(p.ptmom1, b1, 1.0);
      hPtMother_PHENIX->Fill(p.ptmom2, b2, 1.0);
      hPtMother_BRweight_PHENIX->Fill(p.ptmom1, b1, p.weight_br);
      hPtMother_BRweight_PHENIX->Fill(p.ptmom2, b2, p.weight_br);
    }

    if (qaPairs[i].pass_star) {
      hMee_STAR->Fill(p.pair_mass, p.weight_br);
      hMeeVsProcess_STAR->Fill(p.pair_mass, p.srcbin, p.weight_br);
      hPtMother_STAR->Fill(p.ptmom1, b1, 1.0);
      hPtMother_STAR->Fill(p.ptmom2, b2, 1.0);
      hPtMother_BRweight_STAR->Fill(p.ptmom1, b1, p.weight_br);
      hPtMother_BRweight_STAR->Fill(p.ptmom2, b2, p.weight_br);
    }
  }

  // Fill charm-hadron QA independent of dielectron acceptance.
  for (size_t i = 0; i < qaCharm.size(); ++i) {
    const CharmHadronInfo& h = qaCharm[i];
    const int b = parentBin(h.pdg);
    hPtCharm_all->Fill(h.pt, b, 1.0);
    hPtCharm_BRweight->Fill(h.pt, b, h.br_e);
    hPtCharmVsProcess->Fill(h.pt, h.srcbin, 1.0);
  }

  // ------------------------------------------------------------------
  // Fill final tree with smart-rounded relative BR replication.
  // ------------------------------------------------------------------
  TRandom3 rng(12345);

  long long nPairsRead = 0;
  long long nTreeEntries = 0;
  double sumWeightRel = 0.0;

  for (std::map<int, PairInfo>::const_iterator it = prodPairs.begin();
       it != prodPairs.end(); ++it) {

    const int pair_id = it->first;
    const PairInfo& p = it->second;

    std::map<int, std::vector<TrackInfo> >::const_iterator jt =
      tracks.find(pair_id);

    if (jt == tracks.end()) continue;
    if (jt->second.size() < 2) continue;

    nPairsRead++;
    sumWeightRel += p.weight_rel;

    const int nrep = smartRound(p.weight_rel, rng);
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

  fout->cd();
  tree->Write();

  hMee_all->Write();
  hMee_PHENIX->Write();
  hMee_STAR->Write();
  hMeeVsProcess_all->Write();
  hMeeVsProcess_PHENIX->Write();
  hMeeVsProcess_STAR->Write();
  hPtMother_all->Write();
  hPtMother_PHENIX->Write();
  hPtMother_STAR->Write();
  hPtMother_BRweight_all->Write();
  hPtMother_BRweight_PHENIX->Write();
  hPtMother_BRweight_STAR->Write();
  hPtCharm_all->Write();
  hPtCharm_BRweight->Write();
  hPtCharmVsProcess->Write();

  TParameter<double>("n_production_pairs_read", (double)nPairsRead).Write();
  TParameter<double>("n_tree_entries", (double)nTreeEntries).Write();
  TParameter<double>("sum_production_weight_rel", sumWeightRel).Write();
  TParameter<double>("n_qa_pairs_read", (double)qaPairs.size()).Write();
  TParameter<double>("n_qa_charm_hadrons_read", (double)qaCharm.size()).Write();

  saveSummaryText(summaryFile);
  saveSummaryNumbers(summaryFile);

  fout->Close();

  std::cout << "Production pairs read: " << nPairsRead << std::endl;
  std::cout << "Tree entries:          " << nTreeEntries << std::endl;
  std::cout << "QA pairs read:         " << qaPairs.size() << std::endl;
  std::cout << "QA charm hadrons read: " << qaCharm.size() << std::endl;
  std::cout << "Wrote:                 " << outFile << std::endl;

  return 0;
}
