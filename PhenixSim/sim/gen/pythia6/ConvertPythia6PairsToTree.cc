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

struct RunInfo {
  int job_id;
  long long target_pairs;
  long long generated_events;
  long long accepted_pairs;
  double sqrt_s_GeV;
  double sum_pair_BR_weights;
  double sum_relative_rep_weights;
  double br_d0_to_e;
  double br_d0_to_e_squared;
  double pythia_xsec_mb;
  double pythia_xsec_96_mb;
  double estimated_pair_xsec_mb;

  RunInfo()
    : job_id(-1),
      target_pairs(0),
      generated_events(0),
      accepted_pairs(0),
      sqrt_s_GeV(0.0),
      sum_pair_BR_weights(0.0),
      sum_relative_rep_weights(0.0),
      br_d0_to_e(0.0),
      br_d0_to_e_squared(0.0),
      pythia_xsec_mb(0.0),
      pythia_xsec_96_mb(0.0),
      estimated_pair_xsec_mb(0.0)
  {}
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

static bool isSTARGroundStateCharm(int pdg)
{
  // STAR Fig. 16 comparison: add ground-state charm hadrons
  // D0, D+, Ds, Lambda_c in |y| < 1.
  const int a = std::abs(pdg);
  return (a == 421 || a == 411 || a == 431 || a == 4122);
}

static bool isD0ForSTAR(int pdg)
{
  // STAR data conversion: D0 and anti-D0 divided by c -> D0 = 0.565.
  return std::abs(pdg) == 421;
}

static bool isDstarChargedForSTAR(int pdg)
{
  // STAR data conversion: D*+ and D*- divided by c -> D*+ = 0.224.
  // PDG: D*+ = 413, D*- = -413.
  return std::abs(pdg) == 413;
}

static TH1D* makeSTARStyleCcbarCrossSection(const TH1D* hCounts,
                                            const char* name,
                                            const char* title,
                                            const RunInfo& runInfo,
                                            double scaleExtra)
{
  TH1D* h = (TH1D*)hCounts->Clone(name);
  h->SetTitle(title);
  h->Reset();

  const double twoPi = 2.0 * std::acos(-1.0);
  const double deltaY = 2.0; // |y| < 1
  const double sigmaMB = runInfo.pythia_xsec_mb;
  const double nGen = (double)runInfo.generated_events;

  for (int b = 1; b <= hCounts->GetNbinsX(); ++b) {
    const double n = hCounts->GetBinContent(b);
    const double pt = hCounts->GetBinCenter(b);
    const double dpt = hCounts->GetBinWidth(b);

    if (nGen <= 0.0 || pt <= 0.0 || dpt <= 0.0) continue;

    // STAR-style invariant differential cross section:
    // (1 / (2*pi*pT)) * d^2 sigma / (dpT dy)
    // using generated PYTHIA MB cross section and generated-event count.
    const double denom = twoPi * pt * dpt * deltaY;
    const double val = scaleExtra * sigmaMB * n / nGen / denom;
    const double err = scaleExtra * sigmaMB * std::sqrt(n) / nGen / denom;

    h->SetBinContent(b, val);
    h->SetBinError(b, err);
  }

  return h;
}


static inline double deltaPhi0Pi(double phi1, double phi2)
{
  const double pi = std::acos(-1.0);
  double dphi = std::fabs(phi1 - phi2);
  while (dphi > 2.0*pi) dphi -= 2.0*pi;
  if (dphi > pi) dphi = 2.0*pi - dphi;
  return dphi;
}

static inline double totalMomentumFromPtEta(double pt, double eta)
{
  return pt * std::cosh(eta);
}

static double muOverElectronBR(int parentPdg)
{
  // The PYTHIA6 generator was run with charm hadrons forced to electron channels,
  // so PairInfo::weight_br is BR(Hc->e)*BR(Hcbar->e).
  // For the PHENIX dimuon-like QA proxy we convert this to a muon BR weight.
  // By default, use lepton universality: BR_mu / BR_e = 1.
  // If you want species-dependent PDG ratios, change them here.
  const int a = std::abs(parentPdg);
  if (a == 421)  return 1.0; // D0
  if (a == 411)  return 1.0; // D+
  if (a == 431)  return 1.0; // Ds
  if (a == 4122) return 1.0; // Lambda_c
  if (a == 4132) return 1.0; // Xi_c0
  if (a == 4232) return 1.0; // Xi_c+
  if (a == 4332) return 1.0; // Omega_c0
  return 1.0;
}

static TH1D* makeDSigmaDPhi1D(const TH1D* hCounts,
                              const char* name,
                              const char* title,
                              const RunInfo& runInfo)
{
  TH1D* h = (TH1D*)hCounts->Clone(name);
  h->SetTitle(title);
  h->Reset();

  const double sigmaMB = runInfo.pythia_xsec_mb;
  const double nGen = (double)runInfo.generated_events;

  for (int b = 1; b <= hCounts->GetNbinsX(); ++b) {
    const double n = hCounts->GetBinContent(b);
    const double e = hCounts->GetBinError(b);
    const double dphi = hCounts->GetBinWidth(b);
    if (nGen <= 0.0 || dphi <= 0.0) continue;

    // d sigma / d DeltaPhi [mb/rad]
    h->SetBinContent(b, sigmaMB * n / nGen / dphi);
    h->SetBinError(b, sigmaMB * e / nGen / dphi);
  }

  return h;
}

static TH2D* makeDSigmaDPhiProcess2D(const TH2D* hCounts,
                                     const char* name,
                                     const char* title,
                                     const RunInfo& runInfo)
{
  TH2D* h = (TH2D*)hCounts->Clone(name);
  h->SetTitle(title);
  h->Reset();

  const double sigmaMB = runInfo.pythia_xsec_mb;
  const double nGen = (double)runInfo.generated_events;

  for (int bx = 1; bx <= hCounts->GetNbinsX(); ++bx) {
    const double dphi = hCounts->GetXaxis()->GetBinWidth(bx);
    if (nGen <= 0.0 || dphi <= 0.0) continue;

    for (int by = 1; by <= hCounts->GetNbinsY(); ++by) {
      const double n = hCounts->GetBinContent(bx, by);
      const double e = hCounts->GetBinError(bx, by);

      // d sigma / d DeltaPhi [mb/rad] in each PYTHIA6 process bin.
      h->SetBinContent(bx, by, sigmaMB * n / nGen / dphi);
      h->SetBinError(bx, by, sigmaMB * e / nGen / dphi);
    }
  }

  return h;
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


static bool readSummaryInfo(const char* filename, RunInfo& info)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open summary file: " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;

    std::istringstream ss(line);
    std::string key;
    double val = 0.0;

    ss >> key;
    if (!ss) continue;
    if (key.size() > 0 && key[0] == '#') continue;

    ss >> val;
    if (!ss) continue;

    if (key == "job_id") {
      info.job_id = (int)val;
    }
    else if (key == "target_pairs") {
      info.target_pairs = (long long)val;
    }
    else if (key == "generated_events") {
      info.generated_events = (long long)val;
    }
    else if (key == "accepted_pairs") {
      info.accepted_pairs = (long long)val;
    }
    else if (key == "sqrt_s_GeV") {
      info.sqrt_s_GeV = val;
    }
    else if (key == "sum_pair_BR_weights") {
      info.sum_pair_BR_weights = val;
    }
    else if (key == "sum_relative_rep_weights") {
      info.sum_relative_rep_weights = val;
    }
    else if (key == "BRD0_to_e") {
      info.br_d0_to_e = val;
    }
    else if (key == "BRD0_to_e_squared") {
      info.br_d0_to_e_squared = val;
    }
    else if (key == "pythia_XSEC_0_3_mb") {
      info.pythia_xsec_mb = val;
    }
    else if (key == "pythia_XSEC_96_3_mb") {
      info.pythia_xsec_96_mb = val;
    }
    else if (key == "estimated_pair_cross_section_mb") {
      info.estimated_pair_xsec_mb = val;
    }
  }

  return true;
}

static void writeMergeableRunSummary(const RunInfo& runInfo,
                                     long long nPairsRead,
                                     long long nTreeEntries,
                                     long long nQAPairs,
                                     long long nQACharmHadrons,
                                     double sumProductionWeightRel)
{
  // These TParameters are mergeable with hadd.
  // After merging many jobs:
  //   pythia_generated_events is the total generated events.
  //   pythia_xsec_mb_sum / pythia_nfiles is the average PYTHIA xsec.
  TParameter<Long64_t>("pythia_nfiles", 1).Write();
  TParameter<Long64_t>("pythia_job_id_sum", (Long64_t)runInfo.job_id).Write();
  TParameter<Long64_t>("pythia_target_pairs", (Long64_t)runInfo.target_pairs).Write();
  TParameter<Long64_t>("pythia_generated_events", (Long64_t)runInfo.generated_events).Write();
  TParameter<Long64_t>("pythia_accepted_pairs", (Long64_t)runInfo.accepted_pairs).Write();

  TParameter<Long64_t>("n_production_pairs_read", (Long64_t)nPairsRead).Write();
  TParameter<Long64_t>("n_tree_entries", (Long64_t)nTreeEntries).Write();
  TParameter<Long64_t>("n_qa_pairs_read", (Long64_t)nQAPairs).Write();
  TParameter<Long64_t>("n_qa_charm_hadrons_read", (Long64_t)nQACharmHadrons).Write();

  TParameter<double>("pythia_sqrt_s_GeV_sum", runInfo.sqrt_s_GeV).Write();
  TParameter<double>("pythia_sum_pair_BR_weights", runInfo.sum_pair_BR_weights).Write();
  TParameter<double>("pythia_sum_relative_rep_weights", runInfo.sum_relative_rep_weights).Write();
  TParameter<double>("sum_production_weight_rel", sumProductionWeightRel).Write();

  TParameter<double>("pythia_xsec_mb_sum", runInfo.pythia_xsec_mb).Write();
  TParameter<double>("pythia_xsec_96_mb_sum", runInfo.pythia_xsec_96_mb).Write();
  TParameter<double>("pythia_estimated_pair_xsec_mb_sum", runInfo.estimated_pair_xsec_mb).Write();
  TParameter<double>("pythia_brd0_to_e_sum", runInfo.br_d0_to_e).Write();
  TParameter<double>("pythia_brd0_to_e_squared_sum", runInfo.br_d0_to_e_squared).Write();

  // Browser-friendly histograms with the same information.
  // These also merge by summing bins in hadd.
  TH1D* hRunCounts = new TH1D(
    "hRunSummaryCounts",
    "Run summary counts;quantity;sum over files",
    9, 0.5, 9.5
  );

  hRunCounts->GetXaxis()->SetBinLabel(1, "nfiles");
  hRunCounts->GetXaxis()->SetBinLabel(2, "target_pairs");
  hRunCounts->GetXaxis()->SetBinLabel(3, "generated_events");
  hRunCounts->GetXaxis()->SetBinLabel(4, "accepted_pairs");
  hRunCounts->GetXaxis()->SetBinLabel(5, "production_pairs_read");
  hRunCounts->GetXaxis()->SetBinLabel(6, "tree_entries");
  hRunCounts->GetXaxis()->SetBinLabel(7, "qa_pairs_read");
  hRunCounts->GetXaxis()->SetBinLabel(8, "qa_charm_hadrons_read");
  hRunCounts->GetXaxis()->SetBinLabel(9, "job_id_sum");

  hRunCounts->SetBinContent(1, 1.0);
  hRunCounts->SetBinContent(2, (double)runInfo.target_pairs);
  hRunCounts->SetBinContent(3, (double)runInfo.generated_events);
  hRunCounts->SetBinContent(4, (double)runInfo.accepted_pairs);
  hRunCounts->SetBinContent(5, (double)nPairsRead);
  hRunCounts->SetBinContent(6, (double)nTreeEntries);
  hRunCounts->SetBinContent(7, (double)nQAPairs);
  hRunCounts->SetBinContent(8, (double)nQACharmHadrons);
  hRunCounts->SetBinContent(9, (double)runInfo.job_id);
  hRunCounts->Write();

  TH1D* hRunWeights = new TH1D(
    "hRunSummaryWeights",
    "Run summary weights and cross sections;quantity;sum over files",
    9, 0.5, 9.5
  );

  hRunWeights->GetXaxis()->SetBinLabel(1, "sqrt_s_GeV_sum");
  hRunWeights->GetXaxis()->SetBinLabel(2, "sum_pair_BR_weights");
  hRunWeights->GetXaxis()->SetBinLabel(3, "sum_relative_rep_weights");
  hRunWeights->GetXaxis()->SetBinLabel(4, "sum_production_weight_rel");
  hRunWeights->GetXaxis()->SetBinLabel(5, "pythia_xsec_mb_sum");
  hRunWeights->GetXaxis()->SetBinLabel(6, "pythia_xsec_96_mb_sum");
  hRunWeights->GetXaxis()->SetBinLabel(7, "estimated_pair_xsec_mb_sum");
  hRunWeights->GetXaxis()->SetBinLabel(8, "BRD0_to_e_sum");
  hRunWeights->GetXaxis()->SetBinLabel(9, "BRD0_to_e_squared_sum");

  hRunWeights->SetBinContent(1, runInfo.sqrt_s_GeV);
  hRunWeights->SetBinContent(2, runInfo.sum_pair_BR_weights);
  hRunWeights->SetBinContent(3, runInfo.sum_relative_rep_weights);
  hRunWeights->SetBinContent(4, sumProductionWeightRel);
  hRunWeights->SetBinContent(5, runInfo.pythia_xsec_mb);
  hRunWeights->SetBinContent(6, runInfo.pythia_xsec_96_mb);
  hRunWeights->SetBinContent(7, runInfo.estimated_pair_xsec_mb);
  hRunWeights->SetBinContent(8, runInfo.br_d0_to_e);
  hRunWeights->SetBinContent(9, runInfo.br_d0_to_e_squared);
  hRunWeights->Write();

  delete hRunCounts;
  delete hRunWeights;
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

  RunInfo runInfo;
  readSummaryInfo(summaryFile, runInfo);

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


  // ------------------------------------------------------------------
  // PHENIX dimuon-like correlation proxy using forced dielectrons.
  // Cuts follow the charm panel of the PHENIX forward dimuon plot,
  // but applied to e+e- from open charm:
  //   1.5 < m_ee < 2.5 GeV
  //   p_e > 3 GeV/c  where p = pT*cosh(eta)
  //   1.2 < |eta_e| < 2.2
  // The electron-BR histogram uses the generator BR weights directly.
  // The muon-BR histogram rescales by BR_mu/BR_e, currently set to 1
  // in muOverElectronBR(parentPdg).
  // ------------------------------------------------------------------
  TH1D* hDielectronDphi_forwardIM_eBR_counts = new TH1D(
    "hDielectronDphi_forwardIM_eBR_counts",
    "Open-charm e^{+}e^{-}, PHENIX dimuon-like cuts, e-BR weight;#Delta#phi_{ee} [rad];BR_{e} weighted pairs",
    32, 0.0, std::acos(-1.0)
  );

  TH1D* hDielectronDphi_forwardIM_muBR_counts = new TH1D(
    "hDielectronDphi_forwardIM_muBR_counts",
    "Open-charm e^{+}e^{-}, PHENIX dimuon-like cuts, #mu-BR weight;#Delta#phi_{ee} [rad];BR_{#mu} weighted pairs",
    32, 0.0, std::acos(-1.0)
  );

  TH2D* hDielectronDphiVsProcess_forwardIM_eBR_counts = new TH2D(
    "hDielectronDphiVsProcess_forwardIM_eBR_counts",
    "Open-charm e^{+}e^{-} #Delta#phi vs PYTHIA6 process, e-BR weight;#Delta#phi_{ee} [rad];process",
    32, 0.0, std::acos(-1.0),
    8, 0.5, 8.5
  );

  TH2D* hDielectronDphiVsProcess_forwardIM_muBR_counts = new TH2D(
    "hDielectronDphiVsProcess_forwardIM_muBR_counts",
    "Open-charm e^{+}e^{-} #Delta#phi vs PYTHIA6 process, #mu-BR weight;#Delta#phi_{ee} [rad];process",
    32, 0.0, std::acos(-1.0),
    8, 0.5, 8.5
  );

  TH1D* hDielectronForwardIMInfo = new TH1D(
    "hDielectronForwardIMInfo",
    "Forward intermediate-mass dielectron QA;;value",
    6, 0.5, 6.5
  );

  hDielectronForwardIMInfo->GetXaxis()->SetBinLabel(1, "QA pairs");
  hDielectronForwardIMInfo->GetXaxis()->SetBinLabel(2, "pass fwd IM");
  hDielectronForwardIMInfo->GetXaxis()->SetBinLabel(3, "sum eBR all QA");
  hDielectronForwardIMInfo->GetXaxis()->SetBinLabel(4, "sum eBR fwd IM");
  hDielectronForwardIMInfo->GetXaxis()->SetBinLabel(5, "sum muBR fwd IM");
  hDielectronForwardIMInfo->GetXaxis()->SetBinLabel(6, "sigma fwd IM mb");
  hDielectronForwardIMInfo->LabelsOption("v", "X");
  hDielectronForwardIMInfo->SetStats(0);

  hDielectronDphi_forwardIM_eBR_counts->Sumw2();
  hDielectronDphi_forwardIM_muBR_counts->Sumw2();
  hDielectronDphiVsProcess_forwardIM_eBR_counts->Sumw2();
  hDielectronDphiVsProcess_forwardIM_muBR_counts->Sumw2();

  TH1D* hSTARCharmCounts_y1 = new TH1D(
    "hSTARCharmCounts_y1",
    "STAR-style ground-state charm hadron counts, |y|<1;p_{T} [GeV/c];counts",
    60, 0.0, 6.0
  );

  TH1D* hSTARD0Counts_y1 = new TH1D(
    "hSTARD0Counts_y1",
    "D^{0}+#bar{D}^{0} counts, |y|<1;p_{T} [GeV/c];counts",
    60, 0.0, 6.0
  );

  TH1D* hSTARDstarCounts_y1 = new TH1D(
    "hSTARDstarCounts_y1",
    "D^{*+}+D^{*-} counts, |y|<1;p_{T} [GeV/c];counts",
    60, 0.0, 6.0
  );

  labelProcessAxis(hMeeVsProcess_all);
  labelProcessAxis(hMeeVsProcess_PHENIX);
  labelProcessAxis(hMeeVsProcess_STAR);
  labelProcessAxis(hPtCharmVsProcess);
  labelProcessAxis(hDielectronDphiVsProcess_forwardIM_eBR_counts);
  labelProcessAxis(hDielectronDphiVsProcess_forwardIM_muBR_counts);

  labelParentAxis(hPtMother_all);
  labelParentAxis(hPtMother_PHENIX);
  labelParentAxis(hPtMother_STAR);
  labelParentAxis(hPtMother_BRweight_all);
  labelParentAxis(hPtMother_BRweight_PHENIX);
  labelParentAxis(hPtMother_BRweight_STAR);
  labelParentAxis(hPtCharm_all);
  labelParentAxis(hPtCharm_BRweight);

  // Fill pair QA.
  long long nForwardIMPairs = 0;
  double sumQAWeightAll = 0.0;
  double sumForwardIMElectronBRWeight = 0.0;
  double sumForwardIMMuonBRWeight = 0.0;

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

    sumQAWeightAll += p.weight_br;

    const double pAbs1 = totalMomentumFromPtEta(p.pt1, p.eta1);
    const double pAbs2 = totalMomentumFromPtEta(p.pt2, p.eta2);

    const bool passForwardIM =
      (p.pair_mass > 1.5 && p.pair_mass < 2.5 &&
       pAbs1 > 3.0 && pAbs2 > 3.0 &&
       std::fabs(p.eta1) > 1.2 && std::fabs(p.eta1) < 2.2 &&
       std::fabs(p.eta2) > 1.2 && std::fabs(p.eta2) < 2.2);

    if (passForwardIM) {
      const double dphi = deltaPhi0Pi(p.phi1, p.phi2);
      const double muWeight = p.weight_br *
        muOverElectronBR(p.parent1) * muOverElectronBR(p.parent2);

      hDielectronDphi_forwardIM_eBR_counts->Fill(dphi, p.weight_br);
      hDielectronDphi_forwardIM_muBR_counts->Fill(dphi, muWeight);
      hDielectronDphiVsProcess_forwardIM_eBR_counts->Fill(dphi, p.srcbin, p.weight_br);
      hDielectronDphiVsProcess_forwardIM_muBR_counts->Fill(dphi, p.srcbin, muWeight);

      nForwardIMPairs++;
      sumForwardIMElectronBRWeight += p.weight_br;
      sumForwardIMMuonBRWeight += muWeight;
    }

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


  TH1D* hDielectronDphi_forwardIM_eBR_dSigma_dPhi = makeDSigmaDPhi1D(
    hDielectronDphi_forwardIM_eBR_counts,
    "hDielectronDphi_forwardIM_eBR_dSigma_dPhi",
    "Open-charm e^{+}e^{-}, PHENIX dimuon-like cuts;#Delta#phi_{ee} [rad];d#sigma/d#Delta#phi [mb/rad]",
    runInfo
  );

  TH1D* hDielectronDphi_forwardIM_muBR_dSigma_dPhi = makeDSigmaDPhi1D(
    hDielectronDphi_forwardIM_muBR_counts,
    "hDielectronDphi_forwardIM_muBR_dSigma_dPhi",
    "Open-charm e^{+}e^{-} proxy for #mu^{+}#mu^{-}, PHENIX dimuon-like cuts;#Delta#phi_{ee} [rad];d#sigma/d#Delta#phi [mb/rad]",
    runInfo
  );

  TH2D* hDielectronDphiVsProcess_forwardIM_eBR_dSigma_dPhi = makeDSigmaDPhiProcess2D(
    hDielectronDphiVsProcess_forwardIM_eBR_counts,
    "hDielectronDphiVsProcess_forwardIM_eBR_dSigma_dPhi",
    "Open-charm e^{+}e^{-}, d#sigma/d#Delta#phi vs process;#Delta#phi_{ee} [rad];process",
    runInfo
  );
  labelProcessAxis(hDielectronDphiVsProcess_forwardIM_eBR_dSigma_dPhi);

  TH2D* hDielectronDphiVsProcess_forwardIM_muBR_dSigma_dPhi = makeDSigmaDPhiProcess2D(
    hDielectronDphiVsProcess_forwardIM_muBR_counts,
    "hDielectronDphiVsProcess_forwardIM_muBR_dSigma_dPhi",
    "Open-charm e^{+}e^{-} proxy for #mu^{+}#mu^{-}, d#sigma/d#Delta#phi vs process;#Delta#phi_{ee} [rad];process",
    runInfo
  );
  labelProcessAxis(hDielectronDphiVsProcess_forwardIM_muBR_dSigma_dPhi);

  const double sigmaForwardIMMuonBR =
    (runInfo.generated_events > 0) ?
    runInfo.pythia_xsec_mb * sumForwardIMMuonBRWeight / (double)runInfo.generated_events : 0.0;

  hDielectronForwardIMInfo->SetBinContent(1, (double)qaPairs.size());
  hDielectronForwardIMInfo->SetBinContent(2, (double)nForwardIMPairs);
  hDielectronForwardIMInfo->SetBinContent(3, sumQAWeightAll);
  hDielectronForwardIMInfo->SetBinContent(4, sumForwardIMElectronBRWeight);
  hDielectronForwardIMInfo->SetBinContent(5, sumForwardIMMuonBRWeight);
  hDielectronForwardIMInfo->SetBinContent(6, sigmaForwardIMMuonBR);

  // Fill charm-hadron QA independent of dielectron acceptance.
  for (size_t i = 0; i < qaCharm.size(); ++i) {
    const CharmHadronInfo& h = qaCharm[i];
    const int b = parentBin(h.pdg);
    hPtCharm_all->Fill(h.pt, b, 1.0);
    hPtCharm_BRweight->Fill(h.pt, b, h.br_e);
    hPtCharmVsProcess->Fill(h.pt, h.srcbin, 1.0);

    if (std::abs(h.y) < 1.0) {
      if (isSTARGroundStateCharm(h.pdg)) {
        hSTARCharmCounts_y1->Fill(h.pt);
      }
      if (isD0ForSTAR(h.pdg)) {
        hSTARD0Counts_y1->Fill(h.pt);
      }
      if (isDstarChargedForSTAR(h.pdg)) {
        hSTARDstarCounts_y1->Fill(h.pt);
      }
    }
  }

  TH1D* hSTAR_ccbar_d2sigma_y1 = makeSTARStyleCcbarCrossSection(
    hSTARCharmCounts_y1,
    "hSTAR_ccbar_d2sigma_y1",
    "STAR-style PYTHIA charm spectrum, |y|<1;p_{T} [GeV/c];(d^{2}#sigma)/(2#pi p_{T} dp_{T} dy) [mb/(GeV/c)^{2}]",
    runInfo,
    1.0
  );

  TH1D* hSTAR_ccbar_d2sigma_y1_half = makeSTARStyleCcbarCrossSection(
    hSTARCharmCounts_y1,
    "hSTAR_ccbar_d2sigma_y1_half",
    "Same as hSTAR_ccbar_d2sigma_y1 but divided by 2;p_{T} [GeV/c];(d^{2}#sigma)/(2#pi p_{T} dp_{T} dy) [mb/(GeV/c)^{2}]",
    runInfo,
    0.5
  );

  // STAR-data-style conversions.
  // D0 points in the STAR paper were divided by c -> D0 = 0.565.
  // D* points were divided by c -> D*+ = 0.224.
  TH1D* hSTAR_ccbar_from_D0_y1 = makeSTARStyleCcbarCrossSection(
    hSTARD0Counts_y1,
    "hSTAR_ccbar_from_D0_y1",
    "STAR-style c#bar{c} from D^{0}/0.565, |y|<1;p_{T} [GeV/c];(d^{2}#sigma)/(2#pi p_{T} dp_{T} dy) [mb/(GeV/c)^{2}]",
    runInfo,
    1.0/0.565/2
  );

  TH1D* hSTAR_ccbar_from_Dstar_y1 = makeSTARStyleCcbarCrossSection(
    hSTARDstarCounts_y1,
    "hSTAR_ccbar_from_Dstar_y1",
    "STAR-style c#bar{c} from D^{*#pm}/0.224, |y|<1;p_{T} [GeV/c];(d^{2}#sigma)/(2#pi p_{T} dp_{T} dy) [mb/(GeV/c)^{2}]",
    runInfo,
    1.0/0.224/2
  );

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
  hDielectronDphi_forwardIM_eBR_counts->Write();
  hDielectronDphi_forwardIM_muBR_counts->Write();
  hDielectronDphi_forwardIM_eBR_dSigma_dPhi->Write();
  hDielectronDphi_forwardIM_muBR_dSigma_dPhi->Write();
  hDielectronDphiVsProcess_forwardIM_eBR_counts->Write();
  hDielectronDphiVsProcess_forwardIM_muBR_counts->Write();
  hDielectronDphiVsProcess_forwardIM_eBR_dSigma_dPhi->Write();
  hDielectronDphiVsProcess_forwardIM_muBR_dSigma_dPhi->Write();
  hDielectronForwardIMInfo->Write();
  TParameter<Long64_t>("n_qa_forwardIM_pairs", (Long64_t)nForwardIMPairs).Write();
  TParameter<double>("sum_qa_forwardIM_eBR_weight", sumForwardIMElectronBRWeight).Write();
  TParameter<double>("sum_qa_forwardIM_muBR_weight", sumForwardIMMuonBRWeight).Write();
  TParameter<double>("sigma_qa_forwardIM_muBR_mb", sigmaForwardIMMuonBR).Write();
  hSTARCharmCounts_y1->Write();
  hSTARD0Counts_y1->Write();
  hSTARDstarCounts_y1->Write();
  hSTAR_ccbar_d2sigma_y1->Write();
  hSTAR_ccbar_d2sigma_y1_half->Write();
  hSTAR_ccbar_from_D0_y1->Write();
  hSTAR_ccbar_from_Dstar_y1->Write();

  writeMergeableRunSummary(runInfo,
                           nPairsRead,
                           nTreeEntries,
                           (long long)qaPairs.size(),
                           (long long)qaCharm.size(),
                           sumWeightRel);

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
