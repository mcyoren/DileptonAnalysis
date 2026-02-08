// WriteTreeoutputccbarSoft.cc
//
// SoftQCD + charm-only (no bottom) via UserHook
// Force charm hadrons -> e channels only
// For each event: pick one OS ee pair in gate (|y|<0.5, pT>0.2)
// Replicate event N times where N ~ BR(H1)*BR(H2)/BR(D0)^2 using unbiased smart rounding
// Write TTree ONLY when gate has >=2 electrons and replication count > 0
//
// Extra QA requested:
// - histogram of all BRs used (per PDG, bin labels)
// - mee_int: integer round-up scheme: 1 for D0D0 else ceil(wpair)
// - mee_wpair: weighted by wpair = BRprod/BRD0^2
// - mee_brprod: weighted by BRprod = BR(H1)*BR(H2)
//
// Build:
//   g++ -O2 -std=c++17 WriteTreeoutputccbarSoft.cc -o WriteTreeoutputccbarSoft \
//     $(root-config --cflags --libs) $(pythia8-config --cxxflags --libs)
//
// Run:
//   ./WriteTreeoutputccbarSoft <seed> <tree_entries_to_write>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <chrono>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TParameter.h>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// -----------------------------
// TTree container (your vector format)
// -----------------------------
struct MyEvent {
  int ntracks = 0;
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
    pid.clear(); mass.clear(); energy.clear();
    px.clear(); py.clear(); pz.clear();
    vx.clear(); vy.clear(); vz.clear();
  }
};

// -----------------------------
// UserHook: keep events with charm, veto any bottom
// -----------------------------
class CharmOnlyNoBottomHook : public UserHooks {
public:
  long long nCalls = 0;
  long long nVeto  = 0;

  bool canVetoPartonLevel() override { return true; }

  bool doVetoPartonLevel(const Event& event) override {
    ++nCalls;

    bool hasCharm  = false;
    bool hasBottom = false;

    for (int i = 0; i < event.size(); ++i) {
      int id = std::abs(event[i].id());
      if (id == 4) hasCharm = true;
      if (id == 5) hasBottom = true;
      if (hasCharm && hasBottom) break;
    }

    const bool veto = (!hasCharm) || hasBottom;
    if (veto) ++nVeto;
    return veto;
  }
};


// -----------------------------
// Helpers
// -----------------------------
static inline bool isElectron(int pdg) { return (pdg == 11 || pdg == -11); }
static inline int chargeFromPdg(int pdg) { return (pdg > 0) ? -1 : +1; } // e-:11 -> -1 ; e+:-11 -> +1

struct CandE {
  int idx = -1;      // index in pythia.event
  int pdg = 0;
  int q = 0;
  int parent = 0;    // charm parent PDG id (abs)
  TLorentzVector p4;
  double pt = 0;
  double y  = 0;
  double eta = 0;
};

static inline bool passGateY05Pt02(const CandE& e) {
  return (e.pt > 0.2 && std::fabs(e.y) < 0.5);
}

// climb mothers until you hit an allowed charm parent, or stop
static int findCharmParent(const Event& ev, int idx,
                           const std::set<int>& charmHadronSet,
                           int maxSteps = 60)
{
  int cur = idx;
  for (int step = 0; step < maxSteps; ++step) {
    if (cur <= 0 || cur >= ev.size()) return 0;
    const int m1 = ev[cur].mother1();
    const int m2 = ev[cur].mother2();
    const int mom = (m1 > 0 ? m1 : (m2 > 0 ? m2 : 0));
    if (mom <= 0 || mom >= ev.size()) return 0;

    const int idAbs = std::abs(ev[mom].id());
    if (charmHadronSet.count(idAbs)) return idAbs;

    cur = mom;
  }
  return 0;
}

// choose ONE OS pair: pick highest-pt e then best opposite-sign partner by pt
static bool pickBestOSPair(const std::vector<CandE>& v, int& iBest, int& jBest)
{
  iBest = -1; jBest = -1;
  if (v.size() < 2) return false;

  std::vector<int> idx(v.size());
  for (size_t i = 0; i < v.size(); ++i) idx[i] = (int)i;

  std::sort(idx.begin(), idx.end(),
            [&](int a, int b){ return v[a].pt > v[b].pt; });

  for (size_t ia = 0; ia < idx.size(); ++ia) {
    const int i = idx[ia];
    int bestj = -1;
    double bestpt = -1;

    for (size_t jb = ia + 1; jb < idx.size(); ++jb) {
      const int j = idx[jb];
      if (v[i].q * v[j].q >= 0) continue;
      if (v[j].pt > bestpt) { bestpt = v[j].pt; bestj = j; }
    }

    if (bestj >= 0) { iBest = i; jBest = bestj; return true; }
  }
  return false;
}

// unbiased "smart rounding": floor(w) + Bernoulli(frac)
static int smartRound(double w, TRandom3& rng)
{
  if (w <= 0) return 0;
  const int n = (int)std::floor(w);
  const double frac = w - n;
  int out = n;
  if (rng.Uniform() < frac) out += 1;
  return out;
}

static inline double wrapPhi(double a) {
  const double PI = TMath::Pi();
  const double TWOPI = 2.0*PI;
  while (a <= -PI) a += TWOPI;
  while (a >   PI) a -= TWOPI;
  return a;
}

static inline bool inPhenixArm(double phi) {
  const double PI = TMath::Pi();
  const double phi_west_low = -3*PI/16.0, phi_west_up = 5*PI/16.0;
  const double phi_east_low = 11*PI/16.0, phi_east_up = 19*PI/16.0;

  // your code shifts to [-pi/2, 3pi/2)
  if (phi < -PI/2.0) phi += 2.0*PI;

  const bool west = (phi > phi_west_low && phi < phi_west_up);
  const bool east = (phi > phi_east_low && phi < phi_east_up);
  return (west || east);
}

static inline bool passPhenixPhiAcc(double px, double py, int q, double pt) {
  if (pt <= 0) return false;

  const double PI = TMath::Pi();
  const double k_DC   = 0.206; // rad GeV/c
  const double k_RICH = 0.309; // rad GeV/c

  double phi = std::atan2(py, px);

  double phi_rich = phi + q * k_RICH / pt;
  double phi_dc   = phi + q * k_DC   / pt;

  if (phi_rich < -PI/2.0) phi_rich += 2.0*PI;
  if (phi_dc   < -PI/2.0) phi_dc   += 2.0*PI;

  return (inPhenixArm(phi_rich) && inPhenixArm(phi_dc));
}


int main(int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <seed> <tree_entries_to_write>\n";
    return 1;
  }

  const int seed = std::stoi(argv[1]);
  const int targetTreeEntries = std::stoi(argv[2]);

  std::cout << "seed=" << seed << " targetTreeEntries=" << targetTreeEntries << "\n";
  auto t0 = std::chrono::high_resolution_clock::now();

  TRandom3 rng(seed + 12345);

  // ----------------------------
  // Pythia setup
  // ----------------------------
  Pythia pythia;

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");
  pythia.readString("SoftQCD:inelastic = on");

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 100000000");

  // STAR-like tune bits you had
  //pythia.readString("PDF:pSet = 17");
  //pythia.readString("MultipartonInteractions:ecmRef = 200");
  //pythia.readString("MultipartonInteractions:bprofile = 2");
  //pythia.readString("MultipartonInteractions:pT0Ref = 1.40");
  //pythia.readString("MultipartonInteractions:ecmPow = 0.135");
  //pythia.readString("MultipartonInteractions:coreRadius = 0.56");
  //pythia.readString("MultipartonInteractions:coreFraction = 0.78");
  //pythia.readString("ColourReconnection:range = 5.4");

  std::shared_ptr<CharmOnlyNoBottomHook> hook = std::make_shared<CharmOnlyNoBottomHook>();
  pythia.setUserHooksPtr(hook);

  // ----------------------------
  // Charm parents to force into e± channels
  // Keep 421 (D0) as reference!
  // ----------------------------
  const std::vector<int> charmParents = {
    421, 411, 431, 4122,            // D0, D+, Ds, Lambda_c
    10421,10411, 423,413,           // excited D
    10423,10413, 20423,20413,
    425,415,
    10431,433,10433,20433,435
  };

  std::set<int> charmParentSet;
  for (int id : charmParents) charmParentSet.insert(std::abs(id));

  for (int pdgID : charmParents) {
    const int idAbs = std::abs(pdgID);
    pythia.readString(std::to_string(idAbs) + ":onMode = off");
    pythia.readString(std::to_string(idAbs) + ":onIfAny = 11 -11");
  }

  pythia.init();

  // Compute BR(H->e±+X) for each forced parent
  std::map<int,double> brE;
  for (int pdgID : charmParents) {
    const int idAbs = std::abs(pdgID);
    double totalBR = 0.0;

    auto entryPtr = pythia.particleData.particleDataEntryPtr(idAbs);
    if (entryPtr) {
      for (int i = 0; i < entryPtr->sizeChannels(); ++i) {
        if (entryPtr->channel(i).onMode()) totalBR += entryPtr->channel(i).bRatio();
      }
    }
    brE[idAbs] = totalBR;
  }

  const double brD0 = (brE.count(421) ? brE[421] : 0.0);
  if (brD0 <= 0) {
    std::cerr << "ERROR: BR(D0->e) computed as " << brD0 << " (<=0)\n";
    return 2;
  }
  const double brD0sq = brD0 * brD0;

  std::cout << "BR(D0->e) = " << std::setprecision(10) << brD0
            << "  BR(D0)^2 = " << brD0sq << "\n";

  // ----------------------------
  // Output ROOT
  // ----------------------------
  TFile* fout = new TFile("tree_out.root", "RECREATE");

  // Tree
  TTree* tree = new TTree("T","RECREATE");
  MyEvent myevent;

  tree->Branch("ntracks", &myevent.ntracks, "ntracks/I");
  tree->Branch("pid",    &myevent.pid);
  tree->Branch("mass",   &myevent.mass);
  tree->Branch("energy", &myevent.energy);
  tree->Branch("px",     &myevent.px);
  tree->Branch("py",     &myevent.py);
  tree->Branch("pz",     &myevent.pz);
  tree->Branch("vx",     &myevent.vx);
  tree->Branch("vy",     &myevent.vy);
  tree->Branch("vz",     &myevent.vz);

  // BR histogram (bin labels = PDG)
  TH1D* hBR = new TH1D("br_used", "BR(H#rightarrow e^{#pm}+X) used;parent PDG;BR",
                       (int)charmParents.size(), 0.5, (double)charmParents.size() + 0.5);
  for (size_t i = 0; i < charmParents.size(); ++i) {
    const int idAbs = std::abs(charmParents[i]);
    hBR->GetXaxis()->SetBinLabel((int)i+1, std::to_string(idAbs).c_str());
    hBR->SetBinContent((int)i+1, brE[idAbs]);
  }

  // Mass spectra QA
  TH1D* hMee_int    = new TH1D("mee_int",
    "m_{ee} with integer round-up: 1 for D0D0 else ceil(BRprod/BRD0^{2});m_{ee} (GeV/c^{2});counts", 300, 0, 6);

  TH1D* hMee_wpair  = new TH1D("mee_wpair",
    "m_{ee} weighted by wpair=BRprod/BRD0^{2};m_{ee};weighted counts", 300, 0, 6);

  TH1D* hMee_brprod = new TH1D("mee_brprod",
    "m_{ee} weighted by BRprod=BR(H1)BR(H2);m_{ee};weighted counts", 300, 0, 6);

  // This is what your replicated tree represents (smart rounding)
  TH1D* hMee_tree   = new TH1D("mee_tree",
    "m_{ee} from replicated TTree entries (smart rounding);m_{ee};counts", 300, 0, 6);

  TH1D* hWpair = new TH1D("w_pair",
    "w_{pair}=BRprod/BRD0^{2};w_{pair};counts", 200, 0, 5);

  // Counters (QA bookkeeping)
  TH1D* hCounts = new TH1D("counts", "Counters;bin;value", 6, -0.5, 5.5);
  // 0 tried (pythia.next success)
  // 1 events with >=2 gated electrons
  // 2 events with OS pair found
  // 3 sum_nrep_smart (accumulated)
  // 4 tree entries written
  // 5 reserved

  // Parent ID frequency for selected OS pair (2 fills per accepted pair)
  // --- parent category map: PDG -> bin index (1..N), plus OTHER bin
  std::vector<int> parentCats = charmParents;         // same list you force
  std::sort(parentCats.begin(), parentCats.end());
  parentCats.erase(std::unique(parentCats.begin(), parentCats.end()), parentCats.end());

  std::map<int,int> parentToBin;
  for (size_t i = 0; i < parentCats.size(); ++i) parentToBin[parentCats[i]] = (int)i + 1;

  const int BIN_OTHER = (int)parentCats.size() + 1;

  TH1D* hParentPair = new TH1D("parent_pdg_pair",
  "Charm parent (selected ee pair); category; counts",
  BIN_OTHER, 0.5, BIN_OTHER + 0.5);

  TH1D* hParentAll = new TH1D("parent_pdg_all",
    "Charm parent (all gated electrons); category; counts",
    BIN_OTHER, 0.5, BIN_OTHER + 0.5);

  // set bin labels
  for (size_t i = 0; i < parentCats.size(); ++i) {
    hParentPair->GetXaxis()->SetBinLabel((int)i+1, std::to_string(parentCats[i]).c_str());
    hParentAll ->GetXaxis()->SetBinLabel((int)i+1, std::to_string(parentCats[i]).c_str());
  }
  hParentPair->GetXaxis()->SetBinLabel(BIN_OTHER, "OTHER");
  hParentAll ->GetXaxis()->SetBinLabel(BIN_OTHER, "OTHER");

  // nicer look
  hParentPair->LabelsOption("v","X");
  hParentAll ->LabelsOption("v","X");

  auto fillParent = [&](TH1D* h, int pdgAbs) {
  auto it = parentToBin.find(pdgAbs);
  if (it != parentToBin.end()) h->Fill(it->second);
  else h->Fill(BIN_OTHER);
  };

  TH1D* hPtAllCharmE = new TH1D("pt_all",
    "e^{#pm} p_{T} (all charm e);p_{T} (GeV/c);counts", 100, 0, 10);

  TH1D* hMee_exact2        = new TH1D("mee_exact2",
    "m_{ee} (OS), exactly 2 charm e (no acceptance);m_{ee} (GeV/c^{2});counts", 300, 0, 6);

  TH1D* hMee_star          = new TH1D("mee_star",
    "m_{ee} (OS), STAR |y|<1 p_{T}>0.2;m_{ee};counts", 300, 0, 6);

  TH1D* hMee_phenix_eta05  = new TH1D("mee_phenix_eta05",
    "m_{ee} (OS), PHENIX sim |#eta|<0.5 p_{T}>0.2;m_{ee};counts", 300, 0, 6);

  TH1D* hMee_phenix_y035   = new TH1D("mee_phenix_y035",
    "m_{ee} (OS), PHENIX |y|<0.35 p_{T}>0.2;m_{ee};counts", 300, 0, 6);

  TH1D* hMee_phenix_phiacc = new TH1D("mee_phenix_phiacc",
    "m_{ee} (OS), PHENIX |y|<0.35 p_{T}>0.2 + DC/RICH #phi cuts;m_{ee};counts", 300, 0, 6);

  TH1D* hPt_star        = new TH1D("pt_star",
    "e^{#pm} p_{T} (STAR |y|<1);p_{T};counts", 100, 0, 10);

  TH1D* hPt_phenix_eta  = new TH1D("pt_phenix_eta",
    "e^{#pm} p_{T} (PHENIX sim |#eta|<0.5);p_{T};counts", 100, 0, 10);

  TH1D* hPt_phenix_y    = new TH1D("pt_phenix_y",
    "e^{#pm} p_{T} (PHENIX |y|<0.35);p_{T};counts", 100, 0, 10);

  TH1D* hPt_phenix_phi  = new TH1D("pt_phenix_phi",
    "e^{#pm} p_{T} (PHENIX |y|<0.35 + DC/RICH #phi cuts);p_{T};counts", 100, 0, 10);

  TH1D* hPTHat = new TH1D("pTHat",
    "pTHat Distribution;#hat{p}_{T} (GeV/c);counts", 100, 0, 10);



  // store constants too
  TParameter<double>* pBRD0   = new TParameter<double>("BRD0", brD0);
  TParameter<double>* pBRD0sq = new TParameter<double>("BRD0sq", brD0sq);

  // ----------------------------
  // Generation loop: keep going until tree entries == targetTreeEntries
  // ----------------------------
  long long tried = 0;
  long long overall = 0;
  long long writtenEntries = 0;
  long long sumNrepSmart = 0;

  // also track how well the “int round-up” compares:
  long long sumIntWeight = 0;
  double sumWpair = 0.0;
  double sumBRprod = 0.0;

  // Main event loop: keep generating until the target tree size is reached
  while (writtenEntries < targetTreeEntries) {

    overall++;
    if (!pythia.next()) continue;
    ++tried;
    hCounts->Fill(0);

    const int nparticles = pythia.event.size();

    hPTHat->Fill(pythia.info.pTHat());

    std::vector<CandE> allCharm;  allCharm.reserve(16);
    std::vector<CandE> star;      star.reserve(16);
    std::vector<CandE> phen_eta;  phen_eta.reserve(16);
    std::vector<CandE> phen_y;    phen_y.reserve(16);
    std::vector<CandE> phen_phi;  phen_phi.reserve(16);

    // collect gated electrons from allowed charm parents
    std::vector<CandE> gated; // this is the TREE/BR gate: |y|<0.5, pT>0.2
    gated.reserve(8);

    // Loop over all final-state particles and collect charm-origin electrons
    for (int j = 0; j < nparticles; ++j) {
      if (!pythia.event[j].isFinal()) continue;
    
      const int pdg = pythia.event[j].id();
      if (!isElectron(pdg)) continue;
    
      const double Px = pythia.event[j].px();
      const double Py = pythia.event[j].py();
      const double Pz = pythia.event[j].pz();
      const double E  = pythia.event[j].e();
    
      TLorentzVector p4(Px, Py, Pz, E);
      const double pt = p4.Pt();
    
      CandE c;
      c.idx = j;
      c.pdg = pdg;
      c.q   = chargeFromPdg(pdg);
      c.p4  = p4;
      c.pt  = pt;
      c.y   = p4.Rapidity();
      c.eta = p4.Eta();
    
      // Parent finding: climb the history to the nearest allowed charm hadron
      const int parent = findCharmParent(pythia.event, j, charmParentSet);
      if (parent == 0) continue;
    
      c.parent = parent;
    
      // ---- "all charm electrons" QA (no acceptance cuts yet)
      allCharm.push_back(c);
      hPtAllCharmE->Fill(pt);
      if (pt < 0.2) continue;
      fillParent(hParentAll, std::abs(parent));
    
      // ---- STAR-like (|y|<1, pT>0.2)
      if (std::fabs(c.y) < 1.0) {
        star.push_back(c);
        hPt_star->Fill(pt);
      }
    
      // ---- PHENIX sim acceptance (|eta|<0.5)
      if (std::fabs(c.eta) < 0.5) {
        phen_eta.push_back(c);
        hPt_phenix_eta->Fill(pt);
      }
    
      // ---- PHENIX physics acceptance (|y|<0.35)
      if (std::fabs(c.y) < 0.35) {
        phen_y.push_back(c);
        hPt_phenix_y->Fill(pt);
      
        // ---- PHENIX phi acceptance (DC+RICH bending)
        if (passPhenixPhiAcc(Px, Py, c.q, pt)) {
          phen_phi.push_back(c);
          hPt_phenix_phi->Fill(pt);
        }
      }
    
      // ---- TREE/BR gate for weighted replication (|y|<0.5, pT>0.2)
      if (std::fabs(c.y) < 0.5) {
        gated.push_back(c);
      }
    }

    // exactly 2 charm electrons (no further acc), OS pair
    if ((int)allCharm.size() == 2) {
      int ia=-1, ib=-1;
      if (pickBestOSPair(allCharm, ia, ib)) {
        const CandE& e1 = allCharm[ia];
        const CandE& e2 = allCharm[ib];
        const double mee = (e1.p4 + e2.p4).M();

        const int H1 = e1.parent;
        const int H2 = e2.parent;

        const double br1 = brE.count(H1) ? brE[H1] : 0.0;
        const double br2 = brE.count(H2) ? brE[H2] : 0.0;
        if (br1 <= 0 || br2 <= 0) continue;

        const double brprod = br1 * br2;
        const double wpair  = brprod / brD0sq;
        hMee_exact2->Fill(mee, wpair);
      }
    }

    // STAR
    {
      int ia=-1, ib=-1;
      if (pickBestOSPair(star, ia, ib)) {
        const CandE& e1 = star[ia];
        const CandE& e2 = star[ib];
        const double mee = (e1.p4 + e2.p4).M();

        const int H1 = e1.parent;
        const int H2 = e2.parent;

        const double br1 = brE.count(H1) ? brE[H1] : 0.0;
        const double br2 = brE.count(H2) ? brE[H2] : 0.0;
        if (br1 <= 0 || br2 <= 0) continue;

        const double brprod = br1 * br2;
        const double wpair  = brprod / brD0sq;
        hMee_star->Fill(mee, wpair);
      }
    }

    // PHENIX sim eta
    {
      int ia=-1, ib=-1;
      if (pickBestOSPair(phen_eta, ia, ib)) {
        const CandE& e1 = phen_eta[ia];
        const CandE& e2 = phen_eta[ib];
        const double mee = (e1.p4 + e2.p4).M();

        const int H1 = e1.parent;
        const int H2 = e2.parent;

        const double br1 = brE.count(H1) ? brE[H1] : 0.0;
        const double br2 = brE.count(H2) ? brE[H2] : 0.0;
        if (br1 <= 0 || br2 <= 0) continue;

        const double brprod = br1 * br2;
        const double wpair  = brprod / brD0sq;
        hMee_phenix_eta05->Fill(mee, wpair);
      }
    }

    // PHENIX y
    {
      int ia=-1, ib=-1;
      if (pickBestOSPair(phen_y, ia, ib)) {
        const CandE& e1 = phen_y[ia];
        const CandE& e2 = phen_y[ib];
        const double mee = (e1.p4 + e2.p4).M();

        const int H1 = e1.parent;
        const int H2 = e2.parent;

        const double br1 = brE.count(H1) ? brE[H1] : 0.0;
        const double br2 = brE.count(H2) ? brE[H2] : 0.0;
        if (br1 <= 0 || br2 <= 0) continue;

        const double brprod = br1 * br2;
        const double wpair  = brprod / brD0sq;
        hMee_phenix_y035->Fill(mee, wpair);
      }
    }

    // PHENIX phi-acc
    {
      int ia=-1, ib=-1;
      if (pickBestOSPair(phen_phi, ia, ib)) {
        const CandE& e1 = phen_phi[ia];
        const CandE& e2 = phen_phi[ib];
        const double mee = (e1.p4 + e2.p4).M();

        const int H1 = e1.parent;
        const int H2 = e2.parent;

        const double br1 = brE.count(H1) ? brE[H1] : 0.0;
        const double br2 = brE.count(H2) ? brE[H2] : 0.0;
        if (br1 <= 0 || br2 <= 0) continue;

        const double brprod = br1 * br2;
        const double wpair  = brprod / brD0sq;
        hMee_phenix_phiacc->Fill(mee, wpair);
      }
    }


    if ((int)gated.size() < 2) continue;
    hCounts->Fill(1);

    // pick best OS pair
    int ia=-1, ib=-1;
    if (!pickBestOSPair(gated, ia, ib)) continue;
    hCounts->Fill(2);

    const CandE& e1 = gated[ia];
    const CandE& e2 = gated[ib];
    const double mee = (e1.p4 + e2.p4).M();

    const int H1 = e1.parent;
    const int H2 = e2.parent;

    const double br1 = brE.count(H1) ? brE[H1] : 0.0;
    const double br2 = brE.count(H2) ? brE[H2] : 0.0;
    if (br1 <= 0 || br2 <= 0) continue;

    fillParent(hParentPair, H1);
    fillParent(hParentPair, H2);

    const double brprod = br1 * br2;
    const double wpair  = brprod / brD0sq;

    sumWpair  += wpair;
    sumBRprod += brprod;

    hWpair->Fill(wpair);
    hMee_wpair->Fill(mee, wpair);
    hMee_brprod->Fill(mee, brprod);

    // ---- integer round-up scheme for QA (as you requested)
    int wInt = 0;
    if (H1 == 421 && H2 == 421) wInt = 1;               // D0D0 exactly -> 1
    else                        wInt = (int)std::ceil(wpair); // other -> ceil
    if (wInt < 0) wInt = 0;

    sumIntWeight += wInt;
    for (int k = 0; k < wInt; ++k) hMee_int->Fill(mee);

    // ---- smart rounding replication for actual unweighted tree production
    const int nrep = smartRound(wpair, rng);
    sumNrepSmart += nrep;

    // update counter with accumulated nrep
    // (use SetBinContent since Fill would add 1, not nrep)
    hCounts->SetBinContent(4, (double)sumNrepSmart);

    if (nrep <= 0) continue; // accept–reject equivalent

    // prepare tree event with exactly the two electrons of the selected pair
    myevent.set_to_null();
    myevent.ntracks = 2;

    auto pushTrack = [&](const CandE& e) {
      myevent.pid.push_back(e.pdg);
      myevent.mass.push_back(0.000511);
      myevent.energy.push_back(e.p4.E());
      myevent.px.push_back(e.p4.Px());
      myevent.py.push_back(e.p4.Py());
      myevent.pz.push_back(e.p4.Pz());
      myevent.vx.push_back(pythia.event[e.idx].xProd());
      myevent.vy.push_back(pythia.event[e.idx].yProd());
      myevent.vz.push_back(pythia.event[e.idx].zProd());
    };

    pushTrack(e1);
    pushTrack(e2);

    // replicate fills
    for (int r = 0; r < nrep; ++r) {
      if (writtenEntries >= targetTreeEntries) break;
      tree->Fill();
      hMee_tree->Fill(mee);
      ++writtenEntries;
    }

    // store written entries counter in bin at x=4 (bin index 5)
    hCounts->SetBinContent(5, (double)writtenEntries);

    if (writtenEntries % 10 == 0) {
      std::cout << "written=" << writtenEntries
                << " tried=" << tried
                << " <wpair>~" << (sumWpair / std::max(1LL, (long long)hCounts->GetBinContent(3)))
                << "\n";
    }
  }

  // store tried as parameter
  TParameter<long long>* pTried = new TParameter<long long>("Ntried_pythiaNextOK", tried);
  TParameter<long long>* pOverall = new TParameter<long long>("Noverall_pythiaNextOK", hook->nCalls);

  // additional bookkeeping parameters
  TParameter<double>* pSumWpair   = new TParameter<double>("sum_wpair", sumWpair);
  TParameter<double>* pSumBRprod  = new TParameter<double>("sum_BRprod", sumBRprod);
  TParameter<long long>* pSumInt  = new TParameter<long long>("sum_intWeight", sumIntWeight);
  TParameter<long long>* pSumNrep = new TParameter<long long>("sum_nrepSmart", sumNrepSmart);

  // “total norm” you mentioned: Nev / BRD0^2
  // Here Nev = writtenEntries (tree entries), but you can also use tried or events-with-pair.
  TParameter<double>* pNormEntriesOverBRD0sq =
    new TParameter<double>("Ntree_over_BRD0sq", (double)writtenEntries / brD0sq);

  // cross section
  TParameter<double>* pSigmaGen = new TParameter<double>("sigmaGen_mb", pythia.info.sigmaGen());

  // Write everything
  fout->cd();
  tree->Write();

  hBR->Write();
  hWpair->Write();

  hMee_int->Write();
  hMee_wpair->Write();
  hMee_brprod->Write();
  hMee_tree->Write();

  hCounts->Write();
  hParentAll->Write();
  hParentPair->Write();

  hPtAllCharmE->Write();

  hMee_exact2->Write();
  hMee_star->Write();
  hMee_phenix_eta05->Write();
  hMee_phenix_y035->Write();
  hMee_phenix_phiacc->Write();
  
  hPt_star->Write();
  hPt_phenix_eta->Write();
  hPt_phenix_y->Write();
  hPt_phenix_phi->Write();
  
  hPTHat->Write();


  pBRD0->Write();
  pBRD0sq->Write();
  pOverall->Write();
  pTried->Write();
  pSumWpair->Write();
  pSumBRprod->Write();
  pSumInt->Write();
  pSumNrep->Write();
  pNormEntriesOverBRD0sq->Write();
  pSigmaGen->Write();

  fout->Close();

  std::cout << "Done.\n";
  std::cout << "overall=" << hook->nCalls<< " "<< overall << "\n";
  std::cout << "tried=" << tried << "\n";
  std::cout << "treeEntries=" << writtenEntries << "\n";
  std::cout << "BRD0=" << std::setprecision(10) << brD0 << "  BRD0^2=" << brD0sq << "\n";
  std::cout << "sum_wpair=" << std::setprecision(10) << sumWpair << "\n";
  std::cout << "sum_BRprod=" << std::setprecision(10) << sumBRprod << "\n";
  std::cout << "sum_intWeight=" << sumIntWeight << "\n";
  std::cout << "sum_nrepSmart=" << sumNrepSmart << "\n";
  std::cout << "sigmaGen(mb)=" << std::setprecision(8) << pythia.info.sigmaGen() << "\n";

  auto t1 = std::chrono::high_resolution_clock::now();
  auto dt = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  std::cout << "runtime(s)=" << dt << "\n";

  std::cout << "Hook calls = " << hook->nCalls
          << " vetoes = " << hook->nVeto
          << " accept = " << (hook->nCalls - hook->nVeto)
          << std::endl;


  return 0;
}
