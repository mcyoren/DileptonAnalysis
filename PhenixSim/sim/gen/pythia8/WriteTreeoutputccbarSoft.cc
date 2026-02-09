// WriteTreeoutputccbarSoft.cc
//
// SoftQCD + charm-only (no bottom) via UserHook
// Force charm hadrons -> e channels only
// For each event: pick one OS ee pair in gate (|y|<0.5, pT>0.2)
// Replicate event N times where N ~ BR(H1)*BR(H2)/BR(D0)^2 using unbiased smart rounding
// Write TTree ONLY when gate has >=2 electrons and replication count > 0
//
// ADDITION:
// - classify ccbar "source" into 5 bins (event-level tag, computed once per event):
//     (a) s-channel Flavor Creation      : qqbar -> ccbar
//     (b) t-channel Flavor Creation      : gg -> ccbar
//     (c) Flavor Excitation (your def)   : gg -> ccbar + g  (cc and extra g share same 2 mothers = incoming gg)
//     (d) Gluon Splitting                : g -> ccbar (a gluon has both c and cbar as daughters)
//     (e) OTHER
//
// - ALL mass & pT histograms become 2D: (observable vs source).
// - STAR/PHENIX single-e pT are weighted by BR(parent)/BR(D0) (and by pythia event weight).
// - Pair observables (m_ee, pT_ee) are weighted by wpair=BR(H1)BR(H2)/BR(D0)^2 (and by event weight).
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
#include <TH2D.h>
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



static inline void sort2(int& a, int& b) { if (a > b) std::swap(a,b); }

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

static inline bool inPhenixArm(double phi) {
  const double PI = TMath::Pi();
  const double phi_west_low = -3*PI/16.0, phi_west_up = 5*PI/16.0;
  const double phi_east_low = 11*PI/16.0, phi_east_up = 19*PI/16.0;

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
// ---- helpers to iterate daughters of a particle i ----
static std::vector<int> getDaughters(const Event& ev, int i)
{
  std::vector<int> out;
  int d1 = ev[i].daughter1();
  int d2 = ev[i].daughter2();

  // Pythia convention: 0 means none
  if (d1 <= 0 || d2 <= 0) return out;

  // IMPORTANT: if d1 > d2 it's not a contiguous daughter range.
  // Do NOT swap. Treat as unusable for simple "range daughter" logic.
  if (d1 > d2) return out;

  out.reserve(d2 - d1 + 1);
  for (int d = d1; d <= d2; ++d) out.push_back(d);
  return out;
}


static bool getIncomingPartonsIdx(const Event& ev, int& iA, int& iB)
{
  iA = -1; iB = -1;
  for (int i = 0; i < ev.size(); ++i) {
    int st = ev[i].status();
    if (st == -21 || st == -22) {
      if (iA < 0) iA = i;
      else { iB = i; return true; }
    }
  }
  return false;
}

static bool hasDirectCCbarKids(const Event& ev,
                                      const std::vector<std::vector<int>>& kids,
                                      int g, int& cIdx, int& cbIdx)
{
  cIdx = -1; cbIdx = -1;

  for (int ch : kids[g]) {
    const int id = ev[ch].id();
    if (id != 4 && id != -4) continue;

    // reject incoming legs (PDF history can make incoming charm look like a "daughter")
    const int st = ev[ch].status();
    if (st == -21 || st == -22) continue;

    if (id == 4)  cIdx  = ch;
    if (id == -4) cbIdx = ch;

    if (cIdx >= 0 && cbIdx >= 0) return true;
  }
  return false;
}

static inline bool mothersAreIncomingPair(const Event& ev, int i, int in1, int in2)
{
  int a = ev[i].mother1();
  int b = ev[i].mother2();
  if (a <= 0 || b <= 0) return false;
  sort2(a,b);
  int x=in1, y=in2; sort2(x,y);
  return (a==x && b==y);
}

static void getHardOutgoing(const Event& ev, int in1, int in2, std::vector<int>& out)
{
  out.clear();
  for (int i=0;i<ev.size();++i) {
    if (ev[i].status() != -23) continue;
    if (!mothersAreIncomingPair(ev, i, in1, in2)) continue;
    out.push_back(i);
  }
}


static void printParticleOneLine(const Event& ev, int i, const char* tag="")
{
  std::cout
    << tag
    << " idx=" << i
    << " id=" << ev[i].id()
    << " st=" << ev[i].status()
    << " m=(" << ev[i].mother1() << "," << ev[i].mother2() << ")"
    << " d=(" << ev[i].daughter1() << "," << ev[i].daughter2() << ")"
    << " pT=" << std::fixed << std::setprecision(3) << ev[i].pT()
    << " y="  << std::fixed << std::setprecision(3) << ev[i].y()
    << "\n";
}


// Build mother -> children list (for event record indices)
static void buildChildrenMap(const Event& ev, std::vector<std::vector<int>>& kids)
{
  kids.assign(ev.size(), std::vector<int>());
  for (int i = 0; i < ev.size(); ++i) {
    int m1 = ev[i].mother1();
    int m2 = ev[i].mother2();
    if (m1 > 0 && m1 < ev.size()) kids[m1].push_back(i);
    if (m2 > 0 && m2 < ev.size() && m2 != m1) kids[m2].push_back(i);
  }
}

static void debugPrintIncomingAndCharmTopology(const Pythia& pythia, int iev, int maxPrintEvents=100)
{
  if (iev >= maxPrintEvents) return;

  const Event& ev = pythia.event;

  int in1=-1, in2=-1;
  if (!getIncomingPartonsIdx(ev, in1, in2)) return;

  std::cout << "\n================ EVENT " << iev
            << "  weight=" << pythia.info.weight()
            << "  pTHat=" << pythia.info.pTHat()
            << " ================\n";

  std::cout << "Incoming partons:\n";
  printParticleOneLine(ev, in1, "  [IN] ");
  printParticleOneLine(ev, in2, "  [IN] ");

  // Build child map once (robust; does not rely on daughter1..daughter2 ranges)
  std::vector<std::vector<int>> kids;
  buildChildrenMap(ev, kids);

  // Print hard outgoing legs: status -23 and mothers are the incoming pair
  std::cout << "Hard outgoing (-23, mothers=incoming pair):\n";
  std::vector<int> hard;
  getHardOutgoing(ev, in1, in2, hard);
  if (hard.empty()) std::cout << "  (none)\n";
  for (int i : hard) printParticleOneLine(ev, i, "  [H] ");

  // Print ONLY hard gluons that directly split to c+cbar (true GS candidates)
  std::cout << "Hard gluons with direct (c,cbar) kids (GS candidates):\n";
  bool found = false;
  for (int i : hard) {
    if (ev[i].id() != 21) continue;          // only gluons
    int c=-1, cb=-1;
    if (!hasDirectCCbarKids(ev, kids, i, c, cb)) continue; // kids-based direct test
    found = true;
    printParticleOneLine(ev, i, "  [gH] ");
    printParticleOneLine(ev, c,  "     ");
    printParticleOneLine(ev, cb, "     ");
  }
  if (!found) std::cout << "  (none)\n";

  // (Optional) If you still want to see ISR radiator "fake splitters", print them separately:
  std::cout << "Non-hard non-incoming gluons with direct (c,cbar) kids (often ISR bookkeeping):\n";
  found = false;
  for (int i = 0; i < ev.size(); ++i) {
    if (ev[i].id() != 21) continue;
    if (i == in1 || i == in2) continue;
    if (ev[i].status() == -23) continue;    // skip hard; already printed above
    int c=-1, cb=-1;
    if (!hasDirectCCbarKids(ev, kids, i, c, cb)) continue;
    found = true;
    printParticleOneLine(ev, i, "  [gISR] ");
    printParticleOneLine(ev, c,  "        ");
    printParticleOneLine(ev, cb, "        ");
  }
  std::cout << "Any non-incoming gluon with strict direct (c,cbar) kids (shower GS candidates):\n";
  for (int i=0;i<ev.size();++i) {
    if (ev[i].id()!=21) continue;
    if (ev[i].status()==-21 || ev[i].status()==-22) continue; // not incoming
    int c=-1, cb=-1;
    if (!hasDirectCCbarKids(ev, kids, i, c, cb)) continue;
    found = true;
    printParticleOneLine(ev, i, "  [gSplit] ");
    printParticleOneLine(ev, c,  "     ");
    printParticleOneLine(ev, cb, "     ");
  }
  if (!found) std::cout << "  (none)\n";
  
}

// ----------------------------------------------------------------------
// ccbar source classification (5 bins)
// ----------------------------------------------------------------------
enum CcSource {
  SRC_S_FLAVOR_CREATION = 0, // (a) qqbar -> cc
  SRC_T_FLAVOR_CREATION = 1, // (b) gg -> cc
  SRC_FLAVOR_EXCITATION = 2, // (c) gg -> cc + g  (your requested def)
  SRC_GLUON_SPLITTING   = 3, // (d) g -> cc
  SRC_OTHER             = 4  // (e)
};

static inline const char* srcName(int s) {
  switch(s) {
    case SRC_S_FLAVOR_CREATION: return "S_FC (qqbar->cc)";
    case SRC_T_FLAVOR_CREATION: return "T_FC (gg->cc)";
    case SRC_FLAVOR_EXCITATION: return "FE (gg->ccg)";
    case SRC_GLUON_SPLITTING:   return "GS (g->cc)";
    case SRC_OTHER:             return "OTHER";
    default:                    return "UNKNOWN";
  }
}

// Check if a gluon has both c and cbar as daughters (g -> c cbar) in the record
static bool findGluonSplittingPair(const Event& ev,
                                  const std::vector<std::vector<int>>& kids,
                                  int& outG, int& outC, int& outCbar)
{
  outG = outC = outCbar = -1;
  for (int g = 0; g < ev.size(); ++g) {
    if (ev[g].id() != 21) continue;
    int foundC = -1, foundCb = -1;
    for (int ch : kids[g]) {
      int id = ev[ch].id();
      if (id != 4 && id != -4) continue;
      int st = ev[ch].status();
      if (st == -21 || st == -22) continue; // reject incoming
      if (id == 4) foundC = ch;
      if (id == -4) foundCb = ch;
    }
  }
  return false;
}

// Find a c and cbar that share the same two mothers (unordered) and return them
static bool findPairSameTwoMothers(const Event& ev, int& outC, int& outCbar, int& m1, int& m2)
{
  outC = outCbar = -1;
  m1 = m2 = -1;

  // store c candidates by their (sorted) mother pair
  struct Key { int a,b; };
  struct KeyLess { bool operator()(const Key& x, const Key& y) const {
    if (x.a != y.a) return x.a < y.a;
    return x.b < y.b;
  }};

  std::map<Key, std::vector<int>, KeyLess> cByMoms;
  std::map<Key, std::vector<int>, KeyLess> cbByMoms;

  for (int i = 0; i < ev.size(); ++i) {
    const int id = ev[i].id();
    if (id != 4 && id != -4) continue;
    int a = ev[i].mother1();
    int b = ev[i].mother2();
    if (a <= 0 || b <= 0) continue;
    sort2(a,b);
    Key k{a,b};
    if (id == 4)  cByMoms[k].push_back(i);
    if (id == -4) cbByMoms[k].push_back(i);
  }

  // prefer pairs whose mothers are incoming partons
  int inA=-1, inB=-1;
  const bool hasIn = getIncomingPartonsIdx(ev, inA, inB);
  int siA=inA, siB=inB;
  if (hasIn) sort2(siA, siB);

  if (hasIn) {
    Key kin{siA, siB};
    auto itc  = cByMoms.find(kin);
    auto itcb = cbByMoms.find(kin);
    if (itc != cByMoms.end() && itcb != cbByMoms.end() &&
        !itc->second.empty() && !itcb->second.empty()) {
      outC    = itc->second.front();
      outCbar = itcb->second.front();
      m1 = kin.a; m2 = kin.b;
      return true;
    }
  }

  // otherwise take any available mother-pair match
  for (auto& kv : cByMoms) {
    auto itcb = cbByMoms.find(kv.first);
    if (itcb == cbByMoms.end()) continue;
    if (kv.second.empty() || itcb->second.empty()) continue;

    outC    = kv.second.front();
    outCbar = itcb->second.front();
    m1 = kv.first.a; m2 = kv.first.b;
    return true;
  }

  return false;
}

// FE check: existence of extra gluon with the same two mothers as the ccbar pair mothers
static bool hasExtraGluonSameMothers(const Event& ev, int mom1, int mom2, int cIdx, int cbIdx)
{
  int a = mom1, b = mom2;
  sort2(a,b);
  for (int i = 0; i < ev.size(); ++i) {
    if (i == cIdx || i == cbIdx) continue;
    if (ev[i].id() != 21) continue;
    int m1 = ev[i].mother1();
    int m2 = ev[i].mother2();
    if (m1 <= 0 || m2 <= 0) continue;
    sort2(m1,m2);
    if (m1 == a && m2 == b) return true;
  }
  return false;
}

// Event-level classification: compute once per event.
static int classifyEventCcSource(const Pythia& pythia)
{
  const Event& ev = pythia.event;

  int in1=-1, in2=-1;
  if (!getIncomingPartonsIdx(ev, in1, in2)) return SRC_OTHER;

  const int idIn1 = ev[in1].id();
  const int idIn2 = ev[in2].id();
  const int aIn1  = std::abs(idIn1);
  const int aIn2  = std::abs(idIn2);

  // (c) Flavor excitation: incoming charm present (cg -> cg, cq -> cq, etc.)
  // This matches what your debug print shows for many events.
  if (aIn1 == 4 || aIn2 == 4) return SRC_FLAVOR_EXCITATION;

  // Build child map once
  std::vector<std::vector<int>> kids;
  buildChildrenMap(ev, kids);

  // Hard outgoing legs
  std::vector<int> hard;
  getHardOutgoing(ev, in1, in2, hard);

  // ---- (a)/(b) Flavor Creation: hard outgoing contains BOTH c and cbar
  bool hardHasC = false, hardHasCb = false;
  for (int i : hard) {
    if (ev[i].id() == 4)  hardHasC  = true;
    if (ev[i].id() == -4) hardHasCb = true;
  }
  if (hardHasC && hardHasCb) {
    // qqbar -> ccbar
    if (idIn1 == -idIn2 && (aIn1>=1 && aIn1<=5)) return SRC_S_FLAVOR_CREATION;
    // gg -> ccbar
    if (aIn1 == 21 && aIn2 == 21) return SRC_T_FLAVOR_CREATION;
    return SRC_OTHER;
  }

  // ---- (d) Gluon splitting: a HARD outgoing gluon (-23) has direct ccbar kids
  for (int i : hard) {
    if (ev[i].id() != 21) continue;
    int c=-1, cb=-1;
    if (hasDirectCCbarKids(ev, kids, i, c, cb)) {
      return SRC_GLUON_SPLITTING;
    }
  }

  // Optional: if you still want a "gg -> ccbar + g" category,
  // you need a hard 2->3 record (rare in this soft setup). You can attempt:
  // If incoming gg AND hard contains a charm OR anticharm AND a hard gluon, call FE.
  if (aIn1 == 21 && aIn2 == 21) {
    bool hardHasOneCharm = false, hardHasG = false;
    for (int i : hard) {
      if (std::abs(ev[i].id()) == 4) hardHasOneCharm = true;
      if (ev[i].id() == 21)          hardHasG = true;
    }
    if (hardHasOneCharm && hardHasG) return SRC_FLAVOR_EXCITATION;
  }
  // ---- (c') FE from ISR excitation: g -> ccbar where exactly one charm is a beam leg
  auto isBeamLegCharm = [&](int st){
    return (st == -31 || st == -33);
  };

  for (int g = 0; g < ev.size(); ++g) {
    if (ev[g].id() != 21) continue;
    if (ev[g].status() == -21 || ev[g].status() == -22) continue; // skip incoming gluons

    int c=-1, cb=-1;
    if (!hasDirectCCbarKids(ev, kids, g, c, cb)) continue;

    const bool cBeam  = isBeamLegCharm(ev[c].status());
    const bool cbBeam = isBeamLegCharm(ev[cb].status());

    // XOR: exactly one is beam-leg charm => excitation
    if (cBeam ^ cbBeam) return SRC_FLAVOR_EXCITATION;
  }

  // ---- (d') GS from shower: g -> ccbar where both charm legs are timelike (common: -51)
  for (int g = 0; g < ev.size(); ++g) {
    if (ev[g].id() != 21) continue;
    if (ev[g].status() == -21 || ev[g].status() == -22) continue;

    int c=-1, cb=-1;
    if (!hasDirectCCbarKids(ev, kids, g, c, cb)) continue;

    if (ev[c].status() == -51 && ev[cb].status() == -51)
      return SRC_GLUON_SPLITTING;
  }

  return SRC_OTHER;
}

// -----------------------------
// main
// -----------------------------
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

  std::shared_ptr<CharmOnlyNoBottomHook> hook = std::make_shared<CharmOnlyNoBottomHook>();
  pythia.setUserHooksPtr(hook);

  // ----------------------------
  // Charm parents to force into e± channels
  // Keep 421 (D0) as reference!
  // ----------------------------
  const std::vector<int> charmParents = {
    421, 411, 431, 4122,
    10421,10411, 423,413,
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

  // Counters
  TH1D* hCounts = new TH1D("counts", "Counters;bin;value", 6, -0.5, 5.5);
  // 0 tried, 1 >=2 gated e, 2 OS pair, 3 sum_nrep, 4 written, 5 reserved

  // Pair weight distribution
  TH1D* hWpair = new TH1D("w_pair", "w_{pair}=BRprod/BRD0^{2};w_{pair};counts", 200, 0, 5);

  // Event source counts (sanity!)
  TH1D* hSrcCount = new TH1D("src_count", "Event source tag;source;events", 5, -0.5, 4.5);
  for (int b=1;b<=5;++b) hSrcCount->GetXaxis()->SetBinLabel(b, srcName(b-1));
  hSrcCount->LabelsOption("v","X");

  // ----------------------------
  // ALL physics histograms are 2D: x=observable, y=source
  // ----------------------------
  auto make2D = [&](const char* name, const char* title,
                    int nx, double xlo, double xhi) -> TH2D* {
    TH2D* h = new TH2D(name, title, nx, xlo, xhi, 5, -0.5, 4.5);
    for (int b=1;b<=5;++b) h->GetYaxis()->SetBinLabel(b, srcName(b-1));
    h->LabelsOption("v","Y");
    return h;
  };

  // m_ee in various acceptances (weighted by wpair*eventWeight)
  TH2D* hMee_exact2_src       = make2D("mee_exact2_src",       "m_{ee} (OS), exactly 2 charm e (no acc);m_{ee};source", 300,0,6);
  TH2D* hMee_star_src         = make2D("mee_star_src",         "m_{ee} (OS), STAR |y|<1 p_{T}>0.2; m_{ee};source",     300,0,6);
  TH2D* hMee_phenix_eta05_src = make2D("mee_phenix_eta05_src", "m_{ee} (OS), |#eta|<0.5 p_{T}>0.2; m_{ee};source",     300,0,6);
  TH2D* hMee_phenix_y035_src  = make2D("mee_phenix_y035_src",  "m_{ee} (OS), |y|<0.35 p_{T}>0.2; m_{ee};source",        300,0,6);
  TH2D* hMee_phenix_phiacc_src= make2D("mee_phenix_phiacc_src","m_{ee} (OS), |y|<0.35 p_{T}>0.2 + #phi acc; m_{ee};source",300,0,6);

  // Main gate: mee, pair pT, and "tree replicated" mee
  TH2D* hMee_gate_wpair_src = make2D("mee_gate_wpair_src", "m_{ee} (gate) weighted by wpair; m_{ee};source", 300,0,6);
  TH2D* hMee_gate_tree_src  = make2D("mee_gate_tree_src",  "m_{ee} (gate) tree replications (unweighted); m_{ee};source", 300,0,6);
  TH2D* hPairPt_gate_src    = make2D("pairpt_gate_src",    "p_{T}^{ee} (gate) weighted by wpair; p_{T}^{ee};source", 120,0,12);

  // Single-e pT spectra (weighted by BR(parent)/BR(D0))
  TH2D* hPt_star_src         = make2D("pt_star_src",         "e p_{T} (STAR |y|<1) weighted by BR/BR(D0); p_{T};source", 100,0,10);
  TH2D* hPt_phenix_eta_src   = make2D("pt_phenix_eta_src",   "e p_{T} (|#eta|<0.5) weighted by BR/BR(D0); p_{T};source", 100,0,10);
  TH2D* hPt_phenix_y_src     = make2D("pt_phenix_y_src",     "e p_{T} (|y|<0.35) weighted by BR/BR(D0); p_{T};source",   100,0,10);
  TH2D* hPt_phenix_phi_src   = make2D("pt_phenix_phi_src",   "e p_{T} (|y|<0.35 + #phi acc) weighted by BR/BR(D0); p_{T};source",100,0,10);

  // pTHat (2D too, just to be consistent)
  TH2D* hPTHat_src = make2D("pTHat_src", "pTHat distribution;#hat{p}_{T};source", 100,0,10);

  // store constants
  TParameter<double>* pBRD0   = new TParameter<double>("BRD0", brD0);
  TParameter<double>* pBRD0sq = new TParameter<double>("BRD0sq", brD0sq);

  // ----------------------------
  // Loop
  // ----------------------------
  long long tried = 0;
  long long overall = 0;
  long long writtenEntries = 0;
  long long sumNrepSmart = 0;

  double sumWpair = 0.0;

  while (writtenEntries < targetTreeEntries) {

    overall++;
    if (!pythia.next()) continue;
    ++tried;
    hCounts->Fill(0);

    const double wEvt = pythia.info.weight(); // "pp weight" safeguard
    const int srcEvt = classifyEventCcSource(pythia);
    hSrcCount->Fill(srcEvt, wEvt);
    
    static int iev_debug = 0;
    debugPrintIncomingAndCharmTopology(pythia, iev_debug, 25);
    if(iev_debug < 25) std::cout<<"EVENT CLASS IS "<<srcName(srcEvt)<<"\n";
    ++iev_debug;


    hPTHat_src->Fill(pythia.info.pTHat(), srcEvt, wEvt);

    const int nparticles = pythia.event.size();

    std::vector<CandE> allCharm;  allCharm.reserve(16);
    std::vector<CandE> star;      star.reserve(16);
    std::vector<CandE> phen_eta;  phen_eta.reserve(16);
    std::vector<CandE> phen_y;    phen_y.reserve(16);
    std::vector<CandE> phen_phi;  phen_phi.reserve(16);
    std::vector<CandE> gated;     gated.reserve(8);

    // Collect electrons from charm parents
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

      const int parent = findCharmParent(pythia.event, j, charmParentSet);
      if (parent == 0) continue;
      c.parent = parent;

      // For single-e pT spectra: weight by BR(parent)/BR(D0) and by event weight
      const double br1 = brE.count(parent) ? brE[parent] : 0.0;
      if (br1 <= 0) continue;
      const double wSingle = wEvt * (br1 / brD0);

      allCharm.push_back(c);

      // apply pT>0.2 for your acceptance lists (as before)
      if (pt < 0.2) continue;

      if (std::fabs(c.y) < 1.0) {
        star.push_back(c);
        hPt_star_src->Fill(pt, srcEvt, wSingle);
      }

      if (std::fabs(c.eta) < 0.5) {
        phen_eta.push_back(c);
        hPt_phenix_eta_src->Fill(pt, srcEvt, wSingle);
      }

      if (std::fabs(c.y) < 0.35) {
        phen_y.push_back(c);
        hPt_phenix_y_src->Fill(pt, srcEvt, wSingle);

        if (passPhenixPhiAcc(Px, Py, c.q, pt)) {
          phen_phi.push_back(c);
          hPt_phenix_phi_src->Fill(pt, srcEvt, wSingle);
        }
      }

      if (std::fabs(c.y) < 0.5) {
        gated.push_back(c);
      }
    }

    // exactly 2 charm electrons (no acc), OS pair -> mee_exact2_src (pair weight)
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
        if (br1 > 0 && br2 > 0) {
          const double wpair = (br1*br2) / brD0sq;
          hMee_exact2_src->Fill(mee, srcEvt, wEvt*wpair);
        }
      }
    }

    // STAR mee
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
        if (br1 > 0 && br2 > 0) {
          const double wpair = (br1*br2) / brD0sq;
          hMee_star_src->Fill(mee, srcEvt, wEvt*wpair);
        }
      }
    }

    // PHENIX eta mee
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
        if (br1 > 0 && br2 > 0) {
          const double wpair = (br1*br2) / brD0sq;
          hMee_phenix_eta05_src->Fill(mee, srcEvt, wEvt*wpair);
        }
      }
    }

    // PHENIX y mee
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
        if (br1 > 0 && br2 > 0) {
          const double wpair = (br1*br2) / brD0sq;
          hMee_phenix_y035_src->Fill(mee, srcEvt, wEvt*wpair);
        }
      }
    }

    // PHENIX phi mee
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
        if (br1 > 0 && br2 > 0) {
          const double wpair = (br1*br2) / brD0sq;
          hMee_phenix_phiacc_src->Fill(mee, srcEvt, wEvt*wpair);
        }
      }
    }

    // MAIN gate for tree + dedicated gate histos
    if ((int)gated.size() < 2) continue;
    hCounts->Fill(1);

    int ia=-1, ib=-1;
    if (!pickBestOSPair(gated, ia, ib)) continue;
    hCounts->Fill(2);

    const CandE& e1 = gated[ia];
    const CandE& e2 = gated[ib];
    const double mee = (e1.p4 + e2.p4).M();
    const double pairPt = (e1.p4 + e2.p4).Pt();

    const int H1 = e1.parent;
    const int H2 = e2.parent;

    const double br1 = brE.count(H1) ? brE[H1] : 0.0;
    const double br2 = brE.count(H2) ? brE[H2] : 0.0;
    if (br1 <= 0 || br2 <= 0) continue;

    const double brprod = br1 * br2;
    const double wpair  = brprod / brD0sq;

    hWpair->Fill(wpair, wEvt);
    sumWpair += wpair;

    hMee_gate_wpair_src->Fill(mee, srcEvt, wEvt*wpair);
    hPairPt_gate_src->Fill(pairPt, srcEvt, wEvt*wpair);

    // replicate according to wpair (BR scaling), include event weight only in hist weights (tree is unweighted)
    const int nrep = smartRound(wpair, rng);
    sumNrepSmart += nrep;
    hCounts->SetBinContent(4, (double)sumNrepSmart);

    if (nrep <= 0) continue;

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

    for (int r = 0; r < nrep; ++r) {
      if (writtenEntries >= targetTreeEntries) break;
      tree->Fill();
      hMee_gate_tree_src->Fill(mee, srcEvt, 1.0); // tree entries are unweighted counts
      ++writtenEntries;
    }

    hCounts->SetBinContent(5, (double)writtenEntries);

    if (writtenEntries % 50 == 0) {
      std::cout << "written=" << writtenEntries
                << " tried=" << tried
                << " src=" << srcName(srcEvt)
                << "\n";
    }
  }

  // Parameters
  TParameter<long long>* pTried = new TParameter<long long>("Ntried_pythiaNextOK", tried);
  TParameter<long long>* pOverall = new TParameter<long long>("Noverall_hookCalls", hook->nCalls);
  TParameter<long long>* pSumNrep = new TParameter<long long>("sum_nrepSmart", sumNrepSmart);
  TParameter<double>* pSumWpair   = new TParameter<double>("sum_wpair", sumWpair);
  TParameter<double>* pSigmaGen   = new TParameter<double>("sigmaGen_mb", pythia.info.sigmaGen());

  // Write
  fout->cd();
  tree->Write();

  hBR->Write();
  hCounts->Write();
  hWpair->Write();
  hSrcCount->Write();

  hMee_exact2_src->Write();
  hMee_star_src->Write();
  hMee_phenix_eta05_src->Write();
  hMee_phenix_y035_src->Write();
  hMee_phenix_phiacc_src->Write();

  hMee_gate_wpair_src->Write();
  hMee_gate_tree_src->Write();
  hPairPt_gate_src->Write();

  hPt_star_src->Write();
  hPt_phenix_eta_src->Write();
  hPt_phenix_y_src->Write();
  hPt_phenix_phi_src->Write();

  hPTHat_src->Write();

  pBRD0->Write();
  pBRD0sq->Write();
  pOverall->Write();
  pTried->Write();
  pSumNrep->Write();
  pSumWpair->Write();
  pSigmaGen->Write();

  fout->Close();

  std::cout << "Done.\n";
  std::cout << "overall(hook calls)=" << hook->nCalls << " overall(loop)=" << overall << "\n";
  std::cout << "tried=" << tried << "\n";
  std::cout << "treeEntries=" << writtenEntries << "\n";
  std::cout << "BRD0=" << std::setprecision(10) << brD0 << "  BRD0^2=" << brD0sq << "\n";
  std::cout << "sum_wpair=" << std::setprecision(10) << sumWpair << "\n";
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
