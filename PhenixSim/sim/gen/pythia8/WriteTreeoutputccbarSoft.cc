// charm_softqcd_hook_pairs.cc
//
// Build a charm-only (no-bottom) SoftQCD sample using a UserHook,
// write a TTree in your vector-branch format,
// BUT only write the tree (and OSCAR block) when there are at least 2 electrons
// passing: |y| < 0.5 and pT > 0.2.
//
// Also writes several mee histograms for different selections + single-e pT histos,
// plus bookkeeping counters and sigmaGen.
//
// Compile (example):
//   g++ -O2 -std=c++17 charm_softqcd_hook_pairs.cc $(pythia8-config --cxxflags --libs) \
//       $(root-config --cflags --libs) -o charm_softqcd_hook_pairs
//
// Run:
//   ./charm_softqcd_hook_pairs 12345 5000
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// -----------------------------
// Data container for TTree (your format)
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

// -----------------------------
// UserHook: keep events with charm and veto any bottom
// -----------------------------
class CharmOnlyNoBottomHook : public UserHooks {
public:
  bool canVetoPartonLevel() override { return true; }

  bool doVetoPartonLevel(const Event& event) override {
    bool hasCharm  = false;
    bool hasBottom = false;

    for (int i = 0; i < event.size(); ++i) {
      const int id = std::abs(event[i].id());
      if (id == 4) hasCharm = true;
      if (id == 5) hasBottom = true;
      if (hasCharm && hasBottom) break;
    }

    // veto if no charm OR any bottom is present
    if (!hasCharm) return true;
    if (hasBottom) return true;
    return false;
  }
};

// -----------------------------
// Helpers
// -----------------------------
static inline bool isElectron(int pdg) { return (pdg == 11 || pdg == -11); }
static inline int chargeFromPdg(int pdg) { return (pdg > 0) ? -1 : +1; } // e-:11 -> -1; e+:-11 -> +1

static inline double wrapToPhiRange(double phi) {
  // match your snippet: keep in [-pi/2, 3pi/2)
  const double PI = TMath::Pi();
  while (phi < -PI/2) phi += 2.0*PI;
  while (phi >= 3.0*PI/2) phi -= 2.0*PI;
  return phi;
}
static inline bool inArm(double phi,
                         double phi_west_low, double phi_west_up,
                         double phi_east_low, double phi_east_up) {
  return ( (phi > phi_west_low && phi < phi_west_up) ||
           (phi > phi_east_low && phi < phi_east_up) );
}

struct ETrack {
  int pdg = 0;
  int q = 0;
  TLorentzVector p4;
  double pt = 0;
  double eta = 0;
  double y = 0;
};

static inline bool passTreeGate(const ETrack& t) {
  // REQUIRED by user for writing the tree:
  // at least 2 electrons must satisfy |y|<0.5 and pT>0.2
  return (t.pt > 0.2 && std::fabs(t.y) < 0.5);
}
static inline bool passSTAR(const ETrack& t) {
  return (t.pt > 0.2 && std::fabs(t.y) < 1.0);
}
static inline bool passPHENIX_simEta(const ETrack& t) {
  return (t.pt > 0.2 && std::fabs(t.eta) < 0.5);
}
static inline bool passPHENIX_y035(const ETrack& t) {
  return (t.pt > 0.2 && std::fabs(t.y) < 0.35);
}
static inline bool passPHENIX_phiDC_RICH(const ETrack& t) {
  // your constants
  const float phi_west_low = -3*TMath::Pi()/16, phi_west_up = 5*TMath::Pi()/16;
  const float phi_east_low = 11*TMath::Pi()/16, phi_east_up = 19*TMath::Pi()/16;
  const float k_DC   = 0.206f; // rad GeV/c
  const float k_RICH = 0.309f; // rad GeV/c

  const float phi = TMath::ATan2((float)t.p4.Py(), (float)t.p4.Px());
  const float pt  = (float)t.pt;
  const int   q   = t.q;

  // Use k_RICH for RICH and k_DC for DC
  float phi_rich = phi + q * k_RICH / pt;
  float phi_dc   = phi + q * k_DC   / pt;

  phi_rich = (float)wrapToPhiRange(phi_rich);
  phi_dc   = (float)wrapToPhiRange(phi_dc);

  if (!inArm(phi_rich, phi_west_low, phi_west_up, phi_east_low, phi_east_up)) return false;
  if (!inArm(phi_dc,   phi_west_low, phi_west_up, phi_east_low, phi_east_up)) return false;

  return true;
}

static inline void fillOSPairsMass(const std::vector<ETrack>& v, TH1D* h) {
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = i + 1; j < v.size(); ++j) {
      if (v[i].q * v[j].q >= 0) continue; // unlike-sign only
      h->Fill((v[i].p4 + v[j].p4).M());
    }
  }
}

int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <seed> <nevents_to_write>\n";
    return 1;
  }

  const std::string str_seed = argv[1];
  const int nevents_to_write = std::stoi(argv[2]);

  const bool IsWriteOscar = true;
  const int  nevt_max_try = 400000; // attempts cap (SoftQCD + charm-only hook can be sparse)
  const double pi = TMath::ACos(-1.0);

  std::cout << "lets begin\n";
  auto start = std::chrono::high_resolution_clock::now();

  // ----------------------------
  // Pythia setup
  // ----------------------------
  Pythia pythia;

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");

  // As requested: SoftQCD
  pythia.readString("SoftQCD:inelastic = on");

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + str_seed);
  pythia.readString("Next:numberCount = 100000000");

  // STAR-like tune bits you had
  pythia.readString("PDF:pSet = 17");
  pythia.readString("MultipartonInteractions:ecmRef = 200");
  pythia.readString("MultipartonInteractions:bprofile = 2");
  pythia.readString("MultipartonInteractions:pT0Ref = 1.40");
  pythia.readString("MultipartonInteractions:ecmPow = 0.135");
  pythia.readString("MultipartonInteractions:coreRadius = 0.56");
  pythia.readString("MultipartonInteractions:coreFraction = 0.78");
  pythia.readString("ColourReconnection:range = 5.4");

  // Charm-only hook (veto bottom)
  std::shared_ptr<CharmOnlyNoBottomHook> hook = std::make_shared<CharmOnlyNoBottomHook>();
  pythia.setUserHooksPtr(hook);

  pythia.init();
  std::cout << "pythia initialized (SoftQCD inelastic + charm-only hook)\n";

  // ----------------------------
  // ROOT output
  // ----------------------------
  TTree* tree = new TTree("T", "RECREATE");
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

  // -----------------------------
  // Histograms
  // -----------------------------
  TH1D* hPtAllCharmE = new TH1D("pt_all", "e^{#pm} p_{T} (all charm e);p_{T} (GeV/c);counts", 100, 0, 10);

  TH1D* hMee_exact2        = new TH1D("mee_exact2",        "m_{ee} (OS), exactly 2 charm e (tree gate not required);m_{ee} (GeV/c^{2});counts", 300, 0, 6);
  TH1D* hMee_star          = new TH1D("mee_star",          "m_{ee} (OS), STAR |y|<1 pT>0.2;m_{ee};counts",                                     300, 0, 6);
  TH1D* hMee_phenix_eta05  = new TH1D("mee_phenix_eta05",  "m_{ee} (OS), PHENIX sim |#eta|<0.5 pT>0.2;m_{ee};counts",                           300, 0, 6);
  TH1D* hMee_phenix_y035   = new TH1D("mee_phenix_y035",   "m_{ee} (OS), PHENIX |y|<0.35 pT>0.2;m_{ee};counts",                                 300, 0, 6);
  TH1D* hMee_phenix_phiacc = new TH1D("mee_phenix_phiacc", "m_{ee} (OS), PHENIX |y|<0.35 pT>0.2 + DC/RICH #phi cuts;m_{ee};counts",             300, 0, 6);

  TH1D* hPt_star        = new TH1D("pt_star",        "e^{#pm} p_{T} (STAR |y|<1);p_{T};counts",                  100, 0, 10);
  TH1D* hPt_phenix_eta  = new TH1D("pt_phenix_eta",  "e^{#pm} p_{T} (PHENIX sim |#eta|<0.5);p_{T};counts",       100, 0, 10);
  TH1D* hPt_phenix_y    = new TH1D("pt_phenix_y",    "e^{#pm} p_{T} (PHENIX |y|<0.35);p_{T};counts",             100, 0, 10);
  TH1D* hPt_phenix_phi  = new TH1D("pt_phenix_phi",  "e^{#pm} p_{T} (PHENIX |y|<0.35 + DC/RICH #phi cuts);p_{T};counts", 100, 0, 10);

  TH1D* hPTHat = new TH1D("pTHat", "pTHat Distribution;#hat{p}_{T} (GeV/c);counts", 100, 0, 10);
  TH1D* hNorm  = new TH1D("norm",  "Event counter;bin;counts", 4, -0.5, 3.5);
  // bin0: tried (pythia.next success)
  // bin1: passed tree gate (>=2 e with |y|<0.5 pT>0.2)
  // bin2: written (tree filled)
  // bin3: reserved

  TH1D* hSigmaGen = new TH1D("sigmaGen_mb", "Pythia sigmaGen; ; mb", 1, 0, 1);

  // ----------------------------
  // OSCAR output (optional)
  // ----------------------------
  std::ofstream file;
  if (IsWriteOscar) {
    file.open("oscar.particles.dat");
    file << "# OSC1999A\n";
    file << "# final_id_p_x\n";
    file << "# SimName 1.0\n";
    file << "#\n";
    file << "# Some comments...\n";
  }

  // ----------------------------
  // Event loop
  // ----------------------------
  int written = 0;

  for (int iev = 0; iev < nevt_max_try; ++iev) {

    if (written >= nevents_to_write) break;

    if (!pythia.next()) continue;
    hNorm->Fill(0);

    const int nparticles = pythia.event.size();

    // Collect all final-state charm electrons (mother in D-family list)
    std::vector<ETrack> allE;
    allE.reserve(16);

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
      if (pt <= 0.0) continue; // protect k_DC/pt etc

      // Basic pt cut (common baseline)
      if (pt < 0.2) continue;

      // Mother safety
      const int mother_index = pythia.event[j].mother1();
      if (mother_index <= 0 || mother_index >= nparticles) continue;
      const int mother_id = std::abs(pythia.event[mother_index].id());

      // D-family list (same as your original)
      const bool mother_D_meson =
          mother_id == 411   || mother_id == 421   || mother_id == 10411 || mother_id == 10421 ||
          mother_id == 413   || mother_id == 423   || mother_id == 10413 || mother_id == 10423 ||
          mother_id == 20413 || mother_id == 20423 || mother_id == 415   || mother_id == 425   ||
          mother_id == 431   || mother_id == 10431 || mother_id == 433   || mother_id == 10433 ||
          mother_id == 20433 || mother_id == 435;

      if (!mother_D_meson) continue;

      ETrack t;
      t.pdg = pdg;
      t.q   = chargeFromPdg(pdg);
      t.p4  = p4;
      t.pt  = pt;
      t.eta = p4.Eta();
      t.y   = p4.Rapidity();

      allE.push_back(t);

      // single-e pT (all charm e passing pt>0.2)
      hPtAllCharmE->Fill(t.pt);
    }

    // Fill "exact2" mee (OS only) independent of tree gate
    if (allE.size() == 2 && allE[0].q * allE[1].q < 0) {
      hMee_exact2->Fill((allE[0].p4 + allE[1].p4).M());
    }

    // Build selection lists for pair+single histos
    std::vector<ETrack> vSTAR, vPHX_eta, vPHX_y, vPHX_phi;
    vSTAR.reserve(allE.size());
    vPHX_eta.reserve(allE.size());
    vPHX_y.reserve(allE.size());
    vPHX_phi.reserve(allE.size());

    for (const auto& t : allE) {
      if (passSTAR(t)) { vSTAR.push_back(t); hPt_star->Fill(t.pt); }
      if (passPHENIX_simEta(t)) { vPHX_eta.push_back(t); hPt_phenix_eta->Fill(t.pt); }
      if (passPHENIX_y035(t)) { vPHX_y.push_back(t); hPt_phenix_y->Fill(t.pt); }
      if (passPHENIX_y035(t) && passPHENIX_phiDC_RICH(t)) { vPHX_phi.push_back(t); hPt_phenix_phi->Fill(t.pt); }
    }

    fillOSPairsMass(vSTAR,   hMee_star);
    fillOSPairsMass(vPHX_eta,hMee_phenix_eta05);
    fillOSPairsMass(vPHX_y,  hMee_phenix_y035);
    fillOSPairsMass(vPHX_phi,hMee_phenix_phiacc);

    // ----------------------------
    // TREE GATE: write tree ONLY if at least 2 electrons satisfy |y|<0.5, pT>0.2
    // ----------------------------
    int n_gate = 0;
    for (const auto& t : allE) if (passTreeGate(t)) ++n_gate;

    if (n_gate < 2) continue;
    hNorm->Fill(1);

    // Fill tree with ONLY the gated electrons (|y|<0.5, pT>0.2) to match the requirement
    myevent.set_to_null();
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
      if (pt < 0.2) continue;

      const double y = p4.Rapidity();
      if (std::fabs(y) >= 0.5) continue;

      // Mother safety
      const int mother_index = pythia.event[j].mother1();
      if (mother_index <= 0 || mother_index >= nparticles) continue;
      const int mother_id = std::abs(pythia.event[mother_index].id());

      const bool mother_D_meson =
          mother_id == 411   || mother_id == 421   || mother_id == 10411 || mother_id == 10421 ||
          mother_id == 413   || mother_id == 423   || mother_id == 10413 || mother_id == 10423 ||
          mother_id == 20413 || mother_id == 20423 || mother_id == 415   || mother_id == 425   ||
          mother_id == 431   || mother_id == 10431 || mother_id == 433   || mother_id == 10433 ||
          mother_id == 20433 || mother_id == 435;

      if (!mother_D_meson) continue;

      myevent.pid.push_back(pdg);
      myevent.mass.push_back(0.000511);
      myevent.energy.push_back(E);
      myevent.px.push_back(Px);
      myevent.py.push_back(Py);
      myevent.pz.push_back(Pz);
      myevent.vx.push_back(pythia.event[j].xProd());
      myevent.vy.push_back(pythia.event[j].yProd());
      myevent.vz.push_back(pythia.event[j].zProd());
    }

    myevent.ntracks = (int)myevent.pid.size();
    if (myevent.ntracks < 2) continue; // just in case (should not happen)

    hPTHat->Fill(pythia.info.pTHat());

    tree->Fill();
    hNorm->Fill(2);

    if (IsWriteOscar && myevent.ntracks > 0) {
      file << 0 << "\t" << myevent.ntracks << "\n";
      for (int i = 0; i < myevent.ntracks; ++i) {
        file << (i+1) << "\t"
             << myevent.pid[i] << "\t"
             << 0 << "\t"
             << myevent.px[i] << "\t"
             << myevent.py[i] << "\t"
             << myevent.pz[i] << "\t"
             << myevent.energy[i] << "\t"
             << myevent.mass[i] << "\t"
             << myevent.vx[i]*std::pow(10,12) << "\t"
             << myevent.vy[i]*std::pow(10,12) << "\t"
             << myevent.vz[i]*std::pow(10,12) << "\t"
             << 0 << "\n";
      }
      file << 0 << "\t" << 0 << "\n";
    }

    ++written;
    if (written % 1 == 0) std::cout << written << "\tCompleted\n";
  }

  // cross section
  hSigmaGen->SetBinContent(1, pythia.info.sigmaGen()); // mb
  std::cout << std::setprecision(6) << pythia.info.sigmaGen() << "\n";

  // (Optional) convert hPtAllCharmE to invariant-yield-like (1/(2π pT dpT) dN/dpT)
  TH1D* pt_spectra = (TH1D*)hPtAllCharmE->Clone("pt_spectra");
  pt_spectra->Reset("ICESM");
  const int nb = hPtAllCharmE->GetNbinsX();
  for (int ib = 1; ib <= nb; ++ib) {
    const double c = hPtAllCharmE->GetBinContent(ib);
    const double e = hPtAllCharmE->GetBinError(ib);
    const double w = hPtAllCharmE->GetBinWidth(ib);
    const double x = hPtAllCharmE->GetBinCenter(ib);
    const double denom = 2.0 * pi * w * x;
    if (denom > 0) {
      pt_spectra->SetBinContent(ib, c / denom);
      pt_spectra->SetBinError(ib,   e / denom);
    }
  }

  // Write output
  TFile* fout = new TFile("tree_out.root", "RECREATE");
  fout->cd();

  tree->Write();
  hPtAllCharmE->Write();
  pt_spectra->Write();

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
  hNorm->Write();
  hSigmaGen->Write();

  fout->Close();

  if (IsWriteOscar) file.close();

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << duration.count() << "\n";

  return 0;
}
