// WriteLambda.cc
//
// Very simple pp 200 GeV SoftQCD generator.
// Writes pT distributions of:
//   Lambda + anti-Lambda
//   K0S
//   Lambda_c + anti-Lambda_c
//   D0 + anti-D0
//
// Same run structure as your previous code:
//
// Build:
//   g++ -O2 -std=c++17 WriteLambda.cc -o WriteLambda \
//     $(root-config --cflags --libs) $(pythia8-config --cxxflags --libs)
//
// Run:
//   ./WriteLambda <seed> <events_to_analyze>
//
// Example:
//   ./WriteLambda 12345 1000000

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>

#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// ------------------------------------------------------------
// Build mother -> children list
// ------------------------------------------------------------
static void buildChildrenMap(const Event& ev, std::vector<std::vector<int>>& kids)
{
  kids.assign(ev.size(), std::vector<int>());

  for (int i = 0; i < ev.size(); ++i) {
    int m1 = ev[i].mother1();
    int m2 = ev[i].mother2();

    if (m1 > 0 && m1 < ev.size()) {
      kids[m1].push_back(i);
    }

    if (m2 > 0 && m2 < ev.size() && m2 != m1) {
      kids[m2].push_back(i);
    }
  }
}

// ------------------------------------------------------------
// Avoid double-counting copied particles in the event record.
// We fill only the last copy: no daughter with the same abs(PDG).
// ------------------------------------------------------------
static bool hasSameAbsDaughter(const Event& ev,
                               const std::vector<std::vector<int>>& kids,
                               int i)
{
  const int idAbs = std::abs(ev[i].id());

  for (int d : kids[i]) {
    if (d <= 0 || d >= ev.size()) continue;

    if (std::abs(ev[d].id()) == idAbs) {
      return true;
    }
  }

  return false;
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
int main(int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <seed> <events_to_analyze>\n";
    return 1;
  }

  const int seed = std::stoi(argv[1]);
  const long long targetEvents = std::stoll(argv[2]);

  std::cout << "seed = " << seed
            << "  events_to_analyze = " << targetEvents << "\n";

  auto t0 = std::chrono::high_resolution_clock::now();

  // ------------------------------------------------------------
  // Pythia setup: same pp 200 GeV SoftQCD setting
  // ------------------------------------------------------------
  Pythia pythia;

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200");

  pythia.readString("SoftQCD:inelastic = on");

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 1000000");

 
  //Pythia8 tune for STAR (2110.09447 )
  pythia.readString("PDF:pSet = 17");
  pythia.readString("MultipartonInteractions:ecmRef = 200");
  pythia.readString("MultipartonInteractions:bprofile = 2");

  pythia.readString("MultipartonInteractions:pT0Ref = 1.40"); //Gaussian kT term
  pythia.readString("MultipartonInteractions:ecmPow = 0.135");
  pythia.readString("MultipartonInteractions:coreRadius = 0.56");
  pythia.readString("MultipartonInteractions:coreFraction = 0.78");
  pythia.readString("ColourReconnection:range = 5.4");

  //pythia.readString("ColourReconnection:reconnect = on");
  //pythia.readString("ColourReconnection:mode = 1");
  
  // PYTHIA 8 Monash base
//pythia.readString("Tune:pp = 14");
//// QCD-based colour reconnection / CR-BLC Mode 2 style
//pythia.readString("BeamRemnants:remnantMode = 1");
//pythia.readString("ColourReconnection:reconnect = on");
//pythia.readString("ColourReconnection:mode = 1");
//pythia.readString("ColourReconnection:allowDoubleJunRem = off");
//pythia.readString("ColourReconnection:m0 = 0.3");
//pythia.readString("ColourReconnection:allowJunctions = on");
//pythia.readString("ColourReconnection:junctionCorrection = 1.20");
//// This is the important Mode-2 causality/time-dilation part
//pythia.readString("ColourReconnection:timeDilationMode = 2");
//pythia.readString("ColourReconnection:timeDilationPar = 0.18");

  pythia.init();

  // ------------------------------------------------------------
  // Output file
  // ------------------------------------------------------------
  
  TFile* fout = new TFile("tree_out.root", "RECREATE");

  // Inclusive pT histograms
  TH1D* hPtLambda = new TH1D(
      "pt_lambda",
      "#Lambda + #bar{#Lambda};p_{T} [GeV/c];counts",
      200, 0.0, 20.0);

  TH1D* hPtK0S = new TH1D(
      "pt_k0s",
      "K^{0}_{S};p_{T} [GeV/c];counts",
      200, 0.0, 20.0);

  TH1D* hPtLambdaC = new TH1D(
      "pt_lambdac",
      "#Lambda_{c}^{+} + #bar{#Lambda}_{c}^{-};p_{T} [GeV/c];counts",
      200, 0.0, 20.0);

  TH1D* hPtD0 = new TH1D(
      "pt_d0",
      "D^{0} + #bar{D}^{0};p_{T} [GeV/c];counts",
      200, 0.0, 20.0);

  // Optional weighted versions using pythia.info.weight()
  TH1D* hPtLambda_w = new TH1D(
      "pt_lambda_weighted",
      "#Lambda + #bar{#Lambda};p_{T} [GeV/c];weighted counts",
      200, 0.0, 20.0);

  TH1D* hPtK0S_w = new TH1D(
      "pt_k0s_weighted",
      "K^{0}_{S};p_{T} [GeV/c];weighted counts",
      200, 0.0, 20.0);

  TH1D* hPtLambdaC_w = new TH1D(
      "pt_lambdac_weighted",
      "#Lambda_{c}^{+} + #bar{#Lambda}_{c}^{-};p_{T} [GeV/c];weighted counts",
      200, 0.0, 20.0);

  TH1D* hPtD0_w = new TH1D(
      "pt_d0_weighted",
      "D^{0} + #bar{D}^{0};p_{T} [GeV/c];weighted counts",
      200, 0.0, 20.0);

  TH1D* hCounts = new TH1D(
      "counts",
      "Counters;bin;value",
      10, -0.5, 9.5);

  hCounts->GetXaxis()->SetBinLabel(1, "pythia.next ok");
  hCounts->GetXaxis()->SetBinLabel(2, "Lambda");
  hCounts->GetXaxis()->SetBinLabel(3, "K0S");
  hCounts->GetXaxis()->SetBinLabel(4, "Lambda_c");
  hCounts->GetXaxis()->SetBinLabel(5, "D0");

  // ------------------------------------------------------------
  // Event loop
  // ------------------------------------------------------------
  long long nAccepted = 0;
  long long nTried = 0;

  while (nAccepted < targetEvents) {
    ++nTried;

    if (!pythia.next()) {
      continue;
    }

    ++nAccepted;
    hCounts->Fill(0);

    const double wEvt = pythia.info.weight();

    std::vector<std::vector<int>> kids;
    buildChildrenMap(pythia.event, kids);

    for (int i = 0; i < pythia.event.size(); ++i) {
      const int id = pythia.event[i].id();
      const int idAbs = std::abs(id);

      // Fill only last copy to avoid double counting.
      if (hasSameAbsDaughter(pythia.event, kids, i)) {
        continue;
      }

      const double pt = pythia.event[i].pT();

      // Lambda + anti-Lambda
      if (idAbs == 3122) {
        hPtLambda->Fill(pt);
        hPtLambda_w->Fill(pt, wEvt);
        hCounts->Fill(1);
      }

      // K0S
      if (id == 310) {
        hPtK0S->Fill(pt);
        hPtK0S_w->Fill(pt, wEvt);
        hCounts->Fill(2);
      }

      // Lambda_c + anti-Lambda_c
      if (idAbs == 4122) {
        hPtLambdaC->Fill(pt);
        hPtLambdaC_w->Fill(pt, wEvt);
        hCounts->Fill(3);
      }

      // D0 + anti-D0
      if (idAbs == 421) {
        hPtD0->Fill(pt);
        hPtD0_w->Fill(pt, wEvt);
        hCounts->Fill(4);
      }
    }

    if (nAccepted % 100000 == 0) {
      std::cout << "accepted events = " << nAccepted
                << "  tried = " << nTried << "\n";
    }
  }

  // ------------------------------------------------------------
  // Store run info
  // ------------------------------------------------------------
  TParameter<long long>* pNaccepted =
      new TParameter<long long>("Naccepted", nAccepted);

  TParameter<long long>* pNtried =
      new TParameter<long long>("Ntried", nTried);

  TParameter<double>* pSigmaGen =
      new TParameter<double>("sigmaGen_mb", pythia.info.sigmaGen());

  // ------------------------------------------------------------
  // Write output
  // ------------------------------------------------------------

  std::cout << "Integral Lambda   = " << hPtLambda->Integral() << "\n";
  std::cout << "Integral K0S      = " << hPtK0S->Integral() << "\n";
  std::cout << "Integral Lambda_c = " << hPtLambdaC->Integral() << "\n";
  std::cout << "Integral D0       = " << hPtD0->Integral() << "\n";
  fout->cd();

  hPtLambda->Write();
  hPtK0S->Write();
  hPtLambdaC->Write();
  hPtD0->Write();

  hPtLambda_w->Write();
  hPtK0S_w->Write();
  hPtLambdaC_w->Write();
  hPtD0_w->Write();

  hCounts->Write();

  pNaccepted->Write();
  pNtried->Write();
  pSigmaGen->Write();

  fout->Close();

  // ------------------------------------------------------------
  // Print summary
  // ------------------------------------------------------------
  std::cout << "Done.\n";
  std::cout << "accepted events = " << nAccepted << "\n";
  std::cout << "tried events    = " << nTried << "\n";
  std::cout << "sigmaGen(mb)    = " << std::setprecision(8)
            << pythia.info.sigmaGen() << "\n";

  auto t1 = std::chrono::high_resolution_clock::now();
  auto dt = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();

  std::cout << "runtime(s) = " << dt << "\n";

  return 0;
}