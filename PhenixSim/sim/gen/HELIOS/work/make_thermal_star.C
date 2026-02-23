// make_thermal_to_tree.C
// ROOT macro-style code that:
//  1) reads vertexes from a text file (x y z t per line, in cm -> converts to fm)
//  2) generates a thermal "parent" 4-vector in rapidity range y in [-1,1]
//  3) decays parent -> e+ e- isotropically with your TwoBodyDecay
//  4) writes output to a ROOT TTree in the same “WriteEvent/WriteTrack” style as your HELIOS example
//
// NOTE (important bugfix): your original TwoBodyDecay created TRandom3(0) *inside* the function,
// so every call repeats the same random numbers. Below we make the RNG static (one RNG for all decays).

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>

using namespace std;

// If you want the same WriteEvent/WriteTrack format as in your example:
#include "../source/HELIOSLibrary/HELIOSLibrary.h"

static const double kPi = 3.14159265358979323846;

// ---------------------------------------------------------------------
// Isotropic two body decay with a persistent RNG (bugfix vs your version)
// ---------------------------------------------------------------------
void TwoBodyDecay(const TLorentzVector& parent,
                  TLorentzVector& daughter1,
                  TLorentzVector& daughter2)
{
  static TRandom3 randy(0); // persistent RNG, seeded once

  const double mass = parent.M();
  const double m1   = daughter1.M();
  const double m2   = daughter2.M();

  const TVector3 boost = parent.BoostVector();

  // p* in rest frame
  const double term1 = (mass*mass - (m1+m2)*(m1+m2));
  const double term2 = (mass*mass - (m1-m2)*(m1-m2));
  const double p = (term1 > 0 && term2 > 0) ? std::sqrt(term1*term2)/(2.0*mass) : 0.0;

  const double phi   = randy.Uniform(0., 2.0*kPi);
  const double z     = randy.Uniform(-1., 1.);
  const double theta = std::acos(z);

  const double pz = p * std::cos(theta);
  const double px = p * std::sin(theta) * std::cos(phi);
  const double py = p * std::sin(theta) * std::sin(phi);

  double E1 = std::sqrt(px*px + py*py + pz*pz + m1*m1);
  double E2 = std::sqrt(px*px + py*py + pz*pz + m2*m2);

  daughter1.SetPxPyPzE( px,  py,  pz, E1);
  daughter2.SetPxPyPzE(-px, -py, -pz, E2);

  daughter1.Boost(boost);
  daughter2.Boost(boost);
}

// ---------------------------------------------------------------------
// Your helper RNG class (kept mostly as-is, but we add a "flat rapidity" generator)
// ---------------------------------------------------------------------
class TMyRandom
{
public:
  TMyRandom(int seed = 987632621)
  {
    fourVec = new TLorentzVector();
    unif = new TRandom3(seed);
    fhagdorn = new TF1("fhagdorn", "[0]*x/(pow(exp(-[1]*x-[2]*x*x)+x/[3],[4]))", 0.0, 100.0);
    fpow = new TF1("fpow", "exp(-x/[0])", 0.0, 100.0);
  }

  ~TMyRandom()
  {
    delete fourVec;
    delete unif;
    delete fhagdorn;
    delete fpow;
  }

  double Rndm() { return unif->Rndm(); }

  double GetPowLaw(double n, double lpt, double upt)
  {
    fpow->SetRange(lpt, upt);
    fpow->SetParameter(0, n);
    return fpow->GetRandom();
  }

  double GetHagdorn(double lpt, double upt)
  {
    fhagdorn->SetRange(lpt, upt);
    double a = 0.31602, b = 0.124106, c = 31.0164, p0 = 0.805548, n = 8.3051;
    fhagdorn->SetParameters(c, a, b, p0, n);
    return fhagdorn->GetRandom();
  }

  // NEW: flat rapidity in [y_low, y_high]
  // pt shape choice:
  //   if shape > 0   -> "power law" in mT with scale=shape (matches your GetPowLaw usage)
  //   if shape == 0  -> flat pt
  //   if shape < 0   -> Hagedorn
  //
  // NOTE: This preserves your original convention, even though the naming is confusing.
  TLorentzVector* GetFMomFlatYShapePT(double y_low, double y_high,
                                      double shape,
                                      double lpt, double upt,
                                      double m0)
  {
    const double y = unif->Uniform(y_low, y_high);

    double pt = 0.0;
    if (shape > 0)
    {
      // sample mt with exp(-mt/shape) then convert to pt
      const double lmt = std::sqrt(m0*m0 + lpt*lpt);
      const double umt = std::sqrt(m0*m0 + upt*upt);
      const double mt  = GetPowLaw(shape, lmt, umt);
      pt = (mt > m0) ? std::sqrt(mt*mt - m0*m0) : 0.0;
    }
    else if (shape == 0)
    {
      pt = lpt + (upt - lpt) * unif->Rndm();
    }
    else
    {
      pt = GetHagdorn(lpt, upt);
    }

    const double phi = unif->Rndm() * 2.0*kPi - kPi;
    const double mt  = std::sqrt(m0*m0 + pt*pt);
    const double pl  = mt * std::sinh(y);
    const double E   = std::sqrt(mt*mt + pl*pl);

    fourVec->SetPxPyPzE(pt*std::cos(phi), pt*std::sin(phi), pl, E);
    return fourVec;
  }

private:
  TLorentzVector *fourVec = nullptr;
  TRandom3 *unif = nullptr;
  TF1 *fhagdorn = nullptr;
  TF1 *fpow = nullptr;
};

// ---------------------------------------------------------------------
// Seed helper
// ---------------------------------------------------------------------
unsigned int get_seed()
{
  unsigned int fSeed = 0;
  ifstream devrandom("/dev/urandom", ios::binary);
  if (devrandom.good())
  {
    devrandom.read((char *)&fSeed, sizeof(fSeed));
    devrandom.close();
    fSeed = fSeed % 900000000u;
    cout << "Got seed from /dev/urandom: " << fSeed << endl;
  }
  else
  {
    cout << "WARNING: /dev/urandom not available, using seed=0" << endl;
    fSeed = 0;
  }
  return fSeed;
}

// ---------------------------------------------------------------------
// HELIOS-style helper (copied from your example)
// ---------------------------------------------------------------------
void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist,
              float px, float py, float pz, float en, float mass,
              float xpos, float ypos, float zpos,
              double br, int parent_index, long double weight)
{
  newTrack.SetFinal(isFinal);
  newTrack.SetNum(num);
  newTrack.SetID(id);
  newTrack.SetIst(ist);
  newTrack.SetPx(px);
  newTrack.SetPy(py);
  newTrack.SetPz(pz);
  newTrack.SetEnergy(en);
  newTrack.SetMass(mass);
  newTrack.SetXpos(xpos);
  newTrack.SetYpos(ypos);
  newTrack.SetZpos(zpos);
  newTrack.SetBranch(br);
  newTrack.SetParentIndex(parent_index);
  newTrack.SetWeight(weight);
}

// ---------------------------------------------------------------------
// Main function: thermal production -> TTree (WriteEvent style)
// ---------------------------------------------------------------------
void make_thermal_to_tree(const char* vertex_filepath =
                            "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
                          TString outroot = "thermal.root",
                          const int nevt = 1000000,
                          const float pt_min = 0.,
                          const float pt_max = 10.0,
                          const float y_min  = -1.0,
                          const float y_max  =  1.0,
                          const float temp   = 0.3) // GeV
{
  // IDs for e- and e+
  const int    pid_electron  =  11;
  const int    pid_positron  = -11;
  const double me = 0.000511;

  // Vertex units: your file is in cm, OSCAR wanted fm. Keep your scale.
  const double scale = 1e13; // cm -> fm

  // Read vertexes (expects 4 numbers per event: x y z t, or at least x y z)
  vector<double> vertexes;
  {
    ifstream myfile(vertex_filepath);
    string line;
    if (!myfile.is_open())
    {
      cout << "ERROR: cannot open vertex file: " << vertex_filepath << endl;
      return;
    }
    while (getline(myfile, line))
    {
      if (line.empty()) continue;
      string s;
      stringstream ss(line);
      while (getline(ss, s, ' '))
      {
        if (s.empty()) continue;
        vertexes.push_back(atof(s.c_str()) * scale); // store in fm
      }
    }
    myfile.close();
  }

  const int nVtxEvents = (int)vertexes.size() / 4;
  if (nVtxEvents <= 0)
  {
    cout << "ERROR: vertex file did not provide 4 columns per event." << endl;
    return;
  }
  if (nevt > nVtxEvents)
  {
    cout << "WARNING: nevt=" << nevt << " > vertex events=" << nVtxEvents
         << ". Will loop in cycle" << endl;
  }

  // RNG
  const unsigned int seed = get_seed();
  TMyRandom myrand(seed);
  gRandom->SetSeed(seed);

  // Thermal mass distribution (same spirit as your original code)
  TF1 fThermalMass("fThermalMass", "pow(x,1.5)*exp(-x/[0])", 0.0, 10.0);
  fThermalMass.SetParameter(0, temp);

  // Output ROOT file + TTree (HELIOS style)
  TFile* f = new TFile(outroot, "RECREATE");

  WriteEvent MyEvent;
  WriteTrack MyTrack;

  TTree* T = new TTree("T","Thermal parent -> e+e- (WriteEvent)");
  T->Branch("MyEvent", &MyEvent);

  // Diagnostics (optional)
  TH1D h_zvtx("h_zvtx","zvtx [cm]", 400, -20, 20);
  TH1D h_parent_pt("h_parent_pt","parent pT [GeV]", 100, 0, 10);
  TH1D h_parent_y("h_parent_y","parent y", 120, -3, 3);
  TH1D h_parent_m("h_parent_m","parent mass [GeV]", 120, 0, 6);
  TH1D h_dau_pt("h_dau_pt","daughter pT [GeV]", 100, 0, 10);

  // Event loop
  const int nloop = nevt; // you can set this to a smaller number for testing

  for (int ievt = 0; ievt < nloop; ievt++)
  {
    if (ievt % 100000 == 0) cout << "event " << ievt << "/" << nloop << endl;

    MyEvent.ClearEvent();
    int nentry = 0;

    // Vertex in fm
    const double vx = vertexes[4*(ievt%nVtxEvents) + 0];
    const double vy = vertexes[4*(ievt%nVtxEvents) + 1];
    const double vz = vertexes[4*(ievt%nVtxEvents) + 2];
    // const double vt = vertexes[4*(ievt%nVtxEvents) + 3]; // unused here

    h_zvtx.Fill(vz/scale);

    // -------------------- "thermal parent" --------------------
    const double m0 = fThermalMass.GetRandom(); // sample mass
    TLorentzVector parent;
    {
      // We generate parent directly with flat rapidity y in [-1,1]
      // and a simple thermal-ish mt exponential using "shape=temp".
      TLorentzVector* tmp = myrand.GetFMomFlatYShapePT(y_min, y_max,
                                                      temp,        // shape parameter
                                                      pt_min, pt_max,
                                                      m0);
      parent = *tmp;
    }

    h_parent_pt.Fill(parent.Pt());
    h_parent_y.Fill(parent.Rapidity());
    h_parent_m.Fill(parent.M());

    // Store parent as non-final entry (like your HELIOS example)
    // br/branch: use 0 to mean "2-body ee"
    const int mother_index = nentry;
    setTrack(MyTrack,
             /*isFinal=*/0,
             /*num=*/nentry,
             /*id=*/0,          // if you want a PDG id here, set it (e.g. 999999) or a real meson id
             /*ist=*/0,
             parent.Px(), parent.Py(), parent.Pz(), parent.E(), parent.M(),
             vx, vy, vz,
             /*br=*/0.0,
             /*parent_index=*/-999,
             /*weight=*/1.0L);
    MyEvent.AddEntry(MyTrack);
    nentry++;

    // -------------------- decay to e+ e- --------------------
    TLorentzVector eminus, eplus;
    eminus.SetPtEtaPhiM(0., 0., 0., me);
    eplus .SetPtEtaPhiM(0., 0., 0., me);

    TwoBodyDecay(parent, eminus, eplus);

    // daughter 0: e-
    h_dau_pt.Fill(eminus.Pt());
    setTrack(MyTrack,
             /*isFinal=*/1,
             /*num=*/nentry,
             /*id=*/pid_electron,
             /*ist=*/0,
             eminus.Px(), eminus.Py(), eminus.Pz(), eminus.E(), eminus.M(),
             vx, vy, vz,
             /*br=*/1.0,
             /*parent_index=*/mother_index,
             /*weight=*/1.0L);
    MyEvent.AddEntry(MyTrack);
    nentry++;

    // daughter 1: e+
    h_dau_pt.Fill(eplus.Pt());
    setTrack(MyTrack,
             /*isFinal=*/1,
             /*num=*/nentry,
             /*id=*/pid_positron,
             /*ist=*/0,
             eplus.Px(), eplus.Py(), eplus.Pz(), eplus.E(), eplus.M(),
             vx, vy, vz,
             /*br=*/1.0,
             /*parent_index=*/mother_index,
             /*weight=*/1.0L);
    MyEvent.AddEntry(MyTrack);
    nentry++;

    MyEvent.SetNStable(2); // two final daughters
    T->Fill();

    MyEvent.ClearEvent();
  }

  f->cd();
  T->Write();
  h_zvtx.Write();
  h_parent_pt.Write();
  h_parent_y.Write();
  h_parent_m.Write();
  h_dau_pt.Write();
  f->Close();

  cout << "Wrote: " << outroot << " with " << nloop << " events" << endl;
}