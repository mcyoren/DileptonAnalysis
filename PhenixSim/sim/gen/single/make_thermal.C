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
#include <TString.h>
#include <TMath.h>

using namespace std;

void TwoBodyDecay(TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2)
{
  //
  // Isotropic two body decay of parent into two particles
  //
  // input:  parent 4 vector
  // output: daughter1 &2 4 vector
  //
  // Axel Drees 10/17/2018 - checked to work for 2 photon decays
  //            8/29/2021 - adapted from standalone function to this class
  //
  TRandom3 randy = TRandom3(0); // Random Generator
  Double_t pi = 3.14159;        // define pi

  Double_t px, py, pz, E;                // define generic 4 vector component
  Double_t mass = parent.M();            // mass of parent particle
  Double_t m1 = daughter1.M();           // mass of decay particle 1
  Double_t m2 = daughter2.M();           // mass of decay particle 2
  TVector3 boost = parent.BoostVector(); // boost vector of parent particle

  //
  // decay in rest frame
  //
  Double_t p = sqrt((mass * mass - (m1 + m2) * (m1 + m2)) * (mass * mass - (m1 - m2) * (m1 - m2))) / 2 / mass;
  Double_t phi = randy.Uniform(0., 2 * pi);        // phi random from 0. to 2pi
  Double_t z = randy.Uniform(-1., 1.);             // z used to calculate theta
  Double_t theta = acos(z);                        //     so that each angle in 4pi is equaly likely
  pz = p * cos(theta);                             // longitudinal momentum
  px = p * sin(theta) * cos(phi);                  // transverse momentum x component
  py = p * sin(theta) * sin(phi);                  // transverse momentum y component
  E = sqrt(px * px + py * py + pz * pz + m1 * m1); // energy of decay particle 1
  daughter1.SetPxPyPzE(px, py, pz, E);             // set 4 vector of 1 decay particle
  E = sqrt(px * px + py * py + pz * pz + m2 * m2); // energy of decay particle 2
  daughter2.SetPxPyPzE(-px, -py, -pz, E);          // set 4 vector of 2 decay particle in opposite direction

  // boost to parent frame
  daughter1.Boost(boost); // boost with parent momentum
  daughter2.Boost(boost); // boost with parent momentum
}

class TMyRandom
{
public:
  TMyRandom(int seed = 987632621)
  {
    fourVec = new TLorentzVector();
    unif = new TRandom3(seed);
    fhagdorn = new TF1("fhagdorn", "[0]*x/(pow(exp(-[1]*x-[2]*x*x)+x/[3],[4]))", 0.0, 100.0);
    fpow = new TF1("fpow", "pow(x,[0])", 0.0, 100.0);
    pi = 3.14159;
  }

  double Rndm() { return unif->Rndm(); }

  double Gaus(double mean, double sigma) { return unif->Gaus(mean, sigma); }

  double GetPowLaw(double n, double lpt, double upt)
  {
    fpow->SetRange(lpt, upt);
    fpow->SetParameter(0, -n);
    return fpow->GetRandom();
  }

  double GetHagdorn(double lpt, double upt)
  {
    fhagdorn->SetRange(lpt, upt);
    double a = 0.31602, b = 0.124106, c = 31.0164, p0 = 0.805548, n = 8.3051;
    fhagdorn->SetParameters(c, a, b, p0, n);
    return fhagdorn->GetRandom();
  }

  TLorentzVector *GetFMomGaussYPowPT(double ySig, double rapwin, double n, double lpt, double upt, double m0)
  {
    do
      y = unif->Gaus(0, ySig);
    while (fabs(y) > rapwin);
    if (n > 0)
      pt = GetPowLaw(n, lpt, upt);
    else if (n == 0)
      pt = lpt + (upt - lpt) * Rndm();
    else
      pt = GetHagdorn(lpt, upt);
    phi = unif->Rndm() * 2 * pi - pi;
    mt = sqrt(m0 * m0 + pt * pt);
    pl = mt * sinh(y);
    fourVec->SetPxPyPzE(pt * cos(phi), pt * sin(phi), pl, sqrt(mt * mt + pl * pl));
    return fourVec;
  }

private:
  TLorentzVector *fourVec;
  TRandom3 *unif;
  TF1 *fhagdorn;
  TF1 *fpow;
  double phi, y, pt, pl, mt;
  double pi;
};

unsigned int get_seed()
{
  unsigned int fSeed;

  ifstream devrandom;
  devrandom.open("/dev/urandom", ios::binary);
  devrandom.read((char *)&fSeed, sizeof(fSeed));
  devrandom.close();
  if (fSeed != -1)
  {
    cout << " Got seed from /dev/urandom" << endl;
    fSeed = fSeed % 900000000;
  }
  else
    fSeed = 0;
  cout << "seed is " << fSeed << endl;
  return fSeed;
}

void make_thermal(char filepath[200] = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
                  TString fout = "oscar.particles.dat", const int nevt = 10000, const float pt_min = 0.5, const float pt_max = 5,
                  const double n = -1, const int id = 0)
{
  // n: <0 hagdorn (mb HeAu), =0 flat, >0 power law

  // the only key part above is the
  // OSC1999A which specifies the format
  // the rest of the file is listed by
  // events
  // "0 2" (e.g. for two particles)
  // then a list for each particle
  // idpart id ist px,py,pz,E,
  // mass, x,y,z,t
  // and then "0 0"

  // ------- INPUT CARD --------
  double Tpi = 3.14159;
  const int NPARTICLES = 2;
  double scale = 1e13; // cm to fm conversion

  // vertex staff
  std::ifstream myfile(filepath);
  vector<double> vertexes;
  string line;
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      string s;
      std::stringstream ss(line);
      while (getline(ss, s, ' '))
      {
        vertexes.push_back(atof(s.c_str()) * scale);
      }
    }
    myfile.close();
  }

  for (int i = 0; i < (int)vertexes.size() / 4; i++)
  {
    if (false)
      std::cout << vertexes[4 * i] << "     " << vertexes[4 * i + 1] << "     " << vertexes[4 * i + 2] << " " << vertexes[4 * i + 3] << " " << std::endl;
  }

  // ---------------------------

  double rapwidth = 10., rapwin = .5;
  Double_t IDs[2][3] =
      {
          {11, 0, 0.000511},
          {-11, 0, 0.000511}};
  // Double_t IDphoton[3] = {22, 0, 0.};

  int seed = get_seed();
  TMyRandom *myrand = new TMyRandom(seed);
  gRandom->SetSeed(seed);
  cout << "double check seed " << gRandom->Uniform() << endl;

  // Write out the file
  ofstream fileout(fout);
  fileout << "# OSC1999A" << endl;
  fileout << "# final_id_p_x" << endl;
  fileout << "# SimName 1.0" << endl;
  fileout << "#" << endl;
  fileout << "# Some comments..." << endl;

  TF1 *fGaus = new TF1("fGaus", "gaus", -10, 10);
  fGaus->SetParameters(1.27043e+04, 1.77852e+00, 1.27954e+01); // for CuAu

  TH1D h1d_zvtx("h1d_zvtx", "", 400, -20, 20);
  TH1D h1d_pt("h1d_pt", "", 50, 0, 5);
  TH1D h1d_mass("h1d_mass", "", 50, 0, 5);
  TH1D h1d_pt_daughter("h1d_pt_daughter", "", 50, 0, 5);

  TF1 *fThermal = new TF1("mythermalmass", "pow(x,1.5)*exp(-x/[0])", 0.5, 10);
  fThermal->SetParameter(0, 0.3);

  // Get Zvertex dist from data
  //TH1D *zvtx_dat;

  //TFile *f1 = new TFile("heau200_bbcz.root", "READ");
  //zvtx_dat = (TH1D *)f1->Get("h2d_bbcz_px");
  //zvtx_dat->AddDirectory(0);

  for (int ievt = 0; ievt < nevt; ievt++)
  {
    //  Number of particles in event
    fileout << "0 " << NPARTICLES << endl;
    // Try to generate flat in zvtx (-30; 30) and reweight after
    double vertex = vertexes[4 * ievt + 2];
    h1d_zvtx.Fill(vertex / scale);

    //-------------------  PI0 --------------------
    TLorentzVector *thermal = new TLorentzVector();
    thermal = myrand->GetFMomGaussYPowPT(rapwidth, rapwin, n, pt_min, pt_max, fThermal->GetRandom());
    h1d_pt.Fill(thermal->Pt());
    h1d_mass.Fill(thermal->M());

    TLorentzVector *electron = new TLorentzVector();
    TLorentzVector *positron = new TLorentzVector();
    electron->SetPtEtaPhiM(0., 0., 0., IDs[0][2]);
    positron->SetPtEtaPhiM(0., 0., 0., IDs[1][2]);
    TwoBodyDecay(*thermal, *electron, *positron);

    for (int np = 0; np < NPARTICLES; np++)
    {
      TLorentzVector *particle;
      if (np == 0)
        particle = electron;
      else
        particle = positron;

      h1d_pt_daughter.Fill(particle->Pt());
      // Particle index
      Int_t index = np;
      // note that these positions are in femtometers !
      Double_t xpos = vertexes[4 * ievt];
      Double_t ypos = vertexes[4 * ievt + 1];
      Double_t zpos = vertexes[4 * ievt + 2];
      Double_t time = 0.0;
      // idpart id ist px,py,pz,E, mass, x,y,z,t

      fileout << index << " " << IDs[np][0] << " " << IDs[np][1] << " " << particle->Px() << " " << particle->Py() << " " << particle->Pz() << " " << particle->E() << " " << particle->M() << " " << xpos << " " << ypos << " " << zpos << " " << time << endl;
    } // End of particle
    fileout << "0 0" << endl;
  } // End of events

  fileout.close();

  TFile fdiag("diag.root", "recreate");
  h1d_zvtx.Write();
  h1d_pt.Write();
  h1d_pt_daughter.Write();
  h1d_mass.Write();
  fdiag.Write();
  fdiag.Close();
}