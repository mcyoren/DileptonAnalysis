#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TRandom3.h>

using namespace std;
const int MaxTracks = 20;

static int stochastic_round(double x, TRandom3& rng)
{
  // If you want D0-only events to stay 1, enforce minimum 1.
  if (x <= 1.0) return 1;

  int n0 = (int)std::floor(x);
  double frac = x - (double)n0;
  int n = n0 + (rng.Uniform() < frac ? 1 : 0);
  if (n < 1) n = 1;
  return n;
}

void WriteROOT2OscarPythia(
    const TString filepath = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt",
    const TString infile   = "single_pi0_HELIOS_1B.root",
    const TString output   = "oscar.txt"
) {
    TFile* input = TFile::Open(infile, "READ");
    if (!input || input->IsZombie()) {
        cerr << "Error: can't open " << infile << endl;
        return;
    }

    // Vertex file reading
    const double scale = 1e13;
    const double scale_ptyhia = 1e13;
    ifstream vfile(filepath.Data());
    if (!vfile.is_open()) {
        cerr << "Error: can't open vertex file " << filepath << endl;
        input->Close();
        return;
    }

    vector<double> vertexes;
    string line;
    while (getline(vfile, line)) {
        stringstream ss(line);
        double val;
        while (ss >> val) vertexes.push_back(val * scale);
    }
    vfile.close();

    ofstream file(output.Data());
    if (!file.is_open()) {
        cerr << "Error: can't open " << output << endl;
        input->Close();
        return;
    }

    TTree* T = dynamic_cast<TTree*>(input->Get("T"));
    if (!T) {
        cerr << "Error: TTree 'T' not found in file " << infile << endl;
        file.close();
        input->Close();
        return;
    }

    int ntracks = 0;
    double zvtx = 0;

    int pid[MaxTracks] = {0}, isCharm[MaxTracks] = {0}, isBottom[MaxTracks] = {0};
    double px[MaxTracks] = {0}, py[MaxTracks] = {0}, pz[MaxTracks] = {0};
    double energy[MaxTracks] = {0}, vx[MaxTracks] = {0}, vy[MaxTracks] = {0}, vz[MaxTracks] = {0};

    // IMPORTANT: weight as per-track BR factor
    double wgt[MaxTracks] = {1.0};

    T->SetBranchAddress("ntracks", &ntracks);
    T->SetBranchAddress("zvtx", &zvtx);
    T->SetBranchAddress("pid", pid);
    T->SetBranchAddress("isCharm", isCharm);
    T->SetBranchAddress("isBottom", isBottom);
    T->SetBranchAddress("px", px);
    T->SetBranchAddress("py", py);
    T->SetBranchAddress("pz", pz);
    T->SetBranchAddress("energy", energy);
    T->SetBranchAddress("vx", vx);
    T->SetBranchAddress("vy", vy);   // FIX: was &vy
    T->SetBranchAddress("vz", vz);   // FIX: was &vz
    T->SetBranchAddress("weight", wgt); // CHANGED: array, not scalar

    // OSCAR header
    file << "# OSC1999A" << endl
         << "# final_id_p_x" << endl
         << "# SimName 1.0" << endl
         << "# " << endl
         << "# Some comments..." << endl << endl;

    TRandom3 rng(12345);

    Long64_t nEntries = T->GetEntries();
    int nev_written = 0;
    for (Long64_t ievt = 0; ievt < nEntries; ++ievt)
    {
        T->GetEntry(ievt);

        // Product of per-electron weights (BR factors)
        double wprod = 1.0;
        for (int i = 0; i < ntracks; ++i) {
            wprod *= (wgt[i] > 0.0 ? wgt[i] / 0.068 : 1.0);
        }

        // Number of times to write this event (rounded in expectation)
        int nRep = stochastic_round(wprod, rng);

        for (int irep = 0; irep < nRep; ++irep)
        {
            file << 0 << "\t" << ntracks << endl;

            size_t iv = (size_t)nev_written * 4;
            for (int i = 0; i < ntracks; ++i)
            {
                double vx_global = vx[i] * scale_ptyhia + ((iv + 0) < vertexes.size() ? vertexes[iv + 0] : 0.0);
                double vy_global = vy[i] * scale_ptyhia + ((iv + 1) < vertexes.size() ? vertexes[iv + 1] : 0.0);
                double vz_global = (vz[i] - zvtx) * scale_ptyhia + ((iv + 2) < vertexes.size() ? vertexes[iv + 2] : 0.0);

                // keep your bottom hack
                double px_out = px[i];
                double py_out = py[i];
                if (isBottom[i] == 1) {
                    px_out = 0.1;
                    py_out = 0.1;
                }

                file << i + 1 << "\t"
                     << pid[i] << "\t"
                     << 0 << "\t"
                     << px_out << "\t"
                     << py_out << "\t"
                     << pz[i] << "\t"
                     << energy[i] << "\t"
                     << 0.000511 << "\t"
                     << vx_global << "\t"
                     << vy_global << "\t"
                     << vz_global << "\t"
                     << 0 << endl;
            }

            file << 0 << "\t" << 0 << endl;
            nev_written++;
            //stopping as soon as 10000 events are written
            if(nev_written >= 10000) 
            {
                std::cout<<"\033[1;31m"<<"Reached 10000 events, stopping writing more at "<<ievt<<"\033[0m"<<std::endl;
                break;
            }
        }
        //stopping as soon as 10000 events are written
        if(nev_written >= 10000)
        {
            std::cout<<"\033[1;31m"<<"Reached 10000 events, stopping writing more at "<<ievt<<"\033[0m"<<std::endl;
            break;
        }
    }

    file.close();
    input->Close();
    cout << "Written OSCAR file '" << output << "' with BR-replicated events from "
         << nEntries << " entries." << endl;
}
