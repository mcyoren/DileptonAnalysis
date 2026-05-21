#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

using namespace std;
const int MaxTracks = 20;

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
	const double scale_ptyhia = 1e13; // Scale for Pythia output
    ifstream vfile(filepath);
    if (!vfile.is_open()) {
        cerr << "Error: can't open vertex file " << filepath << endl;
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

    ofstream file(output);
    if (!file.is_open()) {
        cerr << "Error: can't open " << output << endl;
        return;
    }

    TTree* T = dynamic_cast<TTree*>(input->Get("T"));
    if (!T) {
        cerr << "Error: TTree 'T' not found in file " << infile << endl;
        return;
    }

    int ntracks = 0;
    double zvtx = 0;
    double weight = 1.0;
    int pid[MaxTracks] = {0}, isCharm[MaxTracks] = {0}, isBottom[MaxTracks] = {0};
    double px[MaxTracks] = {0}, py[MaxTracks] = {0}, pz[MaxTracks] = {0};
    double energy[MaxTracks] = {0}, vx[MaxTracks] = {0}, vy[MaxTracks] = {0}, vz[MaxTracks] = {0};

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
    T->SetBranchAddress("vy", &vy);
    T->SetBranchAddress("vz", &vz);
    T->SetBranchAddress("weight", &weight);

    // OSCAR header
    file << "# OSC1999A" << endl
         << "# final_id_p_x" << endl
         << "# SimName 1.0" << endl
         << "# " << endl
         << "# Some comments..." << endl << endl;

	Long64_t nEntries = T->GetEntries();
	for (int i = 0; i < 2; ++i)
	{
		for (Long64_t ievt = 0; ievt < nEntries; ++ievt)
		{
			T->GetEntry(ievt);
			file << 0 << "\t" << ntracks << endl;

			size_t iv = ievt * 4;
			for (int i = 0; i < ntracks; ++i)
			{
				double vx_global = vx[i] * scale_ptyhia + ((iv + 0) < vertexes.size() ? vertexes[iv + 0] : 0);
				double vy_global = vy[i] * scale_ptyhia + ((iv + 1) < vertexes.size() ? vertexes[iv + 1] : 0);
				double vz_global = (vz[i] - zvtx) * scale_ptyhia + ((iv + 2) < vertexes.size() ? vertexes[iv + 2] : 0);

                if (isBottom[i]==1) {
                    px[i]=0.1;
                    py[i]=0.1;
				}

				file << i + 1 << "\t"
					 << pid[i] << "\t"
					 << 0 << "\t"
					 << px[i] << "\t"
					 << py[i] << "\t"
					 << pz[i] << "\t"
					 << energy[i] << "\t"
					 << 0.000511 << "\t"
					 << vx_global << "\t"
					 << vy_global << "\t"
					 << vz_global << "\t"
					 << 0 << endl;
			}

			file << 0 << "\t" << 0 << endl;
		}
	}

	file.close();
    input->Close();
    cout << "Written OSCAR file '" << output << "' with " << nEntries << " events." << endl;
}
