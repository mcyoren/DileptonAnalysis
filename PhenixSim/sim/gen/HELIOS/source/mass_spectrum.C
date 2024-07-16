#include "TFile.h"
#include "TTree.h"
#include "HELIOSLibrary/HELIOSLibrary.h"

const int nevents = 20000000;
const int centbin = 1;

TH2D* h2d_ee_FG[centbin] = {NULL};
TH2D* h2d_ee_FG_like[centbin] = {NULL};
TH1D* h1d_ee_FG_temp[centbin] = {NULL};
TH1D* h1d_ee_FG[centbin] = {NULL};

TH1D* h1d_ee_FG_like_temp[centbin] = {NULL};
TH1D* h1d_ee_FG_like[centbin] = {NULL};

static const double Me = 0.000510998918;
static const double Me2 = Me*Me;
static const double c = 299792458;

void init_hists(){

  for(int icent = 0; icent < centbin; icent++){

        h2d_ee_FG[icent] = new TH2D(Form("h2d_ee_FG_cent%d", icent), "ee mass distribution FG", 450, 0, 4.5 , 300, 0, 10);
        h2d_ee_FG[icent]->Sumw2();

        h1d_ee_FG[icent] = new TH1D(Form("h1d_ee_FG_cent%d", icent), "ee mass distribution FG", 450, 0, 4.5);
        h1d_ee_FG[icent]->Sumw2();

        h2d_ee_FG_like[icent] = new TH2D(Form("h2d_ee_FG_like_cent%d", icent), "ee mass distribution Like FG", 450, 0, 4.5 , 300, 0, 10);
        h2d_ee_FG_like[icent]->Sumw2();

        h1d_ee_FG_like[icent] = new TH1D(Form("h1d_ee_FG_like_cent%d", icent), "ee mass distribution Like FG", 450, 0, 4.5);
        h1d_ee_FG_like[icent]->Sumw2();

    }
}

double getPT(double px1, double py1, double pz1, double px2, double py2, double pz2)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  pair = p1 + p2;

  return pair.Pt();
}

double getMass(double px1, double py1, double pz1, double px2, double py2, double pz2)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  pair = p1 + p2;

  return pair.M();
}


void mass_spectrum(){

    init_hists();

    TFile* fjpsi = new TFile("helios_jpsi.root", "READ");
    TTree* tjpsi = (TTree*)fjpsi->Get("T");
    WriteEvent* jpsi_event = 0;
    TBranch* br_jpsi  = tjpsi->GetBranch("MyEvent");
    br_jpsi->SetAddress(&jpsi_event);

    for(int ievt=0; ievt<nevents; ievt++){

        if(ievt%1000000==0) cout << "At Tree Entry = " << ievt << endl;

        br_jpsi->GetEntry(ievt);
        int nstable = jpsi_event->GetNStable();

        for(int i = 0; i < nstable; i++){

            WriteTrack trk_A1 = jpsi_event->GetWriteTrack(i+1);
            double px_A1 = trk_A1.GetPx();
            double py_A1 = trk_A1.GetPy();
            double pz_A1 = trk_A1.GetPz();
            double ecore_A1 = sqrt(px_A1*px_A1 + py_A1*py_A1 + pz_A1*pz_A1 + Me2);
            int charge_A1 = -999;
            if(trk_A1.GetID() == 11) charge_A1 = -1;
            else charge_A1 = 1;

            for(int j=i+1; j<nstable; j++){

                WriteTrack trk_A2 = jpsi_event->GetWriteTrack(j+1);
                double px_A2 = trk_A2.GetPx();
                double py_A2 = trk_A2.GetPy();
                double pz_A2 = trk_A2.GetPz();
                double ecore_A2 = sqrt(px_A2*px_A2 + py_A2*py_A2 + pz_A2*pz_A2 + Me2);
                int charge_A2 = -999;
                if(trk_A2.GetID() == 11) charge_A2 = -1;
                else charge_A2 = 1;

                //if(!(ecore_A1 > 1.4 || ecore_A2 > 1.4)) continue;

                double mass_ee_FG  = getMass(px_A1, py_A1, pz_A1, px_A2, py_A2, pz_A2);
                double pt_ee_FG = getPT(px_A1, py_A1, pz_A1, px_A2, py_A2, pz_A2);

                double wt = 1.;

                if (charge_A1 != charge_A2) h2d_ee_FG[0]->Fill(mass_ee_FG, pt_ee_FG, wt);
                else h2d_ee_FG_like[0]->Fill(mass_ee_FG, pt_ee_FG, wt);

            }

        }

    }
        
    for(int icent=0; icent<centbin; icent++){

        h2d_ee_FG[icent]->GetYaxis()->SetRangeUser(0, 10);
        h1d_ee_FG_temp[icent] = (TH1D*)h2d_ee_FG[icent]->ProjectionX();

        int nbinsX = h1d_ee_FG_temp[icent]->GetNbinsX();

        for(int ibin=1; ibin < nbinsX+1; ibin++){

            double binwidth = h1d_ee_FG_temp[icent]->GetBinWidth(ibin);
            double content = h1d_ee_FG_temp[icent]->GetBinContent(ibin);
            double err = h1d_ee_FG_temp[icent]->GetBinError(ibin);
            h1d_ee_FG[icent]->SetBinContent(ibin, content/binwidth);
            h1d_ee_FG[icent]->SetBinError(ibin, err/binwidth);

        }
    }

    for(int icent=0; icent<centbin; icent++){

        h2d_ee_FG_like[icent]->GetYaxis()->SetRangeUser(0, 10);
        h1d_ee_FG_like_temp[icent] = (TH1D*)h2d_ee_FG_like[icent]->ProjectionX();

        int nbinsX = h1d_ee_FG_like_temp[icent]->GetNbinsX();

        for(int ibin=1; ibin < nbinsX+1; ibin++){

            double binwidth = h1d_ee_FG_like_temp[icent]->GetBinWidth(ibin);
            double content = h1d_ee_FG_like_temp[icent]->GetBinContent(ibin);
            double err = h1d_ee_FG_like_temp[icent]->GetBinError(ibin);
            h1d_ee_FG_like[icent]->SetBinContent(ibin, content/binwidth);
            h1d_ee_FG_like[icent]->SetBinError(ibin, err/binwidth);

        }
    }

    cout<<"Writing TTrees to file"<<endl;

    TFile* fout = new TFile("helios_mass_spectrum_jpsi.root","RECREATE");
    fout->cd();

    for(int icent =0; icent < centbin; icent++){

        h2d_ee_FG[icent]->Write();
        h1d_ee_FG[icent]->Write();
        h2d_ee_FG_like[icent]->Write();
        h1d_ee_FG_like[icent]->Write();
    }

    fout->Close();



}