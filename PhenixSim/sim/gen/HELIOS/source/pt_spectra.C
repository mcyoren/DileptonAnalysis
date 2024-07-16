#include "TH1.h"

const int nFiles = 9;
TH1D* hist_temp[nFiles] = {NULL};
TH1D* hist[nFiles] = {NULL};

const char* fileNames[nFiles] = {"pi0_gee", "eta_gee", "etap_gee", "omega_ee", "omega_pi0ee", "phi_ee", "phi_etaee", "jpsi","psip"};
const int color[nFiles] = {kRed, kBlack, kBlue, kGreen+2, kYellow+2, kCyan, kMagenta, kOrange+2, kRed+2};


void pt_spectra(){

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();

    for(int ifile= 0; ifile<nFiles; ifile++){

        TFile* fsim = new TFile(Form("/Users/vassudoomra/Desktop/DileptonAnalysis/HELIOS/helios_%s.root",fileNames[ifile]),"READ");
        hist_temp[ifile] = (TH1D*)fsim->Get("h_input_pt");
        
        hist[ifile] = (TH1D*)hist_temp[ifile]->Clone();
        hist[ifile]->Reset("ICESM");
        hist[ifile]->SetName(Form("pt_%s",fileNames[ifile]));

        int nbins = hist[ifile]->GetNbinsX();

        for(int ibin=1; ibin<nbins+1; ibin++){

            double bin_content = hist_temp[ifile]->GetBinContent(ibin)*(1/hist_temp[ifile]->GetBinWidth(ibin))*(1/hist[ifile]->GetBinCenter(ibin));
            double bin_error = hist_temp[ifile]->GetBinError(ibin)*(1/hist_temp[ifile]->GetBinWidth(ibin))*(1/hist[ifile]->GetBinCenter(ibin));

            hist[ifile]->SetBinContent(ibin, bin_content);
            hist[ifile]->SetBinError(ibin, bin_error);


        }

        hist[ifile]->Scale((double)(1./20e6)*(1./2*3.14)*42.);

        if(ifile==0) { hist[ifile]->Draw(); hist[ifile]->SetMarkerColor(color[ifile]); hist[ifile]->SetMarkerStyle(20); hist[ifile]->SetMarkerSize(0.5); }
        else { hist[ifile]->Draw("psame"); hist[ifile]->SetMarkerColor(color[ifile]); hist[ifile]->SetMarkerStyle(20); hist[ifile]->SetMarkerSize(0.5); }

    }


}