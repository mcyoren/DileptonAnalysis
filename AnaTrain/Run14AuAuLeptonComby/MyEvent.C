#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"

#include "MyEvent.h"

ClassImp(MyDileptonAnalysis::MyEvent)
ClassImp(MyDileptonAnalysis::MyTrack)
ClassImp(MyDileptonAnalysis::MyElectron)
ClassImp(MyDileptonAnalysis::MyHadron)
ClassImp(MyDileptonAnalysis::MyVTXHit)

using namespace MyDileptonAnalysis;
using namespace std;

namespace MyDileptonAnalysis
{
    void MyTrack::SetPrimes(const float bbcz, const float svxz, const int rg_beamoffset, const int rungroup)
    {
        float alpha_offset = -999;
        float phi_offset = -999;

        if (this->GetArm() == 0)
        {
            alpha_offset = (fXoffset[0][rg_beamoffset] / 220) * TMath::Sin(this->GetPhiDC()) + (fYoffset[0][rg_beamoffset] / 220) * TMath::Cos(this->GetPhiDC());
            phi_offset = dilep_par0_phi[0] * TMath::Sin(this->GetPhi0()) + dilep_par1_phi[0] * TMath::Cos(this->GetPhi0()) + dilep_par2_phi[0];
        }
        else
        {
            alpha_offset = (fXoffset[1][rg_beamoffset] / 220) * TMath::Sin(this->GetPhiDC()) + (fYoffset[1][rg_beamoffset] / 220) * TMath::Cos(this->GetPhiDC());
            phi_offset = dilep_par0_phi[1] * TMath::Sin(this->GetPhi0()) + dilep_par1_phi[1] * TMath::Cos(this->GetPhi0()) + dilep_par2_phi[1];
        }
        this->SetAlphaPrime(this->GetAlpha() - alpha_offset);

        if (this->GetArm() == 0)
            this->SetPhi0Prime(this->GetPhi0() - phi_offset + res_rot_east);
        else
            this->SetPhi0Prime(this->GetPhi0() - phi_offset + res_rot_west);
        // set Phi0 to right value
        this->SetPhi0(this->GetPhi0() - 2.0195 * alpha_offset);

        this->SetPtPrime(this->GetPt() * fabs(this->GetAlpha() / this->GetAlphaPrime()) * mscale);

        if (this->GetAlpha() * this->GetAlphaPrime() < 0)
            this->SetQPrime(-this->GetCharge());
        else
            this->SetQPrime(this->GetCharge());

        const int DCArm = this->GetArm();

        const float new_the0 = this->GetThe0() - ((bbcz - svxz) / 220) * TMath::Sin(this->GetThe0());
        const float theta_offset = dilep_par0_theta[DCArm][rungroup] * TMath::Sin(new_the0) + dilep_par1_theta[DCArm][rungroup] * TMath::Cos(new_the0) + dilep_par2_theta[DCArm][rungroup];
        this->SetThe0Prime(new_the0 - theta_offset);
    }

    void MyVTXHit::SetPolars(const float xvtx, const float yvtx, const float zvtx)
    {
        rhit = sqrt(xhit * xhit + yhit * yhit);
        TVector3 hitpoint;
        hitpoint.SetXYZ(xhit - xvtx, yhit - yvtx, zhit - zvtx);
        phihit = hitpoint.Phi();
        if (phihit < -TMath::Pi() / 2)
            phihit += 2 * TMath::Pi();
        thetahit = hitpoint.Theta();
    }

    void MyVTXHit::SetiLayerFromR()
    {
        if (layer == 0 || layer == 1)
            ilayer = layer;
        else if (rhit > 10.20 && rhit < 10.70)
            ilayer = 2;
        else if (rhit > 11.60 && rhit < 12.00)
            ilayer = 3;
        else if (rhit > 12.70 && rhit < 13.40)
            ilayer = 4;
        else if (rhit > 15.20 && rhit < 15.80)
            ilayer = 5;
        else if (rhit > 16.50 && rhit < 17.00)
            ilayer = 6;
        else if (rhit > 17.60 && rhit < 18.20)
            ilayer = 7;
        else
            ilayer = -7777;
    }

    void MyEvent::Associate_Hits_to_Leptons(TH2D *hist_br, TH2D *hist_bz, bool test, int is_fill_hsits, TH3D *dphi_hist_el[N_centr*2],
                                            TH3D *sdphi_hist_el[N_centr*2], TH3D *dthe_hist_el[N_centr*2], TH3D *sdthe_hist_el[N_centr*2], TH3D *chi2_ndf[N_centr])
    {
        const int nleptons = this->GetNtrack();
        const int nvtxhits = this->GetNVTXhit();
        const int centrality = this->GetCentrality();
        const int rungroup = this->GetRunNumber();
        if (test)
            is_fill_hsits = 0;
        const int central_bin = (int)centrality / 20;
        if (central_bin > 4 || central_bin < 0)
            return;
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = this->GetEntry(itrk);
            // mytrk = static_cast<MyDileptonAnalysis::MyHadron*>(myhad);
            const float pt = mytrk->GetPtPrime();

            const float thetaprime = mytrk->GetThe0Prime();

            float phi0_trk_proj = mytrk->GetPhi0Prime();
            float the0_trk_proj = mytrk->GetThe0Prime();
            const float pz = mytrk->GetPtPrime() * (TMath::Cos(thetaprime)) / (TMath::Sin(thetaprime));

            float rp = sqrt(this->GetPreciseX() * this->GetPreciseX() + this->GetPreciseY() * this->GetPreciseY());
            float zp = this->GetPreciseZ();

            float dilep_phi_projection[total_vtx_layers];
            float dilep_the_projection[total_vtx_layers];
            for (int ii = 0; ii < total_vtx_layers; ii++)
            {
                dilep_phi_projection[ii] = -999;
                dilep_the_projection[ii] = -999;
            }
            for (int p = 1; p < N_steps; p++)
            {
                for (int l = 0; l < total_vtx_layers; l++)
                {
                    if (fabs(rp - radii[l]) < step_size && dilep_phi_projection[l] < -900)
                    {
                        dilep_phi_projection[l] = phi0_trk_proj;
                        dilep_the_projection[l] = the0_trk_proj;
                    }
                }

                const int rbin = hist_bz->GetXaxis()->FindBin(rp);
                const int zbin = hist_bz->GetYaxis()->FindBin(zp);

                const float bz = hist_bz->GetBinContent(rbin, zbin) / 10000; // Conversion from Gaus to Tesla

                const float delta_phi0 = (mytrk->GetChargePrime() * 0.3 * step_size * bz) / (2 * mytrk->GetPtPrime() * 100);
                phi0_trk_proj += delta_phi0;

                const float bradial = hist_br->GetBinContent(rbin, zbin) / 10000; // Conversion from Gaus to Tesla
                // Bend in the z direction does not depend upon the charge.

                const float delta_the0 = 0.3 * bradial * (step_size * TMath::Tan(pi / 2 - the0_trk_proj)) / (2 * pz * 100);

                if (thetaprime > pi / 2)
                    the0_trk_proj -= delta_the0;
                else
                    the0_trk_proj += delta_the0;

                zp += step_size * TMath::Tan(pi / 2 - the0_trk_proj);
                rp += step_size;
            }
            const unsigned int charge_bin = (1 - mytrk->GetChargePrime()) / 2;
            const float quality = mytrk->GetTrkQuality();
            if (quality == 63 || quality == 31 || quality == 51)
                mytrk->SetisERT(1);
            else
                mytrk->SetisERT(0);

            float min[nvtx_layers] = {100, 100, 100, 100};
            
            std::vector<int> numbers[4];

            for (int iter_layer = 3; iter_layer >= 0; iter_layer--)
            {
                int n_iter_hits = 1;
                if(iter_layer==2) n_iter_hits += mytrk->GetHitCounter(iter_layer+1);
                if(iter_layer==1 && mytrk->GetHitCounter(2)+mytrk->GetHitCounter(3)==0) break;
                if(iter_layer==1) n_iter_hits = mytrk->GetHitCounter(2)+mytrk->GetHitCounter(3);/////crunh
                if(iter_layer==0 && mytrk->GetHitCounter(1)==0) break;
                if(iter_layer==0) n_iter_hits = mytrk->GetHitCounter(1);
                for (int iassociatedhit = 0; iassociatedhit < n_iter_hits; iassociatedhit++)
                {
                    float dphi_previous_layer = 0;
                    float dthe_previous_layer = 0;
                    if(iassociatedhit>0&&iter_layer>1)
                    {
                        dphi_previous_layer = mytrk->GetdPhi(iter_layer+1, iassociatedhit-1);
                        dthe_previous_layer = mytrk->GetdThe(iter_layer+1, iassociatedhit-1);
                        numbers[3].push_back(iassociatedhit);
                    }
                    if(iter_layer<2 && !(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) )
                    {
                        dphi_previous_layer = mytrk->GetdPhi(iter_layer+1, iassociatedhit);
                        dthe_previous_layer = mytrk->GetdThe(iter_layer+1, iassociatedhit);
                    }
                    if(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2))
                    {
                        dphi_previous_layer = mytrk->GetdPhi(iter_layer+2, iassociatedhit-mytrk->GetHitCounter(2));
                        dthe_previous_layer = mytrk->GetdThe(iter_layer+2, iassociatedhit-mytrk->GetHitCounter(2));
                    }

                    int iter_num = 0;
                    for (int ihit = 0; ihit < nvtxhits; ihit++)
                    {
                        MyDileptonAnalysis::MyVTXHit *vtxhit = this->GetVTXHitEntry(ihit);

                        const int layer = vtxhit->GetLayer();

                        if (layer != iter_layer)
                            continue;

                        const int ilayer = vtxhit->GetiLayer();
                        if (ilayer < 0)
                            continue;

                        const float phi_hit = vtxhit->GetPhiHit();
                        const float theta_hit = vtxhit->GetTheHit();

                        float sigma_phi_value = mytrk->get_sigma_phi_data(rungroup, central_bin, layer);
                        float mean_phi_value = mytrk->get_mean_phi_data(rungroup, central_bin, layer);
                        float sigma_theta_value = mytrk->get_sigma_theta_data(rungroup, central_bin, layer);
                        float mean_theta_value = mytrk->get_mean_theta_data(rungroup, central_bin, layer);

                        if(iter_layer<2||iassociatedhit>0)
                        {
                            sigma_phi_value /= 2;
                            sigma_theta_value /= 1+iter_layer/2;
                            mean_phi_value = dphi_previous_layer*1.15;
                            mean_theta_value += dthe_previous_layer;
                        }

                        const float dphi = (dilep_phi_projection[ilayer] - phi_hit);
                        const float dthe = (dilep_the_projection[ilayer] - theta_hit);
                        const float sdphi = (dphi - mean_phi_value) / sigma_phi_value;
                        const float sdthe = (dthe - mean_theta_value) / sigma_theta_value;

                        const float diff = sqrt(pow(sdphi, 2) + pow(sdthe, 2));

                        const float sigma = 2.0;

                        bool SignTrack = true;
                        if (fabs(sdphi) < sigma && fabs(sdthe) < sigma)
                        {
                            if (diff < min[layer])
                            {
                                min[layer] = diff;
                                mytrk->SetMinDist(diff, layer);
                                mytrk->SetMinDphi(dphi * mytrk->GetChargePrime(), layer);
                                mytrk->SetMinDthe(dthe * mytrk->GetChargePrime(), layer);
                                mytrk->SetMinsDphi(sdphi * mytrk->GetChargePrime(), layer);
                                mytrk->SetMinsDthe(sdthe * mytrk->GetChargePrime(), layer);
                                mytrk->SetHitIndex(ihit, layer);
                            }
                            //vtxhit->AddAssociatedTrack(itrk, diff);
                            mytrk->AddHitCounter(layer);
                            mytrk->SetdPhidThe(iter_layer,dphi,dthe,diff,ihit);
                            iter_num++;
                            if(iter_layer==2 && iassociatedhit >0) numbers[2].push_back(iter_num*10  +numbers[3][iassociatedhit-1]);
                            if(iter_layer==2 && iassociatedhit==0) numbers[2].push_back(iter_num*10  );
                            if(iter_layer==1 && iassociatedhit <  mytrk->GetHitCounter(2)) numbers[1].push_back(iter_num*100 +numbers[2][iassociatedhit]);
                            if(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) numbers[1].push_back(iter_num*100 +numbers[3][iassociatedhit-mytrk->GetHitCounter(2)]);
                            if(iter_layer==0 ) numbers[0].push_back(iter_num*1000+numbers[1][iassociatedhit]);

                        } // end of association
                        else
                        {
                            if (vtxhit->N_AssociatedTracks() > 0)
                                SignTrack = false;
                        }
                        int layer_bin = iter_layer;
                        if(iter_layer==2 && iassociatedhit==0) layer_bin = 4;
                        if (fabs(sdthe) < sigma && SignTrack && is_fill_hsits)
                        {
                            dphi_hist_el[2*layer_bin+charge_bin]->Fill(dphi, dphi_previous_layer, pt);
                            sdphi_hist_el[2*layer_bin+charge_bin]->Fill(sdphi, dphi_previous_layer/sigma_phi_value, pt);
                        }
                        if (fabs(sdphi) < sigma && SignTrack && is_fill_hsits)
                        {
                            dthe_hist_el[2*layer_bin+charge_bin]->Fill(dthe, dthe_previous_layer, pt);
                            sdthe_hist_el[2*layer_bin+charge_bin]->Fill(sdthe, dthe_previous_layer/sigma_theta_value, pt);
                        }
                    } // enf of hit loop
                }
            }
            float min_chi2=100.;
            int final_number = 0;
            for (unsigned int inum = 0; inum < numbers[0].size(); inum++)
            {
                float chi2 = 0.;
                for (int ii1 = 0; ii1 < ((numbers[0][inum] / 1000) >0); ii1++)
                {   
                    int nn = 1;
                    chi2 += mytrk->GetDist(0, ii1);
                    for (int ii2 = 0; ii2 < ((numbers[0][inum] / 100 %10) >0); ii2++)
                    {
                        chi2 += mytrk->GetDist(1, ii2);
                        nn++;
                        for (int ii3 = 0; ii3 < ((numbers[0][inum] / 10 %10) > 0); ii3++)
                        {
                            chi2 += mytrk->GetDist(2, ii3);
                            nn++;
                        }
                        for (int ii4 = 0; ii4 < ((numbers[0][inum] % 10) >0); ii4++)
                        {
                            chi2 += mytrk->GetDist(3, ii4);
                            nn++;
                        }
                        chi2 /= nn;
                        if(chi2<min_chi2) {min_chi2=chi2;final_number=numbers[0][inum];} 
                        chi2_ndf[central_bin]->Fill(chi2, numbers[0][inum] / 1000+numbers[0][inum] / 100 %10+numbers[0][inum] / 10 %10 + numbers[0][inum] % 10 - 3, pt);
                    }
                }
            }
            chi2_ndf[central_bin]->Fill(min_chi2, 19, pt);
                    
            if(min_chi2<3)
            {
                mytrk->SetHitIndex(mytrk->GetHits(0,(int) final_number/1000-0.5), 0);
                mytrk->SetHitIndex(mytrk->GetHits(1,(int) final_number/100 %10-0.5), 1);
                if (mytrk->GetHitCounter(2)>0) {mytrk->SetHitIndex(mytrk->GetHits(2,(int) final_number/10 %10 -0.5 ), 2);mytrk->SetHitCounter(2,1);}
                if (mytrk->GetHitCounter(3)>0) {mytrk->SetHitIndex(mytrk->GetHits(3,(int) final_number%10 -0.5), 3);mytrk->SetHitCounter(3,1);}
                mytrk->SetHitCounter(0,1);mytrk->SetHitCounter(1,1);
            }else{
                mytrk->SetHitCounter(0,0);
            }

        }     // enf of e loop
    }         // end

    void MyEvent::Associate_Hits_to_Hadrons(TH2D *hist_br, TH2D *hist_bz, int is_fill_hsits, TH3D *dphi_hist[N_centr],
                                            TH3D *sdphi_hist[N_centr], TH3D *dthe_hist[N_centr], TH3D *sdthe_hist[N_centr])
    {
        const int nhadron = this->GetNhadron();
        const int nvtxhits = this->GetNVTXhit();
        const int centrality = this->GetCentrality();
        const int rungroup = this->GetRunNumber();
        const int central_bin = (int)centrality / 20;
        if (central_bin > 4 || central_bin < 0)
            return;
        for (int itrk = 0; itrk < nhadron; itrk++)
        {
            MyDileptonAnalysis::MyHadron *mytrk = this->GetHadronEntry(itrk);
            const float pt = mytrk->GetPtPrime();

            const float thetaprime = mytrk->GetThe0Prime();

            float phi0_trk_proj = mytrk->GetPhi0Prime();
            float the0_trk_proj = mytrk->GetThe0Prime();
            const float pz = mytrk->GetPtPrime() * (TMath::Cos(thetaprime)) / (TMath::Sin(thetaprime));

            float rp = sqrt(this->GetPreciseX() * this->GetPreciseX() + this->GetPreciseY() * this->GetPreciseY());
            float zp = this->GetPreciseZ();

            float dilep_phi_projection[total_vtx_layers];
            float dilep_the_projection[total_vtx_layers];
            for (int ii = 0; ii < total_vtx_layers; ii++)
            {
                dilep_phi_projection[ii] = -999;
                dilep_the_projection[ii] = -999;
            }
            for (int p = 1; p < N_steps; p++)
            {
                for (int l = 0; l < total_vtx_layers; l++)
                {
                    if (fabs(rp - radii[l]) < step_size && dilep_phi_projection[l] < -900)
                    {
                        dilep_phi_projection[l] = phi0_trk_proj;
                        dilep_the_projection[l] = the0_trk_proj;
                    }
                }

                const int rbin = hist_bz->GetXaxis()->FindBin(rp);
                const int zbin = hist_bz->GetYaxis()->FindBin(zp);

                const float bz = hist_bz->GetBinContent(rbin, zbin) / 10000;

                const float delta_phi0 = (mytrk->GetChargePrime() * 0.3 * step_size * bz) / (2 * mytrk->GetPtPrime() * 100);
                phi0_trk_proj += delta_phi0;

                const float bradial = hist_br->GetBinContent(rbin, zbin) / 10000;

                const float delta_the0 = 0.3 * bradial * (step_size * TMath::Tan(pi / 2 - the0_trk_proj)) / (2 * pz * 100);

                if (thetaprime > pi / 2)
                    the0_trk_proj -= delta_the0;
                else
                    the0_trk_proj += delta_the0;

                zp += step_size * TMath::Tan(pi / 2 - the0_trk_proj);
                rp += step_size;
            }
            const unsigned int charge_bin = (1 - mytrk->GetChargePrime()) / 2;
            const float quality = mytrk->GetTrkQuality();
            if (quality == 63 || quality == 31 || quality == 51)
                mytrk->SetisERT(1);
            else
                mytrk->SetisERT(0);

            float min[nvtx_layers] = {100, 100, 100, 100};
            for (int ihit = 0; ihit < nvtxhits; ihit++)
            {
                MyDileptonAnalysis::MyVTXHit *vtxhit = this->GetVTXHitEntry(ihit);

                const int layer = vtxhit->GetLayer();

                if (layer < 0 || layer > 3)
                    continue;

                int ilayer = vtxhit->GetiLayer();
                if (ilayer < 0)
                    continue;

                const float phi_hit = vtxhit->GetPhiHit();
                const float theta_hit = vtxhit->GetTheHit();

                const float dphi = (dilep_phi_projection[ilayer] - phi_hit);
                const float dthe = (dilep_the_projection[ilayer] - theta_hit);

                if (abs(dphi) > 0.2 || abs(dthe) > 0.2)
                    continue;

                const float sigma_phi_value = mytrk->get_sigma_phi_data(rungroup, central_bin, layer);
                const float mean_phi_value = mytrk->get_mean_phi_data(rungroup, central_bin, layer);
                const float sigma_theta_value = mytrk->get_sigma_theta_data(rungroup, central_bin, layer);
                const float mean_theta_value = mytrk->get_mean_theta_data(rungroup, central_bin, layer);

                const float sdphi = (dphi - mean_phi_value) / sigma_phi_value;
                const float sdthe = (dthe - mean_theta_value) / sigma_theta_value;

                if (abs(sdphi) > 2 && abs(sdthe) > 2)
                    continue;

                const float diff = sqrt(pow(sdphi, 2) + pow(sdthe, 2));

                bool SignTrack = true;
                if (abs(sdphi) < 2.0 && abs(sdthe) < 2.0)
                {
                    int N_AssociatedTracks = vtxhit->N_AssociatedTracks();
                    for (int iasstrack = 0; iasstrack < N_AssociatedTracks; iasstrack++)
                    {
                        const int CompetitorID = vtxhit->GetAssociatedTrack(iasstrack);
                        MyDileptonAnalysis::MyHadron *CompetitorTrack = this->GetHadronEntry(CompetitorID);
                        if (ihit == CompetitorTrack->GetHitIndex(layer))
                        {
                            if (vtxhit->GetAssociatedTrackDistance(CompetitorID) < diff)
                            {
                                SignTrack = false;
                            }
                            else
                            {
                                if (CompetitorTrack->GetHitCounter(layer) == 1 && mytrk->GetisERT() <= CompetitorTrack->GetisERT())
                                {
                                    SignTrack = false; // maybe delete
                                }
                                else
                                {
                                    vtxhit->DeleteAssociatedTrack(CompetitorID);
                                    CompetitorTrack->RemoveHitCounter(layer);
                                }
                            }
                        }
                    }
                    if (diff < min[layer] && SignTrack)
                    {
                        min[layer] = diff;
                        mytrk->SetMinDist(diff, layer);
                        mytrk->SetHitIndex(ihit, layer);
                        mytrk->AddHitCounter(layer);
                    }
                    vtxhit->AddAssociatedTrack(itrk, diff);
                }
                else
                {
                    if (vtxhit->N_AssociatedTracks() > 0)
                        SignTrack = false;
                }
                if (abs(sdthe) < 2.0 && is_fill_hsits)
                {
                    dphi_hist[central_bin]->Fill(dphi, charge_bin + 2 * layer, pt);
                    sdphi_hist[central_bin]->Fill(sdphi, charge_bin + 2 * layer, pt);
                }
                if (abs(sdphi) < 2.0 && is_fill_hsits)
                {
                    dthe_hist[central_bin]->Fill(dthe, charge_bin + 2 * layer, pt);
                    sdthe_hist[central_bin]->Fill(sdthe, charge_bin + 2 * layer, pt);
                }
            } // enf of hit loop
        }     // end of hadron loop
    }         // end

    void MyEvent::ClearEvent()
    {
        evtno = -999;
        BBCcharge = -999;
        BBCchargeN = -999;
        BBCchargeS = -999;
        BBCtimeN = -999;
        BBCtimeS = -999;
        zvertex = -999;
        centrality = -999;

        preciseX = -999;
        preciseY = -999;
        preciseZ = -999;
        run_number = -999;

        TrackList.clear();
        HadronList.clear();
        VTXHitList.clear();
        ElecCandList.clear();
    }

    int MyEvent::GetRunGroup(int in_run_number)
    {
        for (int irun = 0; irun < N_rg_beam_offset; irun++)
        {
            if (in_run_number >= RunBoarders[irun] && in_run_number < RunBoarders[irun + 1])
                return irun;
        }
        return 0;
    }

    void MyEvent::SetDCA(const unsigned int itrk, const int layer2)
    {
        MyDileptonAnalysis::MyElectron *mytrk = this->GetEntry(itrk);
        if (!mytrk)
            return;

        const float R = mytrk->GetPtPrime() / (0.003 * 0.90);
        const float R2 = R * R;

        MyVTXHit *vtxhit = this->GetVTXHitEntry(mytrk->GetHitIndex(0));
        if (!vtxhit)
            return;

        const float x1 = vtxhit->GetXHit();
        const float y1 = vtxhit->GetYHit();

        vtxhit = this->GetVTXHitEntry(mytrk->GetHitIndex(layer2));
        if (!vtxhit)
            return;

        const float x2 = vtxhit->GetXHit();
        const float y2 = vtxhit->GetYHit();

        const float dx = x2 - x1;
        const float dy = y2 - y1;

        const float dx2 = dx * dx;
        const float dy2 = dy * dy;

        const int arm = (mytrk->GetArm() - 0.5) * 2;
        const int charge = mytrk->GetChargePrime();

        if (R2 < dx2 + dy2)
            return;

        const float delta = arm * charge * 0.5 * sqrt((4 * R2 - dx2 - dy2) / (1 + dy2 / dx2));

        const float x3 = 0.5 * (x1 + x2) - delta * dy / dx;
        const float y3 = 0.5 * (y1 + y2) + delta;

        const float X_circle = x3 - this->GetPreciseX();
        const float Y_circle = y3 - this->GetPreciseY();

        const float L = sqrt(X_circle * X_circle + Y_circle * Y_circle);
        const float dca = (L - R) * 10000; // In Micro Meters

        mytrk->SetDCA(dca);
        const float sdca = dca / (sigma_DCA[0][0][layer2 - 1][0] + sigma_DCA[0][0][layer2 - 1][1] * exp(sigma_DCA[0][0][layer2 - 1][2] * mytrk->GetPtPrime()));
        mytrk->SetsDCA(sdca);

        mytrk->SetDCAX(X_circle * dca / L);
        mytrk->SetDCAY(Y_circle * dca / L);
    }

    void MyEvent::SetDCA2(const unsigned int itrk, const int layer3)
    {
        MyDileptonAnalysis::MyElectron *mytrk = this->GetEntry(itrk);
        if (!mytrk)
            return;

        MyVTXHit *vtxhit = this->GetVTXHitEntry(mytrk->GetHitIndex(0));
        if (!vtxhit)
            return;

        const float x1 = vtxhit->GetXHit();
        const float y1 = vtxhit->GetYHit();

        vtxhit = this->GetVTXHitEntry(mytrk->GetHitIndex(1));
        if (!vtxhit)
            return;

        const float x2 = vtxhit->GetXHit();
        const float y2 = vtxhit->GetYHit();

        vtxhit = this->GetVTXHitEntry(mytrk->GetHitIndex(layer3));
        if (!vtxhit)
            return;

        const float x3 = vtxhit->GetXHit();
        const float y3 = vtxhit->GetYHit();

        const float A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;

        if (A == 0)
            return;
        const float B = (SQR(x1) + SQR(y1)) * (y3 - y2) + (SQR(x2) + SQR(y2)) * (y1 - y3) + (SQR(x3) + SQR(y3)) * (y2 - y1);
        const float C = (SQR(x1) + SQR(y1)) * (x2 - x3) + (SQR(x2) + SQR(y2)) * (x3 - x1) + (SQR(x3) + SQR(y3)) * (x1 - x2);
        const float D = (SQR(x1) + SQR(y1)) * (x3 * y2 - x2 * y3) + (SQR(x2) + SQR(y2)) * (x1 * y3 - x3 * y1) + (SQR(x3) + SQR(y3)) * (x2 * y1 - x1 * y2);

        const float xc = -B / 2 / A;
        const float yc = -C / 2 / A;

        const float X_circle = xc - this->GetPreciseX();
        const float Y_circle = yc - this->GetPreciseY();

        const float R = sqrt((SQR(B) + SQR(C) - 4 * A * D) / 4 / SQR(A));

        const float L = sqrt(X_circle * X_circle + Y_circle * Y_circle);
        const float dca = (L - R) * 10000; // In Micro Meters

        mytrk->SetDCA2(dca);

        mytrk->SetDCAX2(X_circle * dca / L);
        mytrk->SetDCAY2(Y_circle * dca / L);

        const float reconstructed_pt = R * 0.003 * 0.90;
        mytrk->SetReconPT(reconstructed_pt);

        const int arm = (mytrk->GetArm() - 0.5) * 2;

        const int final_charge = (y2 - 0.5 * y1 - 0.5 * y3) / fabs(y2 - 0.5 * y1 - 0.5 * y3) * arm;
        mytrk->SetQ(final_charge);

        /*float xx[3]={x1,x2,x3}, yy[3]={y1,y2,y3};
        TGraph* graph = new TGraph(3,xx,yy);
        TF1 *func;
        if(x1 > 0) func = new TF1("func","pol2", 0, 20);
        else func = new TF1("func","pol2", -20, 0);
        graph->Fit(func,"QR");
        const float par1 = func->GetParameter(1);
        const float par2 = func->GetParameter(2);
        const float slope = par1 + 2*par2*this->GetPreciseX();
        float phi0_new_method = TMath::ATan(slope);
        if((x1 < 0 && y1 > 0) || (x1<0 && y1<0)) phi0_new_method += pi;
        mytrk->SetPhi0(phi0_new_method);*/
    }

    void MyEvent::ReshuffleElectrons()
    {
        int n_electrons = this->GetNtrack();
        for (int itrk = 0; itrk < n_electrons; itrk++)
        {
            MyDileptonAnalysis::MyElectron mytrk = *this->GetEntry(itrk);
            bool do_reshuf = false;

            if (mytrk.GetHitCounter(0) < 1 || mytrk.GetHitCounter(1) < 1 ||
                (mytrk.GetHitCounter(2) < 1 && mytrk.GetHitCounter(3) < 1))
                do_reshuf = true;

            if (do_reshuf)
            {
                this->RemoveTrackEntry(itrk);
                this->AddElecCand(&mytrk);
                n_electrons--;
                itrk--;
                continue;
            }
        }
    }

    void MyEventContainer::Reveal_Hadron()
    {
        const int central_bin = (int)event->GetCentrality() / 20;
        
        event->ReshuffleElectrons();
        const int nleptons = event->GetNtrack();
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            const float pt = mytrk->GetPtPrime();
            const float z1 = mytrk->GetCrkz();
            const float phi1 = mytrk->GetCrkphi();

            const float dcenter_phi_sigma = 0.013;
            const float dcenter_z_sigma = 5.0;

            const int neleccand[2] = {(int)event->GetNeleccand(), nleptons};
            const int first_iter[2] = {0, itrk + 1};
            for (int iway = 0; iway < 2; iway++)
            {
                for (int jtrk = first_iter[iway]; jtrk < neleccand[iway]; jtrk++)
                {
                    MyDileptonAnalysis::MyElectron *myeleccand;
                    if (iway == 0)
                        myeleccand = event->GetElecCand(jtrk);
                    else
                        myeleccand = event->GetEntry(jtrk);

                    const float z2 = myeleccand->GetCrkz();
                    if (z2 < -999)
                        continue;
                    const float phi2 = myeleccand->GetCrkphi();
                    if (phi2 < -999)
                        continue;
                    const int charge_bin = fabs(mytrk->GetChargePrime() - myeleccand->GetChargePrime()) * 5 / 2;

                    const float dcenter_z = (z1 - z2) / dcenter_z_sigma;
                    const float dcenter_phi = (phi1 - phi2) / dcenter_phi_sigma;

                    const float dzed = mytrk->GetZDC() - myeleccand->GetZDC();
                    const float dphi = mytrk->GetPhiDC() - myeleccand->GetPhiDC();
                    const float dalpha = mytrk->GetAlphaPrime() - myeleccand->GetAlphaPrime();

                    float dep1 = mytrk->GetDep();
                    if (dep1 < 0)
                        dep1 = -dep1 * 2.5;
                    float dep2 = myeleccand->GetDep();
                    if (dep2 < 0)
                        dep2 = -dep2 * 2.5;
                    const int EP1 = mytrk->GetEcore() / mytrk->GetPtot() > 0.8;
                    const int EP2 = myeleccand->GetEcore() / myeleccand->GetPtot() > 0.8 || mytrk->GetEcore() < -999;

                    const float dcenter_r = sqrt(dcenter_phi * dcenter_phi + dcenter_z * dcenter_z);
                    const float pt_in_hist = 0.5 * (pt + myeleccand->GetPtPrime());
                    if (fabs(dphi - (0.13 * dalpha)) < 0.015)
                        el_had_dz->Fill(dzed, pt_in_hist, central_bin + charge_bin);
                    el_had_dphi->Fill((dphi - (0.04 * dalpha)) / 0.005, pt_in_hist, central_bin + charge_bin);
                    el_had_dr->Fill(dcenter_r, pt_in_hist, central_bin + charge_bin);

                    if ((fabs(dzed) < 6.0 && fabs(dphi - (0.13 * dalpha)) < 0.015) || fabs(dphi - (0.04 * dalpha)) < 0.015 ||
                        (fabs(dphi - (-0.065 * dalpha)) < 0.015))
                    {
                        if (mytrk->GetGhost() == 0 && (dep1 > dep2 || EP1 < EP2))
                            mytrk->SetGhost(3);
                        if (myeleccand->GetGhost() == 0 && (dep1 < dep2 || EP1 > EP2))
                            myeleccand->SetGhost(3);
                    }
                    if ((fabs(dcenter_z) < 0.01 && fabs(dcenter_phi) < 5) || (fabs(dcenter_phi) < 0.01 && fabs(dcenter_z) < 5))
                    {
                        if (dep1 > dep2 || EP1 < EP2)
                            mytrk->SetGhost(4);
                        if (dep1 < dep2 || EP1 > EP2)
                            myeleccand->SetGhost(4);
                    }

                    if (dcenter_r < 4)
                    {
                        if (dep1 > dep2 || EP1 < EP2)
                            mytrk->SetGhost(1);
                        if (dep1 < dep2 || EP1 > EP2)
                            myeleccand->SetGhost(1);
                        if (dcenter_r < 3 && myeleccand->GetProb() > 0.05 && myeleccand->GetTOFDPHI() > 0 && myeleccand->GetTOFDZ() < 20)
                        {
                            if (dep1 > dep2 || EP1 < EP2)
                                mytrk->SetGhost(2);
                            if (dep1 < dep2 || EP1 > EP2)
                                myeleccand->SetGhost(2);
                        }
                    }
                }
            }
            if (fabs(mytrk->GetEmcdphi()) > 0.05 || fabs(mytrk->GetEmcdz()) > 25)
                continue;
            const int ghost = 5 * mytrk->GetGhost();
            const float EP_new = mytrk->GetEcore() / (mytrk->GetPtot() * mytrk->GetReconPT() / mytrk->GetPtPrime());
            ep_hist->Fill(mytrk->GetEcore() / mytrk->GetPtot(), pt, central_bin + ghost);
            ep_hist_el->Fill(EP_new, pt, central_bin + ghost);
            if (mytrk->GetEcore() / mytrk->GetPtot() > 0.8 && mytrk->GetDep() < 2)
            {
                emc_dphi_el->Fill(mytrk->GetEmcdphi(), pt, central_bin + ghost);
                emc_dz_el->Fill(mytrk->GetEmcdz(), pt, central_bin + ghost);
            }
        }
    }

    void MyEventContainer::CheckVeto()
    {   
        const int centr_bin = event->GetCentrality()/20;
        if(centr_bin<0||centr_bin>4) return;
        for (int ielectron = 0; ielectron < event->GetNtrack(); ielectron++)
        {
            int count = 0;
            MyDileptonAnalysis::MyElectron *electron = event->GetEntry(ielectron);
            const int charge_bin = (1 - electron->GetChargePrime()) / 2;
            const float pt = electron->GetPtPrime();
            if (pt < 0.2)
                continue;
            if (electron->GetHitCounter(0) < 1 || electron->GetHitCounter(1) < 1 ||
                (electron->GetHitCounter(2) < 1 && electron->GetHitCounter(3) < 1))
                continue;
            for (int ilayer = 0; ilayer < 4; ilayer++)
            {
                if(electron->GetHitCounter(ilayer)<1) continue;
                int id_hit = electron->GetHitIndex(ilayer);

                MyDileptonAnalysis::MyVTXHit *hit_orig = event->GetVTXHitEntry(id_hit);

                const float phi_orig = hit_orig->GetPhiHit();
                const float the_orig = hit_orig->GetTheHit();
                
                const int nvtxhits = event->GetNVTXhit();

                for (int ihit = 0; ihit < nvtxhits; ihit++)
                {
                    if(ihit==id_hit) continue;
                    MyDileptonAnalysis::MyVTXHit *vtxhit = event->GetVTXHitEntry(ihit);
                    if(vtxhit->GetLayer()!=ilayer) continue;

                    int layer = vtxhit->GetLayer()-1;
                    if(layer<0) layer=0;

                    const float phi_hit = vtxhit->GetPhiHit();
                    const float theta_hit = vtxhit->GetTheHit();

                    const float dphi = (phi_hit - phi_orig);
                    const float dthe = (theta_hit - the_orig);

                    if (abs(dphi) > 0.05 || abs(dthe) > 0.05 || abs(dphi) < 0.0001)
                        continue;

                    int dphi_index = 1;
                    if (dphi < 0)
                        dphi_index = 0;

                    const float mean = veto_window_mean_par0[layer][charge_bin][dphi_index] + veto_window_mean_par1[layer][charge_bin][dphi_index]*exp(veto_window_mean_par2[layer][charge_bin][dphi_index] * pt);
                    const float sigma = veto_window_sigma_par0[layer][charge_bin][dphi_index] + veto_window_sigma_par1[layer][charge_bin][dphi_index]*exp(veto_window_sigma_par2[layer][charge_bin][dphi_index] * pt);
                    
                     if (abs(dthe) < 0.002 ) veto_hist[centr_bin]->    Fill(dphi,ilayer+4*dphi_index+8*charge_bin,pt);
                    if (abs(dphi) > 0.001 )veto_hist_the[centr_bin]->Fill(dthe,ilayer+4*dphi_index+8*charge_bin,pt);
                    if(fabs(dthe) < 0.002) sveto_hist[centr_bin]->    Fill(fabs(dphi - mean) / sigma,ilayer+4*dphi_index+8*charge_bin,pt);
                    if (fabs(dphi) < 0.05 && fabs(dthe) < 0.001*ilayer)
                    {
                        count++;
                    }
                }
            }
            couter_veto_hist->Fill(count,centr_bin);
        }
    }

    void MyEventContainer::GetHistsFromFile(const std::string loc)
    {
        infile = TFile::Open(loc.c_str(), "read");
        if (!infile)
            std::cout << "NO FILE" << std::endl;
        else
            std::cout << "File opened at " << loc << std::endl;
        hist_br = (TH2D *)infile->Get("hist_br");
        hist_bz = (TH2D *)infile->Get("hist_bz");
    }
    void MyEventContainer::CreateOutFileAndInitHists(std::string outfilename, const int fill_ell, const int fill_had, const int fill_tree, const int fill_dphi,
                                                     const int fill_DCA, const int fill_track_QA, const int fill_reveal, const int fill_true_DCA, const int check_veto)
    {
        outfilename = "my-" + outfilename;
        const int compress = 9;
        if (fill_ell || fill_had || fill_tree || fill_dphi || fill_DCA || fill_track_QA || fill_reveal || fill_true_DCA || check_veto)
            outfile = new TFile(outfilename.c_str(), "RECREATE", outfilename.c_str(), compress);

        if (fill_ell)
        {
            INIT_HISTOS(3, dphi_hist_el,  2*N_centr, 100, -0.1, 0.1, 100, -0.1, 0.1, 50, 0, 5);
            INIT_HISTOS(3, dthe_hist_el,  2*N_centr, 100, -0.1, 0.1, 100, -0.1, 0.1, 50, 0, 5);
            INIT_HISTOS(3, sdphi_hist_el, 2*N_centr, 100, -10, 10, 100, -10, 10, 50, 0, 5);
            INIT_HISTOS(3, sdthe_hist_el, 2*N_centr, 100, -10, 10, 100, -10, 10, 50, 0, 5);
            INIT_HISTOS(3, chi2_ndf, N_centr,      50, 0, 10,  20, 0, 20, 25, 0, 5);
            is_fill_hsits = 1;
        }

        if (fill_had)
        {
            INIT_HISTOS(3, dphi_hist, N_centr, 100, -0.1, 0.1, 8, 0, 8, 50, 0, 5);
            INIT_HISTOS(3, dthe_hist, N_centr, 100, -0.1, 0.1, 8, 0, 8, 50, 0, 5);
            INIT_HISTOS(3, sdphi_hist, N_centr, 100, -10, 10, 8, 0, 8, 50, 0, 5);
            INIT_HISTOS(3, sdthe_hist, N_centr, 100, -10, 10, 8, 0, 8, 50, 0, 5);
            is_fill_hadron_hsits = 1;
        }
        if (fill_tree)
        {
            tree = new TTree("tree", "tree");
            tree->Branch("MyEvent", &ev);
            is_fill_tree = 1;
        }
        if (fill_dphi)
        {
            INIT_HISTOS(3, d_dphi_hist, N_centr, 240, -0.12, 0.12, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, d_dthe_hist, N_centr, 240, -0.12, 0.12, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, DCA_hist, N_centr, 240, -1000, 1000, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sd_dphi_hist, N_centr, 100, -10, 10, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sd_dthe_hist, N_centr, 100, -10, 10, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sDCA_hist, N_centr, 100, -20, 20, 6, 0, 6, 28, 0.2, 3);
            is_fill_dphi_hist = 1;
        }
        if (fill_DCA)
        {
            INIT_HISTOS(3, DCA_2D_hist, N_centr, 100, -1000, 1000, 100, -1000, 1000, 12, 0, 12);
            INIT_HISTOS(3, sDCA_2D_hist, N_centr, 100, -1000, 1000, 100, -1000, 1000, 12, 0, 12);
            INIT_HISTOS(3, DCA2_2D_hist, N_centr, 100, -1000, 1000, 100, -1000, 1000, 12, 0, 12);
            INIT_HISTOS(3, sDCA2_2D_hist, N_centr, 100, -1000, 1000, 100, -1000, 1000, 12, 0, 12);
            is_fill_DCA_hist = 1;
        }
        if (fill_track_QA)
        {
            INIT_HIST(3, temc, 50, -50, 50, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, ttof, 50, -50, 50, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, n0_hist, 10, 0, 10, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, ep_hist, 50, 0, 1.5, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, prob_hist, 50, 0, 1, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, disp_hist, 5, 0, 5, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, chi2npe0_hist, 50, 0, 10, 18, 0.2, 2.0, 5, 0, 5);

            INIT_HIST(3, el_had_dphi, 100, -0.05, 0.05, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, el_had_dz, 100, -50, 50, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, n0_hist_el, 10, 0, 10, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, ep_hist_el, 50, 0, 1.5, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, prob_hist_el, 50, 0, 1, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, disp_hist_el, 5, 0, 5, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, chi2npe0_hist_el, 50, 0, 10, 18, 0.2, 2.0, 5, 0, 5);

            is_fill_track_QA = 1;
        }
        if (fill_reveal)
        {
            INIT_HIST(3, el_had_dphi, 200, -25, 25, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, el_had_dz, 200, -25, 25, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, el_had_dr, 200, 0, 40, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, ep_hist, 25, 0, 1.5, 24, 0.2, 5.0, 25, 0, 25);
            INIT_HIST(3, ep_hist_el, 25, 0, 1.5, 24, 0.2, 5.0, 25, 0, 25);
            INIT_HIST(3, emc_dphi_el, 25, -0.05, 0.05, 24, 0.2, 5.0, 25, 0, 25);
            INIT_HIST(3, emc_dz_el, 25, -50, 50, 24, 0.2, 5.0, 25, 0, 25);
            is_fill_reveal = 1;
        }
        if (fill_true_DCA)
        {
            INIT_HIST(3, DCPT_ReconPT, 50, 0, 5, 50, 0, 5, 10, 0, 10);
            INIT_HIST(3, sDCPT_ReconPT, 500, -0.05, 0.05, 500, -0.05, 0.05, 10, 0, 10);
            INIT_HISTOS(3, DCA12_hist, N_centr, 100, -2000, 2000, 100, -2000, 2000, 12, 0, 12);
            INIT_HISTOS(3, DCA2_hist, N_centr, 200, -4000, 4000, 4, 2, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sDCA2_hist, N_centr, 200, -4000, 4000, 4, 2, 6, 28, 0.2, 3);
            INIT_HIST(3, charge_hist, 8, 0, 8, 24, 0.2, 5.0, 5, 0, 5); // 2bit syst origQ+1, newQ+4, arm+2
            is_fill_DCA2_hist = 1;
        }
        if(check_veto)
        {
            INIT_HISTOS(3, veto_hist,     N_centr, 200, -0.05, 0.05, 16,0,16, 28, 0.2, 3);
            INIT_HISTOS(3, veto_hist_the, N_centr, 200, -0.02, 0.02, 16,0,16, 28, 0.2, 3);
            INIT_HISTOS(3, sveto_hist,    N_centr, 200, 0, 20,       16,0,16, 28, 0.2, 3);
            INIT_HIST(2, couter_veto_hist, 4, 0, 4, 5, 0, 5);
            is_check_veto = 1;
        }
    }
    void MyEventContainer::FillDphiHists()
    {
        const int central_bin = (int)event->GetCentrality() / 20;
        const int nleptons = event->GetNtrack();
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            const float pt = mytrk->GetPtPrime();
            if (pt < 0.2)
                continue;
            if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 ||
                (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1))
                continue;
            if (mytrk->GetMinDphi(0) < 90 && mytrk->GetMinDphi(1) < 90 && (mytrk->GetMinDphi(2) < 90 || mytrk->GetMinDphi(3) < 90))
            {
                const int charge_bin = (1 - mytrk->GetChargePrime()) / 2;
                const float d_dphi0 = mytrk->GetMinDphi(1) - mytrk->GetMinDphi(0);
                const float sd_dphi0 = (d_dphi0 - mean_d_dphi[0][0][charge_bin][0] - mean_d_dphi[0][0][charge_bin][1] / sqrt(pt) - mean_d_dphi[0][0][charge_bin][2] / pt) /
                                       (sigma_d_dphi[0][0][charge_bin][0] + sigma_d_dphi[0][0][charge_bin][1] * exp(sigma_d_dphi[0][0][charge_bin][2] * mytrk->GetPtPrime()));
                for (int ilayer = 3; ilayer > 0; ilayer--)
                {
                    if (mytrk->GetMinDphi(ilayer) < 90)
                    {
                        const int ibin = charge_bin + ilayer * 2 - 2;
                        const float d_dphi = mytrk->GetMinDphi(ilayer) - mytrk->GetMinDphi(ilayer - 1);
                        const float d_dthe = mytrk->GetMinDthe(ilayer) - mytrk->GetMinDthe(ilayer - 1);

                        const float sd_dphi = (d_dphi - mean_d_dphi[0][0][ibin][0] - mean_d_dphi[0][0][ibin][1] / sqrt(pt) - mean_d_dphi[0][0][ibin][2] / pt) /
                                              (sigma_d_dphi[0][0][ibin][0] + sigma_d_dphi[0][0][ibin][1] * exp(sigma_d_dphi[0][0][ibin][2] * mytrk->GetPtPrime()));
                        const float sd_dthe = (d_dthe - mean_d_dthe) /
                                              (sigma_d_dthe[0][0][ibin][0] + sigma_d_dthe[0][0][ibin][1] * exp(sigma_d_dthe[0][0][ibin][2] * mytrk->GetPtPrime()));
                        event->SetDCA(itrk, ilayer);

                        if (is_fill_DCA_hist and ilayer == 1)
                        {
                            float third_bin_input = 0.1 + (1 - mytrk->GetChargePrime()) * 1.5 + 6 * mytrk->GetArm();
                            if (pt > 0.6)
                                third_bin_input += 1;
                            if (pt > 1.2)
                                third_bin_input += 1;
                            DCA_2D_hist[central_bin]->Fill(mytrk->GetDCAX(), mytrk->GetDCAY(), third_bin_input);
                            if (sd_dphi0 < 2 || sd_dphi0 > 8)
                                sDCA_2D_hist[central_bin]->Fill(mytrk->GetDCAX(), mytrk->GetDCAY(), third_bin_input);
                        }
                        if (is_fill_dphi_hist)
                        {
                            const float DCA = mytrk->GetDCA();
                            const float sDCA = mytrk->GetsDCA();
                            d_dphi_hist[central_bin]->Fill(d_dphi, ibin, pt);
                            d_dthe_hist[central_bin]->Fill(d_dthe, ibin, pt);
                            DCA_hist[central_bin]->Fill(DCA, ibin, pt);
                            if ((sd_dphi0 < 2 || sd_dphi0 > 8) || ilayer == 1)
                                sd_dphi_hist[central_bin]->Fill(sd_dphi, ibin, pt);
                            if (sd_dphi0 < 2 || sd_dphi0 > 8)
                                sd_dthe_hist[central_bin]->Fill(sd_dthe, ibin, pt);
                            if (sd_dphi0 < 2 || sd_dphi0 > 8)
                                sDCA_hist[central_bin]->Fill(sDCA, ibin, pt);
                        }

                        mytrk->SetMinDphi(sd_dphi, ilayer);
                        // mytrk->SetMinDthe(sd_dthe,ilayer);
                        if (ilayer > 1)
                        {
                            event->SetDCA2(itrk, ilayer);
                            DCPT_ReconPT->Fill(mytrk->GetReconPT(), pt, central_bin + 5 * (ilayer - 2));
                            //if (sd_dphi0 < 2 || sd_dphi0 > 8)
                                sDCPT_ReconPT->Fill(mytrk->GetPhi0()-mytrk->GetPhi0Prime(), 0., central_bin + 5 * (ilayer - 2));
                            float third_bin_input = 0.1 + (1 - mytrk->GetChargePrime()) * 1.5 + 6 * mytrk->GetArm();
                            if (is_fill_DCA_hist)
                                DCA2_2D_hist[central_bin]->Fill(mytrk->GetDCAX2(), mytrk->GetDCAY2(), third_bin_input);
                            if (is_fill_DCA_hist && (sd_dphi0 < 2 || sd_dphi0 > 8))
                                sDCA2_2D_hist[central_bin]->Fill(mytrk->GetDCAX2(), mytrk->GetDCAY2(), third_bin_input);
                            const float DCA2 = mytrk->GetDCA2();
                            DCA2_hist[central_bin]->Fill(DCA2, ibin, pt);
                            if (sd_dphi0 < 2 || sd_dphi0 > 8)
                                sDCA2_hist[central_bin]->Fill(DCA2, ibin, pt);
                            DCA12_hist[central_bin]->Fill(mytrk->GetDCA2(), mytrk->GetDCA(), third_bin_input);
                        }
                        if (ilayer == 1)
                            charge_hist->Fill(abs(charge_bin - (1 - mytrk->GetCharge()) / 2) * 4 + (1 - mytrk->GetCharge()) + mytrk->GetArm() + 0.1, pt, central_bin);
                    }
                }
            }
            if (is_fill_track_QA)
            {
                temc->Fill(mytrk->GetEmcTOF(), pt, central_bin);
                ttof->Fill(mytrk->GetTOFE() - event->GetBBCtimeN(), pt, central_bin);
                n0_hist_el->Fill(mytrk->GetN0(), pt, central_bin);
                ep_hist_el->Fill(mytrk->GetEcore() / mytrk->GetPtot(), pt, central_bin);
                prob_hist_el->Fill(mytrk->GetProb(), pt, central_bin);
                disp_hist_el->Fill(mytrk->GetDisp(), pt, central_bin);
                chi2npe0_hist_el->Fill(mytrk->GetChi2() / mytrk->GetNpe0(), pt, central_bin);
            }
        }
        if (is_fill_track_QA)
        {
            const int nharons = event->GetNhadron();
            for (int itrk = 0; itrk < nharons; itrk++)
            {
                MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);

                const float pt = mytrk->GetPtPrime();
                if (pt < 0.2 || pt > 2 || mytrk->GetEcore() < -99 || mytrk->GetNpe0() < -99 || mytrk->GetTrkQuality() != 63)
                    continue;
                temc->Fill(mytrk->GetEmcTOF(), pt, central_bin);
                ttof->Fill(mytrk->GetTOFE() - event->GetBBCtimeN(), pt, central_bin);
                if (mytrk->GetDep() < -2 && mytrk->GetDep() > -99)
                    n0_hist->Fill(mytrk->GetN0(), pt, central_bin);
                if (mytrk->GetN0() < 1)
                    ep_hist->Fill(mytrk->GetEcore() / mytrk->GetPtot(), pt, central_bin);
                if (mytrk->GetDep() > -2 || mytrk->GetEcore() / mytrk->GetPtot() > 0.8)
                    continue;
                disp_hist->Fill(mytrk->GetDisp(), pt, central_bin);
                chi2npe0_hist->Fill(mytrk->GetChi2() / mytrk->GetNpe0(), pt, central_bin);
                if (mytrk->GetN0() > 1)
                    continue;
                prob_hist->Fill(mytrk->GetProb(), pt, central_bin);
            }
        }
    }
    void MyEventContainer::WriteOutFile()
    {
        std::cout << "Start writing hists to My outfile" << std::endl;
        infile->Close();
        if (is_fill_tree || is_fill_hadron_hsits || is_fill_hsits || is_fill_dphi_hist || is_fill_DCA_hist || is_fill_track_QA
        || is_fill_reveal || is_fill_DCA2_hist||is_check_veto)
        {
            outfile->cd();
            outfile->Write();
            outfile->Close();
        }
        std::cout << "Hists were written to My outfile" << std::endl;
    }
}
