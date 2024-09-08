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
ClassImp(MyDileptonAnalysis::MyPair)
ClassImp(MyDileptonAnalysis::MyGenTrack)

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

        this->SetPtPrime(this->GetPt() * TMath::Abs(this->GetAlpha() / this->GetAlphaPrime()) * mscale);

        if (this->GetAlpha() * this->GetAlphaPrime() < 0)
            this->SetQPrime(-this->GetCharge());
        else
            this->SetQPrime(this->GetCharge());

        const int DCArm = this->GetArm();

        const float new_the0 = this->GetThe0() - ((bbcz - svxz) / 220) * TMath::Sin(this->GetThe0());
        const float theta_offset = dilep_par0_theta[DCArm][rungroup] * TMath::Sin(new_the0) + dilep_par1_theta[DCArm][rungroup] * TMath::Cos(new_the0) + dilep_par2_theta[DCArm][rungroup];
        this->SetThe0Prime(new_the0 - theta_offset);
    }

    void MyTrack::ResetPrimes(const float bbcz, const float svxz, const int rungroup)
    {
        ///reversing prev corrections
        float phi_offset = -999;
    
        this->SetPhi0(this->GetPhi0() + 2.0195 * (this->GetAlpha() - this->GetAlphaPrime()));

        if (this->GetArm() == 0)
        {
            phi_offset = dilep_par0_phi[0] * TMath::Sin(this->GetPhi0()) + dilep_par1_phi[0] * TMath::Cos(this->GetPhi0()) + dilep_par2_phi[0];
        }
        else
        {
            phi_offset = dilep_par0_phi[1] * TMath::Sin(this->GetPhi0()) + dilep_par1_phi[1] * TMath::Cos(this->GetPhi0()) + dilep_par2_phi[1];
        }

        if (this->GetArm() == 0)
            this->SetPhi0(this->GetPhi0Prime() + phi_offset - res_rot_east);
        else
            this->SetPhi0(this->GetPhi0Prime() + phi_offset - res_rot_west);

        // set Phi0 to right value
        const float alpha_offset = this->GetAlpha() - this->GetAlphaPrime();
        this->SetPhi0(this->GetPhi0() - 2.0195 * alpha_offset*0 );

        ////new correction for phi ant the offset between VTX and DC
        const int DCArm = this->GetArm();
        const int charge = (1-this->GetChargePrime())/2;

        //const float new_phi_offset = phi_offset_params[rungroup][DCArm][0] * TMath::Sin(this->GetPhi0()) + 
        //phi_offset_params[rungroup][DCArm][1] * TMath::Cos(this->GetPhi0()) + phi_offset_params[rungroup][DCArm][2];
        //const float new_phi_offset = ToT_offset[DCArm] + ((fXoffset[DCArm][rungroup] - (fVTXXoffset[rungroup]+1.*(1-2*DCArm))) / 220) * TMath::Sin(this->GetPhiDC()-ToT_offset[DCArm]) +
        //((fYoffset[DCArm][rungroup] - (fVTXYoffset[rungroup]) )/ 220) * TMath::Cos(this->GetPhiDC()-ToT_offset[DCArm]);  
        //const float new_phi_offset = ToT_offset[DCArm] - 2.0195 * 1.005 * (1-2*this->GetArm()) * ( ( fVTXXoffset[rungroup] / 220) * TMath::Sin(this->GetPhiDC()-ToT_offset[DCArm]) -
        //( fVTXYoffset[rungroup] / 220) * TMath::Cos(this->GetPhiDC()-ToT_offset[DCArm]) );
        //const float new_phi_offset = ToT_offset[DCArm] + 2.0195 * ((fXoffset[DCArm][rungroup] - (fVTXXoffset[rungroup])) / 220) * TMath::Sin(this->GetPhiDC()-ToT_offset[DCArm]) +
        //((fYoffset[DCArm][rungroup] - (fVTXYoffset[rungroup]) )/ 220) * TMath::Cos(this->GetPhiDC()-ToT_offset[DCArm])
        //- 0*2.0195 * ( 1. *  DCArm / 220) * TMath::Sin(this->GetPhi0()-ToT_offset[DCArm]);
        //const float new_phi_offset = phi_offset_param0[rungroup][DCArm][charge][corr_layer] * TMath::Sin(this->GetPhi0()) + 
        //phi_offset_param1[rungroup][DCArm][charge][corr_layer] * TMath::Cos(this->GetPhi0()) + phi_offset_param2[rungroup][DCArm][charge][corr_layer];
        //const float new_phi_offset = phi_offset_params[rungroup][DCArm][charge][0] * TMath::Sin(this->GetPhi0()) + 
        //phi_offset_params[rungroup][DCArm][charge][1] * TMath::Cos(this->GetPhi0()) + phi_offset_params[rungroup][DCArm][charge][2];
        //const float new_the0 = this->GetThe0() - ((bbcz - svxz) / 220) * TMath::Sin(this->GetThe0());
        //const float theta_offset = the_offset_params[rungroup][DCArm][0] * TMath::Sin(new_the0) + 
        //the_offset_params[rungroup][DCArm][1] * TMath::Cos(new_the0) + the_offset_params[rungroup][DCArm][2];
        //const float theta_offset = the_offset_param0[rungroup][DCArm][charge][corr_layer] * TMath::Sin(new_the0) + 
        //the_offset_param1[rungroup][DCArm][charge][corr_layer] * TMath::Cos(new_the0) + the_offset_param2[rungroup][DCArm][charge][corr_layer];

        const float new_the0 = this->GetThe0() - ((bbcz - svxz) / 220) * TMath::Sin(this->GetThe0());

        const float theta_offset = the_offset_params[rungroup][DCArm][charge][0] * TMath::Sin(new_the0) + 
        the_offset_params[rungroup][DCArm][charge][1] * TMath::Cos(new_the0) + the_offset_params[rungroup][DCArm][charge][2];

        this->SetThe0Prime( new_the0 - theta_offset);

        const float new_phi_the_offset = phi_the_offset_params[rungroup][DCArm][charge][0] * TMath::Sin(this->GetThe0Prime()) + 
        phi_the_offset_params[rungroup][DCArm][charge][1] * TMath::Cos(this->GetThe0Prime()) + phi_the_offset_params[rungroup][DCArm][charge][2];

        this->SetPhi0(this->GetPhi0() -  new_phi_the_offset);        

        const float new_phi_offset = phi_offset_params[rungroup][DCArm][charge][0] * TMath::Sin(this->GetPhiDC()) + 
        phi_offset_params[rungroup][DCArm][charge][1] * TMath::Cos(this->GetPhiDC()) + phi_offset_params[rungroup][DCArm][charge][2];

        this->SetPhi0Prime(this->GetPhi0() -  new_phi_offset - new_phi_the_offset);
    }

    float MyVTXHit::GetPhiHit(const float xvtx, const float yvtx, const float zvtx) 
    {
        TVector2 hitpoint;
        hitpoint.Set(xhit - xvtx, yhit - yvtx);
        float phihit = hitpoint.Phi();
        if (phihit > 3 * TMath::Pi() / 2)
            phihit -= 2 * TMath::Pi();
        return phihit;
    }

    float MyVTXHit::GetTheHit(const float xvtx, const float yvtx, const float zvtx) 
    {
        TVector3 hitpoint;
        hitpoint.SetXYZ(xhit - xvtx, yhit - yvtx, zhit - zvtx);
        const float thetahit = hitpoint.Theta();
        return thetahit;
    }

    void MyVTXHit::SetiLayerFromR()
    {
        const int layer = this->GetLayer();
        const float rhit = sqrt(xhit * xhit + yhit * yhit);
        int ilayer = 0;
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
        this->SetiLayer(ilayer);
    }

    void MyEventContainer::Associate_Hits_to_Leptons(float sigma, float sigma_veto, float sigma_inner, bool not_fill)
    {
        const int nleptons = event->GetNtrack();
        const int nvtxhits = event->GetNVTXhit();
        const int centrality = event->GetCentrality();
        const int rungroup = event->GetRunNumber();
        const int is_fill_hsits_local = is_fill_hsits && !not_fill;

        const int central_bin = (int)centrality / 20;
        if (central_bin > 4 || central_bin < 0)
            return;
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            mytrk->ZeroHitCounters();
            mytrk->ClearNumberVectors();
            // mytrk = static_cast<MyDileptonAnalysis::MyHadron*>(myhad);
            const float pt = mytrk->GetPtPrime();
            
            if(is_fill_hsits_local)
            {
                charge_recover_hist->Fill(mytrk->GetCharge(),0.,pt);
                charge_recover_hist->Fill(mytrk->GetChargePrime(),1,pt);
                charge_recover_hist->Fill(mytrk->GetPhiDC()>mytrk->GetPhi0()?1:-1,2,pt);
                charge_recover_hist->Fill(mytrk->GetMcId()<3?1:-1,3,pt);
            }

            const float thetaprime = mytrk->GetThe0Prime();

            float phi0_trk_proj = mytrk->GetPhi0Prime();
            float the0_trk_proj = mytrk->GetThe0Prime();
            const float pz = mytrk->GetPtPrime() * (TMath::Cos(thetaprime)) / (TMath::Sin(thetaprime));

            float rp = sqrt(event->GetPreciseX() * event->GetPreciseX() + event->GetPreciseY() * event->GetPreciseY());
            float xp = event->GetPreciseX();
            float yp = event->GetPreciseY();
            float zp = event->GetPreciseZ();

            float phi_now = phi0_trk_proj;
            float the_now = the0_trk_proj;

            float dilep_phi_projection[total_vtx_layers];
            float dilep_the_projection[total_vtx_layers];
            for (int ii = 0; ii < total_vtx_layers; ii++)
            {
                dilep_phi_projection[ii] = -999;
                dilep_the_projection[ii] = -999;
            }
            for (int p = 1; p < N_steps; p++)
            {
                rp = sqrt(SQR(xp) + SQR(yp) );

                for (int l = 0; l < total_vtx_layers; l++)
                {
                    if (TMath::Abs(rp - radii[l]) < step_size && dilep_phi_projection[l] < -900)
                    {
                        dilep_phi_projection[l] = phi0_trk_proj;
                        dilep_the_projection[l] = the0_trk_proj;
                    }
                }

                const int rbin = hist_bz->GetXaxis()->FindBin(rp);
                const int zbin = hist_bz->GetYaxis()->FindBin(zp);

                const float bz = hist_bz->GetBinContent(rbin, zbin) / 10000;

                const float delta_phi0 = (mytrk->GetChargePrime() * 0.3 * step_size * bz) / (2 * mytrk->GetPtPrime() * 100 );
                phi0_trk_proj += delta_phi0;
                phi_now += 2*delta_phi0;

                const float bradial = hist_br->GetBinContent(rbin, zbin) / 10000;

                const float delta_the0 = 0.3 * bradial * (step_size * TMath::Tan(pi / 2 - the_now)) / (2 * pz * 100 );

                if (thetaprime > pi / 2)
                    {the0_trk_proj -= delta_the0; the_now -= 2*delta_the0;}
                else
                    {the0_trk_proj += delta_the0; the_now += 2*delta_the0;}

                zp += step_size * TMath::Tan(pi / 2 - the_now);
                xp += step_size * TMath::Cos(phi_now);
                yp += step_size * TMath::Sin(phi_now);
                //rp += step_size;
            }
            const unsigned int charge_bin = (1 - mytrk->GetChargePrime()) / 2;

            float min[nvtx_layers] = {100, 100, 100, 100};
            
            std::vector<long long> numbers[4];
            long long iter_nums[4] = {0,0,0,0};

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
                    float sdphi_previous_layer = 0;
                    float sdthe_previous_layer = 0;

                    if(iassociatedhit>0&&iter_layer>1)
                    {
                        dphi_previous_layer = mytrk->GetdPhi(iter_layer+1, iassociatedhit-1);
                        dthe_previous_layer = mytrk->GetdThe(iter_layer+1, iassociatedhit-1);
                        sdphi_previous_layer = mytrk->GetsdPhi(iter_layer+1, iassociatedhit-1);
                        sdthe_previous_layer = mytrk->GetsdThe(iter_layer+1, iassociatedhit-1);
                        numbers[3].push_back(iassociatedhit);
                    }
                    if(iter_layer<2 && !(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) )
                    {
                        dphi_previous_layer = mytrk->GetdPhi(iter_layer+1, iassociatedhit);
                        dthe_previous_layer = mytrk->GetdThe(iter_layer+1, iassociatedhit);
                        sdphi_previous_layer = mytrk->GetsdPhi(iter_layer+1, iassociatedhit);
                        sdthe_previous_layer = mytrk->GetsdThe(iter_layer+1, iassociatedhit);
                    }
                    if(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2))
                    {
                        dphi_previous_layer = mytrk->GetdPhi(iter_layer+2, iassociatedhit-mytrk->GetHitCounter(2));
                        dthe_previous_layer = mytrk->GetdThe(iter_layer+2, iassociatedhit-mytrk->GetHitCounter(2));
                        sdphi_previous_layer = mytrk->GetsdPhi(iter_layer+2, iassociatedhit-mytrk->GetHitCounter(2));
                        sdthe_previous_layer = mytrk->GetsdThe(iter_layer+2, iassociatedhit-mytrk->GetHitCounter(2));
                    }

                    for (int ihit = 0; ihit < nvtxhits; ihit++)
                    {
                        MyDileptonAnalysis::MyVTXHit *vtxhit = event->GetVTXHitEntry(ihit);

                        const int layer = vtxhit->GetLayer();

                        if (layer != iter_layer)
                            continue;

                        const int ilayer = vtxhit->GetiLayer();
                        if (ilayer < 0)
                            continue;

                        const float phi_hit = vtxhit->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                        const float theta_hit = vtxhit->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                        
                        const float dphi = (dilep_phi_projection[ilayer] - phi_hit);
                        const float dthe = (dilep_the_projection[ilayer] - theta_hit);
                        if ((TMath::Abs(dphi-2*pi)<0.1||TMath::Abs(dphi+2*pi)<0.1)&&mytrk->GetPtPrime()>0.4) std::cout<<  mytrk->GetChargePrime() << " "<<  mytrk->GetPtPrime() << " "<<dphi<< " "<<dphi-2*pi<< " "<<dphi+2*pi<<std::endl;
                        if (abs(dphi) > 0.1 || abs(dthe) > 0.1) continue;

                        if(vtxhit->GetLadder()>49)vtxhit->SetLadder(vtxhit->GetLadder()-50);

                        float sigma_phi_value = mytrk->get_sigma_phi_data(0*rungroup, central_bin, layer);
                        float mean_phi_value = mytrk->get_mean_phi_data(0*rungroup, central_bin, layer);
                        float sigma_theta_value = mytrk->get_sigma_theta_data(0*rungroup, central_bin, layer);
                        float mean_theta_value = mytrk->get_mean_theta_data(0*rungroup, central_bin, layer);

                        int cycle_layer = layer;
                        if((iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) || iter_layer==2) cycle_layer++;
                        if(iter_layer<2||iassociatedhit>0)
                        {
                            sigma_phi_value   = mytrk->get_dynamic_sigma_phi_data  (0, cycle_layer, dphi_previous_layer);
                            mean_phi_value    = mytrk->get_dynamic_mean_phi_data   (0, cycle_layer, dphi_previous_layer);
                            sigma_theta_value = mytrk->get_dynamic_sigma_theta_data(0, cycle_layer, dthe_previous_layer);
                            mean_theta_value  = mytrk->get_dynamic_mean_theta_data (0, cycle_layer, dthe_previous_layer);
                        }

                        const float sdphi = (dphi - mean_phi_value) / sigma_phi_value;// - mytrk->get_dynamic_smean_phi_data(0, cycle_layer, dphi_previous_layer);
                        const float sdthe = (dthe - mean_theta_value) / sigma_theta_value;

                        const float diff = sqrt(std::pow(sdphi, 2) + std::pow(sdthe, 2));


                        bool SignTrack = true;
                        float sigma_veto0 = sigma, sigma_inner0=sigma;
                        if(layer == 0) {sigma_veto0=sigma_veto;sigma_inner0=sigma_inner;}
                        if ( sdphi*mytrk->GetChargePrime()>-sigma_veto0 && sdphi*mytrk->GetChargePrime() < sigma_inner0 && TMath::Abs(sdthe) < sigma)
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
                            mytrk->SetdPhidThe(iter_layer,dphi,dthe,sdphi,sdthe,diff,ihit);
                            iter_nums[layer]++;
                            if((iter_nums[layer]>99&&layer>=2)||iter_nums[layer]>999) 
                            {
                                std::cout<<layer<<" "<<iter_nums[layer]<< " " <<mytrk->GetChargePrime()<<std::endl;
                                return;
                            }
                            if(iter_layer==2 && iassociatedhit >0) numbers[2].push_back(iter_nums[layer]*100  +numbers[3][iassociatedhit-1]);
                            if(iter_layer==2 && iassociatedhit==0) numbers[2].push_back(iter_nums[layer]*100  );
                            if(iter_layer==1 && iassociatedhit <  mytrk->GetHitCounter(2)) numbers[1].push_back(iter_nums[layer]*10000 +numbers[2][iassociatedhit]);
                            if(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) numbers[1].push_back(iter_nums[layer]*10000 +numbers[3][iassociatedhit-mytrk->GetHitCounter(2)]);
                            if(iter_layer==0 ) numbers[0].push_back(iter_nums[layer]*10000000+numbers[1][iassociatedhit]);

                        } // end of association
                        else
                        {
                            if (vtxhit->N_AssociatedTracks() > 0)
                                SignTrack = false;
                        }
                        int in_arg = 1*mytrk->GetArm()+4*layer+2*charge_bin;
                        if( (layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) || (layer==2 && iassociatedhit>0) ) in_arg+=4;
                        if(iter_layer>1 && iassociatedhit==0) in_arg+=8;

                        if (TMath::Abs(sdthe) < sigma && SignTrack && is_fill_hsits_local)
                        {
                            dphi_hist_el_dynamic[in_arg]->Fill(dphi, dphi_previous_layer, pt);
                            sdphi_hist_el_dynamic[in_arg]->Fill(sdphi, sdphi_previous_layer, pt);
                        }
                        if (sdphi*mytrk->GetChargePrime()>-sigma_veto && sdphi*mytrk->GetChargePrime() < sigma && SignTrack && is_fill_hsits_local)
                        {
                            dthe_hist_el_dynamic[in_arg]->Fill(dthe, dthe_previous_layer, pt);
                            sdthe_hist_el_dynamic[in_arg]->Fill(sdthe, sdthe_previous_layer, pt);
                        }
                    } // enf of hit loop
                } // end of hits in prev layer 
            } //ens of layers   
            if(is_fill_hsits_local)
            {
                charge_recover_hist->Fill(mytrk->GetCharge(),4.,pt);
                charge_recover_hist->Fill(mytrk->GetChargePrime(),5,pt);
                charge_recover_hist->Fill(mytrk->GetPhiDC()>mytrk->GetPhi0()?1:-1,6,pt);
                charge_recover_hist->Fill(mytrk->GetMcId()<3?1:-1,7,pt);
            }
            float min_chi2=1000000.;
            long long final_number = 0;
            for (unsigned int inum = 0; inum < numbers[0].size(); inum++)
            {
                const int inum0 = numbers[0][inum] / 10000000-1;
                const int inum1 = numbers[0][inum] / 10000 %1000-1;
                const int inum2 = numbers[0][inum] / 100 %100-1;
                const int inum3 = numbers[0][inum] %100-1;
                    
                if(inum0>=0 && inum1>=0)
                {   
                    float chi2 = 0;
                    float recon_pt = 0;
                    mytrk->SetHitIndex(mytrk->GetHits(0,inum0), 0);
                    mytrk->SetHitIndex(mytrk->GetHits(1,inum1), 1);
                    if (inum2>=0) mytrk->SetHitIndex(mytrk->GetHits(2,inum2), 2);
                    if (inum3>=0) mytrk->SetHitIndex(mytrk->GetHits(3,inum3), 3);
                    if( inum3>=0) 
                    {
                        event->SetDCA2(itrk,3);
                        recon_pt += mytrk->GetReconPT();
                    }
                    if( inum2>=0) 
                    {
                        event->SetDCA2(itrk,2);
                        recon_pt += mytrk->GetReconPT();
                        if (inum3>=0) recon_pt/=2;
                    }
                    chi2 = TMath::Abs(recon_pt-pt)/pt*30/(2+(int)(inum2>=0)+(int)(inum3>=0));
                    if(chi2<min_chi2) {min_chi2=chi2;final_number=numbers[0][inum];} 
                    
                    if (is_fill_hsits_local&&numbers[0].size()<10) chi2_ndf[central_bin]->Fill(chi2, numbers[0].size(), pt);
                }
            }
            if(is_fill_hsits_local) chi2_ndf[central_bin]->Fill(min_chi2, 19, pt);
            mytrk->SetHitCounter(3,0);mytrk->SetHitCounter(2,0);
            if(min_chi2<800000)
            {

                const int inum0 = final_number / 10000000-1;
                const int inum1 = final_number / 10000 %1000-1;
                const int inum2 = final_number / 100 %100-1;
                const int inum3 = final_number %100-1;
                mytrk->SetHitIndex(mytrk->GetHits(0,inum0), 0);
                mytrk->SetHitIndex(mytrk->GetHits(1,inum1), 1);
                if (inum2>=0) 
                {
                    mytrk->SetHitIndex(mytrk->GetHits(2,inum2 ), 2);
                    mytrk->SetHitCounter(2,1);
                    mytrk->SetMinsDphi(mytrk->GetsdPhi(2, inum2) * mytrk->GetChargePrime(), 2);
                    mytrk->SetMinsDthe(mytrk->GetsdThe(2, inum2) * mytrk->GetChargePrime(), 2);
                }
                if (inum3>=0) 
                {
                    mytrk->SetHitIndex(mytrk->GetHits(3,inum3 ), 3);
                    mytrk->SetHitCounter(3,1);
                    mytrk->SetMinsDphi(mytrk->GetsdPhi(3, inum3) * mytrk->GetChargePrime(), 3);
                    mytrk->SetMinsDthe(mytrk->GetsdThe(3, inum3) * mytrk->GetChargePrime(), 3);
                }
                mytrk->SetHitCounter(0,1);mytrk->SetHitCounter(1,1);
                mytrk->SetMinsDphi(mytrk->GetsdPhi(0, inum0) * mytrk->GetChargePrime(), 0);
                mytrk->SetMinsDthe(mytrk->GetsdThe(0, inum0) * mytrk->GetChargePrime(), 0);
                mytrk->SetMinsDphi(mytrk->GetsdPhi(1, inum1) * mytrk->GetChargePrime(), 1);
                mytrk->SetMinsDthe(mytrk->GetsdThe(1, inum1) * mytrk->GetChargePrime(), 1);
                event->SetDCA(itrk, 1);
                if (mytrk->GetHitCounter(3)>0)  event->SetDCA2(itrk, 3);
                if (mytrk->GetHitCounter(2)>0)  event->SetDCA2(itrk, 2);

                ////////////////////////////////cheking hit assoc effinceincy in sim////////////////////////
                MyVTXHit *vtxhit0 = event->GetVTXHitEntry(mytrk->GetHitIndex(0));
                MyVTXHit *vtxhit1 = event->GetVTXHitEntry(mytrk->GetHitIndex(1));
                MyVTXHit *vtxhit2 = nullptr,*vtxhit3 = nullptr; 
                const int hit_count = 2+mytrk->GetHitCounter(2)+mytrk->GetHitCounter(3);
                if (mytrk->GetHitCounter(2)>0) vtxhit2 = event->GetVTXHitEntry(mytrk->GetHitIndex(2));
                else                           vtxhit2 = event->GetVTXHitEntry(mytrk->GetHitIndex(3));
                if (hit_count==4) vtxhit3 = event->GetVTXHitEntry(mytrk->GetHitIndex(3));
                MyVTXHit *vtxhits[4] = {vtxhit0,vtxhit1,vtxhit2,vtxhit3};
                int istruehitcounter = (hit_count==4)*4;
                for (int i = 0; i < hit_count; i++)
                {
                    if(!vtxhits[i]) std::cout<<"kek"<<std::endl;
                    if(!vtxhits[i]) continue;
                    if(vtxhits[i]->GetSensor() == 0) istruehitcounter++;
                }        
                if(is_fill_hsits_local) truehithist->Fill(istruehitcounter,mytrk->GetPtPrime(),event->GetCentrality());
                if(is_fill_hsits_local) chi2_ndf[central_bin]->Fill(min_chi2, 10+istruehitcounter, pt);
                if(is_fill_hsits_local&&true) 
                                        truehitsigmahist->Fill(istruehitcounter+(sigma-2)*10,mytrk->GetPtPrime(),event->GetCentrality());
                if(is_fill_hsits_local)
                {
                    charge_recover_hist->Fill(mytrk->GetCharge(),8,pt);
                    charge_recover_hist->Fill(mytrk->GetChargePrime(),9,pt);
                    charge_recover_hist->Fill(mytrk->GetPhiDC()>mytrk->GetPhi0()?1:-1,10,pt);
                    charge_recover_hist->Fill(mytrk->GetMcId()<3?1:-1,11,pt);
                }

            }else{
                mytrk->SetHitCounter(0,0);
                if(is_fill_hsits_local)
                {
                    charge_recover_hist->Fill(mytrk->GetCharge(),12,pt);
                    charge_recover_hist->Fill(mytrk->GetChargePrime(),13,pt);
                    charge_recover_hist->Fill(mytrk->GetPhiDC()>mytrk->GetPhi0()?1:-1,14,pt);
                    charge_recover_hist->Fill(mytrk->GetMcId()<3?1:-1,15,pt);
                }
            }
            //mytrk->ClearNumberVectors();
        }     // enf of e loop
    }         // end

    void MyEventContainer::Associate_Hits_to_Hadrons(float sigma)
    {
        const int nhadron = event->GetNhadron();
        const int nvtxhits = event->GetNVTXhit();
        const int centrality = event->GetCentrality();
        const int rungroup = event->GetRunNumber();
        const int central_bin = (int)centrality / 20;
        if (central_bin > 4 || central_bin < 0)
            return;
        for (int itrk = 0; itrk < nhadron; itrk++)
        {
            MyDileptonAnalysis::MyHadron *mytrk = event->GetHadronEntry(itrk);
            const float pt = mytrk->GetPtPrime();

            const float thetaprime = mytrk->GetThe0Prime();

            float phi0_trk_proj = mytrk->GetPhi0Prime()-0.00*mytrk->GetChargePrime()*mytrk->GetArm();
            float the0_trk_proj = mytrk->GetThe0Prime();
            const float pz = mytrk->GetPtPrime() * (TMath::Cos(thetaprime)) / (TMath::Sin(thetaprime));

            float rp = sqrt(event->GetPreciseX() * event->GetPreciseX() + event->GetPreciseY() * event->GetPreciseY());
            float xp = event->GetPreciseX() + 0*mytrk->GetArm();
            float yp = event->GetPreciseY();
            float zp = event->GetPreciseZ();

            float phi_now = phi0_trk_proj;
            float the_now = the0_trk_proj;

            float dilep_phi_projection[total_vtx_layers];
            float dilep_the_projection[total_vtx_layers];
            for (int ii = 0; ii < total_vtx_layers; ii++)
            {
                dilep_phi_projection[ii] = -999;
                dilep_the_projection[ii] = -999;
            }
            for (int p = 1; p < N_steps; p++)
            {
                rp = sqrt(SQR(xp) + SQR(yp) );

                for (int l = 0; l < total_vtx_layers; l++)
                {
                    if (TMath::Abs(rp - radii[l]) < step_size && dilep_phi_projection[l] < -900)
                    {
                        dilep_phi_projection[l] = phi0_trk_proj;
                        dilep_the_projection[l] = the0_trk_proj;
                    }
                }

                const int rbin = hist_bz->GetXaxis()->FindBin(rp);
                const int zbin = hist_bz->GetYaxis()->FindBin(zp);

                const float bz = hist_bz->GetBinContent(rbin, zbin) / 10000;

                const float delta_phi0 = (mytrk->GetChargePrime() * 0.3 * step_size * bz) / (2 * mytrk->GetPtPrime() * 100 );
                phi0_trk_proj += delta_phi0;
                phi_now += 2*delta_phi0;

                const float bradial = hist_br->GetBinContent(rbin, zbin) / 10000;

                const float delta_the0 = 0.3 * bradial * (step_size * TMath::Tan(pi / 2 - the_now)) / (2 * pz * 100 );

                if (thetaprime > pi / 2)
                    {the0_trk_proj -= delta_the0; the_now -= 2*delta_the0;}
                else
                    {the0_trk_proj += delta_the0; the_now += 2*delta_the0;}

                zp += step_size * TMath::Tan(pi / 2 - the_now);
                xp += step_size * TMath::Cos(phi_now);
                yp += step_size * TMath::Sin(phi_now);
                //rp += step_size;
            }
            const unsigned int charge_bin = (1 - mytrk->GetChargePrime()) / 2;

            float min[nvtx_layers] = {100, 100, 100, 100};
            for (int ihit = 0; ihit < nvtxhits; ihit++)
            {
                MyDileptonAnalysis::MyVTXHit *vtxhit = event->GetVTXHitEntry(ihit);

                const int layer = vtxhit->GetLayer();

                if (layer < 0 || layer > 3)
                    continue;

                int ilayer = vtxhit->GetiLayer();
                if (ilayer < 0)
                    continue;

                const float phi_hit = vtxhit->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float theta_hit = vtxhit->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());

                const float dphi = (dilep_phi_projection[ilayer] - phi_hit);
                const float dthe = (dilep_the_projection[ilayer] - theta_hit);

                if (abs(dphi) > 0.1 || abs(dthe) > 0.1)
                    continue;
                
                if(vtxhit->GetLadder()>49)vtxhit->SetLadder(vtxhit->GetLadder()-50);

                const float sigma_phi_value = mytrk->get_sigma_phi_data(0*rungroup, central_bin, layer);
                const float mean_phi_value = mytrk->get_mean_phi_data(0*rungroup, central_bin, layer);
                const float sigma_theta_value = mytrk->get_sigma_theta_data(0*rungroup, central_bin, layer);
                const float mean_theta_value = mytrk->get_mean_theta_data(0*rungroup, central_bin, layer);

                const float sdphi = (dphi - mean_phi_value) / sigma_phi_value;
                const float sdthe = (dthe - mean_theta_value) / sigma_theta_value;

                if (abs(sdphi) > sigma && abs(sdthe) > sigma)
                    continue;

                const float diff = sqrt(std::pow(sdphi, 2) + std::pow(sdthe, 2));

                bool SignTrack = true;
                if (abs(sdphi) < sigma && abs(sdthe) < sigma)
                {
                    int N_AssociatedTracks = vtxhit->N_AssociatedTracks();
                    for (int iasstrack = 0; iasstrack < N_AssociatedTracks; iasstrack++)
                    {
                        const int CompetitorID = vtxhit->GetAssociatedTrack(iasstrack);
                        MyDileptonAnalysis::MyHadron *CompetitorTrack = event->GetHadronEntry(CompetitorID);
                        if (ihit == CompetitorTrack->GetHitIndex(layer))
                        {
                            if (vtxhit->GetAssociatedTrackDistance(CompetitorID) < diff)
                            {
                                SignTrack = false;
                            }
                            else
                            {
                                if (CompetitorTrack->GetHitCounter(layer) == 1 )
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
                    if(false)vtxhit->AddAssociatedTrack(itrk, diff);
                    if(is_fill_hadron_hsits)
                    {
                        myvtx_hist->Fill(event->GetPreciseX()-vtxhit->GetXHit()+sqrt(SQR(vtxhit->GetXHit()-event->GetPreciseX())
                                        +SQR(vtxhit->GetYHit()-event->GetPreciseY()))*TMath::Cos(mytrk->GetPhi0Prime()),event->GetRunGroup(),0.5);
                        myvtx_hist->Fill(event->GetPreciseY()-vtxhit->GetYHit()+sqrt(SQR(vtxhit->GetXHit()-event->GetPreciseX())
                                        +SQR(vtxhit->GetYHit()-event->GetPreciseY()))*TMath::Sin(mytrk->GetPhi0Prime()),event->GetRunGroup(),1.5);
                        myvtx_hist->Fill(event->GetPreciseZ()-vtxhit->GetZHit()+sqrt(SQR(vtxhit->GetXHit()-event->GetPreciseX())
                                        +SQR(vtxhit->GetYHit()-event->GetPreciseY()))/TMath::Tan(mytrk->GetThe0Prime()),event->GetRunGroup(),2.5);
                    }
                }
                else
                {
                    if (vtxhit->N_AssociatedTracks() > 0)
                        SignTrack = false;
                }
                const int hist_2nd_arg = 2*charge_bin + 1 * mytrk->GetArm() + 4 * layer; // for dphi and dthe hists
                if (abs(sdthe) < sigma && is_fill_hadron_hsits)
                {
                    dphi_hist[central_bin] ->Fill( dphi, hist_2nd_arg, pt);
                    sdphi_hist[central_bin]->Fill(sdphi, hist_2nd_arg, pt);
                    const float dphi0 = dphi + mytrk->GetPhi0() - mytrk->GetPhi0Prime();

                    dphi_phi0_init_hist[layer]->Fill(dphi0, mytrk->GetPhi0(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                    dphi_phi0_corr_hist[layer]->Fill(dphi, mytrk->GetPhi0Prime(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                    dphi_the0_init_hist[layer]->Fill(dphi0, mytrk->GetThe0Prime(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                    dphi_the0_corr_hist[layer]->Fill(dphi, mytrk->GetThe0Prime(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                    dthe_phi0_init_hist[layer]->Fill(dphi0, mytrk->GetPhiDC(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                    dthe_phi0_corr_hist[layer]->Fill(dphi, mytrk->GetPhiDC(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                }
                if (abs(sdphi) < sigma && is_fill_hadron_hsits)
                {
                    dthe_hist[central_bin] ->Fill( dthe, hist_2nd_arg, pt);
                    sdthe_hist[central_bin]->Fill(sdthe, hist_2nd_arg, pt);
                    const float newthe0 = mytrk->GetThe0() - ((event->GetVtxZ() - event->GetPreciseZ()) / 220) * TMath::Sin(mytrk->GetThe0());
                    const float dthe0 = dthe + newthe0 - mytrk->GetThe0Prime();
                    dthe_the0_init_hist[layer]->Fill(dthe0, newthe0, 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
                    dthe_the0_corr_hist[layer]->Fill(dthe, mytrk->GetThe0Prime(), 2*mytrk->GetArm() + charge_bin + 4*event->GetRunGroup());
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
        GenTrackList.clear();
    }

    int MyEvent::GetRunGroup(int in_run_number) const
    {
        if(in_run_number == 0) in_run_number = run_number;
        for (int irun = 0; irun < N_rg_beam_offset; irun++)
        {
            if (in_run_number >= RunBoarders[irun] && in_run_number < RunBoarders[irun + 1])
                return irun;
        }
        return 0;
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
        const float phi3 = vtxhit->GetPhiHit(this->GetPreciseX(),this->GetPreciseY(),this->GetPreciseZ());

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
        //const int dir = +1 * mytrk->GetChargePrime();
        //const float phi = mytrk->GetPhi0Prime();
        //const float vx = X_circle - dir * R * sin(phi);
        //const float vy = Y_circle + dir * R * cos(phi);
        mytrk->SetDCA2(dca);

        mytrk->SetDCAX2(X_circle * dca / L);
        mytrk->SetDCAY2(Y_circle * dca / L);

        const float reconstructed_pt = R * 0.003 * 0.90;
        mytrk->SetReconPT(reconstructed_pt);

        if(TMath::Abs(x1)==TMath::Abs(x2)||TMath::Abs(x2)==TMath::Abs(x3)||TMath::Abs(x1)==TMath::Abs(x3))
        {
            std::cout<<"no parabola"<<std::endl;
            return;
        }
        const float a = 1./(x1-x3)*((y1-y2)/(x1-x2)-(y2-y3)/(x2-x3));
        const float b = (y1-y2)*(x3+x2)/(x1-x2)/(x3-x1)-(y2-y3)*(x1+x2)/(x2-x3)/(x3-x1);
        const float c = y1-b*x1-a*x1*x1;
        const float slope = b + 2*a*this->GetPreciseX();
        float phi0_new_method = TMath::ATan(slope);
        if((x1 < 0 && y1 > 0) || (x1<0 && y1<0)) phi0_new_method += pi;
        if(false)mytrk->SetPhi0Prime(phi0_new_method);
        mytrk->SetMinDist(a,0);
        mytrk->SetMinDist(b,1);
        mytrk->SetMinDist(c,2);
        const float xx1 = this->GetPreciseX() + X_circle * dca / L / 10000.;
        const float yy1 = this->GetPreciseY() + Y_circle * dca / L / 10000.;
        const float xx2 = x1;
        const float yy2 = y1;
        const float xx3 = x2;
        const float yy3 = y2;
        const float a1 = 1./(xx1-xx3)*((yy1-yy2)/(xx1-xx2)-(yy2-yy3)/(xx2-xx3));
        const float b1 = (yy1-yy2)*(xx3+xx2)/(xx1-xx2)/(xx3-xx1)-(yy2-yy3)*(xx1+xx2)/(xx2-xx3)/(xx3-xx1);
        //const float c1 = yy1-b1*xx1-a1*xx1*xx1;
        const float slope1 = b1 + 2*a1*xx1;
        float phi0_new_method1 = TMath::ATan(slope1);
        if((x1 < 0 && y1 > 0) || (x1<0 && y1<0)) phi0_new_method1 += pi;
        if(true)mytrk->SetPhi0(phi0_new_method1);
        const int final_charge = (phi3 - phi0_new_method1)  > 0 ? 1: -1;
        mytrk->SetQ(final_charge);
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

    void MyEventContainer::CleanUpHitList()
    {
        int n_hits = event->GetNVTXhit();
        for (int ihit = 0; ihit < n_hits; ihit++)
        {
            MyDileptonAnalysis::MyVTXHit myhit = *event->GetVTXHitEntry(ihit);

            if (myhit.GetLadder()>49)
            {
                event->RemoveVTXHitEntry(ihit);
                n_hits--;
                ihit--;
                continue;
            }
        }
    }

    void MyEventContainer::CheckVeto()
    {   
        const int centr_bin = event->GetCentrality()/20;
        const int centrality = event->GetCentrality();
        if(centr_bin<0||centr_bin>4) return;
        for (int ielectron = 0; ielectron < event->GetNtrack(); ielectron++)
        {
            int count = 0;
            MyDileptonAnalysis::MyElectron *electron = event->GetEntry(ielectron);
            const int charge_bin = (1 - electron->GetChargePrime()) / 2;
            const float pt = electron->GetPtPrime();
            if (pt < 0.2)
                continue;
            if(is_check_veto)
            {
                counter_assoc_eff_hist->Fill(0.5,pt,centrality);
                for (int isigma = 5; isigma > 1; isigma--)
                {
                    const int ientry = 5-isigma;
                    if(TMath::Abs(electron->GetMinsDphi(3))<isigma||TMath::Abs(electron->GetMinsDphi(2))<isigma) 
                    {
                        counter_assoc_eff_hist->Fill(ientry*7+1.5,pt,centrality);
                        if(TMath::Abs(electron->GetMinsDphi(1))<isigma)
                        {
                            counter_assoc_eff_hist->Fill(ientry*7+2.5,pt,centrality);
                            if(TMath::Abs(electron->GetMinsDphi(0))<isigma)
                            {
                                counter_assoc_eff_hist->Fill(ientry*7+3.5,pt,centrality);
                                if(TMath::Abs(electron->GetMinsDphi(3))<isigma&&TMath::Abs(electron->GetMinsDphi(2))<isigma) 
                                    counter_assoc_eff_hist->Fill(ientry*7+4.5,pt,centrality);
                                if(electron->GetMinsDphi(0)>-2) 
                                    counter_assoc_eff_hist->Fill(ientry*7+5.5,pt,centrality);
                                if(electron->GetMinsDphi(0)>-1) 
                                    counter_assoc_eff_hist->Fill(ientry*7+6.5,pt,centrality);
                                if(electron->GetMinsDphi(0)>-0) 
                                    counter_assoc_eff_hist->Fill(ientry*7+7.5,pt,centrality);
                            }
                        }   
                    }
                }
            }
            
            if (electron->GetHitCounter(0) < 1 || electron->GetHitCounter(1) < 1 ||
                (electron->GetHitCounter(2) < 1 && electron->GetHitCounter(3) < 1))
                continue;

            if(is_check_veto) 
            {
                temc->Fill(electron->GetEmcTOF(),pt,centr_bin);
                ttof->Fill(electron->GetTOFE(),pt,centr_bin);
            }

            std::vector<double> prevphis, prevthes;
            std::vector<double> dphivec0, dthevec0, dphivec1, dthevec1;
            std::vector<double> xvec0, yvec0, xvec1, yvec1;
            float phi00=-999, the00=999;
            for (int ilayer = 0; ilayer < 4; ilayer++)
            {
                MyDileptonAnalysis::MyVTXHit *hit_orig = nullptr;
                int id_hit = -999;
                float phi_orig = -999, the_orig = -999, dphi_this = -999, dthe_this = -999;  
                if(electron->GetHitCounter(ilayer)>0) 
                {
                    id_hit = electron->GetHitIndex(ilayer);
                    hit_orig = event->GetVTXHitEntry(id_hit);
                    phi_orig = hit_orig->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                    the_orig = hit_orig->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                    if(ilayer==0)
                    {
                        phi00 = phi_orig;
                        the00 = the_orig;
                    }
                    dphi_this  = (phi00 - phi_orig)*electron->GetChargePrime();
                    dthe_this  = (the00 - the_orig)*electron->GetChargePrime();
                }else{
                    if(ilayer==2 && electron->GetHitCounter(3)>0)
                    {
                        id_hit = electron->GetHitIndex(3);
                        hit_orig = event->GetVTXHitEntry(id_hit);
                        the_orig = hit_orig->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());      
                    }
                    if(ilayer==3 && electron->GetHitCounter(2)>0)
                    {
                        id_hit = electron->GetHitIndex(2);
                        hit_orig = event->GetVTXHitEntry(id_hit);
                        the_orig = hit_orig->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());      
                    }
                }
                float dphi_prev = ilayer == 0 ? 0 : prevphis[prevphis.size()-1];
                float dthe_prev = ilayer == 0 ? 0 : prevthes[prevphis.size()-1];
                const float sigma_phi_value   = electron->get_dynamic_sigma_phi_data  (0, ilayer==0 ? 0 : ilayer-1, dphi_prev);

                if(ilayer<2)
                {
                    prevphis.push_back(dphi_this);
                    prevthes.push_back(dthe_this);
                }

                const int nvtxhits = event->GetNVTXhit();

                for (int ihit = 0; ihit < nvtxhits; ihit++)
                {
                    if(ihit==id_hit) continue;
                    MyDileptonAnalysis::MyVTXHit *vtxhit = event->GetVTXHitEntry(ihit);
                    if(vtxhit->GetLayer()!=ilayer) continue;

                    const int LiLayer = vtxhit->GetiLayer();

                    const float phi_hit = vtxhit->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                    const float theta_hit = vtxhit->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());

                    const float dphi = (phi_hit - phi_orig)*electron->GetChargePrime();
                    const float dthe = (theta_hit - the_orig)*electron->GetChargePrime();

                    const float dphi0 = -(phi_hit - phi00)*electron->GetChargePrime();
                    const float dthe0 = -(theta_hit - the_orig)*electron->GetChargePrime();

                    const float sdphi = (dphi0 ) / sigma_phi_value;
                    const float sdthe = (dthe0 ) * 500;

                    if (abs(dphi0) > 0.5 || abs(dthe0) > 0.5 || abs(dphi) < 0.0005)
                        continue;

                    int dphi_index = 1;
                    if (dphi < 0)
                        dphi_index = 0;

                    if (true&&ilayer==0&&sdphi>0&&sdphi<12.5*4&&TMath::Abs(sdthe)<7) 
                    {
                        dphivec0.push_back(dphi0);
                        dthevec0.push_back(dthe0);
                        xvec0.push_back(vtxhit->GetXHit());
                        yvec0.push_back(vtxhit->GetYHit());
                    }
                    if (true&&ilayer==1&&sdphi>0&&sdphi<12.5*4&&TMath::Abs(sdthe)<7) 
                    {
                        dphivec1.push_back(dphi0);
                        dthevec1.push_back(dthe0);
                        xvec1.push_back(vtxhit->GetXHit());
                        yvec1.push_back(vtxhit->GetYHit());
                    }
                    
                    bool check_2hit = false;
                    if(is_check_veto) 
                    {
                        if (TMath::Abs(dthe) < 0.01  ) veto_phi_hist[centr_bin]->   Fill(dphi,ilayer+4*dphi_index+8*charge_bin,pt);
                        if (TMath::Abs(dphi) < 0.01 + 0.03*ilayer )veto_the_hist[centr_bin]->Fill(dthe,ilayer+4*dphi_index+8*charge_bin,pt);
                        
                        const int countvec0 = (ilayer>0) ? (int) dphivec0.size() : 0;
                        for (int ivec = 0; ivec < countvec0; ivec++)
                        {
                            if(TMath::Abs(sdthe)<7)veto_sphi_sphi_hist[ilayer-1+5*centr_bin]->Fill(dphi0-dphivec0[ivec],dphivec0[ivec],pt);///need a doctor
                            if(sdphi>0&&sdphi<12.5*(ilayer<2?4:4))veto_sthe_sthe_hist[ilayer-1+5*centr_bin]->Fill(dthe0-dthevec0[ivec]*0.5,  dthevec0[ivec]*0.5,pt);  
                            if (TMath::Abs(dthe0-dthevec0[ivec]*0.5)<0.01&&dphi0-dphivec0[ivec]>0) check_2hit = true;    
                        }
                        const int countvec1 = (ilayer>1) ? (int) dphivec1.size() : 0;
                        for (int ivec = 0; ivec < countvec1; ivec++)
                        {
                            const float keff = (  sqrt ( SQR( vtxhit->GetXHit() - event->GetPreciseX() - radii[0]*cos(phi00)) + SQR( vtxhit->GetYHit() - event->GetPreciseY()- radii[0]*sin(phi00)) ) 
                                                / sqrt ( SQR( xvec1[ivec] - event->GetPreciseX() - radii[0]*cos(phi00))       + SQR( yvec1[ivec] - event->GetPreciseY() - radii[0]*sin(phi00)) ) );
                            if(TMath::Abs(sdthe)<2)veto_sphi_sphi_hist[ilayer+1+5*centr_bin]->Fill(dphi0,dphivec1[ivec]*keff+(radii[LiLayer]-radii[0])/(radii[1]-radii[0])*0,pt);///need a doctor
                            if(sdphi>0&&sdphi<12.5*(ilayer<2?1:2))veto_sthe_sthe_hist[ilayer+1+5*centr_bin]->Fill(dthe0-dthevec1[ivec],  dthevec1[ivec],pt);    
                            if (TMath::Abs(dthe0-dthevec1[ivec])<0.01) check_2hit = true;    
                        }
                        if(TMath::Abs(sdthe)<2)veto_sphi_phi_hist[ilayer+4*centr_bin]->Fill(sdphi,  dphi_prev,pt);
                        if(sdphi>0&&sdphi<12.5*(ilayer<2?1:2))veto_sthe_the_hist[ilayer+4*centr_bin]->Fill(sdthe,  dthe_prev,pt);
                        if(centr_bin>2)
                        {
                            if(TMath::Abs(sdthe)<2)veto_phi_phi_hist[ilayer]->Fill(dphi0,dphi_prev,pt);
                            if(sdphi>0&&sdphi<12.5*(ilayer<2?1:2))veto_the_the_hist[ilayer]->Fill(dthe0,dthe_prev,pt);
                        }
                    }
                    if (ilayer<2  && TMath::Abs(dphi) < 0.04+0.04*exp(-2*pt) && TMath::Abs(dthe) < 0.01)
                    {
                        if(electron->GetGhost()<5) electron->SetGhost(ilayer+1);
                        count++;
                    }
                    if (ilayer>1 && TMath::Abs(dphi) < 0.1  && TMath::Abs(dthe) < 0.01 )
                    {
                        if(electron->GetGhost()<5) electron->SetGhost(ilayer+1);
                        count++;
                    }
                    if(electron->GetGhost()==5  && ilayer>0&& sdphi>0           && sdphi<12.5*(ilayer<2?4:4) && TMath::Abs(sdthe)<7) electron->SetGhost(6);
                    if(electron->GetGhost()<5   && ilayer>0&& sdphi>0           && sdphi<12.5*(ilayer<2?4:4) && TMath::Abs(sdthe)<7) electron->SetGhost(5);
                    if(electron->GetGhost()==7  && ilayer>0&& sdphi>0           && sdphi<12.5*(ilayer<2?2:4) && TMath::Abs(sdthe)<5) electron->SetGhost(8);
                    if(electron->GetGhost()<7   && ilayer>0&& sdphi>0           && sdphi<12.5*(ilayer<2?2:4) && TMath::Abs(sdthe)<5) electron->SetGhost(7);
                    if(electron->GetGhost()<10  && ilayer>=0&& sdphi>0 &&check_2hit  && sdphi<12.5*(ilayer<2?4:4) && TMath::Abs(sdthe)<7) electron->SetGhost(10);
                    if(electron->GetGhost()<15  && ilayer>0&& sdphi>0 &&check_2hit  && sdphi<12.5*(ilayer<2?1:2) && TMath::Abs(sdthe)<5) electron->SetGhost(15);
                    if(electron->GetGhost()<20  && ilayer>=0&& sdphi>0 &&check_2hit  && sdphi<12.5*(ilayer<2?3:3) && TMath::Abs(sdthe)<2) electron->SetGhost(20);
                    if(electron->GetGhost()<25  && ilayer>=0&& sdphi>0          && sdphi<12.5*(ilayer<2?1:2) && TMath::Abs(sdthe)<2) electron->SetGhost(ilayer+25);
                }
            }
            if(is_check_veto) 
            {
                couter_veto_hist->Fill(count,pt,centr_bin);
                veto_type_hist->Fill(electron->GetGhost(),pt,centrality);
                counter_assoc_ghost_hist->Fill(0.,pt,centrality);
                if ((TMath::Abs(electron->GetMinsDphi(3))<2.0||TMath::Abs(electron->GetMinsDphi(2))<2.0)&&(TMath::Abs(electron->GetMinsDphi(1))<2.0)&&TMath::Abs((electron->GetMinsDphi(0))<2.0))
                {
                    int cut = electron->GetGhost()<10 ? 10 : 0;
                    counter_assoc_ghost_hist->Fill(1+cut,pt,centrality);
                    if((electron->GetMinsDphi(2)>-1||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>-1||electron->GetHitCounter(3)<1))
                        counter_assoc_ghost_hist->Fill(2+cut,pt,centrality);
                    if((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1))
                        counter_assoc_ghost_hist->Fill(3+cut,pt,centrality);
                    if(electron->GetMinsDphi(0)>-1)
                        counter_assoc_ghost_hist->Fill(4+cut,pt,centrality);
                    if(electron->GetMinsDphi(0)>0)
                        counter_assoc_ghost_hist->Fill(5+cut,pt,centrality);
                    if((electron->GetMinsDphi(2)>-1||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>-1||electron->GetHitCounter(3)<1)&&electron->GetMinsDphi(0)>-1)
                        counter_assoc_ghost_hist->Fill(6+cut,pt,centrality);
                    if((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1)&&electron->GetMinsDphi(0)>-1)
                        counter_assoc_ghost_hist->Fill(7+cut,pt,centrality);
                    if((electron->GetMinsDphi(2)>-1||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>-1||electron->GetHitCounter(3)<1)&&electron->GetMinsDphi(0)>0)
                        counter_assoc_ghost_hist->Fill(8+cut,pt,centrality);
                    if((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1)&&electron->GetMinsDphi(0)>0)
                        counter_assoc_ghost_hist->Fill(9+cut,pt,centrality);
                    ////-----------------------------------------
                    cut = 20;
                    if(((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1))||electron->GetGhost()<20)
                        counter_assoc_ghost_hist->Fill(1+cut,pt,centrality);
                    if(electron->GetMinsDphi(0)>0 || electron->GetGhost()<20)
                        counter_assoc_ghost_hist->Fill(2+cut,pt,centrality);
                    if(((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1))||electron->GetMinsDphi(0)>0||electron->GetGhost()<20)
                        counter_assoc_ghost_hist->Fill(3+cut,pt,centrality);
                    if(((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1)&&electron->GetMinsDphi(0)>0)||electron->GetGhost()<20)
                        counter_assoc_ghost_hist->Fill(4+cut,pt,centrality);
                    ////-----------------------------------------
                    cut = 25;
                    if(((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1))||electron->GetGhost()<10)
                        counter_assoc_ghost_hist->Fill(1+cut,pt,centrality);
                    if(electron->GetMinsDphi(0)>0 || electron->GetGhost()<10)
                        counter_assoc_ghost_hist->Fill(2+cut,pt,centrality);
                    if(((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1))||electron->GetMinsDphi(0)>0||electron->GetGhost()<10)
                        counter_assoc_ghost_hist->Fill(3+cut,pt,centrality);
                    if(((electron->GetMinsDphi(2)>0||electron->GetHitCounter(2)<1)&&(electron->GetMinsDphi(3)>0||electron->GetHitCounter(3)<1)&&electron->GetMinsDphi(0)>0)||electron->GetGhost()<10)
                        counter_assoc_ghost_hist->Fill(4+cut,pt,centrality);
                }
            }
        }
    }

    void MyEventContainer::FillTrueDCA()
    {
        const int central_bin = (int)event->GetCentrality() / 20;
        const int nleptons = event->GetNtrack();
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            const float pt = mytrk->GetPtPrime();
            if (pt < 0.2)
                continue;
            if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1)
                continue;
            const int charge_bin = (1 - mytrk->GetChargePrime()) / 2;
            event->SetDCA(itrk, 1);
            for (int ilayer = 3; ilayer > 1; ilayer--)
            {
                if (mytrk->GetHitCounter(ilayer) < 1) continue;
                const int ibin = charge_bin + ilayer * 2 - 2;
                
                event->SetDCA2(itrk, ilayer);
                if (mytrk->GetMcId()<-1)
                {
                    DCPT_ReconPT->Fill(mytrk->GetReconPT(), pt, central_bin + 5 * (ilayer - 2));
                    if(mytrk->GetGhost()==0)sDCPT_ReconPT->Fill(mytrk->GetReconPT(), pt, central_bin + 5 * (ilayer - 2));
                }
                float third_bin_input = 0.1 + (1 - mytrk->GetChargePrime()) * 1.5 + 6 * mytrk->GetArm();
                if (is_fill_DCA_hist)
                    DCA2_2D_hist[central_bin]->Fill(mytrk->GetDCAX2(), mytrk->GetDCAY2(), third_bin_input);
                if (is_fill_DCA_hist && mytrk->GetGhost()==0)
                    sDCA2_2D_hist[central_bin]->Fill(mytrk->GetDCAX2(), mytrk->GetDCAY2(), third_bin_input);
                const float DCA2 = mytrk->GetDCA2();
                DCA2_hist[central_bin]->Fill(DCA2, ibin, pt);
                if (mytrk->GetGhost()==0)
                    sDCA2_hist[central_bin]->Fill(DCA2, ibin, pt);
                DCA12_hist[central_bin]->Fill(mytrk->GetDCA2(), mytrk->GetDCA(), third_bin_input);
                phi_hist   ->Fill((mytrk->GetPhi0()-mytrk->GetPhi0Prime())*mytrk->GetChargePrime(), pt, central_bin+10);
                the_hist   ->Fill(mytrk->GetThe0()-mytrk->GetThe0Prime(), pt, central_bin+10);
                if(!(mytrk->GetMcId()>-499 && mytrk->GetMcId()-event->GetEvtNo()!=-8))
                {
                    DCPT_ReconPT ->Fill(pt,                  event->GetBBCcharge(), central_bin + 5 * (ilayer - 2));
                    sDCPT_ReconPT->Fill((mytrk->GetReconPT()>pt)?mytrk->GetReconPT():pt, event->GetBBCcharge(), central_bin + 5 * (ilayer - 2));
                    const float dphinew = 1000*(mytrk->GetPhi0Prime()-mytrk->GetPhi0())*mytrk->GetChargePrime();
                    //if(dphinew<0) std::cout<<"Recharge: "<<dphinew<<" "<< pt << " " <<mytrk->GetChargePrime()<<" "<<event->GetEvtNo()<<std::endl;
                    pt_corr ->Fill(dphinew,event->GetBBCcharge(),  pt);
                    charge_hist->Fill( (( mytrk->GetCharge() - event->GetEvtNo()/TMath::Abs(event->GetEvtNo()) )==0)+ (event->GetEvtNo()>0)*2 +  4*mytrk->GetArm(), pt, central_bin);
                    charge_hist->Fill( (( mytrk->GetChargePrime() - event->GetEvtNo()/TMath::Abs(event->GetEvtNo()) )==0)+ (event->GetEvtNo()>0)*2 +  4*mytrk->GetArm(), pt, central_bin+5);
                    phi_hist   ->Fill((mytrk->GetPhi0()-event->GetBBCchargeN())*mytrk->GetChargePrime(), pt, central_bin);
                    phi_hist   ->Fill((mytrk->GetPhi0Prime()-event->GetBBCchargeN())*mytrk->GetChargePrime(), pt, central_bin+5);
                    the_hist   ->Fill(mytrk->GetThe0()-event->GetBBCchargeS(), pt, central_bin);
                    the_hist   ->Fill(mytrk->GetThe0Prime()-event->GetBBCchargeS(), pt, central_bin+5);
                }
                if(!(mytrk->GetMcId()>-499 && mytrk->GetMcId()-event->GetBBCtimeN()!=13))
                {
                    DCPT_ReconPT ->Fill(pt,                  event->GetBBCtimeS(), central_bin + 5 * (ilayer - 2));
                    sDCPT_ReconPT->Fill((mytrk->GetReconPT()>pt)?mytrk->GetReconPT():pt, event->GetBBCtimeS(), central_bin + 5 * (ilayer - 2));
                    const float dphinew = 1000*(mytrk->GetPhi0Prime()-mytrk->GetPhi0())*mytrk->GetChargePrime();
                    //if(dphinew<0) std::cout<<"Recharge: "<<dphinew<<" "<< pt << " " <<mytrk->GetChargePrime()<<" "<<event->GetBBCtimeN()<<std::endl;
                    pt_corr ->Fill(dphinew, event->GetBBCtimeS(), pt);
                    charge_hist->Fill( (( mytrk->GetCharge() - event->GetBBCtimeN()/TMath::Abs(event->GetBBCtimeN()) )==0)+ (event->GetBBCtimeN()>0)*2 +  4*mytrk->GetArm(), pt, central_bin);
                    charge_hist->Fill( (( mytrk->GetChargePrime() - event->GetBBCtimeN()/TMath::Abs(event->GetBBCtimeN()) )==0)+ (event->GetBBCtimeN()>0)*2 +  4*mytrk->GetArm(), pt, central_bin+5);
                    phi_hist   ->Fill((mytrk->GetPhi0()-event->GetPsi3BBC())*mytrk->GetChargePrime(), pt, central_bin);
                    phi_hist   ->Fill((mytrk->GetPhi0Prime()-event->GetPsi3BBC())*mytrk->GetChargePrime(), pt, central_bin+5);
                    the_hist   ->Fill(mytrk->GetThe0()-event->GetPsi3FVTXA0(), pt, central_bin);
                    the_hist   ->Fill(mytrk->GetThe0Prime()-event->GetPsi3FVTXA0(), pt, central_bin+5);
                }
                if((ilayer==3||mytrk->GetHitCounter(3)<1)&&1000*(mytrk->GetPhi0Prime()-mytrk->GetPhi0())*mytrk->GetChargePrime()<-5) 
                mytrk->SetPtPrime((mytrk->GetReconPT()>pt*1.05)?mytrk->GetReconPT():pt);
            }
        }
    }

    void MyEventContainer::fill_evtbuff_list(const unsigned int pool_depth)
    {
        int icent_mix = event->GetCentrality() / 5;
        if (icent_mix > 3)
            icent_mix = 2 + event->GetCentrality() / 10;
        if (icent_mix > 7)
            icent_mix = 8;
        if(icent_mix<0) icent_mix = 9;
        const int izvtx_mix = (event->GetVtxZ() + 10) / 10;
        int ipsi2_mix = (event->GetPsi2FVTXA0() + pi / 2) / pi * 4;
        if(ipsi2_mix<0) ipsi2_mix = 0;
        if(icent_mix<0||izvtx_mix<0||ipsi2_mix<0) std::cout<<"NULL: "<<icent_mix<<" "<<izvtx_mix<<" "<<ipsi2_mix<<std::endl;
        if(icent_mix>=MIX_CENTBIN||izvtx_mix>=MIX_ZVTXBIN||ipsi2_mix>=MIX_RP2BIN) std::cout<<"MAX: "<<icent_mix<<" "<<izvtx_mix<<" "<<ipsi2_mix<<std::endl;
        evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].push_back(*event);
        if (evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].size() > pool_depth)
        {
            evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].pop_front();
        }
    }

    void MyEventContainer::fill_inv_mass(const unsigned int pool_depth)
    {
        int icent_mix = event->GetCentrality() / 5;
        if (icent_mix > 3)
            icent_mix = 2 + event->GetCentrality() / 10;
        if (icent_mix > 7)
            icent_mix = 8;
        if(icent_mix<0) icent_mix = 9;
        const int izvtx_mix = (event->GetVtxZ() + 10) / 10;
        int ipsi2_mix = (event->GetPsi2FVTXA0() + pi / 2) / pi * 4;
        if(ipsi2_mix<0) ipsi2_mix = 0;

        for (int ielectron = 0; ielectron < event->GetNtrack(); ielectron++)
        {
            MyDileptonAnalysis::MyElectron *newTrack1 = event->GetEntry(ielectron);
            if (newTrack1->GetHitCounter(0) < 1 || newTrack1->GetHitCounter(1) < 1 ||
                (newTrack1->GetHitCounter(2) < 1 && newTrack1->GetHitCounter(3) < 1))
                continue;
            //if (!(((newTrack1->GetMinsDphi(2)>0||newTrack1->GetHitCounter(2)<1)&&(newTrack1->GetMinsDphi(3)>0||newTrack1->GetHitCounter(3)<1)&&newTrack1->GetMinsDphi(0)>0)||newTrack1->GetGhost()<10))
            //    continue;
            //const float a1 = newTrack1->GetMinDist(0);
            //const float b1 = newTrack1->GetMinDist(1);
            //const float c1 = newTrack1->GetMinDist(2);
            for (int jelectron = ielectron+1; jelectron < event->GetNtrack(); jelectron++)
            {
                MyDileptonAnalysis::MyElectron *newTrack2 = event->GetEntry(jelectron);
                if (newTrack2->GetHitCounter(0) < 1 || newTrack2->GetHitCounter(1) < 1 ||
                    (newTrack2->GetHitCounter(2) < 1 && newTrack2->GetHitCounter(3) < 1))
                    continue;
                ////////pair cuts
                //const float phi_pip = newTrack1->GetChargePrime() >  newTrack2->GetChargePrime() ? newTrack1->GetPhiDC() : newTrack2->GetPhiDC();
                //const float phi_pim = newTrack1->GetChargePrime() <= newTrack2->GetChargePrime() ? newTrack1->GetPhiDC() : newTrack2->GetPhiDC();
                const float dcenter_phi = (newTrack1->GetChargePrime() * newTrack1->GetCrkphi()+newTrack2->GetChargePrime() * newTrack2->GetCrkphi()) / 0.013;
                const float dcenter_zed = (newTrack1->GetChargePrime() * newTrack1->GetCrkz()+newTrack2->GetChargePrime() * newTrack2->GetCrkz())/5.;
                if ( (TMath::Abs(dcenter_phi)<2|| newTrack1->GetCrkphi()<-99 || newTrack2->GetCrkphi()<-99 || sqrt(SQR(dcenter_zed)+SQR(dcenter_phi)) < 4 ) &&
                      newTrack1->GetChargePrime() != newTrack2->GetChargePrime() && newTrack1->GetArm()== newTrack2->GetArm()) continue;
                const float dalpha = newTrack1->GetAlphaPrime() - newTrack2->GetAlphaPrime();
                const float dphiDC = newTrack1->GetPhiDC() - newTrack2->GetPhiDC();
                const float dzed = newTrack1->GetZDC() - newTrack2->GetZDC();
                if (TMath::Abs(dzed) < 6.0 && TMath::Abs(dphiDC - (0.13 * dalpha)) < 0.015) continue;
                if (TMath::Abs(dphiDC - (0.04 * dalpha)) < 0.015) continue;
                if (TMath::Abs(dphiDC - (-0.065 * dalpha)) < 0.015) continue;
                //////////////end of pair cuts
                
                const int in_hist = (int) (event->GetCentrality() / 20) + N_centr*((newTrack1->GetChargePrime()+newTrack2->GetChargePrime()+2)/2);
                //const float a2 = newTrack2->GetMinDist(0);
                //const float b2 = newTrack2->GetMinDist(1);
                //const float c2 = newTrack2->GetMinDist(2);

                const float dca0 = sqrt(abs(SQR(newTrack1->GetDCA2()) - SQR(newTrack2->GetDCA2())));
                const float dca1 = abs(newTrack1->GetDCAY2() - newTrack2->GetDCAY2());
                const float dca2 = sqrt( SQR(newTrack1->GetDCAX2() - newTrack2->GetDCAX2()) + SQR(newTrack1->GetDCAY2() - newTrack2->GetDCAY2()) );
                const float dca3 = abs(  abs(newTrack1->GetDCAX2() - newTrack2->GetDCAX2()) + abs(newTrack1->GetDCAY2() - newTrack2->GetDCAY2()) );
                const float dca4 = abs(newTrack1->GetDCAY2() / abs(newTrack1->GetDCAY2()) * abs(newTrack1->GetDCA2()) - newTrack2->GetDCAY2() / abs(newTrack2->GetDCAY2()) * abs(newTrack2->GetDCA2()));

                const float pair_pt = sqrt( SQR(newTrack1->GetPx() + newTrack2->GetPx()) + SQR(newTrack1->GetPy() + newTrack2->GetPy()) );

                const float px1 = newTrack1->GetPx();
                const float py1 = newTrack1->GetPy();
                const float pz1 = newTrack1->GetPz();
                const float px2 = newTrack2->GetPx();
                const float py2 = newTrack2->GetPy();
                const float pz2 = newTrack2->GetPz();
                const float pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
                const float pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
                const float es = sqrt(pm1 + me2) + sqrt(pm2 + me2);
                const float px = px1 + px2;
                const float py = py1 + py2;
                const float pz = pz1 + pz2;

                const float invm = sqrt(es * es - px * px - py * py - pz * pz);

                const TVector3 ee1(px1, py1, pz1);
                const TVector3 ee2(px2, py2, pz2);

                const float dphi = ee1.Angle(ee2);
                inv_mass_dca_fg0[in_hist]->Fill(dca0, invm, pair_pt);
                delt_phi_dca_fg0[in_hist]->Fill(dca0, dphi, pair_pt);
                if (!(newTrack1->GetGhost()<25&&newTrack2->GetGhost()<25))
                    continue;
                inv_mass_dca_fg1[in_hist]->Fill(dca1, invm, pair_pt);    
                delt_phi_dca_fg1[in_hist]->Fill(dca1, dphi, pair_pt);
                if (!(newTrack1->GetGhost()<20&&newTrack2->GetGhost()<20))
                    continue;
                inv_mass_dca_fg2[in_hist]->Fill(dca2, invm, pair_pt);
                delt_phi_dca_fg2[in_hist]->Fill(dca2, dphi, pair_pt);
                if (!(newTrack1->GetGhost()<15&&newTrack2->GetGhost()<15))
                    continue;
                inv_mass_dca_fg3[in_hist]->Fill(dca3, invm, pair_pt);
                delt_phi_dca_fg3[in_hist]->Fill(dca3, dphi, pair_pt);
                if(newTrack1->GetHitIndex(0)==newTrack2->GetHitIndex(0)||newTrack1->GetHitIndex(1)==newTrack2->GetHitIndex(1) || 
                  (newTrack1->GetHitIndex(2)==newTrack2->GetHitIndex(2)&&newTrack1->GetHitIndex(3)==newTrack2->GetHitIndex(3))) continue;
                inv_mass_dca_fg4[in_hist]->Fill(dca4, invm, pair_pt);
                delt_phi_dca_fg4[in_hist]->Fill(dca4, dphi, pair_pt);

            }
            const int N_bg_events = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].size();
            for (int ievent = 0; ievent < N_bg_events; ievent++)
            {
                for (int jelectron = 0; jelectron < evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetNtrack(); jelectron++)
                {
                    MyDileptonAnalysis::MyElectron *newTrack2 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetEntry(jelectron);
                    if (newTrack2->GetHitCounter(0) < 1 || newTrack2->GetHitCounter(1) < 1 ||
                        (newTrack2->GetHitCounter(2) < 1 && newTrack2->GetHitCounter(3) < 1))
                        continue;
                    //if (!(((newTrack2->GetMinsDphi(2)>0||newTrack2->GetHitCounter(2)<1)&&(newTrack2->GetMinsDphi(3)>0||newTrack2->GetHitCounter(3)<1)&&newTrack2->GetMinsDphi(0)>0)||newTrack2->GetGhost()<10))
                    //    continue;
                    // const float a2 = newTrack2->GetMinDist(0);
                    // const float b2 = newTrack2->GetMinDist(1);
                    // const float c2 = newTrack2->GetMinDist(2);
                    ////////pair cuts
                    //const float phi_pip = newTrack1->GetChargePrime() >  newTrack2->GetChargePrime() ? newTrack1->GetPhiDC() : newTrack2->GetPhiDC();
                    //const float phi_pim = newTrack1->GetChargePrime() <= newTrack2->GetChargePrime() ? newTrack1->GetPhiDC() : newTrack2->GetPhiDC();
                    const float dcenter_phi = (newTrack1->GetChargePrime() * newTrack1->GetCrkphi()+newTrack2->GetChargePrime() * newTrack2->GetCrkphi()) / 0.013;
                    const float dcenter_zed = (newTrack1->GetChargePrime() * newTrack1->GetCrkz()+newTrack2->GetChargePrime() * newTrack2->GetCrkz())/5.;
                    if ( (TMath::Abs(dcenter_phi)<2|| newTrack1->GetCrkphi()<-99 || newTrack2->GetCrkphi()<-99 || sqrt(SQR(dcenter_zed)+SQR(dcenter_phi)) < 4 ) &&
                          newTrack1->GetChargePrime() != newTrack2->GetChargePrime() && newTrack1->GetArm()== newTrack2->GetArm()) continue;
                    const float dalpha = newTrack1->GetAlpha() - newTrack2->GetAlpha();
                    const float dphiDC = newTrack1->GetPhiDC() - newTrack2->GetPhiDC();
                    const float dzed = newTrack1->GetZDC() - newTrack2->GetZDC();
                    if (TMath::Abs(dzed) < 6.0 && TMath::Abs(dphiDC - (0.13 * dalpha)) < 0.015) continue;
                    if (TMath::Abs(dphiDC - (0.04 * dalpha)) < 0.015) continue;
                    if (TMath::Abs(dphiDC - (-0.065 * dalpha)) < 0.015) continue;
                    //////////////end of pair cuts
                    const int in_hist = (int) (event->GetCentrality() / 20) + N_centr*((newTrack1->GetChargePrime()+newTrack2->GetChargePrime()+2)/2);

                    const float dca0 = sqrt(abs(SQR(newTrack1->GetDCA2()) - SQR(newTrack2->GetDCA2())));
                    const float dca1 = abs(newTrack1->GetDCAY2() - newTrack2->GetDCAY2());
                    const float dca2 = sqrt( SQR(newTrack1->GetDCAX2() - newTrack2->GetDCAX2()) + SQR(newTrack1->GetDCAY2() - newTrack2->GetDCAY2()) );
                    const float dca3 = abs(  abs(newTrack1->GetDCAX2() - newTrack2->GetDCAX2()) + abs(newTrack1->GetDCAY2() - newTrack2->GetDCAY2()) );
                    const float dca4 = abs(newTrack1->GetDCAY2() / abs(newTrack1->GetDCAY2()) * abs(newTrack1->GetDCA2()) - newTrack2->GetDCAY2() / abs(newTrack2->GetDCAY2()) * abs(newTrack2->GetDCA2()));

                    const float pair_pt = sqrt(SQR(newTrack1->GetPx() + newTrack2->GetPx()) + SQR(newTrack1->GetPy() + newTrack2->GetPy()));

                    const float px1 = newTrack1->GetPx();
                    const float py1 = newTrack1->GetPy();
                    const float pz1 = newTrack1->GetPz();
                    const float px2 = newTrack2->GetPx();
                    const float py2 = newTrack2->GetPy();
                    const float pz2 = newTrack2->GetPz();
                    const float pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
                    const float pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
                    const float es = sqrt(pm1 + me2) + sqrt(pm2 + me2);
                    const float px = px1 + px2;
                    const float py = py1 + py2;
                    const float pz = pz1 + pz2;

                    const float invm = sqrt(es * es - px * px - py * py - pz * pz);

                    const TVector3 ee1(px1, py1, pz1);
                    const TVector3 ee2(px2, py2, pz2);
    
                    const float dphi = ee1.Angle(ee2);
                    inv_mass_dca_bg0[in_hist]->Fill(dca0, invm, pair_pt);
                    delt_phi_dca_bg0[in_hist]->Fill(dca0, dphi, pair_pt);
                    if (!(newTrack1->GetGhost()<25&&newTrack2->GetGhost()<25))
                        continue;
                    inv_mass_dca_bg1[in_hist]->Fill(dca1, invm, pair_pt);    
                    delt_phi_dca_bg1[in_hist]->Fill(dca1, dphi, pair_pt);
                    if (!(newTrack1->GetGhost()<20&&newTrack2->GetGhost()<20))
                        continue;
                    inv_mass_dca_bg2[in_hist]->Fill(dca2, invm, pair_pt);
                    delt_phi_dca_bg2[in_hist]->Fill(dca2, dphi, pair_pt);
                    if (!(newTrack1->GetGhost()<15&&newTrack2->GetGhost()<15))
                        continue;
                    inv_mass_dca_bg3[in_hist]->Fill(dca3, invm, pair_pt);
                    delt_phi_dca_bg3[in_hist]->Fill(dca3, dphi, pair_pt);
                    if(newTrack1->GetHitIndex(0)==newTrack2->GetHitIndex(0)||newTrack1->GetHitIndex(1)==newTrack2->GetHitIndex(1) || 
                      (newTrack1->GetHitIndex(2)==newTrack2->GetHitIndex(2)&&newTrack1->GetHitIndex(3)==newTrack2->GetHitIndex(3))) continue;
                    inv_mass_dca_bg4[in_hist]->Fill(dca4, invm, pair_pt);
                    delt_phi_dca_bg4[in_hist]->Fill(dca4, dphi, pair_pt);
                }
            }
        }
        if(event->GetNtrack()>0)this->fill_evtbuff_list(pool_depth);
    }

    void MyEventContainer::fill_inv_mass_sim()
    {
        for (int ielectron = 0; ielectron < event->GetNgentrack(); ielectron++)
        {
            MyDileptonAnalysis::MyGenTrack* newTrack1 = event->GetGenTrack(ielectron);
            for (int jelectron = ielectron+1; jelectron < event->GetNgentrack(); jelectron++)
            {
                MyDileptonAnalysis::MyGenTrack *newTrack2 = event->GetGenTrack(jelectron);
                
                const int in_hist = (int) (event->GetCentrality() / 20) + 5* ( (newTrack1->GetID()>0 ? 1 : 0) + (newTrack2->GetID()>0 ? 1 : 0) );
                const float dca = sqrt( SQR( newTrack1->GetVx() - newTrack2->GetVx() ) + SQR( newTrack1->GetVy() - newTrack2->GetVy() ) ) / 1e13*10000;

                const float pair_pt = sqrt( SQR(newTrack1->GetPx() + newTrack2->GetPx()) + SQR(newTrack1->GetPy() + newTrack2->GetPy()) );

                const float px1 = newTrack1->GetPx();
                const float py1 = newTrack1->GetPy();
                const float pz1 = newTrack1->GetPz();
                const float px2 = newTrack2->GetPx();
                const float py2 = newTrack2->GetPy();
                const float pz2 = newTrack2->GetPz();
                const float pm1 = px1 * px1 + py1 * py1 + pz1 * pz1;
                const float pm2 = px2 * px2 + py2 * py2 + pz2 * pz2;
                const float es = sqrt(pm1 + me2) + sqrt(pm2 + me2);
                const float px = px1 + px2;
                const float py = py1 + py2;
                const float pz = pz1 + pz2;

                const float invm = sqrt(es * es - px * px - py * py - pz * pz);

                inv_mass_dca_gen[in_hist]->Fill(dca, invm, pair_pt);
            }
        }
    }

    int MyEventContainer::GetNGoodElectrons()
    {
        int n_good_el = 0;
        for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            if (mytrk->GetHitCounter(0) > 0 && mytrk->GetHitCounter(1) > 0 &&
               (mytrk->GetHitCounter(2) > 0 || mytrk->GetHitCounter(3) > 0)) 
               {
                    if(is_fill_tree)el_pt_hist->Fill(mytrk->GetPtPrime(),0.5,event->GetCentrality());
                    n_good_el++;
               }else{
                   if(is_fill_tree) el_pt_hist->Fill(mytrk->GetPtPrime(),1.5,event->GetCentrality());
               }
        }
        return n_good_el;
    }

    void MyEventContainer::correct_beam_offset()
    {
        const int n_had = event->GetNhadron();
        for (int i = 0; i < n_had; i++)
        {
            MyDileptonAnalysis::MyHadron *hadron = event->GetHadronEntry(i);
            const float alpha_offset = - (event->GetPreciseX() / 220) * TMath::Sin(hadron->GetPhiDC()) - (event->GetPreciseY() / 220) * TMath::Cos(hadron->GetPhiDC());
     
            hadron->SetAlphaPrime(hadron->GetAlpha() - alpha_offset);
            // set Phi0 to right value
            hadron->SetPhi0Prime(hadron->GetPhi0Prime() - 2.0195 * alpha_offset);

            hadron->SetPtPrime(hadron->GetPtPrime() * TMath::Abs(hadron->GetAlpha() / hadron->GetAlphaPrime()) );

            if (hadron->GetAlpha() * hadron->GetAlphaPrime() < 0)
                hadron->SetQPrime(-hadron->GetChargePrime());
            else
                hadron->SetQPrime(hadron->GetChargePrime());
        }
        const int n_elec = event->GetNtrack();
        for (int i = 0; i < n_elec; i++)
        {
            MyDileptonAnalysis::MyElectron *hadron = event->GetEntry(i);
            const float alpha_offset = - (event->GetPreciseX() / 220) * TMath::Sin(hadron->GetPhiDC()) - (event->GetPreciseY() / 220) * TMath::Cos(hadron->GetPhiDC());
     
            hadron->SetAlphaPrime(hadron->GetAlpha() - alpha_offset);
            // set Phi0 to right value
            hadron->SetPhi0Prime(hadron->GetPhi0Prime() - 2.0195 * alpha_offset);

            hadron->SetPtPrime(hadron->GetPtPrime() * TMath::Abs(hadron->GetAlpha() / hadron->GetAlphaPrime()) );

            if (hadron->GetAlpha() * hadron->GetAlphaPrime() < 0)
                hadron->SetQPrime(-hadron->GetChargePrime());
            else
                hadron->SetQPrime(hadron->GetChargePrime());
        }
    }

    void MyEventContainer::FillQAHist(const int mc_id)
    {
        el_pt_hist->Fill(event->GetBBCcharge(),0.5,event->GetCentrality());
        const int Nelectrons = event->GetNtrack();
        for (int i = 0; i < Nelectrons; i++)
        {
            MyDileptonAnalysis::MyElectron *electron = event->GetEntry(i);
            if (TMath::Abs(electron->GetMcId()-mc_id)>1 && mc_id>-1) continue;
            el_pt_hist->Fill(electron->GetPtPrime(),1.5,event->GetCentrality());
            el_had_dphi->Fill(electron->GetEmcdphi_e(),electron->GetPtPrime(),event->GetCentrality());
            el_had_dz  ->Fill(electron->GetEmcdz_e()  ,electron->GetPtPrime(),event->GetCentrality());
            if(electron->GetEcore()<0 || electron->GetN0()<0|| electron->GetChi2()<0|| electron->GetNpe0()<0|| 
               electron->GetProb()<0|| electron->GetDisp()<0||electron->GetDisp()>10||electron->GetChi2()/electron->GetNpe0()>20) continue;
            el_pt_hist->Fill(electron->GetPtPrime(),2.5,event->GetCentrality());
            ep_hist->Fill(electron->GetEcore()/electron->GetPtot(),electron->GetPtPrime(),event->GetCentrality());
            n0_hist->Fill(electron->GetN0(),electron->GetPtPrime(),event->GetCentrality());
            prob_hist->Fill(electron->GetProb(),electron->GetPtPrime(),event->GetCentrality());
            disp_hist->Fill(electron->GetDisp(),electron->GetPtPrime(),event->GetCentrality());
            chi2npe0_hist->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetPtPrime(),event->GetCentrality());
            const float Rghost = sqrt(SQR(electron->GetEmcdphi_e())+SQR(electron->GetEmcdz_e()));
            el_had_dr->Fill(Rghost,electron->GetPtPrime(),event->GetCentrality());


            if(electron->GetEcore()/electron->GetPtot()>0.8 &&electron->GetEcore()/electron->GetPtot()<1.2 )
            {
                el_pt_hist->Fill(electron->GetPtPrime(),3.5,event->GetCentrality());
                if(electron->GetProb()>0.01)
                {
                    el_pt_hist->Fill(electron->GetPtPrime(),4.5,event->GetCentrality());
                    if(electron->GetN0()>=2 && electron->GetDisp()<5 && electron->GetChi2()/electron->GetNpe0()<10)
                        el_pt_hist->Fill(electron->GetPtPrime(),5.5,event->GetCentrality());
                    if(electron->GetN0()>=3 && electron->GetDisp()<5 && electron->GetChi2()/electron->GetNpe0()<10)
                        el_pt_hist->Fill(electron->GetPtPrime(),6.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8.)
                        el_pt_hist->Fill(electron->GetPtPrime(),7.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8. && electron->GetChi2()/electron->GetNpe0()<10 )
                        el_pt_hist->Fill(electron->GetPtPrime(),8.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8. && electron->GetChi2()/electron->GetNpe0()<12-electron->GetN0()+electron->GetDisp() && electron->GetN0()>electron->GetDisp() )
                        el_pt_hist->Fill(electron->GetPtPrime(),9.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8. && electron->GetDisp()<4 )
                        el_pt_hist->Fill(electron->GetPtPrime(),10.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8. && electron->GetChi2()/electron->GetNpe0()<10 && electron->GetDisp()<4 )
                        el_pt_hist->Fill(electron->GetPtPrime(),11.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8. && electron->GetChi2()/electron->GetNpe0()<10 && electron->GetDisp()<4 && electron->GetProb()>0.03)
                        el_pt_hist->Fill(electron->GetPtPrime(),12.5,event->GetCentrality());
                    if(electron->GetN0() >= 2 + SQR(electron->GetDisp())/8. && electron->GetChi2()/electron->GetNpe0()<10 && electron->GetDisp()<4 && electron->GetProb()>0.03 && 
                       (Rghost>3 || electron->GetN0()- electron->GetDisp() > ((int)electron->GetTOFDPHI())%10 - electron->GetTOFDPHI()/100 || electron->GetTOFE()>electron->GetChi2()/electron->GetNpe0()))
                        el_pt_hist->Fill(electron->GetPtPrime(),13.5,event->GetCentrality());
                    if(electron->GetN0() >= electron->GetDisp() && electron->GetN0() >= 3 && electron->GetChi2()/electron->GetNpe0()<10 )
                        el_pt_hist->Fill(electron->GetPtPrime(),14.5,event->GetCentrality());
                    
                }
            }
            
            //if(Rghost<3 && electron->GetN0()- electron->GetDisp() < ((int)electron->GetTOFDPHI())%10 - electron->GetTOFDPHI()/100 && electron->GetTOFE()<electron->GetChi2()/electron->GetNpe0() ) continue;
            
            ep_hist_el->Fill(electron->GetEcore()/electron->GetPtot(),electron->GetProb(),electron->GetPtPrime());
            n0_hist_el->Fill(electron->GetN0(),electron->GetDisp(),event->GetCentrality());
            prob_hist_el->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetDisp(),event->GetCentrality());
            disp_hist_el->Fill(electron->GetDisp(),electron->GetNpe0(),event->GetCentrality());
            chi2npe0_hist_el->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetNpe0(),event->GetCentrality());
            rich_prob1->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetN0()-1*electron->GetDisp(),event->GetCentrality());
            rich_prob2->Fill(electron->GetNpe0(),electron->GetN0()-1*electron->GetDisp(),event->GetCentrality());
            rich_prob3->Fill(electron->GetEmcdphi(),electron->GetEmcdz(),event->GetCentrality());


            if(electron->GetMcId()>7 && electron->GetEcore()/electron->GetPtot()>0.8 && electron->GetEcore()/electron->GetPtot()<1.2 &&
               electron->GetProb()>0.01 && electron->GetN0() >= 2 + SQR(electron->GetDisp())/8.&& false) 
               std::cout<<electron->GetPtPrime()<<" "<<electron->GetN0()<<" "<<electron->GetEcore()/electron->GetPtot()<<" "<<electron->GetNpe0()
               <<" "<<electron->GetDisp()<<" "<<electron->GetChi2()/electron->GetNpe0()<<" "<<electron->GetProb()
               <<" "<<electron->GetChi2()/electron->GetNpe0()-6-electron->GetDisp()<<" "<<electron->GetN0()-1*electron->GetDisp()<<std::endl;
            //n0_ep, ep_disp, ep_n0-disp+5, disp_n0, chi2_ep
            //TF1 func = TF1("func","pow(x,[0]-1)*exp(-x/2.)/TMath::Gamma([0])/pow(2,[0])",0,100);
            //func.SetNpx(100000);
            //func.SetParameter(0,electron->GetNpe0()/2.);
            //const float p1 = func.Integral(0,electron->GetChi2()/electron->GetNpe0()*electron->GetN0()/2,1.e-3);
            //func.SetParameter(0,electron->GetNpe0()/2.);
            //const float p2 = func.Integral(0,electron->GetChi2()/electron->GetNpe0()*electron->GetDep()*electron->GetN0(),1.e-3);
        }
    
    }

    void MyEventContainer::GetHistsFromFile(const std::string &loc)
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
                                                     const int fill_DCA, const int fill_track_QA, const int fill_reveal, const int fill_true_DCA, 
                                                     const int check_veto, const int fill_inv_mas)
    {
        outfilename = "my-" + outfilename;
        const int compress = 9;
        if (fill_ell || fill_had || fill_tree || fill_dphi || fill_DCA || fill_track_QA || fill_reveal || fill_true_DCA || check_veto || fill_inv_mas)
            outfile = new TFile(outfilename.c_str(), "RECREATE", outfilename.c_str(), compress);

        if (fill_ell)
        {   
            INIT_HISTOS(3, dphi_hist_el_dynamic,  N_dynamic, 100, -0.1, 0.1, 100, -0.1, 0.1, 50, 0, 5);
            INIT_HISTOS(3, dthe_hist_el_dynamic,  N_dynamic, 100, -0.1, 0.1, 100, -0.1, 0.1, 50, 0, 5);
            INIT_HISTOS(3, sdphi_hist_el_dynamic, N_dynamic, 100,  -10,  10, 100,  -10,  10, 50, 0, 5);
            INIT_HISTOS(3, sdthe_hist_el_dynamic, N_dynamic, 100,  -10,  10, 100,  -10,  10, 50, 0, 5);
            INIT_HISTOS(3, chi2_ndf, N_centr,      50, 0, 10,  20, 0, 20, 25, 0, 5);
            INIT_HISTOS(3, dphi_hist_el,  1, 50, -0.1, 0.1, 8, 0, 8, 5, 0, 5);
            INIT_HISTOS(3, dthe_hist_el,  1, 50, -0.1, 0.1, 8, 0, 8, 5, 0, 5);
            INIT_HISTOS(3, sdphi_hist_el, 1, 50, -10, 10,   8, 0, 8, 5, 0, 5);
            INIT_HISTOS(3, sdthe_hist_el, 1, 50, -10, 10,   8, 0, 8, 5, 0, 5);
            INIT_HIST  (3, truehithist,      10, 0, 10, 50, 0, 5, 10, 0, 100);
            INIT_HIST  (3, truehitsigmahist, 50, 0, 50, 50, 0, 5, 10, 0, 100);
            INIT_HIST  (3, charge_recover_hist, 2, -2, 2, 16, 0, 16, 50, 0, 5);
            if(fill_ell>1)
            {
                INIT_HISTOS(3, dphi_hist_el,  N_centr, 100, -0.1, 0.1, 8, 0, 8, 50, 0, 5);
                INIT_HISTOS(3, dthe_hist_el,  N_centr, 100, -0.1, 0.1, 8, 0, 8, 50, 0, 5);
                INIT_HISTOS(3, sdphi_hist_el, N_centr, 100, -10, 10,   8, 0, 8, 50, 0, 5);
                INIT_HISTOS(3, sdthe_hist_el, N_centr, 100, -10, 10,   8, 0, 8, 50, 0, 5);
                is_fill_hsits = 2;
            }
            is_fill_hsits = 1;
        }

        if (fill_had)
        {
            INIT_HISTOS(3, dphi_hist,  N_centr, 100, -0.1, 0.1, 16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, dthe_hist,  N_centr, 100, -0.1, 0.1, 16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, sdphi_hist, N_centr, 100, -10, 10,   16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, sdthe_hist, N_centr, 100, -10, 10,   16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, dphi_phi0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 52, 0, 52);
            INIT_HISTOS(3, dthe_the0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 52, 0, 52);
            INIT_HISTOS(3, dphi_phi0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 52, 0, 52);
            INIT_HISTOS(3, dthe_the0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 52, 0, 52);
            INIT_HISTOS(3, dthe_phi0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 52, 0, 52);
            INIT_HISTOS(3, dphi_the0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 52, 0, 52);
            INIT_HISTOS(3, dthe_phi0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 52, 0, 52);
            INIT_HISTOS(3, dphi_the0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 52, 0, 52);
            INIT_HIST(3, myvtx_hist, 1000, -5, 5, 8 , 0 ,8, 4, 0 ,4);
            is_fill_hadron_hsits = 1;
        }
        if (fill_tree)
        {
            tree = new TTree("tree", "tree");
            tree->Branch("MyEvent", event);
            INIT_HIST(1, event_hist, 10, 0, 10);
            INIT_HIST(1, centr_hist, 100, 0, 100);
            INIT_HIST(3, el_pt_hist, 50, 0, 5, 2, 0, 2, 100, 0, 100);
            is_fill_tree = 1;
        }
        if (fill_dphi)
        {
            INIT_HISTOS(3, d_dphi_hist,  N_centr, 240, -0.12, 0.12, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, d_dthe_hist,  N_centr, 240, -0.12, 0.12, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, DCA_hist,     N_centr, 240, -1000, 1000, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sd_dphi_hist, N_centr, 100,   -10,   10, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sd_dthe_hist, N_centr, 100,   -10,   10, 6, 0, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sDCA_hist,    N_centr, 100,   -20,   20, 6, 0, 6, 28, 0.2, 3);
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
            //INIT_HIST(3, temc, 50, -50, 50, 18, 0.2, 2.0, 5, 0, 5);
            //INIT_HIST(3, ttof, 50, -50, 50, 18, 0.2, 2.0, 5, 0, 5);
            INIT_HIST(3, ep_hist, 50, 0, 1.5, 50, 0, 5.0, 5, 0, 100);
            INIT_HIST(3, n0_hist, 10, 0, 10, 50, 0, 5.0, 5, 0, 100);
            INIT_HIST(3, prob_hist, 100, 0, 1, 50, 0, 5.0, 5, 0, 100);
            INIT_HIST(3, disp_hist, 50, 0, 10, 50, 0, 5.0, 5, 0, 100);
            INIT_HIST(3, chi2npe0_hist, 100, 0, 20, 50, 0., 5.0, 5, 0, 100);
            //n0_ep, ep_disp, ep_n0-disp+5, disp_n0, chi2_ep
            INIT_HIST(3, ep_hist_el, 30, 0, 1.5,  100, 0, 1, 50, 0., 5.0);
            INIT_HIST(3, n0_hist_el, 10, 0, 10, 50, 0, 10, 10, 0., 100);
            INIT_HIST(3, prob_hist_el, 20, 0, 20, 50, 0, 10, 10, 0., 100);
            INIT_HIST(3, disp_hist_el, 50, 0, 10, 30, 0, 30, 10, 0., 100);
            INIT_HIST(3, chi2npe0_hist_el, 50, 0, 20, 30, 0, 30, 10, 0., 100);
            INIT_HIST(3, rich_prob1, 50, 0, 20, 50, -10, 10, 10, 0., 100);
            INIT_HIST(3, rich_prob2, 30, 0, 30, 50, -10, 10, 10, 0., 100);
            INIT_HIST(3, rich_prob3, 100, -0.05, 0.05, 100, -25, 25, 10, 0., 100);

            INIT_HIST(3, el_had_dphi, 100, -10, 10, 50, 0.0, 5.0, 5, 0, 100);
            INIT_HIST(3, el_had_dz, 100, -10, 10, 50, 0.0, 5.0, 5, 0, 100);
            INIT_HIST(3, el_had_dr, 100, 0, 20, 50, 0., 5.0, 5, 0, 100);
            INIT_HIST(3, el_pt_hist, 50, 0, 5, 15, 0, 15, 100, 0, 100);

            //INIT_HIST(3, el_had_dphi, 100, -0.05, 0.05, 24, 0.2, 5.0, 10, 0, 10);
            //INIT_HIST(3, el_had_dz, 100, -50, 50, 24, 0.2, 5.0, 10, 0, 10);
            //INIT_HIST(3, n0_hist_el, 10, 0, 10, 18, 0.2, 2.0, 5, 0, 5);
            //INIT_HIST(3, ep_hist_el, 50, 0, 1.5, 18, 0.2, 2.0, 5, 0, 5);
            //INIT_HIST(3, prob_hist_el, 50, 0, 1, 18, 0.2, 2.0, 5, 0, 5);
            //INIT_HIST(3, disp_hist_el, 5, 0, 5, 18, 0.2, 2.0, 5, 0, 5);
            //INIT_HIST(3, chi2npe0_hist_el, 50, 0, 10, 18, 0.2, 2.0, 5, 0, 5);

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
            INIT_HIST(3, DCPT_ReconPT, 50, 0, 5, 50,  0,  5, 10, 0, 10);
            INIT_HIST(3, sDCPT_ReconPT,50, 0, 5, 50,  0,  5, 10, 0, 10);
            INIT_HIST(3, pt_corr, 500, -50, 50, 50,  0,  5, 50, 0,  5);
            INIT_HISTOS(3, DCA12_hist, N_centr, 100, -2000, 2000, 100, -2000, 2000, 12, 0, 12);
            INIT_HISTOS(3, DCA2_hist, N_centr, 200, -4000, 4000, 4, 2, 6, 28, 0.2, 3);
            INIT_HISTOS(3, sDCA2_hist, N_centr, 200, -4000, 4000, 4, 2, 6, 28, 0.2, 3);
            INIT_HIST(3, charge_hist,   8,     0,    8, 50, 0., 5.0, 10, 0, 10); // 2bit syst origQ+1, newQ+4, arm+2
            INIT_HIST(3,    phi_hist, 300, -0.15, 0.15, 50, 0., 5.0, 15, 0, 15); // 2bit syst origQ+1, newQ+4, arm+2
            INIT_HIST(3,    the_hist, 300, -0.15, 0.15, 50, 0., 5.0, 15, 0, 15); // 2bit syst origQ+1, newQ+4, arm+2
            is_fill_DCA2_hist = 1;
        }
        if(check_veto)
        {
            INIT_HISTOS(3, veto_phi_hist,        N_centr, 150, -0.15, 0.15,  16,    0,  16, 28, 0.2, 3);
            INIT_HISTOS(3, veto_the_hist,        N_centr, 150, -0.15, 0.15,  16,    0,  16, 28, 0.2, 3);
            INIT_HISTOS(3, veto_sphi_phi_hist, 4*N_centr, 150,   -15,   60, 150,-0.15,0.15, 28, 0.2, 3);
            INIT_HISTOS(3, veto_sthe_the_hist, 4*N_centr, 150,   -15,   15, 150,-0.15,0.15, 28, 0.2, 3);
            INIT_HISTOS(3, veto_sphi_sphi_hist,5*N_centr, 150, -0.15, 0.15, 150,-0.15,0.15, 28, 0.2, 3);
            INIT_HISTOS(3, veto_sthe_sthe_hist,5*N_centr, 150, -0.15, 0.15, 150,-0.15,0.15, 28, 0.2, 3);
            INIT_HISTOS(3, veto_phi_phi_hist,    N_centr, 150, -0.15, 0.15, 150,-0.15,0.15, 28, 0.2, 3);
            INIT_HISTOS(3, veto_the_the_hist,    N_centr, 150, -0.15, 0.15, 150,-0.15,0.15, 28, 0.2, 3);
            INIT_HIST(3, couter_veto_hist,         8, 0, 8, 50, 0, 5, 5, 0, 5);
            INIT_HIST(3, counter_assoc_eff_hist,   30,0,30, 50, 0, 5, 10, 0, 100);
            INIT_HIST(3, counter_assoc_ghost_hist, 30,0,30, 50, 0, 5, 10, 0, 100);
            INIT_HIST(3, veto_type_hist,           30,0,30, 50, 0, 5, 10, 0, 100);
            INIT_HIST(3, temc, 150, -50, 150, 50, 0., 2.5, 5, 0, 5);
            INIT_HIST(3, ttof, 1500, -50, 150, 50, 0., 1.25, 5, 0, 5);
            is_check_veto = 1;
        }
        if(fill_inv_mas)
        {
            INIT_HIST(3, el_had_dphi, 200, -25, 25, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, el_had_dz, 200, -25, 25, 24, 0.2, 5.0, 10, 0, 10);
            INIT_HIST(3, el_had_dr, 200, 0, 40, 24, 0.2, 5.0, 10, 0, 10);

            INIT_HISTOS( 3, inv_mass_dca_fg0, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_fg1, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_fg2, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_fg3, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_fg4, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_bg0, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_bg1, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_bg2, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_bg3, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            INIT_HISTOS( 3, inv_mass_dca_bg4, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            
            INIT_HISTOS( 3, delt_phi_dca_fg0, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_fg1, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_fg2, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_fg3, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_fg4, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_bg0, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_bg1, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_bg2, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_bg3, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);
            INIT_HISTOS( 3, delt_phi_dca_bg4, 3*N_centr, 50, 0, 2000, 63, 0, 3.15, 25, 0, 10);

            INIT_HISTOS( 3, inv_mass_dca_gen, 3*N_centr, 50, 0, 2000, 90, 0, 4.50, 25, 0, 10);
            is_fill_inv_mass = 1;
        }
    }
    
    void MyEventContainer::WriteOutFile()
    {
        std::cout << "Start writing hists to My outfile" << std::endl;
        infile->Close();
        if (is_fill_tree || is_fill_hadron_hsits || is_fill_hsits || is_fill_dphi_hist || is_fill_DCA_hist || is_fill_track_QA
        || is_fill_reveal || is_fill_DCA2_hist||is_check_veto||is_fill_inv_mass)
        {
            outfile->cd();
            outfile->Write();
            outfile->Close();
        }
        std::cout << "Hists were written to My outfile" << std::endl;
    }

    /// @yoren no longer in use

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
                        {        
                            charge_hist->Fill(TMath::Abs(charge_bin - (1 - mytrk->GetCharge()) / 2) * 4 + (1 - mytrk->GetCharge()) + mytrk->GetArm() + 0.1, pt, central_bin);
                        }
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
                    const int charge_bin = TMath::Abs(mytrk->GetChargePrime() - myeleccand->GetChargePrime()) * 5 / 2;

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
                    if (TMath::Abs(dphi - (0.13 * dalpha)) < 0.015)
                        el_had_dz->Fill(dzed, pt_in_hist, central_bin + charge_bin);
                    el_had_dphi->Fill((dphi - (0.04 * dalpha)) / 0.005, pt_in_hist, central_bin + charge_bin);
                    el_had_dr->Fill(dcenter_r, pt_in_hist, central_bin + charge_bin);

                    if ((TMath::Abs(dzed) < 6.0 && TMath::Abs(dphi - (0.13 * dalpha)) < 0.015) || TMath::Abs(dphi - (0.04 * dalpha)) < 0.015 ||
                        (TMath::Abs(dphi - (-0.065 * dalpha)) < 0.015))
                    {
                        if (mytrk->GetGhost() == 0 && (dep1 > dep2 || EP1 < EP2))
                            mytrk->SetGhost(3);
                        if (myeleccand->GetGhost() == 0 && (dep1 < dep2 || EP1 > EP2))
                            myeleccand->SetGhost(3);
                    }
                    if ((TMath::Abs(dcenter_z) < 0.01 && TMath::Abs(dcenter_phi) < 5) || (TMath::Abs(dcenter_phi) < 0.01 && TMath::Abs(dcenter_z) < 5))
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
            if (TMath::Abs(mytrk->GetEmcdphi()) > 0.05 || TMath::Abs(mytrk->GetEmcdz()) > 25)
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
    //void MyEventContainer::fill_evtbuff_list(int icent_mix, int izvtx_mix, int ipsi2_mix)
    //{
    //    if (icent_mix < 0 || izvtx_mix < 0 || ipsi2_mix < 0)
    //        return;
    //    evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].push_back(event);
    //    if (debug_mode)
    //        std::cout << "evt buff info: " << icent_mix << " " << izvtx_mix << " " << ipsi2_mix << " depth now " << evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].size() << std::endl;
    //    if (evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].size() > max_evbuf_depth)
    //    {
    //        evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix].pop_front();
    //    }
    //}

}//end of namespace


