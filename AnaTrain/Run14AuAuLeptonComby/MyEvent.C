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
ClassImp(MyDileptonAnalysis::MyBDTrack)
ClassImp(MyDileptonAnalysis::MyBDTHit)

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
        const float theta_offset = dilep_par0_theta[DCArm][rungroup*0] * TMath::Sin(new_the0) + dilep_par1_theta[DCArm][rungroup*0] * TMath::Cos(new_the0) + dilep_par2_theta[DCArm][rungroup*0];
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

        const float new_phi_the_offset = phi_the_offset_params[rungroup*0][DCArm][charge][0] * TMath::Sin(this->GetThe0Prime()) + 
        phi_the_offset_params[rungroup*0][DCArm][charge][1] * TMath::Cos(this->GetThe0Prime()) + phi_the_offset_params[rungroup*0][DCArm][charge][2];

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

    void MyEventContainer::IdenElectrons()
    {
        for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);

            mytrk->SetMcId(0);
            if (mytrk->GetEcore()/mytrk->GetPtot()<0.6 || mytrk->GetN0()<0 ) continue;
        
            if ( mytrk->GetEcore()/mytrk->GetPtot() > 0.8 && mytrk->GetN0()>=2 && mytrk->GetDisp()<5 )  mytrk->SetMcId(mytrk->GetMcId()+1);
            if ( mytrk->GetEcore()/mytrk->GetPtot() > 0.8 && mytrk->GetN0()>=2 + mytrk->GetDisp()*mytrk->GetDisp() / 8. &&  mytrk->GetChi2()/(mytrk->GetNpe0()+0.1)<10 &&
                 mytrk->GetProb()>0.01 && mytrk->GetDisp() < 4)  mytrk->SetMcId(mytrk->GetMcId()+6);
                 
            const float pt = mytrk->GetPtPrime()>0.4?mytrk->GetPtPrime():0.405;
            double treshlods[4] = {0.017808128514646658,  0.0229414147041921,  0.03, 0.04501756860704552};
            if (mytrk->GetPtPrime()<0.4) for (int i = 0; i < 4; i++) treshlods[i] = treshlods[i] + (0.4 - pt) / 20;
            const double input_x[13]=////['centrality', 'pt', 'e/p', 'n0', 'disp', 'chi2', 'npe0', 'prob', 'disp2', 'chi2/npe0', 'centr+pt', 'e/p*pt', 'n0*pt']
            {
                event->GetCentrality(), pt, mytrk->GetEcore()/mytrk->GetPtot(), (double) mytrk->GetN0(), mytrk->GetDisp(), mytrk->GetChi2(), 
                (double) mytrk->GetNpe0(), mytrk->GetProb(), mytrk->GetN0()-SQR(mytrk->GetDisp()), mytrk->GetChi2()/(mytrk->GetNpe0()+0.001), 
                event->GetCentrality()/20.+pt*2, mytrk->GetN0()+4*pt, 
                1./(TMath::Abs( mytrk->GetEcore()/mytrk->GetPtot()-0.9)+0.25)/(1.25-mytrk->GetProb())+4*mytrk->GetPtPrime()
            };
            //mytrk->SetEmcdphi_e(MyML::GetProb(input_x));
            mytrk->SetMcId(mytrk->GetMcId() + MyML::GeteID(input_x, treshlods));
            //if (mytrk->GetMcId()<100 && event->GetCentrality()>20) mytrk->SetMcId(mytrk->GetMcId()+90);
            if (false) std::cout<<mytrk->GetMcId()<<" "<<event->GetCentrality()<<" "<<mytrk->GetPtPrime()<<" "<<mytrk->GetEcore()/mytrk->GetPtot()
                        <<" "<<mytrk->GetN0()<<" "<<mytrk->GetDisp()<<" "<<mytrk->GetChi2()<<" "<<mytrk->GetNpe0()<<std::endl;
            
            if(false) if (mytrk->GetMcId()>10 && mytrk->GetCrkphi()<-99 && mytrk->GetPtPrime()>0.4 ) std::cout<<MyML::GetProb(input_x)<<" "<<mytrk->GetMcId()<<" "<<event->GetCentrality()<<" "<<mytrk->GetPtPrime()<<" "<<mytrk->GetEcore()/mytrk->GetPtot()
                        <<" "<<mytrk->GetN0()<<" "<<mytrk->GetDisp()<<" "<<mytrk->GetChi2()<<" "<<mytrk->GetNpe0()<<" "<<mytrk->GetCrkphi()<<" "<<mytrk->GetEmcdz()<<std::endl;
            
        }    
    }

    void MyEventContainer::Associate_Hits_to_Leptons(float sigma, float sigma_veto, float sigma_inner, int not_fill,int recover_fg, float sigma_theta)
    {
        const int nleptons = event->GetNtrack();
        const int nvtxhits = event->GetNVTXhit();
        const int centrality = event->GetCentrality();
        const int rungroup = event->GetRunGroup();
        is_fill_hsits = !not_fill;

        BDTracklist.clear();

        const int central_bin = (int)centrality / 20;
        if (central_bin > 4 || central_bin < 0)
            return;
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            MyDileptonAnalysis::MyBDTrack newBDTrack = MyDileptonAnalysis::MyBDTrack();

            mytrk->ZeroHitCounters();
            mytrk->ClearNumberVectors();
            mytrk->SetGhost(0);
            const float pt = mytrk->GetPtPrime();

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
            std::vector<std::vector<double> > hit_counter_jlayer[total_vtx_layers+2];

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
                        
                        const float dphi = (dilep_phi_projection[ilayer] - phi_hit) + 
                            (dilep_phi_projection[ilayer<7?ilayer+1:6] - dilep_phi_projection[ilayer])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer<7?ilayer+1:6] - radii[ilayer]);
                        const float dthe = (dilep_the_projection[ilayer] - theta_hit) + 
                            (dilep_the_projection[ilayer<7?ilayer+1:6] - dilep_the_projection[ilayer])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer<7?ilayer+1:6] - radii[ilayer]);
                        if ((TMath::Abs(dphi-2*pi)<0.1||TMath::Abs(dphi+2*pi)<0.1)&&mytrk->GetPtPrime()>0.4) std::cout<<  mytrk->GetChargePrime() << " "<<  mytrk->GetPtPrime() << " "<<dphi<< " "<<dphi-2*pi<< " "<<dphi+2*pi<<std::endl;
                        if (abs(dphi) > 0.1 || abs(dthe) > 0.1) continue;

                        if(vtxhit->GetLadder()>49)vtxhit->SetLadder(vtxhit->GetLadder()-50);
                        if(recover_fg==2) continue;

                        float sigma_phi_value = mytrk->get_sigma_phi_data(0*rungroup, central_bin, layer);
                        float mean_phi_value = mytrk->get_mean_phi_data(0*rungroup, central_bin, layer);
                        float sigma_theta_value = mytrk->get_sigma_theta_data(0*rungroup, central_bin, layer);
                        float mean_theta_value = mytrk->get_mean_theta_data(0*rungroup, central_bin, layer);

                        int cycle_layer = layer;
                        if((iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) || iter_layer==2) cycle_layer++;
                        if(iter_layer<2||iassociatedhit>0)
                        {
                            sigma_phi_value   = mytrk->get_dynamic_sigma_phi_data  (cycle_layer, dphi_previous_layer);
                            mean_phi_value    = mytrk->get_dynamic_mean_phi_data   (cycle_layer, dphi_previous_layer);
                            sigma_theta_value = mytrk->get_dynamic_sigma_theta_data(cycle_layer, dthe_previous_layer);
                            mean_theta_value  = mytrk->get_dynamic_mean_theta_data (cycle_layer, dthe_previous_layer);
                        }

                        const float sdphi = (dphi - mean_phi_value) / sigma_phi_value;// - mytrk->get_dynamic_smean_phi_data(0, cycle_layer, dphi_previous_layer);
                        const float sdthe = (dthe - mean_theta_value) / sigma_theta_value;

                        const float sdphi_rec = (iter_layer==2&&iassociatedhit)?((dphi -  mytrk->get_mean_phi_data(0, central_bin, layer)) / mytrk->get_sigma_phi_data(0, central_bin, layer)):sdphi;
                        const float sdthe_rec = (iter_layer==2&&iassociatedhit)?((dthe -  mytrk->get_mean_theta_data(0, central_bin, layer)) / mytrk->get_sigma_theta_data(0, central_bin, layer)):sdthe;

                        const float diff = sqrt(std::pow(sdphi, 2) + std::pow(sdthe, 2));


                        bool SignTrack = true;
                        float sigma_veto0 = sigma, sigma_inner0=sigma;
                        if(layer == 0) {sigma_veto0=sigma_veto;sigma_inner0=sigma_inner;}
                        if ( sdphi*mytrk->GetChargePrime()>-sigma_veto0 && sdphi*mytrk->GetChargePrime() < sigma_inner0 && TMath::Abs(sdthe) < sigma_theta)
                        {
                            if(false && vtxhit->GetLadder()<25) vtxhit->SetLadder(25+itrk);/// do no taking all hits in cone
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
                            mytrk->SetdPhidThe(iter_layer,dphi,dthe,sdphi_rec,sdthe_rec,diff,ihit);
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
                            std::vector<double> var_counter_jlayer;  
                            var_counter_jlayer.push_back(sdphi);var_counter_jlayer.push_back(sdthe);var_counter_jlayer.push_back(phi_hit);var_counter_jlayer.push_back(theta_hit);
                            hit_counter_jlayer[ilayer].push_back(var_counter_jlayer);
                            if(layer==2) hit_counter_jlayer[8].push_back(var_counter_jlayer);
                            if(layer==3) hit_counter_jlayer[9].push_back(var_counter_jlayer);
                            //std::cout<<ilayer<<" "<<radii[ilayer]<<" "<<sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))<<" "<<dilep_phi_projection[ilayer]<<" "<<mytrk->GetPtPrime()<<" "<<
                            //(dilep_phi_projection[ilayer]-dilep_phi_projection[ilayer<7?ilayer+1:6])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer]-radii[ilayer<7?ilayer+1:6])*1000<<" "<<
                            //(dilep_phi_projection[ilayer]-dilep_phi_projection[ilayer<6?ilayer+2:5])/(radii[ilayer]-radii[ilayer<6?ilayer+2:5])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])*1000<<std::endl;
                        } // end of association
                        else
                        {
                            if (vtxhit->N_AssociatedTracks() > 0)
                                SignTrack = false;
                        }
                        int in_arg = 1*mytrk->GetArm()+4*layer+2*charge_bin;
                        if( (layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) || (layer==2 && iassociatedhit>0) ) in_arg+=4;
                        if(iter_layer>1 && iassociatedhit==0) in_arg+=8;

                        if (TMath::Abs(sdthe) < sigma && SignTrack && is_fill_hsits && !not_fill)
                        {
                            dphi_hist_el_dynamic[in_arg]->Fill(dphi, dphi_previous_layer, pt);
                            sdphi_hist_el_dynamic[in_arg]->Fill(sdphi, sdphi_previous_layer, pt);
                        }
                        if (sdphi*mytrk->GetChargePrime()>-sigma_veto && sdphi*mytrk->GetChargePrime() < sigma && SignTrack && is_fill_hsits)
                        {
                            dthe_hist_el_dynamic[in_arg]->Fill(dthe, dthe_previous_layer, pt);
                            sdthe_hist_el_dynamic[in_arg]->Fill(sdthe, sdthe_previous_layer, pt);
                        }
                    } // enf of hit loop
                } // end of hits in prev layer 
            } //ens of layers   
            if(is_fill_hsits)
            {
                charge_recover_hist->Fill(mytrk->GetCharge(),4.,pt);
                charge_recover_hist->Fill(mytrk->GetChargePrime(),5,pt);
                charge_recover_hist->Fill(mytrk->GetPhiDC()>mytrk->GetPhi0()?1:-1,6,pt);
                charge_recover_hist->Fill(mytrk->GetMcId()<3?1:-1,7,pt);
            }
            float min_chi2=1000000.;
            long long final_number = 0;
            int is_bdt_set = 0;
            for (unsigned int inum = 0; inum < numbers[0].size(); inum++)
            {
                const int inum0 = numbers[0][inum] / 10000000-1;
                const int inum1 = numbers[0][inum] / 10000 %1000-1;
                const int inum2 = numbers[0][inum] / 100 %100-1;
                const int inum3 = numbers[0][inum] %100-1;
                if(inum0>=0 && inum1>=0)
                {   
                    float chi2 = 800001;
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

                    if(true)
                    {
                        if(!is_bdt_set)
                        {
                            newBDTrack.SetAlpha(mytrk->GetAlpha());
                            newBDTrack.SetArm(mytrk->GetArm());
                            newBDTrack.SetCentrality(event->GetCentrality());
                            newBDTrack.SetCharge(mytrk->GetChargePrime());
                            newBDTrack.SetEcore(mytrk->GetEcore());
                            newBDTrack.SetPhi0(mytrk->GetPhi0Prime());
                            newBDTrack.SetPhiDC(mytrk->GetPhiDC());
                            newBDTrack.SetPt(mytrk->GetPtPrime());
                            newBDTrack.SetThe0(mytrk->GetThe0Prime());
                            newBDTrack.SetZDC(mytrk->GetZDC()-event->GetPreciseZ());
                            is_bdt_set = 1;
                        }
                        
                        MyVTXHit *vtxhit0 = event->GetVTXHitEntry(mytrk->GetHitIndex(0));
                        MyVTXHit *vtxhit1 = event->GetVTXHitEntry(mytrk->GetHitIndex(1));
                        MyVTXHit *vtxhit2 = nullptr,*vtxhit3 = nullptr; 
                        if (inum2>=0) vtxhit2 = event->GetVTXHitEntry(mytrk->GetHitIndex(2));
                        else          vtxhit2 = event->GetVTXHitEntry(mytrk->GetHitIndex(3));
                        if (inum3>=0 && inum2>=0) vtxhit3 = event->GetVTXHitEntry(mytrk->GetHitIndex(3));
                        MyVTXHit *vtxhits[4] = {vtxhit0,vtxhit1,vtxhit2,vtxhit3};
                        const int hit_in_graph = (inum3>=0 && inum2>=0) ? 4 : 3; 

                        MyDileptonAnalysis::MyBDTHit newBDTHit = MyDileptonAnalysis::MyBDTHit();
                        newBDTHit.SetReconPt(mytrk->GetReconPT());
                        newBDTHit.SetReconPhi0(mytrk->GetPhi0());
                        newBDTHit.SetReconThe0(mytrk->GetThe0());
                        for (int jlayer = 0; jlayer < hit_in_graph; jlayer++)
                        {
                            if(!vtxhits[jlayer]) std::cout<<"kek"<<std::endl;
                            if(!vtxhits[jlayer]) continue;
                            newBDTHit.SetIsTrue(jlayer, (int) vtxhits[jlayer]->GetSensor() == 0);
                            if(jlayer>1) newBDTHit.SetOutiLayer(jlayer,vtxhits[jlayer]->GetiLayer());
                        }
                        newBDTHit.Setsdphi(0,mytrk->GetsdPhi(0, inum0) * mytrk->GetChargePrime());
                        newBDTHit.Setsdphi(1,mytrk->GetsdPhi(1, inum1) * mytrk->GetChargePrime());
                        if(inum2>=0)newBDTHit.Setsdphi(2,mytrk->GetsdPhi(2, inum2) * mytrk->GetChargePrime());
                        if(inum3>=0)newBDTHit.Setsdphi(3,mytrk->GetsdPhi(3, inum3) * mytrk->GetChargePrime());
                        newBDTHit.Setsdthe(0,mytrk->GetsdThe(0, inum0) * mytrk->GetChargePrime());
                        newBDTHit.Setsdthe(1,mytrk->GetsdThe(1, inum1) * mytrk->GetChargePrime());
                        if(inum2>=0)newBDTHit.Setsdthe(2,mytrk->GetsdThe(2, inum2) * mytrk->GetChargePrime());
                        if(inum3>=0)newBDTHit.Setsdthe(3,mytrk->GetsdThe(3, inum3) * mytrk->GetChargePrime());

                        const float phi_0layer = vtxhits[0]->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                        const float the_0layer = vtxhits[0]->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());

                        float r_first_min[4][2] = {{99,99},{99,99},{99,99},{99,99}};
                        float r_second_min[4][2] = {{99,99},{99,99},{99,99},{99,99}};
                        for (int isecondhit = 0; isecondhit < nvtxhits; isecondhit++)
                        {   
                            int skip_secondhit_search = 0;
                            for (int jlayer = 0; jlayer < hit_in_graph; jlayer++)
                            {
                                if(isecondhit == vtxhits[jlayer]->GetClustId()) skip_secondhit_search = 1;
                            }
                            if (skip_secondhit_search) continue;

                            MyDileptonAnalysis::MyVTXHit *secondvtxhit = event->GetVTXHitEntry(isecondhit);
                            if(secondvtxhit->GetLadder()>49 ) continue;
                            if(secondvtxhit->GetLadder()>24 && secondvtxhit->GetLadder() != itrk + 25 && !recover_fg ) 
                            {   //check on same track cuts is needed
                                MyDileptonAnalysis::MyElectron *second_trk = event->GetEntry(secondvtxhit->GetLadder() - 25);
                                if (!second_trk) {std::cout<<"no second thrack: smt went wrong"<<std::endl;continue;}
                                if (mytrk->GetChargePrime() != second_trk->GetChargePrime()) 
                                {
                                    const float dcenter_phi = ( mytrk->GetCrkphi() - second_trk->GetCrkphi()) / 0.013;
                                    if ( TMath::Abs(dcenter_phi) > 5) 
                                    continue;
                                }
                            }
                            const int secondhit_layer = secondvtxhit->GetLayer(); 

                            const float phi_second_hit = secondvtxhit->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                            const float the_second_hit = secondvtxhit->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());

                            const int secondhit_ilayer = secondvtxhit->GetiLayer(); 
                            const float dphi0 = (dilep_phi_projection[secondhit_ilayer] - phi_second_hit) + 
                                (dilep_phi_projection[secondhit_ilayer<7?secondhit_ilayer+1:6] - dilep_phi_projection[secondhit_ilayer])*
                                (sqrt(SQR(secondvtxhit->GetXHit())+SQR(secondvtxhit->GetYHit()))-radii[secondhit_ilayer])/
                                (radii[secondhit_ilayer<7?secondhit_ilayer+1:6] - radii[secondhit_ilayer]);
                            const float dthe0 = (dilep_the_projection[secondhit_ilayer] - the_second_hit) + 
                                (dilep_the_projection[secondhit_ilayer<7?secondhit_ilayer+1:6] - dilep_the_projection[secondhit_ilayer])*
                                (sqrt(SQR(secondvtxhit->GetXHit())+SQR(secondvtxhit->GetYHit()))-radii[secondhit_ilayer])/
                                (radii[secondhit_ilayer<7?secondhit_ilayer+1:6] - radii[secondhit_ilayer]);
                            if (abs(dphi0) > 0.1 || abs(dthe0) > 0.1) continue;

                            const float sigma_second_phi = mytrk->get_dynamic_sigma_phi_data  (0, secondhit_layer==0 ? 0 : secondhit_layer-1, 0) * mytrk->GetChargePrime();
                            const float sigma_second_the = mytrk->get_dynamic_sigma_theta_data(0, secondhit_layer==0 ? 0 : secondhit_layer-1, 0);// * mytrk->GetChargePrime();

                            const float sdphi_second_hit = (phi_0layer - phi_second_hit) / sigma_second_phi;
                            const float sdthe_second_hit = (the_0layer - the_second_hit) / sigma_second_the;

                            if(sdphi_second_hit>0)
                            {
                                if ( sqrt( SQR(sdphi_second_hit) + SQR(sdthe_second_hit) ) < r_second_min[secondhit_layer][1])
                                {
                                    if ( sqrt( SQR(sdphi_second_hit) + SQR(sdthe_second_hit) ) < r_second_min[secondhit_layer][0])
                                    {
                                        r_second_min[secondhit_layer][0] = sqrt( SQR(sdphi_second_hit) + SQR(sdthe_second_hit) );
                                        newBDTHit.SetSecondHitPhiL(secondhit_layer,0,sdphi_second_hit);
                                        newBDTHit.SetSecondHitTheL(secondhit_layer,0,sdthe_second_hit);
                                    }else
                                    {
                                        r_second_min[secondhit_layer][1] = sqrt( SQR(sdphi_second_hit) + SQR(sdthe_second_hit) );
                                        newBDTHit.SetSecondHitPhiL(secondhit_layer,1,sdphi_second_hit);
                                        newBDTHit.SetSecondHitTheL(secondhit_layer,1,sdthe_second_hit);
                                    }
                                    
                                }
                            }

                            for (int jlayer = 0; jlayer < hit_in_graph; jlayer++)
                            {
                                if(secondhit_layer == vtxhits[jlayer]->GetLayer()) 
                                {
                                    const float phi_loc_layer = vtxhits[jlayer]->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                                    const float the_loc_layer = vtxhits[jlayer]->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                                    const float sdphi_first_hit = ((phi_loc_layer - phi_second_hit) + 
                                        (dilep_phi_projection[secondhit_ilayer<7?secondhit_ilayer+1:6] - dilep_phi_projection[secondhit_ilayer])
                                        *(sqrt(SQR(vtxhits[jlayer]->GetXHit())+SQR(vtxhits[jlayer]->GetYHit()))-sqrt(SQR(secondvtxhit->GetXHit())+SQR(secondvtxhit->GetYHit())))
                                        /(radii[secondhit_ilayer<7?secondhit_ilayer+1:6] - radii[secondhit_ilayer]))/
                                        sigma_second_phi;
                                    const float sdthe_first_hit = ((the_loc_layer - the_second_hit) + 
                                        (dilep_the_projection[secondhit_ilayer<7?secondhit_ilayer+1:6] - dilep_the_projection[secondhit_ilayer])
                                        *(sqrt(SQR(vtxhits[jlayer]->GetXHit())+SQR(vtxhits[jlayer]->GetYHit()))-sqrt(SQR(secondvtxhit->GetXHit())+SQR(secondvtxhit->GetYHit())))
                                        /(radii[secondhit_ilayer<7?secondhit_ilayer+1:6] - radii[secondhit_ilayer]))/
                                        sigma_second_the;
                                    const float r_local = sqrt( SQR(sdphi_first_hit) + SQR(sdthe_first_hit) );

                                    //if ( r_local < 1e-12 ) secondvtxhit->SetLadder(48);///might be un issue
                                    if ( r_local < 1e-12 ) secondvtxhit->SetZHit(4000);
                                    if ( r_local < 1e-12 ) continue;
                                    if ( r_local < r_first_min[secondhit_layer][1] )
                                    {
                                        if ( r_local < r_first_min[secondhit_layer][0])
                                        {
                                            r_first_min[secondhit_layer][0] = r_local;
                                            newBDTHit.SetSecondHitPhiR(secondhit_layer,0,sdphi_first_hit);
                                            newBDTHit.SetSecondHitTheR(secondhit_layer,0,sdthe_first_hit);
                                        }else
                                        {
                                            r_first_min[secondhit_layer][1] = r_local;
                                            newBDTHit.SetSecondHitPhiR(secondhit_layer,1,sdphi_first_hit);
                                            newBDTHit.SetSecondHitTheR(secondhit_layer,1,sdthe_first_hit);
                                        }

                                    }
                                }
                            }

                        }
                        const double ecore  = newBDTrack.GetEcore()*( 1.+1./TMath::Tan(the0_trk_proj)/TMath::Tan(the0_trk_proj));
                        const double ecore1 = ecore - newBDTHit.GetReconPt();
                        const double ecore2 = ecore1/newBDTrack.GetEcore();////need to be replased with ecore
                        const double reconpt1 = newBDTrack.GetPt() - newBDTHit.GetReconPt();
                        const double reconpt2 = reconpt1/newBDTrack.GetPt();
                        if (newBDTrack.GetPt()<0.4) newBDTrack.SetPt(0.405);
                        const double BDTHitInput[24] = {//['Pt', 'Ecore', 'Centrality', 'reconpt', 'sdthe0', 'SecondHitPhiR0', 'SecondHitTheR0', 'SecondHitPhiR01', 'SecondHitTheR01', 'sdthe1', 'SecondHitPhiR1', 'SecondHitTheR1', 'SecondHitPhiR11', 'SecondHitTheR11', 'sdthe2', 'SecondHitPhiR2', 'SecondHitTheR2', 'SecondHitPhiR21', 'SecondHitTheR21', 'sdthe3', 'SecondHitPhiR3', 'SecondHitTheR3', 'SecondHitPhiR31', 'SecondHitTheR31']
                            newBDTrack.GetPt(), ecore2, (double) newBDTrack.GetCentrality(),  reconpt2,
                            newBDTHit.Getsdthe(0), newBDTHit.GetSecondHitPhiR(0), newBDTHit.GetSecondHitTheR(0), newBDTHit.GetSecondHitPhiR(0,1), newBDTHit.GetSecondHitTheR(0,1),
                            newBDTHit.Getsdthe(1), newBDTHit.GetSecondHitPhiR(1), newBDTHit.GetSecondHitTheR(1), newBDTHit.GetSecondHitPhiR(1,1), newBDTHit.GetSecondHitTheR(1,1),
                            newBDTHit.Getsdthe(2), newBDTHit.GetSecondHitPhiR(2), newBDTHit.GetSecondHitTheR(2), newBDTHit.GetSecondHitPhiR(2,1), newBDTHit.GetSecondHitTheR(2,1),
                            newBDTHit.Getsdthe(3), newBDTHit.GetSecondHitPhiR(3), newBDTHit.GetSecondHitTheR(3), newBDTHit.GetSecondHitPhiR(3,1), newBDTHit.GetSecondHitTheR(3,1)   
                        };
                        const double BDTConvInput[38] = {//['Pt', 'Ecore', 'Centrality', 'Charge', 'reconpt', 'sdphi0', 'SecondHitPhiR0', 'SecondHitPhiL0', 'SecondHitTheR0', 'SecondHitTheL0', 'SecondHitPhiR01', 'SecondHitPhiL01', 'SecondHitTheR01', 'SecondHitTheL01', 'SecondHitPhiR1', 'SecondHitPhiL1', 'SecondHitTheR1', 'SecondHitTheL1', 'SecondHitPhiR11', 'SecondHitPhiL11', 'SecondHitTheR11', 'SecondHitTheL11', 'SecondHitPhiR2', 'SecondHitPhiL2', 'SecondHitTheR2', 'SecondHitTheL2', 'SecondHitPhiR21', 'SecondHitPhiL21', 'SecondHitTheR21', 'SecondHitTheL21', 'SecondHitPhiR3', 'SecondHitPhiL3', 'SecondHitTheR3', 'SecondHitTheL3', 'SecondHitPhiR31', 'SecondHitPhiL31', 'SecondHitTheR31', 'SecondHitTheL31']
                            newBDTrack.GetPt(), ecore2, (double) newBDTrack.GetCentrality(), (double)newBDTrack.GetCharge(),
                            reconpt2, newBDTHit.Getsdphi(0)-newBDTHit.Getsdphi(1),
                            newBDTHit.GetSecondHitPhiR(0), newBDTHit.GetSecondHitPhiL(0), newBDTHit.GetSecondHitTheR(0), newBDTHit.GetSecondHitTheL(0),newBDTHit.GetSecondHitPhiR(0,1), newBDTHit.GetSecondHitPhiL(0,1), newBDTHit.GetSecondHitTheR(0,1), newBDTHit.GetSecondHitTheL(0,1),
                            newBDTHit.GetSecondHitPhiR(1), newBDTHit.GetSecondHitPhiL(1), newBDTHit.GetSecondHitTheR(1), newBDTHit.GetSecondHitTheL(1),newBDTHit.GetSecondHitPhiR(1,1), newBDTHit.GetSecondHitPhiL(1,1), newBDTHit.GetSecondHitTheR(1,1), newBDTHit.GetSecondHitTheL(1,1),
                            newBDTHit.GetSecondHitPhiR(2), newBDTHit.GetSecondHitPhiL(2), newBDTHit.GetSecondHitTheR(2), newBDTHit.GetSecondHitTheL(2),newBDTHit.GetSecondHitPhiR(2,1), newBDTHit.GetSecondHitPhiL(2,1), newBDTHit.GetSecondHitTheR(2,1), newBDTHit.GetSecondHitTheL(2,1),
                            newBDTHit.GetSecondHitPhiR(3), newBDTHit.GetSecondHitPhiL(3), newBDTHit.GetSecondHitTheR(3), newBDTHit.GetSecondHitTheL(3),newBDTHit.GetSecondHitPhiR(3,1), newBDTHit.GetSecondHitPhiL(3,1), newBDTHit.GetSecondHitTheR(3,1), newBDTHit.GetSecondHitTheL(3,1)       
                        };
                        const double probsHIT[4] = {0.7048700405814634, 0.80, 0.8971706093225199, 0.95};
                        //const double probsConv[4] = {0.0227005061149391,0.03696777185169043, 0.06905145040745402, 0.13901432529026608};
                        //const double probsConv[4] = {0.16655781590812105,0.20684925395539883, 0.2773192065090052, 0.3935912730094734};
                        const double probsConv[4] = {0.059511564980906206,0.0672349245703329, 0.07958121918163912, 0.10254988100220457};
                        float kek = 1.;
                        if ( mytrk->GetPtPrime()<0.4 && !(newBDTHit.Getsdthe(2)>-5 && newBDTHit.Getsdthe(3)>-5)) kek = 1000;
                        const double BDTHIT_prob = MyML::GetHitBDTProb(BDTHitInput)/kek;
                        const double BDTConv_prob = MyML::GetConvBDTProb(BDTConvInput);
                        //const double check_conv[38] = {1.1865234375, 0.03924144270263441, 8.0, -1.0, 0.00905609130859375, -0.30126953125, -43.8125, -99.0, 18.1875, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -21.0625, 45.40625, 7.8828125, 14.8203125, -21.0625, 61.25, 19.09375, 18.296875, -25.203125, -99.0, 20.765625, -99.0, -48.875, -99.0, -8.578125, -99.0};
                        //const double check_hit[24] = {0.6347407102584839, -0.06209205680455947, 81.0, -0.11198512464761734, 0.869284451007843, -99.0, -99.0, -99.0, -99.0, 0.5129899382591248, -99.0, -99.0, -99.0, -99.0, 1.3870052099227905, -99.0, -99.0, -99.0, -99.0, -10.0, -99.0, -99.0, -99.0, -99.0}; 
                        //std::cout<<MyML::GetHitBDTProb(check_hit)<<" "<<BDTHIT_prob<<" "<<MyML::GetConvBDTProb(check_conv)<<" "<<BDTConv_prob<<" "<<newBDTHit.GetIsTrue(0)<<" "<<newBDTHit.GetIsTrue(1)<<" "<<newBDTHit.GetIsTrue(2)<<" "<<newBDTHit.GetIsTrue(3)<<std::endl;
                        //std::cout<<BDTHitInput[0]<<" "<<BDTHitInput[1]<<" "<<BDTHitInput[2]<<" "<<BDTHitInput[3]<<" "<<BDTHitInput[4]<<" "<<BDTHitInput[5]<<" "<<std::endl;
                        newBDTHit.SetIsTrue(2, (int) 1000*BDTHIT_prob);
                        newBDTHit.SetIsTrue(3, (int) 1000*BDTConv_prob);
                        newBDTHit.SetIsTrue(0, MyML::GetHitBDT(BDTHIT_prob,probsHIT));
                        newBDTHit.SetIsTrue(1, MyML::GetConvBDT(BDTConv_prob,probsConv));
                        newBDTrack.AddBDTHit(&newBDTHit);
                        chi2 = (1-BDTHIT_prob)*3.0/(1-probsHIT[1]);
                        if(chi2<min_chi2) 
                        {
                            if (recover_fg) mytrk->SetNHits(newBDTHit.GetIsTrue(0));
                            if (recover_fg) mytrk->SetTOFDPHI(newBDTHit.GetIsTrue(1) );
                            if (!recover_fg) mytrk->SetPC3SDZ(newBDTHit.GetIsTrue(0));
                            if (!recover_fg) mytrk->SetPC3SDPHI(newBDTHit.GetIsTrue(1));
                            newBDTrack.SetReconPt(mytrk->GetReconPT());
                            newBDTrack.SetReconPhi0(mytrk->GetPhi0());
                            newBDTrack.SetReconThe0(mytrk->GetThe0());
                        }
                        mytrk->SetGhost(mytrk->GetGhost()+(probsConv[1]>BDTConv_prob?30:0));
                    }
                    //chi2 = TMath::Abs(recon_pt-pt)/pt*30/(2+(int)(inum2>=0)+(int)(inum3>=0));
                    if(chi2<min_chi2) {min_chi2=chi2;final_number=numbers[0][inum];} 
                    
                    if (is_fill_hsits&&numbers[0].size()<10) chi2_ndf[central_bin]->Fill(chi2, numbers[0].size(), pt);
                }
            }
            if(is_fill_hsits) chi2_ndf[central_bin]->Fill(min_chi2, 19, pt);

            if(min_chi2<800000)
            {
                mytrk->SetHitCounter(3,0);mytrk->SetHitCounter(2,0);
                mytrk->SetGhost( (int) (mytrk->GetGhost()*1./newBDTrack.GetNBDThit()));
                BDTracklist.push_back(newBDTrack);

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
                    vtxhits[i]->SetLadder(25+itrk);
                    if(vtxhits[i]->GetSensor() == 0) istruehitcounter++;
                }        
                if(is_fill_hsits) truehithist->Fill(istruehitcounter,mytrk->GetPtPrime(),event->GetCentrality());
                if(is_fill_hsits) chi2_ndf[central_bin]->Fill(min_chi2, 10+istruehitcounter, pt);
                if(is_fill_hsits&&true) 
                                        truehitsigmahist->Fill(istruehitcounter+(sigma-2)*10,mytrk->GetPtPrime(),event->GetCentrality());
                if(is_fill_hsits)
                {
                    charge_recover_hist->Fill(mytrk->GetCharge(),8,pt);
                    charge_recover_hist->Fill(mytrk->GetChargePrime(),9,pt);
                    charge_recover_hist->Fill(mytrk->GetPhiDC()>mytrk->GetPhi0()?1:-1,10,pt);
                    charge_recover_hist->Fill(mytrk->GetMcId()<3?1:-1,11,pt);
                }
                if(is_fill_hsits) 
                {
                    float hit_counter_jlayer_inside[total_vtx_layers*5] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                    int passed_association_in_sigma[3][nvtx_layers] = {{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0}};
                    for (int jlayer = 0; jlayer < total_vtx_layers+2; jlayer++)
                    {
                        for (unsigned int jentry = 0; jentry < hit_counter_jlayer[jlayer].size(); jentry++)
                        {
                            for (unsigned int kentry = jentry+1; kentry < hit_counter_jlayer[jlayer].size(); kentry++)
                            {
                                if ( TMath::Abs( hit_counter_jlayer[jlayer][jentry][2] - hit_counter_jlayer[jlayer][kentry][2] ) < 0.001
                                  && TMath::Abs( hit_counter_jlayer[jlayer][jentry][3] - hit_counter_jlayer[jlayer][kentry][3] ) < 0.001 )
                                    hit_counter_jlayer[jlayer][kentry][0] = 100;
                            } 
                        }                        
                    }
                    for (int jlayer = 0; jlayer < nvtx_layers; jlayer++)
                    {
                        int iarraylayer = jlayer;
                        if(jlayer>1) iarraylayer = jlayer+6;
                        for (unsigned int jentry = 0; jentry < hit_counter_jlayer[iarraylayer].size(); jentry++)
                        {
                            if ( TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][0])<4 && TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][1])<4 ) 
                                passed_association_in_sigma[0][jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][0])<3 && TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][1])<3 ) 
                                passed_association_in_sigma[1][jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][0])<2 && TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][1])<2 ) 
                                passed_association_in_sigma[2][jlayer]++;
                        }                        
                    }
                    for (int jlayer = 0; jlayer < total_vtx_layers+2; jlayer++)
                    {
                        for (unsigned int jentry = 0; jentry < hit_counter_jlayer[jlayer].size(); jentry++)
                        {
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<5 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<5) 
                                hit_counter_jlayer_inside[jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<4 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<4 
                            && passed_association_in_sigma[0][0] && passed_association_in_sigma[0][1] && 
                            (passed_association_in_sigma[0][2] || passed_association_in_sigma[0][3] ) ) hit_counter_jlayer_inside[10+jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<3 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<3 
                            && passed_association_in_sigma[1][0] && passed_association_in_sigma[1][1] && 
                            (passed_association_in_sigma[1][2] || passed_association_in_sigma[1][3] ) ) hit_counter_jlayer_inside[20+jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<2 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<2 
                            && passed_association_in_sigma[2][0] && passed_association_in_sigma[2][1] && 
                            (passed_association_in_sigma[2][2] || passed_association_in_sigma[2][3] ) )  hit_counter_jlayer_inside[30+jlayer]++;
                        }                        
                    }
                    for (int jlayer = 0; jlayer < total_vtx_layers*5; jlayer++)
                    {
                        ilayerhitshist[central_bin]->Fill(hit_counter_jlayer_inside[jlayer], jlayer, pt);
                    }
                }

            }else{
                mytrk->SetHitCounter(0,0);
                if(is_fill_hsits)
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
        const int rungroup = event->GetRunGroup();
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

                const float dphi0 = (dilep_phi_projection[ilayer] - phi_hit);
                const float dthe0 = (dilep_the_projection[ilayer] - theta_hit);
                const float dphi = dphi0 +
                    (dilep_phi_projection[ilayer<7?ilayer+1:6] - dilep_phi_projection[ilayer])*
                    (sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/
                    (radii[ilayer<7?ilayer+1:6] - radii[ilayer]);
                const float dthe = dthe0 + 
                    (dilep_the_projection[ilayer<7?ilayer+1:6] - dilep_the_projection[ilayer])*
                    (sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/
                    (radii[ilayer<7?ilayer+1:6] - radii[ilayer]);

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
                    if(false)
                    {
                        myvtx_hist[central_bin]->Fill(event->GetPreciseX()-vtxhit->GetXHit()+sqrt(SQR(vtxhit->GetXHit()-event->GetPreciseX())
                                        +SQR(vtxhit->GetYHit()-event->GetPreciseY()))*TMath::Cos(mytrk->GetPhi0Prime()),event->GetRunGroup(),0.5);
                        myvtx_hist[central_bin]->Fill(event->GetPreciseY()-vtxhit->GetYHit()+sqrt(SQR(vtxhit->GetXHit()-event->GetPreciseX())
                                        +SQR(vtxhit->GetYHit()-event->GetPreciseY()))*TMath::Sin(mytrk->GetPhi0Prime()),event->GetRunGroup(),1.5);
                        myvtx_hist[central_bin]->Fill(event->GetPreciseZ()-vtxhit->GetZHit()+sqrt(SQR(vtxhit->GetXHit()-event->GetPreciseX())
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
                    if (pt<1.5) continue;
                    const float dphi0 = dphi + mytrk->GetPhi0() - mytrk->GetPhi0Prime();
                    dphi_phi0_init_hist[layer]->Fill(dphi0, mytrk->GetPhi0(), 2*mytrk->GetArm() + charge_bin);
                    dphi_phi0_corr_hist[layer]->Fill(dphi, mytrk->GetPhi0Prime(), 2*mytrk->GetArm() + charge_bin);
                    dphi_the0_init_hist[layer]->Fill(dphi0, mytrk->GetThe0Prime(), 2*mytrk->GetArm() + charge_bin);
                    dphi_the0_corr_hist[layer]->Fill(dphi, mytrk->GetThe0Prime(), 2*mytrk->GetArm() + charge_bin);
                    dthe_phi0_init_hist[layer]->Fill(dphi0, mytrk->GetPhiDC(), 2*mytrk->GetArm() + charge_bin);
                    dthe_phi0_corr_hist[layer]->Fill(dphi, mytrk->GetPhiDC(), 2*mytrk->GetArm() + charge_bin);
                }
                if (abs(sdphi) < sigma && is_fill_hadron_hsits)
                {
                    dthe_hist[central_bin] ->Fill( dthe, hist_2nd_arg, pt);
                    sdthe_hist[central_bin]->Fill(sdthe, hist_2nd_arg, pt);
                    if (pt<1.5) continue;
                    const float newthe0 = mytrk->GetThe0() - ((event->GetVtxZ() - event->GetPreciseZ()) / 220) * TMath::Sin(mytrk->GetThe0());
                    const float dthe0 = dthe + newthe0 - mytrk->GetThe0Prime();
                    dthe_the0_init_hist[layer]->Fill(dthe0, newthe0, 2*mytrk->GetArm() + charge_bin);
                    dthe_the0_corr_hist[layer]->Fill(dthe, mytrk->GetThe0Prime(), 2*mytrk->GetArm() + charge_bin);
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

    void MyEvent::SetDCA2(const unsigned int itrk, const int layer3, const int type)
    {
        MyDileptonAnalysis::MyElectron *mytrk;
        if(type) mytrk = this->GetElecCand(itrk);
        else mytrk = this->GetEntry(itrk);
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
            if (pt < 0.4)
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

    void MyEventContainer::FillTrueDCA(const float weight)
    {
        const int central_bin = (int)event->GetCentrality() / 20;
        const int nleptons = event->GetNtrack();
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            if (mytrk->GetNHits()<10) continue;
            const float pt = mytrk->GetPtPrime();
            DCA12_hist[central_bin]->Fill(mytrk->GetDCAX2(),mytrk->GetDCAY2(),pt);
            const float eConv = std::log10(mytrk->GetNHits()+1)*5 + std::log10(mytrk->GetTOFDPHI()+1);
            if (mytrk->GetChargePrime()>0)
                DCA2_hist[central_bin]->Fill(mytrk->GetDCA2(),pt,eConv,weight);
            else
                sDCA2_hist[central_bin]->Fill(mytrk->GetDCA2(),pt,eConv,weight);

            DCPT_ReconPT->Fill(mytrk->GetReconPT(),pt,event->GetCentrality()+100*( mytrk->GetChargePrime() > 0 ? 0 : 1));


            int hit_assocaition = 0;
            if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<3) ||
                   (TMath::Abs(mytrk->GetMinsDphi(2))<3) ) && 
                   (TMath::Abs(mytrk->GetMinsDphi(1))<3) && 
                   (mytrk->GetMinsDphi(0))>-3 ) ) hit_assocaition=1;
            if (hit_assocaition==1 && mytrk->GetMinsDphi(0)<3 ) hit_assocaition=2;
            if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<2 && TMath::Abs(mytrk->GetMinsDthe(3))<2) ||
                   (TMath::Abs(mytrk->GetMinsDphi(2))<2 && TMath::Abs(mytrk->GetMinsDthe(2))<2) ) && 
                   (TMath::Abs(mytrk->GetMinsDphi(1))<2) && 
                   (mytrk->GetMinsDphi(0)> -1 ) )) hit_assocaition=3;
            if (hit_assocaition==3 && mytrk->GetMinsDphi(0)<2 ) hit_assocaition=4;
            int conv_reject = 0;
            if ( ((int)mytrk->GetEmcdphi_e())%100==0) conv_reject=5;
            if ( ((int)mytrk->GetEmcdphi_e())%100<3 && ((int)mytrk->GetEmcdphi_e())/100<3) conv_reject=10;
            if ( ((int)mytrk->GetEmcdphi_e())%100<1 && ((int)mytrk->GetEmcdphi_e())/100<3) conv_reject=15;
            if ( ((int)mytrk->GetEmcdphi_e())%100<1 && ((int)mytrk->GetEmcdphi_e())/100<1) conv_reject=20;

            const int hist_in = hit_assocaition + conv_reject;

            if (mytrk->GetChargePrime()>0)
                DCA_2D_hist[central_bin] ->Fill(mytrk->GetDCA2(),pt,hist_in,weight);
            else
                sDCA_2D_hist[central_bin]->Fill(mytrk->GetDCA2(),pt,hist_in,weight);
            
        }
    }

    void MyEventContainer::FillVTXAcceptance()
    {
        const int nvtx_hits = event->GetNVTXhit();
        for (int ihit = 0; ihit < nvtx_hits; ihit++)
        {
            MyDileptonAnalysis::MyVTXHit *vtxhit = event->GetVTXHitEntry(ihit);
            const float phi = vtxhit->GetPhiHit(0,0,0);
            const float zhit = vtxhit->GetZHit();
            const int layer = vtxhit->GetiLayer();

            vtx_accaptance_hist->Fill(phi,zhit,layer);
            if (vtxhit->GetSensor() > 1) vtx_deadmaps_hist->Fill(phi,zhit,layer);
        }
    }

    void MyEventContainer::FillFlow(const float psi_BBCS, const float psi_BBCN, const float psi_FVTXS, const float psi_FVTXN)
    {
        const int centrality = event->GetCentrality();
        const float psi_BBC = event->GetPsi2BBC();
        const float psi_FVTX = event->GetPsi2FVTXA0();
        if (psi_BBCS <-99||psi_BBCN <-99||psi_BBC <-99||
            psi_FVTXS<-99||psi_FVTXN<-99||psi_FVTX<-99) return;
        BBC_psi_hist->Fill(psi_BBCS,   centrality, 0);
        BBC_psi_hist->Fill(psi_BBCN,   centrality, 1);
        BBC_psi_hist->Fill(psi_BBC,    centrality, 2);
        FVTX_psi_hist->Fill(psi_FVTXS, centrality, 0);
        FVTX_psi_hist->Fill(psi_FVTXN, centrality, 1);
        FVTX_psi_hist->Fill(psi_FVTX,  centrality, 2);
        const float cos_BBC_FVTXS   = cos( 2 * ( psi_BBC - psi_FVTXS ) );
        const float cos_BBC_FVTXN   = cos( 2 * ( psi_BBC - psi_FVTXN ) );
        const float cos_FVTXS_FVTXN = cos( 2 * ( psi_FVTXS - psi_FVTXN ) );
        const float cos_FVTX_BBCS = cos( 2 * ( psi_FVTX - psi_BBCS ) );
        const float cos_FVTX_BBCN = cos( 2 * ( psi_FVTX - psi_BBCN ) );
        const float cos_BBCS_BBCN = cos( 2 * ( psi_BBCS - psi_BBCN ) );
        cos_BBC_hist->Fill(cos_BBC_FVTXS,   centrality, 0);
        cos_BBC_hist->Fill(cos_BBC_FVTXN,   centrality, 1);
        cos_BBC_hist->Fill(cos_FVTXS_FVTXN, centrality, 2);
        cos_FVTX_hist->Fill(cos_FVTX_BBCS,   centrality, 0);
        cos_FVTX_hist->Fill(cos_FVTX_BBCN,   centrality, 1);
        cos_FVTX_hist->Fill(cos_BBCS_BBCN,   centrality, 2);
        
        const int nleptons = event->GetNtrack();
        for (int itrk = 0; itrk < nleptons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            if ( mytrk->GetGhost()>15 || mytrk->GetIsConv()>=1 ) continue;
            if ( mytrk->GetN0()<2 || mytrk->GetDisp()>5 || mytrk->GetChi2()/(mytrk->GetNpe0()+0.1)>10 || TMath::Abs(mytrk->GetEmcTOF())>5)  continue;
            if (mytrk->GetHitCounter(0) < 1) continue;
            const float pt = mytrk->GetPtPrime();
            const float e_phi = mytrk->GetPhi0() < TMath::Pi() ? mytrk->GetPhi0() : mytrk->GetPhi0() - TMath::Pi();
            if(TMath::Abs(e_phi)>TMath::Pi()||TMath::Abs(psi_BBC)>TMath::Pi()||TMath::Abs(psi_FVTX)>TMath::Pi() ) 
                std::cout<<e_phi<<" "<<psi_BBC<<" "<<psi_FVTX<<std::endl;
            const float v2_bbc  = cos ( 2 * ( e_phi - psi_BBC  ) );
            const float v2_fvtx = cos ( 2 * ( e_phi - psi_FVTX ) );
            v2_BBC_hist ->Fill(v2_bbc, pt,centrality+50*(1-mytrk->GetChargePrime()));          
            v2_FVTX_hist->Fill(v2_fvtx,pt,centrality+50*(1-mytrk->GetChargePrime()));
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

    void MyEventContainer::fill_inv_mass(const float weight, const unsigned int pool_depth)
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
                const float dcenter_phi = (newTrack1->GetCrkphi() - newTrack2->GetCrkphi()) / 0.013;
                //const float dcenter_zed = (newTrack1->GetCrkz() - newTrack2->GetCrkz())/5.;
                if ( (TMath::Abs(dcenter_phi)<5|| newTrack1->GetCrkphi()<-99 || newTrack2->GetCrkphi()<-99 ) ) continue;//|| sqrt(SQR(dcenter_zed)+SQR(dcenter_phi)) < 10 
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

                const float dca = sqrt( SQR(newTrack1->GetDCAX2() - newTrack2->GetDCAX2()) + SQR(newTrack1->GetDCAY2() - newTrack2->GetDCAY2()) );

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
                inv_mass_dca_fg0[in_hist]->Fill(dca, invm, pair_pt, weight);
                delt_phi_dca_fg0[in_hist]->Fill(dca, dphi, pair_pt, weight);
                  
                MyDileptonAnalysis::MyVTXHit *vtxhit10 = event->GetVTXHitEntry(newTrack1->GetHitIndex(0));
                MyDileptonAnalysis::MyVTXHit *vtxhit11 = event->GetVTXHitEntry(newTrack1->GetHitIndex(1));
                MyDileptonAnalysis::MyVTXHit *vtxhit12 = nullptr; 
                if (newTrack1->GetHitCounter(2)>0) vtxhit12 = event->GetVTXHitEntry(newTrack1->GetHitIndex(2));
                else                               vtxhit12 = event->GetVTXHitEntry(newTrack1->GetHitIndex(3));
                MyDileptonAnalysis::MyVTXHit *vtxhit20 = event->GetVTXHitEntry(newTrack2->GetHitIndex(0));
                MyDileptonAnalysis::MyVTXHit *vtxhit21 = event->GetVTXHitEntry(newTrack2->GetHitIndex(1));
                MyDileptonAnalysis::MyVTXHit *vtxhit22 = nullptr; 
                if (newTrack2->GetHitCounter(2)>0) vtxhit22 = event->GetVTXHitEntry(newTrack2->GetHitIndex(2));
                else                               vtxhit22 = event->GetVTXHitEntry(newTrack2->GetHitIndex(3));

                const float phi11 = vtxhit10->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float phi12 = vtxhit11->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float phi13 = vtxhit12->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float the11 = vtxhit10->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float the12 = vtxhit11->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float the13 = vtxhit12->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float phi21 = vtxhit20->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float phi22 = vtxhit21->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float phi23 = vtxhit22->GetPhiHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float the21 = vtxhit20->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float the22 = vtxhit21->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());
                const float the23 = vtxhit22->GetTheHit(event->GetPreciseX(),event->GetPreciseY(),event->GetPreciseZ());

                if (!(newTrack1->GetMcId()>900&&newTrack2->GetMcId()>900))
                    continue;
                if (!(newTrack1->GetTOFDPHI()>900&&newTrack2->GetTOFDPHI()>900))
                    continue;
                if (event->GetCentrality()<20 && !(newTrack1->GetMcId()>9000&&newTrack2->GetMcId()>9000))
                    continue;
                inv_mass_dca_fg1[in_hist]->Fill(dca, invm, pair_pt, weight);    
                delt_phi_dca_fg1[in_hist]->Fill(dca, dphi, pair_pt, weight);

                if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.02 ) continue;
                if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.01 ) continue;
                if ( TMath::Abs(phi13-phi23)<0.001 && TMath::Abs(the13-the23)<0.01 ) continue;

                inv_mass_dca_fg2[in_hist]->Fill(dca, invm, pair_pt, weight);
                delt_phi_dca_fg2[in_hist]->Fill(dca, dphi, pair_pt, weight);
                if (!(newTrack1->GetTOFDPHI()>9000&&newTrack2->GetTOFDPHI()>9000))
                    continue;
                inv_mass_dca_fg3[in_hist]->Fill(dca, invm, pair_pt, weight);
                delt_phi_dca_fg3[in_hist]->Fill(dca, dphi, pair_pt, weight);
                if (!(newTrack1->GetNHits()>9000&&newTrack2->GetNHits()>9000))
                    continue;
                inv_mass_dca_fg4[in_hist]->Fill(dca, invm, pair_pt, weight);
                delt_phi_dca_fg4[in_hist]->Fill(dca, dphi, pair_pt, weight);

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
                    if( TMath::Abs ( event->GetPreciseZ() - evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetPreciseZ() ) > 4 ) continue;
                    const float psi_pip = event->GetPsi2BBC();
                    const float psi_pim = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetPsi2BBC();
                    if ( TMath::Abs(psi_pip-psi_pim)>TMath::Pi()/8 && TMath::Abs(psi_pip-psi_pim)<7*TMath::Pi()/8 )
                        continue;
                    //if (!(((newTrack2->GetMinsDphi(2)>0||newTrack2->GetHitCounter(2)<1)&&(newTrack2->GetMinsDphi(3)>0||newTrack2->GetHitCounter(3)<1)&&newTrack2->GetMinsDphi(0)>0)||newTrack2->GetGhost()<10))
                    //    continue;
                    // const float a2 = newTrack2->GetMinDist(0);
                    // const float b2 = newTrack2->GetMinDist(1);
                    // const float c2 = newTrack2->GetMinDist(2);
                    ////////pair cuts
                    //const float phi_pip = newTrack1->GetChargePrime() >  newTrack2->GetChargePrime() ? newTrack1->GetPhiDC() : newTrack2->GetPhiDC();
                    //const float phi_pim = newTrack1->GetChargePrime() <= newTrack2->GetChargePrime() ? newTrack1->GetPhiDC() : newTrack2->GetPhiDC();
                    const float dcenter_phi = (newTrack1->GetCrkphi() - newTrack2->GetCrkphi()) / 0.013;
                    //const float dcenter_zed = (newTrack1->GetCrkz() - newTrack2->GetCrkz())/5.;
                    if ( (TMath::Abs(dcenter_phi)<5|| newTrack1->GetCrkphi()<-99 || newTrack2->GetCrkphi()<-99 ) ) continue;//|| sqrt(SQR(dcenter_zed)+SQR(dcenter_phi)) < 10 
                    const float dalpha = newTrack1->GetAlpha() - newTrack2->GetAlpha();
                    const float dphiDC = newTrack1->GetPhiDC() - newTrack2->GetPhiDC();
                    const float dzed = newTrack1->GetZDC() - newTrack2->GetZDC();
                    if (TMath::Abs(dzed) < 6.0 && TMath::Abs(dphiDC - (0.13 * dalpha)) < 0.015) continue;
                    if (TMath::Abs(dphiDC - (0.04 * dalpha)) < 0.015) continue;
                    if (TMath::Abs(dphiDC - (-0.065 * dalpha)) < 0.015) continue;
                    //////////////end of pair cuts
                    const int in_hist = (int) (event->GetCentrality() / 20) + N_centr*((newTrack1->GetChargePrime()+newTrack2->GetChargePrime()+2)/2);

                    const float dca = sqrt( SQR(newTrack1->GetDCAX2() - newTrack2->GetDCAX2()) + SQR(newTrack1->GetDCAY2() - newTrack2->GetDCAY2()) );

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
                    inv_mass_dca_bg0[in_hist]->Fill(dca, invm, pair_pt, weight);
                    delt_phi_dca_bg0[in_hist]->Fill(dca, dphi, pair_pt, weight);
                    
                    MyDileptonAnalysis::MyVTXHit *vtxhit10 = event->GetVTXHitEntry(newTrack1->GetHitIndex(0));
                    MyDileptonAnalysis::MyVTXHit *vtxhit11 = event->GetVTXHitEntry(newTrack1->GetHitIndex(1));
                    MyDileptonAnalysis::MyVTXHit *vtxhit12 = nullptr; 
                    if (newTrack1->GetHitCounter(2)>0) vtxhit12 = event->GetVTXHitEntry(newTrack1->GetHitIndex(2));
                    else                               vtxhit12 = event->GetVTXHitEntry(newTrack1->GetHitIndex(3));
                    MyDileptonAnalysis::MyVTXHit *vtxhit20 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetVTXHitEntry(newTrack2->GetHitIndex(0));
                    MyDileptonAnalysis::MyVTXHit *vtxhit21 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetVTXHitEntry(newTrack2->GetHitIndex(1));
                    MyDileptonAnalysis::MyVTXHit *vtxhit22 = nullptr; 
                    if (newTrack2->GetHitCounter(2)>0) vtxhit22 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetVTXHitEntry(newTrack2->GetHitIndex(2));
                    else                               vtxhit22 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetVTXHitEntry(newTrack2->GetHitIndex(3));

                    const float ex1 = event->GetPreciseX();
                    const float ey1 = event->GetPreciseY();
                    const float ez1 = event->GetPreciseZ(); 
                    const float ex2 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetPreciseX();
                    const float ey2 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetPreciseY();
                    const float ez2 = evtbuff_list[icent_mix][izvtx_mix][ipsi2_mix][ievent].GetPreciseZ(); 

                    const float phi11 = vtxhit10->GetPhiHit(ex1,ey1,ez1);
                    const float phi12 = vtxhit11->GetPhiHit(ex1,ey1,ez1);
                    const float phi13 = vtxhit12->GetPhiHit(ex1,ey1,ez1);
                    const float the11 = vtxhit10->GetTheHit(ex1,ey1,ez1);
                    const float the12 = vtxhit11->GetTheHit(ex1,ey1,ez1);
                    const float the13 = vtxhit12->GetTheHit(ex1,ey1,ez1);
                    const float phi21 = vtxhit20->GetPhiHit(ex2,ey2,ez2);
                    const float phi22 = vtxhit21->GetPhiHit(ex2,ey2,ez2);
                    const float phi23 = vtxhit22->GetPhiHit(ex2,ey2,ez2);
                    const float the21 = vtxhit20->GetTheHit(ex2,ey2,ez2);
                    const float the22 = vtxhit21->GetTheHit(ex2,ey2,ez2);
                    const float the23 = vtxhit22->GetTheHit(ex2,ey2,ez2);

                    if (!(newTrack1->GetMcId()>900&&newTrack2->GetMcId()>900))
                        continue;
                    if (!(newTrack1->GetTOFDPHI()>900&&newTrack2->GetTOFDPHI()>900))
                        continue;
                    if (event->GetCentrality()<20 && !(newTrack1->GetMcId()>9000&&newTrack2->GetMcId()>9000))
                        continue;
                    inv_mass_dca_bg1[in_hist]->Fill(dca, invm, pair_pt, weight);    
                    delt_phi_dca_bg1[in_hist]->Fill(dca, dphi, pair_pt, weight);
    
                    if ( TMath::Abs(phi11-phi21)<0.002 && TMath::Abs(the11-the21)<0.02 ) continue;
                    if ( TMath::Abs(phi12-phi22)<0.001 && TMath::Abs(the12-the22)<0.01 ) continue;
                    if ( TMath::Abs(phi13-phi23)<0.001 && TMath::Abs(the13-the23)<0.01 ) continue;
    
                    inv_mass_dca_bg2[in_hist]->Fill(dca, invm, pair_pt, weight);
                    delt_phi_dca_bg2[in_hist]->Fill(dca, dphi, pair_pt, weight);
                    if (!(newTrack1->GetTOFDPHI()>9000&&newTrack2->GetTOFDPHI()>9000))
                        continue;
                    inv_mass_dca_bg3[in_hist]->Fill(dca, invm, pair_pt, weight);
                    delt_phi_dca_bg3[in_hist]->Fill(dca, dphi, pair_pt, weight);
                    if (!(newTrack1->GetNHits()>9000&&newTrack2->GetNHits()>9000))
                        continue;
                    inv_mass_dca_bg4[in_hist]->Fill(dca, invm, pair_pt, weight);
                    delt_phi_dca_bg4[in_hist]->Fill(dca, dphi, pair_pt, weight);
                }
            }
        }
        if(event->GetNtrack()>0)this->fill_evtbuff_list(pool_depth);
    }

    void MyEventContainer::fill_inv_mass_sim(const float weight)
    {
        for (int ielectron = 0; ielectron < event->GetNgentrack(); ielectron++)
        {
            MyDileptonAnalysis::MyGenTrack* newTrack1 = event->GetGenTrack(ielectron);

            if (newTrack1->GetID() == 22) continue;
            for (int jelectron = ielectron+1; jelectron < event->GetNgentrack(); jelectron++)
            {
                MyDileptonAnalysis::MyGenTrack *newTrack2 = event->GetGenTrack(jelectron);
                if (newTrack2->GetID() == 22) continue;
                
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

                inv_mass_dca_gen[in_hist]->Fill(dca, invm, pair_pt, weight);
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
                    n_good_el++;
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
            MyDileptonAnalysis::MyElectron *electron = event->GetEntry(i);
            const float alpha_offset = - (event->GetPreciseX() / 220) * TMath::Sin(electron->GetPhiDC()) - (event->GetPreciseY() / 220) * TMath::Cos(electron->GetPhiDC());
     
            electron->SetAlphaPrime(electron->GetAlpha() - alpha_offset);
            // set Phi0 to right value
            electron->SetPhi0Prime(electron->GetPhi0Prime() - 2.0195 * alpha_offset);

            electron->SetPtPrime(electron->GetPtPrime() * TMath::Abs(electron->GetAlpha() / electron->GetAlphaPrime()) );

            if (electron->GetAlpha() * electron->GetAlphaPrime() < 0)
                electron->SetQPrime(-electron->GetChargePrime());
            else
                electron->SetQPrime(electron->GetChargePrime());

            if(false && electron->GetMcId()==2 && electron->GetChargePrime()==-1 && electron->GetPtPrime()>0.4) std::cout << "Positron: " << electron->GetPtPrime() << " " << electron->GetChargePrime() << std::endl;
            if(false && electron->GetMcId()==3 && electron->GetChargePrime()==+1 && electron->GetPtPrime()>0.4) std::cout << "Electron: " << electron->GetPtPrime() << " " << electron->GetChargePrime() << std::endl;
            
        }
    }

    int MyEventContainer::CircleIntersection(float cx0, float cy0, float r0, float cx1, float cy1, float r1, std::pair<float, float> &p1, std::pair<float, float> &p2)
    {
        float dx = cx1 - cx0;
        float dy = cy1 - cy0;
        float d = sqrt(dx * dx + dy * dy);

        if (d > r0 + r1 || d < fabs(r0 - r1) || d == 0)
            return 0;

        float a = (r0 * r0 - r1 * r1 + d * d) / (2 * d);
        float h = sqrt(r0 * r0 - a * a);

        float xm = cx0 + a * dx / d;
        float ym = cy0 + a * dy / d;

        float rx = -dy * (h / d);
        float ry = dx * (h / d);

        p1.first = xm + rx;
        p1.second = ym + ry;
        p2.first = xm - rx;
        p2.second = ym - ry;

        return 1;
    }

    float MyEventContainer::compute_weighted_median(std::vector<std::pair<float, float> > &value_weight_pairs)
    {
        // Sort by value
        std::sort(value_weight_pairs.begin(), value_weight_pairs.end());

        float total_weight = 0.0;
        for (size_t i = 0; i < value_weight_pairs.size(); ++i)
            total_weight += value_weight_pairs[i].second;

        float cumulative = 0.0;
        for (size_t i = 0; i < value_weight_pairs.size(); ++i)
        {
            cumulative += value_weight_pairs[i].second;
            if (cumulative >= 0.5 * total_weight)
                return value_weight_pairs[i].first;
        }

        return value_weight_pairs.empty() ? 0.0 : value_weight_pairs.back().first;
    }

    void MyEventContainer::VertexReFinder(int fill_hist, int verbosity)
    {

        //const float m_dphi[2] = {0.005 / sqrt(12.), 0.008 / sqrt(12.)};
        /// 0.005  : pixel size of pixel detector in phi-direction
        /// 0.008  : pixel size of stripixel detector in phi-direction
        //const float m_dthe[2] = {0.0425 / sqrt(12.), 0.1 / sqrt(12.)};
        /// 0.0425 : pixel size of pixel detector in z-direction
        /// 0.1    : pixel size of stripixel detector in z-direction

        const float sdphi =  0.05;// + 0*m_dphi[0];
        const float sdthe =  0.01;// + 0*m_dthe[0];
        const float sddphi = 0.005;
        const float sddthe = 0.01;

        std::vector<std::pair<float, float> > track_vertices; // (x at y_beam, y at x_beam)
        std::vector<float> weights_x;
        std::vector<float> weights_y;    

        float beam_x = event->GetPreciseX();//0.322;
        float beam_y = event->GetPreciseY();//0.038;
        float beam_z = event->GetPreciseZ();

        const int nvtx_hits = event->GetNVTXhit();

        int layer_hits[4]= {0, 0, 0, 0};
    
        for (int ihit = 0; ihit < nvtx_hits; ihit++)
        {
            MyDileptonAnalysis::MyVTXHit *myhit1 = event->GetVTXHitEntry(ihit);
            layer_hits[myhit1->GetLayer()]++;
            for (int jhit = ihit + 1; jhit < nvtx_hits*0; jhit++)
            {
                MyDileptonAnalysis::MyVTXHit *myhit2 = event->GetVTXHitEntry(jhit);
                if ( myhit1->GetLayer()==myhit2->GetLayer() && TMath::Abs(myhit1->GetPhi() - myhit2->GetPhi()) < 0.0005 && TMath::Abs(myhit1->GetTheHit() - myhit2->GetTheHit()) < 0.0005)     
                {
                    myhit2->SetZHit(-200);                    
                    if (verbosity) std::cout << "\033[31mdouble hit at layer" << myhit1->GetLayer() << "\033[0m" << std::endl;
                }
            }
        }

        if (verbosity>10) std::cout << "layer_hits:  " << layer_hits[0] << " " << layer_hits[1] << " " << layer_hits[2] << " " << layer_hits[3] << std::endl;

        std::vector<int> used_hits(nvtx_hits, 0);
        std::vector<std::vector<float> > photons, all_circles; 
        
        for (int ihit = 0; ihit < nvtx_hits; ++ihit)
        {
            MyDileptonAnalysis::MyVTXHit *hit0 = event->GetVTXHitEntry(ihit);
            if (hit0->GetLayer() < 2)
                continue;
            if (hit0->GetZHit() < -100)
                continue;
            
            std::vector<std::vector<float> > circle_params;  // stores (cx, cy, R)
            for (int jhit = 0; jhit < nvtx_hits; ++jhit)
            {
                MyDileptonAnalysis::MyVTXHit *hit1 = event->GetVTXHitEntry(jhit);
                if (hit1->GetLayer() != 1)
                    continue;
                if (hit1->GetZHit() < -100)
                    continue;

                float phi0 = hit0->GetPhiHit(beam_x, beam_y, beam_z);
                float the0 = hit0->GetTheHit(beam_x, beam_y, beam_z);
                float phi1 = hit1->GetPhiHit(beam_x, beam_y, beam_z);
                float the1 = hit1->GetTheHit(beam_x, beam_y, beam_z);

                float dphi = phi1 - phi0;
                float dtheta = the1 - the0;

                if (fabs(dphi) > sdphi || fabs(dtheta) > sdthe )
                    continue;
                
                float x0 = hit0->GetXHit();
                float y0 = hit0->GetYHit();
                float x1 = hit1->GetXHit();
                float y1 = hit1->GetYHit();
                
                float x12 = x1 - x0;
                float y12 = y1 - y0;
                float x2b = beam_x - x0;
                float y2b = beam_y - y0;
                float det = x12 * y2b - y12 * x2b;
                if (fabs(det) < 1e-9) continue;

                float A = beam_x * beam_x + beam_y * beam_y;
                float B = x0 * x0 + y0 * y0;
                float C = x1 * x1 + y1 * y1;

                float cx0 = (A * (y0 - y1) + B * (y1 - beam_y) + C * (beam_y - y0)) / (2 * det);
                float cy0 = (A * (x1 - x0) + B * (beam_x - x1) + C * (x0 - beam_x)) / (2 * det);
                float R0 = sqrt((beam_x - cx0) * (beam_x - cx0) + (beam_y - cy0) * (beam_y - cy0));
                float pt0 = R0 * (0.003 * 0.9); 

                for (int khit = 0; khit < nvtx_hits; ++khit)
                {
                    MyDileptonAnalysis::MyVTXHit *hit2 = event->GetVTXHitEntry(khit);
                    if (hit2->GetLayer() != 0 )
                        continue;
                    if (hit1->GetZHit() < -100)
                        continue;

                    float phi2 = hit2->GetPhiHit(beam_x, beam_y, beam_z);
                    float the2 = hit2->GetTheHit(beam_x, beam_y, beam_z);
                    float x2 = hit2->GetXHit();
                    float y2 = hit2->GetYHit();

                    float r0  = sqrt( SQR(x2) + SQR(y2) );
                    float r1  = sqrt( SQR(x1) + SQR(y1) );
                    float r01 = r0 - r1;
                    float r12 = r1 - sqrt( SQR(x0) + SQR(y0) );
                    //float z01 = hit1->GetZHit()-hit0->GetZHit();
                    //float z02 = hit2->GetZHit()-hit0->GetZHit();

                    float dphi1 = phi2 - ( phi1 + dphi * r01 / r12 );
                    float dthe1 = the2 - the1;// - ( the0 + dtheta * z02 / z01 );

                    if (fabs(dthe1) < sddthe && fill_hist)
                        vtx_dphi_dphi_hist->Fill(dphi1, dphi, pt0);
                    if (fabs(dphi1) < sddphi && fill_hist)
                        vtx_dthe_dthe_hist->Fill(dthe1, dtheta, pt0);
                        
                    if (fabs(dphi1) > sddphi || fabs(dthe1) > sddthe )
                        continue;

                    if(false)
                    {
                        float r_hit = sqrt(x2 * x2 + y2 * y2);
                        std::pair<float, float> inter1, inter2;
                        if (!CircleIntersection(cx0, cy0, R0, 0.0, 0.0, r_hit, inter1, inter2)) continue;
                        
                        // Pick intersection point closest to hit2
                        float dist1 = (x2 - inter1.first)*(x2 - inter1.first) + (y2 - inter1.second)*(y2 - inter1.second);
                        float dist2 = (x2 - inter2.first)*(x2 - inter2.first) + (y2 - inter2.second)*(y2 - inter2.second);
                        float px = (dist1 < dist2) ? inter1.first : inter2.first;
                        float py = (dist1 < dist2) ? inter1.second : inter2.second;
                        
                        // New projected angle on the track
                        float phi_proj = atan2(py - beam_y, px - beam_x);

                        float dphir = phi2 - phi_proj;       

                        if (fabs(dphir) > sddphi) continue;
                    }

                    float x12 = x1 - x0;
                    float y12 = y1 - y0;
                    float x23 = x2 - x1;
                    float y23 = y2 - y1;
                    float det = x12 * y23 - y12 * x23;
                    if (fabs(det) < 1e-6)
                        continue; // collinear points

                    float A = x0 * x0 + y0 * y0;
                    float B = x1 * x1 + y1 * y1;
                    float C = x2 * x2 + y2 * y2;

                    float cx = (A * (y1 - y2) + B * (y2 - y0) + C * (y0 - y1)) / (2 * det);
                    float cy = (A * (x2 - x1) + B * (x0 - x2) + C * (x1 - x0)) / (2 * det);
                    float R = sqrt((x0 - cx) * (x0 - cx) + (y0 - cy) * (y0 - cy));

                    std::vector<float> best_circle(4);
                    best_circle[0] = cx;
                    best_circle[1] = cy;
                    best_circle[2] = ((phi2 - phi0) > 0 ? 1 : -1) * R;
                    best_circle[3] = TMath::Abs(dphi1);

                    if (false)
                    {
                        float phi_proj0 = atan2( y0 - cy, x0 - cx);
                        float phi0_c = atan2(y0 - beam_y, x0 - beam_x);
                        while ( phi_proj0 - phi0_c > 3*TMath::Pi()/8 )
                            phi_proj0 -=  TMath::Pi()/2;
                        
                        while ( phi_proj0 - phi0_c < -3*TMath::Pi()/8 )
                            phi_proj0 +=  TMath::Pi()/2;
                        //std::cout << "phi_proj0 = " << phi_proj0 << " phi0_c = " << phi0_c << std::endl;
                        best_circle[3] = phi_proj0;
                        best_circle[4] = the0;
                    }

                    if (fill_hist)
                        phi_the_pt_hist->Fill(hit2->GetPhi()-dphi1*r0/r01, hit2->GetTheHit(), best_circle[2] * (0.003 * 0.9));
    
                    circle_params.push_back(best_circle);
                }
            }
            
            if(fill_hist)
                hits_vtx_ntracks->Fill(circle_params.size(), event->GetCentrality());
            
            int charge_flip = 0;
            
            if (false)
            {
                for (size_t i = 0; i < circle_params.size(); ++i)
                {   
                    int charge = (circle_params[i][2] > 0) ? 1 : -1;
                    for (size_t j = i+1; j < circle_params.size(); ++j)
                    {
                        if (charge * circle_params[j][2] < 0)
                        {
                            const double pt1 = TMath::Abs(circle_params[i][2] * (0.003 * 0.9));
                            const double pt2 = TMath::Abs(circle_params[j][2] * (0.003 * 0.9));
                            if(pt1 > 1.5 || pt2 > 1.5) continue;

                            const double px1 = pt1 * cos(circle_params[i][3]);
                            const double py1 = pt1 * sin(circle_params[i][3]);
                            const double pz1 = pt1 / tan(circle_params[i][4]);
                            const double px2 = pt2 * cos(circle_params[j][3]);
                            const double py2 = pt2 * sin(circle_params[j][3]);
                            const double pz2 = pt2 / tan(circle_params[j][4]);

                            const double pt_pair = sqrt(SQR(px1 + px2) + SQR(py1 + py2));
                            const double p_pair = sqrt(SQR(pt_pair) + SQR(pz1 + pz2));
                            const double es = sqrt(SQR(px1) + SQR(py1) + SQR(pz1) + SQR(0.000511)) + sqrt(SQR(px2) + SQR(py2) + SQR(pz2) + SQR(0.000511));
                            double m2 = SQR(es) - SQR(p_pair);
                            if (m2 < 0) m2 = 0;  // Avoid NaN due to numerical precision
                            const double invm = sqrt(m2);

                            if(fill_hist)
                                conv_photon_mass_hist->Fill(invm, pt_pair, event->GetCentrality());

                            if (verbosity>100) std::cout << "\033[31mcharge flip to Q = " << charge << "\033[0m with" << invm << " " << pt_pair << std::endl;

                            std::vector<float> photon(3);
                            photon[0] = px1 + px2;
                            photon[1] = py1 + py2;
                            photon[2] = pz1 + pz2;
                            photons.push_back(photon);
                        }
                    }
                }
            }

            if (circle_params.size() > 1)
            {
                //std::cout << "\033[31mcircle_params.size() = " << circle_params.size() << "\033[0m" << std::endl;
                float min_R = 999999;
                int min_index = -1;
                
                int charge = 0;
                for (size_t i = 0; i < circle_params.size(); ++i)
                {   
                    if ( charge*circle_params[i][2] < 0 )
                    {
                        //std::cout << "\033[31mcharge flip to Q = " << charge << "\033[0m" << std::endl;
                        charge_flip = 1;
                    }
                    charge = (circle_params[i][2] > 0) ? 1 : -1;
                    if (TMath::Abs(circle_params[i][3]) < min_R)
                    {
                        min_R = TMath::Abs(circle_params[i][3]);
                        min_index = i;
                    }
                }

                std::vector<std::vector<float> > keep_one;
                if (min_index >= 0)
                {
                    keep_one.push_back(circle_params[min_index]);
                }
                circle_params = keep_one;
                if(false) circle_params.pop_back();
            }

            if (charge_flip) continue;
            if (circle_params.size()) all_circles.push_back(circle_params[0]);

            for (size_t icircle = 0; icircle < circle_params.size(); ++icircle)
            {
                float cx = circle_params[icircle][0];
                float cy = circle_params[icircle][1];
                float R = TMath::Abs( circle_params[icircle][2] );

                float pt = R * (0.003 * 0.9);
                // Evaluate at y = beam_y to find x// Always project to beam_y and beam_x without rejecting

                //float dx = sqrt((cx - beam_x) * (cx - beam_x) + (cy - beam_y) * (cy - beam_y));
                float angle_beam = atan2(beam_y - cy, beam_x - cx);

                float x_proj = cx + R * cos(angle_beam);
                float y_proj = cy + R * sin(angle_beam);

                float weight = 1.0 / (0.115 / pt + 0.555); // / (1e-6 + 0.1 * pt / (1 + pt)); // rescattering-based weight

                float w_x = fabs(cos(angle_beam)) * weight;
                float w_y = fabs(sin(angle_beam)) * weight;

                track_vertices.push_back(std::make_pair(x_proj, y_proj));
                weights_x.push_back(w_x);
                weights_y.push_back(w_y);

                const float dca = sqrt(SQR(cx - beam_x) + SQR(cy - beam_y)) - R;
                const float dca_beam = sqrt(SQR(cx - 0.322) + SQR(cy - 0.038)) - R;

                if (fill_hist)
                    hist_dca_x->Fill(dca, pt, event->GetCentrality(), 1);
                if (fill_hist)
                    hist_dca_y->Fill(dca_beam, pt, event->GetCentrality(), 1);
            }
        }

        if (fill_hist)
            hits_vtx_ntracks->Fill(track_vertices.size(), event->GetCentrality());

        if (false && photons.size() > 1)
        {
            for (size_t i = 0; i < photons.size(); ++i)
            {
                for (size_t j = i + 1; i < photons.size(); i++)
                {
                    const double px1 = photons[i][0];
                    const double py1 = photons[i][1];
                    const double pz1 = photons[i][2];
                    const double px2 = photons[j][0];
                    const double py2 = photons[j][1];
                    const double pz2 = photons[j][2];

                    const double pt_pair = sqrt(SQR(px1 + px2) + SQR(py1 + py2));
                    const double p_pair = sqrt(SQR(pt_pair) + SQR(pz1 + pz2));
                    const double es = sqrt(SQR(px1) + SQR(py1) + SQR(pz1)) + sqrt(SQR(px2) + SQR(py2) + SQR(pz2));
                    double m2 = SQR(es) - SQR(p_pair);
                    if (m2 < 0)
                        m2 = 0; // Avoid NaN due to numerical precision
                    const double invm = sqrt(m2);

                    if (fill_hist)
                        pi0_mass_hist->Fill(invm, pt_pair, event->GetCentrality());

                    if (invm > 0.09 && invm < 0.2 && verbosity>100)
                        std::cout << "\033[32mpi0 with" << invm << " " << pt_pair << "\033[0m" << std::endl;
                }
            }
        }

        if (!track_vertices.empty())
        {

            // Prepare value-weight pairs
            std::vector<std::pair<float, float> > x_val_wt;
            std::vector<std::pair<float, float> > y_val_wt;

            for (size_t i = 0; i < track_vertices.size(); ++i)
            {
                x_val_wt.push_back(std::make_pair(track_vertices[i].first, weights_x[i]));
                y_val_wt.push_back(std::make_pair(track_vertices[i].second, weights_y[i]));
            }

            // Compute weighted medians
            float median_x = compute_weighted_median(x_val_wt);
            float median_y = compute_weighted_median(y_val_wt);

            // Optionally print and fill histograms
            if (verbosity)
            {
                std::cout << "\033[34mWeighted median VTX: " << median_x << " " << median_y << "\033[0m"<< std::endl;
            }

            float sumx = 0.0, sumy = 0.0;
            float sumwx = 0.0, sumwy = 0.0;
    
            for (size_t i = 0; i < track_vertices.size(); ++i)
            {
                sumx += track_vertices[i].first * weights_x[i];
                sumwx += weights_x[i];
                sumy += track_vertices[i].second * weights_y[i];
                sumwy += weights_y[i];
            }
    
            float vx = (sumwx > 0) ? sumx / sumwx : beam_x;
            float vy = (sumwy > 0) ? sumy / sumwy : beam_y;

            if (verbosity) std::cout << "prevous vertex: " << event->GetPreciseX() << " "  << event->GetPreciseY() 
            << ";  \033[32mnew VTX: " << vx << " " << vy << "\033[0m" << " using " <<track_vertices.size() << " tracks at "<< (int) event->GetCentrality() << std::endl;

            if (fill_hist)
                hist_vtx_x->Fill(median_x - event->GetPreciseX(), median_x, event->GetCentrality());
            if (fill_hist)
                hist_vtx_y->Fill(median_y - event->GetPreciseY(), median_y, event->GetCentrality());

                //event->SetPreciseX(median_x);
                //event->SetPreciseY(median_y);
        }
        else
        {
            if (verbosity) std::cout << "\033[31mNo VTX found!\033[0m"<< " at "<< (int) event->GetCentrality() << " using " <<event->GetPreciseX()<<" "<<event->GetPreciseY()
            << " hits in layers: " << layer_hits[0] << " " << layer_hits[1] << " " << layer_hits[2] << " " << layer_hits[3]  << std::endl;
        }
        if (track_vertices.size() > 10)
        {

            float vx1 = 0, vy1 = 0, sumwx1 = 0, sumwy1 = 0;
            float vx2 = 0, vy2 = 0, sumwx2 = 0, sumwy2 = 0;

            for (size_t i = 0; i < track_vertices.size(); ++i)
            {
                if ( i % 2 == 0 )
                {
                    vx1 += track_vertices[i].first * weights_x[i];
                    sumwx1 += weights_x[i];
                    vy1 += track_vertices[i].second * weights_y[i];
                    sumwy1 += weights_y[i];
                }
                else
                {
                    vx2 += track_vertices[i].first * weights_x[i];
                    sumwx2 += weights_x[i];
                    vy2 += track_vertices[i].second * weights_y[i];
                    sumwy2 += weights_y[i];
                }
            }

            if (sumwx1 > 0 && sumwx2 > 0 && sumwy1 > 0 && sumwy2 > 0)
            {
                float diffx = vx1 / sumwx1 - vx2 / sumwx2;
                float diffy = vy1 / sumwy1 - vy2 / sumwy2;

                if (false)
                    hist_vtx_delta_x->Fill(diffx,event->GetPreciseX(), event->GetCentrality());
                if (false)
                    hist_vtx_delta_y->Fill(diffy,event->GetPreciseY(), event->GetCentrality());
            }
        }

        if (false)
        {
            std::vector<std::pair<float, float> > intersections; // (x at y_beam, y at x_beam)
            std::vector<float> weights_inter_x;
            std::vector<float> weights_inter_y;
            for (size_t icircle = 0; icircle < all_circles.size(); ++icircle)
            {
                float cx1 = all_circles[icircle][0];
                float cy1 = all_circles[icircle][1];
                float R1 = TMath::Abs(all_circles[icircle][2]);
                float dca1 = TMath::Abs(sqrt((cx1 - beam_x) * (cx1 - beam_x) + (cy1 - beam_y) * (cy1 - beam_y)) - R1);
                if (dca1 > 0.015)
                    continue; // reject this circle if DCA > 150 μm
                float pt1 = R1 * (0.003 * 0.9);

                for (size_t jcircle = icircle + 1; jcircle < all_circles.size(); ++jcircle)
                {
                    float cx2 = all_circles[jcircle][0];
                    float cy2 = all_circles[jcircle][1];
                    float R2 = TMath::Abs(all_circles[jcircle][2]);
                    float dca2 = TMath::Abs(sqrt((cx2 - beam_x) * (cx2 - beam_x) + (cy2 - beam_y) * (cy2 - beam_y)) - R2);
                    if (dca2 > 0.015)
                        continue; // reject this circle if DCA > 150 μm
                    float pt2 = R2 * (0.003 * 0.9);

                    std::pair<float, float> inter1, inter2;
                    if (!CircleIntersection(cx1, cy1, R1, cx2, cy2, R2, inter1, inter2))
                        continue;

                    float d1 = SQR(inter1.first - beam_x) + SQR(inter1.second - beam_y);
                    float d2 = SQR(inter2.first - beam_x) + SQR(inter2.second - beam_y);
                    float x_proj = (d1 < d2) ? inter1.first : inter2.first;
                    float y_proj = (d1 < d2) ? inter1.second : inter2.second;

                    // Tangents at intersection point
                    float tx1 = -(y_proj - cy1);
                    float ty1 = (x_proj - cx1);
                    float tx2 = -(y_proj - cy2);
                    float ty2 = (x_proj - cx2);

                    float norm1 = sqrt(tx1 * tx1 + ty1 * ty1);
                    float norm2 = sqrt(tx2 * tx2 + ty2 * ty2);
                    if (norm1 < 1e-6 || norm2 < 1e-6)
                        continue;

                    // Normalize and get cosine of angle
                    float dot = (tx1 * tx2 + ty1 * ty2) / (norm1 * norm2);
                    dot = std::max(-1.0f, std::min(1.0f, dot)); // Clamp for safety
                    // float delta_alpha = acos(dot);
                    float geom_weight = sqrt(1. - SQR(dot)); // fabs(sin(delta_alpha)); // ~0 if parallel, ~1 if orthogonal

                    // Beam direction angle
                    float beam_angle1 = atan2(cy1 - beam_y, cx1 - beam_x);
                    float beam_angle2 = atan2(cy2 - beam_y, cx2 - beam_x);
                    float w_x = sqrt(fabs(cos(beam_angle1) * cos(beam_angle2)) * SQR(geom_weight) / (0.115 / pt1 + 0.555) / (0.115 / pt2 + 0.555));
                    float w_y = sqrt(fabs(sin(beam_angle1) * sin(beam_angle2)) * SQR(geom_weight) / (0.115 / pt1 + 0.555) / (0.115 / pt2 + 0.555));

                    intersections.push_back(std::make_pair(x_proj, y_proj));
                    weights_inter_x.push_back(w_x);
                    weights_inter_y.push_back(w_y);

                    // now you can store or use (x_proj, y_proj) and weights (w_x, w_y)
                }
            }

            if (!intersections.empty())
            {
                float sumx = 0.0, sumy = 0.0;
                float sumwx = 0.0, sumwy = 0.0;

                for (size_t i = 0; i < intersections.size(); ++i)
                {
                    sumx += intersections[i].first * weights_inter_x[i];
                    sumwx += weights_inter_x[i];
                    sumy += intersections[i].second * weights_inter_y[i];
                    sumwy += weights_inter_y[i];
                }

                float vx = (sumwx > 0) ? sumx / sumwx : beam_x;
                float vy = (sumwy > 0) ? sumy / sumwy : beam_y;

                if (verbosity)
                    std::cout << "\033[36mnovel VTX: " << vx << " " << vy << "\033[0m" << " using " << track_vertices.size() << " tracks at " << (int)event->GetCentrality() << std::endl;

                // if (fill_hist)
                //     hist_vtx_x->Fill(vx - event->GetPreciseX(), vx, event->GetCentrality());
                // if (fill_hist)
                //     hist_vtx_y->Fill(vy - event->GetPreciseY(), vy, event->GetCentrality());

                // event->SetPreciseX(vx);
                // event->SetPreciseY(vy);
            }
        }
    
    }

    void MyEventContainer::CorrectVTXOffset(int keff)
    {
        for (int ihit = 0; ihit < event->GetNVTXhit(); ++ihit)
        {
            MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
            const int layer = myhit->GetLayer();
            if (layer < 1) continue; // skip striplex
            for (int ilayer = 0; ilayer < layer; ilayer++)
            {
                const double x = myhit->GetXHit();
                const double y = myhit->GetYHit();
                const double r = sqrt(SQR(x) + SQR(y)); 
                int z_bin = (myhit->GetZHit()+12)/3; // 24 cm range, 3 cm per bin
                if (z_bin < 0) z_bin = 0;
                if (z_bin > 7) z_bin = 7; // 24 cm range
                const int arm = x>0 ? 1 : 0; // east is 0; west is 1
                const int ii = arm + ilayer*2;
                const double phi = myhit->GetPhiHit(0,0,0);
                int phi_bin = 0;
                for (int iphi = 0; iphi < 6; ++iphi)
                {
                    if (phi > phi_vtx_bins[arm][iphi][0] && phi <= phi_vtx_bins[arm][iphi][1])
                    {
                        phi_bin = iphi;
                        break;
                    }
                }
                double delta_x = 0; 
                double delta_y = 0; 
                if(keff>0)
                {
                    delta_x = VTX_seg_alignment[ii][z_bin][phi_bin][0] / 10000 * r + VTX_seg_alignment[ii][z_bin][phi_bin][2]/10000 * r * TMath::Sin(phi);
                    delta_y = VTX_seg_alignment[ii][z_bin][phi_bin][1] / 10000 * r + VTX_seg_alignment[ii][z_bin][phi_bin][2]/10000 * r * TMath::Cos(phi);
                }else
                {
                    delta_x = VTX_seg_alignment1[ii][z_bin][phi_bin][0] / 10000 * r + VTX_seg_alignment1[ii][z_bin][phi_bin][2]/10000 * r * TMath::Sin(phi);
                    delta_y = VTX_seg_alignment1[ii][z_bin][phi_bin][1] / 10000 * r + VTX_seg_alignment1[ii][z_bin][phi_bin][2]/10000 * r * TMath::Cos(phi);
                }

                //const double dphi = VTX_seg_alignment[ii][z_bin][phi_bin][0] / 10000 * TMath::Sin(phi) + VTX_seg_alignment[ii][z_bin][phi_bin][1] / 10000 * TMath::Cos(phi) + VTX_seg_alignment[ii][z_bin][phi_bin][2] / 10000;
                //const double delta_x = dphi * r * TMath::Sin(phi);
                //const double delta_y = dphi * r * TMath::Cos(phi); 
                myhit->SetXHit(x + keff*delta_x);
                myhit->SetYHit(y - keff*delta_y);
            }
        }
        for (int ihit = 0; ihit < event->GetNVTXhit(); ++ihit)
        {
            MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
            const float z = myhit->GetZHit();
            const float x = myhit->GetXHit();
            const float y = myhit->GetYHit();
            const float phi = myhit->GetPhiHit(0, 0, 0);
            int z_bin = (z + 12) / 3; // 24 cm range, 3 cm per bin
            if (z_bin < 0) z_bin = 0;
            if (z_bin > 7) z_bin = 7; // 24 cm range
            const int phi_bin = (phi + 1.5) / 0.3;
            const float offset = offsets_dca[z_bin][phi_bin];
            const float offsets_x = offset*cos(phi+TMath::Pi()/2);
            const float offsets_y = offset*sin(phi+TMath::Pi()/2);
            myhit->SetXHit(x - keff*offsets_x);
            myhit->SetYHit(y - keff*offsets_y);
        }
        if (false)
        {
            for (int ihit = 0; ihit < event->GetNVTXhit(); ++ihit)
            {
                MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
            
                const float y = myhit->GetYHit();
            
                const float z = myhit->GetZHit();
                int z_bin = (z + 12) / 3; // 24 cm range, 3 cm per bin
                if (z_bin < 0) z_bin = 0;
                if (z_bin > 7) z_bin = 7; // 24 cm range
                const float phi = myhit->GetPhiHit(0, 0, 0);
                const int phi_bin = (phi < 1.5) ? 0 : 1; 
                const float offsets_y = offsets_vtx_arms_y[z_bin][phi_bin];
                myhit->SetYHit(y - keff*offsets_y);
            }
            for (int ihit = 0; ihit < event->GetNVTXhit()*0; ++ihit)
            {
                MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
                int layer = myhit->GetLayer();
                if (layer == 0 ) continue; // skip most inner
                for (int ilayer = 0; ilayer < layer; ilayer++)
                {
                    const double z = myhit->GetZHit();
                    const double x = myhit->GetXHit();
                    const double y = myhit->GetYHit();
                    const double r = sqrt(SQR(x) + SQR(y));
                    const int arm = x>0 ? 1 : 0; // east is 0; west is 1
                    const int z_bin = z>0 ? 1 : 0; // south is 0; north is 1
                    const double delta_x =  (vtx_pixel_offsets[arm][ilayer][z_bin][0][0] + vtx_pixel_offsets[arm][ilayer][z_bin][0][1] * z) / 10000 * r;
                    const double delta_y = -(vtx_pixel_offsets[arm][ilayer][z_bin][1][0] + vtx_pixel_offsets[arm][ilayer][z_bin][1][1] * z) / 10000 * r;
                    myhit->SetXHit(x - keff*delta_x);
                    myhit->SetYHit(y - keff*delta_y);
                }
            }
            for (int ihit = 0; ihit < event->GetNVTXhit()*0; ++ihit)
            {
                MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
                const int layer = myhit->GetLayer();
                if (layer > 2) continue; // skip striplex
                for (int ilayer = layer; ilayer < 3; ilayer++)
                {
                    const double x = myhit->GetXHit();
                    const double y = myhit->GetYHit();
                    const double r = sqrt(SQR(x) + SQR(y)); 
                    int z_bin = (myhit->GetZHit()+12)/3; // 24 cm range, 3 cm per bin
                    if (z_bin < 0) z_bin = 0;
                    if (z_bin > 7) z_bin = 7; // 24 cm range
                    const int arm = x>0 ? 1 : 0; // east is 0; west is 1
                    const int ii = arm + ilayer*2;
                    const double phi = myhit->GetPhiHit(0,0,0);
                    int phi_bin = 0;
                    for (int iphi = 0; iphi < 6; ++iphi)
                    {
                        if (phi > phi_vtx_bins[arm][iphi][0] && phi <= phi_vtx_bins[arm][iphi][1])
                        {
                            phi_bin = iphi;
                            break;
                        }
                    }
                    const double delta_x =  VTX_seg_alignment[ii][z_bin][phi_bin][0] / 10000 * r + VTX_seg_alignment[ii][z_bin][phi_bin][2]/10000 * r * TMath::Sin(phi);
                    const double delta_y =  VTX_seg_alignment[ii][z_bin][phi_bin][1] / 10000 * r + VTX_seg_alignment[ii][z_bin][phi_bin][2]/10000 * r * TMath::Cos(phi);
                    //const double dphi = VTX_seg_alignment[ii][z_bin][phi_bin][0] / 10000 * TMath::Sin(phi) + VTX_seg_alignment[ii][z_bin][phi_bin][1] / 10000 * TMath::Cos(phi) + VTX_seg_alignment[ii][z_bin][phi_bin][2] / 10000;
                    //const double delta_x = dphi * r * TMath::Sin(phi);
                    //const double delta_y = dphi * r * TMath::Cos(phi); 
                    myhit->SetXHit(x - keff*delta_x);
                    myhit->SetYHit(y + keff*delta_y);
                }
            }

            for (int ihit = 0; ihit < event->GetNVTXhit()*0; ++ihit)
            {
                MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
                const float z = myhit->GetZHit();
                const float x = myhit->GetXHit();
                const float y = myhit->GetYHit();
                const float phi = myhit->GetPhiHit(0, 0, 0);
                int z_bin = (z + 12) / 3; // 24 cm range, 3 cm per bin
                if (z_bin < 0) z_bin = 0;
                if (z_bin > 7) z_bin = 7; // 24 cm range
                const int phi_bin = (phi + 1.5) / 0.3;
                const float offsets_x = offsets_dca_x[z_bin][phi_bin];
                const float offsets_y = offsets_dca_y[z_bin][phi_bin];
                myhit->SetXHit(x - keff*offsets_x);
                myhit->SetYHit(y - keff*offsets_y);
            }
            for (int ihit = 0; ihit < event->GetNVTXhit()*0; ++ihit)
            {
                MyDileptonAnalysis::MyVTXHit *myhit = event->GetVTXHitEntry(ihit);
                const float x = myhit->GetXHit();
                const float y = myhit->GetYHit();
                const float random_gaus_offset_x = 0.0010 * (gRandom->Gaus(0, 1));
                const float random_gaus_offset_y = 0.0010 * (gRandom->Gaus(0, 1));
                myhit->SetXHit(x + random_gaus_offset_x);
                myhit->SetYHit(y + random_gaus_offset_y);
            }
        }
    }

    // Function to scan vertex Z and find the best estimate using circular track matching
    void MyEventContainer::VertexXYScan(const float run_beam_x = 0.328, const float run_beam_y = 0.038, int fill_hist = 0, int verbosity = 0)
    {
        float sdphi = 0.05;      // +/- rad
        float sdthe = 0.01;      // +/- rad
        float sddphi = 0.005;    // +/- rad
        float sddthe = 0.01;     // +/- rad
        float x_range = 0.05;    // +/- cm
        float y_range = 0.05;    // +/- cm
        float step_size = 0.005; // cm = 25 microns
        int min_track_count = 3;
        float beam_x = run_beam_x;
        float beam_y = run_beam_y;
        float beam_z = event->GetPreciseZ();
        float best_x = beam_x;
        float best_y = beam_y;
        int max_tracks = 0;
        const float fraction = 0.75; // fraction of max-min tracks to use for vertex calculation

        const int nvtx_hits = event->GetNVTXhit();

        std::vector<std::vector<int> > hits_chains; // stores (cx, cy, R)
        for (int ihit = 0; ihit < nvtx_hits; ++ihit)
        {
            MyDileptonAnalysis::MyVTXHit *hit0 = event->GetVTXHitEntry(ihit);
            if (hit0->GetLayer() < 2)
                continue;
            for (int jhit = 0; jhit < nvtx_hits; ++jhit)
            {
                MyDileptonAnalysis::MyVTXHit *hit1 = event->GetVTXHitEntry(jhit);
                if (hit1->GetLayer() != 1)
                    continue;

                float phi0 = hit0->GetPhiHit(beam_x, beam_y, beam_z);
                float the0 = hit0->GetTheHit(beam_x, beam_y, beam_z);
                float phi1 = hit1->GetPhiHit(beam_x, beam_y, beam_z);
                float the1 = hit1->GetTheHit(beam_x, beam_y, beam_z);

                float dphi = phi1 - phi0;
                float dtheta = the1 - the0;

                if (fabs(dphi) > sdphi + 0.01 || fabs(dtheta) > sdthe)
                    continue;

                float x0 = hit0->GetXHit();
                float y0 = hit0->GetYHit();
                float x1 = hit1->GetXHit();
                float y1 = hit1->GetYHit();

                for (int khit = 0; khit < nvtx_hits; ++khit)
                {
                    MyDileptonAnalysis::MyVTXHit *hit2 = event->GetVTXHitEntry(khit);
                    if (hit2->GetLayer() != 0)
                        continue;

                    float phi2 = hit2->GetPhiHit(beam_x, beam_y, beam_z);
                    float the2 = hit2->GetTheHit(beam_x, beam_y, beam_z);
                    float x2 = hit2->GetXHit();
                    float y2 = hit2->GetYHit();

                    float r0 = sqrt(SQR(x2) + SQR(y2));
                    float r1 = sqrt(SQR(x1) + SQR(y1));
                    float r01 = r0 - r1;
                    float r12 = r1 - sqrt(SQR(x0) + SQR(y0));

                    float dphi1 = phi2 - (phi1 + dphi * r01 / r12);
                    float dthe1 = the2 - the1; // - ( the0 + dtheta * z02 / z01 );

                    if (fabs(dphi1) > sddphi + 0.015 || fabs(dthe1) > sddthe)
                        continue;

                    std::vector<int> iter_hit(3);
                    iter_hit[0] = ihit;
                    iter_hit[1] = jhit;
                    iter_hit[2] = khit;
                    hits_chains.push_back(iter_hit);
                }
            }
        }
        
        std::map<int, std::vector< std::pair<int, int> > > khit_to_chain;

        for (size_t i = 0; i < hits_chains.size(); ++i)
        {
            const std::vector<int> &chain = hits_chains[i];
            if (chain.size() != 3)
                continue; // safety check

            int ihit = chain[0];
            int jhit = chain[1];
            int khit = chain[2];

            khit_to_chain[khit].push_back(std::make_pair(ihit, jhit));
        }

        if (verbosity)
            std::cout << "\n\033[32mFound " << hits_chains.size() << " " << khit_to_chain.size() << std::endl;

        std::vector<std::pair<float, float> > track_vertices, track_vertices_west, track_vertices_east, track_vertices_notused, track_vertices_pions; // (x at y_beam, y at x_beam)
        std::vector<int> n_tracks_vec, n_tracks_vec_west, n_tracks_vec_east, n_tracks_vec_notused, n_tracks_vec_pions;
        for (float xvtx = beam_x - x_range; xvtx <= beam_x + x_range; xvtx += step_size)
        {
            for (float yvtx = beam_y - y_range; yvtx <= beam_y + y_range; yvtx += step_size)
            {
                int n_tracks = 0, n_wtracks = 0;
                int n_tracks_east = 0, n_tracks_west = 0, n_tracks_notused = 0, n_tracks_pions = 0;
                std::vector<float> weight_of_tracks;
                std::vector<int> used_outer_hits;
                for (std::map<int, std::vector<std::pair<int, int> > >::iterator it = khit_to_chain.begin();
                     it != khit_to_chain.end(); ++it)
                {
                    const unsigned int khit = it->first;
                    const std::vector<std::pair<int, int> > &layer12_hits = it->second;
                    std::vector<std::vector<float> > circle_params; // stores (cx, cy, R)
                    for (unsigned int ichain = 0; ichain < layer12_hits.size(); ++ichain)
                    {
                        const unsigned int ihit = layer12_hits[ichain].first;
                        const unsigned int jhit = layer12_hits[ichain].second;
                        if (used_outer_hits.size() > 0 && ( std::find(used_outer_hits.begin(), used_outer_hits.end(), ihit) != used_outer_hits.end() ||
                            std::find(used_outer_hits.begin(), used_outer_hits.end(), jhit) != used_outer_hits.end() ))
                            continue; // skip if this khit is already used in another chain

                        MyDileptonAnalysis::MyVTXHit *hit0 = event->GetVTXHitEntry(ihit);
                        if (hit0->GetLayer() < 2)
                            continue;

                        MyDileptonAnalysis::MyVTXHit *hit1 = event->GetVTXHitEntry(jhit);
                        if (hit1->GetLayer() != 1)
                            continue;

                        float phi0 = hit0->GetPhiHit(xvtx, yvtx, beam_z);
                        float the0 = hit0->GetTheHit(xvtx, yvtx, beam_z);
                        float phi1 = hit1->GetPhiHit(xvtx, yvtx, beam_z);
                        float the1 = hit1->GetTheHit(xvtx, yvtx, beam_z);

                        float dphi = phi1 - phi0;
                        float dtheta = the1 - the0;

                        if (fabs(dphi) > sdphi || fabs(dtheta) > sdthe)
                            continue;

                        float x0 = hit0->GetXHit();
                        float y0 = hit0->GetYHit();
                        float x1 = hit1->GetXHit();
                        float y1 = hit1->GetYHit();

                        MyDileptonAnalysis::MyVTXHit *hit2 = event->GetVTXHitEntry(khit);
                        if (hit2->GetLayer() != 0)
                            continue;

                        float phi2 = hit2->GetPhiHit(xvtx, yvtx, beam_z);
                        float the2 = hit2->GetTheHit(xvtx, yvtx, beam_z);
                        float x2 = hit2->GetXHit();
                        float y2 = hit2->GetYHit();

                        float r0 = sqrt(SQR(x2) + SQR(y2));
                        float r1 = sqrt(SQR(x1) + SQR(y1));
                        float r01 = r0 - r1;
                        float r12 = r1 - sqrt(SQR(x0) + SQR(y0));

                        float dphi1 = phi2 - (phi1 + dphi * r01 / r12);
                        float dthe1 = the2 - the1; // - ( the0 + dtheta * z02 / z01 );

                        if (fabs(dphi1) > sddphi || fabs(dthe1) > sddthe)
                            continue;

                        float x12 = x1 - x0;
                        float y12 = y1 - y0;
                        float x23 = x2 - x1;
                        float y23 = y2 - y1;
                        float det = x12 * y23 - y12 * x23;
                        if (fabs(det) < 1e-6)
                            continue; // collinear points

                        float A = x0 * x0 + y0 * y0;
                        float B = x1 * x1 + y1 * y1;
                        float C = x2 * x2 + y2 * y2;

                        float cx = (A * (y1 - y2) + B * (y2 - y0) + C * (y0 - y1)) / (2 * det);
                        float cy = (A * (x2 - x1) + B * (x0 - x2) + C * (x1 - x0)) / (2 * det);
                        float R = sqrt((x0 - cx) * (x0 - cx) + (y0 - cy) * (y0 - cy));

                        std::vector<float> best_circle(7);
                        best_circle[0] = cx;
                        best_circle[1] = cy;
                        best_circle[2] = ((phi2 - phi0) > 0 ? 1 : -1) * R;
                        best_circle[3] = (hit2->GetXHit() > 0 ? 1:-1) * TMath::Abs(dphi1);
                        best_circle[4] = (((hit0->GetLadder()>24&&hit0->GetLadder()<48)) || (hit1->GetLadder()>24&&hit1->GetLadder()<48) || (hit2->GetLadder()>24&&hit2->GetLadder()<48)) ? 1 : 0;
                        if(((hit0->GetLadder()>43&&hit0->GetLadder()<48)) || (hit1->GetLadder()>43&&hit1->GetLadder()<48) || (hit2->GetLadder()>43&&hit2->GetLadder()<48)) best_circle[4] += 2; //
                        best_circle[5] = ihit; // store ihit
                        best_circle[6] = jhit; // store jhit

                        circle_params.push_back(best_circle);

                    } // end of loop over chains

                    if (true) //using only one circle and no conversions:
                    {
                        int charge_flip = 0;
                        if (circle_params.size() > 1)
                        {
                            // std::cout << "\033[31mcircle_params.size() = " << circle_params.size() << "\033[0m" << std::endl;
                            float min_R = 999999;
                            int min_index = -1;

                            int charge = 0;
                            for (size_t i = 0; i < circle_params.size(); ++i)
                            {
                                if (charge * circle_params[i][2] < 0)
                                {
                                    // std::cout << "\033[31mcharge flip to Q = " << charge << "\033[0m" << std::endl;
                                    charge_flip = 1;
                                }
                                charge = (circle_params[i][2] > 0) ? 1 : -1;
                                if (TMath::Abs(circle_params[i][3]) < min_R)
                                {
                                    min_R = TMath::Abs(circle_params[i][3]);
                                    min_index = i;
                                }
                            }
                            if (true) 
                            {
                                std::vector<std::vector<float> > keep_one;
                                if (min_index >= 0)
                                {
                                    keep_one.push_back(circle_params[min_index]);
                                }
                                circle_params = keep_one;
                            }
                            
                        }
                        if (charge_flip)
                            circle_params.clear(); // skip this track if charge flip detected
                    }// end of charge flip

                    for (size_t icircle = 0; icircle < circle_params.size(); ++icircle)
                    {
                        const double pt = TMath::Abs(circle_params[0][2] * (0.003 * 0.9));
                        double dca = sqrt(SQR(circle_params[0][0] - xvtx) + SQR(circle_params[0][1] - yvtx)) - TMath::Abs(circle_params[0][2]) - 0.0671 * pow(pt, -0.062) + 0.0690;
                        const double rescattering_resolution =  step_size/sqrt(2) + 0.0130;// + 0*TMath::Sqrt( 0.18 / pt / pt + 0.99 ) / 100; // rescattering
                        //if (pt<0.2||dca>0.05) continue; // reject tracks with pt < 0.2 GeV/c or DCA > 500 μm
                        const int local_weight = 1. / TMath::Sqrt( SQR(rescattering_resolution) + SQR(dca) ); // weight based on DCA and pt
                        n_wtracks += local_weight;
                        //if( circle_params[0][3]>0 ) n_tracks_east+=local_weight;//(int) (pt*1000)
                        //else n_tracks_west+=local_weight;
                        if(!circle_params[0][4]%2==1) n_tracks_notused+=local_weight;
                        if (circle_params[0][4] < 2) n_tracks_pions += local_weight; // pions
                        weight_of_tracks.push_back(local_weight);

                        n_tracks++;
                        used_outer_hits.push_back((int) circle_params[0][5]); // store ihit
                        used_outer_hits.push_back((int) circle_params[0][6]); // store jhit
                    }
                }
                n_tracks = n_wtracks;
                n_tracks_vec.push_back(n_tracks);
                track_vertices.push_back(std::make_pair(xvtx, yvtx));
                if(weight_of_tracks.size()>1) //splitting vecotr in 2 randomly
                {
                    std::srand(unsigned(time(0))); 
                    std::random_shuffle(weight_of_tracks.begin(), weight_of_tracks.end());
                    n_tracks_west = 0; n_tracks_east = 0; 
                    for (size_t i = 0; i < weight_of_tracks.size(); ++i)
                    {
                        if (i % 2 == 0)
                            n_tracks_west += weight_of_tracks[i];
                        else
                            n_tracks_east += weight_of_tracks[i];
                    }
                }
                
                if(n_tracks_west)
                {
                    track_vertices_west.push_back(std::make_pair(xvtx, yvtx));
                    n_tracks_vec_west.push_back(n_tracks_west);
                }
                if (n_tracks_east)
                {
                    track_vertices_east.push_back(std::make_pair(xvtx, yvtx));
                    n_tracks_vec_east.push_back(n_tracks_east);
                }
                if(n_tracks_notused)
                {
                    track_vertices_notused.push_back(std::make_pair(xvtx, yvtx));
                    n_tracks_vec_notused.push_back(n_tracks_notused);
                }
                if(n_tracks_pions)
                {
                    track_vertices_pions.push_back(std::make_pair(xvtx, yvtx));
                    n_tracks_vec_pions.push_back(n_tracks_pions);
                }
                if (n_tracks > max_tracks)
                {
                    max_tracks = n_tracks;
                    best_x = xvtx;
                    best_y = yvtx;
                }
            }// end of loop over yvtx
        }// end of loop over xvtx
        int min_tracks = *std::min_element(n_tracks_vec.begin(), n_tracks_vec.end());
        int track_treshold = (max_tracks - min_tracks) * fraction + min_tracks;
        if (n_tracks_vec.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            for (size_t i = 0; i < n_tracks_vec.size(); ++i)
            {
                if (n_tracks_vec[i] > track_treshold)
                {
                    {
                        hist_vtx_grid_xy->Fill(track_vertices[i].first, track_vertices[i].second, event->GetCentrality(), n_tracks_vec[i]);
                    }
                    const double weight = n_tracks_vec[i] - track_treshold;
                    sum_weights += weight;
                    sum_x += track_vertices[i].first * weight;
                    sum_y += track_vertices[i].second * weight;
                }
            }
            if (fill_hist && TMath::Abs(event->GetPreciseX()-beam_x)<0.0025 && TMath::Abs(event->GetPreciseY()-beam_y)<0.0025)
            {
                for (size_t i = 0; i < n_tracks_vec.size(); ++i)
                {
                    hist_vtx_grid_xy->Fill(track_vertices[i].first, track_vertices[i].second, event->GetCentrality(), n_tracks_vec[i]);
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x = new_vx;
            best_y = new_vy;
        }
        float best_x_west = 1;
        float best_y_west = 1;
        float best_x_east = -1;
        float best_y_east = -1;
        float best_x_notused = -2;
        float best_y_notused = -2;
        float best_x_pions = -3;
        float best_y_pions = -3;
        if(n_tracks_vec_east.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            int max_tracks_east = *std::max_element(n_tracks_vec_east.begin(), n_tracks_vec_east.end());
            int min_tracks_east = *std::min_element(n_tracks_vec_east.begin(), n_tracks_vec_east.end());
            int track_treshold_east = (max_tracks_east - min_tracks_east) * fraction + min_tracks_east;
            for (size_t i = 0; i < n_tracks_vec_east.size(); ++i)
            {
                if (n_tracks_vec_east[i] > track_treshold_east)
                {
                    const double weight = n_tracks_vec_east[i] - track_treshold_east;
                    sum_weights += weight;
                    sum_x += track_vertices_east[i].first * weight;
                    sum_y += track_vertices_east[i].second * weight;
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x_east = new_vx;
            best_y_east = new_vy;
        }
        if(n_tracks_vec_west.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            int max_tracks_west = *std::max_element(n_tracks_vec_west.begin(), n_tracks_vec_west.end());
            int min_tracks_west = *std::min_element(n_tracks_vec_west.begin(), n_tracks_vec_west.end());
            int track_treshold_west = (max_tracks_west - min_tracks_west) * fraction + min_tracks_west;
            for (size_t i = 0; i < n_tracks_vec_west.size(); ++i)
            {
                if (n_tracks_vec_west[i] > track_treshold_west)
                {
                    const double weight = n_tracks_vec_west[i] - track_treshold_west;
                    sum_weights += weight;
                    sum_x += track_vertices_west[i].first * weight;
                    sum_y += track_vertices_west[i].second * weight;
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x_west = new_vx;
            best_y_west = new_vy;
        }
        if(n_tracks_vec_notused.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            int max_tracks_notused = *std::max_element(n_tracks_vec_notused.begin(), n_tracks_vec_notused.end());
            int min_tracks_notused = *std::min_element(n_tracks_vec_notused.begin(), n_tracks_vec_notused.end());
            int track_treshold_notused = (max_tracks_notused - min_tracks_notused) * fraction + min_tracks_notused;
            for (size_t i = 0; i < n_tracks_vec_notused.size(); ++i)
            {
                if (n_tracks_vec_notused[i] > track_treshold_notused)
                {
                    const double weight = n_tracks_vec_notused[i] - track_treshold_notused;
                    sum_weights += weight;
                    sum_x += track_vertices_notused[i].first * weight;
                    sum_y += track_vertices_notused[i].second * weight;
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x_notused = new_vx;
            best_y_notused = new_vy;
        }
        if(n_tracks_vec_pions.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            int max_tracks_pions = *std::max_element(n_tracks_vec_pions.begin(), n_tracks_vec_pions.end());
            int min_tracks_pions = *std::min_element(n_tracks_vec_pions.begin(), n_tracks_vec_pions.end());
            int track_treshold_pions = (max_tracks_pions - min_tracks_pions) * fraction + min_tracks_pions;
            for (size_t i = 0; i < n_tracks_vec_pions.size(); ++i)
            {
                if (n_tracks_vec_pions[i] > track_treshold_pions)
                {
                    const double weight = n_tracks_vec_pions[i] - track_treshold_pions;
                    sum_weights += weight;
                    sum_x += track_vertices_pions[i].first * weight;
                    sum_y += track_vertices_pions[i].second * weight;
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x_pions = new_vx;
            best_y_pions = new_vy;
        }

        if (verbosity)
            std::cout << "\033[32mBest X and Y vertex = " << best_x << " " << best_y << " with " << max_tracks << " " << min_tracks << " " << track_treshold << " tracks\033[0m" << " previous " << event->GetPreciseX() << " " << event->GetPreciseY() << std::endl;
        if (verbosity)
            std::cout << "\033[32mBest X and Y vertex (east-west) = " << best_x_east<< " "<< best_x_west << " " << best_y_east<<" "<<best_y_west << " tracks\033[0m" << std::endl;
        if (verbosity)
            std::cout << "\033[32mBest X and Y vertex (notused) = " << best_x_notused<< " "<< best_y_notused << " " << best_x_pions << " " << best_y_pions << " tracks\033[0m" << std::endl;
        if (max_tracks < min_track_count)
            return; // Not enough hits to form a vertex

        if (fill_hist)
        {
            hist_vtx_delta_x->Fill(best_x_west-best_x_east, best_x, event->GetCentrality());
            hist_vtx_delta_y->Fill(best_y_west-best_y_east, best_y, event->GetCentrality());
            hist_vtx_x->Fill(best_x_notused - event->GetPreciseX(), best_x, event->GetCentrality());
            hist_vtx_y->Fill(best_y_notused - event->GetPreciseY(), best_y, event->GetCentrality());
            hist_vtx_grid_xy->Fill(best_x_notused, best_y_notused, event->GetCentrality());
            hist_vtx_z->Fill(best_x_notused, best_y_notused, beam_z);
            hist_vtx_delta_x_reuse->Fill(best_x_notused-best_x, best_x, event->GetCentrality());
            hist_vtx_delta_y_reuse->Fill(best_y_notused-best_y, best_y, event->GetCentrality());
            if(event->GetNeleccand()==1)
            {
                hist_vtx_delta_x_pion->Fill(best_x_pions-best_x, best_x, event->GetCentrality());
                hist_vtx_delta_y_pion->Fill(best_y_pions-best_y, best_y, event->GetCentrality());
            }
        }

        event->SetBBCchargeN(best_x_pions);
        event->SetBBCchargeS(best_y_pions);
        if(n_tracks_vec_notused.size()==0) return; // if we have not used tracks, we do not set the vertex
        event->SetPreciseX(best_x_notused);
        event->SetPreciseY(best_y_notused);

        if (fill_hist)
        {   
            int n_used_tracks = 0, n_used_tracks_notused = 0; 
            for (std::map<int, std::vector<std::pair<int, int> > >::iterator it = khit_to_chain.begin();
                 it != khit_to_chain.end(); ++it)
            {
                const unsigned int khit = it->first;
                const std::vector<std::pair<int, int> > &layer12_hits = it->second;
                std::vector<std::vector<float> > circle_params; // stores (cx, cy, R)
                for (unsigned int ichain = 0; ichain < layer12_hits.size(); ++ichain)
                {
                    const unsigned int ihit = layer12_hits[ichain].first;
                    const unsigned int jhit = layer12_hits[ichain].second;

                    MyDileptonAnalysis::MyVTXHit *hit0 = event->GetVTXHitEntry(ihit);
                    if (hit0->GetLayer() < 2)
                        continue;

                    MyDileptonAnalysis::MyVTXHit *hit1 = event->GetVTXHitEntry(jhit);
                    if (hit1->GetLayer() != 1)
                        continue;

                    float phi0 = hit0->GetPhiHit(best_x_notused, best_y_notused, beam_z);
                    float the0 = hit0->GetTheHit(best_x_notused, best_y_notused, beam_z);
                    float phi1 = hit1->GetPhiHit(best_x_notused, best_y_notused, beam_z);
                    float the1 = hit1->GetTheHit(best_x_notused, best_y_notused, beam_z);

                    float dphi = phi1 - phi0;
                    float dtheta = the1 - the0;

                    float x0 = hit0->GetXHit();
                    float y0 = hit0->GetYHit();
                    float x1 = hit1->GetXHit();
                    float y1 = hit1->GetYHit();

                    MyDileptonAnalysis::MyVTXHit *hit2 = event->GetVTXHitEntry(khit);
                    if (hit2->GetLayer() != 0)
                        continue;

                    float phi2 = hit2->GetPhiHit(best_x_notused, best_y_notused, beam_z);
                    float the2 = hit2->GetTheHit(best_x_notused, best_y_notused, beam_z);
                    float x2 = hit2->GetXHit();
                    float y2 = hit2->GetYHit();

                    float r0 = sqrt(SQR(x2) + SQR(y2));
                    float r1 = sqrt(SQR(x1) + SQR(y1));
                    float r01 = r0 - r1;
                    float r12 = r1 - sqrt(SQR(x0) + SQR(y0));

                    float dphi1 = phi2 - (phi1 + dphi * r01 / r12);
                    float dthe1 = the2 - the1; // - ( the0 + dtheta * z02 / z01 );

                    float x12 = x1 - x0;
                    float y12 = y1 - y0;
                    float x23 = x2 - x1;
                    float y23 = y2 - y1;
                    float det = x12 * y23 - y12 * x23;
                    if (fabs(det) < 1e-6)
                        continue; // collinear points

                    float A = x0 * x0 + y0 * y0;
                    float B = x1 * x1 + y1 * y1;
                    float C = x2 * x2 + y2 * y2;

                    float cx = (A * (y1 - y2) + B * (y2 - y0) + C * (y0 - y1)) / (2 * det);
                    float cy = (A * (x2 - x1) + B * (x0 - x2) + C * (x1 - x0)) / (2 * det);
                    float R = sqrt((x0 - cx) * (x0 - cx) + (y0 - cy) * (y0 - cy));
                    float pt = R * (0.003 * 0.9);

                    float angle_beam = atan2(best_y_notused - cy, best_x_notused - cx);
                    float x_proj = cx + R * cos(angle_beam);
                    float y_proj = cy + R * sin(angle_beam);
                    float dca = sqrt(SQR(cx - best_x_notused) + SQR(cy - best_y_notused)) - R;//  - 0.0671 * pow(pt, -0.062) + 0.0690;
                    float dca_x = best_x_notused - x_proj;
                    float dca_y = best_y_notused - y_proj;
                    
                    if (fabs(dthe1) < sddthe && fill_hist)
                        vtx_dphi_dphi_hist->Fill(dphi1, dphi, pt);
                    if (fabs(dphi1) < sddphi && fill_hist)
                        vtx_dthe_dthe_hist->Fill(dthe1, dtheta, pt);
                        
                    if (fabs(dphi) > sdphi || fabs(dtheta) > sdthe)
                        continue;
                    if (fabs(dphi1) > sddphi || fabs(dthe1) > sddthe )
                        continue;
                    n_used_tracks++;
                    if(!(((hit0->GetLadder()>24&&hit0->GetLadder()<48)) || (hit1->GetLadder()>24&&hit1->GetLadder()<48) || (hit2->GetLadder()>24&&hit2->GetLadder()<48))) n_used_tracks_notused++;

                    if ( TMath::Abs(hit2->GetZHit())>12) continue; // skip if z is too large add centrality
                    int z_bin =  ( (int) ((hit2->GetZHit()+12)/3) )*20 + ( (int) ((phi2 + 1.5) / 0.3) ) + 160*( (int) (event->GetCentrality()>20 ? 1:0) ); // 0-20 for negative z, 20-40 for positive z
                    hist_dca_x->Fill(dca_x, pt, z_bin);
                    hist_dca_y->Fill(dca_y, pt, z_bin);
                    const int charge = phi0 > phi2 ? 1 : -1;
                    if (hit2->GetZHit()>0)  conv_photon_mass_hist->Fill(dca*charge, phi2, pt*charge);
                    else pi0_mass_hist->Fill(dca*charge, phi2, pt*charge);
                }
            }
            if(verbosity) std::cout << "\033[32mUsed tracks = " << n_used_tracks << " " << n_used_tracks_notused << "\033[0m" << std::endl;
            hits_vtx_ntracks->Fill(n_used_tracks, event->GetCentrality());
            hits_vtx_ntracks_ofnotusedhits->Fill(n_used_tracks_notused, event->GetCentrality());
        }
    }

    void MyEventContainer::ConversionFinder(int fill_hist = 0, int verbosity = 0)
    {
        const int Ntracks = event->GetNtrack();
        const int central_bin = event->GetCentrality()/20;
        for (int itrk = 0; itrk < Ntracks; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            const int charge = mytrk->GetChargePrime();
            const float pt = mytrk->GetPtPrime();
            //if(charge==1)  continue;
            if(!mytrk->GetHitCounter(0)) continue;
            const int layer0_hit_id = mytrk->GetHitIndex(0);
            MyDileptonAnalysis::MyVTXHit *layer0_hit = event->GetVTXHitEntry(layer0_hit_id);
            const float phi00 = mytrk->GetPhi0();
            const float phi0 = layer0_hit->GetPhiHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ())+0.001;
            const float the0 = layer0_hit->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
            const float x0 = layer0_hit->GetXHit();
            const float y0 = layer0_hit->GetYHit();
            const float r0 = sqrt(SQR(x0-event->GetPreciseX()) + SQR(y0-event->GetPreciseY())); 
            const float dphi_dr = (phi0 - phi00) / r0;
            int is_conversion = 0, is_dalitz = 0;
            for (int ihit2 = 0; ihit2 < event->GetNVTXhit(); ihit2++)
            {
                MyDileptonAnalysis::MyVTXHit *layer1_hit = event->GetVTXHitEntry(ihit2);
                if (layer1_hit->GetLayer() == 0 || layer1_hit->GetLayer() == 3)
                    continue;
                if ((layer1_hit->GetLadder()>24&&layer1_hit->GetLadder()<48)) continue;
                const int inner_layer = layer1_hit->GetLayer();
                const float phi1 = layer1_hit->GetPhiHit(x0, y0, event->GetPreciseZ());
                const float the1 = layer1_hit->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                const float dphi = phi1 - phi0; ///dphi = |f(R,pt)|*charge => dphi*charge -> psotive for track but negative for second tracks
                const float cdphi = charge*dphi;
                const float dthe = the1 - the0;

                if (cdphi < 0 && cdphi> -0.1 && TMath::Abs(dthe) < 0.03 )///make dthe in layer 1 to 0.001
                {
                    for (int ihit3 = 0; ihit3 < event->GetNVTXhit(); ihit3++)
                    {
                        MyDileptonAnalysis::MyVTXHit *layer2_hit = event->GetVTXHitEntry(ihit3);
                        if (layer2_hit->GetLayer() <= inner_layer)
                            continue;
                        if ((layer2_hit->GetLadder()>24&&layer2_hit->GetLadder()<48)) continue;
                        const float phi2 = layer2_hit->GetPhiHit(x0, y0, event->GetPreciseZ());
                        const float the2 = layer2_hit->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                        const float r10 = sqrt(SQR(layer1_hit->GetXHit()-x0) + SQR(layer1_hit->GetYHit()-y0));
                        const float r20 = sqrt(SQR(layer2_hit->GetXHit()-x0) + SQR(layer2_hit->GetYHit()-y0));
                        const float dphi1 = phi2 - phi0;
                        if(dphi1*charge > 0) continue;
                        const float dphihh = phi2 - (phi0 + dphi_dr*r20);
                        if(TMath::Abs(dphihh) < 0.005) continue;

                        const float dphi2 = dphi1 - dphi * r20 / r10;
                        const float dthe2 = the2 - the1;
                        if (fill_hist)
                        {
                            const int layer_bin = (inner_layer>1 ? 1 : 0) + (layer2_hit->GetLayer()-2) + (central_bin>2 ? 3:0);
                            if(TMath::Abs(dthe2) < 0.01) hist_conv_phi_phi[layer_bin]->Fill(dphi2, dphi, pt);
                            if(TMath::Abs(dphi2) < 0.10) hist_conv_the_the[layer_bin]->Fill(dthe2, dthe, pt);
                        }
                        if ( TMath::Abs(dphi2) < 0.1 && TMath::Abs(dthe2) < 0.01 )
                        {
                            is_conversion++;
                            if (verbosity)
                            {
                                std::cout << "\033[32mFound conversion\033[0m" << std::endl;
                                std::cout << "\033[32m" << dphi << " " << dphi2 << " " << dthe2 << " " << inner_layer << " " << layer2_hit->GetLayer() << " " << mytrk->GetPtPrime() << " " << mytrk->GetTOFDPHI() << "\033[0m" << std::endl;
                            }
                        }
                    }
                }
            }
            for (int ihit0 = 0; ihit0 < event->GetNVTXhit(); ihit0++)
            {
                MyDileptonAnalysis::MyVTXHit *layer0_hit2 = event->GetVTXHitEntry(ihit0);
                if (layer0_hit2->GetLayer() == 3)
                    continue;
                if ((layer0_hit2->GetLadder()>24&&layer0_hit2->GetLadder()<48)) continue;
                const int inner_layer = layer0_hit2->GetLayer();
                const float phi1 = layer0_hit2->GetPhiHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                const float the1 = layer0_hit2->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                const float dphi = phi1 - phi00; ///dphi = |f(R,pt)|*charge => dphi*charge -> positive for track but negative for second tracks
                const float cdphi = charge*dphi;
                const float dthe = the1 - the0;
                if (cdphi < -0.002 && cdphi > -0.1 && TMath::Abs(dthe) < 0.03 )
                {
                    for (int ihit1 = 0; ihit1 < event->GetNVTXhit(); ihit1++)
                    {
                        MyDileptonAnalysis::MyVTXHit *layer1_hit = event->GetVTXHitEntry(ihit1);
                        if (layer1_hit->GetLayer() <= inner_layer)
                            continue;
                        if ((layer1_hit->GetLadder()>24&&layer1_hit->GetLadder()<48)) continue;
                        const float phi2 = layer1_hit->GetPhiHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                        const float dphi1 = phi2 - phi00;
                        if(dphi1*charge>-0.002||(phi2-phi1)*charge>-0.002) continue;
                        const float the2 = layer1_hit->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                        const float r10 = sqrt(SQR(layer0_hit2->GetXHit()-event->GetPreciseX()) + SQR(layer0_hit2->GetYHit()-event->GetPreciseY()));
                        const float r20 = sqrt(SQR(layer1_hit ->GetXHit()-event->GetPreciseX()) + SQR(layer1_hit ->GetYHit()-event->GetPreciseY()));
                        const float dphi1hh = phi2 - (phi00 + dphi_dr*r20);
                        if(TMath::Abs(dphi1hh) < 0.005) continue;
                        const float dphi2 = dphi1 - dphi * r20 / r10;
                        const float dthe2 = the2 - the1;
                        if (fill_hist)
                        {
                            const int layer_bin = (inner_layer==0?0:(inner_layer==1?2:3)) + (layer1_hit->GetLayer()-1) + (central_bin>2 ? 6:0);
                            if(TMath::Abs(dthe2) < 0.01) hist_daltz_phi_phi[layer_bin]->Fill(dphi2, dphi, pt);
                            if(TMath::Abs(dphi2) < 0.10) hist_daltz_the_the[layer_bin]->Fill(dthe2, dthe, pt);
                        }
                        if ( TMath::Abs(dphi2) < 0.02 && TMath::Abs(dthe2) < 0.01 )
                        {
                            is_dalitz++;
                            if (verbosity)
                            {
                                std::cout << "\033[32mFound Dalitz\033[0m" << std::endl;
                                std::cout << "\033[32m" << dphi << " " << dphi2 << " " << dthe2 << " " << inner_layer << " " << layer1_hit->GetLayer() << " " << mytrk->GetPtPrime() << " " << mytrk->GetTOFDPHI() << "\033[0m" << std::endl;
                            }
                        }
                    }
                }
            }
            for (int ihit0 = 0; ihit0 < event->GetNVTXhit(); ihit0++)
            {
                MyDileptonAnalysis::MyVTXHit *layer0_hit2 = event->GetVTXHitEntry(ihit0);
                if (layer0_hit2->GetLayer() != 0)
                    continue;
                //if ((layer0_hit2->GetLadder()>24&&layer0_hit2->GetLadder()<48)) continue;
                const int inner_layer = layer0_hit2->GetLayer();
                const float phi1 = layer0_hit2->GetPhiHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                const float the1 = layer0_hit2->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                const float dphi = phi1 - phi00; ///dphi = |f(R,pt)|*charge => dphi*charge -> positive for track but negative for second tracks
                const float cdphi = charge*dphi;
                const float dthe = the1 - the0;
                if (cdphi < 0.1 && cdphi > -0.2 && TMath::Abs(dthe) < 0.1 )
                {
                    for (int ihit1 = 0; ihit1 < event->GetNVTXhit(); ihit1++)
                    {
                        MyDileptonAnalysis::MyVTXHit *layer1_hit = event->GetVTXHitEntry(ihit1);
                        if (layer1_hit->GetLayer() == 0 || layer1_hit->GetLayer() == 3)
                            continue;
                        const int layer2 = layer1_hit->GetLayer();
                        if ((layer1_hit->GetLadder()>24&&layer1_hit->GetLadder()<48)) continue;
                        const float phi2 = layer1_hit->GetPhiHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                        const float dphi1 = phi2 - phi1;
                        const float the2 = layer1_hit->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                        const float dthe2 = the2 - the1;
                        if(dphi1*charge > -0.002 || dphi1*charge < -0.2 || TMath::Abs(dthe2)>0.01 ) continue;
                        for (int ihit2 = 0; ihit2 < event->GetNVTXhit(); ihit2++)
                        {
                            MyDileptonAnalysis::MyVTXHit *layer2_hit = event->GetVTXHitEntry(ihit2);
                            if (layer2_hit->GetLayer() <= layer2)
                                continue;
                            if ((layer2_hit->GetLadder()>24&&layer2_hit->GetLadder()<48)) continue;
                            const float phi3 = layer2_hit->GetPhiHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                            const float dphi2 = phi3 - phi1;
                            if (dphi2 * charge > -0.002)
                                continue;
                            const float the3 = layer2_hit->GetTheHit(event->GetPreciseX(), event->GetPreciseY(), event->GetPreciseZ());
                            const float r00 = sqrt(SQR(layer0_hit2->GetXHit() - event->GetPreciseX()) + SQR(layer0_hit2->GetYHit() - event->GetPreciseY()));
                            const float r10 = sqrt(SQR(layer1_hit->GetXHit() - event->GetPreciseX()) + SQR(layer1_hit->GetYHit() - event->GetPreciseY()))-r00;
                            const float r20 = sqrt(SQR(layer2_hit->GetXHit() - event->GetPreciseX()) + SQR(layer2_hit->GetYHit() - event->GetPreciseY()))-r00;
                            const float dphi1hh = phi2 - (phi00 + dphi_dr * r20);
                            if (TMath::Abs(dphi1hh) < 0.005)
                                continue;
                            const float dphi3 = dphi2 - dphi1 * r20 / r10;
                            const float dthe3 = the3 - the2;
                            if (false)
                            {
                                const int layer_bin = (inner_layer == 0 ? 0 : (inner_layer == 1 ? 2 : 3)) + (layer1_hit->GetLayer() - 1) + (central_bin > 2 ? 6 : 0);
                                if (TMath::Abs(dthe3) < 0.01)
                                    hist_daltz_phi_phi[layer_bin]->Fill(dphi3, dphi, pt);
                                if (TMath::Abs(dphi3) < 0.10)
                                    hist_daltz_the_the[layer_bin]->Fill(dthe3, dthe, pt);
                            }
                            if (TMath::Abs(dphi3) < 0.02 && TMath::Abs(dthe3) < 0.01)
                            {
                                is_dalitz++;
                                if (verbosity)
                                {
                                    std::cout << "\033[32mFound Dalitz\033[0m" << std::endl;
                                    std::cout << "\033[32m" << dphi << " " << dphi2 << " " << dthe2 << " " << inner_layer << " " << layer1_hit->GetLayer() << " " << mytrk->GetPtPrime() << " " << mytrk->GetTOFDPHI() << "\033[0m" << std::endl;
                                }
                            }
                        }
                    }
                }
            }
            if(mytrk->GetMinsDphi(0)+mytrk->GetMinsDphi(1) < -2) is_conversion += 10;
            mytrk->SetEmcdphi_e(is_dalitz*100+is_conversion);
            if(fill_hist)
            {
                sdphi_conv_hist[central_bin]->Fill(mytrk->GetMinsDphi(0)+mytrk->GetMinsDphi(1), mytrk->GetMinsDphi(0), pt);
                if (is_conversion>1) 
                    sdphi_real_conv_hist[central_bin]->Fill(mytrk->GetMinsDphi(0)+mytrk->GetMinsDphi(1), mytrk->GetMinsDphi(0), pt);
                hist_is_dalitz_conv[central_bin]->Fill(is_conversion, is_dalitz, pt);
                hist_is_ml_conv[central_bin]->Fill(is_conversion, TMath::Log10 ( mytrk->GetTOFDPHI()>1 ?  mytrk->GetTOFDPHI() : 1 ), pt);
            }
        }
    }

    void MyEventContainer::FillQAHist(const int mc_id)
    {
        const int Nelectrons = event->GetNtrack();
        for (int i = 0; i < Nelectrons; i++)
        {
            MyDileptonAnalysis::MyElectron *electron = event->GetEntry(i);
            int charge_centr_bin = event->GetCentrality() + 50 * (1 - electron->GetChargePrime());

            if (electron->GetEcore()/electron->GetPtot()<0.6 || electron->GetN0()<0 ) continue;
            //if ( electron->GetEcore()/electron->GetPtot() < 0.5 || electron->GetN0() < 0 )
            //    continue;
            //if (fabs(electron->GetEmcTOF())>5 && mc_id != -99) continue; 

            int central_bin = (int) event->GetCentrality() / 20;
            central_bin += N_centr * ( 1 - electron->GetChargePrime() ) / 2;

            //const float eConv = std::log10(electron->GetNHits()+1)*5 + std::log10(electron->GetTOFDPHI()+1);
            int hit_assocaition = 0;
            if ( (((TMath::Abs(electron->GetMinsDphi(3))<3) ||
                   (TMath::Abs(electron->GetMinsDphi(2))<3) ) && 
                   (TMath::Abs(electron->GetMinsDphi(1))<3) && 
                   (electron->GetMinsDphi(0))>-3 ) ) hit_assocaition=1;
            if (hit_assocaition==1 && electron->GetMinsDphi(0)<3 ) hit_assocaition=2;
            if ( (((TMath::Abs(electron->GetMinsDphi(3))<2 && TMath::Abs(electron->GetMinsDthe(3))<2) ||
                   (TMath::Abs(electron->GetMinsDphi(2))<2 && TMath::Abs(electron->GetMinsDthe(2))<2) ) && 
                   (TMath::Abs(electron->GetMinsDphi(1))<2) && 
                   (electron->GetMinsDphi(0)> -1 ) )) hit_assocaition=3;
            if (hit_assocaition==3 && electron->GetMinsDphi(0)<2 ) hit_assocaition=4;
            int conv_reject = 0;
            if ( ((int)electron->GetEmcdphi_e())%100==0) conv_reject=5;
            if ( ((int)electron->GetEmcdphi_e())%100<3 && ((int)electron->GetEmcdphi_e())/100<3) conv_reject=10;
            if ( ((int)electron->GetEmcdphi_e())%100<1 && ((int)electron->GetEmcdphi_e())/100<3) conv_reject=15;
            if ( ((int)electron->GetEmcdphi_e())%100<1 && ((int)electron->GetEmcdphi_e())/100<1) conv_reject=20;
            const float eConv = hit_assocaition + conv_reject;

            if (electron->GetMcId()%10>0 ) 
                el_pt_hist[central_bin]->Fill(electron->GetPtPrime(),0.,eConv);
            if (electron->GetMcId()%10>5 ) 
                el_pt_hist[central_bin]->Fill(electron->GetPtPrime(),1.,eConv);
            if (electron->GetMcId()>90) 
                el_pt_hist[central_bin]->Fill(electron->GetPtPrime(),2.,eConv);
            if (electron->GetMcId()>900) 
                el_pt_hist[central_bin]->Fill(electron->GetPtPrime(),3.,eConv);
            if (electron->GetMcId()>9000 && electron->GetMcId()%10>5) 
                el_pt_hist[central_bin]->Fill(electron->GetPtPrime(),4.,eConv);

            //if (electron->GetMcId()<100) continue;///figuring out how bdt actually works

            ep_hist->Fill(electron->GetEcore()/electron->GetPtot(),electron->GetPtPrime(), charge_centr_bin);
            n0_hist->Fill(electron->GetN0(),electron->GetPtPrime(), charge_centr_bin);
            prob_hist->Fill(electron->GetProb(),electron->GetPtPrime(), charge_centr_bin);
            disp_hist->Fill(electron->GetDisp(),electron->GetPtPrime(), charge_centr_bin);
            chi2npe0_hist->Fill(electron->GetChi2()/(electron->GetNpe0()+0.1),electron->GetPtPrime(), charge_centr_bin);

            if (electron->GetMcId()>99 && electron->GetMcId()%10>5)///figuring out how bdt actually works
            {

                temc->Fill(electron->GetEmcTOF(),electron->GetPtPrime(),charge_centr_bin);
                ttof->Fill(electron->GetTOFE()*0.01,electron->GetPtPrime(),charge_centr_bin);
    
                ep_hist_el->Fill(electron->GetEcore()/electron->GetPtot(),electron->GetProb(),electron->GetPtPrime());
                n0_hist_el->Fill(electron->GetN0(),electron->GetDisp(),event->GetCentrality());
                prob_hist_el->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetDisp(),event->GetCentrality());
                disp_hist_el->Fill(electron->GetDisp(),electron->GetNpe0(),event->GetCentrality());
                chi2npe0_hist_el->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetNpe0(),event->GetCentrality());
                rich_prob1->Fill(electron->GetChi2()/electron->GetNpe0(),electron->GetN0()-1*electron->GetDisp(),event->GetCentrality());
                rich_prob2->Fill(electron->GetNpe0(),electron->GetN0()-1*electron->GetDisp(),event->GetCentrality());
                rich_prob3->Fill(electron->GetEmcdphi(),electron->GetEmcdz(),event->GetCentrality());

                const float Rghost = sqrt(SQR(electron->GetEmcdphi())+SQR(electron->GetEmcdz()));
                el_had_dr->Fill(Rghost,electron->GetPtPrime(),charge_centr_bin);
                el_had_dphi->Fill(electron->GetEmcdphi(),electron->GetPtPrime(),charge_centr_bin);
                el_had_dz  ->Fill(electron->GetEmcdz()  ,electron->GetPtPrime(),charge_centr_bin);
            }

            if(electron->GetEmcTOF()>0.4 && electron->GetEmcTOF()<1.4) 
            {
                charge_centr_bin += 400;
                
                ep_hist->Fill(electron->GetEcore()/electron->GetPtot(),electron->GetPtPrime(), charge_centr_bin);
                n0_hist->Fill(electron->GetN0(),electron->GetPtPrime(), charge_centr_bin);
                prob_hist->Fill(electron->GetProb(),electron->GetPtPrime(), charge_centr_bin);
                disp_hist->Fill(electron->GetDisp(),electron->GetPtPrime(), charge_centr_bin);
                chi2npe0_hist->Fill(electron->GetChi2()/(electron->GetNpe0()+0.1),electron->GetPtPrime(), charge_centr_bin);

                const float Rghost = sqrt(SQR(electron->GetEmcdphi())+SQR(electron->GetEmcdz()));
                el_had_dr->Fill(Rghost,electron->GetPtPrime(),charge_centr_bin);
                el_had_dphi->Fill(electron->GetEmcdphi(),electron->GetPtPrime(),charge_centr_bin);
                el_had_dz  ->Fill(electron->GetEmcdz()  ,electron->GetPtPrime(),charge_centr_bin);
            }
        }
    }

    void MyEventContainer::FillQAHistPreAssoc()
    {
        const int Nelectrons = event->GetNtrack();
        for (int i = 0; i < Nelectrons; i++)
        {
            MyDileptonAnalysis::MyElectron *electron = event->GetEntry(i);
            int charge_centr_bin = event->GetCentrality() + 50 * (1 - electron->GetChargePrime())+200;

            if (electron->GetEcore()/electron->GetPtot()<0.6 || electron->GetN0()<0 ) continue;

            temc->Fill(electron->GetEmcTOF(),electron->GetPtPrime(),charge_centr_bin+200);
            ttof->Fill(electron->GetTOFE()*0.01,electron->GetPtPrime(),charge_centr_bin+200);
            
            if (electron->GetN0()>= 2 +SQR(electron->GetDisp())/8. && electron->GetDisp()<4 && electron->GetChi2()/(electron->GetNpe0()+0.1)<10 && electron->GetProb()>0.01 && electron->GetChi2()>0) 
                ep_hist->Fill(electron->GetEcore()/electron->GetPtot(),electron->GetPtPrime(), charge_centr_bin);
            if (electron->GetEcore()/electron->GetPtot()>0.8 && electron->GetDisp()<4 && electron->GetChi2()/(electron->GetNpe0()+0.1)<10 && electron->GetProb()>0.01 && electron->GetChi2()>0) 
                n0_hist->Fill(electron->GetN0(),electron->GetPtPrime(), charge_centr_bin);
            if (electron->GetN0()>= 2 +SQR(electron->GetDisp())/8. && electron->GetEcore()/electron->GetPtot()>0.8 && electron->GetDisp()<4 && electron->GetChi2()/(electron->GetNpe0()+0.1)<10 && electron->GetChi2()>0) 
                prob_hist->Fill(electron->GetProb(),electron->GetPtPrime(), charge_centr_bin);
            if (electron->GetN0()>= 4. && electron->GetEcore()/electron->GetPtot()>0.8 && electron->GetChi2()/(electron->GetNpe0()+0.1)<10 && electron->GetProb()>0.01 && electron->GetChi2()>0) 
                disp_hist->Fill(electron->GetDisp(),electron->GetPtPrime(), charge_centr_bin);
            if (electron->GetN0()>= 2 +SQR(electron->GetDisp())/8. && electron->GetEcore()/electron->GetPtot()>0.8 && electron->GetDisp()<4 && electron->GetProb()>0.01 && electron->GetChi2()>0) 
                chi2npe0_hist->Fill(electron->GetChi2()/(electron->GetNpe0()+0.1),electron->GetPtPrime(), charge_centr_bin);

            if (electron->GetMcId()<1000 || electron->GetProb()<0.1) continue;///figuring out how bdt actually works

            temc->Fill(electron->GetEmcTOF(),electron->GetPtPrime(),charge_centr_bin);
            ttof->Fill(electron->GetTOFE()*0.01,electron->GetPtPrime(),charge_centr_bin);

            const float Rghost = sqrt(SQR(electron->GetEmcdphi())+SQR(electron->GetEmcdz()));
            el_had_dr->Fill(Rghost,electron->GetPtPrime(),charge_centr_bin);
            el_had_dphi->Fill(electron->GetTOFDPHI(),electron->GetPtPrime(),charge_centr_bin);
            el_had_dz  ->Fill(electron->GetTOFDZ()  ,electron->GetPtPrime(),charge_centr_bin);
        }
    }

    int MyEventContainer::isGhostEvent()
    {   
        int n_ghosts = 0;
        int pair_counter[3][33] = {{0},{0},{0}};
        for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk1 = event->GetEntry(itrk);
            
            for (int jtrk = itrk+1; jtrk < event->GetNtrack(); jtrk++)
            {
                MyDileptonAnalysis::MyElectron *mytrk2 = event->GetEntry(jtrk);
                //if (mytrk1->GetChargePrime() != mytrk2->GetChargePrime()) continue;
                bool isGhost = false;
                const float dcenter_phi = ( mytrk1->GetCrkphi() - mytrk2->GetCrkphi()) / 0.013;
                if ( TMath::Abs(dcenter_phi) < 5) 
                    isGhost = true;

                const float phi_pip = mytrk1->GetPhiDC();
                const float phi_pim = mytrk2->GetPhiDC();

                const float alpha_pip = mytrk1->GetAlphaPrime();
                const float alpha_pim = mytrk2->GetAlphaPrime();

                const float zed_pip = mytrk1->GetZDC();
                const float zed_pim = mytrk2->GetZDC();

                const float dalpha = alpha_pip - alpha_pim;
                const float dphi = phi_pip - phi_pim;
                const float dzed = zed_pip - zed_pim;
                 // pc1
                if (TMath::Abs(dzed) < 6.0 && TMath::Abs(dphi - (0.13 * dalpha)) < 0.015) 
                    isGhost = true;
                // x1x2_1
                if (TMath::Abs(dphi - (0.04 * dalpha)) < 0.015) 
                    isGhost = true;
                // x1x2_2
                if (TMath::Abs(dphi - (-0.065 * dalpha)) < 0.015)
                    isGhost = true;

                const int Nsect1 = mytrk1->GetSect()+mytrk1->GetArm()*4; 
                const int Ysect1 = mytrk1->GetYsect(); 
                const int Zsect1 = mytrk1->GetZsect(); 
                const int Nsect2 = mytrk2->GetSect()+mytrk2->GetArm()*4; 
                const int Ysect2 = mytrk2->GetYsect(); 
                const int Zsect2 = mytrk2->GetZsect();
                
                if ( Nsect1 == Nsect2 && TMath::Abs( Ysect1 - Ysect2 ) < 3 && TMath::Abs( Zsect1 - Zsect2) <3 )
                    isGhost = true;
                
                if (isGhost) {
                    mytrk1->SetGhost(-9);
                    mytrk2->SetGhost(-9);
                    //if (mytrk1->GetMcId() > mytrk2->GetMcId() ) mytrk2->SetGhost(-9);
                    //else if (mytrk1->GetMcId() < mytrk2->GetMcId() ) mytrk1->SetGhost(-9);
                    //else {
                    //    if (mytrk1->GetPtPrime() > mytrk2->GetPtPrime() ) mytrk2->SetGhost(-9);
                    //    else mytrk1->SetGhost(-9);
                    //}
                }
                if(isGhost) n_ghosts++;
                else{
                    const int charge_bin = 1 + (mytrk1->GetChargePrime() + mytrk2->GetChargePrime())/2;
                    if ( mytrk1->GetHitCounter(0)>0 && mytrk2->GetHitCounter(0)>0 ) pair_counter[charge_bin][0]++;
                    if( mytrk1->GetTOFDPHI()>9      && mytrk2->GetTOFDPHI()>9    && mytrk1->GetGhost()!=-9  && mytrk2->GetGhost()!=-9  ) pair_counter[charge_bin][1]++;
                    if( mytrk1->GetTOFDPHI()>99     && mytrk2->GetTOFDPHI()>99   && mytrk1->GetGhost()!=-9  && mytrk2->GetGhost()!=-9  ) pair_counter[charge_bin][2]++;
                    if( mytrk1->GetTOFDPHI()>999    && mytrk2->GetTOFDPHI()>999  && mytrk1->GetGhost()!=-9  && mytrk2->GetGhost()!=-9  ) pair_counter[charge_bin][3]++;
                    if( mytrk1->GetTOFDPHI()>9999   && mytrk2->GetTOFDPHI()>9999 && mytrk1->GetGhost()!=-9  && mytrk2->GetGhost()!=-9  ) pair_counter[charge_bin][4]++;
                    if( mytrk1->GetTOFDPHI()>9      && mytrk2->GetTOFDPHI()>9    && mytrk1->GetNHits()>9    && mytrk2->GetNHits()>9    ) pair_counter[charge_bin][5]++;
                    if( mytrk1->GetTOFDPHI()>9      && mytrk2->GetTOFDPHI()>9    && mytrk1->GetNHits()>99   && mytrk2->GetNHits()>99   ) pair_counter[charge_bin][6]++;
                    if( mytrk1->GetTOFDPHI()>9      && mytrk2->GetTOFDPHI()>9    && mytrk1->GetNHits()>999  && mytrk2->GetNHits()>999  ) pair_counter[charge_bin][7]++;
                    if( mytrk1->GetTOFDPHI()>9      && mytrk2->GetTOFDPHI()>9    && mytrk1->GetNHits()>9999 && mytrk2->GetNHits()>9999 ) pair_counter[charge_bin][8]++;
                    if( mytrk1->GetTOFDPHI()>99     && mytrk2->GetTOFDPHI()>99   && mytrk1->GetNHits()>9    && mytrk2->GetNHits()>9    ) pair_counter[charge_bin][9]++;
                    if( mytrk1->GetTOFDPHI()>99     && mytrk2->GetTOFDPHI()>99   && mytrk1->GetNHits()>99   && mytrk2->GetNHits()>99   ) pair_counter[charge_bin][10]++;
                    if( mytrk1->GetTOFDPHI()>99     && mytrk2->GetTOFDPHI()>99   && mytrk1->GetNHits()>999  && mytrk2->GetNHits()>999  ) pair_counter[charge_bin][11]++;
                    if( mytrk1->GetTOFDPHI()>99     && mytrk2->GetTOFDPHI()>99   && mytrk1->GetNHits()>9999 && mytrk2->GetNHits()>9999 ) pair_counter[charge_bin][12]++;
                    if( mytrk1->GetTOFDPHI()>999    && mytrk2->GetTOFDPHI()>999  && mytrk1->GetNHits()>9    && mytrk2->GetNHits()>9    ) pair_counter[charge_bin][13]++;
                    if( mytrk1->GetTOFDPHI()>999    && mytrk2->GetTOFDPHI()>999  && mytrk1->GetNHits()>99   && mytrk2->GetNHits()>99   ) pair_counter[charge_bin][14]++;
                    if( mytrk1->GetTOFDPHI()>999    && mytrk2->GetTOFDPHI()>999  && mytrk1->GetNHits()>999  && mytrk2->GetNHits()>999  ) pair_counter[charge_bin][15]++;
                    if( mytrk1->GetTOFDPHI()>999    && mytrk2->GetTOFDPHI()>999  && mytrk1->GetNHits()>9999 && mytrk2->GetNHits()>9999 ) pair_counter[charge_bin][16]++;
                    if( mytrk1->GetTOFDPHI()>9999   && mytrk2->GetTOFDPHI()>9999 && mytrk1->GetNHits()>9    && mytrk2->GetNHits()>9    ) pair_counter[charge_bin][17]++;
                    if( mytrk1->GetTOFDPHI()>9999   && mytrk2->GetTOFDPHI()>9999 && mytrk1->GetNHits()>99   && mytrk2->GetNHits()>99   ) pair_counter[charge_bin][18]++;
                    if( mytrk1->GetTOFDPHI()>9999   && mytrk2->GetTOFDPHI()>9999 && mytrk1->GetNHits()>999  && mytrk2->GetNHits()>999  ) pair_counter[charge_bin][19]++;
                    if( mytrk1->GetTOFDPHI()>9999   && mytrk2->GetTOFDPHI()>9999 && mytrk1->GetNHits()>9999 && mytrk2->GetNHits()>9999 ) pair_counter[charge_bin][20]++;
                    if( mytrk1->GetTOFDPHI()>9      && mytrk2->GetTOFDPHI()>9      ) pair_counter[charge_bin][21]++;
                    if( mytrk1->GetTOFDPHI()>99     && mytrk2->GetTOFDPHI()>99     ) pair_counter[charge_bin][22]++;
                    if( mytrk1->GetTOFDPHI()>999    && mytrk2->GetTOFDPHI()>999    ) pair_counter[charge_bin][23]++;
                    if( mytrk1->GetTOFDPHI()>9999   && mytrk2->GetTOFDPHI()>9999   ) pair_counter[charge_bin][24]++;
                    if( mytrk1->GetNHits()>9        && mytrk2->GetNHits()>9        ) pair_counter[charge_bin][25]++;
                    if( mytrk1->GetNHits()>99       && mytrk2->GetNHits()>99       ) pair_counter[charge_bin][26]++;
                    if( mytrk1->GetNHits()>999      && mytrk2->GetNHits()>999      ) pair_counter[charge_bin][27]++;
                    if( mytrk1->GetNHits()>9999     && mytrk2->GetNHits()>9999     ) pair_counter[charge_bin][28]++;
                    if( mytrk1->GetHitCounter(0) == 0 && mytrk1->GetHitCounter(1) == 0  && mytrk2->GetTOFDPHI()>999 ) pair_counter[charge_bin][29]++;
                    if( mytrk2->GetHitCounter(0) == 0 && mytrk2->GetHitCounter(1) == 0  && mytrk1->GetTOFDPHI()>999 ) pair_counter[charge_bin][30]++;
                    if( mytrk1->GetHitCounter(0) == 0 && mytrk1->GetHitCounter(1) == 0 && mytrk2->GetHitCounter(0) == 0 && mytrk2->GetHitCounter(1) == 0 ) 
                    {
                        pair_counter[charge_bin][32]++;
                        if( (mytrk1->GetHitCounter(2) > 0  || mytrk1->GetHitCounter(3) > 0) &&  (mytrk2->GetHitCounter(2) > 0  || mytrk2->GetHitCounter(3) > 0)  ) pair_counter[charge_bin][31]++;
                    }
                }
                //if(isGhost) return 1;
            }
        }
        if(event->GetNtrack()>1)
        {
            for (int i = 0; i < 33; i++)
            {
                counter_assoc_ghost_hist->Fill(pair_counter[0][i], i, event->GetCentrality());
                counter_assoc_ghost_hist->Fill(pair_counter[1][i], i, event->GetCentrality()+100);
                counter_assoc_ghost_hist->Fill(pair_counter[2][i], i, event->GetCentrality()+200);
                if (n_ghosts)
                {
                    counter_assoc_ghost_hist->Fill(pair_counter[0][i], i, event->GetCentrality()+300);
                    counter_assoc_ghost_hist->Fill(pair_counter[1][i], i, event->GetCentrality()+400);
                    counter_assoc_ghost_hist->Fill(pair_counter[2][i], i, event->GetCentrality()+500);
                }
            }
        }
        return n_ghosts;
    }

    void MyEventContainer::ResetRecoverFGVars()
    {   
        for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk1 = event->GetEntry(itrk);
            if(mytrk1->GetChargePrime() != -1 ) continue;
            //if(mytrk1->GetMcId() < 100 ) continue;
            
            for (int jtrk = 0; jtrk < event->GetNtrack(); jtrk++)
            {
                MyDileptonAnalysis::MyElectron *mytrk2 = event->GetEntry(jtrk);
                if(mytrk2->GetChargePrime() != 1 ) continue;
                //if(mytrk2->GetMcId() < 100 ) continue;

                const float min_hadorn_reject = mytrk1->GetMcId() < mytrk2->GetMcId() ? mytrk1->GetMcId() : mytrk2->GetMcId();
                //std::cout<<"before: "<<mytrk1->GetNHits()<<" "<<mytrk2->GetNHits()<<" "<<mytrk1->GetTOFDPHI()<<" "<<mytrk2->GetTOFDPHI()<<std::endl;
                //std::cout<<"before: "<<(int)mytrk1->GetPC3SDZ()<<" "<<(int)mytrk2->GetPC3SDZ()<<" "<<mytrk1->GetPC3SDPHI()<<" "<<mytrk2->GetPC3SDPHI()<<std::endl;
                const float min_nhits_reco = mytrk1->GetPC3SDZ() < mytrk2->GetPC3SDZ() ? mytrk1->GetPC3SDZ() : mytrk2->GetPC3SDZ();
                if(min_nhits_reco > mytrk1->GetNHits()) { mytrk1->SetNHits(min_nhits_reco); mytrk1->SetMcId(min_hadorn_reject);}
                if(min_nhits_reco > mytrk2->GetNHits()) { mytrk2->SetNHits(min_nhits_reco); mytrk2->SetMcId(min_hadorn_reject);}


                const float min_ghost_reco = mytrk1->GetPC3SDPHI() < mytrk2->GetPC3SDPHI() ? mytrk1->GetPC3SDPHI() : mytrk2->GetPC3SDPHI();
                if(min_ghost_reco > mytrk1->GetTOFDPHI()) { mytrk1->SetTOFDPHI(min_ghost_reco); mytrk1->SetMcId(min_hadorn_reject);}
                if(min_ghost_reco > mytrk2->GetTOFDPHI()) { mytrk2->SetTOFDPHI(min_ghost_reco); mytrk2->SetMcId(min_hadorn_reject);}
                
                //std::cout<<"after:  "<<mytrk1->GetNHits()<<" "<<mytrk2->GetNHits()<<" "<<mytrk1->GetTOFDPHI()<<" "<<mytrk2->GetTOFDPHI()<<std::endl;
            }
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
                                                     const int fill_DCA, const int fill_track_QA, const int fill_flow, const int fill_true_DCA, 
                                                     const int check_veto, const int fill_inv_mas, const int fill_vertex_reco, const int do_conv_dalitz_finder)
    {
        outfilename = "my-" + outfilename;
        const int compress = 9;
        if (fill_ell || fill_had || fill_tree || fill_dphi || fill_DCA || fill_track_QA || fill_flow || fill_true_DCA || check_veto || fill_inv_mas || fill_vertex_reco)
            outfile = new TFile(outfilename.c_str(), "RECREATE", outfilename.c_str(), compress);

        if (fill_ell)
        {   
            INIT_HISTOS(3, dphi_hist_el_dynamic,  N_dynamic, 100, -0.1, 0.1, 100, -0.1, 0.1, 50, 0, 5);
            INIT_HISTOS(3, dthe_hist_el_dynamic,  N_dynamic, 100, -0.1, 0.1, 100, -0.1, 0.1, 50, 0, 5);
            INIT_HISTOS(3, sdphi_hist_el_dynamic, N_dynamic, 100,  -10,  10, 100,  -10,  10, 50, 0, 5);
            INIT_HISTOS(3, sdthe_hist_el_dynamic, N_dynamic, 100,  -10,  10, 100,  -10,  10, 50, 0, 5);
            //INIT_HISTOS(3, dphi_hist_el_dynamic,  N_dynamic, 100, -0.025, 0.025, 60, -1.5, 4.5, 40, -20, 20);
            //INIT_HISTOS(3, dthe_hist_el_dynamic,  N_dynamic, 100, -0.025, 0.025, 60, -1.5, 4.5, 40, -20, 20);
            //INIT_HISTOS(3, sdphi_hist_el_dynamic, N_dynamic, 100, -0.025, 0.025, 60, -1.5, 4.5, 40, -20, 20);
            //INIT_HISTOS(3, sdthe_hist_el_dynamic, N_dynamic, 100, -0.025, 0.025, 60, -1.5, 4.5, 40, -20, 20);
            INIT_HISTOS(3, chi2_ndf, N_centr,      50, 0, 10,  20, 0, 20, 25, 0, 5);
            INIT_HISTOS(3, ilayerhitshist, N_centr,50, -0.5, 49.5, 40, 0, 40, 50, 0, 5);
            INIT_HISTOS(3, dphi_hist_el,  1, 50, -0.1, 0.1, 8, 0, 8, 5, 0, 5);
            INIT_HISTOS(3, dthe_hist_el,  1, 50, -0.1, 0.1, 8, 0, 8, 5, 0, 5);
            INIT_HISTOS(3, sdphi_hist_el, 1, 50, -10, 10,   8, 0, 8, 5, 0, 5);
            INIT_HISTOS(3, sdthe_hist_el, 1, 50, -10, 10,   8, 0, 8, 5, 0, 5);
            INIT_HIST  (3, truehithist,      10, 0, 10, 50, 0, 5, 10, 0, 100);
            INIT_HIST  (3, truehitsigmahist, 50, 0, 50, 50, 0, 5, 10, 0, 100);
            INIT_HIST  (3, charge_recover_hist, 2, -2, 2, 16, 0, 16, 50, 0, 5);
            if(fill_ell==2)
            {
                INIT_HISTOS(3, dphi_hist_el,  N_centr, 100, -0.1, 0.1, 8, 0, 8, 50, 0, 5);
                INIT_HISTOS(3, dthe_hist_el,  N_centr, 100, -0.1, 0.1, 8, 0, 8, 50, 0, 5);
                INIT_HISTOS(3, sdphi_hist_el, N_centr, 100, -10, 10,   8, 0, 8, 50, 0, 5);
                INIT_HISTOS(3, sdthe_hist_el, N_centr, 100, -10, 10,   8, 0, 8, 50, 0, 5);
                is_fill_hsits = 2;
            }
            if(fill_ell<3)is_fill_hsits = 1;
        }

        if (fill_had)
        {
            INIT_HISTOS(3, dphi_hist,  N_centr, 100, -0.1, 0.1, 16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, dthe_hist,  N_centr, 100, -0.1, 0.1, 16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, sdphi_hist, N_centr, 100, -10, 10,   16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, sdthe_hist, N_centr, 100, -10, 10,   16, 0, 16, 50, 0, 5);
            INIT_HISTOS(3, dphi_phi0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 4, 0, 4);
            INIT_HISTOS(3, dthe_the0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 4, 0, 4);
            INIT_HISTOS(3, dphi_phi0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 4, 0, 4);
            INIT_HISTOS(3, dthe_the0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 4, 0, 4);
            INIT_HISTOS(3, dthe_phi0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 4, 0, 4);
            INIT_HISTOS(3, dphi_the0_init_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 4, 0, 4);
            INIT_HISTOS(3, dthe_phi0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, -1.57, 4.71, 4, 0, 4);
            INIT_HISTOS(3, dphi_the0_corr_hist,  nvtx_layers, 400, -0.05, 0.05, 120, 0.785, 2.36, 4, 0, 4);
            //INIT_HISTOS(3, myvtx_hist, N_centr, 1000, -1, 1, 12, 0 ,12, 4, 0 ,4);
            is_fill_hadron_hsits = 1;
        }
        if (fill_tree)
        {
            tree = new TTree("tree", "tree");
            tree->Branch("MyEvent", event);
            INIT_HIST(1, event_hist, 10, 0, 10);
            INIT_HIST(2, centr_hist, 100, 0, 100, 2000, 0, 2000);
            INIT_HISTOS(3, el_pt_hist, N_centr*2, 50, 0, 5, 2, 0, 2, 100, 0, 100);
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
            INIT_HIST(3, ep_hist, 50, 0, 1.5, 50, 0, 5.0, 30, 0, 600);
            INIT_HIST(3, n0_hist, 10, 0, 10, 50, 0, 5.0, 30, 0, 600);
            INIT_HIST(3, prob_hist, 100, 0, 1, 50, 0, 5.0, 30, 0, 600);
            INIT_HIST(3, disp_hist, 50, 0, 10, 50, 0, 5.0, 30, 0, 600);
            INIT_HIST(3, chi2npe0_hist, 100, 0, 20, 50, 0., 5.0, 30, 0, 600);

            INIT_HIST(3, el_had_dphi, 100, -0.1, 0.1, 50, 0.0, 5.0, 30, 0, 600);
            INIT_HIST(3, el_had_dz, 100, -50, 50, 50, 0.0, 5.0, 30, 0, 600);
            INIT_HIST(3, el_had_dr, 100, 0, 20, 50, 0., 5.0, 30, 0, 600);
            INIT_HIST(3, temc, 1000, -0.5, 1.5, 200, 0., 10, 30, 0, 600);
            INIT_HIST(3, ttof, 1000, -0.5, 1.5, 50, 0., 5, 30, 0, 600);

            INIT_HISTOS(3, el_pt_hist, N_centr*2, 100, 0, 10, 5, 0, 5, 25, 0, 25);
            INIT_HIST(3, BDT_eID_hist, 1000, 0, 1, 50, 0, 5.0, 40, 0, 400);

            INIT_HIST(3, ep_hist_el, 30, 0, 1.5,  100, 0, 1, 50, 0., 5.0);
            INIT_HIST(3, n0_hist_el, 10, 0, 10, 50, 0, 10, 10, 0., 100);
            INIT_HIST(3, prob_hist_el, 20, 0, 20, 50, 0, 10, 10, 0., 100);
            INIT_HIST(3, disp_hist_el, 50, 0, 10, 30, 0, 30, 10, 0., 100);
            INIT_HIST(3, chi2npe0_hist_el, 50, 0, 20, 30, 0, 30, 10, 0., 100);
            INIT_HIST(3, rich_prob1, 50, 0, 20, 50, -10, 10, 10, 0., 100);
            INIT_HIST(3, rich_prob2, 30, 0, 30, 50, -10, 10, 10, 0., 100);
            INIT_HIST(3, rich_prob3, 100, -0.05, 0.05, 100, -25, 25, 10, 0., 100);
            is_fill_track_QA = 1;
        }
        if(fill_flow)
        {
            INIT_HIST(3, BBC_psi_hist ,  64, -1.6, 1.6, 20, 0, 100, 3, 0, 3);
            INIT_HIST(3, FVTX_psi_hist,  64, -1.6, 1.6, 20, 0, 100, 3, 0, 3);
            INIT_HIST(3, cos_BBC_hist , 200, -1.0, 1.0, 20, 0, 100, 3, 0, 3);
            INIT_HIST(3, cos_FVTX_hist, 200, -1.0, 1.0, 20, 0, 100, 3, 0, 3);
            INIT_HIST(3, v2_BBC_hist ,  200, -1.0, 1.0, 50, 0, 5, 40, 0, 200);
            INIT_HIST(3, v2_FVTX_hist,  200, -1.0, 1.0, 50, 0, 5, 40, 0, 200);
        }
        if (fill_true_DCA)
        {
            INIT_HIST(1, event_hist, 20, 0, 20);
            INIT_HIST(2, centr_hist, 100, 0, 100, 2000, 0, 2000);
            INIT_HIST(3, counter_assoc_ghost_hist, 50,0,50, 33, 0, 33, 60, 0, 600);
            INIT_HISTOS(3, myvtx_hist, N_centr, 100, 0, 1, 100, -0.5 ,0.5, 20, -10 ,10);

            INIT_HIST(3, DCPT_ReconPT, 50, 0, 5, 50,  0,  5, 10, 0, 200);

            INIT_HISTOS(3, DCA12_hist, N_centr, 50, -500, 500, 100, -5, 5, 160, 0, 160);
            INIT_HISTOS(3, DCA2_hist, N_centr, 200, -4000, 4000, 50, 0, 5, 25, 0, 25);
            INIT_HISTOS(3, sDCA2_hist, N_centr, 200, -4000, 4000, 50, 0, 5, 25, 0, 25);
            INIT_HISTOS(3,  DCA_2D_hist, N_centr, 200, -4000, 4000, 50, 0, 5, 25, 0, 25);
            INIT_HISTOS(3, sDCA_2D_hist, N_centr, 200, -4000, 4000, 50, 0, 5, 25, 0, 25);

            INIT_HIST(3, vtx_accaptance_hist, 300, -1.5, 4.5, 240 , -12, 12, 8, 0 ,8 );
            INIT_HIST(3, vtx_deadmaps_hist,   30, -1.5, 4.5, 24 , -12, 12, 8, 0 ,8 );

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
            //INIT_HIST(3, temc, 150, -50, 150, 50, 0., 2.5, 5, 0, 5);
            INIT_HIST(3, ttof, 1500, -50, 150, 50, 0., 1.25, 5, 0, 5);
            is_check_veto = 1;
        }
        if(fill_inv_mas)
        {
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
        if(fill_vertex_reco)
        {
            do_vertex_reco = 1;

            INIT_HIST( 3, hist_dca_x, 200,  -0.05, 0.05, 50,   0.0, 5.0, 320, 0, 320);
            INIT_HIST( 3, hist_dca_y, 200,  -0.05, 0.05, 50,   0.0, 5.0, 320, 0, 320);
            INIT_HIST( 3, hist_vtx_x, 200,  -0.05, 0.05, 100,  0.1, 0.6, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_y, 200,  -0.05, 0.05, 100,  0.0, 0.2, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_z, 100,   0.00, 0.60, 100,  0.0, 0.2, 20, -10, 10);
            INIT_HIST( 3, hist_vtx_grid_xy, 100,   0.28, 0.38, 100,  -0.01, 0.09, 10, 0, 100);
            INIT_HIST( 3, vtx_dphi_dphi_hist,     100, -0.05, 0.05, 100, -0.05, 0.05, 50, 0, 5);
            INIT_HIST( 3, vtx_dthe_dthe_hist,     100, -0.05, 0.05, 100, -0.05, 0.05, 50, 0, 5);
            INIT_HIST( 3, vtx_dca_pion_hist, 200, -0.10, 0.10, 200, 0, 0.02, 250, 0, 25);
            INIT_HIST( 3, phi_the_pt_hist, 100, -3.14/2, 3.14*3/2, 100, 0.6, 2.6, 50, -5, 5);
            INIT_HIST( 3, conv_photon_mass_hist, 200, -0.05, 0.05, 60, -1.5, 4.5, 100, -5.0, 5.0);
            INIT_HIST( 3, pi0_mass_hist,         200, -0.05, 0.05, 60, -1.5, 4.5, 100, -5.0, 5.0);
            INIT_HIST( 2, hits_vtx_ntracks, 1000, 0, 1000, 10, 0, 100);
            INIT_HIST( 2, hits_vtx_ntracks_ofnotusedhits, 1000, 0, 1000, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_delta_x, 200, -0.05, 0.05, 100, 0.1, 0.6, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_delta_y, 200, -0.05, 0.05, 100, 0.0, 0.2, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_delta_x_reuse, 200, -0.05, 0.05, 100, 0.1, 0.6, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_delta_y_reuse, 200, -0.05, 0.05, 100, 0.0, 0.2, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_delta_x_pion, 200, -0.05, 0.05, 100, 0.1, 0.6, 10, 0, 100);
            INIT_HIST( 3, hist_vtx_delta_y_pion, 200, -0.05, 0.05, 100, 0.0, 0.2, 10, 0, 100);
        }
        if(do_conv_dalitz_finder==2)
        {
            
            INIT_HISTOS(3, hist_conv_phi_phi,    6,       100, -0.1, 0.1, 20, -0.1, 0.1, 25, 0, 5); 
            INIT_HISTOS(3, hist_conv_the_the,    6,       100, -0.1, 0.1, 20, -0.1, 0.1, 25, 0, 5); 
            INIT_HISTOS(3, hist_daltz_phi_phi,   12,      100, -0.1, 0.1, 20, -0.1, 0.1, 25, 0, 5); 
            INIT_HISTOS(3, hist_daltz_the_the,   12,      100, -0.1, 0.1, 20, -0.1, 0.1, 25, 0, 5); 
            INIT_HISTOS(3, hist_is_dalitz_conv,  N_centr, 20, 0, 20, 20, 0, 20, 50, 0, 5); 
            INIT_HISTOS(3, hist_is_ml_conv,      N_centr, 20, 0, 20, 10, 0, 10, 50, 0, 5); 
            INIT_HISTOS(3, sdphi_conv_hist,      N_centr,  50,  -10,  10, 50, -5.0, 5.0, 50, 0, 5); 
            INIT_HISTOS(3, sdphi_real_conv_hist, N_centr,  50,  -10,  10, 50, -5.0, 5.0, 50, 0, 5);
        }
    }
    
    void MyEventContainer::WriteOutFile()
    {
        std::cout << "Start writing hists to My outfile" << std::endl;
        infile->Close();
        if (is_fill_tree || is_fill_hadron_hsits || is_fill_hsits || is_fill_dphi_hist || is_fill_DCA_hist || is_fill_track_QA
        || is_fill_flow || is_fill_DCA2_hist||is_check_veto||is_fill_inv_mass||do_vertex_reco)
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

    void MyEventContainer::Associate_Hits_to_Leptons_OLD(float sigma, float sigma_veto, float sigma_inner, int not_fill)
    {
        const int nleptons = event->GetNtrack();
        const int nvtxhits = event->GetNVTXhit();
        const int centrality = event->GetCentrality();
        const int rungroup = event->GetRunGroup();

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

            int is_fill_hsits_local = is_fill_hsits * (1-not_fill);  /////to be removed
            if (not_fill<0)
            {
                if( mytrk->GetCrkphi()<-99)
                    is_fill_hsits_local = 1;
                else is_fill_hsits_local = 0;
            }
            
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
            std::vector<std::vector<double> > hit_counter_jlayer[total_vtx_layers+2];

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
                            sigma_phi_value   = mytrk->get_dynamic_sigma_phi_data  (cycle_layer, dphi_previous_layer);
                            mean_phi_value    = mytrk->get_dynamic_mean_phi_data   (cycle_layer, dphi_previous_layer);
                            sigma_theta_value = mytrk->get_dynamic_sigma_theta_data(cycle_layer, dthe_previous_layer);
                            mean_theta_value  = mytrk->get_dynamic_mean_theta_data (cycle_layer, dthe_previous_layer);
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
                            std::vector<double> var_counter_jlayer;  
                            var_counter_jlayer.push_back(sdphi);var_counter_jlayer.push_back(sdthe);var_counter_jlayer.push_back(phi_hit);var_counter_jlayer.push_back(theta_hit);
                            hit_counter_jlayer[ilayer].push_back(var_counter_jlayer);
                            if(layer==2) hit_counter_jlayer[8].push_back(var_counter_jlayer);
                            if(layer==3) hit_counter_jlayer[9].push_back(var_counter_jlayer);
                            //std::cout<<ilayer<<" "<<radii[ilayer]<<" "<<sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))<<" "<<dilep_phi_projection[ilayer]<<" "<<mytrk->GetPtPrime()<<" "<<
                            //(dilep_phi_projection[ilayer]-dilep_phi_projection[ilayer<7?ilayer+1:6])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer]-radii[ilayer<7?ilayer+1:6])*1000<<" "<<
                            //(dilep_phi_projection[ilayer]-dilep_phi_projection[ilayer<6?ilayer+2:5])/(radii[ilayer]-radii[ilayer<6?ilayer+2:5])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])*1000<<std::endl;
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
                if(is_fill_hsits_local) 
                {
                    float hit_counter_jlayer_inside[total_vtx_layers*5] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                    int passed_association_in_sigma[3][nvtx_layers] = {{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0}};
                    for (int jlayer = 0; jlayer < total_vtx_layers+2; jlayer++)
                    {
                        for (unsigned int jentry = 0; jentry < hit_counter_jlayer[jlayer].size(); jentry++)
                        {
                            for (unsigned int kentry = jentry+1; kentry < hit_counter_jlayer[jlayer].size(); kentry++)
                            {
                                if ( TMath::Abs( hit_counter_jlayer[jlayer][jentry][2] - hit_counter_jlayer[jlayer][kentry][2] ) < 0.001
                                  && TMath::Abs( hit_counter_jlayer[jlayer][jentry][3] - hit_counter_jlayer[jlayer][kentry][3] ) < 0.001 )
                                    hit_counter_jlayer[jlayer][kentry][0] = 100;
                            } 
                        }                        
                    }
                    for (int jlayer = 0; jlayer < nvtx_layers; jlayer++)
                    {
                        int iarraylayer = jlayer;
                        if(jlayer>1) iarraylayer = jlayer+6;
                        for (unsigned int jentry = 0; jentry < hit_counter_jlayer[iarraylayer].size(); jentry++)
                        {
                            if ( TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][0])<4 && TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][1])<4 ) 
                                passed_association_in_sigma[0][jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][0])<3 && TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][1])<3 ) 
                                passed_association_in_sigma[1][jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][0])<2 && TMath::Abs(hit_counter_jlayer[iarraylayer][jentry][1])<2 ) 
                                passed_association_in_sigma[2][jlayer]++;
                        }                        
                    }
                    for (int jlayer = 0; jlayer < total_vtx_layers+2; jlayer++)
                    {
                        for (unsigned int jentry = 0; jentry < hit_counter_jlayer[jlayer].size(); jentry++)
                        {
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<5 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<5) 
                                hit_counter_jlayer_inside[jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<4 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<4 
                            && passed_association_in_sigma[0][0] && passed_association_in_sigma[0][1] && 
                            (passed_association_in_sigma[0][2] || passed_association_in_sigma[0][3] ) ) hit_counter_jlayer_inside[10+jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<3 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<3 
                            && passed_association_in_sigma[1][0] && passed_association_in_sigma[1][1] && 
                            (passed_association_in_sigma[1][2] || passed_association_in_sigma[1][3] ) ) hit_counter_jlayer_inside[20+jlayer]++;
                            if ( TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<2 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<2 
                            && passed_association_in_sigma[2][0] && passed_association_in_sigma[2][1] && 
                            (passed_association_in_sigma[2][2] || passed_association_in_sigma[2][3] ) )  hit_counter_jlayer_inside[30+jlayer]++;
                            //if ( jlayer == 0 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][0])<2 && TMath::Abs(hit_counter_jlayer[jlayer][jentry][1])<2 
                            //&& passed_association_in_sigma[2][0]>1 && passed_association_in_sigma[2][1] && 
                            //(passed_association_in_sigma[2][2] || passed_association_in_sigma[2][3] ) ) 
                            //    std::cout<<mytrk->GetPtPrime()<<" "<<hit_counter_jlayer[jlayer][jentry][0]<<" "<<hit_counter_jlayer[jlayer][jentry][1]<<" "
                            //    <<hit_counter_jlayer[jlayer][jentry][2]<<" "<<hit_counter_jlayer[jlayer][jentry][3]<<std::endl;

                        }                        
                    }
                    for (int jlayer = 0; jlayer < total_vtx_layers*5; jlayer++)
                    {
                        ilayerhitshist[central_bin]->Fill(hit_counter_jlayer_inside[jlayer], jlayer, pt);
                    }
                }
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
    
    void MyEventContainer::Associate_Hits_to_Hadrons_Dynamic(float sigma = 5, float vertex_x = -999, float vertex_y = -999)
    {
        const int nhadrons = event->GetNeleccand();
        const int nvtxhits = event->GetNVTXhit();
        const int centrality = event->GetCentrality();
        const int central_bin = (int)centrality / 20;
        const int rungroup = event->GetRunGroup();
        const int is_fill_hsits_hdphi = vertex_x>-99?1:0;
        if (central_bin > 4 || central_bin < 0)
            return;

        const float init_vtx_x = vertex_x>-99? vertex_x:event->GetPreciseX();
        const float init_vtx_y = vertex_y>-99? vertex_y:event->GetPreciseY();

        if(vertex_x>-99 && vertex_y>-99)
        {
            if(nhadrons==1){
                const float pt = event->GetElecCand(0)->GetPtPrime() + central_bin*5;
                const float dca = event->GetElecCand(0)->GetDCA2()/10000;
                const float delta_r_vtx = sqrt(SQR(init_vtx_x - event->GetPreciseX()) + SQR(init_vtx_y - event->GetPreciseY()));
                vtx_dca_pion_hist->Fill(dca, delta_r_vtx, pt);
            }
        }
        
        for (int itrk = 0; itrk < nhadrons; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetElecCand(itrk);

            mytrk->ZeroHitCounters();
            mytrk->ClearNumberVectors();
            const float pt = mytrk->GetPtPrime();

            const float thetaprime = mytrk->GetThe0Prime();

            float phi0_trk_proj = mytrk->GetPhi0Prime();
            float the0_trk_proj = mytrk->GetThe0Prime();
            const float pz = mytrk->GetPtPrime() * (TMath::Cos(thetaprime)) / (TMath::Sin(thetaprime));

            float rp = sqrt(init_vtx_x * init_vtx_x + init_vtx_y * init_vtx_y);
            float xp = init_vtx_x;
            float yp = init_vtx_y;
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
            std::vector<std::vector<double> > hit_counter_jlayer[total_vtx_layers+2];

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

                        const float phi_hit = vtxhit->GetPhiHit(init_vtx_x,init_vtx_y,event->GetPreciseZ());
                        const float theta_hit = vtxhit->GetTheHit(init_vtx_x,init_vtx_y,event->GetPreciseZ());
                        
                        const float dphi = (dilep_phi_projection[ilayer] - phi_hit) + 
                            (dilep_phi_projection[ilayer<7?ilayer+1:6] - dilep_phi_projection[ilayer])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer<7?ilayer+1:6] - radii[ilayer]);
                        const float dthe = (dilep_the_projection[ilayer] - theta_hit) + 
                            (dilep_the_projection[ilayer<7?ilayer+1:6] - dilep_the_projection[ilayer])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer<7?ilayer+1:6] - radii[ilayer]);
                        if (abs(dphi) > 0.1 || abs(dthe) > 0.1) continue;

                        float sigma_phi_value = mytrk->get_sigma_phi_data(0*rungroup, central_bin, layer);
                        float mean_phi_value = mytrk->get_mean_phi_data(0*rungroup, central_bin, layer);
                        float sigma_theta_value = mytrk->get_sigma_theta_data(0*rungroup, central_bin, layer);
                        float mean_theta_value = mytrk->get_mean_theta_data(0*rungroup, central_bin, layer);

                        int cycle_layer = layer;
                        if((iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) || iter_layer==2) cycle_layer++;
                        if(iter_layer<2||iassociatedhit>0)
                        {
                            sigma_phi_value   = mytrk->get_dynamic_sigma_phi_data  (cycle_layer, dphi_previous_layer);
                            mean_phi_value    = mytrk->get_dynamic_mean_phi_data   (cycle_layer, dphi_previous_layer);
                            sigma_theta_value = mytrk->get_dynamic_sigma_theta_data(cycle_layer, dthe_previous_layer);
                            mean_theta_value  = mytrk->get_dynamic_mean_theta_data (cycle_layer, dthe_previous_layer);
                        }

                        const float sdphi = (dphi - mean_phi_value) / sigma_phi_value;// - mytrk->get_dynamic_smean_phi_data(0, cycle_layer, dphi_previous_layer);
                        const float sdthe = (dthe - mean_theta_value) / sigma_theta_value;

                        const float diff = sqrt(std::pow(sdphi, 2) + std::pow(sdthe, 2));


                        bool SignTrack = true;
                        if ( sdphi*mytrk->GetChargePrime()>-sigma && sdphi*mytrk->GetChargePrime() < sigma && TMath::Abs(sdthe) < 2)
                        {
                            vtxhit->SetLadder(44);
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
                                mytrk->SetHitCounter(0,1);
                                return;
                            }
                            if(iter_layer==2 && iassociatedhit >0) numbers[2].push_back(iter_nums[layer]*100  +numbers[3][iassociatedhit-1]);
                            if(iter_layer==2 && iassociatedhit==0) numbers[2].push_back(iter_nums[layer]*100  );
                            if(iter_layer==1 && iassociatedhit <  mytrk->GetHitCounter(2)) numbers[1].push_back(iter_nums[layer]*10000 +numbers[2][iassociatedhit]);
                            if(iter_layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) numbers[1].push_back(iter_nums[layer]*10000 +numbers[3][iassociatedhit-mytrk->GetHitCounter(2)]);
                            if(iter_layer==0 ) numbers[0].push_back(iter_nums[layer]*10000000+numbers[1][iassociatedhit]);
                            std::vector<double> var_counter_jlayer;  
                            var_counter_jlayer.push_back(sdphi);var_counter_jlayer.push_back(sdthe);var_counter_jlayer.push_back(phi_hit);var_counter_jlayer.push_back(theta_hit);
                            hit_counter_jlayer[ilayer].push_back(var_counter_jlayer);
                            if(layer==2) hit_counter_jlayer[8].push_back(var_counter_jlayer);
                            if(layer==3) hit_counter_jlayer[9].push_back(var_counter_jlayer);
                            //std::cout<<ilayer<<" "<<radii[ilayer]<<" "<<sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))<<" "<<dilep_phi_projection[ilayer]<<" "<<mytrk->GetPtPrime()<<" "<<
                            //(dilep_phi_projection[ilayer]-dilep_phi_projection[ilayer<7?ilayer+1:6])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])/(radii[ilayer]-radii[ilayer<7?ilayer+1:6])*1000<<" "<<
                            //(dilep_phi_projection[ilayer]-dilep_phi_projection[ilayer<6?ilayer+2:5])/(radii[ilayer]-radii[ilayer<6?ilayer+2:5])*(sqrt(SQR(vtxhit->GetXHit())+SQR(vtxhit->GetYHit()))-radii[ilayer])*1000<<std::endl;
                        } // end of association
                        else
                        {
                            if (vtxhit->N_AssociatedTracks() > 0)
                                SignTrack = false;
                        }
                        int in_arg = 1*mytrk->GetArm()+4*layer+2*charge_bin;
                        if( (layer==1 && iassociatedhit >= mytrk->GetHitCounter(2)) || (layer==2 && iassociatedhit>0) ) in_arg+=4;
                        if(iter_layer>1 && iassociatedhit==0) in_arg+=8;

                        if (TMath::Abs(sdthe) < 2 && SignTrack && is_fill_hsits_hdphi )
                        {
                            //dphi_hist_el_dynamic[in_arg]->Fill(dphi-dphi_previous_layer, vtxhit->GetPhiHit(0,0,event->GetPreciseZ()), vtxhit->GetZHit());
                            //if(pt>1.4)sdphi_hist_el_dynamic[in_arg]->Fill(dphi-dphi_previous_layer, vtxhit->GetPhiHit(0,0,event->GetPreciseZ()), vtxhit->GetZHit());
                            dphi_hist_el_dynamic[in_arg]->Fill(dphi, dphi_previous_layer, pt);
                            sdphi_hist_el_dynamic[in_arg]->Fill(sdphi, sdphi_previous_layer, pt);
                        }
                        if (sdphi*mytrk->GetChargePrime()>-2 && sdphi*mytrk->GetChargePrime() < 2 && SignTrack && is_fill_hsits_hdphi)//TMath::Abs(sdthe) < sigma && 
                        {
                            //dthe_hist_el_dynamic[in_arg]->Fill(dthe-dthe_previous_layer, vtxhit->GetTheHit(0,0,event->GetPreciseZ()), vtxhit->GetZHit());
                            //if(pt>1.4)sdthe_hist_el_dynamic[in_arg]->Fill(dthe-dthe_previous_layer, vtxhit->GetTheHit(0,0,event->GetPreciseZ()), vtxhit->GetZHit());
                            dthe_hist_el_dynamic[in_arg]->Fill(dthe, dthe_previous_layer, pt);
                            sdthe_hist_el_dynamic[in_arg]->Fill(sdthe, sdthe_previous_layer, pt);
                        }
                    } // enf of hit loop
                } // end of hits in prev layer 
            } //ens of layers   
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
                    float chi2 = 800001;
                    float recon_pt = 0;
                    mytrk->SetHitIndex(mytrk->GetHits(0,inum0), 0);
                    mytrk->SetHitIndex(mytrk->GetHits(1,inum1), 1);
                    if (inum2>=0) mytrk->SetHitIndex(mytrk->GetHits(2,inum2), 2);
                    if (inum3>=0) mytrk->SetHitIndex(mytrk->GetHits(3,inum3), 3);
                    if( inum3>=0) 
                    {
                        event->SetDCA2(itrk,3,1);
                        recon_pt += mytrk->GetReconPT();
                    }
                    if( inum2>=0) 
                    {
                        event->SetDCA2(itrk,2,1);
                        recon_pt += mytrk->GetReconPT();
                        if (inum3>=0) recon_pt/=2;
                    }
                    
                    chi2 = SQR(recon_pt-pt)/pt*30/(2+(int)(inum2>=0)+(int)(inum3>=0));
                    if(chi2<min_chi2) {min_chi2=chi2;final_number=numbers[0][inum];} 
                    
                    if (numbers[0].size()<10 && is_fill_hsits_hdphi) chi2_ndf[central_bin]->Fill(chi2, numbers[0].size(), pt);
                }
            }
            if(is_fill_hsits_hdphi) chi2_ndf[central_bin]->Fill(min_chi2, 19, pt);
            mytrk->SetHitCounter(0,0);
            mytrk->SetHitCounter(1,0);
            if(min_chi2<800000)
            {
                mytrk->SetHitCounter(3,0);mytrk->SetHitCounter(2,0);

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
                //event->SetDCA2(itrk, 1);
                if (mytrk->GetHitCounter(3)>0)  event->SetDCA2(itrk, 3, 1);
                if (mytrk->GetHitCounter(2)>0)  event->SetDCA2(itrk, 2, 1);
               

            }else{
                mytrk->SetHitCounter(0,0);
            }
            //mytrk->ClearNumberVectors();
        }     // enf of e loop
        int i = 0;
        while (i < event->GetNeleccand())
        {
            MyDileptonAnalysis::MyElectron myel = *event->GetElecCand(i);
            if (myel.GetHitCounter(0)==0)
            {
                event->RemoveElecCandEntry(i);
            }
            else{
                i++; 
            }
        }
    }         // end

    void MyEventContainer::FillTrueDCAHadrons(const float weight)
    {
        const int central_bin = (int)event->GetCentrality() / 20;
        const int nhadorns = event->GetNeleccand();
        for (int itrk = 0; itrk < nhadorns; itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetElecCand(itrk);
            const float pt = mytrk->GetPtPrime();
            int hit_assocaition = 0;
            if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<3) ||
                   (TMath::Abs(mytrk->GetMinsDphi(2))<3) ) && 
                   (TMath::Abs(mytrk->GetMinsDphi(1))<3) && 
                   (mytrk->GetMinsDphi(0))>-3 ) ) hit_assocaition=5;
            if (hit_assocaition==5 && mytrk->GetMinsDphi(0)<3 ) hit_assocaition=10;
            if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<2 && TMath::Abs(mytrk->GetMinsDthe(3))<2) ||
                   (TMath::Abs(mytrk->GetMinsDphi(2))<2 && TMath::Abs(mytrk->GetMinsDthe(2))<2) ) && 
                   (TMath::Abs(mytrk->GetMinsDphi(1))<2) && 
                   (mytrk->GetMinsDphi(0)> -2 ) )) hit_assocaition=15;
            if (hit_assocaition==15 && mytrk->GetMinsDphi(0)<2 ) hit_assocaition=20;
            int conv_reject = 0;
            if ( TMath::Abs(mytrk->GetMinsDphi(3))<2 || TMath::Abs(mytrk->GetMinsDphi(2))<2) conv_reject=1;
            if ( TMath::Abs(mytrk->GetMinsDphi(2))<2 && TMath::Abs(mytrk->GetMinsDphi(2))<2) conv_reject=2;
            if ( TMath::Abs(mytrk->GetMinsDphi(1))<2) conv_reject=3;
            if ( TMath::Abs(mytrk->GetMinsDphi(0))<2) conv_reject=4;

            const int hist_in = hit_assocaition + conv_reject;

            if (mytrk->GetChargePrime()>0)
                DCA2_hist[central_bin] ->Fill(mytrk->GetDCA2(),pt,hist_in,weight);
            else
                sDCA2_hist[central_bin]->Fill(mytrk->GetDCA2(),pt,hist_in,weight);
            
            if( TMath::Abs(event->GetPreciseZ()) > 12 ) continue; // skip events with large z vertex
            int z_bin =  ( (int) ((event->GetPreciseZ()+12)/3) )*20 + ( (int) ((mytrk->GetPhi0Prime() + 1.5) / 0.3) ); 
            DCA12_hist[central_bin]->Fill(mytrk->GetDCA2()*mytrk->GetChargePrime(),pt*mytrk->GetChargePrime(),z_bin,weight);

            DCPT_ReconPT->Fill(mytrk->GetReconPT(),pt,event->GetCentrality()+100*( mytrk->GetChargePrime() > 0 ? 0 : 1));
            
        }
        for (int itrk = 0; itrk < event->GetNtrack(); itrk++)
        {
            MyDileptonAnalysis::MyElectron *mytrk = event->GetEntry(itrk);
            const float pt = mytrk->GetPtPrime();
            int hit_assocaition = 0;
            if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<3) ||
                   (TMath::Abs(mytrk->GetMinsDphi(2))<3) ) && 
                   (TMath::Abs(mytrk->GetMinsDphi(1))<3) && 
                   (mytrk->GetMinsDphi(0))>-3 ) ) hit_assocaition=1;
            if (hit_assocaition==1 && mytrk->GetMinsDphi(0)<3 ) hit_assocaition=2;
            if ( (((TMath::Abs(mytrk->GetMinsDphi(3))<2 && TMath::Abs(mytrk->GetMinsDthe(3))<2) ||
                   (TMath::Abs(mytrk->GetMinsDphi(2))<2 && TMath::Abs(mytrk->GetMinsDthe(2))<2) ) && 
                   (TMath::Abs(mytrk->GetMinsDphi(1))<2) && 
                   (mytrk->GetMinsDphi(0)> -1 ) )) hit_assocaition=3;
            if (hit_assocaition==3 && mytrk->GetMinsDphi(0)<2 ) hit_assocaition=4;
            int conv_reject = 0;
            if ( ((int)mytrk->GetEmcdphi_e())%100==0) conv_reject=5;
            if ( ((int)mytrk->GetEmcdphi_e())%100<3 && ((int)mytrk->GetEmcdphi_e())/100<3) conv_reject=10;
            if ( ((int)mytrk->GetEmcdphi_e())%100<1 && ((int)mytrk->GetEmcdphi_e())/100<3) conv_reject=15;
            if ( ((int)mytrk->GetEmcdphi_e())%100<1 && ((int)mytrk->GetEmcdphi_e())/100<1) conv_reject=20;

            const int hist_in = hit_assocaition + conv_reject;

            if (mytrk->GetChargePrime()>0)
                DCA_2D_hist[central_bin] ->Fill(mytrk->GetDCA2(),pt,hist_in,weight);
            else
                sDCA_2D_hist[central_bin]->Fill(mytrk->GetDCA2(),pt,hist_in,weight);
        }
    }

    void MyEventContainer::VertexXYScanDC(const float run_beam_x = 0.328, const float run_beam_y = 0.038, int fill_hist = 0, int verbosity = 0)
    {
        if(event->GetCentrality()>56) return; // skip low centrality events 
        float fraction = 0.75;
        float sigma = 5;     
        float x_range = 0.05;    // +/- cm
        float y_range = 0.05;    // +/- cm
        float step_size = 0.005; // cm = 50 microns
        float prev_x = event->GetPreciseX();
        float prev_y = event->GetPreciseY();
        float beam_x = run_beam_x;
        float beam_y = run_beam_y;
        float beam_z = event->GetPreciseZ();
        float best_x = beam_x;
        float best_y = beam_y;// ...existing code...
        int init_n =event->GetNeleccand();
        int i = 0;
        while (i < event->GetNeleccand())
        {
            MyDileptonAnalysis::MyElectron myel = *event->GetElecCand(i);
            if (myel.GetN0() >= 0 || myel.GetPtPrime()<0.4 || myel.GetEcore()/myel.GetPtot()>0.6 || myel.GetEcore()<0)
            {
                event->RemoveElecCandEntry(i);
            }
            else{
                i++; 
            }
        }
        Associate_Hits_to_Hadrons_Dynamic(10);
        i = 0;
        while (i < event->GetNeleccand())
        {
            MyDileptonAnalysis::MyElectron myel = *event->GetElecCand(i);
            if (!myel.GetHitCounter(0))
            {
                event->RemoveElecCandEntry(i);
            }
            else{
                i++; 
            }
        }
        
        const int Ntracks = event->GetNeleccand();
        if(verbosity)
            std::cout << "\033[35mVertexXYScanDC: Ntracks = " << Ntracks << " of " << init_n << " at vertex "<< prev_x<< " "<< prev_y << " at centrality " << event->GetCentrality() << "\033[0m" << std::endl;
        std::vector<std::pair<float, float> > track_vertices, track_vertices_west, track_vertices_east; // (x at y_beam, y at x_beam)
        std::vector<int> n_tracks_vec, n_tracks_vec_west, n_tracks_vec_east;
        int maximum_tracks = 0;
        for (float xvtx = beam_x - x_range; xvtx <= beam_x + x_range; xvtx += step_size)
        {
            for (float yvtx = beam_y - y_range; yvtx <= beam_y + y_range; yvtx += step_size)
            {
                event->SetPreciseX(xvtx);
                event->SetPreciseY(yvtx);
                Associate_Hits_to_Hadrons_Dynamic(sigma);
                int n_tracks = 0, n_wtracks = 0;
                int n_tracks_east = 0, n_tracks_west = 0;
                for (int itrk = 0; itrk < Ntracks; itrk++)
                {
                    MyDileptonAnalysis::MyElectron *mytrk = event->GetElecCand(itrk);
                    if (!mytrk->GetHitCounter(0)) continue;
                    double dca = mytrk->GetDCA2()/10000.;
                    const double rescattering_resolution =  step_size/sqrt(2) + 0.0130;// + 0*TMath::Sqrt( 0.18 / pt / pt + 0.99 ) / 100; // rescattering
                    if(TMath::Abs(dca)>0.9) continue;
                    n_tracks++;
                    //if (pt<0.2||dca>0.05) continue; // reject tracks with pt < 0.2 GeV/c or DCA > 500 μm
                    const int local_weight = 1. / TMath::Sqrt( SQR(rescattering_resolution) + SQR(dca) ); // weight based on DCA and pt
                    n_wtracks += local_weight;
                    if( n_tracks%2==0 ) n_tracks_east+=local_weight;//(int) (pt*1000)
                    else n_tracks_west+=local_weight;
                }
                if (n_tracks > maximum_tracks) maximum_tracks = n_tracks;
                n_tracks_vec.push_back(n_wtracks);
                track_vertices.push_back(std::make_pair(xvtx, yvtx));
                if(n_tracks_west)
                {
                    track_vertices_west.push_back(std::make_pair(xvtx, yvtx));
                    n_tracks_vec_west.push_back(n_tracks_west);
                }
                if (n_tracks_east)
                {
                    track_vertices_east.push_back(std::make_pair(xvtx, yvtx));
                    n_tracks_vec_east.push_back(n_tracks_east);
                }
            }// end of loop over yvtx
        }// end of loop over xvtx
        int max_tracks = *std::max_element(n_tracks_vec.begin(), n_tracks_vec.end());
        int min_tracks = *std::min_element(n_tracks_vec.begin(), n_tracks_vec.end());
        int track_treshold = (max_tracks - min_tracks) * fraction + min_tracks;
        if (n_tracks_vec.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            for (size_t i = 0; i < n_tracks_vec.size(); ++i)
            {
                if (n_tracks_vec[i] > track_treshold)
                {
                    {
                        hist_vtx_grid_xy->Fill(track_vertices[i].first, track_vertices[i].second, event->GetCentrality(), n_tracks_vec[i]);
                    }
                    const double weight = n_tracks_vec[i] - track_treshold;
                    sum_weights += weight;
                    sum_x += track_vertices[i].first * weight;
                    sum_y += track_vertices[i].second * weight;
                }
            }
            if (fill_hist && TMath::Abs(event->GetPreciseX()-beam_x)<0.0025 && TMath::Abs(event->GetPreciseY()-beam_y)<0.0025)
            {
                for (size_t i = 0; i < n_tracks_vec.size(); ++i)
                {
                    hist_vtx_grid_xy->Fill(track_vertices[i].first, track_vertices[i].second, event->GetCentrality(), n_tracks_vec[i]);
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x = new_vx;
            best_y = new_vy;
        }
        float best_x_west = 1;
        float best_y_west = 1;
        float best_x_east = -1;
        float best_y_east = -1;
        if(n_tracks_vec_east.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            int max_tracks_east = *std::max_element(n_tracks_vec_east.begin(), n_tracks_vec_east.end());
            int min_tracks_east = *std::min_element(n_tracks_vec_east.begin(), n_tracks_vec_east.end());
            int track_treshold_east = (max_tracks_east - min_tracks_east) * fraction + min_tracks_east;
            for (size_t i = 0; i < n_tracks_vec_east.size(); ++i)
            {
                if (n_tracks_vec_east[i] > track_treshold_east)
                {
                    const double weight = n_tracks_vec_east[i] - track_treshold_east;
                    sum_weights += weight;
                    sum_x += track_vertices_east[i].first * weight;
                    sum_y += track_vertices_east[i].second * weight;
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x_east = new_vx;
            best_y_east = new_vy;
        }
        if(n_tracks_vec_west.size())
        {
            double sum_weights = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;

            int max_tracks_west = *std::max_element(n_tracks_vec_west.begin(), n_tracks_vec_west.end());
            int min_tracks_west = *std::min_element(n_tracks_vec_west.begin(), n_tracks_vec_west.end());
            int track_treshold_west = (max_tracks_west - min_tracks_west) * fraction + min_tracks_west;
            for (size_t i = 0; i < n_tracks_vec_west.size(); ++i)
            {
                if (n_tracks_vec_west[i] > track_treshold_west)
                {
                    const double weight = n_tracks_vec_west[i] - track_treshold_west;
                    sum_weights += weight;
                    sum_x += track_vertices_west[i].first * weight;
                    sum_y += track_vertices_west[i].second * weight;
                }
            }
            double new_vx = sum_x / sum_weights;
            double new_vy = sum_y / sum_weights;
            best_x_west = new_vx;
            best_y_west = new_vy;
        }

        if (fill_hist)
        {
            hist_vtx_delta_x->Fill(best_x_west-best_x_east, best_x, event->GetCentrality());
            hist_vtx_delta_y->Fill(best_y_west-best_y_east, best_y, event->GetCentrality());
            hist_vtx_x->Fill(best_x - prev_x, best_x, event->GetCentrality());
            hist_vtx_y->Fill(best_y - prev_y, best_y, event->GetCentrality());
            hist_vtx_grid_xy->Fill(best_x, best_y, event->GetCentrality());
            hist_vtx_z->Fill(best_x, best_y, beam_z);
            hist_vtx_delta_x_reuse->Fill(best_x-best_x, best_x, event->GetCentrality());
            hist_vtx_delta_y_reuse->Fill(best_y-best_y, best_y, event->GetCentrality());
        }
        std::cout << "\033[32mBest X and Y vertex = " << best_x << " " << best_y << " with " << maximum_tracks << " "<< max_tracks << " " << min_tracks << " " << track_treshold << " tracks\033[0m" << std::endl;
        std::cout << "\033[32mBest X and Y vertex (east-west) = " << best_x_east<< " "<< best_x_west << " " << best_y_east<<" "<<best_y_west << " tracks\033[0m" << std::endl;

        event->SetPreciseX(best_x);
        event->SetPreciseY(best_y);
    }


}