#include "InvMass.h"
#include "TMath.h"
#include "TF1.h"

double Tsallis(const double pt, const double mass = 0.135, const double n = 8, const double T = 0.15, const double dNdy = 1)
{
  if (n < 1) return 1;
  const double mT = sqrt(pt*pt+mass*mass);
  return  dNdy * pt * (n - 1) * (n - 2) / (n*T + mass*(n - 1)) / (n*T + mass) * pow ( (n*T + mT) / (n*T + mass) , -n) ; //1./(2*TMath::Pi())
}

double Hagedorn(const double pt, const double dNdy = 95.7, const double mass = 0.135, const double A = 504.5, const double a = 0.5169, const double b = 0.1626, const double p0 = 0.7366, const double n = 8.274)
{
  if (mass < 0) return 1;
  const double mT = sqrt(pt*pt + mass*mass - 0.018218683); ///mpi0^2
  return  dNdy / dNdy_pp[0] * pt * A * pow( exp(-a * mT - b *mT * mT ) + mT / p0, -n) ;///2 * TMath::Pi() *
}
Double_t Hagedorn_Yield_Func(Double_t *x, Double_t *p) {
  //
  // Hagedorn function in 1/pt dN/pt from ppg088 
  //
    Double_t pt   = x[0];
    Double_t mass = p[0];
    Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.134*0.134);
    Double_t A    = p[1];
    Double_t a    = p[2];
    Double_t b    = p[3];
    Double_t p0   = p[4];
    Double_t n    = p[5];
   
    Double_t value = 2 * TMath::Pi() * pt *A* pow( exp(-a*mt-b*mt*mt)+mt/p0 , -n);
    return value;
  }
  
TF1 *fHagedorn(const Char_t *name, Double_t lowerlim,  Double_t upperlim,  const double dNdy = 95.7, const double mass = 0.135, const double A = 504.5, const double a = 0.5169, const double b = 0.1626, const double p0 = 0.7366, const double n = 8.274)
{
  TF1 *fHagedorn = new TF1(name, Hagedorn_Yield_Func, lowerlim, upperlim, 6);
  fHagedorn->SetParameters(mass, A, a, b, p0, n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
}

void InvMass(const TString inname = inFile[0],  int itread = 0, int ntreads = 1, int N_max = 1000000, const int part = 0)
{
  std::cout<<"start"<<std::endl;

  TF1* hadron_yield[5];TF1 *generated_spectrum = nullptr;
  for (int i = 0; i < 5; i++)
  {
    TString name = Form("hagedorr_yield_%d", i);
    hadron_yield[i] = fHagedorn(name, 0.01, 25, dNdy_pp[part], m_pp[part]);
    //if (part == 8 || part == 9 )  hadron_yield[i] = fHagedorn(name, 0.01, 20, dNdy_pp[part], m_pp[part], 377, 0.356, 0.068, 0.7,8.25);
    std::cout<<" integral " << i << " " << hadron_yield[i]->Integral(0.01,20) <<std::endl;
    //hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) / hadron_yield[i]->Integral(0.4,10) * dNdy_pp[part]/42.2*Ncolls[i] * Br_to_ee[part] );
    hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) * Ncolls[i] / 257. * Br_to_ee[part] );
    if (part<2||part==8) 
    {
      TString name1 = Form("hagedorr_new_yield_%d", i);
      hadron_yield[i] = fHagedorn(name1, 0.01, 25, dNdy_pp[part], m_pp[part], hag_params_try[i][0], hag_params_try[i][1], hag_params_try[i][2], hag_params_try[i][3], hag_params_try[i][4]);
      hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) * Br_to_ee[part] );
      //hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) * pi0_multiplicity_AuAu[i] / pi0_multiplicity_AuAu[5] * 257. / Ncolls[i] );
    }
    if (part == 8 || part == 9 ) hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) / hadron_yield[i]->Integral(0.01,10) * dNdy_pp[part] * Ncolls[i] * Br_to_ee[part] * RAA_jpsi[i] );
    if (part < 8 && part >0)
      hadron_yield[i]->SetParameter(1, hadron_yield[i]->GetParameter(1) * hith_pt_ratio[part] / (hadron_yield[i]->Eval(5.0) / Br_to_ee[part] / pi0_high_pt_ratio[i] ) );
    if (part == 0 && true)
    {
      generated_spectrum = fHagedorn(name, 0.01, 20, dNdy_pp[part], m_pp[part], 377, 0.356, 0.068, 0.7,8.25);
      generated_spectrum->SetParameter(1, generated_spectrum->GetParameter(1) / generated_spectrum->Integral(0.4,10) *10 );
    }
    std::cout<<" integral " << i << " " << hadron_yield[i]->Integral(0.01,25)/Br_to_ee[part] << " " << dNdy_pp[part]/42.2*Ncolls[i] <<" "<<Br_to_ee[part] << " " << hadron_yield[i]->Eval(5.0)/Br_to_ee[part]<<std::endl;
  }
  TF1 *trk_mult_centrality = new TF1("trk_mult_centrality", "100*exp([0]*x)", 0, 100);
  trk_mult_centrality->SetParameter(0, -0.0411197);
  const int mult_centr_bins[] = {0,20,40,60,80,93};
  float average_mult_per_centr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < 5; i++) {
    average_mult_per_centr[i] = trk_mult_centrality->Integral(mult_centr_bins[i], mult_centr_bins[i+1])/(mult_centr_bins[i+1]-mult_centr_bins[i]);
  }
  TFile *input = new TFile(inname, "READ");
  if (!(input))
  {
    std::cout << "no input file" << std::endl;
    return;
  }
  const int associate_hits = 1;
  const int remove_hadron_hits = 0;
  const int fill_QA_hadron_hists = 0;
  const int fill_QA_lepton_hists = 0;//1
  const int fill_TTree = 0;
  const int fill_d_dphi_hists = 0;
  const int fill_DCA_hists = 0;
  const int use_d_dphi_DCA = 0;
  const int do_track_QA = 0;
  const int do_electron_QA = 3;
  const int do_ident_electrons = 2;//2
  const int do_reveal_hadron = 0;
  const int Use_ident = 0;
  const int fill_true_DCA = 1;//1
  const int check_veto = 0;
  const int fill_inv_mass = 1;//0
  const int fill_inv_mass_sim = 1;//0
  const int fill_vtx_accaptance = 0;
  const int uncorrect_vtx_offset = 1;
  const int do_vertex_reco = 1;
  const int do_conv_dalitz_finder = 2;//2
  const int is_hagedorn = 0;
  const int do_pt_reco = 0;

  char outname[200];
  sprintf(outname,"kek_%d.root",itread);
  TTree *T = (TTree *)input->Get("tree");
  TBranch *br = T->GetBranch("MyEvent");
  MyDileptonAnalysis::MyEventContainer *event_container = new MyDileptonAnalysis::MyEventContainer();
  event_container->InitEvent();
  event_container->GetHistsFromFile("../../ee_QA/AnaTrain/Run14AuAuLeptonComby/field_map.root");
  event_container->CreateOutFileAndInitHists(outname,fill_QA_lepton_hists,fill_QA_hadron_hists,fill_TTree,fill_d_dphi_hists,
                                               fill_DCA_hists, do_track_QA+do_electron_QA, do_reveal_hadron, fill_true_DCA, check_veto, fill_inv_mass, do_vertex_reco, do_conv_dalitz_finder);
  MyDileptonAnalysis::MyEvent *myevent = 0;                                            
  //event = 0;
  br->SetAddress(&myevent);

  TH3D *hist_pt_orig = new TH3D("hist_pt_orig","hist_pt_orig",50,0,10,5,0,100,2,0,2);
  TH2D *hist_pt_mother = new TH2D("hist_pt_mother","hist_pt_mother",500,0,25,5,0,100);
  TH2D *hist_pt_mother_weight = new TH2D("hist_pt_mother_weight","hist_pt_mother_weight",500,0,25,100,0,100);
  TH3D *hist_dca_mother_3d = new TH3D("hist_dca_3d", "hist_dca_3d", 200, -500, 500, 50, 0, 10, 10, 0, 100);
  TH3D *n_hits_hist = new TH3D("n_hits_hist","n_hits_hist",1000,0,1000,100,0,100, 4, 0, 4);
  TH3D *hist_pt_reco = new TH3D("hist_pt_reco","hist_pt_reco",50,0,5,4,0,4,5,0,5);
  TH3D *hist_sing_pt = new TH3D("hist_sing_pt","hist_sing_pt",50,-0.5,0.5,100,0,10,5,0,100);
  TH3D *hist_sing_e = new TH3D("hist_sing_e","hist_sing_e",50,-0.5,0.5,100,0,10,5,0,100);
  TH3D *hist_sing_pt_reco = new TH3D("hist_sing_pt_reco","hist_sing_pt_reco",50,-0.5,0.5, 100, -10,10 , 50,0,5);

  const double ccbar_yield_ratio_fit_params[4] = { 7.213072954495502 , -5.168881501393766 , 1.427529678167659 , -0.15216750641396687 };
  TF1 *ccbar_yield_ratio_fit_pol3 = new TF1("ccbar_yield_ratio_fit_pol3", "pol3", 0, 4.5);
  ccbar_yield_ratio_fit_pol3->SetParameters(ccbar_yield_ratio_fit_params);

  const float phi_west_min = -3*TMath::Pi()/16, phi_west_max = 5*TMath::Pi()/16;
  const float phi_east_min = 11*TMath::Pi()/16, phi_east_max = 19*TMath::Pi()/16;
  const float k_DC = 0.206, k_RICH = 0.309; // rad GeV/c

  std::cout << "Trees read!" << std::endl;

  const int nevt = T->GetEntries();
  const int beggin = nevt * (itread - 1) / ntreads;
  const int endish = nevt * itread / ntreads;
  
  TF1 f("f","gaus",-5,5);
  f.SetParameter(0,1);
  f.SetParameter(1,0);
  f.SetParameter(2,1);
  const float vtx_res_x[10] = {11.573013950655056, 15.215734440960366, 21.771775161191638, 33.458415120576696, 53.807666598553, 60.86003298271064, 34.027076037074465, 58.177590775412874, 87.71163737802117};
  const float vtx_res_y[10] = {7.398283944195731, 9.657671603622507, 13.770485615750985, 21.207408416267437, 34.80991174557141, 40.69221130423994, 24.078456294746875, 42.90731740850407, 67.81600271377448};

  for (int ievent = beggin; ievent < endish; ievent++)
  {
    if ((ievent -beggin) % 50000 == 0)
    std::cout << "ithread, iEvent, N_events: " << itread<< ",  " << ievent -beggin<< " / " << nevt/ntreads << std::endl;
    myevent->ClearEvent();
    event_container->ClearEvent();
    br->GetEntry(ievent);
    if (ievent - beggin > N_max)
      break;
    //if (myevent->GetCentrality() > 30) continue;
    if(false)
    {
      int not_in_acceptance = 0;
      for (int itrk = 0; itrk < myevent->GetNgentrack(); itrk++)
      {
        MyDileptonAnalysis::MyGenTrack *mygentrk = myevent->GetGenTrack(itrk);
        if (mygentrk->GetID()==22) continue;
        if (mygentrk->GetID()>0)
          hist_pt_orig->Fill(mygentrk->GetPt(),myevent->GetCentrality(),0);
        else
          hist_pt_orig->Fill(mygentrk->GetPt(),myevent->GetCentrality(),1);
          for (int jtrk = itrk + 1; jtrk < myevent->GetNgentrack(); jtrk++)
        {
          MyDileptonAnalysis::MyGenTrack *mygentrk2 = myevent->GetGenTrack(jtrk);
          const double pt1 = mygentrk->GetPt();
          const double pt2 = mygentrk2->GetPt();
          if(pt1<0.2 || pt2<0.2) not_in_acceptance++;
          const double eta1 = 0.5 * log( ( sqrt(pt1*pt1 + mygentrk->GetPz()*mygentrk->GetPz()) + mygentrk->GetPz() ) / ( sqrt(pt1*pt1 + mygentrk->GetPz()*mygentrk->GetPz()) - mygentrk->GetPz() ) );
          const double eta2 = 0.5 * log( ( sqrt(pt2*pt2 + mygentrk2->GetPz()*mygentrk2->GetPz()) + mygentrk2->GetPz() ) / ( sqrt(pt2*pt2 + mygentrk2->GetPz()*mygentrk2->GetPz()) - mygentrk2->GetPz() ) );
          if (TMath::Abs(eta1)>0.35 || TMath::Abs(eta2)>0.35) not_in_acceptance++;
          int q1 = mygentrk->GetID()<0 ? 1 : -1;
          int q2 = mygentrk2->GetID()<0 ? 1 : -1;
          double phi1 = TMath::ATan2( mygentrk->GetPy(),  mygentrk->GetPx());
          double phi2 = TMath::ATan2(mygentrk2->GetPy(), mygentrk2->GetPx());
          double phi1_DC = phi1 + k_DC * q1 / pt1 ;
          double phi2_DC = phi2 + k_DC * q2 / pt2 ;
          if (phi1_DC < -TMath::Pi()/2) phi1_DC += 2*TMath::Pi();
          if (phi2_DC < -TMath::Pi()/2) phi2_DC += 2*TMath::Pi();
          if ( !((phi1_DC > phi_west_min && phi1_DC < phi_west_max) || (phi1_DC > phi_east_min && phi1_DC < phi_east_max)) ) not_in_acceptance++;
          if ( !((phi2_DC > phi_west_min && phi2_DC < phi_west_max) || (phi2_DC > phi_east_min && phi2_DC < phi_east_max)) ) not_in_acceptance++;
          double phi1_RICH = phi1 + k_RICH * q1 / pt1 ;
          double phi2_RICH = phi2 + k_RICH * q2 / pt2 ;
          if (phi1_RICH < -TMath::Pi()/2) phi1_RICH += 2*TMath::Pi();
          if (phi2_RICH < -TMath::Pi()/2) phi2_RICH += 2*TMath::Pi();
          if ( !((phi1_RICH > phi_west_min && phi1_RICH < phi_west_max) || (phi1_RICH > phi_east_min && phi1_RICH < phi_east_max)) ) not_in_acceptance++;
          if ( !((phi2_RICH > phi_west_min && phi2_RICH < phi_west_max) || (phi2_RICH > phi_east_min && phi2_RICH < phi_east_max)) ) not_in_acceptance++;
        }
      }
      if(not_in_acceptance) continue;
    }
    double weight = Ncolls[(int)myevent->GetCentrality()/20] * trk_mult_centrality->Eval(myevent->GetCentrality()) / average_mult_per_centr[(int)myevent->GetCentrality()/20];
    //if (part == 8 || part == 9 ) weight = RAA_jpsi[(int)myevent->GetCentrality()/20];
    if(part==10)
    {
      int bbbar_in_ccbar = 0;
      for(int igentrk=0; igentrk<myevent->GetNgentrack(); igentrk++)
      {
        MyDileptonAnalysis::MyGenTrack *mygentrk = myevent->GetGenTrack(igentrk);
        if (TMath::Abs((mygentrk->GetPx()-0.1))<0.01 && TMath::Abs((mygentrk->GetPy()-0.1))<0.01) bbbar_in_ccbar=1;
      }
      //if (bbbar_in_ccbar==1) std::cout<<"found bbbar in ccbar event "<<ievent<<std::endl;
      if (bbbar_in_ccbar==1) continue;
    }
    if (myevent->GetNgentrack()==2)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk1 = myevent->GetGenTrack(0);
      MyDileptonAnalysis::MyGenTrack *mygentrk2 = myevent->GetGenTrack(1);
      const double px = mygentrk1->GetPx()+mygentrk2->GetPx();
      const double py = mygentrk1->GetPy()+mygentrk2->GetPy();
      const double pt = sqrt(px*px+py*py);
      //weight = Tsallis(pt, m_pp[part], n_pp[part], T_pp[part], dNdy_pp[part]/42.2*Ncolls[(int)myevent->GetCentrality()/20]);
      //weight = Hagedorn(pt, dNdy_pp[part]*Br_to_ee[part], m_pp[part]);
      if ( part < 10) weight = hadron_yield[(int)myevent->GetCentrality()/20]->Eval(pt);
      if ( part == 0 && is_hagedorn) weight /= generated_spectrum->Eval(pt);
      weight *= trk_mult_centrality->Eval(myevent->GetCentrality()) / average_mult_per_centr[(int)myevent->GetCentrality()/20];
      hist_pt_mother->Fill(pt,myevent->GetCentrality());
      hist_pt_mother_weight->Fill(pt,myevent->GetCentrality(),weight);
      if(true && part==10)
      {
        const float E1 = sqrt(mygentrk1->GetPx()*mygentrk1->GetPx()+mygentrk1->GetPy()*mygentrk1->GetPy()+mygentrk1->GetPz()*mygentrk1->GetPz());
        const float E2 = sqrt(mygentrk2->GetPx()*mygentrk2->GetPx()+mygentrk2->GetPy()*mygentrk2->GetPy()+mygentrk2->GetPz()*mygentrk2->GetPz());
        float mass = sqrt( (E1+E2)*(E1+E2) - (px*px+py*py+(mygentrk1->GetPz()+mygentrk2->GetPz())*(mygentrk1->GetPz()+mygentrk2->GetPz()) ) );
        mass = mass < 0.5 ? 0.5 : mass >3.0 ? 3.0 : mass;
        weight *= ccbar_yield_ratio_fit_pol3->Eval(mass);
      }
    }
    if (myevent->GetNgentrack()==3)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk1 = myevent->GetGenTrack(0);
      MyDileptonAnalysis::MyGenTrack *mygentrk2 = myevent->GetGenTrack(1);
      MyDileptonAnalysis::MyGenTrack *mygentrk3 = myevent->GetGenTrack(2);
      const double px = mygentrk1->GetPx()+mygentrk2->GetPx()+mygentrk3->GetPx();
      const double py = mygentrk1->GetPy()+mygentrk2->GetPy()+mygentrk3->GetPy();
      const double pt = sqrt(px*px+py*py);
      //weight = Hagedorn(pt, dNdy_pp[part]*Br_to_ee[part], m_pp[part]);
      //weight = Tsallis(pt, m_pp[part], n_pp[part], T_pp[part], dNdy_pp[part]/42.2*Ncolls[(int)myevent->GetCentrality()/20]);
      if ( part < 10) weight = hadron_yield[(int)myevent->GetCentrality()/20]->Eval(pt);
      if ( part == 0 && is_hagedorn) weight /= generated_spectrum->Eval(pt);
      weight *= trk_mult_centrality->Eval(myevent->GetCentrality()) / average_mult_per_centr[(int)myevent->GetCentrality()/20];
      hist_pt_mother->Fill(pt,myevent->GetCentrality());
      hist_pt_mother_weight->Fill(pt,myevent->GetCentrality(),weight);
    }
    if (myevent->GetNgentrack() == 1)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk1 = myevent->GetGenTrack(0);
      const double px = mygentrk1->GetPx();
      const double py = mygentrk1->GetPy();
      const double pt = sqrt(px * px + py * py);
      // weight = Tsallis(pt, m_pp[part], n_pp[part], T_pp[part], dNdy_pp[part]/42.2*Ncolls[(int)myevent->GetCentrality()/20]);
      // weight = Hagedorn(pt, dNdy_pp[part]*Br_to_ee[part], m_pp[part]);
      if (part < 10)
        weight = hadron_yield[(int)myevent->GetCentrality() / 20]->Eval(pt)* trk_mult_centrality->Eval(myevent->GetCentrality()) / average_mult_per_centr[(int)myevent->GetCentrality()/20];
      hist_pt_mother->Fill(pt, myevent->GetCentrality());
      hist_pt_mother_weight->Fill(pt, myevent->GetCentrality(), weight);
    }
    //int energy_loss = 1;
    //if(myevent->GetNtrack()>1) continue;
    //else{
    //  for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    //  {
    //    MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
    //    const float pt = mytrk->GetPtPrime();
    //    for (int i = 0; i < myevent->GetNgentrack(); i++)
    //    {
    //      MyDileptonAnalysis::MyGenTrack *mygentrk = myevent->GetGenTrack(i);
    //      if(mygentrk->GetID()==22) continue;
    //      const float pt_gen = mygentrk->GetPt();
    //      if (pt_gen*0.8<pt) energy_loss = 0;
    //    }
    //  }
    //} 
    //if(energy_loss==0) continue;
    int contain_long_tracks = 0;
    for (int itrk = 0; itrk < myevent->GetNgentrack()*0; itrk++)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk = myevent->GetGenTrack(itrk);
      const int charge = mygentrk->GetID() > 0 ? -1 : 1;
      const double vx = mygentrk->GetVx()*1e-9 - myevent->GetPreciseX()*1e4;
      const double vy = mygentrk->GetVy()*1e-9 - myevent->GetPreciseY()*1e4;
      const double v = sqrt(vx*vx + vy*vy);
      const double phi_vtx = TMath::ATan2(vy, vx);
      const double px = mygentrk->GetPx();
      const double py = mygentrk->GetPy();
      const double pt = sqrt(px*px + py*py+1e-6);
      const double r = pt / 0.003 / 0.90;
      const double phi_pt = TMath::ATan2(py, px);
      const double dca_l = v * TMath::Sin(phi_vtx - phi_pt);
      const double dca_t = SQR(v * TMath::Cos(phi_vtx - phi_pt))/r;
      const double dca = sqrt(dca_l*dca_l*0 + dca_t*dca_t );
      if (dca>400) contain_long_tracks++;
      hist_dca_mother_3d->Fill(dca, mygentrk->GetPt(), myevent->GetCentrality());
    }
    //if (!contain_long_tracks) continue;
    const int icentrality = (int)myevent->GetCentrality()/10;
    const float event_vertex_sigma[3] = {vtx_res_x[icentrality] / 10000, vtx_res_y[icentrality] / 10000, 200. / 10000};
    myevent->SetPreciseX(myevent->GetPreciseX() + f.GetRandom() * event_vertex_sigma[0]);
    myevent->SetPreciseY(myevent->GetPreciseY() + f.GetRandom() * event_vertex_sigma[1]);
    myevent->SetPreciseZ(myevent->GetPreciseZ() + f.GetRandom() * event_vertex_sigma[2]);
    myevent->SetVtxZ(myevent->GetPreciseZ());
    //if(myevent->GetCentrality()>30) continue;
    if(fill_inv_mass_sim) event_container->fill_inv_mass_sim(weight);
    if (fill_inv_mass && myevent->GetNtrack()<2) continue;
    if(myevent->GetNtrack()<1 && !do_track_QA) continue;
    
    for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
      
      if (mytrk->GetArm()<2)mytrk->SetPtPrime(mytrk->GetPtPrime()*0.985);//0.97
      //const double pt_smear = 1.+0.025*f.GetRandom();/// 0.7+1*pt
      //mytrk->SetPt(mytrk->GetPt()*pt_smear);
      //mytrk->SetPtPrime(mytrk->GetPtPrime()*pt_smear);
      if(mytrk->GetEcore()<0.1) continue;
      const int sector = mytrk->GetArm() == 1 ? mytrk->GetSect() : 7 - mytrk->GetSect();//mytrk->GetSect();
      const double e_smear = scale[sector] *  ( 1 + f.GetRandom() * sqrt( pow ( smear_c1[sector] , 2 ) + pow( smear_c2[sector] / sqrt(mytrk->GetEcore()), 2) ) );          // (1.+0.025*f.GetRandom());//// 8.1/sqrt(E)+2.1
      mytrk->SetEcore(mytrk->GetEcore()*e_smear);

      //const float dphi_smear = - 0.001 * mytrk->GetChargePrime() + 0.0005 + 0.0025 /mytrk->GetPtPrime() * f.GetRandom();
      //mytrk->SetPhi0Prime(mytrk->GetPhi0Prime() + dphi_smear);
      //mytrk->SetPhi0(mytrk->GetPhi0Prime() + dphi_smear);
      //mytrk->SetPtPrime(mytrk->GetPtPrime() * mytrk->GetAlphaPrime() / (mytrk->GetAlphaPrime() + dphi_smear ));
    }

    int n_electrons = myevent->GetNtrack();
    if (associate_hits && do_ident_electrons)
    {
      for (int itrk = 0; itrk < n_electrons; itrk++)
      {
        MyDileptonAnalysis::MyElectron mytrk = *myevent->GetEntry(itrk);
        bool skip = false; 
          if(mytrk.GetTrkQuality()!=63 && mytrk.GetTrkQuality()!=51 && mytrk.GetTrkQuality()!=31)
            skip = true;
          if (mytrk.GetPtPrime()>10 || mytrk.GetPtPrime() < 0.3)
            skip = true;
          if (mytrk.GetCrkphi()<-99)
            skip = true;
          if (mytrk.GetEcore()/mytrk.GetPtot() < 0.1)///0.6 podgon
            skip = true;
          if (mytrk.GetN0() < 0 )
            skip = true;
          if (fabs(mytrk.GetEmcdphi())>0.05 || fabs(mytrk.GetEmcdz())>25 )
            skip = true;
      
          if( skip ){
            myevent->RemoveTrackEntry(itrk);
            //event->AddElecCand(&mytrk);
            n_electrons--;
            itrk--;
            continue;
          }
          //std::cout<<mytrk.GetPt()<<" "<<mytrk.GetPtPrime()<<" "<<mytrk.GetReconPT()<<std::endl;
      }
    }
    if (do_ident_electrons == 2)
    {
      event_container->SetEvent(myevent);
      event_container->IdenElectrons(1);
      event_container->FillQAHistPreAssoc(weight);
      n_electrons = myevent->GetNtrack();
      for (int itrk = 0; itrk < n_electrons; itrk++)
      {
        MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
        //if (mytrk->GetPtPrime() > 4.4 && mytrk->GetProb() > 0.8 && mytrk->GetEcore() / mytrk->GetPtot() > 0.5 &&
        //    mytrk->GetEcore() / mytrk->GetPtot() < 1.4 && TMath::Abs(mytrk->GetEmcdphi()) < 0.005 * 3 && TMath::Abs(mytrk->GetEmcdz()) < 10 && mytrk->GetN0() > 1 && mytrk->GetDisp() < 5)
          //if (mytrk->GetPtPrime() > 4.4 && (mytrk->GetEcore() / mytrk->GetPtot() < 0.5 || TMath::Abs(mytrk->GetEmcdz()-1) > 8) )
          if(mytrk->GetPtPrime() > 4.6 && (mytrk->GetEcore()/mytrk->GetPtot()<0.6 || TMath::Abs(mytrk->GetEmcdz()-1) > 8  || mytrk->GetProb() < 0.5
             || mytrk->GetEcore()/mytrk->GetPtot()>1.4 ||TMath::Abs(mytrk->GetEmcdphi()) > 0.005*3 || mytrk->GetN0()<1 || mytrk->GetDisp()>5 ) )
          mytrk->SetMcId(0); // 1000 is the electron id
      }
      for (int itrk = 0; itrk < n_electrons; itrk++)
      {
        MyDileptonAnalysis::MyElectron mytrk = *myevent->GetEntry(itrk);
        if ( (myevent->GetCentrality()<40 && mytrk.GetMcId() + mytrk.GetTrkInfo() < 10000) || (myevent->GetCentrality()>=40 && mytrk.GetTrkInfo() < 100) )
      //if ( mytrk.GetMcId()<100 || (myevent->GetCentrality()<40 && mytrk.GetMcId()<1000) || (myevent->GetCentrality()<20 && mytrk.GetMcId()<1000) || (mytrk.GetProb()<0.1 && myevent->GetCentrality()<40) )
        //if (mytrk.GetMcId() < 100 || (myevent->GetCentrality() < 40 && mytrk.GetMcId() < 1000) || (myevent->GetCentrality() < 20 && mytrk.GetMcId() < 1000) ||
        ////if( mytrk.GetMcId() % 10 < 6 ||
        //    (mytrk.GetPtPrime() < 0.4 && (fabs(mytrk.GetEmcdphi()) > 0.02 || fabs(mytrk.GetEmcdz()) > 8 || mytrk.GetDisp() > 3 || mytrk.GetMcId() % 10 < 6))) // adding regualr electron cuts|| mytrk.GetEcore()<0.3 || mytrk.GetEcore()/mytrk.GetPtot()<0.8
        {
          myevent->RemoveTrackEntry(itrk);
          // event->AddElecCand(&mytrk);
          n_electrons--;
          itrk--;
          continue;
        }
      }
    }
    if (fill_vtx_accaptance) {
      event_container->SetEvent(myevent);
      event_container->FillVTXAcceptance();
      continue;
    }
    if (fill_inv_mass && myevent->GetNtrack()<2) continue;
    if(myevent->GetNtrack()<1 && !do_track_QA) continue;
    
    int n_hits = myevent->GetNVTXhit();
    int h_hits_mc[2] = {0, 0};
    int n_embeded_hits[4] = {0, 0, 0, 0};
    for (int ihit = 0; ihit < n_hits; ihit++)
    {
      MyDileptonAnalysis::MyVTXHit myhit = *myevent->GetVTXHitEntry(ihit);
    if(myhit.GetLadder()>24)std::cout<<myhit.GetLadder()<<std::endl;
    //if (myhit.GetSensor() == 0) std::cout << "0 sesnsor:  " << myhit.GetLayer() << " " << myhit.GetPhiHit() << " " << myhit.GetZHit() << " " <<myhit.GetLadder()<<std::endl; 
    //if (myhit.GetSensor() == 1) std::cout << "1 sesnsor:  " << myhit.GetLayer() << " " << myhit.GetPhiHit() << " " << myhit.GetZHit() << " " <<myhit.GetLadder()<<std::endl;
     //if (myhit.GetSensor() == 0) h_hits_mc[0]++;
     //if (myhit.GetSensor() == 1) h_hits_mc[1]++;
     if (myhit.GetSensor() == 1)
     {
         myevent->RemoveVTXHitEntry(ihit);
         n_hits--;
         ihit--;
         continue;
     }
    if(myhit.GetSensor() > 1)
    {
      n_embeded_hits[myhit.GetLayer()]++;
    }
     if (true||!fill_inv_mass) continue;
     if (myhit.GetLayer()==3 && myhit.GetPhiHit(0,0,0)>0.85&& myhit.GetPhiHit(0,0,0)<1.08)
     {
         myevent->RemoveVTXHitEntry(ihit);
         n_hits--;
         ihit--;
         continue;
     }
     if (myhit.GetLayer()==2 && myhit.GetPhiHit(0,0,0)>3.17&& myhit.GetPhiHit(0,0,0)<3.45)
     {
         myevent->RemoveVTXHitEntry(ihit);
         n_hits--;
         ihit--;
         continue;
     }
    }
    //if(h_hits_mc[0]+5<h_hits_mc[1]) std::cout<<"Warning: hits in sensor 0 less than in 1 sensor: " << h_hits_mc[0] << " " << h_hits_mc[1] << std::endl;
    for (int ilayer = 0; ilayer < 4; ilayer++) {
        n_hits_hist->Fill(n_embeded_hits[ilayer], myevent->GetCentrality(), ilayer);
    }

    if(do_ident_electrons<2) 
      event_container->SetEvent(myevent);

    event_container->correct_beam_offset();
    bool skip_event = false;
    for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    {
        MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
        mytrk->ZeroHitCounters();
        if (mytrk->GetMcId()!=9) skip_event = true;
    }
    //if(skip_event) continue;

    if (uncorrect_vtx_offset) event_container->CorrectVTXOffset(-1);
    for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    {
      MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
      const float dphi_smear = - 0.001 * mytrk->GetChargePrime()*0 + 0.0005*0 + 0.005 /mytrk->GetPtPrime() * f.GetRandom();
      mytrk->SetPhi0Prime(mytrk->GetPhi0Prime() + dphi_smear);
      mytrk->SetPhi0(mytrk->GetPhi0Prime());
      mytrk->SetDCAX(mytrk->GetPhi0());
      mytrk->SetPtPrime(mytrk->GetPtPrime() * mytrk->GetAlphaPrime() / (mytrk->GetAlphaPrime() + dphi_smear / 2  ));
      mytrk->SetAlphaPrime(mytrk->GetAlphaPrime() + dphi_smear / 2);
    }
    
    if(do_ident_electrons==1) event_container->IdenElectrons(1);
    if(do_electron_QA && do_ident_electrons<2)event_container->FillQAHistPreAssoc(weight);
    //if(do_vertex_reco) event_container->VertexXYScan(myevent->GetPreciseX(), myevent->GetPreciseY(), 1,0);
    if (false) {
      event_container->Associate_Hits_to_Leptons(5,5,5,1,2,3);
      event_container->Associate_Hits_to_Leptons(5,5,5,!fill_QA_lepton_hists,0,3,5,weight);
      for (int itrk = 0; itrk < n_electrons; itrk++)
      {
          MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
          if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
             (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 ))
          {
              myevent->RemoveTrackEntry(itrk);
              n_electrons--;
              itrk--;
              continue;
          }
          if(mytrk->GetPtPrime()<0.4)
          {
            if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || (mytrk->GetHitCounter(2) < 1 || mytrk->GetHitCounter(3) < 1))
              {
                myevent->RemoveTrackEntry(itrk);
                n_electrons--;
                itrk--;
                continue;
              }
          }
      }
      event_container->FillTrueDCA(weight);
      continue;
    }
    if(associate_hits) event_container->Associate_Hits_to_Leptons(5,5,5,1,2,3);
    if(associate_hits ) event_container->Associate_Hits_to_Leptons(5,5,5,1,1,3);
    if(associate_hits && event_container->GetNGoodElectrons()<1) continue;
    if(associate_hits && fill_inv_mass && event_container->GetNGoodElectrons()<2  ) continue;
    if(do_pt_reco) event_container->CorrectPtForEventOffset(myevent->GetPreciseX(), myevent->GetPreciseY(), 0);
    if(associate_hits ) event_container->Associate_Hits_to_Leptons(5,5,5,1,1,3);
    if(do_pt_reco) event_container->CorrectPtForEventOffset(myevent->GetPreciseX(), myevent->GetPreciseY(), -2);
    if(false) event_container->Associate_Hits_to_Leptons_OLD(20,20,20);
    if(do_conv_dalitz_finder) event_container->ConversionFinder((int)(do_conv_dalitz_finder==2),0,0,weight);
    if(associate_hits && fill_inv_mass && event_container->GetNGoodElectrons()<1  ) continue;
    if(associate_hits && event_container->GetNGoodElectrons()>=1) event_container->Associate_Hits_to_Leptons(5,5,5,!fill_QA_lepton_hists,0,3,5,weight);
    if(do_track_QA) event_container->FillQAHist(in_id, weight);

    
    if(remove_hadron_hits) 
    {
        event_container->Associate_Hits_to_Hadrons();
        event_container->Associate_Hits_to_Leptons();
    }

    n_electrons = myevent->GetNtrack();
    if (associate_hits)
    {
      for (int itrk = 0; itrk < n_electrons; itrk++)
      {
          MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
          if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || 
             (mytrk->GetHitCounter(2) < 1 && mytrk->GetHitCounter(3) < 1 ))
          {
              myevent->RemoveTrackEntry(itrk);
              n_electrons--;
              itrk--;
              continue;
          }
          if(mytrk->GetPtPrime()<0.4)
          {
            if (mytrk->GetHitCounter(0) < 1 || mytrk->GetHitCounter(1) < 1 || (mytrk->GetHitCounter(2) < 1 || mytrk->GetHitCounter(3) < 1))
              {
                myevent->RemoveTrackEntry(itrk);
                n_electrons--;
                itrk--;
                continue;
              }
          }
      }
    }

    if(false) event_container->CorrectPtForEventOffset(myevent->GetPreciseX(), myevent->GetPreciseY());
    if(associate_hits && event_container->GetNGoodElectrons()<1) continue;
    if(associate_hits && event_container->GetNGoodElectrons()>1) event_container->ResetRecoverFGVars(); 
    
    if(true)
    {
      std::vector<std::vector<float> > clusters;
      for(int ihadron = 0; ihadron < myevent->GetNhadron(); ihadron++)
      {
        MyDileptonAnalysis::MyHadron *myhadron = myevent->GetHadronEntry(ihadron);
        //if (myhadron->GetMcId()==1)continue;
        if( myhadron->GetProb()<0.1) continue;
        if( myhadron->GetEcore()<0.3) continue;
        clusters.push_back({myhadron->GetPhi0(), myhadron->GetThe0(), myhadron->GetZDC(), myhadron->GetEcore(), (float) myhadron->GetSect()});
      }
      event_container->Find_Bremsstrahlung(clusters, weight);
    }

    for (int itrk = 0; itrk < myevent->GetNgentrack(); itrk++)
    {
      MyDileptonAnalysis::MyGenTrack *mygentrk = myevent->GetGenTrack(itrk);
      if (mygentrk->GetID()==22) continue;
      if (mygentrk->GetID()>0)
      {
        for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
        {
          MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
          if(mytrk->GetChargePrime()>0) continue;
          hist_sing_pt->Fill((mytrk->GetPtPrime()-mygentrk->GetPt())/mygentrk->GetPt(), mygentrk->GetPt(), myevent->GetCentrality(),weight);
          hist_sing_e->Fill((mytrk->GetEcore()*-mygentrk->GetPt())/mygentrk->GetPt(), mygentrk->GetPt(), myevent->GetCentrality(),weight);
          float dphi = (mytrk->GetMinsDphi(1)+mytrk->GetMinsDphi(0));
          if(mytrk->GetDep())hist_sing_pt_reco->Fill((mytrk->GetPtPrime()-mygentrk->GetPt())/mygentrk->GetPt(), dphi, mygentrk->GetPt(),weight);
        }
      }else{
        for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
        {
          MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
          if(mytrk->GetChargePrime()<0) continue;
          hist_sing_pt->Fill((mytrk->GetPtPrime()-mygentrk->GetPt())/mygentrk->GetPt(), mygentrk->GetPt(), myevent->GetCentrality(),weight);
          hist_sing_e->Fill((mytrk->GetEcore()*TMath::Abs(TMath::Sin(mytrk->GetThe0Prime()))-mygentrk->GetPt())/mygentrk->GetPt(), mygentrk->GetPt(), myevent->GetCentrality(),weight);
          float dphi = (mytrk->GetMinsDphi(1)+mytrk->GetMinsDphi(0));
          if(mytrk->GetDep())hist_sing_pt_reco->Fill((mytrk->GetPtPrime()-mygentrk->GetPt())/mygentrk->GetPt(), dphi, mygentrk->GetPt(),weight);
        }
      }
    }

    //for (int itrk = 0; itrk < myevent->GetNtrack(); itrk++)
    //{
    //  MyDileptonAnalysis::MyElectron *mytrk = myevent->GetEntry(itrk);
    //  mytrk->SetPhi0(mytrk->GetPhi0() + 0.005*(mytrk->GetArm() == 0 ? 1 : -1));
    //}
//

    if(use_d_dphi_DCA)  event_container->FillDphiHists();
    if(do_reveal_hadron) event_container->Reveal_Hadron();
    //event_container->CheckVeto();
    if(fill_true_DCA) event_container->FillTrueDCA(weight);
    if(fill_TTree) event_container->FillTree();
    if(do_electron_QA)event_container->FillQAHist(in_id, weight);
    //myevent->ReshuffleElectrons();
    if(associate_hits && fill_inv_mass && event_container->GetNGoodElectrons()<2  ) continue;
    if(fill_inv_mass) event_container->fill_inv_mass(weight);
    myevent->ClearEvent();
    event_container->ClearEvent();
    //event->ClearEvent();
  }
  hist_pt_orig->Write();
  hist_pt_mother->Write();
  hist_pt_mother_weight->Write();
  hist_dca_mother_3d->Write();
  n_hits_hist->Write();
  hist_pt_reco->Write();
  hist_sing_pt->Write();
  hist_sing_e->Write();
  hist_sing_pt_reco->Write();
  event_container->WriteOutFile();
}


///backup in unout