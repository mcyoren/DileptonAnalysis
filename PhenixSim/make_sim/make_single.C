class TMyRandom {

    public:
        TMyRandom(int seed=987632621){
            fourVec = new TLorentzVector();
            unif    = new TRandom3(seed);
            fhagdorn = new TF1("fhagdorn","[0]*x/(pow(exp(-[1]*x-[2]*x*x)+x/[3],[4]))",0.0,100.0);
            fpow = new TF1("fpow","pow(x,[0])",0.0,100.0);
            pi = 3.14159;
        }

        double Rndm(){return unif->Rndm();}

        double Gaus(double mean, double sigma){ return unif->Gaus(mean, sigma); }

        double GetPowLaw(double n, double lpt, double upt){
            fpow->SetRange(lpt,upt);
            fpow->SetParameter(0,-n);
            return fpow->GetRandom();
            // return pow(pow(lpt,1-n)+unif->Rndm()*(pow(upt,1-n) - pow(lpt,1-n)), 1/(1-n) );
        }

        double GetHagdorn(double lpt, double upt)
        {
          fhagdorn->SetRange(lpt,upt);
          //double a = 0.5169, b = 0.1626, c = 504.5, p0 = 0.7366, n = 8.274;
          double a = 0.31602, b = 0.124106, c = 31.0164, p0 = 0.805548, n = 8.3051;
          fhagdorn->SetParameters(c,a,b,p0,n);
          return fhagdorn->GetRandom();
        }

        TLorentzVector *GetFMomGaussYPowPT(double ySig, double rapwin, double n, double lpt, double upt, double m0){
            do y  = unif->Gaus(0, ySig); while(fabs(y)>rapwin);
            if (n>0) pt = GetPowLaw(n, lpt, upt);
            else if (n==0) pt = lpt+(upt-lpt)*Rndm();
            else pt = GetHagdorn(lpt,upt);
            phi= unif->Rndm()*2*pi-pi; //-pi; ???
            mt = sqrt(m0*m0+pt*pt);
            pl = mt*sinh(y);
            fourVec->SetPxPyPzE(pt*cos(phi), pt*sin(phi), pl, sqrt(mt*mt+pl*pl) );
            return fourVec;
        }
        
    private:
        TLorentzVector *fourVec;
        TRandom3 *unif;
        TF1* fgaus;
        TF1 *fhagdorn;
        double phi,y,pt,pl,mt;
        double phiPh,yPh,ptPh,plPh;
        double pi;
        double gphi,z,gtheta,px,py,pz;

};

unsigned int get_seed()
{
    unsigned int fSeed;

    ifstream devrandom;
    devrandom.open("/dev/urandom",ios::binary);
    devrandom.read((char*)&fSeed,sizeof(fSeed));
    devrandom.close();
    if ( fSeed != -1 )
    {
      cout << " Got seed from /dev/urandom" << endl;
      fSeed = abs(fSeed)%900000000;
    }
    else fSeed = 0;
    cout << "seed is " << fSeed << endl;
    return fSeed;
}

void make_single(TString fout = "oscar.particles.dat", const int nevt = 10000, const float pt_min = 0.0, const float pt_max = 20, const double n = 0) {

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
const int NPARTICLES = 1;
double scale = 1e13; //cm to fm conversion
// ---------------------------

double rapwidth=10.,rapwin=.5;
Double_t IDpi0[3]     = { 111, 0, 0.1349766};
Double_t IDpiplus[3]  = { 211, 0, 0.13957018};
Double_t IDpiminus[3] = {-211, 0, 0.13957018};
Double_t IDkaonplus[3]  = { 321, 0, 0.493677};
Double_t IDkaonminus[3] = {-321, 0, 0.493677};
Double_t IDprotonplus[3]  = { 2212, 0, 0.938272};
Double_t IDprotonminus[3] = {-2212, 0, 0.938272};
//Double_t IDphoton[3] = {22, 0, 0.};

int seed = get_seed();
TMyRandom* myrand = new TMyRandom(seed);
gRandom->SetSeed(seed);
cout<<"double check seed "<<gRandom->Uniform()<<endl;

// Write out the file
ofstream fileout( fout );
fileout << "# OSC1999A" << endl;
fileout << "# final_id_p_x" << endl;
fileout << "# SimName 1.0" << endl;
fileout << "#" << endl;
fileout << "# Some comments..." << endl;

TF1* fGaus = new TF1("fGaus","gaus",-10,10);
 fGaus->SetParameters(1.27043e+04, 1.77852e+00, 1.27954e+01); // for CuAu

TH1D h1d_zvtx("h1d_zvtx","",400,-20,20);
TH1D h1d_pi0pt("h1d_pi0pt","",600,0,15);

// Get Zvertex dist from data
TH1D* zvtx_dat;

TFile* f1 = new TFile("heau200_bbcz.root", "READ");
zvtx_dat = (TH1D*)f1->Get("h2d_bbcz_px");
zvtx_dat->AddDirectory(0);

for(int ievt=0; ievt<nevt; ievt++) {

    //  Number of particles in event
    fileout << "0 " << NPARTICLES <<endl;

  for(int np=0; np<NPARTICLES; np++) {

	// Try to generate flat in zvtx (-30; 30) and reweight after
    double vertex = scale*fGaus->GetRandom();
    //double vertex = scale*zvtx_dat->GetRandom();
    //double vertex = scale*(60*myrand->Rndm()-30);
    h1d_zvtx.Fill(vertex/scale);

    //-------------------  PI0 --------------------
    TLorentzVector *PiZero = new TLorentzVector();
    PiZero = myrand->GetFMomGaussYPowPT(rapwidth, rapwin, n, pt_min, pt_max, IDpiplus[2]);
    h1d_pi0pt.Fill(PiZero->Pt());
    
    // Particle index
    Int_t index = np;
    // note that these positions are in femtometers !
    Double_t xpos=0.0; Double_t ypos=0.0; Double_t zpos=vertex; Double_t time=0.0;
    // idpart id ist px,py,pz,E, mass, x,y,z,t

    //-------------------  PI0 --------------------
    fileout << index << " " << IDpiplus[0] << " " << IDpiplus[1] << " " << PiZero->Px() << " " << PiZero->Py() << " " << PiZero->Pz() << " " << PiZero->E() << " " <<
  IDpiplus[2] << " " << xpos << " " << ypos << " " << zpos << " " << time << endl;
	
  }//End of particle
  fileout << "0 0" << endl;
}//End of events

fileout.close();

TFile fdiag("diag.root","recreate");
h1d_zvtx.Write();
h1d_pi0pt.Write();
fdiag.Write();
fdiag.Close();
}
