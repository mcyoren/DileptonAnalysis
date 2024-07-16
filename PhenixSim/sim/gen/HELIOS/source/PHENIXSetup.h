/////// PHENIXSetup.h //////////////////////////////////////////////////////////////////////////////////////////
//
// Collection of constants defining PHENIX Setup 
// used in Class PHENIX 
//
// Axel Drees 9/14/2021
//  
#ifndef PHENIXSetup_h
#define PHENIXSetup_h

#define FIELD 0                   // Magnetic field 0 = ++; 1 = +/- ; 2 = + field configuration
#define RUN 14                    // run number

#include "PDG.h"

// Detector Identifiers 
const int arm_west = 1;
const int arm_east = 2;

// central arm acceptance 
// physical phi and eta acceptance of arms
const double P_eta = 0.35;     // rapidity range
const double P_theta_min = 2*atan(exp(-P_eta));
const double P_theta_max = 2*atan(exp(P_eta));
Double_t pi = 3.14159;                              // define pi
const double P_phi_west_max =  11./16.*pi;
const double P_phi_west_min =   3./16.*pi;
const double P_phi_east_max =  -3./16.*pi; 
const double P_phi_east_min = -11./16.*pi;

// radial location of detector in meter
const double P_R_DC = 2.20;          
const double P_R_EMCal = 5.00;

// geometry of 8 calorimeter sectors
// index scheme for sector number
//
//               E3   /5   4\   W3
//               E2  |6     3|  W2
//               E1  |7     2|  W1
//               E0   \8   1/   W0
//
// center of sector in unites of pi
// the phi coordinate is choosen to have phi=0 up
//
const double P_centerPhi[8] = {1./2.+2./15.9,  1./2.,          1./2.-2./15.9,  1./2.-4./15.9, 
	                          -1./2.+4./15.9, -1./2.+2./15.9, -1./2.,         -1./2.-2/16.36 };
//const double P_centerPhi[8] = {1./2.+2./16.,  1./2.,          1./2.-2./16.,  1./2.-4./16., 
//	                          -1./2.+4./16., -1./2.+2./16.,  -1./2.,        -1./2.-2./16. };
// half width of sector in units of pi
const double P_dPhi[8] = {1./16.26,1./16.26,1./16.26,1./16.26,
	                                1./16.26,1./16.26,1./17.48 ,1./17.48  };   
// dito excluding outer towers 
const double P_fdPhi[8] = {1./17.21,1./17.21,1./17.21,1./17.21,
	                                 1./17.21,1./17.21,1./18.24 ,1./18.24  };   

#if FIELD == 0
// magnetic field: the full bent of a particle traversing the magnetic field is approximately
// 	 315 mrad/pt measured in GeV for the ++ field configuration. This is split between alpha 
//   and phi_0 - phi_DC about 1/3 to 2/3
// 	                                
const double  P_kDC = 0.307;                 // this is the full bent angle in the field at the DC
const double  P_kBfield = 1.023*P_kDC;       // bent in full field (i.e. out to 3 meters) ~ .322      
const double  P_delta = 0.35;                // fraction of field bent measured by alpha - also valid for conversions in VTX          
// phi component of magnetic field bent at DC and RICH (assume RICH~EMCAL)
const double P_kDC_phi = (1-P_delta)*P_kDC;  // = 0.210
const double P_kRICH_phi = P_kBfield;        // is this needed ?? 
const double P_kEMCal_phi = 0.263;           // needs to be checked
const double P_kDCtoEMC = 0.0097;            // residual bent from DC to EMC 
const double P_deltaDCtoEMC = 0.081;         // same as P_delta 
// 
const double P_kVTX[4] = {0.300,0.294,0.280,0.254};  // field bend for track startin at vtx B0, B1, B2, B3

#endif


// DC resolution
const double sigma_alpha_ms  = 0.0007;   // rad GeV
const double sigma_alpha_det = 0.0011;   // rad
const double sigma_p_c1      = sigma_alpha_ms/P_kDC/P_delta;  // = 0.67%
const double sigma_p_c2      = sigma_alpha_det/P_kDC/P_delta; // = 1.05%
const double sigma_phi_ms    = 0.0043;   // rad GeV - all material including all VTX layers
const double sigma_eta_det   = 0.00068;  // rad

// EMCal PbSc resolution
const double sigma_E_PbSc_c1 = 0.081;     // %
const double sigma_E_PbSc_c2 = 0.021;     // %
const double sigma_x_PbSc_c1 = 0.0059;    // mm
const double sigma_x_PbSc_c2 = 0.0016;    // mm
const double sigma_x_PbSc_d  = 0.8*0.021;  // mm

// EMCal PbSc resolution
const double sigma_E_PbGl_c1 = 0.059;     // %
const double sigma_E_PbGl_c2 = 0.008;     // %
const double sigma_x_PbGl_c1 = 0.0084;    // mm
const double sigma_x_PbGl_c2 = 0.0002;    // mm
const double sigma_x_PbGl_d  = 0.0;       // mm

//
// live detector areas
#if RUN == 14
// EMCals Au-Au run 14
const double P_SectorLive[8] = {0.1242, 0.113, 0.1551, 0.1956, 0.1273, 0.1493, 0.2398, 0.3431};
//const double P_SectorLive[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
const double EMCProbEff = 0.1;
#endif


// detector material used for Bremsstrahlung and conversions
#if RUN > 10
const double VTX_X0[4] = {0.013,0.013,0.052,0.052}; 
const double VTX_R[4]  = {2.5,5.,10.,20.}; 
#else
const double VTX_X0[4] = {0.,0.,0.,0.};
const double VTX_R[4] = {0.,0.,0.,0.};
#endif

#endif
