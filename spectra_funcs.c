#ifndef spectra_funcs_h
#define spectra_funcs_h

#include <stdio.h>
#include <math.h>
#include <gsl_interp2d.h>
#include <gsl_integration.h>
#include <gsl_spline.h>

//extern gsl_integration_workspace * w;
//#pragma omp threadprivate(w)
//gsl_integration_workspace * w;
extern double p;
extern double c_cmsm1;
extern double pc_cm;

double m_p = 0.9383; //GeV
double m_pi = 0.1350; //GeV
double m_e = 0.511e-3; //GeV
#define T_p_th (2. * m_pi + pow(m_pi, 2)/(2. * m_p)) /* GeV */

/* T_min T_max - minimum and maximum CR energy to integrate over */
double T_p_low = 1.; //GeV
double T_p_high = 1e6; //GeV
double mb_cm2 = 1e-27; //mbarn in cm^2




double alpha = 1./137.;
double hbar_GeVs = 6.582e-25;
double kb_GeVKm1 = 8.617e-14;
double arad_GeVcmm3Km4 = 4.722e-12; //GeV cm^-3 K^-4
double Lsol_GeVsm1 = 2.3882e36; //GeV s^-1

/*
typedef struct {
  double sigma_g;
  double h_pc;
  double * fcal;
  double tau_eff;
  double D_cm2sm1;
  } galaxy;
*/

typedef struct fdata_in {
  double E_g;
  size_t nx;
  size_t ny;
  double * xa;
  double * ya;
  double * za;
  const gsl_interp2d_type *T;
  gsl_interp2d *interp;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
} fd_in;



/* All energies in GeV */

/* Yu et al. fit to Manga data */
double sigma_gas_Yu( double SFR_Msolyrm1 ){return 32.063 * pow( SFR_Msolyrm1, 0.096 );}
/* Inverse KS derived from PHIBSS2 data */
double Sigma_gas_Msolpcm2_iKS( double Sigma_SFR_Msolyrm1pcm2 ){return 2.511e8 * pow( Sigma_SFR_Msolyrm1pcm2, 0.936 );}

double s( double T_p ){return 2. * m_p * ( T_p + 2. * m_p );}
double gam_CM( double T_p ){return ( T_p + 2. * m_p )/( sqrt( s( T_p ) ) );}
double E_pi_CM( double T_p ){return ( s( T_p ) - 4. * pow(m_p, 2) + pow(m_pi, 2) )/( 2. * sqrt( s( T_p ) ) );}
double P_pi_CM( double T_p ){return sqrt( pow(E_pi_CM( T_p ), 2) - pow(m_pi, 2) );}
double beta_CM( double T_p ){return sqrt( 1. - pow(gam_CM( T_p ), -2) );}
double E_pi_LAB_max( double T_p ){return gam_CM( T_p ) * ( E_pi_CM( T_p ) + P_pi_CM( T_p ) * beta_CM( T_p ) );}



double A_max( double T_p ){
  double b0 = 5.9;
  double b1( double T_p ){if (T_p < 5.){return 9.53;} else {return 9.13;}}
  double b2( double T_p ){if (T_p < 5.){return 0.52;} else {return 0.35;}}
  double b3( double T_p ){if (T_p < 5.){return 0.054;} else {return 9.7e-3;}}
  double theta_p = T_p/m_p;

  double sigma_1pi_mb( double T_p ){
    double sigma_0 = 7.66e-3;
    double M_res = 1.1883;
    double Gam_res = 0.2264;
    double eta = sqrt( pow( s( T_p ) - pow(m_pi, 2) - 4.*pow(m_p, 2), 2 ) - 16.*pow(m_pi, 2)*pow(m_p, 2) )/( 2. *  m_pi * sqrt( s( T_p ) ) );
    double gam_s = sqrt( pow(M_res,2) * (pow(M_res,2) + pow(Gam_res,2)) );
    double K = (sqrt(8.) * M_res * Gam_res * gam_s)/( M_PI * sqrt( pow(M_res,2) + gam_s ) );
    double f_BW = (m_p * K)/(pow( pow( sqrt( s( T_p ) ) - m_p, 2) - pow(M_res,2), 2 ) + pow(M_res,2) * pow(Gam_res,2));
    return sigma_0 * pow(eta,1.95) * ( 1. + eta + pow(eta,5) ) * pow(f_BW,1.86);
    }

  double sigma_2pi_mb( double T_p ){return 5.7 /( 1. + exp( -9.3 * ( T_p - 1.4) ) );}

  double sigma_inel_mb( double T_p ){
    if (T_p >= T_p_th){
      return ( 30.7 - 0.96 * log( T_p/T_p_th ) + 0.18 * pow(log( T_p/T_p_th ),2) ) * pow( (1. - pow( (T_p_th/T_p), 1.9 ) ), 3 );
      }
    else {return 0.;}
    }

  double n_pi_low( double T_p ){
    double Q_p = (T_p - T_p_th)/m_p;
    return -6e-3 + 0.237 * Q_p - 0.023 * pow( Q_p, 2 );
    }

  double n_pi_high( double T_p ){
    double out( double T_p, double a1, double a2, double a3, double a4, double a5 ){
      double xi_p = (T_p - 3.)/m_p;
      return a1 * pow( xi_p, a4 ) * ( 1. + exp( -a2 * pow( xi_p, a5 ) ) ) * ( 1. - exp( -a3 * pow( xi_p, 0.25 ) ) );
      }
    if (T_p >= 5.){
      double a1 = 0.728;
      double a2 = 0.596;
      double a3 = 0.491;
      double a4 = 0.2503;
      double a5 = 0.117;
      return out( T_p, a1, a2, a3, a4, a5 );
      }
    else if (T_p > 50.){
      double a1 = 0.652;
      double a2 = 0.0016;
      double a3 = 0.488;
      double a4 = 0.1928;
      double a5 = 0.483;
      return out( T_p, a1, a2, a3, a4, a5 );
      }
    else if (T_p > 100.){
      double a1 = 5.436;
      double a2 = 0.254;
      double a3 = 0.072;
      double a4 = 0.075;
      double a5 = 0.166;
      return out( T_p, a1, a2, a3, a4, a5 );
      }
    else{
      return 0.;
      }
    }

  double sigma_pi_mb( double T_p ){
   if (T_p < T_p_th){return 0.;}
    else if (T_p < 2.){return sigma_1pi_mb( T_p ) + sigma_2pi_mb( T_p );}
    else if (T_p < 5.){return sigma_inel_mb( T_p ) * n_pi_low( T_p );}
    else {return sigma_inel_mb( T_p ) * n_pi_high( T_p );}
    }

  double G( double T_p ){
    double T_p_0 = 1e3;
    return 1. + log(fmax( 1., sigma_inel_mb( T_p )/sigma_inel_mb( T_p_0 ) ));
    }

  double eps( double T_p ){
    double eps_c = 1.37;
    double eps_1 = 0.29;
    double eps_2 = 0.1;
    double sigma_pp_R_mb = 31.4;
    return eps_c + (eps_1 + eps_2) * (sigma_pp_R_mb * G( T_p ))/sigma_inel_mb( T_p );
    }

  double Amax;

  if (T_p < T_p_th){Amax = 0.;}
  if (T_p < 1.){Amax = b0 * sigma_pi_mb( T_p )/E_pi_LAB_max( T_p ) * mb_cm2;}
  else {Amax = b1( T_p ) * pow(theta_p, -b2( T_p )) * exp( b3( T_p ) * pow(log(theta_p), 2) ) * sigma_pi_mb( T_p )/m_p * mb_cm2;}

  return eps( T_p ) * Amax;
  }






double dsig_dEg( double T_p, double E_gam ){
  double gam_pi_LAB =  E_pi_LAB_max( T_p )/m_pi;
  double beta_pi_LAB = sqrt( 1. - pow(gam_pi_LAB, -2) );
  double E_gam_max = m_pi/2. * gam_pi_LAB * ( 1. + beta_pi_LAB );
  double E_gam_min = m_pi/2. * gam_pi_LAB * ( 1. - beta_pi_LAB );

  double Y_gam = E_gam + pow(m_pi, 2)/(4. * E_gam);
  double Y_gam_max = E_gam_max + pow(m_pi, 2)/(4. * E_gam_max);  
  double X_gam = ( Y_gam - m_pi )/( Y_gam_max - m_pi );

  double alpha( double T_p ){if ( T_p < T_p_th ){return 0.;} if ( T_p <= 20.){return 1.0;} else {return 0.5;}}

  double kappa( double T_p ){return 3.29 - 0.2 * pow(T_p/m_p, -1.5);}

  double mu( double T_p ){double q = (T_p - 1.)/m_p; return 5./4. * pow(q, 5./4.) * exp( -5./4. * q );}

  double beta( double T_p ){
    if ( T_p < T_p_th ){return 0.;}
    if ( T_p <= 1.){return kappa( T_p );}
    else if ( T_p <= 4.){return mu( T_p ) + 2.45;}
    else if ( T_p <= 20.){return 1.5 * mu( T_p ) + 4.95;}
    else if ( T_p <= 100.){return 4.2;}
    else {return 4.9;}
    }

  double gamma( double T_p ){
    if ( T_p <= 1.){return 0.;}
    else if ( T_p <= 4.){return mu( T_p ) + 1.45;}
    else if ( T_p <= 20.){return mu( T_p ) + 1.5;}
    else {return 1.;}
    }

  double lambda = 3.;
  double C = lambda * m_pi/Y_gam_max;

  double F( double T_p, double E_gam ){
    if (pow(X_gam, alpha( T_p )) <= 1.){return pow( (1. - pow(X_gam, alpha( T_p ))) , beta( T_p ) ) / pow( (1. + X_gam/C) , gamma( T_p ));}
    else{return 0.;}
    }

  return A_max( T_p ) * F( T_p, E_gam );

  }

//injected number density per energy and time [ L^-3 E^-1 T^-1 ]
// J is dndE * beta_p * c
double J( double T_p, double C ){
  double theta_p = T_p/m_p;
  return C * (p-1.)/m_p * pow( (theta_p + 1.) , -p ) * c_cmsm1;
//  return C * (p-1.)/m_p * pow( sqrt( pow(T_p,2) + 2.*m_p*T_p )/m_p , -p ) * c_cmsm1 * sqrt( 1.-pow(m_p,2)/(pow(T_p,2) + 2.*m_p*T_p + pow(m_p,2)) ) ;
  }

//double n_gam_tot( double C ){
//  return C * ( pow( (T_p_low/m_p + 1.) , 1.-p ) - pow( (T_p_high/m_p + 1.) , 1.-p ) );
//  }

//double EdndE( double T_p, double C ){
//  return C * ( pow( (T_p_low/m_p + 1.) , 1.-p );
//  }

//double tau_gg( double E_gam, double h_pc, double Phi_gam, double C ){
//  double c_pcsm1 = c_cmsm1/pc_cm;
//  return ( h_pc * J( E_gam, C ) )/( 50. * Phi_gam * E_gam ); //* c_pcsm1 );
//  }



//double J_star( double T_p, double D, double f_cal_E ){
//  double P = sqrt( pow(T_p,2) + 2. * T_p * m_p );
//  return 0.4/pow(P,1.7) * ( (T_p + m_p)/P ) * sqrt( 1. - pow(1. + T_p/m_p, -2) ) * f_cal_E;
//  }

double Phi( double E_gam, double n_H, double C,  gsl_spline spline, gsl_interp_accel acc ){
  double result;
  double abserr;
  double f(double T_p, void* p){
    return dsig_dEg( T_p, E_gam ) * J( T_p, C ) * gsl_spline_eval( &spline, T_p + m_p, &acc );
    }
  gsl_function F;
  F.function = &f;

  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);

  gsl_integration_qag( &F , T_p_low, T_p_high, 0., 1e-4, 100000, GSL_INTEG_GAUSS61, w, &result, &abserr );

  gsl_integration_workspace_free(w);

  return 4. * M_PI * n_H * result;
  }


//Breith-Wheeler x-section
double sigma_gg( double E1, double E2 ){
  if ( E1 * E2 >= pow(m_e, 2) ){
    double beta_hat = sqrt( 1. - pow(m_e, 2)/( E1 * E2 ) );
    return M_PI/2. * pow( (alpha*hbar_GeVs*c_cmsm1)/m_e, 2 ) * (1.-pow(beta_hat,2)) * (2.*beta_hat * (pow(beta_hat,2)-2.) + 
           (3.-pow(beta_hat,4))*log((1.+beta_hat)/(1.-beta_hat)));
    }
  else { return 0.; } 
    }

double T_planck_gal_K( double Sigma_star_Msolpcm2 ){
  return pow( 102.731 * Sigma_star_Msolpcm2/(2.*arad_GeVcmm3Km4*c_cmsm1*pow(pc_cm,2)) * Lsol_GeVsm1, 0.25);
  }

//number density per energy bin Planck
double nE_gam_Planck( double E_gam_GeV, double Sigma_star_Msolpcm2 ){
  return  (8.*M_PI*pow(E_gam_GeV,2))/pow(hbar_GeVs * 2.*M_PI * c_cmsm1, 3) * 1./(exp(E_gam_GeV/(kb_GeVKm1 * T_planck_gal_K( Sigma_star_Msolpcm2 )))-1.);
  }

double tau_pair( double E_gam_GeV, double h_pc, double Sigma_star_Msolpcm2 ){
  double result;
  double abserr;
  double f(double E, void* p){
//    return nE_gam_Planck( E, Sigma_star_Msolpcm2 );
    return sigma_gg( E, E_gam_GeV ) * nE_gam_Planck( E, Sigma_star_Msolpcm2 );
    }
  gsl_function F;
  F.function = &f;

  double E_low = pow(m_e, 2)/E_gam_GeV;
  double E_high = fmax(10. * kb_GeVKm1 * T_planck_gal_K( Sigma_star_Msolpcm2 ), E_low); // 1e-6;

//  gsl_integration_qag( &F , E_low, E_high, 0., 1e-10, 100000, GSL_INTEG_GAUSS61, w, &result, &abserr );

//  return result * h_pc * pc_cm;

  return 4.123 * h_pc/100. *(arad_GeVcmm3Km4 * pow(T_planck_gal_K( Sigma_star_Msolpcm2 ),4))/1e-6 * 1e-11/E_low;
  }









#endif
