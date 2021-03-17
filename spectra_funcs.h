#ifndef spectra_funcs_h
#define spectra_funcs_h

#include <stdio.h>
#include <math.h>
#include <gsl_interp2d.h>
#include <gsl_integration.h>
#include <gsl_spline.h>
#include <cubature.h>
#include "share/math_funcs.h"
#include "cosmo_funcs.h"

//extern gsl_integration_workspace * w;
//#pragma omp threadprivate(w)
//gsl_integration_workspace * w;
extern double q;
extern double c_cmsm1;
extern double pc_cm;
extern double E_0;

extern double Emin_GeV;
extern double Emax_GeV;

extern double Emin_e_GeV;
extern double Emax_e_GeV;

extern double Delta_x[11];
extern double Phi_1H[11];
extern double Phi_2H[11];


/* T_min T_max - minimum and maximum CR energy to integrate over */

double E_cutoff = 1e8;
double E_cutoff_e = 1e5;

double T_p_low = 1.e-3; //GeV
double T_p_high = 1.e8; //1e6 //GeV

double T_p_low_norm = 0.; //0.9383; //GeV
double T_p_high_norm = 1e8; //GeV


double e__esu = 4.80320451e-10; //statC or esu

double m_p = 0.93827208816; //GeV
double m_pi0 = 0.1349766; //GeV
double m_piC = 0.13957018; //GeV
double m_mu = 0.1056583755; //GeV
double m_e = 0.51099895000e-3; //GeV
double m_e__GeV = 0.51099895000e-3; //GeV
#define T_p_th (2. * m_pi0 + pow(m_pi0, 2)/(2. * m_p)) /* GeV */ //0.2797
double alpha = 1./137.;
double hbar_GeVs = 6.582e-25;
double h_GeVs = 4.135668e-24;
double h__Js = 6.62607015e-34;
double kb_GeVKm1 = 8.617e-14;
double arad_GeVcmm3Km4 = 4.722e-12; //GeV cm^-3 K^-4
double Lsol_GeVsm1 = 2.3882e36; //GeV s^-1
double sigma_T_mb = 665.2487158; //mbarn
double mb_cm2 = 1e-27; //cm^2
double GeV_erg = 1.6021773e-3; //erg

//Bremsstrahlung Peretti 2019 following Stecker 2971
double sigma_brems_mb = 34.;

double K_pi = 0.17;

struct Phi_out {
  double Phi;
  double Phi_fcal1;
  };

struct gsl_spline_obj {
  gsl_spline *spline;
  gsl_interp_accel *acc;
  };

struct gsl_spline_obj gsl_so( double *xdata, double *ydata, size_t n_data ){
  struct gsl_spline_obj gsl_so;
  gsl_so.acc = gsl_interp_accel_alloc();
  gsl_so.spline = gsl_spline_alloc(gsl_interp_linear, n_data);
  gsl_spline_init( gsl_so.spline, xdata, ydata, n_data);
  return gsl_so;
  }

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

struct hcub_data {
  double E_gam;
  double B_G;
  double T_dust__K;
  double T_optBB__K;
  int n_radcomp;
  double *E_peaks__GeV;
  double *urads;
  double urad;
  double E_peak__GeV;
  double uradFIR;
  double uradopt;
  double uradtot;

  double z;
  double h;
  double Sigstar;
  double C;
  double C_e;
  double n_H;
  double E_cut;
  double E_cut_e;
  struct gsl_spline_obj gsl_so;
  gsl_spline spline;
  gsl_interp_accel acc;
  };






/* All energies in GeV */

/* Yu et al. fit to Manga data */
double sigma_gas_Yu( double SFR_Msolyrm1 ){return 32.063 * pow( SFR_Msolyrm1, 0.096 );}
/* Inverse KS derived from PHIBSS2 data */
//double Sigma_gas_Msolpcm2_iKS( double Sigma_SFR_Msolyrm1pcm2 ){return 2.511e8 * pow( Sigma_SFR_Msolyrm1pcm2, 0.936 );}
/* Bigiel */
//double Sigma_gas_Msolpcm2_iKS( double Sigma_SFR_Msolyrm1pcm2 ){return pow( Sigma_SFR_Msolyrm1pcm2/7.94e-3 * 1e6, 1 ) * 10.;}
/* proper KS */
double Sigma_gas_Msolpcm2_iKS( double Sigma_SFR_Msolyrm1pcm2 ){return pow( Sigma_SFR_Msolyrm1pcm2/2.5e-4 * 1e6, 1/1.4 );}

double Sigma_gas_Shi_Msolpcm2_iKS( double Sigma_SFR_Msolyrm1pcm2, double Sigma_star_Msolpcm2 ){
  return pow( 10., 10.28 ) * Sigma_SFR_Msolyrm1pcm2 * pow( Sigma_star_Msolpcm2, -0.48 ) ;
  }


double s( double T_p ){return 2. * m_p * ( T_p + 2. * m_p );}
double gam_CM( double T_p ){return ( T_p + 2. * m_p )/( sqrt( s( T_p ) ) );}
double E_pi_CM( double T_p ){return ( s( T_p ) - 4. * pow(m_p, 2) + pow(m_pi0, 2) )/( 2. * sqrt( s( T_p ) ) );}
double P_pi_CM( double T_p ){return sqrt( pow(E_pi_CM( T_p ), 2) - pow(m_pi0, 2) );}
double beta_CM( double T_p ){return sqrt( 1. - pow(gam_CM( T_p ), -2) );}
double E_pi_LAB_max( double T_p ){return gam_CM( T_p ) * ( E_pi_CM( T_p ) + P_pi_CM( T_p ) * beta_CM( T_p ) );}

double sigma_tot_mb( double T_p ){
  double sig_inel = 0.;
  double sig_el = 0.;

  if (T_p >= T_p_th){
    sig_inel = ( 30.7 - 0.96 * log( T_p/T_p_th ) + 0.18 * pow(log( T_p/T_p_th ),2) ) * pow( (1. - pow( (T_p_th/T_p), 1.9 ) ), 3 );}
  else {sig_inel = 0.;}

  double p = sqrt( pow(T_p,2) + 2. * T_p * m_p );
  if (p < 0.8){ sig_el =  23.5 + 1000. * pow( p - 0.7 ,4) - 23.6 + 1250./(0.8+50.) - 4. * pow( 0.8 - 1.3, 2); }
  else if (p <= 2.){ sig_el = 1250./(p+50.) - 4. * pow( p - 1.3, 2); }
  else if (p > 2.){ sig_el = 77./(p+1.5) - 22. + 1250./(2.+50.) - 4. * pow( 2. - 1.3, 2); }

  return sig_inel + sig_el;
  }


double sigma_inel_mb( double T_p ){
  if (T_p >= T_p_th){
    return ( 30.7 - 0.96 * log( T_p/T_p_th ) + 0.18 * pow(log( T_p/T_p_th ),2) ) * pow( (1. - pow( (T_p_th/T_p), 1.9 ) ), 3 );
    }
  else {return 0.;}
  }


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
    double eta = sqrt( pow( s( T_p ) - pow(m_pi0, 2) - 4.*pow(m_p, 2), 2 ) - 16.*pow(m_pi0, 2)*pow(m_p, 2) )/( 2. *  m_pi0 * sqrt( s( T_p ) ) );
    double gam_s = sqrt( pow(M_res,2) * (pow(M_res,2) + pow(Gam_res,2)) );
    double K = (sqrt(8.) * M_res * Gam_res * gam_s)/( M_PI * sqrt( pow(M_res,2) + gam_s ) );
    double f_BW = (m_p * K)/(pow( pow( sqrt( s( T_p ) ) - m_p, 2) - pow(M_res,2), 2 ) + pow(M_res,2) * pow(Gam_res,2));
    return sigma_0 * pow(eta,1.95) * ( 1. + eta + pow(eta,5) ) * pow(f_BW,1.86);
    }

  double sigma_2pi_mb( double T_p ){return 5.7 /( 1. + exp( -9.3 * ( T_p - 1.4) ) );}

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

  if ( T_p > m_pi0/K_pi ){return eps( T_p ) * Amax;}
  else{return Amax;}
  }

double dsig_dEg( double T_p, double E_gam ){
  double gam_pi0_LAB =  E_pi_LAB_max( T_p )/m_pi0;
  double beta_pi_LAB = sqrt( 1. - pow(gam_pi0_LAB, -2) );
  double E_gam_max = m_pi0/2. * gam_pi0_LAB * ( 1. + beta_pi_LAB );
  double E_gam_min = m_pi0/2. * gam_pi0_LAB * ( 1. - beta_pi_LAB );

  double Y_gam = E_gam + pow(m_pi0, 2)/(4. * E_gam);
  double Y_gam_max = E_gam_max + pow(m_pi0, 2)/(4. * E_gam_max);  
  double X_gam = ( Y_gam - m_pi0 )/( Y_gam_max - m_pi0 );

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
  double C = lambda * m_pi0/Y_gam_max;

  double F( double T_p, double E_gam ){
    if (pow(X_gam, alpha( T_p )) <= 1.){
      if ( T_p <= 1.){return pow( (1. - X_gam) , beta( T_p ) );}
      else {return pow( (1. - pow(X_gam, alpha( T_p ))) , beta( T_p ) ) / pow( (1. + X_gam/C) , gamma( T_p ));}
      }
    else{return 0.;}
    }

  if (T_p < T_p_th){return 0.;}

  return A_max( T_p ) * F( T_p, E_gam );

  }


double J( double T_p, double C, double E_cut ){
  double E_norm = m_p;
  double beta = 1.0;
  double E = T_p + m_p;
  double p = sqrt( pow(T_p,2) + 2.*m_p*T_p );
//  return C * (p-1.)/E_norm * pow( E/E_norm , -p ) * exp(-pow(E/E_cut, beta));
//  return C * (p-1.)/E_norm * (T_p+m_p)/sqrt( pow(T_p,2) + 2.*m_p*T_p ) * pow( sqrt( pow(T_p,2) + 2.*m_p*T_p )/E_norm , -p );
  //power law in momentum
  return C * pow( p, -q ) * p/E * exp(-pow(E/E_cut, beta));
  }

double J_e( double T_e, double C_e, double E_cut_e ){
  double E_norm = m_e__GeV;
  double beta = 1.0;
  double E_e = T_e + m_e__GeV;
  double p_e = sqrt( pow(T_e,2) + 2.*m_e*T_e );
//  return C_e * (p-1.)/E_norm * pow( E_e/E_norm , -p ) * exp(-pow(E_e/E_cut_e, beta));
  return C_e * pow( p_e, -q ) * E_e/p_e * exp(-pow(E_e/E_cut_e, beta));
  }



double C_norm_E( double E_cut ){
  double result;
  double abserr;

  int f( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      fval[j] = x[j*ndim+0] * J( x[j*ndim+0], 1., fdata_in.E_cut );
      }
    return 0;
    }

  double xmin[1] = { T_p_low_norm };
  double xmax[1] = { T_p_high_norm };

  struct hcub_data fdata;
  fdata.E_cut = E_cut;
  struct hcub_data *fdata_ptr = &fdata;

  hcubature_v( 1, f, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &result, &abserr );

  return result;
  }

double C_norm_E_e( double E_cut_e ){
  double result;
  double abserr;

  int f( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      fval[j] = x[j*ndim+0] * J_e( x[j*ndim+0], 1., fdata_in.E_cut );
      }
    return 0;
    }

  double xmin[1] = { 0. };
  double xmax[1] = { T_p_high_norm };

  struct hcub_data fdata;
  fdata.E_cut = E_cut_e;
  struct hcub_data *fdata_ptr = &fdata;

  hcubature_v( 1, f, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &result, &abserr );

  return result;
  }

double C_norm_E_e_1GeV( double E_cut_e ){
  double result;
  double abserr;

  int f( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      fval[j] = x[j*ndim+0] * J_e( x[j*ndim+0], 1., fdata_in.E_cut );
      }
    return 0;
    }

  double xmin[1] = { 1. };
  double xmax[1] = { T_p_high_norm };

  struct hcub_data fdata;
  fdata.E_cut = E_cut_e;
  struct hcub_data *fdata_ptr = &fdata;

  hcubature_v( 1, f, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &result, &abserr );

  return result;
  }

struct Phi_out Phi( double E_gam, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc ){
  struct Phi_out Phi_out;

  double result;
  double abserr;

  int f( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    double beta_p[npts], f_cal[npts];
    for (j = 0; j < npts; ++j)
      {
      beta_p[j] = sqrt( 1. - pow(m_p,2)/pow( x[j*ndim+0] + m_p, 2) );
      f_cal[j] = gsl_spline_eval( &fdata_in.spline, x[j*ndim+0], &fdata_in.acc );
      fval[j] = dsig_dEg( x[j*ndim+0], fdata_in.E_gam ) * J( x[j*ndim+0], fdata_in.C, fdata_in.E_cut ) * c_cmsm1 * beta_p[j] * f_cal[j] * fdata_in.n_H;
      }
    return 0;
    }

  int f_fcal1( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    double beta_p[npts], f_cal[npts];

    for (j = 0; j < npts; ++j)
      {
      beta_p[j] = sqrt( 1. - pow(m_p,2)/pow( x[j*ndim+0] + m_p, 2) );
      fval[j] = dsig_dEg( x[j*ndim+0], fdata_in.E_gam ) * J( x[j*ndim+0], fdata_in.C, fdata_in.E_cut ) * c_cmsm1 * beta_p[j] * fdata_in.n_H;
      }
    return 0;
    }


  double xmin[1] = { T_p_low };
  double xmax[1] = { T_p_high };

  struct hcub_data fdata;
  fdata.E_gam = E_gam;
  fdata.n_H = n_H;
  fdata.C = C;
  fdata.E_cut = E_cut;
  fdata.spline = spline;
  fdata.acc = acc;

  struct hcub_data *fdata_ptr = &fdata;


  hcubature_v( 1, f, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &result, &abserr );
  Phi_out.Phi = result;
  hcubature_v( 1, f_fcal1, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &result, &abserr );
  Phi_out.Phi_fcal1 = result;

  return Phi_out;
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
//  return pow( 102.731 * Sigma_star_Msolpcm2/(2.*arad_GeVcmm3Km4*c_cmsm1*pow(pc_cm,2)) * Lsol_GeVsm1, 0.25);
  return pow( 3673.35 * Sigma_star_Msolpcm2  * Lsol_GeVsm1 /(2.*arad_GeVcmm3Km4*c_cmsm1*pow(pc_cm,2)), 0.25);
  }

//number density per energy bin Planck
double nE_gam_Planck( double E_gam_GeV, double Sigma_star_Msolpcm2 ){
  return  (8.*M_PI*pow(E_gam_GeV,2))/pow(hbar_GeVs * 2.*M_PI * c_cmsm1, 3) * 1./(exp(E_gam_GeV/(kb_GeVKm1 * T_planck_gal_K( Sigma_star_Msolpcm2 )))-1.);
  }

//Razzaque, Meszaros, Zhang
double u_CIMB( double T_K ){return 6.49394 * pow(kb_GeVKm1 * T_K,4)/( pow(M_PI,2) * pow(hbar_GeVs * c_cmsm1,3) );} 
double tau_pair( double E_gam_GeV, double z, double h_pc, double Sigma_star_Msolpcm2 )
  {
  double EIR_low = pow(m_e, 2)/E_gam_GeV;

  double E_gam_min = pow(m_e,2)/(2.*E_gam_GeV);
  double res_out = fmax( u_CIMB(T_planck_gal_K( Sigma_star_Msolpcm2 ))/E_gam_min, u_CIMB(2.725 * (1+z))/E_gam_min ) *
         3./8. * sigma_T_mb * mb_cm2 * h_pc * pc_cm;
  return res_out;
  }

double dnde_gam_Planck( double E_gam_GeV, double T_K )
  {
  return  (8.*M_PI*pow(E_gam_GeV,2))/pow(h_GeVs * c_cmsm1, 3) * 1./( exp( E_gam_GeV/(kb_GeVKm1 * T_K) ) - 1. );
  }


//Magnelli 2018: Tdust - sSFR relation
double Tdust_K( double z, double SFR_Msolyrm1, double M_star_Msol )
  {
  return 98.*pow(1.+z,-0.065) + 6.9*log10(SFR_Msolyrm1/M_star_Msol);
  }
//Ackermann 2021 (Kennicutt) FIR-SFR
double LFIR_Lsol( double SFR_Msolyrm1 )
  {
  return SFR_Msolyrm1/( 1.7e-10 * 0.79 );
  }

//Persic 2008
double uFIR_GeVsm1cmm3( double L_FIR_Lsol, double R_e_kpc, double h_pc )
  {
  return L_FIR_Lsol*Lsol_GeVsm1/( c_cmsm1 * (4.*M_PI*pow(R_e_kpc*1e3*pc_cm,2)+2.*sqrt(2.)*R_e_kpc*1e3*pc_cm*2.*h_pc*pc_cm) );
  }

double E0_GeV(void){return 2.e12*h_GeVs;}
double Cdil_dBB( double u_FIR_GeVsm1cmm3, double E0_GeV, double T_dust_K)
  {
  return u_FIR_GeVsm1cmm3/24.8863 * E0_GeV*pow(h_GeVs*c_cmsm1,3)/(8.*M_PI*pow(kb_GeVKm1*T_dust_K,5));
  }

double Cdil_BB( double u_rad_GeVsm1cmm3, double T_dust_K)
  {
  return u_rad_GeVsm1cmm3/(arad_GeVcmm3Km4*pow(kb_GeVKm1*T_dust_K,4)); //???? Units are fucked
  }


//Obtained from a fit to DustPedia galaxy data excluding E and S0 type galaxies
//Fraction of luminosity absorbed and reemitted by dust
double fabsdust( double sSFR__yrm1, double Mstar__Msol )
  {
  double res = pow( 10., 0.51175779 * log10(sSFR__yrm1) + 0.30915893 * log10(Mstar__Msol) + 1.45388753 );
  return fmax( 0., fmin( 1., res ) );
  }
//bolometric luminosity fitted excluding E and S0 galaxies in DustPedia
double Lbol__Lsol( double M_star__Msol, double SFR__Msolyrm1 )
  {
  return pow( 10., 0.67450795 * log10(M_star__Msol) + 0.26664641 * log10(SFR__Msolyrm1) + 3.53945016 );
  }
double Labs__Lsol( double SFR__Msolyrm1 )
  {
  return fmax( pow( 10., 0.9617451 * log10(SFR__Msolyrm1) + 9.723582 ), LFIR_Lsol( SFR__Msolyrm1 ) );
  }

double Lunatt__Lsol( double M_star__Msol )
  {
  return pow( 10., 0.8214778 * log10(M_star__Msol) + 1.899166 );
  }

double urad__GeVcmm3( double L__Lsol, double R_e__kpc, double h__pc )
  {
  return L__Lsol*Lsol_GeVsm1/( c_cmsm1 * 2.*M_PI*pow(R_e__kpc*1e3*pc_cm,2) );
//  return L__Lsol*Lsol_GeVsm1/( 2.*M_PI*pow(R_e__kpc*1e3*pc_cm,2) * 2*h__pc*pc_cm ) * h__pc*pc_cm/c_cmsm1;
  }

double dndE_dBB_GeVm1cmm3( double E_gam_GeV, double E0_GeV, double C_dil, double T_dust_K )
  {
  return C_dil*8.*M_PI*pow(E_gam_GeV,2)/pow(h_GeVs*c_cmsm1,3)*(E_gam_GeV/E0_GeV)/(exp(E_gam_GeV/(kb_GeVKm1*T_dust_K))-1.);
  }

double dndE_BB_GeVm1cmm3( double E_gam_GeV, double C_dil, double T_star_K )
  {
  return C_dil*8.*M_PI*pow(E_gam_GeV,2)/pow(h_GeVs*c_cmsm1,3)/(exp(E_gam_GeV/(kb_GeVKm1*T_star_K))-1.);
  }

double tau_gg_gal( double E_gam_GeV, double z, double SFR_Msolyrm1, double M_star_Msol, double r_e_kpc, double h_pc )
  {
  //Magnelli 2018: Tdust - sSFR relation
  double T_dust_K = Tdust_K( z, SFR_Msolyrm1, M_star_Msol );
  //Ackermann 2021 (Kennicutt) FIR-SFR
  double L_FIR_Lsol = LFIR_Lsol( SFR_Msolyrm1 );

  //Persic 2008
//  double u_FIR_GeVsm1cmm3 = uFIR_GeVsm1cmm3( L_FIR_Lsol, R_e_kpc, h_pc );


  double f_abs = fabsdust( SFR_Msolyrm1/M_star_Msol, M_star_Msol );

  double L_bol__Lsol = Lbol__Lsol( M_star_Msol, SFR_Msolyrm1 );

//  double u_rad__GeVcmm3 = urad__GeVcmm3( L_bol__Lsol, R_e_kpc, h_pc );
//  double u_FIR_GeVcmm3 = f_abs * u_rad__GeVcmm3;
//  double u_BB_GeVcmm3 = (1.-f_abs) * u_rad__GeVcmm3;

  double u_FIR__GeVcmm3 = urad__GeVcmm3( Labs__Lsol( SFR_Msolyrm1 ), r_e_kpc, h_pc );
  double u_opt__GeVcmm3 = urad__GeVcmm3( Lunatt__Lsol( M_star_Msol ), r_e_kpc, h_pc );

  double C_dil_dBB = Cdil_dBB( u_FIR__GeVcmm3, E0_GeV(), T_dust_K );

  //make some assumption of average stellar surface temperature for BB contribution.
//  double T_star__K = 6000.; 

//  double C_dil_BB = Cdil_BB( u_BB_GeVsm1cmm3, T_star__K );

  double E_gam_field = pow(m_e,2)/(2.*E_gam_GeV);
//  double res_out = fmax( dndE_dBB_GeVm1cmm3( E_gam_field, E0_GeV(), C_dil, T_dust_K ), dnde_gam_Planck( E_gam_field, 2.725 * (1+z) ) )*E_gam_field * 3./8. * sigma_T_mb * mb_cm2 * h_pc * pc_cm;

  //modify this as it should be the sum of the IR+CMB contributed number densities rather than just the max of them.
  double res_out = (dndE_dBB_GeVm1cmm3( E_gam_field, E0_GeV(), C_dil_dBB, T_dust_K ) + dnde_gam_Planck( E_gam_field, 2.725 * (1+z) ) )*E_gam_field * 3./8. * sigma_T_mb * mb_cm2 * h_pc * pc_cm;

//printf( "%s %e %e %e %e \n", "tau_gg", res_out, T_dust_K, dndE_dBB_GeVm1cmm3( E_gam_field, E0_GeV, C_dil, T_dust_K ), dnde_gam_Planck( E_gam_field, 2.725 * (1+z) ) );
  return res_out;
  }


//Neutrinos
//Kelner et al. (2009)
//Eq (62)
double F_e( double x, double E_p_TeV ){
  double L = log( E_p_TeV );
  double B_e = 1./( 69.5 + 2.65 * L + 0.3 * pow( L, 2 ) );
  double beta_e = 1./pow( 0.201 + 0.062 * L + 0.00042 * pow( L, 2 ) , 0.25 );
  double k_e = ( 0.279 + 0.141 * L + 0.0172 * pow( L, 2 ) )/( 0.3 + pow( 2.3 + L, 2 ) );
  return B_e * pow( 1. + k_e * pow( log(x), 2 ) , 3 )/( x * ( 1. + 0.3/pow( x, beta_e ) ) ) * pow( -log(x), 5 );
  }

//Eq (66)
double F_numu( double x, double E_p_TeV ){
  double y = x/0.427;
  double L = log( E_p_TeV );
  double B_nu = 1.75 + 0.204 * L + 0.010 * pow( L, 2 );
  double beta_nu = 1./( 1.67 + 0.111 * L + 0.0038 * pow( L, 2 ) );
  double k_nu = 1.07 - 0.086 * L + 0.002 * pow( L, 2 );
  return B_nu * log(y)/y * pow( ( 1.-pow(y, beta_nu) )/( 1. + k_nu * pow(y, beta_nu) * (1.-pow(y, beta_nu)) ) , 4)  * 
         ( 1./log(y) - (4.*beta_nu*pow(y, beta_nu))/(1.-pow(y, beta_nu)) - 
         (4. * k_nu * beta_nu * pow(y, beta_nu)*(1. - 2. * pow(y, beta_nu)))/(1. + k_nu * pow(y, beta_nu)*(1.-pow(y, beta_nu))));
  }

//Eq.(71)
double Phi_nu( double E_nu_GeV, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc ){
  double result;
  double abserr;
  double F_nu(double T_p, void* p){
    double E_p_GeV = T_p + m_p;
    double beta_p = sqrt( 1. - pow(m_p,2)/pow(T_p+m_p,2) );
    double x = E_nu_GeV/E_p_GeV;
    double f_cal = gsl_spline_eval( &spline, T_p, &acc );
    double Fmu;
    double Fe;
    if (x < 0.430) {
      Fmu = F_numu( x, E_p_GeV/1e3 );
      Fe = F_e( x, E_p_GeV/1e3 );
      }
    else if (x < 0.99) {
      Fmu = 0.;
      Fe = F_e( x, E_p_GeV/1e3 );
      }
    else {
      Fmu = 0.;
      Fe = 0.;
      }
    double L = log( E_p_GeV/1e3 );
    double sig_inel = sigma_inel_mb( T_p ); //34.3 + 1.88 * L + 0.25 * pow(L,2);
    return n_H * c_cmsm1 * beta_p * J( T_p, C, E_cut ) * sig_inel/E_p_GeV * mb_cm2 * ( Fmu + 2. * Fe ) * f_cal * 2.;
    }
  gsl_function F;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);
  F.function = &F_nu;
  double T_p_high_nu = T_p_high;
  double T_p_low_nu = fmax( E_nu_GeV - m_p, 1e2); // fmax(E_nu_GeV - m_p, 36.3 - m_p);
  if (T_p_low_nu < T_p_high_nu){
    gsl_integration_qag( &F , T_p_low_nu, T_p_high_nu, 0., 1e-4, 100000, GSL_INTEG_GAUSS61, w, &result, &abserr );
    } 
  else {
    result = 0.;
    }
  gsl_integration_workspace_free(w);
  return result;
  }

//Peretti 2019 summarised method with separata Pion production

double Heaviside( double arg ){
  if (arg < 0){return 0.;}
  else {return 1.;}
  }

//Peretti 2019
double q_pi( double E_pi_GeV, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc ){
  double T_p = E_pi_GeV/K_pi;
  double beta_p = sqrt( 1. - pow(m_p,2)/pow(T_p+m_p,2) );
  double f_cal = gsl_spline_eval( &spline, T_p, &acc );
  return n_H * c_cmsm1 * beta_p * J( T_p, C, E_cut )/K_pi * sigma_inel_mb( T_p ) * mb_cm2 * f_cal;
  }
//double q_pi_fcal1( double E_pi_GeV, double n_H, double C, gsl_spline spline, gsl_interp_accel acc ){
//  double T_p = E_pi_GeV/K_pi - m_p;
//  double beta_p = sqrt( 1. - pow(m_p,2)/pow(T_p+m_p,2) );
//  return n_H * c_cmsm1 * beta_p * J( T_p, C )/K_pi * sigma_inel_mb( T_p ) * mb_cm2;
//  }

double q_pi_norm( double E_pi_GeV, double E_cut ){
  double T_p = E_pi_GeV/K_pi - m_p;
  double beta_p = sqrt( 1. - pow(m_p,2)/pow(T_p+m_p,2) );
  return c_cmsm1 * beta_p * J( T_p, 1., E_cut )/K_pi * sigma_inel_mb( T_p ) * mb_cm2;
  }


//Peretti 2019 and Kelner 2006
double f_nu_mu2( double x ) //This is also f_e
  {
  double g_nu_mu( double x )
    {
    double lambda = 1. - pow(m_mu/m_piC,2);
    double r = 1. - lambda;
    return (9. * pow(x,2) - 6. * log(x) - 4. * pow(x,3) - 5.) * (3. - 2. * r) / (9. * pow(1.-r,2));
    }
  double h_nu_mu1( double x )
    {
    double lambda = 1. - pow(m_mu/m_piC,2);
    double r = 1. - lambda;
    return (9. * pow(r,2) - 6. * log(r) - 4. * pow(r,3) - 5.) * (3. - 2. * r) / (9. * pow(1.-r,2));
    }
  double h_nu_mu2( double x )
    {
    double lambda = 1. - pow(m_mu/m_piC,2);
    double r = 1. - lambda;
    return (9. * (r + x) - 4. * (pow(r,2) + r*x + pow(x,2))) * (1. + 2. * r) * (r - x) / (9. * pow(r,2));
    }
  double lambda = 1. - pow(m_mu/m_piC,2);
  double r = 1. - lambda;
  return g_nu_mu(x) * Heaviside(x-r) + (h_nu_mu1(x) + h_nu_mu2(x)) * Heaviside(r-x);
  }


double f_nu_e( double x )
  {
  double g_nu_e( double x )
    {
    //There is an erratum for this eqn in the paper.
    double lambda = 1. - pow(m_mu/m_piC,2);
    double r = 1. - lambda;
    return 2. * ( (1.-x) * (6. * pow(1.-x,2) + r * (5. + 5. * x - 4. * pow(x,2))) + 6. * r * log(x))/(3. * pow(1.-r,2));
//    return 2. * (1.-x) * ( (6. * pow(1.-x,2) + r * (5. + 5. * x - 4. * pow(x,2))) + 6. * r * log(x))/(3. * pow(1.-r,2));
    }
  double h_nu_e1( double x )
    {
    double lambda = 1. - pow(m_mu/m_piC,2);
    double r = 1. - lambda;
    return 2. * ((1.-r) * (6. - 7. * r + 11. * pow(r,2) - 4. * pow(r,3)) + 6. * r * log(r))/(3. * pow(1.-r,2));
    }
  double h_nu_e2( double x )
    {
    double lambda = 1. - pow(m_mu/m_piC,2);
    double r = 1. - lambda;
    return 2. * (r - x) * (7. * pow(r,2) - 4. * pow(r,3) + 7. * x * r - 4. * x * pow(r,2) - 2. * pow(x,2) - 4. * pow(x,2) * r)/(3. * pow(r,2));
    }
  double lambda = 1. - pow(m_mu/m_piC,2);
  double r = 1. - lambda;
  return g_nu_e(x) * Heaviside(x-r) + (h_nu_e1(x) + h_nu_e2(x)) * Heaviside(r-x);
  }

/*
double q_nu( double E_nu_GeV, double n_H, double C, gsl_spline spline, gsl_interp_accel acc ){
  double res_numu2_nue, res_numu1;
  double abserr_numu2_nue, abserr_numu1;
  double x_min, x_max;
  double F_numu2_nue(double x, void* p){
    return 2. *  ( f_nu_e(x) + f_nu_mu2(x) ) * q_pi( E_nu_GeV/x , n_H, C, spline, acc )/x;
    }
  double F_numu1(double x, void* p){
    double lambda = 1. - pow(m_mu/m_piC,2);
    return 2./lambda * q_pi( E_nu_GeV/x , n_H, C, spline, acc )/x;
    }

  gsl_function F;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);

  F.function = &F_numu2_nue;
  x_min = E_nu_GeV/( (T_p_high + m_p) * K_pi);
  x_max = fmin( E_nu_GeV/( (T_p_low + m_p) * K_pi), 1. );
  if (x_min < x_max){
    gsl_integration_qag( &F , x_min, x_max, 0., 1e-3, 100000, GSL_INTEG_GAUSS61, w, &res_numu2_nue, &abserr_numu2_nue );
    }
  else {
    res_numu2_nue = 0.;
    abserr_numu2_nue = 0.;
    }


  F.function = &F_numu1;
  double lambda = 1. - pow(m_mu/m_piC,2);
  x_min = E_nu_GeV/( (T_p_high + m_p) * K_pi);
  x_max = fmin( E_nu_GeV/( (T_p_low + m_p) * K_pi), lambda );
  if (x_min < x_max){
    gsl_integration_qag( &F , x_min, x_max, 0., 1e-4, 100000, GSL_INTEG_GAUSS61, w, &res_numu1, &abserr_numu1 );
    }
  else {
    res_numu1 = 0.;
    abserr_numu1 = 0.;
    }

  gsl_integration_workspace_free(w);
  return res_numu2_nue + res_numu1;
  }
*/

//Neutrino spectrum function
double q_nu( double E_nu_GeV, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc ){
  double res_numu2_nue, res_numu1;
  double abserr_numu2_nue, abserr_numu1;
  double x_min, x_max;


  int F_numu2_nue( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      fval[j] = 2. *  ( f_nu_e(x[j*ndim+0]) + f_nu_mu2(x[j*ndim+0]) ) * q_pi( fdata_in.E_gam/x[j*ndim+0] , fdata_in.n_H, fdata_in.C, 
                fdata_in.E_cut, fdata_in.spline, fdata_in.acc )/x[j*ndim+0];
      }
    return 0;
    }

  double F_numu1(double x, void* p){
    double lambda = 1. - pow(m_mu/m_piC,2);
    return 2./lambda * q_pi( E_nu_GeV/x , n_H, C, E_cut, spline, acc )/x;
    }

  double xmin[1] = { E_nu_GeV/( T_p_high * K_pi) };
  double xmax[1] = { 1. };

  struct hcub_data fdata;
  fdata.E_gam = E_nu_GeV;
  fdata.n_H = n_H;
  fdata.C = C;
  fdata.E_cut = E_cut;
  fdata.spline = spline;
  fdata.acc = acc;

  struct hcub_data *fdata_ptr = &fdata;

  if (xmin[0] < xmax[0]){
    hcubature_v( 1, F_numu2_nue, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res_numu2_nue, &abserr_numu2_nue );
    }
  else {
    res_numu2_nue = 0.;
    abserr_numu2_nue = 0.;
    }

  gsl_function F;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);

  F.function = &F_numu1;
  double lambda = 1. - pow(m_mu/m_piC,2);
  x_min = E_nu_GeV/( T_p_high * K_pi);
  x_max = fmin( E_nu_GeV/( T_p_low * K_pi), lambda );
  if (x_min < x_max){
    gsl_integration_qag( &F , x_min, x_max, 0., 1e-4, 100000, GSL_INTEG_GAUSS61, w, &res_numu1, &abserr_numu1 );
    }
  else {
    res_numu1 = 0.;
    abserr_numu1 = 0.;
    }

  gsl_integration_workspace_free(w);
  return res_numu2_nue + res_numu1;
  }

//Electron spectrum function
double q_e( double E_e_GeV, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc ){
  double res_e = 0.;
  double abserr_e;

  int F_e( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      fval[j] = fmax( 2. * f_nu_mu2(x[j*ndim+0])  * q_pi( fdata_in.E_gam/x[j*ndim+0] , fdata_in.n_H, fdata_in.C, 
                fdata_in.E_cut, fdata_in.spline, fdata_in.acc )/x[j*ndim+0], 0.);
      }
    return 0;
    }

  double xmin[1] = { E_e_GeV/( T_p_high * K_pi) };
  double xmax[1] = { 1. };

  struct hcub_data fdata;
  fdata.E_gam = E_e_GeV;
  fdata.n_H = n_H;
  fdata.C = C;
  fdata.E_cut = E_cut;
  fdata.spline = spline;
  fdata.acc = acc;

  struct hcub_data *fdata_ptr = &fdata;

  if (xmin[0] < xmax[0]){
    hcubature_v( 1, F_e, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res_e, &abserr_e );
    }
  else {
    res_e = 0.;
    abserr_e = 0.;
    }

  return res_e;
  }





double E_dBBpeak__GeV( double T_dBB__K ){ return 2.82144*kb_GeVKm1*T_dBB__K;}
double E_BBpeak__GeV( double T_BB__K ){ return 1.59362*kb_GeVKm1*T_BB__K;}





double tauBS( double E_e__GeV, double n_H__cmm3 )
  {
  //Schlickeiser 2002
  double fac_Delta( double Delta ){
    int i;
    for (i = 0; i < 10; i++){
      if ( Delta > Delta_x[i] && Delta < Delta_x[i+1] ){
        return  (Delta - Delta_x[i])/(Delta_x[i+1] - Delta_x[i]) * (Phi_1H[i+1] - Phi_1H[i]) + Phi_1H[i];
        }
      }
    printf("Error in tauBS");
    return 1.;
    }

  double gamma = E_e__GeV/m_e__GeV;
  if (gamma < 15.){
    return 2.*M_PI/(3. * alpha * c_cmsm1*sigma_T_mb * mb_cm2 * n_H__cmm3 * 2. * (log(gamma) + log(2.) - 1./3.));
    }
  else{
    double Delta = 1./( 4. * alpha * gamma );
    return 8.*M_PI/(3.9 * alpha * c_cmsm1*sigma_T_mb * mb_cm2 * n_H__cmm3 * fac_Delta( Delta ));
    }
  }

double tausync( double E_e__GeV, double B__G )
  {
  //Ghisellini 2013
  double gamma = E_e__GeV/m_e__GeV;
  double beta = sqrt(1. - 1./pow(gamma,2));
  return (6.*M_PI * m_e__GeV*GeV_erg)/(sigma_T_mb * mb_cm2 * c_cmsm1 * pow(B__G,2) * pow(beta,2) * gamma);
  }

double tauIC( double E_e__GeV, double urad__GeVcmm3, double E_peak__GeV )
  {
  double gamma = E_e__GeV/m_e__GeV;
  double beta = sqrt(1. - 1./pow(gamma,2));
//  double gam_K = 0.27*m_e__GeV/E_peak__GeV;
//  return (3. * pow(m_e__GeV,2))/( 4. * sigma_T_mb * mb_cm2 * c_cmsm1 * urad_GeVcmm3 * pow(beta,2) * E_e );
  //Klein Nishina correction Schlickeiser 2010
//  return 3. * m_e__GeV/( 4. * sigma_T_mb * mb_cm2 * c_cmsm1 * urad_GeVcmm3 * pow(beta,2) * (E_e/m_e__GeV) ) * (pow(E_e/m_e__GeV,2) + pow(gam_K,2))/pow(gam_K,2);

  //KN according to Fang et al. 2020
  double x = 4.*gamma*E_peak__GeV/m_e__GeV;
  double Y( double x ){
    if (x <= 1.5e-3){return pow(M_PI,4)/15.;}
    else if (x < 150){return exp( -3.996e-2*pow(log(x),0) - 9.1e-1*pow(log(x),1) - 1.197e-1*pow(log(x),2) + 3.305e-3*pow(log(x),3) + 1.044e-3*pow(log(x),4) - 7.013e-5*pow(log(x),5) - 9.618e-6*pow(log(x),6));}
    else {return 3./4. * pow(M_PI/x,2) * (log(x) - 1.9805);}
    }
  return pow(M_PI,4) * m_e__GeV/( 20. * sigma_T_mb * mb_cm2 * c_cmsm1 * urad__GeVcmm3 * pow(beta,2) * gamma * Y(x) );
  }


double tauion( double E_e__GeV, double n_H__cmm3 )
  {
  //Schlickeiser 2002
  double gamma = E_e__GeV/m_e__GeV;
  double beta = sqrt(1. - 1./pow(gamma,2));
  return gamma/(9./4.*c_cmsm1*sigma_T_mb*mb_cm2*n_H__cmm3*(log(gamma)+2./3.*log(m_e__GeV*1.e9/15.0)));
  //Bethe-Bloch formula
//  return 2.*E_e*beta/(3.*c_cmsm1*sigma_T_mb*mb_cm2*m_e*n_H*(log(2*m_e__GeV * pow(beta,2)*1.e9*pow(E_e/m_e,2)/(15.0))-pow(beta,2)));
  }

double taulosse_s( double E_e__GeV, double n_H__cmm3, double B__G, int n_radcomp, double urad__GeVcmm3[], double Epeak__GeV[] )
  {
  int i;
  double tICm1 = 0.;
  for (i = 0; i < n_radcomp; i++){
    tICm1 += 1./tauIC( E_e__GeV, urad__GeVcmm3[i], Epeak__GeV[i] );
    }
  return 1./(1./tauBS( E_e__GeV, n_H__cmm3 ) + 1./tausync( E_e__GeV, B__G ) + 1./tauion( E_e__GeV, n_H__cmm3 ) + tICm1 );
  }




double eps_sync( double E_gam__GeV, double z, double B__G, double SFR__Msolyrm1, double M_star__Msol, double r_e__kpc, double h__pc, double n_H, double T_dust__K, struct gsl_spline_obj qe_so ){

  double E_larmor__GeV = h_GeVs*e__esu*B__G*c_cmsm1/(2.*M_PI*m_e__GeV*GeV_erg);
  double gamma = pow( E_gam__GeV/E_larmor__GeV, 1./2. );
  double E_syn__GeV = pow( gamma, 2) * E_larmor__GeV;
  double E_e__GeV = gamma*m_e__GeV;

  double Psyn_GeVsm1( double gamma, double B__G ){
    double beta = sqrt(1. - 1./pow(gamma,2));
    return sigma_T_mb * mb_cm2 * c_cmsm1 * pow(B__G,2) * pow(gamma*beta,2)/(6.*M_PI*GeV_erg);
    }


  if (E_e__GeV >= Emin_e_GeV && E_e__GeV <= Emax_e_GeV){

    int n_radcomp = 5;
    double Epeaks[n_radcomp];
    double urads[n_radcomp];
  
//    double u_FIR__GeVcmm3 = urad__GeVcmm3( Labs__Lsol( SFR__Msolyrm1 ), r_e__kpc, h__pc );
//    double u_opt__GeVcmm3 = urad__GeVcmm3( Lunatt__Lsol( M_star__Msol ), r_e__kpc, h__pc );

    double u_FIR__GeVcmm3 = pow(B__G,2)/(8.*M_PI*GeV_erg);//urad__GeVcmm3( Labs__Lsol( SFR__Msolyrm1 ), r_e__kpc, h__pc );
    double u_opt__GeVcmm3 = 0.; //urad__GeVcmm3( Lunatt__Lsol( M_star__Msol ), r_e__kpc, h__pc );

    double f_urad = Labs__Lsol( SFR__Msolyrm1 )/(Labs__Lsol( SFR__Msolyrm1 ) + Lunatt__Lsol( M_star__Msol ));

    //Dust diluted black body
    Epeaks[0] = E_dBBpeak__GeV( T_dust__K );
    urads[0] = u_FIR__GeVcmm3;
    Epeaks[1] = E_BBpeak__GeV( 2.725 * (1.+z) );
    urads[1] = arad_GeVcmm3Km4*pow( 2.725 * (1.+z), 4 );
    //Draine radiation fields p.121 black (grey) bodies
    Epeaks[2] = E_BBpeak__GeV( 3000. );
    urads[2] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
    Epeaks[3] = E_BBpeak__GeV( 4000. );
    urads[3] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
    Epeaks[4] = E_BBpeak__GeV( 7500. );
    urads[4] = u_opt__GeVcmm3 * f_urad;


    double tau_loss = taulosse_s( E_e__GeV, n_H, B__G, n_radcomp, urads, Epeaks );

    return Psyn_GeVsm1( gamma, B__G ) * gsl_spline_eval( qe_so.spline, E_e__GeV, qe_so.acc ) * tau_loss *
           1./2. * E_e__GeV/E_gam__GeV/E_gam__GeV;

    }
  else {return 0.;}
  }


//Stecker 1971
double eps_brems( double E_gam_GeV, double z, double B_G, double SFR__Msolyrm1, double M_star__Msol, double r_e__kpc, double h__pc, double n_H, double T_dust__K, struct gsl_spline_obj qe_so ){
  double res = 0.;
  double abserr;

  //Schlickeiser 2002 based on Gould 1969
  double sigma_BS_mb( double E_gam__GeV, double E_e__GeV ){
    int i;
    double Delta = E_gam__GeV * m_e__GeV/(4.*alpha*E_e__GeV*(E_e__GeV - E_gam__GeV));
    double Phi_1 = 0.;
    double Phi_2 = 0.;
    if (Delta > 2){
      Phi_1 = 4.*(log( 2.*E_e__GeV/m_e__GeV * (E_e__GeV - E_gam__GeV)/E_gam__GeV) - 0.5);
      Phi_2 = Phi_1;
      }
    else {
      for (i = 0; i < 10; i++){
        if ( Delta > Delta_x[i] && Delta < Delta_x[i+1] ){
          Phi_1 =  (Delta - Delta_x[i])/(Delta_x[i+1] - Delta_x[i]) * (Phi_1H[i+1] - Phi_1H[i]) + Phi_1H[i];
          Phi_2 =  (Delta - Delta_x[i])/(Delta_x[i+1] - Delta_x[i]) * (Phi_2H[i+1] - Phi_2H[i]) + Phi_2H[i];
          break;
          }
        }
      }
    return 3./(8.*M_PI) * alpha * sigma_T_mb * ((1.+pow( 1.- E_gam__GeV/E_e__GeV, 2 ) ) * Phi_1 - 2./3. * (1.- E_gam__GeV/E_e__GeV) * Phi_2);
    }

  int F_brems( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    double tau_loss[npts];
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      tau_loss[j] = taulosse_s( x[j*ndim+0] + m_e__GeV, fdata_in.n_H, fdata_in.B_G, fdata_in.n_radcomp, fdata_in.urads, fdata_in.E_peaks__GeV );
      fval[j] = sigma_BS_mb( fdata_in.E_gam, x[j*ndim+0] + m_e__GeV ) * mb_cm2 * c_cmsm1 * fdata_in.n_H/fdata_in.E_gam * gsl_spline_eval( fdata_in.gsl_so.spline, x[j*ndim+0], fdata_in.gsl_so.acc ) * tau_loss[j];
      }
    return 0;
    }

/*
  int F_brems( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    double tau_loss[npts];
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      tau_loss[j] = taulosse_s( x[j*ndim+0], fdata_in.n_H, fdata_in.B_G, fdata_in.n_radcomp, fdata_in.urads, fdata_in.E_peaks__GeV );
      fval[j] = sigma_brems_mb * mb_cm2 * c_cmsm1 * fdata_in.n_H/fdata_in.E_gam * gsl_spline_eval( fdata_in.gsl_so.spline, x[j*ndim+0], fdata_in.gsl_so.acc ) * tau_loss[j];
      }
    return 0;
    }
*/

  double xmin[1] = { fmax( Emin_e_GeV, E_gam_GeV ) };
  double xmax[1] = { Emax_e_GeV };

  struct hcub_data fdata;

  int n_radcomp = 5;
  double Epeaks[n_radcomp];
  double urads[n_radcomp];

//  double u_FIR__GeVcmm3 = urad__GeVcmm3( Labs__Lsol( SFR__Msolyrm1 ), r_e__kpc, h__pc );
//  double u_opt__GeVcmm3 = urad__GeVcmm3( Lunatt__Lsol( M_star__Msol ), r_e__kpc, h__pc );

  double u_FIR__GeVcmm3 = pow(B_G,2)/(8.*M_PI*GeV_erg);//urad__GeVcmm3( Labs__Lsol( SFR__Msolyrm1 ), r_e__kpc, h__pc );
  double u_opt__GeVcmm3 = 0.; //urad__GeVcmm3( Lunatt__Lsol( M_star__Msol ), r_e__kpc, h__pc );

  double f_urad = Labs__Lsol( SFR__Msolyrm1 )/(Labs__Lsol( SFR__Msolyrm1 ) + Lunatt__Lsol( M_star__Msol ));

  //Dust diluted black body
  Epeaks[0] = E_dBBpeak__GeV( T_dust__K );
  urads[0] = u_FIR__GeVcmm3;
  Epeaks[1] = E_BBpeak__GeV( 2.725 * (1.+z) );
  urads[1] = arad_GeVcmm3Km4*pow( 2.725 * (1.+z), 4 );
  //Draine radiation fields p.121 black (grey) bodies
  Epeaks[2] = E_BBpeak__GeV( 3000. );
  urads[2] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
  Epeaks[3] = E_BBpeak__GeV( 4000. );
  urads[3] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
  Epeaks[4] = E_BBpeak__GeV( 7500. );
  urads[4] = u_opt__GeVcmm3 * f_urad;
  
  fdata.n_radcomp = n_radcomp;
  fdata.E_peaks__GeV = Epeaks;
  fdata.urads = urads;

  fdata.B_G = B_G;
  fdata.E_gam = E_gam_GeV;
  fdata.n_H = n_H;
  fdata.gsl_so = qe_so;

  struct hcub_data *fdata_ptr = &fdata;

  if (xmin[0] < xmax[0]){
    hcubature_v( 1, F_brems, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );
    }

  return res;
  }





double eps_IC( double E_gam_GeV, double z, double B_G, double SFR__Msolyrm1, double M_star__Msol, double r_e__kpc, double h__pc, double n_H, double T_dust__K, struct gsl_spline_obj qe_so ){
  double res = 0.;
  double resout = 0.;
  double abserr;
  int i;

  //Jones 1968 - Blumenthal and Gould 1970
  double d2NdtdEgam_sm1GeVm1cm3( double E_gam__GeV, double E_peak_GeV, double E_e__GeV )
    {
    double G( double q, double Gam )
      {
      return 2.*q*log(q) + (1.+2.*q)*(1-q) + pow(Gam*q,2)*(1.-q)/(2.*(1.+Gam*q));
      }
    double gamma = E_e__GeV/m_e__GeV;
    double Gam = 4.*E_peak_GeV*gamma/m_e__GeV;
    double q = E_gam__GeV/(Gam*(E_e__GeV - E_gam_GeV));
    return 3./4. * sigma_T_mb * mb_cm2 * c_cmsm1 * pow(gamma,-2)/E_peak_GeV * G( q, Gam );
    }

  int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    double tau_loss[npts];
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);

    for (j = 0; j < npts; ++j)
      {
      tau_loss[j] = taulosse_s( x[j*ndim+0] + m_e__GeV, fdata_in.n_H, fdata_in.B_G, fdata_in.n_radcomp, fdata_in.urads, fdata_in.E_peaks__GeV );
      fval[j] = gsl_spline_eval( fdata_in.gsl_so.spline, x[j*ndim+0], fdata_in.gsl_so.acc ) * tau_loss[j] * 
                (d2NdtdEgam_sm1GeVm1cm3( fdata_in.E_gam, fdata_in.E_peak__GeV, x[j*ndim+0] + m_e__GeV ) * fdata_in.urad/fdata_in.E_peak__GeV);
      }
    return 0;
    }

  struct hcub_data fdata;


  int n_radcomp = 5;
  double Epeaks[n_radcomp];
  double urads[n_radcomp];

  double u_FIR__GeVcmm3 = pow(B_G,2)/(8.*M_PI*GeV_erg);//urad__GeVcmm3( Labs__Lsol( SFR__Msolyrm1 ), r_e__kpc, h__pc );
  double u_opt__GeVcmm3 = 0.; //urad__GeVcmm3( Lunatt__Lsol( M_star__Msol ), r_e__kpc, h__pc );

//  double u_FIR__GeVcmm3 = urad__GeVcmm3( Labs__Lsol( SFR__Msolyrm1 ), r_e__kpc, h__pc );
//  double u_opt__GeVcmm3 = urad__GeVcmm3( Lunatt__Lsol( M_star__Msol ), r_e__kpc, h__pc );

  double f_urad = Labs__Lsol( SFR__Msolyrm1 )/(Labs__Lsol( SFR__Msolyrm1 ) + Lunatt__Lsol( M_star__Msol ));

  //Dust diluted black body
  Epeaks[0] = E_dBBpeak__GeV( T_dust__K );
  urads[0] = u_FIR__GeVcmm3;
  Epeaks[1] = E_BBpeak__GeV( 2.725 * (1.+z) );
  urads[1] = arad_GeVcmm3Km4*pow( 2.725 * (1.+z), 4 );
  //Draine radiation fields p.121 black (grey) bodies
  Epeaks[2] = E_BBpeak__GeV( 3000. );
  urads[2] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
  Epeaks[3] = E_BBpeak__GeV( 4000. );
  urads[3] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
  Epeaks[4] = E_BBpeak__GeV( 7500. );
  urads[4] = u_opt__GeVcmm3 * f_urad;

  
  fdata.n_radcomp = 2; //n_radcomp;
  fdata.E_peaks__GeV = Epeaks;
  fdata.urads = urads;

  fdata.B_G = B_G;
  fdata.E_gam = E_gam_GeV;
  fdata.n_H = n_H;
  fdata.gsl_so = qe_so;

  double xmin[1], xmax[1];
  struct hcub_data *fdata_ptr = &fdata;

//  for (i = 0; i < n_radcomp; i++){
  for (i = 0; i < 2; i++){
    fdata.urad = urads[i];
    fdata.E_peak__GeV = Epeaks[i];
    xmin[0] = fmax( Emin_e_GeV, E_gam_GeV/2.*(1.+pow(1.+pow(m_e__GeV,2)/(E_gam_GeV*Epeaks[i]),1./2.)) );
    xmax[0] = Emax_e_GeV;

    if (xmin[0] < xmax[0]){
      hcubature_v( 1, F_IC, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );
      resout += res;
      }
    }

  return resout;
  }



double q_nu_norm_fcal1( double E_nu_GeV, double E_cut ){
  double res_numu2_nue, res_numu1;
  double abserr_numu2_nue, abserr_numu1;
  double x_min, x_max;


  int F_numu2_nue( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct hcub_data fdata_in = *((struct hcub_data *)fdata);
    for (j = 0; j < npts; ++j)
      {
      fval[j] = 2. *  ( f_nu_e(x[j*ndim+0]) + f_nu_mu2(x[j*ndim+0]) ) * q_pi_norm( fdata_in.E_gam/x[j*ndim+0], fdata_in.E_cut )/x[j*ndim+0];
      }
    return 0;
    }

  double F_numu1(double x, void* p){
    double lambda = 1. - pow(m_mu/m_piC,2);
    return 2./lambda * q_pi_norm( E_nu_GeV/x, E_cut )/x;
    }

  double xmin[1] = { E_nu_GeV/( (T_p_high + m_p) * K_pi) };
  double xmax[1] = { fmin( E_nu_GeV/( (T_p_low + m_p) * K_pi), 1. ) };

  struct hcub_data fdata;
  fdata.E_gam = E_nu_GeV;
  fdata.E_cut = E_cut;

  struct hcub_data *fdata_ptr = &fdata;

  if (xmin[0] < xmax[0]){
    hcubature_v( 1, F_numu2_nue, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res_numu2_nue, &abserr_numu2_nue );
    }
  else {
    res_numu2_nue = 0.;
    abserr_numu2_nue = 0.;
    }

  gsl_function F;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);

  F.function = &F_numu1;
  double lambda = 1. - pow(m_mu/m_piC,2);
  x_min = E_nu_GeV/( (T_p_high + m_p) * K_pi);
  x_max = fmin( E_nu_GeV/( (T_p_low + m_p) * K_pi), lambda );
  if (x_min < x_max){
    gsl_integration_qag( &F , x_min, x_max, 0., 1e-4, 100000, GSL_INTEG_GAUSS61, w, &res_numu1, &abserr_numu1 );
    }
  else {
    res_numu1 = 0.;
    abserr_numu1 = 0.;
    }

  gsl_integration_workspace_free(w);
  return res_numu2_nue + res_numu1;
  }


/* ############################################################################################################################## */

//Klein Nishina limit cross-section Blumenthal and Gould 1970 originally from Jones (1968)
// Photon E_gam produced on upuscaterring of photon with energy eps by electron of energy E_e
double sigma_KS_mb( double E_gam, double eps, double E_e )
  {
  double Gam_e = 4. * eps * E_e/pow(m_e,2);
  double q = E_gam/(Gam_e * (E_e - E_gam));
  return 3. * sigma_T_mb * pow(m_e,2)/(4. * pow(E_e,2) * eps) * ( 2. * q * log(q) + (1.+2.*q)*(1.-q) + 0.5 * pow(Gam_e*q,2)/(1.+Gam_e*q) * (1.-q) );
  }




double signn( double E_gam, double z, double eps, double E_e, double x, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc )
  {
  double T_CMB = 2.725 * (1. + z);
  return 2. * c_cmsm1 * sigma_KS_mb( E_gam, eps, E_e ) * mb_cm2 * dnde_gam_Planck( eps, T_CMB ) * f_nu_mu2(x) * q_pi( E_e/x, n_H, C, E_cut, spline, acc )/x;
  }







int f_int( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval ) /* x[ndim] are eps, E_e, x */
  {
  unsigned j;
  struct hcub_data fdata_in = *((struct hcub_data *)fdata);

  double T_CMB = 2.725 * (1.+fdata_in.z);

  for (j = 0; j < npts; ++j)
    {
    fval[j] = 2. * c_cmsm1 * sigma_KS_mb( fdata_in.E_gam, x[j*ndim+0], x[j*ndim+1] ) * mb_cm2 * 
              dnde_gam_Planck( x[j*ndim+0], T_CMB ) * f_nu_mu2(x[j*ndim+2] * x[j*ndim+1]) * 
              q_pi( 1./x[j*ndim+2] , fdata_in.n_H, fdata_in.C, fdata_in.E_cut, fdata_in.spline, fdata_in.acc )/(x[j*ndim+2] * x[j*ndim+1]);

    }
  return 0;
  }

int f_int_q( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval ) /* x[ndim] are eps, q, y */
  {
  unsigned j;
  struct hcub_data fdata_in = *((struct hcub_data *)fdata);

  double T_CMB = 2.725 * (1.+fdata_in.z);
  double q[npts], eps[npts], y[npts], E_e[npts];
  double sigma_KS[npts];   
  double Gam_e[npts];

  for (j = 0; j < npts; ++j)
    {
    eps[j] = x[j*ndim+0];
    q[j] = x[j*ndim+1];
    y[j] = x[j*ndim+2];
    E_e[j] = fdata_in.E_gam/2. * ( 1. + sqrt( 1. + pow(m_e,2)/(fdata_in.E_gam * eps[j] * q[j]) ) );
    Gam_e[j] = 4. * eps[j]/pow(m_e,2) * E_e[j];
    sigma_KS[j] = 3. * sigma_T_mb * mb_cm2 * pow(m_e,2)/(4. * eps[j]) * pow(Gam_e[j],2) * q[j]/pow(1.+Gam_e[j]*q[j],3) *
                  ( 2. * q[j] * log(q[j]) + (1.+2.*q[j])*(1.-q[j]) + 0.5 * pow(Gam_e[j]*q[j],2)/(1.+Gam_e[j]*q[j]) * (1.-q[j]) );

    fval[j] = 2. * c_cmsm1 * sigma_KS[j] * mb_cm2 * dnde_gam_Planck( eps[j], T_CMB ) * f_nu_mu2(y[j] * E_e[j]) * 
              q_pi( 1./y[j] , fdata_in.n_H, fdata_in.C, fdata_in.E_cut, fdata_in.spline, fdata_in.acc )/(y[j] * E_e[j]);
    }
  return 0;
  }


double q_gam_IC( double E_gam_GeV, double z, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc )  /* x[ndim] are eps, E_e, x */
  {

  double xmin[3] = { 0. , 0. , 1./( (T_p_high + m_p) * K_pi) };
  double xmax[3] = { 10 * 6.62545338e-13 * (1.+z) , 1. , fmin( 1./( (T_p_low + m_p) * K_pi), 1. ) };

  struct hcub_data fdata;
  fdata.E_gam = E_gam_GeV;
  fdata.z = z;
  fdata.C = C;
  fdata.E_cut = E_cut;
  fdata.n_H = n_H;
  fdata.spline = spline;
  fdata.acc = acc;

  struct hcub_data *fdata_ptr = &fdata;
  double val, err;

  hcubature_v( 1, f_int_q, fdata_ptr, 3, xmin, xmax, 100000, 0., 1e-3, ERROR_INDIVIDUAL, &val, &err );
  return val;
  }

//New function calc electron spectra first the integrate over them by inerpolating


int f_int_q2( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval ) /* x[ndim] are eps, q, y */
  {
  unsigned j;
  struct hcub_data fdata_in = *((struct hcub_data *)fdata);

  double T_CMB = 2.725 * (1.+fdata_in.z);
  double q[npts], eps[npts], E_e[npts];
  double sigma_KS[npts];   
  double Gam_e[npts];

  for (j = 0; j < npts; ++j)
    {
    eps[j] = x[j*ndim+0];
    q[j] = x[j*ndim+1];

    E_e[j] = fdata_in.E_gam/2. * ( 1. + sqrt( 1. + pow(m_e,2)/(fdata_in.E_gam * eps[j] * q[j]) ) );
    Gam_e[j] = 4. * eps[j]/pow(m_e,2) * E_e[j];

    sigma_KS[j] = 3. * sigma_T_mb * mb_cm2 * pow(m_e,4)/(16. * pow(eps[j],2) * pow(E_e[j],2)) * 1./( pow(q[j],2) * sqrt(1+ pow(m_e,2)/( eps[j] * q[j] * fdata_in.E_gam ) ) ) *
                  ( 2. * q[j] * log(q[j]) + (1.+2.*q[j])*(1.-q[j]) + 0.5 * pow(Gam_e[j]*q[j],2)/(1.+Gam_e[j]*q[j]) * (1.-q[j]) );

    if ( E_e[j] >= T_p_low && E_e[j] <= T_p_high ){
//    fval[j] = c_cmsm1 * sigma_KS[j] * dnde_gam_Planck( eps[j], T_CMB ) * gsl_spline_eval( &fdata_in.spline, E_e[j], &fdata_in.acc ) * fdata_in.h * pc_cm/c_cmsm1;

    fval[j] = c_cmsm1 * sigma_KS[j] * nE_gam_Planck( eps[j], fdata_in.Sigstar ) * gsl_spline_eval( &fdata_in.spline, E_e[j], &fdata_in.acc ) * fdata_in.h * pc_cm/c_cmsm1;


    if (fval[j] < 0){printf("Here\n");}
    }
    else { fval[j] = 0.;}
    }
  return 0;
  }


double q_gam_IC2( double E_gam_GeV, double z, double h_pc, double Sigstar, double n_H, double C, double E_cut, gsl_spline spline, gsl_interp_accel acc )  /* x[ndim] are eps, E_e, x */
  {
  int n_Esteps = 101.;
  double E_e_GeV[n_Esteps];
  logspace_array( n_Esteps, T_p_low, T_p_high, E_e_GeV);
  int i;

  double dndE_e[n_Esteps];

  for (i = 0; i < n_Esteps; i++){
    dndE_e[i] =  q_e( E_e_GeV[i], n_H, C, E_cut, spline, acc );
    }
  gsl_interp_accel *accEe;
  gsl_spline *splineEe;
  accEe = gsl_interp_accel_alloc();
  splineEe = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
  gsl_spline_init(splineEe, E_e_GeV, dndE_e, n_Esteps);


  double T_ISM = T_planck_gal_K( Sigstar );
  double eps_peak = 2.431e-13 * T_ISM;

//  printf("%le %le\n", eps_peak*10e9, T_ISM);

  double xmin[2] = { eps_peak/10. , 0. };
  double xmax[2] = { eps_peak*10. , 1. };

  struct hcub_data fdata;
  fdata.E_gam = E_gam_GeV;
  fdata.E_cut = E_cut;
  fdata.z = z;
  fdata.h = h_pc;
  fdata.Sigstar = Sigstar;
  fdata.spline = *splineEe;
  fdata.acc = *accEe;

  struct hcub_data *fdata_ptr = &fdata;
  double val, err;

  hcubature_v( 1, f_int_q2, fdata_ptr, 2, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &val, &err );
  return val;
  }

double N_gam_tot_obs( gsl_spline spline, gsl_interp_accel acc, double Emin, double Emax ){
  struct Ngt_int {
    gsl_spline spline;
    gsl_interp_accel acc;
    };
  struct Ngt_int fdata;
  fdata.spline = spline;
  fdata.acc = acc;

  int f_int_Ngt( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct Ngt_int fdata_in = *((struct Ngt_int *)fdata);
    for (j = 0; j < npts; ++j){fval[j] = gsl_spline_eval( &fdata_in.spline, x[j*ndim+0], &fdata_in.acc );}
    return 0;
    }

//  double xmin[1] = { Emin_GeV };
// Set this to 50 GeV+ to compare with the Ackermann work
  double xmin[1] = { Emin };
  double xmax[1] = { Emax };

  struct Ngt_int *fdata_ptr = &fdata;
  double val, err;

  hcubature_v( 1, f_int_Ngt, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &val, &err );
  return val;
  }


// Input is the spetrum 1/E
double E_gam_tot( gsl_spline spline, gsl_interp_accel acc, double Emin, double Emax ){
  struct Ngt_int {
    gsl_spline spline;
    gsl_interp_accel acc;
    };
  struct Ngt_int fdata;
  fdata.spline = spline;
  fdata.acc = acc;

  int f_int_Ngt( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct Ngt_int fdata_in = *((struct Ngt_int *)fdata);
    for (j = 0; j < npts; ++j){fval[j] = gsl_spline_eval( &fdata_in.spline, x[j*ndim+0], &fdata_in.acc ) * x[j*ndim+0];}
    return 0;
    }

  double xmin[1] = { Emin };
  double xmax[1] = { Emax };

  struct Ngt_int *fdata_ptr = &fdata;
  double val, err;

  hcubature_v( 1, f_int_Ngt, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &val, &err );
  return val;
  }

// Input is the spetrum 1/E
double N_gam_tot( gsl_spline spline, gsl_interp_accel acc, double Emin, double Emax ){
  struct Ngt_int {
    gsl_spline spline;
    gsl_interp_accel acc;
    };
  struct Ngt_int fdata;
  fdata.spline = spline;
  fdata.acc = acc;

  int f_int_Ngt( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct Ngt_int fdata_in = *((struct Ngt_int *)fdata);
    for (j = 0; j < npts; ++j){fval[j] = gsl_spline_eval( &fdata_in.spline, x[j*ndim+0], &fdata_in.acc );}
    return 0;
    }

  double xmin[1] = { Emin };
  double xmax[1] = { Emax };

  struct Ngt_int *fdata_ptr = &fdata;
  double val, err;

  hcubature_v( 1, f_int_Ngt, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &val, &err );
  return val;
  }




#endif
