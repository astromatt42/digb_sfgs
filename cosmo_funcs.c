#ifndef cosmo_funcs_h
#define cosmo_funcs_h
#include <math.h>
#include <gsl_integration.h>
#include "../share/cosmo_params.h"
#include "../share/astro_const.h"

double E_z( double z ){
  return sqrt(om_M * pow((1.+z),3) + om_k * pow((1.+z),2) + om_lam);
  }

/* Comoving distance [MPc] */
double d_c_MPc( double z_low, double z_high ){
  if (z_low < 0.) return 0.;
  if (z_high <= z_low) return 0.; 
  double result;
  double abserr;
  double f(double z, void* p){return 1./E_z( z );}
  gsl_function F;
  F.function = &f;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);
  gsl_integration_qag( &F , z_low, z_high, 0., 1e-8, 10000, GSL_INTEG_GAUSS61, w, &result, &abserr );
  gsl_integration_workspace_free(w);
  return result * D_H_MPc;
  }

/* Transverse comoving distance [MPc] */
double d_m_MPc( double z ){
  if (z < 0.) return 0.;
  if (om_k == 0.){return d_c_MPc( 0., z );}
  else if (om_k > 0.){return D_H_MPc/sqrt(om_k) * sinh( sqrt(om_k) * d_c_MPc( 0., z )/D_H_MPc );}
  else if (om_k < 0.){return D_H_MPc/sqrt(fabs(om_k)) * sin( sqrt(fabs(om_k)) * d_c_MPc( 0., z )/D_H_MPc );}
  return 0.;
  }

/* Luminosity distance [MPc] */
double d_l_MPc( double z ){
  if (z <= 0.) return 0.;
  return d_m_MPc( z ) * (1. + z);
  }

/* Angular diameter distance [MPc] */
double d_a_MPc( double z ){
  if (z < 0.) return 0.;
  return d_m_MPc( z ) / (1. + z);
  }

/* Comoving volume element [MPc^3 d\Omega^-1 dz^-1] */
double dV_c( double z ){
  if (z < 0.) return 0.;
  return D_H_MPc * pow( (1. + z) * d_a_MPc( z ), 2 )/E_z( z ) ;
  }

/* Comoving volume between z_low and z_high [MPc^3 d\Omega^-1] */
double Vc_dOm_MPc3srm1( double z_low, double z_high ){
  if ( (z_low < 0.) || (z_high <= z_low) ) return 0.;
  double result;
  double abserr;
  double f(double z, void* p){return dV_c( z );}
  gsl_function F;
  F.function = &f;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);
  gsl_integration_qag( &F , z_low, z_high, 0., 1e-8, 10000, GSL_INTEG_GAUSS61, w, &result, &abserr );
  gsl_integration_workspace_free(w);
  return result;
  }





#endif
 
