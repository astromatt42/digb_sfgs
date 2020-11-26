#ifndef EBL_funcs_h
#define EBL_funcs_h

//Cascades
// Berezniski and Kalashev 2016
extern double m_e;
extern double Emin_GeV;
extern double Emax_GeV;

double eps_EBL( double z ){return 0.6695 * (1.+z)/1e9;} //GeV
double eps_CMB( double z ){return 6.627e-4 * (1.+z)/1e9;} //GeV 2.795 K CMB

double eps_gam( double z ){return pow(m_e,2)/eps_EBL( z );}
double eps_X( double z ){return pow(m_e,2) * eps_CMB( z )/( 3. * pow(eps_EBL( z ),2) );}

double dndE_gam_casc( double E_gam, double z, double C, fd_in fdata_in )
  {
  double tauEBL = fmax( 0., gsl_interp2d_eval_extrap( fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, 
                                            fmin(E_gam/(1.+z) * 1e9, 1e15), z, fdata_in.xacc, fdata_in.yacc ) );
//  double tausrcgg = tau_pair( (1.+z[i]) * E_GeV[j], h_pc[i], Sig_star_Msolpcm2[i] )
  if ( E_gam <= eps_X(z) )
//    {return C * pow(E_gam/eps_X(z), -3./2.) * exp(-tauEBL);}
    {return C * pow(E_gam/eps_X(z), -3./2.);}
  else if ( E_gam > eps_X(z) )
//    {return C * pow(E_gam/eps_X(z), -2.) * exp(-tauEBL);}
    {return C * pow(E_gam/eps_X(z), -2.);}
  else
    {printf("Error in dndE_gam_casc"); return 0.;}
  }


double norm_casc_C1( double z, fd_in fdata_in )
  {
  double result;
  double abserr;
  double res_out = 0.;
  double f(double E_gam, void* p){return dndE_gam_casc( E_gam, z, 1., fdata_in ) * E_gam;}
  gsl_function F;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);
  F.function = &f;
  gsl_integration_qag( &F , 0., 1e6, 0., 1e-4, 100000, GSL_INTEG_GAUSS61, w, &result, &abserr );
//  gsl_integration_qag( &F , 0., eps_gam( z ), 0., 1e-8, 100000, GSL_INTEG_GAUSS61, w, &result, &abserr );
  gsl_integration_workspace_free(w);
  return 1./result;
  }

/*
double norm_casc_C( double spec_gam[], double E_GeV[], int n_Esteps, double z, fd_in fdata_in)
  {
  double result;
  double abserr;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
  double f(double E_gam, void* p)
    {
    double tauEBL = fmax( 0., gsl_interp2d_eval_extrap( fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, 
                                              fmin(E_gam/(1.+z) * 1e9, 1e15), z, fdata_in.xacc, fdata_in.yacc ) );
    double spec = gsl_spline_eval( spline, E_gam, acc );
    return (1.-exp(-tauEBL)) * spec * E_gam;
    }
  gsl_function F;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(100000);
  gsl_spline_init(spline, E_GeV, spec_gam, n_Esteps);
  F.function = &f;
//  printf("Test1 \n");

  gsl_integration_qag( &F ,  Emin_GeV, Emax_GeV, 0., 1e-4, 10000, GSL_INTEG_GAUSS61, w, &result, &abserr );
//./run   printf("Test2 \n"); 
  gsl_integration_workspace_free(w);
  return result * norm_casc_C1( z, fdata_in );
  }
*/


double norm_casc_C( double spec_gam[], double E_GeV[], int n_Esteps, double z, fd_in fdata_in)
  {

  struct int_data {
    double z;
    fd_in fdata_in;
    gsl_spline spline;
    gsl_interp_accel acc;
    };

  int f_Ccasc_int( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
    unsigned j;
    struct int_data fdata_in = *((struct int_data *)fdata);
    double tauEBL[npts];   
    double spec[npts];

    for (j = 0; j < npts; ++j)
      {
      tauEBL[j] = fmax( 0., gsl_interp2d_eval_extrap( fdata_in.fdata_in.interp, fdata_in.fdata_in.xa, fdata_in.fdata_in.ya, 
                  fdata_in.fdata_in.za, fmin(x[j*ndim+0]/(1.+fdata_in.z) * 1e9, 1e15), fdata_in.z, fdata_in.fdata_in.xacc,    
                  fdata_in.fdata_in.yacc ) );

      spec[j] = gsl_spline_eval( &fdata_in.spline, x[j*ndim+0], &fdata_in.acc );
      fval[j] = (1.-exp(-tauEBL[j])) * spec[j] * x[j*ndim+0];
      }
    return 0;
    }

  gsl_interp_accel *acc;
  gsl_spline *spline;
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
  gsl_spline_init(spline, E_GeV, spec_gam, n_Esteps);

  struct int_data fdata;
  fdata.z = z;
  fdata.fdata_in = fdata_in;
  fdata.spline = *spline;
  fdata.acc = *acc;

  struct int_data *fdata_ptr = &fdata;
  double val, err;

  double xmin[1] = { Emin_GeV };
  double xmax[1] = { Emax_GeV };

  hcubature_v( 1, f_Ccasc_int, fdata_ptr, 1, xmin, xmax, 100000, 0., 1e-4, ERROR_INDIVIDUAL, &val, &err );
  return val * norm_casc_C1( z, fdata_in );
  }


double tauEBL_Peretti( double E_gam_GeV, double z )
  {
  return 95./1100. * ( pow( pow(E_gam_GeV/1e3,-2.7)/2.1 + pow(E_gam_GeV/1e3,-0.31)/0.34 , -1. ) +
         pow( pow(E_gam_GeV/12e3,-3.1)/0.47 + pow(E_gam_GeV/40e3,-0.8)/20. , -1. ) + 7. * pow(E_gam_GeV/100e3,7.8) ) * z/0.003;
  }














#endif
