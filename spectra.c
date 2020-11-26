#include <stdio.h>
#include <math.h>
#include <gsl_interp2d.h>
#include <gsl_integration.h>
#include <gsl_sf_hyperg.h>
#include <gsl_spline.h>
#include "spectra_funcs.h"
#include "EBL_funcs.h"
#include "share/math_funcs.h"
#include "cosmo_funcs.h"

double p = 2.20; //2.23
double c_cmsm1 = 2.99792458e10;
double pc_cm = 3.0857e18;

double Msol_kg = 1.988e30;
double erg_GeV = 624.151;
double m_p_GeV = 0.9383;

double Emin_GeV = 1e-2;
double Emax_GeV = 1e8;

double E_0 = 0.9383;


int main(){
  unsigned long int i, j;

/*################################################################################################################################*/

/* Read tau EBL from file Dominguez data */

  size_t nx_E, ny_z;
//  FILE *tau_in = fopen( "input/tau_Eg_z_Dominguez.txt", "r" );
  FILE *tau_in = fopen( "input/tau_Eg_z_Franceschini.txt", "r" );
//  FILE *tau_in = fopen( "input/tau_Eg_z_Gilmore.txt", "r" );
  fscanf( tau_in , "%*[^\n]\n");
  fscanf( tau_in , "%lu \n", &ny_z );
  fscanf( tau_in , "%*[^\n]\n");
  fscanf( tau_in , "%lu \n", &nx_E );
  fscanf( tau_in , "%*[^\n]\n");

  double xa_E[nx_E];
  double ya_z[ny_z]; 
  double za_tau[nx_E * ny_z];
  double za_tau_err[nx_E * ny_z];

  /* File IO */
  /* Careful with \n characters in za_tau read-in */
  for (i=0;i<ny_z;i++){fscanf( tau_in ,"%lf", &ya_z[i] );}
  fscanf( tau_in , "\n");
  fscanf( tau_in , "%*[^\n]\n");
  for (i=0;i<nx_E;i++){fscanf( tau_in ,"%le", &xa_E[i] );}
  fscanf( tau_in , "\n");
  fscanf( tau_in , "%*[^\n]\n");
  for ( i = 0 ; i < (nx_E * ny_z) ; i++){fscanf( tau_in ,"%le", &za_tau[i] );}
  fclose(tau_in);


  /* Assign fdata for interpolation for data & error - energy in eV */
  fd_in fdata_in;
  fdata_in.nx = nx_E;
  fdata_in.ny = ny_z;
  fdata_in.xa = (double*) malloc(nx_E * sizeof(double));
  fdata_in.ya = (double*) malloc(ny_z * sizeof(double));
  fdata_in.za = (double*) malloc( (nx_E * ny_z) * sizeof(double));
  for (i=0; i < nx_E; i++) fdata_in.xa[i] = xa_E[i];
  for (i=0; i < ny_z; i++) fdata_in.ya[i] = ya_z[i];
  for (i=0; i < (ny_z * nx_E); i++) fdata_in.za[i] = za_tau[i];
  fdata_in.T = gsl_interp2d_bilinear;
  fdata_in.interp = gsl_interp2d_alloc(fdata_in.T, fdata_in.nx, fdata_in.ny);
  fdata_in.xacc = gsl_interp_accel_alloc();
  fdata_in.yacc = gsl_interp_accel_alloc();
  gsl_interp2d_init(fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, fdata_in.nx, fdata_in.ny);

  fd_in fdata_in_err;
  fdata_in_err.nx = nx_E;
  fdata_in_err.ny = ny_z;
  fdata_in_err.xa = (double*) malloc(nx_E * sizeof(double));
  fdata_in_err.ya = (double*) malloc(ny_z * sizeof(double));
  fdata_in_err.za = (double*) malloc( (nx_E * ny_z) * sizeof(double));
  for (i=0; i < nx_E; i++) fdata_in_err.xa[i] = xa_E[i];
  for (i=0; i < ny_z; i++) fdata_in_err.ya[i] = ya_z[i];
  for (i=0; i < (ny_z * nx_E); i++) fdata_in_err.za[i] = za_tau_err[i];
  fdata_in_err.T = gsl_interp2d_bilinear;
  fdata_in_err.interp = gsl_interp2d_alloc(fdata_in.T, fdata_in.nx, fdata_in.ny);
  fdata_in_err.xacc = gsl_interp_accel_alloc();
  fdata_in_err.yacc = gsl_interp_accel_alloc();
  gsl_interp2d_init(fdata_in_err.interp, fdata_in_err.xa, fdata_in_err.ya, fdata_in_err.za, fdata_in_err.nx, fdata_in_err.ny);


/*################################################################################################################################*/



  /* Read in the gals */
  /* Field z logM re_light re_mass SFR */
  /* Field is GOODS-S: 1, GOODS-N: 2, COSMOS: 3, UDS: 4 */

  FILE *gals_in = fopen( "input/cat_GOODS-S.txt" , "r" );
//  FILE *gals_in = fopen( "input/cat_test.txt" , "r" );
//  FILE *gals_in = fopen( "../gal_data/cat_test.txt" , "r" );
  fscanf( gals_in , "%*[^\n]\n");
  unsigned long int n_gal;
  fscanf( gals_in , "%lu\n" , &n_gal);
  printf("Reading %lu galaxies...\n" , n_gal);
  fscanf( gals_in , "%*[^\n]\n");
  fscanf( gals_in , "%*[^\n]\n");
  unsigned int field[n_gal];
  double z[n_gal], logM_Msol[n_gal], re_light_kpc[n_gal], re_mass_kpc[n_gal], SFR_Msolyrm1[n_gal];
  for (i = 0; i < n_gal; i++){
    fscanf( gals_in , "%u %le %le %le %le %le\n" , &field[i], &z[i], &logM_Msol[i], &re_light_kpc[i], &re_mass_kpc[i], &SFR_Msolyrm1[i] ); 
    }
  fclose(gals_in);






/*################################################################################################################################*/

/* Calculate calorimetry fraction */

  double Sig_star_Msolpcm2[n_gal], Sig_SFR_Msolyrm1pcm2[n_gal], Sig_gas_Msolpcm2[n_gal], A_re_pc2[n_gal], sig_gas_kmsm1[n_gal];
  double h_pc[n_gal], n_H_cmm3[n_gal], n_cmm3[n_gal], Phi0[n_gal];

  for (i = 0; i < n_gal; i++){


     A_re_pc2[i] = M_PI * pow( re_light_kpc[i] * 1e3, 2 );


//NGC253
////A_re_pc2[5] = SFR_Msolyrm1[5]/2. * 1e6 / 10.;
//ARP220
////A_re_pc2[2] = SFR_Msolyrm1[2]/2. * 1e6 / 489.78;
//M31
//A_re_pc2[0] = SFR_Msolyrm1[0]/2. * 1e6 / 7.41e-4;
//SMC
//A_re_pc2[3] = SFR_Msolyrm1[3]/2. * 1e6 / 1.0e-3; //Bolatto 2011 Fig 3.
//LMC
//A_re_pc2[1] = SFR_Msolyrm1[1]/2. * 1e6 / 7.41e-4;


    Sig_star_Msolpcm2[i] = pow( 10., logM_Msol[i] ) / ( 2. * A_re_pc2[i] );
    Sig_SFR_Msolyrm1pcm2[i] = SFR_Msolyrm1[i] / ( 2. * A_re_pc2[i] );
//    Sig_gas_Msolpcm2[i] = Sigma_gas_Msolpcm2_iKS( Sig_SFR_Msolyrm1pcm2[i] );
    Sig_gas_Msolpcm2[i] = Sigma_gas_Shi_Msolpcm2_iKS( Sig_SFR_Msolyrm1pcm2[i], Sig_star_Msolpcm2[i] );
    sig_gas_kmsm1[i] = sigma_gas_Yu( SFR_Msolyrm1[i] );
    }

  double chi = 1e-4, M_A = 2., beta = 0.25;
  double m_H_kg = 1.67e-27, sigma_pp_cm2 = 40e-27, mu_H = 1.4, mu_p = 1.17;

  /* Number of stars that go SN per solar mass of stars formed - C2003 IMF */
  double n_SN_Msolm1 = 1.321680e-2;
  /* fraction of energy that goes into CRs */
  double f_EtoCR = 0.1;
  /* Energy of each SN */
  double E_SN_erg = 1e51;

  int n_Esteps = 101.;

  double E_GeV[n_Esteps];
//  double E_GeV_nu[n_Esteps];
//  logspace_array( n_Esteps, Emin_GeV, Emax_GeV, E_GeV);
  logspace_array( n_Esteps, Emin_GeV, Emax_GeV, E_GeV);
  FILE *Ebins_out = fopen( "output/out_Ebins.txt", "w+" );
  for (i = 0; i < n_Esteps; i++){fprintf( Ebins_out, "%le ", E_GeV[i] );}
  fclose(Ebins_out);

  //need to allocate in heap rather than stack due to size
//  double fcal[n_gal][n_Esteps];
//  double tau_eff[n_gal][n_Esteps];
//  double D_cm2sm1[n_gal][n_Esteps];
  
  double **fcal = malloc(sizeof *fcal * n_gal);
  if (fcal){for (i = 0; i < n_gal; i++){fcal[i] = malloc(sizeof *fcal[i] * n_Esteps);}}
  double **tau_eff = malloc(sizeof *tau_eff * n_gal);
  if (tau_eff){for (i = 0; i < n_gal; i++){tau_eff[i] = malloc(sizeof *tau_eff[i] * n_Esteps);}}
  double **D_cm2sm1 = malloc(sizeof *D_cm2sm1 * n_gal);
  if (D_cm2sm1){for (i = 0; i < n_gal; i++){D_cm2sm1[i] = malloc(sizeof *D_cm2sm1[i] * n_Esteps);}}

  double G_h = 4.302e-3, eta_pp = 0.5;
  double u_LA, v_Ai, t_loss_s, E_SNCR_GeVsm1, C[n_gal] ;
  double v_st, L_A, Gam_0, D0;

  double CnormE = C_norm_E();

  for (i = 0; i < n_gal; i++){
    h_pc[i] = pow( sig_gas_kmsm1[i], 2 )/( M_PI * G_h * ( Sig_star_Msolpcm2[i] + Sig_gas_Msolpcm2[i] ) ); //or 12 instead of pi

    n_H_cmm3[i] = Sig_gas_Msolpcm2[i]/( mu_H * m_H_kg * 2. * h_pc[i] ) * Msol_kg/pow( pc_cm, 3 );

    n_cmm3[i] = n_H_cmm3[i] * mu_H/mu_p;

    u_LA = sig_gas_kmsm1[i]/sqrt(2.);
    v_Ai = 1000. * ( u_LA/10. )/( sqrt(chi/1e-4) * M_A );
    L_A = h_pc[i]/pow(M_A,3);

    D0 = v_Ai * L_A * 1e5 * pc_cm;

    t_loss_s = 1./(1./(1./( n_cmm3[i] * sigma_pp_cm2 * eta_pp * c_cmsm1 )) + 1./(pow(10 * h_pc[i] * pc_cm,2)/D0) );
    E_SNCR_GeVsm1 = SFR_Msolyrm1[i] * n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg_GeV/yr_s;

    C[i] = E_SNCR_GeVsm1 * t_loss_s/( CnormE * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) );




    for (j = 0; j < n_Esteps; j++){

      v_st = fmin( v_Ai * (1. + 2.3e-3 * pow( sqrt(pow(E_GeV[j],2) + 2. * m_p_GeV * E_GeV[j]) , p-1.) * pow(n_H_cmm3[i]/1e3, 1.5) * (chi/1e-4) * M_A/( u_LA/10. * C[i]/2e-7 )), c_cmsm1/1e5);
      D_cm2sm1[i][j] = v_st * L_A * 1e5 * pc_cm;


      tau_eff[i][j] = 9.9 * Sig_gas_Msolpcm2[i]/1e3 * h_pc[i]/1e2 * 1e27/D_cm2sm1[i][j];

      Gam_0 = 41.2 * h_pc[i]/1e2 * v_st/1e3 * 1e27/D_cm2sm1[i][j];
      fcal[i][j] = 1. - 1./( gsl_sf_hyperg_0F1( beta/(beta+1.) , tau_eff[i][j]/pow(beta+1.,2) ) + 
                   tau_eff[i][j]/Gam_0 * gsl_sf_hyperg_0F1( (beta+2.)/(beta+1.) , tau_eff[i][j]/pow(beta+1., 2)) );
      }
    //C here is used for the calcultion with dsig/dE, so need t_loss = t_col
    t_loss_s = 1./( n_cmm3[i] * sigma_pp_cm2 * eta_pp * c_cmsm1 );
    C[i] = E_SNCR_GeVsm1 * t_loss_s/( CnormE * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) );

    }

  printf("Done calculating galaxies.\n");


/*################################################################################################################################*/
/*################################################################################################################################*/

/* Calculate spectra */

  double **specs = malloc(sizeof *specs * n_gal);
  if (specs){for (i = 0; i < n_gal; i++){specs[i] = malloc(sizeof *specs[i] * n_Esteps);}}
  double **specs_fcal1 = malloc(sizeof *specs_fcal1 * n_gal);
  if (specs_fcal1){for (i = 0; i < n_gal; i++){specs_fcal1[i] = malloc(sizeof *specs_fcal1[i] * n_Esteps);}}
  double **specs_nu = malloc(sizeof *specs_nu * n_gal);
  if (specs_nu){for (i = 0; i < n_gal; i++){specs_nu[i] = malloc(sizeof *specs_nu[i] * n_Esteps);}}
  double **tau_gg = malloc(sizeof *tau_gg * n_gal);
  if (tau_gg){for (i = 0; i < n_gal; i++){tau_gg[i] = malloc(sizeof *tau_gg[i] * n_Esteps);}}
  double **tau_EBL = malloc(sizeof *tau_EBL * n_gal);
  if (tau_EBL){for (i = 0; i < n_gal; i++){tau_EBL[i] = malloc(sizeof *tau_EBL[i] * n_Esteps);}}

  double **specs_IC = malloc(sizeof *specs_IC * n_gal);
  if (specs_IC){for (i = 0; i < n_gal; i++){specs_IC[i] = malloc(sizeof *specs_IC[i] * n_Esteps);}}
  double **specs_casc_obs = malloc(sizeof *specs_casc_obs * n_gal);
  if (specs_casc_obs){for (i = 0; i < n_gal; i++){specs_casc_obs[i] = malloc(sizeof *specs_casc_obs[i] * n_Esteps);}}
  double **specs_casc = malloc(sizeof *specs_casc * n_gal);
  if (specs_casc){for (i = 0; i < n_gal; i++){specs_casc[i] = malloc(sizeof *specs_casc[i] * n_Esteps);}}

  double **specs_obs = malloc(sizeof *specs_obs * n_gal);
  if (specs_obs){for (i = 0; i < n_gal; i++){specs_obs[i] = malloc(sizeof *specs_obs[i] * n_Esteps);}}

  double **specs_N_emit = malloc(sizeof *specs_N_emit * n_gal);
  if (specs_N_emit){for (i = 0; i < n_gal; i++){specs_N_emit[i] = malloc(sizeof *specs_N_emit[i] * n_Esteps);}}

/* Spline for calorimetry fraction */
  gsl_interp_accel *acc;
  gsl_spline *spline;

  double mod, vol;
  struct Phi_out Phi_out;


  double C_gam;
  double spec_emit[n_Esteps];

  printf("%s \n", "Calculating spectra...");

  #pragma omp parallel for schedule(dynamic) private(j, spline, acc, Phi_out, mod, vol, spec_emit, C_gam )
  for (i = 0; i < n_gal; i++){
    if (field[i] == 1){
      acc = gsl_interp_accel_alloc();
      spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
      printf("%lu \n", i);
      gsl_spline_init(spline, E_GeV, fcal[i], n_Esteps);
      vol = 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3);
      mod = vol * pow((1.+z[i]),2)/( 4. * M_PI * pow(d_l_MPc( z[i] ) * Mpc_cm, 2) );
      for (j = 0; j < n_Esteps; j++){
        Phi_out = Phi( (1.+z[i]) * E_GeV[j], n_H_cmm3[i], C[i], *spline, *acc );
        tau_gg[i][j] = fmax(0., tau_pair( (1.+z[i]) * E_GeV[j], z[i], h_pc[i], Sig_star_Msolpcm2[i] ));
        tau_EBL[i][j] = fmax(0., gsl_interp2d_eval_extrap( fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, 
                                           fmin(E_GeV[j] * 1e9, 1e15), z[i], fdata_in.xacc, fdata_in.yacc ));

        specs[i][j] = Phi_out.Phi * mod;
        specs_obs[i][j] = Phi_out.Phi * mod * exp(-tau_gg[i][j]) * exp(-tau_EBL[i][j]);
        specs_fcal1[i][j] = Phi_out.Phi_fcal1 * mod;
        specs_nu[i][j] = q_nu( (1.+z[i]) * E_GeV[j], n_H_cmm3[i], C[i], *spline, *acc ) * mod;

        //call the function again this time without redshift
        Phi_out = Phi( E_GeV[j], n_H_cmm3[i], C[i], *spline, *acc );
        //number of photons emitted
        specs_N_emit[i][j] = Phi_out.Phi * vol * exp(-tau_gg[i][j]);
        //for cascade
        spec_emit[j] = Phi_out.Phi * exp(-tau_gg[i][j]);
        }
      C_gam = norm_casc_C( spec_emit, E_GeV, n_Esteps, z[i], fdata_in );
      for (j = 0; j < n_Esteps; j++){
        specs_casc_obs[i][j] = dndE_gam_casc( (1.+z[i]) * E_GeV[j], z[i], C_gam, fdata_in ) * mod;
        specs_casc[i][j] = dndE_gam_casc( E_GeV[j], z[i], C_gam, fdata_in ) * vol;
        }
  
      }
    }

/*################################################################################################################################*/
// Calculate total emission observed

  double N_gam_emit_obs[n_gal];

  printf("%s \n", "Calculating photon luminosity...");

  #pragma omp parallel for schedule(dynamic) private( j, spline, acc, spec_emit )
  for (i = 0; i < n_gal; i++){
    if (field[i] == 1){
      for (j = 0; j < n_Esteps; j++){spec_emit[j] = specs_obs[i][j] + specs_casc_obs[i][j];}
      acc = gsl_interp_accel_alloc();
      spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
      gsl_spline_init(spline, E_GeV, spec_emit, n_Esteps);
      N_gam_emit_obs[i] = N_gam_tot_obs( *spline, *acc, 1., 100. );
    }
  }

/*################################################################################################################################*/

  double L_gam_emit[n_gal];
  double N_gam_emit[n_gal];

  printf("%s \n", "Calculating total luminosity...");

  #pragma omp parallel for schedule(dynamic) private( j, spline, acc, spec_emit )
  for (i = 0; i < n_gal; i++){
    if (field[i] == 1){
      for (j = 0; j < n_Esteps; j++){spec_emit[j] = specs_N_emit[i][j] + specs_casc[i][j];}
      acc = gsl_interp_accel_alloc();
      spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
      gsl_spline_init(spline, E_GeV, spec_emit, n_Esteps);
      L_gam_emit[i] = E_gam_tot( *spline, *acc, 0.1, 100. );
      N_gam_emit[i] = N_gam_tot( *spline, *acc, 0.1, 100. );
    }
  }


/*################################################################################################################################*/

  double taugg, tauEBL;

  printf("%s \n", "Writing ouput files...");

/* Output file */
  FILE *specs_out = fopen( "output/out_gal_specs.txt", "w+" );
  FILE *specs_out_fcal1 = fopen( "output/out_gal_specs_fcal1.txt", "w+" );
  FILE *specs_out_nt = fopen( "output/out_gal_specs_notaugg.txt", "w+" );
  FILE *specs_out_ne = fopen( "output/out_gal_specs_noEBL.txt", "w+" );
  FILE *specs_out_nt_ne = fopen( "output/out_gal_specs_notaugg_noEBL.txt", "w+" );
  FILE *specs_out_nu = fopen( "output/out_gal_specs_nu.txt", "w+" );
  FILE *specs_out_IC = fopen( "output/out_gal_specs_IC.txt", "w+" );
  FILE *specs_out_cascade = fopen( "output/out_gal_specs_cascade.txt", "w+" );
  FILE *fcal_out = fopen( "output/out_gal_fcal.txt", "w+" );

  FILE *Ngtobs_out = fopen( "output/out_Ngt_obs.txt", "w+" );
  FILE *Lgt_out = fopen( "output/out_Lgt.txt", "w+" );
  FILE *Ngt_out = fopen( "output/out_Ngt.txt", "w+" );

  for (i = 0; i < n_gal; i++){
    if (field[i] == 1){
      fprintf( Ngtobs_out, "%e \n", N_gam_emit_obs[i] );
      fprintf( Lgt_out, "%e \n", L_gam_emit[i] );
      fprintf( Ngt_out, "%e \n", N_gam_emit[i] );
      for (j = 0; j < n_Esteps; j++){
        taugg = tau_gg[i][j];
        tauEBL = tau_EBL[i][j];
        fprintf( specs_out, "%e ", specs[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_ne, "%e ", specs[i][j] * exp(-taugg) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_nt, "%e ", specs[i][j] * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_nu, "%e ", specs_nu[i][j] * pow( E_GeV[j], 2 ));
        fprintf( specs_out_IC, "%e ", specs_IC[i][j]  * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_cascade, "%e ", specs_casc_obs[i][j]  * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_nt_ne, "%e ", specs[i][j] * pow( E_GeV[j], 2 ));
        fprintf( specs_out_fcal1, "%e ", specs_fcal1[i][j] * pow( E_GeV[j], 2 ));
        fprintf( fcal_out, "%e ", fcal[i][j] );
        }
      fprintf( specs_out, "\n" );
      fprintf( specs_out_nt, "\n" );
      fprintf( specs_out_ne, "\n" );
      fprintf( specs_out_nu, "\n" );
      fprintf( specs_out_IC, "\n" );
      fprintf( specs_out_cascade, "\n" );
      fprintf( specs_out_nt_ne, "\n" );
      fprintf( specs_out_fcal1, "\n" );
      fprintf( fcal_out, "\n" );
      }
    }  

  fclose(specs_out);
  fclose(specs_out_nt);
  fclose(specs_out_ne);
  fclose(specs_out_nu);
  fclose(specs_out_IC);
  fclose(specs_out_cascade);
  fclose(specs_out_nt_ne);
  fclose(specs_out_fcal1);
  fclose(fcal_out);
  fclose(Ngtobs_out);
  fclose(Lgt_out);
  fclose(Ngt_out);

  printf("%s \n", "Done");

/*################################################################################################################################*/

  return 0;
  }
