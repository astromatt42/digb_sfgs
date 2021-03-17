#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl_interp2d.h>
#include <gsl_integration.h>
#include <gsl_sf_hyperg.h>
#include <gsl_spline.h>
#include "spectra_funcs.h"
#include "EBL_funcs.h"
#include "share/math_funcs.h"
#include "cosmo_funcs.h"

double q = 2.2; //2.23
double c_cmsm1 = 2.99792458e10;
double pc_cm = 3.0857e18;

double Msol_kg = 1.988e30;
double erg_GeV = 624.151;
double m_p_GeV = 0.9383;

double Emin_GeV = 1e-3;
double Emax_GeV = 1e8;

double Emin_e_GeV = 1.e-6; //m_e_GeV
double Emax_e_GeV = 1.e8;


double E_0 = 0.9383;


double Delta_x[11] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10. };
double Phi_1H[11] = { 45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65 };
double Phi_2H[11] = { 44.46, 44.38, 44.24, 43.65, 42.49, 40.19, 34.93, 29.78, 24.34, 17.28, 12.41 };


int main( int argc, char *argv[] ){
  unsigned long int i, j, k;


//  printf("Normalisation full/1GeV+: %e \n", C_norm_E_e( E_cutoff )/C_norm_E_e( E_cutoff_e ));


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

//  FILE *gals_in = fopen( "input/cat_GOODS-S.txt" , "r" );
//  FILE *gals_in = fopen( "input/cat_test.txt" , "r" );

  char infp[] = "input/";
  char infile[strlen(argv[1]) + strlen(infp) + 1];
  snprintf(infile, strlen(argv[1]) + strlen(infp) + 1, "%s%s", infp, argv[1]);
  FILE *gals_in = fopen( infile , "r" );

//  char outfp[] =

  fscanf( gals_in , "%*[^\n]\n");
  unsigned long int n_gal;
  fscanf( gals_in , "%lu\n" , &n_gal);
  printf("Reading %lu galaxies from file %s...\n" , n_gal, infile);
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
  double h_pc[n_gal], n_H_cmm3[n_gal], n_cmm3[n_gal], Phi0[n_gal], Mstar__Msol[n_gal];

  #pragma omp parallel for schedule(dynamic)
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

    Mstar__Msol[i] = pow( 10., logM_Msol[i] );
    Sig_star_Msolpcm2[i] = pow( 10., logM_Msol[i] ) / ( 2. * A_re_pc2[i] );
    Sig_SFR_Msolyrm1pcm2[i] = SFR_Msolyrm1[i] / ( 2. * A_re_pc2[i] );
//    Sig_gas_Msolpcm2[i] = Sigma_gas_Msolpcm2_iKS( Sig_SFR_Msolyrm1pcm2[i] );
    Sig_gas_Msolpcm2[i] = Sigma_gas_Shi_Msolpcm2_iKS( Sig_SFR_Msolyrm1pcm2[i], Sig_star_Msolpcm2[i] );
    sig_gas_kmsm1[i] = sigma_gas_Yu( SFR_Msolyrm1[i] );
    }

  double chi = 1e-4, M_A = 2.0, beta = 0.25;
  double m_H_kg = 1.67e-27, sigma_pp_cm2 = 40e-27, mu_H = 1.4, mu_p = 1.17;

  /* Number of stars that go SN per solar mass of stars formed - C2003 IMF */
  double n_SN_Msolm1 = 1.321680e-2;
  /* fraction of energy that goes into CRs */
  double f_EtoCR = 0.1;
  /* Energy of each SN */
  double E_SN_erg = 1e51;

  int n_Esteps = 101;

  double E_GeV[n_Esteps];

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

  double **CR_spec = malloc(sizeof *CR_spec * n_gal);
  if (CR_spec){for (i = 0; i < n_gal; i++){CR_spec[i] = malloc(sizeof *CR_spec[i] * n_Esteps);}}

  double **CR_spec_e1 = malloc(sizeof *CR_spec_e1 * n_gal);
  if (CR_spec_e1){for (i = 0; i < n_gal; i++){CR_spec_e1[i] = malloc(sizeof *CR_spec_e1[i] * n_Esteps);}}



  double G_h = 4.302e-3, eta_pp = 0.5;
  double u_LA, v_Ai, t_loss_s, E_SNCR_GeVsm1, C[n_gal], Ce_Esm1[n_gal];
  double v_st, L_A, Gam_0, D0, B_G[n_gal];

  double DcascKgv, DcascKrn, DfastTTD, t_casc, r_L, Dflrw, r_G_cm, l_outer, LdampA, Ldampf, eps_f, l_MA1, c_s, mu_cff;

  double CnormE[n_gal], CnormE_e[n_gal], E_cut[n_gal];

  #pragma omp parallel for schedule(dynamic) private(j, u_LA, v_Ai, L_A, l_outer, D0, t_loss_s, E_SNCR_GeVsm1, LdampA, Ldampf, v_st, r_G_cm, l_MA1, DcascKgv, eps_f, c_s, mu_cff, DfastTTD, Gam_0 )
  for (i = 0; i < n_gal; i++){
    h_pc[i] = pow( sig_gas_kmsm1[i], 2 )/( M_PI * G_h * ( Sig_star_Msolpcm2[i] + Sig_gas_Msolpcm2[i] ) );

    n_H_cmm3[i] = Sig_gas_Msolpcm2[i]/( mu_H * m_H_kg * 2. * h_pc[i] ) * Msol_kg/pow( pc_cm, 3 );

    n_cmm3[i] = n_H_cmm3[i] * mu_H/mu_p;

    u_LA = sig_gas_kmsm1[i]/sqrt(2.);
    v_Ai = 1000. * ( u_LA/10. )/( sqrt(chi/1e-4) * M_A );
    L_A = h_pc[i]/pow(M_A,3);

    B_G[i] = sqrt(4.* M_PI * chi * n_cmm3[i] * mu_p * m_H_kg * 1e3) * v_Ai * 1e5;
    l_outer = h_pc[i] * pc_cm; 

    //work out the cutoff energy for each galaxy by assuming Hillas criterion with injection at h and Bfield amplified by 10?.
    E_cut[i] = h_pc[i] * pc_cm * B_G[i]/1e4/330.;

//    CnormE[i] = C_norm_E(E_cut[i]);
    CnormE[i] = C_norm_E(E_cutoff);

    D0 = v_Ai * L_A * 1e5 * pc_cm;

    t_loss_s = 1./(1./(1./( n_cmm3[i] * sigma_pp_cm2 * eta_pp * c_cmsm1 )) + 1./(pow( h_pc[i] * pc_cm,2)/D0) );
    E_SNCR_GeVsm1 = SFR_Msolyrm1[i] * n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg_GeV/yr_s;

    C[i] = E_SNCR_GeVsm1 * t_loss_s/( CnormE[i] * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) );

//Output these for paper
//printf( "Phi = %e \n", n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg_GeV/yr_s/C_norm_E(E_cutoff)/pow(m_p_GeV,q) );
//printf( "integrated over E = %e \n", n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg_GeV/yr_s );

    Ce_Esm1[i] = 0.2*E_SNCR_GeVsm1/( C_norm_E_e(E_cutoff_e) * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) );


    LdampA = 0.0011 * pow( h_pc[i]/100., -1./2.) * pow( u_LA/10./ ( (chi/1e-4) * (n_H_cmm3[i]/1e3) ), 3./2.);
    Ldampf = 1.7 * pow( u_LA/10., 2./3.) * pow( sqrt(h_pc[i]/100.)/( (chi/1e-4) * (n_H_cmm3[i]/1e3) ), 1./3.) / pow(M_A, 2);


    for (j = 0; j < n_Esteps; j++){
      v_st = fmin( v_Ai * (1. + 2.3e-3 * pow( sqrt(pow(E_GeV[j],2) + 2. * m_p_GeV * E_GeV[j]) , q-1.) * pow(n_H_cmm3[i]/1e3, 1.5) * (chi/1e-4) * M_A/( u_LA/10. * C[i]/2e-7 )), c_cmsm1/1e5);
      D_cm2sm1[i][j] = v_st * L_A * 1e5 * pc_cm;

      r_G_cm = 330. * E_GeV[j] * pow( B_G[i]/1e4, -1 );

      l_MA1 = l_outer/pow( M_A, 3 ); 

      if ( r_G_cm > l_MA1 ){
        //Kolmogorov ISO-K41 weak isotropic cascade
        DcascKgv = c_cmsm1/3 * pow( M_A, -2 ) * l_outer * pow( r_G_cm/l_outer, 1./3.);
        }
      else {
        //Xu Lazarian 2017 l_parallel + Hopkins 2020
        DcascKgv = c_cmsm1/3 * l_MA1 ; //* pow( r_G_cm/l_MA1, 1./2.);
        }

      //Fraction of turbulent power in fast modes
      eps_f = 0.5;
      //Kraichnan - make sure to adjust damping scale to fast modes
      DcascKrn = c_cmsm1 * pow( M_A, -2 ) * l_outer * pow( r_G_cm/l_outer, 1./2.); //pow( r_G_cm * l_outer, 1./2.) * c_cmsm1/(3. * eps_f);
//printf( "%e %e %e \n", D_cm2sm1[i][j], DcascKgv, DcascKrn );

      c_s = 15.;

      mu_cff = fmin(0.99, pow( 14./M_PI * pow(1.15, 2) * pow(r_G_cm/pc_cm/h_pc[i], 1./2.), 2./11. ));

      DfastTTD = 28./(5.*M_PI) * pow(1.15, -2) * pow(r_G_cm/pc_cm/h_pc[i], -1./2.) * r_G_cm * c_cmsm1 * (4.-sqrt(mu_cff)*(5.-pow(mu_cff,2)));

//      DfastTTD = pow(c_cmsm1,2)/4. * 1./(sqrt(2.)/4 * pow(M_PI,3./2.) * 1./(16.*M_PI) * pow(h_pc[i], -1./2.) * pow(u_LA/sqrt(3)/c_s, 2) * c_cmsm1/pc_cm * ( pow( fmax( r_G_cm/pc_cm, Ldampf ), -1./2) - pow(h_pc[i], -1./2) ));


//printf( "%e %e %e \n", E_GeV[j], D_cm2sm1[i][j], DfastTTD );

        if (r_G_cm/pc_cm/LdampA * 2.*M_PI * 100. > 1.){

//          D_cm2sm1[i][j] = fmin(DfastTTD,D_cm2sm1[i][j]);
          }

/*
      if ( Ldampf > LdampA ){
        if (r_G_cm/pc_cm/Ldampf * 2.*M_PI > 1.){
          D_cm2sm1[i][j] = fmin(DcascKrn,D_cm2sm1[i][j]);
          }
        else if (r_G_cm/pc_cm/LdampA * 2.*M_PI > 1.){
          D_cm2sm1[i][j] = fmin(DcascKgv,D_cm2sm1[i][j]);
          }
        }
      else {
        if (r_G_cm/pc_cm/Ldampf * 2.*M_PI > 1.){
          D_cm2sm1[i][j] = DcascKrn;
          }
        }

*/

//      if (r_G_cm/pc_cm/LdampA * 2.*M_PI > 1.){
//      if (r_G_cm/pc_cm/Ldampf * 2.*M_PI > 1.){
//        D_cm2sm1[i][j] = Dcasc;
//        }

//      printf( "%e %e %e \n", E_GeV[j], D_cm2sm1[i][j], Dcasc );
//      printf( "%e %e \n", E_GeV[j], r_G_cm/pc_cm/LdampA );

      tau_eff[i][j] = 9.9 * Sig_gas_Msolpcm2[i]/1e3 * h_pc[i]/1e2 * 1e27/D_cm2sm1[i][j];

      Gam_0 = 41.2 * h_pc[i]/1e2 * v_st/1e3 * 1e27/D_cm2sm1[i][j];
      fcal[i][j] = 1. - 1./( gsl_sf_hyperg_0F1( beta/(beta+1.) , tau_eff[i][j]/pow(beta+1.,2) ) + 
                   tau_eff[i][j]/Gam_0 * gsl_sf_hyperg_0F1( (beta+2.)/(beta+1.) , tau_eff[i][j]/pow(beta+1., 2)) );
      }
    //C here is used for the calcultion with dsig/dE, so need t_loss = t_col
    t_loss_s = 1./( n_cmm3[i] * sigma_pp_cm2 * eta_pp * c_cmsm1 );
    C[i] = E_SNCR_GeVsm1 * t_loss_s/( CnormE[i] * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) );

    //Calculate the internal cosmic ray spectrum with loss times.
    for (j = 0; j < n_Esteps; j++){
//      t_loss_s = 1./(1./(1./( n_cmm3[i] * sigma_pp_cm2 * eta_pp * c_cmsm1 )) + 1./(pow(10 * h_pc[i] * pc_cm,2)/D_cm2sm1[i][j]) );
//      CR_spec[i][j] = J( E_GeV[j], E_SNCR_GeVsm1 * t_loss_s/( CnormE[i] * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) ), E_cut[i] );
      CR_spec[i][j] = J( E_GeV[j], E_SNCR_GeVsm1/( CnormE[i] * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3) ), E_cutoff ) * fcal[i][j];
      CR_spec_e1[i][j] = J_e( E_GeV[j], E_SNCR_GeVsm1/10./( C_norm_E_e(E_cutoff_e) ), E_cutoff_e );
      }

//    if (i==7){
//      for (j = 0; j < n_Esteps; j++){
//        printf( "%e %e \n", D_cm2sm1[i][j], E_GeV[j] );
//        }
//      }



    }

  /* print out the energy densities */
  FILE *Udens_out = fopen( "output/Udens_out.txt", "w+" );
  fprintf( Udens_out, "%s\n", "U_mag U_FIR U_star Sig_gas" );
  for (i = 0; i < n_gal; i++){
    fprintf( Udens_out, "%e %e %e %e \n", pow(B_G[i],2)/(8*M_PI*GeV_erg), urad__GeVcmm3( Labs__Lsol( SFR_Msolyrm1[i] ), re_light_kpc[i], h_pc[i] ), urad__GeVcmm3( Lunatt__Lsol( pow(10.,logM_Msol[i]) ), re_light_kpc[i], h_pc[i] ), Sig_gas_Msolpcm2[i] * Msol_kg*1.e3/pow(pc_cm,2) );
    }

//  for (i = 0; i < n_gal; i++){
//    printf( "%e %e \n", urad__GeVcmm3( Labs__Lsol( SFR_Msolyrm1[i] ), re_light_kpc[i], h_pc[i] ) + urad__GeVcmm3( Lunatt__Lsol( pow(10.,logM_Msol[i]) ), re_light_kpc[i], h_pc[i] ), uFIR_GeVsm1cmm3( LFIR_Lsol( SFR_Msolyrm1[i] ), re_light_kpc[i], h_pc[i] ) );
//    }

  printf("Done calculating galaxies.\n");

  //print CR spectra
  FILE *CR_specs_out = fopen( "output/CR_specs.txt", "w+" );
  for (j = 0; j < n_Esteps; j++){
    fprintf( CR_specs_out, "%e ", E_GeV[j] );
    }
  fprintf( CR_specs_out, "\n" );
  for (i = 0; i < n_gal; i++){
    for (j = 0; j < n_Esteps; j++){
      fprintf( CR_specs_out, "%e ", CR_spec[i][j] * pow(E_GeV[j],2) * c_cmsm1/(4.*M_PI) );
      }
    fprintf( CR_specs_out, "\n" );
    }
  fclose(CR_specs_out);


  FILE *galdata_out = fopen( "output/gal_data.txt", "w+" );
  fprintf( galdata_out, "h_pc nH_cmm3 B_G sigmag_kmsm1 Are_pc2\n" );
  for (i = 0; i < n_gal; i++){
    fprintf( galdata_out, "%e %e %e %e %e\n", h_pc[i], n_H_cmm3[i], B_G[i], sig_gas_kmsm1[i], A_re_pc2[i] );
    }
  fclose(galdata_out);


  unsigned long int n_radcomp = 5;
  double Epeaks[n_radcomp];
  double urads[n_radcomp];
  
  double u_FIR__GeVcmm3, u_opt__GeVcmm3, f_urad;

  double T_dust__K;

//    double tau_loss = taulosse_s( E_e__GeV, n_H, B__G, n_radcomp, urads, Epeaks );


  double tICm1;

  //print loss times
  FILE *tau_loss_out = fopen( "output/tau_loss_leptons.txt", "w+" );
  for (i = 0; i < n_gal; i++){

    u_FIR__GeVcmm3 = urad__GeVcmm3( Labs__Lsol( SFR_Msolyrm1[i] ), re_light_kpc[i], h_pc[i] );
    u_opt__GeVcmm3 = urad__GeVcmm3( Lunatt__Lsol( pow(10.,logM_Msol[i]) ), re_light_kpc[i], h_pc[i] );

    f_urad = Labs__Lsol( SFR_Msolyrm1[i] )/(Labs__Lsol( SFR_Msolyrm1[i] ) + Lunatt__Lsol( pow(10.,logM_Msol[i]) ));

    T_dust__K = Tdust_K( z[i], SFR_Msolyrm1[i], pow(10., logM_Msol[i]) );

    //Dust diluted black body
    Epeaks[0] = E_dBBpeak__GeV( T_dust__K );
    urads[0] = u_FIR__GeVcmm3;
    //Draine radiation fields p.121 black (grey) bodies
    Epeaks[1] = E_BBpeak__GeV( 3000. );
    urads[1] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
    Epeaks[2] = E_BBpeak__GeV( 4000. );
    urads[2] = u_opt__GeVcmm3 * (1. - f_urad) * 0.5;
    Epeaks[3] = E_BBpeak__GeV( 7500. );
    urads[3] = u_opt__GeVcmm3 * f_urad;
    Epeaks[4] = E_BBpeak__GeV( 2.725 * (1.+z[i]) );
    urads[4] = arad_GeVcmm3Km4*pow( 2.725 * (1.+z[i]), 4 );

    for (j = 0; j < n_Esteps; j++){
      fprintf( tau_loss_out, "%e ", tauBS( E_GeV[j], n_H_cmm3[i] ) );
      }
    fprintf( tau_loss_out, "\n" );
    for (j = 0; j < n_Esteps; j++){
      fprintf( tau_loss_out, "%e ", tausync( E_GeV[j], B_G[i] ) );
      }
    fprintf( tau_loss_out, "\n" );
    for (j = 0; j < n_Esteps; j++){
      tICm1 = 0.;

      for (k = 0; k < n_radcomp; k++){
        tICm1 += 1./tauIC( E_GeV[j], urads[k], Epeaks[k] );
        }
      fprintf( tau_loss_out, "%e ", 1./tICm1 );
      }
    fprintf( tau_loss_out, "\n" );
    for (j = 0; j < n_Esteps; j++){
      fprintf( tau_loss_out, "%e ", tauion( E_GeV[j], n_H_cmm3[i] ) );
      }
    fprintf( tau_loss_out, "\n" );
  }


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


  double **specs_casc_obs = malloc(sizeof *specs_casc_obs * n_gal);
  if (specs_casc_obs){for (i = 0; i < n_gal; i++){specs_casc_obs[i] = malloc(sizeof *specs_casc_obs[i] * n_Esteps);}}
  double **specs_casc = malloc(sizeof *specs_casc * n_gal);
  if (specs_casc){for (i = 0; i < n_gal; i++){specs_casc[i] = malloc(sizeof *specs_casc[i] * n_Esteps);}}

  double **specs_brems_1st = malloc(sizeof *specs_brems_1st * n_gal);
  if (specs_brems_1st){for (i = 0; i < n_gal; i++){specs_brems_1st[i] = malloc(sizeof *specs_brems_1st[i] * n_Esteps);}}

  double **specs_brems_2nd = malloc(sizeof *specs_brems_2nd * n_gal);
  if (specs_brems_2nd){for (i = 0; i < n_gal; i++){specs_brems_2nd[i] = malloc(sizeof *specs_brems_2nd[i] * n_Esteps);}}

  double **specs_sync_1st = malloc(sizeof *specs_sync_1st * n_gal);
  if (specs_sync_1st){for (i = 0; i < n_gal; i++){specs_sync_1st[i] = malloc(sizeof *specs_sync_1st[i] * n_Esteps);}}

  double **specs_sync_2nd = malloc(sizeof *specs_sync_2nd * n_gal);
  if (specs_sync_2nd){for (i = 0; i < n_gal; i++){specs_sync_2nd[i] = malloc(sizeof *specs_sync_2nd[i] * n_Esteps);}}

  double **specs_IC_1st = malloc(sizeof *specs_IC_1st * n_gal);
  if (specs_IC_1st){for (i = 0; i < n_gal; i++){specs_IC_1st[i] = malloc(sizeof *specs_IC_1st[i] * n_Esteps);}}

  double **specs_IC_2nd = malloc(sizeof *specs_IC_2nd * n_gal);
  if (specs_IC_2nd){for (i = 0; i < n_gal; i++){specs_IC_2nd[i] = malloc(sizeof *specs_IC_2nd[i] * n_Esteps);}}

  double **specs_obs = malloc(sizeof *specs_obs * n_gal);
  if (specs_obs){for (i = 0; i < n_gal; i++){specs_obs[i] = malloc(sizeof *specs_obs[i] * n_Esteps);}}

  double **specs_L_emit = malloc(sizeof *specs_L_emit * n_gal);
  if (specs_L_emit){for (i = 0; i < n_gal; i++){specs_L_emit[i] = malloc(sizeof *specs_L_emit[i] * n_Esteps);}}

  double **q_e_2nd = malloc(sizeof *q_e_2nd * n_gal);
  if (q_e_2nd){for (i = 0; i < n_gal; i++){q_e_2nd[i] = malloc(sizeof *q_e_2nd[i] * n_Esteps);}}

  double **q_e_1st = malloc(sizeof *q_e_1st * n_gal);
  if (q_e_1st){for (i = 0; i < n_gal; i++){q_e_1st[i] = malloc(sizeof *q_e_1st[i] * n_Esteps);}}


/* Spline for calorimetry fraction */
  gsl_interp_accel *acc;
  gsl_spline *spline;

  double mod, vol;
  struct Phi_out Phi_out;


  double C_gam;
  double spec_emit[n_Esteps];

  int n_Esteps_e = 1001;
  double E_GeV_e[n_Esteps_e];
  double q_e_1st_nEe[n_Esteps_e];


  logspace_array( n_Esteps_e, Emin_e_GeV, Emax_e_GeV, E_GeV_e);

  struct gsl_spline_obj qe1st_so, qe2nd_so;

  double Lradio[n_gal];

  double sb1, sb2, sIC1, sIC2, ss1, ss2;

  printf("%s \n", "Calculating spectra...");


  #pragma omp parallel for schedule(dynamic) private(j, spline, acc, Phi_out, mod, vol, spec_emit, C_gam, qe2nd_so, qe1st_so, q_e_1st_nEe, T_dust__K, sb1, sb2, sIC1, sIC2, ss1, ss2 )
  for (i = 0; i < n_gal; i++){
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
    printf("%lu \n", i);
    gsl_spline_init(spline, E_GeV, fcal[i], n_Esteps);
    vol = 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3);
    mod = vol * pow((1.+z[i]),2)/( 4. * M_PI * pow(d_l_MPc( z[i] ) * Mpc_cm, 2) );

    T_dust__K = Tdust_K( z[i], SFR_Msolyrm1[i], pow(10., logM_Msol[i]) );



    for (j = 0; j < n_Esteps_e; j++){
      /* Compute the spectra for the secondary electrons and then interpolate on them */
      q_e_1st_nEe[j] = J_e( E_GeV_e[j], Ce_Esm1[i], E_cutoff_e );
      }

    qe1st_so = gsl_so( E_GeV_e, q_e_1st_nEe, n_Esteps_e );

    for (j = 0; j < n_Esteps; j++){
      /* Compute the spectra for the secondary electrons and then interpolate on them */
      q_e_2nd[i][j] = q_e( E_GeV[j] + m_e__GeV, n_H_cmm3[i], C[i], E_cutoff, *spline, *acc );
      q_e_1st[i][j] = gsl_spline_eval( qe1st_so.spline, E_GeV[j], qe1st_so.acc );
      }


    qe2nd_so = gsl_so( E_GeV, q_e_2nd[i], n_Esteps );

    Lradio[i] = (eps_sync( 1.49e9 * h_GeVs, z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so ) +
                 eps_sync( 1.49e9 * h_GeVs, z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so )) * 
                 1.49e9 * h__Js * h_GeVs * 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3);


    for (j = 0; j < n_Esteps; j++){

      Phi_out = Phi( (1.+z[i]) * E_GeV[j], n_H_cmm3[i], C[i], E_cutoff, *spline, *acc );

      tau_gg[i][j] = fmax(0., tau_gg_gal( (1.+z[i]) * E_GeV[j], z[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i] ));
      tau_EBL[i][j] = fmax(0., gsl_interp2d_eval_extrap( fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, fmin(E_GeV[j] * 1e9, 1e15), z[i], fdata_in.xacc, fdata_in.yacc ));

      specs[i][j] = Phi_out.Phi * mod;



      specs_brems_1st[i][j] = eps_brems( (1.+z[i]) * E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so ) * mod;
      specs_brems_2nd[i][j] = eps_brems( (1.+z[i]) * E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so ) * mod;

      specs_IC_1st[i][j] = eps_IC( (1.+z[i]) * E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so ) * mod;
      specs_IC_2nd[i][j] = eps_IC( (1.+z[i]) * E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so ) * mod;

      specs_sync_1st[i][j] = eps_sync( (1.+z[i]) * E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so ) * mod;
      specs_sync_2nd[i][j] = eps_sync( (1.+z[i]) * E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so ) * mod;

      specs_fcal1[i][j] = Phi_out.Phi_fcal1 * mod;
      specs_nu[i][j] = q_nu( (1.+z[i]) * E_GeV[j], n_H_cmm3[i], C[i], E_cutoff, *spline, *acc ) * mod;


      specs_obs[i][j] = (Phi_out.Phi + specs_brems_1st[i][j] + specs_brems_2nd[i][j] + specs_IC_1st[i][j] + specs_IC_2nd[i][j] + 
                         specs_sync_1st[i][j] + specs_sync_2nd[i][j])* mod * exp(-tau_gg[i][j]) * exp(-tau_EBL[i][j]);

      //call the function again this time without redshift
      Phi_out = Phi( E_GeV[j], n_H_cmm3[i], C[i], E_cutoff, *spline, *acc );
      sb1 = eps_brems( E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so );
      sb2 = eps_brems( E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so );

      sIC1 = eps_IC( E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so );
      sIC2 = eps_IC( E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so );

      ss1 = eps_sync( E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe1st_so );
      ss2 = eps_sync( E_GeV[j], z[i], B_G[i], SFR_Msolyrm1[i], Mstar__Msol[i], re_light_kpc[i], h_pc[i], n_H_cmm3[i], T_dust__K, qe2nd_so );
      //total energy emitted
      specs_L_emit[i][j] = ( Phi_out.Phi + sb1 + sb2 +sIC1 + sIC2 + ss1 + ss2 ) * vol * exp(-tau_gg[i][j]);
      //for cascade
      spec_emit[j] = ( Phi_out.Phi + sb1 + sb2 +sIC1 + sIC2 + ss1 + ss2 ) * exp(-tau_gg[i][j]);
      }
    C_gam = norm_casc_C( spec_emit, E_GeV, n_Esteps, z[i], fdata_in );
    for (j = 0; j < n_Esteps; j++){
      specs_casc_obs[i][j] = dndE_gam_casc( (1.+z[i]) * E_GeV[j], z[i], C_gam, fdata_in ) * mod;
      specs_casc[i][j] = dndE_gam_casc( E_GeV[j], z[i], C_gam, fdata_in ) * vol;
      }
    }

  FILE *CR_specs_out_2 = fopen( "output/CR_specs_2.txt", "w+" );
  for (i = 0; i < n_gal; i++){

    vol = 2. * A_re_pc2[i] * 2. * h_pc[i] * pow(pc_cm, 3);

    for (j = 0; j < n_Esteps; j++){
      fprintf( CR_specs_out_2, "%e ", CR_spec[i][j] * vol );
      }
    fprintf( CR_specs_out_2, "\n" );

    for (j = 0; j < n_Esteps; j++){
      fprintf( CR_specs_out_2, "%e ", q_e_1st[i][j] * vol );
      }
    fprintf( CR_specs_out_2, "\n" );

    for (j = 0; j < n_Esteps; j++){
      fprintf( CR_specs_out_2, "%e ", q_e_2nd[i][j] * vol );
      }
    fprintf( CR_specs_out_2, "\n" );
    }
  fclose(CR_specs_out_2);


  FILE *FIR_radio_out = fopen( "output/FIR_radio.txt", "w+" );
  for (i = 0; i < n_gal; i++){
    fprintf( FIR_radio_out, "%e %e\n", Lradio[i], Labs__Lsol( SFR_Msolyrm1[i] ) );
    }
  fclose(FIR_radio_out);

/*################################################################################################################################*/
// Calculate total emission observed

  double N_gam_emit_obs[n_gal];

  printf("%s \n", "Calculating photon luminosity...");

  #pragma omp parallel for schedule(dynamic) private( j, spline, acc, spec_emit )
  for (i = 0; i < n_gal; i++){
    for (j = 0; j < n_Esteps; j++){spec_emit[j] = specs_obs[i][j] + specs_casc_obs[i][j];}
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
    gsl_spline_init(spline, E_GeV, spec_emit, n_Esteps);
    N_gam_emit_obs[i] = N_gam_tot_obs( *spline, *acc, 1., 100. );
    }

/*################################################################################################################################*/

  double L_gam_emit[n_gal];
  double N_gam_emit[n_gal];

  printf("%s \n", "Calculating total luminosity...");

  #pragma omp parallel for schedule(dynamic) private( j, spline, acc, spec_emit )
  for (i = 0; i < n_gal; i++){
    for (j = 0; j < n_Esteps; j++){spec_emit[j] = specs_L_emit[i][j] + specs_casc[i][j];}
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_linear, n_Esteps);
    gsl_spline_init(spline, E_GeV, spec_emit, n_Esteps);
    L_gam_emit[i] = E_gam_tot( *spline, *acc, 0.1, 100. );
    N_gam_emit[i] = N_gam_tot( *spline, *acc, 1., 100. );
    }


/*################################################################################################################################*/

  double taugg, tauEBL;

  printf("%s \n", "Writing ouput files...");

/* Output file */
  FILE *specs_out = fopen( "output/out_gal_specs.txt", "w+" );
  FILE *specs_out_fcal1 = fopen( "output/out_gal_specs_fcal1.txt", "w+" );
  FILE *specs_out_fcal1_tau = fopen( "output/out_gal_specs_fcal1_tau.txt", "w+" );
  FILE *specs_out_nt = fopen( "output/out_gal_specs_notaugg.txt", "w+" );
  FILE *specs_out_ne = fopen( "output/out_gal_specs_noEBL.txt", "w+" );
  FILE *specs_out_nt_ne = fopen( "output/out_gal_specs_notaugg_noEBL.txt", "w+" );
  FILE *specs_out_nu = fopen( "output/out_gal_specs_nu.txt", "w+" );
  FILE *specs_out_brems1 = fopen( "output/out_gal_specs_brems1.txt", "w+" );
  FILE *specs_out_brems2 = fopen( "output/out_gal_specs_brems2.txt", "w+" );
  FILE *specs_out_IC1 = fopen( "output/out_gal_specs_IC1.txt", "w+" );
  FILE *specs_out_IC2 = fopen( "output/out_gal_specs_IC2.txt", "w+" );
  FILE *specs_out_sync1 = fopen( "output/out_gal_specs_sync1.txt", "w+" );
  FILE *specs_out_sync2 = fopen( "output/out_gal_specs_sync2.txt", "w+" );
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
        fprintf( specs_out_brems1, "%e ", specs_brems_1st[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_brems2, "%e ", specs_brems_2nd[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_IC1, "%e ", specs_IC_1st[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_IC2, "%e ", specs_IC_2nd[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_sync1, "%e ", specs_sync_1st[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_sync2, "%e ", specs_sync_2nd[i][j] * exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_cascade, "%e ", specs_casc_obs[i][j]  * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_nt_ne, "%e ", (specs[i][j] + specs_brems_1st[i][j] + specs_brems_2nd[i][j] + specs_IC_1st[i][j] + 
                                          specs_IC_2nd[i][j] + specs_sync_1st[i][j] + specs_sync_2nd[i][j]) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_fcal1, "%e ", (specs_fcal1[i][j] + specs_brems_1st[i][j] + specs_brems_2nd[i][j] + specs_IC_1st[i][j] + 
                                          specs_IC_2nd[i][j] + specs_sync_1st[i][j] + specs_sync_2nd[i][j]) * pow( E_GeV[j], 2 ));
        fprintf( specs_out_fcal1_tau, "%e ", (specs_fcal1[i][j] + specs_brems_1st[i][j] + specs_brems_2nd[i][j] + 
                                              specs_IC_1st[i][j] + specs_IC_2nd[i][j] + specs_sync_1st[i][j] + specs_sync_2nd[i][j]) * 
                                              exp(-taugg) * exp(-tauEBL) * pow( E_GeV[j], 2 ));
        fprintf( fcal_out, "%e ", fcal[i][j] );
        }
      fprintf( specs_out, "\n" );
      fprintf( specs_out_nt, "\n" );
      fprintf( specs_out_ne, "\n" );
      fprintf( specs_out_nu, "\n" );
      fprintf( specs_out_brems1, "\n" );
      fprintf( specs_out_brems2, "\n" );
      fprintf( specs_out_IC1, "\n" );
      fprintf( specs_out_IC2, "\n" );
      fprintf( specs_out_sync1, "\n" );
      fprintf( specs_out_sync2, "\n" );
      fprintf( specs_out_cascade, "\n" );
      fprintf( specs_out_nt_ne, "\n" );
      fprintf( specs_out_fcal1, "\n" );
      fprintf( specs_out_fcal1_tau, "\n" );
      fprintf( fcal_out, "\n" );

      }
    }  

  fclose(specs_out);
  fclose(specs_out_nt);
  fclose(specs_out_ne);
  fclose(specs_out_nu);
  fclose(specs_out_brems1);
  fclose(specs_out_brems2);
  fclose(specs_out_IC1);
  fclose(specs_out_IC2);
  fclose(specs_out_sync1);
  fclose(specs_out_sync2);
  fclose(specs_out_cascade);
  fclose(specs_out_nt_ne);
  fclose(specs_out_fcal1);
  fclose(specs_out_fcal1_tau);
  fclose(fcal_out);
  fclose(Ngtobs_out);
  fclose(Lgt_out);
  fclose(Ngt_out);

  printf("%s \n", "Done");

/*################################################################################################################################*/

  return 0;
  }
