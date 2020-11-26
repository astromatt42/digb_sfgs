#ifndef astro_const_h
#define astro_const_h

#include <math.h>

const double pi = M_PI;

const double c_mks = 299792458.; /* m s^-1 */
const double c_cgs = 29979245800.; /* cm s^-1 */
const double eV_J = 1.602176634e-19; /* J */
const double eV_erg = 1.602176634e-12; /* erg */


const double k_b_mks = 1.3806503e-23; /* m^2 kg s^-2 K^-1 */
const double k_b_eVK = 8.617333262145e-5; /* eV K^-1 */
const double k_b_cgs = 1.380649e-16; /* erg K^-1 */

const double h_mks = 6.62607015e-34; /* J s */
const double h_eVs = 4.135667696e-15; /* eV s */

const double E_e_eV = 510998.95; /* eV */
const double sigma_T_cgs = 6.6524587158e-25; /* cm^2 */


const double H_0 = 70.; /* km s^-1 MPc^-1 */

const double MPc_m = 3.0857e22; /* m */
const double Mpc_cm = 3.0857e24; /* cm */
const double yr_s = 365.25 * 24. * 60. * 60.; /* s */

const double foe = 1e51; /* erg */

#define D_H_MPc ((c_mks/1000.)/H_0) /* MPc */

#define D_H_cm (D_H_MPc * MPc_cm) /* cm */




#endif
