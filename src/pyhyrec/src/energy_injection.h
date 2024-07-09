/* Structure with all energy injection parameters */
/* If adding a new energy injection process 
   make sure to add relevant parameters here */

#ifndef __ENERGY_INJECTION__
#define __ENERGY_INJECTION__


/****** CONSTANTS IN CGS + EV UNIT SYSTEM *******/

#define EI   13.598286071938324        /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */
#define mH   0.93878299831e9           /* Hydrogen atom mass in eV/c^2 */
#define _RHO_CRITICAL_ (double) (10536.723728399931) /* critical density on the universe in units of h^2 * eV / clight^2 / cm^3 */
#define _MPC_TO_CM_ (double) (3.085677581491367e+24) /* conversion factor between Mpc and cm */
#define _MU_0_ (double) (2.2401590387281866e-14) // clight^2 / eV * cm * nG^2 * s^2
#define _C_LIGHT_ (double) (2.99792458e+10) // cm / s

//------------------------------//

typedef struct {
  double h;
  double T0;
  double Omega_b, Omega_cb, Omega_k;
  double w0, wa; 
  double Neff; 
  double Nmnu; 
  double mnu1, mnu2, mnu3;             
  double YHe;
  double fsR, meR;
} INPUT_COSMOPARAMS;

typedef struct {
  double pann;
  double pann_halo; 
  double ann_z, ann_zmax, ann_zmin, ann_var;
  double ann_z_halo;
  int on_the_spot;
  double Mpbh, fpbh;
  double decay;
  double sigmaB_PMF, nB_PMF;
  double sigmaA_PMF;
  double smooth_z_PMF;
  int heat_channel_PMF;
} INPUT_INJ_PARAMS;


//------------------------------//

typedef struct {

  double odmh2;                 /* Omega_dm h^2 */
  
  double pann, pann_halo;       /* DM annihilation parameter in the smooth background and in haloes */
                                /* Units of pann and pann_halo are cm^3/s/GeV */

  double ann_z, ann_zmax, ann_zmin, ann_var; /* Parameters for the variation of pann(z) */
  double ann_z_halo;                         /* Characteristic redshift for annihilation in haloes */
  
  double decay;

  double Mpbh, fpbh;           /* Mass and fraction of DM made of primordial black holes */

  int on_the_spot;            /* if set to 1 assume energy deposition rate = injection rate */
                              /* Otherwise solves for deposition given injection with simple recipe */

  double ion, exclya, dEdtdV_heat;   /* Adding the possibility to have a heating decorrelated from ion and exclya */
  
  double sigmaB_PMF, nB_PMF;      /* adding the possibility for Primordial Magnetic Field Heating (sigmaB_PMF in nG)*/
  double sigmaA_PMF;              /* characteristic amplitude of the PMF on Alfven's scale and Jean's scale */
  double smooth_z_PMF;            /* smoothing scale for PMF energy injection around recombination*/
  
  int heat_channel_PMF;           /* 0 for both ambipolar diffusion and turbulences, 1 for turbulences, 2 for ambipolar diffusion*/

} INJ_PARAMS;


/* Structure for HYREC-2 internal parameters */

typedef struct {
  double h;                                         /* Hubble constant */
  double T0;                                        /* CMB temperature today in K*/
  double obh2, ocbh2, odeh2, okh2, orh2, onuh2;     /* density parameters */
  double w0, wa;                                    /* Dark energy equation of state parameters */
  double Neff;                                      /* total effective number of neutrinos (massive + massless) */
  double Nur;                                       /* number of massless neutrinos */
  double Nmnu;                                      /* number of massive neutrinos */
  double mnu[3];                                    /* neutrino masses */
  double fHe;                                       /* Helium fraction by number */
  double nH0;                                       /* density of hydrogen today in cm^{-3} [Changed from m^{-3} in February 2015] */
  double YHe;                                       /* Helium fraction */
  double fsR, meR;                                  /* fine-structure constant alpha/alpha(today)
                                                       and me/me(today) (Added April 2012)*/
  double dlna, nz;

  INJ_PARAMS *inj_params;                           /* Structure containing all Energy-injection parameters */

} REC_COSMOPARAMS;

double dEdtdV_inj(double z, double xe, double Tgas, INJ_PARAMS *params);
void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
		       double nH, double xH, double H, REC_COSMOPARAMS *params, double *dEdtdV_dep, 
           double *dEdtdV_ion, double *dEdtdV_exclya, double *dEdtdV_heat);
double decay_rate_pmf_turbulences(double z, double tdti, double nB); /* decay rate of PMF turbulences*/
double dEdtdV_heat_turbulences_pmf(double z, double H, double obh2, double ocbh2, double sigmaA, double sigmaB, double nB, double smooth_z);
double dEdtdV_heat_ambipolar_pmf(double z, double xe, double Tgas, double obh2, double sigmaA, double sigmaB, double nB, double smooth_z);
double compute_dEdtdV_heat_turbulences_pmf(double z, double H, INPUT_COSMOPARAMS cosmo_params, INPUT_INJ_PARAMS inj_params);
double compute_dEdtdV_heat_ambipolar_pmf(double z, double xe, double Tgas, INPUT_COSMOPARAMS cosmo_params, INPUT_INJ_PARAMS inj_params);


#endif
