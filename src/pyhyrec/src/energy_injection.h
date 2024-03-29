/* Structure with all energy injection parameters */
/* If adding a new energy injection process 
   make sure to add relevant parameters here */

#ifndef __ENERGY_INJECTION__
#define __ENERGY_INJECTION__


/****** CONSTANTS IN CGS + EV UNIT SYSTEM *******/

#define EI   13.598286071938324        /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */


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

#endif
