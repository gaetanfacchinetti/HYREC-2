/* Structure with all energy injection parameters */
/* If adding a new energy injection process 
   make sure to add relevant parameters here */

#ifndef __ENERGY_INJECTION__
#define __ENERGY_INJECTION__

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


double dEdtdV_inj(double z, double xe, double Tgas, INJ_PARAMS *params);
double decay_rate_pmf_turbulences(double z, double tdti, double nB); /* decay rate of PMF turbulences*/

#endif
