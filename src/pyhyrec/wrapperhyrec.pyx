
import warnings
import numpy as np
import os

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    _filename = 'pyhyrec/' + filename.split("/")[-1]
    return ' %s:%s: %s (%s)\n' % (_filename, lineno, message, category.__name__)

warnings.formatwarning = warning_on_one_line

cdef extern from "src/history.h":

    double test_cython(double x, double y)

    ctypedef struct INPUT_INJ_PARAMS:

        double pann,               #/* DM annihilation parameter in the smooth background and in haloes */ 
        double pann_halo           #/* Units of pann and pann_halo are cm^3/s/GeV */
        double ann_z, ann_zmax, ann_zmin, ann_var #/* Parameters for the variation of pann(z) */
        double ann_z_halo                         #/* Characteristic redshift for annihilation in haloes */   
        int on_the_spot             #/* if set to 1 assume energy deposition rate = injection rate */ /* Otherwise solves for deposition given injection with simple recipe */
        double Mpbh, fpbh           #/* Mass and fraction of DM made of primordial black holes */
        double decay

    ctypedef struct INPUT_COSMOPARAMS:

        double h                                         # Hubble constant 
        double T0                                        #/* CMB temperature today in K*/
        double Omega_b, Omega_cb, Omega_k                # Abundances of baryons, matter and curvature
        double w0, wa                                    #/* Dark energy equation of state parameters */
        double Neff                                      #/* total effective number of neutrinos (massive + massless) */
        double Nmnu                                      #/* number of massive neutrinos */
        double mnu1, mnu2, mnu3                          #/* neutrino masses */
        double YHe                                       #/* Helium fraction */
        double fsR, meR                                  #/* fine-structure constant alpha/alpha(today) and me/me(today) (Added April 2012)*/

    ctypedef struct HYREC_DATA:
        int error
        char * error_message
        pass

    HYREC_DATA * run_hyrec(INPUT_COSMOPARAMS cosmo, INPUT_INJ_PARAMS inj_params, double zmax, double zmin)
    void hyrec_free(HYREC_DATA * data)
    double hyrec_xe(double z, HYREC_DATA * data)
    double hyrec_Tm(double z, HYREC_DATA * data)
    

def call_test_cython(double x, double y) :
    return test_cython(x, y)

def call_run_hyrec(INPUT_COSMOPARAMS cosmo_params, INPUT_INJ_PARAMS inj_params, double zmax = 8000.0, double zmin = 0.0, int nz = 8000):

    location = os.path.join( os.path.dirname(os.path.realpath(__file__)), 'data')

    f = open(os.path.join(location, 'Alpha_inf.dat'), 'r')
    print(f.read())

    data = run_hyrec(cosmo_params, inj_params, zmax, zmin)
    
    ## If something went wrong we print the error message
    if data.error == 1 :
        print(str(data.error_message.decode('utf-8')))

    z_array = np.linspace(np.max([zmin, 1.0]), zmax, nz)
    xe_array = np.zeros(nz)
    Tm_array = np.zeros(nz)
    
    for iz, z in enumerate(z_array):
        xe_array[iz] = hyrec_xe(z, data)
        Tm_array[iz] = hyrec_Tm(z, data)
    
    # Free the memory at the end
    hyrec_free(data)
    
    return z_array, xe_array, Tm_array

    

def init_INPUT_INJ_PARAMS(double pann, double pann_halo, 
                        double ann_z, double ann_zmax, double ann_zmin, double ann_var, 
                        double ann_z_halo, double decay, int on_the_spot,
                        double Mpbh, double fpbh):
    
    cdef INPUT_INJ_PARAMS inj_params
    
    inj_params.pann  = pann
    inj_params.pann_halo = pann_halo
    inj_params.ann_z = ann_z
    inj_params.ann_zmax = ann_zmax
    inj_params.ann_zmin = ann_zmin
    inj_params.ann_var = ann_var
    inj_params.ann_z_halo = ann_z_halo
    inj_params.decay = decay
    inj_params.Mpbh = Mpbh
    inj_params.fpbh = fpbh
    inj_params.on_the_spot = on_the_spot
    
    return inj_params


def init_INPUT_COSMOPARAMS(double h, double T0, 
        double Omega_b, double Omega_cb, double Omega_k, 
        double w0, double wa, double Neff, double Nmnu,
        double mnu1, double mnu2, double mnu3,
        double YHe, double fsR, double meR):

    cdef INPUT_COSMOPARAMS cosmo

    cosmo.h = h
    cosmo.T0 = T0
    cosmo.Omega_b = Omega_b
    cosmo.Omega_cb = Omega_cb
    cosmo.Omega_k = Omega_k
    cosmo.w0 = w0
    cosmo.wa = wa
    cosmo.Neff = Neff
    cosmo.Nmnu = Nmnu
    cosmo.mnu1 = mnu1
    cosmo.mnu2 = mnu2
    cosmo.mnu3 = mnu3
    cosmo.YHe = YHe
    cosmo.fsR = fsR
    cosmo.meR = meR

    return cosmo;
