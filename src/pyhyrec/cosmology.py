########################
# Some simple functions implemented to compute simple cosmological quantities
# from the result of HYREC runs - without invoking involved codes like CLASS
########################

import numpy as np
from scipy import integrate

from .wrapperhyrec import compute_hubble_rate, call_run_hyrec
from .params import HyRecInjectionParams, HyRecCosmoParams


# define usefull units and conversion rates
_MPC_TO_M_      = 3.08567758128e+22
_MSUN_TO_KG_    = 2.0e+30
_KG_TO_EV_      = 5.60958860380445e+35
_KM_TO_MPC_     = 1.0/_MPC_TO_M_ * 1e+3
_M_TO_MPC_     =  1.0/_MPC_TO_M_
_MASS_HYDROGEN_ = 0.93878299831e+9 # in eV
_MASS_PROTON_   = 0.938272e+9 # in eV
_MASS_HELIUM_   = 3.728400e+9 # in eV
_SIGMA_THOMSON_ = 6.6524616e-29 # in m^2
_C_LIGHT_       = 299792458 # in m / s


def hubble_rate(z, cosmo = HyRecCosmoParams()):
    """
        hubble rate in s^{-1} directly computed from the C-code

    Parameters:
    -----------
    - z: np.ndarray, float
        redshift
    - cosmo: HyRecCosmoParams, optional
        cosmology

    Returns:
    --------
    hubble rate in s^{-1}
    """

    _IS_A_LIST = True
    if not isinstance(z, (list, np.ndarray)):
        _IS_A_LIST = False
        z = np.array([z])

    res = np.zeros(len(z))

    for iz, _z in enumerate(z):
        res[iz] = compute_hubble_rate(_z, cosmo(), HyRecInjectionParams()())

    if _IS_A_LIST == True:
        return res
    
    return res[0]


def hubble_factor(z, cosmo=HyRecCosmoParams()):
    """
        hubble factor H(z)/H_0 (dimensionless)

    Parameters:
    -----------
    - z: np.ndarray, float
        redshift
    - cosmo: HyRecCosmoParams, optional
        cosmology

    Returns:
    --------
    hubble factor H(z)/H0
    """
    return hubble_rate(z, cosmo) / 3.2407792896393e-18 / cosmo.h


def rho_baryons(cosmo = HyRecCosmoParams()):
    return 2.7754e+11 * cosmo.Omega_b * cosmo.h**2 * _MSUN_TO_KG_ * _KG_TO_EV_ / _MPC_TO_M_**3 # in eV / m^3

def n_baryons(cosmo = HyRecCosmoParams()):
    return rho_baryons(cosmo) / _MASS_PROTON_ / (1 + cosmo.YHe / 4 * (_MASS_HELIUM_/_MASS_HYDROGEN_ -1)) # in 1/m^3

def rho_gamma(cosmo = HyRecCosmoParams()):
    omega_gamma = 4.48162687719e-7 * cosmo.T0**4 / cosmo.h**2
    return 2.7754e+11 * omega_gamma * cosmo.h**2 * _MSUN_TO_KG_ * _KG_TO_EV_ / _MPC_TO_M_**3 # in eV / m^3

def rho_radiation(cosmo):
    neff_array  = [3.046, 2.0328, 1.0196, 0.00641]
    m_neutrinos = [cosmo.mnu1, cosmo.mnu2, cosmo.mnu3]
    n_ur = neff_array[np.count_nonzero(m_neutrinos)]
    omega_r = 4.48162687719e-7 * cosmo.T0**4 * (1. + 0.227107317660239 * n_ur) / cosmo.h**2

    return 2.7754e+11 * omega_r * cosmo.h**2 * _MSUN_TO_KG_ * _KG_TO_EV_ / _MPC_TO_M_**3 # in eV / m^3

def optical_depth(z, xe, cosmo = HyRecCosmoParams()):
    """
        optical depth assuming no reionization (dimensionless)

    Parameters:
    -----------
    - z: np.ndarray
        redshift
    - xe: np.ndarray
        free electron fraction
    - cosmo: HyRecCosmoParams, optional
        cosmology

    Returns:
    --------
    optical depth for each value of z
    """
    
    e_z = hubble_factor(z)

    # fast trapezoid integration scheme
    integrand = xe * (1+z)**2 / e_z
    trapz = (integrand[1:] + integrand[:-1])/2.0
    dz = np.diff(z)

    res = np.zeros(len(z))
    for i in range(len(z)-1):
        res[i+1] = res[i] + trapz[i] * dz[i]

    pref = _C_LIGHT_ * _SIGMA_THOMSON_ * n_baryons(cosmo) / (cosmo.h * 3.2407792896393e-18)
    return pref * res


def visibility_function(z, xe, cosmo):
    pref = _C_LIGHT_ * _SIGMA_THOMSON_ * n_baryons(cosmo) / (cosmo.h * 3.2407792896393e-18)
    return pref * (1 + z)**2 / hubble_factor(z) * xe * np.exp(-optical_depth(z, xe, cosmo))

def z_rec(z, xe, cosmo):
    vis = visibility_function(z, xe, cosmo)
    return z[np.argmax(vis)]

def compute_z_rec(cosmo = HyRecCosmoParams()):
    z, xe, _ = call_run_hyrec(cosmo(),  HyRecInjectionParams()(), zmax = 8000, zmin = 500, nz = 40000)
    return z_rec(z, xe, cosmo)
    
def acoustic_damping_scale(z, xe, cosmo):
    
    # get the value of the recombination redshift
    z_reco = z_rec(z, xe, cosmo)
    iz_min = np.argmin(np.abs(z - z_reco)) 
    
    # restrict the redshift and free electron fraction 
    # to the interesting interval
    z_int  = z[iz_min:]
    xe_int = xe[iz_min:]

    # get the baryon and photon densities
    rho_b  = rho_baryons(cosmo) * (1+z_int)**3
    rho_g  = rho_gamma(cosmo) * (1+z_int)**4
    r_arr = 3./4. * rho_b / rho_g

    # evaluate the Thomson scattering length
    l_scatt = 1.0/(_SIGMA_THOMSON_ * xe_int * n_baryons(cosmo) * (1+z_int)**3 ) * _M_TO_MPC_ # in Mpc

    # hubble rate in 1/Mpc
    h_z = hubble_rate(z_int) /_M_TO_MPC_ / _C_LIGHT_ # in 1/Mpc

    # damping length square
    ld2 = integrate.trapezoid((1+z_int)/(1+ r_arr)/h_z * l_scatt * (16.0/15.0 + r_arr**2/(1+r_arr)), z_int)/6.0
    
    return np.sqrt(1/ld2)

def compute_acoustic_damping_scale(cosmo = HyRecCosmoParams()):
    z, xe, _ = call_run_hyrec(cosmo(),  HyRecInjectionParams()(), zmax = 10000, zmin = 500, nz = 40000)
    return acoustic_damping_scale(z, xe, cosmo)