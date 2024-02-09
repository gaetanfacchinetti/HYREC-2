import pyhyrec as pyhy
import numpy as np

# Test that we can define INJ_PARAMS and COSMOPARAMS objects
def test_load():
    pyhy.init_INPUT_INJ_PARAMS(0., 0., 0., 0., 0., 0., 0., 0, 0., 1., 0.)
    pyhy.init_INPUT_COSMOPARAMS(6.735837e-01, 2.7255, 0.0494142797907188, 0.31242079216478097, 0., -1, 0, 3.046, 1.0, 0.06, 0., 0., 0.245, 1., 1.)

# Test that we can run the HYREC C code
def test_call_run():
    inj_params = pyhy.init_INPUT_INJ_PARAMS(0., 0., 0., 0., 0., 0., 0., 0, 0., 1., 0.)
    cosmo_params = pyhy.init_INPUT_COSMOPARAMS(6.735837e-01, 2.7255, 0.0494142797907188, 0.31242079216478097, 0., -1, 0, 3.046, 1.0, 0.06, 0., 0., 0.245, 1., 1.)
    pyhy.call_run_hyrec(cosmo_params, inj_params)

# Check here that the result given with the defaut options match with the default output
def test_default_output():
    
    inj_params = pyhy.init_INPUT_INJ_PARAMS(0., 0., 0., 0., 0., 0., 0., 0, 0., 1., 0.)
    cosmo_params = pyhy.init_INPUT_COSMOPARAMS(6.735837e-01, 2.7255, 0.0494142797907188, 0.31242079216478097, 0., -1, 0, 3.046, 1.0, 0.06, 0., 0., 0.245, 1., 1.)

    _, xe, Tm = pyhy.call_run_hyrec(cosmo_params, inj_params)

    data = np.loadtxt("./tests/output_xe.dat")

    xe_c = np.flip(data[:, 1])
    Tm_c = np.flip(data[:, 2])

    assert all(abs(xe-xe_c)/xe_c < 1e-6) 
    assert all(abs(Tm-Tm_c)/Tm_c < 1e-6) 





