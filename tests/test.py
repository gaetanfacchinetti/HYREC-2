import pyhyrec as pyhy

def test_default_output():
    
    inj_params = pyhy.init_INPUT_INJ_PARAMS(0., 0., 0., 0., 0., 0., 0., 0, 0., 1., 0.)
    cosmo_params = pyhy.init_INPUT_COSMOPARAMS(6.735837e-01, 2.7255, 0.0494142797907188, 0.31242079216478097, 0., -1, 0, 3.046, 1.0, 0.06, 0., 0., 0.245, 1., 1.)

    z, xe, Tm = pyhy.call_run_hyrec(cosmo_params, inj_params)
