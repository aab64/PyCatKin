from pkg_resources import resource_filename, Requirement
from pycatkin.functions.load_input import read_from_input_file
import os
import numpy as np


def test_2(tmpdir):
    """Regression test for volcano example.

    """
    print('Regression test for volcano example.')
    print('------------------------------------')

    print('(1/2) Loading input file')
    infile = resource_filename(Requirement.parse("pycatkin"), 'examples/COOxVolcano/input.json')
    assert(os.path.isfile(infile))
    sim_system = read_from_input_file(input_path=infile)

    # Define binding energies
    ECO = -1.0
    EO = -1.0

    # Standard entropies (taken from Atkins, in eV/K)
    SCOg = 2.0487e-3
    SO2g = 2.1261e-3

    # Note the temperature and pressure
    T = sim_system.params['temperature']

    # (a) Set CO adsorption energy and entropy
    sim_system.reactions['CO_ads'].dErxn_user = ECO
    sim_system.reactions['CO_ads'].dGrxn_user = ECO + SCOg * T

    # (b) Set O adsorption energy and entropy
    sim_system.reactions['2O_ads'].dErxn_user = 2.0 * EO
    sim_system.reactions['2O_ads'].dGrxn_user = 2.0 * EO + SO2g * T

    # (c) Set O2 adsorption energy and entropy
    EO2 = sim_system.states['sO2'].get_potential_energy()
    sim_system.reactions['O2_ads'].dErxn_user = EO2
    sim_system.reactions['O2_ads'].dGrxn_user = EO2 + SO2g * T

    # (d) Set CO oxidation barrier
    ETS_CO_ox = sim_system.states['SRTS_ox'].get_potential_energy()
    sim_system.reactions['CO_ox'].dEa_fwd_user = np.max((ETS_CO_ox - (ECO + EO), 0.0))

    # (e) Set O2 dissociation barrier
    ETS_O2_2O = sim_system.states['SRTS_O2'].get_potential_energy()
    sim_system.reactions['O2_2O'].dEa_fwd_user = np.max((ETS_O2_2O - EO2, 0.0))

    print('(2/2) Computing the activity')
    activity = sim_system.activity(tof_terms=['CO_ox'])
    assert(abs(activity - (-1.563)) <= 1e-3)
