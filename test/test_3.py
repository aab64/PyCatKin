from pkg_resources import resource_filename, Requirement
from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import run_temperatures
import os
import pandas as pd


def test_3(tmpdir):
    """Regression test for reactor example.

    """
    print('Regression test for reactor example.')
    print('------------------------------------')

    print('(1/4) Loading input file')
    infile = resource_filename(Requirement.parse("pycatkin"), 'examples/COOxReactor/input_Pd111.json')
    assert(os.path.isfile(infile))
    sim_system = read_from_input_file(input_path=infile)

    print('(2/4) Reconfiguring paths for test')
    for s in sim_system.snames:
        if sim_system.states[s].path:
            sim_system.states[s].path = resource_filename(Requirement.parse("pycatkin"),
                                                          'examples/COOxReactor/' +
                                                          sim_system.states[s].path)

    print('(3/4) Drawing states with ASE')
    for s in sim_system.snames:
        if sim_system.states[s].state_type != 'TS':
            sim_system.states[s].save_pdb(path=tmpdir)

    print('(4/4) Running ODE and SS solvers')
    run_temperatures(sim_system=sim_system,
                     temperatures=[523],
                     steady_state_solve=True,
                     save_results=True,
                     csv_path=tmpdir)
    assert(os.path.isfile(tmpdir + 'pressures_vs_temperature.csv'))
    df = pd.read_csv(filepath_or_buffer=tmpdir + 'pressures_vs_temperature.csv')
    pCOin = sim_system.params['inflow_state']['CO']
    pCOout = df['pCO (bar)'].values
    xCO = 100.0 * (1.0 - pCOout / pCOin)
    assert(abs(xCO - 51.143) <= 1e-3)
