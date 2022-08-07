from pkg_resources import resource_filename, Requirement
from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import run, run_temperatures, run_energy_span_temperatures
from pycatkin.functions.presets import save_state_energies, save_energies, save_energies_temperatures
import os
import numpy as np
import pandas as pd


def test_1(tmpdir):
    """Regression test for DMTM examples.

    """
    print('Regression test for DMTM examples.')
    print('---------------------------------')

    print('(1/7) Loading input file')
    infile = resource_filename(Requirement.parse("pycatkin"), 'examples/DMTM/input.json')
    assert(os.path.isfile(infile))
    sim_system = read_from_input_file(input_path=infile)

    print('(2/7) Reconfiguring paths for test')
    for s in sim_system.snames:
        if sim_system.states[s].path:
            sim_system.states[s].path = resource_filename(Requirement.parse("pycatkin"),
                                                          'examples/DMTM/' +
                                                          sim_system.states[s].path)
        if sim_system.states[s].vibs_path:
            sim_system.states[s].vibs_path = resource_filename(Requirement.parse("pycatkin"),
                                                               'examples/DMTM/' +
                                                               sim_system.states[s].vibs_path)

    print('(3/7) Checking system loaded correctly')
    assert(sim_system is not None)
    assert(sim_system.reactions is not None)
    assert(sim_system.states is not None)
    assert(sim_system.reactor is not None)
    assert(sim_system.energy_landscapes is not None)

    print('(4/7) Running ODE solver')
    run(sim_system=sim_system)
    assert(abs(1 - np.sum(sim_system.solution[-1][sim_system.adsorbate_indices])) <= 1e-6)
    assert(np.max(sim_system.solution[-1][sim_system.adsorbate_indices]) > 0.999)
    assert(sim_system.snames[[i for i in sim_system.adsorbate_indices
                              if sim_system.solution[-1][i] ==
                              np.max(sim_system.solution[-1][sim_system.adsorbate_indices])][0]] == 'sCH3OH')

    print('(5/7) Running DRC calculations')
    tof_terms = ['r5', 'r9']
    temperatures = np.linspace(start=400, stop=800, num=2, endpoint=True)
    run_temperatures(sim_system=sim_system,
                     temperatures=temperatures,
                     tof_terms=tof_terms,
                     steady_state_solve=True,
                     save_results=True,
                     csv_path=tmpdir)
    assert(os.path.isfile(tmpdir + 'drcs_vs_temperature.csv'))
    df = pd.read_csv(filepath_or_buffer=tmpdir + 'drcs_vs_temperature.csv')
    assert([i for i in df.columns[1::] if df[i][0] == max(df.T[0][1::])][0] == 'r9')

    print('(6/7) Running energy span calculations')
    run_energy_span_temperatures(sim_system=sim_system,
                                 temperatures=temperatures,
                                 save_results=True,
                                 csv_path=tmpdir)
    assert(os.path.isfile(tmpdir + 'energy_span_summary_full_pes.csv'))
    df = pd.read_csv(filepath_or_buffer=tmpdir + 'energy_span_summary_full_pes.csv')
    assert(df['TDI'][0] == 'sCH3OH')
    assert(df['TDI'][1] == 's2OCH4')
    assert(df['TDTS'][0] == 'TS6')
    assert(df['TDTS'][1] == 'TS3')

    print('(7/7) Saving energies')
    save_state_energies(sim_system=sim_system,
                        csv_path=tmpdir)
    assert(os.path.isfile(tmpdir + 'state_energies_800.0K_1.0bar.csv'))
    df = pd.read_csv(filepath_or_buffer=tmpdir + 'state_energies_800.0K_1.0bar.csv')
    assert(abs(max(df['Free (eV)']) - (-7.864)) <= 1e-3)
    assert(abs(max(df['Vibrational (eV)']) - 1.142) <= 1e-3)
    assert(abs(min(df['Rotational (eV)']) - (-1.259)) <= 1e-3)
    assert(abs(min(df['Translational (eV)']) - (-0.659)) <= 1e-3)

    save_energies(sim_system=sim_system,
                  csv_path=tmpdir)
    assert(os.path.isfile(tmpdir + 'reaction_energies_and_barriers_800.0K_1.0bar.csv'))
    df = pd.read_csv(filepath_or_buffer=tmpdir + 'reaction_energies_and_barriers_800.0K_1.0bar.csv')
    assert(abs(max(df['dEr (J/mol)']) - 220788.916) <= 1e-3)
    assert(abs(max(df['dGr (J/mol)']) - 66358.978) <= 1e-3)
    assert(abs(max(df['dEa (J/mol)']) - 138934.617) <= 1e-3)
    assert(abs(max(df['dGa (J/mol)']) - 230155.396) <= 1e-3)
