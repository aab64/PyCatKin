from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import draw_states, run_temperatures, plot_data_simple
import os
import numpy as np
import pandas as pd

fig, ax = None, None
if not os.path.isdir('figures'):
    os.mkdir('figures')
if not os.path.isdir('outputs'):
    os.mkdir('outputs')

# Load input files
sim_system_Au = read_from_input_file(input_path='input_AuPd.json')
sim_system_Pd = read_from_input_file(input_path='input_Pd111.json')

# View the states using ASE
draw_states(sim_system=sim_system_Au,
            fig_path='figures/AuPd/')  # rotation='-90x'

# Save the states in proteindatabank (.pdb) format
for s in sim_system_Pd.snames:
    if sim_system_Pd.states[s].state_type != 'TS':
        sim_system_Pd.states[s].save_pdb(path='figures/Pd111/')

# Run simulations for a range of temperatures
temperatures = np.linspace(start=423, stop=623, num=20, endpoint=True)
for sysname, sim_system in [['AuPd', sim_system_Au], ['Pd111', sim_system_Pd]]:
    run_temperatures(sim_system=sim_system,
                     temperatures=temperatures,
                     steady_state_solve=True,
                     plot_results=False,
                     save_results=True,
                     fig_path='figures/%s/' % sysname,
                     csv_path='outputs/%s/' % sysname)

    # Read in one of the output files
    df = pd.read_csv(filepath_or_buffer='outputs/%s/pressures_vs_temperature.csv' % sysname)

    # Compute conversion across the reactor
    pCOin = sim_system_Pd.params['inflow_state']['CO']
    pCOout = df['pCO (bar)'].values
    xCO = 100.0 * (1.0 - pCOout / pCOin)

    fig, ax = plot_data_simple(fig=fig,
                               ax=ax,
                               xdata=temperatures,
                               ydata=xCO,
                               xlabel='Temperature (K)',
                               ylabel='Conversion (%)',
                               label=sysname,
                               addlegend=True,
                               color='teal' if sysname == 'Pd111' else 'salmon',
                               fig_path='figures/',
                               fig_name='conversion')
