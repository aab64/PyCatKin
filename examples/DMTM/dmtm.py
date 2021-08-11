from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import draw_energy_landscapes, compare_energy_landscapes
from pycatkin.functions.presets import run, run_temperatures, run_energy_span_temperatures
from pycatkin.functions.presets import save_state_energies, save_energies, save_energies_temperatures
import copy
import numpy as np

# Load input file
sim_system = read_from_input_file()

# Draw energy landscape
draw_energy_landscapes(sim_system=sim_system,
                       etype='electronic',
                       show_labels=True,
                       fig_path='figures/')

sim_system.params['temperature'] = 450
draw_energy_landscapes(sim_system=sim_system,
                       fig_path='figures/')

# Compare energy landscapes of two systems
sim_system2 = copy.deepcopy(sim_system)
sim_system2.params['temperature'] = 650

sim_systems = {'450 K': sim_system,
               '650 K': sim_system2}

compare_energy_landscapes(sim_systems=sim_systems,
                          legend_location='upper right',
                          show_labels=True,
                          fig_path='figures/')

# Run microkinetics simulation
run(sim_system=sim_system,
    plot_results=True,
    save_results=True,
    fig_path='figures/',
    csv_path='outputs/')

# Run simulations for a range of temperatures
temperatures = np.linspace(start=400, stop=800, num=17, endpoint=True)
run_temperatures(sim_system=sim_system,
                 temperatures=temperatures,
                 steady_state_solve=True,
                 plot_results=True,
                 save_results=True,
                 fig_path='figures/',
                 csv_path='outputs/')

# Run simulations and compute DRC
tof_terms = ['r5', 'r9']
run_temperatures(sim_system=sim_system,
                 temperatures=temperatures,
                 tof_terms=tof_terms,
                 steady_state_solve=True,
                 plot_results=True,
                 save_results=True,
                 fig_path='figures/',
                 csv_path='outputs/')

# Run energy span model for a range of temperatures
run_energy_span_temperatures(sim_system=sim_system,
                             temperatures=temperatures,
                             save_results=True,
                             csv_path='outputs/')

# Save state energies
save_state_energies(sim_system=sim_system,
                    csv_path='outputs/')

# Save reaction energies and barriers
save_energies(sim_system=sim_system,
              csv_path='outputs/')

# Save reaction energies and barriers for a range of temperatures
save_energies_temperatures(sim_system=sim_system,
                           temperatures=temperatures,
                           csv_path='outputs/')
