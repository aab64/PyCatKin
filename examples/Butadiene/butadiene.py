from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import draw_energy_landscapes, compare_energy_landscapes
import matplotlib.pyplot as plt

# Load input file
sim_system = read_from_input_file()

# Draw energy landscapes for p123, p124, p156, BuOH, and doped steps (dehydrogenation, dehydration, aldol, Prins)
# draw_energy_landscapes(sim_system=sim_system,
#                        etype='free',
#                        show_labels=True,
#                        eunits='kcal/mol',
#                        fig_path='all_fes/')

# Compare energy landscapes for doped steps (dehydrogenation, dehydration, aldol, Prins)
for step in ['dehydrogenation', 'dehydration', 'aldol', 'prins']:
    compare_energy_landscapes(sim_system,
                              landscapes=[i for i in sim_system.energy_landscapes.keys() if step in i],
                              etype='electronic',
                              eunits='eV',
                              legend_location='upper left',
                              show_labels=True,
                              fig_path=step + '/',
                              cmap=plt.get_cmap("Dark2"))
