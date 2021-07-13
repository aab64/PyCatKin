from pycatkin.functions.load_input import *
from pycatkin.functions.presets import *

input_file = 'input.json'

print('--------------------')
print('System: CO oxidation')
print('--------------------')

sys = read_from_input_file(input_path=input_file)

# draw_states(sim_system=sys)
run(sim_system=sys, plot_results=True, save_results=False)
# run_temperatures(sim_system=sys,
#                  steady_state_solve=False,
#                  tof_terms=['r5', 'r9'],
#                  temperatures=np.linspace(start=400, stop=800, num=20, endpoint=True),
#                  plot_results=True, save_results=False)
