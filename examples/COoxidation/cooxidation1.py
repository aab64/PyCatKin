from pycatkin.functions.load_input import *

print('--------------------')
print('System: CO oxidation')
print('--------------------')

input_file = 'input.json'
sys = read_from_input_file(input_path=input_file)
sys.solve_odes()
sys.plot_transient()
sys.find_steady(store_steady=True)
sys.reaction_terms(sys.full_steady)

