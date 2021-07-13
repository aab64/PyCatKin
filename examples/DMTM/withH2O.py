from pycatkin.functions.load_input import *
from pycatkin.functions.presets import *
from ase.io import read


def load_from_contcar(state_path):
    contcar_path = state_path + '/CONTCAR'
    assert (os.path.isfile(contcar_path))
    atoms = read(contcar_path)
    inertia = atoms.get_moments_of_inertia()
    outcar_path = state_path + '/OUTCAR'
    assert (os.path.isfile(outcar_path))
    atoms = read(outcar_path, format='vasp-out')
    mass = sum(atoms.get_masses())
    return atoms, mass, inertia


# Location of mechanism file
input_file = 'input_withH2O.json'

# Location of outcars and frequencies
metal = 'cu/'
base_data_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Metals/'
ads_opt_dir = base_data_dir + 'energies/withH2O/' + metal
ads_vib_dir = base_data_dir + 'vibrations/withH2O/' + metal
gas_opt_dir = base_data_dir + 'energies/molecules/'
gas_vib_dir = base_data_dir + 'vibrations/molecules/'

# Location of results files and images
base_out_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Metals/'
results_dir = base_out_dir + 'results/withH2O/' + metal
figures_dir = base_out_dir + 'images/withH2O/' + metal

print('---------------------')
print('System: DMTM with H2O')
print('---------------------')

sys = read_from_input_file(input_path=input_file)

for s in sys.snames:
    opt_dir = gas_opt_dir if sys.states[s].state_type == 'gas' else ads_opt_dir
    vib_dir = gas_vib_dir if sys.states[s].state_type == 'gas' else ads_vib_dir
    sys.states[s].path = opt_dir + sys.states[s].name
    sys.states[s].vibs_path = vib_dir + sys.states[s].name
    read_from_alternate = None
    if sys.states[s].state_type == 'gas':
        read_from_alternate = {'get_atoms': lambda state_path=opt_dir + sys.states[s].name: load_from_contcar(
            state_path)}
    sys.states[s].read_from_alternate = read_from_alternate

# draw_states(sim_system=sys)
# run(sim_system=sys, plot_results=True, save_results=False)
run_temperatures(sim_system=sys,
                 steady_state_solve=False,
                 tof_terms=['r5', 'r9'],
                 temperatures=np.linspace(start=400, stop=800, num=20, endpoint=True),
                 plot_results=True, save_results=False)



