from pycatkin.classes.state import *
from pycatkin.classes.energy import *
from pycatkin.classes.reaction import *
from pycatkin.classes.system import *
from pycatkin.classes.reactor import *
from pycatkin.functions.profiling import *
import glob
import copy
from ase.data.pubchem import pubchem_atoms_search
from ase.io import read
import pandas as pd


def load_acetaldehyde_from_pubchem(outcar_path):
    aa = pubchem_atoms_search(name='acetaldehyde')
    inertia = copy.copy(aa.get_moments_of_inertia())
    vasp_atoms = read(outcar_path, format='vasp-out')
    mass = sum(vasp_atoms.get_masses())
    return vasp_atoms, mass, inertia


def load_other_frequencies(freq_path):
    if verbose:
        print('Checking OUTCAR for frequencies')
    assert (freq_path is not None)
    freq_path = freq_path.split('\\')[0] + '/frequencies/' + freq_path.split('\\')[1] + '.dat'
    assert (os.path.isfile(freq_path))
    freq = []
    i_freq = []
    with open(freq_path, 'r') as fd:
        lines = fd.readlines()
    for line in lines:
        data = line.split()
        if 'THz' in data:
            fHz = float(data[-4]) * 1.0e12
            if 'f/i=' not in data and 'f/i' not in data:
                freq.append(fHz)
            else:
                i_freq.append(fHz)
    return freq, i_freq


def load_ethyl_acetate_energy(energy_path):
    energy_path = energy_path.split('\\')[0] + '/energies/' + energy_path.split('\\')[1] + '.dat'
    if os.path.isfile(energy_path):
        with open(energy_path, 'r') as fd:
            lines = fd.readlines()
            Gelec_aa = float(lines[0].split()[0])
    return Gelec_aa


# Conditions
p = 1.01325e5  # Pressure (Pa)
T = 723  # Temperature (K)
Ts = list(np.linspace(start=523, stop=923, num=17, endpoint=True))  # Temperatures (K)
xethanol = np.linspace(start=0.002, stop=0.2, num=100)
Asite = (13e-10 * 15e-10) / 4  # Site area (m2)
# times = np.logspace(start=-12, stop=3, num=int(0.5e3)) * 24 * 3.6  # Times (s)
times = pd.read_csv('original_times.csv').values[1::, 0]
use_jacobian = True  # Use Jacobian to solve SS and ODEs
verbose = False  # Print messages
rtol = 1.0e-6
atol = 1.0e-8
xtol = 1.0e-16
ftol = 1.0e-16

runMK = True
runES = True
runSweep = False
saveResults = True
profile = False

# Location of outcars and frequencies
adsdir = 'D:/Users/Astrid/Documents/Chalmers/Data/Butadiene/DFT data/'
folders = ['dehydrogenation', 'aldol condensation',
           'MPV red', 'aldol mpv',
           'dehydration', 'prins', 'etoxetoh',
           'adsorption', 'H2', 'gas',
           'butanol formation', 'ethyl acetate']

# Location of results files and images
results_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Butadiene/Results/'
figures_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Butadiene/Images/'

print('-----------------')
print('System: butadiene')
print('-----------------')

# Load DFT data
print('Configuring states from DFT states...')

states = []
for folder in folders:
    paths = glob.glob(adsdir + folder + '/*')
    for path in paths:
        if not os.path.isdir(path):
            if 'ethyl acetate' in path and '7Eiv' not in path:
                read_from_alternate = dict()
                read_from_alternate['get_electronic_energy'] = lambda ea_path=path: load_ethyl_acetate_energy(ea_path)
                read_from_alternate['get_vibrations'] = lambda ea_path=path: load_other_frequencies(ea_path)
            elif 'acetaldehyde' in path:
                read_from_alternate = dict()
                read_from_alternate['get_atoms'] = lambda aa_path=path: load_acetaldehyde_from_pubchem(aa_path)
            elif 'butanol formation' in path:
                read_from_alternate = dict()
                read_from_alternate['get_vibrations'] = lambda bu_path=path: load_other_frequencies(bu_path)
            else:
                read_from_alternate = None
            if 'gas' in path:
                if 'H2' in path or 'H2O' in path:
                    states += [State(state_type='gas', path=path, sigma=2,
                                     read_from_alternate=read_from_alternate, truncate_freq=False)]
                else:
                    states += [State(state_type='gas', path=path, sigma=1,
                                     read_from_alternate=read_from_alternate, truncate_freq=False)]
            else:
                states += [State(state_type='adsorbate', path=path,
                                 read_from_alternate=read_from_alternate, truncate_freq=False)]

# Zip states into a dictionary
snames = [s.name for s in states]
states = dict(zip(snames, states))

print('Done.')

# for s in ['3Ci', '3Cii', '3Ciii']:
#     states[s].view_atoms(path=figures_dir)
#     states[s].save_pdb(path=figures_dir)
#     states[s].calc_electronic_energy(verbose=verbose)
#
# for k in ['7Eiv']:  # states.keys()
#     if states[k].state_type == 'adsorbate':
#         for i in range(len(states[k].atoms)):
#             if states[k].atoms.positions[i, 0] > 12:
#                 states[k].atoms.positions[i, :] -= states[k].atoms.get_cell()[0, :]
#             if states[k].atoms.positions[i, -1] > 10:
#                 states[k].atoms.positions[i, :] -= states[k].atoms.get_cell()[-1, :]
#     states[k].view_atoms(path=figures_dir + 'ASE/')
#     states[k].save_pdb(path=figures_dir + 'xyz/')

for s in snames:
    states[s].get_free_energy(T=T, p=p)
    if states[s].i_freq is not None:
        if len(states[s].i_freq) == 1:
            # print(s)
            states[s].state_type = 'TS'

# Base reactions
reactions = dict()
initial_states = ['1A', '2A', '2F', '2J', '2L', '3A', '3D', '3F',
                  '4A', '4D', '4D', '4F', '4I', '5A',
                  '6A', '6C', '6E', '6G', '7A', '9D',
                  '8A', '3Ci', '7Ei']  # , '3Civ'
transition_states = ['1B', '2B', '2G', '2K', '2M', '3B', '3E', '3Fi',
                     '4B', None, '4E', '4G', '4J', '5B',
                     '6B', '6D', '6F', '6Gi', '7D', None,
                     '8B', '3Cii', '7Eii']  # , '3Cv'
final_states = ['1C', '2C', '2H', '2L', '2N', '3C', '3F', '3G',
                '4C', '4Ca', '4F', '4H', '4K', '5C',
                '6C', '6E', '6G', '6H', '7E', '9C',
                '8C', '3Ciii', '7Eiii']  # , '3Cvi'
for r in range(len(initial_states)):
    rname = initial_states[r] + '-' + final_states[r]
    reactions[rname] = Reaction(name=rname,
                                area=Asite,
                                reac_type='Arrhenius',
                                reactants=[states[initial_states[r]]],
                                products=[states[final_states[r]]],
                                TS=[states[transition_states[r]]] if transition_states[r] else None)

ads_initial_states = [['surface', 'ethanol'],
                      ['surface', 'H2O'],
                      ['surface', 'acetaldehyde'],
                      ['2O', 'crotonaldehyde'],
                      ['7Eiv', 'ethylacetate']]
ads_final_states = ['1A', '9B', '10B', '2N', '7Eiii']
for r in range(len(ads_initial_states)):
    rname = ads_initial_states[r][1] + '-' + ads_final_states[r]
    reactions[rname] = Reaction(name=rname,
                                area=Asite,
                                reac_type='adsorption',
                                reactants=[states[i] for i in ads_initial_states[r]],
                                products=[states[ads_final_states[r]]])

print('Reaction energies and barriers:')
print('-------------------------------')
for r in reactions.keys():
    print(r + ' & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g' %
          (reactions[r].get_reaction_barriers(T=623, p=p)[0] / kcaltoJ,
           reactions[r].get_reaction_barriers(T=723, p=p)[0] / kcaltoJ,
           reactions[r].get_reaction_barriers(T=823, p=p)[0] / kcaltoJ,
           reactions[r].get_reaction_barriers(T=623, p=p)[1] / kcaltoJ,
           reactions[r].get_reaction_barriers(T=723, p=p)[1] / kcaltoJ,
           reactions[r].get_reaction_barriers(T=823, p=p)[1] / kcaltoJ,
           reactions[r].get_reaction_energy(T=623, p=p) / kcaltoJ,
           reactions[r].get_reaction_energy(T=723, p=p) / kcaltoJ,
           reactions[r].get_reaction_energy(T=823, p=p) / kcaltoJ))

ref_free = states['surface'].get_free_energy(T=723, p=p) + 3.0 * states['ethanol'].get_free_energy(T=723, p=p)
ref_elec = states['surface'].Gelec + 3.0 * states['ethanol'].Gelec

newstates = ['3F', '3Fi', '3G',
             '6G', '6Gi', '6H',
             '4Ca',
             '3Ci', '3Cii', '3Ciii',
             # '3Civ', '3Cv', '3Cvi', '3Cvii',
             '7Ei', '7Eii', '7Eiii', '7Eiv',
             '8A', '8B', '8C',
             '9A', '9B', '9C', '9D',
             '10A', '10B']

gasstates = [['H2', 'H2', 'H2O', 'acetaldehyde'],  # 3F
             ['H2', 'H2', 'H2O', 'acetaldehyde'],  # 3Fi
             ['H2', 'H2', 'H2O', 'acetaldehyde'],  # 3G
             ['ethanol', 'H2', 'H2O'],  # 6G
             ['ethanol', 'H2', 'H2O'],  # 6Gi
             ['ethanol', 'H2', 'H2O'],  # 6H
             ['H2', 'H2', 'acetaldehyde'],  # 4Ca
             ['H2', 'acetaldehyde'],  # 3Ci
             ['H2', 'acetaldehyde'],  # 3Cii
             ['H2', 'acetaldehyde'],  # 3Ciii
             # ['H2', 'H2O', 'acetaldehyde'],  # 3Civ
             # ['H2', 'H2O', 'acetaldehyde'],  # 3Cv
             # ['H2', 'H2O', 'acetaldehyde'],  # 3Cvi
             # ['H2', 'H2O', 'acetaldehyde', 'butanol'],  # 3Cvi
             ['ethanol', 'H2'],  # 7Ei
             ['ethanol', 'H2'],  # 7Eii
             ['ethanol', 'H2'],  # 7Eiii
             ['ethanol', 'H2', 'ethylacetate'],  # 7Eiv
             ['ethanol', 'ethanol', 'ethanol'],  # 8A
             ['ethanol', 'ethanol', 'ethanol'],  # 8B
             ['ethanol', 'ethanol', 'ethanol'],  # 8C
             ['ethanol', 'ethanol', 'ethanol', 'H2O'],  # 9A
             ['ethanol', 'ethanol', 'ethanol'],  # 9B
             ['ethanol', 'ethanol', 'ethanol'],  # 9C
             ['ethanol', 'ethanol', 'ethanol'],  # 9D
             ['ethanol', 'ethanol', 'ethanol', 'acetaldehyde'],  # 10A
             ['ethanol', 'ethanol', 'ethanol']]  # 10B

print('State energies:')
print('---------------')
for i, s in enumerate(newstates):
    Gfree = np.sum([states[j].get_free_energy(T=T, p=p) for j in gasstates[i]]) - ref_free
    Gelec = np.sum([states[j].Gelec for j in gasstates[i]]) - ref_elec
    if s not in ['9A', '10A']:
        Gfree += states[s].get_free_energy(T=T, p=p)
        Gelec += states[s].Gelec
    else:
        Gfree += states['surface'].get_free_energy(T=T, p=p)
        Gelec += states['surface'].Gelec
    print(s + ' & %.3g & %.3g' %
          (Gelec * eVtokJ * 1.0e3 / kcaltoJ,
           Gfree * eVtokJ * 1.0e3 / kcaltoJ))

# Energy landscapes
# =================

p123 = dict()
p123[0] = [states['surface'], states['ethanol'], states['ethanol'], states['ethanol']]
p123[1] = [states['1A'], states['ethanol'], states['ethanol']]
p123[2] = [states['1B'], states['ethanol'], states['ethanol']]
p123[3] = [states['1C'], states['ethanol'], states['ethanol']]
p123[4] = [states['2A'], states['ethanol'], states['ethanol'], states['H2']]
p123[5] = [states['2B'], states['ethanol'], states['ethanol'], states['H2']]
p123[6] = [states['2C'], states['ethanol'], states['ethanol'], states['H2']]
p123[7] = [states['2D'], states['ethanol'], states['ethanol'], states['H2']]
p123[8] = [states['2E'], states['ethanol'], states['ethanol'], states['H2']]
p123[9] = [states['2F'], states['ethanol'], states['H2'], states['H2']]
p123[10] = [states['2G'], states['ethanol'], states['H2'], states['H2']]
p123[11] = [states['2H'], states['ethanol'], states['H2'], states['H2']]
p123[12] = [states['2I'], states['ethanol'], states['H2'], states['H2']]
p123[13] = [states['2J'], states['ethanol'], states['H2'], states['H2']]
p123[14] = [states['2K'], states['ethanol'], states['H2'], states['H2']]
p123[15] = [states['2L'], states['ethanol'], states['H2'], states['H2']]
p123[16] = [states['2M'], states['ethanol'], states['H2'], states['H2']]
p123[17] = [states['2N'], states['ethanol'], states['H2'], states['H2']]
p123[18] = [states['3A'], states['H2'], states['H2'], states['H2O']]
p123[19] = [states['3B'], states['H2'], states['H2'], states['H2O']]
p123[20] = [states['3C'], states['H2'], states['H2'], states['H2O']]
p123[21] = [states['3D'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
p123[22] = [states['3E'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
p123[23] = [states['3F'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
p123[24] = [states['3Fi'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
p123[25] = [states['3G'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
l123 = [i[0].name if i[0].name != 'surface' else '' for i in p123.values()]
energy_p123 = Energy(minima=p123, labels=l123, name='p123')

p124 = dict()
p124[0] = [states['surface'], states['ethanol'], states['ethanol'], states['ethanol']]
p124[1] = [states['1A'], states['ethanol'], states['ethanol']]
p124[2] = [states['1B'], states['ethanol'], states['ethanol']]
p124[3] = [states['1C'], states['ethanol'], states['ethanol']]
p124[4] = [states['2A'], states['ethanol'], states['ethanol'], states['H2']]
p124[5] = [states['2B'], states['ethanol'], states['ethanol'], states['H2']]
p124[6] = [states['2C'], states['ethanol'], states['ethanol'], states['H2']]
p124[7] = [states['2D'], states['ethanol'], states['ethanol'], states['H2']]
p124[8] = [states['2E'], states['ethanol'], states['ethanol'], states['H2']]
p124[9] = [states['2F'], states['ethanol'], states['H2'], states['H2']]
p124[10] = [states['2G'], states['ethanol'], states['H2'], states['H2']]
p124[11] = [states['2H'], states['ethanol'], states['H2'], states['H2']]
p124[12] = [states['2I'], states['ethanol'], states['H2'], states['H2']]
p124[13] = [states['2J'], states['ethanol'], states['H2'], states['H2']]
p124[14] = [states['4A'], states['H2'], states['H2']]
p124[15] = [states['4B'], states['H2'], states['H2']]
p124[16] = [states['4C'], states['H2'], states['H2']]
p124[17] = [states['4D'], states['H2'], states['H2'], states['acetaldehyde']]
p124[18] = [states['4E'], states['H2'], states['H2'], states['acetaldehyde']]
p124[19] = [states['4F'], states['H2'], states['H2'], states['acetaldehyde']]
p124[20] = [states['4G'], states['H2'], states['H2'], states['acetaldehyde']]
p124[21] = [states['4H'], states['H2'], states['H2'], states['acetaldehyde']]
p124[22] = [states['4I'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
p124[23] = [states['4J'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
p124[24] = [states['4K'], states['H2'], states['H2'], states['H2O'], states['acetaldehyde']]
l124 = [i[0].name if i[0].name != 'surface' else '' for i in p124.values()]
energy_p124 = Energy(minima=p124, labels=l124, name='p124')

p156 = dict()
p156[0] = [states['surface'], states['surface'], states['ethanol'], states['ethanol'], states['ethanol']]
p156[1] = [states['1A'], states['ethanol'], states['ethanol'], states['surface']]
p156[2] = [states['1B'], states['ethanol'], states['ethanol'], states['surface']]
p156[3] = [states['1C'], states['ethanol'], states['ethanol'], states['surface']]
p156[4] = [states['5A'], states['1C'], states['ethanol']]
p156[5] = [states['5B'], states['1C'], states['ethanol']]
p156[6] = [states['5C'], states['1C'], states['ethanol']]
p156[7] = [states['6A'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[8] = [states['6B'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[9] = [states['6C'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[10] = [states['6D'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[11] = [states['6E'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[12] = [states['6F'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[13] = [states['6G'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[14] = [states['6Gi'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[15] = [states['6H'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
l156 = [i[0].name if i[0].name != 'surface' else '' for i in p156.values()]
energy_p156 = Energy(minima=p156, labels=l156, name='p156')

pBuOH = dict()
pBuOH[0] = [states['surface'], states['ethanol'], states['ethanol'], states['ethanol']]
pBuOH[1] = [states['3C'], states['H2'], states['H2'], states['H2O']]
pBuOH[2] = [states['3Ci'], states['H2'], states['acetaldehyde']]
pBuOH[3] = [states['3Cii'], states['H2'], states['acetaldehyde']]
pBuOH[4] = [states['3Ciii'], states['H2'], states['acetaldehyde']]
lBuOH = [i[0].name if i[0].name != 'surface' else '' for i in pBuOH.values()]
energy_pBuOH = Energy(minima=pBuOH, labels=lBuOH, name='p3b')

# MK model
# ========
gas_states = ['H2', 'H2O', 'C2H4', 'C4H6', 'C2H4O', 'C2H5OH', 'C4H6O', 'C4H10O', 'C4H8O2']
gas_masses = [2.016, 18.02, 28.05, 54.09, 44.05, 46.07, 70.09, 74.12, 88.11]
gas_fracts = [0.01, 1.0e-8, 0.01, 1.0e-8, 1.0e-8, 0.02, 1.0e-8, 1.0e-8, 1.0e-8]
gas_sigmas = [2 if i == 'H2' or i == 'H2O' else 1 for i in gas_states]
ads_states = ['H*', 'O*', 'OH*',
              'C2H3O*', 'C2H4O*', 'C2H5O*',
              'C4H6O_1*', 'C4H6O_2*', 'C4H6O_3*',
              'C4H7O_1*', 'C4H7O_2*', 'C4H7O_3*',
              'C4H8O_1*',
              'C4H6O2*', 'C4H7O2*',
              'C4H8O2_1*', 'C4H8O2_2*',
              'C4H9O2_1*', 'C4H9O2_2*',
              'C4H8O_2*', 'C4H8O_3*',
              'C4H9O*',
              'C4H8O2_3*',
              'C4H10O*']
surf_state = ['*']

# MK states
mk_states = dict()
for gind, g in enumerate(gas_states):
    mk_states[g] = State(name=g, state_type='gas',
                         sigma=gas_sigmas[gind], mass=gas_masses[gind])
for aind, a in enumerate(ads_states):
    mk_states[a] = State(name=a, state_type='adsorbate')
mk_states['*'] = State(name='*', state_type='surface')

# MK reactions
mk_reactants = [['C2H5O*', 'H*'],
                ['C2H4O*', '*'],
                ['C2H3O*', 'C2H4O*'],
                ['C4H7O2*', '*'],
                ['C4H6O2*', '*'],
                ['C4H6O_1*', 'C2H5O*'],
                ['C4H7O_1*', '*'],
                ['C4H6O_2*'],
                ['C2H5O*', 'C4H7O2*'],
                ['C4H9O2_1*', '*'],  # ['C4H8O2_1*', 'H*'],
                ['C4H9O2_1*', '*'],
                ['C4H8O2_2*', '*'],
                ['C4H7O_2*', '*'],
                ['C2H5O*', 'H*'],
                ['C2H4', 'C2H4O*', '*'],
                ['C4H8O_1*', '*'],
                ['C4H7O_3*', '*'],
                ['C4H6O_3*'],
                ['C2H4O*', 'C2H5O*'],
                ['O*', 'H*'],  # ['OH*', '*'],
                ['C2H5OH', '*', '*'],
                ['H2', '*', '*'],
                ['H2O', '*', '*'],
                ['C2H4O', '*'],
                ['C4H6O', '*'],
                ['C4H7O_1*', 'H*'],
                ['C4H8O_2*'],
                ['C4H8O_3*', 'H*'],
                ['C4H10O', '*', '*'],
                # ['C4H7O_1*', 'H*', 'H*'],
                # ['C4H9O*', 'H*'],
                # ['C4H10O', '*'],
                ['C4H9O2_2*', '*'],
                ['C4H8O2', '*']
                ]

mk_products = [['H2', 'C2H4O*', '*'],
               ['C2H3O*', 'H*'],
               ['C4H7O2*', '*'],
               ['C4H6O2*', 'H*'],
               ['C4H6O_1*', 'O*'],
               ['C4H7O_1*', 'C2H4O*'],
               ['C4H6O_2*', 'H*'],
               ['C4H6', 'O*'],
               ['C4H8O2_1*', 'C2H4O*'],
               ['C4H8O2_1*', 'H*'],  # ['C4H9O2_1*', '*'],
               ['C4H8O2_2*', 'H*'],
               ['C4H7O_2*', 'OH*'],
               ['C4H6', 'O*', 'H*'],
               ['C2H4', 'OH*', 'H*'],
               ['C4H8O_1*', '*'],
               ['C4H7O_3*', 'H*'],
               ['C4H6O_3*', 'H*'],
               ['C4H6', 'O*'],
               ['C4H9O2_2*', '*'],
               ['OH*', '*'],  # ['O*', 'H*'],
               ['C2H5O*', 'H*'],
               ['H*', 'H*'],
               ['OH*', 'H*'],
               ['C2H4O*'],
               ['C4H6O_1*'],
               ['C4H8O_2*', '*'],
               ['C4H8O_3*'],
               ['C4H9O*', '*'],
               ['C4H9O*', 'H*'],
               # ['C4H9O*', '*', '*'],
               # ['C4H10O*', '*'],
               # ['C4H10O*'],
               ['C4H8O2_3*', 'H*'],
               ['C4H8O2_3*']
               ]

mk_rnames = ['1A-1C',
             '2A-2C',
             '2F-2H',
             '2J-2L',
             '2L-2N',
             '3A-3C',
             '3D-3F',
             '3F-3G',
             '4A-4C',
             '4D-4Ca',
             '4D-4F',
             '4F-4H',
             '4I-4K',
             '5A-5C',
             '6A-6C',
             '6C-6E',
             '6E-6G',
             '6G-6H',
             '7A-7E',
             '9D-9C',
             'ethanol-1A',
             '8A-8C',
             'H2O-9B',
             'acetaldehyde-10B',
             'crotonaldehyde-2N',
             '12A-12C',
             '12C-12E',
             '12E-12G',
             'butanol-12G',
             # '3Ci-3Ciii',
             # '3Civ-3Cvi',
             # 'butanol-3Cvi',
             '7Ei-7Eiii',
             'ethylacetate-7Eiii'
             ]


# User defined reactions (free energies at 673 K from Souza et al. (2020)
souza_reactions = ['12A-12C', '12C-12E', '12E-12G', 'butanol-12G']
souza_energies = [[0.61, 0.71, -0.1], [1.06, 1.11, -0.05], [0.20, 0.89, -0.69], [None, None, 0.46]]

mk_reactions = dict()
for r in range(len(mk_rnames)):
    rname = mk_rnames[r]
    if rname not in souza_reactions:
        rtype = reactions[rname].reac_type if (rname != '8A-8C' and rname != '14A-14C') else 'adsorption'
        mk_reactions[rname] = ReactionDerivedReaction(name=rname,
                                                      area=Asite,
                                                      reac_type=rtype,
                                                      base_reaction=reactions[rname],
                                                      reactants=[mk_states[i] for i in mk_reactants[r]],
                                                      products=[mk_states[i] for i in mk_products[r]])
    else:
        rtype = 'Arrhenius' if r != 'butanol-12G' else 'adsorption'
        mk_reactions[rname] = UserDefinedReaction(reac_type=rtype,
                                                  reactants=[mk_states[i]
                                                             for i in mk_reactants[mk_rnames.index(rname)]],
                                                  products=[mk_states[i]
                                                            for i in mk_products[mk_rnames.index(rname)]],
                                                  area=Asite,
                                                  name=rname,
                                                  dGrxn_user=souza_energies[souza_reactions.index(rname)][2],
                                                  dGa_fwd_user=souza_energies[souza_reactions.index(rname)][0])

# Run MK model
mk_results = dict()
mk_sweep = dict()
if runMK:
    print('Solving MK model...')
    mk_tests = ['all4',
                'p123', 'p124', 'p156',
                'byproducts_Souza',
                # 'byproducts_House'
                ]
    test_inds = [list(np.arange(0, 25)),
                 list(np.arange(0, 8)) + list(np.arange(19, 25)),
                 list(np.arange(0, 3)) + list(np.arange(8, 13)) + list(np.arange(19, 25)),
                 [0] + list(np.arange(13, 18)) + list(np.arange(19, 25)),
                 list(np.arange(0, 29))  # + list(np.arange(32, 34)),
                 # list(np.arange(0, 25)) + list(np.arange(29, 34))
                 ]
    mk_reactor = InfiniteDilutionReactor()

    for ind, test_case in enumerate(mk_tests):
        print(test_case)
        mk_results[test_case] = dict()
        test_rxns = dict(zip([list(mk_reactions.keys())[i] for i in test_inds[ind]],
                             [list(mk_reactions.values())[i] for i in test_inds[ind]]))
        mk_sys = System(reactor=mk_reactor, reactions=test_rxns)

        for Ti in Ts:
            print('%1.0f K' % Ti)
            start_state = dict()
            for gind, g in enumerate(gas_states):
                if g in mk_sys.snames:
                    if g != 'C2H5OH':
                        start_state[g] = gas_fracts[gind] * 0.02 * p / bartoPa
                    else:
                        start_state[g] = gas_fracts[gind] * p / bartoPa
            start_state['*'] = 1.0

            mk_sys.set_parameters(times=times, start_state=start_state, T=Ti, p=p,
                                  use_jacobian=use_jacobian, verbose=True,
                                  atol=atol, rtol=rtol, xtol=None, ftol=1.0e-8)
            mk_sys.solve_odes()
            mk_sys.reaction_terms(y=mk_sys.solution[-1])

            tof = 0.0
            for rterm in ['3F-3G', '4I-4K', '6G-6H']:
                if rterm in mk_sys.reactions.keys():
                    i = list(mk_sys.reactions.keys()).index(rterm)
                    tof += mk_sys.rates[i, 0] - mk_sys.rates[i, 1]
            bu_tof = 0.0
            for rterm in ['butanol-3Cvi', 'butanol-12G']:
                if rterm in mk_sys.reactions.keys():
                    i = list(mk_sys.reactions.keys()).index(rterm)
                    bu_tof += mk_sys.rates[i, 1] - mk_sys.rates[i, 0]
            ea_tof = 0.0
            if 'ethylacetate-7Eiii' in mk_sys.reactions.keys():
                i = list(mk_sys.reactions.keys()).index(rterm)
                ea_tof += mk_sys.rates[i, 1] - mk_sys.rates[i, 0]

            rates = np.zeros(len(mk_reactions))
            rates[test_inds[ind]] = mk_sys.rates[:, 0] - mk_sys.rates[:, 1]

            transient_rates = np.zeros((len(mk_sys.times), 2))
            transient_rates[:, 0] = mk_sys.times
            for t in range(len(mk_sys.times)):
                mk_sys.reaction_terms(y=mk_sys.solution[t])
                for rterm in ['3F-3G', '4I-4K', '6G-6H']:
                    if rterm in mk_sys.reactions.keys():
                        i = list(mk_sys.reactions.keys()).index(rterm)
                        transient_rates[t, 1] += mk_sys.rates[i, 0] - mk_sys.rates[i, 1]

            odesol = np.zeros((len(mk_sys.times) + 1, len(ads_states) + 1))
            odesol[0:-1, 0] = mk_sys.times
            odesol[-1, 0] = np.inf
            final_fval = np.zeros(len(ads_states))
            final_diff = np.zeros(len(ads_states))

            mk_sys.find_steady(store_steady=True)
            # mk_sys.reaction_terms(y=mk_sys.full_steady)

            for sind, s in enumerate(ads_states):
                if s in mk_sys.snames:
                    odesol[0:-1, sind + 1] = mk_sys.solution[:, mk_sys.snames.index(s)]
                    odesol[-1, sind + 1] = mk_sys.full_steady[mk_sys.snames.index(s)]
                    final_fval[sind] = mk_sys.species_odes(y=mk_sys.full_steady)[mk_sys.snames.index(s)]
                    final_diff[sind] = odesol[-2, sind + 1] - odesol[-1, sind + 1]

            mk_results[test_case][Ti] = {'tof': tof, 'rates': rates, 'transient_rates': transient_rates,
                                         'cover': odesol[-1, 1::], 'odesol': odesol[0:-1, :],
                                         'final_fval': pow(np.linalg.norm(final_fval), 2),
                                         'final_diff': pow(np.linalg.norm(final_diff), 2),
                                         'ea_tof': ea_tof, 'bu_tof': bu_tof}

            print('||fvals|| = %1.3e, ||xdiff|| = %1.3e' %
                  (np.linalg.norm(final_fval), np.linalg.norm(final_diff)))

            if test_case == 'all4' and runSweep:
                mk_sweep[Ti] = dict()
                for xind, x in enumerate(xethanol):
                    start_state['C2H5OH'] = x * p / bartoPa
                    mk_sys.set_parameters(times=times, start_state=start_state, T=Ti, p=p,
                                          use_jacobian=use_jacobian, verbose=verbose,
                                          atol=atol, rtol=rtol)
                    mk_sys.solve_odes()
                    # mk_sys.find_steady(store_steady=True)
                    # mk_sys.reaction_terms(y=mk_sys.full_steady)
                    mk_sys.reaction_terms(y=mk_sys.solution[-1])
                    mk_sweep[Ti][x] = np.zeros(len(mk_reactions))
                    mk_sweep[Ti][x][test_inds[ind]] = mk_sys.rates[:, 0] - mk_sys.rates[:, 1]
    print('Done.')

# Run ES model
es_results = {'p123': dict(), 'p124': dict(), 'p156': dict()}  # , 'p174': dict()
if runES:
    print('Investigating energy landscapes...')
    energy_p123.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol', path=figures_dir, show_labels=True)
    energy_p124.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol', path=figures_dir, show_labels=True)
    energy_p156.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol', path=figures_dir, show_labels=True)
    energy_pBuOH.draw_energy_landscape(T=T, p=p, etype='electronic', eunits='kcal/mol',
                                       path=figures_dir, show_labels=True, legend_location='lower left')
    energy_pBuOH.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol',
                                       path=figures_dir, show_labels=True, legend_location='upper left')

    for Ti in Ts:
        tof, Espan, TDTS, TDI, xTDTS, xTDI, lTi, lIj = energy_p123.evaluate_energy_span_model(T=Ti, p=p, verbose=False,
                                                                                              etype='free')
        es_results['p123'][Ti] = {'tof': tof, 'Espan': Espan, 'TDTS': TDTS, 'TDI': TDI,
                                  'xTDTS': xTDTS, 'xTDI': xTDI, 'lTi': lTi, 'lIj': lIj}
        tof, Espan, TDTS, TDI, xTDTS, xTDI, lTi, lIj = energy_p124.evaluate_energy_span_model(T=Ti, p=p, verbose=False,
                                                                                              etype='free')
        es_results['p124'][Ti] = {'tof': tof, 'Espan': Espan, 'TDTS': TDTS, 'TDI': TDI,
                                  'xTDTS': xTDTS, 'xTDI': xTDI, 'lTi': lTi, 'lIj': lIj}
        tof, Espan, TDTS, TDI, xTDTS, xTDI, lTi, lIj = energy_p156.evaluate_energy_span_model(T=Ti, p=p, verbose=False,
                                                                                              etype='free')
        es_results['p156'][Ti] = {'tof': tof, 'Espan': Espan, 'TDTS': TDTS, 'TDI': TDI,
                                  'xTDTS': xTDTS, 'xTDI': xTDI, 'lTi': lTi, 'lIj': lIj}
    print('Done.')

if saveResults:
    # MK base case
    mk_dir = results_dir + 'MK/base/'
    ending = '_1perH2_1perC2H4_traceother_2kPaethanol_24hr/'
    alt_mk_dir = results_dir + 'MK/sweep/'
    for k in mk_results.keys():
        file_dir = mk_dir + k + ending
        cover = np.zeros((len(Ts), len(ads_states) + 1))
        rates = np.zeros((len(Ts), len(mk_reactions) + 1))
        ss_vals = np.zeros((len(Ts), 3))
        cover[:, 0] = Ts
        rates[:, 0] = Ts
        ss_vals[:, 0] = Ts
        for Tind, Ti in enumerate(Ts):
            df = pd.DataFrame(mk_results[k][Ti]['odesol'], columns=['Time (s)'] + ads_states)
            df.to_csv(path_or_buf=(file_dir + 'odesol_%1.0fK.csv' % Ti), sep=',', header=True, index=False)
            df = pd.DataFrame(mk_results[k][Ti]['transient_rates'], columns=['Time (s)', 'BD TOF (1/s)'])
            df.to_csv(path_or_buf=(file_dir + 'transient/tof_%1.0fK.csv' % Ti), sep=',', header=True, index=False)
            cover[Tind, 1::] = mk_results[k][Ti]['cover']
            rates[Tind, 1::] = mk_results[k][Ti]['rates']
            ss_vals[Tind, 1] = mk_results[k][Ti]['final_fval']
            ss_vals[Tind, 2] = mk_results[k][Ti]['final_diff']
        df = pd.DataFrame(cover, columns=['Temperature (K)'] + ads_states)
        df.to_csv(path_or_buf=file_dir + 'cover.csv', sep=',', header=True, index=False)
        df = pd.DataFrame(rates, columns=['Temperature (K)'] + list(mk_reactions.keys()))
        df.to_csv(path_or_buf=file_dir + 'rates.csv', sep=',', header=True, index=False)
        df = pd.DataFrame(ss_vals, columns=['Temperature (K)', 'Objective function', 'Norm(xdiff)'])
        df.to_csv(path_or_buf=file_dir + 'steady_states_analysis.csv', sep=',', header=True, index=False,
                  float_format='%1.3e')

        if k == 'all4' and runSweep:
            for xind, x in enumerate(xethanol):
                alt_ending = '_0.02kPaH2_0.02kPaC2H4_traceother_%1.3fkPaethanol_24hr/' % x
                alt_file_dir = alt_mk_dir + k + alt_ending
                if not os.path.isdir(alt_file_dir):
                    os.mkdir(path=alt_file_dir)
                rates = np.zeros((len(Ts), len(mk_reactions) + 1))
                rates[:, 0] = Ts
                for Tind, Ti in enumerate(Ts):
                    rates[Tind, 1::] = mk_sweep[Ti][x]
                df = pd.DataFrame(rates, columns=['Temperature (K)'] + list(mk_reactions.keys()))
                df.to_csv(path_or_buf=alt_file_dir + 'rates.csv', sep=',', header=True, index=False)

    # Energy
    en_dir = results_dir + 'Energy/'
    for Ti in Ts:
        evals = np.zeros((len(mk_reactions), 3))
        for kind, k in enumerate(mk_reactions.keys()):
            evals[kind, 0:-1] = mk_reactions[k].get_reaction_barriers(T=Ti, p=p, verbose=False)
            evals[kind, -1] = mk_reactions[k].get_reaction_energy(T=Ti, p=p, verbose=False)
        evals *= (1.0 / kcaltoJ)
        df = pd.DataFrame(evals, columns=['dGa_f (kcal/mol)', 'dGa_r (kcal/mol)', 'dGr (kcal/mol)'],
                          index=list(mk_reactions.keys()))
        df.to_csv(path_or_buf=(en_dir + 'constants_%1.0fK.csv' % Ti), sep=',', header=True, index=True)

    # ES base case
    es_dir = results_dir + 'ES/base/'
    for p in es_results.keys():
        esones = np.empty((len(Ts), 5), dtype=object)
        esTDTS = np.zeros((len(Ts), len(es_results[p][Ts[0]]['xTDTS']) + 1))
        esTDI = np.zeros((len(Ts), len(es_results[p][Ts[0]]['xTDI']) + 1))

        for Tind, Ti in enumerate(Ts):
            esones[Tind, :] = [str(Ti), str(es_results[p][Ti]['tof']), str(es_results[p][Ti]['Espan']),
                               es_results[p][Ti]['TDTS'], es_results[p][Ti]['TDI']]
            esTDTS[Tind, :] = [Ti] + es_results[p][Ti]['xTDTS']
            esTDI[Tind, :] = [Ti] + es_results[p][Ti]['xTDI']

        df = pd.DataFrame(esones, columns=['Temperature (K)', 'TOF (1/s)', 'Espan (kcal/mol)', 'TDTS', 'TDI'])
        df.to_csv(path_or_buf=es_dir + p + '_summary.csv', sep=',', header=True, index=False)
        df = pd.DataFrame(esTDTS, columns=['Temperature (K)'] + es_results[p][Ts[0]]['lTi'])
        df.to_csv(path_or_buf=es_dir + p + '_xtof_TDTS.csv', sep=',', header=True, index=False)
        df = pd.DataFrame(esTDI, columns=['Temperature (K)'] + es_results[p][Ts[0]]['lIj'])
        df.to_csv(path_or_buf=es_dir + p + '_xtof_TDI.csv', sep=',', header=True, index=False)

# TOF figure
if runES or runMK:
    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    if runES:
        ax.plot(Ts, [es_results['p123'][Ti]['tof'] for Ti in Ts], color='darkorchid', label='p123')
        ax.plot(Ts, [es_results['p124'][Ti]['tof'] for Ti in Ts], color='dodgerblue', label='p124')
        ax.plot(Ts, [es_results['p156'][Ti]['tof'] for Ti in Ts], color='darkorange', label='p156')
    if runMK:
        ax.plot(Ts, [mk_results['all4'][Ti]['tof'] for Ti in Ts], color='black', label='MK')
    ax.set(yscale='log', xlabel='Temperature (K)', ylabel='TOF (1/s)', ylim=(1e-8, 1e2), xlim=(523, 923))
    ax.grid()
    fig.tight_layout()
    fig.show()

if runMK:
    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    for k in ['all4', 'byproducts_Souza']:
        ls = '-' if k == 'all4' else '--' if k == 'byproducts_souza' else ':'
        ax.plot(Ts, [mk_results[k][Ti]['tof'] for Ti in Ts],
                color='black', label=k, linestyle=ls, alpha=0.6)
        ax.plot(Ts, [mk_results[k][Ti]['bu_tof'] for Ti in Ts],
                color='mediumvioletred', label=k, linestyle=ls, alpha=0.6)
        ax.plot(Ts, [mk_results[k][Ti]['ea_tof'] for Ti in Ts],
                color='darkslateblue', label=k, linestyle=ls, alpha=0.6)
    ax.set(xlabel='Temperature (K)', ylabel='TOF (1/s)', yscale='log',
           ylim=(1e-12, 1e0), xlim=(523, 923),
           yticks=np.logspace(start=-12, stop=0, endpoint=True, num=9))
    fig.tight_layout()
    fig.show()

# Profiling for testing purposes (unrelated to manuscript)
if profile:
    print('Profiling...')
    mk_reactor = InfiniteDilutionReactor()
    mk_sys = System(reactor=mk_reactor, reactions=mk_reactions)
    start_state = dict()
    for gind, g in enumerate(gas_states):
        if g in mk_sys.snames:
            if gind != 5:
                start_state[g] = gas_fracts[gind] * 0.02 * p / bartoPa
            else:
                start_state[g] = gas_fracts[gind] * p / bartoPa
    start_state['*'] = 1.0
    print('- without Jacobian:')
    mk_sys.set_parameters(times=times, start_state=start_state, T=T, p=p,
                          use_jacobian=False, verbose=False,
                          atol=atol, rtol=rtol, xtol=xtol)
    # draw_call_graph(mk_sys.solve_odes, path=figures_dir + 'pycallgraph_nojac.png')
    # run_timed(mk_sys.solve_odes)
    run_cprofiler('mk_sys.solve_odes()')
    print('Number of steps: %1.0f' % len(mk_sys.times))

    print('\n* * *\n')
    print('- with Jacobian:')
    mk_sys.set_parameters(times=times, start_state=start_state, T=T, p=p,
                          use_jacobian=True, verbose=False,
                          atol=atol, rtol=rtol, xtol=xtol)
    # draw_call_graph(mk_sys.solve_odes, path=figures_dir + 'pycallgraph.png')
    # run_timed(mk_sys.solve_odes)
    run_cprofiler('mk_sys.solve_odes()')
    print('Number of steps: %1.0f' % len(mk_sys.times))
    print('Done.')
