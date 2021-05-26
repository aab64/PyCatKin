from microkinetics.classes.state import *
from microkinetics.classes.energy import *
from microkinetics.classes.reaction import *
from microkinetics.classes.system import *
from microkinetics.classes.reactor import *
import glob
import copy
from ase.data.pubchem import pubchem_atoms_search
from ase.io import read


def load_acetaldehyde_from_pubchem(outcar_path):
    aa = pubchem_atoms_search(name='acetaldehyde')
    inertia = copy.copy(aa.get_moments_of_inertia())
    vasp_atoms = ase.io.read(outcar_path, format='vasp-out')
    mass = sum(vasp_atoms.get_masses())
    return vasp_atoms, mass, inertia


def load_ethylacetate_frequencies(freq_path):
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
        else:
            break
    return freq, i_freq


def load_ethyl_acetate_energy(energy_path):
    energy_path = energy_path.split('\\')[0] + '/energies/' + energy_path.split('\\')[1] + '.dat'
    if os.path.isfile(energy_path):
        with open(energy_path, 'r') as fd:
            lines = fd.readlines()
            Gelec = float(lines[0].split()[0])
    return Gelec


# Conditions
p = 1.01325e5  # Pressure (Pa)
T = 723  # Temperature (K)
Ts = list(np.linspace(start=523, stop=923,
                      num=17, endpoint=True))  # Temperatures (K)
Asite = (13e-10 * 15e-10) / 4  # Site area (m2)
times = np.logspace(start=-14, stop=3, num=int(1e5)) * 24 * 3.6  # Times (s)
use_jacobian = True  # Use Jacobian to solve SS and ODEs
verbose = True  # Print messages

# Location of outcars and frequencies
adsdir = 'D:/Users/Astrid/Documents/Chalmers/Data/Butadiene/DFT data/'
folders = ['dehydrogenation', 'aldol condensation',
           'MPV red', 'aldol mpv',
           'dehydration', 'prins', 'etoxetoh',
           'adsorption', 'H2', 'gas',
           'butanol_test', 'ethyl acetate']

# Location of results files and images
# results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/butadiene/results/'
figures_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Butadiene/images/'

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
                ea_path = path
                read_from_alternate = dict()
                read_from_alternate['get_electronic_energy'] = lambda: load_ethyl_acetate_energy(ea_path)
                read_from_alternate['get_vibrations'] = lambda: load_ethylacetate_frequencies(ea_path)
            elif 'acetaldehyde' in path:
                aa_path = path
                read_from_alternate = dict()
                read_from_alternate['get_atoms'] = lambda: load_acetaldehyde_from_pubchem(aa_path)
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

# for s in snames:
#     states[s].view_atoms(path=figures_dir)
#     states[s].save_pdb(path=figures_dir)
#     states[s].calc_electronic_energy(verbose=verbose)

# for Ts in [623, 723, 823]:
#     print('At %1.0f K' % Ts)
#     for k in states.keys():
#         print(k)
#         print(states[k].get_free_energy(T=Ts, p=p, verbose=verbose))
#         print('--')
#
# (states['7Eiii'].get_free_energy(T=T, p=p) - (states['7Eiv'].get_free_energy(T=T, p=p) +
#                                               states['ethylacetate'].get_free_energy(T=T, p=p))) * eVtokcal
#
# for k in states.keys():
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
# for s in snames:
#     print(s)
#     print('%1.1f' % (states[s].Gvibr * eVtokcal))

# Base reactions
reactions = dict()
initial_states = ['1A', '2A', '2F', '2J', '2L', '3A', '3D', '3F',
                  '4A', '4D', '4D', '4F', '4I', '5A',
                  '6A', '6C', '6E', '6G', '7A', '9D',
                  '8A', '13A', '14A', '7Ei']
final_states = ['1C', '2C', '2H', '2L', '2N', '3C', '3F', '3G',
                '4C', '4Ca', '4F', '4H', '4K', '5C',
                '6C', '6E', '6G', '6H', '7E', '9C',
                '8C', '13C', '14C', '7Eiii']
transition_states = ['1B', '2B', '2G', '2K', '2M', '3B', '3E', '3Fi',
                     '4B', None, '4E', '4G', '4J', '5B',
                     '6B', '6D', '6F', '6Gi', '7D', None,
                     '8B', '13B', '14B', '7Eii']
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

# MK model
# ========
gas_states = ['H2', 'H2O', 'C2H4', 'C4H6', 'C2H4O', 'C2H5OH', 'C4H6O', 'C4H10O', 'C4H8O2']
gas_sigmas = [2 if i == 'H2' or i == 'H2O' else 1 for i in gas_states]
gas_masses = [2.016, 18.02, 28.054, 54.09, 44.05, 46.07, 70.09, 74.12, 72.10, 88.11]
ads_states = ['H*', 'O*', 'OH*',
              'C2H3O*', 'C2H4O*', 'C2H5O*',
              'C4H6O_1*', 'C4H6O_2*', 'C4H6O_3*',
              'C4H7O_1*', 'C4H7O_2*', 'C4H7O_3*',
              'C4H8O_1*', 'C4H8O_2*', 'C4H8O_3*',
              'C4H9O*',
              'C4H6O2*', 'C4H7O2*',
              'C4H8O2_1*', 'C4H8O2_2*', 'C4H8O2_3*',
              'C4H9O2_1*', 'C4H9O2_2*']
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
                # ['C4H7O_1*', 'H*'],
                # ['C4H8O_2*'],
                # ['C4H8O_3*', 'H*'],
                # ['C4H10O', '*', '*'],
                # ['C4H7O_1*', 'H*', 'H*'],
                # ['C4H10O', '*', '*'],
                # ['C4H9O2_2*', '*'],
                # ['C4H8O2', '*']
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
               # ['C4H8O_2*', '*'],
               # ['C4H8O_3*'],
               # ['C4H9O*', '*'],
               # ['C4H9O*', 'H*'],
               # ['C4H9O*', '*', '*'],
               # ['C4H9O*', 'H*'],
               # ['C4H8O2_3*', 'H*'],
               # ['C4H8O2_3*']
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
             # '12A-12C',
             # '12C-12E',
             # '12E-12G',
             # '0-12G',
             # '13A-13C',
             # '14A-14C',
             # '7Ei-7Eiii',
             # 'ethylacetate-7Eiii'
             ]

mk_reactions = dict()
for r in range(len(mk_rnames)):
    rname = mk_rnames[r]
    rtype = reactions[rname].reac_type if (rname != '8A-8C' and rname != '14A-14C') else 'adsorption'
    mk_reactions[rname] = ReactionDerivedReaction(name=rname,
                                                  area=Asite,
                                                  reac_type=rtype,
                                                  base_reaction=reactions[rname],
                                                  reactants=[mk_states[i] for i in mk_reactants[r]],
                                                  products=[mk_states[i] for i in mk_products[r]])

# MK system
xgas = [0.01, 1.0e-8, 0.01, 1.0e-8, 1.0e-8, 0.02, 1.0e-8, 1.0e-8, 1.0e-8]
mk_reactor = InfiniteDilutionReactor()
mk_sys = System(reactor=mk_reactor, reactions=mk_reactions)
start_state = dict()
for gind, g in enumerate(gas_states):
    if gind < len(gas_states) - 2:
        if gind != 5:
            start_state[g] = xgas[gind] * 0.02 * p / bartoPa
        else:
            start_state[g] = xgas[gind] * p / bartoPa
start_state['*'] = 1.0


mk_tofs = []
for Ti in Ts:
    mk_sys.set_parameters(times=times, start_state=start_state, T=Ti, p=p,
                          use_jacobian=use_jacobian, verbose=verbose)
    mk_sys.solve_odes()
    mk_sys.find_steady(plot_comparison=True, store_steady=True)
    mk_sys.reaction_terms(y=mk_sys.full_steady)
    # mk_sys.plot_transient()
    mk_tofs.append(sum([mk_sys.rates[i, 0] - mk_sys.rates[i, 1]
                        for i in [list(mk_sys.reactions.keys()).index('3F-3G'),
                                  list(mk_sys.reactions.keys()).index('4I-4K'),
                                  list(mk_sys.reactions.keys()).index('6G-6H')]]))

# ES model
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

energy_p123 = Energy(minima=p123)
energy_p123.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol', path=figures_dir)

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

energy_p124 = Energy(minima=p124)
energy_p124.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol', path=figures_dir)

p156 = dict()
p156[0] = [states['surface'], states['surface'], states['ethanol'], states['ethanol'], states['ethanol']]
p156[1] = [states['1A'], states['ethanol'], states['ethanol'], states['surface']]
p156[2] = [states['1B'], states['ethanol'], states['ethanol'], states['surface']]
p156[3] = [states['1C'], states['ethanol'], states['ethanol'], states['surface']]
p156[4] = [states['1C'], states['5A'], states['ethanol']]
p156[5] = [states['1C'], states['5B'], states['ethanol']]
p156[6] = [states['1C'], states['5C'], states['ethanol']]
p156[7] = [states['6A'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[8] = [states['6B'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[9] = [states['6C'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[10] = [states['6D'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[11] = [states['6E'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[12] = [states['6F'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[13] = [states['6G'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[14] = [states['6Gi'], states['ethanol'], states['H2'], states['H2O'], states['surface']]
p156[15] = [states['6H'], states['ethanol'], states['H2'], states['H2O'], states['surface']]

energy_p156 = Energy(minima=p156)
energy_p156.draw_energy_landscape(T=T, p=p, etype='free', eunits='kcal/mol', path=figures_dir)

tofs = []
fig, ax = plt.subplots(figsize=(3.2, 3.2))
for Ti in Ts:
    tof, num_i, num_j, lTi, lIj = energy_p123.evaluate_energy_span_model(T=Ti, p=p, verbose=False, etype='free')
    tofs.append(tof)
    print('ES:')
    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (Ti, tof, 1.0e-3 *
                                                               np.log((h * tof * 2.0) / (kB * Ti)) * (-R * Ti)))
ax.plot(Ts, tofs, color='darkorchid', label='p123')

tofs = []
for Ti in Ts:
    tof, num_i, num_j, lTi, lIj = energy_p124.evaluate_energy_span_model(T=Ti, p=p, verbose=False, etype='free')
    tofs.append(tof)
    print('ES:')
    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (Ti, tof, 1.0e-3 *
                                                               np.log((h * tof * 2.0) / (kB * Ti)) * (-R * Ti)))
ax.plot(Ts, tofs, color='dodgerblue', label='p124')

tofs = []
for Ti in Ts:
    tof, num_i, num_j, lTi, lIj = energy_p156.evaluate_energy_span_model(T=Ti, p=p, verbose=False, etype='free')
    tofs.append(tof)
    print('ES:')
    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (Ti, tof, 1.0e-3 *
                                                               np.log((h * tof * 2.0) / (kB * Ti)) * (-R * Ti)))
ax.plot(Ts, tofs, color='darkorange', label='p156')
ax.plot(Ts, mk_tofs, color='black', label='MK')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='TOF (1/s)', ylim=(1e-8, 1e2), xlim=(523, 923))
ax.grid()
fig.tight_layout()
fig.show()


