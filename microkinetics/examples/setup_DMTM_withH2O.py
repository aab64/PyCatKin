from microkinetics.classes.state import State
from microkinetics.classes.reaction import Reaction
from microkinetics.classes.system import System
from microkinetics.classes.reactor import *
from microkinetics.classes.energy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
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


font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['lines.linewidth'] = 1.5

# Conditions
p = 1.0e5  # Pressure (Pa)
T = 448.15  # Temperature (K)
Apore = 3.8e-10 ** 2  # Pore area (m2) - taken from wikipedia for SSZ-13
verbose = False  # Print messages
use_jacobian = False  # Use Jacobian to solve SS and ODEs
savexyz = False  # Save xyz files (not used)
savefig = True

# Location of outcars and frequencies
ads_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/energies/withH2O/'
ads_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/vibrations/withH2O/'
gas_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/energies/molecules/'
gas_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/vibrations/molecules/'

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/methanol/results/withH2O/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/methanol/images/'

print('-----------------------')
print('System 2: DMTM with H2O')
print('-----------------------')

# Load states
print('Configuring states...')

# Gases
states = [State(state_type='gas', path='o2', sigma=2)]
states += [State(state_type='gas', path='ch4', sigma=12)]
states += [State(state_type='gas', path='ch3oh', sigma=1)]
states += [State(state_type='gas', path='h2o', sigma=2)]
# Adsorbates
states += [State(state_type='adsorbate', path='2CuH2O')]
states += [State(state_type='adsorbate', path='Cu-pairH2O')]
states += [State(state_type='adsorbate', path='CuO2CuH2O')]
# states += [State(state_type='adsorbate', path='CuOO-CuH2O')]
states += [State(state_type='adsorbate', path='CuOOCuH2O')]
states += [State(state_type='adsorbate', path='s2Och4H2O')]
states += [State(state_type='adsorbate', path='sOsCH3OHH2O')]
states += [State(state_type='adsorbate', path='sOH2O')]
states += [State(state_type='adsorbate', path='sOch4H2O')]
states += [State(state_type='adsorbate', path='sOHsCH3H2O')]
states += [State(state_type='adsorbate', path='sCH3OHH2O')]
states += [State(state_type='adsorbate', path='sH2O')]
# Transition states
states += [State(state_type='TS', path='ts1H2O')]
states += [State(state_type='TS', path='ts2H2O')]
states += [State(state_type='TS', path='ts3H2O')]
states += [State(state_type='TS', path='ts4H2O')]
states += [State(state_type='TS', path='ts5H2O')]

snames = [s.name for s in states]
states = dict(zip(snames, states))

for s in snames:
    opt_dir = gas_opt_dir if states[s].state_type == 'gas' else ads_opt_dir
    vib_dir = gas_vib_dir if states[s].state_type == 'gas' else ads_vib_dir
    states[s].path = opt_dir + states[s].name
    states[s].vibs_path = vib_dir + states[s].name
    read_from_alternate = None
    if states[s].state_type == 'gas':
        read_from_alternate = {'get_atoms': lambda state_path=opt_dir + states[s].name: load_from_contcar(state_path)}
    states[s].read_from_alternate = read_from_alternate

frac = 0.67
for s in states.keys():
    if states[s].state_type != 'gas':
        states[s].gasdata = {'fraction': [frac], 'state': [states['h2o']]}
states['s2Och4H2O'].gasdata['fraction'].append(frac)
states['s2Och4H2O'].gasdata['state'].append(states['ch4'])
states['sOch4H2O'].gasdata['fraction'].append(frac)
states['sOch4H2O'].gasdata['state'].append(states['ch4'])

print('Done.')

# Reactions
print('Configuring reactions...')

reactions = [Reaction(reac_type='Arrhenius',
                      reactants=[states['2CuH2O']],
                      products=[states['Cu-pairH2O']],
                      TS=[states['ts1H2O']],
                      area=Apore,
                      name='r0')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['Cu-pairH2O'], states['o2']],
                       products=[states['CuO2CuH2O']],
                       TS=None,
                       area=Apore,
                       name='r1')]
# reactions += [Reaction(reac_type='Arrhenius',
#                        reactants=[states['CuO2CuH2O']],
#                        products=[states['CuOO-CuH2O']],
#                        TS=None,
#                        area=Apore,
#                        name='r2')]
# reactions += [Reaction(reac_type='Arrhenius',
#                        reactants=[states['CuOO-CuH2O']],
#                        products=[states['CuOOCuH2O']],
#                        TS=None,
#                        area=Apore,
#                        name='r3')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['CuO2CuH2O']],
                       products=[states['CuOOCuH2O']],
                       TS=None,
                       area=Apore,
                       name='r2')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CuOOCuH2O'], states['ch4']],
                       products=[states['s2Och4H2O']],
                       TS=None,
                       area=Apore,
                       name='r3')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['s2Och4H2O']],
                       products=[states['sOsCH3OHH2O']],
                       TS=[states['ts2H2O']],
                       area=Apore,
                       name='r4')]
reactions += [Reaction(reac_type='desorption',
                       reactants=[states['sOsCH3OHH2O']],
                       products=[states['sOH2O'], states['ch3oh']],
                       TS=None,
                       area=Apore,
                       name='r5')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['sOH2O'], states['ch4']],
                       products=[states['sOch4H2O']],
                       TS=None,
                       area=Apore,
                       name='r6')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sOch4H2O']],
                       products=[states['sOHsCH3H2O']],
                       TS=[states['ts3H2O']],
                       area=Apore,
                       name='r7')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sOHsCH3H2O']],
                       products=[states['sCH3OHH2O']],
                       TS=[states['ts4H2O']],
                       area=Apore,
                       name='r8')]
reactions += [Reaction(reac_type='desorption',
                       reactants=[states['sCH3OHH2O']],
                       products=[states['sH2O'], states['ch3oh']],
                       TS=None,
                       area=Apore,
                       name='r9')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sH2O']],
                       products=[states['Cu-pairH2O']],
                       TS=[states['ts5H2O']],
                       area=Apore,
                       name='r10')]

rnames = [r.name for r in reactions]
reactions = dict(zip(rnames, reactions))

print('Done.')

# System
print('Configuring system...')

reactor = InfiniteDilutionReactor()
sys = System(reactions, reactor)

print('Done.')

# Solve
print('Solving ODEs...')

# Initial conditions (2CH4 + O2 => 2CH3OH)
start_state = dict()
start_state['2CuH2O'] = 1.0
start_state['o2'] = 0.10
start_state['ch4'] = 0.02
start_state['ch3oh'] = 1.0e-11

Ts = np.linspace(start=400, stop=800, num=20, endpoint=True)
final_rates = np.zeros((len(Ts), len(reactions)))
final_cover = np.zeros((len(Ts), len([s for s in snames if states[s].state_type == 'adsorbate'])))
dorc = dict()
ipath = figures_dir + 'withH2O/' if savefig else None
for Tind, T in enumerate(Ts):

    times = np.logspace(start=int(np.log10(1e-16)), stop=int(np.log10(1e6)), num=int(1e6))
    sys.set_parameters(times=times, start_state=start_state, T=T, p=p, inflow_state=None,
                       use_jacobian=use_jacobian, verbose=verbose, xtol=1e-8)
    sys.solve_odes()
    sys.find_steady(store_steady=True)
    sys.reaction_terms(sys.full_steady)
    final_rates[Tind, :] = sys.rates[:, 0] - sys.rates[:, 1]
    final_cover[Tind, :] = sys.full_steady[sys.adsorbate_indices]
    tof = final_rates[Tind, 5] + final_rates[Tind, 9]

    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, tof, 1.0e-3 *
                                                               np.log((h * tof) / (kB * T)) * (-R * T)))

    tof_terms = ['r5', 'r9']
    dorc[T] = sys.degree_of_rate_control(tof_terms, ss_solve=False, eps=1.0e-3)

    print('Done.')


rates_output = pd.DataFrame(data=final_rates, index=Ts, columns=reactions.keys())
cover_output = pd.DataFrame(data=final_cover, index=Ts, columns=[s for s in sys.snames if states[s].state_type ==
                                                                 'adsorbate'])
dorc_output = pd.DataFrame(data=dorc)
rates_output.to_csv(results_dir + 'rates.csv')
cover_output.to_csv(results_dir + 'cover.csv')
dorc_output.to_csv(results_dir + 'dorc.csv')

fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot(Ts, final_rates[:, 6], 'o-', color='teal', label='first')
ax.plot(Ts, final_rates[:, 10], '+-', color='tomato', label='second')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)', title=r'With H$_2$O')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withH2O/tof_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(6.6, 3.3))
for i in range(final_cover.shape[1]):
    ax.plot(Ts, final_cover[:, i], label=sys.snames[sys.adsorbate_indices[i]])
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Coverage')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withH2O/cover_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(6.6, 3.3))
for i in range(final_rates.shape[1]):
    ax.plot(Ts, final_rates[:, i], label=list(reactions.keys())[i])
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withH2O/rates_temperature.png', dpi=300)

cmap = plt.get_cmap("Accent", len(rnames))
fig, ax = plt.subplots(figsize=(3.3, 3.3))
for rind, r in enumerate(rnames):
    if np.max([dorc[i][r] for i in Ts]) > 1e-2:
        ax.plot(Ts, [dorc[i][r] for i in Ts], '-', label=r, color=cmap.colors[rind, :])
ax.set(xlabel='Temperature (K)', ylabel=r'$\chi$ CH$_3$OH formation', title=r'With H$_2$O')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withH2O/dorc_temperature.png', dpi=300)

cmap = plt.get_cmap("Accent", len(sys.adsorbate_indices))

fig, ax = plt.subplots()
for r in ['r3', 'r4', 'r5', 'r6', 'r7', 'r9']:
    ax.plot(Ts, [reactions[r].get_reaction_energy(T=T, p=p) / eVtokJ / 1000 for T in Ts], label=r)
    if r in ['r4', 'r7']:
        ax.plot(Ts, [reactions[r].get_reaction_barriers(T=T, p=p)[0] / eVtokJ / 1000 for T in Ts], label=r+'TS')
ax.legend()
ax.grid()
ax.set(xlabel='Temperature (K)', ylabel='dG (eV)', title='With H2O')
fig.tight_layout()

minima = dict()
minima[0] = [states['2CuH2O'], states['o2'], states['ch4'], states['ch4']]
minima[1] = [states['ts1H2O'], states['o2'], states['ch4'], states['ch4']]
minima[2] = [states['Cu-pairH2O'], states['o2'], states['ch4'], states['ch4']]
minima[3] = [states['CuO2CuH2O'], states['ch4'], states['ch4']]
# minima[4] = [states['CuOO-CuH2O'], states['ch4'], states['ch4']]
minima[4] = [states['CuOOCuH2O'], states['ch4'], states['ch4']]
minima[5] = [states['s2Och4H2O'], states['ch4']]
minima[6] = [states['ts2H2O'], states['ch4']]
minima[7] = [states['sOsCH3OHH2O'], states['ch4']]
minima[8] = [states['sOH2O'], states['ch4'], states['ch3oh']]
minima[9] = [states['sOch4H2O'], states['ch3oh']]
minima[10] = [states['ts3H2O'], states['ch3oh']]
minima[11] = [states['sOHsCH3H2O'], states['ch3oh']]
minima[12] = [states['ts4H2O'], states['ch3oh']]
minima[13] = [states['sCH3OHH2O'], states['ch3oh']]
minima[14] = [states['sH2O'], states['ch3oh'], states['ch3oh']]
minima[15] = [states['ts5H2O'], states['ch3oh'], states['ch3oh']]
minima[16] = [states['Cu-pairH2O'], states['ch3oh'], states['ch3oh']]
minima[17] = [states['ts1H2O'], states['ch3oh'], states['ch3oh']]
minima[18] = [states['2CuH2O'], states['ch3oh'], states['ch3oh']]

energy = Energy(minima=minima)
energy.draw_energy_landscape(T=450, p=p, verbose=False, etype='electronic')
energy.draw_energy_landscape(T=450, p=p, verbose=False, etype='free')

fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot(Ts, final_rates[:, 5] + final_rates[:, 9], 'o-', color='tomato', label='MK')
for T in Ts:
    tof, num_i, num_j, lTi, lIj = energy.evaluate_energy_span_model(T=T, p=p, verbose=False, etype='free')
    ax.plot(T, 2.0 * tof, '*-', color='teal', label='')
    print('ES:')
    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, tof, 1.0e-3 *
                                                               np.log((h * tof * 2.0) / (kB * T)) * (-R * T)))
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)', title=r'With H$_2$O')
ax.legend(('MK', 'ES'), frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withH2O/tof_cmp.png', dpi=300)
