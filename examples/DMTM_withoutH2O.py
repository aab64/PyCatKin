from pycatkin.classes.state import State
from pycatkin.classes.reaction import Reaction
from pycatkin.classes.system import System
from pycatkin.classes.reactor import *
from pycatkin.classes.energy import *
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
ads_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/energies/withoutH2O/'
ads_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/vibrations/withoutH2O/'
gas_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/energies/molecules/'
gas_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/vibrations/molecules/'

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Cu/results/withoutH2O/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Cu/images/'

print('--------------------------')
print('System 1: DMTM without H2O')
print('--------------------------')

# Load states
print('Configuring states...')

# Adsorbates
states = [State(state_type='adsorbate', name='2cu')]
states += [State(state_type='adsorbate', name='Cu-pair')]
states += [State(state_type='adsorbate', name='CuO2Cu')]
# states += [State(state_type='adsorbate', name='CuOO-Cu')]
states += [State(state_type='adsorbate', name='CuOOCu')]
states += [State(state_type='adsorbate', name='s2Och4')]
states += [State(state_type='adsorbate', name='sOsCH3OH')]
states += [State(state_type='adsorbate', name='sO')]
states += [State(state_type='adsorbate', name='sOch4')]
states += [State(state_type='adsorbate', name='sOHsCH3')]
states += [State(state_type='adsorbate', name='sCH3OH')]
states += [State(state_type='adsorbate', name='s')]
# Transition states
states += [State(state_type='TS', name='ts1')]
states += [State(state_type='TS', name='ts2')]
states += [State(state_type='TS', name='ts3')]
states += [State(state_type='TS', name='ts4')]
states += [State(state_type='TS', name='ts5')]
states += [State(state_type='TS', name='ts6')]
# Gases
states += [State(state_type='gas', name='o2', sigma=2)]
states += [State(state_type='gas', name='ch4', sigma=12)]
states += [State(state_type='gas', name='ch3oh', sigma=1)]
states += [State(state_type='gas', name='h2o', sigma=2)]

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

# Add gas phase entropy for hindered molecules
frac = 0.67
states['s2Och4'].gasdata = {'fraction': [frac], 'state': [states['ch4']]}
states['sOch4'].gasdata = {'fraction': [frac], 'state': [states['ch4']]}

print('Done.')

# Reactions
print('Configuring reactions...')

reactions = [Reaction(reac_type='Arrhenius',
                      reactants=[states['2cu']],
                      products=[states['Cu-pair']],
                      TS=[states['ts1']],
                      area=Apore,
                      name='r0')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['Cu-pair'], states['o2']],
                       products=[states['CuO2Cu']],
                       TS=None,
                       area=Apore,
                       name='r1')]
# reactions += [Reaction(reac_type='Arrhenius',
#                        reactants=[states['CuO2Cu']],
#                        products=[states['CuOO-Cu']],
#                        TS=[states['ts2']],
#                        area=Apore,
#                        name='r2')]
# reactions += [Reaction(reac_type='Arrhenius',
#                        reactants=[states['CuOO-Cu']],
#                        products=[states['CuOOCu']],
#                        TS=None,
#                        area=Apore,
#                        name='r3')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['CuO2Cu']],
                       products=[states['CuOOCu']],
                       TS=[states['ts2']],
                       area=Apore,
                       name='r2')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CuOOCu'], states['ch4']],
                       products=[states['s2Och4']],
                       TS=None,
                       area=Apore,
                       name='r3')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['s2Och4']],
                       products=[states['sOsCH3OH']],
                       TS=[states['ts3']],
                       area=Apore,
                       name='r4')]
reactions += [Reaction(reac_type='desorption',
                       reactants=[states['sOsCH3OH']],
                       products=[states['sO'], states['ch3oh']],
                       TS=None,
                       area=Apore,
                       name='r5')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['sO'], states['ch4']],
                       products=[states['sOch4']],
                       TS=None,
                       area=Apore,
                       name='r6')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sOch4']],
                       products=[states['sOHsCH3']],
                       TS=[states['ts4']],
                       area=Apore,
                       name='r7')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sOHsCH3']],
                       products=[states['sCH3OH']],
                       TS=[states['ts5']],
                       area=Apore,
                       name='r8')]
reactions += [Reaction(reac_type='desorption',
                       reactants=[states['sCH3OH']],
                       products=[states['s'], states['ch3oh']],
                       TS=None,
                       area=Apore,
                       name='r9')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['s']],
                       products=[states['Cu-pair']],
                       TS=[states['ts6']],
                       area=Apore,
                       name='r10')]

rnames = [r.name for r in reactions]
reactions = dict(zip(rnames, reactions))

print('Done.')

# System
print('Configuring system...')

reactor = InfiniteDilutionReactor()
sys = System(reactions=reactions, reactor=reactor)

print('Done.')

# Solve
print('Solving ODEs...')

# Initial conditions
start_state = dict()
start_state['2cu'] = 1.0
start_state['o2'] = 0.10
start_state['ch4'] = 0.02
start_state['ch3oh'] = 1.0e-11

Ts = np.linspace(start=400, stop=800, num=20, endpoint=True)
dorc = dict()
final_rates = np.zeros((len(Ts), len(reactions)))
final_cover = np.zeros((len(Ts), len([s for s in snames if states[s].state_type == 'adsorbate'])))
ipath = figures_dir + 'withoutH2O/' if savefig else None
for Tind, T in enumerate(Ts):

    times = np.logspace(start=int(np.log10(1e-12)), stop=int(np.log10(1e12)), num=int(1e6))
    sys.set_parameters(times=times, start_state=start_state, inflow_state=None, T=T, p=p,
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

rates_output = pd.DataFrame(data=final_rates, index=Ts, columns=list(reactions.keys()))
cover_output = pd.DataFrame(data=final_cover, index=Ts, columns=[s for s in sys.snames
                                                                 if states[s].state_type == 'adsorbate'])
dorc_output = pd.DataFrame(data=dorc)
rates_output.to_csv(results_dir + 'rates.csv')
cover_output.to_csv(results_dir + 'cover.csv')
dorc_output.to_csv(results_dir + 'dorc.csv')

fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot(Ts, final_rates[:, 6], 'o-', color='teal', label='first')
ax.plot(Ts, final_rates[:, 10], '+-', color='tomato', label='second')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)', title=r'Without H$_2$O')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withoutH2O/tof_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(6.6, 3.3))
for i in range(final_cover.shape[1]):
    ax.plot(Ts, final_cover[:, i], label=sys.snames[sys.adsorbate_indices[i]])
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Coverage')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withoutH2O/cover_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(6.6, 3.3))
for i in range(final_rates.shape[1]):
    ax.plot(Ts, final_rates[:, i], label=list(reactions.keys())[i])
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withoutH2O/rates_temperature.png', dpi=300)

cmap = plt.get_cmap("Accent", len(rnames))
fig, ax = plt.subplots(figsize=(3.3, 3.3))
for rind, r in enumerate(rnames):
    if np.max([dorc[i][r] for i in Ts]) > 1e-2:
        ax.plot(Ts, [dorc[i][r] for i in Ts], '-', label=r, color=cmap.colors[rind, :])
ax.set(xlabel='Temperature (K)', ylabel=r'$\chi$ CH$_3$OH formation', title=r'Without H$_2$O')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withoutH2O/dorc_temperature.png', dpi=300)

cmap = plt.get_cmap("Accent", len(sys.adsorbate_indices))

fig, ax = plt.subplots()
for r in ['r3', 'r4', 'r5', 'r6', 'r7', 'r9']:
    ax.plot(Ts, [reactions[r].get_reaction_energy(T=T, p=p) / eVtokJ / 1000 for T in Ts], label=r)
    if r in ['r4', 'r7']:
        ax.plot(Ts, [reactions[r].get_reaction_barriers(T=T, p=p)[0] / eVtokJ / 1000 for T in Ts], label=r+'TS')
ax.legend()
ax.grid()
ax.set(xlabel='Temperature (K)', ylabel='dG (eV)', title='Without H2O')
fig.tight_layout()

clrs = ['tomato', 'dodgerblue', 'slategrey', 'teal']
fig, ax = plt.subplots(figsize=(3.33, 3.33))
for rind, r in enumerate(['r6', 'r7', 'r8', 'r9']):
    ax.plot(Ts, [reactions[r].get_reaction_energy(T=T, p=p) / eVtokJ / 1000 for T in Ts], label=r,
            color=clrs[rind])
    if r in ['r7', 'r8']:
        ax.plot(Ts, [reactions[r].get_reaction_barriers(T=T, p=p)[0] / eVtokJ / 1000 for T in Ts], ':', label='',
                color=clrs[rind])
ax.legend(frameon=False)
ax.set(xlabel='Temperature (K)', ylabel='dG (eV)', title='Without H2O')
fig.tight_layout()
plt.savefig(figures_dir + 'withoutH2O/dGs.png', format='png', dpi=300)

minima = dict()
minima[0] = [states['2cu'], states['o2'], states['ch4'], states['ch4']]
minima[1] = [states['ts1'], states['o2'], states['ch4'], states['ch4']]
minima[2] = [states['Cu-pair'], states['o2'], states['ch4'], states['ch4']]
minima[3] = [states['CuO2Cu'], states['ch4'], states['ch4']]
minima[4] = [states['ts2'], states['ch4'], states['ch4']]
# minima[5] = [states['CuOO-Cu'], states['ch4'], states['ch4']]
minima[5] = [states['CuOOCu'], states['ch4'], states['ch4']]
minima[6] = [states['s2Och4'], states['ch4']]
minima[7] = [states['ts3'], states['ch4']]
minima[8] = [states['sOsCH3OH'], states['ch4']]
minima[9] = [states['sO'], states['ch4'], states['ch3oh']]
minima[10] = [states['sOch4'], states['ch3oh']]
minima[11] = [states['ts4'], states['ch3oh']]
minima[12] = [states['sOHsCH3'], states['ch3oh']]
minima[13] = [states['ts5'], states['ch3oh']]
minima[14] = [states['sCH3OH'], states['ch3oh']]
minima[15] = [states['s'], states['ch3oh'], states['ch3oh']]
minima[16] = [states['ts6'], states['ch3oh'], states['ch3oh']]
minima[17] = [states['Cu-pair'], states['ch3oh'], states['ch3oh']]
minima[18] = [states['ts1'], states['ch3oh'], states['ch3oh']]
minima[19] = [states['2cu'], states['ch3oh'], states['ch3oh']]

peslabs = [i[0].name for i in minima.values()]

energy = Energy(minima=minima, labels=peslabs)
energy.construct_energy_landscape(T=450, p=p, verbose=False)
energy.draw_energy_landscape(T=450, p=p, verbose=False, etype='electronic')

fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot(Ts, final_rates[:, 5] + final_rates[:, 9], 'o-', color='tomato', label='MK')
for T in Ts:
    tof, Espan, TDTS, TDI, num_i, num_j, lTi, lIj = energy.evaluate_energy_span_model(T=T, p=p, etype='free',
                                                                                      verbose=False)
    ax.plot(T, 2.0 * tof, '*-', color='teal', label='')

ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)', title=r'Without H$_2$O')
ax.legend(('MK', 'ES'), frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'withoutH2O/tof_cmp.png', dpi=300)
