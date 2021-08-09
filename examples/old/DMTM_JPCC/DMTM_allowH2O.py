from pycatkin.classes.state import State
from pycatkin.classes.reaction import Reaction
from pycatkin.classes.system import System
from pycatkin.classes.reactor import *
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
savefig = True

# Load states
print('Configuring states...')

# Location of outcars and frequencies
gas_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/energies/molecules/'
gas_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/vibrations/molecules/'

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Cu/results/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Cu/images/allowH2O/'

print('------------------------')
print('System 3: DMTM with both')
print('------------------------')

print('- System 1: without H2O')

# Adsorbates
states = [State(state_type='adsorbate', name='2cu')]
states += [State(state_type='adsorbate', name='Cu-pair')]
states += [State(state_type='adsorbate', name='CuO2Cu')]
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

print('- System 2: with H2O')

# Adsorbates
states += [State(state_type='adsorbate', path='2CuH2O')]
states += [State(state_type='adsorbate', path='Cu-pairH2O')]
states += [State(state_type='adsorbate', path='CuO2CuH2O')]
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

# Add path info and some gas entropy contributions
frac = 0.67
for s in states.keys():
    if 'H2O' in states[s].name:
        ads_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/energies/withH2O/'
        ads_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/vibrations/withH2O/'
        if states[s].state_type != 'gas':
            states[s].gasdata = {'fraction': [frac], 'state': [states['h2o']]}
    else:
        ads_opt_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/energies/withoutH2O/'
        ads_vib_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Cu/vibrations/withoutH2O/'
    opt_dir = gas_opt_dir if states[s].state_type == 'gas' else ads_opt_dir
    vib_dir = gas_vib_dir if states[s].state_type == 'gas' else ads_vib_dir
    states[s].path = opt_dir + states[s].name
    states[s].vibs_path = vib_dir + states[s].name
    read_from_alternate = None
    if states[s].state_type == 'gas':
        read_from_alternate = {'get_atoms': lambda state_path=opt_dir + states[s].name: load_from_contcar(state_path)}
    states[s].read_from_alternate = read_from_alternate
states['s2Och4'].gasdata = {'fraction': [frac], 'state': [states['ch4']]}
states['sOch4'].gasdata = {'fraction': [frac], 'state': [states['ch4']]}
states['s2Och4H2O'].gasdata['fraction'].append(frac)
states['s2Och4H2O'].gasdata['state'].append(states['ch4'])
states['sOch4H2O'].gasdata['fraction'].append(frac)
states['sOch4H2O'].gasdata['state'].append(states['ch4'])

print('Done.')

# Reactions
print('Configuring reactions...')

# Without H2O
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

# With H2O
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['2cu'], states['h2o']],
                       products=[states['2CuH2O']],
                       TS=None,
                       area=Apore,
                       name='rXH2O')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['2CuH2O']],
                       products=[states['Cu-pairH2O']],
                       TS=[states['ts1H2O']],
                       area=Apore,
                       name='r0H2O')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['Cu-pairH2O'], states['o2']],
                       products=[states['CuO2CuH2O']],
                       TS=None,
                       area=Apore,
                       name='r1H2O')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['CuO2CuH2O']],
                       products=[states['CuOOCuH2O']],
                       TS=None,
                       area=Apore,
                       name='r2H2O')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CuOOCuH2O'], states['ch4']],
                       products=[states['s2Och4H2O']],
                       TS=None,
                       area=Apore,
                       name='r3H2O')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['s2Och4H2O']],
                       products=[states['sOsCH3OHH2O']],
                       TS=[states['ts2H2O']],
                       area=Apore,
                       name='r4H2O')]
reactions += [Reaction(reac_type='desorption',
                       reactants=[states['sOsCH3OHH2O']],
                       products=[states['sOH2O'], states['ch3oh']],
                       TS=None,
                       area=Apore,
                       name='r5H2O')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['sOH2O'], states['ch4']],
                       products=[states['sOch4H2O']],
                       TS=None,
                       area=Apore,
                       name='r6H2O')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sOch4H2O']],
                       products=[states['sOHsCH3H2O']],
                       TS=[states['ts3H2O']],
                       area=Apore,
                       name='r7H2O')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sOHsCH3H2O']],
                       products=[states['sCH3OHH2O']],
                       TS=[states['ts4H2O']],
                       area=Apore,
                       name='r8H2O')]
reactions += [Reaction(reac_type='desorption',
                       reactants=[states['sCH3OHH2O']],
                       products=[states['sH2O'], states['ch3oh']],
                       TS=None,
                       area=Apore,
                       name='r9H2O')]
reactions += [Reaction(reac_type='Arrhenius',
                       reactants=[states['sH2O']],
                       products=[states['Cu-pairH2O']],
                       TS=[states['ts5H2O']],
                       area=Apore,
                       name='r10H2O')]

rnames = [r.name for r in reactions]
reactions = dict(zip(rnames, reactions))

print('Done.')

# System
print('Configuring system...')

reactor = InfiniteDilutionReactor()
sys = System()

for s in states.keys():
    sys.add_state(state=states[s])
for r in reactions.keys():
    sys.add_reaction(reaction=reactions[r])
sys.add_reactor(reactor=reactor)
sys.names_to_indices()

print('Done.')

# Solve
print('Solving ODEs...')

# Initial conditions
start_state = dict()
start_state['2cu'] = 1.0
start_state['o2'] = 0.10
start_state['ch4'] = 0.02
start_state['ch3oh'] = 1.0e-11
start_state['h2o'] = 1.0e-11

Xs = np.logspace(start=-13, stop=-1, num=20, endpoint=True)
Ts = np.linspace(start=400, stop=800, num=20, endpoint=True)

dorc = dict()
rmap = np.zeros((len(Ts), len(Xs)))
final_rates = np.zeros((len(Ts), len(reactions)))
final_cover = np.zeros((len(Ts), len([s for s in snames if states[s].state_type == 'adsorbate'])))
ipath = figures_dir + 'allowH2O/' if savefig else None
for Tind, T in enumerate(Ts):
    for Xind, start_state['h2o'] in enumerate(Xs):

        times = np.logspace(start=int(np.log10(1e-12)), stop=int(np.log10(1e6)), num=int(1e4))
        sys.set_parameters(times=times, start_state=start_state, inflow_state=None, T=T, p=p,
                           use_jacobian=use_jacobian, verbose=verbose, xtol=1e-8)
        sys.solve_odes()
        # sys.find_steady(store_steady=True)
        # sys.reaction_terms(sys.full_steady)

        final_rates[Tind, :] = sys.rates[:, 0] - sys.rates[:, 1]
        # final_cover[Tind, :] = sys.full_steady[sys.adsorbate_indices]
        final_cover[Tind, :] = sys.solution[-1][sys.adsorbate_indices]
        rmap[Tind, Xind] = final_rates[Tind, 5] + final_rates[Tind, 9] + final_rates[Tind, 17] + final_rates[Tind, 21]

        tof = rmap[Tind, Xind]
        print('T = %1.0f K, xH2O = %1.2e, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, Xs[Xind], tof,
                                                                                 1.0e-3 * np.log((h * tof) / (kB * T)) *
                                                                                 (-R * T)))
        if Xind == len(Xs) - 1:
            tof_terms = ['r5', 'r9', 'r17', 'r21']
            dorc[T] = sys.degree_of_rate_control(tof_terms, eps=1.0e-3, ss_solve=False)

        print('Xind %1.0f/%1.0f done.' % (Xind + 1, len(Xs)))

    print('Tind %1.0f/%1.0f done.' % (Tind + 1, len(Ts)))

rmap_output = pd.DataFrame(data=rmap, index=Ts, columns=Xs)
rmap_output.to_csv(results_dir + 'rmap.csv')

fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot(Ts, final_rates[:, 5], 'o-', color='teal', label='1st')
ax.plot(Ts, final_rates[:, 9], '+-', color='tomato', label='2nd')
ax.plot(Ts, final_rates[:, 17], 's-', color='darkslategrey', label='1st [H2O]')
ax.plot(Ts, final_rates[:, 21], 'x-', color='royalblue', label='2nd [H2O]')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'tof_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot(Ts, final_rates[:, 5] + final_rates[:, 9] + final_rates[:, 17] + final_rates[:, 21], 'o-', color='k')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Rate (1/s)')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'tof_temperature.png', dpi=300)

cmap = plt.get_cmap("Accent", len(rnames))
fig, ax = plt.subplots(figsize=(3.3, 3.3))
for rind, r in enumerate(rnames):
    if np.max([dorc[i][r] for i in Ts]) > 1e-1:
        ax.plot(Ts, [dorc[i][r] for i in Ts], '-', label=r, color=cmap.colors[rind, :])
ax.set(xlabel='Temperature (K)', ylabel=r'$\chi$ CH$_3$OH formation')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'dorc_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(6.6, 3.3))
for i in range(final_cover.shape[1]):
    ax.plot(Ts, final_cover[:, i], label=sys.snames[sys.adsorbate_indices[i]])
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='Coverage')
ax.legend(frameon=False, loc='best')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'cover_temperature.png', dpi=300)

fig, ax = plt.subplots(figsize=(3.3, 3.3))
XX, TT = np.meshgrid(Xs, Ts)
CS = ax.contourf(XX, TT, np.log10(rmap), cmap=plt.cm.coolwarm)
ax.set(xscale='log', ylabel='Temperature (K)', xlabel='H$_2$O pressure (bar)', xlim=(1e-13, 1e-1))
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel(r'${\log}_{10}(\mathsf{Rate})$')
fig.tight_layout()
if savefig:
    plt.savefig(figures_dir + 'rmap.svg', dpi=300, format='svg')

elen_path = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/pycatkin/examples/DMTM/wetdata/energy/'
vibs_path = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/pycatkin/examples/DMTM/wetdata/vibrations/'
for s in sys.snames:
    if s is 'h2o':
        sys.states[s].save_energy(path=elen_path)
        sys.states[s].save_vibrations(vibs_path=vibs_path)
        if sys.states[s].state_type == 'gas':
            print(s)
            print(sys.states[s].mass)
            print(sys.states[s].shape)
            print(sys.states[s].inertia)
