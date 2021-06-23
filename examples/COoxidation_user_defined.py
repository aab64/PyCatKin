from pycatkin.classes.state import *
from pycatkin.classes.reaction import *
from pycatkin.classes.system import *
from pycatkin.classes.reactor import *
from pycatkin.classes.scaling import *
import numpy as np

# Conditions
p = 1.0e5  # Pressure (Pa)
Ts = list(np.linspace(start=423, stop=623, num=20, endpoint=True))  # Temperature (K)
use_jacobian = True  # Use Jacobian in SS and ODE solvers
verbose = False  # Print messages

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/results/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/images/'

print('--------------------')
print('System: CO oxidation')
print('--------------------')

# Load states
print('Configuring states...')

# Adsorbates
states = [State(state_type='adsorbate', name='CO.Pd111')]
states += [State(state_type='adsorbate', name='O.Pd111')]
states += [State(state_type='adsorbate', name='COonAu')]
states += [State(state_type='adsorbate', name='OonAu')]
states += [State(state_type='adsorbate', name='COonPd')]
states += [State(state_type='adsorbate', name='OonPd')]
states += [State(state_type='adsorbate', name='COonAu_alloy')]
states += [State(state_type='adsorbate', name='COonPd_alloy')]
states += [State(state_type='adsorbate', name='OonAu_alloy')]
states += [State(state_type='adsorbate', name='OonPd_alloy')]

# Surfaces
states += [State(state_type='surface', name='Clean')]

# Gases
states += [State(state_type='gas', name='O2', sigma=2, mass=32)]
states += [State(state_type='gas', name='CO', sigma=1, mass=28)]
states += [State(state_type='gas', name='CO2', sigma=2, mass=44)]

# Zip states into a dictionary
snames = [s.name for s in states]
states = dict(zip(snames, states))

print('Done.')

# Reactions
print('Configuring reactions...')

# Electronic energies from DFT
dEads_CO_Pd = -1.18
dEads_O_Pd = -1.14
dEads_CO_alloy = -1.34
dEads_O_alloy = -0.49
dEads_CO_AuTop = -0.093
dEads_O_AuTop = 0.35
dEads_CO_PdTop = -1.14
dEads_O_PdTop = -1.15

# Standard entropies (Atkins, J/molK)
SCOg = 197.67 * 1.0e-3 / eVtokJ
SO2g = 205.138 * 1.0e-3 / eVtokJ

# Free energies
dGads_CO_Pd = dict(zip(([0] + Ts), [dEads_CO_Pd + SCOg * T for T in ([0] + Ts)]))
dGads_O2_Pd = dict(zip(([0] + Ts), [dEads_O_Pd * 2.0 + SO2g * T for T in ([0] + Ts)]))
dGads_CO_alloy = dict(zip(([0] + Ts), [dEads_CO_alloy + SCOg * T for T in ([0] + Ts)]))
dGads_O2_alloy = dict(zip(([0] + Ts), [dEads_O_alloy * 2.0 + SO2g * T for T in ([0] + Ts)]))
dGads_CO_AuTop = dict(zip(([0] + Ts), [dEads_CO_AuTop + SCOg * T for T in ([0] + Ts)]))
dGads_O2_AuTop = dict(zip(([0] + Ts), [dEads_O_AuTop * 2.0 + SO2g * T for T in ([0] + Ts)]))
dGads_CO_PdTop = dict(zip(([0] + Ts), [dEads_CO_PdTop + SCOg * T for T in ([0] + Ts)]))
dGads_O2_PdTop = dict(zip(([0] + Ts), [dEads_O_PdTop * 2.0 + SO2g * T for T in ([0] + Ts)]))

# Site properties

# Pd
rPd = 2.0e-10  # Pd-Pd radius (m)
aPd = np.pi * rPd ** 2  # Pd site area (m2)
Acat_Pd = 0.062 * 2e-4  # Total catalyst area (m2)

# AuPd
rAu = 2.0e-10  # Au-Au radius (m)
aAu = np.pi * rAu ** 2  # Au site area (m2)
Acat_Au = 0.069 * 2e-4  # Total catalyst area (m2)

reactions = [UserDefinedReaction(reac_type='adsorption',
                                 reactants=[states['CO'], states['Clean']],
                                 products=[states['CO.Pd111']],
                                 TS=None,
                                 area=aPd,
                                 name='CO_ads_Pd',
                                 dErxn_user=dEads_CO_Pd,
                                 dGrxn_user=dGads_CO_Pd)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2'], states['Clean'], states['Clean']],
                                  products=[states['O.Pd111'], states['O.Pd111']],  # Two O formed
                                  TS=None,
                                  area=aPd,
                                  name='O2_ads_Pd',
                                  dErxn_user=dEads_O_Pd * 2.0,
                                  dGrxn_user=dGads_O2_Pd)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['CO'], states['Clean']],
                                  products=[states['COonAu']],
                                  TS=None,
                                  area=aAu,
                                  name='CO_ads_AuTop',
                                  dErxn_user=dEads_CO_AuTop,
                                  dGrxn_user=dGads_CO_AuTop)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2'], states['Clean'], states['Clean']],
                                  products=[states['OonAu'], states['OonAu']],
                                  TS=None,
                                  area=aAu,
                                  name='O2_ads_AuTop',
                                  dErxn_user=dEads_O_AuTop * 2.0,
                                  dGrxn_user=dGads_O2_AuTop)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['CO'], states['Clean']],
                                  products=[states['COonPd']],
                                  TS=None,
                                  area=aPd,
                                  name='CO_ads_PdTop',
                                  dErxn_user=dEads_CO_PdTop,
                                  dGrxn_user=dGads_CO_PdTop)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2'], states['Clean'], states['Clean']],
                                  products=[states['OonPd'], states['OonPd']],
                                  TS=None,
                                  area=aPd,
                                  name='O2_ads_PdTop',
                                  dErxn_user=dEads_O_PdTop * 2.0,
                                  dGrxn_user=dGads_O2_PdTop)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['CO'], states['Clean']],
                                  products=[states['COonAu_alloy']],
                                  TS=None,
                                  area=aAu,
                                  name='CO_ads_alloy_Au',
                                  dErxn_user=dEads_CO_alloy,
                                  dGrxn_user=dGads_CO_alloy,
                                  scaling=0.5)]  # Only allow adsorption on half the sites
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2'], states['Clean'], states['Clean']],
                                  products=[states['OonPd_alloy'], states['OonPd_alloy']],
                                  TS=None,
                                  area=aAu,
                                  name='O2_ads_alloy_Pd',
                                  dErxn_user=dEads_O_alloy * 2.0,
                                  dGrxn_user=dGads_O2_alloy,
                                  scaling=0.5)]  # Only allow adsorption on half the sites

# Zip reactions into a dictionary
rnames = [r.name for r in reactions]
reactions = dict(zip(rnames, reactions))

# Transition states - energies from Falsig et al. (2008) scaling relation for CO oxidation
falsig_scaling_coeffs = dict({'gradient': 0.7, 'intercept': 0.02})
snames += ['SRTS.Pd', 'SRTS.AuTop', 'SRTS.PdTop', 'SRTS.alloy']
states['SRTS.Pd'] = Scaling(state_type='TS',
                            name='SRTS.Pd',
                            scaling_coeffs=falsig_scaling_coeffs,
                            scaling_reactions={'CO': {'reaction': reactions['CO_ads_Pd'], 'multiplicity': 1.0},
                                               'O': {'reaction': reactions['O2_ads_Pd'], 'multiplicity': 0.5}},
                            dereference=False)
states['SRTS.AuTop'] = Scaling(state_type='TS',
                               name='SRTS.AuTop',
                               scaling_coeffs=falsig_scaling_coeffs,
                               scaling_reactions={'CO': {'reaction': reactions['CO_ads_AuTop'], 'multiplicity': 1.0},
                                                  'O': {'reaction': reactions['O2_ads_AuTop'], 'multiplicity': 0.5}},
                               dereference=False)
states['SRTS.PdTop'] = Scaling(state_type='TS',
                               name='SRTS.PdTop',
                               scaling_coeffs=falsig_scaling_coeffs,
                               scaling_reactions={'CO': {'reaction': reactions['CO_ads_PdTop'], 'multiplicity': 1.0},
                                                  'O': {'reaction': reactions['O2_ads_PdTop'], 'multiplicity': 0.5}},
                               dereference=False)
states['SRTS.alloy'] = Scaling(state_type='TS',
                               name='SRTS.alloy',
                               scaling_coeffs=falsig_scaling_coeffs,
                               scaling_reactions={'CO': {'reaction': reactions['CO_ads_alloy_Au'], 'multiplicity': 1.0},
                                                  'O': {'reaction': reactions['O2_ads_alloy_Pd'], 'multiplicity': 0.5}},
                               dereference=False)

# Add reactions involving transition states
dEa_CO_ox_Pd = np.max((states['SRTS.Pd'].get_free_energy(T=0, p=p) - dEads_CO_Pd - dEads_O_Pd, 0.0))
dEa_CO_ox_AuTop = np.max((states['SRTS.AuTop'].get_free_energy(T=0, p=p) - dEads_CO_AuTop - dEads_O_AuTop, 0.0))
dEa_CO_ox_PdTop = np.max((states['SRTS.PdTop'].get_free_energy(T=0, p=p) - dEads_CO_PdTop - dEads_O_PdTop, 0.0))
dEa_CO_ox_alloy = np.max((states['SRTS.alloy'].get_free_energy(T=0, p=p) - dEads_CO_alloy - dEads_O_alloy, 0.0))

dGa_CO_ox_Pd = dict(zip(([0] + Ts), [states['SRTS.Pd'].get_free_energy(T=T, p=p) -
                                     (dGads_CO_Pd[T] + 0.5 * dGads_O2_Pd[T])
                                     for T in ([0] + Ts)]))
dGa_CO_ox_AuTop = dict(zip(([0] + Ts), [states['SRTS.AuTop'].get_free_energy(T=T, p=p) -
                                        (dGads_CO_AuTop[T] + 0.5 * dGads_O2_AuTop[T])
                                        for T in ([0] + Ts)]))
dGa_CO_ox_PdTop = dict(zip(([0] + Ts), [states['SRTS.PdTop'].get_free_energy(T=T, p=p) -
                                        (dGads_CO_PdTop[T] + 0.5 * dGads_O2_PdTop[T])
                                        for T in ([0] + Ts)]))
dGa_CO_ox_alloy = dict(zip(([0] + Ts), [states['SRTS.alloy'].get_free_energy(T=T, p=p) -
                                        (dGads_CO_alloy[T] + 0.5 * dGads_O2_alloy[T])
                                        for T in ([0] + Ts)]))

rnames += ['CO_ox_Pd', 'CO_ox_AuTop', 'CO_ox_PdTop', 'CO_ox_alloy']
reactions['CO_ox_Pd'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                            reactants=[states['CO.Pd111'], states['O.Pd111']],
                                            products=[states['CO2'], states['Clean'], states['Clean']],
                                            TS=[states['SRTS.Pd']],
                                            area=aPd,
                                            name='CO_ox_Pd',
                                            dEa_fwd_user=dEa_CO_ox_Pd,
                                            dGa_fwd_user=dGa_CO_ox_Pd,
                                            scaling=1.0)
reactions['CO_ox_AuTop'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                               reactants=[states['COonAu'], states['OonAu']],
                                               products=[states['CO2'], states['Clean'], states['Clean']],
                                               TS=[states['SRTS.AuTop']],
                                               area=aAu,
                                               name='CO_ox_AuTop',
                                               dEa_fwd_user=dEa_CO_ox_AuTop,
                                               dGa_fwd_user=dGa_CO_ox_AuTop,
                                               scaling=1.0)
reactions['CO_ox_PdTop'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                               reactants=[states['COonPd'], states['OonPd']],
                                               products=[states['CO2'], states['Clean'], states['Clean']],
                                               TS=[states['SRTS.PdTop']],
                                               area=aPd,
                                               name='CO_ox_PdTop',
                                               dEa_fwd_user=dEa_CO_ox_PdTop,
                                               dGa_fwd_user=dGa_CO_ox_PdTop,
                                               scaling=1.0)
reactions['CO_ox_alloy'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                               reactants=[states['COonAu_alloy'], states['OonPd_alloy']],
                                               products=[states['CO2'], states['Clean'], states['Clean']],
                                               TS=[states['SRTS.alloy']],
                                               area=aAu,
                                               name='CO_ox_alloy',
                                               dEa_fwd_user=dEa_CO_ox_alloy,
                                               dGa_fwd_user=dGa_CO_ox_alloy,
                                               scaling=1.0)

print('Done.')

# Inflow conditions
inflow_state = dict()
inflow_state['O2'] = 0.08 * p / bartoPa
inflow_state['CO'] = 0.02 * p / bartoPa

# Start conditions
start_state = dict()
start_state['Clean'] = 1.0

flow_rate = (2.4 / 60.0) * 1.0e-6

# Systems
reactionsets = [['CO_ads_Pd', 'O2_ads_Pd', 'CO_ox_Pd'],
                ['CO_ads_AuTop', 'O2_ads_AuTop', 'CO_ox_AuTop'],
                ['CO_ads_PdTop', 'O2_ads_PdTop', 'CO_ox_PdTop'],
                ['CO_ads_alloy_Au', 'O2_ads_alloy_Pd', 'CO_ox_alloy']]
cat_config = ['Pd111', 'AuTop', 'PdTop', 'alloy']
XCO = dict()
ssCO = dict()
ssO = dict()
for v, considered_reactions in enumerate(reactionsets):
    print('Configuring system ' + str(v) + '...')
    if cat_config[v] != 'AuTop':
        reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Pd / 3250)
    else:
        reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Au / 1250)
    sys = System(reactions={key: reactions[key] for key in considered_reactions}, reactor=reactor)

    print('Done.')

    # Solve
    print('Solving ODEs...')

    XCO[cat_config[v]] = []
    ssCO[cat_config[v]] = np.zeros(len(Ts))
    ssO[cat_config[v]] = np.zeros(len(Ts))

    for Tind, T in enumerate(Ts):

        times = np.logspace(start=-14, stop=3, num=int(1e4)) * 3.6
        sys.set_parameters(times=times, inflow_state=inflow_state, start_state=start_state, T=T, p=p,
                           use_jacobian=use_jacobian, verbose=verbose)
        sys.solve_odes()

        # sys.write_results(T=T, p=p, path=results_dir + cat_config[v] + '/')
        # sys.plot_transient(T=T, p=p, path=None)

        XCO[cat_config[v]].append(100.0 - 100.0 * sys.solution[-1][sys.snames.index('CO')] / inflow_state['CO'])
        i = [i for i, s in enumerate(sys.snames) if 'CO' in s and s != 'CO' and s != 'CO2'][0]
        ssCO[cat_config[v]][Tind] = sys.solution[-1][i]
        j = [i for i, s in enumerate(sys.snames) if '2' not in s and 'C' not in s][0]
        ssO[cat_config[v]][Tind] = sys.solution[-1][j]
        print('Tind %1.0f/%1.0f done.' % (Tind + 1, len(Ts)))

Xexp_Pd = [0, 0, 0, 0, 1, 2, 3, 4.5, 6, 38, 46, 52, 52, 52, 52, 52]
Xexp_Au = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 3, 5, 7, 9, 10, 10, 10]
Texp_Pd = np.linspace(start=150, stop=300, endpoint=True, num=16)
Texp_Au = np.linspace(start=150, stop=330, endpoint=True, num=19)

clrs = ['teal', 'tomato', 'royalblue', 'mediumorchid']

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.33, 3.33), sharey='row')
for v, c in enumerate(cat_config):
    if v == 0:
        ax.plot(Texp_Pd, Xexp_Pd, 's', color='k', label='Pd (exp.)')
        ax.plot(Texp_Au, Xexp_Au, '+', color='k', label='PdAu (exp.)')
    ax.plot([T - 273 for T in Ts], XCO[c], color=clrs[v], label=c)
    ax.set(xlabel='Temperature (K)')
ax.set(ylabel='Conversion (%)', ylim=(-1, 101))
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figures_dir + 'conversion_temperature_cmp.png', format='png', dpi=300)
fig.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.33, 3.33), sharey='row')
for v, c in enumerate(cat_config):
    ax.plot([T - 273 for T in Ts], ssCO[c], label=c, color=clrs[v])
    ax.plot([T - 273 for T in Ts], ssO[c], linestyle=':', color=clrs[v], label='')
ax.set(ylabel='Coverage', ylim=(-0.1, 1.1), xlabel='Temperature (K)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figures_dir + 'coverage_temperature_no_rxn.png', format='png', dpi=300)
fig.show()

