from pycatkin.classes.state import *
from pycatkin.classes.reaction import *
from pycatkin.classes.system import *
from pycatkin.classes.reactor import *
from pycatkin.classes.scaling import *
import numpy as np
import pandas as pd

# Conditions
p = 1.0e5  # Pressure (Pa)
Ts = list(np.linspace(start=423, stop=623, num=20, endpoint=True))  # Temperature (K)
times = np.logspace(start=-14, stop=3, num=int(1e4)) * 3.6  # Times (s)
xtol = 1e-12  # SS solver tolerance
use_jacobian = True  # Print messages
verbose = False  # Print messages

# Location of outcars and frequencies
adsdir = 'D:/Users/Astrid/Documents/Chalmers/Data/CO oxidation/Pd111_PdAu_alloys/RPBE/'
gasdir = 'D:/Users/Astrid/Documents/Chalmers/Data/CO oxidation/Pd111_PdAu_alloys/RPBE/'

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/results/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/images/'

print('--------------------')
print('System: CO oxidation')
print('--------------------')

# Load states
print('Configuring states...')

# Adsorbates
states = [State(state_type='adsorbate', path=adsdir + 'CO.Pd111')]
states += [State(state_type='adsorbate', path=adsdir + 'O.Pd111')]
states += [State(state_type='adsorbate', path=adsdir + 'OrderedPdAu-alloy/AuOnTop/COonAu')]
states += [State(state_type='adsorbate', path=adsdir + 'OrderedPdAu-alloy/AuOnTop/OonAu')]
states += [State(state_type='adsorbate', path=adsdir + 'OrderedPdAu-alloy/PdOnTop/COonPd')]
states += [State(state_type='adsorbate', path=adsdir + 'OrderedPdAu-alloy/PdOnTop/OonPd')]
states += [State(state_type='adsorbate', path=adsdir + 'PdAu-alloy/COonAu', name='COonAu_alloy')]
states += [State(state_type='adsorbate', path=adsdir + 'PdAu-alloy/COonPd', name='COonPd_alloy')]
states += [State(state_type='adsorbate', path=adsdir + 'PdAu-alloy/OonAu', name='OonAu_alloy')]
states += [State(state_type='adsorbate', path=adsdir + 'PdAu-alloy/OonPd', name='OonPd_alloy')]

# Gases
states += [State(state_type='gas', path=gasdir + 'O2', sigma=2)]
states += [State(state_type='gas', path=gasdir + 'CO', sigma=1)]
states += [State(state_type='gas', path=gasdir + 'CO2', sigma=2)]

# Surfaces
states += [State(state_type='surface', path=adsdir + 'Clean')]
states += [State(state_type='surface', path=adsdir + 'OrderedPdAu-alloy/AuOnTop/Clean', name='Clean_AuTop')]
states += [State(state_type='surface', path=adsdir + 'OrderedPdAu-alloy/PdOnTop/Clean', name='Clean_PdTop')]
states += [State(state_type='surface', path=adsdir + 'PdAu-alloy/Clean', name='Clean_alloy')]

# Combined adsorbates
# states += [State(state_type='adsorbate', path=adsdir + 'CO_O.Pd111/newstart')]
# states += [State(state_type='adsorbate', path=adsdir + 'OrderedPdAu-alloy/AuOnTop/CO_OonAu')]
# states += [State(state_type='adsorbate', path=adsdir + 'OrderedPdAu-alloy/PdOnTop/CO_OonPd')]
# states += [State(state_type='adsorbate', path=adsdir + 'PdAu-alloy/CO_O_COonAu')]
# states += [State(state_type='adsorbate', path=adsdir + 'PdAu-alloy/CO_O_COonPd')]

# Zip states into a dictionary
snames = [s.name for s in states]
states = dict(zip(snames, states))

# All frequencies are in a sub-folder called freqs
for s in snames:
    if states[s].path is not None:
        states[s].vibs_path = states[s].path + '/nowave'

print('Done.')

# Reactions
print('Configuring reactions...')

# Site properties

# Pd
rPd = 2.0e-10  # Pd-Pd radius (m)
aPd = np.pi * rPd ** 2  # Pd site area (m2)
Acat_Pd = 0.062 * 2e-4  # Total catalyst area (m2)

# AuPd
rAu = 2.0e-10  # Au-Au radius (m)
aAu = np.pi * rAu ** 2  # Au site area (m2)
Acat_Au = 0.069 * 2e-4  # Total catalyst area (m2)

reactions = [Reaction(reac_type='adsorption',
                      reactants=[states['CO'], states['Clean']],
                      products=[states['CO.Pd111']],
                      TS=None,
                      area=aPd,
                      name='CO_ads_Pd')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean'], states['Clean']],
                       products=[states['O.Pd111'], states['O.Pd111']],  # Two O formed
                       TS=None,
                       area=aPd,
                       name='O2_ads_Pd')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CO'], states['Clean_AuTop']],
                       products=[states['COonAu']],
                       TS=None,
                       area=aAu,
                       name='CO_ads_AuTop')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean_AuTop'], states['Clean_AuTop']],
                       products=[states['OonAu'], states['OonAu']],
                       TS=None,
                       area=aAu,
                       name='O2_ads_AuTop')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CO'], states['Clean_PdTop']],
                       products=[states['COonPd']],
                       TS=None,
                       area=aPd,
                       name='CO_ads_PdTop')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean_PdTop'], states['Clean_PdTop']],
                       products=[states['OonPd'], states['OonPd']],
                       TS=None,
                       area=aPd,
                       name='O2_ads_PdTop')]
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CO'], states['Clean_alloy']],
                       products=[states['COonAu_alloy']],
                       TS=None,
                       area=aAu,
                       name='CO_ads_alloy_Au',
                       scaling=0.5)]  # Only allow adsorption on half the sites
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CO'], states['Clean_alloy']],
                       products=[states['COonPd_alloy']],
                       TS=None,
                       area=aAu,
                       name='CO_ads_alloy_Pd',
                       scaling=0.5)]  # Only allow adsorption on half the sites
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean_alloy'], states['Clean_alloy']],
                       products=[states['OonAu_alloy'], states['OonAu_alloy']],
                       TS=None,
                       area=aPd,
                       name='O2_ads_alloy_Au',
                       scaling=0.5)]  # Only allow adsorption on half the sites
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean_alloy'], states['Clean_alloy']],
                       products=[states['OonPd_alloy'], states['OonPd_alloy']],
                       TS=None,
                       area=aPd,
                       name='O2_ads_alloy_Pd',
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
                            dereference=True)
states['SRTS.AuTop'] = Scaling(state_type='TS',
                               name='SRTS.AuTop',
                               scaling_coeffs=falsig_scaling_coeffs,
                               scaling_reactions={'CO': {'reaction': reactions['CO_ads_AuTop'], 'multiplicity': 1.0},
                                                  'O': {'reaction': reactions['O2_ads_AuTop'], 'multiplicity': 0.5}},
                               dereference=True)
states['SRTS.PdTop'] = Scaling(state_type='TS',
                               name='SRTS.PdTop',
                               scaling_coeffs=falsig_scaling_coeffs,
                               scaling_reactions={'CO': {'reaction': reactions['CO_ads_PdTop'], 'multiplicity': 1.0},
                                                  'O': {'reaction': reactions['O2_ads_PdTop'], 'multiplicity': 0.5}},
                               dereference=True)
states['SRTS.alloy'] = Scaling(state_type='TS',
                               name='SRTS.alloy',
                               scaling_coeffs=falsig_scaling_coeffs,
                               scaling_reactions={'CO': {'reaction': reactions['CO_ads_alloy_Au'], 'multiplicity': 1.0},
                                                  'O': {'reaction': reactions['O2_ads_alloy_Pd'], 'multiplicity': 0.5}},
                               dereference=True)

# Add reactions involving transition states
rnames += ['CO_ox_Pd', 'CO_ox_AuTop', 'CO_ox_PdTop', 'CO_ox_alloy']
reactions['CO_ox_Pd'] = Reaction(reac_type='scaling',
                                 reversible=False,
                                 reactants=[states['CO.Pd111'],
                                            states['O.Pd111']],
                                 products=[states['CO2'], states['Clean'], states['Clean']],
                                 TS=[states['SRTS.Pd']],
                                 area=(4.0 * aPd),
                                 name='CO_ox_Pd',
                                 scaling=1.0)
reactions['CO_ox_AuTop'] = Reaction(reac_type='scaling',
                                    reversible=False,
                                    reactants=[states['COonAu'], states['OonAu']],
                                    products=[states['CO2'], states['Clean_AuTop'], states['Clean_AuTop']],
                                    TS=[states['SRTS.AuTop']],
                                    area=(4.0 * aAu),
                                    name='CO_ox_AuTop',
                                    scaling=1.0)
reactions['CO_ox_PdTop'] = Reaction(reac_type='scaling', reversible=False,
                                    reactants=[states['COonPd'], states['OonPd']],
                                    products=[states['CO2'], states['Clean_PdTop'], states['Clean_PdTop']],
                                    TS=[states['SRTS.PdTop']],
                                    area=(4.0 * aPd),
                                    name='CO_ox_PdTop',
                                    scaling=1.0)
reactions['CO_ox_alloy'] = Reaction(reac_type='scaling', reversible=False,
                                    reactants=[states['COonAu_alloy'], states['OonPd_alloy']],
                                    products=[states['CO2'], states['Clean_alloy'], states['Clean_alloy']],
                                    TS=[states['SRTS.alloy']],
                                    area=(4.0 * aAu),
                                    name='CO_ox_alloy',
                                    scaling=1.0)

print('Done.')

# Inflow conditions
inflow_state = dict()
inflow_state['O2'] = 0.08 * p / bartoPa
inflow_state['CO'] = 0.02 * p / bartoPa

# Systems
reactionsets = [['CO_ads_Pd', 'O2_ads_Pd', 'CO_ox_Pd'],
                ['CO_ads_alloy_Au', 'O2_ads_alloy_Pd', 'CO_ox_alloy'],
                ['CO_ads_PdTop', 'O2_ads_PdTop', 'CO_ox_PdTop'],
                ['CO_ads_AuTop', 'O2_ads_AuTop', 'CO_ox_AuTop']]
cat_config = ['Pd(111)', 'Mixed', 'Pd top', 'Au top']
XCO = dict()
ssCO = dict()
ssO = dict()
rCO = dict()
dorc = dict()
Ea = dict()
final_rates = dict()
final_cover = dict()

for v, considered_reactions in enumerate(reactionsets):
    print('Configuring system ' + str(v) + '...')
    if v == 0:
        reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Pd / 3250)
    elif v == 1:
        reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=0.5 * (Acat_Au / 1250 + Acat_Pd / 3250))
    elif v < 3:
        reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Pd / 3250)
    else:
        reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Au / 1250)
    sys = System(reactions={key: reactions[key] for key in considered_reactions}, reactor=reactor)

    cleanname = [[reacs[i].name for i in range(len(reacs)) if 'Clean' in reacs[i].name]
                 for reacs in [reactions[k].reactants for k in considered_reactions]][0][0]

    # Initial conditions
    start_state = dict()
    start_state[cleanname] = 1

    print('Done.')

    # Solve
    print('Solving ODEs...')

    XCO[cat_config[v]] = []
    ssCO[cat_config[v]] = np.zeros(len(Ts))
    ssO[cat_config[v]] = np.zeros(len(Ts))
    rCO[cat_config[v]] = np.zeros(len(Ts))
    final_rates[cat_config[v]] = np.zeros(len(Ts))
    final_cover[cat_config[v]] = np.zeros((len(Ts), 3))
    dorc[cat_config[v]] = dict()
    Ea[cat_config[v]] = np.zeros(len(Ts))

    for Tind, T in enumerate(Ts):

        if Tind > 1 or v != 3:
            sys.set_parameters(times=times, inflow_state=inflow_state, start_state=start_state, T=T, p=p,
                               xtol=xtol, use_jacobian=use_jacobian, verbose=verbose)
        else:
            sys.set_parameters(times=times, inflow_state=inflow_state, start_state=start_state, T=T, p=p,
                               xtol=xtol, rtol=1e-10, atol=1e-12, use_jacobian=use_jacobian, verbose=verbose)
        sys.solve_odes()

        val = 100.0 - 100.0 * sys.solution[-1][sys.snames.index('CO')] / inflow_state['CO']
        XCO[cat_config[v]].append(val if val > 1.0e-9 else 0.0)
        i = [i for i, s in enumerate(sys.snames) if 'CO' in s and s != 'CO' and s != 'CO2'][0]
        ssCO[cat_config[v]][Tind] = sys.solution[-1][i]
        j = [i for i, s in enumerate(sys.snames) if '2' not in s and 'C' not in s][0]
        ssO[cat_config[v]][Tind] = sys.solution[-1][j]
        rCO[cat_config[v]][Tind] = sys.rates[-1][0]

        sys.find_steady(store_steady=True)
        sys.reaction_terms(sys.full_steady)
        final_rates[cat_config[v]][Tind] = sys.rates[-1, 0]
        final_cover[cat_config[v]][Tind] = sys.full_steady[sys.adsorbate_indices]
        tof = final_rates[cat_config[v]][Tind]
        Ea[cat_config[v]][Tind] = 1.0e-3 * np.log((h * tof) / (kB * T)) * (-R * T)
        print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, tof, Ea[cat_config[v]][Tind]))

test = pd.DataFrame(data=XCO, index=Ts)
test.to_csv(results_dir + 'conversion_data.csv')

Xexp_Pd = [0, 0, 0, 0, 1, 2, 3, 4.5, 6, 38, 46, 52, 52, 52, 52, 52]
Xexp_Au = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 3, 5, 7, 9, 10, 10, 10]
Texp_Pd = np.linspace(start=150, stop=300, endpoint=True, num=16)
Texp_Au = np.linspace(start=150, stop=330, endpoint=True, num=19)

clrs = ['teal', 'mediumorchid', 'royalblue', 'tomato', 'mediumorchid']

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.2, 3.2), sharey='row')
for v, c in enumerate(cat_config):
    if v < 4:
        ax.plot([T - 273 for T in Ts], XCO[c], label=c, color=clrs[v])
    else:
        ax.plot([T - 273 for T in Ts], XCO[c], label='', color=clrs[v], linestyle=':')
    ax.set(xlabel=r'Temperature ($^\circ\mathsf{ C}$)')
ax.set(ylabel='Conversion (%)', ylim=(-2, 61), xlim=(140, 360))
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figures_dir + 'conversion_temperature.png', format='png', dpi=600)
fig.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.2, 3.2), sharey='row')
for v, c in enumerate(cat_config):
    if v < 4:
        ax.plot([T - 273 for T in Ts], rCO[c], label=c, color=clrs[v])
    else:
        ax.plot([T - 273 for T in Ts], rCO[c], label='', color=clrs[v], linestyle=':')
    ax.set(xlabel=r'Temperature ($^\circ\mathsf{ C}$)')
ax.set(ylabel='TOF (1/s)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figures_dir + 'tof_temperature_cmp.png', format='png', dpi=600)
fig.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.2, 3.2), sharey='row')
for v, c in enumerate(cat_config):
    ax.plot([T - 273 for T in Ts], ssCO[c], label=c, color=clrs[v])
    ax.plot([T - 273 for T in Ts], ssO[c], linestyle=':', label='', color=clrs[v])
ax.set(ylabel='Coverage', ylim=(-0.05, 1.05),
       xlabel=r'Temperature ($^\circ\mathsf{ C}$)', xlim=(140, 360))
ax.legend(frameon=False)
ax.text(160, 0.85, 'O', color='mediumorchid', ha='left', va='top')
ax.text(160, 0.07, 'CO', color='mediumorchid', ha='left', va='bottom')
fig.tight_layout()
fig.savefig(figures_dir + 'coverage_temperature.png', format='png', dpi=600)
fig.show()

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6.4, 3.2))
for v, s in enumerate([states['Clean_alloy'], states['Clean_PdTop'], states['Clean_AuTop']]):
    ax[0].bar(v, s.Gelec - states['Clean'].Gelec, label=cat_config[v], facecolor=clrs[v])
    ax[1].plot([T - 273 for T in Ts],
               [s.get_free_energy(T, p=p) - states['Clean'].get_free_energy(T, p=p) for T in Ts],
               label=cat_config[v + 1], color=clrs[v])
    ax[0].text(v, s.Gelec - states['Clean'].Gelec - 0.5, '%1.1f eV' % (s.Gelec - states['Clean'].Gelec),
               color='white', fontweight='bold', rotation=90, va='top', ha='center')
ax[0].set(ylabel='Electronic energy (eV)', xlabel=r'Surface model',
          xticks=np.arange(3), xticklabels=cat_config[1::])
ax[1].set(ylabel='Free energy (eV)', xlabel=r'Temperature ($^\circ\mathsf{ C}$)')
ax[1].legend(frameon=False)
fig.tight_layout()
fig.savefig(figures_dir + 'surface_energies.png', format='png', dpi=600)
fig.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.2, 3.2), sharey='row')
for v, c in enumerate(cat_config):
    if v == 0:
        ax.plot(Texp_Pd, Xexp_Pd, 'o', label='Exp',
                markerfacecolor='white', color='teal')
        ax.plot(Texp_Au, Xexp_Au, 'o', label='',
                markerfacecolor='white', color='tomato')
        ax.plot([T - 273 for T in Ts], XCO[c], label='Model', color=clrs[v])
    elif v == 3:
        ax.plot([T - 273 for T in Ts], XCO[c], label='', color=clrs[v])
    ax.set(xlabel=r'Temperature ($^\circ\mathsf{ C}$)', xlim=(140, 360),
           ylabel='Conversion (%)', ylim=(-2, 61))
ax.text(325, 51, 'Pd', color='teal', ha='left', va='top')
ax.text(325, 7, 'AuPd', color='tomato', ha='left', va='top')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figures_dir + 'conversion_temperature_cmp.png', format='png', dpi=600)
fig.show()

fig, ax = plt.subplots(figsize=(3.2, 3.2))
ax.plot(Texp_Pd, Xexp_Pd, '-s', label='Pd',
        color=(0, 108 / 255, 137 / 255))
ax.plot(Texp_Au, Xexp_Au, '-o', label='AuPd',
        markerfacecolor=(1, 230 / 255, 41 / 255), color=(0, 108 / 255, 137 / 255))
ax.set(xlabel=r'Temperature ($^\circ\mathsf{ C}$)', xlim=(140, 355),
       ylabel='Conversion (%)', ylim=(-2, 61))
ax.text(302, 49, 'Pd', ha='left', va='top')
ax.text(302, 6, 'PdAu', ha='left', va='top')
fig.tight_layout()
fig.savefig(figures_dir + 'exp_conversion.png', format='png', dpi=600)
fig.show()
