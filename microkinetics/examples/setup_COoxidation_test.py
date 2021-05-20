from microkinetics.classes.state import *
from microkinetics.classes.reaction import *
from microkinetics.classes.system import *
from microkinetics.classes.reactor import *
from microkinetics.classes.scaling import *
from microkinetics.classes.energy import *
import numpy as np
import time

# Conditions
p = 1.01325e5  # Pressure (Pa)
Ts = list(np.linspace(start=423, stop=623, num=20, endpoint=True))  # Temperature (K)
times = np.logspace(start=-10, stop=3, num=int(1e4))  # Times (s)
withflow = False  # Include reactor model
use_jacobian = True  # Use Jacobian to solve SS and ODEs
verbose = True  # Print messages

# Location of outcars and frequencies
adsdir = 'D:/Users/Astrid/Documents/Chalmers/Data/CO oxidation/Pd111_PdAu_alloys/RPBE/'
gasdir = 'D:/Users/Astrid/Documents/Chalmers/Data/CO oxidation/Pd111_PdAu_alloys/RPBE/'

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/results/toy/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/images/toy/'

print('--------------------')
print('System: CO oxidation')
print('--------------------')

# Load states
print('Configuring states...')

# Adsorbates
states = [State(state_type='adsorbate', path=adsdir + 'CO.Pd111')]
states += [State(state_type='adsorbate', path=adsdir + 'O.Pd111')]

# Gases
states += [State(state_type='gas', path=gasdir + 'O2', sigma=2)]
states += [State(state_type='gas', path=gasdir + 'CO', sigma=1)]
states += [State(state_type='gas', path=gasdir + 'CO2', sigma=2)]

# Surfaces
states += [State(state_type='surface', path=adsdir + 'Clean')]

# Zip states into a dictionary
snames = [s.name for s in states]
states = dict(zip(snames, states))

# All frequencies are in a sub-folder called freqs
for s in snames:
    if states[s].path is not None:
        states[s].vibs_path = states[s].path + '/nowave'

print('Done.')

# Site properties

# Pd
rPd = 2.0e-10  # Pd-Pd radius (m)
aPd = np.pi * rPd ** 2  # Pd site area (m2)
Acat_Pd = 0.062 * 2e-4 / 3250  # Total catalyst area (m2)

# Reactions
print('Configuring reactions...')

reactions = [Reaction(reac_type='adsorption', name='CO_ads_Pd', area=aPd,
                      reactants=[states['CO'], states['Clean']], TS=None,
                      products=[states['CO.Pd111']])]
reactions += [Reaction(reac_type='adsorption', name='O2_ads_Pd', area=aPd,
                       reactants=[states['O2'], states['Clean'], states['Clean']],
                       TS=None,
                       products=[states['O.Pd111'], states['O.Pd111']])]  # Two O formed

# Zip reactions into a dictionary
rnames = [r.name for r in reactions]
reactions = dict(zip(rnames, reactions))

# Transition states - energies from Falsig et al. (2008) scaling relation for CO oxidation
falsig_scaling_coeffs = dict({'gradient': 0.7, 'intercept': 0.02})
snames += ['SRTS.Pd']
states['SRTS.Pd'] = Scaling(state_type='TS',
                            name='SRTS.Pd',
                            scaling_coeffs=falsig_scaling_coeffs,
                            scaling_reactions={'CO': {'reaction': reactions['CO_ads_Pd'], 'multiplicity': 1.0},
                                               'O': {'reaction': reactions['O2_ads_Pd'], 'multiplicity': 0.5}},
                            dereference=True)

# Add reactions involving transition states
rnames += ['CO_ox_Pd']
reactions['CO_ox_Pd'] = Reaction(reac_type='scaling',
                                 reversible=False,
                                 reactants=[states['CO.Pd111'],
                                            states['O.Pd111']],
                                 products=[states['CO2'], states['Clean'], states['Clean']],
                                 TS=[states['SRTS.Pd']],
                                 area=aPd,
                                 name='CO_ox_Pd',
                                 scaling=1.0)

print('Done.')

# Reactor
print('Configuring reactor...')

if not withflow:
    # InfiniteDilutionReactor
    start_state = dict()
    start_state['O2'] = 0.8 * 0.02 * p / bartoPa
    start_state['CO'] = 0.2 * 0.02 * p / bartoPa
    start_state['Clean'] = 1.0
    inflow_state = None
    reactor = InfiniteDilutionReactor()
else:
    # CSTR
    start_state = None
    inflow_state = dict()
    inflow_state['O2'] = 0.8 * p / bartoPa
    inflow_state['CO'] = 0.2 * p / bartoPa
    start_state['Clean'] = 1.0
    reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Pd)

sys = System(reactions=reactions, reactor=reactor)

XCO = []
ssCO = np.zeros(len(Ts))
ssO = np.zeros(len(Ts))

final_rates = np.zeros((len(Ts), len(reactions)))
final_cover = np.zeros((len(Ts), len([s for s in snames
                                      if states[s].state_type == 'adsorbate'
                                      or states[s].state_type == 'surface'])))
dorc = dict()

print('Done.')

# Solve
print('Solving ODEs...')

for Tind, T in enumerate(Ts):
    sys.set_parameters(times=times, start_state=start_state, inflow_state=inflow_state, T=T, p=p,
                       xtol=1.0e-12, use_jacobian=use_jacobian, verbose=verbose)
    start_time = time.time()
    sys.solve_odes()
    print("--- %s seconds ---" % (time.time() - start_time))
    # sys.plot_transient(T=T, p=p, path=figures_dir)
    # sys.write_results(T=T, p=p, path=results_dir)

    val = 100.0 - 100.0 * sys.solution[-1][sys.snames.index('CO')] / inflow_state['CO'] if inflow_state else 0.0
    XCO.append(val if val > 1.0e-9 else 0.0)
    ssCO[Tind] = sys.solution[-1][sys.snames.index('CO.Pd111')]
    ssO[Tind] = sys.solution[-1][sys.snames.index('O.Pd111')]

    sys.find_steady(store_steady=True)  # , path=figures_dir
    sys.reaction_terms(sys.full_steady)

    final_rates[Tind, :] = sys.rates[:, 0] - sys.rates[:, 1]
    final_cover[Tind, :] = sys.full_steady[sys.adsorbate_indices]
    tof = final_rates[Tind, -1]

    tof_terms = ['CO_ox_Pd']
    dorc[T] = sys.degree_of_rate_control(tof_terms, ss_solve=False, eps=1.0e-3)

    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, tof, 1.0e-3 *
                                                               np.log((h * tof) / (kB * T)) * (-R * T)))


fig, ax = plt.subplots(figsize=(3.2, 3.2))
ax.plot([T - 273 for T in Ts], XCO, color='slategrey')
ax.set(xlabel=r'Temperature ($^\circ\mathsf{ C}$)')
ax.set(ylabel='CO conversion (%)')
fig.tight_layout()
fig.savefig(figures_dir + 'conversion_temperature.png', format='png', dpi=300)
fig.show()

fig, ax = plt.subplots(figsize=(3.2, 3.2))
ax.plot(Ts, final_rates[:, -1], color='teal')
ax.set(xlabel=r'Temperature (K)')
ax.set(ylabel='TOF (1/s)', yscale='log')
fig.tight_layout()
fig.savefig(figures_dir + 'tof_temperature.png', format='png', dpi=300)
fig.show()

fig, ax = plt.subplots(figsize=(3.2, 3.2))
ax.plot(Ts, ssCO, label='CO', color='slategrey')
ax.plot(Ts, ssO, label='O', color='darkorchid')
ax.set(ylabel='Coverage', ylim=(-0.01, 1.01), xlabel=r'Temperature (K)')
ax.legend(frameon=False, loc='center right')
fig.tight_layout()
fig.savefig(figures_dir + 'coverage_temperature.png', format='png', dpi=300)
fig.show()

fig, ax = plt.subplots(figsize=(3.2, 3.2))
ax.plot(Ts, [dorc[T]['CO_ads_Pd'] for T in Ts], label='CO ads', color='slategrey')
ax.plot(Ts, [dorc[T]['O2_ads_Pd'] for T in Ts], label='O ads', color='darkorchid')
ax.plot(Ts, [dorc[T]['CO_ox_Pd'] for T in Ts], label='CO ox', color='teal')
ax.plot(Ts, [sum([dorc[T][k] for k in dorc[T].keys()]) for T in Ts], label='', color='k', linestyle=':')
ax.set(ylabel='DORC', xlabel=r'Temperature (K)')
ax.legend(frameon=False, loc='center right')
fig.tight_layout()
fig.savefig(figures_dir + 'dorc_temperature.png', format='png', dpi=300)
fig.show()

minima = dict()
minima[0] = [states['Clean'], states['Clean'], states['Clean'], states['Clean'],
             states['O2'], states['CO'], states['CO']]
minima[1] = [states['Clean'], states['Clean'], states['Clean'],
             states['CO.Pd111'], states['O2'], states['CO']]
minima[2] = [states['Clean'], states['Clean'], states['O2'],
             states['CO.Pd111'], states['CO.Pd111']]
minima[3] = [states['CO.Pd111'], states['CO.Pd111'],
             states['O.Pd111'], states['O.Pd111']]
minima[4] = [states['SRTS.Pd'], states['CO.Pd111'], states['O.Pd111']]
minima[5] = [states['CO2'], states['CO.Pd111'], states['O.Pd111'],
             states['Clean'], states['Clean']]
minima[6] = [states['SRTS.Pd'], states['CO2'], states['Clean'], states['Clean']]
minima[7] = [states['CO2'], states['CO2'],
             states['Clean'], states['Clean'], states['Clean'], states['Clean']]

energy = Energy(minima)
energy.draw_energy_landscape(T=500, p=p, etype='electronic', path=figures_dir)
energy.draw_energy_landscape(T=500, p=p, etype='free', path=figures_dir)

tofs = []
fig, ax = plt.subplots(figsize=(3.2, 3.2))
for T in Ts:
    tof, num_i, num_j, lTi, lIj = energy.evaluate_energy_span_model(T=T, p=p, verbose=False, etype='free')
    tofs.append(tof)
    print('ES:')
    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, tof, 1.0e-3 *
                                                               np.log((h * tof * 2.0) / (kB * T)) * (-R * T)))
ax.plot(Ts, tofs, color='darkorchid', label='ES v1')
ax.plot(Ts, final_rates[:, -1], '-', color='teal', label='MK')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='TOF (1/s)')
ax.legend(frameon=False, loc='center right')
fig.tight_layout()
fig.savefig(figures_dir + 'EStof_temperature.png', format='png', dpi=300)
fig.show()

minima = dict()
minima[0] = [states['Clean'], states['Clean'], states['Clean'], states['O2'], states['CO'], states['CO']]
minima[1] = [states['Clean'], states['Clean'], states['CO.Pd111'], states['O2'], states['CO']]
minima[2] = [states['CO.Pd111'], states['O.Pd111'], states['O.Pd111'], states['CO']]
minima[3] = [states['SRTS.Pd'], states['O.Pd111'], states['CO']]
minima[4] = [states['CO2'], states['O.Pd111'], states['Clean'], states['Clean'], states['CO']]
minima[5] = [states['CO2'], states['O.Pd111'], states['CO.Pd111'], states['Clean']]
minima[6] = [states['CO2'], states['SRTS.Pd'], states['Clean']]
minima[7] = [states['CO2'], states['CO2'], states['Clean'], states['Clean'], states['Clean']]

energy = Energy(minima)
energy.draw_energy_landscape(T=500, p=p, etype='electronic', path=figures_dir)
energy.draw_energy_landscape(T=500, p=p, etype='free', path=figures_dir)

tofs = []
# fig, ax = plt.subplots(figsize=(3.2, 3.2))
for T in Ts:
    tof, num_i, num_j, lTi, lIj = energy.evaluate_energy_span_model(T=T, p=p, verbose=False, etype='free')
    tofs.append(tof)
    print('ES:')
    print('T = %1.0f K, TOF = %1.2e 1/s, Ea = %1.2f kJ/mol' % (T, tof, 1.0e-3 *
                                                               np.log((h * tof * 2.0) / (kB * T)) * (-R * T)))
ax.plot(Ts, tofs, '--', color='orchid', label='ES v2')
ax.set(yscale='log', xlabel='Temperature (K)', ylabel='TOF (1/s)')
ax.legend(frameon=False, loc='center right')
fig.tight_layout()
fig.savefig(figures_dir + 'EStof_temperature.png', format='png', dpi=300)
fig.show()

