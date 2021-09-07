from pycatkin.classes.state import *
from pycatkin.classes.reaction import *
from pycatkin.classes.system import *
from pycatkin.classes.reactor import *
from pycatkin.classes.scaling import *
import numpy as np
import time

# Conditions
p = 1.0e5  # Pressure (Pa)
T = 600
times = np.logspace(start=-14, stop=2, num=int(1e3))
# times = np.concatenate((np.zeros(1), np.logspace(start=-14, stop=3, num=int(1e5))))
withflow = False  # Choose reactor model
use_jacobian = True  # Use Jacbian for ODE and SS solvers
verbose = False  # Print messages

# Location of results files and images
results_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/results/surface/'
figures_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/COoxidation/images/surface/'

print('----------------------------')
print('System: CO oxidation volcano')
print('----------------------------')

# Catalyst properties
rPd = 1.0e-10  # Pd-Pd radius (m)
aPd = np.pi * rPd ** 2  # Pd site area (m2)
Acat_Pd = 0.062 * 2e-4 / 1250  # Total catalyst area (m2)

# Standard entropies (Atkins, J/molK)
SCOg = 197.67 * 1.0e-3 / eVtokJ
SO2g = 205.138 * 1.0e-3 / eVtokJ

if not withflow:
    # Conditions - IDR
    inflow_state = None
    start_state = dict()
    start_state['O2'] = 0.33 * p / bartoPa
    start_state['CO'] = 0.67 * p / bartoPa
    start_state['*'] = 1
    reactor = InfiniteDilutionReactor()
else:
    # Conditions - CSTR
    start_state = dict()
    start_state['*'] = 1
    inflow_state = dict()
    inflow_state['O2'] = 0.08 * p / bartoPa
    inflow_state['CO'] = 0.02 * p / bartoPa
    reactor = CSTReactor(residence_time=4.5, volume=180.0e-9, catalyst_area=Acat_Pd)

# Systems
npts = 10
besCO = np.linspace(start=-2.5, stop=0.5, num=npts, endpoint=True)
besO = np.linspace(start=-2.5, stop=0.5, num=npts, endpoint=True)
act = np.zeros((npts, npts))
tofs = np.zeros((npts, npts))

# User defined states
states = [State(state_type='adsorbate', name='CO*')]
states += [State(state_type='adsorbate', name='O*')]
states += [State(state_type='adsorbate', name='O2*')]
states += [State(state_type='surface', name='*')]
states += [State(state_type='gas', name='O2', sigma=2, mass=32)]
states += [State(state_type='gas', name='CO', sigma=1, mass=28)]
states += [State(state_type='gas', name='CO2', sigma=2, mass=44)]
snames = [s.name for s in states]
states = dict(zip(snames, states))

# Transition states - energies from Falsig et al. (2008) scaling relation for CO oxidation
falsig_scaling_coeffs = dict({'gradient': 0.7, 'intercept': 0.02})
falsig_scaling_coeffs_OTS = dict({'gradient': 1.39, 'intercept': 1.56})
falsig_scaling_coeffs_O2O = dict({'gradient': 0.89, 'intercept': 0.17})

start_time = time.time()

for iCO, bCO in enumerate(besCO):
    for iO, bO in enumerate(besO):
        bO2 = 0.89 * bO + 0.17
        dGads_CO = dict(zip(([0, T]), [bCO + SCOg * Ti for Ti in [0, T]]))
        dGads_O2 = dict(zip(([0, T]), [bO2 + SO2g * Ti for Ti in [0, T]]))
        dGads_O = dict(zip(([0, T]), [bO * 2.0 + SO2g * Ti for Ti in [0, T]]))

        # Reactions
        reactions = [UserDefinedReaction(reac_type='adsorption',
                                         reactants=[states['CO'], states['*']],
                                         products=[states['CO*']],
                                         TS=None,
                                         area=aPd,
                                         name='CO_ads',
                                         dErxn_user=bCO,
                                         dGrxn_user=dGads_CO)]
        reactions += [UserDefinedReaction(reac_type='adsorption',
                                          reactants=[states['O2'], states['*']],
                                          products=[states['O2*']],
                                          TS=None,
                                          area=aPd,
                                          name='O2_ads',
                                          dErxn_user=bO2,
                                          dGrxn_user=dGads_O2)]
        # Zip reactions into a dictionary
        rnames = [r.name for r in reactions]
        reactions = dict(zip(rnames, reactions))
        # Add scaling relation
        tmpreaction = dict()
        tmpreaction['O_ads'] = UserDefinedReaction(reac_type='adsorption',
                                                   reactants=[states['O2'], states['*'], states['*']],
                                                   products=[states['O*'], states['O*']],
                                                   TS=None,
                                                   area=aPd,
                                                   name='O_ads',
                                                   dErxn_user=bO * 2.0,
                                                   dGrxn_user=dGads_O)
        snames += ['SRTS_O']
        states['SRTS_O'] = Scaling(state_type='TS',
                                   name='SRTS_O',
                                   scaling_coeffs=falsig_scaling_coeffs_OTS,
                                   scaling_reactions={'O': {'reaction': tmpreaction['O_ads'], 'multiplicity': 0.5}},
                                   dereference=False)
        dEa_O2O = np.max((states['SRTS_O'].get_free_energy(T=0, p=p) - bO2, 0.0))

        rnames += ['O2_diss']
        reactions['O2_diss'] = UserDefinedReaction(reac_type='scaling',
                                                   reversible=False,
                                                   reactants=[states['O2*'], states['*']],
                                                   products=[states['O*'], states['O*']],
                                                   TS=[states['SRTS_O']],
                                                   area=aPd,
                                                   name='O2_diss',
                                                   dEa_fwd_user=dEa_O2O,
                                                   dGa_fwd_user=dEa_O2O)

        snames += ['SRTS']
        states['SRTS'] = Scaling(state_type='TS',
                                 name='SRTS',
                                 scaling_coeffs=falsig_scaling_coeffs,
                                 scaling_reactions={'CO': {'reaction': reactions['CO_ads'], 'multiplicity': 1.0},
                                                    'O': {'reaction': tmpreaction['O_ads'], 'multiplicity': 0.5}},
                                 dereference=False)
        dEa_COox = np.max((states['SRTS'].get_free_energy(T=0, p=p) - bCO - bO, 0.0))
        dGa_COox = dict(zip(([0, T]), [states['SRTS'].get_free_energy(T=Ti, p=p) - (dGads_CO[Ti] + 0.5 * dGads_O[Ti])
                                       for Ti in [0, T]]))
        rnames += ['COox']
        reactions['COox'] = UserDefinedReaction(reac_type='scaling',
                                                reversible=False,
                                                reactants=[states['CO*'], states['O*']],
                                                products=[states['CO2'], states['*'], states['*']],
                                                TS=[states['SRTS']],
                                                area=aPd,
                                                name='COox',
                                                dEa_fwd_user=dEa_COox,
                                                dGa_fwd_user=dGa_COox)

        print('Energy comp: %1.2f eV, %1.2f eV, %1.2f eV' % (dEa_COox, dGa_COox[0], dGa_COox[T]))

        # Let's get started...
        print('Configuring system with bCO=%1.2f eV and bO=%1.2f eV...' % (bCO, bO))
        sys = System()
        for s in states.keys():
            sys.add_state(state=states[s])
        for r in reactions.keys():
            sys.add_reaction(reaction=reactions[r])
        sys.add_reactor(reactor=reactor)
        sys.names_to_indices()

        # Solve
        print('Solving ODEs...')
        sys.set_parameters(times=times, start_state=start_state, inflow_state=inflow_state, T=T, p=p,
                           use_jacobian=use_jacobian, verbose=verbose)
        sys.solve_odes()

        for r in sys.reactions.keys():
            print(r)
            try:
                print('dGr = %1.2f' % (sys.reactions[r].get_reaction_energy(T=600, p=1e5, etype='electronic') * 1e-3 / eVtokJ))
            except:
                print('irreversible')
            try:
                print('dGa = %1.2f' % (sys.reactions[r].get_reaction_barriers(T=600, p=1e5, etype='electronic')[0] * 1e-3 / eVtokJ))
            except:
                print('barrierless')

        print('Finding steady-state solution...')
        # sys.find_steady(store_steady=True)
        # sys.reaction_terms(sys.full_steady)

#         tof = sys.rates[-1, 0]
#         if tof < 1.0e-40:
#             tof = 1.0e-40
#         tofs[iCO, iO] = tof
#         act[iCO, iO] = 1.0e-3 * np.log((h * tof) / (kB * T)) * (R * T) / eVtokJ
#         print('Done.')
#
# print("--- %s seconds ---" % (time.time() - start_time))
#
# ending = ''
# if isinstance(reactor, CSTReactor):
#     ending = '_cstr'
#
# fig, ax = plt.subplots(figsize=(4.3, 3.2))
# CS = ax.contourf(besO, besCO, act, levels=25, cmap=plt.cm.RdYlBu_r)
# cbar = fig.colorbar(CS)
# cbar.ax.set_ylabel('Activity (eV)')
# ax.set(ylabel=r'$E_{\mathsf{CO}}$ (eV)', xlabel=r'$E_{\mathsf{O}}$ (eV)',
#        xlim=(-2.5, 0.5), ylim=(-2.5, 0.5),
#        xticks=(-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5),
#        yticks=(-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5))
# fig.tight_layout()
# fig.savefig(figures_dir + 'be_activity' + ending + '.png', format='png', dpi=600)
# fig.show()
