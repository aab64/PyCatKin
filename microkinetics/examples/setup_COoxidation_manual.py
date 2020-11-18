from microkinetics.classes.state import *
from microkinetics.classes.reaction import *
from microkinetics.classes.system import *
from microkinetics.classes.reactor import *
from microkinetics.classes.scaling import *

# Conditions
p = 1.0e5  # Pressure (Pa)
verbose = True  # Print messages
savexyz = False  # Save xyz files (not invoked atm)

# Location of outcars and frequencies
adsdir = 'D:/Users/Astrid/Dropbox/COox.Pd111-Astrid/RPBE/'
gasdir = 'D:/Users/Astrid/Dropbox/COox.Pd111-Astrid/RPBE/'

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

# Gases
states += [State(state_type='gas', name='O2', sigma=2, mass=32)]  # Specify mass for adsorption
states += [State(state_type='gas', name='CO', sigma=1, mass=28)]
states += [State(state_type='gas', name='CO2', sigma=2)]

# Zip states into a dictionary
snames = [s.name for s in states]
states = dict(zip(snames, states))

print('Done.')

# Reactions
print('Configuring reactions...')

# Electronic energies from DFT
dGads_CO_Pd = -1.18
dGads_O_Pd = -1.15
dGads_CO_alloy = -1.34
dGads_O_alloy = -0.49
dGads_CO_AuTop = -0.093
dGads_O_AuTop = 0.35
dGads_CO_PdTop = -1.14
dGads_O_PdTop = -1.15

# Site properties

# Pd
rPd = 1.63e-10  # Pd vdw radius (m)
aPd = np.pi * rPd ** 2  # Pd site area (m2)
Acat_Pd = 0.062 * 2e-4  # Total catalyst area (m2)

# AuPd
rAu = 1.66e-10  # Au vdw radius (m)
aAu = np.pi * rAu ** 2  # Au site area (m2)
Acat_Au = 0.069 * 2e-4  # Total catalyst area (m2)

# Temporary changes: fix scaling too!
# rAu = 1.63e-10
# aAu = np.pi * rPd ** 2
# Acat_Au = 0.062 * 2e-4

reactions = [UserDefinedReaction(reac_type='adsorption',
                                 reactants=[states['CO']],
                                 products=[states['CO.Pd111']],
                                 TS=None,
                                 area=aPd,
                                 name='CO_ads_Pd',
                                 dGrxn_user=dGads_CO_Pd)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2']],
                                  products=[states['O.Pd111'], states['O.Pd111']],  # Two O formed
                                  TS=None,
                                  area=aPd,
                                  name='O2_ads_Pd',
                                  dGrxn_user=dGads_O_Pd * 2.0)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['CO']],
                                  products=[states['COonAu']],
                                  TS=None,
                                  area=aAu,
                                  name='CO_ads_AuTop',
                                  dGrxn_user=dGads_CO_AuTop)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2']],
                                  products=[states['OonAu'], states['OonAu']],
                                  TS=None,
                                  area=aAu,
                                  name='O2_ads_AuTop',
                                  dGrxn_user=dGads_O_AuTop * 2.0)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['CO']],
                                  products=[states['COonPd']],
                                  TS=None,
                                  area=aPd,
                                  name='CO_ads_PdTop',
                                  dGrxn_user=dGads_CO_PdTop)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2']],
                                  products=[states['OonPd'], states['OonPd']],
                                  TS=None,
                                  area=aPd,
                                  name='O2_ads_PdTop',
                                  dGrxn_user=dGads_O_PdTop * 2.0)]
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['CO']],
                                  products=[states['COonAu_alloy']],
                                  TS=None,
                                  area=aAu,
                                  name='CO_ads_alloy_Au',
                                  dGrxn_user=dGads_CO_alloy,
                                  scaling=0.5)]  # Only allow adsorption on half the sites
reactions += [UserDefinedReaction(reac_type='adsorption',
                                  reactants=[states['O2']],
                                  products=[states['OonPd_alloy'], states['OonPd_alloy']],
                                  TS=None,
                                  area=aPd,
                                  name='O2_ads_alloy_Pd',
                                  dGrxn_user=dGads_O_alloy * 2.0,
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
dGa_CO_ox_Pd = states['SRTS.Pd'].get_free_energy(T=0, p=p) - dGads_CO_Pd - dGads_O_Pd
dGa_CO_ox_AuTop = states['SRTS.AuTop'].get_free_energy(T=0, p=p) - dGads_CO_AuTop - dGads_O_AuTop
dGa_CO_ox_PdTop = states['SRTS.PdTop'].get_free_energy(T=0, p=p) - dGads_CO_PdTop - dGads_O_PdTop
dGa_CO_ox_alloy = states['SRTS.alloy'].get_free_energy(T=0, p=p) - dGads_CO_alloy - dGads_O_alloy

rnames += ['CO_ox_Pd', 'CO_ox_AuTop', 'CO_ox_PdTop', 'CO_ox_alloy']
reactions['CO_ox_Pd'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                            reactants=[states['CO.Pd111'], states['O.Pd111']],
                                            products=[states['CO2']],
                                            TS=states['SRTS.Pd'],
                                            area=aPd,
                                            name='CO_ox_Pd',
                                            dGa_fwd_user=dGa_CO_ox_Pd)
reactions['CO_ox_AuTop'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                               reactants=[states['COonAu'], states['OonAu']],
                                               products=[states['CO2']],
                                               TS=states['SRTS.AuTop'],
                                               area=aAu,
                                               name='CO_ox_AuTop',
                                               dGa_fwd_user=dGa_CO_ox_AuTop)
reactions['CO_ox_PdTop'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                               reactants=[states['COonPd'], states['OonPd']],
                                               products=[states['CO2']],
                                               TS=states['SRTS.PdTop'],
                                               area=aPd,
                                               name='CO_ox_PdTop',
                                               dGa_fwd_user=dGa_CO_ox_PdTop)
reactions['CO_ox_alloy'] = UserDefinedReaction(reac_type='scaling', reversible=False,
                                               reactants=[states['COonAu_alloy'], states['OonPd_alloy']],
                                               products=[states['CO2']],
                                               TS=states['SRTS.alloy'],
                                               area=0.5 * (aPd + aAu),
                                               name='CO_ox_alloy',
                                               dGa_fwd_user=dGa_CO_ox_alloy)

print('Done.')

# Inflow conditions
inflow_state = dict()
inflow_state['O2'] = 0.08 * p / bartoPa
inflow_state['CO'] = 0.02 * p / bartoPa

# Temperatures
Ts = np.linspace(start=373, stop=673, num=61, endpoint=True)

# Systems
reactionsets = [['CO_ads_Pd', 'O2_ads_Pd', 'CO_ox_Pd'],
                ['CO_ads_AuTop', 'O2_ads_AuTop', 'CO_ox_AuTop'],
                ['CO_ads_PdTop', 'O2_ads_PdTop', 'CO_ox_PdTop'],
                ['CO_ads_alloy_Au', 'O2_ads_alloy_Pd', 'CO_ox_alloy']]
cat_config = ['Pd111', 'AuTop', 'PdTop', 'alloy']
for v, considered_reactions in enumerate(reactionsets):
    print('Configuring system ' + str(v) + '...')
    if v == 0:
        reactor = CSTReactor(flow_rate=(2.4 / 60) * 1.0e-6, volume=180.0e-9, catalyst_area=Acat_Pd)
    else:
        reactor = CSTReactor(flow_rate=(2.4 / 60) * 1.0e-6, volume=180.0e-9, catalyst_area=Acat_Au)
    sys = System(reactions={key: reactions[key] for key in considered_reactions}, reactor=reactor)

    print('Done.')

    # Solve
    print('Solving ODEs...')

    for Tind, T in enumerate(Ts):

        times = np.logspace(start=-16, stop=3, num=int(1e2)) * 2.4
        sys.run_odeint(times=times, inflow_state=inflow_state, T=T, p=p, verbose=verbose)
        print('Done.')

        print('Post-processing...')
        sys.write_results(T=T, p=p, path=results_dir + cat_config[v] + '/')
        # sys.plot_transient(T=T, p=p, path=None)
