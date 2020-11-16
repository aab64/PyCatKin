from state import *
from reaction import *
from system import *
from reactor import *
from scaling import *
import numpy as np

# Conditions
p = 1.0e5  # Pressure (Pa)
verbose = True  # Print messages
savexyz = False  # Save xyz files (not invoked atm)

# Location of outcars and frequencies
adsdir = 'D:/Users/Astrid/Documents/Chalmers/Data/CO oxidation/Pd111_PdAu_alloys/RPBE/'
gasdir = 'D:/Users/Astrid/Documents/Chalmers/Data/CO oxidation/Pd111_PdAu_alloys/RPBE/'

# Location of results files and images
results_dir = 'results/COoxidation/'
figures_dir = 'Images/COoxidation/'

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
states += [State(state_type='gas', name='CO2', sigma=2)]

# Surfaces
states += [State(state_type='surface', path=adsdir + 'Clean')]
states += [State(state_type='surface', path=adsdir + 'OrderedPdAu-alloy/AuOnTop/Clean', name='Clean_AuTop')]
states += [State(state_type='surface', path=adsdir + 'OrderedPdAu-alloy/PdOnTop/Clean', name='Clean_PdTop')]
states += [State(state_type='surface', path=adsdir + 'PdAu-alloy/Clean', name='Clean_alloy')]

# Zip states into a dictionary
snames = [s.name for s in states]
states = dict(zip(snames, states))

# All frequencies are in a sub-folder called freqs
for s in snames:
    if states[s].path is not None:
        states[s].vibs_path = states[s].path + '/freqs_ase'

print('Done.')

# Reactions
print('Configuring reactions...')

# Site properties

# Pd
rPd = 1.63e-10  # Pd vdw radius (m)
aPd = np.pi * rPd ** 2  # Pd site area (m2)
Acat_Pd = 0.062 * 2e-4  # Total catalyst area (m2)

# AuPd
rAu = 1.66e-10  # Au vdw radius (m)
aAu = np.pi * rAu ** 2  # Au site area (m2)
Acat_Au = 0.069 * 2e-4  # Total catalyst area (m2)

# Temporary (also fix scaling!)
rAu = 1.63e-10
aAu = np.pi * rPd ** 2
Acat_Au = 0.062 * 2e-4

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
                       scaling=1.0)]  # Only allow adsorption on half the sites
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['CO'], states['Clean_alloy']],
                       products=[states['COonPd_alloy']],
                       TS=None,
                       area=aAu,
                       name='CO_ads_alloy_Pd',
                       scaling=1.0)]  # Only allow adsorption on half the sites
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean_alloy'], states['Clean_alloy']],
                       products=[states['OonAu_alloy'], states['OonAu_alloy']],
                       TS=None,
                       area=aPd,
                       name='O2_ads_alloy_Au',
                       scaling=1.0)]  # Only allow adsorption on half the sites
reactions += [Reaction(reac_type='adsorption',
                       reactants=[states['O2'], states['Clean_alloy'], states['Clean_alloy']],
                       products=[states['OonPd_alloy'], states['OonPd_alloy']],
                       TS=None,
                       area=aPd,
                       name='O2_ads_alloy_Pd',
                       scaling=1.0)]  # Only allow adsorption on half the sites

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
reactions['CO_ox_Pd'] = Reaction(reac_type='scaling', reversible=False,
                                 reactants=[states['CO.Pd111'], states['O.Pd111']],
                                 products=[states['CO2']],
                                 TS=[states['SRTS.Pd']],
                                 area=aPd,
                                 name='CO_ox_Pd')
reactions['CO_ox_AuTop'] = Reaction(reac_type='scaling', reversible=False,
                                    reactants=[states['COonAu'], states['OonAu']],
                                    products=[states['CO2']],
                                    TS=[states['SRTS.AuTop']],
                                    area=aAu,
                                    name='CO_ox_AuTop')
reactions['CO_ox_PdTop'] = Reaction(reac_type='scaling', reversible=False,
                                    reactants=[states['COonPd'], states['OonPd']],
                                    products=[states['CO2']],
                                    TS=[states['SRTS.PdTop']],
                                    area=aPd,
                                    name='CO_ox_PdTop')
reactions['CO_ox_alloy'] = Reaction(reac_type='scaling', reversible=False,
                                    reactants=[states['COonAu_alloy'], states['OonPd_alloy']],
                                    products=[states['CO2']],
                                    TS=[states['SRTS.alloy']],
                                    area=0.5*(aPd + aAu),
                                    name='CO_ox_alloy')

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
