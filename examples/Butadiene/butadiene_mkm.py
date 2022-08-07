from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import get_tof_for_given_reactions
import matplotlib.pyplot as plt
import copy
import numpy as np

# Load DFT mechanisms input file
dft_system = read_from_input_file()

# Load microkinetic model input file
mkm_system = read_from_input_file(input_path='input_mkm.json',
                                  base_system=dft_system)

# Define pathways to consider
adsorption_reactions = ['9D-9C', 'ethanol-1A', '8A-8C', 'H2O-9B', 'acetaldehyde-10B', 'crotonaldehyde-2N']
p123 = ['1A-1C', '2A-2C', '2F-2H', '2J-2L', '2L-2N', '3A-3C', '3D-3F', '3F-3G'] + adsorption_reactions
p124 = ['1A-1C', '2A-2C', '2F-2H', '4A-4C', '4D-4Ca', '4D-4F', '4F-4H', '4I-4K'] + adsorption_reactions
p156 = ['1A-1C', '5A-5C', '6A-6C', '6C-6E', '6E-6G', '6G-6H'] + adsorption_reactions
souza_butanol = ['12A-12C', '12C-12E', '12E-12G', 'butanol-12G']
jonas_butanol = ['3Ci-3Ciii', '3Civ-3Cvi', 'butanol-3Cvi']
p123_p124_p156 = list(set(p123 + p124 + p156))
with_souza_byproducts = p123_p124_p156 + souza_butanol + ['7Ei-7Eiii', 'ethylacetate-7Eiii']
with_jonas_byproducts = p123_p124_p156 + jonas_butanol + ['7Ei-7Eiii', 'ethylacetate-7Eiii']

# Define temperature range (K)
Ts = list(np.linspace(start=523, stop=923, num=17, endpoint=True))

# Create a dictionary to store the results
mkm_results = dict()

# For each set of reactions, solve the ODEs and get the TOF
for test_case, pathways in zip(['p123_p124_p156', 'p123', 'p124', 'p156',
                                'with_souza_byproducts', 'with_jonas_byproducts'],
                               [p123_p124_p156, p123, p124, p156,
                                with_souza_byproducts, with_jonas_byproducts]):
    mkm_results[test_case] = dict()

    # Make a copy of the system and delete reactions not considered
    sim_system = copy.deepcopy(mkm_system)
    reactions_to_discard = []
    for reaction in sim_system.reactions.keys():
        if reaction not in pathways:
            reactions_to_discard.append(reaction)
    for reaction in reactions_to_discard:
        del sim_system.reactions[reaction]
    sim_system.names_to_indices()

    # Loop over temperatures and solve
    for T in Ts:
        sim_system.params['temperature'] = T

        # Solve the ODEs
        sim_system.solve_odes()

        # Get the TOFs for butadiene and byproducts
        bd_tof = get_tof_for_given_reactions(sim_system=sim_system,
                                             tof_terms=['3F-3G', '4I-4K', '6G-6H'])
        bu_tof = get_tof_for_given_reactions(sim_system=sim_system,
                                             tof_terms=['butanol-3Cvi', 'butanol-12G']) * -1.0
        ea_tof = get_tof_for_given_reactions(sim_system=sim_system,
                                             tof_terms=['ethylacetate-7Eiii']) * -1.0
        eth_tof = get_tof_for_given_reactions(sim_system=sim_system,
                                              tof_terms=['ethanol-1A', 'ethanol-1A_Cu', 'ethanol-1A_Zn'])

        # Store the final rates
        sim_system.reaction_terms(y=sim_system.solution[-1])
        rates = sim_system.rates[:, 0] - sim_system.rates[:, 1]

        # Compute the steady state solution and check it is a steady state (small fval)
        sim_system.find_steady(store_steady=True)
        final_fval = sim_system.species_odes(y=sim_system.full_steady)[sim_system.adsorbate_indices]

        # Save the results
        mkm_results[test_case][T] = {'tof': bd_tof,
                                     'rates': rates,
                                     'cover': sim_system.full_steady[sim_system.adsorbate_indices],
                                     'odesol': np.concatenate((sim_system.times.reshape(-1, 1),
                                                               sim_system.solution[:, sim_system.adsorbate_indices]),
                                                              axis=1),
                                     'final_fval': pow(np.linalg.norm(final_fval), 2),
                                     'ea_tof': ea_tof,
                                     'bu_tof': bu_tof,
                                     'eth_tof': eth_tof
                                     }

# Plot the butadiene turnover frequency for each pathway
fig, ax = plt.subplots(figsize=(3.2, 3.2))
for ti, test_case in enumerate(['p123_p124_p156', 'p123', 'p124', 'p156']):
    cl = 'k' if ti == 0 else 'purple' if ti == 1 else 'dodgerblue' if ti == 2 else 'orange'
    ax.plot(Ts, [mkm_results[test_case][T]['tof'] for T in Ts],
            label=test_case,
            color=cl)
ax.set(xlabel='Temperature (K)',
       ylabel='TOF (1/s)',
       xlim=(523, 923),
       ylim=(1e-8, 1e2),
       yscale='log')
ax.legend()
ax.grid()
fig.tight_layout()
fig.show()

# Plot the butadiene and byproduct turnover frequencies
fig, ax = plt.subplots(figsize=(3.2, 3.2))
for ki, k in enumerate(['p123_p124_p156', 'with_souza_byproducts', 'with_jonas_byproducts']):
    ls = '-' if ki == 0 else '--' if ki == 1 else ':'
    ax.plot(Ts, [mkm_results[k][T]['tof'] for T in Ts],
            color='black',
            label='base case' if ki == 0 else 'Souza steps' if ki == 1 else 'Jonas steps',
            linestyle=ls,
            alpha=0.6)
    ax.plot(Ts, [mkm_results[k][T]['bu_tof'] for T in Ts],
            color='mediumvioletred',
            linestyle=ls,
            alpha=0.6)
ax.text(865, 1e-2, 'BD',
        color='k', ha='right', va='bottom', weight='bold')
ax.text(865, 2e-7, 'BuOH',
        color='mediumvioletred', ha='right', va='bottom', weight='bold')
ax.set(xlabel='Temperature (K)',
       ylabel='TOF (1/s)',
       xlim=(523, 923),
       ylim=(1e-12, 1e0),
       yscale='log',
       yticks=np.logspace(start=-12, stop=0, endpoint=True, num=9))
ax.legend()
fig.tight_layout()
fig.show()
