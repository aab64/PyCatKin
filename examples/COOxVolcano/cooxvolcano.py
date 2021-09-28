from pycatkin.functions.load_input import read_from_input_file
import numpy as np
import matplotlib.pyplot as plt
import os

# Load input file
sim_system = read_from_input_file()

# Define a range of binding energies
be = np.linspace(start=-2.5, stop=0.5, num=10, endpoint=True)

# Standard entropies (taken from Atkins, in eV/K)
SCOg = 2.0487e-3
SO2g = 2.1261e-3

# Note the temperature
T = sim_system.params['temperature']

# Loop over energies and compute the activity
activity = np.zeros((len(be), len(be)))

for iCO, ECO in enumerate(be):
    for iO, EO in enumerate(be):

        print('* Binding energy of CO is %1.2f, binding energy of O is %1.2f' % (ECO, EO))

        # (a) Set CO adsorption energy and entropy
        sim_system.reactions['CO_ads'].dErxn_user = ECO
        sim_system.reactions['CO_ads'].dGrxn_user = ECO + SCOg * T

        # (b) Set O adsorption energy and entropy
        sim_system.reactions['2O_ads'].dErxn_user = 2.0 * EO
        sim_system.reactions['2O_ads'].dGrxn_user = 2.0 * EO + SO2g * T

        # (c) Set O2 adsorption energy and entropy
        EO2 = sim_system.states['sO2'].get_potential_energy()
        sim_system.reactions['O2_ads'].dErxn_user = EO2
        sim_system.reactions['O2_ads'].dGrxn_user = EO2 + SO2g * T

        # (d) Set CO oxidation barrier
        ETS_CO_ox = sim_system.states['SRTS_ox'].get_potential_energy()
        sim_system.reactions['CO_ox'].dEa_fwd_user = np.max((ETS_CO_ox - (ECO + EO), 0.0))

        # (e) Set O2 dissociation barrier
        ETS_O2_2O = sim_system.states['SRTS_O2'].get_potential_energy()
        sim_system.reactions['O2_2O'].dEa_fwd_user = np.max((ETS_O2_2O - EO2, 0.0))

        # Now compute and save the activity
        activity[iCO, iO] = sim_system.activity(tof_terms=['CO_ox'])


if not os.path.isdir('figures'):
    os.mkdir('figures')

fig, ax = plt.subplots(figsize=(4, 3))
CS = ax.contourf(be, be, activity, levels=25, cmap=plt.get_cmap("RdYlBu_r"))
fig.colorbar(CS).ax.set_ylabel('Activity (eV)')
ax.set(xlabel=r'$E_{\mathsf{O}}$ (eV)', ylabel=r'$E_{\mathsf{CO}}$ (eV)')
fig.tight_layout()
fig.savefig('figures/activity.png', format='png', dpi=600)
