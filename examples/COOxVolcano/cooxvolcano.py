from pycatkin.functions.load_input import *
from pycatkin.functions.presets import *

# Load input file
sim_system = read_from_input_file()

npts = 10
be = np.linspace(start=-2.5, stop=0.5, num=npts, endpoint=True)
act = np.zeros((npts, npts))

# Standard entropies (taken from Atkins, in kJ/molK)
SCOg = 0.19767 / eVtokJ
SO2g = 0.205138 / eVtokJ

for iCO, bCO in enumerate(be):
    for iO, bO in enumerate(be):

        # Set CO adsorption energy and entropy
        sim_system.reactions['CO_ads'].dErxn_user = bCO
        sim_system.reactions['CO_ads'].dGrxn_user = bCO + SCOg * sim_system.params['temperature']

        # Set O adsorption energy and entropy
        sim_system.reactions['2O_ads'].dErxn_user = 2.0 * bO
        sim_system.reactions['2O_ads'].dGrxn_user = 2.0 * bO + SO2g * sim_system.params['temperature']

        # Add adsorption entropy change for gases
        sim_system.states['sO2'].set_energy_modifier(modifier=SO2g * sim_system.params['temperature'])
        sim_system.states['SRTS_O2'].set_energy_modifier(modifier=SO2g * sim_system.params['temperature'])
        sim_system.states['SRTS_ox'].set_energy_modifier(modifier=(0.5 * SO2g * sim_system.params['temperature'] +
                                                                   SCOg * sim_system.params['temperature']))

        # Set O2 adsorption free energy
        sim_system.reactions['O2_ads'].dGrxn_user = sim_system.states['sO2'].get_free_energy(
            T=sim_system.params['temperature'],
            p=sim_system.params['pressure'])

        # Set CO oxidation free energy barrier
        sim_system.reactions['CO_ox'].dGa_fwd_user = sim_system.states['SRTS_ox'].get_free_energy(
            T=sim_system.params['temperature'],
            p=sim_system.params['pressure']) - (sim_system.reactions['CO_ads'].dGrxn_user +
                                                0.5 * sim_system.reactions['2O_ads'].dGrxn_user)

        # Set O2 dissociation free energy barrier
        sim_system.reactions['O2_2O'].dGa_fwd_user = sim_system.states['SRTS_O2'].get_free_energy(
            T=sim_system.params['temperature'],
            p=sim_system.params['pressure']) - sim_system.reactions['O2_ads'].dGrxn_user

        # Now compute and save the activity
        act[iCO, iO] = sim_system.activity(tof_terms=['CO_ox'])

fig, ax = plt.subplots(figsize=(4, 3))
CS = ax.contourf(be, be, act, levels=25, cmap=plt.get_cmap("RdYlBu_r"))
fig.colorbar(CS).ax.set_ylabel('Activity (eV)')
ax.set(xlabel=r'$E_{\mathsf{O}}$ (eV)', ylabel=r'$E_{\mathsf{CO}}$ (eV)')
fig.tight_layout()
fig.savefig('figures/activity.png', format='png', dpi=600)
