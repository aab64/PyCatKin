from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import *
from ase.io import read
import matplotlib.pyplot as plt

font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['lines.linewidth'] = 1.5


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


# Location of mechanism file
input_file = 'input_withoutH2O.json'

sim_systems = dict()
all_states = []
all_metals = ['ag', 'au', 'aucu', 'aupd', 'co', 'cu', 'fe', 'ni', 'pd', 'pdau', 'pdcu']
metals = ['ag', 'au', 'aucu', 'aupd', 'co', 'cu', 'fe', 'ni', 'pd', 'pdcu']  # 'pdau',
temperatures = np.linspace(start=400, stop=800, num=5, endpoint=True)
for metal in metals:
    print(metal)
    # Location of outcars and frequencies
    base_data_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Metals/'
    ads_opt_dir = base_data_dir + 'energies/withoutH2O/%s/' % metal
    ads_vib_dir = base_data_dir + 'vibrations/withoutH2O/%s/' % metal
    gas_opt_dir = base_data_dir + 'energies/molecules/'
    gas_vib_dir = base_data_dir + 'vibrations/molecules/'

    # Location of results files and images
    base_out_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Metals/'
    results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
    figures_dir = base_out_dir + 'images/withoutH2O/%s/' % metal

    print('------------------------')
    print('System: DMTM without H2O')
    print('------------------------')

    sys = read_from_input_file(input_path=metal + '/' + input_file)

    for s in sys.snames:
        opt_dir = gas_opt_dir if sys.states[s].state_type == 'gas' else ads_opt_dir
        vib_dir = gas_vib_dir if sys.states[s].state_type == 'gas' else ads_vib_dir
        sys.states[s].path = opt_dir + sys.states[s].name
        sys.states[s].vibs_path = vib_dir + sys.states[s].name
        if not os.path.isdir(sys.states[s].vibs_path):
            print(sys.states[s].vibs_path)
            sys.states[s].vibs_path = base_data_dir + '../DMTM_Cu/vibrations/tmp-old-Cu/withoutH2O/' + sys.states[s].name
        read_from_alternate = None
        if sys.states[s].state_type == 'gas':
            read_from_alternate = {'get_atoms': lambda state_path=opt_dir + sys.states[s].name: load_from_contcar(
                state_path)}
        sys.states[s].read_from_alternate = read_from_alternate
        if s not in all_states:
            all_states.append(s)

    sys.params['times'][-1] = 1.0e4
    sim_systems[metal] = sys
    # sys.params['ftol'] = 1.0e-15
    # sys.params['xtol'] = 1.0e-15
    # draw_states(sim_system=sys, fig_path=figures_dir)
    # run(sim_system=sys, plot_results=True, save_results=True,
    #     fig_path=figures_dir, csv_path=results_dir)
    # run_temperatures(sim_system=sys,
    #                  steady_state_solve=True,
    #                  tof_terms=None,  # ['r5', 'r9']
    #                  temperatures=temperatures,
    #                  plot_results=False, save_results=True,
    #                  fig_path=figures_dir, csv_path=results_dir)
    # draw_energy_landscapes(sim_system=sys, etype='electronic', show_labels=True, fig_path=figures_dir)
    # draw_energy_landscapes(sim_system=sys, etype='free', show_labels=True, fig_path=figures_dir)
    # run_energy_span_temperatures(sim_system=sys,
    #                              temperatures=temperatures,
    #                              etype='free',
    #                              save_results=True,
    #                              csv_path=results_dir)
    # sys.params['temperature'] = 500
    # save_state_energies(sim_system=sys,
    #                     csv_path=results_dir)
    # save_energies(sim_system=sys,
    #               csv_path=results_dir)
    # save_energies_temperatures(sim_system=sys,
    #                            temperatures=temperatures,
    #                            csv_path=results_dir)
    # compare_energy_landscapes(sim_systems, etype='free', eunits='eV', legend_location='lower right',
    #                           show_labels=True, fig_path=figures_dir)

all_states.sort()
cmap = plt.get_cmap("Spectral", len(all_metals))
fig, ax = plt.subplots(ncols=2, figsize=(6.4, 3.2), sharey=True)
for mind, metal in enumerate(all_metals):
    if metal in metals:
        results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
        df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
        ax[0].plot(df['Temperature (K)'].values, df['r5'].values + df['r9'].values,
                   color=cmap(metals.index(metal)), label=metal)
        df = pd.read_csv(results_dir + 'energy_span_summary_full_pes.csv')
        ax[1].plot(df['Temperature (K)'].values, 2 * df['TOF (1/s)'].values,
                   color=cmap(metals.index(metal)), label=metal, linestyle='-')
ax[0].plot(temperatures, np.ones(len(temperatures)), ':', color='lightslategrey')
ax[1].plot(temperatures, np.ones(len(temperatures)), ':', color='lightslategrey')
ax[0].text(temperatures[0], 1.0, 'TOF = 1 s$^{-1}$', ha='left', va='bottom', color='lightslategrey')
ax[0].set(xlabel='Temperature (K)',
          ylabel='TOF (1/s)', yscale='log', ylim=(1e-12, 1e2),
          title='MK model')
ax[1].set(xlabel='Temperature (K)',
          yscale='log', ylim=(1e-12, 1e2),
          title='ES model')
ax[0].legend(loc='lower center', ncol=3)
fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.savefig(base_out_dir + 'images/withoutH2O/tof_vs_temperature_vs_metals.png', format='png', dpi=600)

cmap = plt.get_cmap("Spectral", len(all_metals))
fig, ax = plt.subplots(figsize=(3.2, 3.2))
for mind, metal in enumerate(all_metals):
    if metal in metals:
        results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
        df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
        ax.plot(df['Temperature (K)'].values, df['r5'].values + df['r9'].values,
                color=cmap(all_metals.index(metal)), label=metal)
        df = pd.read_csv(results_dir + 'energy_span_summary_full_pes.csv')
        ax.plot(df['Temperature (K)'].values, 2 * df['TOF (1/s)'].values,
                color=cmap(all_metals.index(metal)), label='', linestyle=':', alpha=0.7)
ax.plot(temperatures, np.ones(len(temperatures)), '--', color='lightslategrey')
ax.text(temperatures[0], 1.0, 'TOF = 1 s$^{-1}$', ha='left', va='bottom', color='lightslategrey')
ax.set(xlabel='Temperature (K)',
       ylabel='TOF (1/s)', yscale='log', ylim=(1e-24, 1e2))
ax.legend(loc='lower center', ncol=3, mode='expand')
fig.tight_layout()

fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex=True, sharey=True)
locs = [[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1], [2, 2], [2, 3]]
labs = []
for mind, metal in enumerate(all_metals):
    if metal in metals:
        results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
        df = pd.read_csv(results_dir + 'energy_span_xTDTS_full_pes.csv')
        cmap = plt.get_cmap("Spectral", len(all_states) + 2)
        for sind, s in enumerate(df.columns[1::]):
            if s not in all_states:
                all_states.append(s)
            if max(df[s].values) > 0.01:
                ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values, df[s].values,
                                                      color=cmap(all_states.index(s)), label=s)
                if s not in labs:
                    labs.append(s)
        # ax[locs[mind][0]][locs[mind][1]].legend(loc='best')
    if locs[mind][0] == 2:
        ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    if locs[mind][1] == 0:
        ax[locs[mind][0]][locs[mind][1]].set(ylabel='xTDTS')
    ax[locs[mind][0]][locs[mind][1]].set(title=metal)
for s in labs:
    ax[-1, -1].plot((temperatures[0], temperatures[-1]), (-1, -1),
                    label=s, color=cmap(all_states.index(s)))
ax[-1, -1].legend(loc='center', frameon=False, mode='expand', ncol=2, fontsize=5)
ax[-1, -1].set(ylim=(-0.05, 1.05))
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labelleft=False)
ax[-1, -1].axis('off')
fig.tight_layout()
fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.savefig(base_out_dir + 'images/withoutH2O/xTDTS_vs_temperature_vs_metals.png', format='png', dpi=600)

fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex=True, sharey=True)
locs = [[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1], [2, 2], [2, 3]]
labs = []
for mind, metal in enumerate(all_metals):
    if metal in metals:
        results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
        df = pd.read_csv(results_dir + 'energy_span_summary_full_pes.csv')
        if metal == 'pd' or metal == 'ag':
            print(metal)
            print(df)
        df = pd.read_csv(results_dir + 'energy_span_xTDI_full_pes.csv')
        cmap = plt.get_cmap("Spectral", len(all_states))
        for sind, s in enumerate(df.columns[1::]):
            if s not in all_states:
                all_states.append(s)
            if max(df[s].values) > 0.01:
                ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values, df[s].values,
                                                      color=cmap(all_states.index(s)), label=s)
                if s not in labs:
                    labs.append(s)
        # ax[locs[mind][0]][locs[mind][1]].legend(loc='best')
    if locs[mind][0] == 2:
        ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    if locs[mind][1] == 0:
        ax[locs[mind][0]][locs[mind][1]].set(ylabel='xTDI')
    ax[locs[mind][0]][locs[mind][1]].set(title=metal)
for s in labs:
    ax[-1, -1].plot((temperatures[0], temperatures[-1]), (-1, -1),
                    label=s, color=cmap(all_states.index(s)))
ax[-1, -1].legend(loc='center', frameon=False, ncol=2, fontsize=5, mode='expand')
ax[-1, -1].set(ylim=(-0.05, 1.05))
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labelleft=False)
ax[-1, -1].axis('off')
fig.tight_layout()
fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.savefig(base_out_dir + 'images/withoutH2O/xTDI_vs_temperature_vs_metals.png', format='png', dpi=600)

fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex=True, sharey=True)
locs = [[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1], [2, 2], [2, 3]]
labs = []
for mind, metal in enumerate(all_metals):
    if metal in metals:
        results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
        df = pd.read_csv(results_dir + 'coverages_vs_temperature.csv')
        cmap = plt.get_cmap("Spectral", len(all_states) + 2)
        for sind, s in enumerate(df.columns[1::]):
            if s not in all_states:
                all_states.append(s)
            if max(df[s].values) > 0.01:
                ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values, df[s].values,
                                                      color=cmap(all_states.index(s)), label=s)
                if s not in labs:
                    labs.append(s)
        # ax[locs[mind][0]][locs[mind][1]].legend(loc='best')
    if locs[mind][0] == 2:
        ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    if locs[mind][1] == 0:
        ax[locs[mind][0]][locs[mind][1]].set(ylabel='Coverage')
    ax[locs[mind][0]][locs[mind][1]].set(title=metal)
for s in labs:
    ax[-1, -1].plot((temperatures[0], temperatures[-1]), (-1, -1),
                    label=s, color=cmap(all_states.index(s)))
ax[-1, -1].legend(loc='center', frameon=False, mode='expand', ncol=2, fontsize=5)
ax[-1, -1].set(ylim=(-0.05, 1.05))
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labelleft=False)
ax[-1, -1].axis('off')
fig.tight_layout()
fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.savefig(base_out_dir + 'images/withoutH2O/coverage_vs_temperature_vs_metals.png', format='png', dpi=600)

##

cmap = plt.get_cmap("Spectral", len(all_metals))
fig, ax = plt.subplots(ncols=2, figsize=(6.4, 3.2), sharey=True)
evals = []
eavals = []
gvals = []
gavals = []
svals = []
gsvals = []
rvals = []
ind = 1
for mind, metal in enumerate(all_metals):
    if metal in metals:
        results_dir = base_out_dir + 'results/withoutH2O/%s/' % metal
        df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
        rvals += [np.log10(df['r5'][ind] + df['r9'][ind])]

        df = pd.read_csv(results_dir + 'reaction_energies_and_barriers_r5.csv')
        evals += [df['dEr (J/mol)'][ind] * 1.0e-3 / eVtokJ]
        gvals += [df['dGr (J/mol)'][ind] * 1.0e-3 / eVtokJ]
        df = pd.read_csv(results_dir + 'reaction_energies_and_barriers_r4.csv')
        eavals += [df['dEa (J/mol)'][ind] * 1.0e-3 / eVtokJ]
        gavals += [df['dGa (J/mol)'][ind] * 1.0e-3 / eVtokJ]
        df = pd.read_csv(results_dir + 'state_energies_500.0K_1.0bar.csv')
        svals += [df['Electronic (eV)'][df[df['State'] == 'sO'].index.item()] +
                  df['Electronic (eV)'][df[df['State'] == 'ch3oh'].index.item()] -
                  df['Electronic (eV)'][df[df['State'] == '2s'].index.item()] -
                  df['Electronic (eV)'][df[df['State'] == 'o2'].index.item()] -
                  df['Electronic (eV)'][df[df['State'] == 'ch4'].index.item()]]
        gsvals += [df['Free (eV)'][df[df['State'] == 'sO'].index.item()] +
                   df['Free (eV)'][df[df['State'] == 'ch3oh'].index.item()] -
                   df['Free (eV)'][df[df['State'] == '2s'].index.item()] -
                   df['Free (eV)'][df[df['State'] == 'o2'].index.item()] -
                   df['Free (eV)'][df[df['State'] == 'ch4'].index.item()]]

ax[0].plot(evals, rvals, 'o', color='k')
ax[1].plot(gvals, rvals, 'o', color='k')

ax[0].set(xlabel='Reaction energy, r5 (eV)',
          ylabel='log(TOF) (1/s)',
          title='Electronic')
ax[1].set(xlabel='Reaction energy, r5 (eV)',
          title='Free')
fig.tight_layout()
fig.subplots_adjust(wspace=0)
# plt.savefig(base_out_dir + 'images/withoutH2O/tof_vs_temperature_vs_metals.png', format='png', dpi=600)

fig, ax = plt.subplots(ncols=2, figsize=(6.4, 3.2), sharey=True)
ax[0].plot(evals, eavals, 'o', color='k')
ax[1].plot(gvals, gavals, 'o', color='k')

ax[0].set(xlabel='Reaction energy, r5 (eV)',
          ylabel='Reaction barrier, r4 (eV)',
          title='Electronic')
ax[1].set(xlabel='Reaction energy, r5 (eV)',
          title='Free')
fig.tight_layout()
fig.subplots_adjust(wspace=0)

fig, ax = plt.subplots(ncols=1, figsize=(3.2, 3.2))
ax.plot(svals, rvals, 'o', color='k')
ax.set(xlabel='dE, sO (eV)',
       ylabel='TOF (1/s)')
fig.tight_layout()
fig.subplots_adjust(wspace=0)

fig, ax = plt.subplots(ncols=1, figsize=(3.2, 3.2))
ax.plot(gsvals, rvals, 'o', color='k')
ax.set(xlabel='dG, sO (eV)',
       ylabel='TOF (1/s)')
fig.tight_layout()
fig.subplots_adjust(wspace=0)

fig, ax = plt.subplots(ncols=1, figsize=(3.2, 3.2))
ax.plot(svals, eavals, 'o', color='k')
ax.set(xlabel='dE, sO (eV)',
       ylabel='dGa, r4 (eV)')
fig.tight_layout()
fig.subplots_adjust(wspace=0)
