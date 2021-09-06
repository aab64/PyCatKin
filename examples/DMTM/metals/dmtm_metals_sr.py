from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import *
import copy

font = {'family': 'sans-serif', 'weight': 'bold', 'size': 8}
mpl.rc('font', **font)
mpl.rcParams['font.sans-serif'] = "Ariel"
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['xtick.bottom'] = False
mpl.rcParams['xtick.top'] = False
mpl.rcParams['ytick.left'] = False
mpl.rcParams['ytick.right'] = False
mpl.rcParams['xtick.color'] = (87 / 255, 87 / 255, 87 / 255)
mpl.rcParams['ytick.color'] = (87 / 255, 87 / 255, 87 / 255)
plt.rcParams.update({'figure.max_open_warning': 40})

base_out_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Metals/'
all_metals = ['ag', 'au', 'aucu', 'aupd', 'co', 'cu', 'fe', 'ni', 'pd', 'pdcu', 'zn']
meclrs = [(170, 170, 164),
          (230, 198, 56),
          (255, 133, 43),
          (28, 151, 66),
          (126, 154, 241),
          (251, 73, 53),
          (167, 96, 64),
          (138, 142, 142),
          (52, 95, 209),
          (161, 74, 181),
          (190, 215, 212)]
meclrs = [np.divide(m, 255) for m in meclrs]
gclr = (87 / 255, 87 / 255, 87 / 255)

minima = [
    ["2s", "o2", "ch4", "ch4"],
    ["ts1", "o2", "ch4", "ch4"],
    ["s-pair", "o2", "ch4", "ch4"],
    ["sO2s", "ch4", "ch4"],
    ["ts2", "ch4", "ch4"],
    ["sOOs", "ch4", "ch4"],
    ["s2Och4", "ch4"],
    ["rad1", "ch4"],
    ["sOsCH3OH", "ch4"],
    ["sO", "ch4", "ch3oh"],
    ["sOch4", "ch3oh"],
    ["rad2", "ch3oh"],
    ["sOHsCH3", "ch3oh"],
    ["ts5", "ch3oh"],
    ["sCH3OH", "ch3oh"],
    ["s", "ch3oh", "ch3oh"],
    ["ts6", "ch3oh", "ch3oh"],
    ["s-pair.1", "ch3oh", "ch3oh"]]

sim_systems = dict()

fig1, ax1 = plt.subplots(figsize=(3.2, 3.2))
fig2, ax2 = plt.subplots(ncols=3, nrows=2, figsize=(3.2, 3.2), sharey='row', sharex='col')
fig3, ax3 = plt.subplots(ncols=3, figsize=(6.4, 3.2), sharey='row')

for Ti, temperature in enumerate([500, 650, 800]):
    for study in ['dry', 'wet']:

        sim_systems[study] = dict()

        results_dir = base_out_dir + 'results/%s/' % study
        figures_dir = base_out_dir + 'images/%s/' % study
        input_file = 'input_%s_sr_subset.json' % study
        sys = read_from_input_file(input_path=input_file)
        sys.params['temperature'] = temperature

        gas_entropies = dict()
        ges = dict()
        for gas in ['o2_mk', 'ch4_mk', 'ch3oh_mk']:
            sys.states[gas].get_free_energy(T=sys.params['temperature'],
                                            p=sys.params['pressure'],
                                            verbose=sys.params['verbose'])
            gas_entropies[gas] = sys.states[gas].Gtran + sys.states[gas].Grota  # + sys.states[gas].Gvibr
            ges[gas] = sys.states[gas].Gtran + sys.states[gas].Grota + sys.states[gas].Gvibr

        for m in minima:
            modifier = sum([gas_entropies[g + '_mk'] for g in m[1::]])
            modifier -= sum([gas_entropies[g + '_mk'] for g in minima[0][1::]])
            if 'Och4' in m[0]:
                modifier += (gas_entropies['ch4_mk'] * 0.67)
            sys.states[m[0]].set_energy_modifier(modifier=modifier)

        bsOs = np.zeros(len(all_metals))
        # bsOs = np.linspace(start=-6, stop=0, num=50, endpoint=True)

        tofs = np.zeros(len(bsOs))
        acts = np.zeros(len(bsOs))
        # for bind, metal in enumerate(bsOs):
        #     bsO = metal
        for bind, metal in enumerate(all_metals):

            sim_systems[study][metal] = dict()
            sim_systems[study][metal]['DFT'] = read_from_input_file(input_path='%s/input_%s.json' % (metal, study))
            sim_systems[study][metal]['DFT'].params['temperature'] = temperature

            ads_opt_dir = metal + '/data' if study == 'dry' else '/wetdata' + '/energies/%s/%s/' % (study, metal)
            gas_opt_dir = metal + '/data/energies/%s/%s/' % (study, metal)
            ads_vib_dir = metal + '/data' if study == 'dry' else '/wetdata' + '/vibrations/%s/%s/' % (study, metal)
            gas_vib_dir = metal + '/data/vibrations/%s/%s/' % (study, metal)

            for s in sim_systems[study][metal]['DFT'].snames:
                opt_dir = gas_opt_dir if sim_systems[study][metal]['DFT'].states[s].state_type == 'gas' else ads_opt_dir
                vib_dir = gas_vib_dir if sim_systems[study][metal]['DFT'].states[s].state_type == 'gas' else ads_vib_dir
                sim_systems[study][metal]['DFT'].states[s].path = opt_dir + s + '_energy.dat'
                sim_systems[study][metal]['DFT'].states[s].vibs_path = vib_dir + s + '_frequencies.dat'

            for s in ['o2', 'ch4', 'ch3oh']:
                sim_systems[study][metal]['DFT'].states[s].sigma = sys.states[s].sigma
                sim_systems[study][metal]['DFT'].states[s].mass = sys.states[s].mass
                sim_systems[study][metal]['DFT'].states[s].inertia = sys.states[s].inertia

            for s in ['sO', 'o2', 'ch4', 'ch3oh', '2s']:
                if study == 'wet' and s in ['sO', '2s']:
                    s += '-h2o'
                sim_systems[study][metal]['DFT'].states[s].get_free_energy(T=sys.params['temperature'],
                                                                           p=sys.params['pressure'])
            bsO = (sim_systems[study][metal]['DFT'].states['sO' if study == 'dry' else 'sO-h2o'].Gelec +
                   sim_systems[study][metal]['DFT'].states['ch4'].Gelec +
                   sim_systems[study][metal]['DFT'].states['ch3oh'].Gelec) - (
                sim_systems[study][metal]['DFT'].states['2s' if study == 'dry' else '2s-h2o'].Gelec +
                sim_systems[study][metal]['DFT'].states['o2'].Gelec +
                sim_systems[study][metal]['DFT'].states['ch4'].Gelec +
                sim_systems[study][metal]['DFT'].states['ch4'].Gelec)

            print('%s: %1.2f' % (metal, bsO))

            sim_systems[study][metal]['scaling (gas T+R)'] = copy.deepcopy(sys)
            sim_systems[study][metal]['scaling (gas T+R)'].states['sO'].Gelec = bsO
            sim_systems[study][metal]['scaling (gas T+R)'].states['sO'].freq = np.array((0, 0))
            sim_systems[study][metal]['scaling (gas T+R)'].reactions['rsO'].dErxn_user = bsO

            sim_systems[study][metal]['scaling (gas V+T+R)'] = copy.deepcopy(sys)
            sim_systems[study][metal]['scaling (gas V+T+R)'].states['sO'].Gelec = bsO
            sim_systems[study][metal]['scaling (gas V+T+R)'].states['sO'].freq = np.array((0, 0))
            sim_systems[study][metal]['scaling (gas V+T+R)'].reactions['rsO'].dErxn_user = bsO
            for m in minima:
                modifier = sum([ges[g + '_mk'] for g in m[1::]])
                modifier -= sum([ges[g + '_mk'] for g in minima[0][1::]])
                if 'Och4' in m[0]:
                    modifier += (ges['ch4_mk'] * 0.67)
                sim_systems[study][metal]['scaling (gas V+T+R)'].states[m[0]].set_energy_modifier(modifier=modifier)

            def cmap(ind):
                if not ind:
                    return meclrs[bind]
                elif ind == 1:
                    return gclr
                else:
                    return 'sienna'

            if metal != 'zn':
                sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima[0:-2]
                sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels[0:-2]
            else:
                sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima[0:-1]
                sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels[0:-1]
            # compare_energy_landscapes(sim_systems=sim_systems[study][metal],
            #                           etype='electronic',
            #                           legend_location='best',
            #                           show_labels=True,
            #                           fig_path=figures_dir + '%s/' % metal,
            #                           cmap=cmap)
            compare_energy_landscapes(sim_systems=sim_systems[study][metal],
                                      etype='free',
                                      legend_location='best',
                                      show_labels=True,
                                      fig_path=figures_dir + '%s/' % metal,
                                      cmap=cmap)

            run(sim_system=sim_systems[study][metal]['scaling (gas T+R)'],
                plot_results=False,
                save_results=False,
                steady_state_solve=True,
                fig_path=figures_dir + 'sr/%s/' % metal,
                csv_path=results_dir + 'sr/%s/' % metal)
            sim_systems[study][metal]['scaling (gas T+R)'].reaction_terms(y=sim_systems[study][metal]['scaling (gas T+R)'].full_steady)
            i5 = list(sim_systems[study][metal]['scaling (gas T+R)'].reactions.keys()).index('r5_rdr')
            i9 = list(sim_systems[study][metal]['scaling (gas T+R)'].reactions.keys()).index('r9_rdr')
            tofs[bind] = (sim_systems[study][metal]['scaling (gas T+R)'].rates[i5, 0] -
                          sim_systems[study][metal]['scaling (gas T+R)'].rates[i5, 1]) + (
                    sim_systems[study][metal]['scaling (gas T+R)'].rates[i9, 0] -
                    sim_systems[study][metal]['scaling (gas T+R)'].rates[i9, 1])
            acts[bind] = (np.log((h * tofs[bind]) / (kB * sys.params['temperature'])) *
                          (R * sys.params['temperature'])) * 1.0e-3 / eVtokJ
            bsOs[bind] = bsO

            ###

        inds = np.argsort(bsOs)
        ax1.plot(bsOs[inds], acts[inds],
                 marker='o' if study == 'wet' else 'v',
                 color='k' if study == 'wet' else 'k',
                 linestyle=':' if study == 'wet' else '-',
                 label=study, markevery=[])
        if study == 'wet':
            ax3[Ti].plot(bsOs[inds], tofs[inds],
                         marker='o' if study == 'wet' else 'o',
                         color='k' if study == 'wet' else 'k',
                         linestyle='-' if study == 'wet' else '-',
                         label=study, markevery=[])
        ax2[0 if study == 'wet' else 1][Ti].plot(bsOs[inds], tofs[inds],
                                                 marker='o' if study == 'wet' else 'v',
                                                 color='k' if study == 'wet' else 'k',
                                                 linestyle=':' if study == 'wet' else '-',
                                                 label=study, markevery=[])

    mEsOs = []
    for mind, metal in enumerate(all_metals):
        for sind, study in enumerate(['dry', 'wet']):
            results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
            descriptor = 'sO' if study == 'dry' else 'sO-h2o'
            df = pd.read_csv(results_dir + 'full_pes_energy_landscape_%1.1fK_1.0bar.csv' % sys.params['temperature'])
            EsO = df['Electronic (eV)'][list(df['State']).index('sO')]
            mEsOs.append(EsO)
            df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
            Tind = list(df['Temperature (K)']).index(sys.params['temperature'])
            ax1.plot(EsO,
                     (np.log((h * (df['r5'][Tind] + df['r9'][Tind])) / (kB * sys.params['temperature'])) *
                      (R * sys.params['temperature'])) * 1.0e-3 / eVtokJ,
                     color=meclrs[all_metals.index(metal)],
                     marker='o' if study == 'wet' else 'v')
            if study == 'wet':
                ax3[Ti].plot(EsO, df['r5'][Tind] + df['r9'][Tind],
                             color=meclrs[all_metals.index(metal)],
                             marker='o' if study == 'wet' else 'o')
            ax2[0 if study == 'wet' else 1][Ti].plot(EsO, df['r5'][Tind] + df['r9'][Tind],
                                                     color=meclrs[all_metals.index(metal)],
                                                     marker='o' if study == 'wet' else 'v')

ax1.set(xlabel='E$_{*O}$ (eV)', ylabel='Activity')
ax2[0][0].set(ylabel='TOF (1/s)', yscale='log', ylim=(2e-16, 2e2))
ax2[1][0].set(xlabel='E$_{*O}$ (eV)', ylabel='TOF (1/s)', yscale='log', xticks=(-4, -2, 0),
              ylim=(2e-16, 2e2),
              xlim=(-6.5, 0.5))
ax2[1][1].set(xlabel='E$_{*O}$ (eV)', xticks=(-4, -2, 0), xlim=(-6.5, 0.5))
ax2[1][2].set(xlabel='E$_{*O}$ (eV)', xticks=(-4, -2, 0), xlim=(-6.5, 0.5))

ax3[0].set(xlabel='E$_{*O}$ (eV)', ylabel='TOF (1/s)', yscale='log', ylim=(1e-14, 1e2))
ax1.legend(loc='best', frameon=False)
# ax2[0][0].legend(loc='best', frameon=False)
ax3[0].legend(loc='best', frameon=False)

ax2[0][0].set(title='500 K')
ax2[0][1].set(title='650 K')
ax2[0][2].set(title='800 K')

fig1.tight_layout()
fig2.tight_layout()
fig2.subplots_adjust(wspace=0, hspace=0)
# fig1.savefig(base_out_dir + 'images/scaling_activity.png', format='png', dpi=600)
fig2.savefig(base_out_dir + 'images/scaling_tof_novibs.png', format='png', dpi=600)
fig2.savefig(base_out_dir + 'images/scaling_tof_novibs.eps', format='eps', dpi=600)
