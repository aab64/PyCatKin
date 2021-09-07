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
fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(3.2, 3.2), sharey='row', sharex='col')

for Ti, temperature in enumerate([500, 650, 800]):
    for study in ['dry', 'wet']:
        sim_systems[study] = dict()
        results_dir = base_out_dir + 'results/%s/' % study
        figures_dir = base_out_dir + 'images/%s/' % study
        input_file = 'input_%s_sr_subset.json' % study
        sys = read_from_input_file(input_path=input_file)
        sys.params['temperature'] = temperature

        gas_entropies = dict()
        for gas in ['o2_mk', 'ch4_mk', 'ch3oh_mk']:
            sys.states[gas].get_free_energy(T=sys.params['temperature'],
                                            p=sys.params['pressure'],
                                            verbose=sys.params['verbose'])
            gas_entropies[gas] = sys.states[gas].Gtran + sys.states[gas].Grota + sys.states[gas].Gvibr

        for m in minima:
            modifier = sum([gas_entropies[g + '_mk'] for g in m[1::]])
            modifier -= sum([gas_entropies[g + '_mk'] for g in minima[0][1::]])
            if 'Och4' in m[0]:
                modifier += ((gas_entropies['ch4_mk'] - sys.states['ch4_mk'].Gvibr) * 0.67)
            sys.states[m[0]].set_energy_modifier(modifier=modifier)

        # bsOs = np.zeros(len(all_metals))
        bsOs = np.linspace(start=-5, stop=0, num=50, endpoint=True)
        tofs = np.zeros(len(bsOs))

        for bind, metal in enumerate(bsOs):
            bsO = metal
        # for bind, metal in enumerate(all_metals):

            sim_systems[study][metal] = dict()
            # sim_systems[study][metal]['DFT'] = read_from_input_file(input_path='%s/input_%s.json' % (metal, study))
            # sim_systems[study][metal]['DFT'].params['temperature'] = temperature

            # for s in ['sO', 'o2', 'ch4', 'ch3oh', '2s']:
            #     if study == 'wet' and s in ['sO', '2s']:
            #         s += '-h2o'
            #     sim_systems[study][metal]['DFT'].states[s].get_free_energy(T=sys.params['temperature'],
            #                                                                p=sys.params['pressure'])
            # bsO = (sim_systems[study][metal]['DFT'].states['sO' if study == 'dry' else 'sO-h2o'].Gelec +
            #        sim_systems[study][metal]['DFT'].states['ch4'].Gelec +
            #        sim_systems[study][metal]['DFT'].states['ch3oh'].Gelec) - (
            #     sim_systems[study][metal]['DFT'].states['2s' if study == 'dry' else '2s-h2o'].Gelec +
            #     sim_systems[study][metal]['DFT'].states['o2'].Gelec +
            #     sim_systems[study][metal]['DFT'].states['ch4'].Gelec +
            #     sim_systems[study][metal]['DFT'].states['ch4'].Gelec)
            print('%s: %1.2f' % (metal, bsO))

            sim_systems[study][metal]['scaling'] = copy.deepcopy(sys)
            sim_systems[study][metal]['scaling'].states['sO'].Gelec = bsO
            # sim_systems[study][metal]['scaling'].states['sO'].freq = np.array((0, 0))
            sim_systems[study][metal]['scaling'].reactions['rsO'].dErxn_user = bsO
            #
            # def cmap(ind):
            #     if not ind:
            #         return meclrs[bind]
            #     else:
            #         return gclr
            #
            # if metal != 'zn':
            #     sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima[0:-2]
            #     sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels[0:-2]
            # else:
            #     sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima[0:-1]
            #     sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels[0:-1]
            # compare_energy_landscapes(sim_systems=sim_systems[study][metal],
            #                           etype='electronic',
            #                           legend_location='best',
            #                           show_labels=True,
            #                           fig_path=figures_dir + '%s/' % metal,
            #                           cmap=cmap)
            # compare_energy_landscapes(sim_systems=sim_systems[study][metal],
            #                           etype='free',
            #                           legend_location='best',
            #                           show_labels=True,
            #                           fig_path=figures_dir + '%s/' % metal,
            #                           cmap=cmap)

            run(sim_system=sim_systems[study][metal]['scaling'],
                plot_results=False,
                save_results=False,
                steady_state_solve=True,
                fig_path=figures_dir + 'sr/%s/' % metal,
                csv_path=results_dir + 'sr/%s/' % metal)
            
            sim_systems[study][metal]['scaling'].reaction_terms(y=sim_systems[study][metal]['scaling'].full_steady)
            i5 = list(sim_systems[study][metal]['scaling'].reactions.keys()).index('r5_rdr')
            i9 = list(sim_systems[study][metal]['scaling'].reactions.keys()).index('r9_rdr')
            tofs[bind] = (sim_systems[study][metal]['scaling'].rates[i5, 0] -
                          sim_systems[study][metal]['scaling'].rates[i5, 1]) + (
                    sim_systems[study][metal]['scaling'].rates[i9, 0] -
                    sim_systems[study][metal]['scaling'].rates[i9, 1])
            bsOs[bind] = bsO
        inds = np.argsort(bsOs)

        ax[0 if study == 'wet' else 1][Ti].plot(bsOs[inds], tofs[inds],
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
            ax[0 if study == 'wet' else 1][Ti].plot(EsO, df['r5'][Tind] + df['r9'][Tind],
                                                    color=meclrs[all_metals.index(metal)],
                                                    marker='o' if study == 'wet' else 'v')

ax[0][0].set(ylabel='TOF (1/s)', yscale='log', ylim=(1e-16, 1e3))
ax[1][0].set(ylabel='TOF (1/s)', yscale='log', ylim=(1e-16, 1e3),
             xlabel='E$_{*O}$ (eV)', xticks=(-6, -4, -2), xlim=(-7, 0))
ax[1][1].set(xlabel='E$_{*O}$ (eV)', xticks=(-6, -4, -2), xlim=(-7, 0))
ax[1][2].set(xlabel='E$_{*O}$ (eV)', xticks=(-6, -4, -2), xlim=(-7, 0))
# ax[0][0].legend(loc='best', frameon=False)
ax[0][0].set(title='500 K')
ax[0][1].set(title='650 K')
ax[0][2].set(title='800 K')

fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
fig.savefig(base_out_dir + 'images/scaling_tof.png', format='png', dpi=600)
fig.savefig(base_out_dir + 'images/scaling_tof.eps', format='eps', dpi=600)
