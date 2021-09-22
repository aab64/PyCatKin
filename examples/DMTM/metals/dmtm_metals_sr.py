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

base_out_dir = '../../../../methanol/DMTM_Metals/'
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

plotcase = 'landscapes'

if plotcase == 'volcano':
    bsOs = np.linspace(start=-6, stop=0, num=50, endpoint=True)
elif plotcase == 'drc':
    bsOs = [-3.8, -2, -1.8, -0.8]
else:
    bsOs = np.zeros(len(all_metals))

tofdict = dict()
sim_systems = dict()

for study in ['dry', 'wet']:
    tofdict[study] = np.zeros((3, len(bsOs)))
    for Ti, temperature in enumerate([500, 650, 800]):
        sim_systems[study] = dict()
        results_dir = base_out_dir + 'results/%s/' % study
        figures_dir = base_out_dir + 'images/%s/' % study
        input_file = 'input_%s_sr.json' % study
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

        if plotcase != 'landscapes':
            for bind, bsO in enumerate(bsOs):
                sim_systems[study][bsO] = dict()
                sys.states['sO'].Gelec = bsO
                sys.reactions['rsO'].dErxn_user = bsO

                if plotcase == 'volcano':
                    run(sim_system=sys,
                        plot_results=False,
                        save_results=False,
                        steady_state_solve=True,
                        fig_path=figures_dir + 'sr/%s/' % bsO,
                        csv_path=results_dir + 'sr/%s/' % bsO)
                    sys.reaction_terms(y=sys.full_steady)
                    i5 = list(sys.reactions.keys()).index('r5_rdr')
                    i9 = list(sys.reactions.keys()).index('r9_rdr')
                    tofdict[study][Ti, bind] = (sys.rates[i5, 0] -
                                                sys.rates[i5, 1]) + (sys.rates[i9, 0] -
                                                                     sys.rates[i9, 1])
                elif plotcase == 'drc':
                    if not os.path.isdir(results_dir + 'sr/%s/' % bsO):
                        os.mkdir(results_dir + 'sr/%s/' % bsO)
                        if not os.path.isdir(results_dir + 'sr/%s/%1.0f/' % (bsO, temperature)):
                            os.mkdir(results_dir + 'sr/%s/%1.0f/' % (bsO, temperature))

                    run_temperatures(sim_system=sys,
                                     temperatures=np.linspace(start=temperature,
                                                              stop=temperature,
                                                              num=1,
                                                              endpoint=True),
                                     plot_results=False,
                                     save_results=True,
                                     steady_state_solve=True,
                                     eps=1.0e-3,
                                     tof_terms=['r5_rdr', 'r9_rdr'],
                                     csv_path=results_dir + 'sr/%s/%1.0f/' % (bsO, temperature))
        else:
            for mind, metal in enumerate(all_metals):
                if not os.path.isdir(figures_dir + 'sr/'):
                    os.mkdir(figures_dir + 'sr/')
                if not os.path.isdir(figures_dir + 'sr/%s/' % metal):
                    os.mkdir(figures_dir + 'sr/%s/' % metal)

                sim_systems[study][metal] = dict()
                sim_systems[study][metal]['DFT'] = read_from_input_file(input_path='%s/input_%s.json' % (metal, study))
                sim_systems[study][metal]['DFT'].params['temperature'] = temperature

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
                bsOs[mind] = bsO

                sys.states['sO'].Gelec = bsO
                sys.reactions['rsO'].dErxn_user = bsO
                sim_systems[study][metal]['scaling'] = copy.deepcopy(sys)

                def cmap(ind):
                    if not ind:
                        return meclrs[mind]
                    else:
                        return gclr

                if metal != 'zn':
                    mins = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima[0:-2]
                    labs = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels[0:-2]
                else:
                    mins = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima[0:-1]
                    labs = sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels[0:-1]
                sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].minima = mins
                sim_systems[study][metal]['DFT'].energy_landscapes['full_pes'].labels = labs

                compare_energy_landscapes(sim_systems=sim_systems[study][metal],
                                          etype='electronic',
                                          legend_location='best',
                                          show_labels=True,
                                          fig_path=figures_dir + 'sr/%s/' % metal,
                                          cmap=cmap)
                compare_energy_landscapes(sim_systems=sim_systems[study][metal],
                                          etype='free',
                                          legend_location='best',
                                          show_labels=True,
                                          fig_path=figures_dir + 'sr/%s/' % metal,
                                          cmap=cmap)

    if plotcase == 'volcano':
        df = pd.DataFrame(tofdict[study],
                          columns=bsOs,
                          index=[500, 650, 800])
        results_dir = base_out_dir + 'results/%s/' % study
        df.to_csv(results_dir + 'sr/volcano_tofs.csv')

if plotcase == 'volcano':
    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(3.2, 3.2), sharey='all', sharex='all')
    for Ti, temperature in enumerate([500, 650, 800]):
        for sind, study in enumerate(['dry', 'wet']):
            ax[0 if study == 'wet' else 1][Ti].plot([b for bi, b in enumerate(bsOs) if tofdict[study][Ti, bi] > 0],
                                                    [t for t in tofdict[study][Ti, :] if t > 0],
                                                    marker='o' if study == 'wet' else 'v',
                                                    color='k' if study == 'wet' else 'k',
                                                    linestyle=':' if study == 'wet' else '-',
                                                    label=study, markevery=[])

            for mind, metal in enumerate(all_metals):
                results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
                descriptor = 'sO' if study == 'dry' else 'sO-h2o'
                df = pd.read_csv(results_dir + 'full_pes_energy_landscape_%1.1fK_1.0bar.csv' % temperature)
                EsO = df['Electronic (eV)'][list(df['State']).index('sO')]
                df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
                Tind = list(df['Temperature (K)']).index(temperature)
                ax[0 if study == 'wet' else 1][Ti].plot(EsO, df['r5'][Tind] + df['r9'][Tind],
                                                        color=meclrs[all_metals.index(metal)],
                                                        marker='o' if study == 'wet' else 'v')

    ax[0][0].set(ylabel='TOF (1/s)', yscale='log', ylim=(1e-16, 1e5), yticks=(1e-15, 1e-10, 1e-5, 1e0))
    ax[1][0].set(ylabel='TOF (1/s)', yscale='log', ylim=(1e-16, 1e5), yticks=(1e-15, 1e-10, 1e-5, 1e0),
                 xlabel='E$_{*O}$ (eV)', xticks=(-6, -4, -2), xlim=(-7, 0))
    ax[1][1].set(xlabel='E$_{*O}$ (eV)', xticks=(-6, -4, -2), xlim=(-7, 0), yticks=(1e-15, 1e-10, 1e-5, 1e0))
    ax[1][2].set(xlabel='E$_{*O}$ (eV)', xticks=(-6, -4, -2), xlim=(-7, 0),
                 yticks=(1e-15, 1e-10, 1e-5, 1e0), ylim=(1e-16, 1e5))
    ax[0][0].set(title='500 K')
    ax[0][1].set(title='650 K')
    ax[0][2].set(title='800 K')
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(base_out_dir + 'images/scaling_tof.png', format='png', dpi=600)
    fig.savefig(base_out_dir + 'images/scaling_tof.eps', format='eps', dpi=600)
elif plotcase == 'drc':
    cmap = plt.get_cmap("tab20", 12)
    for study in ['dry', 'wet']:
        results_dir = base_out_dir + 'results/%s/' % study
        figures_dir = base_out_dir + 'images/%s/' % study
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.4, 3.2), sharey='row')
        for ti, temperature in enumerate([500, 650, 800]):
            be = [-3.8, -1.8, -0.8] if study == 'dry' else [-3.8, -2, -0.8]
            for mi, metal in enumerate(be):
                df = pd.read_csv(results_dir + 'sr/%s/%1.0f/drcs_vs_temperature.csv' % (metal, temperature))
                rnames = [r for r in df.columns[1::].values if r != 'rsO']
                drcs = [df[c].values[0] for c in rnames]
                for xi, x in enumerate(np.linspace(start=mi, stop=mi + 1 - 0.5 / (len(rnames) + 1), num=len(rnames))):
                    ax[ti].bar(x + 0.6 / (len(rnames) + 1), drcs[xi], width=0.9 / (len(rnames) + 1), color=cmap(xi),
                               label=rnames[xi].split('_')[0] if mi == 0 else '')
                if mi < 2:
                    ax[ti].plot((mi + 1, mi + 1), (-1, 1.2), ':', color=gclr)
            ax[ti].set(title='%1.0f K' % temperature, xlabel='*O energy [eV]', ylim=(0, 1.5), xlim=(0, 3),
                       xticks=(0.5, 1.5, 2.5), xticklabels=be)
        ax[0].set(ylabel='DRC', ylim=(0, 1.2))
        fig.tight_layout()
        fig.savefig(base_out_dir + 'images/drc_%s.eps' % study, format='eps', dpi=600)
        ax[0].set(ylabel='DRC', ylim=(0, 1.5))
        ax[1].legend(ncol=3, mode='expand', loc='upper center')
        fig.tight_layout()
        fig.savefig(base_out_dir + 'images/drc_%s.png' % study, format='png', dpi=600)

    for study in ['dry', 'wet']:
        results_dir = base_out_dir + 'results/%s/' % study
        figures_dir = base_out_dir + 'images/%s/' % study
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.4, 3.2), sharey='row')
        for ti, temperature in enumerate([500, 650, 800]):
            be = [-3.8, -1.8, -0.8] if study == 'dry' else [-3.8, -2, -0.8]
            labs = []
            for mi, metal in enumerate(be):
                df = pd.read_csv(results_dir + 'sr/%s/%1.0f/coverages_vs_temperature.csv' % (metal, temperature))
                rnames = [r for r in df.columns[1::].values if '0' not in r]
                drcs = [df[c].values[0] for c in rnames]
                for xi, x in enumerate(np.linspace(start=mi, stop=mi + 1 - 0.5 / (len(rnames) + 1), num=len(rnames))):
                    if drcs[xi] > 0:
                        ax[ti].bar(x + 0.6 / (len(rnames) + 1), drcs[xi], width=0.9 / (len(rnames) + 1), color=cmap(xi),
                                   label=rnames[xi].split('_')[0] if rnames[xi].split('_')[0] not in labs else '')
                        if rnames[xi].split('_')[0] not in labs:
                            labs += [rnames[xi].split('_')[0]]
                if mi < 2:
                    ax[ti].plot((mi + 1, mi + 1), (-1, 1.2), ':', color=gclr)
            ax[ti].set(title='%1.0f K' % temperature, xlabel='*O energy [eV]', ylim=(0, 1.5), xlim=(0, 3),
                       xticks=(0.5, 1.5, 2.5), xticklabels=be)
        ax[0].set(ylabel='Coverage', ylim=(0, 1.2))
        fig.tight_layout()
        fig.savefig(base_out_dir + 'images/coverage_%s.eps' % study, format='eps', dpi=600)
        ax[0].legend(ncol=1, loc='upper center')
        ax[1].legend(ncol=1, loc='upper center')
        ax[2].legend(ncol=1, loc='upper center')
        ax[0].set(ylabel='Coverage', ylim=(0, 1.5))
        fig.tight_layout()
        fig.savefig(base_out_dir + 'images/coverage_%s.png' % study, format='png', dpi=600)
