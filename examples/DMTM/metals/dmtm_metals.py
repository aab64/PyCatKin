from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import *
from ase.io import read
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
from scipy.stats import t

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
temperatures = np.linspace(start=400, stop=800, num=17, endpoint=True)

sim_systems = dict()

for study in ['dry', 'wet']:
    all_states = []
    sim_systems[study] = dict()
    for metal in all_metals:
        print(metal)
        # Location of files and images
        input_file = 'input_%s.json' % study
        # base_data_dir = 'D:/Users/Astrid/Documents/Chalmers/Data/Methanol/DMTM_Metals/'
        # gas_opt_dir = base_data_dir + 'energies/molecules/'
        # gas_vib_dir = base_data_dir + 'vibrations/molecules/'
        # ads_opt_dir = base_data_dir + 'energies/%s/%s/' % (study, metal)
        # ads_vib_dir = base_data_dir + 'vibrations/%s/%s/' % (study, metal)
        base_out_dir = 'D:/Users/Astrid/Dropbox/Chalmers/Simulations/microkinetics/Methanol/DMTM_Metals/'
        results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
        figures_dir = base_out_dir + 'images/%s/%s/' % (study, metal)
        sim_systems[study][metal] = read_from_input_file(input_path=metal + '/' + input_file)

        for s in sim_systems[study][metal].snames:
            # opt_dir = gas_opt_dir if sim_systems[study][metal].states[s].state_type == 'gas' else ads_opt_dir
            # vib_dir = gas_vib_dir if sim_systems[study][metal].states[s].state_type == 'gas' else ads_vib_dir
            # sim_systems[study][metal].states[s].path = opt_dir + sim_systems[study][metal].states[s].name
            # sim_systems[study][metal].states[s].vibs_path = vib_dir + sim_systems[study][metal].states[s].name
            # if not os.path.isdir(sim_systems[study][metal].states[s].vibs_path):
            #     print(sim_systems[study][metal].states[s].vibs_path)
            #     sim_systems[study][metal].states[s].vibs_path = (base_data_dir +
            #                                                      'vibrations/tmp-old-Cu/%s/' % study +
            #                                                      sim_systems[study][metal].states[s].name)
            # read_from_alternate = None
            # if sim_systems[study][metal].states[s].state_type == 'gas':
            #     read_from_alternate = {'get_atoms': lambda state_path=opt_dir + sim_systems[study][metal].states[s].name: load_from_contcar(
            #         state_path)}
            # sim_systems[study][metal].states[s].read_from_alternate = read_from_alternate

            if s not in all_states:
                all_states.append(s)

        run_temperatures(sim_system=sim_systems[study][metal],
                         temperatures=temperatures,
                         steady_state_solve=True,
                         tof_terms=None,
                         plot_results=False,
                         save_results=True,
                         fig_path=figures_dir,
                         csv_path=results_dir)
        # run_energy_span_temperatures(sim_system=sim_systems[study][metal],
        #                              temperatures=temperatures,
        #                              etype='free',
        #                              save_results=True,
        #                              csv_path=results_dir)
        #
        # save_energies_temperatures(sim_system=sim_systems[study][metal],
        #                            temperatures=temperatures,
        #                            csv_path=results_dir)

        sim_systems[study][metal].params['temperature'] = 500
        # draw_energy_landscapes(sim_system=sim_systems[study][metal],
        #                        etype='free',
        #                        show_labels=True,
        #                        fig_path=figures_dir)
        # save_state_energies(sim_system=sim_systems[study][metal],
        #                     csv_path=results_dir)
        # save_energies(sim_system=sim_systems[study][metal],
        #               csv_path=results_dir)
        save_pes_energies(sim_system=sim_systems[study][metal],
                          csv_path=results_dir)

    all_states.sort()

    # # TOF MK vs ES
    # fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex='all', sharey='all')
    # locs = [[0, 0], [0, 1], [0, 2], [0, 3],
    #         [1, 0], [1, 1], [1, 2], [1, 3],
    #         [2, 0], [2, 1], [2, 2], [2, 3]]
    # for mind, metal in enumerate(all_metals):
    #     results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
    #     df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
    #     ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values,
    #                                           df['r5'].values + df['r9'].values,
    #                                           color=meclrs[all_metals.index(metal)],
    #                                           label='MK')
    #     df = pd.read_csv(results_dir + 'energy_span_summary_full_pes.csv')
    #     ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values,
    #                                           2 * df['TOF (1/s)'].values,
    #                                           color=meclrs[all_metals.index(metal)],
    #                                           label='ES',
    #                                           linestyle=':')
    #     ax[locs[mind][0]][locs[mind][1]].set(yscale='log', ylim=(1e-20, 1e2))
    #     if locs[mind][0] == 2:
    #         ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    #     else:
    #         ax[locs[mind][0]][locs[mind][1]].tick_params(axis='x', which='both',
    #                                                      bottom=False, top=False, labelbottom=False)
    #     if locs[mind][1] == 0:
    #         ax[locs[mind][0]][locs[mind][1]].set(ylabel='TOF (1/s)')
    #     else:
    #         ax[locs[mind][0]][locs[mind][1]].tick_params(axis='y', which='both',
    #                                                      left=False, right=False, labelleft=False)
    #     ax[locs[mind][0]][locs[mind][1]].set(title=metal)
    # for s in ['MK', 'ES']:
    #     ax[-1, -1].plot((temperatures[0], temperatures[-1]), (0, 0),
    #                     label=s,
    #                     color='k',
    #                     linestyle='-' if s == 'MK' else ':')
    # ax[-1, -1].legend(loc='center', frameon=False)
    # ax[-1, -1].tick_params(axis='both', which='both',
    #                        bottom=False, top=False, labelbottom=False,
    #                        left=False, right=False, labelleft=False)
    # ax[-1, -1].axis('off')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0.2, hspace=0.2)
    # plt.savefig(base_out_dir + 'images/%s/TOF_vs_temperature_each_metal.png' % study, format='png', dpi=600)
    #
    # # TDTS
    # fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex='all', sharey='all')
    # labs = []
    # cmap = plt.get_cmap("Spectral", len(all_states) + 2)
    # for mind, metal in enumerate(all_metals):
    #     results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
    #     df = pd.read_csv(results_dir + 'energy_span_xTDTS_full_pes.csv')
    #     for sind, s in enumerate(df.columns[1::]):
    #         if s not in all_states:
    #             all_states.append(s)
    #         if max(df[s].values) > 0.01:
    #             ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values, df[s].values,
    #                                                   color=cmap(all_states.index(s)), label=s)
    #             if s not in labs:
    #                 labs.append(s)
    #     if locs[mind][0] == 2:
    #         ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    #     if locs[mind][1] == 0:
    #         ax[locs[mind][0]][locs[mind][1]].set(ylabel='xTDTS')
    #     ax[locs[mind][0]][locs[mind][1]].set(title=metal)
    # for s in labs:
    #     ax[-1, -1].plot((temperatures[0], temperatures[-1]), (-1, -1),
    #                     label=s.split('-h2o')[0],
    #                     color=cmap(all_states.index(s)))
    # ax[-1, -1].legend(loc='center', frameon=False, mode='expand', ncol=2, fontsize=5)
    # ax[-1, -1].set(ylim=(-0.05, 1.05))
    # plt.tick_params(axis='both', which='both',
    #                 bottom=False, top=False, labelbottom=False,
    #                 left=False, right=False, labelleft=False)
    # ax[-1, -1].axis('off')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0.2, hspace=0.2)
    # plt.savefig(base_out_dir + 'images/%s/xTDTS_vs_temperature_vs_metals.png' % study, format='png', dpi=600)
    #
    # # TDI
    # fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex='all', sharey='all')
    # labs = []
    # for mind, metal in enumerate(all_metals):
    #     results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
    #     df = pd.read_csv(results_dir + 'energy_span_summary_full_pes.csv')
    #     if metal == 'pd' or metal == 'ag':
    #         print(metal)
    #         print(df)
    #     df = pd.read_csv(results_dir + 'energy_span_xTDI_full_pes.csv')
    #     for sind, s in enumerate(df.columns[1::]):
    #         if s not in all_states:
    #             all_states.append(s)
    #         if max(df[s].values) > 0.01:
    #             ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values, df[s].values,
    #                                                   color=cmap(all_states.index(s)), label=s)
    #             if s not in labs:
    #                 labs.append(s)
    #     if locs[mind][0] == 2:
    #         ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    #     if locs[mind][1] == 0:
    #         ax[locs[mind][0]][locs[mind][1]].set(ylabel='xTDI')
    #     ax[locs[mind][0]][locs[mind][1]].set(title=metal)
    # for s in labs:
    #     ax[-1, -1].plot((temperatures[0], temperatures[-1]), (-1, -1),
    #                     label=s.split('-h2o')[0],
    #                     color=cmap(all_states.index(s)))
    # ax[-1, -1].legend(loc='center', frameon=False, ncol=2, fontsize=5, mode='expand')
    # ax[-1, -1].set(ylim=(-0.05, 1.05))
    # plt.tick_params(axis='both', which='both',
    #                 bottom=False, top=False, labelbottom=False,
    #                 left=False, right=False, labelleft=False)
    # ax[-1, -1].axis('off')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0.2, hspace=0.2)
    # plt.savefig(base_out_dir + 'images/%s/xTDI_vs_temperature_vs_metals.png' % study, format='png', dpi=600)
    #
    # # Coverage
    # fig, ax = plt.subplots(figsize=(6.4, 4.8), ncols=4, nrows=3, sharex='all', sharey='all')
    # labs = []
    # for mind, metal in enumerate(all_metals):
    #     results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
    #     df = pd.read_csv(results_dir + 'coverages_vs_temperature.csv')
    #     for sind, s in enumerate(df.columns[1::]):
    #         if s not in all_states:
    #             all_states.append(s)
    #         if max(df[s].values) > 0.01:
    #             ax[locs[mind][0]][locs[mind][1]].plot(df['Temperature (K)'].values, df[s].values,
    #                                                   color=cmap(all_states.index(s)), label=s)
    #             if s not in labs:
    #                 labs.append(s)
    #     if locs[mind][0] == 2:
    #         ax[locs[mind][0]][locs[mind][1]].set(xlabel='Temperature (K)')
    #     if locs[mind][1] == 0:
    #         ax[locs[mind][0]][locs[mind][1]].set(ylabel='Coverage')
    #     ax[locs[mind][0]][locs[mind][1]].set(title=metal)
    # for s in labs:
    #     ax[-1, -1].plot((temperatures[0], temperatures[-1]), (-1, -1),
    #                     label=s.split('-h2o')[0],
    #                     color=cmap(all_states.index(s)))
    # ax[-1, -1].legend(loc='center', frameon=False, mode='expand', ncol=2, fontsize=5)
    # ax[-1, -1].set(ylim=(-0.05, 1.05))
    # plt.tick_params(axis='both', which='both',
    #                 bottom=False, top=False, labelbottom=False,
    #                 left=False, right=False, labelleft=False)
    # ax[-1, -1].axis('off')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0.2, hspace=0.2)
    # plt.savefig(base_out_dir + 'images/%s/coverage_vs_temperature_vs_metals.png' % study, format='png', dpi=600)

    # Scaling relation
    all_states = ["2s", "ts1", "s-pair", "sO2s", "ts2", "sOOs", "s2Och4", "rad1", "sOsCH3OH", "sO",
                  "sOch4", "rad2", "sOHsCH3", "ts5", "sCH3OH", "s", "ts6"]
    state_names = ["2Me", "TS1", "Me-pair", "MeO$_2$Me", "TS2", "MeOOMe", "*2OCH$_4$", "TS3", "*O*CH$_3$OH", "*O",
                   "*OCH$_4$", "TS4", "*OH*CH$_3$", "TS5", "*CH$_3$OH", "*", "TS6"]
    descriptor = 'sO'
    Ttopplot = 500.0

    locs = [[0, 0], [0, 1], [0, 2], [0, 3],
            [1, 0], [1, 1], [1, 2], [1, 3],
            [2, 0], [2, 1], [2, 2], [2, 3],
            [3, 0], [3, 1], [3, 2], [3, 3]]

    xvals = dict()
    yvals = dict()
    mvals = dict()
    dfvals = dict({'Gradient': dict(), 'Std err': dict(), '95% CI': dict(),
                   'Intercept': dict(), 'R2': dict(), 'MAE': dict(), 'MAX': dict()})

    subset = ['ag', 'zn']

    for mind, metal in enumerate([m for m in all_metals if m not in subset]):  #
        results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
        for s in all_states:
            if s in sim_systems[study][metal].energy_landscapes['full_pes'].labels and s != descriptor:
                if s not in xvals.keys():
                    xvals[s] = []
                if s not in yvals.keys():
                    yvals[s] = []
                if s not in mvals.keys():
                    mvals[s] = []
                mvals[s].append(metal)
                df = pd.read_csv(results_dir + 'full_pes_energy_landscape_500.0K_1.0bar.csv')
                xvals[s].append(df['Electronic (eV)'][list(df['State']).index(descriptor)])
                yvals[s].append(df['Electronic (eV)'][list(df['State']).index(s)])

    xmin = -7  # min([min(xvals[k]) for k in xvals.keys()]) - 0.5
    xmax = -1  # max([max(xvals[k]) for k in xvals.keys()]) + 0.5
    ymin = -8  # min([min(yvals[k]) for k in yvals.keys()]) - 0.5
    ymax = 2  # max([max(yvals[k]) for k in yvals.keys()]) + 0.5

    fig, ax = plt.subplots(ncols=4, nrows=4, figsize=(6.4, 6.4), sharex='all', sharey='row')
    for k in xvals.keys():
        kind = all_states.index(k)
        if kind > all_states.index(descriptor):
            kind -= 1
        mkr = 'o'
        lsy = ':'
        for mind, m in enumerate(mvals[k]):
            ax[locs[kind][0]][locs[kind][1]].plot(xvals[k][mind], yvals[k][mind],
                                                  marker=mkr, color=meclrs[all_metals.index(m)])
        if locs[kind][1] == 0:
            ax[locs[kind][0]][locs[kind][1]].set(ylabel='dE$_{state}$ (eV)')
            ax[locs[kind][0]][locs[kind][1]].tick_params(axis='y', which='both',
                                                         left=False, right=False)
        else:
            ax[locs[kind][0]][locs[kind][1]].tick_params(axis='y', which='both',
                                                         left=False, right=False, labelleft=False)
        if locs[kind][0] == 3:
            ax[locs[kind][0]][locs[kind][1]].set(xlabel='dE$_{*O}$ (eV)')
            ax[locs[kind][0]][locs[kind][1]].tick_params(axis='x', which='both',
                                                         bottom=False, top=False)
        else:
            ax[locs[kind][0]][locs[kind][1]].tick_params(axis='x', which='both',
                                                         bottom=False, top=False, labelbottom=False)

        xtofit = np.array([xvals[k][i] for i in range(len(xvals[k]))])
        ytofit = np.array([yvals[k][i] for i in range(len(yvals[k]))])
        slope, intercept, r, p, se = linregress(xtofit, ytofit)
        yfit = slope * xtofit + intercept

        ax[locs[kind][0]][locs[kind][1]].plot(np.array((min(xmin, xmin), max(xmax, xmax))),
                                              np.array((min(xmin, xmin), max(xmax, xmax))) * slope + intercept,
                                              linestyle=lsy, color=gclr, linewidth=2)

        tinv = lambda p, dof: abs(t.ppf(p / 2, dof))
        ts = tinv(0.05, len(xtofit) - 2)

        dfvals['Gradient'][k] = slope
        dfvals['95% CI'][k] = ts * se
        dfvals['Intercept'][k] = intercept
        dfvals['R2'][k] = r**2
        dfvals['MAE'][k] = np.mean(abs(ytofit - yfit))
        dfvals['MAX'][k] = np.max(abs(ytofit - yfit))

        # ax[locs[kind][0]][locs[kind][1]].text(xmax - 0.25,
        #                                       ymin + 0.25,
        #                                       state_names[all_states.index(k)],
        #                                       ha='right', va='bottom', color=gclr)

        # ax[locs[kind][0]][locs[kind][1]].text(xmax - 0.25,
        #                                       ymin + 0.25,
        #                                       '$R^2=%1.2f$' % (sol.rvalue**2),
        #                                       ha='right', va='bottom', color=gclr)
        # if 's-pair.1' not in k:
        #     ax[locs[kind][0]][locs[kind][1]].text(max(xmax, ymax) - 0.25,
        #                                           min(xmin, ymin) + 0.25,
        #                                           '$%1.2f, %1.2f$' % (sol.slope, sol.intercept),
        #                                           ha='right', va='bottom', color=gclr)
        # else:
        #     ax[locs[kind][0]][locs[kind][1]].text(max(xmax, ymax) - 0.25,
        #                                           min(xmin, ymin) + 0.75,
        #                                           '$%1.2f, %1.2f$' % (sol.slope, sol.intercept),
        #                                           ha='right', va='bottom', color=gclr)

        ax[locs[kind][0]][locs[kind][1]].set(xlim=(min(xmin, xmin), max(xmax, xmax)),
                                             ylim=(min(ymin, ymin), max(ymax, ymax)),
                                             xticks=(-6, -4, -2), yticks=(-6, -4, -2, 0))
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    print('%1.2f, %1.2f, %1.2f, %1.2f' % (xmin, xmax, ymin, ymax))
    plt.savefig(base_out_dir + 'images/%s/SR_dE%s_vs_dEstate_T=%1.0fK.png' %
                (study, descriptor, Ttopplot),
                format='png', dpi=600)
    plt.savefig(base_out_dir + 'images/%s/SR_dE%s_vs_dEstate_T=%1.0fK.eps' %
                (study, descriptor, Ttopplot),
                format='eps', dpi=600)
    df = pd.DataFrame(dfvals)
    df.to_csv(path_or_buf=base_out_dir + 'results/%s/scaling.csv' % study)

# TOF comparison wet/dry
mnames = ['Ag$_2$', 'Au$_2$', 'AuCu$_2$', 'AuPd$_2$', 'Co$_2$',
          'Cu$_2$', 'Fe$_2$', 'Ni$_2$', 'Pd$_2$', 'PdCu$_2$', 'Zn$_2$']

fig, ax = plt.subplots(ncols=3, figsize=(6.4, 3.2), sharey='all',
                       gridspec_kw={'width_ratios': [1.2, 1.2, 1]})

ax[0].plot(temperatures, np.ones(len(temperatures)), color=gclr, linewidth=1.5)
ax[1].plot(temperatures, np.ones(len(temperatures)), color=gclr, linewidth=1.5)
ax[0].text(405, 1.0, 'TOF = 1 s$^{-1}$', ha='left', va='bottom', color=gclr)

for sind, study in enumerate(['dry', 'wet']):
    for mind, metal in enumerate(all_metals):
        results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
        df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
        if metal not in ['aucu', 'aupd', 'pdcu']:
            ax[0].plot(df['Temperature (K)'].values,
                       df['r5'].values + df['r9'].values,
                       color=meclrs[all_metals.index(metal)],
                       label=mnames[mind] if study == 'dry' else 'w /H_2O',
                       linestyle='-' if study == 'dry' else ':',
                       linewidth=2)
        else:
            ax[1].plot(df['Temperature (K)'].values,
                       df['r5'].values + df['r9'].values,
                       color=meclrs[all_metals.index(metal)],
                       label=mnames[mind] if study == 'dry' else 'w /H_2O',
                       linestyle='-' if study == 'dry' else ':',
                       linewidth=2)
        # ax[2].plot(df['Temperature (K)'].values,
        #            1e-16 * np.ones(17),
        #            color=meclrs[all_metals.index(metal)],
        #            label=mnames[mind] if study == 'dry' else 'w /H$_2$O',
        #            linestyle='-' if study == 'dry' else ':',
        #            linewidth=2)
# leg = ax[2].legend(loc='center', frameon=False,
#                    ncol=2, mode='expand', labelspacing=0.8)
ax[0].set(xlabel='Temperature (K)', xlim=(400, 800), xticks=(450, 550, 650, 750),
          ylabel='TOF (1/s)', yscale='log', ylim=(1e-16, 2e2),
          title='Metals')
ax[1].set(xlabel='Temperature (K)', xlim=(400, 800), xticks=(450, 550, 650, 750),
          yscale='log', ylim=(1e-16, 2e2),
          title='Alloys')
ax[2].axis('off')
ax[0].tick_params(axis='y', which='right',
                  right=False, labelright=False)
ax[1].tick_params(axis='y', which='both',
                  left=False, right=False, labelleft=False)
fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.savefig(base_out_dir + 'images/tof_vs_temperature_vs_metals.png', format='png', dpi=600)
plt.savefig(base_out_dir + 'images/tof_vs_temperature_vs_metals.eps', format='eps', dpi=600)
