from pycatkin.functions.load_input import read_from_input_file
from pycatkin.functions.presets import *
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

base_out_dir = '../../../../methanol/DMTM_Metals/'

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
all_metals = ['ag', 'au', 'aucu', 'aupd', 'co', 'cu', 'fe', 'ni', 'pd', 'pdcu', 'zn']
metal_names = ['Ag', 'Au', 'AuCu', 'AuPd', 'Co', 'Cu', 'Fe', 'Ni', 'Pd', 'PdCu', 'Zn']
all_states = ["2s", "ts1", "s-pair", "sO2s", "ts2", "sOOs", "s2Och4", "rad1", "sOsCH3OH", "sO",
              "sOch4", "rad2", "sOHsCH3", "ts5", "sCH3OH", "s", "ts6"]
state_names = ["2Me", "TS1", "Me-pair", "MeO$_2$Me", "TS2", "MeOOMe", "*2OCH$_4$", "TS3", "*O*CH$_3$OH", "*O",
               "*OCH$_4$", "TS4", "*OH*CH$_3$", "TS5", "*CH$_3$OH", "*", "TS6"]

temperatures = np.linspace(start=400, stop=800, num=17, endpoint=True)
sim_systems = dict()
for study in ['dry', 'wet']:
    sim_systems[study] = dict()
    for metal in all_metals:
        print(metal)
        results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
        figures_dir = base_out_dir + 'images/%s/%s/' % (study, metal)
        sim_systems[study][metal] = read_from_input_file(input_path='%s/input_%s.json' % (metal, study))

        run_temperatures(sim_system=sim_systems[study][metal],
                         temperatures=temperatures,
                         steady_state_solve=True,
                         plot_results=True,
                         save_results=True,
                         fig_path=figures_dir,
                         csv_path=results_dir)

        for temperature in [500, 650, 800]:
            sim_systems[study][metal].params['temperature'] = temperature
            save_pes_energies(sim_system=sim_systems[study][metal],
                              csv_path=results_dir)

    # Plot scaling relations
    descriptor = 'sO'
    locs = [[0, 0], [0, 1], [0, 2], [0, 3],
            [1, 0], [1, 1], [1, 2], [1, 3],
            [2, 0], [2, 1], [2, 2], [2, 3],
            [3, 0], [3, 1], [3, 2], [3, 3]]
    xvals = dict()
    yvals = dict()
    mvals = dict()
    dfvals = dict({'Gradient': dict(), '95% CI': dict(), 'Intercept': dict(),
                   'R2': dict(), 'MAE': dict(), 'MAX': dict()})
    for mind, metal in enumerate([m for m in all_metals if m not in ['ag', 'zn']]):
        results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
        df = pd.read_csv(results_dir + 'full_pes_energy_landscape_500.0K_1.0bar.csv')
        for s in all_states:
            if s in df['State'].values and s != descriptor:
                if s not in xvals.keys():
                    xvals[s] = []
                if s not in yvals.keys():
                    yvals[s] = []
                if s not in mvals.keys():
                    mvals[s] = []
                mvals[s].append(metal)
                xvals[s].append(df['Electronic (eV)'][list(df['State']).index(descriptor)])
                yvals[s].append(df['Electronic (eV)'][list(df['State']).index(s)])
    xmin = -7
    xmax = -1
    ymin = -8
    ymax = 2
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
        ax[locs[kind][0]][locs[kind][1]].plot(np.array((xmin, xmax)), np.array((xmin, xmax)) * slope + intercept,
                                              linestyle=lsy, color=gclr, linewidth=2)
        tinv = lambda p, dof: abs(t.ppf(p / 2, dof))
        ts = tinv(0.05, len(xtofit) - 2)
        dfvals['Gradient'][k] = slope
        dfvals['95% CI'][k] = ts * se
        dfvals['Intercept'][k] = intercept
        dfvals['R2'][k] = r**2
        dfvals['MAE'][k] = np.mean(abs(ytofit - yfit))
        dfvals['MAX'][k] = np.max(abs(ytofit - yfit))
        ax[locs[kind][0]][locs[kind][1]].text(xmax - 0.1, ymin, '$y=%1.2fx+%1.2f$' % (slope, intercept),
                                              ha='right', va='bottom', color=gclr)
        ax[locs[kind][0]][locs[kind][1]].text(xmax - 0.1, ymin + 0.75, state_names[kind],
                                              ha='right', va='bottom', color=gclr)

        ax[locs[kind][0]][locs[kind][1]].set(xlim=(xmin, xmax), ylim=(ymin, ymax),
                                             xticks=(-6, -4, -2), yticks=(-6, -4, -2, 0))
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(base_out_dir + 'images/%s/SR_dE%s_vs_dEstate.png' % (study, descriptor), format='png', dpi=600)
    plt.savefig(base_out_dir + 'images/%s/SR_dE%s_vs_dEstate.eps' % (study, descriptor), format='eps', dpi=600)
    df = pd.DataFrame(dfvals)
    df.to_csv(path_or_buf=base_out_dir + 'results/%s/scaling.csv' % study)

# Plot TOF comparison wet/dry
fig, ax = plt.subplots(ncols=3, figsize=(6.4, 3.2), sharey='all',
                       gridspec_kw={'width_ratios': [1.2, 1.2, 1]})
ax[0].plot(temperatures, np.ones(len(temperatures)), color=gclr, linewidth=1.5)
ax[1].plot(temperatures, np.ones(len(temperatures)), color=gclr, linewidth=1.5)
ax[1].text(405, 1.0, 'TOF = 1 s$^{-1}$', ha='left', va='bottom', color=gclr)
for sind, study in enumerate(['dry', 'wet']):
    for mind, metal in enumerate(all_metals):
        results_dir = base_out_dir + 'results/%s/%s/' % (study, metal)
        df = pd.read_csv(results_dir + 'rates_vs_temperature.csv')
        if metal not in ['aucu', 'aupd', 'pdcu']:
            ax[0].plot(df['Temperature (K)'].values,
                       df['r5'].values + df['r9'].values,
                       color=meclrs[all_metals.index(metal)],
                       label=metal_names[mind] if study == 'dry' else 'w /H_2O',
                       linestyle='-' if study == 'dry' else ':',
                       linewidth=2)
        else:
            ax[1].plot(df['Temperature (K)'].values,
                       df['r5'].values + df['r9'].values,
                       color=meclrs[all_metals.index(metal)],
                       label=metal_names[mind] if study == 'dry' else 'w /H_2O',
                       linestyle='-' if study == 'dry' else ':',
                       linewidth=2)
        ax[2].plot(df['Temperature (K)'].values,
                   1e-20 * np.ones(17),
                   color=meclrs[all_metals.index(metal)],
                   label=metal_names[mind] if study == 'dry' else 'w /H$_2$O',
                   linestyle='-' if study == 'dry' else ':',
                   linewidth=2)
ylims = (1e-16, 1e5)
leg = ax[2].legend(loc='center', frameon=False,
                   ncol=2, mode='expand', labelspacing=0.8)
ax[0].set(xlabel='Temperature (K)', xlim=(400, 800), xticks=(450, 550, 650, 750),
          ylabel='TOF (1/s)', yscale='log', ylim=ylims, title='Metals')
ax[1].set(xlabel='Temperature (K)', xlim=(400, 800), xticks=(450, 550, 650, 750),
          yscale='log', ylim=ylims, title='Alloys')
ax[2].axis('off')
ax[0].tick_params(axis='y', which='right',
                  right=False, labelright=False)
ax[1].tick_params(axis='y', which='both',
                  left=False, right=False, labelleft=False)
fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.savefig(base_out_dir + 'images/tof_vs_temperature_vs_metals.png', format='png', dpi=600)
plt.savefig(base_out_dir + 'images/tof_vs_temperature_vs_metals.eps', format='eps', dpi=600)

# for me in all_metals:
#     print(me)
#     for study in ['dry', 'wet']:
#         print(study)
#         base_dir = '../../../../../../../../Documents/Chalmers/Data/Methanol/DMTM_Metals/vibrations/%s/%s/' % (study,
#                                                                                                                me)
#         for s in all_states:
#             if study == 'wet':
#                 s += '-h2o'
#             if s in sim_systems[study][me].snames and sim_systems[study][me].states[s].state_type != 'gas':
#                 sim_systems[study][me].states[s].freq = None
#                 sim_systems[study][me].states[s].i_freq = None
#                 sim_systems[study][me].states[s].vibs_path = base_dir + s
#                 if (os.path.isdir(sim_systems[study][me].states[s].vibs_path)):
#                     sim_systems[study][me].states[s].freq_source = None
#                     sim_systems[study][me].states[s].get_free_energy(T=500, p=1e5, verbose=False)
#                     print(s + ': %1.0f' % len(sim_systems[study][me].states[s].freq))
#                     sim_systems[study][me].states[s].save_vibrations(vibs_path='%s/%sdata/vibrations/' % (me,
#                                                                                                           study
#                                                                                                           if study == 'wet'
#                                                                                                           else ''))
#
# for me in all_metals:
#     print(me)
#     for study in ['dry', 'wet']:
#         print(study)
#         base_dir = '../../../../../../../../Documents/Chalmers/Data/Methanol/DMTM_Metals/energies/%s/%s/' % (study, me)
#         for s in all_states:
#             if study == 'wet':
#                 s += '-h2o'
#             if s in sim_systems[study][me].snames and sim_systems[study][me].states[s].state_type != 'gas':
#                 if not (study == 'dry' and me == 'co' and s == 'ts6'):
#                     if sim_systems[study][me].states[s].Gelec is not None:
#                         print(s + ': %1.1f' % sim_systems[study][me].states[s].Gelec)
#                     sim_systems[study][me].states[s].Gelec = None
#                     sim_systems[study][me].states[s].path = base_dir + s + '/OUTCAR'
#                     if (os.path.isfile(sim_systems[study][me].states[s].path)):
#                         sim_systems[study][me].states[s].energy_source = None
#                         sim_systems[study][me].states[s].get_free_energy(T=500, p=1e5, verbose=False)
#                         print(s + ': %1.1f' % sim_systems[study][me].states[s].Gelec)
#                         sim_systems[study][me].states[s].save_energy(path='%s/%sdata/energies/' % (me,
#                                                                                                    study
#                                                                                                    if study == 'wet'
#                                                                                                    else ''))

# for study in ['dry', 'wet']:
#     locs = [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5],
#             [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5],
#             [2, 0], [2, 1], [2, 2], [2, 3], [2, 4], [2, 5]]
#     fig, ax = plt.subplots(ncols=6, nrows=3, figsize=(8, 4), sharex='all', sharey='all')
#     for sind, s in enumerate(all_states):
#         for mind, metal in enumerate(all_metals):
#             if s in sim_systems[study][metal].energy_landscapes['full_pes'].labels:
#                 key = s + '-h2o' if study == 'wet' else s
#                 barval = sim_systems[study][metal].states[key].Gvibr
#             else:
#                 barval = 0.0
#             ax[locs[sind][0]][locs[sind][1]].bar(mind + 1, barval, color=meclrs[all_metals.index(metal)])
#         if locs[sind][1] == 0:
#             ax[locs[sind][0]][locs[sind][1]].tick_params(axis='y', which='both',
#                                                          left=False, right=False)
#             ax[locs[sind][0]][locs[sind][1]].set(ylabel='Gvibr [eV]')
#         else:
#             ax[locs[sind][0]][locs[sind][1]].tick_params(axis='y', which='both',
#                                                          left=False, right=False, labelleft=False)
#         if locs[sind][0] == 2:
#             ax[locs[sind][0]][locs[sind][1]].set(xticks=np.arange(1, len(all_metals) + 1),
#                                                  xticklabels=metal_names)
#             ax[locs[sind][0]][locs[sind][1]].tick_params(axis='x', which='both', rotation=90,
#                                                          bottom=False, top=False)
#         else:
#             ax[locs[sind][0]][locs[sind][1]].tick_params(axis='x', which='both',
#                                                          bottom=False, top=False, labelbottom=False)
#         ax[locs[sind][0]][locs[sind][1]].text(6, 0.9, state_names[sind],
#                                               ha='center', va='top', color=gclr)
#
#     ax[-1][-1].set(xticks=np.arange(1, len(all_metals) + 1), xticklabels=metal_names,
#                    xlim=(0, 12), ylim=(-1.2, 1))
#     ax[-1][-1].tick_params(axis='x', which='both', rotation=90, bottom=False, top=False)
#
#     fig.tight_layout()
#     fig.subplots_adjust(wspace=0, hspace=0)
#     plt.savefig(base_out_dir + 'images/%s/Gvibr.png' % study, format='png', dpi=600)
#     plt.savefig(base_out_dir + 'images/%s/Gvibr.eps' % study, format='eps', dpi=600)
