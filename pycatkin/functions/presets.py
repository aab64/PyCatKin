import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
import os

font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['lines.linewidth'] = 1.5


def run(sim_system, plot_results=False, save_results=False,
        fig_path=None, csv_path=''):
    """Runs the ODE solver and optionally plots/saves
    the results.

    """
    sim_system.solve_odes()
    if plot_results:
        sim_system.plot_transient(path=fig_path)
    if save_results:
        sim_system.write_results(path=csv_path)


def run_temperatures(sim_system, temperatures, steady_state_solve=False, tof_terms=None, eps=1.0e-3,
                     plot_results=False, save_results=False, fig_path=None, csv_path=''):
    """Runs the ODE solver for a range of temperatures
    and optionally plots/saves the results.

    """

    rates = np.zeros((len(temperatures), len(sim_system.reactions)))
    final = np.zeros((len(temperatures), len(sim_system.snames)))
    drcs = dict()
    for Tind, T in enumerate(temperatures):
        print('* %1.0f K' % T)
        sim_system.params['temperature'] = T
        run(sim_system=sim_system, plot_results=plot_results, save_results=save_results,
            fig_path=fig_path, csv_path=csv_path)
        if steady_state_solve:
            sim_system.find_steady(store_steady=True)
            final[Tind, :] = sim_system.full_steady
        else:
            final[Tind, :] = sim_system.solution[-1]
        sim_system.reaction_terms(final[Tind, :])
        rates[Tind, :] = sim_system.rates[:, 0] - sim_system.rates[:, 1]
        if tof_terms is not None:
            drcs[T] = sim_system.degree_of_rate_control(tof_terms, eps=eps)

    if plot_results:
        if fig_path:
            if not os.path.isdir(fig_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(fig_path)

        cmap = plt.get_cmap("Spectral", len(sim_system.adsorbate_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(sim_system.snames):
            if i in sim_system.adsorbate_indices:
                ax.plot(temperatures, final[:, i], label=sname,
                        color=cmap(sim_system.adsorbate_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel='Temperature (K)',
               ylabel='Coverage', ylim=(-0.1, 1.1))
        fig.tight_layout()
        if fig_path:
            plt.savefig(fig_path + 'coverages_vs_temperature.png', format='png', dpi=300)

        cmap = plt.get_cmap("Accent", len(sim_system.gas_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(sim_system.snames):
            if i in sim_system.gas_indices:
                ax.plot(temperatures, final[:, i], label=sname,
                        color=cmap(sim_system.gas_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel='Temperature (K)',
               ylabel='Pressure (bar)')
        fig.tight_layout()
        if fig_path:
            plt.savefig(fig_path + 'pressures_vs_temperature.png', format='png', dpi=300)

        cmap = plt.get_cmap("Accent", len(sim_system.reactions))
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        for i, rname in enumerate(sim_system.reactions.keys()):
            ax.plot(temperatures, rates[:, i], label=rname, color=cmap(i))
        ax.legend(loc='best', frameon=False, ncol=4)
        ax.set(xlabel='Temperature (K)',
               ylabel='Rate (1/s)', yscale='log')
        fig.tight_layout()
        if fig_path:
            plt.savefig(fig_path + 'surfrates_vs_temperature.png', format='png', dpi=300)

        if tof_terms is not None:
            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            for rind, rname in enumerate(sim_system.reactions.keys()):
                ax.plot(temperatures, [drcs[i][rname] for i in temperatures], label=rname, color=cmap(rind))
            ax.set(xlabel='Temperature (K)',
                   ylabel='Degree of rate control')
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
            if fig_path:
                plt.savefig(fig_path + 'drc_vs_temperature.png', format='png', dpi=300)

            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            ax.plot(temperatures,
                    np.sum(rates[:, [list(sim_system.reactions.keys()).index(i) for i in tof_terms]], axis=1),
                    '-+', color='k')
            ax.set(xlabel='Temperature (K)',
                   ylabel='Rate (1/s)', yscale='log')
            if fig_path:
                plt.savefig(fig_path + 'tof_vs_temperature.png', format='png', dpi=300)

    if save_results:
        if csv_path is not '':
            if not os.path.isdir(csv_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(csv_path)

        rfile = csv_path + 'rates_vs_temperature.csv'
        cfile = csv_path + 'coverages_vs_temperature.csv'
        pfile = csv_path + 'pressures__vs_temperature.csv'

        rheader = ['Temperature (K)'] + list(sim_system.reactions.keys())
        cheader = ['Temperature (K)'] + [s for i, s in enumerate(sim_system.snames)
                                         if i in sim_system.adsorbate_indices]
        pheader = ['Temperature (K)'] + [s for i, s in enumerate(sim_system.snames)
                                         if i in sim_system.gas_indices]

        df = pd.DataFrame(np.concatenate((np.reshape(temperatures, (len(temperatures), 1)),
                                          rates), axis=1), columns=rheader)
        df.to_csv(path_or_buf=rfile, sep=',', header=True, index=False)
        df = pd.DataFrame(np.concatenate((np.reshape(temperatures, (len(temperatures), 1)),
                                          final[:, sim_system.adsorbate_indices]), axis=1), columns=cheader)
        df.to_csv(path_or_buf=cfile, sep=',', header=True, index=False)
        df = pd.DataFrame(np.concatenate((np.reshape(temperatures, (len(temperatures), 1)),
                                          final[:, sim_system.gas_indices]), axis=1), columns=pheader)
        df.to_csv(path_or_buf=pfile, sep=',', header=True, index=False)

        if tof_terms is not None:
            dfile = csv_path + 'drcs_vs_temperature.csv'
            dheader = ['Temperature (K)'] + list(sim_system.reactions.keys())
            vals = np.zeros((len(temperatures), len(list(sim_system.reactions.keys())) + 1))
            vals[:, 0] = temperatures
            for Tind, T in enumerate(temperatures):
                vals[Tind, 1::] = np.array(list(drcs[T].values()))
            df = pd.DataFrame(vals, columns=dheader)
            df.to_csv(path_or_buf=dfile, sep=',', header=True, index=False)


def draw_states(sim_system, rotation='', fig_path=None):
    """Uses ASE to draw the atoms objects for each state
    and optionally saves them.

    """

    if fig_path:
        if not os.path.isdir(fig_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(fig_path)
    for s in sim_system.snames:
        sim_system.states[s].view_atoms(rotation=rotation, path=fig_path)


def draw_energy_landscapes(sim_system, etype='free', eunits='eV', legend_location='upper right',
                           show_labels=False, fig_path=None):
    """Draws the energy landscapes using the parameters
    saved to sim_system.params and optionally saves them.

    """

    if fig_path:
        if not os.path.isdir(fig_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(fig_path)
    for k in sim_system.energy_landscapes.keys():
        sim_system.energy_landscapes[k].draw_energy_landscape(T=sim_system.params['temperature'],
                                                              p=sim_system.params['pressure'],
                                                              verbose=sim_system.params['verbose'],
                                                              etype=etype, eunits=eunits,
                                                              legend_location=legend_location,
                                                              path=fig_path, show_labels=show_labels)


def run_energy_span_temperatures(sim_system, temperatures, etype='free', save_results=False, csv_path=''):
    """Runs the energy span model for a range of temperatures
    and optionally saves the results.

    """

    if save_results:
        if csv_path is not '':
            if not os.path.isdir(csv_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(csv_path)

    for k in sim_system.energy_landscapes.keys():
        print('* Landscape %s:' % k)
        esm = dict()
        for Tind, T in enumerate(temperatures):
            sim_system.params['temperature'] = T

            esm[T] = sim_system.energy_landscapes[k].evaluate_energy_span_model(T=sim_system.params['temperature'],
                                                                                p=sim_system.params['pressure'],
                                                                                verbose=sim_system.params['verbose'],
                                                                                etype=etype)
        if save_results:
            df = pd.DataFrame(data=[[T] + list(esm[T][0:4]) for T in temperatures],
                              columns=['Temperature (K)', 'TOF (1/s)', 'Espan (eV)', 'TDTS', 'TDI'])
            df.to_csv(path_or_buf=csv_path + 'energy_span_summary_%s.csv' % k, sep=',', header=True, index=False)
            df = pd.DataFrame(data=[[T] + esm[T][4] for T in temperatures],
                              columns=['Temperature (K)'] + esm[temperatures[0]][6])
            df.to_csv(path_or_buf=csv_path + 'energy_span_xTDTS_%s.csv' % k, sep=',', header=True, index=False)
            df = pd.DataFrame(data=[[T] + esm[T][5] for T in temperatures],
                              columns=['Temperature (K)'] + esm[temperatures[0]][7])
            df.to_csv(path_or_buf=csv_path + 'energy_span_xTDI_%s.csv' % k, sep=',', header=True, index=False)
