from pycatkin.constants.physical_constants import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
import os

font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['lines.linewidth'] = 1.5


def run(sim_system, steady_state_solve=False, plot_results=False, save_results=False,
        fig_path=None, csv_path=''):
    """Runs the ODE solver and optionally plots/saves
    the results.

    """
    sim_system.solve_odes()
    if plot_results:
        sim_system.plot_transient(path=fig_path)
    if save_results:
        sim_system.write_results(path=csv_path)
    if steady_state_solve:
        sim_system.find_steady(store_steady=True)


def run_temperatures(sim_system, temperatures, steady_state_solve=False, tof_terms=None, eps=5.0e-2,
                     plot_results=False, save_results=False, plot_transient=False, save_transient=False,
                     fig_path=None, csv_path=''):
    """Runs the ODE solver for a range of temperatures
    and optionally plots/saves the results.

    """

    rates = np.zeros((len(temperatures), len(sim_system.reactions)))
    final = np.zeros((len(temperatures), len(sim_system.snames)))
    drcs = dict()
    print('Running simulations for T in [%1.1f K, %1.1f K]...' % (temperatures[0], temperatures[-1]))
    for Tind, T in enumerate(temperatures):
        sim_system.params['temperature'] = T
        run(sim_system=sim_system, plot_results=plot_transient, save_results=save_transient,
            fig_path=fig_path, csv_path=csv_path)
        final_time = sim_system.params['times'][-1]
        if steady_state_solve:
            # while (max(abs(sim_system.species_odes(sim_system.solution[-1])[sim_system.dynamic_indices])) > 1.0e-10
            #        and sim_system.params['times'][-1] < 1.0e11):
            #     sim_system.params['times'][-1] = sim_system.params['times'][-1] ** 2.0
            #     print('System not steady, increasing final time to %1.2e s' % sim_system.params['times'][-1])
            #     run(sim_system=sim_system)
            sim_system.find_steady(store_steady=True)
            final[Tind, :] = sim_system.full_steady
            # final[Tind, :] = sim_system.solution[-1]
            sim_system.params['times'][-1] = final_time
        else:
            final[Tind, :] = sim_system.solution[-1]
        sim_system.reaction_terms(final[Tind, :])
        rates[Tind, :] = sim_system.rates[:, 0] - sim_system.rates[:, 1]
        if tof_terms is not None:
            drcs[T] = sim_system.degree_of_rate_control(tof_terms, eps=eps)
        print('* %1.0f K done' % T)

    if plot_results:
        if fig_path is not None and fig_path != '':
            if not os.path.isdir(fig_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(fig_path)

        cmap = plt.get_cmap("tab20", len(sim_system.adsorbate_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(sim_system.snames):
            if i in sim_system.adsorbate_indices and max(final[:, i]) > 0.01:
                ax.plot(temperatures, final[:, i], label=sname,
                        color=cmap(sim_system.adsorbate_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel='Temperature (K)',
               ylabel='Coverage', ylim=(-0.1, 1.1))
        fig.tight_layout()
        if fig_path is not None:
            plt.savefig(fig_path + 'coverages_vs_temperature.png', format='png', dpi=600)

        cmap = plt.get_cmap("tab20", len(sim_system.gas_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(sim_system.snames):
            if i in sim_system.gas_indices:
                ax.plot(temperatures, final[:, i], label=sname,
                        color=cmap(sim_system.gas_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel='Temperature (K)',
               ylabel='Pressure (bar)')
        fig.tight_layout()
        if fig_path is not None:
            plt.savefig(fig_path + 'pressures_vs_temperature.png', format='png', dpi=600)

        cmap = plt.get_cmap("tab20", len(sim_system.reactions))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, rname in enumerate(sim_system.reactions.keys()):
            ax.plot(temperatures, rates[:, i], label=rname, color=cmap(i))
        ax.legend(loc='best', frameon=False, ncol=1)
        yvals = ax.get_ylim()
        ax.set(xlabel='Temperature (K)',
               ylabel='Rate (1/s)', yscale='log', ylim=(max(1e-10, yvals[0]), yvals[1]))
        fig.tight_layout()
        if fig_path is not None:
            plt.savefig(fig_path + 'surfrates_vs_temperature.png', format='png', dpi=600)

        if tof_terms is not None:
            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            for rind, rname in enumerate(sim_system.reactions.keys()):
                drc = [drcs[i][rname] for i in temperatures]
                if max([abs(d) for d in drc]) > 0.01:
                    ax.plot(temperatures, drc, label=rname, color=cmap(rind))
            ax.set(xlabel='Temperature (K)',
                   ylabel='Degree of rate control')
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
            if fig_path is not None:
                plt.savefig(fig_path + 'drc_vs_temperature.png', format='png', dpi=600)

            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            ax.plot(temperatures,
                    np.sum(rates[:, [list(sim_system.reactions.keys()).index(i) for i in tof_terms]], axis=1),
                    color='k')
            ax.set(xlabel='Temperature (K)',
                   ylabel='TOF (1/s)', yscale='log')
            fig.tight_layout()
            if fig_path is not None:
                plt.savefig(fig_path + 'tof_vs_temperature.png', format='png', dpi=600)

    if save_results:
        if csv_path != '':
            if not os.path.isdir(csv_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(csv_path)

        rfile = csv_path + 'rates_vs_temperature.csv'
        cfile = csv_path + 'coverages_vs_temperature.csv'
        pfile = csv_path + 'pressures_vs_temperature.csv'

        rheader = ['Temperature (K)'] + list(sim_system.reactions.keys())
        cheader = ['Temperature (K)'] + [s for i, s in enumerate(sim_system.snames)
                                         if i in sim_system.adsorbate_indices]
        pheader = ['Temperature (K)'] + ['p' + s + ' (bar)' for i, s in enumerate(sim_system.snames)
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


def run_parameters(sim_system, parameters, params_name, steady_state_solve=False, tof_terms=None, eps=5.0e-2,
                   plot_results=False, save_results=False, plot_transient=False, save_transient=False,
                   fig_path=None, csv_path=''):
    """Runs the ODE solver for a range of parameter values
    and optionally plots/saves the results.

    """

    rates = np.zeros((len(parameters), len(sim_system.reactions)))
    final = np.zeros((len(parameters), len(sim_system.snames)))
    drcs = dict()
    print('Running simulations for %s in [%1.3f, %1.3f]...' % (params_name, parameters[0], parameters[-1]))
    for pind, param in enumerate(parameters):
        if 'start_state' in params_name:
            sim_system.params['start_state'][params_name.split('start_state_')[1]] = param
        elif 'inflow_state' in params_name:
            sim_system.params['inflow_state'][params_name.split('inflow_state_')[1]] = param
        else:
            sim_system.params[params_name] = param
        run(sim_system=sim_system, plot_results=plot_transient, save_results=save_transient,
            fig_path=fig_path, csv_path=csv_path)
        final_time = sim_system.params['times'][-1]
        if steady_state_solve:
            sim_system.find_steady(store_steady=True)
            final[pind, :] = sim_system.full_steady
            sim_system.params['times'][-1] = final_time
        else:
            final[pind, :] = sim_system.solution[-1]
        sim_system.reaction_terms(final[pind, :])
        rates[pind, :] = sim_system.rates[:, 0] - sim_system.rates[:, 1]
        if tof_terms is not None:
            drcs[param] = sim_system.degree_of_rate_control(tof_terms, eps=eps)
        print('* %1.3f done' % param)

    if plot_results:
        if fig_path is not None and fig_path != '':
            if not os.path.isdir(fig_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(fig_path)

        cmap = plt.get_cmap("tab20", len(sim_system.adsorbate_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(sim_system.snames):
            if i in sim_system.adsorbate_indices and max(final[:, i]) > 0.01:
                ax.plot(parameters, final[:, i], label=sname,
                        color=cmap(sim_system.adsorbate_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel=params_name,
               ylabel='Coverage', ylim=(-0.1, 1.1))
        fig.tight_layout()
        if fig_path is not None:
            plt.savefig(fig_path + 'coverages_vs_%s.png' % params_name, format='png', dpi=600)

        cmap = plt.get_cmap("tab20", len(sim_system.gas_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(sim_system.snames):
            if i in sim_system.gas_indices:
                ax.plot(parameters, final[:, i], label=sname,
                        color=cmap(sim_system.gas_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel=params_name,
               ylabel='Pressure (bar)')
        fig.tight_layout()
        if fig_path is not None:
            plt.savefig(fig_path + 'pressures_vs_%s.png' % params_name, format='png', dpi=600)

        cmap = plt.get_cmap("tab20", len(sim_system.reactions))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, rname in enumerate(sim_system.reactions.keys()):
            ax.plot(parameters, rates[:, i], label=rname, color=cmap(i))
        ax.legend(loc='best', frameon=False, ncol=1)
        yvals = ax.get_ylim()
        ax.set(xlabel=params_name,
               ylabel='Rate (1/s)', yscale='log', ylim=(max(1e-10, yvals[0]), yvals[1]))
        fig.tight_layout()
        if fig_path is not None:
            plt.savefig(fig_path + 'surfrates_vs_%s.png' % params_name, format='png', dpi=600)

        if tof_terms is not None:
            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            for rind, rname in enumerate(sim_system.reactions.keys()):
                drc = [drcs[i][rname] for i in parameters]
                if max([abs(d) for d in drc]) > 0.01:
                    ax.plot(parameters, drc, label=rname, color=cmap(rind))
            ax.set(xlabel=params_name,
                   ylabel='Degree of rate control')
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
            if fig_path is not None:
                plt.savefig(fig_path + 'drc_vs_%s.png' % params_name, format='png', dpi=600)

            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            ax.plot(parameters,
                    np.sum(rates[:, [list(sim_system.reactions.keys()).index(i) for i in tof_terms]], axis=1),
                    color='k')
            ax.set(xlabel=params_name,
                   ylabel='TOF (1/s)', yscale='log')
            fig.tight_layout()
            if fig_path is not None:
                plt.savefig(fig_path + 'tof_vs_%s.png' % params_name, format='png', dpi=600)

    if save_results:
        if csv_path != '':
            if not os.path.isdir(csv_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(csv_path)

        rfile = csv_path + 'rates_vs_%s.csv' % params_name
        cfile = csv_path + 'coverages_vs_%s.csv' % params_name
        pfile = csv_path + 'pressures_vs_%s.csv' % params_name

        rheader = [params_name] + list(sim_system.reactions.keys())
        cheader = [params_name] + [s for i, s in enumerate(sim_system.snames)
                                   if i in sim_system.adsorbate_indices]
        pheader = [params_name] + ['p' + s + ' (bar)' for i, s in enumerate(sim_system.snames)
                                   if i in sim_system.gas_indices]

        df = pd.DataFrame(np.concatenate((np.reshape(parameters, (len(parameters), 1)),
                                          rates), axis=1), columns=rheader)
        df.to_csv(path_or_buf=rfile, sep=',', header=True, index=False)
        df = pd.DataFrame(np.concatenate((np.reshape(parameters, (len(parameters), 1)),
                                          final[:, sim_system.adsorbate_indices]), axis=1), columns=cheader)
        df.to_csv(path_or_buf=cfile, sep=',', header=True, index=False)
        df = pd.DataFrame(np.concatenate((np.reshape(parameters, (len(parameters), 1)),
                                          final[:, sim_system.gas_indices]), axis=1), columns=pheader)
        df.to_csv(path_or_buf=pfile, sep=',', header=True, index=False)

        if tof_terms is not None:
            dfile = csv_path + 'drcs_vs_%s.csv' % params_name
            dheader = [params_name] + list(sim_system.reactions.keys())
            vals = np.zeros((len(parameters), len(list(sim_system.reactions.keys())) + 1))
            vals[:, 0] = parameters
            for pind, param in enumerate(parameters):
                vals[pind, 1::] = np.array(list(drcs[param].values()))
            df = pd.DataFrame(vals, columns=dheader)
            df.to_csv(path_or_buf=dfile, sep=',', header=True, index=False)


def draw_states(sim_system, rotation='', fig_path=None):
    """Uses ASE to draw the atoms objects for each state
    and optionally saves them.

    """

    if fig_path is not None and fig_path != '':
        if not os.path.isdir(fig_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(fig_path)
    for s in sim_system.snames:
        if not isinstance(sim_system.states[s], Scaling):
            sim_system.states[s].view_atoms(rotation=rotation, path=fig_path)


def draw_energy_landscapes(sim_system, etype='free', eunits='eV', legend_location='upper right',
                           show_labels=False, fig_path=None):
    """Draws the energy landscapes using the parameters
    saved to sim_system.params and optionally saves them.

    """

    if fig_path is not None and fig_path != '':
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
        if csv_path != '':
            if not os.path.isdir(csv_path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(csv_path)

    for k in sim_system.energy_landscapes.keys():
        print('Landscape %s:' % k)
        print('-----------------')
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


def save_energies(sim_system, csv_path=''):
    """Save the reaction energies for the current temperature.

    """

    if csv_path != '':
        if not os.path.isdir(csv_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(csv_path)

    T = sim_system.params['temperature']
    p = sim_system.params['pressure']
    v = sim_system.params['verbose']

    evals = dict()
    print('Saving reaction energies...')
    for r in sim_system.reactions.keys():

        evals[r] = [sim_system.reactions[r].get_reaction_energy(T=T, p=p, verbose=v, etype='electronic')]
        evals[r].append(sim_system.reactions[r].get_reaction_energy(T=T, p=p, verbose=v, etype='free'))
        evals[r].append(sim_system.reactions[r].get_reaction_barriers(T=T, p=p, verbose=v, etype='electronic')[0])
        evals[r].append(sim_system.reactions[r].get_reaction_barriers(T=T, p=p, verbose=v, etype='free')[0])

        print('* Reaction %s done' % r)

    df = pd.DataFrame(data=[[r] + evals[r] for r in sim_system.reactions.keys()],
                      columns=['Reaction', 'dEr (J/mol)', 'dGr (J/mol)', 'dEa (J/mol)', 'dGa (J/mol)'])
    df.to_csv(path_or_buf=csv_path + 'reaction_energies_and_barriers_%1.1fK_%1.1fbar.csv' % (T, p / bartoPa),
              sep=',', header=True, index=False)


def save_energies_temperatures(sim_system, temperatures, csv_path=''):
    """Save the reaction energies for a range of temperatures.

    """

    if csv_path != '':
        if not os.path.isdir(csv_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(csv_path)

    p = sim_system.params['pressure']
    v = sim_system.params['verbose']

    print('Saving reaction energies...')
    for r in sim_system.reactions.keys():
        evals = dict()
        for Tind, T in enumerate(temperatures):
            sim_system.params['temperature'] = T

            evals[T] = [T]
            evals[T].append(sim_system.reactions[r].get_reaction_energy(T=T, p=p, verbose=v, etype='electronic'))
            evals[T].append(sim_system.reactions[r].get_reaction_energy(T=T, p=p, verbose=v, etype='free'))
            evals[T].append(sim_system.reactions[r].get_reaction_barriers(T=T, p=p, verbose=v, etype='electronic')[0])
            evals[T].append(sim_system.reactions[r].get_reaction_barriers(T=T, p=p, verbose=v, etype='free')[0])

        df = pd.DataFrame(data=[evals[T] for T in temperatures],
                          columns=['Temperature (K)', 'dEr (J/mol)', 'dGr (J/mol)', 'dEa (J/mol)', 'dGa (J/mol)'])
        df.to_csv(path_or_buf=csv_path + 'reaction_energies_and_barriers_%s.csv' % r,
                  sep=',', header=True, index=False)
        print('* Reaction %s done' % r)


def save_state_energies(sim_system, csv_path=''):
    """Save the state energies for the current temperature.

    """

    if csv_path != '':
        if not os.path.isdir(csv_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(csv_path)

    T = sim_system.params['temperature']
    p = sim_system.params['pressure']
    v = sim_system.params['verbose']

    evals = dict()
    print('Saving state energies...')
    for s in sim_system.snames:

        evals[s] = [sim_system.states[s].get_free_energy(T=T, p=p, verbose=v)]
        evals[s] += [sim_system.states[s].Gelec if sim_system.states[s].Gelec is not None else None]
        evals[s] += [sim_system.states[s].Gvibr if sim_system.states[s].Gvibr is not None else None]
        evals[s] += [sim_system.states[s].Grota if sim_system.states[s].Grota is not None else None]
        evals[s] += [sim_system.states[s].Gtran if sim_system.states[s].Gtran is not None else None]

        print('* State %s done' % s)

    df = pd.DataFrame(data=[[s] + evals[s] for s in sim_system.snames],
                      columns=['State', 'Free (eV)', 'Electronic (eV)', 'Vibrational (eV)', 'Translational (eV)',
                               'Rotational (eV)'])
    df.to_csv(path_or_buf=csv_path + 'state_energies_%1.1fK_%1.1fbar.csv' % (T, p / bartoPa),
              sep=',', header=True, index=False)


def save_pes_energies(sim_system, csv_path=''):
    """Save the state energies for the current temperature.

    """

    if csv_path != '':
        if not os.path.isdir(csv_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(csv_path)

    T = sim_system.params['temperature']
    p = sim_system.params['pressure']
    v = sim_system.params['verbose']

    evals = dict()
    print('Saving state energies...')
    for k in sim_system.energy_landscapes.keys():
        sim_system.energy_landscapes[k].construct_energy_landscape(T=T, p=p, verbose=v)
        for s in sim_system.energy_landscapes[k].energy_landscape['free'].keys():
            evals[s] = [sim_system.energy_landscapes[k].energy_landscape['free'][s],
                        sim_system.energy_landscapes[k].energy_landscape['electronic'][s]]
        df = pd.DataFrame(data=[[sim_system.energy_landscapes[k].labels[s]] + evals[s] for s in evals.keys()],
                          columns=['State', 'Free (eV)', 'Electronic (eV)'])
        df.to_csv(path_or_buf=csv_path + str(k) + '_energy_landscape_%1.1fK_%1.1fbar.csv' % (T, p / bartoPa),
                  sep=',', header=True, index=False)


def compare_energy_landscapes(sim_systems, etype='free', eunits='eV', legend_location=None,
                              show_labels=False, fig_path=None, cmap=None):
    """Draws the energy landscapes for a number of systems using the parameters
    saved to sim_system.params and optionally saves them.

    """

    if fig_path is not None and fig_path != '':
        if not os.path.isdir(fig_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(fig_path)

    fig, ax = plt.subplots(figsize=(10, 4))
    if cmap is None:
        cmap = plt.get_cmap("tab20", len(sim_systems))

    for sind, sim_system in enumerate(sim_systems.values()):
        for k in sim_system.energy_landscapes.keys():
            fig, ax = sim_system.energy_landscapes[k].draw_energy_landscape_simple(T=sim_system.params['temperature'],
                                                                                   p=sim_system.params['pressure'],
                                                                                   verbose=sim_system.params['verbose'],
                                                                                   fig=fig, ax=ax,
                                                                                   linecolor=cmap(sind),
                                                                                   etype=etype, eunits=eunits,
                                                                                   show_labels=show_labels)
    if legend_location is not None:
        yvals = ax.get_ylim()
        xvals = ax.get_xlim()

        for sind, sname in enumerate(sim_systems.keys()):
            ax.plot(xvals, (yvals[0] - 1e6, yvals[0] - 1e6), color=cmap(sind), label=sname)
        ax.set(xlim=xvals, ylim=(yvals[0] - 0.05 * abs(yvals[0]), yvals[1] + 0.05 * abs(yvals[1])))
        ax.legend(loc=legend_location)

    if fig_path is not None:
        fig.savefig(fig_path + etype + '_energy_landscapes.png', format='png', dpi=600)


def plot_data_simple(fig=None, ax=None, xdata=None, ydata=None, label=None,
                     linestyle='-', color='k',
                     xlabel=None, ylabel=None, title=None, addlegend=False, legendloc='best',
                     fig_path=None, fig_name='figure'):
    """Generic function to plot a set of data and save the figure.

     """

    if fig_path is not None and fig_path != '':
        if not os.path.isdir(fig_path):
            print('Directory does not exist. Will try creating it...')
            os.mkdir(fig_path)
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
    ax.plot(xdata, ydata, linestyle, color=color, label=label)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    if addlegend:
        ax.legend(loc=legendloc, frameon=False)
    fig.tight_layout()

    if fig_path is not None:
        fig.savefig(fig_path + fig_name + '.png', format='png', dpi=600)

    return fig, ax

