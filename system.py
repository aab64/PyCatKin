import numpy as np
from scipy.integrate import odeint
from physical_constants import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


class System:

    def __init__(self, reactions, reactor=None):
        self.reactions = reactions
        self.reactor = reactor
        self.snames = None
        self.species_map = None
        self.adsorbate_indices = None
        self.gas_indices = None
        self.rate_constants = None
        self.conditions = None
        self.rates = None
        self.solution = None
        self.times = None
        self.names_to_indices()

    def names_to_indices(self):
        """Compiles species names and assigns indicies.

        """
        self.snames = []
        for r in self.reactions.keys():
            self.snames += [i.name for i in self.reactions[r].reactants if i.name not in self.snames]
            self.snames += [i.name for i in self.reactions[r].products if i.name not in self.snames]
        self.snames = list(set(self.snames))

        self.species_map = dict()

        for r in self.reactions.keys():
            yreac = [self.snames.index(i.name) for i in self.reactions[r].reactants if i.state_type == 'adsorbate']
            preac = [self.snames.index(i.name) for i in self.reactions[r].reactants if i.state_type == 'gas']
            yprod = [self.snames.index(i.name) for i in self.reactions[r].products if i.state_type == 'adsorbate']
            pprod = [self.snames.index(i.name) for i in self.reactions[r].products if i.state_type == 'gas']
            self.species_map[r] = {'yreac': yreac, 'yprod': yprod, 'preac': preac, 'pprod': pprod,
                                   'site_density': 1.0 / self.reactions[r].area, 'scaling': self.reactions[r].scaling}
            if self.adsorbate_indices is None:
                if yreac or yprod:
                    self.adsorbate_indices = []
                    self.adsorbate_indices += yreac
                    self.adsorbate_indices += yprod
            else:
                self.adsorbate_indices += yreac
                self.adsorbate_indices += yprod
            if self.gas_indices is None:
                if preac or pprod:
                    self.gas_indices = []
                    self.gas_indices += preac
                    self.gas_indices += pprod
            else:
                self.gas_indices += preac
                self.gas_indices += pprod

        if self.adsorbate_indices is not None:
            self.adsorbate_indices = list(set(self.adsorbate_indices))
            self.reactor.is_adsorbate = [1 if i in self.adsorbate_indices else 0 for i in range(len(self.snames))]
        else:
            self.reactor.is_adsorbate = np.zeros(len(self.snames))
        if self.gas_indices is not None:
            self.gas_indices = list(set(self.gas_indices))
            self.reactor.is_gas = [1 if i in self.gas_indices else 0 for i in range(len(self.snames))]
        else:
            self.reactor.is_gas = np.zeros(len(self.snames))

    def check_rate_constants(self, T, p):
        """Check if the rate constants have been calculated
        and are updated to the current conditions.

        """
        update = True
        if self.conditions is None or self.rate_constants is None:
            self.conditions = dict()
            self.conditions['temperature'] = T
            self.conditions['pressure'] = p
            self.rate_constants = dict()
        elif self.conditions['temperature'] != T or self.conditions['pressure'] != p:
            self.conditions['temperature'] = T
            self.conditions['pressure'] = p
        else:
            update = False
        if update:
            for r in self.reactions.keys():
                self.reactions[r].calc_rate_constants(T, p)
                self.rate_constants[r] = {'kfwd': self.reactions[r].kfwd,
                                          'krev': self.reactions[r].krev}

    def reaction_terms(self, y, T, p):
        """Constructs forward and reverse reaction rate terms
        by multiplying rate constants by reactant species coverages.

        """
        self.check_rate_constants(T=T, p=p)
        self.rates = np.zeros((len(self.species_map), 2))
        yfree = 1 - sum(y[self.adsorbate_indices])
        for rind, r in enumerate(self.species_map.keys()):
            self.rates[rind, 0] = self.rate_constants[r]['kfwd']
            self.rates[rind, 1] = self.rate_constants[r]['krev']
            freebal = len(self.species_map[r]['yreac']) - len(self.species_map[r]['yprod'])
            for i in self.species_map[r]['yreac']:  # Forward rate species coverages
                self.rates[rind, 0] *= y[i]
            for i in self.species_map[r]['yprod']:  # Reverse rate species coverages
                self.rates[rind, 1] *= y[i]
            for i in self.species_map[r]['preac']:  # Forward rate species pressures
                self.rates[rind, 0] *= (y[i] * bartoPa)
            for i in self.species_map[r]['pprod']:  # Reverse rate species pressures
                self.rates[rind, 1] *= (y[i] * bartoPa)
            for i in range(0, freebal):
                self.rates[rind, 1] *= yfree
            for i in range(freebal, 0):
                self.rates[rind, 0] *= yfree

    def species_odes(self, y, T, p):
        """Constructs species ODEs for adsorbate coverage from reaction terms.

        Returns array of species ODEs."""
        # Reaction rate terms
        self.reaction_terms(y=y, T=T, p=p)
        net_rates = self.rates[:, 0] - self.rates[:, 1]

        # Construct ODE for each species
        dy_dt = np.zeros(len(y))
        for rind, rinfo in enumerate(self.species_map.values()):
            for sp in rinfo['yreac']:  # Species is consumed
                dy_dt[sp] -= net_rates[rind] * rinfo['scaling']
            for sp in rinfo['yprod']:  # Species is formed
                dy_dt[sp] += net_rates[rind] * rinfo['scaling']
            for sp in rinfo['preac']:
                dy_dt[sp] -= net_rates[rind] * rinfo['scaling'] * rinfo['site_density']
            for sp in rinfo['pprod']:
                dy_dt[sp] += net_rates[rind] * rinfo['scaling'] * rinfo['site_density']
        return dy_dt

    def run_odeint(self, times=None, start_state=None, T=None, p=None, inflow_state=None,
                   rtol=1e-6, atol=1e-8, verbose=False):
        """Wrapper for ODE integrator.

        Returns times and solution variable."""
        # Set initial coverages to zero if not specified
        yinit = np.zeros(len(self.snames))
        print(self.snames)
        print(len(self.snames))
        print('-------------------------')
        if start_state is not None:
            for s in start_state.keys():
                yinit[self.snames.index(s)] = start_state[s]

        # Set inflow mole fractions to zero if not specified
        yinflow = np.zeros(len(self.snames))
        print(self.snames)
        print(len(self.snames))
        print('-------------------------')
        if inflow_state is not None:
            for s in inflow_state.keys():
                yinflow[self.snames.index(s)] = inflow_state[s]

        if verbose:
            print('=========\nInitial conditions:\n')
            for s, sname in enumerate(self.snames):
                print('%15s : %1.2e' % (sname, yinit[s]))
            if yinflow.any():
                print('=========\nInflow conditions:\n')
                for s, sname in enumerate(self.snames):
                    if s in self.gas_indices:
                        print('%15s : %1.2e' % (sname, yinflow[s]))

        # Choose solution times if not specified
        if times is None:
            times = np.logspace(np.log10(1e-14), np.log10(1e14), num=60)

        # Integrate species ODEs for given times
        y = odeint(self.reactor.rhs(self.species_odes), yinit, times, rtol=rtol, atol=atol, args=(T, p, yinflow))

        if verbose:
            print('=========\nFinal conditions:\n')
            for s, sname in enumerate(self.snames):
                print('%15s : %9.2e' % (sname, y[-1][s]))

        self.times = times
        self.solution = y

    def write_results(self, T, p, path=''):
        """Write reaction rates, coverages and pressures to file.
        Reaction rates computed at temperature T and pressure p.
        File written to current directory unless otherwise specified

        """

        if path is not '':
            if not os.path.isdir(path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(path)

        rfile = path + 'surfrates_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.txt'
        cfile = path + 'coverages_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.txt'
        pfile = path + 'pressures_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.txt'

        with open(rfile, 'w') as file:
            header = ', '.join(['Time (s)'] + [(r.name + '_fwd, ' + r.name + '_rev')
                                               for r in self.reactions.values()]) + '\n'
            file.write(header)
        with open(cfile, 'w') as file:
            header = ', '.join(['Time (s)'] + [s for i, s in enumerate(self.snames)
                                               if i in self.adsorbate_indices]) + '\n'
            file.write(header)
        with open(pfile, 'w') as file:
            header = ', '.join(['Time (s)'] + [s for i, s in enumerate(self.snames)
                                               if i in self.gas_indices]) + '\n'
            file.write(header)

        for t in range(len(self.times)):
            with open(rfile, 'a') as file:
                self.reaction_terms(y=self.solution[t, :], T=T, p=p)
                line = ', '.join([str(self.times[t])] + [(str(r[0]) + ', ' + str(r[1])) for r in self.rates]) + '\n'
                file.write(line)
            with open(cfile, 'a') as file:
                line = ', '.join([str(self.times[t])] + [str(self.solution[t, s]) for s in range(len(self.snames))
                                                         if s in self.adsorbate_indices]) + '\n'
                file.write(line)
            with open(pfile, 'a') as file:
                line = ', '.join([str(self.times[t])] + [str(self.solution[t, s]) for s in range(len(self.snames))
                                                         if s in self.gas_indices]) + '\n'
                file.write(line)

    def plot_transient(self, T, p, path=None):
        """Plot transient rates, coverages and pressures.
        If path is specified, figures are saved to path

        """
        font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
        plt.rc('font', **font)
        mpl.rcParams['lines.markersize'] = 6
        mpl.rcParams['lines.linewidth'] = 1.5

        if path is not None:
            if not os.path.isdir(path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(path)

        rates = np.zeros((len(self.times), len(self.reactions) * 2))
        for t in range(len(self.times)):
            self.reaction_terms(y=self.solution[t, :], T=T, p=p)
            for i in range(len(self.reactions)):
                rates[t, 2 * i] = self.rates[i, 0]
                rates[t, 2 * i + 1] = self.rates[i, 1]

        fig, ax = plt.subplots(figsize=(3.33, 3.33))
        for i, sname in enumerate(self.snames):
            if i in self.adsorbate_indices:
                ax.plot(self.times, self.solution[:, i], label=sname)
        ax.legend(loc='upper center', frameon=False, ncol=2)
        ax.set(xlabel='Time (s)', ylabel='Coverage', title=('%1.1f K' % T))
        fig.tight_layout()
        if path is not None:
            figname = path + 'coverages_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=300)

        fig, ax = plt.subplots(figsize=(3.33, 3.33))
        for i, sname in enumerate(self.snames):
            if i in self.gas_indices:
                ax.plot(self.times, self.solution[:, i], label=sname)
        ax.legend(loc='upper center', frameon=False, ncol=2)
        ax.set(xlabel='Time (s)', ylabel='Pressure (bar)', title=('%1.1f K' % T))
        fig.tight_layout()
        if path is not None:
            figname = path + 'pressures_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=300)

        fig, ax = plt.subplots(figsize=(3.33, 3.33))
        for i, rname in enumerate([r for rname in self.reactions.keys() for r in [rname + '_fwd', rname + '_rev']]):
            ax.plot(self.times, rates[:, i], label=rname)
        ax.legend(loc='upper center', frameon=False, ncol=2)
        ax.set(xlabel='Time (s)', ylabel='Rate (1/s)', title=('%1.1f K' % T))
        fig.tight_layout()
        if path is not None:
            figname = path + 'surfrates_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=300)
