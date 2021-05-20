import numpy as np
from scipy.integrate import ode, solve_ivp
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import copy
from scipy.optimize import fsolve
from microkinetics.classes.reactor import *


class System:

    def __init__(self, reactions, reactor=None):
        self.reactions = reactions
        self.reactor = reactor
        self.snames = None
        self.species_map = None
        self.adsorbate_indices = None
        self.gas_indices = None
        self.dynamic_indices = None
        self.rate_constants = None
        self.conditions = None
        self.rates = None
        self.times = None
        self.solution = None
        self.full_steady = None
        self.params = None
        self.names_to_indices()
        self.set_parameters()

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
            yreac = [self.snames.index(i.name) for i in self.reactions[r].reactants if i.state_type == 'adsorbate' or
                     i.state_type == 'surface']
            preac = [self.snames.index(i.name) for i in self.reactions[r].reactants if i.state_type == 'gas']
            yprod = [self.snames.index(i.name) for i in self.reactions[r].products if i.state_type == 'adsorbate' or
                     i.state_type == 'surface']
            pprod = [self.snames.index(i.name) for i in self.reactions[r].products if i.state_type == 'gas']
            self.species_map[r] = {'yreac': yreac, 'yprod': yprod, 'preac': preac, 'pprod': pprod,
                                   'site_density': 1.0 / self.reactions[r].area, 'scaling': self.reactions[r].scaling,
                                   'perturbation': 0.0}
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
        self.dynamic_indices = self.reactor.get_dynamic_indices(self.adsorbate_indices, self.gas_indices)

    def set_parameters(self, times=None, start_state=None, inflow_state=None, T=293.15, p=101325.0,
                       rtol=1e-8, atol=1e-10, xtol=1e-10, use_jacobian=True, verbose=False):
        """Store simulation conditions, solver tolerances and verbosity.

        """
        self.params = dict()
        self.params['times'] = copy.copy(times)
        self.params['start_state'] = copy.copy(start_state)
        self.params['inflow_state'] = copy.copy(inflow_state)
        self.params['temperature'] = T
        self.params['pressure'] = p
        self.params['rtol'] = rtol
        self.params['atol'] = atol
        self.params['xtol'] = xtol
        self.params['jacobian'] = use_jacobian
        self.params['verbose'] = verbose

    def check_rate_constants(self):
        """Check if the rate constants have been calculated
        and are updated to the current conditions.

        """
        update = True
        if self.conditions is None or self.rate_constants is None:
            self.conditions = dict()
            self.conditions['temperature'] = self.params['temperature']
            self.conditions['pressure'] = self.params['pressure']
            self.rate_constants = dict()
        elif (self.conditions['temperature'] != self.params['temperature']) or (self.conditions['pressure'] !=
                                                                                self.params['pressure']):
            self.conditions['temperature'] = self.params['temperature']
            self.conditions['pressure'] = self.params['pressure']
        else:
            update = False
        if update:
            for r in self.reactions.keys():
                self.reactions[r].calc_rate_constants(T=self.params['temperature'], p=self.params['pressure'],
                                                      verbose=self.params['verbose'])
                self.rate_constants[r] = {'kfwd': self.reactions[r].kfwd,
                                          'krev': self.reactions[r].krev}

    def reaction_terms(self, y):
        """Constructs forward and reverse reaction rate terms
        by multiplying rate constants by reactant coverages/pressures.

        """
        self.check_rate_constants()
        self.rates = np.zeros((len(self.species_map), 2))

        ny = max(y.shape)
        y = y.reshape((ny, 1))

        # yfree = 1 - sum(y[self.adsorbate_indices, 0])

        for rind, r in enumerate(self.species_map.keys()):
            self.rates[rind, 0] = self.rate_constants[r]['kfwd'] + self.species_map[r]['perturbation']
            self.rates[rind, 1] = self.rate_constants[r]['krev'] * (1.0 + self.species_map[r]['perturbation'] /
                                                                    self.rate_constants[r]['kfwd'])
            # freebal = len(self.species_map[r]['yreac']) - len(self.species_map[r]['yprod'])
            for i in self.species_map[r]['yreac']:  # Forward rate species coverages
                self.rates[rind, 0] *= y[i, 0]
            for i in self.species_map[r]['yprod']:  # Reverse rate species coverages
                self.rates[rind, 1] *= y[i, 0]
            for i in self.species_map[r]['preac']:  # Forward rate species pressures
                self.rates[rind, 0] *= (y[i, 0] * bartoPa)
            for i in self.species_map[r]['pprod']:  # Reverse rate species pressures
                self.rates[rind, 1] *= (y[i, 0] * bartoPa)
            # for i in range(0, freebal):
            #     self.rates[rind, 1] *= yfree
            # for i in range(freebal, 0):
            #     self.rates[rind, 0] *= yfree

    def species_odes(self, y):
        """Constructs ODEs for adsorbate coverages from reaction rates.

        Returns array of species ODEs."""
        # Reaction rate terms
        self.reaction_terms(y=y)
        net_rates = self.rates[:, 0] - self.rates[:, 1]

        ny = max(y.shape)

        # Construct ODE for each species
        dy_dt = np.zeros(ny)
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

    def reaction_derivatives(self, y):
        """Constructs derivative of reactions wrt each species
        by multiplying rate constants by reactant coverages/pressures.

        returns ndarray of derivatives
        """
        self.check_rate_constants()

        ny = max(y.shape)
        y = y.reshape((ny, 1))
        dr_dtheta = np.zeros((len(self.species_map), ny))

        # yfree = 1 - sum(y[self.adsorbate_indices, 0])
        # assert(abs(yfree) <= 1e-16)  # Implementation assumes free sites are a state

        for rind, r in enumerate(self.species_map.keys()):
            kfwd = self.rate_constants[r]['kfwd'] + self.species_map[r]['perturbation']
            krev = self.rate_constants[r]['krev'] * (1.0 + self.species_map[r]['perturbation'] /
                                                     self.rate_constants[r]['kfwd'])

            yfwd = np.prod([y[i, 0] for i in self.species_map[r]['yreac']])
            yrev = np.prod([y[i, 0] for i in self.species_map[r]['yprod']])
            pfwd = np.prod([y[i, 0] * bartoPa for i in self.species_map[r]['preac']])
            prev = np.prod([y[i, 0] * bartoPa for i in self.species_map[r]['pprod']])
            for ind, i in enumerate(self.species_map[r]['yreac']):  # Forward rate species coverages
                dr_dtheta[rind, i] += kfwd * pfwd * np.prod([y[self.species_map[r]['yreac'][j], 0] for j in
                                                             range(len(self.species_map[r]['yreac'])) if j != ind])
            for ind, i in enumerate(self.species_map[r]['yprod']):  # Reverse rate species coverages
                dr_dtheta[rind, i] -= krev * prev * np.prod([y[self.species_map[r]['yprod'][j], 0] for j in
                                                             range(len(self.species_map[r]['yprod'])) if j != ind])
            for ind, i in enumerate(self.species_map[r]['preac']):  # Forward rate species pressures
                dr_dtheta[rind, i] += kfwd * yfwd * np.prod([y[self.species_map[r]['preac'][j], 0] * bartoPa for j in
                                                             range(len(self.species_map[r]['preac'])) if j != ind])
            for ind, i in enumerate(self.species_map[r]['pprod']):  # Reverse rate species pressures
                dr_dtheta[rind, i] -= krev * yrev * np.prod([y[self.species_map[r]['pprod'][j], 0] * bartoPa for j in
                                                             range(len(self.species_map[r]['pprod'])) if j != ind])
        return dr_dtheta

    def species_jacobian(self, y):
        """Constructs derivatives of species ODEs for adsorbate coverages.

        Returns Jacobian with shape (ny x ny)."""
        # Reaction rate derivatives
        dr_dtheta = self.reaction_derivatives(y=y)

        ny = max(y.shape)

        # Construct ODE for each species
        jac = np.zeros((ny, ny))
        for rind, rinfo in enumerate(self.species_map.values()):
            for sp1 in range(ny):
                for sp2 in rinfo['yreac']:  # Species is consumed
                    jac[sp2, sp1] -= dr_dtheta[rind, sp1] * rinfo['scaling']
                for sp2 in rinfo['yprod']:  # Species is formed
                    jac[sp2, sp1] += dr_dtheta[rind, sp1] * rinfo['scaling']
                for sp2 in rinfo['preac']:
                    jac[sp2, sp1] -= dr_dtheta[rind, sp1] * rinfo['scaling'] * rinfo['site_density']
                for sp2 in rinfo['pprod']:
                    jac[sp2, sp1] += dr_dtheta[rind, sp1] * rinfo['scaling'] * rinfo['site_density']
        return jac

    def solve_odes(self):
        """Wrapper for ODE integrator.

        Returns times and solution variable."""
        self.conditions = None  # Force rate constants to be recalculated

        # Set initial coverages to zero if not specified
        yinit = np.zeros(len(self.snames))
        if self.params['start_state'] is not None:
            for s in self.params['start_state'].keys():
                yinit[self.snames.index(s)] = self.params['start_state'][s]

        # Set inflow mole fractions to zero if not specified
        yinflow = np.zeros(len(self.snames))
        if self.params['inflow_state'] is not None:
            for s in self.params['inflow_state'].keys():
                yinflow[self.snames.index(s)] = self.params['inflow_state'][s]

        if self.params['verbose']:
            print('=========\nInitial conditions:\n')
            for s, sname in enumerate(self.snames):
                print('%15s : %1.2e' % (sname, yinit[s]))
            if yinflow.any():
                print('=========\nInflow conditions:\n')
                for s, sname in enumerate(self.snames):
                    if s in self.gas_indices:
                        print('%15s : %1.2e' % (sname, yinflow[s]))

        # Choose solution times if not specified
        times = self.params['times']
        if times is None:
            times = np.logspace(np.log10(1e-14), np.log10(1e3), num=1000)

        # sol = solve_ivp(fun=self.reactor.rhs(self.species_odes),
        #                 t_span=(self.params['times'][0], self.params['times'][-1]), y0=yinit, method='LSODA',
        #                 args=(self.params['temperature'], yinflow),
        #                 rtol=self.params['rtol'], atol=self.params['atol'])
        # if self.params['verbose']:
        #     print(sol.message)
        # self.times = sol.t
        # self.solution = np.transpose(sol.y)

        y = np.zeros((len(times) + 1, len(yinit)))
        y[0, :] = yinit

        if self.params['jacobian']:
            jacfun = self.reactor.jacobian(self.species_jacobian)
        else:
            jacfun = None

        r = ode(self.reactor.rhs(self.species_odes), jacfun).set_integrator('lsoda', method='bdf',
                                                                            rtol=self.params['rtol'],
                                                                            atol=self.params['atol'])
        r.set_initial_value(yinit, 0.0).set_f_params(self.params['temperature'], yinflow)

        if self.params['jacobian']:
            r.set_jac_params(self.params['temperature'])

        for i, t in enumerate(times):
            if r.successful() and r.t < times[-1]:
                r.integrate(t)
                y[i + 1, :] = r.y

        self.times = np.concatenate((np.zeros(1), times))
        self.solution = y

        if self.params['verbose']:
            print('=========\nFinal conditions:\n')
            for s, sname in enumerate(self.snames):
                print('%15s : %9.2e' % (sname, self.solution[-1][s]))

    def find_steady(self, store_steady=False, plot_comparison=False, path=None):
        """Solve for the steady state solution

        Returns steady state solution."""

        self.conditions = None  # Force rate constants to be recalculated

        if self.solution is not None:
            y_guess = copy.copy(self.solution[-1, self.dynamic_indices])
            full_steady = copy.copy(self.solution[-1, :])
        else:
            y_guess = np.zeros(len(self.dynamic_indices))
            full_steady = np.zeros(len(self.adsorbate_indices) + len(self.gas_indices))

        # Set inflow mole fractions to zero if not specified
        yinflow = None
        if self.params['inflow_state']:
            yinflow = np.zeros(len(self.snames))
            for s in self.params['inflow_state'].keys():
                yinflow[self.snames.index(s)] = self.params['inflow_state'][s]

        def func(y):
            full_steady[self.dynamic_indices] = y
            return self.reactor.rhs(self.species_odes)(t=0, y=full_steady, T=self.params['temperature'],
                                                       inflow_state=yinflow)[self.dynamic_indices]

        if self.params['jacobian']:
            def jacfun(y):
                full_steady[self.dynamic_indices] = y
                full_jacobian = self.reactor.jacobian(self.species_jacobian)(t=0, y=full_steady,
                                                                             T=self.params['temperature'])
                return np.array([[full_jacobian[i1, i2] for i1 in self.dynamic_indices]
                                 for i2 in self.dynamic_indices])
        else:
            jacfun = None

        y_steady, info, ier, mesg = fsolve(func=func, x0=y_guess, fprime=jacfun, full_output=True,
                                           xtol=self.params['xtol'])
        full_steady[self.dynamic_indices] = y_steady

        if store_steady:
            self.full_steady = full_steady

        if self.params['verbose']:
            print('Results of steady state search...')
            print('- At %1.0f K: %s, %1i' % (self.params['temperature'], mesg, ier))
            print('- Max abs function value at steady state: %1.3e' % np.max(np.abs(func(y_steady))))
            print('- Max abs difference: %1.3e' % (np.max(np.abs(y_steady - y_guess))))

        if plot_comparison:
            font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
            plt.rc('font', **font)
            mpl.rcParams['lines.markersize'] = 6
            mpl.rcParams['lines.linewidth'] = 1.5
            cmap = plt.get_cmap("Spectral", len(self.dynamic_indices))
            fig, ax = plt.subplots(figsize=(3.2, 3.2))
            for i in self.dynamic_indices:
                if np.max(self.solution[:, i]) > 1.0e-6:
                    ax.plot(self.times, self.solution[:, i], label=self.snames[i],
                            color=cmap(self.dynamic_indices.index(i)))
                    ax.plot(self.times, [full_steady[i] for t in self.times], label='', linestyle=':',
                            color=cmap(self.dynamic_indices.index(i)))
            ax.legend(frameon=False, loc='center right')
            ax.set(xlabel='Time (s)', ylabel='Coverage', title=(r'$T=%1.0f$ K' % self.params['temperature']),
                   ylim=(1e-6, 1e1), xscale='log', yscale='log')
            fig.tight_layout()
            if path:
                fig.savefig((path + 'SS_vs_transience_%1.1fK.png') % self.params['temperature'], format='png', dpi=300)

        return full_steady

    def run_and_return_tof(self, tof_terms, ss_solve=False):
        """Integrate or solve for the steady state and compute the TOF by summing steps in tof_terms

        Returns array of xi_i terms for each step i."""

        if ss_solve:
            full_steady = self.find_steady()
        else:
            self.solve_odes()
            full_steady = self.solution[-1, :]

        self.reaction_terms(full_steady)

        tof = 0.0
        for rind, r in enumerate(self.species_map.keys()):
            if r in tof_terms:
                tof += self.rates[rind, 0] - self.rates[rind, 1]  # assumes forward directions
        return tof

    def degree_of_rate_control(self, tof_terms, ss_solve=False, eps=1.0e-3):
        """Calculate the degree of rate control xi_i

        Returns array of xi_i terms for each step i."""

        self.conditions = None  # Force rate constants to be recalculated

        r0 = self.run_and_return_tof(tof_terms=tof_terms, ss_solve=ss_solve)
        xi = dict()
        if self.params['verbose']:
            print('Checking degree of rate control...')
        for r in self.reactions.keys():
            self.species_map[r]['perturbation'] = eps * self.rate_constants[r]['kfwd']
            xi_r = self.run_and_return_tof(tof_terms=tof_terms, ss_solve=ss_solve)
            self.species_map[r]['perturbation'] = -eps * self.rate_constants[r]['kfwd']
            xi_r -= self.run_and_return_tof(tof_terms=tof_terms, ss_solve=ss_solve)
            xi_r *= (self.rate_constants[r]['kfwd']) / (2.0 * eps * self.rate_constants[r]['kfwd'] * r0)
            xi[r] = xi_r
            self.species_map[r]['perturbation'] = 0.0
            if self.params['verbose']:
                print(r + ': done.')
        return xi

    def write_results(self, path=''):
        """Write reaction rates, coverages and pressures to file.
        Reaction rates computed at temperature T and pressure p.
        File written to current directory unless otherwise specified

        """

        if path != '':
            if not os.path.isdir(path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(path)

        T = self.params['temperature']
        p = self.params['pressure']

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
                self.reaction_terms(y=self.solution[t, :])
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

    def plot_transient(self, path=None):
        """Plot transient rates, coverages and pressures.
        If path is specified, figures are saved to path

        """
        font = {'family': 'sans-serif', 'weight': 'normal', 'size': 8}
        plt.rc('font', **font)
        mpl.rcParams['lines.markersize'] = 6
        mpl.rcParams['lines.linewidth'] = 1.5

        T = self.params['temperature']
        p = self.params['pressure']

        if path:
            if not os.path.isdir(path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(path)

        rates = np.zeros((len(self.times), len(self.reactions) * 2))
        for t in range(len(self.times)):
            self.reaction_terms(y=self.solution[t, :], T=T, p=p)
            for i in range(len(self.reactions)):
                rates[t, 2 * i] = self.rates[i, 0]
                rates[t, 2 * i + 1] = self.rates[i, 1]

        cmap = plt.get_cmap("Spectral", len(self.adsorbate_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(self.snames):
            if i in self.adsorbate_indices:
                ax.plot(self.times / 3600, self.solution[:, i], label=sname,
                        color=cmap(self.adsorbate_indices.index(i)))
        ax.legend(loc='center right', frameon=False, ncol=1)
        ax.set(xlabel='Time (hr)', ylabel='Coverage', title=(r'$T=%1.1f$ K' % T), xscale='log', ylim=(-0.1, 1.1))
        fig.tight_layout()
        if path:
            figname = path + 'coverages_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=300)

        cmap = plt.get_cmap("Accent", len(self.gas_indices))  # Spectral
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(self.snames):
            if i in self.gas_indices:
                ax.plot(self.times / 3600, self.solution[:, i], label=sname, color=cmap(self.gas_indices.index(i)))
        ax.legend(loc='center right', frameon=False, ncol=1)
        ax.set(xlabel='Time (hr)', ylabel='Pressure (bar)', title=('%1.1f K' % T), xscale='log')
        fig.tight_layout()
        if path:
            figname = path + 'pressures_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=300)

        cmap = plt.get_cmap("Accent", len(self.reactions) * 2)
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        for i, rname in enumerate([r for rname in self.reactions.keys() for r in [rname + '_fwd', rname + '_rev']]):
            ax.plot(self.times / 3600, rates[:, i], label=rname, color=cmap(i))
        ax.legend(loc='lower center', frameon=False, ncol=4)
        ax.set(xlabel='Time (hr)', ylabel='Rate (1/s)', title=('%1.1f K' % T), yscale='log')
        fig.tight_layout()
        if path:
            figname = path + 'surfrates_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=300)

