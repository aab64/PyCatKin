from pycatkin.classes.reactor import *
import os
import copy
import pickle
from scipy.integrate import solve_ivp, ode
from scipy.optimize import least_squares
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd


class System:

    def __init__(self, path_to_pickle=None):
        """Initialises System class.
        System class stores the states, reactions and reactor.
        It is used to construct the reaction rates and species ODEs
        and solve the transient or steady-state dynamics.
        If path_to_pickle is defined, the pickled object is loaded.

        """

        if path_to_pickle:
            assert (os.path.isfile(path_to_pickle))
            newself = pickle.load(open(path_to_pickle, 'rb'))
            assert (isinstance(newself, System))
            for att in newself.__dict__.keys():
                setattr(self, att, getattr(newself, att))
        else:
            self.states = None
            self.reactions = None
            self.reactor = None
            self.energy_landscapes = None
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
            self.set_parameters()

    def add_state(self, state):
        """Adds a state to the dictionary of states
        and adds its name to the list of names.
        Checks state names are unique.

        """

        if self.states is None:
            self.states = dict()
        if self.snames is None:
            self.snames = []
        if self.params['verbose']:
            print('Adding state %s.' % state.name)
        self.states[state.name] = state
        if state.name in self.snames:
            raise ValueError('Found two copies of state %s. State names must be unique!' % state.name)
        else:
            self.snames = sorted(self.snames + [state.name])

    def add_reaction(self, reaction):
        """Adds a reaction to the dictionary of reactions.

        """

        if self.reactions is None:
            self.reactions = dict()
        if self.params['verbose']:
            print('Adding reaction %s.' % reaction.name)
        self.reactions[reaction.name] = reaction

    def add_reactor(self, reactor):
        """Adds a reactor.

        """

        if self.params['verbose']:
            print('Adding the reactor.')
        self.reactor = reactor

    def add_energy_landscape(self, energy_landscape):
        """Adds an energy landscape to the dictionary of reactions.

        """

        if self.energy_landscapes is None:
            self.energy_landscapes = dict()
        if self.params['verbose']:
            print('Adding energy landscape %s.' % energy_landscape.name)
        self.energy_landscapes[energy_landscape.name] = energy_landscape

    def names_to_indices(self):
        """Assigns indicies corresponding to the species for
        easier access to elements of the solution vector.

        """

        self.species_map = dict()
        for r in self.reactions.keys():
            yreac = [self.snames.index(i.name) for i in self.reactions[r].reactants
                     if i.state_type == 'adsorbate' or i.state_type == 'surface']
            preac = [self.snames.index(i.name) for i in self.reactions[r].reactants
                     if i.state_type == 'gas']
            yprod = [self.snames.index(i.name) for i in self.reactions[r].products
                     if i.state_type == 'adsorbate' or i.state_type == 'surface']
            pprod = [self.snames.index(i.name) for i in self.reactions[r].products
                     if i.state_type == 'gas']
            self.species_map[r] = {'yreac': yreac,
                                   'yprod': yprod,
                                   'preac': preac,
                                   'pprod': pprod,
                                   'site_density': 1.0 / self.reactions[r].area if self.reactions[r].area else 0.0,
                                   'scaling': self.reactions[r].scaling,
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
            is_adsorbate = [1 if i in self.adsorbate_indices else
                            0 for i in range(len(self.snames))]
        else:
            is_adsorbate = np.zeros(len(self.snames))
        if self.gas_indices is not None:
            self.gas_indices = list(set(self.gas_indices))
            is_gas = [1 if i in self.gas_indices else
                      0 for i in range(len(self.snames))]
        else:
            is_gas = np.zeros(len(self.snames))
        self.reactor.set_indices(is_adsorbate=is_adsorbate, is_gas=is_gas)
        self.dynamic_indices = self.reactor.get_dynamic_indices(self.adsorbate_indices, self.gas_indices)

    def set_parameters(self, times=None, start_state=None, inflow_state=None, T=293.15, p=101325.0,
                       use_jacobian=True, ode_solver='solve_ivp', nsteps=1e4, rtol=1e-8, atol=1e-10,
                       xtol=1e-8, ftol=1e-8, verbose=False):
        """Store simulation conditions, solver tolerances and verbosity.

        """

        self.params = dict()
        self.params['times'] = copy.deepcopy(times)
        self.params['start_state'] = copy.deepcopy(start_state)
        self.params['inflow_state'] = copy.deepcopy(inflow_state)
        self.params['temperature'] = T
        self.params['pressure'] = p
        self.params['rtol'] = rtol
        self.params['atol'] = atol
        self.params['xtol'] = xtol
        self.params['ftol'] = ftol
        self.params['jacobian'] = use_jacobian
        self.params['nsteps'] = int(nsteps)
        self.params['ode_solver'] = ode_solver
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
                self.reactions[r].calc_rate_constants(T=self.params['temperature'],
                                                      p=self.params['pressure'],
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

        for rind, r in enumerate(self.species_map.keys()):
            self.rates[rind, 0] = self.rate_constants[r]['kfwd'] + self.species_map[r]['perturbation']
            self.rates[rind, 1] = self.rate_constants[r]['krev'] * (1.0 + self.species_map[r]['perturbation'] /
                                                                    self.rate_constants[r]['kfwd'])
            for i in self.species_map[r]['yreac']:  # Forward rate species coverages
                self.rates[rind, 0] *= y[i, 0]
            for i in self.species_map[r]['yprod']:  # Reverse rate species coverages
                self.rates[rind, 1] *= y[i, 0]
            for i in self.species_map[r]['preac']:  # Forward rate species pressures
                self.rates[rind, 0] *= (y[i, 0] * bartoPa)
            for i in self.species_map[r]['pprod']:  # Reverse rate species pressures
                self.rates[rind, 1] *= (y[i, 0] * bartoPa)

    def species_odes(self, y):
        """Constructs ODEs for adsorbate coverages and pressures
        from the reaction rates.

        Returns array of species ODEs."""

        self.reaction_terms(y=y)
        net_rates = self.rates[:, 0] - self.rates[:, 1]

        ny = max(y.shape)
        dy_dt = np.zeros(ny)

        for rind, rinfo in enumerate(self.species_map.values()):
            for sp in rinfo['yreac']:  # Species consumed
                dy_dt[sp] -= net_rates[rind] * rinfo['scaling']
            for sp in rinfo['yprod']:  # Species formed
                dy_dt[sp] += net_rates[rind] * rinfo['scaling']
            for sp in rinfo['preac']:
                dy_dt[sp] -= net_rates[rind] * rinfo['scaling'] * rinfo['site_density']
            for sp in rinfo['pprod']:
                dy_dt[sp] += net_rates[rind] * rinfo['scaling'] * rinfo['site_density']
        return dy_dt

    def reaction_derivatives(self, y):
        """Constructs derivative of reactions wrt each species
        by multiplying rate constants by reactant coverages/pressures.

        Returns an (Nr x Ns) array of derivatives."""

        self.check_rate_constants()

        ny = max(y.shape)
        y = y.reshape((ny, 1))
        dr_dtheta = np.zeros((len(self.species_map), ny))

        def prodfun(reac, vartype, species):
            val = 1.0
            scaling = 1.0
            nsp = len(self.species_map[reac][vartype])
            for j in range(nsp):
                if j != species:
                    val *= y[self.species_map[reac][vartype][j], 0]
                    if vartype in ['preac', 'pprod']:
                        scaling = bartoPa
            return val * scaling

        for rind, r in enumerate(self.species_map.keys()):
            kfwd = self.rate_constants[r]['kfwd'] + self.species_map[r]['perturbation']
            krev = self.rate_constants[r]['krev'] * (1.0 + self.species_map[r]['perturbation'] /
                                                     self.rate_constants[r]['kfwd'])

            yfwd = prodfun(reac=r, vartype='yreac', species=None)
            yrev = prodfun(reac=r, vartype='yprod', species=None)
            pfwd = prodfun(reac=r, vartype='preac', species=None)
            prev = prodfun(reac=r, vartype='pprod', species=None)

            for ind, i in enumerate(self.species_map[r]['yreac']):
                dr_dtheta[rind, i] += kfwd * pfwd * prodfun(reac=r, vartype='yreac', species=ind)
            for ind, i in enumerate(self.species_map[r]['yprod']):
                dr_dtheta[rind, i] -= krev * prev * prodfun(reac=r, vartype='yprod', species=ind)
            for ind, i in enumerate(self.species_map[r]['preac']):
                dr_dtheta[rind, i] += kfwd * yfwd * prodfun(reac=r, vartype='preac', species=ind)
            for ind, i in enumerate(self.species_map[r]['pprod']):
                dr_dtheta[rind, i] -= krev * yrev * prodfun(reac=r, vartype='pprod', species=ind)
        return dr_dtheta

    def species_jacobian(self, y):
        """Constructs derivatives of species ODEs
        for adsorbate coverages and pressures.

        Returns Jacobian with shape (Ns x Ns)."""

        dr_dtheta = self.reaction_derivatives(y=y)

        ny = max(y.shape)
        jac = np.zeros((ny, ny))
        for rind, rinfo in enumerate(self.species_map.values()):
            for sp1 in range(ny):
                for sp2 in rinfo['yreac']:  # Species consumed
                    jac[sp2, sp1] -= dr_dtheta[rind, sp1] * rinfo['scaling']
                for sp2 in rinfo['yprod']:  # Species formed
                    jac[sp2, sp1] += dr_dtheta[rind, sp1] * rinfo['scaling']
                for sp2 in rinfo['preac']:
                    jac[sp2, sp1] -= dr_dtheta[rind, sp1] * rinfo['scaling'] * rinfo['site_density']
                for sp2 in rinfo['pprod']:
                    jac[sp2, sp1] += dr_dtheta[rind, sp1] * rinfo['scaling'] * rinfo['site_density']
        return jac

    def solve_odes(self):
        """Wrapper for ODE integrator.

        """

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

        solfun = lambda tval, yval: self.reactor.rhs(self.species_odes)(t=tval, y=yval, T=self.params['temperature'],
                                                                        inflow_state=yinflow)
        jacfun = lambda tval, yval: self.reactor.jacobian(self.species_jacobian)(t=tval, y=yval,
                                                                                 T=self.params['temperature'])

        # Create ODE solver
        if self.params['ode_solver'] == 'solve_ivp':
            sol = solve_ivp(fun=solfun, jac=jacfun if self.params['jacobian'] else None,
                            t_span=(self.params['times'][0], self.params['times'][-1]),
                            y0=yinit, method='BDF',
                            rtol=self.params['rtol'], atol=self.params['atol'])
            if self.params['verbose']:
                print(sol.message)
            self.times = sol.t
            self.solution = np.transpose(sol.y)
        elif self.params['ode_solver'] == 'ode':
            sol = ode(f=solfun, jac=jacfun if self.params['jacobian'] else None)
            sol.set_integrator('lsoda', method='bdf', rtol=self.params['rtol'], atol=self.params['atol'])
            sol.set_initial_value(yinit, self.params['times'][0])
            self.times = np.concatenate((np.zeros(1),
                                         np.logspace(start=np.log10(self.params['times'][0]
                                                                    if self.params['times'][0] else 1.0e-8),
                                                     stop=np.log10(self.params['times'][-1]),
                                                     num=self.params['nsteps'],
                                                     endpoint=True)))
            self.solution = np.zeros((self.params['nsteps'] + 1,
                                      len(self.snames)))
            self.solution[0, :] = yinit
            i = 1
            while sol.successful() and i <= self.params['nsteps']:
                sol.integrate(self.times[i])
                self.solution[i, :] = sol.y
                i += 1
        else:
            raise RuntimeError('Unknown ODE solver specified. Please use solve_ivp or ode, or add a new option here.')

        if self.params['verbose']:
            print('=========\nFinal conditions:\n')
            for s, sname in enumerate(self.snames):
                print('%15s : %9.2e' % (sname, self.solution[-1][s]))

    def find_steady(self, store_steady=False, plot_comparison=False, path=None):
        """Solve for the steady state solution

        Returns the steady state solution."""

        self.conditions = None  # Force rate constants to be recalculated

        # Establish an initial guess
        if self.solution is not None:
            y_guess = copy.deepcopy(self.solution[-1, self.dynamic_indices])
            full_steady = copy.deepcopy(self.solution[-1, :])
        else:
            y_guess = np.zeros(len(self.dynamic_indices))
            full_steady = np.zeros(len(self.adsorbate_indices) + len(self.gas_indices))

        # Set inflow mole fractions to zero if not specified
        yinflow = np.zeros(len(self.snames))
        if self.params['inflow_state']:
            yinflow = np.zeros(len(self.snames))
            for s in self.params['inflow_state'].keys():
                yinflow[self.snames.index(s)] = self.params['inflow_state'][s]

        # Function to solve
        # Adds variable y to full solution vector and passes to function
        # Returns only variable parts of function
        def func(y):
            full_steady[self.dynamic_indices] = y
            return self.reactor.rhs(self.species_odes)(t=0, y=full_steady, T=self.params['temperature'],
                                                       inflow_state=yinflow)[self.dynamic_indices]

        # Jacobian to use
        if self.params['jacobian']:
            def jacfun(y):
                full_steady[self.dynamic_indices] = y
                full_jacobian = self.reactor.jacobian(self.species_jacobian)(t=0, y=full_steady,
                                                                             T=self.params['temperature'])
                return np.array([[full_jacobian[i1, i2] for i1 in self.dynamic_indices]
                                 for i2 in self.dynamic_indices])
        else:
            jacfun = '3-point'

        sol = least_squares(fun=func, x0=y_guess, jac=jacfun, method='trf',
                            xtol=self.params['xtol'], ftol=self.params['ftol'],
                            max_nfev=np.max((int(1e4), 100 * len(y_guess))))
        y_steady = sol.x
        mesg = sol.message
        ier = sol.nfev

        full_steady[self.dynamic_indices] = y_steady

        if store_steady:
            self.full_steady = full_steady

        if self.params['verbose']:
            print('Results of steady state search...')
            print('- At %1.0f K: %s, %1i' % (self.params['temperature'], mesg, ier))
            print('- Cost function value at steady state: %.3g' % sol.cost)
            print(sol.fun)
            print('- Norm of function value at steady state: %.3g' % np.linalg.norm(func(y_steady)))
            print('- Norm of guess minus steady state: %.3g' % np.linalg.norm(y_guess - y_steady))

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
                    ax.plot(self.times, [full_steady[i] for t in self.times], label='',
                            color=cmap(self.dynamic_indices.index(i)), linestyle=':')
            ax.legend(frameon=False, loc='center right')
            ax.set(xlabel='Time (s)', xscale='log',
                   ylabel='Coverage', yscale='log', ylim=(1e-6, 1e1),
                   title=(r'$T=%1.0f$ K' % self.params['temperature']))
            fig.tight_layout()
            if path:
                fig.savefig((path + 'SS_vs_transience_%1.1fK.png') % self.params['temperature'],
                            format='png', dpi=300)

        return full_steady

    def run_and_return_tof(self, tof_terms, ss_solve=False):
        """Integrate or solve for the steady state and
        compute the TOF by summing steps in tof_terms

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
                tof += self.rates[rind, 0] - self.rates[rind, 1]
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

    def activity(self, tof_terms, ss_solve=False):
        """Calculate the activity from the TOF

        Returns the activity."""

        self.conditions = None  # Force rate constants to be recalculated

        tof = self.run_and_return_tof(tof_terms=tof_terms, ss_solve=ss_solve)

        activity = (np.log((h * tof) / (kB * self.params['temperature'])) *
                    (R * self.params['temperature'])) * 1.0e-3 / eVtokJ

        return activity

    def write_results(self, path=''):
        """Write reaction rates, coverages and pressures to file.
        File written to current directory unless otherwise specified.

        """

        if path != '':
            if not os.path.isdir(path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(path)

        T = self.params['temperature']
        p = self.params['pressure']

        rfile = path + 'rates_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.csv'
        cfile = path + 'coverages_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.csv'
        pfile = path + 'pressures_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.csv'

        rheader = ['Time (s)'] + [j for k in [i.split(',')
                                              for i in [(r.name + '_fwd,' + r.name + '_rev')
                                                        for r in self.reactions.values()]]
                                  for j in k]
        cheader = ['Time (s)'] + [s for i, s in enumerate(self.snames) if i in self.adsorbate_indices]
        pheader = ['Time (s)'] + [s for i, s in enumerate(self.snames) if i in self.gas_indices]

        rmat = np.zeros((len(self.times), 2 * self.rates.shape[0]))
        for t in range(len(self.times)):
            self.reaction_terms(y=self.solution[t, :])
            rmat[t, :] = self.rates.flatten()

        times = self.times.reshape(len(self.times), 1)

        df = pd.DataFrame(np.concatenate((times, rmat), axis=1), columns=rheader)
        df.to_csv(path_or_buf=rfile, sep=',', header=True, index=False)
        df = pd.DataFrame(np.concatenate((times, self.solution[:, self.adsorbate_indices]), axis=1), columns=cheader)
        df.to_csv(path_or_buf=cfile, sep=',', header=True, index=False)
        df = pd.DataFrame(np.concatenate((times, self.solution[:, self.gas_indices]), axis=1), columns=pheader)
        df.to_csv(path_or_buf=pfile, sep=',', header=True, index=False)

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

        if path is not None and path != '':
            if not os.path.isdir(path):
                print('Directory does not exist. Will try creating it...')
                os.mkdir(path)

        rates = np.zeros((len(self.times), len(self.reactions) * 2))
        for t in range(len(self.times)):
            self.reaction_terms(y=self.solution[t, :])
            for i in range(len(self.reactions)):
                rates[t, 2 * i] = self.rates[i, 0]
                rates[t, 2 * i + 1] = self.rates[i, 1]

        cmap = plt.get_cmap("tab20", len(self.adsorbate_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(self.snames):
            if i in self.adsorbate_indices and max(self.solution[:, i]) > 0.01:
                ax.plot(self.times / 3600, self.solution[:, i], label=sname,
                        color=cmap(self.adsorbate_indices.index(i)))
        ax.legend(loc='best', frameon=False, ncol=1)
        ax.set(xlabel='Time (hr)', xscale='log',
               ylabel='Coverage', ylim=(-0.1, 1.1),
               title=(r'$T=%1.1f$ K' % T))
        fig.tight_layout()
        if path is not None:
            figname = path + 'coverages_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=600)

        cmap = plt.get_cmap("tab20", len(self.gas_indices))
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        for i, sname in enumerate(self.snames):
            if i in self.gas_indices:
                ax.plot(self.times / 3600, self.solution[:, i], label=sname,
                        color=cmap(self.gas_indices.index(i)))
        ax.legend(loc='center right', frameon=False, ncol=1)
        ax.set(xlabel='Time (hr)', xscale='log',
               ylabel='Pressure (bar)',
               title=('T = %1.1f K' % T))
        fig.tight_layout()
        if path is not None:
            figname = path + 'pressures_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=600)

        cmap = plt.get_cmap("tab20", len(self.reactions) * 2)
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        for i, rname in enumerate([r for rname in self.reactions.keys()
                                   for r in [rname + '_fwd', rname + '_rev']]):
            ax.plot(self.times / 3600, rates[:, i], label=rname, color=cmap(i))
        ax.legend(loc='lower center', frameon=False, ncol=4)
        yvals = ax.get_ylim()
        ax.set(xlabel='Time (hr)', xscale='log',
               ylabel='Rate (1/s)', yscale='log', ylim=(max(1e-10, yvals[0]), yvals[1]),
               title=('T = %1.1f K' % T))
        fig.tight_layout()
        if path is not None:
            figname = path + 'surfrates_' + ('%1.1f' % T) + 'K_' + ('%1.1f' % (p / bartoPa)) + 'bar.png'
            plt.savefig(figname, format='png', dpi=600)

    def save_pickle(self, path=None):
        """Save the system as a pickle object.

        """

        path = path if path is not None else ''
        pickle.dump(self, open(path + 'system' + '.pckl', 'wb'))
