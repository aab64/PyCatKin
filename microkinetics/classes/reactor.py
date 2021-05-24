from microkinetics.constants.physical_constants import *
import copy
import pickle
import os
import numpy as np


class Reactor:

    def __init__(self, name='reactor', volume=None, catalyst_area=None,
                 residence_time=None, flow_rate=None, path_to_pickle=None):
        """Initialises generic Reactor class.
        If path_to_pickle is defined, the pickled object is loaded.

        """

        if path_to_pickle:
            assert (os.path.isfile(path_to_pickle))
            newself = pickle.load(open(path_to_pickle, 'rb'))
            assert (isinstance(newself, Reactor))
            for att in newself.__dict__.keys():
                setattr(self, att, getattr(newself, att))
        else:
            self.name = name
            self.volume = volume
            self.catalyst_area = catalyst_area
            self.residence_time = residence_time
            self.flow_rate = flow_rate
            self.scaling = None
            self.is_adsorbate = None
            self.is_gas = None

    def set_scaling(self, T):
        """Scaling from per site per second to per Pa per second.
        site_density will be added when reaction rates are computed.
        P = nRT/V, n = (site_density * catalyst_area) / NA, kB = R / NA

        """

        self.scaling = kB * T * self.catalyst_area / self.volume

    def rhs(self, adsorbate_kinetics):
        """Construct the ODE right hand side (rhs).
        Multiplies the species ODEs by 1 if they are required or 0 if
        they are not required (e.g., pressure boundary conditions).

        Returns a function handle."""

        return lambda y: np.multiply(adsorbate_kinetics(y), self.is_adsorbate)

    def jacobian(self, adsorbate_jacobian):
        """Construct the Jacobian function.
        Multiplies the rows of the Jacobian by 1 if they are required or
        0 if they are not required (e.g., pressure boundary conditions).

        Returns a function handle."""

        return lambda y: np.multiply(adsorbate_jacobian(y),
                                     np.transpose(np.tile(self.is_adsorbate,
                                                          (len(self.is_adsorbate), 1))))

    def get_dynamic_indices(self, adsorbate_indices, gas_indices):
        """Returns which indicies in the solution vector vary with time
        (e.g., not pressure boundary conditions for ID reactors).

        """

        return copy.copy(adsorbate_indices)

    def save_pickle(self, path=None):
        """Save the reactor as a pickle object.

        """

        path = path if path else ''
        pickle.dump(self, open(path + 'reactor_' + self.name + '.pckl', 'wb'))


class InfiniteDilutionReactor(Reactor):

    def rhs(self, adsorbate_kinetics):
        """Construct the ODE right hand side (rhs).
        Multiplies the species ODEs by 1 if they are required or 0 if
        they are not required (e.g., pressure boundary conditions).

        Returns a function handle."""

        def combined(t, y, T, inflow_state):
            return np.multiply(adsorbate_kinetics(y=y), self.is_adsorbate)

        return combined

    def jacobian(self, adsorbate_jacobian):
        """Construct the Jacobian function.
        Multiplies the rows of the Jacobian by 1 if they are required or
        0 if they are not required (e.g., pressure boundary conditions).

        Returns a function handle."""

        def combined(t, y, T):
            return np.multiply(adsorbate_jacobian(y=y),
                               np.transpose(np.tile(self.is_adsorbate, (len(self.is_adsorbate), 1))))

        return combined

    def get_dynamic_indices(self, adsorbate_indices, gas_indices):
        """Returns which indicies in the solution vector have transient change.

        """

        return copy.copy(adsorbate_indices)


class CSTReactor(Reactor):

    def __init__(self, name='reactor', volume=None, catalyst_area=None,
                 residence_time=None, flow_rate=None):
        """Initialises Continuously Stirred Tank (CST) Reactor class.

        """

        super(CSTReactor, self).__init__(residence_time=residence_time, flow_rate=flow_rate, volume=volume,
                                         catalyst_area=catalyst_area, name=name)
        if self.residence_time is None:
            assert(self.flow_rate is not None and
                   self.volume is not None)
            print('Computing residence time from flow rate and volume, assuming SI units...')
            self.residence_time = self.volume / self.flow_rate

    def rhs(self, adsorbate_kinetics):
        """Construct the ODE right hand side (rhs).
        Multiplies the species ODEs by 1 if they are surface equations
        or sigma if they are gas equations. For gas equations, adds flow.

        Returns a function handle."""

        def combined(t, y, T, inflow_state):
            ny = max(y.shape)
            y = y.reshape((ny, 1))
            self.set_scaling(T=T)
            scaling = [1 if i else (self.scaling / bartoPa)
                       for i in self.is_adsorbate]
            flow = np.array([0 if not self.is_gas[i] else
                             (inflow_state[i] - y[i, 0]) / self.residence_time
                             for i in range(len(self.is_gas))])
            return np.multiply(adsorbate_kinetics(y=y), np.array(scaling)) + flow

        return combined

    def jacobian(self, adsorbate_jacobian):
        """Construct the Jacobian function.
        Multiplies the rows of the Jacobian by 1 if they are surface equations
        or sigma if they are gas equations. For rows corresponding to gases,
        adds flow derivative (-1/tau).

        Returns a function handle."""

        def combined(t, y, T):
            ny = max(y.shape)
            y = y.reshape((ny, 1))
            self.set_scaling(T=T)
            scaling = [1 if i else (self.scaling / bartoPa)
                       for i in self.is_adsorbate]
            flow = np.array([0 if not self.is_gas[i] else
                             -1.0 / self.residence_time
                             for i in range(len(self.is_gas))])
            return np.multiply(adsorbate_jacobian(y=y),
                               np.transpose(np.tile(scaling, (len(scaling), 1)))) + np.diag(flow)

        return combined

    def get_dynamic_indices(self, adsorbate_indices, gas_indices):
        """Returns which indicies in the solution
        vector experience transient changes.

        """

        return copy.copy(adsorbate_indices) + copy.copy(gas_indices)
