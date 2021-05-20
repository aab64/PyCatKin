import copy
import numpy as np
import pickle
from microkinetics.constants.physical_constants import *


class Reactor:

    def __init__(self, residence_time=None, flow_rate=None, volume=None, catalyst_area=None, name=None):
        """Initialises generic Reactor class.

        """
        if name is None:
            name = 'reactor'
        self.residence_time = residence_time
        self.flow_rate = flow_rate
        self.volume = volume
        self.catalyst_area = catalyst_area
        self.name = name
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
        return lambda y: np.multiply(adsorbate_kinetics(y), self.is_adsorbate)

    def jacobian(self, adsorbate_jacobian):
        return lambda y: np.multiply(adsorbate_jacobian(y),
                                     np.transpose(np.tile(self.is_adsorbate, (len(self.is_adsorbate), 1))))

    def get_dynamic_indices(self, adsorbate_indices, gas_indices):
        """Returns which indicies in the solution vector have transient change.

        """
        return copy.copy(adsorbate_indices)

    def save_pickle(self, path=None):
        """Save the reactor as a pickle object.

        """
        path = '' if path is None else path
        pickle.dump(self, open(path + self.name + '.pckl', 'wb'))


class InfiniteDilutionReactor(Reactor):

    def rhs(self, adsorbate_kinetics):
        def combined(t, y, T, inflow_state):
            return np.multiply(adsorbate_kinetics(y=y), self.is_adsorbate)
        return combined

    def jacobian(self, adsorbate_jacobian):
        def combined(t, y, T):
            return np.multiply(adsorbate_jacobian(y=y),
                               np.transpose(np.tile(self.is_adsorbate, (len(self.is_adsorbate), 1))))
        return combined

    def get_dynamic_indices(self, adsorbate_indices, gas_indices):
        """Returns which indicies in the solution vector have transient change.

        """
        return copy.copy(adsorbate_indices)


class CSTReactor(Reactor):

    def __init__(self, residence_time=None, flow_rate=None, volume=None, catalyst_area=None, name=None):
        """Initialises Continuously Stirred Tank (CTR) Reactor class.

        """
        super(CSTReactor, self).__init__(residence_time=residence_time, flow_rate=flow_rate, volume=volume,
                                         catalyst_area=catalyst_area, name=name)
        if self.residence_time is None:
            assert(self.flow_rate is not None and self.volume is not None)
            print('Computing residence time from flow rate and volume, assuming SI units...')
            self.residence_time = self.volume / self.flow_rate

    def rhs(self, adsorbate_kinetics):
        def combined(t, y, T, inflow_state):
            ny = max(y.shape)
            y = y.reshape((ny, 1))
            self.set_scaling(T=T)
            scaling = [1 if i else (self.scaling / bartoPa) for i in self.is_adsorbate]
            flow = np.array([0 if not self.is_gas[i] else
                             (inflow_state[i] - y[i, 0]) / self.residence_time
                             for i in range(len(self.is_gas))])
            return np.multiply(adsorbate_kinetics(y=y), np.array(scaling)) + flow
        return combined

    def jacobian(self, adsorbate_jacobian):
        def combined(t, y, T):
            ny = max(y.shape)
            y = y.reshape((ny, 1))
            self.set_scaling(T=T)
            scaling = [1 if i else (self.scaling / bartoPa) for i in self.is_adsorbate]
            flow = np.array([0 if not self.is_gas[i] else
                             -1.0 / self.residence_time
                             for i in range(len(self.is_gas))])
            return np.multiply(adsorbate_jacobian(y=y),
                               np.transpose(np.tile(scaling, (len(scaling), 1)))) + np.diag(flow)
        return combined

    def get_dynamic_indices(self, adsorbate_indices, gas_indices):
        """Returns which indicies in the solution vector have transient change.

        """
        return copy.copy(adsorbate_indices) + copy.copy(gas_indices)
