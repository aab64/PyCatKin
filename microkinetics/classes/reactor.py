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
        """Calculates the scaling from reaction per site per second to per Pa per second.
        Reaction specific site_density will be added when reaction rates are computed.
        P = nRT/V, n = (site_density * catalyst_area) / NA, kB = R / NA

        """
        self.scaling = kB * T * self.catalyst_area / self.volume

    def rhs(self, adsorbate_kinetics):
        return lambda y: np.dot(adsorbate_kinetics(y), self.is_adsorbate)

    def save_pickle(self, path=None):
        """Save the reactor as a pickle object.

        """
        path = '' if path is None else path
        pickle.dump(self, open(path + self.name + '.pckl', 'wb'))


class InfiniteDilutionReactor(Reactor):

    def rhs(self, adsorbate_kinetics):
        def combined(y, t, T, p, inflow_state, verbose=False):
            return np.multiply(adsorbate_kinetics(y=y, T=T, p=p, verbose=verbose), self.is_adsorbate)
        return combined


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
        def combined(y, t, T, p, inflow_state, verbose=False):
            self.set_scaling(T=T)
            scaling = [1 if i else (self.scaling / bartoPa) for i in self.is_adsorbate]
            flow = np.array([0 if not self.is_gas[i] else (inflow_state[i] - y[i]) / self.residence_time
                             for i in range(len(self.is_gas))])
            return np.multiply(adsorbate_kinetics(y=y, T=T, p=p, verbose=verbose), np.array(scaling)) + flow
        return combined
