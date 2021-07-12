from pycatkin.classes.reaction import *
import copy
import numpy as np


class Uncertainty:

    def __init__(self, sys=None, mu=0.0, sigma=0.01, nruns=1):
        """Initialises Uncertainty class.
        Uncertainty class stores the base case system object, and
        performs experiments be perturbing its parameters to quantify
        the impact of uncertainty in the parameters on the predictions.

        """

        self.sys = copy.deepcopy(sys)
        self.mu = mu
        self.sigma = sigma
        self.nruns = nruns
        self.noisy_sys = None
        self.state_noises = None

    def get_noise(self, noise_type='white'):
        """Samples from either a Gaussian or uniform distribution.
        Default is white noise with mean zero and variance sigma^2.

        Returns a noise value.
        """
        if noise_type == 'white':
            return np.random.normal(loc=self.mu, scale=self.sigma, size=None)
        elif noise_type == 'uniform':
            return np.random.uniform()
        return 0.0

    def get_correlated_state_noises(self):
        """Samples a white noise value to add to
        adsorbate energies. Multiplies this by a uniform
        variate to perturb each transition state energy.

        Returns a dictionary of state names and noises.
        """
        noise = self.get_noise(noise_type='white')
        state_noises = dict()
        for rname in self.sys.reactions.keys():
            intermediates = []
            transition_states = []
            if isinstance(self.sys.reactions[rname], ReactionDerivedReaction):
                intermediates = self.sys.reactions[rname].base_reaction.reactants + \
                                self.sys.reactions[rname].base_reaction.products
                if self.sys.reactions[rname].base_reaction.TS:
                    transition_states = self.sys.reactions[rname].base_reaction.TS
            else:
                intermediates = self.sys.reactions[rname].reactants + \
                                self.sys.reactions[rname].products
                if self.sys.reactions[rname].TS:
                    transition_states = self.sys.reactions[rname].TS
            for reac in intermediates:
                if reac.state_type == 'adsorbate':
                    if reac.name not in state_noises.keys():
                        state_noises[reac.name] = noise
            for reac in transition_states:
                if reac.name not in state_noises.keys():
                    frac = self.get_noise(noise_type='uniform')
                    state_noises[reac.name] = (noise * frac)
        return state_noises

    def set_correlated_state_noises(self, state_noises):
        """Obtains a dictionary of state names and noises.
        Copies the system object and sets the perturbations
        for state energies in the new noisy system.

        Returns the new noisy system.
        """
        noisy_sys = copy.deepcopy(self.sys)
        for rname in noisy_sys.reactions.keys():
            if isinstance(noisy_sys.reactions[rname], ReactionDerivedReaction):
                for reac in noisy_sys.reactions[rname].base_reaction.reactants + \
                            noisy_sys.reactions[rname].base_reaction.products:
                    if reac.state_type == 'adsorbate':
                        noise = state_noises[reac.name]
                        reac.set_energy_modifier(noise)
                if noisy_sys.reactions[rname].base_reaction.TS:
                    for reac in noisy_sys.reactions[rname].base_reaction.TS:
                        noise = state_noises[reac.name]
                        reac.set_energy_modifier(noise)
            else:
                for reac in noisy_sys.reactions[rname].reactants + \
                            noisy_sys.reactions[rname].products:
                    if reac.state_type == 'adsorbate':
                        noise = state_noises[reac.name]
                        reac.set_energy_modifier(noise)
                if noisy_sys.reactions[rname].TS:
                    for reac in noisy_sys.reactions[rname].TS:
                        noise = state_noises[reac.name]
                        reac.set_energy_modifier(noise)
        return noisy_sys

    def get_noisy_sys_samples(self):
        """Runs the ODE solver for the original system and
         generates nruns copies, with noisy state energies.
         Runs the ODE solver for the noisy systems. Saves the
         perturbations used in each run.

        """
        self.sys.solve_odes()
        self.noisy_sys = dict()
        self.state_noises = dict()
        self.noisy_sys[0] = copy.deepcopy(self.sys)
        for run in range(1, self.nruns + 1):
            self.state_noises[run] = copy.deepcopy(self.get_correlated_state_noises())
            self.noisy_sys[run] = copy.deepcopy(self.set_correlated_state_noises(state_noises=self.state_noises[run]))
            self.noisy_sys[run].solve_odes()
        self.state_noises[0] = dict(zip(self.state_noises[1].keys(), np.zeros(len(self.state_noises[1].keys()))))

    def get_mean_property_value(self, property_handle):
        """Uses the function property_handle to compute the
        property values for all noisy systems.

        Returns the values, their mean and standard deviation.
        """
        self.get_noisy_sys_samples()
        property_values = np.array([property_handle(self.noisy_sys[i]) for i in self.noisy_sys.keys()])
        mean_property_value = np.mean(property_values[1::])
        std_property_value = np.std(property_values[1::])
        return property_values, mean_property_value, std_property_value
