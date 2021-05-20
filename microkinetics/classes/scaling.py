from microkinetics.constants.physical_constants import *


class Scaling:

    def __init__(self, state_type, name=None, sigma=None, add_to_energy=None,
                 scaling_coeffs=None, scaling_reactions=None, dereference=False):
        """Initialises scaling relation class.
        This should be enough like a state that it can supply the reaction class with its free energy
        as a function of conditions.

        """
        self.state_type = state_type
        self.name = name
        self.scaling_coeffs = scaling_coeffs
        self.scaling_reactions = scaling_reactions
        self.dereference = dereference
        self.add_to_energy = add_to_energy
        self.path = None
        self.atoms = None
        self.freq = None
        self.i_freq = None
        self.Gelec = None
        self.Gtran = None
        self.Gvibr = None
        self.Grota = None
        self.Gfree = None
        self.Gzpe = None
        self.mass = None
        self.inertia = None
        self.shape = None
        self.sigma = sigma
        if self.state_type == 'gas':
            assert(self.sigma is not None)

    def calc_free_energy(self, T, p, verbose=False):
        """Calculates free energy.

        Saves value in eV."""
        assert(self.scaling_reactions is not None)
        assert(self.scaling_coeffs is not None)

        self.Gelec = self.scaling_coeffs['intercept']
        self.Gfree = 0.0
        for r in self.scaling_reactions.values():
            dEIS = r['reaction'].get_reaction_energy(T=T, p=p, verbose=verbose, etype='electronic') / (eVtokJ * 1.0e3)
            dGIS = r['reaction'].get_reaction_energy(T=T, p=p, verbose=verbose, etype='free') / (eVtokJ * 1.0e3)
            if self.dereference:
                ref_EIS = sum([reac.Gelec for reac in r['reaction'].reactants])
                ref_GIS = sum([reac.get_free_energy(T=T, p=p, verbose=verbose) for reac in r['reaction'].reactants])
            else:
                ref_EIS = 0.0
                ref_GIS = 0.0
            self.Gelec += r['multiplicity'] * (self.scaling_coeffs['gradient'] * dEIS + ref_EIS)
            self.Gfree += r['multiplicity'] * (-ref_EIS - dEIS + dGIS + ref_GIS)
        if self.add_to_energy:
            self.Gelec += self.add_to_energy
        self.Gfree += self.Gelec

        if verbose:
            print((self.name + ' elec: %1.2f eV') % self.Gelec)
            print((self.name + ' free: %1.2f eV') % self.Gfree)

    def get_free_energy(self, T, p, verbose=False):
        """Returns the free energy in eV.

        """
        self.calc_free_energy(T=T, p=p, verbose=verbose)
        return self.Gfree

    def set_energy_modifier(self, modifier, verbose=False):
        """Sets modifier to the electronic energy.

        Updates stored value in eV."""
        self.add_to_energy = modifier
        if self.Gelec:
            self.Gelec += self.add_to_energy
