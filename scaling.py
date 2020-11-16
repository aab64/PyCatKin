from physical_constants import eVtokJ


class Scaling:

    def __init__(self, state_type, name=None, sigma=None,
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

        self.Gfree = self.scaling_coeffs['intercept']
        for r in self.scaling_reactions.values():
            energy = r['reaction'].get_reaction_energy(T=T, p=p, verbose=verbose) / (eVtokJ * 1.0e3)
            if self.dereference:
                ref_energy = sum([reac.get_free_energy(T=T, p=p, verbose=verbose) for reac in r['reaction'].reactants])
            else:
                ref_energy = 0.0
            self.Gfree += r['multiplicity'] * (self.scaling_coeffs['gradient'] * energy + ref_energy)

        if verbose:
            print((self.name + ': %1.2f eV') % self.Gfree)

    def get_free_energy(self, T, p, verbose=False):
        """Returns the free energy in eV.

        """
        self.calc_free_energy(T=T, p=p, verbose=verbose)
        return self.Gfree
