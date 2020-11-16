import pickle
from microkinetics.functions.rate_constants import *


class Reaction:

    def __init__(self, reac_type, reversible=True, reactants=None, products=None, TS=None,
                 area=1.0e-19, name='reaction', scaling=1.0):
        """Initialises Reaction class.
        Reaction class stores the states involved in the reaction,
        the rate constants, reaction energy and barrier.

        """
        self.reac_type = reac_type
        self.reversible = reversible
        self.reactants = reactants
        self.products = products
        self.TS = TS
        self.area = area
        self.name = name
        self.scaling = scaling
        self.kfwd = None
        self.krev = None
        self.Keq = None
        self.dGrxn = None
        self.dGa_fwd = None
        self.dGa_rev = None

    def calc_reaction_energy(self, T, p, verbose=False):
        """Computes reaction energies and barriers in J/mol.

        """
        Greac = sum([i.get_free_energy(T=T, p=p, verbose=verbose) for i in self.reactants])
        if self.reversible:
            Gprod = sum([i.get_free_energy(T=T, p=p, verbose=verbose) for i in self.products])
            self.dGrxn = (Gprod - Greac) * eVtokJ * 1.0e3
        if self.TS is not None:
            GTS = sum([i.get_free_energy(T=T, p=p, verbose=verbose) for i in self.TS])
            self.dGa_fwd = (GTS - Greac) * eVtokJ * 1.0e3
            if self.reversible:
                self.dGa_rev = (GTS - Gprod) * eVtokJ * 1.0e3
        else:
            self.dGa_fwd = 0.0
            self.dGa_rev = 0.0
        if verbose:
            print('---------------------')
            print(self.name)
            print('reactants:')
            for i in self.reactants:
                print('* ' + i.name + ', ' + i.state_type)
            print('products:')
            for i in self.products:
                print('* ' + i.name + ', ' + i.state_type)
            if self.TS is not None:
                for i in self.TS:
                    print('* ' + i.name + ', ' + i.state_type)
            print('dGfwd: % 1.2f (kJ/mol)' % (self.dGa_fwd * 1.0e-3))
            if self.reversible:
                print('dGrev: % 1.2f (kJ/mol)' % (self.dGa_rev * 1.0e-3))
                print('dGrxn: % 1.2f (kJ/mol)' % (self.dGrxn * 1.0e-3))
            print('---------------------')

    def calc_rate_constants(self, T, p, verbose=False):
        """Computes reaction rate constants.

        """
        self.calc_reaction_energy(T=T, p=p, verbose=verbose)
        if self.reac_type == 'adsorption':
            gassp = [s for s in self.reactants if s.state_type == 'gas']
            assert(len(gassp) == 1)
            self.kfwd = kads(T=T, mass=gassp[0].mass, area=self.area)
            if self.reversible:
                self.Keq = keq_therm(T=T, rxn_en=self.dGrxn)
                self.krev = k_from_eq_rel(kknown=self.kfwd, Keq=self.Keq, direction='forward')
            else:
                self.krev = 0.0
        elif self.reac_type == 'desorption':
            gassp = [s for s in self.products if s.state_type == 'gas']
            assert(len(gassp) == 1)
            if self.reversible:
                self.krev = kads(T=T, mass=gassp[0].mass, area=self.area)
                self.Keq = keq_therm(T=T, rxn_en=self.dGrxn)
                self.kfwd = k_from_eq_rel(kknown=self.krev, Keq=self.Keq, direction='reverse')
            else:
                print('Inconsistent definition of desorption from equilibrium relation!')
                self.Keq = None
                self.krev = 0.0
        else:
            self.kfwd = karr(T=T, prefac=prefactor(T), barrier=self.dGa_fwd)
            if self.reversible:
                self.Keq = keq_therm(T=T, rxn_en=self.dGrxn)
                self.krev = k_from_eq_rel(kknown=self.kfwd, Keq=self.Keq, direction='forward')
            else:
                self.krev = 0.0

    def get_reaction_energy(self, T, p, verbose=False):
        """Returns the reaction energy in J/mol.

        """
        self.calc_reaction_energy(T=T, p=p, verbose=verbose)
        return self.dGrxn

    def get_reaction_barriers(self, T, p, verbose=False):
        """Returns the reaction barriers in J/mol.

        """
        self.calc_reaction_energy(T=T, p=p, verbose=verbose)
        return self.dGa_fwd, self.dGa_rev

    def save_pickle(self, path=None):
        """Save the reaction as a pickle object.

        """
        path = '' if path is None else path
        pickle.dump(self, open(path + self.name + '.pckl', 'wb'))


class UserDefinedReaction(Reaction):

    def __init__(self, reac_type, reversible=True, reactants=None, products=None, TS=None,
                 area=1.0e-19, name='reaction', scaling=1.0,
                 dGrxn_user=None, dGa_fwd_user=None, dGa_rev_user=None):
        """Initialises UserDefinedReaction class with energies specified by the user.
        Reaction class stores the states involved in the reaction,
        the rate constants, reaction energy and barrier.

        """
        super(UserDefinedReaction, self).__init__(reac_type=reac_type, reversible=reversible, reactants=reactants,
                                                  products=products, TS=TS, area=area, name=name, scaling=scaling)
        self.dGrxn_user = dGrxn_user
        self.dGa_fwd_user = dGa_fwd_user
        self.dGa_rev_user = dGa_rev_user

    def calc_reaction_energy(self, T, p, verbose=False):
        """Computes reaction energies and barriers in J/mol.

        """
        if self.reversible:
            if isinstance(self.dGrxn_user, dict):
                self.dGrxn = self.dGrxn_user[p][T] * eVtokJ * 1.0e3
            else:
                assert (self.dGrxn_user is not None)
                self.dGrxn = self.dGrxn_user * eVtokJ * 1.0e3

        if self.dGa_fwd_user is not None:
            if isinstance(self.dGa_fwd_user, dict):
                self.dGa_fwd = self.dGa_fwd_user[p][T] * eVtokJ * 1.0e3
            else:
                self.dGa_fwd = self.dGa_fwd_user * eVtokJ * 1.0e3
            if self.reversible:
                self.dGa_rev = (self.dGa_fwd - self.dGrxn) * eVtokJ * 1.0e3
        else:
            self.dGa_fwd = 0.0
            self.dGa_rev = 0.0

        if verbose:
            print('---------------------')
            print(self.name)
            print('reactants:')
            for i in self.reactants:
                print('* ' + i.name + ', ' + i.state_type)
            print('products:')
            for i in self.products:
                print('* ' + i.name + ', ' + i.state_type)
            if self.TS is not None:
                for i in self.TS:
                    print('* ' + i.name + ', ' + i.state_type)
            print('dGfwd: % 1.2f (kJ/mol)' % (self.dGa_fwd * 1.0e-3))
            if self.reversible:
                print('dGrev: % 1.2f (kJ/mol)' % (self.dGa_rev * 1.0e-3))
                print('dGrxn: % 1.2f (kJ/mol)' % (self.dGrxn * 1.0e-3))
            print('---------------------')
