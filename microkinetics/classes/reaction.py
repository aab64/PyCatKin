from microkinetics.functions.rate_constants import *
import pickle
import os


class Reaction:

    def __init__(self, name='reaction', reac_type=None, reversible=True,
                 reactants=None, products=None, TS=None,
                 area=1.0e-19, scaling=1.0, path_to_pickle=None):
        """Initialises Reaction class.
        Reaction class stores the states involved in the reaction,
        the rate constants, reaction energy and barrier.
        If path_to_pickle is defined, the pickled object is loaded.

        """

        if path_to_pickle:
            assert (os.path.isfile(path_to_pickle))
            newself = pickle.load(open(path_to_pickle, 'rb'))
            assert (isinstance(newself, Reaction))
            for att in newself.__dict__.keys():
                setattr(self, att, getattr(newself, att))
        else:
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
            self.dErxn = None
            self.dEa_fwd = None
            self.dEa_rev = None

    def calc_reaction_energy(self, T, p, verbose=False):
        """Computes reaction energies and barriers in J/mol.

        """

        Greac = sum([i.get_free_energy(T=T, p=p, verbose=verbose)
                     for i in self.reactants])
        Ereac = sum([i.Gelec for i in self.reactants])
        if self.reversible:
            Gprod = sum([i.get_free_energy(T=T, p=p, verbose=verbose)
                         for i in self.products])
            Eprod = sum([i.Gelec for i in self.products])
            self.dGrxn = (Gprod - Greac) * eVtokJ * 1.0e3
            self.dErxn = (Eprod - Ereac) * eVtokJ * 1.0e3
        if self.TS is not None:
            GTS = sum([i.get_free_energy(T=T, p=p, verbose=verbose)
                       for i in self.TS])
            ETS = sum([i.Gelec for i in self.TS])
            self.dGa_fwd = (GTS - Greac) * eVtokJ * 1.0e3
            self.dEa_fwd = (ETS - Ereac) * eVtokJ * 1.0e3
            if self.reversible:
                self.dGa_rev = (GTS - Gprod) * eVtokJ * 1.0e3
                self.dEa_rev = (ETS - Eprod) * eVtokJ * 1.0e3
        else:
            self.dGa_fwd = 0.0
            self.dGa_rev = 0.0
            self.dEa_fwd = 0.0
            self.dEa_rev = 0.0

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
            print('dEfwd: % 1.2f (kJ/mol)' % (self.dEa_fwd * 1.0e-3))
            if self.reversible:
                print('dGrev: % 1.2f (kJ/mol)' % (self.dGa_rev * 1.0e-3))
                print('dGrxn: % 1.2f (kJ/mol)' % (self.dGrxn * 1.0e-3))
                print('dErev: % 1.2f (kJ/mol)' % (self.dEa_rev * 1.0e-3))
                print('dErxn: % 1.2f (kJ/mol)' % (self.dErxn * 1.0e-3))
            print('---------------------')

    def calc_rate_constants(self, T, p, verbose=False):
        """Computes reaction rate constants.

        """

        self.calc_reaction_energy(T=T, p=p, verbose=verbose)

        if self.reac_type == 'adsorption':
            gassp = [s for s in self.reactants
                     if s.state_type == 'gas']
            assert(len(gassp) == 1)
            self.kfwd = kads(T=T, mass=gassp[0].mass, area=self.area)
            if self.dGa_fwd:
                if verbose:
                    print('Assuming activated adsorption has Arrhenius rate constant!')
                self.kfwd = karr(T=T, prefac=self.kfwd, barrier=np.max((self.dGa_fwd, 0.0)))
            if self.reversible:
                self.Keq = keq_therm(T=T, rxn_en=self.dGrxn)
                self.krev = k_from_eq_rel(kknown=self.kfwd, Keq=self.Keq, direction='forward')
            else:
                self.krev = 0.0
        elif self.reac_type == 'desorption':
            gassp = [s for s in self.products
                     if s.state_type == 'gas']
            assert(len(gassp) == 1)
            if self.reversible:
                self.krev = kads(T=T, mass=gassp[0].mass, area=self.area)
                if self.dGa_rev:
                    if verbose:
                        print('Assuming activated adsorption has Arrhenius rate constant!')
                    self.krev = karr(T=T, prefac=self.krev, barrier=np.max((self.dGa_rev, 0.0)))
                self.Keq = keq_therm(T=T, rxn_en=self.dGrxn)
                self.kfwd = k_from_eq_rel(kknown=self.krev, Keq=self.Keq, direction='reverse')
            else:
                print('Inconsistent definition of desorption from equilibrium relation!')
                self.Keq = None
                self.krev = 0.0
        else:
            self.kfwd = karr(T=T, prefac=prefactor(T), barrier=np.max((self.dGa_fwd, 0.0)))
            if self.reversible:
                self.Keq = keq_therm(T=T, rxn_en=self.dGrxn)
                self.krev = k_from_eq_rel(kknown=self.kfwd, Keq=self.Keq, direction='forward')
            else:
                self.krev = 0.0

    def get_reaction_energy(self, T, p, verbose=False, etype='free'):
        """Returns the reaction energy in J/mol.

        """

        self.calc_reaction_energy(T=T, p=p, verbose=verbose)

        if etype == 'electronic':
            return self.dErxn
        return self.dGrxn

    def get_reaction_barriers(self, T, p, verbose=False, etype='free'):
        """Returns the reaction barriers in J/mol.

        """

        self.calc_reaction_energy(T=T, p=p, verbose=verbose)

        if etype == 'electronic':
            return self.dEa_fwd, self.dEa_rev
        return self.dGa_fwd, self.dGa_rev

    def save_pickle(self, path=None):
        """Save the reaction as a pickle object.

        """

        path = path if path else ''
        pickle.dump(self, open(path + 'reaction_' + self.name + '.pckl', 'wb'))


class UserDefinedReaction(Reaction):

    def __init__(self, reac_type, reversible=True, reactants=None, products=None, TS=None,
                 area=1.0e-19, name='reaction', scaling=1.0,
                 dErxn_user=None, dEa_fwd_user=None, dEa_rev_user=None,
                 dGrxn_user=None, dGa_fwd_user=None, dGa_rev_user=None):
        """Initialises UserDefinedReaction class
        in which energies are specified by the user.

        """

        super(UserDefinedReaction, self).__init__(reac_type=reac_type, reversible=reversible, reactants=reactants,
                                                  products=products, TS=TS, area=area, name=name, scaling=scaling)
        self.dErxn_user = dErxn_user
        self.dEa_fwd_user = dEa_fwd_user
        self.dEa_rev_user = dEa_rev_user
        self.dGrxn_user = dGrxn_user
        self.dGa_fwd_user = dGa_fwd_user
        self.dGa_rev_user = dGa_rev_user

    def calc_reaction_energy(self, T, p, verbose=False):
        """Computes reaction energies and barriers in J/mol.

        """

        if self.reversible:
            if isinstance(self.dErxn_user, dict):
                self.dErxn = self.dErxn_user[T] * eVtokJ * 1.0e3
            else:
                if self.dErxn_user is not None:
                    self.dErxn = self.dErxn_user * eVtokJ * 1.0e3
            if isinstance(self.dGrxn_user, dict):
                self.dGrxn = self.dGrxn_user[T] * eVtokJ * 1.0e3
            else:
                if self.dGrxn_user is not None:
                    self.dGrxn = self.dGrxn_user * eVtokJ * 1.0e3
            if self.dErxn is None:
                assert(self.dGrxn is not None)
                self.dErxn = self.dGrxn
            if self.dGrxn is None:
                assert(self.dErxn is not None)
                self.dGrxn = self.dErxn

        if self.dEa_fwd_user is not None:
            if isinstance(self.dEa_fwd_user, dict):
                self.dEa_fwd = self.dEa_fwd_user[T] * eVtokJ * 1.0e3
            else:
                self.dEa_fwd = self.dEa_fwd_user * eVtokJ * 1.0e3
            if self.reversible:
                self.dEa_rev = (self.dEa_fwd - self.dErxn)
        else:
            self.dEa_fwd = 0.0
            self.dEa_rev = 0.0

        if self.dGa_fwd_user is not None:
            if isinstance(self.dGa_fwd_user, dict):
                self.dGa_fwd = self.dGa_fwd_user[T] * eVtokJ * 1.0e3
            else:
                self.dGa_fwd = self.dGa_fwd_user * eVtokJ * 1.0e3
            if self.reversible:
                self.dGa_rev = (self.dGa_fwd - self.dGrxn)
        else:
            self.dGa_fwd = 0.0
            self.dGa_rev = 0.0

        if self.dEa_fwd is None and self.dGa_fwd is not None:
            self.dEa_fwd = self.dGa_fwd
            self.dEa_rev = self.dGa_rev
        elif self.dEa_fwd is not None and self.dGa_fwd is None:
            self.dGa_fwd = self.dEa_fwd
            self.dGa_rev = self.dEa_rev

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
            print('dEfwd: % 1.2f (kJ/mol)' % (self.dEa_fwd * 1.0e-3))
            if self.reversible:
                print('dGrev: % 1.2f (kJ/mol)' % (self.dGa_rev * 1.0e-3))
                print('dGrxn: % 1.2f (kJ/mol)' % (self.dGrxn * 1.0e-3))
                print('dErev: % 1.2f (kJ/mol)' % (self.dEa_rev * 1.0e-3))
                print('dErxn: % 1.2f (kJ/mol)' % (self.dErxn * 1.0e-3))
            print('---------------------')


class ReactionDerivedReaction(Reaction):

    def __init__(self, reac_type, reversible=True, reactants=None, products=None, TS=None,
                 area=1.0e-19, name='reaction', scaling=1.0, base_reaction=None):
        """Initialises ReactionDerivedReaction class
        in which energies are specified by a different reaction.

        """

        super(ReactionDerivedReaction, self).__init__(reac_type=reac_type, reversible=reversible, reactants=reactants,
                                                      products=products, TS=TS, area=area, name=name, scaling=scaling)
        assert (base_reaction is not None)
        self.base_reaction = base_reaction
        if reversible != base_reaction.reversible:
            revcheck = lambda s: 'reversible' if s else 'irreversible'
            print('Warning! ReactionDerivedReaction is defined as %s but base Reaction is %s. Will use %s' %
                  (revcheck(reversible), revcheck(base_reaction.reversible), revcheck(base_reaction.reversible)))
            self.reversible = base_reaction.reversible

    def calc_reaction_energy(self, T, p, verbose=False):
        """Computes reaction energies and barriers in J/mol.

        """

        Greac = sum([i.get_free_energy(T=T, p=p, verbose=verbose)
                     for i in self.base_reaction.reactants])
        Ereac = sum([i.Gelec for i in self.base_reaction.reactants])
        if self.base_reaction.reversible:
            Gprod = sum([i.get_free_energy(T=T, p=p, verbose=verbose)
                         for i in self.base_reaction.products])
            Eprod = sum([i.Gelec for i in self.base_reaction.products])
            self.dGrxn = (Gprod - Greac) * eVtokJ * 1.0e3
            self.dErxn = (Eprod - Ereac) * eVtokJ * 1.0e3
        if self.base_reaction.TS is not None:
            GTS = sum([i.get_free_energy(T=T, p=p, verbose=verbose)
                       for i in self.base_reaction.TS])
            ETS = sum([i.Gelec for i in self.base_reaction.TS])
            self.dGa_fwd = (GTS - Greac) * eVtokJ * 1.0e3
            self.dEa_fwd = (ETS - Ereac) * eVtokJ * 1.0e3
            if self.base_reaction.reversible:
                self.dGa_rev = (GTS - Gprod) * eVtokJ * 1.0e3
                self.dEa_rev = (ETS - Eprod) * eVtokJ * 1.0e3
        else:
            self.dGa_fwd = 0.0
            self.dGa_rev = 0.0
            self.dEa_fwd = 0.0
            self.dEa_rev = 0.0

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
            print('dEfwd: % 1.2f (kJ/mol)' % (self.dEa_fwd * 1.0e-3))
            if self.reversible:
                print('dGrev: % 1.2f (kJ/mol)' % (self.dGa_rev * 1.0e-3))
                print('dGrxn: % 1.2f (kJ/mol)' % (self.dGrxn * 1.0e-3))
                print('dErev: % 1.2f (kJ/mol)' % (self.dEa_rev * 1.0e-3))
                print('dErxn: % 1.2f (kJ/mol)' % (self.dErxn * 1.0e-3))
            print('---------------------')
