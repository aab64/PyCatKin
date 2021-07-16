from pycatkin.constants.physical_constants import *
import copy
import pickle
import os
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


class Energy:

    def __init__(self, name='landscape', minima=None, labels=None, path_to_pickle=None):
        """Initialises Energy class.
        Energy class stores the states in the energy landscape,
        and computes energy span model predictions.
        If path_to_pickle is defined, the pickled object is loaded.

        """

        if path_to_pickle:
            assert (os.path.isfile(path_to_pickle))
            newself = pickle.load(open(path_to_pickle, 'rb'))
            assert (isinstance(newself, Energy))
            for att in newself.__dict__.keys():
                setattr(self, att, getattr(newself, att))
        else:
            self.name = name
            self.minima = minima
            if labels is not None:
                self.labels = labels
            else:
                self.labels = [i[0].name for i in minima]
            self.energy_landscape = None
            if self.minima is None:
                print('No states loaded.')
            if self.labels is not None:
                assert(len(self.labels) == len(self.minima))

    def construct_energy_landscape(self, T, p, verbose=False):
        """Records free and electronic energies of minima and
        transition states on the energy landscape
        relative to the first entry in minima.
        
        """

        self.energy_landscape = dict({'free': {},
                                      'electronic': {},
                                      'isTS': {},
                                      'T': T,
                                      'p': p})

        ref_free = sum([s.get_free_energy(T=T, p=p, verbose=verbose) for s in self.minima[0]])
        ref_elec = sum([s.Gelec for s in self.minima[0]])

        for sind, state in enumerate(self.minima):
            self.energy_landscape['free'][sind] = sum([s.get_free_energy(T=T, p=p, verbose=verbose)
                                                       for s in self.minima[sind]]) - ref_free
            self.energy_landscape['electronic'][sind] = sum([s.Gelec for s in self.minima[sind]]) - ref_elec
            self.energy_landscape['isTS'][sind] = 1 if True in [i.state_type == 'TS'
                                                                for i in self.minima[sind]] else 0

    def draw_energy_landscape(self, T, p, etype='free', eunits='eV', legend_location='upper right',
                              verbose=False, path=None, show_labels=False):
        """Records free and electronic energies of minima and
        transition states on the energy landscape
        relative to the first entry in minima.

        """

        if self.energy_landscape is None:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)
        elif self.energy_landscape['T'] != T or self.energy_landscape['p'] != p:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)

        if show_labels:
            assert(self.labels is not None)

        fmt = '%.3g'
        if eunits == 'eV':
            conv = 1.0
        elif eunits == 'kcal/mol':
            conv = eVtokcal
        elif eunits == 'kJ/mol':
            conv = eVtokJ
        elif eunits == 'J/mol':
            conv = eVtokJ * 1.0e3
        else:
            print('Specified conversion not possible, using eV')
            conv = 1.0
            eunits = 'eV'

        fig, ax = plt.subplots(figsize=(10, 4))
        xpoints = []
        ypoints = []
        for i in range(len(self.energy_landscape[etype].keys())):
            toadd = 0.25 if self.energy_landscape['isTS'][i] else 0.25
            if not self.energy_landscape['isTS'][i]:
                xpoints += [i - toadd,
                            i + toadd]
                ypoints += [self.energy_landscape[etype][i] * conv,
                            self.energy_landscape[etype][i] * conv]
            else:
                xs = [i - 1 + toadd,
                      i]
                ys = [self.energy_landscape[etype][i - 1],
                      self.energy_landscape[etype][i]]
                spl = CubicSpline(xs, ys, bc_type='clamped')
                xint = np.linspace(start=xs[0], stop=xs[-1], num=100)
                yint = spl(xint)
                xpoints += [x for x in xint]
                ypoints += [y * conv for y in yint]
                xs = [i,
                      i + 1 - toadd]
                ys = [self.energy_landscape[etype][i],
                      self.energy_landscape[etype][i + 1]]
                spl = CubicSpline(xs, ys, bc_type='clamped')
                xint = np.linspace(start=xs[0], stop=xs[-1], num=100)
                yint = spl(xint)
                xpoints += [x for x in xint]
                ypoints += [y * conv for y in yint]

        ax.plot(xpoints, ypoints, '-', color='black')
        label_TS = True
        label_I = True
        for k in self.energy_landscape[etype].keys():
            if self.energy_landscape['isTS'][k] == 1:
                ax.plot(k, self.energy_landscape[etype][k] * conv, 's',
                        label=('Transition state' if label_TS else ''),
                        color='tomato')
                label_TS = False
            else:
                ax.plot(k, self.energy_landscape[etype][k] * conv, 's',
                        label=('Intermediate' if label_I else ''),
                        color='darkturquoise')
                label_I = False
            ax.text(k, self.energy_landscape[etype][k] * conv + 0.2 * conv,
                    fmt % (self.energy_landscape[etype][k] * conv),
                    ha='center')
            if show_labels:
                ax.text(k, self.energy_landscape[etype][k] * conv - 0.2 * conv,
                        self.labels[k],
                        ha='center', va='top')
        ax.legend(loc=legend_location)
        ax.set(xlabel='Reaction coordinate',
               xlim=(-1, len(self.energy_landscape[etype].keys())),
               xticks=range(len(self.energy_landscape[etype].keys())),
               ylabel='Relative ' + etype + ' energy (' + eunits + ')',
               ylim=(ax.get_ylim()[0] - 0.25 * conv, ax.get_ylim()[1] + 0.25 * conv))
        plt.tick_params(axis='x', which='both',
                        bottom=False, top=False, labelbottom=False)
        fig.tight_layout()
        if path:
            fig.savefig(path + etype + '_energy_landscape.png',
                        format='png', dpi=300)

    def evaluate_energy_span_model(self, T, p, etype='free', verbose=False, opath=None):
        """Energy span calculations.

        Returns turnover frequency (tof), relative contributions
        of transition states (num_i) and intermediates (num_j)
        and labels of both respectively (lTi, lIj)."""

        if self.energy_landscape is None:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)
        elif self.energy_landscape['T'] != T or self.energy_landscape['p'] != p:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)

        nTi = len([self.energy_landscape[etype][s]
                   for s in self.energy_landscape[etype].keys()
                   if self.energy_landscape['isTS'][s] == 1])
        nIj = len([self.energy_landscape[etype][s]
                   for s in self.energy_landscape[etype].keys()
                   if self.energy_landscape['isTS'][s] == 0]) - 1

        drxn = self.energy_landscape[etype][max(list(self.energy_landscape[etype].keys()))] * eVtokJ * 1.0e3
        print('dGrxn = %1.2f eV' % (drxn * 1.0e-3 / eVtokJ))

        XTOFTi = np.zeros((nTi, nIj))
        ctri = 0
        ctrj = 0
        for i in range(nTi + nIj):
            if self.energy_landscape['isTS'][i]:
                Ti = self.energy_landscape[etype][i] * eVtokJ * 1.0e3
                for j in range(nTi + nIj):
                    if not self.energy_landscape['isTS'][j]:
                        Ij = self.energy_landscape[etype][j] * eVtokJ * 1.0e3
                        dGij = drxn if i >= j else 0.0
                        XTOFTi[ctri, ctrj] = Ti - Ij - dGij
                        ctrj += 1
                ctri += 1
                ctrj = 0

        den = sum(sum(np.exp(XTOFTi / (R * T))))
        num_i = [sum([(np.exp(vals / (R * T)) / den)
                      for vals in XTOFTi[i, :]])
                 for i in range(nTi)]
        num_j = [sum([(np.exp(vals / (R * T)) / den)
                      for vals in XTOFTi[:, j]])
                 for j in range(nIj)]

        iTDTS = [i for i in range(len(num_i)) if num_i[i] == max(num_i)][0]
        iTDTS = [k for k in self.energy_landscape['isTS'].keys()
                 if self.energy_landscape['isTS'][k] == 1][iTDTS]
        iTDI = [j for j in range(len(num_j)) if num_j[j] == max(num_j)][0]
        iTDI = [k for k in self.energy_landscape['isTS'].keys()
                 if self.energy_landscape['isTS'][k] == 0][iTDI]

        TDTS = self.labels[iTDTS]
        TDI = self.labels[iTDI]

        tof = (kB * T / h) * np.exp((-drxn / (R * T)) - 1.0) / den

        lTi = [self.labels[lab] for lab in self.energy_landscape['isTS'].keys()
               if self.energy_landscape['isTS'][lab] == 1]
        lIj = [self.labels[lab] for lab in self.energy_landscape['isTS'].keys()
               if self.energy_landscape['isTS'][lab] == 0][0:-1]

        Espan = self.energy_landscape[etype][iTDTS] - self.energy_landscape[etype][iTDI]
        Eapp = np.log((h * tof) / (kB * T)) * (-R * T) * 1.0e-3
        print('Energy span model results (%1.0f K): ' % T)
        print('* TOF = % .3g 1/s' % tof)
        print('* Espan = %.3g eV = %.3g kcal/mol = %.3g kJ/mol' %
              (Espan, Espan * eVtokcal, Espan * eVtokJ))
        print('* TDTS is %s.' % TDTS)
        print('* TDI is %s.' % TDI)
        print('* dGrxn = %.3g eV = %.3g kcal/mol = %.3g kJ/mol' %
              (drxn * 1.0e-3 / eVtokJ, drxn / kcaltoJ, drxn * 1.0e-3))
        print('* Eapp = %.3g eV = %.3g kcal/mol = %.3g kJ/mol' %
              (Eapp / eVtokJ, Eapp * 1.0e3 / kcaltoJ, Eapp))

        if opath:
            with open(opath, 'w') as tfile:
                tfile.write(str(tof) + '\n')
                tfile.write(', '.join([str(i) for i in num_i] + ['\n']))
                tfile.write(', '.join([str(j) for j in num_j] + ['\n']))

        return tof, Espan, TDTS, TDI, num_i, num_j, lTi, lIj

    def save_pickle(self, path=None):
        """Save the energy landscape as a pickle object.

        """

        path = path if path else ''
        pickle.dump(self, open(path + 'energy_' + self.name + '.pckl', 'wb'))
