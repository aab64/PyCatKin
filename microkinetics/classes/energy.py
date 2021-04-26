import copy
import numpy as np
import matplotlib.pyplot as plt
from microkinetics.constants.physical_constants import *


class Energy:
    def __init__(self, minima=None):
        """Initialises Energy class.
        Energy class stores the states involved in the energy landscape,
        and computes energy spans.

        """
        self.minima = copy.copy(minima)
        if self.minima is None:
            print('No states loaded.')
        self.energy_landscape = None

    def construct_energy_landscape(self, T, p, verbose=False):
        """Records free and electronic energies of minima and
        transition states on the energy landscape
        relative to the first entry in minima.
        
        """
        
        self.energy_landscape = dict({'free': {}, 'electronic': {}, 'isTS': {}, 'T': T, 'p': p})

        ref_free = sum([s.get_free_energy(T=T, p=p, verbose=verbose) for s in self.minima[0]])
        ref_elec = sum([s.Gelec for s in self.minima[0]])

        for sind, state in self.minima.items():
            self.energy_landscape['free'][sind] = sum([s.get_free_energy(T=T, p=p, verbose=verbose)
                                                       for s in self.minima[sind]]) - ref_free
            self.energy_landscape['electronic'][sind] = sum([s.Gelec for s in self.minima[sind]]) - ref_elec
            self.energy_landscape['isTS'][sind] = 1 if True in [i.state_type == 'TS' for i in self.minima[sind]] else 0

    def draw_energy_landscape(self, T, p, etype='free', eunits='eV', verbose=False, path=None):
        """Records free and electronic energies of minima and
        transition states on the energy landscape
        relative to the first entry in minima.

        """

        if self.energy_landscape is None:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)
        elif self.energy_landscape['T'] != T or self.energy_landscape['p'] != p:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)

        from scipy.interpolate import CubicSpline

        fmt = '%1.2f'
        if eunits == 'eV':
            conv = 1.0
        elif eunits == 'kcal/mol':
            conv = eVtokcal
            fmt = '%1.1f'
        elif eunits == 'kJ/mol':
            conv = eVtokJ
            fmt = '%1.1f'
        elif eunits == 'J/mol':
            conv = eVtokJ * 1.0e3
            fmt = '%1.2e'
        else:
            print('Specified conversion not possible, using eV')
            conv = 1.0
            eunits = 'eV'

        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        xpoints = []
        ypoints = []
        for i in range(len(self.energy_landscape[etype].keys())):
            toadd = 0.25 if self.energy_landscape['isTS'][i] else 0.25
            if not self.energy_landscape['isTS'][i]:
                xpoints += [i - toadd, i + toadd]
                ypoints += [self.energy_landscape[etype][i] * conv, self.energy_landscape[etype][i] * conv]
            else:
                xs = [i - 1 + toadd, i]
                ys = [self.energy_landscape[etype][i - 1], self.energy_landscape[etype][i]]
                spl = CubicSpline(xs, ys, bc_type='clamped')
                xint = np.linspace(start=xs[0], stop=xs[-1], num=100)
                yint = spl(xint)
                xpoints += [x for x in xint]
                ypoints += [y * conv for y in yint]
                xs = [i, i + 1 - toadd]
                ys = [self.energy_landscape[etype][i], self.energy_landscape[etype][i + 1]]
                spl = CubicSpline(xs, ys, bc_type='clamped')
                xint = np.linspace(start=xs[0], stop=xs[-1], num=100)
                yint = spl(xint)
                xpoints += [x for x in xint]
                ypoints += [y * conv for y in yint]

        ax.plot(xpoints, ypoints, '-', color='grey')
        label_TS = True
        label_I = True
        for k in self.energy_landscape[etype].keys():
            if self.energy_landscape['isTS'][k] == 1:
                ax.plot(k, self.energy_landscape[etype][k] * conv, 's', color='dodgerblue', label=('Transition state' if
                                                                                                   label_TS else ''))
                label_TS = False
            else:
                ax.plot(k, self.energy_landscape[etype][k] * conv, 's', color='tomato', label=('Intermediate' if
                                                                                               label_I else ''))
                label_I = False
            ax.text(k, self.energy_landscape[etype][k] * conv + 0.2 * conv, fmt % (self.energy_landscape[etype][k] *
                                                                                   conv), ha='center')
        ax.legend()
        ax.set(xlabel='Reaction coordinate', ylabel='Relative ' + etype + ' energy (' + eunits + ')',
               xticks=range(len(self.energy_landscape[etype].keys())),
               ylim=(ax.get_ylim()[0] - 0.25 * conv, ax.get_ylim()[1] + 0.25 * conv),
               xlim=(-1, len(self.energy_landscape[etype].keys())))
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        fig.tight_layout()
        if path is not None:
            fig.savefig(path + etype + '_energy_landscape.png', format='png', dpi=300)

    def evaluate_energy_span_model(self, T, p, etype='free', verbose=False, opath=None):
        """Energy span calculations.

        Returns turnover frequency (tof), relative contributions
        of transition states (num_i) and intermediates (num_j)
        and labels of both respectively (lTi, lIj)."""

        if self.energy_landscape is None:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)
        elif self.energy_landscape['T'] != T or self.energy_landscape['p'] != p:
            self.construct_energy_landscape(T=T, p=p, verbose=verbose)

        Ti = [self.energy_landscape[etype][s] for s in self.energy_landscape[etype].keys()
              if self.energy_landscape['isTS'][s] == 1]
        Ij = [self.energy_landscape[etype][s] for s in self.energy_landscape[etype].keys()
              if self.energy_landscape['isTS'][s] == 0]

        nTi = len(Ti)
        nIj = len(Ij)

        drxn = self.energy_landscape[etype][max(list(self.energy_landscape[etype].keys()))] * eVtokJ * 1.0e3
        print('dGrxn = %1.2f eV' % (drxn * 1.0e-3 / eVtokJ))

        XTOFTi = np.zeros((nTi, nIj))
        ctri = 0
        ctrj = 0
        for i in range(len(self.energy_landscape['free'].keys())):
            if self.energy_landscape['isTS'][i]:
                Ti = self.energy_landscape[etype][i] * eVtokJ * 1.0e3
                for j in range(len(self.energy_landscape['free'].keys())):
                    if not self.energy_landscape['isTS'][j]:
                        Ij = self.energy_landscape[etype][j] * eVtokJ * 1.0e3
                        if i >= j:
                            dGij = drxn
                        else:
                            dGij = 0.0
                        XTOFTi[ctri, ctrj] = Ti - Ij - dGij
                        ctrj += 1
                ctri += 1
                ctrj = 0

        den = sum(sum(np.exp(XTOFTi / (R * T))))
        num_i = [sum([(np.exp(vals / (R * T)) / den) for vals in XTOFTi[i, :]]) for i in range(nTi)]
        num_j = [sum([(np.exp(vals / (R * T)) / den) for vals in XTOFTi[:, j]]) for j in range(nIj)]

        tof = (kB * T / h) * np.exp((-drxn / (R * T)) - 1.0) / den

        lTi = [lab for lab in self.energy_landscape['isTS'].keys() if self.energy_landscape['isTS'][lab] == 1]
        lIj = [lab for lab in self.energy_landscape['isTS'].keys() if self.energy_landscape['isTS'][lab] == 0]

        print('dEmax = %1.2f eV' % (np.max(XTOFTi) * 1.0e-3 / eVtokJ))

        if opath:
            with open(opath, 'w') as tfile:
                tfile.write(str(tof) + '\n')
                tfile.write(', '.join([str(i) for i in num_i] + ['\n']))
                tfile.write(', '.join([str(j) for j in num_j] + ['\n']))

        return tof, num_i, num_j, lTi, lIj


# sigma = 0.3 # eV
# nsamples = 100

# for smpl in range(nsamples):
    # print('Noisy sample ' + str(smpl) + '...\n')
    # g = np.random.normal(loc=0.0, scale=sigma, size=None)
    # pert_energy_dict = copy.copy(energy_dict)
    # for s in pert_energy_dict['Free'].keys():
        # if s not in ['reference', 'gas', 'surface']:
            # if pert_energy_dict['isTS'][s]:
                # u = np.random.uniform()
                # add = u * g
            # else:
                # add = g
            # pert_energy_dict['Free'][s] += np.random.normal(loc=0.0, scale=sigma, size=None)
