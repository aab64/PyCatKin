import numpy as np
import ase.io
import os
import pickle
from microkinetics.constants.physical_constants import *
from ase.visualize import view


class State:

    def __init__(self, state_type, path=None, vibs_path=None, name=None, sigma=None, mass=None, gasdata=None):
        """Initialises State class.
        State class stores the species name and atoms object,
        the electronic energy and vibrational frequencies from DFT,
        the degrees of freedom contributions to the free energy.

        """
        if name is None:
            name = os.path.basename(path)
        self.state_type = state_type
        self.path = path
        self.name = name
        self.vibs_path = vibs_path
        self.mass = mass
        self.gasdata = gasdata
        self.atoms = None
        self.freq = None
        self.i_freq = None
        self.Gelec = None
        self.Gtran = None
        self.Gvibr = None
        self.Grota = None
        self.Gfree = None
        self.Gzpe = None
        self.inertia = None
        self.shape = None
        self.sigma = sigma
        if self.state_type == 'gas':
            assert(self.sigma is not None)

    def get_atoms(self):
        """Reads OUTCAR file from path and extracts atoms object.

        """
        assert(self.path is not None)
        outcar_path = self.path + '/OUTCAR'
        assert(os.path.isfile(outcar_path))
        self.atoms = ase.io.read(outcar_path, format='vasp-out')
        self.mass = sum(self.atoms.get_masses())
        if self.state_type == 'gas':
            self.inertia = self.atoms.get_moments_of_inertia()
            # Truncate inertial components likely resulting from low precision
            inertia_cutoff = 1.0e-12
            self.inertia = np.array([i if i > inertia_cutoff else 0.0 for i in self.inertia])
            self.shape = len([i for i in self.inertia if i > 0.0])
            # Check if there are at least 2 non-zero components
            if self.shape < 2:
                print('Too many components of the moments of inertia are zero. Trying geometry in CONTCAR...')
                outcar_path = self.path + '/CONTCAR'
                assert(os.path.isfile(outcar_path))
                atoms = ase.io.read(outcar_path)
                self.inertia = atoms.get_moments_of_inertia()
                self.inertia = np.array([i if i > inertia_cutoff else 0.0 for i in self.inertia])
                self.shape = len([i for i in self.inertia if i > 0.0])
                assert(self.shape > 1)

    def get_vibrations(self, verbose=False):
        """Reads vibrations file from path and extracts frequencies.

        """
        if self.vibs_path is not None:
            freq_path = self.vibs_path + '/log.vib'
        else:
            assert(self.path is not None)
            freq_path = self.path + '/log.vib'

        if os.path.isfile(freq_path):
            with open(freq_path, 'r') as fd:
                lines = fd.readlines()
            initat = 0
            endat = 0
            for lind, line in enumerate(lines):
                if '#' in line:
                    initat = lind + 2
                    endat = 0
                if lind > initat and not endat and '---' in line:
                    endat = lind - 1
            freq = [(float(line.strip().split()[1]) * 1e-3 / (h * JtoeV)) for line in
                    lines[initat:endat + 1] if 'i' not in line]
            i_freq = [(float(line.strip().split()[1].split('i')[0]) * 1e-3 / (h * JtoeV)) for line in
                      lines[initat:endat + 1] if 'i' in line]
        else:
            if verbose:
                print('Checking OUTCAR for frequencies')
            assert(self.path is not None)
            freq_path = self.path + '/OUTCAR'
            assert(os.path.isfile(freq_path))
            freq = []
            i_freq = []
            firstcopy = None
            with open(freq_path, 'r') as fd:
                lines = fd.readlines()
            for line in lines:
                data = line.split()
                if 'THz' in data:
                    if firstcopy != data[0]:
                        fHz = float(data[-8]) * 1e12
                        if 'f/i=' not in data:
                            freq.append(fHz)
                        else:
                            i_freq.append(fHz)
                        if firstcopy is None:
                            firstcopy = data[0]
                    else:
                        break
        # Truncate small freqs
        for f in range(len(freq)):
            if (freq[f] * h * JtoeV * 1e3) < 12.4:
                freq[f] = 12.4 * 1e-3 / (h * JtoeV)
                if verbose:
                    print('Truncating small freq %1.2f to 12.4 meV' % (freq[f] * h * JtoeV * 1e3))
        # Check correct DOF
        n_freq = len(freq)
        n_dof = len(freq) + len(i_freq)  # 3 * N_moving_atoms
        if self.state_type == 'gas':
            n_dof -= 3  # Translational DOF
        if n_freq < n_dof:
            if verbose:
                print('Incorrect number of frequencies! n_dof = %1.0f and n_freq = %1.0f' % (n_dof, n_freq))
                print('Adding %1.0f extra frequencies of 12.4 meV...' % (n_dof - n_freq))
            freq += [12.4 * 1e-3 / (h * JtoeV) for f in range(n_dof - n_freq)]
        self.freq = np.array(sorted(freq, reverse=True))
        self.i_freq = np.array(i_freq)

    def calc_electronic_energy(self, verbose=False):
        """Calculates electronic energy.

        Saves value in eV."""
        if self.atoms is None:
            self.get_atoms()
        self.Gelec = self.atoms.get_potential_energy(force_consistent=True)

    def calc_zpe(self, verbose=False):
        """Calculates zero point energy.

        Saves value in eV."""
        if self.freq is None:
            self.get_vibrations(verbose=verbose)
        self.Gzpe = 0.5 * h * sum(self.freq) * JtoeV

    def calc_vibrational_contrib(self, T, verbose=False):
        """Calculates vibrational contribution to free energy.

        Saves value in eV."""
        if self.freq is None:
            self.get_vibrations(verbose=verbose)
        if self.atoms is None:
            self.get_atoms()

        # Truncate some modes if required
        if self.state_type == 'gas':
            ntrunc = self.shape
        else:
            ntrunc = 0
        nfreqs = self.freq.shape[0] - ntrunc
        use_freq = self.freq[0:nfreqs]

        self.Gvibr = (0.5 * h * sum(use_freq) + kB * T * sum(np.log(1 - np.exp(-use_freq * h / (kB * T))))) * JtoeV

    def calc_translational_contrib(self, T, p, verbose=False):
        """Calculates translational contribution to free energy.

        Saves value in eV."""
        if self.atoms is None:
            self.get_atoms()
        if self.state_type == 'gas':
            self.Gtran = (-kB * T * np.log((kB * T / p) * pow(2 * np.pi * (self.mass * amutokg) * kB * T / (h ** 2),
                                                              1.5))) * JtoeV
        else:
            self.Gtran = 0.0
        if self.gasdata is not None:
            for s in range(len(self.gasdata['fraction'])):
                self.gasdata['state'][s].calc_translational_contrib(T=T, p=p, verbose=verbose)
                self.Gtran += self.gasdata['fraction'][s] * self.gasdata['state'][s].Gtran

    def calc_rotational_contrib(self, T, verbose=False):
        """Calculates rotational contribution to free energy
        accounting for linear/non-linear molecule.

        Saves value in eV."""
        if self.state_type == 'gas':
            I = self.inertia * amuA2tokgm2
            if self.shape == 2:
                I = np.sqrt(np.prod([I[i] for i in range(len(I)) if I[i] != 0]))
                self.Grota = (-kB * T * np.log(8 * np.pi * np.pi * kB * T * I / (self.sigma * h ** 2))) * JtoeV
            else:
                self.Grota = (-kB * T * np.log((np.sqrt(np.pi) / self.sigma) * pow(
                    8 * np.pi * np.pi * kB * T / (h ** 2), 1.5) * np.sqrt(np.prod(I)))) * JtoeV
        else:
            self.Grota = 0.0
        if self.gasdata is not None:
            for s in range(len(self.gasdata['fraction'])):
                self.gasdata['state'][s].calc_rotational_contrib(T=T, verbose=verbose)
                self.Grota += self.gasdata['fraction'][s] * self.gasdata['state'][s].Grota

    def calc_free_energy(self, T, p, verbose=False):
        """Calculates free energy. 

        Saves value in eV."""
        self.calc_electronic_energy()
        self.calc_vibrational_contrib(T=T, verbose=verbose)
        self.calc_translational_contrib(T=T, p=p, verbose=verbose)
        self.calc_rotational_contrib(T=T, verbose=verbose)

        self.Gfree = self.Gelec + self.Gtran + self.Grota + self.Gvibr

        if verbose:
            print((self.name + ': %1.2f eV') % self.Gfree)

    def get_free_energy(self, T, p, verbose=False):
        """Returns the free energy in eV.

        """
        self.calc_free_energy(T=T, p=p, verbose=verbose)
        return self.Gfree

    def save_pdb(self, path='Images/xyz/'):
        """Saves the atoms object as a pdb structure file.

        """
        if self.atoms is None:
            self.get_atoms()
        assert(os.path.isdir(path))
        ase.io.write(path + self.name + '.pdb', self.atoms, format='proteindatabank')

    def save_pickle(self, path=None):
        """Save the state as a pickle object.

        """
        path = '' if path is None else path
        pickle.dump(self, open(path + self.name + '.pckl', 'wb'))

    def view_atoms(self, path=None, rotation=None):
        """Views the atoms object using the ASE visualizer.
        If path is not None, saves as a png.

        """
        if self.atoms is None:
            self.get_atoms()
        view(self.atoms)
        if path is not None:
            ase.io.write(path + self.name + '.png', self.atoms, rotation=rotation)
