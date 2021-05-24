from microkinetics.constants.physical_constants import *
import numpy as np


def karr(T, prefac, barrier):
    """Calculates reaction rate constant from Arrhenius expression.

    Returns rate constant."""

    k = prefac * np.exp(-barrier / (R * T))

    return k


def kads(T, mass, area):
    """Calculates adsorption rate constant from collision theory.

    Returns rate constant."""

    k = area / np.sqrt(2.0 * np.pi * (mass * amutokg) * kB * T)

    return k


def kdes(T, mass, area, sigma, theta, des_en):
    """Calculates desorption rate constant from collision theory.

    Returns rate constant."""

    k = ((kB ** 2) * area * 2.0 * np.pi * (mass * amutokg) * (T ** 3)) / (
            (h ** 3) * sigma * theta) * np.exp(-des_en / (R * T))

    return k


def keq_kin(ka, kd):
    """Calculates equilibrium constant from kinetics.

    Returns equilibrium rate constant."""

    k = ka / kd

    return k


def keq_therm(T, rxn_en):
    """Calculates equilibrium constant from thermodynamics.

    Returns equilibrium rate constant."""

    k = np.exp(-rxn_en / (R * T))

    return k


def k_from_eq_rel(kknown, Keq, direction='forward'):
    """Calculates unknown forward/reverse rate constant from equilibrium relation.

    Returns unknown rate constant."""

    if direction == 'forward':
        k = kknown / Keq
    else:
        k = kknown * Keq

    return k


def prefactor(T):
    """Calculates prefactor from transition state theory.

    Returns prefactor."""

    prefac = kB * T / h

    return prefac
