from matplotlib import pyplot as plt
import itertools
from typing import List
from operator import mul
from functools import reduce


class Element:
    def __init__(self, masses: List[float], abundances: List[float]):
        assert len(masses) == len(abundances)

        self.masses = masses
        self.abundances = abundances


class Atom:
    def __init__(self, element: Element):
        self.element = element


class Bond:
    def __init__(self, start: Atom, end: Atom):
        self.start = start
        self.end = end


def estimate_relative_frequencies(bonds: List[Bond]):
    atoms = []
    masses = []
    abundances = []

    for bond in bonds:
        start = bond.start
        end = bond.end

        if start not in atoms:
            atoms.append(start)
            masses.append(start.element.masses)
            abundances.append(start.element.abundances)

        if end not in atoms:
            atoms.append(end)
            masses.append(end.element.masses)
            abundances.append(end.element.abundances)

    isotope_permutation_count = reduce(mul, [len(mass) for mass in masses], 1)
    permutation = [0]*len(atoms)

    mzs = {}
    for i in range(isotope_permutation_count):
        if i != 0:
            for j in range(len(permutation)):
                permutation[j] += 1
                if permutation[j] != len(masses[j]):
                    break
                permutation[j] = 0

        mass = sum([mass[permutation[j]] for j, mass in enumerate(masses)])
        abundance = reduce(mul, [abundance[permutation[j]] for j, abundance in enumerate(abundances)], 1)

        if mzs.get(mass) is None:
            mzs[mass] = abundance
        else:
            mzs[mass] += abundance

    xs = range(1, max(mzs.keys())+1)

    fig, ax = plt.subplots()

    ax.bar(xs, [0 if mzs.get(i) is None else mzs.get(i) for i in xs])
    fig.show()


carbon = Element([12, 13], [0.989, 0.011])
chlorine = Element([35, 37], [0.75, 0.25])
hydrogen = Element([1], [1])

middle_carbon = Atom(carbon)

bonds = [
    Bond(middle_carbon, Atom(hydrogen)),
    Bond(middle_carbon, Atom(chlorine)),
    Bond(middle_carbon, Atom(chlorine)),
    Bond(middle_carbon, Atom(chlorine))
]

estimate_relative_frequencies(bonds)
