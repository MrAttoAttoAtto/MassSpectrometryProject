from matplotlib import pyplot as plt
import itertools
from typing import List
from operator import mul
from functools import reduce
import numpy as np


class Element:
    def __init__(self, masses: List[int], abundances: List[float]):
        assert len(masses) == len(abundances)

        self.masses = masses
        self.abundances = abundances


class Atom:
    def __init__(self, element: Element):
        self.element = element


class Bond:
    def __init__(self, start: Atom, end: Atom, strength: float):
        self.start = start
        self.end = end
        self.strength = strength/50


def estimate_relative_frequencies(bonds: List[Bond], molecule_name=None, factor: float = 1):
    atoms = []
    masses = []
    abundances = []
    connections = {}

    for bond in bonds:
        start = bond.start
        end = bond.end

        if start not in atoms:
            atoms.append(start)
            masses.append(start.element.masses)
            abundances.append(start.element.abundances)
            connections[start] = []

        connections[start].append(bond)

        if end not in atoms:
            atoms.append(end)
            masses.append(end.element.masses)
            abundances.append(end.element.abundances)
            connections[end] = []

        connections[end].append(bond)

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
        abundance = reduce(mul, [abundance[permutation[j]] for j, abundance in enumerate(abundances)], 1) * factor

        if mzs.get(mass) is None:
            mzs[mass] = abundance
        else:
            mzs[mass] += abundance

    for i in range(len(bonds)):
        if factor == 1:
            print(i)
        old_bond = bonds[i]

        j = 0
        left_molecule = [old_bond]
        checked_atoms = []
        while True:
            if j == len(left_molecule):
                break
            bond_to_check = left_molecule[j]

            atom_to_check = bond_to_check.start
            if atom_to_check not in checked_atoms:
                for bond in connections[atom_to_check]:
                    if bond not in left_molecule:
                        left_molecule.append(bond)

                checked_atoms.append(bond_to_check.start)

            if j != 0:
                atom_to_check = bond_to_check.end
                if atom_to_check not in checked_atoms:
                    for bond in connections[atom_to_check]:
                        if bond not in left_molecule:
                            left_molecule.append(bond)

                    checked_atoms.append(bond_to_check.end)

            j += 1

        del left_molecule[0]

        j = 0
        right_molecule = [old_bond]
        checked_atoms = []
        while True:
            if j == len(right_molecule):
                break
            bond_to_check = right_molecule[j]

            atom_to_check = bond_to_check.end
            if atom_to_check not in checked_atoms:
                for bond in connections[atom_to_check]:
                    if bond not in right_molecule:
                        right_molecule.append(bond)

                checked_atoms.append(bond_to_check.end)

            if j != 0:
                atom_to_check = bond_to_check.start
                if atom_to_check not in checked_atoms:
                    for bond in connections[atom_to_check]:
                        if bond not in right_molecule:
                            right_molecule.append(bond)

                    checked_atoms.append(bond_to_check.start)

            j += 1

        del right_molecule[0]

        left_mzs = {}
        if len(left_molecule) == 0:
            if 1 not in old_bond.start.element.masses:
                left_mzs = dict(zip(old_bond.start.element.masses, np.array(old_bond.start.element.abundances)*factor*0.1))
        else:
            left_mzs = estimate_relative_frequencies(left_molecule, None, factor / old_bond.strength)

        mzs = {key: mzs.get(key, 0) + left_mzs.get(key, 0)
               for key in set(mzs) | set(left_mzs)}

        right_mzs = {}
        if len(right_molecule) == 0:
            if 1 not in old_bond.end.element.masses:
                right_mzs = dict(zip(old_bond.end.element.masses, np.array(old_bond.end.element.abundances)*factor*0.1))
        else:
            right_mzs = estimate_relative_frequencies(right_molecule, None, factor / old_bond.strength)

        mzs = {key: mzs.get(key, 0) + right_mzs.get(key, 0)
                  for key in set(mzs) | set(right_mzs)}

    if factor != 1:
        return mzs

    print(mzs)

    xs = range(1, max(mzs.keys())+1)

    fig, ax = plt.subplots(figsize=(20, 4))

    ys = [0 if mzs.get(i) is None else mzs.get(i) for i in xs]
    ax.bar(xs, ys, width=0.3, color="red")

    last_elevated = False
    for i, y in enumerate(ys):
        if y == 0:
            continue

        if abs(ys[i-1] - y) < 0.1 and not last_elevated:
            y_coord = y + 0.07
            last_elevated = True
        else:
            y_coord = y
            last_elevated = False

        ax.text(i+1, y_coord+max(ys)/35, str(i+1), color='black', fontsize=7, fontweight='bold', ha="center")

    ax.set(ylim=(0, max(ys)*1.125))

    ax.set_ylabel("Relative Intensity")
    ax.set_xlabel("m/z")

    ax.set_title(f"Mass Spectrum{' for ' + molecule_name if molecule_name is not None else ''}")

    plt.show()
    #fig.savefig("C:\\Users\\<username>\\Downloads\\cool.png")


carbon = Element([12, 13], [0.989, 0.011])
chlorine = Element([35, 37], [0.75, 0.25])
hydrogen = Element([1], [1])
oxygen = Element([16], [1])

'''
left_carbon = Atom(carbon)
right_carbon = Atom(carbon)

left_oxygen = Atom(oxygen)
right_oxygen = Atom(oxygen)

bonds = [
    Bond(left_oxygen, Atom(hydrogen)),
    Bond(left_oxygen, left_carbon),
    Bond(left_carbon, Atom(hydrogen)),
    Bond(left_carbon, Atom(hydrogen)),
    Bond(left_carbon, right_carbon),
    Bond(right_oxygen, Atom(hydrogen)),
    Bond(right_oxygen, right_carbon),
    Bond(right_carbon, Atom(hydrogen)),
    Bond(right_carbon, Atom(hydrogen))
]'''

'''
left_carbon = Atom(carbon)
right_carbon = Atom(carbon)

bonds = [
    Bond(left_carbon, Atom(hydrogen), 413),
    Bond(left_carbon, Atom(hydrogen), 413),
    Bond(left_carbon, Atom(hydrogen), 413),
    Bond(left_carbon, right_carbon, 346),
    Bond(right_carbon, Atom(hydrogen), 413),
    Bond(right_carbon, Atom(hydrogen), 413),
    Bond(right_carbon, Atom(hydrogen), 413)
]

estimate_relative_frequencies(bonds, "Ethane")
'''

right_carbon = Atom(carbon)

bonds = [
    Bond(right_carbon, Atom(chlorine), 346),
    Bond(right_carbon, Atom(chlorine), 346),
    Bond(right_carbon, Atom(chlorine), 346),
    Bond(right_carbon, Atom(hydrogen), 413)
]

estimate_relative_frequencies(bonds, "Chloroform")
