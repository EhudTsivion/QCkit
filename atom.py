# -*- coding: utf-8 -*-

# This is a free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ASE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

import numpy as np
import QCkit.physical_constants as constants


class Atom(object):
    def __init__(self, sym=None, coords=[], coord_units="bohr"):
        """

        default constructor of atom object

        :param sym: chemical symbol of atom
        :param coords: cartesian coordinates of the atom, in units
        :param coord_units: the units in which coordinates are provided, can
                be bohr or angstroms, angs, Angs etc.. anthing with first letter 'A'
        """

        # make sure chemical symbol has the correct
        # letter-case form
        sym = correct_symbol_case(sym)
        self.symbol = str(sym)

        # set mass
        self.mass = mass_for_sym(sym)

        # verify correctness of input
        if len(coords) == 3:

            # if units are given in Angstroms, convert to bohr
            if coord_units[0].lower() == "a":
                self._xyz = np.array(coords).astype(np.float64) * constants.angstrom_to_bohr

            else:
                self._xyz = np.array(coords).astype(np.float64)

        else:
            raise ValueError("error reading the coordinates: \"{}\"".format(coords))

    @property
    def coords(self):
        return self._xyz

    @coords.setter
    def coords(self, new_coord):
        self._xyz = np.copy(new_coord)

    @property
    def atomic_number(self):
        return number_for_sym(self.symbol)

    def covalent_radius(self, bond_order=1):
        return covalent_radii[self.atomic_number - 1, bond_order - 1]

    def __str__(self):
        return "{} {} {:< 10.8f} {:< 10.8f} {:< 10.8f}".format(self.symbol, self.mass, self._xyz[0], self._xyz[1],
                                                               self._xyz[2])


chemical_symbols = ['X', 'H', 'He', 'Li', 'Be',
                    'B', 'C', 'N', 'O', 'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P', 'S', 'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V', 'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I', 'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U', 'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

# this should have the same size as "chemical_symbols" list
atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium']

# this should have the same size as "chemical_symbols" list
# the atomic masses are given in atomic mass units
atomic_masses = np.array([
    0.00000,  # X
    1.00794,  # H
    4.00260,  # He
    6.94100,  # Li
    9.01218,  # Be
    10.81100,  # B
    12.01100,  # C
    14.00670,  # N
    15.99940,  # O
    18.99840,  # F
    20.17970,  # Ne
    22.98977,  # Na
    24.30500,  # Mg
    26.98154,  # Al
    28.08550,  # Si
    30.97376,  # P
    32.06600,  # S
    35.45270,  # Cl
    39.94800,  # Ar
    39.09830,  # K
    40.07800,  # Ca
    44.95590,  # Sc
    47.88000,  # Ti
    50.94150,  # V
    51.99600,  # Cr
    54.93800,  # Mn
    55.84700,  # Fe
    58.93320,  # Co
    58.69340,  # Ni
    63.54600,  # Cu
    65.39000,  # Zn
    69.72300,  # Ga
    72.61000,  # Ge
    74.92160,  # As
    78.96000,  # Se
    79.90400,  # Br
    83.80000,  # Kr
    85.46780,  # Rb
    87.62000,  # Sr
    88.90590,  # Y
    91.22400,  # Zr
    92.90640,  # Nb
    95.94000,  # Mo
    np.nan,  # Tc
    101.07000,  # Ru
    102.90550,  # Rh
    106.42000,  # Pd
    107.86800,  # Ag
    112.41000,  # Cd
    114.82000,  # In
    118.71000,  # Sn
    121.75700,  # Sb
    127.60000,  # Te
    126.90450,  # I
    131.29000,  # Xe
    132.90540,  # Cs
    137.33000,  # Ba
    138.90550,  # La
    140.12000,  # Ce
    140.90770,  # Pr
    144.24000,  # Nd
    np.nan,  # Pm
    150.36000,  # Sm
    151.96500,  # Eu
    157.25000,  # Gd
    158.92530,  # Tb
    162.50000,  # Dy
    164.93030,  # Ho
    167.26000,  # Er
    168.93420,  # Tm
    173.04000,  # Yb
    174.96700,  # Lu
    178.49000,  # Hf
    180.94790,  # Ta
    183.85000,  # W
    186.20700,  # Re
    190.20000,  # Os
    192.22000,  # Ir
    195.08000,  # Pt
    196.96650,  # Au
    200.59000,  # Hg
    204.38300,  # Tl
    207.20000,  # Pb
    208.98040,  # Bi
    np.nan,  # Po
    np.nan,  # At
    np.nan,  # Rn
    np.nan,  # Fr
    226.02540,  # Ra
    np.nan,  # Ac
    232.03810,  # Th
    231.03590,  # Pa
    238.02900,  # U
    237.04820,  # Np
    np.nan,  # Pu
    np.nan,  # Am
    np.nan,  # Cm
    np.nan,  # Bk
    np.nan,  # Cf
    np.nan,  # Es
    np.nan,  # Fm
    np.nan,  # Md
    np.nan,  # No
    np.nan]).astype(np.float)  # Lw

# an array of covalent radii of atoms
#
# 1st column -> single bond
# 2nd column -> double bond
# 3rd column -> triple bond
# 
# adapted from http://en.wikipedia.org/w/index.php?title=Covalent_radius&oldid=623288198
#
# based on
# P. Pyykkö, M. Atsumi (2009). "Molecular Single-Bond Covalent Radii for Elements 1-118".
# Chemistry: A European Journal 15: 186–197. doi:10.1002/chem.200800987
# 
# P. Pyykkö, M. Atsumi (2009). "Molecular Double-Bond Covalent Radii for Elements Li–E112". 
# Chemistry: A European Journal 15 (46): 12770–12779. doi:10.1002/chem.200901472
#
# P. Pyykkö, S. Riedel, M. Patzschke (2005). "Triple-Bond Covalent Radii". 
# Chemistry: A European Journal 11 (12): 3511–3520. doi:10.1002/chem.200401299
#
# Bohr units are used
covalent_radii = np.array([
    [0.605, np.nan, np.nan],  # H
    [0.869, np.nan, np.nan],  # He
    [2.513, 2.343, np.nan],  # Li
    [1.928, 1.701, 1.606],  # Be
    [1.606, 1.474, 1.379],  # B
    [1.417, 1.266, 1.134],  # C
    [1.342, 1.134, 1.020],  # N
    [1.191, 1.077, 1.002],  # O
    [1.209, 1.115, 1.002],  # F
    [1.266, 1.814, np.nan],  # Ne
    [2.929, 3.024, np.nan],  # Na
    [2.627, 2.494, 2.400],  # Mg
    [2.381, 2.135, 2.098],  # Al
    [2.192, 2.022, 1.928],  # Si
    [2.098, 1.928, 1.776],  # P
    [1.946, 1.776, 1.795],  # S
    [1.871, 1.795, 1.757],  # Cl
    [1.814, 2.022, 1.814],  # Ar
    [3.704, 3.647, np.nan],  # K
    [3.231, 2.778, 2.513],  # Ca
    [2.797, 2.192, 2.154],  # Sc
    [2.570, 2.211, 2.041],  # Ti
    [2.532, 2.116, 2.003],  # V
    [2.305, 2.098, 1.946],  # Cr
    [2.249, 1.984, 1.946],  # Mn
    [2.192, 2.060, 1.928],  # Fe
    [2.098, 1.946, 1.814],  # Co
    [2.079, 1.909, 1.909],  # Ni
    [2.116, 2.173, 2.268],  # Cu
    [2.230, 2.268, np.nan],  # Zn
    [2.343, 2.211, 2.287],  # Ga
    [2.287, 2.098, 2.154],  # Ge
    [2.287, 2.154, 2.003],  # As
    [2.192, 2.022, 2.022],  # Se
    [2.154, 2.060, 2.079],  # Br
    [2.211, 2.287, 2.041],  # Kr
    [3.968, 3.817, np.nan],  # Rb
    [3.496, 2.967, 2.627],  # Sr
    [3.080, 2.457, 2.343],  # Y
    [2.910, 2.400, 2.287],  # Zr
    [2.778, 2.362, 2.192],  # Nb
    [2.608, 2.287, 2.135],  # Mo
    [2.419, 2.268, 2.079],  # Tc
    [2.362, 2.154, 1.946],  # Ru
    [2.362, 2.079, 2.003],  # Rh
    [2.268, 2.211, 2.116],  # Pd
    [2.419, 2.627, 2.589],  # Ag
    [2.570, 2.721, np.nan],  # Cd
    [2.683, 2.570, 2.759],  # In
    [2.646, 2.457, 2.494],  # Sn
    [2.646, 2.513, 2.400],  # Sb
    [2.570, 2.419, 2.287],  # Te
    [2.513, 2.438, 2.362],  # I
    [2.476, 2.551, 2.305],  # Xe
    [4.384, 3.950, np.nan],  # Cs
    [3.704, 3.042, 2.816],  # Ba
    [3.402, 2.627, 2.627],  # La
    [3.080, 2.589, 2.476],  # Ce
    [3.326, 2.608, 2.419],  # Pr
    [3.288, 2.589, np.nan],  # Nd
    [3.269, 2.551, np.nan],  # Pm
    [3.250, 2.532, np.nan],  # Sm
    [3.175, 2.532, np.nan],  # Eu
    [3.194, 2.551, 2.494],  # Gd
    [3.175, 2.551, np.nan],  # Tb
    [3.156, 2.513, np.nan],  # Dy
    [3.137, 2.513, np.nan],  # Ho
    [3.118, 2.513, np.nan],  # Er
    [3.099, 2.476, np.nan],  # Tm
    [3.213, 2.438, np.nan],  # Yb
    [3.061, 2.476, 2.476],  # Lu
    [2.872, 2.419, 2.305],  # Hf
    [2.759, 2.381, 2.249],  # Ta
    [2.589, 2.268, 2.173],  # W
    [2.476, 2.249, 2.079],  # Re
    [2.438, 2.192, 2.060],  # Os
    [2.305, 2.173, 2.022],  # Ir
    [2.324, 2.116, 2.079],  # Pt
    [2.343, 2.287, 2.324],  # Au
    [2.513, 2.683, np.nan],  # Hg
    [2.721, 2.683, 2.835],  # Tl
    [2.721, 2.551, 2.589],  # Pb
    [2.853, 2.665, 2.551],  # Bi
    [2.740, 2.551, 2.438],  # Po
    [2.778, 2.608, 2.608],  # At
    [2.683, 2.740, 2.513],  # Rn
    [4.214, 4.120, np.nan],  # Fr
    [3.798, 3.269, 3.005],  # Ra
    [3.515, 2.891, 2.646],  # Ac
    [3.307, 2.702, 2.570],  # Th
    [3.194, 2.608, 2.438],  # Pa
    [3.213, 2.532, 2.230],  # U
    [3.231, 2.570, 2.192],  # Np
    [3.250, 2.551, np.nan],  # Pu
    [3.137, 2.551, np.nan],  # Am
    [3.137, 2.570, np.nan],  # Cm
    [3.175, 2.627, np.nan],  # Bk
    [3.175, 2.646, np.nan],  # Cf
    [3.118, 2.646, np.nan],  # Es
    [3.156, np.nan, np.nan],  # Fm
    [3.269, 2.627, np.nan],  # Md
    [3.326, np.nan, np.nan],  # No
    [3.042, 2.665, np.nan],  # Lr
    [2.967, 2.646, 2.476],  # Rf
    [2.816, 2.570, 2.381],  # Db
    [2.702, 2.419, 2.287],  # Sg
    [2.665, 2.419, 2.249],  # Bh
    [2.532, 2.362, 2.230],  # Hs
    [2.438, 2.362, 2.135],  # Mt
    [2.419, 2.192, 2.116],  # Ds
    [2.287, 2.192, 2.230],  # Rg
    [2.305, 2.589, 2.457],  # Cn
    [2.570, np.nan, np.nan],  # Uut
    [2.702, np.nan, np.nan],  # Fl
    [3.061, np.nan, np.nan],  # Uup
    [3.307, np.nan, np.nan],  # Lv
    [3.118, np.nan, np.nan],  # Uus
    [2.967, np.nan, np.nan]]).astype(np.float)  # Uuo


def correct_symbol_case(chem_symbol):
    """
    correct the chemical symbol to be in the form:

    Upperlower
    example:

    correct_symbol_case('he') -> 'He'
    correct_symbol_case('HE') -> 'He'
    correct_symbol_case('He') -> 'He'

    correct_symbol_case('C') -> 'C'
    correct_symbol_case('c') -> 'c'


    :param chem_symbol: the chemical symbole to correct
    :return: a correct case chemical symbol
    """

    # make sure chemical symbol is not longer than 2.
    sym_length = len(chem_symbol)

    if sym_length > 2:
        raise Exception('Found error while reading the chemical symbol \"{}\"'.format(chem_symbol))

    elif sym_length == 2:
        return chem_symbol[0].upper() + chem_symbol[1].lower()

    elif sym_length == 1:
        return chem_symbol.upper()


def mass_for_sym(symbol):
    """

    :param symbol:
    :return: the atomic mass of the chemical symbol
    """

    if len(chemical_symbols) != len(atomic_masses):
        raise Exception('\"chem_symbol\" list and \"atomic_masses\" numpy array are not '
                        'the same size.')

    try:
        index = chemical_symbols.index(symbol)

    except ValueError():
        print('Chemical symbol \"{}\" does not exist'.format(symbol))
        return 0

    mass = constants.dict_of_atomic_masses[symbol]

    return mass


def number_for_sym(symbol):
    """

    :param symbol:
    :return: the atomic number of the chemical symbol
    """

    try:
        number = chemical_symbols.index(symbol)

    except ValueError:
        print("{} is not a symbol of an atom".format(symbol))

    return number
