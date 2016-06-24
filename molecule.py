import numpy as np
import subprocess as sp
import tempfile
import copy
import os
import logging
import shutil

import QCkit.physical_constants as const
from QCkit.atom import Atom

log = logging.getLogger()


class Molecule(object):
    """

        A class for holding the information and operations
         related to collections of atoms, such as molecule,
         dimer or any complex.

    """

    """
    some notes to the programmer, about atom numbering:

    like everything in python, the atoms are numbered starting 0
    however, to make things more friendly to the user, he refers
    to the atoms as starting from 1. The numbers entered by the user
    can then be converted to the usual python convention.

    for instance, if there are 15 atoms, and the user wishes to split
    the molecule after the 10'th atom, he would use:

    molecule.split(10) (1 .. numbering)

    and not molecule.split(9) (0 .. numbering)
    """

    def __init__(self, atoms=None):

        # all the *_list variables are managed by getter and setters
        if type(atoms) == list or not atoms:
            self._atoms = atoms
        else:
            raise TypeError('Wrong error type when setting molecule')

        self.charge = 0
        self.multiplicity = 1

    def add_atom(self, new_atom):
        """

        :param new_atom:
        :return: adds the atom to the atom list of the molecule
         """

        if not self._atoms:
            self._atoms = []

        self._atoms.append(copy.deepcopy(new_atom))

    def rotate(self, atom_num1, atom_num2, angle=0, atom_list=None, radians=False):
        """

        counter-clockwise rotation of the molecule around the
        axis formed by connecting atom1 and atom2

        :param angle: how much to rotate around the axis formed
        by atoms 1 and 2. Degrees is the defaults

        :param radians: False is angle is given in degrees, True if
        given in radians

        :param atom_list: a list of atoms to be rotated

        :return: Changes the internal state (more specifically: the location)
        of the location of the atoms that make the molecule
        """

        if not atom_list:
            raise ValueError("No atoms to be rotated are found")

        else:
            # correct user input to start from 0 (not user 1..)
            atom_list = np.array(atom_list, dtype=int)
            atom_list -= 1

        max_num = max(atom_num1, atom_num2)
        if max_num > self.atom_count:
            raise ValueError("atom index {} exceeds the number of atoms in molecule."
                             .format(max_num))

        min_num = min(atom_num1, atom_num2)
        if min_num < 1:
            raise ValueError("atom index {} is smaller than 1".format(min_num))

        log.debug('rotate atoms around the axis which goes through: '
                  '\natom {}: {}\natom {}: {}\nby {} degrees/radians'.format(atom_num1, self.atoms[atom_num1],
                                                                             atom_num2, self.atoms[atom_num1],
                                                                             angle))

        # correct the atom numbering to start from 0
        # since python uses C numbering, which starts at 0 and not 1
        atom_num1 -= 1
        atom_num2 -= 1

        if not radians:
            # convert to radians
            angle = angle * np.pi / 180

        coordinates = self.positions

        # get rotation vector, and normalize
        rot_vec = coordinates[atom_num1] - coordinates[atom_num2]
        rot_vec /= np.linalg.norm(rot_vec)

        # shift the coordinate system, such that (0,0,0) will be located
        # on one of the atoms (atom_num1 or atom_num2) which define
        # the rotation angle
        coord_shift = np.copy(coordinates[atom_num1])
        coordinates -= coord_shift

        # form the rotation matrix
        # see: http://goo.gl/lUCF3P for wikipedia article
        rot_mat = np.zeros((3, 3), dtype=np.float)

        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)

        rot_mat[0, 0] = cos_theta + rot_vec[0] ** 2 * (1 - cos_theta)
        rot_mat[0, 1] = rot_vec[0] * rot_vec[1] * (1 - cos_theta) - rot_vec[2] * sin_theta
        rot_mat[0, 2] = rot_vec[0] * rot_vec[2] * (1 - cos_theta) + rot_vec[1] * sin_theta

        rot_mat[1, 0] = rot_vec[1] * rot_vec[0] * (1 - cos_theta) + rot_vec[2] * sin_theta
        rot_mat[1, 1] = cos_theta + rot_vec[1] ** 2 * (1 - cos_theta)
        rot_mat[1, 2] = rot_vec[1] * rot_vec[2] * (1 - cos_theta) - rot_vec[0] * sin_theta

        rot_mat[2, 0] = rot_vec[2] * rot_vec[0] * (1 - cos_theta) - rot_vec[1] * sin_theta
        rot_mat[2, 1] = rot_vec[2] * rot_vec[1] * (1 - cos_theta) + rot_vec[0] * sin_theta
        rot_mat[2, 2] = cos_theta + rot_vec[2] ** 2 * (1 - cos_theta)

        coordinates[atom_list] = np.dot(rot_mat, coordinates[atom_list].T).T

        # shift the coordinate system back to original center
        coordinates += coord_shift

        # update the locations of the molecule object
        self.positions = coordinates

    def center_on_com(self):
        """

        Change the coordinates of the molecule, to be centered on the center-of-mass

        """

        # shift coordinates to center of mass
        self.positions -= self.center_of_mass

    def move_atoms(self, atom_num1, atom_num2, distance=0, atom_list=None, units="bohr"):
        """

        move atoms along the vector formed by connecting atom1 and atom2

        :param distance: The distance in bohr

        :param atom_list: a list of atoms to be rotated

        :param units: can be "Bohr" or "Angstrom", only the first letter is important, not case sensitive

        :return: Changes the internal state (more specifically: the location)
        of the location of the atoms in atoms_list
        """

        if not atom_list:
            raise ValueError("No atoms to be moved are stated")

        else:
            # correct user input to start from 0 (not user 1..)
            atom_list = np.array(atom_list, dtype=int)
            atom_list -= 1

        max_num = max(atom_num1, atom_num2)
        if max_num > self.atom_count:
            raise ValueError("atom index {} exceeds the number of atoms in molecule."
                             .format(max_num))

        min_num = min(atom_num1, atom_num2)
        if min_num < 1:
            raise ValueError("atom index {} is smaller than 1".format(min_num))

        log.debug('move atoms along the axis which goes through: '
                  "\natom {}: {}\natom {}: {}\nby {} Bohr/Angstroms".format(atom_num1, self.atoms[atom_num1],
                                                                            atom_num2, self.atoms[atom_num1],
                                                                            distance))

        # correct the atom numbering to start from 0
        # since python uses C numbering, which starts at 0 and not 1
        atom_num1 -= 1
        atom_num2 -= 1

        if units[0].lower() == 'a':
            # convert to Bohr
            distance = distance * const.angstrom_to_bohr

        elif units[0].lower() == 'b':
            # do nothing
            pass

        else:
            raise ValueError("Unit type \"{}\" is not identified".format(units))

        coordinates = self.positions

        # get rotation vector, and normalize
        move_vec = coordinates[atom_num1] - coordinates[atom_num2]
        move_vec /= np.linalg.norm(move_vec)
        move_vec *= distance

        coords = self.positions
        coords[atom_list] += move_vec

        # update the locations of the molecule object
        self.positions = coords

    @property
    def mass_array(self):
        """

        Return a numpy float array of the masses of the molecules atoms
        with the same ordering.

        If N is the number of atoms, the size of the array is N x 1

        """
        if not self.atoms:
            raise Exception('There are no atoms set for this molecule')

        mass_list = []

        for atom in self.atoms:
            mass_list.append(atom.mass)

        return np.array(mass_list).astype(float)

    @property
    def positions(self):
        """
        Return a numpy float array of the positions of the molecules atoms
        with the same ordering.

        If N is the number of atoms, the size of the array is N x 3
        """

        position_list = []

        for atom in self.atoms:
            position_list.append(atom.coords)

        positions = np.array(position_list).astype(float)

        return positions

    @positions.setter
    def positions(self, new_positions):
        if not new_positions.shape[0] == self.atom_count:
            raise IndexError("number new atomic positions do not match number of atoms")

        counter = 0

        for atom in self.atoms:
            atom.coords = new_positions[counter]
            counter += 1

    @property
    def mass(self):
        """
        :return: The weight of the molecule in atomic units
        """

        return self.mass_array.sum()

    @property
    def center_of_mass(self):
        """

            Get the center of mass
            return a vector of the form [x_position, y_position, z_position]

            units are Bohr (atomic)

        """

        com = np.dot(self.mass_array, self.positions) / self.mass_array.sum()

        return com

    @property
    def inertia_tensor(self):
        """Get the moments of inertia along the principal axes.

        Units of the moments of inertia are amu*bohr**2.

        :returns a 1 x 3 vector of momenta values, a 3 x 3 matrix of which
        the columns are the principle axes of inertia.

        """
        com = self.center_of_mass
        positions = self.positions
        positions -= com  # translate center of mass to origin
        masses = self.mass_array

        # initialize elements of the inertial tensor
        i_xx = 0.0
        i_yy = 0.0
        i_zz = 0.0
        i_xy = 0.0
        i_xz = 0.0
        i_yz = 0.0

        for i in range(self.atom_count):
            x, y, z = positions[i]
            m = masses[i]

            i_xx += m * (y ** 2 + z ** 2)
            i_yy += m * (x ** 2 + z ** 2)
            i_zz += m * (x ** 2 + y ** 2)
            i_xy += -m * x * y
            i_xz += -m * x * z
            i_yz += -m * y * z

        inertia_tensor = np.array([[i_xx, i_xy, i_xz],
                                   [i_xy, i_yy, i_yz],
                                   [i_xz, i_yz, i_zz]])

        moments, principle_axes = np.linalg.eigh(inertia_tensor)

        sort_perm = moments.argsort()

        # sort them
        moments.sort()
        principle_axes = principle_axes[sort_perm]

        return moments, principle_axes

    @property
    def principal_axes(self):
        """
        Get the moments of inertia along the principal axes.

        Units of the moments of inertia are amu*bohr**2.

        :returns a 3 x 3 matrix of which the columns are the principle axes of inertia.
        """

        return self.inertia_tensor[1]

    def split(self, atom_number):
        """
        splits the molecule into several molecules

        :param atom_number: a list number of atoms, which will  be the spliting
        points of the molecule. The molecule is split after the giver atom.

        for instance, if molecule M has 18 atoms
        M.split([5, 6, 12]) would return in the following:
        [M1, M2, M3, M4]

        M1 contains the atoms 0..5
        M2 contains the atom 6
        M3 contains the atoms 7..12
        M4 contains the atoms 13..17

        :return: list of molecules
        """

        if not isinstance(atom_number, list):
            atom_number = [atom_number]

        try:
            point_list = np.array(atom_number, dtype=int)

        except ValueError:
            print('there are non-integer values')
            return

        point_list = point_list[point_list > 0]
        point_list = np.unique(point_list)
        point_list = np.sort(point_list)

        if point_list[-1] >= self.atom_count:
            raise ValueError('Cannot split molecule at atom number {} - '
                             'not enough atoms'.format(point_list[-1]))

        log.debug('splitting the molecule at atoms number: {}'.format(point_list))

        molecule_counter = 0

        molecule_list = list()

        molecule_list.append(Molecule())

        for i in range(self.atom_count):

            molecule_list[molecule_counter].add_atom(self.atoms[i])

            if i in point_list - 1:  # - 1 to correct numbering to start from 1
                # and not from 0 like C lang
                molecule_list.append(Molecule())
                molecule_counter += 1

        return molecule_list

    @property
    def symmetry_number(self):
        return int(1)

    @property
    def atom_count(self):
        return len(self._atoms)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, new_atoms):
        self._atoms = new_atoms

    @property
    def is_linear(self):
        """
        Returns True value if molecule is linear
         as determined by the Inertia moments
        """

        moments, axes = self.inertia_tensor

        if moments[0] == 0.0 and (moments[1] - moments[2]) < self.mom_thresh:
            linear = True

        else:
            linear = False

        return linear

    def to_jmol(self):
        """
        use JMol to visualize molecule
        """

        # an xyz file that contains all the information about the molecule
        jmol_script = tempfile.NamedTemporaryFile(suffix='.xyz', prefix="QCkitVis_", delete=False)

        jmol_script.write(self.to_xyz_format())

        jmol_script.flush()

        jmol_script.close()

        if shutil.which('jmol'):

            sp.call(['jmol', "{}".format(jmol_script.name)])

        else:

            try:

                jmol_home = os.environ["JMOLHOME"]

                sp.call(["java", "-jar", jmol_home + "\Jmol.jar", "-i", "{}".format(jmol_script.name)])

            except KeyError:
                print("JMOLHOME environment variable is not defined. (The thing that points to the root directory of JMOL)")

            except OSError:
                print('cannot find JMOL executable')

    def to_molden_format(self):
        """
        write the molecular data in xyz format

        :return: the content of an xyz formatted file
        """
        output = ""

        output += "[Molden Format]\n"
        output += "[Atoms] (AU)\n".format(self.atom_count)

        counter = 1

        for atom in self.atoms:
            output += "{:>5} {:>5} {:>5} {:10f} {:10f} {:10f}\n" \
                .format(atom.symbol,
                        counter,
                        atom.atomic_number,
                        atom.coords[0],
                        atom.coords[1],
                        atom.coords[2])
            counter += 1

        # also, add normal modes data
        if self.normal_modes:
            output += "[FREQ]\n"
            for freq in self.normal_modes.frequency_array:
                output += "{}\n".format(freq)

            output += "[FR-COORD]\n"

            for atom in self.atoms:
                output += "{:>5} {: 10f} {: 10f} {: 10f}\n" \
                    .format(atom.symbol,
                            atom.coords[0],
                            atom.coords[1],
                            atom.coords[2])
                counter += 1

            output += "[FR-NORM-COORD]\n"
            mode_counter = 0

            for mode in self.normal_modes:
                mode_counter += 1
                output += "vibration {}\n".format(mode_counter)
                output += mode.motion_print

            log.info('Warning - intensities are not really implemented')
            output += "[INT]\n"
            for mode in self.normal_modes:
                output += "1.0\n".format(mode_counter)

        return output

    def to_xyz_format(self):
        """
        Generate an xyz format text with molecular specifications

        :return: xyz format text with molecular specifications
        """

        xyz_output = ""

        xyz_output += "{:>5}\n{}\n".format(self.atom_count, "pyQChem output")

        xyz_output += self.xyz()

        return str.encode(xyz_output)

    def xyz(self, units="angstrom"):
        """
        obtain the coordinate of the atoms

        :return: a list of symbols and coordinates (?)
        """

        xyz_text = ""
        units = units.lower()

        if units[0].lower() == "b":
            factor = 1
            units = "Bohr"

        elif units[0].lower() == "a":
            factor = const.bohr_to_angstrom
            units = "Angstrom"

        else:
            raise ValueError("units requested \"{}\" not recognized".format(units))

        log.info('writing geometry in {} units'.format(units))

        for atom in self.atoms:
            xyz_text += "{:>5} {:< 15.10f} {:< 15.10f} {:< 15.10f}\n" \
                .format(atom.symbol,
                        atom.coords[0] * factor,
                        atom.coords[1] * factor,
                        atom.coords[2] * factor)

        return xyz_text

    @property
    def normal_modes(self):
        if self._normal_modes:
            return self._normal_modes
        else:
            raise Exception('No normal-modes data was assigned for this molecule')

    def to_file(self, file_name="pyQChem_output", file_format="xyz"):
        """
        writes a file with molecular specification

        :param file_name: the name of the file to write

        :param file_format: the file format. Currently supported formats are:
                xyz and Molden

        :return: None
        """

        file_format = file_format.lower()

        supported_formats = ["xyz", "molden"]

        if file_format not in supported_formats:
            raise NotImplementedError("format {} is not supported "
                                      "(yet?)".format(file_format))

        if file_format == "xyz":
            file_text = self.to_xyz()
            file_name += ".xyz"

        elif file_format == "molden":
            file_text = self.to_molden_format()
            file_name += ".molden"

        else:
            raise NotImplementedError("format {} is not supported "
                                      "(yet?)".format(file_format))

        f = open(file_name, 'w')
        f.write(file_text)
        f.close()

        log.info('write the file {} to directory {}'.format(file_name, os.getcwd()))

    def get_fragment(self, atom_numbers):
        """
        This functions returns a Molecule object, which is composed
        from atoms of the parent molecule

        Data is copied

        :return: a molecule object
        """

        fragment = Molecule(verbosity=self.verbosity)

        for i in atom_numbers:
            fragment.add_atom(copy.copy(self.atoms[i]))

        if self.normal_modes:

            fragment._normal_modes = copy.copy(self.normal_modes)
            for nm in fragment._normal_modes:
                # remove unwanted data from original mode data
                nm.select_coordinates(atom_numbers)

        return fragment

    def __str__(self):
        output = ""

        if not self.atoms:
            raise Exception('There are no atoms in this Molecule object')

        else:
            for atom in self.atoms:
                output += "{:4} {:.2e} {: e} {: e} {: e}\n" \
                    .format(atom.symbol,
                            atom.mass,
                            atom.coords[0],
                            atom.coords[1],
                            atom.coords[2])

        return output

    def __add__(self, other_molecule):
        """
        an operator to add two molecule together

        :return: a molecule object with the atoms of the two molecules
        joined
        """

        joined_molecule = Molecule()

        if self.atoms:
            for atom in self.atoms:
                joined_molecule.add_atom(atom)

        if other_molecule.atoms:
            for atom in other_molecule.atoms:
                joined_molecule.add_atom(atom)

        return joined_molecule

    def __getitem__(self, item):
        raise Exception("__getitem__ not implemented!")


def join(molecule_list):
    new_mol = Molecule()

    for mol in molecule_list:
        new_mol += mol

    return new_mol


def from_molden(file_name, molden_job_num=1, verbosity=1):
    """
    parse a molecule object out of a molden format

    :param: file_name: the file containing the molden format

    :param molden_job_num: the number of molden output. Relevant
                        if a job contains several "molden format"
                        groups, use the output from molden_job.

    :param verbosity: how much output: 0 none, 1 some (default),
                                            2 (or higher) debug.

    :return: a molecule object
    """

    parsed_molden = molden_parser.MoldenIO(file_name,
                                           molden_job_num=molden_job_num,
                                           verbosity=verbosity)

    parsed_molden.molecule.charge = 0

    parsed_molden.molecule.multiplicity = 0

    return parsed_molden.molecule


def from_xyz(file_name):
    """
    generate a molecule object from xyz formatted file

    :return: a molecule object
    """

    new_molecule = Molecule()

    f = open(file_name)

    log.info('parsing xyz file {}'.format(file_name))

    text = f.readlines()

    read_flag = False
    atom_counter = 0
    atom_number = -1
    skip_line = False

    for line in text:
        line = line.strip().split()

        if atom_counter - atom_number == 0:
            break

        elif skip_line:
            skip_line = False
            continue

        elif not read_flag:
            if len(line) == 1:
                atom_number = int(line[0])
                read_flag = True
                skip_line = True
                log.info('xyz file contains {} atoms'.format(atom_number))

        elif read_flag:
            atom_counter += 1
            new_atom = Atom(sym=line[0], coords=[line[1], line[2], line[3]],
                            coord_units="angs")
            new_molecule.add_atom(new_atom)

    new_molecule.charge = 0

    new_molecule.multiplicity = 1

    return new_molecule


if __name__ == "__main__":
    pass
