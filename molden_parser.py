import numpy as np
import os
import logging as log

import QCkit.molecule as molecule
from QCkit.atom import Atom


# Normal modes stuff is commented out


class MoldenIO:
    """

    a parser for molden format output

    """

    def __init__(self, file_name, molden_job_num=1):
        """

        :param file_name: the name of the file to parse

        :param molden_job_num: the number of molden output. Relevant
                            if a job contains several "molden format"
                            groups, use the output from molden_job.

        :returns an initialized Molden object
        """

        self._energy = None

        # talk to the user
        log.info('\n{:*^40}'.format(' Molden Format parsing '))
        log.info('Reading file \"{}\"'.format(file_name))

        # can the file be found?
        if not os.path.isfile(file_name):
            raise IOError('\"{}\" input file not found'.format(file_name))

        # open output file and read content
        f = open(file_name, 'r')
        text = f.read()
        f.close()

        # is there a "molden format" section in it?
        molden_count = text.count("[Molden Format]")

        # talk to the user
        log.info('Detected {} [Molden Format] groups. '
                 'User asked number {}'.format(molden_count, molden_job_num))

        # if the output is complex, there can be several [Molden Format]
        # groups in the file. Or there can be none at all.
        if molden_count == 0:
            print('[Molden Format] section not found in input file')
            return

        while molden_count < molden_job_num and molden_job_num >= 0:
            print("There are only {0} molden formats, please enter a number 1 .. {0}:".format(molden_count))
            try:
                molden_job_num = int(input())

            except ValueError:
                print('This is not an integer number')

        # parse the relevant [Molden Format] section
        molden_format = text.split("[Molden Format]")[molden_job_num]

        # parse molecular data
        # molecular data is stored in a special "Molecule" object
        if "[Atoms]" in molden_format:
            log.info('Detected Atomic data. '
                     'Parsing geometry')
            self._atoms = get_atoms(molden_format)
            log.debug(self._atoms)

        # parse frequencies. frequencies are simply a numpy array of floats
        if "[FREQ]" in molden_format:
            log.info('Detected frequency data. '
                     'Parsing frequencies (cm -1 units)')
            self._frequencies = get_frequencies(molden_format)
            log.debug(self._frequencies)

            # parse normal mode data.
            # normal modes are stored in a special "Motions" object
            # if "[FR-NORM-COORD]" in molden_format:
            #     log.info('Detected normal modes. '
            #                             'Parsing molecular motions')
            #     self._normal_modes = get_normal_modes(molden_format, self._frequencies)
            #     log.debug(self._normal_modes)

        if "[GEOCONV]" in molden_format:
            log.info('Detected optimization data. ')
            self._energy = get_last_energy(molden_format)
            log.debug('energy is {} Hartree'.format(self._energy))

    @property
    def molecule(self):
        return self._atoms

    @property
    def atoms(self):
        return self._atoms

    @property
    def energy(self):
        """
        return the energy of the molecule.
        """

        return self._energy

        # @property
        # def normal_modes(self):
        #     return self._normal_modes


def get_atoms(molden_format):
    """
    extraction of the [Atoms] section from a molden output

    :return a molecule object
    """

    # first, locate the [Atoms] sections
    atoms_start = molden_format.find('[Atoms]')

    # new molecule object
    molec = molecule.Molecule()

    if atoms_start == -1:
        raise Exception('[Atoms] section no found in [Molden Format]')

    # start scanning from the start of [Atoms]
    lines = molden_format[atoms_start + 1:].splitlines()

    # remove redundant first line (contains [Atom] (angs))
    units = lines.pop(0)

    # what are the units?
    if "(Angs)" in units:
        units = "ang"
    else:
        units = "bohr"
        # raise Exception('Geometry is stated in Bohr units is not currently supported\n'
        #                 'but is very easy to add!')

    for line in lines:

        # make sure that this is an atomic specifications,
        # the only verification is the number of terms in
        # the line
        #
        # and if not, than continue

        line = line.split()
        if len(line) != 6:
            break
        elif not line[0].isalnum() and \
                line[2].isdigit() and \
                isfloat(line[3]) and \
                isfloat(line[4]) and \
                isfloat(line[5]):
            break
        else:
            # do not use the atomic mass specification
            # as given by the molden output. It is not
            # accurate enough.
            molec.add_atom(Atom(line[0], line[3:6], coord_units=units))

    # if units == "bohr":
    #     mol.center_of_mass

    return molec


# def get_normal_modes(molden_format, frequency_data):
#     """
#     Procedure for extraction the vibrational data from
#     molden format output
#
#     in the HMA analysis, these are not vibrations,
#     but motions, hence their name.
#
#     :param molden_format:  the output file
#     :return:
#     """
#
#
#     # first, locate the [Atoms] sections
#     motions_start = molden_format.find('[FR-NORM-COORD]')
#
#     if motions_start == -1:
#         raise Exception('[FR-NORM-COORD] section no found in [Molden Format]')
#
#     # start scanning from the start of [FR-NORM-COORD]
#     lines = molden_format[motions_start+1:].splitlines()
#
#     # remove redundant first line, which contains the [FR-NORM-COORD] header
#     lines.pop(0)
#
#     # remove another line, that contains "vibration 1"
#     # this will be useful, because the word "vibration
#     # will be used to identify  a new parsed motion
#     lines.pop(0)
#
#     coords_list = []
#
#     modes_counter = -1
#
#     for line in lines:
#
#         # make sure that this is an atomic specifications
#         line = line.split()
#
#         if len(line) == 3 and isfloat(line[0]) and isfloat(line[1]) and isfloat(line[2]):
#             coords_list.append(line)
#
#         elif "vibration" in line:
#             modes_counter += 1
#             m = NMode(coords_list, frequency_data[modes_counter])
#             coords_list = []
#             molecular_motions.add_mode(m)
#
#         else:
#             modes_counter += 1
#             m = NMode(coords_list, frequency_data[modes_counter])
#             molecular_motions.add_mode(m)
#             break
#
#     if modes_counter + 1 != len(frequency_data):
#         raise Exception("[Freq] data is not the same size as [FR-NORM-COORD] data")
#
#     return molecular_motions


def get_frequencies(molden_format):
    """
    a procedure for extraction the vibrational frequencies

    :param molden_format:  the output file
    :return: a numpy array of frequencies (cm -1 units)
    """

    freq_list = []

    # first, locate the [FREQ] sections
    freq_start = molden_format.find('[FREQ]')

    if freq_start == -1:
        raise Exception('[FREQ] section no found in [Molden Format]')

    # start scanning from the start of [FREQ]
    lines = molden_format[freq_start + 1:].splitlines()

    # remove redundant first line, which contains the [FREQ] header
    lines.pop(0)

    for line in lines:

        # make sure that this is an FREQ specifications,
        # the only verification is the number of terms in
        # the line
        #
        # and if not, continue

        line = line.split()
        if len(line) == 1 and isfloat(line[0]):
            freq_list.append(line[0])
        else:
            break

    freq_array = np.array(freq_list).astype(np.float)

    return freq_array


def get_last_energy(molden_format):
    """
    a procedure for extraction the the last energy,
    in case of geometry optimization

    :param molden_format:  the output file
    :return: the energy (Hartree units)
    """

    energy = None

    # first, locate the "energy" sections
    energy_start = molden_format.find('energy')

    if energy_start == -1:
        raise Exception('\"energy\" section no found in [Molden Format]')

    # start scanning from the start of "energy"
    lines = molden_format[energy_start + 1:].splitlines()

    # remove redundant first line, which contains the "energy" header
    lines.pop(0)

    for line in lines:

        line = line.strip()
        if line == "max-force":
            break

        else:
            energy = float(line)

    return energy


def isfloat(text):
    """
    check whether the text string represent a float number

    :param text: a string
    :return: True if the text represents a float number, False otherwise
    """

    # in case this is a negative number
    if text[0] == "-":
        text = text[1:]

    if text.isdigit():
        return True

    elif "." in text:
        text = text.split(".")
        if len(text) > 2:
            return False
        else:
            if text[0].isdigit() and text[1].isdigit():
                return True


if __name__ == "__main__":
    mol = MoldenIO("./example_outputs/acetylene_opt_freq.qchem")
    print(mol.energy)
