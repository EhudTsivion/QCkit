import numpy as np


class TrjProcessor:
    """
    This class processes the TRJ files
    which contain the results from the TPD simulation

    The main products of this process are:

    1. Graphs of the distance between the metal and the hydrogen atoms
    2. Detection of desorption events

    How it works:

    (a) Parse the trj file into a list of individual data sets
    (b) The H2 molecule and metal ion are detected by looking
        at the first XYZ data set.
    (c) The relative position with respect to the metal is calculated
    (d) Detection of desorption event.
    (e) Optionally: print report into PDF file

    """

    def __init__(self, trj_file):

        # these are adjusted in the 'parse_trj' method
        self.atom_count = None
        self.data_length = None
        self.h2_metal_distance = None
        self.temperature_vec = None
        self.time_vec = None


        self.data_set = self.parse_trj(trj_file)
        self.metal_index = self.find_metal_index()
        self.hydrogens = self.find_hydrogen_indices()
        self.temperature_vec = self.get_temperature_vec()
        self.time_vec = self.get_time_vec()


        self.get_distance_from_metal()


    def parse_trj(self, trj_file):
        """

        :param      trj_file: A name of a file to be parsed
        :return:    a list the contains atomic positions,
                    target temperature, instantaneous temperature and
                    step number.
        """

        data_set = []

        with open(trj_file, 'r') as f:
            content = f.readlines()

        self.atom_count = int(content[0])

        line_number = len(content)

        # get the total number of xyz data points
        # first line is atom number second is comment
        # third id atomic coordinates
        xyz_num = int(line_number / (self.atom_count + 2))

        data = dict()

        # iterate over XYZ data sets
        for i in range(xyz_num):

            # get rid of first line
            content.pop(0)

            # next line contains data about
            # simulation time and temperature
            line = content.pop(0).split()

            target_temperature = line[1]
            simulation_time = line[-2]

            # next lines contain xyz data
            coords = list()
            atoms = list()

            for j in range(self.atom_count):
                line = content.pop(0).split()
                coords.append(np.array(line[1:4], dtype=float))
                atoms.append(str(line[0]))

            data = {'temperature': target_temperature,
                    'time': simulation_time,
                    'atoms': atoms,
                    'coordinates': coords}

            data_set.append(data)

            self.data_length = len(data_set)

        return data_set

    def find_metal_index(self):
        """

        :return: the index (from 0 .. Natoms-1) of the
          metal ion, Ca or Mg are currently implemented
        """

        index = None

        for i in range(self.atom_count):

            # inspect only the first data point
            atom = self.data_set[0]['atoms']

            try:
                index = atom.index('Mg')

            except ValueError:
                pass

            try:
                index = atom.index('Ca')

            except ValueError:
                pass

        if index:
            return index

        else:
            raise Exception("Cannot fund metal (Ca or Mg - UppercaseLowecase) ion in TRJ file.\n"
                            "Something must be very bad.")

    def find_hydrogen_indices(self):
        """

        :return:
        """

        indices = list()

        # inspect only the first data point
        atoms = self.data_set[0]['atoms']

        for i in range(self.atom_count):

            atom = atoms[i]

            if atom.lower() == "h":

                for j in range(i + 1, self.atom_count):

                    atom2 = atoms[j]

                    if atom2.lower() == 'h':

                        pos_1 = self.data_set[0]['coordinates'][i]
                        pos_2 = self.data_set[0]['coordinates'][j]

                        distance = np.linalg.norm(pos_1 - pos_2)

                        if 0.5 < distance < 1:

                            indices.append([i, j])

        if not indices:
            raise Exception("Can't find hydrogen molecules. Something must be very wrong")

        return indices

    def get_temperature_vec(self):

        temperature_vec = np.zeros(self.data_length, dtype=np.float)

        for i in range(self.data_length):

            temperature_vec[i] = self.data_set[i]['temperature']

        return temperature_vec

    def get_time_vec(self):

        time_vec = np.zeros(self.data_length, dtype=np.float)

        for i in range(self.data_length):
            time_vec[i] = self.data_set[i]['time']

        return time_vec

    def get_distance_from_metal(self):
        """
        calculate the distance from the center of the hydrogen molecule
        to the metal ion.

        to do this with maximum efficiency, we an array for positions
        of the center of mass of each hydrogen molecule and for the metal
        ion.

        :return:    numpy array size N-hydrogens X number_of_simulation steps
                    containing the distance.

        """

        # an array to hold the center of mass of the hydrogen
        # molecules

        h2_com_xyz = np.zeros((self.data_length, self.hydrogen_number, 3), dtype=np.float)
        metal_xyz = np.zeros((self.data_length, 3), dtype=np.float)
        distance_vec = np.zeros((self.data_length, self.hydrogen_number), dtype=np.float)

        for i in range(self.data_length):

            metal_xyz[i] = self.data_set[i]['coordinates'][self.metal_index]

            for j in range(self.hydrogen_number):

                a = self.data_set[i]['coordinates'][self.hydrogens[j][0]]
                b = self.data_set[i]['coordinates'][self.hydrogens[j][1]]

                h2_com_xyz[i][j] = np.true_divide(a + b, 2)

        for i in range(self.hydrogen_number):

            distance_vec[:, i] = np.linalg.norm(h2_com_xyz[:, i] - metal_xyz[:], axis=1)

        self.h2_metal_distance = distance_vec

        return distance_vec

    @property
    def hydrogen_number(self):
        return len(self.hydrogens)

    def plot_h2_metal_distance(self):
        import matplotlib.pyplot as plt

        """
        plot
        """

        x = np.arange(self.data_length)

        try:
            plt.plot(self.time_vec, self.h2_metal_distance)
            plt.show()

        except ValueError:
            print("The H2 to Metal distance wasn't calculated")

        return None

if __name__ == "__main__":
    trj_data = TrjProcessor('../examples/example_trj.trj')
    trj_data.plot_h2_metal_distance()
