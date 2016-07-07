import numpy as np
import json
import os

from scipy.stats import gamma, norm
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from QCkit.thermalDesorption import tdpCurve


class TrjProcessor:
    """
    This class processes the TRJ files
    which contain the results from the TPD simulation

    The main products of this process is a json
    file which contained all the important information from
    the trajectory

    How it works:

    (a) Parse the trj file into a list of individual data sets
    (b) The H2 molecule and metal ion are detected by looking
        at the first XYZ data set.
    (c) The relative position with respect to the metal is calculated
    (d) the data is dumped into a json file

    once a json data file exists, it can be processes by the
    JsonData object.

    The reason to make this split, is that processing of the TRJ
    files take a lot of time, while processing of the JSON files
    is quick
    """

    def __init__(self, trj_file):

        # these are adjusted in the 'parse_trj' method
        self.atom_count = None
        self.data_length = None
        self.h2_metal_distance = None
        self.temperature_vec = None
        self.time_vec = None
        self.metal_xyz = None  # position of the metal ion
        self.h2_com_xyz = None  # positions of the H2 center-of-mass

        self.file_name = trj_file
        self.data_set = self.parse_trj(self.file_name)
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
        """
        collected the time data into a vector

        in case of a restart run
        it will correct the time
        assuming time-step similar to
        original run

        :return:
        """

        time_vec = np.zeros(self.data_length, dtype=np.float)

        for i in range(0, self.data_length):
            time_vec[i] = self.data_set[i]['time']

        # validates that time actually advances
        # this is important for parsing restart jobs
        # were time goes back

        flag = False

        for i in range(2, self.data_length):

            diff = time_vec[i - 1] - time_vec[i - 2]

            if time_vec[i] - time_vec[i - 1] < 0:  # this is an indicator of new time-keeping due to restart
                time_vec[i] = '{:0.2f}'.format(time_vec[i - 1] + diff)
                flag = True

        if flag:
            print('found possible restart job')

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

        self.h2_com_xyz = np.zeros((self.data_length, self.hydrogen_number, 3), dtype=np.float)
        self.metal_xyz = np.zeros((self.data_length, 3), dtype=np.float)
        distance_vec = np.zeros((self.data_length, self.hydrogen_number), dtype=np.float)

        for i in range(self.data_length):

            self.metal_xyz[i] = self.data_set[i]['coordinates'][self.metal_index]

            for j in range(self.hydrogen_number):
                a = self.data_set[i]['coordinates'][self.hydrogens[j][0]]
                b = self.data_set[i]['coordinates'][self.hydrogens[j][1]]

                self.h2_com_xyz[i][j] = np.true_divide(a + b, 2)

        for i in range(self.hydrogen_number):
            distance_vec[:, i] = np.linalg.norm(self.h2_com_xyz[:, i] - self.metal_xyz[:], axis=1)

        self.h2_metal_distance = distance_vec

        return distance_vec

    @property
    def hydrogen_number(self):
        return len(self.hydrogens)

    def dump_data(self):
        """
        writes the analyzed trajectory data as json file

        The data regarding the distance between the H2
        molecules and the metal ion is dumped. Its shape is:
        (N, M) where N is the number of H2 molecules and M is the
        number of steps in the simulation

        the information about the time of the simulation as

        """

        dump_file_name = self.file_name.replace('.trj', '.json')

        dump_data = dict()

        dump_data['simulation_temperature'] = self.temperature_vec.tolist()
        dump_data['simulation_time'] = self.time_vec.tolist()
        dump_data['h2_metal_distance'] = self.h2_metal_distance.tolist()

        with open(dump_file_name, 'w') as f:
            json.dump(dump_data, f)

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


class JsonData:
    """
    Reads data from a single json file, containing trajectory data,
    and processes the information.

    """

    def __init__(self, json_file):

        with open(json_file, 'r') as f:
            content = f.read()

        self.data = json.loads(content)

    def get_desorption_temperature(self, thresh=5):
        """
        :param thresh: the distance, in A units, between the hydrogen molecule and the metal ion
        which is considered as a threshold for detachment.

        :return: a vector of detachment times the trajectory
        """

        counter = 0

        # find time of detachment

        found = False

        time = None

        for distance in self.data['h2_metal_distance']:

            # check distance for each individual h2 molecule
            for h2_molecule in distance:

                if h2_molecule >= thresh:
                    time = self.data['simulation_temperature'][counter]
                    found = True

            counter += 1

            if found:
                break

        return time

    def get_max_temp(self):
        """

        :return: return the maximal temperature that the simulation has reached
        """

        return self.data['simulation_temperature'][-1]


class DirectoryInformation:
    """
    A class for working with all the information in a given directory

    A directory, or folder, contains many json files.
    The goal of this class is to work will all the information in a
    particular directory.

    """

    def __init__(self, directory='.'):
        """
        Initialize by getting all the information in the directory

        """

        self.info = list()

        print('reading data from all JSON files in directory')

        current_dir = os.getcwd()

        os.chdir(directory)

        for f in os.listdir(directory):

            if f.endswith('.json'):
                self.info.append(JsonData(f))

        if not self.info:
            raise IOError('Did not find any valid json pared information')

        os.chdir(current_dir)

    def get_desorption_distribution(self, thresh=5, plot=False, bins=15):
        """
        Obtain the temperatures at which the molecules have desorbed
        from the metal center, by moving farther than a certain threshold
        distance.

        You can optionally plot a histogram of the distribution

        :param thresh:
        :param plot:
        :param bins:
        :return:
        """

        detachment_temperatures = list()

        for trj_info in self.info:

            time = trj_info.get_desorption_temperature(thresh=thresh)

            if time:
                detachment_temperatures.append(trj_info.get_desorption_temperature(thresh=thresh))

        detachment_temperatures = np.array(detachment_temperatures, dtype=np.float)

        if plot:
            plt.hist(detachment_temperatures,
                     range=(0, detachment_temperatures.max()),
                     bins=bins)

            plt.show()

        return detachment_temperatures

    def analyze_traj_max_time(self, hist_bins=15, plot=True):
        """
        Obtain the temperatures at which the molecules have desorbed
        from the metal center, by moving farther than a certain threshold
        distance.

        You can optionally plot a histogram of the distribution

        :param thresh:
        :param plot:
        :param bins:
        :return:
        """

        # the maximal temperature achieved by all trjectories
        max_temperatures = list()

        for trj_info in self.info:

            max_temp = trj_info.get_max_temp()

            max_temperatures.append(max_temp)

        max_temperatures = np.array(max_temperatures, dtype=np.float)

        if plot:
            plt.hist(max_temperatures,
                     range=(0, max_temperatures.max()),
                     bins=hist_bins)

            plt.show()

        return max_temperatures

    def fit_gamma_distribution(self, desorption_thresh=5, plot=False, bins=15, normed=True):

        #  first get the detachment time
        dt = self.get_desorption_distribution(thresh=desorption_thresh)

        # you need to invert the data - to be able to fit gamma
        # dt = dt.max() - dt

        fit_alpha, fit_loc, fit_beta = gamma.fit(dt)

        x = np.linspace(0, dt.max(), 100)

        pdf_fitted = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)

        # normalize
        # pdf_fitted = pdf_fitted / pdf_fitted.max()

        # this is the maximum of the distribution
        mode = x[pdf_fitted.argmax()]

        if plot:
            x = np.linspace(0, dt.max(), 100)

            plt.hist(dt, bins=bins, normed=normed)
            plt.plot(x, pdf_fitted)
            plt.show()

        return fit_alpha, fit_loc, fit_beta

    def test_all_gamma(self, thresh_min=3.5, thresh_max=6.5):

        for thresh in np.arange(thresh_min, thresh_max, 0.1):
            print(thresh, self.fit_gamma_distribution(thresh=thresh))

        return None

    def fit_normal_distribution(self, desorption_thresh=5, plot=False, bins=15, normed=True):
        """
        Fit a normal distribution to the results

        You can optionally plot the fit of the distribution

        :parameter desorption_thresh: the distance between the metal and the
        hydrogen molecule from which the molecule is considered as desorbed
        :parameter plot: a flag to turn on or off the plotting of the results
        :parameter bins: the number of bins used to plot the histogram
        :parameter normed: whether to normalize the results
        :return:
        """

        #  first get the detachment time
        dt = self.get_desorption_distribution(thresh=desorption_thresh)

        fit_loc, fit_scale = norm.fit(dt)

        if plot:
            x = np.linspace(0, dt.max(), 100)

            pdf_fitted = norm.pdf(x, fit_loc, fit_scale)

            # normalize
            pdf_fitted = pdf_fitted  # / pdf_fitted.max()

            plt.hist(dt, bins=bins, normed=normed)
            plt.plot(x, pdf_fitted, color='r')
            plt.show()

        return fit_loc, fit_scale

    def scan_desorption_threshold(self, thresh_min=3.5, thresh_max=6.5):

        for thresh in np.arange(thresh_min, thresh_max, 0.1):
            print(thresh, self.fit_normal_distribution(desorption_thresh=thresh))

        return None

    def extract_enthalpy_entropy(self,
                                 guess_enthalpy,
                                 guess_entropy,
                                 heating_rate,
                                 desorption_thresh=5.5):
        """

        :param heating_rate: The heating rate in K fs-1
        :param guess_enthalpy:
        :param desorption_thresh:
        :param guess_entropy:
        :return:
        """

        loc, scale = self.fit_normal_distribution(desorption_thresh=desorption_thresh)

        center = loc
        fwhm = scale * 2.355

        # Nelder-Mead is the only only which seems to work
        optimized_values = minimize(error2, np.array([guess_enthalpy, guess_entropy]),
                                    args=(center, fwhm, heating_rate),
                                    method='Nelder-Mead')

        results = {'enthalpy': optimized_values['x'][0],
                   'entropy': optimized_values['x'][1],
                   'error': optimized_values['fun']}

        return results


def extract_enthalpy_entropy_two_info(guess_enthalpy,
                                      guess_entropy,
                                      data_dir_1,
                                      data_dir_2,
                                      heating_rate1,
                                      heating_rate2,
                                      desorption_thresh=5.5):
    """

    :param heating_rate: The heating rate in K fs-1
    :return:
    """

    center1_target, width1_target = DirectoryInformation(data_dir_1).fit_normal_distribution(
        desorption_thresh=desorption_thresh)

    center2_target, width2_target = DirectoryInformation(data_dir_2).fit_normal_distribution(
        desorption_thresh=desorption_thresh)

    # Nelder-Mead is the only only which seems to work (?)
    optimized_values = minimize(error_twoRates, np.array([guess_enthalpy, guess_entropy]),
                                args=(center1_target,
                                      center2_target,
                                      width1_target * 2.355,  # convert width to FWHM
                                      width2_target * 2.355,  # convert width to FWHM
                                      heating_rate1,
                                      heating_rate2),
                                method='Nelder-Mead',
                                tol=0.0005,
                                options={'maxiter': 200})

    results = {'enthalpy': optimized_values['x'][0],
               'entropy': optimized_values['x'][1],
               'error': optimized_values['fun'],
               'success': optimized_values['success']}

    if not results['success']:
        raise Exception('Fitting procedure failed')

    opt_center1 = tdpCurve.TpdCurve(enthalpy_atRoomT=results['enthalpy'],
                                    entropy=results['entropy'],
                                    heating_rate=heating_rate1).max

    opt_center2 = tdpCurve.TpdCurve(enthalpy_atRoomT=results['enthalpy'],
                                    entropy=results['entropy'],
                                    heating_rate=heating_rate2).max

    print('Center of slow heating curve is: {:0.0f}, optimized {:0.0f}'.format(center1_target, opt_center1))
    print('Center of fast heating curve is: {:0.0f}, optimized {:0.0f}'.format(center2_target, opt_center2))

    return results


def error_twoRates(parameters,
                   center1_target,
                   center2_target,
                   width1_target,
                   width2_target,
                   heating_rate1,
                   heating_rate2):
    enthalpy = parameters[0]

    entropy = parameters[1]

    tdp1 = tdpCurve.TpdCurve(enthalpy_atRoomT=enthalpy,
                             entropy=entropy,
                             heating_rate=heating_rate1)

    tdp2 = tdpCurve.TpdCurve(enthalpy_atRoomT=enthalpy,
                             entropy=entropy,
                             heating_rate=heating_rate2)

    print(tdp1.fwhm)

    error = np.abs(center1_target - tdp1.max) + np.abs(center2_target - tdp2.max)

    return error


def error2(parameters, target_max, target_fwhm, heating_rate):
    enthalpy = parameters[0]

    entropy = parameters[1]

    tdp = tdpCurve.TpdCurve(enthalpy_atRoomT=enthalpy, entropy=entropy, heating_rate=heating_rate)

    error = np.abs(target_max - tdp.max) + np.abs(target_fwhm - tdp.fwhm)

    return error


def process_all_trj():
    """
    parses all TRJ files in a directory and dump
    the information into JSON files

    """

    for f in os.listdir('.'):

        if f.endswith('.trj'):
            print('Parsing and dumping data from: {}'.format(f))
            trjp = TrjProcessor(f)
            trjp.dump_data()


if __name__ == "__main__":
    thresh = 10

    di = DirectoryInformation('C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2')
    print(di.fit_normal_distribution(plot=True, desorption_thresh=thresh))

    di = DirectoryInformation(
        'C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2_cont')
    print(di.fit_normal_distribution(plot=True, desorption_thresh=thresh))
