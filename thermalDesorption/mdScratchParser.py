import path
import os
import numpy as np


class MdScratchParser:
    """
    A class for retrieving information from MD runs

    """

    def __init__(self, job_name):

        self.job_name = job_name

        try:
            scratch_dir = os.environ["QCSCRATCH"]

        except KeyError:

            print("Cannot find $QCSCRATCH directory")

        job_dir = scratch_dir + "/" + self.job_name

        if not os.path.isdir(job_dir):
            print(job_dir)

            raise NotADirectoryError("Cannot locate specific job scratch dir for job {}".format(job_name))

        self.aimd_scratch = job_dir + "/AIMD"

        if not os.listdir(self.aimd_scratch):
            raise FileNotFoundError("AMID run scratch files could not be found. Make sure -save options was used")

    def get_velocities(self):
        """
        get information about velocities from the scratch info
        """

        with open(self.aimd_scratch + "/NucVeloc", 'r') as f:
            content = str(f.read()).splitlines()

        # get last velocity
        # remove time by poping
        line = content[-1]
        line = line.split()
        line.pop(0)

        velocity = ""

        while line:
            velocity += "{:< e} {:< e} {:< e}\n".format(float(line.pop(0)),
                                                        float(line.pop(0)),
                                                        float(line.pop(0)))

        return velocity

    def get_positions(self):
        """
        get last nuclear positions from amid scratch
        units are Angstrom

        :return: an N x 3 numpy array with positions, which is required for the
        position setter defined in the Molecule object
        """

        with open(self.aimd_scratch + "/NucCarts", 'r') as f:
            content = str(f.read()).splitlines()

        # get last position
        # remove time by poping
        line = content[-1]
        line = line.split()
        line.pop(0)

        positions = []

        while line:
            positions.append([np.float64(line.pop(0)), np.float64(line.pop(0)), np.float64(line.pop(0))])

        return np.array(positions, dtype=np.float64)

    def append_trj_data(self, first_step=False):
        """
        collect data from aimd scratch and create trj file

        """

        # create trjectory file
        with open(self.aimd_scratch + "/View.xyz", 'r') as f:

            trj_data = str(f.read())

        with open(self.job_name + ".trj", 'a') as f:

            f.write(trj_data)

        # create T and V file
        with open(self.aimd_scratch + "/TandV", 'r') as f:
            tandv_data = f.read().splitlines()

        # in case of not first write, remove first line
        if not first_step:
            tandv_data.pop(0)

        tandv_data = ' '.join(tandv_data) + '\n'

        with open(self.job_name + ".TandV", 'a') as f:
            f.write(str(tandv_data))




