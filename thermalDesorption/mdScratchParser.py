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

    @property
    def trajectory(self, first_step=False):
        """
        collect data from aimd scratch and create trj file

        """

        # create trjectory file
        with open(self.aimd_scratch + "/View.xyz", 'r') as f:
            trj_data = str(f.read())

        return trj_data

    @property
    def get_temp_and_potential(self):
        """
        collect data from aimd scratch and create trj file

        """

        results = {'time': list(),
                   'potentialE': list(),
                   'temperature': list()}

        # create trajectory file
        with open(self.aimd_scratch + "/TandV", 'r') as f:
            tandv_data = str(f.read())

        # remove header line
        tandv_data = tandv_data.splitlines()
        tandv_data.pop(0)

        for line in tandv_data:

            line = line.split()

            results['time'].append(float(line[0]))
            results['potentialE'].append(float(line[1]))
            results['temperature'].append(float(line[2]))

        results['time'] = np.array(results['time'], dtype=np.float)
        results['potentialE'] = np.array(results['potentialE'], dtype=np.float)
        results['temperature'] = np.array(results['temperature'], dtype=np.float)

        return results

    @property
    def finished_ok(self):

        for line in self.content:

            if "Have a nice day" in line:

                return True

            else:
                return False
