from QCkit.molecule import Molecule
from QCkit.atom import Atom
from QCkit.thermalDesorption import simulatedTDP


class TPDrestart:
    def __init__(self,
                 pos_vel_file,
                 basis,  # basis set to use
                 exchange,  # exchange to use (such as B3LYP etc.)
                 high_temperature,  # maximum temperature of the simulation
                 temp_advance,  # change in temperature each simulation step
                 aimd_step_duration,  # time step of the MD simulation, in atomic units = 0.0242 fs
                 aimd_num_steps,  # number of steps for each AIMD run (one run of each temperature)
                 thermostat_timescale,  # friction of system with the thermostat head-bath
                 thermostat="langevin",  # type of thermostat
                 threads=None):  # number of openmp threads to use. If none then is OMP_NUM_THREADS
        """

        :return:
        """

        # read content of restart file
        with open(pos_vel_file, 'r') as f:
            content = f.read()

        # extract job name
        job_name = pos_vel_file.replace('.last_pos_vel', '')

        # get the last position and velocity
        content = content.split('Ended Q-Chem aimd run with the following positions and velocities:\n')
        last_temperature = int(content[-2].splitlines()[-1].split()[2])

        if last_temperature + temp_advance >= high_temperature:
            raise Exception('Old job started with temperature already higher then new target')

        content = content[-1].splitlines()

        molecule = Molecule()
        velocities = ""

        for line in content:
            line = line.split()
            molecule.add_atom(Atom(sym=str(line[0]),
                                   coords=(line[1:4]),
                                   coord_units="ang"))

            velocities += "{} {} {}\n".format(line[4], line[5], line[6])

        tpd_simulation = simulatedTDP.TPD(molecule=molecule,
                                          tpd_job_name=job_name,
                                          high_temperature=high_temperature,
                                          low_temperature=last_temperature + temp_advance,
                                          temp_advance=temp_advance,
                                          aimd_step_duration=aimd_step_duration,
                                          aimd_num_steps=aimd_num_steps,
                                          basis=basis,
                                          threads=threads,
                                          exchange=exchange,
                                          thermostat=thermostat,
                                          thermostat_timescale=thermostat_timescale,
                                          velocities=velocities,
                                          restart=True)


if __name__ == "__main__":
    TPDrestart(pos_vel_file='./example.last_pos_vel',
               basis="6-31g*",
               exchange="b97-d3",
               threads=str(4),
               high_temperature=124,
               temp_advance=2,  # advance temperature by 2 K
               thermostat="langevin",
               thermostat_timescale=120,  # fs
               aimd_step_duration=60,  # atomic units
               aimd_num_steps=10)
