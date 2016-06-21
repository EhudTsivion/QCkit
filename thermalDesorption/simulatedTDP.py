import datetime
import random
import logging as log

from QCkit.thermalDesorption.mdjob import MDjob
from QCkit import physical_constants
from QCkit.thermalDesorption.mdScratchParser import MdScratchParser


class TPD:
    """
    A class for simulation of temperature programmed desorption experiments
    using molecular dynamics
    """

    def __init__(self,
                 molecule,          # the molecule object, contains all the information about the molecule
                 basis,             # basis set to use
                 exchange,          # exchange to use (such as B3LYP etc.)
                 high_temperature,  # maximum temperature of the simulation
                 low_temperature,   # starting temperature of the simulation
                 temp_advance,      # change in temperature each simulation step
                 time_step,         # time step of the MD simulation, in atomic units = 0.0242 fs
                 aimd_steps,        # number of steps for each AIMD run (one run of each temperature)
                 thermostat_timescale,      # friction of system with the thermostat head-bath
                 thermostat="langevin",     # type of thermostat
                 threads=None,          # number of openmp threads to use. If none then is OMP_NUM_THREADS
                 tpd_job_name=None):    # name of the job. appears in all related output files.

        self.low_temperature = low_temperature
        self.high_temperature = high_temperature
        self.temp_advance = temp_advance
        self.molecule = molecule
        self.time_step = time_step
        self.aimd_steps = aimd_steps
        self.basis = basis
        self.threads = threads
        self.exchange = exchange
        self.thermostat = thermostat
        self.thermostat_timescale = thermostat_timescale
        self.tpd_job_name = tpd_job_name

        self.current_temp = self.low_temperature

        if not tpd_job_name:

            self.tpd_job_name = "tpd-QChem" + datetime.datetime.now().strftime('-%Y-%m-%d-') \
                                + str(random.randint(10000, 99999))

        else:
            # always append a random number, because several
            # jobs with same name can run in parallel
            self.tpd_job_name = '{}-{}'.format(tpd_job_name, str(random.randint(10000, 99999)))

        log.basicConfig(filename="{}.log".format(self.tpd_job_name),
                        filemode='w',
                        level='INFO',
                        format='')

        self.run()

        log.info("******** Simulation ended")

    def run(self):

        rems = {"aimd_thermostat": self.thermostat,
                "aimd_time_step": self.time_step,
                "aimd_temp": self.low_temperature,
                "aimd_steps": self.aimd_steps,
                "aimd_init_veloc": "thermal",
                "aimd_print": "1",
                "max_scf_cycles": "200",
                "aimd_langevin_timescale": self.thermostat_timescale}

        # I think this MDjob thing is
        # probably unnecessary
        job = MDjob(molecule=self.molecule,
                    job_type="aimd",
                    basis=self.basis,
                    threads=self.threads,
                    exchange=self.exchange,
                    molden_format="False",
                    job_name=self.tpd_job_name,
                    rems=rems)

        log.info("\n\n{:*^30}".format("new TPD simulation"))
        log.info("\ncomputational details:")
        log.info("basis: {}, exchange: {}".format(self.basis, self.exchange))
        log.info("starting temperature {}, ending temperature {} K".format(self.low_temperature,
                                                                           self.high_temperature))
        log.info("advancing with steps of {} K".format(self.temp_advance))

        log.info("using {} thermostat".format(job.rems["aimd_thermostat"]))

        first_run = True

        simulation_time = 0  # keep record of the time

        for temp in range(self.low_temperature,
                          self.high_temperature,
                          self.temp_advance):

            job.change_temperature(temp)

            self.current_temp = temp

            log.info("\n++ starting a new Q-Chem AIMD job at {} K".format(self.current_temp))

            job.run()

            if not job.failed:
                log.info('Q-Chem finished {} K run successfully'.format(self.current_temp))

            else:
                log.info("Q-Chem job failed at temperature {} K".format(self.current_temp))
                break

            log.info('current simulation duration is {:0.2} fs\n'.format(float(simulation_time)))

            # create a new scratch parser to parse the scratch
            # from $QCSCRATCH/AIMD
            scr_parser = MdScratchParser(self.tpd_job_name)

            job.set_velocities(scr_parser.get_velocities())

            job.molecule.positions = scr_parser.get_positions() * physical_constants.angstrom_to_bohr

            self.scratch_generator(job.temperatures_list,
                                   job.trajectory,
                                   scr_parser.get_velocities(),
                                   temp,
                                   simulation_time)

            simulation_time += self.time_step * self.aimd_steps / physical_constants.atomic_unit_of_time_to_femtosec

            if first_run:
                first_run = False
                job.rm_rem("aimd_init_veloc")

    def scratch_generator(self, inst_temperature, trj_data, velocities, temperature, sim_time):

        # append TRJ data to file
        with open(self.tpd_job_name + ".trj", 'a') as f:

            # counts the steps of the simulation
            step_counter = 0

            # counts the lines of xyz data which includes:
            # first line is number of atoms
            # second line comment
            # rest of the lines xyz of atoms
            # total of Natoms + 2 lines ( +1 if you count from 0)
            lines_counter = 0

            for line in trj_data.splitlines():

                if step_counter >= 1:

                    if lines_counter == 0:
                        lines_counter += 1
                        f.write(line + '\n')

                    elif lines_counter == 1:
                        lines_counter += 1
                        f.write('T {} K, time {:0.2f} fs\n'.format(temperature, sim_time +
                                                                   step_counter *
                                                                   self.time_step /
                                                                   physical_constants.atomic_unit_of_time_to_femtosec))

                    elif lines_counter % (self.molecule.atom_count + 1) == 0:
                        lines_counter = 0
                        step_counter += 1
                        f.write(line + '\n')

                    else:
                        lines_counter += 1
                        f.write(line + '\n')

                # the first xyz data is just a starting
                # single-point calculation which is of no interest
                # and is therefore discarded
                elif step_counter == 0 and lines_counter < self.molecule.atom_count + 1:
                    lines_counter += 1

                elif step_counter == 0 and lines_counter == self.molecule.atom_count + 1:
                    lines_counter = 0
                    step_counter += 1

        # generate file with temperatures
        with open(self.tpd_job_name + ".tempra", 'a') as f:
            f.write('Target Temperature {} K\n'.format(self.current_temp))
            f.write(inst_temperature)

        # since we're interested only in the last couple of lines
        # the data is now parsed into lines
        velocities = velocities.splitlines()
        trj_data = trj_data.splitlines()

        with open(self.tpd_job_name + ".last_pos_vel", 'a') as f:

            f.write('Target Temperature {} K\n'.format(self.current_temp))
            f.write('Ended Q-Chem aimd run with the following positions and velocities:\n')

            for i in reversed(range(1, self.molecule.atom_count + 1)):
                f.write('{}          {}\n'.format(trj_data[-1 * i], velocities[-1 * i]))
