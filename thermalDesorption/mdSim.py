import datetime
import random
import logging as log
from matplotlib import pyplot as plt
import numpy as np

from QCkit.atom import Atom
from QCkit.molecule import Molecule
from QCkit import physical_constants
from QCkit.thermalDesorption.mdjob import MDjob
from QCkit.thermalDesorption.mdScratchParser import MdScratchParser


class MD:
    """
    A class for simulation of temperature programmed desorption experiments
    using molecular dynamics
    """

    def __init__(self,
                 molecule,  # the molecule object, contains all the information about the molecule
                 basis,  # basis set to use
                 exchange,  # exchange to use (such as B3LYP etc.)
                 temperature,  # target temperature simulation
                 aimd_num_steps,  # number of MD steps
                 thermostat_timescale,  # friction of system with the thermostat head-bath

                 # time step of the MD simulation, in atomic units = 0.0242 fs, default is 0.0242 fs
                 aimd_step_duration=10,

                 thermostat="langevin",  # type of thermostat
                 threads=None,  # number of openmp threads to use. If none then is OMP_NUM_THREADS
                 job_name=None,
                 velocities=None,
                 restart=False,
                 extra_rems=None):  # name of the job. appears in all related output files.

        self.temperature = temperature
        self.molecule = molecule
        self.time_step = aimd_step_duration
        self.aimd_steps = aimd_num_steps
        self.basis = basis
        self.threads = threads
        self.exchange = exchange
        self.thermostat = thermostat
        self.thermostat_timescale = thermostat_timescale
        self.job_name = job_name
        self.velocities = velocities
        self.extra_rems = extra_rems

        if not job_name:

            self.job_name = "mdrun-QChem" + datetime.datetime.now().strftime('-%Y-%m-%d-') \
                            + str(random.randint(1000000, 9999999))

        else:

            if not restart and job_name != 'test':
                # always append a random number, because several
                # jobs with same name can run in parallel
                self.job_name = '{}-{}'.format(job_name, str(random.randint(1000000, 9999999)))

            else:
                # all output is appended
                self.job_name = job_name

        if not type(self.time_step) == int:

            print('non integer simulation time step detected, attempt to round')

            self.time_step = round(self.time_step)


        log.basicConfig(filename="{}.log".format(self.job_name),
                        filemode='a',
                        level='INFO',
                        format='')

        if restart:
            log.info("\n\n{:*^30}".format("RESTART"))

        self.run(restart)

        log.info("******** Simulation ended")

    def run(self, restart=False):

        print('called run')

        rems = {"aimd_thermostat": self.thermostat,
                "aimd_time_step": self.time_step,
                "aimd_temp": self.temperature,
                "aimd_steps": self.aimd_steps,
                "aimd_init_veloc": "thermal",
                "aimd_print": "1",
                "max_scf_cycles": "200",
                "aimd_langevin_timescale": self.thermostat_timescale}

        # if extra_rems are specified
        # add them to these rems
        if self.extra_rems:
            rems.update(self.extra_rems)

        # I think this MDjob thing is
        # probably unnecessary
        job = MDjob(molecule=self.molecule,
                    job_type='aimd',
                    basis=self.basis,
                    threads=self.threads,
                    exchange=self.exchange,
                    molden_format='False',
                    job_name=self.job_name,
                    rems=rems)

        log.info('\n\n{:*^30}'.format('new MD simulation'))
        log.info('\ncomputational details:')
        log.info('basis: {}, exchange: {}'.format(self.basis, self.exchange))
        log.info('Target temperature {}'.format(self.temperature))
        log.info('using {} thermostat'.format(job.rems["aimd_thermostat"]))
        log.info('thermostat frequency: '.format(job.rems['aimd_langevin_timescale']))

        if restart:
            job.set_velocities(self.velocities)
            job.rm_rem('aimd_init_veloc')
            first_run = False

        else:
            first_run = True

        job.run()

        if not job.failed:
            log.info('Q-Chem finished MD run successfully')

        else:
            log.info('Q-Chem job failed')

        # create a new scratch parser to parse the scratch
        # from $QCSCRATCH/AIMD
        scr_parser = MdScratchParser(self.job_name)

        # update position for next AIMD run
        job.molecule.positions = scr_parser.get_positions() * physical_constants.angstrom_to_bohr

        self.scratch_generator(job.temperatures_list,
                               job.trajectory,
                               scr_parser.get_velocities(),
                               self.temperature)

        # plot interesting quantities

        tandv = scr_parser.get_temp_and_potential

        plt.title('job: {}'.format(self.job_name))

        fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

        ax0.plot(tandv['time'], tandv['temperature'])
        ax1.plot(tandv['time'], tandv['potentialE'])

        ax0.set_title('temperature')
        ax1.set_title('potential energy')

        #
        # plt.plot(x_axis, job.potential_energy_list)
        #
        # plt.savefig(self.job_name + 'temperature.png')
        #
        # plt.
        #
        plt.savefig(self.job_name + '_temperature.png')

    def scratch_generator(self, inst_temperature, trj_data, velocities, temperature):

        # append TRJ data to file
        with open(self.job_name + ".trj", 'a') as f:

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
                        f.write('T {} K, time {:0.2f} fs\n'.format(inst_temperature[step_counter],
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
        with open(self.job_name + ".tempra", 'a') as f:
            f.write('Target Temperature {} K\n'.format(self.temperature))
            f.write(inst_temperature)

        # since we're interested only in the last couple of lines
        # the data is now parsed into lines
        velocities = velocities.splitlines()
        trj_data = trj_data.splitlines()

        with open(self.job_name + ".last_pos_vel", 'a') as f:

            f.write('Target Temperature {} K\n'.format(self.temperature))
            f.write('Ended Q-Chem aimd run with the following positions and velocities:\n')

            for i in reversed(range(1, self.molecule.atom_count + 1)):
                f.write('{}          {}\n'.format(trj_data[-1 * i], velocities[-1 * i]))


if __name__ == "__main__":
    pass
