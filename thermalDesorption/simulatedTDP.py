import datetime
import random
import logging as log

from QCkit.thermalDesorption.mdjob import MDjob
from QCkit import physical_constants
from QCkit.outputParser import OutputParser
from QCkit.thermalDesorption.mdScratchParser import MdScratchParser


class TPD:
    """
    A class for simulation of temperature programmed desorption experiments
    using molecular dynamics
    """

    def __init__(self,
                 molecule,
                 basis,
                 exchange,
                 high_T,
                 low_T,
                 temp_advance,
                 time_step,
                 aimd_steps,
                 thermostat_timescale,
                 thermostat="langevin",
                 threads=None,
                 tpd_job_name=None):

        self.low_temperature = low_T
        self.high_temperature = high_T
        self.temp_advance = temp_advance
        self.molecule = molecule
        self.time_step = time_step
        self.aimd_steps = aimd_steps
        self.basis = basis
        self.threads = threads
        self.exchange = exchange
        self.thermostat = thermostat
        self.thermostat_timescale = thermostat_timescale

        self.current_temp = self.low_temperature

        if not tpd_job_name:

            self.tpd_job_name = "tpd-qchem-" + datetime.datetime.now().strftime('%Y-%m-%d-') \
                                + str(random.randint(10000, 99999))

        else:
            # always append date and some random number
            self.tpd_job_name = tpd_job_name + datetime.datetime.now().strftime('%Y-%m-%d-') \
                                + str(random.randint(10000, 99999))

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

        for temp in range(self.low_temperature,
                          self.high_temperature,
                          self.temp_advance):

            job.change_temperature(temp)

            self.current_temp = temp

            log.info("\n++ starting a new Q-Chem AIMD job at {} K".format(self.current_temp))

            job.run()

            if not job.failed:
                log.info('Q-Chem done with {}'.format(self.current_temp))

            else:
                log.info("Q-Chem job failed at temperature {} K".format(self.current_temp))
                break

            # create a new scratch parser to parse the scratch
            # from $QCSCRATCH/AIMD
            scr_parser = MdScratchParser(self.tpd_job_name)

            job.set_velocities(scr_parser.get_velocities())

            job.molecule.positions = scr_parser.get_positions() * physical_constants.angstrom_to_bohr

            # self.temperature_file_generator(job.temperatures_list)
            #
            # self.trj_file_generator(job.trajectory)

            self.scratch_generator(job.temperatures_list, job.trajectory, scr_parser.get_velocities())

            if first_run:
                first_run = False
                job.rm_rem("aimd_init_veloc")

    def scratch_generator(self, temperature_data, trj_data, velocities):

        # generate TRJ file
        with open(self.tpd_job_name + ".trj", 'a') as f:
            # f.write('Target Temperature {} K\n'.format(self.current_temp))
            f.write(trj_data)

        # generate file with temperatures
        with open(self.tpd_job_name + ".tempra", 'a') as f:
            f.write('Target Temperature {} K\n'.format(self.current_temp))
            f.write(temperature_data)

        # since we're interested only in the last couple of lines
        # the data is now parsed into lines
        velocities = velocities.splitlines()
        trj_data = trj_data.splitlines()

        with open(self.tpd_job_name + ".last_pos_vel", 'a') as f:

            f.write('Target Temperature {} K\n'.format(self.current_temp))
            f.write('Ended Q-Chem aimd run with the following positions and velocities:\n')

            for i in reversed(range(1, self.molecule.atom_count + 1)):
                f.write('{}          {}\n'.format(trj_data[-1*i], velocities[-1*i]))