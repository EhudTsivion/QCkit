from QCkit.thermalDesorption.mdjob import MDjob
import datetime
import random
import logging as log
import physical_constants
from outputParser import OutputParser
from thermalDesorption.mdScratchParser import MdScratchParser


class TPD:

    """
    A class for simulation of temperature programmed desorption experiments
    using molecular dynamics
    """

    def __init__(self,
                 molecule,
                 basis,
                 exchange,
                 threads,
                 high_T,
                 low_T,
                 temp_advance,
                 thermostat,
                 nh_length,
                 nh_timescale,
                 time_step,
                 aimd_steps,
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

        self.nose_hoover_length = nh_length

        self.nose_hoover_timescale = nh_timescale

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

    def run(self):

        if not self.thermostat.lower() == "nose_hoover":
            raise Exception("only nose-hoover thermostat is currently supported")

        # first run starts at the lowest temperature
        # and uses thermal initial velocities guess

        job = MDjob(molecule=self.molecule,
                    job_type="aimd",
                    basis=self.basis,
                    threads=self.threads,
                    exchange=self.exchange,
                    molden_format="False",
                    job_name=self.tpd_job_name,
                    rems={"aimd_thermostat": self.thermostat,
                          "aimd_time_step": self.time_step,
                          "aimd_temp": self.low_temperature,
                          "aimd_steps": self.aimd_steps,
                          "aimd_init_veloc": "thermal",
                          "aimd_print": "1",
                          "nose_hoover_length": self.nose_hoover_length,
                          "nose_hoover_timescale": self.nose_hoover_timescale} )

        log.info("\n\n{:*^30}".format("new TPD simulation"))
        log.info("\ncomputational details:")
        log.info("basis: {}, exchange: {}".format(self.basis, self.exchange))
        log.info("starting temperature {}, ending temperature {} K".format(self.low_temperature,
                                                                           self.high_temperature))
        log.info("advancing with steps of {} K".format(self.temp_advance))

        job.run()

        scr_parser = MdScratchParser(self.tpd_job_name)

        # create the trjectory
        scr_parser.append_trj_data(first_step=True)

        self.temprature_file_generator(new=True,
                                       tempra_list=OutputParser(self.tpd_job_name + ".qchem").
                                       get_temperatures())

        # remove aimd_init_veloc rem
        job.rm_rem("aimd_init_veloc")

        # add read guess rem
        job.add_rem("scf_guess", "read")

        # subsequent runs use the last geometry and last
        # velocities as obtained by

        for temp in range(self.low_temperature + self.temp_advance,
                          self.high_temperature,
                          self.temp_advance):

            log.info("Temperature increased to {}".format(temp))

            job.set_velocities(scr_parser.get_velocities())

            job.molecule.positions = scr_parser.get_positions()*physical_constants.angstrom_to_bohr

            job.run()

            self.temprature_file_generator(tempra_list=OutputParser(self.tpd_job_name + ".qchem").
                                           get_temperatures())

            scr_parser.append_trj_data()

    def temprature_file_generator(self, tempra_list, new=False):

        if new:

            status = 'w'

        else:

            status = 'a'

        with open(self.tpd_job_name + ".tempra", status) as f:
            f.write(tempra_list)







