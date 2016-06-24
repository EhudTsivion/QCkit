import QCkit.qcjob as qcjob


class MDjob(qcjob.QCjob):

    def set_velocities(self, velocities):
        """
        append velocities data to end of current job_file
        """

        # add velocities
        self.other_things = {"velocity": velocities}

    def change_temperature(self, new_temperature):
        """
        change the temperature of the job
        """

        # remove old
        self.rm_rem("aimd_temp")

        # add new temperature:
        self.add_rem("aimd_temp", str(new_temperature))

    @property
    def temperatures_list(self):
        """
        Return a list of all instantaneous temperatures
        during the simulation
        """

        return self.output_parser.temperature_list

    @property
    def trajectory(self):
        """
        Return a list of all instantaneous temperatures
        during the simulation
        """

        return self.aimd_scrach_parser.trajectory
