import qcjob


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






