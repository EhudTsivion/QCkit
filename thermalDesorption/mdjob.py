import qcjob


class MDjob(qcjob.QCjob):

    def set_velocities(self, velocities):
        """
        append velocities data to end of current job_file
        """

        # add velocities
        self.other_things = {"velocity": velocities}


