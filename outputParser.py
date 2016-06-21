class OutputParser:
    """
    class for parsing Q-Chem output
    For the time being, only properties of interest are implemented
    """

    def __init__(self, output_file):

        with open(output_file, 'r') as f:
            self.content = f.readlines()

    @property
    def temperature_list(self):

        """
        extract the temperatures from each md simulation step
        It's in Kelvins

        :return:
        """

        temp_list = []

        for line in self.content:

            if "Instantaneous Temperature" in line:
                temp_list.append(str(line.split()[-2]))

        return "\n".join(temp_list) + "\n"

    @property
    def job_failed(self):

        """

        :return: True if job didn't finish properly, otherwise False
        """

        flag = True

        for line in self.content:

            if "Thank you very much for using Q-Chem.  Have a nice day" in line:
                flag = False

        return flag
