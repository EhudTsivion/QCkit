

class OutputParser:
    """
    class for parsing Q-Chem output
    For the time being, only properties of interest are implemented
    """

    def __init__(self, output_file):

        with open(output_file, 'r') as f:

            self.content = f.readlines()

    def get_temperatures(self):

        """
        extract the temperatures from each md simulation step
        It's in Kelvins

        :return:
        """

        temp_array = []

        for line in self.content:

            if "Instantaneous Temperature" in line:
                temp_array.append(str(line.split()[-2]))

        return "\n".join(temp_array) + "\n"
