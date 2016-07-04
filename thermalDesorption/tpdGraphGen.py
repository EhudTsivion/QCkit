import os
import numpy as np

import matplotlib.pylab as plt
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams

from scipy.stats import norm

from QCkit.thermalDesorption.trjProcessor import DirectoryInformation


class TpdGraphGenerator:
    def __init__(self):

        self.source_dir_list = list()

    def add_source_dir(self, dir_name):
        """
        Add an information source for generating the graphs
        the source is a directory which contains JSON data files

        :param dir_name: the name of the directory to add
        """

        if os.path.isdir(dir_name):

            self.source_dir_list.append(dir_name)

        else:

            raise NotADirectoryError('This is not a directory')

        pass

    def gen_figure(self, desorption_thresh=5.5):

        if not self.source_dir_list:
            raise Exception('No information directory was provided')

        max_temp = 0

        fit_loc_array = list()

        fit_scale_array = list()

        information_array = list()

        max_temp = 0

        for info in self.source_dir_list:

            di = DirectoryInformation(info)

            dd = di.get_desorption_distribution(thresh=desorption_thresh)

            information_array.append(dd)

            if max_temp < dd.max():

                max_temp = dd.max()

            fit_loc, fit_scale = norm.fit(dd)

            fit_loc_array.append(fit_loc)

            fit_scale_array.append(fit_scale)

        max_graph_temp = max_temp + 50

        temp_axis = np.linspace(0, max_graph_temp, 200)

        counter = 0

        font = FontProperties()

        for info in information_array:

            pdf_fitted = norm.pdf(temp_axis,
                                  fit_loc_array[counter],
                                  fit_scale_array[counter])

            counter += 1

            plt.fill(temp_axis, pdf_fitted, alpha=0.6)

            # plt.hist(info, bins=15, normed=True, alpha=0.3)

        plt.xticks(np.arange(0, max_graph_temp, 100))

        plt.yticks(np.linspace(0, plt.axes().get_ylim()[1], 3))

        # plt.show()

        plt.savefig('tpd_graph.png', dpi=300, figsize=(3.33, 3), fontsize=14)


if __name__ == "__main__":
    graphs = TpdGraphGenerator()

    rcParams.update({'font.size': 22, 'xtick.major.size': 10})

    graphs.add_source_dir(
        'C:/Users/Udi-BRIX/Dropbox/abinitio/multi_h2/dynamics/production/tcatMg/tcatMg2H2_good_results')

    graphs.add_source_dir(
        'C:/Users/Udi-BRIX/Dropbox/abinitio/multi_h2/dynamics/production/tcatMg/tcatMg1H2')

    graphs.gen_figure()
