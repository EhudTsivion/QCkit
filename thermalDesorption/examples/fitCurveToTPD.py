from QCkit.thermalDesorption import trjProcessor

"""
This examples reads data from a directory which contain TPD data
and fits a normal distributions for the time which the molecules were
desorbed.

"""

di = trjProcessor.DirectoryInformation(
    'C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2_weaker2_friction')

# this returns the loc and scale of the normal distribution
curve_fit = di.fit_normal_distribution(plot=True,
                                   desorption_thresh=6.50,
                                   bins=30,
                                   normed=True)
print(curve_fit)
