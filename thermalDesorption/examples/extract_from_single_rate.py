from QCkit.thermalDesorption import trjProcessor
import numpy as np

"""
This examples reads data from two separate directories
which contain TPD data of the same system, heated in different rates
and extracts the thermodynamic properties of entropy (S) and enthalpy (H)

"""

di = trjProcessor.DirectoryInformation(
    'C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2_slow_step')

opt=False

for h in np.arange(7, 12, 0.5):
    for s in np.arange(0.01, 0.05, 0.01):
        optimization_values = di.extract_enthalpy_entropy(desorption_thresh=4.5,
                                                          guess_enthalpy=h,
                                                          guess_entropy=s,
                                                          heating_rate=0.022222222)

        print(optimization_values)
        if optimization_values['error'] < 1:
            opt = True
            print('opt!')
            break


    if opt:
        break


di.fit_normal_distribution(plot=True, normed=False)

print(optimization_values)

print('enthalpy -ΔH: ', optimization_values['enthalpy'])
print('entropy -ΔS: ', optimization_values['entropy'])
print('error: ', optimization_values['error'])
