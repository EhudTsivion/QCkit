from QCkit.thermalDesorption import trjProcessor

"""
This examples reads data from two separate directories
which contain TPD data of the same system, heated in different rates
and extracts the thermodynamic properties of entropy (S) and enthalpy (H)

"""

optimization_values = trjProcessor.extract_enthalpy_entropy_two_info(
    desorption_thresh=5.6,
    guess_enthalpy=2.5,
    guess_entropy=0.026,
    data_dir_1='C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2',
    data_dir_2='C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2_fast',
    heating_rate1=0.022222222,
    heating_rate2=0.022222222 * 5 / 2)

print('enthalpy ΔH: ', optimization_values['enthalpy'])
print('entropy ΔS: ', optimization_values['entropy'])
print('error: ', optimization_values['error'])
