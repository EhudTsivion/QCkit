from QCkit.thermalDesorption import trjProcessor

"""
This examples reads data from two separate directories
which contain TPD data of the same system, heated in different rates
and extracts the thermodynamic properties of entropy (S) and enthalpy (H)

"""

optimized_values = trjProcessor.extract_enthalpy_entropy_two_info(
    desorption_thresh=5.5,
    guess_enthalpy=6,
    guess_entropy=0.01,
    data_dir_1='C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2_strong_friction',
    data_dir_2='C:/Users/Udi-BRIX\Dropbox/abinitio\multi_h2\dynamics\production/tcatMg/tcatMg2H2_fast',
    heating_rate1=0.0267,
    heating_rate2=0.0267 * 5 / 2)

print(optimized_values)
