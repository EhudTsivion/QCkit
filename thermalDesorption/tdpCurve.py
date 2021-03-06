import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from QCkit import physical_constants as const


class TpdCurve:
    """
    A class for generation of temperature programmed desorption
    curves. These curves are the flux measured by the mass spectrometer
    of molecules emitted from the surface as the temperature is increased

    Take note, that the enthalpy supplied as input is the enthalpy at 298.15 K
    or room temperature enthalpy. The changes in the pV=RT expression are
    taken into account.

    We assumes that contributions to enthalpy from internal vibrations is constant
    in temperature.

    """

    def __init__(self,
                 entropy,
                 enthalpy_atRoomT,
                 heating_rate,
                 initial_temperature=0,
                 final_temperature=500,
                 grid_size=1000):
        """

        Calculate the normalized desorption curve as a function of temperature

        :param entropy: The entropy in kJ mol-1 K-1 units
        :param enthalpy_atRoomT: the enthalpy in kJ mol-1 units, where T = 298.15K or 25 C
        :param heating_rate: the rate of heating in K femtosec-1
        :param initial_temperature: The initial temperature of the simulation
        :param final_temperature: The final temperature of the simulation
        :param grid_size: the number of points used to evaluate the curve
        """

        self.entropy = entropy
        self.enthalpy = enthalpy_atRoomT

        self.occupation_curve = np.zeros(grid_size, dtype=np.float)
        self.desorption_curve = np.zeros(grid_size, dtype=np.float)
        self.temperatures = np.arange(grid_size)

        # the starting point is that all sites are occupied
        self.occupation_curve[0] = 1

        # convert units to sec
        heating_rate /= 1E-15

        # the duration of each simulation step
        step_duration = (final_temperature - initial_temperature) / grid_size / heating_rate

        enthalpy_zeroT = enthalpy_atRoomT - const.molar_gas_constant * 298.15 / 1000

        # create an array of temperatures
        self.temperatures = np.linspace(initial_temperature, final_temperature, grid_size)

        # create array of enthalpy
        enthalpy = np.linspace(initial_temperature, final_temperature,
                                     grid_size) * const.molar_gas_constant / 1000 + enthalpy_zeroT

        for i in np.arange(1, grid_size):
            temperature = self.temperatures[i]

            rate = rate_tst(temperature, enthalpy=enthalpy[i], entropy=entropy)

            self.occupation_curve[i] = self.occupation_curve[i - 1] * concentration(rate, step_duration)

        self.desorption_curve = np.gradient(-1 * self.occupation_curve)

        # now it is calculated
        self.desorption_curve /= self.desorption_curve.max()

        # temperature of maximum desorption
        self.max = self.temperatures[np.argmax(self.desorption_curve)]

        # this variable holds all the temperatures where he desorption
        # curve is greater than 0.5
        temperature_range = self.temperatures[np.where(self.desorption_curve > 0.5)]

        # width at half maximum, in temperature
        self.fwhm = temperature_range[-1] - temperature_range[0]


def rate_tst(temperature, enthalpy, entropy):
    """
    get the reaction rate according to transition-state-theory

    :param temperature: the temperature in K
    :param enthalpy: the enthalpy in kJ mol-1
    :param entropy: the entropy in kJ mol-1 K-1

    :return: the rate in units of sec-1
    """

    gas_constant = np.float(const.molar_gas_constant / 1000)  # in kJ/mol

    # check the units here, the're OK
    pre_exp = const.Boltzmann_constant * temperature / const.Planck_constant

    rate = pre_exp * np.exp(entropy / gas_constant) * np.exp(
        (-1 * enthalpy) / (gas_constant * temperature))  # divide by 1000 for kj/mol

    return np.float(rate)


def concentration(rate, delta_time):
    return np.exp(-1 * rate * delta_time)


if __name__ == "__main__":
    factor = 1

    curve1 = TpdCurve(entropy=0.043, enthalpy_atRoomT=10.5, initial_temperature=0, final_temperature=400,
                      heating_rate=0.0223)
    plt.plot(curve1.temperatures, curve1.desorption_curve, 'b')

    curve2 = TpdCurve(entropy=0.023, enthalpy_atRoomT=10.5, initial_temperature=0, final_temperature=400,
                      heating_rate=0.0223 * factor)
    plt.plot(curve2.temperatures, curve2.desorption_curve, 'r')

    curve3 = TpdCurve(entropy=0.013, enthalpy_atRoomT=10.5, initial_temperature=0, final_temperature=400,
                      heating_rate=0.0223 / factor)
    plt.plot(curve3.temperatures, curve3.desorption_curve, 'g')

    plt.show()
