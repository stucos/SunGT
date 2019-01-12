"""
Classes for calculating antenna G/T from measurements of the sun

Using th work in paper:
SSC18-X-05
Determination of Earth Station Antenna G/T Using the Sun or the Moon as an RF Source

"""

from scipy.constants import Boltzmann as k
from scipy.constants import speed_of_light as c
from numpy import log10, pi, square, power, mean, array
from datetime import datetime
from WebData import WebData

__author__ = 'Stuart Cossar'

class SunGT(object):
    """

    """
    def __init__(self, measurement_frequency, beam_correction_factor, atmospheric_attenuation, measurement_time, flux_indices=None):
        """

        :param noise_delta: Delta of source noise power density to cold sky power density, linear no units
        :param flux_indices: a two member dict for flux indicies for two frequencies {freq_1:fi_1, frq_2:fi_2}
        :param wavelength: wavelength in metres
        :param beam_correction_factor: Beam Correction Factor, linear no units
        :param atmospheric_attenuation: Atmospheric Attenuation at elevation angle, linear no units
        """
        self.noise_delta = None
        # flux_density: Solar or Lunar Flux Density, expressed in Solar Flux Units. can be entered manually or retrieved automatically
        self.fd = None
        self.flux_indices = flux_indices
        self.wavelength = self.calculate_wavelength(measurement_frequency)
        self.bcf = beam_correction_factor
        self.atmospheric_atten = atmospheric_attenuation
        self.measurement_frequency = measurement_frequency
        self.measurement_time = measurement_time

    def g_over_t(self):
        """
        Calculate and return G/T - equation (1)
        :return: G/T
        """
        gt = 10*log10((8*pi*k*(self.noise_delta-1))/(self.fd*square(self.wavelength)*self.bcf*self.atmospheric_atten))
        return gt

    def calculate_noise_delta(self, source_power, cold_sky_power):
        """
        Calculate the noise delta - equation (3) and (4)
        The value y that is required for G/T calculation, represents the difference of received power from the RF
        source versus cold sky. This value is obtained by the same means regardless of the RF emissions source
        chosen. The antenna should be pointed at the RF emissions source and a received power level recorded in
        decibels (source_power). Then the antenna should be pointed to cold sky (at the same elevation angle) and again
        the received power level recorded in decibels (cold_sky_power). The difference between the two is delta_power
        which is taken out of the decibel scale to obtain a linear value.

        for best results, delta power should be > 1dB

        :param source_power: Power in dB received by the antenna system
        :param cold_sky_power: Power in dB received by the antenna system when pointed at cold sky, at the same
        elevation as source_power
        :return: noise_delta
        """

        delta_power = float(source_power) - float(cold_sky_power)
        self.noise_delta = power(10, (delta_power/10))
        return  self.noise_delta

    def calculate_wavelength(self, measurement_frequency):
        """
        Calculate wavelength from frequency (in Hz) and speed of light - equation (5)

        :param measurement_frequency: in Hz
        :return: wavelength in meters
        """

        self.wavelength = c / measurement_frequency
        return self.wavelength

    def calculate_solar_flux_density(self):
        """
        Calculate the flux density - equations (7) and (8)
        :param method:
        :return:
        """
        # check to see if we have specified a flux density, if not, get data from noaa
        if self.flux_indices is None:
            self.flux_indices = self.get_flux_indices_noaa()

        flux_freq_1 = float(list(self.flux_indices.keys())[0])
        flux_indices_1 = float(list(self.flux_indices.values())[0])
        flux_freq_2 = float(list(self.flux_indices.keys())[1])
        flux_indices_2 = float(list(self.flux_indices.values())[1])

        # calculate Interpolation Exponent
        interpolation_exponent = (log10((self.measurement_frequency/flux_freq_2)))/(log10(flux_freq_1/flux_freq_2))

        self.sfd = flux_indices_2*power((flux_indices_1/flux_indices_2), interpolation_exponent)

        print("Solar Flux Density: %s" % self.sfd)
        return self.sfd

    def get_flux_indices_noaa(self):
        """
        The industry standard method for obtaining a value for
        solar flux density when a direct measurement is not
        possible is to interpolate a value for your desired
        frequency from the table of measured values published
        by the National Oceanic and Atmospheric
        Administration (NOAA) Space Weather Prediction
        Center. NOAA has stations scattered across the Earth
        that measure solar flux density at local noon at 9 different
        frequencies from 245 MHz to 15400 MHz
        :return: dict for flux indices for the two frequencies {freq_1:fi_1, frq_2:fi_2}
        """
        # get all NOAA data
        noaa_data = WebData().solar_radio_flux_data_noaa()

        # parse the data to get frequencies below and above our test frequency
        freqs_lower = []
        freqs_higher = []
        for frequency, flux_data_list in noaa_data['data_by_frequency'].items():
            if float(frequency) <= self.measurement_frequency:
                freqs_lower.append(int(frequency))
            elif float(frequency) >= self.measurement_frequency:
                freqs_higher.append(int(frequency))

        # get max lower freq and min higher freq to get the two either side of ours
        freq_1 = max(freqs_lower)
        freq_2 = min(freqs_higher)

        # get the flux data for the two frequencies for the time closest to our measurement time
        flux_time = min(noaa_data['time_list'], key=lambda x: abs(x - datetime.strptime(self.measurement_time, '%Y %b %d %H:%M')))
        # get the index of the time from the times list. as there are multiple meaurments at a single time, we need to get them all
        time_element_list = [i for i, e in enumerate(noaa_data['time_list']) if e == flux_time]

        # get the data for the times in the list for the frequency
        flux_values_1 = []
        flux_values_2 = []

        for sample_time in time_element_list:
            flux_value_1 = noaa_data['data_by_frequency'][str(freq_1)][sample_time]
            flux_value_2 = noaa_data['data_by_frequency'][str(freq_2)][sample_time]
            # do not add if the data is -1
            if int(flux_value_1) != -1:
                flux_values_1.append(float(flux_value_1))
            if int(flux_value_2) != -1:
                flux_values_2.append(float(flux_value_2))

        print("Flux values for %s MHz: %s @ %s" % (freq_1, flux_values_1, [noaa_data['time_list'][x] for x in time_element_list]))
        print("Flux values for %s MHz: %s @ %s" % (freq_2, flux_values_2, [noaa_data['time_list'][x] for x in time_element_list]))

        return {freq_1: mean(array(flux_values_1)), freq_2: mean(array(flux_values_2))}

if __name__ == "__main__":

    measurement_frequency = 8000
    beam_correction_factor = 0
    atmospheric_attenuation = 0
    measurement_time = '2019 Jan 12 12:23'

    gt_obj = SunGT(measurement_frequency, beam_correction_factor, atmospheric_attenuation, measurement_time)
    gt_obj.calculate_solar_flux_density()