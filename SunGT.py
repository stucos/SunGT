"""
Classes for calculating antenna G/T from measurements of the sun

Using th work in paper:
SSC18-X-05
Determination of Earth Station Antenna G/T Using the Sun or the Moon as an RF Source

"""

from scipy.constants import Boltzmann as k
from scipy.constants import speed_of_light as c
from numpy import log10, pi, square, power, mean, array, exp, sin
from datetime import datetime
from WebData import WebData

__author__ = 'Stuart Cossar'

class SunGT(object):
    """

    """
    def __init__(self):
        """

        :param noise_delta: Delta of source noise power density to cold sky power density, linear no units
        :param flux_indices: a two member dict for flux indicies for two frequencies {freq_1:fi_1, frq_2:fi_2}
        :param wavelength: wavelength in metres
        """
        """
        The apparent angular diameter used in equation 14 is the
        angular diameter of the object from an observer on the
        Earth’s surface. This is also known as optical HPBW.
        For the Sun, this value is typically defined as 0.525
        degrees.
        """
        self.hpbw = 0.525

    def g_over_t(self):
        """
        Calculate and return G/T - equation (1)
        :return: G/T
        """
        gt = 10*log10((8*pi*k*(self.noise_delta-1))/(self.fd*square(self.wavelength)*self.bcf*self.atmospheric_atten))
        return gt

    def get_noise_delta(self, source_power, cold_sky_power):
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

    def get_wavelength(self, measurement_frequency):
        """
        Calculate wavelength from frequency (in Hz) and speed of light - equation (5)

        :param measurement_frequency: in Hz
        :return: wavelength in meters
        """

        self.wavelength = c / measurement_frequency
        return self.wavelength

    def get_solar_flux_density(self, measurement_frequency, flux_indices):
        """
        Calculate the flux density - equations (7) and (8)
        :param method:
        :return:
        """
        # check to see if we have specified a flux density, if not, get data from noaa
        if flux_indices is None:
            flux_indices = self.get_flux_indices_noaa()

        flux_freq_1 = float(list(flux_indices.keys())[0])
        flux_indices_1 = float(list(flux_indices.values())[0])
        flux_freq_2 = float(list(flux_indices.keys())[1])
        flux_indices_2 = float(list(flux_indices.values())[1])

        # calculate Interpolation Exponent
        interpolation_exponent = (log10((measurement_frequency/flux_freq_2)))/(log10(flux_freq_1/flux_freq_2))

        self.sfd = flux_indices_2*power((flux_indices_1/flux_indices_2), interpolation_exponent)

        return self.sfd

    def get_beamwidth(self, wavelength, antenna_diameter):
        """
        Equation (15)

        :return: Float - Beamwidth in degrees
        """
        self.beamwidth = 68.0 * (wavelength / antenna_diameter)
        return self.beamwidth

    def get_effective_rf_diameter(self, measurement_frequency, sun_hpbw=0.525):
        """
        equation (14)
        Effective RF Diameter is a function of Beamwidth.

        :return:
        """
        self.effective_rf_diameter = sun_hpbw * (1.24 - (0.162 * log10((measurement_frequency/1000000000))))
        return self.effective_rf_diameter

    def get_beam_correction_factor(self, effective_rf_diameter, beamwidth):
        """
        equation (13)
        # TODO: answer does not match example

        :return:
        """
        self.beam_correction_factor = (1.0-(exp(-((effective_rf_diameter/beamwidth)**2)*log10(2.0))))/(((effective_rf_diameter/beamwidth)**2)*log10(2.0))

        return self.beam_correction_factor

    def get_atmospheric_attenuation(self):
        """
        Atmospheric attenuation is defined as the reduction in
        power density of an electromagnetic wave as it
        propagates through space, specifically referring to losses
        due to atmospheric composition, elevation angle, clouds,
        rain, and barometric pressure.

        The expressions used to calculate losses due to
        atmospheric attenuation are the same whether using the
        Sun or the Moon as an RF source and are based on Annex
        2 of ITU-R P.676-11 [8].

        :return:
        """
        pass

    def get_slant_path_attenuation(self, specific_atten_oxygen, equiv_height_dry_air, specific_atten_water_vapor,
                               equiv_height_water_vapor, elevation_angle):
        """
        Equation (16)

        :param specific_atten_oxygen:
        :param equiv_height_dry_air:
        :param specific_atten_water_vapor:
        :param equiv_height_water_vapor:
        :param elevation_angle:
        :return:
        """

        self.slant_path_attenuation = ((specific_atten_oxygen*equiv_height_dry_air) +
                                       (specific_atten_water_vapor*equiv_height_water_vapor))/sin(elevation_angle)
        return self.slant_path_attenuation

    def get_dry_air_specific_attenuation(self, measurement_frequency, i_oxygen_line_strength, oxygen_line_shape_factor,
                                         dry_air_continuum):
        """
        equation (17)

        :param measurement_frequency:
        :param i_oxygen_line_strength:
        :param oxygen_line_shape_factor:
        :param dry_air_continuum:
        :return:
        """

        self.dry_air_specific_attenuation = 0.1820*measurement_frequency*(i_oxygen_line_strength*oxygen_line_shape_factor + dry_air_continuum)
        return self.dry_air_specific_attenuation

    def get_water_vapor_specific_attenuation(self, measurement_frequency, i_water_vapor_line_strength,
                                             water_vapor_line_shape_factor):
        """
        equation (18)

        :param measurement_frequency:
        :param i_oxygen_line_strength:
        :param oxygen_line_shape_factor:
        :param dry_air_continuum:
        :return:
        """

        self.water_vapor_specific_attenuation = 0.1820*measurement_frequency*(i_water_vapor_line_strength*water_vapor_line_shape_factor)
        return self.water_vapor_specific_attenuation

    def get_oxygen_line_strength(self, itu_attenuation_1, itu_attenuation_2, dry_air_pressure, temperature_k):
        """
        equation (19-1)

        :param itu_attenuation_1:
        :param itu_attenuation_2:
        :param dry_air_pressure:
        :param temperature_k:
        :return:
        """

        self.oxygen_line_strength = itu_attenuation_1 * power(10, -7)*power((dry_air_pressure*(300.0/temperature_k)), 3)*(itu_attenuation_2*(1-(300.0/temperature_k)))
        return self.oxygen_line_strength

    def get_water_line_strength(self, itu_attenuation_1, itu_attenuation_2, water_vapor_pressure, temperature_k):
        """
        equation (19-2)
        :param itu_attenuation_1:
        :param itu_attenuation_2:
        :param water_vapor_pressure:
        :param temperature_k:
        :return:
        "
        self.water_line_strength = itu_attenuation_1 * power(10, -1) * power(
            (water_vapor_pressure * (300.0 / temperature_k)), 3.5) * (itu_attenuation_2 * (1 - (300.0 / temperature_k)))
        return self.water_line_strength

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

        # print("Flux values for %s MHz: %s @ %s" % (freq_1, flux_values_1, [noaa_data['time_list'][x] for x in time_element_list]))
        # print("Flux values for %s MHz: %s @ %s" % (freq_2, flux_values_2, [noaa_data['time_list'][x] for x in time_element_list]))

        return {freq_1: mean(array(flux_values_1)), freq_2: mean(array(flux_values_2))}

if __name__ == "__main__":

    antenna_diameter = 3.66 # meters
    measurement_frequency = 8200000000 # Hz
    measurement_time = '2019 Jan 12 12:23'
    flux_indices = {4995.0: 109, 8800.0: 235}

    gt_obj = SunGT()
    print('Solar FLux Density: %s SFU' % gt_obj.get_solar_flux_density(measurement_frequency, flux_indices))
    wavelength = gt_obj.get_wavelength(measurement_frequency)
    print("Wavelength: %s meters" % wavelength)
    beamwidth = gt_obj.get_beamwidth(wavelength, antenna_diameter)
    print("Beamwidth: %s degrees" % beamwidth)
    sun_effective_rf_diameter = gt_obj.get_effective_rf_diameter(measurement_frequency)
    print("Sun Effective RF diameter: %s degrees" % sun_effective_rf_diameter)
    beam_correction_factor = gt_obj.get_beam_correction_factor(sun_effective_rf_diameter, beamwidth)
    print("Beam correction factor: %s" % beam_correction_factor)
