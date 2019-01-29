"""
Classes for calculating antenna G/T from measurements of the sun

Using th work in paper:
SSC18-X-05
Determination of Earth Station Antenna G/T Using the Sun or the Moon as an RF Source

"""

import itur
from scipy.constants import Boltzmann as k
from scipy.constants import speed_of_light as c
from numpy import log10, pi,power, mean, array, exp, log
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

    def g_over_t(self, noise_delta, solar_flux_density, wavelength, beam_corr_factor, atmospheric_atten):
        """
        Calculate and return G/T - equation (1)
        :return: G/T
        """
        sfd = solar_flux_density*10**-22 # convert to correct units
        atmospheric_attenuation = power(10, (atmospheric_atten/10)) # turn to linear from dB

        gt = 10*log10(((8*pi*k*(noise_delta-1))/(sfd*(wavelength**2)*beam_corr_factor*atmospheric_attenuation)))

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

    def get_solar_flux_density(self, measurement_frequency, flux_indices, test_date, test_time):
        """
        Calculate the flux density - equations (7) and (8)
        :param method:
        :return:
        """
        # check to see if freq is above what noaa provide. if so, use quiet sun approximation
        if measurement_frequency > 15400000000:
            print("Measurement frequency above NOAA observations, using least-square fitting equation")
            solar_flux_density = gt_obj.extrapolate_solar_flux(measurement_frequency)
            return solar_flux_density

        # convert Hz to MHz
        measurement_frequency = measurement_frequency/1000000
        # check to see if we have specified a flux density, if not, get data from noaa
        if flux_indices is None:
            flux_indices = self.get_flux_indices_noaa(test_date, test_time, measurement_frequency)

        flux_freq_1 = float(list(flux_indices.keys())[0])
        flux_indices_1 = float(list(flux_indices.values())[0])
        flux_freq_2 = float(list(flux_indices.keys())[1])
        flux_indices_2 = float(list(flux_indices.values())[1])

        # print("SFD retrieved from Web: %s Mhz: %s, %s MHz: %s" % (flux_freq_1, flux_indices_1, flux_freq_2, flux_indices_2))
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
        # 68 or 70????
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
        numerator_arr = array((-(effective_rf_diameter/beamwidth)**2)*log(2.0))
        numerator = 1.0-exp(numerator_arr)
        denominator = ((effective_rf_diameter/beamwidth)**2)*log(2.0)

        self.beam_correction_factor = numerator / denominator

        # self.beam_correction_factor = 1 + 0.38*((effective_rf_diameter/beamwidth)**2)

        return self.beam_correction_factor

    def get_atmospheric_attenuation(self, lat, lon, measurement_frequency, elevation, antenna_diameter,
                                    unavailability=0.1, clear_skies=True):
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
        """
        Ag: Gasious Attenuation 
        Ac: Cloud attenuation 
        Ar: rain attenuation 
        As: Scintilation attenuation 
        Att: total attenuation
        """
        Ag, Ac, Ar, As, Att = itur.atmospheric_attenuation_slant_path(lat, lon,
                                                                      (measurement_frequency/1000000000)*itur.u.GHz,
                                                                      elevation,
                                                                      unavailability,
                                                                      antenna_diameter*itur.u.m,
                                                                      return_contributions=True)
        # if clear skies (default) return only the gasious attenuation
        if clear_skies:
            self.atmospheric_atten = Ag.value
            return Ag.value
        self.atmospheric_atten = Att.value
        return Att.value

    def extrapolate_solar_flux(self, frequency):
        """
        from "10-60Ghz G/T measurments using the sun as a source - a preliminary study" William C. Daywitt
        https://www.govinfo.gov/content/pkg/GOVPUB-C13-53a55ea34f3ca8aaacf289a9caa0bee6/pdf/GOVPUB-C13-53a55ea34f3ca8aaacf289a9caa0bee6.pdf
        equation 4

        :param frequency:
        :return:
        """
        frequency = frequency/1000000000.0
        flux_log10 = 1.20 + 1.10*(log10(frequency)) + 0.179*(log10(frequency))**2
        flux = power(10, flux_log10)

        return flux

    def get_flux_indices_noaa(self, test_date, test_time, measurement_frequency):
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
        noaa_data = WebData().solar_radio_flux_data_noaa(test_date)

        # parse the data to get frequencies below and above our test frequency
        freqs_lower = []
        freqs_higher = []
        for frequency, flux_data_list in noaa_data['data_by_frequency'].items():
            if float(frequency) <= measurement_frequency:
                freqs_lower.append(int(frequency))
            elif float(frequency) >= measurement_frequency:
                freqs_higher.append(int(frequency))

        # get max lower freq and min higher freq to get the two either side of ours
        freq_1 = max(freqs_lower)
        freq_2 = min(freqs_higher)

        # get the flux data for the two frequencies for the time closest to our measurement time
        # print(noaa_data['time_list'])
        try:
            flux_time = min(noaa_data['time_list'], key=lambda x: abs(x - datetime.strptime('%s %s' % (test_date, test_time), '%Y %b %d %H:%M')))
        except ValueError:
            print("Invalid test time")
            exit(1)
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

        print("Flux values (from NOAA web) for %s MHz: %s @ %s" % (freq_1, flux_values_1[0], [noaa_data['time_list'][x] for x in time_element_list][0].strftime('%Y-%m-%d %H:%M')))
        print("Flux values (from NOAA web) for %s MHz: %s @ %s" % (freq_2, flux_values_2[0], [noaa_data['time_list'][x] for x in time_element_list][0].strftime('%Y-%m-%d %H:%M')))

        return {freq_1: mean(array(flux_values_1)), freq_2: mean(array(flux_values_2))}

    def calculate_g_t(self, measurement_frequency, flux_indices, measurement_date, measurement_time,
                      source_power, cold_sky_power, antenna_diameter, lat, lon, elevation):
        """

        :param measurement_frequency:
        :param flux_indices:
        :param measurement_date:
        :param measurement_time:
        :param source_power:
        :param cold_sky_power:
        :param antenna_diameter:
        :param lat:
        :param lon:
        :param elevation:
        :return:
        """
        solar_flux_density = self.get_solar_flux_density(measurement_frequency, flux_indices, measurement_date,
                                                           measurement_time)
        print('Solar FLux Density: %s SFU' % solar_flux_density)
        noise_delta = self.get_noise_delta(source_power, cold_sky_power)
        print("Noise delta: %s" % noise_delta)
        wavelength = self.get_wavelength(measurement_frequency)
        print("Wavelength: %s meters" % wavelength)
        beamwidth = self.get_beamwidth(wavelength, antenna_diameter)
        print("Beamwidth: %s degrees" % beamwidth)
        sun_effective_rf_diameter = self.get_effective_rf_diameter(measurement_frequency)
        print("Sun Effective RF diameter: %s degrees" % sun_effective_rf_diameter)
        beam_correction_factor = self.get_beam_correction_factor(sun_effective_rf_diameter, beamwidth)
        print("Beam correction factor: %s" % beam_correction_factor)
        slant_path_attenuation = self.get_atmospheric_attenuation(lat, lon, measurement_frequency, elevation,
                                                                    antenna_diameter)
        print("Atmospheric attenuation: %s dB" % slant_path_attenuation)
        g_over_t = self.g_over_t(noise_delta, solar_flux_density, wavelength, beam_correction_factor,
                                   slant_path_attenuation)
        print("G/T: %s" % g_over_t)

if __name__ == "__main__":

    antenna_diameter = 2.2  # meters
    measurement_frequency = 19300000000  # Hz
    measurement_date = '2019 Jan 2'
    measurement_time = '17:00'
    #flux_indices = {8800.0: 216, 15400: 525}
    flux_indices = None
    lat = 32.2909615
    lon = 34.865974
    elevation = 38  # degrees
    source_power = -43.72  # dBm
    cold_sky_power = -55.836 # dBm

    gt_obj = SunGT()
    g_t = gt_obj.calculate_g_t(measurement_frequency, flux_indices, measurement_date, measurement_time, source_power,
                               cold_sky_power, antenna_diameter, lat, lon, elevation)
