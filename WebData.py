import requests
import re
from datetime import datetime

__author__ = 'Stuart Cossar'

class WebData(object):
    """
    Class of methods to get data from the web such as solar radio flux:
    http://services.swpc.noaa.gov/text/current-space-weather-indices.txt
    """
    def __init__(self):
        self.noaa_solar_radio_flux_txt = 'http://services.swpc.noaa.gov/text/current-space-weather-indices.txt'
        #self.noaa_solar_radio_flux_txt = ''

    def get_noa_web_data(self):
        """

        :return:
        """
        try:
            srf_data = requests.get(self.noaa_solar_radio_flux_txt).text
            return srf_data
        except:
            return False

    def noaa_to_dict(self, noaa_srf_data_file):
        """
        Convert the noaa text file to python_dict for ease of use
        :param noaa_srf_data_file: The text file from noaa
        :return: Dict
        """
        # get date: ":Solar_Radio_Flux: (\d+)\W([a-z, A-Z]+\W(\d)+)"
        re_pattern = re.compile(':Solar_Radio_Flux:(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)')
        # turn data in to list of the rows
        noaa_regex_groups_raw = re.findall(re_pattern, noaa_srf_data_file)[0]
        # print(noaa_regex_groups)
        # strip '#' from list
        noaa_regex_groups = [s.strip('#') for s in noaa_regex_groups_raw]
        # parse the data from the rows
        date_of_data = noaa_regex_groups[0]
        sites = noaa_regex_groups[1].split()
        times = noaa_regex_groups[2].split()
        freq_1_row = noaa_regex_groups[3].split()
        freq_2_row = noaa_regex_groups[4].split()
        freq_3_row = noaa_regex_groups[5].split()
        freq_4_row = noaa_regex_groups[6].split()
        freq_5_row = noaa_regex_groups[7].split()
        freq_6_row = noaa_regex_groups[8].split()
        freq_7_row = noaa_regex_groups[9].split()
        freq_8_row = noaa_regex_groups[10].split()
        freq_9_row = noaa_regex_groups[11].split()

        noaa_dict = {
            'date': datetime.strptime(date_of_data, ' %Y %b %d').date(),
            'data_by_frequency': {
                freq_1_row[0]: freq_1_row[1:],
                freq_2_row[0]: freq_2_row[1:],
                freq_3_row[0]: freq_3_row[1:],
                freq_4_row[0]: freq_4_row[1:],
                freq_5_row[0]: freq_5_row[1:],
                freq_6_row[0]: freq_6_row[1:],
                freq_7_row[0]: freq_7_row[1:],
                freq_8_row[0]: freq_8_row[1:],
                freq_9_row[0]: freq_9_row[1:],
            },
            'time_list': [datetime.strptime(date_of_data + ' ' + str(x[:2]) + ':' + str(x[2:]), ' %Y %b %d %H:%M') for x in times],
            'measurement_site_list': sites
        }

        # print(noaa_dict)

        return noaa_dict

    def solar_radio_flux_data_noaa(self):
        """
        """
        srf_data = self.get_noa_web_data()
        parse_noaa_data = self.noaa_to_dict(srf_data)
        return parse_noaa_data

if __name__ == "__main__":
    print(WebData().solar_radio_flux_data_noaa())