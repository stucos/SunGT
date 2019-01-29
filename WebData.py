from ftplib import FTP
import re
from datetime import datetime
from io import BytesIO

__author__ = 'Stuart Cossar'

class WebData(object):
    """
    Class of methods to get data from the web such as solar radio flux:
    http://services.swpc.noaa.gov/text/current-space-weather-indices.txt
    """
    def __init__(self):
        # self.noaa_solar_radio_flux_txt = 'http://services.swpc.noaa.gov/text/current-space-weather-indices.txt'
        self.noaa_solar_radio_flux_dir = 'ftp.swpc.noaa.gov'

    def get_noa_web_data(self):
        """

        :return:
        """
        # try:
            #ftp = FTP(self.noaa_solar_radio_flux_dir)
        ftp = FTP('ftp.swpc.noaa.gov')
        ftp.login()
        ftp.cwd('/pub/lists/radio/')
        #ftp.retrlines('LIST')
        r = BytesIO()
        ftp.retrbinary('RETR 45day_rad.txt', r.write)
        #print(r.getvalue())
        return r.getvalue().decode("utf-8")

        # except:
        #     return False

    def noaa_to_dict(self, noaa_srf_data_file, test_date):
        """
        Convert the noaa text file to python_dict for ease of use
        :param noaa_srf_data_file: The text file from noaa
        :return: Dict
        """

        re_pattern = re.compile('%s\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)\n(.+)' % test_date)
        # turn data in to list of the rows
        try:
            noaa_regex_groups = re.findall(re_pattern, noaa_srf_data_file)[0]
        except IndexError:
            print("Could not get NOAA data for given date")
            exit(0)

        # print(noaa_regex_groups)
        # parse the data from the rows
        date_of_data = test_date
        sites = ['Learmonth', 'San Vito', 'Sag Hill', 'Penticton', 'Penticton', 'Palehua', 'Penticton']
        times = ['0500',  '1200', '1700',  '1700',  '2000', '2300', '2300']
        freq_1_row = noaa_regex_groups[0].split()
        freq_2_row = noaa_regex_groups[1].split()
        freq_3_row = noaa_regex_groups[2].split()
        freq_4_row = noaa_regex_groups[3].split()
        freq_5_row = noaa_regex_groups[4].split()
        freq_6_row = noaa_regex_groups[5].split()
        freq_7_row = noaa_regex_groups[6].split()
        freq_8_row = noaa_regex_groups[7].split()
        freq_9_row = noaa_regex_groups[8].split()

        noaa_dict = {
            'date': datetime.strptime(date_of_data, '%Y %b %d').date(),
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
            'time_list': [datetime.strptime('%s %s' % (date_of_data, x), '%Y %b %d %H%M') for x in times],
            'measurement_site_list': sites
        }

        # print(noaa_dict)

        return noaa_dict

    def solar_radio_flux_data_noaa(self, test_date):
        """
        """
        srf_data = self.get_noa_web_data()
        parse_noaa_data = self.noaa_to_dict(srf_data, test_date)
        return parse_noaa_data

if __name__ == "__main__":
    print(WebData().solar_radio_flux_data_noaa('2019 Jan 2'))