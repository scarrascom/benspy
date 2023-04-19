#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 08:22:03 2020

@author: ascarrasco
"""

import argparse
import os

from obspy import read, UTCDateTime

def inputs():
    desc = 'Merge data from old 2-minutes files and new Centaur data.'
    parser = argparse.ArgumentParser(description=desc, 
                            formatter_class=argparse.RawTextHelpFormatter)
    def_permanent = '/mnt/Projects/SDSarchive'
    def_seedlink = '/mnt/SC-Share/seiscomp/var/lib/archive'
    parser.add_argument('date', type=str, 
              help='Day of interest to retrieve data [YYYYMMDD]')
    parser.add_argument('stationid', type=str,
                        help='ID of new station [e.g. B4.HNK..EH] (without component!)')
    parser.add_argument('-pa', '--parchive', type=str, default=def_permanent,
                        help=f'Path to the SDS directory where continuous waveforms will be stored permanently\n\
[Default: {def_permanent}]')
    parser.add_argument('-rta', '--rtarchive', type=str, default=def_seedlink, 
                        help=f'Path to the real-time seedlink archive.\n\
[Default: {def_seedlink}]')
    arg = parser.parse_args()
    return arg

if __name__=="__main__":
    inputs = inputs()
    station_id = inputs.stationid
    date = inputs.date
    permanent_archive = inputs.parchive
    seedlink_archive = inputs.rtarchive
    
    net, sta, loc, cha = station_id.split('.')
    this_date = UTCDateTime.strptime(date, '%Y%m%d')
    
    for comp in 'ZNE':
        print(f'Merging {station_id}{comp}')
        this_file = f'{net}.{sta}.{loc}.{cha}{comp}.D.{this_date.year}.{this_date.julday:03}'
        old_path = os.path.join(permanent_archive, f'{this_date.year}', 
                                net, sta, loc, cha+comp+'.D', this_file)
        new_path = os.path.join(seedlink_archive, f'{this_date.year}', 
                                net, sta, loc, cha+comp+'.D', this_file)
        new_stream = read(old_path) + read(new_path)
        print(f'Saving {old_path}')
        new_stream.write(old_path, format='MSEED', reclen=512, encoding='STEIM2')
        print(f'Saving {new_path}')
        new_stream.write(new_path, format='MSEED', reclen=512, encoding='STEIM2')
