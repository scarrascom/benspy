#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Script for archiving GD data

 Filename : bens_gdwaves.py
  Purpose : Retrieve custom data from GD files (event) to archive it.
   Author : Sebastian Carrasco
    Email : acarrasc@uni-koeln.de
"""

import argparse
import json
import os
import re
import sys

import numpy as np
import benspy_toarchive as bw

from glob import glob
from obspy.core import read, UTCDateTime, Stream

def receive():
    #######################
    desc = "Get waveforms of GD network and store them in an SDSarchive (event files only)."
    parser = argparse.ArgumentParser(description=desc, 
                                formatter_class=argparse.RawTextHelpFormatter)
    def_acq = '/mnt/SeisAn/Seismo/WAV'
    def_archive = '/mnt/Projects/EVTarchive'
    def_stas = 'all'
    def_nosta = 'NONE'
    def_dict = '/mnt/Station/Inventory/json/gd_stations.json'
    parser.add_argument('tini_date', type=str, 
                        help='Date of starting time in format YYYY-MM-DD [e.g. 2019-03-20] at 00:00')
    parser.add_argument('tend_date', type=str, 
                        help='Date of ending time in format YYYY-MM-DD [e.g. 2019-03-21] at 00:00')
    parser.add_argument('-q', '--acquisition', type=str, 
                        default=def_acq,
                        help=f"Directory where the original GD data is stored (usually WAV)\n \
[Default: {def_acq}]")
    parser.add_argument('-a', '--archive', type=str, 
                        default=def_archive,
                        help=f"Directory where the waveform data will be stored\n \
[Default: {def_archive}]")
    parser.add_argument('-s', '--stations', type=str, 
                        default=def_stas,
                        help=f"Store only these stations, comma-separated [e.g. GSH,ACN,WBS, etc]\n \
[Default: {def_stas}]")
    parser.add_argument('-ns', '--nostations', type=str, 
                        default=def_nosta,
                        help=f"Stations to not be stored, comma separated [e.g. BHE,ACN]\n \
[Default: {def_nosta}]")
    parser.add_argument('-ed', '--extradict', type=str, 
                        default=def_dict,
                        help=f"Dictionary with station specifications\n \
[Default: {def_dict}]")
    arg = parser.parse_args()
    return arg

if __name__ == "__main__":
    ## Receiving inputs
    inputs = receive()
    acq_dir = inputs.acquisition
    dict_stations = json.load(open(inputs.extradict))
    tini_date = inputs.tini_date
    tend_date = inputs.tend_date
    sdsarchi = inputs.archive
    stationsl = [ s.upper() for s in (inputs.stations).split(',') ]
    networksl = ['ALL']
    nostatnsl = [ ns.upper() for ns in (inputs.nostations).split(',') ]
    re_date = re.compile('.*-.*-.')
    dict_station = bw.crea_newdict(dict_stations, stationsl, networksl, nostatnsl)
    print("These are the selected stations:")
    print(dict_station.keys())
    ## Checking given initial and end time
    if not re_date.match(tini_date) or not re_date.match(tend_date):
        sys.exit('Date must be in YY-MM-DD format')
    
    ## Get the dates and search for the files
    dict_out = {'format': 'MSEED', 
                'reclen' : 512, 
                'encoding' : 'STEIM1', 
                'dataquality' : 'D',
                'byteorder' : '>'}
    
    starttime = UTCDateTime(tini_date)
    endtime = UTCDateTime(tend_date)
    day2sec = 3600*24
    auxendday = starttime + day2sec
    sub_wavgd = os.path.join(acq_dir, 'BNS__/*/*/20*GD*_0??')
    sub_wavgsh = os.path.join(acq_dir, 'BNS__/*/*/20*GSH*_0??')
    main_wavgd = os.path.join(acq_dir, '20*GD*_???')
    main_wavgsh = os.path.join(acq_dir, '20*GSH*_???')
    all_data = glob(main_wavgd) + glob(main_wavgsh) + glob(sub_wavgd) + glob(sub_wavgsh)
    while auxendday<=endtime:
        year = str(starttime.year)
        julday = str(starttime.julday).zfill(3)
        this_day = starttime.strftime('%Y-%m-%d')
        pre_day = (starttime-1).strftime('%Y-%m-%d')
        after_day = auxendday.strftime('%Y-%m-%d')
        print("Retrieving data for day %s" % (this_day))
        data_date = [ dates for dates in all_data 
                         if dates.split('/')[-1][:10]==this_day or 
                             dates.split('/')[-1][:10]==pre_day or 
                             dates.split('/')[-1][:10]==after_day ]
        day_stream = Stream()
        print(data_date)
        if not data_date:
            print(f'No data on {this_day}')
            starttime = auxendday
            auxendday = starttime + day2sec
            continue
        for files in data_date:
            day_stream += read(files)
        day_stream.trim(starttime, auxendday)
        day_stream.merge()
        dstream = day_stream.split().copy()
        avail_sta = np.unique([ tra.stats.station for tra in dstream ])
        for sta_name in avail_sta:
            if sta_name not in dict_station.keys():
                continue
            param_list = dict_station[sta_name]
            sta_stream = dstream.select(station=sta_name, channel='??[ZNE12]').copy()
            if not sta_stream:
                continue
            print(f'Station {sta_name}')
            for params in param_list:
                network = params[1]
                channel = params[2]
                locid = params[3]
                for tra in sta_stream:
                    tra.stats.network = network
                    tra.stats.station = sta_name
                    tra.stats.location = locid
                    newchan = channel + tra.stats.channel[-1]
                    tra.stats.channel = newchan
                    
                for comp in 'ZNE12':
                    compst = sta_stream.select(channel=f'??{comp}').copy()
                    if not compst:
                        continue
                    nett = compst[0].stats.network
                    staa = compst[0].stats.station
                    locc = compst[0].stats.location
                    chann = compst[0].stats.channel
                    chanid = f'{chann}.D'
                    fname = f'{nett}.{staa}.{locc}.{chanid}.{year}.{julday}'
                    fdir = os.path.join(sdsarchi, f'{year}/{nett}/{staa}/{chanid}')
                    os.makedirs(fdir, exist_ok=True)
                    fdirname = os.path.join(fdir, fname)
                    compst.write(fdirname, **dict_out)
        starttime = auxendday
        auxendday = starttime + day2sec
