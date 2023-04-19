#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:59:16 2021

@author: seismo
"""

import argparse
import json
import os
import re
import sys

import benspy_toarchive as bw
import pandas as pd

from obspy import read, UTCDateTime

def receive():
    #######################
    desc = "Check the times of the data gaps (missing or corrupted data)."
    parser = argparse.ArgumentParser(description=desc, 
                                formatter_class=argparse.RawTextHelpFormatter)
    def_acq = '/mnt/ACQ/work'
    def_sta = 'all'
    def_nets = 'all'
    def_nost = 'NONE'
    def_extd = '/mnt/Station/station/Inventory/json/all_bns_stations.json'
    def_outd = '/mnt/Station/data/gap_check'
    parser.add_argument('tini_date', type=str, 
                        help='Date of starting time in format YYYY-MM-DD [e.g. 2019-03-20] at 00:00')
    parser.add_argument('tend_date', type=str, 
                        help='Date of ending time in format YYYY-MM-DD [e.g. 2019-03-21] at 00:00')
    parser.add_argument('-q', '--acquisition', type=str, 
                        default=def_acq,
                        help=f"Directory where the 2-min files are currently stored\n\
[Default: {def_acq}]")
    parser.add_argument('-s', '--stations', type=str, 
                        default=def_sta,
                        help=f"Store only these stations, comma-separated [e.g. STB,HOB,BNS]\n\
[Default: {def_sta}]")
    parser.add_argument('-n', '--networks', type=str, 
                        default=def_nets,
                        help=f"Store only these networks, comma-separated [bns,sefo,rhb]\n\
[Default: {def_nets}]")
    parser.add_argument('-ns', '--nostations', type=str, 
                        default=def_nost,
                        help=f"Stations to not be stored, comma separated [e.g. NAST,ROD]\n\
[Default: {def_nost}]")
    parser.add_argument('-ed', '--extradict', type=str, 
                        default=def_extd,
                        help=f"Use different dictionary for station specifications\n\
[Default: {def_extd}]")
    parser.add_argument('-o', '--outdir', type=str, 
                        default=def_outd, 
                        help=f'Directory where to save the output files\n\
[Default: {def_outd}].')
    arg = parser.parse_args()
    return arg


def check_data(_this_time,
                _station_,
                _acq_dir_,
                _extn, 
                get_arc):
    """
    Function to check whether the data file does exist or is corrupted.
    
    : _this_time : UTCDateTime object containing the starting time of stream
    : _station_ :  String object, it corresponds to the station name.
    : _acq_dir_ :  String object, it corresponds to the root archive directory
    : _extn : 2-elements tuple containing primary and secondary extension of files
              (usually 'STA___003.gz' and 'STA___003')
    """
    
    file_suffix = _this_time.strftime('%Y-%m-%d-%H%M') + '-00S'
    arc_path = get_arc(_acq_dir_, _this_time, _station_)
    try:
        path_file = os.path.join(arc_path, file_suffix + '.' +_extn[0])
        auxst = read(path_file)
        var = 'O'
    except (Exception, FileNotFoundError, TypeError, AssertionError) as error:
        if type(error)==TypeError:
            var = 'C'
        elif type(error)==FileNotFoundError:
            try:
                aux_path = os.path.join(arc_path, file_suffix + '.' +_extn[1])
                auxst = read(aux_path)
                var = 'O'
            except (TypeError, FileNotFoundError) as new_err:
                if type(new_err)==TypeError:
                    var = 'C'
                elif type(new_err)==FileNotFoundError:
                    var = 'M'
    return var

if __name__=="__main__":
    ## Receiving inputs
    inputs = receive()
    acq_dir = inputs.acquisition
    dict_stations = json.load(open(inputs.extradict))
    tini_date = inputs.tini_date
    tend_date = inputs.tend_date
    outdir = inputs.outdir
    stationsl = [ s.upper() for s in (inputs.stations).split(',') ]
    networksl = [ n.lower() for n in (inputs.networks).split(',') ]
    nostatnsl = [ ns.upper() for ns in (inputs.nostations).split(',') ]
    re_date = re.compile('.*-.*-.')
    # Create new dictionary with desired stations or networks
    dict_station = bw.crea_newdict(dict_stations, stationsl, networksl, nostatnsl)
    print("These are the selected stations:")
    print(dict_station.keys())
    ## Checking given initial and end time
    if not re_date.match(tini_date) or not re_date.match(tend_date):
        sys.exit('Date must be in YY-MM-DD format')
    
    dict_inmode = {'evt': [bw.get_extension_evt, bw.get_archive_evt], 
                   'cont': [bw.get_extension_cont, bw.get_archive_cont]}
    get_extension, get_archive = dict_inmode['cont']
    
    ## Get the dates and search for the files
    starttime = UTCDateTime(tini_date)
    endtime = UTCDateTime(tend_date)
    day2sec = 3600*24
    auxendday = starttime + day2sec
    while auxendday<=endtime:
        mdict = {}
        year = str(starttime.year)
        julday = str(starttime.julday).zfill(3)
        this_day = starttime.strftime('%Y-%m-%d')
        print("Retrieving data for day %s" % (this_day))
        aux_list_times = bw.get2min_array(starttime, auxendday)
        for auxt in aux_list_times:
            str_auxt = auxt.strftime('%Y-%m-%d-%H%M-%SS')
            print(f'Checking time: {str_auxt}')
            mdict[str_auxt] = {}
            for int_station_name, param_list in dict_station.items():
                for params in param_list:
                    network = params[1]
                    channel = params[2]
                    loccode = params[3]
                    group = params[4]
                    station = params[5].upper()
                    acq_sta = os.path.join(acq_dir, params[0])
                    extension = get_extension(group, station, year)
                    gap_type = check_data(auxt, station, acq_sta, extension, get_archive)
                    sta_code = f'{network}.{station}.{loccode}.{channel}'
                    mdict[str_auxt] = {station: gap_type, **mdict[str_auxt]}
        
        dfcorr = pd.DataFrame.from_dict(mdict).T
        for station in dfcorr:
            ## Export the corrupted data times
            thisdf = dfcorr[station]
            bool_corr = thisdf=='C'
            if any(bool_corr):
                daystr = starttime.strftime('%Y%m%d')
                fname = f'{daystr}_{station}.corr'
                path_out = os.path.join(outdir, fname)
                with open(path_out, 'w') as f:
                    for item in thisdf[bool_corr].index:
                        f.write("%s\n" % item)
        starttime = auxendday
        auxendday = starttime + day2sec
