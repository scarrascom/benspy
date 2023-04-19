#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 15:39:13 2020

Script for copying files from seiscomp archive to SDSarchive and create drums

@author: seismo
"""

import argparse
import json
import os
import subprocess
import sys
import tqdm

from obspy import UTCDateTime
#from obspy.clients.fdsn import Client
from obspy.clients.filesystem.sds import Client

def receive():
    #######################
    desc   = "Create drum files for external stations from external server for visualization in Bensberg Drum Viewer"
    parser = argparse.ArgumentParser(description=desc, 
                                formatter_class=argparse.RawTextHelpFormatter)
    def_day = 'yesterday'
    def_comp = 'Z'
    def_stas = 'AHRW,TNS,GSH,LAGB,MEM,BUG'
    def_dict = '/mnt/Station/station/Inventory/json/drumplot_ext.json'
    def_outd = '/mnt/Station/data'
    parser.add_argument('-d', '--day', type=str, 
                        default=def_day, 
                        help=f"Day to create the drum files in format YYYY-MM-DD also possible 'today' or 'yesterday'\n\
[Default: {def_day}].")
    parser.add_argument('-c', '--comp', type=str,
                        default=def_comp,
                        help=f"Channel to use to create drum files\n\
[Default: {def_comp}].")
    parser.add_argument('-s', '--stations', type=str, 
                        default=def_stas,
                        help=f"Create drum files for these stations only, comma-separated\n\
[Default: {def_stas}].")
    parser.add_argument('-ed', '--extradict', type=str, 
                        default=def_dict,
                        help=f"Use different dictionary for station specifications\n\
[Default: {def_dict}].")
    parser.add_argument('-o', '--outputdir', type=str, 
                        default=def_outd,
                        help=f"Output directory where store the drum files\n\
[Default: {def_outd}].")
    arg = parser.parse_args()
    return arg

def data_drum(stream, details, ddate):
    '''
    Compute and create drum list from data stream during one specific day.
    Parameters
    ----------
    stream : Obspy stream object.
        Data stream including data for the desired day.
    details : 3x1 list or numpy array.
        List containing station name, minimum and maximum frequency to filter.
    date : UTCDateTime object.
        Desired day to compute the drum file.

    Returns
    -------
    drum_list : Python list.
        DESCRIPTION.
    '''
    
    DAY2SEC = 86400
    TWOMINSEC = 120
    LAPSE = 1-1e-5 # Infinitesimaly shorter than 1 s.
    internal_name, fmin, fmax, _ = details
    
    stream.detrend('linear')
    stream.filter('bandpass', freqmin=fmin, freqmax=fmax, zerophase=True)

    drum_list = [ [] for gr in range(DAY2SEC//TWOMINSEC) ]
    startt = UTCDateTime(ddate.year, ddate.month, ddate.day)
    sbar = tqdm.tqdm(range(DAY2SEC))
    for second in sbar:
        sbar.set_description('Seconds')
        sliced_st = stream.slice(startt+second, startt+second+LAPSE).copy()
        if sliced_st:
            aux_max = -1e20
            aux_min = 1e20
            for tr in sliced_st:
                aux_max = tr.data.max() if tr.data.max()>=aux_max else aux_max
                aux_min = tr.data.min() if tr.data.min()<=aux_min else aux_min
        else:
            aux_max = aux_min = 0
            sbar.set_description('No data')
        index = second // TWOMINSEC
        drum_list[index].append(int(aux_max))
        drum_list[index].append(int(aux_min))

    return drum_list


def data2file(drum_data, drum_path, startt, station):
    '''
    Save drum data to drum ASCII file.

    Parameters
    ----------
    drum_data : Nx1 List or Numpy Array.
        List object containing min/max per second data per two minutes data.
    drum_path : String.
        Path to where the desired file will be saved.
    startt : UTCDateTime object.
        Starting time at midnight of the desired day.

    Returns
    -------
    None.

    '''
    
    TWOMINSEC = 120    
    COLHEAD = 73    # Length of header
    EVERY = 6       # 6 values per row
    if os.path.exists(drum_path):
        os.remove(drum_path)
    os.makedirs(os.path.dirname(drum_path), exist_ok=True)
    drum_file = open(drum_path, 'a')
    for k, segment in enumerate(drum_data):
        if not all(v==0 for v in segment): # Do not save if full of zeroes.
            aux_start = startt + k*TWOMINSEC
            day_header = f'{aux_start.year:04} {aux_start.month:02} {aux_start.day:02}'
            hour_header = f'{aux_start.hour:02} {aux_start.minute:02}'
            header = f' {station:<4}{day_header} {hour_header} {TWOMINSEC*2} 1'
            full_head = f'{header:<{COLHEAD}}' + '\n'
            drum_file.write(full_head)
            subsegs = [ segment[i:i+EVERY] for i in range(0, len(segment), EVERY) ]
            for subseg in subsegs:
                data_line = ''.join([ f'{val:>12}' for val in subseg]) + '\n'
                drum_file.write(data_line)
    drum_file.close()

def crea_newdict(base_dict, _stal, _netl, _nostal):
    '''
    Create new dictionary based on (un)selected stations/groups
    Parameters
    ----------
    base_dict : Dictionary formatted as: { 'NET.STA.LOCID.CHAN' : 
                                           ['Name', fmin, fmax, 'group']}
                Base dictionary containing information
    _stal : List
        Selected stations to be included (default is ALL)
    _netl : List
        Selected group of stations to be included (default is ALL)
    _nostal : List
        Selected stations to not be included (default is NONE)

    Returns
    -------
    new_dict : Dictionary
        Dictionary with selected stations

    '''
    if _stal[0].upper()=='ALL' and _netl[0].upper()=='ALL':
        new_dict = base_dict
    elif _netl[0].upper()=='ALL' and _stal[0].upper()!='ALL':
        new_dict = { key:value for key, value in base_dict.items() 
                     if key.split('.')[1] in _stal }
    elif _netl[0].upper()!='ALL' and _stal[0]=='ALL':
        new_dict = { key:value for key, value in base_dict.items() 
                    if value[3] in _netl and key.split('.')[1] not in _nostal }
    else:
        new_dict = { key:base_dict[key] for key, values in base_dict.items() 
                      if ( values[0] in _stal or values[3] in _netl ) and values[0] not in _nostal } 
    return new_dict

def input2date(date):
    ## Setting the date for the drum plot creation
    if date.upper() in ('TODAY', 'T', 'NOW'):
        target = UTCDateTime.now()
    elif date.upper() in ('YESTERDAY', 'YES', 'Y'):
        target = UTCDateTime.now() - 3600*24
    else:
        try:
            target = UTCDateTime(date)
        except:
            sys.exit('Time string is wrong. Format is YYYY-MM-DD (or YYYYMMDD)')
    return target


def ext_drum(day, dictio, comp, output_dir):
    sds_arc = '/mnt/SC-Share/seiscomp/var/lib/archive'
    cl = Client(sds_arc)
    for key, details in dictio.items():
        print(key)
        network, station, location, chan = key.split('.')
        channel = f'{chan}{comp}'
        start = UTCDateTime(day.year, day.month, day.day)
        aux_start = start - 5*60
        aux_endt = start + 24*3600 + 5*60
        st = cl.get_waveforms(network, station, location,
                              channel, aux_start, aux_endt)
        if not st:
            print('No waveform data for this station')
            continue
        drum_list = data_drum(st, details, day)
        file_name = f'DR_{start.year:04}{start.month:02}{start.day:02}.{station}'
        path_drum_file = os.path.join(output_dir, details[3],
                                      'drum', f'{day.month:02}', file_name)
        data2file(drum_list, path_drum_file, start, station)
        unix2dos = f'sudo unix2dos -n {path_drum_file} {path_drum_file}'
        os.system(unix2dos)
  
if __name__ == "__main__":
    args = receive()
    date = args.day
    stations = args.stations
    component = args.comp
    output_drum = args.outputdir
    sta_dictfile = args.extradict
    stationsl = [ s.upper() for s in (args.stations).split(',') ]
    ## Initialize parameters
    sta_dict_drum = json.load(open(sta_dictfile))
    new_dict = crea_newdict(sta_dict_drum, stationsl, 'ALL', '')
    target_day = input2date(date)
    ## Create drum files if desired
    ext_drum(target_day, new_dict, component, output_drum)
