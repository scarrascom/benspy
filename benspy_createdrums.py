#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""

Filename : bens_createdrums.py
 Purpose : Create drum files from waveforms, to be opened by DrumViewer
  Author : Sebastian Carrasco
   Email : acarrasc@uni-koeln.de

"""

import argparse
import json
import os
import shutil
import sys
import tqdm
import numpy as np

from glob import glob

from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client

def receive():
    #######################
    desc   = "Create drum files from SDS archived mseed data for visualization in Bensberg Drum Viewer"
    #def_archive = '/home/seismo/archive'
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_day = 'yesterday'
    def_comp = 'Z'
    def_stas = 'all'
    def_nets = 'all'
    def_nost = 'NONE'
    def_dict = '/mnt/Station/station/Inventory/json/drumplot.json'
    def_outd = '/mnt/Station/data'
    #def_resamp = 50
    parser = argparse.ArgumentParser(description=desc, 
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a', '--archive', type=str, 
                        default=def_archive,
                        help=f"SDS archive where the waveform data is being stored\n\
[Default: {def_archive}]")
    parser.add_argument('-d', '--day', type=str, 
                        default=def_day, 
                        help=f"Day to create the drum files in format YYYY-MM-DD, also possible 'today' or 'yesterday'\n\
[Default: {def_day}].")
    parser.add_argument('-c', '--comp', type=str,
                        default=def_comp,
                        help=f"Channel to use to create drum files\n\
[Default: {def_comp}]")
    parser.add_argument('-s', '--stations', type=str, 
                        default=def_stas,
                        help=f"Create drum files for these stations only, comma-separated [e.g. STB,HOBG,BNS]\n\
[Default: {def_stas}].")
    parser.add_argument('-n', '--networks', type=str, 
                        default=def_nets,
                        help=f"Create drums files for these networks only, comma-separated [bens,sefonib,rhb]\n\
[Default: {def_nets}].")
    parser.add_argument('-ns', '--nostations', type=str, 
                        default=def_nost,
                        help=f"Do not store these stations, comma separated [e.g. NAST,RODG]\n\
[Default: {def_nost}].")
    parser.add_argument('-ed', '--extradict', type=str, 
                        default=def_dict,
                        help=f"Dictionary to use for station specifications\n\
[Default: {def_dict}]")
    parser.add_argument('-o', '--outputdir', type=str, 
                        default=def_outd,
                        help=f"Output directory where store the drum files\n\
[Default: {def_outd}].")
    #parser.add_argument('-r', '--resample', type=float,
    #                    default=def_resamp,
    #                    help=f"Sampling rate to downsample the data\n\
#[Default: {def_resamp} Hz]")
#    parser.add_argument('-nc', '--nocopy', action='store_true', 
#                        help='Do not copy the files to the permanent archive.')
#    parser.add_argument('-nd', '--nodrum', action='store_true', 
#                        help='Do not create drum files.')
#    parser.add_argument('-soh', '--state', action='store_true', 
#                        help='Copy only state of health files.')
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
    LAPSE = 1-1e-5          # Infinitesimaly shorter than 1 s.
    internal_name, fmin, fmax, _, _ = details
    
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
                     if (key.split('.')[1] in _stal or value[0].upper() in _stal)}
    elif _netl[0].upper()!='ALL' and _stal[0].upper()=='ALL':
        new_dict = { key:value for key, value in base_dict.items() 
                    if value[3].upper() in _netl and key.split('.')[1] not in _nostal }
    else:
        new_dict = { key:base_dict[key] for key, values in base_dict.items() 
                     if ( values[0].upper() in _stal or values[3].upper() in _netl ) and values[0].upper() not in _nostal }
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
            sys.exit('Time string is wrong. Format should be YYYY-MM-DD.')
    return target

"""
def copyday(tday, dictio, soh=False,
            perm_archive='/mnt/Projects/SDSarchive',
            base_archive='/home/seismo/seiscomp/var/lib/archive/'):
    year = tday.year
    jday = tday.julday
    for key, details in dictio.items():
        network, station, location, channel = key.split('.')
        soh_code = 'D0' if soh else ''
        search = glob(os.path.join(base_archive, f'{year}', network, 
                                   station, '*', f'*{soh_code}*{jday}'))
        for f in search:
            extension = f.split(base_archive)[-1]
            dest_file = os.path.join(perm_archive, extension)
            print(dest_file)
            os.makedirs(os.path.dirname(dest_file), exist_ok=True)
            shutil.copy(f, dest_file)
"""

def create_drum(day, dictio, comp, archive, output_dir):
    sdsc = Client(archive)
    for key, details in dictio.items():
        print(key)
        network, station, location, chan = key.split('.')
        channel = f'{chan}{comp}'
        start = UTCDateTime(day.year, day.month, day.day)
        aux_start = start - 10*60
        aux_endt = start + 24*3600 + 10*60
        ## TOFIX: InternalMSEEDERROR MESSAGES BECAUSE OF SPIKES
        try:
            st = sdsc.get_waveforms(network, station, location, channel, 
                                    aux_start, aux_endt, merge=None)
        except:
            ## READ THE WAVEFORM FILE IN TWO-MIN PACKETS
            from obspy import Stream
            aux_st = Stream()
            for j in range(int((aux_endt - aux_start)/120)+1):
                aux_t0 = aux_start + j*120
                aux_tf = aux_start + (j+1)*120
                try:
                    aux_st += sdsc.get_waveforms(network, station, location, channel, 
                                                 aux_t0, aux_tf, merge=None)
                except:
                    print(f'Corrupt data or no data at {aux_t0}')
            
            aux_st.merge()
            st = aux_st.split()
            
        if st:
            drum_list = data_drum(st, details, day)
            if details[0] == 'BA09':
                station = 'BA09'
            file_name = f'DR_{start.year:04}{start.month:02}{start.day:02}.{station}'
            path_drum_file = os.path.join(output_dir, details[3], 
                                          'drum', f'{day.month:02}', file_name) 
            data2file(drum_list, path_drum_file, start, station)
            unix2dos = f'sudo unix2dos -n {path_drum_file} {path_drum_file}'
            print(unix2dos)
            os.system(unix2dos)
        else:
            print('There is no waveform file for this stream.')
            continue

if __name__ == "__main__":
    args = receive()
    date = args.day
    archive = args.archive
    stations = args.stations
    component = args.comp
    output_drum = args.outputdir
    sta_dictfile = args.extradict
    stationsl = [ s.upper() for s in (args.stations).split(',') ]
    networksl = [ n.upper() for n in (args.networks).split(',') ]
    nostatnsl = [ ns.upper() for ns in (args.nostations).split(',') ]
    ## Boolean vars
    #nocopy = args.nocopy
    #nodrum = args.nodrum
    #soh = args.state
    ## Initialize parameters
    sta_dict_drum = json.load(open(sta_dictfile))
    new_dict = crea_newdict(sta_dict_drum, stationsl, networksl, nostatnsl)
    print(new_dict)
    target_day = input2date(date)
    #if not nocopy:
    #    copyday(target_day, new_dict, soh, perm_archive='/mnt/Projects/SDSarchive')
    ## Create drum files if desired
    #if not nodrum:
    create_drum(target_day, new_dict, component, archive, output_drum)
