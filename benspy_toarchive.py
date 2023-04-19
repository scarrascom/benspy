#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
#-----------------------------------
# Filename : distribute_waveforms.py
#  Purpose : Read waveform data from BNS archive 
#            directories and save into SDS archive format (compatible with SeisComp3)
#   Author : Sebastian Carrasco
#    Email : acarrasc@uni-koeln.de
#
#-----------------------------------

import numpy as np
import argparse
import re
import os
import json
import sys

from obspy.core import read, UTCDateTime, Stream

def receive():
    #######################
    script_name = __file__
    desc = f"Get waveforms of BENS stations and store them in an SDS archive.\n\
For the discrete-event case: \n \
> {script_name} tini tdate -q /mnt/Station/data -a /mnt/Projects/EVTarchive -ed /mnt/Station/station/Inventory/json/all_bns_stations_evt.json -in evt"
    def_acq = '/mnt/ACQ/work'
    def_archive = '/mnt/Projects/SDSarchive'
    def_stas = 'all'
    def_nets = 'all'
    def_nosta = 'NONE'
    def_dict = '/mnt/Station/Inventory/json/all_bns_stations.json'
    def_mode = 'cont'
    parser = argparse.ArgumentParser(description=desc, 
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tini_date', type=str, 
                        help='Date of starting time in format YYYY-MM-DD [e.g. 2019-03-20] at 00:00')
    parser.add_argument('tend_date', type=str, 
                        help='Date of ending time in format YYYY-MM-DD [e.g. 2019-03-21] at 00:00')
    parser.add_argument('-q', '--acquisition', type=str, 
                        default=def_acq,
                        help=f"Directory where the 2-min files are currently stored\n\
[Default: {def_acq}]")
    parser.add_argument('-a', '--archive', type=str, 
                        default=def_archive,
                        help=f"Directory where the waveform data will be stored [SDS]\n\
[Default: {def_archive}]")
    parser.add_argument('-s', '--stations', type=str,
                        default=def_stas,
                        help=f"Store only these stations, comma-separated [e.g. STB,HOB,BNS or all]\n\
[Default: {def_stas}]")
    parser.add_argument('-n', '--networks', type=str,
                        default=def_nets,
                        help=f"Store only these networks, comma-separated [bns,sefo,rhb,all]\n\
[Default: {def_nets}]")
    parser.add_argument('-ns', '--nostations', type=str, 
                        default=def_nosta,
                        help=f"Do not store these stations, comma separated [e.g. NAST,ROD]\n\
[Default: {def_nosta}]")
    parser.add_argument('-ed', '--extradict', type=str, 
                        default=def_dict,
                        help=f"Use this dictionary for station specifications\n\
[Default: {def_dict}]")
    parser.add_argument('-in', '--inputmode', type=str, 
                        default=def_mode, 
                        help=f"Choose between continuous (cont) or event (evt) data mode.\n\
[Default: {def_mode}]")
    arg = parser.parse_args()
    return arg


def get_archive_cont(acquisition_, 
                     _time, 
                     _statn):
    now = UTCDateTime().replace(minute=0, second=0, microsecond=0)
    path_time = os.path.join(str(_time.year), str(_time.month).zfill(2))
    hours_dt = (now - _time.replace(minute=0, second=0, microsecond=0)) / 3600.
    bns_sta  = ('BGG', 'BNS', 'DRE', 'HIL', 'HOB', 'JUE', 'KLL', 
                'KOE', 'LAU', 'NAST', 'ROD',  'STB', 'SYW7')
    sefo_sta = ('BA01', 'BA02', 'BA03', 'BA04', 'BA05', 'BA06', 'BA07', 'BA08', 'BA09', 'BA10', 'BA11', 'BA13', 'BA16')
    # Archive depends on the selected and current date (arc: only last 99 days)
    if _statn in bns_sta:
        if hours_dt <= 23:
            archive = ''
        elif 23 < hours_dt <= 2400:
            archive = 'arc/'        
        elif hours_dt > 2400:
            archive = os.path.join('archiev', path_time)
    elif _statn in sefo_sta:
        if hours_dt <= 23:
            archive = ''
        elif 23 < hours_dt:
            archive = 'arc/'
    else:
        archive = ''
    return os.path.join(acquisition_, _statn , archive)

def get_archive_evt(acquisition_, 
                    _time, 
                    _statn):
    evt_year = str(_time.year)
    evt_month = str(_time.month).zfill(2)
    if evt_year=='2023':
        archive = os.path.join('wav')
    else:
        archive = os.path.join(evt_year, 'wav', evt_month)
    return os.path.join(acquisition_, archive)

def get_extension_cont(group_, 
                       station_, 
                       year_):
    ## After 2020-08-13T10:00, BA05 started to receive only 3 channels
    #if station_=='BA05':
    #    nchans = '?'
    #else:
    #    nchans = 3
    nchans = 3
    if station_ in ('BA21','BA22'):
        nchans = 4
    main_ext = station_ + ('______00%s'%nchans)[len(station_):]
    if group_=='bns' or station_ in ('BA04', 'BA09', 'BA12', 'BA14'):
        extension = [ main_ext + '', main_ext + '.gz']
    else:
        extension = [ main_ext + '.gz', main_ext + '']
    return extension

def get_extension_evt(group_, 
                      station_, 
                      year_):
    nchans = 3
#    if station_=='BA05':
#        nchans = 4
#    else:
#        nchans = 3
    if station_ in ('BA21', 'BA22'):
        nchans = 4
    main_ext = station_ + ('______00%s'%nchans)[len(station_):]
    if year_=='2021':
        extension = [ main_ext + '', main_ext + '.gz']
    else:
        extension = [ main_ext + '.gz', main_ext + '']
    return extension

def _interp_fix(minsec, 
                maxsec, 
                this_npts,
                nominal_samples,
                delta_case_var, 
                delta_nomi_var, 
                data2interp, 
                raw_trace):
    """
    Documentation of this function
    """   
    tdelta_aux = np.linspace(minsec, maxsec, this_npts+1)
    this_delta = np.take(tdelta_aux, range(len(tdelta_aux)+delta_case_var))
    
    ndelta_aux = np.linspace(minsec, maxsec, nominal_samples+1)
    delta_nom = np.take(ndelta_aux, range(len(ndelta_aux)+delta_nomi_var))
    
    interp_data = np.interp(delta_nom, this_delta, data2interp)
    resam2int = [round(val) for val in interp_data]
    raw_trace.data = np.asarray(resam2int).astype(np.int32)
    raw_trace.stats.sampling_rate = nominal_samples/(maxsec-minsec)

def interpolate_function(nom_sr, 
                         nominal_samp, 
                         raw_trace, 
                         next_stream,
                         maxisec=120):
    """
    Documentation of this function
    """
    this_npts = raw_trace.stats.npts
    this_channel = raw_trace.stats.channel
    if this_npts==nominal_samp:
        raw_trace.stats.sampling_rate = nom_sr
    elif this_npts>nominal_samp :
        delta_case = -1
        delta_nomi = -1
        print("There are more samples than expected! Downsampling...")
        _interp_fix(0, maxisec, this_npts, nominal_samp, 
                           delta_case, delta_nomi, raw_trace.data, raw_trace)
    elif this_npts<nominal_samp:
        print("There are less samples than expected! Interpolating...")
        try:
            delta_case = 0
            delta_nomi = -1
            to_interp = np.concatenate([raw_trace.data, 
                            next_stream.select(channel=this_channel)[0].data[0:1]])
            _interp_fix(0, maxisec, this_npts, nominal_samp, 
                           delta_case, delta_nomi, to_interp, raw_trace)
        except IndexError:
            print("Stream is empty, interpolating the trace itself")
            delta_case = -1
            delta_nomi = -2
            _interp_fix(0, maxisec, this_npts, nominal_samp, 
                           delta_case, delta_nomi, raw_trace.data, raw_trace)

def check_anomalous_record(_stream, 
                           _acqui_dir, 
                           _exte, get_arc):
    """
    Interpolate and resample all the anomalous records with N++ or N--
    samples (N is expected npts)
    
    : _stream : Stream object containing the 2-min traces.
    : acqui_dir : Full path to the directory where waveforms are stored
    : exte : 2-elements tuple containing primary and secondary extension of files
             (usually 'STA___003.gz' and 'STA___003')
    """
    min2sec = 120
    _nom_sr = round(_stream[0].stats.sampling_rate)
    nominal_samp = int(_nom_sr*min2sec)
    if _stream[0].stats.npts<nominal_samp:
        print("There are less samples than expected, need to get next trace!")
        next_stream = get_this_stream( _stream[0].stats.starttime + min2sec, 
                                    _stream[0].stats.station, _acqui_dir, 
                                    _exte, get_arc)
    else:
        next_stream = Stream()
    
    for _trace in _stream:
        print(_trace)
        interpolate_function(_nom_sr, nominal_samp, _trace, 
                             next_stream,  min2sec)
    return _stream

def get_this_stream(_this_time,
                    _station_,
                    _acq_dir_,
                    _extn, 
                    get_arc):
    """
    Function to retrieve file stream for one specific date.
    
    : _this_time : UTCDateTime object containing the starting time of stream
    : _station_ :  String object, it corresponds to the station name.
    : _acq_dir_ :  String object, it corresponds to the root archive directory
    : _extn : 2-elements tuple containing primary and secondary extension of files
              (usually 'STA___003.gz' and 'STA___003')
    """
    
    file_suffix = _this_time.strftime('%Y-%m-%d-%H%M') + '-00S'
    arc_path = get_arc(_acq_dir_, _this_time, _station_)
    print( "Searching in directory: %s" % arc_path )
    print("#############################################")
    try:
        print("Reading file %s.%s" % (file_suffix, _extn[0]))
        auxst = read(os.path.join(arc_path, file_suffix + '.' +_extn[0]))
    except (Exception, IOError, TypeError, AssertionError) as error:
        print(error)
        try:
            print("Trying with %s extension" % _extn[1])
            auxst = read(os.path.join(arc_path, file_suffix + '.' +_extn[1]))
        except (Exception, IOError, TypeError, AssertionError) as new_error:
            print(new_error)
            print("There are no waveform files for this time, \
                    returning empty stream")
            auxst = Stream()
    return auxst

def get_stream(_station, 
               _aux_list_times,
               _acq_dir,
               _ext, 
               getarc):
    st  = Stream()
    for this_time in _aux_list_times:
        this_stream = get_this_stream(this_time, _station, _acq_dir, 
                                      _ext, getarc)
        if this_stream:
            checked_stream = check_anomalous_record(this_stream, _acq_dir, 
                                                    _ext, getarc)
            st += checked_stream
    return st

def get2min_time(_utctime):
    aux_time = _utctime
    aux_time.second = 0
    aux_time.microsecond = 0
    if aux_time.minute % 2:
        aux_time.minute = aux_time.minute-1
    return aux_time

def get2min_array(_starttime, _endtime):
    every = 120.                       # Every 2 minutes
    aux_start = get2min_time(_starttime)
    aux_end = get2min_time(_endtime)
    n2minbins = int( ( aux_end - aux_start )/ every )
    bins2min = [ aux_start + val * every for val in range(n2minbins) ]
    return bins2min

def save_stream(network, station, loccode, 
                            channel, year, julday, sdsarchi, fstream):
    dict_out = {'reclen' : 512, 
                'encoding' : 'STEIM2', 
                'dataquality' : 'D',
                'byteorder' : '>'}
    for comp in 'ZNE':
        chanid = channel + comp
        name_file = ("%s.%s.%s.%s.D.%s.%s"
                    % (network, station, loccode, chanid, year, julday))
        dir_path = os.path.join(sdsarchi, year, network, station, chanid+".D")
        full_path = os.path.join(dir_path, name_file)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        if fstream.select(channel=chanid):
            print(full_path)
            fstream.select(channel=chanid).write(full_path, 
                        format='MSEED', **dict_out)
        else:
            print(f"There are no waveforms for this day for station {station} and component {chanid}")

def crea_newdict(base_dict, _stal, _netl, _nostal):
    if _stal[0].upper()=='ALL' and _netl[0].upper()=='ALL':
        new_dict = base_dict
    elif _netl[0].upper()=='ALL' and _stal[0].upper()!='ALL':
        new_dict = { x:base_dict[x] for x in _stal }
    elif _netl[0].upper()!='ALL' and _stal[0].upper()=='ALL':
        new_dict = { key:base_dict[key] for x in _netl 
                            for key, values in base_dict.items() 
                            for value in values
                                if value[4].lower()==x and key not in _nostal }
    else:
        new_dict = { key:base_dict[key] for x in _netl 
                        for key, values in base_dict.items() 
                        for value in values
                            if ( value[4].lower()==x or key in _stal ) and key not in _nostal } 
    return new_dict

if __name__ == "__main__":
    ## Receiving inputs
    inputs = receive()
    acq_dir = inputs.acquisition
    dict_stations = json.load(open(inputs.extradict))
    tini_date = inputs.tini_date
    tend_date = inputs.tend_date
    inmode = inputs.inputmode
    sdsarchi = inputs.archive
    stationsl = [ s.upper() for s in (inputs.stations).split(',') ]
    networksl = [ n.lower() for n in (inputs.networks).split(',') ]
    nostatnsl = [ ns.upper() for ns in (inputs.nostations).split(',') ]
    re_date = re.compile('.*-.*-.')
    # Create new dictionary with desired stations or networks
    dict_station = crea_newdict(dict_stations, stationsl, networksl, nostatnsl)
    print("These are the selected stations:")
    print(dict_station.keys())
    ## Checking given initial and end time
    if not re_date.match(tini_date) or not re_date.match(tend_date):
        sys.exit('Date must be in YY-MM-DD format')
    
    dict_inmode = {'evt': [get_extension_evt, get_archive_evt], 
                   'cont': [get_extension_cont, get_archive_cont]}
    get_extension, get_archive = dict_inmode[inmode]
    ## Get the dates and search for the files
    starttime = UTCDateTime(tini_date)
    endtime = UTCDateTime(tend_date)
    day2sec = 3600*24
    auxendday = starttime + day2sec
    while auxendday<=endtime:
        year = str(starttime.year)
        julday = str(starttime.julday).zfill(3)
        this_day = starttime.strftime('%Y-%m-%d')
        print("Retrieving data for day %s" % (this_day))
        aux_list_times = get2min_array(starttime, auxendday)
        for int_station_name, param_list in dict_station.items():
            for params in param_list:
                network = params[1]
                channel = params[2]
                loccode = params[3]
                group = params[4]
                station = params[5].upper()
                acq_sta = os.path.join(acq_dir, params[0])
                print(acq_sta)
                extension = get_extension(group, station, year)
                station_stream = get_stream(station, aux_list_times, 
                                            acq_sta, extension, get_archive)
#                print(station_stream.__str__(extended=True))
                try:
                    station_stream.merge(method=0)
                except:
                    print('Likely different sampling rates, not merging the traces...')
                    pass
                nomasked_stream = station_stream.split()
                ## Rewrite the stream with new nomenclature (according to IRIS)
                bool_ba23 = False
                for trace in nomasked_stream:
                    trace.stats.station = int_station_name
                    trace.stats.network = network                    
                    trace.stats.channel = channel + trace.stats.channel[-1]
                    trace.stats.location = loccode
                    if trace.stats.station=='BA21' and trace.stats.channel=='DNL':
                        trace.stats.station = 'BA23'
                        trace.stats.channel = 'DNN'
                        trace.data = trace.data*-1
                        bool_ba23 = True
                    elif trace.stats.station=='BA22' and trace.stats.channel=='DNT':
                        trace.stats.station = 'BA23'
                        trace.stats.channel = 'DNE'
                        bool_ba23 = True
                print("Writing output file")
                print(nomasked_stream)
                ## If BA23 is present, save it first.
                if bool_ba23:
                    print('BA23 is here!')
                    ba23st = nomasked_stream.select(station='BA23')
                    for tt in ba23st:
                        nomasked_stream.remove(tt)
                    save_stream(network, 'BA23', loccode, channel, 
                                year, julday, sdsarchi, ba23st)
                ## Saving the original station
                save_stream(network, int_station_name, loccode, channel, 
                            year, julday, sdsarchi, nomasked_stream)
        starttime = auxendday
        auxendday = starttime + day2sec
