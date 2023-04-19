#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 23:35:03 2020

@author: ascarrasco
"""

import os
import sys
import argparse

import benspy_getwaves as get_wvfm
import numpy as np
import pandas as pd

from obspy import UTCDateTime, Trace, Stream

def inputs():
    desc = 'Get waveform data from drum selection for a certain date.'
    parser = argparse.ArgumentParser(description=desc, 
                            formatter_class=argparse.RawTextHelpFormatter)
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_selectdir = '/mnt/Station/data/select'
    def_dirout = '/mnt/SeisAn/Seismo/WOR'
    def_stations = 'ALL'
    def_nets = 'SAME'
    def_nost = 'NONE'
    parser.add_argument('date', type=str, 
                        help='Day of interest to retrieve data [YYYY-MM-DD]')
    parser.add_argument('group', type=str,
                        help='Drum-select file to use, from DrumViewer selection [bens,sefo,rhb]')
    parser.add_argument('-a', '--archive', type=str, 
                        default=def_archive,
                        help=f'Path to the SDS directory where continuous waveforms are stored\n\
[Default: {def_archive}].')
    parser.add_argument('-sd', '--selectdir', type=str, 
                        default=def_selectdir, 
                        help=f'Path to the directory containing the drum-select files\n\
[Default: {def_selectdir}].')
    parser.add_argument('-s', '--stations', type=str,
                        default=def_stations,
                        help=f"Stations to retrieve, comma-separated [e.g. STB,HOBG,BNS]\n\
[Default: {def_stations}].")
    parser.add_argument('-n', '--networks', type=str,
                        default=def_nets,
                        help="Retrieve this(ese) network(s) only, comma-separated [bens,sefo,rhb,lgb,gd,eifel]\n\
[Default: same as select group].")
    parser.add_argument('-ns', '--nostations', type=str,
                        default=def_nost,
                        help=f"Do not retrieve these stations, comma separated [e.g. NAST,RODG]\n\
[Default: {def_nost}].")
    parser.add_argument('-o', '--outdir', type=str, 
                        default=def_dirout, 
                        help=f'Path to output directory where to save the waveform data\n\
[Default: {def_dirout}].')
    parser.add_argument('-ch', '--choose', action='store_true',
                        help='Choose specific event(s) from drum list (needs user input)')
    parser.add_argument('-kw', '--keepwaves', action='store_true',
                        help='Keep (do not remove) old waveform files.')
    arg = parser.parse_args()
    return arg

def aux_stream(times):   
    auxdata = np.asarray([0,1])
    st = Stream()
    for time in times:
        tr = Trace(data=auxdata, header={'delta': 60, 'npts':2, 'starttime': time})
        st.append(tr)
    
    st.merge()
    finalst = st.split()
    
    return finalst

if __name__=="__main__":
    inputs = inputs()
    date = inputs.date
    selgroup = inputs.group
    archive = inputs.archive
    select = inputs.selectdir
    outdir = inputs.outdir
    choose = inputs.choose
    keep_waves = inputs.keepwaves
    ch_stations = [ s.upper() for s in (inputs.stations).split(',') ]
    networksl = [ n.upper() for n in (inputs.networks).split(',') ]
    nostatnsl = [ ns.upper() for ns in (inputs.nostations).split(',') ]
    
    date_str = ''.join(date.split('-'))
    select_file = f'{date_str}{selgroup.lower()}.select'
    selectf = os.path.join(select, select_file)
    
    if not os.path.isfile(selectf):
        sys.exit(f'File {selectf} was not found')
    
    filedata = np.loadtxt(open(selectf), dtype='object')
    date_str = '%Y-%m-%d-%H-%M-%SS'
    times = np.asarray([ UTCDateTime.strptime(datte, date_str) 
                            for datte in filedata[:,1] ])
    new_stream = aux_stream(times)
    if networksl[0]=='SAME':
        networksl = [selgroup.upper()]

    stations, nets = get_wvfm.crea_dictsta(ch_stations, networksl, ['NONE'])
    
    dictout = {'encoding': 'STEIM2', 'record_length': 512, 
                'dataquality': 'D', 'byteorder': '<' }
    
    nevents = len(new_stream)
    filenr = 'filenr.lis'
    filenrlis = os.path.join(outdir, filenr)
    
    ## Remove previous filenr.lis and waveforme files
    if os.path.exists(filenrlis):
        try:
            data = pd.read_csv(filenrlis, header=None)
            if not keep_waves:
                print('Trying to remove old waveform files...')
                for waves_name in data[0]:
                    fname = waves_name.split('  ')[-1]
                    path_wave = os.path.join(outdir, fname)
                    if os.path.exists(path_wave):
                        os.remove(path_wave)
            print('Removing previous filenr.lis and data (if existed)')
        except:
            print('filenr.lis is likely empty')
            pass
        os.remove(filenrlis)
    
    if choose:
        for k, tra in enumerate(new_stream):
            startt = tra.stats.starttime.strftime('%Y-%m-%d %H:%M:%S')
            print(f' #{k+1:3d} | {startt}')
        inputnm = input('Choose the event(s) you want to retrieve (comma-separated, ranges with - , e.g., 2-4,7)\n')
        auxnum = []
        for nums in inputnm.split(','):
            ranges = nums.split('-')
            if len(ranges)>1:
                auxnum.append(range(int(ranges[0])-1, int(ranges[-1]), 1))
            else:
                auxnum.append(int(ranges[0])-1)
        numevs = np.hstack(auxnum)
    else:
        numevs = np.arange(len(new_stream))

    with open(filenrlis, "w") as flis:
        j = 0
        for k, tra in enumerate(new_stream):
            ti = tra.stats.starttime
            tf = tra.stats.endtime + 60
            strti = ti.strftime('%Y-%m-%d-%H%M')
            if k in numevs:
                print(f'### {strti} ({k+1}/{nevents}) ###')
                this_st = get_wvfm._get_data(ti, tf, archive, 
                                             stations, 'ZNE', merge=True)
                filename = get_wvfm.out_put(ti, this_st, 'MSEED', dictout, nets, 'full', outdir)
                fname = filename.split('/')[-1]
                flis.write(f' #{j+1:3d}  {fname}\n')
                j += 1
    os.system(f'more {filenrlis}')
