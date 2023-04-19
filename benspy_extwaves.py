#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 23:35:03 2020

@author: ascarrasco
"""

import os
import argparse
import sys
import tqdm

import numpy as np

from obspy import UTCDateTime, Trace, Stream
from obspy.clients.fdsn import Client
from obspy.clients.filesystem.sds import Client as sdsClient
from obspy.clients.fdsn.client import FDSNNoDataException

def inputs():
    desc = 'Retrieve EXTERNAL waveform data from drum select file according to a\
            given date.'
    parser = argparse.ArgumentParser(description=desc, 
                            formatter_class=argparse.RawTextHelpFormatter)
    select_dir = '/mnt/Station/data/select'
    def_dirout = '/mnt/SeisAn/Seismo/WAV'
    def_group = 'bens'
    def_stations = 'ALL'
    parser.add_argument('date', type=str, 
              help='Day of interest to retrieve data [YYYY-MM-DD]')
    parser.add_argument('-g', '--group', type=str, default=def_group,
                        help=f'Drums to use [bens,sefo,rhb]\n\
[Default: {def_group}]')
    parser.add_argument('-s', '--stations', type=str, default=def_stations,
                        help=f'Stations to retrieve, comma-separated [e.g. MEM, AHRW]\n\
[Default: {def_stations}]')
    parser.add_argument('-ch', '--choose', action='store_true', 
                        help='Choose specific event(s) from drum list (needs input)')
    parser.add_argument('-se', '--select', type=str, default=select_dir, 
                        help=f'Path to the directory containing the drum-select files.\n\
[Default: {select_dir}]')
    parser.add_argument('-o', '--outdir', type=str, default=def_dirout, 
                        help=f'Path to output directory where to save the waveform data.\n\
[Default: {def_dirout}]')
    arg = parser.parse_args()
    return arg

def aux_stream(times):   
    auxdata = np.asarray([0,1])
    st = Stream()
    for time in times:
        tr = Trace(data=auxdata, header={'delta': 60, 'npts': 2, 
                                         'starttime': time})
        st.append(tr)
    
    st.merge()
    finalst = st.split()
    
    return finalst

if __name__=="__main__":
    inputs = inputs()
    date = inputs.date
    netw = inputs.group
    select = inputs.select
    outdir = inputs.outdir
    choose = inputs.choose
    str_stations = inputs.stations
    
    stations = [ sta.upper() for sta in str_stations.split(',')]
    print('Stations')
    print(stations)
        
    date_str = ''.join(date.split('-'))
    select_file = f'{date_str}{netw.lower()}.select'
    selectf = os.path.join(select, select_file)
    
    if not os.path.isfile(selectf):
        sys.exit(f'File {selectf} was not found')
    
    filedata = np.loadtxt(open(selectf), dtype='object')
    times = np.asarray([ UTCDateTime.strptime(date, '%Y-%m-%d-%H-%M-%SS') 
                            for date in filedata[:,1] ])
    
    new_stream = aux_stream(times)
    
    if choose:
        for k, tra in enumerate(new_stream):
            startt = tra.stats.starttime.strftime('%Y-%m-%d %H:%M:%S')
            print(f' #{k+1:3d} | {startt}')
        inputnm = input('Choose the event(s) you want to retrieve (comma-separated, e.g., 2,4,7)\n')
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

    raspi_server = 'https://data.raspberryshake.org'
    sta_dict = {'KNMI': {'NL.MAME..HH', 'NL.ROLD..HH'}, 
                raspi_server: {'AM.R5B44.00.EH', # Ex-Monschau
                               'AM.R653E.00.SH', # Mayen
                               'AM.RA52C.00.EH', # Ex-Bad Breisig
                               'AM.RFA17.00.SH', # W Koenigswinter
                               'AM.R8795.00.EH', # Bornheim
                               'AM.R3B38.00.EH', # W Taunusstein
                               'AM.R28DF.00.EH', # Bad Breisig
                               'AM.R79BC.00.EH', # Nastaetten
                               'AM.R1055.00.EH', # W Schlangenbad
                               'AM.R00A4.00.EH', # Sittard (NL)
                               'AM.R5D58.00.EH', # Alsdorf
                               'AM.R4C26.00.EH', # Wuerselen
                               'AM.R9335.00.EH', # Rott
                               'AM.R81B0.00.EH'} # Bonn
                }
    
    be_stas = ['BE.HOU..HH', 'BE.MRG..HH', 'BE.MEM..HH', 'BE.BEBN..HH', 
               'BE.RCHB..HH', 'BE.OPTB..HH']
    
    bgr_stas = ['GR.AHRW..HH', 'GR.BUG..HH', 'GR.KAST..HH', 
                'GR.MILB..HH', 'GR.TNS..HH']
    
    hs_stas = ['HS.GWBC..HH', 'HS.GWBD..HH', 'HS.GWBE..HH', 
               'HS.GWBO..HH', 'HS.WBFO..HH']
    
    lgb_stas = ['LE.MUEZ..HH', 'LE.LAGB..HH', 'LE.BIW..EH', 'LE.GLOK..HH']
    
    lux_stas = ['LU.KLB..HH', 'LU.VIA..HH', 'LU.WILW..HH']
    
    gd_stas = ['NH.ACN.00.EH', 'NH.BHE.00.SH', 'NH.BHE.00.EH', 'NH.ENT.00.SH', 
               'NH.ENTS.00.EH', 'NH.GSH.00.SH', 'NH.GSH.00.EH', 'NH.HES.00.EH', 
               'NH.HES.00.SH', 'NH.JCK.01.EH', 'NH.OLF.00.EH', 'NH.OLFT.00.EH', 
               'NH.PLH.00.SH', 'NH.PLH.00.EH', 'NH.RWB.00.EH', 'NH.SOR.00.SH', 
               'NH.SORT.00.EH', 'NH.TDN.00.SH', 'NH.TDN.00.EH', 'NH.WBS.00.EH', 
               'NH.XAN.00.EH', 'NH.XANT.00.EH']   
   
    nl_stas = ['NL.TERZ.00.HH', 'NL.TERZ.01.HH', 'NL.HRKB..BH', 'NL.VKB..HH', 
               'NL.OPLO.01.HH', 'NL.HGN.02.HH']
    
    rub_stas = ['RN.HMES..EH',  'RN.BPFI..EH',  'RN.ZER1..EH',
                'YD.LUN39..EH', 'YD.WER38..EH', 'YD.HAM29..EH']
    
    wei_stas = ['ZB.CB21..HH', 'ZB.CB22..HH', 'ZB.TB16..HH', 'ZB.TB17..HH', 
                'ZB.TB18..HH', 'ZB.TB19..HH', 'ZB.TB20..HH']
    
    bns_slink = be_stas + bgr_stas + hs_stas + lgb_stas + lux_stas + gd_stas +\
                nl_stas + rub_stas + wei_stas
    
    
    sds = sdsClient('/mnt/SC-Share/seiscomp/var/lib/archive')

    dictout = {'encoding': 'STEIM2', 'record_length': 512, 
               'dataquality': 'D', 'byteorder': '<' }
    
    nevents = len(new_stream)
    time_tqdm = tqdm.tqdm(new_stream)
    for k, tra in enumerate(time_tqdm):
        ti = tra.stats.starttime
        tf = tra.stats.endtime + 60
        strti = ti.strftime('%Y-%m-%d-%H%M')
        if k in numevs:
            print(f'### {strti} ({k+1}/{nevents}) ###')
            this_st = Stream()
            for server in sta_dict.keys():
                try:
                    cl = Client(server)
                except:
                    print('Server error')
                    continue
                for station in sta_dict[server]:
                    net, sta, loc, cha = station.split('.')
                    if sta in stations or stations[0]=='ALL':
                        print(station)
                        try:
                            st = cl.get_waveforms(net, sta, loc, f'{cha}?', 
                                                  starttime=ti, endtime=tf)
                            st.merge(method=1, fill_value='latest', interpolation_samples=0)
                            this_st += st
                        except:
                            print('No data or error!!')
            ## Retrieve stations stored in local SDSarchive (from Wiechert)
            for sl_sta in bns_slink:
                net, sta, loc, cha = sl_sta.split('.')
                if sta in stations or stations[0]=='ALL':
                    print(sl_sta)
                    slst = sds.get_waveforms(network=net, station=sta, location=loc, channel=f'{cha}[ZNE12]', 
                                             starttime=ti, endtime=tf)
                    if len(slst)==0:
                        print('No data or error!!')
                    this_st += slst
            try:
                filename = os.path.join(outdir, f'{strti}.EXT.mseed')
                this_st.write(filename, **dictout)
                print(filename)
            except:
                print(f'Error while saving file {strti} (File is empty?)')
