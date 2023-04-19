#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 18:42:54 2022

@author: ascarrasco
"""
import argparse
import os
import pytz
import sys

from obspy import read, UTCDateTime
from obspy.clients.filesystem.sds import Client

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def receive():
    #######################
    desc   = "Create daily drum plots for one specific date [0-24 h]"
    parser = argparse.ArgumentParser(description=desc,
                                formatter_class=argparse.RawTextHelpFormatter)
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_day = 'yesterday'
    def_interval = 20
    def_outd = '/mnt/Station/data/drumplots'
    def_stations = 'ALL'
    parser.add_argument('-a', '--archive', type=str,
                        default=def_archive,
                        help=f"Directory where the waveform data is stored\n\
[DEF {def_archive}]")
    parser.add_argument('-d', '--day', type=str, 
                        default=def_day,
                        help=f"Day to create the drum plots in format YYYYMMDD, also possible 'today' or 'yesterday'.\n\
[DEF {def_day}]")
    parser.add_argument('-i', '--inval', type=float,
                        default=def_interval,
                        help=f"Interval length for each row, in minutes\n\
[DEF {def_interval}]")
    parser.add_argument('-o', '--outputdir', type=str,
                        default=def_outd,
                        help=f"Output directory where store the drum files.\n\
[DEF {def_outd}]")
    parser.add_argument('-s', '--stations', type=str, 
                        default=def_stations,
                        help=f"Create drumplots for specific stations only.\n\
[DEF {def_stations}]")
    parser.add_argument('-nu', '--noupload', action='store_true',
                        help='Do not transfer/upload drum plots to website.')
    parser.add_argument('-nc', '--nocopy', action='store_true', 
                        help='Do not copy to the station version file (without date)')
    arg = parser.parse_args()
    return arg

def input2date(date):
    ## Setting the date for the drum plot creation
    if date.upper() in ('TODAY', 'T', 'NOW'):
        auxt = UTCDateTime.now()
        target = UTCDateTime(auxt.year, auxt.month, auxt.day)
    elif date.upper() in ('YESTERDAY', 'YES', 'Y'):
        auxt = UTCDateTime.now() - 3600*24
        target = UTCDateTime(auxt.year, auxt.month, auxt.day)
    else:
        try:
            target = UTCDateTime(date)
        except:
            sys.exit("Time string is wrong, format must be:\n\
                     - today/t/now\n\
                     - yesterday/yes/y or\n\
                     - Any date with format YYYYMMDD")
    return target

if __name__ == "__main__":
    args = receive()
    date = args.day
    archive = args.archive
    noup = args.noupload
    nocopy = args.nocopy
    out_dir = args.outputdir
    inval = args.inval
    stations = [ s.upper() for s in (args.stations).split(',') ]

    t0 = input2date(date)
    tf = t0 + 3600*24
    todayst = t0.strftime('%d.%m.%Y')
    
    sds = Client(archive)
    
    # Plotting settings
    xres, yres = 800, 1050
    xdim = 8
    ydim = xdim*yres/xres
    ftsize = 12
    colors = ['red', 'blue', 'black']
    dpi = 100

    dict_amp = {
                'BA02': { 'id': 'BQ.BA02..EHZ', 'amp': 160000000*1e-6},
                'BNS':  { 'id': 'BQ.BNS..HHZ',  'amp': 301720000*2e-6},
                'BGG':  { 'id': 'BQ.BGG..HHZ',  'amp': 481000000*5e-7},
                'DREG': { 'id': 'BQ.DREG..HHZ', 'amp': 600000000*1e-6},
                'HILG': { 'id': 'BQ.HILG..HHZ', 'amp': 301720000*9e-7},
                'HOBG': { 'id': 'BQ.HOBG..HHZ', 'amp': 600000000*1e-6},
                'JUE':  { 'id': 'BQ.JUE..EHZ',  'amp': 160000000*1e-5},
                'KLL':  { 'id': 'BQ.KLL..EHZ',  'amp': 160000000*1e-6},
                'KOE':  { 'id': 'BQ.KOE..HHZ',  'amp': 301720000*5e-7},
                'LAUG': { 'id': 'BQ.LAUG..EHZ', 'amp': 160000000*4e-7},
                'NAST': { 'id': 'BQ.NAST..HHZ', 'amp': 301720000*1e-6},
                'RODG': { 'id': 'BQ.RODG..EHZ', 'amp': 59600000*1e-4},
                'BD14': { 'id': 'B4.BD14..EHZ',  'amp': 160000000*5e-5},
                'BD15': { 'id': 'B4.BD15..EHZ',  'amp': 160000000*5e-5},
                'BD17': { 'id': 'B4.BD17..EHZ',  'amp': 160000000*5e-5},
                'BER':  { 'id': 'B4.BER..EHZ',  'amp': 160000000*5e-5},
                'FUN':  { 'id': 'B4.FUN..EHZ',  'amp': 160000000*5e-5},
                'HMB':  { 'id': 'B4.HMB..EHZ',  'amp': 160000000*5e-5},
                'HNK':  { 'id': 'B4.HNK..EHZ',  'amp': 160000000*5e-5},
                'MIL':  { 'id': 'B4.MIL..EHZ',  'amp': 160000000*5e-5},
                'ROE':  { 'id': 'B4.ROE..EHZ',  'amp': 160000000*5e-5},
                'SIN':  { 'id': 'B4.SIN..EHZ',  'amp': 160000000*5e-5},
                }
    
    new_dict = {}
    if 'ALL' in stations:
        new_dict = dict_amp
    else:
        for sta in stations:
            if sta not in dict_amp.keys():
                sys.exit(f'Station {sta} is not in the available drumplots')
            new_dict[sta] = dict_amp[sta]

    ## Pameters for uploading to website
    pwdtxt = '/home/seismo/scripte/ebsbg1pwd.txt'
    address = 'ebsbg1@dialog.rrz.uni-koeln.de:/afs/.rrz.uni-koeln.de/common/info/www/docs/math-nat-fak/geologie/seismo/seismogramme/'
    
    local = pytz.timezone('Europe/Berlin')
    for stan, stdic in new_dict.items():
        stid = stdic['id']
        print(f'Dayplot for {stid} ({t0.day:02d}/{t0.month:02d}/{t0.year:04d})')
        now = UTCDateTime()
        nowst = now.strftime('%d.%m.%Y %H:%M UTC')
        out_today = t0.strftime('%Y%m%d')
        net, sta, loc, chan = stid.split('.')
        st = sds.get_waveforms(net, sta, loc, chan, 
                               starttime=t0-5*60, endtime=tf+5*60, merge=False)
        if not st:
            print(f'No waveforms for station {sta}')
            continue
        
        st.detrend('linear')
        st.filter('bandpass', freqmin=0.4, freqmax=8)
        
        vsr = stdic['amp']
        fig = plt.figure(figsize=(xdim, ydim))
        st.plot(type='dayplot', size=(xres, yres), interval=inval, 
                starttime=t0, endtime=tf, title='', fig=fig,
                subplots_adjust_top=0.95, subplots_adjust_bottom=0.05, 
                subplots_adjust_left=0.12, subplots_adjust_right=0.97, 
                dpi=dpi, tick_format='%H:%M', color=colors, show_y_UTC_label=False,
                vertical_scaling_range=vsr, linewidth=0.5)
        
        ax = fig.get_axes()[0]

        hh = int(local.utcoffset(now.datetime).seconds/3600)
        ax.set_ylabel(f'UTC (local time = UTC + {hh:02d}:00)', fontweight='bold', fontsize=8)
        ax.set_xlabel('Minutes', fontweight='bold')
        ax.xaxis.get_label().set_fontsize(ftsize)
        ax.tick_params(axis='both', which='major', labelsize=ftsize)
        ax.xaxis.set_minor_locator(MultipleLocator(40))
        ax.grid(axis='x', which='minor', linewidth=0.5, linestyle='--', alpha=0.5)
        
        ax.yaxis.get_label().set_fontsize(ftsize)
        ax.yaxis.get_label().set_fontweight('bold')
        
        ax.set_title(f'{st[0].id}, {todayst}\n' + u'Erdbebenstation Bensberg \u00A9', 
                     loc='left', fontsize=ftsize)
        ax.set_title(f'Last update: {nowst}', loc='right', fontsize=10, 
                     fontstyle='italic')
        
        if sta in ('DREG', 'HILG', 'HOBG', 'LAUG', 'RODG'):
            sta = sta[:-1]

        fname = f'{out_today}_{sta}'
        full_path = os.path.join(out_dir, f'{fname}.png')
        new_path_sta = os.path.join(out_dir, f'{sta}.png')
        fig.savefig(full_path, dpi=dpi)
        plt.close(fig)
        if not nocopy:
            os.system(f'cp {full_path} {new_path_sta}')
        
        if not noup:
            up_date = f'sshpass -f {pwdtxt} scp {full_path} {address}'
            up_sta = f'sshpass -f {pwdtxt} scp {new_path_sta} {address}'
            os.system(up_date)
            os.system(up_sta)
