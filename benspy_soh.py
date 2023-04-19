#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-

import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from obspy import UTCDateTime
from datetime import datetime, timedelta

from obspy.clients.filesystem.sds import Client

from benspy_getwaves import crea_dictsta

def receive():
    #######################
    desc = "Create State of Health (soh) plots for stations transmitted via Seedlink"
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawTextHelpFormatter)
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_outd = '/mnt/Station/data/soh'
    def_soh = 'VDT,GNS,VEI,LCQ'
    parser.add_argument('ini_file', type=str, 
                        help='File with channel codes to plot')
    parser.add_argument('tini_date', type=str,
                        help='Starting date in format YYYY-MM-DD [e.g. 2019-03-20]')
    parser.add_argument('tend_date', type=str,
                        help='Ending date in format YYYY-MM-DD [e.g. 2019-03-21]')
    parser.add_argument('-a', '--archive', type=str, 
                        default=def_archive,
                        help=f"Directory where the waveform data is stored [SC-Share]\n\
[Default: {def_archive}]")
    parser.add_argument('-d', '--dirout', type=str, 
                        default=def_outd,
                        help=f"Output directory where to save the waveform files\n\
[Default: {def_outd}]")
    parser.add_argument('-s', '--soh', type=str,
                        default=def_soh,
                        help=f'SOH channels to plot, comma separated\n\
[Default: {def_soh}]')
    arg = parser.parse_args()
    return arg

def run_sensor(sen_sor, days_list):
    net, sta, loc, cha = sen_sor.split('.')
    ns = f'{net}.{sta}'
    dict_sensor = { sen_sor: { 'times': [], 'data': []}}
    for day in days_list:
        print(day)
        utcday = UTCDateTime(day)
        tiaux = utcday - 60*30
        tfaux = utcday + 24*3600 + 60*30
        try:	
            thissds = sdsa
            st = thissds.get_waveforms(net, sta, loc, cha, 
                                       tiaux, tfaux, merge=-1)
            if st:
                st.trim(utcday, utcday+24*3600)
                for tra in st:
                    dict_sensor[sen_sor]['times'].append(
                        np.asarray([ str(tra.stats.starttime + dt) for dt in tra.times() ]))
                    dict_sensor[sen_sor]['data'].append(tra.data)
        except FileNotFoundError:	
            print("There are no waveforms for this sensor at this day")
            continue
    ## Flatten the times and data lists:
    if dict_sensor[sen_sor]['times']:
        dict_sensor[sen_sor]['times'] = np.hstack(dict_sensor[sen_sor]['times'])
        dict_sensor[sen_sor]['data'] = np.hstack(dict_sensor[sen_sor]['data'])
    return(dict_sensor)

if __name__=="__main__":
    inputs = receive()
    fname = inputs.ini_file
    tstart = inputs.tini_date
    tfin = inputs.tend_date
    archive = inputs.archive
    outdir = inputs.dirout
    soh_chs = inputs.soh
    
    dpis = 200
    sdsa = Client(archive)
    
    dataname = fname.split('/')[-1].split('.')[0]
    sensors1 = np.loadtxt(fname, dtype=str)
    
    t1 = datetime.strptime(tstart, '%Y-%m-%d')
    t2 = datetime.strptime(tfin, '%Y-%m-%d')
    list_days = [ t1 + timedelta(days=x) for x in range(0, (t2-t1).days)]
    
    soh_chans = soh_chs.split(',')
    soh_dict = {
                'GAN': {'label': 'GNSS antenna status',  'units': 'N',      'amp': 1 }, 
                'GEL': {'label': 'Elevation',            'units': 'm',      'amp': 1e6 },
                'GLO': {'label': 'Longitude',            'units': '$^o$',   'amp': 1e6 }, 
                'GLA': {'label': 'Latitude',             'units': '$^o$',   'amp': 1e6 },
                'GNS': {'label': 'GNSS satellites',      'units': 'N',      'amp': 1 }, 
                'GPL': {'label': 'GNSS PLL status',      'units': 'N',      'amp': 1 },
                'GST': {'label': 'GNSS status',          'units': 'N',      'amp': 1 },
                'LCE': {'label': 'Clock phase error',    'units': 's',      'amp': 1e6},
                'LCQ': {'label': 'Clock quality',        'units': '$\%$',   'amp': 1},
                'VCO': {'label': 'Control voltage',      'units': 'cts',    'amp': 1},
                'VDT': {'label': 'Digitizer temperature','units': '$^o$C',  'amp': 1e3},
                'VEC': {'label': 'Digitizer current',    'units': 'A',      'amp': 1e3},
                'VEI': {'label': 'Input voltage',        'units': 'V',      'amp': 1e3}
                }

    ## Check first if the SOH channels are available
    for soh_ch in soh_chans:
        if soh_ch not in soh_dict.keys():
            sys.exit(f'Channel {soh_ch} is not in the pre-defined list')

    for soh_ch in soh_chans:
        print('#############################')
        print(f'Plotting channel {soh_ch}')
        soh_amp = soh_dict[soh_ch]['amp']
        soh_label = soh_dict[soh_ch]['label']
        soh_units = soh_dict[soh_ch]['units']
        soh_title = f'{soh_ch} - {soh_label} [{soh_units}]'
        sensor_dict = {}
        sensors = [ f'{sen}.D0.{soh_ch}' for sen in sensors1 ]
        for sensor in sensors:        
            auxd_sensor = run_sensor(sensor, list_days)
            sensor_dict = {**sensor_dict, **auxd_sensor}
        
        total_days = (t2-t1).days    
        print(sensor_dict.keys())
        sensor_dict = dict(sorted(sensor_dict.items(), key=lambda x: x[0]))
        str_tini = ''.join(tstart.split('-'))
        str_tend = ''.join(tfin.split('-'))
        
        ## Creating Figures
        total_sensor = len(sensor_dict.keys())
        indexs = np.arange(total_sensor)   
        color_rec = '#f1f1f1'
        
        ## TODO: One figure for each station, if required
        ## One Figure for each SOH channel
        fig, axs = plt.subplots(figsize=(0.9*total_sensor, 1.2*total_sensor), 
                                nrows=total_sensor, sharex=True, sharey=True)
        gridmaj = dict(which='major', linestyle='--', alpha=0.8, color='gray', zorder=2)    
        gridmin = dict(which='minor', axis='y', linestyle='--', alpha=0.5, color='gray', zorder=2)
        
        for inde, sensor in enumerate(sensor_dict.keys()):
            sensor_data = sensor_dict[sensor]
            stimes = np.asarray([pd.to_datetime(date) for date in sensor_data['times']])
            thisax = axs[inde]
            coldots = 'orange'
            thisax.set_facecolor(color_rec)
            if (inde%2==1):
                coldots = 'royalblue'
            thisax.scatter(stimes, np.asarray(sensor_data['data'])/soh_amp, 
                        s=40, c=coldots, alpha=0.5, zorder=2)
            thisax.grid(**gridmaj)
            thisax.grid(**gridmin)
            ylab = '.'.join(sensor.split('.')[:2])
            thisax.set_ylabel(ylab, fontweight='bold', rotation=0, ha='right')
        
        if total_days<10:
            locator = mdates.AutoDateLocator(minticks=10, maxticks=10)
        else:
            locator = mdates.AutoDateLocator(minticks=10, maxticks=30)
        formatter = mdates.ConciseDateFormatter(locator)
        ax0 = axs[0]
        ax0.xaxis.set_major_locator(locator)
        ax0.xaxis.set_major_formatter(formatter)
             
        ax0.set_xlabel('Date')
        ax0.set_xlim(t1, t2)
                  
        if total_days==1:
            fig_name = f'SOH_{str_tini}.{soh_ch}.{dataname}'
            tini_title = UTCDateTime(tstart).strftime('%d-%b-%Y')
        else:
            fig_name = f'SOH_{str_tini}_{str_tend}.{soh_ch}.{dataname}'
        
        ax0.set_title(soh_title, y=1)
        fig.autofmt_xdate(bottom=0.2, rotation=30, ha='right')
        fig.tight_layout()
        
        out_fname = f'{outdir}/{fig_name}'
        sensors = np.asarray([ key for key in sensor_dict.keys() ])    
        print(sensors)
        fig.savefig(f'{out_fname}.png', dpi=dpis)
