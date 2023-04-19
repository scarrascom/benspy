#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-

import argparse
import tqdm

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import numpy as np

from obspy import UTCDateTime
from datetime import datetime, timedelta

from obspy.clients.filesystem.sds import Client
from matplotlib.patches import Patch

def receive():
    #######################
    desc = "Create data availability figures from custom archive"
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawTextHelpFormatter)
    default_cont_bens = '/mnt/SC-Share/seiscomp/var/lib/archive'
    default_out = '/mnt/Station/data/data_avail'
    default_critdays = 180
    parser.add_argument('ini_file', type=str, 
                        help=f'File with channel codes to plot (see {default_out}/BENS.dat as example)')
    parser.add_argument('tini_date', type=str,
                        help='Starting date in format YYYY-MM-DD [e.g. 2019-03-20]')
    parser.add_argument('tend_date', type=str,
                        help='Ending date in format YYYY-MM-DD [e.g. 2019-03-21]')
    parser.add_argument('-a', '--archive', type=str, 
                        default=default_cont_bens,
                        help=f"Directory where the waveform data is stored\n\
[Default: {default_cont_bens}]")
    parser.add_argument('-d', '--dirout', type=str, 
                        default=default_out,
                        help=f"Output directory where to save the waveform files\n\
[Default: {default_out}]")
    parser.add_argument('-nt', '--notxt', action='store_true', 
                        help='Do not create txt files with availability summary.')
    parser.add_argument('-cd', '--critdays', type=int,
                        default=default_critdays,
                        help=f"Maximum number of days for which the data is read continuously\n\
[Default: {default_critdays}]")
    arg = parser.parse_args()
    return arg

"""
def run_sensor(sen_sor, days_list, _sdsa):
    net, sta, loc, cha = sen_sor.split('.')
    dict_sensor = { sen_sor: {'times': [], 'sumon': 0, 'sumoff': 0, 
                              'onper': 0, 'offper': 0}}
    nslc = f'{net}.{sta}.{loc}.{cha[:-1]}'
    daybar = tqdm.tqdm(days_list)
    for day in daybar:
        utcday = UTCDateTime(day)
        strday = utcday.strftime('%Y-%m-%d')
        daybar.set_description(strday)
        tiaux = utcday - 60*30
        tfaux = utcday + 24*3600 + 60*30
        try:	
            st = _sdsa.get_waveforms(net, sta, loc, cha, 
                                     tiaux, tfaux, merge=-1)          
            st.trim(utcday, utcday+24*3600)
            for tra in st:
                tini = tra.stats.starttime
                tend = tra.stats.endtime
                dict_sensor[sen_sor]['times'].extend([str(tini), str(tend)])
        except FileNotFoundError:	
            print("There are no waveforms for this sensor on this day")
            continue
    return(dict_sensor)
"""
def run_sensor(sen_sor, days, _sdsa, dict_in={}):
    net, sta, loc, cha = sen_sor.split('.')
    if not dict_in:
        dict_sensor = { sen_sor: {'times': [], 'sumon': 0, 'sumoff': 0,
                                  'onper': 0, 'offper': 0}}
    else:
        dict_sensor = dict_in
    nslc = f'{net}.{sta}.{loc}.{cha[:-1]}'
    iniday, endday = days
    t0aux = UTCDateTime(iniday)
    tfaux = UTCDateTime(endday)
    t0str = t0aux.strftime('%Y-%m-%d %H:%M:%S')
    tfstr = tfaux.strftime('%Y-%m-%d %H:%M:%S')
    try:
        print(f'Getting data between {t0str} and {tfstr}...')
        st = _sdsa.get_waveforms(net, sta, loc, cha, t0aux, tfaux)
        st.trim(t0aux, tfaux)
        for tra in st:
            tini = tra.stats.starttime
            tend = tra.stats.endtime
            dict_sensor[sen_sor]['times'].extend([str(tini), str(tend)])
        del(st)
    except FileNotFoundError:
        print("There are no waveforms for this sensor on this time window")
    return(dict_sensor)

def split_run_sensor(sen_sor, days, _sdsa, ninters=4):
    """
    Use this function when the time window is too long (e.g., > 1 year)
    """
    tini, tfin = days
    tot_days = (tfin - tini).days
    time_span = int(tot_days/ninters)+1
    days_list = [ tini + timedelta(days=int(x*time_span)) for x in range(ninters+1) ]
    if days_list[-1] > tfin:
        days_list[-1] = tfin

    dict_sensor = { sen_sor: {'times': [], 'sumon': 0, 'sumoff': 0,
                              'onper': 0, 'offper': 0}}
    for d, end_day in enumerate(days_list[1:]):
        ini_day = days_list[d]
        aux_sensor = run_sensor(sen_sor, (ini_day, end_day), _sdsa, dict_in=dict_sensor)
        dict_sensor = aux_sensor

    return(dict_sensor)


def get_stats(dict_sensor, tini, tfin):
    total_secs = (tfin - tini).total_seconds()
    tot_days = (tfin-tini).days
    auxdel = []
    for sensor, sensor_data in dict_sensor.items():
        print(sensor)
        times = sensor_data['times']
        if times:
            dt1 = UTCDateTime(times[0]) - UTCDateTime(t1)
            dt2 = UTCDateTime(t2) - UTCDateTime(times[-1])
            sum_on = 0
            sum_off = 0
            if dt1 > 0:
                sum_off = dt1
            for j, ttime in enumerate(times[:-1]):
                if j%2==0:
                    sum_on += (UTCDateTime(times[j+1]) - UTCDateTime(ttime))
                else:
                    sum_off += (UTCDateTime(times[j+1]) - UTCDateTime(ttime))
            if dt2>0:
                sum_off += dt2
        else:
            auxdel.append(sensor)
            sum_on = 0.0
            sum_off = total_secs
        sensor_data['sumon'] = round(sum_on, 2)
        sensor_data['sumoff'] = round(sum_off, 2)
        sensor_data['onper'] = round(sum_on/total_secs, 3)
        sensor_data['offper'] = round(sum_off/total_secs, 3)

    return(dict_sensor, tot_days)

if __name__=="__main__":
    inputs = receive()
    fname = inputs.ini_file
    tstart = inputs.tini_date
    tfin = inputs.tend_date
    archive = inputs.archive
    outdir = inputs.dirout
    notxt = inputs.notxt
    critdays = inputs.critdays
    
    sdsa = Client(archive)
    dataname = fname.split('/')[-1].split('.')[0]
    sensors = np.loadtxt(fname, dtype=str)
    
    t1 = datetime.strptime(tstart, '%Y-%m-%d')
    t2 = datetime.strptime(tfin, '%Y-%m-%d')
    tot_days = (t2-t1).days

    sensor_dict = {}

    if tot_days<=critdays:
        for sensor in sensors:
            print(sensor)
            aux_sensor_dict = run_sensor(sensor, (t1, t2), sdsa)
            sensor_dict = {**sensor_dict, **aux_sensor_dict}
    else:
        print('Splitting the time window')
        for sensor in sensors:
            print(sensor)
            aux_sensor_dict = split_run_sensor(sensor, (t1, t2), sdsa)
            sensor_dict = {**sensor_dict, **aux_sensor_dict}
    
    aux_sen_dict, total_days = get_stats(sensor_dict, t1, t2)

    print(aux_sen_dict.keys())
    sensor_dict = dict(sorted(aux_sen_dict.items(), key=lambda x: x[0]))
    str_tini = ''.join(tstart.split('-'))
    str_tend = ''.join(tfin.split('-'))
       
    ## Creating Figures
    time_on = np.asarray([ data['sumon'] for key, data in sensor_dict.items() ])
    time_off = np.asarray([ data['sumoff'] for key, data in sensor_dict.items() ])
    pert_on = np.asarray([ 100*data['onper'] for key, data in sensor_dict.items() ])
    pert_off = np.asarray([ 100*data['offper'] for key, data in sensor_dict.items() ])
    total_sensor = len(sensor_dict.keys())
    indexs = np.arange(total_sensor)
    
    ## Plot the percentage of data in the right-hand side plot
    bar_width = 0.35
    opacity   = 0.4
    
    fig = plt.figure(figsize=(15, total_sensor+1))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1 = plt.subplot(gs[1])
    
    rects1 = ax1.barh(indexs+bar_width, pert_on, bar_width, alpha=opacity, color='b', label='Data')
    
    axoff = ax1.twiny()
    rects2 = axoff.barh(indexs+bar_width, -pert_off, bar_width, 
                        alpha=opacity, color='r', label='No data')
    axoff.set_xlim(-100, 0)
    axoff.set_xticklabels([])
    
    for k, val in enumerate(pert_off):
        axoff.text(-1, indexs[k]+bar_width, f'{val:.1f}', ha='right', va='center', 
                   transform=axoff.transData)
        ax1.text(1, indexs[k]+bar_width, f'{pert_on[k]:.1f}', ha='left', va='center', 
                   transform=ax1.transData)
    
    ax1.grid(True)
    ax1.set_xlabel('Time percentage [%]')
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, total_sensor)
    ax1.set_yticks(indexs + bar_width)
    ax1.set_yticklabels([])
    
    legend_elements = [ Patch(facecolor='b', edgecolor='b', label='Data', alpha=opacity),
                        Patch(facecolor='r', edgecolor='r', label='No data', alpha=opacity)]
    
    ax1.legend(handles=legend_elements, bbox_to_anchor=(0., 1.00, 1., .102), 
               loc=3, ncol=2, mode="expand", borderaxespad=0.)
    
    ## Plot segments of data in the left-hand side plot
    ax2 = plt.subplot(gs[0])
    ax2.grid(True)
    for inde, sensor in enumerate(sensor_dict.keys()):
        sensor_data = sensor_dict[sensor]
        stimes = sensor_data['times']
        nt = len(stimes)
        print(f"This sensor has {nt} times")
        for i in range(nt):
            if i%2==0:
                t1i = UTCDateTime(stimes[i])
                t11 = datetime(t1i.year, t1i.month, t1i.day, t1i.hour, t1i.minute, t1i.second)
                t2i = UTCDateTime(stimes[i+1])
                t22 = datetime(t2i.year, t2i.month, t2i.day, t2i.hour, t2i.minute, 	t2i.second)
                ax2.hlines(y=inde+bar_width, xmin=t11, xmax=t22, 
                           linewidth=20, color='b')
    
    locator = mdates.AutoDateLocator(minticks=10, maxticks=31)
    formatter = mdates.ConciseDateFormatter(locator)
    ax2.xaxis.set_major_locator(locator)
    ax2.xaxis.set_major_formatter(formatter)    
       
    ax2.set_ylim(0, total_sensor)
    ax2.set_yticks(indexs + bar_width)
    ax2.set_yticklabels(sensor_dict.keys(), fontweight='bold')
    ax2.set_xlabel('Date')
    ax2.set_xlim(t1, t2)
    
    if total_days==1:
        fig_name = f'Dataavail_{str_tini}.{dataname}'
        tini_title = UTCDateTime(tstart).strftime('%d-%b-%Y')
        title = f'Data availability on {tini_title}'
    else:
        fig_name = f'Dataavail_{str_tini}_{str_tend}.{dataname}'
        tini_title = UTCDateTime(tstart).strftime('%d-%b-%Y')
        tend_title = UTCDateTime(tfin).strftime('%d-%b-%Y')
        title = f'Data availability between {tini_title} and {tend_title} ({total_days} days)'

    ax2.set_title(title, y=1)
    
    fig.autofmt_xdate(bottom=0.2, rotation=30, ha='right')
    fig.tight_layout()
    
    # pickle.dump(fig, open(f'{fig_name}.fig.pkl', 'wb'))
    out_fname = f'{outdir}/{fig_name}'
    sensors = np.asarray([ key for key in sensor_dict.keys() ])

    print(sensors)
    print(pert_on)
    if not notxt:
        np.savetxt(f'{out_fname}.txt', np.transpose([sensors, pert_on]), fmt='%s')
    fig.savefig(f'{out_fname}.png', dpi=400)
