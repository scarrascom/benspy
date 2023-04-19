#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
#-----------------------------------
# Filename : ppsd_sta.py
#  Purpose : Compute noise levels for different range frequencies for a given 
#            time range, specific station.
#   Author : Sebastian Carrasco (based on T. Lecocq script)
#    Email : acarrasc@uni-koeln.de
#
#-----------------------------------

import argparse
import matplotlib
import os
import sys
import seismospec
import tqdm
import warnings

import numpy as np
import pandas as pd

from glob import glob
from obspy import UTCDateTime, read, read_inventory
from obspy.clients.filesystem.sds import Client as sdsClient
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx

matplotlib.rcParams['pdf.fonttype'] = 42

def receive():
    #######################
    desc = "Create PPSDs and noise plots for custom station and date ranges"
    parser = argparse.ArgumentParser(description=desc,
                            formatter_class=argparse.RawTextHelpFormatter)
    def_inih = 0
    def_endh = 24
    def_interval = 20
    def_overlap = 0.5
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_inventory = '/mnt/Station/Inventory/LOCAL/STATIONXML/FULL_BENS.xml'
    def_step_octaves = 0.00125
    def_smooth_octaves = 0.125
    def_daytype = 'wf'
    parser.add_argument('stid', type=str, 
                        help='Full station-channel id [NET.STA.LOC.CHAN, e.g., BQ.HOBG..HHZ]')
    parser.add_argument('tini_date', type=str,
            help='Starting date in format YYYY-MM-DD [at 00:00:00]')
    parser.add_argument('tend_date', type=str,
            help='Ending date in format YYYY-MM-DD [until 23:59:59]')
    parser.add_argument('-ih', '--inihour', type=float, default=def_inih,
                        help=f'Starting hour interval\n\
[Default: {def_inih}]')
    parser.add_argument('-eh', '--endhour', type=float, default=def_endh,
                        help=f'Ending hour interval\n\
[Default: {def_endh}]')
    parser.add_argument('-i', '--interval', type=float, default=def_interval, 
                        help=f'Segment intervals in minutes\n\
[Default: {def_interval}]')
    parser.add_argument('-o', '--overlap', type=float, default=def_overlap, 
                        help=f'Overlapping factor\n\
[Default: {def_overlap}]')
    parser.add_argument('-a', '--archive', type=str, default=def_archive,
                        help=f"Directory where the waveform data is stored\n\
[Default: {def_archive}]")
    parser.add_argument('-inv', '--inventory', type=str, default=def_inventory,
                        help=f"Inventory file to be used, as STATIONXML file\n\
[Default: {def_inventory}]")
    parser.add_argument('-sm', '--smooth', type=float, default=def_smooth_octaves, 
                        help=f'Smoothing parameters - determines over what period/freq range the PSD is smoothed around T0/f0, given in fractions of octaves\n\
[Default: {def_smooth_octaves}]')
    parser.add_argument('-st', '--step', type=float, default=def_step_octaves, 
                        help=f'Step length on frequency axis in fraction of octaves\n\
[Default: {def_step_octaves}]')
    parser.add_argument('-dt', '--daytype', type=str, default=def_daytype, 
                        help=f'Type of days to calculate PPSD [wd: weekday, wn: weekend, wf: weekfull]\n\
[Default: {def_daytype}]')
    parser.add_argument('-kz', '--keepnpz', action='store_true',
                        help='Store the individual PPSDs as npz files.')
    arg = parser.parse_args()
    return arg

if __name__=="__main__":
    inputs = receive()
    stationid = inputs.stid
    t0 = inputs.tini_date
    tf = inputs.tend_date
    hini = inputs.inihour
    hend =  inputs.endhour
    mins = inputs.interval
    my_archive = inputs.archive
    inv_path = inputs.inventory
    ovlap = inputs.overlap
    smooth = inputs.smooth
    step = inputs.step
    daytype = inputs.daytype
    keep_npz = inputs.keepnpz

    network, station, location, channel = stationid.split('.')

    start = UTCDateTime(t0)
    end = UTCDateTime(tf)
    ppsd_len = mins*60

    dataset  = "Noise"
    time_zone = "Europe/Berlin"
    ## Check whether the instrument is accelerometer or not
    if channel[1]=='N':
        print('It seems to be a strong-motion sensor!')
        special_hand = 'ringlaser'
    else:
        special_hand = None

    data_provider = "BENS"
    logo = None

    # Collect the seismic waveform data
    datelist   = pd.date_range(start.datetime, end.datetime, freq="D")
    my_inventory = read_inventory(inv_path)
    c = sdsClient(sds_root=my_archive)
    
    nslc = stationid
    nslc = nslc.replace("*", "").replace("?", "")
    pbar = tqdm.tqdm(datelist)
    for day in pbar:
        datestr = day.strftime("%Y-%m-%d")
        if (daytype in ('wd', 'weekday')) and day.weekday()<5:
            pass
        elif (daytype in ('wn', 'weekend')) and day.weekday()>4:
            pass
        elif daytype in ('wf', 'weekfull'):
            pass
        else:
            print(f'Day {datestr} is not in the criteria [{daytype}]')
            continue
        fn = "{}_{}_{}.mseed".format(dataset, datestr, nslc)
        if day != UTCDateTime().datetime and os.path.isfile(fn):
            continue
        else:
            pbar.set_description("Fetching %s" % fn)
            st = c.get_waveforms(network, station, location, channel, 
                                UTCDateTime(day) + hini*3600, 
                                UTCDateTime(day) + hend*3600)
            if st:
                st.write(fn)
            else:
                pbar.set_description("No data in Archive storage for %s" % fn)
                continue

    # Compute PPSDs using custom parameters
    force_reprocess = False
    pbar = tqdm.tqdm(datelist)
    for day in pbar:
        datestr = day.strftime("%Y-%m-%d")
        if (daytype in ('wd', 'weekday')) and day.weekday()<5:
            pass
        elif (daytype in ('wn', 'weekend')) and day.weekday()>4:
            pass
        elif daytype in ('wf', 'weekfull'):
            pass
        else:
            print(f'Day {datestr} is not in the criteria [{daytype}]')
            continue
        fn_in = "{}_{}_{}.mseed".format(dataset, datestr, nslc)
        pbar.set_description("Processing %s" % fn_in)
        if not os.path.isfile(fn_in):
            continue
        stall = read(fn_in, headonly=True)
        for mseedid in list(set([tr.id for tr in stall])):
            fn_out = "{}_{}_{}.npz".format(dataset, datestr, mseedid)
            if os.path.isfile(fn_out) and not force_reprocess:
                continue
            st = read(fn_in, sourcename=mseedid)
            resp = my_inventory.select(network, station, location, channel, 
                                       starttime=UTCDateTime(day), endtime=UTCDateTime(day)+3600*24)
            if channel[1]=='N':
                resp = {'sensitivity': resp.get_response(stationid, datetime=UTCDateTime(day)).instrument_sensitivity.value }
            ppsd = PPSD(st[0].stats, metadata=resp,
                        ppsd_length=ppsd_len, overlap=ovlap,
                        special_handling=special_hand,
                        period_smoothing_width_octaves=smooth,
                        period_step_octaves=step,
                        period_limits=(0.008, 120),
                        db_bins=(-200, 20, 0.25))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ppsd.add(st)
            ppsd.save_npz(fn_out[:-4])
            del st, ppsd
        os.remove(fn_in)
        del stall

    # Reload daily PSDs from the disk and create a single PPSD object
    ppsds = {}
    pbar  = tqdm.tqdm(datelist)
    for day in pbar:
        datestr = day.strftime("%Y-%m-%d")
        fn_pattern = "{}_{}_{}.npz".format(dataset, datestr, nslc)
        pbar.set_description("Reading %s" % fn_pattern)
        for fn in glob(fn_pattern):
            mseedid = fn.replace(".npz", "").split("_")[-1]
            if mseedid not in ppsds:
                ppsds[mseedid] = PPSD.load_npz(fn)#, allow_pickle=True)
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ppsds[mseedid].add_npz(fn)#, allow_pickle=True)

    # Standard plots
    [ ppsd.plot(filename="{}_PPSD.png".format(nslc), max_percentage=10,
                        show=False, xaxis_frequency=True, cmap=pqlx)
                        for mseedid, ppsd in ppsds.items() ]
    #[ ppsd.plot_temporal(0.10) for mseedid, ppsd in ppsds.items() ]
    [ ppsd.plot_spectrogram(filename="{}_spectrogram.png".format(nslc), 
                        clim=(-160,-100), show=False)  
                        for mseedid, ppsd in ppsds.items() ]

    # Define frequency bands of interest:
    freqs = [(0.1,1.0),(0.8, 2.0),(1.0,20.0),(4.0,14.0),(4.0,20.0),(2.0,10.0)]

    displacement_RMS = {}
    for mseedid, ppsd in tqdm.tqdm(ppsds.items()):
        for per in [15, 50, 85]:
            periods, psdvalues = ppsd.get_percentile(per)
            psdvalues[psdvalues>0] = np.nan
            np.savetxt('%s_PPSD_p%s.csv' % (mseedid, per), np.transpose([periods, psdvalues]), delimiter=',', fmt='%.5e')
        ind_times = pd.DatetimeIndex([ d.datetime for d in ppsd.current_times_used ])
        data = pd.DataFrame(ppsd.psd_values, index=ind_times, 
                            columns=1./ppsd.period_bin_centers)
        data = data.sort_index(axis=1)
        displacement_RMS[mseedid] = seismospec.df_rms(data, freqs, output="DISP")
        displacement_RMS[mseedid].to_csv("%s.csv" % mseedid)

    ## Custom plot for a single frequency band
    args = {'band':"0.8-2.0",       # might be None or commented ("4.0-14.0" per default) or any of the tupples in freqs
            'time_zone':time_zone,   # required for clockplots
            'logo':logo,             # might be None or commented
            'save':'./',             # might be None or commented or a path 
            'show': False,
            'unit': 'nm'
           }

    seismospec.plot(displacement_RMS, type='timeseries', **args)
    if not keep_npz:
        print('Removing npz files...')
        os.system(f'rm *{stationid}.npz')

    # Weekday / Time of day Analysis
    #seismospec.plot(displacement_RMS, type='clockplots', **args)

    # Noise distribution over time of the day
    #seismospec.plot(displacement_RMS, type='clockmaps', **args)

    # Noise distribution over time of the day
    #seismospec.plot(displacement_RMS, type='dailyplots', **args)

    #seismospec.plot(displacement_RMS, type='gridmaps', **args)
