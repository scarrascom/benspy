#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 10:22:02 2022

@author: seismo

Script for plotting the seismogram of a teleseismic earthquake recorded at station DREG

"""

import argparse
import sys

import numpy as np

from obspy import UTCDateTime, read_inventory
from obspy.geodetics.base import gps2dist_azimuth
from obspy.clients.filesystem.sds import Client as sdsClient
from obspy.taup import TauPyModel

import matplotlib.pyplot as plt

def receive():
    #######################
    desc = "Draw seismogram of a station for custom earthquakes"
    parser = argparse.ArgumentParser(description=desc,
                            formatter_class=argparse.RawTextHelpFormatter)
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_fmin = 0.02
    def_fmax = 0.1
    def_unit = 'DISP'
    def_length = 4
    parser.add_argument('ot', type=str,
              help='Origin time of the earthquake, in format YYYY-MM-DDTHH:MM:SS [e.g. 2019-03-20T06:15:49]')
    parser.add_argument('lat', type=float,
              help='Latitude of the hypocenter [e.g. 55.364]')
    parser.add_argument('lon', type=float,
              help='Longitude of the hypocenter [e.g. -157.888]')
    parser.add_argument('depth', type=float,
              help='Depth of the hypocenter in km [e.g. 50]')
    parser.add_argument('mag', type=float,
              help='Magnitude of the earthquake [e.g. 8.2]')
    parser.add_argument('stid', type=str,
              help='Station ID [net.sta.loc.channel]')
    parser.add_argument('ref', type=str,
              help='Reference location for the earthquake [e.g. Alaska]')
    parser.add_argument('-a', '--archive', type=str, default=def_archive,
              help=f"Directory where the waveform data is stored\n\
[Default: {def_archive}]")
    parser.add_argument('-fm', '--fmin', type=float, default=def_fmin,
              help=f"Minimum frequency to bandpass-filter the data\n\
[Default: {def_fmin} Hz]")
    parser.add_argument('-fx', '--fmax', type=float, default=def_fmax,
              help=f"Maximum frequency to bandpass-filter the data\n\
[Default: {def_fmax} Hz]")
    parser.add_argument('-u', '--units', type=str, default=def_unit,
              help=f"Output units\n\
[Default: {def_unit}]")
    parser.add_argument('-l', '--length', type=float, default=def_length,
              help=f'This number controls the duration of the seismogram, larger is longer time windows\n\
[Default: {def_length}]')
    arg = parser.parse_args()
    return arg

## Earthquake parameters
# Origin time of earthquake
if __name__=="__main__":
    inputs = receive()
    ot = inputs.ot
    lat = inputs.lat
    lon =  inputs.lon
    zkm = inputs.depth
    mag = inputs.mag
    stid = inputs.stid
    fmin = inputs.fmin
    fmax = inputs.fmax
    ref = inputs.ref
    units = inputs.units.upper()
    length = inputs.length
    archive = inputs.archive
    
    if units=='DISP':
        labelu = 'Bodenverschiebung [mm]'
    elif units=='VEL':
        labelu = 'Bodenschwinggeschwindigkeit [mm/s]'
    elif units=='ACC':
        labelu = 'Bodenbeschleunigung [mm/s/s]'
    else:
        sys.exit('Units must be DISP, VEL or ACC')
    
    t0 = UTCDateTime(ot)
    sdscl = sdsClient(archive)
    
    net, sta, loc, chan = stid.split('.')
    dl = read_inventory(f'/mnt/Station/Inventory/LOCAL/STATIONXML/{net}/{net}_{sta}.xml')
    
    lat0 = dl.get_coordinates(stid, datetime=t0)['latitude']
    lon0 = dl.get_coordinates(stid, datetime=t0)['longitude']
    z0 = dl.get_coordinates(stid, datetime=t0)['elevation']
    refname = dl.select(network=net, station=sta, location=loc, 
                        channel=chan, time=t0)[0][0].get_contents()['stations'][0]
    
    ## Get the distance between the station and hypocenter
    dm, _, _ = gps2dist_azimuth(lat0, lon0, lat, lon)
    dkm = dm/1000                       # Great-circle distance in km
    
    ## Compute the arrival times of the main waves (plist)
    ## Ray Tracing
    print('Computing the arrival times of the main waves...')
    model = TauPyModel(model='iasp91')
    Re = 6371.0                         # Radius of the Earth in kilometers
    ddeg = 360*dkm/(2*np.pi*Re)         # Distance in degrees
    print(f'Located at {dkm:.1f} km or {ddeg:.1f} degrees distance')
    plist = ['P', 'Pdiff', 'PKP', 'PKIKP', 'PP', 'Sdiff', 'S', 'SKS', 'SS']
    arrivals = model.get_travel_times(source_depth_in_km=zkm, 
                                      distance_in_degree=ddeg,
                                      phase_list=plist)
    
    print('Plotting the raypaths...')
    figr = plt.figure()
    plotarr = model.get_ray_paths(source_depth_in_km=zkm, 
                                  distance_in_degree=ddeg, phase_list=plist)
    type = 'cartesian'
    if ddeg>5:
        ptype = 'spherical'
    
    axr = plotarr.plot_rays(fig=figr, plot_type=ptype, show=False)
    
    markers = { line.get_label(): line.get_color() for line in axr.lines }
    
    ## Get the time of the earliest and latest arrivals and estimate a proper time window
    print(f'Getting the waveform data from station {stid}')
    early = arrivals[0].time
    late = arrivals[-1].time
    dt = late - early
    
    t1 = t0 + early - dt
    t2 = t0 + early + length*dt
    
    st = sdscl.get_waveforms(net, sta, loc, chan, t1, t2)
    factor = 1e3        # From meters to milimeters
    
    st.detrend('linear')
    st.taper(0.05)
    st.filter('bandpass', freqmin=fmin, freqmax=fmax, zerophase=True, corners=4)
    pre_filt = [0.9*fmin, fmin, fmax, 1.1*fmax]

    st.remove_response(inventory=dl, output='DEF', 
                       water_level=None, pre_filt=pre_filt)
    respp = dl.get_response(stid, datetime=t0)
    in_units = respp.instrument_sensitivity.input_units.upper()
    if in_units in ('M/S**2', 'M/S/S', 'M/(S**2)', 'M/SEC**2', 'M/(SEC**2)'):
        ## If input data is in acceleration units
        if units in ('VEL', 'DISP'):
            # Integrate data once
            print('Integrating once...')
            st.detrend('linear')
            st.taper(0.05)
            st.integrate()
            if units=='DISP':
                # Integrate data twice
                print('Integrating twice...')
                st.detrend('linear')
                st.taper(0.05)
                st.integrate()
    elif in_units in ('M/S', 'M/SEC'):        
        if units=='ACC':
            # Differentiate once
            print('Differentiating once...')
            st.detrend('linear')
            st.taper(0.05)
            st.differentiate()
        elif units=='DISP':
            # Integrate once
            print('Integrating once...')
            st.detrend('linear')
            st.taper(0.05)
            st.integrate()
    elif in_units=='M':
        if units in ('VEL', 'ACC'):
            print('Differentiating once...')
            st.detrend('linear')
            st.taper(0.05)
            st.differentiate()
            if units=='ACC':
                print('Differentiating twice...')
                st.detrend('linear')
                st.taper(0.05)
                st.differentiate()
    
    data = st[0].data*factor
        
    if (t2-t1) < 120:
        labelt = 's'
        fact = 1
    elif (t2-t1)>=120 and (t2-t1)<3600*2:
        labelt = 'min'
        fact = 60
    else:
        labelt = 'h'
        fact = 3600
    
    times = ((t1-t0) + st[0].times())/fact
    print('And finally plotting the seismogram...')
    fig, ax = plt.subplots(figsize=(12,6))
    ax.plot(times, data, linewidth=1.2, color='royalblue')
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5, color='gray')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel(labelu, fontweight='bold')
    ax.set_xlabel(f'Zeit nach Erdbeben [{labelt}]', fontweight='bold')

#    if t1 < t0:
    ax.set_xlim(0,)
    
    fig.subplots_adjust(top=0.981, bottom=0.073, left=0.1, 
                        right=0.988, hspace=0.2, wspace=0.2)
    
    t_text = t0.strftime('%d/%m/%Y um %H:%M:%S UTC')
    
    dict_comp = {'Z': 'vertikalen', 
                 'N': 'horizontalen (N-S)', 
                 'E': 'horizontalen (O-W)'}
    
    komp = dict_comp[chan[-1]]
    
    ax.text(0.02, 0.99, f'Seismogramm der {komp} Komponente der', 
            ha='left', va='top', transform=ax.transAxes, 
            fontsize=11, fontweight='bold')
    ax.text(0.02, 0.96, f'Station {refname}', 
            ha='left', va='top', transform=ax.transAxes, 
            fontsize=11, fontweight='bold')
    text_ort = f'Erdbeben {ref}, ~{dkm:.0f} km entfernt'
    ax.text(0.02, 0.92, text_ort, 
            ha='left', va='top', transform=ax.transAxes, 
            fontsize=11, style='italic')
    text_chars = f'Magnitude {mag}, Tiefe: {zkm} km'
    ax.text(0.02, 0.88, text_chars, 
            ha='left', va='top', transform=ax.transAxes, 
            fontsize=11, style='italic')
    ax.text(0.02, 0.84, t_text, 
            ha='left', va='top', transform=ax.transAxes, 
            fontsize=11, style='italic')
    ax.text(0.02, 0.80, f'Filter: {fmin:g}-{fmax:g} Hz', 
            ha='left', va='top', transform=ax.transAxes, 
            fontsize=11)
    
    day_eq = t0.strftime('%Y%m%d')
    fname = f'{day_eq}_{ref}_M{mag}_{sta}{chan}_{units}'
    
    fig.savefig(f'{fname}.pdf')
    fig.savefig(f'{fname}.png', dpi=300)
    figr.savefig(f'{fname}_RayPaths.png', dpi=300)
    figr.savefig(f'{fname}_RayPaths.pdf')
    
    for idx, phase in enumerate(plist):
        for arrs in arrivals:
            if arrs.name==phase:
                print(arrs.name, arrs.time)
                ax.axvline(arrs.time/fact, 0.3, 0.7, color=markers[phase], lw=2,
                           label=f'{phase}-Welle')
                # ax.text(x=arrs.time/fact, y=0.7, s=f'{phase}-Welle' , va='bottom', 
                #         transform=ax.get_xaxis_transform(), color=markers[phase], 
                #         rotation=60)
    
    ax.legend(loc='lower left', fontsize=12)
    
    fig.savefig(f'{fname}_mitWelle.png', dpi=300)
    fig.savefig(f'{fname}_mitWelle.pdf')
