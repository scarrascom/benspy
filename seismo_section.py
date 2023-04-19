#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
"""
Script for creating seismogram section, based on waveform files of each event

@author: Sebastián Carrasco
 e-mail: acarrasc@uni-koeln.de
"""

import argparse
import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import obspy.io.nordic.core as nordics

from glob import glob

from matplotlib.ticker import FuncFormatter, MultipleLocator

from obspy.io.nordic.utils import _get_line_tags
from obspy import Stream, read, read_inventory, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth

warnings.filterwarnings("ignore")

logfmt = FuncFormatter(lambda y, _: '{:g}'.format(y))

def inputs():
    desc = 'Create seismic section for specific events (from S-files with format DD*.SYYYYMM)'
    parser = argparse.ArgumentParser(description=desc, 
                            formatter_class=argparse.RawTextHelpFormatter)
    def_sdir = '/mnt/SeisAn/Seismo/REA/BNS__'
    def_wav = '/mnt/SeisAn/Seismo/WAV'
    def_invdir = '/mnt/Station/Inventory/LOCAL/STATIONXML'
    def_fmin = 2
    def_fmax = 10
    def_dmax = 70
    def_dmin = 0
    def_tini = 0
    def_tlen = 30
    def_nosta = 'NONE'
    def_nets = 'ALL'
    def_nonets = 'NONE'
    def_scale = 2
    def_outd = os.getcwd()
    def_disttype = 'epi'
    def_comp = 'Z'
    parser.add_argument('date', type=str,
                        help='Day of the event(s) to create seismic section plots [YYYY-MM-DD].')
    parser.add_argument('-sd', '--sfiledir', type=str, 
                        default=def_sdir,
                        help=f"Path of the WAV directory where the waveform data are stored\n\
[Default: {def_sdir}]")
    parser.add_argument('-wd', '--wavdir', type=str, 
                        default=def_wav,
                        help=f"Path of the WAV directory where the waveform data are stored\n\
[Default: {def_wav}]")
    parser.add_argument('-n', '--networks', type=str, 
                        default=def_nets,
                        help=f"Plot only these networks [BQN,B4,NH,LE,6E,GR,etc]\n\
[Default: {def_nets}]")
    parser.add_argument('-sc', '--scale', type=float, 
                        default=def_scale,
                        help=f"Vertical amplification of signal\n\
[Default: {def_scale}]")
    parser.add_argument('-id', '--invdir', type=str, 
                        default=def_invdir,
                        help=f"Path to the Inventories directory\n\
[Default: {def_invdir}]")
    parser.add_argument('-fn', '--fmin', type=float, 
                        default=def_fmin,
                        help=f"Minimum frequency to bandpass-filter the data\n\
[Default: {def_fmin} Hz]")
    parser.add_argument('-fx', '--fmax', type=float,
                        default=def_fmax,
                        help=f"Maximum frequency to bandpass-filter the data\n\
[Default: {def_fmax} Hz]")
    parser.add_argument('-dn', '--dmin', type=float, 
                        default=def_dmin,
                        help=f"Minimum distance from which the data is plotted\n\
[Default: {def_dmin} km]")
    parser.add_argument('-dx', '--dmax', type=float, 
                        default=def_dmax,
                        help=f"Maximum distance up to which data is plotted\n\
[Default: {def_dmax} km]")
    parser.add_argument('-ti', '--tini', type=float, 
                        default=def_tini,
                        help=f"Starting time of plotting\n\
[Default: {def_tini} s after Origin Time]")
    parser.add_argument('-tl', '--tlen', type=float, 
                        default=def_tlen,
                        help=f"Length of plotting time window\n\
[Default: {def_tlen} s]")
    parser.add_argument('-ns', '--nostations', type=str,  
                        default=def_nosta,
                        help=f"Do not plot these stations, comma separated [e.g. NAST,RODG]\n\
[Default: {def_nosta}]")
    parser.add_argument('-nn', '--nonets', type=str,  
                        default=def_nosta,
                        help=f"Do not plot these networks, comma separated [NH,GR,BQN,etc]\n\
[Default: {def_nonets}]")
    parser.add_argument('-do', '--dirout', type=str, 
                        default=def_outd, 
                        help=f"Output directory where to save the waveform files\n\
[Default: current directory = {def_outd}]")
    parser.add_argument('-dt', '--disttype', type=str, 
                        default=def_disttype, 
                        help=f"Type of distance (epi) or (hypo)central\n\
[Default: {def_disttype}]")
    parser.add_argument('-c', '--comp', type=str, 
                        default=def_comp, 
                        help=f"Component to plot (Z,N,E,R,T)\n\
[Default: {def_comp}]")
    parser.add_argument('-tx', '--texts', action='store_true',
                        help='Write station name on the right side of the plot')
    parser.add_argument('-ch', '--choose', action='store_true',
                        help='Choose specific event(s) from list (needs user input)')
    arg = parser.parse_args()
    return arg


if __name__=="__main__":
    in_puts = inputs()
    tdate = in_puts.date
    sfiles = in_puts.sfiledir
    wavdir = in_puts.wavdir
    invdir = in_puts.invdir
    fmin = in_puts.fmin
    fmax = in_puts.fmax
    dmin = in_puts.dmin
    dmax = in_puts.dmax
    tini = in_puts.tini
    tlen = in_puts.tlen
    scale = in_puts.scale
    disttype = in_puts.disttype
    texts = in_puts.texts
    choose = in_puts.choose
    comp = in_puts.comp
    nets = [ n.upper() for n in (in_puts.networks).split(',') ]
    nostatns = [ ns.upper() for ns in (in_puts.nostations).split(',') ]
    nonets = [ ns.upper() for ns in (in_puts.nonets).split(',') ]
    dirout = in_puts.dirout
    m2km = 1e3
    alfa = 0.5
    lw = 0.4
    
    dayutc = UTCDateTime(tdate)
    year = f'{dayutc.year:04}'
    month = f'{dayutc.month:02}'
    day = f'{dayutc.day:02}'
    
    gen_sfile = f'{day}-*.S{year}{month}'
    here_dir = os.path.join(sfiles, gen_sfile)
    deep_dir = os.path.join(sfiles, year, month, gen_sfile)
    events_list = np.asarray(glob(here_dir) + glob(deep_dir))
    
    ## Reading one event
    if choose:
        for s, sfile in enumerate(events_list):
            print(f' #{s+1:3d} | {sfile}')
        inputnm = input('Choose the event(s) you want to plot (comma-separated, ranges with - , e.g., 2-4,7)\n')
        auxnum = []
        for nums in inputnm.split(','):
            ranges = nums.split('-')
            if len(ranges)>1:
                auxnum.append(range(int(ranges[0])-1, int(ranges[-1]), 1))
            else:
                auxnum.append(int(ranges[0])-1)
        numevs = np.hstack(auxnum)
    else:
        numevs = np.arange(len(events_list))

    chosen_sfiles = events_list[numevs]
    for sfile in chosen_sfiles:
        print(f'# {sfile}')
        tagged_lines = _get_line_tags(open(sfile, 'r'))
        
        event = nordics.readheader(sfile)
        nordics._read_comments(tagged_lines, event)
        nordics._read_picks(tagged_lines, event)
        
        lat0 = event.origins[0].latitude
        lon0 = event.origins[0].longitude
        dep0 = event.origins[0].depth
        if lat0==None or lon0==None or dep0==None:
            print('ERROR: Missing parameters to define hypocentral location (lat, lon or depth).')
            continue
        
        t0 =  event.origins[0].time
        mag = event.magnitudes[0].mag
        magtype = event.magnitudes[0].magnitude_type
        str_t0 = t0.strftime('%d/%m/%Y %H:%M:%S')
        print(f'Latitude: {lat0:.3f}')
        print(f'Longitude: {lon0:.3f}')
        print(f'Depth: {dep0/m2km:.1f} km')
        print(f'Magnitude: {mag} {magtype}')
        print(str_t0)
        
        try:
            with open(sfile) as myFile:
                for num, line in enumerate(myFile):
                    if 'ACTION' in line:
                        nl_action = num
            for tlines in tagged_lines['3']:
                if tlines[1]==nl_action-1:
                    place = tlines[0].split('3')[0].strip()
        except IndexError:
            print('Please provide a reference town/city')
            continue
    
        ## Working with the waveforms
        wave_names = nordics.readwavename(sfile)
        full_st = Stream()
        
        for wvfile in wave_names:
            try:
                print(f'Reading {wvfile}')
                wvfile_path = os.path.join(wavdir, wvfile)
                if not os.path.isfile(wvfile_path):
                    wvfile_path = os.path.join(wavdir, 'BNS__', year, month, wvfile)
                full_st += read(wvfile_path)
            except:
                print('File is not available or not readable')

        old_stid = {'\x00\x00.ABH.\x00\x00.HH' : 'LE.ABH..HH',
                    '\x00\x00.BIW.\x00\x00.EH' : 'LE.BIW..EH',
                    '\x00\x00.FACH.\x00\x00.HH': 'LE.FACH..HH',
                    '\x00\x00.GLOK.\x00\x00.HH': 'LE.GLOK..HH',
                    '\x00\x00.LAGB.\x00\x00.HH': 'LE.LAGB..HH',
                    '\x00\x00.NICK.\x00\x00.HH': 'LE.NICK..HH',
                    '\x00\x00.OCHT.\x00\x00.EH': 'LE.OCHT..EH',
                    '\x00\x00.PYRM.\x00\x00.HH': 'LE.PYRM..HH',
                    '\x00\x00.RIVT.\x00\x00.EH': 'LE.RIVT..EH',
                    '.ACN..S ' : 'NH.ACN..EH',
                    '.BHE..S ' : 'NH.BHE..EH',
                    '.ENT..S ' : 'NH.ENT..EH',
                    '.GM1..S ' : 'NH.GM1..EH',
                    '.GM2..S ' : 'NH.GM2..EH',
                    '.GSH..S ' : 'NH.GSH..EH',
                    '.HES..S ' : 'NH.HES..EH',
                    '.JCK..S ' : 'NH.JCK..EH',
                    '.LOH..S ' : 'NH.LOH..EH',
                    '.OLF..S ' : 'NH.OLF..EH',
                    '.PLH..S ' : 'NH.PLH..EH',
                    '.RWB..S ' : 'NH.RWB..EH',
                    '.SOR..S ' : 'NH.SOR..EH',
                    '.TDN..S ' : 'NH.TDN..EH',
                    '.WBS..S ' : 'NH.WBS..EH',
                    '.XAN..S ' : 'NH.XAN..EH',
                    '.BA01..A ': 'BQ.BA01..DN',
                    '.BA02..A ': 'BQ.BA02..DN',
                    '.BA03..A ': 'BQ.BA03..DN',
                    '.BA04..A ': 'BQ.BA04..DN',
                    '.BA05..A ': 'BQ.BA05..DN',
                    '.BA06..A ': 'BQ.BA06..DN',
                    '.BA07..A ': 'BQ.BA07..DN',
                    '.BA08..A ': 'BQ.BA08..DN',
                    '.BA09..A ': 'BQ.BA09..DN',
                    '.BA10..A ': 'BQ.BA10..DN',
                    '.BA11..A ': 'BQ.BA11..DN',
                    '.BA12..A ': 'BQ.BA12..DN',
                    '.BA13..A ': 'BQ.BA13..DN',
                    '.BA14..A ': 'BQ.BA14..DN',
                    '.BA15..A ': 'BQ.BA15..DN',
                    '.BA16..A ': 'BQ.BA16..DN',
                    '.BA17..A ': 'BQ.BA17..DN',
                    '.BA18..A ': 'BQ.BA18..DN',
                    '.BA19..A ': 'BQ.BA19..DN',
                    '.BA20..A ': 'BQ.BA20..DN',
                    '.BA20..S ': 'BQ.BA20..DN',
                    '.BA21..A ': 'BQ.BA21..DN',
                    '.BA21..A L': 'BQ.BA23..DNN',
                    '.BA22..A T': 'BQ.BA23..DNE',
                    '.BA22..A ': 'BQ.BA22..DN',
                    '.BA23..A ': 'BQ.BA23..DN',
                    '.BA24..A ': 'BQ.BA24..DN',
                    '.BA26..A ': 'BQ.BA26..DN',
                    '.BA30..A ': 'BQ.BA30..DN',
                    '.BGG..S ' : 'BQ.BGG..EH',
                    '.BNS..S ' : 'BQ.BNS..EH',
                    '.BT1Q..S ': 'BQ.BT1Q..EH',
                    '.DRE..B ' : 'BQ.DREG..HH',
                    '.HIL..S ' : 'BQ.HILG..EH',
                    '.HOB..S ' : 'BQ.HOBG..EH',
                    '.JUE..S ' : 'BQ.JUE..EH',
                    '.KLL..S ' : 'BQ.KLL..EH',
                    '.KOE..S ' : 'BQ.KOE..EH',
                    '.LAU..S ' : 'BQ.LAUG..EH',
                    '.NAST..S ': 'BQ.NAST..EH',
                    '.ROD..S ' : 'BQ.RODG..EH',
                    '.STB..S ' : 'BQ.STB..EH',
                    '.SYW7..B ': 'BQ.HOBG..HH',
                    '.BD13..S ': 'B4.BD13..EH',
                    '.BD14..S ': 'B4.BD14..EH',
                    '.BD15..S ': 'B4.BD15..EH',
                    '.BD16..S ': 'B4.BD16..EH',
                    '.BD17..S ': 'B4.BD17..EH',
                    '.BHM..S ' : 'B4.BHM..EH',
                    '.BER..S ' : 'B4.BER..EH',
                    '.BOR..S ' : 'B4.BOR..EH',
                    '.FUN..S ' : 'B4.FUN..EH',
                    '.GMA..S ' : 'B4.GMA..EH',
                    '.HMB..S ' : 'B4.HMB..EH', 
                    '.HNK..S ' : 'B4.HNK..EH',
                    '.MIL..S ' : 'B4.MIL..EH',
                    '.PES..S ' : 'B4.PES..EH',
                    '.ROE..S ' : 'B4.ROE..EH',
                    '.SIN..S ' : 'B4.SIN..EH',
                    '.TGD..S ' : 'B4.TGD..EH',
                    '.AHRW..BH': 'GR.AHRW..BH',
                    '.AHRW..HH': 'GR.AHRW..HH',
                    '.BFO..BH' : 'GR.BFO..BH',
                    '.BFO..HH' : 'GR.BFO..HH',
                    '.BUG..BH' : 'GR.BUG..BH',
                    '.BUG..HH' : 'GR.BUG..HH',
                    '.TNS..HH' : 'GR.TNS..HH',
                    '.TNS..BH' : 'GR.TNS..BH',
                    '.BAVN..BH': 'RN.BAVN..BH',
                    '.BAVN..HH': 'RN.BAVN..HH',
                    'UC.VIA..HH': 'LU.VIA..HH',
                    'UC.KLB..HH': 'LU.KLB..HH',
                    }

        for tr in full_st:
            ## Checking if BA23 is on the waveform files
            if tr.id=='.BA21..A L':
                tr.id = old_stid['.BA21..A L']
            elif tr.id=='.BA22..A T':
                tr.id = old_stid['.BA22..A T']
            nslc = tr.id[:-1]
            if nslc in old_stid.keys():
                tr.id = old_stid[nslc] + tr.id[-1]
        
        print('Computing distances...')
        if comp in 'ZNE':
            sel_st = full_st.select(channel=f'??{comp}').sort(keys=['network','station', 'channel'])
            for tr in sel_st:
                net, stan, loc, chan = tr.id.split('.')
                full_inv = os.path.join(invdir, net, f'{net}_{stan}.xml')
                if net!='' and os.path.isfile(full_inv):
                    dl = read_inventory(full_inv)
                    stadl = dl.select(station=stan)
                    thislat = stadl[0][0].latitude
                    thislon = stadl[0][0].longitude
                    epidist, az, baz = gps2dist_azimuth(lat0, lon0, thislat, thislon)
                    tr.stats.baz = baz
                    if disttype.lower()=='epi':
                        tr.stats.distance = epidist/m2km
                    elif disttype.lower()=='hypo':
                        tr.stats.distance = np.sqrt(epidist**2 + dep0**2)/m2km
                    else:
                        print('ERROR: Wrong type of distance!')
                        continue
                else:
                    print(f'## Not using channel: {tr.id} ##')
                    sel_st.remove(tr)
        elif comp in 'RT':
            sel_st = full_st.sort(keys=['network','station', 'channel'])
            nslcs = []
            for tt in sel_st:
                nslcs.append(tt.id[:-1])
            uni_nslcs = np.unique(nslcs)
            for stid in uni_nslcs:
                net, stan, loc, chan = stid.split('.')
                this_sel = sel_st.select(network=net, station=stan, location=loc, channel=f'{chan}?')
                full_inv = os.path.join(invdir, net, f'{net}_{stan}.xml')
                if net!='' and os.path.isfile(full_inv):
                    dl = read_inventory(full_inv)
                    stadl = dl.select(station=stan, time=t0)
                    thislat = stadl[0][0].latitude
                    thislon = stadl[0][0].longitude
                    epidist, az, baz = gps2dist_azimuth(lat0, lon0, thislat, thislon)
                    try:
                        this_sel.rotate('->ZNE', inventory=stadl)
                        this_sel.rotate('NE->RT', baz, inventory=stadl)
                        new_sel = this_sel.select(channel=f'??{comp}')
                        if new_sel:
                            for tr in new_sel: 
                                tr.stats.baz = baz
                                if disttype.lower()=='epi':
                                    tr.stats.distance = epidist/m2km
                                elif disttype.lower()=='hypo':
                                    tr.stats.distance = np.sqrt(epidist**2 + dep0**2)/m2km
                                else:
                                    print('Wrong type of distance!')
                                    continue
                        else:
                            print(f'## Not plotting channel: {stid} (no horizontals?)##')
                    except:
                        print(f'## Cannot rotate {stid}')
                        continue
                else:
                    print(f'## Not plotting channel: {tr.id} (likely no inventory)##')
            sel_st = sel_st.select(channel=f'??{comp}')
        else:
            print('Choose one of the following components: Z, N, E, R, T')
            continue
        
        if len(sel_st)==0:
            print('There are no traces to plot for this event and criteria')
            continue

        sefo_list = ['BQ.BA01..DN', 'BQ.BA02..DN', 'BQ.BA03..DN', 'BQ.BA04..DN', 'BQ.BA05..DN',
                     'BQ.BA06..DN', 'BQ.BA07..DN', 'BQ.BA08..DN', 'BQ.BA09..DN', 'BQ.BA10..DN', 
                     'BQ.BA11..DN', 'BQ.BA12..DN', 'BQ.BA13..DN', 'BQ.BA14..DN', 'BQ.BA15..DN',
                     'BQ.BA16..DN', 'BQ.BA17..DN', 'BQ.BA18..DN', 'BQ.BA19..DN', 'BQ.BA20..DN',
                     'BQ.BA21..DN', 'BQ.BA22..DN', 'BQ.BA23..DN', 'BQ.BA24..DN', 'BQ.BA25..DN',
                     'BQ.BA26..DN', 'BQ.BA27..DN', 'BQ.BA28..DN', 'BQ.BA29..DN', 'BQ.BA30..DN']

        print('Plotting...')
        fig, ax = plt.subplots(figsize=(6, 9))
        sel_st.detrend('linear')
        sel_st.taper(0.1)
        sel_st.filter('bandpass', freqmin=fmin, freqmax=fmax, zerophase=True)
        
        data_cols = pd.read_csv('/mnt/Station/scripts/color_networks.csv', header=None)
        network_cols = { net: data_cols[1][n] for n, net in enumerate(data_cols[0])}
        
        sel_st.sort(keys=['distance'])
        for tra in sel_st:
            this_net = tra.stats.network
            this_sta = tra.stats.station
            if tra.stats.network=='BQ':
                if tra.id[:-1] in sefo_list:
                    this_net = this_net + 'N'
                else:
                    this_net = this_net + 'H'
            nonet_bool = np.logical_or(this_net in nonets, this_net[:2] in nonets)
            if tra.stats.distance <= dmax and tra.stats.distance >= dmin and\
                this_sta not in nostatns and not nonet_bool:
                label = tra.stats.network
                color = network_cols[tra.stats.network]
                if this_net in nets or nets[0]=='ALL' or this_net[:2] in nets:
                    print(f'{label}.{this_sta} - {tra.stats.distance:.1f} km')
                    sefos_bool = any( nn in nets for nn in ['ALL', 'BQ', 'BQN'])
                    if sefos_bool and tra.stats.station in sefo_list:
                        extralab = ' (N)'
                        newalpha = alfa
                    else:
                        extralab = ''
                        newalpha = alfa*2
                    ax.plot(tra.stats.starttime - t0 + tra.times(), 
                            scale*tra.data/max(abs(tra.data)) + tra.stats.distance,
                            color=color, alpha=newalpha, 
                            label=f'{label}{extralab}', linewidth=lw)
                    if texts:
                        ax.text(1.005, tra.stats.distance, f'{label}.{this_sta}', 
                                transform=ax.get_yaxis_transform(), va='center',
                                color=color, alpha=newalpha, fontsize=6)
            
        ax.set_xlim(tini, tini+tlen)
        if tlen<=30:
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_minor_locator(MultipleLocator(1))
        elif tlen > 30 and tlen <=60:
            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
        elif tlen > 60 and tlen <=200:
            ax.xaxis.set_major_locator(MultipleLocator(20))
            ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.xaxis.set_major_formatter(logfmt)
        ax.set_xlabel('Time after OT [s]')
        
        ax.set_ylim(dmin, dmax)
        ax.set_ylabel(f'{disttype.capitalize()}. distance [km]')
        
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        leg = ax.legend(by_label.values(), by_label.keys(), loc='upper left', 
                        framealpha=0.9, fontsize=8)
        
        # Set the linewidth of each legend object
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3)
        
        griddmaj = dict(which='minor', linestyle='dotted', alpha=0.3, 
                        color='gray', zorder=1)
        griddmin = dict(which='major', linestyle='dashed', alpha=0.5, 
                        color='gray', zorder=1)
        ax.grid(visible=True, **griddmaj)
        ax.grid(visible=True, **griddmin)
        ax.tick_params(direction='in', which='both')
        
        if lat0>0:
            str_lat = f'{abs(lat0):.2f}N'
        else:
            str_lat = f'{abs(lat0):.2f}S'
        
        if lon0>0:
            str_lon = f'{abs(lon0):.2f}E'
        else:
            str_lat = f'{abs(lon0):.2f}W'
        
        text_section = f"Time: {str_t0} UTC\n\
Epicenter: {str_lat}, {str_lon}\n\
Ref: {place}\n\
Depth: {dep0/m2km:.1f} km\n\
Magnitude : {mag} {magtype}"
        
        tbox = dict(boxstyle='round', ec='lightgray', fc='white', alpha=alfa)
        ax.text(0.98, 0.01, text_section, fontsize=9, ha='right', va='bottom',
                bbox=tbox, transform=ax.transAxes, ma='left', linespacing=1.5)
        
        ax.set_title(u'Erdbebenstation Bensberg, Uni Koeln \u00A9',
                     loc='left', fontsize=9)
        ax.set_title(f'Filter: {fmin:g} - {fmax:g} Hz | Comp: {comp}',
                     loc='right', fontsize=9)
        
        fig.tight_layout()
        date = t0.strftime('%Y-%m-%d_%H%M')
        
        foutname = os.path.join(dirout, f'{date}_prof')
        fig.savefig(f'{foutname}.png', dpi=200)
        #fig.savefig(f'{foutname}.pdf')