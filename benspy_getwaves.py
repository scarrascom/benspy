#!/home/seismo/venvs/py3seismo/bin/python
# -*- coding: utf-8 -*-
#-----------------------------------
# Filename : bens_getwvfm.py
#  Purpose : Retrieve custom data from BENS (continuous or event) archive and
#            save it into waveform files (mseed, sac, ascii).
#   Author : Sebastian Carrasco
#    Email : acarrasc@uni-koeln.de
#
#-----------------------------------

import argparse
import json
import os
import re
import sys
import tqdm

from obspy import UTCDateTime, Stream
from obspy.clients.filesystem.sds import Client as sdsClient

def check_inputs(inputs):
    tini_date = inputs.tini_date
    tini_hour = inputs.tini_hour
    tend_date = inputs.tend_date
    tend_hour = inputs.tend_hour
    comps = (inputs.components).split(',')
    modeout = inputs.modeout
    formatout = inputs.formatout
    re_date = re.compile('.*-.*-.')
    re_hour = re.compile('.*:.*:.')
    re_hors = re.compile('.*:.*')
    
    ## Checking given initial and end time
    if not re_date.match(tini_date) or not re_date.match(tend_date):
        sys.exit('Date must be in YY-MM-DD format')
        
    if not re_hour.match(tini_hour) and not re_hors.match(tini_hour):
        sys.exit('Starting time must be in HH:MM[:SS] format')
        
    if not re_hour.match(tend_hour) and not re_hors.match(tend_hour):
        sys.exit('Ending time must be in HH:MM[:SS] format')
    
    ## Checking output mode
    if modeout not in ('full', 'sta', 'com', 'comsta'):
        sys.exit('Output mode must be one of the following options: \
                             full, sta, cha, chasta')
    
    used_comps = [ c.upper() in 'ZNE12' for c in comps ]
    if not all(used_comps):
        sys.exit('Components must be Z, N and/or E.')
    
    ## Checking output format
    if formatout.upper()=='MSEED':
        dict_out = {'encoding': 'STEIM2', 'record_length': 512, 
                    'dataquality': 'D',     'byteorder': '>' }
    elif formatout.upper()=='SAC':
        dict_out = { 'byteorder': '>' }
    elif formatout.upper()=='SLIST':
        dict_out = {}
    else:
        sys.exit('The desired output format is not available, choose among SAC, MSEED or SLIST')
    ## Return the dictionary with the kwargs for output
    return dict_out

def receive():
    #######################
    name_script = __file__
    desc = f"Retrieve waveform data for a set of stations belonging to the BENS network for a given time range.\n\
Examples:\n\
1. Retrieve waveform data from all strong- and weak-motion stations in one-single waveform file:\n\
   > {name_script} 2020-05-01 14:30:43 2020-05-01 14:44\n\
2. Retrieve waveform data from stations BGG, HOBG, DREG and BA01 and save data into per-component files:\n\
   > {name_script} 2019-12-10 23:54 2019-12-11 00:10:35 -s BGG,HOBG,DREG,BA01 -m com\n\
3. Retrieve waveform data from all Sefonib stations except BA09 and provide one miniseed file per station:\n\
   > {name_script} 2020-03-20 15:55 2020-03-20 16:30:25 -g sefo -ns BA09 -m sta\n\
4. Retrieve waveform data using the discrete per-event archive:\n\
   > {name_script} 2011-07-19 05:56 2011-07-19 06:01 -a /mnt/Projects/EVTarchive -ed /mnt/Station/Inventory/json/all_bns_stations_evt.json -s BGG"
    parser = argparse.ArgumentParser(description=desc, 
                            formatter_class=argparse.RawTextHelpFormatter)
    def_archive = '/mnt/SC-Share/seiscomp/var/lib/archive'
    def_json = '/mnt/Station/Inventory/json/all_bns_stations.json'
    def_sta = 'ALL'
    def_group = 'ALL'
    def_nosta = 'NONE'
    def_comps = 'Z,N,E,1,2'
    def_mode = 'full'
    def_outd = os.getcwd()
    def_format =  'MSEED'
    parser.add_argument('tini_date', type=str, 
                        help='Date of starting time in format YYYY-MM-DD [e.g. 2019-03-20]')
    parser.add_argument('tini_hour', type=str, 
                        help='Hour of starting time in format HH:MM[:SS] [e.g. 20:30:54]')
    parser.add_argument('tend_date', type=str, 
                        help='Date of ending time in format YYYY-MM-DD [e.g. 2019-03-21]')
    parser.add_argument('tend_hour', type=str, 
                        help='Hour of starting time in format HH:MM[:SS] [e.g. 01:02:53]')
    parser.add_argument('-a', '--archive', type=str, 
                        default=def_archive, 
                        help=f"Directory where the waveform data is stored\n\
[Default: {def_archive}]")
    parser.add_argument('-s', '--stations', type=str, 
                        default=def_sta,
                        help=f"Stations to retrieve data, comma separated\n\
[Default: {def_sta}]")
    parser.add_argument('-g', '--group', type=str, 
                        default=def_group,
                        help=f"Group of stations to retrieve [BENS, SEFO, RHB, GD, LGB, EIFEL]\n\
[Default: {def_group}]")
    parser.add_argument('-ns', '--nostations', type=str,  
                        default=def_nosta,
                        help=f"Do not store these stations, comma separated [e.g. NAST,ROD]\n\
[Default: {def_nosta}]")
    parser.add_argument('-c', '--components', type=str, 
                        default=def_comps,
                        help=f'Components to retrieve, comma separated [e.g. Z,N,E]\n\
[Default: {def_comps}]')
    help_modeout = f"Output mode of waveforms. Choose among:\n\
    -full: full (One single file with all the data)\n\
    -com: component (One file per component including all the stations)\n\
    -sta: station (One file per station including all the components)\n\
    -comsta: component and station (One file per component and per station)\n\
[Default: {def_mode}]"
    parser.add_argument('-m', '--modeout', type=str, 
                        default='full',
                        help=help_modeout)
    parser.add_argument('-d', '--dirout', type=str, 
                        default=def_outd, 
                        help=f"Output directory where to save the waveform files\n\
[Default: current directory = {def_outd}]")
    parser.add_argument('-f', '--formatout', type=str, 
                        default=def_format, 
                        help=f"Output format of the waveform data among [MSEED, SAC or ASCII]\n\
[Default: {def_format}]")
    parser.add_argument('-mg', '--merge', action='store_true',
                        help='Merge gappy data by filling the gaps (latest value)')
    arg = parser.parse_args()
    return arg

def crea_dictsta(_stal, _groups, _nostal):
    sefo_group = ['BA01', 'BA02', 'BA03', 'BA04', 'BA05', 'BA06', 
                  'BA08', 'BA09', 'BA10', 'BA11', 'BA12', 'BA13', 
                  'BA14', 'BA15', 'BA16', 'BA17', 'BA18', 'BA19', 
                  'BA20', 'BA21', 'BA22', 'BA23', 'BA26', 'BA30']
    bens_group = ['BNS', 'BGG', 'DREG', 'HILG', 'HOBG', 'JUE', 'KLL',  
                  'KOE', 'LAUG', 'NAST', 'RODG', 'STB', 'SYW7', 
                  'BT2A', 'BT2B', 'BT2C', 'BT2D', 'BT2E', 'BA02H']
    rhb_group = ['BD14', 'BD15', 'BD17', 'BER', 'FUN', 
                 'HMB', 'HNK', 'MIL', 'ROE', 'SIN']
    gd_group = ['ABH', 'ACN', 'BHE', 'ENT', 'FACH', 
                'GMA', 'GM1', 'GM2', 'GSH', 'GWBC', 
                'HES', 'JCK', 'LOH', 'OLF', 'PLH', 
                'RWB', 'SOR', 'TDN', 'WBS', 'XAN']
    lgb_group = ['ABH', 'BEDO', 'BEUR', 'BIW', 'FACH', 'GLOK', 'KOGO', 
                 'MUEZ', 'NICK', 'OMED', 'PYRM', 'RIVT']
    eifel_group = ['TA000', 'TA001', 'TA002', 'TA003', 'TA004', 'TA005', 'TA006', 'TA007', 'TA008', 'TA009', 
                   'TA011', 'TA012', 'TA013', 'TA014', 'TA018', 'TA019', 'TA020', 'TA022', 'TA024', 'TA025', 
                   'TA026', 'TA027', 'TA028', 'TA029', 'TA030', 'TA031', 'TA032', 'TA033', 'TA034', 'TA035', 
                   'TA036', 'TA037', 'TA038', 'TA039', 'TA040', 'TA041', 'TA042', 'TA043', 'TA044', 'TA045', 
                   'TA046', 'TA047', 'TA048', 'TA049', 'TA050', 'TA052', 'TA053', 'TA054', 'TA056', 'TA057', 
                   'TA058', 'TA059', 'TA060', 'TA061', 'TA062', 'TA063', 'TA064', 'TA065', 'TA066', 'TA067', 
                   'TA068', 'TA070', 'TA072', 'TA073', 'TA074', 'TA075', 'TA076', 'TA077', 'TA087', 'TA088', 
                   'TA089', 'TA090', 'TA097', 'TA098', 'TA113', 'TA114', 'TA120', 'TA125', 'TA152', 'TA166', 
                   'TA167', 'TA168', 'TA169', 'TA170', 'TA171', 'TA174', 'TA176', 'TA177', 'TA226', 'TA247', 
                   'TA253', 'TA254', 'TA255', 'TA256', 'TA257', 'TA259', 'TA260', 'TA261', 'TA284', 'TA289', 
                   'TA290', 'TA302', 'TA317', 'TA335', 'TA341', 'TA410', 'TB015', 'TB016', 'TB017', 'TB078', 
                   'TB079', 'TB080', 'TB081', 'TB082', 'TB083', 'TB085', 'TB092', 'TB093', 'TB094', 'TB095', 
                   'TB096', 'TB099', 'TB100', 'TB103', 'TB104', 'TB105', 'TB108', 'TB110', 'TB112', 'TB115', 
                   'TB117', 'TB118', 'TB119', 'TB123', 'TB124', 'TB126', 'TB127', 'TB145', 'TB148', 'TB151', 
                   'TB153', 'TB154', 'TB157', 'TB158', 'TB163', 'TB172', 'TB175', 'TB208', 'TB209', 'TB211', 
                   'TB232', 'TB236', 'TB243', 'TB251', 'TB262', 'TB263', 'TB264', 'TB265', 'TB268', 'TB282', 
                   'TB288', 'TB318', 'TB337', 'TB356', 'TB357', 'TB386', 'TB390', 'TB412', 'TC101', 'TC102', 
                   'TC111', 'TC140', 'TC141', 'TC150', 'TC160', 'TC161', 'TC181', 'TC213', 'TC219', 'TC233', 
                   'TC238', 'TC269', 'TC274', 'TC287', 'TC297', 'TC303', 'TC312', 'TC332', 'TC339', 'TC387', 
                   'TD137', 'TD178', 'TD188', 'TD191', 'TD194', 'TD201', 'TD205', 'TD215', 'TD216', 'TD235', 
                   'TD237', 'TD273', 'TD291', 'TD292', 'TD368', 'TE131', 'TE134', 'TE138', 'TE182', 'TE183', 
                   'TE207', 'TE231', 'TE252', 'TE278', 'TE281', 'TE342', 'TE380', 'TF129', 'TF180', 'TF298', 
                   'TF338', 'TG128', 'TG130', 'TG192', 'TG223', 'TG280', 'TG314', 'TG329', 'TH132', 'TH189', 
                   'TH315', 'TI221', 'TJ240', 'TX347', 'TX348', 'TX349', 'TX350', 'TX351', 'TX352', 'TX353']
    dflt = 'ALL'
    dict_groups = {'SEFO': sefo_group, 
                   'BENS': bens_group,
                   'RHB': rhb_group,
                   'GD': gd_group,
                   'LGB': lgb_group,
                   'EIFEL': eifel_group,
                    dflt: sefo_group + bens_group + rhb_group + gd_group + lgb_group + eifel_group}
    new_dict = []
    ## TODO: IMPROVE THIS FUNCTION
    if dflt in _groups:
        if dflt in _stal:
            new_dict = new_dict + dict_groups[dflt]
        else:
            new_dict = new_dict + _stal
    else:
        for group in _groups:
            if dflt in _stal:
                new_dict = new_dict + dict_groups[group]
            else:
                for sta in dict_groups[group]:
                    if sta in _stal:
                        new_dict.append(sta)
    
    for nostation in _nostal:
        if nostation in new_dict:
            new_dict.remove(nostation)

    print(new_dict)        
    sta_list = list(dict.fromkeys(new_dict))
    nets = []
    for key, val in dict_groups.items():
        if any([sta in val for sta in sta_list]):
            nets.append(key)
    
    nets.remove('ALL')
    if len(nets)==5:
        nets = ['ALL']
    
    return sta_list, nets

def out_put(t0, st, formout, dictout, nets, modeout='full', outdir='.'):
    found_stations = list(dict.fromkeys([ tr.stats.station for tr in st ]))
    nchannels = list(dict.fromkeys([ tr.id  for tr in st ]))
    ncha = len(nchannels)
    if not nchannels:
        sys.exit('There are no waveforms for the indicated component, \
                         station and time interval.')
    else:
        suffix_out = f'________{ncha:03}'

    preffix_out = t0.strftime('%Y-%m-%d-%H%M-%SS')
    print("Writing output file")

    if modeout=='full':
        net_join = ''.join(nets)
        name_file = f'{preffix_out}.{net_join}{suffix_out[len(net_join):]}'
        outfile = os.path.join(outdir, name_file)
        st.write(outfile, format=formout.upper(), **dictout)
        print(outfile)
        return outfile
    elif modeout=='com':
        net_join = ''.join(nets)
        for comp in 'ZNEXY12':
            name_file = f'{preffix_out}{comp}.{net_join}{suffix_out[len(net_join):]}'
            outfile = os.path.join(outdir, name_file)
            sel_st = st.select(channel='??'+comp).copy()
            if len(sel_st)>0:
                sel_st.write(outfile, format=formout.upper(), **dictout)
                print(outfile)
    elif modeout=='sta':
        for sta in found_stations:
            name_file = f'{preffix_out}.{sta}.{formout.lower()}'
            outfile = os.path.join(outdir, name_file)
            sel_st = st.select(station=sta)
            if len(sel_st)>0:
                sel_st.write(outfile, format=formout.upper(), **dictout)
                print(outfile)
    elif modeout=='cs':
        for sta in stations:
            for comp in 'ZNEXY12':
                name_file = f'{preffix_out}.{sta}.{comp}.{formout.lower()}'
                outfile = os.path.join(outdir, name_file)
                sel_st = st.select(station=sta, channel='??'+comp).copy()
                if len(sel_st)>0:
                    sel_st.write(outfile, format=formout.upper(), **dictout)
                    print(outfile)
    else:
        sys.exit("Something went wrong, the desired output mode is likely incorrect.")

def _get_data(t0, tf, acq_dir, dict_station, components, 
              prio_chans=['HH', 'EH', 'BH', 'SH', 'DN', 'HN', 'DP'],  ## The first channel with the data available will be retrieved
              merge=False):
    ## Set the client variables for retrieving the data
    cthis = sdsClient(sds_root=acq_dir)
    st_aux = Stream()
    sta_tqdm = tqdm.tqdm(dict_station)
    for station in sta_tqdm:
        sta_tqdm.set_description(station)
        network = '??'
        location = '*'
        prech0 = ''
        statn = station.upper()
        if statn=='BA09' and t0 > UTCDateTime('2022-02-14T07:50:00'):
            network = 'BQ'
            statn = 'BNS'
            prech0 = 'DN'
        elif statn=='BA02H':
            network = 'BQ'
            statn = 'BA02'
            prech0 = 'EH'
        elif statn=='BA02':
            network = 'BQ'
            statn = 'BA02'
            prech0 = 'DN'
        for prech in prio_chans:
            if prech0:
                prech = prech0
            channel = f'{prech}[ZNE12]'
            aux_st = cthis.get_waveforms(network, statn, location, channel, t0, tf, merge=-1)
            st_aux += aux_st
            ## If there are data available for this channel, then don't include any more channels (first in the list are priorities)
            if len(aux_st)>0:
                break
    if merge:
        ## We just decide to use the latest sample value if yes (so Seisan can handle the waveform data and no jumps are observed)
        st = st_aux.merge(method=0, fill_value='latest')
    else:
        st_aux.merge(method=0)
        st = st_aux.split()
    return st

if __name__=="__main__":
    ## Receiving inputs
    inputs = receive()
    acq_dir = inputs.archive
    tini_date = inputs.tini_date
    tini_hour = inputs.tini_hour
    tend_date = inputs.tend_date
    tend_hour = inputs.tend_hour
    formatout = inputs.formatout
    modeout = inputs.modeout
    dirout = inputs.dirout
    merge = inputs.merge
    stations = [ s.upper() for s in (inputs.stations).split(',') ]
    networks = [ g.upper() for g in (inputs.group).split(',') ]
    nostatns = [ ns.upper() for ns in (inputs.nostations).split(',') ]
    components = [ c.upper() for c in (inputs.components).split(',') ]
    
    ## Check the inputs and retrieve the dictionary 
    dict_out = check_inputs(inputs)
    dict_station, nets = crea_dictsta(stations, networks, nostatns)
    ## Get the dates and search for the files
    starttime = UTCDateTime("%sT%s" % (tini_date, tini_hour))
    endtime = UTCDateTime("%sT%s" % (tend_date, tend_hour))

    st = _get_data(starttime, endtime, acq_dir, 
                   dict_station, components, merge=merge)
    ## Output data
    out_put(starttime, st, formatout, dict_out, nets, modeout, dirout)
