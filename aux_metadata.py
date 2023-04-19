#!/home/seismo/venvs/py3seismo/bin/python

import argparse
import base64
import hashlib
import os
import sys

from glob import glob

from obspy import read_inventory, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.io.xseed import Parser

def receive():
    #######################
    desc = "Download metadata from SMP website and save it into different formats"
    def_repo = 'bensberg-nrw'
    parser = argparse.ArgumentParser(description=desc, 
                            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r', '--repo', type=str, 
                        default=def_repo,
                        help=f'Repository to download [bensberg-nrw or bensberg-eida]\n\
[Default: {def_repo}]')
    arg = parser.parse_args()
    return arg

def hashSMPAuth(username, passphrase):
    s = hashlib.sha256()
    busername = username.encode()
    bpass = passphrase.encode()
    bserver = "smp.gempa.de".encode()
    s.update(busername)
    s.update(bserver)
    s.update(bpass)
    key = username + ':' + s.hexdigest()
    return base64.b64encode(key.encode())

if __name__=="__main__":
    inputs = receive()
    xml2dless = '/mnt/Station/scripts/stationxml-seed-converter.jar'
    repo = inputs.repo
    xmlfile = 'FULL_BENS'
    username = 'ascarrasco'
    password = 'YOsoytebo1.G'
    tokensmp = hashSMPAuth(username, password).decode('utf-8')

    ## Obspy can read up to SC3ML v0.9
    version = '0.12'
    if repo=='bensberg-nrw':
        mdata_dir = '/mnt/Station/Inventory/LOCAL'
    elif repo=='bensberg-eida':
        mdata_dir = '/mnt/Station/Inventory/EIDA'
    else:
        print('ERROR: repository must be (only) one of the following:')
        print('- bensberg-nrw')
        print('- bensberg-eida')
        sys.exit()

    root_dl = os.path.join(mdata_dir)
    wadress = f'https://smp.gempa.de/api/{username}/{repo}/download?format=SCXML/{version}'

    pathrawxml = os.path.join(root_dl, 'SCXML', f'v{version}')
    os.makedirs(pathrawxml, exist_ok=True)
    outxml = os.path.join(pathrawxml, f'{xmlfile}.scxml')
    command = f'wget -O {outxml} --header="Authorization: SMP {tokensmp}" {wadress}'
    print(command)
    os.system(command)

    print('Exporting full inventory as STATIONXML')
    inv_scml = read_inventory(outxml, format='SC3ML')
    file_path_xml = os.path.join(root_dl, 'STATIONXML')
    os.makedirs(file_path_xml, exist_ok=True)
    out_xmlfull = os.path.join(file_path_xml, f'{xmlfile}.xml')

    ## Retrieving RASPISHAKE inventory to create full inventory
    client = Client('https://data.raspberryshake.org')
    minlat = 49.9
    maxlat = 51.57
    minlon = 5.39
    maxlon = 8.11
    today = UTCDateTime('2019-01-01T00:00:00')
    print('Downloading inventory Raspberry Shakes')
    inv_rasp = client.get_stations(minlatitude=minlat, maxlatitude=maxlat, endafter=today,
                                   minlongitude=minlon, maxlongitude=maxlon, 
                                   level='response')

    inv_dl = inv_scml + inv_rasp
    print(f'Saving full inventory as STATIONXML ({xmlfile})')
    inv_dl.write(out_xmlfull, format='STATIONXML')

    print('Exporting full inventory as SEED')
    file_path_dless = os.path.join(root_dl, 'SEED')
    os.makedirs(file_path_dless, exist_ok=True)
    out_dlessfull = os.path.join(file_path_dless, f'{xmlfile}.seed')
    todless = f'java -jar {xml2dless} --input {out_xmlfull} --output {out_dlessfull}'
    os.system(todless)

    print('Exporting RESP files for all the stations')
    inv_parser = Parser(out_dlessfull)
    inv_parser.write_resp(folder=root_dl)

    print('Move all the RESP files')
    resp_files = glob(os.path.join(root_dl, 'RESP.*'))
    networks = list(dict.fromkeys([ resp.split('.')[1] for resp in resp_files ]))

    for nett in networks:
        print(f'- Network {nett}')
        path_resp = os.path.join(root_dl, 'RESP', nett)
        move_resp = os.path.join(root_dl, f'RESP.{nett}.*')
        os.makedirs(path_resp, exist_ok=True)
        comm = f'mv {move_resp} {path_resp}/.'
        os.system(comm)
        # Download the SC3ML file
        wadress = f'https://smp.gempa.de/api/{username}/{repo}/network/{nett}/download?format=SCXML/{version}'
        path_scxml = os.path.join(root_dl, 'SC3ML', f'v{version}', nett)
        os.makedirs(path_scxml, exist_ok=True)
        netxml = os.path.join(path_scxml, f'FULL_{nett}.scxml')
        command = f'wget -O {netxml} --header="Authorization: SMP {tokensmp}" {wadress}'
        os.system(command)

    if repo=='bensberg-nrw':
        print('Copying new RESP files to SeisAn CAL directory')
        cal_dir = '/mnt/SeisAn/Seismo/CAL/.'
        resp_nets = ['AM', 'BE', 'B4', 'BQ', 'LU', 'NH', 'YD', 'ZB']
        for rnet in resp_nets:
            regex_copy = os.path.join(root_dl, 'RESP', rnet, '*')
            os.system(f'cp {regex_copy} {cal_dir}')

    dict_conts = inv_dl.get_contents()
    allnets = dict_conts['networks']

    for net in allnets:
        print('##########################################')
        print(f'Saving XML and SEED files for Network {net}')
        print('##########################################')
        aux_dl = inv_dl.select(network=net)
        preffix = f'FULL_{net}'
        print(root_dl)
        file_path_xml = os.path.join(root_dl, 'STATIONXML', net)
        os.makedirs(file_path_xml, exist_ok=True)
        out_xmlfull = os.path.join(file_path_xml, f'{preffix}.xml')
        print(f'Saving {out_xmlfull}')
        aux_dl.write(out_xmlfull, format='STATIONXML')
        ## Convert from STATIONXML to SEED format
        file_path_dless = os.path.join(root_dl, 'SEED', net)
        os.makedirs(file_path_dless, exist_ok=True)
        out_dlessfull = os.path.join(file_path_dless, f'{preffix}.seed')
        todless = f'java -jar {xml2dless} --input {out_xmlfull} --output {out_dlessfull}'
        os.system(todless)
        file_path_resp = os.path.join(root_dl, 'RESP', net)
        os.makedirs(file_path_resp, exist_ok=True)
        print('Saving the stations:')
        list_sta = aux_dl.get_contents()['stations']
        stations = list(dict.fromkeys([ fullsta.split(' ')[0].split('.')[1] 
                                       for fullsta in list_sta]))
        for s, sta in enumerate(stations):
            print(f'{s+1}. {sta}')
            sta_dl = aux_dl.select(network=net, station=sta)
            preffix = f'{net}_{sta}'
            ## Save as STATIONXML file
            out_xmlsta = os.path.join(file_path_xml, f'{preffix}.xml')
            sta_dl.write(out_xmlsta, format='STATIONXML')
            ## Convert from STATIONXML to SEED format
            out_dlesssta = os.path.join(file_path_dless, f'{preffix}.seed')
            todless = f'java -jar {xml2dless} --input {out_xmlsta} --output {out_dlesssta}'
            os.system(todless)

    ### Get stations to share via EIDA, with custom starting times
    print('Saving metadata file for EIDA')
    eida_dict = {'BQ.BA01': '2019-07-21T14:00:00',
                 'BQ.BA02': '2020-03-18T12:00:00',
                 'BQ.BA03': '2019-07-21T14:00:00',
                 'BQ.BA04': '2019-07-21T14:00:00',
                 'BQ.BA05': '2020-08-05T00:00:00',
                 'BQ.BA06': '2019-07-21T14:00:00',
                 'BQ.BA08': '2019-07-21T14:00:00',
                 'BQ.BA09': '2019-11-12T07:30:00',
                 'BQ.BA10': '2019-07-21T14:00:00',
                 'BQ.BA11': '2020-07-08T00:00:00',
                 'BQ.BA12': '2023-03-06T10:15:00',
                 'BQ.BA13': '2019-12-21T00:00:00',
                 'BQ.BA16': '2019-10-28T13:00:00',
                 'BQ.BA17': '2023-02-15T00:00:00',
                 'BQ.BA18': '2019-11-01T00:00:00',
                 'BQ.BA26': '2023-02-24T10:30:00',
                 'BQ.BA30': '2020-12-01T00:00:00',
                 'BQ.BGG' : '2019-10-09T14:00:00',
                 'BQ.BNS' : '2019-08-01T00:00:00',
                 'BQ.DREG': '2016-01-31T08:56:00',
                 'BQ.HILG': '2019-08-01T00:00:00',
                 'BQ.HOBG': '2019-08-01T00:00:00',
                 'BQ.JUE' : '2019-08-01T00:00:00',
                 'BQ.KLL' : '2016-01-31T09:00:00',
                 'BQ.KOE' : '2019-08-01T00:00:00',
                 'BQ.LAUG': '2019-08-01T00:00:00',
                 'BQ.NAST': '2019-08-01T00:00:00',
                 'BQ.RODG': '2019-08-01T00:00:00',
                 'BQ.STB' : '2016-01-01T00:00:00'}

    eida_keys = list(eida_dict.keys())
    eidast = eida_keys[0]
    print('#### EIDA ####')
    print(eidast)
    net, sta = eidast.split('.')
    ini_time = UTCDateTime(eida_dict[eidast])
    eida_inv = inv_scml.select(network=net, station=sta, starttime=ini_time)
    for nets in eida_inv.networks:
        for sta in nets.stations:
            if sta.start_date < ini_time:
                sta.start_date = ini_time
            for chan in sta.channels:
                if chan.start_date < ini_time:
                    chan.start_date = ini_time

    eida_meta = '/mnt/Station/Inventory/EIDA/STATIONXML'
    for eidast in eida_keys[1:]:
        print(eidast)
        net, sta = eidast.split('.')
        sttime = eida_dict[eidast]
        ini_time = UTCDateTime(sttime)
        auxinv = inv_scml.select(network=net, station=sta, starttime=ini_time)
        eida_sta = os.path.join(eida_meta, f'{net}/{net}_{sta}.xml')
        for nets in auxinv.networks:
            if nets.start_date < ini_time:
                nets.start_date = ini_time
            for sta in nets.stations:
                if sta.start_date < ini_time:
                    sta.start_date = ini_time
                for chan in sta.channels:
                    if chan.start_date < ini_time:
                        chan.start_date = ini_time
                eida_inv.networks[0].stations.append(sta)
        auxinv.write(eida_sta, format='STATIONXML')

    eida_inv[0].start_date = UTCDateTime('2016-01-01T00:00:00')
    eida_full = os.path.join(eida_meta, 'FULL_BENS.xml')
    print(eida_full)
    eida_inv.write(eida_full, format='STATIONXML')
