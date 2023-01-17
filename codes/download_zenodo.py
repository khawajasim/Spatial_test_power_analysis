#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 13:42:37 2023

@author: khawaja
"""

import requests
import sys
import json
import os
import hashlib


def check_hash(filename, checksum):
    algorithm, value = checksum.split(':')
    if not os.path.exists(filename):
        return value, 'invalid'
    h = hashlib.new(algorithm)
    with open(filename, 'rb') as f:
        while True:
            data = f.read(4096)
            if not data:
                break
            h.update(data)
    digest = h.hexdigest()
    return value, digest


def download_file(url, filename):
    progress_bar_length = 72
    block_size = 1024
    r = requests.get(url, stream=True)
    try:
        total_size = int(r.headers.get('content-length'))
    except TypeError:
        total_size = None
    download_size = 0
    if total_size:
        print(f'Downloading file with size of {total_size / block_size:.3f} kB')
    else:
        print(f'Downloading file with unknown size')
    with open(filename, 'wb') as f:
        for data in r.iter_content(chunk_size=block_size):
            download_size += len(data)
            f.write(data)
            # for progress bar
            if total_size:
                progress = int(progress_bar_length*download_size/total_size)
                sys.stdout.write('\r[{}{}] {:.1f}%'.format('█'*progress, '.' * (progress_bar_length-progress),
                    100*download_size/total_size))
                sys.stdout.flush()
        sys.stdout.write('\n')
        

link = "https://zenodo.org/api/records/7153351"    
filename =  '../data'

r = requests.get(link)

download_url = [f['links']['self'] for f in r.json()['files']]
# filenames = [f['key'] for f in r.json()['files']]
filenames = [(f['key'], f['checksum']) for f in r.json()['files']]


for (fname, checksum), url in zip(filenames, download_url):
    print(url)
    print(fname)  
    print(checksum)
    # os.makedirs(dir_map[fname])
    folder_for_download =  os.path.join('../', fname)
    #full_path = os.path.join(dir_map[fname], fname)
    download_file(url, folder_for_download)
    
    
    value, digest = check_hash(folder_for_download, checksum)
    if value != digest:
       print("Error: Checksum does not match.")
       sys.exit(-1)
    
    # r = requests.get(url)
