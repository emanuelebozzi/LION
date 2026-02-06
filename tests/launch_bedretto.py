#!/usr/bin/env python

"""
To call LOKI you need this file and make it executable (or call python)

"""

from loki.loki import Loki

db_path = '/home/emanuele/data/emanuele/loki-das/Traveltimes/bedretto'
#data_path = '/home/emanuele/data/emanuele/loki-das/bedretto/Bedretto_DAS/M0B_snippets_mseed_denoised'
#data_path = '/home/emanuele/data/emanuele/loki-das/bedretto/test_hybrid_network_3_02_26'
data_path = '/home/emanuele/data/emanuele/loki-das/bedretto/processed_events_denoised'
output_path = '/home/emanuele/data/emanuele/loki-das/bedretto/output_first_50_events_magnitude_geophones_only_denoised'
hdr_filename = 'header_long_bedretto.hdr'
geometry_filename_fiber = 'DAS_geom.txt'
geometry_filename_stat = 'hybrid_geom_geophones_only.dat'
inputs = {}
inputs['tshortp_min_sta'] = 0.005 #0.1
inputs['tshortp_max_sta'] = 0.02 #0.2
inputs['tshorts_min_sta'] = 0.005 #0.15
inputs['tshorts_max_sta'] = 0.02 #0.25
inputs['tshortp_min_fiber'] = 0.005 #0.1
inputs['tshortp_max_fiber'] = 0.02 #0.3
inputs['tshorts_min_fiber'] = 0.005 #0.15
inputs['tshorts_max_fiber'] = 0.02 #0.35
inputs['slrat'] = 2
inputs['npr'] = 30
inputs['ntrial'] = 4
inputs['derivative'] = True
inputs['normalize'] = True
inputs['vfunc'] = 'erg' #'tkeo'
inputs['hfunc'] = 'tkeo' #tkeo
inputs['model'] = 'layer'
inputs['epsilon'] = 0.001
#inputs['freq']=[2,15]
precision = 'single'
comp = ['E', 'N', 'Z']
inputs['extension_sta'] = '*'
inputs['extension_das'] = 'CANDAS2_2023-01-07_10-48-10.h5'
inputs['delta_das'] = 1000

# =========  Call


l1 = Loki(data_path, output_path, db_path, hdr_filename, geometry_filename_fiber, geometry_filename_stat, mode='locator')
l1.location( comp, precision, **inputs)
