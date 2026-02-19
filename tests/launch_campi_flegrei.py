#!/usr/bin/env python

"""
To call LOKI you need this file and make it executable (or call python)

"""

from loki.loki import Loki

db_path = '/home/emanuele/data/emanuele/loki-das/Traveltimes/campi-flegrei'
data_path = '/home/emanuele/data/emanuele/loki-das/campi-flegrei/test_vlp_12_30'
output_path = '/home/emanuele/data/emanuele/loki-das/campi-flegrei/output_vlp_12_30_no_sta_lta'
hdr_filename = 'header_long_campi_flegrei.hdr'
geometry_filename_fiber = 'DAS_geom.txt'
geometry_filename_stat = 'DAS_geom.txt'
inputs = {}
#inputs['tshortp_min_sta'] = 1 #0.1
#inputs['tshortp_max_sta'] = 1 #0.2
#inputs['tshorts_min_sta'] = 1 #0.15
#inputs['tshorts_max_sta'] = 1 #0.25
#inputs['tshortp_min_fiber'] = 1 #0.1
#inputs['tshortp_max_fiber'] = 1 #0.3
#inputs['tshorts_min_fiber'] = 1 #0.15
#inputs['tshorts_max_fiber'] = 1 #0.35
#inputs['slrat'] = 2
inputs['npr'] = 36
inputs['ntrial'] = 1
inputs['derivative'] = True
inputs['normalize'] = True
inputs['vfunc'] = 'tkeo'
inputs['hfunc'] = 'tkeo'
inputs['model'] = 'layer'
inputs['epsilon'] = 0.001
#inputs['freq']=[2,15]
precision = 'single'
comp = ['E', 'N', 'Z']
inputs['extension_sta'] = '*'
inputs['extension_das'] = 'CANDAS2_2023-01-07_10-48-10.h5'
inputs['delta_das'] = 200

# =========  Call


l1 = Loki(data_path, output_path, db_path, hdr_filename, geometry_filename_fiber, geometry_filename_stat, mode='locator')
l1.location( comp, precision, **inputs)
