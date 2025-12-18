#!/usr/bin/env python

"""
To call LOKI you need this file and make it executable (or call python)

"""

from loki.loki import Loki

db_path = '/home/emanuele/data/emanuele/loki-das/Traveltimes/bedretto'
data_path = '/home/emanuele/data/emanuele/loki-das/bedretto/m_new/hybrid'
output_path = '/home/emanuele/data/emanuele/loki-das/bedretto/output_n_new_prova'
hdr_filename = 'header_long_bedretto.hdr'
geometry_filename_fiber = 'DAS_geom.txt'
geometry_filename_stat = 'hybrid_geom.dat'
inputs = {}
inputs['tshortp_min_sta'] = 0.004 #0.1
inputs['tshortp_max_sta'] = 0.004#0.2
inputs['tshorts_min_sta'] = 0.002 #0.15
inputs['tshorts_max_sta'] = 0.004 #0.25
inputs['tshortp_min_fiber'] = 0.002 #0.1
inputs['tshortp_max_fiber'] = 0.004 #0.3
inputs['tshorts_min_fiber'] = 0.002 #0.15
inputs['tshorts_max_fiber'] = 0.004 #0.35
inputs['slrat'] = 2
inputs['npr'] = 36
inputs['ntrial'] = 5
inputs['derivative'] = True
inputs['normalize'] = True
inputs['vfunc'] = 'erg' #'tkeo'
inputs['hfunc'] = 'pca' #tkeo
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
