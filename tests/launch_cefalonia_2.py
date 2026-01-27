#!/usr/bin/env python

"""
To call LOKI you need this file and make it executable (or call python)

"""

from loki.loki import Loki

db_path = '/home/emanuele/data/emanuele/loki-das/Traveltimes/cefalonia'
data_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/orion_das_high_freq_complete_weak_less_subsection_also_non_selected_2_gilbert'
output_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/output_new_hybrid_3'
hdr_filename = 'header_long_cefalonia.hdr'
geometry_filename_fiber = 'channels.dat'
geometry_filename_stat = 'stations_cefalonia_total.dat'
inputs = {}
inputs['tshortp_min_sta'] = 0.1 #0.1
inputs['tshortp_max_sta'] = 0.2 #0.2
inputs['tshorts_min_sta'] = 0.1 #0.15
inputs['tshorts_max_sta'] = 0.2 #0.25
inputs['tshortp_min_fiber'] = 0.1 #0.1
inputs['tshortp_max_fiber'] = 0.2 #0.3
inputs['tshorts_min_fiber'] = 0.1 #0.15
inputs['tshorts_max_fiber'] = 0.2 #0.35
inputs['slrat'] = 2
inputs['npr'] = 15
inputs['ntrial'] = 4
inputs['derivative'] = True
inputs['normalize'] = True
inputs['vfunc'] = 'erg'
inputs['hfunc'] = 'pca'
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
