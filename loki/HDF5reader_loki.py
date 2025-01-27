# -*- coding: utf-8 -*-
"""
HDF5reader v1.0 - Halliburton HDF5 raw file reader for Python. 

Copyright (C) 2024 Halliburton

Author: A. Ellmauthaler
Last modified: 03/19/2024
"""

import h5py
import os
import datetime
import numpy as np
import math



""" 
Get RAW data from HDF5 files

Input:
******
dirName:       Parent directory holding HDF5 files 
channRange:    Channel indices to load (slice or None)
               None => load all available channels (default)
sampleRange:   Sample indices to read. When reading multiple files, the index counter continuous into subsequent files
               None => load all available time samples (default)
               int   => Load this number of samples, starting from first sample.
               slice => Range of indices, i.e. samples=slice(start,stop)  
integrate:     Integrate along time axis of phase data (default)
               If False, the output will be the time differential of the phase (delta phase)

Output:
*******
data:          Numpy data matrix
                Rows: Channels
                Columns: Time
metaData:      Recording Meta data (adjusted to reflect extracted channel and sample range)            
"""
def getData(dirName, dataFileList, channRange=None, sampleRange=None, integrate=True):

    # Get list of HDF5 files
    #dataFileList = []
    #for root, directories, files in os.walk(dirName):
    #        for file in files:
    #            if file.endswith(".h5"):
    #                dataFileList.append(os.path.join(root, file).replace("\\","/"))
    #dataFileList.sort()
    #print('\n\nFound {0} HDF5 files\n\n'.format(len(dataFileList)))

    # Retrieve meta data
    metaData = getMetaData(dirName)

    # Extract user-defined channel range and adjust recording meta data information accordingly
    if channRange is None:
        channRange = slice(0,metaData['numChanns'])
    elif channRange.stop <= metaData['numChanns']:
        metaData['numChanns'] = channRange.stop-channRange.start
        metaData['channelMDs'] = metaData['channelMDs'][channRange]
    else:
        print('Invalid Channel Range\n\n')
        return(None, None)

    ## Extract user-defined sample range and adjust recording meta data information accordingly
    #if sampleRange is None:
    #    sampleRange = slice(0,metaData['framesPerFile']*len(dataFileList))
    #elif isinstance(sampleRange, int):
    #    sampleRange = slice(0, sampleRange)
        
    #if sampleRange.stop <= metaData['recDuration']*metaData['fs_out']:
    #    metaData['timeStamps'] = metaData['timeStamps'][sampleRange]
    #    metaData['recDuration'] =(sampleRange.stop-sampleRange.start)/metaData['fs_out']
    #    metaData['recStartTime'] += datetime.timedelta(microseconds=(sampleRange.start/metaData['fs_out'])*10**6)
    #else:
    #    print('Invalid Sample Range\n\n')
    #    return(None, None)
    #
    ## Get data
    #print('Extracting data ... \n')
    #dataFileList = dataFileList[(sampleRange.start // metaData['framesPerFile']):(math.ceil(sampleRange.stop / metaData['framesPerFile']))]
    
    #print(dataFileList)
    sampleRange=sampleRange = slice(0,metaData['framesPerFile']*len(dataFileList))   # 19/11/2024 
    metaData['timeStamps'] = metaData['timeStamps'][sampleRange]
    metaData['recDuration'] =(sampleRange.stop-sampleRange.start)/metaData['fs_out']
    metaData['recStartTime']
    
    data = np.zeros((channRange.stop-channRange.start, len(dataFileList)*metaData['framesPerFile']))
    n = 0
    for fileName in dataFileList:
        with h5py.File(fileName, 'r') as f:
            data[:, n:n+metaData['framesPerFile']] = f['Acquisition']['Raw[0]']['RawData'][channRange, :]
            n += metaData['framesPerFile']
            print(fileName)
    # Extract sample interval of interest
    data = data[:,(sampleRange.start-(sampleRange.start // metaData['framesPerFile'])*metaData['framesPerFile']):(sampleRange.stop-(sampleRange.start // metaData['framesPerFile'])*metaData['framesPerFile'])]
    # Normalize the data
    data = data/metaData['normFact']

    # If integrate is True, convert delta phase to relative phase
    if integrate is True:
        data = np.cumsum(data, axis=1)

    print('Done\n\n')

    return(data, metaData)
                  

# Get critical meta data from HDF5 file or directory holding HDF5 files
def getMetaData(name):

    if os.path.isdir(name):
        # Get list of HDF5 files
        dataFileList = []
        for root, directories, files in os.walk(name):
            for file in files:
                if file.endswith(".h5"):
                    dataFileList.append(os.path.join(root, file).replace("\\","/"))
        file = h5py.File(dataFileList[0], 'r')
    else:
        file = h5py.File(name, 'r')
         
    # Initialize dictionary holding metadata information
    metaData = {}
    # Populate dictionary with metadata
    metaData['fs_out'] = file['Acquisition']['Raw[0]'].attrs['OutputDataRate'] # Output Sampling Frequency
    metaData['fs_orig'] = file['Acquisition'].attrs['PulseRate'] # Original Sampling Frequency
    metaData['gaugeLength'] = file['Acquisition'].attrs['GaugeLength'] # Gauge Length
    metaData['pulseWidth'] = file['Acquisition'].attrs['PulseWidth'] # Pulse Width
    metaData['dx_orig'] = file['Acquisition']['Custom'].attrs['OriginalSpatialSamplingInterval'] # Original Spatial Sampling Interval
    metaData['dx_out'] = file['Acquisition'].attrs['SpatialSamplingInterval'] # Output Spatial Sampling Interval
    metaData['interrogatorSerial'] = file['Acquisition']['Custom'].attrs['InterrogatorSerial'].tostring().decode('UTF-8') # DAS Interrogator Serial Number
    metaData['ituChannels'] = file['Acquisition']['Custom'].attrs['ITUChannels'] # DAS Interrogator Serial Number
    metaData['lengthUnit'] = file['Acquisition'].attrs['SpatialSamplingInterval.uom'].tostring().decode('UTF-8') # Length Unit (meter/feet)
    metaData['dataType'] = 'Differential Phase' # Data Type
    metaData['normFact'] = 10430 # Data Normalization Factor
    metaData['numChanns'] = file['Acquisition'].attrs['NumberOfLoci'] # Channel Number
    metaData['framesPerFile'] = int(file['Acquisition']['Raw[0]']['RawData'].attrs['Count'] / metaData['numChanns']) # Number of frames per file
    metaData['channelMDs'] = file['Acquisition']['FacilityCalibration[0]']['Calibration[1]']['LocusDepthPoint']['FacilityLength'] # Measured Depth of each channel
    metaData['recStartTime'] = datetime.datetime.strptime(file['Acquisition'].attrs['MeasurementStartTime'].tostring().decode('UTF-8'), '%Y-%m-%dT%H:%M:%S.%f%z') # Recording Start Time
    # Recording duration
    if os.path.isdir(name):
        metaData['recDuration'] = int(len(dataFileList) * metaData['framesPerFile']/metaData['fs_out'])
    else:
        metaData['recDuration'] = int(metaData['framesPerFile']/metaData['fs_out'])
    # Timestamp for each sample in file
    metaData['timeStamps'] = np.arange(0,metaData['recDuration'], 1/metaData['fs_out'])

    # Close HDF5 file
    file.close()

    return(metaData)


# Print critical meta data
def printMetaData(metaData):
    print('HDF5 Meta Data:')
    print('***************')
    print('DAS Interrogator: ' + str(metaData['interrogatorSerial']))
    print('ITU Channels: ' + str(metaData['ituChannels']))
    print('Recording Start Time: ' + str(metaData['recStartTime']))
    print('Recording Duration: {0}h {1}m {2}s'.format(metaData['recDuration']//3600, (metaData['recDuration'] % 3600) // 60, metaData['recDuration'] % 60))
    print('Sampling Frequency - Original: ' + str(metaData['fs_orig']) + 'Hz')
    print('Sampling Frequency - Output: ' + str(metaData['fs_out']) + 'Hz')
    print('Spatial Sampling Interval - Original: {0:.2f}{1}'.format(metaData['dx_orig'], metaData['lengthUnit'])) 
    print('Spatial Sampling Interval - Output: {0:.2f}{1}'.format(metaData['dx_out'],metaData['lengthUnit']))
    print('Gauge Length: ' + str(metaData['gaugeLength']) + 'm')
    print('Pulse Width ' + str(metaData['pulseWidth']) + 'ns')
    print('Unit of Length: ' + str(metaData['lengthUnit']))
    print('Data Type: ' + str(metaData['dataType']))
    print('Frames per File: ' + str(metaData['framesPerFile']))
    print('Number of Channels: ' + str(metaData['numChanns']))
    print('Depth Range: {0:.2f}{1} - {2:.2f}{1}'.format(metaData['channelMDs'][0], metaData['lengthUnit'], metaData['channelMDs'][-1]))
    print('\n\n')


# Function that dumps all meta data contained in HDF5 file or directory holding HDF5 files to user-defined text file
def dumpHDF5HeadertoFile(name, outFileName):

    if os.path.isdir(name):
        # Get list of HDF5 files
        dataFileList = []
        for root, directories, files in os.walk(name):
            for file in files:
                if file.endswith(".h5"):
                    dataFileList.append(os.path.join(root, file).replace("\\","/"))
        file = h5py.File(dataFileList[0], 'r')
    else:
            file = h5py.File(name, 'r')

    # Get groups from HDF5 file
    groups = []
    file.visit(groups.append)

    with open(outFileName, "w") as outFile:
        # Get attributed from HDF5 file
        for group in groups:
            dset = file[group]
            outFile.write(str(group) + '\n')
            for attrs in dset.attrs:
                outFile.write(str(attrs) + ': ' + str(dset.attrs[attrs]) + '\n')
            outFile.write('\n\n')
