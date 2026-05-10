##########################################
# utilities for the analysis of MES data
# part of this file is temporary, to be used for testing and development before integration in scpy
#
# 1. signal to noise utilities for synthetic data (snr, make_noise, add_noise)
# 2. pure spectra generation
# 3. concentration evolution simulation for various reaction mechanisms
# 4. PSD utilities
# 5. CP decomposition
# 6. specific readers of publicated data
# 7. animated gif plotter
#
# %%
import numpy as np
import traitlets as tr
from numpy.random import RandomState
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import re
from datetime import datetime, timedelta
import os
import tensorly as tl

from typing import Union, Sequence

import spectrochempy as scp
from spectrochempy import NDDataset, Coord, CoordSet
from spectrochempy.analysis._base._analysisbase import AnalysisConfigurable

num = Union[int, float]



def read_raman_raw(filename):
    times = []
    data = []

    with open(filename, 'r') as file:
        # return file.readline().strip()
        ramanshifts  = ([float(x) for x in file.readline().strip().split()])

        while True:
            line = file.readline()
            if line != '':
                items = line.strip().split()
                times.append(float(items[0]))
                data.append([float(x) for x in items[1:]])
            else:
                break

    data = np.array(data)
    times = np.array(times) * 24 * 3600
    out = NDDataset(data)
    out.name = filename.split('/')[-1]
    out.name = os.path.basename(filename)
    out.x = Coord(ramanshifts, title='Raman Shift', units='1/cm')
    out.y = Coord(times, title='timestamp', units='s')
    out.title = 'Raman intensity'
    out.units = 'a.u.'

    return out



def read_xas(dir_path):

    dates = []
    timestamps = []
    data = []
    energies = []

    # read files in the directory
    for filename in os.listdir(dir_path):
        if filename.endswith('.txt'):

            data_ = []
            meta_ = {}

            with open(dir_path+'/'+filename, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('#'):
                        key_value = re.split(':|=', line[1:], 1  )
                        if len(key_value) == 2:
                            meta_[key_value[0].strip()] = key_value[1].strip()
                    else:
                        data_.append([float(x) for x in line.split()])

                hdf_file = meta_['Data from HDF File']
                date_start = datetime.strptime(meta_['Time at start'], '%Y-%m-%d %H:%M:%S.%f')

                date = date_start + timedelta(seconds=float(meta_['Time from start (seconds)']))
                dates.append(date)
                timestamps.append(date.timestamp())
                data_ = np.array(data_)
                data.append(data_[:,1])
                energies.append(data_[:,0])

    # check that all energies are the same
    energies = np.array(energies)
    if np.count_nonzero(energies - energies[0]):
        raise ValueError('Inconsistent energy values in the files')


    data = np.array(data)

    out = NDDataset(data)
    out.name = dir_path.split('/')[-2]
    out.title = 'absorbance'
    out.units = 'absorbance'
    out.x = Coord(energies[0], title='energy', units='eV')
    out.y = Coord(timestamps, title='timestamp', units='s')

    return out

