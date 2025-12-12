#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
@author: john-marshall@northwestern.edu
## helper functions used in loading output of the MALTAB CNMFE files
# create_fluorescence_time_delta
# normalize
# return_spatial_info
# create_contour_layouts
# com

adapted from 
"Helper functions to support 'run_cnmfe_matlab.py'
  Please see that file for more details.
  Created on Nov 1 2017
  @author: tamachado@stanford.edu"

"""

from past.utils import old_div
from matplotlib import pyplot as plt
from scipy import stats
import scipy.sparse as sparse
import numpy as np
import miniscope_analysis as ma
from tqdm import tqdm
import scipy.io as sio
import pandas as pd
import scipy.spatial.distance as dist
import itertools
import statsmodels.formula.api as smf
import math
import dlc_utils

## helper functions used in loading output of the MALTAB CNMFE files
# create_fluorescence_time_delta
# normalize
# return_spatial_info
# create_contour_layouts
# com


def normalize(trace, percentile=True):
    """ Normalize a fluorescence trace by its max or its 99th percentile. """
    trace = trace - np.min(trace)
    if np.percentile(trace, 99) > 0:
        if percentile:
            trace = trace / np.percentile(trace, 99)
        else:
            trace = trace / np.max(trace)
    return trace


def com(A, d1, d2):
    """Calculation of the center of mass for spatial components

       From Caiman: https://github.com/flatironinstitute/CaImAn
       @author: agiovann

     Inputs:
     ------
     A:   np.ndarray
          matrix of spatial components (d x K)

     d1:  int
          number of pixels in x-direction

     d2:  int
          number of pixels in y-direction

     Output:
     -------
     cm:  np.ndarray
          center of mass for spatial components (K x 2)
    """
    nr = np.shape(A)[-1]
    Coor = dict()
    Coor['x'] = np.kron(np.ones((d2, 1)), np.expand_dims(list(range(d1)),
                        axis=1))
    Coor['y'] = np.kron(np.expand_dims(list(range(d2)), axis=1),
                        np.ones((d1, 1)))
    cm = np.zeros((nr, 2))        # vector for center of mass
    cm[:, 0] = old_div(np.dot(Coor['x'].T, A), A.sum(axis=0))
    cm[:, 1] = old_div(np.dot(Coor['y'].T, A), A.sum(axis=0))

    return cm

def return_spatial_info(path_to_cnmfe, spatial_threshold, dims=(752, 480)):
  CNMFE_results = sio.loadmat(path_to_cnmfe)
  spatial_components=np.array(CNMFE_results['A'].todense())
  d1 = dims[0]
  d2 = dims[1]
  coms = com(spatial_components, d1, d2)
  com_df = pd.DataFrame(coms, columns=['y', 'x'], index=[int(index) for index in np.linspace(1, len(coms), len(coms))])
  return(com_df, spatial_components)



def create_contour_layouts(spatial_components, dims=(752, 480)):
  # return dict with info for plotting
  x, y = np.mgrid[0:dims[0]:1, 0:dims[1]:1]
  cell_contours = {}
  for_dims = {}
  to_plot = (0, np.shape(spatial_components)[1])
  for i in range(to_plot[0], to_plot[1]):
    Bvec = spatial_components[:, i].flatten()
    #normalize contours to 1
    Bvec /= np.max(Bvec)
    thr = 0.6
  # rehape to dimensions of image
    Bmat = np.reshape(Bvec, (dims), order='F')
    cell_contours[i+1] = Bmat
    for_dims[i+1] = Bvec
  return(cell_contours, for_dims)

def create_fluorescence_time_delta(fluoresence_data):
    fluorescence = pd.DataFrame(np.transpose(z_score_CNMFE(fluoresence_data)),
      columns=[int(cell_num) for cell_num in np.linspace(1, len(fluoresence_data), len(fluoresence_data))])
    fluorescence['msCamFrame'] = fluorescence.index.values
    fluorescence = fluorescence.set_index(pd.to_timedelta(np.linspace(0, (len(fluorescence)-1*(1/20)), len(fluorescence)), unit='s'), drop=False)
    return(fluorescence)
