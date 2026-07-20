import sys
import os

import numpy as np

class Normalizer:

    #stats_tuple is over the entire training dataset
    def __init__(self, tensor_tuple, stats_tuple):
        #np arrays
        coords, target, edge_attr, features = tensor_tuple
        #scalars
        coords_mean, target_mean, target_std, euclid_mean, euclid_std, degree_mean, degree_std = stats_tuple

        self.CenterCoordinates(coords, coords_mean)
        self.Standardize(target, target_mean, target_std, 0)
        self.Standardize(edge_attr, euclid_mean, euclid_std, 4)
        self.Standardize(features, degree_mean, degree_std, 14)

    #out_mean = np.empty((1, 3), dtype=np.float32) 
    def CenterCoordinates(self, coords, coords_mean):
        coords -= coords_mean
    
#standardized over entire dataset
#applied to: target, euclidean edge distance, atom degree
    def Standardize(self, data, mean, std, col):

        is_row = data.ndim == 1

        if (is_row): 
            data[col] -= mean
            if (std == 0): std += 1e-8
            data[col] /= std
        
        else:
            data[:,col] -= mean
            if (std == 0): std += 1e-8

            data[:,col] /= std

        

        

    




