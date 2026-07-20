import sys
import os
import numpy as np

#individual instances for this class for each coords, target, edge_attr, features

class Statistics:
    def __init__(self): 
        self.count = 0
        self.mean = 0
        self.M2 = 0
    
    def WelfordAlgorithm(self, col):
        
        for curr_val in col:
            self.count+=1
            delta1 = curr_val - self.mean

            self.mean = self.mean + (curr_val - self.mean)/self.count

            delta2 = curr_val - self.mean

            self.M2 = self.M2 + delta1*delta2

    def CalculatePopulationVariance(self):
        if self.count == 0:
            return 0.0
        return self.M2/self.count

    def CalculatePopulationStd(self):
        return np.sqrt(self.CalculatePopulationVariance())

