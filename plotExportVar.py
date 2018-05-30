# -*- coding: utf-8 -*-
"""
Created on Wed May 30 23:16:57 2018

@author: Thomas Heidebrecht
"""

import numpy as np 
import matplotlib.pyplot as plt

def plotExportVar(path):
    var = np.loadtxt(path)
    plt.plot(var)
    
    
if __name__ == '__main__':
    plotExportVar('../EXPORT.txt')