# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import time, math, warnings, sys, os
import numpy as np
import itertools as it
try:
    from textRepr import * 
    import visualization
    from abaqusConstants import *
    import abaqus
    print('Richtig bei Gui')
except:  
    from odbAccess import *
    print('Richtig bei Command')
from Modul_Schaedigung import *


debug = 0


    
if __name__ == '__main__' and debug == 0:
# Wert reinfolge für aufruf aus Konsole
# E-Modul, k, Zugfestigkeit, odb-Name, Set-Name
# Für Konsolen aufruf
#abaqus python Schaediungsparameter.py 210000 0.5 10000 03_Tension_Static OneElement
#abaqus python Schaediungsparameter.py 210000 0.5 10000 05_Validierung Pobe Set-2
    
    rawArgList = sys.argv
    interactiveFlag = 0
    Para = {}
    E = float(rawArgList[1])
    k =  float(rawArgList[2])
    S_yield = float(rawArgList[3])
    odb_name = rawArgList[4]
    part_name = rawArgList[5]
    print('Commandline-Mode:')
    print('--------------------')
    print('E: [MPa]',E)
    print('Zugfestigkeit [MPa] ',S_yield)
    print('k [-]:',k)
    print('Odb-Name:',odb_name)
    print('Part-Name:',part_name)
    print('interactiveFlag:',interactiveFlag)
    print('--------------------')
    if len(rawArgList) == 7:
        print('I am here')
        setName = rawArgList[6]
        makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_name,interactiveFlag,setName.upper())
        
    else:
        makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_name,interactiveFlag)    
    
