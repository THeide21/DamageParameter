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
#abaqus python Schaediungsparameter.py 210000 0.5 1000 03_Tension_Static OneElement
#abaqus python Schaediungsparameter.py 210000 0.5 1000 05_Validierung Pobe Set-2
#abaqus python Schaediungsparameter.py 204900 0.5 1000 07_Validierung_Coarse_Mat Pobe
#abaqus python Schaediungsparameter.py 204900 0.3 400 LCF-Strain2 Step-Static Gleich
    rawArgList = sys.argv
    interactiveFlag = 0
    Para = {}
    E = float(rawArgList[1])
    k =  float(rawArgList[2])
    S_yield = float(rawArgList[3])
    odb_name = rawArgList[4]
    part_name = rawArgList[5]
    StepName = rawArgList[6]

    print('Commandline-Mode:')
    print('--------------------')
    print('E: [MPa]',E)
    print('Zugfestigkeit [MPa] ',S_yield)
    print('k [-]:',k)
    print('Odb-Name:',odb_name)
    print('Part-Name:',part_name)
    print('Step-Name:',StepName)
    print('interactiveFlag:',interactiveFlag)
    print('--------------------')
    if len(rawArgList) == 8:
        setName = rawArgList[7]
        makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_name,StepName,interactiveFlag,setName.upper())
        
    else:
        makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_name,StepName,interactiveFlag)    
    
