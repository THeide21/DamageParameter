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
# Für Konsolen aufruf:
#abaqus python Schaediungsparameter.py 204900 0 450 oe_elastic_Tension Part-1 Step-1
#abaqus python Schaediungsparameter.py 204900 0.5 400 LCF-Strain2_jan Probe Step-Static Set-Mitte2
#abaqus python Schaediungsparameter.py 204900 0.5 2270 LCF-Strain2  Probe Step-Static DEBUG
#abaqus python Schaediungsparameter.py 204900 0.5 2270 McClafin_Test_Setup-A1_Mitte3 Tube Step-load Set-Mitte3
#abaqus python Schaediungsparameter.py 20450.000000 0.500000 450.000000 oe_8_Elements PolyElement Step-1 -0.113000 -610.000000 1311.000000 0.875500 0.500000
#abaqus python Schaediungsparameter.py 204900 0.0 2270 McClafin_Test_Setup-A1_Mitte3 Tube Step-load -0.113000 -610.000000 1311.000000 0.875500 0.500000 0 10 all Set-Mitte3
#abaqus python Schaediungsparameter.py 204900 0.5 2270 McClafin_Test_Setup-D1_Mitte1 Tube Step-load -0.113000 -610.000000 1311.000000 0.875500 0.500000 0 10 all Set-Mitte1
    rawArgList = sys.argv
    interactiveFlag = 1
    Para = {}
    E = float(rawArgList[1])
    k =  float(rawArgList[2])
    S_yield = float(rawArgList[3])
    odb_name = rawArgList[4]
    part_name = rawArgList[5]
    StepName = rawArgList[6]
    b = float(rawArgList[7])
    c = float(rawArgList[8])
    sigma_F = float(rawArgList[9])
    epsilion_F = float(rawArgList[10])
    portion = float(rawArgList[11])
    t_begin = float(rawArgList[12])
    t_end = float(rawArgList[13])
    frame_flag = rawArgList[14]
    print('Commandline-Mode:')
    print('--------------------')
    print('E: [MPa]',E)
    print('Zugfestigkeit [MPa] ',S_yield)
    print('k [-]:',k)
    print('Odb-Name:',odb_name)
    print('Part-Name:',part_name)
    print('Step-Name:',StepName)
    print('interactiveFlag:',interactiveFlag)
    print('Parameter für Lebensdauer:')
    print('b:',b)
    print('c:',c)
    print('sigma_F',sigma_F)
    print('epsilion_F:',epsilion_F)
    print('portion',portion)
    print('t_begin',t_begin)
    print('t_end',t_end)
    print('frame_flag',frame_flag)
    print('--------------------')
    if len(rawArgList) == 16:
        setName = rawArgList[15]
        # odb_name,E,k,S_yield,part_Name,StepName,t_begin,t_end,b,c,sigma_F,epsilion_F,portion,interactiveFlag=0,setName = None,frame_flag='all'
        makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_name,StepName,t_begin,t_end,b,c,sigma_F,epsilion_F,portion,interactiveFlag,setName.upper(),frame_flag)
        
    else:
        setName = None
        makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_name,StepName,t_begin,t_end,b,c,sigma_F,epsilion_F,portion,interactiveFlag,setName,frame_flag)    
    
