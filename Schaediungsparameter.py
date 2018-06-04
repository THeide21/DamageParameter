# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#+-------------------------------+
#|   block of user's constants   |
#+-------------------------------+
#odbPathName = 'C:\Users\ThomasHeidebrecht\Documents\Abaqus\Benchmark_Coarse.odb'
#oldStepName = 'Step-1'


#+-------------------------------+
#|   TO-DO:   |
#+-------------------------------+
# 1. Function witch creat Set of 
#   Elements by pretending Set-Name. If no set is defined: default= all Elements)
#2. GUI
#
#
#
#
#
#
#+----------------------+
#|   block of modules   |
#+----------------------+
#	

import visualization
from abaqusConstants import *
import math
import abaqus
from odbAccess import *
import time
import numpy as np
import warnings

def tic():
    # Homemade version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    if 'startTime_for_tictoc' in globals():
        print( "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")
        
        
#
#
#+------------------------+
#|   block of functions   |
#+------------------------+

def OPENodb(name_ODB,odbPathName): 
    'Opens ODB-File as readable'
    odb =session.openOdb(name=name_ODB, path=odbPathName, readOnly=False)     
    return odb

#+------------------------+

def getFrames(step_name,odb):
    'return Frames from a Step'
    frames = odb.steps[step_name].frames
    return frames
#+------------------------+

def getSteps(odb):
    'return Frames from a Step'
    steps = odb.steps
    return steps

#+------------------------+
 
def getEIDS(instance,SetName = None): 
    if SetName == None:
        elements = instance.elements
    else:
        elements = instance.elementSets[SetName]
    eIDS=map(lambda element:element.label,elements)
    return eIDS

#+------------------------+
 
def getGlobalValues(frame,flag):
    'return of the whole Stressfield'
    var = frame.fieldOutputs[flag]
    return var

#+------------------------+
 
def getValueONelement(frame,flag,eID):
    'return of the values on ELement'
    var = frame.fieldOutputs[flag].values[eID].data
    return var

#+------------------------+
 
def VectorToTensor(Vec,flag):
    'Transform Stress or Strain Vector in Voigt-Notation to numpy-Tensor'
    warnings.simplefilter('error', UserWarning)
    Tensor = np.zeros([3, 3])
    Tensor[0,0] = Vec[0] 
    Tensor[1,1] = Vec[1]
    Tensor[2,2] = Vec[2]
    if flag == 'S':
        Tensor[0,1] =Vec[3]#1,2
        Tensor[0,2] =Vec[4]#1,3
        Tensor[1,0] =Vec[3]#2,1
        Tensor[1,2] =Vec[5]#2,3
        Tensor[2,0] =Vec[4]#3,1
        Tensor[2,1] =Vec[5]#3,2
    elif  flag == 'LE':
        Tensor[0,1] =Vec[3]*0.5#1,2
        Tensor[0,2] =Vec[4]*0.5#1,3
        Tensor[1,0] =Vec[3]*0.5#2,1
        Tensor[1,2] =Vec[5]*0.5#2,3
        Tensor[2,0] =Vec[4]*0.5#3,1
        Tensor[2,1] =Vec[5]*0.5#3,2
    else:
        warnings.warn('Flag is unknown. Transformation. Conversion is performed without transforming the Vectorcomponents')
        Tensor[0,1] =Vec[3]#1,2
        Tensor[0,2] =Vec[4]#1,3
        Tensor[1,0] =Vec[3]#2,1
        Tensor[1,2] =Vec[5]#2,3
        Tensor[2,0] =Vec[4]#3,1
        Tensor[2,1] =Vec[5]#3,2
        print(Tensor)
    return Tensor
    
#+------------------------+
 
def getValueHistory(odb,flag,eID):
    'INPUt: odb,flag,eID'
    'Returns the time course of LE or S (flag) of the Elemente with the eID'
    histoValue = []
    frames = getFrames('Step-1',odb)
    histoValue=map(lambda temp:getValueONelement(temp,flag,eID),frames)
    return histoValue     

#+------------------------+

def getMaxEigVal(TENSOR):
    'Calculat the Eigevalues of a Tensor'
    EigVal=np.linalg.eigvals(TENSOR)
    maxValue=EigVal.max()
    return maxValue
 
#+------------------------+
  
def getMinEigVal(TENSOR):
    'Calculat the Eigevalues of a Tensor'
    EigVal=np.linalg.eigvals(TENSOR)
    minValue=EigVal.min()
    return minValue    

#+------------------------+
    
def calculateSWT(E,S_max,E_min,E_max):
    'Calculate the Smith-Wattson-Tropper-Parameter'
    'Input:E,S_max,E_min,E_max'
    SWF=np.sqrt((E_max-E_min)/2*S_max*E)
    return SWF

#+------------------------+

def calculateFS(E,k,S_max,S_yield,E_min,E_max):
    'Calculate the Fatemi-Socie-Paratmeter'
    'Input E,k,S_max,S_yield,E_min,E_max'
    FS = 0.5*(E_max - E_min)*(1 + k*S_max/S_yield) 
    return FS

#+-------------------------+
def getTimeMinMax(histoValue,flag):
    Max = np.max(map(getMaxEigVal,map(lambda temp : VectorToTensor(temp,flag),histoValue)))
    Min = np.min(map(getMinEigVal,map(lambda temp : VectorToTensor(temp,flag),histoValue)))
    return Min,Max
#+------------------------+

def calcuParameter(eIDs,histoLE,histoS,Para):
    'Calculate the Fatemi-Socie-Paratmeter and Smith-Wattson-Tropper-Parameter'
    'Returns the field for Countor plotting' 
    'Input Dir with the Parameter for Fatemi-Socie-Paratmeter and Smith-Wattson-Tropper-Parameter (E,K,S_yield)'
#    minLE,maxLE = 
#    minS,maxS = 
    FS = map(lambda temp:calculateSWT(Para,temp[2],temp[1],temp[0]),[minLE,maxLE,maxS])

#+------------------------+

    
def getDataForAreaOfIntrest(odb,instance_name,El_SetName=None):
    'generats the fieldOutput for the  with the pretended Function'
    'Input:frames,elements,func'
    'Return:eIDS,histS,histoLE'
    frames = getFrames('Step-1',odb)
    instance = odb.rootAssembly.instances[instance_name]
    eIDs=getEIDS(instance,El_SetName)
    histoS = map(lambda temp: getValueHistory(odb,'S',temp),range(len(eIDs)))
    histoLE = map(lambda temp: getValueHistory(odb,'LE',temp),range(len(eIDs)))
    return eIDs,histoS,histoLE

#+------------------------+
     
def ScalarNewFieldOutput(frame,Name,Field):
    'Genarete a new Fieldoutput for element'
    FieldOut =  frame.FieldOutput(name=Name,description=Name,type=SCALAR)
    
#+------------------------+
    
def exportVariable(Var,fileName='EXPORT.txt'):
    'export a variable in a txt-file'
    obj = open(fileName,'w')
    for var in Var:
        obj.write('%6.5f\n' % var)

#+------------------------+




#Run DEBUG Mode
odb=OPENodb('TEST','Shear_OneElement.odb')
frames = getFrames('Step-1',odb)
eIDS,histoS,histoLE = getDataForAreaOfIntrest(odb,odb.rootAssembly.instances['ONEELEMENT-1'].name)

Para = {'E':200000,'k':0.5,'S_yield':1200}





#histoValue = getValueHistory(odb,'S',0)
#print('LE-Wert von Element:1')
#print(histoValue)


#HistMaxEigVal = map(getMaxEigVal,map(lambda temp : VectorToTensor(temp,'S'),histoValue))
#exportVariable(HistMaxEigVal,fileName='EXPORT.txt')