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

    
    
import time
import numpy as np
import warnings
import math
import itertools as it


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

# Calculates Rotation Matrix given euler angles.
def eulerAnglesToRotationMatrix(theta) :
    ''
    'From https://www.learnopencv.com/rotation-matrix-to-euler-angles/'
     
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
         
         
                     
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
                     
                     
    R = np.dot(R_z, np.dot( R_y, R_x ))
 
    return R

#+------------------------+

def getRotatioField(T_min,T_max,n):
    'Return all Ratationsmatrx for all angel(radian)'
    'input: Intial Angel: T_0, Maximal Angel T_max, Delta T'
        #print('New Theta Field will genarated an exported')
    T = np.linspace(T_min,T_max,n)
        # For numpy Version '1.13.3'
    #theta = np.array(np.meshgrid(T, T, T)).T.reshape(-1,3) # Matrix, which contains all combines of theta from 0,2*Pi
    theta = list(it.combinations_with_replacement(T,3))
    R = map(eulerAnglesToRotationMatrix,theta)
    return R

#+------------------------+

def rotT_pv(T, g):     # @pv.'s soln
    'Rotation of Tensor nach https://stackoverflow.com/questions/4962606/fast-tensor-rotation-with-numpy'
    return np.einsum('ac,bd,ab->cd', g, g, T)
    
#+------------------------+
 
def getValueHistory(odb,flag,eID):
    'INPUt: odb,flag,eID'
    'Returns the time course of LE or S (flag) of the Elemente with the eID'
    histoValue = []
    frames = getFrames('Step-1',odb)
    histoValue = map(lambda temp : VectorToTensor(temp,flag),map(lambda temp:getValueONelement(temp,flag,eID),frames))
    #histoValue=map(lambda temp : VectorToTensor(temp,flag),histoVec)
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
    temp = (E_max-E_min)/(2*S_max*E)
    if temp > 0:
        SWF=np.sqrt(temp)
    return SWF

#+------------------------+

def calculateFS(E,k,S_yield,S_max,E_min,E_max):
    'Calculate the Fatemi-Socie-Paratmeter'
    'Input E,k,S_max,S_yield,E_min,E_max'
    FS = 0.5*(E_max - E_min)*(1 + k*S_max/S_yield) 
    return FS

#+-------------------------+
    
def getTimeMin(histoValue,flag):
    Min = np.min(map(getMinEigVal,map(lambda temp : VectorToTensor(temp,flag),histoValue)))
    return Min

#+------------------------+

def getTimeMax(histoValue,flag):
    Max = np.max(map(getMaxEigVal,map(lambda temp : VectorToTensor(temp,flag),histoValue)))
    return Max

#+-------------------------+
    
def getTimeMinMax(histoValue,flag):
    Max = np.max(map(getMaxEigVal,map(lambda temp : VectorToTensor(temp,flag),histoValue)))
    Min = np.min(map(getMinEigVal,map(lambda temp : VectorToTensor(temp,flag),histoValue)))
    return Min,Max

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
    
def calcuParameter(histoLE,histoS,Para):
    'Calculate the Fatemi-Socie-Paratmeter and Smith-Wattson-Tropper-Parameter'
    'Returns the field for Countor plotting' 
    'Input Dir with the Parameter for Fatemi-Socie-Paratmeter and Smith-Wattson-Tropper-Parameter (E,K,S_yield)'
    #minLE = map(lambda temp: getTimeMin(temp,'LE'),histoLE)
    #maxLE = map(lambda temp: getTimeMax(temp,'LE'),histoLE)
    #maxS = map(lambda temp: getTimeMax(temp,'S'),histoS)
    maxS,pos = getMAXbyRotation(histoS)
    SWT =[]
    FS = []
    for i in range(len(minLE)):
        SWT.append((calculateSWT(Para['E'],maxS[i],minLE[i],maxLE[i]),))
        FS.append((calculateFS(Para['E'],Para['k'],Para['S_yield'],maxS[i],minLE[i],maxLE[i]),))
    return FS,SWT

#+------------------------+
    
def ScalarNewFieldOutput(odb,instance,frame,Name,eIDs,Field):
    'Genarete a new Fieldoutput for element'
    'Attention: OBD will be saved and closed. It has to be opend again'
    'Input:odb,instance,frame,Name and eIDs and Field as tuple'
    try:
        FieldOut =  frame.FieldOutput(name=Name,description=Name,type=SCALAR)
        FieldOut.addData(position=INTEGRATION_POINT,instance=instance,labels=eIDS,data=Field)
    except:
        print('Variable %s already exist' % Name)
    odb.save()
    odb.close()    
    
#+------------------------+

def getMAXbyRotation(histo):
    'Returns the MaxValue, the Postition, the Rotationmatrix'
    'and the frame of an time-depended Value for an element by rotating the Tensor'
    'Input: timedependen Valau (histo)'
    R_field = getRotatioField(0,2*math.pi,200)
    #MaxValaues = np.zeros([len(histo),1])
    #pos = np.zeros([len(histo),1])
    #R = np.zeros([len(histo),1])
    MaxValues = []
    pos = []
    R = []
    for frame in range(len(histo)):
        print(frame)
        print
        if TensorIsZero(histo[frame]) == True:
            MaxValues.append(0)
            pos.append([0,0])
            R.append(0)
            print('Nur Nullen')
        else:
            Tensor_roted  = map(lambda temp:rotT_pv(histo[frame], temp),R_field)
            hilfsTensor = map(getMaxTension,Tensor_roted)
            MaxValaue_roted = map(lambda temp:temp[0], hilfsTensor)
            pos_roted = map(lambda temp:temp[1], hilfsTensor)
            R_id = np.argmax(MaxValaue_roted)
            MaxValues.append(MaxValaue_roted[R_id])
            pos.append(pos_roted[R_id])
            R.append(R_field[R_id])#
            print('Vektor wird rotiert')
    frame_id = np.argmax(MaxValues)
    print('---------------------------------------')
    print('Max Values')
    print(MaxValues)
    print('---------------------------------------')
    #print(pos)
    #print(R)
    return MaxValues[frame_id],pos[frame_id],R[frame_id],frame_id

#+------------------------+

def TensorIsZero(Tensor):
    if np.count_nonzero(Tensor) != 9:
        temp = True
    else:
        temp = False
    return temp

#+------------------------+
def getMaxTension(Tensor):
    'Returns the maximum tension'
    tension = np.zeros([3,1])
    tension[0] = Tensor[0,0]
    tension[1] = Tensor[1,1]
    tension[2] = Tensor[2,2]
    S_max = np.max(tension)
    pos = np.argwhere(Tensor == S_max)
    return S_max,pos

    

#+------------------------+
    
def exportVariable(Var,fileName='EXPORT.txt'):
    'export a variable in a txt-file'
    obj = open(fileName,'w')
    for var in Var:
        obj.write('%6.5f\n' % var)
#+------------------------+


debug = 1

if debug ==1:
    print('Start')
    from textRepr import * 
    import visualization
    from abaqusConstants import *
    import math
    import abaqus
    from odbAccess import *
    #from odbAccess import *
#Run DEBUG Mode
    #odb=OPENodb('TEST','Shear_OneElement.odb')
    #odb=OPENodb('TEST','Benchmark_Fine.odb')
    odb=OPENodb('TEST','Lochscheibe_Fine.odb')
    frames = getFrames('Step-1',odb)
    eIDS,histoS,histoLE = getDataForAreaOfIntrest(odb,odb.rootAssembly.instances['LOCHSCHEIBE_3D-1'].name)
    Para = {'E':200000,'k':0.5,'S_yield':1200}
    histo = histoS[-1]
    MaxValue,pos,R,frame_id = getMAXbyRotation(histo)
    print(MaxValue)
    print(pos)
    print(R)
    print(frame_id)
    #FS,SWT = calcuParameter(histoLE,histoS,Para)
    
#    ScalarNewFieldOutput(odb,odb.rootAssembly.instances['LOCHSCHEIBE_3D-1'],frames[-1],'SWT',tuple(eIDS),tuple(SWT))
#    odb=OPENodb('TEST','Benchmark_Coarse.odb')
#    ScalarNewFieldOutput(odb,odb.rootAssembly.instances['LOCHSCHEIBE_3D-1'],frames[-1],'FS',tuple(eIDS),tuple(FS))
#    odb=OPENodb('TEST','Benchmark_Coarse.odb')



#histoValue = getValueHistory(odb,'S',0)
#print('LE-Wert von Element:1')
#print(histoValue)


#HistMaxEigVal = map(getMaxEigVal,map(lambda temp : VectorToTensor(temp,'S'),histoValue))
#exportVariable(HistMaxEigVal,fileName='EXPORT.txt')