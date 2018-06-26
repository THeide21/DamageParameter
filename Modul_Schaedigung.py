# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 11:38:17 2018

@author: ThomasHeidebrecht
"""
import time, math, warnings, sys, os
import numpy as np
import scipy as sc
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

def OPENodb(name_ODB,odbPathName,interactiveFlag): 
    'Opens ODB-File as readable'
    if interactiveFlag == 1:
        print('session')
        odb =session.openOdb(name=name_ODB, path=odbPathName, readOnly=False)
        
    else:
        print('openOdb')
        odb  = openOdb(path=odbPathName, readOnly=False)
        
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
        SET = instance.elementSets[SetName]
        elements=SET.elements
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
    var = frame.fieldOutputs[flag].values[eID-1].data
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
    elif  flag == 'LE' or flag == 'E':
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
 
def getValueHistory(odb,stepName,flag,eID):
    'INPUt: odb,flag,eID'
    'Returns the time course of LE or S (flag) of the Elemente with the eID'
    histoValue = []
    frames = getFrames(stepName,odb)
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
    
def getDataForAreaOfIntrest(odb,instance_name,stepName,El_SetName=None):
    'generats the fieldOutput for the  with the pretended Function'
    'Input:odb,instance_name,stepName,El_SetName=None'
    'Return:eIDS,histS,histoLE'
    frames = getFrames(stepName,odb)
    instance = odb.rootAssembly.instances[instance_name]
    eIDs=getEIDS(instance,El_SetName) 
    histoS = map(lambda temp: getValueHistory(odb,stepName,'S',temp),eIDs)
    try:
        histoLE = map(lambda temp: getValueHistory(odb,stepName,'LE',temp),eIDs)
    except:
        histoLE = map(lambda temp: getValueHistory(odb,stepName,'E',temp),eIDs)
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
    FieldOut =  frame.FieldOutput(name=Name,description=Name,type=SCALAR)
    FieldOut.addData(position=CENTROID,instance=instance,labels=eIDs,data=Field)
    odb.save()
    odb.close()    

#+------------------------+
def VectorNewFieldOutput(odb,instance,frame,Name,eIDs,Vec):
    FieldOut =  frame.FieldOutput(name=Name,description=Name,type=VECTOR,validInvariants=(MAGNITUDE,))
    FieldOut.addData(position=CENTROID,instance=instance,labels=eIDs,data=Vec)
    odb.save()
    odb.close()    
      
#+------------------------+

def getMAXbyRotation(histo):
    'Returns the MaxValue, the Postition, the Rotationmatrix'
    'and the frame of an time-depended Value for an element by rotating the Tensor'
    'Input: timedependen Valau (histo)'
    R_field = getRotatioField(-math.pi,math.pi,20)
    #MaxValaues = np.zeros([len(histo),1])
    #pos = np.zeros([len(histo),1])
    #R = np.zeros([len(histo),1])
    MaxValues = []
    pos = []
    R = []
    for frame in range(len(histo)):
        if TensorIsZero(histo[frame]) == True:
            MaxValues.append(0)
            pos.append([0,0])
            R.append(np.eye((3)))
        else:
            Tensor_roted  = map(lambda temp:rotT_pv(histo[frame], temp),R_field)
            hilfsTensor = map(getMaxTension,Tensor_roted)
            MaxValaue_roted = map(lambda temp:temp[0], hilfsTensor)
            pos_roted = map(lambda temp:temp[1], hilfsTensor)
            R_id = np.argmax(MaxValaue_roted)
            MaxValues.append(MaxValaue_roted[R_id])
            pos.append(pos_roted[R_id])
            R.append(R_field[R_id])#
    frame_id = np.argmax(MaxValues)
    #print(MaxValues)
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
    if tension[0] == S_max:
        pos = [0,0]
    elif tension[1] == S_max:
        pos = [1,1]
    elif tension[2] == S_max:
        pos = [2,2]
    else:
         print('Fehler ! in getMaxTension')       
    return S_max,pos

    
#+------------------------+

def getMinTension(Tensor):
    'Returns the minium tension'
    tension = np.zeros([3,1])
    tension[0] = Tensor[0,0]
    tension[1] = Tensor[1,1]
    tension[2] = Tensor[2,2]
    S_min = np.min(tension)
    pos = np.argwhere(Tensor == S_min)
    if tension[0] == S_min:
        pos = [0,0]
    elif tension[1] == S_min:
        pos = [1,1]
    elif tension[2] == S_min:
        pos = [2,2]
    else:
         print('Fehler ! in getMaxTension')      
    return S_min,pos

#+------------------------+

def getStrainForSWT(histo,R,pos): 
    hist_roted = map(lambda temp: rotT_pv(temp, R),histo)#
    histo_component = map(lambda temp: temp[pos[0],pos[1]],hist_roted)
    histo_min = np.min(histo_component)
    histo_max = np.max(histo_component)
    return histo_min, histo_max

 
#+------------------------+

def getStrainForFS(histo,R,pos): 
    hist_roted = map(lambda temp: rotT_pv(temp, R),histo)#
    if pos == [0,0]:
        histo_component_1 = map(lambda temp: temp[1,0],hist_roted)
        histo_component_2 = map(lambda temp: temp[2,0],hist_roted)
    elif pos == [1,1]:
        histo_component_1 = map(lambda temp: temp[0,1],hist_roted)
        histo_component_2 = map(lambda temp: temp[2,1],hist_roted)
    elif pos == [2,2]:
        histo_component_1 = map(lambda temp: temp[0,2],hist_roted)
        histo_component_2 = map(lambda temp: temp[1,2],hist_roted)
    else:
        print('Error in getStrainForFS')
    histo_min = np.min(np.min(histo_component_1),np.min(histo_component_2))
    histo_max = np.max(np.max(histo_component_1),np.max(histo_component_2))
    return histo_min, histo_max

#+------------------------+
    
def creatPlane(pos,R):
    'creat Plane'
    if pos == [0,0]:
        norm = R[:,0]
    elif pos ==[1,1]:
        norm = R[:,1]
    elif pos ==[2,2]:
        norm = R[:,2]
    else:
        print('Fehler in creatPlane')
    return tuple(norm)
        
#+------------------------+
    
def calculateSWT(histoS,histoLE,E):
    'Calculate the Smith-Wattson-Tropper-Parameter'
    S_max,pos,R,frame_id = getMAXbyRotation(histoS)
    E_SWT_min,E_SWT_max=getStrainForSWT(histoLE,R,pos)
    temp = (E_SWT_max-E_SWT_min)/2*S_max*E
    if temp > 0:
        SWT=np.sqrt(temp)
    else:
        SWT = 0
    SWT_vec=creatPlane(pos,R)
    return SWT,SWT_vec
        
#+------------------------+
    
def calculateFS(E,k,S_yield,S_max,E_amp):
    'Calculate the Fatemi-Socie-Paratmeter'
    FS = 0.5*(E_amp)*(1 + k*S_max/S_yield) 
    return FS

#-------------------------+

def getAmplitude(histo):
    'return the Amplitude of an timedepented-Variable'
    Max = np.max(histo)
    Min = np.min(histo)
    return Max-Min

#-------------------------+
    
def getMaxFSbyRotation(histoS,histoLE,E,S_yield,k):
    R_field = getRotatioField(-math.pi,math.pi,10)
    MaxValues = []
    pos = []
    for R in R_field:
        S_roted  = map(lambda temp: rotT_pv(temp, R),histoS)
        LE_roted  = map(lambda temp: rotT_pv(temp, R),histoLE) 
        n_1_max =np.max(map(lambda Tensor:Tensor[0,0],S_roted))  # Normalspannungen
        n_2_max=np.max(map(lambda Tensor:Tensor[1,1],S_roted)) # Normalspannungen
        n_3_max=np.max(map(lambda Tensor:Tensor[2,2],S_roted)) # Normalspannungen
        E_S_12 = getAmplitude(map(lambda temp:temp[0,1],LE_roted))    # Scherung 1,2
        E_S_13= getAmplitude(map(lambda temp:temp[0,2],LE_roted))    # Scherung 1,3
        E_S_23 = getAmplitude(map(lambda temp:temp[1,2],LE_roted))    # Scherung 2,3
        FS_temp = np.array([[calculateFS(E,k,S_yield,n_1_max,E_S_13),calculateFS(E,k,S_yield,n_1_max,E_S_12)],[calculateFS(E,k,S_yield,n_2_max,E_S_12),calculateFS(E,k,S_yield,n_2_max,E_S_23)],[calculateFS(E,k,S_yield,n_3_max,E_S_13),calculateFS(E,k,S_yield,n_3_max,E_S_23)]])
        FS_arg = np.unravel_index(np.argmax(FS_temp, axis=None), FS_temp.shape) 
        MaxValues.append(FS_temp[FS_arg]) 
        if FS_arg[0] == 0:
            pos.append([0,0])
        elif FS_arg[0] == 1:
            pos.append([1,1])
        elif FS_arg[0] == 2:
            pos.append([2,2])
        else:
            print('Fehler in getMaxFSbyRotation')
    arg = np.argmax(MaxValues)
    FS = MaxValues[arg]
    FS_vec = creatPlane(pos[arg],R_field[arg])
    return FS,FS_vec
#+------------------------+
def calcFatLifeSWT(SWT_list,eIDS,b,c,sigma_F,epsilion_F,E,portion = 0.05):
    if portion == 0 or len(SWT_list)==1:
        eIDS_short = eIDS
        SWT_short= SWT_list
    else:
        arg = np.argsort(np.array(SWT_list))
        eIDS_sort = np.array(eIDS)[arg]
        SWT_sort = np.array(SWT_list)[arg]
        eIDS_short = eIDS_sort[int(len(SWT_list)*(1-portion)):]
        SWT_short = SWT_sort[int(len(SWT_list)*(1-portion)):]
    N=[]
    for SWT in SWT_short:
        data = (SWT,sigma_F,epsilion_F,E,b,c)
        N.append(sc.optimize.fsolve(SWT_eq,1,data)[0])
    return eIDS_short,N

#+------------------------+
def SWT_eq(N,*data):
    SWT,sigma_F,epsilion_F,E,b,c = data
    return np.power(np.power(sigma_F,2)*np.power(2*N,2*b)+sigma_F*epsilion_F*E*np.power(2*N,(b+c)),0.5)-SWT
#+------------------------+
    
def makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_Name,StepName,interactiveFlag=0,setName = None):
    'This Function creats new Counter-Plots in ODB of Smith-Watson-Tropper and Fatemi-Socie -Paramater'
#    print('In Funktion:')
#    print('odb_name,E,k,S_yield,part_Name,interactiveFlag,setName')
#    print(odb_name,E,k,S_yield,part_Name,interactiveFlag,setName)
    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    frames = getFrames(StepName,odb)
    eIDS,histoS,histoLE = getDataForAreaOfIntrest(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'].name,StepName,setName)
    FS = []
    FS_vec = []
    SWT = []
    SWT_vec = []
    message = 'Berechnete Parameter'
    OB = 'Element'
    total = len(eIDS)
    print('Berechnung der Parameter')
    for jj in range(len(eIDS)):
        #print('%d von %d Elementen' % (jj,len(eIDS)))
        # Berchnung von SWT
        SWT_temp=calculateSWT(histoS[jj],histoLE[jj],E)
        SWT.append((SWT_temp[0],))
        SWT_vec.append(SWT_temp[1])
        #Berechnun von FS
        FS_temp = getMaxFSbyRotation(histoS[jj],histoLE[jj],E,S_yield,k)
        FS.append((FS_temp[0],))
        FS_vec.append(FS_temp[1])
        print(SWT_temp[0])
        print(SWT_temp[1])
        print(FS_temp[0])
        print(FS_temp[1])
    print('Counter Plotting')
#    ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'SWT',tuple(eIDS),tuple(SWT))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
#    VectorNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'SWT_Vec',tuple(eIDS),tuple(SWT_vec))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
#    ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'FS',tuple(eIDS),tuple(FS))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
#    VectorNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'FS_Vec',tuple(eIDS),tuple(FS_vec))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    #Ploting life-Time nach SWT
    # Einlesen SWT über GUI Später; 
    b = -0.113
    c = -610
    sigma_F = 1311
    epsilion_F = 0.8755
    portion = 1
    eIDs_SWT,N_SWT = calcFatLifeSWT(SWT,eIDS,b,c,sigma_F,epsilion_F,E,)
    ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'N_SWT',tuple(eIDs_SWT),tuple(N_SWT))
    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    #
#    #----------------------------------------------------------------------------------------------------------------------------
    try:
        ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'SWT',tuple(eIDS),tuple(SWT))
        odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    except:
        print('SWT exestiert')
    #----------------------------------------------------------------------------------------------------------------------------
    try:
        VectorNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'SWT_Vec',tuple(eIDS),tuple(SWT_vec))
        odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    except:
        print('SWT-Vec exestiert')
    #----------------------------------------------------------------------------------------------------------------------------   
        
    try:
        ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'FS',tuple(eIDS),tuple(FS))
        odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    except:
        print('FS exestiert')
    #----------------------------------------------------------------------------------------------------------------------------
    try:
        VectorNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'FS_Vec',tuple(eIDS),tuple(FS_vec))
        odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    except:
        print('FS-Vec exestiert')     
#    #----------------------------------------------------------------------------------------------------------------------------
     

    

#+------------------------+
    
def exportVariable(Var,fileName='EXPORT.txt'):
    'export a variable in a txt-file'
    obj = open(fileName,'w')
    for var in Var:
        obj.write('%6.5f\n' % var)
#+------------------------+