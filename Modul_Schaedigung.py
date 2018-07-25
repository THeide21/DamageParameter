# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 11:38:17 2018

@author: ThomasHeidebrecht
"""
import time, math, warnings, sys, os
import numpy as np
#import scipy as sc
import itertools as it
import logging
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
    if interactiveFlag == 0:
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
 
def getEIDS(odb,instance_name,El_SetName = None):
    instance = odb.rootAssembly.instances[instance_name]    
    if El_SetName == None:
        elements = instance.elements
    else:
        SET = instance.elementSets[El_SetName]
        print('Set wird auf Assembly-Ebene geholt')
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


# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6
 
 
# Calculates rotation matrix to euler angles
# The result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
    
def rotationMatrixToEulerAngles(R) :
    assert(isRotationMatrix(R))
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
    return np.array([x, y, z])

#+------------------------+

def getRotatioField(T_min,T_max,n):
    'Return all Ratationsmatrx for all angel(radian)'
    'input: Intial Angel: T_0, Maximal Angel T_max, Delta T'
        #print('New Theta Field will genarated an exported')
    T = np.linspace(T_min,T_max,n)
    theta = list(it.combinations_with_replacement(T,10))
    #theta = [[0,0,0]]#,[0,0,0],[0,0.7853981634,0],[0,0,0.7853981634]]
    R = map(eulerAnglesToRotationMatrix,theta)
    return R,theta

#+------------------------+

def rotT_pv(T, g):     # @pv.'s soln
    'Rotation of Tensor nach https://stackoverflow.com/questions/4962606/fast-tensor-rotation-with-numpy'
    return np.einsum('ac,bd,ab->cd', g, g, T)
    
#+------------------------+
#odb,stepName,'S',temp,t,_begin,t_end,frame_flag
def getValueHistory(odb,stepName,flag,eID,t_begin,t_end,frame_flag='all'):
    'INPUt: odb,flag,eID'
    'Returns the time course of LE or S (flag) of the Elemente with the eID'
    frames = getFrames(stepName,odb)
#    try:
    if frame_flag == 'all':
        histoValue = map(lambda temp : VectorToTensor(temp,flag),map(lambda temp:getValueONelement(temp,flag,eID),frames))
    else:
        histoValue = []
        frame_begin=getNearstFrame(t_begin,frames)
        frame_end=getNearstFrame(t_end,frames)
        print(frame_begin)
        print(frame_end)
        if frame_flag == 'interval':
            frame_interval = range(frame_begin,frame_end+1)
        else:
            frame_interval = [frame_begin,frame_end]
        for frame_id in frame_interval: 
            histoValue.append(VectorToTensor(getValueONelement(frames[frame_id],flag,eID),flag))           
#    except:
#        print('Fehler in getValueHistory')
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
    
def getDataForAreaOfIntrest(odb,instance_name,stepName,t_begin,t_end,frame_flag='all',El_SetName=None):
    'generats the fieldOutput for the  with the pretended Function'
    'Input:odb,instance_name,stepName,El_SetName=None'
    'Return:eIDS,histS,histoLE'
    frames = getFrames(stepName,odb)
    eIDs=getEIDS(odb,instance_name,El_SetName) 
    histoS = map(lambda temp: getValueHistory(odb,stepName,'S',temp,t_begin,t_end,frame_flag),eIDs)
    try:
        histoLE = map(lambda temp: getValueHistory(odb,stepName,'LE',temp,t_begin,t_end,frame_flag),eIDs)
    except:
        histoLE = map(lambda temp: getValueHistory(odb,stepName,'E',temp,t_begin,t_end,frame_flag),eIDs)
    return eIDs,histoS,histoLE

#+------------------------+
    
def calcuParameter(histoLE,histoS,Para):
    'Calculate the Fatemi-Socie-Paratmeter and Smith-Wattson-Tropper-Parameter'
    'Returns the field for Countor plotting' 
    'Input Dir with the Parameter for Fatemi-Socie-Paratmeter and Smith-Wattson-Tropper-Parameter (E,K,S_yield)'
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
    R_field,theta = getRotatioField(-math.pi,math.pi,20)
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

def derivative(f):
    'https://www.daniweb.com/programming/software-development/code/444930/newton-s-method-example-python'
    def compute(x, dx):
        return (f(x+dx) - f(x))/dx
    return compute

#+------------------------+
def newtons_method(f, x, dx=0.000001, tolerance=0.000001):
    '''f is the function f(x)'''
    'https://www.daniweb.com/programming/software-development/code/444930/newton-s-method-example-python'
    df = derivative(f)
    while True:
        x1 = x - f(x)/df(x, dx)
        t = abs(x1 - x)
        if t < tolerance:
            break
        x = x1
    return x

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
    FS = (E_amp)*(1 + k*S_max/S_yield)
    return FS

#-------------------------+

def getAmplitude(histo):
    'return the Amplitude of an timedepented-Variable'
    Max = np.max(histo)
    Min = np.min(histo)
    print('Min: %f \n Max: %f' % (Min,Max))
    return Max-Min

#-------------------------+
    
def getMaxFSbyRotation(histoS,histoLE,E,S_yield,k):
    R_field,theta = getRotatioField(0,math.pi,10)
    MaxValues = []
    pos = []
    for R in R_field:
        S_roted  = map(lambda temp: rotT_pv(temp, R),histoS)
        LE_roted  = map(lambda temp: rotT_pv(temp, R),histoLE)
        n_1_max =np.max(map(lambda Tensor:Tensor[0,0],S_roted))  # Normalspannungen
        n_2_max=np.max(map(lambda Tensor:Tensor[1,1],S_roted)) # Normalspannungen
        n_3_max=np.max(map(lambda Tensor:Tensor[2,2],S_roted)) # Normalspannungen
        print('1-2')
        E_S_12 = getAmplitude(map(lambda temp:temp[0,1],LE_roted))    # Scherung 1,2
        print('1-3')
        E_S_13= getAmplitude(map(lambda temp:temp[0,2],LE_roted))    # Scherung 1,3
        print('2_3')
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
    print('FS max: %f' % FS)
    FS_vec = creatPlane(pos[arg],R_field[arg])
    print(rotationMatrixToEulerAngles(R_field[arg]))
    return FS,FS_vec

#+------------------------+
    
def reduceEIDs(eIDS,variable_List,portion):
    print('reduceEIDs:')
    print(eIDS)
    print(variable_List)
    if len(variable_List)==1:
        eIDS_short = eIDS
        short_variable= variable_List
    else:
        arg = np.argsort(np.array(variable_List))
        print(arg)
        eIDS_sort = np.array(eIDS)[arg]
        variable_sort = np.array(variable_List)[arg]
        print('eiDs:----------------------')
        print(eIDS_sort)
        print('SWT: -----------------------')
        print(variable_sort)
        if portion == 0: 
            print('Worst case Element ')
            eIDS_short = eIDS_sort[-1]
            short_variable = variable_sort[-1]
        else:
            eIDS_short = eIDS_sort[int(len(variable_List)*(1-portion)):]
            short_variable = variable_sort[int(len(variable_List)*(1-portion)):]
    return eIDS_short[0],short_variable[0]

#*------------------------+
def calcFatLifeSWT(SWT_list,eIDS,b,c,sigma_F,epsilion_F,E,portion = 0.05):
    eIDS_short,SWT_short = reduceEIDs(eIDS,SWT_list,portion)
    N=[]
    for SWT in SWT_short:
        N.append(newtons_method(lambda N_temp: SWT_eq(N_temp,SWT,sigma_F,epsilion_F,E,b,c), 1))
    return eIDS_short,N

#+------------------------+
    
def SWT_eq(N,SWT,sigma_F,epsilion_F,E,b,c):  
    return np.power(np.power(sigma_F,2)*np.power(2*N,2*b)+sigma_F*epsilion_F*E*np.power(2*N,(b+c)),0.5)-SWT

#+------------------------+
    
def FS_eq(N,SWT,tau_F,gemma_F,G,b,c):
    return tau_F

#+------------------------+
    
def getNearstFrame(t,frames):
    'Returns the Index of the nearst Fra,e'
    return np.argmin(map(lambda frame: np.abs(frame.frameValue - t),frames))

#+------------------------+
#+------------------------+
#+------------------------+   

def makeCounterPlotof_FS_SWT(odb_name,E,k,S_yield,part_Name,StepName,t_begin,t_end,b,c,sigma_F,epsilion_F,portion,interactiveFlag=0,setName = None,frame_flag='all'):
    'This Function creats new Counter-Plots in ODB of Smith-Watson-Tropper and Fatemi-Socie -Paramater'
    sys.stdout = open(odb_name+'.txt', 'w')
    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    frames = getFrames(StepName,odb) 
    eIDS,histoS,histoLE = getDataForAreaOfIntrest(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'].name,StepName,t_begin,t_end,frame_flag,setName)
    FS = []
    FS_vec = []
    SWT = []
    SWT_vec = []
    print('Berechnung der Parameter')
    for jj in range(len(eIDS)):
         print('Elementlabel %d von %d' % (eIDS[jj],len(eIDS)-1 ))
        #Berchnung von SWT
         SWT_temp=calculateSWT(histoS[jj],histoLE[jj],E)
         SWT.append((SWT_temp[0],))
         SWT_vec.append(SWT_temp[1])
        #Berechnun von FS
         FS_temp = getMaxFSbyRotation(histoS[jj],histoLE[jj],E,S_yield,k)
         FS.append((FS_temp[0],))
         FS_vec.append(FS_temp[1])

	# eIDs_SWT,N_SWT = ca lcFatLifeSWT(SWT,eIDS,b,c,sigma_F,epsilion_F,E,portion)  
    # #try: 
	# ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'N_SWT',tuple(eIDs_SWT),tuple(N_SWT))
    # odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
    # #except:
     # #   print('Lebensdauer bereits berechnet ')
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
        
#    ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'SWT',tuple(eIDS),tuple(SWT))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
#    VectorNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'SWT_Vec',tuple(eIDS),tuple(SWT_vec))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
#    ScalarNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'FS',tuple(eIDS),tuple(FS))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)
#    VectorNewFieldOutput(odb,odb.rootAssembly.instances[part_Name.upper()+'-1'],frames[-1],'FS_Vec',tuple(eIDS),tuple(FS_vec))
#    odb=OPENodb(odb_name,odb_name+'.odb',interactiveFlag)