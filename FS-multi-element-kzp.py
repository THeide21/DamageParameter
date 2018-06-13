#
#+----------------------+
#|   block of modules   |
#+----------------------+
#	
from odbAccess import *
from abaqusConstants import *
#from Numeric import array
import math
import abaqus
import time
#
# To measure timing
#' http://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
def tic():
    # Homemade version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    if 'startTime_for_tictoc' in globals():
        print( "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")
#+-------------------------------+
#|   block of user's constants   |
#+-------------------------------+
tic()
# strain_label = "E" or "LE"
strain_label = 'LE'
#
# Drehung um z-Achse
ALPHA = [0,10,20,25,30,35,40,45,50,60,70,80,90,100,110,120,125,130,131,135,140,150,160,165,170,175]
#ALPHA = [30]
BETA = [90] 
#
odbPathName = 'C:\Users\ThomasHeidebrecht\Documents\Abaqus\Lochscheibe_Fine.odb'
oldStepName = 'Step-1'
#
#
#+------------------------+
#|   block of functions   |
#+------------------------+
#
def TransformField(C, alpha, beta, flag):
    t = [[math.cos(alpha)*math.sin(beta), math.sin(alpha)*math.sin(beta), math.cos(beta)],\
         [-math.sin(alpha), math.cos(alpha),0],\
         [-math.cos(alpha)*math.cos(beta), -math.sin(alpha)*math.cos(beta),math.sin(beta)]]
    #
    B = [[t[0][0]*t[0][0], t[0][1]*t[0][1], t[0][2]*t[0][2], 2*t[0][0]*t[0][1], 2*t[0][0]*t[0][2], 2*t[0][2]*t[0][1]],\
         [t[1][0]*t[1][0], t[1][1]*t[1][1], t[1][2]*t[1][2], 2*t[1][0]*t[1][1], 2*t[1][0]*t[1][2], 2*t[1][2]*t[1][1]],\
         [t[2][0]*t[2][0], t[2][1]*t[2][1], t[2][2]*t[2][2], 2*t[2][0]*t[2][1], 2*t[2][0]*t[2][2], 2*t[2][2]*t[2][1]],\
         [t[0][0]*t[1][0], t[0][1]*t[1][1], t[0][2]*t[1][2], t[0][0]*t[1][1]+t[0][1]*t[1][0], t[0][2]*t[1][0]+t[0][0]*t[1][2], t[0][1]*t[1][2]+t[0][2]*t[1][1]],\
         [t[0][0]*t[2][0], t[0][1]*t[2][1], t[0][2]*t[2][2], t[0][0]*t[2][1]+t[0][1]*t[2][0], t[0][2]*t[2][0]+t[0][0]*t[2][2], t[0][2]*t[2][1]+t[0][1]*t[2][2]],\
         [t[1][0]*t[2][0], t[1][1]*t[2][1], t[1][2]*t[2][2], t[1][0]*t[2][1]+t[1][1]*t[2][0], t[1][2]*t[2][0]+t[1][0]*t[2][2], t[1][1]*t[2][2]+t[1][2]*t[2][1]]]
    #
    A = [0.0,0.0,0.0,0.0,0.0,0.0]
    G = [0.0,0.0,0.0,0.0,0.0,0.0]
    #
    for i in [0,1,2]:
        G[i] = C[i]
    for i in [3,4,5]:
        if flag == 'strain':
           G[i] = C[i]/2.0
        else:
           G[i] = C[i]        
    #
    for i in [0,1,2,3,4,5]:
        for j in [0,1,2,3,4,5]:
            A[i] = A[i]+B[i][j]*G[j]      
    #
    if flag == 'strain':
       for i in [3,4,5]:
           A[i] = A[i]*2.0
    #    
    return A 
#
#
def FindCycle(frames, Step):
    frames_of_Cycle = []
    for i in frames:
        frames_of_Cycle.append(Step.frames[i])
    return frames_of_Cycle
#
#
def newSWT(ALPHA,BETA,frames):
    SWT_Table = [] 
    #
    element_range = range(0,len(frames[0].fieldOutputs['S'].values))
    for element_id in element_range:
        oldStress_Table = []
        oldStrain_Table = []
        for frame in frames:
            stressElement = frame.fieldOutputs['S'].values[element_id]
            strainElement = frame.fieldOutputs[strain_label].values[element_id]
            stressTemp = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            strainTemp = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            component_range = range(0,len(stressElement.data))
            for component in component_range:
                stressTemp[component] = stressElement.data[component]
                strainTemp[component] = strainElement.data[component]
            #
            oldStress_Table.append(stressTemp)
            oldStrain_Table.append(strainTemp)
        #
        for beta in BETA:
            for alpha in ALPHA:
                newStress_Table = []
                newStrain_Table = []
                for i in range(0,len(oldStress_Table)):
                    newStress_Table.append(TransformField(oldStress_Table[i],abaqus.degreeToRadian(alpha),abaqus.degreeToRadian(beta),''))
                    newStrain_Table.append(TransformField(oldStrain_Table[i],abaqus.degreeToRadian(alpha),abaqus.degreeToRadian(beta),'strain'))
                #
                Sx = []
                Ex = []
                Temp_range = range(0,len(newStress_Table))
                for i in Temp_range:
                    Sx.append(newStress_Table[i][0])
                    Ex.append(newStrain_Table[i][0])
                #
                Sx.sort()
                Ex.sort()
                Smax = Sx[-1]
                Emax = Ex[-1]
                Emin = Ex[0]
                SWT = Smax*(Emax-Emin)/2.0*2.00e5/1000000
                SWT_Table.append([frame.fieldOutputs['S'].values[element_id].elementLabel,\
                                  frame.fieldOutputs['S'].values[element_id].integrationPoint,\
                                  alpha,beta,SWT]) 
				#
			#	
        #
    #
    return SWT_Table
#
#
def export_results(data, fileName='results.txt'):
    file = open(''+fileName,'w+')
    for line in range(0,len(data)):
        file.write(str(data[line][0])+', '+str(data[line][1])+', '+str(data[line][2])+', '+str(data[line][3])+', '+format(data[line][4], '1.6f')+chr(10))
    #
    file.close()
    return
#
#
def FatemiSocie(ALPHA,BETA,frames,k,Syield):
    #k = 0.4
    #Syield = 500000000 #[MPa]
    FS_Table = [] 
    #
    element_range = range(0,len(frames[0].fieldOutputs['S'].values))
    for element_id in element_range:
        oldStress_Table = []
        oldStrain_Table = []
        for frame in frames:
            stressElement = frame.fieldOutputs['S'].values[element_id]
            strainElement = frame.fieldOutputs[strain_label].values[element_id]
            stressTemp = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            strainTemp = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            component_range = range(0,len(stressElement.data))
            for component in component_range:
                stressTemp[component] = stressElement.data[component]
                strainTemp[component] = strainElement.data[component]
            #
            oldStress_Table.append(stressTemp)
            oldStrain_Table.append(strainTemp)        
        #
        for beta in BETA:
            for alpha in ALPHA:
                newStress_Table = []
                newStrain_Table = []
                for i in range(0,len(oldStress_Table)):
                    newStress_Table.append(TransformField(oldStress_Table[i],abaqus.degreeToRadian(alpha),abaqus.degreeToRadian(beta),''))
                    newStrain_Table.append(TransformField(oldStrain_Table[i],abaqus.degreeToRadian(alpha),abaqus.degreeToRadian(beta),'strain'))                
                #
                Sx = []
                Exy = []
                Temp_range = range(0,len(newStress_Table))
                for i in Temp_range:
                    Sx.append(newStress_Table[i][0])
                    Exy.append(newStrain_Table[i][3])
                #
                Sx.sort()
                Exy.sort()
                Smax = Sx[-1]
                Emax = Exy[-1]
                Emin = Exy[0]
                FS = 0.5*(Emax - Emin)*(1 + k*Smax/Syield)
                FS_Table.append([frame.fieldOutputs['S'].values[element_id].elementLabel,\
                                 frame.fieldOutputs['S'].values[element_id].integrationPoint,\
                                 alpha,beta,FS])
				#
			#	
        #
    #
    return FS_Table
#
#
def findMax(fatigue_Parameter, maxFP):
    maxFP.append(fatigue_Parameter[0])
    for i in range(1,len(fatigue_Parameter)):
        if fatigue_Parameter[i][0] == maxFP[-1][0]:
           if fatigue_Parameter[i][4] > maxFP[-1][4]:
              maxFP[-1] = fatigue_Parameter[i]
        else:
           maxFP.append(fatigue_Parameter[i]) 
#
#
#
#
#+----------------+
#|   Main block   |
#+----------------+
#
oldOdb = openOdb(path=odbPathName)
oldRootAssembly = oldOdb.rootAssembly
oldStep = oldOdb.steps[oldStepName]
#Select frames 
#frames_id = [72,73,74,75,76,77,78,79,80]
#frames_id = [712,713,714,715,716,717,718,719,720]
frames_id = [16,17,18,19,20]
frames = FindCycle(frames_id, oldStep)
#
maxFP = []
#
#fatigue_Parameter = newSWT(ALPHA,BETA,frames) # Quadrat von SWT wird berechnet !!!
#findMax(fatigue_Parameter, maxFP)
#export_results(fatigue_Parameter, fileName='SWT-results.txt')
#print fatigue_Parameter
#
k = 0.4
Syield = 356000000 #[MPa]
fatigue_Parameter = FatemiSocie(ALPHA,BETA,frames,k,Syield)
findMax(fatigue_Parameter, maxFP)
export_results(fatigue_Parameter, fileName='FS-results.txt')
print(fatigue_Parameter)
#
#print maxFP
from operator import itemgetter
#print (sorted(maxFP, key=itemgetter(4)))
#
print('FERTIG')
#
toc()