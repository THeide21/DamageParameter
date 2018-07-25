# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 13:00:16 2018

@author: ThomasHeidebrecht
"""

def startupGUI(odb_name,part_name,E,k,S_yield,StepName,b,c,sigma_F,epsilion_F,portion,t_begin,t_end,frame_flag,SetName=None):
    print('odb_name,part_name,E,k,S_yield,StepName,SetName')
    print(odb_name,part_name,E,k,S_yield,StepName,SetName)
    import os
    import subprocess as sp
    if SetName == None or SetName == 'None':
        print('abaqus python Schaediungsparameter.py %f %f %f %s %s %s %f %f %f %f %f %f %f %s'%(E,k,S_yield,odb_name,part_name,StepName,b,c,sigma_F,epsilion_F,portion,t_begin,t_end,frame_flag))
        var = sp.call('abaqus python Schaediungsparameter.py %f %f %f %s %s %s %f %f %f %f %f %f %f %s'%(E,k,S_yield,odb_name,part_name,StepName,b,c,sigma_F,epsilion_F,portion,t_begin,t_end,frame_flag),shell=True)
    else:
        print('abaqus python Schaediungsparameter.py %f %f %f %s %s %s %f %f %f %f %f %f %f %s %s'%(E,k,S_yield,odb_name,part_name,StepName,t_begin,t_end,frame_flag,SetName))
        var = sp.call('abaqus python Schaediungsparameter.py %f %f %f %s %s %s %f %f %f %f %f %f %f %s %s'%(E,k,S_yield,odb_name,part_name,StepName,b,c,sigma_F,epsilion_F,portion,t_begin,t_end,frame_flag,SetName),shell=True)
    if var == 1: 
        print('Auswertung könnte nicht durchgeführt werden')
    else:
        print('Auswertung wurde durchgeführt')
    return 