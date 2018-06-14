# Do not edit this file or it may not load correctly
# if you try to open it with the RSG Dialog Builder.

# Note: thisDir is defined by the Activator class when
#       this file gets exec'd

from rsg.rsgGui import *
from abaqusConstants import INTEGER, FLOAT
dialogBox = RsgDialog(title='Berechnung von SWT und FS-Sch\xe4digungsparameter', kernelModule='Schaediungsparameter', kernelFunction='makeCounterPlotof_FS_SWT', includeApplyBtn=False, includeSeparator=True, okBtnText='OK', applyBtnText='Apply', execDir=thisDir)
RsgGroupBox(name='GroupBox_2', p='DialogBox', text='Smith-Watson-Topper', layout='0')
RsgIcon(p='GroupBox_2', fileName=r'SWT.PNG')
RsgGroupBox(name='GroupBox_3', p='DialogBox', text='Fatemi-Socie', layout='0')
RsgIcon(p='GroupBox_3', fileName=r'FS.PNG')
RsgGroupBox(name='GroupBox_1', p='DialogBox', text='Materialparameter', layout='0')
RsgTextField(p='GroupBox_1', fieldType='Float', ncols=12, labelText='E-Modul [GPa]', keyword='E', default='')
RsgTextField(p='GroupBox_1', fieldType='Float', ncols=12, labelText='Streckgrenze [MPa]', keyword='S_yield', default='')
RsgSlider(p='GroupBox_1', text='k', minLabelText='0', maxLabelText='1', valueType=FLOAT, minValue=0, maxValue=1, decimalPlaces=2, showValue=True, width=200, keyword='k', default=0.5)
RsgSeparator(p='GroupBox_1')
RsgTextField(p='GroupBox_1', fieldType='String', ncols=12, labelText='Odb-Name', keyword='odb_name', default='')
RsgTextField(p='GroupBox_1', fieldType='String', ncols=12, labelText='Part-Name', keyword='part_Name', default='')
dialogBox.show()