from abaqusConstants import *
from abaqusGui import *
from kernelAccess import mdb, session
import os

thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)


###########################################################################
# Class definition
###########################################################################

class GUI_ParamteterDB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #

        AFXDataDialog.__init__(self, form, 'Smith-Watson-Topper und  Fatemi-Socie  Sch\xe4digungsparameter',
            self.OK|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
            

        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
            
        GroupBox_1 = FXGroupBox(p=self, text='Parameter', opts=FRAME_GROOVE|LAYOUT_FILL_X)
        AFXTextField(p=self, ncols=20, labelText='E [GPa]', tgt=form.E_paraKw, sel=0)
        if isinstance(self, FXHorizontalFrame):
            FXVerticalSeparator(p=self, x=0, y=0, w=0, h=0, pl=2, pr=2, pt=2, pb=2)
        else:
            FXHorizontalSeparator(p=self, x=0, y=0, w=0, h=0, pl=2, pr=2, pt=2, pb=2)
        AFXTextField(p=self, ncols=12, labelText='Streckgrenze [MPa]', tgt=form.S_yield_paraKw, sel=0)
        slider = AFXSlider(self, form.k_parameterKw, 0, 
            AFXSLIDER_INSIDE_BAR|AFXSLIDER_SHOW_VALUE|LAYOUT_FIX_WIDTH, 0,0,200,0)
        slider.setTitleLabelText('k')
        slider.setTitleLabelJustify(JUSTIFY_CENTER_X)
        slider.setMinLabelText('0')
        slider.setMaxLabelText('1')
        slider.setDecimalPlaces(2)
        slider.setRange(0, 100)
