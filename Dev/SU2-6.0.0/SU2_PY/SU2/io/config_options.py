## \file config_options.py
#  \brief python package for config
#  \author T. Lukaczyk, F. Palacios
#  \version 6.0.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from ..util import ordered_bunch

class OptionError(Exception):
    pass

class Option(object):
    
    def __init__(self):
        self.val = ""

    def __get__(self):
        return self.val

    def __set__(self,newval):
        self.val = newval

#: class Option

class MathProblem(Option):

    def __init__(self,*args,**kwarg):
        super(MathProblem,self).__init__(*args,**kwarg)
        self.validoptions = ['DIRECT','CONTINUOUS_ADJOINT','LINEARIZED']

    def __set__(self,newval):
        if not self.newval in self.validoptions:
            raise OptionError("Invalid option. Valid options are: %s"%self.validoptions)
        super(MathProblem,self).__set__(newval)

#: class MathProblem

class DEFINITION_DV(ordered_bunch):
    """ SU2.io.config.DEFINITION_DV()
    
        List of design variables (Design variables are separated by semicolons)
        - HICKS_HENNE ( 1, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc )
        - SURFACE_BUMP ( 2, Scale | Mark. List | x_Start, x_End, x_Loc )
        - NACA_4DIGITS ( 4, Scale | Mark. List |  1st digit, 2nd digit, 3rd and 4th digit )
        - TRANSLATION ( 5, Scale | Mark. List | x_Disp, y_Disp, z_Disp )
        - ROTATION ( 6, Scale | Mark. List | x_Axis, y_Axis, z_Axis, x_Turn, y_Turn, z_Turn )
        - FFD_CONTROL_POINT ( 7, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind, k_Ind, x_Mov, y_Mov, z_Mov )
        - FFD_TWIST_ANGLE ( 9, Scale | Mark. List | FFD_Box_ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
        - FFD_ROTATION ( 10, Scale | Mark. List | FFD_Box_ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
        - FFD_CAMBER ( 11, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind )
        - FFD_THICKNESS ( 12, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind )
        - FFD_CONTROL_POINT_2D ( 15, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind, x_Mov, y_Mov )
        - FFD_CAMBER_2D ( 16, Scale | Mark. List | FFD_Box_ID, i_Ind )
        - FFD_THICKNESS_2D ( 17, Scale | Mark. List | FFD_Box_ID, i_Ind )
        
    """
    
    def __init__(self,*args,**kwarg):
        ordered_bunch.__init__(self)
        self.KIND   = []
        self.SCALE  = []
        self.MARKER = []
        self.FFDTAG = []
        self.PARAM  = []
        self.update(ordered_bunch(*args,**kwarg))
    
    def append(self,new_dv):
        self.KIND.  append(new_dv['KIND'])
        self.SCALE. append(new_dv['SCALE'])
        self.MARKER.append(new_dv['MARKER'])
        self.FFDTAG.append(new_dv['FFDTAG'])
        self.PARAM. append(new_dv['PARAM'])
    
    def extend(self,new_dvs):
        assert isinstance(new_dvs,DEFINITION_DV) , 'input must be of type DEFINITION_DV'
        self.KIND.  extend(new_dvs['KIND'])
        self.SCALE. extend(new_dvs['SCALE'])
        self.MARKER.extend(new_dvs['MARKER'])
        self.FFDTAG.extend(new_dvs['FFDTAG'])
        self.PARAM. extend(new_dvs['PARAM'])

#: class DEFINITION_DV

class DV_KIND(ordered_bunch):
  """ SU2.io.config.DV_KIND()
    
    List of design variables (Design variables are separated by semicolons)
    - HICKS_HENNE ( 1, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc )
    - SURFACE_BUMP ( 2, Scale | Mark. List | x_Start, x_End, x_Loc )
    - NACA_4DIGITS ( 4, Scale | Mark. List |  1st digit, 2nd digit, 3rd and 4th digit )
    - TRANSLATION ( 5, Scale | Mark. List | x_Disp, y_Disp, z_Disp )
    - ROTATION ( 6, Scale | Mark. List | x_Axis, y_Axis, z_Axis, x_Turn, y_Turn, z_Turn )
    - FFD_CONTROL_POINT ( 7, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind, k_Ind, x_Mov, y_Mov, z_Mov )
    - FFD_TWIST_ANGLE ( 9, Scale | Mark. List | FFD_Box_ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
    - FFD_ROTATION ( 10, Scale | Mark. List | FFD_Box_ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
    - FFD_CAMBER ( 11, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind )
    - FFD_THICKNESS ( 12, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind )
    - FFD_CONTROL_POINT_2D ( 15, Scale | Mark. List | FFD_Box_ID, i_Ind, j_Ind, x_Mov, y_Mov )
    - FFD_CAMBER_2D ( 16, Scale | Mark. List | FFD_Box_ID, i_Ind )
    - FFD_THICKNESS_2D ( 17, Scale | Mark. List | FFD_Box_ID, i_Ind )
    
    """
  
  def __init__(self,*args,**kwarg):
    ordered_bunch.__init__(self)
    self.FFDTAG = []
    self.PARAM  = []
    self.update(ordered_bunch(*args,**kwarg))
  
  def append(self,new_dv):
    self.FFDTAG.append(new_dv['FFDTAG'])
    self.PARAM. append(new_dv['PARAM'])
  
  def extend(self,new_dvs):
    assert isinstance(new_dvs,DV_KIND) , 'input must be of type DV_KIND'
    self.FFDTAG.extend(new_dvs['FFDTAG'])
    self.PARAM. extend(new_dvs['PARAM'])

#: class DV_KIND
