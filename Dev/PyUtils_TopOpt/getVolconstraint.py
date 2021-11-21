from __future__ import print_function
import numpy as np
import readwriteData as rwd

def getVolconstraint(x, grad,ccxversion,inp,volfrac,rmin,penalty, volume_scaling):
	volConstraint=rwd.getVolconstr() # g = v-vf*vmax
	total_volume = rwd.getTotalVol() # maximum volume of solid domain, vmax
	#scaling=1

	scaling = volume_scaling/(volfrac*total_volume) # to normalize as g : V/(VF*Vmax)-1<=0
	if grad.size > 0:
		a=rwd.getCurrentVolsens()
		grad[:]=a*scaling
	return volConstraint*scaling
