from __future__ import print_function
import numpy as np
import readwriteData as rwd

def getVolconstraint(x, grad,ccxversion,inp,volfrac,rmin,penalty):
	volConstraint=rwd.getVolconstr()
	scaling=1
	if grad.size > 0:
		a=rwd.getCurrentVolsens()
		grad[:]=a*scaling
	return volConstraint*scaling
