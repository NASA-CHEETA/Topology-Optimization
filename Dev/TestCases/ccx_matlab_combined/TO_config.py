#!/usr/bin/env python3

from __future__ import print_function
import nlopt
import numpy as np
from topopt_functions import callrunccx, getVolconstraint
import readwriteData as rwd

#################################################################################################3

def TopOpt(restartFlag):

    ccxversion="caltop_gHVrc"
    inp= "cant3d"
    passivefilename = False #'passiveElemList.nam'   # elements excluded from design doamin, put a randon integer if empty


    dimension= 9000
    volfrac=0.15    
    rmin=3    
    p=3

    compliance_scaling = 1
    volume_scaling = 1

    callrunccx.beta = 0.1# Beta
    callrunccx.beta_factor = 4 # Beta multiplier
    callrunccx.betamax = 40
    callrunccx.beta_iter = 20 # Iteration for each beta update
    callrunccx.beta_iter_init = 0
    
    callrunccx.eta = 0.5
    callrunccx.betamaxcounter = 0

    callrunccx.nnz = 300 # Number of nonzeros in each row of density filter assumed
    callrunccx.counter=0    #Overall iteration counter
    
    callrunccx.q =3     #Filter order start
    callrunccx.q_max=3   #Max filter order
    callrunccx.q_iter = 30  #Iteration number to update the fitler order
    callrunccx.q_factor = 1.5 #Factor for filter order update

    topopt_maxIter = 200

    plotIterationWiseFlag = False # True if post process for each TO iteration 
    


    #########################################################################################
    #########################################################################################
    opt = nlopt.opt(nlopt.LD_MMA, dimension)

    outputCompliance="Obj_iterations.dat"

    lower=0.00001*np.ones(dimension)
    upper=0.99999*np.ones(dimension)

    if restartFlag:
        x0 = rwd.rho_column("density.dat")
    else:
        x0 =volfrac*np.ones(dimension)

    # Initialize passive elements to 1
    if isinstance(passivefilename, str):
        print('Setting passive elements')
        passiveEl = rwd.readElemList(passivefilename)
        x0 = rwd.setElemVal(x0, passiveEl,0.99999* np.ones(np.size(passiveEl))) 

    # DV Box bounds
    opt.set_lower_bounds(lower)
    opt.set_upper_bounds(upper)

    # DV Box bounds
    f=lambda x,grad: callrunccx(x,grad,ccxversion,inp,volfrac,rmin,p,outputCompliance, plotIterationWiseFlag, compliance_scaling,passivefilename)


    ff = open(outputCompliance, "w")
    ff.write("#Contains all compliance in each iteration")
    ff.close()


    opt.set_min_objective(f)

    # Volume constraint
    opt.add_inequality_constraint(lambda x, grad: getVolconstraint(x,grad,ccxversion,inp,volfrac,rmin,p, volume_scaling,passivefilename), 1e-8)
    #opt.add_inequality_constraint(lambda x, grad: myconstraint(x,grad, -1, 1), 1e-8)
    #opt.set_ftol_rel(1e-4)
    opt.set_xtol_rel(1e-3)
    opt.set_maxeval(topopt_maxIter)
    opt.set_param("inner_maxeval",10)
    x = opt.optimize(x0)
    minf = opt.last_optimum_value()
    #print('optimum at ', x)
    print('minimum value = ', minf)
    print('result code = ', opt.last_optimize_result())
    print('nevals = ', opt.get_numevals())
    #print('initial step =', opt.get_initial_step(x0))

    return 0


if __name__=='__main__':

    restartFlag = False
    topopt = TopOpt(restartFlag)
