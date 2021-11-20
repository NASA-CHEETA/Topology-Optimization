import RunMDA

Nnodes = 7946
Nelems = 31005

E   =  2.5E5
nu  =  0.35
rho =  100

n = 4

for n in range(1,4):
    Flag = RunMDA.MDA(Nnodes, Nelems, E, nu, rho, n)