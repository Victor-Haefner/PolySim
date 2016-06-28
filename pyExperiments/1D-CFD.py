# source: http://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/01_Step_1.ipynb
# Prof. Barba, 12 steps to Navier-Stokes

import numpy                       
from matplotlib import pyplot      
import time, sys

nx = 41 		# system size
dx = 2.0/(nx-1) # system scale
nt = 25    		# timesteps simulated
dt = .025  		# timestep
c = 1      		# wavespeed

# setup system
u = numpy.ones(nx)
u[.5/dx : 1.0/dx+1]=2
pyplot.plot(numpy.linspace(0,2,nx), u); # draw initial state

for n in range(nt): # for each timestep
	un = u.copy()
	for i in range(1,nx): u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1]) # compute new state
	pyplot.plot(numpy.linspace(0,2,nx), u); # draw current state
	
pyplot.show()
