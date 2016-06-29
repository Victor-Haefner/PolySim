# source: http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/
# Prof. Barba, 12 steps to Navier-Stokes

import numpy                       
from matplotlib import pyplot      
import time, sys

nx = 161 		# system size
dx = 2.0/(nx-1) # system scale
nt = 125    		# timesteps simulated
dt = .0001  		# timestep
c = 40      		# wavespeed
vis = 0.1		# viscuosity

def initWave(): # setup system
	u = numpy.ones(nx)
	u[.9/dx : 1.0/dx+1]=2
	pyplot.plot(numpy.linspace(0,2,nx), u); # draw initial state
	return u

def propagate(u, nonlinear):
	for n in range(nt): # for each timestep
		un = u.copy()
		for i in range(1,nx-1):
			diffusion = vis*dt/dx/dx*(un[i+1]-2*un[i]+un[i-1])
			convection = c*dt/dx*(un[i]-un[i-1])
			if nonlinear: convection *= un[i]
			u[i] = un[i] - convection + diffusion # compute new state, based on Burgers equation
		pyplot.plot(numpy.linspace(0,2,nx), u); # draw current state
		

u = initWave()
propagate(u, False) # linear

u = initWave()
propagate(u, True) # non linear


pyplot.show()
