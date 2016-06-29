# source: http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/
# Prof. Barba, 12 steps to Navier-Stokes

import numpy   
from mpl_toolkits.mplot3d import Axes3D                    
from matplotlib import pyplot      
import time, sys

nx = 81 		# system size
ny = 81 		# system size
dx = 2.0/(nx-1) # system scale
dy = 2.0/(ny-1) # system scale
nt = 125    		# timesteps simulated
sigma = 0.2
dt = sigma*dx  		# timestep
c = 1      		# wavespeed
vis = 0.01		# viscuosity

# plot
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

def plot(u):
	fig = pyplot.figure(figsize=(11,7), dpi=100)          ##the figsize parameter can be used to produce different sized images
	ax = fig.gca(projection='3d')                      
	X, Y = numpy.meshgrid(x,y)                            
	surf = ax.plot_surface(X,Y,u[:])

def initWave(): # setup system
	u = numpy.ones((nx,ny))
	u[.5/dy:1.0/dy+1,.5/dx:1.0/dx+1]=2
	#plot(u) # draw initial state
	return u

def propagate(u, nonlinear):
	for n in range(nt): # for each timestep
		un = u.copy()
		
		
		row, col = u.shape
		for j in range(1, row-1):
			for i in range(1, col-1):
		
				#diffusion = vis*dt/dx/dx*(un[i+1]-2*un[i]+un[i-1])
				#convection = c*dt/dx*(un[i]-un[i-1])
				#u[i] = un[i] - convection + diffusion # compute new state, based on Burgers equation
				
				convection = c*dt/dx*(un[j,i] - un[j,i-1]) + c*dt/dy*(un[j,i]-un[j-1,i])
				diffusion = vis*dt/dx/dx*(un[j,i+1]-2*un[j,i]+un[j,i-1]) + vis*dt/dy/dy*(un[j+1,i]-2*un[j,i]+un[j-1,i])
				if nonlinear: convection *= un[i]
				if nonlinear: diffusion *= un[i]
				
				u[j,i] = un[j,i] - convection + diffusion
            
			
		

u = initWave()
propagate(u, False) # linear
plot(u) # draw current state
pyplot.show()
