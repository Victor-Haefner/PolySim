# source: http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/
# Prof. Barba, 12 steps to Navier-Stokes

import numpy   
from mpl_toolkits.mplot3d import Axes3D                    
from matplotlib import cm, pyplot      
import time, sys

nx = 41 		# system size
ny = 41 		# system size
dx = 2.0/(nx-1) # system scale
dy = 2.0/(ny-1) # system scale
nt = 200   		# timesteps simulated
nt2 = 50		# substeps
dt = 0.001  		# timestep
c = 1      		# wavespeed
nu = 0.1		# viscuosity
rho = 0.1
F = 0

# plot
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)                    
X, Y = numpy.meshgrid(x,y)    

def plot(u,v,p):
	fig = pyplot.figure(figsize=(11,7), dpi=100)          ##the figsize parameter can be used to produce different sized images
	pyplot.contourf(X,Y,p,alpha=0.5)    ###plnttong the pressure field as a contour
	pyplot.colorbar()
	#pyplot.contour(X,Y,p)               ###plotting the pressure field outlines
	pyplot.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
	pyplot.xlabel('X')
	pyplot.ylabel('Y')

def initGrid(): # setup system
	u = numpy.zeros((nx,ny))
	v = numpy.zeros((nx,ny))
	p = numpy.zeros((nx,ny))
	#u[.5/dy:1.0/dy+1,.5/dx:1.0/dx+1]=2
	#v[.5/dy:1.0/dy+1,.5/dx:1.0/dx+1]=2
	#p[.5/dy:1.0/dy+1,.5/dx:1.0/dx+1]=2
	return u,v,p

def buildUpB(b, rho, dt, u, v, dx, dy):
    
    b[1:-1,1:-1]=rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))-\
                      ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2-\
                      2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))-\
                      ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)

    return b
      
def pressure(p, dx, dy, b, rho, dt, u, v):
	pn = numpy.empty_like(p)
	pn = p.copy()

	b = buildUpB(b, rho, dt, u, v, dx, dy)

	for q in range(nt2):
		pn = p.copy()
		p[1:-1,1:-1] = ((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2+(pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)/\
						(2*(dx**2+dy**2)) -\
						dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]

		p[:,-1] =p[:,-2] ##dp/dy = 0 at x = 2
		p[0,:] = p[1,:]  ##dp/dy = 0 at y = 0
		p[:,0]=p[:,1]    ##dp/dx = 0 at x = 0
		p[-1,:]=0        ##p = 0 at y = 2
		
	return p, b
    
def pressure2(p, dx, dy, b, rho, dt, un, vn):
	row, col = u.shape
	dx2 = dx**2
	dy2 = dy**2
	K = 2*(dx2+dy2)
	P2 = rho*dx2*dy2/K
	
	b = numpy.zeros((ny, nx))
	for j in range(1, row-1):
		for i in range(1, col-1):
			
			Dui = un[j,i+1]-un[j,i-1]
			Dvj = vn[j+1,i]-vn[j-1,i]
			
			P3 = (Dui/dx + Dvj/dy)/dt/2
			P4 = Dui*Dui/(4*dx2)
			P5 = (un[j+1,i]-un[j-1,i])*(vn[j,i+1]-vn[j,i-1])/(2*dx*dy)
			P6 = Dvj*Dvj/(4*dy2)
			
			b[j,i] = P3-P4-P5-P6
	
	for q in range(nt2): # poisson 
		pn = p.copy()
		
		for j in range(1, row-1):
			for i in range(1, col-1):				
				P1 = ((pn[j,i+1]+pn[j,i-1])*dy2 + (pn[j+1,i]+pn[j-1,i])*dx2)/K
				p[j,i] = P1 - P2*b[j,i]
	
		p[:,-1] =p[:,-2] ##dp/dy = 0 at x = 2
		p[0,:] = p[1,:]  ##dp/dy = 0 at y = 0
		p[:,0]=p[:,1]    ##dp/dx = 0 at x = 0
		p[-1,:]=0        ##p = 0 at y = 2
				
	return p, None
    
def flow(un, vn, rho, nu, dt):
	u[1:-1,1:-1] = un[1:-1,1:-1]-\
					un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
					vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
					dt/(2*rho*dx)*(p[1:-1,2:]-p[1:-1,0:-2])+\
					nu*(dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+\
					dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])) + dt*F

	v[1:-1,1:-1] = vn[1:-1,1:-1]-\
					un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
					vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
					dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])+\
					nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+\
					(dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])))
					
	
	u[0,:] = 0
	u[:,0] = 0
	u[:,-1] = 0
	u[-1,:] = 1    #set velocity on cavity lid equal to 1
	v[0,:] = 0
	v[-1,:]=0
	v[:,0] = 0
	v[:,-1] = 0
					
	return u,v
        
def flow2(un, vn, rho, nu, dt):
	row, col = u.shape
	for j in range(1, row-1):
		for i in range(1, col-1):
			
			u_convection = un[j,i]*dt*(un[j,i] - un[j,i-1])/dx + vn[j,i]*dt*(un[j,i]-un[j-1,i])/dy
			v_convection = un[j,i]*dt*(vn[j,i] - vn[j,i-1])/dx + vn[j,i]*dt*(vn[j,i]-vn[j-1,i])/dy
			u_diffusion = nu*dt*( (un[j,i+1]-2*un[j,i]+un[j,i-1])/dx/dx + (un[j+1,i]-2*un[j,i]+un[j-1,i])/dy/dy )
			v_diffusion = nu*dt*( (vn[j,i+1]-2*vn[j,i]+vn[j,i-1])/dx/dx + (vn[j+1,i]-2*vn[j,i]+vn[j-1,i])/dy/dy )
			u_pressure = dt/(2*rho*dx)*(p[j,i+1]-p[j,i-1])
			v_pressure = dt/(2*rho*dy)*(p[j+1,i]-p[j-1,i])
			u_source = dt*F
			v_source = 0
			
			u[j,i] = un[j,i] - u_convection - u_pressure + u_diffusion + u_source
			v[j,i] = vn[j,i] - v_convection - v_pressure + v_diffusion + v_source
			

	u[0,:] = 0
	u[:,0] = 0
	u[:,-1] = 0
	u[-1,:] = 1    #set velocity on cavity lid equal to 1
	v[0,:] = 0
	v[-1,:]=0
	v[:,0] = 0
	v[:,-1] = 0
			
	return u,v
	
def compute(nt, u, v, dt, dx, dy, p, rho, nu):    
	un = numpy.empty_like(u)
	vn = numpy.empty_like(v)
	
	b = numpy.zeros((ny, nx))

	for n in range(nt):
		print 'compute', n
		
		un = u.copy()
		vn = v.copy()

		p, b = pressure(p, dx, dy, b, rho, dt, un, vn)
		u,v = flow(un, vn, rho, nu, dt)

	return u, v, p

u,v,p = initGrid()
u,v,p = compute(nt, u, v, dt, dx, dy, p, rho, nu)

plot(u,v,p) # draw current state
pyplot.show()
