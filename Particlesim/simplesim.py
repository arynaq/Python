import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import time

N = 2												#Number of particles
R = 400												#Box width
pD= R/100											#Particle diameter
mu= 9*(10**-4)										#Viscosity of medium
alpha_chance = 0.998								#Chance to radiate from core

r = np.random.randint(0, R, (N,2))					#Position vector
v = np.zeros((N,2))									#Velocity vector
a = np.array([0,0])		  	  						#Forces
cM = np.sum(r,axis=0)/N								#Center of mass


#Set up animation frames
plt.ion()
Nline,  = plt.plot([],'o')
NTline, = plt.plot([],'o')
cMline, = plt.plot([],'o', markersize=10, alpha="0.7")							
plt.axis([0,R,0,R])


while True:
	#Update plot
	Nline.set_ydata(r[::2,1])
	Nline.set_xdata(r[::2,0])
	NTline.set_ydata(r[1::2,1])
	NTline.set_xdata(r[1::2,0])
	cMline.set_xdata(cM[0])
	cMline.set_ydata(cM[1])
	plt.draw()

	#Collision tests
	r_hit_x0 = np.where(r[:,0]<0)					#Hit floor?
	r_hit_x1 = np.where(r[:,0]>R)					#Hit roof?
	r_hit_LR = np.where(r[:,1]<0)					#Left wall?
	r_hit_RR = np.where(r[:,1]>R)					#Right wall?
	

	#Stop at walls
	r[r_hit_x0,0] = 0
	r[r_hit_x1,0] = R
	r[r_hit_LR,1] = 0
	r[r_hit_RR,1] = R
	
	#Reverse velocities
	v[r_hit_x0,0] = -v[r_hit_x0,0]
	v[r_hit_x1,0] = -v[r_hit_x1,0]
	v[r_hit_LR,1] = -v[r_hit_LR,1]
	v[r_hit_RR,1] = -v[r_hit_RR,1]

	
	#Distance checks
	D = squareform(pdist(r))
	ind1, ind2 	= np.where(D < pD)
	c_1,c_2  	= np.where(D < 30*R)
	unique 		= (ind1 < ind2)
	c_unique 	= (c_1<c_2)
	
	ind1 = ind1[unique]
	ind2 = ind2[unique]
	c_1 = c_1[c_unique]
	c_2 = c_2[c_unique]
		
	
	#Columb
	for i1,i2 in zip(c_1,c_2):
		diff = r[i1,:] - r[i2,:]
		lngth= np.sqrt(np.sum(diff**2))
		sign = np.array([1,1])**[i1%2,i2%2]
		
		if lngth<pD:
			lngth = pD
			
		v[i1,:] = v[i1,:]  - sign[0]*5*(diff/lngth)
		v[i2,:] = v[i2,:]  - sign[1]*5*(-diff/lngth)
		
	#Collisions and radiation
	for i1, i2 in zip(ind1, ind2):
		v[i1,:]*= -1
		v[i2,:]*= -1
		
		
	
	#Physics
	friction = 0 #-3*np.pi*mu*v
	cM = np.sum(r,axis=0)/N
	v += 0.01*(a+friction)
	r += 0.01*v
	
		
		

	
	

