import numpy as np

#1. Constants
lam = 1 #Wavelength: 1cm
arcmin = (1.0/60.0)*(np.pi/180.0) #1 arc minute in radians

#2. Load Data
vis = np.loadtxt("Visibilities.csv",delimiter=',') #visibilities ( (i,j,A,phi) )
pos = np.loadtxt("AntennaPositions.csv",delimiter=',') #positions ( (x,y) )

#3. Create derived data
n = len(pos) #store total number of antennae
B = np.zeros( (n,n) ) #Will store baseline vector (u,v) at each index (i,j)
for i in range(n):
    for j in range(n):       
        #Get antennae positions
        xi,yi = pos[i]
        xj,yj = pos[j]
        #Get wavelength-calibrated baseline vector
        u = (xi - xj)/lam
        v = (yi - yj)/lam
        #Store in array
        B[i,j] = (u,v)

#4. Methods

def gauss(x,y,sig=arcmin):
    return 
    
