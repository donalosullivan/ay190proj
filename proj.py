import numpy as np

#1. Constants
lam = 1 #Wavelength: 1cm
arcmin = (1.0/60.0)*(np.pi/180.0) #1 arc minute in radians

#2. Load Data
vis = np.loadtxt("Visibilities.csv",delimiter=',') #visibilities ( (i,j,A,phi) )
pos = np.loadtxt("AntennaPositions.csv",delimiter=',') #positions ( (x,y) )

#3. Create derived data
n = np.shape(pos)[0] #store total number of antennae
nbas=n*(n-1)/2. #calculate the number of independent baselines
u=np.zeros((n,n))
v=np.zeros((n,n))
for i in range(n):
    for j in range(n):       
        #Get antennae positions
        xi,yi = pos[i][1],pos[i][2]
        xj,yj = pos[j][1],pos[j][2]
        #Get wavelength-calibrated baseline vector
        u[i][j] = (xi - xj)/lam
        v[i][j] = (yi - yj)/lam

uvvis=np.zeros(np.shape(vis))
for i in range(np.shape(uvvis)[0]):
    a,b,A,phi=vis[i]
    uvvis[i][0]=u[a-1,b-1]
    uvvis[i][1]=v[a-1,b-1]
    uvvis[i][2]=A
    uvvis[i][3]=phi

#4. Methods

def gauss(x,y,sig=arcmin):
    return (1.0/sig)*(2*np.pi)**(-0.5)*np.exp( (1.0/(2*
    
