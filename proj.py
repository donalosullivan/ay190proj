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

#3.2 Angular Ranges (l,m)
lmax,mmax = 100*arcsec,100*arcsec #Span to +/- 100'' 
res = 1e3 #Resolution
l = np.linspace(-lmax,lmax,res) #RA from zenith
m = np.linspace(-mmax,mmax,res) #DEC from zenith

#4. Method definitions

#Gaussian function for Antenna Beam A(l,m)
def A(l,m,sig=arcmin):
    a = 1.0/(sig*np.sqrt(2*pi))
    b = (1.0/(2*sig**2))
    return a*np.exp( b*( l**2 - m**2 ) )

#RHS function for Discrete FT: Sum over (u,v)
def DFT_rhs(uvvis,l,m):
    rhs = 0
    for (u,v,A,phi) in uvvis:
        V = A*np.exp(im*phi)
        rhs += V*np.exp(-2*pi*im*(u*l + v*m) )

#Returns intensity array I(l,m) using DFT_rhs
def DFT(uvvis,L,M):
    I = np.zeros( (len(L),len(M)) )
    for l in L:
        for m in M:
            a = np.sqrt( 1 - l**2 - m**2 )
            I[l,m] = ( a/A(l,m) )*DFT_rhs(uvvis,l,m)
    return I      
    
    
