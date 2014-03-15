import numpy as np
import matplotlib.pyplot as mpl

#1. Constants
lam = 1 #Wavelength: 1cm
arcmin = (1.0/60.0)*(np.pi/180.0) #1 arc minute in radians
arcsec = arcmin/60.0
im = 1j
pi = np.pi

#2. Load Data
vis = np.loadtxt("Visibilities.csv",delimiter=',') #visibilities ( (i,j,A,phi) )
pos = np.loadtxt("AntennaPositions.csv",delimiter=',') #positions ( (x,y) )

#3. Create derived data
n = pos.shape[0] #store total number of antennae
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
        
        
#Create array of form (u,v,A,phi)
uvvis=np.zeros(np.shape(vis))
for i in range(np.shape(uvvis)[0]):
    a,b,A,phi=vis[i]
    uvvis[i][0]=u[a-1,b-1]
    uvvis[i][1]=v[a-1,b-1]
    uvvis[i][2]=A
    uvvis[i][3]=phi

print uvvis.shape

#3.2 Angular Ranges (l,m)
lmax,mmax = 100*arcsec,100*arcsec #Span to +/- 100'' 
res = 1e2 #Resolution
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
    return rhs
    
#Returns intensity array I(l,m) using DFT_rhs
def DFT(uvvis,L,M):
    I = np.zeros( (len(L),len(M)) )
    for l in L:
        for m in M:
            a = np.sqrt( 1 - l**2 - m**2 )
            I[l,m] = a*DFT_rhs(uvvis,l,m)/A(l,m)
    return I      

#Return indices of rows in uvvis which only use N closest/farthest antennae
def get_selection(pos,vis,N,orderby='asc'):

    #Create array of (i,dist from origin)
    r = np.zeros( (pos.shape[0],2) )
    r[:,0] = pos[:,0]
    r[:,1] = np.sqrt( pos[:,1]**2 + pos[:,2]**2 )
    
    #Sort by distance
    sorted_indices = r[:,1].argsort()
    r_sorted = r[sorted_indices]
    
    #Extract indices of N closest/farthest antennae
    antennae = np.zeros(N)
    if orderby=='asc': antennae = r_sorted[:N][:,0]
    elif orderby=='desc': antennae = r_sorted[-N:][:,0]
    else:
        print "Error in parameters: orderby must be either 'asc' or 'desc'"
        return []
        
    #Create dictionary to quickly reference if an antenna is one of N chosen    
    ant_dic = {}
    for x in antennae: ant_dic[x]=True
    
    #Run through 'vis' array and pick rows where both antenna are in ant_dic
    rows = []
    for ind,(i,j,A,phi) in enumerate(vis):
        if ant_dic.has_key(i) and ant_dic.has_key(j):
            rows.append(ind)

    #Return the indices of these rows 
    return rows
    
#Input parameters  
order = 'asc'
N=10

#Get rows to use from uvvis
rows = get_selection(pos,vis,N,order) 
L = len(rows)

#Create cropped selection of uvvis
uvvis_cropped = np.zeros( (L,vis.shape[1]) )
for i in range(L):
    uvvis_cropped[i] = uvvis[rows[i]]

I = DFT(uvvis_cropped,l,m) #Get intensity array using UVIS    

mpl.figure()
mpl.pcolor(I)
mpl.show()
    
