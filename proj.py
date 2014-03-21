import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as mpl
import time
from scipy.optimize import curve_fit as fit

################
#1. Constants
################
lam = 1 #Wavelength: 1cm
arcmin = (1.0/60.0)*(np.pi/180.0) #1 arc minute in radians
arcsec = arcmin/60.0
im = 1j
pi = np.pi
meter = 100

################
#2. Load Data
################
vis = np.loadtxt("Visibilities.csv",delimiter=',') #visibilities ( (i,j,A,phi) )
pos = np.loadtxt("AntennaPositions.csv",delimiter=',') #positions ( (x,y) )

##########################
#3. Create derived data
##########################
n = pos.shape[0] #Store total number of antennae
nbas=n*(n-1)/2. #Calculate the number of independent baselines
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
uvvis=np.zeros((2*np.shape(vis)[0],np.shape(vis)[1]))
for i in range(np.shape(uvvis)[0]/2):
    a,b,amp,phi=vis[i]
#The first half of the array is for baselines (1,2), (1,3), ...
    uvvis[i][0]=u[a-1,b-1]
    uvvis[i][1]=v[a-1,b-1]
    uvvis[i][2]=amp
    uvvis[i][3]=phi
#The second half is for baselines (2,1), (3,1), ...
#The (u,v) values of these have opposite sign, as does the phase    
    uvvis[i+np.shape(uvvis)[0]/2]=-u[a-1,b-1]
    uvvis[i+np.shape(uvvis)[0]/2][1]=-v[a-1,b-1]
    uvvis[i+np.shape(uvvis)[0]/2][2]=amp
    uvvis[i+np.shape(uvvis)[0]/2][3]=-phi

#Angular Ranges (l,m)
lmax,mmax = 100*arcsec,100*arcsec #Span to +/- 100'' 
res = 1e2 #Resolution
l = np.linspace(-lmax,lmax,res) #RA from zenith
m = np.linspace(-mmax,mmax,res) #DEC from zenith

##########################
#4. Method definitions
##########################

#Gaussian function for Antenna Beam A(l,m)
def A(l,m,sig=arcmin):
    a = 1.0/(sig*np.sqrt(2*pi))
    b = (1.0/(2*sig**2))
    return a*np.exp( b*(- l**2 - m**2 ) )

#RHS function for Discrete FT: Sum over (u,v)
def DFT_rhs(uvvis,_l,_m):
    rhs = 0.0
    for (u,v,A,phi) in uvvis:
        rhs += A*np.exp(im*phi)*np.exp(2*pi*im*(u*_l + v*_m) )
	rhs += A*np.exp(-im*phi)*np.exp(-2*pi*im*(u*_l + v*_m) )
    return rhs
    
#Returns intensity array I(l,m) using DFT_rhs
def DFT(uvvis,L,M):
    I = np.zeros( (len(L),len(M)) ) + im*np.zeros( (len(L),len(M)) ) 
    for i,_l in enumerate(L):
        for j,_m in enumerate(M):
            a = np.sqrt( 1 - _l**2 - _m**2 )
            I[i,j] = (a/A(_l,_m))*DFT_rhs(uvvis,_l,_m)
    return I  

#Takes a measured (u,v) point and locates the nearest (u,v) point on an evenly spaced grid
def find_nearest_gridpoint(x,xlist):
	pos=bisect_left(xlist,x)
	if pos == 0:
        	return 0
        if pos == len(xlist):
        	return len(xlist)-1
        before = xlist[pos - 1]
        after = xlist[pos]
        if after - x < x - before:
        	return pos
        else:
        	return pos-1

#Create an evenly spaced grid of visibilities to use in an inverse fft
def uv_grid(uvvis,ugrid,vgrid):
    #visibilities will be complex
	gridded_visibilities = np.zeros((len(ugrid),len(vgrid)))+im*np.zeros((len(ugrid),len(vgrid)))
	for i in range(np.shape(uvvis)[0]):
		#Take the measurements,
		umeas=uvvis[i][0]
		vmeas=uvvis[i][1]
		amp=uvvis[i][2]
		phi=uvvis[i][3]
		#Find the nearest gridpoint,
		j=find_nearest_gridpoint(umeas,ugrid)
		k=find_nearest_gridpoint(vmeas,vgrid)
		#and add the visibility to that gridpoint:
		gridded_visibilities[j][k]+=amp*np.exp(im*phi)
	return gridded_visibilities

    
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

    antennae_int = np.array( [ int(antennae[i]) for i in range(N) ] ) #Cast to ints
    #Return the indices of these rows 
    return rows,antennae_int



##############################################
#MAIN 1: DFT
##############################################
def part_one(N,order,l,m,lmax,mmax):

    if order=='asc': imstring="closest"
    elif order=='desc': imstring="farthest"
    else:
        print "Error: order must be 'asc' or 'desc'"  

    #Get uvvis rows and indices of chosen antennae
    rows,antennae = get_selection(pos,vis,N,order) 

    #Plot positions of antennae being used
    mpl.figure()
    mpl.plot( pos[:,1]/meter, pos[:,2]/meter, 'k.',label="inactive antennae")
    mpl.plot(pos[antennae-1,1]/meter,pos[antennae-1,2]/meter,'ro',label="active antennae")
    mpl.xlabel('$x [m]$',fontsize=20)
    mpl.ylabel('$y [m]$',fontsize=20)
    mpl.legend()
    mpl.savefig("%i_%santennae.pdf" % (N,imstring))

    #Create cropped selection of uvvis using only the selected antennae
    uvvis_cropped = np.zeros( (len(rows),vis.shape[1]) )
    for i in range(len(rows)): uvvis_cropped[i] = uvvis[rows[i]]

    #Get intensity using cropped array
    #t0 = time.clock()
    dftI = DFT(uvvis_cropped,l,m)
    #t = time.clock()-t0 

    #Plot results
    lmax=lmax/arcsec
    mmax=mmax/arcsec
    mpl.figure()
    mpl.imshow(np.abs(dftI),extent=[-lmax,lmax,-mmax,mmax])
    mpl.xlabel('$l (arcseconds)$',fontsize=20)
    mpl.ylabel('$m (arcseconds)$',fontsize=20)
    mpl.savefig("DFT_image_%i%s.pdf" % (N,imstring))
    mpl.show()
    
    return t

#TIMING: UNCOMMENT TO USE
#def quad(x,a,b,c): return a*x**2 + b*x + c #define a linear function for fitting
#NN = np.array( [2,4,6,8,10,20] ) #Values of N to time
#TT = np.zeros(len(NN))
#for i,N in enumerate(NN):
#    TT[i] = part_one(N,'asc',l,m,lmax,mmax) #Store time taken to run
#    print N,TT[i]   
#a,b,c = fit( quad, NN, TT )[0] #Fit quadratic to data
#NN_sm = np.linspace(NN[0],NN[-1],100) #Create smoothed domain
#mpl.figure()
#mpl.plot(NN,TT,'ko')
#mpl.plot(NN_sm,quad(NN_sm,a,b,c),'r-',label=r"$%.2fx^{2} + %.2fx + %.2f$" % (a,b,c))
#mpl.xlabel("$t$")
#mpl.ylabel("$N$")
#mpl.legend()
#mpl.savefig("DFT_image_timing.pdf")
#mpl.show()


#############################################
#MAIN 2: FFT
#############################################

#Create a grid of u and v values to fill
ugrid = np.linspace(-60000,60000,res)
vgrid = np.linspace(-60000,60000,res)

#Fill the grid
Vgrid=uv_grid(uvvis,ugrid,vgrid)

#Calculate intensity (which could have a small imaginary part due to numerical error)
fftI=np.zeros(np.shape(Vgrid))+np.zeros(np.shape(Vgrid))*im
for i in range(len(l)):
	for j in range(len(m)):
		fftI[i][j]=(np.fft.fftshift(np.fft.ifft2(Vgrid))[i][j]*(1-l[i]**2-m[j]**2)**0.5)/(A(l[i],m[j]))
		

#Plot the results (use absolute value of intensities just in case there is a small imagniary part)
mpl.imshow(np.abs(fftI),extent=[-lmax,lmax,-mmax,mmax])
mpl.xlabel('$l$ (arcseconds)',fontsize=20)
mpl.ylabel('$m$ (arcseconds)',fontsize=20)
mpl.savefig('FFT_Image.pdf')
mpl.show()


