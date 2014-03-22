#######################
# PROJECT DESCRIPTION
#######################

This 'Ay190: Computational Astrophysics' project takes visibilities from a radio interferometry array and constructs an image of the source using both a Discrete Fourier Transform method and NumPy's Fast Fourier Transform method.

#############
#INPUT FILES
#############
AntennaPositions.csv - This is a text file that contains the position (of each antenna) measured in centimeters. This first column gives the antenna number, the second column gives the x-position of the antenna, and the final column gives the y-position of the antenna.

Visibilities.csv - This file gives the visibilities measured by the radio interferometer. These visibilities have no thermal noise and are perfectly phase-calibrated (but they are not flux-calibrated so the amplitude is meaningless). The first two columns define the antennas that the visibility was measured between (for example, 2 and 17 indicates that the visibility was measured between antennas 2 and 17). The final two columns give the amplitude and phase of the visibility respectively (the complex visibility is then amplitude*exp(i*phase) ).

These observations have been made at a wavelength of 1cm (with an infinitesimally small bandpass filter so that there is no bandwidth smearing). Additionally, the source is located directly overhead at the zenith (so there are no projection effects), and the beam of each antenna is Gaussian with sigma=1 arcminute. The integration time is too small for the rotation of the Earth to have any impact on the source's position in the sky.

##################
#DATA STRUCTURES
##################

#'pos' - Antenna Positions
    pos( (i,x,y) )
        i is the antenna number, x and y are cartesian coordinates
        -->i,x,y = pos(0) would unpack the number and position of the first antenna

#'vis' - Positions and Visibilities       
    vis( (i,j,A,phi) )
        i and j are indices of the two antennae in 'pos', 
        A is amplitude,
        phi is phase
    
#'uvvis' - Baselines and Visibilities
    uvvis( (u,v,A,phi) )
        for each row in 'vis', 'uvvis' stores:
            u,v - the baseline vectors (x2-x1) and (y2-y1) of antennae i and j
            A,phi - the same visibility as in 'vis'

#'l' and 'm'
    These 1D arrays are grids of RA and DEC ranging from -100'' to +100''
    Centered on zenith


###################
#METHODS OVERVIEW (more detailed comments in code and description in writeup)
###################

#Gaussian function for Antenna Beam A(l,m)
def A(l,m,sig=arcmin)

#RHS function for Discrete FT: Sum over (u,v)
def DFT_rhs(uvvis,_l,_m)
 
#Returns intensity array I(l,m) using DFT_rhs
def DFT(uvvis,L,M)

#Takes a measured (u,v) point and locates the nearest (u,v) point on an evenly spaced grid
def find_nearest_gridpoint(x,xlist)

#Create an evenly spaced grid of visibilities to use in an inverse fft
def uv_grid(uvvis,ugrid,vgrid)
 
#Return indices of rows in uvvis which only use N closest/farthest antennae
def get_selection(pos,vis,N,orderby='asc')





