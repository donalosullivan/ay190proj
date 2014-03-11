Naming Convetions:

raw data:
    Antenna Positions: load into 2d numpy array called 'pos': 
        pos( (x,y) )
            such that x,y = pos(0) would unpack the position of the first antenna
            
    Visibilities: load into 2d numpy array called 'vis': 
        vis( (i,j,A,phi) )
            i and j are indices of antennae in 'pos', 
            A is amplitude,
            phi is phase
        
        indexes are such that vis[0] gives first row in text file

lalalalalalalalala
