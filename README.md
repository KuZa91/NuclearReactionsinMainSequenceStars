# NuclearReactionsinMainSequenceStars
Software for the simulation of the evolution of the light elements abundances due to nuclear reactions regarding pp and cno cycle in main sequence stars (Sun Like).
It's possible to simulate the evolution in abundance of main sequence star elements untill a end time tend has been reached or untill the H abundance go under a certain fixed value.
The abundances of pp and cno cycle will be automatically graphied using gnuplot at the end of the simulation plus the program will create a abundance.txt file in which there will be all abundances in function of time.
The evolution is simulated using Explicit Euler method, for fiducial values try :

T = 1.7e7
rho = 1.6e2
x(4He) =0.28
x(3h3) = 0.
Z = 0.02
X(12C) = 3.4714e-3
x(14N) = 1.0652e-3
x(16O) = 9.6702e-3
