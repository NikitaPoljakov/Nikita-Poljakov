import scipy as sp
import matplotlib
import matplotlib.pyplot as py
import math
import numpy as np

#Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'

# [Starter code for Problem 8]

# ------------------------------------------------------------
# Solving the bulk with an inhomogeneous iterative method
# ------------------------------------------------------------

#Define some useful objects
class Latticesite:
    def __init__(self, siteNumber, epsilon_w, pot):
        self.index = siteNumber
        self.coordinate = []
        self.NNs = []
        self.NNNs = []
        
        if pot == True:
            self.potential = potential(np.floor(siteNumber/Ly), epsilon_w)
        else:
            self.potential = 0.0
            
        self.density_current = 0.0
        self.density_previous = 0.0
    
    def update(self):
        '''Update the density. Remember to account for divergences under iteration'''
        if self.density_current<1.0:
            self.density_previous = self.density_current
        elif self.density_current>1.0:
            self.density_previous = 1.0
        else:
            self.density_previous = 0.0

def iterate(sites, k, mu, beta, epsilon):
    '''Perform a single iteration of eq. 46 for the particle at site k, given the sitelist sites'''
    nDens = 0.0
    nnDens = 0.0
    for neighbor in sites[k].NNs:
        nDens = nDens + sites[neighbor].density_previous
    for neighbor in sites[k].NNNs:
        nnDens = nnDens + sites[neighbor].density_previous
    return (1-sites[k].density_previous)*sp.exp(beta*(mu + epsilon * nDens + 0.25 * epsilon * nnDens - sites[k].potential)) #Here we assume T,mu,V are all in units of the interaction strength 

def potential(y, epsilon_w):
    
    if y >= 1:
        return - epsilon_w * y**(-3)
    
    else:
        return 10**10

#The lattice is a (Lx x Ly square lattice)
Lx = 10
Ly = 10

#Initialize the lattice
def initialize(Lx, Ly, epsilon_w, potential):

    def getIndex(Lx,nSites):
        '''Converts from coordinates (x,y) to lattice index'''
        return lambda x,y:(Lx*x + y)

    nSites = Lx*Ly
    sites = [Latticesite(k, epsilon_w, potential) for k in range(nSites)]  
    pos = getIndex(Lx,Ly)      
    for site in sites:
        x,y = [site.index//Lx,site.index%Lx]                                #Get x and y coordinates and from those the coordinates of the neighbors
        nns = set([((x+1)%Lx,y),((x-1)%Lx,y),(x,(y+1)%Ly),(x,(y-1)%Ly)])    
        nnns = set([( (x+1)%Lx,(y+1)%Ly ),( (x-1)%Lx, (y-1)%Ly),((x-1)%Lx,(y+1)%Ly),((x+1)%Lx,(y-1)%Ly)])
        site.NNs = [pos(x[0],x[1]) for x in nns]           #Store the neighbor indices as instance variables
        site.NNNs = [pos(x[0],x[1]) for x in nnns]
        site.density_previous = 0.2                        #Initialize the system in the low density limit
    return sites

#Now we iterate the solver until the density is converged
def run(mu, T, epsilon, epsilon_w, Lx,Ly, potential, cTol=10**-8, mixing=0.1, iterMax=1000, show=True):
    'Calculates the density profile at a given mu,T'
    sites = initialize(Lx,Ly, epsilon_w, potential)
    convergence = 0.1
    iteration = 0.0
    while (convergence>cTol) and (iteration<iterMax):
        for k in range(len(sites)):
            sites[k].density_current = sites[k].density_previous*(1 - mixing) + mixing*iterate(sites,k,mu,1/T, epsilon) #Calculate new state of the system from the old state of the system
        two_norm = sum([(site.density_current-site.density_previous)**2 for site in sites])
        convergence = math.sqrt(two_norm)
        iteration = iteration + 1
        for site in sites:
            site.update()

    'Can then return an image of the density profile'
    z = []
    for site in sites:
        z.append(site.density_previous)
    Z = sp.array(z)
    Z = sp.reshape(Z,(Lx,Ly))
    Z = np.transpose(Z)    
    
    return Lx * Z[0]

#Run a few examples
epsilon = 1
epsilon_w = [1.33, 1.42, 1.5]
beta = 1.2
y1 = []
y2 = []
y3 = []

x = np.linspace(-0.5, 0, 50)

for i in x:
    fluidDensities = run((i/beta) - 2.5 * epsilon, 1 / beta, epsilon, epsilon_w[0], Lx, Ly, potential = True)
    bulkGasDensities = run((i/beta) - 2.5 * epsilon, 1 / beta, epsilon, epsilon_w[0], Lx, Ly, potential = False)
    y1.append(np.sum(fluidDensities - bulkGasDensities))
    
    fluidDensities = run((i/beta) - 2.5 * epsilon, 1 / beta, epsilon, epsilon_w[1], Lx, Ly, potential = True)
    bulkGasDensities = run((i/beta) - 2.5 * epsilon, 1 / beta, epsilon, epsilon_w[1], Lx, Ly, potential = False)
    y2.append(np.sum(fluidDensities - bulkGasDensities))
    
    fluidDensities = run((i/beta) - 2.5 * epsilon, 1 / beta, epsilon, epsilon_w[2], Lx, Ly, potential = True)
    bulkGasDensities = run((i/beta) - 2.5 * epsilon, 1 / beta, epsilon, epsilon_w[2], Lx, Ly, potential = False)
    y3.append(np.sum(fluidDensities - bulkGasDensities))

py.style.use("dark_background")
py.plot(x, y3, color = 'b', label = "$βƐ_w$ = 1.8")
py.plot(x, y2, color = 'c', label = "$βƐ_w$ = 1.7")
py.plot(x, y1, color = 'm', label = "$βƐ_w$ = 1.6")
py.xlabel("$β$($μ$ - $μ_{coex}$)", size = 20)
py.ylabel('$Γ$', size = 20)
py.legend(prop={'size': 15})
py.tight_layout()
py.savefig('Problem12.pdf')
py.show()