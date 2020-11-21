"""
Finding density when the upper boundary for the Gibbs-Bogoliubov inequality is minimal
"""
import matplotlib.pyplot as plt
import numpy as np

rho = np.linspace(0, 1, 200) # Densities
y1, y2, y3 = np.zeros(len(rho)), np.zeros(len(rho)), np.zeros(len(rho))
xAxis = np.zeros(len(rho))

# Part (i): beta * epsilon = 1 with µ/epsilon = −3.0

def graph(beta, epsilon, mu):

    for i in range(len(rho)):
        y1[i] = rho[i] - (1 - rho[i]) * np.exp(beta*(mu[0] + 5 * epsilon * rho[i]))
        y2[i] = rho[i] - (1 - rho[i]) * np.exp(beta*(mu[1] + 5 * epsilon * rho[i]))
        y3[i] = rho[i] - (1 - rho[i]) * np.exp(beta*(mu[2] + 5 * epsilon * rho[i]))
        
    plt.style.use('dark_background')    
    plt.plot(rho, y1, '-', color = 'g', label = "µ = " + str(mu[0]))
    plt.plot(rho, y2, '-', color = 'c', label = "µ = " + str(mu[1]))
    plt.plot(rho, y3, '-', color = 'b', label = "µ = " + str(mu[2]))
    plt.plot(rho, xAxis, '-', color = 'r')
    #plt.title("βε = "+ str(beta * epsilon) + ", µ = " + str(mu))
    plt.xlabel("ρ", size = 15)
    plt.ylabel("dΩ / dρ", size = 15)
    plt.legend(prop={'size': 15})

    
    # Finding the roots
    index1 = np.argmin(np.absolute(y1))
    index2 = np.argmin(np.absolute(y2))
    index3 = np.argmin(np.absolute(y3))
    root1 = rho[index1]
    root2 = rho[index2]
    root3 = rho[index3]
    
    # Adding the roots to the plot
    plt.plot(root1, 0, color = 'g', marker = 'd')
    plt.plot(root2, 0, color = 'c', marker = 'd')
    plt.plot(root3, 0, color = 'b', marker = 'd')
    
    # Add multiple roots to the beta = 1 case manually case because I am too lazy
    if beta == 1:
        plt.plot(0.5, 0, color = 'c', marker = 'd')
        plt.plot(0.85, 0, color = 'c', marker = 'd')
    
    plt.savefig("Minimising the Upper Bound for beta = " + str(beta) + ".pdf", bbox_inches='tight')
    plt.show()
        
    print("The roots are at " + str(root1) + ", " + str(root2) + ", " + str(root3))
    
    return 0

graph(1, 1, [-3, -2.5, -2])
graph(0.67, 1, [-3, -2.5, -2])


    