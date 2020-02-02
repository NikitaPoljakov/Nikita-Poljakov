import numpy as np
import matplotlib.pyplot as plt
from math import log10, floor
from scipy.optimize import curve_fit
from scipy.stats import linregress

# Rounding function
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

# Function to fit
def func(x, a, b):
    return a * x + b

#Define constants
masses = [0, 0.0135030, 0.035900, 0.0253356]
delta_masses = [0, 10**(-4), 10**(-3), 6*10**(-4)]
lengths = [0.27, 0.743, 0.51]
lengths_err = [0.01, 0.005, 0.01]
v = 0.01285
delta_v = 0.0003
r = 0.01
delta_r = 0.0005 #Error in diameter / 2
p = 101325
delta_p = 0.025 * 9.81 / (np.pi * r**2)
g = 9.81 

#Plotting and Fitting Data
xdata = masses[1:4]
xerr = delta_masses[1:4] 
ydata = lengths
yerr = lengths_err
y2data = []

# Plotting lengths vs masses
plt.plot(xdata, ydata, '.')
plt.errorbar(xdata, ydata, yerr, xerr,'r.')

# Fitting values to func (defined above) and errors, printing r-value
popt, pcov = curve_fit(func, xdata, ydata)
print("R-value = " + str(linregress(xdata, ydata)[2]))

# Creating linear fit data points
for i in range(3):
    y2data.append(func(xdata[i], popt[0], popt[1]))

#Gamma Value
gamma = 2*g*v/(popt[0]*p*(np.pi*r**2)**2) 
gamma_err = gamma * np.sqrt((delta_v/v)**2 + (pcov[0][0]/popt[0])**2 + (4*delta_r/r)**2 + (delta_p/p)**2)

# Plotting
plt.plot(xdata, y2data)
plt.title('Rinkel Plot for Air')
plt.ylabel('Drop Length (m)')
plt.xlabel('Mass (kg)')
plt.text(0.01, 0.07, "The gamma value is " + str(round(gamma, len(str(round_to_1(gamma_err))) - 2)) + u"\u00B1" + str(round_to_1(gamma_err)))
plt.savefig('Rinkel Air.pdf')
plt.show()


