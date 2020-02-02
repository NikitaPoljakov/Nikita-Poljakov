import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import log10, floor
from scipy.stats import linregress

# Rounding function
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

# Determining peaks and their times
def Peaks(filename, start = 0, stop = -1):
    
    # File selection
    file = open(filename + '.txt', 'r')

    # Creating arrays
    time, pressure, rel_time, delta_time = [], [], [], []

    # Reading data and filling in pressure
    n = 0
    for line in file:
        if n >= 8:
            column = line.split()
            time.append(column[2])
            pressure.append(float(column[4]))
            n += 1
            
        elif n < 8:
            n+= 1

    file.close()
            
    # Generating time array
    for m in range(len(time)):
        col2 = time[m].split(':')
        rel_time.append(int(col2[0])*3600 + int(col2[1])*60 + float(col2[2]))
        delta_time.append(rel_time[m] - rel_time[0]) 

    # Finding peaks
    x = pressure[start:stop]
    peaks, _ = find_peaks(x, height=1, distance = 50)

    xdata = [delta_time[p] for p in peaks] 
    ydata = [x[p] for p in peaks]
    
    return [xdata, ydata]                                                             

# Function to fit for envelope
def func1(x, a, b, c):
    return a + b * np.exp(-c * x)

# Function to fit for linear fit to obtain b
def func2(x, d, e):
    return d * x + e

# Fitting and plotting, printing the exponent factor (b/2m)
def FitnPlot(filename):
    info = Peaks(filename)
    xdata = info[0]
    ydata = info[1]
    y2data = []
    popt, pcov = curve_fit(func1, xdata, ydata)
    print(str(popt[2]) + '\u00B1' + str(pcov[2][2]))
    plt.plot(xdata, ydata, 'b.')

    for i in range(len(xdata)):
        y2data.append(func1(xdata[i], *popt))
    
    plt.plot(xdata, y2data, 'r-')
    plt.title('')
    plt.ylabel('Relative Pressure')
    plt.xlabel('Time (s)')
    plt.savefig(filename + ' envelope.pdf')
    plt.show()
    
    return [popt[2], pcov[2][2]]

# Define constants
v = 0.01285
delta_v = 0.0001
r = 0.01
area = np.pi * r**2
delta_r = 0.0005 #Error in diameter / 2
p = 101325
delta_p = 0.025 * 9.81 / (np.pi * r**2)
delta_masses = [0, 10**(-4), 6*10**(-4), 7*10**(-4)]
masses = [0, 0.0135028, 0.0358996, 0.0253353]

for i in range(1, 4):
    for j in range(1, 4):
        x = FitnPlot('Air Ball ' + str(i) + ' Run ' + str(j))

#Plotting and Fitting Data
xdata = []
xerr = []
for i in range(1, 4):
    xdata.append(1/masses[i])
    xerr.append(delta_masses[i]/masses[i]**2)
ydata = [0.33, 0.18, 0.23]
yerr = [0.01, 0.01, 0.01]
y2data = []

# Plotting lengths vs masses
plt.plot(xdata, ydata, '.')
plt.errorbar(xdata, ydata, yerr, xerr,'r.')

# Fitting values to func (defined above) and errors
popt, pcov = curve_fit(func2, xdata, ydata)
print("R-value = " + str(linregress(xdata, ydata)))

# Creating linear fit data points
for i in range(3):
    y2data.append(func2(xdata[i], popt[0], popt[1]))

# Plotting
plt.plot(xdata, y2data)
plt.title('')
plt.ylabel('Exponent of the Envelope (s^(-1))')
plt.xlabel('Mass^-1 (kg^(-1))')
plt.savefig('Damping Coefficient.pdf')
plt.show()

print(2 * popt[0], round_to_1(2*pcov[0][0]))











