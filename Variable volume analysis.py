import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from math import log10, floor
from scipy.optimize import curve_fit
from scipy.stats import linregress

#_________________________ Defining Function __________________________________

def Data(filename, volumes, start = 0, stop = -1):
    
    # File selection
    file = open(filename + '.txt', 'r')
    
    # Creating arrays
    time, pressure, rel_time, delta_time, periods = [], [], [], [], []
    
    # Reading data and filling in pressure
    n = 0 # n is used to skip the inital information lines
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

    
    # Fnding and appending pressure peaks to peaks
    x = pressure[start:stop]
    peaks, _ = find_peaks(x, height=-10, distance = 1)
    peak_times = [delta_time[p] for p in peaks] 
    No_peaks = len(peak_times)

    # Plotting data
    plt.plot(delta_time[start:stop], pressure[start:stop])
    plt.title('')
    plt.xlabel('Time (s)')
    plt.ylabel('Relative Pressure')
    plt.savefig(filename + '.pdf')
    
    # Calculating periods
    for i in range(No_peaks-1):
        periods.append(peak_times[i+1] - peak_times[i])

    # Statistics
    mean = (peak_times[-1] - peak_times[0])/len(periods)
    sq_sum = 0
    for i in range(len(periods)):
        sq_sum += (periods[i] - mean)**2
        stan_div = np.sqrt(sq_sum/len(periods)) 
        err_mean.append(stan_div)

    # Outputs
    plt.show()
    print('Number of peaks is ', No_peaks)
    print('Number of periods is: ', len(periods))
    print(periods)
    print('Mean period is: ', mean, u"\u00B1", stan_div)
                                                                  
    means.append(mean)

# Rounding function
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

# Function to fit
def func(x, a, b):
    return a * x + b

#_____________________________________ Analysis _______________________________

# Define constants
r = 0.0171 #Change this after we find the value of r
delta_r = 0.001 #Error in diameter / 2
p = 101325
delta_p = 0.1 * 9.81 / (np.pi * r**2) # Improve upon this approximation
mass = 0.0952506
delta_mass = 0.0004

# Running for all files
means = []
err_mean = []
volumes = [] #volumes in m**3
volumesname = [45, 50, 54, 60, 64, 69, 74, 79, 85, 88, 90, 91, 93, 94, 100] #volumes in ml used to call the files

for i in range(len(volumesname)):
    volumes.append(volumesname[i])

for i in range(len(volumes)):
    Data(str(volumesname[i]), volumes[i])   
       
# Assigning variables for analysis of data fitting
xdata = []

for i in range(len(means)):
    xdata.append((means[i]*1000)**2)
 
ydata = volumes
y2data = [] # Linear fit y data according to model
yerr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
xerr = []

for i in range(len(volumes)):
    xerr.append(2 * err_mean[i] * means[i] * 10**6) # Here is the place where a mistake was made, should multiplythe error by period. I multiplied by 10^6 because I converted to ms

# Plotting masses versus periods^2
plt.plot(xdata, ydata, '.')
plt.errorbar(xdata, ydata, yerr, xerr,'r.')

# Fitting values to func (defined above) and errors
popt, pcov = curve_fit(func, xdata, ydata)
print("R-value = " + str(linregress(xdata, ydata)))

#y-values for linear fit
for i in range(15):
    y2data.append(func(xdata[i], popt[0], popt[1]))

plt.plot(xdata, y2data)
plt.title('')
plt.xlabel('Period^2 (ms^2)')
plt.ylabel('Volume (ml)')
plt.savefig('VarVolFinal.pdf')
plt.show()






    
    