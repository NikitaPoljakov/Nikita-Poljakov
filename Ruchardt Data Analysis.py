import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from math import log10, floor
from scipy.optimize import curve_fit
from scipy.stats import linregress

def Data(filename, masses, delta_masses, start = 0, stop = -1):
    
    # File selection
    file = open(filename + '.txt', 'r')
    
    # Creating arrays
    time, pressure, rel_time, delta_time = [], [], [], []
    periods = []
    
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

    
    # I don't understand
    x = pressure[start:stop]
    
    peaks, _ = find_peaks(x, height=0, distance = 60)
    #print(peaks)
    peak_times = [delta_time[p] for p in peaks] 
    No_peaks = len(peak_times)

    # Plotting data
    plt.plot(delta_time[start:stop], pressure[start:stop])
    plt.title('Graph for ' + filename)
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
    
    mean_mean.append(mean)
    


# Rounding function
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

# Function to fit
def func(x, a, b):
    return a * x + b

#_____________________________________ Analysis _______________________________


# Define constants
delta_masses = [0, 10**(-1), 6*10**(-1), 7*10**(-1)]
masses = [0, 0.0135028, 0.0358996, 0.0253353]
lengths = [0.223, 0.622, 0.433]
lengths_err = [0.005, 0.005, 0.005]
v = 0.01285
delta_v = 0.0001
r = 0.01
delta_r = 0.0005 #Error in diameter / 2
p = 101325
delta_p = 0.025 * 9.81 / (np.pi * r**2)
g = 9.81

mean_mean, err_mean = [], []

for i in range(1, 4):
    for j in range(1, 4):
            Data('Air Ball' + ' ' + str(i) + ' ' + 'Run' + ' ' + str(j), masses[i], delta_masses[i]) 

# Creating Variables for mean periods & errors for each ball
ball1, ball2, ball3 = 0, 0, 0
err1, err2, err3 = 0, 0, 0

# Summing mean periods
for i in range(len(mean_mean)):
    if i <= 2:
        ball1 += mean_mean[i]
    elif i <= 5:
        ball2 += mean_mean[i]
    else:
        ball3 += mean_mean[i]

# Finding mean for each ball
ball1 = ball1/3
ball2 = ball2/3
ball3 = ball3/3

# Error analysis
err1 = 1/3 * ball1 * np.sqrt(err_mean[0]**2 + err_mean[1]**2 + err_mean[2]**2)
err2 = 1/3 * ball2 * np.sqrt(err_mean[3]**2 + err_mean[4]**2 + err_mean[5]**2)
err3 = 1/3 * ball3 * np.sqrt(err_mean[6]**2 + err_mean[7]**2 + err_mean[8]**2)

print(ball1, u"\u00B1", err1)
print(ball2, u"\u00B1", err2)
print(ball3, u"\u00B1", err3)      
        
xdata = [ball1**2, ball2**2, ball3**2]
ydata = masses[1:4]
for i in range(len(ydata)):
    ydata[i] *= 1000 
y2data = []
yerr = delta_masses[1:4]
xerr = [2 * err1, 2 * err2, 2 * err3]

plt.plot(xdata, ydata, '.')
plt.errorbar(xdata, ydata, yerr, xerr, 'r.')

popt, pcov = curve_fit(func, xdata, ydata)
print("R-value = " + str(linregress(xdata, ydata)))

for i in range(3):
    y2data.append(func(xdata[i], popt[0], popt[1]))

#Gamma value
gamma = popt[0]*4*v/(1000*p*r**4)
gamma_err = gamma * np.sqrt((0.05)**2 + (0.01)**2 + (4*0.05)**2 + (0.01)**2)


plt.plot(xdata, y2data)
plt.title('Final Graph for Ruchardt with Air')
plt.xlabel('Period^2 (s^2)')
plt.ylabel('Masses (g)')
plt.savefig('Ruchardt Air Day 1.pdf')
plt.show()

print(gamma, gamma_err)
  