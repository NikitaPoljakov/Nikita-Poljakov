import numpy as np
import matplotlib.pyplot as plt

N = 1000 # Number of points
f = 1 # focal length
x_circle = np.zeros(N) # x-coordinates of the circle
y_upper_circle = np.zeros(N) # y-coordinates of the upper semicircle
y_lower_circle = np.zeros(N) # y-coordinates of the lower semicircle

print("Welcome to the optics simulator where the real image of a circle through a convex lens of focal distance 1 is computed. ")
x_0 = float(input("Select the x-coordinate of the circle in the range (-5, 0): "))
y_0 = float(input("Select the y-coordinate of the circle in the range (-5, 5): "))
R = float(input("Select the radius of the circle in the range (0, 5): "))

# Get x-coordinates of the circle
x = - R + x_0 # Initial x-coordinate
for i in range(N):
    x_circle[i] = x
    x += 2.01 * R / N
    
    y_upper_circle[i] = y_0 + np.sqrt(float(R**2 - (x_circle[i] - x_0)**2)) # Get y-coordinates of the upper semicircle
    y_lower_circle[i] = y_0 - np.sqrt(float(R**2 - (x_circle[i] - x_0)**2)) # Get y-coordinates of the lower semicircle

# Lens coordinates
x_lens = [0, 0]
y_lens = [-5, 5]

# Optical axis coordinates
x_axis = [-5, 5]
y_axis = [0, 0]

# Get Image coordinates
x_new = np.zeros(N)
y_upper_new = np.zeros(N)
y_lower_new = np.zeros(N)

for i in range(N):
    x_new[i] = 1 / ((1 / f) - (1 / x_circle[i]))
    y_upper_new[i] = y_upper_circle[i] / ((x_circle[i] / f) - 1) 
    y_lower_new[i] = y_lower_circle[i] / ((x_circle[i] / f) - 1) 

plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.plot(x_circle, y_upper_circle, 'm')
plt.plot(x_circle, y_lower_circle, 'm')
plt.plot(x_lens, y_lens, 'b-', linewidth = 4)
plt.plot(x_axis, y_axis, 'b:')
plt.plot(x_new, y_upper_new, 'r')
plt.plot(x_new, y_lower_new, 'r')
# Plot foci
plt.plot(-1, 0, 'bo')
plt.plot(1, 0, 'bo')
plt.show()


