# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 08:16:35 2021

@author: ljubo
"""

#from matplotlib import ticker
#from scipy.optimize import curve_fit
import numpy as np
#import scipy.integrate as integral
import math as sp
import random as rd
import matplotlib.pyplot as plt
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D


# Constants 

G = 6.67408*10**(-11); #[m^3*kg^-1*s^-2]
R = 6371*10**3; #[m] Earth's radius and distance unit for initial coordinates


m_A = 5.972*10**24 #[kg] Mass of the asteroid belt
m_B = 5.972*10**24 #[kg] Mass of the moon
m_C = 5.972*10**24 #[kg] Mass of the earth

m_proj = 6.3 #[kg] Mass of the projectile


# Coordinates

x_A = 8*R;
y_A = 3*R;

x_B = 8*R;
y_B = 8*R;

x_C = 2*R;
y_C = 5*R;

x_0 = 5*R;
y_0 = 0.5*R;

x_fin = 5*R;
y_fin = 9.5*R;


# Define function f(x,y)

def force(x,y):
    d_A_magn = sp.sqrt((x_A-x)**2+(y_A-y)**2);
    d_B_magn = sp.sqrt((x_B-x)**2+(y_B-y)**2);
    d_C_magn = sp.sqrt((x_C-x)**2+(y_C-y)**2);
    
    d_Ax = x_A - x;
    d_Ay = y_A - y;

    d_Bx = x_B - x;
    d_By = y_B - y;

    d_Cx = x_C - x;
    d_Cy = y_C - y;
    
    F_Ax = G*m_A*m_proj*d_Ax/(d_A_magn**3);
    F_Ay = G*m_A*m_proj*d_Ay/(d_A_magn**3);
    
    F_Bx = G*m_B*m_proj*d_Bx/(d_B_magn**3);
    F_By = G*m_B*m_proj*d_By/(d_B_magn**3);
    
    F_Cx = G*m_C*m_proj*d_Cx/(d_C_magn**3);
    F_Cy = G*m_C*m_proj*d_Cy/(d_C_magn**3);
    
    F_x = F_Ax + F_Bx + F_Cx;
    F_y = F_Ay + F_By + F_Cy;
    
    F = sp.sqrt(F_x**2 + F_y**2);
    
    return F;


# Define random variable storages
x_storage = [];
y_storage = [];
f_storage = [];


# Misc variables
r_x = 0; # Random x coordinate
r_y = 0; # Random y coordinate
f_temp = 0;
i = 0;



# Find random points and store
while i<10000:
    r_x = rd.uniform(0*R, 10*R);
    r_y = rd.uniform(0*R, 10*R);
    
    x_storage.append(r_x);
    y_storage.append(r_y);
    
    f_temp = force(r_x, r_y);
    f_storage.append(f_temp);
    
    i = i+1;
  
    
# Find coordinates of the minimum and plot a circle
min_f = min(f_storage);
min_index_f = f_storage.index(min_f);

min_fx = x_storage[min_index_f];
min_fy = y_storage[min_index_f];



# Draw A:
x_A_circ = [];
y_A_circ = [];
angles = np.linspace(0,2*sp.pi, 1000);
for j in range(1000):
    x_A_circ.append(x_A+sp.cos(angles[j])*1*R);  # Radius of A is 0.3*R
    y_A_circ.append(y_A+sp.sin(angles[j])*1*R);
    

# Draw B:
x_B_circ = [];
y_B_circ = [];
angles = np.linspace(0,2*sp.pi, 1000);
for j in range(1000):
    x_B_circ.append(x_B+sp.cos(angles[j])*1*R);
    y_B_circ.append(y_B+sp.sin(angles[j])*1*R);      
    

# Draw C:
x_C_circ = [];
y_C_circ = [];
angles = np.linspace(0,2*sp.pi, 1000);
for j in range(1000):
    x_C_circ.append(x_C+sp.cos(angles[j])*1*R);
    y_C_circ.append(y_C+sp.sin(angles[j])*1*R);    
    
    
    
# Draw the initial position:
x_0_circ = [];
y_0_circ = [];
angles = np.linspace(0,2*sp.pi, 1000);
for j in range(1000):
    x_0_circ.append(x_0+sp.cos(angles[j])*0.4*R);
    y_0_circ.append(y_0+sp.sin(angles[j])*0.4*R);     
    
    
# Draw the final position:
x_fin_circ = [];
y_fin_circ = [];
angles = np.linspace(0,2*sp.pi, 1000);
for j in range(1000):
    x_fin_circ.append(x_fin+sp.cos(angles[j])*0.4*R);
    y_fin_circ.append(y_fin+sp.sin(angles[j])*0.4*R);   

# Draw the minimum's position:
x_min_circ = [];
y_min_circ = [];
angles = np.linspace(0,2*sp.pi, 1000);
for j in range(1000):
    x_min_circ.append(min_fx+sp.cos(angles[j])*0.4*R);
    y_min_circ.append(min_fy+sp.sin(angles[j])*0.4*R);  
    
    
# Plotting

"""
plt.xlim([0, 10])
plt.ylim([0, 10])
plt.plot(x_min_circ, y_min_circ);
plt.plot(x_fin_circ, y_fin_circ);
plt.plot(x_0_circ, y_0_circ);
plt.plot(x_A_circ,y_A_circ);
plt.plot(x_B_circ,y_B_circ);
plt.plot(x_C_circ,y_C_circ);

plt.axis('equal')
plt.show()
"""

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(x_0_circ, y_0_circ)
ax.plot(x_fin_circ, y_fin_circ)
ax.plot(x_min_circ, y_min_circ, 'k--', label="{c_x:.2f},{c_y:.2f}\n F = {c_f:.2f} N".format(c_x = min_fx/R, c_y = min_fy/R, c_f = min_f))
ax.plot(x_A_circ,y_A_circ)
ax.plot(x_B_circ,y_B_circ)
ax.plot(x_C_circ,y_C_circ)
ax.set_xlim(xmin=0.0, xmax=10*R)
ax.set_ylim(ymin=0.0, ymax=10*R)
ax.set_aspect('equal');
ax.set_xlabel("x-coordinate")
ax.set_ylabel("y-coordinate")
legend = ax.legend(loc='upper left', shadow=True, fontsize='medium')
