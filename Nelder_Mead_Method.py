# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:01:57 2021

@author: ljubo
"""

#from matplotlib import ticker
#from scipy.optimize import curve_fit
#import numpy as np
#import scipy.integrate as integral
import math as sp
#import matplotlib.pyplot as plt
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


# Misc variables
length = 15;
F_temp = 0;
x_storage = [];
y_storage = [];
min_vertex = [];




# Third variable:

F_tot = [];



def vec_magn(a0,a1,b0,b1):
    magn = sp.sqrt((a1-a0)**2+(b1-b0)**2);
    return magn;



def Nelder_Mead(u0_x,u0_y,v0_x,v0_y,w0_x,w0_y):
    i = 0;
    while i < length:
        u_temp = force(u0_x,u0_y);
        v_temp = force(v0_x,v0_y);
        w_temp = force(w0_x,w0_y);
        if u_temp > v_temp and u_temp > w_temp:
            # Find midpoint on the opposite side
            m_x = 0.5*(v0_x+w0_x);
            m_y = 0.5*(v0_y+w0_y);
            
            # Find vector between worst point(u) and midpoint and take the opposite(-1)
            mu_x = (u0_x - m_x)*(-1);
            mu_y = (u0_y - m_y)*(-1);
        
            # Find new better point r
            r_x = m_x + mu_x;
            r_y = m_y + mu_y;
            
            # Evaluate at r
            r_temp =  force(r_x,r_y);
            print("Step U1 passed");
            if r_temp < u_temp:
                if r_temp < w_temp and r_temp < v_temp:
                    e_x = r_x + mu_x;
                    e_y = r_y + mu_y;
                    e_temp = force(e_x, e_y)
                    
                    u_temp = e_temp;
                    u0_x = e_x;
                    u0_y = e_y;
                    
                    min_vertex.append(u_temp);
                    i = i+1;
                    print("Step U1.1.1 passed");
                else:
                    u_temp = r_temp;
                    u0_x = r_x;
                    u0_y = r_y;
                    
                    min_vertex.append(u_temp);
                    i =i+1;
                    print("Step U1.1.2 passed");
                    
            elif r_temp > u_temp and r_temp > v_temp:
                ci_x = u0_x + 0.5*mu_x;
                ci_y = u0_y + 0.5*mu_y;
                
                c0_x = m_x + 0.5*mu_x;
                c0_y = m_y + 0.5*mu_y;
                
                ci_temp = force(ci_x, ci_y);
                c0_temp = force(c0_x, c0_y);
                print("Step U1.2 passed");
                if ci_temp < c0_temp:
                    u_temp = ci_temp;
                    u0_x = ci_x;
                    u0_y = ci_y;
                    
                    min_vertex.append(u_temp);
                    i = i+1;
                    print("Step U1.2.1 passed");
                else:
                    u_temp = c0_temp;
                    u0_x = c0_x;
                    u0_y = c0_y;
                    
                    min_vertex.append(u_temp);
                    i = i+1;
                    print("Step U1.2.2 passed");
            else: 
                w0_x = 0.5*(u0_x+w0_x);
                w0_y = 0.5*(u0_y+w0_y);
                w_temp = force(w0_x, w0_y);
                
                v0_x = 0.5*(u0_x+v0_x);
                v0_y = 0.5*(u0_y+v0_y);
                v_temp = force(v0_x, v0_y);
                
                min_vertex.append(u_temp);
                i = i+1;
                print("Step U1.3 passed");
                
        if v_temp > u_temp and v_temp > w_temp:
            # Find midpoint on the opposite side
            m_x = 0.5*(u0_x+w0_x);
            m_y = 0.5*(u0_y+w0_y);
            
            # Find vector between worst point(v) and midpoint and take the opposite(-1)
            mu_x = (v0_x - m_x)*(-1);
            mu_y = (v0_y - m_y)*(-1);
            
            # Find new better point r
            r_x = m_x + mu_x;
            r_y = m_y + mu_y;
            
            # Evaluate at r
            r_temp =  force(r_x,r_y);
            print("Step V1 passed");
            if r_temp < v_temp:
                if r_temp < w_temp and r_temp < u_temp:
                    e_x = r_x + mu_x;
                    e_y = r_y + mu_y;
                    e_temp = force(e_x, e_y)
                    
                    v_temp = e_temp;
                    v0_x = e_x;
                    v0_y = e_y;
                    
                    min_vertex.append(v_temp);
                    i = i+1;
                    print("Step V1.1.1 passed");
                else:
                    v_temp = r_temp;
                    v0_x = r_x;
                    v0_y = r_y;
                    
                    min_vertex.append(v_temp);
                    i = i+1;
                    print("Step V1.1.2 passed");
            elif r_temp > v_temp and r_temp > u_temp:
                ci_x = v0_x + 0.5*mu_x;
                ci_y = v0_y + 0.5*mu_y;
                
                c0_x = m_x + 0.5*mu_x;
                c0_y = m_y + 0.5*mu_y;
                
                ci_temp = force(ci_x, ci_y);
                c0_temp = force(c0_x, c0_y);
                print("Step V1.2 passed");
                if ci_temp < c0_temp:
                    v_temp = ci_temp;
                    v0_x = ci_x;
                    v0_y = ci_y;
                    
                    min_vertex.append(v_temp);
                    i = i+1;
                    print("Step V1.2.1 passed");
                    print(v0_x/R,v0_y/R)
                else:
                    v_temp = c0_temp;
                    v0_x = c0_x;
                    v0_y = c0_y;
                    
                    min_vertex.append(v_temp);
                    i = i+1;
                    print("Step V1.2.2 passed");
            else: 
                w0_x = 0.5*(v0_x+w0_x);
                w0_y = 0.5*(v0_y+w0_y);
                w_temp = force(w0_x, w0_y);
                
                u0_x = 0.5*(v0_x+u0_x);
                u0_y = 0.5*(v0_y+u0_y);
                u_temp = force(u0_x, u0_y);
                
                min_vertex.append(v_temp);
                i = i+1;
                print("Step V1.3 passed");
                
                
        if w_temp > v_temp and w_temp > u_temp:
            # Find midpoint on the opposite side
            m_x = 0.5*(v0_x+u0_x);
            m_y = 0.5*(v0_y+u0_y);
            
            # Find vector between worst point(u) and midpoint and take the opposite(-1)
            mu_x = (w0_x - m_x)*(-1);
            mu_y = (w0_y - m_y)*(-1);
            
            # Find new better point r
            r_x = m_x + mu_x;
            r_y = m_y + mu_y;
            
            # Evaluate at r
            r_temp =  force(r_x,r_y);
            print("Step W1 passed");
            if r_temp < w_temp:
                if r_temp < u_temp and r_temp < v_temp:
                    e_x = r_x + mu_x;
                    e_y = r_y + mu_y;
                    e_temp = force(e_x, e_y)
                    
                    w_temp = e_temp;
                    w0_x = e_x;
                    w0_y = e_y;
                    
                    min_vertex.append(w_temp);
                    i = i+1;
                    print("Step W1.1.1 passed");
                else:
                    w_temp = r_temp;
                    w0_x = r_x;
                    w0_y = r_y;
                    
                    min_vertex.append(w_temp);
                    i = i+1;
                    print("Step W1.1.2 passed");
                    
            elif r_temp > w_temp and r_temp > v_temp:
                ci_x = w0_x + 0.5*mu_x;
                ci_y = w0_y + 0.5*mu_y;
                
                c0_x = m_x + 0.5*mu_x;
                c0_y = m_y + 0.5*mu_y;
                
                ci_temp = force(ci_x, ci_y);
                c0_temp = force(c0_x, c0_y);
                print("Step W1.2 passed");
                if ci_temp < c0_temp:
                    w_temp = ci_temp;
                    w0_x = ci_x;
                    w0_y = ci_y;
                    
                    min_vertex.append(w_temp);
                    i = i+1;
                    print("Step W1.2.1 passed");
                else:
                    w_temp = c0_temp;
                    w0_x = c0_x;
                    w0_y = c0_y;
                    
                    min_vertex.append(w_temp);
                    i = i+1;
                    print("Step W1.2.2 passed");
            else: 
                u0_x = 0.5*(u0_x+w0_x);
                u0_y = 0.5*(u0_y+w0_y);
                u_temp = force(u0_x, u0_y);
                
                v0_x = 0.5*(w0_x+v0_x);
                v0_y = 0.5*(w0_y+v0_y);
                v_temp = force(v0_x, v0_y);
                
                min_vertex.append(w_temp);
                i = i+1;
                print("Step W1.3 passed");
    print(min_vertex);            
    return min_vertex[length-1];           

#u0_x = 6*R;
#u0_y = 0*R; 

#v0_x = 4*R;
#v0_y = 0*R;

#w0_x = 5*R;
#w0_y = 1*R;

print(Nelder_Mead(5.5*R, 4.5*R, 4.5*R, 2.5*R, 5*R, 3*R))
print("The force is:",force(3.75*R,2.25*R));


