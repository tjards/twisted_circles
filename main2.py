#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 07:35:05 2023

@author: tjards

Note: Plot stuff here for exploration of various curves 
    https://www.desmos.com/calculator/65ps3jwzqn


"""

#%% import stuff
#-------------
from sympy import *
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import random

#%% useful functions

# Hamilton product
# ----------------
def hamilton_product(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return Matrix([w, x, y, z])

# quaternion conjugate
# --------------------
def quatjugate(q):
    w, x, y, z = q
    return Matrix([w, -x, -y, -z])

# quaternion rotation (p by q)
# ----------------------------
def quatrotate(q, p):
    
    if p.shape[0] != 4:
        print('error: express point as 4x1 quaternion')
        print('comment: typically, just appending a zero in the first row is sufficient')
    
    # rotation of point, p 
    # by quaterion, q
    rotated = hamilton_product(hamilton_product(q, p),quatjugate(q))
    return rotated


#%% parameters
# ------------
d, d2, r, R, psi, t = symbols('d d2 r R t psi')
a, b = symbols('a b', nonzero=True)             
t = symbols('t', nonzero=True)

#%% Reference: circle with radius d
# ----------------------------------
w = 0
x = d*cos(t)
y = d*sin(t)
z = 0
eqn_cir = Matrix([w,x,y,z])

#%% Inputs: Define what we're solving for
# ---------------------------------------
u = Matrix([a,b,0,0])           # given rotation about x-axis
RHS = quatrotate(u,eqn_cir)     # this is the Right Hand Side equation

#%% Targets: various curves with known parametric equations
SOLS = list()

# 2D // lemniscate of Gerono
# -------------------
w = 0
x = d*cos(t)
y = d*sin(t)*cos(t)
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8ger = Matrix([w,x,y,z])
EQNS = eqn_8ger - RHS
SOLS.append(nonlinsolve([EQNS[1], EQNS[2]], [a, b], S.Reals).simplify())

# 2D // lemniscate of Bernoulli
# -----------------------
# d is half-width: distance from crossing point to horiz extreme (tunable)
d1 = 2*d
w = 0
x = d1*cos(t)/(1+sin(t)*sin(t))
y = d1*sin(t)*cos(t)/(1+sin(t)*sin(t))
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8ber = Matrix([w,x,y,z])
EQNS = eqn_8ber - RHS
SOLS.append(nonlinsolve([EQNS[1], EQNS[2]], [a, b], S.Reals).simplify())

# 2D // dumbell curve (sextic curve), bowtie
# https://mathcurve.com/courbes2d.gb/doublegouttedeau/doublegouttedeau.shtml
# -----------------------------------
d2 = d
w = 0
x = d*cos(t)
y = ((d*d)/d2)*cos(t)*cos(t)*sin(t)
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8bow = Matrix([w,x,y,z])
EQNS = eqn_8bow - RHS
SOLS.append(nonlinsolve([EQNS[1], EQNS[2]], [a, b], S.Reals).simplify())


#%% Pull out and ensure unit quaternion (norm =1)
# ------------------------------------------
QUATS = list()

j = 0 # some of the above will have multiple solutions; choose first

for i in range(0,len(SOLS)):
    j = random.randint(0, 3)
    quaternion_terms = SOLS[i].args[j].simplify()
    isunit = sqrt(quaternion_terms[0]*quaternion_terms[0]+quaternion_terms[1]*quaternion_terms[1]).simplify()
    if isunit == 1:
        print("Solution ",i," is a unit quaternion")
    else:
        print("Solution ",i," has norm ", isunit)
        quaternion_terms = Matrix([quaternion_terms[0],quaternion_terms[1]])
        quaternion_terms = quaternion_terms/quaternion_terms.norm()
        print("... Solution ",i," has been normized ")
    QUATS.append(quaternion_terms)
    
        
#%% Print stuff
# -------------
radius = 4 # radius of circle to use
plots = list()
eqn_cir = eqn_cir.subs(d, radius)

for i in range(0,len(QUATS)):
    QED = quatrotate(Matrix([QUATS[i][0],QUATS[i][1],0,0]),eqn_cir)  
    QED = QED.subs(d, radius)
    plots.append(QED)       

#2D
plot_2D = plot_parametric(
                (eqn_cir[1], eqn_cir[2]),
                (plots[0][1], plots[0][2]),
                (plots[1][1], plots[1][2]),
                (plots[2][1], plots[2][2]),
                (t,-np.pi,np.pi), 
                legend = True,
                title = 'Projection onto x-y plane')           
plot_2D[0].line_color = 'black'
plot_2D[0].label = 'Circle'
plot_2D[1].line_color = 'blue'
plot_2D[1].label = 'Gerono'
plot_2D[2].line_color = 'red'
plot_2D[2].label = 'Bernoulli'
plot_2D[3].line_color = 'green'
plot_2D[3].label = 'Dumbell'
plot_2D.show()

#3D    
from sympy.plotting import plot3d_parametric_line   
plot_3D = plot3d_parametric_line(
                (eqn_cir[1], eqn_cir[2], eqn_cir[1]*0.01, (t,-np.pi,np.pi)),
                (plots[0][1], plots[0][2], plots[0][3], (t,-np.pi,np.pi)),
                (plots[1][1], plots[1][2], plots[1][3],(t,-np.pi,np.pi)),
                (plots[2][1], plots[2][2], plots[2][3],(t,-np.pi,np.pi)),
                (eqn_cir[1], eqn_cir[1]*0.01, eqn_cir[2], (t,-np.pi,np.pi)),
                legend = True)           
plot_3D[0].line_color = 'black'
plot_3D[4].line_color = 'black'
plot_3D[0].label = 'Circle'
plot_3D[4].label = ''
plot_3D[1].line_color = 'blue'
plot_3D[1].label = 'Gerono'
plot_3D[2].line_color = 'red'
plot_3D[2].label = 'Bernoulli'
plot_3D[3].line_color = 'green'
plot_3D[3].label = 'Dumbell'
plot_3D.show()







        
        
        