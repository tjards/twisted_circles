#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 07:35:05 2023

@author: tjards

This program takes the parametric equations for arbitrary simple closed 
curves and finds the corresponding functional rotation that "deforms" 
a circle into a 3D curve for which the 2D projection is this simple
closed curve.

The rotation is accomplished using a functional unit quaternion 
which can be shown to be a 1-parameter homeomorphism. 

The problem is indeterminate, so there are multiple solutions; we randomly
select for the purpose of plotting the resultant curves.

Points are expressed as "pure" quaternions (i.e. p=(0,x,y,z))
Rotation are expressed as "unit" quaternions (i.e. norm(q) = 1)
 
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
        print('comment: for true quaternions, just appending a zero in the first row is sufficient')
    
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

#%% Inputs: Define what we're solving for (i.e. the rotation)
# -----------------------------------------------------------
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


# 2D // Lemniscatic arch (note the phase shift)
# ---------------------------------------------
w = 0
x = d*cos(t)
y = d*sin(t)*cos(t+pi/3)
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8boo = Matrix([w,x,y,z])
EQNS = eqn_8boo - RHS
SOLS.append(nonlinsolve([EQNS[1], EQNS[2]], [a, b], S.Reals).simplify())





#%% Pull out and ensure unit quaternion (norm = 1)
# ------------------------------------------
QUATS = list()

j         = 0 # some of the above will have multiple solutions; choose first
normalize = 0 # enforce unit quaternions?


for i in range(0,len(SOLS)):
    #j = random.randint(0, 3)
    quaternion_terms = SOLS[i].args[j].simplify()
    isunit = sqrt(quaternion_terms[0]*quaternion_terms[0]+quaternion_terms[1]*quaternion_terms[1]).simplify()
    if isunit == 1:
        print("Solution ",i," is a unit quaternion")
    else:
        print("Solution ",i," has norm ", isunit)
        quaternion_terms = Matrix([quaternion_terms[0],quaternion_terms[1]])
        if normalize == 1:
            quaternion_terms = quaternion_terms/quaternion_terms.norm()
            print("... Solution ",i," has been normized")
        else:
            print("... Solution ",i," has been not been normized")
    QUATS.append(quaternion_terms)
    
#%% Print stuff
# -------------
radius = 3 # radius of circle to use
plots = list()
eqn_cir = eqn_cir.subs(d, radius)

for i in range(0,len(QUATS)):
    QED = quatrotate(Matrix([QUATS[i][0],QUATS[i][1],0,0]),eqn_cir)  
    QED = QED.subs(d, radius)
    plots.append(QED)       

# 2D projection (x-y)
plot_2D = plot_parametric(
                (eqn_cir[1], eqn_cir[2]),
                (plots[0][1], plots[0][2]),
                (plots[1][1], plots[1][2]),
                (plots[2][1], plots[2][2]),
                (plots[3][1], plots[3][2]),
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
plot_2D[4].line_color = 'magenta'
plot_2D[4].label = 'Lemniscatic arch'
plot_2D.show()


# 2D projection
plot_2D = plot_parametric(
                (eqn_cir[1], eqn_cir[2]),
                (plots[0][2], plots[0][3]),
                (plots[1][2], plots[1][3]),
                (plots[2][2], plots[2][3]),
                (plots[3][2], plots[3][3]),
                (t,-np.pi,np.pi), 
                legend = True,
                title = 'Projection onto y-z plane')           
plot_2D[0].line_color = 'black'
plot_2D[0].label = 'Circle'
plot_2D[1].line_color = 'blue'
plot_2D[1].label = 'Gerono'
plot_2D[2].line_color = 'red'
plot_2D[2].label = 'Bernoulli'
plot_2D[3].line_color = 'green'
plot_2D[3].label = 'Dumbell'
plot_2D[4].line_color = 'magenta'
plot_2D[4].label = 'Lemniscatic arch'
plot_2D.show()


#3D plot   
from sympy.plotting import plot3d_parametric_line   
plot_3D = plot3d_parametric_line(
                (eqn_cir[1], eqn_cir[2], eqn_cir[1]*0.01, (t,-np.pi,np.pi)),
                (plots[0][1], plots[0][2], plots[0][3], (t,-np.pi,np.pi)),
                (plots[1][1], plots[1][2], plots[1][3],(t,-np.pi,np.pi)),
                (plots[2][1], plots[2][2], plots[2][3],(t,-np.pi,np.pi)),
                (plots[3][1], plots[3][2], plots[3][3],(t,-np.pi,np.pi)),
                (eqn_cir[1], eqn_cir[1]*0.01, eqn_cir[2], (t,-np.pi,np.pi)),
                title = 'Curves formed from deformed circles',
                legend = True)           
plot_3D[0].line_color = 'black'
plot_3D[5].line_color = 'black'
plot_3D[0].label = 'Circle'
plot_3D[5].label = ''
plot_3D[1].line_color = 'blue'
plot_3D[1].label = 'Gerono'
plot_3D[2].line_color = 'red'
plot_3D[2].label = 'Bernoulli'
plot_3D[3].line_color = 'green'
plot_3D[3].label = 'Dumbell'
plot_3D[4].line_color = 'magenta'
plot_3D[4].label = 'Lemniscatic arch'
plot_3D.show()

#%% Animate one curve
# -------------------
pick = 3 # which curve to pick from the list?

#plot the 3-d trajectory
ax = plt.figure().add_subplot(projection='3d')
ax.set_xlim3d([-5, 5])
ax.set_ylim3d([-5, 5])
ax.set_zlim3d([-5, 5])

# initial state (cartesian)
x = np.array([0,radius,0,0])
for i in np.arange(0,2*np.pi,0.1):
    quaternion_terms = QUATS[pick]
    x_prev = x
    eqn_cir_ = eqn_cir.subs({d:radius, t:i})
    quaternion_terms_ = quaternion_terms.subs({d:radius, t:i})
    pos = quatrotate(np.array([quaternion_terms_[0],quaternion_terms_[1],0,0]),eqn_cir_)
    x = pos
    plt.pause(0.001)
    ax.plot([x_prev[1],x[1]],[x_prev[2], x[2]],[x_prev[3], x[3]],'-b.')
    
#ax.title('Curve')
plt.show()






        
        
        