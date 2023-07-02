#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 07:35:05 2023

@author: tjards

play here:
    
    https://www.desmos.com/calculator/65ps3jwzqn


"""

# import stuff
#-------------
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
#import time


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
    
    

#%% parametric equations for various curves
# -----------------------------------------

# variables
var('a b d r R t psi t z', positive = True)

# hypotrochoidal curves
# ----------------------
psi = (R-r)/r
w = 0
x = (R-r)*cos(t)+d*cos(psi*t)
y = (R-r)*sin(t)-d*sin(psi*t)
z = 0
eqn_hyp = Matrix([w,x,y,z])
#pprint(eqn_hyp)

# lemniscate of Gerono
# -------------------
w = 0
x = d*cos(t)
y = d*sin(t)*cos(t)
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8_1 = Matrix([w,x,y,z])
#pprint(eqn_8_1)

# lemniscate of Bernoulli
# -----------------------
c2 = d   # half-width: distance from crossing point to horiz extreme (tunable)
w = 0
x = c2*cos(t)/(1+sin(t)*sin(t))
y = c2*sin(t)*cos(t)/(1+sin(t)*sin(t))
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8_2 = Matrix([w,x,y,z])
#pprint(eqn_8_2)

# Vivians curve
# https://mathworld.wolfram.com/VivianisCurve.html
# -------------
c3 = d
w = 0
x = c3*(1+cos(t))
y = c3*sin(t)
z = 2*c3*sin(t/2)
eqn_8_3 = Matrix([w,x,y,z])

# dumbell curve (sextic curve), bowtie
# https://mathcurve.com/courbes2d.gb/doublegouttedeau/doublegouttedeau.shtml
# -----------------------------------
c4 = d
c5 = d
w = 0
x = c4*cos(t)
y = ((c4*c4)/c5)*cos(t)*cos(t)*sin(t)
z = 0
eqn_8_4 = Matrix([w,x,y,z])

#Cayley's sextet
# --------------
w = 0
x = cos(t)*cos(t)*cos(t)*cos(3*t)
y = cos(t)*cos(t)*cos(t)*sin(3*t)
z = 0
eqn_8_5 = Matrix([w,x,y,z])


# # lemniscate of Booth (for later exploration)

# #circle
# ca = d
# cb = d
# #oval
# ca = d   
# cb = d/sqrt(2)
# # centered conic
# ca = sqrt(d)   # there are conditions for these (oval, indented oval, lemniscate, isolated circles, figure-8)
# cb = sqrt(d)
# w = 0
# x = cb*ca*cb*cos(t)/(cb*cb*cos(t)*cos(t) + ca*ca*sin(t)*sin(t))
# y = ca*ca*cb*sin(t)/(cb*cb*cos(t)*cos(t) + ca*ca*sin(t)*sin(t))
# z = 0 
# eqn_8_3 = Matrix([w,x,y,z])
# pprint(eqn_8_3)




#%% circle with radius d
# ----------------------
w = 0
x = d*cos(t)
y = d*sin(t)
z = 0
eqn_cir = Matrix([w,x,y,z])
#pprint(eqn_cir)

#%% unknown parameters (for which we will solve)
u = Matrix([a,b,0,0])  # given rotation about x and y
#u2 = quatjugate(u1)
#pprint(u)
#pprint(u2)




'''
Unnecessary complicated thing:
Let us solve for a, b that generate a hypotrochoidal curves as
    a function of projection of quaternion rotations around x,y

'''
#%%
RHS = quatrotate(u,eqn_cir)
LHS = eqn_8_5
#pprint(LHS)
#print('=')
#pprint(RHS)

EQNS = LHS-RHS
#pprint(EQNS)

#%%
solutions = nonlinsolve([EQNS[1], EQNS[2]], [a, b])

#equations = Eq(EQNS[1],EQNS[2],EQNS[3])
#solutions = solve(equations, a)
#solutions.args[0] # when not a list
#pprint(solutions)
pprint(solutions.args[0].simplify())

#%% plot stuff

# plt.figure()
# plt.axis([-5,5,-5,5])
# r1   = 1
# #plt_type = 'circle'
# #plt_type = 'gerono'
# #plt_type = 'bernoulli'
# plt_type = 'booth'

# if plt_type == 'booth':
    
#     # for oval: r/sqrt(2) < r2 < r*sqrt(2)
#     #r2  = r*sqrt(2)
#     r2  = r1/np.sqrt(2)

#     # for circle: r2 = r1
#     r2  = r1


# xs = 0
# ys = 0

# #for rs in np.arange(-2*r1/np.sqrt(2),r1*np.sqrt(2)+2*r1/np.sqrt(2),0.25*r1/np.sqrt(2)):
# #for rs in np.arange(-2*r1^2,2*r1^3,0.25*r1/np.sqrt(2)):
# #for rs in np.arange(r1*sqrt(2),r1/sqrt(2)-r1/sqrt(2),-0.2):
# for rs in [0,0.1,0.2,0.5,1,1.5,2]:
    
#     r2 = 1
#     r1 = rs
    
#     #if r2 >= - 0.2 and r2 <= 0.2:
        
#         #continue

#     for t in np.arange(0,2*np.pi,0.1):
        
#         xs_prev = xs
#         ys_prev = ys
        
#         if plt_type == 'circle':
#             xs = r1*cos(t)
#             ys = r1*sin(t)
        
#         if plt_type == 'gerono':
#             xs = r1*cos(t)
#             ys = r1*sin(t)*cos(t)
        
#         if plt_type == 'bernoulli':
#             xs = r1*cos(t)/(1+sin(t)*sin(t))
#             ys = r1*sin(t)*cos(t)/(1+sin(t)*sin(t))
            
#         if plt_type == 'booth':
#             xs = r2*r1*r2*cos(t)/(r2*r2*cos(t)*cos(t) + r1*r1*sin(t)*sin(t))
#             ys = r2*r1*r2*sin(t)/(r2*r2*cos(t)*cos(t) + r1*r1*sin(t)*sin(t))
        
#         plt.axis([-5,5,-5,5])
#         plt.text(-1.5,-1.2,'r1: '+str(r1))
#         plt.text(-1.5,-1.5,'r2: '+str(r2))
#         plt.pause(0.01)
#         plt.plot([xs_prev,xs],[ys_prev, ys],'-bo')

        
#     plt.clf()
# plt.show()











# test examples
# --------
# var('w x y z')
# q1 = Matrix([w, x, y, z])
# var('a b c d')
# q2 = Matrix([a, b, c, d])
# product = hamilton_product(q1, q2)
# pprint(product)


    



