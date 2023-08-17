#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 07:35:05 2023

@author: tjards

Plot stuff here for exploration:
    
    https://www.desmos.com/calculator/65ps3jwzqn


"""

#%% import stuff
#-------------
from sympy import *
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
#import time


#%% useful functions
# ------------------

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
#var('a b d d2 r R t psi t z', positive = True)
# let is use the "d" parameter as the radius of reference circle
#a, b, d, d2, r, R, t, psi, t = symbols('a b d d2 r R t psi t')
d, d2, r, R, psi, t = symbols('d d2 r R t psi')
a, b = symbols('a b', nonzero=True)
t = symbols('t', nonzero=True)


#%% circle with radius d
# ----------------------
w = 0
x = d*cos(t)
y = d*sin(t)
z = 0
eqn_cir = Matrix([w,x,y,z])
#pprint(eqn_cir)

# 2D // lemniscate of Gerono
# -------------------
w = 0
x = d*cos(t)
y = d*sin(t)*cos(t)
z = 0.5*d*sin(t)*sin(t)
eqn_8ger = Matrix([w,x,y,z])
#pprint(eqn_8_1)

# 2D // lemniscate of Bernoulli
# -----------------------
# d is half-width: distance from crossing point to horiz extreme (tunable)
d1 = d
w = 0
x = d1*cos(t)/(1+sin(t)*sin(t))
y = d1*sin(t)*cos(t)/(1+sin(t)*sin(t))
z = 0 #0.5*d*sin(t)*sin(t)
eqn_8ber = Matrix([w,x,y,z])
#pprint(eqn_8_2)

# 2D // dumbell curve (sextic curve), bowtie
# https://mathcurve.com/courbes2d.gb/doublegouttedeau/doublegouttedeau.shtml
# -----------------------------------
d2 = d
w = 0
x = d*cos(t)
y = ((d*d)/d2)*cos(t)*cos(t)*sin(t)
z = 0#0.5*d*sin(t)*sin(t)
eqn_8bow = Matrix([w,x,y,z])

# 2D // Cayley's sextet
# --------------
w = 0
x = d*cos(t)*cos(t)*cos(t)*cos(3*t)
y = d*cos(t)*cos(t)*cos(t)*sin(3*t)
z = 0
eqn_cay = Matrix([w,x,y,z])

#%% Define what we're solving for
u = Matrix([a,b,0,0])  # given rotation about x and y
#pprint(u)


#%% Define the equations
RHS = quatrotate(u,eqn_cir)

# travis - this is where you enter the curve:
LHS = eqn_8bow

#pprint(LHS)
#print('=')
#pprint(RHS)

EQNS = LHS-RHS
#pprint(EQNS)

#%% Solve for u
solutions = nonlinsolve([EQNS[1], EQNS[2]], [a, b], S.Reals)
#solutions = nonlinsolve([EQNS[1], EQNS[2], EQNS[3]], [a, b])
#pprint(solutions.args[0].simplify())


#%% Try it
radius = 5
quaternion_terms = solutions.args[0].simplify()

#try randome ones 
#quaternion_terms = Matrix([cos(t/2),sin(t/2),0,0])  #<--- known Gerono
#quaternion_terms = Matrix([cos(t),sin(t)*cos(t),0,0]) #< --- flat bowtie type
#quaternion_terms = Matrix([cos(t/2),sin(t/2)*cos(2*t),0,0]) #< --- funky #1



#quaternion_terms = Matrix([cos(t/2),sin(t/2),0,0])
#quaternion_terms = Matrix([cos(t/2),sin(t/2),0,0])
#quaternion_terms = Matrix([cos(t),1*sin(t)*cos(t),0,0]) #<-- bow tie
#quaternion_terms = Matrix([cos(t),0.5*sin(2*t),0,0]) #<-- bow tie?
#quaternion_terms = Matrix([cos(2*t),sin(t/2)*cos(t/2),0,0]) #<-- cool weird one
#quaternion_terms = Matrix([cos(t),sin(t/2)*cos(t),0,0]) #<-- cool weird one

# confirm unit
isunit = sqrt(quaternion_terms[0]*quaternion_terms[0]+quaternion_terms[1]*quaternion_terms[1]).simplify()
if isunit == 1:
    print('confirmed unit quaterion')
else:
    print('warning: not a unit quaterion')
pprint(quaternion_terms)

#%%

QED = quatrotate(Matrix([quaternion_terms[0],quaternion_terms[1],0,0]),eqn_cir)
#QED = quatrotate(Matrix([quaternion_terms[0],quaternion_terms[1],quaternion_terms[2],0]),eqn_cir)
#pprint(QED)
QED = QED.subs(d, radius)
# plot the 2-d projection
plot_parametric((QED[1], QED[2]),(t,-np.pi,np.pi))
    
# plot the 3-d trajectory
ax = plt.figure().add_subplot(projection='3d')
ax.set_xlim3d([-5, 5])
ax.set_ylim3d([-5, 5])
ax.set_zlim3d([-5, 5])

# initial state (cartesian)
x = np.array([0,radius,0,0])
# for all angles
for i in np.arange(0,2*np.pi,0.1):
#for i in np.arange(-np.pi,np.pi,0.1):
    
    x_prev = x
    eqn_cir_ = eqn_cir.subs({d:radius, t:i})
    quaternion_terms_ = quaternion_terms.subs(t, i)
    pos = quatrotate(np.array([quaternion_terms_[0],quaternion_terms_[1],0,0]),eqn_cir_)
    #pos = quatrotate(np.array([quaternion_terms_[0]*quaternion_terms_[0],quaternion_terms_[1]*quaternion_terms_[1],0,0]),eqn_cir_)
    x = pos
    plt.pause(0.001)
    ax.plot([x_prev[1],x[1]],[x_prev[2], x[2]],[x_prev[3], x[3]],'-b.')
    
#ax.title('Curve')
plt.show()

 
# #%% confirm unit
# isunit = sqrt(solutions.args[0][0]*solutions.args[0][0]+solutions.args[0][1]*solutions.args[0][1]).simplify()
# if isunit == 1:
#     print('confirmed unit quaterion')
# else:
#     print('warning: not a unit quaterion')
    
#%% try by alternate means (TBD)
# Given p2 = q * p * q^(-1), solve for q

# NOTE: this doesn't work yet; there is an error somewhere

alternate = 0

if alternate == 1:
    
    #p, p2, v, costheta, angle, q = symbols('p p2 v costheta angle q', nonzero=True, real=True)

    #rotation axis (3D vector v) is given by the cross product of vectors p and p2, normalized to length 1:
    p = Matrix([eqn_cir[1],eqn_cir[2],eqn_cir[3]])
    p2 = Matrix([LHS[1],LHS[2],LHS[3]]) 
    
    v = p.cross(p2)
    v = v.norm()

    costheta = (p.dot(p2))
    costheta = costheta/(norm_p*norm_v)
    angle = acos(costheta / (p.norm() * p2.norm()))
    q = Matrix([cos(angle/2), v[0]*sin(angle/2), v[1]*sin(angle/2), v[2]*sin(angle/2)])

    pprint(q)


#The rotation angle θ can be calculated from the dot product of the vector parts of p and p2:



#From the rotation axis v and angle θ, we can construct the rotation quaternion q:

#%%



# basic = quatrotate(np.array([a,b,0,0]),eqn_cir)
# pprint(basic)


# proove_bow = Matrix([cos(t),sin(t)*cos(t),0,0]) #<-- bow tie
# proove_bow_ = quatrotate(proove_bow,eqn_cir)
# pprint(proove_bow_)

#%% Solve for a,b directly from the points

# p2 = symbols('p2')

# # this is the desired equations (pure quaternion)
# p2 = Matrix(eqn_8ger[1:4])

# # this is the original equation (pure quaternion)
# p1 = Matrix(eqn_cir[1:4])

# # compute the angle (the norms should be equal to d?, but keep for now)
# theta = p1.dot(p2)/(p1.norm() * p2.norm())
# #theta = p1.dot(p2)/d

# # compute the axis of rotation
# r = p1.cross(p2)
# rhat = r/r.norm()

# # q_direct
# q1 = cos(theta/2).simplify()
# q2 = rhat[0]*sin(theta/2).simplify()
# q3 = rhat[1]*sin(theta/2).simplify()
# q4 = rhat[2]*sin(theta/2).simplify()

# print("a = ")
# pprint(q1)
# print("b = ")
# pprint(q2)

# a1, b1, a2, b2, t, t1, t2 = symbols('a1 b1 a2 b2 t t1, t2', nonzero=True)
# eqn1 = basic.subs({a:a1, b:b1, t:t1})
# eqn2 = basic.subs({a:a2, b:b2, t:t2})
# solve = nonlinsolve(eqn1-eqn2, [t1, t2], S.Reals)
# #bounded_solve = [sol.subs('n', n) for sol in solve for n in range(0, 1)]
# #pprint(bounded_solve)
#%%
#pprint(solve.args[1])
#%%
#test = quatrotate(eqn_8bow,eqn_cir)
#test.simplify()
#test2 = test.subs({t:np.pi})
#pprint(test2)

#%% 
#test = (d*d)*cos(t)*cos(t) + d*d*cos(t)*cos(t)*sin(t)*sin(t) + d*d*sin(t)*sin(t)*sin(t)*sin(t)
#test = cos(t)*cos(t) + cos(t)*sin(t)*sin(t)
#test = p1.cross(p2)
#test = test.norm().simplify()


#pprint(test)

#%% For Later

# ===
# 3 D 
# ===

# # Vivians curve (this looks like my lemni arch?)
# # https://mathworld.wolfram.com/VivianisCurve.html
# # -------------
# w = 0
# x = d*(1+cos(t))
# y = d*sin(t)
# z = 2*d*sin(t/2)
# eqn_viv = Matrix([w,x,y,z])

# # lemniscate of Booth (for later exploration)
# http://www.iosrjournals.org/iosr-jm/papers/Vol6-issue2/G0624353.pdf

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


    
# LEGACY
# # hypotrochoidal curves
# # ----------------------
# psi = (R-r)/r
# w = 0
# x = (R-r)*cos(t)+d*cos(psi*t)
# y = (R-r)*sin(t)-d*sin(psi*t)
# z = 0
# eqn_hyp = Matrix([w,x,y,z])
# #pprint(eqn_hyp)


