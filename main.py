#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 07:35:05 2023

@author: tjards
"""

# import stuff
from sympy import *

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
    
    

# define quatertions we want to manipulate
# ---------------------------------------


var('a b d r R t psi t z', positive = True)

# parametric equations for hypotrochoidal curves
psi = (R-r)/r
w = 0
x = (R-r)*cos(t)+d*cos(psi*t)
y = (R-r)*sin(t)-d*sin(psi*t)
z = 0
eqn_hyp = Matrix([w,x,y,z])
#pprint(eqn_hyp)

# figure-8
# --------
w = 0
x = d*cos(t)
y = d*sin(t)*cos(t)
z = 0.5*d*sin(t)*sin(t)
eqn_8 = Matrix([w,x,y,z])
pprint(eqn_8)

# equation of circle with radius d
w = 0
x = d*cos(t)
y = d*sin(t)
z = 0
eqn_cir = Matrix([w,x,y,z])
pprint(eqn_cir)

# unknowns
u1 = Matrix([a,b,0,0])
u2 = quatjugate(u1)
pprint(u1)
#pprint(u2)




'''
Unnecessary complicated thing:
Let us solve for a, b that generate a hypotrochoidal curves as
    a function of projection of quaternion rotations around x,y

'''
#%%
RHS = quatrotate(u1,eqn_cir)
LHS = eqn_8
pprint(LHS)
print('=')
pprint(RHS)

#EQNS = LHS-RHS
#equations = Eq(EQNS[1],EQNS[2],EQNS[3])
#solutions = solve(equations, a)
#pprint(solutions)







# test examples
# --------
# var('w x y z')
# q1 = Matrix([w, x, y, z])
# var('a b c d')
# q2 = Matrix([a, b, c, d])
# product = hamilton_product(q1, q2)
# pprint(product)


    



