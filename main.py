#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 07:35:05 2023

@author: tjards
"""

from sympy import *

# define quatertions we want to manipulate
# ---------------------------------------

#q1 = np.array([w, x, y, z])





a       = Symbol('x')
x,y,z   = symbols('p,q,r')


print(a + a + 1)
print(x)




# multiplication between quaternions
# -----------------------------------
# def quat_mult(q1, q2):
#     w1, x1, y1, z1 = q1
#     w2, x2, y2, z2 = q2
#     w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
#     x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
#     y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
#     z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
#     return np.array([w, x, y, z])