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


# define quatertions we want to manipulate
# ---------------------------------------
# define our variables
var('w x y z')
q1 = Matrix([w, x, y, z])

var('a b c d')
q2 = Matrix([a, b, c, d])

product = hamilton_product(q1, q2)


pprint(product) # pprint is "pretty print"


    



