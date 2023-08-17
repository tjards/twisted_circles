#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 22:05:40 2023

@author: tjards
"""
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(-2*np.pi, 2*np.pi, 1000)  # Define the range of t values

x = (-np.sqrt(2) * np.sqrt(np.cos(t)**2 + 1)) / 2
y = (-np.sqrt(2) * np.sqrt(-(np.cos(t) - 1)*(np.cos(t) + 1))) / 2

plt.plot(x, y)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Parametric Plot')
plt.grid(True)
plt.show()