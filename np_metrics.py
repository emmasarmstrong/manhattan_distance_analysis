#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:16:56 2022

@author: emmaarmstrong
"""

### METRICS USED TO COMBINE S_ALPHA AND S_BETA INTO A SINGLE MANHATTAN DISTANCE VALUE
### TAXICAB, EUCLIDEAN AND CHEBYSHEV METRICS ARE USED

from math import sqrt
import numpy as np

def taxicab(a,b): #list or array of values to be combined
    return a+b

def euclidean(a,b): #list or array of values to be combined
    return np.sqrt(a**2 + b**2)

def chebyshev(a,b): #list or array of values to be combined
    return np.maximum(a,b)