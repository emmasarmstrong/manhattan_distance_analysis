#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:29:55 2022

@author: emmaarmstrong
"""

### CALCULATING THE AVERAGE (MEAN) NUMBER OF CALCITE AND ARAGONITE IDENTIFIED IN THE SYSTEM
### PER FRAME, PER ION 
### CLUSTER IS CONSIDERED CALCITE- OR ARAGONITE-LIKE IF THE ION'S STRUCTURE QUANTIFIED WITH
### THE MANHATTAN DISTANCE LIES WITHIN A SET NUMBER OF STANDARD DEVIATIONS OF THE REFERENCE VALUES
### METRIC VARIABLE DETERMINES WHICH METHOD TO COMBINE S_ALPHA AND S_BETA IS USED
### 0 - SEPARATE VALUES
### 1 - TAXICAB
### 2 - EUCLIDEAN
### 3 - CHEBYSHEV

#modules required
from np_manhattan_distance_values import np_manhattan_distances
from np_metrics import taxicab,euclidean,chebyshev
import numpy as np

def np_average_polymorphs(filepaths,nca,nframes,sd,metric): #CV-dist file, number of calcium ions in system, s.d. cutoff and metric used
    #s.d. distances calculated from manhattan distance values in CV_dist
    sd_alpha_cal,sd_beta_cal,sd_alpha_ara,sd_beta_ara = np_manhattan_distances(filepaths,nca,nframes)
    sd_alpha_cal = sd_alpha_cal.flatten()
    sd_beta_cal = sd_beta_cal.flatten()
    sd_alpha_ara = sd_alpha_ara.flatten()
    sd_beta_ara = sd_beta_ara.flatten()
    #frames is actually the total number of timeframes multiplied by the number of Ca ions
    frames = len(sd_alpha_ara)

    #if values lie within the cutoff, counted as calcite or aragonite like in the frame
    if metric == 0:
        sd_cal = len((np.where(np.logical_and(sd_alpha_cal < sd, sd_beta_cal < sd)))[0])
        sd_ara = len((np.where(np.logical_and(sd_alpha_ara < sd, sd_beta_ara < sd)))[0])
    elif metric == 1:
        cal_t = taxicab(sd_alpha_cal,sd_beta_cal)
        ara_t = taxicab(sd_alpha_ara,sd_beta_ara)
        sd_cal = len(cal_t[cal_t<sd])
        sd_ara = len(ara_t[ara_t<sd])
    elif metric == 2:
        cal_e = euclidean(sd_alpha_cal,sd_beta_cal)
        ara_e = euclidean(sd_alpha_ara,sd_beta_ara)
        sd_cal = len(cal_e[cal_e<sd])
        sd_ara = len(ara_e[ara_e<sd])
    elif metric == 3:
        cal_c = chebyshev(sd_alpha_cal, sd_beta_cal)
        ara_c = chebyshev(sd_alpha_ara, sd_beta_ara)
        sd_cal = len(cal_c[cal_c<sd])
        sd_ara = len(ara_c[ara_c<sd])
    #mean values calculated using total number of times clusters are identified as polymorph-like 
    #dived by the number of frames in CVdist 
    ave_cal = sd_cal/frames
    ave_ara = sd_ara/frames
    return ave_cal, ave_ara