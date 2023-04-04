#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 15:19:19 2022

@author: emmaarmstrong
"""
import numpy as np

def np_manhattan_distances(filepaths,nca,nframes):  #CV-dist file, number of calcium ions in system and number of frames
    
    #mean values of s_alpha and s_beta for calcite and aragonite
    alpha_cal_mean = 0.19
    beta_cal_mean = 0.15
    alpha_ara_mean = 0.17
    beta_ara_mean = 0.13
    
    #s.d. values for s_alpha and s_beta for calcite and aragonite 
    alpha_cal_sd = 0.062
    beta_cal_sd = 0.056
    alpha_ara_sd = 0.055
    beta_ara_sd = 0.045

    nfiles = len(filepaths)
    
    s_alpha_cal = np.zeros([int(nfiles*nframes-nfiles+1),nca]) #empty list for Ca-C MD value with reference to calcite
    s_beta_cal = np.zeros([int(nfiles*nframes-nfiles+1),nca]) #empty list for C-Ca-C MD value with reference to calcite
    s_alpha_ara = np.zeros([int(nfiles*nframes-nfiles+1),nca]) #empty list for Ca-C MD value with reference to aragonite
    s_beta_ara = np.zeros([int(nfiles*nframes-nfiles+1),nca]) #empty list for C-Ca-C MD value with reference to aragonite
    
    #reading in s_alpha and s_beta values from CV-dist file for each Ca ion in system
    file = 0
    frame = 0
    for f in filepaths:
        file += 1
        #open single file and extracting Manhattan distance values
        cvdist = open(f,"r")
        for line in cvdist:
            ls = line.split()
            if ls[1] == 'Frame': #identifying frames in files
                frame += 1
                #finding manhattan distance values for the given number of Ca ions
                if frame/file != nframes or frame == nfiles*nframes:
                    for i in range(nca): #not including the last frame to avoid repetition of values                     
                        s = next(cvdist).split()
                    
                        s_alpha_ara[frame-1][i] = float(s[5])
                        s_alpha_cal[frame-1][i] = float(s[6])
                        s_beta_ara[frame-1][i] = float(s[7])
                        s_beta_cal[frame-1][i] = float(s[8])
    
    sd_alpha_cal = np.abs(s_alpha_cal - alpha_cal_mean)/alpha_cal_sd
    sd_beta_cal = np.abs(s_beta_cal - beta_cal_mean)/beta_cal_sd
    sd_alpha_ara = np.abs(s_alpha_ara - alpha_ara_mean)/alpha_ara_sd
    sd_beta_ara = np.abs(s_beta_ara - beta_ara_mean)/beta_ara_sd
            
    return sd_alpha_cal,sd_beta_cal,sd_alpha_ara,sd_beta_ara