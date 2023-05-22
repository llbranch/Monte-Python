"""
    Used this script to fix photon_interaction and re-emmition angle in the main class file
    Copyright Liam Branch, Michelle Pico..., Robert Johnson
    UCSC 2023
"""
import numpy as np
import matplotlib.pyplot as plt
import random
# import pandas as pd
# from scipy.signal import find_peaks, peak_prominences, peak_widths
# import sys

def normalize(x):
    x /= np.linalg.norm(x)
    return x

x = normalize([-1,0,0]) # x n-vec
y = normalize([0,1,0]) # y n-vec
z_top = normalize([0,0,-1]) # z_top n-vec
z_bot = normalize([0,0,1]) # z_bot n-vec
n = z_bot

new_u = np.zeros((100,3))
for i in range(len(new_u)):
    theta_new = random.uniform(-np.pi/2,np.pi/2)       # new theta direction of photon TODO: generate cos theta uniform distribution
    phi_new = random.uniform(-np.pi, np.pi)           # new phi   direction of photon
    new_u[i] = normalize(np.array([np.sin(phi_new)*np.cos(theta_new),np.sin(phi_new)*np.sin(theta_new),np.cos(phi_new)]))
    change_factor = 0.95
    new_u[i] = normalize(change_factor*new_u[i] + n)
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111,projection='3d')
ax.quiver(1,1,1,n[0],n[1],n[2], color='k', length=0.07)
ax.quiver(np.ones(100),np.ones(100),np.ones(100),new_u[:,0],new_u[:,1],new_u[:,2], alpha=0.5, length=0.05)
plt.show()