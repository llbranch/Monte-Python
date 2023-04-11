#############################
# AESOPE-Lite Monte Carlo
# Created by Liam Branch and Robert Johnson
# Copyright UCSC 2023
#############################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
from tqdm import tqdm
from time import sleep
import math
import scipy
import random
 
#############################
# CONSTANTS
#############################
c = 0.0299792 # Speed of Light in cm / ps
q = 1.60217663e-19 # charge of electron columbs
t_rise = 800 #ps 
T3z=0 #cm is the bottom of T3
T1z=33.782 #cm is the bottom of T1
T4z=-28.07297 #cm is the bottom of T4
T1_radius = 13 #cm
T4_radius = 18 #cm 
xPMT4=9.5*np.cos(110)*2.54
yPMT4=9.5*np.sin(110)*2.54
xPMT1=8.*np.cos(np.radians(-45))*2.54
yPMT1=8.*np.sin(np.radians(-45))*2.54
PMT1_radius = 4.6/2 #cm 
PMT4_radius = 4.6/2 #cm 
# x2PMT1=8.*np.cos(np.radians(-53.72))*2.54 For test
# y2PMT1=8.*np.sin(np.radians(-53.72))*2.54 For test
xPMT4=9.5*np.cos(110)*2.54
yPMT4=9.5*np.sin(110)*2.54
n_dynodes = 8
V = np.linspace(150,850,n_dynodes)
# V = [150,300,350,600,750,850]
n_incident_photons = 10000
E_per_electron = 20
Nmax = 100000
QE = 0.23

#############################
# HELPER FUNCTIONS
#############################

# FIND SIGNIFICANT DIGIT POWER OF 10
def round_to_sig(x):
    return -int(np.floor(np.log10(np.abs(x))))

# NORMALIZE A VECTOR
def normalize(x):
    x /= np.linalg.norm(x)
    return x

# REDEF NORM FOR READABILITY
def mag(x):
    return np.linalg.norm(x)

# DISTANCE 2-DIM CIRCLE WITH LINE SEGMENT
# t = -D . ∆ ± √(D . ∆)^2 - |D|^2(|∆|^2 - R^2)
#     over |D|^2
# ARGUMENTS : # 3d directional vector, 3d point, center of scintillator, radius of scintillator, use corner circle boolean
def distance_circle(u, o, center, radius, quadrant=False): 
    P = o
    D = u*-1 if np.dot(u,P) < 0 else u # does a normalized vector in 3d equate to normalized vector in 2d?
    C = center
    R = radius
    bigDelta = np.array(P)-np.array(C)

    magDsq = mag(D)**2
    magDeltasq = mag(bigDelta)**2
    DdotDelta = np.dot(D,bigDelta)
    if DdotDelta**2 - magDsq * (magDeltasq - R**2) < 0:
        return 100 # some large value that won't be chosen because photon has no intersection with circle
    sqrt_term = np.sqrt(DdotDelta**2 - magDsq * (magDeltasq - R**2))/magDsq
    b_term = -DdotDelta/magDsq
    rootA = b_term - sqrt_term
    rootB = b_term + sqrt_term
    if quadrant is not False: # if in corner don't use the other 3/4ths of the circle to find distance only 4th quadrant part
        return np.abs(rootA) if np.abs(rootA) > np.abs(rootB) else np.abs(rootB)
    return np.abs(rootA) if (rootA < 0) & (np.dot(u,P) < 0) else np.abs(rootB)

# ARGUMENTS : 3d directional vector, 3d point, z positions of planes bottom and top, plane dimension number
def distance_plane(u, o, plane, dim):                                     
    P = o
    if dim==2:
        d_plane = plane[0] if u[dim] < 0 else plane[1]                    # make sure direction matches location of plane 
    else:
        d_plane = plane
    return np.abs((d_plane - P[dim])/u[dim])


# SOLVE FOR DISTANCE LOGIC FUNCTION
def distance_solver(u, o, center, radius, plane_z, corner_center, corner_radius, pmt_center, pmt_radius):
    dcircle = distance_circle(u,o,center,radius)                          # checks distance to circle boundary
    dplane_z = distance_plane(u,o,plane_z,dim=2)                          # checks distance to z boundary in general scint
    dist = dplane_z if dcircle > dplane_z else dcircle
    temp_o = o+dist*u
    PMT_cond = False
    if (temp_o[0] > 0) & (temp_o[1] < 0) & ((temp_o[0]**2+temp_o[1]**2) >= radius**2-1):
        dplanex = distance_plane(u,o,radius,dim=0)                        # checks distance to x boundary
        dplaney = distance_plane(u,o,-radius,dim=1)                       # checks distance to y boundary
        dplanez = distance_plane(u,o,plane_z,dim=2)                       # checks distance to z boundary inside light guide
        dcorner = distance_circle(u,o,corner_center, corner_radius, True) # checks distance to corner boundary
        light_guide_dist = np.min([dplanex,dplaney,dplanez,dcorner])
        temp_o = o+(light_guide_dist)*u                                   # resuse this variable
                                                                          # if close to z = zero and within PMT circle
        if (temp_o[2] < (plane_z[0]+0.01)) & (((temp_o[0]-pmt_center[0])**2+(temp_o[1]-pmt_center[1])**2) <= pmt_radius**2): 
            PMT_cond = True
        return light_guide_dist, PMT_cond
    else:
        return dist, PMT_cond

# PSEUDOCODE FOR EACH PHOTON INTERACTION WITH BOUNDARY
    # if random number X_1 < mean ( Reflectance s_polarization + Reflectance s_polarization ):
        # Reflect
    # else if random number X_2 < absorbption into scintillator boundary probability:
        # Absorbed and exit current particle simulation
    # else if not absorbed:
        # assume photon transmitted through boundary, 
        # absorbed by white paint and reemmitted back 
        # into scintillator with random direction given by random angles Phi_3, Theta_3
        # with constraint of z coordinate entering
def photon_interaction(o, u, n, notabsorbed, scint_plane):
    u_i = u
    u_r = u - 2*np.dot(u, n)*n                      # u_new = u - 2 (u . n)*n
                                                    # CONSTRAINT n_1 <= n_2
    n_1 = 1.000293                                  # Sample index of refraction of air
    n_2 = 1.85                                      # 1.85 for NaI
    v = u*-1 if np.dot(u,n) < 0 else u
    theta = np.arcsin(mag(np.cross(v,n))/(mag(u)*mag(n)))
    inside_sqrt = ((n_1/n_2)*np.sin(theta))**2
    sqrt_term = np.sqrt(1 - inside_sqrt)            # cos(theta)_transmission
    Rs = np.abs((n_1*np.cos(theta) - n_2*sqrt_term)/(n_1*np.cos(theta) + n_2*sqrt_term))**2
    Rp = np.abs((n_1*sqrt_term - n_2*np.cos(theta))/(n_1*sqrt_term + n_2*np.cos(theta)))**2
                                                    # Determine probability of reflectance
    if np.random.random() < ((Rs+Rp)/2):            # if random chance is high enough reflect !
        return normalize(u_r), True                 # return full internal reflection and not absorbed is True
                                                    # else photon is transmitted to white paint
    elif np.random.random() < 0.80:                 # does it get absorbed? change probability when you get more data
        return normalize(u_r), False                # not absorbed is False
    else:                                           # no it didn't get absorbed!
        theta_new = random.uniform(0,2*np.pi)       # new theta direction of photon
        phi_new = random.uniform(0,np.pi)           # new phi   direction of photon
        new_u = normalize(np.array([np.sin(phi_new)*np.cos(theta_new),np.sin(phi_new)*np.sin(theta_new),np.cos(phi_new)]))
        change_factor = np.random.random()-0.5
        u_r = u_r + change_factor*new_u
        return normalize(u_r), True                 # new small change in direction (should be random), and not absorbed is True

# Calculate n vector for all planes and surfaces in apparatus
def n_vec_calculate(o, scint_plane, light_guide_planes, corner_center, corner_radius):
    if o[2] == scint_plane[0]:                                      # bottom of scint
        return np.array([0,0,+1])
    elif o[2] == scint_plane[1]:                                    # top of scint
        return np.array([0,0,-1])
    elif o[0] == light_guide_planes[0]:                             # y plane of light guide 
        return np.array([0,+1,0])
    elif o[1] == light_guide_planes[1]:                             # x plane of light guide
        return np.array([-1,0,0])
    elif (o[0] >= corner_center[0]) & (o[1] <= corner_center[1]):   # in corner
        return normalize(o-corner_center)
    else:                                                           # in main scintillator
        return normalize(o-np.array([0,0,0]))


#############################
# SIMULATION FUNCTIONS
#############################

# PSEUDOCODE FOR PARTICLE GENERATION
    # Generate random position in circle and random direction in allowed cone
    # Walk
    # while z of particle > lowest z point of T4
    #     if point is outside of scintillator
    #          then step to next scintillator boundary
    #     for each scintillator:
    #          if point is insde of scintillator_i
    #                Generate photons if random number X_1 < Pr(scintillate)
    #                walk mean free path length and store position if still within scintillator

def particle_path(t, phi_range_deg, T1_z, T1_width, T4_z, T4_width, T1_radius, T4_radius, mean_free_path, photons_per_E, prob_scint):
    theta = random.uniform(0,2*np.pi)                                                 # random theta in circle above T1
    phi = random.uniform(np.pi-phi_range_deg*np.pi/180/2,np.pi+phi_range_deg*np.pi/180/2) # phi angle pointing in -k given phi range
    maxdist = np.random.random()*T1_radius/2                                          # half the radius of T1
    round_const = round_to_sig(mean_free_path)
    o = np.float64((maxdist*np.cos(theta), maxdist*np.sin(theta), T1_z+T1_width+2))   # x, y, top of T1_z+2
    u = np.array((np.cos(theta)*np.sin(phi),np.sin(theta)*np.sin(phi),np.cos(phi)),dtype=np.float64)
    print(f"u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
    photons = [0]                                                                     # begin photon array
    points = [o]                                                                      # current point 
    times = [t]                                                                       # current t 
    z = points[-1][2]                                                                 # current z 
    z_1 = (z+mean_free_path*u[2]).round(round_const)                                  # next z step
    inside_scint = False
    while z_1 >= T4_z:
        if not inside_scint:
            distT1 = np.abs((T1_z+T1_width - z)/u[2])
            distT4 = np.abs((T4_z+T4_width - z)/u[2])
            dist = distT4 if z_1 < T1_z else distT1
            t +=  dist/c                                                               # calculate time in ps passed
            times.append(t)
            points.append(points[-1] + dist*u)
            phot = np.random.poisson(photons_per_E)
            if np.random.random() < prob_scint: photons.append(phot)
            else: photons.append(0)
            z = points[-1][2]
            z_1 = (z+mean_free_path*u[2]).round(round_const)
        for Tbottom,Ttop in [(T1_z,T1_z+T1_width),(T4_z,T4_z+T4_width)]:
            inside_scint = (z_1 <= (Ttop)) & (z_1 >= Tbottom)
            if inside_scint:
                while inside_scint:
                    t +=  mean_free_path/c
                    times.append(t)
                    points.append(points[-1] + mean_free_path*u)
                    phot = np.random.poisson(photons_per_E)
                    if np.random.random() < prob_scint: photons.append(phot)
                    else: photons.append(0)
                    z = points[-1][2]
                    z_1 = (z+mean_free_path*u[2]).round(round_const)
                    inside_scint = (z_1 <= (Ttop)) & (z_1 >= Tbottom)
    print(f"time elapsed in scintillators: {np.abs(times[1]-times[-1]):.2f}ps total scintillation points: {len(points[1:])}")
    return np.array(times, dtype=np.float64)[1:], np.array(points, dtype=np.float64)[1:], np.array(photons[1:], dtype=np.float64)


def scintillator_monte_carlo(o, notabsorbed, scint_radius, scint_plane, scint_width, light_guide_planes, pmt_center, pmt_radius, N_max, dt):
    track_history = np.zeros((N_max+1,7))         # x, y history of Photon
    corner_radius = scint_width*4
    corner_center = [scint_radius-corner_radius,-scint_radius+corner_radius,scint_plane[0]]
    theta = random.uniform(0,2*np.pi)             # first theta direction of photon
    phi = random.uniform(0,np.pi)                 # first phi   direction of photon
    PMT_hit_condition = False
    u = np.array([np.sin(phi)*np.cos(theta),np.sin(phi)*np.sin(theta),np.cos(phi)]) # first direction unit vector
    track_history[0,:] = [o[0],o[1],o[2],u[0],u[1],u[2],notabsorbed]
    i = 1
    while (i < N_max+1) & (not PMT_hit_condition) & (notabsorbed is True):
        ds, PMT_hit_condition = distance_solver(u, o, np.array([0,0,scint_plane[0]]),scint_radius, scint_plane, corner_center, corner_radius, pmt_center, pmt_radius)
        x, y, z = o+ds*u
        o = np.array([x, y, np.abs(z) if np.abs(z-scint_plane).any() < 1e-5 else z])
        dt += np.abs(ds)/c                        # time taken in ps traveling in direction theta
#         print(f"step {i}: ds={ds:.2f}cm dt={dt:.2f}ps Absorbed?={not notabsorbed} xyz =({x:.2f},{y:.2f},{z:.2f}) u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
        n = n_vec_calculate(o, scint_plane, light_guide_planes, corner_center, corner_radius)
        u, notabsorbed = photon_interaction(o, u, n, notabsorbed, scint_plane)
        track_history[i] = [x,y,z,u[0],u[1],u[2],notabsorbed]
        i+=1
    if i < N_max+1:
        track_history = track_history[:i,:]
    return PMT_hit_condition, dt, track_history

# PMT SIMULATION
def photoElectrons(photons, QE): # Main monte carlo
    pe = 0. 
    for i in range(int(photons)):
        if np.random.random()<QE:
            pe+=1
    return pe
def photontoElectrons(photons, V, QE, N, E_per_electron):
    e = photoElectrons(photons, QE)
    for dynode in range(N-1):
        delta_voltage = V[dynode+1]-V[dynode]
        e += np.random.poisson(e*delta_voltage/E_per_electron)
    return e

#############################
# RUN SIMULATION 
#############################

# FIND PARTICLE PATH
times, points, photons = particle_path(t=0, phi_range_deg=40, T1_z=T1z, T1_width=0.5, 
                                       T4_z=T4z, T4_width=1, T1_radius=13, T4_radius=18, mean_free_path=0.00024, 
                                       photons_per_E=10, prob_scint=0.8)
N = np.sum(photons)
print("Photons generated", N)

# SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
T1_input_times = []
T4_input_times = []
pmt_hits = 0
with tqdm(total=N, file=sys.stdout) as pbar:
    for i,(point_i,time_i) in enumerate(zip(points,times)):
        # IF IN T1
        if point_i[2] >= T1z:
            for photon in range(photons[i].astype(int)):
                hit_PMT, travel_time, check = scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=T1_radius, 
                                                  scint_plane=np.array([T1z,T1z+0.5]), scint_width=0.5, 
                                                  light_guide_planes=[T1_radius,-T1_radius], 
                                                  pmt_center=[T1_radius-4*0.5,-T1_radius+4*0.5,T1z], pmt_radius=PMT1_radius,
                                                  N_max=8, dt=time_i)
                if hit_PMT: 
                    T1_input_times.append(travel_time)
                    pmt_hits +=1
                pbar.set_description(f'hits: {pmt_hits}')
                pbar.update(1)
        else:
        # ELSE IN T4
            for photon in range(photons[i].astype(int)):
                hit_PMT, travel_time, _ = scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=T4_radius, 
                                                  scint_plane=np.array([T4z,T4z+1]), scint_width=1, 
                                                  light_guide_planes=[T4_radius,-T4_radius], 
                                                  pmt_center=[T4_radius-4*0.5,-T4_radius+4*0.5,T4z], pmt_radius=PMT4_radius,
                                                  N_max=8, dt=time_i)
                if hit_PMT: 
                    T4_input_times.append(travel_time)
                    pmt_hits +=1
                pbar.set_description(f'hits: {pmt_hits}')
                pbar.update(1)
# PRINT RESULTS
print("HITS on T1",len(T1_input_times),"\n",T1_input_times)
print("HITS on T4",len(T4_input_times),"\n",T4_input_times)

# BEGIN SIMULATING PMT PULSE
pmt_electron_travel_time = 16000 # approx 16 ns
signals = []
output_times = []
QE = 1 # for testing purposes
for i,t in enumerate(T1_input_times):
    pmtSignal_i = photontoElectrons(1, V, QE, n_dynodes, E_per_electron)
    output_times.append(pmt_electron_travel_time+t)
    signals.append(pmtSignal_i)
for i,t in enumerate(T4_input_times):
    pmtSignal_i = photontoElectrons(1, V, QE, n_dynodes, E_per_electron)
    output_times.append(pmt_electron_travel_time+t)
    signals.append(pmtSignal_i)

# CONVERTION Electron count to Current
signals *= q / t_rise 

# OUTPUT FORMATTING
print("OUTPUT TIMES", output_times)
print("SIGNALS", signals)
print("Exporing...", end='')
pd.DataFrame(rows=[output_times,signals], columns=['time','current']).to_csv('monte_carlo_output.txt')
print("Done!")