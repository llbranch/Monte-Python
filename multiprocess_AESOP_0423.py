#############################
# AESOP-Lite Monte Carlo
# Multiprocess Edited Version 1 (Better)
# Created by Liam Branch and Robert Johnson
# Copyright UCSC 2023
#############################

import numpy as np
import pandas as pd
from tqdm import tqdm
import random
from time import perf_counter
from datetime import timedelta, datetime
from multiprocessing import Pool, cpu_count, freeze_support
 
class Simulation:
    def __init__(self):
        #############################
        # CONSTANTS
        #############################
        self.c = 0.0299792 # Speed of Light in cm / ps
        self.q = 1.60217663e-19 # charge of electron columbs
        # CONSTRAINT n_1 <= n_2
        self.n_1 = 1.000293 # Sample index of refraction of air
        self.n_2 = 1.58 # 1.58 for EJ-200
        self.T3z=0 #cm is the bottom of T3
        self.T1z=33.782 #cm is the bottom of T1
        self.T4z=-28.07297 #cm is the bottom of T4
        self.T1_radius = 13 #cm
        self.T4_radius = 18 #cm 
        self.T1_width = 0.5 #cm
        self.T4_width = 1 #cm
        self.xPMT4=9.5*np.cos(np.radians(110))*2.54
        self.yPMT4=9.5*np.sin(np.radians(110))*2.54
        self.xPMT1=8.*np.cos(np.radians(-45))*2.54 # x2PMT1=8.*np.cos(np.radians(-53.72))*2.54 For test
        self.yPMT1=8.*np.sin(np.radians(-45))*2.54 # y2PMT1=8.*np.sin(np.radians(-53.72))*2.54 For test
        self.PMT1_radius = 4.6/2 #cm need to change this to 46 milimeters or 0.046 cm
        self.PMT4_radius = 4.6/2 #cm 
        
        self.n_dynodes = 8
        self.V = np.linspace(150,850,self.n_dynodes)
        # self.V = [150,300,350,600,750,850]
        self.E_per_electron = 20
        self.QE = 1 #0.23
        self.t_initial = 0 #ps
        self.particle_init_angle_range = 40 #degrees
        self.particle_gen_area = self.T1_radius
        self.particle_gen_z = self.T1z+self.T1_width + 2 #cm
        self.mean_free_path_scints = 8e-5 #cm or 80 micro meters
        self.photons_produced_per_MeV = 10 # true value is closer to 10000 per 1MeV
        self.pr_of_scintillation = 0.8
        self.max_simulated_reflections = 40
        self.pmt_electron_travel_time = 0 # approx 16 ns
        self.artificial_gain = 1 # gain factor
        self.pr_absorption = 0.1 # probability of boundary absorbing
        self.seperation_time = 1e5 # ps 
        self.output_bin_width = 100 # ps
        self.num_particles = 1 # default Muons
        self.CMOS_thresh = 1.5 # V
        self.reemission_angle_factor = 0.9 # range [0,1] --> cone from [-pi,pi]
        
        # Introduction Print Statement
        print("######################################################")
        print("Generated Apparatus Simulation with following defaults")
        print("######################################################")
        print("PARTICLE: Mean Free Path =", self.mean_free_path_scints, "cm")
        print("PARTICLE: Time Seperation between sequential Particles if simulation more than 1 =",self.seperation_time)
        print("SCINT:    Probability of Scintillaton =", self.pr_of_scintillation)
        print("PMT:      Quantum Efficiency is set to", self.QE, "by default to keep more pulses")
        print("PMT:      Energy per Photoelectron is set to", self.E_per_electron, "by best estimation")
        print("PMT:      Artificial Gain on Output Current =", self.artificial_gain)
        print("OUTPUT:   Binning Width for PWL output file =", self.output_bin_width, "ps")
        print("\nRun with .run() function given optional arguments below")
        print("integer n particles, 'delta_t' =", self.seperation_time, "ps particle time seperation")
        
    #############################
    # HELPER FUNCTIONS
    #############################

    # FIND SIGNIFICANT DIGIT POWER OF 10
    def round_to_sig(self, x):
        return -int(np.floor(np.log10(np.abs(x))))

    # NORMALIZE A VECTOR
    def normalize(self, x):
        x /= np.linalg.norm(x)
        return x

    # REDEF NORM FOR READABILITY
    def mag(self, x):
        return np.linalg.norm(x)
    
    # LIGHT GUIDE CONDITION
    def lg_condition(self, corner_pt, scint_corner, scint_num):
        ret = (corner_pt[0] > 0) & (corner_pt[0] < scint_corner[0]) & (corner_pt[1] < 0) & (corner_pt[1] > scint_corner[1])
        if scint_num == 4:
            ret = (corner_pt[0] > 0) & (corner_pt[0] < scint_corner[0]) & (corner_pt[1] < 0) & (corner_pt[1] > scint_corner[1])
        return ret

    # SCINT RADIUS CONDITION
    def scint_condition(self, corner_pt, scint_radius, scint_num):
        ret = np.sqrt(np.sum(corner_pt[0:2]**2)) < self.T1_radius
        if scint_num == 4:
            ret = np.sqrt(np.sum(corner_pt[0:2]**2)) < self.T4_radius
        return ret

    # DISTANCE 2-DIM CIRCLE WITH LINE SEGMENT
    # t = -D . ∆ ± √(D . ∆)^2 - |D|^2(|∆|^2 - R^2)
    #     over |D|^2
    # ARGUMENTS : # 3d directional vector, 3d point, center of scintillator, radius of scintillator, use corner circle boolean
    def distance_circle(self, u, o, center, radius, quadrant=False): 
        cond = np.dot(u,o) < 0
        D = u*-1 if cond else u # does a normalized vector in 3d equate to normalized vector in 2d?
        bigDelta = np.array(o)-np.array(center)
        magDsq = self.mag(D)**2
        magDeltasq = self.mag(bigDelta)**2
        DdotDelta = np.dot(D,bigDelta)
        if DdotDelta**2 - magDsq * (magDeltasq - radius**2) < 0:
            return 100 # some large value that won't be chosen because photon has no intersection with circle
        sqrt_term = np.sqrt(DdotDelta**2 - magDsq * (magDeltasq - radius**2))/magDsq
        b_term = -DdotDelta/magDsq
        rootA = b_term - sqrt_term
        rootB = b_term + sqrt_term
        if quadrant is not False: # if in corner don't use the other 3/4ths of the circle to find distance only 2nd or 4th quadrant part
            return np.abs(rootA) if np.abs(rootA) > np.abs(rootB) else np.abs(rootB)
        return np.abs(rootA) if (rootA < 0) & cond else np.abs(rootB)

    # ARGUMENTS : 3d directional vector, 3d point, z positions of planes bottom and top, plane dimension number
    def distance_plane(self, u, o, plane, dim):                                     
        if dim==2:
            d_plane = plane[0] if u[dim] < 0 else plane[1]                    # make sure direction matches location of plane 
        else:
            d_plane = plane
        return np.abs((d_plane - o[dim])/u[dim])


    # SOLVE FOR DISTANCE LOGIC FUNCTION
    def distance_solver(self, u, o, center, radius, plane_z, corner_center, corner_radius, pmt_center, pmt_radius):
        dcircle = self.distance_circle(u,o,center,radius)                          # checks distance to circle boundary
        dplane_z = self.distance_plane(u,o,plane_z,dim=2)                          # checks distance to z boundary in general scint
        dist = dplane_z if dcircle > dplane_z else dcircle
        temp_o = o+dist*u
        PMT_cond = False
        if (pmt_center[0] > 0) & (temp_o[0] > 0) & (temp_o[1] < 0) & ((temp_o[0]**2+temp_o[1]**2) >= radius**2-1):
            dplanex = self.distance_plane(u,o,radius,dim=0)                        # checks distance to +x boundary
            dplaney = self.distance_plane(u,o,-radius,dim=1)                       # checks distance to -y boundary
            dplanez = self.distance_plane(u,o,plane_z,dim=2)                       # checks distance to z boundary inside light guide
            dcorner = self.distance_circle(u,o,corner_center, corner_radius, True) # checks distance to corner boundary
            light_guide_dist = np.min([dplanex,dplaney,dplanez,dcorner])
            temp_o = o+(light_guide_dist)*u                                   # resuse this variable
                                                                            # if close to z = zero and within PMT circle
            if (temp_o[2] < (plane_z[0]+0.01)) & (((temp_o[0]-pmt_center[0])**2+(temp_o[1]-pmt_center[1])**2) <= pmt_radius**2): 
                PMT_cond = True
            return light_guide_dist, PMT_cond
        elif (pmt_center[0] < 0) & (temp_o[0] < 0) & (temp_o[1] > 0) & ((temp_o[0]**2+temp_o[1]**2) >= radius**2-1):
            dplanex = self.distance_plane(u,o,-radius,dim=0)                       # checks distance to -x boundary
            dplaney = self.distance_plane(u,o,radius,dim=1)                        # checks distance to +y boundary
            dplanez = self.distance_plane(u,o,plane_z,dim=2)                       # checks distance to z boundary inside light guide
            dcorner = self.distance_circle(u,o,corner_center, corner_radius, True) # checks distance to corner boundary
            light_guide_dist = np.min([dplanex,dplaney,dplanez,dcorner])
            temp_o = o+(light_guide_dist)*u                                   # resuse this variable
                                                                            # if close to z = zero and within PMT circle
            if (temp_o[2] < (plane_z[0]+0.01)) & (((temp_o[0]-pmt_center[0])**2+(temp_o[1]-pmt_center[1])**2) <= pmt_radius**2): 
                PMT_cond = True
            return light_guide_dist, PMT_cond
        else:
            return dist, PMT_cond

    # PSEUDOCODE FOR EACH PHOTON INTERACTION WITH BOUNDARY
        # if random number X_1 < mean ( Reflectance s_polarization + Reflectance p_polarization ):
            # Reflect
        # else if random number X_2 < absorbption into scintillator boundary probability:
            # Absorbed and exit current particle simulation
        # else if not absorbed:
            # assume photon transmitted through boundary, 
            # absorbed by white paint and reemmitted back 
            # into scintillator with random direction given by random angles Phi_3, Theta_3
            # with constraint of z coordinate entering
    def photon_interaction(self, u, n):
        u_r = u - 2*np.dot(u, n)*n                              # u_new = u - 2 (u . n)*n
        v = u*-1 if np.dot(u,n) < 0 else u
        theta = np.arcsin(self.mag(np.cross(v,n))/(self.mag(u)*self.mag(n)))
        inside_sqrt = ((self.n_1/self.n_2)*np.sin(theta))**2
        sqrt_term = np.sqrt(1 - inside_sqrt)                    # cos(theta)_transmission
        Rs = np.abs((self.n_1*np.cos(theta) - self.n_2*sqrt_term)/(self.n_1*np.cos(theta) + self.n_2*sqrt_term))**2
        Rp = np.abs((self.n_1*sqrt_term - self.n_2*np.cos(theta))/(self.n_1*sqrt_term + self.n_2*np.cos(theta)))**2
                                                                # Determine probability of reflectance
        if np.random.random() < ((Rs+Rp)/2):                    # if random chance is high enough reflect !
            return self.normalize(u_r), True                        # return full internal reflection and not absorbed is True
                                                                # else photon is transmitted to white paint
        elif np.random.random() < self.pr_absorption:               # does it get absorbed? change probability when you get more data
            return self.normalize(u_r), False                       # not absorbed is False
        else:                                                   # no it didn't get absorbed!
            theta_new = random.uniform(-np.pi/2,np.pi/2)            # new theta direction of photon
            phi_new = random.uniform(-np.pi, np.pi)                 # new phi   direction of photon
            new_u = self.normalize(np.array([np.sin(phi_new)*np.cos(theta_new),np.sin(phi_new)*np.sin(theta_new),np.cos(phi_new)]))
            u_r = self.reemission_angle_factor*new_u + n
            return self.normalize(u_r), True                        # new small change in direction (should be random), and not absorbed is True

    # Calculate n vector for all planes and surfaces in apparatus
    def n_vec_calculate(self, o, scint_plane, light_guide_planes, corner_center, corner_radius):
        if o[2] == scint_plane[0]:                                      # bottom of scint
            return np.array([0,0,+1])
        elif o[2] == scint_plane[1]:                                    # top of scint
            return np.array([0,0,-1])
        elif o[0] == light_guide_planes[0]:                             # y plane of light guide 
            return np.array([0,light_guide_planes[0]/abs(light_guide_planes[0]),0])
        elif o[1] == light_guide_planes[1]:                             # x plane of light guide
            return np.array([light_guide_planes[1]/abs(light_guide_planes[1]),0,0])
        elif (o[0] >= corner_center[0]) & (o[1] <= corner_center[1]):   # in corner
            return self.normalize(o-corner_center)
        else:                                                           # in main scintillator
            return self.normalize(o-np.array([0,0,0]))


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

    def particle_path(self, t, phi_range_deg, T1_z, T1_width, T4_z, T4_width, T1_radius, T4_radius, T1_corner, T4_corner, mean_free_path, photons_per_E, prob_scint):
        theta = random.uniform(0,2*np.pi)                                                     # random theta in circle above T1
        phi = random.uniform(np.pi-phi_range_deg*np.pi/180/2,np.pi+phi_range_deg*np.pi/180/2) # phi angle pointing in -k given phi range
        maxdist = np.random.random()*self.particle_gen_area                                   # radius of generation
        round_const = self.round_to_sig(mean_free_path)
        o = np.float64((maxdist*np.cos(theta), maxdist*np.sin(theta), self.particle_gen_z))   # x, y, z of new particle
        u = np.array((np.cos(theta)*np.sin(phi),np.sin(theta)*np.sin(phi),np.cos(phi)),dtype=np.float64)
        # print(f"u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
        photons = [0]                                                                         # begin photon array
        points = [o]                                                                          # current point 
        times = [t]                                                                           # current t 
        cur_o = points[-1]                                                                    # current z 
        next_o = (cur_o+mean_free_path*u).round(round_const)                                  # next z step
        inside_scint = False
        missed = 0
        while next_o[2] >= T4_z:
            if not inside_scint:
                if missed:
                    theta = random.uniform(0,2*np.pi)                                                 # reset random theta in circle above T1
                    phi = random.uniform(np.pi-phi_range_deg*np.pi/180/2,np.pi+phi_range_deg*np.pi/180/2) # reset phi angle pointing in -k given phi range
                    maxdist = np.random.random()*T1_radius/2                                          # reset random point inside half the radius of T1
                    round_const = self.round_to_sig(mean_free_path)
                    o = np.float64((maxdist*np.cos(theta), maxdist*np.sin(theta), T1_z+T1_width+2))   # reset x, y, top of T1_z+2
                    u = np.array((np.cos(theta)*np.sin(phi),np.sin(theta)*np.sin(phi),np.cos(phi)),dtype=np.float64) # reset u direction
                    photons.clear(); points.clear(); times.clear()
                    photons = [0]                                                                     # reset photon array
                    points = [o]                                                                      # reset current point 
                    times = [t]                                                                       # reset current t 
                    cur_o = points[-1]                                                                # reset current z 
                    next_o = (cur_o+mean_free_path*u).round(round_const)                              # reset next z step
                    missed = False
                distT1 = np.abs((T1_z+T1_width - cur_o[2])/u[2])
                distT4 = np.abs((T4_z+T4_width - cur_o[2])/u[2])
                dist = distT4 if next_o[2] < T1_z else distT1
                check = (cur_o+dist*u).round(round_const)
                inside_T1 = self.scint_condition(check, T1_radius, 1) | self.lg_condition(check, T1_corner, 1)
                inside_T4 = self.scint_condition(check, T4_radius, 4) | self.lg_condition(check, T4_corner, 4)
                scint_cond = inside_T4 if check[2] < T1_z else inside_T1
                # print(f"inside_T1={inside_T1} inside_T4={inside_T4}")
                # print("outer whileloop", scint_cond, next_o, dist, T4_z)          
                if scint_cond:
                    t +=  dist/self.c                                              # calculate time in ps passed
                    times.append(t)
                    points.append(points[-1]+dist*u+mean_free_path*u)
                    phot = np.random.poisson(photons_per_E)
                    if np.random.random() < prob_scint: photons.append(phot)
                    else: photons.append(0)
                    cur_o = points[-1]                                             # current point 
                    next_o = (cur_o+mean_free_path*u).round(round_const)           # next point
                    # print("z",cur_o[2],"z_1",next_o[2])
                    inside_scint = True
                else:                                                              # missed a scintillator / lightguide so throw away and restart
                    # print("missed!")
                    missed = True
                    inside_scint = False
                    continue
            for Tbottom,Ttop,Tradius,Tcorner,num in [(T1_z,T1_z+T1_width,T1_radius,T1_corner,1),(T4_z,T4_z+T4_width,T4_radius,T4_corner,4)]:
                inside_scint = (next_o[2] <= (Ttop)) & (next_o[2] >= Tbottom) & (self.scint_condition(next_o, Tradius, num) | self.lg_condition(next_o, Tcorner, num))
                while inside_scint:
                    # print("inner whileloop", inside_scint)
                    t += mean_free_path/self.c
                    times.append(t)
                    points.append(cur_o+mean_free_path*u)
                    phot = np.random.poisson(photons_per_E)
                    if np.random.random() < prob_scint: photons.append(phot)
                    else: photons.append(0)
                    cur_o = points[-1]                                             # current point 
                    next_o = (cur_o+mean_free_path*u).round(round_const)           # next point
                    inside_scint = (next_o[2] <= (Ttop)) & (next_o[2] >= Tbottom) & (self.scint_condition(next_o, Tradius, num) | self.lg_condition(next_o, Tcorner, num))
        # print(f"time elapsed in scintillators: {np.abs(times[1]-times[-1]):.2f}ps total scintillation points: {len(points[1:])}")
        return np.array(times, dtype=np.float64)[1:], np.array(points, dtype=np.float64)[1:], np.array(photons[1:], dtype=np.float64)


    def scintillator_monte_carlo(self, o, notabsorbed, scint_radius, scint_plane, scint_width, light_guide_planes, pmt_center, pmt_radius, N_max, t, keepdata):
        if keepdata: track_history = np.zeros((N_max+1,7))         # x, y history of Photon
        corner_radius = scint_width*4
        signs = light_guide_planes/np.abs(light_guide_planes)
        corner_center = [signs[0]*(scint_radius-corner_radius),signs[1]*(scint_radius-corner_radius),scint_plane[0]]
        endpoint_dist = self.mag(o-pmt_center)
        theta = random.uniform(0,2*np.pi)             # first theta direction of photon
        phi = random.uniform(0,np.pi)                 # first phi   direction of photon
        PMT_hit_condition = False
        total_dist = 0
        u = np.array([np.sin(phi)*np.cos(theta),np.sin(phi)*np.sin(theta),np.cos(phi)]) # first direction unit vector
        if keepdata: track_history[0,:] = [o[0],o[1],o[2],u[0],u[1],u[2],notabsorbed]
        i = 0
        while (i < N_max) & (not PMT_hit_condition) & (notabsorbed is True):
            ds, PMT_hit_condition = self.distance_solver(u, o, np.array([0,0,scint_plane[0]]),scint_radius, scint_plane, corner_center, corner_radius, pmt_center, pmt_radius)
            x, y, z = o+ds*u
            total_dist += self.mag(ds*u[0:2])
            o = np.array([x, y, np.abs(z) if np.abs(z-scint_plane).any() < 1e-5 else z])
            t += np.abs(ds)/self.c                        # time taken in ps traveling in direction theta
    #         print(f"step {i}: ds={ds:.2f}cm dt={dt:.2f}ps Absorbed?={not notabsorbed} xyz =({x:.2f},{y:.2f},{z:.2f}) u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
            n = self.n_vec_calculate(o, scint_plane, light_guide_planes, corner_center, corner_radius)
            u, notabsorbed = self.photon_interaction(u, n)
            if keepdata: track_history[i+1] = [x,y,z,u[0],u[1],u[2],notabsorbed]
            i+=1
        if keepdata & (i < N_max):
            track_history = track_history[:i+1,:]
        if keepdata:
            return PMT_hit_condition, t, track_history
        else:
            return PMT_hit_condition, t, total_dist, endpoint_dist, i

    # PMT SIMULATION
    def photoElectrons(self, photons): # Main monte carlo
        pe = 0. 
        for i in range(int(photons)):
            if np.random.random()<self.QE:
                pe+=1
        return pe
    def photontoElectrons(self, photons):
        e = self.photoElectrons(photons)
        for dynode in range(self.n_dynodes-1):
            delta_voltage = self.V[dynode+1]-self.V[dynode]
            e += np.random.poisson(e*delta_voltage/self.E_per_electron)
        return e

    #############################
    # RUN SIMULATION 
    #############################
    def particle_task(self, mult):
        return self.particle_path(t=self.t_initial+self.seperation_time*mult, phi_range_deg=self.particle_init_angle_range, T1_z=self.T1z, T1_width=self.T1_width, 
                                                T4_z=self.T4z, T4_width=self.T4_width, T1_radius=self.T1_radius, T4_radius=self.T4_radius, T1_corner=[self.T4_radius,-self.T4_radius],
                                                T4_corner=[self.T1_radius,-self.T1_radius], mean_free_path=self.mean_free_path_scints, 
                                                photons_per_E=self.photons_produced_per_MeV, prob_scint=self.pr_of_scintillation)
    def scint_taskT1(self, xpoint, ypoint, zpoint, time_i):
        point_i = np.hstack((xpoint,ypoint,zpoint))
        return self.scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=self.T1_radius, 
                                                        scint_plane=np.array([self.T1z,self.T1z+self.T1_width]), scint_width=self.T1_width, 
                                                        light_guide_planes=[self.T1_radius,-self.T1_radius], 
                                                        pmt_center=[self.T1_radius-4*0.5,-self.T1_radius+4*0.5,self.T1z], pmt_radius=self.PMT1_radius,
                                                        N_max=self.max_simulated_reflections, t=time_i, keepdata=False)
    def scint_taskT4(self, xpoint, ypoint, zpoint, time_i):
        point_i = np.hstack((xpoint,ypoint,zpoint))
        return self.scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=self.T4_radius, 
                                                        scint_plane=np.array([self.T4z,self.T4z+self.T4_width]), scint_width=self.T4_width, 
                                                        light_guide_planes=[-self.T4_radius,+self.T4_radius], 
                                                        pmt_center=[-self.T4_radius+4*0.5,self.T4_radius-4*0.5,self.T4z], pmt_radius=self.PMT4_radius,
                                                        N_max=self.max_simulated_reflections, t=time_i, keepdata=False)
    def days_hours_minutes(self, td):
        # days 
        # hours - 24*#days 
        # minutes - 60*#hours
        # seconds - 60*#minutes
        return td.days, td.seconds//3600-(24*td.days), (td.seconds//60)%60-(60*(td.seconds//3600)), td.total_seconds()-(60*((td.seconds//60)%60))

        """Run simulation with default 1 particle or arg[0] as number of particles and a time seperation of 'delta_t'=1e-5"""
    def run(self, *arg, **kwargs):
        freeze_support()
        if arg:
            self.num_particles = int(arg[0])
            print(f"Generating {self.num_particles} particles now...")
        else:
            self.num_particles = 1
            print(f"Generating {self.num_particles} particle now...")
        self.seperation_time = kwargs.get('delta_t', self.seperation_time) # in ps
        logstarttime = perf_counter()
        # FIND PARTICLE PATH
        times = np.zeros(1); points = np.zeros((1,3)); photons = np.zeros(1)
        with Pool(processes=cpu_count()-1) as pool:
            res = pool.map(self.particle_task, range(self.num_particles))
            for (time_i, point_i, photon_i) in res:
                times = np.concatenate((times, time_i), axis=0)
                points = np.concatenate((points,point_i), axis=0)
                photons = np.concatenate((photons,photon_i), axis=0)
        logendparticle = perf_counter()
        N = np.sum(photons)
        print("Photons generated", N)

        # SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
        # Gather TOF data
        T1_input_times = []
        T4_input_times = []
        # Gather Extra Data for analysis
        self.T1_prop_dist = []
        self.T4_prop_dist = []
        self.T1_endpoint_dist = []
        self.T4_endpoint_dist = []
        self.T1_prop_times = []
        self.T4_prop_times = []
        self.T1_interactions = []
        self.T4_interactions = []
        T1points = (points[1:])[points[1:,2] >= self.T1z]
        T1times = (times[1:])[points[1:,2] >= self.T1z]
        T1photons = (photons[1:])[points[1:,2] >= self.T1z]
        T4points = (points[1:])[points[1:,2] < self.T1z]
        T4times = (times[1:])[points[1:,2] < self.T1z]
        T4photons = (photons[1:])[points[1:,2] < self.T1z]
        print(f"Photons in T1: {np.sum(T1photons)} and Photons in T4: {np.sum(T4photons)}")
        del times; del points; del photons; # remove copies
        logstartphoton = perf_counter()
        with Pool(processes=cpu_count()) as pool:
            print("T1 Photon Propagation working...", end='\r')
            T1res = pool.starmap(self.scint_taskT1, np.repeat(np.c_[T1points,T1times],T1photons.astype(int), axis=0))
            print("T1 Photon Propagation done.     ")
            print("T4 Photon Propagation working...", end='\r')
            T4res = pool.starmap(self.scint_taskT4, np.repeat(np.c_[T4points,T4times],T4photons.astype(int), axis=0))
            print("T4 Photon Propagation done.     ")
            print("Unzipping reuslts into arrays...")
            for (T1hit_PMT, T1travel_time, T1tot_dist, T1endpt, T1bounces) in T1res:
                if T1hit_PMT:
                    T1_input_times.append(T1travel_time)
                    self.T1_prop_dist.append(T1tot_dist)
                    self.T1_endpoint_dist.append(T1endpt)
                    self.T1_prop_times.append(T1tot_dist/self.c)
                    self.T1_interactions.append(T1bounces)
            for (T4hit_PMT, T4travel_time, T4tot_dist, T4endpt, T4bounces) in T4res:
                if T4hit_PMT:
                    T4_input_times.append(T4travel_time)
                    self.T4_prop_dist.append(T4tot_dist)
                    self.T4_endpoint_dist.append(T4endpt)
                    self.T4_prop_times.append(T4tot_dist/self.c)
                    self.T4_interactions.append(T4bounces)
        logendtime = perf_counter()
        # PRINT RESULTS
        print("TIME ANALYSIS:")
        pgtime = self.days_hours_minutes(timedelta(seconds=logendparticle-logstarttime))
        phtime = self.days_hours_minutes(timedelta(seconds=logendtime-logstartphoton))
        ttime = self.days_hours_minutes(timedelta(seconds=logendtime-logstarttime))
        print(f"Generation of Particles     {pgtime[0]}days {pgtime[1]}hrs {pgtime[2]}mins {pgtime[3]:.5f}s")
        print(f"Simulation of Photon Travel {phtime[0]}days {phtime[1]}hrs {phtime[2]}mins {phtime[3]:.5f}s")
        print(f"Total Time Elapsed:         {ttime[0]}days {ttime[1]}hrs {ttime[2]}mins {ttime[3]:.5f}s")
        print("RESULTS SUMMARY:")
        print("HITS on T1",len(T1_input_times))
        print("RATIO T1   total photons", np.sum(T1photons), "total incident photons", len(T1_input_times), f"ratio={np.sum(T1photons)/len(T1_input_times):.2f}")
        print("HITS on T4",len(T4_input_times))
        print("RATIO T4   total photons ", np.sum(T4photons),"total incident photons", len(T4_input_times), f"ratio={np.sum(T4photons)/len(T4_input_times):.2f}")
        print("DISTANCE: ")
        del T1points; del T1times; del T1photons; del T4points; del T4times; del T4photons; # remove unused variables
        # print(T4_input_times)
        # BEGIN SIMULATING PMT PULSE
        signals_channelT1 = []
        signals_channelT4 = []
        output_times_channelT1 = []
        output_times_channelT4 = []
        signals = []
        output_times = []
        for t in T1_input_times:
            pmtSignal_i = self.photontoElectrons(1)
            output_times.append(self.pmt_electron_travel_time+t)
            output_times_channelT1.append(self.pmt_electron_travel_time+t)
            signals.append(pmtSignal_i)
            signals_channelT1.append(pmtSignal_i)
        for t in T4_input_times:
            pmtSignal_i = self.photontoElectrons(1)
            output_times.append(self.pmt_electron_travel_time+t)
            output_times_channelT4.append(self.pmt_electron_travel_time+t)
            signals.append(pmtSignal_i)
            signals_channelT4.append(pmtSignal_i)

        # CONVERTION Electron count to Current and save in array
        self.signals = np.array(signals) * self.q / 1e-12 * self.artificial_gain # divided by 1ps 
        self.output_times = np.array(output_times)
        self.signals_channelT1 = np.array(signals_channelT1) * self.q / 1e-12 * self.artificial_gain
        self.signals_channelT4 = np.array(signals_channelT4) * self.q / 1e-12 * self.artificial_gain * 0.6 # factor to limit pulses to 50miliamps and stop contant comparator firing. however, current should be smaller from Quantum Efficiency and current should be larger from 3kV potential difference across PMT dynodes instead of current 1kV potential difference
        self.output_times_channelT1 = np.array(output_times_channelT1)
        self.output_times_channelT4 = np.array(output_times_channelT4)
    
    # Output function
    def to_csv(self, **kwargs):
        output_extra = kwargs.get('extra_data_only', False)
        output_both = kwargs.get('output_both', False)
        # OUTPUT FORMATTING
        if output_extra or output_both:
            print("Exporting Extra Data...")
            dft1 = pd.DataFrame({'T1_prop_dist':self.T1_prop_dist,'T1_endpoint_dist':self.T1_endpoint_dist, 'T1_prop_times':self.T1_prop_times, 'T1_interactions':self.T1_interactions})
            dft4 = pd.DataFrame({'T4_prop_dist':self.T4_prop_dist,'T4_endpoint_dist':self.T4_endpoint_dist, 'T4_prop_times':self.T4_prop_times, 'T4_interactions':self.T4_interactions})
            dft1.to_csv('monte_carlo_extradata'+str(self.num_particles)+'chT1_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
            dft4.to_csv('monte_carlo_extradata'+str(self.num_particles)+'chT4_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
            if not output_both:
                return
        print("Exporing to 2 channels...")
        # for each channel
        for time,signal,ch in zip([self.output_times_channelT1,self.output_times_channelT4],[self.signals_channelT1,self.signals_channelT4],[1,4]):
            # sort by time
            df = pd.DataFrame({'time':time,'current':signal}).sort_values(by=['time'])          # return sorted dataframe
            fill_data = []                                                                      # declare empty array
            # begin padding data at time 1/5th bin width before first time stamp
            fill_data.append([df['time'].iloc[0]-self.output_bin_width/5,0])                    # add zero at beginning
            for i in range(len(time)-1):                                                        # for each time index
                if abs(df['time'].iloc[i]-df['time'].iloc[i+1]) > self.output_bin_width:        # if dt between signals is greater than minimum bin width
                    fill_data.append([df['time'].iloc[i]+self.output_bin_width/5,0])            # place zero after current signal
                    fill_data.append([df['time'].iloc[i+1]-self.output_bin_width/5,0])          # place zero before next signal
            fill_data.append([df['time'].iloc[-1]+self.output_bin_width/5,0])                   # add zero at end
            fill_data = np.array(fill_data)
            fill = pd.DataFrame(fill_data, columns=['time','current'])
            df = pd.concat([fill, df], ignore_index=True).sort_values(by=['time']).reset_index(drop=True)
            df['time'] = df['time']/1e12
            df = df[['time', 'current']] # need this for LTSpice PWL current input file to work
            df.to_csv('monte_carlo_input'+str(self.num_particles)+'ch'+str(ch)+'_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt', float_format='%.13f', header=False, index=False, sep=' ') # PWL file formatting
        print("Done!")

    def fix_artificial_gain(self, T1factor, T4factor, filedate=None, filenum=None):
        import os
        date = datetime.now().strftime('%m_%d_%Y') # Defaults
        num_part = self.num_particles
        if filedate is not None: # Take input given correct format
            date = filedate
        if filenum is not None: # Take input given correct integer number
            num_part = int(filenum)
        filename_ch1 = os.path.abspath('monte_carlo_input'+str(num_part)+'ch1_'+str(date)+'.txt')
        filename_ch4 = os.path.abspath('monte_carlo_input'+str(num_part)+'ch4_'+str(date)+'.txt')
        print("Loading files...")
        dfch1 = pd.read_csv(filename_ch1, names=['time', 'current'], sep=' ')
        dfch4 = pd.read_csv(filename_ch4, names=['time', 'current'], sep=' ')
        print("Fixing files...")
        dfch1['current'] = dfch1['current']*T1factor
        dfch4['current'] = dfch4['current']*T4factor
        dfch1.to_csv(filename_ch1, float_format='%.13f', header=False, index=False, sep=' ')
        dfch4.to_csv(filename_ch4, float_format='%.13f', header=False, index=False, sep=' ')
        print("Fixed! with T1 factor,", T1factor, " and T4 factor,", T4factor)

    def load_extradata(self, filename=None, filenum=None):
        import os
        # Default
        date = datetime.now().strftime('%m_%d_%Y')
        num = self.num_particles
        if filenum is not None:
            num = int(filenum)
        fileT1 = 'monte_carlo_extradata'+str(num)+'chT1_'+str(date)+'.txt'
        fileT4 = 'monte_carlo_extradata'+str(num)+'chT4_'+str(date)+'.txt'
        if filename is not None: # if special name use it
            if filename[-16] == '1':
                fileT1 = filename
                fileT4 = filename[:-16]+'4'+filename[-15:]
            else:
                fileT1 = filename[:-16]+'1'+filename[-15:]
                fileT4 = filename
        # Check for path errors
        print(fileT1)
        print(fileT4)
        valid = True
        if os.path.exists(fileT1) is False:
            print("Path to T1 result file below doesn't exist!")
            print(fileT1)
            valid = False
        if os.path.exists(fileT4) is False:
            print("Path to T4 result file below doesn't exist!")
            print(fileT4)
            valid = False
        if not valid:
            print("Please try again with correct path to result file")
            return
        # Load data
        T1extra_data = pd.read_csv(fileT1)
        T4extra_data = pd.read_csv(fileT4)
        # If require to replace, replace
        self.T1_prop_dist = T1extra_data['T1_prop_dist']
        self.T4_prop_dist = T4extra_data['T4_prop_dist']
        self.T1_endpoint_dist = T1extra_data['T1_endpoint_dist']
        self.T4_endpoint_dist = T4extra_data['T4_endpoint_dist']
        self.T1_prop_times = T1extra_data['T1_prop_times']
        self.T4_prop_times = T4extra_data['T4_prop_times']
        self.T1_interactions = T1extra_data['T1_interactions']
        self.T4_interactions = T4extra_data['T4_interactions']
        print("Saved to class fields")
    """
      Alternate ToF Method Assuming Seperation Width and Known # of Particles
      Uses rising edge time to compare for Time-of-Flight calculations like CMOS chip
    """
    def time_at_thresh(self, rawtime, rawVoltage, num, thresh, ch):
        out = []
        grad = np.gradient(rawVoltage) # find gradients on all data
        limit = (grad > 0) & (grad >= 0.1) # positive slope bigger than 0.1
        dtime = rawtime[limit]
        dtvoltage = rawVoltage[limit]
        tdiff = np.diff(dtime)
        condition = tdiff > self.seperation_time/1e12/10 # check if next point is 10ns away
        count = 0
        start_index = 0
        for i in range(len(tdiff)): # for number of particles we expect
            if count > num-1:
                return np.array(out)
            if condition[i]: # if condition is true at index i
                times = dtime[start_index:i+1] # take snippet of time from starting index flag to index i
                Voltages = dtvoltage[start_index:i+1]
                start_index = i+1 # reset flag to next position
                if len(times) < 1 or len(Voltages) < 1: # if no particle here then skip
                    continue
                m, b = np.polyfit(times,Voltages, deg=1) # find linear fit
                # if deg=1 returns two params slope m and y-intercept b
                # now use slope and intercept to solve for x value given our mid value y
                # y = mx + b  --> x = (y - b) / m
                out.append( (thresh - b) / m )
                count +=1 # count particle!
            
        print("Ch",ch,"counted",count,"particles!")
        if count < num: 
            print(f"Note: Counted less particles than the expected {num}")
            print("Check LTSpice that all particles were simulated.")
        return np.array(out)

    def ToF_finalize(self, tofch1, tofch4, time_limit=10e-9):
        print("ch1 length", len(tofch1))
        print("ch4 length", len(tofch4))
        out = []
        i = 0; j = 0
        lim = time_limit if time_limit is not None else 10e-9
        while (i < len(tofch1)) and (j < len(tofch4)):
            check = tofch4[j] - tofch1[i]
            if abs(check) < lim:
                out.append(abs(check))
                i+=1
                j+=1
            elif (i+1 < len(tofch1)) and (j+1 < len(tofch4)) and check < 0:
                j+=1
            else: #check < 0
                i+=1
        print("finished calculating,", len(out), "particles")
        self.FinalToF = np.array(out)
        # return np.array(out)
    
    """LTSpice Command to Analyze, Simulate and Calculate TOF"""
    def ltspice(self, filedate=None, filenum=None):
        print("\n##################################")
        print("Running LTSpice on each channel...")
        print("###################################\n")
        import os
        from PyLTSpice import SimCommander, RawRead # Use version 3.1 by pip3 install PyLTSpice==3.1
        # Make the .net file (netlist) by opening file first then saving a seperate text file
        LTC = SimCommander("PHAReduced_sim.net", parallel_sims=cpu_count()-1)
        # When running this file, two LTSpice libaries must be in same folder location:
        # LTC1.lib and LTC7.lib
        ch1ToF = 0 # declare for scope
        ch4ToF = 0 # declare for scope
        # Save the filenames of the inputs to LTSpice (Need to be Same Day and # of particles)
        # 'monte_carlo_input<X>ch1_<MM>_<DD>_<YYYY>.txt' is the format where X is # of particles if manual input is desired
        # filename_ch1 = os.path.abspath('monte_carlo_input<X>ch1_<MM>_<DD>_<YYYY>.txt')
        date = datetime.now().strftime('%m_%d_%Y') # Defaults
        num_part = self.num_particles
        if filedate is not None: # Take input given correct format
            date = filedate
        if filenum is not None: # Take input given correct integer number
            num_part = int(filenum)
        filename_ch1 = os.path.abspath('monte_carlo_input'+str(num_part)+'ch1_'+str(date)+'.txt')
        filename_ch4 = os.path.abspath('monte_carlo_input'+str(num_part)+'ch4_'+str(date)+'.txt')
        for filename,strname in zip([filename_ch1,filename_ch4],['ch1','ch4']):
            print('PWL file='+str(filename))
            LTC.set_element_model('I1', 'PWL file='+str(filename))
            LTC.add_instructions("; Simulation settings", f".tran 0 {num_part}.5u 0 0.002u") # fix this to adjust for time seperation
            # print(LTC.get_component_info('I1')) # to check if correctly set
            LTC.run(run_filename=f'PHAReduced_{strname}.net')
            LTC.wait_completion()
            print('Successful/Total Simulations: ' + str(LTC.okSim) + '/' + str(LTC.runno))
            # Now read the output
            LTR = RawRead(f"PHAReduced_{strname}.raw")
            # print(LTR.get_trace_names()) # check the outputs
            # print(LTR.get_raw_property()) # what properies does the simulation have
            # get trace gets the output waveform, get wave retrieves data from waveform object
            t = LTR.get_trace('time').get_wave() 
            compOut = LTR.get_trace('V(compout)').get_wave()
            # input ndarrays into DataFrame and fix weird negative time values
            df = pd.DataFrame({'t':np.abs(t),'V':compOut}).sort_values(by='t')
            df.to_csv('output'+str(num_part)+strname+'_'+str(date)+'.txt', header=False, index=False)
            # implement csv creation!
            # Clean up extra files
            os.remove(f"PHAReduced_{strname}.log")
            os.remove(f"PHAReduced_{strname}.op.raw")
            os.remove(f"PHAReduced_{strname}.raw")
            os.remove(f"PHAReduced_{strname}.net")
            # Remove LTSpice Object
            del LTR
       
    """ToF load LTSpice output function and call time_at_thresh and ToF_finalize""" 
    def calc_ToF(self, filedate=None, filenum=None):
        import os
        # Make final calulcation all time of flight data
        date = datetime.now().strftime('%m_%d_%Y') # Defaults
        num_part = self.num_particles
        if filedate is not None: # Take input given correct format
            date = filedate
        if filenum is not None: # Take input given correct integer number
            num_part = int(filenum)
        filename_ch1 = os.path.abspath('output'+str(num_part)+'ch1_'+str(date)+'.txt')
        filename_ch4 = os.path.abspath('output'+str(num_part)+'ch4_'+str(date)+'.txt')
        ch1 = pd.read_csv(filename_ch1, names=['t', 'V'], sep=',')
        ch4 = pd.read_csv(filename_ch4, names=['t', 'V'], sep=',')
        ch1ToF = self.time_at_thresh(ch1['t'],ch1['V'], num_part, self.CMOS_thresh, 1)
        ch4ToF = self.time_at_thresh(ch4['t'],ch4['V'], num_part, self.CMOS_thresh, 4)
        self.ToF_finalize(ch1ToF,ch4ToF) # Calculated correct time of flight
        print(pd.DataFrame(self.FinalToF).describe())

    """ToF save result data to a csv file"""
    def save_ToF(self, filename=None):
        # Default
        date = datetime.now().strftime('%m_%d_%Y')
        num_total = self.num_particles
        counted = len(self.FinalToF)
        file = 'result_'+str(counted)+'_of_'+str(num_total)+'_'+str(date)+'.txt'
        if filename is not None: # if special name use it
            file = filename
        # Output using DataFrame format and column title
        pd.DataFrame({'Time-of-Flight [s]':self.FinalToF}).to_csv(file, index=False)
    
    """
    Loads ToF onto FinalToF ndarray and will append if there already exists an array
    Unless replace=True 
    """
    def load_ToF(self, count, filename, replace=True):
        import os
        # Default
        date = datetime.now().strftime('%m_%d_%Y')
        num_total = self.num_particles
        counted = count
        file = 'result_'+str(counted)+'_of_'+str(num_total)+'_'+str(date)+'.txt'
        if filename is not None: # if special name use it
            file = filename
        # Check for path errors
        if os.path.exists(file) is False:
            print("Path to result file below doesn't exist!")
            print(file)
            print("Please try again with correct path to result file")
        # Load data
        new_ToF_data = pd.read_csv(file)
        print(new_ToF_data)
        # If require to replace, replace
        if replace is not False:
            self.FinalToF = new_ToF_data.to_numpy()
        else:
            self.FinalToF = np.append(self.FinalToF, new_ToF_data.to_numpy())
    
    def plotToF(self):
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit
        # Setup least squares loss for guassian function with 4 parameters (amplitude, mean, stddev, bias)
        fitfunc  = lambda x, p0, p1, p2, p3: p0*np.exp(-0.5*((x-p1)/p2)**2)+p3
        errfunc  = lambda p, x, y: (y - fitfunc(p, x))
        init  = [200e0, 2e-9, 1e-9, 0e0]
        plt.title(f'TOF, Total Points {len(self.FinalToF)}')
        up_bound = np.quantile(self.FinalToF, 0.95)
        low_bound = np.quantile(self.FinalToF, 0.05)
        # Plot histogram
        h, xedges, _ = plt.hist(self.FinalToF, bins=np.linspace(low_bound-1e-9,up_bound+1e-9,125), histtype='step', edgecolor='k', linewidth=2)
        xdata = (xedges[:-1] + xedges[1:]) / 2
        # Run least squares fit
        p, cov = curve_fit(fitfunc, xdata, h, p0=init)
        dp = np.sqrt(np.diag(cov))
        # find FWHM
        fwhm = 2*np.sqrt(2*np.log(2))*p[2]
        # Plot fitted guassian
        plotlabel = '$\mu=$ '+f'{p[1]/1e-9:.2f}'+'$\pm$'+f'{dp[1]/1e-9:.2f}ns\n'+'$\sigma=$ '+f'{p[2]/1e-9:.2f}'+'$\pm$'+f'{dp[2]/1e-9:.2f}ns\n'+f'FWHM = {fwhm/1e-9:.2f}ns'
        plt.plot(xedges,fitfunc(xedges, p[0], p[1], p[2], p[3]), label=plotlabel, ls='--')
        plt.grid()
        plt.legend()
        plt.show()

    #############################
    # DEBUG AND PLOTTING
    #############################

    def data_for_cylinder_along_z(self, center_x,center_y,radius,top_z,pos_z, start_theta, end_theta):
        z = np.linspace(pos_z, pos_z+top_z, 20)
        theta = np.linspace(start_theta, end_theta, 20)
        theta_grid, z_grid=np.meshgrid(theta, z)
        x_grid = radius*np.cos(theta_grid) + center_x
        y_grid = radius*np.sin(theta_grid) + center_y
        return x_grid,y_grid,z_grid
    
    def light_guide_corner(self, border_point, height_z, pos_z, border_radius, scint): 
        z = np.linspace(pos_z, height_z, 20)
        theta = np.linspace(0, -np.pi/2, 20)
        if scint == 4:
            theta = np.linspace(np.pi, np.pi/2, 20)
        theta_grid, Z = np.meshgrid(theta, z)
        X = border_radius*np.cos(theta_grid) + border_point[0]-border_radius
        Y = border_radius*np.sin(theta_grid) + border_point[1]+border_radius
        if scint == 4:
            X = border_radius*np.cos(theta_grid) + border_point[0]+border_radius
            Y = border_radius*np.sin(theta_grid) + border_point[1]-border_radius
        return X,Y,Z
    
    # Input: 3d photon position, input_time | Optional: int index of genreated scintillation points, int max iteration for PMT hit search
    def plot_scint(self, scint, photon_pos, input_time, *arg):
        # Handle arguments
        Tz = 0; T_width = 0; T_radius = 0; PMT_radius = 0; sign = 1 # intialize main variables
        if scint == 1:
            Tz = self.T1z; T_width = self.T1_width; PMTx = self.xPMT1; PMTy = self.yPMT1 
            T_radius = self.T1_radius; PMT_radius = self.PMT1_radius
        elif scint == 4:
            Tz = self.T4z; T_width = self.T4_width; PMTx = self.xPMT4; PMTy = self.yPMT4 
            T_radius = self.T4_radius; PMT_radius = self.PMT4_radius; sign = -1
        else:
            print("Incorrect scint number. Use 1 or 4")
            return
        rep = 1
        if arg and arg[0] == True:
            print("Generate random particle!")
            (time_i, point_i, _) = self.particle_task(1)
            time_i = time_i[(point_i[:,2] >= Tz) & (point_i[:,2] < (Tz+T_width))]
            point_i = point_i[(point_i[:,2] >= Tz) & (point_i[:,2] < (Tz+T_width))]
            print(f"Array Length: {len(point_i):d}")
            i = 0
            if (len(arg) > 1) and (arg[1] < len(point_i)):
                print(f"Chosen index: {arg[1]:d}")
                i = int(arg[1])
            if (len(arg) > 2) and (arg[2] < 20000):
                print(f"Max iter for PMT hit search: {arg[2]:d}")
                rep = int(arg[2])
            photon_pos = point_i[i]; input_time = time_i[i]
            print(f"pos=({photon_pos[0]:.2f},{photon_pos[1]:.2f},{photon_pos[2]:.2f}) at t={input_time:.2f}ps")
        j = 0
        hit_PMT = 0
        tracks = []
        # Search for PMT_hit (or up to limit rep)
        while j < rep:
            hit_PMT, _, tracks = self.scintillator_monte_carlo(photon_pos, notabsorbed=True, scint_radius=T_radius, 
                                                            scint_plane=np.array([Tz,Tz+T_width]), scint_width=T_width, 
                                                            light_guide_planes=[sign*T_radius,sign*-T_radius], 
                                                            pmt_center=[sign*(T_radius-4*0.5),sign*(-T_radius+4*0.5),Tz], pmt_radius=PMT_radius,
                                                            N_max=self.max_simulated_reflections, t=input_time, keepdata=True)
            j += 1
            if hit_PMT: # hit_PMT original condition
                break
        print("Returned iteration",j)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111,projection='3d')
        self.plot_full_apparatus(ax, scint)
        for x1,y1,z1,u1,v1,w1,l in zip(tracks[:-1,0],tracks[:-1,1],tracks[:-1,2], tracks[:-1,3], tracks[:-1,4],tracks[:-1,5],np.linalg.norm(tracks[:-1,3:6],ord=2, axis=1)):
            ax.quiver(x1, y1, z1, u1, v1, w1, length=l*0.5, edgecolor='k', facecolor='black')
        ax.plot(tracks[:,0],tracks[:,1],tracks[:,2], alpha=0.5, color='C1', marker='.')
        ax.text(photon_pos[0],photon_pos[1],photon_pos[2], f'hitPMT1? {hit_PMT == 1}')
        plt.show()

    def plot_full_apparatus(self, ax, *arg):
        # Create plot data
        plotT1 = True
        plotT4 = True
        lower_z_limit = self.T4z
        upper_z_limit = self.T1z+self.T1_width
        xybounds = self.T4_radius
        if len(arg) > 0:
            if arg[0] == 1: # if first scint is selected
                lower_z_limit = self.T1z
                xybounds = self.T1_radius
                plotT4 = False
            if arg[0] == 4: # if fourth scint is selected
                plotT1 = False
                upper_z_limit = self.T4z+self.T4_width
        
        if plotT1:
            XT1, YT1, ZT1 = self.data_for_cylinder_along_z(0,0,self.T1_radius,self.T1_width,self.T1z, 0, 3*np.pi/2)
            Xlg1T1, Zlg1T1 = np.meshgrid(np.linspace(0,self.T1_radius-self.T1_width*4,10),np.linspace(self.T1z,self.T1z+self.T1_width,10))
            Ylg1T1 = -np.ones((10,10))*self.T1_radius
            Ylg2T1, Zlg2T1 = np.meshgrid(np.linspace(0,-self.T1_radius+self.T1_width*4,10),np.linspace(self.T1z,self.T1z+self.T1_width,10))
            Xlg2T1 = np.ones((10,10))*self.T1_radius
            PMTxyzT1 = [+self.T1_radius-4*0.5,-self.T1_radius+4*0.5,self.T1z]
            PMT_XT1, PMT_YT1, PMT_ZT1 = self.data_for_cylinder_along_z(PMTxyzT1[0],PMTxyzT1[1],self.PMT1_radius,0.5,PMTxyzT1[2]-0.5, 0, 2*np.pi)
            SphXT1, SphYT1, SphZT1 = self.light_guide_corner([self.T1_radius,-self.T1_radius], self.T1_width+self.T1z, self.T1z, self.T1_width*4, 1)
            # Plot data 
            ax.plot_wireframe(Xlg1T1, Ylg1T1, Zlg1T1, alpha=0.2, color='C0')
            ax.plot_wireframe(Xlg2T1, Ylg2T1, Zlg2T1, alpha=0.2, color='C0')
            ax.plot_wireframe(XT1, YT1, ZT1, alpha=0.2)
            ax.plot_wireframe(PMT_XT1, PMT_YT1, PMT_ZT1, alpha=0.2, color='purple')
            ax.plot_wireframe(SphXT1, SphYT1, SphZT1, alpha=0.2, color='C0')
            ax.scatter(self.xPMT1, self.yPMT1, self.T1z, color='red', marker='o')
        if plotT4:
            XT4, YT4, ZT4 = self.data_for_cylinder_along_z(0,0,self.T4_radius,self.T4_width,self.T4z, -np.pi, np.pi/2)
            Xlg1T4, Zlg1T4 = np.meshgrid(np.linspace(0,-self.T4_radius+self.T4_width*4,10),np.linspace(self.T4z,self.T4z+self.T4_width,10))
            Ylg1T4 = +np.ones((10,10))*self.T4_radius
            Ylg2T4, Zlg2T4 = np.meshgrid(np.linspace(0,+self.T4_radius-self.T4_width*4,10),np.linspace(self.T4z,self.T4z+self.T4_width,10))
            Xlg2T4 = -np.ones((10,10))*self.T4_radius
            PMTxyzT4 = [-self.T4_radius+4*0.5,+self.T4_radius-4*0.5,self.T4z]
            PMT_XT4, PMT_YT4, PMT_ZT4 = self.data_for_cylinder_along_z(PMTxyzT4[0],PMTxyzT4[1],self.PMT4_radius,0.5,PMTxyzT4[2]-0.5, 0, 2*np.pi)
            SphXT4, SphYT4, SphZT4 = self.light_guide_corner([-self.T4_radius,self.T4_radius], self.T4_width+self.T4z, self.T4z, self.T4_width*4, 4)
            # Plot data 
            ax.plot_wireframe(Xlg1T4, Ylg1T4, Zlg1T4, alpha=0.2, color='C0')
            ax.plot_wireframe(Xlg2T4, Ylg2T4, Zlg2T4, alpha=0.2, color='C0')
            ax.plot_wireframe(XT4, YT4, ZT4, alpha=0.2)
            ax.plot_wireframe(PMT_XT4, PMT_YT4, PMT_ZT4, alpha=0.2, color='purple')
            ax.plot_wireframe(SphXT4, SphYT4, SphZT4, alpha=0.2, color='C0')
            ax.scatter(self.xPMT4, self.yPMT4, self.T4z, color='red', marker='o')
        # Fix axes
        ax.grid(True)
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.set_zlabel('z [cm]')
        ax.set_zlim([lower_z_limit-0.5,upper_z_limit+0.5])
        ax.set_ylim([-xybounds,xybounds])
        ax.set_xlim([-xybounds,xybounds])
        # ax.view_init(elev=90, azim=-90, roll=0)

    def plot_particle_dist(self, *arg):
        rep = 1
        if arg:
            if (len(arg) > 0) and (arg[0] < 200):
                print(f"Particles to simulate scintillation: {arg[0]:d}")
                rep = int(arg[0])
        print(f"Generating {rep} random particles!")
        points = np.zeros((1,3))
        lines = np.zeros((rep+1,2,3))
        for i in range(rep):
            _, points_i, _ = self.particle_task(0)
            points = np.concatenate((points,points_i), axis=0)
            lines[i,0,:] = points_i[0]
            lines[i,1,:] = points_i[-1]
        print(f"Particles shape: {points.shape}")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111,projection='3d')
        self.plot_full_apparatus(ax)
        ax.scatter(points[::100,0],points[::100,1],points[::100,2], color='green', marker='o', s=0.5)
        for i in range(rep):
            ax.plot(lines[i,:,0],lines[i,:,1],lines[i,:,2], color='C3', alpha=0.5)
        ax.set_title(f'Sample Particle Distribution: {rep} particles')
        plt.show()
        
    def plot_gen_PE(self, rep):
        times = np.zeros(1); points = np.zeros((1,3)); photons = np.zeros(1); particle_ids = np.zeros(1) 
        with Pool(processes=cpu_count()-1) as pool:
            res = pool.map(self.particle_task, range(rep))
            for i,(time_i, point_i, photon_i) in enumerate(res):
                print(time_i.shape, point_i.shape, photon_i.shape, np.repeat(i, len(time_i)).shape)
                particle_ids = np.concatenate((particle_ids,np.repeat(i, len(time_i))),axis=0)
                times = np.concatenate((times, time_i), axis=0)
                points = np.concatenate((points,point_i), axis=0)
                photons = np.concatenate((photons,photon_i), axis=0)
        N = np.sum(photons)
        print("Photons generated", N)
        # SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
        T1_input_times = []
        T4_input_times = []
        T1_prop_dist = []
        T4_prop_dist = []
        T1_output_id = []
        T4_output_id = []
        pmt_hitsT1 = 0
        pmt_hitsT4 = 0
        T1points = (points[1:])[points[1:,2] >= self.T1z]
        T1times = (times[1:])[points[1:,2] >= self.T1z]
        T1photons = (photons[1:])[points[1:,2] >= self.T1z]
        T1part_ids = (particle_ids[1:])[points[1:,2] >= self.T1z]
        T4points = (points[1:])[points[1:,2] < self.T1z]
        T4times = (times[1:])[points[1:,2] < self.T1z]
        T4photons = (photons[1:])[points[1:,2] < self.T1z]
        T4part_ids = (particle_ids[1:])[points[1:,2] < self.T1z]
        T1part_ids = np.repeat(T1part_ids, T1photons.astype(int),axis=0)
        T4part_ids = np.repeat(T4part_ids, T4photons.astype(int),axis=0)
        del times; del points; del photons;
        print(f"Photons in T1: {np.sum(T1photons)} and Photons in T4: {np.sum(T4photons)}")
        with Pool(processes=cpu_count()) as pool:
            T1res = pool.starmap(self.scint_taskT1, tqdm(np.repeat(np.c_[T1points,T1times],T1photons.astype(int), axis=0),total=np.sum(T1photons)))
            T4res = pool.starmap(self.scint_taskT4, tqdm(np.repeat(np.c_[T4points,T4times],T4photons.astype(int), axis=0),total=np.sum(T4photons)))
            for (T1hit_PMT, T1travel_time, T1tot_dist, _), T1part_id in zip(T1res, T1part_ids):
                if T1hit_PMT:
                    T1_input_times.append(T1travel_time)
                    T1_prop_dist.append(T1tot_dist)
                    T1_output_id.append(T1part_id)
                    pmt_hitsT1 +=1
            for (T4hit_PMT, T4travel_time, T4tot_dist, _), T4part_id in zip(T4res, T4part_ids):
                if T4hit_PMT:
                    T4_input_times.append(T4travel_time)
                    T4_prop_dist.append(T4tot_dist)
                    T4_output_id.append(T4part_id)
                    pmt_hitsT4 +=1
        
        _, T1cnts = np.unique(T1part_ids, return_counts=True)
        _, T4cnts = np.unique(T4part_ids, return_counts=True)
        _, T1output_cnts = np.unique(T1_output_id, return_counts=True)
        _, T4output_cnts = np.unique(T4_output_id, return_counts=True)
        import matplotlib.pyplot as plt
        fig0, ax0 = plt.subplots()
        ax0.set_title(f'Photon Generation distribution T1: {rep} particles')
        ax0.hist(T1cnts, bins=np.linspace(int(np.mean(T1cnts)-np.std(T1cnts)*2),int(np.mean(T1cnts)+np.std(T1cnts)*2),50), alpha=0.5, label='T1 Count')
        plt.show()
        fig1, ax1 = plt.subplots()
        ax1.set_title(f'Photon Generation distribution T4: {rep} particles')
        ax1.hist(T4cnts, bins=np.linspace(int(np.mean(T4cnts)-np.std(T4cnts)*2),int(np.mean(T4cnts)+np.std(T4cnts)*2),50), alpha=0.5,label='T1 Count')
        ax1.legend()
        plt.show()
        fig2, ax2 = plt.subplots()
        ax2.set_title(f'Incident Photon distribution: {min(len(T1cnts),len(T4cnts))}/{rep} particles')
        photonbins = np.linspace(int(np.mean(T1output_cnts)-np.std(T1output_cnts)*2),int(np.mean(T4output_cnts)+np.std(T4output_cnts)),50)
        ax2.hist(T1output_cnts, bins=photonbins, alpha=0.5, label='T4 Count')
        ax2.hist(T4output_cnts, bins=photonbins, alpha=0.5, label='T4 Count')
        ax2.legend()
        plt.show()
            
    def plot_xydistance_distr(self):
        if hasattr(self, 'T1_prop_dist') & hasattr(self, 'T4_prop_dist') & hasattr(self, 'T1_interactions') & hasattr(self, 'T4_interactions'):
            import matplotlib.pyplot as plt
            # import seaborn as sns
            fig0, ax0 = plt.subplots(1,2)
            fig0.set_size_inches(12,6)
            bins = np.linspace(0,40,125)
            ax0[0].hist(self.T1_prop_dist, bins=bins, histtype='step', edgecolor='k')
            ax0[0].axvline(self.T1_radius, color='C2', label='T1 radius')
            ax0[0].set_xlabel('Distance [cm]')
            ax0[0].set_ylabel('Entries')
            ax0[0].set_title('T1 Propagation XY Distance Distribution')
            ax0[0].legend()
            ax0[0].grid()
            bins = np.linspace(0,80,125)
            ax0[1].hist(self.T4_prop_dist, bins=bins, histtype='step', edgecolor='k')
            ax0[1].axvline(self.T4_radius, color='C2', label='T4 radius')
            ax0[1].set_xlabel('Distance [cm]')
            ax0[1].set_ylabel('Entries')
            ax0[1].set_title('T4 Propagation XY Distance Distribution')
            ax0[1].legend()
            ax0[1].grid()
            plt.show()
            fig1, ax1 = plt.subplots(1,2)
            fig1.set_size_inches(13,6)
            # ax1[0].scatter(self.T1_prop_dist, self.T1_interactions, s=10, alpha=0.5, facecolors='none', edgecolors='C0')
            h0 = ax1[0].hist2d(self.T1_prop_dist, self.T1_interactions, bins=[250,40], range=[[0,np.mean(self.T1_prop_dist)+np.std(self.T1_prop_dist)*5],[0,max(self.T1_interactions)]], cmin = 1)
            fig1.colorbar(h0[3], ax=ax1[0])
            # sns.kdeplot(self.T1_prop_dist, self.T1_interactions, n_levels=250, cbar=True, shade_lowest=False, cmap='viridis', ax=ax1[0])
            ax1[0].axvline(self.T1_radius, color='C2', label='T1 radius')
            ax1[0].set_xlabel('Distance [cm]')
            ax1[0].set_ylabel('# Interactions with Boundary')
            ax1[0].set_title('T1 Propagation XY Distance vs. Interactions / Reflections')
            ax1[0].legend()
            ax1[0].grid()
            # ax1[1].scatter(self.T4_prop_dist, self.T4_interactions, s=10, alpha=0.3, facecolors='none', edgecolors='C0')
            h1 = ax1[1].hist2d(self.T4_prop_dist, self.T4_interactions, bins=[250,40], range=[[0,np.mean(self.T4_prop_dist)+np.std(self.T4_prop_dist)*5],[0,max(self.T4_interactions)]], cmin = 1)
            fig1.colorbar(h1[3], ax=ax1[1])
            # sns.kdeplot(self.T4_prop_dist, self.T4_interactions, n_levels=250, cbar=True, shade_lowest=False, cmap='viridis', ax=ax1[1])
            ax1[1].axvline(self.T4_radius, color='C2',label='T4 radius')
            ax1[1].set_xlabel('Distance [cm]')
            ax1[1].set_ylabel('# Interactions with Boundary')
            ax1[1].set_title('T4 Propagation XY Distance vs. Interactions / Reflections')
            ax1[1].legend()
            ax1[1].grid()
            plt.show()
        else:
            print("Need to run simulation first with sim.run()")

    def plot_distPMT_proptime(self):
        if hasattr(self, 'T1_endpoint_dist') & hasattr(self, 'T1_prop_times') & hasattr(self, 'T4_endpoint_dist') & hasattr(self, 'T4_prop_times'):
            import matplotlib.pyplot as plt
            fig0, ax0 = plt.subplots(1,2)
            fig0.set_size_inches(12,6)
            T1proptime = np.array(self.T1_prop_times)
            T4proptime = np.array(self.T4_prop_times)
            limT1 = T1proptime < np.mean(T1proptime)+np.std(T1proptime)*3
            limT4 = T4proptime < np.mean(T4proptime)+np.std(T4proptime)*3
            # ax0[0].scatter(self.T1_prop_times, self.T1_endpoint_dist, s=1.5)
            h0 = ax0[0].hist2d(self.T1_prop_times, self.T1_endpoint_dist, bins=[500,250], range=[[0,np.mean(T1proptime)+np.std(T1proptime)*5],[5,max(self.T1_endpoint_dist)]], cmin = 1)
            fig0.colorbar(h0[3], ax=ax0[0])
            ax0[0].plot(T1proptime[limT1], T1proptime[limT1]*self.c+self.PMT1_radius, color='C3', label=f'c={self.c:.5f}cm/ps slope line')
            ax0[0].set_xlabel('time [ps]')
            ax0[0].set_ylabel('Distance [cm]')
            ax0[0].set_title(f'T1 Distance From Track to PMT vs. Propagation Time')#, outliers={np.count_nonzero(~limT1)}')
            ax0[0].legend()
            ax0[0].grid()
            # ax0[1].scatter(self.T4_prop_times, self.T4_endpoint_dist, s=1.5)
            h1 = ax0[1].hist2d(self.T4_prop_times, self.T4_endpoint_dist, bins=[1000,500], range=[[0,np.mean(T4proptime)+np.std(T4proptime)*6],[2.5,max(self.T4_endpoint_dist)]], cmin = 1)
            fig0.colorbar(h1[3], ax=ax0[1])
            ax0[1].plot(T4proptime[limT4], T4proptime[limT4]*self.c+self.PMT4_radius, color='C3', label=f'c={self.c:.5f}cm/ps slope line')
            ax0[1].set_xlabel('time [ps]')
            ax0[1].set_ylabel('Distance [cm]')
            ax0[1].set_title(f'T4 Distance From Track to PMT vs. Propagation Time')#, outliers={np.count_nonzero(~limT4)}')
            ax0[1].legend(loc='upper right')
            ax0[1].grid()
            plt.show()
        else:
            print("Need to run simulation first with sim.run()")



if __name__ == '__main__':
    sim = Simulation()
    # sim.plot_scint(1,[0,0,0],1,True,100,2000)
    # sim.plot_scint(4,[0,0,0],1,True,100,2000)
    # sim.plot_particle_dist(100)
    # sim.max_simulated_reflections = 60
    # sim.run(1)
    sim.num_particles = 4000
    # sim.to_csv(output_both=True)
    sim.load_extradata(filename='monte_carlo_extradata4000chT1_07_01_2023.txt')
    sim.plot_xydistance_distr()
    sim.plot_distPMT_proptime()
    # sim.plot_gen_PE(10)
    # sim.to_csv()
    # sim.ltspice(filedate='07_01_2023',filenum=4000)
    # sim.calc_ToF(filedate='07_01_2023',filenum=4000)
    # sim.save_ToF()
    sim.load_ToF(3858, filename='result_3858_of_4000_07_02_2023.txt')
    sim.plotToF()

