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
        self.t_rise = 800 #ps 
        self.T3z=0 #cm is the bottom of T3
        self.T1z=33.782 #cm is the bottom of T1
        self.T4z=-28.07297 #cm is the bottom of T4
        self.T1_radius = 13 #cm
        self.T4_radius = 18 #cm 
        self.xPMT4=9.5*np.cos(np.radians(110))*2.54
        self.yPMT4=9.5*np.sin(np.radians(110))*2.54
        self.xPMT1=8.*np.cos(np.radians(-45))*2.54
        self.yPMT1=8.*np.sin(np.radians(-45))*2.54
        self.PMT1_radius = 4.6/2 #cm 
        self.PMT4_radius = 4.6/2 #cm 
        # x2PMT1=8.*np.cos(np.radians(-53.72))*2.54 For test
        # y2PMT1=8.*np.sin(np.radians(-53.72))*2.54 For test
        self.n_dynodes = 8
        self.V = np.linspace(150,850,self.n_dynodes)
        # self.V = [150,300,350,600,750,850]
        self.E_per_electron = 20
        self.QE = 1 #0.23
        self.t_initial = 0 #ps
        self.particle_init_angle_range = 40 #degrees
        self.T1_width = 0.5 #cm
        self.T4_width = 1 #cm
        self.mean_free_path_scints = 0.00024 #cm
        self.photons_produced_per_MeV = 10 # electrons
        self.pr_of_scintillation = 0.8
        self.max_simulated_reflections = 16
        self.pmt_electron_travel_time = 0 # approx 16 ns
        self.artificial_gain = 1 # gain factor
        self.pr_absorption = 0.1 # probability of white paint absorbing
        self.seperation_time = 1e2 # ps 
        self.output_bin_width = 100 # ps
        self.num_particles = 1 # Muons
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
    def quadrant_sign(self, corner_pt):
        negatives = 0
        if corner_pt[0] < 0:
            negatives = (negatives + 1) % 2
        if corner_pt[1] < 0:
            negatives = (negatives + 1) % 2
        return -1*negatives

    # DISTANCE 2-DIM CIRCLE WITH LINE SEGMENT
    # t = -D . ∆ ± √(D . ∆)^2 - |D|^2(|∆|^2 - R^2)
    #     over |D|^2
    # ARGUMENTS : # 3d directional vector, 3d point, center of scintillator, radius of scintillator, use corner circle boolean
    def distance_circle(self, u, o, center, radius, quadrant=False): 
        P = o
        D = u*-1 if np.dot(u,P) < 0 else u # does a normalized vector in 3d equate to normalized vector in 2d?
        C = center
        R = radius
        bigDelta = np.array(P)-np.array(C)

        magDsq = self.mag(D)**2
        magDeltasq = self.mag(bigDelta)**2
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
    def distance_plane(self, u, o, plane, dim):                                     
        P = o
        if dim==2:
            d_plane = plane[0] if u[dim] < 0 else plane[1]                    # make sure direction matches location of plane 
        else:
            d_plane = plane
        return np.abs((d_plane - P[dim])/u[dim])


    # SOLVE FOR DISTANCE LOGIC FUNCTION
    def distance_solver(self, u, o, center, radius, plane_z, corner_center, corner_radius, pmt_center, pmt_radius):
        dcircle = self.distance_circle(u,o,center,radius)                          # checks distance to circle boundary
        dplane_z = self.distance_plane(u,o,plane_z,dim=2)                          # checks distance to z boundary in general scint
        dist = dplane_z if dcircle > dplane_z else dcircle
        temp_o = o+dist*u
        PMT_cond = False
        if (temp_o[0] > 0) & (temp_o[1] < 0) & ((temp_o[0]**2+temp_o[1]**2) >= radius**2-1):
            dplanex = self.distance_plane(u,o,radius,dim=0)                        # checks distance to x boundary
            dplaney = self.distance_plane(u,o,-radius,dim=1)                       # checks distance to y boundary
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
    def photon_interaction(self, o, u, n, notabsorbed, scint_plane):
        u_r = u - 2*np.dot(u, n)*n                      # u_new = u - 2 (u . n)*n
        v = u*-1 if np.dot(u,n) < 0 else u
        theta = np.arcsin(self.mag(np.cross(v,n))/(self.mag(u)*self.mag(n)))
        inside_sqrt = ((self.n_1/self.n_2)*np.sin(theta))**2
        sqrt_term = np.sqrt(1 - inside_sqrt)            # cos(theta)_transmission
        Rs = np.abs((self.n_1*np.cos(theta) - self.n_2*sqrt_term)/(self.n_1*np.cos(theta) + self.n_2*sqrt_term))**2
        Rp = np.abs((self.n_1*sqrt_term - self.n_2*np.cos(theta))/(self.n_1*sqrt_term + self.n_2*np.cos(theta)))**2
                                                        # Determine probability of reflectance
        if np.random.random() < ((Rs+Rp)/2):            # if random chance is high enough reflect !
            return self.normalize(u_r), True                 # return full internal reflection and not absorbed is True
                                                        # else photon is transmitted to white paint
        elif np.random.random() < self.pr_absorption:                 # does it get absorbed? change probability when you get more data
            return self.normalize(u_r), False                # not absorbed is False
        else:                                           # no it didn't get absorbed!
            theta_new = random.uniform(-np.pi/2,np.pi/2)       # new theta direction of photon
            phi_new = random.uniform(-np.pi, np.pi)           # new phi   direction of photon
            new_u = self.normalize(np.array([np.sin(phi_new)*np.cos(theta_new),np.sin(phi_new)*np.sin(theta_new),np.cos(phi_new)]))
            u_r = self.reemission_angle_factor*new_u + n
            return self.normalize(u_r), True            # new small change in direction (should be random), and not absorbed is True

    # Calculate n vector for all planes and surfaces in apparatus
    def n_vec_calculate(self, o, scint_plane, light_guide_planes, corner_center, corner_radius):
        if o[2] == scint_plane[0]:                                      # bottom of scint
            return np.array([0,0,+1])
        elif o[2] == scint_plane[1]:                                    # top of scint
            return np.array([0,0,-1])
        elif o[0] == light_guide_planes[0]:                             # y plane of light guide 
            return np.array([0,+1,0])
        elif o[1] == light_guide_planes[1]:                             # x plane of light guide
            return np.array([-1,0,0])
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
        theta = random.uniform(0,2*np.pi)                                                 # random theta in circle above T1
        phi = random.uniform(np.pi-phi_range_deg*np.pi/180/2,np.pi+phi_range_deg*np.pi/180/2) # phi angle pointing in -k given phi range
        maxdist = np.random.random()*T1_radius/2                                          # half the radius of T1
        round_const = self.round_to_sig(mean_free_path)
        o = np.float64((maxdist*np.cos(theta), maxdist*np.sin(theta), T1_z+T1_width+2))   # x, y, top of T1_z+2
        u = np.array((np.cos(theta)*np.sin(phi),np.sin(theta)*np.sin(phi),np.cos(phi)),dtype=np.float64)
        # print(f"u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
        photons = [0]                                                                     # begin photon array
        points = [o]                                                                      # current point 
        times = [t]                                                                       # current t 
        cur_o = points[-1]                                                                 # current z 
        next_o = (cur_o+mean_free_path*u).round(round_const)                                  # next z step
        inside_scint = False
        while next_o[2] >= T4_z:
            if not inside_scint:
                distT1 = np.abs((T1_z+T1_width - cur_o[2])/u[2])
                distT4 = np.abs((T4_z+T4_width - cur_o[2])/u[2])
                dist = distT4 if next_o[2] < T1_z else distT1
                next_o = (cur_o+dist*u).round(round_const)
                inside_T1 = ((np.sqrt(np.sum(next_o[0:2]**2)) < T1_radius) | ((next_o[0] < T1_corner[0]) & (next_o[1] > T1_corner[1])))
                inside_T4 = ((np.sqrt(np.sum(next_o[0:2]**2)) < T4_radius) | ((next_o[0] < T4_corner[0]) & (next_o[1] > T4_corner[1])))
                inside_scint = inside_T4 if next_o[2] < T1_z else inside_T1
                if not inside_scint:
                    continue
                t +=  dist/self.c                                                               # calculate time in ps passed
                times.append(t)
                points.append(next_o)
                phot = np.random.poisson(photons_per_E)
                if np.random.random() < prob_scint: photons.append(phot)
                else: photons.append(0)
                cur_o = points[-1]                                                                 # current z 
                next_o = (cur_o+mean_free_path*u).round(round_const)  # make this just the whole vector calculation
            for Tbottom,Ttop,Tradius,Tcorner in [(T1_z,T1_z+T1_width,T1_radius,T1_corner),(T4_z,T4_z+T4_width,T4_radius,T4_corner)]:
                inside_scint = (next_o[2] <= (Ttop)) & (next_o[2] >= Tbottom) & ((np.sqrt(np.sum(next_o[0:2]**2)) < Tradius) | ((next_o[0] < Tcorner[0]) & (next_o[1] > Tcorner[1])))
                while inside_scint:
                    t +=  mean_free_path/self.c
                    times.append(t)
                    points.append(cur_o+mean_free_path*u)
                    phot = np.random.poisson(photons_per_E)
                    if np.random.random() < prob_scint: photons.append(phot)
                    else: photons.append(0)
                    cur_o = points[-1]                                                                 # current z 
                    next_o = (cur_o+mean_free_path*u).round(round_const)
                    inside_scint = (next_o[2] <= (Ttop)) & (next_o[2] >= Tbottom) & ((np.sqrt(np.sum(next_o[0:2]**2)) < Tradius) | ((next_o[0] < Tcorner[0]) & (next_o[1] > Tcorner[1])))
        # print(f"time elapsed in scintillators: {np.abs(times[1]-times[-1]):.2f}ps total scintillation points: {len(points[1:])}")
        return np.array(times, dtype=np.float64)[1:], np.array(points, dtype=np.float64)[1:], np.array(photons[1:], dtype=np.float64)


    def scintillator_monte_carlo(self, o, notabsorbed, scint_radius, scint_plane, scint_width, light_guide_planes, pmt_center, pmt_radius, N_max, dt, keepdata):
        if keepdata: track_history = np.zeros((N_max+1,7))         # x, y history of Photon
        corner_radius = scint_width*4
        corner_center = [scint_radius-corner_radius,-scint_radius+corner_radius,scint_plane[0]]
        theta = random.uniform(0,2*np.pi)             # first theta direction of photon
        phi = random.uniform(0,np.pi)                 # first phi   direction of photon
        PMT_hit_condition = False
        u = np.array([np.sin(phi)*np.cos(theta),np.sin(phi)*np.sin(theta),np.cos(phi)]) # first direction unit vector
        if keepdata: track_history[0,:] = [o[0],o[1],o[2],u[0],u[1],u[2],notabsorbed]
        i = 1
        while (i < N_max+1) & (not PMT_hit_condition) & (notabsorbed is True):
            ds, PMT_hit_condition = self.distance_solver(u, o, np.array([0,0,scint_plane[0]]),scint_radius, scint_plane, corner_center, corner_radius, pmt_center, pmt_radius)
            x, y, z = o+ds*u
            o = np.array([x, y, np.abs(z) if np.abs(z-scint_plane).any() < 1e-5 else z])
            dt += np.abs(ds)/self.c                        # time taken in ps traveling in direction theta
    #         print(f"step {i}: ds={ds:.2f}cm dt={dt:.2f}ps Absorbed?={not notabsorbed} xyz =({x:.2f},{y:.2f},{z:.2f}) u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
            n = self.n_vec_calculate(o, scint_plane, light_guide_planes, corner_center, corner_radius)
            u, notabsorbed = self.photon_interaction(o, u, n, notabsorbed, scint_plane)
            if keepdata: track_history[i] = [x,y,z,u[0],u[1],u[2],notabsorbed]
            i+=1
        if keepdata & (i < N_max+1):
            track_history = track_history[:i,:]
        if keepdata:
            return PMT_hit_condition, dt, track_history
        else:
            return PMT_hit_condition, dt, 1 # placeholder

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
                                                        N_max=self.max_simulated_reflections, dt=time_i, keepdata=False)
    def scint_taskT4(self, xpoint, ypoint, zpoint, time_i):
        point_i = np.hstack((xpoint,ypoint,zpoint))
        return self.scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=self.T4_radius, 
                                                        scint_plane=np.array([self.T4z,self.T4z+self.T4_width]), scint_width=self.T4_width, 
                                                        light_guide_planes=[self.T4_radius,-self.T4_radius], 
                                                        pmt_center=[self.T4_radius-4*0.5,-self.T4_radius+4*0.5,self.T4z], pmt_radius=self.PMT4_radius,
                                                        N_max=self.max_simulated_reflections, dt=time_i, keepdata=False)
    def days_hours_minutes(self, td):
        # days 
        # hours - 24*#days 
        # minutes - 60*#hours
        # seconds - 60*#minutes
        return td.days, td.seconds//3600-(24*td.days), (td.seconds//60)%60-(60*(td.seconds//3600)), td.total_seconds()-(60*((td.seconds//60)%60))
    def run(self, *arg, **kwargs):
        freeze_support()
        if arg:
            self.num_particles = int(arg[0])
            print(f"Generating {self.num_particles} particles now...")
        else:
            self.num_particles = 1
            print(f"Generating {self.num_particles} particle now...")
        self.seperation_time = kwargs.get('delta_t', 1e6) # in ps
        logstarttime = perf_counter()
        # FIND PARTICLE PATH
        times = np.zeros(1); points = np.zeros((1,3)); photons = np.zeros(1)
        with Pool(processes=cpu_count()-1) as pool:
            res = pool.map(self.particle_task, range(self.num_particles))
            for (time_i, point_i, photon_i) in res:
                times = np.append(times, time_i, axis=0)
                points = np.append(points,point_i, axis=0)
                photons = np.append(photons,photon_i)
        logendparticle = perf_counter()
        N = np.sum(photons)
        print("Photons generated", N)

        # SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
        T1_input_times = []
        T4_input_times = []
        pmt_hits = 0
        T1points = (points[1:])[points[1:,2] >= self.T1z]
        T1times = (times[1:])[points[1:,2] >= self.T1z]
        T1photons = (photons[1:])[points[1:,2] >= self.T1z]
        T4points = (points[1:])[points[1:,2] < self.T1z]
        T4times = (times[1:])[points[1:,2] < self.T1z]
        T4photons = (photons[1:])[points[1:,2] < self.T1z]
        print(f"Photons in T1: {np.sum(T1photons)} and Photons in T4: {np.sum(T4photons)}")
        logstartphoton = perf_counter()
        with Pool(processes=cpu_count()) as pool:
            T1res = pool.starmap(self.scint_taskT1, tqdm(np.repeat(np.c_[T1points,T1times],T1photons.astype(int), axis=0),total=np.sum(T1photons)))
            T4res = pool.starmap(self.scint_taskT4, tqdm(np.repeat(np.c_[T4points,T4times],T4photons.astype(int), axis=0),total=np.sum(T4photons)))
            for (T1hit_PMT, T1travel_time, _) in T1res:
                if T1hit_PMT:
                    T1_input_times.append(T1travel_time)
                    pmt_hits +=1
            for (T4hit_PMT, T4travel_time, _) in T4res:
                if T4hit_PMT:
                    T4_input_times.append(T4travel_time)
                    pmt_hits +=1
        logendtime = perf_counter()
        # PRINT RESULTS
        print("TIME ANALYSIS:")
        pgtime = self.days_hours_minutes(timedelta(seconds=logendparticle-logstarttime))
        phtime = self.days_hours_minutes(timedelta(seconds=logendtime-logstartphoton))
        ttime = self.days_hours_minutes(timedelta(seconds=logendtime-logstarttime))
        print(f"Generation of Particles     {pgtime[0]}days {pgtime[1]}hrs {pgtime[2]}mins {pgtime[3]}s")
        print(f"Simulation of Photon Travel {phtime[0]}days {phtime[1]}hrs {phtime[2]}mins {phtime[3]}s")
        print(f"Total Time Elapsed:         {ttime[0]}days {ttime[1]}hrs {ttime[2]}mins {ttime[3]}s")
        print("RESULTS:")
        print("HITS on T1",len(T1_input_times))
        # print(T1_input_times)
        print("HITS on T4",len(T4_input_times))
        # print(T4_input_times)
        # BEGIN SIMULATING PMT PULSE
        signals_channelT1 = []
        signals_channelT4 = []
        output_times_channelT1 = []
        output_times_channelT4 = []
        signals = []
        output_times = []
        for i,t in enumerate(T1_input_times):
            pmtSignal_i = self.photontoElectrons(1)
            output_times.append(self.pmt_electron_travel_time+t)
            output_times_channelT1.append(self.pmt_electron_travel_time+t)
            signals.append(pmtSignal_i)
            signals_channelT1.append(pmtSignal_i)
        for i,t in enumerate(T4_input_times):
            pmtSignal_i = self.photontoElectrons(1)
            output_times.append(self.pmt_electron_travel_time+t)
            output_times_channelT4.append(self.pmt_electron_travel_time+t)
            signals.append(pmtSignal_i)
            signals_channelT4.append(pmtSignal_i)

        # CONVERTION Electron count to Current and save in array
        self.signals = np.array(signals) * self.q / 1e-12 * self.artificial_gain # divided by 1ps 
        self.output_times = np.array(output_times)
        self.signals_channelT1 = np.array(signals_channelT1) * self.q / 1e-12 * self.artificial_gain
        self.signals_channelT4 = np.array(signals_channelT4) * self.q / 1e-12 * self.artificial_gain
        self.output_times_channelT1 = np.array(output_times_channelT1)
        self.output_times_channelT4 = np.array(output_times_channelT4)
    
    # Output function
    def to_csv(self):
        
        # OUTPUT FORMATTING
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
            df.to_csv('monte_carlo_input'+str(self.num_particles)+'ch'+str(ch)+'_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt', float_format='%.13f', header=False, index=False, sep=' ')
            print(df)
        print("Done!")
    
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
        condition = tdiff > 10e-9 # check if next point is 10ns away
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
        LTC = SimCommander("PHAReduced_sim.net")
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
        ch1ToF = self.time_at_thresh(ch1['t'],ch1['V'], self.num_particles, self.CMOS_thresh, 1)
        ch4ToF = self.time_at_thresh(ch4['t'],ch4['V'], self.num_particles, self.CMOS_thresh, 4)
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
    def load_ToF(self, filename, replace=True):
        import os
        # Default
        date = datetime.now().strftime('%m_%d_%Y')
        num_total = self.num_particles
        counted = self.num_particles
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
        plt.title(f'TOF, Total Points {len(self.FinalToF)}')
        plt.hist(self.FinalToF, bins=np.linspace(0,5e-9,125), histtype='step', edgecolor='k', linewidth=2)
        plt.axvline(np.mean(self.FinalToF), label=f'$\mu$={np.mean(self.FinalToF):.2e}', color='C1')
        plt.grid()
        plt.legend()
        plt.show()
        return
    #############################
    # DEBUG AND PLOTTING
    #############################

    def data_for_cylinder_along_z(self, center_x,center_y,radius,height_z,pos_z, start_theta, end_theta):
        z = np.linspace(pos_z, pos_z+height_z, 20)
        theta = np.linspace(start_theta, end_theta, 20)
        theta_grid, z_grid=np.meshgrid(theta, z)
        x_grid = radius*np.cos(theta_grid) + center_x
        y_grid = radius*np.sin(theta_grid) + center_y
        return x_grid,y_grid,z_grid
    
    # Input: 3d photon position, input_time | Optional: int index of genreated scintillation points, int max iteration for PMT hit search
    def plot_scint(self, scint, photon_pos, input_time, *arg):
        # Handle arguments
        Tz = 0; T_width = 0; T_radius = 0; PMT_radius = 0 # intialize main variables
        if scint == 1:
            Tz = self.T1z; T_width = self.T1_width; PMTx = self.xPMT1; PMTy = self.yPMT1 
            T_radius = self.T1_radius; PMT_radius = self.PMT1_radius
        elif scint == 4:
            Tz = self.T4z; T_width = self.T4_width; PMTx = self.xPMT4; PMTy = self.yPMT4 
            T_radius = self.T4_radius; PMT_radius = self.PMT4_radius
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
        tracks = []
        # Search for PMT_hit (or up to limit rep)
        while j < rep:
            hit_PMT, _, tracks = self.scintillator_monte_carlo(photon_pos, notabsorbed=True, scint_radius=T_radius, 
                                                            scint_plane=np.array([Tz,Tz+T_width]), scint_width=T_width, 
                                                            light_guide_planes=[T_radius,-T_radius], 
                                                            pmt_center=[T_radius-4*0.5,-T_radius+4*0.5,Tz], pmt_radius=PMT_radius,
                                                            N_max=self.max_simulated_reflections, dt=input_time, keepdata=True)
            j += 1
            if np.sqrt(np.sum(photon_pos[0:2]**2)) > T_radius: # hit_PMT original condition
                break
        print("Returned iteration",j)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111,projection='3d')
        # Create plot data
        X, Y, Z = self.data_for_cylinder_along_z(0,0,T_radius,T_width,Tz, 0, 3*np.pi/2)
        Xlg1, Zlg1 = np.meshgrid(np.linspace(0,T_radius,10),np.linspace(Tz,Tz+T_width,10))
        Ylg1 = -np.ones((10,10))*T_radius
        Ylg2, Zlg2 = np.meshgrid(np.linspace(0,-T_radius,10),np.linspace(Tz,Tz+T_width,10))
        Xlg2 = np.ones((10,10))*T_radius
        PMTxyz = [T_radius-4*0.5,-T_radius+4*0.5,Tz]
        PMT_X, PMT_Y, PMT_Z = self.data_for_cylinder_along_z(PMTxyz[0],PMTxyz[1],PMT_radius,0.5,PMTxyz[2]-0.5, 0, 2*np.pi)
        # Plot data
        ax.plot_wireframe(Xlg1, Ylg1, Zlg1, alpha=0.2, color='C0')
        ax.plot_wireframe(Xlg2, Ylg2, Zlg2, alpha=0.2, color='C0')
        ax.plot_wireframe(X, Y, Z, alpha=0.2)
        ax.plot_wireframe(PMT_X, PMT_Y, PMT_Z, alpha=0.2, color='purple')
        ax.scatter(photon_pos[0],photon_pos[1],photon_pos[2], color='red', marker='o')
        ax.scatter(PMTx, PMTy, Tz, color='red', marker='o')
        ax.quiver(tracks[:,0],tracks[:,1],tracks[:,2], tracks[:,3], tracks[:,4],tracks[:,5], length=1, edgecolor='k', facecolor='black', linewidth=1.5)
        ax.plot(tracks[:,0],tracks[:,1],tracks[:,2], alpha=0.5, color='C1', marker='.')
        ax.text(photon_pos[0],photon_pos[1],photon_pos[2], f'hitPMT1? {hit_PMT}')
        # Fix axes
        ax.grid(True)
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.set_zlabel('z [cm]')
        ax.set_zlim([Tz-0.5,Tz+T_width+0.5])
        ax.set_ylim([-T_radius,T_radius])
        ax.set_xlim([-T_radius,T_radius])
        # ax.view_init(elev=90, azim=-90, roll=0)
        plt.show()




if __name__ == '__main__':
    sim = Simulation()
    # sim.plot_scint(1,[0,0,0],1,True,100,1000)
    sim.plot_scint(4,[0,0,0],1,True,100,1000)
    # sim.artificial_gain = 2 # was 3
    # sim.num_particles = 5
    # sim.run(5)
    # sim.to_csv()
    # sim.ltspice(filedate='05_07_2023',filenum=20000)
    # sim.calc_ToF(filedate='05_07_2023',filenum=20000)
    # sim.save_ToF()
    # sim.load_ToF(filename='result_3867_of_8000_05_17_2023.txt')
    # sim.plotToF()