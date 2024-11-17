#############################
# AESOP-Lite Monte Carlo (2024)
# Simulation of Particle to PMT
# Created by Liam Branch, Michelle Picardo, Robert Johnson
# UCSC Alumni + UCSC Professor
#############################

# import numpy as np
import utility as np
import h5py
from pandas import DataFrame, read_csv, concat
from tqdm import tqdm
from random import uniform
from time import perf_counter_ns
from datetime import timedelta, datetime
from multiprocessing import Pool, cpu_count, freeze_support, Manager, set_start_method


class Simulation:
    
    def __init__(self):
        #############################
        # CONSTANTS
        #############################
        self.c = 0.0299792/1.58 # Speed of Light in cm / ps
        self.q = 1.60217663e-19 # charge of electron columbs
        # CONSTRAINT n_1 <= n_2
        # SIMULATOR FIELDS
        self.n_1 = 1.000293 # Sample index of refraction of air
        self.n_2 = 1.58 # 1.58 for EJ-200
        self.T3z=0 #cm is the bottom of T3
        self.T1z=33.782 #cm is the bottom of T1
        self.T4z=-28.07297 #cm is the bottom of T4
        self.T1_radius = 13 #cm
        self.T4_radius = 18 #cm
        self.T1_width = 0.5 #cm
        self.T4_width = 1 #cm
        self.T1top = self.T1z+self.T1_width
        self.T4top = self.T4z+self.T4_width
        self.T1_corner_radius = self.T1_width*4
        self.T4_corner_radius = self.T4_width*4
        self.T1_corner_center = [self.T1_radius-self.T1_corner_radius,self.T1_corner_radius-self.T1_radius, self.T1_radius]
        self.T4_corner_center = [self.T4_corner_radius-self.T4_radius, self.T4_radius-self.T4_corner_radius, self.T4_radius]
        self.PMT1_center = [self.T1_radius-4*0.5,-self.T1_radius+4*0.5,self.T1z]
        self.PMT4_center = [-self.T4_radius+4*0.5,self.T4_radius-4*0.5,self.T4z]
        self.xPMT4=9.5*np.cos(np.radians(110))*2.54
        self.yPMT4=9.5*np.sin(np.radians(110))*2.54
        self.xPMT1=8.*np.cos(np.radians(-45))*2.54 # x2PMT1=8.*np.cos(np.radians(-53.72))*2.54 For test
        self.yPMT1=8.*np.sin(np.radians(-45))*2.54 # y2PMT1=8.*np.sin(np.radians(-53.72))*2.54 For test
        self.PMT1_radius = 4.6/2 #cm need to change this to 46 milimeters or 0.046 cm
        self.PMT4_radius = 4.6/2 #cm 
        # PMT SIGNAL GENERATION FIELDS 
        self.n_dynodes = 8
        self.V = np.linspace(150,850,self.n_dynodes)
        # self.V = [150,300,350,600,750,850]
        self.E_per_electron = 20
        self.QE = 1 #0.23
        self.sigma_smoothing = 400 #ps
        self.t_initial = 0 #ps
        self.particle_init_angle_range = 40 #degrees
        self.particle_gen_area = self.T1_radius
        self.particle_gen_z = self.T1z+self.T1_width + 2 #cm
        self.mean_free_path_scints = 24e-5 #cm or 80 micro meters
        self.photons_produced_per_MeV = 10 # true value is closer to 10000 per 1MeV
        self.pr_of_scintillation = 0.8
        self.max_simulated_reflections = 40
        self.pmt_electron_travel_time = 0 # approx 16 ns
        self.artificial_gain = 1 # gain factor
        self.max_pmt_current_output = 80e-3 # mA
        self.pr_absorption = 0.1 # probability of boundary absorbing
        self.seperation_time = 1e5 # ps 
        self.output_bin_width = 100 # ps
        self.num_particles = 1 # default Muons
        self.CMOS_thresh = 1.5 # V for rising edge detector
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
        np.divide(x, np.norm(x), out = x) # removing temp array 
        return x
    
    # LIGHT GUIDE CONDITION
    def lg_condition(self, corner_pt, scint_corner, scint_num):
        # removed temp variable ret
        if scint_num == 4:
            return (corner_pt[0] > 0) & (corner_pt[0] < scint_corner[0]) & (corner_pt[1] < 0) & (corner_pt[1] > scint_corner[1])
        return (corner_pt[0] > 0) & (corner_pt[0] < scint_corner[0]) & (corner_pt[1] < 0) & (corner_pt[1] > scint_corner[1])

    # SCINT RADIUS CONDITION
    def scint_condition(self, corner_pt, scint_radius, scint_num):
        # reduced computation 
        radius_sq = corner_pt[0]**2 + corner_pt[1]**2
        if scint_num == 4:
            return np.sqrt(radius_sq) < self.T4_radius
        return np.sqrt(radius_sq) < self.T1_radius

    # DISTANCE 2-DIM CIRCLE WITH LINE SEGMENT
    # t = -D . ∆ ± √(D . ∆)^2 - |D|^2(|∆|^2 - R^2)
    #     over |D|^2
    # ARGUMENTS : # 3d directional vector, 3d point, center of scintillator, radius of scintillator, use corner circle boolean
    def distance_circle(self, u, o, center, radius, quadrant=False): 
        # Calculate the dot prod. of u and o 
        cond = np.dot(u,o) < 0
        # Calculate the normalized direction vector D 
        D = -u if cond else u # does a normalized vector in 3d equate to normalized vector in 2d?- the process is the same and it will collapse to 2d if the 3d vector has a zero as a component 
        # Calculate the components of bigDelta, removed array creation 
        bigDelta_x = o[0] - center[0]
        bigDelta_y = o[1] - center[1]
        # Calculate the squared mag of D aka the dot of D aka squared norm
        magDsq = np.dot(D,D)
        # Calculate the squared mag. of bigDelta
        magDeltasq = bigDelta_x**2 + bigDelta_y**2
        # Calculate the dot prod. of D and bigDelta
        DdotDelta = D[0]*bigDelta_x + D[1]*bigDelta_y
        # Calculate discriminant
        discriminant = DdotDelta**2 - magDsq * (magDeltasq - radius**2)

        if discriminant < 0:
            # Release memory for intermediate variables
            D = None; bigDelta_x = None; bigDelta_y = None; magDsq = None;
            magDeltasq = None; DdotDelta = None;
            return 100 # some large value that won't be chosen because photon has no intersection with circle
        
        sqrt_term = np.sqrt(discriminant)/magDsq
        b_term = -DdotDelta/magDsq
        rootA = b_term - sqrt_term
        rootB = b_term + sqrt_term

        if quadrant is not False: # if in corner don't use the other 3/4ths of the circle to find distance only 2nd or 4th quadrant part
            # Release memory for intermediate variables
            D = None; bigDelta_x = None; bigDelta_y = None; magDsq = None;
            magDeltasq = None; DdotDelta = None;
            return np.abs(rootA) if np.abs(rootA) > np.abs(rootB) else np.abs(rootB)
        
        # Release memory for intermediate variables
        D = None; bigDelta_x = None; bigDelta_y = None; magDsq = None;
        magDeltasq = None; DdotDelta = None;
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
            # release
            dcircle = None; dplane_z = None; dist = None; temp_o = None; dplanex = None; dplaney = None; dplanez = None; 
            dcorner = None;
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
            # release
            dcircle = None; dplane_z = None; dist = None; temp_o = None; dplanex = None; dplaney = None; dplanez = None; 
            dcorner = None;
            return light_guide_dist, PMT_cond
        else:
            # release
            dcircle = None; dplane_z = None; temp_o = None; dplanex = None; dplaney = None; dplanez = None; 
            dcorner = None;
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
        theta = np.arcsin(np.norm(np.cross(v,n))/(np.norm(u)*np.norm(n)))
        inside_sqrt = ((self.n_1/self.n_2)*np.sin(theta))**2
        sqrt_term = np.sqrt(1 - inside_sqrt)                    # cos(theta)_transmission
        Rs = np.abs((self.n_1*np.cos(theta) - self.n_2*sqrt_term)/(self.n_1*np.cos(theta) + self.n_2*sqrt_term))**2
        Rp = np.abs((self.n_1*sqrt_term - self.n_2*np.cos(theta))/(self.n_1*sqrt_term + self.n_2*np.cos(theta)))**2
                                                                # Determine probability of reflectance
        if np.random() < ((Rs+Rp)/2):                    # if random chance is high enough reflect !
            # release
            v = None; theta = None; inside_sqrt = None; sqrt_term = None; Rs = None; Rp = None;             
            return self.normalize(u_r), True                        # return full internal reflection and not absorbed is True
                                                                # else photon is transmitted to white paint
        elif np.random() < self.pr_absorption:               # does it get absorbed? change probability when you get more data
            # release
            v = None; theta = None; inside_sqrt = None; sqrt_term = None; Rs = None; Rp = None;
            return self.normalize(u_r), False                       # not absorbed is False
        else:                                                   # no it didn't get absorbed!
            theta_new = uniform(-np.pi/2,np.pi/2)            # new theta direction of photon
            phi_new = uniform(-np.pi, np.pi)                 # new phi   direction of photon
            new_u = self.normalize(np.array([np.sin(phi_new)*np.cos(theta_new),np.sin(phi_new)*np.sin(theta_new),np.cos(phi_new)]))
            u_r = self.reemission_angle_factor*new_u + n
            # release
            v = None; theta = None; inside_sqrt = None; sqrt_term = None; Rs = None; Rp = None;
            theta_new = None; phi_new = None;
            return self.normalize(u_r), True                        # new small change in direction (should be random), and not absorbed is True

    # Predefined unit vectors as NumPy arrays - to prevent recreation of np.arrays
    unit_z = np.array([0, 0, 1])
    unit_neg_z = np.array([0, 0, -1])
    # Calculate n vector for all planes and surfaces in apparatus
    def n_vec_calculate(self, o, scint_plane, light_guide_planes, corner_center, corner_radius=None):
        if o[2] == scint_plane[0]:                                      # bottom of scint
            return self.unit_z
        elif o[2] == scint_plane[1]:                                    # top of scint
            return self.unit_neg_z
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
        theta = uniform(0,2*np.pi)                                                     # random theta in circle above T1
        phi = uniform(np.pi-phi_range_deg*np.pi/180/2,np.pi+phi_range_deg*np.pi/180/2) # phi angle pointing in -k given phi range
        maxdist = np.random()*self.particle_gen_area                                   # radius of generation
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
                    theta = uniform(0,2*np.pi)                                                 # reset random theta in circle above T1
                    phi = uniform(np.pi-phi_range_deg*np.pi/180/2,np.pi+phi_range_deg*np.pi/180/2) # reset phi angle pointing in -k given phi range
                    maxdist = np.random()*T1_radius/2                                          # reset random point inside half the radius of T1
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
                    phot = np.poisson(photons_per_E)
                    if np.random() < prob_scint: photons.append(phot)
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
                    phot = np.poisson(photons_per_E)
                    if np.random() < prob_scint: photons.append(phot)
                    else: photons.append(0)
                    cur_o = points[-1]                                             # current point 
                    next_o = (cur_o+mean_free_path*u).round(round_const)           # next point
                    inside_scint = (next_o[2] <= (Ttop)) & (next_o[2] >= Tbottom) & (self.scint_condition(next_o, Tradius, num) | self.lg_condition(next_o, Tcorner, num))
        # Release
        theta = None; phi = None; maxdist = None; round_const = None; o = None; u = None; cur_o = None; next_o = None; inside_scint = None; missed = None; distT1 = None; distT4 = None; 
        dist = None; check = None; inside_T1 = None; inside_T4 = None; scint_cond = None; t = None; phot = None; 
        # write to file compressed arrays float 64s
        return np.array(times, dtype=np.float64)[1:], np.array(points, dtype=np.float64)[1:], np.array(photons[1:], dtype=np.float64)

    # @profile(precision=4)
    def scintillator_monte_carlo(self, o, notabsorbed, scint_radius, scint_plane, light_guide_planes, pmt_center, pmt_radius, corner_center, corner_radius, N_max, t, keepdata):
        if keepdata: track_history = np.zeros((N_max+1,7))         # x, y history of Photon
        endpoint_dist = np.norm(o-pmt_center)
        theta = uniform(0,2*np.pi)             # first theta direction of photon
        phi = uniform(0,np.pi)                 # first phi   direction of photon
        PMT_hit_condition = False
        total_dist = 0
        dt = 0
        u = np.array([np.sin(phi)*np.cos(theta),np.sin(phi)*np.sin(theta),np.cos(phi)]) # first direction unit vector
        if keepdata: track_history[0,:] = [o[0],o[1],o[2],u[0],u[1],u[2],notabsorbed]
        i = 0
        while (i < N_max) & (not PMT_hit_condition) & (notabsorbed is True):
            ds, PMT_hit_condition = self.distance_solver(u, o, np.array([0,0,scint_plane[0]]),scint_radius, scint_plane, corner_center, corner_radius, pmt_center, pmt_radius)
            x, y, z = o+ds*u
            total_dist += np.norm(ds*u[0:2])
            o = np.array([x, y, np.abs(z) if np.abs(z-scint_plane).any() < 1e-5 else z])
            dt += np.abs(ds)/self.c                        # time taken in ps traveling in direction theta
    #         print(f"step {i}: ds={ds:.2f}cm dt={dt:.2f}ps Absorbed?={not notabsorbed} xyz =({x:.2f},{y:.2f},{z:.2f}) u=({u[0]:.2f},{u[1]:.2f},{u[2]:.2f})")
            n = self.n_vec_calculate(o, scint_plane, light_guide_planes, corner_center, corner_radius)
            u, notabsorbed = self.photon_interaction(u, n)
            if keepdata: track_history[i+1] = [x,y,z,u[0],u[1],u[2],notabsorbed]
            i+=1
        if keepdata:
            if (i < N_max):
                track_history = track_history[:i+1,:]
            # release
            endpoint_dist = None; theta = None; phi = None; total_dist = None; u = None; ds = None; x = None; y = None; z = None; 
            o = None; n = None; notabsorbed = None;
            return PMT_hit_condition, (t+dt), track_history
        else:
            # release
            endpoint_dist = None; theta = None; phi = None; total_dist = None; u = None; ds = None; x = None; y = None; z = None; 
            o = None; n = None; notabsorbed = None;
            return PMT_hit_condition, (t+dt), total_dist, endpoint_dist, i, dt

    # PMT SIMULATION
    def photontoElectrons(self, photons):
        e = 0.
        for i in range(int(photons)):
            if np.random()<self.QE: # Main Monte Carlo 
                e+=1
        for dynode in range(self.n_dynodes-1):
            delta_voltage = self.V[dynode+1]-self.V[dynode]
            e += np.poisson(e*delta_voltage/self.E_per_electron)
        return e

    #############################
    # RUN SIMULATION 
    #############################
    def particle_task(self, mult):
        return self.particle_path(t=self.t_initial+self.seperation_time*mult, phi_range_deg=self.particle_init_angle_range, T1_z=self.T1z, T1_width=self.T1_width, 
                                                T4_z=self.T4z, T4_width=self.T4_width, T1_radius=self.T1_radius, T4_radius=self.T4_radius, T1_corner=[self.T4_radius,-self.T4_radius],
                                                T4_corner=[self.T1_radius,-self.T1_radius], mean_free_path=self.mean_free_path_scints, 
                                                photons_per_E=self.photons_produced_per_MeV, prob_scint=self.pr_of_scintillation)
    def scint_taskT1(self, point_i, time_i):
        return self.scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=self.T1_radius, 
                                                        scint_plane=np.array([self.T1z,self.T1top]),  
                                                        light_guide_planes=[self.T1_radius,-self.T1_radius], 
                                                        pmt_center=self.PMT1_center, pmt_radius=self.PMT1_radius, corner_center=self.T1_corner_center,
                                                        corner_radius=self.T1_corner_radius, N_max=self.max_simulated_reflections, t=time_i, keepdata=False)
    def scint_taskT4(self, point_i, time_i):
        return self.scintillator_monte_carlo(point_i, notabsorbed=True, scint_radius=self.T4_radius, 
                                                        scint_plane=np.array([self.T4z,self.T4top]),
                                                        light_guide_planes=[-self.T4_radius,+self.T4_radius], 
                                                        pmt_center=self.PMT4_center, pmt_radius=self.PMT4_radius, corner_center=self.T4_corner_center,
                                                        corner_radius=self.T4_corner_radius, N_max=self.max_simulated_reflections, t=time_i, keepdata=False)
   
    def run_worker_T1(self, i, q):
        with h5py.File('temp.hdf5', 'r') as f:
            data = f['T1']
            point = data['points'][i]
            time = data['times'][i]
            particle_id = data['particleID'][i]
        # Move the Q out of the scope of the file 
        hit, travel_time, prop_dist, endpt_dist, prop_time, interactions = self.scint_taskT1(point, time)
        q.put([hit, travel_time, prop_dist, endpt_dist, prop_time, interactions, particle_id])
        hit = None; travel_time = None; prop_dist = None; endpt_dist = None; prop_time = None; interactions = None;
        
    def run_worker_T4(self, i, q):
        with h5py.File('temp.hdf5', 'r') as f:
            data = f['T4']
            point = data['points'][i]
            time = data['times'][i]
            particle_id = data['particleID'][i]
        # Move the Q out of the scope of the file 
        hit, travel_time, prop_dist, endpt_dist, prop_time, interactions = self.scint_taskT4(point, time)
        q.put([hit, travel_time, prop_dist, endpt_dist, prop_time, interactions, particle_id])
        hit = None; travel_time = None; prop_dist = None; endpt_dist = None; prop_time = None; interactions = None;

    # def listener(self, q, filename):
    #     '''listens for messages on the q, writes to file. '''
    #     f = h5py.File(f'{str(filename)}.hdf5', 'w')
    #     not_created = True
    #     while 1:
    #         new_data = q.get()
    #         if new_data == 'kill':
    #             f['data'].resize((f['data'].attrs['n_photons']), axis=0)
    #             f.close()
    #             print("Queue", filename, "finished!")
    #             break
    #         if not_created:
    #             ds = f.create_dataset('data', data=new_data, dtype='float64', compression="gzip", chunks=True, shape=(1,7), maxshape=(None,7))
    #             ds.attrs['n_photons'] = 0
    #             not_created = False
    #         else:
    #             if f['data'].attrs['n_photons'] == f['data'].shape[0]:
    #                 # if out of space add 10 rows
    #                 f['data'].resize((f['data'].shape[0] + 10), axis=0)
    #             # add data regardless and increase counter
    #             f['data'][f['data'].attrs['n_photons'],:] = new_data
    #             f['data'].attrs['n_photons'] += 1
    def listener(self, q, filename, chunk_size=3000):
        '''listens for messages on the q, writes to file. '''
        f = h5py.File(f'{str(filename)}.hdf5', 'w')
        not_created = True
        chunk = []
        while 1:
            new_data = q.get()
            if new_data == 'kill':
                if chunk:
                    if not_created:
                        ds = f.create_dataset('data', data=chunk, dtype='float64', compression="gzip", chunks=True, shape=(len(chunk), 7), maxshape=(None, 7))
                        ds.attrs['n_photons'] = len(chunk)
                        not_created = False
                    else:
                        if f['data'].attrs['n_photons'] + len(chunk) > f['data'].shape[0]:
                            # Resize the dataset if needed
                            f['data'].resize((f['data'].attrs['n_photons'] + len(chunk)), axis=0)
                        f['data'][f['data'].attrs['n_photons']:(f['data'].attrs['n_photons'] + len(chunk)), :] = chunk
                        f['data'].attrs['n_photons'] += len(chunk)
                f.close()
                print("Queue", filename, "finished!")
                break

            # Accumulate data in the chunk list until the chunk size is reached
            chunk.append(new_data)
            if len(chunk) >= chunk_size:
                if not_created:
                    ds = f.create_dataset('data', data=chunk, dtype='float64', compression="gzip", chunks=True, shape=(len(chunk), 7), maxshape=(None, 7))
                    ds.attrs['n_photons'] = len(chunk)
                    not_created = False
                else:
                    if f['data'].attrs['n_photons'] + len(chunk) > f['data'].shape[0]:
                        # Resize the dataset if needed
                        f['data'].resize((f['data'].attrs['n_photons'] + len(chunk)), axis=0)
                    f['data'][f['data'].attrs['n_photons']:(f['data'].attrs['n_photons'] + len(chunk)), :] = chunk
                    f['data'].attrs['n_photons'] += len(chunk)
                chunk = []


    # @profile(precision=4)
    def run(self, *arg, **kwargs):
        """Run simulation with default 1 particle or arg[0] as number of particles and a time seperation of 'delta_t'=1e-5"""
        import gc
        freeze_support()
        if arg:
            self.num_particles = int(arg[0])
            print(f"Generating {self.num_particles} particles now...")
        else:
            self.num_particles = 1
            print(f"Generating {self.num_particles} particle now...")
        self.seperation_time = kwargs.get('delta_t', self.seperation_time) # in ps
        logstarttime = perf_counter_ns()
        # FIND PARTICLE PATH
        times = []
        points = []
        photons = []
        particleID = []
        i = 0
        with Pool(processes=cpu_count()-1) as pool:
            res = pool.map(self.particle_task, range(self.num_particles))
            for (time_i, point_i, photon_i) in res:
                i = 0
                times.extend(time_i)
                points.extend(point_i)
                photons.extend(photon_i)
                particleID.extend(np.repeat(i, len(time_i))) # particle it belongs to
                i += 1
        logendparticle = perf_counter_ns()
        N = np.sum(photons)
        print("Photons generated", N)
        times = np.asarray(times); points = np.asarray(points); photons = np.asarray(photons); particleID = np.asarray(particleID)
        T1_count = np.sum(photons[points[:,2] >= self.T1z]).astype(int)
        T4_count = np.sum(photons[points[:,2] < self.T1z]).astype(int)
        with h5py.File('temp.hdf5', 'w') as f:
            print(f"Photons in T1: {T1_count} and Photons in T4: {T4_count}")
            t1 = f.create_group("T1")
            t1.create_dataset("times", data=np.repeat(times[points[:,2] >= self.T1z], photons[points[:,2] >= self.T1z].astype(int), axis=0), dtype=np.float64)
            t1.create_dataset("points", data=np.repeat(points[points[:,2] >= self.T1z], photons[points[:,2] >= self.T1z].astype(int), axis=0), dtype=np.float64)
            t1.create_dataset("particleID", data=np.repeat(particleID[points[:,2] >= self.T1z], photons[points[:,2] >= self.T1z].astype(int), axis=0), dtype=np.float64)
            t4 = f.create_group("T4")
            t4.create_dataset("times", data=np.repeat(times[points[:,2] < self.T1z],photons[points[:,2] < self.T1z].astype(int), axis=0), dtype=np.float64)
            t4.create_dataset("points", data=np.repeat(points[points[:,2] < self.T1z],photons[points[:,2] < self.T1z].astype(int), axis=0), dtype=np.float64)
            t4.create_dataset("particleID", data=np.repeat(particleID[points[:,2] < self.T1z],photons[points[:,2] < self.T1z].astype(int), axis=0), dtype=np.float64)
        T1_total = len(times[points[:,2] >= self.T1z]) * T1_count
        T4_total = len(times[points[:,2] < self.T1z]) * T4_count
        del times; del points; del photons; del particleID
        gc.collect()
        logstartphoton = perf_counter_ns()
        
        # New write and collect executor
        manager = Manager()
        q1 = manager.Queue()
        q4 = manager.Queue()
        with Pool(processes=cpu_count() -1 ) as pool:
            #put listeners to work first
            print("Created file t1_data.hdf5 with size", T1_total)
            watcher_t1 = pool.apply_async(self.listener, (q1, 't1_data'))
            print("Created file t4_data.hdf5 with size", T4_total)
            watcher_t4 = pool.apply_async(self.listener, (q4, 't4_data'))

            #fire off workers
            jobs = []
            print("T1 Photon Propagation working...")
            for i in range(T1_count-1):
                job = pool.apply_async(self.run_worker_T1, (i,q1))
                jobs.append(job)
            
            print("T4 Photon Propagation working...")
            for i in range(T4_count-1):
                job = pool.apply_async(self.run_worker_T4, (i,q4))
                jobs.append(job)
            
            # Collect results
            for j in tqdm(jobs):
                j.get()

            print("Done!")

            # once collected kill the queues
            q1.put('kill')
            q4.put('kill')

            watcher_t1.get()
            watcher_t4.get()

        logendtime = perf_counter_ns()
        # LOAD RESULTS
        f_t1 = h5py.File('t1_data.hdf5', 'r')
        f_t4 = h5py.File('t4_data.hdf5', 'r')
        
        # PRINT RESULTS
        print("TIME ANALYSIS:")
        pgtime = timedelta(seconds=(logendparticle-logstarttime)/1e9)
        phtime = timedelta(seconds=(logendtime-logstartphoton)/1e9)
        ttime = timedelta(seconds=(logendtime-logstarttime)/1e9)
        print(f"Generation of Particles     {str(pgtime)}")
        print(f"Simulation of Photon Travel {str(phtime)}")
        print(f"Total Time Elapsed:         {str(ttime)}")
        print("RESULTS SUMMARY:")
        print("HITS on T1", np.sum(f_t1['data'][:,0]))
        if np.sum(f_t1['data'][:,0]) > 0:
            print("RATIO T1   total photons", T1_count, "total incident photons", np.sum(f_t1['data'][:,0]), f"ratio={T1_count/np.sum(f_t1['data'][:,0]):.2f}")
        print("HITS on T4", np.sum(f_t4['data'][:,0]))
        if np.sum(f_t4['data'][:,0]) > 0:
            print("RATIO T4   total photons ", T4_count,"total incident photons", np.sum(f_t4['data'][:,0]), f"ratio={T1_count/np.sum(f_t4['data'][:,0]):.2f}")
        # BEGIN SIMULATING PMT PULSE
        signals_channelT1 = []
        signals_channelT4 = []
        output_times_channelT1 = []
        output_times_channelT4 = []
        signals = []
        for t in f_t1['data'][(f_t1['data'][:,0] == 1)][:,1]:
            pmtSignal_i = self.photontoElectrons(1)
            output_times_channelT1.append(self.pmt_electron_travel_time+t)
            signals.append(pmtSignal_i)
            signals_channelT1.append(pmtSignal_i)
        for t in f_t4['data'][(f_t4['data'][:,0] == 1)][:,1]:
            pmtSignal_i = self.photontoElectrons(1)
            output_times_channelT4.append(self.pmt_electron_travel_time+t)
            signals.append(pmtSignal_i)
            signals_channelT4.append(pmtSignal_i)

        # CONVERTION Electron count to Current and save in array
        self.signals = np.array(signals) * self.q / 1e-12 * self.artificial_gain # divided by 1ps 
        self.signals_channelT1 = np.array(signals_channelT1) * self.q / 1e-12 * self.artificial_gain
        self.signals_channelT4 = np.array(signals_channelT4) * self.q / 1e-12 * self.artificial_gain * 0.6 # factor to limit pulses to 50miliamps and stop contant comparator firing. however, current should be smaller from Quantum Efficiency and current should be larger from 3kV potential difference across PMT dynodes instead of current 1kV potential difference
        self.output_times_channelT1 = np.array(output_times_channelT1)
        self.output_times_channelT4 = np.array(output_times_channelT4)
        print(self.output_times_channelT1)
        print(self.output_times_channelT4)
    # Output function
    def to_csv(self, **kwargs):
        from scipy.stats import norm
        output_extra = kwargs.get('extra_data_only', False)
        output_both = kwargs.get('output_both', False)
        # OUTPUT FORMATTING
        if output_extra or output_both:
            # data index lookup:
            # 0  ,     1      ,    2      ,    3      ,    4    ,      5     ,     6                  
            # hit, travel_time, prop_dist, endpt_dist, prop_time, interactions, particleID
            file_t1 = h5py.File('t1_data.hdf5', 'r')['data']
            file_t4 = h5py.File('t4_data.hdf5', 'r')['data']
            f_t1 = np.array([data for data in file_t1 if data[0]])
            f_t4 = np.array([data for data in file_t4 if data[0]])
            # print(f_t1.shape, f_t4.shape)
            # print(self.output_times_channelT1.shape, self.output_times_channelT4.shape)
            
            print("Exporting Extra Data...")
            if f_t1.shape[0] > 0:
                dft1 = DataFrame({'T1_part_ids':f_t1[:,6],'time':f_t1[:,1],'T1_prop_dist':f_t1[:,2],'T1_endpoint_dist':f_t1[:,3], 'T1_prop_times':f_t1[:,4], 'T1_interactions':f_t1[:,5]})
                dft1.to_csv('monte_carlo_extradata'+str(self.num_particles)+'chT1_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
            else:
                print("WARN: Not enough PMT hits in T1! (< 1)")
            if f_t4.shape[0] > 0:
                dft4 = DataFrame({'T4_part_ids':f_t4[:,6],'time':f_t4[:,1],'T4_prop_dist':f_t4[:,2],'T4_endpoint_dist':f_t4[:,3], 'T4_prop_times':f_t4[:,4], 'T4_interactions':f_t4[:,5]})
                dft4.to_csv('monte_carlo_extradata'+str(self.num_particles)+'chT4_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
            else:
                print("WARN: Not enough PMT hits in T4! (< 1)")
            if not output_both:
                return
        print("Exporing to 2 channels...")
        # for each channel
        for time,signal,ch in zip([self.output_times_channelT1,self.output_times_channelT4],[self.signals_channelT1,self.signals_channelT4],[1,4]):

            # from io import StringIO
            # from csv import writer 
            # output = StringIO()
            # csv_writer = writer(output)
            
            print("Smoothing Signals...")
            t_binned = [0.] # left edges of bins
            y_binned = [0.]
            for i,y in enumerate(signal):
                # print(f"i={i},t[{i}]={time[i]} y[{i}]={y}")
                lower_bound = max(time[i]-2*self.sigma_smoothing,0) # 2 sigma away backward
                upper_bound = min(time[i]+2*self.sigma_smoothing,max(time)+2*self.sigma_smoothing) # 2 sigma away forward
                # MAKE NEW DATA CENTERED AROUND PULSE
                if lower_bound < max(t_binned): # if already binned
                    lower_bound = t_binned[np.digitize(lower_bound, t_binned)]+self.output_bin_width/2
                cur_x = np.arange(lower_bound,upper_bound,self.output_bin_width)+self.output_bin_width/2
                # print(f"cur_x from {lower_bound}-->{upper_bound}", cur_x)
                # ADD DATA IF NEEDED
                for x in cur_x:
                    if x > max(t_binned): 
                        t_binned.append(x)
                        y_binned.append(0)
                    elif (np.digitize(x, t_binned)-1 > 0) and (np.digitize(x, t_binned) < len(t_binned)):
                        index = np.digitize(x, t_binned)
                        if abs(t_binned[index]-t_binned[index-1]) > self.output_bin_width:
                            t_binned.insert(index, x) # check if need -1 or just np.digitize()
                            y_binned.insert(index, 0) # check 
                # GET INDICIES
                index_lower = [i for i,t in enumerate(t_binned) if t >= lower_bound][0] # left edge in time binned
                index_upper = [i for i,t in enumerate(t_binned) if t <= upper_bound][-1] # right edge in time binned
                # GAUSSIAN SMOOTH
                gaussian = norm.pdf(t_binned[index_lower:index_upper], loc=time[i], scale=self.sigma_smoothing)*self.sigma_smoothing*y/4
                # ADD TO CORRECT BINS
                for i,y_add in enumerate(gaussian):
                    if y_binned[index_lower+i]+y_add < self.max_pmt_current_output:
                        y_binned[index_lower+i] += y_add
                    else:
                        y_binned[index_lower+i] = self.max_pmt_current_output

            df = DataFrame({'time':t_binned,'current':y_binned}).sort_values(by=['time'])
            print("Formatting PWL dataframe...")
            fill_data = []                                                                      # declare empty array
            # begin padding data at time 1/5th bin width before first time stamp
            fill_data.append([df['time'].iloc[0]-self.output_bin_width/5,0])                    # add zero at beginning
            for i in range(len(df['time'])-1):                                                        # for each time index
                if abs(df['time'].iloc[i]-df['time'].iloc[i+1]) > self.output_bin_width:        # if dt between signals is greater than minimum bin width
                    fill_data.append([df['time'].iloc[i]+self.output_bin_width/5,0])            # place zero after current signal
                    fill_data.append([df['time'].iloc[i+1]-self.output_bin_width/5,0])          # place zero before next signal
            fill_data.append([df['time'].iloc[-1]+self.output_bin_width/5,0])                   # add zero at end
            fill_data = np.array(fill_data)
            fill = DataFrame(fill_data, columns=['time','current'])
            df = concat([fill, df], ignore_index=True).sort_values(by=['time']).reset_index(drop=True)
            df['time'] = df['time']/1e12
            df = df[['time', 'current']] # need this for LTSpice PWL current input file to work
            df.to_csv('monte_carlo_input'+str(self.num_particles)+'ch'+str(ch)+'_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt', float_format='%.13f', header=False, index=False, sep=' ') # PWL file formatting
        print("Done!")