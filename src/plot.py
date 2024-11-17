#############################
# AESOP-Lite Monte Carlo (2024)
# Plotter
# Created by Liam Branch, Michelle Picardo, Robert Johnson
# UCSC Alumni + UCSC Professor
#############################

# import numpy as np
# import utility as np
import h5py
# from pandas import DataFrame, read_csv, concat
from tqdm import tqdm
from random import uniform
from time import perf_counter_ns
from datetime import timedelta, datetime
from multiprocessing import Pool, cpu_count, freeze_support, Manager, set_start_method
#############################
# DEBUG AND PLOTTING
#############################

class Plotter:
    def __init__(self, Simulation):
        #############################
        # CONSTANTS
        #############################
        self.c = Simulation.c
        self.q = Simulation.q
        self.n_1 = Simulation.n_1
        self.n_2 = Simulation.n_2
        self.T1z= Simulation.T1z
        self.T4z= Simulation.T4z
        self.T1_radius = Simulation.T1_radius #
        self.T4_radius = Simulation.T4_radius #
        self.T1_width = Simulation.T1_width
        self.T4_width = Simulation.T4_width
        self.xPMT4=Simulation.xPMT4
        self.yPMT4=Simulation.yPMT4
        self.xPMT1=Simulation.xPMT1
        self.yPMT1=Simulation.yPMT1
        self.PMT1_radius = Simulation.PMT1_radius #
        self.PMT4_radius = Simulation.PMT4_radius #
        self.n_dynodes = Simulation.n_dynodes
        self.V = Simulation.V
        # self.V = [150,300,350,600,750,850]
        self.E_per_electron = Simulation.E_per_electron
        self.QE = Simulation.QE
        self.t_initial = Simulation.t_initial
        self.particle_init_angle_range = Simulation.particle_init_angle_range
        self.particle_gen_area = Simulation.particle_gen_area
        self.particle_gen_z = Simulation.particle_gen_z
        self.mean_free_path_scints = Simulation.mean_free_path_scints
        self.photons_produced_per_MeV = Simulation.photons_produced_per_MeV
        self.pr_of_scintillation = Simulation.pr_of_scintillation
        self.max_simulated_reflections = Simulation.max_simulated_reflections
        self.pmt_electron_travel_time = Simulation.pmt_electron_travel_time
        self.artificial_gain = Simulation.artificial_gain
        self.pr_absorption = Simulation.pr_absorption
        self.seperation_time = Simulation.seperation_time
        self.output_bin_width = Simulation.output_bin_width
        self.num_particles = Simulation.num_particles
        self.CMOS_thresh = Simulation.CMOS_thresh
        self.reemission_angle_factor = Simulation.reemission_angle_factor

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
        T1extra_data = read_csv(fileT1)
        T4extra_data = read_csv(fileT4)
        # If require to replace, replace
        if 'T1_part_ids' in T1extra_data and 'T4_part_ids' in T4extra_data:
            self.T1_part_ids = T1extra_data['T1_part_ids']
            self.T4_part_ids = T4extra_data['T4_part_ids']
        if 'time' in T1extra_data and 'time' in T4extra_data:
            self.T1_extratimes = T1extra_data['time']
            self.T4_extratimes = T4extra_data['time']
        self.T1_prop_dist = T1extra_data['T1_prop_dist']
        self.T4_prop_dist = T4extra_data['T4_prop_dist']
        self.T1_endpoint_dist = T1extra_data['T1_endpoint_dist']
        self.T4_endpoint_dist = T4extra_data['T4_endpoint_dist']
        self.T1_prop_times = T1extra_data['T1_prop_times']
        self.T4_prop_times = T4extra_data['T4_prop_times']
        self.T1_interactions = T1extra_data['T1_interactions']
        self.T4_interactions = T4extra_data['T4_interactions']
        print("Saved to class fields")

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
            h0 = ax1[0].hist2d(self.T1_prop_times, self.T1_interactions, bins=[250,40], range=[[0,np.mean(self.T1_prop_times)+np.std(self.T1_prop_times)*5],[0,max(self.T1_interactions)]], cmin = 1)
            fig1.colorbar(h0[3], ax=ax1[0])
            # sns.kdeplot(self.T1_prop_dist, self.T1_interactions, n_levels=250, cbar=True, shade_lowest=False, cmap='viridis', ax=ax1[0])
            # ax1[0].axvline(self.T1_radius, color='C2', label='T1 radius')
            ax1[0].set_xlabel('Time [ps]')
            ax1[0].set_ylabel('# Interactions with Boundary')
            ax1[0].set_title('T1 Propagation Times vs. Interactions / Reflections')
            ax1[0].legend()
            ax1[0].grid()
            # ax1[1].scatter(self.T4_prop_dist, self.T4_interactions, s=10, alpha=0.3, facecolors='none', edgecolors='C0')
            h1 = ax1[1].hist2d(self.T4_prop_times, self.T4_interactions, bins=[250,40], range=[[0,np.mean(self.T4_prop_times)+np.std(self.T4_prop_times)*5],[0,max(self.T4_interactions)]], cmin = 1)
            fig1.colorbar(h1[3], ax=ax1[1])
            # sns.kdeplot(self.T4_prop_dist, self.T4_interactions, n_levels=250, cbar=True, shade_lowest=False, cmap='viridis', ax=ax1[1])
            # ax1[1].axvline(self.T4_radius, color='C2',label='T4 radius')
            ax1[1].set_xlabel('Time [ps]')
            ax1[1].set_ylabel('# Interactions with Boundary')
            ax1[1].set_title('T4 Propagation Times vs. Interactions / Reflections')
            ax1[1].legend()
            ax1[1].grid()
            plt.show()
        else:
            print("Need to run simulation first with sim.run()")

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
        for x1,y1,z1,u1,v1,w1,l in zip(tracks[:-1,0],tracks[:-1,1],tracks[:-1,2], tracks[:-1,3], tracks[:-1,4],tracks[:-1,5],np.norm(tracks[:-1,3:6],ord=2, axis=1)):
            ax.quiver(x1, y1, z1, u1, v1, w1, length=l*0.5, edgecolor='k', facecolor='black')
        ax.plot(tracks[:,0],tracks[:,1],tracks[:,2], alpha=0.5, color='C1', marker='.')
        ax.text(photon_pos[0],photon_pos[1],photon_pos[2], f'hitPMT1? {hit_PMT == 1}')
        plt.show()
    
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
        new_ToF_data = read_csv(file)
        print(new_ToF_data)
        # If require to replace, replace
        if replace is not False:
            self.FinalToF = new_ToF_data.to_numpy()
        else:
            self.FinalToF = np.append(self.FinalToF, new_ToF_data.to_numpy())
    
    def plotToF(self):
        if not hasattr(self, 'FinalToF'):
            print("Need to run calc_ToF or load_ToF result file")
            return
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit
        # Setup least squares loss for guassian function with 4 parameters (amplitude, mean, stddev, bias)
        fitfunc  = lambda x, p0, p1, p2, p3: p0*np.exp(-0.5*((x-p1)/p2)**2)+p3
        errfunc  = lambda p, x, y: (y - fitfunc(p, x))
        init  = [200e0, 2e-9, 1e-9, 0e0]
        # plt.title(f'TOF, Total Points {len(self.FinalToF)} {self.max_simulated_reflections} maximum reflections')
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
        plt.xlabel('time [s]', fontsize=14)
        plt.ylabel('Counts', fontsize=14)
        plt.grid()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=14)
        plt.show()

    def correct_tof(self):
        if not hasattr(self, 'FinalToF'): # assume final tof found
            print("Need to import ToF result file!")
            return
        if not hasattr(self, 'T1_prop_times') or not hasattr(self, 'T4_prop_times'): # assume extradata loaded
            print("Need to import ToF extradata files!")
            return
        self.plotToF()
        import matplotlib.pyplot as plt
        fig0, ax0 = plt.subplots()
        up_bound = np.quantile(self.FinalToF, 0.95) # should be already defined in plotToF()
        low_bound = np.quantile(self.FinalToF, 0.05)
        # Plot histogram
        # DEFINE PMT4dist, PMT1dist
        # correction = self.FinalToF + self.n_2(PMT4dist - PMT1dist) 
        correction = np.zeros(len(self.FinalToF))
        for i in range(len(self.FinalToF)):
            propagation_mean = np.mean(self.T1_prop_times[self.T1_part_ids == i])/1e12
            correction[i] = self.FinalToF[i]-propagation_mean
        ax0.hist(self.FinalToF, bins=np.linspace(low_bound-1e-9,up_bound+1e-9,125), histtype='step', edgecolor='r', linewidth=2)
        ax0.hist(correction, bins=np.linspace(low_bound-1e-9,up_bound+1e-9,125), histtype='step', edgecolor='g', linewidth=2)
        ax0.set_xlabel('time [s]', fontsize=14)
        ax0.set_ylabel('Counts', fontsize=14)
        ax0.grid()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax0.legend(fontsize=14)
        plt.show()

