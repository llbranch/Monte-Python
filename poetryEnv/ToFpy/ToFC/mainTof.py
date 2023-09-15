#############################
# Testing .close/.join
# without: 
#############################

# Installs Needed 
import numpy as np #analysis
import pandas as pd #analysis
from tqdm import tqdm #progressbar
# Native 
import random #generator 
from time import perf_counter #timing 
from datetime import timedelta, datetime #timing
from memory_profiler import profile, LogFile #memory checks

# Testing Multiprocessing 
from multiprocessing import Pool, cpu_count, freeze_support
# from concurrent.futures import ProcessPoolExecutor

# Native
import cProfile # tuna 
import pstats
import sys


# Module Pull
from modTof import *


#######################################  P O O L #############################################################################
"""Run simulation with default 1 particle or arg[0] as number of particles and a time seperation of 'delta_t'=1e-5"""
fp = open("memory_profiler_test.log", "w+")
@profile(precision=4, stream=fp)
# @profile(precision =4)
def run(simClass, *arg, **kwargs):
    import gc
    freeze_support() # best practice 
    if arg:
        simClass.num_particles = int(arg[0])
        print(f"Generating {simClass.num_particles} particles now...")
    else:
        simClass.num_particles = 1
        print(f"Generating {simClass.num_particles} particle now...")
    simClass.seperation_time = kwargs.get('delta_t', simClass.seperation_time) # in ps
    logstarttime = perf_counter()
    # FIND PARTICLE PATH
    times = []
    points = []
    photons = []
    particleID = []
    i = 0
    with Pool(processes=cpu_count(), maxtasksperchild=2) as pool:
        res = pool.map(simClass.particle_task, range(simClass.num_particles))

        pool.close()
        pool.join()


        # print(f'SIZE OF RES: ', sys.getsizeof(res))
        # start = perf_counter()

        for (time_i, point_i, photon_i) in res:
            i = 0
            times.extend(time_i)
            points.extend(point_i)
            photons.extend(photon_i)
            particleID.extend(np.repeat(i, len(time_i))) # particle it belongs to
            i += 1

        # end = perf_counter() -start

    logendparticle = perf_counter()
    N = np.sum(photons)
    print("Photons generated", N)
    times = np.asarray(times); points = np.asarray(points); photons = np.asarray(photons); particleID = np.asarray(particleID)
    del(N)
    # print(f'SIZE OF TIMES, POINTS, ETC: ', sys.getsizeof(times))
    # print('Time to unpack times,points, etc.', end)


    # RETURNS A FILE
    # SPLIT HERE
    # RUN #2
    

    # SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
    # Gather TOF data
    T1_input_times = []
    T4_input_times = []
    # Gather Extra Data for analysis
    simClass.T1_prop_dist = []
    simClass.T4_prop_dist = []
    simClass.T1_endpoint_dist = []
    simClass.T4_endpoint_dist = []
    simClass.T1_prop_times = []
    simClass.T4_prop_times = []
    simClass.T1_interactions = []
    simClass.T4_interactions = []
    simClass.T1_part_ids = []
    simClass.T4_part_ids = []
    T1points = (points[:])[points[:,2] >= simClass.T1z]
    T1times = (times[:])[points[:,2] >= simClass.T1z]
    T1photons = (photons[:])[points[:,2] >= simClass.T1z]
    T1part_ids = (particleID[:])[points[:,2] >= simClass.T1z]
    T1part_ids = np.repeat(T1part_ids, T1photons.astype(int), axis=0) # big id bank
    T4points = (points[:])[points[:,2] < simClass.T1z]
    T4times = (times[:])[points[:,2] < simClass.T1z]
    T4photons = (photons[:])[points[:,2] < simClass.T1z]
    T4part_ids = (particleID[:])[points[:,2] < simClass.T1z]
    T4part_ids = np.repeat(T4part_ids, T4photons.astype(int), axis=0) # big id bank
    print(f"Photons in T1: {np.sum(T1photons)} and Photons in T4: {np.sum(T4photons)}")
    del times; del points; del photons; # remove copies
    # gc.collect()

    logstartphoton = perf_counter()

    # check this link https://stackoverflow.com/questions/14749897/python-multiprocessing-memory-usage
    with Pool(processes=cpu_count(), maxtasksperchild=2) as pool: # this way of making the pool causes all the data to copy! 
        print("T1 Photon Propagation working...")

        # start = perf_counter()
        T1res = pool.starmap(simClass.scint_taskT1, np.repeat(np.c_[T1points,T1times],T1photons.astype(int), axis=0))
        # end = perf_counter() - start
        # print("Time to process T1:", end)

        print("Done!")


        print("T4 Photon Propagation working...")
        # start = perf_counter()
        T4res = pool.starmap(simClass.scint_taskT4, np.repeat(np.c_[T4points,T4times],T4photons.astype(int), axis=0))
        # end = perf_counter() - start
        # print("Time to process T4:", end)

        print("Done!")
        print("Unzipping reuslts into arrays...")

        pool.close()
        pool.join()

        # start = perf_counter()
        for (T1hit_PMT, T1travel_time, T1tot_dist, T1endpt, T1bounces, T1prop),T1part_id in zip(T1res, T1part_ids): # check if moving starmap here helps
            if T1hit_PMT:
                T1_input_times.append(T1travel_time)
                simClass.T1_prop_dist.append(T1tot_dist)
                simClass.T1_endpoint_dist.append(T1endpt)
                simClass.T1_prop_times.append(T1prop)
                simClass.T1_interactions.append(T1bounces)
                simClass.T1_part_ids.append(T1part_id)
        # end = perf_counter() - start
        # print("Time to unpack T1:", end)

        # start = perf_counter()
        for (T4hit_PMT, T4travel_time, T4tot_dist, T4endpt, T4bounces, T4prop),T4part_id in zip(T4res, T4part_ids): # check if moving starmap here helps
            if T4hit_PMT:
                T4_input_times.append(T4travel_time)
                simClass.T4_prop_dist.append(T4tot_dist)
                simClass.T4_endpoint_dist.append(T4endpt)
                simClass.T4_prop_times.append(T4prop)
                simClass.T4_interactions.append(T4bounces)
                simClass.T4_part_ids.append(T4part_id)
        # end = perf_counter() - start
        # print("Time to unpack T4:", end)

    logendtime = perf_counter()
    # PRINT RESULTS
    print("TIME ANALYSIS:")
    pgtime = timedelta(seconds=logendparticle-logstarttime)
    phtime = timedelta(seconds=logendtime-logstartphoton)
    ttime = timedelta(seconds=logendtime-logstarttime)
    print(f"Generation of Particles     {str(pgtime)}")
    print(f"Simulation of Photon Travel {str(phtime)}")
    print(f"Total Time Elapsed:         {str(ttime)}")
    print("RESULTS SUMMARY:")
    print("HITS on T1",len(T1_input_times))
    print("RATIO T1   total photons", np.sum(T1photons), "total incident photons", len(T1_input_times), f"ratio={np.sum(T1photons)/len(T1_input_times):.2f}")
    print("HITS on T4",len(T4_input_times))
    print("RATIO T4   total photons ", np.sum(T4photons),"total incident photons", len(T4_input_times), f"ratio={np.sum(T4photons)/len(T4_input_times):.2f}")
    print("DISTANCE: ")
    del T1points; del T1times; del T1photons; del T4points; del T4times; del T4photons; # remove unused variables
    # gc.collect()

    # print(T4_input_times)
    # BEGIN SIMULATING PMT PULSE
    signals_channelT1 = []
    signals_channelT4 = []
    output_times_channelT1 = []
    output_times_channelT4 = []
    signals = []
    output_times = []
    for t in T1_input_times:
        pmtSignal_i = simClass.photontoElectrons(1)
        output_times.append(simClass.pmt_electron_travel_time+t)
        output_times_channelT1.append(simClass.pmt_electron_travel_time+t)
        signals.append(pmtSignal_i)
        signals_channelT1.append(pmtSignal_i)
    for t in T4_input_times:
        pmtSignal_i = simClass.photontoElectrons(1)
        output_times.append(simClass.pmt_electron_travel_time+t)
        output_times_channelT4.append(simClass.pmt_electron_travel_time+t)
        signals.append(pmtSignal_i)
        signals_channelT4.append(pmtSignal_i)

    # CONVERTION Electron count to Current and save in array
    simClass.signals = np.array(signals) * simClass.q / 1e-12 * simClass.artificial_gain # divided by 1ps 
    simClass.output_times = np.array(output_times)
    simClass.signals_channelT1 = np.array(signals_channelT1) * simClass.q / 1e-12 * simClass.artificial_gain
    simClass.signals_channelT4 = np.array(signals_channelT4) * simClass.q / 1e-12 * simClass.artificial_gain * 0.6 # factor to limit pulses to 50miliamps and stop contant comparator firing. however, current should be smaller from Quantum Efficiency and current should be larger from 3kV potential difference across PMT dynodes instead of current 1kV potential difference
    simClass.output_times_channelT1 = np.array(output_times_channelT1)
    simClass.output_times_channelT4 = np.array(output_times_channelT4)

    # Output function
    # def to_csv(simClass, **kwargs):
    from scipy.stats import norm # type: ignore 
    output_extra = kwargs.get('extra_data_only', False)
    output_both = kwargs.get('output_both', False)
    # OUTPUT FORMATTING
    if output_extra or output_both:
        print("Exporting Extra Data...")
        dft1 = pd.DataFrame({'T1_part_ids':simClass.T1_part_ids,'time':simClass.output_times_channelT1,'T1_prop_dist':simClass.T1_prop_dist,'T1_endpoint_dist':simClass.T1_endpoint_dist, 'T1_prop_times':simClass.T1_prop_times, 'T1_interactions':simClass.T1_interactions})
        dft4 = pd.DataFrame({'T4_part_ids':simClass.T4_part_ids,'time':simClass.output_times_channelT4,'T4_prop_dist':simClass.T4_prop_dist,'T4_endpoint_dist':simClass.T4_endpoint_dist, 'T4_prop_times':simClass.T4_prop_times, 'T4_interactions':simClass.T4_interactions})
        dft1.to_csv('monte_carlo_extradata'+str(simClass.num_particles)+'chT1_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
        dft4.to_csv('monte_carlo_extradata'+str(simClass.num_particles)+'chT4_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
        if not output_both:
            return
    print("Exporing to 2 channels...")
    # for each channel
    for time,signal,ch in zip([simClass.output_times_channelT1,simClass.output_times_channelT4],[simClass.signals_channelT1,simClass.signals_channelT4],[1,4]):

        # from io import StringIO
        # from csv import writer 
        # output = StringIO()
        # csv_writer = writer(output)
        
        print("Smoothing Signals...")
        t_binned = [0.] # left edges of bins
        y_binned = [0.]
        for i,y in enumerate(signal):
            # print(f"i={i},t[{i}]={time[i]} y[{i}]={y}")
            lower_bound = max(time[i]-2*simClass.sigma_smoothing,0) # 2 sigma away backward
            upper_bound = min(time[i]+2*simClass.sigma_smoothing,max(time)+2*simClass.sigma_smoothing) # 2 sigma away forward
            # MAKE NEW DATA CENTERED AROUND PULSE
            if lower_bound < max(t_binned): # if already binned
                lower_bound = t_binned[np.digitize(lower_bound, t_binned)]+simClass.output_bin_width/2
            cur_x = np.arange(lower_bound,upper_bound,simClass.output_bin_width)+simClass.output_bin_width/2
            # print(f"cur_x from {lower_bound}-->{upper_bound}", cur_x)
            # ADD DATA IF NEEDED
            for x in cur_x:
                if x > max(t_binned): 
                    t_binned.append(x)
                    y_binned.append(0)
                elif (np.digitize(x, t_binned)-1 > 0) and (np.digitize(x, t_binned) < len(t_binned)):
                    index = np.digitize(x, t_binned)
                    if abs(t_binned[index]-t_binned[index-1]) > simClass.output_bin_width:
                        t_binned.insert(index, x) # check if need -1 or just np.digitize()
                        y_binned.insert(index, 0) # check 
            # GET INDICIES
            index_lower = [i for i,t in enumerate(t_binned) if t >= lower_bound][0] # left edge in time binned
            index_upper = [i for i,t in enumerate(t_binned) if t <= upper_bound][-1] # right edge in time binned
            # GAUSSIAN SMOOTH
            gaussian = norm.pdf(t_binned[index_lower:index_upper], loc=time[i], scale=simClass.sigma_smoothing)*simClass.sigma_smoothing*y/4
            # ADD TO CORRECT BINS
            for i,y_add in enumerate(gaussian):
                if y_binned[index_lower+i]+y_add < simClass.max_pmt_current_output:
                    y_binned[index_lower+i] += y_add
                else:
                    y_binned[index_lower+i] = simClass.max_pmt_current_output

        df = pd.DataFrame({'time':t_binned,'current':y_binned}).sort_values(by=['time'])
        print("Formatting PWL dataframe...")
        fill_data = []                                                                      # declare empty array
        # begin padding data at time 1/5th bin width before first time stamp
        fill_data.append([df['time'].iloc[0]-simClass.output_bin_width/5,0])                    # add zero at beginning
        for i in range(len(df['time'])-1):                                                        # for each time index
            if abs(df['time'].iloc[i]-df['time'].iloc[i+1]) > simClass.output_bin_width:        # if dt between signals is greater than minimum bin width
                fill_data.append([df['time'].iloc[i]+simClass.output_bin_width/5,0])            # place zero after current signal
                fill_data.append([df['time'].iloc[i+1]-simClass.output_bin_width/5,0])          # place zero before next signal
        fill_data.append([df['time'].iloc[-1]+simClass.output_bin_width/5,0])                   # add zero at end
        fill_data = np.array(fill_data)
        fill = pd.DataFrame(fill_data, columns=['time','current'])
        df = pd.concat([fill, df], ignore_index=True).sort_values(by=['time']).reset_index(drop=True)
        df['time'] = df['time']/1e12
        df = df[['time', 'current']] # need this for LTSpice PWL current input file to work
        df.to_csv('monte_carlo_input'+str(simClass.num_particles)+'ch'+str(ch)+'_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt', float_format='%.13f', header=False, index=False, sep=' ') # PWL file formatting
    print("Done!")
#######################################  P O O L #############################################################################



if __name__ == '__main__':

    with cProfile.Profile() as profile:

        sim = Simulation()

        sim.max_simulated_reflections = 8

        run(sim,1)

    results = pstats.Stats(profile)
    results.sort_stats(pstats.SortKey.TIME)
    # results.print_stats()
    results.dump_stats("report_time.prof")
    sys.stdout = LogFile('memory_profile_log',reportIncrementFlag=False)


    # sim = Simulation()
    # sim.max_simulated_reflections = 8
    # run(sim,1)
