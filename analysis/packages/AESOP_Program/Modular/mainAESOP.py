from modAESOP import *
from modAESOP import Simulation as sim




# @profile(precision=4)
def run(sim, *arg, **kwargs):
    """Run simulation with default 1 particle or arg[0] as number of particles and a time seperation of 'delta_t'=1e-5"""
    import gc
    freeze_support()
    sim.intro() # new 
    if arg:
        sim.num_particles = int(arg[0])
        print(f"Generating {sim.num_particles} particles now...")
    else:
        sim.num_particles = 1
        print(f"Generating {sim.num_particles} particle now...")
    sim.seperation_time = kwargs.get('delta_t', sim.seperation_time) # in ps
    logstarttime = perf_counter_ns()
    # FIND PARTICLE PATH
    times = []
    points = []
    photons = []
    particleID = []
    i = 0
    with Pool(processes=cpu_count()-1) as pool:
        res = pool.map(sim.particle_task, range(sim.num_particles))
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
    T1_count = np.sum(photons[points[:,2] >= sim.T1z]).astype(int)
    T4_count = np.sum(photons[points[:,2] < sim.T1z]).astype(int)
    with h5py.File('temp.hdf5', 'w') as f:
        print(f"Photons in T1: {T1_count} and Photons in T4: {T4_count}")
        t1 = f.create_group("T1")
        t1.create_dataset("times", data=np.repeat(times[points[:,2] >= sim.T1z], photons[points[:,2] >= sim.T1z].astype(int), axis=0), dtype=np.float64)
        t1.create_dataset("points", data=np.repeat(points[points[:,2] >= sim.T1z], photons[points[:,2] >= sim.T1z].astype(int), axis=0), dtype=np.float64)
        t1.create_dataset("particleID", data=np.repeat(particleID[points[:,2] >= sim.T1z], photons[points[:,2] >= sim.T1z].astype(int), axis=0), dtype=np.float64)
        t4 = f.create_group("T4")
        t4.create_dataset("times", data=np.repeat(times[points[:,2] < sim.T1z],photons[points[:,2] < sim.T1z].astype(int), axis=0), dtype=np.float64)
        t4.create_dataset("points", data=np.repeat(points[points[:,2] < sim.T1z],photons[points[:,2] < sim.T1z].astype(int), axis=0), dtype=np.float64)
        t4.create_dataset("particleID", data=np.repeat(particleID[points[:,2] < sim.T1z],photons[points[:,2] < sim.T1z].astype(int), axis=0), dtype=np.float64)
    T1_total = len(times[points[:,2] >= sim.T1z]) * T1_count
    T4_total = len(times[points[:,2] < sim.T1z]) * T4_count
    del times; del points; del photons; del particleID
    gc.collect()
    logstartphoton = perf_counter_ns()
    
    # New write and collect executor
    manager = Manager()
    q1 = manager.Queue()
    q4 = manager.Queue()
    with Pool(processes=cpu_count()-1) as pool:
        #put listeners to work first
        print("Created file t1_data.hdf5 with size", T1_total)
        watcher_t1 = pool.apply_async(sim.listener, (q1, 't1_data'))
        print("Created file t4_data.hdf5 with size", T4_total)
        watcher_t4 = pool.apply_async(sim.listener, (q4, 't4_data'))

        chunk_size = 100
        #fire off workers
        jobs_T1 = []
        print("T1 Photon Propagation working...")
        for i in range(T1_count-1):
            job = pool.apply_async(sim.run_worker_T1, (i,q1))
            jobs_T1.append(job)

        jobs_T4 = []
        print("T4 Photon Propagation working...")
        for i in range(T4_count-1):
            job = pool.apply_async(sim.run_worker_T4, (i,q4))
            jobs_T4.append(job)
        
        # Collect results for T1
        for j in tqdm(jobs_T1):
            j.get()

        # Collect results for T4
        for j in tqdm(jobs_T4):
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
        pmtSignal_i = sim.photontoElectrons(1)
        output_times_channelT1.append(sim.pmt_electron_travel_time+t)
        signals.append(pmtSignal_i)
        signals_channelT1.append(pmtSignal_i)
    for t in f_t4['data'][(f_t4['data'][:,0] == 1)][:,1]:
        pmtSignal_i = sim.photontoElectrons(1)
        output_times_channelT4.append(sim.pmt_electron_travel_time+t)
        signals.append(pmtSignal_i)
        signals_channelT4.append(pmtSignal_i)

    # CONVERTION Electron count to Current and save in array
    sim.signals = np.array(signals) * sim.q / 1e-12 * sim.artificial_gain # divided by 1ps 
    sim.signals_channelT1 = np.array(signals_channelT1) * sim.q / 1e-12 * sim.artificial_gain
    sim.signals_channelT4 = np.array(signals_channelT4) * sim.q / 1e-12 * sim.artificial_gain * 0.6 # factor to limit pulses to 50miliamps and stop contant comparator firing. however, current should be smaller from Quantum Efficiency and current should be larger from 3kV potential difference across PMT dynodes instead of current 1kV potential difference
    sim.output_times_channelT1 = np.array(output_times_channelT1)
    sim.output_times_channelT4 = np.array(output_times_channelT4)
    print(sim.output_times_channelT1)
    print(sim.output_times_channelT4)

if __name__ == "__main__": 
    ##########################################
    # DECLARE SIMULATION AND PLOTTER CLASSES
    ##########################################
    sim = Simulation()
    plot = plotter(sim)

    #####################
    # RUN SIMULATION 
    #####################
    sim.max_simulated_reflections = 8
    sim.mean_free_path_scints = 0.0001

    # # sim.mean_free_path_scints = 0.00024 # cm -> 2.4 micrometers
    # # sim.num_particles = 4000
    run(sim,1)
    sim.to_csv(output_both=True)

    ###############################################################
    # RUN LTSPICE AND CALCULATE TIME OF FLIGHT --> SAVE TO FILE
    ###############################################################
    # sim.ltspice(filedate='07_12_2023',filenum=1)
    # sim.calc_ToF(filedate='07_12_2023',filenum=1)
    # sim.save_ToF()

    #########################################
    # LOAD CORRECTED MODEL AND PLOT EXTRA DATA
    #########################################
    plot.load_extradata(filename='monte_carlo_extradata1chT1_02_17_2024.txt')
    plot.plot_xydistance_distr()
    plot.plot_distPMT_proptime()
    # plot.load_ToF(1, filename='result_1_of_1_07_12_2023.txt')
    # plot.correct_tof()
    # sim.load_ToF(3858, filename='result_3858_of_4000_07_02_2023.txt')
    # sim.plotToF()