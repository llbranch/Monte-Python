# Testing changes to Main AESOP file 
- 8 reflections 
- 1 particle 
- run times will always vary due to the CPU allocating resources on it's own 

Stopped timing analysis due to the many variables involved. If conditions are not kept the same then the run times are affected. 
Watching Youtube while running the sim added ~5 seconds to the run time. 

> Run Tests with cProfile ; needs tuna package to view it online

1. add to code 
```
import cProfile
import pstats
.
.
.

if __name__ == '__main__':

    with cProfile.Profile() as profile:

        sim = Simulation()

        sim.max_simulated_reflections = 8
        run(sim,1)

    results = pstats.Stats(profile)
    results.sort_stats(pstats.SortKey.TIME)
    # results.print_stats()
    results.dump_stats("report_time.prof")

```
2. terminal commands to run sim and call time analysis 
```
poetry run python mainTof.py 
poetry run tuna report_time.prof 
```

> View Memory with memory_profiler; will slow down code 
1. add this code before the function of interest 
```
fp = open("report_mem.log", "w+")
@profile(precision=4, stream=fp)
.
.

if __name__ == '__main__':
.
.
sys.stdout = LogFile('memory_profile_log',reportIncrementFlag=False)
```

> Saving plot 
```
poetry run mprof plot -o plot.png --backend agg
```

## Timer
See jupyter file to see results of perf_counter vs perf_counter_ns

### Run1: no changes to `AESOP_OG.py`: carbon copy of Liam's code 
Output: Program did not finish 
```
######################################################
Generated Apparatus Simulation with following defaults
######################################################
PARTICLE: Mean Free Path = 0.00024 cm
PARTICLE: Time Seperation between sequential Particles if simulation more than 1 = 100000.0
SCINT:    Probability of Scintillaton = 0.8
PMT:      Quantum Efficiency is set to 1 by default to keep more pulses
PMT:      Energy per Photoelectron is set to 20 by best estimation
PMT:      Artificial Gain on Output Current = 1
OUTPUT:   Binning Width for PWL output file = 100 ps

Run with .run() function given optional arguments below
integer n particles, 'delta_t' = 100000.0 ps particle time seperation
Generating 1 particles now...
Photons generated 49398.0
Photons in T1: 16226.0 and Photons in T4: 33172.0
T1 Photon Propagation working...
Done!
T4 Photon Propagation working...
Done!
Unzipping reuslts into arrays...
TIME ANALYSIS:
Generation of Particles     0:00:00.187369
Simulation of Photon Travel 0:00:12.363645
Total Time Elapsed:         0:00:12.584832
RESULTS SUMMARY:
HITS on T1 66
RATIO T1   total photons 16226.0 total incident photons 66 ratio=245.85
HITS on T4 3
RATIO T4   total photons  33172.0 total incident photons 3 ratio=11057.33
DISTANCE: 
```

### Run2: Moved run() out of module file 
- run() out of the module file and running it on `manTof.py`
- added print statements 
```
######################################################
Generated Apparatus Simulation with following defaults
######################################################
PARTICLE: Mean Free Path = 0.00024 cm
PARTICLE: Time Seperation between sequential Particles if simulation more than 1 = 100000.0
SCINT:    Probability of Scintillaton = 0.8
PMT:      Quantum Efficiency is set to 1 by default to keep more pulses
PMT:      Energy per Photoelectron is set to 20 by best estimation
PMT:      Artificial Gain on Output Current = 1
OUTPUT:   Binning Width for PWL output file = 100 ps

Run with .run() function given optional arguments below
integer n particles, 'delta_t' = 100000.0 ps particle time seperation
Generating 1 particles now...
SIZE OF RES:  64
Photons generated 50044.0
SIZE OF TIMES, POINTS, ETC:  50272
Time to unpack times,points, etc. 0.0043006999912904575 <------------
Photons in T1: 16557.0 and Photons in T4: 33487.0
T1 Photon Propagation working...
Time to process T1: 4.643296800000826  <------------
Done!
T4 Photon Propagation working...
Time to process T4: 12.528069000007235  <------------
Done!
Unzipping reuslts into arrays...
Time to unpack T1: 0.003748599992832169  <------------
Time to unpack T4: 0.008279799993033521  <------------
TIME ANALYSIS:
Generation of Particles     0:00:00.204723
Simulation of Photon Travel 0:00:17.216911
Total Time Elapsed:         0:00:17.465147
RESULTS SUMMARY:
HITS on T1 13
RATIO T1   total photons 16557.0 total incident photons 13 ratio=1273.62
HITS on T4 7
RATIO T4   total photons  33487.0 total incident photons 7 ratio=4783.86
DISTANCE: 
Exporing to 2 channels...
Smoothing Signals...
Formatting PWL dataframe...
Smoothing Signals...
Formatting PWL dataframe...
Done!
```
Ran this 10 times analyze 
```
Time to unpack times,points, etc. 0.0038453000015579164
Time to process T1: 4.717084000003524
Time to process T4: 13.080929999996442
Time to unpack T1: 0.0016912999999476597
Time to unpack T4: 0.004023399989819154
Generation of Particles     0:00:00.207846
Simulation of Photon Travel 0:00:17.828522
Total Time Elapsed:         0:00:18.078770

Time to unpack times,points, etc. 0.005897300012293272
Time to process T1: 5.369369199994253
Time to process T4: 12.483720100004575
Time to unpack T1: 0.0017503999988548458
Time to unpack T4: 0.003987299991422333
Generation of Particles     0:00:00.216202
Simulation of Photon Travel 0:00:17.905974
Total Time Elapsed:         0:00:18.172509

Time to unpack times,points, etc. 0.004000300003099255
Time to process T1: 6.334556400004658
Time to process T4: 15.196734099998139
Time to unpack T1: 0.0017188000056194142
Time to unpack T4: 0.00336990000505466
Generation of Particles     0:00:00.216428
Simulation of Photon Travel 0:00:21.561585
Total Time Elapsed:         0:00:21.817211

Time to unpack times,points, etc. 0.004335399993578903
Time to process T1: 7.510844600008568
Time to process T4: 14.187703799994779
Time to unpack T1: 0.0019486000092001632
Time to unpack T4: 0.0037727999879280105
Generation of Particles     0:00:00.207782
Simulation of Photon Travel 0:00:21.731866
Total Time Elapsed:         0:00:21.981767

Time to unpack times,points, etc. 0.005189599993173033
Time to process T1: 8.015519499997026
Time to process T4: 17.222090399998706
Time to unpack T1: 0.0022756999969715253
Time to unpack T4: 0.004099199999473058
Generation of Particles     0:00:00.321412
Simulation of Photon Travel 0:00:25.267383
Total Time Elapsed:         0:00:25.626843

Time to unpack times,points, etc. 0.005598600007942878
Time to process T1: 4.67165709999972
Time to process T4: 12.166835600000923
Time to unpack T1: 0.0015801999979885295
Time to unpack T4: 0.00347729999339208
Generation of Particles     0:00:00.194377
Simulation of Photon Travel 0:00:16.868972
Total Time Elapsed:         0:00:17.098705

Time to unpack times,points, etc. 0.0041745000053197145
Time to process T1: 4.864355199999409
Time to process T4: 11.663345499997376
Time to unpack T1: 0.002212599996710196
Time to unpack T4: 0.00328990000707563
Generation of Particles     0:00:00.163343
Simulation of Photon Travel 0:00:16.557039
Total Time Elapsed:         0:00:16.759118

Time to unpack times,points, etc. 0.0034118000039597973
Time to process T1: 4.2409252999932505
Time to process T4: 10.348038600001018
Time to unpack T1: 0.00176090000604745
Time to unpack T4: 0.003521299993735738
Generation of Particles     0:00:00.184783
Simulation of Photon Travel 0:00:14.617020
Total Time Elapsed:         0:00:14.836530

Time to unpack times,points, etc. 0.003861199991661124
Time to process T1: 4.486565100000007
Time to process T4: 10.617457700005616
Time to unpack T1: 0.0022634000051766634
Time to unpack T4: 0.005091400002129376
Generation of Particles     0:00:00.181520
Simulation of Photon Travel 0:00:15.139724
Total Time Elapsed:         0:00:15.360047

Time to unpack times,points, etc. 0.003699699998833239
Photons in T1: 17355.0 and Photons in T4: 34868.0
Time to process T1: 4.750455300003523
Time to process T4: 11.839638699995703
Time to unpack T1: 0.0016458000027341768
Time to unpack T4: 0.0033072999940486625
Generation of Particles     0:00:00.174950
Simulation of Photon Travel 0:00:16.616631
Total Time Elapsed:         0:00:16.825960

Time to unpack times,points, etc. 0.0037311000050976872
Time to process T1: 5.100151799997548
Time to process T4: 13.393178099999204
Time to unpack T1: 0.001559099997393787
Time to unpack T4: 0.003504599997540936
Generation of Particles     0:00:00.186976
Simulation of Photon Travel 0:00:18.532617
Total Time Elapsed:         0:00:18.761063

Time to unpack times,points, etc. 0.0035319000016897917
Time to process T1: 4.485526199991
Time to process T4: 10.23125379999692
Time to unpack T1: 0.0017469000013079494
Time to unpack T4: 0.0035051000013481826
Generation of Particles     0:00:00.195727
Simulation of Photon Travel 0:00:14.746720
Total Time Elapsed:         0:00:14.980612

Time to unpack times,points, etc. 0.00518190000730101
Time to process T1: 4.084145600005286
Time to process T4: 12.34758060000604
Time to unpack T1: 0.00147790000482928
Time to unpack T4: 0.0030190999968908727
Generation of Particles     0:00:00.216756
Simulation of Photon Travel 0:00:16.459638
Total Time Elapsed:         0:00:16.713639

Time to unpack times,points, etc. 0.004287999996449798
Time to process T1: 4.322206700002425
Time to process T4: 13.715005199992447
Time to unpack T1: 0.002017199993133545
Time to unpack T4: 0.0038554999919142574
Generation of Particles     0:00:00.200458
Simulation of Photon Travel 0:00:18.068496
Total Time Elapsed:         0:00:18.310706

Time to unpack times,points, etc. 0.003803699990385212
Time to process T1: 4.74256250000326
Time to process T4: 13.440829300001496
Time to unpack T1: 0.0018959999870276079
Time to unpack T4: 0.003311099993879907
Generation of Particles     0:00:00.245326
Simulation of Photon Travel 0:00:18.213929
Total Time Elapsed:         0:00:18.494478
```

### Run3: run() in `mainTof.py` added `.close()` and `.join()`
- no other changes 

```
######################################################
Generated Apparatus Simulation with following defaults
######################################################
PARTICLE: Mean Free Path = 0.00024 cm
PARTICLE: Time Seperation between sequential Particles if simulation more than 1 = 100000.0
SCINT:    Probability of Scintillaton = 0.8
PMT:      Quantum Efficiency is set to 1 by default to keep more pulses
PMT:      Energy per Photoelectron is set to 20 by best estimation
PMT:      Artificial Gain on Output Current = 1
OUTPUT:   Binning Width for PWL output file = 100 ps

Run with .run() function given optional arguments below
integer n particles, 'delta_t' = 100000.0 ps particle time seperation
Generating 1 particles now...
SIZE OF RES:  64
Photons generated 50303.0
SIZE OF TIMES, POINTS, ETC:  50168
Time to unpack times,points, etc. 0.002401699995971285 <------------
Photons in T1: 16492.0 and Photons in T4: 33811.0
T1 Photon Propagation working...
Time to process T1: 4.139595000000554 <------------
Done!
T4 Photon Propagation working...
Time to process T4: 11.396087999994052 <------------
Done!
Unzipping reuslts into arrays...
Time to unpack T1: 0.0015951999957906082 <------------
Time to unpack T4: 0.00310689999605529 <------------
TIME ANALYSIS:
Generation of Particles     0:00:00.268368
Simulation of Photon Travel 0:00:15.559308
Total Time Elapsed:         0:00:15.860324
RESULTS SUMMARY:
HITS on T1 15
RATIO T1   total photons 16492.0 total incident photons 15 ratio=1099.47
HITS on T4 8
RATIO T4   total photons  33811.0 total incident photons 8 ratio=4226.38
DISTANCE: 
Exporing to 2 channels...
Smoothing Signals...
Formatting PWL dataframe...
Smoothing Signals...
Formatting PWL dataframe...
Done!
```

15 samples .. I was watching a youtube video... and the times were affected 
```
Time to unpack times,points, etc. 0.002509400001144968
Time to process T1: 5.02223000000231
Time to process T4: 9.664573399990331
Time to unpack T1: 0.001724300003843382
Time to unpack T4: 0.003609500010497868
Generation of Particles     0:00:00.295424
Simulation of Photon Travel 0:00:14.719148
Total Time Elapsed:         0:00:15.055746

Time to unpack times,points, etc. 0.0032100999960675836
Time to process T1: 5.241305999996257
Time to process T4: 15.809703200007789
Time to unpack T1: 0.0018291999876964837
Time to unpack T4: 0.0049673000030452386
Generation of Particles     0:00:00.225565
Simulation of Photon Travel 0:00:21.084554
Total Time Elapsed:         0:00:21.345556

Time to unpack times,points, etc. 0.0025191000022459775
Time to process T1: 6.416496500009089
Time to process T4: 12.862061399995582
Time to unpack T1: 0.003260099998442456
Time to unpack T4: 0.005013599991798401
Generation of Particles     0:00:00.239213
Simulation of Photon Travel 0:00:19.311876
Total Time Elapsed:         0:00:19.595070

Time to unpack times,points, etc. 0.002636700010043569
Time to process T1: 6.341107499989448
Time to process T4: 13.13832719999482
Time to unpack T1: 0.001910100007080473
Time to unpack T4: 0.004655299999285489
Generation of Particles     0:00:00.242716
Simulation of Photon Travel 0:00:19.520187
Total Time Elapsed:         0:00:19.817901

Time to unpack times,points, etc. 0.0031235000060405582
Time to process T1: 7.6858969000022626
Time to process T4: 14.937267899993458
Time to unpack T1: 0.0026584000006550923
Time to unpack T4: 0.0038923999964026734
Generation of Particles     0:00:00.192505
Simulation of Photon Travel 0:00:22.655905
Total Time Elapsed:         0:00:22.898012

Time to unpack times,points, etc. 0.0037701999972341582
Time to process T1: 5.952577399992151
Time to process T4: 15.900621200009482
Time to unpack T1: 0.0015655000024707988
Time to unpack T4: 0.0030743000097572803
Generation of Particles     0:00:00.271378
Simulation of Photon Travel 0:00:21.891316
Total Time Elapsed:         0:00:22.219520

Time to unpack times,points, etc. 0.0021852000063518062
Time to process T1: 7.045508799987147
Time to process T4: 14.018052900006296
Time to unpack T1: 0.0017857999919215217
Time to unpack T4: 0.0036008000024594367
Generation of Particles     0:00:00.200034
Simulation of Photon Travel 0:00:21.092211
Total Time Elapsed:         0:00:21.330255

Time to unpack times,points, etc. 0.0025425999920116737
Time to process T1: 7.542093600000953
Time to process T4: 16.048452699993504
Time to unpack T1: 0.0028106999961892143
Time to unpack T4: 0.0035709999938262627
Generation of Particles     0:00:00.237167
Simulation of Photon Travel 0:00:23.628461
Total Time Elapsed:         0:00:23.902074

Time to unpack times,points, etc. 0.0025880000030156225
Time to process T1: 7.641546899991226
Time to process T4: 15.177921700000297
Time to unpack T1: 0.001638799993088469
Time to unpack T4: 0.0044071000011172146
Generation of Particles     0:00:00.213527
Simulation of Photon Travel 0:00:22.856056
Total Time Elapsed:         0:00:23.110998

Time to unpack times,points, etc. 0.0030282999941846356
Time to process T1: 8.571980499997153
Time to process T4: 15.865097300003981
Time to unpack T1: 0.0022413999977288768
Time to unpack T4: 0.003618100003222935
Generation of Particles     0:00:00.213530
Simulation of Photon Travel 0:00:24.464594
Total Time Elapsed:         0:00:24.721597

Time to unpack times,points, etc. 0.002366999993682839
Time to process T1: 7.022463300003437
Time to process T4: 15.416343199991388
Time to unpack T1: 0.002511899991077371
Time to unpack T4: 0.004090000002179295
Generation of Particles     0:00:00.224899
Simulation of Photon Travel 0:00:22.473050
Total Time Elapsed:         0:00:22.738384

Time to unpack times,points, etc. 0.002204500007792376
Time to process T1: 3.8341020999941975
Time to process T4: 9.162641999995685
Time to unpack T1: 0.001779099999112077
Time to unpack T4: 0.003414499995415099
Generation of Particles     0:00:00.164971
Simulation of Photon Travel 0:00:13.022813
Total Time Elapsed:         0:00:13.220030

```

## Memory 
### Run1: Not using `.close` and `.join`

```
Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    26  89.0586 MiB  89.0586 MiB           1   @profile(precision=4)
    27                                         def run(simClass, *arg, **kwargs):
    28  89.0586 MiB   0.0000 MiB           1       import gc
    29  89.0586 MiB   0.0000 MiB           1       freeze_support() # best practice 
    30  89.0586 MiB   0.0000 MiB           1       if arg:
    31  89.0586 MiB   0.0000 MiB           1           simClass.num_particles = int(arg[0])
    32  89.0586 MiB   0.0000 MiB           1           print(f"Generating {simClass.num_particles} particles now...")
    33                                             else:
    34                                                 simClass.num_particles = 1
    35                                                 print(f"Generating {simClass.num_particles} particle now...")
    36  89.0586 MiB   0.0000 MiB           1       simClass.seperation_time = kwargs.get('delta_t', simClass.seperation_time) # in ps
    37  89.0586 MiB   0.0000 MiB           1       logstarttime = perf_counter()
    38                                             # FIND PARTICLE PATH
    39  89.0586 MiB   0.0000 MiB           1       times = []
    40  89.0586 MiB   0.0000 MiB           1       points = []
    41  89.0586 MiB   0.0000 MiB           1       photons = []
    42  89.0586 MiB   0.0000 MiB           1       particleID = []
    43  89.0586 MiB   0.0000 MiB           1       i = 0
    44  87.0820 MiB  -1.9766 MiB           1       with Pool(processes=cpu_count()-1) as pool:
    45  87.5742 MiB   0.4922 MiB           1           res = pool.map(simClass.particle_task, range(simClass.num_particles))
    46                                                 # pool.close()
    47                                                 # pool.join()
    48  87.5742 MiB   0.0000 MiB           1           import sys
    49  87.5742 MiB   0.0000 MiB           1           print(f'SIZE OF RES: ', sys.getsizeof(res))
    50  87.5742 MiB   0.0000 MiB           1           start = perf_counter()
    51  88.8320 MiB   0.0000 MiB           2           for (time_i, point_i, photon_i) in res:
    52  87.5742 MiB   0.0000 MiB           1               i = 0
    53  87.5977 MiB   0.0234 MiB           1               times.extend(time_i)
    54  88.4453 MiB   0.8477 MiB           1               points.extend(point_i)
    55  88.6836 MiB   0.2383 MiB           1               photons.extend(photon_i)
    56  88.8320 MiB   0.1484 MiB           1               particleID.extend(np.repeat(i, len(time_i))) # particle it belongs to
    57  88.8320 MiB   0.0000 MiB           1               i += 1
    58  89.3906 MiB   0.5586 MiB           1           end = perf_counter() -start
    59  89.3906 MiB   0.0000 MiB           1       logendparticle = perf_counter()
    60  89.3906 MiB   0.0000 MiB           1       N = np.sum(photons)
    61  89.3906 MiB   0.0000 MiB           1       print("Photons generated", N)
    62  89.5117 MiB   0.1211 MiB           1       times = np.asarray(times); points = np.asarray(points); photons = np.asarray(photons); particleID = np.asarray(particleID)
    63  89.5117 MiB   0.0000 MiB           1       print(f'SIZE OF TIMES, POINTS, ETC: ', sys.getsizeof(times))
    64  89.5117 MiB   0.0000 MiB           1       print('Time to unpack times,points, etc.', end)
    65                                             # RETURNS A FILE
    66                                             # SPLIT HERE
    67                                             # RUN #2
    68                                             
    69                                         
    70                                             # SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
    71                                             # Gather TOF data
    72  89.5117 MiB   0.0000 MiB           1       T1_input_times = []
    73  89.5117 MiB   0.0000 MiB           1       T4_input_times = []
    74                                             # Gather Extra Data for analysis
    75  89.5117 MiB   0.0000 MiB           1       simClass.T1_prop_dist = []
    76  89.5117 MiB   0.0000 MiB           1       simClass.T4_prop_dist = []
    77  89.5117 MiB   0.0000 MiB           1       simClass.T1_endpoint_dist = []
    78  89.5117 MiB   0.0000 MiB           1       simClass.T4_endpoint_dist = []
    79  89.5117 MiB   0.0000 MiB           1       simClass.T1_prop_times = []
    80  89.5117 MiB   0.0000 MiB           1       simClass.T4_prop_times = []
    81  89.5117 MiB   0.0000 MiB           1       simClass.T1_interactions = []
    82  89.5117 MiB   0.0000 MiB           1       simClass.T4_interactions = []
    83  89.5117 MiB   0.0000 MiB           1       simClass.T1_part_ids = []
    84  89.5117 MiB   0.0000 MiB           1       simClass.T4_part_ids = []
    85  89.5117 MiB   0.0000 MiB           1       T1points = (points[:])[points[:,2] >= simClass.T1z]
    86  89.5117 MiB   0.0000 MiB           1       T1times = (times[:])[points[:,2] >= simClass.T1z]
    87  89.5117 MiB   0.0000 MiB           1       T1photons = (photons[:])[points[:,2] >= simClass.T1z]
    88  89.5117 MiB   0.0000 MiB           1       T1part_ids = (particleID[:])[points[:,2] >= simClass.T1z]
    89  89.5117 MiB   0.0000 MiB           1       T1part_ids = np.repeat(T1part_ids, T1photons.astype(int), axis=0) # big id bank
    90  89.5117 MiB   0.0000 MiB           1       T4points = (points[:])[points[:,2] < simClass.T1z]
    91  89.5117 MiB   0.0000 MiB           1       T4times = (times[:])[points[:,2] < simClass.T1z]
    92  89.5117 MiB   0.0000 MiB           1       T4photons = (photons[:])[points[:,2] < simClass.T1z]
    93  89.5117 MiB   0.0000 MiB           1       T4part_ids = (particleID[:])[points[:,2] < simClass.T1z]
    94  89.9492 MiB   0.4375 MiB           1       T4part_ids = np.repeat(T4part_ids, T4photons.astype(int), axis=0) # big id bank
    95  89.9492 MiB   0.0000 MiB           1       print(f"Photons in T1: {np.sum(T1photons)} and Photons in T4: {np.sum(T4photons)}")
    96  89.9492 MiB   0.0000 MiB           1       del times; del points; del photons; # remove copies
    97  89.9688 MiB   0.0195 MiB           1       gc.collect()
    98  89.9688 MiB   0.0000 MiB           1       logstartphoton = perf_counter()
    99                                         
   100                                             # check this link https://stackoverflow.com/questions/14749897/python-multiprocessing-memory-usage
   101  89.9531 MiB  -0.0156 MiB           1       with Pool(processes=cpu_count()) as pool: # this way of making the pool causes all the data to copy! 
   102  89.9531 MiB   0.0000 MiB           1           print("T1 Photon Propagation working...")
   103  89.9531 MiB   0.0000 MiB           1           start = perf_counter()
   104  93.4219 MiB   3.4688 MiB           1           T1res = pool.starmap(simClass.scint_taskT1, np.repeat(np.c_[T1points,T1times],T1photons.astype(int), axis=0))
   105  93.4219 MiB   0.0000 MiB           1           end = perf_counter() - start
   106  93.4219 MiB   0.0000 MiB           1           print("Time to process T1:", end)
   107  93.4219 MiB   0.0000 MiB           1           print("Done!")
   108                                         
   109                                         
   110  93.4219 MiB   0.0000 MiB           1           print("T4 Photon Propagation working...")
   111  93.4219 MiB   0.0000 MiB           1           start = perf_counter()
   112 102.7500 MiB   9.3281 MiB           1           T4res = pool.starmap(simClass.scint_taskT4, np.repeat(np.c_[T4points,T4times],T4photons.astype(int), axis=0))
   113 102.7500 MiB   0.0000 MiB           1           end = perf_counter() - start
   114 102.7500 MiB   0.0000 MiB           1           print("Time to process T4:", end)
   115 102.7500 MiB   0.0000 MiB           1           print("Done!")
   116 102.7500 MiB   0.0000 MiB           1           print("Unzipping reuslts into arrays...")
   117                                         
   118                                                 # pool.close()
   119                                                 # pool.join()
   120                                         
   121 102.7500 MiB   0.0000 MiB           1           start = perf_counter()
   122 102.7500 MiB   0.0000 MiB       16708           for (T1hit_PMT, T1travel_time, T1tot_dist, T1endpt, T1bounces, T1prop),T1part_id in zip(T1res, T1part_ids): # check if moving starmap here helps
   123 102.7500 MiB   0.0000 MiB       16707               if T1hit_PMT:
   124 102.7500 MiB   0.0000 MiB          11                   T1_input_times.append(T1travel_time)
   125 102.7500 MiB   0.0000 MiB          11                   simClass.T1_prop_dist.append(T1tot_dist)
   126 102.7500 MiB   0.0000 MiB          11                   simClass.T1_endpoint_dist.append(T1endpt)
   127 102.7500 MiB   0.0000 MiB          11                   simClass.T1_prop_times.append(T1prop)
   128 102.7500 MiB   0.0000 MiB          11                   simClass.T1_interactions.append(T1bounces)
   129 102.7500 MiB   0.0000 MiB          11                   simClass.T1_part_ids.append(T1part_id)
   130 102.7500 MiB   0.0000 MiB           1           end = perf_counter() - start
   131 102.7500 MiB   0.0000 MiB           1           print("Time to unpack T1:", end)
   132                                         
   133 102.7500 MiB   0.0000 MiB           1           start = perf_counter()
   134 102.7500 MiB   0.0000 MiB       33861           for (T4hit_PMT, T4travel_time, T4tot_dist, T4endpt, T4bounces, T4prop),T4part_id in zip(T4res, T4part_ids): # check if moving starmap here helps
   135 102.7500 MiB   0.0000 MiB       33860               if T4hit_PMT:
   136 102.7500 MiB   0.0000 MiB           2                   T4_input_times.append(T4travel_time)
   137 102.7500 MiB   0.0000 MiB           2                   simClass.T4_prop_dist.append(T4tot_dist)
   138 102.7500 MiB   0.0000 MiB           2                   simClass.T4_endpoint_dist.append(T4endpt)
   139 102.7500 MiB   0.0000 MiB           2                   simClass.T4_prop_times.append(T4prop)
   140 102.7500 MiB   0.0000 MiB           2                   simClass.T4_interactions.append(T4bounces)
   141 102.7500 MiB   0.0000 MiB           2                   simClass.T4_part_ids.append(T4part_id)
   142 102.7500 MiB   0.0000 MiB           1           end = perf_counter() - start
   143 102.9453 MiB   0.1953 MiB           1           print("Time to unpack T4:", end)
   144 102.9453 MiB   0.0000 MiB           1       logendtime = perf_counter()
   145                                             # PRINT RESULTS
   146 102.9453 MiB   0.0000 MiB           1       print("TIME ANALYSIS:")
   147 102.9453 MiB   0.0000 MiB           1       pgtime = timedelta(seconds=logendparticle-logstarttime)
   148 102.9453 MiB   0.0000 MiB           1       phtime = timedelta(seconds=logendtime-logstartphoton)
   149 102.9453 MiB   0.0000 MiB           1       ttime = timedelta(seconds=logendtime-logstarttime)
   150 102.9453 MiB   0.0000 MiB           1       print(f"Generation of Particles     {str(pgtime)}")
   151 102.9453 MiB   0.0000 MiB           1       print(f"Simulation of Photon Travel {str(phtime)}")
   152 102.9453 MiB   0.0000 MiB           1       print(f"Total Time Elapsed:         {str(ttime)}")
   153 102.9453 MiB   0.0000 MiB           1       print("RESULTS SUMMARY:")
   154 102.9453 MiB   0.0000 MiB           1       print("HITS on T1",len(T1_input_times))
   155 103.0195 MiB   0.0742 MiB           1       print("RATIO T1   total photons", np.sum(T1photons), "total incident photons", len(T1_input_times), f"ratio={np.sum(T1photons)/len(T1_input_times):.2f}")
   156 103.0195 MiB   0.0000 MiB           1       print("HITS on T4",len(T4_input_times))
   157 103.0195 MiB   0.0000 MiB           1       print("RATIO T4   total photons ", np.sum(T4photons),"total incident photons", len(T4_input_times), f"ratio={np.sum(T4photons)/len(T4_input_times):.2f}")
   158 103.0195 MiB   0.0000 MiB           1       print("DISTANCE: ")
   159 103.0195 MiB   0.0000 MiB           1       del T1points; del T1times; del T1photons; del T4points; del T4times; del T4photons; # remove unused variables
   160 103.0195 MiB   0.0000 MiB           1       gc.collect()
   161                                             # print(T4_input_times)
   162                                             # BEGIN SIMULATING PMT PULSE
   163 103.0195 MiB   0.0000 MiB           1       signals_channelT1 = []
   164 103.0195 MiB   0.0000 MiB           1       signals_channelT4 = []
   165 103.0195 MiB   0.0000 MiB           1       output_times_channelT1 = []
   166 103.0195 MiB   0.0000 MiB           1       output_times_channelT4 = []
   167 103.0195 MiB   0.0000 MiB           1       signals = []
   168 103.0195 MiB   0.0000 MiB           1       output_times = []
   169 103.0195 MiB   0.0000 MiB          12       for t in T1_input_times:
   170 103.0195 MiB   0.0000 MiB          11           pmtSignal_i = simClass.photontoElectrons(1)
   171 103.0195 MiB   0.0000 MiB          11           output_times.append(simClass.pmt_electron_travel_time+t)
   172 103.0195 MiB   0.0000 MiB          11           output_times_channelT1.append(simClass.pmt_electron_travel_time+t)
   173 103.0195 MiB   0.0000 MiB          11           signals.append(pmtSignal_i)
   174 103.0195 MiB   0.0000 MiB          11           signals_channelT1.append(pmtSignal_i)
   175 103.0195 MiB   0.0000 MiB           3       for t in T4_input_times:
   176 103.0195 MiB   0.0000 MiB           2           pmtSignal_i = simClass.photontoElectrons(1)
   177 103.0195 MiB   0.0000 MiB           2           output_times.append(simClass.pmt_electron_travel_time+t)
   178 103.0195 MiB   0.0000 MiB           2           output_times_channelT4.append(simClass.pmt_electron_travel_time+t)
   179 103.0195 MiB   0.0000 MiB           2           signals.append(pmtSignal_i)
   180 103.0195 MiB   0.0000 MiB           2           signals_channelT4.append(pmtSignal_i)
   181                                         
   182                                             # CONVERTION Electron count to Current and save in array
   183 103.0195 MiB   0.0000 MiB           1       simClass.signals = np.array(signals) * simClass.q / 1e-12 * simClass.artificial_gain # divided by 1ps 
   184 103.0195 MiB   0.0000 MiB           1       simClass.output_times = np.array(output_times)
   185 103.0195 MiB   0.0000 MiB           1       simClass.signals_channelT1 = np.array(signals_channelT1) * simClass.q / 1e-12 * simClass.artificial_gain
   186 103.0195 MiB   0.0000 MiB           1       simClass.signals_channelT4 = np.array(signals_channelT4) * simClass.q / 1e-12 * simClass.artificial_gain * 0.6 # factor to limit pulses to 50miliamps and stop contant comparator firing. however, current should be smaller from Quantum Efficiency and current should be larger from 3kV potential difference across PMT dynodes instead of current 1kV potential difference
   187 103.0195 MiB   0.0000 MiB           1       simClass.output_times_channelT1 = np.array(output_times_channelT1)
   188 103.0195 MiB   0.0000 MiB           1       simClass.output_times_channelT4 = np.array(output_times_channelT4)
   189                                         
   190                                         # Output function
   191                                         # def to_csv(simClass, **kwargs):
   192 139.2773 MiB  36.2578 MiB           1       from scipy.stats import norm # type: ignore 
   193 139.2773 MiB   0.0000 MiB           1       output_extra = kwargs.get('extra_data_only', False)
   194 139.2773 MiB   0.0000 MiB           1       output_both = kwargs.get('output_both', False)
   195                                             # OUTPUT FORMATTING
   196 139.2773 MiB   0.0000 MiB           1       if output_extra or output_both:
   197                                                 print("Exporting Extra Data...")
   198                                                 dft1 = pd.DataFrame({'T1_part_ids':simClass.T1_part_ids,'time':simClass.output_times_channelT1,'T1_prop_dist':simClass.T1_prop_dist,'T1_endpoint_dist':simClass.T1_endpoint_dist, 'T1_prop_times':simClass.T1_prop_times, 'T1_interactions':simClass.T1_interactions})
   199                                                 dft4 = pd.DataFrame({'T4_part_ids':simClass.T4_part_ids,'time':simClass.output_times_channelT4,'T4_prop_dist':simClass.T4_prop_dist,'T4_endpoint_dist':simClass.T4_endpoint_dist, 'T4_prop_times':simClass.T4_prop_times, 'T4_interactions':simClass.T4_interactions})
   200                                                 dft1.to_csv('monte_carlo_extradata'+str(simClass.num_particles)+'chT1_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
   201                                                 dft4.to_csv('monte_carlo_extradata'+str(simClass.num_particles)+'chT4_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
   202                                                 if not output_both:
   203                                                     return
   204 139.2773 MiB   0.0000 MiB           1       print("Exporing to 2 channels...")
   205                                             # for each channel
   206 140.9414 MiB   0.0000 MiB           3       for time,signal,ch in zip([simClass.output_times_channelT1,simClass.output_times_channelT4],[simClass.signals_channelT1,simClass.signals_channelT4],[1,4]):
   207                                         
   208                                                 # from io import StringIO
   209                                                 # from csv import writer 
   210                                                 # output = StringIO()
   211                                                 # csv_writer = writer(output)
   212                                                 
   213 140.9414 MiB   0.0000 MiB           2           print("Smoothing Signals...")
   214 140.9414 MiB   0.0000 MiB           2           t_binned = [0.] # left edges of bins
   215 140.9414 MiB   0.0000 MiB           2           y_binned = [0.]
   216 140.9414 MiB   0.0000 MiB          15           for i,y in enumerate(signal):
   217                                                     # print(f"i={i},t[{i}]={time[i]} y[{i}]={y}")
   218 140.9414 MiB   0.0000 MiB          13               lower_bound = max(time[i]-2*simClass.sigma_smoothing,0) # 2 sigma away backward
   219 140.9414 MiB   0.0000 MiB          13               upper_bound = min(time[i]+2*simClass.sigma_smoothing,max(time)+2*simClass.sigma_smoothing) # 2 sigma away forward
   220                                                     # MAKE NEW DATA CENTERED AROUND PULSE
   221 140.9414 MiB   0.0000 MiB          13               if lower_bound < max(t_binned): # if already binned
   222 140.9414 MiB   0.0000 MiB          11                   lower_bound = t_binned[np.digitize(lower_bound, t_binned)]+simClass.output_bin_width/2
   223 140.9414 MiB   0.0000 MiB          13               cur_x = np.arange(lower_bound,upper_bound,simClass.output_bin_width)+simClass.output_bin_width/2
   224                                                     # print(f"cur_x from {lower_bound}-->{upper_bound}", cur_x)
   225                                                     # ADD DATA IF NEEDED
   226 140.9414 MiB   0.0000 MiB         213               for x in cur_x:
   227 140.9414 MiB   0.0000 MiB         200                   if x > max(t_binned): 
   228 140.9414 MiB   0.0000 MiB          35                       t_binned.append(x)
   229 140.9414 MiB   0.0000 MiB          35                       y_binned.append(0)
   230 140.9414 MiB   0.2969 MiB         165                   elif (np.digitize(x, t_binned)-1 > 0) and (np.digitize(x, t_binned) < len(t_binned)):
   231 140.9414 MiB   0.0000 MiB         163                       index = np.digitize(x, t_binned)
   232 140.9414 MiB   0.0000 MiB         163                       if abs(t_binned[index]-t_binned[index-1]) > simClass.output_bin_width:
   233 139.5742 MiB   0.0000 MiB          10                           t_binned.insert(index, x) # check if need -1 or just np.digitize()
   234 139.5742 MiB   0.0000 MiB          10                           y_binned.insert(index, 0) # check 
   235                                                     # GET INDICIES
   236 140.9414 MiB   0.0000 MiB         336               index_lower = [i for i,t in enumerate(t_binned) if t >= lower_bound][0] # left edge in time binned
   237 140.9414 MiB   0.0000 MiB         336               index_upper = [i for i,t in enumerate(t_binned) if t <= upper_bound][-1] # right edge in time binned
   238                                                     # GAUSSIAN SMOOTH
   239 140.9414 MiB   0.0000 MiB          13               gaussian = norm.pdf(t_binned[index_lower:index_upper], loc=time[i], scale=simClass.sigma_smoothing)*simClass.sigma_smoothing*y/4
   240                                                     # ADD TO CORRECT BINS
   241 140.9414 MiB   0.0000 MiB         248               for i,y_add in enumerate(gaussian):
   242 140.9414 MiB   0.0000 MiB         235                   if y_binned[index_lower+i]+y_add < simClass.max_pmt_current_output:
   243 140.9414 MiB   0.0000 MiB         235                       y_binned[index_lower+i] += y_add
   244                                                         else:
   245                                                             y_binned[index_lower+i] = simClass.max_pmt_current_output
   246                                         
   247 140.9414 MiB   1.0000 MiB           2           df = pd.DataFrame({'time':t_binned,'current':y_binned}).sort_values(by=['time'])
   248 140.9414 MiB   0.0000 MiB           2           print("Formatting PWL dataframe...")
   249 140.9414 MiB   0.0000 MiB           2           fill_data = []                                                                      # declare empty array
   250                                                 # begin padding data at time 1/5th bin width before first time stamp
   251 140.9414 MiB   0.0000 MiB           2           fill_data.append([df['time'].iloc[0]-simClass.output_bin_width/5,0])                    # add zero at beginning
   252 140.9414 MiB   0.0000 MiB          47           for i in range(len(df['time'])-1):                                                        # for each time index
   253 140.9414 MiB   0.0000 MiB          45               if abs(df['time'].iloc[i]-df['time'].iloc[i+1]) > simClass.output_bin_width:        # if dt between signals is greater than minimum bin width
   254 140.9414 MiB   0.0000 MiB           3                   fill_data.append([df['time'].iloc[i]+simClass.output_bin_width/5,0])            # place zero after current signal
   255 140.9414 MiB   0.0000 MiB           3                   fill_data.append([df['time'].iloc[i+1]-simClass.output_bin_width/5,0])          # place zero before next signal
   256 140.9414 MiB   0.0000 MiB           2           fill_data.append([df['time'].iloc[-1]+simClass.output_bin_width/5,0])                   # add zero at end
   257 140.9414 MiB   0.0000 MiB           2           fill_data = np.array(fill_data)
   258 140.9414 MiB   0.0000 MiB           2           fill = pd.DataFrame(fill_data, columns=['time','current'])
   259 140.9414 MiB   0.1719 MiB           2           df = pd.concat([fill, df], ignore_index=True).sort_values(by=['time']).reset_index(drop=True)
   260 140.9414 MiB   0.0977 MiB           2           df['time'] = df['time']/1e12
   261 140.9414 MiB   0.0000 MiB           2           df = df[['time', 'current']] # need this for LTSpice PWL current input file to work
   262 140.9414 MiB   0.0977 MiB           2           df.to_csv('monte_carlo_input'+str(simClass.num_particles)+'ch'+str(ch)+'_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt', float_format='%.13f', header=False, index=False, sep=' ') # PWL file formatting
   263 140.9414 MiB   0.0000 MiB           1       print("Done!")

```

### Run2: with `.close` and `.join`

```
Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    26     91.3 MiB     91.3 MiB           1   @profile
    27                                         # @profile(precision=4)
    28                                         def run(simClass, *arg, **kwargs):
    29     91.3 MiB      0.0 MiB           1       import gc
    30     91.3 MiB      0.0 MiB           1       freeze_support() # best practice 
    31     91.3 MiB      0.0 MiB           1       if arg:
    32     91.3 MiB      0.0 MiB           1           simClass.num_particles = int(arg[0])
    33     91.3 MiB      0.0 MiB           1           print(f"Generating {simClass.num_particles} particles now...")
    34                                             else:
    35                                                 simClass.num_particles = 1
    36                                                 print(f"Generating {simClass.num_particles} particle now...")
    37     91.3 MiB      0.0 MiB           1       simClass.seperation_time = kwargs.get('delta_t', simClass.seperation_time) # in ps
    38     91.3 MiB      0.0 MiB           1       logstarttime = perf_counter()
    39                                             # FIND PARTICLE PATH
    40     91.3 MiB      0.0 MiB           1       times = []
    41     91.3 MiB      0.0 MiB           1       points = []
    42     91.3 MiB      0.0 MiB           1       photons = []
    43     91.3 MiB      0.0 MiB           1       particleID = []
    44     91.3 MiB      0.0 MiB           1       i = 0
    45     91.6 MiB      0.3 MiB           1       with Pool(processes=cpu_count()-1) as pool:
    46     92.1 MiB      0.5 MiB           1           res = pool.map(simClass.particle_task, range(simClass.num_particles))
    47     92.1 MiB      0.0 MiB           1           pool.close()
    48     92.5 MiB      0.4 MiB           1           pool.join()
    49     92.5 MiB      0.0 MiB           1           import sys
    50     92.5 MiB      0.0 MiB           1           print(f'SIZE OF RES: ', sys.getsizeof(res))
    51     92.5 MiB      0.0 MiB           1           start = perf_counter()
    52     93.8 MiB      0.0 MiB           2           for (time_i, point_i, photon_i) in res:
    53     92.5 MiB      0.0 MiB           1               i = 0
    54     92.7 MiB      0.2 MiB           1               times.extend(time_i)
    55     93.3 MiB      0.6 MiB           1               points.extend(point_i)
    56     93.4 MiB      0.1 MiB           1               photons.extend(photon_i)
    57     93.8 MiB      0.4 MiB           1               particleID.extend(np.repeat(i, len(time_i))) # particle it belongs to
    58     93.8 MiB      0.0 MiB           1               i += 1
    59     93.8 MiB      0.0 MiB           1           end = perf_counter() -start
    60     93.8 MiB      0.0 MiB           1       logendparticle = perf_counter()
    61     93.8 MiB      0.0 MiB           1       N = np.sum(photons)
    62     93.8 MiB      0.0 MiB           1       print("Photons generated", N)
    63     94.1 MiB      0.3 MiB           1       times = np.asarray(times); points = np.asarray(points); photons = np.asarray(photons); particleID = np.asarray(particleID)
    64     94.1 MiB      0.0 MiB           1       print(f'SIZE OF TIMES, POINTS, ETC: ', sys.getsizeof(times))
    65     94.1 MiB      0.0 MiB           1       print('Time to unpack times,points, etc.', end)
    66                                             # RETURNS A FILE
    67                                             # SPLIT HERE
    68                                             # RUN #2
    69                                             
    70                                         
    71                                             # SIMULATE EACH PHOTON PATH IN BOTH SCINTILLATORS
    72                                             # Gather TOF data
    73     94.1 MiB      0.0 MiB           1       T1_input_times = []
    74     94.1 MiB      0.0 MiB           1       T4_input_times = []
    75                                             # Gather Extra Data for analysis
    76     94.1 MiB      0.0 MiB           1       simClass.T1_prop_dist = []
    77     94.1 MiB      0.0 MiB           1       simClass.T4_prop_dist = []
    78     94.1 MiB      0.0 MiB           1       simClass.T1_endpoint_dist = []
    79     94.1 MiB      0.0 MiB           1       simClass.T4_endpoint_dist = []
    80     94.1 MiB      0.0 MiB           1       simClass.T1_prop_times = []
    81     94.1 MiB      0.0 MiB           1       simClass.T4_prop_times = []
    82     94.1 MiB      0.0 MiB           1       simClass.T1_interactions = []
    83     94.1 MiB      0.0 MiB           1       simClass.T4_interactions = []
    84     94.1 MiB      0.0 MiB           1       simClass.T1_part_ids = []
    85     94.1 MiB      0.0 MiB           1       simClass.T4_part_ids = []
    86     94.1 MiB      0.0 MiB           1       T1points = (points[:])[points[:,2] >= simClass.T1z]
    87     94.1 MiB      0.0 MiB           1       T1times = (times[:])[points[:,2] >= simClass.T1z]
    88     94.1 MiB      0.0 MiB           1       T1photons = (photons[:])[points[:,2] >= simClass.T1z]
    89     94.1 MiB      0.0 MiB           1       T1part_ids = (particleID[:])[points[:,2] >= simClass.T1z]
    90     94.1 MiB      0.0 MiB           1       T1part_ids = np.repeat(T1part_ids, T1photons.astype(int), axis=0) # big id bank
    91     94.1 MiB      0.0 MiB           1       T4points = (points[:])[points[:,2] < simClass.T1z]
    92     94.1 MiB      0.0 MiB           1       T4times = (times[:])[points[:,2] < simClass.T1z]
    93     94.1 MiB      0.0 MiB           1       T4photons = (photons[:])[points[:,2] < simClass.T1z]
    94     94.1 MiB      0.0 MiB           1       T4part_ids = (particleID[:])[points[:,2] < simClass.T1z]
    95     94.6 MiB      0.5 MiB           1       T4part_ids = np.repeat(T4part_ids, T4photons.astype(int), axis=0) # big id bank
    96     94.6 MiB      0.0 MiB           1       print(f"Photons in T1: {np.sum(T1photons)} and Photons in T4: {np.sum(T4photons)}")
    97     94.6 MiB      0.0 MiB           1       del times; del points; del photons; # remove copies
    98     94.6 MiB      0.0 MiB           1       gc.collect()
    99     94.6 MiB      0.0 MiB           1       logstartphoton = perf_counter()
   100                                         
   101                                             # check this link https://stackoverflow.com/questions/14749897/python-multiprocessing-memory-usage
   102     94.6 MiB     -0.0 MiB           1       with Pool(processes=cpu_count()) as pool: # this way of making the pool causes all the data to copy! 
   103     94.6 MiB      0.0 MiB           1           print("T1 Photon Propagation working...")
   104     94.6 MiB      0.0 MiB           1           start = perf_counter()
   105     98.1 MiB      3.5 MiB           1           T1res = pool.starmap(simClass.scint_taskT1, np.repeat(np.c_[T1points,T1times],T1photons.astype(int), axis=0))
   106     98.1 MiB      0.0 MiB           1           end = perf_counter() - start
   107     98.1 MiB      0.0 MiB           1           print("Time to process T1:", end)
   108     98.1 MiB      0.0 MiB           1           print("Done!")
   109                                         
   110                                         
   111     98.1 MiB      0.0 MiB           1           print("T4 Photon Propagation working...")
   112     98.1 MiB      0.0 MiB           1           start = perf_counter()
   113    107.2 MiB      9.2 MiB           1           T4res = pool.starmap(simClass.scint_taskT4, np.repeat(np.c_[T4points,T4times],T4photons.astype(int), axis=0))
   114    107.2 MiB      0.0 MiB           1           end = perf_counter() - start
   115    107.2 MiB      0.0 MiB           1           print("Time to process T4:", end)
   116    107.2 MiB      0.0 MiB           1           print("Done!")
   117    107.2 MiB      0.0 MiB           1           print("Unzipping reuslts into arrays...")
   118                                         
   119    107.5 MiB      0.3 MiB           1           pool.close()
   120    107.7 MiB      0.2 MiB           1           pool.join()
   121                                         
   122    107.7 MiB      0.0 MiB           1           start = perf_counter()
   123    107.7 MiB      0.0 MiB       16585           for (T1hit_PMT, T1travel_time, T1tot_dist, T1endpt, T1bounces, T1prop),T1part_id in zip(T1res, T1part_ids): # check if moving starmap here helps
   124    107.7 MiB      0.0 MiB       16584               if T1hit_PMT:
   125    107.7 MiB      0.0 MiB           5                   T1_input_times.append(T1travel_time)
   126    107.7 MiB      0.0 MiB           5                   simClass.T1_prop_dist.append(T1tot_dist)
   127    107.7 MiB      0.0 MiB           5                   simClass.T1_endpoint_dist.append(T1endpt)
   128    107.7 MiB      0.0 MiB           5                   simClass.T1_prop_times.append(T1prop)
   129    107.7 MiB      0.0 MiB           5                   simClass.T1_interactions.append(T1bounces)
   130    107.7 MiB      0.0 MiB           5                   simClass.T1_part_ids.append(T1part_id)
   131    107.7 MiB      0.0 MiB           1           end = perf_counter() - start
   132    107.7 MiB      0.0 MiB           1           print("Time to unpack T1:", end)
   133                                         
   134    107.7 MiB      0.0 MiB           1           start = perf_counter()
   135    107.7 MiB      0.0 MiB       33968           for (T4hit_PMT, T4travel_time, T4tot_dist, T4endpt, T4bounces, T4prop),T4part_id in zip(T4res, T4part_ids): # check if moving starmap here helps
   136    107.7 MiB      0.0 MiB       33967               if T4hit_PMT:
   137    107.7 MiB      0.0 MiB          13                   T4_input_times.append(T4travel_time)
   138    107.7 MiB      0.0 MiB          13                   simClass.T4_prop_dist.append(T4tot_dist)
   139    107.7 MiB      0.0 MiB          13                   simClass.T4_endpoint_dist.append(T4endpt)
   140    107.7 MiB      0.0 MiB          13                   simClass.T4_prop_times.append(T4prop)
   141    107.7 MiB      0.0 MiB          13                   simClass.T4_interactions.append(T4bounces)
   142    107.7 MiB      0.0 MiB          13                   simClass.T4_part_ids.append(T4part_id)
   143    107.7 MiB      0.0 MiB           1           end = perf_counter() - start
   144    107.7 MiB      0.0 MiB           1           print("Time to unpack T4:", end)
   145    107.7 MiB      0.0 MiB           1       logendtime = perf_counter()
   146                                             # PRINT RESULTS
   147    107.7 MiB      0.0 MiB           1       print("TIME ANALYSIS:")
   148    107.7 MiB      0.0 MiB           1       pgtime = timedelta(seconds=logendparticle-logstarttime)
   149    107.7 MiB      0.0 MiB           1       phtime = timedelta(seconds=logendtime-logstartphoton)
   150    107.7 MiB      0.0 MiB           1       ttime = timedelta(seconds=logendtime-logstarttime)
   151    107.7 MiB      0.0 MiB           1       print(f"Generation of Particles     {str(pgtime)}")
   152    107.7 MiB      0.0 MiB           1       print(f"Simulation of Photon Travel {str(phtime)}")
   153    107.7 MiB      0.0 MiB           1       print(f"Total Time Elapsed:         {str(ttime)}")
   154    107.7 MiB      0.0 MiB           1       print("RESULTS SUMMARY:")
   155    107.7 MiB      0.0 MiB           1       print("HITS on T1",len(T1_input_times))
   156    107.7 MiB      0.0 MiB           1       print("RATIO T1   total photons", np.sum(T1photons), "total incident photons", len(T1_input_times), f"ratio={np.sum(T1photons)/len(T1_input_times):.2f}")
   157    107.7 MiB      0.0 MiB           1       print("HITS on T4",len(T4_input_times))
   158    107.7 MiB      0.0 MiB           1       print("RATIO T4   total photons ", np.sum(T4photons),"total incident photons", len(T4_input_times), f"ratio={np.sum(T4photons)/len(T4_input_times):.2f}")
   159    107.7 MiB      0.0 MiB           1       print("DISTANCE: ")
   160    107.7 MiB      0.0 MiB           1       del T1points; del T1times; del T1photons; del T4points; del T4times; del T4photons; # remove unused variables
   161    107.7 MiB      0.0 MiB           1       gc.collect()
   162                                             # print(T4_input_times)
   163                                             # BEGIN SIMULATING PMT PULSE
   164    107.7 MiB      0.0 MiB           1       signals_channelT1 = []
   165    107.7 MiB      0.0 MiB           1       signals_channelT4 = []
   166    107.7 MiB      0.0 MiB           1       output_times_channelT1 = []
   167    107.7 MiB      0.0 MiB           1       output_times_channelT4 = []
   168    107.7 MiB      0.0 MiB           1       signals = []
   169    107.7 MiB      0.0 MiB           1       output_times = []
   170    107.7 MiB      0.0 MiB           6       for t in T1_input_times:
   171    107.7 MiB      0.0 MiB           5           pmtSignal_i = simClass.photontoElectrons(1)
   172    107.7 MiB      0.0 MiB           5           output_times.append(simClass.pmt_electron_travel_time+t)
   173    107.7 MiB      0.0 MiB           5           output_times_channelT1.append(simClass.pmt_electron_travel_time+t)
   174    107.7 MiB      0.0 MiB           5           signals.append(pmtSignal_i)
   175    107.7 MiB      0.0 MiB           5           signals_channelT1.append(pmtSignal_i)
   176    107.7 MiB      0.0 MiB          14       for t in T4_input_times:
   177    107.7 MiB      0.0 MiB          13           pmtSignal_i = simClass.photontoElectrons(1)
   178    107.7 MiB      0.0 MiB          13           output_times.append(simClass.pmt_electron_travel_time+t)
   179    107.7 MiB      0.0 MiB          13           output_times_channelT4.append(simClass.pmt_electron_travel_time+t)
   180    107.7 MiB      0.0 MiB          13           signals.append(pmtSignal_i)
   181    107.7 MiB      0.0 MiB          13           signals_channelT4.append(pmtSignal_i)
   182                                         
   183                                             # CONVERTION Electron count to Current and save in array
   184    107.7 MiB      0.0 MiB           1       simClass.signals = np.array(signals) * simClass.q / 1e-12 * simClass.artificial_gain # divided by 1ps 
   185    107.7 MiB      0.0 MiB           1       simClass.output_times = np.array(output_times)
   186    107.7 MiB      0.0 MiB           1       simClass.signals_channelT1 = np.array(signals_channelT1) * simClass.q / 1e-12 * simClass.artificial_gain
   187    107.7 MiB      0.0 MiB           1       simClass.signals_channelT4 = np.array(signals_channelT4) * simClass.q / 1e-12 * simClass.artificial_gain * 0.6 # factor to limit pulses to 50miliamps and stop contant comparator firing. however, current should be smaller from Quantum Efficiency and current should be larger from 3kV potential difference across PMT dynodes instead of current 1kV potential difference
   188    107.7 MiB      0.0 MiB           1       simClass.output_times_channelT1 = np.array(output_times_channelT1)
   189    107.7 MiB      0.0 MiB           1       simClass.output_times_channelT4 = np.array(output_times_channelT4)
   190                                         
   191                                         # Output function
   192                                         # def to_csv(simClass, **kwargs):
   193    142.3 MiB     34.6 MiB           1       from scipy.stats import norm # type: ignore 
   194    142.3 MiB      0.0 MiB           1       output_extra = kwargs.get('extra_data_only', False)
   195    142.3 MiB      0.0 MiB           1       output_both = kwargs.get('output_both', False)
   196                                             # OUTPUT FORMATTING
   197    142.3 MiB      0.0 MiB           1       if output_extra or output_both:
   198                                                 print("Exporting Extra Data...")
   199                                                 dft1 = pd.DataFrame({'T1_part_ids':simClass.T1_part_ids,'time':simClass.output_times_channelT1,'T1_prop_dist':simClass.T1_prop_dist,'T1_endpoint_dist':simClass.T1_endpoint_dist, 'T1_prop_times':simClass.T1_prop_times, 'T1_interactions':simClass.T1_interactions})
   200                                                 dft4 = pd.DataFrame({'T4_part_ids':simClass.T4_part_ids,'time':simClass.output_times_channelT4,'T4_prop_dist':simClass.T4_prop_dist,'T4_endpoint_dist':simClass.T4_endpoint_dist, 'T4_prop_times':simClass.T4_prop_times, 'T4_interactions':simClass.T4_interactions})
   201                                                 dft1.to_csv('monte_carlo_extradata'+str(simClass.num_particles)+'chT1_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
   202                                                 dft4.to_csv('monte_carlo_extradata'+str(simClass.num_particles)+'chT4_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt') # default sep=','
   203                                                 if not output_both:
   204                                                     return
   205    142.3 MiB      0.0 MiB           1       print("Exporing to 2 channels...")
   206                                             # for each channel
   207    144.2 MiB      0.0 MiB           3       for time,signal,ch in zip([simClass.output_times_channelT1,simClass.output_times_channelT4],[simClass.signals_channelT1,simClass.signals_channelT4],[1,4]):
   208                                         
   209                                                 # from io import StringIO
   210                                                 # from csv import writer 
   211                                                 # output = StringIO()
   212                                                 # csv_writer = writer(output)
   213                                                 
   214    144.2 MiB      0.0 MiB           2           print("Smoothing Signals...")
   215    144.2 MiB      0.0 MiB           2           t_binned = [0.] # left edges of bins
   216    144.2 MiB      0.0 MiB           2           y_binned = [0.]
   217    144.2 MiB      0.0 MiB          20           for i,y in enumerate(signal):
   218                                                     # print(f"i={i},t[{i}]={time[i]} y[{i}]={y}")
   219    144.2 MiB      0.0 MiB          18               lower_bound = max(time[i]-2*simClass.sigma_smoothing,0) # 2 sigma away backward
   220    144.2 MiB      0.0 MiB          18               upper_bound = min(time[i]+2*simClass.sigma_smoothing,max(time)+2*simClass.sigma_smoothing) # 2 sigma away forward
   221                                                     # MAKE NEW DATA CENTERED AROUND PULSE
   222    144.2 MiB      0.0 MiB          18               if lower_bound < max(t_binned): # if already binned
   223    144.2 MiB      0.0 MiB          16                   lower_bound = t_binned[np.digitize(lower_bound, t_binned)]+simClass.output_bin_width/2
   224    144.2 MiB      0.0 MiB          18               cur_x = np.arange(lower_bound,upper_bound,simClass.output_bin_width)+simClass.output_bin_width/2
   225                                                     # print(f"cur_x from {lower_bound}-->{upper_bound}", cur_x)
   226                                                     # ADD DATA IF NEEDED
   227    144.2 MiB      0.0 MiB         293               for x in cur_x:
   228    144.2 MiB      0.0 MiB         275                   if x > max(t_binned): 
   229    144.2 MiB      0.0 MiB          35                       t_binned.append(x)
   230    144.2 MiB      0.0 MiB          35                       y_binned.append(0)
   231    144.2 MiB      0.0 MiB         240                   elif (np.digitize(x, t_binned)-1 > 0) and (np.digitize(x, t_binned) < len(t_binned)):
   232    144.2 MiB      0.0 MiB         237                       index = np.digitize(x, t_binned)
   233    144.2 MiB      0.0 MiB         237                       if abs(t_binned[index]-t_binned[index-1]) > simClass.output_bin_width:
   234                                                                 t_binned.insert(index, x) # check if need -1 or just np.digitize()
   235                                                                 y_binned.insert(index, 0) # check 
   236                                                     # GET INDICIES
   237    144.2 MiB      0.0 MiB         392               index_lower = [i for i,t in enumerate(t_binned) if t >= lower_bound][0] # left edge in time binned
   238    144.2 MiB      0.0 MiB         392               index_upper = [i for i,t in enumerate(t_binned) if t <= upper_bound][-1] # right edge in time binned
   239                                                     # GAUSSIAN SMOOTH
   240    144.2 MiB      0.0 MiB          18               gaussian = norm.pdf(t_binned[index_lower:index_upper], loc=time[i], scale=simClass.sigma_smoothing)*simClass.sigma_smoothing*y/4
   241                                                     # ADD TO CORRECT BINS
   242    144.2 MiB      0.0 MiB         266               for i,y_add in enumerate(gaussian):
   243    144.2 MiB      0.0 MiB         248                   if y_binned[index_lower+i]+y_add < simClass.max_pmt_current_output:
   244    144.2 MiB      0.0 MiB         248                       y_binned[index_lower+i] += y_add
   245                                                         else:
   246                                                             y_binned[index_lower+i] = simClass.max_pmt_current_output
   247                                         
   248    144.2 MiB      1.1 MiB           2           df = pd.DataFrame({'time':t_binned,'current':y_binned}).sort_values(by=['time'])
   249    144.2 MiB      0.0 MiB           2           print("Formatting PWL dataframe...")
   250    144.2 MiB      0.0 MiB           2           fill_data = []                                                                      # declare empty array
   251                                                 # begin padding data at time 1/5th bin width before first time stamp
   252    144.2 MiB      0.0 MiB           2           fill_data.append([df['time'].iloc[0]-simClass.output_bin_width/5,0])                    # add zero at beginning
   253    144.2 MiB      0.0 MiB          37           for i in range(len(df['time'])-1):                                                        # for each time index
   254    144.2 MiB      0.0 MiB          35               if abs(df['time'].iloc[i]-df['time'].iloc[i+1]) > simClass.output_bin_width:        # if dt between signals is greater than minimum bin width
   255    144.2 MiB      0.0 MiB           2                   fill_data.append([df['time'].iloc[i]+simClass.output_bin_width/5,0])            # place zero after current signal
   256    144.2 MiB      0.0 MiB           2                   fill_data.append([df['time'].iloc[i+1]-simClass.output_bin_width/5,0])          # place zero before next signal
   257    144.2 MiB      0.0 MiB           2           fill_data.append([df['time'].iloc[-1]+simClass.output_bin_width/5,0])                   # add zero at end
   258    144.2 MiB      0.0 MiB           2           fill_data = np.array(fill_data)
   259    144.2 MiB      0.0 MiB           2           fill = pd.DataFrame(fill_data, columns=['time','current'])
   260    144.2 MiB      0.3 MiB           2           df = pd.concat([fill, df], ignore_index=True).sort_values(by=['time']).reset_index(drop=True)
   261    144.2 MiB      0.3 MiB           2           df['time'] = df['time']/1e12
   262    144.2 MiB      0.0 MiB           2           df = df[['time', 'current']] # need this for LTSpice PWL current input file to work
   263    144.2 MiB      0.2 MiB           2           df.to_csv('monte_carlo_input'+str(simClass.num_particles)+'ch'+str(ch)+'_'+str(datetime.now().strftime('%m_%d_%Y'))+'.txt', float_format='%.13f', header=False, index=False, sep=' ') # PWL file formatting
   264    144.2 MiB      0.0 MiB           1       print("Done!")
```

## Changing the Pool commands 
- cpu count * 2 
- max tasks per child = 2 
- running with cProfile instead of memory_profiler 
    - checks the time it takes to run a program
- discovered multiprocessing copies over the overhead needed per process aka the packages for each interpreter
