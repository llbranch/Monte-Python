# Outputs from runs: 1 particle 8 reflections 
1. Ran in 1 .py file no changes

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
Photons generated 50704.0
Photons in T1: 16877.0 and Photons in T4: 33827.0
T1 Photon Propagation working...
Done!
T4 Photon Propagation working...
Done!
Unzipping reuslts into arrays...
TIME ANALYSIS:
Generation of Particles     0:00:00.185221
Simulation of Photon Travel 0:00:11.085438
Total Time Elapsed:         0:00:11.332092
RESULTS SUMMARY:
HITS on T1 9
RATIO T1   total photons 16877.0 total incident photons 9 ratio=1875.22
HITS on T4 11
RATIO T4   total photons  33827.0 total incident photons 11 ratio=3075.18
DISTANCE: 
```

2. Ran with mainTof.py and modTof.py; no changes to classes 
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
Photons generated 50322.0
Photons in T1: 16928.0 and Photons in T4: 33394.0
T1 Photon Propagation working...
Done!
T4 Photon Propagation working...
Done!
Unzipping reuslts into arrays...
TIME ANALYSIS:
Generation of Particles     0:00:00.236504
Simulation of Photon Travel 0:00:11.302593
Total Time Elapsed:         0:00:11.572675
RESULTS SUMMARY:
HITS on T1 18
RATIO T1   total photons 16928.0 total incident photons 18 ratio=940.44
HITS on T4 16
RATIO T4   total photons  33394.0 total incident photons 16 ratio=2087.12
DISTANCE: 
Exporing to 2 channels...
Smoothing Signals...
Formatting PWL dataframe...
Smoothing Signals...
Formatting PWL dataframe...
Done!
```
3. ran with simple compile; modTof.py is compiled without typesetting 
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
Photons generated 50234.0
Photons in T1: 16908.0 and Photons in T4: 33326.0
T1 Photon Propagation working...
Done!
T4 Photon Propagation working...
Done!
Unzipping reuslts into arrays...
TIME ANALYSIS:
Generation of Particles     0:00:00.181830
Simulation of Photon Travel 0:00:11.204276
Total Time Elapsed:         0:00:11.421281
RESULTS SUMMARY:
HITS on T1 14
RATIO T1   total photons 16908.0 total incident photons 14 ratio=1207.71
HITS on T4 5
RATIO T4   total photons  33326.0 total incident photons 5 ratio=6665.20
DISTANCE: 
Exporing to 2 channels...
Smoothing Signals...
Formatting PWL dataframe...
Smoothing Signals...
Formatting PWL dataframe...
Done!
```