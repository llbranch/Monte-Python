# Testing by Michelle P 9/4/23
# Keep Multiprocessing 
#   add .join and .close 

# Installs Needed 
import numpy as np #analysis
# np calls: 
#       - cos
#       - sin
#       - radians
#       - linspace
#       - floor
#       - log10
#       - abs
#       - linalg.norm
#       - sqrt
#       - sum
#       - dot
#       - array
#       - min
#       - cross
#       - arcsin
#       - random.random
#       - pi
#       - float64
#       - random.poisson
#       - zeros
#       - repeat
#       - asarray
#       - digitize
#       - arange
#       - gradient
#       - diff
#       - polyfit
#       - meshgrid
#       - ones
#       - std
#       - mean
#       - count_nonzero
#       - concatenate
#       - c_
#       - unique
#       - append
#       - exp
#       - quantile
#       - diag
#       - log
from pandas import DataFrame, read_csv, concat #analysis
import h5py
from tqdm import tqdm #progressbar

# Native 
from random import uniform #generator 
from time import time, perf_counter_ns #timing 
from datetime import timedelta, datetime #timing
# Native + memory
from memory_profiler import profile #memory checks
# Testing Multiprocessing 
from multiprocessing import Pool, cpu_count, freeze_support, Manager, Queue
from concurrent.futures import ProcessPoolExecutor

# Vscode highlights and helps navigate the regression in the program 


# Class Simulation: self.function()/self.variable() called within itself
### No self calls 
#   init(self) : constants -
#   round_to_sig(self,float) : int -
#   normalize(self, float): float-
#   lg_condition(self, x1, x2, x3): y-
#   distance_circle(self,x1,x2,x3,x4,x5): ?-
#   distance_plane(self, x1,x2,x3,x4): ?-

### Calls self Constants ONLY
#   scint_condition(self, x1,x2,x3): y-
#   photontoElectrons(self, x): float ? -

### Calls self.other_functions
#   distance_solver(self, x1,x2,x3,x4,x5,x6,x7,x8,x9): ?,Bool-
#   - distance_circle(self,x1,x2,x3,x4,x5): ?-
#   - distance_plane(self, x1,x2,x3,x4): ?-

#   photon_interaction(self, x1,x2): ?,Bool-
#   - normalize(self, float): float-
#   - self.constants -

#   n_vec_calculate(self, x1,x2,x3,x4,x5): array ? 
#   - normalize(self, float): float-

#   particle_path(self, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x12,x13): array ?- 
#   - round_to_sig(self,float) : int -
#   - self.constants -
#   - scint_condition(self, x1,x2,x3): y-
#   - lg_condition(self, x1, x2, x3): y-


#   scintillator_monte_carlo(self,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x12 ): bool, constants ... -
#   - distance_solver(self, x1,x2,x3,x4,x5,x6,x7,x8,x9): ?,Bool-
#   - n_vec_calculate(self, x1,x2,x3,x4,x5): array ? -
#   - photon_interaction(self, x1,x2): ?,Bool-

#   particle_task(self, x1,x2,x3): void --> (self, multi)
#   - particle_path(self, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x12,x13): array ? 

#   scint_taskT1(self, x1,x2,x3,x4): bool, constants... --> (self, point_i, time_i)
#   - scintillator_monte_carlo(self,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x12 ): bool, constants ... 

#   scint_taskT4(self, x1,x2,x3,x4): bool, constants... --> (self, point_i, time_i)
#   - scintillator_monte_carlo(self,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x12 ): bool, constants ... 

#   run_worker_T1(self, i,q): ?
#   - scint_taskT1(self, x1,x2,x3,x4)

#   run_worker_T4(self, i,q): ?
#   - scint_taskT4(self, x1,x2,x3,x4)

#   listener(self, q, filename): 
#   - h5py commands

#   Pool  ###################################################
#   run(self, x1*,x2*): void --> (self, *arg, **kwargs)
#   - self.constants 
#   - particle_task(self, x1,x2,x3): void --> (self, multi)
#   - listener(self, q, filename): 
#   - run_worker_T1(self, i,q): ?
#   - run_worker_T4(self, i,q): ?
#   - photontoElectrons(self, x): float ? 
# ###########################################################

# Plotting and other functions without other function calls 

#   to_csv(self, *kwargs): cases 2(?, void)
#   time_at_thresh(self, x1,x2,x3,x4,x5): array ? 
#   ToF_finalize(self, x1,x2,x3): void
#   ltspice(self, x1,x2,x3): void 
#   calc_ToF(self, x1,x2): void
#   save_ToF(self, x1): void 

# ------------------ end of 1st class -----------------


# Data for saving 
# propagation time
# propagation dist
# pmt dist prop time
# end point dist
# interaction 
# extra times 