#############################
# AESOP-Lite Monte Carlo (2024)
# Plotter
# Created by Liam Branch, Michelle Picardo, Robert Johnson
# UCSC Alumni + UCSC Professor
#############################

import numpy as np
# import utility as np
import h5py
# from pandas import DataFrame, read_csv, concat
from tqdm import tqdm
from random import uniform
from time import perf_counter_ns
from datetime import timedelta, datetime
from multiprocessing import Pool, cpu_count, freeze_support, Manager, set_start_method

# Import Main Modules
from simulation import Simulation 
from ltspice_wrapper import LTSpiceWrapper
from plot import Plotter

if __name__ == '__main__':

    ##########################################
    # DECLARE SIMULATION AND PLOTTER CLASSES
    ##########################################
    set_start_method('fork')
    sim = Simulation()
    # plot = plotter(sim)

    #####################
    # RUN SIMULATION 
    #####################
    sim.max_simulated_reflections = 8
    sim.mean_free_path_scints = 0.01
    # sim.mean_free_path_scints = 0.00024 # cm -> 2.4 micrometers
    # sim.num_particles = 4000
    sim.run(1)
    # sim.to_csv(output_both=True)
    # print("hi! run this command to view memory w.o imports")
    # print("mprof run -M -C <file>")
    # print("mprof plot -s")

    ###############################################################
    # RUN LTSPICE AND CALCULATE TIME OF FLIGHT --> SAVE TO FILE
    ###############################################################
    # sim.ltspice(filedate='07_12_2023',filenum=1)
    # sim.calc_ToF(filedate='07_12_2023',filenum=1)
    # sim.save_ToF()

    #########################################
    # LOAD CORRECTED MODEL AND PLOT EXTRA DATA
    #########################################
    # plot.load_extradata(filename='monte_carlo_extradata4000chT1_07_11_2023.txt')
    # plot.plot_xydistance_distr()
    # plot.plot_distPMT_proptime()
    # plot.load_ToF(1, filename='result_1_of_1_07_12_2023.txt')
    # plot.correct_tof()
    # sim.load_ToF(3858, filename='result_3858_of_4000_07_02_2023.txt')
    # sim.plotToF()

