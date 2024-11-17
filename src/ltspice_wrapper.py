#############################
# AESOP-Lite Monte Carlo (2024)
# LT Spice Run and Parsing
# Created by Liam Branch, Michelle Picardo, Robert Johnson
# UCSC Alumni + UCSC Professor
#############################
class LTSpiceWrapper():
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

    """LTSpice Command to Analyze, Simulate and Calculate TOF"""
    def ltspice(self, filedate=None, filenum=None, filesep=None):
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
        sep_part = self.seperation_time
        if filedate is not None: # Take input given correct format
            date = filedate
        if filenum is not None: # Take input given correct integer number
            num_part = int(filenum)
        if filesep is not None:
            sep_part = int(filesep)
        filename_ch1 = os.path.abspath('monte_carlo_input'+str(num_part)+'ch1_'+str(date)+'.txt')
        filename_ch4 = os.path.abspath('monte_carlo_input'+str(num_part)+'ch4_'+str(date)+'.txt')
        for filename,strname in zip([filename_ch1,filename_ch4],['ch1','ch4']):
            print('PWL file='+str(filename))
            LTC.set_element_model('I1', 'PWL file='+str(filename))
            LTC.add_instructions("; Simulation settings", f".tran 0 {int(round(num_part*sep_part/1e6,0))}.5u 0 0.002u") # fix this to adjust for time seperation
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
            df = DataFrame({'t':np.abs(t),'V':compOut}).sort_values(by='t')
            df.to_csv('output'+str(num_part)+strname+'_'+str(date)+'.txt', header=False, index=False)
            # implement csv creation!
            # Clean up extra files
            os.remove(f"PHAReduced_{strname}.log")
            os.remove(f"PHAReduced_{strname}.op.raw")
            os.remove(f"PHAReduced_{strname}.raw")
            os.remove(f"PHAReduced_{strname}.net")
            # Remove LTSpice Object
            del LTR

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
        first_particle = True
        condition_0 = (grad > 0) & (grad < 0.001)
        condition = tdiff > self.seperation_time/1e12/10 # check if next point is 10ns away
        count = 0
        start_index = 0
        for i in range(len(tdiff)): # for number of particles we expect
            if count > num-1:
                return np.array(out)
            if condition[i] or condition_0[i]: # if condition is true at index i
                if first_particle:
                    condition_0 = False
                    first_particle = False
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
        ch1 = read_csv(filename_ch1, names=['t', 'V'], sep=',')
        ch4 = read_csv(filename_ch4, names=['t', 'V'], sep=',')
        ch1ToF = self.time_at_thresh(ch1['t'],ch1['V'], num_part, self.CMOS_thresh, 1)
        ch4ToF = self.time_at_thresh(ch4['t'],ch4['V'], num_part, self.CMOS_thresh, 4)
        self.ToF_finalize(ch1ToF,ch4ToF) # Calculated correct time of flight
        print(DataFrame(self.FinalToF).describe())

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
        DataFrame({'Time-of-Flight [s]':self.FinalToF}).to_csv(file, index=False)