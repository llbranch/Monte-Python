# AESOP-Lite Monte Carlo Simulation
### By Liam Branch & Robert Johnson
Simulating the propagation of scintillated photons in a circular plastic scintillator into a photomultiplier tube function to generate signal data. This signal data will then be used in a SPICE simulation program to observe any possible unwanted behavior.

## How to use SLURM and the UCSC Hummingbird Cluster to run the Simulation
All information on Hummingbird and `SLURM` job scheduling for HPC clusters can be found [here](https://hummingbird.ucsc.edu/).

First we sign in to Hummingbird via the command line:
```bash
$ ssh <userid>@hb.ucsc.edu
$ yes
$ <enter gold password>
```
Use the `Hummingbird_Run/monte-carlo-tasks.mpi` file in this repo to schedule the properties of your task to submit. An bash script file of this format is below:
```bash
#!/bin/bash
#SBATCH --job-name=monte-carlo-sim    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END FAIL, ALL)
#SBATCH --mail-user=<userid>@ucsc.edu # Where to send mail	
#SBATCH --nodes=1                     # Use one node
#SBATCH --ntasks=1                    # Run a single task	
#SBATCH --cpus-per-task=16            # Number of CPU cores we use per task
#SBATCH --output=parallel_%j.out      # Standard output and error log
#SBATCH --mem=10G                     # 10GB minimum memory

# Note: module load of gnu is optional (its the default)
# You may need other modules, however, to run your program. 
# Note: use module avail to check available modules on Hummingbird
module load python-3.6.2 # load python from available slurm modules
python3 multiprocess_AESOP_0423.py # run script with main function defined
```
Now we edit the `monte-carlo-tasks.mpi` to send to the your correct `<userid>@ucsc.edu` email and set `--cpus-per-task=16` to the max per the node you are trying to use. There are four public options for nodes (called partitions or queues) students can use. 
- **128x24** --> 128GB of memory - 24 CPU Cores 
- **256x44** --> 256Gb of memory - 44 CPU Cores
- **96x24gpu** --> 96Gb of memory - 24 CPU Cores - 14336 GPU Cores - 64Gb of HBM2 memory (4 NVIDIA Tesla P100 GPUs)
- Instruction (for classes)

This means to get the max available CPU Cores for our script to run quickly we want max `24` or max `44` depending on the node you specify. To specify a specific nodes other than the 128x24 add another line like below to your bash script:
```bash
#SBATCH -p 256x44
```
Make sure your `.mpi` file is in the same directory as your `multiprocess_AESOP_0423.py` file. Finally to run this bash script we're going to use `sbatch` command rather than `srun` so our task runs in the background rather than interactively. 
```bash
sbatch monte-carlo-tasks.mpi
```
You can check the status of the job using `squeue` or by visiting the [official Hummingbird website](https://hummingbird.ucsc.edu/current-usage/)