# AESOP-Lite Monte Carlo Simulation
### By Liam Branch & Robert Johnson
Simulating the propagation of scintillated photons in a circular plastic scintillator into a photomultiplier tube function to generate signal data. This signal data will then be used in a SPICE simulation program to observe any possible unwanted behavior.

## How to use SLURM and the UCSC Hummingbird Cluster to run the Simulation
Use the `Hummingbird_Run/monte-carlo-tasks.mpi` file to schedule the properties of your task to submit. We're going to use `sbatch` command rather than `srun` so our task runs in the background rather than interactively. 
```
ssh <youruserid>@hb.ucsc.edu
yes
<enter gold password>
```
Now we edit the `monte-carlo-tasks.mpi` to send to your `<userid>@ucsc.edu` email and set `--cpus-per-task=16` to the max per the node you are trying to use. Options are max `24` or max `44` depending on the node you use. To change nodes make another similar line like below:
```
#SBATCH -p 256x44
```


You can check the status of the job using `squeue` or by visiting the [official Hummingbird website](https://hummingbird.ucsc.edu/current-usage/)