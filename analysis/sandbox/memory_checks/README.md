# How to use
The following is specific to running files in this directory.
## Check memory 
`memory_profiler`\
To save memory data to a `.log` include the following: 
```py
from memory_profiler import profile, LogFile
import sys  # used with memory_profiler

# set IO stream for profiler
fp = open("mem_output/memory_profiler.log", "w+")

# add decorator to desired functions
@profile(precision=4, stream=fp)

if __name__ == "__main__": 
    sys.stdout = LogFile('memory_profile_log',reportIncrementFlag=True)

```
Run the following in the command line 
```bash
# The following can be run on a python file WITHOUT decorators

$ poetry run mprof run -M <file.py> # shows individual memory plots (forked processes)

# or 

$ poetry run mprof run -M -C <file.py> # shows sum and individual (forked processes)

# and

$ mprof plot -o plot.png # this plots the data

# or

# this plots the data & trend line <0 is a potential memory leak
$ mprof plot -s -o plot.png 

```

## Check timing 
`cProfile`\
To check the timing of a file include the following: 
```py
# Add this within the file 
from cProfile import Profile
from pstats import SortKey, Stats # organizes the output of Profile

# Functions 


if __name__ == "__main__": 
    with Profile() as profile:
        # execution 
    
    # store statistical results of profile
    results = Stats(profile)
    # removed the extraneous path from all the module names
    results.strip_dirs() 
    # sorts
    results.sort_stats(SortKey.CALLS)
    # store with .prof
    results.dump_stats("time_output/report_time.prof")
```
To view the time analysis run this cmd. Note, package `snakeviz` is required.
```bash
$ poetry run python <file.py>

$ snakeviz time_output/report_time.prof
```