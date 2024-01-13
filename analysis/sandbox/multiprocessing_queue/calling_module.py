""" Demonstrates the ability to call a module and run a process """

import multiprocessing
from mod_worker_function import worker_function
# call function 

if __name__ == '__main__': 
    jobs = []
    for i in range(5): 
        process_i = multiprocessing.Process(target=worker_function)
        jobs.append(process_i)
        process_i.start()

    