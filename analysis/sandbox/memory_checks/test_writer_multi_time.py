# test writer
# loc: https://github.com/llbranch/Monte-Python/blob/master/test_writer.py


# multiprocessing 
from multiprocessing import Queue, Pool, cpu_count, Manager, current_process
# file system
import h5py
# time check
from cProfile import Profile # run command: snakeviz report_time.prof 
from pstats import SortKey, Stats # used with cProfile

from sudo_particle_path import create_pp, array_size


'''
Questions to answer: 
    1. How many workers? 
    2. How many listeners? 
    3. What are we using to verify the program is working? 

Comments: 
    1. Assuming we have a way to verify the work is - working, 
       why not store in binary? 
# '''

            
def worker(in_file, item_array, q):
    '''
    io reading
        - reads (?) files with particle path data 
        - passes that data to the scintillator pool
        - the data from the worker should be sent to a q
    '''
    name = current_process().name
    # open an existing file to read 
    with h5py.File(in_file, 'r') as f:

        # keep the list in RAM 
        # data = list(f.keys()) 
        # for d in data:
        #     q.put(f[d][i]) # file[key][item_array]: for key in keys put item_array

        for key in f.keys(): 
            q.put(f[key][item_array])
            # q.put(f[key][()]) # prints the entire array

    # print(name)

            
def listener(out_file, q):
    '''
    io writing 
        - listens for messages on the q, writes to file. 
    '''
    name = current_process().name

    count = 0
    with open(out_file, 'w') as f:
        while 1:
            m = q.get()
            # print(m)
            if m == 'kill':
                f.write('killed')
                break
            if m != 'END':
                count+=1
                # f.write(name + '-')
                # f.write(str(count)+': '+str(m) + '\n')
                f.write(str(m) + '\n')

            # f.flush()
                





if __name__ == '__main__':

    with Profile() as profile:

        pp_file = create_pp(print_keys=False) # returns the name of the file
        # check_contents(pp_file)

        manager = Manager()
        q = manager.Queue()

        with Pool(cpu_count() - 1) as pool: 

            #put listener to work first
            watcher = pool.apply_async(listener, ('out.txt', q,))

            #fire off workers, why is upper_bound the length of the item_array? 
            jobs = []
            for i in range(array_size):
                job = pool.apply_async(worker, (pp_file, i, q))    
                jobs.append(job)
            
            for j in jobs:
                j.get()
            
            q.put('kill')

            pool.close()
            pool.join()
            print("it worked")
            print(f"Expected # of lines {array_size} x 12: ", array_size*12 )

    results = Stats(profile)
    results.strip_dirs()
    results.sort_stats(SortKey.CALLS)
    results.dump_stats("time_output/report_time.prof")



