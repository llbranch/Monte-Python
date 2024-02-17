# test writer
# loc: https://github.com/llbranch/Monte-Python/blob/master/test_writer.py
from multiprocessing import Pool, cpu_count, Manager
import h5py
from sudo_particle_path import create_pp, array_size
# https://superfastpython.com/multiprocessing-pool-share-with-workers/

# Q's are FIFO 

def worker(in_file, item_array, q): 
    '''
    io reading
        - read an item from a key in the file
        - Rather than create a sudo function SCINT I'm using read 
            to gather the parameters from the file 
    '''
# def worker(args):
#     in_file, item_array, q = args
    
    with h5py.File(in_file, 'r') as f: 
        for key in f.keys(): 
            q.put(f[key][item_array])

    # print('Worker: Done', flush=True)


def listener(out_file, q):
    '''
    io writing 
        - listens for messages on the q, writes to file. 
    '''
    with open(out_file, 'w') as f: 
        while 1: 
            item = q.get(timeout=.0003)
            if item == 'kill':
                f.write('killed')
                break
            if item is not None:
                f.write(str(item) + '\n')

# callback function
def custom_callback(result):
	print(f'Got result')

# error callback function
def handler(error):
    print(error, flush=True)





if __name__ == '__main__':

    print('Create sudo particle path data with data size: ', array_size)
    pp_file = create_pp(print_keys=False) # returns the name of the file
    
    # Required & with manager does not work
    manager = Manager()
    q = manager.Queue()
    # q2 = manager.Queue()

    with Pool(processes=cpu_count() + 2) as pool: 
        
        # start watching the Q
        watcher1 = pool.apply_async(listener, (f'out1.txt', q), callback=custom_callback, error_callback=handler)
        watcher2 = pool.apply_async(listener, (f'out2.txt', q), callback=custom_callback, error_callback=handler)
        watcher2 = pool.apply_async(listener, (f'out3.txt', q), callback=custom_callback, error_callback=handler)
        watcher2 = pool.apply_async(listener, (f'out4.txt', q), callback=custom_callback, error_callback=handler)
        watcher2 = pool.apply_async(listener, (f'out5.txt', q), callback=custom_callback, error_callback=handler)
  
        # Version 1
        # Use imap_unordered to iterate over results as they become available
        # result_iterator = pool.imap_unordered(worker, [(pp_file, i, q) for i in range(array_size)])
        # result_iterator = list(result_iterator)
    
        jobs = []
        for i in range(array_size):
            job = pool.apply_async(worker, (pp_file, i, q))
            jobs.append(job)
        
        for j in jobs:
            job.get()



        pool.close()
        pool.join()
        q.put('kill')
        

    print("it worked")
    print(f"Expected # of lines {array_size} x 12: ", array_size*12 )
