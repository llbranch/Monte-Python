# test writer
# loc: https://github.com/llbranch/Monte-Python/blob/master/test_writer.py

# multiprocessing 
from multiprocessing import Pool, cpu_count, Manager, current_process
# file system
import h5py
# calculations
from numpy import arange

from sudo_particle_path import create_pp, array_size


     
    
            
def worker(in_file, item_array, q):
    '''
    io reading
        - read an item from a key in the file
        - Rather than create a sudo function SCINT I'm using read 
            to gather the parameters from the file 
    '''
    name = current_process().name
    # open an existing file to read 
    with h5py.File(in_file, 'r') as f:

        # Old Version: 
        # data = list(f.keys()) 
        # for d in data:
        #     q.put(f[d][i]) # file[key][item_array]: for key in keys put item_array

        # My version
        for key in f.keys(): 
            q.put(f[key][item_array], block=False)
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
            m = q.get_nowait()
            # m = q.get(block = False)
            # print(m)
            if m == 'kill':
                f.write('killed')
                break
            if m != 'END':
                count+=1
                # f.write(name + '-')
                # f.write(str(count)+': '+str(m) + '\n')
                f.write(str(m) + '\n')

                


if __name__ == '__main__':

    print('Create sudo particle path data with data size: ', array_size)
    pp_file = create_pp(print_keys=True) # returns the name of the file
    # check_contents(pp_file)
    Manager = Manager()
    q = Manager.Queue()

    with Pool(cpu_count() - 1) as pool: 

        #put listener to work first
        watcher = pool.apply_async(listener, ('out.txt', q,))

        #fire off workers
        # this is bound by the length of the data arrays
        jobs = []
        for i in range(array_size):
            job = pool.apply_async(worker, (pp_file, i, q))
            jobs.append(job)


        for j in jobs:
            j.get()

        pool.close()
        pool.join()
        # print(dir(pool.apply_async(worker, (pp_file, i, q))))
        q.put('kill')


        print("it worked")
        print(f"Expected # of lines {array_size} x 12: ", array_size*12 )



