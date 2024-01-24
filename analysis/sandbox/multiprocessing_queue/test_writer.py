# test writer
from multiprocessing import Queue, Pool, cpu_count, Manager
import h5py
from numpy import linspace,repeat,arange,array
from sudo_particle_path import create_pp, array_size



def worker(in_file, i, q):
    with h5py.File(in_file, 'r') as f:
        data = list(f.keys())
        for d in data:
            q.put(f[d][i])
    # print("put", data)

def listener(out_file, q):
    '''listens for messages on the q, writes to file. '''
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
                # f.write(str(count)+': '+str(m) + '\n')
                f.write(str(m) + '\n')
            # f.flush()



if __name__ == '__main__':
    print('Create sudo particle path data with data size: ', array_size)
    pp_file = create_pp(print_keys=True) # returns the name of the file

    manager = Manager()
    q = manager.Queue()
    pool = Pool(cpu_count())

    #put listener to work first
    watcher = pool.apply_async(listener, ('out1.txt', q))

    #fire off workers
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