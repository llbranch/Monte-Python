# test writer
# loc: https://github.com/llbranch/Monte-Python/blob/master/test_writer.py

from multiprocessing import Queue, Pool, cpu_count, Manager
import h5py
from numpy import linspace,repeat,arange,array
result_queue = Queue()
filename = 'out.txt'
in_file = 'in.hdf5'
def worker(i, q):
    with h5py.File(in_file, 'r') as f:
        data = list(f.keys())
        for d in data:
            q.put(f[d][i])
    # print("put", data)

def listener(q):
    '''listens for messages on the q, writes to file. '''
    count = 0
    with open(filename, 'w') as f:
        while 1:
            m = q.get()
            print(m)
            if m == 'kill':
                f.write('killed')
                break
            if m != 'END':
                count+=1
                f.write(str(count)+': '+str(m) + '\n')
            # f.flush()



if __name__ == '__main__':

    manager = Manager()
    q = manager.Queue()
    pool = Pool(cpu_count() + 2)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,))

    #fire off workers
    jobs = []
    for i in range(100):
        job = pool.apply_async(worker, (i, q))    
        jobs.append(job)
    
    for j in jobs:
        j.get()
    
    q.put('kill')

    pool.close()
    pool.join()
    
    
    # with Pool(processes=cpu_count()-1) as pool:
    #     pool.starmap(processHugeData, (range(10),repeat(result_queue, 10)))
    print("it worked")