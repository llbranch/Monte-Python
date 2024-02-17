from multiprocessing import Process, Queue, Manager, Pool, cpu_count
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
    # open an existing file to read 
    with h5py.File(in_file, 'r') as f:

        # My version
        # for key in f.keys(): 
        #     q.put(f[key][item_array])
            # q.put(f[key][()]) # prints the entire array


        # # LAZY WAY Process and put data into the queue without loading the entire array
        # q.put(process_data(f[key][item_array]))


def listener(out_file, q):
    '''
    io writing 
        - listens for messages on the q, writes to file. 
    '''

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
                # for item in m:
                #     f.write(str(item) + '\n')       


if __name__ == '__main__':

    print('Create sudo particle path data with data size: ', array_size)
    pp_file = create_pp(print_keys=True)  # returns the name of the file
    # check_contents(pp_file)

    manager = Manager()
    queues = [manager.Queue() for _ in range(12)]

    processes = []

    # Create and start listener processes
    for i in range(12):
        listener_process = Process(target=listener, args=(f'data/out_{i}.txt', queues[i],))
        listener_process.start()
        processes.append(listener_process)

  
    

    with Pool(cpu_count() - 1) as pool:

        # Fire off workers
        for i in range(array_size):
            for j in range(12):
                pool.apply_async(worker, (pp_file, i, queues[j]))

        # Wait for all workers to finish
        # pool.close()
        # pool.join()

        # Signal workers to stop
        for j in range(12):
            queues[j].put('kill')

        print("it worked")

        # Wait for listener processes to finish
        for listener_process in processes:
            listener_process.join()

    print(f"Expected # of lines {array_size} x 12: ", array_size * 12)