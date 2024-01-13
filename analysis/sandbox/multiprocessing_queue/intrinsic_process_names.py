""" Demonstrates that each process has an intrinsic name """

import multiprocessing
import time

def my_worker(): 
    name = multiprocessing.current_process().name
    print(name, 'Starting...')
    time.sleep(2)
    print(name, 'Ending')
    # no return 

def my_service():
    name = multiprocessing.current_process().name
    print(name, 'Starting...')
    time.sleep(2)
    print(name, 'Ending')
    # no return 

if __name__ == '__main__': 

    # name a service, though here worker and service are identical
    service = multiprocessing.Process(name = 'My assigned service', 
                                      target=my_service,)
    
    # name a worker
    worker_1 = multiprocessing.Process(name='New worker 1',
                                       target= my_worker,)
    
    # Use default name
    worker_2 = multiprocessing.Process(
        # default name 
        target=my_worker,)

    service.start()
    worker_1.start()
    worker_2.start()