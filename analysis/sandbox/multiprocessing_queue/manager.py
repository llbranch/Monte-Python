# SuperFastPython.com
# example of sharing a process pool among processes
from multiprocessing import Manager
 
# error callback function
def handler(error):
    print(error, flush=True)
 
# task executed in a worker process
def task(pool):
    # report a message
    print(f'Pool Details: {pool}')
 
# protect the entry point
if __name__ == '__main__':
    # create a manager
    with Manager() as manager:
        # create and configure the process pool
        with manager.Pool() as pool:
            # issue a task to the process pool
            pool.apply_async(task, args=(pool,), error_callback=handler)
            # close the pool
            pool.close()
            # wait for all issued tasks to complete
            pool.join()