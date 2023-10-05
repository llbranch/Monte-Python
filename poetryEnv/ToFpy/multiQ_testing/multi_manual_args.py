"""Takes arguments for multiprocessing"""
import multiprocessing

def worker_function_v2(num): 
    """ argument worker function """
    print('Worker function: ', num) # no returns 

if __name__ == '__main__': 
    jobs = [] 
    for i in range(5):
        process_i = multiprocessing.Process(target=worker_function_v2, args=(i,))
        jobs.append(process_i)
        process_i.start()
    
    # process_i.join() # include this line to reconnect jobs before proceeding to next line in the code
    
    print(len(jobs))